# The algorithm

This is the walkthrough. [`DESIGN.md`](DESIGN.md) covers the decisions —
number types, API shape, what is deliberately not here. Read this one first.

## What a straight skeleton is

Take a polygon. Slide every edge inward, perpendicular to itself, all at the
same speed, keeping each edge straight and keeping the whole thing connected.
The polygon shrinks. Its corners trace out paths. Those paths are the straight
skeleton.

```
    +---------------------------+          +---------------------------+
    |                           |          | \                       / |
    |                           |          |   \                   /   |
    |                           |          |     \_______________/     |
    |                           |    =>    |     /               \     |
    |                           |          |   /                   \   |
    |                           |          | /                       \ |
    +---------------------------+          +---------------------------+

          the input                    the shrinking corners' traces:
                                            the straight skeleton
```

The shrinking outline is called the **wavefront**. Watch it for a rectangle:

```
    +-------------------------+     t=0    the wavefront is the input
    |  +-------------------+  |     t=1    every edge has moved in by 1
    |  |  +-------------+  |  |     t=2
    |  |  |             |  |  |
    |  |  |  o-------o  |  |  |     t=5    the top and bottom edges collide;
    |  |  |             |  |  |            all that is left is a segment
    |  |  +-------------+  |  |
    |  +-------------------+  |
    +-------------------------+

    That final segment is the roof's ridge, and it is part of the skeleton.
```

## Why the wavefront is worth simulating

Between two consecutive edges sits a wavefront **vertex**. As the two edges
slide, that vertex has to stay on both of them at once, so it slides along the
angle bisector — at speed `1 / sin(θ/2)`, faster than the edges themselves when
the corner is sharp.

Crucially, **the vertex moves in a straight line at a constant velocity** for as
long as its two edges are its two edges. So its whole path is determined by one
vector, and the simulation only has to notice the moments where the wavefront's
*structure* changes. Those moments are **events**, and between them nothing
interesting happens. That is what makes this tractable: the simulation hops from
event to event rather than stepping through time.

## The three events

### Edge event — an edge shrinks to nothing

Two adjacent vertices converge and meet. The edge between them vanishes, and
they fuse into one vertex carrying the two outer edges.

```
        \             /                          \       /
         \           /                            \     /
          \_________/          =>                  \   /
          /    e    \                               \ /
         /           \                               o     e is gone; the two
                                                           neighbours have met
```

This is where a triangle's three corners meet at its incenter.

### Split event — a reflex vertex hits an opposing edge

A polygon with a notch has a **reflex** vertex, one whose interior angle exceeds
180°. Reflex vertices move *outward* along their bisector, into the material, and
can run into a wavefront edge on the far side. When that happens the wavefront
tears into two independent loops.

```
      +---------------------------+                +--------------+------------+
      |                           |                |              |            |
      |         reflex v          |                |              |            |
      |            \              |       =>       |              |            |
      +-------+     \    +--------+                +-------+      |   +--------+
              |      \   |                                 |      |   |
              |       \  |                                 |      |   |

      v drives into the far edge            the loop splits; each half now
                                            shrinks on its own
```

Only reflex vertices can split. Convex ones always move into shrinking material
and cannot reach an opposing edge before their own neighbours do. That is why
`skeleton()` on a convex polygon never searches for split events at all.

**Holes work by exactly this mechanism.** A hole is a separate wavefront loop
that *expands* into the material. Its vertices are all reflex when seen from the
material's side (`Polygon::is_reflex` is written to say so), so they split.
When the splitting vertex and the edge it splits are in *different* loops, the
same relinking that tears one loop in two instead **merges** two loops into one.
No special case; the code does not even know which happened.

### Speed change — an edge hits its distance limit

Only in `skeleton_constrained`. See "The one idea" below.

## The one idea

Every edge is a moving line. At time `t`, edge `e`'s line is

```
    { x : normal_e · x = c_e + offset_e(t) }
```

A wavefront vertex must sit on both of its edges' lines at once. So its velocity
`v` is whatever satisfies both:

```
    normal_left  · v = speed_left
    normal_right · v = speed_right
```

Two equations, two unknowns, one 2x2 solve (`Sim::velocity_of`). And that single
solve is the whole crate:

| speed_left | speed_right | what falls out |
|---|---|---|
| 1 | 1 | the classic angle bisector, at `1/sin(θ/2)` — the plain skeleton |
| 1 | 0 | the vertex **slides along** the stopped edge |
| 0 | 0 | the vertex **freezes** |

So a per-edge distance limit is not a special mode. It is just "this edge's
speed drops from 1 to 0 at `t = limit`". `skeleton` and `skeleton_constrained`
run the same simulation; the former simply never schedules a speed change.

This is the *weighted* straight skeleton, restricted to weights in `{0, 1}`.

## Finding the events

**Edge collapse** is one-dimensional, which is the trick worth knowing. Both
endpoints of a wavefront edge lie on that edge's line for the whole simulation,
so there is no need for a 2D intersection test. Project their separation onto
the edge's own direction and ask when it reaches zero. That is linear in `t`, so
it is one divide, with no special cases (`Sim::edge_event`).

**Split** is the expensive one, and how you ask the question decides whether the
whole algorithm works or melts down.

The naive version asks the whole question at once: *when does this reflex vertex
land on a live stretch of some other edge?* That needs the current endpoints of
every candidate edge — and those endpoints belong to a part of the wavefront
arbitrarily far away. It couples everything to everything. See below.

The version here asks a deliberately weaker question (`split_lower_bound`):

> when does this vertex reach some edge's **moving line** — never mind whether it
> lands on a live stretch of it?

That is *only a lower bound* on the true split time. It is never late, though,
which is the only property needed for popping events in time order to stay
correct. And it is cheap to keep true: an edge's wavefront slides along its own
offset track and never leaves it, so the answer depends only on the vertex's
trajectory and the edge's original line. **Nothing happening elsewhere can
change it.**

The real question is settled later, in `handle_split_event`, when the event is
popped. By then `now == t`, every earlier event has been processed, and the
wavefront's shape at `t` is not a forecast but settled fact — so
`live_stretch_at` simply *looks*. If the vertex came down off the end of every
live stretch, no split happens: the edge is struck off (a vertex travels in a
straight line, so it meets that line once and the question is closed for good)
and the next candidate is queued.

The scan is `O(n)` per reflex vertex, and it is the only non-constant step left.

## Staleness, and the meltdown it hid

Events are queued, never removed; obsolete ones are recognised when popped. Two
independent stamps do it, and conflating them loses events:

- `gen`, bumped when a vertex **moves or is relinked** — invalidates events
  computed *from* it.
- `evt`, bumped when a vertex is **rescheduled** — supersedes only its own
  previous event, disturbing nothing else.

The critical property is that an event's `refs` — the vertices its timing was
computed from — stay **O(1)**: the owner and its one neighbour. So `touch`
marks a vertex and its `prev`, and that is provably complete, because exactly
two events are computed from any vertex's geometry: its own, and its
predecessor's edge event watching the two of them converge.

It is worth seeing what happens when that is not true, because this crate did it
the other way first and the result was spectacular. Stamping a split event with
the two endpoints of the edge it was aimed at made a vertex's event depend on a
part of the wavefront arbitrarily far away. Then:

1. any event invalidated events all over the polygon;
2. every reschedule re-registered more dependencies;
3. which made the next event invalidate even more.

It fed back on itself. A 132-vertex comb took 124ms; a 260-vertex one asked for
**27GB** and died. Measured growth: **n^5.5**. The file you are reading claimed
`O(n^2 log n)` at the time — the claim was theory, and nobody had run it.

The fix was not better bookkeeping. It was noticing that a split's *timing never
depended on those endpoints in the first place*, which is what `split_lower_bound`
is. The same comb now runs 1028 vertices in 8.9ms.

A second, smaller one hid behind it: the struck-off list was scanned linearly
*inside* the per-edge loop, so one scan cost `O(n * rejections)`. Sorting it and
binary-searching took a random star from `n^3.0` to `n^2.2`.

Both were bookkeeping, not geometry. That is the argument for measuring.

## The degeneracies that actually bite

Textbook descriptions quietly assume events are distinct and generic. They are
not, on input as ordinary as a square.

### Simultaneous events: the vertex event

A square's four corners reach the centre **at the same instant**. Handle that as
a cascade of two-vertex merges and the first merge leaves a vertex whose two
edges are the square's *opposite* sides — antiparallel, so the velocity solve has
no solution.

The fix is to stop pretending: when several consecutive vertices arrive at one
point together, gather the whole coincident run and retire it in one event
(`Sim::coincident_chain`). One node, one set of arcs, no impossible leftovers.

### Needles

A rectangle's long sides collide head-on. So do the two sides of any strip of
material — a hole sitting 2*d* from a wall pinches the strip between them shut.
The vertex left behind has antiparallel edges. Since it lies on both of their
lines, and two antiparallel lines through a point are the *same* line, its edges
have collided and the material between them is gone. The wavefront has folded
back on itself:

```
    prev  o<---------------------o m          prev  o
          o--------------------->'                  |
    next                                            |   the strip is gone;
                              =>                    |   what remains is
    the two edges now lie on top                    |   one skeleton arc
    of one another                            next  o
```

This one *must* be handled explicitly. The folded edges are parallel, so no edge
event can ever fire on them again: leave it alone and the simulation stalls with
the loop still live, silently dropping every arc it had left to trace.

`Sim::resolve_needle` zips it up. The overlap — from the vertex to whichever
neighbour is nearer — is exactly one skeleton arc, bisecting the two edges that
collided. For a rectangle, that arc *is* the ridge. Emit it, retire the arm,
splice the rest back, and repeat, because zipping one needle shut routinely
exposes the next.

A two-vertex loop stalls for the same reason and is resolved the same way.

## Complexity, and the meltdown that hid in it

Measured, with `cargo run --release --example bench`:

| input | scaling | why |
|---|---|---|
| convex | ~n^1.3 | no reflex vertices, so `split_lower_bound` is never called |
| reflex-heavy | ~n^1.9 to n^2.1 | `O(1)` events, each `O(n)` to schedule |

Space is `O(n)`.

The event count is **linear** — about 5n pops for a comb, 10n for a star — and
each event reschedules `O(1)` vertices. The quadratic term is entirely
`split_lower_bound`'s scan over every edge.

That is worth stating carefully, because an earlier version of this file claimed
`O(n^2 log n)` and the real figure was `n^5.5`, on its way to asking for 27GB at
260 vertices. The claim was theory; nobody had measured. Two things were wrong,
and both were in the *bookkeeping*, not the geometry:

- Split events depended on the endpoints of the edge they were aimed at, so any
  event invalidated events all over the polygon, each reschedule registered more
  dependencies, and it fed back on itself. Fixed by
  [`split_lower_bound`](#finding-the-events): the timing never depended on those
  endpoints in the first place.
- The reject list was scanned linearly *inside* the per-edge loop, making one
  scan `O(n * rejections)`. Sorting it and binary-searching took a star from
  `n^3.0` to `n^2.2`.

Beating `O(n^2)` in the worst case needs the motorcycle graph. Reflex vertices
launch "motorcycles" along their bisectors; split events are exactly where
motorcycles crash, and — crucially — motorcycle traces are **static rays**, so
they can go in a spatial index, which moving wavefront edges cannot. Given the
motorcycle graph, the skeleton follows in `O(n log n)` (Cheng–Vigneron,
Huber–Held). It is a different algorithm, and a much larger one; CGAL and
Surfer2 ship an `O(n^2)` worst case too.

## Where it is not the medial axis

Worth internalising, because it is the most common misconception about straight
skeletons and it will bite you when reading the output.

Both bisect their input. They differ in *what* they bisect:

- a **straight skeleton** bisects edges' infinite **supporting lines**. Every arc
  is straight — that is the whole point, and the source of the name.
- a **medial axis** bisects the nearest **features**, and so grows *parabolic*
  arcs around reflex vertices, where the nearest feature is a point.

They agree exactly when the polygon is convex. Around reflex corners they do not:

```
    In the plus-shape, the centre is at offset 5 — the four arms' walls have all
    travelled 5 to reach it. But the *nearest input feature* to that centre is a
    reflex corner, 7.07 away. The reflex vertex swept out along its bisector at
    1/sin(45°) = 1.41x the wavefront's speed, and 5 * 1.41 = 7.07.

    offset != distance-to-boundary. Never assume it does.
```

This is why `Node::sources` is documented as *"the edges whose faces meet here"*
rather than *"the closest edges"*. The claim it makes is exact, and it is the one
that is actually useful: it is what assigns a roof panel to its wall.

## Faces, and why roofs are free

Each input edge owns one **face**: everything its wavefront swept. The faces tile
the polygon, and `Skeleton::face` recovers one by walking the arcs that name that
edge as a source.

Now lift every node to `z = offset`. Each face becomes **planar**, because every
point on it is `offset` from the same line and height is a linear function of
that distance. So:

- one face = one flat roof panel,
- one arc = a hip, valley, or ridge where two panels meet,
- `offset * pitch` = height.

The `roof` example does nothing but read this off, and asserts each panel really
is planar as it goes.

## Further reading

- Aichholzer, Aurenhammer, Alberts, Gärtner, *A Novel Type of Skeleton for
  Polygons* (1995) — the original.
- Aichholzer, Aurenhammer, *Straight Skeletons for General Polygonal Figures in
  the Plane* (1996) — holes and the general case.
- Felkel, Obdržálek, *Straight Skeleton Implementation* (1998) — the SLAV
  formulation this crate's structure resembles. Note that its split-event
  handling is known to be incomplete on some inputs.
- Huber, Held, *Theoretical and Practical Results on Straight Skeletons of
  Planar Straight-Line Graphs* (2011) — on the degeneracies that actually matter.
