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

### What stops when everything stops

The table's last row has a consequence worth following. A plain skeleton always
finishes by its wavefront shrinking away to nothing — every edge is closing on
another, so something always collapses eventually. Freeze every edge around a
loop, though, and nothing is closing on anything: the loop stops, and stays
exactly where it is, forever.

So the simulation runs out of events with the loop still standing. That is not a
stall to be fixed — it is the answer. What is standing is the input polygon
offset inward by the limit, and `Sim::collect_residual` reads it straight off the
`prev`/`next` links once the queue is empty. Each surviving vertex's `node` is
already the node its arc stopped at, so there is nothing to compute.

It is the only way a vertex can outlive the queue, which is what makes the
reading unambiguous: `velocity_of` returns zero only when *both* of a vertex's
edges have stopped. (A needle's antiparallel pair also sits still, but
`resolve_needle` retires those on the spot rather than leaving them behind.)

See [`DESIGN.md`](DESIGN.md#the-residual-wavefront-is-half-the-answer) for why
this comes back as its own type rather than as more arcs.

## Finding the events

**Edge collapse** is one-dimensional, which is the trick worth knowing. Both
endpoints of a wavefront edge lie on that edge's line for the whole simulation,
so there is no need for a 2D intersection test. Project their separation onto
the edge's own direction and ask when it reaches zero. That is linear in `t`, so
it is one divide, with no special cases (`Sim::edge_event`).

**Split** is the expensive one, and how you ask the question decides whether the
whole algorithm works.

The obvious version asks the whole question at once: *when does this reflex
vertex land on a live stretch of some other edge?* That needs the current
endpoints of every candidate edge — and those endpoints belong to a part of the
wavefront arbitrarily far away, which couples everything to everything.

The version here asks a deliberately weaker question (`Sim::scan_for_split`):

> when does this vertex reach some edge's **moving line** — never mind whether it
> lands on a live stretch of it?

That is *only a lower bound* on the true split time. It is never late, though,
which is the only property needed for popping events in time order to stay
correct. And it is cheap to keep true: an edge's wavefront slides along its own
offset track and never leaves it, so the answer depends only on the vertex's
trajectory and the edge's original line. **Nothing happening elsewhere can
change it.** That is what keeps an event's dependencies `O(1)` — see
"Staleness" below.

The real question is settled later, in `Sim::handle_split_event`, when the event
is popped. By then `now == t`, every earlier event has been processed, and the
wavefront's shape at `t` is not a forecast but settled fact — so
`live_stretch_at` simply *looks*. If the vertex came down off the end of every
live stretch, no split happens: the edge is struck off (a vertex travels in a
straight line, so it meets that line once and the question is closed for good)
and the next candidate is taken.

### The bound is time-invariant, and the scan exploits it

Write the distance from the vertex to edge `e`'s moving line at time `u` as
`d(u) = d(t0) + closing * (u - t0)`, with `closing` constant along a fixed
trajectory. The crossing time is then

```
    t = u + d(u) / -closing = t0 - d(t0) / closing
```

and the `u` cancels. **The answer does not depend on when it was asked.** So the
scan is done once per trajectory and its result kept in a `SplitCache`, rather
than recomputed each time a neighbour moves and forces a reschedule. Only
`handle_speed_change` can falsify it, by changing an edge's `closing` — and it
invalidates every vertex, not just the ones it moves, because a bound is a race
between a vertex and a *target* edge.

The cache keeps the earliest `SPLIT_FANOUT` (8) candidates, not just one.
Rejections are not the rare case — on star-like input about three quarters of all
split events are rejections, averaging ~3 per reflex vertex — and each one would
otherwise cost a fresh scan to replace the candidate it struck off. Since the
candidates are consumed in increasing time order, the next one is simply the next
in the list.

The scan itself is two passes, because they want opposite things. Computing a
crossing time is the same handful of arithmetic for every edge with no reason to
branch, so the first pass does exactly that over `EdgeLines` — the edges'
normals, offsets and limits transposed into one flat array each — and writes
`INFINITY` where there is no crossing. No `continue`, nothing to trip the
vectoriser. Picking the earliest few is all branching and no arithmetic, so it
gets its own pass over the resulting contiguous `f32`s.

It is still `O(n)` per scan, and it is the only non-constant step left in the
simulation.

## Staleness

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

Split events are not in that set at all, which is the point of asking the weaker
question above: a split's timing is computed from the target edge's supporting
line, and no vertex over there can invalidate it.

### Rescheduling is deduplicated

A single event reaches the same vertex by several routes — as the merged vertex,
as a neighbour's `prev`, and as an explicit mark. Rescheduling it once per route
queues an event per route and immediately strands all but the last. So marks go
through `Sim::mark`, which is idempotent within an event, and every handler
defers scheduling to the drain rather than scheduling inline.

Left undeduplicated this is not a rounding error: it was **77–89% of all popped
events**, including on convex input, which has no split events to blame.

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

### Simultaneous events have no defined resolution

Where several events genuinely coincide — the interior of a comb, whose
identical evenly-spaced teeth make exact ties the norm — which one is processed
first is not defined, and different orderings produce **different but equally
valid** skeletons. The arcs are the same length and meet the same edges; which
node ends up carrying which sources can differ.

So do not diff two skeletons node by node and expect a match. What is guaranteed
is what `tests/common/check_invariants` checks: every node is `offset` from every
edge it names, every arc bisects its pair along its whole length, boundary nodes
have degree 1, and every face closes.

## What it costs

Measured, with `cargo run --release --example bench`, which also times
`Polygon::new` — validation has its own complexity, and leaving it off the clock
would report the crate as faster than any caller can actually get a skeleton.

| input | 1024 vertices | 3200 vertices | scaling |
|---|---|---|---|
| convex | 1.2 ms (at 828) | — (coordinate cap) | ~n^1.2 |
| comb, half reflex | 2.0 ms | 12 ms | ~n^1.3 rising to ~n^1.7 |
| random star, half reflex | 2.9 ms | 16 ms | ~n^1.4 rising to ~n^1.6 |

Space is `O(n)`.

The event count is **linear** — about 3n pops for a comb, 6n for a star — and
each event reschedules `O(1)` vertices. The quadratic term is entirely
`scan_for_split`'s pass over every edge, which is why the exponent climbs with
`n`: the constant is now small enough that the linear work dominates at these
sizes, and the quadratic term only takes over later. It does still take over.

`Polygon::new` is separately `O(n^2)` in the worst case; see
[`DESIGN.md`](DESIGN.md#validation-is-strict-and-up-front).

## Beating the quadratic: the motorcycle graph

The `O(n)` split scan is the only thing standing between this and a
sub-quadratic algorithm, and the known way past it is the **motorcycle graph**.

Every reflex vertex launches a "motorcycle" from its position along its bisector,
at its wavefront speed. Motorcycles leave a trace and crash when they hit a wall
or an earlier trace. The resulting arrangement is the motorcycle graph, and the
theorem that makes it interesting (Eppstein–Erickson; Cheng–Vigneron) is that
**every arc a reflex wavefront vertex traces is part of it**. So the graph
contains the hard part of the skeleton, and given it the rest follows in
`O(n log n)`.

Why it can go faster is worth being precise about, because it is the same reason
this crate cannot. Motorcycle traces are **static rays**. Static things can go in
a spatial index. The moving lines `scan_for_split` searches cannot — a line is
infinite, so no bounding volume excludes it, and the "line" it is racing is at a
different offset every instant. That is the whole trade, and it is why the scan
looks at all `n` edges rather than a local few.

It is not a drop-in, and three things make it a different project rather than an
optimisation:

- **The graph is its own kinetic simulation**, with its own degeneracies —
  simultaneous crashes, and the chicken-and-egg where a motorcycle crashes into
  the trace of one that itself crashes earlier and so never laid that trace.
- **It would not serve `skeleton_constrained`.** The theory is developed for the
  unweighted skeleton; this crate's edges have speeds in `{0, 1}`, and the
  reflex-trace-containment theorem is not something to assume carries over.
- **It does not help the measured bottleneck cheaply.** The tempting shortcut —
  compute the graph, use each motorcycle's crash time to bound its vertex's
  search — does not work. Rejected candidates have crossing times *before* the
  true split, so a bound at the crash excludes none of them. The graph pays off
  only if you adopt the whole construction and stop searching for splits at all.

CGAL and Surfer2 both ship an `O(n^2)` worst case too.

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
- Eppstein, Erickson, *Raising Roofs, Crashing Cycles, and Playing Pool* (1999) —
  the motorcycle graph.
- Cheng, Vigneron, *Motorcycle Graphs and Straight Skeletons* (2002) — the
  sub-quadratic construction.
- Huber, Held, *Theoretical and Practical Results on Straight Skeletons of
  Planar Straight-Line Graphs* (2011) — on the degeneracies that actually matter.
