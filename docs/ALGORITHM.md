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

**Split** is the expensive one. A reflex vertex is tested against every wavefront
edge: solve for when it reaches that edge's moving line, then check the hit
actually lands **within the live segment**, not merely on its infinite supporting
line. That last check is what separates a real split from a bisector sailing past
outside the edge. It is `O(n)` per reflex vertex, and it is where the quadratic
term comes from.

## Staleness, and the bug it hides

Events are queued, never removed; obsolete ones are recognised when popped. That
needs care, and getting it subtly wrong produced the nastiest bug in this
crate's history.

A split event's timing is computed from **three** vertices: the reflex vertex,
and the two endpoints of the edge it is aimed at. Move any of them and the timing
is worthless. But the two endpoints belong to a completely different part of the
wavefront — so an ordinary edge event over *there* silently invalidates an event
owned by a vertex over *here*.

Reschedule only the vertices you touched, and that vertex is never rescheduled.
It never fires again. Its wavefront sails straight through the boundary, and you
get skeleton nodes sitting outside the polygon.

Nor can it be fixed by checking staleness lazily on pop: the stale event may sit
*later* in the queue than the true one, so by the time it surfaces, the moment to
act is long past.

So `Sim` keeps a reverse index, `dependents[j]` — every vertex whose queued event
was computed from `j`'s geometry. Touching `j` marks all of them dirty, and the
simulation refuses to advance until every dirty vertex has a fresh event. Two
distinct stamps keep the two concerns apart:

- `gen`, bumped when a vertex **moves or is relinked** — invalidates events
  computed *from* it, anywhere in the wavefront.
- `evt`, bumped when a vertex is **rescheduled** — supersedes only its own
  previous event, disturbing nothing else.

Conflating them either loses events or cascades forever.

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

## Complexity

| | time | why |
|---|---|---|
| Convex polygon | `O(n^2)` | no reflex vertices, so no split search; the `n^2` is the input's self-intersection check |
| General | `O(n^2 log n)` typical | each reflex vertex scans the wavefront |
| Constrained | `+ O(n^2)` per distinct limit | every vertex's velocity can bend when an edge stops |

Space is `O(n)`.

Sub-quadratic algorithms exist (Eppstein–Erickson, Cheng–Vigneron). This crate
does not use one. The ordering is **correct > understandable > fast**, and the
straightforward search is dramatically easier to convince yourself of — which,
given how many of the bugs above hid in the *bookkeeping* rather than the
geometry, was the right trade.

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
