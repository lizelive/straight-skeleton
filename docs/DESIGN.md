# Design notes

Why the crate is shaped the way it is. [`ALGORITHM.md`](ALGORITHM.md) covers how
the skeleton is actually computed; this covers the decisions around it.

The stated priority order is **correct > understandable > fast**, and it is not
decoration — it decided most of what follows, and it is the reason to reach for
this crate or not.

## Number types

### The rule

**`i16` in and out. `i32` and `f32` in between. Nothing wider, anywhere.**

No `f64`, and no `i64` either, in any part of the algorithm. The point is
portability to hardware where `f64` is slow or missing, and a type you only use
"internally" is still a type the hardware has to have. The one exception is
`ring_area2`, which is validation — see below.

That rule is not free. It is paid for with exactly one bit of coordinate range.

### The cap, and the single expression that sets it

Coordinates are capped at `-16384..=16383` (`Point::MIN_COORD`/`MAX_COORD`),
one bit narrower than `i16`. `Polygon` rejects anything outside.

Everything traces back to the orientation determinant:

```
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
```

Its width is set by the largest coordinate difference `d`, as `2 * d^2`:

| coordinates | largest `d` | `2 * d^2` | in `i32`? |
|---|---|---|---|
| full `i16` | 65_535 | 8_589_672_450 | **overflows** — and wraps to the *wrong side* |
| capped | 32_767 | 2_147_352_578 | fits, with 131_069 to spare |

So one bit is exactly the price of an **exact** `i32` predicate: no epsilon, no
rounding, no overflow, for every input `Polygon` accepts. And it is exactly one
bit — the capped worst case uses 99.994% of `i32`. Two bits would be waste; zero
bits silently reports the wrong side.

### The counterexamples are real, not rhetorical

Both failure modes are pinned by tests, against triples found by **search**
rather than asserted:

- `beyond_the_cap_i32_would_report_the_wrong_side` — at `(21203,-24650)`,
  `(-22519,1049)`, `(26449,26335)`, the true determinant is `-2_363_983_124` and
  `i32` wraps to `+1_930_984_172`. Sign flipped: it reports *left* where the
  truth is *right*.
- `f32_would_report_a_real_turn_as_collinear_even_inside_the_cap` — at
  `(14176,-12146)`, `(-9937,5341)`, `(4434,-5081)`, all **within** the cap, the
  true determinant is `9` and `f32` returns `0`. A genuine turn, reported as
  collinear.

The second is why the cap does not rescue `f32` for predicates. `f32` has 24
mantissa bits; these products need up to 31. No range short of absurd fixes
that, so the predicate is `i32` and the simulation is `f32`, and they are
different tools for different jobs.

(An earlier draft of this file asserted that `f32` rounds a full-scale
determinant to zero. That turned out to be **false** when checked — it returns
65536 against a true 65535: wrong, but not zero. The tests now encode what was
actually verified. Search, don't assert.)

### The simulation is `f32` — and what that costs

Skeleton nodes are irrational in general (rotate a 3-4-5 triangle and its
incenter stops being rational), so there is no lattice to compute *on*. The
simulation is `f32`, and the cap is what leaves it enough room: `0.002` of
absolute resolution at the worst corner of the coordinate space, against `0.004`
uncapped.

**This is a real robustness trade, and worth being straight about.** `f64` at
the same coordinates resolves `~1e-11`. Two skeleton features on an integer
lattice of size `R` can genuinely be `~1/R^2` apart, which at `R = 16384` is
`~4e-9` — far below what `f32` can see. So:

- For ordinary input, `f32` is fine, and the whole test suite passes at the cap.
- For adversarially near-degenerate input, `f32` cannot distinguish what `f64`
  could. The tolerances (`EPS = 1e-4`, `MERGE_EPS = 1e-2`) are set to absorb
  `f32` noise, and two features closer than `MERGE_EPS` will be fused.

What makes that survivable is that robustness here is mostly *structural*, not
numeric: needle zipping and node deduplication are what handle degeneracy, and
both are tolerance-based by design rather than precision-based.

If you need worst-case robustness at full `i16` range, the honest answer is
double-float arithmetic — a pair of `f32`s giving ~48 mantissa bits, the standard
technique on hardware without `f64`. That would restore the margin at roughly
10-20x the arithmetic cost, and it is the obvious next step if the cap ever bites.

### The one exception: `ring_area2`

`ring_area2` keeps an `i64` accumulator. Unlike `orient2d` it *sums* triangles,
and a ring that doubles back can wind around a region more than once, so the
running total is bounded by the vertex count rather than the coordinate box. One
bit of range cannot fix that; only a wider accumulator can.

It is a fair exception because it is **validation**: it runs once, on the host,
when a `Polygon` is built, and never during the simulation. The skeleton itself
is `i32` and `f32` throughout.

### Summary

| Where | Type | Why |
|---|---|---|
| Public input and output | `i16`, capped to `±16384` | one bit buys exact `i32` predicates |
| Predicates | `i32` | exact within the cap; `f32` gets it wrong even inside it |
| Simulation interior | `f32` | nodes are irrational; the cap leaves enough resolution |
| `Node::exact`, offsets, limits, roof heights | `f32` | what the simulation actually produced |
| `ring_area2` only | `i64` | sums can wind; validation-only, never in the simulation |

### So is it GPU-ready?

Honest answer: the *arithmetic* now is — `i32` and `f32` only. The *algorithm*
is not, and no choice of number type would make it so. It is a sequential event
simulation over a priority queue and a linked structure rewritten at every event.
That is inherently serial. A parallel straight skeleton is a motorcycle-graph
construction, not this.

What the rule buys is real regardless: the whole thing runs on hardware without
`f64`, results feed a GPU pipeline with no widening pass, and the predicate is
exact rather than approximately-usually-right.

## `i16` output, and rounding

`Node::position` is the nearest lattice point, rounding half away from zero and
**saturating** rather than wrapping at the `i16` bounds. Saturation is a
safety net, not a path anything should take: a straight skeleton lies within its
input's convex hull, so a node can only exceed the coordinate range through
floating-point error at the very edge of the space.

`Node::exact` carries the unrounded `f32` alongside. Two distinct nodes can round
to the same lattice point on a small enough polygon; when that matters, use
`exact`.

## No required dependencies

The crate compiles with zero dependencies, `std` or `no_std`. That cost exactly
one thing: `sqrt`, which is `std`-only and which `no_std` builds normally take
from `libm`.

`math::sqrt_soft` is a Newton–Raphson refinement over the classic
exponent-halving bit trick — about 20 lines. It is **always compiled**, even
under `std`, specifically so it can be differentially tested against the hardware
instruction. The tests sweep a wide range and assert agreement within 1 ULP, so
turning `std` on or off cannot change which branch the algorithm takes.

That is the only transcendental the algorithm needs; everything else is `+ - * /`.

`alloc` is required. The event queue and wavefront arena grow with the input and
there is no sensible fixed bound.

## The API

### Flattened rings, and free provenance

`Polygon` stores all rings in one `Vec<Point>` with ring boundaries alongside.
This buys the identity that carries the whole traceability story:

> **Edge `i` starts at vertex `i`.**

So `EdgeId` and `VertexId` are the same number, converting between an edge and
its start vertex is free, and the boundary nodes come out in vertex order —
which is what makes `Skeleton::boundary_node` an index rather than a search.

### Interior on the left, always

Outer rings are normalised counter-clockwise and holes clockwise, so the
polygon's interior is to the **left of every directed edge without exception**.
Rings supplied the other way round are reversed for you.

This is the invariant that lets the wavefront treat holes and the outer boundary
identically. There is no `if is_hole` anywhere in the simulation. A hole's
corners simply come out reflex, and reflex vertices split — so holes work by the
mechanism that already existed.

### Validation is strict, and up front

`Polygon::new` rejects rather than repairs: too few vertices, repeated vertices,
zero area, crossing edges, 180° spikes, holes outside the outer ring. Every error
names the offending ring, and where meaningful the vertex or edge.

This is deliberate. A straight skeleton of an invalid polygon is not
ill-conditioned, it is *undefined* — and quietly returning a plausible-looking
graph for a self-intersecting input is worse than an error. Validating once at
construction is also what lets the simulation carry no defensive checks: by the
time `skeleton()` runs, zero-length edges and spikes cannot exist.

The one repair performed is winding normalisation, because it is unambiguous and
every caller wants it.

The self-intersection check is naive all-pairs `O(n^2)`. It is comfortably
cheaper than the skeleton itself, and keeping it obvious is worth more than the
constant factor.

### `sources`, not `closest`

The brief asked that it be trivial to find the input features nearest an output
feature. What the crate provides is **exact and free** — `Arc::sources` is two
`EdgeId`s sitting in the struct — but it is worth being precise about what it
means, because the obvious reading is subtly wrong.

A straight skeleton is **not** the medial axis. Its arcs bisect edges'
*supporting lines*; a medial axis bisects the nearest *features*. They coincide
exactly when the polygon is convex, and diverge around reflex corners. The
plus-shape's centre is at offset 5, but its nearest input feature is a reflex
corner 7.07 away.

So `sources` means **"the input edges whose faces meet here"**. That is what is
true, and it is the notion that is actually wanted: it is what assigns a roof
panel to its wall. If you need genuinely-nearest features on a non-convex
polygon, you want a medial axis, which is a different structure with curved arcs.

This distinction is documented on `Node::sources`, `Arc`, and the crate root,
because it is the thing most likely to surprise someone reading the output.

### Constrained: one mechanism, not two

`skeleton_constrained` does not branch off a separate code path. Both entry
points run the same **weighted** wavefront; a per-edge limit is just "this edge's
speed drops to 0 at `t = limit`". Passing all-`INFINITY` limits reproduces
`skeleton()` **exactly**, which is asserted rather than assumed
(`infinite_limits_reproduce_the_plain_skeleton` compares nodes and arcs for
equality).

Two consequences fall out of the semantics and are documented on the API:

- A constrained skeleton is **disconnected** when limits bind. Once every edge
  stops, what is left is disjoint stubs reaching in from the boundary. That is
  the point of the transform, not a defect.
- `offset` stops being a distance and becomes the wavefront's **time**. An edge
  that stopped at `limit` stays `limit` away however long the simulation runs, so
  the distance to a source edge `e` is `min(offset, limit_e)`.

One configuration is refused: two **collinear neighbouring** edges given
different limits. One line stops while the other, parallel to it, keeps going,
and the vertex between them has nowhere to be — the wavefront would have to tear
open. `IncompatibleCollinearLimits` says so rather than inventing a shape.

## Testing

The suite is built so it can actually fail when the implementation is wrong.

**Invariants are re-derived from the definition**, in `tests/common`, taking the
polygon and the skeleton and checking what must be true of any straight skeleton
— sources are equidistant from their supporting lines, arcs bisect their pair at
the *midpoint* (endpoints would pass trivially), boundary nodes have degree
exactly 1, the graph is connected, offsets never exceed the boundary distance.

**Expected geometry is derived by hand**, not read off the implementation: the
9-12-15 triangle's incenter is at `(3,3)` because its inradius is
`(9 + 12 - 15)/2 = 3`; a regular *n*-gon's peak offset is its apothem
`r·cos(π/n)`; a 20x10 rectangle's ridge runs `(5,5)` to `(15,5)`.

**Symmetries are tested**, since they catch whole classes of ordering bug:
translation invariance, invariance under which vertex is listed first, invariance
under winding direction. The starting-vertex test is what caught a duplicate-node
bug that every other test missed.

**The degeneracies have their own tests**, because they are where the real bugs
were: rectangle ridges, simultaneous four-corner collapse, collinear
straight-through vertices, holes placed to pinch strips shut symmetrically,
slivers, coordinates at the `i16` extremes.

**The examples are tests.** `roof` asserts every panel is genuinely planar, which
holds only if the faces, node positions, and offsets are all correct *together* —
a single misplaced node buckles its panel and trips it.

## What is not here

- **A sub-quadratic worst case.** Measured scaling is ~n^1.3 convex and ~n^1.9
  to n^2.1 with reflex vertices; the quadratic term is `split_lower_bound`'s scan
  over every edge. Beating it needs the **motorcycle graph**: reflex vertices
  launch motorcycles along their bisectors, split events are exactly where they
  crash, and motorcycle traces are *static rays* — so they can go in a spatial
  index, which moving wavefront edges cannot. Given the graph, the skeleton
  follows in `O(n log n)` (Cheng–Vigneron, Huber–Held). It is a different
  algorithm and a much larger one; CGAL and Surfer2 also ship `O(n^2)` worst
  cases. See ALGORITHM.md.
- **Weights other than 0 and 1.** The machinery is general — `velocity_of` solves
  for arbitrary speeds — but arbitrary weights raise degeneracies that are not
  tested, so the API does not expose them.
- **The residual wavefront.** A constrained skeleton leaves an offset polygon
  behind; it is not returned, since a wavefront edge is *parallel* to an input
  edge rather than bisecting two, and putting it in `arcs` would break what an
  `Arc` means. Worth adding as its own type.
- **Medial axis.** Different structure, curved arcs. See above.
- **Polygons over 65534 vertices.** `VertexId` is a `u16`, which keeps `Arc` at
  16 bytes. `TooManyVertices` says so.
