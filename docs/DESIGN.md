# Design notes

Why the crate is shaped the way it is. [`ALGORITHM.md`](ALGORITHM.md) covers how
the skeleton is actually computed; this covers the decisions around it.

## The priority order

**correct > fast > understandable.**

It is not decoration — it decided most of what follows, and it is the reason to
reach for this crate or not.

Read the `>` as strict. Where speed and clarity conflict, speed wins and the
clarity is bought back with a comment explaining what the code is doing and why
the obvious version is not there: `scan_for_split` is two passes rather than one
readable loop, the edges are stored transposed as well as as structs, and the
wavefront carries an eight-slot cache that a naive implementation would not need.
None of that is free to read, and all of it is measured.

Where **correctness** and speed conflict, though, it is not a trade — correctness
wins outright, and there is no amount of speed that buys it back:

- `Polygon::check_simple` prunes, but only by tests that are exact (two segments
  disjoint in x cannot cross). It stays `O(n^2)` in the worst case rather than
  adopt a sweep-line algorithm whose degeneracy handling is a known source of
  subtle wrongness. See [below](#validation-is-strict-and-up-front).
- The predicates are exact `i32`, never `f32`, and the coordinate range is capped
  a bit short to keep them so. See [below](#the-cap-and-the-single-expression-that-sets-it).
- Every optimisation in the crate's history has been checked against the
  behaviour it replaced, not just against the test suite. See
  [Testing](#testing).

The order is also why the numbers below are measured rather than derived. A
priority you do not measure is a preference.

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

Honest answer: the *arithmetic* is — `i32` and `f32` only. The *algorithm* is
not, and no choice of number type would make it so. It is a sequential event
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
two things, both small.

`sqrt` is `std`-only, and `no_std` builds normally take it from `libm`.
`math::sqrt_soft` is a Newton–Raphson refinement over the classic
exponent-halving bit trick — about 20 lines. It is **always compiled**, even
under `std`, specifically so it can be differentially tested against the hardware
instruction. The tests sweep a wide range and assert agreement within 1 ULP, so
turning `std` on or off cannot change which branch the algorithm takes.

`f32::floor` is `std`-only too, and `math::floor_i32` covers the one place that
needs it (the node grid's cell lookup). It is four lines and tested against
`f32::floor` over a sweep.

Everything else the algorithm needs is `+ - * /`.

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

**The self-intersection check is where the priority order earns its keep.** The
obvious implementation is a naive all-pairs loop, and the obvious justification
is that validation is far cheaper than the skeleton anyway. It is not: all-pairs
costs 73ms against the skeleton's 13ms on a 3200-vertex comb — five times more
than the thing it feeds, on the critical path of every caller.

So it is a sweep along x, holding open only the edges whose x-range still
overlaps the sweep line. Two segments with disjoint x-ranges cannot cross, so the
skipped pairs are exactly the pairs that could not have failed; a y-overlap test
then drops most of the rest before the exact predicates run. That takes the comb
to 0.13ms.

It is still `O(n^2)` in the worst case, and that is a deliberate stop. A polygon
whose edges all span the full width genuinely has `n^2` pairs to test — a star of
long spokes radiating from a centre only improves from 67ms to 17ms. Beating
*that* needs a real sweep-line intersection algorithm, whose event ordering
around vertical segments, shared endpoints and collinear overlaps is a well-known
source of subtle wrongness. An exact prune that is sometimes no help beats an
asymptotically better algorithm that is sometimes incorrect. `correct > fast`.

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

- A constrained skeleton's **arcs** are disconnected when limits bind. Once every
  edge stops, what the arcs are left as is disjoint stubs reaching in from the
  boundary. That is the point of the transform, not a defect.
- `offset` stops being a distance and becomes the wavefront's **time**. An edge
  that stopped at `limit` stays `limit` away however long the simulation runs, so
  the distance to a source edge `e` is `min(offset, limit_e)`.

### The residual wavefront is half the answer

Those stubs are not the whole result, and treating them as such is the mistake
the API used to invite. A plain skeleton's wavefront shrinks away to nothing —
that is what it means for it to be finished. A constrained one's need not: once
every edge around a loop has stopped, the loop stops too, and stays there. What
it stays as is the input offset inward by the limit, and `Skeleton::residual`
returns it.

That is most of the picture, not a footnote. An L-shape with every wall stopped
at 20 has six stub arcs and a six-sided flat; the flat is the interesting part,
and it looks like the input seen from the inside because that is what it is.

It is **not** made of `Arc`s, and that is the one real design decision here. An
arc bisects the supporting lines of exactly two input edges — that is what makes
`Arc::sources` mean anything, and the whole provenance story rests on it. A
residual segment is *parallel* to one input edge and belongs to it alone. Putting
one in `arcs` would quietly break the invariant every consumer of `sources`
relies on, to save a type. So it gets the type, with the shape it actually has:
each segment names the one edge it came from.

The loops inherit the input's winding, so the interior stays on the left of every
segment — the outer loop counter-clockwise, a loop around a surviving hole
clockwise. A hole's residual *grows*: its wavefront expands into the material, so
an 80x50 hole limited at 10 comes back 100x70 while the outer boundary shrinks.

`Skeleton::face` closes across it, so a constrained skeleton's faces are closed
regions like any other's. Together with the residual they still tile the polygon
exactly — the faces are what the wavefront swept, the residual is what it never
reached — and the test that checks it has to sum the residual's areas *signed*,
because a loop around a surviving hole is a hole in the unswept region rather
than more of it.

## Roof styles are a height function, not four algorithms

`Roof` reads a skeleton off. What it does *not* get from the skeleton is height,
because the skeleton has none: it is the roof's **plan**, saying where the hips,
valleys and ridges run, and that is the same whatever the roof is shaped like.

So the styles are one variable's worth of difference. Height is a function of
`Node::offset` alone, and that function is the `Profile`:

| style | profile | skeleton |
|---|---|---|
| hip | one pitch | plain |
| mansard | two pitches with a break | plain |
| truncated hip | one pitch | uniform limit |
| truncated mansard | two pitches with a break | uniform limit |

A **mansard**'s break sits at a constant offset, which is a constant height, so
it comes out as a level kerb all the way round — which is what a real mansard
has. Offset is affine in position across a face, so the break's level set is a
straight line and cutting there leaves both halves flat. That cut is the only
real work: a panel spanning the break would be *bent*, so each face is split in
two and the corners the cut introduces are the roof's only vertices that stand
over no skeleton node.

Those corners are shared between the two panels either side of every arc, keyed
by the unordered node pair. Minting one per panel would put two vertices at the
same point and leave a crack down every hip.

A **truncated** roof's flat is the residual raised to the limit's height. It
needs no cutting whatever the profile, being level already.

### Only uniform limits have a roof

This is a real constraint rather than a missing feature, and it is worth being
precise about. Height is a function of `offset`, and that only works while offset
means *distance from the wall*. On a plain skeleton it always does. On a
constrained one `offset` is the wavefront's **time**, which is the same thing
only until something stops early — an edge that halted at 3 stays 3 from its face
however long the clock runs on.

With one uniform limit nothing stops early: every edge stops together, at the
top, and the roof is a hip roof truncated to a flat. With uneven limits one
wall's panel would want to end lower than its neighbour's and the surface between
them would have to tear. There is no roof, so `RoofError::UnevenLimits` says so
rather than returning a plausible-looking wrong one.

`Roof` cannot measure distances — it never sees the polygon — so it tests the
condition on the skeleton instead: every `LimitReached` node must sit at
`max_offset`. Uniform limits satisfy that; any edge stopping early does not. It
is conservative by construction, which is the right direction for it to fail in.

Which edge you limit matters in a way that is easy to miss. On a 40x20 rectangle
the ridge sits at offset 10 because the two *long* edges meet there, so limiting a
long edge lowers the ridge — but limiting a *short* edge cannot lower it at all.
What it does instead is **lengthen** it: the short wall's corners stop bisecting
once it freezes and slide straight along it, meeting further out than they
otherwise would.

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

**Faces are walked**, which is the check that pins the *combinatorics* rather
than the geometry. A face only closes if the arcs naming its edge form exactly
one loop, so a node placed on the wrong side of a degenerate tie fails it even
though every individual node is still equidistant from the edges it names.

**Expected geometry is derived by hand**, not read off the implementation: the
9-12-15 triangle's incenter is at `(3,3)` because its inradius is
`(9 + 12 - 15)/2 = 3`; a regular *n*-gon's peak offset is its apothem
`r·cos(π/n)`; a 20x10 rectangle's ridge runs `(5,5)` to `(15,5)`; limiting a
40x20 rectangle's short wall at 3 kinks its corners at `(37,3)` and `(37,17)`.

**Symmetries are tested**, since they catch whole classes of ordering bug:
translation invariance, invariance under which vertex is listed first, invariance
under winding direction.

**The degeneracies have their own tests**, because they are where the real bugs
are: rectangle ridges, simultaneous four-corner collapse, collinear
straight-through vertices, holes placed to pinch strips shut symmetrically,
slivers, coordinates at the `i16` extremes.

**The shapes the benchmark uses are tested too**, at full size. They have
hundreds of reflex vertices, splits and needles, and are by far the most
demanding input the crate sees. A skeleton that is fast and wrong still passes a
benchmark.

**Fast paths are checked against the slow ones they replaced.**
`sweep_agrees_with_all_pairs_on_random_rings` runs `check_simple`'s sweep and the
all-pairs loop it replaced over several thousand random rings on a deliberately
tiny coordinate grid — small enough that crossings, collinear overlaps and shared
endpoints are common rather than rare — and asserts the verdicts match. The
claim being made is about every input, not about the shapes someone thought to
write a test for. `edge_state_and_edge_lines_agree` does the same job for the two
views of an edge's speed limit.

**The examples are tests.** `roof` asserts every panel is genuinely planar, which
holds only if the faces, node positions, and offsets are all correct *together* —
a single misplaced node buckles its panel and trips it.

**There is a harness for whole-skeleton diffs.** `examples/snapshot.rs` dumps
every node and arc of a corpus of shapes with positions as raw bit patterns;
`examples/compare.rs` diffs two such dumps *geometrically*, matching nodes by
their sources and reporting the worst distance it had to bridge. A textual diff
is useless here — the simulation is `f32`, so any change to the order arithmetic
happens in moves results by an ULP — but "worst drift 8e-4 across every shape"
is exactly the statement an optimisation needs to make.

## What is not here

- **A sub-quadratic worst case.** Measured scaling is ~n^1.2 convex, and ~n^1.3
  rising to ~n^1.7 with reflex vertices as `n` grows; the quadratic term is
  `scan_for_split`'s pass over every edge. Beating it needs the **motorcycle
  graph**, which is a different algorithm, would not obviously serve
  `skeleton_constrained`, and cannot be bolted on as a pruner. ALGORITHM.md
  works through why. CGAL and Surfer2 also ship `O(n^2)` worst cases.
- **A defined resolution of simultaneous events.** Where events genuinely
  coincide, different orderings give different but equally valid skeletons, and
  the crate does not promise which. See ALGORITHM.md.
- **Weights other than 0 and 1.** The machinery is general — `velocity_of` solves
  for arbitrary speeds — but arbitrary weights raise degeneracies that are not
  tested, so the API does not expose them.
- **Medial axis.** Different structure, curved arcs. See above.
- **Polygons over 65534 vertices.** `VertexId` is a `u16`, which keeps `Arc` at
  16 bytes. `TooManyVertices` says so.
