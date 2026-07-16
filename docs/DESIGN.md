# Design notes

Why the crate is shaped the way it is. [`ALGORITHM.md`](ALGORITHM.md) covers how
the skeleton is actually computed; this covers the decisions around it.

The stated priority order is **correct > understandable > fast**, and it is not
decoration — it decided most of what follows, and it is the reason to reach for
this crate or not.

## Number types

### The constraint, and where it breaks

The brief asked for `i16` in and out, and for `i32`/`f32` internals so the
algorithm could be ported to a GPU. The first is met exactly. The second is met
almost everywhere, and **deliberately broken in two places**. Both are worth
understanding, because both are forced.

### Predicates are `i64`, and they have to be

An `i16` coordinate spans `-32768..=32767`. A *difference* of two spans
`-65535..=65535` — **17 bits**. The orientation predicate multiplies two
differences and subtracts:

```
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
```

Each product reaches `65535^2 ≈ 2^32`, and the difference of two of them needs
**35 bits**.

- `i32` holds 31 bits plus sign. It **overflows and wraps**, and a wrapped
  determinant does not just lose precision — it reports the *wrong side*. There
  is a test pinning a real triple where `i32` says "left" when the truth is
  "right": `i32_would_overflow_and_report_the_wrong_side`.
- `f32` holds 24 mantissa bits, so it rounds each product by up to `2^8 = 256`.
  Any true determinant below roughly 512 can vanish entirely. There is a test
  pinning a real triple where `f32` reports "collinear" for a genuine turn:
  `f32_would_report_a_real_turn_as_collinear`.

`i64` covers 35 bits with room to spare, so every predicate in `predicates` is
**exact**: no epsilon, no rounding, no overflow, for every `i16` input. Given
that `correct` is first on the list, this was not a close call.

Both tests were written by searching for genuine counterexamples rather than
asserting a plausible-sounding claim — an earlier draft of this file asserted
`f32` would round a full-scale determinant to zero, and that turned out to be
false. The tests now encode what was actually verified.

### The simulation runs in `f64`

Two independent reasons:

1. **Skeleton nodes are irrational.** A 3-4-5 triangle's incenter is rational,
   but rotate it and it is not. There is no lattice to compute *on*; the output
   is rounded to `i16` at the boundary and nowhere else.
2. **`f32` is too coarse across the `i16` range.** At coordinates near 32767,
   `f32` resolves about `0.004`. Event times are computed from divisions of
   accumulated quantities, and errors compound across events. The tolerances
   this crate relies on (`1e-7`, `1e-6`) do not exist in `f32` at that
   magnitude.

`f64` resolves about `1e-11` there — six orders of margin under the merge
tolerance and eleven under one lattice unit.

### So what about the GPU?

Honest answer: **this crate is not GPU-ready, and would not be even in `f32`.**
It is a sequential event simulation over a priority queue and a linked structure
that is rewritten at every event. That is inherently serial; the number types are
not what stands in the way.

What is preserved is the part that transfers: the *public interface* is `i16` and
`f32` throughout (`Point`, `Node::exact`, `Node::offset`, the per-edge limits),
so results feed a GPU pipeline without a widening pass. If you want a parallel
straight skeleton, the literature to start from is motorcycle-graph
constructions, not this.

Calling that out is more useful than quietly shipping an `f32` predicate that is
wrong on a few thousand inputs per million.

### Summary

| Where | Type | Why |
|---|---|---|
| Public input and output | `i16` | as specified |
| `Node::exact`, `Node::offset`, limits | `f32` | narrow, and lossless enough at the boundary |
| Predicates | `i64` | 35 bits needed; exactness is non-negotiable |
| Simulation interior | `f64` | nodes are irrational; `f32` is too coarse at scale |

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

- **Sub-quadratic algorithms.** Priorities said otherwise. See ALGORITHM.md.
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
