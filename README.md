# straight-skeleton

[![CI](https://github.com/lizelive/straight-skeleton/actions/workflows/ci.yml/badge.svg)](https://github.com/lizelive/straight-skeleton/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/straight-skeleton.svg)](https://crates.io/crates/straight-skeleton)
[![docs.rs](https://img.shields.io/docsrs/straight-skeleton)](https://docs.rs/straight-skeleton)
[![licence](https://img.shields.io/crates/l/straight-skeleton.svg)](LICENSE)

The **straight skeleton** of a polygon, with holes, on the `i16` integer
lattice. No required dependencies. `no_std`. No `unsafe`.

Shrink a polygon by sliding every edge inward at the same speed, keeping the
edges straight. The corners trace out a tree of straight segments — that is the
straight skeleton.

```
    +---------------------------+          +---------------------------+
    |                           |          | \                       / |
    |                           |          |   \                   /   |
    |                           |    =>    |     \_______________/     |
    |                           |          |     /               \     |
    |                           |          |   /                   \   |
    +---------------------------+          +---------------------------+
             input                          its straight skeleton
```

Use it to find a polygon's medial ridge, generate mitred offsets, or raise a hip
roof over a floor plan.

## Install

```toml
[dependencies]
straight-skeleton = "0.1"
```

## Use

```rust
use straight_skeleton::{skeleton, Point, Polygon};

let square = Polygon::from_outer(&[
    Point::new(0, 0),
    Point::new(10, 0),
    Point::new(10, 10),
    Point::new(0, 10),
])?;

let skel = skeleton(&square)?;

// A square's skeleton is an X: four corners meeting at the centre.
assert_eq!(skel.arc_count(), 4);
let centre = skel.nodes().iter().find(|n| !n.is_boundary()).unwrap();
assert_eq!(centre.position, Point::new(5, 5));
assert_eq!(centre.offset, 5.0);
# Ok::<(), Box<dyn std::error::Error>>(())
```

## What you get

**Tracing output back to input is a field access.** Every arc separates the faces
of exactly two input edges and carries both ids — not a nearest-neighbour search
bolted on afterwards, but what an arc *is*:

```rust
for arc in skel.arcs() {
    let [e0, e1] = arc.sources;   // the two input edges this arc came from
}
```

**Per-edge distance limits.** `skeleton_constrained` caps how far each edge may
travel, individually. An edge that hits its limit stops; its neighbours slide
*along* it rather than over it:

```rust
use straight_skeleton::skeleton_constrained;

// Stop every edge after 3 units.
let skel = skeleton_constrained(&square, &[3.0; 4])?;
assert!(skel.max_offset() <= 3.0);

// Or hold just one edge back while the rest run free.
let skel = skeleton_constrained(&square, &[3.0, f32::INFINITY, f32::INFINITY, f32::INFINITY])?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

This is not a bolted-on mode: both entry points run the same *weighted*
wavefront, and all-`INFINITY` limits reproduce `skeleton()` exactly.

**Faces, so roofs are nearly free.** Each input edge owns one face; lift its
nodes to `z = offset` and the face is planar. One face, one roof panel.

**Holes, with no special case.** A hole is a wavefront loop that expands instead
of shrinking. Its corners come out reflex, and reflex corners already split the
wavefront, so holes work by the mechanism that was already there.

## Examples

```sh
cargo run --example svg     # draws a gallery to target/svg/index.html
cargo run --example roof    # writes hip roofs as .obj to target/roofs/
```

`svg` renders each shape with its skeleton, colouring arcs by their source edge
so the provenance is visible. `roof` lifts skeletons into 3D and asserts every
panel really is planar.

## Coordinates

`i16` in and out, on the lattice. Two places deliberately use wider arithmetic,
and both are forced:

- **Predicates are `i64`.** For `i16` inputs the orientation determinant needs
  **35 bits**. `i32` overflows and reports the *wrong side*; `f32`'s 24-bit
  mantissa reports a genuine turn as *collinear*. Both failures are pinned by
  tests against real counterexamples. `i64` makes every predicate exact — no
  epsilons, no rounding, no overflow.
- **The simulation is `f64`.** Skeleton nodes are irrational in general, so
  there is no lattice to compute *on*; `f32` resolves only ~0.004 at the far end
  of the `i16` range, which is coarser than the algorithm's tolerances.

Node positions round to the lattice at the boundary, and `Node::exact` keeps the
unrounded `f32`. See [`docs/DESIGN.md`](docs/DESIGN.md) for the width analysis.

## Not the medial axis

Worth knowing before you read the output. Both bisect their input, but a
**straight skeleton** bisects edges' infinite *supporting lines* — which is what
keeps every arc straight — while a **medial axis** bisects the nearest
*features*, growing parabolic arcs around reflex vertices.

They agree exactly when the polygon is convex. Around a reflex corner they do
not: the plus-shape's centre is at offset 5, but its nearest input feature is a
reflex corner 7.07 away. So `offset` is not "distance to the boundary", and
`sources` means "the edges whose faces meet here" rather than "the closest
edges". That is the useful notion anyway — it is what assigns a roof panel to
its wall.

## Features

No required dependencies. Everything is opt-in.

| Feature | Default | Effect |
|---|---|---|
| `std` | yes | `std::error::Error` impls, hardware `sqrt` |
| `serde` | no | `Serialize`/`Deserialize` on the public types |
| `geo-types` | no | conversions to and from `geo_types` |
| `glam` | no | conversions to and from `glam` vectors |
| `mint` | no | conversions to and from `mint` vectors |
| `num-traits` | no | generic numeric conversions |

### `no_std`

```toml
straight-skeleton = { version = "0.1", default-features = false }
```

Needs `alloc`, nothing else. The only `std` maths involved is `sqrt`, and the
crate carries its own rather than depend on `libm` — tested to agree with the
hardware instruction to within 1 ULP, so the feature flag cannot change results.

## Performance

`O(n^2 log n)` typical, `O(n^2)` for convex input, `O(n)` space. The quadratic
term is the split-event search.

Sub-quadratic algorithms exist. This crate does not use one: the priority is
**correct > understandable > fast**, in that order, and most of the real bugs
here hid in event bookkeeping rather than geometry — which is an argument for
code you can check by reading. If you need a straight skeleton of a
hundred-thousand-vertex polygon, this is the wrong crate.

## Documentation

- [API docs](https://docs.rs/straight-skeleton)
- [`docs/ALGORITHM.md`](docs/ALGORITHM.md) — how it works, with diagrams: the
  wavefront, the three events, and the degeneracies that actually bite.
- [`docs/DESIGN.md`](docs/DESIGN.md) — why it is shaped this way: number types,
  API decisions, and what is deliberately missing.

## Licence

GPL-2.0-or-later. See [LICENSE](LICENSE).
