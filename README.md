# straight-skeleton

[![CI](https://github.com/lizelive/straight-skeleton/actions/workflows/ci.yml/badge.svg)](https://github.com/lizelive/straight-skeleton/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/straight-skeleton.svg)](https://crates.io/crates/straight-skeleton)
[![docs.rs](https://img.shields.io/docsrs/straight-skeleton)](https://docs.rs/straight-skeleton)
[![licence](https://img.shields.io/crates/l/straight-skeleton.svg)](LICENSE)

The **straight skeleton** of a polygon, with holes, on the `i16` integer
lattice. Computed entirely in `i32` and `f32` — **no `f64`, no `i64`**. No
required dependencies. `no_std`. No `unsafe`.

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

**Roofs, in `i16`, built in.** Each input edge owns one skeleton face; lift its
nodes to `z = offset * pitch` and the face is a flat roof panel. `Roof` does
that for you, and every panel carries the wall it rises from:

```rust
use straight_skeleton::Roof;

let roof = Roof::new(&skel, 0.5)?;          // pitch: rise over run
assert_eq!(roof.panels().len(), 4);          // one panel per wall
let apex: i16 = roof.ridge_height();         // i16, like everything else

for panel in roof.panels() {
    let wall = panel.wall;                   // which wall this panel rises from
}
# Ok::<(), Box<dyn std::error::Error>>(())
```

Panels are planar by construction, not by fitting: every point on one is
`offset` from the same wall's line, and height is linear in that distance.

**Holes, with no special case.** A hole is a wavefront loop that expands instead
of shrinking. Its corners come out reflex, and reflex corners already split the
wavefront, so holes work by the mechanism that was already there.

## Examples

```sh
cargo run --example svg     # draws a gallery to target/svg/index.html
cargo run --example roof    # writes hip roofs as .obj to target/roofs/
```

`svg` renders each shape with its skeleton, colouring arcs by their source edge
so the provenance is visible. `roof` uses the library's `Roof` type and only
adds the OBJ serialisation, since choosing a file format is the caller's
business rather than the crate's.

## Coordinates: `i32` and `f32`, no `f64`

`i16` in and out. Everything in between is `i32` and `f32` — nothing wider, so
the arithmetic ports to hardware where `f64` is slow or absent.

That costs **one bit of range**: coordinates are capped at `-16384..=16383`
(`Point::MIN_COORD`/`MAX_COORD`), and `Polygon` rejects anything outside. One
expression sets the cap — the orientation determinant, which needs `2 * d^2` for
the largest coordinate difference `d`:

| coordinates | `2 * d^2` | in `i32`? |
|---|---|---|
| full `i16` | 8,589,672,450 | **overflows**, reporting the *wrong side* |
| capped | 2,147,352,578 | fits, with 131,069 to spare |

So one bit buys **exact** predicates: no epsilon, no rounding, no overflow.
`f32` cannot do that job at any range — the tests pin a real triple *inside* the
cap where it calls a genuine turn collinear, and one outside where `i32` flips
sign. Both were found by search, not asserted.

The simulation is `f32`. Skeleton nodes are irrational in general, so there is
no lattice to compute *on*; positions round back to it at the boundary and
`Node::exact` keeps the unrounded value. The cap is also what leaves `f32`
enough absolute resolution — ~0.002 at worst — to work in.

[`docs/DESIGN.md`](docs/DESIGN.md) has the full analysis, including what this
costs in robustness.

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

Measured with `cargo run --release --example bench`, which reports the empirical
growth exponent rather than a claim:

| input | 1024 vertices | scaling |
|---|---|---|
| convex | 1.1 ms | ~n^1.3 |
| comb (512 reflex) | 8.9 ms | ~n^1.9 |
| random star (512 reflex) | 13 ms | ~n^2.1 |

**Convex input is sub-quadratic** — it has no reflex vertices, so it never
searches for split events at all. With reflex vertices the search is `O(n)` per
reflex vertex, giving `O(n^2)` overall, and the measurements agree.

Beating `O(n^2)` in the worst case needs the motorcycle-graph machinery of
Eppstein–Erickson or Cheng–Vigneron. That is a different algorithm and a much
larger one; practical implementations (CGAL, Surfer2) also ship an `O(n^2)`
worst case. `docs/DESIGN.md` explains what it would take and why it is not here.

Space is `O(n)`.

## Documentation

- [API docs](https://docs.rs/straight-skeleton)
- [`docs/ALGORITHM.md`](docs/ALGORITHM.md) — how it works, with diagrams: the
  wavefront, the three events, and the degeneracies that actually bite.
- [`docs/DESIGN.md`](docs/DESIGN.md) — why it is shaped this way: number types,
  API decisions, and what is deliberately missing.

## Licence

GPL-2.0-or-later. See [LICENSE](LICENSE).
