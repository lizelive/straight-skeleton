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

When limits bind hard enough, the wavefront **stops** instead of collapsing —
and what it stops as is the other half of the answer. The arcs are stubs
reaching in from the boundary; `residual()` is the outline they stop on, which
is the input offset inward by the limit. It looks like the outside seen from the
inside, because that is exactly what it is.

```rust
// The 10x10 square, every edge stopped at 3: a 4x4 square is left standing.
let skel = skeleton_constrained(&square, &[3.0; 4])?;
let flat = &skel.residual()[0];

let mut corners: Vec<Point> = flat.nodes.iter().map(|&n| skel.node(n).position).collect();
corners.sort();
assert_eq!(corners, vec![
    Point::new(3, 3), Point::new(3, 7), Point::new(7, 3), Point::new(7, 7),
]);

// A plain skeleton has none: its wavefront always shrinks away to nothing.
assert!(skeleton(&square)?.residual().is_empty());
# Ok::<(), Box<dyn std::error::Error>>(())
```

Each segment names the one wall it came from and runs parallel to it — which is
why it is not an `Arc`: an arc bisects *two* edges, and putting these in `arcs`
would break what `sources` means.

**Roofs, in `i16`, built in.** Each input edge owns one skeleton face; lift its
nodes to `z = offset * pitch` and the face is a flat roof panel. `Roof` does
that for you, and every panel carries the wall it rises from:

```rust
use straight_skeleton::Roof;

let roof = Roof::new(&skel, 0.5)?;          // pitch: rise over run
assert_eq!(roof.panels().len(), 4);          // one panel per wall
let apex: i16 = roof.ridge_height();         // i16, like everything else

for panel in roof.panels() {
    let wall = panel.wall();                 // which wall this panel rises from
}
# Ok::<(), Box<dyn std::error::Error>>(())
```

Panels are planar by construction, not by fitting: every point on one is
`offset` from the same wall's line, and height is linear in that distance.

**Mansard and truncated roofs, off the same skeleton.** A skeleton is a roof's
*plan* — where the hips, valleys and ridges run — and that does not depend on how
high anything is. So the style is one variable: the `Profile`, which says how
height grows with distance from the eaves.

```rust
use straight_skeleton::{skeleton_constrained, Roof};

// Steep to a kerb at 10, then shallow. Two panels per wall, split along a
// level break line; the corners that split introduces are shared between
// neighbours, so the mesh stays watertight.
let mansard = Roof::mansard(&skel, 2.0, 10.0, 0.25)?;

// Stop every wall at 15 and the apex is cut off, leaving a flat on top —
// which is the constrained skeleton's residual, raised.
let truncated = Roof::new(&skeleton_constrained(&square, &[15.0; 4])?, 0.5)?;
assert_eq!(truncated.flat().count(), 1);

// The two compose: steep, then shallow, then flat.
let both = Roof::mansard(&skeleton_constrained(&square, &[15.0; 4])?, 2.0, 10.0, 0.25)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

A limit of **zero** does the opposite job, and gives you gables. The wall never
moves, so it sweeps nothing: its panel is the degenerate face stood on end — a
vertical gable — and the ridge runs out to it rather than hipping away. A mansard
with gable ends is a *gambrel*, the barn roof:

```rust
// A 240x90 hall with both short ends frozen. The ridge now runs the whole 240.
let gabled = skeleton_constrained(&hall, &[f32::INFINITY, 0.0, f32::INFINITY, 0.0])?;
let gambrel = Roof::mansard(&gabled, 2.0, 10.0, 0.25)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

`cargo run --example roof` writes every style of every plan as OBJ.

Limits must be uniform **or zero** for a roof, and that is geometry rather than a
gap. An edge stopping partway still has a sloping panel, which would end lower
than its neighbour's and tear the surface between them — so
`RoofError::UnevenLimits` says so rather than returning a plausible wrong answer.
A zero-limit wall has no sloping panel to be inconsistent about, which is exactly
why gables are allowed where half-measures are not.

**Holes, with no special case.** A hole is a wavefront loop that expands instead
of shrinking. Its corners come out reflex, and reflex corners already split the
wavefront, so holes work by the mechanism that was already there.

## Examples

```sh
cargo run --example svg     # draws a gallery to target/svg/index.html
cargo run --example roof    # writes roofs as .obj to target/roofs/
```

`svg` renders each shape with its skeleton, colouring arcs by their source edge
so the provenance is visible, and shading the residual where a constrained shape
has one. `roof` writes every plan four ways — hip, mansard, truncated, truncated
mansard — off one skeleton each; it uses the library's `Roof` type and adds only
the OBJ serialisation and the triangulation that goes with it, since choosing a
file format and a mesh representation are the caller's business rather than the
crate's.

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
growth exponent rather than a claim, and times `Polygon::new` alongside the
skeleton — validation has its own complexity, and leaving it off the clock would
report the crate as faster than any caller can actually get a skeleton.

| input | 1024 vertices | 3200 vertices | scaling |
|---|---|---|---|
| convex | 1.2 ms (at 828) | — (hits the coordinate cap) | ~n^1.2 |
| comb (half reflex) | 2.0 ms | 12 ms | ~n^1.3, rising to ~n^1.7 |
| random star (half reflex) | 2.9 ms | 16 ms | ~n^1.4, rising to ~n^1.6 |

The event count is linear and each event reschedules `O(1)` vertices. The one
non-constant step is the split search, which scans every edge for each reflex
vertex — so the worst case is `O(n^2)`, and the rising exponent is that term
gradually taking over as `n` grows. **Convex input never runs it at all**, having
no reflex vertices to search from.

Beating `O(n^2)` needs the motorcycle graph. It is a genuinely different
algorithm rather than an optimisation, it would not obviously serve
`skeleton_constrained`, and — for the reason `docs/ALGORITHM.md` works through —
it cannot be bolted on as a mere pruner. CGAL and Surfer2 ship an `O(n^2)` worst
case too.

Space is `O(n)`.

## Documentation

- [API docs](https://docs.rs/straight-skeleton)
- [`docs/ALGORITHM.md`](docs/ALGORITHM.md) — how it works, with diagrams: the
  wavefront, the three events, and the degeneracies that actually bite.
- [`docs/DESIGN.md`](docs/DESIGN.md) — why it is shaped this way: number types,
  API decisions, and what is deliberately missing.

## Licence

GPL-2.0-or-later. See [LICENSE](LICENSE).
