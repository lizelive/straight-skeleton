//! The **straight skeleton** of a polygon, on the `i16` integer lattice.
//!
//! Shrink a polygon by sliding every edge inward at the same speed, keeping the
//! edges straight and letting them stay connected. The corners trace out a tree
//! of straight line segments. That tree is the straight skeleton, and this
//! crate computes it.
//!
//! ```text
//!    +-----------------------+        +-----------------------+
//!    |                       |        |\                     /|
//!    |                       |        |  \                 /  |
//!    |                       |        |    \_____________/    |
//!    |                       |   ->   |    /             \    |
//!    |                       |        |  /                 \  |
//!    |                       |        |/                     \|
//!    +-----------------------+        +-----------------------+
//!         input polygon                   its straight skeleton
//! ```
//!
//! Straight skeletons are how you find a polygon's medial ridge, generate
//! mitred offset curves, or raise a hip roof over a floor plan (each skeleton
//! node's distance from the boundary *is* the roof height there — see the
//! `roof` example).
//!
//! # Quick start
//!
//! ```
//! use straight_skeleton::{skeleton, Point, Polygon};
//!
//! // A 10x10 square.
//! let square = Polygon::from_outer(&[
//!     Point::new(0, 0),
//!     Point::new(10, 0),
//!     Point::new(10, 10),
//!     Point::new(0, 10),
//! ])?;
//!
//! let skel = skeleton(&square)?;
//!
//! // Its skeleton is an X: the four corners meet at the centre.
//! assert_eq!(skel.arc_count(), 4);
//! let centre = skel.nodes().iter().find(|n| !n.is_boundary()).unwrap();
//! assert_eq!(centre.position, Point::new(5, 5));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Tracing output back to input
//!
//! Every skeleton [`Arc`] separates the faces of **exactly two** input edges,
//! and carries those two ids in [`Arc::sources`]. This is not a
//! nearest-neighbour search bolted on afterwards — it is what an arc *is*, so
//! the lookup is a field access:
//!
//! ```
//! use straight_skeleton::{skeleton, Point, Polygon};
//!
//! let square = Polygon::from_outer(&[
//!     Point::new(0, 0), Point::new(10, 0), Point::new(10, 10), Point::new(0, 10),
//! ])?;
//! let skel = skeleton(&square)?;
//!
//! for arc in skel.arcs() {
//!     let [e0, e1] = arc.sources;
//!     // Every point on this arc is equidistant from the supporting lines of
//!     // input edges e0 and e1.
//!     assert_ne!(e0, e1);
//! }
//!
//! // Nodes carry the same information, with 3+ sources where arcs meet.
//! let centre = skel.nodes().iter().find(|n| !n.is_boundary()).unwrap();
//! assert_eq!(centre.sources.len(), 4); // equidistant from all four edges
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Each input edge also owns one skeleton **[face]**, the region its wavefront
//! swept. The faces tile the polygon, and every face is planar once lifted to
//! `z = offset` — which is exactly why this builds roofs.
//!
//! One caveat, and it is a real one rather than an implementation wrinkle: a
//! straight skeleton is **not** the medial axis. It bisects edges' infinite
//! *supporting lines*, which is what keeps every arc straight; a medial axis
//! bisects the nearest *features*, and grows parabolas around reflex vertices.
//! The two agree exactly when the polygon is convex, and part company around
//! reflex corners. So `sources` means "the edges whose faces meet here", which
//! is the notion you actually want. [`Node::sources`] works through an example.
//!
//! [face]: Skeleton::face
//!
//! # Constrained skeletons
//!
//! [`skeleton_constrained`] caps how far each edge is allowed to travel,
//! **individually**. An edge that hits its limit simply stops, and its
//! neighbours slide along it instead of over it. Use it to truncate a roof to a
//! given eave-to-ridge rise, or to build a variable-width offset.
//!
//! ```
//! use straight_skeleton::{skeleton_constrained, Point, Polygon};
//!
//! let square = Polygon::from_outer(&[
//!     Point::new(0, 0), Point::new(20, 0), Point::new(20, 20), Point::new(0, 20),
//! ])?;
//!
//! // Stop every edge after travelling 3 units, well before the centre at 10.
//! let limits = [3.0; 4];
//! let skel = skeleton_constrained(&square, &limits)?;
//!
//! // Nothing gets further from the boundary than the limit allows.
//! assert!(skel.max_offset() <= 3.0 + 1e-4);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Coordinates, and where the arithmetic lives
//!
//! Input and output coordinates are [`Point`]s of `i16`. The crate prefers
//! narrow types so the algorithm can be ported to a GPU, but two places make a
//! deliberate, documented exception:
//!
//! - **Predicates use `i64`** ([`predicates`]). For `i16` inputs the
//!   orientation determinant needs 35 bits — it does not fit `i32`, and `f32`'s
//!   24-bit mantissa cannot hold it either. `i64` makes every predicate
//!   *exact*: no epsilons, no rounding, no overflow.
//! - **The simulation runs in `f64`.** Skeleton nodes are irrational in
//!   general, so they cannot be computed on the lattice, and `f32` loses too
//!   much precision across the full `i16` range.
//!
//! Node positions are rounded to the lattice at the boundary; [`Node::exact`]
//! keeps the unrounded value as `f32`. `docs/DESIGN.md` works through the
//! width analysis.
//!
//! # Feature flags
//!
//! The crate has **no required dependencies**. Everything below is opt-in.
//!
//! | Feature | Default | Effect |
//! |---|---|---|
//! | `std` | yes | `std::error::Error` impls, hardware `sqrt`. Disable for `no_std`. |
//! | `serde` | no | `Serialize`/`Deserialize` on the public types. |
//! | `geo-types` | no | Conversions to and from `geo_types`. |
//! | `glam` | no | Conversions to and from `glam` vectors. |
//! | `mint` | no | Conversions to and from `mint` vectors. |
//! | `num-traits` | no | Generic numeric conversions. |
//!
//! # `no_std`
//!
//! Disable default features. The crate needs [`alloc`] but nothing else — the
//! only `std` maths it uses is `sqrt`, which it carries its own implementation
//! of rather than take a dependency on `libm`.
//!
//! ```toml
//! straight-skeleton = { version = "0.1", default-features = false }
//! ```
//!
//! [`alloc`]: https://doc.rust-lang.org/alloc/
//! [`Arc::sources`]: crate::Arc::sources
//! [`Node::exact`]: crate::Node::exact

#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]

extern crate alloc;

#[cfg(any(feature = "std", test))]
extern crate std;

mod math;
mod point;
mod polygon;
mod skeleton;
mod wavefront;

pub mod predicates;

#[cfg(any(
    feature = "geo-types",
    feature = "glam",
    feature = "mint",
    feature = "num-traits"
))]
mod interop;

pub use point::Point;
pub use polygon::{EdgeId, Polygon, PolygonError, RingId, VertexId};
pub use skeleton::{Arc, ArcId, Node, NodeId, NodeKind, Skeleton};
pub use wavefront::SkeletonError;

/// Computes the straight skeleton of a polygon.
///
/// Every edge slides inward at unit speed until the polygon has shrunk to
/// nothing. The paths its vertices trace form the skeleton.
///
/// # Errors
///
/// Returns [`SkeletonError`] only for inputs the simulation cannot resolve.
/// A [`Polygon`] is already validated at construction, so for the unconstrained
/// transform this is not expected to fail — if it does, it is a bug.
///
/// # Complexity
///
/// `O(n^2 log n)` time and `O(n)` space, for `n` input vertices. The quadratic
/// term is the split-event search, which tests each reflex vertex against every
/// wavefront edge. Sub-quadratic algorithms exist; this crate ranks correctness
/// and clarity above speed, and the straightforward search is much easier to
/// verify. Convex polygons have no reflex vertices and so run in `O(n log n)`.
///
/// # Examples
///
/// ```
/// use straight_skeleton::{skeleton, Point, Polygon};
///
/// // A triangle's skeleton is three arcs meeting at the incenter.
/// let tri = Polygon::from_outer(&[
///     Point::new(0, 0),
///     Point::new(12, 0),
///     Point::new(0, 9),
/// ])?;
/// let skel = skeleton(&tri)?;
/// assert_eq!(skel.arc_count(), 3);
///
/// // The 9-12-15 triangle has inradius 3, so the incenter sits at (3, 3).
/// let incenter = skel.nodes().iter().find(|n| !n.is_boundary()).unwrap();
/// assert_eq!(incenter.position, Point::new(3, 3));
/// assert!((incenter.offset - 3.0).abs() < 1e-4);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn skeleton(polygon: &Polygon) -> Result<Skeleton, SkeletonError> {
    wavefront::compute(polygon, None)
}

/// Computes a straight skeleton in which each edge stops after travelling a
/// given distance.
///
/// `limits[i]` is the furthest edge `EdgeId(i)` may travel; use
/// [`f32::INFINITY`] to leave an edge unconstrained. Limits are per-edge and
/// need not agree with one another.
///
/// # How it differs from [`skeleton`]
///
/// An edge that reaches its limit stops dead. Its neighbours do not: they keep
/// moving, and the vertices joining them to the stopped edge slide *along* it.
/// So the skeleton's arcs bend at the moment an adjacent edge stops, and the
/// node at that bend is marked [`NodeKind::LimitReached`].
///
/// Passing all-[`f32::INFINITY`] limits reproduces [`skeleton`] exactly: both
/// functions run the same weighted wavefront, and an unlimited edge is one
/// whose speed never drops.
///
/// # Errors
///
/// - [`SkeletonError::LimitCountMismatch`] if `limits.len()` is not
///   [`Polygon::edge_count`].
/// - [`SkeletonError::InvalidLimit`] if a limit is negative or NaN.
/// - [`SkeletonError::IncompatibleCollinearLimits`] if two collinear
///   neighbouring edges are given different limits, which would tear the
///   wavefront apart.
///
/// # Examples
///
/// ```
/// use straight_skeleton::{skeleton, skeleton_constrained, Point, Polygon};
///
/// let square = Polygon::from_outer(&[
///     Point::new(0, 0), Point::new(20, 0), Point::new(20, 20), Point::new(0, 20),
/// ])?;
///
/// // Unlimited on every edge is exactly the plain skeleton.
/// let unlimited = skeleton_constrained(&square, &[f32::INFINITY; 4])?;
/// assert_eq!(unlimited.arcs().len(), skeleton(&square)?.arcs().len());
///
/// // Limits need not be uniform: hold one edge back while the rest advance.
/// let mixed = skeleton_constrained(&square, &[2.0, f32::INFINITY, f32::INFINITY, f32::INFINITY])?;
/// assert!(mixed.node_count() > 0);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn skeleton_constrained(polygon: &Polygon, limits: &[f32]) -> Result<Skeleton, SkeletonError> {
    wavefront::compute(polygon, Some(limits))
}
