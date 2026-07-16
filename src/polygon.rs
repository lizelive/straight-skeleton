//! Input polygons, their identifiers, and validation.

use alloc::vec::Vec;
use core::fmt;

use crate::predicates::{is_ccw, orient2d, ring_area2, segments_properly_cross, Orientation};
use crate::Point;

/// Identifies an input vertex of a [`Polygon`].
///
/// Vertices are numbered across the whole polygon, outer ring first, then each
/// hole in order. See [`Polygon`] for the numbering guarantee.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct VertexId(pub u16);

/// Identifies an input edge of a [`Polygon`].
///
/// # The `EdgeId` / `VertexId` correspondence
///
/// Edge `i` is the edge that **starts** at vertex `i` and ends at the next
/// vertex of the same ring (wrapping at the ring's end). So `EdgeId(i)` and
/// `VertexId(i)` always share a number, and converting between an edge and its
/// start vertex is free. This is the whole reason the crate stores rings
/// flattened.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EdgeId(pub u16);

impl EdgeId {
    /// The vertex this edge starts at.
    #[inline]
    pub const fn start_vertex(self) -> VertexId {
        VertexId(self.0)
    }
}

impl VertexId {
    /// The edge that starts at this vertex.
    #[inline]
    pub const fn outgoing_edge(self) -> EdgeId {
        EdgeId(self.0)
    }
}

impl fmt::Display for VertexId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "v{}", self.0)
    }
}

impl fmt::Display for EdgeId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "e{}", self.0)
    }
}

/// Identifies a ring of a [`Polygon`]. Ring 0 is always the outer boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RingId(pub u16);

/// Why a [`Polygon`] could not be built.
///
/// Every variant names the ring (and where meaningful the vertex or edge)
/// responsible, so the caller can point at the offending input rather than
/// guessing.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum PolygonError {
    /// A polygon needs an outer ring; none was supplied.
    NoOuterRing,
    /// A ring had fewer than three distinct vertices.
    TooFewVertices {
        /// The offending ring.
        ring: RingId,
        /// How many vertices it had.
        count: usize,
    },
    /// A ring repeated a vertex back to back, giving a zero-length edge.
    RepeatedVertex {
        /// The offending ring.
        ring: RingId,
        /// Index of the repeat within the ring.
        index: usize,
        /// The duplicated point.
        point: Point,
    },
    /// A ring encloses no area (every vertex is collinear).
    DegenerateRing {
        /// The offending ring.
        ring: RingId,
    },
    /// A ring doubles back on itself through 180°, forming a zero-width spike.
    ///
    /// The wavefront vertex at such a corner would have to move infinitely
    /// fast, so the skeleton is undefined there. Nudge the spike tip sideways
    /// by one unit, or drop it.
    Spike {
        /// The offending ring.
        ring: RingId,
        /// The spike tip.
        vertex: VertexId,
    },
    /// Two edges of the polygon cross. Simple polygons only.
    SelfIntersection {
        /// One of the crossing edges.
        a: EdgeId,
        /// The other crossing edge.
        b: EdgeId,
    },
    /// A hole is not contained in the outer ring.
    HoleOutsideOuter {
        /// The offending hole.
        ring: RingId,
    },
    /// The polygon has more vertices than [`EdgeId`] can number.
    TooManyVertices {
        /// How many vertices were supplied.
        count: usize,
        /// The most that can be numbered.
        max: usize,
    },
}

impl fmt::Display for PolygonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PolygonError::NoOuterRing => write!(f, "polygon has no outer ring"),
            PolygonError::TooFewVertices { ring, count } => write!(
                f,
                "ring {} has {count} vertices; at least 3 are required",
                ring.0
            ),
            PolygonError::RepeatedVertex { ring, index, point } => write!(
                f,
                "ring {} repeats vertex ({}, {}) at index {index}, giving a zero-length edge",
                ring.0, point.x, point.y
            ),
            PolygonError::DegenerateRing { ring } => {
                write!(f, "ring {} encloses no area", ring.0)
            }
            PolygonError::Spike { ring, vertex } => write!(
                f,
                "ring {} doubles back through 180° at vertex {}, forming a zero-width spike",
                ring.0, vertex.0
            ),
            PolygonError::SelfIntersection { a, b } => {
                write!(
                    f,
                    "edges {} and {} cross; the polygon must be simple",
                    a.0, b.0
                )
            }
            PolygonError::HoleOutsideOuter { ring } => {
                write!(f, "hole {} is not contained in the outer ring", ring.0)
            }
            PolygonError::TooManyVertices { count, max } => {
                write!(f, "polygon has {count} vertices; the maximum is {max}")
            }
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for PolygonError {}

/// A simple polygon, optionally with holes, on the `i16` lattice.
///
/// # Invariants
///
/// A `Polygon` can only be constructed through [`Polygon::new`] or
/// [`Polygon::from_outer`], which enforce that:
///
/// - Ring 0 is the outer boundary; rings `1..` are holes.
/// - The outer ring winds **counter-clockwise** and holes wind **clockwise**,
///   so the polygon's interior is always on the *left* of every directed edge.
///   Rings supplied the other way round are reversed automatically.
/// - Every ring has at least 3 vertices, no repeated consecutive vertices, and
///   encloses a non-zero area.
/// - No two edges cross, and no vertex is a zero-width spike.
/// - Every hole lies inside the outer ring.
///
/// The uniform "interior on the left" rule is what lets the wavefront treat
/// outer boundary and holes identically — see `docs/ALGORITHM.md`.
///
/// # Vertex and edge numbering
///
/// Vertices are numbered `0..n` across all rings, outer ring first. Edge `i`
/// starts at vertex `i`; see [`EdgeId`].
///
/// # Examples
///
/// ```
/// use straight_skeleton::{Point, Polygon};
///
/// // A square. Winding is fixed up for you.
/// let square = Polygon::from_outer(&[
///     Point::new(0, 0),
///     Point::new(10, 0),
///     Point::new(10, 10),
///     Point::new(0, 10),
/// ])?;
/// assert_eq!(square.vertex_count(), 4);
/// assert_eq!(square.ring_count(), 1);
///
/// // A square with a square hole.
/// let with_hole = Polygon::new(
///     &[Point::new(0, 0), Point::new(30, 0), Point::new(30, 30), Point::new(0, 30)],
///     &[vec![
///         Point::new(10, 10),
///         Point::new(20, 10),
///         Point::new(20, 20),
///         Point::new(10, 20),
///     ]],
/// )?;
/// assert_eq!(with_hole.ring_count(), 2);
/// assert_eq!(with_hole.vertex_count(), 8);
/// # Ok::<(), straight_skeleton::PolygonError>(())
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Polygon {
    /// All vertices, outer ring first, then holes in order.
    verts: Vec<Point>,
    /// `ring_starts[i]..ring_starts[i + 1]` is ring `i`'s slice of `verts`.
    /// Has `ring_count + 1` entries; the last is `verts.len()`.
    ring_starts: Vec<u16>,
}

impl Polygon {
    /// The largest number of vertices a polygon may have.
    ///
    /// Bounded by [`VertexId`]'s `u16`, minus one so that "one past the end"
    /// indices cannot overflow.
    pub const MAX_VERTICES: usize = u16::MAX as usize - 1;

    /// Builds a polygon from an outer ring and a list of holes.
    ///
    /// Ring winding is normalised for you: the outer ring is made
    /// counter-clockwise and holes clockwise, reversing any ring given the
    /// other way round.
    ///
    /// # Errors
    ///
    /// Returns a [`PolygonError`] naming the offending ring if the input is not
    /// a simple polygon with holes. See [`Polygon`]'s invariants for the full
    /// list of checks.
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::{Point, Polygon, PolygonError};
    ///
    /// // A ring that crosses itself is rejected, not silently accepted.
    /// let crossed = Polygon::from_outer(&[
    ///     Point::new(0, 0),
    ///     Point::new(10, 10),
    ///     Point::new(10, 0),
    ///     Point::new(0, 4),
    /// ]);
    /// assert!(matches!(crossed, Err(PolygonError::SelfIntersection { .. })));
    /// ```
    pub fn new(outer: &[Point], holes: &[Vec<Point>]) -> Result<Self, PolygonError> {
        if outer.is_empty() {
            return Err(PolygonError::NoOuterRing);
        }

        let total: usize = outer.len() + holes.iter().map(|h| h.len()).sum::<usize>();
        if total > Self::MAX_VERTICES {
            return Err(PolygonError::TooManyVertices {
                count: total,
                max: Self::MAX_VERTICES,
            });
        }

        let mut verts: Vec<Point> = Vec::with_capacity(total);
        let mut ring_starts: Vec<u16> = Vec::with_capacity(holes.len() + 2);
        ring_starts.push(0);

        // The outer ring must be CCW and holes CW, so that the interior lies to
        // the left of every directed edge without exception.
        push_ring(&mut verts, &mut ring_starts, outer, RingId(0), true)?;
        for (i, hole) in holes.iter().enumerate() {
            let id = RingId((i + 1) as u16);
            push_ring(&mut verts, &mut ring_starts, hole, id, false)?;
        }

        let poly = Polygon { verts, ring_starts };
        poly.check_no_spikes()?;
        poly.check_simple()?;
        poly.check_holes_inside()?;
        Ok(poly)
    }

    /// Builds a polygon with no holes.
    ///
    /// # Errors
    ///
    /// As [`Polygon::new`].
    pub fn from_outer(outer: &[Point]) -> Result<Self, PolygonError> {
        Polygon::new(outer, &[])
    }

    /// Total number of vertices across all rings.
    #[inline]
    pub fn vertex_count(&self) -> usize {
        self.verts.len()
    }

    /// Number of rings: 1 (outer) plus one per hole.
    #[inline]
    pub fn ring_count(&self) -> usize {
        self.ring_starts.len() - 1
    }

    /// Number of holes.
    #[inline]
    pub fn hole_count(&self) -> usize {
        self.ring_count() - 1
    }

    /// All vertices, outer ring first, then holes in order.
    #[inline]
    pub fn vertices(&self) -> &[Point] {
        &self.verts
    }

    /// The position of a vertex.
    ///
    /// # Panics
    ///
    /// Panics if `v` does not belong to this polygon.
    #[inline]
    pub fn vertex(&self, v: VertexId) -> Point {
        self.verts[v.0 as usize]
    }

    /// The vertices of one ring, in order.
    ///
    /// # Panics
    ///
    /// Panics if `ring` does not belong to this polygon.
    #[inline]
    pub fn ring(&self, ring: RingId) -> &[Point] {
        let lo = self.ring_starts[ring.0 as usize] as usize;
        let hi = self.ring_starts[ring.0 as usize + 1] as usize;
        &self.verts[lo..hi]
    }

    /// Which ring a vertex belongs to.
    ///
    /// # Panics
    ///
    /// Panics if `v` does not belong to this polygon.
    pub fn ring_of(&self, v: VertexId) -> RingId {
        let i = v.0;
        // Rings are contiguous and ordered, so the ring is the last start <= i.
        let idx = self
            .ring_starts
            .partition_point(|&start| start <= i)
            .saturating_sub(1);
        debug_assert!(idx < self.ring_count());
        RingId(idx as u16)
    }

    /// Iterates every ring's vertices in order, outer ring first.
    pub fn rings(&self) -> impl Iterator<Item = &[Point]> + '_ {
        (0..self.ring_count()).map(move |i| self.ring(RingId(i as u16)))
    }

    /// The vertex following `v` within its ring, wrapping at the ring's end.
    pub fn next_vertex(&self, v: VertexId) -> VertexId {
        let ring = self.ring_of(v);
        let lo = self.ring_starts[ring.0 as usize];
        let hi = self.ring_starts[ring.0 as usize + 1];
        if v.0 + 1 == hi {
            VertexId(lo)
        } else {
            VertexId(v.0 + 1)
        }
    }

    /// The vertex preceding `v` within its ring, wrapping at the ring's start.
    pub fn prev_vertex(&self, v: VertexId) -> VertexId {
        let ring = self.ring_of(v);
        let lo = self.ring_starts[ring.0 as usize];
        let hi = self.ring_starts[ring.0 as usize + 1];
        if v.0 == lo {
            VertexId(hi - 1)
        } else {
            VertexId(v.0 - 1)
        }
    }

    /// The endpoints of an edge, in direction order.
    ///
    /// The polygon's interior lies to the **left** of `start -> end`.
    ///
    /// # Panics
    ///
    /// Panics if `e` does not belong to this polygon.
    #[inline]
    pub fn edge(&self, e: EdgeId) -> (Point, Point) {
        let start = e.start_vertex();
        (self.vertex(start), self.vertex(self.next_vertex(start)))
    }

    /// Total number of edges, which equals the number of vertices.
    #[inline]
    pub fn edge_count(&self) -> usize {
        self.verts.len()
    }

    /// Iterates every edge id.
    pub fn edge_ids(&self) -> impl Iterator<Item = EdgeId> + '_ {
        (0..self.edge_count() as u16).map(EdgeId)
    }

    /// Iterates every vertex id.
    pub fn vertex_ids(&self) -> impl Iterator<Item = VertexId> + '_ {
        (0..self.vertex_count() as u16).map(VertexId)
    }

    /// Whether the interior angle at `v` exceeds 180°, i.e. `v` is a reflex
    /// ("notch") corner.
    ///
    /// Reflex vertices are the only ones that can trigger split events, so this
    /// drives the algorithm's main branch.
    pub fn is_reflex(&self, v: VertexId) -> bool {
        let prev = self.vertex(self.prev_vertex(v));
        let cur = self.vertex(v);
        let next = self.vertex(self.next_vertex(v));
        // Interior is on the left, so a left turn (CCW) is convex.
        orient2d(prev, cur, next) == Orientation::Clockwise
    }

    /// Twice the signed area of the polygon: the outer ring's area minus every
    /// hole's. Always positive for a valid polygon.
    pub fn signed_area2(&self) -> i64 {
        self.rings().map(ring_area2).sum()
    }

    /// Rejects vertices where the ring reverses through exactly 180°.
    ///
    /// Such a corner is a zero-width spike: its wavefront vertex would need
    /// infinite speed, so no finite skeleton exists.
    fn check_no_spikes(&self) -> Result<(), PolygonError> {
        for v in self.vertex_ids() {
            let prev = self.vertex(self.prev_vertex(v));
            let cur = self.vertex(v);
            let next = self.vertex(self.next_vertex(v));

            if orient2d(prev, cur, next) != Orientation::Collinear {
                continue;
            }
            // Collinear at `cur`: either a straight-through vertex (fine, the
            // wavefront just translates) or a 180° reversal (a spike). They are
            // told apart by the sign of the dot product of the two edges.
            let inc = (cur.x as i64 - prev.x as i64, cur.y as i64 - prev.y as i64);
            let out = (next.x as i64 - cur.x as i64, next.y as i64 - cur.y as i64);
            if inc.0 * out.0 + inc.1 * out.1 < 0 {
                return Err(PolygonError::Spike {
                    ring: self.ring_of(v),
                    vertex: v,
                });
            }
        }
        Ok(())
    }

    /// Rejects polygons whose edges cross.
    ///
    /// This is the naive all-pairs test, O(n^2). It is comfortably the
    /// cheapest part of building a skeleton (which is itself O(n^2 log n)), and
    /// keeping it obvious is worth more than the constant factor — consistent
    /// with the crate's correct > understandable > fast ordering.
    fn check_simple(&self) -> Result<(), PolygonError> {
        let n = self.edge_count();
        for i in 0..n {
            let a = EdgeId(i as u16);
            let (a1, a2) = self.edge(a);
            for j in (i + 1)..n {
                let b = EdgeId(j as u16);
                let (b1, b2) = self.edge(b);

                if segments_properly_cross(a1, a2, b1, b2) {
                    return Err(PolygonError::SelfIntersection { a, b });
                }

                // A proper crossing test deliberately ignores touching, since
                // consecutive edges must share a vertex. But two *non*-adjacent
                // edges touching is still an invalid pinch, so check for it.
                if !self.edges_are_adjacent(a, b) && self.edges_touch(a, b) {
                    return Err(PolygonError::SelfIntersection { a, b });
                }
            }
        }
        Ok(())
    }

    /// Whether two edges share a vertex by construction (consecutive in a ring).
    fn edges_are_adjacent(&self, a: EdgeId, b: EdgeId) -> bool {
        let a_start = a.start_vertex();
        let b_start = b.start_vertex();
        self.next_vertex(a_start) == b_start || self.next_vertex(b_start) == a_start
    }

    /// Whether two edges share any point at all.
    fn edges_touch(&self, a: EdgeId, b: EdgeId) -> bool {
        use crate::predicates::point_on_segment;
        let (a1, a2) = self.edge(a);
        let (b1, b2) = self.edge(b);
        point_on_segment(a1, b1, b2)
            || point_on_segment(a2, b1, b2)
            || point_on_segment(b1, a1, a2)
            || point_on_segment(b2, a1, a2)
    }

    /// Rejects holes that escape the outer ring.
    ///
    /// Holes overlapping each other, or poking out of the outer ring, would
    /// already have been caught by [`Polygon::check_simple`] as a crossing.
    /// What remains is a hole entirely *outside* the outer ring, which crosses
    /// nothing — so one containment test per hole closes the gap.
    fn check_holes_inside(&self) -> Result<(), PolygonError> {
        let outer = self.ring(RingId(0));
        for h in 1..self.ring_count() {
            let ring = RingId(h as u16);
            let probe = self.ring(ring)[0];
            if !point_in_ring(probe, outer) {
                return Err(PolygonError::HoleOutsideOuter { ring });
            }
        }
        Ok(())
    }
}

/// Normalises and appends one ring, enforcing the per-ring invariants.
fn push_ring(
    verts: &mut Vec<Point>,
    ring_starts: &mut Vec<u16>,
    ring: &[Point],
    id: RingId,
    want_ccw: bool,
) -> Result<(), PolygonError> {
    // Tolerate the common convention of repeating the first point to close the
    // ring; the crate's own representation leaves rings implicitly closed.
    let ring = match ring {
        [first, mid @ .., last] if first == last && !mid.is_empty() => &ring[..ring.len() - 1],
        _ => ring,
    };

    if ring.len() < 3 {
        return Err(PolygonError::TooFewVertices {
            ring: id,
            count: ring.len(),
        });
    }

    for i in 0..ring.len() {
        let next = (i + 1) % ring.len();
        if ring[i] == ring[next] {
            return Err(PolygonError::RepeatedVertex {
                ring: id,
                index: next,
                point: ring[i],
            });
        }
    }

    let area2 = ring_area2(ring);
    if area2 == 0 {
        return Err(PolygonError::DegenerateRing { ring: id });
    }

    let start = verts.len();
    verts.extend_from_slice(ring);
    if is_ccw(ring) != want_ccw {
        verts[start..].reverse();
    }
    ring_starts.push(verts.len() as u16);
    Ok(())
}

/// Exact point-in-ring test by crossing number.
///
/// Uses only `i64` predicates, so it is exact for every `i16` input. Points
/// exactly on the boundary are reported as inside.
fn point_in_ring(p: Point, ring: &[Point]) -> bool {
    let n = ring.len();
    let mut inside = false;
    for i in 0..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];

        if crate::predicates::point_on_segment(p, a, b) {
            return true;
        }

        // Cast a ray in +x and count crossings. The half-open rule
        // (a.y <= p.y < b.y) counts each crossing exactly once, so vertices
        // touched by the ray don't get double-counted.
        let crosses = (a.y > p.y) != (b.y > p.y);
        if crosses {
            // Is the crossing strictly right of p? Compare exactly, without
            // dividing: the sign of the orientation, flipped when the edge
            // points downward.
            let side = orient2d(a, b, p);
            let upward = b.y > a.y;
            let right_of_p = if upward {
                side == Orientation::Clockwise
            } else {
                side == Orientation::CounterClockwise
            };
            if right_of_p {
                inside = !inside;
            }
        }
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn square(size: i16) -> Vec<Point> {
        vec![
            Point::new(0, 0),
            Point::new(size, 0),
            Point::new(size, size),
            Point::new(0, size),
        ]
    }

    #[test]
    fn builds_a_square() {
        let p = Polygon::from_outer(&square(10)).unwrap();
        assert_eq!(p.vertex_count(), 4);
        assert_eq!(p.edge_count(), 4);
        assert_eq!(p.ring_count(), 1);
        assert_eq!(p.hole_count(), 0);
        assert_eq!(p.signed_area2(), 200);
    }

    #[test]
    fn normalises_outer_ring_to_ccw() {
        let mut cw = square(10);
        cw.reverse();
        let p = Polygon::from_outer(&cw).unwrap();
        assert!(is_ccw(p.ring(RingId(0))), "outer ring must end up CCW");
        assert!(p.signed_area2() > 0);
    }

    #[test]
    fn normalises_holes_to_cw() {
        let hole_ccw = vec![
            Point::new(10, 10),
            Point::new(20, 10),
            Point::new(20, 20),
            Point::new(10, 20),
        ];
        let p = Polygon::new(&square(30), &[hole_ccw]).unwrap();
        assert!(!is_ccw(p.ring(RingId(1))), "hole must end up CW");
        // Outer 30x30 = 900, hole 10x10 = 100. Twice the difference is 1600.
        assert_eq!(p.signed_area2(), 1600);
    }

    #[test]
    fn accepts_explicitly_closed_rings() {
        let mut closed = square(10);
        closed.push(closed[0]);
        let p = Polygon::from_outer(&closed).unwrap();
        assert_eq!(p.vertex_count(), 4, "the repeated closing point is dropped");
    }

    #[test]
    fn edge_and_vertex_ids_correspond() {
        let p = Polygon::from_outer(&square(10)).unwrap();
        for v in p.vertex_ids() {
            assert_eq!(v.outgoing_edge().start_vertex(), v);
            let (start, _) = p.edge(v.outgoing_edge());
            assert_eq!(start, p.vertex(v));
        }
    }

    #[test]
    fn edges_wrap_within_their_ring() {
        let p = Polygon::new(
            &square(30),
            &[vec![
                Point::new(10, 10),
                Point::new(10, 20),
                Point::new(20, 20),
                Point::new(20, 10),
            ]],
        )
        .unwrap();

        // The outer ring's last edge closes back to the outer ring's first
        // vertex, not into the hole.
        assert_eq!(p.next_vertex(VertexId(3)), VertexId(0));
        // The hole's last edge closes back to the hole's first vertex.
        assert_eq!(p.next_vertex(VertexId(7)), VertexId(4));
        assert_eq!(p.prev_vertex(VertexId(4)), VertexId(7));
    }

    #[test]
    fn ring_of_maps_vertices_correctly() {
        let p = Polygon::new(
            &square(30),
            &[vec![
                Point::new(10, 10),
                Point::new(10, 20),
                Point::new(20, 20),
                Point::new(20, 10),
            ]],
        )
        .unwrap();
        for v in 0..4 {
            assert_eq!(p.ring_of(VertexId(v)), RingId(0));
        }
        for v in 4..8 {
            assert_eq!(p.ring_of(VertexId(v)), RingId(1));
        }
    }

    #[test]
    fn rejects_too_few_vertices() {
        let e = Polygon::from_outer(&[Point::new(0, 0), Point::new(1, 1)]).unwrap_err();
        assert!(matches!(e, PolygonError::TooFewVertices { count: 2, .. }));
    }

    #[test]
    fn rejects_empty_outer_ring() {
        assert_eq!(
            Polygon::from_outer(&[]).unwrap_err(),
            PolygonError::NoOuterRing
        );
    }

    #[test]
    fn rejects_repeated_vertices() {
        let e = Polygon::from_outer(&[
            Point::new(0, 0),
            Point::new(10, 0),
            Point::new(10, 0),
            Point::new(10, 10),
        ])
        .unwrap_err();
        assert!(matches!(e, PolygonError::RepeatedVertex { .. }));
    }

    #[test]
    fn rejects_collinear_ring() {
        let e = Polygon::from_outer(&[Point::new(0, 0), Point::new(5, 0), Point::new(9, 0)])
            .unwrap_err();
        assert!(matches!(e, PolygonError::DegenerateRing { .. }));
    }

    #[test]
    fn rejects_self_intersecting_ring() {
        // Asymmetric, so it has non-zero area and must be caught by the
        // crossing test rather than incidentally by the area test.
        let e = Polygon::from_outer(&[
            Point::new(0, 0),
            Point::new(10, 10),
            Point::new(10, 0),
            Point::new(0, 4),
        ])
        .unwrap_err();
        assert!(
            matches!(e, PolygonError::SelfIntersection { .. }),
            "got {e:?}"
        );
    }

    #[test]
    fn rejects_symmetric_bowtie() {
        // A symmetric bowtie's two lobes cancel exactly, so it trips the
        // zero-area check before the crossing check ever runs. Either
        // rejection is correct; what matters is that it does not build.
        let e = Polygon::from_outer(&[
            Point::new(0, 0),
            Point::new(10, 10),
            Point::new(10, 0),
            Point::new(0, 10),
        ])
        .unwrap_err();
        assert_eq!(e, PolygonError::DegenerateRing { ring: RingId(0) });
    }

    #[test]
    fn rejects_zero_width_spike() {
        // Out to (20, 5) and straight back along the same line.
        let e = Polygon::from_outer(&[
            Point::new(0, 0),
            Point::new(10, 0),
            Point::new(10, 5),
            Point::new(20, 5),
            Point::new(10, 5),
            Point::new(0, 10),
        ])
        .unwrap_err();
        assert!(
            matches!(
                e,
                PolygonError::Spike { .. } | PolygonError::SelfIntersection { .. }
            ),
            "got {e:?}"
        );
    }

    #[test]
    fn accepts_straight_through_vertices() {
        // A collinear vertex mid-edge is not a spike; the wavefront handles it.
        let p = Polygon::from_outer(&[
            Point::new(0, 0),
            Point::new(5, 0), // straight-through
            Point::new(10, 0),
            Point::new(10, 10),
            Point::new(0, 10),
        ])
        .unwrap();
        assert_eq!(p.vertex_count(), 5);
    }

    #[test]
    fn rejects_hole_outside_outer() {
        let far_hole = vec![
            Point::new(100, 100),
            Point::new(110, 100),
            Point::new(110, 110),
            Point::new(100, 110),
        ];
        let e = Polygon::new(&square(30), &[far_hole]).unwrap_err();
        assert!(matches!(e, PolygonError::HoleOutsideOuter { .. }));
    }

    #[test]
    fn rejects_overlapping_holes() {
        let a = vec![
            Point::new(5, 5),
            Point::new(15, 5),
            Point::new(15, 15),
            Point::new(5, 15),
        ];
        let b = vec![
            Point::new(10, 10),
            Point::new(20, 10),
            Point::new(20, 20),
            Point::new(10, 20),
        ];
        assert!(Polygon::new(&square(30), &[a, b]).is_err());
    }

    #[test]
    fn rejects_hole_touching_outer_ring() {
        // A hole whose vertex lands on the outer boundary pinches the interior.
        let touching = vec![
            Point::new(0, 10),
            Point::new(10, 10),
            Point::new(10, 20),
            Point::new(0, 20),
        ];
        assert!(Polygon::new(&square(30), &[touching]).is_err());
    }

    #[test]
    fn detects_reflex_vertices() {
        // An L-shape: exactly one reflex corner, at the inner elbow.
        let l = Polygon::from_outer(&[
            Point::new(0, 0),
            Point::new(20, 0),
            Point::new(20, 10),
            Point::new(10, 10), // reflex elbow
            Point::new(10, 20),
            Point::new(0, 20),
        ])
        .unwrap();

        let reflex: Vec<_> = l.vertex_ids().filter(|&v| l.is_reflex(v)).collect();
        assert_eq!(reflex, vec![VertexId(3)]);
    }

    #[test]
    fn convex_polygons_have_no_reflex_vertices() {
        let p = Polygon::from_outer(&square(10)).unwrap();
        assert!(p.vertex_ids().all(|v| !p.is_reflex(v)));
    }

    #[test]
    fn hole_vertices_are_reflex_from_the_interiors_view() {
        // A hole's convex-looking corners bulge *into* the material, so under
        // the interior-on-the-left rule they are reflex. This is what makes
        // holes generate split events.
        let p = Polygon::new(
            &square(30),
            &[vec![
                Point::new(10, 10),
                Point::new(20, 10),
                Point::new(20, 20),
                Point::new(10, 20),
            ]],
        )
        .unwrap();
        for v in 4..8 {
            assert!(p.is_reflex(VertexId(v)), "hole vertex {v} should be reflex");
        }
    }

    #[test]
    fn point_in_ring_basics() {
        let sq = square(10);
        assert!(point_in_ring(Point::new(5, 5), &sq));
        assert!(!point_in_ring(Point::new(15, 5), &sq));
        assert!(!point_in_ring(Point::new(-1, 5), &sq));
        // Boundary counts as inside.
        assert!(point_in_ring(Point::new(0, 5), &sq));
        assert!(point_in_ring(Point::new(0, 0), &sq));
    }

    #[test]
    fn point_in_ring_handles_rays_through_vertices() {
        // A diamond: the +x ray from (0, 0) passes exactly through vertex
        // (10, 0), which a naive crossing count would tally twice.
        let diamond = vec![
            Point::new(10, 0),
            Point::new(20, 10),
            Point::new(10, 20),
            Point::new(0, 10),
        ];
        assert!(!point_in_ring(Point::new(-5, 0), &diamond));
        assert!(point_in_ring(Point::new(10, 10), &diamond));
        assert!(!point_in_ring(Point::new(10, 25), &diamond));
    }

    #[test]
    fn display_for_errors_names_the_ring() {
        let e = PolygonError::HoleOutsideOuter { ring: RingId(2) };
        assert!(alloc::format!("{e}").contains('2'));
    }
}
