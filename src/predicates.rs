//! Exact geometric predicates, in `i32`.
//!
//! # Why the coordinate cap exists
//!
//! [`orient2d`] multiplies two coordinate *differences* and subtracts:
//!
//! ```text
//!     (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
//! ```
//!
//! so the width it needs is set by the largest difference `d`, as `2 * d^2`.
//! That single expression is what fixes the crate's coordinate range:
//!
//! | coordinates | largest `d` | `2 * d^2` | in `i32`? |
//! |---|---|---|---|
//! | full `i16`, `-32768..=32767` | 65_535 | 8_589_672_450 | **overflows**, and wraps to the *wrong side* |
//! | capped, `-16384..=16383` | 32_767 | 2_147_352_578 | fits, with 131_069 to spare |
//!
//! Giving up one bit of range therefore buys **exact** predicates in `i32`: no
//! epsilon, no rounding, no overflow, for every input a [`Polygon`] will
//! accept. The tests in this module pin down both halves of that table —
//! including real triples where `i32` at the full range, and `f32` at any
//! range, each get the answer wrong.
//!
//! `f32` is not an option here at any range: its mantissa holds 24 bits, and
//! these products need up to 31.
//!
//! [`Polygon`]: crate::Polygon

use crate::Point;

/// Which side of a directed line a point falls on.
///
/// Returned by [`orient2d`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orientation {
    /// `c` lies strictly left of the directed line `a -> b`
    /// (a counter-clockwise turn).
    CounterClockwise,
    /// `c` lies exactly on the line through `a` and `b`.
    Collinear,
    /// `c` lies strictly right of the directed line `a -> b`
    /// (a clockwise turn).
    Clockwise,
}

/// Twice the signed area of triangle `(a, b, c)`, computed exactly in `i32`.
///
/// Positive when `a -> b -> c` turns counter-clockwise, negative when it turns
/// clockwise, and zero exactly when the three points are collinear.
///
/// Exact — no error term, no epsilon — for every point within the
/// [coordinate cap](crate::Point::MAX_COORD), which is precisely the range
/// where it cannot overflow. See the [module docs](self) for the arithmetic.
///
/// # Panics
///
/// In debug builds, if a coordinate is outside the cap, where the result would
/// silently wrap. [`Polygon`] rejects such points, so this cannot be reached
/// through the normal API.
///
/// [`Polygon`]: crate::Polygon
///
/// # Examples
///
/// ```
/// use straight_skeleton::predicates::signed_area2;
/// use straight_skeleton::Point;
///
/// let a = Point::new(0, 0);
/// let b = Point::new(4, 0);
/// // The unit triangle has area 1/2, so twice its area is 1.
/// assert_eq!(signed_area2(Point::new(0, 0), Point::new(1, 0), Point::new(0, 1)), 1);
/// // Counter-clockwise is positive, clockwise negative, collinear zero.
/// assert!(signed_area2(a, b, Point::new(0, 3)) > 0);
/// assert!(signed_area2(a, b, Point::new(0, -3)) < 0);
/// assert_eq!(signed_area2(a, b, Point::new(9, 0)), 0);
/// ```
#[inline]
pub fn signed_area2(a: Point, b: Point, c: Point) -> i32 {
    debug_assert!(
        a.in_range() && b.in_range() && c.in_range(),
        "signed_area2 is only exact within the coordinate cap"
    );
    let abx = b.x as i32 - a.x as i32;
    let aby = b.y as i32 - a.y as i32;
    let acx = c.x as i32 - a.x as i32;
    let acy = c.y as i32 - a.y as i32;
    abx * acy - aby * acx
}

/// Exact orientation of the point triple `(a, b, c)`.
///
/// This is the sign of [`signed_area2`], named for readability at call sites.
///
/// # Examples
///
/// ```
/// use straight_skeleton::predicates::{orient2d, Orientation};
/// use straight_skeleton::Point;
///
/// let a = Point::new(0, 0);
/// let b = Point::new(4, 0);
/// assert_eq!(orient2d(a, b, Point::new(2, 1)), Orientation::CounterClockwise);
/// assert_eq!(orient2d(a, b, Point::new(2, -1)), Orientation::Clockwise);
/// assert_eq!(orient2d(a, b, Point::new(2, 0)), Orientation::Collinear);
/// ```
#[inline]
pub fn orient2d(a: Point, b: Point, c: Point) -> Orientation {
    match signed_area2(a, b, c) {
        d if d > 0 => Orientation::CounterClockwise,
        d if d < 0 => Orientation::Clockwise,
        _ => Orientation::Collinear,
    }
}

/// Twice the signed area enclosed by a closed ring, computed exactly.
///
/// Positive for a counter-clockwise ring, negative for a clockwise one, and
/// zero for a ring enclosing no area (all vertices collinear, or a ring that
/// doubles back on itself exactly).
///
/// # Why this one is `i64`
///
/// Unlike [`orient2d`], this *sums* triangles, and a ring that doubles back can
/// wind around a region more than once, so the running total is bounded by the
/// vertex count rather than by the coordinate box. One bit of range cannot fix
/// that; only a wider accumulator can.
///
/// It is a fair exception to make, because this is validation — it runs once,
/// on the host, when a [`Polygon`] is built, and never during the simulation.
/// The skeleton itself uses only `i32` and `f32`.
///
/// [`Polygon`]: crate::Polygon
///
/// # Examples
///
/// ```
/// use straight_skeleton::predicates::ring_area2;
/// use straight_skeleton::Point;
///
/// let ccw = [Point::new(0, 0), Point::new(10, 0), Point::new(10, 10), Point::new(0, 10)];
/// assert_eq!(ring_area2(&ccw), 200); // twice the 10x10 area
///
/// let mut cw = ccw;
/// cw.reverse();
/// assert_eq!(ring_area2(&cw), -200);
/// ```
pub fn ring_area2(ring: &[Point]) -> i64 {
    if ring.len() < 3 {
        return 0;
    }
    // The shoelace formula about the ring's first vertex. Anchoring at a real
    // vertex, rather than at the origin, keeps every individual term inside the
    // exact i32 range that `signed_area2` needs.
    let origin = ring[0];
    let mut acc: i64 = 0;
    for w in ring.windows(2) {
        acc += signed_area2(origin, w[0], w[1]) as i64;
    }
    acc
}

/// Whether a ring winds counter-clockwise (encloses positive area).
///
/// # Examples
///
/// ```
/// use straight_skeleton::predicates::is_ccw;
/// use straight_skeleton::Point;
///
/// let ring = [Point::new(0, 0), Point::new(4, 0), Point::new(4, 4), Point::new(0, 4)];
/// assert!(is_ccw(&ring));
/// ```
#[inline]
pub fn is_ccw(ring: &[Point]) -> bool {
    ring_area2(ring) > 0
}

/// Whether closed segments `a1-a2` and `b1-b2` properly cross.
///
/// "Properly" means they intersect at a single point interior to both segments.
/// Segments that merely touch at an endpoint, or that overlap collinearly, are
/// **not** proper crossings and return `false`. This is exactly the test needed
/// for detecting self-intersecting rings, where shared endpoints between
/// consecutive edges are expected and legal.
#[inline]
pub fn segments_properly_cross(a1: Point, a2: Point, b1: Point, b2: Point) -> bool {
    let d1 = orient2d(a1, a2, b1);
    let d2 = orient2d(a1, a2, b2);
    let d3 = orient2d(b1, b2, a1);
    let d4 = orient2d(b1, b2, a2);

    // A strict crossing requires each segment to straddle the other's line.
    // Any Collinear result means a touch or an overlap, not a proper crossing.
    d1 != d2
        && d3 != d4
        && d1 != Orientation::Collinear
        && d2 != Orientation::Collinear
        && d3 != Orientation::Collinear
        && d4 != Orientation::Collinear
}

/// Whether point `p` lies on the closed segment `a-b`.
#[inline]
pub fn point_on_segment(p: Point, a: Point, b: Point) -> bool {
    if orient2d(a, b, p) != Orientation::Collinear {
        return false;
    }
    // Collinear: p is on the segment iff it is inside the bounding box.
    p.x >= a.x.min(b.x) && p.x <= a.x.max(b.x) && p.y >= a.y.min(b.y) && p.y <= a.y.max(b.y)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orientation_basics() {
        let a = Point::new(0, 0);
        let b = Point::new(1, 0);
        assert_eq!(
            orient2d(a, b, Point::new(0, 1)),
            Orientation::CounterClockwise
        );
        assert_eq!(orient2d(a, b, Point::new(0, -1)), Orientation::Clockwise);
        assert_eq!(orient2d(a, b, Point::new(2, 0)), Orientation::Collinear);
    }

    #[test]
    fn orientation_is_antisymmetric_under_swap() {
        let a = Point::new(-5, 3);
        let b = Point::new(7, -2);
        let c = Point::new(1, 9);
        assert_eq!(signed_area2(a, b, c), -signed_area2(a, c, b));
        assert_eq!(signed_area2(a, b, c), signed_area2(b, c, a));
        assert_eq!(signed_area2(a, b, c), signed_area2(c, a, b));
    }

    /// The cap is exactly as wide as `i32` allows, and no wider: the largest
    /// triangle it admits lands just inside the type.
    #[test]
    fn orientation_is_exact_at_the_coordinate_cap() {
        let a = Point::new(Point::MIN_COORD, Point::MIN_COORD);
        let b = Point::new(Point::MAX_COORD, Point::MIN_COORD);
        let c = Point::new(Point::MIN_COORD, Point::MAX_COORD);

        // Twice the area of the largest right triangle that fits: 32767^2.
        assert_eq!(signed_area2(a, b, c), 32_767 * 32_767);
        assert_eq!(orient2d(a, b, c), Orientation::CounterClockwise);

        // The worst case is twice that, and it still fits -- which is what
        // makes the cap exactly one bit, rather than two.
        assert!(2i64 * 32_767 * 32_767 <= i32::MAX as i64);
        assert!(
            2i64 * 65_535 * 65_535 > i32::MAX as i64,
            "the uncapped range would not fit, which is why the cap exists"
        );
    }

    #[test]
    fn orientation_detects_one_unit_of_non_collinearity_at_full_scale() {
        let a = Point::new(Point::MIN_COORD, Point::MIN_COORD);
        let b = Point::new(Point::MAX_COORD, Point::MAX_COORD);
        // Exactly on the line a->b.
        assert_eq!(orient2d(a, b, Point::new(0, 0)), Orientation::Collinear);
        // One unit off the line — must still be detected, not rounded away.
        assert_eq!(
            orient2d(a, b, Point::new(0, 1)),
            Orientation::CounterClockwise
        );
        assert_eq!(orient2d(a, b, Point::new(0, -1)), Orientation::Clockwise);
    }

    /// A real triple, found by search, **inside the coordinate cap**, where the
    /// naive `f32` determinant collapses a genuine turn into a false
    /// "collinear".
    ///
    /// This is why capping the range does not rescue `f32`, and so why the
    /// predicate is `i32` instead. `f32` rounds each product by up to `2^6`
    /// here, so any true determinant below that can vanish entirely — and
    /// near-collinear triples are exactly what degenerate input looks like.
    #[test]
    fn f32_would_report_a_real_turn_as_collinear_even_inside_the_cap() {
        let (a, b, c) = (
            Point::new(14176, -12146),
            Point::new(-9937, 5341),
            Point::new(4434, -5081),
        );
        assert!(
            a.in_range() && b.in_range() && c.in_range(),
            "the point is that capping the range does not save f32"
        );

        // The truth: c is (just barely) left of a->b.
        assert_eq!(signed_area2(a, b, c), 9);
        assert_eq!(orient2d(a, b, c), Orientation::CounterClockwise);

        let naive_f32 = {
            let (ax, ay) = (a.x as f32, a.y as f32);
            let (bx, by) = (b.x as f32, b.y as f32);
            let (cx, cy) = (c.x as f32, c.y as f32);
            (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        };
        assert_eq!(naive_f32, 0.0, "f32 loses the turn entirely");
    }

    /// A real triple, found by search, showing why the cap is not optional:
    /// **beyond** it the `i32` determinant overflows and flips sign, reporting
    /// a left turn where the truth is a right turn. Silent, and catastrophic
    /// for any algorithm that branches on it.
    ///
    /// These points are outside [`Point::MAX_COORD`], so a [`Polygon`] rejects
    /// them and [`signed_area2`] is never asked about them. This test is what
    /// justifies that rejection.
    ///
    /// [`Polygon`]: crate::Polygon
    #[test]
    fn beyond_the_cap_i32_would_report_the_wrong_side() {
        let (a, b, c) = (
            Point::new(21203, -24650),
            Point::new(-22519, 1049),
            Point::new(26449, 26335),
        );
        assert!(
            !a.in_range() && !b.in_range() && !c.in_range(),
            "the whole point is that these are out of range"
        );

        // The truth, computed in i64.
        let exact = |p: Point, q: Point, r: Point| -> i64 {
            (q.x as i64 - p.x as i64) * (r.y as i64 - p.y as i64)
                - (q.y as i64 - p.y as i64) * (r.x as i64 - p.x as i64)
        };
        assert_eq!(exact(a, b, c), -2_363_983_124);

        // The same expression in i32 wraps clean around to positive.
        let naive_i32 = {
            let (ax, ay) = (a.x as i32, a.y as i32);
            let (bx, by) = (b.x as i32, b.y as i32);
            let (cx, cy) = (c.x as i32, c.y as i32);
            (bx - ax)
                .wrapping_mul(cy - ay)
                .wrapping_sub((by - ay).wrapping_mul(cx - ax))
        };
        assert_eq!(naive_i32, 1_930_984_172);
        assert!(
            naive_i32 > 0 && exact(a, b, c) < 0,
            "i32 reports the opposite side once the cap is exceeded"
        );
    }

    #[test]
    fn ring_area_signs_and_magnitude() {
        let square = [
            Point::new(0, 0),
            Point::new(10, 0),
            Point::new(10, 10),
            Point::new(0, 10),
        ];
        assert_eq!(ring_area2(&square), 200);
        assert!(is_ccw(&square));

        let mut rev = square;
        rev.reverse();
        assert_eq!(ring_area2(&rev), -200);
        assert!(!is_ccw(&rev));
    }

    #[test]
    fn ring_area_is_translation_invariant() {
        let tri = [Point::new(0, 0), Point::new(30, 0), Point::new(0, 40)];
        let shifted: [Point; 3] =
            core::array::from_fn(|i| Point::new(tri[i].x - 1000, tri[i].y + 500));
        assert_eq!(ring_area2(&tri), ring_area2(&shifted));
        assert_eq!(ring_area2(&tri), 1200);
    }

    #[test]
    fn ring_area_of_degenerate_rings_is_zero() {
        assert_eq!(ring_area2(&[]), 0);
        assert_eq!(ring_area2(&[Point::new(1, 1)]), 0);
        assert_eq!(ring_area2(&[Point::new(1, 1), Point::new(2, 2)]), 0);
        // Collinear "ring" encloses nothing.
        assert_eq!(
            ring_area2(&[Point::new(0, 0), Point::new(5, 0), Point::new(9, 0)]),
            0
        );
    }

    #[test]
    fn ring_area_does_not_overflow_at_full_extent() {
        let big = [
            Point::new(Point::MIN_COORD, Point::MIN_COORD),
            Point::new(Point::MAX_COORD, Point::MIN_COORD),
            Point::new(Point::MAX_COORD, Point::MAX_COORD),
            Point::new(Point::MIN_COORD, Point::MAX_COORD),
        ];
        assert_eq!(ring_area2(&big), 2 * 32_767i64 * 32_767i64);
    }

    #[test]
    fn proper_crossing_detection() {
        // A clean X.
        assert!(segments_properly_cross(
            Point::new(0, 0),
            Point::new(10, 10),
            Point::new(0, 10),
            Point::new(10, 0),
        ));
        // Disjoint.
        assert!(!segments_properly_cross(
            Point::new(0, 0),
            Point::new(1, 1),
            Point::new(5, 5),
            Point::new(6, 6),
        ));
    }

    #[test]
    fn touching_endpoints_are_not_proper_crossings() {
        // Consecutive polygon edges share a vertex; that must stay legal.
        assert!(!segments_properly_cross(
            Point::new(0, 0),
            Point::new(5, 0),
            Point::new(5, 0),
            Point::new(5, 5),
        ));
        // T-junction: endpoint lands in the other segment's interior.
        assert!(!segments_properly_cross(
            Point::new(0, 0),
            Point::new(10, 0),
            Point::new(5, 0),
            Point::new(5, 5),
        ));
    }

    #[test]
    fn collinear_overlap_is_not_a_proper_crossing() {
        assert!(!segments_properly_cross(
            Point::new(0, 0),
            Point::new(10, 0),
            Point::new(5, 0),
            Point::new(15, 0),
        ));
    }

    #[test]
    fn point_on_segment_cases() {
        let a = Point::new(0, 0);
        let b = Point::new(10, 5);
        assert!(point_on_segment(Point::new(4, 2), a, b));
        assert!(point_on_segment(a, a, b), "endpoints count as on-segment");
        assert!(point_on_segment(b, a, b));
        // Collinear with the line but beyond the segment.
        assert!(!point_on_segment(Point::new(20, 10), a, b));
        assert!(!point_on_segment(Point::new(-2, -1), a, b));
        // Off the line entirely.
        assert!(!point_on_segment(Point::new(4, 3), a, b));
    }
}
