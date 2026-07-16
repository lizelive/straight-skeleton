//! Exact geometric predicates over `i16` inputs.
//!
//! # Why `i64`
//!
//! The crate prefers `i32`/`f32` arithmetic so the algorithm can be ported to a
//! GPU, but the orientation predicate is a deliberate, documented exception —
//! and the width analysis shows why it has to be.
//!
//! For `i16` coordinates in `-32768..=32767`, a coordinate *difference* spans
//! `-65535..=65535`, needing **17 bits**. [`orient2d`] multiplies two such
//! differences, giving a product up to `65535^2 ≈ 2^32`, then subtracts two of
//! them — a result needing **35 bits**. That does not fit in `i32` (31 bits plus
//! sign), and `f32` cannot represent it exactly either, since its mantissa
//! holds only 24 bits.
//!
//! `i64` covers the full 35-bit range with room to spare, so every predicate
//! here is **exact with no rounding, no epsilons, and no overflow** — a
//! property the tests in this module pin down at the coordinate extremes.
//!
//! This affects predicates only. Everything else in the crate stays within the
//! narrower types where it can, and the public API is `i16` throughout.

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

/// Twice the signed area of triangle `(a, b, c)`, computed exactly.
///
/// Positive when `a -> b -> c` turns counter-clockwise, negative when it turns
/// clockwise, and zero exactly when the three points are collinear.
///
/// Because the inputs are `i16` and the accumulator is `i64`, this is exact for
/// **every** input — there is no error term and no epsilon to tune.
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
pub fn signed_area2(a: Point, b: Point, c: Point) -> i64 {
    let abx = b.x as i64 - a.x as i64;
    let aby = b.y as i64 - a.y as i64;
    let acx = c.x as i64 - a.x as i64;
    let acy = c.y as i64 - a.y as i64;
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
/// The accumulator cannot overflow for any ring representable in `i16`: the
/// enclosed area is bounded by the `65535 x 65535` coordinate box, so twice the
/// area is at most about `2^33`.
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
    // vertex (rather than the origin) keeps every term within the 35-bit bound
    // established above, so the i64 sum has no chance of overflowing.
    let origin = ring[0];
    let mut acc: i64 = 0;
    for w in ring.windows(2) {
        acc += signed_area2(origin, w[0], w[1]);
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

    /// The entire reason predicates use `i64`: this case needs 35 bits and
    /// would silently wrap in `i32` or round in `f32`.
    #[test]
    fn orientation_is_exact_at_coordinate_extremes() {
        let a = Point::new(i16::MIN, i16::MIN);
        let b = Point::new(i16::MAX, i16::MIN);
        let c = Point::new(i16::MIN, i16::MAX);
        // Twice the area of this huge right triangle: 65535 * 65535.
        assert_eq!(signed_area2(a, b, c), 65_535i64 * 65_535i64);
        assert_eq!(orient2d(a, b, c), Orientation::CounterClockwise);

        // Confirm the magnitude genuinely exceeds what i32 could hold.
        assert!(65_535i64 * 65_535i64 > i32::MAX as i64);
    }

    #[test]
    fn orientation_detects_one_unit_of_non_collinearity_at_full_scale() {
        let a = Point::new(-32768, -32768);
        let b = Point::new(32767, 32767);
        // Exactly on the line a->b.
        assert_eq!(orient2d(a, b, Point::new(0, 0)), Orientation::Collinear);
        // One unit off the line — must still be detected, not rounded away.
        assert_eq!(
            orient2d(a, b, Point::new(0, 1)),
            Orientation::CounterClockwise
        );
        assert_eq!(orient2d(a, b, Point::new(0, -1)), Orientation::Clockwise);
    }

    /// A real triple, found by search, where the naive `f32` determinant
    /// collapses a genuine turn into a false "collinear".
    ///
    /// `f32` computes each product to within about `2^(32-24) = 256`, so any
    /// true determinant smaller than roughly 512 can be rounded away entirely.
    /// This is why the predicate is not `f32`, and it is not a corner case:
    /// near-collinear triples are exactly what a skeleton's degenerate inputs
    /// look like.
    #[test]
    fn f32_would_report_a_real_turn_as_collinear() {
        let (a, b, c) = (
            Point::new(614, 17634),
            Point::new(13803, 1088),
            Point::new(10216, 5588),
        );

        // The truth: c is (just barely) right of a->b.
        assert_eq!(signed_area2(a, b, c), -2);
        assert_eq!(orient2d(a, b, c), Orientation::Clockwise);

        let naive_f32 = {
            let (ax, ay) = (a.x as f32, a.y as f32);
            let (bx, by) = (b.x as f32, b.y as f32);
            let (cx, cy) = (c.x as f32, c.y as f32);
            (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        };
        assert_eq!(naive_f32, 0.0, "f32 loses the turn entirely");
    }

    /// A real triple, found by search, where the naive `i32` determinant
    /// **overflows and flips sign** — reporting a left turn where the truth is
    /// a right turn. Silent, and catastrophic for any algorithm that branches
    /// on it.
    #[test]
    fn i32_would_overflow_and_report_the_wrong_side() {
        let (a, b, c) = (
            Point::new(21203, -24650),
            Point::new(-22519, 1049),
            Point::new(26449, 26335),
        );

        // The truth, in i64.
        assert_eq!(signed_area2(a, b, c), -2_363_983_124);
        assert_eq!(orient2d(a, b, c), Orientation::Clockwise);

        // The same expression in i32 wraps around.
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
            naive_i32 > 0 && signed_area2(a, b, c) < 0,
            "i32 reports the opposite side"
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
            Point::new(i16::MIN, i16::MIN),
            Point::new(i16::MAX, i16::MIN),
            Point::new(i16::MAX, i16::MAX),
            Point::new(i16::MIN, i16::MAX),
        ];
        assert_eq!(ring_area2(&big), 2 * 65_535i64 * 65_535i64);
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
