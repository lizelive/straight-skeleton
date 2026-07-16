//! The integer point type that forms the crate's public boundary.

use crate::math::Vec2;

/// A point on the integer lattice.
///
/// `i16` is the crate's coordinate type for **both input and output**. The
/// usable range is therefore `-32768..=32767` on each axis, and the longest
/// representable distance is `sqrt(65535^2 + 65535^2) ≈ 92_681`.
///
/// # Rounding
///
/// A straight skeleton's interior nodes generally land on *irrational*
/// coordinates even when every input vertex is an integer — the classic
/// example is the incenter of a 3-4-5 triangle. The algorithm therefore works
/// internally in `f64` and rounds only at the boundary, so a [`Node`]'s
/// `position` is the nearest lattice point to its true location. When you need
/// the unrounded value, every node also carries [`Node::exact`].
///
/// [`Node`]: crate::Node
/// [`Node::exact`]: crate::Node::exact
///
/// # Examples
///
/// ```
/// use straight_skeleton::Point;
///
/// let p = Point::new(3, 4);
/// assert_eq!(p.x, 3);
/// assert_eq!(p.y, 4);
///
/// // Points convert from the obvious tuple and array forms.
/// assert_eq!(Point::from((3, 4)), p);
/// assert_eq!(Point::from([3, 4]), p);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Point {
    /// Horizontal coordinate.
    pub x: i16,
    /// Vertical coordinate.
    pub y: i16,
}

impl Point {
    /// The origin, `(0, 0)`.
    pub const ORIGIN: Point = Point { x: 0, y: 0 };

    /// Constructs a point from its coordinates.
    #[inline]
    pub const fn new(x: i16, y: i16) -> Self {
        Point { x, y }
    }

    /// Widens to the internal `f64` working space.
    #[inline]
    pub(crate) fn to_vec2(self) -> Vec2 {
        Vec2::new(self.x as f64, self.y as f64)
    }

    /// Rounds an internal `f64` position back to the lattice, saturating at the
    /// `i16` bounds rather than wrapping.
    ///
    /// Saturation is a safety net, not an expected path: a straight skeleton
    /// lies inside the convex hull of its input, so a node can only exceed the
    /// input's range through floating-point error at the very edge of the
    /// coordinate space.
    #[inline]
    pub(crate) fn from_vec2_rounded(v: Vec2) -> Self {
        Point {
            x: round_to_i16(v.x),
            y: round_to_i16(v.y),
        }
    }

    /// Squared distance to `other`.
    ///
    /// Exact for every pair of points except those spanning nearly the whole
    /// coordinate space: the maximum possible value is `2 * 65535^2`, which
    /// overflows `u32` by a hair, so the result saturates at [`u32::MAX`]
    /// rather than wrapping. Comparisons stay monotone in all other cases.
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::Point;
    ///
    /// assert_eq!(Point::new(0, 0).distance_squared(Point::new(3, 4)), 25);
    /// ```
    #[inline]
    pub fn distance_squared(self, other: Point) -> u32 {
        let dx = (self.x as i32 - other.x as i32).unsigned_abs();
        let dy = (self.y as i32 - other.y as i32).unsigned_abs();
        // dx, dy <= 65_535, so dx^2 + dy^2 <= 8_589_803_970 > u32::MAX.
        // Saturate rather than wrap; callers comparing distances get a
        // monotone answer everywhere except the extreme corner case.
        (dx * dx).saturating_add(dy * dy)
    }
}

/// Rounds half-away-from-zero and saturates into `i16`.
#[inline]
fn round_to_i16(v: f64) -> i16 {
    if v.is_nan() {
        return 0;
    }
    let r = round_half_away_from_zero(v);
    if r <= i16::MIN as f64 {
        i16::MIN
    } else if r >= i16::MAX as f64 {
        i16::MAX
    } else {
        r as i16
    }
}

/// `f64::round` is unavailable in `no_std`, so we spell it out.
#[inline]
fn round_half_away_from_zero(v: f64) -> f64 {
    // `as i64` truncates toward zero; nudging by 0.5 in the sign direction
    // turns that into round-half-away-from-zero for the magnitudes we see.
    if v >= 0.0 {
        let t = (v + 0.5) as i64 as f64;
        // Guard the exact-half-below case introduced by the nudge.
        if t - v > 0.5 {
            t - 1.0
        } else {
            t
        }
    } else {
        let t = (v - 0.5) as i64 as f64;
        if v - t > 0.5 {
            t + 1.0
        } else {
            t
        }
    }
}

impl From<(i16, i16)> for Point {
    #[inline]
    fn from((x, y): (i16, i16)) -> Self {
        Point::new(x, y)
    }
}

impl From<[i16; 2]> for Point {
    #[inline]
    fn from([x, y]: [i16; 2]) -> Self {
        Point::new(x, y)
    }
}

impl From<Point> for (i16, i16) {
    #[inline]
    fn from(p: Point) -> Self {
        (p.x, p.y)
    }
}

impl From<Point> for [i16; 2] {
    #[inline]
    fn from(p: Point) -> Self {
        [p.x, p.y]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rounds_to_nearest() {
        assert_eq!(round_to_i16(0.4), 0);
        assert_eq!(round_to_i16(0.5), 1);
        assert_eq!(round_to_i16(0.6), 1);
        assert_eq!(round_to_i16(-0.4), 0);
        assert_eq!(round_to_i16(-0.5), -1);
        assert_eq!(round_to_i16(-0.6), -1);
        assert_eq!(round_to_i16(1.5), 2);
        assert_eq!(round_to_i16(2.5), 3);
        assert_eq!(round_to_i16(-2.5), -3);
    }

    #[test]
    fn rounding_saturates_instead_of_wrapping() {
        assert_eq!(round_to_i16(1e9), i16::MAX);
        assert_eq!(round_to_i16(-1e9), i16::MIN);
        assert_eq!(round_to_i16(f64::INFINITY), i16::MAX);
        assert_eq!(round_to_i16(f64::NEG_INFINITY), i16::MIN);
        assert_eq!(round_to_i16(f64::NAN), 0);
        assert_eq!(round_to_i16(32767.4), i16::MAX);
        assert_eq!(round_to_i16(-32768.4), i16::MIN);
    }

    #[test]
    fn distance_squared_is_exact_in_range() {
        assert_eq!(Point::new(0, 0).distance_squared(Point::new(3, 4)), 25);
        assert_eq!(Point::new(-3, -4).distance_squared(Point::new(0, 0)), 25);
        assert_eq!(Point::new(5, 5).distance_squared(Point::new(5, 5)), 0);
    }

    #[test]
    fn distance_squared_saturates_at_the_extreme_corner() {
        // The full diagonal overflows u32; saturating keeps it monotone.
        let d = Point::new(i16::MIN, i16::MIN).distance_squared(Point::new(i16::MAX, i16::MAX));
        assert_eq!(d, u32::MAX);
    }

    #[test]
    fn conversions_round_trip() {
        let p = Point::new(-7, 12);
        assert_eq!(Point::from(<(i16, i16)>::from(p)), p);
        assert_eq!(Point::from(<[i16; 2]>::from(p)), p);
    }

    #[test]
    fn to_vec2_is_lossless() {
        for v in [i16::MIN, -1, 0, 1, i16::MAX] {
            let p = Point::new(v, v);
            assert_eq!(Point::from_vec2_rounded(p.to_vec2()), p);
        }
    }
}
