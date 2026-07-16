//! The integer point type that forms the crate's public boundary.

use crate::math::Vec2;

/// A point on the integer lattice.
///
/// `i16` is the crate's coordinate type for **both input and output**.
///
/// # The coordinate cap
///
/// Coordinates are restricted to [`Point::MIN_COORD`]`..=`[`Point::MAX_COORD`],
/// i.e. `-16384..=16383` — half of what `i16` could hold. [`Polygon`] rejects
/// anything outside it.
///
/// That one bit is what pays for everything else. It is exactly the bit that
/// lets the crate compute the whole skeleton in `i32` and `f32`, with no `f64`
/// and no `i64` in the algorithm:
///
/// - The orientation determinant of three points needs `2 * d^2` where `d` is
///   the largest coordinate difference. At the full `i16` range that is
///   `8_589_672_450` — it overflows `i32` and reports the *wrong side*. Capped,
///   it is `2_147_352_578`, which fits `i32` with 131_069 to spare, making
///   every predicate exact. See [`crate::predicates`].
/// - `f32` resolves `0.002` at the cap, against `0.004` at the full range,
///   which is what leaves the simulation enough room to work in.
///
/// [`Polygon`]: crate::Polygon
///
/// # Rounding
///
/// A straight skeleton's interior nodes generally land on *irrational*
/// coordinates even when every input vertex is an integer, so there is no
/// lattice to compute *on*. The algorithm works internally in `f32` and rounds
/// only at the boundary, so a [`Node`]'s `position` is the nearest lattice
/// point to its true location. When you need the unrounded value, every node
/// also carries [`Node::exact`].
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

    /// The most negative coordinate a [`Polygon`] may use, `-16384`.
    ///
    /// See [`Point`] for why the range is capped below what `i16` could hold.
    ///
    /// [`Polygon`]: crate::Polygon
    pub const MIN_COORD: i16 = -16384;

    /// The largest coordinate a [`Polygon`](crate::Polygon) may use, `16383`.
    ///
    /// See [`Point`] for why the range is capped below what `i16` could hold.
    pub const MAX_COORD: i16 = 16383;

    /// Whether both coordinates are within
    /// [`MIN_COORD`](Point::MIN_COORD)`..=`[`MAX_COORD`](Point::MAX_COORD).
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::Point;
    ///
    /// assert!(Point::new(16383, -16384).in_range());
    /// assert!(!Point::new(16384, 0).in_range());
    /// assert!(!Point::new(0, i16::MIN).in_range());
    /// ```
    #[inline]
    pub const fn in_range(self) -> bool {
        self.x >= Self::MIN_COORD
            && self.x <= Self::MAX_COORD
            && self.y >= Self::MIN_COORD
            && self.y <= Self::MAX_COORD
    }

    /// Constructs a point from its coordinates.
    #[inline]
    pub const fn new(x: i16, y: i16) -> Self {
        Point { x, y }
    }

    /// Widens to the internal `f32` working space.
    ///
    /// Exact: `i16` needs 16 bits and `f32`'s mantissa holds 24.
    #[inline]
    pub(crate) fn to_vec2(self) -> Vec2 {
        Vec2::new(self.x as f32, self.y as f32)
    }

    /// Rounds an internal `f32` position back to the lattice, saturating at the
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

    /// Squared distance to `other`, computed exactly in `i32`.
    ///
    /// Exact for every in-range pair: the largest possible value is
    /// `2 * 32767^2 = 2_147_352_578`, which fits comfortably.
    ///
    /// Points outside the [coordinate cap](Point::MAX_COORD) saturate rather
    /// than wrap, which keeps comparisons monotone; a [`Polygon`] cannot
    /// contain such a point anyway.
    ///
    /// [`Polygon`]: crate::Polygon
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
        dx.saturating_mul(dx).saturating_add(dy.saturating_mul(dy))
    }
}

/// Rounds half-away-from-zero and saturates into `i16`.
#[inline]
fn round_to_i16(v: f32) -> i16 {
    if v.is_nan() {
        return 0;
    }
    let r = round_half_away_from_zero(v);
    if r <= i16::MIN as f32 {
        i16::MIN
    } else if r >= i16::MAX as f32 {
        i16::MAX
    } else {
        r as i16
    }
}

/// `f32::round` is unavailable in `no_std`, so we spell it out.
#[inline]
pub(crate) fn round_half_away_from_zero(v: f32) -> f32 {
    // `as i32` truncates toward zero; nudging by 0.5 in the sign direction
    // turns that into round-half-away-from-zero for the magnitudes we see.
    if v >= 0.0 {
        let t = (v + 0.5) as i32 as f32;
        // Guard the exact-half-below case introduced by the nudge.
        if t - v > 0.5 {
            t - 1.0
        } else {
            t
        }
    } else {
        let t = (v - 0.5) as i32 as f32;
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
        assert_eq!(round_to_i16(f32::INFINITY), i16::MAX);
        assert_eq!(round_to_i16(f32::NEG_INFINITY), i16::MIN);
        assert_eq!(round_to_i16(f32::NAN), 0);
        assert_eq!(round_to_i16(32767.4), i16::MAX);
        assert_eq!(round_to_i16(-32768.4), i16::MIN);
    }

    #[test]
    fn coordinate_cap_is_where_the_arithmetic_says() {
        // One bit below i16, which is what makes i32 predicates exact.
        assert_eq!(Point::MIN_COORD, -16384);
        assert_eq!(Point::MAX_COORD, 16383);

        assert!(Point::new(0, 0).in_range());
        assert!(Point::new(16383, -16384).in_range());
        assert!(!Point::new(16384, 0).in_range());
        assert!(!Point::new(0, -16385).in_range());
        assert!(!Point::new(i16::MAX, i16::MIN).in_range());
    }

    #[test]
    fn distance_squared_is_exact_across_the_whole_capped_range() {
        let d = Point::new(Point::MIN_COORD, Point::MIN_COORD)
            .distance_squared(Point::new(Point::MAX_COORD, Point::MAX_COORD));
        // 2 * 32767^2, exact rather than saturated.
        assert_eq!(d, 2 * 32767 * 32767);
    }

    #[test]
    fn distance_squared_is_exact_in_range() {
        assert_eq!(Point::new(0, 0).distance_squared(Point::new(3, 4)), 25);
        assert_eq!(Point::new(-3, -4).distance_squared(Point::new(0, 0)), 25);
        assert_eq!(Point::new(5, 5).distance_squared(Point::new(5, 5)), 0);
    }

    #[test]
    fn distance_squared_saturates_beyond_the_cap_rather_than_wrapping() {
        // The full i16 diagonal needs 2 * 65535^2 = 8_589_672_450, which does
        // not fit u32. Such points cannot reach a Polygon, but must saturate
        // rather than wrap if handed here directly.
        let d = Point::new(i16::MIN, i16::MIN).distance_squared(Point::new(i16::MAX, i16::MAX));
        assert_eq!(d, u32::MAX);
        assert!(2u64 * 65535 * 65535 > u32::MAX as u64);
    }

    #[test]
    fn conversions_round_trip() {
        let p = Point::new(-7, 12);
        assert_eq!(Point::from(<(i16, i16)>::from(p)), p);
        assert_eq!(Point::from(<[i16; 2]>::from(p)), p);
    }

    #[test]
    fn to_vec2_is_lossless() {
        // i16 needs 16 bits; f32's mantissa holds 24, so this cannot round.
        for v in [i16::MIN, -1, 0, 1, i16::MAX] {
            let p = Point::new(v, v);
            assert_eq!(Point::from_vec2_rounded(p.to_vec2()), p);
        }
    }
}
