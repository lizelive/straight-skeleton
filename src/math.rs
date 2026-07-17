//! Internal floating-point scaffolding.
//!
//! # No `f64`, anywhere
//!
//! The whole crate computes in `f32`. Nothing here, or anywhere else, uses
//! `f64`: the algorithm is meant to be portable to hardware where `f64` is
//! slow or absent, and a type you only use "internally" is still a type the
//! hardware has to have.
//!
//! That is affordable only because the coordinate range is capped at
//! [`Point::MAX_COORD`], which is what buys `f32` enough absolute resolution to
//! work in. `docs/DESIGN.md` has the arithmetic.
//!
//! # No dependencies either
//!
//! `no_std` with **zero required dependencies** rules out both `std`'s
//! `f32::sqrt` and the `libm` crate. A square root is the only transcendental
//! the algorithm needs — solely to normalise edge normals, see
//! [`crate::wavefront`] — so we carry a small implementation rather than take a
//! dependency for one function.
//!
//! With the `std` feature on we defer to the hardware instruction, which is
//! both faster and correctly rounded.
//!
//! [`Point::MAX_COORD`]: crate::Point::MAX_COORD

use core::ops::{Add, AddAssign, Mul, Neg, Sub};

/// Square root of a non-negative `f32`.
///
/// With the `std` feature this is `f32::sqrt` (correctly rounded, typically a
/// single CPU instruction). Without it, this is a Newton–Raphson refinement
/// seeded by the classic exponent-halving bit trick.
///
/// The `no_std` path agrees with the `std` path to within 1 ULP for all finite
/// non-negative inputs, which is verified by a test in this module.
#[inline]
pub fn sqrt(x: f32) -> f32 {
    #[cfg(feature = "std")]
    {
        x.sqrt()
    }
    #[cfg(not(feature = "std"))]
    {
        sqrt_soft(x)
    }
}

/// Dependency-free `sqrt`, always compiled so it can be differentially tested
/// against the `std` implementation.
#[inline]
#[allow(dead_code)]
pub(crate) fn sqrt_soft(x: f32) -> f32 {
    if x.is_nan() || x < 0.0 {
        return f32::NAN;
    }
    if x == 0.0 || x == f32::INFINITY {
        // Preserves the sign of zero, matching `f32::sqrt`.
        return x;
    }

    // Denormals break the exponent trick below, so scale them up into the
    // normal range, take the root there, and scale the answer back by the
    // square root of the same factor.
    if x < f32::MIN_POSITIVE {
        return sqrt_soft(x * 16_777_216.0) * (1.0 / 4096.0);
    }

    // Seed: halving the biased exponent approximates halving the exponent,
    // which approximates a square root to within a factor of ~2.
    let bits = x.to_bits();
    let mut y = f32::from_bits((bits >> 1) + 0x1fc0_0000);

    // Newton–Raphson on f(y) = y^2 - x. Each step doubles the correct digits;
    // the seed is good to ~5 bits, so 4 steps saturate f32's 24-bit mantissa.
    for _ in 0..4 {
        y = 0.5 * (y + x / y);
    }
    y
}

/// `x.floor() as i32`, without `std`.
///
/// `f32::floor` lives in `std`, and the crate takes no dependency on `libm`, so
/// the one place that needs a floor carries it here. `as i32` truncates toward
/// zero, which is the floor only for non-negative input; the correction is one
/// step down for negatives that were not already integral.
///
/// Saturates for input outside `i32`, which is what `as i32` already does.
#[inline]
pub fn floor_i32(x: f32) -> i32 {
    let t = x as i32;
    if x < t as f32 {
        t - 1
    } else {
        t
    }
}

/// A 2D vector in the algorithm's internal `f32` working space.
///
/// Public input and output coordinates are `i16` (see [`crate::Point`]); this
/// type exists only between them. See `docs/DESIGN.md` for why `f32` is enough,
/// and what the coordinate cap has to do with it.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec2 {
    /// Horizontal component.
    pub x: f32,
    /// Vertical component.
    pub y: f32,
}

impl Vec2 {
    /// The zero vector.
    pub const ZERO: Vec2 = Vec2 { x: 0.0, y: 0.0 };

    /// Constructs a vector from its components.
    #[inline]
    pub const fn new(x: f32, y: f32) -> Self {
        Vec2 { x, y }
    }

    /// Dot product `self · other`.
    #[inline]
    pub fn dot(self, other: Vec2) -> f32 {
        self.x * other.x + self.y * other.y
    }

    /// 2D cross product `self × other`, i.e. the signed area of the
    /// parallelogram they span. Positive when `other` is counter-clockwise
    /// from `self`.
    #[inline]
    pub fn cross(self, other: Vec2) -> f32 {
        self.x * other.y - self.y * other.x
    }

    /// Euclidean length.
    #[inline]
    pub fn length(self) -> f32 {
        sqrt(self.dot(self))
    }

    /// Squared Euclidean length. Prefer this over [`Vec2::length`] when
    /// comparing magnitudes, since it avoids a square root.
    #[inline]
    pub fn length_squared(self) -> f32 {
        self.dot(self)
    }

    /// Returns `self` scaled to unit length, or `None` if `self` is too short
    /// for the direction to be meaningful.
    #[inline]
    pub fn normalize(self) -> Option<Vec2> {
        let len = self.length();
        if len > 0.0 && len.is_finite() {
            Some(Vec2::new(self.x / len, self.y / len))
        } else {
            None
        }
    }

    /// Rotates 90° counter-clockwise: `(x, y) -> (-y, x)`.
    ///
    /// For an edge pointing along `d`, this yields the normal facing the
    /// polygon interior under the crate's CCW-outer-ring convention.
    #[inline]
    pub fn perp(self) -> Vec2 {
        Vec2::new(-self.y, self.x)
    }

    /// True when both components are finite.
    #[inline]
    #[allow(dead_code)]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }
}

impl Add for Vec2 {
    type Output = Vec2;
    #[inline]
    fn add(self, rhs: Vec2) -> Vec2 {
        Vec2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl AddAssign for Vec2 {
    #[inline]
    fn add_assign(&mut self, rhs: Vec2) {
        *self = *self + rhs;
    }
}

impl Sub for Vec2 {
    type Output = Vec2;
    #[inline]
    fn sub(self, rhs: Vec2) -> Vec2 {
        Vec2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl Mul<f32> for Vec2 {
    type Output = Vec2;
    #[inline]
    fn mul(self, rhs: f32) -> Vec2 {
        Vec2::new(self.x * rhs, self.y * rhs)
    }
}

impl Neg for Vec2 {
    type Output = Vec2;
    #[inline]
    fn neg(self) -> Vec2 {
        Vec2::new(-self.x, -self.y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The `no_std` square root must track the hardware one closely enough that
    /// swapping features cannot change which branch the algorithm takes.
    #[test]
    fn soft_sqrt_matches_std_within_one_ulp() {
        let cases = [
            1.0,
            2.0,
            4.0,
            1e-30,
            1e30,
            0.5,
            32767.0,
            32767.0 * 32767.0 * 2.0,
            f32::MIN_POSITIVE,
            f32::MAX,
            3.0000001,
        ];
        for x in cases {
            let soft = sqrt_soft(x);
            let hard = x.sqrt();
            let ulps = (soft.to_bits() as i64 - hard.to_bits() as i64).abs();
            assert!(ulps <= 1, "sqrt({x}): soft={soft} hard={hard} ulps={ulps}");
        }
    }

    #[test]
    fn soft_sqrt_matches_std_over_a_sweep() {
        // Walk a wide range of magnitudes and mantissas.
        let mut x = 1e-12f32;
        while x < 1e12 {
            let ulps = (sqrt_soft(x).to_bits() as i64 - x.sqrt().to_bits() as i64).abs();
            assert!(ulps <= 1, "sqrt({x}) differed by {ulps} ulps");
            x *= 1.000_137;
        }
    }

    /// Denormals defeat the exponent-halving seed, so they take a scaled path.
    #[test]
    fn soft_sqrt_handles_denormals() {
        for x in [1e-40f32, 1e-44, f32::MIN_POSITIVE / 2.0, f32::from_bits(1)] {
            let ulps = (sqrt_soft(x).to_bits() as i64 - x.sqrt().to_bits() as i64).abs();
            assert!(ulps <= 1, "sqrt({x:e}) differed by {ulps} ulps");
        }
    }

    #[test]
    fn soft_sqrt_edge_values() {
        assert_eq!(sqrt_soft(0.0), 0.0);
        assert!(sqrt_soft(0.0).is_sign_positive());
        assert!(sqrt_soft(-0.0).is_sign_negative(), "must preserve -0.0");
        assert!(sqrt_soft(-1.0).is_nan());
        assert!(sqrt_soft(f32::NAN).is_nan());
        assert_eq!(sqrt_soft(f32::INFINITY), f32::INFINITY);
    }

    #[test]
    fn floor_i32_matches_std_floor() {
        for x in [
            0.0f32,
            1.0,
            1.5,
            -1.0,
            -1.5,
            -0.5,
            0.5,
            1e6,
            -1e6,
            16384.0,
            -16384.0,
            1_638_400.0,
            -1_638_400.0,
            -0.0,
        ] {
            assert_eq!(floor_i32(x), x.floor() as i32, "floor_i32({x})");
        }
    }

    #[test]
    fn floor_i32_over_a_sweep() {
        let mut x = -3000.0f32;
        while x < 3000.0 {
            assert_eq!(floor_i32(x), x.floor() as i32, "floor_i32({x})");
            x += 0.37;
        }
    }

    #[test]
    fn vec2_perp_faces_left() {
        // +x edge direction -> +y normal, i.e. interior-left for CCW rings.
        assert_eq!(Vec2::new(1.0, 0.0).perp(), Vec2::new(0.0, 1.0));
    }

    #[test]
    fn vec2_cross_sign_is_ccw_positive() {
        assert!(Vec2::new(1.0, 0.0).cross(Vec2::new(0.0, 1.0)) > 0.0);
        assert!(Vec2::new(0.0, 1.0).cross(Vec2::new(1.0, 0.0)) < 0.0);
    }

    #[test]
    fn vec2_normalize_rejects_degenerate() {
        assert!(Vec2::ZERO.normalize().is_none());
        let n = Vec2::new(3.0, 4.0).normalize().unwrap();
        assert!((n.length() - 1.0).abs() < 1e-6);
    }
}
