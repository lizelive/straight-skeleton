//! Optional conversions to and from other ecosystem crates.
//!
//! Every integration here is opt-in and off by default, so the crate keeps its
//! zero required dependencies. Each is gated behind a feature named after the
//! crate it bridges to.
//!
//! # Conversions are lossy in one direction only
//!
//! Going *out* of this crate — `Point -> anything` — is always exact: an `i16`
//! fits in every target type without rounding.
//!
//! Coming *in* is not, because a point on some other crate's `f32`/`f64` plane
//! need not land on the `i16` lattice at all. So inbound conversions from
//! floating-point types are [`TryFrom`], and they fail rather than silently
//! round or wrap. Integer inbound conversions are infallible where the width
//! allows, and `TryFrom` where it does not.

use crate::Point;

/// Why a coordinate could not be brought onto the `i16` lattice.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum CoordError {
    /// The value lies outside `i16::MIN..=i16::MAX`.
    OutOfRange,
    /// The value is not an integer, so it is not on the lattice.
    NotAnInteger,
    /// The value is NaN or infinite.
    NotFinite,
}

impl core::fmt::Display for CoordError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CoordError::OutOfRange => write!(f, "coordinate is outside the i16 range"),
            CoordError::NotAnInteger => write!(f, "coordinate is not an integer"),
            CoordError::NotFinite => write!(f, "coordinate is not finite"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for CoordError {}

/// Narrows one floating-point coordinate onto the lattice, exactly or not at
/// all.
fn lattice(v: f64) -> Result<i16, CoordError> {
    if !v.is_finite() {
        return Err(CoordError::NotFinite);
    }
    // `as i64` truncates, so comparing back catches any fractional part.
    let t = v as i64;
    if t as f64 != v {
        return Err(CoordError::NotAnInteger);
    }
    if t < i16::MIN as i64 || t > i16::MAX as i64 {
        return Err(CoordError::OutOfRange);
    }
    Ok(t as i16)
}

// --- geo-types --------------------------------------------------------------

#[cfg(feature = "geo-types")]
#[cfg_attr(docsrs, doc(cfg(feature = "geo-types")))]
mod geo_types_impl {
    use super::*;

    impl From<Point> for geo_types::Coord<i16> {
        #[inline]
        fn from(p: Point) -> Self {
            geo_types::Coord { x: p.x, y: p.y }
        }
    }

    impl From<geo_types::Coord<i16>> for Point {
        #[inline]
        fn from(c: geo_types::Coord<i16>) -> Self {
            Point::new(c.x, c.y)
        }
    }

    impl From<Point> for geo_types::Point<i16> {
        #[inline]
        fn from(p: Point) -> Self {
            geo_types::Point::new(p.x, p.y)
        }
    }

    impl From<geo_types::Point<i16>> for Point {
        #[inline]
        fn from(p: geo_types::Point<i16>) -> Self {
            Point::new(p.x(), p.y())
        }
    }

    /// Widening to `f64` is exact, which is what makes this `From` rather than
    /// `TryFrom`: `f64` has 53 mantissa bits and `i16` needs 16.
    impl From<Point> for geo_types::Coord<f64> {
        #[inline]
        fn from(p: Point) -> Self {
            geo_types::Coord {
                x: p.x as f64,
                y: p.y as f64,
            }
        }
    }

    impl TryFrom<geo_types::Coord<f64>> for Point {
        type Error = CoordError;
        #[inline]
        fn try_from(c: geo_types::Coord<f64>) -> Result<Self, CoordError> {
            Ok(Point::new(lattice(c.x)?, lattice(c.y)?))
        }
    }
}

// --- glam -------------------------------------------------------------------

#[cfg(feature = "glam")]
#[cfg_attr(docsrs, doc(cfg(feature = "glam")))]
mod glam_impl {
    use super::*;

    impl From<Point> for glam::I16Vec2 {
        #[inline]
        fn from(p: Point) -> Self {
            glam::I16Vec2::new(p.x, p.y)
        }
    }

    impl From<glam::I16Vec2> for Point {
        #[inline]
        fn from(v: glam::I16Vec2) -> Self {
            Point::new(v.x, v.y)
        }
    }

    impl From<Point> for glam::IVec2 {
        #[inline]
        fn from(p: Point) -> Self {
            glam::IVec2::new(p.x as i32, p.y as i32)
        }
    }

    impl TryFrom<glam::IVec2> for Point {
        type Error = CoordError;
        #[inline]
        fn try_from(v: glam::IVec2) -> Result<Self, CoordError> {
            let x = i16::try_from(v.x).map_err(|_| CoordError::OutOfRange)?;
            let y = i16::try_from(v.y).map_err(|_| CoordError::OutOfRange)?;
            Ok(Point::new(x, y))
        }
    }

    /// `i16 -> f32` is exact: 16 bits fits comfortably in a 24-bit mantissa.
    impl From<Point> for glam::Vec2 {
        #[inline]
        fn from(p: Point) -> Self {
            glam::Vec2::new(p.x as f32, p.y as f32)
        }
    }

    impl TryFrom<glam::Vec2> for Point {
        type Error = CoordError;
        #[inline]
        fn try_from(v: glam::Vec2) -> Result<Self, CoordError> {
            Ok(Point::new(lattice(v.x as f64)?, lattice(v.y as f64)?))
        }
    }
}

// --- mint -------------------------------------------------------------------

#[cfg(feature = "mint")]
#[cfg_attr(docsrs, doc(cfg(feature = "mint")))]
mod mint_impl {
    use super::*;

    impl From<Point> for mint::Point2<i16> {
        #[inline]
        fn from(p: Point) -> Self {
            mint::Point2 { x: p.x, y: p.y }
        }
    }

    impl From<mint::Point2<i16>> for Point {
        #[inline]
        fn from(p: mint::Point2<i16>) -> Self {
            Point::new(p.x, p.y)
        }
    }

    impl From<Point> for mint::Vector2<i16> {
        #[inline]
        fn from(p: Point) -> Self {
            mint::Vector2 { x: p.x, y: p.y }
        }
    }

    impl From<mint::Vector2<i16>> for Point {
        #[inline]
        fn from(v: mint::Vector2<i16>) -> Self {
            Point::new(v.x, v.y)
        }
    }

    impl From<Point> for mint::Point2<f32> {
        #[inline]
        fn from(p: Point) -> Self {
            mint::Point2 {
                x: p.x as f32,
                y: p.y as f32,
            }
        }
    }

    impl TryFrom<mint::Point2<f32>> for Point {
        type Error = CoordError;
        #[inline]
        fn try_from(p: mint::Point2<f32>) -> Result<Self, CoordError> {
            Ok(Point::new(lattice(p.x as f64)?, lattice(p.y as f64)?))
        }
    }
}

// --- num-traits -------------------------------------------------------------

#[cfg(feature = "num-traits")]
#[cfg_attr(docsrs, doc(cfg(feature = "num-traits")))]
mod num_traits_impl {
    use super::*;
    use num_traits::{NumCast, ToPrimitive};

    impl Point {
        /// Converts to any numeric type `num-traits` can cast to.
        ///
        /// Always succeeds for the usual targets, since every `i16` fits.
        ///
        /// # Examples
        ///
        /// ```
        /// use straight_skeleton::Point;
        ///
        /// let p = Point::new(3, -4);
        /// assert_eq!(p.cast::<f32>(), Some((3.0, -4.0)));
        /// assert_eq!(p.cast::<i64>(), Some((3, -4)));
        /// // -4 has no unsigned representation.
        /// assert_eq!(p.cast::<u8>(), None);
        /// ```
        pub fn cast<T: NumCast>(self) -> Option<(T, T)> {
            Some((T::from(self.x)?, T::from(self.y)?))
        }

        /// Builds a point from any numeric type, if both coordinates land
        /// **exactly** on the `i16` lattice.
        ///
        /// Returns `None` rather than rounding. Note this is stricter than
        /// `num_traits`' own `to_i16`, which truncates `0.5` to `0`; a silently
        /// moved vertex is the last thing a geometry crate should hand you.
        ///
        /// # Examples
        ///
        /// ```
        /// use straight_skeleton::Point;
        ///
        /// assert_eq!(Point::try_cast(3.0f64, -4.0f64), Some(Point::new(3, -4)));
        /// assert_eq!(Point::try_cast(1e9f64, 0.0f64), None);   // out of range
        /// assert_eq!(Point::try_cast(0.5f64, 0.0f64), None);   // off the lattice
        /// assert_eq!(Point::try_cast(3i64, -4i64), Some(Point::new(3, -4)));
        /// ```
        pub fn try_cast<T: ToPrimitive>(x: T, y: T) -> Option<Point> {
            // Via f64 rather than `to_i16`, so that `lattice` can reject a
            // fractional value instead of truncating it. Any input too large
            // for f64 to hold exactly is far outside i16 anyway, so it is
            // rejected as out of range regardless.
            Some(Point::new(
                lattice(x.to_f64()?).ok()?,
                lattice(y.to_f64()?).ok()?,
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lattice_accepts_exact_integers() {
        assert_eq!(lattice(0.0), Ok(0));
        assert_eq!(lattice(-32768.0), Ok(i16::MIN));
        assert_eq!(lattice(32767.0), Ok(i16::MAX));
        assert_eq!(lattice(-0.0), Ok(0));
    }

    #[test]
    fn lattice_rejects_rather_than_rounds() {
        assert_eq!(lattice(0.5), Err(CoordError::NotAnInteger));
        assert_eq!(lattice(-1.25), Err(CoordError::NotAnInteger));
    }

    #[test]
    fn lattice_rejects_rather_than_wraps() {
        assert_eq!(lattice(32768.0), Err(CoordError::OutOfRange));
        assert_eq!(lattice(-32769.0), Err(CoordError::OutOfRange));
        assert_eq!(lattice(1e18), Err(CoordError::OutOfRange));
    }

    #[test]
    fn lattice_rejects_nonfinite() {
        assert_eq!(lattice(f64::NAN), Err(CoordError::NotFinite));
        assert_eq!(lattice(f64::INFINITY), Err(CoordError::NotFinite));
        assert_eq!(lattice(f64::NEG_INFINITY), Err(CoordError::NotFinite));
    }

    #[cfg(feature = "geo-types")]
    #[test]
    fn geo_types_round_trips() {
        let p = Point::new(-7, 12);
        let c: geo_types::Coord<i16> = p.into();
        assert_eq!(Point::from(c), p);

        let gp: geo_types::Point<i16> = p.into();
        assert_eq!(Point::from(gp), p);

        let f: geo_types::Coord<f64> = p.into();
        assert_eq!(Point::try_from(f), Ok(p));

        assert_eq!(
            Point::try_from(geo_types::Coord { x: 0.5f64, y: 0.0 }),
            Err(CoordError::NotAnInteger)
        );
    }

    #[cfg(feature = "glam")]
    #[test]
    fn glam_round_trips() {
        let p = Point::new(-7, 12);
        assert_eq!(Point::from(glam::I16Vec2::from(p)), p);
        assert_eq!(Point::try_from(glam::IVec2::from(p)), Ok(p));
        assert_eq!(Point::try_from(glam::Vec2::from(p)), Ok(p));

        assert_eq!(
            Point::try_from(glam::IVec2::new(70_000, 0)),
            Err(CoordError::OutOfRange)
        );
        assert_eq!(
            Point::try_from(glam::Vec2::new(0.5, 0.0)),
            Err(CoordError::NotAnInteger)
        );
    }

    #[cfg(feature = "mint")]
    #[test]
    fn mint_round_trips() {
        let p = Point::new(-7, 12);
        assert_eq!(Point::from(mint::Point2::from(p)), p);
        assert_eq!(Point::from(mint::Vector2::from(p)), p);
        let f: mint::Point2<f32> = p.into();
        assert_eq!(Point::try_from(f), Ok(p));
    }

    #[cfg(feature = "num-traits")]
    #[test]
    fn num_traits_casts() {
        let p = Point::new(3, -4);
        assert_eq!(p.cast::<f64>(), Some((3.0, -4.0)));
        assert_eq!(p.cast::<i32>(), Some((3, -4)));
        assert_eq!(p.cast::<u16>(), None, "-4 is not representable");
        assert_eq!(Point::try_cast(3i32, -4i32), Some(p));
        assert_eq!(Point::try_cast(99_999i32, 0i32), None);
        // Stricter than num-traits' own to_i16, which would truncate to 0.
        assert_eq!(Point::try_cast(0.5f64, 0.0f64), None);
        assert_eq!(Point::try_cast(-0.5f32, 0.0f32), None);
        assert_eq!(Point::try_cast(f64::NAN, 0.0), None);
    }
}
