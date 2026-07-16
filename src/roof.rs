//! Hip roofs, raised from a straight skeleton.
//!
//! This is the classic application of a straight skeleton, and it is almost too
//! neat. Picture a roof being built by raising the walls' eaves inward at a
//! constant slope. At any moment the still-unroofed floor area is exactly the
//! shrinking wavefront, and the height reached is exactly how far it has
//! travelled. So the roof is not *computed* from the skeleton — it is *read off*
//! it:
//!
//! - every skeleton [`Node`] is a roof vertex, at height `offset * pitch`;
//! - every skeleton [face] is one flat roof panel, rising from the wall that
//!   face belongs to;
//! - every skeleton [`Arc`] is a hip, valley, or ridge, where two panels meet.
//!
//! ```text
//!      plan (a 2:1 rectangle)              roof, seen from the side
//!
//!    +-------------------+                        ______
//!    | \               / |                       /|    |\
//!    |   \___________/   |          =>          / |    | \
//!    |   /           \   |                     /  |    |  \
//!    | /               \ |                    /___|____|___\
//!    +-------------------+
//!         the skeleton                  the ridge is the skeleton's ridge,
//!      (four hips, one ridge)           at height ridge_offset * pitch
//! ```
//!
//! [`Node`]: crate::Node
//! [`Arc`]: crate::Arc
//! [face]: crate::Skeleton::face
//!
//! # Examples
//!
//! ```
//! use straight_skeleton::{skeleton, Point, Point3, Polygon, Roof};
//!
//! // A 120 x 80 floor plan.
//! let plan = Polygon::from_outer(&[
//!     Point::new(0, 0),
//!     Point::new(120, 0),
//!     Point::new(120, 80),
//!     Point::new(0, 80),
//! ])?;
//!
//! let roof = Roof::new(&skeleton(&plan)?, 0.5)?;
//!
//! // One panel per wall.
//! assert_eq!(roof.panels().len(), 4);
//!
//! // The ridge sits over the middle of the short side, 40 in, so at a pitch
//! // of 0.5 it is 20 high.
//! assert_eq!(roof.ridge_height(), 20);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use alloc::vec::Vec;
use core::fmt;

use crate::point::round_half_away_from_zero;
use crate::polygon::EdgeId;
use crate::skeleton::{NodeId, Skeleton};

/// A point on the 3D integer lattice.
///
/// The `z` axis is up. Like [`Point`], every coordinate is `i16`, so a roof is
/// `i16` in and `i16` out just as the rest of the crate is.
///
/// [`Point`]: crate::Point
///
/// # Examples
///
/// ```
/// use straight_skeleton::Point3;
///
/// let p = Point3::new(1, 2, 3);
/// assert_eq!(p.z, 3);
/// assert_eq!(<[i16; 3]>::from(p), [1, 2, 3]);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Point3 {
    /// Horizontal coordinate.
    pub x: i16,
    /// Depth coordinate.
    pub y: i16,
    /// Height above the eaves.
    pub z: i16,
}

impl Point3 {
    /// The origin.
    pub const ORIGIN: Point3 = Point3 { x: 0, y: 0, z: 0 };

    /// Constructs a point from its coordinates.
    #[inline]
    pub const fn new(x: i16, y: i16, z: i16) -> Self {
        Point3 { x, y, z }
    }
}

impl From<Point3> for [i16; 3] {
    #[inline]
    fn from(p: Point3) -> Self {
        [p.x, p.y, p.z]
    }
}

impl From<[i16; 3]> for Point3 {
    #[inline]
    fn from([x, y, z]: [i16; 3]) -> Self {
        Point3::new(x, y, z)
    }
}

impl From<Point3> for (i16, i16, i16) {
    #[inline]
    fn from(p: Point3) -> Self {
        (p.x, p.y, p.z)
    }
}

/// Why a roof could not be raised.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum RoofError {
    /// The pitch was negative, NaN, or infinite.
    InvalidPitch {
        /// The value supplied.
        pitch: f32,
    },
    /// A wall's face could not be walked, so its panel has no outline.
    ///
    /// The usual cause is handing in a **constrained** skeleton: limits
    /// truncate it into disconnected stubs, so its faces are not closed
    /// regions and there is no panel to raise. Roofs need a plain
    /// [`skeleton`].
    ///
    /// [`skeleton`]: crate::skeleton
    UnwalkableFace {
        /// The wall whose panel could not be outlined.
        wall: EdgeId,
    },
    /// A vertex would stand higher than `i16` can hold.
    ///
    /// The plan is too wide for this pitch. Lower the pitch, or scale the plan
    /// down. This is reported rather than saturated: a silently flattened ridge
    /// is a wrong roof, not an approximate one.
    HeightOverflow {
        /// The offending node.
        node: NodeId,
        /// The height it wanted.
        height: f32,
    },
}

impl fmt::Display for RoofError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RoofError::InvalidPitch { pitch } => {
                write!(f, "pitch {pitch} is not a finite, non-negative number")
            }
            RoofError::UnwalkableFace { wall } => write!(
                f,
                "could not walk the face of wall {}; roofs need a plain skeleton, \
                 not a constrained one",
                wall.0
            ),
            RoofError::HeightOverflow { node, height } => write!(
                f,
                "node {} would stand {height} high, which overflows i16; \
                 lower the pitch or scale the plan down",
                node.0
            ),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for RoofError {}

/// One corner of a roof.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RoofVertex {
    /// The corner's position, rounded to the integer lattice.
    ///
    /// `x` and `y` come from the skeleton node; `z` is `offset * pitch`.
    pub position: Point3,
    /// The corner's unrounded position.
    ///
    /// Rounding `z` to the lattice tilts each panel very slightly, so a roof
    /// built from [`RoofVertex::position`] is planar only to within half a
    /// unit. Use this when exact planarity matters — a renderer computing
    /// normals, say. It mirrors [`Node::exact`], for the same reason.
    ///
    /// [`Node::exact`]: crate::Node::exact
    pub exact: [f32; 3],
    /// The skeleton node this corner was raised from.
    pub node: NodeId,
}

/// One flat plane of roof, rising from a single wall.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Panel {
    /// The wall this panel rises from.
    ///
    /// This is the traceability a straight skeleton gives for free: no search,
    /// no nearest-neighbour query. The panel *is* the region that wall's
    /// wavefront swept.
    pub wall: EdgeId,
    /// The panel's corners, in order, as indices into [`Roof::verts`].
    ///
    /// The first two are the wall's own endpoints, so the outline starts along
    /// the eave and works back over the roof.
    pub corners: Vec<NodeId>,
}

/// A hip roof over a floor plan.
///
/// Build one with [`Roof::new`]. Vertices are indexed by [`NodeId`]: roof
/// vertex `i` stands directly over skeleton node `i`, so the two structures
/// stay in step and provenance survives into 3D.
///
/// # Why every panel is flat
///
/// A panel is the region swept by one wall's wavefront, so every point on it is
/// `offset` away from that wall's supporting line. Height is `offset * pitch` —
/// a linear function of that distance — so the panel is a plane. That is
/// guaranteed by construction rather than fitted, and the crate's tests check
/// it holds.
///
/// # Examples
///
/// ```
/// use straight_skeleton::{skeleton, Point, Polygon, Roof};
///
/// // An L-shaped house.
/// let plan = Polygon::from_outer(&[
///     Point::new(0, 0),
///     Point::new(160, 0),
///     Point::new(160, 70),
///     Point::new(70, 70),
///     Point::new(70, 150),
///     Point::new(0, 150),
/// ])?;
///
/// let roof = Roof::new(&skeleton(&plan)?, 0.6)?;
///
/// // Six walls, six panels.
/// assert_eq!(roof.panels().len(), 6);
///
/// // Every panel knows which wall it rises from.
/// for (i, panel) in roof.panels().iter().enumerate() {
///     assert_eq!(panel.wall.0 as usize, i);
/// }
///
/// // Eaves sit at zero; nothing is below them.
/// assert!(roof.verts().iter().all(|v| v.position.z >= 0));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Roof {
    verts: Vec<RoofVertex>,
    panels: Vec<Panel>,
    pitch: f32,
}

impl Roof {
    /// Raises a roof over a skeleton.
    ///
    /// `pitch` is rise over run: 1.0 gives 45°, 0.5 gives a shallower roof
    /// half as tall, 0.0 gives a flat one.
    ///
    /// # Errors
    ///
    /// - [`RoofError::InvalidPitch`] if `pitch` is negative, NaN, or infinite.
    /// - [`RoofError::UnwalkableFace`] if a wall's face is not a closed region.
    ///   In practice this means the skeleton came from
    ///   [`skeleton_constrained`], which truncates it into disconnected stubs;
    ///   roofs need a plain [`skeleton`].
    /// - [`RoofError::HeightOverflow`] if the plan is too wide for the pitch to
    ///   fit in `i16`.
    ///
    /// [`skeleton`]: crate::skeleton
    /// [`skeleton_constrained`]: crate::skeleton_constrained
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::{skeleton, skeleton_constrained, Point, Polygon, Roof, RoofError};
    ///
    /// let plan = Polygon::from_outer(&[
    ///     Point::new(0, 0), Point::new(80, 0), Point::new(80, 80), Point::new(0, 80),
    /// ])?;
    /// let skel = skeleton(&plan)?;
    ///
    /// // A square plan gives a pyramid: its apex is 40 in from every wall.
    /// assert_eq!(Roof::new(&skel, 1.0)?.ridge_height(), 40);
    /// assert_eq!(Roof::new(&skel, 0.5)?.ridge_height(), 20);
    /// assert_eq!(Roof::new(&skel, 0.0)?.ridge_height(), 0);
    ///
    /// // A pitch that would push the apex past i16 is refused, not clamped.
    /// assert!(matches!(
    ///     Roof::new(&skel, 1000.0),
    ///     Err(RoofError::HeightOverflow { .. })
    /// ));
    ///
    /// // Nor can a truncated skeleton be roofed: its faces are not closed.
    /// let truncated = skeleton_constrained(&plan, &[5.0; 4])?;
    /// assert!(matches!(
    ///     Roof::new(&truncated, 1.0),
    ///     Err(RoofError::UnwalkableFace { .. })
    /// ));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(skeleton: &Skeleton, pitch: f32) -> Result<Roof, RoofError> {
        if !pitch.is_finite() || pitch < 0.0 {
            return Err(RoofError::InvalidPitch { pitch });
        }

        // A skeleton node at offset d becomes a roof vertex at height d * pitch:
        // the wavefront's travel *is* the run, so the rise follows directly.
        let mut verts = Vec::with_capacity(skeleton.node_count());
        for (i, node) in skeleton.nodes().iter().enumerate() {
            let id = NodeId(i as u32);
            let height = node.offset * pitch;
            let z = lattice_height(height).ok_or(RoofError::HeightOverflow { node: id, height })?;
            verts.push(RoofVertex {
                position: Point3::new(node.position.x, node.position.y, z),
                exact: [node.exact[0], node.exact[1], height],
                node: id,
            });
        }

        let mut panels = Vec::with_capacity(skeleton.input_edge_count());
        for i in 0..skeleton.input_edge_count() as u16 {
            let wall = EdgeId(i);
            let corners = skeleton
                .face(wall)
                .ok_or(RoofError::UnwalkableFace { wall })?;
            panels.push(Panel { wall, corners });
        }

        Ok(Roof {
            verts,
            panels,
            pitch,
        })
    }

    /// Every corner of the roof, indexed by [`NodeId`].
    #[inline]
    pub fn verts(&self) -> &[RoofVertex] {
        &self.verts
    }

    /// Every panel, one per wall, in wall order.
    #[inline]
    pub fn panels(&self) -> &[Panel] {
        &self.panels
    }

    /// The pitch this roof was raised at.
    #[inline]
    pub fn pitch(&self) -> f32 {
        self.pitch
    }

    /// The corner standing over a given skeleton node.
    ///
    /// # Panics
    ///
    /// Panics if `n` does not belong to the skeleton this roof was built from.
    #[inline]
    pub fn vertex(&self, n: NodeId) -> &RoofVertex {
        &self.verts[n.0 as usize]
    }

    /// The panel rising from a given wall.
    ///
    /// # Panics
    ///
    /// Panics if `wall` does not belong to the skeleton this roof was built
    /// from.
    #[inline]
    pub fn panel(&self, wall: EdgeId) -> &Panel {
        &self.panels[wall.0 as usize]
    }

    /// The height of the highest point: the ridge, or a pyramid's apex.
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::{skeleton, Point, Polygon, Roof};
    ///
    /// // A 120-wide, 80-deep plan. The ridge runs down the middle of the long
    /// // axis, 40 in from each long wall, so at pitch 0.5 it stands 20 high.
    /// let plan = Polygon::from_outer(&[
    ///     Point::new(0, 0), Point::new(120, 0), Point::new(120, 80), Point::new(0, 80),
    /// ])?;
    /// assert_eq!(Roof::new(&skeleton(&plan)?, 0.5)?.ridge_height(), 20);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn ridge_height(&self) -> i16 {
        self.verts.iter().map(|v| v.position.z).max().unwrap_or(0)
    }

    /// The corners of one panel, as positions.
    ///
    /// # Panics
    ///
    /// Panics if `wall` does not belong to the skeleton this roof was built
    /// from.
    pub fn panel_outline(&self, wall: EdgeId) -> Vec<Point3> {
        self.panel(wall)
            .corners
            .iter()
            .map(|&n| self.vertex(n).position)
            .collect()
    }
}

/// Rounds a height onto the lattice, or `None` if it will not fit.
///
/// Unlike [`crate::Point`]'s rounding, this refuses rather than saturates. A
/// saturated `x`/`y` is a node nudged by a hair at the very edge of the
/// coordinate space; a saturated `z` is a roof with its ridge sliced flat, and
/// the caller needs to know.
fn lattice_height(h: f32) -> Option<i16> {
    if !h.is_finite() {
        return None;
    }
    let r = round_half_away_from_zero(h);
    if r < i16::MIN as f32 || r > i16::MAX as f32 {
        return None;
    }
    Some(r as i16)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lattice_height_rounds_to_nearest() {
        assert_eq!(lattice_height(0.0), Some(0));
        assert_eq!(lattice_height(2.4), Some(2));
        assert_eq!(lattice_height(2.5), Some(3));
        assert_eq!(lattice_height(2.6), Some(3));
    }

    #[test]
    fn lattice_height_refuses_rather_than_saturating() {
        assert_eq!(lattice_height(32767.0), Some(i16::MAX));
        assert_eq!(lattice_height(32768.0), None);
        assert_eq!(lattice_height(1e9), None);
        assert_eq!(lattice_height(f32::INFINITY), None);
        assert_eq!(lattice_height(f32::NAN), None);
    }

    #[test]
    fn point3_conversions() {
        let p = Point3::new(1, -2, 3);
        assert_eq!(<[i16; 3]>::from(p), [1, -2, 3]);
        assert_eq!(Point3::from([1, -2, 3]), p);
        assert_eq!(<(i16, i16, i16)>::from(p), (1, -2, 3));
        assert_eq!(Point3::ORIGIN, Point3::new(0, 0, 0));
    }
}
