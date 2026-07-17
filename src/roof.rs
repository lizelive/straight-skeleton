//! Roofs, raised from a straight skeleton.
//!
//! This is the classic application of a straight skeleton, and it is almost too
//! neat. Picture a roof being built by raising the walls' eaves inward at a
//! constant slope. At any moment the still-unroofed floor area is exactly the
//! shrinking wavefront, and the height reached is exactly how far it has
//! travelled. So the roof is not *computed* from the skeleton — it is *read off*
//! it:
//!
//! - every skeleton [`Node`] is a roof vertex, at height [`Profile::height_at`]
//!   its offset;
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
//! # The skeleton is the plan, not the roof
//!
//! Worth separating, because it is what makes everything below cheap. A
//! straight skeleton says where the hips, valleys and ridges *run*. It says
//! nothing about how high anything is — height is a function of one variable,
//! [`Node::offset`], and that function is the [`Profile`].
//!
//! So the roof styles here are not different algorithms. They are the same
//! skeleton read with a different profile:
//!
//! ```text
//!     hip                mansard              truncated         truncated
//!                                                                mansard
//!        /\               ___                  _____              _____
//!       /  \             /   \                /     \            /     \
//!      /    \           |     |              /       \          |       |
//!     /______\          |_____|             /_________\         |_______|
//!
//!     one pitch       steep, then          stopped short,     all three at
//!                      shallow             leaving a flat        once
//! ```
//!
//! The last two need a [`skeleton_constrained`], whose wavefront stops rather
//! than collapsing; the flat is that [residual] raised to the limit's height.
//! See [`Roof::with_profile`].
//!
//! [`Node`]: crate::Node
//! [`Node::offset`]: crate::Node::offset
//! [`Arc`]: crate::Arc
//! [face]: crate::Skeleton::face
//! [`skeleton_constrained`]: crate::skeleton_constrained
//! [residual]: crate::Skeleton::residual
//!
//! # Examples
//!
//! ```
//! use straight_skeleton::{skeleton, skeleton_constrained, Point, Polygon, Roof};
//!
//! // A 120 x 80 floor plan.
//! let plan = Polygon::from_outer(&[
//!     Point::new(0, 0),
//!     Point::new(120, 0),
//!     Point::new(120, 80),
//!     Point::new(0, 80),
//! ])?;
//! let skel = skeleton(&plan)?;
//!
//! // A hip roof: one panel per wall. The ridge sits over the middle of the
//! // short side, 40 in, so at a pitch of 0.5 it is 20 high.
//! let hip = Roof::new(&skel, 0.5)?;
//! assert_eq!(hip.panels().len(), 4);
//! assert_eq!(hip.ridge_height(), 20);
//!
//! // A mansard over the *same* skeleton: steep to 10, then shallow. Two
//! // panels per wall now, because the break cuts each one in half.
//! let mansard = Roof::mansard(&skel, 2.0, 10.0, 0.25)?;
//! assert_eq!(mansard.panels().len(), 8);
//! assert_eq!(mansard.ridge_height(), 28);   // 10 * 2 + 30 * 0.25
//!
//! // Stop every wall at 15 and the apex is cut off, leaving a flat on top.
//! let truncated = Roof::new(&skeleton_constrained(&plan, &[15.0; 4])?, 0.5)?;
//! assert_eq!(truncated.flat().count(), 1);
//! assert_eq!(truncated.ridge_height(), 8);  // 15 * 0.5, rounded
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use alloc::collections::BTreeMap;
use alloc::vec::Vec;
use core::fmt;

use crate::math::Vec2;
use crate::point::round_half_away_from_zero;
use crate::polygon::EdgeId;
use crate::skeleton::{NodeId, NodeKind, Skeleton};
use crate::Point;

/// Tolerance for deciding which side of a profile's break a corner sits on, in
/// coordinate units.
///
/// Comfortably above `f32`'s resolution at the top of the coordinate range
/// (~0.002) and far below one lattice unit, so a corner is only called "on the
/// break" when it genuinely is.
const BREAK_EPS: f32 = 1e-2;

/// How a roof's height grows with distance from the eaves.
///
/// A straight skeleton gives the *plan* of a roof: where the hips, valleys and
/// ridges run. It does not say how high anything is — that is this. Height is a
/// function of one variable, [`Node::offset`], because a skeleton node's offset
/// is how far the wavefront had travelled to reach it, which is exactly the run
/// the roof has had to rise over.
///
/// So changing the roof *style* does not change the skeleton at all. A mansard
/// over a plan has the same hips and ridge as a hip roof over it; only `z`
/// differs.
///
/// [`Node::offset`]: crate::Node::offset
///
/// # Examples
///
/// ```
/// use straight_skeleton::Profile;
///
/// // A hip: one slope all the way up.
/// let hip = Profile::Hip { pitch: 0.5 };
/// assert_eq!(hip.height_at(10.0), 5.0);
///
/// // A mansard: steep to 10, then shallow.
/// let mansard = Profile::Mansard {
///     lower_pitch: 2.0,
///     break_offset: 10.0,
///     upper_pitch: 0.25,
/// };
/// assert_eq!(mansard.height_at(10.0), 20.0);       // the break
/// assert_eq!(mansard.height_at(30.0), 25.0);       // 20 + 20 * 0.25
/// ```
/// Exhaustive on purpose, unlike the crate's error types: this is a value
/// callers are meant to match on, and a roof style they cannot handle is better
/// caught by the compiler than by a `_` arm quietly treating it as something
/// else.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Profile {
    /// One slope all the way to the ridge: `z = offset * pitch`.
    ///
    /// The classic hip roof.
    Hip {
        /// Rise over run. 1.0 gives 45°, 0.0 gives a flat roof.
        pitch: f32,
    },
    /// Two slopes with a break between them: a **mansard**.
    ///
    /// Steep from the eaves up to `break_offset`, shallow from there on. That
    /// is what a mansard is for — the steep lower slope buys headroom in the
    /// storey inside it, and the shallow upper one keeps the whole thing from
    /// becoming absurdly tall.
    ///
    /// The break is at a constant *offset*, which is a constant *height*
    /// (`break_offset * lower_pitch`), so it comes out as a level line all the
    /// way round the roof — the kerb a real mansard has.
    ///
    /// Nothing stops `upper_pitch` being the steeper of the two; the type does
    /// not police taste. Equal pitches reduce to a [`Profile::Hip`].
    Mansard {
        /// Rise over run below the break.
        lower_pitch: f32,
        /// How far from the eaves the pitch changes, in plan units.
        break_offset: f32,
        /// Rise over run above the break.
        upper_pitch: f32,
    },
}

impl Profile {
    /// The height at a given distance from the eaves.
    #[inline]
    pub fn height_at(self, offset: f32) -> f32 {
        match self {
            Profile::Hip { pitch } => offset * pitch,
            Profile::Mansard {
                lower_pitch,
                break_offset,
                upper_pitch,
            } => {
                if offset <= break_offset {
                    offset * lower_pitch
                } else {
                    break_offset * lower_pitch + (offset - break_offset) * upper_pitch
                }
            }
        }
    }

    /// The offset the slope changes at, if it changes at all.
    ///
    /// This is where panels have to be cut in two: a panel spanning it would be
    /// bent rather than flat.
    #[inline]
    fn break_offset(self) -> Option<f32> {
        match self {
            Profile::Hip { .. } => None,
            Profile::Mansard { break_offset, .. } => Some(break_offset),
        }
    }

    fn validate(self) -> Result<(), RoofError> {
        let bad = |p: f32| !p.is_finite() || p < 0.0;
        match self {
            Profile::Hip { pitch } => {
                if bad(pitch) {
                    return Err(RoofError::InvalidPitch { pitch });
                }
            }
            Profile::Mansard {
                lower_pitch,
                break_offset,
                upper_pitch,
            } => {
                if bad(lower_pitch) {
                    return Err(RoofError::InvalidPitch { pitch: lower_pitch });
                }
                if bad(upper_pitch) {
                    return Err(RoofError::InvalidPitch { pitch: upper_pitch });
                }
                if !break_offset.is_finite() || break_offset < 0.0 {
                    return Err(RoofError::InvalidBreak { break_offset });
                }
            }
        }
        Ok(())
    }
}

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
    /// Should not happen for a skeleton of a valid polygon, constrained or not.
    UnwalkableFace {
        /// The wall whose panel could not be outlined.
        wall: EdgeId,
    },
    /// A [`Profile::Mansard`]'s break was negative, NaN, or infinite.
    InvalidBreak {
        /// The value supplied.
        break_offset: f32,
    },
    /// The skeleton's edges stopped **partway**, at different distances, so
    /// there is no roof over it.
    ///
    /// A roof's height is a function of [`Node::offset`], and that only works
    /// while offset means *distance from the wall*. On a plain skeleton it
    /// always does. On a [`skeleton_constrained`] it is really the wavefront's
    /// **time**, which is the same thing only until something stops early: an
    /// edge that halted at 3 is 3 away from its face forever, however long the
    /// clock runs on.
    ///
    /// What is allowed follows from that, and it is more than it first looks:
    ///
    /// - **No limits**, or limits that never bind. A plain hip roof.
    /// - **One uniform limit.** Every edge stops together, at the top, so
    ///   nothing stops *early* and the roof is a hip truncated to a flat.
    /// - **A limit of zero**, mixed freely with either of the above. A wall that
    ///   never moves sweeps nothing, so it has no sloping panel to be
    ///   inconsistent about: its face is degenerate in plan, and stands up as a
    ///   vertical **gable**. The neighbouring walls' corners slide along it, so
    ///   the ridge runs out to it. See [`PanelKind::Slope`].
    ///
    /// What is refused is an edge stopping at some distance *in between*, while
    /// others go further. That edge has a real sloping panel, and its far
    /// corners would sit at their offset's height while being only `limit` from
    /// the wall — so the panel would want to end lower than its neighbour's, and
    /// the surface between them would have to tear. There is no roof to return,
    /// so this says so rather than inventing one.
    ///
    /// [`Node::offset`]: crate::Node::offset
    /// [`skeleton_constrained`]: crate::skeleton_constrained
    UnevenLimits {
        /// A node where an edge stopped before the rest of the roof did.
        node: NodeId,
        /// The offset it stopped at.
        stopped_at: f32,
        /// The offset the rest of the roof reaches.
        reaches: f32,
    },
    /// A [`Profile::Mansard`]'s break left a wall's face in pieces that could
    /// not be closed back up into coherent panels.
    ///
    /// The break is a straight line across a face, so it normally splits it into
    /// a lower piece and an upper one — and a face whose wavefront was split
    /// apart by two reflex vertices into several lower and upper pieces is
    /// handled too, each piece emitted as its own panel. This is the safety net
    /// for the geometry that is left: a face the break meets in a way that does
    /// not resolve into closed pieces, which should not arise for the face of a
    /// valid polygon. Refused rather than returned bent or self-touching.
    BreakSplitsPanel {
        /// The wall whose panel the break could not cleanly cut.
        wall: EdgeId,
        /// How many loose ends the break left on the face's outline. An odd
        /// count is the tell — a closed outline must cross a line an even number
        /// of times.
        crossings: usize,
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
            RoofError::UnwalkableFace { wall } => {
                write!(f, "could not walk the face of wall {}", wall.0)
            }
            RoofError::InvalidBreak { break_offset } => write!(
                f,
                "mansard break {break_offset} is not a finite, non-negative distance"
            ),
            RoofError::UnevenLimits {
                node,
                stopped_at,
                reaches,
            } => write!(
                f,
                "an edge stopped at {stopped_at} (node {}) while the rest of the roof \
                 reaches {reaches}; a roof needs every limit to be the same, or none",
                node.0
            ),
            RoofError::BreakSplitsPanel { wall, crossings } => write!(
                f,
                "the mansard break left wall {}'s face with {crossings} loose ends, \
                 so it could not be cut into coherent panels",
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

/// Identifies a corner of a [`Roof`].
///
/// The first [`Skeleton::node_count`] of these stand directly over the skeleton
/// nodes of the same number, so `RoofVertexId(i)` is `NodeId(i)` raised — which
/// is what keeps provenance alive into 3D. A [`Profile::Mansard`] appends more
/// beyond that, where its break line cuts across a panel; those stand over no
/// skeleton node, and their [`RoofVertex::node`] is `None`.
///
/// [`Skeleton::node_count`]: crate::Skeleton::node_count
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RoofVertexId(pub u32);

/// One corner of a roof.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RoofVertex {
    /// The corner's position, rounded to the integer lattice.
    ///
    /// `z` is [`Profile::height_at`] the corner's offset.
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
    ///
    /// `None` for a corner a [`Profile::Mansard`]'s break line introduced,
    /// which lies partway along a skeleton arc rather than at either end of it.
    pub node: Option<NodeId>,
}

/// What part of a roof a [`Panel`] is.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PanelKind {
    /// A panel rising from one wall.
    ///
    /// # Gables
    ///
    /// Usually sloping, but not always. A wall given a limit of **zero** never
    /// moves, so it sweeps no plan area at all, and its face is degenerate —
    /// every corner sits on the wall's own line. Stood up at `z =
    /// height_at(offset)` that degenerate face becomes a vertical triangle: a
    /// **gable**. The neighbouring walls' corners slide along the frozen wall
    /// rather than over it, so the ridge runs out to it instead of hipping away.
    ///
    /// It is the same panel with the same rule applied, so it gets no special
    /// variant. A consumer that wants to tell the two apart can: a gable's
    /// footprint has zero area, and nothing else's does. The `roof` example
    /// names its OBJ groups that way.
    Slope {
        /// The wall it rises from.
        ///
        /// This is the traceability a straight skeleton gives for free: no
        /// search, no nearest-neighbour query. The panel *is* the region that
        /// wall's wavefront swept.
        wall: EdgeId,
        /// Which band of the [`Profile`] it belongs to, counting from the
        /// eaves. A [`Profile::Hip`] only ever has band 0; a
        /// [`Profile::Mansard`] has band 0 below the break and band 1 above.
        band: u8,
    },
    /// The flat a truncated roof stops at, standing over the skeleton's
    /// [residual].
    ///
    /// One of these per residual loop, wound like the input's rings:
    /// counter-clockwise for a flat's outline, clockwise for a hole in one. So
    /// a flat with a hole in it — a courtyard's, say — is two `Flat` panels, and
    /// the winding is what says which is which.
    ///
    /// That is deliberately the raw material rather than a finished mesh. A
    /// consumer that cannot express holes has to triangulate, and *how* is its
    /// own business: the crate hands over the loops and their winding, which is
    /// the part only it knows. The `roof` example shows one way, in about thirty
    /// lines.
    ///
    /// [residual]: crate::Skeleton::residual
    Flat,
}

/// One flat plane of roof.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Panel {
    /// What part of the roof this is.
    pub kind: PanelKind,
    /// The panel's corners, in order, as indices into [`Roof::verts`].
    pub corners: Vec<RoofVertexId>,
}

impl Panel {
    /// The wall this panel rises from, or `None` if it is a flat.
    #[inline]
    pub fn wall(&self) -> Option<EdgeId> {
        match self.kind {
            PanelKind::Slope { wall, .. } => Some(wall),
            PanelKind::Flat => None,
        }
    }

    /// Whether this is the flat rather than a slope.
    #[inline]
    pub fn is_flat(&self) -> bool {
        matches!(self.kind, PanelKind::Flat)
    }
}

/// A roof over a floor plan.
///
/// Build one with [`Roof::new`] for a hip, [`Roof::mansard`] for two pitches, or
/// [`Roof::with_profile`] for either — and for the truncated versions of both,
/// over a constrained skeleton.
///
/// Vertices are indexed by [`RoofVertexId`], and the first
/// [`Skeleton::node_count`] of those stand directly over the skeleton nodes of
/// the same number, so the two structures stay in step and provenance survives
/// into 3D.
///
/// # Why every panel is flat
///
/// A panel is the region swept by one wall's wavefront, so every point on it is
/// `offset` away from that wall's supporting line — an affine function of
/// position. Each band of a [`Profile`] is affine in that distance, so the
/// composition is affine too, and the panel is a plane.
///
/// That is why a [`Profile::Mansard`]'s break has to *cut* panels rather than
/// merely bend them: the profile is only affine within a band, so a panel
/// spanning the break would be a fold, not a plane. Cut at the break, both
/// halves are planes again.
///
/// Planarity is therefore guaranteed by construction rather than fitted, and the
/// crate's tests re-derive every corner's height from its wall to check it.
///
/// [`Skeleton::node_count`]: crate::Skeleton::node_count
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
///     assert_eq!(panel.wall().unwrap().0 as usize, i);
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
    profile: Profile,
}

impl Roof {
    /// Raises a hip roof over a skeleton.
    ///
    /// `pitch` is rise over run: 1.0 gives 45°, 0.5 gives a shallower roof
    /// half as tall, 0.0 gives a flat one.
    ///
    /// Shorthand for [`Roof::with_profile`] with a [`Profile::Hip`]. Use
    /// [`Roof::mansard`] for two pitches with a break between them, and see
    /// [`Roof::with_profile`] for what a **constrained** skeleton does here.
    ///
    /// # Errors
    ///
    /// As [`Roof::with_profile`].
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
    /// // Stopping every wall at 10 cuts the apex off, leaving a flat.
    /// let truncated = skeleton_constrained(&plan, &[10.0; 4])?;
    /// let roof = Roof::new(&truncated, 1.0)?;
    /// assert_eq!(roof.ridge_height(), 10);
    /// assert_eq!(roof.flat().count(), 1);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(skeleton: &Skeleton, pitch: f32) -> Result<Roof, RoofError> {
        Roof::with_profile(skeleton, Profile::Hip { pitch })
    }

    /// Raises a **mansard** roof: steep to `break_offset`, shallow above it.
    ///
    /// Shorthand for [`Roof::with_profile`] with a [`Profile::Mansard`]. See
    /// there for what a mansard is and why the skeleton underneath is the same
    /// one a hip roof uses.
    ///
    /// # Errors
    ///
    /// As [`Roof::with_profile`].
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::{skeleton, Point, Polygon, Roof};
    ///
    /// // A 120 x 80 plan. Its ridge is 40 in from the long walls.
    /// let plan = Polygon::from_outer(&[
    ///     Point::new(0, 0), Point::new(120, 0), Point::new(120, 80), Point::new(0, 80),
    /// ])?;
    /// let skel = skeleton(&plan)?;
    ///
    /// // Steep (2:1) for the first 10, then shallow (1:4) to the ridge.
    /// let roof = Roof::mansard(&skel, 2.0, 10.0, 0.25)?;
    ///
    /// // 10 * 2 = 20 at the kerb, then 30 more of run at 0.25 = 7.5 -> 28.
    /// assert_eq!(roof.ridge_height(), 28);
    ///
    /// // Each of the four walls now carries two panels rather than one: the
    /// // steep skirt, and the shallow slope above it.
    /// assert_eq!(roof.panels().len(), 8);
    /// assert_eq!(roof.panels_of(straight_skeleton::EdgeId(0)).count(), 2);
    ///
    /// // A hip roof of the same plan is much taller for the same lower pitch.
    /// assert_eq!(Roof::new(&skel, 2.0)?.ridge_height(), 80);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn mansard(
        skeleton: &Skeleton,
        lower_pitch: f32,
        break_offset: f32,
        upper_pitch: f32,
    ) -> Result<Roof, RoofError> {
        Roof::with_profile(
            skeleton,
            Profile::Mansard {
                lower_pitch,
                break_offset,
                upper_pitch,
            },
        )
    }

    /// Raises a roof over a skeleton, with any [`Profile`].
    ///
    /// # Constrained skeletons
    ///
    /// A [`skeleton_constrained`] with one **uniform** limit gives a *truncated*
    /// roof: the slopes rise to the limit and stop, and the [residual] the
    /// wavefront stopped as becomes a [`PanelKind::Flat`] on top. That works
    /// with any profile, so a truncated mansard is steep, then shallow, then
    /// flat.
    ///
    /// **Uneven** limits are refused — see [`RoofError::UnevenLimits`], which
    /// explains why there is no such roof rather than merely no implementation.
    ///
    /// # Errors
    ///
    /// - [`RoofError::InvalidPitch`] if a pitch is negative, NaN, or infinite.
    /// - [`RoofError::InvalidBreak`] likewise for a mansard's break.
    /// - [`RoofError::UnevenLimits`] if the skeleton's edges stopped at
    ///   different distances.
    /// - [`RoofError::BreakSplitsPanel`] if a mansard's break leaves a face in
    ///   pieces that cannot be closed back up (a safety net that should not fire
    ///   for the face of a valid polygon).
    /// - [`RoofError::HeightOverflow`] if the plan is too wide for the pitch to
    ///   fit in `i16`.
    /// - [`RoofError::UnwalkableFace`] if a wall's face is not a closed region,
    ///   which should not happen for a skeleton of a valid polygon.
    ///
    /// [`skeleton_constrained`]: crate::skeleton_constrained
    /// [residual]: crate::Skeleton::residual
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::{skeleton_constrained, PanelKind, Point, Polygon, Profile, Roof};
    ///
    /// let plan = Polygon::from_outer(&[
    ///     Point::new(0, 0), Point::new(100, 0), Point::new(100, 100), Point::new(0, 100),
    /// ])?;
    ///
    /// // Every wall stopped at 20: a hip roof with its apex cut off.
    /// let skel = skeleton_constrained(&plan, &[20.0; 4])?;
    /// let roof = Roof::with_profile(&skel, Profile::Hip { pitch: 0.5 })?;
    ///
    /// // Four slopes, and the flat they stop at.
    /// assert_eq!(roof.panels().len(), 5);
    /// assert_eq!(roof.flat().count(), 1);
    ///
    /// // The flat stands at 20 * 0.5, and that is the top of the roof.
    /// assert_eq!(roof.ridge_height(), 10);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn with_profile(skeleton: &Skeleton, profile: Profile) -> Result<Roof, RoofError> {
        profile.validate()?;
        check_limits_are_even(skeleton)?;

        let mut build = Build {
            skel: skeleton,
            profile,
            verts: Vec::with_capacity(skeleton.node_count()),
            cuts: BTreeMap::new(),
        };

        // Every skeleton node becomes a roof vertex, in step, so that
        // `RoofVertexId(i)` is `NodeId(i)` raised. The wavefront's travel *is*
        // the run, so the rise follows from the offset alone.
        for (i, node) in skeleton.nodes().iter().enumerate() {
            let id = NodeId(i as u32);
            let height = profile.height_at(node.offset);
            let z = lattice_height(height).ok_or(RoofError::HeightOverflow { node: id, height })?;
            build.verts.push(RoofVertex {
                position: Point3::new(node.position.x, node.position.y, z),
                exact: [node.exact[0], node.exact[1], height],
                node: Some(id),
            });
        }

        let mut panels = Vec::with_capacity(skeleton.input_edge_count());
        for i in 0..skeleton.input_edge_count() as u16 {
            let wall = EdgeId(i);
            let face = skeleton
                .face(wall)
                .ok_or(RoofError::UnwalkableFace { wall })?;
            build.push_slopes(wall, &face, &mut panels)?;
        }

        // The flat, where the wavefront stopped rather than collapsing. It sits
        // at one height throughout, so it needs no cutting whatever the profile.
        for loop_ in skeleton.residual() {
            panels.push(Panel {
                kind: PanelKind::Flat,
                corners: loop_.nodes.iter().map(|&n| RoofVertexId(n.0)).collect(),
            });
        }

        Ok(Roof {
            verts: build.verts,
            panels,
            profile,
        })
    }

    /// Every corner of the roof, indexed by [`RoofVertexId`].
    #[inline]
    pub fn verts(&self) -> &[RoofVertex] {
        &self.verts
    }

    /// Every panel: the slopes in wall order, then any flats.
    #[inline]
    pub fn panels(&self) -> &[Panel] {
        &self.panels
    }

    /// The profile this roof was raised with.
    #[inline]
    pub fn profile(&self) -> Profile {
        self.profile
    }

    /// A corner of the roof.
    ///
    /// # Panics
    ///
    /// Panics if `v` does not belong to this roof.
    #[inline]
    pub fn vertex(&self, v: RoofVertexId) -> &RoofVertex {
        &self.verts[v.0 as usize]
    }

    /// The corner standing over a given skeleton node.
    ///
    /// # Panics
    ///
    /// Panics if `n` does not belong to the skeleton this roof was built from.
    #[inline]
    pub fn vertex_at(&self, n: NodeId) -> &RoofVertex {
        &self.verts[n.0 as usize]
    }

    /// The panels rising from a given wall, from the eaves up.
    ///
    /// A [`Profile::Hip`] gives exactly one; a [`Profile::Mansard`] gives two
    /// where its break crosses the panel, and one where it does not reach.
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::{skeleton, EdgeId, Point, Polygon, Roof};
    ///
    /// let plan = Polygon::from_outer(&[
    ///     Point::new(0, 0), Point::new(120, 0), Point::new(120, 80), Point::new(0, 80),
    /// ])?;
    /// let skel = skeleton(&plan)?;
    ///
    /// assert_eq!(Roof::new(&skel, 0.5)?.panels_of(EdgeId(0)).count(), 1);
    /// assert_eq!(Roof::mansard(&skel, 2.0, 10.0, 0.25)?.panels_of(EdgeId(0)).count(), 2);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn panels_of(&self, wall: EdgeId) -> impl Iterator<Item = &Panel> + '_ {
        self.panels.iter().filter(move |p| p.wall() == Some(wall))
    }

    /// The flat panels, if this roof is truncated. Empty otherwise.
    ///
    /// More than one only when the flat has a hole in it — see
    /// [`PanelKind::Flat`].
    pub fn flat(&self) -> impl Iterator<Item = &Panel> + '_ {
        self.panels.iter().filter(|p| p.is_flat())
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
    /// Panics if `panel` does not belong to this roof.
    pub fn outline(&self, panel: &Panel) -> Vec<Point3> {
        panel
            .corners
            .iter()
            .map(|&v| self.vertex(v).position)
            .collect()
    }
}

/// Every edge either never stopped, or stopped level with the top of the roof.
///
/// See [`RoofError::UnevenLimits`] for why anything else has no roof at all.
fn check_limits_are_even(skeleton: &Skeleton) -> Result<(), RoofError> {
    let reaches = skeleton.max_offset();
    for (i, node) in skeleton.nodes().iter().enumerate() {
        if node.kind == NodeKind::LimitReached && node.offset < reaches - BREAK_EPS {
            return Err(RoofError::UnevenLimits {
                node: NodeId(i as u32),
                stopped_at: node.offset,
                reaches,
            });
        }
    }
    Ok(())
}

/// Builds a roof's vertices, minting new ones where a profile's break cuts an
/// arc.
struct Build<'a> {
    skel: &'a Skeleton,
    profile: Profile,
    verts: Vec<RoofVertex>,
    /// Break-line corners, keyed by the **unordered** pair of skeleton nodes
    /// whose arc they sit on.
    ///
    /// Two panels share every arc, and both cut it at the same place. Minting a
    /// corner per panel would leave two vertices at one point, and a roof with a
    /// seam down every hip that a renderer would show as a crack. Keyed this way
    /// they get the same one.
    cuts: BTreeMap<(u32, u32), RoofVertexId>,
}

impl Build<'_> {
    /// Which side of the break a node's offset falls: -1 below, +1 above, 0 on.
    fn side(&self, n: NodeId, at: f32) -> i8 {
        let o = self.skel.node(n).offset;
        if o < at - BREAK_EPS {
            -1
        } else if o > at + BREAK_EPS {
            1
        } else {
            0
        }
    }

    /// The corner where the break crosses the arc between two nodes.
    ///
    /// Offset is an affine function of position, so the crossing is a simple
    /// interpolation — and the break's level set across a whole panel is a
    /// straight line, which is why cutting there leaves both halves flat.
    fn cut_between(&mut self, a: NodeId, b: NodeId, at: f32) -> Result<RoofVertexId, RoofError> {
        // Canonical order, so the two panels sharing this arc compute the very
        // same point rather than two that differ in the last bit.
        let (lo, hi) = if a.0 < b.0 { (a, b) } else { (b, a) };
        if let Some(&v) = self.cuts.get(&(lo.0, hi.0)) {
            return Ok(v);
        }

        let (nlo, nhi) = (self.skel.node(lo), self.skel.node(hi));
        let span = nhi.offset - nlo.offset;
        // Guarded, though a caller only asks when the ends straddle the break,
        // which needs them at least 2 * BREAK_EPS apart.
        let t = if span.abs() < f32::EPSILON {
            0.0
        } else {
            (at - nlo.offset) / span
        };
        let x = nlo.exact[0] + (nhi.exact[0] - nlo.exact[0]) * t;
        let y = nlo.exact[1] + (nhi.exact[1] - nlo.exact[1]) * t;
        let height = self.profile.height_at(at);
        let z = lattice_height(height).ok_or(RoofError::HeightOverflow { node: lo, height })?;

        let p = Point::from_vec2_rounded(Vec2::new(x, y));
        let id = RoofVertexId(self.verts.len() as u32);
        self.verts.push(RoofVertex {
            position: Point3::new(p.x, p.y, z),
            exact: [x, y, height],
            node: None,
        });
        self.cuts.insert((lo.0, hi.0), id);
        Ok(id)
    }

    /// Turns one wall's face into its panels: one if the break misses it, or the
    /// several the break cuts it into.
    ///
    /// The break is a level set of `offset`, and `offset` is affine across a
    /// single face (it is the distance to the wall's own line), so the break is a
    /// straight line **within this face**. Splitting a face is therefore
    /// splitting a simple polygon by a line — see [`Build::split_at_break`] for
    /// the general case, which cuts a face into as many lower and upper pieces as
    /// the line leaves it in.
    fn push_slopes(
        &mut self,
        wall: EdgeId,
        face: &[NodeId],
        out: &mut Vec<Panel>,
    ) -> Result<(), RoofError> {
        let Some(at) = self.profile.break_offset() else {
            out.push(Panel {
                kind: PanelKind::Slope { wall, band: 0 },
                corners: face.iter().map(|&n| RoofVertexId(n.0)).collect(),
            });
            return Ok(());
        };

        let sides: Vec<i8> = face.iter().map(|&n| self.side(n, at)).collect();

        // Wholly on one side: the break does not reach this panel, or clears it
        // entirely. Either way it stays in one piece, corners untouched. Corners
        // sitting *on* the break (side 0) do not cut it — they are the kerb the
        // panel runs up to, not a crossing of it.
        if sides.iter().all(|&s| s <= 0) || sides.iter().all(|&s| s >= 0) {
            let band = if sides.iter().any(|&s| s > 0) { 1 } else { 0 };
            out.push(Panel {
                kind: PanelKind::Slope { wall, band },
                corners: face.iter().map(|&n| RoofVertexId(n.0)).collect(),
            });
            return Ok(());
        }

        self.split_at_break(wall, face, &sides, at, out)
    }

    /// Splits a face the break genuinely crosses into its lower and upper pieces.
    ///
    /// A face can end up in more than two pieces: a wall whose wavefront was
    /// split apart by two reflex vertices has a face that dips below the break
    /// and back more than once, so the break cuts it into several lower pieces
    /// and several upper ones. This handles the general even-crossing partition,
    /// of which the ordinary "one lower, one upper" split is the two-crossing
    /// case. See [`Build::trace_band`] for the tracing itself.
    fn split_at_break(
        &mut self,
        wall: EdgeId,
        face: &[NodeId],
        sides: &[i8],
        at: f32,
        out: &mut Vec<Panel>,
    ) -> Result<(), RoofError> {
        // Order break-line corners along the break, which is parallel to the
        // wall. A face's first two nodes are the wall's own endpoints, so their
        // difference is the wall's direction — projecting onto it sorts corners
        // as they lie along the break.
        let a0 = self.skel.node(face[0]).exact;
        let a1 = self.skel.node(face[1]).exact;
        let (dx, dy) = (a1[0] - a0[0], a1[1] - a0[1]);
        let along = |v: RoofVertexId, verts: &[RoofVertex]| {
            let e = verts[v.0 as usize].exact;
            (e[0] - a0[0]) * dx + (e[1] - a0[1]) * dy
        };

        // The augmented boundary: the face's own corners, with a new one spliced
        // in wherever an edge crosses the break strictly (a `-1`/`+1` pair). That
        // makes every crossing of the break pass through a corner *on* it, so
        // each piece is bounded by real corners rather than mid-edge points.
        let mut bound: Vec<Bv> = Vec::with_capacity(face.len() + 4);
        let n = face.len();
        for i in 0..n {
            let v = RoofVertexId(face[i].0);
            let key = along(v, &self.verts);
            bound.push(Bv {
                id: v,
                side: sides[i],
                key,
            });

            let j = (i + 1) % n;
            if sides[i] * sides[j] < 0 {
                let c = self.cut_between(face[i], face[j], at)?;
                let key = along(c, &self.verts);
                bound.push(Bv {
                    id: c,
                    side: 0,
                    key,
                });
            }
        }

        let lower = self.trace_band(&bound, true, wall)?;
        let upper = self.trace_band(&bound, false, wall)?;

        // A sliver with fewer than three corners encloses nothing; dropping it
        // is what keeps a break laid exactly on a node from emitting a
        // degenerate panel. Lower pieces first, so a wall reads eaves-up.
        for (band, pieces) in [(0u8, lower), (1u8, upper)] {
            for corners in pieces {
                if corners.len() >= 3 {
                    out.push(Panel {
                        kind: PanelKind::Slope { wall, band },
                        corners,
                    });
                }
            }
        }
        Ok(())
    }

    /// Traces the closed pieces of one band of a cut face.
    ///
    /// The band's boundary is made of two kinds of edge: stretches of the face's
    /// own outline that lie on this side of the break, and chords *along* the
    /// break that close those stretches back up. The outline stretches are known
    /// outright; the chords are the puzzle, and they are recovered from where the
    /// outline meets the break.
    ///
    /// The break is a straight line, so the band meets it along a set of disjoint
    /// intervals. Each interval's two ends are an outline stretch ending on the
    /// break and another starting off it again — so, sorted along the break, the
    /// loose ends pair up consecutively, and joining each pair is the chord.
    fn trace_band(
        &self,
        bound: &[Bv],
        lower: bool,
        wall: EdgeId,
    ) -> Result<Vec<Vec<RoofVertexId>>, RoofError> {
        let m = bound.len();
        let in_band = |s: i8| if lower { s <= 0 } else { s >= 0 };

        // `next[i]` is the corner following corner `i` around this band. Seed it
        // from the outline: an edge is this band's boundary when both its ends
        // are on this side (and it is not a stretch lying *along* the break,
        // which is a chord's job, not an edge's).
        let mut next = alloc::vec![usize::MAX; m];
        let mut has_prev = alloc::vec![false; m];
        for i in 0..m {
            let j = (i + 1) % m;
            let (si, sj) = (bound[i].side, bound[j].side);
            if si == 0 && sj == 0 {
                continue;
            }
            if in_band(si) && in_band(sj) {
                next[i] = j;
                has_prev[j] = true;
            }
        }

        // The loose ends on the break: corners on it that the outline leaves by
        // (need an outgoing chord) or arrives at (need an incoming one).
        let mut needs_out: Vec<usize> = Vec::new();
        let mut needs_in: Vec<usize> = Vec::new();
        for i in 0..m {
            if bound[i].side != 0 {
                continue;
            }
            match (has_prev[i], next[i] != usize::MAX) {
                (true, false) => needs_out.push(i),
                (false, true) => needs_in.push(i),
                // Fully joined already (a corner the break merely touches), or
                // not part of this band at all.
                _ => {}
            }
        }

        // Sorted along the break, a loose end that needs a chord out and one that
        // needs a chord in alternate, and each adjacent pair bounds one interval
        // the band covers. Join them.
        let mut ends: Vec<usize> = needs_out.iter().chain(&needs_in).copied().collect();
        ends.sort_by(|&a, &b| bound[a].key.total_cmp(&bound[b].key));
        if ends.len() % 2 != 0 {
            return Err(RoofError::BreakSplitsPanel {
                wall,
                crossings: ends.len(),
            });
        }
        for pair in ends.chunks_exact(2) {
            let (a, b) = (pair[0], pair[1]);
            // One end must be an exit, the other an entry; anything else is
            // geometry this cannot make coherent pieces of.
            let out_in = match (next[a] == usize::MAX, next[b] == usize::MAX) {
                (true, false) => Some((a, b)),
                (false, true) => Some((b, a)),
                _ => None,
            };
            let Some((from, to)) = out_in else {
                return Err(RoofError::BreakSplitsPanel {
                    wall,
                    crossings: ends.len(),
                });
            };
            next[from] = to;
        }

        // Every loose end is joined now, so following `next` walks closed loops.
        let mut pieces = Vec::new();
        let mut seen = alloc::vec![false; m];
        for start in 0..m {
            if seen[start] || next[start] == usize::MAX {
                continue;
            }
            let mut corners = Vec::new();
            let mut cur = start;
            while !seen[cur] {
                seen[cur] = true;
                corners.push(bound[cur].id);
                cur = next[cur];
                if cur == usize::MAX {
                    return Err(RoofError::BreakSplitsPanel {
                        wall,
                        crossings: ends.len(),
                    });
                }
            }
            pieces.push(corners);
        }
        Ok(pieces)
    }
}

/// One corner of a face's augmented boundary while a break is being traced
/// through it: a face corner, or one the break spliced onto an edge.
struct Bv {
    /// The roof vertex this corner is.
    id: RoofVertexId,
    /// Which side of the break it falls: -1 below, +1 above, 0 on.
    side: i8,
    /// Its position projected along the break, for ordering the on-break corners.
    key: f32,
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
