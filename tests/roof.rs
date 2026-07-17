//! Integration tests for [`Roof`].
//!
//! The load-bearing check is `assert_panels_match_profile`: every corner is
//! re-derived from the distance to the wall its panel claims to rise from. It
//! holds only if the skeleton's faces, node positions and offsets are all
//! correct *together*, so a single misplaced node trips it — and being derived
//! rather than fitted, it catches a panel raised over the wrong wall too.

mod common;

use common::*;
use straight_skeleton::{
    skeleton, skeleton_constrained, EdgeId, PanelKind, Point, Point3, Polygon, Profile, Roof,
    RoofError,
};

const PITCH: f32 = 0.6;

/// Every corner of every panel stands exactly where the profile says, measured
/// from that panel's **own wall**.
///
/// This is the load-bearing check, and it is stronger than testing planarity.
/// Fitting a plane to a panel's points would pass for any flat panel, however
/// wrongly placed. Re-deriving each corner's height from the distance to the
/// wall the panel claims to rise from means a panel raised over the *wrong*
/// wall fails, a node placed at the wrong offset fails, and a mansard whose
/// break landed at the wrong height fails.
///
/// Planarity then follows rather than being asserted: distance-to-a-line is
/// affine in position, and each band of a profile is affine in distance, so a
/// panel that satisfies this is a plane by construction.
///
/// Runs on `exact`, not `position`: rounding z to the lattice tilts each panel
/// by up to half a unit, which is documented on `RoofVertex::exact`.
fn assert_panels_match_profile(plan: &Polygon, roof: &Roof) {
    for panel in roof.panels() {
        match panel.kind {
            PanelKind::Slope { wall, .. } => {
                let runs: Vec<f64> = panel
                    .corners
                    .iter()
                    .map(|&c| {
                        let v = roof.vertex(c);
                        signed_dist_to_edge_line(plan, wall, [v.exact[0] as f64, v.exact[1] as f64])
                    })
                    .collect();

                // A wall frozen at a limit of zero sweeps nothing, so every
                // corner of its panel sits on the wall's own line and the panel
                // stands vertical: a gable. There is no height-to-derive-from-run
                // there, since the run is zero all the way up. Its geometry is
                // pinned by `a_zero_limit_wall_becomes_a_gable` instead, and by
                // the fact that its corners are shared with the sloping panels
                // either side, which *are* checked here.
                if runs.iter().all(|r| r.abs() < 1e-2) {
                    continue;
                }

                for (&c, &run) in panel.corners.iter().zip(&runs) {
                    let want = roof.profile().height_at(run as f32) as f64;
                    let got = roof.vertex(c).exact[2] as f64;
                    assert!(
                        (got - want).abs() < 1e-2,
                        "panel over wall {wall} has a corner {run} from that wall, \
                         so the profile puts it at {want}, but it stands at {got}"
                    );
                }
            }
            // The flat is level, by definition. Nothing to derive from a wall,
            // so the claim is simply that every corner is at one height.
            PanelKind::Flat => {
                let z0 = roof.vertex(panel.corners[0]).exact[2] as f64;
                for &c in &panel.corners {
                    let z = roof.vertex(c).exact[2] as f64;
                    assert!(
                        (z - z0).abs() < 1e-2,
                        "the flat is not level: {z} at one corner, {z0} at another"
                    );
                }
            }
        }
    }
}

fn plans() -> Vec<(&'static str, Polygon)> {
    vec![
        ("square", Polygon::from_outer(&rect(80, 80)).unwrap()),
        ("rectangle", Polygon::from_outer(&rect(120, 80)).unwrap()),
        ("l-shape", Polygon::from_outer(&l_shape()).unwrap()),
        ("plus", Polygon::from_outer(&plus_shape()).unwrap()),
        (
            "hexagon",
            Polygon::from_outer(&ngon(6, 100.0, 0.0, 0.0)).unwrap(),
        ),
        (
            "courtyard",
            Polygon::new(
                &rect(200, 160),
                &[vec![
                    Point::new(70, 55),
                    Point::new(130, 55),
                    Point::new(130, 105),
                    Point::new(70, 105),
                ]],
            )
            .unwrap(),
        ),
    ]
}

#[test]
fn every_panel_is_planar() {
    for (name, plan) in plans() {
        let roof = Roof::new(&skeleton(&plan).unwrap(), PITCH).unwrap();
        assert_panels_match_profile(&plan, &roof);
        assert_eq!(
            roof.panels().len(),
            plan.edge_count(),
            "{name}: one panel per wall"
        );
    }
}

#[test]
fn panels_stay_planar_across_pitches() {
    let plan = Polygon::from_outer(&l_shape()).unwrap();
    let skel = skeleton(&plan).unwrap();
    for pitch in [0.0f32, 0.1, 0.25, 0.5, 1.0, 2.0] {
        let roof = Roof::new(&skel, pitch).unwrap();
        assert_panels_match_profile(&plan, &roof);
        assert_eq!(roof.profile(), Profile::Hip { pitch });
    }
}

/// A square plan gives a pyramid whose apex is half the width in from any wall.
#[test]
fn a_square_plan_gives_a_pyramid() {
    let plan = Polygon::from_outer(&rect(80, 80)).unwrap();
    let roof = Roof::new(&skeleton(&plan).unwrap(), 1.0).unwrap();

    // Four eaves at zero, one apex.
    let eaves: Vec<_> = roof.verts().iter().filter(|v| v.position.z == 0).collect();
    assert_eq!(eaves.len(), 4);
    assert_eq!(roof.verts().len(), 5);

    // The apex is 40 in from each wall, so at pitch 1.0 it stands 40 high.
    assert_eq!(roof.ridge_height(), 40);
    let apex = roof.verts().iter().find(|v| v.position.z > 0).unwrap();
    assert_eq!(apex.position, Point3::new(40, 40, 40));
}

/// A 2:1 plan gives a ridge, not a point — the case that catches naive
/// implementations.
#[test]
fn a_rectangular_plan_gives_a_ridge() {
    let plan = Polygon::from_outer(&rect(120, 80)).unwrap();
    let roof = Roof::new(&skeleton(&plan).unwrap(), 0.5).unwrap();

    // The ridge runs along the middle of the long axis, 40 in from each long
    // wall, so at pitch 0.5 it stands 20 high.
    assert_eq!(roof.ridge_height(), 20);

    let mut ridge: Vec<Point3> = roof
        .verts()
        .iter()
        .filter(|v| v.position.z == 20)
        .map(|v| v.position)
        .collect();
    ridge.sort();
    assert_eq!(
        ridge,
        vec![Point3::new(40, 40, 20), Point3::new(80, 40, 20)]
    );
}

#[test]
fn pitch_scales_the_height_linearly() {
    let plan = Polygon::from_outer(&rect(80, 80)).unwrap();
    let skel = skeleton(&plan).unwrap();

    assert_eq!(Roof::new(&skel, 0.0).unwrap().ridge_height(), 0);
    assert_eq!(Roof::new(&skel, 0.25).unwrap().ridge_height(), 10);
    assert_eq!(Roof::new(&skel, 0.5).unwrap().ridge_height(), 20);
    assert_eq!(Roof::new(&skel, 1.0).unwrap().ridge_height(), 40);
    assert_eq!(Roof::new(&skel, 2.0).unwrap().ridge_height(), 80);
}

#[test]
fn a_flat_roof_is_flat() {
    let plan = Polygon::from_outer(&l_shape()).unwrap();
    let roof = Roof::new(&skeleton(&plan).unwrap(), 0.0).unwrap();
    assert!(roof.verts().iter().all(|v| v.position.z == 0));
    assert_eq!(roof.ridge_height(), 0);
}

/// A hip roof's vertices line up with the skeleton's nodes one for one, so they
/// stand directly over them. This is what keeps provenance alive into 3D.
#[test]
fn vertices_line_up_with_skeleton_nodes() {
    let plan = Polygon::from_outer(&plus_shape()).unwrap();
    let skel = skeleton(&plan).unwrap();
    let roof = Roof::new(&skel, PITCH).unwrap();

    assert_eq!(roof.verts().len(), skel.node_count());
    for n in skel.node_ids() {
        let v = roof.vertex_at(n);
        assert_eq!(v.node, Some(n));
        assert_eq!(v.position.x, skel.node(n).position.x);
        assert_eq!(v.position.y, skel.node(n).position.y);
        // Height is the node's offset scaled by the pitch.
        assert!((v.exact[2] - skel.node(n).offset * PITCH).abs() < 1e-4);
    }
}

#[test]
fn every_panel_knows_its_wall() {
    let plan = Polygon::from_outer(&l_shape()).unwrap();
    let roof = Roof::new(&skeleton(&plan).unwrap(), PITCH).unwrap();

    for (i, panel) in roof.panels().iter().enumerate() {
        let wall = EdgeId(i as u16);
        assert_eq!(panel.wall(), Some(wall));
        assert_eq!(panel.kind, PanelKind::Slope { wall, band: 0 });
        assert!(!panel.is_flat());
        // A hip roof gives each wall exactly one panel.
        assert_eq!(roof.panels_of(wall).count(), 1);
        // A panel's outline starts along its own eave.
        let outline = roof.outline(panel);
        assert_eq!(outline.len(), panel.corners.len());
        assert_eq!(outline[0].z, 0, "a panel starts on its eave");
        assert_eq!(outline[1].z, 0);
    }
    assert_eq!(roof.flat().count(), 0, "a plain hip roof has no flat");
}

/// The panels tile the plan, so their footprints must sum to its area.
#[test]
fn panel_footprints_tile_the_plan() {
    for (name, plan) in plans() {
        let roof = Roof::new(&skeleton(&plan).unwrap(), PITCH).unwrap();

        let total: f64 = roof
            .panels()
            .iter()
            .map(|panel| footprint_area(&roof, panel).abs())
            .sum();

        let want = plan.signed_area2() as f64 / 2.0;
        assert!(
            (total - want).abs() < 0.05 * want,
            "{name}: panels cover {total} but the plan is {want}"
        );
    }
}

// --- Error handling ---------------------------------------------------------

#[test]
fn rejects_an_invalid_pitch() {
    let plan = Polygon::from_outer(&rect(80, 80)).unwrap();
    let skel = skeleton(&plan).unwrap();

    for bad in [-1.0f32, -0.001, f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
        assert!(
            matches!(Roof::new(&skel, bad), Err(RoofError::InvalidPitch { .. })),
            "pitch {bad} should be rejected"
        );
    }
}

/// Overflow is refused rather than clamped: a silently flattened ridge is a
/// wrong roof, not an approximate one.
#[test]
fn refuses_to_overflow_rather_than_flattening_the_ridge() {
    // The widest plan the coordinate cap allows.
    let plan = Polygon::from_outer(&rect(16000, 16000)).unwrap();
    let skel = skeleton(&plan).unwrap();

    // The apex is 8000 in, so a pitch of 5 wants 40000 — past i16.
    assert!(matches!(
        Roof::new(&skel, 5.0),
        Err(RoofError::HeightOverflow { .. })
    ));

    // A pitch that does fit is fine. Heights may exceed the *coordinate* cap,
    // which only constrains the plan: z is not fed back in as input.
    let roof = Roof::new(&skel, 4.0).unwrap();
    assert_eq!(roof.ridge_height(), 32000);
}

/// **Uneven** limits have no roof at all, and that is a fact about geometry
/// rather than a gap in the implementation.
///
/// Height is a function of `offset`, which only means "distance from the wall"
/// while nothing has stopped early. An edge that halted at 3 stays 3 from its
/// face however long the clock runs, so its panel would want to end lower than
/// its neighbour's and the surface between them would have to tear.
#[test]
fn refuses_a_skeleton_whose_limits_are_uneven() {
    let plan = Polygon::from_outer(&rect(80, 80)).unwrap();

    // One wall held back while the rest run on.
    let uneven =
        skeleton_constrained(&plan, &[5.0, f32::INFINITY, f32::INFINITY, f32::INFINITY]).unwrap();
    assert!(matches!(
        Roof::new(&uneven, 1.0),
        Err(RoofError::UnevenLimits { .. })
    ));

    // Two different finite limits are just as torn.
    let mixed = skeleton_constrained(&plan, &[5.0, 20.0, 20.0, 20.0]).unwrap();
    assert!(matches!(
        Roof::new(&mixed, 1.0),
        Err(RoofError::UnevenLimits { .. })
    ));

    // Limits that never bind leave the skeleton intact, so that roofs fine.
    let intact = skeleton_constrained(&plan, &[f32::INFINITY; 4]).unwrap();
    assert_eq!(Roof::new(&intact, 1.0).unwrap().ridge_height(), 40);
}

#[test]
fn rejects_an_invalid_mansard_break() {
    let plan = Polygon::from_outer(&rect(80, 80)).unwrap();
    let skel = skeleton(&plan).unwrap();

    for bad in [-1.0f32, f32::NAN, f32::INFINITY] {
        assert!(
            matches!(
                Roof::mansard(&skel, 1.0, bad, 0.5),
                Err(RoofError::InvalidBreak { .. })
            ),
            "break {bad} should be rejected"
        );
    }
    // Both pitches are still policed.
    assert!(matches!(
        Roof::mansard(&skel, -1.0, 10.0, 0.5),
        Err(RoofError::InvalidPitch { .. })
    ));
    assert!(matches!(
        Roof::mansard(&skel, 1.0, 10.0, f32::NAN),
        Err(RoofError::InvalidPitch { .. })
    ));
}

#[test]
fn error_messages_are_useful() {
    let plan = Polygon::from_outer(&rect(80, 80)).unwrap();
    let skel = skeleton(&plan).unwrap();
    let msg = Roof::new(&skel, -1.0).unwrap_err().to_string();
    assert!(msg.contains("-1"), "unhelpful: {msg}");

    let uneven =
        skeleton_constrained(&plan, &[5.0, f32::INFINITY, f32::INFINITY, f32::INFINITY]).unwrap();
    let msg = Roof::new(&uneven, 1.0).unwrap_err().to_string();
    assert!(msg.contains("same"), "unhelpful: {msg}");

    let msg = Roof::mansard(&skel, 1.0, -3.0, 0.5)
        .unwrap_err()
        .to_string();
    assert!(msg.contains("-3"), "unhelpful: {msg}");
}

// --- Truncated roofs: a uniform limit gives a flat top ----------------------

/// A uniform limit is a hip roof with its apex cut off, and the flat it stops
/// at is the skeleton's residual raised to the limit's height.
#[test]
fn a_uniform_limit_gives_a_flat_top() {
    let plan = Polygon::from_outer(&rect(80, 80)).unwrap();
    // The apex is 40 in; stopping at 10 cuts it well short.
    let skel = skeleton_constrained(&plan, &[10.0; 4]).unwrap();
    let roof = Roof::new(&skel, 1.0).unwrap();

    assert_panels_match_profile(&plan, &roof);

    // Four slopes and one flat.
    assert_eq!(roof.panels().len(), 5);
    assert_eq!(roof.flat().count(), 1);
    for w in 0..4u16 {
        assert_eq!(roof.panels_of(EdgeId(w)).count(), 1);
    }

    // At pitch 1.0 the flat stands at the limit, and is the top of the roof.
    assert_eq!(roof.ridge_height(), 10);

    // It is the 60x60 square left in the middle of the 80x80 plan.
    let flat = roof.flat().next().unwrap();
    let mut corners: Vec<Point3> = roof.outline(flat);
    corners.sort();
    assert_eq!(
        corners,
        vec![
            Point3::new(10, 10, 10),
            Point3::new(10, 70, 10),
            Point3::new(70, 10, 10),
            Point3::new(70, 70, 10),
        ]
    );
}

#[test]
fn truncated_roofs_stay_planar_over_every_plan() {
    for (name, plan) in plans() {
        // Half way to the ridge, so the limit is sure to bite on every plan
        // whatever its size. A fixed limit would silently clear the ridge of the
        // smaller ones and test nothing.
        let ridge = skeleton(&plan).unwrap().max_offset();
        let limit = ridge / 2.0;

        let skel = skeleton_constrained(&plan, &vec![limit; plan.edge_count()]).unwrap();
        let roof = Roof::new(&skel, 0.75).unwrap();

        assert_panels_match_profile(&plan, &roof);
        assert!(
            roof.flat().count() >= 1,
            "{name}: a limit of {limit} is half the ridge at {ridge}, so it must bite"
        );
        // Nothing rises above the flat.
        let want = (limit * 0.75).round() as i16;
        assert_eq!(roof.ridge_height(), want, "{name}: {limit} * 0.75");
    }
}

/// A courtyard's hole survives the limit, so the flat has a hole in it: two
/// `Flat` panels, wound opposite ways.
#[test]
fn a_flat_with_a_hole_in_it_is_two_panels() {
    let plan = Polygon::new(
        &rect(200, 160),
        &[vec![
            Point::new(70, 55),
            Point::new(130, 55),
            Point::new(130, 105),
            Point::new(70, 105),
        ]],
    )
    .unwrap();
    let skel = skeleton_constrained(&plan, &vec![8.0; plan.edge_count()]).unwrap();
    let roof = Roof::new(&skel, 1.0).unwrap();

    assert_panels_match_profile(&plan, &roof);
    assert_eq!(roof.flat().count(), 2, "an outline and a hole in it");

    // Opposite windings say which is which.
    let signed: Vec<f64> = roof.flat().map(|p| footprint_area(&roof, p)).collect();
    assert_eq!(signed.iter().filter(|a| **a > 0.0).count(), 1);
    assert_eq!(signed.iter().filter(|a| **a < 0.0).count(), 1);
}

// --- Mansards ---------------------------------------------------------------

/// The headline: a mansard is steep, then shallow, and the break is a level line
/// all the way round.
#[test]
fn a_mansard_breaks_where_it_is_told_to() {
    // 120 x 80: the ridge is 40 in from the long walls.
    let plan = Polygon::from_outer(&rect(120, 80)).unwrap();
    let skel = skeleton(&plan).unwrap();
    let roof = Roof::mansard(&skel, 2.0, 10.0, 0.25).unwrap();

    assert_panels_match_profile(&plan, &roof);

    // Each wall now carries two panels: the steep skirt and the shallow slope.
    assert_eq!(roof.panels().len(), 8);
    for w in 0..4u16 {
        let bands: Vec<u8> = roof
            .panels_of(EdgeId(w))
            .map(|p| match p.kind {
                PanelKind::Slope { band, .. } => band,
                PanelKind::Flat => unreachable!(),
            })
            .collect();
        assert_eq!(
            bands,
            vec![0, 1],
            "wall {w}: a steep band then a shallow one"
        );
    }

    // Derived by hand: 10 of run at 2.0 is 20 at the kerb, then the remaining
    // 30 of run at 0.25 adds 7.5 -> 27.5, which rounds to 28.
    assert_eq!(roof.ridge_height(), 28);

    // The kerb is level: every corner the break introduced stands at 20.
    let kerb: Vec<&straight_skeleton::RoofVertex> =
        roof.verts().iter().filter(|v| v.node.is_none()).collect();
    assert!(!kerb.is_empty(), "the break must introduce corners");
    for v in kerb {
        assert!(
            (v.exact[2] - 20.0).abs() < 1e-3,
            "the kerb should be level at 20, found {}",
            v.exact[2]
        );
    }
}

/// Equal pitches make a mansard indistinguishable from a hip roof, which is the
/// sharpest test that the break itself introduces no geometry of its own.
#[test]
fn a_mansard_with_equal_pitches_is_a_hip_roof() {
    for (name, plan) in plans() {
        let skel = skeleton(&plan).unwrap();
        let hip = Roof::new(&skel, 0.5).unwrap();
        let mansard = Roof::mansard(&skel, 0.5, 6.0, 0.5).unwrap();

        assert_eq!(hip.ridge_height(), mansard.ridge_height(), "{name}");
        assert_panels_match_profile(&plan, &mansard);

        // The footprints still tile the plan, whatever the break did to them.
        let a: f64 = hip.panels().iter().map(|p| footprint_area(&hip, p)).sum();
        let b: f64 = mansard
            .panels()
            .iter()
            .map(|p| footprint_area(&mansard, p))
            .sum();
        assert!((a - b).abs() < 0.05 * a.abs(), "{name}: {a} vs {b}");
    }
}

/// A break beyond the ridge never bites, so only the lower pitch is ever used.
#[test]
fn a_break_out_of_reach_leaves_a_plain_hip_roof() {
    let plan = Polygon::from_outer(&rect(80, 80)).unwrap();
    let skel = skeleton(&plan).unwrap();

    // The apex is 40 in; break at 100.
    let roof = Roof::mansard(&skel, 1.0, 100.0, 0.1).unwrap();
    assert_eq!(roof.panels().len(), 4, "no panel is cut");
    assert_eq!(roof.ridge_height(), 40, "the upper pitch never applies");
    assert!(roof.verts().iter().all(|v| v.node.is_some()));
    assert_panels_match_profile(&plan, &roof);
}

/// A mansard over a uniform limit: steep, then shallow, then flat.
#[test]
fn a_constrained_mansard_is_steep_then_shallow_then_flat() {
    let plan = Polygon::from_outer(&rect(120, 80)).unwrap();
    // Stop at 20, which is past the break at 10 but short of the ridge at 40.
    let skel = skeleton_constrained(&plan, &[20.0; 4]).unwrap();
    let roof = Roof::mansard(&skel, 2.0, 10.0, 0.25).unwrap();

    assert_panels_match_profile(&plan, &roof);

    // Two bands per wall, plus the flat they stop at.
    assert_eq!(roof.panels().len(), 4 * 2 + 1);
    assert_eq!(roof.flat().count(), 1);

    // 10 at 2.0 = 20 at the kerb; then 10 more of run at 0.25 = 22.5 -> 23.
    assert_eq!(roof.ridge_height(), 23);
    for &c in &roof.flat().next().unwrap().corners {
        assert!((roof.vertex(c).exact[2] - 22.5).abs() < 1e-3);
    }
}

/// A limit *below* the break stops the roof before the shallow band ever
/// starts: the result is a truncated steep roof with no upper band at all.
#[test]
fn a_limit_below_the_break_never_reaches_the_shallow_band() {
    let plan = Polygon::from_outer(&rect(120, 80)).unwrap();
    let skel = skeleton_constrained(&plan, &[6.0; 4]).unwrap();
    let roof = Roof::mansard(&skel, 2.0, 10.0, 0.25).unwrap();

    assert_panels_match_profile(&plan, &roof);
    assert_eq!(
        roof.panels().len(),
        4 + 1,
        "one band per wall, plus the flat"
    );
    assert!(roof.panels().iter().all(|p| p.kind
        != PanelKind::Slope {
            wall: EdgeId(0),
            band: 1
        }));
    assert_eq!(roof.ridge_height(), 12, "6 * 2.0, all in the steep band");
}

// --- Gables: a zero limit freezes a wall ------------------------------------

/// A limit of **zero** is how you gable a wall, and it composes with everything
/// else here. It is also the case that makes `UnevenLimits` subtler than "all
/// the limits must match": zero mixes with unlimited freely, because a wall that
/// never moves has no sloping panel to be inconsistent about.
///
/// Geometry derived by hand. A 240 x 90 hall, both 90-long ends frozen: the long
/// walls still close in on each other and meet along `y = 45` at offset 45. With
/// nothing eating in from the ends, that ridge runs the whole 240 — out to
/// `(0, 45)` and `(240, 45)`, which sit *on* the frozen walls. Each frozen wall
/// is then left with a vertical triangle standing on it: eaves at both corners,
/// apex at the ridge. That is a gable.
#[test]
fn a_zero_limit_wall_becomes_a_gable() {
    let plan = Polygon::from_outer(&rect(240, 90)).unwrap();
    // e1 is the right end (x = 240), e3 the left (x = 0).
    let limits = [f32::INFINITY, 0.0, f32::INFINITY, 0.0];
    let skel = skeleton_constrained(&plan, &limits).unwrap();
    let roof = Roof::new(&skel, 1.0).unwrap();

    assert_panels_match_profile(&plan, &roof);
    assert_eq!(roof.ridge_height(), 45, "the ridge is half the 90 span");

    // The ridge reaches the frozen walls rather than hipping away from them.
    let mut ridge: Vec<Point3> = roof
        .verts()
        .iter()
        .filter(|v| v.position.z == 45)
        .map(|v| v.position)
        .collect();
    ridge.sort();
    assert_eq!(
        ridge,
        vec![Point3::new(0, 45, 45), Point3::new(240, 45, 45)],
        "a gabled ridge runs right out to both end walls"
    );

    // Each frozen wall's panel is a vertical triangle standing on it.
    for (wall, x) in [(EdgeId(1), 240i16), (EdgeId(3), 0)] {
        let gable = roof.panels_of(wall).next().unwrap();
        let mut corners = roof.outline(gable);
        corners.sort();
        assert_eq!(
            corners,
            vec![
                Point3::new(x, 0, 0),
                Point3::new(x, 45, 45),
                Point3::new(x, 90, 0)
            ],
            "wall {wall} should be a gable triangle standing at x = {x}"
        );
        // Vertical means no plan area at all.
        assert!(
            footprint_area(&roof, gable).abs() < 1e-3,
            "a gable sweeps nothing, so its footprint is degenerate"
        );
    }

    // The unfrozen walls still get real sloping panels, running the full length.
    for wall in [EdgeId(0), EdgeId(2)] {
        let slope = roof.panels_of(wall).next().unwrap();
        assert!((footprint_area(&roof, slope).abs() - 240.0 * 45.0).abs() < 1.0);
    }
}

/// A **gambrel** is a mansard with gable ends — the barn roof — and it needs
/// both constraints at once: the mansard's break and the gables' zero limits.
///
/// The break cuts the gable too, which is exactly right: a gambrel's end wall is
/// the classic kinked pentagon, not a triangle.
#[test]
fn a_gambrel_is_a_mansard_with_gable_ends() {
    let plan = Polygon::from_outer(&rect(240, 90)).unwrap();
    let limits = [f32::INFINITY, 0.0, f32::INFINITY, 0.0];
    let skel = skeleton_constrained(&plan, &limits).unwrap();
    // Steep to 10, then shallow.
    let roof = Roof::mansard(&skel, 2.0, 10.0, 0.25).unwrap();

    assert_panels_match_profile(&plan, &roof);

    // 10 of run at 2.0 is 20 at the kerb; the remaining 35 at 0.25 adds 8.75.
    assert_eq!(roof.ridge_height(), 29);

    // The end wall comes out as a pentagon, split into a band either side of the
    // kerb: eaves at 0, kink at 20, apex at 28.75.
    let gable: Vec<&straight_skeleton::Panel> = roof.panels_of(EdgeId(1)).collect();
    assert_eq!(gable.len(), 2, "the break cuts the gable in two");
    let mut zs: Vec<i16> = gable
        .iter()
        .flat_map(|p| roof.outline(p))
        .map(|p| p.z)
        .collect();
    zs.sort();
    zs.dedup();
    assert_eq!(zs, vec![0, 20, 29], "eaves, kerb, apex");

    // Still vertical, still standing on x = 240.
    for p in gable {
        assert!(roof.outline(p).iter().all(|c| c.x == 240));
        assert!(footprint_area(&roof, p).abs() < 1e-3);
    }
}

/// Twice the signed footprint area of a panel, halved: the shoelace formula on
/// its x/y, ignoring height.
fn footprint_area(roof: &Roof, panel: &straight_skeleton::Panel) -> f64 {
    let pts: Vec<[f64; 2]> = panel
        .corners
        .iter()
        .map(|&c| {
            let e = roof.vertex(c).exact;
            [e[0] as f64, e[1] as f64]
        })
        .collect();
    let mut a = 0.0;
    for i in 0..pts.len() {
        let p = pts[i];
        let q = pts[(i + 1) % pts.len()];
        a += p[0] * q[1] - q[0] * p[1];
    }
    a / 2.0
}
