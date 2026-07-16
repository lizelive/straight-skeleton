//! Integration tests for [`Roof`].
//!
//! The load-bearing test here is planarity: it holds only if the skeleton's
//! faces, node positions, and offsets are all correct *together*, so a single
//! misplaced node buckles its panel and trips it. It used to live in the `roof`
//! example; it belongs here now that the geometry does.

mod common;

use common::*;
use straight_skeleton::{
    skeleton, skeleton_constrained, EdgeId, Point, Point3, Polygon, Roof, RoofError,
};

const PITCH: f32 = 0.6;

/// Every panel must be flat, or it is not a roof you could build.
///
/// The panel's plane is pinned by its own wall: height must be
/// `pitch * distance-from-that-wall's-line` at every corner. Checking it that
/// way, rather than fitting a plane to the points, means a panel raised over
/// the *wrong* wall fails too.
///
/// This runs on `exact`, not `position`: rounding z to the lattice tilts each
/// panel by up to half a unit, which is documented on `RoofVertex::exact`.
fn assert_panels_are_planar(plan: &Polygon, roof: &Roof, pitch: f32) {
    for panel in roof.panels() {
        let mut slope: Option<f64> = None;
        for &c in &panel.corners {
            let v = roof.vertex(c);
            let run =
                signed_dist_to_edge_line(plan, panel.wall, [v.exact[0] as f64, v.exact[1] as f64]);
            let rise = v.exact[2] as f64;

            if run < 1e-6 {
                // At the eave, where the panel meets its own wall.
                assert!(
                    rise.abs() < 1e-3,
                    "panel over wall {} sits {rise} high on its own eave",
                    panel.wall,
                );
                continue;
            }
            let s = rise / run;
            match slope {
                None => slope = Some(s),
                Some(want) => assert!(
                    (s - want).abs() < 1e-3,
                    "panel over wall {} is not planar: slope {s} at one corner, {want} at another",
                    panel.wall,
                ),
            }
        }
        // Whatever slope the panel settled on, it must be the pitch asked for.
        if let Some(s) = slope {
            assert!(
                (s - pitch as f64).abs() < 1e-3,
                "panel over wall {} rises at {s}, not the requested pitch {pitch}",
                panel.wall,
            );
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
        assert_panels_are_planar(&plan, &roof, PITCH);
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
        assert_panels_are_planar(&plan, &roof, pitch);
        assert_eq!(roof.pitch(), pitch);
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

/// Roof vertices are indexed by `NodeId`, so they stand directly over their
/// skeleton nodes. This is what keeps provenance alive into 3D.
#[test]
fn vertices_line_up_with_skeleton_nodes() {
    let plan = Polygon::from_outer(&plus_shape()).unwrap();
    let skel = skeleton(&plan).unwrap();
    let roof = Roof::new(&skel, PITCH).unwrap();

    assert_eq!(roof.verts().len(), skel.node_count());
    for n in skel.node_ids() {
        let v = roof.vertex(n);
        assert_eq!(v.node, n);
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
        assert_eq!(panel.wall, EdgeId(i as u16));
        assert_eq!(roof.panel(EdgeId(i as u16)), panel);
        // A panel's outline starts along its own eave.
        let outline = roof.panel_outline(EdgeId(i as u16));
        assert_eq!(outline.len(), panel.corners.len());
        assert_eq!(outline[0].z, 0, "a panel starts on its eave");
        assert_eq!(outline[1].z, 0);
    }
}

/// The panels tile the plan, so their footprints must sum to its area.
#[test]
fn panel_footprints_tile_the_plan() {
    for (name, plan) in plans() {
        let roof = Roof::new(&skeleton(&plan).unwrap(), PITCH).unwrap();

        let total: f64 = roof
            .panels()
            .iter()
            .map(|panel| {
                let pts: Vec<[f64; 2]> = panel
                    .corners
                    .iter()
                    .map(|&n| {
                        let e = roof.vertex(n).exact;
                        [e[0] as f64, e[1] as f64]
                    })
                    .collect();
                let mut a = 0.0;
                for i in 0..pts.len() {
                    let p = pts[i];
                    let q = pts[(i + 1) % pts.len()];
                    a += p[0] * q[1] - q[0] * p[1];
                }
                (a / 2.0).abs()
            })
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

/// A constrained skeleton is truncated into disconnected stubs, so its faces
/// are not closed and there is no panel to raise. Saying so beats returning
/// something broken.
#[test]
fn refuses_to_roof_a_constrained_skeleton() {
    let plan = Polygon::from_outer(&rect(80, 80)).unwrap();
    let truncated = skeleton_constrained(&plan, &[5.0; 4]).unwrap();

    assert!(matches!(
        Roof::new(&truncated, 1.0),
        Err(RoofError::UnwalkableFace { .. })
    ));

    // Unlimited limits leave the skeleton intact, so that roofs fine.
    let intact = skeleton_constrained(&plan, &[f32::INFINITY; 4]).unwrap();
    assert_eq!(Roof::new(&intact, 1.0).unwrap().ridge_height(), 40);
}

#[test]
fn error_messages_are_useful() {
    let plan = Polygon::from_outer(&rect(80, 80)).unwrap();
    let skel = skeleton(&plan).unwrap();
    let msg = Roof::new(&skel, -1.0).unwrap_err().to_string();
    assert!(msg.contains("-1"), "unhelpful: {msg}");
}
