//! Integration tests for the per-edge distance-limited skeleton.

mod common;

use common::*;
use straight_skeleton::{
    skeleton, skeleton_constrained, EdgeId, NodeKind, Point, Polygon, SkeletonError,
};

const TOL: f64 = 1e-2;

/// The checks that still apply once limits truncate the skeleton.
///
/// Some of the unconstrained invariants do not survive here, and each for a
/// real reason rather than convenience. Once limits bind, `offset` stops being
/// a distance at all: it is the wavefront's **time**, and an edge that stopped
/// early is nearer than that time suggests.
///
/// - **Connectivity** goes: a constrained skeleton is deliberately
///   disconnected. When every edge stops, the arcs are left as disjoint stubs
///   reaching in from the boundary. That is the point of the transform.
/// - **Source distance** changes: a node is `min(offset, limit)` from a source
///   edge's line, not `offset`. This is a *stronger* check than the
///   unconstrained one, since it pins down where each edge stopped too.
/// - **The boundary-distance bound** goes: a node at time 17 can sit 3 from an
///   edge that stopped at 3, so the wavefront legitimately ends up closer to
///   the boundary than its time.
fn check_constrained(poly: &Polygon, skel: &straight_skeleton::Skeleton, limits: &[f32], tol: f64) {
    check_boundary_nodes(poly, skel);
    check_arc_orientation(skel);
    check_no_nans(skel);
    check_sources_respect_limits(poly, skel, limits, tol);
}

/// Every source edge's line is `min(offset, limit)` from the node.
fn check_sources_respect_limits(
    poly: &Polygon,
    skel: &straight_skeleton::Skeleton,
    limits: &[f32],
    tol: f64,
) {
    for (i, node) in skel.nodes().iter().enumerate() {
        let p = [node.exact[0] as f64, node.exact[1] as f64];
        assert!(node.sources.len() >= 2, "node {i} has too few sources");
        for &e in &node.sources {
            let want = (node.offset as f64).min(limits[e.0 as usize] as f64);
            let got = signed_dist_to_edge_line(poly, e, p);
            assert!(
                (got - want).abs() < tol,
                "node {i} at {p:?} (offset {}) claims edge {e}, whose limit is {}: \
                 its line should be {want} away but is {got}",
                node.offset,
                limits[e.0 as usize],
            );
        }
    }
}

fn square(size: i16) -> Vec<Point> {
    vec![
        Point::new(0, 0),
        Point::new(size, 0),
        Point::new(size, size),
        Point::new(0, size),
    ]
}

/// The headline promise: an unlimited constrained skeleton *is* the plain one.
/// Both run the same weighted wavefront, so this must hold exactly.
#[test]
fn infinite_limits_reproduce_the_plain_skeleton() {
    for pts in [square(20), l_shape(), plus_shape(), rect(30, 10)] {
        let poly = Polygon::from_outer(&pts).unwrap();
        let plain = skeleton(&poly).unwrap();
        let limits = vec![f32::INFINITY; poly.edge_count()];
        let constrained = skeleton_constrained(&poly, &limits).unwrap();

        assert_eq!(plain.node_count(), constrained.node_count());
        assert_eq!(plain.arc_count(), constrained.arc_count());
        assert_eq!(plain.nodes(), constrained.nodes());
        assert_eq!(plain.arcs(), constrained.arcs());
    }
}

/// A limit beyond the ridge cannot bind, so it changes nothing.
#[test]
fn a_limit_larger_than_the_ridge_changes_nothing() {
    let poly = Polygon::from_outer(&square(20)).unwrap();
    let plain = skeleton(&poly).unwrap();

    // The ridge is at 10; 50 is far out of reach.
    let constrained = skeleton_constrained(&poly, &[50.0; 4]).unwrap();
    assert_eq!(plain.nodes(), constrained.nodes());
    assert_eq!(plain.arcs(), constrained.arcs());
    assert!(constrained
        .nodes()
        .iter()
        .all(|n| n.kind != NodeKind::LimitReached));
}

#[test]
fn a_uniform_limit_truncates_the_skeleton() {
    let poly = Polygon::from_outer(&square(20)).unwrap();
    let skel = skeleton_constrained(&poly, &[3.0; 4]).unwrap();

    assert!(
        skel.max_offset() <= 3.0 + 1e-4,
        "max offset {} exceeds the limit",
        skel.max_offset()
    );
    // Four corners run inward and stop.
    assert_eq!(skel.arc_count(), 4);
    assert_eq!(
        skel.nodes()
            .iter()
            .filter(|n| n.kind == NodeKind::LimitReached)
            .count(),
        4
    );
    check_constrained(&poly, &skel, &[3.0; 4], TOL);
}

/// The corner of a square runs along the 45-degree diagonal, so stopping the
/// wavefront at distance `d` leaves the node at exactly `(d, d)`.
#[test]
fn a_truncated_corner_lands_where_the_geometry_says() {
    let poly = Polygon::from_outer(&square(40)).unwrap();
    let skel = skeleton_constrained(&poly, &[6.0; 4]).unwrap();

    let mut stops: Vec<Point> = skel
        .nodes()
        .iter()
        .filter(|n| n.kind == NodeKind::LimitReached)
        .map(|n| n.position)
        .collect();
    stops.sort();
    assert_eq!(
        stops,
        vec![
            Point::new(6, 6),
            Point::new(6, 34),
            Point::new(34, 6),
            Point::new(34, 34)
        ]
    );
    for n in skel.nodes().iter().filter(|n| !n.is_boundary()) {
        assert!((n.offset - 6.0).abs() < 1e-3);
    }
}

#[test]
fn a_limit_of_zero_leaves_only_the_boundary() {
    let poly = Polygon::from_outer(&square(20)).unwrap();
    let skel = skeleton_constrained(&poly, &[0.0; 4]).unwrap();

    assert_eq!(skel.max_offset(), 0.0);
    // Nothing moved, so nothing was traced.
    assert_eq!(skel.arc_count(), 0);
    assert_eq!(skel.node_count(), 4);
    check_boundary_nodes(&poly, &skel);
}

/// Limits are per edge, so holding one edge back must bend the arcs that touch
/// it and leave the rest alone.
#[test]
fn limits_apply_per_edge() {
    let poly = Polygon::from_outer(&rect(40, 20)).unwrap();
    // Only the bottom edge is limited.
    let limits = [3.0, f32::INFINITY, f32::INFINITY, f32::INFINITY];
    let skel = skeleton_constrained(&poly, &limits).unwrap();

    // The unlimited edges still meet well above the limited edge's stop.
    assert!(
        skel.max_offset() > 3.0,
        "unlimited edges should keep going past 3, got {}",
        skel.max_offset()
    );
    // The bottom edge stopping bends the two arcs that ride on it.
    assert!(skel
        .nodes()
        .iter()
        .any(|n| n.kind == NodeKind::LimitReached));
    check_constrained(&poly, &skel, &limits, TOL);
}

/// Once an edge stops, its neighbours keep moving and slide *along* it, rather
/// than over it. So nothing on that edge's side gets further from it than the
/// limit allowed.
#[test]
fn a_stopped_edge_holds_its_offset() {
    let poly = Polygon::from_outer(&rect(60, 40)).unwrap();
    let limit = 4.0;
    let skel =
        skeleton_constrained(&poly, &[limit, f32::INFINITY, f32::INFINITY, f32::INFINITY]).unwrap();

    // Edge 0 is the bottom, y = 0. Nothing claiming it as a source may sit
    // further from it than its limit.
    for (i, n) in skel.nodes().iter().enumerate() {
        if n.sources.contains(&EdgeId(0)) && !n.is_boundary() {
            assert!(
                n.exact[1] <= limit + 1e-3,
                "node {i} claims the stopped bottom edge but is at y = {}",
                n.exact[1]
            );
        }
    }
}

#[test]
fn different_limits_on_different_edges() {
    let poly = Polygon::from_outer(&square(40)).unwrap();
    let limits = [2.0, 5.0, 9.0, 14.0];
    let skel = skeleton_constrained(&poly, &limits).unwrap();
    check_constrained(&poly, &skel, &limits, TOL);
    assert!(skel.max_offset() <= 14.0 + 1e-3);
}

#[test]
fn constrained_skeleton_of_a_polygon_with_a_hole() {
    let poly = Polygon::new(
        &square(40),
        &[vec![
            Point::new(15, 15),
            Point::new(25, 15),
            Point::new(25, 25),
            Point::new(15, 25),
        ]],
    )
    .unwrap();
    let limits = vec![3.0; poly.edge_count()];
    let skel = skeleton_constrained(&poly, &limits).unwrap();

    assert!(skel.max_offset() <= 3.0 + 1e-3);
    check_constrained(&poly, &skel, &limits, TOL);
}

#[test]
fn constrained_l_shape_is_wellformed() {
    let poly = Polygon::from_outer(&l_shape()).unwrap();
    let limits = vec![5.0; poly.edge_count()];
    let skel = skeleton_constrained(&poly, &limits).unwrap();
    assert!(skel.max_offset() <= 5.0 + 1e-3);
    check_constrained(&poly, &skel, &limits, TOL);
}

/// Increasing a limit can only ever reveal more skeleton, never less.
#[test]
fn raising_the_limit_is_monotone() {
    let poly = Polygon::from_outer(&l_shape()).unwrap();
    let n = poly.edge_count();

    let mut last = 0.0f32;
    for d in [1.0f32, 2.0, 4.0, 8.0, 10.0] {
        let skel = skeleton_constrained(&poly, &vec![d; n]).unwrap();
        assert!(
            skel.max_offset() >= last - 1e-3,
            "limit {d} reached {} but limit before reached {last}",
            skel.max_offset()
        );
        assert!(skel.max_offset() <= d + 1e-3);
        last = skel.max_offset();
    }
}

// --- Error handling ---------------------------------------------------------

#[test]
fn rejects_a_wrong_number_of_limits() {
    let poly = Polygon::from_outer(&square(10)).unwrap();
    assert_eq!(
        skeleton_constrained(&poly, &[1.0, 2.0]).unwrap_err(),
        SkeletonError::LimitCountMismatch {
            got: 2,
            expected: 4
        }
    );
    assert!(skeleton_constrained(&poly, &[]).is_err());
}

#[test]
fn rejects_negative_and_nan_limits() {
    let poly = Polygon::from_outer(&square(10)).unwrap();

    let e = skeleton_constrained(&poly, &[1.0, -1.0, 1.0, 1.0]).unwrap_err();
    assert!(matches!(
        e,
        SkeletonError::InvalidLimit {
            edge: EdgeId(1),
            ..
        }
    ));

    let e = skeleton_constrained(&poly, &[1.0, 1.0, f32::NAN, 1.0]).unwrap_err();
    assert!(matches!(
        e,
        SkeletonError::InvalidLimit {
            edge: EdgeId(2),
            ..
        }
    ));
}

/// Two collinear edges given different limits would tear the wavefront: one
/// line stops while the other, parallel to it, keeps going, and the vertex
/// between them has nowhere to be. The crate refuses rather than invent
/// something.
#[test]
fn rejects_incompatible_limits_on_collinear_edges() {
    // The bottom edge is split in two by a collinear mid-edge vertex.
    let poly = Polygon::from_outer(&[
        Point::new(0, 0),
        Point::new(10, 0), // collinear
        Point::new(20, 0),
        Point::new(20, 20),
        Point::new(0, 20),
    ])
    .unwrap();

    // Edges 0 and 1 are the two collinear halves. Give them different limits.
    let e = skeleton_constrained(
        &poly,
        &[2.0, 7.0, f32::INFINITY, f32::INFINITY, f32::INFINITY],
    )
    .unwrap_err();
    assert!(
        matches!(e, SkeletonError::IncompatibleCollinearLimits { .. }),
        "got {e:?}"
    );

    // The same limit on both is fine.
    assert!(skeleton_constrained(
        &poly,
        &[2.0, 2.0, f32::INFINITY, f32::INFINITY, f32::INFINITY]
    )
    .is_ok());
}

#[test]
fn error_messages_name_the_problem() {
    let poly = Polygon::from_outer(&square(10)).unwrap();
    let e = skeleton_constrained(&poly, &[1.0, 2.0]).unwrap_err();
    let msg = e.to_string();
    assert!(msg.contains('2') && msg.contains('4'), "unhelpful: {msg}");
}
