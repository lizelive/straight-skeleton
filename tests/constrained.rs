//! Integration tests for the per-edge distance-limited skeleton.

mod common;

use common::*;
use straight_skeleton::{
    skeleton, skeleton_constrained, EdgeId, NodeId, NodeKind, Point, Polygon, SkeletonError,
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

/// The counterpart to [`limits_apply_per_edge`], holding back a **short** wall
/// rather than a long one — which is a different situation, not the same one
/// rotated.
///
/// A 40x20 rectangle's ridge runs along the long axis at offset 10, and it is
/// the two *long* edges meeting that puts it there. So limiting a long edge
/// lowers the ridge, but limiting a short edge cannot: the long edges still meet
/// at offset 10 regardless. What it does instead is **lengthen** the ridge. The
/// short wall's two corners stop bisecting once it freezes and slide straight
/// along it, so they meet further out than they otherwise would have.
#[test]
fn limiting_a_short_edge_lengthens_the_ridge_rather_than_lowering_it() {
    let poly = Polygon::from_outer(&rect(40, 20)).unwrap();
    // Edge 1 is the right-hand short wall, x = 40.
    let limits = [f32::INFINITY, 3.0, f32::INFINITY, f32::INFINITY];
    let skel = skeleton_constrained(&poly, &limits).unwrap();
    check_constrained(&poly, &skel, &limits, TOL);

    // All derived from the shape, not read off the output. The wall stops with
    // its line at x = 37. Its corners have risen at 45 degrees to (37, 3) and
    // (37, 17) by then; from there each must stay on both its long edge's line
    // (y = t, still moving) and the frozen x = 37, so each slides vertically at
    // unit speed.
    for want in [Point::new(37, 3), Point::new(37, 17)] {
        assert!(
            skel.nodes()
                .iter()
                .any(|n| n.position == want && n.kind == NodeKind::LimitReached),
            "expected the wall to stop and its corner to kink at {want:?}"
        );
    }

    // Those two corners are 14 apart and closing at 2, so they meet at (37, 10)
    // at t = 10 -- which is also when the long edges meet at y = 10. The ridge
    // therefore runs (10, 10) to (37, 10), reaching 7 further right than the
    // unconstrained (10, 10)-(30, 10).
    for want in [Point::new(37, 10), Point::new(10, 10)] {
        assert!(
            skel.nodes().iter().any(|n| n.position == want),
            "expected a ridge end at {want:?}"
        );
    }

    // And the ridge is no lower, because the edges that set its height were
    // never limited. This is the whole contrast with `limits_apply_per_edge`.
    let plain = skeleton(&poly).unwrap();
    assert!(
        (skel.max_offset() - plain.max_offset()).abs() < 1e-4,
        "limiting a short edge must not lower the ridge: {} vs {}",
        skel.max_offset(),
        plain.max_offset()
    );
}

/// A uniform limit leaves the input polygon, offset inward by the limit. That
/// residual outline is half of what a constrained skeleton *is* — the arcs are
/// the stubs reaching in, and this is what they stop on.
#[test]
fn a_uniform_limit_leaves_the_polygon_offset_inward() {
    // The L-shape, every wall stopped at 20. Offsetting each wall inward by 20
    // and re-intersecting gives a smaller L with the same six corners: the
    // reflex corner at (100, 100) moves to (80, 80), *outward* along its
    // bisector, because reflex corners run backwards into the material.
    let poly = Polygon::from_outer(&[
        Point::new(0, 0),
        Point::new(200, 0),
        Point::new(200, 100),
        Point::new(100, 100),
        Point::new(100, 200),
        Point::new(0, 200),
    ])
    .unwrap();
    let skel = skeleton_constrained(&poly, &[20.0; 6]).unwrap();

    assert_eq!(skel.residual().len(), 1, "one loop, around the outer ring");
    let flat = &skel.residual()[0];
    assert_eq!(flat.len(), 6);

    let mut corners: Vec<Point> = flat.nodes.iter().map(|&n| skel.node(n).position).collect();
    corners.sort();
    assert_eq!(
        corners,
        vec![
            Point::new(20, 20),
            Point::new(20, 180),
            Point::new(80, 80),
            Point::new(80, 180),
            Point::new(180, 20),
            Point::new(180, 80),
        ]
    );

    // Every corner is a node where the wavefront stopped, at the limit.
    for &n in &flat.nodes {
        assert_eq!(skel.node(n).kind, NodeKind::LimitReached);
        assert!((skel.node(n).offset - 20.0).abs() < 1e-4);
    }

    check_residual_is_parallel_to_its_walls(&poly, &skel);
    check_constrained(&poly, &skel, &[20.0; 6], TOL);
}

/// The faces of a constrained skeleton close, by running along the residual
/// where no arc bounds them — and together with the residual they still tile the
/// polygon exactly.
///
/// This is the invariant that says the residual really is the missing boundary
/// rather than a separate object bolted alongside: the faces cover everything
/// the wavefront swept, the residual covers what it never reached, and there is
/// nothing left over.
#[test]
fn constrained_faces_and_the_residual_tile_the_polygon() {
    let l = vec![
        Point::new(0, 0),
        Point::new(200, 0),
        Point::new(200, 100),
        Point::new(100, 100),
        Point::new(100, 200),
        Point::new(0, 200),
    ];
    let cases: Vec<(Polygon, Vec<f32>)> = vec![
        (Polygon::from_outer(&square(100)).unwrap(), vec![20.0; 4]),
        (Polygon::from_outer(&rect(200, 100)).unwrap(), vec![20.0; 4]),
        (Polygon::from_outer(&l).unwrap(), vec![20.0; 6]),
        // A limit that binds on only part of the shape: the narrow arm's
        // wavefront collapses before 20, the wide part survives to it.
        (Polygon::from_outer(&l).unwrap(), vec![45.0; 6]),
        (
            Polygon::new(
                &rect(200, 150),
                &[vec![
                    Point::new(60, 50),
                    Point::new(140, 50),
                    Point::new(140, 100),
                    Point::new(60, 100),
                ]],
            )
            .unwrap(),
            vec![10.0; 8],
        ),
    ];

    for (poly, limits) in cases {
        let skel = skeleton_constrained(&poly, &limits).unwrap();
        let faces = skel
            .faces()
            .unwrap_or_else(|| panic!("a constrained face failed to close at limit {}", limits[0]));
        assert_eq!(faces.len(), poly.edge_count());

        // Faces are disjoint regions, so their areas add. The residual loops are
        // not: a loop around a surviving hole is a hole *in* the region the
        // wavefront never reached, so it must subtract. Summing the signed areas
        // does exactly that, which is the winding convention doing real work
        // rather than being a tidy-looking promise.
        let swept: f64 = faces.iter().map(|f| ring_area(&skel, f).abs()).sum();
        let unswept: f64 = skel
            .residual()
            .iter()
            .map(|l| residual_area(&skel, l))
            .sum();

        let want = poly.signed_area2() as f64 / 2.0;
        let total = swept + unswept;
        assert!(
            (total - want).abs() < 0.02 * want,
            "swept {swept} + unswept {unswept} = {total}, but the polygon \
             encloses {want} (limit {})",
            limits[0]
        );
    }
}

/// Twice the signed area of a ring of nodes, halved: the shoelace formula.
fn ring_area(skel: &straight_skeleton::Skeleton, ring: &[NodeId]) -> f64 {
    let mut a = 0.0;
    for i in 0..ring.len() {
        let p = skel.node(ring[i]).exact;
        let q = skel.node(ring[(i + 1) % ring.len()]).exact;
        a += p[0] as f64 * q[1] as f64 - q[0] as f64 * p[1] as f64;
    }
    a / 2.0
}

/// A plain skeleton has no residual at all: its wavefront always collapses to
/// nothing, which is what it means for it to have finished.
#[test]
fn an_unconstrained_skeleton_has_no_residual() {
    let cases: Vec<Vec<Point>> = vec![
        rect(10, 10),
        rect(40, 10),
        l_shape(),
        plus_shape(),
        ngon(7, 300.0, 0.0, 0.0),
    ];
    for pts in cases {
        let poly = Polygon::from_outer(&pts).unwrap();
        assert!(
            skeleton(&poly).unwrap().residual().is_empty(),
            "a plain skeleton's wavefront must shrink away entirely"
        );
        // Limits that never bind must reproduce that, since the wavefront still
        // collapses before it reaches them.
        let far = vec![1e6; poly.edge_count()];
        assert!(skeleton_constrained(&poly, &far)
            .unwrap()
            .residual()
            .is_empty());
    }
}

/// A hole that survives leaves its own residual loop, wound the other way.
#[test]
fn a_surviving_hole_leaves_its_own_residual_loop() {
    // A 200x150 rectangle with a 80x50 hole in the middle. At a limit of 10
    // neither the outer wavefront nor the hole's has met anything.
    let poly = Polygon::new(
        &[
            Point::new(0, 0),
            Point::new(200, 0),
            Point::new(200, 150),
            Point::new(0, 150),
        ],
        &[vec![
            Point::new(60, 50),
            Point::new(140, 50),
            Point::new(140, 100),
            Point::new(60, 100),
        ]],
    )
    .unwrap();
    let n = poly.edge_count();
    let skel = skeleton_constrained(&poly, &vec![10.0; n]).unwrap();

    assert_eq!(skel.residual().len(), 2, "outer outline, plus the hole's");
    check_residual_is_parallel_to_its_walls(&poly, &skel);

    // The hole's wavefront expands into the material, so its loop grows: the
    // 80x50 hole becomes 100x70. The outer one shrinks, 200x150 to 180x130.
    let areas: Vec<f64> = skel
        .residual()
        .iter()
        .map(|l| residual_area(&skel, l).abs())
        .collect();
    let mut sorted = areas.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert!(
        (sorted[0] - 100.0 * 70.0).abs() < 1.0,
        "hole should offset outward to 100x70, got area {}",
        sorted[0]
    );
    assert!(
        (sorted[1] - 180.0 * 130.0).abs() < 1.0,
        "outer should offset inward to 180x130, got area {}",
        sorted[1]
    );

    // Winding is inherited from the input, so the interior stays on the left of
    // every segment: outer counter-clockwise (positive), hole clockwise.
    let signed: Vec<f64> = skel
        .residual()
        .iter()
        .map(|l| residual_area(&skel, l))
        .collect();
    assert_eq!(
        signed.iter().filter(|a| **a > 0.0).count(),
        1,
        "exactly one loop winds counter-clockwise: the outer one"
    );
    assert_eq!(signed.iter().filter(|a| **a < 0.0).count(), 1);
}

/// Twice the signed area of a residual loop; positive when counter-clockwise.
fn residual_area(skel: &straight_skeleton::Skeleton, l: &straight_skeleton::ResidualLoop) -> f64 {
    let mut a = 0.0;
    for (from, to, _) in l.segments() {
        let p = skel.node(from).exact;
        let q = skel.node(to).exact;
        a += p[0] as f64 * q[1] as f64 - q[0] as f64 * p[1] as f64;
    }
    a / 2.0
}

/// Each residual segment is **parallel** to the wall it names, and `min(offset,
/// limit)` from that wall's supporting line.
///
/// This is the property that says the residual really is the input offset
/// inward, and it is also the reason these segments are not [`Arc`]s: an arc
/// bisects two edges, whereas this runs alongside exactly one.
fn check_residual_is_parallel_to_its_walls(poly: &Polygon, skel: &straight_skeleton::Skeleton) {
    for l in skel.residual() {
        assert_eq!(l.nodes.len(), l.edges.len(), "one source edge per segment");
        for (from, to, e) in l.segments() {
            let (a, b) = poly.edge(e);
            let wall = [(b.x - a.x) as f64, (b.y - a.y) as f64];
            let p = skel.node(from).exact;
            let q = skel.node(to).exact;
            let seg = [(q[0] - p[0]) as f64, (q[1] - p[1]) as f64];

            // Parallel, and pointing the same way as the wall.
            let cross = wall[0] * seg[1] - wall[1] * seg[0];
            let dot = wall[0] * seg[0] + wall[1] * seg[1];
            let len = (seg[0] * seg[0] + seg[1] * seg[1]).sqrt();
            assert!(
                cross.abs() < 1e-2 * len.max(1.0),
                "residual segment {p:?}->{q:?} is not parallel to its wall {e}"
            );
            assert!(dot > 0.0, "residual segment {p:?}->{q:?} runs against {e}");

            // Both ends sit the same distance out from that wall's line.
            for end in [from, to] {
                let x = skel.node(end).exact;
                let d = signed_dist_to_edge_line(poly, e, [x[0] as f64, x[1] as f64]);
                assert!(
                    (d - skel.node(end).offset as f64).abs() < TOL,
                    "residual corner is {d} from wall {e}, not its offset {}",
                    skel.node(end).offset
                );
            }
        }
    }
}

/// A limit of zero on a **single** wall: that wall never moves at all, while
/// its neighbours run free.
///
/// This is the sharpest test of the "an edge that stops is just an edge at
/// speed zero" model, because the stopped edge never has a moving phase to be
/// wrong about. Its two corners must slide straight **along** it rather than
/// away from it, so both arcs leaving that wall stay pinned to `y = 0`.
#[test]
fn a_single_wall_with_a_zero_limit_never_moves() {
    let poly = Polygon::from_outer(&rect(40, 20)).unwrap();
    // Edge 0 is the bottom wall, y = 0. It is frozen from the start.
    let limits = [0.0, f32::INFINITY, f32::INFINITY, f32::INFINITY];
    let skel = skeleton_constrained(&poly, &limits).unwrap();

    check_constrained(&poly, &skel, &limits, TOL);

    // The other three walls are unconstrained, so the wavefront still runs.
    assert!(
        skel.max_offset() > 0.0,
        "only one wall was limited; the rest must keep going"
    );

    // Every node touching the frozen wall must lie *on* it.
    let touching: Vec<_> = skel
        .nodes()
        .iter()
        .filter(|n| n.sources.contains(&EdgeId(0)))
        .collect();
    assert!(
        touching.len() >= 2,
        "the frozen wall still bounds part of the skeleton"
    );
    for n in &touching {
        assert!(
            n.exact[1].abs() < 1e-3,
            "node at {:?} claims the frozen wall but has left it",
            n.exact
        );
    }

    // Its two corners slide along it, inward from (0,0) and (40,0), rather
    // than lifting off along a 45-degree bisector.
    let mut on_wall: Vec<i16> = touching
        .iter()
        .filter(|n| !n.is_boundary())
        .map(|n| n.position.x)
        .collect();
    on_wall.sort();
    assert!(
        on_wall.iter().all(|&x| x > 0 && x < 40),
        "corners should slide inward along the wall, got x = {on_wall:?}"
    );
}

/// Zero on *every* wall is the degenerate limit: nothing moves, so nothing is
/// traced.
#[test]
fn a_zero_limit_on_every_wall_traces_nothing() {
    let poly = Polygon::from_outer(&l_shape()).unwrap();
    let limits = vec![0.0; poly.edge_count()];
    let skel = skeleton_constrained(&poly, &limits).unwrap();

    assert_eq!(skel.max_offset(), 0.0);
    assert_eq!(skel.arc_count(), 0);
    assert_eq!(skel.node_count(), poly.vertex_count());
    check_boundary_nodes(&poly, &skel);
}

/// A zero limit on a wall of a polygon that has reflex corners, so the frozen
/// wall has to interact with split events rather than just its own neighbours.
#[test]
fn a_zero_limit_wall_on_a_reflex_polygon() {
    let poly = Polygon::from_outer(&l_shape()).unwrap();
    let mut limits = vec![f32::INFINITY; poly.edge_count()];
    limits[0] = 0.0;
    let skel = skeleton_constrained(&poly, &limits).unwrap();

    check_constrained(&poly, &skel, &limits, TOL);
    for n in skel
        .nodes()
        .iter()
        .filter(|n| n.sources.contains(&EdgeId(0)))
    {
        // Edge 0 of the L runs along y = 0.
        assert!(
            n.exact[1].abs() < 1e-3,
            "node at {:?} claims the frozen wall but has left it",
            n.exact
        );
    }
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
