//! Shared checking utilities for the integration tests.
//!
//! The point of this module is that the assertions here are derived from the
//! *definition* of a straight skeleton, independently of how the crate
//! computes one. `check_invariants` is the workhorse: any skeleton, of any
//! polygon, must satisfy it.

#![allow(dead_code)]

use straight_skeleton::{Point, Polygon, Skeleton};

/// Exact-ish distance from a point to a segment, in `f64`.
pub fn dist_point_segment(p: [f64; 2], a: Point, b: Point) -> f64 {
    let (ax, ay) = (a.x as f64, a.y as f64);
    let (bx, by) = (b.x as f64, b.y as f64);
    let (dx, dy) = (bx - ax, by - ay);
    let len2 = dx * dx + dy * dy;
    let t = if len2 == 0.0 {
        0.0
    } else {
        (((p[0] - ax) * dx + (p[1] - ay) * dy) / len2).clamp(0.0, 1.0)
    };
    let (cx, cy) = (ax + t * dx, ay + t * dy);
    ((p[0] - cx).powi(2) + (p[1] - cy).powi(2)).sqrt()
}

/// The distance from a point to the polygon's boundary — the minimum over
/// every edge.
pub fn dist_to_boundary(poly: &Polygon, p: [f64; 2]) -> f64 {
    poly.edge_ids()
        .map(|e| {
            let (a, b) = poly.edge(e);
            dist_point_segment(p, a, b)
        })
        .fold(f64::INFINITY, f64::min)
}

/// Signed distance from a point to an edge's *supporting line*, positive on the
/// polygon's interior side.
///
/// This, not distance-to-segment, is the quantity a straight skeleton is built
/// from. See [`check_sources_are_equidistant`] for why that distinction is the
/// whole ballgame.
pub fn signed_dist_to_edge_line(poly: &Polygon, e: straight_skeleton::EdgeId, p: [f64; 2]) -> f64 {
    let (a, b) = poly.edge(e);
    let (ax, ay) = (a.x as f64, a.y as f64);
    let (bx, by) = (b.x as f64, b.y as f64);
    let (dx, dy) = (bx - ax, by - ay);
    let len = (dx * dx + dy * dy).sqrt();
    // Interior is to the left of a->b, so the inward normal is perp(dir).
    let (nx, ny) = (-dy / len, dx / len);
    nx * (p[0] - ax) + ny * (p[1] - ay)
}

/// Checks every property a straight skeleton must have, from first principles.
///
/// This is deliberately independent of the wavefront implementation: it takes
/// the polygon and the skeleton and re-derives what must be true.
pub fn check_invariants(poly: &Polygon, skel: &Skeleton, tol: f64) {
    check_boundary_nodes(poly, skel);
    check_sources_are_equidistant(poly, skel, tol);
    check_offset_does_not_exceed_boundary_distance(poly, skel, tol);
    check_arcs_are_bisectors(poly, skel, tol);
    check_arc_orientation(skel);
    check_graph_is_connected_and_wellformed(skel);
    check_no_nans(skel);
}

/// There is exactly one boundary node per input vertex, at that vertex, at
/// offset 0.
pub fn check_boundary_nodes(poly: &Polygon, skel: &Skeleton) {
    let mut seen = vec![false; poly.vertex_count()];
    for node in skel.nodes() {
        if let Some(v) = node.input_vertex() {
            assert!(
                !seen[v.0 as usize],
                "input vertex {v} has more than one boundary node"
            );
            seen[v.0 as usize] = true;
            assert_eq!(
                node.position,
                poly.vertex(v),
                "boundary node for {v} is not at that vertex"
            );
            assert_eq!(
                node.offset, 0.0,
                "boundary node for {v} has non-zero offset"
            );
        }
    }
    for (i, s) in seen.iter().enumerate() {
        assert!(s, "input vertex {i} has no boundary node");
    }
}

/// Every edge a node names as a source is exactly `offset` from it, measured to
/// that edge's **supporting line**.
///
/// This is the strongest check in the suite: it ties the simulation's notion of
/// time back to real geometry, so a node placed at the wrong spot or stamped
/// with the wrong time cannot pass.
///
/// It must be the *line*, not the segment, because a straight skeleton is not
/// the medial axis. Both bisect their input, but a straight skeleton only ever
/// bisects edges' infinite supporting lines, which is what keeps every arc
/// straight. A medial axis bisects the nearest *features*, and so grows
/// parabolic arcs around reflex vertices. Where they differ, this crate is
/// computing the straight skeleton, by definition.
pub fn check_sources_are_equidistant(poly: &Polygon, skel: &Skeleton, tol: f64) {
    for (i, node) in skel.nodes().iter().enumerate() {
        let p = [node.exact[0] as f64, node.exact[1] as f64];
        assert!(
            node.sources.len() >= 2,
            "node {i} lists {} sources; every node is equidistant from at least 2 edges",
            node.sources.len()
        );
        for &e in &node.sources {
            let d = signed_dist_to_edge_line(poly, e, p);
            assert!(
                (d - node.offset as f64).abs() < tol,
                "node {i} at {p:?} claims edge {e} as a source at offset {}, \
                 but its supporting line is {d} away",
                node.offset
            );
        }
    }
}

/// A skeleton node is never *closer* to the boundary than the wavefront had
/// travelled when it was made.
///
/// Equality holds wherever the polygon is locally convex. Near a reflex corner
/// the node can be strictly further away, because the reflex vertex sweeps out
/// along its bisector faster than its edges advance — the plus-shape's centre
/// sits at offset 5 but a full 7.07 from the nearest reflex corner. So this is
/// an inequality, and asserting equality here is exactly the mistake of
/// confusing this with a medial axis.
pub fn check_offset_does_not_exceed_boundary_distance(poly: &Polygon, skel: &Skeleton, tol: f64) {
    for (i, node) in skel.nodes().iter().enumerate() {
        let p = [node.exact[0] as f64, node.exact[1] as f64];
        let d = dist_to_boundary(poly, p);
        assert!(
            node.offset as f64 <= d + tol,
            "node {i} at {p:?} claims offset {} but is only {d} from the boundary",
            node.offset
        );
    }
}

/// For a **convex** polygon the straight skeleton and the medial axis coincide,
/// so there every node's offset is exactly its distance to the boundary, and
/// its sources really are its nearest edges.
///
/// Only valid for convex input; see
/// [`check_offset_does_not_exceed_boundary_distance`].
pub fn check_convex_offsets_are_boundary_distances(poly: &Polygon, skel: &Skeleton, tol: f64) {
    assert!(
        poly.vertex_ids().all(|v| !poly.is_reflex(v)),
        "this check is only meaningful for convex polygons"
    );
    for (i, node) in skel.nodes().iter().enumerate() {
        let p = [node.exact[0] as f64, node.exact[1] as f64];
        let d = dist_to_boundary(poly, p);
        assert!(
            (d - node.offset as f64).abs() < tol,
            "convex polygon: node {i} at {p:?} claims offset {} but is {d} from the boundary",
            node.offset
        );
    }
}

/// An arc's midpoint is equidistant from the supporting lines of both of the
/// arc's source edges, and that shared distance is the midpoint's interpolated
/// offset.
///
/// Sampling the midpoint rather than an endpoint matters: endpoints are shared
/// with other arcs and would pass trivially, whereas the midpoint is only
/// equidistant if the arc genuinely bisects the pair along its whole length.
pub fn check_arcs_are_bisectors(poly: &Polygon, skel: &Skeleton, tol: f64) {
    for a in skel.arc_ids() {
        let arc = skel.arc(a);
        let n0 = skel.node(arc.nodes[0]);
        let n1 = skel.node(arc.nodes[1]);
        let mid = [
            (n0.exact[0] as f64 + n1.exact[0] as f64) * 0.5,
            (n0.exact[1] as f64 + n1.exact[1] as f64) * 0.5,
        ];

        let [e0, e1] = arc.sources;
        assert_ne!(e0, e1, "arc {a:?} bisects an edge with itself");

        let d0 = signed_dist_to_edge_line(poly, e0, mid);
        let d1 = signed_dist_to_edge_line(poly, e1, mid);
        assert!(
            (d0 - d1).abs() < tol,
            "arc {a:?} midpoint {mid:?} is {d0} from source {e0}'s line \
             but {d1} from source {e1}'s line"
        );

        // Offset is linear along an arc, so the midpoint's offset is the mean
        // of its endpoints'. That shared bisected distance must equal it.
        let want = (n0.offset as f64 + n1.offset as f64) * 0.5;
        assert!(
            (d0 - want).abs() < tol,
            "arc {a:?} midpoint bisects its sources at {d0}, but its \
             interpolated offset is {want}"
        );
    }
}

/// Arcs point away from the boundary: `nodes[0]` is the lower-offset end.
pub fn check_arc_orientation(skel: &Skeleton) {
    for a in skel.arc_ids() {
        let arc = skel.arc(a);
        let lo = skel.node(arc.nodes[0]).offset;
        let hi = skel.node(arc.nodes[1]).offset;
        assert!(lo <= hi + 1e-4, "arc {a:?} runs downhill: {lo} -> {hi}");
    }
}

/// The graph is sane: no self-loops, every boundary node has degree 1, and the
/// whole thing is connected.
pub fn check_graph_is_connected_and_wellformed(skel: &Skeleton) {
    if skel.node_count() == 0 {
        return;
    }
    for a in skel.arc_ids() {
        let arc = skel.arc(a);
        assert_ne!(arc.nodes[0], arc.nodes[1], "arc {a:?} is a self-loop");
    }

    for n in skel.node_ids() {
        let deg = skel.arcs_at(n).len();
        if skel.node(n).is_boundary() {
            assert_eq!(
                deg, 1,
                "boundary node {n:?} has degree {deg}; each input vertex emits exactly one arc"
            );
        } else {
            assert!(deg >= 1, "interior node {n:?} is isolated");
        }
    }

    // Flood fill from node 0.
    let mut seen = vec![false; skel.node_count()];
    let mut stack = vec![0usize];
    seen[0] = true;
    while let Some(i) = stack.pop() {
        for &a in skel.arcs_at(straight_skeleton::NodeId(i as u32)) {
            let arc = skel.arc(a);
            for &m in &arc.nodes {
                if !seen[m.0 as usize] {
                    seen[m.0 as usize] = true;
                    stack.push(m.0 as usize);
                }
            }
        }
    }
    let unreached: Vec<_> = seen
        .iter()
        .enumerate()
        .filter(|(_, &s)| !s)
        .map(|(i, _)| i)
        .collect();
    assert!(
        unreached.is_empty(),
        "skeleton is disconnected; nodes {unreached:?} are unreachable"
    );
}

/// No node position or offset is NaN or infinite.
pub fn check_no_nans(skel: &Skeleton) {
    for (i, n) in skel.nodes().iter().enumerate() {
        assert!(n.exact[0].is_finite(), "node {i} has non-finite x");
        assert!(n.exact[1].is_finite(), "node {i} has non-finite y");
        assert!(n.offset.is_finite(), "node {i} has non-finite offset");
        assert!(
            n.offset >= -1e-4,
            "node {i} has negative offset {}",
            n.offset
        );
    }
}

/// A rectangle, CCW from the origin.
pub fn rect(w: i16, h: i16) -> Vec<Point> {
    vec![
        Point::new(0, 0),
        Point::new(w, 0),
        Point::new(w, h),
        Point::new(0, h),
    ]
}

/// A regular-ish polygon inscribed in a circle, CCW.
pub fn ngon(n: usize, r: f64, cx: f64, cy: f64) -> Vec<Point> {
    (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64) / (n as f64);
            Point::new(
                (cx + r * a.cos()).round() as i16,
                (cy + r * a.sin()).round() as i16,
            )
        })
        .collect()
}

/// An L-shaped polygon with exactly one reflex vertex.
pub fn l_shape() -> Vec<Point> {
    vec![
        Point::new(0, 0),
        Point::new(40, 0),
        Point::new(40, 20),
        Point::new(20, 20),
        Point::new(20, 40),
        Point::new(0, 40),
    ]
}

/// A plus/cross shape with four reflex vertices.
pub fn plus_shape() -> Vec<Point> {
    vec![
        Point::new(10, 0),
        Point::new(20, 0),
        Point::new(20, 10),
        Point::new(30, 10),
        Point::new(30, 20),
        Point::new(20, 20),
        Point::new(20, 30),
        Point::new(10, 30),
        Point::new(10, 20),
        Point::new(0, 20),
        Point::new(0, 10),
        Point::new(10, 10),
    ]
}
