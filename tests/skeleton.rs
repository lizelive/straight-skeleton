//! Integration tests for the unconstrained straight skeleton.
//!
//! Expected geometry here is derived by hand from the shape, not read off the
//! implementation, so these tests can actually catch the implementation being
//! wrong.

mod common;

use common::*;
use straight_skeleton::{skeleton, NodeKind, Point, Polygon};

/// A tolerance in coordinate units. Node positions are narrowed to `f32` on
/// output, which costs about 1e-3 of absolute precision at coordinates in the
/// tens of thousands, so checks are held to a little looser than that.
const TOL: f64 = 1e-2;

#[test]
fn square_skeleton_is_an_x() {
    let poly = Polygon::from_outer(&rect(10, 10)).unwrap();
    let skel = skeleton(&poly).unwrap();

    // Four corners, one centre.
    assert_eq!(skel.node_count(), 5);
    assert_eq!(skel.arc_count(), 4);

    let centre: Vec<_> = skel.nodes().iter().filter(|n| !n.is_boundary()).collect();
    assert_eq!(centre.len(), 1);
    assert_eq!(centre[0].position, Point::new(5, 5));
    assert!((centre[0].offset - 5.0).abs() < TOL as f32);
    // The centre is equidistant from all four sides.
    assert_eq!(centre[0].sources.len(), 4);

    check_invariants(&poly, &skel, TOL);
}

/// The canonical case the naive pairwise merge gets wrong: a rectangle's
/// skeleton is a ridge, not a point.
#[test]
fn rectangle_skeleton_is_a_ridge() {
    let poly = Polygon::from_outer(&rect(20, 10)).unwrap();
    let skel = skeleton(&poly).unwrap();

    // Four corners plus the two ridge ends.
    assert_eq!(skel.node_count(), 6);
    // Four corner arcs plus the ridge itself.
    assert_eq!(skel.arc_count(), 5);

    let mut ridge: Vec<_> = skel
        .nodes()
        .iter()
        .filter(|n| !n.is_boundary())
        .map(|n| n.position)
        .collect();
    ridge.sort();
    assert_eq!(ridge, vec![Point::new(5, 5), Point::new(15, 5)]);

    // The ridge sits at half the rectangle's short side.
    for n in skel.nodes().iter().filter(|n| !n.is_boundary()) {
        assert!((n.offset - 5.0).abs() < TOL as f32);
    }
    assert!((skel.max_offset() - 5.0).abs() < TOL as f32);

    check_invariants(&poly, &skel, TOL);
}

#[test]
fn tall_rectangle_ridge_is_vertical() {
    let poly = Polygon::from_outer(&rect(10, 20)).unwrap();
    let skel = skeleton(&poly).unwrap();

    let mut ridge: Vec<_> = skel
        .nodes()
        .iter()
        .filter(|n| !n.is_boundary())
        .map(|n| n.position)
        .collect();
    ridge.sort();
    assert_eq!(ridge, vec![Point::new(5, 5), Point::new(5, 15)]);

    check_invariants(&poly, &skel, TOL);
}

#[test]
fn right_triangle_meets_at_the_incenter() {
    // The 9-12-15 right triangle. Its inradius is (9 + 12 - 15) / 2 = 3,
    // so the incenter is 3 from each leg: (3, 3).
    let poly =
        Polygon::from_outer(&[Point::new(0, 0), Point::new(12, 0), Point::new(0, 9)]).unwrap();
    let skel = skeleton(&poly).unwrap();

    assert_eq!(skel.node_count(), 4);
    assert_eq!(skel.arc_count(), 3);

    let inc: Vec<_> = skel.nodes().iter().filter(|n| !n.is_boundary()).collect();
    assert_eq!(inc.len(), 1);
    assert_eq!(inc[0].position, Point::new(3, 3));
    assert!((inc[0].offset - 3.0).abs() < TOL as f32);
    assert_eq!(inc[0].sources.len(), 3);

    check_invariants(&poly, &skel, TOL);
}

#[test]
fn any_convex_polygon_has_no_split_events() {
    for n in [3usize, 4, 5, 6, 7, 8, 12, 20] {
        let poly = Polygon::from_outer(&ngon(n, 100.0, 0.0, 0.0)).unwrap();
        let skel = skeleton(&poly).unwrap();

        assert!(
            skel.nodes()
                .iter()
                .all(|nd| nd.kind != NodeKind::SplitEvent),
            "{n}-gon is convex, so it cannot have a split event"
        );
        check_invariants(&poly, &skel, TOL);
    }
}

/// On convex input the straight skeleton *is* the medial axis, so the stronger
/// claim holds there: a node's offset is exactly its distance to the boundary.
/// (It does not hold in general — see the plus-shape.)
#[test]
fn on_convex_input_offset_is_exactly_the_boundary_distance() {
    let cases: Vec<Vec<Point>> = vec![
        rect(10, 10),
        rect(40, 10),
        ngon(5, 300.0, 0.0, 0.0),
        ngon(7, 900.0, 0.0, 0.0),
        vec![Point::new(0, 0), Point::new(12, 0), Point::new(0, 9)],
    ];
    for pts in cases {
        let poly = Polygon::from_outer(&pts).unwrap();
        let skel = skeleton(&poly).unwrap();
        check_convex_offsets_are_boundary_distances(&poly, &skel, 0.05);
    }
}

/// Every input edge owns exactly one face, and the faces tile the polygon: their
/// areas must sum to the polygon's own.
#[test]
fn faces_tile_the_polygon() {
    let cases: Vec<Polygon> = vec![
        Polygon::from_outer(&rect(10, 10)).unwrap(),
        Polygon::from_outer(&rect(40, 10)).unwrap(),
        Polygon::from_outer(&l_shape()).unwrap(),
        Polygon::from_outer(&plus_shape()).unwrap(),
        Polygon::from_outer(&ngon(6, 200.0, 0.0, 0.0)).unwrap(),
        Polygon::new(
            &rect(30, 30),
            &[vec![
                Point::new(10, 10),
                Point::new(20, 10),
                Point::new(20, 20),
                Point::new(10, 20),
            ]],
        )
        .unwrap(),
    ];

    for poly in cases {
        let skel = skeleton(&poly).unwrap();
        let faces = skel.faces().expect("every edge must have a walkable face");
        assert_eq!(faces.len(), poly.edge_count());

        let total: f64 = faces
            .iter()
            .map(|face| {
                let pts: Vec<[f64; 2]> = face
                    .iter()
                    .map(|&n| {
                        let e = skel.node(n).exact;
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

        // `signed_area2` is twice the polygon's area.
        let want = poly.signed_area2() as f64 / 2.0;
        assert!(
            (total - want).abs() < 0.05 * want.max(1.0),
            "faces cover {total} but the polygon encloses {want}"
        );
    }
}

#[test]
fn regular_ngon_skeleton_converges_on_the_center() {
    // A regular polygon's skeleton is a star: every corner runs to the centre.
    for n in [4usize, 5, 6, 8] {
        let poly = Polygon::from_outer(&ngon(n, 1000.0, 0.0, 0.0)).unwrap();
        let skel = skeleton(&poly).unwrap();

        let peak = skel
            .nodes()
            .iter()
            .max_by(|a, b| a.offset.partial_cmp(&b.offset).unwrap())
            .unwrap();

        // The apothem: distance from centre to each side.
        let apothem = 1000.0 * (std::f64::consts::PI / n as f64).cos();
        assert!(
            (peak.offset as f64 - apothem).abs() < 2.0,
            "{n}-gon: peak offset {} should be near the apothem {apothem}",
            peak.offset
        );
        assert!(peak.position.x.abs() <= 2 && peak.position.y.abs() <= 2);

        check_invariants(&poly, &skel, TOL.max(0.5));
    }
}

/// The first shape with a reflex vertex, so the first that can produce a split
/// event.
#[test]
fn l_shape_is_wellformed() {
    let poly = Polygon::from_outer(&l_shape()).unwrap();
    let skel = skeleton(&poly).unwrap();

    assert_eq!(skel.nodes().iter().filter(|n| n.is_boundary()).count(), 6);
    check_invariants(&poly, &skel, TOL);

    // The reflex elbow at (20, 20) is on the boundary, so its node is there.
    let elbow = skel
        .nodes()
        .iter()
        .find(|n| n.position == Point::new(20, 20) && n.is_boundary());
    assert!(elbow.is_some(), "the reflex elbow needs a boundary node");
}

#[test]
fn plus_shape_is_wellformed() {
    let poly = Polygon::from_outer(&plus_shape()).unwrap();
    let skel = skeleton(&poly).unwrap();

    assert_eq!(skel.nodes().iter().filter(|n| n.is_boundary()).count(), 12);
    check_invariants(&poly, &skel, TOL);
}

/// The skeleton must not care where the polygon sits or which vertex is first.
#[test]
fn skeleton_is_invariant_under_translation() {
    let base = Polygon::from_outer(&l_shape()).unwrap();
    let base_skel = skeleton(&base).unwrap();

    for (dx, dy) in [(100i16, 50i16), (-300, 200), (1000, -1000)] {
        let moved: Vec<Point> = l_shape()
            .iter()
            .map(|p| Point::new(p.x + dx, p.y + dy))
            .collect();
        let poly = Polygon::from_outer(&moved).unwrap();
        let skel = skeleton(&poly).unwrap();

        assert_eq!(skel.node_count(), base_skel.node_count());
        assert_eq!(skel.arc_count(), base_skel.arc_count());

        let mut got: Vec<_> = skel
            .nodes()
            .iter()
            .map(|n| (n.position.x - dx, n.position.y - dy))
            .collect();
        let mut want: Vec<_> = base_skel
            .nodes()
            .iter()
            .map(|n| (n.position.x, n.position.y))
            .collect();
        got.sort();
        want.sort();
        assert_eq!(
            got, want,
            "translating by ({dx}, {dy}) changed the skeleton"
        );
    }
}

#[test]
fn skeleton_is_invariant_under_starting_vertex() {
    let base = Polygon::from_outer(&l_shape()).unwrap();
    let base_skel = skeleton(&base).unwrap();
    let mut want: Vec<_> = base_skel
        .nodes()
        .iter()
        .map(|n| (n.position.x, n.position.y))
        .collect();
    want.sort();

    for rot in 1..6 {
        let mut pts = l_shape();
        pts.rotate_left(rot);
        let poly = Polygon::from_outer(&pts).unwrap();
        let skel = skeleton(&poly).unwrap();

        let mut got: Vec<_> = skel
            .nodes()
            .iter()
            .map(|n| (n.position.x, n.position.y))
            .collect();
        got.sort();
        assert_eq!(
            got, want,
            "rotating the input by {rot} changed the skeleton"
        );
    }
}

#[test]
fn skeleton_is_invariant_under_winding_direction() {
    // Polygon::new normalises winding, so a reversed ring must give the same
    // skeleton geometry.
    let mut reversed = l_shape();
    reversed.reverse();

    let a = skeleton(&Polygon::from_outer(&l_shape()).unwrap()).unwrap();
    let b = skeleton(&Polygon::from_outer(&reversed).unwrap()).unwrap();

    let mut pa: Vec<_> = a
        .nodes()
        .iter()
        .map(|n| (n.position.x, n.position.y))
        .collect();
    let mut pb: Vec<_> = b
        .nodes()
        .iter()
        .map(|n| (n.position.x, n.position.y))
        .collect();
    pa.sort();
    pb.sort();
    assert_eq!(pa, pb);
}

#[test]
fn scaling_the_polygon_scales_the_offsets() {
    let small = Polygon::from_outer(&rect(10, 10)).unwrap();
    let large = Polygon::from_outer(&rect(100, 100)).unwrap();

    let s = skeleton(&small).unwrap();
    let l = skeleton(&large).unwrap();

    assert_eq!(s.node_count(), l.node_count());
    assert!((s.max_offset() * 10.0 - l.max_offset()).abs() < 1e-2);
}

#[test]
fn square_with_a_square_hole() {
    let poly = Polygon::new(
        &rect(30, 30),
        &[vec![
            Point::new(10, 10),
            Point::new(20, 10),
            Point::new(20, 20),
            Point::new(10, 20),
        ]],
    )
    .unwrap();
    let skel = skeleton(&poly).unwrap();

    // Eight boundary nodes: four outer corners, four hole corners.
    assert_eq!(skel.nodes().iter().filter(|n| n.is_boundary()).count(), 8);

    // A hole forces split events: the wavefront from the outer ring must meet
    // the wavefront from the hole, which tears loops apart.
    assert!(
        skel.nodes().iter().any(|n| n.kind == NodeKind::SplitEvent),
        "a hole must produce at least one split event"
    );

    // The material is 10 units thick all round, so nothing exceeds an offset
    // of 5 by much.
    assert!(skel.max_offset() < 5.5, "max offset {}", skel.max_offset());

    check_invariants(&poly, &skel, TOL);
}

#[test]
fn rectangle_with_an_offcentre_hole() {
    let poly = Polygon::new(
        &rect(60, 40),
        &[vec![
            Point::new(10, 10),
            Point::new(25, 10),
            Point::new(25, 25),
            Point::new(10, 25),
        ]],
    )
    .unwrap();
    let skel = skeleton(&poly).unwrap();
    check_invariants(&poly, &skel, TOL);
}

#[test]
fn multiple_holes() {
    let poly = Polygon::new(
        &rect(80, 40),
        &[
            vec![
                Point::new(10, 10),
                Point::new(20, 10),
                Point::new(20, 30),
                Point::new(10, 30),
            ],
            vec![
                Point::new(50, 10),
                Point::new(70, 10),
                Point::new(70, 30),
                Point::new(50, 30),
            ],
        ],
    )
    .unwrap();
    let skel = skeleton(&poly).unwrap();

    assert_eq!(skel.nodes().iter().filter(|n| n.is_boundary()).count(), 12);
    check_invariants(&poly, &skel, TOL);
}

#[test]
fn triangular_hole_in_a_square() {
    let poly = Polygon::new(
        &rect(40, 40),
        &[vec![
            Point::new(15, 15),
            Point::new(25, 15),
            Point::new(20, 25),
        ]],
    )
    .unwrap();
    let skel = skeleton(&poly).unwrap();
    check_invariants(&poly, &skel, TOL);
}

#[test]
fn straight_through_vertices_are_harmless() {
    // Collinear vertices along an edge are legal and must not change the
    // skeleton's shape, only add boundary nodes.
    let plain = Polygon::from_outer(&rect(20, 10)).unwrap();
    let subdivided = Polygon::from_outer(&[
        Point::new(0, 0),
        Point::new(10, 0), // collinear, mid-edge
        Point::new(20, 0),
        Point::new(20, 10),
        Point::new(10, 10), // collinear, mid-edge
        Point::new(0, 10),
    ])
    .unwrap();

    let a = skeleton(&plain).unwrap();
    let b = skeleton(&subdivided).unwrap();
    check_invariants(&subdivided, &b, TOL);

    // Same ridge height either way.
    assert!((a.max_offset() - b.max_offset()).abs() < 1e-3);
    // Each extra input vertex adds a boundary node and its arc.
    assert_eq!(b.nodes().iter().filter(|n| n.is_boundary()).count(), 6);
}

#[test]
fn very_thin_sliver_still_terminates() {
    let poly = Polygon::from_outer(&[
        Point::new(0, 0),
        Point::new(1000, 0),
        Point::new(1000, 1),
        Point::new(0, 1),
    ])
    .unwrap();
    let skel = skeleton(&poly).unwrap();
    assert!((skel.max_offset() - 0.5).abs() < 1e-3);
    check_invariants(&poly, &skel, TOL);
}

/// The far corner of the usable coordinate space, where `f32` has the least
/// absolute resolution (about 0.002) and the arithmetic is under most strain.
#[test]
fn coordinates_at_the_coordinate_cap() {
    let (lo, hi) = (Point::MIN_COORD, Point::MAX_COORD);
    let poly = Polygon::from_outer(&[
        Point::new(lo, lo),
        Point::new(hi, lo),
        Point::new(hi, hi),
        Point::new(lo, hi),
    ])
    .unwrap();
    let skel = skeleton(&poly).unwrap();

    assert_eq!(skel.node_count(), 5);
    let centre = skel.nodes().iter().find(|n| !n.is_boundary()).unwrap();
    assert!(
        centre.position.x.abs() <= 1 && centre.position.y.abs() <= 1,
        "centre landed at {:?}",
        centre.position
    );
    // Half the 32767-wide square.
    assert!(
        (centre.offset - 16383.5).abs() < 1.0,
        "offset {}",
        centre.offset
    );

    // A generous tolerance: f32 resolves only ~0.002 out here, which is the
    // price of the cap and is documented as such.
    check_invariants(&poly, &skel, 0.5);
}

/// A polygon whose coordinates exceed the cap is rejected rather than computed
/// with wrapped predicates.
#[test]
fn coordinates_beyond_the_cap_are_rejected() {
    let e = Polygon::from_outer(&[
        Point::new(-32000, -32000),
        Point::new(32000, -32000),
        Point::new(32000, 32000),
        Point::new(-32000, 32000),
    ])
    .unwrap_err();
    assert!(
        matches!(
            e,
            straight_skeleton::PolygonError::CoordinateOutOfRange { .. }
        ),
        "got {e:?}"
    );
}

#[test]
fn skeleton_of_a_star() {
    // Alternating radii give five sharp points and five deep reflex notches:
    // a workout for split events.
    let mut pts = Vec::new();
    for i in 0..10 {
        let a = std::f64::consts::TAU * (i as f64) / 10.0;
        let r = if i % 2 == 0 { 200.0 } else { 80.0 };
        pts.push(Point::new(
            (r * a.cos()).round() as i16,
            (r * a.sin()).round() as i16,
        ));
    }
    let poly = Polygon::from_outer(&pts).unwrap();
    let skel = skeleton(&poly).unwrap();

    assert_eq!(skel.nodes().iter().filter(|n| n.is_boundary()).count(), 10);
    check_invariants(&poly, &skel, 0.1);
}

#[test]
fn comb_shape_with_many_reflex_vertices() {
    // A comb: teeth pointing up, deep notches between them.
    let mut pts = vec![Point::new(0, 0)];
    for i in 0..5 {
        let x = i * 20;
        pts.push(Point::new(x + 5, 0));
        pts.push(Point::new(x + 5, 30));
        pts.push(Point::new(x + 15, 30));
        pts.push(Point::new(x + 15, 0));
    }
    pts.push(Point::new(100, 0));
    pts.push(Point::new(100, -10));
    pts.push(Point::new(0, -10));

    let poly = Polygon::from_outer(&pts).unwrap();
    let skel = skeleton(&poly).unwrap();
    check_invariants(&poly, &skel, 0.1);
}

#[test]
fn each_boundary_node_emits_exactly_one_arc() {
    for pts in [rect(10, 10), rect(30, 10), l_shape(), plus_shape()] {
        let poly = Polygon::from_outer(&pts).unwrap();
        let skel = skeleton(&poly).unwrap();
        for n in skel.node_ids() {
            if skel.node(n).is_boundary() {
                assert_eq!(skel.arcs_at(n).len(), 1);
            }
        }
    }
}

#[test]
fn arc_sources_always_name_two_distinct_edges() {
    for pts in [rect(20, 10), l_shape(), plus_shape()] {
        let poly = Polygon::from_outer(&pts).unwrap();
        let skel = skeleton(&poly).unwrap();
        for arc in skel.arcs() {
            assert_ne!(arc.sources[0], arc.sources[1]);
            assert!((arc.sources[0].0 as usize) < poly.edge_count());
            assert!((arc.sources[1].0 as usize) < poly.edge_count());
        }
    }
}
