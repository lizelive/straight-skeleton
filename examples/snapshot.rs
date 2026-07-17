//! Dumps skeletons of a fixed corpus of shapes, for diffing across changes.
//!
//! ```text
//! cargo run --release --example snapshot > before.txt
//! ```
//!
//! An optimisation is only allowed to make the crate faster, never to change
//! what it computes. This prints every node and arc of every shape, so the two
//! runs can be compared byte for byte.

use std::fmt::Write as _;
use straight_skeleton::{skeleton, skeleton_constrained, Point, Polygon};

fn main() {
    let mut out = String::new();
    for (name, poly, limits) in corpus() {
        let skel = match &limits {
            None => skeleton(&poly),
            Some(l) => skeleton_constrained(&poly, l),
        };
        let Ok(skel) = skel else {
            let _ = writeln!(out, "{name}: ERROR {:?}", skel.err().unwrap());
            continue;
        };
        let _ = writeln!(
            out,
            "{name}: {} nodes, {} arcs",
            skel.node_count(),
            skel.arc_count()
        );
        for n in skel.nodes() {
            // Bit patterns, not formatted decimals: an optimisation that shifts
            // a result by one ULP has still changed the answer, and rounding to
            // a few decimals would hide it.
            let _ = writeln!(
                out,
                "  N {:?} exact=({:08x},{:08x}) off={:08x} kind={:?} src={:?}",
                n.position,
                n.exact[0].to_bits(),
                n.exact[1].to_bits(),
                n.offset.to_bits(),
                n.kind,
                n.sources
            );
        }
        for a in skel.arcs() {
            let _ = writeln!(out, "  A {:?} src={:?}", a.nodes, a.sources);
        }
        for (i, l) in skel.residual().iter().enumerate() {
            let _ = writeln!(out, "  R{i} nodes={:?} edges={:?}", l.nodes, l.edges);
        }
    }
    print!("{out}");
}

type Case = (String, Polygon, Option<Vec<f32>>);

/// Every shape the test suite and the gallery exercise, plus the bench shapes at
/// a size big enough to hit the interesting paths.
fn corpus() -> Vec<Case> {
    let mut g: Vec<Case> = Vec::new();

    let square = vec![
        Point::new(0, 0),
        Point::new(100, 0),
        Point::new(100, 100),
        Point::new(0, 100),
    ];
    g.push(("square".into(), Polygon::from_outer(&square).unwrap(), None));

    let rect = vec![
        Point::new(0, 0),
        Point::new(200, 0),
        Point::new(200, 100),
        Point::new(0, 100),
    ];
    g.push((
        "rectangle".into(),
        Polygon::from_outer(&rect).unwrap(),
        None,
    ));

    g.push((
        "triangle".into(),
        Polygon::from_outer(&[Point::new(0, 0), Point::new(120, 0), Point::new(0, 90)]).unwrap(),
        None,
    ));

    g.push((
        "l-shape".into(),
        Polygon::from_outer(&[
            Point::new(0, 0),
            Point::new(200, 0),
            Point::new(200, 100),
            Point::new(100, 100),
            Point::new(100, 200),
            Point::new(0, 200),
        ])
        .unwrap(),
        None,
    ));

    g.push((
        "plus".into(),
        Polygon::from_outer(&[
            Point::new(50, 0),
            Point::new(100, 0),
            Point::new(100, 50),
            Point::new(150, 50),
            Point::new(150, 100),
            Point::new(100, 100),
            Point::new(100, 150),
            Point::new(50, 150),
            Point::new(50, 100),
            Point::new(0, 100),
            Point::new(0, 50),
            Point::new(50, 50),
        ])
        .unwrap(),
        None,
    ));

    let mut star = Vec::new();
    for i in 0..10 {
        let a = std::f64::consts::TAU * (i as f64) / 10.0 - std::f64::consts::FRAC_PI_2;
        let r = if i % 2 == 0 { 120.0 } else { 48.0 };
        star.push(Point::new(
            (120.0 + r * a.cos()).round() as i16,
            (120.0 + r * a.sin()).round() as i16,
        ));
    }
    g.push(("star".into(), Polygon::from_outer(&star).unwrap(), None));

    g.push((
        "rect-with-hole".into(),
        Polygon::new(
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
        .unwrap(),
        None,
    ));

    g.push((
        "two-holes".into(),
        Polygon::new(
            &[
                Point::new(0, 0),
                Point::new(240, 0),
                Point::new(240, 120),
                Point::new(0, 120),
            ],
            &[
                vec![
                    Point::new(30, 30),
                    Point::new(80, 30),
                    Point::new(80, 90),
                    Point::new(30, 90),
                ],
                vec![
                    Point::new(150, 40),
                    Point::new(210, 40),
                    Point::new(180, 90),
                ],
            ],
        )
        .unwrap(),
        None,
    ));

    // Constrained cases: the speed-change path and the sliding vertices it makes.
    let p = Polygon::from_outer(&square).unwrap();
    let n = p.edge_count();
    g.push(("square-limit-15".into(), p, Some(vec![15.0; n])));

    let p = Polygon::from_outer(&rect).unwrap();
    let mut limits = vec![f32::INFINITY; p.edge_count()];
    limits[0] = 10.0;
    g.push(("rect-limit-long-edge".into(), p, Some(limits)));

    let p = Polygon::from_outer(&rect).unwrap();
    let mut limits = vec![f32::INFINITY; p.edge_count()];
    limits[1] = 10.0;
    g.push(("rect-limit-short-edge".into(), p, Some(limits)));

    // Limits that bind on every wall, so the wavefront stops rather than
    // collapsing and the result is mostly residual.
    let p = Polygon::from_outer(&[
        Point::new(0, 0),
        Point::new(200, 0),
        Point::new(200, 100),
        Point::new(100, 100),
        Point::new(100, 200),
        Point::new(0, 200),
    ])
    .unwrap();
    let n = p.edge_count();
    g.push(("l-shape-limit-20".into(), p, Some(vec![20.0; n])));

    // Two residual loops, wound opposite ways.
    let p = Polygon::new(
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
    let n = p.edge_count();
    g.push(("rect-with-hole-limit-10".into(), p, Some(vec![10.0; n])));

    // The bench shapes, at a size that exercises many splits and needles.
    g.push((
        "comb-132".into(),
        Polygon::from_outer(&comb(132)).unwrap(),
        None,
    ));
    g.push((
        "star-128".into(),
        Polygon::from_outer(&rand_star(128)).unwrap(),
        None,
    ));

    g
}

fn comb(n: usize) -> Vec<Point> {
    let teeth = (n / 4).max(1);
    let mut pts = vec![Point::new(0, 0)];
    for i in 0..teeth {
        let x = (i as i16) * 20;
        pts.push(Point::new(x + 5, 0));
        pts.push(Point::new(x + 5, 300));
        pts.push(Point::new(x + 15, 300));
        pts.push(Point::new(x + 15, 0));
    }
    pts.push(Point::new(teeth as i16 * 20, 0));
    pts.push(Point::new(teeth as i16 * 20, -40));
    pts.push(Point::new(0, -40));
    pts
}

fn rand_star(n: usize) -> Vec<Point> {
    let mut rng = 0x2545_F491_4F6C_DD1Du64;
    let mut next = || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        (rng >> 33) as f64 / (1u64 << 31) as f64
    };
    (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64) / (n as f64);
            let r = if i % 2 == 0 {
                3000.0 + next() * 1000.0
            } else {
                1200.0 + next() * 400.0
            };
            Point::new((r * a.cos()) as i16, (r * a.sin()) as i16)
        })
        .collect()
}
