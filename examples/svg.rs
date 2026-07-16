//! Renders straight skeletons to SVG so you can look at them.
//!
//! ```text
//! cargo run --example svg
//! cargo run --example svg -- --out my_dir
//! ```
//!
//! Writes one SVG per shape into `target/svg/` (or `--out`), plus an
//! `index.html` that shows them all on one page. Open that file in a browser.
//!
//! Each drawing shows:
//!
//! - the **input polygon** as a heavy dark outline,
//! - the **skeleton arcs**, coloured by which input edge they belong to, so you
//!   can see the provenance from [`Arc::sources`] directly,
//! - the **skeleton faces** as translucent fills — one face per input edge,
//! - **nodes**, sized by kind, labelled with their offset.

use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use straight_skeleton::{
    skeleton, skeleton_constrained, EdgeId, NodeKind, Point, Polygon, Skeleton,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut out = PathBuf::from("target/svg");
    let args: Vec<String> = std::env::args().collect();
    if let Some(i) = args.iter().position(|a| a == "--out") {
        out = PathBuf::from(args.get(i + 1).ok_or("--out needs a directory")?);
    }
    fs::create_dir_all(&out)?;

    let mut written: Vec<(String, String)> = Vec::new();

    for (name, poly, limits) in gallery() {
        let skel = match &limits {
            Some(l) => skeleton_constrained(&poly, l)?,
            None => skeleton(&poly)?,
        };

        let svg = render(&poly, &skel);
        let file = format!("{name}.svg");
        fs::write(out.join(&file), &svg)?;

        let note = format!(
            "{} vertices, {} nodes, {} arcs, ridge {:.2}{}",
            poly.vertex_count(),
            skel.node_count(),
            skel.arc_count(),
            skel.max_offset(),
            if limits.is_some() {
                " (constrained)"
            } else {
                ""
            },
        );
        println!("{:<28} {}", file, note);
        written.push((file, format!("{name} — {note}")));
    }

    let index = out.join("index.html");
    fs::write(&index, index_html(&written))?;
    println!("\nOpen {} to view them all.", display(&index));
    Ok(())
}

/// The shapes to draw, each with optional per-edge distance limits.
#[allow(clippy::type_complexity)]
fn gallery() -> Vec<(String, Polygon, Option<Vec<f32>>)> {
    let mut g: Vec<(String, Polygon, Option<Vec<f32>>)> = Vec::new();

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
        "rectangle-ridge".into(),
        Polygon::from_outer(&rect).unwrap(),
        None,
    ));

    g.push((
        "triangle".into(),
        Polygon::from_outer(&[Point::new(0, 0), Point::new(160, 0), Point::new(40, 120)]).unwrap(),
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

    // A five-pointed star: sharp tips and deep reflex notches.
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

    let with_hole = Polygon::new(
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
    g.push(("rect-with-hole".into(), with_hole, None));

    let two_holes = Polygon::new(
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
    .unwrap();
    g.push(("two-holes".into(), two_holes, None));

    // Constrained: every edge stops at 15, truncating the skeleton.
    let p = Polygon::from_outer(&square).unwrap();
    let n = p.edge_count();
    g.push(("square-limit-15".into(), p, Some(vec![15.0; n])));

    // Constrained, non-uniformly: the bottom edge is held back at 10 while the
    // rest run free.
    let p = Polygon::from_outer(&rect).unwrap();
    let mut limits = vec![f32::INFINITY; p.edge_count()];
    limits[0] = 10.0;
    g.push(("rectangle-limit-one-edge".into(), p, Some(limits)));

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

    g
}

/// A stable, readable colour per input edge.
fn edge_colour(e: EdgeId, total: usize) -> String {
    // Even spacing round the hue wheel, at fixed saturation and lightness, so
    // adjacent faces are always distinguishable.
    let hue = 360.0 * (e.0 as f64) / (total.max(1) as f64);
    format!("hsl({hue:.0} 70% 45%)")
}

fn render(poly: &Polygon, skel: &Skeleton) -> String {
    let (min_x, min_y, max_x, max_y) = bounds(poly);
    let pad = 16.0;
    let w = (max_x - min_x) as f64 + 2.0 * pad;
    let h = (max_y - min_y) as f64 + 2.0 * pad;
    let ox = min_x as f64 - pad;
    let oy = min_y as f64 - pad;
    let ne = poly.edge_count();

    let mut s = String::new();
    // SVG's y axis points down; flip it so the drawing matches the maths.
    let _ = write!(
        s,
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="{ox} {oy} {w} {h}" width="{w}" height="{h}">
<g transform="translate(0 {ty}) scale(1 -1)">
<rect x="{ox}" y="{oy}" width="{w}" height="{h}" fill="#fbfbfd"/>"##,
        ty = 2.0 * oy + h,
    );

    // Faces, translucent, one per input edge.
    if let Some(faces) = skel.faces() {
        for (i, face) in faces.iter().enumerate() {
            let pts: Vec<String> = face
                .iter()
                .map(|&n| {
                    let p = skel.node(n).exact;
                    format!("{:.2},{:.2}", p[0], p[1])
                })
                .collect();
            let _ = writeln!(
                s,
                r##"<polygon points="{}" fill="{}" fill-opacity="0.14" stroke="none"/>"##,
                pts.join(" "),
                edge_colour(EdgeId(i as u16), ne),
            );
        }
    }

    // Skeleton arcs, coloured by their first source edge.
    for a in skel.arc_ids() {
        let arc = skel.arc(a);
        let p0 = skel.node(arc.nodes[0]).exact;
        let p1 = skel.node(arc.nodes[1]).exact;
        let _ = writeln!(
            s,
            r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{}" stroke-width="1.1" stroke-linecap="round"/>"##,
            p0[0],
            p0[1],
            p1[0],
            p1[1],
            edge_colour(arc.sources[0], ne),
        );
    }

    // The input polygon, on top, heavy.
    for ring in poly.rings() {
        let pts: Vec<String> = ring.iter().map(|p| format!("{},{}", p.x, p.y)).collect();
        let _ = writeln!(
            s,
            r##"<polygon points="{}" fill="none" stroke="#16161d" stroke-width="1.8"/>"##,
            pts.join(" ")
        );
    }

    // Nodes.
    for n in skel.node_ids() {
        let node = skel.node(n);
        let (r, fill) = match node.kind {
            NodeKind::Boundary(_) => (1.6, "#16161d"),
            NodeKind::EdgeEvent => (2.0, "#2f6fed"),
            NodeKind::SplitEvent => (2.6, "#e0463c"),
            NodeKind::LimitReached => (2.6, "#0f9d58"),
        };
        let _ = writeln!(
            s,
            r##"<circle cx="{:.2}" cy="{:.2}" r="{r}" fill="{fill}"/>"##,
            node.exact[0], node.exact[1],
        );
    }

    s.push_str("</g>\n");

    // Legend, drawn outside the flipped group so the text is not upside down.
    let _ = writeln!(
        s,
        r##"<g font-family="ui-monospace,monospace" font-size="7" fill="#55555f">
<circle cx="{lx}" cy="{l0}" r="2" fill="#2f6fed"/><text x="{tx}" y="{t0}">edge event</text>
<circle cx="{lx}" cy="{l1}" r="2.6" fill="#e0463c"/><text x="{tx}" y="{t1}">split event</text>
<circle cx="{lx}" cy="{l2}" r="2.6" fill="#0f9d58"/><text x="{tx}" y="{t2}">limit reached</text>
</g>
</svg>"##,
        lx = ox + 6.0,
        tx = ox + 12.0,
        l0 = oy + 8.0,
        t0 = oy + 10.5,
        l1 = oy + 18.0,
        t1 = oy + 20.5,
        l2 = oy + 28.0,
        t2 = oy + 30.5,
    );
    s
}

fn bounds(poly: &Polygon) -> (i16, i16, i16, i16) {
    let vs = poly.vertices();
    let min_x = vs.iter().map(|p| p.x).min().unwrap_or(0);
    let min_y = vs.iter().map(|p| p.y).min().unwrap_or(0);
    let max_x = vs.iter().map(|p| p.x).max().unwrap_or(0);
    let max_y = vs.iter().map(|p| p.y).max().unwrap_or(0);
    (min_x, min_y, max_x, max_y)
}

fn index_html(items: &[(String, String)]) -> String {
    let mut s = String::from(
        "<!doctype html><meta charset=utf-8><title>straight-skeleton gallery</title>\n\
         <style>body{font:14px system-ui;margin:2rem;background:#fff;color:#16161d}\
         .g{display:flex;flex-wrap:wrap;gap:1.5rem}\
         figure{margin:0}img{border:1px solid #e5e5ea;border-radius:6px;display:block}\
         figcaption{font:12px ui-monospace,monospace;color:#55555f;margin-top:.4rem;max-width:34ch}\
         </style>\n<h1>straight-skeleton</h1><div class=g>\n",
    );
    for (file, note) in items {
        let _ = writeln!(
            s,
            "<figure><img src=\"{file}\" alt=\"{file}\"><figcaption>{note}</figcaption></figure>\n"
        );
    }
    s.push_str("</div>\n");
    s
}

fn display(p: &Path) -> String {
    p.display().to_string().replace('\\', "/")
}
