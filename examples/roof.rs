//! Raises roofs over floor plans and writes them out as Wavefront OBJ.
//!
//! ```text
//! cargo run --example roof
//! cargo run --example roof -- --out my_dir --pitch 0.8
//! ```
//!
//! Open the `.obj` files in Blender, MeshLab, or any 3D viewer.
//!
//! Each plan is roofed four ways, which is the point of the example: **the
//! skeleton is the same every time**. A straight skeleton gives a roof's *plan*
//! — where the hips, valleys and ridges run — and that does not depend on how
//! high anything is. So all four styles read off one skeleton, and differ only
//! in their [`Profile`]:
//!
//! | file | what it is |
//! |---|---|
//! | `*-hip.obj` | the classic: one slope all the way to the ridge |
//! | `*-mansard.obj` | steep to the kerb, then shallow — two panels per wall |
//! | `*-truncated.obj` | a hip stopped short, leaving a flat on top |
//! | `*-truncated-mansard.obj` | steep, then shallow, then flat |
//! | `*-gambrel.obj` | a mansard with **gable** ends: the barn roof |
//!
//! The last three need a **constrained** skeleton, and the two constraints do
//! opposite jobs. A **uniform** limit stops the whole wavefront before it
//! collapses, and the flat on top is that [residual] raised to the limit's
//! height. A **zero** limit on one wall stops that wall alone: it sweeps
//! nothing, so its panel is the degenerate face stood on end — a gable — and
//! the ridge runs out to it rather than hipping away.
//!
//! Only plans with a ridge to run out to get a gambrel, so the square and the
//! courtyard skip it.
//!
//! The geometry itself lives in the library — see [`Roof`]. All this example
//! adds is the OBJ serialisation, since choosing a file format is the caller's
//! business, not the crate's, and the triangulation that goes with it: a
//! courtyard's flat has a hole in it, and OBJ has no way to say so. See
//! [`flat_triangles`], which is here rather than in the crate for the same
//! reason — choosing a mesh representation is downstream of computing one.
//!
//! [residual]: straight_skeleton::Skeleton::residual

use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use straight_skeleton::{
    skeleton, skeleton_constrained, Panel, PanelKind, Point, Polygon, Profile, Roof, RoofVertexId,
};

/// How steep the roof is: rise over run. A pitch of 1.0 gives 45 degrees.
const DEFAULT_PITCH: f32 = 0.6;

/// How far up a mansard's steep lower slope runs before the pitch breaks.
const MANSARD_BREAK: f32 = 12.0;

/// The steep lower pitch of a mansard. Much steeper than the hip's, which is
/// the whole point: it buys headroom in the storey inside the roof.
const MANSARD_LOWER: f32 = 2.2;

/// The shallow upper pitch of a mansard, which keeps it from getting silly.
const MANSARD_UPPER: f32 = 0.22;

/// How far the wavefront runs before a truncated roof stops it, leaving a flat.
const TRUNCATE_AT: f32 = 22.0;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let mut out = PathBuf::from("target/roofs");
    if let Some(i) = args.iter().position(|a| a == "--out") {
        out = PathBuf::from(args.get(i + 1).ok_or("--out needs a directory")?);
    }
    let pitch: f32 = match args.iter().position(|a| a == "--pitch") {
        Some(i) => args.get(i + 1).ok_or("--pitch needs a number")?.parse()?,
        None => DEFAULT_PITCH,
    };
    fs::create_dir_all(&out)?;

    let mansard = Profile::Mansard {
        lower_pitch: MANSARD_LOWER,
        break_offset: MANSARD_BREAK,
        upper_pitch: MANSARD_UPPER,
    };

    println!(
        "{:<34} {:>6} {:>8} {:>6} {:>7}",
        "file", "panels", "vertices", "flats", "ridge"
    );

    for (name, plan, gable_walls) in plans() {
        // One skeleton per plan. Every style below reads off this same one.
        let skel = skeleton(&plan)?;

        // ...and one constrained skeleton, whose wavefront stops at the limit
        // instead of collapsing. A uniform limit is what makes it a roof: with
        // uneven ones the surface would have to tear, and `Roof` says so rather
        // than inventing something. See `RoofError::UnevenLimits`.
        let stopped = skeleton_constrained(&plan, &vec![TRUNCATE_AT; plan.edge_count()])?;

        let mut styles: Vec<(&str, Roof)> = vec![
            ("hip", Roof::with_profile(&skel, Profile::Hip { pitch })?),
            ("mansard", Roof::with_profile(&skel, mansard)?),
            (
                "truncated",
                Roof::with_profile(&stopped, Profile::Hip { pitch })?,
            ),
            ("truncated-mansard", Roof::with_profile(&stopped, mansard)?),
        ];

        // A **gambrel**: a mansard with gable ends, which is the classic barn
        // roof. Both halves of that come free.
        //
        // A limit of **zero** is what makes a gable. The wall never moves, so it
        // sweeps nothing, so its face is degenerate in plan — and standing that
        // degenerate face up at `z = offset * pitch` gives the vertical triangle
        // that *is* the gable. The ridge then runs right out to the wall instead
        // of hipping away from it, because the neighbouring walls' corners slide
        // **along** the frozen wall rather than over it.
        //
        // Mixing zero limits with unlimited ones is allowed where mixing two
        // *finite* limits is not, and the difference is not arbitrary: a wall
        // that never moves has no sloping panel to be inconsistent about. See
        // `RoofError::UnevenLimits`.
        if !gable_walls.is_empty() {
            let mut limits = vec![f32::INFINITY; plan.edge_count()];
            for &w in gable_walls {
                limits[w as usize] = 0.0;
            }
            let gabled = skeleton_constrained(&plan, &limits)?;
            styles.push(("gambrel", Roof::with_profile(&gabled, mansard)?));
        }

        for (style, roof) in styles {
            let stem = format!("{name}-{style}");
            let path = out.join(format!("{stem}.obj"));
            fs::write(&path, to_obj(&roof, &stem))?;
            println!(
                "{:<34} {:>6} {:>8} {:>6} {:>7}",
                format!("{stem}.obj"),
                roof.panels().len(),
                roof.verts().len(),
                roof.flat().count(),
                roof.ridge_height(),
            );
        }
    }

    println!("\nWrote OBJ files to {}", out.display());
    Ok(())
}

/// Twice the signed area of a panel's footprint, halved: the shoelace formula
/// on its `x`/`y`. Positive when the loop winds counter-clockwise.
fn footprint_area(roof: &Roof, panel: &Panel) -> f64 {
    let mut a = 0.0;
    for i in 0..panel.corners.len() {
        let p = roof.vertex(panel.corners[i]).exact;
        let q = roof
            .vertex(panel.corners[(i + 1) % panel.corners.len()])
            .exact;
        a += p[0] as f64 * q[1] as f64 - q[0] as f64 * p[1] as f64;
    }
    a / 2.0
}

/// Whether a point falls inside a panel's footprint, by crossing number.
///
/// Only used to decide which outline a hole belongs to, and the two are nested
/// with real clearance between them, so this does not need to be careful about
/// points exactly on the boundary.
fn footprint_contains(roof: &Roof, panel: &Panel, p: [f32; 3]) -> bool {
    let mut inside = false;
    for i in 0..panel.corners.len() {
        let a = roof.vertex(panel.corners[i]).exact;
        let b = roof
            .vertex(panel.corners[(i + 1) % panel.corners.len()])
            .exact;
        if (a[1] > p[1]) != (b[1] > p[1])
            && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]
        {
            inside = !inside;
        }
    }
    inside
}

/// Cuts the roof's flats into triangles.
///
/// # Why this is here at all
///
/// OBJ has no way to say "this face has a hole in it", and a courtyard's flat
/// has exactly that: the wavefront stops before it reaches the middle, so what
/// is left is a ring. `Roof::flat` returns it the way the crate returns
/// everything with holes — as several loops, an outline wound counter-clockwise
/// and a hole wound clockwise, like a `Polygon`'s rings. Emitted naively that
/// becomes a solid slab paved straight over the courtyard.
///
/// # Why it is not in the library
///
/// Triangulating is choosing a mesh representation, which is the caller's
/// business in the same way that choosing a file format is. The crate hands over
/// the loops and their winding, which is the part only it knows; turning those
/// into whatever a consumer's renderer wants is downstream of that. So this
/// leans on `earcutr` — a well-tested ear clipper — rather than the crate
/// growing a triangulator, or this example growing a half-tested one.
fn flat_triangles(roof: &Roof) -> Vec<[RoofVertexId; 3]> {
    let flats: Vec<&Panel> = roof.flat().collect();
    // Winding is what says which is which: an outline encloses area, a hole
    // un-encloses it.
    let (outlines, holes): (Vec<&Panel>, Vec<&Panel>) =
        flats.iter().partition(|p| footprint_area(roof, p) > 0.0);

    let mut tris = Vec::new();
    for outline in outlines {
        // A flat can have several holes, and a roof several flats, so pair them
        // up by containment rather than assuming there is only ever one of each.
        let mine: Vec<&&Panel> = holes
            .iter()
            .filter(|h| footprint_contains(roof, outline, roof.vertex(h.corners[0]).exact))
            .collect();

        let mut coords: Vec<f64> = Vec::new();
        let mut ids: Vec<RoofVertexId> = Vec::new();
        let mut hole_starts: Vec<usize> = Vec::new();

        let push = |panel: &Panel, coords: &mut Vec<f64>, ids: &mut Vec<RoofVertexId>| {
            for &c in &panel.corners {
                let e = roof.vertex(c).exact;
                coords.push(e[0] as f64);
                coords.push(e[1] as f64);
                ids.push(c);
            }
        };
        push(outline, &mut coords, &mut ids);
        for h in &mine {
            hole_starts.push(ids.len());
            push(h, &mut coords, &mut ids);
        }

        // A flat is level, so triangulating its footprint in 2D triangulates the
        // flat itself; the ids carry the height back.
        let flat_tris = earcutr::earcut(&coords, &hole_starts, 2)
            .expect("a flat is a simple polygon with simple holes, which earcut always handles");
        for t in flat_tris.chunks(3) {
            tris.push([ids[t[0]], ids[t[1]], ids[t[2]]]);
        }

        // The triangles must cover the outline minus its holes, and no more.
        // Signed areas make that one subtraction: an outline encloses, a hole
        // un-encloses. This is what would catch a hole paired to the wrong
        // outline — which is the only judgement call in here, and the one that
        // would otherwise show up as a courtyard paved over.
        let want: f64 = footprint_area(roof, outline)
            + mine.iter().map(|h| footprint_area(roof, h)).sum::<f64>();
        let got: f64 = tris
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|v| roof.vertex(v).exact);
                let ab = [(b[0] - a[0]) as f64, (b[1] - a[1]) as f64];
                let ac = [(c[0] - a[0]) as f64, (c[1] - a[1]) as f64];
                (ab[0] * ac[1] - ab[1] * ac[0]) / 2.0
            })
            .sum();
        assert!(
            (got - want).abs() < 1e-3 * want.abs().max(1.0),
            "the triangles cover {got} but the flat is {want}"
        );
    }
    tris
}

/// Serialises a roof to Wavefront OBJ.
///
/// Sloping panels go out as single faces: OBJ handles n-gons, and a panel is
/// flat by construction, so there is nothing to gain by cutting them up. The
/// flats are triangulated, because they can have holes and OBJ cannot say so —
/// see [`flat_triangles`].
fn to_obj(roof: &Roof, name: &str) -> String {
    let mut s = format!("# {name} - generated by the straight-skeleton crate\n");
    let _ = match roof.profile() {
        Profile::Hip { pitch } => writeln!(s, "# hip profile, pitch {pitch}"),
        Profile::Mansard {
            lower_pitch,
            break_offset,
            upper_pitch,
        } => writeln!(
            s,
            "# mansard profile: pitch {lower_pitch} to a kerb at offset \
             {break_offset}, then {upper_pitch}"
        ),
    };
    let flats = roof.flat().count();
    if flats > 0 {
        let _ = writeln!(s, "# truncated: {flats} flat loop(s) on top, triangulated");
    }
    let _ = writeln!(s, "o {name}");

    // Vertices go out in `RoofVertexId` order, so a panel's corner ids are
    // already the right OBJ indices, offset by one.
    for v in roof.verts() {
        let p = v.position;
        let _ = writeln!(s, "v {} {} {}", p.x, p.y, p.z);
    }

    // OBJ indices are 1-based.
    let obj = |c: &RoofVertexId| (c.0 + 1).to_string();

    for panel in roof.panels() {
        // Name each group after where it came from — the provenance the skeleton
        // hands over for free.
        let PanelKind::Slope { wall, band } = panel.kind else {
            continue; // the flats go out below, as triangles
        };
        // A wall given a zero limit never moves, so it sweeps no plan area at
        // all: its panel is the degenerate face standing on edge, which is a
        // gable. Same `Slope` panel, same `z = height_at(offset)`, but worth
        // naming for what it is rather than calling it a slope that happens to
        // be vertical.
        let _ = if footprint_area(roof, panel).abs() < 1e-3 {
            writeln!(s, "g gable_{}_band_{band}", wall.0)
        } else {
            writeln!(s, "g wall_{}_band_{band}", wall.0)
        };
        let corners = panel.corners.iter().map(obj).collect::<Vec<_>>().join(" ");
        let _ = writeln!(s, "f {corners}");
    }

    let tris = flat_triangles(roof);
    if !tris.is_empty() {
        let _ = writeln!(s, "g flat");
        for t in tris {
            let _ = writeln!(s, "f {} {} {}", obj(&t[0]), obj(&t[1]), obj(&t[2]));
        }
    }
    s
}

/// Floor plans to roof, each with the walls to gable — see the `gambrel` style
/// in `main`, and leave it empty to skip that style for a plan.
///
/// The plans are deliberately larger than the crate's test shapes: a mansard's
/// kerb and a truncated roof's flat both need room to sit between the eaves and
/// the ridge, and on a small plan the limits would clear the ridge and never
/// bite.
fn plans() -> Vec<(&'static str, Polygon, &'static [u16])> {
    vec![
        (
            "simple",
            Polygon::from_outer(&[
                Point::new(0, 0),
                Point::new(160, 0),
                Point::new(160, 110),
                Point::new(0, 110),
            ])
            .unwrap(),
            // The two 110-long ends, either side of the ridge that runs the
            // 160-long way. Gabling them turns a hipped rectangle into the
            // plainest possible barn.
            &[1, 3],
        ),
        (
            "pyramid",
            Polygon::from_outer(&[
                Point::new(0, 0),
                Point::new(120, 0),
                Point::new(120, 120),
                Point::new(0, 120),
            ])
            .unwrap(),
            // A square has no ridge to run out to a gable — its skeleton meets
            // at a point — so there is no end wall to pick, and no gambrel.
            &[],
        ),
        (
            "l-shaped-house",
            Polygon::from_outer(&[
                Point::new(0, 0),
                Point::new(220, 0),
                Point::new(220, 100),
                Point::new(100, 100),
                Point::new(100, 210),
                Point::new(0, 210),
            ])
            .unwrap(),
            // e1 is the right end of the long wing, e4 the top of the short one:
            // the two walls each of the L's ridges runs out to.
            &[1, 4],
        ),
        (
            "t-shaped-house",
            Polygon::from_outer(&[
                Point::new(0, 0),   // e0: 240 along the bottom
                Point::new(240, 0), // e1: the crossbar's right end, 90
                Point::new(240, 90),
                Point::new(165, 90),  // e2: the right shoulder, 75
                Point::new(165, 200), // e3: the stem's right side, 110
                Point::new(75, 200),  // e4: the stem's end, 90
                Point::new(75, 90),   // e5: the stem's left side, 110
                Point::new(0, 90),    // e6: the left shoulder, 75
            ])
            .unwrap(),
            // The three short end walls, one at each end of a ridge: the
            // crossbar's ridge runs east-west and ends at e1 and e7, the stem's
            // runs north-south and ends at e4.
            //
            // Not the *shortest* walls — the shoulders e2 and e6 are 75 to these
            // 90. But a shoulder faces into the crossbar rather than capping a
            // ridge, so freezing one gives a lean-to rather than a gable, and
            // drags the crossbar's roof up to offset 90 instead of meeting at 45.
            &[1, 4, 7],
        ),
        (
            "courtyard",
            Polygon::new(
                &[
                    Point::new(0, 0),
                    Point::new(260, 0),
                    Point::new(260, 200),
                    Point::new(0, 200),
                ],
                &[vec![
                    Point::new(90, 70),
                    Point::new(170, 70),
                    Point::new(170, 130),
                    Point::new(90, 130),
                ]],
            )
            .unwrap(),
            // A courtyard's roof is a closed ring with no ridge end anywhere, so
            // again nothing to gable.
            &[],
        ),
    ]
}
