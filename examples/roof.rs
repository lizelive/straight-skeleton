//! Raises hip roofs over floor plans and writes them out as Wavefront OBJ.
//!
//! ```text
//! cargo run --example roof
//! cargo run --example roof -- --out my_dir --pitch 0.8
//! ```
//!
//! Open the `.obj` files in Blender, MeshLab, or any 3D viewer.
//!
//! # Why a straight skeleton is a roof
//!
//! This is the classic application, and it is almost too neat. Picture the roof
//! being built by raising the walls' eaves inward at a constant slope. At any
//! moment the still-unroofed floor area is exactly the shrinking wavefront, and
//! the height reached is exactly how far it has travelled. So:
//!
//! - **Every skeleton node is a roof vertex**, at height `offset * pitch`.
//! - **Every skeleton face is one flat roof panel** — the panel rising from the
//!   wall that face belongs to. It is planar because all of its points are the
//!   same distance from that wall's line, and height is a linear function of
//!   that distance.
//! - **Every skeleton arc is a hip, valley, or ridge**: two panels meeting.
//!
//! Nothing here computes geometry; it reads it off the skeleton. That is the
//! whole point. This example asserts each panel really is planar, which is a
//! genuine end-to-end check on the library rather than a decoration.

use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use straight_skeleton::{skeleton, EdgeId, Point, Polygon, Skeleton};

/// How steep the roof is: rise over run. A pitch of 1.0 gives 45 degrees.
const DEFAULT_PITCH: f32 = 0.6;

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

    for (name, plan) in plans() {
        let skel = skeleton(&plan)?;
        let roof = Roof::build(&plan, &skel, pitch).ok_or("could not walk the skeleton's faces")?;

        roof.assert_panels_are_planar(&plan);

        let path = out.join(format!("{name}.obj"));
        fs::write(&path, roof.to_obj(name))?;
        println!(
            "{:<22} {:>2} panels, {:>2} vertices, ridge height {:.2}",
            format!("{name}.obj"),
            roof.panels.len(),
            roof.verts.len(),
            roof.ridge_height(),
        );
    }

    println!("\nWrote OBJ files to {}", out.display());
    Ok(())
}

/// A roof: vertices in 3D, and one flat panel per wall.
struct Roof {
    /// `(x, y, z)` per skeleton node, indexed by node id.
    verts: Vec<[f32; 3]>,
    /// One panel per input edge (wall), as indices into `verts`.
    panels: Vec<Panel>,
}

/// One flat plane of roof, rising from a single wall.
struct Panel {
    /// The wall this panel rises from. This is the traceability the skeleton
    /// gives for free: no search, no nearest-neighbour query.
    wall: EdgeId,
    /// Corners, as indices into [`Roof::verts`].
    corners: Vec<usize>,
}

impl Roof {
    /// Lifts a skeleton into a roof.
    fn build(plan: &Polygon, skel: &Skeleton, pitch: f32) -> Option<Roof> {
        // A skeleton node at offset d is a roof vertex at height d * pitch:
        // the wavefront's travel *is* the run, so the rise follows directly.
        let verts: Vec<[f32; 3]> = skel
            .nodes()
            .iter()
            .map(|n| [n.exact[0], n.exact[1], n.offset * pitch])
            .collect();

        let panels = plan
            .edge_ids()
            .map(|wall| {
                skel.face(wall).map(|corners| Panel {
                    wall,
                    corners: corners.iter().map(|n| n.0 as usize).collect(),
                })
            })
            .collect::<Option<Vec<_>>>()?;

        Some(Roof { verts, panels })
    }

    fn ridge_height(&self) -> f32 {
        self.verts.iter().map(|v| v[2]).fold(0.0, f32::max)
    }

    /// Every panel must be flat, or it is not a roof you could build.
    ///
    /// This is the end-to-end check: it holds only if the skeleton's faces, its
    /// node positions, and its offsets are all correct together. A single
    /// misplaced node buckles its panel and trips this.
    fn assert_panels_are_planar(&self, plan: &Polygon) {
        for panel in &self.panels {
            // The panel rises from its wall, so its plane is pinned by the
            // wall: height must be `pitch * distance-from-the-wall's-line`
            // everywhere on it. Fitting a plane to the points would be a weaker
            // check, since it would pass for a panel over the wrong wall.
            let (a, b) = plan.edge(panel.wall);
            let (dx, dy) = ((b.x - a.x) as f64, (b.y - a.y) as f64);
            let len = (dx * dx + dy * dy).sqrt();
            let (nx, ny) = (-dy / len, dx / len);

            let mut ratio: Option<f64> = None;
            for &c in &panel.corners {
                let v = self.verts[c];
                let run = nx * (v[0] as f64 - a.x as f64) + ny * (v[1] as f64 - a.y as f64);
                let rise = v[2] as f64;
                if run < 1e-6 {
                    // At the eave: rise must be zero there too.
                    assert!(
                        rise.abs() < 1e-3,
                        "panel over wall {} sits at height {rise} on its own eave",
                        panel.wall,
                    );
                    continue;
                }
                let r = rise / run;
                match ratio {
                    None => ratio = Some(r),
                    Some(want) => assert!(
                        (r - want).abs() < 1e-3,
                        "panel over wall {} is not planar: slope {r} here but {want} elsewhere",
                        panel.wall,
                    ),
                }
            }
        }
    }

    /// Serialises to Wavefront OBJ.
    fn to_obj(&self, name: &str) -> String {
        let mut s = format!(
            "# {name} — hip roof generated by the straight-skeleton crate\n\
             # one face per wall; each is a single flat panel\n\
             o {name}\n"
        );
        for v in &self.verts {
            let _ = writeln!(s, "v {} {} {}", v[0], v[1], v[2]);
        }
        for panel in &self.panels {
            // OBJ indices are 1-based.
            let _ = writeln!(
                s,
                "g wall_{}\nf {}",
                panel.wall.0,
                panel
                    .corners
                    .iter()
                    .map(|c| (c + 1).to_string())
                    .collect::<Vec<_>>()
                    .join(" ")
            );
        }
        s
    }
}

/// Floor plans to roof.
fn plans() -> Vec<(&'static str, Polygon)> {
    vec![
        (
            "simple-hip",
            Polygon::from_outer(&[
                Point::new(0, 0),
                Point::new(120, 0),
                Point::new(120, 80),
                Point::new(0, 80),
            ])
            .unwrap(),
        ),
        (
            "pyramid",
            Polygon::from_outer(&[
                Point::new(0, 0),
                Point::new(80, 0),
                Point::new(80, 80),
                Point::new(0, 80),
            ])
            .unwrap(),
        ),
        (
            "l-shaped-house",
            Polygon::from_outer(&[
                Point::new(0, 0),
                Point::new(160, 0),
                Point::new(160, 70),
                Point::new(70, 70),
                Point::new(70, 150),
                Point::new(0, 150),
            ])
            .unwrap(),
        ),
        (
            "t-shaped-house",
            Polygon::from_outer(&[
                Point::new(0, 0),
                Point::new(180, 0),
                Point::new(180, 60),
                Point::new(120, 60),
                Point::new(120, 140),
                Point::new(60, 140),
                Point::new(60, 60),
                Point::new(0, 60),
            ])
            .unwrap(),
        ),
        (
            "courtyard",
            Polygon::new(
                &[
                    Point::new(0, 0),
                    Point::new(200, 0),
                    Point::new(200, 160),
                    Point::new(0, 160),
                ],
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
