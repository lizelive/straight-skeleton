//! Times the skeleton over growing inputs, to measure how it actually scales.
//!
//! ```text
//! cargo run --release --example bench
//! ```
//!
//! Prints the empirical growth exponent: if doubling `n` multiplies the time by
//! `2^k`, then the algorithm behaves like `O(n^k)` over that range. That number,
//! not the theory, is what says whether an optimisation worked.

use std::time::Instant;
use straight_skeleton::{skeleton, Point, Polygon};

fn main() {
    println!("{:>7}  {:>7}  {:>12}  {:>6}", "n", "reflex", "time", "n^k");

    for shape in ["zigzag-comb", "random-star", "convex-ngon"] {
        println!("\n--- {shape}");
        let mut prev: Option<(usize, f64)> = None;

        for &n in &[16usize, 32, 64, 128, 256, 512, 1024] {
            let pts = match shape {
                "zigzag-comb" => comb(n),
                "random-star" => star(n),
                _ => convex(n),
            };
            let Ok(poly) = Polygon::from_outer(&pts) else {
                continue;
            };
            let reflex = poly.vertex_ids().filter(|&v| poly.is_reflex(v)).count();

            let t0 = Instant::now();
            let skel = skeleton(&poly).unwrap();
            let dt = t0.elapsed().as_secs_f64();
            std::hint::black_box(&skel);

            // Empirical exponent against the previous size.
            let k = match prev {
                Some((pn, pt)) if pt > 0.0 && dt > 0.0 => {
                    let s = format!("{:.2}", (dt / pt).log2() / ((n as f64 / pn as f64).log2()));
                    s
                }
                _ => "-".to_string(),
            };
            println!(
                "{:>7}  {:>7}  {:>10.2}ms  {:>6}",
                poly.vertex_count(),
                reflex,
                dt * 1000.0,
                k
            );
            prev = Some((n, dt));

            if dt > 20.0 {
                println!("        (stopping: too slow to keep doubling)");
                break;
            }
        }
    }
}

/// A comb: `n/4` teeth, so about `n/4` reflex vertices. The classic bad case.
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

/// A star with alternating radii: half its vertices are reflex.
fn star(n: usize) -> Vec<Point> {
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

/// A genuinely convex polygon: no reflex vertices, so no split search at all.
///
/// Sampling a circle and rounding to the lattice is not enough — at a few
/// hundred vertices the rounding alone bends some corners the wrong way, and
/// the "convex" case quietly grows reflex vertices and stops measuring what it
/// claims to. So the sampled points go through a convex hull, which by
/// construction cannot.
fn convex(n: usize) -> Vec<Point> {
    let mut pts: Vec<Point> = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64) / (n as f64);
            Point::new((14000.0 * a.cos()) as i16, (14000.0 * a.sin()) as i16)
        })
        .collect();
    pts.sort_by_key(|p| (p.x, p.y));
    pts.dedup();
    hull(&pts)
}

/// Andrew's monotone chain, using the crate's own exact predicate.
fn hull(sorted: &[Point]) -> Vec<Point> {
    use straight_skeleton::predicates::{orient2d, Orientation};
    if sorted.len() < 3 {
        return sorted.to_vec();
    }
    let build = |it: &mut dyn Iterator<Item = &Point>| -> Vec<Point> {
        let mut out: Vec<Point> = Vec::new();
        for &p in it {
            while out.len() >= 2
                && orient2d(out[out.len() - 2], out[out.len() - 1], p)
                    != Orientation::CounterClockwise
            {
                out.pop();
            }
            out.push(p);
        }
        out.pop();
        out
    };
    let mut lower = build(&mut sorted.iter());
    let upper = build(&mut sorted.iter().rev());
    lower.extend(upper);
    lower
}
