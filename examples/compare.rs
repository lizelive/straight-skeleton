//! Diffs two snapshot files geometrically rather than textually.
//!
//! ```text
//! cargo run --release --example compare -- before.txt after.txt
//! ```
//!
//! The simulation runs in `f32`, so any change to the order arithmetic happens
//! in moves results by an ULP or two, and a textual diff drowns in that. This
//! matches each node in one file to the nearest node in the other carrying the
//! same sources, and reports the largest gap it had to bridge. A reordering or
//! an ULP of drift comes out as a tiny number; a genuinely different skeleton
//! does not.

use std::collections::HashMap;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: compare <before.txt> <after.txt>");
        std::process::exit(2);
    }
    let a = parse(&args[1]);
    let b = parse(&args[2]);

    let mut worst_overall = 0.0f64;
    let mut failed = false;

    let mut names: Vec<&String> = a.keys().collect();
    names.sort();
    for name in names {
        let na = &a[name];
        let Some(nb) = b.get(name) else {
            println!("{name}: MISSING from second file");
            failed = true;
            continue;
        };
        if na.len() != nb.len() {
            println!("{name}: NODE COUNT {} -> {}", na.len(), nb.len());
            failed = true;
            continue;
        }

        // Greedy nearest match among nodes with identical sources. Sources are a
        // combinatorial fact, not a numeric one, so they must agree exactly;
        // only the positions are allowed to wobble.
        let mut used = vec![false; nb.len()];
        let mut worst = 0.0f64;
        let mut unmatched = 0;
        for x in na {
            let mut best: Option<(f64, usize)> = None;
            for (j, y) in nb.iter().enumerate() {
                if used[j] || x.sources != y.sources {
                    continue;
                }
                let d = ((x.pos.0 - y.pos.0).powi(2) + (x.pos.1 - y.pos.1).powi(2)).sqrt();
                let d = d.max((x.offset - y.offset).abs());
                if best.map_or(true, |(bd, _)| d < bd) {
                    best = Some((d, j));
                }
            }
            match best {
                Some((d, j)) => {
                    used[j] = true;
                    worst = worst.max(d);
                }
                None => unmatched += 1,
            }
        }
        worst_overall = worst_overall.max(worst);
        if unmatched > 0 {
            println!("{name}: {unmatched} node(s) with no counterpart at all");
            failed = true;
        } else {
            println!(
                "{name}: {} nodes matched, worst drift {worst:.2e}",
                na.len()
            );
        }
    }

    println!("---");
    println!("worst drift across every shape: {worst_overall:.3e}");
    if failed {
        println!("VERDICT: shapes differ structurally");
        std::process::exit(1);
    }
    // One lattice unit is the resolution the public API reports positions at, so
    // drift orders of magnitude below it cannot change what a caller sees.
    if worst_overall < 1e-2 {
        println!("VERDICT: identical up to f32 noise");
    } else {
        println!("VERDICT: drift exceeds f32 noise - investigate");
        std::process::exit(1);
    }
}

struct Node {
    pos: (f64, f64),
    offset: f64,
    sources: String,
}

fn parse(path: &str) -> HashMap<String, Vec<Node>> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut out: HashMap<String, Vec<Node>> = HashMap::new();
    let mut cur = String::new();
    for line in text.lines() {
        let line = line.trim_end();
        if let Some(colon) = line.find(": ") {
            if !line.starts_with(' ') && line.contains(" nodes, ") {
                cur = line[..colon].to_string();
                out.insert(cur.clone(), Vec::new());
                continue;
            }
        }
        let t = line.trim_start();
        if let Some(rest) = t.strip_prefix("N ") {
            let exact = rest
                .split("exact=(")
                .nth(1)
                .and_then(|s| s.split(')').next())
                .expect("exact field");
            let (xh, yh) = exact.split_once(',').expect("two components");
            let off = rest
                .split("off=")
                .nth(1)
                .and_then(|s| s.split(' ').next())
                .expect("off field");
            let src = rest.split("src=").nth(1).expect("src field").to_string();
            out.get_mut(&cur).expect("node before header").push(Node {
                pos: (bits(xh) as f64, bits(yh) as f64),
                offset: bits(off) as f64,
                sources: src,
            });
        }
    }
    out
}

fn bits(hex: &str) -> f32 {
    f32::from_bits(u32::from_str_radix(hex.trim(), 16).expect("hex float"))
}
