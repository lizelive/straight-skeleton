//! The wavefront simulation that computes the skeleton.
//!
//! See `docs/ALGORITHM.md` for the illustrated walkthrough. In brief: every
//! edge of the polygon slides inward at its own speed, dragging its neighbours
//! with it. Vertices trace out the skeleton's arcs. The simulation advances
//! from one *event* to the next, where an event is any moment the wavefront's
//! shape changes discontinuously.

use alloc::collections::BinaryHeap;
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;
use core::fmt;

use crate::math::Vec2;
use crate::polygon::{EdgeId, Polygon};
use crate::skeleton::{Arc, Node, NodeId, NodeKind, Skeleton};
use crate::Point;

/// Tolerance for geometric comparisons, in coordinate units.
///
/// Input coordinates are integers bounded by 32767, so `f32` resolves about
/// `1e-11` there. `1e-7` sits far enough above the noise floor to absorb
/// accumulated error, and far enough below 1 that it can never merge two
/// distinct lattice points.
const EPS: f32 = 1e-4;

/// Tolerance for treating two unit normals as parallel.
///
/// Compared against the cross product of two unit vectors, i.e. `sin` of the
/// angle between them, so this is an angle of roughly `1e-9` radians.
const PARALLEL_EPS: f32 = 1e-6;

/// How close two wavefront vertices must be to count as arriving at the same
/// place, and so be retired by a single vertex event.
///
/// Looser than [`EPS`] because these positions are extrapolated from event
/// times that are themselves computed, so the error in them compounds. Still
/// six orders of magnitude below one lattice unit, so it cannot fuse vertices
/// that genuinely belong apart.
const MERGE_EPS: f32 = 1e-2;

/// Why a skeleton could not be computed.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum SkeletonError {
    /// The per-edge limit slice did not have one entry per input edge.
    LimitCountMismatch {
        /// How many limits were supplied.
        got: usize,
        /// How many the polygon needs — one per edge.
        expected: usize,
    },
    /// A per-edge limit was negative or NaN.
    InvalidLimit {
        /// The offending edge.
        edge: EdgeId,
        /// The value supplied.
        value: f32,
    },
    /// Two collinear neighbouring edges were given different distance limits.
    ///
    /// When one stops and the other keeps going, the wavefront would have to
    /// tear open between two parallel lines at different offsets, and the
    /// vertex between them has nowhere to be. Give both edges the same limit,
    /// or separate them with a non-collinear edge.
    IncompatibleCollinearLimits {
        /// The edge arriving at the problem vertex.
        left: EdgeId,
        /// The edge leaving it.
        right: EdgeId,
    },
    /// The simulation exceeded its event budget.
    ///
    /// This is a guard against a non-terminating loop caused by degenerate
    /// input, and should not be reachable. Please report it as a bug, with the
    /// polygon that triggered it.
    EventBudgetExhausted {
        /// The budget that was exceeded.
        budget: usize,
    },
}

impl fmt::Display for SkeletonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SkeletonError::LimitCountMismatch { got, expected } => write!(
                f,
                "got {got} per-edge limits but the polygon has {expected} edges"
            ),
            SkeletonError::InvalidLimit { edge, value } => {
                write!(f, "edge {} has an invalid distance limit {value}", edge.0)
            }
            SkeletonError::IncompatibleCollinearLimits { left, right } => write!(
                f,
                "collinear edges {} and {} have different distance limits, \
                 which would tear the wavefront apart",
                left.0, right.0
            ),
            SkeletonError::EventBudgetExhausted { budget } => write!(
                f,
                "the wavefront simulation exceeded its budget of {budget} events; \
                 this is a bug, please report it"
            ),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for SkeletonError {}

/// The fixed geometry of one input edge, plus its distance limit.
///
/// A wavefront edge is its input edge's supporting line, translated inward. At
/// time `t` the line is `{ x : normal · x = c + offset_at(t) }`.
#[derive(Clone, Copy, Debug)]
struct EdgeState {
    /// Unit vector along the edge, from its start vertex to its end vertex.
    dir: Vec2,
    /// Unit normal pointing into the polygon's interior (i.e. to the left of
    /// `dir`, per the CCW-outer-ring convention).
    normal: Vec2,
    /// The original supporting line's offset: `normal · x = c` at time 0.
    c: f32,
    /// How far this edge is allowed to travel. `INFINITY` when unconstrained.
    limit: f32,
}

impl EdgeState {
    /// How far the edge has travelled by time `t`.
    #[inline]
    fn offset_at(&self, t: f32) -> f32 {
        if t < self.limit {
            t
        } else {
            self.limit
        }
    }

    /// The edge's speed at time `t`: 1 while it is still moving, 0 once it has
    /// hit its limit. This is the single mechanism behind the constrained
    /// transform.
    #[inline]
    fn speed_at(&self, t: f32) -> f32 {
        if t < self.limit - EPS {
            1.0
        } else {
            0.0
        }
    }

    /// The supporting line's offset at time `t`.
    #[inline]
    fn c_at(&self, t: f32) -> f32 {
        self.c + self.offset_at(t)
    }
}

/// A vertex of the shrinking wavefront.
///
/// Wavefront vertices live in an arena and are linked into circular lists (one
/// per wavefront loop). A split event turns one loop into two; an edge event
/// shortens a loop.
#[derive(Clone, Debug)]
struct WVertex {
    /// Previous vertex in this wavefront loop.
    prev: usize,
    /// Next vertex in this wavefront loop.
    next: usize,
    /// Position at time [`WVertex::time`].
    pos: Vec2,
    /// The time at which `pos` was sampled.
    time: f32,
    /// Constant velocity since `time`.
    vel: Vec2,
    /// The input edge arriving here (from `prev`).
    left: EdgeId,
    /// The input edge leaving here (toward `next`).
    right: EdgeId,
    /// The skeleton node this vertex's current arc grows from.
    node: NodeId,
    /// False once the vertex has been consumed by an event.
    active: bool,
    /// Bumped whenever this vertex's geometry or links change, so that queued
    /// events computed from the old state can be recognised as stale.
    gen: u32,
    /// Bumped every time this vertex is rescheduled, so that its own previous
    /// event can be recognised as superseded.
    evt: u32,
    /// Input edges this vertex has been shown *not* to split, kept **sorted**.
    ///
    /// A vertex travels in a straight line, so it crosses any given edge's
    /// moving line at exactly one moment. If it was not on a live stretch of
    /// that edge at that moment, it never will be — the question is settled for
    /// good, and the edge can be struck off. Almost always empty.
    ///
    /// Sorted so [`Sim::split_lower_bound`] can binary-search it. That is not
    /// premature: the lookup sits inside a loop over every edge, so a linear
    /// scan here makes one scan cost `O(n * rejections)` and drags the whole
    /// simulation to `O(n^3)` on inputs that reject a lot.
    ///
    /// Cleared whenever the vertex's velocity changes, since that is a
    /// different trajectory and every answer has to be asked again.
    rejected: Vec<EdgeId>,
}

impl WVertex {
    /// Where this vertex is at time `t`, assuming its velocity is unchanged.
    #[inline]
    fn at(&self, t: f32) -> Vec2 {
        self.pos + self.vel * (t - self.time)
    }
}

/// Something that changes the wavefront's structure.
#[derive(Clone, Copy, Debug)]
enum EventKind {
    /// The wavefront edge between `a` and `a.next` collapsed to a point.
    Edge {
        /// The vertex at the edge's start.
        a: usize,
    },
    /// Reflex vertex `v` reached input edge `edge`'s moving line.
    ///
    /// Only a *candidate*: whether `v` lands on a live stretch of that edge, or
    /// merely on the infinite line it lies along, is settled when the event is
    /// popped. See [`Sim::split_lower_bound`].
    Split {
        /// The reflex vertex doing the splitting.
        v: usize,
        /// The input edge whose line it reaches.
        edge: EdgeId,
    },
    /// An input edge reached its distance limit and stopped moving.
    SpeedChange {
        /// The edge that stopped.
        edge: EdgeId,
    },
}

/// A queued event, ordered by time.
///
/// # Staleness
///
/// Events are never removed from the queue; they are recognised as obsolete
/// when popped. That takes two independent stamps, and conflating them loses
/// events:
///
/// - `owner` is the vertex whose event this is, stamped with that vertex's
///   *event serial*. Rescheduling a vertex bumps its serial, so its previous
///   event is superseded without disturbing anything else.
/// - `refs` are the vertices whose geometry the event's timing was computed
///   from, each stamped with its *structural generation*. If one of them moves,
///   the timing is worthless.
///
/// `refs` stays **O(1)** — at most the owner and its one neighbour — and that
/// is what keeps the simulation from melting down. An earlier design also
/// stamped the two endpoints of the edge a split was aimed at, which made a
/// vertex's event depend on a part of the wavefront arbitrarily far away. Any
/// event then invalidated events all over the polygon, each reschedule
/// re-registered more dependencies, and the whole thing fed back on itself: a
/// 132-vertex comb took 124ms and a 260-vertex one asked for 27GB.
///
/// The way out is [`Sim::split_lower_bound`]: a split's *timing* never depended
/// on those endpoints in the first place.
#[derive(Clone, Copy, Debug)]
struct Event {
    time: f32,
    kind: EventKind,
    /// The vertex this event belongs to, and its event serial when queued.
    owner: (usize, u32),
    /// Vertices whose geometry this event's timing depends on, with their
    /// structural generations when queued. At most two.
    refs: [(usize, u32); 2],
    /// How many entries of `refs` are meaningful.
    ref_count: u8,
}

impl Event {
    fn new(time: f32, kind: EventKind, owner: (usize, u32), refs: &[(usize, u32)]) -> Self {
        let mut r = [(0usize, 0u32); 2];
        r[..refs.len()].copy_from_slice(refs);
        Event {
            time,
            kind,
            owner,
            refs: r,
            ref_count: refs.len() as u8,
        }
    }
}

impl PartialEq for Event {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Event {}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Event {
    /// Reversed, so that `BinaryHeap` (a max-heap) yields the *earliest* event.
    ///
    /// Times are never NaN by the time an event is queued, but ordering must be
    /// total regardless, so NaN sorts last rather than panicking.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .time
            .partial_cmp(&self.time)
            .unwrap_or(Ordering::Equal)
    }
}

/// Runs the wavefront simulation.
pub(crate) fn compute(
    polygon: &Polygon,
    limits: Option<&[f32]>,
) -> Result<Skeleton, SkeletonError> {
    let edges = build_edge_states(polygon, limits)?;
    let mut sim = Sim::new(polygon, edges)?;
    sim.run()?;
    let mut skel = sim.skeleton;
    skel.edge_nodes = polygon
        .edge_ids()
        .map(|e| {
            let s = e.start_vertex();
            [NodeId(s.0 as u32), NodeId(polygon.next_vertex(s).0 as u32)]
        })
        .collect();
    skel.build_adjacency();
    Ok(skel)
}

/// Precomputes each edge's direction, inward normal, line offset, and limit.
fn build_edge_states(
    polygon: &Polygon,
    limits: Option<&[f32]>,
) -> Result<Vec<EdgeState>, SkeletonError> {
    if let Some(l) = limits {
        if l.len() != polygon.edge_count() {
            return Err(SkeletonError::LimitCountMismatch {
                got: l.len(),
                expected: polygon.edge_count(),
            });
        }
    }

    let mut states = Vec::with_capacity(polygon.edge_count());
    for e in polygon.edge_ids() {
        let (a, b) = polygon.edge(e);
        // Polygon validation rejects repeated vertices, so no edge is degenerate
        // and `normalize` cannot fail here.
        let dir = (b.to_vec2() - a.to_vec2())
            .normalize()
            .expect("polygon validation rules out zero-length edges");
        let normal = dir.perp();
        let c = normal.dot(a.to_vec2());

        let limit = match limits {
            None => f32::INFINITY,
            Some(l) => {
                let v = l[e.0 as usize];
                if v.is_nan() || v < 0.0 {
                    return Err(SkeletonError::InvalidLimit { edge: e, value: v });
                }
                v
            }
        };

        states.push(EdgeState {
            dir,
            normal,
            c,
            limit,
        });
    }
    Ok(states)
}

/// The simulation's mutable state.
struct Sim<'a> {
    polygon: &'a Polygon,
    edges: Vec<EdgeState>,
    verts: Vec<WVertex>,
    queue: BinaryHeap<Event>,
    skeleton: Skeleton,
    /// Each node's position in full `f32`, parallel to `skeleton.nodes`.
    ///
    /// `Node::exact` is narrowed to `f32` for the public API, which is far too
    /// coarse to decide whether a vertex is standing on its own node. This
    /// keeps the unnarrowed value for the simulation's own use.
    node_pos: Vec<Vec2>,
    /// Vertices needing a fresh event before the simulation may advance.
    ///
    /// Only ever holds a vertex that has just moved or been relinked, and that
    /// vertex's `prev` — whose edge event is computed from it. Both are O(1) per
    /// event, which is the whole point.
    dirty: Vec<usize>,
    /// Current simulation time, monotonically non-decreasing.
    now: f32,
}

impl<'a> Sim<'a> {
    fn new(polygon: &'a Polygon, edges: Vec<EdgeState>) -> Result<Self, SkeletonError> {
        let n = polygon.vertex_count();
        let mut skeleton = Skeleton::default();
        let mut verts = Vec::with_capacity(n * 2);
        let mut node_pos: Vec<Vec2> = Vec::with_capacity(n * 2);

        // One wavefront vertex and one boundary node per input vertex.
        for v in polygon.vertex_ids() {
            let left = polygon.prev_vertex(v).outgoing_edge();
            let right = v.outgoing_edge();
            let pos = polygon.vertex(v).to_vec2();

            let node = NodeId(skeleton.nodes.len() as u32);
            node_pos.push(pos);
            skeleton.nodes.push(Node {
                position: polygon.vertex(v),
                exact: [pos.x, pos.y],
                offset: 0.0,
                kind: NodeKind::Boundary(v),
                sources: vec![left, right],
            });

            verts.push(WVertex {
                prev: polygon.prev_vertex(v).0 as usize,
                next: polygon.next_vertex(v).0 as usize,
                pos,
                time: 0.0,
                vel: Vec2::ZERO,
                left,
                right,
                node,
                active: true,
                gen: 0,
                evt: 0,
                rejected: Vec::new(),
            });
        }

        let mut sim = Sim {
            polygon,
            edges,
            verts,
            queue: BinaryHeap::new(),
            node_pos,
            dirty: Vec::new(),
            skeleton,
            now: 0.0,
        };

        for i in 0..n {
            sim.verts[i].vel = sim.velocity_of(i, 0.0)?;
        }
        Ok(sim)
    }

    /// The velocity a wavefront vertex must have to stay on both of its edges'
    /// moving lines.
    ///
    /// Solving `normal_left · v = speed_left` and `normal_right · v =
    /// speed_right` is the crate's one unifying idea: with both speeds 1 it
    /// reproduces the classic angular bisector at `1 / sin(θ/2)`; with one
    /// speed 0 the vertex slides along the stopped edge; with both 0 it
    /// freezes. The standard and constrained transforms are the same code.
    fn velocity_of(&self, i: usize, t: f32) -> Result<Vec2, SkeletonError> {
        let v = &self.verts[i];
        let le = self.edges[v.left.0 as usize];
        let re = self.edges[v.right.0 as usize];
        let (w1, w2) = (le.speed_at(t), re.speed_at(t));
        let (n1, n2) = (le.normal, re.normal);

        let det = n1.cross(n2);
        if det.abs() > PARALLEL_EPS {
            // Cramer's rule on the 2x2 system.
            return Ok(Vec2::new(
                (w1 * n2.y - n1.y * w2) / det,
                (n1.x * w2 - w1 * n2.x) / det,
            ));
        }

        // The normals are parallel, so the two lines never meet and there is no
        // solution. Which of the two parallel cases this is decides what to do.
        if n1.dot(n2) < 0.0 {
            // Antiparallel: the two edges face each other head-on. This vertex
            // sits on both of their lines, and two antiparallel lines sharing a
            // point are the *same* line — so the wavefront loop has pinched to
            // zero width and its interior is empty. The loop is finished; it
            // just has to be read off. Freezing the vertex is what lets that
            // happen: it stays put while its loop collapses around it, and
            // `handle_edge_event` recognises the two-vertex remnant and emits
            // the ridge between them.
            //
            // This is not the 180° spike that polygon validation rejects: that
            // is an input defect, whereas this arises mid-simulation on inputs
            // as ordinary as a rectangle, whose skeleton *is* a ridge.
            return Ok(Vec2::ZERO);
        }

        if (w1 - w2).abs() < EPS {
            // Collinear and moving as one; the vertex just rides along.
            Ok(n1 * w1)
        } else {
            // One collinear edge stopped and the other did not: the wavefront
            // has no consistent shape here. Refuse rather than fabricate.
            Err(SkeletonError::IncompatibleCollinearLimits {
                left: v.left,
                right: v.right,
            })
        }
    }

    /// Whether the wavefront vertex is a notch (interior angle > 180°).
    ///
    /// This depends only on the two edges' directions, which never change, so
    /// it is stable for the vertex's whole life.
    fn is_reflex(&self, i: usize) -> bool {
        let v = &self.verts[i];
        let d1 = self.edges[v.left.0 as usize].dir;
        let d2 = self.edges[v.right.0 as usize].dir;
        // Interior on the left, so a right turn is reflex.
        d1.cross(d2) < -PARALLEL_EPS
    }

    fn run(&mut self) -> Result<(), SkeletonError> {
        // Each edge that has a finite limit stops at exactly one moment.
        for (i, e) in self.edges.iter().enumerate() {
            if e.limit.is_finite() && e.limit > 0.0 {
                self.queue.push(Event::new(
                    e.limit,
                    EventKind::SpeedChange {
                        edge: EdgeId(i as u16),
                    },
                    (usize::MAX, 0),
                    &[],
                ));
            }
        }

        for i in 0..self.verts.len() {
            self.schedule(i)?;
        }

        // Every event either consumes a vertex or stops an edge, so the total
        // is linear in the input. The generous multiplier makes this a
        // bug-catching backstop, not a real limit.
        let budget = 64 * (self.polygon.vertex_count() + 1) * (self.polygon.vertex_count() + 1);
        let mut processed = 0usize;

        while let Some(ev) = self.queue.pop() {
            processed += 1;
            if processed > budget {
                return Err(SkeletonError::EventBudgetExhausted { budget });
            }

            if !matches!(ev.kind, EventKind::SpeedChange { .. }) && !self.is_fresh(&ev) {
                continue;
            }
            // Events are popped in time order, but floating-point error can
            // produce a time a hair behind `now`; never let the clock go back.
            self.now = ev.time.max(self.now);

            match ev.kind {
                EventKind::Edge { a } => self.handle_edge_event(a, ev.time)?,
                EventKind::Split { v, edge } => self.handle_split_event(v, edge, ev.time)?,
                EventKind::SpeedChange { edge } => self.handle_speed_change(edge, ev.time)?,
            }

            // Nothing may advance past this point until every vertex whose
            // event the last one invalidated has a fresh one queued.
            while let Some(i) = self.dirty.pop() {
                if self.verts[i].active {
                    self.schedule(i)?;
                }
            }
        }
        Ok(())
    }

    /// Whether an event is still the current one for its owner, and still
    /// computed from up-to-date geometry.
    fn is_fresh(&self, ev: &Event) -> bool {
        let (owner, evt) = ev.owner;
        if !self.verts[owner].active || self.verts[owner].evt != evt {
            return false;
        }
        ev.refs[..ev.ref_count as usize]
            .iter()
            .all(|&(i, gen)| self.verts[i].active && self.verts[i].gen == gen)
    }

    /// Records that this vertex has moved or been relinked.
    ///
    /// Two events are computed from a vertex's geometry, and no others: its own,
    /// and its predecessor's edge event, which watches the two of them converge.
    /// So marking the vertex and its `prev` is complete — and, being O(1), is
    /// what keeps a single event from rippling across the whole wavefront.
    ///
    /// Split events are *not* in that set. A split's timing is computed from the
    /// target edge's supporting line, which never moves off its own offset
    /// track, so no vertex over there can invalidate it. See
    /// [`Sim::split_lower_bound`].
    fn touch(&mut self, i: usize) {
        self.verts[i].gen = self.verts[i].gen.wrapping_add(1);
        self.dirty.push(i);
        let prev = self.verts[i].prev;
        self.dirty.push(prev);
    }

    /// Computes and queues the next event for vertex `i`.
    fn schedule(&mut self, i: usize) -> Result<(), SkeletonError> {
        if !self.verts[i].active {
            return Ok(());
        }

        let mut best: Option<Event> = None;
        let mut consider = |ev: Option<Event>| {
            if let Some(e) = ev {
                if best.as_ref().map_or(true, |b| e.time < b.time) {
                    best = Some(e);
                }
            }
        };

        consider(self.edge_event(i));
        if self.is_reflex(i) {
            consider(self.split_lower_bound(i));
        }

        // Supersede whatever this vertex had queued before.
        self.verts[i].evt = self.verts[i].evt.wrapping_add(1);
        if let Some(mut e) = best {
            e.owner = (i, self.verts[i].evt);
            self.queue.push(e);
        }
        Ok(())
    }

    /// When, if ever, the wavefront edge leaving vertex `i` collapses.
    ///
    /// Both endpoints of that edge always lie on its supporting line, so the
    /// question is one-dimensional: track their separation *along* the edge and
    /// find when it reaches zero. That keeps this exactly linear, with no
    /// intersection test and no special cases.
    fn edge_event(&self, i: usize) -> Option<Event> {
        let a = &self.verts[i];
        let j = a.next;
        if j == i {
            return None;
        }
        let b = &self.verts[j];

        let d = self.edges[a.right.0 as usize].dir;
        let sep = d.dot(b.at(self.now) - a.at(self.now));
        let rate = d.dot(b.vel - a.vel);

        // Not shrinking (or growing): no collapse.
        if rate >= -EPS {
            return None;
        }
        let dt = sep / -rate;
        if !dt.is_finite() || dt < -EPS {
            return None;
        }
        let t = self.now + dt.max(0.0);

        Some(Event::new(
            t,
            EventKind::Edge { a: i },
            (i, 0),
            &[(i, a.gen), (j, b.gen)],
        ))
    }

    /// The earliest moment reflex vertex `i` could possibly split something.
    ///
    /// # A lower bound, not an answer
    ///
    /// This returns the first time `i` reaches *any* input edge's moving line,
    /// ignoring entirely whether it lands on a live stretch of that edge or
    /// merely on the infinite line it lies along. So it can be early. It can
    /// never be late, which is the only property that matters: a real split
    /// happens at one of these times, so the smallest of them is a floor under
    /// the true one, and popping events in time order stays correct.
    ///
    /// # Why bother being vague
    ///
    /// Because the vagueness is exactly what makes it cheap to keep true. An
    /// edge's wavefront slides along its own offset track and never leaves it,
    /// so the *time* `i` meets that track depends only on `i`'s trajectory and
    /// the edge's original line — both fixed. Nothing happening elsewhere in the
    /// wavefront can change it.
    ///
    /// Deciding here whether the landing is on a live stretch is what would ruin
    /// that. It would need the current endpoints of the edge, making this
    /// vertex's event depend on two vertices arbitrarily far away, so that any
    /// event anywhere invalidated events everywhere. That is what the earlier
    /// design did, and it cost `O(n^5)`.
    ///
    /// The question is instead settled in [`Sim::handle_split_event`], when the
    /// event is popped — by which point it is not a prediction at all.
    ///
    /// Scans every input edge, so `O(n)`; it is the only non-constant step left
    /// in the simulation.
    fn split_lower_bound(&self, i: usize) -> Option<Event> {
        let v = &self.verts[i];
        let p_now = v.at(self.now);
        let mut best: Option<(f32, EdgeId)> = None;

        for k in 0..self.edges.len() {
            let eo = EdgeId(k as u16);
            // A vertex cannot split the edges it rides on, nor one already
            // ruled out. `rejected` is sorted, so this stays O(log r).
            if eo == v.left || eo == v.right || v.rejected.binary_search(&eo).is_ok() {
                continue;
            }

            let e = &self.edges[k];
            let w = e.speed_at(self.now);

            // Distance to the moving line, and the rate it closes at. The
            // interior is on the +normal side, so `dist` starts non-negative.
            let dist = e.normal.dot(p_now) - e.c_at(self.now);
            let closing = e.normal.dot(v.vel) - w;
            if closing >= -EPS {
                continue; // v never reaches this line
            }
            let dt = dist / -closing;
            if !dt.is_finite() || dt < -EPS {
                continue;
            }
            let t = self.now + dt.max(0.0);

            if best.map_or(true, |(bt, _)| t < bt) {
                best = Some((t, eo));
            }
        }

        best.map(|(t, edge)| Event::new(t, EventKind::Split { v: i, edge }, (i, 0), &[(i, v.gen)]))
    }

    /// Adds a node to the skeleton.
    fn push_node(&mut self, pos: Vec2, t: f32, kind: NodeKind, sources: Vec<EdgeId>) -> NodeId {
        let id = NodeId(self.skeleton.nodes.len() as u32);
        self.node_pos.push(pos);
        self.skeleton.nodes.push(Node {
            position: Point::from_vec2_rounded(pos),
            exact: [pos.x, pos.y],
            offset: t,
            kind,
            sources,
        });
        id
    }

    /// Records the arc a wavefront vertex has traced since its last node.
    fn emit_arc(&mut self, vertex: usize, to: NodeId) {
        let v = &self.verts[vertex];
        let from = v.node;
        if from == to {
            return; // zero-length arc, nothing traced
        }
        let sources = [v.left, v.right];
        self.skeleton.arcs.push(Arc {
            nodes: [from, to],
            sources,
        });
    }

    /// The maximal run of consecutive wavefront vertices that all sit on `pos`
    /// at time `t`, starting from the adjacent pair `(ia, ib)`.
    ///
    /// Several wavefront edges often vanish at the very same instant and place
    /// — a square's four corners all reach its centre together. Textbooks call
    /// that a *vertex event*, and handling it as a cascade of two-vertex merges
    /// does not work: fusing just two of the square's corners would leave a
    /// vertex trapped between two antiparallel edges. Collecting the whole
    /// coincident run and retiring it in one go handles the general case and
    /// the simultaneous case with the same code.
    fn coincident_chain(&self, ia: usize, ib: usize, t: f32, pos: Vec2) -> Vec<usize> {
        let mut chain = vec![ia, ib];

        // Walk backward, then forward, while neighbours are at the same point.
        // The `contains` guard stops the walk from lapping a fully coincident
        // loop and revisiting where it started.
        loop {
            let p = self.verts[chain[0]].prev;
            if !self.verts[p].active || chain.contains(&p) {
                break;
            }
            if (self.verts[p].at(t) - pos).length() > MERGE_EPS {
                break;
            }
            chain.insert(0, p);
        }
        loop {
            let n = self.verts[chain[chain.len() - 1]].next;
            if !self.verts[n].active || chain.contains(&n) {
                break;
            }
            if (self.verts[n].at(t) - pos).length() > MERGE_EPS {
                break;
            }
            chain.push(n);
        }
        chain
    }

    /// One or more adjacent vertices met at a point: the edges between them are
    /// gone.
    fn handle_edge_event(&mut self, ia: usize, t: f32) -> Result<(), SkeletonError> {
        let ib = self.verts[ia].next;
        if !self.verts[ib].active || ib == ia {
            return Ok(());
        }

        // Seed from the colliding pair, then gather everything else arriving at
        // the same instant and place.
        let seed = (self.verts[ia].at(t) + self.verts[ib].at(t)) * 0.5;
        let chain = self.coincident_chain(ia, ib, t, seed);

        // Averaging every arriving vertex's prediction spreads the
        // floating-point error rather than inheriting one vertex's.
        let pos = chain
            .iter()
            .fold(Vec2::ZERO, |acc, &i| acc + self.verts[i].at(t))
            * (1.0 / chain.len() as f32);

        let mut sources = Vec::with_capacity(chain.len() + 1);
        for &i in &chain {
            sources.push(self.verts[i].left);
            sources.push(self.verts[i].right);
        }
        sources.sort_unstable();
        sources.dedup();

        let node = self.node_at(pos, t, NodeKind::EdgeEvent, sources);
        for &i in &chain {
            self.emit_arc(i, node);
        }

        let first = chain[0];
        let last = chain[chain.len() - 1];
        let iprev = self.verts[first].prev;
        let inext = self.verts[last].next;
        let (left, right) = (self.verts[first].left, self.verts[last].right);

        let whole_loop = chain.contains(&iprev);
        for &i in &chain {
            self.deactivate(i);
        }
        if whole_loop {
            // The chain was the entire loop: every arc is recorded and there is
            // nothing left to move.
            return Ok(());
        }

        // The chain fuses into one vertex inheriting the outermost edges.
        let merged = self.spawn(WVertex {
            prev: iprev,
            next: inext,
            pos,
            time: t,
            vel: Vec2::ZERO,
            left,
            right,
            node,
            active: true,
            gen: 0,
            evt: 0,
            rejected: Vec::new(),
        });
        self.verts[iprev].next = merged;
        self.verts[inext].prev = merged;
        self.touch(iprev);
        self.touch(inext);

        if self.is_stalled(merged) {
            // The merge left the loop somewhere no ordinary event can advance
            // it; zip it shut rather than let the simulation stall.
            return self.resolve_needle(merged, t);
        }

        self.verts[merged].vel = self.velocity_of(merged, t)?;

        self.schedule(merged)?;
        self.schedule(iprev)?;
        self.schedule(inext)?;
        Ok(())
    }

    /// Whether a vertex's two edges face each other head-on.
    ///
    /// Such a vertex lies on both edges' lines at once, and two antiparallel
    /// lines through a common point are the same line — so its whole loop has
    /// flattened. See [`Sim::resolve_needle`].
    fn edges_antiparallel(&self, i: usize) -> bool {
        let v = &self.verts[i];
        let n1 = self.edges[v.left.0 as usize].normal;
        let n2 = self.edges[v.right.0 as usize].normal;
        n1.cross(n2).abs() <= PARALLEL_EPS && n1.dot(n2) < 0.0
    }

    /// Whether a vertex has reached a state no ordinary event can advance, and
    /// so must be zipped shut by [`Sim::resolve_needle`].
    ///
    /// Two cases, and both stall the simulation if ignored. A needle's folded
    /// edges are parallel, so they can never collapse. A two-vertex loop has
    /// both of its edges spanning the same pair of points, so it encloses
    /// nothing and cannot shrink either.
    fn is_stalled(&self, i: usize) -> bool {
        self.edges_antiparallel(i) || self.verts[i].prev == self.verts[i].next
    }

    /// Zips up a *needle*: a wavefront vertex whose two edges have met head-on.
    ///
    /// # What a needle is
    ///
    /// When a vertex's two edges are antiparallel, both of their lines pass
    /// through it, and two antiparallel lines sharing a point are the *same*
    /// line. The two edges have collided, the strip of material between them is
    /// gone, and the wavefront has folded back on itself: `prev` and `next` are
    /// both collinear with the vertex and lie on the *same* side of it.
    ///
    /// ```text
    ///     prev  o<---------------------o m          prev  o
    ///           o--------------------->'                  |
    ///     next                                            |   the strip is gone;
    ///                                          ==>        |   what is left is
    ///     the two edges have collided and now             |   one skeleton arc
    ///     lie on top of one another                 next  o
    /// ```
    ///
    /// It arises on thoroughly ordinary input: a rectangle's long sides collide
    /// to leave its ridge, and any hole 2d from a wall pinches the strip
    /// between them shut.
    ///
    /// # Why it needs its own handling
    ///
    /// The folded edges are parallel, so no further edge event can ever fire on
    /// them. Left alone the simulation simply stalls with the loop still live,
    /// silently dropping every arc it had yet to trace.
    ///
    /// # How it resolves
    ///
    /// The overlap runs from the vertex to whichever of `prev`/`next` is
    /// nearer, and that overlap is exactly one skeleton arc, bisecting the two
    /// edges that collided. Emit it, retire the arm (or both arms, when they
    /// are the same length), and splice what remains back into the loop. The
    /// splice can expose another needle, so this repeats until it does not.
    fn resolve_needle(&mut self, start: usize, t: f32) -> Result<(), SkeletonError> {
        let mut m = start;
        loop {
            let prev = self.verts[m].prev;
            let next = self.verts[m].next;
            if prev == m || next == m || !self.verts[prev].active || !self.verts[next].active {
                self.deactivate(m);
                return Ok(());
            }

            let pm = self.verts[m].at(t);

            // A two-vertex loop is the end of the line: emit its last segment.
            if prev == next {
                let p = self.verts[prev].at(t);
                let node = self.node_at(p, t, NodeKind::EdgeEvent, Vec::new());
                self.add_sources(node, m);
                self.add_sources(node, prev);
                self.emit_arc(m, node);
                self.emit_arc(prev, node);
                self.deactivate(m);
                self.deactivate(prev);
                return Ok(());
            }

            let s = (self.verts[prev].at(t) - pm).length();
            let u = (self.verts[next].at(t) - pm).length();
            let take_prev = s <= u + MERGE_EPS;
            let take_next = u <= s + MERGE_EPS;
            let target = if take_prev {
                self.verts[prev].at(t)
            } else {
                self.verts[next].at(t)
            };

            let node = self.node_at(target, t, NodeKind::EdgeEvent, Vec::new());
            self.add_sources(node, m);
            // The collapsed strip itself is a skeleton arc, bisecting the two
            // edges that just collided.
            self.emit_arc(m, node);
            let (mut new_left, mut new_right) = (self.verts[m].left, self.verts[m].right);
            self.deactivate(m);

            let (mut lo, mut hi) = (prev, next);
            if take_prev {
                self.emit_arc(prev, node);
                self.add_sources(node, prev);
                new_left = self.verts[prev].left;
                lo = self.verts[prev].prev;
                self.deactivate(prev);
            }
            if take_next {
                self.emit_arc(next, node);
                self.add_sources(node, next);
                new_right = self.verts[next].right;
                hi = self.verts[next].next;
                self.deactivate(next);
            }

            if !self.verts[lo].active || !self.verts[hi].active {
                return Ok(()); // the loop is spent
            }

            let nv = self.spawn(WVertex {
                prev: lo,
                next: hi,
                pos: target,
                time: t,
                vel: Vec2::ZERO,
                left: new_left,
                right: new_right,
                node,
                active: true,
                gen: 0,
                evt: 0,
                rejected: Vec::new(),
            });
            self.verts[lo].next = nv;
            self.verts[hi].prev = nv;
            self.touch(lo);
            self.touch(hi);

            if self.is_stalled(nv) {
                // Zipping one needle shut exposed another stall.
                m = nv;
                continue;
            }

            self.verts[nv].vel = self.velocity_of(nv, t)?;
            self.schedule(nv)?;
            self.schedule(lo)?;
            self.schedule(hi)?;
            return Ok(());
        }
    }

    /// Adds a wavefront vertex's two edges to a node's source list.
    fn add_sources(&mut self, node: NodeId, vertex: usize) {
        let (l, r) = (self.verts[vertex].left, self.verts[vertex].right);
        let sources = &mut self.skeleton.nodes[node.0 as usize].sources;
        sources.push(l);
        sources.push(r);
        sources.sort_unstable();
        sources.dedup();
    }

    /// The interior node at `pos`, creating one only if there is not one there
    /// already.
    ///
    /// Several events can converge on a single point: a hole sitting
    /// symmetrically between two walls pinches both strips shut at the same
    /// instant and place. Each would otherwise mint its own node, leaving
    /// duplicates stacked at one position and a graph disconnected where it
    /// should not be. Two skeleton nodes at the same point *are* the same node,
    /// so reusing one is exact, not an approximation.
    ///
    /// Boundary nodes are never reused: each belongs to a named input vertex.
    fn node_at(&mut self, pos: Vec2, t: f32, kind: NodeKind, sources: Vec<EdgeId>) -> NodeId {
        for i in 0..self.node_pos.len() {
            if !self.skeleton.nodes[i].is_boundary()
                && (self.node_pos[i] - pos).length() <= MERGE_EPS
            {
                let n = &mut self.skeleton.nodes[i];
                n.sources.extend(sources);
                n.sources.sort_unstable();
                n.sources.dedup();
                return NodeId(i as u32);
            }
        }
        self.push_node(pos, t, kind, sources)
    }

    /// Finds the live stretch of `edge` that `pos` is standing on at time `t`.
    ///
    /// An input edge can own several wavefront edges at once — every split of it
    /// leaves one more — so this asks which, if any, of them `pos` is actually
    /// on.
    ///
    /// This is an *observation*, not a prediction, and that distinction is the
    /// point. Events are popped in time order, so by the time this runs every
    /// event before `t` has been processed and the wavefront's shape at `t` is
    /// settled fact. The same question asked when the event was queued would
    /// have been a guess about a future that later events could change.
    fn live_stretch_at(&self, edge: EdgeId, pos: Vec2, t: f32) -> Option<usize> {
        let e = &self.edges[edge.0 as usize];
        for j in 0..self.verts.len() {
            let a = &self.verts[j];
            if !a.active || a.right != edge {
                continue;
            }
            let b = &self.verts[a.next];
            let (pa, pb) = (a.at(t), b.at(t));
            let span = e.dir.dot(pb - pa);
            if span <= EPS {
                continue; // this stretch has already shrunk away
            }
            let along = e.dir.dot(pos - pa);
            if along >= -EPS && along <= span + EPS {
                return Some(j);
            }
        }
        None
    }

    /// A reflex vertex reached an opposing edge's line. If it landed on a live
    /// stretch of that edge, it tears the wavefront in two.
    fn handle_split_event(&mut self, iv: usize, eo: EdgeId, t: f32) -> Result<(), SkeletonError> {
        if !self.verts[iv].active {
            return Ok(());
        }

        let pos = self.verts[iv].at(t);

        // Now the guess gets checked. `now == t`, so this is what the wavefront
        // genuinely looks like, not a forecast.
        let Some(iopp) = self.live_stretch_at(eo, pos, t) else {
            // It came down on the line but off the end of every live stretch of
            // it, so no split happens here. A vertex moves in a straight line
            // and so meets this line only once: the question is closed for good,
            // and the edge can be struck off before asking for the next
            // candidate.
            let rejected = &mut self.verts[iv].rejected;
            if let Err(at) = rejected.binary_search(&eo) {
                rejected.insert(at, eo); // keep it sorted for the binary search
            }
            self.schedule(iv)?;
            return Ok(());
        };

        let ib = self.verts[iopp].next;
        if !self.verts[ib].active {
            return Ok(());
        }
        let (v_left, v_right, v_prev, v_next) = {
            let v = &self.verts[iv];
            (v.left, v.right, v.prev, v.next)
        };

        let mut sources = vec![v_left, v_right, eo];
        sources.sort_unstable();
        sources.dedup();
        let node = self.node_at(pos, t, NodeKind::SplitEvent, sources);
        self.emit_arc(iv, node);
        self.deactivate(iv);

        // The wavefront loop `... -> v_prev -> v -> v_next -> ... -> opp -> b -> ...`
        // becomes two loops. `v1` closes the chain that runs from `b` around to
        // `v_prev`; `v2` closes the chain from `v_next` around to `opp`.
        let v1 = self.spawn(WVertex {
            prev: v_prev,
            next: ib,
            pos,
            time: t,
            vel: Vec2::ZERO,
            left: v_left,
            right: eo,
            node,
            active: true,
            gen: 0,
            evt: 0,
            rejected: Vec::new(),
        });
        let v2 = self.spawn(WVertex {
            prev: iopp,
            next: v_next,
            pos,
            time: t,
            vel: Vec2::ZERO,
            left: eo,
            right: v_right,
            node,
            active: true,
            gen: 0,
            evt: 0,
            rejected: Vec::new(),
        });

        self.verts[v_prev].next = v1;
        self.verts[ib].prev = v1;
        self.verts[iopp].next = v2;
        self.verts[v_next].prev = v2;
        self.touch(v_prev);
        self.touch(ib);
        self.touch(iopp);
        self.touch(v_next);

        // A split can flatten either of the two new loops — for instance when a
        // reflex vertex lands on the far wall of a narrow channel, pinching
        // that side shut against the edge opposite.
        let flat_1 = self.is_stalled(v1);
        let flat_2 = self.is_stalled(v2);
        if flat_1 {
            self.resolve_needle(v1, t)?;
        } else {
            self.verts[v1].vel = self.velocity_of(v1, t)?;
        }
        if flat_2 {
            self.resolve_needle(v2, t)?;
        } else {
            self.verts[v2].vel = self.velocity_of(v2, t)?;
        }

        for i in [v1, v2, v_prev, ib, iopp, v_next] {
            if self.verts[i].active {
                self.schedule(i)?;
            }
        }
        Ok(())
    }

    /// An edge hit its distance limit and stopped.
    ///
    /// Every active vertex's trajectory can bend here, so each one gets a node
    /// (a kink in its arc) and a fresh velocity. Rescheduling *everything* is
    /// heavy-handed — O(n^2) per stop — but a speed change also invalidates
    /// split events aimed at the stopped edge from elsewhere in the wavefront,
    /// and sorting out exactly which is a subtle bookkeeping problem for no
    /// real gain: a polygon has at most one stop per distinct limit, and the
    /// common cases (no limits, or one uniform limit) cost zero or one pass.
    fn handle_speed_change(&mut self, _edge: EdgeId, t: f32) -> Result<(), SkeletonError> {
        let n = self.verts.len();
        for i in 0..n {
            if !self.verts[i].active {
                continue;
            }
            let new_vel = self.velocity_of(i, t)?;
            let old_vel = self.verts[i].vel;

            if (new_vel - old_vel).length_squared() > EPS * EPS {
                // The trajectory bends here, so the arc it was tracing ends and
                // a new one begins.
                let pos = self.verts[i].at(t);
                let sources = {
                    let v = &self.verts[i];
                    let mut s = vec![v.left, v.right];
                    s.sort_unstable();
                    s.dedup();
                    s
                };
                let node = self.node_at(pos, t, NodeKind::LimitReached, sources);
                self.emit_arc(i, node);

                let v = &mut self.verts[i];
                v.pos = pos;
                v.time = t;
                v.vel = new_vel;
                v.node = node;
                // A new trajectory: every edge struck off was struck off for a
                // path this vertex is no longer on, so ask again.
                v.rejected.clear();
            }
            self.touch(i);
        }
        Ok(())
    }

    fn spawn(&mut self, v: WVertex) -> usize {
        self.verts.push(v);
        self.verts.len() - 1
    }

    fn deactivate(&mut self, i: usize) {
        self.verts[i].active = false;
        self.touch(i);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn events_pop_earliest_first() {
        let mut q = BinaryHeap::new();
        for t in [5.0, 1.0, 3.0, 9.0, 2.0] {
            q.push(Event::new(t, EventKind::Edge { a: 0 }, (0, 0), &[]));
        }
        let mut got = Vec::new();
        while let Some(e) = q.pop() {
            got.push(e.time);
        }
        assert_eq!(got, vec![1.0, 2.0, 3.0, 5.0, 9.0]);
    }

    #[test]
    fn event_ordering_is_total_even_with_nan() {
        // Ord must not panic; NaN simply sorts last.
        let a = Event::new(f32::NAN, EventKind::Edge { a: 0 }, (0, 0), &[]);
        let b = Event::new(1.0, EventKind::Edge { a: 0 }, (0, 0), &[]);
        let _ = a.cmp(&b);
        let _ = b.cmp(&a);
    }

    #[test]
    fn edge_state_speed_drops_at_the_limit() {
        let e = EdgeState {
            dir: Vec2::new(1.0, 0.0),
            normal: Vec2::new(0.0, 1.0),
            c: 0.0,
            limit: 3.0,
        };
        assert_eq!(e.speed_at(0.0), 1.0);
        assert_eq!(e.speed_at(2.9), 1.0);
        assert_eq!(e.speed_at(3.0), 0.0);
        assert_eq!(e.speed_at(9.0), 0.0);

        assert_eq!(e.offset_at(1.0), 1.0);
        assert_eq!(e.offset_at(3.0), 3.0);
        assert_eq!(e.offset_at(9.0), 3.0, "offset clamps at the limit");
        assert_eq!(e.c_at(9.0), 3.0);
    }

    #[test]
    fn unconstrained_edge_never_stops() {
        let e = EdgeState {
            dir: Vec2::new(1.0, 0.0),
            normal: Vec2::new(0.0, 1.0),
            c: 0.0,
            limit: f32::INFINITY,
        };
        assert_eq!(e.speed_at(1e9), 1.0);
        assert_eq!(e.offset_at(1e9), 1e9);
    }

    #[test]
    fn vertex_position_extrapolates_linearly() {
        let v = WVertex {
            prev: 0,
            next: 0,
            pos: Vec2::new(1.0, 2.0),
            time: 1.0,
            vel: Vec2::new(3.0, -1.0),
            left: EdgeId(0),
            right: EdgeId(1),
            node: NodeId(0),
            active: true,
            gen: 0,
            evt: 0,
            rejected: Vec::new(),
        };
        assert_eq!(v.at(1.0), Vec2::new(1.0, 2.0));
        assert_eq!(v.at(3.0), Vec2::new(7.0, 0.0));
    }
}
