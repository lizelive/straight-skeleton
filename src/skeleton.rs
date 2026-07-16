//! The computed straight skeleton and its provenance queries.

use alloc::vec;
use alloc::vec::Vec;

use crate::polygon::{EdgeId, VertexId};
use crate::Point;

/// Identifies a [`Node`] of a [`Skeleton`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct NodeId(pub u32);

/// Identifies an [`Arc`] of a [`Skeleton`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ArcId(pub u32);

/// What produced a [`Node`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum NodeKind {
    /// A node sitting on the input boundary, at offset 0. There is exactly one
    /// of these per input vertex, and it carries that vertex's id.
    Boundary(VertexId),
    /// An interior node, created when the shrinking wavefront collapsed an edge
    /// to nothing and two skeleton arcs met.
    EdgeEvent,
    /// An interior node, created when a reflex vertex ran into an opposing
    /// wavefront edge and tore the wavefront into two loops.
    SplitEvent,
    /// A node where the wavefront stopped because every incident edge had
    /// reached its per-edge distance limit.
    ///
    /// Only produced by [`skeleton_constrained`]; a plain skeleton never has
    /// these.
    ///
    /// [`skeleton_constrained`]: crate::skeleton_constrained
    LimitReached,
}

/// A vertex of the skeleton graph.
///
/// # Examples
///
/// ```
/// use straight_skeleton::{skeleton, NodeKind, Point, Polygon};
///
/// let square = Polygon::from_outer(&[
///     Point::new(0, 0), Point::new(10, 0), Point::new(10, 10), Point::new(0, 10),
/// ])?;
/// let skel = skeleton(&square)?;
///
/// // Every input vertex gets a boundary node at offset 0.
/// let boundary: Vec<_> = skel.nodes().iter().filter(|n| n.is_boundary()).collect();
/// assert_eq!(boundary.len(), 4);
/// assert!(boundary.iter().all(|n| n.offset == 0.0));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Node {
    /// The node's position, rounded to the integer lattice.
    ///
    /// Skeleton nodes are generally irrational even for integer input, so this
    /// is the nearest lattice point. Use [`Node::exact`] when the rounding
    /// matters.
    pub position: Point,
    /// The node's unrounded position.
    ///
    /// The algorithm computes in `f32` throughout, so this is the value it
    /// actually arrived at, not a narrowing of something wider.
    pub exact: [f32; 2],
    /// How far the wavefront had travelled when this node was created.
    ///
    /// For a plain [`skeleton`], this is the node's distance to the supporting
    /// line of each of its [`Node::sources`]. It is the node's height on a
    /// roof, and the offset at which the node appears on an offset curve.
    ///
    /// For a [`skeleton_constrained`], it is the wavefront's **time**, which is
    /// no longer the same thing: an edge that stopped at `limit` stays `limit`
    /// away however long the wavefront runs on. The distance to a source edge
    /// `e`'s line is `min(offset, limit_e)`.
    ///
    /// [`skeleton`]: crate::skeleton
    /// [`skeleton_constrained`]: crate::skeleton_constrained
    pub offset: f32,
    /// What produced this node.
    pub kind: NodeKind,
    /// The input edges whose wavefronts arrived here together.
    ///
    /// Always at least 2 entries, and 3 or more where several skeleton arcs
    /// meet. Each one's supporting line is [`Node::offset`] away (see that
    /// field for the constrained case).
    ///
    /// # This is not quite "nearest"
    ///
    /// For a **convex** polygon these really are the nearest input edges, since
    /// there the straight skeleton coincides with the medial axis.
    ///
    /// Elsewhere they may not be, and the difference is the definition of a
    /// straight skeleton rather than a wrinkle in this implementation. A
    /// straight skeleton bisects edges' infinite **supporting lines**, which is
    /// what keeps every arc straight. A medial axis bisects the nearest
    /// **features**, and so grows parabolic arcs around reflex vertices. Around
    /// a reflex corner the two part company: a plus-shape's centre is at offset
    /// 5 from the four arms' walls, but the nearest input feature is a reflex
    /// corner 7.07 away.
    ///
    /// So read `sources` as *"the input edges whose faces meet here"* — which is
    /// the useful notion anyway, and the one a roof needs. See [`Arc::sources`].
    pub sources: Vec<EdgeId>,
}

impl Node {
    /// Whether this node lies on the input boundary.
    #[inline]
    pub fn is_boundary(&self) -> bool {
        matches!(self.kind, NodeKind::Boundary(_))
    }

    /// The input vertex this node sits on, if it is a boundary node.
    #[inline]
    pub fn input_vertex(&self) -> Option<VertexId> {
        match self.kind {
            NodeKind::Boundary(v) => Some(v),
            _ => None,
        }
    }
}

/// An edge of the skeleton graph: a straight segment traced by one wavefront
/// vertex as it moved.
///
/// # Provenance
///
/// Every arc separates the [faces] of **exactly two input edges** — the two in
/// [`Arc::sources`] — and every point along it is equidistant from those two
/// edges' supporting lines. This is not an approximation or a post-hoc lookup:
/// it is what an arc *is*. A wavefront vertex exists precisely where two
/// shrinking edges meet, and it carries those two edge ids with it as it sweeps
/// the arc out.
///
/// So tracing an arc back to the input it came from is a field access, at zero
/// cost. See [`Node::sources`] for the one caveat: on a non-convex polygon
/// these are the edges whose faces meet here, which is not always the same as
/// the Euclidean-nearest edges.
///
/// [faces]: Skeleton::face
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Arc {
    /// The arc's two endpoints. `nodes[0]` is always the lower-offset end, so
    /// arcs point "uphill", away from the boundary.
    pub nodes: [NodeId; 2],
    /// The two input edges whose faces this arc separates.
    pub sources: [EdgeId; 2],
}

impl Arc {
    /// The endpoint nearer the input boundary.
    #[inline]
    pub fn lower(&self) -> NodeId {
        self.nodes[0]
    }

    /// The endpoint further from the input boundary.
    #[inline]
    pub fn upper(&self) -> NodeId {
        self.nodes[1]
    }

    /// The endpoint that is not `n`, or `None` if `n` is not an endpoint.
    #[inline]
    pub fn other(&self, n: NodeId) -> Option<NodeId> {
        if self.nodes[0] == n {
            Some(self.nodes[1])
        } else if self.nodes[1] == n {
            Some(self.nodes[0])
        } else {
            None
        }
    }
}

/// The straight skeleton of a [`Polygon`].
///
/// A skeleton is a planar graph of [`Node`]s joined by [`Arc`]s. Build one with
/// [`skeleton`] or [`skeleton_constrained`].
///
/// [`Polygon`]: crate::Polygon
/// [`skeleton`]: crate::skeleton
/// [`skeleton_constrained`]: crate::skeleton_constrained
///
/// # Examples
///
/// ```
/// use straight_skeleton::{skeleton, Point, Polygon};
///
/// // A 10x10 square's skeleton is an X: four boundary nodes, one centre node,
/// // four arcs running corner to centre.
/// let square = Polygon::from_outer(&[
///     Point::new(0, 0), Point::new(10, 0), Point::new(10, 10), Point::new(0, 10),
/// ])?;
/// let skel = skeleton(&square)?;
///
/// assert_eq!(skel.node_count(), 5);
/// assert_eq!(skel.arc_count(), 4);
///
/// // The interior node is the centre, at offset 5 (half the width).
/// let centre = skel.nodes().iter().find(|n| !n.is_boundary()).unwrap();
/// assert_eq!(centre.position, Point::new(5, 5));
/// assert!((centre.offset - 5.0).abs() < 1e-4);
///
/// // ...and it is equidistant from all four input edges.
/// assert_eq!(centre.sources.len(), 4);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Skeleton {
    pub(crate) nodes: Vec<Node>,
    pub(crate) arcs: Vec<Arc>,
    /// `node_arcs[i]` lists the arcs incident to node `i`. Built once at the
    /// end of the algorithm so that traversal queries are O(degree).
    pub(crate) node_arcs: Vec<Vec<ArcId>>,
    /// `edge_nodes[i]` is the pair of boundary nodes at input edge `i`'s start
    /// and end vertices, which is where [`Skeleton::face`] begins its walk.
    pub(crate) edge_nodes: Vec<[NodeId; 2]>,
}

impl Skeleton {
    /// All nodes.
    #[inline]
    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    /// All arcs.
    #[inline]
    pub fn arcs(&self) -> &[Arc] {
        &self.arcs
    }

    /// Number of nodes.
    #[inline]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of arcs.
    #[inline]
    pub fn arc_count(&self) -> usize {
        self.arcs.len()
    }

    /// Looks up a node.
    ///
    /// # Panics
    ///
    /// Panics if `n` does not belong to this skeleton.
    #[inline]
    pub fn node(&self, n: NodeId) -> &Node {
        &self.nodes[n.0 as usize]
    }

    /// Looks up an arc.
    ///
    /// # Panics
    ///
    /// Panics if `a` does not belong to this skeleton.
    #[inline]
    pub fn arc(&self, a: ArcId) -> &Arc {
        &self.arcs[a.0 as usize]
    }

    /// Iterates node ids.
    pub fn node_ids(&self) -> impl Iterator<Item = NodeId> + '_ {
        (0..self.nodes.len() as u32).map(NodeId)
    }

    /// Iterates arc ids.
    pub fn arc_ids(&self) -> impl Iterator<Item = ArcId> + '_ {
        (0..self.arcs.len() as u32).map(ArcId)
    }

    /// The arcs incident to a node.
    ///
    /// # Panics
    ///
    /// Panics if `n` does not belong to this skeleton.
    #[inline]
    pub fn arcs_at(&self, n: NodeId) -> &[ArcId] {
        &self.node_arcs[n.0 as usize]
    }

    /// The arc's two endpoints, as positions.
    ///
    /// # Panics
    ///
    /// Panics if `a` does not belong to this skeleton.
    #[inline]
    pub fn arc_segment(&self, a: ArcId) -> (Point, Point) {
        let arc = self.arc(a);
        (
            self.node(arc.nodes[0]).position,
            self.node(arc.nodes[1]).position,
        )
    }

    /// The two input edges a given arc came from.
    ///
    /// Every point along the arc is equidistant from these two edges'
    /// supporting lines, and the arc separates their two faces. See [`Arc`] for
    /// why this is exact rather than a search, and [`Node::sources`] for why
    /// "came from" is more accurate than "closest to".
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::{skeleton, Point, Polygon};
    ///
    /// let square = Polygon::from_outer(&[
    ///     Point::new(0, 0), Point::new(10, 0), Point::new(10, 10), Point::new(0, 10),
    /// ])?;
    /// let skel = skeleton(&square)?;
    ///
    /// // Each of the square's four arcs bisects two adjacent input edges.
    /// for a in skel.arc_ids() {
    ///     let [e0, e1] = skel.closest_input_edges(a);
    ///     assert_ne!(e0, e1);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if `a` does not belong to this skeleton.
    #[inline]
    pub fn closest_input_edges(&self, a: ArcId) -> [EdgeId; 2] {
        self.arc(a).sources
    }

    /// The input edges a given node came from.
    ///
    /// See [`Node::sources`], which this returns.
    ///
    /// # Panics
    ///
    /// Panics if `n` does not belong to this skeleton.
    #[inline]
    pub fn closest_input_edges_to_node(&self, n: NodeId) -> &[EdgeId] {
        &self.node(n).sources
    }

    /// The largest offset reached by any node: the radius of the largest disc
    /// that fits inside the polygon.
    ///
    /// For a roof, this is the ridge height. Returns 0 for an empty skeleton.
    pub fn max_offset(&self) -> f32 {
        self.nodes.iter().map(|n| n.offset).fold(0.0, f32::max)
    }

    /// The boundary node sitting on a given input vertex.
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::{skeleton, Point, Polygon, VertexId};
    ///
    /// let square = Polygon::from_outer(&[
    ///     Point::new(0, 0), Point::new(10, 0), Point::new(10, 10), Point::new(0, 10),
    /// ])?;
    /// let skel = skeleton(&square)?;
    ///
    /// let n = skel.boundary_node(VertexId(2)).unwrap();
    /// assert_eq!(skel.node(n).position, Point::new(10, 10));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn boundary_node(&self, v: VertexId) -> Option<NodeId> {
        // The algorithm emits one boundary node per input vertex, in order,
        // before any interior node, so the ids line up.
        let id = NodeId(v.0 as u32);
        match self.nodes.get(v.0 as usize)?.kind {
            NodeKind::Boundary(w) if w == v => Some(id),
            _ => None,
        }
    }

    /// The **face** of an input edge: the closed region the wavefront of that
    /// edge swept out, as a loop of nodes.
    ///
    /// Every input edge has exactly one face, and the faces tile the polygon.
    /// The returned loop starts with the edge's own two endpoints, then follows
    /// the skeleton arcs that bound the face back around. Each face is planar
    /// when nodes are lifted to `z = offset`, which is exactly why a straight
    /// skeleton builds roofs: **one face is one roof plane**. See the `roof`
    /// example.
    ///
    /// Returns `None` if the face cannot be walked, which should not happen for
    /// a skeleton of a valid polygon.
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::{skeleton, EdgeId, Point, Polygon};
    ///
    /// // Each of a square's four edges has a triangular face running to the
    /// // centre.
    /// let square = Polygon::from_outer(&[
    ///     Point::new(0, 0), Point::new(10, 0), Point::new(10, 10), Point::new(0, 10),
    /// ])?;
    /// let skel = skeleton(&square)?;
    ///
    /// let face = skel.face(EdgeId(0)).unwrap();
    /// assert_eq!(face.len(), 3);
    ///
    /// // A rectangle's long edges get quadrilateral faces, because the ridge
    /// // gives them a fourth corner.
    /// let rect = Polygon::from_outer(&[
    ///     Point::new(0, 0), Point::new(20, 0), Point::new(20, 10), Point::new(0, 10),
    /// ])?;
    /// let skel = skeleton(&rect)?;
    /// assert_eq!(skel.face(EdgeId(0)).unwrap().len(), 4); // long edge
    /// assert_eq!(skel.face(EdgeId(1)).unwrap().len(), 3); // short edge
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn face(&self, e: EdgeId) -> Option<Vec<NodeId>> {
        let [start, end] = *self.edge_nodes.get(e.0 as usize)?;
        let mut loop_ = vec![start, end];

        // Walk from the edge's far endpoint back to its near one, taking only
        // arcs that this edge is a source of — those are precisely the arcs
        // bounding its face.
        let mut cur = end;
        let mut came_from: Option<ArcId> = None;
        loop {
            let next_arc = self
                .arcs_at(cur)
                .iter()
                .copied()
                .find(|&a| Some(a) != came_from && self.arc(a).sources.contains(&e))?;
            let other = self.arc(next_arc).other(cur)?;
            if other == start {
                return Some(loop_);
            }
            loop_.push(other);
            came_from = Some(next_arc);
            cur = other;
            // A face cannot have more corners than the skeleton has nodes.
            if loop_.len() > self.nodes.len() + 2 {
                return None;
            }
        }
    }

    /// How many input edges the polygon had.
    ///
    /// Each one owns exactly one [`face`](Skeleton::face), so this is also the
    /// number of faces.
    ///
    /// # Examples
    ///
    /// ```
    /// use straight_skeleton::{skeleton, Point, Polygon};
    ///
    /// let square = Polygon::from_outer(&[
    ///     Point::new(0, 0), Point::new(10, 0), Point::new(10, 10), Point::new(0, 10),
    /// ])?;
    /// assert_eq!(skeleton(&square)?.input_edge_count(), 4);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[inline]
    pub fn input_edge_count(&self) -> usize {
        self.edge_nodes.len()
    }

    /// Every input edge's face, in edge order.
    ///
    /// Returns `None` if any face cannot be walked.
    pub fn faces(&self) -> Option<Vec<Vec<NodeId>>> {
        (0..self.edge_nodes.len() as u16)
            .map(|i| self.face(EdgeId(i)))
            .collect()
    }

    /// Rebuilds the node-to-arc adjacency. Called once when the algorithm
    /// finishes.
    pub(crate) fn build_adjacency(&mut self) {
        self.node_arcs.clear();
        self.node_arcs.resize(self.nodes.len(), Vec::new());
        for (i, arc) in self.arcs.iter().enumerate() {
            let id = ArcId(i as u32);
            self.node_arcs[arc.nodes[0].0 as usize].push(id);
            // A degenerate zero-length arc would otherwise be listed twice.
            if arc.nodes[0] != arc.nodes[1] {
                self.node_arcs[arc.nodes[1].0 as usize].push(id);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn node(kind: NodeKind, offset: f32) -> Node {
        Node {
            position: Point::ORIGIN,
            exact: [0.0, 0.0],
            offset,
            kind,
            sources: vec![EdgeId(0), EdgeId(1)],
        }
    }

    #[test]
    fn arc_other_endpoint() {
        let a = Arc {
            nodes: [NodeId(1), NodeId(2)],
            sources: [EdgeId(0), EdgeId(1)],
        };
        assert_eq!(a.other(NodeId(1)), Some(NodeId(2)));
        assert_eq!(a.other(NodeId(2)), Some(NodeId(1)));
        assert_eq!(a.other(NodeId(3)), None);
        assert_eq!(a.lower(), NodeId(1));
        assert_eq!(a.upper(), NodeId(2));
    }

    #[test]
    fn node_kind_helpers() {
        let b = node(NodeKind::Boundary(VertexId(7)), 0.0);
        assert!(b.is_boundary());
        assert_eq!(b.input_vertex(), Some(VertexId(7)));

        let i = node(NodeKind::EdgeEvent, 3.0);
        assert!(!i.is_boundary());
        assert_eq!(i.input_vertex(), None);
    }

    #[test]
    fn adjacency_lists_every_incident_arc() {
        let mut s = Skeleton {
            nodes: vec![
                node(NodeKind::Boundary(VertexId(0)), 0.0),
                node(NodeKind::Boundary(VertexId(1)), 0.0),
                node(NodeKind::EdgeEvent, 5.0),
            ],
            arcs: vec![
                Arc {
                    nodes: [NodeId(0), NodeId(2)],
                    sources: [EdgeId(0), EdgeId(1)],
                },
                Arc {
                    nodes: [NodeId(1), NodeId(2)],
                    sources: [EdgeId(1), EdgeId(2)],
                },
            ],
            node_arcs: Vec::new(),
            edge_nodes: Vec::new(),
        };
        s.build_adjacency();

        assert_eq!(s.arcs_at(NodeId(0)), &[ArcId(0)]);
        assert_eq!(s.arcs_at(NodeId(1)), &[ArcId(1)]);
        assert_eq!(s.arcs_at(NodeId(2)), &[ArcId(0), ArcId(1)]);
    }

    #[test]
    fn max_offset_of_empty_skeleton_is_zero() {
        assert_eq!(Skeleton::default().max_offset(), 0.0);
    }

    #[test]
    fn max_offset_finds_the_ridge() {
        let mut s = Skeleton {
            nodes: vec![
                node(NodeKind::Boundary(VertexId(0)), 0.0),
                node(NodeKind::EdgeEvent, 5.0),
                node(NodeKind::EdgeEvent, 2.0),
            ],
            arcs: vec![],
            node_arcs: Vec::new(),
            edge_nodes: Vec::new(),
        };
        s.build_adjacency();
        assert_eq!(s.max_offset(), 5.0);
    }
}
