#ifndef _HAPPENS_AFTER_GRAPH_H_
#define _HAPPENS_AFTER_GRAPH_H_

#include <set>
#include <vector>
#include <unordered_map>

#include "Event.h"
#include "Trace.h"
#include "Annotations.h"
#include "Basis.h"
#include "2SAT.h"
#include "CPid.h"

class HappensExactlyBefore;

class HappensAfterGraph {
  const PositiveAnnotation& annotation;
  const HappensExactlyBefore& happens_before;
  const Basis& basis;
  unsigned max_search_index = 0;
  bool blocked_before_sat = false;

public:
  HappensAfterGraph(const PositiveAnnotation& annot,
                    const HappensExactlyBefore& happens_before,
                    const Basis& b)
  : annotation(annot), happens_before(happens_before), basis(b)
  {
    _constructGraph();
    reachable_nodes.reserve(nodes.size());
    rev_reachable_nodes.reserve(nodes.size());
  }

  HappensAfterGraph(const HappensAfterGraph& oth,
                    const HappensExactlyBefore& hb);

  HappensAfterGraph(HappensAfterGraph&& oth);
  HappensAfterGraph& operator=(HappensAfterGraph&& oth);
  ~HappensAfterGraph();

  const PositiveAnnotation& getAnnotation() const { return annotation; }
  const Basis& getBasis() const { return basis; }
  bool wasBlocked() const { return blocked_before_sat; }

  class Node {
    using EdgesT = std::set<Node *>;
    EdgesT successors;
    EdgesT predecessors;

    // for searching cycles and SCC
    // FIXME: move this out of the node to auxiliary mappings
    unsigned dfsid, lowpt;
    bool is_on_stack;
    // unique id of a node,
    // serves as the index to the vector of nodes
    unsigned id = 0;
    // indices into basis
    unsigned process_id = 0;
    unsigned event_id = 0;

  public:
    const DCEvent *event;
    Node(const DCEvent *ev)
      : dfsid(0), lowpt(0), is_on_stack(0), event(ev) {}

    void set_on_stack(bool a) { is_on_stack = a; }
    void set_dfsid(unsigned id) { dfsid = id; }
    void set_lowpt(unsigned lp) { lowpt = lp; }

    unsigned get_dfsid() const { return dfsid; }
    unsigned get_lowpt() const { return lowpt; }
    bool on_stack() const { return is_on_stack; }

    void addSuccessor(Node *nd) {
      assert(nd != this && "We do not want self-loops");
      assert(nd && "Given node is nullptr");
      successors.insert(nd);
      nd->predecessors.insert(this);
    }

    void removeSuccessor(Node *nd) {
      bool a = successors.erase(nd);
      assert(a && "Didn't have this successor");
      bool b = nd->predecessors.erase(this);
      assert(a == b && "Had successor, but no predecessor");
    }

    EdgesT& getSuccessors() { return successors; }
    const EdgesT& getSuccessors() const { return successors; }
    EdgesT& getPredecessors() { return predecessors; }
    const EdgesT& getPredecessors() const { return predecessors; }

    friend class HappensAfterGraph;
  };

private:
  // a cache for reachable nodes (instead of doing the transitive closure)
  std::unordered_map<Node *, std::set<Node *>> reachable_nodes;
  // and the reverse relation
  std::unordered_map<Node *, std::set<Node *>> rev_reachable_nodes;

  PositiveAnnotation
  getPastConeAnnotation(const PositiveAnnotation& O,
                        std::vector<const Node *> &&tp) const;
public:
  std::set<Node *>& getNodesReachableFromNoLazy(Node *n);
  std::set<Node *>& getNodesReachingToNoLazy(Node *n);

  std::set<Node *>& getNodesReachableFrom(Node *n) {
    auto& ret = reachable_nodes[n];
    if (!ret.empty() || n->successors.empty())
      return ret;

    return getNodesReachableFromNoLazy(n);
  }

  std::set<Node *>& getNodesReachingTo(Node *n) {
    auto& ret = rev_reachable_nodes[n];
    if (!ret.empty() || n->predecessors.empty())
      return ret;

    return getNodesReachingToNoLazy(n);
  }

  std::set<Node *>& getNodesReachableFrom(const DCEvent *e) {
    Node *tmp = getNode(e);
    assert(tmp && "Invalid event");
    return getNodesReachableFrom(tmp);
  }

  std::set<Node *>& getNodesReachableFromNoLazy(const DCEvent *e) {
    Node *tmp = getNode(e);
    assert(tmp && "Invalid event");
    return getNodesReachableFromNoLazy(tmp);
  }

  const Node *getNode(const DCEvent *e) const {
    auto it = nodes.find(e);
    if (it == nodes.end())
      return nullptr;

    return it->second;
  }

  Node *getNode(const DCEvent *e) {
    auto it = nodes.find(e);
    if (it == nodes.end())
      return nullptr;

    return it->second;
  }

  Node *getNode(const AnnotationKeyT& iid) const {
      return mapAnnotToNode(iid);
  };

  unsigned addNode(Node *nd) {
      nd->id = nodes.size();
      nodes.emplace(nd->event, nd);
      return nd->id;
  }

  // add edges between reads and write from the annotation
  void addAnnotation(const PositiveAnnotation& annotation);

  PositiveAnnotation
  getPastConeAnnotation(const PositiveAnnotation& O,
                        const std::initializer_list<const DCEvent *>& l) const;

  PositiveAnnotation
  getPastConeAnnotation(const PositiveAnnotation& O,
                        const std::initializer_list<const DCIID *>& l) const;

  // check whether ev2 is reachable from ev1
  bool canReach(const DCEvent *ev1, const DCEvent *ev2) const;
  bool canReach(const DCEvent& ev1, const DCEvent& ev2) const {
    return canReach(&ev1, &ev2);
  }

  bool hasEdge(const DCEvent *e1, const DCEvent *e2) {
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    assert(n1 && n2 && "Do not have such node");
    return n1->successors.count(n2) == 1;
  }

  bool hasEdge(Node *n1, Node *n2) {
    assert(n1 && n2 && "Do not have such node");
    return n1->successors.count(n2) == 1;
  }

  size_t size() const { return nodes.size(); }

  bool isAcyclic();
  bool checkGraph();
  // for two threads we do not need the 2SAT
  bool addNecessaryEdges();
  void makeTransitiveClosure();
  void makeTransitiveClosureSimple();

  Basis getRestrictedBasis(const PositiveAnnotation& annot) const;

  // for debugging
  bool isTransitivelyClosed();

  // linearize the PO induced by this graph
  // to a trace
  std::vector<DCEvent> linearize();

  void to_dot(const char *edge_params=nullptr) const;
  void dump() const;

  class nodes_iterator {
    std::unordered_map<const DCEvent *, Node *>::iterator it;
    nodes_iterator(std::unordered_map<const DCEvent *, Node *>::iterator&& i)
    : it(i) {}

    friend class HappensAfterGraph;
  public:
    nodes_iterator& operator++() { ++it; return *this;}
    nodes_iterator operator++(int) { auto tmp = *this; ++it; return tmp;}
    bool operator==(const nodes_iterator& oth) { return it == oth.it; }
    bool operator!=(const nodes_iterator& oth) { return it != oth.it; }
    Node& operator*() { return *it->second; }
  };

  nodes_iterator begin() { return nodes_iterator(nodes.begin()); }
  nodes_iterator end() { return nodes_iterator(nodes.end()); }

private:
  void _constructGraph();

  std::unordered_map<const DCEvent *, Node *> nodes;

  // create a mapping from instructions to nodes.
  // We could have mapping from (CPid, Instruction) -> Node,
  // but the comparsion of CPid would be too unefficient for our purposes
  // right now...
  std::unordered_map<const llvm::Instruction *, std::vector<Node *>> instr_to_node;
  Node *initial_node = nullptr;

  Node *mapAnnotToNode(const AnnotationKeyT& iid) const;
};

#endif // _HAPPENS_AFTER_GRAPH_H_
