#ifndef _HAPPENS_AFTER_GRAPH_VCLOCK_H_
#define _HAPPENS_AFTER_GRAPH_VCLOCK_H_

#include <vector>
#include <unordered_map>

#include "Basis.h"
#include "Event.h"
#include "Trace.h"
#include "Annotations.h"
#include "2SAT.h"
#include "CPid.h"

class HappensExactlyBefore;
class HappensExactlyBefore;

namespace HAG {
// node represents a node in basis and keeps the information
// about happens-after relation between the nodes
struct Node {
  // these are indices into basis.
  // If they are INT_MAX, INT_MAX, then this is the initial node
  unsigned process_id;
  unsigned event_id;

  std::vector<int> successors;
  const Basis& basis;

  Node(const Basis& basis, unsigned pid, unsigned evid)
    : process_id(pid), event_id(evid), basis(basis)
  {
    // initialize the successors to INT_MAX, which means
    // there is no edge to those threads
    successors.resize(basis.size(), INT_MAX);
  }

  Node(const Node& oth) = default;
  Node& operator=(const Node& oth) = default;

  const std::vector<int>& getSuccessors() { return successors; }

  // FIXME: move this function to HappensAfterGraphVclock
  // and get rid of the basis attribute
  const DCEvent *getEvent() const {
    if (process_id == INT_MAX)
      return nullptr;

    assert(process_id < basis.size());
    assert(event_id < basis[process_id].size());
    return basis[process_id][event_id];
  }

  void addSuccessor(Node *nd) {
    assert(nd != this && "We do not want self-loops");
    assert(nd->process_id < successors.size());
    successors[nd->process_id]
      = std::min(static_cast<unsigned>(successors[nd->process_id]), nd->event_id);
    assert(successors[nd->process_id] != INT_MAX);
  }
};
} // namespace HAG

class HappensAfterGraphVclock : private Processes<HAG::Node *> {
  using Node = HAG::Node;
  // every event has id that is the index into the trace,
  // let's use this id as a key in a map
  using MappingT = std::unordered_map<const DCEvent *, Node *>;

  const PositiveAnnotation& annotation;
  const HappensExactlyBefore& happens_before;
  const Basis& basis;
  bool blocked_before_sat = false;

  void reserveMemory() {
    if (basis.size() == 0)
      return;

    // this should be more efficient than to
    // dynamically adjust the hash table
    unsigned reserve_size = 0;
    for (auto&b : basis)
      reserve_size += b.size();
    instr_to_node.reserve(reserve_size);
    mapping.reserve(reserve_size);
  }

  void createInitialNode() {
    initial_node = new Node(basis, INT_MAX, INT_MAX);
    // the initial node is before all other events
    for (int& x : initial_node->successors)
      x = 0;
    mapping.emplace(nullptr, initial_node);
  }

public:
  HappensAfterGraphVclock(const PositiveAnnotation& annot,
                          const HappensExactlyBefore& happens_before,
                          const Basis& b)
  : annotation(annot), happens_before(happens_before), basis(b)
  {
    reserveMemory();
    createInitialNode();

    if (basis.size() > 0)
      _constructGraph();
  }

  HappensAfterGraphVclock(HappensAfterGraphVclock&& oth);
  HappensAfterGraphVclock& operator=(HappensAfterGraphVclock&& oth);

  HappensAfterGraphVclock(const HappensAfterGraphVclock& oth);

  ~HappensAfterGraphVclock();

  bool operator==(const HappensAfterGraphVclock&) const;

  bool wasBlocked() const { return blocked_before_sat; }

  // returns {true, order} if the graph is acyclic and
  // {false, _} if there's a cycle
  std::pair<bool, std::vector<Node *>> computeTopoOrder();

  std::pair<PositiveAnnotation, std::set<const DCEvent *>>
  getUnorderedEvents() const;

  Basis createSubbasis(std::set<const DCEvent *>& new_nodes) const;

  // create a subgraph induced by the undecided edges
  // -- we can create a formula only for this subgraph
  HappensAfterGraphVclock createSubgraph(const PositiveAnnotation& annot,
                                         const Basis& new_basis) const;

  const Node *getNode(const DCEvent *e) const {
    auto it = mapping.find(e);
    if (it == mapping.end())
      return nullptr;

    return it->second;
  }

  Node *getNode(const DCEvent *e) {
    auto it = mapping.find(e);
    if (it == mapping.end())
      return nullptr;

    return it->second;
  }

  bool hasEdge(const DCEvent *e1, const DCEvent *e2) const {
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    assert(n1 && n2 && "Do not have such node");
    return hasEdge(n1, n2);
  }

  bool hasEdge(const Node *n1, const Node *n2) const {
    // no node can have an edge to init
    if (n2->process_id == INT_MAX)
      return false;
    return (n1->successors[n2->process_id] != INT_MAX
            && (unsigned) n1->successors[n2->process_id] <= n2->event_id);
  }

  // has edge n1->n2 or n2->n1
  bool hasConnected(const Node *n1, const Node *n2) const {
      return hasEdge(n1, n2) || hasEdge(n2, n1);
  }

  bool hasConnected(const DCEvent *e1, const DCEvent *e2) const {
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    assert(n1 && n2 && "Do not have such node");

    return hasConnected(n1, n2);
  }

  void addEdge(const DCEvent *e1, const DCEvent *e2) {
    assert(e2 && "An edge to init would create a cycle");
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    assert(n1 && n2 && "Do not have such node");
    addEdge(n1, n2);
  }

  // add an edge, return false is the edge creates a cycle
  void addEdge(Node *n1, Node *n2) {
    n1->addSuccessor(n2);
  }
  size_t size() const { return mapping.size(); }

  bool isAcyclic() {
    return makeTransitiveClosure();
  }

  bool addNecessaryEdges();
  bool addNecessaryEdgesFromAnnotation();
  bool addNecessaryEdgesFromAnnotationClauses();
  bool addNecessaryEdgesFromHB();
  // return false if the graph is acyclic
  bool makeTransitiveClosure();

  bool addTransitivity(Events2SAT& sat, const DCEvent *e1, const DCEvent *e2);
  bool addAnnotationClauses(Events2SAT& sat);

  // create the 2SAT formula from the graph.
  // Returns false if some contradiction was detected
  // already while creating
  bool createFormula(Events2SAT& sat);
  // add new edges in such manner that the
  // partial order induced by this graph
  // realizes the annotation
  // return false if the annotation is unrealizable
  bool realize();

  // linearize the PO induced by this graph to a trace
  std::vector<DCEvent> linearize();

  void to_dot(const char *edge_params=nullptr) const;
  void dump() const;

  class nodes_iterator {
    MappingT::iterator it;
    nodes_iterator(MappingT::iterator&& i) : it(i) {}

    friend class HappensAfterGraphVclock;
  public:
    nodes_iterator& operator++() { ++it; return *this;}
    nodes_iterator operator++(int) { auto tmp = *this; ++it; return tmp;}

    bool operator==(const nodes_iterator& oth) const { return it == oth.it; }
    bool operator!=(const nodes_iterator& oth) const { return it != oth.it; }
    Node& operator*() { return *it->second; }
    const Node& operator*() const { return *it->second; }
  };

  // inherit from nodes_iterator or make it via templates and typedefs/inheritance
  class nodes_const_iterator {
    MappingT::const_iterator it;
    nodes_const_iterator(MappingT::const_iterator&& i) : it(i) {}

    friend class HappensAfterGraphVclock;
  public:
    nodes_const_iterator& operator++() { ++it; return *this;}
    nodes_const_iterator operator++(int) { auto tmp = *this; ++it; return tmp;}

    bool operator==(const nodes_const_iterator& oth) const { return it == oth.it; }
    bool operator!=(const nodes_const_iterator& oth) const { return it != oth.it; }
    const Node& operator*() const { return *it->second; }
  };

  nodes_iterator begin() { return nodes_iterator(mapping.begin()); }
  nodes_iterator end() { return nodes_iterator(mapping.end()); }

  nodes_const_iterator begin() const { return nodes_const_iterator(mapping.begin()); }
  nodes_const_iterator end() const { return nodes_const_iterator(mapping.end()); }

  struct edge_iterator {
    HappensAfterGraphVclock *processes;
    // which node we use?
    Node *nd;
    // iterator to the current event
    HappensAfterGraphVclock::events_iterator event_it;

    edge_iterator(HappensAfterGraphVclock *h, Node *nd) :
      processes(h), nd(nd), event_it(processes->events_end())
    {}

    edge_iterator(const edge_iterator&) = default;
    edge_iterator& operator=(const edge_iterator&) = default;

    Node *operator*() const {
      return *event_it;
    }

    bool operator==(const edge_iterator& oth) const {
      return event_it == oth.event_it;
    }

    bool operator!=(const edge_iterator& oth) const {
      return event_it != oth.event_it;
    }
  };

  class succ_iterator : public edge_iterator {
   // the end iterator of current thread
    HappensAfterGraphVclock::events_iterator thread_end;

    void set_to(unsigned thr, unsigned ev) {
      event_it = processes->getIterator(thr, ev); 
      thread_end = event_it.nextProcess();
    }

    int find_next_thread(int current_thread) {
        assert(current_thread >= 0);
        for (unsigned i = current_thread; i < nd->successors.size(); ++i) {
          // find the first successor
          if (nd->successors[i] != INT_MAX) {
            return i;
          }
        }

        return -1;
    }

    void set_next(int current_thread) {
      int nxt = find_next_thread(current_thread);
      if (nxt != -1) {
        set_to(nxt, nd->successors[nxt]);
        assert(event_it.getProcessID() < processes->size());
        assert(event_it.getEventID() < (*processes)[event_it.getProcessID()].size());
      } else {
        // we haven't found anything,
        // so set the iterator to the end
        thread_end = event_it = processes->events_end();
      }
    }

    succ_iterator(HappensAfterGraphVclock *h, Node *nd, bool end)
      : edge_iterator(h, nd), thread_end(processes->events_end())
    {}

    succ_iterator(HappensAfterGraphVclock *h, Node *nd)
      : edge_iterator(h, nd), thread_end(processes->events_end()) {
        set_next(0);
      }

    friend class HappensAfterGraphVclock;

    public:
      succ_iterator(const succ_iterator& oth) = default;
      succ_iterator& operator=(const succ_iterator&) = default;

      succ_iterator& operator++() {
        assert(event_it != processes->events_end());
        assert(event_it != thread_end);

        ++event_it;

        if (event_it == thread_end) {
          // event_it is now already in the next thread,
          // so start searching from getProcessID()
          set_next(event_it.getProcessID());
        }

        return *this;
      }

      succ_iterator operator++(int) {
        auto tmp = *this;
        operator++();
        return tmp;
      }
  };

  // iterator over predecessors of a node in
  // a graph that is transitively closed
  class pred_iterator : public edge_iterator {
    void set_to(unsigned thr, unsigned ev) {
      event_it = processes->getIterator(thr, ev); 
    }

    int find_next_thread(int current_thread) {
      assert(current_thread >= -1);
      for (unsigned i = current_thread + 1; i < nd->successors.size(); ++i) {
        unsigned tmp = (*processes)[i][0]->successors[nd->process_id];
        if (tmp != INT_MAX && tmp <= nd->event_id) {
          return i;
        }
      }

      return -1;
    }

    void set_next(int current_thread) {
        int nxt = find_next_thread(current_thread);
        if (nxt != -1) {
          set_to(nxt, 0);
        } else {
          // we haven't found anything,
          // so set the iterator to the end
          event_it = processes->events_end();
        }
    }

    pred_iterator(HappensAfterGraphVclock *h, Node *nd, bool end)
     : edge_iterator(h, nd) {}

    pred_iterator(HappensAfterGraphVclock *h, Node *nd)
      : edge_iterator(h, nd) {
      if (nd->process_id != INT_MAX)
        set_next(-1);
      // else keep it on end, since init can not have any predecessors
    }

    friend class HappensAfterGraphVclock;

    public:
      pred_iterator(const pred_iterator& oth) = default;
      pred_iterator& operator=(const pred_iterator&) = default;

      pred_iterator& operator++() {
        assert(event_it != processes->events_end());

        ++event_it;

        if (event_it != processes->events_end()) {
          Node *cur = *event_it;
          if (cur == nd || (unsigned)cur->successors[nd->process_id] > nd->event_id) {
            set_next(event_it.getProcessID());
          }
        }

        return *this;
      }

      pred_iterator operator++(int) {
        auto tmp = *this;
        operator++();
        return tmp;
      }
  };

  succ_iterator succ_begin(Node *nd) {
    return succ_iterator(this, nd);
  }

  succ_iterator succ_end(Node *nd) {
    return succ_iterator(this, nd, true /* end */);
  }

  // works correctly only in transitively closed graph!
  pred_iterator pred_begin(Node *nd) {
    return pred_iterator(this, nd);
  }

  // works correctly only in transitively closed graph!
  pred_iterator pred_end(Node *nd) {
    return pred_iterator(this, nd, true /* end */);
  }

private:
  // wrapper around add implication:
  // representing this: (e11, e21) => (e21, e22)
  // if fst or snd, resp. is set to false, it means that
  // the first part, second part resp. is negated,
  // e.g. for fst = false:
  // !(e11, e21) => (e21, e22)
  //
  // This wrapper removes the need to have antisymmetry clauses,
  // becase when we have X12 then instead of variable X21
  // we will use !X12. We keep only the variables for which
  // e1 < e2.
  //
  // Example: for addImplication(1, 2, 3, 1, true, true) we will
  // call addImplication(1, 2, 1, 3, true, false)
  // (instead of X31 we will use NOT X13)
  bool addImplication(Events2SAT& sat,
                      const DCEvent *e11, const DCEvent *e12,
                      const DCEvent *e21, const DCEvent *e22,
                      bool fst = true, bool snd = true) {
      if (e12 < e11) {
          // swap the events
          const DCEvent *tmp = e12;
          e12 = e11;
          e11 = tmp;

          // negate the variable
          fst = !fst;
      }

      if (e22 < e21) {
          // swap the events
          const DCEvent *tmp = e22;
          e22 = e21;
          e21 = tmp;

          // negate the variable
          snd = !snd;
      }

      // can we decide it already?
      if (fst ? hasEdge(e11, e12) : hasEdge(e12, e11))
          return sat.addAssignment(e21, e22, snd);
      else
        return sat.addImplication(e11, e12, e21, e22, fst, snd);
  }

  bool addAssignment(Events2SAT& sat,
                     const DCEvent *e1, const DCEvent *e2, bool val) {
      if (e2 < e1) {
          // swap the events
          const DCEvent *tmp = e2;
          e2 = e1;
          e1 = tmp;

          // negate the variable
          val = !val;
      }

      return sat.addAssignment(e1, e2, val);
  }

  MappingT mapping;
  // create a mapping from instructions to nodes.
  // We could have mapping from (CPid, Instruction) -> Node,
  // but the comparsion of CPid would be too unefficient for our purposes
  // right now...
  std::unordered_map<const llvm::Instruction *, std::vector<Node *>> instr_to_node;
  Node *initial_node = nullptr;

  void _constructGraph();
  Node *mapAnnotToNode(const AnnotationKeyT& iid) const;
};

#endif // _HAPPENS_AFTER_GRAPH_H_
