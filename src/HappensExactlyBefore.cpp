#include <fstream>
#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/raw_os_ostream.h>

#include "HappensAfterGraph.h"
#include "Basis.h"
#include "HappensExactlyBeforeGen.h"
#include "HappensExactlyBefore.h"
#include "DCTBHelpers.h" // isRead, etc.

Basis::events_iterator getFirstConflRW(const Basis::events_iterator& from,
                                          const DCEvent *ev,
                                          bool excluding_from = true) {
    auto it = from;
    auto end = from;
    end.shiftProcess();

    if (excluding_from)
      ++it;

    for (; it != end; ++it) {
        if (!isRead(**it) && !isWrite(**it))
          continue;

        if ((*it)->ml.overlaps(ev->ml))
          return it;
    }

    return from.end();
}

Basis::events_iterator getNextConflRW(const Basis::events_iterator& from) {
    return getFirstConflRW(from, *from);
}

static std::vector<std::pair<HappensAfterGraph::Node *, HappensAfterGraph::Node *>>
updateHappensAfterGraph(HappensAfterGraph& G,
                        const std::pair<const DCEvent *,
                                        const DCEvent *>& edge) {
    HappensAfterGraph::Node *n1 = G.getNode(edge.first);
    HappensAfterGraph::Node *n2 = G.getNode(edge.second);
    assert(n1 && n2);
    std::vector<std::pair<HappensAfterGraph::Node *,
                          HappensAfterGraph::Node *>> new_edges;

    // the graph already has this edge, no work for us
    assert(!G.hasEdge(n2, n1) && "Trying to add a forbidden edge");
    if(G.hasEdge(n1, n2))
      return new_edges;

    // the graph is transitively closed, so we can just iterate over
    // predecessors and successors
    for (auto pred : n1->getPredecessors()) {
      if (!G.hasEdge(pred, n2)) {
        new_edges.push_back({pred, n2});
        pred->addSuccessor(n2);
      }

      for (auto succ : n2->getSuccessors()) {
        if (!G.hasEdge(pred, succ)) {
          new_edges.push_back({pred, succ});
          pred->addSuccessor(succ);
        }
      }
    }

    for (auto succ : n2->getSuccessors()) {
      if (!G.hasEdge(n1, succ)) {
        new_edges.push_back({n1, succ});
        n1->addSuccessor(succ);
      }
    }

    new_edges.push_back({n1, n2});
    n1->addSuccessor(n2);

    assert(new_edges.size() > 0 || edge.first == edge.second);
    assert(G.isTransitivelyClosed());
    return new_edges;
}

static void rollbackHappensAfterGraph(HappensAfterGraph& G,
                                      const std::vector<std::pair<HappensAfterGraph::Node *,
                                                                  HappensAfterGraph::Node *>>& edges) {
  for (auto& edge : edges) {
    edge.first->removeSuccessor(edge.second);
  }
}

Basis::events_iterator
HappensExactlyBeforeGen::needsHB(HappensExactlyBefore& happens_before,
                                 const Basis::events_iterator& start) {
  Basis::events_iterator it = start;
  for (auto end = basis.events_end(); it != end; ++it) {
    // skip root process
    if (basis.isTopologyRoot(it.getProcessID())) {
        it.skipProcess();
        continue;
    }

    assert(!basis.isTopologyRoot(it.getProcessID()));

    const DCEvent *event = *it;
    if (!isWrite(*event) && !isRead(*event))
      continue;

    const HappensExactlyBefore::ValueT *HB
        = happens_before.get({event->cpid, event->instruction, event->order});
    // if have the happens-before relation for all other threads,
    // (except the root) we can skip this write/read
    assert(!HB || HB->size() <= basis.size() - 2);
    if (HB && HB->size() == basis.size() - 2)
      continue;

    // if we got here, we found a write/read that needs HB
    return it;
  }

  assert(it == basis.events_end());
  return it;
}

#if DUMP_HB_TREE
void HappensExactlyBeforeGen::_generate(HappensExactlyBefore& happens_before,
                                        const Basis::events_iterator& start,
                                        unsigned tree_parent_node)
#else
void HappensExactlyBeforeGen::_generate(HappensExactlyBefore& happens_before,
                                        const Basis::events_iterator& start)
#endif
{
  assert(basis.size() > 2);
  assert(G.isAcyclic());

  // Find the event that needs the happens-before relation.
  // If there is none such event, we are done.
  Basis::events_iterator it1 = needsHB(happens_before, start);
  if (it1 == basis.events_end()) {
#if DUMP_HB_TREE
    tree.addNode(tree_parent_node, happens_before, "DONE");
#endif
    generated_hb.push_back(std::move(happens_before));
    return;
  }

  const DCEvent *event = *it1;
  assert(isWrite(*event) || isRead(*event));
  const HappensExactlyBefore::ValueT *HB
    = happens_before.get({event->cpid, event->instruction, event->order});

  // we have incomplete HB for the 'event' (it1), now find
  // the thread for which we are missing the information
  for (auto it2 = basis.events_begin(), end2 = basis.events_end();
       it2 != end2; ++it2) {
    // do not add the relation in the same process,
    // it makes no sense. Also do not make happens-before when one process
    // is root, we are observation-equivalent on such processes.
    if (it1.getProcessID() == it2.getProcessID()
        || basis.isTopologyRoot(it2.getProcessID())) {
      it2.skipProcess();
      continue;
    }

    // we already have the happens before
    // relation for this thread
    if (HB && HB->find((*it2)->cpid) != HB->end()) {
      it2.skipProcess();
      continue;
    }

    // generate the happens-before relation for the thread
    // into which belongs the event it2
#if DUMP_HB_TREE
    _generate(happens_before, it1, it2, tree_parent_node);
#else
    _generate(happens_before, it1, it2);
#endif
    return;
  }

  assert(0 && "Not reachable");
}

static const llvm::Instruction *hbGetConstrains(const Basis& basis,
                                                const Basis::events_iterator& it1,
                                                const Basis::events_iterator& it2,
                                                HappensExactlyBefore& happens_before)
{
  // If we have a conflicting event in this thread above this one,
  // then this read/write is constrained by the HB of the above one
  unsigned i = it1.getProcessID();
  unsigned j = it1.getEventID();
  const HappensExactlyBefore::ValueT *constrains = nullptr;

  while (j-- > 0) {
    if (basis[i][j]->ml.overlaps((*it1)->ml)) {
      constrains = happens_before.get(*basis[i][j]);
      assert(!constrains || constrains->size() == basis.size() - 2);
      break;
    }
  }

  if (constrains) {
      auto C = constrains->find((*it2)->cpid);
      if (C != constrains->end())
        return C->second.first;
  }

  return nullptr;
}

#if DUMP_HB_TREE
void HappensExactlyBeforeGen::_generate(HappensExactlyBefore& happens_before,
                                        const Basis::events_iterator& it1,
                                        Basis::events_iterator& it2,
                                        unsigned tree_parent_node)
#else
void HappensExactlyBeforeGen::_generate(HappensExactlyBefore& happens_before,
                                        const Basis::events_iterator& it1,
                                        Basis::events_iterator& it2)
#endif
{
  assert(it1.getProcessID() != it2.getProcessID());
  assert(it2.getEventID() == 0);

  // we do not have happens-before relation for this thread?
  // then we must take all the events that are lower then
  // HB of the previous conflicting event in this thread and try them.
  // go_from_instruction will tell us whether we should go through
  // the whole thread or whether we can start somewhere in the middle.
  const llvm::Instruction *go_from_instruction
    = hbGetConstrains(basis, it1, it2, happens_before);

  // find the end of this thread, we'll need that
  // when going through the events
  auto process2_end = it2.nextProcess();
  const DCEvent *event = *it1;

  if (go_from_instruction == nullptr) {
    // we should go from the beginning of the thread,
    // so try also adding that the event sees the initial event

    // if some (conflicting) event from the other thread sees the init event,
    // then there's no point in adding this, it will be unralizable
    assert(it2.getEventID() == 0);
    auto confl = getFirstConflRW(it2, event, false /* excluding from */);
    bool add_initial = true;
    if (confl != basis.events_end()) {
      assert(isRead(**confl) || isWrite(**confl));
      // I'm trying to add the add event->confl, so check,
      // whether we do not have the oposite edge already
      // (event sees init means that there will be and edge from
      // init->event and event->confl. Beacuse the graph is
      // transitively closed, than whenever confl sees something
      // above event -- which would create a loop -- then there
      // is also the edge confl->event, so it is enough to check
      // only for this one)
      if (G.hasEdge(*confl, event)) {
#if DUMP_HB_TREE
        // XXX DEBUGGING
        HappensExactlyBefore new_happens_before = happens_before;
        new_happens_before.add(*event, (*it2)->cpid, nullptr, 0);
        tree.addNode(tree_parent_node, new_happens_before, "BLOCKED (239)");
#endif
        add_initial = false;
      }
    }

    if (add_initial) {
        HappensExactlyBefore new_happens_before = happens_before;
        std::vector<std::pair<HappensAfterGraph::Node *,
                              HappensAfterGraph::Node *>> new_edges;

        // NOTE: do not use initial_event as such, it has invalid cpid,
        // we need correct initial event for every pid
        assert(!basis.isTopologyRoot(it2.getProcessID()));
        new_happens_before.add(*event, (*it2)->cpid, nullptr, 0);

#if DUMP_HB_TREE
        // XXX: DEBUGGING
        unsigned tnd = tree.addNode(tree_parent_node, new_happens_before);
        _generate(new_happens_before, it1, tnd);
#else
        _generate(new_happens_before, it1);
#endif
        rollbackHappensAfterGraph(G, new_edges);
    }
  }/*
   XXX: caused a bug in tests/test1-false.c. Probably something wrong
   in hbGetConstrains(basis, it1, it2, happens_before);
  else {
    // we should go from some particular instruction, so find it
    for (; it2 != process2_end; ++it2) {
      if ((*it2)->instruction != go_from_instruction)
        continue;
    }

    assert(it1 != process2_end && "Should go from instr, but did find it");
  }
  */

  // we now want to go from current position to the end of the process
  for (; it2 != process2_end; ++it2) {
    const DCEvent *event2 = *it2;
    bool is_read2 = isRead(*event2);

    if (!isWrite(*event2) && !is_read2)
      continue;

    // two reads are independent
    if (isRead(**it1) && is_read2)
      continue;

    // if those events are not overlapping, we don't care
    if (!event->ml.overlaps(event2->ml))
      continue;

    assert(it1.getProcessID() != it2.getProcessID());
    assert(!basis.isTopologyRoot(it2.getProcessID()));

    // we found a conflicting event.
    // We're about to add 'event (sees) event2', so check
    // whether we do not have 'event2 (sees) event'.
    if (auto V = happens_before.get(*event2)) {
      const auto& vit = V->find(event->cpid);
      if (vit != V->end() &&
          vit->second == HappensExactlyBefore::InstrIdT(event->instruction, event->order)) {
        ++blocked_hb;
        // we can bail out, because no other event that follows event2 can not
        // be seen by event in this case

#if DUMP_HB_TREE
        // XXX DEBUGGING
        HappensExactlyBefore new_happens_before = happens_before;
        new_happens_before.add(*event, *event2);
        tree.addNode(tree_parent_node, new_happens_before, "BLOCKED (for good)");
#endif
        return;
      }
    }

    // if we have the edge (path) event->event2, then this will be unrealizable
    if (G.hasEdge(event, event2)) {
      ++blocked_hb;

#if DUMP_HB_TREE
      // XXX DEBUGGING
      HappensExactlyBefore new_happens_before = happens_before;
      new_happens_before.add(*event, *event2);
      tree.addNode(tree_parent_node, new_happens_before, "BLOCKED (285)");
#endif
      continue;
    }

    std::vector<std::pair<HappensAfterGraph::Node *,
                          HappensAfterGraph::Node *>> new_edges;
    auto confl = getNextConflRW(it2);
    // there's some confl event
    if (confl != basis.events_end()) {
      assert(isRead(**confl) || isWrite(**confl));

      if (G.hasEdge(*confl, event)) {
        ++blocked_hb;

#if DUMP_HB_TREE
        // XXX DEBUGGING
        HappensExactlyBefore new_happens_before = happens_before;
        new_happens_before.add(*event, *event2);
        tree.addNode(tree_parent_node, new_happens_before, "BLOCKED (298)");
#endif
        continue;
      }

      new_edges = updateHappensAfterGraph(G, {event, *confl});
    }

    // new edges that we add to the happens-after graph
    auto new_edges2  = updateHappensAfterGraph(G, {event2, event});

    HappensExactlyBefore new_happens_before = happens_before;
    // event sees event2
    new_happens_before.add(*event, *event2);

    // call it recursively to generate the relation
    // for the rest or writes
#if DUMP_HB_TREE
    unsigned tnd = tree.addNode(tree_parent_node, new_happens_before);
    _generate(new_happens_before, it1, tnd);
#else
    _generate(new_happens_before, it1);
#endif

    rollbackHappensAfterGraph(G, new_edges);
    rollbackHappensAfterGraph(G, new_edges2);
    assert(G.isTransitivelyClosed());
  }
}

extern llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const DCIID& iid);

llvm::raw_ostream& operator<<(llvm::raw_ostream& out,
                              const HappensExactlyBefore& hb) {
  out << "Happens-Before- {\n";
  for (auto& pr : hb) {
    out << pr.first << " =>\n  {\n";
    for (auto &it : pr.second) {
      out << "  " << it.first << "-" << it.second.second;
      if (it.second.first)
        out << *it.second.first << ")\n";
      else
        out << " initial)\n";
    }
    out << "  }\n";
  }
  out << "}\n";

  return out;
}

void HappensExactlyBefore::dump() const {
    llvm::errs() << *this;
}

#ifdef DUMP_HB_TREE
void HBTree::dump() const {
  if (nodes.empty())
      return;

  static unsigned num = 0;
  std::string name = "hb" + std::to_string(num++) + ".dot";
  std::ofstream _annot_debug(name);
  llvm::raw_os_ostream annot_debug(_annot_debug);

  annot_debug << "digraph HappensExactlyBefore {\n";

  // nodes
  for (auto& nd : nodes) {
    // the node description
    annot_debug << "NODE" << nd.id << " [label=\"";
    annot_debug << nd.id << "\\n";
    annot_debug << nd.hb;
    if (nd.msg)
        annot_debug << "\\n(" << nd.msg << ")";
    annot_debug << "\"";

    if (!nd.msg && nd.successors.empty())
        annot_debug << " style=filled fillcolor=green";

    annot_debug << "]\n";
  }

  // edges
  for (auto& nd : nodes) {
    for (unsigned succ : nd.successors) {
      // a new edge
      annot_debug << "NODE" << nd.id << "->NODE" << succ << "\n";
    }
  }

  annot_debug << "}\n";
  annot_debug.flush();
  _annot_debug.close();
}
#endif

