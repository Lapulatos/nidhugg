#include <unordered_map>
#include <stack>

#include "DCTBHelpers.h" // isWrite
#include "HappensAfterGraphVclock.h"
#include "2SAT.h"
#include "SCC.h"

HappensAfterGraphVclock::~HappensAfterGraphVclock()
{
  for (auto& it : *this)
    delete &it;
}

HappensAfterGraphVclock::HappensAfterGraphVclock(HappensAfterGraphVclock&& oth)
: annotation(std::move(oth.annotation)),
  happens_before(std::move(oth.happens_before)),
  basis(oth.basis),
  blocked_before_sat(oth.blocked_before_sat),
  mapping(std::move(oth.mapping)),
  instr_to_node(std::move(oth.instr_to_node)),
  initial_node(oth.initial_node)
{
}

HappensAfterGraphVclock::HappensAfterGraphVclock(const HappensAfterGraphVclock& oth)
: annotation(oth.annotation),
  happens_before(oth.happens_before),
  basis(oth.basis),
  blocked_before_sat(oth.blocked_before_sat)
{
  mapping.reserve(oth.size());
  instr_to_node.reserve(instr_to_node.size());

  for (auto& it : oth.mapping) {
    // copy the node (successors are indices into basis,
    // so the successors will be mapped fine)
    assert(it.second);
    Node *new_nd = new Node(*it.second);

    if (it.first && it.first->instruction)
      instr_to_node[it.first->instruction].emplace_back(new_nd);

    mapping.emplace(it.first, new_nd);
  }

  // also set the initial_node
  initial_node = mapping[nullptr];
  assert(initial_node);
  assert(*this == oth);
}

HappensAfterGraphVclock& HappensAfterGraphVclock::operator=(HappensAfterGraphVclock&& oth) {
  if (&oth != this)
    *this = std::move(oth);
  return *this;
}

bool HappensAfterGraphVclock::operator==(const HappensAfterGraphVclock& oth) const {
  if (mapping.size() != oth.mapping.size())
    return false;

  for (auto& it : *this) {
    const Node *oth_nd = oth.getNode(it.getEvent());
    if (!oth_nd  || (it.successors != oth_nd->successors))
      return false;
  }

  return true;
}

HappensAfterGraphVclock::Node *
HappensAfterGraphVclock::mapAnnotToNode(const AnnotationKeyT& iid) const
{
  if (iid.instruction == nullptr) {
    assert(initial_node);
    return initial_node;
  }

  // find the read and write nodes for this annotation pair
  auto it = instr_to_node.find(iid.instruction);
  if (it == instr_to_node.end())
    return nullptr;

  for (Node *nd : it->second) {
    if (nd->getEvent()->order == iid.order &&
        nd->getEvent()->cpid == iid.cpid)
      return nd;
  }

  return nullptr;
}

// get events (along with reduced annotation made from these events)
// for which we do not know the order in the graph
std::pair<PositiveAnnotation, std::set<const DCEvent *>>
HappensAfterGraphVclock::getUnorderedEvents() const {
  PositiveAnnotation annot;
  std::set<const DCEvent *> new_nodes;

  // map annotation to nodes beforehead, because otherwise we will try to look
  // them up  it for every node
  std::vector<std::pair<Node *, Node *>> annot_nodes;
  annot_nodes.reserve(annotation.size());

  for (auto& rw_pair : annotation) {
    Node *write = mapAnnotToNode(rw_pair.second);
    assert(write);
    // we ignore initial annotation, because
    // if the write is the initial write, then we added all the edges
    // directly to the graph in addNeceassaryEdges, so do not consider the
    // initial write here
    if (!write->getEvent())
      continue;

    Node *read = mapAnnotToNode(rw_pair.first);
    assert(read);
    annot_nodes.emplace_back(read, write);
  }

  for (auto& it : *this) {
    const DCEvent *ev = it.getEvent();
    // no ev means the initial node which we want to skip
    if (!ev)
      continue;

    if (!ev->may_conflict || !isWrite(*ev))
      continue;

    // found a write, so check whether it is conflicting
    // with some read and write in the annotation
    for (auto& rw_pair : annot_nodes) {
      Node *write = rw_pair.second;
      assert(write->getEvent());

      if (write->getEvent() == ev)
        continue;

      // is the event conflicting with the annotation?
      if (!ev->ml.overlaps(write->getEvent()->ml))
        continue;

      Node *read = rw_pair.first;
      if (read->getEvent() == ev)
        continue;

      // here we need the transitive closure -- so that a path in a graph
      // is also an edge and we can do this check.
      if (!hasConnected(ev, read->getEvent())
          || !hasConnected(write->getEvent(), ev)) {
        // this is our guy
        new_nodes.insert(ev);
        new_nodes.insert(read->getEvent());
        new_nodes.insert(write->getEvent());
        annot.add(*rw_pair.first->getEvent(), *rw_pair.second->getEvent());
      }
    }
  }

  return {annot, new_nodes};
}

Basis HappensAfterGraphVclock::createSubbasis(std::set<const DCEvent *>& new_nodes) const {
  // we have a set of nodes that will be in our graph,
  // so now create the induced subraph, along with the subbasis
  // and annotation
  Basis::ProcessesT events;
  events.reserve(basis.size());
  // create the initial process
  events.emplace_back();
  for (auto& proc : basis) {
    // create a new process when we added something to the
    // previous one (otherwise reuse it)
    if (!events.back().empty())
      events.emplace_back();

    for (const DCEvent *ev : proc) {
      if (new_nodes.count(ev) > 0)
        events.back().push_back(ev);
    }
  }

  // remove the empty processes
  while (!events.empty() && events.back().empty())
    events.pop_back();

  return Basis(std::move(events));
}

// create a subgraph induced by the undecided edges
// we suppose that the graph is transitively closed
HappensAfterGraphVclock HappensAfterGraphVclock::createSubgraph(const PositiveAnnotation& annot,
                                                                const Basis& new_basis) const {
  // Now we can create the subgraph with the new annotation and new basis
  // - happens_before is going to be empty
  HappensAfterGraphVclock subg(annot, {}, new_basis);

  // the only thing that left is to add the induced edges by these nodes
  // FIXME: this could be done more smartly
  for (Node& nd1 : subg) {
    for (Node& nd2 : subg) {
      if (&nd1 == &nd2)
        continue;

      // if the original graph has a path between the nodes add the edge
      // to new nodes too
      // (because of the transitive closure we can check only for the edge)
      if (hasEdge(nd1.getEvent(), nd2.getEvent()))
        nd1.addSuccessor(&nd2);
    }
  }

  return subg;
}

bool HappensAfterGraphVclock::addTransitivity(Events2SAT& sat, const DCEvent *e1, const DCEvent *e2)
{
  Node *n1 = getNode(e1);
  Node *n2 = getNode(e2);
  assert(n1);
  assert(n2);

  // the procss of the node 'n1' and the index
  auto& source_prcs = basis[n1->process_id];
  unsigned idx = n1->event_id;

  const auto& succs = n2->getSuccessors(); // where am I now?
  // processes to explore
  std::set<unsigned> to_process;
  for (int i = 0, e = basis.size(); i < e; ++i) {
    if (succs[i] != INT_MAX) {
      to_process.insert(i);
    }
  }

  // we will take every event x upward from n1 and try to find conflicting events
  // reachable from n2 -- if there will be a thread where we do not find nothing,
  // then shift x one event up and repeat
  while (idx-- > 0) {
    Node *source_nd = getNode(source_prcs[idx]);
    assert(source_nd);

    // check whether we are lucky already with the first node (n2),
    // if so, directly add the edge and bail out, since we found
    // something that covers all threads
    if (source_nd->getEvent()->ml.overlaps(e2->ml)) {
      if (!addImplication(sat, e1, e2, source_nd->getEvent(), e2))
        return false;

      return true;
    }

    for (int i : to_process) {
      std::set<unsigned> to_erase;
      for (int idx2 = succs[i], e = basis[i].size(); idx2 < e; ++idx2) {
        // current node
        Node *target_nd = getNode(basis[i][idx2]);
        assert(target_nd);

        // check that the node I am currently at is conflicting or not
        if (target_nd->getEvent()->ml.overlaps(source_nd->getEvent()->ml)) {
          if (!addImplication(sat, e1, e2,
                              source_nd->getEvent(), target_nd->getEvent()))
            return false;

          // stop processing this thread
          to_erase.insert(i);
        }
      }

      for (int x : to_erase)
        to_process.erase(x);

      if (to_process.empty())
        return true;
    }
  }

  return true;
}

bool HappensAfterGraphVclock::addAnnotationClauses(Events2SAT& sat) {
  // map annotation to nodes beforehead, because otherwise we will try to look
  // them up  it for every node
  std::vector<std::pair<Node *, Node *>> annot_nodes;
  annot_nodes.reserve(annotation.size());

  for (auto& rw_pair : annotation) {
    Node *write = mapAnnotToNode(rw_pair.second);
    assert(write);
    // we ignore initial annotation, because
    // if the write is the initial write, then we added all the edges
    // directly to the graph in addNeceassaryEdges, so do not consider the
    // initial write here
    if (!write->getEvent())
      continue;

    Node *read = mapAnnotToNode(rw_pair.first);
    assert(read);
    annot_nodes.emplace_back(read, write);
  }

  // iterate over all nodes in the graph and find writes
  for (auto& it : *this) {
    const DCEvent *ev = it.getEvent();

    if (!ev || !ev->may_conflict)
      continue;

    if (!isWrite(*ev))
      continue;

    // add the annotation clauses
    for (auto& rw_pair : annot_nodes) {
      Node *write = rw_pair.second;

      if (ev == write->getEvent())
        continue;

      // the event is conflicting?
      if (!ev->ml.overlaps(write->getEvent()->ml))
        continue;

      Node *read = rw_pair.first;

      // this clause do not need to be added for the initial event,
      // since there the ev->read and ev->write edges are already there
      if (!addImplication(sat, ev, read->getEvent(), ev, write->getEvent()))
        return false;

      if (!addImplication(sat, write->getEvent(), ev, read->getEvent(), ev))
        return false;
    }
  }

  return true;
}

bool HappensAfterGraphVclock::createFormula(Events2SAT& sat)
{
  std::vector<std::pair<Node *, Node *>> Vc_confl;

  Vc_confl.reserve(size());
  // reserve some space for the formula beforehead,
  // so that we get rid of many reallocations
  sat.reserve(size());

  const Basis::ProcessT *root = &basis.getTopologyRoot();

  // gather the events that are in our topology
  // - those are the events about which's order
  // we need to decide
  for (const DCEvent *ev1 : *root) {
    // for every two processes where one of them
    // is the root or that are the same
    Node *nd1 = getNode(ev1);
    assert(nd1);

    for (auto& process2 : basis) {
      if (&process2 == root)
        continue;

      for (const DCEvent *ev2 : process2) {
        assert(ev2);
        assert(ev1 != ev2);

        Node *nd2 = getNode(ev2);
        assert(nd2);
        if (ev1->ml.overlaps(ev2->ml) && !hasConnected(nd1, nd2)) {
          // are the event conflicting and we do not know their order yet?
          Vc_confl.push_back({nd1, nd2});
        }
      }
    }
  }

  // annotation clauses
  if (!addAnnotationClauses(sat))
    return false;

  // add transitivity and fact clausees
  for (auto& it : Vc_confl) {
    const DCEvent *e1 = it.first->getEvent();
    const DCEvent *e2 = it.second->getEvent();

    if (!addTransitivity(sat, e1, e2))
      return false;
  }
  // the Vc_confl contains only one ward, we need to fill
  // it into a symmetric relation
  for (auto& it : Vc_confl) {
    const DCEvent *e2 = it.first->getEvent();
    const DCEvent *e1 = it.second->getEvent();

    if (!addTransitivity(sat, e1, e2))
      return false;
  }

  return true;
}

bool HappensAfterGraphVclock::realize()
{
  // here we also check whether the graph is acyclic
  if (!addNecessaryEdges()) {
    blocked_before_sat = true;
    return false;
  }

  auto tmp = getUnorderedEvents();
  PositiveAnnotation& annot = tmp.first;
  std::set<const DCEvent *>& events = tmp.second;

  if (events.empty()) {
    assert(annot.empty());
    // we're done! no unordered events
    assert(isAcyclic() && "Bug in formula or 2SAT");
    return true;
  }

  Basis new_basis = createSubbasis(events);

  // set the basis' root
  int new_tr = new_basis.cpidToProc(basis.getTopologyRootCPid());
  if (new_tr == -1) {
    // we do not have any event in the topology root, which means
    // that we should have everything decided already
    return true;
  } else
    new_basis.setTopologyRoot(new_tr);

  HappensAfterGraphVclock subg = createSubgraph(annot, new_basis);
  if (subg.size() == 0) {
    // we have nothing to decide, bail out
    return true;
  }

  // create the 2SAT formula
  Events2SAT sat;
  if (!subg.createFormula(sat))
    return false;

  // solve the 2SAT formula
  if (!sat.solve())
    return false;

  // the formula is satisfiable,
  // add the edges that we computed
  for (const auto& var : sat) {
    int valuation = std::get<2>(var);
    assert(valuation != -1 && "2SAT: unvaluated variable");

    const DCEvent *e1 = std::get<0>(var);
    const DCEvent *e2 = std::get<1>(var);

#ifdef NDEBUG
    if ((e1 && e1->cpid == basis.getTopologyRootCPid())
        || (e2 && e2->cpid == basis.getTopologyRootCPid()))
        assert((sat.getValuation(e2, e1) == false) && "Order is not antisymmetric");
#endif

    Node *n1 = getNode(e1);
    Node *n2 = getNode(e2);
    assert(n1);
    assert(n2);

    if (valuation == 1)
      n1->addSuccessor(n2);
    else
      // if there is not an edge between these two,
      // then there is the opposite edge
      n2->addSuccessor(n1);
  }

  // after we added all the edges when the
  // 2SAT formula is satisfiable, the graph
  // must be acyclic
  assert(isAcyclic() && "Bug in formula or 2SAT");

  return true;
}

bool HappensAfterGraphVclock::makeTransitiveClosure()
{
  auto order = computeTopoOrder();
  if (order.first == false)
    return false;

  assert(order.second.size() == size());
  for (unsigned i = 0, e = size(); i < e; ++i) {
    Node *cur = order.second[i];
    // we need a copy of the original successor,
    // so that we know the immediate successors even
    // after some modification of the successors
    auto orig_succs = cur->successors;

    // iterate over immediate successors
    for (unsigned idx = 0, end_idx = orig_succs.size();
         idx < end_idx; ++idx) {

      // no successor in the thread?
      if (orig_succs[idx] == INT_MAX)
        continue;

      // get the immediate successor of this node
      assert(idx < processes.size());
      assert(idx < orig_succs.size());
      assert((unsigned)orig_succs[idx] < processes[idx].size());

      Node *succ = processes[idx][orig_succs[idx]];

      // check where we can get from the immediate successor
      // and adjust our successors
      assert(succ->successors.size() == cur->successors.size());
      for (unsigned idx2 = 0, end2 = succ->successors.size();
           idx2 < end2; ++idx2) {
        cur->successors[idx2]
          = std::min(cur->successors[idx2], succ->successors[idx2]);
      }
    }
  }

  return true;
}

// add edges that are implied by the annotation
bool HappensAfterGraphVclock::addNecessaryEdgesFromAnnotation() {
  for (auto& it : annotation) {
    Node *r = mapAnnotToNode(it.first);
    Node *w = mapAnnotToNode(it.second);
    assert(r);
    assert(w);

    if (r->process_id == w->process_id) {
      // we checked whether such annotation is valid
      // before (the 'w' must be the closest one in
      // this thread)
      continue;
    }

    if (w->process_id == INT_MAX) {
      // if r should see the initial write, then r must be before
      // every other confl. write in every thread
      // (it is enough to add it to the first confl. event in the thread,
      // the rest is implied by the transitivity)
      for (auto& proc : processes) {
        for (Node *nd : proc) {
          const DCEvent *ev = nd->getEvent();
          if (isWrite(*ev) && ev->ml.overlaps(r->getEvent()->ml)) {
            r->addSuccessor(nd);
            break;
          }
        }
      }
    } else {
      assert(w->event_id != INT_MAX);
      // we already have the edge w->r,
      // we can directly add edges x->w where
      // x is the first confl. write above r
      // and r->y where y is the first confl.
      // write after w in its thread
      unsigned j = r->event_id;
      while (j > 0) {
        const DCEvent *ev = processes[r->process_id][j - 1]->getEvent();
        assert(ev);
        if (isWrite(*ev) && ev->ml.overlaps(w->getEvent()->ml)) {
          processes[r->process_id][j - 1]->addSuccessor(w);
          break;
        }
        --j;
      }

      j = w->event_id + 1;
      for (unsigned end = processes[w->process_id].size(); j < end; ++j) {
        const DCEvent *ev = processes[w->process_id][j]->getEvent();
        if (isWrite(*ev) && ev->ml.overlaps(w->getEvent()->ml)) {
          r->addSuccessor(processes[w->process_id][j]);
          break;
        }
      }
    }
  }

  return true;
}

// add edges that are implied by the annotation
bool HappensAfterGraphVclock::addNecessaryEdgesFromAnnotationClauses() {
  // map annotation to nodes beforehead, because otherwise we will try to look
  // them up  it for every node
  std::vector<std::pair<Node *, Node *>> annot_nodes;
  annot_nodes.reserve(annotation.size());

  for (auto& rw_pair : annotation) {
    Node *write = mapAnnotToNode(rw_pair.second);
    assert(write);
    // we ignore initial annotation, because
    // if the write is the initial write, then we added all the edges
    // directly to the graph in addNeceassaryEdges, so do not consider the
    // initial write here
    if (!write->getEvent())
      continue;

    Node *read = mapAnnotToNode(rw_pair.first);
    assert(read);
    annot_nodes.emplace_back(read, write);
  }

  for (auto& it : *this) {
    const DCEvent *ev = it.getEvent();
    if (!ev || !ev->may_conflict)
      continue;

     if (!isWrite(*ev))
       continue;

     // add the annotation clauses
     for (auto& rw_pair : annot_nodes) {
       Node *write = rw_pair.second;
       assert(write->getEvent());

       if (!write->getEvent())
          continue;
       // no ev means the initial node, which belong to the
       // initial event that is always before anything,
       // so we do not consider it
       if (ev == write->getEvent())
         continue;

      // the event is conflicting?
      if (!ev->ml.overlaps(write->getEvent()->ml))
        continue;

      Node *read = rw_pair.first;
      Node *evnd = getNode(ev);
      assert(evnd);
      if (hasEdge(ev, read->getEvent()))
          evnd->addSuccessor(write);
      if (hasEdge(write->getEvent(), ev))
        read->addSuccessor(evnd);
     }
   }

  return true;
}

// add edges that are implied by the happens-exactly-before relation
bool HappensAfterGraphVclock::addNecessaryEdgesFromHB()
{
  for (auto&it : happens_before) {
    Node *ev = mapAnnotToNode(it.first);
    assert(ev);
    assert(ev != initial_node);

    for (auto& it2 : it.second) {
      Node *see_ev
        = mapAnnotToNode({it2.first, it2.second.first, it2.second.second});
      // see ev must be before ev
      assert(see_ev != ev);
      see_ev->addSuccessor(ev);

      // it2 is {CPid, {instruction, order}}
      unsigned int i, j;
      if (see_ev == initial_node) {
          j = 0;
          // find the right process according to CPid
          i = 0;
          for (auto& proc : processes) {
            if (proc[0]->getEvent()->cpid == it2.first)
              break;
            ++i;
          }
          // we must have found the process
          assert(i != processes.size());
      } else {
          assert(see_ev->process_id != INT_MAX && see_ev->event_id != INT_MAX);
          i = see_ev->process_id;
          j = see_ev->event_id + 1; // move to the next event in the thread
      }
      // find the first conflicting ev after see_ev
      // in its thread
      for (unsigned e = processes[i].size(); j < e; ++j) {
        Node *nd = processes[i][j];
        // this can not be the initial ev
        assert(nd->getEvent());
        const DCEvent *ndev = nd->getEvent();

        // check this because of pthread locks,
        // for which mutex_init is conflicting with lock/unlock
        if (!isRead(*ndev) && !isWrite(*ndev))
          continue;

        if (ndev->ml.overlaps(ev->getEvent()->ml)) {
            assert(ev != nd);
            ev->addSuccessor(nd);
            // found, done
            break;
        }
      }
    }
  }

  return true;
}

bool HappensAfterGraphVclock::addNecessaryEdges()
{
  if (!addNecessaryEdgesFromAnnotation())
    return false;

  if (!addNecessaryEdgesFromHB())
    return false;

  // NOTE: We need to do a transitive closure before trying to
  // create the formula, because with the happens-before relation
  // it is no longer true that 2SAT can decide everything.
  // The 2SAT will only work on happens-after graph that is
  // transitively closed and that has the edges we added
  // at the beggining of this method. Moreover, while doing
  // the transitive closure, we compute the topological order
  // which tells us whether the graph is already acyclic or not
  if (!makeTransitiveClosure())
    return false;

  // when we have the transitive closure,
  // we can check for some edges implied by the annotation
  if (!addNecessaryEdgesFromAnnotationClauses())
    return false;

  // we need to make the transitive closure again...
  // (but on acyclic graph we have it quadratic)
  if (!makeTransitiveClosure())
    return false;

  // FIXME: the transitive closure now introduced new relations
  // that can allow us to decide more edges. These edges
  // should be added to 2SAT as fact clauses, since we got rid
  // of fact clauses. The example is:
  //
  //  w'---|
  //  |    |
  //  |    v
  //  |    R
  //  w<---|
  //
  //  if we do not have decided the edge wR or w'R, then
  //  annotation clauses will add variable Xww', which is
  //  not restricted in any other way (fact clauses are not there,
  //  so it may decide that w->w')
  //
  //  NOTE: is this still true??

  return true;
}

// this works since the graph is acyclic
std::vector<DCEvent> HappensAfterGraphVclock::linearize()
{
  std::vector<DCEvent> trace;
  auto order = computeTopoOrder();
  assert(order.first && "The graph is not acyclic");

  // the order is reversed
  assert(order.second.size() == size());
  for (unsigned i = order.second.size(); i > 0; --i)
    if (order.second[i - 1]->getEvent()) {
      trace.emplace_back(order.second[i - 1]->getEvent()->blank_copy());
      trace.back().id = trace.size() - 1;
    } else
      // the initial event must be the first event
      assert(i == size());

  return trace;
}

void HappensAfterGraphVclock::_constructGraph()
{
  // create nodes and starting points
  assert(!basis.empty());

  // set of nodes that call pthread_create and pthread_join
  std::vector<Node *> spawns;
  std::vector<Node *> joins;
  spawns.reserve(4);
  joins.reserve(4);

  unsigned process_idx = 0;
  for (auto& base : basis) {
    Node *node;
    Node *last = nullptr;

    // create a new process in nodes
    //  - nodes are basically a copy of basis,
    //  in this moment, but with nodes instead of DCEvents
    processes.emplace_back();

    // create the nodes and add the program structure
    unsigned event_idx = 0;
    for (const DCEvent *event : base) {
      assert(event->iid.get_pid() >= 0);
      node = new Node(basis, process_idx, event_idx);
      mapping.emplace(event, node);
      processes.back().push_back(node);

      if (event->instruction)
        instr_to_node[event->instruction].push_back(node);

      if (last)
        last->addSuccessor(node);

      last = node;

      // XXX: what about function pointer calls?
      if (event->instruction) {
        if (is_function_call(event->instruction, "pthread_create"))
          spawns.push_back(node);
        else if (is_function_call(event->instruction, "pthread_join"))
          joins.push_back(node);
      }

      ++event_idx;
    }

    ++process_idx;
  }

  // fill in the annotation edges
  for (auto& rw_pair : annotation) {
    Node *read = mapAnnotToNode(rw_pair.first);
    Node *write = mapAnnotToNode(rw_pair.second);

    assert(read);
    if (!write)
      write = initial_node;

    write->addSuccessor(read);
  }

  // fill in the threads creation happens-before relation
  for (Node *spwn : spawns) {
    // find the first event of the thread in basis
    for (auto& base : basis) {
      assert(!base.empty());
      if (base[0]->cpid == spwn->getEvent()->childs_cpid) {
        spwn->addSuccessor(getNode(base[0]));
        break;
      }
    }
  }

  for (Node *jn : joins) {
    // find the first event of the thread in basis
    for (auto& base : basis) {
      assert(!base.empty());
      if (base[0]->cpid == jn->getEvent()->childs_cpid) {
        getNode(base[base.size() - 1])->addSuccessor(jn);
        break;
      }
    }
  }
}

void HappensAfterGraphVclock::to_dot(const char *edge_params) const
{
  llvm::errs() << "digraph {\n";
  for (auto& it : *this) {
    llvm::errs() << "NODE" << &it << " [label=\"";
    if (!it.getEvent())
      llvm::errs() << "init";
    else
      llvm::errs() << *it.getEvent();

    llvm::errs() << "\"]\n";

  }

  llvm::errs() << "\n";
  for (auto& it : *this) {
    for (unsigned idx = 0; idx < it.successors.size(); ++idx) {
      unsigned eid = it.successors[idx];
      // no edge to this thread
      if (eid == INT_MAX)
        continue;
      auto succ = processes[idx][eid];
      llvm::errs() << "NODE" << &it << " -> NODE" << succ;
      if (edge_params) {
        llvm::errs() << "[" << edge_params << "]\n";
      } else {
        llvm::errs() << "\n";
      }
    }
  }

  llvm::errs() << "}\n";
}

void HappensAfterGraphVclock::dump() const
{
  llvm::errs() << " --- annotation:\n";
  annotation.dump();
  llvm::errs() << " --- HB:\n";
  happens_before.dump();
  llvm::errs() << " --- processes:\n";
  unsigned idx = 0;
  llvm::errs() << "[";
  for (int x : initial_node->successors) {
    if (x == INT_MAX)
      llvm::errs() << "x, ";
    else
      llvm::errs() << x << ", ";
  }
  llvm::errs() << "] ";
  llvm::errs() << " -- init\n";

  for (auto& proc : processes) {
    llvm::errs() << "process " << idx++ << "\n";
    unsigned id = 0;
    for (const Node *nd : proc) {
      llvm::errs() << id++ << " [";
      for (int x : nd->successors) {
        if (x == INT_MAX)
          llvm::errs() << "x, ";
        else
          llvm::errs() << x << ", ";
      }
      llvm::errs() << "] ";
      llvm::errs() << *nd->getEvent() << "\n";
   }
  }
  llvm::errs() << " -----------\n";
}

/* DEBUGGING
#include <sstream>
#include <llvm/Support/raw_os_ostream.h>
std::string HappensAfterGraphVclock::to_string() const
{
  std::ostringstream sstream;
  llvm::raw_os_ostream stream(sstream);

  stream << " --- annotation:\n";
  annotation.dump();
  stream << " --- HB:\n";
  happens_before.dump();
  stream << " --- processes:\n";
  unsigned idx = 0;
  stream << "[";
  for (int x : initial_node->successors) {
    if (x == INT_MAX)
      stream << "x, ";
    else
      stream << x << ", ";
  }
  stream << "] ";
  stream << " -- init\n";

  for (auto& proc : processes) {
    stream << "process " << idx++ << "\n";
    unsigned id = 0;
    for (const Node *nd : proc) {
      stream << id++ << " [";
      for (int x : nd->successors) {
        if (x == INT_MAX)
          stream << "x, ";
        else
          stream << x << ", ";
      }
      stream << "] ";
      stream << *nd->getEvent() << "\n";
   }
  }
  stream << " -----------\n";

  stream.flush();
  return sstream.str();
}
*/

 std::pair<bool, std::vector<HappensAfterGraphVclock::Node *>>
 HappensAfterGraphVclock::computeTopoOrder()
{
  struct info {
    unsigned dfsid = 0, lowpt = 0;
    bool is_on_stack = false;
  };

  // initialize the info about nodes
  info init_info; // info of the initial node
  std::vector<std::vector<info>> infos;
  infos.resize(processes.size());
  unsigned i = 0;
  for (auto& inf : infos)
    inf.resize(processes[i++].size());

  std::function<info&(Node *)> get_info = [&infos, &init_info](Node *nd)-> info& {
    if (nd->process_id == INT_MAX)
      return init_info;
    return infos[nd->process_id][nd->event_id];
  };

  std::stack<Node *> stack;
  std::vector<Node *> ret;
  unsigned index = 0;

  std::function<bool(Node *n)> _compute = [&](Node *n) -> bool {
    ++index;
    auto& inf = get_info(n);
    inf.dfsid = index;
    inf.lowpt = index;
    stack.push(n);
    inf.is_on_stack = true;

    for (auto it = succ_begin(n), end = succ_end(n); it != end; ++it) {
      Node *succ = *it;
      auto& succ_inf = get_info(succ);
      if (succ_inf.dfsid == 0) {
        assert(succ_inf.is_on_stack == false);
        if (!_compute(succ))
          return false;
        inf.lowpt = std::min(inf.lowpt, succ_inf.lowpt);
      } else if (succ_inf.is_on_stack) {
        return false; // cycle
      }
    }

    if (inf.lowpt == inf.dfsid) {
        int its = 0;
        while (get_info(stack.top()).dfsid >= inf.dfsid) {
          ++its;
          Node *w = stack.top();
          stack.pop();

          get_info(w).is_on_stack = false;
          ret.push_back(w);

          if (stack.empty())
              break;
        }
        assert(its == 1);
    }

    return true;
  };

  if (!_compute(initial_node))
    return {false, {}};

  assert(stack.empty());
  return {true, ret};
}

