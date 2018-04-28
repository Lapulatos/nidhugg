#include <unordered_map>

#include "DCTBHelpers.h" // isWrite
#include "HappensAfterGraph.h"
#include "2SAT.h"
#include "SCC.h"

HappensAfterGraph::~HappensAfterGraph()
{
  for (auto& it : nodes)
    delete it.second;
}

HappensAfterGraph::HappensAfterGraph(const HappensAfterGraph& oth,
                                     const HappensExactlyBefore& hb)
: annotation(oth.annotation),
  //happens_before(oth.happens_before),
  happens_before(hb),
  basis(oth.basis),
  max_search_index(0),
  blocked_before_sat(false)
{
  nodes.reserve(oth.nodes.size());
  instr_to_node.reserve(instr_to_node.size());

  // create new copies of nodes and update mappings
  for (auto& it : oth.nodes) {
    Node *new_nd = new Node(it.first);

    if (it.first && it.first->instruction)
      instr_to_node[it.first->instruction].push_back(new_nd);

    nodes.emplace(it.first, new_nd);
  }

  // map successors
  for (auto& it : oth.nodes) {
    Node *new_nd = nodes[it.first];
    assert(new_nd);

    for (Node *succ : it.second->successors) {
      Node *new_succ = nodes[succ->event];
      assert(new_succ);
      new_nd->addSuccessor(new_succ);
    }
  }

  // also set the initial_node
  initial_node = nodes[nullptr];
  assert(initial_node);
}

HappensAfterGraph::HappensAfterGraph(HappensAfterGraph&& oth)
: annotation(std::move(oth.annotation)),
  happens_before(std::move(oth.happens_before)),
  basis(std::move(oth.basis)),
  max_search_index(oth.max_search_index),
  blocked_before_sat(oth.blocked_before_sat),
  nodes(std::move(oth.nodes)),
  instr_to_node(std::move(oth.instr_to_node)),
  initial_node(oth.initial_node)
{
}

HappensAfterGraph& HappensAfterGraph::operator=(HappensAfterGraph&& oth) {
  if (&oth != this)
    *this = std::move(oth);
  return *this;
}

bool HappensAfterGraph::isAcyclic()
{
  AcyclicTopoOrder<Node> detect(max_search_index,
                                true /* only detect cycles */);
  auto start = begin();
  auto en = end();
  bool ret = detect.compute(start, en);
  max_search_index = detect.getIndex();

  return ret;
}

HappensAfterGraph::Node *
HappensAfterGraph::mapAnnotToNode(const AnnotationKeyT& iid) const
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
    if (nd->event->order == iid.order &&
        nd->event->cpid == iid.cpid)
      return nd;
  }

  return nullptr;
}

std::set<HappensAfterGraph::Node *>&
HappensAfterGraph::getNodesReachableFromNoLazy(Node *n)
{
  std::vector<Node *> to_process;
  to_process.push_back(n);
  auto& ret = reachable_nodes[n];
  while(!to_process.empty()) {
    std::vector<Node *> newly_discovered;
    newly_discovered.reserve(to_process.size());

    for (Node *nd : to_process) {
      if (ret.insert(nd).second) {
        for (Node *succ : nd->successors)
          newly_discovered.push_back(succ);
      }
    }

    to_process.swap(newly_discovered);
  }

  // we do not want this to be reflexive
  ret.erase(n);
  return ret;
}

std::set<HappensAfterGraph::Node *>&
HappensAfterGraph::getNodesReachingToNoLazy(Node *n)
{
  std::vector<Node *> to_process;
  to_process.push_back(n);
  auto& ret = rev_reachable_nodes[n];
  while(!to_process.empty()) {
    std::vector<Node *> newly_discovered;
    newly_discovered.reserve(to_process.size());

    for (Node *nd : to_process) {
      if (ret.insert(nd).second) {
        for (Node *pred : nd->predecessors)
          newly_discovered.push_back(pred);
      }
    }

    to_process.swap(newly_discovered);
  }

  // we do not want this to be reflexive
  ret.erase(n);
  return ret;
}

// this code suppose that the graph is acyclic, so it does
// not do a fixpoint (since on acyclic graph the fixpoint
// is reached after one iteration).
void HappensAfterGraph::makeTransitiveClosureSimple()
{
    for (auto& it : nodes) {
        Node *nd = it.second;
        for (Node *s : getNodesReachableFrom(nd))
          nd->addSuccessor(s);
    }
}

void HappensAfterGraph::makeTransitiveClosure()
{
  AcyclicTopoOrder<Node> topo(max_search_index);

  auto start = begin();
  auto en = end();
  bool ret = topo.compute(start, en);
  assert(ret && "The graph is not acyclic");

  // XXX: maybe we should do just simple recursive DFS?
  max_search_index = topo.getIndex();
  auto& ord = topo.getOrder();
  assert(ord.size() == size());
  for (unsigned i = 0, e = ord.size(); i < e; ++i) {
    // make a copy of current successors
    auto succs = ord[i]->getSuccessors();
    for (Node *succ : succs)
      for (Node *succ_of_succ : succ->getSuccessors())
          ord[i]->addSuccessor(succ_of_succ);
  }

  assert(isTransitivelyClosed());
}

bool HappensAfterGraph::isTransitivelyClosed()
{
  for (Node& nd1 : *this) {
    for (Node *succ1 : nd1.getSuccessors()) {
        for (Node *nd2 : succ1->getSuccessors()) {
            // we have edge nd1->succ1 && succ1->nd2
            // so check whether we have nd1->nd2
            if (!hasEdge(&nd1, nd2))
              return false;
        }
    }
  }

  return true;
}

bool HappensAfterGraph::addNecessaryEdges()
{
  // projection of reads and writes from basis
  std::vector<std::vector<Node *>> readswrites;
  // indexes of nodes in the readswrites
  std::unordered_map<Node *, std::pair<unsigned, unsigned>> indexes;
  std::vector<std::pair<Node *, Node *>> new_edges;
  std::map<CPid, unsigned> cpid_to_threads;

  readswrites.reserve((nodes.size() / 2) * 3);
  indexes.reserve((nodes.size() / 2) * 3);
  new_edges.reserve(nodes.size() / 5);

  unsigned i = 0;
  for (auto& base : basis) {
    readswrites.emplace_back();
    cpid_to_threads[base[0]->cpid] = i;

    unsigned j = 0;
    for (const DCEvent *e : base) {
      // make projection to only writes and reads
      if (!isRead(*e) && !isWrite(*e))
        continue;

      Node *n = nodes[e];
      assert(n && "Do not have the node");
      indexes[n] = std::make_pair(i, j);
      readswrites[i].push_back(n);
      ++j;
    }
    ++i;
  }

  for (auto& proc : readswrites) {
    for (Node *n : proc) {
      // if a node has more successors, then
      // it means there's the annotation edge
      // to the other thread, so it must be a read
      // that is in the image of annotation
      if (n->successors.size() < 2)
        continue;

      assert(isWrite(*n->event));
      // for all successors of the write except the
      // one in this process:
      for (Node *read : n->successors) {
        if (read->event->cpid == n->event->cpid)
          continue;

        assert(isRead(*read->event));
        // find the closest conflicting write in this process
        // that happend before this read and add new edge
        // from that write to 'n' (which is the observation write)
        auto it = indexes.find(read);
        assert(it != indexes.end());
        for (int j = it->second.second - 1; j >= 0; --j) {
          Node *confl_w = readswrites[it->second.first][j];
          if (!confl_w->event->instruction ||
              !isWrite(*confl_w->event))
            continue;

          if (confl_w->event->ml.overlaps(read->event->ml)) {
            assert(confl_w != n);
            new_edges.emplace_back(confl_w, n);
            break;
          }
        }

        // find the closest write after 'n' and add an edge
        // from read to it
        it = indexes.find(n);
        assert(it != indexes.end());
        for (int j = it->second.second + 1,
             e = readswrites[it->second.first].size(); j < e ; ++j) {
          Node *confl_w = readswrites[it->second.first][j];
          if (!confl_w->event->instruction ||
              !isWrite(*confl_w->event))
            continue;

          if (confl_w->event->ml.overlaps(read->event->ml)) {
            new_edges.emplace_back(read, confl_w);
            break;
          }
        }
      }
    }
  }

  for (auto&it : happens_before) {
    Node *ev = mapAnnotToNode(it.first);
    assert(ev);
    assert(ev != initial_node);

    for (auto& it2 : it.second) {
        Node *seen_ev
          = mapAnnotToNode({it2.first, it2.second.first, it2.second.second});
        // seen ev must be before ev
        assert(seen_ev != ev);
        new_edges.push_back({seen_ev, ev});

        // it2 is {CPid, instruction}
        int i, j;
        if (seen_ev == initial_node) {
            j = 0;
            i = cpid_to_threads[it2.first];
        } else {
            assert(indexes.count(seen_ev) > 0);
            std::tie(i, j) = indexes[seen_ev];
            ++j; // move to the next event in the thread
        }
        // find the first conflicting ev after see_ev
        // in its thread
        for (int e = readswrites[i].size(); j < e; ++j) {
            Node *nd = readswrites[i][j];
            // this can not be the initial ev
            assert(nd->event);

            if (nd->event->ml.overlaps(ev->event->ml)) {
                assert(ev != nd);
                new_edges.push_back({ev, nd});
                // found, done
                break;
            }
        }
    }
  }

  for (auto& pr : new_edges)
    pr.first->addSuccessor(pr.second);

  // check for acyclicity, since makeTransitiveClosure()
  // works on acyclic graphs (it does not do a fixpoint)
  if (!isAcyclic())
    return false;

  // NOTE: We need to do a transitive closure before trying to
  // create the formula, because with the happens-before relation
  // it is no longer true that 2SAT can decide everything.
  // The 2SAT will only work on happens-after graph that is
  // transitively closed and that has the edges we added
  // at the beggining of this method.
  makeTransitiveClosure();

  return true;
}

// this works since the graph is acyclic
std::vector<DCEvent> HappensAfterGraph::linearize()
{
  std::vector<DCEvent> trace;
  AcyclicTopoOrder<Node> detect(max_search_index);

  auto start = begin();
  auto en = end();
  bool ret = detect.compute(start, en);
  assert(ret && "The graph is not acyclic");

  max_search_index = detect.getIndex();
  // the order is reversed
  auto& ord = detect.getOrder();
  for (unsigned i = ord.size(); i > 0; --i)
    if (ord[i - 1]->event) {
      trace.push_back(ord[i - 1]->event->blank_copy());
      trace.back().id = trace.size() - 1;
    } else
      // the initial event must be the first event
      assert(i == ord.size());

  return trace;
}

// create a new basis that is restricted to the events in the annotation and the
// past cones from the annotation
Basis HappensAfterGraph::getRestrictedBasis(const PositiveAnnotation& annot) const
{
  // a mapping from pid to maximal length of a process
  std::map<CPid, std::pair<unsigned, Node *>> mapping;

  // find the furthest event in every thread that
  // is annotated (and also the node with this event)
  for (const auto& rw : annot) {
    auto& r = mapping[rw.first.cpid];
    Node *rn = mapAnnotToNode(rw.first);
    if (rn->event_id + 1 > r.first) {
      r.second = rn;
      assert(r.second);
      r.first = rn->event_id + 1;
    }

    auto& w = mapping[rw.second.cpid];
    Node *wn = mapAnnotToNode(rw.second);
    if (wn->event_id + 1 > w.first) {
      w.second = wn;
      assert(r.second);
      w.first = wn->event_id + 1;
    }
  }

  // now do the same with past-cones of all these events
  // BFS flags
  std::vector<bool> visited;
  visited.resize(nodes.size(), false);

  std::vector<Node *> to_process;
  to_process.reserve(mapping.size());
  for (auto& it: mapping) {
    Node *cur = it.second.second;
    if (!cur)
      continue;
    visited[cur->id] = true;
    to_process.emplace_back(cur);
  }

  while (!to_process.empty()) {
    std::vector<Node *> new_to_process;
    new_to_process.reserve(to_process.size());
    for (Node *cur : to_process) {
      for (Node *pred : cur->getPredecessors()) {
        if (!visited[pred->id]) {
          // put the node into the queue
          new_to_process.emplace_back(pred);
          visited[pred->id] = true;

          if (!pred->event)
            continue;

          auto& m = mapping[pred->event->cpid];
          // did we discover something new
          if (m.first < pred->event_id + 1) {
            m.first = pred->event_id + 1;
          }
        }
      }
    }

    to_process.swap(new_to_process);
  }

  // create the new basis
  Basis::ProcessesT ret;
  for (auto& process : basis) {
    auto& m = mapping[process[0]->cpid];
    unsigned ev_num = m.first;

    if (ev_num != 0) {
      Basis::ProcessT proc;
      proc.reserve(process.size());
      unsigned int idx = 0;
      for (const DCEvent *ev : process) {
        proc.emplace_back(ev);
        if (++idx == ev_num)
          break;
      }
      ret.emplace_back(std::move(proc));
    }
  }

  return ret;
}

static inline const DCIID evToDCIID(const DCEvent& ev) {
    return DCIID(ev.cpid, ev.instruction, ev.order);
}


PositiveAnnotation
HappensAfterGraph::getPastConeAnnotation(const PositiveAnnotation& A,
                                         const std::initializer_list<const DCEvent *>& l) const
{
  std::vector<const Node *> to_process;
  to_process.reserve(std::max(l.size(), (size_t)8));

  for (const DCEvent *ev : l) {
    if (!ev->instruction) // initial event
      continue;

    const Node *nd = getNode(ev);
    assert(nd);
    to_process.emplace_back(nd);
  }

  return getPastConeAnnotation(A, std::move(to_process));
}

PositiveAnnotation
HappensAfterGraph::getPastConeAnnotation(const PositiveAnnotation& A,
                                         const std::initializer_list<const DCIID *>& l) const
{
  std::vector<const Node *> to_process;
  to_process.reserve(std::max(l.size(), (size_t)8));

  for (const DCIID *ev : l) {
    if (!ev->instruction) // initial event
      continue;

    const Node *nd = mapAnnotToNode(*ev);
    assert(nd);
    to_process.emplace_back(nd);
  }

  return getPastConeAnnotation(A, std::move(to_process));
}

PositiveAnnotation
HappensAfterGraph::getPastConeAnnotation(const PositiveAnnotation& A,
                                         std::vector<const Node *> &&tp) const
{
  PositiveAnnotation ret;
  std::vector<const Node *> to_process = tp;

  // BFS flags
  std::vector<bool> visited;
  visited.resize(nodes.size(), false);

  // mark the reads from given annotation as visited,
  // since the past-cone of such reads is already in the annotation
  // and we do not need to look further
  for (const auto& rw_pair : A) {
    Node *nd = mapAnnotToNode(rw_pair.first);
    assert(nd);
    visited[nd->id] = true;
  }

  while (!to_process.empty()) {
    std::vector<const Node *> new_to_process;
    new_to_process.reserve(to_process.size());

    for (const Node *cur : to_process) {
      for (const Node *pred : cur->getPredecessors()) {
        if (!visited[pred->id]) {
          // put the node into the queue
          new_to_process.emplace_back(pred);
          visited[pred->id] = true;

          if (!pred->event)
            continue;

          if (isRead(*pred->event)) {
            // annotation (the class attribute) is the current
            // observation function
            auto obs = annotation.get(evToDCIID(*pred->event));
            assert(obs);
            ret.add(evToDCIID(*pred->event), *obs);
          }
        }
      }
    }

    to_process.swap(new_to_process);
  }

  return ret;
}

bool HappensAfterGraph::canReach(const DCEvent *ev1, const DCEvent *ev2) const
{
  const Node *ev1_node = getNode(ev1);
  const Node *ev2_node = getNode(ev2);
  assert(ev1_node);
  assert(ev2_node);

  std::vector<const Node *> to_process = { ev1_node };

  // BFS flags
  std::vector<bool> visited;
  visited.resize(nodes.size(), false);
  visited[ev1_node->id] = true;

  while (!to_process.empty()) {
    std::vector<const Node *> new_to_process;
    new_to_process.reserve(to_process.size());

    for (const Node *cur : to_process) {
      for (const Node *succ : cur->getSuccessors()) {
        // the graph is acyclic
        assert(succ != ev1_node);

        if (succ == ev2_node)
          return true;

        if (!visited[succ->id]) {
          // put the node into the queue
          new_to_process.emplace_back(succ);
          visited[succ->id] = true;
        }
      }
    }

    to_process.swap(new_to_process);
  }

  // did not find the ev2_node
  return false;
}

void HappensAfterGraph::addAnnotation(const PositiveAnnotation& annotation)
{
  // fill in the annotation edges
  for (auto& rw_pair : annotation) {
    Node *read = mapAnnotToNode(rw_pair.first);
    Node *write = mapAnnotToNode(rw_pair.second);

    assert(read);
    if (!write)
      write = initial_node;

    write->addSuccessor(read);
  }
}

void HappensAfterGraph::_constructGraph()
{
  // this should be more efficient than to dynamically adjust the hash table
  unsigned reserve_size = 0;
  for (auto&b : basis)
    reserve_size += b.size();
  instr_to_node.reserve(reserve_size);
  nodes.reserve(reserve_size);
  reachable_nodes.reserve(reserve_size);
  rev_reachable_nodes.reserve(reserve_size);

  // create the initial node
  initial_node = new Node(nullptr);
  addNode(initial_node);

  // create nodes and starting points
  assert(!basis.empty());
  // set of nodes that call pthread_create and pthread_join
  std::vector<Node *> spawns;
  std::vector<Node *> joins;

  unsigned idx = 0;
  for (auto& base : basis) {
    Node *node;
    Node *last = nullptr;

    // create the nodes and add the program structure
    unsigned idx2 = 0;
    for (const DCEvent *event : base) {
      assert(event->iid.get_pid() >= 0);
      node = new Node(event);
      node->process_id = idx;
      node->event_id = idx2;
      addNode(node);
      if (event->instruction)
        instr_to_node[event->instruction].push_back(node);

      if (last)
        last->addSuccessor(node);
      else
        // make the first event in the base
        // a successor of the inital node
        initial_node->addSuccessor(node);

      last = node;

      // XXX: what about function pointer calls?
      if (event->instruction) {
        if (is_function_call(event->instruction, "pthread_create"))
          spawns.push_back(node);
        else if (is_function_call(event->instruction, "pthread_join"))
          joins.push_back(node);
      }
      ++idx2;
    }
    ++idx;
  }

  addAnnotation(annotation);

  // fill in the threads creation happens-before relation
  for (Node *spwn : spawns) {
    // find the first event of the thread in basis
    for (auto& base : basis) {
      assert(!base.empty());
      if (base[0]->cpid == spwn->event->childs_cpid) {
        spwn->addSuccessor(nodes[base[0]]);
        break;
      }
    }
  }

  for (Node *jn : joins) {
    // find the first event of the thread in basis
    for (auto& base : basis) {
      assert(!base.empty());
      if (base[0]->cpid == jn->event->childs_cpid) {
        nodes[base[base.size() - 1]]->addSuccessor(jn);
        break;
      }
    }
  }
}

void HappensAfterGraph::to_dot(const char *edge_params) const
{
  llvm::errs() << "digraph {\n";
  for (auto& it : nodes) {
    llvm::errs() << "NODE" << it.second << " [label=\"";
    if (!it.second->event)
      llvm::errs() << "init";
    else
      llvm::errs() << *it.second->event;

    llvm::errs() << "\"]\n";

  }

  llvm::errs() << "\n";
  for (auto& it : nodes) {
    // succ is a pointer
    for (auto succ : it.second->successors) {
      llvm::errs() << "NODE" << it.second << " -> NODE" << succ;
      if (edge_params) {
        llvm::errs() << "[" << edge_params << "]\n";
      } else {
        llvm::errs() << "\n";
      }
    }
  }

  llvm::errs() << "}\n";
}

void HappensAfterGraph::dump() const
{
  llvm::errs() << " --- annotation:\n";
  annotation.dump();
  llvm::errs() << " --- basis:\n";
  unsigned idx = 0;
  for (auto& base : basis) {
    llvm::errs() << "base " << idx++ << "\n";
    for (const DCEvent *e : base)
      llvm::errs() << *e << "\n";
  }
}

