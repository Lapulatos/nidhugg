// this file was taken (but modified) from
// https://github.com/mchalupa/dg/blob/master/src/analysis/SCC.h


#ifndef _SCC_H_
#define  _SCC_H_

#include <stack>
#include <vector>
#include <set>

// implementation of tarjan's algorithm for
// computing strongly connected components
// for a directed graph that has a starting vertex
// from which are all other vertices reachable
template <typename NodeT>
class SCC {
public:
    typedef std::vector<unsigned> SCC_component_t;
    typedef std::vector<SCC_component_t> SCC_t;

    SCC<NodeT>(std::vector<NodeT>& nds, unsigned not_visit = 0)
    : nodes(nds), index(not_visit), NOT_VISITED(not_visit) {}

    // returns a vector of vectors - every inner vector
    // contains the nodes that for a SCC
    SCC_t& compute()
    {
      for (unsigned i = 0; i < nodes.size(); ++i) {
        if (not_visited(nodes[i]))
          _compute(i);
      }

        assert(stack.empty());

        return scc;
    }

    SCC_component_t& operator[](unsigned idx)
    {
        assert(idx < scc.size());
        return scc[idx];
    }

    unsigned getIndex() const { return index; }

private:
    std::vector<NodeT>& nodes;
    std::stack<unsigned> stack;
    unsigned index;
    // it the dfsid is less or equal to this value,
    // then it is considered not to be visited.
    // it is due to the repeated execution of this algorithm
    // on the same nodes
    const unsigned NOT_VISITED;

    // container for the strongly connected components.
    SCC_t scc;

    bool not_visited(NodeT& nd) { return nd.get_dfsid() <= NOT_VISITED; }

    void _compute(unsigned i)
    {
        NodeT& n = nodes[i];
        ++index;
        n.set_dfsid(index);
        n.set_lowpt(index);
        stack.push(i);
        n.set_on_stack(true);

        for (unsigned i : n.getSuccessors()) {
            NodeT& succ = nodes[i];

            if (not_visited(succ)) {
                assert(!succ.on_stack());
                _compute(i);
                n.set_lowpt(std::min(n.get_lowpt(), succ.get_lowpt()));
            } else if (succ.on_stack()) {
                n.set_lowpt(std::min(n.get_lowpt(), succ.get_dfsid()));
            }
        }

        if (n.get_lowpt() == n.get_dfsid()) {
            SCC_component_t component;
            size_t component_num = scc.size();

            unsigned ndfsid = n.get_dfsid();
            while (nodes[stack.top()].get_dfsid() >= ndfsid) {
                unsigned widx = stack.top();
                NodeT& w = nodes[widx];
                stack.pop();

                w.set_on_stack(false);
                component.emplace_back(widx);
                // the numbers scc_id give
                // a reverse topological order
                w.set_scc_id(component_num);

                if (stack.empty())
                    break;
            }

            scc.emplace_back(std::move(component));
        }
    }
};

// this class goes through the graph and detects cycles,
// if there are no cycles, it returns vertices in
// reverse topological order
template <typename NodeT>
class AcyclicTopoOrder {
public:
    AcyclicTopoOrder(unsigned not_visit = 0, bool only_cycle = false)
    : index(not_visit), NOT_VISITED(not_visit), only_detect_cycle(only_cycle) {}

    template <typename IT>
    bool compute(IT& start, IT& end)
    {
      for (; start != end; ++start) {
        if (not_visited(*start))
          if (!_compute(*start))
            return false;
      }

        assert(stack.empty());
        return true;
    }

    const std::vector<NodeT *>& getOrder() const { return scc; }

    unsigned operator[](unsigned idx)
    {
        assert(idx < scc.size());
        return scc[idx];
    }

    unsigned getIndex() const { return index; }

private:
    std::stack<NodeT *> stack;
    unsigned index;
    // it the dfsid is less or equal to this value,
    // then it is considered not to be visited.
    // it is due to the repeated execution of this algorithm
    // on the same nodes
    const unsigned NOT_VISITED;
    bool only_detect_cycle;

    // container for the strongly connected components.
    std::vector<NodeT *> scc;

    bool not_visited(NodeT& nd) { return nd.get_dfsid() <= NOT_VISITED; }

    bool _compute(NodeT& n)
    {
        ++index;
        n.set_dfsid(index);
        n.set_lowpt(index);
        stack.push(&n);
        n.set_on_stack(true);

        for (NodeT *s : n.getSuccessors()) {
            NodeT& succ = *s;

            if (not_visited(succ)) {
                assert(!succ.on_stack());
                if (!_compute(succ))
                  return false;
                n.set_lowpt(std::min(n.get_lowpt(), succ.get_lowpt()));
            } else if (succ.on_stack()) {
              // cycle
              return false;
            }
        }

        if (n.get_lowpt() == n.get_dfsid()) {
            unsigned ndfsid = n.get_dfsid();
            int its = 0;
            while (stack.top()->get_dfsid() >= ndfsid) {
              ++its;
              NodeT *w = stack.top();
              stack.pop();

              w->set_on_stack(false);
              if (!only_detect_cycle)
                scc.push_back(w);

              if (stack.empty())
                  break;
            }
            assert(its == 1);
        }

        return true;
    }
};



#endif //  _SCC_H_
