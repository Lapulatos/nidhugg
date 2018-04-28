#ifndef _2SAT_H_
#define _2SAT_H_

#include <unordered_map>
#include <set>

#include "Event.h"

// this class internally holds the implication graph
// and assignments to other variables
class Events2SAT {
public:
  class Node {
    // successors in the implication graph
    //std::set<unsigned> _implies;
    std::vector<unsigned> _implies;
    //llvm::SmallSet<unsigned, 8> _implies;

    // for the computation of strongly connected components
    // XXX: we could keep that in external arrays...
    unsigned dfsid, lowpt, scc_id;
    bool is_on_stack;

    friend class Events2SAT;

  public:
    // this node represents an edge
    // from e1 to e2 (that is, event e1
    // happens before e2)
    const DCEvent *e1, *e2;

    Node(const DCEvent *a, const DCEvent *b)
    // scc_id is undefined after initialization
    : dfsid(0), lowpt(0), is_on_stack(0),
      e1(a), e2(b) {
        assert(e1 != e2); /* do not want self-loops */
        _implies.reserve(8);
      }

    void set_on_stack(bool a) { is_on_stack = a; }
    void set_dfsid(unsigned id) { dfsid = id; }
    void set_lowpt(unsigned lp) { lowpt = lp; }
    void set_scc_id(unsigned id) { scc_id = id; }

    unsigned get_dfsid() const { return dfsid; }
    unsigned get_lowpt() const { return lowpt; }
    unsigned get_scc_id() const { return scc_id; }
    bool on_stack() const { return is_on_stack; }

    //const std::set<unsigned> getSuccessors() const { return _implies; }
    const std::vector<unsigned>& getSuccessors() const { return _implies; }

    void implies(unsigned idx) {
        // suppose we have only few successors,
        // so it is faster to iterate over sequential
        // array, than over big set
        /* OR: we just ignore multiple-edges,
         * there's not much of them (if any)
        for (unsigned i : _implies) {
          if (i == idx)
            return;
        }
        */

        _implies.push_back(idx);
    }
  };

  struct Info {
    int index = -1;
    int valuation = -1;

    Info() : index(-1), valuation(-1) {}
    Info(int idx, int val = -1) : index(idx), valuation(val) {}
  };

  typedef std::unordered_map<const DCEvent *,
                             std::unordered_map<const DCEvent *, Info>> InfosT;

  void reserve(size_t size) {
      infos.reserve(size);
      nodes.reserve(2*size);
  }

  // add an assignment of an edge e1->e2,
  // return false if the formula becomes provably
  // unsatisfable by adding this valuation
  bool addAssignment(const DCEvent *e1, const DCEvent *e2, bool val)
  {
    // a = a v a = !a -> a
    // NOTE: maybe we could just set a to 'val' and add the
    // edge !a->a (or a->!a for val == false)
    if (!addImplication(e1, e2, e1, e2, !val, val))
      return false;

    Info& info1 = infos[e1][e2];
    assert(info1.valuation == -1 || info1.valuation == val);

    info1.valuation = val;
    return true;
  }

  // (e11,e12) => (e21, e22)
  // val1 and val2 describe whether the variables are negated
  // (true means non negated and false means negated)
  bool addImplication(const DCEvent *e11, const DCEvent *e12,
                      const DCEvent *e21, const DCEvent *e22,
                      bool val1 = true, bool val2 = true)
  {
    Info& info1 = infos[e11][e12];
    Info& info2 = infos[e21][e22];

    // check whether we do not have an unsatisfiable implication
    if (info1.valuation > -1) {
      bool v1 = val1 ? info1.valuation : !info1.valuation;

      if (info2.valuation > -1) {
        bool v2 = val2 ? info2.valuation : !info2.valuation;
        // true => false ?
        if (v1 && !v2)
          return false;
        else
          // we already have the nodes, so we do not need to add
          // them and can return from here
          return true;
      } else if (v1) {
        // we may also direcly set an assignment when we have the
        // valuation of the premise set to true
        info2.valuation = val2 ? 1 : 0;
        // we do not need to add a node into the graph
        // when we derived the value of the second node
        return true;
      }
    } else {
      // we have valuation for implication, but not for the premise?
      // in that case check wether the implication is false,
      // because that means that the premise must be false to in order
      // to satisfy the implication
      if (info2.valuation > -1) {
        bool v2 = val2 ? info2.valuation : !info2.valuation;
        if (!v2) {
          info1.valuation = val1 ? 0 : 1;
          return true;
        }
      }
    }

    // OK, add the corresponding nodes if we do not have them
    if (info1.index == -1) {
      assert(nodes.size() % 2 == 0);

      info1.index = nodes.size();
      nodes.push_back(Node(e11, e12));
      nodes.push_back(Node(e11, e12));
    }

    if (info2.index == -1) {
      assert(nodes.size() % 2 == 0);

      info2.index = nodes.size();
      nodes.push_back(Node(e21, e22));
      nodes.push_back(Node(e21, e22));
    }

    assert(info1.index != -1);
    assert(info2.index != -1);

    // Add the edges to the implications graph.
    // If we have implication (a -> b), it is equivalent
    // to a CNF formula (!a v b), so we must add edges
    // a -> b && !b -> !a
    // Positive literals are on even indices
    // and the negative are on odd indices.
    unsigned idx1 = info1.index;
    if (!val1)
      ++idx1;

    unsigned idx2 = info2.index;
    if (!val2)
      ++idx2;

    assert(idx1 < nodes.size());
    assert(idx2 < nodes.size());
    // a -> b
    nodes[idx1].implies(idx2);

    idx1 = info2.index;
    if (val2) // negate
      ++idx1;

    idx2 = info1.index;
    if (val1) // negate
      ++idx2;

    assert(idx1 < nodes.size());
    assert(idx2 < nodes.size());
    // !b -> !a
    nodes[idx1].implies(idx2);

    return true;
  }

  bool getValuation(const DCEvent *e1, const DCEvent *e2)
  {
    Info& inf = infos[e1][e2];
    assert(inf.valuation != -1);

    return (bool) inf.valuation;
  }

  bool solve();

  const Node *getNode(unsigned idx) const { return &nodes[idx]; }
  const Info* getInfo(const DCEvent *e1, const DCEvent *e2) const
  {
    auto it = infos.find(e1);
    if (it == std::end(infos))
      return nullptr;

    auto it2 = it->second.find(e2);
    if (it2 == std::end(it->second))
      return nullptr;

    return &it2->second;
  }

  const InfosT& getInfos() const { return infos; }

  // this allows to iterate over all variables for which
  // we have some an assignment
  class const_iterator {
    const InfosT* _cont;
    InfosT::const_iterator it;
    InfosT::mapped_type::const_iterator sub_it;

    const_iterator(const InfosT *cont, bool begin)
      :_cont(cont)
    {
      if (begin) {
        it = cont->begin();
        sub_it = it->second.begin();
      } else
        it = cont->end();
    }

    friend class Events2SAT;

  protected:
    bool at_end() const { return it == _cont->end(); }
    void set_next()
    {
      ++sub_it;
      if (sub_it == it->second.end()) {
        ++it;
        if (!at_end())
          sub_it = it->second.begin();
      }
    }

  public:
    // empty
    const_iterator() {};
    // copy
    const_iterator(const const_iterator& oth)
      :_cont(oth._cont), it(oth.it), sub_it(oth.sub_it) {}

    const_iterator& operator++()
    {
      set_next();
      return *this;
    }

    const_iterator operator++(int)
    {
      auto tmp = *this;
      set_next();
      return tmp;
    }

    std::tuple<const DCEvent *, const DCEvent *, int> operator*() const
    {
      return std::tuple<const DCEvent *, const DCEvent *, int>(it->first,
                                                           sub_it->first,
                                                           sub_it->second.valuation);
    }

    bool operator==(const const_iterator& oth) const
    {
      assert(_cont);
      return oth.it == it && (at_end() || oth.sub_it == sub_it);
    }

    bool operator!=(const const_iterator& oth) const
    {
      return !operator==(oth);
    }

  };

  const_iterator begin() const { return const_iterator(&infos, true); }
  const_iterator end() const { return const_iterator(&infos, false); }

  size_t vars_num() const
  {
    size_t sz = 0;
    for (auto& it : infos)
      sz += it.second.size();

    return sz;
  }

  void to_dot() const;

private:
  // An edge (e1, e2) corresponds to a pair of variables
  // X and X' such that X' = not X. Xs are on even indices,
  // X's are on odd indices
  std::vector<Node> nodes;

  // already computed assignments and nodes that we have
  // infos[e1][e2] = (index to nodes, valuation)
  InfosT infos;
};

#endif // _2SAT_H_
