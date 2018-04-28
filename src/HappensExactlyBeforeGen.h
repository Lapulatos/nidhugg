#ifndef _HAPPENS_EXACTLY_BEFORE_GEN_H_
#define _HAPPENS_EXACTLY_BEFORE_GEN_H_

#include <map>
#include <vector>
#include "HappensExactlyBefore.h"
#include "HappensAfterGraph.h"

//#define DUMP_HB_TREE 1

#if DUMP_HB_TREE
struct HBTree {
    struct Node {
        unsigned id = 0;
        const char *msg;
        const HappensExactlyBefore hb;
        std::vector<unsigned> successors;

        Node() = default;
        Node(const HappensExactlyBefore& h, unsigned i = 0, const char *m = nullptr)
            : id(i), msg(m), hb(h) {}
    };

    std::vector<Node> nodes;

    void createRoot() {
        assert(nodes.empty());
        nodes.emplace_back();
    }

    unsigned addNode(unsigned parent_id,
                     const HappensExactlyBefore& h,
                     const char *m = nullptr) {
        assert(nodes.size() >= 1 && "Do not have a root");

        nodes.emplace_back(h, nodes.size(), m);
        nodes[parent_id].successors.push_back(nodes.size() - 1);
        assert(nodes.back().id == nodes.size() - 1);

        return nodes.size() - 1;
    }

    void dump() const;
};
#endif // DUMP_HB_TREE


///
// Generator of the Happens-Exactly-Before relation
//
class HappensExactlyBeforeGen {
  HappensAfterGraph& G;
  Basis& basis;

#if DUMP_HB_TREE
  // for debugging
  HBTree tree;
#endif

  // how many trials we blocked
  unsigned blocked_hb = 0;

  // here we store the generated HB
  std::vector<HappensExactlyBefore> generated_hb;

  Basis::events_iterator
  needsHB(HappensExactlyBefore& happens_before,
          const Basis::events_iterator& start);

#if DUMP_HB_TREE
  void _generate(HappensExactlyBefore& happens_before,
                 const Basis::events_iterator& start,
                 unsigned tree_parent_node);
  void _generate(HappensExactlyBefore& happens_before,
                 const Basis::events_iterator& it1,
                 Basis::events_iterator& it2,
                 unsigned tree_parent_node);
#else
  void _generate(HappensExactlyBefore& happens_before,
                 const Basis::events_iterator& start);
  void _generate(HappensExactlyBefore& happens_before,
                 const Basis::events_iterator& it1,
                 Basis::events_iterator& it2);
#endif

public:
  HappensExactlyBeforeGen(HappensAfterGraph& G,
                          Basis& b)
  : G(G), basis(b) {
#if DUMP_HB_TREE
      // XXX: debugging
      tree.createRoot();
#endif
    }

#if DUMP_HB_TREE
  ~HappensExactlyBeforeGen() {
    tree.dump();
  }
#endif

  bool generate(HappensExactlyBefore& happens_before, const Basis::events_iterator& start) {
#if DUMP_HB_TREE
    _generate(happens_before, start, 0);
#else
    _generate(happens_before, start);
#endif
    return !generated_hb.empty();
  }

  std::vector<HappensExactlyBefore>& generated() { return generated_hb; }
  unsigned blocked() { return blocked_hb; }
};

#endif // _HAPPENS_EXACTLY_BEFORE_H_
