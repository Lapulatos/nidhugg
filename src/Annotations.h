#ifndef _ANNOTATIONS_H_
#define _ANNOTATIONS_H_

#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>

#include "IID.h"
#include "Event.h"

class Basis;

struct DCIID {
    CPid cpid;
    const llvm::Instruction *instruction = nullptr;
    unsigned order = 0;

    DCIID() = default;
    DCIID(const CPid& pid, const llvm::Instruction *i, unsigned ord)
      : cpid(pid), instruction(i), order(ord) {}

    DCIID(const std::tuple<const CPid&, const llvm::Instruction *, unsigned>& tup)
      : cpid(std::get<0>(tup)),
        instruction(std::get<1>(tup)),
        order(std::get<2>(tup)) {}

    DCIID(const DCIID& oth)
      : cpid(oth.cpid),
        instruction(oth.instruction),
        order(oth.order) {}

    DCIID& operator=(const DCIID&) = default;

    DCIID(DCIID&& oth)
      : cpid(std::move(oth.cpid)),
        instruction(oth.instruction),
        order(oth.order) {
    }

    DCIID& operator=(DCIID&& oth) {
        cpid = std::move(oth.cpid);
        instruction = oth.instruction;
        order = oth.order;
        return *this;
    }

    // compare the cpid as the last one, because it is the most expensive
    bool operator<(const DCIID& rhs) const {
        return std::tie(instruction, order, cpid)
            < std::tie(rhs.instruction, rhs.order, rhs.cpid);
    }

    bool operator==(const DCIID& rhs) const {
        return std::tie(instruction, order, cpid)
            == std::tie(rhs.instruction, rhs.order, rhs.cpid);
    }

    bool operator!=(const DCIID& rhs) const {
        return !operator==(rhs);
    }

    bool operator==(const DCEvent& ev) const {
        return std::tie(instruction, order, cpid)
            == std::tie(ev.instruction, ev.order, ev.cpid);
    }

    bool operator!=(const DCEvent& ev) const {
        return !operator==(ev);
    }

    void dump() const;
};

typedef DCIID AnnotationKeyT;
typedef DCIID AnnotationValueT;

// a mapping from Instructions to Instructions
class PositiveAnnotation {
public:
  typedef std::map<AnnotationKeyT, AnnotationValueT> MappingT;
  typedef MappingT::iterator iterator;
  typedef MappingT::const_iterator const_iterator;

  const MappingT& getMapping() const { return mapping; }

  PositiveAnnotation() = default;
  PositiveAnnotation(const PositiveAnnotation&) = default;
  PositiveAnnotation& operator=(const PositiveAnnotation&) = default;
  PositiveAnnotation(PositiveAnnotation&& a) = default;
  PositiveAnnotation& operator=(PositiveAnnotation&& a) = default;

  bool operator==(const PositiveAnnotation& oth) const {
      return mapping == oth.mapping;
  }

  void intersect(const PositiveAnnotation& rhs) {
    MappingT new_mapping;
    auto it1 = mapping.begin();
    auto it2 = rhs.mapping.begin();
    auto last1 = mapping.end();
    auto last2 = rhs.mapping.end();

    while (it1 != last1 && it2 != last2) {
        if (*it1 < *it2) {
            ++it1;
        } else  {
            if (!(*it2 < *it1)) {
                new_mapping.emplace(*it1++);
            }
            ++it2;
        }
    }

    new_mapping.swap(mapping);
  }

  size_t size() const { return mapping.size(); }

  const AnnotationValueT *get(const AnnotationKeyT& k) const
  {
    auto it = mapping.find(k);
    if (it == mapping.end())
      return nullptr;

    return &it->second;
  }

  const AnnotationValueT *get(const DCEvent& e) const
  {
    return get(AnnotationKeyT(e.cpid, e.instruction, e.order));
  }

  // add a annotation (a, b). If there is already an annotation
  // for a, return false.
  bool add(const AnnotationKeyT& a, const AnnotationValueT& b)
  {
    auto it = mapping.find(a);
    if (it != mapping.end())
      return false;

    mapping.emplace_hint(it, a, b);
    return true;
  }

  bool add(AnnotationKeyT&& a, AnnotationValueT&& b)
  {
    auto it = mapping.find(a);
    if (it != mapping.end())
      return false;

    mapping.emplace_hint(it, a, b);
    return true;
  }

  void erase(const AnnotationKeyT& k) {
    mapping.erase(k);
  }

  bool add(const DCEvent& a, const DCEvent& b)
  {
    // b may not have instruction when it is the initial event
    assert(a.instruction && "Read does not have an instruction");
    return add(AnnotationKeyT(a.cpid, a.instruction, a.order),
               AnnotationValueT(b.cpid, b.instruction, b.order));
  }

  bool defines(const AnnotationKeyT& k) const
  {
    return mapping.find(k) != mapping.end();
  }

  bool defines(const DCEvent& ev) const
  {
    return defines(ev.cpid, ev.instruction, ev.order);
  }

  bool defines(const CPid& c, const llvm::Instruction *i, unsigned ord) const
  {
    return defines(AnnotationKeyT(c, i, ord));
  }

  bool empty() const { return mapping.empty(); }

  iterator begin() { return mapping.begin(); }
  const_iterator begin() const { return mapping.begin(); }
  iterator end() { return mapping.end(); }
  const_iterator end() const { return mapping.end(); }

  void dump() const;

  static PositiveAnnotation getObservationFunction(const std::vector<DCEvent>& trace);
  // returns a positive annotation for an event that forces threads to come
  // to this event
  // \param ev     the event in question
  // \param O      observation function
  // \param trace  trace from which to deduce the past code
  static PositiveAnnotation getPastConeAnnotation(const DCEvent& ev,
                                                  const PositiveAnnotation& O,
                                                  Basis& basis);

private:
  MappingT mapping;
};

// a mapping from Instructions to subset of Instructions
class NegativeAnnotation {
public:
  typedef std::pair<AnnotationKeyT, AnnotationValueT> AnnotPair;
  typedef std::set<AnnotPair> AnnotPairSet;
  typedef std::map<AnnotationValueT, AnnotPairSet> AnnotPairSetMap;
  // the negative annotation is a mapping
  // Event -> (Event -> {Event})
  typedef std::map<AnnotationKeyT, AnnotPairSetMap> MappingT;
  typedef MappingT::iterator iterator;
  typedef MappingT::const_iterator const_iterator;

  const MappingT& getMapping() const { return mapping; }

  NegativeAnnotation() = default;

  NegativeAnnotation(const NegativeAnnotation& a)
  : mapping(a.mapping) {}

  NegativeAnnotation& operator=(const NegativeAnnotation& a)
  {
    mapping = a.mapping;
    return *this;
  }

  NegativeAnnotation(NegativeAnnotation&& a)
  : mapping(std::move(a.mapping)) {}

  NegativeAnnotation& operator=(NegativeAnnotation&& a)
  {
    mapping = std::move(a.mapping);
    return *this;
  }

  iterator begin() { return mapping.begin(); }
  const_iterator begin() const { return mapping.begin(); }
  iterator end() { return mapping.end(); }
  const_iterator end() const { return mapping.end(); }

  bool forbids(const PositiveAnnotation& positive_annotation,
               const AnnotationKeyT& r, const AnnotationValueT& w) const
  {
    auto it = mapping.find(r);
    if (it == mapping.end())
      return false;

    auto it2 = it->second.find(w);
    if (it2 == it->second.end())
        return false;

    // find out wheter we have annotated any read in the extra negative annotation set
    bool found_annotated_r = false;
    for (const auto& annot : it2->second) {
      auto rann = positive_annotation.get(annot.first);
      if (rann != nullptr) {
        found_annotated_r = true;
        // if we have annotated read and the read sees
        // the same, we forbid this mutation
        if (*rann == annot.second) {
            return true;
        }
      }
    }

    // we have annotated read, but it sees something different?
    // Then we do not forbid this mutation
    if (found_annotated_r)
      return false;

    return true;
  }

  bool forbids(const PositiveAnnotation& positive_annotation,
               const DCEvent& r, const DCEvent& w) const
  {
    return forbids(positive_annotation,
                    AnnotationKeyT(r.cpid, r.instruction, r.order),
                    AnnotationValueT(w.cpid, w.instruction, w.order));
  }

  bool forbids(const PositiveAnnotation& positive_annotation,
               const PositiveAnnotation& past_cone) const
  {
    // return true if there are some events in the
    // past cone and all pairs from the past-cone
    // are forbidden by the negative annotation
    if (past_cone.empty())
      return false;

    for (auto& rw_pair : past_cone) {
      if (!forbids(positive_annotation, rw_pair.first, rw_pair.second))
        return false;
    }
    return true;
  }

  bool forbidsSome(const PositiveAnnotation& positive_annotation,
                   const PositiveAnnotation& past_cone) const
  {
    // return true if there are some events in the
    // past cone and all pairs from the past-cone
    // are forbidden by the negative annotation
    if (past_cone.empty())
      return false;

    for (auto& rw_pair : past_cone) {
      if (forbids(positive_annotation, rw_pair.first, rw_pair.second))
        return true;
    }
    return false;
  }

  bool empty() const { return mapping.empty(); }

  void add(const AnnotationKeyT& a,
           const AnnotationValueT& b,
           const AnnotPairSet& extraAnnots) {
    auto& S = mapping[a][b];
    S.insert(extraAnnots.begin(), extraAnnots.end());
  }

  void add(const AnnotationKeyT& a,
           const AnnotationValueT& b,
           const PositiveAnnotation& extraAnnots) {
    auto& S = mapping[a][b];
    S.insert(extraAnnots.begin(), extraAnnots.end());
  }

  void add(const AnnotationKeyT& a,
           const AnnotationValueT& b,
           PositiveAnnotation&& extraAnnots) {
    auto& S = mapping[a][b];
    S.insert(extraAnnots.begin(), extraAnnots.end());
  }

  void add(const DCEvent& a, const DCEvent& b,
           const AnnotPairSet& extraAnnots) {
    add(AnnotationKeyT(a.cpid, a.instruction, a.order),
        AnnotationValueT(b.cpid, b.instruction, b.order),
        extraAnnots);
  }

  void add(const DCEvent& a, const DCEvent& b,
           const PositiveAnnotation& extraAnnots) {
    add(AnnotationKeyT(a.cpid, a.instruction, a.order),
        AnnotationValueT(b.cpid, b.instruction, b.order),
        extraAnnots);
  }


  void dump() const;

private:
  MappingT mapping;
};

#include "HappensExactlyBefore.h"

struct AnnotTree {
    struct Node {
        unsigned id = 0;
        const char *msg;
        const PositiveAnnotation annotation;
        const HappensExactlyBefore hb;
        const NegativeAnnotation neg_annot;
        std::vector<unsigned> successors;

        Node() : annotation(), hb() {}
        Node(const PositiveAnnotation& a, unsigned i = 0, const char *m = nullptr)
            : id(i), msg(m), annotation(a), hb() {}

        Node(const PositiveAnnotation& a, const HappensExactlyBefore& h,
             unsigned i = 0, const char *m = nullptr)
            : id(i), msg(m), annotation(a), hb(h) {}

        Node(const PositiveAnnotation& a, const NegativeAnnotation& na,
             unsigned i = 0, const char *m = nullptr)
            : id(i), msg(m), annotation(a), hb(), neg_annot(na) {}

        Node(const PositiveAnnotation& a, const NegativeAnnotation& na,
             const HappensExactlyBefore& h, unsigned i = 0, const char *m = nullptr)
            : id(i), msg(m), annotation(a), hb(h), neg_annot(na) {}
    };

    std::vector<Node> nodes;

    AnnotTree() {nodes.reserve(128);}

    void createRoot() {
        assert(nodes.empty());
        nodes.emplace_back();
    }

    void createRoot(const PositiveAnnotation& ann) {
        assert(nodes.empty());
        nodes.emplace_back(ann, 0, "initial");
    }

    unsigned addNode(unsigned parent_id,
                     const PositiveAnnotation& annot,
                     const HappensExactlyBefore& h,
                     const char *m = nullptr) {
        assert(nodes.size() >= 1 && "Do not have a root");

        nodes.emplace_back(annot, h, nodes.size(), m);
        nodes[parent_id].successors.push_back(nodes.size() - 1);

        assert(nodes.back().id == nodes.size() - 1);

        return nodes.size() - 1;
    }

    unsigned addNode(unsigned parent_id,
                     const PositiveAnnotation& annot,
                     const NegativeAnnotation& na,
                     const char *m = nullptr) {
        assert(nodes.size() >= 1 && "Do not have a root");

        nodes.emplace_back(annot, na, nodes.size(), m);
        nodes[parent_id].successors.push_back(nodes.size() - 1);

        assert(nodes.back().id == nodes.size() - 1);

        return nodes.size() - 1;
    }

    unsigned addNode(unsigned parent_id,
                     const PositiveAnnotation& annot,
                     const NegativeAnnotation& na,
                     const HappensExactlyBefore& h,
                     const char *m = nullptr) {
        assert(nodes.size() >= 1 && "Do not have a root");

        nodes.emplace_back(annot, na, h, nodes.size(), m);
        nodes[parent_id].successors.push_back(nodes.size() - 1);

        assert(nodes.back().id == nodes.size() - 1);

        return nodes.size() - 1;
    }



    void dump() const;
};

#endif // _ANNOTATIONS_H_
