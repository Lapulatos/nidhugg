#ifndef _HAPPENS_EXACTLY_BEFORE_H_
#define _HAPPENS_EXACTLY_BEFORE_H_

#include <map>
#include <vector>

#include "Annotations.h"

// a mapping from Instructions to subset of Instructions
// XXX factor out the common code with annotations to a template
class HappensExactlyBefore {
public:
  // a mapping Events -> (Threads -> Instruction)
  typedef std::pair<const llvm::Instruction *, unsigned> InstrIdT;
  typedef std::map<CPid, InstrIdT> ValueT;
  typedef std::map<AnnotationKeyT, ValueT> MappingT;
  typedef MappingT::iterator iterator;
  typedef MappingT::const_iterator const_iterator;

  const MappingT& getMapping() const { return mapping; }

  HappensExactlyBefore() = default;

  HappensExactlyBefore(const HappensExactlyBefore& a)
  : mapping(a.mapping) {}

  HappensExactlyBefore& operator=(const HappensExactlyBefore& a)
  {
    mapping = a.mapping;
    return *this;
  }

  HappensExactlyBefore(HappensExactlyBefore&& a)
  : mapping(std::move(a.mapping)) {}

  HappensExactlyBefore& operator=(HappensExactlyBefore&& a)
  {
    mapping = std::move(a.mapping);
    return *this;
  }

  const ValueT *get(const AnnotationKeyT& k) const {
    auto it = mapping.find(k);
    if (it == mapping.end())
      return nullptr;

    return &it->second;
  }

  const ValueT *get(const DCEvent& ev) const {
    return get(AnnotationKeyT(ev.cpid, ev.instruction, ev.order));
  }

  iterator begin() { return mapping.begin(); }
  const_iterator begin() const { return mapping.begin(); }
  iterator end() { return mapping.end(); }
  const_iterator end() const { return mapping.end(); }

  bool empty() const { return mapping.empty(); }

  void add(const AnnotationKeyT& a, const AnnotationValueT& b)
  {
    assert(a != b);
    mapping[a].emplace(b.cpid, InstrIdT(b.instruction, b.order));
  }

  void add(const DCEvent& a, const CPid& cpid,
           const llvm::Instruction *I, unsigned ord)
  {
    assert(a.cpid != cpid || a.instruction != I);
    mapping[{a.cpid, a.instruction, a.order}].emplace(cpid, InstrIdT(I, ord));
  }

  void add(const DCEvent& a, const DCEvent& b)
  {
    add(AnnotationKeyT(a.cpid, a.instruction, a.order),
        AnnotationValueT(b.cpid, b.instruction, b.order));
  }

  void dump() const;

private:
  MappingT mapping;
};

#endif // _HAPPENS_EXACTLY_BEFORE_H_
