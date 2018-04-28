#include <fstream>

#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/raw_os_ostream.h>
#include <llvm/IR/Instruction.h>
#include <llvm/IR/Instructions.h>

#include "DCTBHelpers.h" // isWrite
#include "Annotations.h"
#include "IID.h"
#include "CPid.h"
#include "Basis.h"

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const DCEvent& ev)
{
  out << ev.cpid << "-" << ev.order;
  if (ev.instruction)
    out << *ev.instruction;
  else
    out << " <null>";

  return out;
}

void DCEvent::dump() const {
  llvm::errs() << *this << "\n";
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const DCIID& iid)
{
  out << iid.cpid << "-" << iid.order;
  if (iid.instruction)
    out << *iid.instruction;
  else
    out << " <null>";

  return out;
}

void DCIID::dump() const {
  llvm::errs() << *this << "\n";
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const PositiveAnnotation& annot) {
  out << "Annotation+ {\n";
  for (auto& pr : annot) {
    out << "("  << pr.first << ",\n"
                << " "  << pr.second << ")\n";
  }
  out << "}\n";

  return out;
}

void PositiveAnnotation::dump() const {
  llvm::errs() << *this;
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const NegativeAnnotation& annot) {
  out << "Annotation- {\n";
  for (auto& pr : annot) {
    out << pr.first << " =>\n  {\n";
    for (auto &it : pr.second) {
      out << "  " << it.first << "\n";
      out << "  past cone {\n";
      for (auto &it3 : it.second) {
        out << "    " << it3.first <<
               ", "   << it3.second << "\n";
      }
      out << "  }\n";
    }
    out << "  }\n";
  }
  out << "}\n";

  return out;
}

void NegativeAnnotation::dump() const
{
  llvm::errs() << *this;
}

PositiveAnnotation
PositiveAnnotation::getObservationFunction(const std::vector<DCEvent>& trace)
{
  PositiveAnnotation annot;

  // maps memory to its definition
  std::map<ConstMRef, const DCEvent *> mapping;

  // go through events in this process. We assume that every
  // write overwrites the whole memory
  for (unsigned idx = 0, trace_size = trace.size(); idx < trace_size; ++idx) {
    const DCEvent& ev = trace[idx];
    if (!ev.instruction || !ev.may_conflict)
      continue;

    if (isWrite(ev)) {
      mapping[ev.ml] = &ev;
    } else if (isRead(ev)) {
      auto wit = mapping.find(ev.ml);
      if (wit == mapping.end()) {
        // no write to this memory means the initial write
        annot.add(ev, DCEvent(IID<IPid>(0, 0)));
      } else {
        annot.add(ev, *wit->second);
      }
    }
  }

  return annot;
}

PositiveAnnotation
PositiveAnnotation::getPastConeAnnotation(const DCEvent& ev,
                                          const PositiveAnnotation& O, //observation func
                                          Basis& basis)
{
  PositiveAnnotation annot; // return value

  // if this is the initial event, it has no past-cone
  if (ev.instruction == nullptr)
    return annot;

  Basis::events_iterator it = basis.getIterator(ev);
  assert(it != basis.events_end());

  // Check whether there is some predecessor event in this thread
  // If not, we're done, because it has no past-cone
  if (it.getEventID() == 0)
    return annot;

  // move to the predecessor event, because we do not want
  // to add an annotation for ev
  --it;
  assert(it != basis.events_end());

  // mapping for the maximal events in a process. Key is a process
  // and the value is the maximum of indexes of processes that we have
  // to traverse
  std::map<unsigned, unsigned> mapping = {{it.getProcessID(), it.getEventID()}};
  std::set<unsigned> processed = {it.getProcessID()};
  while (!mapping.empty()) {
    auto it = mapping.begin();
    mapping.erase(it);

    unsigned proc = it->first;
    // idx = it->second + 1; --> we will go until idx == 1
    // and access using idx - 1
    for (unsigned idx = it->second + 1; idx > 0; --idx) {
      const DCEvent& ev = *basis[proc][idx - 1];

      if (!isRead(ev))
        continue;

      // we found a read, so add the write it observes
      const AnnotationValueT *wr = O.get(ev);
      assert(wr && "Do not have observation function for a read");
      if (!wr->instruction) { // the initial write
        annot.add({ev.cpid, ev.instruction, ev.order}, {wr->cpid, nullptr, 0});
        // nothing more to be done here
        continue;
      }

      assert(wr->instruction);
      auto wr_ev_it = basis.getIterator(wr->cpid, wr->order);
      assert(wr_ev_it != basis.events_end());
      const DCEvent& wr_ev = **wr_ev_it;
      assert(wr_ev.instruction == wr->instruction);
      annot.add(ev, wr_ev);

      // queue this new write for the backward searching.
      // But only if it is in different process, because we
      // are going through this process now. Also do not do anything
      // if we already processed that process, because... we processed
      // that process :)
      if (wr->cpid == ev.cpid || processed.insert(proc).second)
        continue;

      auto wit2 = mapping.find(wr_ev_it.getProcessID());
      if (wit2 == mapping.end())
        mapping.emplace_hint(wit2,
                             wr_ev_it.getProcessID(),
                             wr_ev_it.getEventID());
      else
        // set from where to search the process
        mapping[wr_ev_it.getProcessID()] = std::max(wit2->second,
                                                    wr_ev_it.getEventID());
    }
  }

  return annot;
}


extern llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const HappensExactlyBefore& hb);

void AnnotTree::dump() const {
  // this runs in every DCTraceBuilder
  // (even when getting basis)
  if (nodes.empty())
      return;

  std::ofstream _annot_debug("annotations.dot");
  llvm::raw_os_ostream annot_debug(_annot_debug);

  annot_debug << "digraph Annotations {\n";

  // nodes
  for (auto& nd : nodes) {
    // the node description
    annot_debug << "NODE" << nd.id << " [label=\"";
    annot_debug << nd.id << "\\n";
    annot_debug << nd.annotation;
    if (!nd.hb.empty()) {
      annot_debug << "------------\\n";
      annot_debug << nd.hb;
      annot_debug << "------------\\n";
    }

    if (!nd.neg_annot.empty()) {
      annot_debug << "------------\\n";
      annot_debug << nd.neg_annot;
      annot_debug << "------------\\n";
    }

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

