/*
 * Copyright (C) 2016-2017 Marek Chalupa
 *
 * This file is part of Nidhugg.
 *
 * Nidhugg is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Nidhugg is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <config.h>
#ifndef __DC_TRACE_BUILDER_H__
#define __DC_TRACE_BUILDER_H__

#include <map>
#include <unordered_map>
#include <deque>
#include <llvm/IR/Instructions.h>
#include <llvm/ExecutionEngine/GenericValue.h>

#include <ctime>

#include "TSOTraceBuilder.h"
#include "VClock.h"
#include "Trace.h"
#include "Annotations.h"
#include "AnnotatedTrace.h"
#include "Basis.h"
#include "HappensExactlyBefore.h"
#include "HappensAfterGraphVclock.h"

static inline void not_supported(const char *msg = nullptr) __attribute__((noreturn));
static void not_supported(const char *msg) {
	if (msg)
        llvm::errs() << "Not supported yet: " << *msg << "\n";
	else
        llvm::errs() << "Not supported yet\n";
    abort();
}

//#define DUMP_ANNOT_TREE 1

class HappensAfterGraph;

class DCTraceBuilder : public TSOTraceBuilder {
#ifdef DUMP_ANNOT_TREE
  // for debugging
  AnnotTree annot_tree;
#endif

  // number of leaves of recursion tree
  unsigned executed_traces = 0;
  unsigned leaves_number = 0;
  unsigned succ_leaves_number = 0;
  // executed instructions in final leaves
  unsigned instr_executed = 0;
  // number of unrealizable traces that were
  // identified to be unrealizable before
  // reducing to SAT
  unsigned blocked_before_sat = 0;
  // number of tried mutations
  unsigned blocked_after_sat = 0;
  unsigned blocked_hb = 0;
  unsigned blocked_past = 0;
  unsigned blocked_past2 = 0;
  unsigned blocked_in_mutation = 0;
  unsigned realize_called = 0;
  unsigned realize_succeeded = 0;
  unsigned realized_using_swap = 0;
  unsigned merged_traces = 0;

  const Configuration &config;
  llvm::Module *M;
  int next_to_schedule;

  // this event works as the global's initialization event
  // We could always put it into the prefix, but then
  // we would need to change the whole code that works
  // with prefix, so it is better here...
  // We need this event when creating well-formed
  // annotations
  const DCEvent initial_event;

  // set to true when we are looking for some
  // maximal-length trace from the point where
  // we set this. Special case is at the beginning,
  // where we set it to true to get some initial
  // (arbitrary) maximal trace
  bool getting_maximal_extension;
  bool getting_initial_trace;
  // should we merge traces
  bool should_merge_traces = false;

  Trace *error_trace = nullptr;

  // this is the currently executed instruction
  const llvm::Instruction *current_inst = nullptr;
  void mayConflict(const ConstMRef *ml = nullptr);

  // explore the given trace -- try mutations and recursively
  // call itself on realized traces. Return an error trace
  // if found any or nullptr.
  void explore(AnnotatedTrace& trace);

  void mergeTraces(AnnotatedTrace& trace,
                   AnnotatedTrace& trace2,
                   const PositiveAnnotation& annot);

  std::vector<DCEvent> extendCurrentTrace();
  // extend current annotated trace (the extension is stored
  // again to trace.trace). If the extended trace finds an error,
  // return the error trace, otherwise return nullptr.
  std::vector<DCEvent> extendTrace(std::vector<DCEvent>&& tr);

  // try adding a new pair to positive annotation,
  // return false if this annotation can not be realized for
  // some reason
  bool addMutation(PositiveAnnotation&, const DCIID&, const DCIID&);
  bool addMergedAnnotation(AnnotatedTrace& trace, HappensAfterGraph& currentPO,
                           PositiveAnnotation& positive_annotation,
                           const DCEvent& ev, const DCEvent& wr_ev);

  bool isInFutureCone(std::vector<DCEvent>& trace, int read_idx, int write_idx) const;
  bool hasAnnotationFreeFutureCone(std::vector<DCEvent>& trace, int idx, int idx2, const PositiveAnnotation& annot) const;
  bool addPastCone(const PositiveAnnotation& cone,
                   PositiveAnnotation& annot,
                   const NegativeAnnotation& negative_annotation);

#ifdef DUMP_ANNOT_TREE
  void addAnnotTreeNode(AnnotatedTrace& trace, const char *msg = nullptr) {
    trace.annot_tree_node
        = annot_tree.addNode(getCurrentTrace().annot_tree_node,
                             trace.positive_annotation,
                             trace.negative_annotation,
                             trace.happens_before,
                             msg);
  }
#endif

  std::vector<DCEvent> swapWithoutSAT(std::vector<DCEvent>& trace, int read_idx, int write_idx);
  std::vector<DCEvent> moveWithoutSAT(std::vector<DCEvent>& trace, int write_idx, int read_idx);

  template <typename T>
  void updateNegativeAnnotation(NegativeAnnotation& neg_annot,
                                const T& ev, const T& wr_ev,
                                PositiveAnnotation& past_annot) {
    neg_annot.add(ev, wr_ev, past_annot);
  }
  template <typename T>
  void updateNegativeAnnotation(NegativeAnnotation& neg_annot,
                                const T& ev, const T& wr_ev,
                                PositiveAnnotation&& past_annot) {
    neg_annot.add(ev, wr_ev, past_annot);
  }

  // check whether we can make the read to see the write
  bool checkMutation(const HappensAfterGraph& currentPO,
                     std::vector<DCEvent>& trace,
                     int read_idx, int write_idx);

  // the trace that is currently being scheduled
  AnnotatedTrace current_trace;

  void update_prefix(unsigned p);
  bool schedule_thread(int *proc, unsigned p, bool dont_update_prefix = false);
  bool schedule_arbitrary_thread(int *proc);
  bool schedule_current_trace(int *proc);

  // stack of traces to be explored
  std::vector<AnnotatedTrace> traces;

  // will over-take the owner ship of resources
  void addTrace(AnnotatedTrace& tr) {
    traces.push_back(std::move(tr));
  }

  AnnotatedTrace& getCurrentTrace() { return current_trace; }
  const AnnotatedTrace& getCurrentTrace() const { return current_trace; }

  // add necessary happens-before edges that are derived from annotations
  void initializeHappensBefore(AnnotatedTrace& trace, Basis& basis);

  bool tryRealizeAnnotations(AnnotatedTrace& trace,
                             HappensExactlyBefore& happens_before,
                             const Basis& basis);

  // try to realize the annotated prefix
  bool realize(AnnotatedTrace& trace, Basis& basis);
  // try to realize the annotation using direct swap of events
  void realizeUsingSwap(AnnotatedTrace& trace,
                        PositiveAnnotation& positive_annotation,
                        PositiveAnnotation& past_annot,
                        int read_idx, int write_idx);
  bool trySwap(AnnotatedTrace& trace,
               PositiveAnnotation& positive_annotation,
               PositiveAnnotation& past_annot,
               int read_idx, int write_idx);
 
  void tryRealizeMutation(HappensAfterGraph& currentPO,
                          AnnotatedTrace& trace,
                          int read_idx,
                          int write_idx);

  bool tryRealizeMutationToInit(int read_idx,
                                AnnotatedTrace& trace,
                                HappensAfterGraph& currentPO,
                                PositiveAnnotation& observation);

  const std::vector<DCEvent>& getPrefix() const { return prefix; }

  void dump_prefix() const;
  void dump_current_trace() const;

public:
  DCTraceBuilder(const Configuration &conf, llvm::Module *m);
  DCTraceBuilder(const Configuration &conf,
                 llvm::Module *m, std::vector<DCEvent>&& tr);

  virtual bool reset();
  virtual bool schedule(int *proc, int *aux, int *alt, bool *dryrun);
  virtual void schedule_next(unsigned idx) { next_to_schedule = 2*idx; }

  virtual void executing_instruction(const llvm::Instruction *Instr) {
    current_inst = Instr;
    ++instr_executed;

    unsigned p = curnode().iid.get_pid();
    curnode().order = threads[p].clock[p];

  }

  /* for debugging */
  //void setValue(int v) { curnode().val = v; }
  virtual bool has_error() const { return (error_trace != nullptr)
                                            || (errors.size() > 0); };

  virtual void atomic_store(const ConstMRef &ml);
  virtual void load(const ConstMRef &ml);
  virtual void refuse_schedule();
  virtual void spawn();
  virtual void join(int tgt_proc);
  virtual void fence();

  virtual void mutex_lock(const ConstMRef &ml);
  //virtual void mutex_lock_fail(const ConstMRef &ml);
  virtual void mutex_trylock(const ConstMRef &ml);
  virtual void mutex_unlock(const ConstMRef &ml);
  virtual void mutex_init(const ConstMRef &ml);
  virtual void mutex_destroy(const ConstMRef &ml);
  virtual void metadata(const llvm::MDNode *md);
  virtual IID<CPid> get_iid() const;
  virtual Trace *get_trace() const;

  virtual void full_memory_conflict() {
    not_supported("full memory conflict");
  }

  virtual void store(const ConstMRef &ml) {
    not_supported("relaxed memory models");
  }

  virtual bool cond_init(const ConstMRef &ml) {
	  // do nothing. If the program just creates
	  // the cond variable but does not use it,
	  // we can work
	return true;
  }

  virtual bool cond_signal(const ConstMRef &ml) {
    not_supported("signaling");
	return false;
  }

  virtual bool cond_broadcast(const ConstMRef &ml) {
    not_supported("signaling");
	return false;
  }

  virtual bool cond_wait(const ConstMRef &cond_ml, const ConstMRef &mutex_ml) {
    not_supported("signaling");
	return false;
  }

  virtual int cond_destroy(const ConstMRef &ml) {
	  // the same as cond_init
	  return 0;
  }

protected:
  /* The fixed prefix of events in the current execution. This may be
   * either the complete sequence of events executed thus far in the
   * execution, or the events executed followed by the subsequent
   * events that are determined in advance to be executed.
   */
  std::vector<DCEvent> prefix;

  DCEvent &curnode() {
    assert(0 <= prefix_idx);
    assert(prefix_idx < int(prefix.size()));
    return prefix[prefix_idx];
  };

  const DCEvent &curnode() const {
    assert(0 <= prefix_idx);
    assert(prefix_idx < int(prefix.size()));
    return prefix[prefix_idx];
  };
};

#endif

