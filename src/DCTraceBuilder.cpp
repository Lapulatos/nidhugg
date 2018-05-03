/* Copyright (C) 2014-2016 Carl Leonardsson
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

#include "Debug.h"
#include "DCTraceBuilder.h"
#include "Annotations.h"
#include "Basis.h"
#include "HappensAfterGraph.h"
#include "HappensAfterGraphVclock.h"
#include "HappensExactlyBeforeGen.h"
#include "DCInterpreter.h"
#include "DCTBHelpers.h"
#include "DPORDriver.h"

#include <llvm/IR/Instructions.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const DCIID& iid);

static inline const DCIID evToDCIID(const DCEvent& ev) {
    return DCIID(ev.cpid, ev.instruction, ev.order);
}

DCTraceBuilder::DCTraceBuilder(const Configuration &conf, llvm::Module *m)
: TSOTraceBuilder(conf), config(conf), M(m),
  next_to_schedule(-1), initial_event(IID<IPid>(0,0)),
  getting_maximal_extension(false), getting_initial_trace(true)
{
}

// start with this trace
DCTraceBuilder::DCTraceBuilder(const Configuration &conf,
                               llvm::Module *m, std::vector<DCEvent>&& trace)
: TSOTraceBuilder(conf), config(conf), M(m),
  next_to_schedule(-1), initial_event(IID<IPid>(0,0)),
  getting_maximal_extension(false), getting_initial_trace(false)
{
  current_trace.trace = std::move(trace);
}

bool DCTraceBuilder::schedule_current_trace(int *proc)
{
  AnnotatedTrace& current_annotated_trace = getCurrentTrace();
  auto& trace = current_annotated_trace.trace;

  assert(!trace.empty());
  assert(prefix_idx < (int) trace.size());

  // the execution copies the the events in trace,
  // so the prefix_idx is also index into the events
  if (prefix_idx == -1) {
    prefix_idx = 0;
  } else {
    // if we already executed the last instruction
    // in the event, continue with next event
    if (trace[prefix_idx].size == prefix[prefix_idx].size) {
      ++prefix_idx;
      // are we done replaying the trace?
      // So get a maximal extension.
      if (prefix_idx == (int) trace.size()) {
        // if we do not have any other event, try to do
        // the maximial extension of this execution:
        // set 'getting_maximal_extension', so that the
        // schedule_tread method will try to schedule
        // any available thread
        getting_maximal_extension = true;
        // decrease back the prefix_idx, so that the
        // shedule_arbitrary_thread can access the last
        // event in prefix - it will bump the prefix back
        // if needed
        --prefix_idx;
        return schedule_arbitrary_thread(proc);
      } else {
        assert(trace.size() > prefix.size());
      }
    } else {
      assert(trace[prefix_idx].size > prefix[prefix_idx].size);
      // bump the size of current event, it will be scheduled
      ++prefix[prefix_idx].size;
    }
  }

  assert((unsigned)prefix_idx < trace.size());
  unsigned p = trace[prefix_idx].iid.get_pid();
  // increase clock
  ++threads[p].clock[p];


  // did we bumped the prefix?
  if (prefix_idx == (int) prefix.size()) {
      prefix.emplace_back(IID<IPid>(IPid(p),threads[p].clock[p]),
                          threads[p].cpid, prefix.size());
      assert(prefix.back().id == prefix.size() - 1);
  }

  bool ret = schedule_thread(proc, p, true /* do not update prefix */);
  assert(ret && "Bug in scheduling");

  return ret;
}

void DCTraceBuilder::update_prefix(unsigned p)
{
    // increase clock
    ++threads[p].clock[p];

    if ((prefix_idx != -1 && !curnode().may_conflict
         && (int) p == curnode().iid.get_pid())) {
      assert(prefix_idx == (int) (prefix.size() - 1));
      ++prefix[prefix_idx].size;
    } else {
      // extend the prefix
      prefix.emplace_back(IID<IPid>(IPid(p),threads[p].clock[p]),
                               threads[p].cpid, prefix.size());
      ++prefix_idx;

      assert(prefix.back().id == prefix.size() - 1);
      assert((unsigned) prefix_idx == prefix.size() - 1);
    }
}

// schedule a next event to be the next event from thread p
bool DCTraceBuilder::schedule_thread(int *proc, unsigned p,
                                     bool dont_update_prefix)
{
  if(threads[p].available && !threads[p].sleeping &&
     (conf.max_search_depth < 0
      || threads[p].clock[p] < conf.max_search_depth)) {
    if (!dont_update_prefix) {
      update_prefix(p);
    } else {
      assert(!getting_maximal_extension && !getting_initial_trace);
    }

    // set the scheduled thread
    *proc = p/2;

    return true;
  }

  return false;
}

bool DCTraceBuilder::schedule_arbitrary_thread(int *proc)
{
  assert(getting_maximal_extension || getting_initial_trace);

  const unsigned sz = threads.size();
  unsigned p;

  for(p = 0; p < sz; p += 2) {
    if (schedule_thread(proc, p))
      return true;
  }

  // we did not schedule anything?
  return false;
}


bool DCTraceBuilder::schedule(int *proc, int *aux, int *alt, bool *dryrun)
{
  *dryrun = false;
  *alt = 0;
  *aux = -1;
  this->dryrun = false;

  // if we are getting some initial trace, just schedule
  // some available thread
  if (getting_maximal_extension || getting_initial_trace)
    return schedule_arbitrary_thread(proc);

  // replay the trace that we have computed
  // and stored into current_trace attribute
  assert(!getting_maximal_extension);
  assert(!getting_initial_trace);
  return schedule_current_trace(proc);
}

bool DCTraceBuilder::addMutation(PositiveAnnotation& annotation,
                                 const DCIID& ev, const DCIID& wr_ev)
{
  assert(ev.instruction);

  // the mapping must be unique (injective) on locks, that is,
  // we do not want an unlock be observed by more locks
  // XXX: we could do it somehow efficiently. Maybe inverse mapping?
  if (is_function_call(ev.instruction, "pthread_mutex_lock")) {
    for (auto& rw_pair : annotation) {
      assert(rw_pair.first.instruction && "read has no instruction");
      if (// is read event lock?
          is_function_call(rw_pair.first.instruction, "pthread_mutex_lock") &&
          // is write event same as the one we want to add?
          rw_pair.second.instruction == wr_ev.instruction &&
          rw_pair.second.cpid == wr_ev.cpid) {
        // fail!
        return false;
      }
    }
  }

  // add a new annotation pair. If we try to add a pair
  // that is already in annotation, the annotation is
  // unrealizable, so we can bail out (the annotation
  // would not be a function E(t) -> E(t))
  return annotation.add(ev, wr_ev);
}

static bool checkAnnotation(std::vector<DCEvent>& trace, const PositiveAnnotation& annot)
{
  // go from right to left and once you find annotated read,
  // all reads from the same process to the left must be
  // annotated too
  std::set<CPid> must_be_annotated;
  for (int idx = trace.size() - 1; idx >= 0; --idx) {
    const DCEvent& ev = trace[idx];
    if (isRead(ev)) {
      bool must = must_be_annotated.count(ev.cpid) > 0;
      bool annotated = annot.defines(ev);

      if (must && !annotated)
        return false;

      if (annotated && !must)
        must_be_annotated.insert(ev.cpid);
    }
  }

  return true;
}

bool DCTraceBuilder::addPastCone(const PositiveAnnotation& past_annot,
                                 PositiveAnnotation& annot,
                                 const NegativeAnnotation& negative_annotation)
{
  for (auto& annot_pair : past_annot) {
    // check whether we already have an annotation for the write
    // annot_pair.first. If so, it must be the same (compatible),
    // otherwise it is unrealizible
    if (auto x = annot.get(annot_pair.first)) {
      if (*x == annot_pair.second) {
        continue;
      } else {
        return false;
      }
    }

    if (!addMutation(annot, annot_pair.first, annot_pair.second)) {
      return false;
    }
  }

  return true;
}

static void rollbackGraph(HappensAfterGraph& G,
                          HappensAfterGraph::Node *read,
                          HappensAfterGraph::Node *old_obs,
                          HappensAfterGraph::Node *new_obs)
{
  new_obs->removeSuccessor(read);
  if (old_obs)
    old_obs->addSuccessor(read);
}

static std::tuple<HappensAfterGraph::Node *,
                  HappensAfterGraph::Node *,
                  HappensAfterGraph::Node *> rearrangeGraph(HappensAfterGraph& currentPO,
                                                            const std::vector<DCEvent>& trace,
                                                            int read_idx, int write_idx) {
  using Node = HappensAfterGraph::Node;

  // get the observation function of current trace
  const PositiveAnnotation& O = currentPO.getAnnotation();
  const DCEvent& ev = trace[read_idx];

  auto evnd = currentPO.getNode(&ev);
  assert(evnd);
  Node *current_observation_node = nullptr;

  // remove the successors now and replace it by the new successor
  // -- but don't remove the successor if the read and the currently observed
  // write are in the same event, because then the successor is also implied
  // by the structure
  auto current_observation = O.get(ev);
  assert(current_observation);
  if (current_observation->cpid != ev.cpid) {
    current_observation_node = currentPO.getNode(*current_observation);
    current_observation_node->removeSuccessor(evnd);
  }

  Node *new_observation
    = write_idx == -1 ? currentPO.getNode(nullptr) : currentPO.getNode(&trace[write_idx]);
  assert(new_observation);
  new_observation->addSuccessor(evnd);

  return std::tuple<Node *, Node *, Node *>(evnd, current_observation_node, new_observation);
}

// create a new negative annotation by copying 'na', but only those parts that
// are not annotated by 'annot'
static NegativeAnnotation pruneNegativeAnnotation(const PositiveAnnotation& annot,
                                                  const NegativeAnnotation& na)
{
    NegativeAnnotation ret;
    for (const auto& it : na) {
      if (annot.defines(it.first))
        continue;

      for (const auto& it2 : it.second) {
        ret.add(it.first, it2.first, it2.second);
      }
    }

    return ret;
}

void DCTraceBuilder::realizeUsingSwap(AnnotatedTrace& trace,
                                      PositiveAnnotation& positive_annotation,
                                      PositiveAnnotation& past_annot,
                                      int read_idx, int write_idx) {
  // if the write that we should observe is after the read,
  // we can do other magic than with the formula
  assert(read_idx < write_idx);
  // we check this in checkMutation()
  assert(!isInFutureCone(trace.trace, read_idx, write_idx));

  const DCEvent& ev = trace.trace[read_idx];
  const DCEvent& wr_ev = trace.trace[write_idx];

  AnnotatedTrace new_trace;
  new_trace.trace = swapWithoutSAT(trace.trace, read_idx, write_idx);

  auto tr = extendTrace(std::move(new_trace.trace));
  if (has_error())
    return;

  new_trace.negative_annotation
    = pruneNegativeAnnotation(positive_annotation,
                              trace.negative_annotation);
  new_trace.positive_annotation = std::move(positive_annotation);
  new_trace.happens_before = trace.happens_before;

#ifdef DUMP_ANNOT_TREE
  addAnnotTreeNode(new_trace, "(swapped)");
#endif

  new_trace.trace = std::move(tr);
  // recursively explore new traces
  // we'll explore new traces after merging some of them,
  // so just store the new trace now
  //explore(new_trace);
  traces.emplace_back(std::move(new_trace));

  updateNegativeAnnotation(trace.negative_annotation,
                           ev, wr_ev, past_annot);
  for (auto& rw_pair : past_annot)
    updateNegativeAnnotation(trace.negative_annotation,
                             rw_pair.first, rw_pair.second, {});
  ++realized_using_swap;
}

bool DCTraceBuilder::trySwap(AnnotatedTrace& trace,
                             PositiveAnnotation& positive_annotation,
                             PositiveAnnotation& past_annot,
                             int read_idx, int write_idx) {
  assert(write_idx < read_idx);

  //  a. If w has no annotated read in its observation future:
  //    a.i.  If w and r are in the same process, we can cut the segment from w
  //          to r in that process and place it after the current observation of r
  //    a.ii. If w and r are in different processes, you can remove the future
  //          cone of w and make it happen exactly before r


  if (!hasAnnotationFreeFutureCone(trace.trace, write_idx, read_idx, positive_annotation)) {
    return false;
  }

  const DCEvent& ev = trace.trace[read_idx];
  const DCEvent& wr_ev = trace.trace[write_idx];

  std::vector<DCEvent> new_events_sequence;
  if (ev.iid.get_pid() != wr_ev.iid.get_pid()) {
    new_events_sequence = moveWithoutSAT(trace.trace, write_idx, read_idx);
  } else {
    llvm::errs() << "SWAPPING read after write\n";
    new_events_sequence.reserve(trace.trace.size());

    // copy the part of new_events_sequence until the write
    for (int i = 0; i < write_idx; ++i) {
      new_events_sequence.push_back(trace.trace[i]);
    }
    // skip write and copy until processes from other
    // events until the read
    std::vector<DCEvent> to_put_before_read;
    for (int i = write_idx + 1; i < read_idx; ++i) {
      if (trace.trace[i].iid.get_pid() == ev.iid.get_pid())
        to_put_before_read.push_back(trace.trace[i]);
      else
        new_events_sequence.push_back(trace.trace[i]);
    }

    // now put here everything from the read's process
    for (auto& ne : to_put_before_read)
      new_events_sequence.push_back(ne);

    // now put there the read and the rest of the trace
    for (int i = read_idx; i < trace.trace.size(); ++i)
      new_events_sequence.push_back(trace.trace[i]);
  }

  AnnotatedTrace new_trace;
  new_trace.trace = std::move(new_events_sequence);

  auto tr = extendTrace(std::move(new_trace.trace));
  if (has_error())
    return true;

  new_trace.negative_annotation
    = pruneNegativeAnnotation(positive_annotation,
                              trace.negative_annotation);
  new_trace.positive_annotation = std::move(positive_annotation);
  new_trace.happens_before = trace.happens_before;

#ifdef DUMP_ANNOT_TREE
  addAnnotTreeNode(new_trace, "(swapped)");
#endif

  new_trace.trace = std::move(tr);
  // recursively explore new traces
  // we'll explore new traces after merging some of them,
  // so just store the new trace now
  //explore(new_trace);
  traces.emplace_back(std::move(new_trace));

  updateNegativeAnnotation(trace.negative_annotation,
                           ev, wr_ev, past_annot);
  for (auto& rw_pair : past_annot)
    updateNegativeAnnotation(trace.negative_annotation,
                             rw_pair.first, rw_pair.second, {});
  ++realized_using_swap;
  return true;
}

void DCTraceBuilder::tryRealizeMutation(HappensAfterGraph& currentPO,
                                        AnnotatedTrace& trace,
                                        int read_idx, int write_idx)
{
  // we could bail out earlier if this would hold
  assert(checkMutation(currentPO, trace.trace, read_idx, write_idx));

  // take the current graph and rearrange the edges according to the
  // read-write pair we are trying to mutate. This graph we can
  // use to derive the past cones of events and so on.
  HappensAfterGraph::Node *evnd, *current_observation_node, *new_observation;
  std::tie(evnd, current_observation_node, new_observation)
    = rearrangeGraph(currentPO, trace.trace, read_idx, write_idx);

  const DCEvent& ev = trace.trace[read_idx];
  const DCEvent& wr_ev = write_idx == -1 ? initial_event : trace.trace[write_idx];
  // create a copy of the positive annotation, we will modify it
  PositiveAnnotation positive_annotation = trace.positive_annotation;

  // add the annotation that we want
  if (!addMutation(positive_annotation,
                   evToDCIID(ev), evToDCIID(wr_ev))) {
    // return back original successors so that we can reuse the graph
    rollbackGraph(currentPO, evnd, current_observation_node, new_observation);
    return;
  }

  /*
  if (!addMergedAnnotation(trace, currentPO, positive_annotation, ev, wr_ev)) {
      rollbackGraph(currentPO, evnd, current_observation_node, new_observation);
      return;
  }
  */

  // get a basis that contains only the annotated events
  Basis basis = currentPO.getRestrictedBasis(positive_annotation);
  assert(!basis.empty());

  // get the past cone of the write and read
  auto past_annot
    = currentPO.getPastConeAnnotation(positive_annotation, {&ev, &wr_ev});
  assert(!past_annot.defines(wr_ev));

  if (!addPastCone(past_annot, positive_annotation,
                   trace.negative_annotation)) {
    rollbackGraph(currentPO, evnd, current_observation_node, new_observation);
    return;
  }

  // if the write that we should observe is after the read,
  // we can swap the events directly without using 2SAT reduction
  // (the cases where we can not do it were blocked in checkMutation())
  if (read_idx < write_idx) {
    realizeUsingSwap(trace, positive_annotation, past_annot, read_idx, write_idx);
    assert(positive_annotation.empty() || has_error()); // we moved the positive annotation
    rollbackGraph(currentPO, evnd, current_observation_node, new_observation);
    return;
  }

  assert(write_idx < read_idx);

  // if the write is before the read, we can still try to swap it
  if (trySwap(trace, positive_annotation, past_annot, read_idx, write_idx)) {
    assert(positive_annotation.empty() || has_error()); // we moved the positive annotation
    rollbackGraph(currentPO, evnd, current_observation_node, new_observation);
    return;
  }

  assert(checkAnnotation(trace.trace, positive_annotation)
         && "Did not annotated all that we need");

  AnnotatedTrace new_trace;
  new_trace.positive_annotation = std::move(positive_annotation);
  new_trace.negative_annotation = trace.negative_annotation;
  new_trace.happens_before = trace.happens_before;

  // realize the new trace (if it is possible)
  if (realize(new_trace, basis)) {
    updateNegativeAnnotation(trace.negative_annotation,
                             ev, wr_ev, past_annot);
    for (auto& rw_pair : past_annot)
      updateNegativeAnnotation(trace.negative_annotation,
                               rw_pair.first, rw_pair.second, {});
    // we queue new traces for execution in the realize method
  }

  rollbackGraph(currentPO, evnd, current_observation_node, new_observation);
}


// is the idx2 in the future cone of idx?
bool DCTraceBuilder::isInFutureCone(std::vector<DCEvent>& trace, int idx, int idx2) const
{
    std::set<ConstMRef> dep_mem;
    // pids to ommit because they are in the future cone
    std::set<unsigned> in_future_cone;
    in_future_cone.insert(trace[idx].iid.get_pid());
    for (int i = idx, e = trace.size(); i < e; ++i) {
      if (i == idx2) {
        return in_future_cone.count(trace[i].iid.get_pid()) > 0;
      } else if (in_future_cone.count(trace[i].iid.get_pid()) > 0) {
        // if this is a write in the future cone,
        // then now add every read that sees that write
        if (isWrite(trace[i]))
          dep_mem.insert(trace[i].ml);
        else if (isPthreadCreate(trace[i])) {
          // if this is a thread creation from an event
          // that is in the future cone, then the new thread
          // is also in the future cone
          for (int ii = 0, ee = threads.size(); ii < ee; ++ii) {
            if (trace[i].childs_cpid == threads[ii].cpid) {
              in_future_cone.insert(ii);
              break;
            }
          }
        }
      } else {
        if (isRead(trace[i])) {
          for (const ConstMRef& mem : dep_mem) {
            if (mem.overlaps(trace[i].ml)) {
              in_future_cone.insert(trace[i].iid.get_pid());
              break;
            }
          }
        } else if (isWrite(trace[i])) {
          // if we already tracked this memory, that means
          // that this write overwrites some previous write.
          // So anything that goes after this write is already
          // not in our future cone
          // NOTE: here we assume that we overwrite the whole variable always
          auto it = dep_mem.find(trace[i].ml);
          if (it != dep_mem.end()) {
              dep_mem.erase(it);
          }
        } else if (isJoin(trace[i])) {
          // If we join a thread that is in the future cone,
          // than this join is also in the future cone.
          for (int x : in_future_cone) {
            if (threads[x].cpid == trace[i].childs_cpid) {
              in_future_cone.insert(trace[i].iid.get_pid());
              break;
            }
          }
        }
      }
    }

    assert(0 && "Unreachable");
    abort();
}

// is the idx2 in the future cone of idx?
bool DCTraceBuilder::hasAnnotationFreeFutureCone(std::vector<DCEvent>& trace, int idx, int read_idx, const PositiveAnnotation& annot) const
{
    // initial write?
    if (idx == -1)
      return false;

    std::set<ConstMRef> dep_mem;
    // pids to ommit because they are in the future cone
    std::set<unsigned> in_future_cone;
    in_future_cone.insert(trace[idx].iid.get_pid());
    for (int i = idx, e = trace.size(); i < e; ++i) {
      if (i == read_idx) {
        //if the read is in the future cone in a different process,
        // than we fail
        if (trace[read_idx].iid.get_pid() != trace[idx].iid.get_pid()
            && in_future_cone.count(trace[i].iid.get_pid()) > 0)
          return false;
        // all reads between idx and read_idx in the future cone
        // are not annotated
        return true;
      } else if (in_future_cone.count(trace[i].iid.get_pid()) > 0) {
        // if this is a write in the future cone,
        // then now add every read that sees that write
        if (isWrite(trace[i])) {
          dep_mem.insert(trace[i].ml);
		    } else if (isRead(trace[i])) {
          if (annot.get(evToDCIID(trace[i])))
            return false;
        } else if (isPthreadCreate(trace[i])) {
          // if this is a thread creation from an event
          // that is in the future cone, then the new thread
          // is also in the future cone
          for (int ii = 0, ee = threads.size(); ii < ee; ++ii) {
            if (trace[i].childs_cpid == threads[ii].cpid) {
              in_future_cone.insert(ii);
              break;
            }
          }
        }
      } else {
        if (isRead(trace[i])) {
          if (annot.get(evToDCIID(trace[i])))
            return false;

          for (const ConstMRef& mem : dep_mem) {
            if (mem.overlaps(trace[i].ml)) {
              in_future_cone.insert(trace[i].iid.get_pid());
              break;
            }
          }
        } else if (isWrite(trace[i])) {
          // if we already tracked this memory, that means
          // that this write overwrites some previous write.
          // So anything that goes after this write is already
          // not in our future cone
          // NOTE: here we assume that we overwrite the whole variable always
          auto it = dep_mem.find(trace[i].ml);
          if (it != dep_mem.end()) {
              dep_mem.erase(it);
          }
        } else if (isJoin(trace[i])) {
          // If we join a thread that is in the future cone,
          // than this join is also in the future cone.
          for (int x : in_future_cone) {
            if (threads[x].cpid == trace[i].childs_cpid) {
              in_future_cone.insert(trace[i].iid.get_pid());
              break;
            }
          }
        }
      }
    }

    return true;
}


std::vector<DCEvent> DCTraceBuilder::swapWithoutSAT(std::vector<DCEvent>& trace, int read_idx, int write_idx)
{
    std::vector<DCEvent> new_trace;
    new_trace.reserve(trace.size());

    std::set<ConstMRef> dep_mem;

    // copy the part of new_trace until the read
    for (int i = 0; i < read_idx; ++i) {
      new_trace.push_back(trace[i]);
    }

    // pids to ommit because they are in the future cone
    std::set<unsigned> to_omit;
    to_omit.insert(trace[read_idx].iid.get_pid());
    for (int i = read_idx, e = trace.size(); i < e; ++i) {
      if (i == write_idx) {
        new_trace.push_back(trace[i]);
        new_trace.back().id = new_trace.size() - 1;
        new_trace.push_back(trace[read_idx]);
        new_trace.back().id = new_trace.size() - 1;
        // we are not done, still need to copy the reset
        // of new_trace (there may be some annotated events
        // that we need in the new_trace too)
      } else if (to_omit.count(trace[i].iid.get_pid()) > 0) {
        // if this is a write in the future cone,
        // then now add every read that sees that write
        if (isWrite(trace[i]))
          dep_mem.insert(trace[i].ml);

        continue;
      } else {
        if (isRead(trace[i])) {
          bool omit = false;
          for (const ConstMRef& mem : dep_mem) {
            if (mem.overlaps(trace[i].ml)) {
              to_omit.insert(trace[i].iid.get_pid());
              omit = true;
              break;
            }
          }
          if (omit)
            continue;
        } else if (isWrite(trace[i])) {
          // if we already tracked this memory, that means
          // that this write overwrites some previous write.
          // So anything that goes after this write is already
          // not in our future cone
          // NOTE: here we assume that we overwrite the whole variable always
          auto it = dep_mem.find(trace[i].ml);
          if (it != dep_mem.end()) {
              dep_mem.erase(it);
          }
        } else if (trace[i].instruction && isJoin(trace[i].instruction)) {
          // if we join a thread that is in the future cone,
          // that this join is also in the future cone
          bool omit = false;
          for (int x : to_omit) {
            if (threads[x].cpid == trace[i].childs_cpid) {
              omit = true;
              break;
            }
          }
          if (omit) {
            // joins are in the future cone
            to_omit.insert(trace[i].iid.get_pid());
            continue;
          }
        }

        new_trace.push_back(trace[i]);
        new_trace.back().id = new_trace.size() - 1;
      }
    }

    return new_trace;
}

std::vector<DCEvent> DCTraceBuilder::moveWithoutSAT(std::vector<DCEvent>& trace, int write_idx, int read_idx)
{
  std::vector<DCEvent> new_trace;
  new_trace.reserve(trace.size());

  std::set<ConstMRef> dep_mem;

  // copy the part of new_trace until the write
  for (int i = 0; i < write_idx; ++i) {
    new_trace.push_back(trace[i]);
  }

  // pids to ommit because they are in the future cone
  std::set<unsigned> to_omit;
  to_omit.insert(trace[write_idx].iid.get_pid());
  for (int i = write_idx, e = trace.size(); i < e; ++i) {
    if (i == read_idx) {
      new_trace.push_back(trace[write_idx]);
      new_trace.back().id = new_trace.size() - 1;
      new_trace.push_back(trace[i]);
      new_trace.back().id = new_trace.size() - 1;
      // we are done, they may not be any annotated events that
      // we should keep in the trace
      break;
    } else if (to_omit.count(trace[i].iid.get_pid()) > 0) {
      // if this is a write in the future cone,
      // then now add every read that sees that write
      if (isWrite(trace[i]))
        dep_mem.insert(trace[i].ml);

      continue;
    } else {
      if (isRead(trace[i])) {
        bool omit = false;
        for (const ConstMRef& mem : dep_mem) {
          if (mem.overlaps(trace[i].ml)) {
            to_omit.insert(trace[i].iid.get_pid());
            omit = true;
            break;
          }
        }
        if (omit)
          continue;
      } else if (isWrite(trace[i])) {
        // if we already tracked this memory, that means
        // that this write overwrites some previous write.
        // So anything that goes after this write is already
        // not in our future cone
        // NOTE: here we assume that we overwrite the whole variable always
        auto it = dep_mem.find(trace[i].ml);
        if (it != dep_mem.end()) {
            dep_mem.erase(it);
        }
      } else if (trace[i].instruction && isJoin(trace[i].instruction)) {
        // if we join a thread that is in the future cone,
        // that this join is also in the future cone
        bool omit = false;
        for (int x : to_omit) {
          if (threads[x].cpid == trace[i].childs_cpid) {
            omit = true;
            break;
          }
        }
        if (omit) {
          // joins are in the future cone
          to_omit.insert(trace[i].iid.get_pid());
          continue;
        }
      }

      new_trace.push_back(trace[i]);
      new_trace.back().id = new_trace.size() - 1;
    }
  }

  return new_trace;
}



bool DCTraceBuilder::checkMutation(const HappensAfterGraph& currentPO,
                                   std::vector<DCEvent>& trace,
                                   int read_idx, int write_idx)
{
  assert(write_idx == -1 || trace[read_idx].ml.overlaps(trace[write_idx].ml));

  // if the write (or any other instruction) is after the read
  // in the same process, we have no chance to make
  // the read observe the write
  // FIXME: isn't this covered in the later cases that
  // check the future cone?
  if (write_idx != -1 &&
      (trace[write_idx].iid.get_pid() == trace[read_idx].iid.get_pid())) {
    if (write_idx > read_idx) {
        ++blocked_past;
        return false;
    }

    // if this is a write in the same thread, check whether
    // it is the closest one, because it is not possible
    // for the read to observe other writes in the same thread
    // than the "closest" one
    assert(write_idx < read_idx);
    bool found_write_in_between = false;
    for (int write_idx2 = write_idx + 1; write_idx2 < read_idx; ++write_idx2) {
        DCEvent& wr_ev2 = trace[write_idx2];
        if (isWrite(wr_ev2) && wr_ev2.ml.overlaps(trace[read_idx].ml)
            && wr_ev2.iid.get_pid() == trace[read_idx].iid.get_pid()) {
            found_write_in_between = true;
            break;
        }
    }

    // do not try annotating this write
    if (found_write_in_between) {
      ++blocked_past2;
      return false;
    }
  }

  // if the write that we should observe is after the read,
  // we can do other magic than with the formula
  if (read_idx < write_idx) {
    if (isInFutureCone(trace, read_idx, write_idx)) {
      // if the write is after the read, but read is in the
      // past cone, then making the read observe the write is
      // impossible
      ++blocked_past;
      return false;
    }
  } else {
    // if there is a conflicting write before that is in the causal past
    // of r and also the write_idx is in the causal past of this write,
    // then the trace is unrealizable
    assert(write_idx < read_idx);
    for (int i = write_idx == -1 ? 0 : write_idx; i < read_idx; ++i) {
      if (!trace[i].may_conflict || !trace[i].instruction || !isWrite(trace[i]))
        continue;
      if (!trace[i].ml.overlaps(trace[read_idx].ml))
        continue;

      // we found conflicting write in between
      int wr_id = trace[i].id;

      // the observed write is in between of 'wr_ev' and 'ev'?
      if (wr_id > write_idx) {
        assert(wr_id < read_idx);
        if (isInFutureCone(trace, wr_id, read_idx)) { // we can not swap wr_id and read_idx
          // and we can not swap write_idx with wr_id
          if (write_idx == -1 || isInFutureCone(trace, write_idx, wr_id)) {
            ++blocked_past2;

            // nothing to do here, this would be unrealizable
            return false;
          }
        }
      }
    }
  }

  return true;
}

bool DCTraceBuilder::tryRealizeMutationToInit(int read_idx,
                                              AnnotatedTrace& trace,
                                              HappensAfterGraph& currentPO,
                                              PositiveAnnotation& observation)
{
  const DCEvent& ev = trace.trace[read_idx];

  if (trace.positive_annotation.defines(initial_event) ||
      trace.negative_annotation.forbids(trace.positive_annotation,
                                         ev, initial_event))
    return false;

  auto obs = observation.get(ev);
  if (obs && obs->instruction == nullptr) {
    // if the read already observes the initial event, skip it
    return false;
  } else if (checkMutation(currentPO, trace.trace, read_idx, -1)) {
    tryRealizeMutation(currentPO, trace,
                       read_idx, -1 /* initial event */);
    return true;
  }

  return false;
}


void DCTraceBuilder::explore(AnnotatedTrace& trace)
{
  // FIXME: make the following three things be computed at once, by one
  // traversal of the trace
  PositiveAnnotation observation = PositiveAnnotation::getObservationFunction(trace.trace);
  // events divided into threads
  Basis basis(trace.trace);
  // current partial order (happens-after graph) given by this trace
  HappensAfterGraph currentPO(observation, {}, basis);

  unsigned executed_traces_old = executed_traces;

  // iterate over reads in the trace from left to right
  unsigned trace_size = trace.trace.size();
  bool tried_realizing_something = false;
  for (unsigned read_idx = 0; read_idx < trace_size; ++read_idx) {
    DCEvent& ev = trace.trace[read_idx];

    // Find only global reads.
    // We are also interested only in the reads that we did not annotated yet
    if (!ev.may_conflict || !isRead(ev)
        || trace.positive_annotation.defines(ev))
      continue;

    // We have found non-annotated read. For every such read,
    // find conflicting writes.

    // (We also need to fake the initial write
    // -- it behaves like a conflicting write for every read)
    tried_realizing_something |=
      tryRealizeMutationToInit(read_idx, trace, currentPO, observation);

    // iterate over the trace and find conflicting writes
    for (unsigned write_idx = 0; write_idx < trace_size; ++write_idx) {
      DCEvent& wr_ev = trace.trace[write_idx];

      if (!wr_ev.may_conflict || !isWrite(wr_ev))
        continue;

      // is the read and write conflicting?
      if (!wr_ev.ml.overlaps(ev.ml))
        continue;

      // if the read already observes the write, skip it
      auto wr_id = observation.get(ev);
      if (wr_id && *wr_id == wr_ev)
        continue;

      // check whether this annotation may be realizable
      if (!checkMutation(currentPO, trace.trace, read_idx, write_idx))
        continue;

      // does the negative annotation forbid this event?
      if (trace.negative_annotation.forbids(trace.positive_annotation,
                                            ev, wr_ev)) {
        continue;
      }

      tried_realizing_something = true;
      tryRealizeMutation(currentPO, trace, read_idx, write_idx);
    }
  }

  // if we haven't tried realize something, then we are in the "recursion leaf",
  // since we have everything annotated - that means that we have explored
  // another class of the observation equivalence
  if (!tried_realizing_something)
    ++succ_leaves_number;

  // if we haven't recurse any further, this is a leaf of the recursion
  if (executed_traces == executed_traces_old)
    ++leaves_number;

  // explore the queued traces
  auto traces2 = std::move(traces);
  assert(traces.empty());
  for (auto& trace : traces2) {
      if (has_error())
          break;
      explore(trace);
  }
}

// compute a new trace from trace and given annotations
bool DCTraceBuilder::realize(AnnotatedTrace& trace, Basis& basis)
{
  assert(trace.trace.empty());
  assert(!trace.positive_annotation.empty());

  if (basis.size() < 3) {
    assert(!trace.has_root);
    // no point in trying happens-exactly before when the architecture
    // is acyclic for sure. This can help us also get possibly better
    // topology root. However, we must set one of the two processes
    // as a root, so that we have some root
    basis.setTopologyRoot(basis[0][0]->cpid);
    return tryRealizeAnnotations(trace, trace.happens_before, basis);
  } else {
    // we need to fix the root now
    if (!trace.has_root) {
      trace.topology_root = basis[1][0]->cpid;
      trace.has_root = true;
    }

    basis.setTopologyRoot(trace.topology_root);
    initializeHappensBefore(trace, basis);

    // if the annotation + happens before is cyclic already at this point,
    // do not try anything
    HappensAfterGraph G(trace.positive_annotation,
                        trace.happens_before, basis);
    if (!G.addNecessaryEdges()) {
#ifdef DUMP_ANNOT_TREE
      assert(trace.annot_tree_node == 0);
      addAnnotTreeNode(trace, "HB: Unrealizable before generating");
#endif
      ++blocked_hb;
      return false;
    }

    // addNecessaryEdges makes transitive closure
    // -- it is necessary, we rely on that while generating
    // happens exactly before.
    assert(G.isTransitivelyClosed());

    // generate the happens-before relation
    HappensExactlyBeforeGen HBGen(G, basis);
    HBGen.generate(trace.happens_before, basis.events_begin());
    blocked_hb += HBGen.blocked();

    bool realized_something = false;
    for (auto& hb :  HBGen.generated()) {
      realized_something |= tryRealizeAnnotations(trace, hb, basis);
    }

    return realized_something;
  }
}

void DCTraceBuilder::initializeHappensBefore(AnnotatedTrace& trace, Basis& basis) {
  // this gives some constrains induced by annotation
  std::vector<std::map<const llvm::Instruction *, unsigned>> instrs_to_index;

  for (auto it = basis.events_begin(), end = basis.events_end(); it != end; ++it) {
    // we want to work only with threads that are non-observatioanlly equivalent
    if (basis.isTopologyRoot(it.getProcessID())) {
      it.skipProcess();
      continue;
    }

    if (!isRead(**it))
      continue;

    const AnnotationValueT *A = trace.positive_annotation.get(**it);
    // if this read is not annotated or it is observed by the initial event
    if (!A || A->instruction == nullptr)
      continue;

    assert(isWrite(*A->instruction));

    // if the read is in the same process as the write that
    // it observes, we can skip it, since we do not add HB
    // for those
    if ((*it)->cpid == A->cpid)
      continue;

    if (basis.isTopologyRoot(A->cpid))
      continue;

    // if we already added HB from this write to some event in this thread,
    // then this is not the top-most read and we need to skip it
    auto HB = trace.happens_before.get(*A);
    if (HB && HB->count((*it)->cpid) != 0)
      continue;

    // We found annotated read that is top-most in its local process.
    // So now find the first conflicting event above
    // this read in this thread and make it happen-before the write from annotation
    // it must be the top-most read from this annotation,
    // because we are going from top to bottom
    unsigned i = it.getProcessID();
    unsigned j = it.getEventID();
    const DCEvent *found_confl = nullptr;
    while (j-- > 0) {
      if (basis[i][j]->ml.overlaps((*it)->ml)) {
        found_confl = basis[i][j];
        // we found the conflicting event, so if it is a read or write,
        // then use it (we check for reads or writes, because for example
        // pthread_mutex_init is conflicting with lock/unlock too)
        if (isWrite(*found_confl) || isRead(*found_confl))
          break;
        else
          // or reset the found event
          found_confl = nullptr;
      }
    }

    if (found_confl) {
      trace.happens_before.add({A->cpid, A->instruction, A->order},
                               {found_confl->cpid, found_confl->instruction, found_confl->order});
    } else {
      trace.happens_before.add({A->cpid, A->instruction, A->order},
                               {basis[i][0]->cpid, nullptr, 0});
    }
  }
}

bool DCTraceBuilder::tryRealizeAnnotations(AnnotatedTrace& trace,
                                           HappensExactlyBefore& happens_before,
                                           const Basis& basis)
{
  assert(!trace.has_root || basis.isTopologyRoot(trace.topology_root));
  assert(trace.has_root || basis.size() < 3);

  HappensAfterGraphVclock G(trace.positive_annotation,
                            happens_before, basis);

  ++realize_called;

  // add the new edges that we need to add in order
  // to be able to realize the positive annotation
  if (!G.realize()) {
    if (G.wasBlocked())
      ++blocked_before_sat;
    else
      ++blocked_after_sat;

#ifdef DUMP_ANNOT_TREE
     trace.annot_tree_node
         = annot_tree.addNode(getCurrentTrace().annot_tree_node,
                              trace.positive_annotation,
                              trace.negative_annotation,
                              happens_before,
                              "UNREALIZABLE");
#endif

    return false;
  }

  ++realize_succeeded;

  // need to do a copy of the trace, since we need the original
  // trace further for generating the traces (due to the happens-before)
  AnnotatedTrace new_trace;
  new_trace.positive_annotation = trace.positive_annotation;
  //new_trace.last_annotated_read = DCIID(trace.last_annotated_read);
  new_trace.negative_annotation
    = pruneNegativeAnnotation(trace.positive_annotation,
                              trace.negative_annotation);
  //new_trace.merged_pairs = trace.merged_pairs;

#ifdef DUMP_ANNOT_TREE
  new_trace.annot_tree_node
    = annot_tree.addNode(getCurrentTrace().annot_tree_node,
                         new_trace.positive_annotation,
                         new_trace.negative_annotation,
                         happens_before);
#endif

  // we won't need this happens before further, so we can move it
  new_trace.happens_before = std::move(happens_before);

  if (trace.has_root) {
    new_trace.has_root = true;
    new_trace.topology_root = trace.topology_root;
  }

  // DONT move the new trace -  it prevents compiler
  // from doing copy ellision
  //new_trace.trace = std::move(G.linearize());
  new_trace.trace = G.linearize();
  assert(!new_trace.trace.empty());

  // queue it for execution
  auto tr = extendTrace(std::move(new_trace.trace));
  if (has_error()) {
    return false;
  }

  new_trace.trace = std::move(tr);
  traces.emplace_back(std::move(new_trace));

  return true;
}

bool DCTraceBuilder::reset()
{
  // we have the initial trace in prefix
  AnnotatedTrace init_trace;

#ifdef DUMP_ANNOT_TREE
  annot_tree.createRoot(PositiveAnnotation::getObservationFunction(prefix));
#endif

  init_trace.trace = std::move(prefix);

  // recursively explore all possible traces
  explore(init_trace);

  // +1 because of the initial trace
  llvm::dbgs() << "Executed traces: " << executed_traces + 1 << "\n";
  llvm::dbgs() << "Explored observation functions num: " << leaves_number + 1
               << " (whole " << succ_leaves_number + 1 << ")\n";
  llvm::dbgs() << "Blocked due to past cone: " << blocked_past << "\n";
  llvm::dbgs() << "Blocked due to past cones (2): " << blocked_past2 << "\n";
  llvm::dbgs() << "Blocked before mutation: " << blocked_in_mutation << "\n";
  llvm::dbgs() << "Blocked happens-before recursions: " << blocked_hb << "\n";
  llvm::dbgs() << "Blocked before SAT reduction (excl. HB): " << blocked_before_sat << "\n";
  llvm::dbgs() << "Blocked in SAT reduction: " << blocked_after_sat << "\n";
  llvm::dbgs() << "Realized using direct swap: " << realized_using_swap << "\n";
  llvm::dbgs() << "Tried realizing annotations: " << realize_called << "("
               << realize_succeeded << " succeeded)\n";

#ifdef DUMP_ANNOT_TREE
    annot_tree.dump();
#endif

  return false;
}

template <typename TraceT>
static void dump_tr(TraceT& trace) {
  unsigned n = 0;
  for (const DCEvent& ev : trace) {
    assert(ev.id == n);
    llvm::errs() << n++ << ". " << ev << "\n";
  }
}

void AnnotatedTrace::dump() const {
    dump_tr(trace);
}

void DCTraceBuilder::dump_prefix() const {
  dump_tr(prefix);
}

void DCTraceBuilder::dump_current_trace() const {
  dump_tr(getCurrentTrace().trace);
}

void DCTraceBuilder::mayConflict(const ConstMRef *ml) {
  auto& curn = curnode();
  curn.may_conflict = true;
  if (ml)
    curn.ml = *ml;
  curn.instruction = current_inst;

  assert(getting_initial_trace || getting_maximal_extension ||
         (prefix[prefix_idx].instruction
            == getCurrentTrace().trace[prefix_idx].instruction));
  assert(getting_initial_trace || getting_maximal_extension ||
         getCurrentTrace().trace[prefix_idx].instruction == current_inst);
}

void DCTraceBuilder::atomic_store(const ConstMRef &ml)
{
  // stores to local memory on stack may not conflict
  assert(llvm::isa<llvm::StoreInst>(current_inst));
  if (!llvm::isa<llvm::AllocaInst>(
        current_inst->getOperand(1)->stripInBoundsOffsets()
      )) {
      mayConflict(&ml);
  }
}

void DCTraceBuilder::load(const ConstMRef &ml)
{
  // loads from stack may not conflict
  assert(llvm::isa<llvm::LoadInst>(current_inst));
  if (!llvm::isa<llvm::AllocaInst>(current_inst->getOperand(0)->stripInBoundsOffsets())) {
    mayConflict(&ml);
  }
}

void DCTraceBuilder::refuse_schedule()
{
  assert(prefix_idx == int(prefix.size())-1);
  assert(!prefix.back().may_conflict);

  // we already incremented p's clock
  unsigned p = prefix[prefix_idx].iid.get_pid();
  --threads[p].clock[p];

  if (curnode().size == 1) {
    prefix.pop_back();
    --prefix_idx;
  } else
    --curnode().size;

  assert(prefix_idx == int(prefix.size())-1);
}

void DCTraceBuilder::spawn()
{
  mayConflict();

  // store the CPid of the new thread
  IPid parent_ipid = curnode().iid.get_pid();
  curnode().childs_cpid = CPS.dry_spawn(threads[parent_ipid].cpid);

  CPid child_cpid = CPS.spawn(threads[parent_ipid].cpid);
  threads.emplace_back(child_cpid, threads[parent_ipid].clock);
  threads.emplace_back(CPS.new_aux(child_cpid), threads[parent_ipid].clock);
  threads.back().available = false; // Empty store buffer
}

void DCTraceBuilder::join(int tgt_proc){
  mayConflict();
  curnode().childs_cpid = threads[2*tgt_proc].cpid;
}

void DCTraceBuilder::mutex_lock(const ConstMRef &ml)
{
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.ref)){
    // Assume static initialization
    mutexes[ml.ref] = Mutex();
  }
  assert(mutexes.count(ml.ref));
  mayConflict(&ml);

  Mutex &mutex = mutexes[ml.ref];
  mutex.last_lock = mutex.last_access = prefix_idx;
  mutex.locked = true;
}

void DCTraceBuilder::mutex_trylock(const ConstMRef &ml){
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.ref)){
    // Assume static initialization
    mutexes[ml.ref] = Mutex();
  }
  assert(mutexes.count(ml.ref));
  mayConflict(&ml);

  Mutex &mutex = mutexes[ml.ref];
  mutex.last_access = prefix_idx;
  if(!mutex.locked){ // Mutex is free
    mutex.last_lock = prefix_idx;
    mutex.locked = true;
  }
}

void DCTraceBuilder::mutex_unlock(const ConstMRef &ml){
  assert(mutexes.count(ml.ref));

  fence();
  Mutex &mutex = mutexes[ml.ref];
  assert(0 <= mutex.last_access);

  mayConflict(&ml);

  mutex.last_access = prefix_idx;
  mutex.locked = false;
}

void DCTraceBuilder::mutex_init(const ConstMRef &ml){
  fence();
  assert(mutexes.count(ml.ref) == 0);
  curnode().may_conflict = true;
  curnode().ml = ml;
  curnode().instruction = current_inst;
  mutexes[ml.ref] = Mutex(prefix_idx);
}

void DCTraceBuilder::mutex_destroy(const ConstMRef &ml){
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.ref)){
    // Assume static initialization
    mutexes[ml.ref] = Mutex();
  }
  assert(mutexes.count(ml.ref));
  curnode().may_conflict = true;
  curnode().ml = ml;
  curnode().instruction = current_inst;

  mutexes.erase(ml.ref);
}

void DCTraceBuilder::metadata(const llvm::MDNode *md) {
  auto& cur = curnode();
  if (cur.md == nullptr)
    cur.md = md;

  last_md = md;
}

void DCTraceBuilder::fence(){
  IPid ipid = curnode().iid.get_pid();
  assert(ipid % 2 == 0);
}

IID<CPid> DCTraceBuilder::get_iid() const {
  IPid pid = curnode().iid.get_pid();
  int idx = curnode().iid.get_index();
  return IID<CPid>(threads[pid].cpid,idx);
}

Trace *DCTraceBuilder::get_trace() const {
  if (error_trace) {
    assert(errors.size() == 0);
    return error_trace;
  }

  std::vector<IID<CPid> > cmp;
  std::vector<const llvm::MDNode*> cmp_md;
  std::vector<Error*> errs;
  if (errors.size() == 0)
    return nullptr;

  for(unsigned i = 0; i < prefix.size(); ++i){
    cmp.push_back(IID<CPid>(threads[prefix[i].iid.get_pid()].cpid,prefix[i].iid.get_index()));
    cmp_md.push_back(prefix[i].md);
  };
  for(unsigned i = 0; i < errors.size(); ++i){
    errs.push_back(errors[i]->clone());
  }
  Trace *t = new IIDSeqTrace(cmp,cmp_md,errs);
  t->set_blocked(false);
  return t;
}

std::vector<DCEvent> DCTraceBuilder::extendCurrentTrace() {
  assert(prefix.empty());

  std::unique_ptr<llvm::ExecutionEngine> EE(DPORDriver::create_execution_engine(M, *this, config));

  // Run main.
  EE->runFunctionAsMain(M->getFunction("main"), {"prog"}, 0);

  // Run static destructors.
  EE->runStaticConstructorsDestructors(true);

  return prefix;
}

std::vector<DCEvent> DCTraceBuilder::extendTrace(std::vector<DCEvent>&& tr) {
    DCTraceBuilder TB(conf, M, std::move(tr));
    auto ret = TB.extendCurrentTrace();

    if (TB.has_error())
      error_trace = TB.get_trace();

    ++executed_traces;

    return ret;
}
