/* Copyright (C) 2014-2017 Carl Leonardsson
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

#ifndef __TSO_PSO_TRACE_BUILDER_H__
#define __TSO_PSO_TRACE_BUILDER_H__

#include "Configuration.h"
#include "MRef.h"
#include "Trace.h"
#include "DetCheckTraceBuilder.h"
#include "Event.h"

#include <string>
#include <vector>

/*  === Thread Identification ===
 *
 * Througout this class, we will use the following convention for
 * identifying threads:
 *
 * A thread is identified by a pair (proc,aux). For real threads, aux
 * == -1. For auxiliary threads, 0 <= aux. The thread (0,-1) is the
 * single original thread in each execution. The thread (proc,-1) is
 * the proc:th thread that was created during the execution. The
 * thread (proc,aux) for 0 <= aux is an auxiliary thread of
 * (proc,-1). The exact interpretation of an auxiliary thread is
 * memory model specific.
 */

class TSOPSOTraceBuilder : public DetCheckTraceBuilder{
public:
  TSOPSOTraceBuilder(const Configuration &conf = Configuration::default_conf)
    : DetCheckTraceBuilder(conf) {};
  virtual ~TSOPSOTraceBuilder() {};
  /* Schedules the next thread.
   *
   * Sets *proc and *aux, such that (*proc,*aux) should be the next
   * thread to run.
   *
   * *alt is set to indicate which alternative execution should be
   * used for the next event, in cases where events have several
   * nondeterministic alternatives for execution. *alt == 0 means the
   * first alternative (the only alternative for deterministic
   * events). Higher values means other alternatives.
   *
   * *dryrun is set to indicate that the next event should only be dry
   * run. Dry running an event means simulating its effects without
   * making any changes to the program state, and without considering
   * the event to have been executed. During dry run, the TraceBuilder
   * will collect information about the event through calls to the
   * methods spawn, store, load, ..., which are used to register
   * aspect that affect the chronological trace. This is used to infer
   * when sleeping threads should be awoken.
   *
   * true is returned on successful scheduling. false is returned if
   * no thread can be scheduled (because all threads are marked
   * unavailable).
   */
  virtual bool schedule(int *proc, int *aux, int *alt, bool *dryrun) = 0;
  /* Reject the last scheduling. Erase the corresponding event and its
   * effects on the chronological trace. Also mark the last scheduled
   * thread as unavailable.
   *
   * Can only be used after a successful call to schedule(_,_,_,d)
   * where *d is false. Multiple refusals without intermediate
   * schedulings are not allowed.
   *
   * Use this for example in cases when the scheduled event turns out
   * to be blocked, waiting for e.g. store buffer flushes or the
   * termination of another thread.
   */
  virtual void refuse_schedule() = 0;
  /* Notify the TraceBuilder that (proc,aux) is available for
   * scheduling.
   */
  virtual void mark_available(int proc, int aux = -1) = 0;
  /* Notify the TraceBuilder that (proc,aux) is unavailable for
   * scheduling. This may be e.g. because it has terminated, or
   * because it is blocked waiting for something.
   */
  virtual void mark_unavailable(int proc, int aux = -1) = 0;
  /* If we are not in a replay, do nothing. Otherwise cancel the
   * replay from here on, so that the computation may continue
   * according to an arbitrary schedule.
   */
  virtual void cancel_replay() = 0;
  /* Associate the currently scheduled event with LLVM "dbg" metadata. */
  virtual void metadata(const llvm::MDNode *md) = 0;

  /*******************************************/
  /*           Registration Methods          */
  /*******************************************/
  /* The following methods are used to register information about the
   * currently scheduled event, which may affect the chronological
   * trace.
   */

  /* The current event spawned a new thread. */
  virtual void spawn() = 0;
  /* Perform a store to ml. */
  virtual void store(const ConstMRef &ml) = 0;
  /* Perform an atomic store to ml.
   *
   * The exact interpretation depends on the memory model. But
   * examples would typically include instructions such as "store
   * atomic, ..., seq_cst" or "cmpxchg ... seq_cst ...".
   *
   * atomic_store performed by auxiliary threads can be used to model
   * memory updates under non-SC memory models.
   */
  virtual void atomic_store(const ConstMRef &ml) = 0;
  /* Perform a load to ml. */
  virtual void load(const ConstMRef &ml) = 0;
  /* Perform an action that conflicts with all memory accesses and all
   * other full memory conflicts.
   *
   * This is used for things that execute as blackboxes (e.g. calls to
   * unknown external functions) in order to capture all potential
   * conflicts that such blackboxes may cause.
   */
  virtual void full_memory_conflict() = 0;
  /* Perform a full memory fence. */
  virtual void fence() = 0;
  /* Perform a pthread_join, with the thread (tgt_proc,-1). */
  virtual void join(int tgt_proc) = 0;
  /* Perform a successful pthread_mutex_lock to the pthread mutex at
   * ml.
   */
  virtual void mutex_lock(const ConstMRef &ml) = 0;
  /* Perform a failed attempt at pthread_mutex_lock to the pthread
   * mutex at ml.
   */
  virtual void mutex_lock_fail(const ConstMRef &ml) = 0;
  /* Perform a pthread_mutex_trylock (successful or failing) to the
   * pthread mutex at ml.
   */
  virtual void mutex_trylock(const ConstMRef &ml) = 0;
  /* Perform a pthread_mutex_unlock to the pthread mutex at ml. */
  virtual void mutex_unlock(const ConstMRef &ml) = 0;
  /* Initialize a pthread mutex at ml. */
  virtual void mutex_init(const ConstMRef &ml) = 0;
  /* Destroy a pthread mutex at ml. */
  virtual void mutex_destroy(const ConstMRef &ml) = 0;
  /* Initialize a pthread condition variable at ml.
   *
   * Returns true on success, false if a pthreads_error has been
   * generated.
   */
  virtual bool cond_init(const ConstMRef &ml) = 0;
  /* Signal on the condition variable at ml.
   *
   * Returns true on success, false if a pthreads_error has been
   * generated.
   */
  virtual bool cond_signal(const ConstMRef &ml) = 0;
  /* Broadcast on the condition variable at ml.
   *
   * Returns true on success, false if a pthreads_error has been
   * generated.
   */
  virtual bool cond_broadcast(const ConstMRef &ml) = 0;
  /* Wait for a condition variable at cond_ml.
   *
   * The mutex at mutex_ml is used as the mutex for pthread_cond_wait.
   *
   * Returns true on success, false if a pthreads_error has been
   * generated.
   */
  virtual bool cond_wait(const ConstMRef &cond_ml, const ConstMRef &mutex_ml) = 0;
  /* Destroy a pthread condition variable at ml.
   *
   * Returns 0 on success, EBUSY if some thread was waiting for the
   * condition variable, another value if a pthreads_error has been
   * generated.
   */
  virtual int cond_destroy(const ConstMRef &ml) = 0;
  /* Notify TraceBuilder that the current event may
   * nondeterministically execute in several alternative ways. The
   * number of ways is given in alt_count.
   */
  virtual void register_alternatives(int alt_count) = 0;
protected:
  /* The fixed prefix of events in the current execution. This may be
   * either the complete sequence of events executed thus far in the
   * execution, or the events executed followed by the subsequent
   * events that are determined in advance to be executed.
   */
  std::vector<Event> prefix;

  /* The index into prefix corresponding to the last event that was
   * scheduled. Has the value -1 when no events have been scheduled.
   */
  int prefix_idx;

  /* An Access is a pair (tp,ml) representing an access to
   * memory. Accesses come in two varieties:
   *
   * (W_ALL_MEMORY,0) is considered as a store access to all memory.
   *
   * (tp,ml) with tp in {R,W}, is a Read or Write access to the byte
   * indicated by the pointer ml.
   */
  class Access{
  public:
    /* The type of memory access. */
    enum Type {R, W, W_ALL_MEMORY, NA} type;
    /* The accessed byte. */
    const void *ml;
    bool operator<(const Access &a) const{
      return type < a.type || (type == a.type && ml < a.ml);
    };
    bool operator==(const Access &a) const{
      return type == a.type && (type == NA || ml == a.ml);
    };
    Access() : type(NA), ml(0) {};
    Access(Type t, const void *m) : type(t), ml(m) {};
  };

  /* A Mutex represents a pthread_mutex_t object.
   */
  class Mutex{
  public:
    Mutex() : last_access(-1), last_lock(-1) {};
    Mutex(int lacc) : last_access(lacc), last_lock(-1) {};
    int last_access;
    int last_lock;
  };
  /* A map containing all pthread mutex objects in the current
   * execution. The key is the position in memory of the actual
   * pthread_mutex_t object.
   */
  std::map<void const*,Mutex> mutexes;

  /* A CondVar represents a pthread_cond_t object. */
  class CondVar{
  public:
    CondVar() : last_signal(-1) {};
    CondVar(int init_idx) : last_signal(init_idx) {};
    /* Index in prefix of the latest call to either of
     * pthread_cond_init, pthread_cond_signal, or
     * pthread_cond_broadcast for this condition variable.
     *
     * -1 if there has been no such call.
     */
    int last_signal;
    /* For each thread which is currently waiting for this condition
     * variable, waiters contains the index into prefix of the event
     * where the thread called pthread_cond_wait and started waiting.
     */
    std::vector<int> waiters;
  };
  /* A map containing all pthread condition variable objects in the
   * current execution. The key is the position in memory of the
   * actual pthread_cond_t object.
   */
  std::map<void const*,CondVar> cond_vars;

  /* The number of threads that have been dry run since the last
   * non-dry run event was scheduled.
   */
  int dry_sleepers;

  /* Are we currently executing an event in dry run mode? */
  bool dryrun;

  /* Are we currently replaying the events given in prefix from the
   * previous execution? Or are we executing new events by arbitrary
   * scheduling?
   */
  bool replay;

  /* The latest value passed to this->metadata(). */
  const llvm::MDNode *last_md;

  /* A ByteInfo object contains information about one byte in
   * memory. In particular, it recalls which events have recently
   * accessed that byte.
   */
  class ByteInfo{
  public:
    ByteInfo() : last_update(-1), last_update_ml(0,1) {};
    /* An index into prefix, to the latest update that accessed this
     * byte. last_update == -1 if there has been no update to this
     * byte.
     */
    int last_update;
    /* The complete memory location (possibly multiple bytes) that was
     * accessed by the last update. Undefined if there has been no
     * update to this byte.
     */
    ConstMRef last_update_ml;
    /* last_read[tid] is the index in prefix of the latest (visible)
     * read of thread tid to this memory location, or -1 if thread tid
     * has not read this memory location.
     *
     * The indexing counts real threads only, so e.g. last_read[1] is
     * the last read of the second real thread.
     *
     * last_read_t is simply a wrapper around a vector, which expands
     * the vector as necessary to accomodate accesses through
     * operator[].
     */
    struct last_read_t {
      std::vector<int> v;
      int operator[](int i) const { return (i < int(v.size()) ? v[i] : -1); };
      int &operator[](int i) {
        if(int(v.size()) <= i){
          v.resize(i+1,-1);
        }
        return v[i];
      };
      std::vector<int>::iterator begin() { return v.begin(); };
      std::vector<int>::const_iterator begin() const { return v.begin(); };
      std::vector<int>::iterator end() { return v.end(); };
      std::vector<int>::const_iterator end() const { return v.end(); };
    } last_read;
  };
  std::map<const void*,ByteInfo> mem;
  /* Index into prefix pointing to the latest full memory conflict.
   * -1 if there has been no full memory conflict.
   */
  int last_full_memory_conflict;

  Event &curnode() {
    assert(0 <= prefix_idx);
    assert(prefix_idx < int(prefix.size()));
    return prefix[prefix_idx];
  };

  const Event &curnode() const {
    assert(0 <= prefix_idx);
    assert(prefix_idx < int(prefix.size()));
    return prefix[prefix_idx];
  };
};

#endif
