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
#ifndef __TSO_TRACE_BUILDER_H__
#define __TSO_TRACE_BUILDER_H__

#include "TSOPSOTraceBuilder.h"
#include "VClock.h"

class TSOTraceBuilder : public TSOPSOTraceBuilder{
public:
  TSOTraceBuilder(const Configuration &conf = Configuration::default_conf);
  virtual ~TSOTraceBuilder();
  virtual bool schedule(int *proc, int *aux, int *alt, bool *dryrun);
  virtual void refuse_schedule();
  virtual void mark_available(int proc, int aux = -1);
  virtual void mark_unavailable(int proc, int aux = -1);
  virtual void cancel_replay();
  virtual bool is_replaying() const;
  virtual void metadata(const llvm::MDNode *md);
  virtual bool sleepset_is_empty() const;
  virtual bool check_for_cycles();
  virtual Trace *get_trace() const;
  virtual bool reset();
  virtual IID<CPid> get_iid() const;

  virtual void debug_print() const ;

  virtual void spawn();
  virtual void store(const ConstMRef &ml);
  virtual void atomic_store(const ConstMRef &ml);
  virtual void load(const ConstMRef &ml);
  virtual void full_memory_conflict();
  virtual void fence();
  virtual void join(int tgt_proc);
  virtual void mutex_lock(const ConstMRef &ml);
  virtual void mutex_lock_fail(const ConstMRef &ml);
  virtual void mutex_trylock(const ConstMRef &ml);
  virtual void mutex_unlock(const ConstMRef &ml);
  virtual void mutex_init(const ConstMRef &ml);
  virtual void mutex_destroy(const ConstMRef &ml);
  virtual bool cond_init(const ConstMRef &ml);
  virtual bool cond_signal(const ConstMRef &ml);
  virtual bool cond_broadcast(const ConstMRef &ml);
  virtual bool cond_wait(const ConstMRef &cond_ml, const ConstMRef &mutex_ml);
  virtual int cond_destroy(const ConstMRef &ml);
  virtual void register_alternatives(int alt_count);
  virtual int estimate_trace_count() const;
protected:
  /* A store pending in a store buffer. */
  class PendingStore{
  public:
    PendingStore(const ConstMRef &ml, const VClock<IPid> &clk, const llvm::MDNode *md)
      : ml(ml), clock(clk), last_rowe(-1), md(md) {};
    /* The memory location that is being written to. */
    ConstMRef ml;
    /* The clock of the store event which produced this store buffer
     * entry.
     */
    VClock<IPid> clock;
    /* An index into prefix to the event of the last load that fetched
     * its value from this store buffer entry by Read-Own-Write-Early.
     *
     * Has the value -1 if there has been no such load.
     */
    int last_rowe;
    /* "dbg" metadata for the store which produced this store buffer
     * entry.
     */
    const llvm::MDNode *md;
  };

  /* Various information about a thread in the current execution.
   */
  class Thread{
  public:
    Thread(const CPid &cpid, const VClock<IPid> &clk)
      : cpid(cpid), available(true), clock(clk), sleeping(false), sleep_full_memory_conflict(false) {};
    CPid cpid;
    /* Is the thread available for scheduling? */
    bool available;
    /* The clock containing all events that have been seen by this
     * thread.
     */
    VClock<IPid> clock;
    /* The store buffer of this thread. The store buffer is kept in
     * the Thread object for the real thread, not for the auxiliary.
     *
     * Newer entries are further to the back.
     */
    std::vector<PendingStore> store_buffer;
    /* True iff this thread is currently in the sleep set. */
    bool sleeping;
    /* sleep_accesses_r is the set of bytes that will be read by the
     * next event to be executed by this thread (as determined by dry
     * running).
     *
     * Empty if !sleeping.
     */
    VecSet<void const *> sleep_accesses_r;
    /* sleep_accesses_w is the set of bytes that will be written by
     * the next event to be executed by this thread (as determined by
     * dry running).
     *
     * Empty if !sleeping.
     */
    VecSet<void const *> sleep_accesses_w;
    /* sleep_full_memory_conflict is set when the next event to be
     * executed by this thread will be a full memory conflict (as
     * determined by dry running).
     */
    bool sleep_full_memory_conflict;
  };
  /* The threads in the current execution, in the order they were
   * created. Threads on even indexes are real, threads on odd indexes
   * i are the auxiliary threads corresponding to the real threads at
   * index i-1.
   */
  std::vector<Thread> threads;
  /* The CPids of threads in the current execution. */
  CPidSystem CPS;

  IPid ipid(int proc, int aux) const {
    assert(-1 <= aux && aux <= 0);
    assert(proc*2+1 < int(threads.size()));
    return aux ? proc*2 : proc*2+1;
  };

  std::string iid_string(const Event &evt) const;
  void add_branch(int i, int j);
  /* Add clocks and branches.
   *
   * All elements e in seen should either be indices into prefix, or
   * be negative. In the latter case they are ignored.
   */
  void see_events(const VecSet<int> &seen);
  /* Traverses prefix to compute the set of threads that were sleeping
   * as the first event of prefix[i] started executing. Returns that
   * set.
   */
  VecSet<IPid> sleep_set_at(int i);
  /* Wake up all threads which are sleeping, waiting for an access
   * (type,ml). */
  void wakeup(Access::Type type, void const *ml);
  /* Returns true iff the thread pid has a pending store to some
   * memory location including the byte ml.
   */
  bool has_pending_store(IPid pid, void const *ml) const;
  /* Helper for check_for_cycles. */
  bool has_cycle(IID<IPid> *loc) const;
  /* Estimate the total number of traces that have the same prefix as
   * the current one, up to the first idx events.
   */
  int estimate_trace_count(int idx) const;
};

#endif

