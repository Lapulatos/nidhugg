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
#ifndef __PSO_TRACE_BUILDER_H__
#define __PSO_TRACE_BUILDER_H__

#include "TSOPSOTraceBuilder.h"

class PSOTraceBuilder : public TSOPSOTraceBuilder{
public:
  PSOTraceBuilder(const Configuration &conf = Configuration::default_conf);
  virtual ~PSOTraceBuilder();
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
  /* A byte of a store pending in a store buffer. */
  class PendingStoreByte{
  public:
    PendingStoreByte(const ConstMRef &ml, const VClock<IPid> &clk, const llvm::MDNode *md)
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
    Thread(int proc, const CPid &cpid, const VClock<IPid> &clk, IPid parent)
      : proc(proc), cpid(cpid), available(true), clock(clk), parent(parent),
        sleeping(false), sleep_full_memory_conflict(false) {};
    /* The process number used to communicate process identities with
     * the interpreter. For auxiliary threads, proc is the process
     * number of the corresponding real thread.
     */
    int proc;
    CPid cpid;
    /* Is the thread available for scheduling? */
    bool available;
    /* The clock containing all events that have been seen by this
     * thread.
     */
    VClock<IPid> clock;
    /* The IPid of the parent thread of this thread. "Parent" is here
     * interpreted in the same way as in the context of CPids.
     */
    IPid parent;
    /* aux_to_byte and byte_to_aux provide the mapping between
     * auxiliary thread indices and the first byte in the memory
     * locations for which that auxiliary thread is responsible.
     */
    std::vector<void const*> aux_to_byte;
    std::map<void const*,int> byte_to_aux;
    /* For each auxiliary thread i of this thread, aux_to_ipid[i] is
     * the corresponding IPid.
     */
    std::vector<IPid> aux_to_ipid;
    /* Each pending store s is split into PendingStoreBytes (see
     * above) and ordered into store buffers. For each byte b in
     * memory, store_buffers[b] contains precisely all
     * PendingStoreBytes of pending stores to that byte. The entries
     * in store_buffers[b] are ordered such that newer entries are
     * further to the back.
     *
     * The store buffer is kept in the Thread object for the real
     * thread, not for the auxiliary.
     */
    std::map<void const*,std::vector<PendingStoreByte> > store_buffers;
    /* For a non-auxiliary thread, aux_clock_sum is the sum of the
     * clocks of all auxiliary threads belonging to this thread.
     */
    VClock<IPid> aux_clock_sum;
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

    bool all_buffers_empty() const{
#ifndef NDEBUG
      /* Empty buffers should be removed from store_buffers. */
      for(auto it = store_buffers.begin(); it != store_buffers.end(); ++it){
        assert(it->second.size());
      }
#endif
      return store_buffers.empty();
    };
  };
  /* The threads in the current execution, in the order they were
   * created. Threads on even indexes are real, threads on odd indexes
   * i are the auxiliary threads corresponding to the real threads at
   * index i-1.
   */
  std::vector<Thread> threads;
  /* A map from the process numbers used in communication with the
   * interpreter to IPids.
   */
  std::vector<IPid> proc_to_ipid;
  /* The CPids of threads in the current execution. */
  CPidSystem CPS;
  /* The set sleepers contains precisely the IPids p such that
   * threads[p].sleeping is true.
   */
  VecSet<IPid> sleepers;
  /* The set available_threads contains precisely the IPids p such
   * that threads[p].cpid.is_auxiliary is false and
   * threads[p].available is true.
   */
  VecSet<IPid> available_threads;
  /* The set available_auxs contains precisely the IPids p such that
   * threads[p].cpid.is_auxiliary is true and threads[p].available is
   * true.
   */
  VecSet<IPid> available_auxs;

  IPid ipid(int proc, int aux) const {
    assert(0 <= proc && proc < int(proc_to_ipid.size()));
    IPid p = proc_to_ipid[proc];
    if(0 <= aux){
      assert(aux < int(threads[p].aux_to_ipid.size()));
      p = threads[p].aux_to_ipid[aux];
    }
    return p;
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
  /* Returns true when it is possible to perform a memory update from
   * the auxiliary thread pid. I.e., there should be a pending store s
   * in the buffers of the real thread corresponding to pid, the first
   * accessed byte of which is the one for which pid is
   * responsible. Furthermore, there should be no other pending stores
   * which overlap with s and which is ordered before s in the store
   * buffers.
   *
   * Pre: threads[pid].cpid.is_auxiliary()
   */
  bool is_aux_at_head(IPid pid) const;
  /* Estimate the total number of traces that have the same prefix as
   * the current one, up to the first idx events.
   */
  int estimate_trace_count(int idx) const;
  /* Same as mark_available, but takes an IPid as thread
   * identifier.
   */
  void mark_available_ipid(IPid pid);
  /* Same as mark_unavailable, but takes an IPid as thread
   * identifier.
   */
  void mark_unavailable_ipid(IPid pid);
};

#endif

