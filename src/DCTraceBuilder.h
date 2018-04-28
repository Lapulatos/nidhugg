/* Copyright (C) 2014-2016 Carl Leonardsson
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

#include "TSOTraceBuilder.h"
#include "VClock.h"

class DCTraceBuilder : public TSOTraceBuilder{
  bool schedule_reply(int *proc, int *aux, int *alt, bool *dryrun);

public:
  DCTraceBuilder(const Configuration &conf = Configuration::default_conf);
  virtual ~DCTraceBuilder();

  virtual bool schedule(int *proc, int *aux, int *alt, bool *dryrun);
/*
  virtual void refuse_schedule();
  virtual void mark_available(int proc, int aux = -1);
  virtual void mark_unavailable(int proc, int aux = -1);
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
*/

};

#endif

