/* Copyright (C) 2014-2016 Carl Leonardsson
 * Copyright (C) 2016 Marek Chalupa

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

#ifndef __EVENT_H__
#define __EVENT_H__

#if defined(HAVE_LLVM_IR_METADATA_H)
#include <llvm/IR/Metadata.h>
#elif defined(HAVE_LLVM_METADATA_H)
#include <llvm/Metadata.h>
#endif

#include "VClock.h"
#include "vecset.h"

/* An identifier for a thread. An index into this->threads.
 *
 * Even indexes are for real threads. Odd indexes i are for
 * auxiliary threads corresponding to the real thread at index i-1.
 */
typedef int IPid;


/* Information about a (short) sequence of consecutive events by the
 * same thread. At most one event in the sequence may have conflicts
 * with other events, and if the sequence has a conflicting event,
 * it must be the first event in the sequence.
 */
class Event{
public:
  Event(const IID<IPid> &iid,
        const VClock<IPid> &clk)
    : iid(iid), origin_iid(iid), size(1), alt(0), md(0), clock(clk),
      may_conflict(false), sleep_branch_trace_count(0) {};

  /* A Branch object is a pair of an IPid p and an alternative index
   * (see Event::alt below) i. It will be tagged on an event in the
   * execution to indicate that if instead of that event, p is allowed
   * to execute (with alternative index i), then a different trace can
   * be produced.
   */
  class Branch{
  public:
    IPid pid;
    int alt;
    bool operator<(const Branch &b) const{
      return pid < b.pid || (pid == b.pid && alt < b.alt);
    };
    bool operator==(const Branch &b) const{
      return pid == b.pid && alt == b.alt;
    };
  };

  /* The identifier for the first event in this event sequence. */
  IID<IPid> iid;
  /* The IID of the program instruction which is the origin of this
   * event. For updates, this is the IID of the corresponding store
   * instruction. For other instructions origin_iid == iid.
   */
  IID<IPid> origin_iid;
  /* The number of events in this sequence. */
  int size;
  /* Some instructions may execute in several alternative ways
   * nondeterministically. (E.g. malloc may succeed or fail
   * nondeterministically if Configuration::malloy_may_fail is set.)
   * Event::alt is the index of the alternative for the first event
   * in this event sequence. The default execution alternative has
   * index 0. All events in this sequence, except the first, are
   * assumed to run their default execution alternative.
   */
  int alt;
  /* Metadata corresponding to the first event in this sequence. */
  const llvm::MDNode *md;
  /* The clock of the first event in this sequence. */
  VClock<IPid> clock;
  /* Is it possible for any event in this sequence to have a
   * conflict with another event?
   */
  bool may_conflict;
  /* Different, yet untried, branches that should be attempted from
   * this position in prefix.
   */
  VecSet<Branch> branch;
  /* The set of threads that go to sleep immediately before this
   * event sequence.
   */
  VecSet<IPid> sleep;
  /* The set of sleeping threads that wake up during or after this
   * event sequence.
   */
  VecSet<IPid> wakeup;
  /* For each previous IID that has been explored at this position
   * with the exact same prefix, some number of traces (both sleep
   * set blocked and otherwise) have been
   * explored. sleep_branch_trace_count is the total number of such
   * explored traces.
   */
  int sleep_branch_trace_count;
};

#endif
