/* Copyright (C) 2014-2016 Carl Leonardsson
 * Copyright (C) 2016 Marek Chalupa
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

DCTraceBuilder::DCTraceBuilder(const Configuration &conf) : TSOTraceBuilder(conf)
{
}

DCTraceBuilder::~DCTraceBuilder()
{
}

bool DCTraceBuilder::schedule_replay(int *proc, int *aux, int *alt, bool *dryrun)
{
  /* Are we done with the current Event? */
  if(0 <= prefix_idx &&
     threads[curnode().iid.get_pid()].clock[curnode().iid.get_pid()] <
     curnode().iid.get_index() + curnode().size - 1){
    /* Continue executing the current Event */
    IPid pid = curnode().iid.get_pid();
    *proc = pid/2;
    *aux = pid % 2 - 1;
    *alt = 0;
    assert(threads[pid].available);
    ++threads[pid].clock[pid];
    return true;
  }else if(prefix_idx + 1 == int(prefix.size())){
    /* We are done replaying. Continue below... */
    replay = false;
    return false;
  }else if(dry_sleepers < int(prefix[prefix_idx+1].sleep.size())){
    /* Before going to the next event, dry run the threads that are
     * being put to sleep.
     */
    IPid pid = prefix[prefix_idx+1].sleep[dry_sleepers];
    ++dry_sleepers;
    threads[pid].sleeping = true;
    *proc = pid/2;
    *aux = pid % 2 - 1;
    *alt = 0;
    *dryrun = true;
    this->dryrun = true;
    return true;
  }else{
    /* Go to the next event. */
    dry_sleepers = 0;
    ++prefix_idx;
    IPid pid = curnode().iid.get_pid();
    *proc = pid/2;
    *aux = pid % 2 - 1;
    *alt = curnode().alt;
    assert(threads[pid].available);
    ++threads[pid].clock[pid];
    curnode().clock = threads[pid].clock;
    return true;
  }
}

// schedule a next event to be the next event from thread p
bool DCTraceBuilder::schedule_thread(int *proc, unsigned p)
{

  if(threads[p].available && !threads[p].sleeping &&
     (conf.max_search_depth < 0
      || threads[p].clock[p] < conf.max_search_depth)) {
    // increase clock
    ++threads[p].clock[p];
    // extend the prefix
    prefix.push_back(Event(IID<IPid>(IPid(p),threads[p].clock[p]),
                           threads[p].clock));
    // set the scheduled thread
    *proc = p/2;

    return true;
  }

  return false;
}

void DCTraceBuilder::squeezeLastEvent()
{
    assert(prefix[prefix.size()-1].branch.empty());
    assert(prefix[prefix.size()-1].wakeup.empty());
    ++prefix[prefix.size()-2].size;
    prefix.pop_back();
    --prefix_idx;
}

void DCTraceBuilder::mergeInvisibleEvents()
{
  if(prefix.size() > 1 &&
     prefix[prefix.size()-1].iid.get_pid() == prefix[prefix.size()-2].iid.get_pid() &&
     !prefix[prefix.size()-1].may_conflict && prefix[prefix.size()-1].sleep.empty())
    squeezeLastEvent();
}

bool DCTraceBuilder::schedule(int *proc, int *aux, int *alt, bool *dryrun)
{
  *dryrun = false;
  *alt = 0;
  this->dryrun = false;

  if(replay){
    if (schedule_replay(proc, aux, alt, dryrun))
      return true;
    // else continue in scheduling,
    // since we have scheduled nothing
  }

  assert(!replay);

  /* Merge previous invisible events if needed */
  mergeInvisibleEvents();

  /* Create a new Event */
  ++prefix_idx;
  assert(prefix_idx == int(prefix.size()));

  /* Find an available thread (auxiliary or real).
   *
   * Prioritize auxiliary before real, and older before younger
   * threads.
   */
  const unsigned sz = threads.size();
  unsigned p;

  // Loop through auxiliary threads
  for(p = 1; p < sz; p += 2){
    if (schedule_thread(proc, p)) {
      // auxiliary thread
      *aux = 0;
      return true;
    }
  }

  // Loop through real threads
  for(p = 0; p < sz; p += 2){
    if (schedule_thread(proc, p)) {
      *aux = -1;
      return true;
    }
  }

  return false; // No available threads
}

