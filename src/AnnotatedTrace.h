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
#ifndef __DC_ANNOTATED_TRACE_H__
#define __DC_ANNOTATED_TRACE_H__

#include "Annotations.h"
#include "HappensExactlyBefore.h"

//#define DUMP_ANNOT_TREE 1

// AnnotatedTrace is a trace along with positive and
// negative annotation of the events
struct AnnotatedTrace {
  PositiveAnnotation positive_annotation;
  NegativeAnnotation negative_annotation;
  // happens before relation for cyclic architectures
  HappensExactlyBefore happens_before;

  // a root of topology for handling acyclic architectures,
  // once we fix it, it must stay the same for all traces derived
  // from this one
  bool has_root = false;
  CPid topology_root;

#ifdef DUMP_ANNOT_TREE
  unsigned annot_tree_node = 0;
#endif

  std::vector<DCEvent> trace;

  AnnotatedTrace() = default;

  AnnotatedTrace(const PositiveAnnotation& pa,
                 const NegativeAnnotation& na,
                 const HappensExactlyBefore& heb)
      : positive_annotation(pa), negative_annotation(na),
        happens_before(heb) {}

  AnnotatedTrace(AnnotatedTrace&& tr) = default;
  AnnotatedTrace& operator=(AnnotatedTrace&& tr) = default;

  void dump() const;
};

#endif

