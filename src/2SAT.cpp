#include <llvm/Support/raw_ostream.h>
#include <llvm/IR/Instruction.h>
#include <llvm/IR/Instructions.h>

#include "2SAT.h"
#include "SCC.h"

bool Events2SAT::solve()
{
  // compute strongly connected components
  SCC<Node> scc_comp(nodes);
  auto& SCCs = scc_comp.compute();

  // process the components in the reverse topological order:
  // (the components are stored in the reverse topological
  // order in the vector)
  for (auto& scc : SCCs) {
    std::vector<Info *> scc_infos;
    scc_infos.reserve(scc.size());
    // gather the information about nodes we have
    // in the SCC
    for (unsigned i : scc) {
      Node& nd = nodes[i];
      scc_infos.push_back(&infos[nd.e1][nd.e2]);
    }

    int scc_valuation = -1;
    for (unsigned i = 0; i < scc.size(); ++i) {
      unsigned ndidx = scc[i];
      Info *ndinf = scc_infos[i];
      bool negative = ndidx % 2;

      // get the index of the positive literal
      // - it is either ndidx or the index before that
      // if ndidx is odd
      unsigned posidx = negative ? ndidx - 1 : ndidx;
      assert(posidx % 2 == 0);
      assert(ndinf->index == (int) posidx);

      // a node is in the same SCC as its complement?
      if (nodes[posidx].scc_id == nodes[posidx + 1].scc_id)
        return false;

      // we already have a valuation of some node
      // in the SCC? Then it must agree with the valuation
      // of all other nodes in the SCC
      if (ndinf->valuation != -1) {
        unsigned ndval = negative ? !ndinf->valuation : ndinf->valuation;
        if (scc_valuation == -1)
          // in the info, we have that the positivie literal is
          // false or true, therefore if the component contains
          // the negative literal, that means that we must negate
          // the valuation to get the right value that is going to
          // be stored in the info
          scc_valuation = ndval;
        else if (scc_valuation != (int) ndval)
          return false;
      }
    }

    // we got here? - then all the nodes are either
    // unevaluated or the evaluation agrees, so this component
    // is fine. Set the evaluation to all the nodes (and thus
    // implicitly to the complements of nodes)
    if (scc_valuation == -1)
      scc_valuation = 1;

    for (unsigned i = 0; i < scc.size(); ++i) {
      unsigned ndidx = scc[i];
      Info *ndinf = scc_infos[i];
      bool negative = ndidx % 2;

      unsigned ndval = negative ? !scc_valuation : scc_valuation;
      assert(ndinf->valuation == -1 || ndinf->valuation == (int) ndval);

      ndinf->valuation = ndval;
    }
  }

  // if we got here, the formula is satisfiable
  return true;
}

void Events2SAT::to_dot() const
{
  llvm::errs() << "digraph {\n";
  for (unsigned i = 0; i < nodes.size(); i += 2) {
  //for (auto I = nodes_begin(), E = nodes_end(); I!=E; ++I) {
  // const Node *pos = std::get<0>(*I);
  // const Node *neg = std::get<1>(*I);
  // const Info *inf = std::get<2>(*I);
    const Node *pos = &nodes[i];
    const Node *neg = &nodes[i+1];
    assert(pos->e1 == neg->e1 && pos->e2 == neg->e2);
    const Info *inf = getInfo(pos->e1, pos->e2);

    llvm::errs() << "NODE" << pos << " [label=\"" << inf->index << "\\n";
    llvm::errs() << "scc: " << pos->get_scc_id() << "\\n";
    if (pos->e1) {
      llvm::errs() << pos->e1->cpid;
      if (pos->e1->instruction)
        llvm::errs() << *pos->e1->instruction;
      else
        llvm::errs() << "(null)";
    } else
      llvm::errs() << "initial";

    llvm::errs() << "\\n-> ";

    if (pos->e2) {
      llvm::errs() << pos->e2->cpid;
      if (pos->e2->instruction)
        llvm::errs() << *pos->e2->instruction;
      else
        llvm::errs() << "(null)";
    } else
      llvm::errs() << "initial";

    llvm::errs() << "\\nval: " << inf->valuation;
    llvm::errs() << "\"]\n";

    llvm::errs() << "NODE" << neg << " [label=\"" << inf->index << " NOT\\n";
    llvm::errs() << "scc: " << neg->get_scc_id() << "\\n";
    if (neg->e1) {
      llvm::errs() << neg->e1->cpid;
      if (neg->e1->instruction)
        llvm::errs() << *neg->e1->instruction;
      else
        llvm::errs() << "(null)";
    } else
      llvm::errs() << "initial";

    llvm::errs() << "\\n-> ";
    if (neg->e2) {
      llvm::errs() << neg->e2->cpid;
      if (neg->e2->instruction)
        llvm::errs() << *neg->e2->instruction;
      else
        llvm::errs() << "(null)";
    } else
      llvm::errs() << "initial";

    if (inf->valuation == -1)
      llvm::errs() << "\\nval: " <<  -1;
    else
      llvm::errs() << "\\nval: " << (int)!inf->valuation;

    llvm::errs() << "\"]\n";
  }

  llvm::errs() << "\n";
  //for (auto I = nodes_begin(), E = nodes_end(); I!=E; ++I) {
  // const Node *pos = std::get<0>(*I);
  // const Node *neg = std::get<1>(*I);
  for (unsigned i = 0; i < nodes.size(); i += 2) {
    const Node *pos = &nodes[i];
    const Node *neg = &nodes[i+1];
    assert(pos->e1 == neg->e1 && pos->e2 == neg->e2);

    for (unsigned s : pos->getSuccessors()) {
      llvm::errs() << "NODE" << pos << " -> NODE" << getNode(s) << "\n";
    }
    for (unsigned s : neg->getSuccessors()) {
      llvm::errs() << "NODE" << neg << " -> NODE" << getNode(s) << "\n";
    }
  }

  llvm::errs() << "}\n";
  llvm::errs().flush();
}

