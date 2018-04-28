#include <set>

#include "Event.h"
#include "Annotations.h"
#include "Basis.h"

void Basis::dump() const
{
  llvm::errs() << "Basis {\n";
  for (auto& process : processes) {
    llvm::errs() << process[0]->cpid << ":\n";
    for (const DCEvent *ev : process) {
      llvm::errs() << "  " << *ev << "\n";
    }
  }
}
