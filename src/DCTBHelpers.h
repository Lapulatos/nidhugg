#ifndef _DCTBHELPERS_H_
#define _DCTBHELPERS_H_

#include <llvm/IR/Instructions.h>
#include <llvm/IR/Function.h>
#include "Event.h"

/* helper functions */
inline bool is_function_call(const llvm::Instruction *I, const char *name)
{
  using namespace llvm;

  const CallInst *CI = dyn_cast<CallInst>(I);
  if(!CI)
    return false;

  const Function *F
    = dyn_cast<Function>(CI->getCalledValue()->stripPointerCasts());
  if(!F || !F->getName().equals(name))
    return false;

  return true;
}

inline bool isLock(const llvm::Instruction& I)
{
  return is_function_call(&I, "pthread_mutex_lock");
}

inline bool isLock(const DCEvent& ev)
{
  return isLock(*ev.instruction);
}

inline bool isUnlock(const llvm::Instruction& I)
{
  return is_function_call(&I, "pthread_mutex_unlock");
}

inline bool isUnlock(const DCEvent& ev)
{
  return isUnlock(*ev.instruction);
}

inline bool isRead(const llvm::Instruction& I)
{
  return llvm::isa<llvm::LoadInst>(&I) || isLock(I);
}

inline bool isRead(const DCEvent& ev)
{
  if (!ev.instruction)
    return false;
  return isRead(*ev.instruction);
}

inline bool isGlobalRead(const llvm::Instruction *I)
{
  // XXX: what about pointers!
  return llvm::isa<llvm::LoadInst>(I) &&
          !llvm::isa<llvm::AllocaInst>(I->getOperand(0)->stripInBoundsOffsets());
}

inline bool isGlobalRead(const DCEvent& ev)
{
  if (!ev.instruction)
    return false;
  return isGlobalRead(ev.instruction);
}

inline bool isWrite(const llvm::Instruction& I)
{
  return llvm::isa<llvm::StoreInst>(&I) || isUnlock(I);
}

inline bool isWrite(const DCEvent& ev)
{
  if (!ev.instruction)
    return false;
  return isWrite(*ev.instruction);
}

inline bool isJoin(const llvm::Instruction *I) {
  return is_function_call(I, "pthread_join");
}

inline bool isJoin(const DCEvent& ev) {
  if (!ev.instruction)
    return false;
  return is_function_call(ev.instruction, "pthread_join");
}

inline bool isPthreadCreate(const llvm::Instruction *I) {
  return is_function_call(I, "pthread_create");
}

inline bool isPthreadCreate(const DCEvent& ev) {
  if (!ev.instruction)
    return false;
  return is_function_call(ev.instruction, "pthread_create");
}

#endif // _DCTBHELPERS_H_

