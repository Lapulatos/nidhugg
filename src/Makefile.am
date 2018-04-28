AUTOMAKE_OPTIONS = -Wno-override ## Make autotools quit complaining about explicit rule for nidhuggc.

noinst_LIBRARIES = libnidhugg.a
libnidhugg_a_SOURCES = \
  AddLibPass.cpp AddLibPass.h \
  BVClock.cpp BVClock.h \
  CheckModule.cpp CheckModule.h \
  Configuration.cpp Configuration.h \
  CPid.cpp CPid.h \
  Debug.cpp Debug.h \
  DetCheckTraceBuilder.cpp DetCheckTraceBuilder.h \
  DPORDriver.cpp DPORDriver.h \
  Execution.cpp \
  ExternalFunctions.cpp \
  FBVClock.cpp FBVClock.h \
  GlobalContext.cpp GlobalContext.h \
  IID.h IID.tcc \
  Interpreter.cpp Interpreter.h \
  LoopUnrollPass.cpp LoopUnrollPass.h \
  MRef.cpp MRef.h \
  nregex.cpp nregex.h \
  POWERExecution.cpp \
  POWERInterpreter.cpp POWERInterpreter.h \
  POWERARMTraceBuilder.cpp POWERARMTraceBuilder.tcc POWERARMTraceBuilder.h POWERARMTraceBuilder_decl.h \
  PSOInterpreter.cpp PSOInterpreter.h \
  PSOTraceBuilder.cpp PSOTraceBuilder.h \
  SigSegvHandler.cpp SigSegvHandler.h \
  SpinAssumePass.cpp SpinAssumePass.h \
  StrModule.cpp StrModule.h \
  Trace.cpp Trace.h Event.h\
  TraceBuilder.cpp TraceBuilder.h \
  Transform.cpp Transform.h \
  TSOInterpreter.cpp TSOInterpreter.h \
  DCInterpreter.cpp DCInterpreter.h \
  TSOPSOTraceBuilder.h \
  TSOTraceBuilder.cpp TSOTraceBuilder.h \
  DCTraceBuilder.cpp DCTraceBuilder.h \
  VClock.cpp VClock.h VClock.tcc \
  Basis.h Basis.cpp \
  Annotations.h Annotations.cpp \
  AnnotatedTrace.h \
  2SAT.h 2SAT.cpp \
  HappensAfterGraph.h HappensAfterGraph.cpp \
  HappensAfterGraphVclock.h HappensAfterGraphVclock.cpp \
  HappensExactlyBefore.h HappensExactlyBefore.cpp \
  vecset.h vecset.tcc
libnidhugg_a_CXXFLAGS = -fno-rtti

EXTRA_PROGRAMS = nidhuggc
bin_PROGRAMS = nidhugg @NIDHUGGCBIN@
nidhugg_SOURCES = main.cpp
nidhugg_LDADD = libnidhugg.a
nidhugg_CXXFLAGS = -fno-rtti
nidhuggc_SOURCES = nidhuggc.py

nidhuggc$(EXEEXT): $(srcdir)/nidhuggc.py
	cat $(srcdir)/nidhuggc.py \
	| sed 's|%%PYTHON%%|@PYTHON@|g' \
	| sed 's|%%CLANG%%|@CLANG@|g' \
	| sed 's|%%CLANGXX%%|@CLANGXX@|g' \
	| sed 's|%%NIDHUGGPATH%%|$(bindir)/nidhugg|g' \
	> nidhuggc
	chmod a+x nidhuggc
