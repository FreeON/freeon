add_library( freeonmodules
  AInv.F90
  Anneal.F90
  AtomPairs.F90
  BasisSetParameters.F90
  BoundingBox.F90
  BraBloks.F90
  CholFact.F90
  Clock.F90
  D1Mach.F
  DCSEvl.F
  DErf.F
  DErfc.F
  DFRepository.F
  DMach.F
  DerivedTypes.F90
  FastMat.F90
  FileCopy.c
  FreeONSleep.c
  Functionals.F90
  GetAddress.c
  GetHostname.c
  GetLimits.c
  GetPID.c
  GetPWD.c
  GetRUsage.c
  GlobalCharacters.F90
  GlobalObjects.F90
  GlobalScalars.F90
  HDF5Wrap.c
  InOut.F90
  Indexing.F90
  InitDS.F
  LinAlg.F90
  LinAlg77.F
  Macros.F90
  MatFunk.F90
  McMurch.F90
  Mechanics.F90
  MemMan.F90
  MondoLogger.F90
  MondoLogger.h
  MondoMPI.F90
  Order.F90
  Order77.F
  PBC.F90
  Parse.F90
  ParsingConstants.F90
  PrettyPrint.F90
  ProcessControl.F90
  Random.c
  SMWrap.c
  Seconds.c
  SetEq.F90
  SetEq77.F
  Slatec.F90
  SpecFun.F90
  TempDirectory.c
  Thresholding.F90
  TimerWrapper.c
  Trap.c
  UnlimitStack.c
  Utilities.F90
  mm_double_NN_1.c
  mm_double_NN_c.c
  mm_double_NT_1.c
  mm_double_NT_c.c
  mm_double_TN_1.c
  mm_double_TN_c.c
  mm_double_TT_1.c
  mm_double_TT_c.c
  phipac_dgemm.c
  )

include_directories( ${CMAKE_BINARY_DIR} )
