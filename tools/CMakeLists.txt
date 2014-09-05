add_executable( printBCSR printBCSR.F90 )
add_executable( convertBCSR
  bcsr.cc
  bcsr.h
  convertBCSR.cc
  index.h
  lapack_interface.h
  logger.cc
  logger.h )
add_executable( plotDecay plotDecay.F90 )
