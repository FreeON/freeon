PROGRAM EllToInclude
   USE DerivedTypes
   USE GlobalScalars
#ifdef NAG
   USE F90_UNIX
#endif
   CHARACTER(LEN=DEFAULT_CHR_LEN) :: IncDir
   CALL GetEnv('MONDO_INC',IncDir)
   OPEN(UNIT=90,FILE=TRIM(IncDir)//'/Ell.m',STATUS='REPLACE')
   WRITE(*,*)TRIM(IncDir)//'Ell.m'
   WRITE(90,*)'BFEll=',BFEll,';'             
   WRITE(90,*)'HGEll=',HGEll,';'             
   WRITE(90,*)'SPEll=',SPEll,';'             
   CLOSE(UNIT=90,STATUS='KEEP')
END PROGRAM EllToInclude

