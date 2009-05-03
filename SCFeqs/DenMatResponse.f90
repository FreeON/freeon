!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
MODULE DenMatResponse
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE LinAlg
  USE DenMatMethods, ONLY: SpectralBounds,SetVarThresh
  USE MondoLogger

  IMPLICIT NONE
  PRIVATE

!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
  PUBLIC  :: DoTC2R
  PUBLIC  :: SetPPrmOld
  PUBLIC  :: SavePPrm
  PUBLIC  :: LoadMatrices
  PUBLIC  :: AllocArray
  PUBLIC  :: DeAllocArray
  PUBLIC  :: FockPrimGuess
  !
!---------------------------------------------------------------------------------
! PRIVATE DECLARATIONS
!---------------------------------------------------------------------------------
  PRIVATE :: PutXFormPrim
  PRIVATE :: TC2R_DMP
  PRIVATE :: CnvrgChckPrim
  PRIVATE :: Order0
  PRIVATE :: Order1_A
  PRIVATE :: Order2_AB,Order2_AA
  PRIVATE :: Order3_ABC,Order3_AAA,Order3_AAC,Order3_ABB
  PRIVATE :: Order4,Order5,Order6,Order7,Order8,Order9,Order10
  !
#ifdef TC2R_EIGENVAL
  INTERFACE DSYEV
     SUBROUTINE DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
       USE GlobalScalars
       CHARACTER(LEN=1), INTENT(IN)    :: JOBZ, UPLO
       INTEGER,          INTENT(IN)    :: LDA,  LWORK, N
       INTEGER,          INTENT(OUT)   :: INFO
       REAL(DOUBLE),     INTENT(INOUT) :: A(LDA,*)
       REAL(DOUBLE),     INTENT(OUT)   :: W(*)
       REAL(DOUBLE),     INTENT(OUT)   :: WORK(*)
     END SUBROUTINE DSYEV
  END INTERFACE
#endif
  !
CONTAINS
  !
  !
  SUBROUTINE DoTC2R(P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1,PPrmOld, &
                    T,Tmp1,Tmp2,Tmp3,RespOrder,SwitchThresh,Args)
    IMPLICIT NONE
    !
    !-------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(DBCSR), INTENT(INOUT) :: P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
    TYPE(DBCSR), INTENT(INOUT) :: PPrmOld,T,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR), INTENT(INOUT) :: P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
    TYPE( BCSR), INTENT(INOUT) :: PPrmOld,T,Tmp1,Tmp2,Tmp3
#endif
    INTEGER    , INTENT(IN   ) :: RespOrder
    REAL(DOUBLE),INTENT(IN   ) :: SwitchThresh
    TYPE(ARGMT)                :: Args
    !-------------------------------------------------------------------
    INTEGER                    :: MM,I
    REAL(DOUBLE)               :: Ne
    CHARACTER(LEN=*), PARAMETER:: Prog='TC2R'
    !-------------------------------------------------------------------
#ifdef TC2R_EIGENVAL
    INTEGER :: LWORK,Info
    TYPE(DBL_VECT)                 :: EigenV,Work
    TYPE(DBL_RNK2)                 :: PPrimA
#endif
    !-------------------------------------------------------------------
    !
    MM=0
    Ne=Half*DBLE(NEl)
    !
    DO I=1,100
       CALL TC2R_DMP(P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                     Tmp1,Tmp2,Tmp3,Ne,MM,RespOrder,SwitchThresh,Args)
       IF(CnvrgChckPrim(Prog,I,Ne,MM,T,P,PPrm1_1,PPrm2_1,PPrm3_1,PPrmOld,Tmp1,Tmp2,RespOrder)) EXIT
    ENDDO
    !
#ifdef TC2R_EIGENVAL
  CALL New(PPrimA,(/NBasF,NBasF/))
  CALL SetEq(PPrimA,PPrm1_1)
  CALL New(EigenV,NBasF)
  CALL SetEq(EigenV,Zero)
  LWORK=MAX(1,3*NBasF+10)
  CALL New(Work,LWork)
  CALL DSYEV('V','U',NBasF,PPrimA%D(1,1),NBasF,EigenV%D(1),Work%D(1),LWORK,Info)
  IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in TC2R. INFO='//TRIM(IntToChar(Info)))
  CALL PLOT_ERROR_WITH_BOXES(EigenV,41,'test'//TRIM(IntToChar(Args%I%I(1))))
  CALL Delete(PPrimA)
  CAll Delete(EigenV)
  CALL Delete(Work)
#endif
  !
  END SUBROUTINE DoTC2R
  !
  !
  SUBROUTINE SetPPrmOld(PPrmOld,PPrm1_1,PPrm2_1,PPrm3_1,RespOrder)
    IMPLICIT NONE
    !
    !-------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(DBCSR), INTENT(INOUT) :: PPrmOld,PPrm1_1,PPrm2_1,PPrm3_1
#else
    TYPE(BCSR ), INTENT(INOUT) :: PPrmOld,PPrm1_1,PPrm2_1,PPrm3_1
#endif
    INTEGER    , INTENT(IN   ) :: RespOrder
    !-------------------------------------------------------------------
    !
    SELECT CASE(RespOrder)
    CASE(1); CALL SetEq(PPrmOld,PPrm1_1)
    CASE(2); CALL SetEq(PPrmOld,PPrm2_1)
    CASE(3); CALL SetEq(PPrmOld,PPrm3_1)
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
  END SUBROUTINE SetPPrmOld
  !
  !
  SUBROUTINE SavePPrm(PPrm1_1,PPrm2_1,PPrm3_1,Tmp1,Tmp2,RespOrder,Args)
    IMPLICIT NONE
    !
    !-------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(DBCSR), INTENT(INOUT) :: PPrm1_1,PPrm2_1,PPrm3_1,Tmp1,Tmp2
#else
    TYPE(BCSR ), INTENT(INOUT) :: PPrm1_1,PPrm2_1,PPrm3_1,Tmp1,Tmp2
#endif
    INTEGER    , INTENT(IN   ) :: RespOrder
    TYPE(ARGMT)                :: Args
    !-------------------------------------------------------------------
    CHARACTER(LEN=*), PARAMETER:: Prog='TC2R'
    !-------------------------------------------------------------------
    !
    SELECT CASE(RespOrder)
    CASE(1); CALL PutXFormPrim(Prog,Args,PPrm1_1,Tmp1,Tmp2)
    CASE(2); CALL PutXFormPrim(Prog,Args,PPrm2_1,Tmp1,Tmp2)
    CASE(3); CALL PutXFormPrim(Prog,Args,PPrm3_1,Tmp1,Tmp2)
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
  END SUBROUTINE SavePPrm
  !
  !
  SUBROUTINE LoadMatrices(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
                          P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                          T,RespOrder,Args)
    IMPLICIT NONE
    !
    !-------------------------------------------------------------------
    TYPE(BCSR ), INTENT(INOUT) :: F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1
#ifdef PARALLEL
    TYPE(DBCSR), INTENT(INOUT) :: P,T
    TYPE(DBCSR), INTENT(INOUT) :: PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
#else
    TYPE(BCSR ), INTENT(INOUT) :: P,T
    TYPE(BCSR ), INTENT(INOUT) :: PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
#endif
    INTEGER    , INTENT(IN   ) :: RespOrder
    TYPE(ARGMT)                :: Args
    !-------------------------------------------------------------------
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: FFile
    CHARACTER(LEN=1)               :: Chr1,Chr2,Chr3
    INTEGER                        :: LastSCFCycle,LastCPSCFCycle
    LOGICAL                        :: EXIST,IsPresent
    !-------------------------------------------------------------------

    type(BCSR) :: S,Z
    TYPE(DBL_RNK2)                       :: BTmp
    integer :: iii,jjj
    !
    Chr1=TRIM(Args%C%C(3)(1:1))
    SELECT CASE(RespOrder)
    CASE(1);Chr2=' '                   ;Chr3=' '
    CASE(2);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=' '
    CASE(3);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=TRIM(Args%C%C(3)(3:3))
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
    ! Load Fock Matrix.
    CALL Get(LastSCFCycle,'lastscfcycle')
    CALL Get(F,TrixFile('OrthoF',Args,LastSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
    !
    ! Load FockPrim Matrices.
    SELECT CASE(RespOrder)
    CASE(1)
       ! Load FockPrime Matrix.
       FFile=TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(3)),Args,0)
       INQUIRE(FILE=FFile,EXIST=IsPresent)
       IF(IsPresent) THEN
          CALL Get(FPrm1_1,TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(3)),Args,0),BCast_O=.FALSE.)
       ELSE
          CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)),Args,0),BCast_O=.FALSE.)
       ENDIF
    CASE(2)
       ! Load FockPrime Matrix.
       IF(Chr1.EQ.Chr2) THEN
          ! PPrm2_1 <-> aa
          ! PPrm1_1 <-> a
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
          CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
       ELSE
          ! PPrm2_1 <-> ab
          ! PPrm1_1 <-> a
          ! PPrm1_2 <-> b
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
          CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:2)))
          CALL Get(FPrm1_2,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:2)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
       ENDIF
       !
       ! Load FockPrim2 Matrix.
       FFile=TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(3)),Args,0)
       INQUIRE(FILE=FFile,EXIST=IsPresent)
       IF(IsPresent) THEN
          CALL Get(FPrm2_1,TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(3)),Args,0),BCast_O=.FALSE.)
       ELSE
          CALL Get(FPrm2_1,TrixFile('OrthoFPrime' //TRIM(Args%C%C(3)),Args,0),BCast_O=.FALSE.)
       ENDIF
       !
    CASE(3)
       ! Load FockPrime Matrix.
       IF(Chr1.EQ.Chr2.AND.Chr1.EQ.Chr3) THEN
          ! PPrm3_1 <-> aaa
          ! PPrm2_1 <-> aa
          ! PPrm1_1 <-> a
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
          CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:2)))
          CALL Get(FPrm2_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:2)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
       ELSEIF(Chr1.EQ.Chr2.AND.Chr1.NE.Chr3) THEN
          ! PPrm3_1 <-> aac
          ! PPrm2_1 <-> aa
          ! PPrm2_2 <-> ac
          ! PPrm1_1 <-> a
          ! PPrm1_3 <-> c
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
          CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(3:3)))
          CALL Get(FPrm1_3,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(3:3)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:2)))
          CALL Get(FPrm2_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:2)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:3)))
          CALL Get(FPrm2_2,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:3)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
       ELSEIF(Chr1.NE.Chr2.AND.Chr2.EQ.Chr3) THEN
          ! PPrm3_1 <-> abb
          ! PPrm2_1 <-> ab
          ! PPrm2_3 <-> bb
          ! PPrm1_1 <-> a
          ! PPrm1_2 <-> b
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
          CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:2)))
          CALL Get(FPrm1_2,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:2)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:2)))
          CALL Get(FPrm2_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:2)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:3)))
          CALL Get(FPrm2_3,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:3)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
       ELSEIF(Chr1.NE.Chr2.AND.Chr1.NE.Chr3.AND.Chr2.NE.Chr3) THEN
          ! PPrm3_1 <-> abc
          ! PPrm2_1 <-> ab
          ! PPrm2_2 <-> ac
          ! PPrm2_3 <-> bc
          ! PPrm1_1 <-> a
          ! PPrm1_2 <-> b
          ! PPrm1_3 <-> c
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
          CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:2)))
          CALL Get(FPrm1_2,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:2)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(3:3)))
          CALL Get(FPrm1_3,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(3:3)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:2)))
          CALL Get(FPrm2_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:2)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1))//TRIM(Args%C%C(3)(3:3)))
          CALL Get(FPrm2_2,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1))//TRIM(Args%C%C(3)(3:3)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
          !
          CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:3)))
          CALL Get(FPrm2_3,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:3)), &
                   Args,LastCPSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
       ELSE
          CALL Halt('Response: unknown symmetry <'//Chr1//Chr2//Chr3//'>.')
       ENDIF
       !
       ! Load FockPrim3 Matrix.
       FFile=TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(3)),Args,0)
       INQUIRE(FILE=FFile,EXIST=IsPresent)
       IF(IsPresent) THEN
          CALL Get(FPrm3_1,TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(3)),Args,0),BCast_O=.FALSE.)
       ELSE
          CALL Get(FPrm3_1,TrixFile('OrthoFPrime' //TRIM(Args%C%C(3)),Args,0),BCast_O=.FALSE.)
       ENDIF
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
    ! Load perturbation (dipole moment at this moment).
    FFile='Ortho'//TRIM(Args%C%C(4))//TRIM(Args%C%C(3)(1:1))
    CALL Get(T,TrixFile(FFile,Args))

!!$    CALL SetEq(BTmp,F)
!!$    open(20)
!!$    do iii=1,size(BTmp%D,1)
!!$       do jjj=1,size(BTmp%D,1)
!!$          write(20,'(E26.15)') BTmp%D(iii,jjj)
!!$       enddo
!!$    enddo
!!$    close(20)
!!$
!!$    select case(TRIM(Args%C%C(3)))
!!$    case('X')
!!$       CALL SetEq(BTmp,T)
!!$       open(21)
!!$       do iii=1,size(BTmp%D,1)
!!$          do jjj=1,size(BTmp%D,1)
!!$             write(21,'(E26.15)') BTmp%D(iii,jjj)
!!$          enddo
!!$       enddo
!!$       close(21)
!!$    case('Y')
!!$       CALL SetEq(BTmp,T)
!!$       open(22)
!!$       do iii=1,size(BTmp%D,1)
!!$          do jjj=1,size(BTmp%D,1)
!!$             write(22,'(E26.15)') BTmp%D(iii,jjj)
!!$          enddo
!!$       enddo
!!$       close(22)
!!$    case('Z')
!!$       CALL SetEq(BTmp,T)
!!$       open(23)
!!$       do iii=1,size(BTmp%D,1)
!!$          do jjj=1,size(BTmp%D,1)
!!$             write(23,'(E26.15)') BTmp%D(iii,jjj)
!!$          enddo
!!$       enddo
!!$       close(23)
!!$    end select
!!$    CALL Get(S,TrixFile('S',Args))
!!$    CALL SetEq(BTmp,S)
!!$    open(20)
!!$    do iii=1,size(BTmp%D,1)
!!$       do jjj=1,size(BTmp%D,1)
!!$          write(20,'(E26.15)') BTmp%D(iii,jjj)
!!$       enddo
!!$    enddo
!!$    close(20)
!!$    CALL Get(Z,TrixFile('Z',Args))
!!$    CALL SetEq(BTmp,Z)
!!$    open(23)
!!$    do iii=1,size(BTmp%D,1)
!!$       do jjj=1,size(BTmp%D,1)
!!$          write(23,'(E26.15)') BTmp%D(iii,jjj)
!!$       enddo
!!$    enddo
!!$    close(23)
    !
  END SUBROUTINE LoadMatrices
  !
  !
  SUBROUTINE AllocArray(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
                        P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                        PPrmOld,T,Tmp1,Tmp2,Tmp3,RespOrder,Args)
    IMPLICIT NONE
    !
    !-------------------------------------------------------------------
    TYPE(BCSR ), INTENT(INOUT) :: F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1
#ifdef PARALLEL
    TYPE(DBCSR), INTENT(INOUT) :: P,PPrmOld,T,Tmp1,Tmp2,Tmp3
    TYPE(DBCSR), INTENT(INOUT) :: PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
#else
    TYPE(BCSR ), INTENT(INOUT) :: P,PPrmOld,T,Tmp1,Tmp2,Tmp3
    TYPE(BCSR ), INTENT(INOUT) :: PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
#endif
    INTEGER    , INTENT(IN   ) :: RespOrder
    TYPE(ARGMT)                :: Args
    !-------------------------------------------------------------------
    CHARACTER(LEN=1)           :: Chr1,Chr2,Chr3
    !-------------------------------------------------------------------
    !
    Chr1=TRIM(Args%C%C(3)(1:1))
    SELECT CASE(RespOrder)
    CASE(1);Chr2=' '                   ;Chr3=' '
    CASE(2);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=' '
    CASE(3);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=TRIM(Args%C%C(3)(3:3))
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
    CALL New(F)
    CALL New(T)
    CALL New(P)
    CALL New(Tmp1)
    CALL New(Tmp2)
    CALL New(Tmp3)
    CALL New(PPrmOld)
    !
    SELECT CASE(RespOrder)
    CASE(1)
       CALL New(FPrm1_1)
       CALL New(PPrm1_1)
    CASE(2)
       CALL New(FPrm1_1)
       CALL New(PPrm1_1)
       CALL New(FPrm2_1)
       CALL New(PPrm2_1)
       IF(Chr1.EQ.Chr2) THEN
       ELSE
          CALL New(FPrm1_2)
          CALL New(PPrm1_2)
       ENDIF
    CASE(3)
       CALL New(FPrm1_1)
       CALL New(PPrm1_1)
       CALL New(FPrm2_1)
       CALL New(PPrm2_1)
       CALL New(FPrm3_1)
       CALL New(PPrm3_1)
       IF(Chr1.EQ.Chr2.AND.Chr1.EQ.Chr3) THEN
       ELSEIF(Chr1.EQ.Chr2.AND.Chr1.NE.Chr3) THEN
          CALL New(FPrm1_3)
          CALL New(PPrm1_3)
          CALL New(FPrm2_2)
          CALL New(PPrm2_2)
       ELSEIF(Chr1.NE.Chr2.AND.Chr2.EQ.Chr3) THEN
          CALL New(FPrm1_2)
          CALL New(PPrm1_2)
          CALL New(FPrm2_3)
          CALL New(PPrm2_3)
       ELSEIF(Chr1.NE.Chr2.AND.Chr1.NE.Chr3.AND.Chr2.NE.Chr3) THEN
          CALL New(FPrm1_2)
          CALL New(PPrm1_2)
          CALL New(FPrm1_3)
          CALL New(PPrm1_3)
          CALL New(FPrm2_2)
          CALL New(PPrm2_2)
          CALL New(FPrm2_3)
          CALL New(PPrm2_3)
       ELSE
          CALL Halt('Response: unknown symmetry <'//Chr1//Chr2//Chr3//'>.')
       ENDIF
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT

  END SUBROUTINE AllocArray
  !
  !
  SUBROUTINE DeAllocArray(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
                          P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                          PPrmOld,T,Tmp1,Tmp2,Tmp3,RespOrder,Args)
    IMPLICIT NONE
    !
    !-------------------------------------------------------------------
    TYPE(BCSR ), INTENT(INOUT) :: F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1
#ifdef PARALLEL
    TYPE(DBCSR), INTENT(INOUT) :: P,PPrmOld,T,Tmp1,Tmp2,Tmp3
    TYPE(DBCSR), INTENT(INOUT) :: PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
#else
    TYPE(BCSR ), INTENT(INOUT) :: P,PPrmOld,T,Tmp1,Tmp2,Tmp3
    TYPE(BCSR ), INTENT(INOUT) :: PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
#endif
    INTEGER    , INTENT(IN   ) :: RespOrder
    TYPE(ARGMT)                :: Args
    !-------------------------------------------------------------------
    CHARACTER(LEN=1)           :: Chr1,Chr2,Chr3
    !-------------------------------------------------------------------
    !
    Chr1=TRIM(Args%C%C(3)(1:1))
    SELECT CASE(RespOrder)
    CASE(1);Chr2=' '                   ;Chr3=' '
    CASE(2);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=' '
    CASE(3);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=TRIM(Args%C%C(3)(3:3))
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
    CALL Delete(P)
    CALL Delete(T)
    CALL Delete(F)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)
    CALL Delete(Tmp3)
    CALL Delete(PPrmOld)
    !
    SELECT CASE(RespOrder)
    CASE(1)
       CALL Delete(PPrm1_1)
       CALL Delete(FPrm1_1)
    CASE(2)
       CALL Delete(PPrm1_1)
       CALL Delete(PPrm2_1)
       CALL Delete(FPrm1_1)
       CALL Delete(FPrm2_1)
       IF(Chr1.NE.Chr2) THEN
          CALL Delete(PPrm1_2)
          CALL Delete(FPrm1_2)
       ENDIF
    CASE(3)
       CALL Delete(PPrm1_1)
       CALL Delete(PPrm2_1)
       CALL Delete(PPrm3_1)
       CALL Delete(FPrm1_1)
       CALL Delete(FPrm2_1)
       CALL Delete(FPrm3_1)
       IF(Chr1.EQ.Chr2.AND.Chr1.EQ.Chr3) THEN
       ELSEIF(Chr1.EQ.Chr2.AND.Chr1.NE.Chr3) THEN
          CALL Delete(PPrm1_3)
          CALL Delete(PPrm2_2)
          CALL Delete(FPrm1_3)
          CALL Delete(FPrm2_2)
       ELSEIF(Chr1.NE.Chr2.AND.Chr2.EQ.Chr3) THEN
          CALL Delete(PPrm1_2)
          CALL Delete(PPrm2_3)
          CALL Delete(FPrm1_2)
          CALL Delete(FPrm2_3)
       ELSEIF(Chr1.NE.Chr2.AND.Chr1.NE.Chr3.AND.Chr2.NE.Chr3) THEN
          CALL Delete(PPrm1_2)
          CALL Delete(PPrm1_3)
          CALL Delete(PPrm2_2)
          CALL Delete(PPrm2_3)
          CALL Delete(FPrm1_2)
          CALL Delete(FPrm1_3)
          CALL Delete(FPrm2_2)
          CALL Delete(FPrm2_3)
       ELSE
          CALL Halt('Response: unknown symmetry <'//Chr1//Chr2//Chr3//'>.')
       ENDIF
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
  END SUBROUTINE DeAllocArray
  !
  !
  SUBROUTINE FockPrimGuess(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
                           P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                           RespOrder,Args)
!H---------------------------------------------------------------------------------
!H SUBROUTINE FockPrimGuess
!H  Set the FockPrim guess for TC2R.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BCSR ), INTENT(INOUT) :: F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1
#ifdef PARALLEL
    TYPE(DBCSR), INTENT(INOUT) :: P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
#else
    TYPE(BCSR ), INTENT(INOUT) :: P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
#endif
    INTEGER    , INTENT(IN   ) :: RespOrder
    TYPE(ARGMT)                :: Args
    !-------------------------------------------------------------------
    REAL(DOUBLE)               :: Fmin,Fmax,DF,Coeff
    CHARACTER(LEN=1)           :: Chr1,Chr2,Chr3
    !-------------------------------------------------------------------
    !
    Chr1=TRIM(Args%C%C(3)(1:1))
    SELECT CASE(RespOrder)
    CASE(1);Chr2=' '                   ;Chr3=' '
    CASE(2);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=' '
    CASE(3);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=TRIM(Args%C%C(3)(3:3))
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
    ! Estimate spectral bounds.
#ifdef PARALLEL
    IF(MyId==ROOT) CALL SpectralBounds(F,Fmin,Fmax)
    CALL BCast(FMin)
    CALL BCast(FMax)
#else
    CALL SpectralBounds(F,Fmin,Fmax)
#endif
    !
    DF=(Fmax-Fmin)
    Coeff=-One/DF
    !
    ! Set up P.
    Call SetEq(P,F)
    CALL Add(P,-Fmax)
    CALL Multiply(P,Coeff)          ! P = (I*F_max-F)/DF = X0
    !
    ! Set up PPrim.
    SELECT CASE(RespOrder)
    CASE(1)
       CALL SetEq(PPrm1_1,FPrm1_1)
       CALL Multiply(PPrm1_1,Coeff)
    CASE(2)
       CALL SetEq(PPrm1_1,FPrm1_1)
       CALL Multiply(PPrm1_1,Coeff)
       !
       IF(Chr1.EQ.Chr2) THEN
          ! PPrm2_1 <-> aa
          ! PPrm1_1 <-> a
       ELSE
          ! PPrm2_1 <-> ab
          ! PPrm1_1 <-> a
          ! PPrm1_2 <-> b
          CALL SetEq(PPrm1_1,FPrm1_1)
          CALL Multiply(PPrm1_1,Coeff)
          !
          CALL SetEq(PPrm1_2,FPrm1_2)
          CALL Multiply(PPrm1_2,Coeff)
       ENDIF
       !
       CALL SetEq(PPrm2_1,FPrm2_1)
       CALL Multiply(PPrm2_1,Coeff)
    CASE(3)
       CALL SetEq(PPrm1_1,FPrm1_1)
       CALL Multiply(PPrm1_1,Coeff)
       !
       CALL SetEq(PPrm2_1,FPrm2_1)
       CALL Multiply(PPrm2_1,Coeff)
       !
       IF(Chr1.EQ.Chr2.AND.Chr1.EQ.Chr3) THEN
          ! PPrm3_1 <-> aaa
          ! PPrm2_1 <-> aa
          ! PPrm1_1 <-> a
       ELSEIF(Chr1.EQ.Chr2.AND.Chr1.NE.Chr3) THEN
          ! PPrm3_1 <-> aac
          ! PPrm2_1 <-> aa
          ! PPrm2_2 <-> ac
          ! PPrm1_1 <-> a
          ! PPrm1_3 <-> c
          CALL SetEq(PPrm1_3,FPrm1_3)
          CALL Multiply(PPrm1_3,Coeff)
          !
          CALL SetEq(PPrm2_2,FPrm2_2)
          CALL Multiply(PPrm2_2,Coeff)
       ELSEIF(Chr1.NE.Chr2.AND.Chr2.EQ.Chr3) THEN
          ! PPrm3_1 <-> abb
          ! PPrm2_1 <-> ab
          ! PPrm2_3 <-> bb
          ! PPrm1_1 <-> a
          ! PPrm1_2 <-> b
          CALL SetEq(PPrm1_2,FPrm1_2)
          CALL Multiply(PPrm1_2,Coeff)
          !
          CALL SetEq(PPrm2_3,FPrm2_3)
          CALL Multiply(PPrm2_3,Coeff)
       ELSEIF(Chr1.NE.Chr2.AND.Chr1.NE.Chr3.AND.Chr2.NE.Chr3) THEN
          ! PPrm3_1 <-> abc
          ! PPrm2_1 <-> ab
          ! PPrm2_2 <-> ac
          ! PPrm2_3 <-> bc
          ! PPrm1_1 <-> a
          ! PPrm1_2 <-> b
          ! PPrm1_3 <-> c
          CALL SetEq(PPrm1_2,FPrm1_2)
          CALL Multiply(PPrm1_2,Coeff)
          !
          CALL SetEq(PPrm1_3,FPrm1_3)
          CALL Multiply(PPrm1_3,Coeff)
          !
          CALL SetEq(PPrm2_2,FPrm2_2)
          CALL Multiply(PPrm2_2,Coeff)
          !
          CALL SetEq(PPrm2_3,FPrm2_3)
          CALL Multiply(PPrm2_3,Coeff)
       ELSE
          CALL Halt('Response: unknown symmetry <'//Chr1//Chr2//Chr3//'>.')
       ENDIF
       CALL SetEq(PPrm3_1,FPrm3_1)
       CALL Multiply(PPrm3_1,Coeff)
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
  END SUBROUTINE FockPrimGuess
  !
  !
  SUBROUTINE PutXFormPrim(Prog,Args,PPrm,Z,Tmp)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PutXFormPrim(Prog,Args,PPrm,Z,Tmp)
!H  Save onto the disc the OrthoPPrim and the PPrim matrices.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(DBCSR)     , INTENT(INOUT) :: PPrm,Z,Tmp
#else
    TYPE(BCSR )     , INTENT(INOUT) :: PPrm,Z,Tmp
#endif
    TYPE(ARGMT)     , INTENT(IN   ) :: Args
    CHARACTER(LEN=*), INTENT(IN   ) :: Prog
    !-------------------------------------------------------------------
    LOGICAL                         :: IsPresent
    !-------------------------------------------------------------------
    !
    ! IO for the orthogonal PPrm
    CALL Put(PPrm,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)),Args,1))
    CALL PChkSum(PPrm,'OrthoDPrime'//TRIM(Args%C%C(3))//'['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint( PPrm,'OrthoDPrime'//TRIM(Args%C%C(3))//'['//TRIM(NxtCycl)//']')
    CALL Plot(   PPrm,'OrthoDPrime'//TRIM(Args%C%C(3))//'_'//TRIM(NxtCycl))
    !
    ! Convert to AO representation
    INQUIRE(FILE=TrixFile('X',Args),EXIST=IsPresent)
    IF(IsPresent)THEN
       CALL Get(Z,TrixFile('X',Args))   ! Z=S^(-1/2)
       CALL Multiply(Z,PPrm,Tmp)
       CALL Multiply(Tmp,Z,PPrm)
    ELSE
       CALL Get(Z,TrixFile('Z',Args))   ! Z=S^(-L)
       CALL Multiply(Z,PPrm,Tmp)
       CALL Get(Z,TrixFile('ZT',Args))
       CALL Multiply(Tmp,Z,PPrm)
    ENDIF
    CALL Filter(Tmp,PPrm)              ! Thresholding
    !
    ! IO for the non-orthogonal PPrm
    CALL Put(Tmp,TrixFile('DPrime'//TRIM(Args%C%C(3)),Args,1))
    CALL PChkSum(Tmp,'DPrime'//TRIM(Args%C%C(3))//'['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint( Tmp,'DPrime'//TRIM(Args%C%C(3))//'['//TRIM(NxtCycl)//']')
    CALL Plot(   Tmp,'DPrime'//TRIM(Args%C%C(3))//'_'//TRIM(NxtCycl))
    !
  END SUBROUTINE PutXFormPrim
  !
  !
  SUBROUTINE TC2R_DMP(P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                      Tmp1,Tmp2,Tmp3,Ne,MM,RespOrder,SwitchThresh,Args)
!H---------------------------------------------------------------------------------
!H SUBROUTINE TC2R_DMP(P,PPrim,Tmp1,Tmp2,Tmp3,Ne,MM)
!H  This routine does the linear, quadratic and cubic TC2Response.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(DBCSR) , INTENT(INOUT) :: P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                                   Tmp1,Tmp2,Tmp3
#else
    TYPE(BCSR ) , INTENT(INOUT) :: P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                                   Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE), INTENT(IN   ) :: Ne,SwitchThresh
    INTEGER     , INTENT(IN   ) :: RespOrder
    INTEGER     , INTENT(INOUT) :: MM
    TYPE(ARGMT)                 :: Args
    !-------------------------------------------------------------------
    REAL(DOUBLE)                :: N
    CHARACTER(LEN=1)            :: Chr1,Chr2,Chr3
    LOGICAL, SAVE               :: Switching=.FALSE.
    LOGICAL, SAVE               :: Switch=.FALSE.
    !-------------------------------------------------------------------
    !
    Chr1=TRIM(Args%C%C(3)(1:1))
    SELECT CASE(RespOrder)
    CASE(1);Chr2=' '                   ;Chr3=' '
    CASE(2);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=' '
    CASE(3);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=TRIM(Args%C%C(3)(3:3))
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
    MM=MM+1
    N=Trace(P)
    !
    IF(ABS(N-Ne).LT.SwitchThresh.OR.Switching) THEN
       Switching=.TRUE.
       IF(Switch) THEN
          N=Ne+1.00D-3
          Switch=.FALSE.
       ELSE
          N=Ne-1.00D-3
          Switch=.TRUE.
       ENDIF
    ENDIF
    !
    SELECT CASE(RespOrder)
    CASE(1)
       ! Da = {Da,X0}
       ! PPrm1_1 <-> a
       CALL Order1_A(P,PPrm1_1,Tmp1,Tmp2,Tmp3,N,Ne)
    CASE(2)
       IF(Chr1.EQ.Chr2) THEN
          ! Daa = {Daa,X0} + Da*Da
          ! PPrm2_1 <-> aa
          ! PPrm1_1 <-> a
          CALL Order2_AA(P,PPrm1_1,PPrm2_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_1,Tmp1,Tmp2,Tmp3,N,Ne)
       ELSE
          ! Dab = {Dab,X0} + 1/2{Da,Db}
          ! PPrm2_1 <-> ab
          ! PPrm1_1 <-> a
          ! PPrm1_2 <-> b
          CALL Order2_AB(P,PPrm1_1,PPrm1_2,PPrm2_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_2,Tmp1,Tmp2,Tmp3,N,Ne)
       ENDIF
    CASE(3)
       IF(Chr1.EQ.Chr2.AND.Chr1.EQ.Chr3) THEN
          ! Daaa = {Dabc,X0} + {Daa,Da}
          ! PPrm3_1 <-> aaa
          ! PPrm2_1 <-> aa
          ! PPrm1_1 <-> a
          CALL Order3_AAA(P,PPrm1_1,PPrm2_1,PPrm3_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order2_AA(P,PPrm1_1,PPrm2_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_1,Tmp1,Tmp2,Tmp3,N,Ne)
       ELSEIF(Chr1.EQ.Chr2.AND.Chr1.NE.Chr3) THEN
          ! Daac = {Daac,X0} + 1/3{Daa,Dc} + 2/3{Dac,Da}
          ! PPrm3_1 <-> aac
          ! PPrm2_1 <-> aa
          ! PPrm2_2 <-> ac
          ! PPrm1_1 <-> a
          ! PPrm1_3 <-> c
          CALL Order3_AAC(P,PPrm1_1,PPrm1_3,PPrm2_1,PPrm2_2,PPrm3_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order2_AA(P,PPrm1_1,PPrm2_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order2_AB(P,PPrm1_1,PPrm1_3,PPrm2_2,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_3,Tmp1,Tmp2,Tmp3,N,Ne)
       ELSEIF(Chr1.NE.Chr2.AND.Chr2.EQ.Chr3) THEN
          ! Dabb = {Dabb,X0} + 2/3{Dab,Db} + 1/3{Dbb,Da}
          ! PPrm3_1 <-> abb
          ! PPrm2_1 <-> ab
          ! PPrm2_3 <-> bb
          ! PPrm1_1 <-> a
          ! PPrm1_2 <-> b
          CALL Order3_ABB(P,PPrm1_1,PPrm1_2,PPrm2_1,PPrm2_3,PPrm3_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order2_AB(P,PPrm1_1,PPrm1_2,PPrm2_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order2_AA(P,PPrm1_2,PPrm2_3,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_2,Tmp1,Tmp2,Tmp3,N,Ne)
       ELSEIF(Chr1.NE.Chr2.AND.Chr1.NE.Chr3.AND.Chr2.NE.Chr3) THEN
          ! Dabc = {Dabc,X0} + 1/3{Dab,Dc} + 1/3{Dac,Db} + 1/3{Dbc,Da}
          ! PPrm3_1 <-> abc
          ! PPrm2_1 <-> ab
          ! PPrm2_2 <-> ac
          ! PPrm2_3 <-> bc
          ! PPrm1_1 <-> a
          ! PPrm1_2 <-> b
          ! PPrm1_3 <-> c
          CALL Order3_ABC(P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order2_AB(P,PPrm1_1,PPrm1_2,PPrm2_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order2_AB(P,PPrm1_1,PPrm1_3,PPrm2_2,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order2_AB(P,PPrm1_2,PPrm1_3,PPrm2_3,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_1,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_2,Tmp1,Tmp2,Tmp3,N,Ne)
          CALL Order1_A(P,PPrm1_3,Tmp1,Tmp2,Tmp3,N,Ne)
       ELSE
          CALL Halt('Response: unknown symmetry <'//Chr1//Chr2//Chr3//'>.')
       ENDIF
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
    CALL Order0(P,Tmp1,Tmp2,Tmp3,N,Ne)
    !
  END SUBROUTINE TC2R_DMP
  !
  !
  LOGICAL FUNCTION CnvrgChckPrim(Prog,NPur,Ne,MM,T,P,PPrim,PPrm2,PPrm3,PPrmOld,Tmp1,Tmp2,RespOrder)
!H---------------------------------------------------------------------------------
!H FUNCTION CnvrgChckPrim
!H  Check the convergence during TC2R steps.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(DBCSR)                     :: T,P,PPrim,PPrm2,PPrm3,PPrmOld,Tmp1,Tmp2
#else
    TYPE(BCSR )                     :: T,P,PPrim,PPrm2,PPrm3,PPrmOld,Tmp1,Tmp2
#endif
    REAL(DOUBLE)    , INTENT(IN   ) :: Ne
    CHARACTER(LEN=*), INTENT(IN   ) :: Prog
    INTEGER         , INTENT(IN   ) :: NPur,RespOrder
    INTEGER         , INTENT(INOUT) :: MM
    !-------------------------------------------------------------------
    REAL(DOUBLE)                    :: Prop,AbsErrProp,RelErrProp
    REAL(DOUBLE)                    :: AveErrProp,PNon0,PPrmNon0
    REAL(DOUBLE)                    :: AbsErrPPrim,FNormErrPPrim,TRacePPrim
    REAL(DOUBLE)    ,SAVE           :: OldProp,OldAEP
    CHARACTER(LEN=2*DEFAULT_CHR_LEN):: Mssg,MssgTag,CnvrgCmmnt
    INTEGER, PARAMETER              :: CyclMin=10
    !-------------------------------------------------------------------
    !
    IF(NPur==0)THEN
       OldProp=BIG_DBL
       OldAEP=BIG_DBL
    ENDIF
#ifdef PARALLEL
    PNon0=DBLE(Reduce(P%NNon0))
    PNon0=100.D0*PNon0/DBLE(NBasF*NBasF)
    SELECT CASE(RespOrder)
    CASE(1)
       PPrmNon0=DBLE(Reduce(PPrim%NNon0))
       PPrmNon0=100.D0*DBLE(PPrmNon0)/DBLE(NBasF*NBasF)
    CASE(2)
       PPrmNon0=DBLE(Reduce(PPrm2%NNon0))
       PPrmNon0=100.D0*DBLE(PPrmNon0)/DBLE(NBasF*NBasF)
    CASE(3)
       PPrmNon0=DBLE(Reduce(PPrm3%NNon0))
       PPrmNon0=100.D0*DBLE(PPrmNon0)/DBLE(NBasF*NBasF)
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
#else
    PNon0=100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
    SELECT CASE(RespOrder)
    CASE(1); PPrmNon0=100.D0*DBLE(PPrim%NNon0)/DBLE(NBasF*NBasF)
    CASE(2); PPrmNon0=100.D0*DBLE(PPrm2%NNon0)/DBLE(NBasF*NBasF)
    CASE(3); PPrmNon0=100.D0*DBLE(PPrm3%NNon0)/DBLE(NBasF*NBasF)
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
#endif
    !
    ! Density matrix errors.
    CALL Multiply(PPrmOld,-One)
    !
    SELECT CASE(RespOrder)
    CASE(1); CALL Add(PPrmOld,PPrim,Tmp1)
    CASE(2); CALL Add(PPrmOld,PPrm2,Tmp1)
    CASE(3); CALL Add(PPrmOld,PPrm3,Tmp1)
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
    AbsErrPPrim=ABS(Max(Tmp1)+1.D-20)
    !
    ! Check for divergence steps.
    IF(AbsErrPPrim.GT.1D+10) &
           CALL Halt('The absolute derivative density matrix error is too big! AbsErrPPrim='&
                      //TRIM(DblToShrtChar(AbsErrPPrim)))
    !
    ! Compute Frobenus norm.
    FNormErrPPrim=FNorm(Tmp1)
    !
    ! Energy errors
#ifdef PARALLEL
    SELECT CASE(RespOrder)
    CASE(1)
       CALL Multiply(PPrim,T,Tmp1)
       Prop=-Two   *Trace(Tmp1)
    CASE(2)
       CALL Multiply(PPrm2,T,Tmp1)
       Prop=-Four  *Trace(Tmp1)
    CASE(3)
       CALL Multiply(PPrm3,T,Tmp1)
       Prop=-12.0d0*Trace(Tmp1)
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
#else
    SELECT CASE(RespOrder)
    CASE(1); Prop=-Two   *Trace(PPrim,T)
    CASE(2); Prop=-Four  *Trace(PPrm2,T)
    CASE(3); Prop=-12.0d0*Trace(PPrm3,T)
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !if(RespOrder==2)CALL PChkSum(PPrm2,'PPrm2_1'//TRIM(Args%C%C(3)),Prog)
    !if(RespOrder==2)write(*,*) 'Prop',Prop
#endif
    !
    ! Compute Abs and Rel error.
    AbsErrProp=ABS(OldProp-Prop)
    ! should check if it is zero
    IF(ABS(Prop).LT.1.00D-8) THEN
       RelErrProp=AbsErrProp
    ELSE
       RelErrProp=AbsErrProp/ABS(Prop)
    ENDIF
    !if(RespOrder==2)write(*,*) 'AbsErrProp',AbsErrProp,'AbsErrPPrim',AbsErrPPrim
    !if(RespOrder==2)write(*,*) 'Thresholds%ETol',Thresholds%ETol,'Thresholds%DTol',Thresholds%DTol
    !
    ! Set convergence check.
    CnvrgChckPrim=.FALSE.
    !
    IF(RelErrProp<SQRT(Thresholds%ETol/100d0).AND.AbsErrPPrim<Thresholds%DTol.AND.MM.GT.CyclMin) THEN
       CnvrgChckPrim=.TRUE.
       SELECT CASE(RespOrder)
       CASE(1); CnvrgCmmnt='Met dE1 goals'
       CASE(2); CnvrgCmmnt='Met dE2 goals'
       CASE(3); CnvrgCmmnt='Met dE3 goals'
       CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
       END SELECT
    !ELSEIF(AbsErrPPrim<Thresholds%DTol/10.0d0.AND.ABS(Prop).LT.1.00D-10) THEN
    !   CnvrgChckPrim=.TRUE.
    !   SELECT CASE(RespOrder)
    !   CASE(1); CnvrgCmmnt='Met dP1 goals'
    !   CASE(2); CnvrgCmmnt='Met dP2 goals'
    !   CASE(3); CnvrgCmmnt='Met dP3 goals'
    !   CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    !   END SELECT
    !ELSEIF(RelErrProp<SQRT(Thresholds%ETol).AND.ABS(Prop).GT.1.00D-10)THEN
    !   IF(AbsErrPPrim>OldAEP)THEN
    !      CnvrgChckPrim=.TRUE.
    !      SELECT CASE(RespOrder)
    !      CASE(1); CnvrgCmmnt='dP1 increase'
    !      CASE(2); CnvrgCmmnt='dP2 increase'
    !      CASE(3); CnvrgCmmnt='dP3 increase'
    !      CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    !      END SELECT
    !   ENDIF
    ENDIF
    !
    ! Updtate previous cycle values
    OldProp=Prop
    OldAEP=AbsErrPPrim
    !
    ! Print convergence stats
    MssgTag = 'Pure '//TRIM(IntToChar(NPur))
    SELECT CASE(RespOrder)
    CASE(1)
      Mssg= 'dE1='//TRIM(DblToShrtChar(RelErrProp))                &
            //', dP1='//TRIM(DblToShrtChar(AbsErrPPrim))           &
            //', %Non0='//TRIM(DblToShrtChar(PPrmNon0))
    CASE(2)
       Mssg='dE2='//TRIM(DblToShrtChar(RelErrProp))                &
            //', dP2='//TRIM(DblToShrtChar(AbsErrPPrim))           &
            //', %Non0='//TRIM(DblToShrtChar(PPrmNon0))
    CASE(3)
       Mssg='dE3='//TRIM(DblToShrtChar(RelErrProp))                &
            //', dP3='//TRIM(DblToShrtChar(AbsErrPPrim))           &
            //', %Non0='//TRIM(DblToShrtChar(PPrmNon0))
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT

    IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
       CALL MondoLog(DEBUG_MAXIMUM, Prog, Mssg, MssgTag)
#ifdef PARALLEL
    ENDIF
#endif
    ENDIF
    !
    ! Set thresholding for next cycle
    CALL SetVarThresh(MM)
    !
    ! Look for convergence
    IF(.NOT.CnvrgChckPrim)THEN
       SELECT CASE(RespOrder)
       CASE(1); CALL SetEq(PPrmOld,PPrim)
       CASE(2); CALL SetEq(PPrmOld,PPrm2)
       CASE(3); CALL SetEq(PPrmOld,PPrm3)
       CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
       END SELECT
       RETURN
    ENDIF
    !
    ! Increment the TC2R counter.
    MM=MM+1
    !
    ! Print summary stats.
    SELECT CASE(RespOrder)
    CASE(1); TracePPrim = Trace(PPrim)
    CASE(2); TracePPrim = Trace(PPrm2)
    CASE(3); TracePPrim = Trace(PPrm3)
    CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
    END SELECT
    !
    IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
       CALL MondoLog(DEBUG_NONE, Prog, CnvrgCmmnt//'Tr{PT}='//TRIM(DblToChar(Prop)))

       Mssg=TRIM(IntToChar(NPur))//' purification steps, '//TRIM(IntToChar(MM))//' matrix multiplies'
       CALL MondoLog(DEBUG_MAXIMUM, Prog, Mssg)

       Mssg='ThrX='//TRIM(DblToShrtChar(Thresholds%Trix))//', %Non0s = '//TRIM(DblToShrtChar(PPrmNon0))
       CALL MondoLog(DEBUG_MAXIMUM, Prog, Mssg)

       SELECT CASE(RespOrder)
       CASE(1)
          Mssg='Max abs errors'//'dE1='//TRIM(DblToShrtChar(RelErrProp))//', '        &
                               //'dP1='//TRIM(DblToShrtChar(AbsErrPPrim))
       CASE(2)
          Mssg='Max abs errors'//'dE2='//TRIM(DblToShrtChar(RelErrProp))//', '        &
                               //'dP2='//TRIM(DblToShrtChar(AbsErrPPrim))
       CASE(3)
          Mssg='Max abs errors'//'dE3='//TRIM(DblToShrtChar(RelErrProp))//', '        &
                               //'dP3='//TRIM(DblToShrtChar(AbsErrPPrim))
       CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
       END SELECT
       CALL MondoLog(DEBUG_MAXIMUM, Prog, Mssg)

       SELECT CASE(RespOrder)
       CASE(1)
          Mssg='Rel dE1='//TRIM(DblToShrtChar(RelErrProp))//', '    &
             //'||dP1||_F='//TRIM(DblToShrtChar(FNormErrPPrim))
       CASE(2)
          Mssg='Rel dE2='//TRIM(DblToShrtChar(RelErrProp))//', '    &
             //'||dP2||_F='//TRIM(DblToShrtChar(FNormErrPPrim))
       CASE(3)
          Mssg='Rel dE3='//TRIM(DblToShrtChar(RelErrProp))//', '    &
             //'||dP3||_F='//TRIM(DblToShrtChar(FNormErrPPrim))
       CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
       END SELECT
       CALL MondoLog(DEBUG_MAXIMUM, Prog, Mssg)

#ifdef PARALLEL
    ENDIF
#endif
    ENDIF
    !
  END FUNCTION CnvrgChckPrim
  !
  !
  SUBROUTINE Order0(P0,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    !
    CALL Multiply(P0,P0,Tmp3)
    IF(N.GE.Ne) THEN
       CALL Filter(P0,Tmp3)
    ELSE
       CALL Multiply(Tmp3,-One)
       CALL Multiply(P0,Two)
       CALL Add(P0,Tmp3,Tmp1)
       CALL Filter(P0,Tmp1)
    ENDIF
  END SUBROUTINE Order0
  !
  !
  SUBROUTINE Order1_A(P0,P1_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! Da = {Da,X0}
    !
    CALL Multiply(P0,P1_1,Tmp1)
    CALL Multiply(P1_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    IF(N.GE.Ne) THEN
       CALL Filter(P1_1,Tmp3)
    ELSE
       CALL Multiply(Tmp3,-One)
       CALL Multiply(P1_1,Two)
       CALL Add(P1_1,Tmp3,Tmp1)
       CALL Filter(P1_1,Tmp1)
    ENDIF
  END SUBROUTINE Order1_A
  !
  !
  SUBROUTINE Order2_AB(P0,P1_1,P1_2,P2_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P1_2,P2_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P1_2,P2_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! Dab = {Dab,X0} + 1/2{Da,Db}
    !
    CALL Multiply(P0,P2_1,Tmp1)
    CALL Multiply(P2_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    CALL Multiply(P1_1,P1_2,Tmp1)
    CALL Multiply(Tmp1,0.5d0)
    CALL Add(Tmp1,Tmp3,Tmp2)
    CALL Multiply(P1_2,P1_1,Tmp1)
    CALL Multiply(Tmp1,0.5d0)
    CALL Add(Tmp1,Tmp2,Tmp3)
    IF(N.GE.Ne) THEN
       CALL Filter(P2_1,Tmp3)
    ELSE
       CALL Multiply(Tmp3,-One)
       CALL Multiply(P2_1,Two)
       CALL Add(P2_1,Tmp3,Tmp1)
       CALL Filter(P2_1,Tmp1)
    ENDIF
  END SUBROUTINE Order2_AB
  !
  !
  SUBROUTINE Order2_AA(P0,P1_1,P2_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P2_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P2_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! Daa = {Daa,X0} + Da*Da
    !
    CALL Multiply(P0,P2_1,Tmp1)
    CALL Multiply(P2_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    CALL Multiply(P1_1,P1_1,Tmp1)
    CALL Add(Tmp1,Tmp3,Tmp2)
    IF(N.GE.Ne) THEN
       CALL Filter(P2_1,Tmp2)
    ELSE
       CALL Multiply(Tmp2,-One)
       CALL Multiply(P2_1,Two)
       CALL Add(P2_1,Tmp2,Tmp1)
       CALL Filter(P2_1,Tmp1)
    ENDIF
  END SUBROUTINE Order2_AA
  !
  !
  SUBROUTINE Order3_ABC(P0,P1_1,P1_2,P1_3,P2_1,P2_2,P2_3,P3_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P1_2,P1_3,P2_1,P2_2,P2_3,P3_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P1_2,P1_3,P2_1,P2_2,P2_3,P3_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! Dabc = {Dabc,X0} + 1/3({Dab,Dc} + {Dac,Db} + {Dbc,Da})
    ! D0*Dabc+Dabc*D0
    CALL Multiply(P0,P3_1,Tmp1)
    CALL Multiply(P3_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! Da*Dbc
    CALL Multiply(P1_1,P2_3,Tmp2)
    CALL Multiply(Tmp2,1.0d0/3.0d0)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! Dbc*Da
    CALL Multiply(P2_3,P1_1,Tmp2)
    CALL Multiply(Tmp2,1.0d0/3.0d0)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! Db*Dac
    CALL Multiply(P1_2,P2_2,Tmp2)
    CALL Multiply(Tmp2,1.0d0/3.0d0)
    CALL Add(Tmp3,Tmp2,Tmp1)
    ! Dac*Db
    CALL Multiply(P2_2,P1_2,Tmp2)
    CALL Multiply(Tmp2,1.0d0/3.0d0)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! Dc*Dab
    CALL Multiply(P1_3,P2_1,Tmp2)
    CALL Multiply(Tmp2,1.0d0/3.0d0)
    CALL Add(Tmp3,Tmp2,Tmp1)
    ! Dab*Dc
    CALL Multiply(P2_1,P1_3,Tmp2)
    CALL Multiply(Tmp2,1.0d0/3.0d0)
    CALL Add(Tmp1,Tmp2,Tmp3)
    IF(N.GE.Ne) THEN
       CALL Filter(P3_1,Tmp3)
    ELSE
       CALL Multiply(Tmp3,-One)
       CALL Multiply(P3_1,Two)
       CALL Add(P3_1,Tmp3,Tmp1)
       CALL Filter(P3_1,Tmp1)
    ENDIF
  END SUBROUTINE Order3_ABC
  !
  !
  SUBROUTINE Order3_AAA(P0,P1_1,P2_1,P3_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! Daaa = {Daaa,X0} + {Daa,Da}
    ! D0*Daaa+Daaa*D0
    CALL Multiply(P0,P3_1,Tmp1)
    CALL Multiply(P3_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! Da*Daa
    CALL Multiply(P1_1,P2_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! Daa*Da
    CALL Multiply(P2_1,P1_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    IF(N.GE.Ne) THEN
       CALL Filter(P3_1,Tmp3)
    ELSE
       CALL Multiply(Tmp3,-One)
       CALL Multiply(P3_1,Two)
       CALL Add(P3_1,Tmp3,Tmp1)
       CALL Filter(P3_1,Tmp1)
    ENDIF
  END SUBROUTINE Order3_AAA
  !
  !
  SUBROUTINE Order3_AAC(P0,P1_1,P1_3,P2_1,P2_2,P3_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P1_3,P2_1,P2_2,P3_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P1_3,P2_1,P2_2,P3_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! Daac = {Daac,X0} + 1/3{Daa,Dc} + 2/3{Dac,Da}
    ! D0*Daac+Daac*D0
    CALL Multiply(P0,P3_1,Tmp1)
    CALL Multiply(P3_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! Da*Dac
    CALL Multiply(P1_1,P2_2,Tmp2)
    CALL Multiply(Tmp2,2.0d0/3.0d0)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! Dac*Da
    CALL Multiply(P2_2,P1_1,Tmp2)
    CALL Multiply(Tmp2,2.0d0/3.0d0)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! Dc*Daa
    CALL Multiply(P1_3,P2_1,Tmp2)
    CALL Multiply(Tmp2,1.0d0/3.0d0)
    CALL Add(Tmp3,Tmp2,Tmp1)
    ! Daa*Dc
    CALL Multiply(P2_1,P1_3,Tmp2)
    CALL Multiply(Tmp2,1.0d0/3.0d0)
    CALL Add(Tmp1,Tmp2,Tmp3)
    IF(N.GE.Ne) THEN
       CALL Filter(P3_1,Tmp3)
    ELSE
       CALL Multiply(Tmp3,-One)
       CALL Multiply(P3_1,Two)
       CALL Add(P3_1,Tmp3,Tmp1)
       CALL Filter(P3_1,Tmp1)
    ENDIF
  END SUBROUTINE Order3_AAC
  !
  !
  SUBROUTINE Order3_ABB(P0,P1_1,P1_2,P2_1,P2_3,P3_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P1_2,P2_1,P2_3,P3_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P1_2,P2_1,P2_3,P3_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! Dabb = {Dabb,X0}+2/3{Dab,Db}+1/3{Dbb,Da}
    ! D0*Dabb+Dabb*D0
    CALL Multiply(P0,P3_1,Tmp1)
    CALL Multiply(P3_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! Db*Dab
    CALL Multiply(P1_2,P2_1,Tmp2)
    CALL Multiply(Tmp2,2.0d0/3.0d0)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! Dab*Db
    CALL Multiply(P2_1,P1_2,Tmp2)
    CALL Multiply(Tmp2,2.0d0/3.0d0)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! Da*Dbb
    CALL Multiply(P1_1,P2_3,Tmp2)
    CALL Multiply(Tmp2,1.0d0/3.0d0)
    CALL Add(Tmp3,Tmp2,Tmp1)
    ! Dbb*Da
    CALL Multiply(P2_3,P1_1,Tmp2)
    CALL Multiply(Tmp2,1.0d0/3.0d0)
    CALL Add(Tmp1,Tmp2,Tmp3)
    IF(N.GE.Ne) THEN
       CALL Filter(P3_1,Tmp3)
    ELSE
       CALL Multiply(Tmp3,-One)
       CALL Multiply(P3_1,Two)
       CALL Add(P3_1,Tmp3,Tmp1)
       CALL Filter(P3_1,Tmp1)
    ENDIF
  END SUBROUTINE Order3_ABB
  !
  !
  SUBROUTINE Order4(P0,P1_1,P2_1,P3_1,P4_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! Daaaa = {Daaaa,X0} + {Daaa,Da} + Daa*Daa
    ! D0*Daaaa+Daaaa*D0
    CALL Multiply(P0,P4_1,Tmp1)
    CALL Multiply(P4_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! Da*Daaa
    CALL Multiply(P1_1,P3_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! Daaa*Da
    CALL Multiply(P3_1,P1_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! Daa*Daa
    CALL Multiply(P2_1,P2_1,Tmp2)
    CALL Add(Tmp3,Tmp2,Tmp1)
    IF(N.GE.Ne) THEN
       CALL Filter(P4_1,Tmp1)
    ELSE
       CALL Multiply(Tmp1,-One)
       CALL Multiply(P4_1,Two)
       CALL Add(P4_1,Tmp1,Tmp2)
       CALL Filter(P3_1,Tmp2)
    ENDIF
  END SUBROUTINE Order4
  !
  !
  SUBROUTINE Order5(P0,P1_1,P2_1,P3_1,P4_1,P5_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! D5 = {D5,X0} + {D4,D1} + {D3,D2}
    ! D0*D5+D5*D0
    CALL Multiply(P0,P5_1,Tmp1)
    CALL Multiply(P5_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D1*D4
    CALL Multiply(P1_1,P4_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D4*D1
    CALL Multiply(P4_1,P1_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D2*D3
    CALL Multiply(P2_1,P3_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D3*D2
    CALL Multiply(P3_1,P2_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    IF(N.GE.Ne) THEN
       CALL Filter(P5_1,Tmp3)
    ELSE
       CALL Multiply(Tmp3,-One)
       CALL Multiply(P5_1,Two)
       CALL Add(P5_1,Tmp3,Tmp2)
       CALL Filter(P5_1,Tmp2)
    ENDIF
  END SUBROUTINE Order5
  !
  !
  SUBROUTINE Order6(P0,P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! D6 = {D6,X0} + {D5,D1} + {D4,D2} + D3*D3
    ! D0*D6+D6*D0
    CALL Multiply(P0,P6_1,Tmp1)
    CALL Multiply(P6_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D1*D5
    CALL Multiply(P1_1,P5_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D5*D1
    CALL Multiply(P5_1,P1_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D2*D4
    CALL Multiply(P2_1,P4_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D4*D2
    CALL Multiply(P4_1,P2_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D3*D3
    CALL Multiply(P3_1,P3_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    IF(N.GE.Ne) THEN
       CALL Filter(P6_1,Tmp1)
    ELSE
       CALL Multiply(Tmp1,-One)
       CALL Multiply(P6_1,Two)
       CALL Add(P6_1,Tmp1,Tmp2)
       CALL Filter(P6_1,Tmp2)
    ENDIF
  END SUBROUTINE Order6
  !
  !
  SUBROUTINE Order7(P0,P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! D7 = {D7,X0} + {D6,D1} + {D5,D2} + {D4,D3}
    ! D0*D7+D7*D0
    CALL Multiply(P0,P7_1,Tmp1)
    CALL Multiply(P7_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D1*D6
    CALL Multiply(P1_1,P6_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D6*D1
    CALL Multiply(P6_1,P1_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D2*D5
    CALL Multiply(P2_1,P5_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D5*D2
    CALL Multiply(P5_1,P2_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D3*D4
    CALL Multiply(P3_1,P4_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D4*D3
    CALL Multiply(P4_1,P3_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    IF(N.GE.Ne) THEN
       CALL Filter(P7_1,Tmp3)
    ELSE
       CALL Multiply(Tmp3,-One)
       CALL Multiply(P7_1,Two)
       CALL Add(P7_1,Tmp3,Tmp2)
       CALL Filter(P7_1,Tmp2)
    ENDIF
  END SUBROUTINE Order7
  !
  !
  SUBROUTINE Order8(P0,P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,P8_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,P8_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,P8_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! D8 = {D8,X0} + {D7,D1} + {D6,D2} + {D5,D3} + D4*D4
    ! D0*D8+D8*D0
    CALL Multiply(P0,P8_1,Tmp1)
    CALL Multiply(P8_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D1*D7
    CALL Multiply(P1_1,P7_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D7*D1
    CALL Multiply(P7_1,P1_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D2*D6
    CALL Multiply(P2_1,P6_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D6*D2
    CALL Multiply(P6_1,P2_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D3*D5
    CALL Multiply(P3_1,P5_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D5*D3
    CALL Multiply(P5_1,P3_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D4*D4
    CALL Multiply(P4_1,P4_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    IF(N.GE.Ne) THEN
       CALL Filter(P8_1,Tmp1)
    ELSE
       CALL Multiply(Tmp1,-One)
       CALL Multiply(P8_1,Two)
       CALL Add(P8_1,Tmp1,Tmp2)
       CALL Filter(P8_1,Tmp2)
    ENDIF
  END SUBROUTINE Order8
  !
  !
  SUBROUTINE Order9(P0,P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,P8_1,P9_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,P8_1,P9_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,P8_1,P9_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! D9 = {D9,X0} + {D8,D1} + {D7,D2} + {D6,D3} + {D5,D4}
    ! D0*D9+D9*D0
    CALL Multiply(P0,P9_1,Tmp1)
    CALL Multiply(P9_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D1*D8
    CALL Multiply(P1_1,P8_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D8*D1
    CALL Multiply(P8_1,P1_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D2*D7
    CALL Multiply(P2_1,P7_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D7*D2
    CALL Multiply(P7_1,P2_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D3*D6
    CALL Multiply(P3_1,P6_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D6*D3
    CALL Multiply(P6_1,P3_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D4*D5
    CALL Multiply(P4_1,P5_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D5*D4
    CALL Multiply(P5_1,P4_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    IF(N.GE.Ne) THEN
       CALL Filter(P9_1,Tmp3)
    ELSE
       CALL Multiply(Tmp3,-One)
       CALL Multiply(P9_1,Two)
       CALL Add(P9_1,Tmp3,Tmp2)
       CALL Filter(P9_1,Tmp2)
    ENDIF
  END SUBROUTINE Order9
  !
  !
  SUBROUTINE Order10(P0,P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,P8_1,P9_1,P10_1,Tmp1,Tmp2,Tmp3,N,Ne)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)  , INTENT(INOUT) :: P0
    TYPE(DBCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,P8_1,P9_1,P10_1,Tmp1,Tmp2,Tmp3
#else
    TYPE( BCSR)  , INTENT(INOUT) :: P0
    TYPE( BCSR)  , INTENT(INOUT) :: P1_1,P2_1,P3_1,P4_1,P5_1,P6_1,P7_1,P8_1,P9_1,P10_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE) , INTENT(IN   ) :: N,Ne
    ! D10 = {D10,X0} + {D9,D1} + {D8,D2} + {D7,D3} + {D6,D4} + D5*D5
    ! D0*D10+D10*D0
    CALL Multiply(P0,P10_1,Tmp1)
    CALL Multiply(P10_1,P0,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D1*D9
    CALL Multiply(P1_1,P9_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D9*D1
    CALL Multiply(P9_1,P1_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D2*D8
    CALL Multiply(P2_1,P8_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D8*D2
    CALL Multiply(P8_1,P2_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D3*D7
    CALL Multiply(P3_1,P7_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D7*D3
    CALL Multiply(P7_1,P3_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D4*D6
    CALL Multiply(P4_1,P6_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    ! D6*D4
    CALL Multiply(P6_1,P4_1,Tmp2)
    CALL Add(Tmp1,Tmp2,Tmp3)
    ! D5*D5
    CALL Multiply(P5_1,P5_1,Tmp2)
    CALL Add(Tmp2,Tmp3,Tmp1)
    IF(N.GE.Ne) THEN
       CALL Filter(P10_1,Tmp1)
    ELSE
       CALL Multiply(Tmp1,-One)
       CALL Multiply(P10_1,Two)
       CALL Add(P10_1,Tmp1,Tmp2)
       CALL Filter(P10_1,Tmp2)
    ENDIF
  END SUBROUTINE Order10
  !
  !
#ifdef TC2R_EIGENVAL
  !
  SUBROUTINE PLOT_ERROR_WITH_BOXES(Eig,NBrBox,Name)
    TYPE(DBL_VECT)               :: Eig
    INTEGER         , INTENT(IN) :: NBrBox
    CHARACTER(LEN=*), INTENT(IN) :: Name
    TYPE(INT_VECT)               :: Count
    REAL(DOUBLE)                 :: MaxE,MinE,DeltaE
    INTEGER                      :: ICount,IBox,IBas,Plt
    !
    CALL New(Count,NBrBox)
    !
    MaxE=MAXVAL(Eig%D(1:NBasF))
    MinE=MINVAL(Eig%D(1:NBasF))
    !
    DeltaE=(MaxE-MinE)/DBLE(NBrBox)
    !
    DO IBox=1,NBrBox
       ICount=0
       DO IBas=1,NBasF
          IF(MinE+DBLE((IBox-1))*DeltaE.LE.Eig%D(IBas).AND. &
                 Eig%D(IBas).LE.MinE+DBLE(IBox)*DeltaE) ICount=ICount+1
       ENDDO
       Count%I(IBox)=ICount
    ENDDO
    !
!    IF(SUM(Count%I).NE.NBasF) STOP 11111
    !
    CALL OpenASCII(TRIM(Name)//'_PLOT_BAR_data',Plt,NewFile_O=.TRUE.)
    DO IBox=1,NBrBox
       WRITE(Plt,*) MinE+(DBLE(IBox-1)+0.5d0)*DeltaE ,Count%I(IBox)
    ENDDO
    CLOSE(Plt)
    !
    CALL OpenASCII(TRIM(Name)//'_PLOT_BAR_PlotMe',Plt,NewFile_O=.TRUE.)
    WRITE(Plt,2)
    WRITE(Plt,3)TRIM(Name)//'.eps'
    WRITE(Plt,6)
    WRITE(Plt,*)'set pointsize 0.1'
!    WRITE(Plt,*)'set logscale y'
    WRITE(Plt,*)"plot '"//TRIM(Name)//"_PLOT_BAR_data' using 1:2 notitle with boxes "
    CLOSE(Plt)
    !
    CALL Delete(Count)
    !
   1   FORMAT(2(1x,I16))
   2   FORMAT('set term  postscript eps  "Times-Roman" 18')
!   2   FORMAT('set term  jpeg transparent')
   3   FORMAT('set output "',A,'"')
   6   FORMAT('set size square ')
   9   FORMAT('plot [0 : ',I12,' ] ',I12,', \\')
  10   FORMAT('                    ',I12,', \\')
  END SUBROUTINE PLOT_ERROR_WITH_BOXES
#endif
  !
  !
END MODULE DenMatResponse
