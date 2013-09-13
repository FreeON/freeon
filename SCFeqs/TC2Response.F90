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
PROGRAM TC2R
!H=================================================================================
!H PROGRAM TC2Response
!H
!H  OPTIONS:
!H  DEBUGING: Use -DTC2R_DBUG.
!H  INFO    : Use -DTC2R_INFO.
!H  EIGENVAL: Use -DTC2R_EIGENVAL.
!H
!H Comment:
!H
!H
!H=================================================================================
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE Parse
  USE Macros
  USE LinAlg
  USE DenMatMethods , ONLY: SpectralBounds,SetVarThresh
  USE DenMatResponse, ONLY: AllocArray,LoadMatrices,FockPrimGuess,SetPPrmOld,DoTC2R,SavePPrm,DeAllocArray
  !
  IMPLICIT NONE
  !
  !------------------------------------------------------------------
  TYPE(BCSR )                    :: F,FPrm1_1,FPrm1_2,FPrm1_3, &
                                    FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1
#ifdef PARALLEL
  TYPE(DBCSR)                    :: T,P,PPrm1_1,PPrm1_2,PPrm1_3, &
                                    PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
  TYPE(DBCSR)                    :: PPrmOld,Tmp1,Tmp2,Tmp3
#else
  TYPE(BCSR )                    :: T,P,PPrm1_1,PPrm1_2,PPrm1_3, &
                                    PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
  TYPE(BCSR )                    :: PPrmOld,Tmp1,Tmp2,Tmp3
#endif
  TYPE(ARGMT)                    :: Args
  !-------------------------------------------------------------------
  INTEGER                        :: RespOrder
  REAL(DOUBLE)                   :: SwitchThresh
  CHARACTER(LEN=*), PARAMETER    :: Prog='TC2R'
  !-------------------------------------------------------------------
  REAL(DOUBLE)    , PARAMETER    :: DEFAULT_SWITCHTHRESH=1.00D-3
  !-------------------------------------------------------------------
  !
  ! Initial setup.
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
  !
  ! Open input
  CALL OpenASCII(InpFile,Inp)
  !
  ! Switching.
  IF(.NOT.OptDblQ(Inp,'TC2RSwitchingThresh',SwitchThresh)) SwitchThresh=DEFAULT_SWITCHTHRESH
  !
  ! Close the input.
  CLOSE(Inp)
  !
  ! Get the response order.
  RespOrder=LEN(TRIM(Args%C%C(3)))
  !
  !-------------------------------------------------------------------
  ! Allocations.
  !-------------------------------------------------------------------
  !
  CALL AllocArray(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
                  P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                  PPrmOld,T,Tmp1,Tmp2,Tmp3,RespOrder,Args)
  !
  !-------------------------------------------------------------------
  ! Loading matrices.
  !-------------------------------------------------------------------
  !
  CALL LoadMatrices(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
                    P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                    T,RespOrder,Args)
  !
  !-------------------------------------------------------------------
  ! Get Density matrix guess.
  !-------------------------------------------------------------------
  !
  CALL FockPrimGuess(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
                     P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                     RespOrder,Args)
  !
  !-------------------------------------------------------------------
  ! Set up the old densiry derivative.
  !-------------------------------------------------------------------
  !
  CALL SetPPrmOld(PPrmOld,PPrm1_1,PPrm2_1,PPrm3_1,RespOrder)
  !
  !-------------------------------------------------------------------
  ! Let's TC2R.
  !-------------------------------------------------------------------
  !
  CALL DoTC2R(P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1,PPrmOld, &
              T,Tmp1,Tmp2,Tmp3,RespOrder,SwitchThresh,Args)
  !
  !-------------------------------------------------------------------
  ! Save the matrices.
  !-------------------------------------------------------------------
  !
  CALL SavePPrm(PPrm1_1,PPrm2_1,PPrm3_1,Tmp1,Tmp2,RespOrder,Args)
  !
  !-------------------------------------------------------------------
  ! Deallocations the matrices.
  !-------------------------------------------------------------------
  !
  CALL DeAllocArray(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
                    P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
                    PPrmOld,T,Tmp1,Tmp2,Tmp3,RespOrder,Args)
  !
  !-------------------------------------------------------------------
  ! We are done.
  !-------------------------------------------------------------------
  !
  CALL ShutDown(Prog)
  !
END PROGRAM TC2R
