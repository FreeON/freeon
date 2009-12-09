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
PROGRAM FockNGrueven
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE SetXYZ
  USE LinAlg
  USE Functionals
  USE McMurchie
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)                    :: F,V,T,J,K,X,Tmp1
#else
  TYPE(BCSR)                     :: F,V,T,J,K,X,Tmp1
#endif
  TYPE(ARGMT)                    :: Args
  TYPE(CRDS)                     :: GM
  TYPE(BSET)                     :: BS
  REAL(DOUBLE)                   :: KScale,E_NukeNuke
  INTEGER                        :: ISCF,I,RespOrder
  CHARACTER(LEN=2)               :: Cycl
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: XFile,DevFile,Mssg,ResponsePostFix
  CHARACTER(LEN=12),PARAMETER    :: Prog='FockNGrueven'
  LOGICAL                        :: Present,ExchangeShift,HasECPs
!------------------------------------------------------------------
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  ISCF=Args%I%I(1)
  Cycl=IntToChar(ISCF)

  CALL Initialize(F)
  CALL Initialize(V)
  CALL Initialize(T)
  CALL Initialize(J)
  CALL Initialize(K)
  CALL Initialize(X)
  CALL Initialize(Tmp1)

  IF(SCFActn=='StartResponse')THEN
     RespOrder=LEN(TRIM(Args%C%C(3)))
     SELECT CASE(RespOrder)
     CASE(1)
        ! At this moment F'=mu.
        ResponsePostFix=TRIM(Args%C%C(4))//TRIM(Args%C%C(5))
        CALL Get(F,TrixFile(ResponsePostFix,Args))            ! T=M_{x} or whatever...
        !
        ! Transform also the dipole matrix for TC2.
        XFile=TrixFile('X',Args)
        INQUIRE(FILE=XFile,EXIST=Present)

        IF(Present)THEN
           CALL Get(X,XFile)                ! X =S^{-1/2}
           CALL Multiply(X,F,Tmp1)          ! T1=S^{-1/2}.F
           CALL Multiply(Tmp1,X,F)          ! T1=S^{-1/2}.F.S^{-1/2}
        ELSE
           CALL Get(X,TrixFile('ZT',Args))  ! X =Z^t=L^{-t}
           CALL Multiply(X,F,Tmp1)          ! T1=Z^t.F_AO
           CALL Get(X,TrixFile('Z',Args))   ! X=Z=L^{-1}
           CALL Multiply(Tmp1,X,F)          ! F=Z^t.F_AO.Z
        ENDIF
        CALL Filter(Tmp1,F)                 ! T1 =F_Orthog=Filter[Z^t.F_AO.Z]
        !
        CALL Put(Tmp1,TrixFile('Ortho'//TRIM(Args%C%C(4))//TRIM(Args%C%C(5)),Args))
        CALL PChkSum(Tmp1,'Ortho'//TRIM(Args%C%C(4))//TRIM(Args%C%C(5)),Prog)
        CALL PPrint(Tmp1,'Ortho'//TRIM(Args%C%C(4))//TRIM(Args%C%C(5)))
        CALL Plot(Tmp1,'Ortho'//TRIM(Args%C%C(4))//TRIM(Args%C%C(5)))
        !
        ! Reload a fresh matrix.
        CALL Get(F,TrixFile(ResponsePostFix,Args))            ! T=M_{x} or whatever...
        !
     CASE(2)
        ! Set F^2 to zero.
        ResponsePostFix=TRIM(Args%C%C(4))//TRIM(Args%C%C(5))
        CALL Get(F,TrixFile(ResponsePostFix,Args))
        CALL Multiply(F,Zero)
     CASE(3)
        ! Set F^3 to zero.
        ResponsePostFix=TRIM(Args%C%C(4))//TRIM(Args%C%C(5))
        CALL Get(F,TrixFile(ResponsePostFix,Args))
        CALL Multiply(F,Zero)
     CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
     END SELECT
  ELSEIF(SCFActn=='FockPrimeBuild') THEN
     RespOrder=LEN(TRIM(Args%C%C(3)))
     ! Start with the Coulomb matrix
     CALL Get(J,TrixFile('JPrime'//TRIM(Args%C%C(3)),Args,0))! J=J_Coulomb[Rho_Total(Nuc+El)]
     ! If we are doing a response calc, then we add in just the multpole matrix or whatever...
     SELECT CASE(RespOrder)
     CASE(1)
        ResponsePostFix=TRIM(Args%C%C(4))//TRIM(Args%C%C(5))
        CALL Get(T,TrixFile(ResponsePostFix,Args))            ! T=M_{x} or whatever...
        CALL Add(J,T,F)                                       ! F=J+T
     CASE(2); CALL SetEq(F,J)                                 ! F=J
     CASE(3); CALL SetEq(F,J)                                 ! F=J
     CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
     END SELECT
  ELSE
     ! Start with the Coulomb matrix
     CALL Get(J,TrixFile('J',Args,0))                      ! J=J_Coulomb[Rho_Total(Nuc+El)]
     ! Doing normal SCF, add in the kinetic energy
     CALL Get(T,TrixFile('T',Args))                        ! T=T_{Kinetic}
     CALL Add(J,T,F)                                       ! F=J+T
  ENDIF
  ! Check to see if we have ECPs
  CALL Get(HasECPs,'hasecps',Tag_O=CurBase)
  IF(HasECPs)THEN
     CALL SetEq(J,F)
     ! Get the pseudopotential matrix U
     CALL Get(T,TrixFile('U',Args))                        ! T=U_{ECP}
     CALL Add(T,J,F)
  ENDIF

  IF(AllocQ(T%Alloc)) CALL Delete(T)
  IF(AllocQ(J%Alloc)) CALL Delete(J)

  IF(SCFActn/='GuessEqCore'.AND.SCFActn.NE.'StartResponse')THEN
     IF(HasHF(ModelChem).AND.HasDFT(ModelChem))THEN
        ! Add in Hartree-Fock exact exchange
        IF(SCFActn=='FockPrimeBuild') THEN
           CALL Halt('Response with DFT does not work yet.')
           CALL Get(K,TrixFile('KPrime'//TRIM(Args%C%C(3)),Args,0))   ! K=K_hf
        ELSE
           CALL Get(K,TrixFile('K',Args,0))                           ! K=K_hf
        ENDIF
        KScale=ExactXScale(ModelChem)
        CALL Multiply(K,KScale)                                       ! K=KScale*K_hf
        CALL Add(F,K,Tmp1)                                            ! F=J+T+K
        ! Now add in exchange-correlation
        IF(SCFActn=='FockPrimeBuild') THEN
           CALL Halt('Response with DFT does not work yet.')
           CALL Get(K,TrixFile('KxcPrime'//TRIM(Args%C%C(3)),Args,0)) ! K=K_xc
        ELSE
           CALL Get(K,TrixFile('Kxc',Args,0))                         ! K=K_xc
        ENDIF
        CALL Add(K,Tmp1,F)                                            ! F=J+T+KShift*K_hf+K_xc
     ELSEIF(HasHF(ModelChem))THEN
        !  Add in Hartree-Fock exact exchange
        IF(SCFActn=='FockPrimeBuild')THEN
           CALL Get(K,TrixFile('KPrime'//TRIM(Args%C%C(3)),Args,0))   ! K=K_hf
        ELSE
           CALL Get(K,TrixFile('K',Args,0))                           ! K=K_hf
        ENDIF
        CALL Add(F,K,Tmp1)                                            ! F=J+T+K_hf
        CALL SetEq(F,Tmp1)
     ELSEIF(HasDFT(ModelChem))THEN
        !  Add in exchange-correlation
        IF(SCFActn=='FockPrimeBuild') THEN
           CALL Halt('Response with DFT does not work yet.')
           CALL Get(K,TrixFile('KxcPrime'//TRIM(Args%C%C(3)),Args,0)) ! K=K_xc
        ELSE
           CALL Get(K,TrixFile('Kxc',Args,0))                         ! K=K_{xc}
        ENDIF
        CALL Add(K,F,Tmp1)                                            ! F=J+T+K_xc
        CALL SetEq(F,Tmp1)
     ENDIF
     CALL Delete(K)
  ENDIF
!
  IF(SCFActn=='FockPrimeBuild'.OR.SCFActn=='StartResponse') THEN
     !I do not need that! (vw)
     CALL Put(F,TrixFile('FPrime'//TRIM(Args%C%C(3)),Args,0))
     CALL PChkSum(F,'FPrime'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( F,'FPrime'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
     !CALL Print_BCSR(F,'FPrime'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']',Unit_O=6)
     CALL Plot(   F,'FPrime'//TRIM(Args%C%C(3))//'_'//TRIM(SCFCycl))
  ELSE
     CALL Put(F,TrixFile('F',Args,0))
     CALL PChkSum(F,'F['//TRIM(Cycl)//']',Prog)
     !  CALL PPrint( F,'F['//TRIM(Cycl)//']',Unit_O=6)
     CALL PPrint( F,'F['//TRIM(Cycl)//']')
     CALL Plot(   F,'F['//TRIM(Cycl)//']')
  ENDIF
! Now transform to an orthogonal representation
  XFile=TrixFile('X',Args)
  INQUIRE(FILE=XFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(X,XFile)                ! X =S^{-1/2}
     CALL Multiply(X,F,Tmp1)          ! T1=S^{-1/2}.F
     CALL Multiply(Tmp1,X,F)          ! T1=S^{-1/2}.F.S^{-1/2}
  ELSE
     CALL Get(X,TrixFile('ZT',Args))  ! X =Z^t=L^{-t}
     CALL Multiply(X,F,Tmp1)          ! T1=Z^t.F_AO
     CALL Get(X,TrixFile('Z',Args))   ! X=Z=L^{-1}
     CALL Multiply(Tmp1,X,F)          ! F=Z^t.F_AO.Z
  ENDIF
  CALL Filter(Tmp1,F)                 ! T1 =F_Orthog=Filter[Z^t.F_AO.Z]
!
  IF(SCFActn=='FockPrimeBuild'.OR.SCFActn=='StartResponse') THEN
     CALL Put(Tmp1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)),Args,0))
     CALL PChkSum(Tmp1,'OrthoFPrime'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( Tmp1,'OrthoFPrime'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
     CALL Plot(   Tmp1,'OrthoFPrime'//TRIM(Args%C%C(3))//'_'//TRIM(SCFCycl))
  ELSE
     CALL Put(Tmp1,TrixFile('OrthoF',Args,0))
     CALL PChkSum(Tmp1,'OrthoF['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint(Tmp1,'OrthoF['//TRIM(SCFCycl)//']')
     CALL Plot(Tmp1,'OrthoF_'//TRIM(SCFCycl))
  ENDIF
! Tidy up
  CALL Delete(F)
  CALL Delete(Tmp1)
  CALL ShutDown(Prog)
END PROGRAM FockNGrueven
