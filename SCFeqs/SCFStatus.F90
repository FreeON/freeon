!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    DETERMINE THE STATUS OF AN SCF CYCLE 
!
PROGRAM SCFStatus
!--------------------------------------- 
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE Functionals
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
  TYPE(ARGMT)                     :: Args
#ifdef PARALLEL
  TYPE(DBCSR)                     :: P,Tmp1,Tmp2
#else
  TYPE(BCSR)                      :: P,Tmp1,Tmp2
#endif
  REAL(DOUBLE)                    :: E_el_tot,E_nuc_tot,E_es_tot,KinE,ExchE,Exc, &
                                     Gap,Etot,DMax,Virial
  CHARACTER(LEN=DEFAULT_CHR_LEN)  :: SCFString,PrevBSetName,CurrBSetName
  CHARACTER(LEN=9),PARAMETER      :: Prog='SCFStatus'  
!---------------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog)
! Allocate some matrices
  CALL New(P)
  CALL New(Tmp1)
  CALL New(Tmp2)
!---------------------------------------------
! Get the density matrix 
  IF(SCFActn=='Switch')THEN 
!    if switching the density matrix, 
!    this is all that is available ...    
     CALL Get(P,TrixFile('D',Args,1))
  ELSE
     CALL Get(P,TrixFile('D',Args,0))
  ENDIF
!---------------------------------------------
! KinE=<T>=Tr{P.T}
  CALL Get(Tmp1,TrixFile('T',Args))
#ifdef PARALLEL
  CALL Multiply(P,Tmp1,Tmp2)
  KinE=Two*Trace(Tmp2)     
#else
  KinE=Two*Trace(P,Tmp1)    
#endif 
  CALL Put(KinE,'Kin',Tag_O=SCFCycl)
!---------------------------------------------
! E_el_tot=<Vee+Vne>=Tr{P.(Vee+Vne)}
  CALL Get(Tmp1,TrixFile('J',Args,0))
#ifdef PARALLEL
  CALL Multiply(P,Tmp1,Tmp2)
  E_el_tot=Trace(Tmp2)     
#else
  E_el_tot=Trace(P,Tmp1)    
#endif
  CALL Put(E_el_tot,'ene+eee',Tag_O=SCFCycl)
!---------------------------------------------
! ExchE=<K>=Tr{P.K}
  IF(HasHF(ModelChem))THEN
     CALL Get(Tmp1,TrixFile('K',Args,0))
#ifdef PARALLEL
     CALL Multiply(P,Tmp1,Tmp2)
     ExchE=ExactXScale(ModelChem)*Trace(Tmp2)     
#else
     ExchE=ExactXScale(ModelChem)*Trace(P,Tmp1)    
#endif
  ELSE
     ExchE=Zero
  ENDIF
  CALL Put(ExchE,'Ex',Tag_O=SCFCycl)
!----------------------------------------------
! Get E_nuc_tot =<Vnn+Vne> 
  CALL Get(E_nuc_tot,'enn+ene',Tag_O=SCFCycl)
!----------------------------------------------
! Get the exchange correlation energy
  CALL Get(ModelChem,'ModelChemistry',CurBase)
  IF(HasDFT(ModelChem))THEN      
     CALL Get(Exc,'Exc',Tag_O=SCFCycl)
  ELSE
     Exc=Zero
  ENDIF
  E_es_tot = E_el_tot+E_nuc_tot
  Etot=KinE+E_es_tot+Exc+ExchE
  CALL Put(Etot,'E_total',Tag_O=SCFCycl)
!----------------------------------------------
! Calculate the Viraial Coefficient
  Virial=E_es_tot/KinE
!--------------------------------------------------------
!Find the largest block of the delta density matrix
!
 DMax=BIG_DBL
 IF(Current(2)==1.AND.Current(1)>1)THEN
    CALL Get(Tmp1,TrixFile('OrthoD',Args,0))
    CALL Get(Tmp2,TrixFile('OrthoD',Args,1))
    CALL Multiply(Tmp1,-One)
    CALL Add(Tmp1,Tmp2,P)
    DMax=Max(P)
 ELSEIF(Current(2)==Previous(2))THEN
    CALL Get(PrevBSetName,'bsetname',Tag_O=PrvBase)
    CALL Get(CurrBSetName,'bsetname',Tag_O=CurBase)
    IF(TRIM(PrevBSetName)==TRIM(CurrBSetName))THEN
       CALL Get(Tmp1,TrixFile('OrthoD',Args,0))!,Stats_O=Previous))
       CALL Get(Tmp2,TrixFile('OrthoD',Args,1))!,Stats_O=Current))
       CALL Multiply(Tmp1,-One)
       CALL Add(Tmp1,Tmp2,P)
       DMax=Max(P)
    ENDIF
 ENDIF
!
 CALL Put(DMax,'DMax',Tag_O='_'//TRIM(CurGeom) &
                          //'_'//TRIM(CurBase) &
                          //'_'//TRIM(SCFCycl))

 CALL Put(Etot,'Etot',Tag_O='_'//TRIM(CurGeom) & 
                          //'_'//TRIM(CurBase) &
                          //'_'//TRIM(SCFCycl))
!
 CALL Get(Gap,'HomoLumoGap')
 SCFString='================== Geometry #'//TRIM(CurGeom) &
                           //', MondoSCF #'//TRIM(SCFCycl)//' ==================== '
#ifdef PARALLEL
 IF(MyId==ROOT)THEN
#endif
    WRITE(*,* )TRIM(SCFString)
    WRITE(*,10) KinE     ,RTRN, &
                E_es_tot ,RTRN, &
                ExchE    ,RTRN, &
                Exc      ,RTRN, &
                Etot     ,RTRN, &
                Gap      ,RTRN
    IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
       CALL OpenASCII(OutFile,Out)
       IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN

!          WRITE(Out,9) SCFCycl                                ,RTRN, &
!                       FRACTION(KinE)     ,EXPONENT(KinE)     ,RTRN, &
!                       FRACTION(E_el_tot) ,EXPONENT(E_el_tot) ,RTRN, &
!                       FRACTION(E_nuc_tot),EXPONENT(E_nuc_tot),RTRN, &
!                       FRACTION(E_es_tot) ,EXPONENT(E_es_tot) ,RTRN, &
!                       FRACTION(ExchE)    ,EXPONENT(ExchE)    ,RTRN, &
!                       FRACTION(Exc)      ,EXPONENT(Exc)      ,RTRN, &
!                       FRACTION(Etot)     ,EXPONENT(Etot)
          9 FORMAT('(*================== MondoSCF[',I3,']=====================*) ',A1,&
                   '  TrKin[',I3,' ] = ',F19.16,'*2^(',I4,');',A1, &
                   '  TrVel[',I3,' ] = ',F19.16,'*2^(',I4,');',A1, &
                   '  TrVnu[',I3,' ] = ',F19.16,'*2^(',I4,');',A1, &
                   '  TrVnu[',I3,' ] = ',F19.16,'*2^(',I4,');',A1, &
                   '  TrKhf[',I3,' ] = ',F19.16,'*2^(',I4,');',A1, & 
                   '  TrKdft[',I3,'] = ',F19.16,'*2^(',I4,');',A1, & 
                   '  ESCF[',I3,  '] = ',F19.16,'*2^(',I4,');')
       ELSE
          WRITE(Out,* )TRIM(SCFString)
          WRITE(Out,10) KinE     ,RTRN, &
                        E_es_tot ,RTRN, &
                        ExchE    ,RTRN, &
                        Exc      ,RTRN, &
                        Etot     ,RTRN, &
                        Gap      ,RTRN
          10 FORMAT('  <T>       = ',D22.16,A1, &
                    '  <V>       = ',D22.16,A1, &
                    '  <HF>      = ',D22.16,A1, &
                    '  <DFT>     = ',D22.16,A1, &
                    '  <SCF>     = ',D22.16,A1, &
                    '  HOMO-LUMO = ',F8.5  ,A1)
       ENDIF
       CLOSE(Out)
    ENDIF
#ifdef PARALLEL
 ENDIF
#endif

 CALL Delete(P)
 CALL Delete(Tmp1)
 CALL Delete(Tmp2)

 CALL ShutDown(Prog)
END PROGRAM SCFStatus


