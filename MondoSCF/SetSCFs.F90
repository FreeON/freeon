MODULE SetSCFs
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE ProcessControl
   USE BasisSetParameters   
   USE InOut
   USE Parse
   USE SetXYZ
   USE PrettyPrint
   USE SCFLocals
   USE McMurchie
   USE Overlay
   USE Order
   IMPLICIT NONE
   CHARACTER(LEN=3) :: CSet
   TYPE(BSet)       :: Base
   CONTAINS 
!======================================================================
!     Set up the HDF InfFile, output to ASCII OutFile
!======================================================================
      SUBROUTINE SetSCF(Ctrl)
         TYPE(SCFControls),INTENT(INOUT) :: Ctrl
         INTEGER                         :: ISet
         TYPE(INT_VECT)                  :: Stat
         CALL OpenHDF(Ctrl%Info)         
         IF(Ctrl%Rest)THEN
            CALL New(Stat,3)
            CALL Get(Stat,'PreviousStatus')
            Ctrl%Previous=Stat%I
            CALL Get(Stat,'CurrentStatus')
            Ctrl%Current=Stat%I
            CALL Delete(Stat)
            CALL Get(Ctrl%NSet,'NumberOfSets')
            CALL Get(Ctrl%NGeom,'NumberOfGeometries')
            DO ISet=1,Ctrl%NSet
               CALL Get(Ctrl%AccL(ISet),  'SCFAccuracy',   Tag_O=IntToChar(ISet))
               CALL Get(Ctrl%Model(ISet) ,'ModelChemistry',Tag_O=IntToChar(ISet))
               CALL Get(Ctrl%Method(ISet),'SCFMethod'     ,Tag_O=IntToChar(ISet))
            ENDDO
         ELSE
            DO ISet=1,Ctrl%NSet
               Ctrl%Current(2)=ISet
               CSet=TRIM(IntToChar(ISet))
               CALL Get(Base,Tag_O=CSet)
!              Matrix limits to HDF InfFile
               CALL SetLimits(Ctrl)
!              Primitive distribution stuff to HDF InfFile
               CALL SetDist(Ctrl)   
#ifdef PARALLEL
!              Compute a more sophisticated domain decomposition
!              CALL SetDomainD(Ctrl)
#endif
               CALL Delete(Base)
            ENDDO
         ENDIF
         CALL CloseHDF()
      END SUBROUTINE SetSCF
!=============================================================
!    Set BCSR matrix limits and global auxiliary arrays
!=============================================================
     SUBROUTINE SetLimits(Ctrl,OverRideNProc_O)
         USE GlobalScalars
         TYPE(SCFControls) :: Ctrl
         INTEGER, OPTIONAL :: OverRideNProc_O
         INTEGER           :: I,K,Acc,BWEstim,MaxBaseNode
         CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: Mssg
!-------------------------------------------------------------
         I=Ctrl%Current(2)
         Acc=Ctrl%AccL(I)
         Ctrl%EErr=1.D0
!        Use assymptotics to set the max matrix dimensions
         MaxAtms=1+NAtoms
         BWEstim=MIN(NAtoms, CEILING( (DBLE(NAtoms) &
                +BandWidth(Acc)*BWDecay(Acc)*DBLE(NAtoms)**2) &
                /(One+BWDecay(Acc)*DBLE(NAtoms)**2) ) ) 
         MaxBlks=1+NAtoms*BWEstim
         NBasF=Base%NBasF
         MaxNon0=1+NBasF*(DBLE(NBasF)*DBLE(BWEstim)/DBLE(NAtoms))
!        Put the limits into the HDF InfFile
         CALL Put(MaxAtms,'maxatms',Tag_O=CSet)
         CALL Put(MaxBlks,'maxblks',Tag_O=CSet)
         CALL Put(MaxNon0,'maxnon0',Tag_O=CSet)
#ifdef PARALLEL
         CALL New(BSiz,NAtoms)
         CALL New(OffS,NAtoms)
         CALL Get(BSiz,'atsiz',Tag_O=CSet)
         CALL Get(OffS,'atoff',Tag_O=CSet)
!-----------------------------------------------------------
!        Compute and set limits of the domain decomposition
!
         CALL Decomp(Ctrl,OverRideNProc_O)
!---------------------------------------------------------
!        Calculate max limits per node
!
         MaxAtmsNode=0
         MaxBaseNode=0
         DO K=0,NPrc-1
            MaxAtmsNode=MAX(MaxAtmsNode,End%I(K)-Beg%I(K)+1)
            MaxBaseNode=MAX(MaxBaseNode,                            & 
                            SUM(BSiz%I(Beg%I(K):End%I(K)))  &
                           )
         ENDDO

         MaxAtmsNode=1+MaxAtmsNode
         MaxBlksNode=2.8D0*(1+MaxAtmsNode*BWEstim)
         MaxNon0Node=2.8D0*(1+MaxBaseNode*(DBLE(MaxBaseNode)*DBLE(BWEstim) &
                                   /DBLE(MaxAtmsNode)))
!-----------------------------------------------------------
!        To HDF with the indecies/limits
!
         CALL Put(MaxAtmsNode,'maxatmsnode',Tag_O=IntToChar(I))
         CALL Put(MaxBlksNode,'maxblksnode',Tag_O=IntToChar(I))
         CALL Put(MaxNon0Node,'maxnon0node',Tag_O=IntToChar(I))
         CALL Put(NPrc,'chknprc')
!---------------------------------------------------------
!        Tidy up
!
         CALL Delete(Beg)
         CALL Delete(End)
         CALL Delete(BSiz)
         CALL Delete(OffS)
#endif
!------------------------------------------------------------
!        Debug statements
!
         IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN 
            CALL OpenASCII(OutFile,Out)
            CALL PrintProtectL(Out)
            Mssg='Accuracy level '//TRIM(IntToChar(Acc))                            & 
                      //' for the '//TRIM(Base%BName)//' basis set:'//Rtrn  &
             //'  TrixNeglect = '//TRIM(DblToShrtChar(TrixNeglect(Acc)))            &
             //', CubeNeglect = '//TRIM(DblToShrtChar(CubeNeglect(Acc)))//','//Rtrn &
             //'  TwoENeglect = '//TRIM(DblToShrtChar(TwoENeglect(Acc)))            &
             //', DistNeglect = '//TRIM(DblToShrtChar(DistNeglect(Acc)))
            WRITE(Out,*)Mssg
            CALL PrintProtectR(Out)
            CLOSE(UNIT=Out,STATUS='KEEP')
         ENDIF
         IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
            CALL OpenASCII(OutFile,Out)
            CALL PrintProtectL(Out)
            Mssg='MaxAtms = '//TRIM(IntToChar(MaxAtms))   & 
             //', MaxBlks = '//TRIM(IntToChar(MaxBlks))   & 
             //', MaxNon0 = '//TRIM(IntToChar(MaxNon0))   &
             //', BndWdth = '//TRIM(IntToChar(BWEstim))
            WRITE(Out,*)Mssg
#ifdef PARALLEL
            Mssg='MaxAtmsNode = '//TRIM(IntToChar(MaxAtmsNode))//', ' & 
              //' MaxBlksNode = '//TRIM(IntToChar(MaxBlksNode))//', ' & 
              //' MaxNon0Node = '//TRIM(IntToChar(MaxNon0Node))
            WRITE(Out,*)Mssg
#endif
            CALL PrintProtectR(Out)
            CLOSE(UNIT=Out,STATUS='KEEP')
         ENDIF
      END SUBROUTINE SetLimits
#ifdef PARALLEL
!============================================================================
!     Brain-dead domain decomposition for basis set Ctrl%ISet.  
!============================================================================      
      SUBROUTINE Decomp(Ctrl,OverRideNProc_O)
         TYPE(SCFControls)                   :: Ctrl
         INTEGER, OPTIONAL                   :: OverRideNProc_O
         TYPE(INT_VECT)                      :: OffSt
         TYPE(CHR_VECT)                      :: Chr,CBeg,CEnd
         INTEGER                             :: NAtsAv,I,K,IPrc,ISet
         INTEGER, PARAMETER                  :: Mns=1,Pls=2
         INTEGER, DIMENSION(Mns:Pls)         :: Beg3,End2
         INTEGER, DIMENSION(Mns:Pls,Mns:Pls) :: End3,Dv
         CHARACTER(LEN=DEFAULT_CHR_LEN)      :: Mssg
         ISet=Ctrl%Current(2)
!---------------------------------------------------------------
!        Parse MPI_INVOKE for NPrc or procs
!
         CALL LineToChars(MPI_FLAGS,Chr)
         NPrc=FAIL
         DO I=1,SIZE(Chr%C)
            IF(TRIM(Chr%C(I))=='-np'.OR. &
               TRIM(Chr%C(I))=='-procs'  )THEN
               NPrc=CharToInt(Chr%C(I+1))
               IF(PRESENT(OverRideNProc_O))THEN
                  NPrc=OverRideNProc_O
                  Chr%C(I+1)=TRIM(IntToChar(NPrc))
                  MPI_FLAGS=' '
                  DO K=1,SIZE(Chr%C)
                     MPI_FLAGS=ADJUSTL(TRIM(MPI_FLAGS))//Blnk//TRIM(Chr%C(K))
                  ENDDO
               ENDIF                   
               EXIT
            ENDIF
         ENDDO
         IF(NPrc==FAIL)CALL Halt(' Failed to find NPrc in' &
                       //' MPI_FLAGS = <'//TRIM(MPI_FLAGS)//'>')
!----------------------------------------------------------------
!        Allocate domain limits
!
         IF(.NOT.AllocQ(Beg%Alloc))CALL New(Beg,NPrc-1,0)
         IF(.NOT.AllocQ(End%Alloc))CALL New(End,NPrc-1,0)
!------------------------------------------------------------------------
!        Look ahead 1 node algorithm for decomposition
!        If NPrc==2^p, should use bisection instead. 
!        -- See ORB.F90 for work in progress...
!
         NBasF=SUM(BSiz%I(1:NAtoms)) 
         Beg%I(0)=1
         End%I(NPrc-1)=NAtoms
         DO IPrc=0,NPrc-2
!---------------------------------------------------------------------
!           Compute running average to section
!
            NAtsAv=DBLE(SUM(BSiz%I(Beg%I(IPrc):NAtoms))) &
                  /DBLE(NPrc-IPrc)
!--------------------------------------------------------------------
!           Forcast Beg and End for (de/inc)rements of (+/- 1) 
!
            DO K=Beg%I(IPrc),NAtoms
               IF(SUM(BSiz%I(Beg%I(IPrc):K))>=NAtsAv)THEN
                  End2(Pls)=K
                  EXIT            
               ENDIF
            ENDDO
            End2(Mns)=End2(Pls)-1
            Beg3(Pls)=End2(Pls)+1
            Beg3(Mns)=End2(Mns)+1
            End3(Pls,Pls)=NAtoms            
            DO K=Beg3(Pls),NAtoms   
               IF(SUM(BSiz%I(Beg3(Pls):K))>=NAtsAv)THEN
                  End3(Pls,Pls)=K
                  EXIT            
               ENDIF 
            ENDDO
            End3(Pls,Mns)=End3(Pls,Pls)-1
            DO K=Beg3(Mns),NAtoms    
               IF(SUM(BSiz%I(Beg3(Mns):K))>=NAtsAv)THEN
                  End3(Mns,Pls)=K
                  EXIT            
               ENDIF
            ENDDO
            End3(Mns,Mns)=End3(Mns,Pls)-1
!-------------------------------------------------------------------------------
!           These are deviations from the running average for each choice
!           of a place to section and the possilbe sectioning in the next
!           iteration
!           
            Dv(Pls,Pls)= &
              ABS(NAtsAv-SUM(Bsiz%I(Beg%I(IPrc):End2(Pls)))) &
             +ABS(NAtsAv-SUM(Bsiz%I(Beg3(Pls):End3(Pls,Pls)))) 
            Dv(Pls,Mns)= &
              ABS(NAtsAv-SUM(Bsiz%I(Beg%I(IPrc):End2(Pls))))   &
             +ABS(NAtsAv-SUM(Bsiz%I(Beg3(Pls):End3(Pls,Mns)))) 
            Dv(Mns,Pls)= &
              ABS(NAtsAv-SUM(Bsiz%I(Beg%I(IPrc):End2(Mns))))   &
             +ABS(NAtsAv-SUM(Bsiz%I(Beg3(Mns):End3(Mns,Pls)))) 
            Dv(Mns,Mns)= &
              ABS(NAtsAv-SUM(Bsiz%I(Beg%I(IPrc):End2(Mns))))   &
             +ABS(NAtsAv-SUM(Bsiz%I(Beg3(Mns):End3(Mns,Mns)))) 
!write(*,*)' 0th ', &
!              ABS(NAtsAv-SUM(Bsiz%I(Beg%I(IPrc):End2(Pls)))),   &
!              ABS(NAtsAv-SUM(Bsiz%I(Beg%I(IPrc):End2(Mns))))   
!write(*,*)' 1st ', &
!             ABS(NAtsAv-SUM(Bsiz%I(Beg3(Pls):End3(Pls,Pls)))) ,&
!             ABS(NAtsAv-SUM(Bsiz%I(Beg3(Pls):End3(Pls,Mns)))) ,&
!             ABS(NAtsAv-SUM(Bsiz%I(Beg3(Mns):End3(Mns,Pls)))) ,&
!             ABS(NAtsAv-SUM(Bsiz%I(Beg3(Mns):End3(Mns,Mns)))) 
!-------------------------------------------------------------------------------
!           Pick the best section based on minimizing the I and I+1 deviation
!           from the average
!          
            IF(MIN(Dv(Pls,Pls),Dv(Pls,Mns))< &
               MIN(Dv(Mns,Pls),Dv(Mns,Mns)))THEN
               End%I(IPrc)=End2(Pls)
               Beg%I(IPrc+1)=End%I(IPrc)+1
!               WRITE(*,12)Dv(Pls,Pls),Dv(Pls,Mns),Dv(Mns,Pls),Dv(Mns,Mns), &
!                          IPrc,Beg%I(IPrc),End%I(IPrc),                    &
!                          SUM(Bsiz%I(Beg%I(IPrc):End%I(IPrc)))
            ELSE
               End%I(IPrc)=MAX(Beg%I(IPrc),End2(Mns))
               Beg%I(IPrc+1)=End%I(IPrc)+1
!               WRITE(*,13)Dv(Pls,Pls),Dv(Pls,Mns),Dv(Mns,Pls),Dv(Mns,Mns), &
!                          IPrc,Beg%I(IPrc),End%I(IPrc) ,                   &
!                          SUM(Bsiz%I(Beg%I(IPrc):End%I(IPrc)))
            ENDIF
         ENDDO
!12 format('I  ++ = ',I3,', +- = ',I3,', -+ = ',I3,', -- = ',I3, &
!          ' IPrc = ',I3,' B = ',I3,' E = ',I3,' Sum = ',I4)
!13 format('II ++ = ',I3,', +- = ',I3,', -+ = ',I3,', -- = ',I3, &
!          ' IPrc = ',I3,' B = ',I3,' E = ',I3,'Sum = ',I4)
!---------------------------------------------------------        
!        Calculate DBCSR matrix RowPt -> GRwPt off-sets
!
         CALL New(OffSt,NPrc-1,0)
         DO I=0,NPrc-1
            OffSt%I(I)=End%I(I)-Beg%I(I)+1
         ENDDO
         DO I=1,NPrc-1
            OffSt%I(I)=OffSt%I(I)+OffSt%I(I-1)
         ENDDO
         DO I=NPrc-1,1,-1
            OffSt%I(I)=OffSt%I(I-1)
         ENDDO
         OffSt%I(0)=0
!-------------------------------------------------------
!        Put the domain boundaries to disk
!
         CALL Put(Beg,'beg',Tag_O=IntToChar(ISet))
         CALL Put(End,'end',Tag_O=IntToChar(ISet))
         CALL Put(OffSt,'dbcsroffsets',Tag_O=IntToChar(ISet))
!-------------------------------------------------------
!        Debug if asked
!
         IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
            CALL OpenASCII(OutFile,Out)           
            CALL PrintProtectL(Out)
            Mssg='MPI Invokation : '//TRIM(MPI_INVOKE)
            WRITE(Out,*)Mssg
            CALL New(CBeg,NPrc-1,0)
            CALL New(CEnd,NPrc-1,0)
            DO I=0,NPrc-1
               CBeg%C(I)=IntToChar(Beg%I(I))
               CEnd%C(I)=IntToChar(End%I(I))
            ENDDO
            WRITE(Out,*)'Atomic partitioning = ',                                 &
                 ('['//TRIM(CBeg%C(K))//'-'//TRIM(CEnd%C(K))//'], ',K=0,NPrc-2),  &
                  '['//TRIM(CBeg%C(NPrc-1))//'-'//TRIM(CEnd%C(NPrc-1)),']'

            DO K=0,NPrc-1
               CBeg%C(K)=IntToChar(SUM(Bsiz%I(Beg%I(K):End%I(K))))
            ENDDO
            WRITE(Out,*)'Basis Functions per node = ',                &
               (TRIM(CBeg%C(K))//', ',K=0,NPrc-2),TRIM(CBeg%C(NPrc-1))
            CALL PrintProtectR(Out)
            CLOSE(UNIT=Out,STATUS='KEEP')
            CALL Delete(CBeg)
            CALL Delete(CEnd)
         ENDIF
!-------------------------------------------------------
!        Tidy up
!
         CALL Delete(Chr)
         CALL Delete(OffSt)
      END SUBROUTINE Decomp
#endif

      SUBROUTINE SetDist(Ctrl)
         TYPE(SCFControls) :: Ctrl
         TYPE(DBL_VECT)    :: DExpt
         TYPE(INT_VECT)    :: Lndex
         TYPE(DBL_RNK6)    :: PrmEst
         INTEGER           :: I,NExpt
         CALL CalcPDist(Base,NExpt,DExpt,Lndex)
         CALL Put(NExpt,'nexpt',Tag_O=CSet)
         CALL Put(DExpt,'dexpt',Tag_O=CSet)
         CALL Put(Lndex,'lndex',Tag_O=CSet)
         CALL Delete(Lndex)
         CALL Delete(Dexpt)
!         CALL BKEst(Base,PrmEst)
!         CALL Put(PrmEst,'best',Tag_O=CSet)
!         CALL Delete(PrmEst)
      END SUBROUTINE SetDist
!-----------------------------------------------------------------------------
!     Pre-compute some primitive distribution stuff
!
      SUBROUTINE CalcPDist(BS,NExpt,DExpt,Lndex)
         TYPE(BSET),     INTENT(INOUT):: BS
         TYPE(DBL_VECT), INTENT(OUT)  :: DExpt
         TYPE(INT_VECT), INTENT(OUT)  :: Lndex
         INTEGER,        INTENT(OUT)  :: NExpt
         TYPE(INT_VECT)               :: ITmp,IPnt
         INTEGER                      :: I,K,KA,KB,CFA,CFB,PFA,PFB,NPB
!---------------------------------------------------------------------------
!        First count the types of primitive distributions(K=1), then compute 
!        their primitive distribution exponents and associated values of 
!        angular symmetry (K=2)
!
         DO K=1,2
            IF(K==1)THEN
               NExpt=1
            ELSE
               CALL New(DExpt,NExpt)
               CALL New(Lndex,NExpt)
               NExpt=1
            ENDIF               
            DO KA=1,BS%NKind
               DO KB=1,KA-1
                  DO CFA=1,BS%NCFnc%I(KA)
                     DO CFB=1,BS%NCFnc%I(KB)
                        DO PFA=1,BS%NPFnc%I(CFA,KA)
                           DO PFB=1,BS%NPFnc%I(CFB,KB)
                              IF(K==2)THEN
                                 DExpt%D(NExpt)=BS%Expnt%D(PFA,CFA,KA)  &
                                               +BS%Expnt%D(PFB,CFB,KB)  
                                 Lndex%I(NExpt)=BS%ASymm%I(2,CFA,KA)    &
                                               +BS%ASymm%I(2,CFB,KB)
                              ENDIF
                              NExpt=NExpt+1          
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            DO KA=1,BS%NKind
               DO CFA=1,BS%NCFnc%I(KA)
                  DO CFB=1,CFA
                     DO PFA=1,BS%NPFnc%I(CFA,KA)
                        NPB=BS%NPFnc%I(CFB,KA)
                        IF(CFA.EQ.CFB)NPB=PFA
                        DO PFB=1,NPB
                            IF(K==2)THEN 
                               DExpt%D(NExpt)=BS%Expnt%D(PFA,CFA,KA)  &
                                             +BS%Expnt%D(PFB,CFB,KA)  
                               Lndex%I(NExpt)=BS%ASymm%I(2,CFA,KA)    &
                                             +BS%ASymm%I(2,CFB,KA)
                            ENDIF
                            NExpt=NExpt+1 
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
!-------------------------------------------------------
!        Add exponent for nuclear delta function
!
         DExpt%D(NExpt)=NuclearExpnt
         Lndex%I(NExpt)=0
!-------------------------------------------------------
!        Sort the exponents in assending order and carry
!        allong angular symetry indecies
!
         CALL New(IPnt,NExpt)
         CALL New(ITmp,NExpt)
         DO I=1,NExpt
            IPnt%I(I)=I
         ENDDO
         CALL Sort(DExpt,IPnt,NExpt,2)   
         DO I=1,NExpt
            ITmp%I(I)=Lndex%I(IPnt%I(I))
         ENDDO
         DO I=1,NExpt
            Lndex%I(I)=ITmp%I(I)
         ENDDO
!--------------------------------
!        Debug
!
!         IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
!            CALL PPrint(DExpt,'DExpt')
!            CALL PPrint(Lndex,'Lndex')
!         ENDIF
!--------------------------------
!        Tidy up
!
         CALL Delete(ITmp)
         CALL Delete(IPnt)
   END SUBROUTINE CalcPDist

      SUBROUTINE BKEst(BS,PrmEst)
      TYPE(BSET),     INTENT(IN)  :: BS
      TYPE(DBL_RNK6), INTENT(OUT) :: PrmEst 
      TYPE(DBL_RNK4)              :: MD,R 
      TYPE(DBL_VECT)              :: AuxR
      REAL(DOUBLE)                :: ZetaA,ZetaB,EtaAB,Omega, &
                                     Upq,Ca,Cb,CC2,Int
      INTEGER                     :: KA,KB,CFA,CFB,PFA,PFB,           &
                                     StartLA,StartLB,StopLA,StopLB,   &
                                     StrideA,StrideB,StrideAB,        &
                                     MaxL,LTot,MaxLA,MaxLB,LMNA,LMNB, &
                                     LA,MA,NA,LB,MB,NB,LAB,MAB,NAB,   &
                                     LCD,MCD,NCD
!-----------------------------------------------------------
!     Allocations
!
      MaxL=4*BS%NASym
      CALL New(PrmEst,(/BS%NPrim,BS%NCtrt,BS%NKind, &
                        BS%NPrim,BS%NCtrt,BS%NKind/))
      CALL New(MD,(/3,BS%NASym,BS%NASym,MaxL/),(/1,0,0,0/))
      CALL New(R,(/MaxL,MaxL,MaxL,MaxL/),(/0,0,0,0/))
      CALL New(AuxR,MaxL,0)
!-----------------------------------------------------------
!     Loops over contracted functions
!
      DO KA=1,BS%NKind
      DO CFA=1,BS%NCFnc%I(KA)
         MaxLA=BS%ASymm%I(2,CFA,KA)
         StartLA=BS%LStrt%I(CFA,KA)
         StopLA =BS%LStop%I(CFA,KA)
         StrideA=StopLA-StartLA+1
         DO KB=1,BS%NKind
         DO CFB=1,BS%NCFnc%I(KB)
            MaxLB=BS%ASymM%I(2,CFB,KB)
            StartLB=BS%LStrt%I(CFB,KB)
            StopLB =BS%LStop%I(CFB,KB)
            StrideB=StopLB-StartLB+1
            StrideAB=StrideA*StrideB
            LTot=2*(MaxLA+MaxLB)
!-----------------------------------------------------------
!           Loop over primitive functions
!
            DO PFA=1,BS%NPFnc%I(CFA,KA)
            DO PFB=1,BS%NPFnc%I(CFB,KB)
               ZetaA=BS%Expnt%D(PFA,CFA,KA)
               ZetaB=BS%Expnt%D(PFB,CFB,KB)
               EtaAB=ZetaA+ZetaB
               Omega=EtaAB*Half
               Upq=Sqrt2Pi5x2*EtaAB**(-FiveHalves)
               CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,EtaAB,MD%D, &
                           Zero,Zero,Zero,Zero,Zero,Zero) 
               CALL AuxInts(MaxL,LTot,AuxR%D,Omega,Zero)
               CALL MD3TRR(MaxL,LTot,R%D,AuxR%D,Upq,Zero,Zero,Zero)
!----------------------------------------------------------------------------
!              (ab|ab)=(Ca*Cb)^2*Upq*Sum_{LMN} Sum_{LpMpNp}( 
!                       (-1)^(Lp+Mp+Np) E^{ab}_{LMN} E^{ab}_{LpMpNp} R_{L+Lp,M+Mp,N+Np} )
!
!              PrmEst=Max(PrmEst,SQRT|(ab|ab)|)
!
               PrmEst%D(PFB,CFB,KB,PFA,CFA,KA)=Zero
               DO LMNA=StartLA,StopLA
                  LA=BS%LxDex%I(LMNA)
                  MA=BS%LyDex%I(LMNA)
                  NA=BS%LzDex%I(LMNA)
                  CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                  DO LMNB=StartLB,StopLB
                     LB=BS%LxDex%I(LMNB)
                     MB=BS%LyDex%I(LMNB)
                     NB=BS%LzDex%I(LMNB)
                     CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                     CC2=(CA*CB)**2
                     Int=Zero 
                     DO LAB=0,LA+LB
                     DO MAB=0,MA+MB
                     DO NAB=0,NA+NB
                        DO LCD=0,LA+LB 
                        DO MCD=0,MA+MB
                        DO NCD=0,NA+NB
                           Int=Int                   &
                              +MD%D(1,LA,LB,LAB)     &
                              *MD%D(2,MA,MB,MAB)     &
                              *MD%D(3,NA,NB,NAB)     &
                              *MD%D(1,LA,LB,LCD)     &
                              *MD%D(2,MA,MB,MCD)     &
                              *MD%D(3,NA,NB,NCD)     &
                              *(-One)**(LCD+MCD+NCD) &
                              *R%D(LAB+LCD,MAB+MCD,NAB+NCD,0)
                        ENDDO
                        ENDDO
                        ENDDO
                     ENDDO
                     ENDDO
                     ENDDO
                     PrmEst%D(PFB,CFB,KB,PFA,CFA,KA)=MAX( &
                     PrmEst%D(PFB,CFB,KB,PFA,CFA,KA),CC2*ABS(Int))
                  ENDDO
               ENDDO
!               PrmEst%D(PFB,CFB,KB,PFA,CFA,KA)=SQRT(Upq* &
!               PrmEst%D(PFB,CFB,KB,PFA,CFA,KA))
               PrmEst%D(PFB,CFB,KB,PFA,CFA,KA)=BIG_DBL
!               WRITE(*,*)' EtaAB = ',EtaAB,' PrmEst = ',PrmEst%D(PFB,CFB,KB,PFA,CFA,KA)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
!-------------------------------------------------------------
!     Tidy up
!
      CALL Delete(MD)
      CALL Delete(R)
      CALL Delete(AuxR)
      END SUBROUTINE BKEst
END MODULE
