  SUBROUTINE RangeOfDensity(D,NameBuf,BfnInd,DB,BS,GM)
    USE DerivedTypes
    USE GlobalScalars
    USE Macros
    USE ONXParameters
    IMPLICIT NONE
#ifdef PARALLEL_ONX
    TYPE(DBCSR),INTENT(IN)        :: D
#else
    TYPE(BCSR),INTENT(IN)         :: D
#endif
    TYPE(INT_VECT),INTENT(INOUT)  :: NameBuf
    TYPE(INT_RNK2),INTENT(INOUT)  :: BfnInd
    TYPE(BSET),INTENT(IN)         :: BS
    TYPE(CRDS),INTENT(IN)         :: GM
    TYPE(DBuf)                    :: DB
    INTEGER                       :: AtA,ShellA,KA,CFA
    INTEGER                       :: AtB
    INTEGER                       :: ri,ci

#ifdef PARALLEL_ONX
    write(*,*)'Natoms',Natoms,' MyID:',MyID,NPrc 
    write(*,*)'Beg%I(MyID)',Beg%I(MyID)
    write(*,*)'End%I(MyID)',End%I(MyID)
#endif

    NameBuf%I(:)=0
#ifdef PARALLEL_ONX
    DO AtA=Beg%I(MyID),End%I(MyID)
       ri=AtA-Beg%I(MyID)+1
!       write(*,*) 'ri',ri
!       write(*,*) 'AtA',AtA
       NameBuf%I(AtA)=1
       DO ci=D%RowPt%I(ri),D%RowPt%I(ri+1)-1
          AtB=D%ColPt%I(ci)
!       write(*,*) 'bound D%ColPt%I',lbound(D%ColPt%I(:)),ubound(D%ColPt%I(:))
!       write(*,*) 'AtB',AtB
          NameBuf%I(AtB)=1
       END DO
    END DO
#else   !vw this is strange, why only this 2 loops? !
    DO AtA=1,NAtoms
       NameBuf%I(AtA)=1
       !write(*,*) 'AtA',AtA
       DO ci=D%RowPt%I(AtA),D%RowPt%I(AtA+1)-1
          AtB=D%ColPt%I(ci)
       !write(*,*) 'bound D%ColPt%I',lbound(D%ColPt%I(:)),ubound(D%ColPt%I(:))
       !write(*,*) 'AtB',AtB
          NameBuf%I(AtB)=1
       END DO
    END DO
!vw---->
!    stop 9999
#endif
!vw<----
    ShellA=0
    DO AtA=1,NAtoms
       KA=GM%AtTyp%I(AtA)
       DO CFA=1,BS%NCFnc%I(KA)
          ShellA=ShellA+1
          BfnInd%I(AtA,CFA)=ShellA
       END DO
    END DO
    DB%NShells=ShellA
!#endif

    !write(*,*) MyId,'NameBuf%Alloc:',NameBuf%Alloc
    !write(*,*) MyId,'NameBuf%I:',NameBuf%I
    !write(*,*) MyId,'Size(NameBuf):',Size(NameBuf%I(:))
    !write(*,*) 'NameBuf%I',NameBuf%I
  END SUBROUTINE RangeOfDensity

