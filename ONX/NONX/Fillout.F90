SUBROUTINE Fillout_BCSR(BS,GM,A)
  USE DerivedTypes
  USE GlobalScalars
  IMPLICIT NONE
  TYPE(BSET),INTENT(IN)    :: BS
  TYPE(CRDS),INTENT(IN)    :: GM
  TYPE(BCSR),INTENT(INOUT) :: A
  INTEGER                  :: AtA,KA,NBFA
  INTEGER                  :: AtB,KB,NBFB
  INTEGER                  :: iPnt1,iPnt2,ci,Ind

#ifdef PARALLEL_ONX     !useful if reduce to root and after transpose !vw
  IF (MyID==ROOT) THEN
#endif
     DO AtA=1,NAtoms
        KA=GM%AtTyp%I(AtA)
        NBFA=BS%BfKnd%I(KA)
        DO ci=A%RowPt%I(AtA),A%RowPt%I(AtA+1)-1
           AtB=A%ColPt%I(ci)
           KB=GM%AtTyp%I(AtB)
           NBFB=BS%BfKnd%I(KB)
           IF (AtA.GE.AtB) THEN
              CALL GetAdrB(AtB,AtA,Ind,A,0)
              iPnt1=A%BlkPt%I(ci)
              iPnt2=A%BlkPt%I(Ind)
              IF (iPnt1.EQ.iPnt2) THEN
                 CALL XPose1C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(ipnt2))
              ELSE
                 CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(ipnt2))
              END IF
           END IF
        END DO
     END DO
#ifdef PARALLEL_ONX
  END IF
#endif
END SUBROUTINE Fillout_BCSR

#ifdef PARALLEL_ONX
SUBROUTINE Fillout_DBCSR(BS,GM,A,NameBuf)
  USE DerivedTypes
  USE GlobalScalars
  USE ONXParameters
  IMPLICIT NONE
  TYPE(BSET),INTENT(IN)     :: BS
  TYPE(CRDS),INTENT(IN)     :: GM
  TYPE(DBCSR),INTENT(INOUT) :: A
  TYPE(INT_VECT),INTENT(IN) :: NameBuf
  INTEGER                   :: AtA,KA,NBFA
  INTEGER                   :: AtB,KB,NBFB
  INTEGER                   :: iPnt1,iPnt2,ci,Ind,ri,rib

  !vw
  write(*,*) MyId,' Nrows',Nrows
  write(*,*) MyId, 'HERE: A%RowPt%I',A%RowPt%I(1),A%RowPt%I(nrows+1)

  !write(*,*) MyId,'NameBuf%Alloc:',NameBuf%Alloc
  !write(*,*) MyId,'NameBuf%I:',NameBuf%I
  !write(*,*) MyId,'Size(NameBuf):',Size(NameBuf%I(:))
  !write(*,*) MyId,'A%RowPt%Alloc',A%RowPt%Alloc
  !write(*,*) MyId,'A%RowPt%I',A%RowPt%I
  write(*,*) MyId,'Size(A%RowPt%I):',Size(A%RowPt%I)
  write(*,*) MyId,'Size(A%ColPt%I):',Size(A%ColPt%I)
  !vw

  DO ri=1,NRows
     AtA=NameBuf%I(ri)
     KA=GM%AtTyp%I(AtA)
     NBFA=BS%BfKnd%I(KA)
     !DO ci=A%RowPt%I(AtA),A%RowPt%I(AtA+1)-1 !<--- problem here ! This is wrong
     DO ci=A%RowPt%I(ri),A%RowPt%I(ri+1)-1 !<--- I am not sure ! Now should be ok!
       
        !write(*,*) 'A%ColPt%I(ata)',A%ColPt%I(ata),' A%ColPt%I(ata+1)',A%ColPt%I(ata+1)-1,ci
        !write(*,*) 'A%ColPt%I(ci)',A%ColPt%I(ci),' A%ColPt%I(ci+1)',A%ColPt%I(ci+1)-1,ci
       
        call checkI1(A%ColPt%I(:),ci,lbound(A%ColPt%I(:)),ubound(A%ColPt%I(:)))

        AtB=A%ColPt%I(ci)
        KB=GM%AtTyp%I(AtB)
        NBFB=BS%BfKnd%I(KB)
        IF (AtA.GE.AtB) THEN
!vw        CALL GetAdrB(AtB,AtA,Ind,0,A%ColPt%I,A%RowPt%I)

! Find the row which correspond to AtB
           rib=Binary_Search(NameBuf%I(:),AtB,NRows)

!vw-->
           CALL GetAdrB(rib,ri,Ind,A,0) !<--- I am not sure ! Now should be ok!
           !CALL GetAdrB(AtB,AtA,Ind,A,0) !<--- problem here ! This is wrong
!vw<--

           call checkI1(A%BlkPt%I(:),ci,lbound(A%BlkPt%I(:)),ubound(A%BlkPt%I(:)))
           call checkI1(A%BlkPt%I(:),ind,lbound(A%BlkPt%I(:)),ubound(A%BlkPt%I(:)))

           iPnt1=A%BlkPt%I(ci)
           iPnt2=A%BlkPt%I(Ind)

           if(myid==2) write(*,'7(A,I3)') 'MyID',MyID,' ri',ri,' AtA',AtA,' rib',rib,' AtB',AtB,' iPnt1',iPnt1,' iPnt2',iPnt2

           call checkD1(A%MTrix%D(:),iPnt1,lbound(A%MTrix%D(:)),ubound(A%MTrix%D(:)))
           call checkD1(A%MTrix%D(:),iPnt2,lbound(A%MTrix%D(:)),ubound(A%MTrix%D(:)))

           call checkD1(A%MTrix%D(:),iPnt1+NBFA*NBFB-1,lbound(A%MTrix%D(:)),ubound(A%MTrix%D(:)))
           call checkD1(A%MTrix%D(:),iPnt2+NBFA*NBFB-1,lbound(A%MTrix%D(:)),ubound(A%MTrix%D(:)))


           if(myid==2)call printD(A%MTrix%D(iPnt1:iPnt1+NBFA*NBFB-1),'old')

           IF (iPnt1.EQ.iPnt2) THEN
              CALL XPose1C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(ipnt2))
           ELSE
              CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(ipnt2))
           END IF

           if(myid==2)call printD(A%MTrix%D(iPnt2:iPnt2+NBFA*NBFB-1),'new')

        END IF
     END DO
  END DO
  
contains

  subroutine printD(D,name)
    real(double), dimension(:), intent(in) :: d
    integer :: j
    character(len=*), intent(in) :: name
    write(*,*) MyID
    write(*,*) trim(name)
    do j=1,size(D)
       write(*,*) D(j)
    enddo
  end subroutine PrintD


  FUNCTION Binary_Search(IVec,IVal,NDim) RESULT(Idx)
    IMPLICIT NONE
    INTEGER :: Idx,IMin,IMax
    INTEGER, INTENT(In) :: IVal,NDim
    INTEGER, DIMENSION(:), INTENT(In) :: IVec
    LOGICAL :: Found
    IF(IVal.lt.IVec(1).or.IVal.gt.IVec(NDim).or.NDim.gt.size(IVec)) THEN
       Idx=-1
       RETURN
    ENDIF
    IMin=1    
    IMax=NDim
    Idx=NDim/2
    Found=.FALSE.
    DO
       IF(IVec(Idx).eq.IVal) THEN
          Found=.TRUE.
          EXIT
       ELSEIF(IVec(Idx).lt.IVal) THEN
          IMin=Idx
          Idx=IMax-(IMax-IMin)/2
          IF(Idx.lt.1) EXIT
       ELSE
          IMax=Idx
          Idx=(IMax-IMin)/2+IMin
          IF(Idx.gt.NDim) EXIT
       ENDIF
    ENDDO
    IF(.NOT.Found) Idx=-1
  END FUNCTION Binary_Search

  subroutine checkI1(I,idx,lb,ub)
    integer,intent(in) :: idx,lb,ub
    integer,dimension(lb:ub),intent(in) :: I
    integer,save :: inbr=1

    if(idx.gt.ub.or.idx.lt.lb) then
       write(*,*) 'idx does not match with array size'
       write(*,*) 'lbound=',lb,' ubound=',ub,' size=',size(I),' idx',idx
       stop
    endif
    inbr=inbr+1
  end subroutine checkI1


  subroutine checkD1(I,idx,lb,ub)
    integer,intent(in) :: idx,lb,ub
    real(double),dimension(lb:ub),intent(in) :: I
    integer,save :: inbr=1

    if(idx.gt.ub.or.idx.lt.lb) then
       write(*,*) 'idx does not match with array size'
       write(*,*) 'lbound=',lb,' ubound=',ub,' size=',size(I),' idx',idx
       stop
    endif
    inbr=inbr+1
  end subroutine checkD1


END SUBROUTINE Fillout_DBCSR
#endif


