SUBROUTINE Scatter(N,NA,NB,IndexA,SB,SubInd,DB,KB,K,ivalue,namebuf)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  USE ONXMemory
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: N,NA,NB,IndexA
  TYPE(DSL),INTENT(IN)       :: SB
  TYPE(INT_RNK2),INTENT(IN)  :: SubInd
  TYPE(DBuf),INTENT(IN)      :: DB
  REAL(DOUBLE),INTENT(IN)    :: KB(N,NA,NB)
#ifdef PARALLEL_ONX
  TYPE(DBCSR),INTENT(INOUT)  :: K

  integer :: ria,rib

#else
  TYPE(BCSR) ,INTENT(INOUT)  :: K
#endif
  INTEGER                    :: I,I0,Ind,IndexB,Ioff
  INTEGER                    :: AtA,NBFA,RS
  INTEGER                    :: AtB,NBFB,CS


!--test--test--test--test--test--test--test--test--
  type(INT_VECT),intent(in) :: namebuf
  integer, dimension(0:20),intent(inout) :: ivalue
  integer :: j,l
  real*8 :: dbl
!--test--test--test--test--test--test--test--test--

  !if(MyID==2) Write(*,*) KB

!  write(*,*) 'I am in Ke 11.4.1',MyID,NPrc !vw


!!$  do i=1,N
!!$     do j=1,NA
!!$        do l=1,NB
!!$           dble=KB(i,j,l)
!!$           if(abs(KB(i,j,l)).lt.10D-10) dble=0.0D+00
!!$           write(*,'(A,4I3,E24.14)') 'MyID',MyID,i,j,l,dble
!!$        enddo
!!$     enddo
!!$  enddo


  DO I=1,N
     I0     = SB%SLDis%I(I)-3 !some magic number!
     IndexB = INT(ABS(DB%DisBuf%D(I0)))
     AtA    = SubInd%I(1,IndexA)
     NBFA   = SubInd%I(2,IndexA)
     RS     = SubInd%I(3,IndexA)
     AtB    = SubInd%I(1,IndexB)
     NBFB   = SubInd%I(2,IndexB)
     CS     = SubInd%I(3,IndexB)
     !write(*,'(A,10I3)') 'ID',IndexA,IndexB,AtA,AtB,NBFA,NBFB,RS,CS,MyID,NPrc
     !    write(*,*) 'I am in Ke 11.4.10',MyID,NPrc !vw
     !write(*,'(A,10I3)') 'ID',IndexA,IndexB,AtA,AtB,RS,CS,MyID

     !CAUTION CAUTION CAUTION CAUTION CAUTION
     !----> I may get problems here !vw

#ifdef PARALLEL_ONX
     !
     !write(*,'5(A,I3)') 'MyID',MyID,' AtA',AtA,' AtB',AtB
     !CALL GetAdrB(ria,rib,Ind,K,0)
     !write(*,*) 'namebuf%I',size(namebuf%I),namebuf%I(:)
     
     CALL GetAdrB(Binary_Search(namebuf%I(:),AtA,NRows),AtB,Ind,K,0)
     !CALL GetAdrB(AtA,AtB,Ind,K,0)
#else
     CALL GetAdrB(AtA,AtB,Ind,K,0)
#endif
     !<----
     !CAUTION CAUTION CAUTION CAUTION CAUTION


     !write(*,*) 'Ind',Ind

!!$     select case(Ind)!case(IndexB)
!!$     case( 0);ivalue( 0)=ivalue( 0)+1
!!$     case( 1);ivalue( 1)=ivalue( 1)+1
!!$     case( 2);ivalue( 2)=ivalue( 2)+1
!!$     case( 3);ivalue( 3)=ivalue( 3)+1
!!$     case( 4);ivalue( 4)=ivalue( 4)+1
!!$     case( 5);ivalue( 5)=ivalue( 5)+1
!!$     case( 6);ivalue( 6)=ivalue( 6)+1
!!$     case( 7);ivalue( 7)=ivalue( 7)+1
!!$     case( 8);ivalue( 8)=ivalue( 8)+1
!!$     case( 9);ivalue( 9)=ivalue( 9)+1
!!$     case(10);ivalue(10)=ivalue(10)+1
!!$     case(11);ivalue(11)=ivalue(11)+1
!!$     case(12);ivalue(12)=ivalue(12)+1
!!$     case(13);ivalue(13)=ivalue(13)+1
!!$     case(14);ivalue(14)=ivalue(14)+1
!!$     case(15);ivalue(15)=ivalue(15)+1
!!$     case(16);ivalue(16)=ivalue(16)+1
!!$     case(17);ivalue(17)=ivalue(17)+1
!!$     case(18);ivalue(18)=ivalue(18)+1
!!$     case(19);ivalue(19)=ivalue(19)+1
!!$     case(20);ivalue(20)=ivalue(20)+1
!!$     case default; write(*,*) 'BUG';stop 999999
!!$     end select
!!$     cycle

     !    write(*,*) 'I am in Ke 11.4.11',MyID,NPrc !vw
     Ioff=K%BlkPt%I(Ind)
     !    write(*,*) 'I am in Ke 11.4.12',MyID,NPrc !vw
     IF (Ind > 0) CALL PutSubBlk(I,N,NBFA,NBFB,NA,NB,RS,CS, &
          &                      K%MTrix%D(Ioff),KB)
     !    write(*,*) 'I am in Ke 11.4.13',MyID,NPrc !vw

!!$     if(MyID==0) write(*,*) 'AtA',AtA,' AtB',AtB
!!$     do l=0,NBFA*NBFB-1
!!$        dbl=K%MTrix%D(Ioff+l)
!!$        if(abs(K%MTrix%D(Ioff+l)).lt.10D-10) dbl=0.0D+00
!!$        if(MyID==0) write(*,'(A,I3,E24.14)') 'MyID',MyID,dbl
!!$     enddo

  END DO

contains

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
END SUBROUTINE Scatter
