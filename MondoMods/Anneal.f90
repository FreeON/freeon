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
! Author: C.K. Gan

! a module to solve the travelling salesman problem with the begin
! and end cities as far as possible.
! Opt_Trav_Band contains the subroutine Anneal to obtain a band matrix S
! by minimising the cost function
! F = (Alpha/DMax)*(Sum_i D_{i,i+1}) + &
!     Sum_{i>j} exp( Gamma* [abs(i-j)*DMax/(2*NCity*Dij)]**NLFac] )
! where i and j are the row and column indices, and Dij is the physical
! distance between atoms visited on the i-th and j-th time respectively.

! abbreviation list ---
! TravO: Traverse Order ! NCity: Number of Cities ! MD: Main Diagonal
! OD: Off Diagonal ! SD: Store Distances

MODULE Opt_Trav_Band
  USE MemMan
  USE PrettyPrint

  IMPLICIT NONE

  TYPE(TIME),SAVE::AnnealTime
  REAL(DOUBLE),PRIVATE::DMax,Beta,RCut,Gamma,Alpha,NLFac,Eta
  INTEGER,PRIVATE::NCity,FNumOpt,N_NEIGHBOR

#define NLFacEq1
#define SD ! #undef SD or #define SD
#if defined(SD)
  TYPE(DBL_VECT)::DistArr
  INTEGER::NumPair
#endif

CONTAINS

!---------------------------------------------------------------------------

  SUBROUTINE Anneal(RCoor,NCityV,TravO,GammaV,NLFacV,N_NeighborV,AlphaV)
    TYPE(DBL_RNK2)::RCoor
    TYPE(DBL_VECT)::XX,YY,ZZ,NeighborDis
    INTEGER::IndexInt,NCityV,LeftCity,RightCity,I,J,K,NOver,&
      NLimit,NSucc,N_NeighborV
    TYPE(INT_VECT)::TravO,ChosenCity,NewTravO,FCity,MCity
    REAL(DOUBLE) ::NLFacV,GammaV,AlphaV,IncCost,TotMDCost1,TotMDCost,&
      IncMDCost,T,MaxIncCost,TotODCost1,TotODCost,Dec,IncODCost,TotCost
    REAL(DOUBLE),EXTERNAL::Random
    CHARACTER(LEN=*),PARAMETER::RunTravOFile='RunTravO',TravOExt='.dat'
    REAL(DOUBLE),PARAMETER::MaxCostFrac=0.50D0,TFacR=0.90D0
    INTEGER,PARAMETER::RunTravOU=49,Temp_Step=150
    LOGICAL::Ans,Okay

    CALL OpenASCII(OutFile,Out)
#ifdef SD
    WRITE(Out,*) 'SD is defined!'
#else
    WRITE(Out,*) 'SD is NOT defined!'
#endif
#ifdef NLFacEq1
    WRITE(Out,*) 'SNLFacEq1 is defined!'
#else
    WRITE(Out,*) 'SNLFacEq1 is NOT defined!'
#endif
    CLOSE(Out,STATUS='KEEP')

    Gamma = GammaV; N_Neighbor = N_NeighborV; Alpha = AlphaV
    NLFac = NLFacV
#ifdef NLFacEq1
    IF(abs(NLFac-1.0D0) > 1.0D-15) THEN
      WRITE(*,*) &
        'ERROR: undefine the compilation switch NLFacEq1 or change NLFac to 1'
      STOP
    ENDIF
#endif
    NCity = NCityV

    CALL OpenASCII(OutFile,Out)
    WRITE(Out,'(A,E8.3,A,E8.3,A,I8)') &
      'Parameters: MaxCostFrac = ',MaxCostFrac,&
      ', TFacR = ', TFacR,',  Temp_Step =', Temp_Step
    WRITE(Out,'(A,I8,A,E8.3,A,E8.3)') &
      'Parameters: NCity = ',NCity,', Alpha = ',Alpha,', Gamma = ', Gamma
    WRITE(Out,'(A,E8.3,A,I10)') 'NLFac = ',NLFac,', N_Neighbor =', N_Neighbor
    CLOSE(Out,STATUS='KEEP')
    IF(NCity <= 4) THEN
      STOP 'Error in Anneal! Too few cities!'
    ENDIF

#if defined(SD)
    ! Num_Pair is wrong if NCity is larger than (sqrt(2^31) = 2**15.5 = 46340
    IF(NCity > 46300) THEN
      STOP 'ERROR: Integer multiplication overflows, reduce NCity!'
    ENDIF
    NumPair = (NCity*(NCity-1))/2
    WRITE(*,*) 'NCity = ',NCity, ', NumPair = ', NumPair
    CALL New(DistArr,NumPair)
    IndexInt = 0
    DO I = 1, NCity
      DO J = 1, I-1
        IndexInt = IndexInt + 1
        IF( ((I-1)*(I-2)/2 + J) /= IndexInt) THEN
          STOP 'ERROR: Indexing problem! Might be due to overflow !!'
        ENDIF
        DistArr%D(IndexInt) = ALen(RCoor%D(1,I),RCoor%D(2,I),RCoor%D(3,I),&
                                 RCoor%D(1,J),RCoor%D(2,J),RCoor%D(3,J))
      ENDDO
    ENDDO
    IF(IndexInt /= NumPair) THEN
      STOP 'ERROR: IndexInt must be the same as NumPair'
    ENDIF
#endif

    CALL GetEndCities1(RCoor,LeftCity,RightCity)
    DO I = 1, NCity
      TravO%I(I) = I
    ENDDO
    TravO%I(1) = LeftCity
    TravO%I(LeftCity) = 1
    TravO%I(NCity) = RightCity
    TravO%I(RightCity) = NCity

    I = TravO%I(1)
    J = TravO%I(NCity)
    DMax = ALen(RCoor%D(1,I),RCoor%D(2,I),RCoor%D(3,I),&
             RCoor%D(1,J),RCoor%D(2,J),RCoor%D(3,J))
    CALL OpenASCII(OutFile,Out)
    WRITE(*,*) 'DMax = ', DMax
    Beta = Gamma*(DMax/(2.0*NCity))**NLFac
    WRITE(*,*) 'Beta = ',Beta
    WRITE(Out,*) 'Beta = ',Beta
    Eta = Alpha/DMax
    WRITE(*,*) 'Eta = ',Eta
    WRITE(Out,*) 'Eta = ',Eta
    CLOSE(Out,STATUS='KEEP')

    FNumOpt = Get_Opt_Num()
    WRITE(*,*) 'NCity = ',NCity, ' FNumOpt = ',FNumOpt

    IF(N_Neighbor > 0 .AND. N_Neighbor < NCity) THEN
      CALL New(NeighborDis,NCity-1)
      CALL AveNeighborDis(RCoor,NeighborDis)
    ENDIF

    IF(N_Neighbor == 0) THEN
      RCut = 0.0D0 ! No off-diagonal cost
    ELSEIF( N_Neighbor > 0 .AND. N_Neighbor < NCity) THEN
      RCut = NeighborDis%D(N_Neighbor)
    ELSEIF(N_Neighbor .GE. NCity) THEN
      RCut = Big_DBL
    ELSE
      WRITE(*,*) 'ERROR: N_Neigbhor must be non-negative !'
    ENDIF

    CALL OpenASCII(OutFile,Out)
    WRITE(Out,'(A,I10,A,E10.5)') 'N_Neighbor = ',N_Neighbor,', RCut = ',RCut
    WRITE(*,'(A,I10,A,E10.5)') 'N_Neighbor = ',N_Neighbor,', RCut = ',RCut

    ! Calculate the length of the total path (main diagonal part)
    TotMDCost = CalTotMDCost(RCoor,TravO)

    ! Calculate the length of the total path (off-diagonal part)
    TotODCost = CalTotODCost(RCoor,TravO)
    TotCost = TotMDCost+TotODCost
    WRITE(Out,*) 'TotMDCost =',TotMDCost,', TotODCost =',TotODCost,', TotCost =',TotCost
    WRITE(*,*)
    CLOSE(Out,STATUS='KEEP')

    CALL New(ChosenCity,6)
    CALL New(XX,6); CALL New(YY,6); CALL New(ZZ,6)
    CALL New(NewTravO,NCity)
    CALL New(FCity,NCity)
    CALL New(MCity,NCity)

    WRITE(*,*) 'Get MaxIncCost...'
    MaxIncCost = -BIG_DBL

    Call Elapsed_Time(AnnealTime,'Init')
    DO I = 1,5000
      Dec = Random()
      IF(Dec < 0.5) THEN
        CALL RevCst(RCoor,TravO,NewTravO,ChosenCity,XX,YY,ZZ,&
          TotMDCost,IncMDCost,TotODCost,IncODCost,FCity,MCity)
        IncCost = IncMDCost + IncODCost
        IF(IncCost > MaxIncCost) THEN
          MaxIncCost = IncCost
        ENDIF
        CALL Revers(TravO,NewTravO,TotMDCost,IncMDCost,TotODCost,IncODCost)

        IF(I < 100) THEN
          TotMDCost1 = CalTotMDCost(RCoor,TravO)
          Okay = Check_Accuracy(TotMDCost,TotMDCost1)
          IF(.NOT. Okay) THEN
            STOP 'TotMDCost and TotMDCost1 differ!!'
          ENDIF
          TotODCost1 = CalTotODCost(RCoor,TravO)
          Okay = Check_Accuracy(TotODCost,TotODCost1)
          IF(.NOT. Okay) THEN
            STOP 'TotODCost and TotODCost1 differ!!'
          ENDIF
        ENDIF
      ELSE
        CALL TrnCst(RCoor,TravO,NewTravO,ChosenCity,&
          XX,YY,ZZ,TotMDCost,IncMDCost,TotODCost,IncODCost,FCity,MCity)
        IncCost = IncMDCost + IncODCost
        IF(IncCost > MaxIncCost) THEN
          MaxIncCost = IncCost
        ENDIF
        CALL Trnspt(TravO,NewTravO,TotMDCost,IncMDCost,TotODCost,IncODCost)

        IF(I < 100) THEN
          TotMDCost1 = CalTotMDCost(RCoor,TravO)
          Okay = Check_Accuracy(TotMDCost,TotMDCost1)
          IF(.NOT. Okay) THEN
            STOP 'TotMDCost and TotMDCost1 differ!!'
          ENDIF
          TotODCost1 = CalTotODCost(RCoor,TravO)
          Okay = Check_Accuracy(TotODCost,TotODCost1)
          IF(.NOT. Okay) THEN
            STOP 'TotODCost and TotODCost1 differ!!'
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    Call Elapsed_Time(AnnealTime,'Accum')
    Call PPrint(AnnealTime,Proc_O='Anneal:5000times took')

    CALL OpenASCII(OutFile,Out)
    WRITE(Out,*) 'NCity,MaxIncCost =',NCity,MaxIncCost
    WRITE(*,*) 'NCity,MaxIncCost =',NCity,MaxIncCost
    CLOSE(Out,STATUS='KEEP')

    CALL OpenASCII(OutFile,Out)
    WRITE(*,*) 'Use the last list for the annealing process.'
    WRITE(Out,*) 'Use the last list for the annealing process.'
    CLOSE(Out,STATUS='KEEP')


    NOver = 100*NCity
    NLimit = 10*NCity
    T = MaxIncCost*MaxCostFrac
    WRITE(*,*) 'Initial T is ', T
    WRITE(*,*) 'Start the annealing process...'
    Call Elapsed_Time(AnnealTime,'Init')

    CALL OpenASCII(RunTravOFile//Trim(IntToChar(0))//TravOExt,RunTravOU,NewFile_O=.TRUE.)
    WRITE(RunTravOU,*) NCity
    WRITE(RunTravOU,*) (TravO%I(I),I=1,NCity)
    CLOSE(RunTravOU)

    DO J = 1, Temp_Step
      NSucc = 0
      DO K = 1,NOver
        Dec = Random()
        IF(Dec < 0.5) THEN
          CALL RevCst(RCoor,TravO,NewTravO,ChosenCity,XX,YY,ZZ,&
            TotMDCost,IncMDCost,TotODCost,IncODCost,FCity,MCity)

          IncCost = IncMDCost + IncODCost
          Ans = (IncCost < 0.0) .OR. (Random().LT.EXP(-IncCost/T))
          IF(Ans) THEN
            NSucc = NSucc + 1
            CALL Revers(TravO,NewTravO,TotMDCost,IncMDCost,TotODCost,IncODCost)
          ENDIF
        ELSE
          CALL TrnCst(RCoor,TravO,NewTravO,ChosenCity,&
            XX,YY,ZZ,TotMDCost,IncMDCost,TotODCost,IncODCost,FCity,MCity)
          IncCost = IncMDCost + IncODCost
          Ans = (IncCost < 0.0) .OR. (Random().LT.EXP(-IncCost/T))
          IF(Ans) THEN
            NSucc = NSucc + 1
            CALL Trnspt(TravO,NewTravO,TotMDCost,IncMDCost,TotODCost,IncODCost)
          ENDIF
        ENDIF
        IF(NSucc .GE. NLimit) Exit
      ENDDO
      CALL OpenASCII(OutFile,Out)
      WRITE(Out,'(A,I5,A,F7.2,A,E10.5)') &
        'Temperature index = ', J, ' (',J*100./Temp_Step,'% ), T = ', T
      WRITE(Out,'(A,I6,A,F7.2,A,I8,A,F7.2,A)') &
        'K-1 = ', K-1,' (',(K-1)*100.0/NOver, '%), NSucc = ',&
         NSucc, ' (',NSucc*100.0/(K-1),'% acceptance)'
      TotCost = TotMDCost+TotODCost
      WRITE(Out,*) 'TotMDCost =',TotMDCost,', TotODCost =',TotODCost,', TotCost =',TotCost
      WRITE(Out,*)
      CLOSE(Out,STATUS='KEEP')
      WRITE(*,'(A,I5,A,F7.2,A,E10.5)') &
        'Temperature index = ', J, ' (',J*100./Temp_Step,'% ), T = ', T
      WRITE(*,'(A,I6,A,F7.2,A,I8,A,F7.2,A)') &
        'K-1 = ', K-1,' (',(K-1)*100.0/NOver, '%), NSucc = ',&
         NSucc, ' (',NSucc*100.0/(K-1),'% acceptance)'
      WRITE(*,*) 'TotMDCost =',TotMDCost,', TotODCost =',TotODCost
      WRITE(*,*)
      CALL OpenASCII(RunTravOFile//Trim(IntToChar(J))//TravOExt,RunTravOU,NewFile_O=.TRUE.)
      WRITE(RunTravOU,*) NCity
      WRITE(RunTravOU,*) (TravO%I(I),I=1,NCity)
      CLOSE(RunTravOU)

      T = T*TFacR
      IF(NSucc == 0) EXIT
    ENDDO
    Call Elapsed_Time(AnnealTime,'Accum')
    Call PPrint(AnnealTime,Proc_O='Anneal:find min took')

    CALL OpenASCII(OutFile,Out)
    WRITE(Out,*) 'Final temperature index J-1 = ', J-1
    WRITE(Out,*) 'Final success number index K-1 = ', K-1
    WRITE(*,*) 'Final temperature index J-1 = ', J-1
    WRITE(*,*) 'Final success number index K-1 = ', K-1

    ! stringent tests for the correctness of the subroutine.
    TotMDCost1 = CalTotMDCost(RCoor,TravO)
    WRITE(Out,*) 'Final main diagonal cost: TotMDCost,TotMDCost1 =',&
      TotMDCost,TotMDCost1
    Okay = Check_Accuracy(TotMDCost,TotMDCost1)
    IF(.NOT. Okay) THEN
      STOP 'TotMDCost and TotMDCost1 differ!!'
    ENDIF
    TotODCost1 = CalTotODCost(RCoor,TravO)
    WRITE(Out,*) 'Final off-diagonal cost: TotODCost,TotODCost1 =',&
      TotODCost,TotODCost1
    Okay = Check_Accuracy(TotODCost,TotODCost1)
    IF(.NOT. Okay) THEN
      STOP 'TotODCost and TotODCost1 differ!!'
    ENDIF
    TotCost = TotMDCost+TotODCost
    WRITE(Out,*) 'TotMDCost =',TotMDCost,', TotODCost =',TotODCost,', TotCost =',TotCost
    CLOSE(Out,STATUS='KEEP')

  END SUBROUTINE Anneal

!---------------------------------------------------------------------------
  FUNCTION Check_Accuracy(Val1,Val2)
  REAL(DOUBLE)::Val1,Val2,MaxV,MinV,Abs1,Abs2
  LOGICAL::Check_Accuracy
  MaxV = Max(Val1,Val2)
  MinV = Min(Val1,Val2)

  IF(MaxV == MinV) THEN
    Check_Accuracy = .TRUE.
  ELSE
    Abs1 = ABS(Val1)
    Abs2 = ABS(Val2)
    MaxV = MAX(Abs1,Abs2)
    IF( ABS( (Val1-Val2)/MaxV) < 1.0D-5) THEN
       Check_Accuracy = .TRUE.
    ELSE
       Check_Accuracy = .FALSE.
    ENDIF
  ENDIF
  END FUNCTION Check_Accuracy

!---------------------------------------------------------------------------
#ifdef SD
  ! I and J are the indices of a city in RCoor, don't swap I and J!
  FUNCTION TableDist(I,J)
    INTEGER::I,J
    REAL(DOUBLE)::TableDist

    IF(I == J) THEN
      STOP 'ERROR: I == J in TableDist!'
    ENDIF
    IF(I > J) THEN
      TableDist = DistArr%D( ((I-2)*(I-1))/2 + J)
    ELSE
      TableDist = DistArr%D( ((J-2)*(J-1))/2 + I)
    ENDIF
    END FUNCTION TableDist
#endif

!---------------------------------------------------------------------
! obtain the average distance for the n-th neighbour
  SUBROUTINE AveNeighborDis(RCoor,Dist)
    TYPE(DBL_RNK2)::RCoor
    INTEGER::I,J,IndexInt
    TYPE(DBL_VECT)::Dist,DisSum

    CALL New(DisSum,NCity-1)
    DO I = 1, NCity-1
      DisSum%D(I) = 0.0
    END DO
    DO I = 1,NCity
      IndexInt = 0
      DO J = 1, NCity
        IF(I == J) CYCLE
        IndexInt = IndexInt + 1
        Dist%D(IndexInt) = ALen(RCoor%D(1,I),RCoor%D(2,I),RCoor%D(3,I),&
                             RCoor%D(1,J),RCoor%D(2,J),RCoor%D(3,J))
      END DO
      IF(IndexInt /= (NCity-1)) THEN
        STOP 'ERROR: Neighbor problem, number of cities is not correct !'
      END IF
      CALL SortDouble(Dist,NCity-1)
      DO J = 1, NCity-1
        DisSum%D(I) = DisSum%D(I) + Dist%D(I)
      END DO
    END DO
    DO I = 1,NCity-1
      Dist%D(I) = DisSum%D(I)/(NCity*1.0)
    ENDDO
    CALL SortDouble(Dist,NCity-1)
    WRITE(*,*) 'after sorting: average neighbour distance:'
    DO I = 1,NCity-1
      WRITE(*,*) Dist%D(I)
    END DO
    WRITE(*,*)
  END SUBROUTINE AveNeighborDis

!---------------------------------------------------------------------
  SUBROUTINE SortDouble(Dis,N)
    TYPE(DBL_VECT)::Dis
    INTEGER::N,I,IP
    REAL(DOUBLE)::Val

    DO I = 2, N
      Val = Dis%D(I); IP = I
      DO
        IF(.NOT. (IP > 1 .AND. Val < Dis%D(IP-1))) EXIT
        Dis%D(IP) = Dis%D(IP-1)
        IP = IP-1
      ENDDO
      Dis%D(IP) = Val
    ENDDO
  END SUBROUTINE SortDouble

!---------------------------------------------------------------------
  ! the cost-difference calculation switchs to Fix-Move Strategy
  ! when FCNum (the number of fixed cities) is greater or equal to Get_Opt_Num
  FUNCTION Get_Opt_Num()
    INTEGER::I,J,Get_Opt_Num
    REAL(DOUBLE)::A,B,Ratio,NewCost,ODCost

    ODCost = NCity*(NCity*1.0-1.)/2.0
    DO I = 1, NCity
      A = I; B = NCity-I
      NewCost = 2.0*(A*B + B*(B-1.0)/2.0) ! factor 2 is for 2 lists
      Ratio = NewCost/ODCost
      IF(Ratio < 1.0) THEN
        J = I
        EXIT
      END IF
    END DO
    ! 0.05 is a value to compensate the array copying overhead, this value
    ! is not optimized.
    Get_Opt_Num = ( J/(NCity*1.0) + 0.05)*NCity
    WRITE(*,*) 'J = ', J, '  Get_Opt_Num = ', Get_Opt_Num
  END FUNCTION Get_Opt_Num

!---------------------------------------------------------------------
  SUBROUTINE TrnCst(RCoor,TravO,NewTravO,ChosenCity,&
    XX,YY,ZZ,TotMDCost,IncMDCost,TotODCost,IncODCost,FCity,MCity)
    TYPE(DBL_RNK2)::RCoor
    TYPE(DBL_VECT)::XX,YY,ZZ
    INTEGER::Val,IP,ArrLen,I,J,RanIndex,Index,FCNum,MCNum
    TYPE(INT_VECT)::TravO,NewTravO,ChosenCity,FCity,MCity
    REAL(DOUBLE)::IncMDCost,TotODCost,TotMDCost,IncODCost
    REAL(DOUBLE),EXTERNAL::Random

    ArrLen = NCity-2 ! the first and the last cities are fixed
    DO
      ChosenCity%I(1) = 2+ArrLen*Random()
      ChosenCity%I(2) = 2+ArrLen*Random()
      ChosenCity%I(3) = 2+ArrLen*Random()

      DO I = 2, 3
        Val = ChosenCity%I(I); IP = I
        DO
          IF(.NOT.(IP > 1 .AND. Val < ChosenCity%I(IP-1))) EXIT
          ChosenCity%I(IP) = ChosenCity%I(IP-1)
          IP = IP-1
        ENDDO
        ChosenCity%I(IP) = Val
      ENDDO

      IF(ChosenCity%I(1) /= ChosenCity%I(2) .AND. &
         ChosenCity%I(2) /= ChosenCity%I(3)) EXIT
    ENDDO
    ChosenCity%I(4) = ChosenCity%I(3)+1
    ChosenCity%I(5) = ChosenCity%I(1)-1
    ChosenCity%I(6) = ChosenCity%I(2)+1

    ! get NewTravO (i.e. new permutation)
    CALL NewTravO_Trnspt(TravO,ChosenCity,NewTravO)

    ! calculate the increase in main diagonal cost
#if !defined(SD)
    DO I = 1, 6
      J = TravO%I(ChosenCity%I(I))
      XX%D(I) = RCoor%D(1,J)
      YY%D(I) = RCoor%D(2,J)
      ZZ%D(I) = RCoor%D(3,J)
    ENDDO
    IncMDCost = Eta*(-ALen(XX%D(2),YY%D(2),ZZ%D(2),XX%D(6),YY%D(6),ZZ%D(6)) &
                       -ALen(XX%D(1),YY%D(1),ZZ%D(1),XX%D(5),YY%D(5),ZZ%D(5)) &
                       -ALen(XX%D(3),YY%D(3),ZZ%D(3),XX%D(4),YY%D(4),ZZ%D(4)) &
                       +ALen(XX%D(1),YY%D(1),ZZ%D(1),XX%D(3),YY%D(3),ZZ%D(3)) &
                       +ALen(XX%D(2),YY%D(2),ZZ%D(2),XX%D(4),YY%D(4),ZZ%D(4)) &
                       +ALen(XX%D(5),YY%D(5),ZZ%D(5),XX%D(6),YY%D(6),ZZ%D(6)))
#else
    I = TravO%I(ChosenCity%I(2))
    J = TravO%I(ChosenCity%I(6))
    IncMDCost = -TableDist(I,J)
    I = TravO%I(ChosenCity%I(1))
    J = TravO%I(ChosenCity%I(5))
    IncMDCost = IncMDCost - TableDist(I,J)
    I = TravO%I(ChosenCity%I(3))
    J = TravO%I(ChosenCity%I(4))
    IncMDCost = IncMDCost - TableDist(I,J)
    I = TravO%I(ChosenCity%I(1))
    J = TravO%I(ChosenCity%I(3))
    IncMDCost = IncMDCost + TableDist(I,J)
    I = TravO%I(ChosenCity%I(2))
    J = TravO%I(ChosenCity%I(4))
    IncMDCost = IncMDCost + TableDist(I,J)
    I = TravO%I(ChosenCity%I(5))
    J = TravO%I(ChosenCity%I(6))
    IncMDCost = IncMDCost + TableDist(I,J)
    IncMDCost = Eta*IncMDCost
#endif

    ! calculate the increase in off-diagonal cost
    MCNum = ChosenCity%I(3)-ChosenCity%I(1)+1
    FCNum = NCity-MCNum

    IF(FCNum < FNumOpt) THEN
      IncODCost = CalTotODCost(RCoor,NewTravO) - TotODCost
    ELSE
      Index = 0
      DO I = 1,ChosenCity%I(5)
        Index = Index + 1
        FCity%I(Index) = I
      END DO
      DO I = ChosenCity%I(4),NCity
        Index = Index + 1
        FCity%I(Index) = I
      ENDDO
      FCNum = Index

      Index = 0
      DO I = ChosenCity%I(1),ChosenCity%I(3)
        Index = Index + 1
        MCity%I(Index) = I
      END DO
      MCNum = Index
      IF( (FCNum + MCNum) /= NCity) THEN
        STOP 'ERROR: the numbers of cities differ!'
      END IF

      IncODCost = FM_ODCostDiff(FCNum,MCNum,FCity,MCity,TravO,NewTravO,RCoor)
    END IF
  END SUBROUTINE TrnCst

!---------------------------------------------------------------------
  FUNCTION City2CityODCost(RCoor,M,N,I,J)
    TYPE(DBL_RNK2)::RCoor
    INTEGER::M,N,I,J
    REAL(DOUBLE)::City2CityODCost,Dij

#ifdef SD
    Dij = TableDist(M,N)
#else
    Dij = ALen(RCoor%D(1,M),RCoor%D(2,M),RCoor%D(3,M),&
            RCoor%D(1,N),RCoor%D(2,N),RCoor%D(3,N))
#endif

    IF(Dij <= RCut) THEN
#ifdef NLFacEq1
      City2CityODCost = EXP( Beta*( ABS(I-J)/Dij))
#else
      City2CityODCost = EXP( Beta*( ABS(I-J)/Dij)**NLFac )
#endif
    ELSE
      City2CityODCost = 0.0d0
    ENDIF
  END FUNCTION City2CityODCost

!---------------------------------------------------------------------
  SUBROUTINE NewTravO_Trnspt(TravO,ChosenCity,NewTravO)
    TYPE(INT_VECT)::TravO,ChosenCity,NewTravO
    INTEGER::NN,I,J,IndexInt

    NN = ChosenCity%I(5)
    IndexInt = 1
    DO I = 1, NN
      NewTravO%I(IndexInt) = TravO%I(I)
      IndexInt = IndexInt + 1
    ENDDO

    NN = ChosenCity%I(3)-ChosenCity%I(6)+1
    J = ChosenCity%I(6)
    DO I = 1,NN
      NewTravO%I(IndexInt) = TravO%I(J+I-1)
      IndexInt = IndexInt + 1
    ENDDO

    NN = ChosenCity%I(2)-ChosenCity%I(1)+1
    J = ChosenCity%I(1)
    DO I = 1,NN
      NewTravO%I(IndexInt) = TravO%I(J+I-1)
      IndexInt = IndexInt + 1
    ENDDO

    NN = NCity-ChosenCity%I(4)+1
    J = ChosenCity%I(4)
    DO I = 1, NN
      NewTravO%I(IndexInt) = TravO%I(J+I-1)
      IndexInt = IndexInt + 1
    ENDDO

    IF(IndexInt /= (NCity+1)) THEN
      STOP 'missing cities !!'
    ENDIF
  END SUBROUTINE NewTravO_Trnspt

!---------------------------------------------------------------------
  SUBROUTINE Trnspt(TravO,NewTravO,TotMDCost,IncMDCost,TotODCost,IncODCost)
    TYPE(INT_VECT)::TravO,NewTravO
    INTEGER::I
    REAL(DOUBLE)::TotMDCost,IncMDCost,TotODCost,IncODCost

    DO I = 1, NCity
      TravO%I(I) = NewTravO%I(I)
    ENDDO
    TotODCost = TotODCost + IncODCost
    TotMDCost = TotMDCost + IncMDCost
  END SUBROUTINE Trnspt

!----------------------------------------------------------------------
  SUBROUTINE RevCst(RCoor,TravO,NewTravO,ChosenCity,XX,YY,ZZ,&
    TotMDCost,IncMDCost,TotODCost,IncODCost,FCity,MCity)
    TYPE(DBL_RNK2)::RCoor
    TYPE(DBL_VECT)::XX,YY,ZZ
    INTEGER::Index,TmpInt,ArrLen,I,J,RanIndex1,RanIndex2,FCNum,MCNum
    TYPE(INT_VECT)::TravO,NewTravO,ChosenCity,FCity,MCity
    REAL(DOUBLE)::TotODCost,IncODCost,TotMDCost,IncMDCost
    REAL(DOUBLE),EXTERNAL::Random

    ArrLen = NCity-2
    DO
      RanIndex1 = 2 + ArrLen*Random()
      RanIndex2 = 2 + ArrLen*Random()
      IF(RanIndex1 /= RanIndex2) EXIT
    ENDDO
    IF(RanIndex1 > RanIndex2) THEN
      TmpInt = RanIndex1
      RanIndex1 = RanIndex2
      RanIndex2 = TmpInt
    ENDIF

    ChosenCity%I(1) = RanIndex1
    ChosenCity%I(2) = RanIndex2

    ChosenCity%I(3) = ChosenCity%I(1)-1
    ChosenCity%I(4) = ChosenCity%I(2)+1

    ! Create a new TravO (i.e. Permutation)
    DO I = 1, NCity
      NewTravO%I(I) = TravO%I(I)
    END DO
    CALL NewTravO_Revers(NewTravO,ChosenCity)

    ! calculate the increase in main diagonal cost
#if !defined(SD)
    DO I = 1, 4
      J = TravO%I(ChosenCity%I(I))
      XX%D(I) = RCoor%D(1,J)
      YY%D(I) = RCoor%D(2,J)
      ZZ%D(I) = RCoor%D(3,J)
    ENDDO
    IncMDCost =  Eta*(-ALen(XX%D(1),YY%D(1),ZZ%D(1),XX%D(3),YY%D(3),ZZ%D(3)) &
                        -ALen(XX%D(2),YY%D(2),ZZ%D(2),XX%D(4),YY%D(4),ZZ%D(4)) &
                        +ALen(XX%D(1),YY%D(1),ZZ%D(1),XX%D(4),YY%D(4),ZZ%D(4)) &
                        +ALen(XX%D(2),YY%D(2),ZZ%D(2),XX%D(3),YY%D(3),ZZ%D(3)))
#else
    I = TravO%I(ChosenCity%I(1))
    J = TravO%I(ChosenCity%I(3))
    IncMDCost = -TableDist(I,J)
    I = TravO%I(ChosenCity%I(2))
    J = TravO%I(ChosenCity%I(4))
    IncMDCost = IncMDCost - TableDist(I,J)
    I = TravO%I(ChosenCity%I(1))
    J = TravO%I(ChosenCity%I(4))
    IncMDCost = IncMDCost + TableDist(I,J)
    I = TravO%I(ChosenCity%I(2))
    J = TravO%I(ChosenCity%I(3))
    IncMDCost = IncMDCost + TableDist(I,J)
    IncMDCost = Eta*IncMDCost
#endif

    ! calculate the increase in off-diagonal cost
    MCNum = ChosenCity%I(2)-ChosenCity%I(1)+1
    FCNum = NCity-MCNum

    IF(FCNum < FNumOpt) THEN
      IncODCost = CalTotODCost(RCoor,NewTravO) - TotODCost
    ELSE
      ! again, calculate the local change only
      Index = 0
      DO I = 1, ChosenCity%I(3)
        Index = Index + 1
        FCity%I(Index) = I
      ENDDO
      DO I = ChosenCity%I(4), NCity
        Index = Index + 1
        FCity%I(Index) = I
      ENDDO
      FCNum = Index

      Index = 0
      DO I = ChosenCity%I(1),ChosenCity%I(2)
        Index = Index + 1
        MCity%I(Index) = I
      ENDDO

      MCNum = Index
      IF( (FCNum + MCNum) /= NCity) THEN
        STOP 'ERROR: the numbers of cities differ!'
      END IF

      IncODCost = FM_ODCostDiff(FCNum,MCNum,FCity,MCity,TravO,NewTravO,RCoor)

    ENDIF
  END SUBROUTINE RevCst

!----------------------------------------------------------------------
! F=> Fixed Cities, M => Moved Cities
  FUNCTION FM_ODCostDiff(FCNum,MCNum,FCity,MCity,TravO,NewTravO,RCoor)
    TYPE(DBL_RNK2)::RCoor
    TYPE(INT_VECT)::FCity,MCity,TravO,NewTravO
    REAL(DOUBLE)::FM_ODCostDiff,OldODCost,NewODCost
    INTEGER::FCNum,MCNum,I,J,II,JJ,M,N,NewM,NewN

    ! now calculate the F and M interaction energy
    OldODCost = 0.0
    NewODCost = 0.0
    DO I = 1, FCNum
      DO J = 1,MCNum
        II = FCity%I(I) ! II is the index in TravO
        JJ = MCity%I(J) ! JJ is the index in TravO
        M = TravO%I(II) ! M is physical index in RCoor
        N = TravO%I(JJ) ! N is a physical index in RCoor
        OldODCost = OldODCost + City2CityODCost(RCoor,M,N,II,JJ)
        NewN = NewTravO%I(JJ)
        NewODCost = NewODCost + City2CityODCost(RCoor,M,NewN,II,JJ)
      END DO
    END DO

    ! self-interacting cost
    DO I = 1,MCNum
      DO J = I+1,MCNum
        II = MCity%I(I) !! II is the index in TravO
        JJ = MCity%I(J) !! JJ is the index in TravO
        M = TravO%I(II) ! M is physical index in RCoor
        N = TravO%I(JJ) ! N is a physical index in RCoor
        OldODCost = OldODCost + City2CityODCost(RCoor,M,N,II,JJ)
        NewM = NewTravO%I(II)
        NewN = NewTravO%I(JJ)
        NewODCost = NewODCost + City2CityODCost(RCoor,NewM,NewN,II,JJ)
      END DO
    END DO

    FM_ODCostDiff = NewODCost - OldODCost
  END FUNCTION FM_ODCostDiff

!----------------------------------------------------------------------
  SUBROUTINE NewTravO_Revers(TravO,ChosenCity)
    INTEGER::NN,I,TmpInt,Left,Right
    TYPE(INT_VECT)::TravO,ChosenCity

    NN = (ChosenCity%I(2)-ChosenCity%I(1)+1)/2
    DO I = 1,NN
      Left = ChosenCity%I(1)+I-1
      Right = ChosenCity%I(2)-I+1
      TmpInt = TravO%I(Left)
      TravO%I(Left) = TravO%I(Right)
      TravO%I(Right) = TmpInt
    ENDDO
  END SUBROUTINE NewTravO_Revers

!----------------------------------------------------------------------
  SUBROUTINE Revers(TravO,NewTravO,TotMDCost,IncMDCost,TotODCost,IncODCost)
    INTEGER::I
    TYPE(INT_VECT)::TravO,NewTravO
    REAL(DOUBLE)::TotODCost,IncODCost,TotMDCost,IncMDCost

    DO I = 1,NCity
      TravO%I(I) = NewTravO%I(I)
    ENDDO
    TotODCost = TotODCost + IncODCost
    TotMDCost = TotMDCost + IncMDCost
  END SUBROUTINE Revers

!----------------------------------------------------------------------
  ! the order is coordinates of 1, and followed by 2
  FUNCTION ALen(x1,y1,z1,x2,y2,z2)
    REAL(DOUBLE)::x1,y1,z1,x2,y2,z2,ALen,dx,dy,dz
    dx = x1-x2; dy = y1-y2; dz = z1-z2
    ALen = SQRT(dx*dx + dy*dy + dz*dz)
  END FUNCTION ALen

!----------------------------------------------------------------------
  FUNCTION CalTotMDCost(RCoor,TravO)
    TYPE(DBL_RNK2)::RCoor
    TYPE(INT_VECT)::TravO
    INTEGER::I,M,N
    REAL(DOUBLE)::CalTotMDCost

    CalTotMDCost = 0.0
    DO I = 1, NCity-1
      M = TravO%I(I)
      N = TravO%I(I+1)
#if !defined(SD)
      CalTotMDCost = CalTotMDCost + &
        ALen(RCoor%D(1,M),RCoor%D(2,M),RCoor%D(3,M),&
             RCoor%D(1,N),RCoor%D(2,N),RCoor%D(3,N))
#else
      CalTotMDCost = CalTotMDCost + TableDist(M,N)
#endif
    END DO
    CalTotMDCost = Eta*CalTotMDCost
  END FUNCTION CalTotMDCost

!----------------------------------------------------------------------
  FUNCTION CalTotODCost(RCoor,TravO)
    TYPE(DBL_RNK2)::RCoor
    TYPE(INT_VECT)::TravO
    INTEGER::I,J,M,N
    REAL(DOUBLE)::CalTotODCost,Dij

    CalTotODCost = 0.0
    DO I = 1, NCity
      M = TravO%I(I)
      DO J = I+1, NCity
        N = TravO%I(J)
       CalTotODCost = CalTotODCost + City2CityODCost(RCoor,M,N,I,J)
      END DO
    END DO
  END FUNCTION CalTotODCost

!----------------------------------------------------------------------
  ! get the end cities using the tight bounding box
  SUBROUTINE GetEndCities1(RCoor,LeftCity,RightCity)
    TYPE(DBL_RNK2)::RCoor
    REAL(DOUBLE)::Dis,MinLDis,MinRDis,X,Y,Z,XLEx,YLEx,ZLEx,XREx,YREx,ZREx
    INTEGER::I,LeftCity,RightCity

    XLEx = Big_DBL;  XREx = -Big_DBL
    YLEx = Big_DBL;  YREx = -Big_DBL
    ZLEx = Big_DBL;  ZREx = -Big_DBL

    DO I = 1, NCity
      X = RCoor%D(1,I); Y = RCoor%D(2,I); Z = RCoor%D(3,I)
      IF(X < XLEx) XLEx = X; IF(X > XREx) XREx = X
      IF(Y < YLEx) YLEx = Y; IF(Y > YREx) YREx = Y
      IF(Z < ZLEx) ZLEx = Z; IF(Z > ZREx) ZREx = Z
    ENDDO
    WRITE(*,*) 'XEx : ',XLEx, XREx
    WRITE(*,*) 'YEx : ',YLEx, YREx
    WRITE(*,*) 'ZEx : ',ZLEx, ZREx

    MinLDis = Big_DBL; MinRDis = Big_DBL

    DO I = 1, NCity
      Dis = ALen(RCoor%D(1,I),RCoor%D(2,I),RCoor%D(3,I),XLEx,YLEx,ZLEx)
      IF(Dis < MinLDis) THEN
        MinLDis = Dis
        LeftCity = I
      ENDIF
      Dis = ALen(RCoor%D(1,I),RCoor%D(2,I),RCoor%D(3,I),XREx,YREx,ZREx)
      IF(Dis < MinRDis) THEN
        MinRDis = Dis
        RightCity = I
      ENDIF
    ENDDO

    WRITE(*,*) 'MinLDis = ',MinLDis
    WRITE(*,*) 'MinRDis = ',MinRDis
    WRITE(*,*) 'LeftCity,RightCity = ',LeftCity,RightCity
  END SUBROUTINE GetEndCities1

!----------------------------------------------------------------------
  ! find out the largest end-to-end distance
  SUBROUTINE GetEndCities2(RCoor,LeftCity,RightCity)
    TYPE(DBL_RNK2)::RCoor
    INTEGER::I,J,LeftCity,RightCity
    REAL(DOUBLE)::DisIJ,MaxDis

    MaxDis = 0.0
    DO I = 1, NCity
      DO J = I+1, NCity
        DisIJ = ALen(RCoor%D(1,I),RCoor%D(2,I),RCoor%D(3,I),&
          RCoor%D(1,J),RCoor%D(2,J),RCoor%D(3,J))
        IF(DisIJ .GE. MaxDis) THEN
          MaxDis = DisIJ
          LeftCity = I
          RightCity = J
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE GetEndCities2

END MODULE Opt_Trav_Band

!---------------------------------------------------------------------------
MODULE AnnealMap
  USE PrettyPrint
  IMPLICIT NONE
  INTEGER,PRIVATE::NCity
CONTAINS
!---------------------------------------------------------------------------
  SUBROUTINE TableAnneal(RCoor,NCityV,TravO)
    TYPE(DBL_RNK2)::RCoor
    INTEGER(INT8),ALLOCATABLE::CityRank(:)
    TYPE(INT_VECT)::TravO
    TYPE(INT_RNK3)::AnnealKey
    INTEGER::Imax,IntVect(3),NCityV,I,J,K,LinDim,LM1,ReadNCity,M,Rank,IndexInt
    REAL(DOUBLE)::Ratio,RMin(3),Diff,MaxDiff

!#define USE_Final_TravO_Dat
#undef USE_Final_TravO_Dat

#ifdef  USE_Final_TravO_Dat
    INTEGER,PARAMETER::ReadU=30
    INTEGER,ALLOCATABLE::ReadTravO(:)
#else
    INCLUDE 'Final_TravO.Inc'
#endif

#ifdef  USE_Final_TravO_Dat
    write(*,*) 'Final_TravO.dat is used.'
#else
    write(*,*) 'Final_TravO.Inc is used.'
#endif

    NCity = NCityV
#ifdef USE_Final_TravO_Dat
    CALL OpenASCII('Final_TravO.dat',ReadU,OldFileQ_O=.TRUE.,Rewind_O=.TRUE.)
    Read(ReadU,*) ReadNCity
    WRITE(*,*) 'TableAnneal, ReadNCity = ',ReadNCity
#else
    ! assign ReadNCity
    ReadNCity = NCity_Const
#endif
    LinDim = NINT(ReadNCity**(1.0d0/3.0d0))
!   WRITE(*,*) 'LinDim = ', LinDim
    IF(LinDim*LinDim*LinDim /= ReadNCity) THEN
      STOP 'Error: Cube root problem in TableAnneal!'
    ENDIF
    LM1 = LinDim-1
!   WRITE(*,*) 'LinDim = ', LinDim, ', Lm1=', LM1
    CALL New(AnnealKey,(/LM1,LM1,LM1/),(/0,0,0/))

#ifdef USE_Final_TravO_Dat
    ALLOCATE(ReadTravO(ReadNCity))
    Read(ReadU,*) (ReadTravO(I),I=1,ReadNCity)
#endif

    IndexInt = 0
    DO I = 0,LM1
      DO J = 0, LM1
        DO K = 0, LM1
          IndexInt = IndexInt + 1
          Rank = 0
          DO M = 1, ReadNCity
            IF(ReadTravO(M) == IndexInt) THEN
              Rank = M
              EXIT
            ENDIF
          ENDDO
          IF(Rank < 1 .OR. Rank > ReadNCity) THEN
            WRITE(*,*) 'IndexInt = ',IndexInt
            WRITE(*,*) 'Rank = ',Rank
            STOP 'Error: Some serious problem in inverse mapping!'
          ENDIF
          AnnealKey%I(K,J,I) = Rank
        ENDDO
      ENDDO
    ENDDO

    DO I = 0,LM1
      DO J = 0, LM1
        DO K = 0, LM1
!         WRITE(*,*) 'AnnealKey',K,J,I,AnnealKey%I(K,J,I)
        ENDDO
      ENDDO
    ENDDO
#ifdef USE_Final_TravO_Dat
    CLOSE(ReadU)
#endif

    RMin(1) = Big_DBL
    RMin(2) = Big_DBL
    RMin(3) = Big_DBL
    DO I = 1, NCity
      DO J = 1,3
        RMin(J) = DMin1(RMin(J),RCoor%D(J,I))
      ENDDO
      TravO%I(I) = I
    ENDDO
    ! any -ve value for MaxDiff is okay, since we deal with the 1st quadrant
    MaxDiff = -1.0D0
    DO I = 1,NCity
      DO J = 1,3
        Diff = RCoor%D(J,I)-RMin(J)
        MaxDiff = DMax1(MaxDiff,Diff)
      ENDDO
    ENDDO
    RATIO = (LinDim*1.D0)/MaxDiff
    IMax = -Big_Int
    ALLOCATE(CityRank(1:NCity))
    DO I = 1,NCity
      DO J = 1, 3
        IntVect(J) = (RCoor%D(J,I)-RMin(J))*Ratio !truncate the fraction part.
        IMax = Max(IMax,IntVect(J))
        IF(IntVect(J) == LinDim) THEN
          IntVect(J) = IntVect(J)-1
        ENDIF
        IF(IntVect(J) < 0 .OR. IntVect(J) .GE. LinDim) THEN
          WRITE(*,*) 'IntVect',j, ' =',IntVect(j)
          STOP 'ERROR: An error has occured! array bound problem'
        ENDIF
      ENDDO
      CityRank(I) = AnnealKey%I(IntVect(1),IntVect(2),IntVect(3))
!     WRITE(*,*) 'CityRank',I, ' =', CityRank(I)
    ENDDO

    CALL I8Sort(CityRank,TravO%I,NCity,2)
    CALL Delete(AnnealKey)
#ifdef USE_Final_TravO_Dat
    DEALLOCATE(ReadTravO)
#endif

  END SUBROUTINE TableAnneal
END MODULE AnnealMap
