      SUBROUTINE GetAdrL(I,J,Ind,JA,IA)
c---------------------------------------------------------------------- 
c Finds the CSR address corresponding to the full matrix index pair.
c Assumes that the collumn index array is not in order (linear search)
c Inputs:
c   I   = row index of full matrix
c   J   = collumn index of full matrix
c   JA  = array of collumn indices of sparse matrix
c   IA  = array of pointers to each row in sparse matrix
c Outputs:
c   Ind = sparse matrix index corresponding to (I,J)
c---------------------------------------------------------------------- 
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER JA(*),IA(*)
      INTEGER I,J,Ind,L
      DO L=IA(I),IA(I+1)-1
         IF(JA(L).EQ.J) GOTO 1000
      ENDDO
      CALL MondoHalt('Matrix element not found in GetAdr')
c      Ind=0
c      RETURN
 1000 CONTINUE
      Ind=L
      RETURN
      END

      SUBROUTINE GetAdrB(I,J,Ind,IType,JA,IA)
c----------------------------------------------------------------------
c Finds the CSR address corresponding to the full matrix index pair.
c Assumes that the collumn index array is in order (binary search)
c Inputs:
c   I     = row index of full matrix
c   J     = collumn index of full matrix
c   IType = 0 halt if element not found
c         = 1 return index 0 if not found
c   JA    = array of collumn indices of sparse matrix
c   IA    = array of pointers to each row in sparse matrix
c Outputs:
c   Ind   = sparse matrix index corresponding to (I,J)
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER I,J,Ind,Mid,IType
      INTEGER IBgn,IEnd
      INTEGER JA(*),IA(*)

      IBgn=IA(I)
      IEnd=IA(I+1)-1
      Ind=0
      IF(IType.EQ.0) THEN
 1000    CONTINUE
         Mid=(IEnd+IBgn)/2
         IF(JA(Mid).EQ.J) GOTO 2000
         IF(IBgn.GE.IEnd) THEN
            WRITE(*,*) "I=",I," J=",J
            CALL MondoHalt('Matrix element not found in GetAdrB')
         ENDIF
         IF(JA(Mid).GT.J) THEN
            IEnd=Mid-1
         ELSE
            IBgn=Mid+1
         ENDIF
         GOTO 1000
 2000    CONTINUE
         Ind=Mid
         RETURN
      ELSEIF(IType.EQ.1) THEN
 3000    CONTINUE
         Mid=(IEnd+IBgn)/2
         IF(JA(Mid).EQ.J) GOTO 4000
         IF(IBgn.GE.IEnd) RETURN
         IF(JA(Mid).GT.J) THEN
            IEnd=Mid-1
         ELSE
            IBgn=Mid+1
         ENDIF
         GOTO 3000
 4000    CONTINUE
         Ind=Mid
         RETURN
      ELSE
         WRITE(*,*) "IType=",IType
         CALL MondoHalt('Illegal IType in GetAdrB')
      ENDIF
      RETURN
      END
