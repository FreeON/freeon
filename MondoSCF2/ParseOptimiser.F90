MODULE ParseOptimizer
  IMPLICIT NONE
!---------------------------------------------------------------------------- 
      SUBROUTINE ParseGrad(Ctrl)
         TYPE(SCFControls)          :: Ctrl
         TYPE(TOLS)                 :: Thrsh ! Thresholds
         REAL(DOUBLE),DIMENSION(3)  :: Sum
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line
         TYPE(CHR_VECT)                 :: Chars
!
!
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Max_Steps = 'Max_Steps'
!
!----------------------------------------------------------------------------
         CALL OpenASCII(OutFile,Out)
         CALL OpenASCII(InpFile,Inp)
!
!        Parse <OPTIONS> for <Grad=>
!     
         IF(OptKeyQ(Inp,GRADIENTS,FORCE))THEN
            Ctrl%Grad=GRAD_ONE_FORCE
            Ctrl%NGeom=1
         ELSEIF(OptKeyQ(Inp,DYNAMICS,MD_VELOCITY))THEN
            write(6,*)'Doing molecular dynamics with velocity verlet'
            Ctrl%Grad=GRAD_MD
            Ctrl%MDC%ALGORITHM=1
            Call ParseMD(Ctrl)
         ELSEIF(OptKeyQ(Inp,OPTIMIZATION,OPT_QUNew))THEN
            IF(OptKeyQ(Inp,OPTIMIZATION,OPT_ONE_BASE))THEN
               Ctrl%Grad=GRAD_QNew_ONE_OPT
            ELSE
               Ctrl%Grad=GRAD_QNew_OPT
            ENDIF
            IF(.NOT.OptIntQ(Inp,Max_Steps,Ctrl%NGeom))THEN
               Ctrl%NGeom=500
            ENDIF
         ELSEIF(OptKeyQ(Inp,OPTIMIZATION,OPT_StpDesc))THEN
               Ctrl%Grad=GRAD_StpDesc_OPT
         ELSEIF(OptKeyQ(Inp,OPTIMIZATION,OPT_DiagHess))THEN
               Ctrl%Grad=GRAD_DiagHess_OPT
         ELSEIF(OptKeyQ(Inp,OPTIMIZATION,OPT_TSTATE))THEN
            Ctrl%Grad=GRAD_TS_SEARCH
            IF(.NOT.OptIntQ(Inp,Max_Steps,Ctrl%NGeom))THEN
               Ctrl%NGeom=500
            ENDIF
         ELSE
            Ctrl%Grad=GRAD_NO_GRAD
         ENDIF
         IF(Ctrl%Grad>=GRAD_QNew_OPT)THEN
            IF(.NOT.OptIntQ(Inp,Max_Steps,Ctrl%NGeom))THEN
               Ctrl%NGeom=500
            ENDIF
         ENDIF



!
! Parse for GDfIIS coordinate type
!
             Ctrl%GeOp%NoGDIIS=.FALSE.
             Ctrl%GeOp%GDIISCoordType=OPT_CartDIIS
         IF(OptKeyQ(Inp,OPTIMIZATION,OPT_CartDIIS)) THEN
            Ctrl%GeOp%GDIISCoordType=OPT_CartDIIS
         ELSE IF(OptKeyQ(Inp,OPTIMIZATION,OPT_IntDIIS)) THEN
            Ctrl%GeOp%GDIISCoordType=OPT_IntDIIS
         ELSE
           IF(OptKeyQ(Inp,OPTIMIZATION,OPT_NoGDIIS)) THEN
             Ctrl%GeOp%NoGDIIS=.TRUE.
           ENDIF
         ENDIF
!
! Parse for projecting out rotations and translations
! from geometry displacements.
!
             Ctrl%GeOp%DoRotOff=.TRUE.
             Ctrl%GeOp%DoTranslOff=.TRUE.
         IF(OptKeyQ(Inp,OPTIMIZATION,OPT_NoRotOff)) THEN
            Ctrl%GeOp%DoRotOff=.FALSE.
         ELSE IF(OptKeyQ(Inp,OPTIMIZATION,OPT_NoTranslOff)) THEN
            Ctrl%GeOp%DoTranslOff=.FALSE.
         ENDIF
!
! Parse for Steepest descent Inverse Hessian
!
              Ctrl%GeOp%StpDescInvH=0.5 !default value
         IF(FindKey(STPDESCINVH,Inp)) THEN
           REWIND(UNIT=Inp)
           DO
              READ(Inp,DEFAULT_CHR_FMT,END=99) Line
              IF(INDEX(Line,STPDESCINVH)/=0) THEN
                CALL New(Chars,2)
                CALL LineToChars(Line,Chars)
                CALL LineToDbls(Chars%C(2),1,Sum(1))
                CALL Delete(Chars)
                Ctrl%GeOp%StpDescInvH=Sum(1)
              ENDIF
           ENDDO
99         CONTINUE
         ENDIF
!
! Parse for coordtype
!
         IF(OptKeyQ(Inp,COORDTYPE,CoordType_PrimInt)) THEN
            Ctrl%GeOp%CoordType=CoordType_PrimInt
         ELSE IF(OptKeyQ(Inp,COORDTYPE,CoordType_Cartesian)) THEN
            Ctrl%GeOp%CoordType=CoordType_Cartesian
         ELSE
            Ctrl%GeOp%CoordType=CoordType_PrimInt
         ENDIF
!
        CLOSE(Out,STATUS='KEEP')
        CLOSE(Inp,STATUS='KEEP')
!
CONTAINS
  
END MODULE ParseOptimizer
