MODULE Moments
  USE DerivedTypes
  USE GlobalScalars
  USE InOut
  USE PrettyPrint
  USE SpecFun
  IMPLICIT NONE
!----------------------------------------------------------
! TYPE SpMoments, Holds Info on Spherical Multipole Moments
!----------------------------------------------------------
  TYPE SpMoments
     INTEGER            :: Alloc
     INTEGER            :: MaxL
     INTEGER            :: MaxLM
     REAL(DOUBLE)       :: CenterX
     REAL(DOUBLE)       :: CenterY
     REAL(DOUBLE)       :: CenterZ    
     TYPE(DBL_VECT)     :: CMMat
     TYPE(DBL_VECT)     :: SMMat
  ENDTYPE SpMoments
!----------------------------------------------------------
! TYPE QMoments, Holds Info on Cartesian Multipole Moments
!----------------------------------------------------------
  TYPE QMoments
     INTEGER            :: Alloc
     INTEGER            :: MaxL
     INTEGER            :: MaxLen
     REAL(DOUBLE)       :: CenterX
     REAL(DOUBLE)       :: CenterY
     REAL(DOUBLE)       :: CenterZ    
     TYPE(DBL_VECT)     :: QMat
  ENDTYPE QMoments
CONTAINS
!========================================================================================
! ALLOCATE  new moments
!========================================================================================
  SUBROUTINE New_SpMoments(A,N_O)
    TYPE(SpMoments)                 :: A
    INTEGER,OPTIONAL                :: N_O
    IF(AllocQ(A%Alloc)) THEN
       CALL Delete_SpMoments(A)       
       IF(PRESENT(N_O)) THEN
          A%MaxL  = N_O
          A%MaxLM = LSP(N_O)
       ELSE
          A%MaxL  = 0        
          A%MaxLM = 0
       ENDIF
       CALL New(A%CMMat,A%MaxLM,0)  
       CALL New(A%SMMat,A%MaxLM,0)     
    ELSE
       A%Alloc=ALLOCATED_TRUE
       IF(PRESENT(N_O)) THEN
          A%MaxL  = N_O
          A%MaxLM = LSP(N_O)
       ELSE
          A%MaxL  = 0
          A%MaxLM = 0          
       ENDIF
       CALL New(A%CMMat,A%MaxLM,0)  
       CALL New(A%SMMat,A%MaxLM,0) 
    ENDIF
  END SUBROUTINE New_SpMoments
!
!
  SUBROUTINE New_QMoments(A,N_O)
    TYPE(QMoments)                  :: A
    INTEGER,OPTIONAL                :: N_O
    IF(AllocQ(A%Alloc)) THEN
       CALL Delete_QMoments(A)       
       IF(PRESENT(N_O)) THEN
          A%MaxL   = N_O
          A%MaxLen = LHGTF(N_O)
       ELSE
          A%MaxL   = 0        
          A%MaxLen = 0
       ENDIF
       CALL New(A%QMat,A%MaxLen,0)     
    ELSE
       A%Alloc=ALLOCATED_TRUE
       IF(PRESENT(N_O)) THEN
          A%MaxL   = N_O
          A%MaxLen = LHGTF(N_O)
       ELSE
          A%MaxL   = 0        
          A%MaxLen = 0
       ENDIF
       CALL New(A%QMat,A%MaxLen,0) 
    ENDIF
  END SUBROUTINE New_QMoments
!========================================================================================
! Delete  moments
!========================================================================================
  SUBROUTINE Delete_SpMoments(A)
    TYPE(SpMoments)                 :: A
    IF(AllocQ(A%Alloc)) THEN
       A%Alloc=ALLOCATED_FALSE
       A%MaxL  = 0
       A%MaxLM = 0
       CALL Delete(A%CMMat)
       CALL Delete(A%SMMat)
    ENDIF
  END SUBROUTINE Delete_SpMoments
!
!
  SUBROUTINE Delete_QMoments(A)
    TYPE(QMoments)                 :: A
    IF(AllocQ(A%Alloc)) THEN
       A%Alloc=ALLOCATED_FALSE
       A%MaxL   = 0
       A%MaxLen = 0
       CALL Delete(A%QMat)
    ENDIF
  END SUBROUTINE Delete_QMoments
!========================================================================================
! Get  moments
!========================================================================================
  SUBROUTINE Put_SpMoments(A,Tag_O)
    TYPE(SpMoments)                 :: A 
    CHARACTER(LEN=*),OPTIONAL       :: Tag_O
    IF(A%Alloc .NE. 0) THEN
       CALL Put(A%MaxL   ,'MaxL')
       CALL Put(A%MaxLM  ,'MaxLM')
       CALL Put(A%CenterX,'CenterX')
       CALL Put(A%CenterY,'CenterY')
       CALL Put(A%CenterZ,'CenterZ')
       CALL Put(A%CMMat,  'CMMat')
       CALL Put(A%SMMat,  'SMMat')
    ELSE
       CALL Halt('Error in Put_SpMoments: Not Allocated')
    ENDIF
  END SUBROUTINE Put_SpMoments
!
!
  SUBROUTINE Put_QMoments(A,Tag_O)
    TYPE(QMoments)                  :: A 
    CHARACTER(LEN=*),OPTIONAL       :: Tag_O
    IF(A%Alloc .NE. 0) THEN
       CALL Put(A%MaxL    ,'MaxL')
       CALL Put(A%MaxLen  ,'MaxLen')
       CALL Put(A%CenterX,'CenterX')
       CALL Put(A%CenterY,'CenterY')
       CALL Put(A%CenterZ,'CenterZ')
       CALL Put(A%QMat,   'QMat')
    ELSE
       CALL Halt('Error in Put_QMoments: Not Allocated')
    ENDIF
  END SUBROUTINE Put_QMoments

!========================================================================================
! Put  moments
!========================================================================================
  SUBROUTINE Get_SpMoments(A,Tag_O)
    TYPE(SpMoments)                 :: A
    CHARACTER(LEN=*),OPTIONAL       :: Tag_O
    CALL Get(A%MaxL   ,'MaxL')
    CALL Get(A%MaxLM  ,'MaxLM')
    CALL Get(A%CenterX,'CenterX')
    CALL Get(A%CenterY,'CenterY')
    CALL Get(A%CenterZ,'CenterZ')
    CALL Get(A%CMMat,  'CMMat')
    CALL Get(A%SMMat,  'SMMat')
  END SUBROUTINE Get_SpMoments
!
!
  SUBROUTINE Get_QMoments(A,Tag_O)
    TYPE(QMoments)                  :: A
    CHARACTER(LEN=*),OPTIONAL       :: Tag_O
    CALL Get(A%MaxL,   'MaxL')
    CALL Get(A%MaxLen, 'MaxLen')
    CALL Get(A%CenterX,'CenterX')
    CALL Get(A%CenterY,'CenterY')
    CALL Get(A%CenterZ,'CenterZ')
    CALL Get(A%QMat,   'QMat')
  END SUBROUTINE Get_QMoments
!========================================================================================
! Print  Sperical moments
!========================================================================================
  SUBROUTINE PPrint_SpMoments(A,Name,FileName_O,Unit_O)
    TYPE(SpMoments)                  :: A
    CHARACTER(LEN=*)                 :: Name   
    CHARACTER(LEN=*),OPTIONAL        :: FileName_O
    INTEGER,OPTIONAL                 :: Unit_O
    INTEGER                          :: OutU
    INTEGER                          :: L,M,LM
!
    IF(.NOT. AllocQ(A%Alloc)) THEN
       CALL Halt(' SpMoments not allocated in PPrint_SpMoments')
    ENDIF
    IF(PRESENT(Unit_O)) THEN
       OutU=Unit_O
    ELSE
       OutU=Out
    ENDIF
    IF(PRESENT(FileName_O) .AND. OutU /= 6) THEN
       CALL OpenASCII(FileName_O,OutU)
    ELSEIF(OutU /= 6) THEN
       CALL OpenASCII(OutFile,OutU)
    ENDIF
!
    WRITE(OutU,5)
    WRITE(OutU,10) Name
    WRITE(OutU,11) A%MaxL,A%MaxLM
    WRITE(OutU,12) A%CenterX,A%CenterY,A%CenterZ
    WRITE(OutU,6)
!
    DO L = 0,A%MaxL
       WRITE(OutU,20) L
       DO M = 0,L
          LM = LTD(L)+M
          WRITE(OutU,21) L,M,A%CMMat%D(LM),L,M,A%SMMat%D(LM)
       ENDDO
    ENDDO
    WRITE(OutU,5)
!
    RETURN
!
5   FORMAT(1x,'=========================================================================')
6   FORMAT(1x,'=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=')
10  FORMAT(1x,A)
11  FORMAT(1x,'MaxL = ',I2,' MaxLM = ',I3)
12  FORMAT(1x,'CellCenter = ('F20.16,','F20.16,',',F20.16,')')
20  FORMAT(1x,'L = ',I2)
21  FORMAT(2x,'CMat('I2,',',I2,') = ',F20.16,2x,'SMat('I2,',',I2,') = ',F20.16)
!
  END SUBROUTINE PPrint_SpMoments
!========================================================================================
! Print  Cartesian moments
!========================================================================================
  SUBROUTINE PPrint_QMoments(A,Name,FileName_O,Unit_O)
    TYPE(QMoments)                   :: A
    CHARACTER(LEN=*)                 :: Name   
    CHARACTER(LEN=*),OPTIONAL        :: FileName_O
    INTEGER,OPTIONAL                 :: Unit_O
    INTEGER                          :: OutU
    INTEGER                          :: L,M,N,LMN
!
    IF(.NOT. AllocQ(A%Alloc)) THEN
       CALL Halt(' QMoments not allocated in PPrint_QMoments')
    ENDIF
    IF(PRESENT(Unit_O)) THEN
       OutU=Unit_O
    ELSE
       OutU=Out
    ENDIF
    IF(PRESENT(FileName_O) .AND. OutU /= 6) THEN
       CALL OpenASCII(FileName_O,OutU)
    ELSEIF(OutU /= 6) THEN
       CALL OpenASCII(OutFile,OutU)
    ENDIF
!
    WRITE(OutU,5)
    WRITE(OutU,10) Name
    WRITE(OutU,11) A%MaxL,A%MaxLen
    WRITE(OutU,12) A%CenterX,A%CenterY,A%CenterZ
    WRITE(OutU,6)
!
    DO L=0,A%MaxL
       DO M=0,A%MaxL-L
          DO N=0,A%MaxL-L-M
             LMN=LMNDex(L,M,N)
             WRITE(OutU,21) L,M,N,A%QMat%D(LMN)
          ENDDO
       ENDDO
    ENDDO
    WRITE(OutU,5)
!
    RETURN
!
5   FORMAT(1x,'=========================================================================')
6   FORMAT(1x,'=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=')
10  FORMAT(1x,A)
11  FORMAT(1x,'MaxL = ',I2,' MaxLen = ',I3)
12  FORMAT(1x,'CellCenter = ('F20.16,','F20.16,',',F20.16,')')
21  FORMAT(2x,'QMat('I2,',',I2,',',I2,') = ',F20.16)
!
  END SUBROUTINE PPrint_QMoments
!========================================================================================
! FT_FSCriptC
!========================================================================================
  FUNCTION FT_FScriptC(L,ExpFac,R)
    INTEGER                    :: L,IFac,Isgn
    REAL(DOUBLE)               :: ExpFac,R,FT_FScriptC,Fac
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.00000000000000000D0 ,1.61199195401646964D0 , &
                    1.39399011673806825D0 ,1.24836057075206815D0 ,1.14180213615392840D0 , &
                    1.05927697236749128D0 ,9.92824332104422014D-1,9.37767832227969827D-1, &
                    8.91149652589504141D-1,8.50992321055495993D-1,8.15915401864279526D-1, &
                    7.84921233738737833D-1,7.57267966537096089D-1,7.32390674279586654D-1, &
                    7.09850372573553877D-1,6.89299957250757043D-1,6.70460791874183787D-1, &
                    6.53106213773379367D-1,6.37049661171763433D-1,6.22135962734544435D-1, &
                    6.08234838286273008D-1,5.95235975457893168D-1,5.83045248969098145D-1, &
                    5.71581781316674666D-1,5.60775631817133548D-1,5.50565960943842365D-1, &
                    5.40899558418448758D-1,5.31729652703953292D-1,5.23014940361368490D-1, &
                    5.14718788772630577D-1,5.06808576734283345D-1,4.99255145565447538D-1, &
                    4.92032339458264325D-1,4.85116618392624313D-1,4.78486730436807416D-1, & 
                    4.72123432945047398D-1,4.66009254246335764D-1,4.60128289044821961D-1, &
                    4.54466022030378857D-1,4.49009175209461305D-1,4.43745575272018022D-1 /) 
!
    IFac = 1+(-1)**L
    IF(IFac == 0) THEN
       FT_FScriptC = Zero
    ELSE
       Isgn = (-1)**(L/2)
       Fac  = ExpFac/DBLE(2*L-1)
       FT_FScriptC = Isgn*(Sfac(L)*R*EXP(-Fac*R*R))**(2*L-1)   
    ENDIF
! 
  END FUNCTION FT_FScriptC
!========================================================================================
! FT_FSCriptS
!========================================================================================
  FUNCTION FT_FScriptS(L,ExpFac,R)
    INTEGER                    :: L,IFac,Isgn
    REAL(DOUBLE)               :: ExpFac,R,FT_FScriptS,Fac
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.00000000000000000D0 ,1.61199195401646964D0 , &
                    1.39399011673806825D0 ,1.24836057075206815D0 ,1.14180213615392840D0 , &
                    1.05927697236749128D0 ,9.92824332104422014D-1,9.37767832227969827D-1, &
                    8.91149652589504141D-1,8.50992321055495993D-1,8.15915401864279526D-1, &
                    7.84921233738737833D-1,7.57267966537096089D-1,7.32390674279586654D-1, &
                    7.09850372573553877D-1,6.89299957250757043D-1,6.70460791874183787D-1, &
                    6.53106213773379367D-1,6.37049661171763433D-1,6.22135962734544435D-1, &
                    6.08234838286273008D-1,5.95235975457893168D-1,5.83045248969098145D-1, &
                    5.71581781316674666D-1,5.60775631817133548D-1,5.50565960943842365D-1, &
                    5.40899558418448758D-1,5.31729652703953292D-1,5.23014940361368490D-1, &
                    5.14718788772630577D-1,5.06808576734283345D-1,4.99255145565447538D-1, &
                    4.92032339458264325D-1,4.85116618392624313D-1,4.78486730436807416D-1, & 
                    4.72123432945047398D-1,4.66009254246335764D-1,4.60128289044821961D-1, &
                    4.54466022030378857D-1,4.49009175209461305D-1,4.43745575272018022D-1 /) 
!
    IFac = 1-(-1)**L
    IF(IFac == 0) THEN
       FT_FScriptS = Zero
    ELSE
       Isgn = (-1)**((L-1)/2)
       Fac  = ExpFac/DBLE(2*L-1)
       FT_FScriptS = Isgn*(Sfac(L)*R*EXP(-Fac*R*R))**(2*L-1)
    ENDIF
!
  END FUNCTION FT_FScriptS
!========================================================================================
! FT_FSCriptC
!========================================================================================
  FUNCTION GScript(L,R)
    INTEGER                    :: L,LS
    REAL(DOUBLE)               :: R,SqrtR,GScript,XSUM
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.0000000000000000D0 ,1.333333333333333333D0 ,&
                       7.3029674334022148D-1,5.3412580601946498D-1,4.2897258944050779D-1,&
                       3.6130260767122454D-1,3.1338173978619282D-1,2.7736675692301273D-1,&
                       2.4916970717917522D-1,2.2642030439700098D-1,2.0763707534036204D-1,&
                       1.9184081786310919D-1,1.7835560611619894D-1,1.6669836456773905D-1,&
                       1.5651386447776487D-1,1.4753458765116469D-1,1.3955494670757578D-1,&
                       1.3241416731800685D-1,1.2598459485542503D-1,1.2016350807025758D-1,&
                       1.1486726370711862D-1,1.1002702835105538D-1,1.0558561445605698D-1,&
                       1.0149509929347097D-1,9.7715008595692035D-2,9.4210913822889116D-2,&
                       9.0953336661706548D-2,8.7916884656620846D-2,8.5079562763772593D-2,&
                       8.2422220248006543D-2,7.9928102738676689D-2,7.7582486742601370D-2,&
                       7.5372379364895518D-2,7.3286270006162048D-2,7.1313923796257465D-2,&
                       6.9446208774383175D-2,6.7674950532187841D-2,6.5992809342867331D-2,&
                       6.4393175806982839D-2,6.2870081828999622D-2,6.1418124351705448D-2 /)
!
    SqrtR      = SQRT(R)
    IF(L == 0) THEN
       GScript = (One-ERF(SqrtR))
    ELSE
       XSUM = SFAC(1)*EXP(-R)
       DO LS = 2,L
          XSUM  = XSUM+(Sfac(LS)*R*EXP(-R/DBLE(LS-1)))**(LS-1)
       ENDDO 
       GScript = (One-ERF(SqrtR))+(SqrtR/SqrtPi)*XSUM
    ENDIF
!
  END FUNCTION GScript
!========================================================================================
! FT_FSCriptS
!========================================================================================
  FUNCTION FScript(L,R)
    INTEGER                    :: L,LS
    REAL(DOUBLE)               :: R,SqrtR,FScript,XSUM
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.0000000000000000D0 ,1.333333333333333333D0 ,&
                       7.3029674334022148D-1,5.3412580601946498D-1,4.2897258944050779D-1,&
                       3.6130260767122454D-1,3.1338173978619282D-1,2.7736675692301273D-1,&
                       2.4916970717917522D-1,2.2642030439700098D-1,2.0763707534036204D-1,&
                       1.9184081786310919D-1,1.7835560611619894D-1,1.6669836456773905D-1,&
                       1.5651386447776487D-1,1.4753458765116469D-1,1.3955494670757578D-1,&
                       1.3241416731800685D-1,1.2598459485542503D-1,1.2016350807025758D-1,&
                       1.1486726370711862D-1,1.1002702835105538D-1,1.0558561445605698D-1,&
                       1.0149509929347097D-1,9.7715008595692035D-2,9.4210913822889116D-2,&
                       9.0953336661706548D-2,8.7916884656620846D-2,8.5079562763772593D-2,&
                       8.2422220248006543D-2,7.9928102738676689D-2,7.7582486742601370D-2,&
                       7.5372379364895518D-2,7.3286270006162048D-2,7.1313923796257465D-2,&
                       6.9446208774383175D-2,6.7674950532187841D-2,6.5992809342867331D-2,&
                       6.4393175806982839D-2,6.2870081828999622D-2,6.1418124351705448D-2 /)
!
    SqrtR      = SQRT(R)
    IF(L == 0) THEN
       FScript = ERF(SqrtR)
    ELSE
       XSUM = SFAC(1)*EXP(-R)
       DO LS = 2,L
          XSUM  = XSUM+(Sfac(LS)*R*EXP(-R/DBLE(LS-1)))**(LS-1)
       ENDDO 
       FScript = ERF(SqrtR)-(SqrtR/SqrtPi)*XSUM
    ENDIF
!
  END FUNCTION FScript
!
END MODULE Moments

