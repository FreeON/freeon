MODULE PFFT
#ifdef PERIODIC  
  USE Derivedtypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE Parse
  USE InOut
  USE Macros
  USE Thresholding
  USE BoundingBox
  USE PoleTree
  USE SpecFun
  USE Globals 
  IMPLICIT NONE
!====================================================================================
! Globals
!====================================================================================
  INTEGER                           :: Dimen,Volume
  REAL(DOUBLE), DIMENSION(0:FFLen2) :: TensorC,TensorS
  CONTAINS
!========================================================================================
! Calculate the PFF
!========================================================================================
    SUBROUTINE PFFTensor(RMIN,Volume,MaxL_O)
      INTEGER,OPTIONAL    :: MaxL_O
      INTEGER             :: MaxL
      LOGICAL             :: HaveTensor
      REAL(DOUBLE)        :: RMIN,Volume 
!   
      IF(PRESENT(MaxL_O)) THEN
         MaxL = MaxL_O
      ELSE
         MaxL = FFEll2
      ENDIF
!
      TensorC = Zero
      TensorS = Zero
      IF(Dimen==0) RETURN
!
!      WRITE(*,*) 'GETTING TENSOR'
      CALL GetTensor(MaxL,HaveTensor)
!
      IF(.NOT. HaveTensor) THEN
!         WRITE(*,*) 'MAKING  TENSOR'
         CALL MakePFFT(RMIN,Volume,MaxL)
!         WRITE(*,*) 'SAVEING TENSOR'
         CALL PutTensor(MaxL)
      ENDIF
!
    END SUBROUTINE PFFTensor
!========================================================================================
! 
!========================================================================================
    SUBROUTINE GetTensor(MaxL,HaveTensor)
      INTEGER              :: MaxL,L,M,LM
      LOGICAL              :: HaveTensor
      CHARACTER(LEN=120)   :: FileName
      CHARACTER(LEN=1)     :: AWX,AWY,AWZ
!
      AWX = '0'
      AWY = '0'
      AWZ = '0'
      IF(GM%AutoW(1)) AWX = '1'
      IF(GM%AutoW(2)) AWY = '1'      
      IF(GM%AutoW(3)) AWZ = '1'
      FileName= TRIM(Args%C%C(1)) // "_Geom#" // TRIM(CurGeom) //    &
                "_AW"   // TRIM(AWX) // TRIM(AWY) // TRIM(AWZ) //    &
                "_NC"   // TRIM(IntToChar(CSMM1%NCells))   //        &    
                "_LM"   // TRIM(IntToChar(MaxL)) //                  &
                ".PFFT"
!
      INQUIRE(FILE=FileName,EXIST=HaveTensor)
      IF(HaveTensor) THEN
         OPEN(UNIT=77,FILE=FileName,FORM='UNFORMATTED',STATUS='OLD')
         DO L = 1,MaxL
            DO M = 0,L
               LM = LTD(L)+M
               READ(77) TensorC(LM),TensorS(LM)
            ENDDO
         ENDDO
         CLOSE(77)
      ENDIF
!
    END SUBROUTINE GetTensor
!========================================================================================
! 
!========================================================================================
    SUBROUTINE PutTensor(MaxL)
      INTEGER             :: MaxL,IMin,JMin,KMin,L,M,LM
      CHARACTER(LEN=120)  :: FileName
      CHARACTER(LEN=1)     :: AWX,AWY,AWZ
!
      AWX = '0'
      AWY = '0'
      AWZ = '0'
      IF(GM%AutoW(1)) AWX = '1'
      IF(GM%AutoW(2)) AWY = '1'      
      IF(GM%AutoW(3)) AWZ = '1'
      FileName= TRIM(Args%C%C(1)) // "_Geom#" // TRIM(CurGeom) //    &
                "_AW"   // TRIM(AWX) // TRIM(AWY) // TRIM(AWZ) //    &
                "_NC"   // TRIM(IntToChar(CSMM1%NCells))   //        &    
                "_LM"   // TRIM(IntToChar(MaxL)) //                  &
                ".PFFT"
!
      OPEN(UNIT=77,FILE=FileName,FORM='UNFORMATTED',STATUS='NEW')
      DO L = 1,MaxL
         DO M = 0,L
            LM = LTD(L)+M
            WRITE(77) TensorC(LM),TensorS(LM)
         ENDDO
      ENDDO
      CLOSE(77)
!
    END SUBROUTINE PutTensor
!========================================================================================
! Calculate the PFFTensor
!========================================================================================
    SUBROUTINE MakePFFT(RMIN,Volume,MaxL)
      INTEGER                           :: MaxL
      INTEGER                           :: I,J,K,L,M,LM,NC
      INTEGER                           :: LSwitch
      REAL(DOUBLE),DIMENSION(3)         :: PQ
      REAL(DOUBLE)                      :: CFac,SFac,BetaSq,Rad,RadSq,ExpFac
      REAL(DOUBLE)                      :: RMIN,RMAX,Accuracy,LenScale,Volume
!
!     Initialize and Zero dimension
!
      TensorC = Zero
      TensorS = Zero
      Accuracy = 1.D-16
      IF(Dimen==0) RETURN
!
!     One Dimension
!
      IF(Dimen==1) THEN
         IF(GM%AutoW(1)) THEN
            CALL IrRegular(MaxL,GM%BoxShape%D(1,1),Zero,Zero)
         ELSEIF(GM%AutoW(2)) THEN
            CALL IrRegular(MaxL,Zero,GM%BoxShape%D(2,2),Zero)
         ELSEIF(GM%AutoW(3)) THEN      
            CALL IrRegular(MaxL,Zero,Zero,GM%BoxShape%D(3,3))
         ENDIF
         NC = (CSMM1%NCells-1)/2
         DO L=1,MaxL
            DO M = 0,L
               LM = LTD(L)
               TensorC(LM) = Cpq(LM)*RZeta(L+1,NC)
            ENDDO
         ENDDO
      ENDIF
!
!     Two and Three Dimension
!
      IF(Dimen == 2 .OR. Dimen == 3) THEN
         LSwitch  = 10
         LenScale = Volume**(One/DBLE(Dimen))
         BetaSq   = (0.5D0/LenScale)**2
         RMAX = RMIN+LenScale*(One/Accuracy)**(One/DBLE(LSwitch))
!           
!        Sum the Real Space
!
         DO
            CALL New_CellSet_Sphere(CSMM2,GM%AutoW,GM%BoxShape%D,RMAX)
            IF(CSMM2%NCells .GT. 400000) THEN
               RMAX = 0.99*RMAX
               CALL Delete_CellSet(CSMM2)
            ELSE
               EXIT
            ENDIF
         ENDDO
!
         DO NC = 1,CSMM2%NCells
            PQ(:) = CSMM2%CellCarts%D(:,NC)
            RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
            IF(.NOT. InCell_CellSet(CSMM1,PQ(1),PQ(2),PQ(3))) THEN
               CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
               DO L = 1,MaxL
                  IF(L .LE. LSwitch) THEN
                     CFac = GScript(L,RadSq)
                  ELSE
                     CFac = One
                  ENDIF
                  DO M = 0,L
                     LM = LTD(L)+M
                     TensorC(LM)=TensorC(LM)+Cpq(LM)*CFac
                     TensorS(LM)=TensorS(LM)+Spq(LM)*CFac
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         CALL Delete_CellSet(CSMM2)
         IF(.TRUE.) RETURN
!
!        Sum the Reciprical Space 
!
         ExpFac = (Pi/BetaSq)**2
         RMAX = SQRT(ABS(LOG(Accuracy/(10.D0**(LSwitch)))/ExpFac))
         DO
            CALL New_CellSet_Sphere(CSMM2,GM%AutoW,GM%InvBoxSh%D,RMAX)
            IF(CSMM2%NCells .LT.27) THEN
               RMAX = 1.01D0*RMAX
            ELSE
               EXIT
            ENDIF
         ENDDO
!
         DO NC = 1,CSMM2%NCells
            PQ(:) = CSMM2%CellCarts%D(:,NC)
            Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
            IF(Rad .GT. 1.D-14) THEN
               CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
               DO L = 1,MaxL
                  IF(L .LE. LSwitch) THEN        
                     CFac = FT_FScriptC(L,ExpFac,Rad)/Volume
                     SFac = FT_FScriptS(L,ExpFac,Rad)/Volume
                  ELSE
                     CFac = Zero
                     SFac = Zero
                  ENDIF
                  DO M = 0,L
                     LM = LTD(L)+M
                     TensorC(LM)=TensorC(LM)+Cpq(LM)*CFac-Spq(LM)*SFac
                     TensorS(LM)=TensorS(LM)+Spq(LM)*CFac+Cpq(LM)*SFac
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         CALL Delete_CellSet(CSMM2)
!
!        Substract the inner boxes
!
         DO NC = 1,CSMM1%NCells
            PQ(:) = CSMM1%CellCarts%D(:,NC)
            RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
            IF(RadSq .GT. 1.D-14) THEN
               CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
               DO L = 1,MaxL
                  IF(L .LE. LSwitch) THEN
                     CFac = FScript(L,RadSq)
                  ELSE
                     CFac = Zero
                  ENDIF
                  DO M = 0,L
                     LM = LTD(L)+M
                     TensorC(LM)=TensorC(LM)-Cpq(LM)*CFac
                     TensorS(LM)=TensorS(LM)-Spq(LM)*CFac
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         CALL Delete_CellSet(CSMM2)
      ENDIF
!
!     Filter out the Zero Elements
!
      DO L = 0,MaxL
         DO M = 0,L
            LM = LTD(L)+M
            IF(ABS(TensorC(LM)) .LT. 1.D-15) TensorC(LM) = Zero
            IF(ABS(TensorS(LM)) .LT. 1.D-15) TensorS(LM) = Zero
         ENDDO
      ENDDO
!
    END SUBROUTINE MakePFFT
!========================================================================================
!   FT_FSCriptC
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
!   FT_FSCriptC
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
!   FT_FSCriptS
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
!========================================================================================
! Rieman Zeta Function
!========================================================================================
    FUNCTION RZeta(N,M)
      INTEGER                    :: N,M,I
      REAL(DOUBLE)               :: RZeta,RSum
      REAL(DOUBLE),DIMENSION(56) :: RZ = (/ 0.00000000000000000D0, 1.64493406684822644D0, & 
                   1.20205690315959429D0, 1.08232323371113819D0, 1.03692775514336993D0, & 
                   1.01734306198444914D0, 1.00834927738192283D0, 1.00407735619794434D0, &
                   1.00200839282608221D0, 1.00099457512781809D0, 1.00049418860411946D0, &
                   1.00024608655330805D0, 1.00012271334757849D0, 1.00006124813505870D0, &
                   1.00003058823630702D0, 1.00001528225940865D0, 1.00000763719763790D0, &
                   1.00000381729326500D0, 1.00000190821271655D0, 1.00000095396203387D0, &
                   1.00000047693298679D0, 1.00000023845050273D0, 1.00000011921992597D0, &
                   1.00000005960818905D0, 1.00000002980350351D0, 1.00000001490155483D0, &
                   1.00000000745071179D0, 1.00000000372533402D0, 1.00000000186265972D0, &
                   1.00000000093132743D0, 1.00000000046566291D0, 1.00000000023283118D0, &
                   1.00000000011641550D0, 1.00000000005820772D0, 1.00000000002910385D0, &
                   1.00000000001455192D0, 1.00000000000727596D0, 1.00000000000363798D0, &
                   1.00000000000181899D0, 1.00000000000090949D0, 1.00000000000045475D0, &
                   1.00000000000022737D0, 1.00000000000011369D0, 1.00000000000005684D0, & 
                   1.00000000000002842D0, 1.00000000000001421D0, 1.00000000000000711D0, &
                   1.00000000000000355D0, 1.00000000000000178D0, 1.00000000000000089D0, &
                   1.00000000000000044D0, 1.00000000000000022D0, 1.00000000000000011D0, &
                   1.00000000000000006D0, 1.00000000000000003D0, 1.00000000000000001D0  /)

!
      IF(N .LE. 56) THEN
         RZeta = RZ(N)
      ELSE
         RZeta = One
      ENDIF
!
      DO I=2,M
         RSum = RSum + One/(DBLE(I)**N)
      ENDDO
      RZeta = RZ(N) - RSum
!
    END FUNCTION RZeta
#endif
  END MODULE PFFT

