MODULE DYNAMO
#ifdef MMech
!
USE DerivedTypes
USE GlobalScalars
Use MemMan
USE InOut
USE Mechanics !Macros, ONLY : HasQM,HasMM,MMOnly 
USE ProcessControl
USE IntCoo
!
   USE definitions
   USE constants
   USE files
   USE io_units
   USE status
   USE sequence
   USE symmetry
   USE elements
   USE parsing
   USE strings
   USE pdb_io 
   USE mm_system
   USE connectivity 
   USE sort 
   USE mm_file_data
   USE mm_file_io 
   USE mm_terms
   USE ATOMS
   USE coordinate_io
   USE constraint   
   USE energy_covalent
   USE energy_non_bonding
IMPLICIT NONE
TYPE(CRDS) :: GM_MM
!
CONTAINS
!
!--------------------------------------------------------
!
   SUBROUTINE EXCL(Natoms,Cur,InfFile,E_LJ_EXCL,E_C_EXCL,Grad_Loc)
!
! Calculate the Coulomb and Lennard-Jones exclusion 
! energies and gradients for MM and QMMM calculations
! MM_Natoms and MM geometry must be used
! Presently, energies are given in KJ/mol, 
! gradients in KJ/mol/angstroem
!
     IMPLICIT NONE
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: InfFile
     REAL(DOUBLE),OPTIONAL :: E_LJ_EXCL,E_C_EXCL
     TYPE(INT_VECT) :: AtmMark
     TYPE(CRDS) :: GM_Loc 
     INTEGER :: I,J,K,L,M,N,Natoms
     INTEGER :: NMAX12,NMAX13,NMAX14,NMax_Excl,NMax_Excl14
     TYPE(DBL_VECT)    :: LJEps14,LJEps,LJRad,Charge14
     TYPE(INT_RNK2)    :: Top_Loc
     CHARACTER(LEN=3)  :: Cur
     REAL(DOUBLE) :: XI,YI,ZI,XJ,YJ,ZJ,QI,QJ,CONVF,CONVF2,QI14,QJ14
     REAL(DOUBLE) :: RI,RJ,R6,DIJ,DIJ2,DIJAu,DIJ2Au,Pref1,Pref2
     REAL(DOUBLE) :: QQI,QQI14,QQJ,QQJ14
     TYPE(DBL_RNK2),OPTIONAL :: Grad_Loc
     TYPE(DBL_VECT) :: DVect     
     INTEGER :: ITop
!
! If QMMM, distinguish QM and MM atoms
!   
!    CALL OpenHDF(InfFile)
!
     GM_Loc%NAtms=Natoms
     CALL New(GM_Loc)
     CALL New(DVect,3)
     CALL New(AtmMark,Natoms)
!
     IF(HasQM()) THEN
       CALL Get(AtmMark,'AtmMark')
     ELSE 
       AtmMark%I(:)=0 !!! all MM atoms
     ENDIF
!
! Get MM geometry and LJ data
! also get 14-charges
! also get Topology matrices for 1-2, 1-3 and 1-4 neighbour atoms
!
     CALL Get(GM_Loc%Carts,'cartesians',Tag_O='GM_MM'//Cur)
     IF(PRESENT(E_LJ_EXCL)) THEN
       CALL New(LJEps,Natoms)
       CALL New(LJEps14,Natoms)
       CALL New(LJRad,Natoms)
       CALL Get(LJRad,'LJRad')
       CALL Get(LJEps,'LJEps')
       CALL Get(LJEps14,'LJEps14')
       E_LJ_EXCL=Zero
       DO I=1,Natoms !!!later rather use scale factors
         IF(DABS(LJEps%D(I))<1.D-5) LJEps14%D(I)=Zero
       ENDDO
     ENDIF
     IF(PRESENT(E_C_EXCL)) THEN
       CALL Get(GM_Loc%AtNum,'atomicnumbers',Tag_O='GM_MM'//Cur)
       CALL New(CHARGE14,Natoms)
       CALL Get(Charge14,'Charge14')
       E_C_EXCL=Zero
       DO I=1,Natoms !!!later rather use scale factors
         IF(DABS(GM_Loc%AtNum%D(I))<1.D-5) Charge14%D(I)=Zero
       ENDDO
     ENDIF
!
! Now calculate Lennard-Jones and Coulomb exclusion energies
! Please note that all MM charges on the QM positions
! and on the link atoms, as well as on the atoms of the 
! broken covalent bonds should be set to zero at this point.
!
     CONVF=e2PerAngstroemToKJPerMol
!
! 12, 13 and 14 topologies merged
!
     DO ITop=1,2
!
     IF(ITop==1) THEN
       CALL Get(NMAX_EXCL,'NMAX_EXCL')
       CALL New(Top_Loc,(/Natoms,NMAX_EXCL+1/)) 
       CALL Get(Top_Loc,'TOP_EXCL')
     ENDIF
!
     IF(ITop==2) THEN
       CALL Get(NMax_Excl14,'NMAX_EXCL14')
       CALL New(Top_Loc,(/Natoms,NMax_Excl14+1/)) 
       CALL Get(Top_Loc,'TOP_EXCL14')
     ENDIF
!
     DO I=1,Natoms
       IF(PRESENT(E_LJ_EXCL)) THEN
         QI=LJEps%D(I)
         QI14=LJEps14%D(I)
         RI=LJRad%D(I)
       ENDIF
       IF(PRESENT(E_C_EXCL)) THEN
         QQI=GM_Loc%AtNum%D(I)
         QQI14=CHARGE14%D(I)
       ENDIF
     DO M=1,Top_Loc%I(I,1)
       J=Top_Loc%I(I,M+1)
     IF(J<=I) CYCLE !!!avoid double counting
     IF(AtmMark%I(I)/=0.AND.AtmMark%I(J)/=0) CYCLE !!! no QM-QM intr.act
       DVect%D(:)=GM_Loc%Carts%D(:,I)-GM_Loc%Carts%D(:,J)
       DVect%D(:)=DVect%D(:)/AngstromsToAU !!!work in angstroems
       DIJ2=DOT_PRODUCT(DVect%D,DVect%D)
       DIJ=SQRT(DIJ2)
       DIJAu=DIJ*AngstromsToAU
       DIJ2Au=DIJAu*DIJAu
       IF(PRESENT(E_LJ_EXCL)) THEN
         QJ=LJEps%D(J)
         QJ14=LJEps14%D(J)
         RJ=LJRad%D(J)
       ENDIF
       IF(PRESENT(E_C_EXCL)) THEN
         QQJ=GM_Loc%AtNum%D(J)
         QQJ14=CHARGE14%D(J)
       ENDIF
         IF(DIJ>0.001D0) THEN
           IF(PRESENT(E_LJ_EXCL)) THEN
             R6=(RI*RJ/DIJ)**6
             IF(ITop==2) THEN
               Pref1=(-QI14*QJ14)*R6
             ELSE
               Pref1=QI*QJ*R6
             ENDIF
               E_LJ_EXCL=E_LJ_EXCL+Pref1*(R6-One)
               IF(PRESENT(Grad_Loc)) THEN
             Pref2=Six*Pref1*(One-Two*R6)/DIJ2 
!since it is about exclusion gradients, they get a negative sign
             Grad_Loc%D(:,I)=Grad_Loc%D(:,I)+(-Pref2*DVect%D(:))
             Grad_Loc%D(:,J)=Grad_Loc%D(:,J)-(-Pref2*DVect%D(:))
               ENDIF
           ENDIF
           IF(PRESENT(E_C_EXCL)) THEN
             IF(ITop==2) THEN
               Pref1=(-QQI14*QQJ14)/DIJ*CONVF
             ELSE
               Pref1=QQI*QQJ/DIJ*CONVF
             ENDIF
               E_C_EXCL=E_C_EXCL+Pref1
               IF(PRESENT(Grad_Loc)) THEN
             Pref2=-Pref1/DIJ2 
!since it is about exclusion gradients, they get a negative sign
             Grad_Loc%D(:,I)=Grad_Loc%D(:,I)+(-Pref2*DVect%D(:))
             Grad_Loc%D(:,J)=Grad_Loc%D(:,J)-(-Pref2*DVect%D(:))
               ENDIF
           ENDIF
         ELSE
       write(*,*) 'd dij= ',i,j,dij,' x1= ',GM_Loc%Carts%D(1:3,i),GM_Loc%Carts%D(1:3,j) 
  Call MondoHalt(QMMM_ERROR,'Atoms are too close to each other in MM set')
         ENDIF
     ENDDO
     ENDDO
!
     CALL Delete(Top_Loc)
     ENDDO ! ITop
!
     CALL Delete(GM_Loc)
     CALL Delete(DVect)
     CALL Delete(AtmMark)
     IF(PRESENT(E_LJ_EXCL)) THEN
       CALL Delete(LJEps)
       CALL Delete(LJEps14)
       CALL Delete(LJRad)
       write(*,*) 'E_LJ_EXCL= ',E_LJ_EXCL
       CALL Put(E_LJ_EXCL,'E_LJ_EXCL')
     ENDIF
     IF(PRESENT(E_C_EXCL)) THEN
       CALL Delete(Charge14)
       write(*,*) 'E_C_EXCL= ',E_C_EXCL
       CALL Put(E_C_EXCL,'E_C_EXCL')
     ENDIF
!
   END SUBROUTINE EXCL
!-------------------------------------------------------------- 
!
   SUBROUTINE ENERGY_LENNARD_JONES(ELJ,ISet,BoxSize,Grad_Loc)
!
   IMPLICIT NONE
   INTEGER :: NBox,Natoms,NX,NY,NZ,I,J
   CHARACTER(LEN=3)   :: ISet
   REAL(DOUBLE) :: BXMIN,BYMIN,BZMIN,BoxSize
   INTEGER :: IX,IY,IZ,IOrd,IOrdD
   INTEGER :: I1,I2,JJ1,JJ2,IXD,IYD,IZD
   REAL(DOUBLE) :: X1,Y1,Z1,X2,Y2,Z2,R1,R2,R12,R6,Q1,Q2
   TYPE(INT_VECT) :: BoxI,BoxJ,AtmMark
   TYPE(CRDS) :: GM_Loc
   TYPE(DBL_VECT) :: LJEps,LJRad,DVect
   TYPE(DBL_RNK2),OPTIONAL :: Grad_Loc
   REAL(DOUBLE) :: ELJ,r12_2,Pref1,Pref2
!
! BoxSize is in angstroems
!
! get MM or QMMM atom coordinates as stored in GM_Loc
! get AtmMark to distinguish between QM and MM atoms
! get LJEps%D : LJ epsilon
! get LJRad%D : LJ radii of atoms
! calculate and get BoxI and BoxJ for a certain Box size, then
! sum up LJ contributions within Boxes and between 
! neighbouring Boxes
!
   CALL Get(GM_Loc,'GM_MM'//ISet)
   Natoms=GM_Loc%Natms
!
     CALL New(AtmMark,Natoms)
   IF(HasQM()) THEN
     CALL Get(AtmMark,'AtmMark')
   ELSE 
     AtmMark%I(:)=0 !!! all MM atoms
   ENDIF
!
   CALL New(DVect,3)
   CALL New(LJEps,Natoms)
   CALL New(LJRad,Natoms)
   CALL Get(LJEps,'LJEps')
   CALL Get(LJRad,'LJRad') 
!
   GM_Loc%Carts%D(:,:)=GM_Loc%Carts%D(:,:)/AngstromsToAU !!! work with angstroems here
!
! Now calculate distributions of atoms in the Boxes
!
   CALL SORT_INTO_Box1(BoxSize,GM_Loc%Carts%D,Natoms,NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
   NBox=NX*NY*NZ
   CALL New(BoxI,NBox+1)
   CALL New(BoxJ,Natoms)
!
   CALL SORT_INTO_Box2(BoxSize,GM_Loc%Carts%D,Natoms,NX,NY,NZ,BXMIN,BYMIN,BZMIN,BoxI,BoxJ)
!
! Calculate LJ energy
!
   ELJ=Zero
!
   DO IZ=1,NZ   
   DO IX=1,NX   
   DO IY=1,NY   
!
! indices of central and neighboring Boxes
!
       IOrd=NX*NY*(IZ-1)+NY*(IX-1)+IY
!
     DO I1=BoxI%I(IOrd),BoxI%I(IOrd+1)-1
       JJ1=BoxJ%I(I1)
        R1=LJRad%D(JJ1)
        Q1=LJEps%D(JJ1)
! second atom may come from central or neigbouring Boxes 
! and must be an MM atom, LJ is not calculated for QM-QM pairs
     DO IZD=-1,1
       IF(IZ+IZD>0 .AND. IZ+IZD<=NZ) THEN
     DO IXD=-1,1
       IF(IX+IXD>0 .AND. IX+IXD<=NX) THEN
     DO IYD=-1,1
       IF(IY+IYD>0 .AND. IY+IYD<=NY) THEN
!
       IOrdD=NX*NY*(IZ-1+IZD)+NY*(IX-1+IXD)+IY+IYD
         DO I2=BoxI%I(IOrdD),BoxI%I(IOrdD+1)-1
           JJ2=BoxJ%I(I2)
       IF(JJ2<=JJ1) CYCLE !!!avoid double counting
       IF(AtmMark%I(JJ1)/=0.AND.AtmMark%I(JJ2)/=0) CYCLE 
           R2=LJRad%D(JJ2)
           Q2=LJEps%D(JJ2)
           DVect%D(:)=GM_Loc%Carts%D(:,JJ1)-GM_Loc%Carts%D(:,JJ2)
           r12_2=DOT_PRODUCT(DVect%D,DVect%D)                
           R12=SQRT(r12_2)
           IF(R12>0.001D0) THEN
! calculate energy contribution
             R6=(R1*R2/R12)**6
             Pref1=R6*Q1*Q2
             ELJ=ELJ+Pref1*(R6-One)
! calculate gradient contribution
             IF(PRESENT(Grad_Loc)) THEN
           Pref2=Six*PREF1*(One-Two*R6)/r12_2
           Grad_Loc%D(:,JJ1)=Grad_Loc%D(:,JJ1)+Pref2*DVect%D(:)
           Grad_Loc%D(:,JJ2)=Grad_Loc%D(:,JJ2)-Pref2*DVect%D(:)
             ENDIF
           ENDIF
         ENDDO
!
       ENDIF
     ENDDO
       ENDIF
     ENDDO
       ENDIF
     ENDDO
!
     ENDDO
!
   ENDDO
   ENDDO
   ENDDO
!
!  CALL CloseHDF()
!
   CALL Delete(GM_Loc)
   CALL Delete(AtmMark)
   CALL Delete(DVect)
   CALL Delete(LJEps)
   CALL Delete(LJRad)
   CALL Delete(BoxI)
   CALL Delete(BoxJ)
!
   END SUBROUTINE ENERGY_LENNARD_JONES
!
!------------------------------------------------------
!
   SUBROUTINE Bond_Energy(EBOND,XYZ,Grad_Loc)  
!
   REAL(DOUBLE)  :: EBOND
   TYPE(DBL_RNK2),OPTIONAL :: Grad_Loc
   INTEGER            :: I,J,IB,NBond
   LOGICAL            :: CalcGrad
   REAL(DOUBLE)       :: D_Force,D_Bond,RIJ
   TYPE(DBL_Vect)     :: DVect,BondFC,BondEQ
   TYPE(INT_VECT) :: Active_Bond
   TYPE(INT_RNK2) :: BondIJ
   REAL(DOUBLE),DIMENSION(:,:) :: XYZ !!! must be in Angstroem !!!
!
   EBOND = Zero  
   CalcGrad=PRESENT(Grad_Loc)
!
   CALL Get(NBond,'MM_NBond')
   IF(NBond==0) RETURN
!
   CALL New(BondFC,NBond)
   CALL New(BondEQ,NBond)
   CALL New(BondIJ,(/2,NBond/))
   CALL New(Active_Bond,NBond)
   CALL Get(BondFC,'BondFC')
   CALL Get(BondEQ,'BondEQ')
   CALL Get(BondIJ,'MM_BondIJ')
   CALL Get(Active_Bond,'Active_Bond')
!
   CALL New(DVect,3)
!
   DO IB = 1,NBond
     IF(Active_Bond%I(IB)/=1) CYCLE
!
      I = BondIJ%I(1,IB)
      J = BondIJ%I(2,IB)
      DVect%D = XYZ(1:3,I)-XYZ(1:3,J)
      RIJ  = SQRT(DOT_PRODUCT(DVect%D,DVect%D))
      D_Bond = RIJ-BondEQ%D(IB)
      D_Force=BondFC%D(IB)*D_Bond
!
! bond energy
!
      EBOND=EBOND+(D_Force*D_Bond)
!
      IF(.NOT.CalcGrad ) CYCLE
!
! bond gradient
!
      D_Force=(Two*D_Force)/RIJ
      Grad_Loc%D(1:3,I)=Grad_Loc%D(1:3,I)+(D_Force*DVect%D)
      Grad_Loc%D(1:3,J)=Grad_Loc%D(1:3,J)-(D_Force*DVect%D)
!
   END DO
!
   CALL Delete(Active_Bond)
   CALL Delete(DVect)
   CALL Delete(BondFC)
   CALL Delete(BondEQ)
   CALL Delete(BondIJ)
!
   END SUBROUTINE Bond_Energy
!
!-------------------------------------------------------------
!
   SUBROUTINE Angle_Energy(EANGLE,XYZ,Grad_Loc)
!
   REAL(DOUBLE) :: EANGLE
   TYPE(DBL_RNK2),OPTIONAL :: Grad_Loc
   TYPE(INT_RNK2) :: AngleIJK
   INTEGER        :: NAngle,IA,I,J,K
   LOGICAL        :: CalcGrad
   REAL(DOUBLE) :: D_Force,D_Angle,Cos_AngleIJK,RIJ,RJK,VAngleIJK
   TYPE(DBL_VECT) :: DVectIJ,DVectJK,FCI,FCJ,FCK
   REAL(DOUBLE), PARAMETER :: DOT_LIMIT = 0.999999D0
   TYPE(INT_VECT) :: ACTIVE_ANGLE
   REAL(DOUBLE),DIMENSION(:,:) :: XYZ !!! must be in Angstroem !!!
!
   EANGLE = Zero   
   CalcGrad=PRESENT(Grad_Loc)
!
   CALL Get(NAngle,'MM_NAngle')
!
   IF(NAngle==0) RETURN
!
   CALL New(ACTIVE_ANGLE,NAngle)
   CALL New(AngleIJK,(/3,NAngle/))
   CALL New(DVectIJ,3)
   CALL New(DVectJK,3)
   CALL New(FCI,3)
   CALL New(FCJ,3)
   CALL New(FCK,3)
   CALL Get(ACTIVE_ANGLE,'ACTIVE_ANGLE')
   CALL Get(AngleIJK,'MM_AngleIJK')
!
   DO IA = 1,NAngle
!
     IF(ACTIVE_ANGLE%I(IA)/=1) CYCLE
!
      I = AngleIJK%I(1,IA)
      J = AngleIJK%I(2,IA)
      K = AngleIJK%I(3,IA)
!
      DVectIJ%D=XYZ(1:3,I)-XYZ(1:3,J)
      DVectJK%D=XYZ(1:3,K)-XYZ(1:3,J)
!
      RIJ=SQRT(DOT_PRODUCT(DVectIJ%D,DVectIJ%D))
      RJK=SQRT(DOT_PRODUCT(DVectJK%D,DVectJK%D))
!
      DVectIJ%D=DVectIJ%D/RIJ
      DVectJK%D=DVectJK%D/RJK
!
      Cos_AngleIJK=DOT_PRODUCT(DVectIJ%D,DVectJK%D)
!
! accuracy of cos(angle) controlled for ACOS function
!
      IF(DABS(Cos_AngleIJK)>0.999999D0) THEN
        Cos_AngleIJK=SIGN(DABS(Cos_AngleIJK),Cos_AngleIJK)
      ENDIF
!
      VAngleIJK=ACOS(Cos_AngleIJK)
!
      D_Angle=VAngleIJK-ANGLES(IA)%EQ
      D_Force=ANGLES(IA)%FC*D_Angle
!
! contributions to energy
!
      EANGLE=EANGLE+(D_Force*D_Angle)
!
      IF(.NOT.CalcGrad) CYCLE
!
      D_Force=-Two*D_Force/SQRT(One-Cos_AngleIJK*Cos_AngleIJK)
!
      FCI%D=(DVectJK%D-Cos_AngleIJK*DVectIJ%D)/RIJ
      FCK%D=(DVectIJ%D-Cos_AngleIJK*DVectJK%D)/RJK
      FCJ%D=-(FCI%D+FCK%D)

      Grad_Loc%D(1:3,I)=Grad_Loc%D(1:3,I)+D_Force*FCI%D
      Grad_Loc%D(1:3,J)=Grad_Loc%D(1:3,J)+D_Force*FCJ%D
      Grad_Loc%D(1:3,K)=Grad_Loc%D(1:3,K)+D_Force*FCK%D
!
   END DO
!
   CALL Delete(ACTIVE_ANGLE)
   CALL Delete(AngleIJK)
   CALL Delete(DVectIJ)
   CALL Delete(DVectJK)
   CALL Delete(FCI)
   CALL Delete(FCJ)
   CALL Delete(FCK)
!
   END SUBROUTINE Angle_Energy
!
!--------------------------------------------------------------------
!
   SUBROUTINE Torsion_Energy(ETorsion,XYZ,Grad_Loc)
!
   REAL(DOUBLE) :: ETorsion
   TYPE(DBL_RNK2),OPTIONAL :: Grad_Loc
   TYPE(INT_VECT) :: Active_Torsion
   TYPE(DBL_VECT) :: TorsionEQ,TorsionFC
   TYPE(INT_VECT) :: TorsionPeriod
   TYPE(INT_RNK2) :: TorsionIJKL
   INTEGER :: NTorsion
   REAL(DOUBLE),DIMENSION(:,:) :: XYZ !!! must be in Angstroem !!!
!
   ETorsion=Zero
   CALL Get(NTorsion,'MM_NTorsion')
!
   IF(NTorsion==0) RETURN 
!
   CALL New(TorsionIJKL,(/4,NTorsion/))
   CALL New(Active_Torsion,NTorsion)
   CALL New(TorsionEQ,NTorsion)
   CALL New(TorsionFC,NTorsion)
   CALL New(TorsionPeriod,NTorsion)
   CALL Get(TorsionIJKL,'MM_TorsionIJKL')
   CALL Get(Active_Torsion,'Active_Torsion')
   CALL Get(TorsionEQ,'TorsionEQ')
   CALL Get(TorsionFC,'TorsionFC')
   CALL Get(TorsionPeriod,'TorsionPeriod')
!
   CALL TorsionalEnergy(NTorsion,TorsionIJKL,ETorsion,XYZ, &
       Grad_Loc,Active_Torsion,TorsionEQ,TorsionFC,TorsionPeriod)
!
   CALL Delete(Active_Torsion)
   CALL Delete(TorsionIJKL)
   CALL Delete(TorsionEQ)
   CALL Delete(TorsionFC)
   CALL Delete(TorsionPeriod)
!
   END SUBROUTINE Torsion_Energy
!
!--------------------------------------------------------------------
!
   SUBROUTINE OutOfPlane_Energy(EOutOfPlane,XYZ,Grad_Loc)
!
   REAL(DOUBLE) :: EOutOfPlane
   TYPE(DBL_RNK2),OPTIONAL :: Grad_Loc
   TYPE(INT_VECT) :: Active_OutOfPlane
   TYPE(INT_RNK2) :: OutOfPlaneIJKL
   TYPE(DBL_VECT) :: OutOfPlaneEQ,OutOfPlaneFC
   TYPE(INT_VECT) :: OutOfPlanePeriod
   INTEGER :: NOutOfPlane
   REAL(DOUBLE),DIMENSION(:,:) :: XYZ !!! must be in Angstroem !!!
!
   EOutOfPlane=Zero
   CALL Get(NOutOfPlane,'MM_NOutOfPlane')
!
   IF(NOutOfPlane==0) RETURN 
!
   CALL New(OutOfPlaneIJKL,(/4,NOutOfPlane/))
   CALL New(Active_OutOfPlane,NOutOfPlane)
   CALL New(OutOfPlaneEQ,NOutOfPlane)
   CALL New(OutOfPlaneFC,NOutOfPlane)
   CALL New(OutOfPlanePeriod,NOutOfPlane)
   CALL Get(OutOfPlaneIJKL,'MM_OutOfPlaneIJKL')
   CALL Get(Active_OutOfPlane,'Active_OutOfPlane')
   CALL Get(OutOfPlaneEQ,'OutOfPlaneEQ')
   CALL Get(OutOfPlaneFC,'OutOfPlaneFC')
   CALL Get(OutOfPlanePeriod,'OutOfPlanePeriod')
!
   CALL TorsionalEnergy(NOutOfPlane,OutOfPlaneIJKL,EOutOfPlane,XYZ, &
        Grad_Loc,Active_OutOfPlane,OutOfPlaneEQ,OutOfPlaneFC,OutOfPlanePeriod)
!
   CALL Delete(Active_OutOfPlane)
   CALL Delete(OutOfPlaneIJKL)
   CALL Delete(OutOfPlaneEQ)
   CALL Delete(OutOfPlaneFC)
   CALL Delete(OutOfPlanePeriod)
!
   END SUBROUTINE OutOfPlane_Energy
!
!-------------------------------------------------------------
!
   SUBROUTINE TorsionalEnergy(NCoord,CoordIJKL,ECoord,XYZ,Grad_Loc,Active_Coord,CoordEQ,CoordFC,CoordPeriod)
!
   IMPLICIT NONE
   REAL(DOUBLE) :: Energy_Coord
   TYPE(DBL_RNK2),OPTIONAL :: Grad_Loc
   TYPE(INT_VECT) :: Active_Coord
   TYPE(DBL_VECT) :: CoordEQ,CoordFC
   TYPE(INT_VECT) :: CoordPeriod
   TYPE(INT_RNK2) :: CoordIJKL
   INTEGER :: NCoord     
   REAL(DOUBLE),DIMENSION(:,:) :: XYZ !!! must be in Angstroem !!!
   INTEGER            :: ICoord,I,IFAC,J,JFAC,K,KFAC,L,LFAC
   LOGICAL            :: CalcGrad
   REAL(DOUBLE) :: CosNPhi,CosPhi,CosPhi2,DCos,DForce
   REAL(DOUBLE) :: IJCOSJK,KLCOSJK,D2Cos,Fact1,Fact2,NFJKL,NFIJK
   REAL(DOUBLE) :: ECoord,DJK
   TYPE(DBL_VECT) :: FIJK,DVectIJ,DVectJK,DVectKL,FJKL
   TYPE(DBL_VECT) :: DGradI,DGradJ,DGradK,DGradL
!
   ECoord=Zero   
   CalcGrad=PRESENT(Grad_Loc)
!
   CALL New(DVectIJ,3)
   CALL New(DVectJK,3)
   CALL New(DVectKL,3)
   CALL New(FIJK,3)
   CALL New(FJKL,3)
!
   CALL New(DGradI,3)
   CALL New(DGradJ,3)
   CALL New(DGradK,3)
   CALL New(DGradL,3)
!
   DO ICoord = 1,NCoord
   IF(Active_Coord%I(ICoord)/=1) CYCLE
!
      I=CoordIJKL%I(1,ICoord)
      J=CoordIJKL%I(2,ICoord)
      K=CoordIJKL%I(3,ICoord)
      L=CoordIJKL%I(4,ICoord)
!
      DVectIJ%D=XYZ(1:3,I)-XYZ(1:3,J)
      DVectJK%D=XYZ(1:3,K)-XYZ(1:3,J)
      DVectKL%D=XYZ(1:3,L)-XYZ(1:3,K)
!
      DJK=SQRT(DOT_PRODUCT(DVectJK%D,DVectJK%D))
!
      DVectJK%D=DVectJK%D/DJK
!
      IJCOSJK=DOT_PRODUCT(DVectIJ%D,DVectJK%D)
      KLCOSJK=DOT_PRODUCT(DVectKL%D,DVectJK%D)
!
      FIJK%D=DVectIJ%D-IJCOSJK*DVectJK%D
      FJKL%D=DVectKL%D-KLCOSJK*DVectJK%D
!
      NFIJK=SQRT(DOT_PRODUCT(FIJK%D,FIJK%D))
      NFJKL=SQRT(DOT_PRODUCT(FJKL%D,FJKL%D))
!
      FIJK%D=FIJK%D/NFIJK
      FJKL%D=FJKL%D/NFJKL
!
      CosPhi=DOT_PRODUCT(FIJK%D,FJKL%D)
!
      CosPhi=SIGN(MIN(ABS(CosPhi),1.D0),CosPhi)
      CosPhi2=CosPhi*CosPhi
!
! Observe periodicity depending on hybridization of J and K
! Use Dynamo's convention
!
      SELECT CASE(CoordPeriod%I(ICoord))
      CASE(0) ; CosNPhi = One   
                DCos    = Zero  
                D2Cos   = Zero  
      CASE(1) ; CosNPhi = CosPhi
                DCos    = One    
                D2Cos   = Zero  
      CASE(2) ; CosNPhi = Two*CosPhi2-One
                DCos    = Four*CosPhi
                D2Cos   = Four  
      CASE(3) ; CosNPhi = (Four*CosPhi2-Three)*CosPhi
                DCos    = 12.D0*CosPhi2-Three
                D2Cos   = 24.D0*CosPhi
      CASE(4) ; CosNPhi =  8.D0*(CosPhi2-One)*CosPhi2+One   
                DCos    = 16.D0*(Two*CosPhi2-One)*CosPhi
                D2Cos   = 16.D0*(Six*CosPhi2-One)
      CASE(5) ; CosNPhi = ((16.D0*CosPhi2-20.D0)*CosPhi2+5.D0)*CosPhi
                DCos    = 20.D0*(Four*CosPhi2-Three)*CosPhi2+Five
                D2Cos   = 40.D0*(8.D0*CosPhi2-Three)*CosPhi
      CASE(6) ; CosNPhi = ((32.D0*CosPhi2-48.D0)*CosPhi2+18.D0) *CosPhi2-One
                DCos    = (192.D0*(CosPhi2-One)*CosPhi2+36.D0) *CosPhi
                D2Cos   =  192.D0*(Five*CosPhi2-Three)*CosPhi2 + 36.D0
!
      END SELECT
!
! Calculate energy
!
ECoord=ECoord+CoordFC%D(ICoord)*(One+CoordEQ%D(ICoord)*CosNPhi)
!
      IF(.NOT.CalcGrad) CYCLE
!
! Calculate force
!
DForce=CoordFC%D(ICoord)*CoordEQ%D(ICoord)*DCos
!
      Fact1=IJCOSJK/DJK
      Fact2=KLCOSJK/DJK
      DGradI%D=(FJKL%D-CosPhi*FIJK%D)/NFIJK
      DGradL%D=(FIJK%D-CosPhi*FJKL%D)/NFJKL
      DGradJ%D= DGradI%D*(Fact1-One)+Fact2*DGradL%D
      DGradK%D=-DGradL%D*(Fact2+One)-Fact1*DGradI%D
!
      Grad_Loc%D(1:3,I)=Grad_Loc%D(1:3,I)+DForce*DGradI%D
      Grad_Loc%D(1:3,J)=Grad_Loc%D(1:3,J)+DForce*DGradJ%D
      Grad_Loc%D(1:3,K)=Grad_Loc%D(1:3,K)+DForce*DGradK%D
      Grad_Loc%D(1:3,L)=Grad_Loc%D(1:3,L)+DForce*DGradL%D
!
   END DO
!
   CALL Delete(DVectIJ)
   CALL Delete(DVectJK)
   CALL Delete(DVectKL)
   CALL Delete(FIJK)
   CALL Delete(FJKL)
!
   CALL Delete(DGradI)
   CALL Delete(DGradJ)
   CALL Delete(DGradK)
   CALL Delete(DGradL)
!
   END SUBROUTINE TorsionalEnergy
!
#endif
END MODULE DYNAMO
