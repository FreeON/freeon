MODULE DYNAMO
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
   USE string 
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
REAL ( KIND = DP ) :: EBOND, VIRIAL, EANGLE, EDIHEDRAL, EIMPROPER
TYPE(CRDS) :: GM_MM
!
CONTAINS
!
!-------------------------------------------------------------- 
!
#ifdef MMech
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
     INTEGER :: I,J,K,L,M,N,Natoms,NMAX12,NMAX13,NMAX14
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
     CALL OpenHDF(InfFile)
!
     GM_Loc%NAtms=Natoms
     CALL New(GM_Loc)
     CALL New(DVect,3)
     CALL New(AtmMark,Natoms)
!
     CALL Get(NMAX12,'NMAX12')
     CALL Get(NMAX13,'NMAX13')
     CALL Get(NMAX14,'NMAX14')
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
     ENDIF
     IF(PRESENT(E_C_EXCL)) THEN
       CALL Get(GM_Loc%AtNum,'atomicnumbers',Tag_O='GM_MM'//Cur)
       CALL New(CHARGE14,Natoms)
       CALL Get(Charge14,'Charge14')
       E_C_EXCL=Zero
     ENDIF
!
! Now calculate Lennard-Jones and Coulomb exclusion energies
! Please note that all MM charges on the QM positions
! and on the link atoms, as well as on the atoms of the 
! broken covalent bonds should be set to zero at this point.
!
     CONVF=e2PerAngstroemToKJPerMol
!
! bond, angle and dihedral terms
!
     DO ITop=1,3
!
     IF(ITop==1) THEN
       CALL New(Top_Loc,(/Natoms,NMAX12/)) 
       CALL Get(Top_Loc,'Top12')
     ENDIF
     IF(ITop==2) THEN
       CALL New(Top_Loc,(/Natoms,NMAX13/)) 
       CALL Get(Top_Loc,'Top13')
     ENDIF
     IF(ITop==3) THEN
       CALL New(Top_Loc,(/Natoms,NMAX14/)) 
       CALL Get(Top_Loc,'Top14')
     ENDIF
!
     DO I=1,Natoms
       IF(AtmMark%I(I)/=0) CYCLE !!!interactions of MM atoms only
       IF(PRESENT(E_LJ_EXCL)) THEN
         QI=LJEps%D(I)
         QI14=LJEps14%D(I)
         RI=LJRad%D(I)
       ENDIF
       IF(PRESENT(E_LJ_EXCL)) THEN
         QQI=GM_Loc%AtNum%D(I)
         QQI14=CHARGE14%D(I)
       ENDIF
     DO M=1,Top_Loc%I(I,1)
       J=Top_Loc%I(I,M+1)
     IF(AtmMark%I(J)==0.AND.J<I) CYCLE !!!avoid double counting
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
       IF(PRESENT(E_LJ_EXCL)) THEN
         QQJ=GM_Loc%AtNum%D(J)
         QQJ14=CHARGE14%D(J)
       ENDIF
         IF(DIJ>0.001D0) THEN
           IF(PRESENT(E_LJ_EXCL)) THEN
             R6=(RI*RJ/DIJ)**6
             IF(ITop==3) THEN
               Pref1=(QI*QJ-QI14*QJ14)*R6
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
             IF(ITop==3) THEN
               Pref1=(QQI*QQJ-QQI14*QQJ14)/DIJ*CONVF
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
       write(*,*) 'd dij= ',i,j,dij
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
     CALL CloseHDF()
!
   END SUBROUTINE EXCL
!-------------------------------------------------------------- 
!
   SUBROUTINE ENERGY_LENNARD_JONES(ELJ,ISet,BoxSize,Grad_Loc)
!
   IMPLICIT NONE
   INTEGER :: ISET,NBox,Natoms,NX,NY,NZ
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
   CALL OpenHDF(InfFile)
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
   CALL Get(GM_Loc,'GM_MM'//TRIM(IntToChar(ISet)))
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
       IF(AtmMark%I(JJ1)/=0) CYCLE !!!interactions of MM atoms only
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
       IF(AtmMark%I(JJ2)==0.AND.JJ2<JJ1) CYCLE !!!avoid double counting
         IF(AtmMark%I(JJ2)==0) THEN
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
! All interactions have been counted twice
!
!  write(*,*) 'elj in dynamo= ',elj
!
   CALL CloseHDF()
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
#endif
!
END MODULE DYNAMO
