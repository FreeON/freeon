MODULE DYNAMO
!
USE DerivedTypes
USE GlobalScalars
Use MemMan
USE InOut
USE Macros, ONLY : HasQM,HasMM,MMOnly 
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
   SUBROUTINE COULOMB_EXCL(Natoms,E_C_EXCL,InfFile,Cur)
!
! Calculate Coulomb exclusion for MM and QMMM calculations
! MM_Natoms and MM geometry must be used
!
     IMPLICIT NONE
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: InfFile
     REAL(DOUBLE) :: E_C_EXCL12,E_C_EXCL13,E_C_EXCL14,E_C_EXCL
     TYPE(INT_VECT) :: ATMMARK
     TYPE(CRDS) :: GM_MM 
     INTEGER :: I,J,K,L,M,N,Natoms,NMAX12,NMAX13,NMAX14
     TYPE(DBL_VECT)    :: CHARGE14
     TYPE(INT_RNK2)    :: TOP12,TOP13,TOP14
     CHARACTER(LEN=3)  :: Cur
     REAL(DOUBLE) :: XI,YI,ZI,XJ,YJ,ZJ,DIJ,QI,QJ,CONVF,QI14,QJ14
!
! If QMMM distinguish QM and MM atoms
!   
     CALL OpenHDF(InfFile)
!
     GM_MM%NAtms=Natoms
     CALL New(GM_MM)
     CALL New(ATMMARK,Natoms)
     CALL New(CHARGE14,Natoms)
!
     CALL Get(NMAX12,'NMAX12')
     CALL Get(NMAX13,'NMAX13')
     CALL Get(NMAX14,'NMAX14')
     CALL New(TOP12,(/Natoms,NMAX12/))
     CALL New(TOP13,(/Natoms,NMAX13/))
     CALL New(TOP14,(/Natoms,NMAX14/))
!
     IF(HasQM()) THEN
       CALL Get(ATMMARK,'ATMMARK')
     ELSE 
       ATMMARK%I(:)=0 !!! all MM atoms
     ENDIF
!
! Get MM geometry and charges
! also get 14-charges
! also get topology matrices for 1-2, 1-3 and 1-4 neighbour atoms
!
!    CALL Get(GM_MM,'GM_MM'//Cur)
     CALL Get(GM_MM%Carts,'cartesians',Tag_O='GM_MM'//Cur)
     CALL Get(GM_MM%AtNum,'atomicnumbers',Tag_O='GM_MM'//Cur)
     CALL Get(CHARGE14,'CHARGE14')
     CALL Get(TOP12,'TOP12')
     CALL Get(TOP13,'TOP13')
     CALL Get(TOP14,'TOP14')
!
     CALL CloseHDF() 
!!
!!test data got from HDF
!      write(*,*) 'coordinates test'
!     do i=1,Natoms
!      write(*,200) I,GM_MM%Carts%D(1:3,I)/AngstromsToAU 
!200   format(I5,3F12.6)
!     enddo
!      write(*,*) 'charges, ch14 eps14'
!     do i=1,Natoms
!      write(*,200) I,GM_MM%AtNum%D(I),CHARGE14%D(I)
!     enddo
!      write(*,*) 'top12 '
!     do i=1,Natoms
!      write(*,300) I,TOP12%I(I,1:NMAX12) 
!     enddo
!      write(*,*) 'top13 '
!     do i=1,Natoms
!      write(*,300) I,TOP13%I(I,1:NMAX13) 
!     enddo
!      write(*,*) 'top14 '
!     do i=1,Natoms
!      write(*,300) I,TOP14%I(I,1:NMAX14) 
!     enddo
!300  format(I5,20I4)
!
! Now calculate Coulomb exclusion energy
! Please note that all MM charges on the QM positions
! and on the link atoms, as well as on the atoms of the 
! broken covalent bonds should be set to zero at this point.
!
     E_C_EXCL12=Zero
     E_C_EXCL13=Zero
     E_C_EXCL14=Zero
     CONVF=e2PerAngstroemToKJPerMol*AngstromsToAU
!
! bond terms
!
     DO I=1,Natoms
       IF(ATMMARK%I(I)==0) THEN
       XI=GM_MM%Carts%D(1,I)
       YI=GM_MM%Carts%D(2,I)
       ZI=GM_MM%Carts%D(3,I)
       QI=GM_MM%AtNum%D(I)
     DO M=1,TOP12%I(I,1)
       J=TOP12%I(I,M+1)
       IF(ATMMARK%I(J)==0) THEN
       XJ=GM_MM%Carts%D(1,J)
       YJ=GM_MM%Carts%D(2,J)
       ZJ=GM_MM%Carts%D(3,J)
       QJ=GM_MM%AtNum%D(J)
       DIJ=SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
         IF(DIJ>0.001D0) THEN
           E_C_EXCL12=E_C_EXCL12+QI*QJ/DIJ
         ELSE
       write(*,*) 'bond dij= ',i,j,dij
  Call MondoHalt(QMMM_ERROR,'Atoms are too close to each other in MM set')
         ENDIF
       ENDIF
     ENDDO
       ENDIF
     ENDDO
!
! angle terms
!
     DO I=1,Natoms
       IF(ATMMARK%I(I)==0) THEN
       XI=GM_MM%Carts%D(1,I)
       YI=GM_MM%Carts%D(2,I)
       ZI=GM_MM%Carts%D(3,I)
       QI=GM_MM%AtNum%D(I)
     DO M=1,TOP13%I(I,1)
       J=TOP13%I(I,M+1)
       IF(ATMMARK%I(J)==0) THEN
       XJ=GM_MM%Carts%D(1,J)
       YJ=GM_MM%Carts%D(2,J)
       ZJ=GM_MM%Carts%D(3,J)
       QJ=GM_MM%AtNum%D(J)
       DIJ=SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
         IF(DIJ>0.001D0) THEN
           E_C_EXCL13=E_C_EXCL13+QI*QJ/DIJ
         ELSE
       write(*,*) 'angle dij= ',i,j,dij
  Call MondoHalt(QMMM_ERROR,'Atoms are too close to each other in MM set')
         ENDIF
       ENDIF
     ENDDO
       ENDIF
     ENDDO
!
! dihedral terms
!
     DO I=1,Natoms
       IF(ATMMARK%I(I)==0) THEN
       XI=GM_MM%Carts%D(1,I)
       YI=GM_MM%Carts%D(2,I)
       ZI=GM_MM%Carts%D(3,I)
       QI=GM_MM%AtNum%D(I)
       QI14=CHARGE14%D(I)
     DO M=1,TOP14%I(I,1)
       J=TOP14%I(I,M+1)
       IF(ATMMARK%I(J)==0) THEN
       XJ=GM_MM%Carts%D(1,J)
       YJ=GM_MM%Carts%D(2,J)
       ZJ=GM_MM%Carts%D(3,J)
       QJ=GM_MM%AtNum%D(J)
       QJ14=CHARGE14%D(J)
       DIJ=SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
         IF(DIJ>0.001D0) THEN
           E_C_EXCL14=E_C_EXCL14+(QI*QJ-QI14*QJ14)/DIJ
         ELSE
       write(*,*) 'angle dij= ',i,j,dij
  Call MondoHalt(QMMM_ERROR,'Atoms are too close to each other in MM set')
         ENDIF
       ENDIF
     ENDDO
       ENDIF
     ENDDO
!
     E_C_EXCL12=E_C_EXCL12*CONVF/Two
     E_C_EXCL13=E_C_EXCL13*CONVF/Two
     E_C_EXCL14=E_C_EXCL14*CONVF/Two
     E_C_EXCL=E_C_EXCL12+E_C_EXCL13+E_C_EXCL14
!
     CALL DELETE(GM_MM)
     CALL DELETE(ATMMARK)
     CALL DELETE(CHARGE14)
     CALL DELETE(TOP12)
     CALL DELETE(TOP13)
     CALL DELETE(TOP14)
     write(*,*) 'E_C_EXCL12= ',E_C_EXCL12
     write(*,*) 'E_C_EXCL13= ',E_C_EXCL13
     write(*,*) 'E_C_EXCL14= ',E_C_EXCL14
     write(*,*) 'E_C_EXCL= ',E_C_EXCL
!
     CALL OpenHDF(InfFile)
       CALL Put(E_C_EXCL12,'E_C_EXCL12')
       CALL Put(E_C_EXCL13,'E_C_EXCL13')
       CALL Put(E_C_EXCL14,'E_C_EXCL14')
       CALL Put(E_C_EXCL,'E_C_EXCL')
     CALL CloseHDF()
!
   END SUBROUTINE COULOMB_EXCL
!
!--------------------------------------------------------
!
   SUBROUTINE LJ_EXCL(Natoms,E_LJ_EXCL,InfFile,Cur)
!
! Calculate the Lennard-Jones exclusion for MM and QMMM calculations
! MM_Natoms and MM geometry must be used
!
     IMPLICIT NONE
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: InfFile
     REAL(DOUBLE) :: E_LJ_EXCL12,E_LJ_EXCL13,E_LJ_EXCL14,E_LJ_EXCL
     TYPE(INT_VECT) :: ATMMARK
     TYPE(CRDS) :: GM_MM 
     INTEGER :: I,J,K,L,M,N,Natoms,NMAX12,NMAX13,NMAX14
     TYPE(DBL_VECT)    :: LJEPS14,LJEPS,LJRAD
     TYPE(INT_RNK2)    :: TOP12,TOP13,TOP14
     CHARACTER(LEN=3)  :: Cur
     REAL(DOUBLE) :: XI,YI,ZI,XJ,YJ,ZJ,DIJ,QI,QJ,CONVF,QI14,QJ14
     REAL(DOUBLE) :: RI,RJ,R6
!
! If QMMM distinguish QM and MM atoms
!   
     CALL OpenHDF(InfFile)
!
     GM_MM%NAtms=Natoms
     CALL New(GM_MM)
     CALL New(ATMMARK,Natoms)
     CALL New(LJEPS,Natoms)
     CALL New(LJEPS14,Natoms)
     CALL New(LJRAD,Natoms)
!
     CALL Get(NMAX12,'NMAX12')
     CALL Get(NMAX13,'NMAX13')
     CALL Get(NMAX14,'NMAX14')
     CALL New(TOP12,(/Natoms,NMAX12/))
     CALL New(TOP13,(/Natoms,NMAX13/))
     CALL New(TOP14,(/Natoms,NMAX14/))
!
     IF(HasQM()) THEN
       CALL Get(ATMMARK,'ATMMARK')
     ELSE 
       ATMMARK%I(:)=0 !!! all MM atoms
     ENDIF
!
! Get MM geometry and charges
! also get 14-charges
! also get topology matrices for 1-2, 1-3 and 1-4 neighbour atoms
!
     CALL Get(GM_MM%Carts,'cartesians',Tag_O='GM_MM'//Cur)
     CALL Get(LJRAD,'LJRAD')
     CALL Get(LJEPS,'LJEPS')
     CALL Get(LJEPS14,'LJEPS14')
     CALL Get(TOP12,'TOP12')
     CALL Get(TOP13,'TOP13')
     CALL Get(TOP14,'TOP14')
!
     CALL CloseHDF() 
!!
!!test data got from HDF
!      write(*,*) 'coordinates test'
!     do i=1,Natoms
!      write(*,200) I,GM_MM%Carts%D(1:3,I)/AngstromsToAU 
!200   format(I5,3F12.6)
!     enddo
!      write(*,*) 'charges, ch14 eps14'
!     do i=1,Natoms
!      write(*,200) I,LJEPS%D(I),LJEPS14%D(I)
!     enddo
!      write(*,*) 'top12 '
!     do i=1,Natoms
!      write(*,300) I,TOP12%I(I,1:NMAX12) 
!     enddo
!      write(*,*) 'top13 '
!     do i=1,Natoms
!      write(*,300) I,TOP13%I(I,1:NMAX13) 
!     enddo
!      write(*,*) 'top14 '
!     do i=1,Natoms
!      write(*,300) I,TOP14%I(I,1:NMAX14) 
!     enddo
!300  format(I5,20I4)
!
! Now calculate Lennard-Jones exclusion energy
! Please note that all MM charges on the QM positions
! and on the link atoms, as well as on the atoms of the 
! broken covalent bonds should be set to zero at this point.
!
     E_LJ_EXCL12=Zero
     E_LJ_EXCL13=Zero
     E_LJ_EXCL14=Zero
     CONVF=e2PerAngstroemToKJPerMol*AngstromsToAU
!
! bond terms
!
     DO I=1,Natoms
       IF(ATMMARK%I(I)==0) THEN
       XI=GM_MM%Carts%D(1,I)
       YI=GM_MM%Carts%D(2,I)
       ZI=GM_MM%Carts%D(3,I)
       QI=LJEPS%D(I)
       RI=LJRAD%D(I)
     DO M=1,TOP12%I(I,1)
       J=TOP12%I(I,M+1)
!      IF(ATMMARK%I(J)==0) THEN !!! the second atom may also be QM
       XJ=GM_MM%Carts%D(1,J)
       YJ=GM_MM%Carts%D(2,J)
       ZJ=GM_MM%Carts%D(3,J)
       QJ=LJEPS%D(J)
       RJ=LJRAD%D(J)
       DIJ=SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)/AngstromsToAU
         IF(DIJ>0.001D0) THEN
           R6=(RI*RJ/DIJ)**6
           E_LJ_EXCL12=E_LJ_EXCL12+QI*QJ*R6*(R6-One)
         ELSE
       write(*,*) 'bond dij= ',i,j,dij
  Call MondoHalt(QMMM_ERROR,'Atoms are too close to each other in MM set')
         ENDIF
!      ENDIF !!!ATMMARK(J)
     ENDDO
       ENDIF
     ENDDO
!
! angle terms
!
     DO I=1,Natoms
       IF(ATMMARK%I(I)==0) THEN
       XI=GM_MM%Carts%D(1,I)
       YI=GM_MM%Carts%D(2,I)
       ZI=GM_MM%Carts%D(3,I)
       QI=LJEPS%D(I)
       RI=LJRAD%D(I)
     DO M=1,TOP13%I(I,1)
       J=TOP13%I(I,M+1)
!      IF(ATMMARK%I(J)==0) THEN !!! the second atom may also be QM
       XJ=GM_MM%Carts%D(1,J)
       YJ=GM_MM%Carts%D(2,J)
       ZJ=GM_MM%Carts%D(3,J)
       QJ=LJEPS%D(J)
       RJ=LJRAD%D(J)
       DIJ=SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)/AngstromsToAU
         IF(DIJ>0.001D0) THEN
           R6=(RI*RJ/DIJ)**6
           E_LJ_EXCL13=E_LJ_EXCL13+QI*QJ*R6*(R6-One)
         ELSE
       write(*,*) 'angle dij= ',i,j,dij
  Call MondoHalt(QMMM_ERROR,'Atoms are too close to each other in MM set')
         ENDIF
!      ENDIF !!!ATMMARK(J)
     ENDDO
       ENDIF
     ENDDO
!
! dihedral terms
!
     DO I=1,Natoms
       IF(ATMMARK%I(I)==0) THEN
       XI=GM_MM%Carts%D(1,I)
       YI=GM_MM%Carts%D(2,I)
       ZI=GM_MM%Carts%D(3,I)
       QI=LJEPS%D(I)
       QI14=LJEPS14%D(I)
       RI=LJRAD%D(I)
     DO M=1,TOP14%I(I,1)
       J=TOP14%I(I,M+1)
!      IF(ATMMARK%I(J)==0) THEN !!! the second atom may also be QM
       XJ=GM_MM%Carts%D(1,J)
       YJ=GM_MM%Carts%D(2,J)
       ZJ=GM_MM%Carts%D(3,J)
       QJ=LJEPS%D(J)
       QJ14=LJEPS14%D(J)
       RJ=LJRAD%D(J)
       DIJ=SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)/AngstromsToAU
         IF(DIJ>0.001D0) THEN
           R6=(RI*RJ/DIJ)**6
!          E_LJ_EXCL14=E_LJ_EXCL14+(QI14*QJ14)*R6*(R6-One)
           E_LJ_EXCL14=E_LJ_EXCL14+(QI*QJ-QI14*QJ14)*R6*(R6-One)
         ELSE
       write(*,*) 'dihedral dij= ',i,j,dij
  Call MondoHalt(QMMM_ERROR,'Atoms are too close to each other in MM set')
         ENDIF
!      ENDIF !!!ATMMARK(J)
     ENDDO
       ENDIF
     ENDDO
!
     E_LJ_EXCL12=E_LJ_EXCL12/Two
     E_LJ_EXCL13=E_LJ_EXCL13/Two
     E_LJ_EXCL14=E_LJ_EXCL14/Two
     E_LJ_EXCL=E_LJ_EXCL12+E_LJ_EXCL13+E_LJ_EXCL14
!
     CALL DELETE(GM_MM)
     CALL DELETE(ATMMARK)
     CALL DELETE(LJEPS)
     CALL DELETE(LJEPS14)
     CALL DELETE(LJRAD)
     CALL DELETE(TOP12)
     CALL DELETE(TOP13)
     CALL DELETE(TOP14)
     write(*,*) 'E_LJ_EXCL12= ',E_LJ_EXCL12
     write(*,*) 'E_LJ_EXCL13= ',E_LJ_EXCL13
     write(*,*) 'E_LJ_EXCL14= ',E_LJ_EXCL14
     write(*,*) 'E_LJ_EXCL= ',E_LJ_EXCL
!
     CALL OpenHDF(InfFile)
       CALL Put(E_LJ_EXCL12,'E_LJ_EXCL12')
       CALL Put(E_LJ_EXCL13,'E_LJ_EXCL13')
       CALL Put(E_LJ_EXCL14,'E_LJ_EXCL14')
       CALL Put(E_LJ_EXCL,'E_LJ_EXCL')
     CALL CloseHDF()
!
   END SUBROUTINE LJ_EXCL
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
   TYPE(INT_VECT) :: BoxI,BoxJ,ATMMARK
   TYPE(CRDS) :: GM_Loc
   TYPE(DBL_VECT) :: LJEPS,LJRAD,DVect
   TYPE(DBL_RNK2),OPTIONAL :: Grad_Loc
   REAL(DOUBLE) :: ELJ,r12_2,Pref1,Pref2
!
   CALL OpenHDF(InfFile)
!
! BoxSize is in angstroems
!
! get MM or QMMM atom coordinates as stored in GM_Loc
! get ATMMARK to distinguish between QM and MM atoms
! get LJEPS%D : LJ epsilon
! get LJRAD%D : LJ radii of atoms
! calculate and get BoxI and BoxJ for a certain Box size, then
! sum up LJ contributions within Boxes and between 
! neighbouring Boxes
!
   CALL Get(GM_Loc,'GM_MM'//TRIM(IntToChar(ISet)))
   Natoms=GM_Loc%Natms
!
     CALL New(ATMMARK,Natoms)
   IF(HasQM()) THEN
     CALL Get(ATMMARK,'ATMMARK')
   ELSE 
     ATMMARK%I(:)=0 !!! all MM atoms
   ENDIF
!
   CALL New(DVect,3)
   CALL New(LJEPS,Natoms)
   CALL New(LJRAD,Natoms)
   CALL Get(LJEPS,'LJEPS')
   CALL Get(LJRAD,'LJRAD') 
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
        R1=LJRAD%D(JJ1)
        Q1=LJEPS%D(JJ1)
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
       IF(JJ1>JJ2) CYCLE !!!! do not count twice
         IF(ATMMARK%I(JJ2)==0) THEN
           R2=LJRAD%D(JJ2)
           Q2=LJEPS%D(JJ2)
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
   CALL Delete(ATMMARK)
   CALL Delete(DVect)
   CALL Delete(LJEPS)
   CALL Delete(LJRAD)
   CALL Delete(BoxI)
   CALL Delete(BoxJ)
!
   END SUBROUTINE ENERGY_LENNARD_JONES
!
#endif
!
END MODULE DYNAMO
