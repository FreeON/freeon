!     Last change:  AA   23 May 2002    2:31 pm
MODULE ChargeUtil
  USE MatrixUtil
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: output,bader,wigner_seitz,read_in,write_max_rho,write_atm_rho
  CONTAINS
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_in(chargefile,lattice,ndim,wdim,ngxf,ngyf,ngzf,nrho,rho,Rcar,Rdir,&
                     ws_charge,vasp)
    CHARACTER(LEN=20),INTENT(IN) :: chargefile
    REAL(q2),INTENT(OUT),DIMENSION(3,3) :: lattice
    REAL(q2),POINTER,DIMENSION(:,:,:) :: rho
    REAL(q2),POINTER,DIMENSION(:,:) :: Rcar,Rdir,ws_charge
    INTEGER,INTENT(OUT) :: ndim,wdim,ngxf,ngyf,ngzf,nrho
    LOGICAL,INTENT(OUT) :: vasp

    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: box,steps,v
    REAL(q2) :: side,vol,bohr2ang=0.529177_q2
    INTEGER :: i
    INTEGER,DIMENSION(110) :: elements
    CHARACTER(LEN=7) :: text

    vasp=.false.

    OPEN(100,FILE=chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
    WRITE(*,'(1A11,1A20)') 'OPEN ... ',chargefile
    READ(100,'(6/,1A7)') text
    REWIND(100)
    IF (text == 'Direct') THEN
      elements=0
      vasp=.true.
      WRITE(*,'(1A23)') 'VASP-STYLE INPUT FILE'
      READ(100,'(/,1F20.16)') side
      READ(100,'(3F13.6)') (lattice(i,1:3) , i=1,3)
      READ(100,'(110I4)') elements
      READ(100,*)
!      READ(100,'(1I4,/)') ndim
      ndim=SUM(elements)
      lattice=side*lattice
      CALL transform_matrix(lattice,B,3,3)
      wdim=ndim
      ALLOCATE(Rcar(ndim,3),Rdir(ndim,3),ws_charge(wdim,4))
      DO i=1,ndim
!   Shouldn't Rdir be multiplied by side?
        READ(100,'(3(2X,1F8.6))') Rdir(i,:)
        CALL matrix_vector(B,Rdir(i,:),v,3,3)
        Rcar(i,:)=v
      END DO
      READ(100,*) 
      READ(100,*) ngxf,ngyf,ngzf
!      READ(100,'(/,3I5)') ngxf,ngyf,ngzf
    ELSE
! temp. readin
      WRITE(*,'(1A27)') 'GAUSSIAN-STYLE INPUT FILE'
      READ(100,'(2/,1I5,3(3X,1F9.6))') ndim,box
      wdim=ndim
      ALLOCATE(Rcar(ndim,3),Rdir(ndim,3),ws_charge(wdim,4))
      READ(100,'(1I5,3X,1F9.6)') ngxf,steps(1)
      READ(100,'(1I5,15X,1F9.6)') ngyf,steps(2)
      READ(100,'(1I5,27X,1F9.6)') ngzf,steps(3)
      box(1)=REAL((ngxf),q2)*steps(1)
      box(2)=REAL((ngyf),q2)*steps(2)
      box(3)=REAL((ngzf),q2)*steps(3)
      lattice=0.0_q2
      DO i=1,3
        lattice(i,i)=box(i)
      END DO
      vol=volume(lattice)
      DO i=1,3
        lattice(i,i)=lattice(i,i)-steps(i)
      END DO
      lattice=lattice*bohr2ang
      CALL transform_matrix(lattice,B,3,3)
      DO i=1,ndim
        READ(100,'(17X,3(3X,1F9.6))') Rdir(i,:)
        Rdir(i,:)=Rdir(i,:)/(box-steps)+0.5_q2
        CALL matrix_vector(B,Rdir(i,:),v,3,3)
        Rcar(i,:)=v
      END DO
!
    END IF
    nrho=ngxf*ngyf*ngzf
    ALLOCATE(rho(ngxf,ngyf,ngzf))
    CALL charge_density(rho,ngxf,ngyf,ngzf,nrho,vasp)
    IF (.not.vasp) rho=rho*vol
!      lattice=lattice*bohr2ang

    WRITE(*,'(1A12,1A20)') 'CLOSE ... ', chargefile
    CLOSE(100)

!    CALL write_max_rho(rho,ngxf,ngyf,ngzf)

!    mag=2
!    CALL make_rhop(rho,ngxf,ngyf,ngzf,rhop,ngxp,ngyp,ngzp,mag)
!    nrho=ngxp*ngyp*ngzp

  RETURN
  END SUBROUTINE read_in

!------------------------------------------------------------------------------------!

  SUBROUTINE charge_density(rho,ngxf,ngyf,ngzf,nrho,vasp)
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf,nrho
    REAL(q2),INTENT(OUT),DIMENSION(ngxf,ngyf,ngzf) :: rho
    LOGICAL,INTENT(IN) :: vasp

    INTEGER :: nx,ny,nz,i,c,j,readstat,endl
    REAL(q2),POINTER,DIMENSION(:) :: r
    CHARACTER(LEN=20) :: frmt

    IF (vasp) THEN
      ALLOCATE(r(5))
      READ(100,'(5E18.11)',IOSTAT=readstat) r
      BACKSPACE(100)
      IF (readstat == 0) THEN
        endl=5
        frmt='(5E18.11)'
      ELSE IF (readstat > 0) THEN
        endl=10
        frmt='(10F8.3)'
        DEALLOCATE(r)
        ALLOCATE(r(endl))
      ELSE 
        WRITE(*,*) 'END OF FILE'
        STOP
      END IF
      i=0
      nz=1
      ny=1
      nx=0
      c=0
      DO
        i=i+1
        r=0.0_q2
        READ(100,frmt) r
        DO j=1,endl
          nx=nx+1
          c=c+1
          IF (nx > ngxf) THEN
            nx=1
            ny=ny+1
            IF (ny > ngyf) THEN
              ny=1
              nz=nz+1
            END IF
          END IF
          rho(nx,ny,nz)=r(j)
          IF (c == nrho) EXIT
        END DO
        IF (c == nrho) EXIT  
      END DO
      
    ELSE
      DO nx=1,ngxf
        DO ny=1,ngyf
          READ(100,'(6E13.5)') (rho(nx,ny,nz) , nz=1,ngzf)
        END DO
      END DO
    END IF

!    IF (vasp) THEN
!      DO nz=1,ngzf
!        write(*,*) 'nz: ',nz
!        DO ny=1,ngyf
!          write(*,*) 'ny: ',ny
!          READ(100,'(5E18.11)') (rho(nx,ny,nz) , nx=1,ngxf)
!          write(*,'(5E18.11)') (rho(nx,ny,nz) , nx=1,ngxf)
!          pause
!        END DO
!      END DO
!    ELSE
!      DO nx=1,ngxf
!        DO ny=1,ngyf
!          READ(100,'(6E13.5)') (rho(nx,ny,nz) , nz=1,ngzf)
!        END DO
!      END DO
!    END IF

  RETURN
  END SUBROUTINE charge_density

!-----------------------------------------------------------------------------------!

  SUBROUTINE bader(bader_charge,bdim,rho,ngxf,ngyf,ngzf,nrho,Rdir,ndim,lattice,     &
  &                bader_achg,bader_atom,bader_dist,max_rho)
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf,nrho,ndim
    REAL(q2),INTENT(IN),DIMENSION(ngxf,ngyf,ngzf) :: rho
    REAL(q2),POINTER,DIMENSION(:,:) :: bader_charge
    REAL(q2),INTENT(IN),DIMENSION(ndim,3) :: Rdir
    REAL(q2),INTENT(IN),DIMENSION(3,3) :: lattice
    INTEGER,INTENT(OUT) :: bdim

!    INTEGER,POINTER,DIMENSION(:,:,:) :: max_rho
    INTEGER,DIMENSION(ngxf,ngyf,ngzf) :: max_rho
    INTEGER,POINTER,DIMENSION(:,:) :: path
    REAL(q2),POINTER,DIMENSION(:,:) :: tmp
    REAL(q2),POINTER,DIMENSION(:) :: bader_dist,bader_achg
    INTEGER,POINTER,DIMENSION(:) :: bader_atom
    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: v,rngf
    INTEGER :: nx,ny,nz,px,py,pz,i,pdim,pnum,known_max,bnum,p,tenths_done

    WRITE(*,'(1A39)') 'CALCULATING BADER CHARGE DISTRIBUTION'
    bdim=64
    pdim=64
    ALLOCATE(bader_charge(bdim,4))
    ALLOCATE(path(pdim,3))
    bader_charge=0.0_q2
    max_rho=0
    bnum=0
    tenths_done=0
    DO nx=1,ngxf
      IF ((nx*10/ngxf) > tenths_done) THEN
        tenths_done=(nx*10/ngxf)
        WRITE(*,'(1X,1I4,1A6)') (tenths_done*10),'% done'
      END IF
      DO ny=1,ngyf
        DO nz=1,ngzf
          px=nx
          py=ny
          pz=nz
          IF(max_rho(px,py,pz) == 0) THEN
            CALL maximize(px,py,pz,rho,ngxf,ngyf,ngzf,path,pdim,pnum,max_rho)
            CALL pbc(px,py,pz,ngxf,ngyf,ngzf)
            known_max=max_rho(px,py,pz)
            IF (known_max == 0) THEN
              bnum=bnum+1
              known_max=bnum
              IF (bnum > bdim) THEN
                ALLOCATE(tmp(bdim,4))
                tmp=bader_charge
                DEALLOCATE(bader_charge)
                bdim=2*bdim
                ALLOCATE(bader_charge(bdim,4))
                bader_charge=0.0_q2
                bader_charge(1:bnum-1,1:4)=tmp
                DEALLOCATE(tmp)
              END IF
              bader_charge(bnum,1)=REAL(px,q2)
              bader_charge(bnum,2)=REAL(py,q2)
              bader_charge(bnum,3)=REAL(pz,q2)
            END IF
            DO p=1,pnum
              max_rho(path(p,1),path(p,2),path(p,3))=known_max
            END DO
          END IF
        END DO
      END DO
    END DO

    DO nx=1,ngxf
      DO ny=1,ngyf
        DO nz=1,ngzf
          bader_charge(max_rho(nx,ny,nz),4)=bader_charge(max_rho(nx,ny,nz),4)+      &
  &                    rho(nx,ny,nz)
        END DO
      END DO
    END DO

    ALLOCATE(tmp(bdim,4))
    tmp=bader_charge
    DEALLOCATE(bader_charge)
    ALLOCATE(bader_charge(bnum,4))
    bader_charge=tmp(1:bnum,:)
    bdim=bnum
    ALLOCATE(bader_atom(bdim),bader_dist(bdim),bader_achg(ndim))

! Don't have this normalization in MONDO
!    bader_charge(:,4)=bader_charge(:,4)/REAL(nrho,q2)

    rngf(1)=REAL(ngxf,q2)
    rngf(2)=REAL(ngyf,q2)
    rngf(3)=REAL(ngzf,q2)
    DO i=1,bdim
      bader_charge(i,1:3)=(bader_charge(i,1:3)-1.0_q2)/rngf
    END DO

    CALL charge2atom(bader_charge,bdim,Rdir,ndim,bader_achg,bader_dist,bader_atom,  &
    &                lattice)

    CALL transform_matrix(lattice,B,3,3)
    DO i=1,bdim
      CALL matrix_vector(B,bader_charge(i,1:3),v,3,3)
      bader_charge(i,1:3)=v
    END DO

!    CALL write_max_rho(max_rho,ngxf,ngyf,ngzf)

  RETURN
  END SUBROUTINE bader

!-----------------------------------------------------------------------------------!

  SUBROUTINE maximize(px,py,pz,rho,ngxf,ngyf,ngzf,path,pdim,pnum,max_rho)
    INTEGER,INTENT(INOUT) :: px,py,pz,pdim,pnum
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf
    REAL(q2),INTENT(IN),DIMENSION(ngxf,ngyf,ngzf) :: rho
    INTEGER,INTENT(IN),DIMENSION(ngxf,ngyf,ngzf) :: max_rho
    INTEGER,POINTER,DIMENSION(:,:) :: path

    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: tmp

    pnum=1
    path(pnum,1:3)=(/px,py,pz/)

    DO
      IF(max_neighbour(px,py,pz,rho,ngxf,ngyf,ngzf)) THEN
        pnum=pnum+1
        IF (pnum > pdim) THEN
          ALLOCATE(tmp(pdim,3))
          tmp=path
          DEALLOCATE(path)
          pdim=2*pdim
          ALLOCATE(path(pdim,3))
          path=0.0_q2
          path(1:pnum-1,:)=tmp
          DEALLOCATE(tmp)
        END IF
	CALL pbc(px,py,pz,ngxf,ngyf,ngzf)
        path(pnum,1:3)=(/px,py,pz/)
        IF(max_rho(px,py,pz) /= 0) EXIT
      ELSE
        EXIT
      END IF
    END DO

  RETURN
  END SUBROUTINE maximize

!-----------------------------------------------------------------------------------!

  FUNCTION max_neighbour(px,py,pz,rho,ngxf,ngyf,ngzf)
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf
    REAL(q2),INTENT(IN),DIMENSION(ngxf,ngyf,ngzf) :: rho
    INTEGER,INTENT(INOUT) :: px,py,pz
    LOGICAL :: max_neighbour

    REAL(q2) :: rho_max,rho_tmp,rho_ctr
    INTEGER :: dx,dy,dz,pxt,pyt,pzt,pxm,pym,pzm
    REAL(q2),DIMENSION(-1:1,-1:1,-1:1),SAVE :: w=RESHAPE((/           &
    &    0.5773502691896_q2,0.7071067811865_q2,0.5773502691896_q2,    &
    &    0.7071067811865_q2,1.0000000000000_q2,0.7071067811865_q2,    &
    &    0.5773502691896_q2,0.7071067811865_q2,0.5773502691896_q2,    &
    &    0.7071067811865_q2,1.0000000000000_q2,0.7071067811865_q2,    &
    &    1.0000000000000_q2,0.0000000000000_q2,1.0000000000000_q2,    &
    &    0.7071067811865_q2,1.0000000000000_q2,0.7071067811865_q2,    &
    &    0.5773502691896_q2,0.7071067811865_q2,0.5773502691896_q2,    &
    &    0.7071067811865_q2,1.0000000000000_q2,0.7071067811865_q2,    &
    &    0.5773502691896_q2,0.7071067811865_q2,0.5773502691896_q2     &
    &    /),(/3,3,3/))

    rho_max=0.0_q2
    pxm=px
    pym=py
    pzm=pz
    rho_ctr=rho_value(px,py,pz,rho,ngxf,ngyf,ngzf)
    DO dx=-1,1
      pxt=px+dx
      DO dy=-1,1
        pyt=py+dy
        DO dz=-1,1
          pzt=pz+dz
          rho_tmp=rho_value(pxt,pyt,pzt,rho,ngxf,ngyf,ngzf)
          rho_tmp=rho_ctr+w(dx,dy,dz)*(rho_tmp-rho_ctr)
          IF (rho_tmp > rho_max) THEN
            rho_max=rho_tmp
            pxm=pxt
            pym=pyt
            pzm=pzt
          END IF
        END DO
      END DO
    END DO

    max_neighbour=((pxm /= px) .or. (pym /= py) .or. (pzm /= pz))
    IF (max_neighbour) THEN
      px=pxm
      py=pym
      pz=pzm
    END IF

  RETURN
  END FUNCTION max_neighbour

!-----------------------------------------------------------------------------------!

  SUBROUTINE wigner_seitz(ws_charge,wdim,rho,ngxf,ngyf,ngzf,Rcar,Rdir,ndim,nrho,vasp)
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf,ndim,nrho
    INTEGER,INTENT(OUT) :: wdim
    REAL(q2),INTENT(IN),DIMENSION(ndim,3) :: Rcar,Rdir
    REAL(q2),INTENT(IN),DIMENSION(ngxf,ngyf,ngzf) :: rho
!    REAL(q2),INTENT(OUT),DIMENSION(wdim,4) :: ws_charge
    REAL(q2),POINTER,DIMENSION(:,:) :: ws_charge
    LOGICAL,INTENT(IN) :: vasp

    REAL(q2),DIMENSION(ndim,3) :: Ratm
    REAL(q2),DIMENSION(3) :: Rcur,dR,ngf,ngf_2
    REAL(q2) :: dist,min_dist,add
    INTEGER :: i,nx,ny,nz,closest,tenths_done

    WRITE(*,'(1A46)') 'CALCULATING WIGNER-SEITZ CHARGE DISTRIBUTION'

    wdim=ndim
    ALLOCATE(ws_charge(wdim,4))

    ngf(1)=REAL(ngxf,q2)
    ngf(2)=REAL(ngyf,q2)
    ngf(3)=REAL(ngzf,q2)
    ngf_2=REAL(ngf,q2)/2.0_q2

    add=0.5_q2
    IF (vasp) add=1.0_q2

    Ratm(:,1)=Rdir(:,1)*ngf(1)+add
    Ratm(:,2)=Rdir(:,2)*ngf(2)+add
    Ratm(:,3)=Rdir(:,3)*ngf(3)+add

    ws_charge=0.0_q2
    tenths_done=0
    DO nx=1,ngxf
      Rcur(1)=REAL(nx,q2)
      IF ((nx*10/ngxf) > tenths_done) THEN
        tenths_done=(nx*10/ngxf)
        WRITE(*,'(1X,1I4,1A6)') (tenths_done*10),'% done'
      END IF
      DO ny=1,ngyf
        Rcur(2)=REAL(ny,q2)
        DO nz=1,ngzf
          Rcur(3)=REAL(nz,q2)
          closest=1
          dR=Rcur-Ratm(1,:)
          CALL dpbc(dR,ngf,ngf_2)
          min_dist=DOT_PRODUCT(dR,dR)
          DO i=2,wdim
            dR=Rcur-Ratm(i,:)
            CALL dpbc(dR,ngf,ngf_2)
	    dist=DOT_PRODUCT(dR,dR)
            IF (dist < min_dist) THEN
              min_dist=dist
              closest=i
            END IF
          END DO
          ws_charge(closest,4)=ws_charge(closest,4)+                    &
                               rho_value(nx,ny,nz,rho,ngxf,ngyf,ngzf)
        END DO
      END DO
    END DO
! Don't have this normalization for MONDO
!    ws_charge(:,4)=ws_charge(:,4)/REAL(nrho,q2)
    ws_charge(:,1:3)=Rcar

  RETURN
  END SUBROUTINE wigner_seitz

!-----------------------------------------------------------------------------------!

  FUNCTION rho_value(px,py,pz,rho,ngxf,ngyf,ngzf)
    INTEGER,INTENT(IN) :: px,py,pz
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf
    REAL(q2),INTENT(IN),DIMENSION(ngxf,ngyf,ngzf) :: rho
    REAL(q2) :: rho_value

    INTEGER :: pxt,pyt,pzt

    pxt=px
    pyt=py
    pzt=pz
    DO
      IF(pxt >= 1) EXIT
      pxt=pxt+ngxf
    END DO
    DO
      IF(pxt <= ngxf) EXIT
      pxt=pxt-ngxf
    END DO
    DO
      IF(pyt >= 1) EXIT
      pyt=pyt+ngyf
    END DO
    DO
      IF(pyt <= ngyf) EXIT
      pyt=pyt-ngyf
    END DO
    DO
      IF(pzt >= 1) EXIT
      pzt=pzt+ngzf
    END DO
    DO
      IF(pzt <= ngzf) EXIT
      pzt=pzt-ngzf
    END DO

    rho_value=rho(pxt,pyt,pzt)

  RETURN
  END FUNCTION rho_value

!-----------------------------------------------------------------------------------!

  SUBROUTINE output(bader_charge,bdim,ws_charge,wdim,ndim,lattice,ngxf,ngyf,ngzf,   &
  &                 bader_achg,bader_dist,bader_atom,bader_tol)
    INTEGER,INTENT(IN) :: bdim,wdim,ndim,ngxf,ngyf,ngzf
    REAL(q2),INTENT(IN),DIMENSION(bdim,4) :: bader_charge
    REAL(q2),INTENT(IN),DIMENSION(wdim,4) :: ws_charge
    REAL(q2),INTENT(IN),DIMENSION(3,3) :: lattice
    REAL(q2),INTENT(IN),DIMENSION(ndim) :: bader_achg,bader_dist
    INTEGER,INTENT(IN),DIMENSION(ndim) :: bader_atom

    INTEGER,DIMENSION(bdim) :: bader_found
    REAL(q2),DIMENSION(3) :: rngf,dxyz,v
    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2) :: sum_achg,bader_tol
    INTEGER :: i,j,bdimsig
    LOGICAL :: accept

    OPEN(100,FILE='ACF.dat',STATUS='replace',ACTION='write')
    WRITE(*,555) '#','X','Y','Z','WS','BADER','%'
    555 FORMAT(2/,3X,1A,10X,1A1,2(11X,1A1),10X,1A2,9X,1A5,8X,1A1)
    WRITE(*,666) '----------------------------------------------------------------',&
  &              '-------------'
    666 FORMAT(1A66,1A13)
    sum_achg=SUM(bader_achg)
    DO i=1,wdim
      WRITE(*,'(1I5,6F12.4)') i,ws_charge(i,1:4),bader_achg(i),100.*bader_achg(i)   &
  &                                                                 /sum_achg
      WRITE(100,'(1I5,6F12.4)') i,ws_charge(i,1:4),bader_achg(i),100.*bader_achg(i) &
  &                                                                   /sum_achg
    END DO
    CLOSE(100)
    OPEN(200,FILE='BCF.dat',STATUS='replace',ACTION='write')
    DO i=1,bdim
        if(Bader_charge(i,4)>Bader_Tol) then
           bdimsig=bdimsig+1
           WRITE(200,777) bdimsig,bader_charge(i,:),bader_atom(i),bader_dist(i)
           777 FORMAT(1I5,4F12.4,1I5,1F12.4)
        endif
    END DO
    CLOSE(200)

    WRITE(*,'(/,1A32,6X,1I8)') 'NUMBER OF BADER MAXIMA FOUND: ',bdim
    WRITE(*,'(1A32,6X,1I8)') '    SIGNIFICANT MAXIMA FOUND: ',bdimsig
    WRITE(*,'(1A23,11X,1F12.5,/)') 'NUMBER OF ELECTRONS: ',SUM(bader_charge(1:bdim,4))


!    OPEN(100,FILE='output.dat',STATUS='replace',ACTION='write')
!
!    CALL transform_matrix(lattice,B,3,3)
!    rngf(1)=REAL(ngxf,q2)
!    rngf(2)=REAL(ngyf,q2)
!    rngf(3)=REAL(ngzf,q2)
!    dxyz=1.0_q2/rngf
!    CALL matrix_vector(B,dxyz,v,3,3)
!    dxyz=v
!
!    bader_found=0
!    WRITE(*,555) '#','X','Y','Z','WS','BADER'
!    555 FORMAT(2/,3X,1A,10X,1A1,2(11X,1A1),10X,1A2,9X,1A5)
!    WRITE(*,666) '-----------------------------------------------------------------'
!    666 FORMAT(1A67)
!    DO i=1,wdim
!      DO j=1,bdim
!        accept=ALL(ABS(bader_charge(j,1:3)-ws_charge(i,1:3)) < dxyz)
!        IF (accept) THEN
!          WRITE(*,'(1I5,5F12.4)') i,ws_charge(i,1:4),bader_charge(j,4)
!          WRITE(100,'(1I5,5F12.4)') i,ws_charge(i,1:4),bader_charge(j,4)
!          bader_found(j)=1
!          EXIT
!        END IF
!      END DO
!      IF (.not.accept) THEN
!        WRITE(*,'(1I5,5F12.4)') i,ws_charge(i,1:4),0.0_q2
!        WRITE(100,'(1I5,5F12.4)') i,ws_charge(i,1:4),0.0_q2
!      END IF
!    END DO
!
!    IF (ANY(bader_found == 0)) THEN
!      i=wdim
!      WRITE(*,'(/,1A47)') 'BADER CHARGE REGIONS NOT ASSOCIATED WITH IONS'
!      DO j=1,bdim
!        IF(bader_found(j) == 0) THEN
!          i=i+1
!          WRITE(*,'(1I5,5F12.4)') i,bader_charge(j,1:3),0.0_q2,bader_charge(j,4)
!          WRITE(100,'(1I5,5F12.4)') -1,bader_charge(j,1:3),0.0_q2,bader_charge(j,4)
!        END IF
!      END DO
!    ELSE
!      WRITE(*,'(/,1A44)') 'ALL BADER CHARGES ARE ASSOCIATED WITH IONS'
!    END IF
!
!    CLOSE(100)

  RETURN
  END SUBROUTINE output

!-----------------------------------------------------------------------------------!

  SUBROUTINE write_max_rho(int_rho,ngxf,ngyf,ngzf)
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf
!    REAL(q2),INTENT(IN),DIMENSION(ngxf,ngyf,ngzf) :: max_rho
    INTEGER,DIMENSION(ngxf,ngyf,ngzf) :: int_rho
    REAL(q2),DIMENSION(ngxf,ngyf,ngzf) :: max_rho
!    INTEGER,INTENT(IN),DIMENSION(ngxf,ngyf,ngzf) :: max_rho
    INTEGER :: nx,ny,nz,i,j,nrho,c
    REAL(q2),DIMENSION(5) :: r

    max_rho=REAL(int_rho)
    OPEN(100,FILE='max_rho.dat',STATUS='replace',ACTION='write')

    write(*,*)'Writing Bader Volumes'

    nrho=ngxf*ngyf*ngzf

    i=0
    nz=1
    ny=1
    nx=0
    c=0
    DO
      i=i+1
      DO j=1,5
        nx=nx+1
        c=c+1
        IF (nx > ngxf) THEN
          nx=1
          ny=ny+1
          IF (ny > ngyf) THEN
            ny=1
            nz=nz+1
          END IF
        END IF
        r(j)=max_rho(nx,ny,nz)
        IF (c == nrho) EXIT
      END DO
      WRITE(100,'(5E18.11)') REAL(r,q2)
      IF (c == nrho) EXIT
    END DO

    IF (j==1) WRITE(100,'(1E18.11)') REAL(r(1),q1)
    IF (j==2) WRITE(100,'(2E18.11)') REAL(r(1:2),q1)
    IF (j==3) WRITE(100,'(3E18.11)') REAL(r(1:3),q1)
    IF (j==4) WRITE(100,'(4E18.11)') REAL(r(1:4),q1)

    CLOSE(100)

  RETURN
  END SUBROUTINE write_max_rho

!-----------------------------------------------------------------------------------!


  SUBROUTINE write_atm_rho(int_rho,ngxf,ngyf,ngzf,lattice,Rdir,ndim,bader_charge,Bdim,bader_tol)
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf,ndim,Bdim
    REAL(q2),DIMENSION(ngxf,ngyf,ngzf) :: max_rho
    INTEGER,INTENT(IN),DIMENSION(ngxf,ngyf,ngzf) :: int_rho
    REAL(q2),POINTER,DIMENSION(:,:) :: bader_charge
    REAL(q2),POINTER,DIMENSION(:,:) :: Rdir
    REAL(q2),INTENT(IN),DIMENSION(3,3) :: lattice
    REAL(q2) :: bader_tol
    INTEGER :: nx,ny,nz,i,j,nrho,c,AtomNum,BaderCur
    REAL(q2),DIMENSION(5) :: r
    CHARACTER(15) :: AtomFileName,AtomNumText

    max_rho=REAL(int_rho)
    AtomNum=0
    DO BaderCur=1,bdim
      if(Bader_charge(BaderCur,4)>Bader_Tol) then
        AtomNum=AtomNum+1
        if(AtomNum<10) then
          write(AtomNumText,'(1A3,I1)') '000',AtomNum
        else if(AtomNum<100) then
          write(AtomNumText,'(1A2,I2)') '00',AtomNum
        else if(AtomNum<1000) then
          write(AtomNumText,'(1A1,I3)') '0',AtomNum
        else
          write(AtomNumText,'(I4)') AtomNum
        endif
!        AtomFileName = "Atom"//Trim(AtomNumText(2:))//".dat"
        AtomFileName = "Atom"//Trim(AtomNumText(1:))//".dat"
!        write(*,*)'AtomFileName ',AtomFileName
        open(100, file=AtomFileName)
        
        write(*,'(1A33,I4')'Writing Bader volume around Atom',AtomNum

        WRITE(100,*)'Bader Volume of Atom',AtomNum
        WRITE(100,*)'1.00'
        WRITE(100,'(3F13.6)') (lattice(i,1:3) , i=1,3)
        WRITE(100,'(1I4,/)') NDim
        WRITE(100,*)'DIRECT'
        DO i=1,ndim
          WRITE(100,'(3(2X,1F8.6))') Rdir(i,:)
!          CALL matrix_vector(B,Rdir(i,:),v,3,3)
!          Rcar(i,:)=v
        END DO
        WRITE(100,*)
        WRITE(100,*) ngxf,ngyf,ngzf

        nrho=ngxf*ngyf*ngzf
        i=0
        nz=1
        ny=1
        nx=0
        c=0
        DO
          i=i+1
          DO j=1,5
            nx=nx+1
            c=c+1
            IF (nx > ngxf) THEN
              nx=1
              ny=ny+1
              IF (ny > ngyf) THEN
                ny=1
                nz=nz+1
              END IF
            END IF
            if(int_rho(nx,ny,nz)==BaderCur) then
              r(j)=100
            else
              r(j)=0
            endif
!            r(j)=max_rho(nx,ny,nz)
            IF (c == nrho) EXIT
          END DO
          WRITE(100,'(5E18.11)') REAL(r,q2)
          IF (c == nrho) EXIT
        END DO
        IF (j==1) WRITE(100,'(1E18.11)') REAL(r(1),q1)
        IF (j==2) WRITE(100,'(2E18.11)') REAL(r(1:2),q1)
        IF (j==3) WRITE(100,'(3E18.11)') REAL(r(1:3),q1)
        IF (j==4) WRITE(100,'(4E18.11)') REAL(r(1:4),q1)

        CLOSE(100)
      endif
    END DO

  RETURN
  END SUBROUTINE write_atm_rho

!-----------------------------------------------------------------------------------!

  SUBROUTINE pbc(px,py,pz,ngxf,ngyf,ngzf)
    INTEGER,INTENT(INOUT) :: px,py,pz
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf

    DO
      IF(px > 0) EXIT
      px=px+ngxf
    END DO
    DO
      IF(px <= ngxf) EXIT
      px=px-ngxf
    END DO
    DO
      IF(py > 0) EXIT
      py=py+ngyf
    END DO
    DO
      IF(py <= ngyf) EXIT
      py=py-ngyf
    END DO
    DO
      IF(pz > 0) EXIT
      pz=pz+ngzf
    END DO
    DO
      IF(pz <= ngzf) EXIT
      pz=pz-ngzf
    END DO

  RETURN
  END SUBROUTINE pbc

!-----------------------------------------------------------------------------------!

  SUBROUTINE dpbc_dir(dR)
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dR

    INTEGER :: i

    DO i=1,3
      DO
        IF(dR(i) > -0.5_q2) EXIT
        dR(i)=dR(i)+1.0_q2
      END DO
      DO
        IF(dR(i) < 0.5_q2) EXIT
        dR(i)=dR(i)-1.0_q2
      END DO
    END DO
  RETURN
  END SUBROUTINE dpbc_dir

!-----------------------------------------------------------------------------------!

  SUBROUTINE dpbc(dR,ngf,ngf_2)
    REAL(q2),INTENT(IN),DIMENSION(3) :: ngf,ngf_2
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dR

    INTEGER :: i

    DO i=1,3
      DO
        IF(dR(i) > -ngf_2(i)) EXIT
        dR(i)=dR(i)+ngf(i)
      END DO
      DO
        IF(dR(i) < ngf_2(i)) EXIT
        dR(i)=dR(i)-ngf(i)
      END DO
    END DO

  RETURN
  END SUBROUTINE dpbc

!-----------------------------------------------------------------------------------!

  SUBROUTINE pbcq2(px,py,pz,rngxf,rngyf,rngzf)
    REAL(q2),INTENT(INOUT) :: px,py,pz
    REAL(q2),INTENT(IN) :: rngxf,rngyf,rngzf

    DO
      IF(px >= 0.0_q2) EXIT
      px=px+rngxf
    END DO
    DO
      IF(px < rngxf) EXIT
      px=px-rngxf
    END DO
    DO
      IF(py >= 0.0_q2) EXIT
      py=py+rngyf
    END DO
    DO
      IF(py < rngyf) EXIT
      py=py-rngyf
    END DO
    DO
      IF(pz >= 0.0_q2) EXIT
      pz=pz+rngzf
    END DO
    DO
      IF(pz < rngzf) EXIT
      pz=pz-rngzf
    END DO

  RETURN
  END SUBROUTINE pbcq2

!-----------------------------------------------------------------------------------!

  SUBROUTINE verticies(px,py,pz,ngxf,ngyf,ngzf,ix,iy,iz,ixm,iym,izm)
    REAL(q2),INTENT(INOUT) :: px,py,pz
    INTEGER,INTENT(IN) :: ngxf,ngyf,ngzf
    INTEGER,INTENT(OUT) :: ix,iy,iz,ixm,iym,izm

    ix=FLOOR(px)
    IF (ix == 0) THEN
      ix=ngxf
      ixm=1
    ELSE
      ixm=ix+1
    END IF
    iy=FLOOR(py)
    IF (iy == 0) THEN
      iy=ngyf
      iym=1
    ELSE
      iym=iy+1
    END IF
    iz=FLOOR(pz)
    IF (iz == 0) THEN
      iz=ngzf
      izm=1
    ELSE
      izm=iz+1
    END IF

  RETURN
  END SUBROUTINE verticies

!-----------------------------------------------------------------------------------!

  FUNCTION volume(h)
    REAL(q2),INTENT(IN),DIMENSION(3,3) :: h
    REAL(q2) :: volume

    volume = h(1,1)*(h(2,2)*h(3,3)-h(2,3)*h(3,2))   &
    &       -h(1,2)*(h(2,1)*h(3,3)-h(3,1)*h(2,3))   &
    &       +h(1,3)*(h(2,1)*h(3,2)-h(3,1)*h(2,2))

  RETURN
  END FUNCTION volume

!-----------------------------------------------------------------------------------!

  SUBROUTINE charge2atom(bader_charge,bdim,Rdir,ndim,bader_achg,bader_dist,         &
  &                      bader_atom,lattice)
    INTEGER,INTENT(IN) :: bdim,ndim
    REAL(q2),INTENT(IN),DIMENSION(bdim,4) :: bader_charge
    REAL(q2),INTENT(IN),DIMENSION(ndim,3) :: Rdir
    REAL(q2),INTENT(OUT),DIMENSION(ndim) :: bader_achg
    REAL(q2),INTENT(OUT),DIMENSION(bdim) :: bader_dist
    INTEGER,INTENT(OUT),DIMENSION(bdim) :: bader_atom
    REAL(q2),INTENT(IN),DIMENSION(3,3) :: lattice

    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: dv,v
    REAL(q2) :: dsq,dminsq
    INTEGER :: i,j,dindex

    bader_achg=0.0_q2
    CALL transform_matrix(lattice,B,3,3)
    DO i=1,bdim
      dv=bader_charge(i,1:3)-Rdir(1,:)
      CALL dpbc_dir(dv)
      CALL matrix_vector(B,dv,v,3,3)
      dminsq=DOT_PRODUCT(v,v)
      dindex=1
      DO j=2,ndim
        dv=bader_charge(i,1:3)-Rdir(j,:)
        CALL dpbc_dir(dv)
        CALL matrix_vector(B,dv,v,3,3)
        dsq=DOT_PRODUCT(v,v)
        IF (dsq < dminsq) THEN
          dminsq=dsq
          dindex=j
        END IF
      END DO
      bader_dist(i)=SQRT(dminsq)
      bader_atom(i)=dindex
      bader_achg(dindex)=bader_achg(dindex)+bader_charge(i,4)
    END DO

  RETURN
  END SUBROUTINE charge2atom

!-----------------------------------------------------------------------------------!

END MODULE ChargeUtil
