SUBROUTINE Int1161(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
 USE DerivedTypes
 USE GlobalScalars
 USE PrettyPrint
 USE InvExp
 INCLUDE "Int1161.Inc"
 REAL(DOUBLE)  :: Cx,Cy,Cz
 REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz,cp1,cp2,cp3
 REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz,cq1,cq2,cq3
 REAL(DOUBLE)  :: EtaQx,EtaQy,EtaQz
 REAL(DOUBLE)  :: r1xZpE,ZpE,Rkk,r1xZ,r1x2Z,rpxZ
 REAL(DOUBLE)  :: PAx,PAy,PAz,PQx,PQy,PQz
 REAL(DOUBLE)  :: QBx,QBy,QBz,WQx,WQy,WQz
 REAL(DOUBLE)  :: Wx,Wy,Wz,WPx,WPy,WPz
 REAL(DOUBLE)  :: T1,T2,T3,T4,TwoT
 REAL(DOUBLE)  :: R1,R2,R3,R4,R5,R6,R7,R8,R9
 REAL(DOUBLE)  :: G0,G1,G2,G3,G4,G5,G6,G7,G8,ET
 REAL(DOUBLE)  :: Ts0,Ts1,r1xE,r1x2E,rpxE,r1x2ZE
 INTEGER       :: I,J,K,ME,MG,l
 INTEGER       :: I0,I1,I2
 IF(N.EQ.0) RETURN
 C=0.0D0
 Cx=DisBufB( 8)
 Cy=DisBufB( 9)
 Cz=DisBufB(10)
 IF(IntCode.EQ. 1010601) THEN
 DO I=1,N
 I0 = SB%SLPrm%I(I)
 I2 = SB%SLDis%I(I)-4
 DO J=1,CKet
 I1=I0+(J-1)*DB%MAXP
 Eta = DB%PrmBuf%D(I1)
 Qx  = DB%PrmBuf%D(I1+1)
 Qy  = DB%PrmBuf%D(I1+2)
 Qz  = DB%PrmBuf%D(I1+3)
 Uq  = DB%PrmBuf%D(I1+4)
 cq1 = DB%PrmBuf%D(I1+5)
 cq2 = DB%PrmBuf%D(I1+6)
 cq3 = DB%PrmBuf%D(I1+7)
 r1xE  = 1.0D0/Eta
 r1x2E = 0.5D0*r1xE
 EtaQx = Eta*Qx
 EtaQy = Eta*Qy
 EtaQz = Eta*Qz
 QBx = Qx-DB%DisBuf%D(I2+7)
 QBy = Qy-DB%DisBuf%D(I2+8)
 QBz = Qz-DB%DisBuf%D(I2+9)
 DO K=1,CBra
 Zeta = PrmBufB(1,K)
 Px   = PrmBufB(2,K)
 Py   = PrmBufB(3,K)
 Pz   = PrmBufB(4,K)
 Up   = PrmBufB(5,K)
 cp1  = PrmBufB(6,K)
 cp2  = PrmBufB(7,K)
 cp3  = PrmBufB(8,K)
 ZpE  = Zeta+Eta
 r1xZpE = 1.0D0/ZpE
 r1x2ZE = 0.5D0*r1xZpE
 r1xZ   = 1.0D0/Zeta
 r1x2Z  = 0.5D0*r1xZ
 rpxZ   = Eta*r1xZpE
 rpxE   = Zeta*r1xZpE
 Rkk  = Up*Uq*DSQRT(r1xZpE)
 PAx = Px-Cx
 PAy = Py-Cy
 PAz = Pz-Cz
 Wx  = (Zeta*Px+EtaQx)*r1xZpE
 Wy  = (Zeta*Py+EtaQy)*r1xZpE
 Wz  = (Zeta*Pz+EtaQz)*r1xZpE
 WPx = Wx-Px
 WPy = Wy-Py
 WPz = Wz-Pz
 WQx = Wx - Qx
 WQy = Wy - Qy
 WQz = Wz - Qz
 PQx = Px - Qx
 PQy = Py - Qy
 PQz = Pz - Qz
 T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
 INCLUDE "RGen2.Inc"
 U(1)=R1                                                          
 U(6)=R2                                                          
 U(5)=R3                                                          
 U(2)=QBx*U(1)+WQx*U(6)                                           
 U(3)=QBy*U(1)+WQy*U(6)                                           
 U(4)=QBz*U(1)+WQz*U(6)                                           
 U(7)=QBx*U(6)+WQx*U(5)                                           
 U(8)=QBy*U(6)+WQy*U(5)                                           
 U(9)=QBz*U(6)+WQz*U(5)                                           
 U(5)=QBx*U(2)+WQx*U(7)+r1x2E*(U(1)-rpxE*U(6))                                             
 U(7)=QBy*U(3)+WQy*U(8)+r1x2E*(U(1)-rpxE*U(6))                                             
 U(10)=QBz*U(4)+WQz*U(9)+r1x2E*(U(1)-rpxE*U(6))                                            
 U(6)=QBx*U(3)+WQx*U(8)                                           
 U(8)=QBx*U(4)+WQx*U(9)                                           
 U(9)=QBy*U(4)+WQy*U(9)                                           
 C(I,1)=C(I,1)+U(5)                                               
 C(I,2)=C(I,2)+U(6)                                               
 C(I,3)=C(I,3)+U(7)                                               
 C(I,4)=C(I,4)+U(8)                                               
 C(I,5)=C(I,5)+U(9)                                               
 C(I,6)=C(I,6)+U(10)                                              
 ENDDO
 ENDDO
 ENDDO
 ELSE
  CALL Halt("Illegal IntCode")
 ENDIF
END SUBROUTINE Int1161
