         DO IGrid=1,NGrid
            RQx=Grid(IGrid,1)-Node%Qx            
            RQy=Grid(IGrid,2)-Node%Qy
            RQz=Grid(IGrid,3)-Node%Qz
            RQ2=RQx*RQx+RQy*RQy+RQz*RQz
            Xpt=EXP(-Node%Zeta*RQ2)
            IF(Xpt>1.D-16)THEN
               TwoZ=Two*Node%Zeta
               LambdaX(0)=Xpt
               LambdaY(0)=One
               LambdaZ(0)=One
               LambdaX(1)=TwoZ*RQx*Xpt
               Lambday(1)=TwoZ*RQy
               LambdaZ(1)=TwoZ*RQz
               DO L=2,Node%Ell+1
                  L1=L-1
                  L2=L-2
                  RL1=DBLE(L1)
                  LambdaX(L)=TwoZ*(RQx*LambdaX(L1)-RL1*LambdaX(L2))
                  LambdaY(L)=TwoZ*(RQy*LambdaY(L1)-RL1*LambdaY(L2))
                  LambdaZ(L)=TwoZ*(RQz*LambdaZ(L1)-RL1*LambdaZ(L2))
               ENDDO
               DO L=0,Node%Ell
                  DO M=0,Node%Ell-L
                     DO N=0,Node%Ell-M-L
                        LMN=LMNDex(L,M,N)
                        RhoV(IGrid,1)=RhoV(IGrid,1)+LambdaX(L  )*LambdaY(M  )*LambdaZ(N  )*Node%Co(LMN)
                        RhoV(IGrid,2)=RhoV(IGrid,2)-LambdaX(L+1)*LambdaY(M  )*LambdaZ(N  )*Node%Co(LMN)
                        RhoV(IGrid,3)=RhoV(IGrid,3)-LambdaX(L  )*LambdaY(M+1)*LambdaZ(N  )*Node%Co(LMN)
                        RhoV(IGrid,4)=RhoV(IGrid,4)-LambdaX(L  )*LambdaY(M  )*LambdaZ(N+1)*Node%Co(LMN)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
