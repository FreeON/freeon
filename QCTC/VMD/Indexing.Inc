c
      LHGTF(L)=(L+1)*(L+2)*(L+3)/6                                    
      IStride(I,J)=LHGTF(J)-LHGTF(I-1)
      LBegin(L)=(L*(L+1)*(L+2))/6+1
      LEnd(L)=LBegin(L)+L*(L+1)/2+L
      LMNDex(L,M,N)=LBegin(L+M+N)+N*(2*(L+M+N)-N+3)/2+M ! Q-Chem ordering
C
      iCo2x3(l)=l*(l+1)*(l+2)/6                                        
      iCo1x2(l)=l*(l+1)/2                                              
      LMNDexM(l,m,n)=iCo2x3(l+m+n)+iCo1x2(m+n)+n+1        ! MONDO ordering
