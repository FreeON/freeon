! Scratch space for integral intermediates
MODULE VScratchB
  USE DerivedTypes
  INTEGER,PARAMETER :: MaxV=1000
  REAL(DOUBLE),DIMENSION(MaxV) :: V
  REAL(DOUBLE) :: ABX,ABY,ABZ,CDX,CDY,CDZ,QCX,QCY,QCZ,R1XZPE,HFXZPE, &
                  R1X2E,R1X2Z,EXZPE,ZXZPE,PAX,PAY,PAZ,WPX,WPY,WPZ,WQX, &
                  WQY,WQZ,SpSpB,FnSpB,SpFnB,SpSpK,FnSpK,SpFnK
END MODULE VScratchB
