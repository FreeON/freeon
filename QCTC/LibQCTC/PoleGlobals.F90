MODULE PoleGlobals
  USE Derivedtypes
  USE GlobalScalars   
  IMPLICIT NONE
  INTEGER,PARAMETER                            :: FFELL=64
  INTEGER,PARAMETER                            :: FFELL2=2*FFELL
  INTEGER,PARAMETER                            :: FFLen=FFEll*(FFEll+3)/2 
  INTEGER,PARAMETER                            :: FFLen2=FFEll2*(FFEll2+3)/2!
  REAL(DOUBLE), DIMENSION(0:2*FFEll2)          :: Factorial
  REAL(DOUBLE), DIMENSION(0:FFEll2)            :: FactOlm0,FactMlm0,Sine,Cosine,CoFact
  REAL(DOUBLE), DIMENSION(0:FFLen2)            :: FactOlm2,FactMlm2,ALegendreP,Spq,Cpq
  REAL(DOUBLE), DIMENSION(0:SPEll+1,0:FFELL)   :: FudgeFactorial
!
  INTEGER,PARAMETER                            :: MaxUEll=4
END MODULE
