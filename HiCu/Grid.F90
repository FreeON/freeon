MODULE CubeGrid
   USE DerivedTypes
   USE GlobalScalars   
   IMPLICIT NONE
!----------------------------------------------------------------------------------
!  Cubature rules
#ifdef RULE3
   INCLUDE 'MMA/CubeRules/Rule3.Inc'
#endif
#ifdef RULE5
!  Currently the rule of choice
   INCLUDE 'MMA/CubeRules/Rule5.Inc'
#endif
#ifdef RULE7
   INCLUDE 'MMA/CubeRules/Rule7.Inc'
#endif
#ifdef RULE9
   INCLUDE 'MMA/CubeRules/Rule9.Inc'
#endif
END MODULE 
