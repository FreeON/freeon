MODULE CubeGrid
   USE DerivedTypes
   USE GlobalScalars   
   IMPLICIT NONE
!----------------------------------------------------------------------------------
!  Cubature rules
#ifdef RULE5
   INCLUDE 'MMA/CubeRules/Rule5.Inc'
#endif
#ifdef RULE7
   INCLUDE 'MMA/CubeRules/Rule7.Inc'
#endif
#ifdef RULEB
   INCLUDE 'MMA/CubeRules/RuleB.Inc'
#endif
#ifdef RULEC
   INCLUDE 'MMA/CubeRules/RuleC.Inc'
#endif
#ifdef RULED
   INCLUDE 'MMA/CubeRules/RuleD.Inc'
#endif
#ifdef RULEE
   INCLUDE 'MMA/CubeRules/RuleE.Inc'
#endif
END MODULE 
