      REAL(DOUBLE), PARAMETER :: Gamma_Grid  =+2.190625000000000D+01
      REAL(DOUBLE), PARAMETER :: Gamma_Switch=+3.200000000000000D+01
      INTEGER,      PARAMETER :: Gamma_Mesh  =700


      REAL(DOUBLE),DIMENSION(0:700) :: & 
                 F0_0,F0_1,F0_2,F0_3,F0_4, & 
                 F1_0,F1_1,F1_2,F1_3,F1_4, & 
                 F2_0,F2_1,F2_2,F2_3,F2_4, & 
                 F3_0,F3_1,F3_2,F3_3,F3_4, & 
                 F4_0,F4_1,F4_2,F4_3,F4_4, & 
                 F5_0,F5_1,F5_2,F5_3,F5_4, & 
                 F6_0,F6_1,F6_2,F6_3,F6_4, & 
                 F7_0,F7_1,F7_2,F7_3,F7_4, & 
                 F8_0,F8_1,F8_2,F8_3,F8_4, & 
                 F9_0,F9_1,F9_2,F9_3,F9_4, & 
                 F10_0,F10_1,F10_2,F10_3,F10_4, & 
                 F11_0,F11_1,F11_2,F11_3,F11_4
      INTEGER :: IF0,IF1,IF2,IF3,IF4,IF5,IF6,IF7,IF8,IF9,IF10,IF11
!      INCLUDE 'F0.Inc' 
!      INCLUDE 'F1.Inc' 
      INCLUDE 'F2.Inc' 
      INCLUDE 'F3.Inc' 
      INCLUDE 'F4.Inc' 
      INCLUDE 'F5.Inc' 
      INCLUDE 'F6.Inc' 
      INCLUDE 'F7.Inc' 
      INCLUDE 'F8.Inc' 
      INCLUDE 'F9.Inc' 
      INCLUDE 'F10.Inc' 
      INCLUDE 'F11.Inc' 
