MODULE BasisSetParameters
   USE GlobalCharacters
   IMPLICIT NONE
   INTEGER, PARAMETER :: MaxKinds=20
   INTEGER, PARAMETER :: MaxCntrx=20
   INTEGER, PARAMETER :: MaxPrmtv=20
   INTEGER, PARAMETER :: MaxASymt= 4
   INTEGER, PARAMETER :: MaxLTyps=10
   INTEGER, PARAMETER :: NSupSets=31

   ! Element 105 is a ghost function with charge = 0. 
  
   CHARACTER(LEN=2),DIMENSION(105):: Ats=(/'h ','he','li','be','b ','c ', &
        'n ','o ','f ','ne','na','mg','al','si','p ','s ','cl','ar','k ', &
        'ca','sc','ti','v ','cr','mn','fe','co','ni','cu','zn','ga','ge', & 
        'as','se','br','kr','rb','sr','y ','zr','nb','mo','tc','ru','rh', &
        'pd','ag','cd','in','sn','sb','te','i ','xe','cs','ba','la','ce', &
        'pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu', &
        'hf','ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po', &
        'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm','bk', &
        'cf','es','fm','md','no','lr','ky','gh'/)  
!
   REAL(DOUBLE),DIMENSION(105)    :: AtsMss = (/ 1.007900D0,    4.0026D0,      6.941D0,     &  
    9.01218D0,    10.81D0,       12.011D0,      14.00670D0,    15.9994D0,     18.998403D0,  & 
   20.179D0,      22.98977D0,    24.305D0,      26.98154D0,    28.0855D0,     30.97376D0,   & 
   32.06D0,       35.453D0,      39.948D0,      39.0983D0,     40.08D0,       44.9559D0,    & 
   47.9D0,        50.9415D0,     51.996D0,      54.938D0,      55.847D0,      58.9332D0,    & 
   58.7D0,        63.546D0,      65.38D0,       69.72D0,       72.59D0,       74.9216D0,    & 
   78.96D0,       79.904D0,      83.8D0,        85.4678D0,     87.62D0,       88.9059D0,    &
   91.22D0,       92.9064D0,     95.94D0,       98.0D0,       101.07D0,      102.9055D0,    &
  106.4D0,       107.868D0,     112.41D0,      114.82D0,      118.69D0,      121.75D0,      &
  127.6D0,       126.9045D0,    131.3D0,       132.9054D0,    137.33D0,      138.9055D0,    &
  140.12D0,      140.9077D0,    144.24D0,      145.0D0,       150.4D0,       151.96D0,      &
  157.25D0,      158.9254D0,    162.5D0,       164.9304D0,    167.26D0,      168.9342D0,    &
  173.04D0,      174.967D0,     178.49D0,      180.9479D0,    183.85D0,      186.207D0,     &
  190.2D0,       192.22D0,      195.09D0,      196.9665D0,    200.59D0,      204.37D0,      &
  207.2D0,       208.9804D0,    209.0D0,       210.0D0,       222.0D0,       223.0197D0,    &
  226.0254D0,    227.0278D0,    232.03804D0,   231.03588D0,   238.029D0,     237.0482D0,    &
  244.0D0,       243.0D0,       247.0D0,       247.0D0,       251.0D0,       252.0D0,       &
  257.0D0,       258.0D0,       259.0D0,       260.0D0,       261.0D0,       0.0D0/)
!
   CHARACTER(LEN=BASESET_CHR_LEN), DIMENSION(2,NSupSets)     :: CSets=RESHAPE( (/           &
!               ALL ELECTRON BASIS SETS
                'sto-2g          ' , 'STO-2G          ',   &  ! 1
                'sto-3g          ' , 'STO-3G          ',   &  ! 2
                'sto-6g          ' , 'STO-6G          ',   &  ! 3
                '3-21g           ' , '3-21G           ',   &  ! 4
                '3-21gsp         ' , '3-21GSP         ',   &  ! 5 
                '3-21g*          ' , '3-21Gs          ',   &  ! 6
                '4-22gsp         ' , '4-22GSP         ',   &  ! 7
                '6-31g           ' , '6-31G           ',   &  ! 8
                '6-31g*          ' , '6-31Gs          ',   &  ! 9
                '6-31g**         ' , '6-31Gss         ',   &  ! 10
                '6-31++g**       ' , '6-31ppGss       ',   &  ! 11
                '6-311g**        ' , '6-311Gss        ',   &  ! 12
                '6-311gbb        ' , '6-311G3df,3pd   ',   &  ! 13
                '6-311++gbb      ' , '6-311ppG3df,3pd ',   &  ! 14
                'mini            ' , 'MINI            ',   &  ! 15
                'midi            ' , 'MIDI            ',   &  ! 16
                'dunninghay_sv   ' , 'DunningHay_SV   ',   &  ! 17
                'dunninghay_svp  ' , 'DunningHay_SVP  ',   &  ! 18
                'dunning_dz      ' , 'Dunning_DZ      ',   &  ! 19
                'dunning_dzp     ' , 'Dunning_DZP     ',   &  ! 20
                'dunning_tz      ' , 'Dunning_TZ      ',   &  ! 21
                'ahlrichs_vdz    ' , 'Ahlrichs_VDZ    ',   &  ! 22
                'ahlrichs_pvdz   ' , 'Ahlrichs_pVDZ   ',   &  ! 23
                'ahlrichs_vtz    ' , 'Ahlrichs_VTZ    ',   &  ! 24
                'ahlrichs_tzv    ' , 'Ahlrichs_TZV    ',   &  ! 25
                'wtbs            ' , 'WTBS            ',   &  ! 26
!               ECP BASIS SETS
                'sbkjc           ' , 'SBKJC_VDZ_ECP   ',   &  ! 27
                'lanl2           ' , 'LANL2_DZ_ECP    ',   &  ! 28
!               USER DEFINED BASIS SETS
                'user1           ' , 'User1           ',   &  ! 29
                'user2           ' , 'User2           ',   &  ! 30
!               Crystal Basis Set
                'Crystal         ' , 'Crystal98       '    &  ! 31
                /), (/2,NSupSets/) )
   CHARACTER(LEN=5),DIMENSION(MaxLTyps) ::  CLTyps = &
         (/'s    ','p    ','d    ','f    ','sp   ', &
           'pd   ','spd  ','df   ','pdf  ','spdf '/)
   INTEGER,DIMENSION(2,MaxLTyps) ::  LTyps = &
         RESHAPE( (/0,0,1,1,2,2,3,3,0,1 ,1,2 ,0,2  ,2,3 ,1,3  ,0,3/),(/2,MaxLTyps/))
   CHARACTER(LEN=5),DIMENSION(10)::ASymmTyps = &
         (/ 'S    ', &
            'Px   ','Py   ','Pz   ', &
            'Dxx  ','Dxy  ','Dyy  ','Dxz  ','Dyz  ','Dzz  '/)
END MODULE

