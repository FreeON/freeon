MODULE BasisSetParameters
   USE GlobalCharacters
   IMPLICIT NONE
   INTEGER, PARAMETER :: MaxKinds=20
   INTEGER, PARAMETER :: MaxCntrx=20
   INTEGER, PARAMETER :: MaxPrmtv=20
   INTEGER, PARAMETER :: MaxASymt= 4
   INTEGER, PARAMETER :: MaxLTyps=10
   INTEGER, PARAMETER :: NSupSets=24
   CHARACTER(LEN=2),DIMENSION(104):: Ats=(/'h ','he','li','be','b ','c ', &
        'n ','o ','f ','ne','na','mg','al','si','p ','s ','cl','ar','k ', &
        'ca','sc','ti','v ','cr','mn','fe','co','ni','cu','zn','ga','ge', & 
        'as','se','br','kr','rb','sr','y ','zr','nb','mo','tc','ru','rh', &
        'pd','ag','cd','in','sn','sb','te','i ','xe','cs','ba','la','ce', &
        'pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu', &
        'hf','ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po', &
        'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm','bk', &
        'cf','es','fm','md','no','lr','ky'/)
   CHARACTER(LEN=BASESET_CHR_LEN), &
     DIMENSION(2,NSupSets)     :: CSets=RESHAPE( (/        &
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
                'mini            ' , 'MINI            ',   &  ! 13
                'midi            ' , 'MIDI            ',   &  ! 14
                'dunninghay_sv   ' , 'DunningHay_SV   ',   &  ! 15
                'dunninghay_svp  ' , 'DunningHay_SVP  ',   &  ! 16
                'dunning_dz      ' , 'Dunning_DZ      ',   &  ! 17
                'dunning_dzp     ' , 'Dunning_DZP     ',   &  ! 18
                'dunning_tz      ' , 'Dunning_TZ      ',   &  ! 19
                'ahlrichs_vdz    ' , 'Ahlrichs_VDZ    ',   &  ! 20
                'ahlrichs_pvdz   ' , 'Ahlrichs_pVDZ   ',   &  ! 21
                'ahlrichs_vtz    ' , 'Ahlrichs_VTZ    ',   &  ! 22
                'ahlrichs_tzv    ' , 'Ahlrichs_TZV    ',   &  ! 23
                'wtbs            ' , 'WTBS            '    &  ! 24
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

