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
     DIMENSION(2,NSupSets)     :: CSets=RESHAPE( (/ &
                'sto-2g          ' , 'STO-2G          ',   &
                'sto-3g          ' , 'STO-3G          ',   &
                'sto-6g          ' , 'STO-6G          ',   &
                '3-21g           ' , '3-21G           ',   &
                '3-21gsp         ' , '3-21GSP         ',   &
                '3-21g*          ' , '3-21Gs          ',   &
                '4-22gsp         ' , '4-22GSP         ',   &
                '6-31g           ' , '6-31G           ',   &
                '6-31g*          ' , '6-31Gs          ',   &
                '6-31g**         ' , '6-31Gss         ',   &
                '6-31++g**       ' , '6-31ppGss       ',   &
                '6-311g**        ' , '6-311Gss        ',   &
                'mini            ' , 'MINI            ',   &
                'midi            ' , 'MIDI            ',   &
                'dunninghay_sv   ' , 'DunningHay_SV   ',   &
                'dunninghay_svp  ' , 'DunningHay_SVP  ',   &
                'dunning_dz      ' , 'Dunning_DZ      ',   &
                'dunning_dzp     ' , 'Dunning_DZP     ',   &
                'dunning_tz      ' , 'Dunning_TZ      ',   &
                'ahlrichs_vdz    ' , 'Ahlrichs_VDZ    ',   &
                'ahlrichs_pvdz   ' , 'Ahlrichs_pVDZ   ',   &
                'ahlrichs_vtz    ' , 'Ahlrichs_VTZ    ',   &
                'ahlrichs_tzv    ' , 'Ahlrichs_TZV    ',   &
                'wtbs            ' , 'WTBS            '    &
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

