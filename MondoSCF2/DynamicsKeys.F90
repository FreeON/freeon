MODULE DynamicsKeys
  CHARACTER(LEN=8),  PARAMETER :: MOLDYNE          ='Dynamics'
  ! const. E/Vol velocity verlet
  CHARACTER(LEN=7),  PARAMETER :: MD_VV            ='VVerlet' 
  INTEGER,           PARAMETER :: MD_SERIAL_VERLET =34235423
  ! Do parallel tempering (requires parallel code)
  CHARACTER(LEN=7),  PARAMETER :: MD_PTEMPER       ='PTemper'
  INTEGER,           PARAMETER :: MD_PARALLEL_REP  =21532454
  CHARACTER(LEN=8),  PARAMETER :: MD_ATOMWRAP      ='AtomWrap'
  CHARACTER(LEN=7),  PARAMETER :: MD_CLOBBER       ='Clobber'
  ! Number of replicas to use (requires parallel code)
  CHARACTER(LEN=8),  PARAMETER :: MD_REPLICAS      ='Replicas'      
  ! restart input 
  CHARACTER(LEN=8),  PARAMETER :: MD_INPUTS        ='MDinputs'      
  ! integration timestep
  CHARACTER(LEN=2),  PARAMETER :: MD_TIME_STEP     ='dT'            
  ! freq. to write coords
  CHARACTER(LEN=7),  PARAMETER :: MD_CRDFREQ       ='CRDfreq'       
  ! freq. to write energies
  CHARACTER(LEN=7),  PARAMETER :: MD_ENEFREQ       ='ENEfreq'       
  ! freq. to write velocities
  CHARACTER(LEN=7),  PARAMETER :: MD_VELFREQ       ='VELfreq'       
  ! freq. to write restart
  CHARACTER(LEN=7),  PARAMETER :: MD_RESFREQ       ='RESfreq'       
  ! starting temp. from Boltzmann distribution (0=don't assign)
  CHARACTER(LEN=5),  PARAMETER :: MD_TEMP0         ='TEMP0'              
  ! set temperature (KELVIN)
  CHARACTER(LEN=4),  PARAMETER :: MD_TEMP          ='TEMP'               
  ! set pressure (ATM)
  CHARACTER(LEN=4),  PARAMETER :: MD_PRES          ='PRES'               
  ! temperature coupling time
  CHARACTER(LEN=4),  PARAMETER :: MD_TTAU          ='TTAU'               
  ! pressure coupling time
  CHARACTER(LEN=4),  PARAMETER :: MD_PTAU          ='PTAU'               
  ! periodically remove CM translation
  CHARACTER(LEN=17), PARAMETER :: MD_REM_TRANS     ='RemoveTranslation'  
  ! periodically remove uniform rotation
  CHARACTER(LEN=14), PARAMETER :: MD_REM_ROTAT     ='RemoveRotation'     
  ! frequency to remove CM motion
  CHARACTER(LEN=9),  PARAMETER :: MD_TRANS_FREQ    ='TRANSfreq'          
  ! frequency to remove rotation
  CHARACTER(LEN=9),  PARAMETER :: MD_ROTAT_FREQ    ='ROTATfreq'          
  ! restart input 
  CHARACTER(9),      PARAMETER :: MD_RESTRT_IN     ='RESTRT_IN'          
  ! restart output 
  CHARACTER(10),     PARAMETER :: MD_RESTRT_OUT    ='RESTRT_OUT'         
  ! coordinate output 
  CHARACTER(7),      PARAMETER :: MD_CRD_OUT       ='CRD_OUT'            
  ! velocity output 
  CHARACTER(7),      PARAMETER :: MD_VEL_OUT       ='VEL_OUT'            
  ! energy output  
  CHARACTER(7),      PARAMETER :: MD_ENE_OUT       ='ENE_OUT'            
END MODULE DynamicsKeys
