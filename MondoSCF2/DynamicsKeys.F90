MODULE DynamicsKeys
  CHARACTER(LEN=*),  PARAMETER :: MOLDYNE          ='Dynamics'
  ! const. E/Vol velocity verlet
  CHARACTER(LEN=*),  PARAMETER :: MD_VV            ='VVerlet' 
  INTEGER,           PARAMETER :: MD_SERIAL_VERLET =34235423
  ! Do parallel tempering (requires parallel code)
  CHARACTER(LEN=*),  PARAMETER :: MD_PTEMPER       ='PTemper'
  INTEGER,           PARAMETER :: MD_PARALLEL_REP  =21532454
  CHARACTER(LEN=*),  PARAMETER :: MD_ATOMWRAP      ='AtomWrap'
  CHARACTER(LEN=*),  PARAMETER :: MD_CLOBBER       ='Clobber'
  ! Number of replicas to use (requires parallel code)
  CHARACTER(LEN=*),  PARAMETER :: MD_REPLICAS      ='Replicas'      
  ! restart input 
  CHARACTER(LEN=*),  PARAMETER :: MD_INPUTS        ='MDinputs'      
  ! integration timestep
  CHARACTER(LEN=*),  PARAMETER :: MD_TIME_STEP     ='dT'            
  ! freq. to write coords
  CHARACTER(LEN=*),  PARAMETER :: MD_CRDFREQ       ='CRDfreq'       
  ! freq. to write energies
  CHARACTER(LEN=*),  PARAMETER :: MD_ENEFREQ       ='ENEfreq'       
  ! freq. to write velocities
  CHARACTER(LEN=*),  PARAMETER :: MD_VELFREQ       ='VELfreq'       
  ! freq. to write restart
  CHARACTER(LEN=*),  PARAMETER :: MD_RESFREQ       ='RESfreq'       
  ! starting temp. from Boltzmann distribution (0=don't assign)
  CHARACTER(LEN=*),  PARAMETER :: MD_TEMP0         ='TEMP0'              
  ! set temperature (KELVIN)
  CHARACTER(LEN=*),  PARAMETER :: MD_TEMP          ='TEMP'               
  ! set pressure (ATM)
  CHARACTER(LEN=*),  PARAMETER :: MD_PRES          ='PRES'               
  ! temperature coupling time
  CHARACTER(LEN=*),  PARAMETER :: MD_TTAU          ='TTAU'               
  ! pressure coupling time
  CHARACTER(LEN=*),  PARAMETER :: MD_PTAU          ='PTAU'               
  ! periodically remove CM translation
  CHARACTER(LEN=*), PARAMETER :: MD_REM_TRANS     ='RemoveTranslation'  
  ! periodically remove uniform rotation
  CHARACTER(LEN=*), PARAMETER :: MD_REM_ROTAT     ='RemoveRotation'     
  ! frequency to remove CM motion
  CHARACTER(LEN=*),  PARAMETER :: MD_TRANS_FREQ    ='TRANSfreq'          
  ! frequency to remove rotation
  CHARACTER(LEN=*),  PARAMETER :: MD_ROTAT_FREQ    ='ROTATfreq'          
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
