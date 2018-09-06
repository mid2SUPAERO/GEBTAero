!******************************************************
!*                                                    *
!* This module handles I/O of the program including   *
!* definition of inputs/outputs, reading inputs from  *
!* a data file, echo inputs, and write outputs        *
!* to data files. To interface with outside           *
!* environment, one needs to provide such             *
!* capability from the outside environment.           *
!*                                                    *
!*             Copyright (c) Wenbin Yu                *
!*         All rights reserved                        *
!*     Mechanical and Aerospace Engineering           *
!*                                                    *
!*           Utah State University                    *
!******************************************************

MODULE IOaero

USE GlobalDataFun 
USE PrescribedCondition
USE TimeFunctionModule
USE InternalData

IMPLICIT NONE

PRIVATE                              ! So everything is private, except declared by public

! functions called by the main program
!---------------------------------------------
PUBLIC Input,WriteError,Output,OutputVtk             

! variables needed by the main program
!----------------------------------------------
PUBLIC error,nkp,nelem,ndof_el,coord,member,pt_condition,material,aerodyn_coef,sol_pt,sol_mb, &
     & frame,mb_condition,distr_fun,curvature,omega_a0,omega_a_tf,v_root_a0,v_root_a_tf,time_function,simu_time, &
	 & EIN,niter,nstep,ncond_mb,ntimefun,analysis_flag,init_cond,nev,eigen_val,eigen_vec_pt,eigen_vec_mb,&
	 & nmemb,ncond_pt,nmate,nframe,ndistrfun,ncurv, aero_flag, grav_flag,nvtk,Velocity_str

! files needed/generated 
!--------------------------------------------------------------------------------
INTEGER,PARAMETER,PRIVATE:: CHAR_LEN=256
INTEGER,PARAMETER:: IN  =10     ! input file: inp_name
CHARACTER(CHAR_LEN)    :: inp_name 

INTEGER,PARAMETER:: EIN =20     ! file for echoing the inputs: inp_name.ech
CHARACTER(CHAR_LEN+3)    :: ech_name 

INTEGER,PARAMETER:: OUT =40     ! file for output: inp_name.out
CHARACTER(CHAR_LEN+3)    :: out_name 
!--------------------------------------------------------------------------------
INTEGER,PARAMETER:: INIT =50     ! file for initial conditions: inp_name.ini
CHARACTER(CHAR_LEN+3)    :: init_name 
!--------------------------------------------------------------------------------


!Private constants
!-----------------------------------------------


! Private variables
!---------------------------------------------------------------------------
INTEGER::nkp           ! number of key points
INTEGER::nelem         ! total number of elements
INTEGER::nmemb         ! number of members
INTEGER::nmate         ! number of cross-sectional properties sets
INTEGER::nframe        ! number of frames 
INTEGER::ncond_pt      ! number of point conditions for concentrated loads and boundary conditions
INTEGER::ndistrfun     ! number of distributed functions 
INTEGER::ncurv         ! number of initial curvatures/twists
INTEGER::analysis_flag ! 0: static analysis; 1: steady state response; 2: transient analysis; 3: eigenvalue analysis 
INTEGER::nev           ! number of frequencies and modeshapes.
INTEGER::aero_flag     ! 0: no aero analasys; 1: stationary aerodynamic, 2: unsteady aerodynamic
INTEGER::grav_flag     ! 0: without gravity; 1: with gravity


!Public integer variables
!---------------------------------------
INTEGER::ncond_mb      ! number of member conditions for distributed loads
INTEGER::ntimefun      ! number of time functions 
INTEGER::niter         ! number of maximum iterations
INTEGER::nstep         ! number of time steps/load steps
INTEGER::nvtk         ! number of the aerodynamic cycle
INTEGER,ALLOCATABLE::member(:,:)   ! member property array: member(nmemb,MEMB_CONST)
INTEGER::ndof_el       ! dofs per element: 12 for static analysis, 18 for dynamic analysis
INTEGER::omega_a_tf(NDIM)    ! time function numbers for the angular velocity of frame a 
INTEGER::v_root_a_tf(NDIM)   ! time function numbers for the velocity of the starting point of the first member


!Public real variables
!--------------------------------------------------------------------------------
REAL(DBL),ALLOCATABLE:: coord(:,:)       ! nodal coordinates: coord(nkp,NDIM)
REAL(DBL),ALLOCATABLE:: material(:,:,:)  ! flexibility matrix: (nmate,12,6)
REAL(DBL),ALLOCATABLE:: aerodyn_coef(:,:)! 2D aerodynamic coefficient : (nmate,:)
REAL(DBL),ALLOCATABLE:: frame(:,:,:)     ! member frames: (nframe,3,3)
REAL(DBL),ALLOCATABLE:: distr_fun(:,:)   ! prescribed functions: (ndistrfun,6)
REAL(DBL),ALLOCATABLE:: curvature(:,:)   ! curvatures: (ncurv,NDIM)
REAL(DBL),ALLOCATABLE:: sol_pt(:,:,:)    ! solutions for points sol_pt(nstep,nkp,NDIM+NDOF_ND)
REAL(DBL),ALLOCATABLE:: sol_mb(:,:,:)    ! solutions for member sol_mb(nstep,nelem,NDIM+ndof_el): nelem: total number of elements
REAL(DBL)            :: simu_time(2)     ! start and end time of the simulation. 
REAL(DBL)            :: omega_a0(NDIM)   ! the magnitude of angular velocity of frame a 
REAL(DBL)            :: v_root_a0(NDIM)  ! the magnitude of linear velocity of the starting point of the first member
REAL(DBL),ALLOCATABLE:: init_cond(:,:)   ! initial conditions: init_cond(nelem,12);
                                         ! init_cond(nelem,1:6) for initial displacements/rotations 
									     ! init_cond(nelem,7:12) for initial velocities
									     ! init_cond(nelem,13:12+NSTATES) Peters finite state parameter at time t+dt
REAL(DBL),ALLOCATABLE::eigen_val(:,:),eigen_vec_pt(:,:,:),eigen_vec_mb(:,:,:) ! arrays for holding eigenvalues and eigenvectors

!Public derived types
!--------------------------------------------------------------------------------
TYPE(PrescriInf),ALLOCATABLE::pt_condition(:) ! prescribed information concentrated at nodes
TYPE(PrescriInf),ALLOCATABLE::mb_condition(:) ! prescribed information distributed along beam members
TYPE(TimeFunction),ALLOCATABLE::time_function(:) ! time functions
!adding
CHARACTER(CHAR_LEN)  :: Velocity_str = ''
INTEGER::arpack,eigenoutput


!Public character variables
!============================================================================
CHARACTER(300)::error         ! a character variable holding  error message
!=========================================================================================

CONTAINS
!=============================


!*************************************************************
!*                                                           *
!* To read and echo problem definition data                  * 
!*															 *
!*************************************************************
SUBROUTINE Input

INTEGER:: i,j,narg,cptArg,length,ind,arpack
INTEGER::tmp_no ! a temporary integer
CHARACTER(CHAR_LEN):: arg,nb_modes_str='',aero_flag_str='',alpha_ac_str='',beta_ac_str='',analysis_str=''
CHARACTER(CHAR_LEN):: niter_str='',nstep_str='',time_str='',nvtk_str='',vz_str='',wy_str='',tfe_str='',tfper_str=''
CHARACTER(CHAR_LEN):: solver_str='',arpack_str='',eigenoutput_str=''
!~ REAL(DBL)          :: Velocity,nb_modes,alpha_ac,beta_ac


! Get all the command argument
narg = command_argument_count()
IF (narg>0) THEN
	DO cptArg = 1, narg
		CALL get_command_argument(cptArg,arg)
		arg=trim(arg)
		length = len_trim(arg)
		ind = index(arg,"=")
		! input_file name
		IF (arg(max(1,length-3):length)==".dat") THEN
			inp_name=arg
		!parameter (velocity,rho,...)
		ELSEIF (ind>0) THEN
		SELECT CASE(arg(1:ind-1))
		CASE("v","V","velocity")
			Velocity_str=arg(ind+1:)
		CASE("m","M","modes")
			nb_modes_str=arg(ind+1:)
		CASE("analysis")
			analysis_str=arg(ind+1:)
		CASE("aero")
			aero_flag_str=arg(ind+1:)
		CASE("alpha_ac")
			alpha_ac_str=arg(ind+1:)
		CASE("beta_ac")
			beta_ac_str=arg(ind+1:)
		CASE("niter")
			niter_str=arg(ind+1:)
		CASE("nstep")
			nstep_str=arg(ind+1:)
		CASE("time")
			time_str=arg(ind+1:)
		CASE("nvtk")
			nvtk_str=arg(ind+1:)
		CASE("vz")
			vz_str=arg(ind+1:)
		CASE("wy")
			wy_str=arg(ind+1:)
		CASE("tfe")
			tfe_str=arg(ind+1:)
		CASE("tfper")
			tfper_str=arg(ind+1:)
		CASE("solver")
			solver_str=arg(ind+1:)
        CASE("arpack")
            arpack_str=arg(ind+1:)
        CASE("eigenoutput")
            eigenoutput_str=arg(ind+1:)
		END SELECT
		!running option
		ELSEIF (arg(1:1)=="-") THEN
		SELECT CASE(arg(2:))
		CASE ("p","python")
			RUNMOD = 1
		CASE ("s","silent")
			RUNMOD = 2
		CASE ("version")
			write(*,*) 'gebtaero version 18.07 developped by French Air Force Academy Research Center'
			RETURN		
		END SELECT
		ENDIF
	ENDDO
ELSE 
	RETURN
ENDIF

! Modify the solver if necessary
IF (solver_str /= '') READ(solver_str,*) SOLVER

IF(TRIM(inp_name)=='') THEN
    error='Please provide an input file name, executing as GEBT input_file_name'
    WRITE(0,*) 'IOaero.f90 : ',error
    ERROR STOP 107
ENDIF
IF(FileOpen(IN, inp_name, 'OLD', 'READ',error))	 RETURN 


! Create a file name for echoing the input data
!--------------------------------------------------
!~ ech_name=TRIM(inp_name) // ".ech" 
ech_name="input.ech" 
IF(FileOpen(EIN,  ech_name,'REPLACE','WRITE',error)) RETURN
CALL TitlePrint(EIN, 'Input file echo of '//inp_name)

! Input and echo analysis control parameters.
!---------------------------------------------------------
READ(IN,*,IOSTAT=in_stat) analysis_flag,aero_flag,grav_flag,niter,nstep, nvtk
IF(IOError('read analysis control parameters',error)) RETURN

CALL TitlePrint(EIN, 'Analysis Control Parameters')
WRITE(EIN,*) "Structural Analysis Type     = ", analysis_flag
WRITE(EIN,*) "Aerodynamic Analysis Type    = ", aero_flag
WRITE(EIN,*) "Gravity forces               = ", grav_flag
WRITE(EIN,*) "Number of Maximum Iterations = ", niter
WRITE(EIN,*) "Number of Time/Load steps    = ", nstep
WRITE(EIN,*) "Number of vtk file           = ", nvtk

!overwrite with console argument

IF (analysis_str /= '') READ(analysis_str,*) analysis_flag
IF (aero_flag_str /= '') READ(aero_flag_str,*) aero_flag
IF (niter_str /= '') READ(niter_str,*) niter
IF (nstep_str /= '') READ(nstep_str,*) nstep
IF (nvtk_str /= '') READ(nvtk_str,*) nvtk
IF (arpack_str /= '') READ(arpack_str,*) arpack
IF (eigenoutput_str /= '') READ(eigenoutput_str,*) eigenoutput

ARPACK_MOD = arpack
EIGEN_OUTPUT = eigenoutput

IF(analysis_flag==0) THEN
	ndof_el=NDOF_ND
ELSE
    ndof_el=18
	READ(IN,*,IOSTAT=in_stat) omega_a0
	IF(IOError('read angular velocity of the global body-attached frame',error)) RETURN
    WRITE(EIN,*)"Angular Velocity of a frame"
    WRITE(EIN,*) '--------------------------------'	
    IF (wy_str /= '') READ(wy_str,*) omega_a0(2)
	CALL WriteVec(EIN,omega_a0)


	READ(IN,*,IOSTAT=in_stat) omega_a_tf
	IF(IOError('read time functions for angular velocity of the global body-attached frame',error)) RETURN
    WRITE(EIN,*)"Time Function for Angular Velocity of a frame"
    WRITE(EIN,*) '--------------------------------'	
	CALL WriteVec(EIN,omega_a_tf)

	
	READ(IN,*,IOSTAT=in_stat) v_root_a0
	IF(IOError('read velocities of starting point of the first member',error)) RETURN
	WRITE(EIN,*)"Initial Velocity of the starting point of the first member"
    WRITE(EIN,*) '--------------------------------'	
    IF (vz_str /= '') READ(vz_str,*) v_root_a0(3)
    CALL WriteVec(EIN,v_root_a0)
    
	READ(IN,*,IOSTAT=in_stat) v_root_a_tf
	IF(IOError('read time function for velocities of starting point of the first member',error)) RETURN
	WRITE(EIN,*)"Time Function for Initial Velocity of the starting point of the first member"
    WRITE(EIN,*) '--------------------------------'	
    CALL WriteVec(EIN,v_root_a_tf)


	IF(analysis_flag==3) THEN
		READ(IN,*,IOSTAT=in_stat) nev
		IF(IOError('read number of frequencies/mode shapes',error)) RETURN
		WRITE(EIN,*) "Number of Frequencies/Mode Shapes= ", nev
		! overwrite with console argument
		IF (nb_modes_str /= '') READ(nb_modes_str,*) nev
	ENDIF
ENDIF

! modification of the length of the system for the Peters aerodynamic model
IF (aero_flag==3) ndof_el = ndof_el+NSTATES

! Input and echo mesh control parameters
!----------------------------------------------------------------
READ(IN,*,IOSTAT=in_stat)  nkp, nmemb, ncond_pt,nmate, nframe,ncond_mb,ndistrfun,ntimefun, &
                        &   ncurv
IF(IOError('read mesh conontrol parameters',error)) RETURN

CALL TitlePrint(EIN, 'Mesh Control Parameters')
WRITE(EIN,*) "Number of Key Points             = ", nkp
WRITE(EIN,*) "Number of Members                = ", nmemb
WRITE(EIN,*) "Number of Point Conditions       = ", ncond_pt
WRITE(EIN,*) "Number of Cross-sections         = ", nmate
WRITE(EIN,*) "Number of Frames                 = ", nframe
WRITE(EIN,*) "Number of Distributed Load Cases = ", ncond_mb
WRITE(EIN,*) "Number of Time Functions         = ", ntimefun
WRITE(EIN,*) "Number of Distributed Functions  = ", ndistrfun
WRITE(EIN,*) "Number of Curvatures/Twist       = ", ncurv

! A small check for the control data
!---------------------------------
IF(nkp<=0 .OR. nmemb<=0 .OR. nmate<=0 .OR. ncond_pt<=0) THEN
   error='nkp, nmemb, nmat, ncond_pt must be greater than 0'
    WRITE(0,*) 'IOaero.f90 : ',error
    ERROR STOP 108
ENDIF  


! Input and echo coordinates for key points
!=============================================================================================
ALLOCATE(coord(nkp,NDIM),STAT=allo_stat)
IF(MemoryError('coord',error)) GOTO 9999
coord=0.0D0

CALL TitlePrint(EIN, 'Key Point Coordinates')

DO i=1,nkp
   READ(IN,*,IOSTAT=in_stat)tmp_no,coord(tmp_no,:)
   IF(IOError('read key point coordinates',error)) GOTO 9999

   WRITE(EIN,*)"Point NO: ",tmp_no
   WRITE(EIN,*) '--------------------------------'	
   CALL WriteVec(EIN,coord(tmp_no,:))
ENDDO   


! Input and echo member properties: we need 7 numbers to describe a member including starting
! and ending points, structural/inertial property set # for the starting and ending points, 
! frame # of the starting point, number of elements divided, geometry property set #
!=============================================================================================
ALLOCATE(member(nmemb,MEMB_CONST),STAT=allo_stat)
IF(MemoryError('member',error)) GOTO 9999
member=0

CALL TitlePrint(EIN, 'Member Definition')
WRITE(EIN,*) '      N1        N2     Sec#1  Sec#2  Frm#   #Elems   Geom#'
WRITE(EIN,*) '----------------------------------------------------------------'	

DO i=1,nmemb  
   READ(IN,*,IOSTAT=in_stat)tmp_no,member(tmp_no,:)
   IF(IOError('read member properties',error)) GOTO 9999
   WRITE(EIN,*)"Member NO: ",tmp_no
   WRITE(EIN,*) '--------------------------------'	
   CALL WriteVec(EIN,member(tmp_no,:))
ENDDO   

! Input and echo point conditions
!=============================================================================================
ALLOCATE(pt_condition(ncond_pt),STAT=allo_stat)
IF(MemoryError('pt_condition',error)) GOTO 9999
pt_condition=InitPI()
 
CALL TitlePrint(EIN, 'Prescribed Point Conditions')
DO i=1,ncond_pt
	pt_condition(i)=InputEchoPrescribedConditions(IN,EIN,error)
ENDDO       


! Input and echo cross-sectional properties including flexibility matrix and mass matrix
!=============================================================================================
ALLOCATE(material(nmate,NSTRN+NSTRN,NSTRN),STAT=allo_stat)
IF(MemoryError('material',error)) GOTO 9999
material=0.0D0

IF(aero_flag/=0.OR.grav_flag/=0) THEN
	ALLOCATE(aerodyn_coef(nmate,8),STAT=allo_stat)
	IF(MemoryError('aerodyn_coef',error)) GOTO 9999
	aerodyn_coef=0.0D0
ENDIF

CALL TitlePrint(EIN, 'Sectional Properties')

DO i=1,nmate
   
   READ(IN,*,IOSTAT=in_stat) tmp_no 
   IF(IOError('read number of sectional properties',error)) GOTO 9999

   WRITE(EIN,*) 'Section No.         =', tmp_no
   WRITE(EIN,*) '--------------------------------'
   CALL TitlePrint(EIN, 'Sectional Flexibility Matrix')
   
   DO j=1,NSTRN
      READ(IN,*,IOSTAT=in_stat)material(tmp_no,j, :) 
	  IF(IOError('read sectional flexibility matrix',error)) GOTO 9999
	  CALL WriteVec(EIN,material(tmp_no,j,:))
   ENDDO

   IF(analysis_flag/=0) THEN
   		CALL TitlePrint(EIN, 'Sectional Mass Matrix')

		DO j=NSTRN+1,NSTRN+NSTRN
			READ(IN,*,IOSTAT=in_stat)material(tmp_no,j, :)  
			IF(IOError('read sectional mass matrix',error)) GOTO 9999
			CALL WriteVec(EIN,material(tmp_no,j,:)) 
		ENDDO
    ENDIF
    
    IF(aero_flag/=0) THEN

   		CALL TitlePrint(EIN, 'Steady Aerodynamic Coefficients')		
		READ(IN,*,IOSTAT=in_stat)aerodyn_coef(tmp_no,1:7)
        !overwrite with console argument
		IF (Velocity_str /= '') READ(Velocity_str,*) aerodyn_coef(tmp_no,1)
		IF (alpha_ac_str /= '') READ(alpha_ac_str,*) aerodyn_coef(tmp_no,5)
		IF (beta_ac_str /= '')  READ(beta_ac_str,*) aerodyn_coef(tmp_no,6)
		WRITE(EIN,*)'Velocity = ',aerodyn_coef(tmp_no,1)
		WRITE(EIN,*)'Rho = ',aerodyn_coef(tmp_no,2)
		WRITE(EIN,*)'chord = ',aerodyn_coef(tmp_no,3)
		WRITE(EIN,*)'a = ',aerodyn_coef(tmp_no,4)
		WRITE(EIN,*)'alpha_ac = ',aerodyn_coef(tmp_no,5)
		WRITE(EIN,*)'beta_ac = ',aerodyn_coef(tmp_no,6)
		WRITE(EIN,*)'Empty = ',aerodyn_coef(tmp_no,7)

		
    ENDIF
    
    IF(grav_flag/=0) THEN
		! X_cg is positf near the trailing edge.
		READ(IN,*,IOSTAT=in_stat)aerodyn_coef(tmp_no,8)
		WRITE(EIN,*)'Xcg = ',aerodyn_coef(tmp_no,8)
    ENDIF

ENDDO       

! Input and echo frame 
!=============================================================================================
IF(nframe>0) THEN
	ALLOCATE(frame(nframe,NDIM,NDIM),STAT=allo_stat)
	IF(MemoryError('frame',error)) GOTO 9999
	frame=0.0D0

	CALL TitlePrint(EIN, 'Member Frames')
	
	DO i=1,nframe
   
		READ(IN,*,IOSTAT=in_stat) tmp_no 
		IF(IOError('read number of member frames',error)) GOTO 9999
   
		WRITE(EIN,*) 'Frame No.         =', tmp_no
		WRITE(EIN,*) '-----------------------------------'	
		
		DO j=1,NDIM
			READ(IN,*,IOSTAT=in_stat)frame(tmp_no,j, :) 
			IF(IOError('read member frames',error)) GOTO 9999
			CALL WriteVec(EIN,frame(tmp_no,j,:)) ! echo
		ENDDO

	ENDDO       
ENDIF

! Input and echo distributed loads
!=============================================================================================
IF (ncond_mb>0) THEN
	ALLOCATE(mb_condition(ncond_mb),STAT=allo_stat)
	IF(MemoryError('mb_condition',error)) GOTO 9999
	mb_condition=InitPI()

	DO i=1,ncond_mb
		CALL TitlePrint(EIN, 'Distributed Loads')
		mb_condition(i)=InputEchoPrescribedConditions(IN,EIN,error)
	ENDDO
ENDIF  

IF(ndistrfun>0) THEN
	ALLOCATE(distr_fun(ndistrfun,NSTRN),STAT=allo_stat)
	IF(MemoryError('distr_fun',error)) GOTO 9999
	distr_fun=0.0D0

    CALL TitlePrint(EIN, 'Distributed Functions Approximated Using Chebychev polynomials')
   		
	DO i=1,ndistrfun
		READ(IN,*,IOSTAT=in_stat) tmp_no
	    IF(IOError('read function number for Chebychev polynomials',error)) GOTO 9999
		
		READ(IN,*,IOSTAT=in_stat)distr_fun(tmp_no, :) 
		IF(IOError('read functions for distributed load',error)) GOTO 9999
	
		WRITE(EIN,*) 'Function No.=', tmp_no
		WRITE(EIN,*) '------------------------------------'	
		CALL WriteVec(EIN,distr_fun(tmp_no,:))
    ENDDO

ENDIF     

! Input and echo initial curvatures/twist
!=============================================================================================
IF (ncurv>0) THEN
	ALLOCATE(curvature(ncurv,NDIM),STAT=allo_stat)
	IF(MemoryError('curvature',error)) GOTO 9999
	curvature=0.0D0

    CALL TitlePrint(EIN, 'Initial Curvatures/Twist')
	
	DO i=1,ncurv
		
		READ(IN,*,IOSTAT=in_stat) tmp_no
		IF(IOError('read # for initial curvatures/twist',error)) GOTO 9999
		 
	    READ(IN,*,IOSTAT=in_stat)curvature(tmp_no, :) 
	    IF(IOError('read initial curvatures/twist',error)) GOTO 9999
	
		WRITE(EIN,*) 'Case No.=', tmp_no
		WRITE(EIN,*) '---------------------------------------'	
		CALL WriteVec(EIN,curvature(tmp_no,:))
   ENDDO
   
ENDIF       


! Input and echo time functions
!=============================================================================================
IF (ntimefun>0.OR.analysis_flag==2) THEN
	WRITE(EIN,*) 
    
	READ(IN,*,IOSTAT=in_stat)simu_time
	IF(IOError('read simulation range',error)) GOTO 9999
	WRITE(EIN,'(A13,'//FMT_REAL//')') "Start Time =", simu_time(1)
    WRITE(EIN,'(A13,'//FMT_REAL//')') "End Time   =", simu_time(2)
    IF (time_str /= '')  READ(time_str,*) simu_time(2)
ENDIF

!~ simu_time(2) = simu_time(2)/(1.1**restart_nb)! modifing nstep depending on the number of restart

IF (ntimefun>0) THEN
	ALLOCATE(time_function(ntimefun),STAT=allo_stat)
	IF(MemoryError('time_function',error)) GOTO 9999
	time_function=InitTF()
   
	DO i=1,ntimefun
	
		READ(IN,*,IOSTAT=in_stat) tmp_no 
		IF(IOError('read number of time functions',error)) GOTO 9999
		IF(tmp_no<1) THEN
			error="time function number cannot be less than 1"
            WRITE(0,*) 'IOaero.f90 : ',error
            ERROR STOP 108
		ENDIF
		WRITE(EIN,*) 'Time Function No.         =', tmp_no
		WRITE(EIN,*) '-----------------------------------'	
		time_function(tmp_no)=InputEchoTimeFunctions(IN,EIN,error)
	    
	    IF (tfe_str /= '')  READ(tfe_str,*) time_function(tmp_no)%te
	    IF (tfper_str /= ''.AND. time_function(tmp_no)%fun_type == 1)  READ(tfper_str,*) time_function(tmp_no)%fun_val(1)

	ENDDO
ENDIF

! For time marching, we also need to read initial conditions, could be provided by 
! users or automatically generatedly from steady state solution
!===============================================================
nelem=SUM(member(:,6))  ! total number of elements for the structure
IF(analysis_flag==2) THEN
	ALLOCATE(init_cond(nelem,12+NSTATES),STAT=allo_stat)  ! storing the initial coniditions
	IF(MemoryError('init_cond',error)) GOTO 9999
	init_cond=0.D0
	! Create a file name for reading initial data
	!--------------------------------------------------
	init_name=TRIM(inp_name) // ".ini" 
	IF(FileOpen(INIT, init_name, 'OLD', 'READ',error))	 RETURN 
	
	WRITE(EIN,*) 
	WRITE(EIN,*) 'Initial conditions'
	WRITE(EIN,*) '-----------------------------------'	
	
	DO i=1,nelem
		READ(INIT,*,IOSTAT=in_stat)init_cond(i, 1:6) 
		IF(IOError('read initial displacements/rotations',error)) GOTO 9999
		CALL WriteVec(EIN,init_cond(i, 1:NSTRN))
	ENDDO

	DO i=1,nelem
		READ(INIT,*,IOSTAT=in_stat)init_cond(i, 7:12) 
		IF(IOError('read initial velocities',error)) GOTO 9999
        CALL WriteVec(EIN,init_cond(i, 7:12))
        
        IF(aero_flag /=0) THEN 	
			init_cond(i,13:) = 0.D0
        ENDIF
	ENDDO		
	
	CLOSE(INIT)
ENDIF


! Creat arrays holding the outputs
!===============================================================
ALLOCATE(sol_pt(nstep,nkp,NDIM+NDOF_ND),STAT=allo_stat)  ! for each point there are 12 variables
IF(MemoryError('sol_pt',error)) GOTO 9999
ALLOCATE(sol_mb(nstep,nelem,2*NDIM+ndof_el+9),STAT=allo_stat) ! for each element, there are 18 variables.
IF(MemoryError('sol_mb',error)) GOTO 9999

IF(analysis_flag==3) THEN
	ALLOCATE(eigen_val(2,nev+1),STAT=allo_stat)
	IF(MemoryError('eigen_val',error)) GOTO 9999
	ALLOCATE(eigen_vec_pt(nev+1,nkp,NDIM+NDOF_ND),STAT=allo_stat)
	IF(MemoryError('eigen_vec_pt',error)) GOTO 9999
	ALLOCATE(eigen_vec_mb(nev+1,nelem,NDIM+ndof_el),STAT=allo_stat)
	IF(MemoryError('eigen_vec_mb',error)) GOTO 9999
ENDIF

IF (RUNMOD == 0) WRITE(*,*) 'The inputs are echoed in ',TRIM(ech_name)
CLOSE(IN)

9999 IF(error/='')THEN
        WRITE(0,*) 'IOaero.f90 : ',error
        ERROR STOP 109
     ENDIF

END SUBROUTINE Input
!*****************************************



!****************************************************************
!*                                                              *
!* To output the results for each point and each division       *
!* within each member                                           *
!*==============================================================*
!* x y z                                                        *
!* ux uy uz thetax thetay thetaz                                *
!* Fx Fy Fz Mx My Mz                                            *
!* x y z are the coordinates in the global coordinate system    *
!* for element the coordinate for the mid-point is used         *
!*															    *
!****************************************************************
SUBROUTINE Output
 
INTEGER::istep,ikp,imemb,j,i,idiv,imodes=1,vec_index,vec_incr,complex_flag

INTEGER::ndiv ! number of divisions
INTEGER::div_no ! current division number
REAL(DBL) :: temp(18+NSTATES),LambdaRatio(nelem+1),EigenVec(nsize)

IF (RUNMOD ==0) THEN

    out_name=TRIM(inp_name) // ".out"
    IF(FileOpen(OUT,  out_name,'REPLACE','WRITE',error)) RETURN

    CALL TitlePrint(OUT, 'The Solution of Internal Variables')   

    IF(nvtk>nstep) nvtk = nstep

    IF(analysis_flag==3) THEN
    WRITE(OUT,*)"Eigenvalue #"
        DO istep=1,nev
            CALL WriteVec(OUT,eigen_val(:,istep))
        ENDDO
        
        DO istep=1,nev
       
            WRITE(OUT,*)"Eigenvalue #", istep
            CALL WriteVec(OUT,eigen_val(:,istep))

            DO ikp=1,nkp
                WRITE(OUT,*)"Point #: ",ikp
                WRITE(OUT,*) '--------------------------------'	
                CALL WriteVec(OUT,eigen_vec_pt(istep,ikp,1:3))
                CALL WriteVec(OUT,eigen_vec_pt(istep,ikp,4:9))
                CALL WriteVec(OUT,eigen_vec_pt(istep,ikp,10:15))
                WRITE(OUT,*) 
            ENDDO
       
            div_no=0

            DO imemb=1,nmemb
                WRITE(OUT,*)"Member #: ",imemb
                WRITE(OUT,*) '--------------------------------'
            
                ndiv=member(imemb,6)
            
                DO j=1,ndiv
                    div_no=div_no+1
                    CALL WriteVec(OUT,eigen_vec_mb(istep,div_no,1:3))
                    CALL WriteVec(OUT,eigen_vec_mb(istep,div_no,4:9))
                    CALL WriteVec(OUT,eigen_vec_mb(istep,div_no,10:15))
                    IF(ndof_el>=18) CALL WriteVec(OUT,eigen_vec_mb(istep,div_no,16:21))
                    IF(ndof_el>18) THEN
                        CALL WriteVec(OUT,eigen_vec_mb(istep,div_no,22:))                        
                    ENDIF    
                    WRITE(OUT,*)
                ENDDO
            ENDDO
        ENDDO
        
    ELSE
        DO istep=1,nstep, FLOOR(REAL(nstep/max(nvtk,1)))
       
            IF(nstep/=1) WRITE(OUT,*)"Step #", istep

           DO ikp=1,nkp
              WRITE(OUT,*)"Point #: ",ikp
              WRITE(OUT,*) '--------------------------------'	
              CALL WriteVec(OUT,sol_pt(istep,ikp,1:3))
              CALL WriteVec(OUT,sol_pt(istep,ikp,4:9))
              CALL WriteVec(OUT,sol_pt(istep,ikp,10:15))
              WRITE(OUT,*)
           ENDDO
           
           div_no=0

           DO imemb=1,nmemb
                WRITE(OUT,*)"Member #: ",imemb
                WRITE(OUT,*) '--------------------------------'
                
                ndiv=member(imemb,6)
                
                DO j=1,ndiv
                    div_no=div_no+1
                    CALL WriteVec(OUT,sol_mb(istep,div_no,1:3))
                    CALL WriteVec(OUT,sol_mb(istep,div_no,4:9))
                    CALL WriteVec(OUT,sol_mb(istep,div_no,10:15))
                    IF(ndof_el>=18) CALL WriteVec(OUT,sol_mb(istep,div_no,16:21))
                    WRITE(OUT,*)
                ENDDO
            ENDDO
            WRITE(OUT,*)
        ENDDO
    ENDIF    
    WRITE(*,*) 'The results can be found in ',TRIM(out_name)

    
ELSEIF (RUNMOD==1) THEN  
  
    IF (analysis_flag ==1) THEN
        WRITE(*,*) "*STATIC"
        WRITE(*,'(1x,6ES25.15)') sol_pt(nstep,1,10:15)
        
    ELSEIF (analysis_flag ==3) THEN
        complex_flag = 0
        IF (aero_flag <3) THEN
            vec_incr = 18
        ELSE
            vec_incr = 18+NSTATES
        ENDIF
        

        DO imodes=1,nev
            IF (EIGEN_OUTPUT == 0) THEN
                EigenVec = 0.D0
                div_no = 0
                vec_index = 0
                DO imemb=1,nmemb
                    ndiv=member(imemb,6)
                    DO j=1,ndiv
                        div_no=div_no+1
                        EigenVec(vec_index+1:vec_index+vec_incr) = eigen_vec_mb(imodes,div_no,4:)
                        vec_index = vec_index+vec_incr
                    ENDDO
                ENDDO
            ENDIF    
            
            IF (imodes < nev) THEN
                IF (eigen_val(2,imodes) == -eigen_val(2,imodes+1) .AND. eigen_val(2,imodes)>0.) THEN
                    complex_flag = 1
                    WRITE(*,*) "*complex"
                    WRITE(*,*) eigen_val(:,imodes)
                    IF (EIGEN_OUTPUT==0) THEN
!~                         WRITE(*,*) 0.
                        WRITE(*,*) EigenVec
                    ENDIF    
                ELSE IF (complex_flag ==1) THEN
                    IF (EIGEN_OUTPUT==0) THEN
                        WRITE(*,*) EigenVec
                    ENDIF   
                    complex_flag = 0
                ELSE
                    IF (ABS(eigen_val(1,imodes))>0.D0) THEN
                        WRITE(*,*) "*real"
                        WRITE(*,*) eigen_val(:,imodes)
                        IF (EIGEN_OUTPUT==0) THEN
!~                             WRITE(*,*) 0.
                            WRITE(*,*) EigenVec
                        ENDIF   
                    ENDIF      
                ENDIF
            ELSE IF (complex_flag ==1) THEN
                IF (EIGEN_OUTPUT==0) THEN
                    WRITE(*,*) EigenVec
                ENDIF   
                complex_flag = 0
            ENDIF    
        ENDDO
!~         out_name=TRIM(inp_name) // ".vec"
!~         IF(FileOpen(OUT,  out_name,'REPLACE','WRITE',error)) RETURN
        
!~         ! provisoire correlation
!~         WRITE(*,*) "*EigenVectors"
!~         DO i = 1,nev      
!~             DO imemb = 1,nmemb
!~                 ndiv=member(imemb,6)
!~                 DO idiv = 1,ndiv
!~                     ! normalisation of each component
!~                     eigen_vec_mb(i,idiv,4:6) = eigen_vec_mb(i,idiv,4:6)/Norm(eigen_vec_mb(i,idiv,4:6))
!~                     eigen_vec_mb(i,idiv,7:9) = eigen_vec_mb(i,idiv,7:9)/Norm(eigen_vec_mb(i,idiv,7:9))
!~                     eigen_vec_mb(i,idiv,10:12) = eigen_vec_mb(i,idiv,10:12)/Norm(eigen_vec_mb(i,idiv,10:12))
!~                     eigen_vec_mb(i,idiv,13:15) = eigen_vec_mb(i,idiv,13:15)/Norm(eigen_vec_mb(i,idiv,13:15))
!~                     eigen_vec_mb(i,idiv,16:18) = eigen_vec_mb(i,idiv,16:18)/Norm(eigen_vec_mb(i,idiv,16:18))
!~                     eigen_vec_mb(i,idiv,19:21) = eigen_vec_mb(i,idiv,19:21)/Norm(eigen_vec_mb(i,idiv,19:21))
!~                     IF (aero_flag == 3) THEN
!~                         eigen_vec_mb(i,idiv,22:) = eigen_vec_mb(i,idiv,22:)/Norm(eigen_vec_mb(i,idiv,22:))
!~                     ENDIF
!~                     eigen_vec_mb(i,idiv,4:) = eigen_vec_mb(i,idiv,4:)/Norm(eigen_vec_mb(i,idiv,4:))    
!~                     WRITE(*,*) eigen_vec_mb(i,idiv,4:)
!~                 ENDDO
!~                 WRITE(*,*) "ENDMODE"
!~             ENDDO    
!~         ENDDO
        
        
    ENDIF
ENDIF
CLOSE(OUT)

! Construct given initial conditions
! Assume the initial displacement/rotations is taken from the last step of the 
! steady state solution and the initial velocities are zero, how these values are given
! should be determined by the end user.
!---------------------------------------------------------------------
IF(analysis_flag==1) THEN 
   	init_name=TRIM(inp_name) // ".ini" 
	IF(FileOpen(INIT, init_name, 'REPLACE','WRITE',error))	 RETURN 
   WRITE(INIT,'(1x,6ES25.15)')(sol_mb(nstep,j,4:9),j=1,div_no)
   WRITE(INIT,'(1x,6ES25.15)') ((0.0D0,i=1,NSTRN),j=1,div_no)
   CLOSE(INIT)
ENDIF


END SUBROUTINE Output
!*****************************************

SUBROUTINE OutputVtk
 
INTEGER                :: istep,imemb,j,ielem,compt,nframe,div_no,ndiv
INTEGER                :: iaero ! current iteration in the aero loop
CHARACTER(CHAR_LEN)    :: iaero_str ! iaero converted to string
CHARACTER(CHAR_LEN)    :: ifreq_str ! ith frequency converted to string
CHARACTER(CHAR_LEN)    :: commande_system !
CHARACTER(CHAR_LEN)    :: file_name!
REAL(DBL)              :: norme,temp(NDIM),eCab(NDIM,NDIM),CT(NDIM,NDIM),bw
! airfoil triad variables
REAL(DBL) :: alpha_ac,beta_ac,xB(3),xflow(3),dir_lift(3),dir_moment(3),dir_airfoil(3)


file_name = TRIM(inp_name) //'_vtk'
! create the folder if it doesnt exist
commande_system='mkdir -p ' // TRIM(file_name)
CALL EXECUTE_COMMAND_LINE(commande_system, wait=.TRUE.)

file_name = TRIM(file_name) // '/' //TRIM(inp_name)
!~ IF (Velocity_str/='') THEN
!~ 		file_name = TRIM(file_name) //TRIM(Velocity_str) //'ms'
!~ ENDIF

! cleaning of previous simulation
commande_system='find '// TRIM(file_name) // '* -delete 2>/dev/null'
CALL EXECUTE_COMMAND_LINE(commande_system, wait=.TRUE.)

compt = 1

!~ IF (analysis_flag==3) nstep=1
!~ IF(nvtk>nstep) nvtk = nstep

!~ DO iaero=1,nstep, FLOOR(REAL(nstep/nvtk))
!~ 	! vtk file creation for each aero iteration
!~ 	WRITE(iaero_str,*)compt
!~ 	out_name = TRIM(file_name) // TRIM(ADJUSTL(iaero_str)) // ".vtk"
	
	
!~ 	IF(FileOpen(OUT,  out_name,'REPLACE','WRITE',error)) RETURN

!~ 	WRITE(OUT,'(A)')'# vtk DataFile Version 3.0'
!~ 	WRITE(OUT,'(A)')out_name
!~ 	WRITE(OUT,'(A)')'ASCII'
!~ 	WRITE(OUT,*)' '
!~ 	WRITE(OUT,'(A)')'DATASET POLYDATA'

	
!~ 	IF (analysis_flag==3) THEN 
	
!~ 	! Points coordinates
!~ 	WRITE(OUT,'(A,I6,A)')'POINTS ',nelem,' float'
!~ 	DO ielem=1,nelem
!~ 		WRITE(OUT,*)eigen_vec_mb(1,ielem,1:3)
!~ 	ENDDO

!~ 	! Segments connections
!~ 	WRITE(OUT,*)' '
!~ 	WRITE(OUT,'(A,I4,I6,I6)')'LINES ',1,nelem+1,nelem
!~ 	DO ielem=1,nelem
!~ 		WRITE(OUT,'(I6)')ielem-1
!~ 	ENDDO

!~ 	WRITE(OUT,*)' '
!~ 	WRITE(OUT,'(A,I6)')'POINT_DATA ',nelem
	
!~ 		DO istep=1,nev
!~ 			IF (eigen_val(2,istep)>=-TOLERANCE) THEN
!~ 				compt=1 ! modal displacement
!~ 				WRITE(ifreq_str,'(E16.4)')eigen_val(2,istep)
!~ 				ifreq_str='disp_'//TRIM(ADJUSTL(ifreq_str))//'Hz'
!~ 				WRITE(OUT,'(A,A,A)')'VECTORS ',TRIM(ifreq_str),' float'
!~ 				! normalisation
!~ 				norme = MAXVAL(ABS(eigen_vec_mb(istep,:,4:6)))
								
!~ 				DO ielem=1,nelem
!~ 					WRITE(OUT,*)eigen_vec_mb(istep,ielem,4:6)/norme
!~ 				ENDDO
!~ 				WRITE(OUT,*)' '
				
!~ 				! modal rotation
!~ 				WRITE(ifreq_str,'(E16.4)')eigen_val(2,istep)
!~ 				ifreq_str='rot_'//TRIM(ADJUSTL(ifreq_str))//'Hz'
!~ 				WRITE(OUT,'(A,A,A)')'VECTORS ',TRIM(ifreq_str),' float'
!~ 				! normalisation
!~ 				norme = MAXVAL(ABS(eigen_vec_mb(istep,:,7:9)))
				
!~ 				DO ielem=1,nelem
!~ 					WRITE(OUT,*) eigen_vec_mb(istep,nelem,7:9)/norme
!~ 				ENDDO
!~ 				WRITE(OUT,*)' '

!~ 		! Orientation of the mid chord
!~ 				WRITE(ifreq_str,'(E16.4)')eigen_val(2,istep)
!~ 				ifreq_str='orientation_'//TRIM(ADJUSTL(ifreq_str))//'Hz'
!~ 				WRITE(OUT,'(A,A,A)')'VECTORS ',TRIM(ifreq_str),' float'
!~                 div_no = 0
!~                 DO imemb=1,nmemb
!~                     ndiv=member(imemb,6)
!~                     IF (nframe > 0) THEN
!~                         eCab = frame(imemb,:,:)
!~                     ELSE
!~                         eCab = I3
!~                     ENDIF
!~                     DO j=1,ndiv
!~                         div_no=div_no+1
!~                         CT = DirCosineTRodrigues(eigen_vec_mb(istep,nelem,7:9))
!~                         bw = aerodyn_coef(member(imemb,3),3)
!~                     WRITE(OUT,*) 2*bw*MATMUL(MATMUL(CT,eCab),(/1.D0,.0D0,0.D0/))
!~                     ENDDO
!~                 ENDDO
!~ 				WRITE(OUT,*)' '
!~ 			ENDIF
!~ 		ENDDO		


IF (analysis_flag==3) THEN
    nstep=nev
    nvtk=nev
ENDIF    
IF(nvtk>nstep) nvtk = nstep

DO iaero=1,nstep, FLOOR(REAL(nstep/nvtk))
	! vtk file creation for each aero iteration
	WRITE(iaero_str,*)compt
	out_name = TRIM(file_name) // TRIM(ADJUSTL(iaero_str)) // ".vtk"
	
	
	IF(FileOpen(OUT,  out_name,'REPLACE','WRITE',error)) RETURN

	WRITE(OUT,'(A)')'# vtk DataFile Version 3.0'
	WRITE(OUT,'(A)')out_name
	WRITE(OUT,'(A)')'ASCII'
	WRITE(OUT,*)' '
	WRITE(OUT,'(A)')'DATASET POLYDATA'

	
	IF (analysis_flag==3) THEN 
	
        ! Points coordinates
        WRITE(OUT,'(A,I6,A)')'POINTS ',nelem,' float'
        DO ielem=1,nelem
            WRITE(OUT,*)eigen_vec_mb(iaero,ielem,1:3)
        ENDDO

        ! Segments connections
        WRITE(OUT,*)' '
        WRITE(OUT,'(A,I4,I6,I6)')'LINES ',1,nelem+1,nelem
        DO ielem=1,nelem
            WRITE(OUT,'(I6)')ielem-1
        ENDDO
        
            ! simulation paramaters
        WRITE(OUT,*)' '
        WRITE(OUT,'(A)')'FIELD FieldData 5 '
        WRITE(OUT,*)'Vinf 1 1 float'
        WRITE(OUT,*)aerodyn_coef(1,1)
        WRITE(OUT,*)'Rho 1 1 float'
        WRITE(OUT,*)aerodyn_coef(1,2)
        WRITE(OUT,*)'Chord 1 1 float'
        WRITE(OUT,*)aerodyn_coef(1,3)
        WRITE(OUT,*)'Alpha_AC 1 1 float'
        WRITE(OUT,*)aerodyn_coef(1,5)
        WRITE(OUT,*)'Beta_AC 1 1 float'
        WRITE(OUT,*)aerodyn_coef(1,6)
        
        WRITE(OUT,*)' '
        WRITE(OUT,'(A,I6)')'POINT_DATA ',nelem
        
!~         WRITE(ifreq_str,'(E16.4)')eigen_val(2,istep)
!~         ifreq_str='disp_'//TRIM(ADJUSTL(ifreq_str))//'Hz'
        WRITE(OUT,'(A,A,A)')'VECTORS ','displacements',' float'
        ! normalisation
        norme = MAXVAL(ABS(eigen_vec_mb(iaero,:,4:6)))
                        
        DO ielem=1,nelem
            WRITE(OUT,*)eigen_vec_mb(iaero,ielem,4:6)/norme
        ENDDO
        WRITE(OUT,*)' '
        
        ! modal rotation
!~         WRITE(ifreq_str,'(E16.4)')eigen_val(2,istep)
!~         ifreq_str='rot_'//TRIM(ADJUSTL(ifreq_str))//'Hz'
        WRITE(OUT,'(A,A,A)')'VECTORS ','rotation',' float'
        ! normalisation
        norme = MAXVAL(ABS(eigen_vec_mb(iaero,:,7:9)))
        
        DO ielem=1,nelem
            WRITE(OUT,*) eigen_vec_mb(iaero,ielem,7:9)/norme
        ENDDO
        WRITE(OUT,*)' '

    ! Orientation of the mid chord
!~         WRITE(ifreq_str,'(E16.4)')eigen_val(2,istep)
!~         ifreq_str='orientation_'//TRIM(ADJUSTL(ifreq_str))//'Hz'
        WRITE(OUT,'(A,A,A)')'VECTORS ','midchord_ori',' float'
        div_no = 0
        DO imemb=1,nmemb
            ndiv=member(imemb,6)
            IF (nframe > 0) THEN
                eCab = frame(imemb,:,:)
            ELSE
                eCab = I3
            ENDIF
            DO j=1,ndiv
                div_no=div_no+1
                CT = DirCosineTRodrigues(eigen_vec_mb(iaero,nelem,7:9))
                bw = aerodyn_coef(member(imemb,3),3)
            WRITE(OUT,*) bw*MATMUL(MATMUL(CT,eCab),(/1.D0,.0D0,0.D0/))
            ENDDO
        ENDDO
        WRITE(OUT,*)' '		
	
    
    ELSE 
	! Points coordinates
	WRITE(OUT,'(A,I6,A)')'POINTS ',nelem,' float'
	DO ielem=1,nelem
		WRITE(OUT,*)sol_mb(iaero,ielem,1:3)
	ENDDO

	! Segments connections
	WRITE(OUT,*)' '
	WRITE(OUT,'(A,I4,I6,I6)')'LINES ',1,nelem+1,nelem
	DO ielem=1,nelem
		WRITE(OUT,'(I6)')ielem-1
	ENDDO
	
    ! simulation paramaters
    WRITE(OUT,*)' '
    WRITE(OUT,'(A)')'FIELD FieldData 5 '
    WRITE(OUT,*)'Vinf 1 1 float'
    WRITE(OUT,*)aerodyn_coef(1,1)
    WRITE(OUT,*)'Rho 1 1 float'
    WRITE(OUT,*)aerodyn_coef(1,2)
    WRITE(OUT,*)'Chord 1 1 float'
    WRITE(OUT,*)aerodyn_coef(1,3)
    WRITE(OUT,*)'Alpha_AC 1 1 float'
    WRITE(OUT,*)aerodyn_coef(1,5)
    WRITE(OUT,*)'Beta_AC 1 1 float'
    WRITE(OUT,*)aerodyn_coef(1,6)
    
    
	WRITE(OUT,*)' '
	WRITE(OUT,'(A,I6)')'POINT_DATA ',nelem
	
	! displacement
		WRITE(OUT,'(A,A,A)')'VECTORS ','displacements',' float'
		DO ielem=1,nelem
			WRITE(OUT,*)sol_mb(iaero,ielem,4:6)
		ENDDO
		WRITE(OUT,*)' '
		
		! members rotations
		WRITE(OUT,'(A,A,A)')'VECTORS ','rotations',' float'
		DO ielem=1,nelem
			WRITE(OUT,*)sol_mb(iaero,ielem,7:9)
		ENDDO

		! members forces
		WRITE(OUT,*)' '
		WRITE(OUT,'(A,A,A)')'VECTORS ','forces',' float'
		DO ielem=1,nelem
			WRITE(OUT,*)sol_mb(iaero,ielem,10:12)
		ENDDO

		! members moments
		WRITE(OUT,*)' '
		WRITE(OUT,'(A,A,A)')'VECTORS ','moments',' float'
		DO ielem=1,nelem
			WRITE(OUT,*)sol_mb(iaero,ielem,13:15)
		ENDDO
		
		IF (analysis_flag>0) THEN
			
			! members linear momentum
			WRITE(OUT,*)' '
			WRITE(OUT,'(A,A,A)')'VECTORS ','lin_momentum',' float'
			DO ielem=1,nelem
				WRITE(OUT,*)sol_mb(iaero,ielem,16:18)
			ENDDO
			
			! members angular momentum
			WRITE(OUT,*)' '
			WRITE(OUT,'(A,A,A)')'VECTORS ','ang_momemtum',' float'
			DO ielem=1,nelem
				WRITE(OUT,*)sol_mb(iaero,ielem,19:21)
			ENDDO
			
		ENDIF
		
		IF(aero_flag>0) THEN
			! Orientation of the mid chord
			WRITE(OUT,*)' '
			WRITE(OUT,'(A,A,A)')'VECTORS ','midchord_ori',' float'
			div_no = 0
			DO imemb=1,nmemb
				ndiv=member(imemb,6)
				IF (nframe > 0) THEN
					eCab = frame(imemb,:,:)
				ELSE
					eCab = I3
				ENDIF
				DO j=1,ndiv
					div_no=div_no+1
					CT = DirCosineTRodrigues(sol_mb(iaero,j,7:9))
					bw = aerodyn_coef(member(imemb,3),3)
				WRITE(OUT,*) bw*MATMUL(MATMUL(CT,eCab),(/1.D0,.0D0,0.D0/))
				ENDDO
			ENDDO
			WRITE(OUT,*)' '
			WRITE(OUT,'(A,A,A)')'VECTORS ','airfoil_ori',' float'
			div_no = 0
			DO imemb=1,nmemb
				ndiv=member(imemb,6)
			
				! angles between flow and aircraft coordinate system 
				alpha_ac = aerodyn_coef(member(imemb,3),5)
				beta_ac = aerodyn_coef(member(imemb,3),6)
				
				! wind direction
				xflow = (/COS(alpha_ac)*COS(beta_ac),-COS(alpha_ac)*SIN(beta_ac),-SIN(alpha_ac)/)

				DO j=1,ndiv
					div_no=div_no+1

					IF (member(imemb,5)>0) THEN
						eCab = frame(imemb,:,:)
					ELSE
						eCab = I3
					ENDIF
					CT = DirCosineTRodrigues(sol_mb(iaero,div_no,7:9))

					eCab=MATMUL(CT,eCab)
					xB=eCab(:,1)
					! direction of the airfoil coordinate system (xflow,dir_moment,dir_lift) in the aircraft reference
					dir_lift = CrossProduct(xB,xflow)
					dir_lift = 1/Norm(dir_lift)*dir_lift
					dir_moment = CrossProduct(xflow,dir_lift)
					dir_moment = 1/Norm(dir_moment)*dir_moment
					dir_airfoil = CrossProduct(dir_lift,dir_moment)
					WRITE(OUT,*) dir_lift
				ENDDO
			ENDDO
			WRITE(OUT,*)' '
		ENDIF
	ENDIF


	compt = compt+1

ENDDO


END SUBROUTINE OutputVtk


END MODULE IOaero
!=========================================================
