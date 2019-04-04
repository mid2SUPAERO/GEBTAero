!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************

!=========================================================
!
!>  This module contains the linear & nonlinear solver interfaced with MUMPS direct solver library
!
!=========================================================
MODULE SolveMumps

USE InternalData
USE PrescribedCondition
USE System

IMPLICIT NONE

INCLUDE 'mpif.h'
INCLUDE 'dmumps_struc.h'
TYPE (DMUMPS_STRUC) mumps_par   !< an array containing the configuration parameters of MUMPS solver
INTEGER IERR    !< error code of the MUMPS solver

PRIVATE       ! So everything is private, except declared by PUBLIC
PUBLIC LinearSolutionMumps,NewtonRaphsonMumps,CTCabPH,InsertElementValues,ExtractElementValues,ExtractSolution

CONTAINS
!===========================================================


!**************************************************************
!*                                                            *                                      
!>  The linear solver is basically the Newton-Raphson with    *
!!  initial guess equal to zero and only uses one iterations  *
!*															  *
!**************************************************************
SUBROUTINE LinearSolutionMumps(ndof_el,memb_info,v_root_a,omega_a,member,error,&
	& ncond_mb,mb_condition,distr_fun,dof_con,x,aero_flag,grav_flag,init_cond)

INTEGER,INTENT(IN)::ndof_el !<#ioaero::ndof_el
INTEGER,INTENT(IN)::aero_flag   !<#ioaero::aero_flag
INTEGER,INTENT(IN)::grav_flag   !<#ioaero::grav_flag
TYPE (MemberInf),INTENT(IN)::memb_info(:)   ! contains the member parameters of the whole structure

REAL(DBL),INTENT(IN)::v_root_a(:) !< linear velocity of frame a
REAL(DBL),INTENT(IN)::omega_a(:)    !< angular velocity of frame a

REAL(DBL),INTENT(IN) ::distr_fun(:,:)   !<#ioaero::distr_fun
INTEGER,INTENT(IN)   ::member(:,:)  !<#ioaero::member
INTEGER,INTENT(IN)   ::ncond_mb !<#ioaero::ncond_mb
TYPE(PrescriInf),INTENT(IN)::mb_condition(:)    !<#ioaero::mb_condition

INTEGER                 ::dof_con(:) !< this array is passed by value
CHARACTER(*),INTENT(OUT)::error !<#ioaero::error

REAL(DBL),INTENT(OUT) ::x(:)    !<The solution vector

REAL(DBL),OPTIONAL,INTENT(IN)::init_cond(:,:)   !<#ioaero::init_cond
REAL(DBL)  ::rhs(nsize) !<RHS vector

x=0.D0

CALL MPI_INIT(IERR)
! Define a communicator for the package.
mumps_par%COMM = MPI_COMM_WORLD
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
CALL DMUMPS(mumps_par)
!~       output parameters
mumps_par%ICNTL(3) = 0

! Assemble the right hand side
!---------------------------------------------------------
CALL AssembleRHS(ndof_el,memb_info,v_root_a,omega_a,member,error,&
& ncond_mb,mb_condition,distr_fun,(dof_con),x,rhs,aero_flag,grav_flag,init_cond)
IF(error/='') RETURN
IF(MAXVAL(ABS(rhs))<TOLERANCE) RETURN

! Assemble the coefficient matrix
!---------------------------------------------------------
CALL AssembleJacobian(ndof_el,1,memb_info,v_root_a,omega_a,member,error,&
& ncond_mb,mb_condition,distr_fun,(dof_con),x,aero_flag,grav_flag,init_cond)
IF(error/='') RETURN

!  Define problem on the host (processor 0)
IF ( mumps_par%MYID .eq. 0 ) THEN
	mumps_par%N=nsize
	mumps_par%NZ=ne
	ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
	ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
	ALLOCATE( mumps_par%A( mumps_par%NZ ) )
	ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
	mumps_par%IRN(:)=irn(1:mumps_par%NZ)
	mumps_par%JCN(:)=jcn(1:mumps_par%NZ)
	mumps_par%A(:)=coef(1:mumps_par%NZ)
	mumps_par%RHS(:)=rhs(1:mumps_par%N)
	! Analysis and factorisation 
	mumps_par%JOB = 6
	IF(ALLOCATED(coef))DEALLOCATE(coef)
	IF(ALLOCATED(irn)) DEALLOCATE (irn)
	IF(ALLOCATED(jcn)) DEALLOCATE (jcn)
END IF

CALL DMUMPS(mumps_par)

! Solve the linear system
!====================================
!~ CALL  LinearSolver((rhs),x,error)
IF ( mumps_par%MYID .eq. 0 ) THEN
	x(:)=mumps_par%RHS(:)
END IF
!  Deallocate user data
IF ( mumps_par%MYID .eq. 0 )THEN
	DEALLOCATE( mumps_par%IRN )
	DEALLOCATE( mumps_par%JCN )
	DEALLOCATE( mumps_par%A )
	DEALLOCATE( mumps_par%RHS )
END IF
! Destroy the instance (deallocate internal data structures)
mumps_par%JOB = -2
CALL DMUMPS(mumps_par)
CALL MPI_FINALIZE(IERR)

END SUBROUTINE LinearSolutionMumps
!***************************************************************


!************************************************************
!*                                                          *                                      
!>  Use Newton-Raphson method to solve the nonlinear system *
!*															*
!************************************************************ 
SUBROUTINE NewtonRaphsonMumps(ndof_el,memb_info,v_root_a,omega_a,member,niter,error,&
	 & ncond_mb,mb_condition,distr_fun,dof_con,x,aero_flag,grav_flag,init_cond)

INTEGER,INTENT(IN)::ndof_el !<#ioaero::ndof_el
INTEGER,INTENT(IN)::aero_flag   !<#ioaero::aero_flag
INTEGER,INTENT(IN)::grav_flag   !<#ioaero::grav_flag
TYPE (MemberInf),INTENT(IN)::memb_info(:)   !< contains the member parameters of the whole structure

REAL(DBL),INTENT(IN)::v_root_a(:) !< linear velocity of frame a
REAL(DBL),INTENT(IN)::omega_a(:)    !< angular velocity of frame a

REAL(DBL),INTENT(IN) ::distr_fun(:,:)  !<#ioaero::distr_fun
INTEGER,INTENT(IN)   ::member(:,:)  !<#ioaero::member
INTEGER,INTENT(IN)   ::niter    !<#ioaero::niter
INTEGER,INTENT(IN)   ::ncond_mb !<#ioaero::ncond_mb
TYPE(PrescriInf),INTENT(IN):: mb_condition(:)   !<#ioaero::mb_condition
REAL(DBL),INTENT(INOUT) ::x(:)  !<The solution vector

INTEGER                 ::dof_con(:) ! this array is passed by value
CHARACTER(*),INTENT(OUT)::error !<#ioaero::error
REAL(DBL)::rhs(nsize)   !<RHS vector
REAL(DBL),OPTIONAL,INTENT(IN)::init_cond(:,:)   !<#ioaero::init_cond


LOGICAL:: check
REAL(DBL):: fmin,fold !< 0.5*rhs.rhs
REAL(DBL),PARAMETER::STPMX=1000.D0 !< scaled maximum step length allowed in line searches
REAL(DBL)::stpmax
REAL(DBL)::gradient(SIZE(x)) !<-rhs.coef
REAL(DBL)::xold(SIZE(x)) !<holding x from previous iteration
REAL(DBL)::dx(SIZE(x)) !<xnew-xold, the increment calculated by N-R method

REAL(DBL),PARAMETER::TOLF=1.0D-8 !< convergence criterion for the function values

INTEGER::i


CALL AssembleRHS(ndof_el,memb_info,v_root_a,omega_a,member,error,&
& ncond_mb,mb_condition,distr_fun,(dof_con),x,rhs,aero_flag,grav_flag,init_cond)
IF(error/='') RETURN

IF(MAXVAL(ABS(rhs))<TOLF) RETURN
    
fmin=0.5D0*DOT_PRODUCT(rhs,rhs)  ! calculate fnew=1/2*F.F
    
stpmax=STPMX*MAX(Norm(x),REAL(SIZE(x)))   ! Calculate step length for line searches

check=.FALSE.	

CALL MPI_INIT(IERR)
! Define a communicator for the package.
mumps_par%COMM = MPI_COMM_WORLD
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
CALL DMUMPS(mumps_par)
!~       output parameters
mumps_par%ICNTL(3) = 0
mumps_par%N=nsize
IF ( mumps_par%MYID .eq. 0 ) THEN
	ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
	mumps_par%RHS(:)=rhs(1:mumps_par%N)
ENDIF

DO i=1,niter
!~    WRITE(*,*) "ITERATION=",i
   IF(DEBUG) WRITE(IOUT,*) "ITERATION=",i

   ! Assemble the nonlinear system
   !---------------------------------------------------
   CALL AssembleJacobian(ndof_el,niter,memb_info,v_root_a,omega_a,member,error,&
	& ncond_mb,mb_condition,distr_fun,(dof_con),x,aero_flag,grav_flag,init_cond)
   IF(error/='') RETURN
   gradient=MATMUL_sparse(-rhs,nsize,ne,irn,jcn,coef) ! compute gradient for the line search
   xold=x
   fold=fmin
	
	IF ( mumps_par%MYID .eq. 0 ) THEN
		mumps_par%NZ=ne
		ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
		ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
		ALLOCATE( mumps_par%A( mumps_par%NZ ) )
		mumps_par%IRN(:)=irn(1:mumps_par%NZ)
		mumps_par%JCN(:)=jcn(1:mumps_par%NZ)
		mumps_par%A(:)=coef(1:mumps_par%NZ)
		IF(ALLOCATED(coef))DEALLOCATE(coef)
		IF(ALLOCATED(irn)) DEALLOCATE (irn)
		IF(ALLOCATED(jcn)) DEALLOCATE (jcn)
	
		! Analysis and factorisation 		
		mumps_par%JOB = 6
	END IF
	
   ! Solve the linear system
   !====================================
	CALL DMUMPS(mumps_par)
!~    CALL  LinearSolver((rhs),dx,error)
	IF ( mumps_par%MYID .eq. 0 ) THEN
		dx(1:nsize)=mumps_par%RHS(1:nsize)
	END IF
   IF(error/='')	 RETURN
	
	!  Deallocate user data
	IF ( mumps_par%MYID .eq. 0 )THEN
		DEALLOCATE( mumps_par%IRN )
		DEALLOCATE( mumps_par%JCN )
		DEALLOCATE( mumps_par%A )
	END IF
	
   CALL LineSearch  ! use line search to find a new x and rhs
	mumps_par%RHS(1:nsize) = rhs(1:nsize)
   IF(MAXVAL(ABS(rhs))<TOLF) THEN
       check=.FALSE.  ! the function values converge to zero
!~        WRITE(*,*) "ITERATIONS=",i
	   EXIT
   ENDIF 

   IF(check) THEN  ! converges on deltaX
      ! Check for gradient of 1/2 F.F zero, i.e. spurious convergence
	  !---------------------------------------------------------------------
	   check=(MAXVAL(ABS(gradient)*MAX(ABS(x),1.0D0)/MAX(fmin,0.5D0*SIZE(x)))<TOLF) 
	   EXIT
   ENDIF   

!~ ! It seems the following test is needed, if it converges on x and the gradient is not equal to zero.  as if rhs=0, then the following test must be true for a real solution
!~ ! if it is not, then, it must converge to a local minimum, one must restart with different nonlinear parameters
   IF(MAXVAL(ABS(x-xold)/MAX(ABS(x),1.0D0))<=TOLF) THEN
		check=.FALSE.
		RETURN
	ENDIF  ! no more corrections found for x	

    IF(i==niter) THEN
        error='The solution does not converge after the maximum number of iterations' !The maximum number of iterations reached, the solution is still not converged. 
        WRITE(0,*) 'SolveMumps.f90 : ',error
        ERROR STOP 122
    ENDIF    
ENDDO
IF(check) THEN
    error="The solution converges to a local minimum. Please restart &
	&  with different number of load steps and number of iterations"
    WRITE(0,*) 'SolveMumps.f90 : ',error
    WRITE(0,*) MAXVAL(ABS(gradient)*MAX(ABS(x),1.0D0)/MAX(fmin,0.5D0*SIZE(x)))
    ERROR STOP 123
ENDIF
        
IF ( mumps_par%MYID .eq. 0 )THEN
	DEALLOCATE( mumps_par%RHS )
	mumps_par%JOB = -2
ENDIF
! Destroy the instance (deallocate internal data structures)

CALL DMUMPS(mumps_par)
CALL MPI_FINALIZE(IERR)

CONTAINS ! The line search subroutine 

!*************************************************************
!*                                                           *                                      
!>  Use line search to improve the convergence of Newton
!! Raphson method, modified from the book: Numerical Recipes*
!*															 *
!************************************************************* 
SUBROUTINE LineSearch

REAL(DBL),PARAMETER:: ALF=1.0D-4 !< A small number to indicate sufficient decrease of the function
REAL(DBL),PARAMETER:: TOLY=1.0D-12 !< a small number to calculate the minimum step size
REAL(DBL):: tmp, slope,alamin, alam,tmplam,rhs1,rhs2,f2,alam2,a,b,disc

check=.FALSE.

tmp=Norm(dx)

IF(tmp>stpmax) dx=dx*stpmax/tmp  ! Scale if attempted step is too big

slope=DOT_PRODUCT(gradient,dx)

!~ IF(slope>=0.0D0) THEN
!~     error='roundoff problem in line search'
!~     WRITE(0,*) 'SolveMumps.f90 : ',error
!~     ERROR STOP 124
!~ ENDIF

alamin=TOLY/MAXVAL(ABS(dx)/MAX(ABS(xold),1.0D0))  ! Compute lambda_min

alam=1.0D0

DO 

  x=xold+alam*dx
  
  CALL AssembleRHS(ndof_el,memb_info,v_root_a,omega_a,member,error,&
	& ncond_mb,mb_condition,distr_fun,(dof_con),x,rhs,aero_flag,grav_flag,init_cond)
  IF(error/='') RETURN

  fmin=0.5D0*DOT_PRODUCT(rhs,rhs)
 
  IF(alam<alamin) THEN  ! CONVERGENCE ON deltaX, the calling program should verify the convergence
    x=xold
	check=.TRUE.
	RETURN
  ELSE IF(fmin<=fold+ALF*alam*slope) THEN
	   RETURN  ! SUFFICIENT function decrease
  ELSE 
     IF(ABS(alam-1.0D0)<TOLERANCE) THEN
	   tmplam=-slope/(2.0D0*(fmin-fold-slope))
	 ELSE  
	   rhs1=fmin-fold-alam*slope
	   rhs2=f2-fold- alam2*slope
	   a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
       b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
       
	   IF(ABS(a)<TOLERANCE)THEN
           tmplam=-slope/(2.0D0*b)
	   ELSE 
	       disc=b*b-3.0D0*a*slope
		   IF(disc<0.0d0) THEN
		       tmplam=0.5D0*alam  ! for imagine roots choose the max allowed
           ELSE 
		       tmplam=(-b+sqrt(disc))/(3.0D0*a)
		   ENDIF
       ENDIF
       IF(tmplam>0.5d0*alam) tmplam=0.5d0*alam 
     ENDIF
  ENDIF
  
  alam2=alam
  f2=fmin
  alam=MAX(tmplam,0.1D0*alam)

ENDDO
END SUBROUTINE LineSearch
!************************************************************ 


END SUBROUTINE NewtonRaphsonMumps
!***************************************************************
!***************************************************************
!*                                                             *
!> The subroutine extracts the solution for each key point 
!! and each member from the solution vector                    *
!*                                                             *
!***************************************************************
SUBROUTINE ExtractSolution(ndof_el,member,coord,memb_info,x,dof_con,sol_pt_i,sol_mb_i)

INTEGER,INTENT(IN)::ndof_el !<#ioaero::ndof_el
INTEGER,INTENT(IN)::member(:,:) !<#ioaero::member
REAL(DBL),INTENT(IN)::coord(:,:)    !<#ioaero::coord
TYPE (MemberInf),INTENT(IN)::memb_info(:)   !< contains the member parameters of the whole structure
REAL(DBL),INTENT(IN)::x(:) !< the solution vector
INTEGER,INTENT(IN)::dof_con(:) !< the connecting condition for key point.

REAL(DBL),INTENT(OUT)      ::sol_pt_i(:,:)     !< solutions for points for ith step
REAL(DBL),INTENT(OUT)      ::sol_mb_i(:,:)     !< solutions for members for ith step


INTEGER::i,j,ncol,pre_dof,free_dof

INTEGER::ndiv !< divisions of the member
INTEGER::val_no !< the position for a value in the member solution vector


! Extract values for the point
!---------------------------------------------	
DO i=1,SIZE(dof_con)

	sol_pt_i(i,1:3)=coord(i,:) ! Extract the coordinates

    ncol=index_kp(i,2)
	
	! for a follower, the prescribed values should be updated
	! by the current orientation calculated from current solution
	!---------------------------------------------------------------------
	IF(ANY(follower_all(i,:)==1)) &  
	& cond_all(i,:)=UpdateFollower(dof_all(i,:),follower_all(i,:),cond_all(i,:),x(ncol:ncol+5))

	IF(dof_con(i)==0) THEN	   ! a boundary point
		DO j=1,NSTRN	    
			pre_dof=dof_all(i,j)
			IF(pre_dof/=0) THEN
				IF(pre_dof>NSTRN) THEN
					free_dof=pre_dof-NSTRN
				ELSE 
					free_dof=pre_dof+NSTRN  
				ENDIF
				sol_pt_i(i,3+pre_dof)=cond_all(i,j)
				sol_pt_i(i,3+free_dof)=x(ncol+j-1)
			ENDIF
		ENDDO   
	ELSE   ! a connection point
		sol_pt_i(i,4:9)=x(ncol:ncol+5)  ! note that forces/moments of connection points are not variables 
	ENDIF
	
ENDDO

val_no=0

! Extract values for the member
!---------------------------------------------	
DO i=1,SIZE(member,1)

	ndiv=memb_info(i)%ndiv

    ncol=index_mb(i,2)
	DO j=1,ndiv
	   val_no=val_no+1
	   sol_mb_i(val_no,1:3)=memb_info(i)%coordinate(j,:)
	   sol_mb_i(val_no,4:3+ndof_el)=x(ncol:ncol+ndof_el-1)
       ncol=ncol+ndof_el
	ENDDO

ENDDO

END SUBROUTINE ExtractSolution
!************************************************************


!***************************************************************
!*                                                             *
!> The subroutine extracts elemental values from               *
!! the solution vector						                   *
!*                                                             *
!***************************************************************
SUBROUTINE ExtractElementValues(ndof_el,member,x,sol_mb_i)

INTEGER,INTENT(IN)::ndof_el !<#ioaero::ndof_el
INTEGER,INTENT(IN)::member(:,:) !<#ioaero::member
REAL(DBL),INTENT(IN)::x(:) !< the solution vector
REAL(DBL),INTENT(OUT)      ::sol_mb_i(:,:)     !< solutions for all the elements for ith step


INTEGER::i,j,ncol

INTEGER::val_no !< the position for a value in the member solution vector

val_no=0

! Extract values for the member
!---------------------------------------------	
DO i=1,SIZE(member,1)

    ncol=index_mb(i,2)
	DO j=1,member(i,6)
	   val_no=val_no+1
	   sol_mb_i(val_no,4:3+ndof_el)=x(ncol:ncol+ndof_el-1)
       ncol=ncol+ndof_el
	ENDDO
ENDDO

END SUBROUTINE ExtractElementValues
!************************************************************


!***************************************************************
!*                                                             *
!> The subroutine insert the elemental values into the solution
!! vector: needed for initial guess for starting
!! time marching: replace the first six valumes for 
!! each element with given initial conditions *
!***************************************************************
SUBROUTINE InsertElementValues(ndof_el,member,x,init_cond)

INTEGER,INTENT(IN)::ndof_el !<#ioaero::ndof_el
INTEGER,INTENT(IN)::member(:,:) !<#ioaero::memeber
REAL(DBL),INTENT(INOUT)::x(:) !< the solution vector

REAL(DBL),INTENT(IN)::init_cond(:,:)    !<#ioaero::init_cond

INTEGER::i,j,ncol

INTEGER::val_no !< the position for a value in the member solution vector


ncol=0
val_no=0

! Insert values into the solution vector
!---------------------------------------------	
DO i=1,SIZE(member,1)

    ncol=index_mb(i,2)
	DO j=1,member(i,6)
	   val_no=val_no+1
	   x(ncol:ncol+5)=init_cond(val_no,1:6)
       ncol=ncol+ndof_el
	ENDDO
ENDDO


END SUBROUTINE InsertElementValues
!************************************************************


!> Transfer the vector PH (linear and angular momenta) form frame B to frame a
FUNCTION CTCabPH(niter,member,memb_info ,sol_mb_i)

INTEGER,INTENT(IN)::niter   !<#ioaero::niter
INTEGER,INTENT(IN)::member(:,:)     !<#ioaero::member
TYPE (MemberInf),INTENT(IN)::memb_info(:)   
REAL(DBL),INTENT(IN) ::sol_mb_i(:,:)     !< solutions for members for ith step
REAL(DBL)            ::CTCabPH(SIZE(sol_mb_i,1),6)

INTEGER::i,j,val_no
REAL(DBL)::eCab(NDIM,NDIM),theta(NDIM),eCTCab(NDIM,NDIM)

val_no=0

DO i=1,SIZE(member,1)

	DO j=1,member(i,6)
		eCab=memb_info(i)%triad(j,:,:)

	   val_no=val_no+1
	   IF(niter/=1) THEN
 			theta  =  sol_mb_i(val_no,4:6)          ! given initial rotations
			eCTCab =MATMUL(DirCosineTRodrigues(theta),eCab)
	   ELSE
		    eCTCab=eCab
	   ENDIF
 	   CTCabPH(val_no,1:3)=MATMUL(eCTCab,sol_mb_i(val_no,13:15))  ! CTCabP
	   CTCabPH(val_no,4:6)=MATMUL(eCTCab,sol_mb_i(val_no,16:18)) ! CTCabH
	ENDDO

ENDDO
END FUNCTION CTCabPH
!************************************************************
END MODULE SolveMumps
!=========================================================
