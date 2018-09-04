!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************
!=========================================================
!
!  This module contains the main routines needed for eigen value
! analysis  
!
!=========================================================
MODULE EigenSlepc

USE InternalData
USE PrescribedCondition
USE System
USE GenericSolver


IMPLICIT NONE

PRIVATE       ! So everything is private, except declared by PUBLIC
PUBLIC EigenSolve

CONTAINS
!===========================================================


!**************************************************************
!*                                                            *                                      
!*  Assemble the linear matrices                              *
!*															  *
!**************************************************************
SUBROUTINE EigenSolve(ndof_el,memb_info,v_root_a,omega_a,member,pt_condition,niter,error,&
           & ncond_mb,mb_condition,distr_fun,dof_con,x,nev,eigen_val,eigen_vec,aero_flag,grav_flag)


INTEGER,INTENT(IN)::ndof_el,aero_flag,grav_flag
TYPE (MemberInf),INTENT(IN)::memb_info(:)

REAL(DBL),INTENT(IN)::v_root_a(:),omega_a(:)

REAL(DBL),INTENT(IN) ::distr_fun(:,:)
INTEGER,INTENT(IN)   ::member(:,:),niter,ncond_mb
TYPE(PrescriInf),INTENT(IN)::pt_condition(:), mb_condition(:) 
INTEGER                 ::dof_con(:) ! this array is passed by value
CHARACTER(*),INTENT(OUT)::error

REAL(DBL),INTENT(IN) ::x(:)

INTEGER,INTENT(INOUT)::nev
REAL(DBL),INTENT(OUT)::eigen_val(:,:) !eigen_val(1,:): real parts; and eigen_val(2,:): imaginary parts of eigenvalue
REAL(DBL),INTENT(OUT)::eigen_vec(:,:) ! the first nev (or nev+1) columns are the eigenvectors

! variables needed for MC19
!----------------------------------------------------------
REAL(DBL)            ::r_mc19(nsize),c_mc19(nsize)
!----------------------------------------------------------

! variables needed for MA28
!----------------------------------------------------------
REAL(DBL),ALLOCATABLE::coef_ma28(:)
INTEGER,ALLOCATABLE::irn_ma28(:),jcn_ma28(:),keep_ma28(:)
INTEGER::ljcn,lirn,neA,neB

! variables needed for ARPACK
!----------------------------------------------------------
REAL(DBL)::er(nev+1) ,ei(nev+1) 
INTEGER:: ncv=0

!----------------------------------------------------------
INTEGER::j
REAL(DBL)::tmpR=0.D0

!~ #if defined(PETSC_USE_FORTRAN_DATATYPES)
!~       type(Mat)      A
!~       type(EPS)      eps
!~ #else
!~       Mat            A
!~       EPS            eps
!~ #endif

! Assemble the coefficient matrix
!---------------------------------------------------------
CALL AssembleJacobian(ndof_el,niter,memb_info,v_root_a,omega_a,member,error,&
	& ncond_mb,mb_condition,distr_fun,(dof_con),x,aero_flag,grav_flag)
IF(error/='') RETURN
neA = ne
ljcn=5*ne
lirn=2*ne
ALLOCATE(coef_ma28(ljcn),STAT=allo_stat)
IF(MemoryError('coef_ma28',error)) GOTO 9999
coef_ma28(1:ne)=coef(1:ne)
DEALLOCATE(coef)

ALLOCATE(irn_ma28(lirn),STAT=allo_stat)
IF(MemoryError('irn_ma28',error)) GOTO 9999
irn_ma28(1:ne)=irn(1:ne)
DEALLOCATE(irn)

ALLOCATE(jcn_ma28(ljcn),STAT=allo_stat)
IF(MemoryError('jcn_ma28',error)) GOTO 9999
jcn_ma28(1:ne)=jcn(1:ne)
DEALLOCATE(jcn)

ALLOCATE(keep_ma28(5*nsize),STAT=allo_stat)
IF(MemoryError('ikeep',error)) GOTO 9999


! Process coef before entering eigen solver
!--------------------------------------------------
!~ CALL ProcessKMatrix(ljcn,lirn,coef_ma28,irn_ma28,jcn_ma28,r_mc19,c_mc19,keep_ma28,error)
!~ IF(error/='') GOTO 9999
!~ DEALLOCATE(irn_ma28)

! Assemble the mass matrix
!---------------------------------------------------------
assemble_flag=1
CALL AssembleJacobian(ndof_el,niter,memb_info,v_root_a,omega_a,member,error,&
	& ncond_mb,mb_condition,distr_fun,(dof_con),x,aero_flag,grav_flag)
IF(error/='') GOTO 9999
! Note ne has been changed to be the nonzero values of the mass matrix

!Using ARPACK to find the needed eigen values and eigen vectors
!-----------------------------------------------------------------------------
neB = ne
ncv=MIN(nsize,2*nev)

!~ CALL ARPACK(r_mc19,c_mc19,coef_ma28,ljcn,jcn_ma28,keep_ma28,nev,ncv,er,ei,eigen_vec,error)
!~ IF(error/='')GOTO 9999
CALL GepSolver(irn_ma28,jcn_ma28,coef_ma28,neA,irn,jcn,coef,neB,ncv)

!DO j=1,nev
!      WRITE(*,*)"Eigenvalues"
!      WRITE (*,'(1x,2G12.4)') er(j),ei(j)
!ENDDO

! The real eigen value is -1/lambda
!--------------------------------------
DO j=1,nev
   IF(ei(j)/=0.0d0) THEN
		tmpR=1.0d0/(er(j)*er(j)+ei(j)*ei(j))
		er(j)=-er(j)*tmpR
		ei(j)= ei(j)*tmpR
	ELSE
	   er(j)=-1.0d0/er(j)
	ENDIF
ENDDO

tmpR=1/(2.0d0*PI)

eigen_val(1,1:nev)=er(1:nev)*tmpR
eigen_val(2,1:nev)=ei(1:nev)*tmpR


DEAllOCATE(keep_ma28,jcn_ma28,coef_ma28)
IF(ALLOCATED(coef)) 	 DEALLOCATE(coef)
IF(ALLOCATED(irn)) 	 DEALLOCATE(irn)
IF(ALLOCATED(jcn)) 	 DEALLOCATE(jcn)

9999 IF(error/='')THEN
		IF(ALLOCATED(keep_ma28)) 	 DEALLOCATE(keep_ma28)
		IF(ALLOCATED(jcn_ma28)) 	 DEALLOCATE(jcn_ma28)
		IF(ALLOCATED(irn_ma28)) 	 DEALLOCATE(irn_ma28)
		IF(ALLOCATED(coef_ma28)) 	 DEALLOCATE(coef_ma28)
	 ENDIF

END SUBROUTINE EigenSolve
!***************************************************************

 
!----------------------------------------------------

!*************************************************************
! Process the K matrix before entering the eigen solver
! including:
!         Scaling
!         Setting default controls
!         Finding pivoting order
!         Factorization
!*************************************************************
SUBROUTINE ProcessKMatrix(ljcn,lirn,coef_ma28,irn_ma28,jcn_ma28,r_mc19,c_mc19,keep_ma28,error)

IMPLICIT NONE
INTERFACE
   SUBROUTINE MC19AD(nsize,ne,coef_ma28,irn_ma28,jcn_ma28,r_mc19,c_mc19,w_mc19)
		USE	GlobalDataFun
		INTEGER::nsize,ne,irn_ma28(ne),jcn_ma28(ne)
		REAL(DBL)::coef_ma28(ne)
		REAL(DBL)::r_mc19(nsize), c_mc19(nsize),w_mc19(5*nsize)
	END SUBROUTINE MC19AD
	 SUBROUTINE MA28AD(nsize,ne,coef_ma28,ljcn,irn_ma28,lirn,jcn_ma28,u,ikeep,iw,r_mc19,iflag)
		USE	GlobalDataFun
		INTEGER::nsize,ne,ljcn,irn_ma28(lirn),lirn,jcn_ma28(ljcn),ikeep(5*nsize),iw(8*nsize),iflag
		REAL(DBL)::coef_ma28(ljcn),u,r_mc19(nsize)
	END SUBROUTINE MA28AD

END INTERFACE

INTEGER,INTENT(IN)::ljcn,lirn
REAL(DBL),INTENT(INOUT) ::coef_ma28(:)
INTEGER,INTENT(INOUT)   ::irn_ma28(:),jcn_ma28(:)
REAL(DBL),INTENT(OUT)   ::r_mc19(:),c_mc19(:)
INTEGER,INTENT(OUT)     ::keep_ma28(:)
CHARACTER(*),INTENT(OUT)::error         ! a character variable holding  error message

REAL(DBL)::w_mc19(5*nsize)  ! work space for MC19AD
INTEGER::iw(8*nsize),iflag  ! work space for	MA28AD
INTEGER::i
REAL(DBL)::u

!Scale input matrices using MC19, coef_ma28 is used as a temp array
!-------------------------------
!~ CALL MC19AD(nsize,ne,coef_ma28(1:ne),irn_ma28(1:ne),jcn_ma28(1:ne),r_mc19,c_mc19,w_mc19)

!~ r_mc19=EXP(r_mc19)
!~ c_mc19=EXP(c_mc19)

!~ DO i=1,ne
!~    coef_ma28(i)= coef_ma28(i)*r_mc19(irn_ma28(i))*c_mc19(jcn_ma28(i))
!~ ENDDO

!------------------------------------------------------------------------------
! DECOMPOSE MATRIX INTO ITS FACTORS, w_MC19 is used as a working array
!------------------------------------------------------------------------------
u=0.1D0
CALL MA28AD(nsize,ne,coef_ma28,ljcn,irn_ma28,lirn,jcn_ma28,u,keep_ma28,iw,w_mc19,iflag)
IF (iflag.lt.0)THEN 
   WRITE(error,'(A,I2)') ' Failure of MA28AD with iflag=', iflag
   RETURN
ENDIF

END SUBROUTINE ProcessKMatrix
!**********************************************************************************


!*************************************************************
! Use ARPACK to solve the eigen value problem, modified from 
! dnsimp: for double precision nonsymmetric problem 
!
!*************************************************************
SUBROUTINE ARPACK(r_mc19,c_mc19,coef_ma28,ljcn,jcn_ma28,keep_ma28,nev,ncv,er,ei,vector,error)

IMPLICIT NONE
REAL(DBL),INTENT(IN)    ::r_mc19(:),c_mc19(:),coef_ma28(:)
INTEGER,INTENT(IN)      ::ljcn,jcn_ma28(:),keep_ma28(:)


INTEGER,INTENT(INOUT)::nev
INTEGER,INTENT(IN)::ncv
CHARACTER(*),INTENT(OUT)::error         ! a character variable holding  error message
REAL(DBL),INTENT(OUT)::er(:),ei(:),vector(:,:)
INTEGER::ido,iparam(11),ipntr(14),lworkl,info
REAL(DBL)::resid(nsize),v(nsize,ncv),workd(3*nsize),workl(3*ncv**2 + 6*ncv)

LOGICAL::select(ncv)
REAL(DBL)::sigmar,sigmai,workev(3*ncv)
INTEGER::ierr

ido=0


iparam(1)=1  ! exact shifts
iparam(3) = 10 ! matrix iterations
iparam(7) = 1  ! mode 1 is used solve A*x = lambda*x in a regular mode

lworkl=3*ncv**2 + 6*ncv

info=0



DO 

	! %---------------------------------------------%
	! | Repeatedly call the routine DNAUPD and take |
	! | actions indicated by parameter IDO until    |
	! | either convergence is indicated or maxitr   |
	! | has been exceeded.                          |
	! %---------------------------------------------%
	CALL dnaupd ( ido, 'I', nsize, 'LM', nev, TOLERANCE, resid, ncv, &
			 &     v, nsize, iparam, ipntr, workd, workl, lworkl,& 
			 &     info )
	IF( info < 0 ) THEN
		WRITE(error,*)' Error with dnaupd, info = ',info,' Check the documentation of dnaupd.'
		RETURN
	ELSEIF (ido == -1 .OR. ido == 1) THEN
	! %-------------------------------------------%
	! | Perform matrix vector multiplication      |
	! |                y <--- Op*x                |
	! | The user should supply his/her own        |
	! | matrix vector multiplication routine here |
	! | that takes workd(ipntr(1)) as the input   |
	! | vector, and return the matrix vector      |
	! | product to workd(ipntr(2)).               | 
	! %-------------------------------------------%

		CALL AW(workd(ipntr(2)),ne,irn,jcn,coef,workd(ipntr(1)), &
			  & r_mc19,c_mc19,coef_ma28,ljcn,jcn_ma28,keep_ma28)
	    
	! %-----------------------------------------%
	! | L O O P   B A C K to call DNAUPD again. |
	! %-----------------------------------------%
		ELSE
			EXIT
	ENDIF 
ENDDO

CALL dneupd (.TRUE., 'A', select, er, ei, v, nsize, sigmar, sigmai, workev, &
           & 'I', nsize, 'LM', nev, TOLERANCE, resid, ncv, v, &
           & nsize, iparam, ipntr, workd, workl, lworkl, ierr )
IF( ierr < 0 ) THEN
	WRITE(error,*)' Error with dneupd, info = ',ierr,' Check the documentation of dneupd.'
	RETURN
ELSE
	nev = iparam(5)
	vector(:,1:nev)=v(:,1:nev)
ENDIF

9999 RETURN
END SUBROUTINE  ARPACK
!******************************************************



!******************************************************
!  Matrix multipy A*W, A is sparse
!			rhs=coef_mass*W
!           coef_ma48*U=rhs
! note, coef_mass and coef_ma48 remains the same during eigen solution
!********************************************************
SUBROUTINE AW(u,nzM,irnM,jcnM,coef_mass,w,r_mc19,c_mc19,coef_ma28,ljcn,jcn_ma28,keep_ma28)
IMPLICIT NONE
INTERFACE
	SUBROUTINE MA28CD(nsize,coef_ma28,ljcn,jcn_ma28,ikeep,rhs,r_mc19,mytype)
		USE	GlobalDataFun
		INTEGER::nsize,ljcn,jcn_ma28(ljcn),ikeep(5*nsize),mytype
		REAL(DBL)::coef_ma28(ljcn),r_mc19(nsize),rhs(nsize)
	END SUBROUTINE MA28CD
END INTERFACE

INTEGER,INTENT(IN)::nzM,irnM(:),jcnM(:),ljcn,jcn_ma28(:),keep_ma28(:)
REAL(DBL),INTENT(IN) ::coef_mass(:),w(nsize),r_mc19(:),c_mc19(:),coef_ma28(:)
REAL(DBL),INTENT(OUT) ::u(nsize)

INTEGER i,j,k

!-------------------------------------
! Compute u=coef_mass*W
!-------------------------------------	   
u=0.0d0
DO k = 1,nzM
   i = irnM(k)
   j = jcnM(k)
   u(i) = u(i) + coef_mass(k)*w(j)
ENDDO

!------------------------------------------------
! scale U correspondingly and store it in rhs
!------------------------------------------------
DO i=1,nsize
	u(i)=u(i)*r_mc19(i)
ENDDO

!------------------------------
! Solve linear system.
!------------------------------
CALL MA28CD(nsize,coef_ma28,ljcn,jcn_ma28,keep_ma28,u,w,1)

DO i=1,nsize
	u(i)=u(i)*c_mc19(i)
ENDDO

     

9999 RETURN      
END SUBROUTINE AW
!******************************************************


END MODULE EigenSlepc
!=========================================================

