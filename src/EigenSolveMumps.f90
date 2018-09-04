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
MODULE EigenMumps

USE InternalData
USE PrescribedCondition
USE System


IMPLICIT NONE

INCLUDE 'mpif.h'
INCLUDE 'dmumps_struc.h'
TYPE (DMUMPS_STRUC) mumps_par
INTEGER IERR


PRIVATE       ! So everything is private, except declared by PUBLIC
PUBLIC EigenSolveMumps

CONTAINS
!===========================================================


!**************************************************************
!*                                                            *                                      
!*  Assemble the linear matrices                              *
!*															  *
!**************************************************************
SUBROUTINE EigenSolveMumps(ndof_el,memb_info,v_root_a,omega_a,member,pt_condition,niter,error,&
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

! Lapack DGEEV parameters
INTEGER,PARAMETER :: LWMAX = 1000
REAL(DBL) :: A(NSTATES,NSTATES),EigAR(NSTATES),EigAI(NSTATES),VL(NSTATES),VR(NSTATES),WORK(LWMAX),U,Chord
INTEGER :: INFO,LWORK
REAL(DBL) :: IFmodes(SIZE(member,1),NSTATES,2) ! Peters induced flow modes ; not to output

! variables needed for ARPACK
!----------------------------------------------------------
REAL(DBL)::er(nev+1) ,ei(nev+1) 
INTEGER:: ncv=0, imodes=1

!----------------------------------------------------------
INTEGER::i,j,k
REAL(DBL)::tmpR=0.D0

CALL MPI_INIT(IERR)
! Define a communicator for the package.

mumps_par%COMM = MPI_COMM_WORLD
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
CALL DMUMPS(mumps_par)
IF (mumps_par%INFOG(1) <0 ) THEN
    WRITE (0,*) 'EigensolveMumps.f90, error in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
    ERROR STOP 200
ENDIF    
!  Control parameters (see Mumps user gudie))
mumps_par%ICNTL(1) = 6
mumps_par%ICNTL(2) = 0
mumps_par%ICNTL(3) = 0
mumps_par%ICNTL(6) = 5
mumps_par%ICNTL(7) = 3
mumps_par%ICNTL(8) = -2
mumps_par%ICNTL(28) = 1
mumps_par%ICNTL(35) = 1

! Analysis and factorisation 
mumps_par%JOB = 4

! Assemble the coefficient matrix
!---------------------------------------------------------
CALL AssembleJacobian(ndof_el,niter,memb_info,v_root_a,omega_a,member,error,&
	& ncond_mb,mb_condition,distr_fun,(dof_con),x,aero_flag,grav_flag)


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

	IF(ALLOCATED(coef)) 	 DEALLOCATE(coef)
	IF(ALLOCATED(irn)) 	 DEALLOCATE(irn)
	IF(ALLOCATED(jcn)) 	 DEALLOCATE(jcn)

	! Assemble the mass matrix
	!---------------------------------------------------------
	assemble_flag=1
	CALL AssembleJacobian(ndof_el,niter,memb_info,v_root_a,omega_a,member,error,&
		& ncond_mb,mb_condition,distr_fun,(dof_con),x,aero_flag,grav_flag)
    assemble_flag=0    
	! Note ne has been changed to be the nonzero values of the mass matrix
ENDIF
	! preprocessing of Mumps solver
	CALL DMUMPS(mumps_par)
    IF (mumps_par%INFOG(1) <0 ) THEN
    WRITE (0,*) 'EigensolveMumps.f90, error in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
        ERROR STOP 200
    ENDIF    
    mumps_par%JOB = 3
	!Using ARPACK to find the needed eigen values and eigen vectors
	!-----------------------------------------------------------------------------
	ncv=MAX(14,MIN(nsize,2*nev))

	CALL ARPACK(nev,ncv,er,ei,eigen_vec,error)

	IF(ALLOCATED(coef)) 	 DEALLOCATE(coef)
	IF(ALLOCATED(irn)) 	 DEALLOCATE(irn)
	IF(ALLOCATED(jcn)) 	 DEALLOCATE(jcn)
    
!~      ! Computation of Peters MatriX A eigenvalue in order eliminate modes related to induced flow states
!~     ! Using Lapack routine DGEEV to compute eigenvalues of Matrix A
!~     IF (aero_flag ==3) THEN
!~         A = MatA(NSTATES)
!~         LWORK = -1
!~         CALL DGEEV ('N','N',NSTATES,A,NSTATES,EigAR,EigAI,VL,NSTATES,VR,NSTATES,WORK,LWORK,INFO)
!~         LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!~         CALL DGEEV ('N','N',NSTATES,A,NSTATES,EigAR,EigAI,VL,NSTATES,VR,NSTATES,WORK,LWORK,INFO)
!~         DO i=1,NSTATES
!~                 DO j=1,NSTATES
!~                 WRITE(*,*) A(i,j)-A(j,i)
!~             ENDDO    
!~         ENDDO    
!~         DO i=1,SIZE(member,1)
!~             U = memb_info(1)%aerodyn_coef(i,1)
!~             U = max(1.D0/max(1e-12,U),U)
!~             Chord = memb_info(1)%aerodyn_coef(i,3)
!~             IFmodes(i,:,1) = PI*Chord/U*EigAR(:)
!~             IFmodes(i,:,2) = PI*Chord/U*EigAI(:)
!~             DO j=1,NSTATES
!~                 tmpR = 1.D0/(IFmodes(i,j,1)**2+IFmodes(i,j,2)**2)
!~                 IFmodes(i,j,1) = - IFmodes(i,j,1)*tmpR
!~                 IFmodes(i,j,2) = - IFmodes(i,j,2)*tmpR
!~             ENDDO
!~         ENDDO    
!~     ENDIF    
    
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
    
!~     IF (aero_flag ==3) THEN
!~         DO j=1,nev
        
!~             ! temporaire : calcul de la norme des valeurs propres avant inversion
!~             WRITE(*,*) er(j),ei(j)

!~             DO i=1,NSTATES
!~                 DO k=1,SIZE(member,1)
!~                     IF(ABS(eigen_val(1,j)-IFmodes(k,i,1))/ABS(eigen_val(1,j))+&
!~                     &ABS(eigen_val(2,j)-IFmodes(k,i,2))/min(1e-12,ABS(eigen_val(2,j)))<1e-3) THEN
!~                     IF(ABS(eigen_val(1,j)-IFmodes(k,i,1))/ABS(eigen_val(1,j))<1e-3) THEN
!~                         IF (ABS(IFmodes(k,i,2))==0.D0) THEN
!~                             eigen_val(:,j) = 0.D0
!~                         ENDIF    
!~                     ENDIF
!~                 ENDDO    
!~             ENDDO    
!~         ENDDO     
!~     ENDIF     
    

	!  Deallocate user data
IF ( mumps_par%MYID .eq. 0 )THEN
	DEALLOCATE( mumps_par%IRN )
	DEALLOCATE( mumps_par%JCN )
	DEALLOCATE( mumps_par%A )
	DEALLOCATE( mumps_par%RHS )
END IF

! Mumps error checking

! Destroy the instance (deallocate internal data structures)
mumps_par%JOB = -2
CALL DMUMPS(mumps_par)
IF (mumps_par%INFOG(1) <0 ) THEN
    WRITE (0,*) 'EigensolveMumps.f90, error in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
    ERROR STOP 200
ENDIF    
CALL MPI_FINALIZE(IERR)

END SUBROUTINE EigenSolveMUMPS

!*************************************************************
SUBROUTINE ARPACK(nev,ncv,er,ei,vector,error)

IMPLICIT NONE

INTEGER,INTENT(INOUT)::nev
INTEGER,INTENT(IN)::ncv
CHARACTER(*),INTENT(OUT)::error         ! a character variable holding  error message
REAL(DBL),INTENT(OUT)::er(:),ei(:),vector(:,:)
INTEGER::ido,iparam(11),ipntr(14),lworkl,info
REAL(DBL)::resid(nsize),v(nsize,ncv),workd(3*nsize),workl(3*ncv**2 + 6*ncv)
CHARACTER(2) :: WichA,WichE

LOGICAL::select(ncv)
REAL(DBL)::sigmar,sigmai,workev(3*ncv)
INTEGER::ierr

ido=0
iparam(1)=1  ! exact shifts
iparam(3) = 10 ! matrix iterations
iparam(7) = 1  ! mode 1 is used solve A*x = lambda*x in a regular mode

lworkl=3*ncv**2 + 6*ncv
info=0


! define the spectral 
SELECT CASE(ARPACK_MOD)
    CASE(1) ! only modes with non-zero frequency
        WichA = 'LI'
        WichE = 'LI'
    CASE(2) ! modes with lowest module (including zero frequency modes)
        WichA = 'LM'
        WichE = 'LM'
    CASE(3) ! only zero frequency modes
        WichA = "LR"
        WichE = "LR"
    CASE(4)
        WichA = "SR"
        WichE = "SR"
    CASE DEFAULT ! mainly non-zero frequency modes including some low damping zero-frequency modes
        WichA = 'LM'
        WichE = 'LI'
END SELECT
DO 

	! %---------------------------------------------%
	! | Repeatedly call the routine DNAUPD and take |
	! | actions indicated by parameter IDO until    |
	! | either convergence is indicated or maxitr   |
	! | has been exceeded.                          |
	! %---------------------------------------------%
	CALL dnaupd ( ido, 'I', nsize, WichA, nev, TOLERANCE, resid, ncv, &
			 &     v, nsize, iparam, ipntr, workd, workl, lworkl,& 
			 &     info )

	IF( info < 0 ) THEN
		WRITE(error,*)' Error with dnaupd, info = ',info,' Check the documentation of dnaupd.'
		WRITE(0,*) 'EigenSolveMumps.f90 : ',error
        ERROR STOP 106
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

		CALL AW(workd(ipntr(2)),ne,irn,jcn,coef,workd(ipntr(1)))
	    
	! %-----------------------------------------%
	! | L O O P   B A C K to call DNAUPD again. |
	! %-----------------------------------------%
		ELSE
			EXIT
	ENDIF 
ENDDO

CALL dneupd (.TRUE., 'A', select, er, ei, v, nsize, sigmar, sigmai, workev, &
           & 'I', nsize, WichE, nev, TOLERANCE, resid, ncv, v, &
           & nsize, iparam, ipntr, workd, workl, lworkl, ierr )
           
!~ write(*,*) iparam(:)           

IF( ierr < 0 ) THEN
	WRITE(error,*)' Error with dneupd, info = ',ierr,' Check the documentation of dneupd.'
    WRITE(0,*) 'EigenSolveMumps.f90 : ',error
    ERROR STOP 106
ELSE
	nev = iparam(5)
	vector(:,1:nev)=v(:,1:nev)
ENDIF

END SUBROUTINE  ARPACK
!******************************************************



!******************************************************
!  Matrix multipy A*W, A is sparse
!			rhs=Mass_Matrix*W
!           Stiffness_Matrix*U=rhs
! note, Mass_Matrix and Stiffness_Matrix remains the same during eigen solution
!********************************************************
SUBROUTINE AW(u,nzM,irnM,jcnM,coef_mass,w)
IMPLICIT NONE

INTEGER,INTENT(IN)::nzM,irnM(:),jcnM(:)
REAL(DBL),INTENT(IN) ::coef_mass(:),w(nsize)
REAL(DBL),INTENT(OUT) ::u(nsize)

INTEGER i,j,k

!-------------------------------------
! Compute RHS=coef_mass*W
!-------------------------------------	   
IF ( mumps_par%MYID .eq. 0 ) THEN
	mumps_par%RHS=0
	DO k = 1,nzM
	   i = irnM(k)
	   j = jcnM(k)
	   mumps_par%RHS(i) = mumps_par%RHS(i) + coef_mass(k)*w(j)
	ENDDO
ENDIF

!  Call package for solution

CALL DMUMPS(mumps_par)
IF (mumps_par%INFOG(1) <0 ) THEN
    WRITE (0,*) 'EigensolveMumps.f90, error in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
    ERROR STOP 200
ENDIF    


!  Solution has been assembled on the host
      IF ( mumps_par%MYID .eq. 0 ) THEN
        u=mumps_par%RHS
      END IF
 
END SUBROUTINE AW

!~ !******************************************************
!~ ! COO format sparse matrix sum operation
!~ ! D = A-sigma*M
!~ ! A,D and M are double complex matrix	
!~ ! return matrix D in COO format
!~ ! *****************************************************
!~ SUBROUTINE CooSum(n,nza,nzm,nzd,sigma,irnA,jcnA,coefA,irnM,jcnM,coefM,irnD,jcnD,coefD)
!~ INTEGER, INTENT(IN) :: n,nza,nzm,irnA(:),jcnA(:),irnM(:),jcnM(:)
!~ INTEGER, INTENT(OUT) ::nzd
!~ INTEGER, INTENT(OUT),ALLOCATABLE :: irnD(:),jcnD(:)
!~ COMPLEX(DBL),INTENT(IN) :: coefA(:),coefM(:)
!~ COMPLEX(DBL),INTENT(IN) :: sigma
!~ COMPLEX(DBL),INTENT(OUT) ,ALLOCATABLE :: coefD(:)

!~ INTEGER nrow,jcol,indexA,indexM,indexD

!~ nzd = nza+nzm+1
!~ ALLOCATE(irnD(nzd))
!~ ALLOCATE(jcnD(nzd))
!~ ALLOCATE(coefD(nzd))

!~ ! initialisation with A matrix
!~ irnD(1:nza) = irnA(1:nza)
!~ jcnD(1:nza) = jcnA(1:nza)
!~ coefD(1:nza) = coefA(1:nza)
!~ indexD = nza+1
!~ ! searching if matrixA has a nonzero value at position (irnM,jcnM)
!~ DO indexM=1,nzm
!~ 	indexA = 1
!~ 	nrow = COUNT(irnA==irnM(indexM),1)
!~ 	IF (nrow>0) THEN
!~ 		indexA = MINLOC(irnA,1,irnA==irnM(indexM))
!~ 		jcol = MINLOC(jcnA(indexA:indexA+nrow-1),1,jcnA(indexA:indexA+nrow-1)==jcnM(indexM))
!~ 	ENDIF

!~ 	IF(jcol/=0) THEN
!~ 		indexA = indexA+jcol-1
!~ 	ELSE
!~ 		indexA = 0
!~ 	ENDIF
!~ 	! if A has a nonzero value, sum at indexA
!~ 	IF (indexA>0) THEN
!~ 		coefD(indexA) = coefD(indexA) - sigma*coefM(indexM)
!~ 	ELSE ! put the coefficient at the end of the matrix D
!~ 		irnD(indexD) = irnM(indexM)
!~ 		jcnD(indexD) = jcnM(indexM)
!~ 		coefD(indexD) = -sigma*coefM(indexM)
!~ 		indexD = indexD+1
!~ 	ENDIF
!~ ENDDO

!~ nzd = indexD-1
!~ irnD=irnD(1:nzd)
!~ jcnD=jcnD(1:nzd)
!~ coefD=coefD(1:nzd)

!~ END SUBROUTINE CooSum

!~ !******************************************************
!~ ! COO format sparse matrix index reordering 
!~ ! by row then by column
!~ ! preprocessing before coo sparse matrix sum subroutine
!~ ! *****************************************************

!~ SUBROUTINE CooSort(n,nz,irn,jcn,coef)
!~ INTEGER,INTENT(IN) :: n,nz
!~ INTEGER,INTENT(INOUT) :: irn(:),jcn(:)
!~ COMPLEX(DBL) coef(:)

!~ INTEGER tempi,tempj, irow,jcol,ncol,indexA,indexC,indexP,crit
!~ COMPLEX(DBL) tempz

!~ ! sorting by row
!~ DO indexA = 1,nz-1
!~ 	crit = MINLOC(irn(indexA:nz),1)
!~ 	IF (crit /= 1) THEN
!~ 		indexP = indexA + crit - 1
!~ 		! circular permutation
!~ 		tempi = irn(indexA)
!~ 		tempj = jcn(indexA)
!~ 		tempz = coef(indexA)
!~ 		irn(indexA) = irn(indexP)
!~ 		jcn(indexA) = jcn(indexP)
!~ 		coef(indexA) = coef(indexP)
!~ 		irn(indexP) = tempi
!~ 		jcn(indexP) = tempj
!~ 		coef(indexP) = tempz
!~ 	ENDIF
!~ ENDDO
!~ indexA = 1
!~ ! sorting by column
!~ DO WHILE (indexA <nz)
!~ 	irow = irn(indexA)
!~ 	ncol = COUNT(irn == irow,1)
!~ 	DO jcol = 1,ncol
!~ 		indexC = indexA+jcol-1
!~ 		crit = MINLOC(jcn(indexC:indexA+ncol-1),1)
!~ 		IF (crit /= 1) THEN
!~ 			indexP = indexC + crit - 1
!~ 			! circular permutation
!~ 			tempi = irn(indexC)
!~ 			tempj = jcn(indexC)
!~ 			tempz = coef(indexC)
!~ 			irn(indexC) = irn(indexP)
!~ 			jcn(indexC) = jcn(indexP)
!~ 			coef(indexC) = coef(indexP)
!~ 			irn(indexP) = tempi
!~ 			jcn(indexP) = tempj
!~ 			coef(indexP) = tempz
!~ 		ENDIF
!~ 	ENDDO
!~ 	indexA=indexA+ncol
!~ ENDDO
!~ END SUBROUTINE CooSort

END MODULE EigenMumps
!=========================================================

