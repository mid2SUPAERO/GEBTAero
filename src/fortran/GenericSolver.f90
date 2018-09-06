!***********************************************************
!*  Copyright (c) Bertrand Kirsch                                                *
!*  French Air Force Reseach Center                                          *
!***********************************************************

!=========================================================
!
!  This module contains generic solvers based on the Slepc and Petsc package
!	The following solver are implemented
!	A generelized eigenvalue problem solver
!
!=========================================================

MODULE GenericSolver


USE InternalData
USE PrescribedCondition
USE System

IMPLICIT NONE

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

PRIVATE
PUBLIC GepSolver

CONTAINS

SUBROUTINE GepSolver(irnA,jcnA,coefA,neA,irnB,jcnB,coefB,neB,ncv)

INTEGER :: neA,neB,ncv
INTEGER,ALLOCATABLE :: irnA(:),jcnA(:),irnB(:),jcnB(:)
REAL(DBL),ALLOCATABLE :: coefA(:),coefB(:)

#if defined(PETSC_USE_FORTRAN_DATATYPES)
	type(Mat)      A,B
	type(EPS)      eps
#else
	Mat            A,B
	EPS            eps
#endif
	EPSType        tname
	PetscInt       n, i, rstart, rend, cstart, cend
	PetscInt       row, col, maxit
	PetscMPIInt    rank
	PetscErrorCode ierr
	PetscBool      flg, terse
	PetscScalar    value,zero,tol

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	zero = 0.D0
	n = nsize
	call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
	call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
	call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-n',n,flg,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Compute the operator matrix that defines the eigensystem, Ax=kx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      call MatCreate(PETSC_COMM_WORLD,A,ierr)
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(A,ierr)
      call MatSetUp(A,ierr)
      
      call MatCreate(PETSC_COMM_WORLD,B,ierr)
      call MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(B,ierr)
      call MatSetUp(B,ierr)

      call MatGetOwnershipRange(A,rstart,rend,ierr)
      call MatGetOwnershipRangeColumn(A,cstart,cend,ierr)
      
	do i=1,nsize
		call MatSetValue(A,i-1,i-1,zero,INSERT_VALUES,ierr)
		call MatSetValue(B,i-1,i-1,zero,INSERT_VALUES,ierr)
	enddo
      do i=1,neA
			row = irnA(i)-1
			col = jcnA(i)-1
			value = coefA(i)
		    if(rstart < row+1)  then
				if (row < rend) then
					if(cstart < col+1) then
						if(col < cend) then
							call MatSetValue(A,row,col,value,INSERT_VALUES,ierr)
						endif
					endif
				endif
			endif
      enddo
		
	do i=1,neB
		row = irnB(i)-1
		col = jcnB(i)-1
		value = coefB(i)
		if(rstart < row+1)  then
			if (row < rend) then
				if(cstart < col+1) then
					if(col < cend) then
						call MatSetValue(B,row,col,value,INSERT_VALUES,ierr)
					endif
				endif
			endif
		endif
      enddo
      

      
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Create the eigensolver and display info
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!     ** Create eigensolver context
	call EPSCreate(PETSC_COMM_WORLD,eps,ierr)

!     ** Set operators. In this case, it is a standard eigenvalue problem
	call EPSSetOperators(eps,B,A,ierr)
!~ 	call EPSSetBalance(eps,EPS_BALANCE_TWOSIDE,10,100,ierr)
!~ 	call EPSSetProblemType(eps,EPS_GNHEP,ierr)

!     ** Set solver parameters at runtime

	call EPSSetFromOptions(eps,ierr)
!~ 	call EPSSetType(eps,EPSGD,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!~ 	call EPSSetTolerances(eps,1e-4,2000,ierr)
	call EPSSolve(eps,ierr) 
	
!     ** Optional: Get some information from the solver and display it
      call EPSGetType(eps,tname,ierr)
      if (rank .eq. 0) then
        write(*,120) tname
      endif
 120  format (' Solution method: ',A)
      call EPSSetDimensions(eps,ncv,min(nsize,2*ncv),nsize,ierr)
       call EPSGetTolerances(eps,tol,maxit,ierr)
      if (rank .eq. 0) then
		write(*,*) tol,maxit
        write(*,130) ncv
      endif
 130  format (' Number of requested eigenvalues:',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!     ** show detailed info unless -terse option is given by user
      call PetscOptionsHasName(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,  &
     &                        '-terse',terse,ierr)
      if (terse) then
        call EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_NULL_OBJECT,ierr)
      else
        call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,           &
     &                   PETSC_VIEWER_ASCII_INFO_DETAIL,ierr)
        call EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call EPSErrorView(eps,EPS_ERROR_RELATIVE,                       &
     &                   PETSC_VIEWER_STDOUT_WORLD,ierr)
        call PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr)
      endif
      call EPSDestroy(eps,ierr)
      call MatDestroy(A,ierr)
      call MatDestroy(B,ierr)

      call SlepcFinalize(ierr)

END SUBROUTINE GepSolver

END MODULE GenericSolver
