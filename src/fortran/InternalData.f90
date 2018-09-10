!***********************************************************
!*                                                         *
!* This module contains the variables needed internally in *
!* the program. Not necessary to be defined in the outside *
!* environment. And these variables will not be passed     *
!* between the dll and outside environment.                *
!*                                                         *
!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************

!> This module contains the variables needed internally in the program. Not necessary to be defined in the outside environment

MODULE InternalData

USE GlobalDataFun


IMPLICIT NONE

LOGICAL,PARAMETER:: DEBUG=.FALSE. ! flag for debugging purpose

!File needed internally for debugging purpose
!===================================================
INTEGER,PARAMETER:: IOUT=30     ! file for debugging: inp_name.deb
CHARACTER(64)    :: deb_name 

!Global integer constants
!============================================================================
INTEGER::           nsize          ! the size of the problem
!=======================================================================

!Global variables
!============================================================================
INTEGER,ALLOCATABLE::dof_all(:,:)           ! prescribed condition for each point, system unique
INTEGER,ALLOCATABLE::follower_all(:,:)      ! follower condition for the prescribed condition for each point,system unique
REAL(DBL),ALLOCATABLE::cond_all(:,:)        ! the prescribed value for each point, load step unique for no follower, 
                                            ! if a follower, we need to update based on current deformation
REAL(DBL),ALLOCATABLE::init_memb(:,:)       ! initial displacement/rotations and time derivates for each member

INTEGER:: init_flag  ! 1-initial step; 2-time marching; 0, otherwise
REAL(DBL) ::two_divide_dt  !< 2/dt
INTEGER   ::assemble_flag=0  !< for the purpose to share the routines between assembly of stiffness matrix and mass matrix
INTEGER,ALLOCATABLE ::index_kp(:,:) !< the starting row and column for each kp
INTEGER,ALLOCATABLE::index_mb(:,:) !< the starting row and column for each member
REAL(DBL)::xyz_pt1(NDIM) !< the coordinate of the starting point of the first member

TYPE MemberInf !< structure containing the caracteritics of a finite element.
	 INTEGER   ::ndiv                       !< number of divisions
	 INTEGER   ::ncol_memb                  !< total number of columns of the member
	 REAL(DBL) ::dL                         !< length of each division 
	 REAL(DBL),ALLOCATABLE::mate(:,:,:)     !< mate(ndiv,12,6): flexibity and mass properties for each division
	 REAL(DBL),ALLOCATABLE::triad(:,:,:)    !< triad(ndiv,3,3): cab for each division, evaluated at the middle point of the division
     REAL(DBL),ALLOCATABLE::Le(:)           !< Le(ndiv): ending arc length for each division
	 REAL(DBL),ALLOCATABLE::coordinate(:,:) !< coordinate(ndiv,3) coordinate for the middle of the division
	 REAL(DBL),ALLOCATABLE::aerodyn_coef(:,:) !< aerodynamic coefficients
END TYPE MemberInf

INTEGER,PARAMETER::nzElemMax=500 ! maximum nonzero entries of the element matrix will be 318

INTEGER::neMax

END MODULE InternalData
!============================================================================
