!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************

!***************************************************************
!*                                                             *
!> This module assemles within a member without considering the particular conditions of the end points 
!Outputs:MemberEqn, MemberJacobian                           *
!***************************************************************

MODULE Member
USE InternalData
USE PrescribedCondition
USE Element

IMPLICIT NONE

PRIVATE       ! So everything is private, except declared by PUBLIC
PUBLIC AssembleMemberRHS,AssembleMemberJacobian,ExtractMemberProperties,ndiv,ncol_memb 


INTEGER::ndiv      !< number of divisions
INTEGER::ncol_memb             !< the total number of columns of the member
!=============================================

CONTAINS

!************************************************************
!*                                                          *                                      
!>  Assemble the equations for each member                  *
!*															*
!************************************************************ 
SUBROUTINE  AssembleMemberRHS(ndof_el,memb_info_i,v_root_a,omega_a,x_memb,rhs_memb,aero_flag,grav_flag,init_cond)

INTEGER,INTENT(IN)::ndof_el !<#ioaero::ndof_el
INTEGER,INTENT(IN)::aero_flag   !<#ioaero::aero_flag
INTEGER,INTENT(IN)::grav_flag   !<#ioaero::grav_flag
TYPE (MemberInf),INTENT(IN)::memb_info_i    !< array containing the characteristics of a beam member

REAL(DBL),INTENT(IN)::v_root_a(:)   !< linear velocity of frame a
REAL(DBL),INTENT(IN)::omega_a(:)     !< angular velocity of frame a
REAL(DBL),INTENT(IN)::x_memb(:) !< solution vector of a beam member

REAL(DBL),INTENT(OUT)::rhs_memb(:)  !<RHS vector of a beam member

REAL(DBL),OPTIONAL,INTENT(IN)::init_cond(:,:)   !<#ioaero::init_cond

INTEGER:: i
INTEGER:: nrl !< the starting row/col for elemental rhs to fill in 
 
rhs_memb=0.0D0

nrl=0

DO i=1,ndiv
   CALL ExtractElementProperties(i,memb_info_i,x_memb(nrl+1:nrl+ndof_el),v_root_a,omega_a,ndof_el,init_memb(:,:),&
																							& aero_flag,grav_flag)
 
! assemble the equations within one element   
!------------------------------------------------------
   rhs_memb(nrl+1:nrl+NDOF_ND+ndof_el) =rhs_memb(nrl+1:nrl+NDOF_ND+ndof_el) +ElemEqn(ndof_el)

   nrl=nrl+ndof_el

ENDDO

END SUBROUTINE AssembleMemberRHS
!****************************************************************


!************************************************************
!*                                                          *                                      
!>  Assemble the Jacobian for each member                   *
!*															*
!************************************************************ 
SUBROUTINE  AssembleMemberJacobian(ndof_el,niter,memb_info_i,v_root_a,omega_a,x_memb,&
                                 & nz_memb,irn_memb,jcn_memb,coef_memb,aero_flag,grav_flag)

INTEGER,INTENT(IN)::ndof_el !<#ioaero::ndof_el
INTEGER,INTENT(IN)::niter   !<#ioaero::niter
INTEGER,INTENT(IN)::aero_flag   !<#ioaero::aero_flag
INTEGER,INTENT(IN)::grav_flag  !<#ioaero::grav_flag
TYPE (MemberInf),INTENT(IN)::memb_info_i    !< array containing the characteristics of a beam member

REAL(DBL),INTENT(IN)::v_root_a(:)    !< linear velocity of frame a
REAL(DBL),INTENT(IN)::omega_a(:)    !< angular velocity of frame a
REAL(DBL),INTENT(IN)::x_memb(:) !< solution vector of a beam member

INTEGER,INTENT(OUT)::nz_memb    !<Number of nonzero value of the member Jacobian
INTEGER,INTENT(OUT)::irn_memb(:)    !< line index of the nonzero coefficient
INTEGER,INTENT(OUT)::jcn_memb(:)    !< column index of the nonzero coefficient
REAL(DBL),INTENT(OUT)::coef_memb(:) !< value of the nonzero coefficient


INTEGER:: i,j
INTEGER:: nrl ! the starting row/col for elemental matrix to fill in 
 
REAL(DBL)::elemJac(NDOF_ND+ndof_el,ndof_el)

coef_memb=0.0D0

nrl=0
nz_memb=0.0D0

DO i=1,ndiv
 
   CALL ExtractElementProperties(i,memb_info_i,x_memb(nrl+1:nrl+ndof_el),v_root_a,omega_a,ndof_el,init_memb(:,:),&
								& aero_flag,grav_flag)

!Assemble the Jacobian matrix
!------------------------------------------------------	  
   IF(assemble_flag/=0) THEN
	   !coef_memb(nrl+1 :nrl+NDOF_ND+ndof_el,nrl+1:nrl+ndof_el)= ElemMass()
	   CALL ElemMass(ndof_el,elemJac)
   ELSE
	   !coef_memb(nrl+1 :nrl+NDOF_ND+ndof_el,nrl+1:nrl+ndof_el)= ElemJacobian(ndof_el,niter)
	   CALL ElemJacobian(ndof_el,niter,elemJac)
   ENDIF
   
   CALL Insert1DElement(nz_memb,elemJac,irn_memb,jcn_memb,coef_memb,nrl,nrl)
    
   nrl=nrl+ndof_el
    
ENDDO

IF(DEBUG)THEN
	WRITE(IOUT,*) 'MEMBER COEFFICIENT MATRIX'
	WRITE(IOUT,*)nz_memb
	DO j=1, nz_memb
		WRITE(IOUT,'(1X,2i5,ES20.12)')irn_memb(j),jcn_memb(j),coef_memb(j)
	ENDDO

ENDIF
 
END SUBROUTINE  AssembleMemberJacobian
!****************************************************************


!************************************************************
!*                                                          *                                      
!>  Extract member properties                               *
!*															*
!************************************************************ 
SUBROUTINE  ExtractMemberProperties(memb_no,memb_info_i,member,ncond_mb,mb_condition,distr_fun,error,init_cond)      

INTEGER,INTENT(IN)   ::memb_no  !<Number of the current member
INTEGER,INTENT(IN)   ::ncond_mb !<#ioaero::ncond_mb
INTEGER,INTENT(IN)   ::member(:,:)  !<#ioaero::member
TYPE (MemberInf),INTENT(IN)::memb_info_i    !< array containing the characteristics of a beam member
REAL(DBL),INTENT(IN) ::distr_fun(:,:)   !<#ioaero::distr_fun
TYPE(PrescriInf),INTENT(IN)::mb_condition(:) !<#ioaero::mb_condition
CHARACTER(*),INTENT(OUT)::error !<#ioaero::error
REAL(DBL),OPTIONAL,INTENT(IN)::init_cond(:,:)   !<#ioaero::init_cond

INTEGER:: ninit

ndiv=memb_info_i%ndiv
ncol_memb=memb_info_i%ncol_memb

! Distributed load information
!---------------------------------------
exist_load=.FALSE.
follower_load=.FALSE. ! Initialize the load flag and follower flag

IF(ncond_mb>0) CALL ExistPI(memb_no,mb_condition,exist_load,follower_load)
IF(exist_load) load=GetDistributedLoad(memb_no,mb_condition,distr_fun)
!------------------------------------------------------------------------------


! Obtain the initial conditions for this member
! for both initial step and time marching, it holds the values from previous time step
!---------------------------------------------
IF(init_flag/=0) THEN  
    ALLOCATE(init_memb(ndiv,12+NSTATES),STAT=allo_stat)
	IF(MemoryError('init_memb',error)) GOTO 9999
	ninit=SUM(member(1:memb_no-1,6))
	init_memb(1:ndiv,:)=init_cond(ninit+1:ninit+ndiv,:)
ENDIF

9999 IF(error/='')THEN
		WRITE(*,*) 'Member.f90 : ',error
        ERROR STOP 107
	 ENDIF

END SUBROUTINE ExtractMemberProperties
!****************************************************************

END MODULE Member
