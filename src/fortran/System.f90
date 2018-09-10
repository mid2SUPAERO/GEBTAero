!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************

!***************************************************************
!*                                                             *
!> This module assemles the system including the coefficient matrix (jacobian matrix) and the right hand side (negative of the equation values)                                     *
!* Outputs:coef, rhs                                           *
!*                                                             *
!***************************************************************
MODULE System

USE InternalData
USE PrescribedCondition
USE Member

IMPLICIT NONE

PRIVATE       ! So everything is private, except declared by PUBLIC

PUBLIC AssembleJacobian,AssembleRHS,ne,irn,jcn,coef

REAL(DBL)::x_pt(NSTRN)         !< the solution from the previous step for a key point.
INTEGER  ::kp_dof(NSTRN)       !< prescribed dof
REAL(DBL)::kp_cond(NSTRN)      !< prescribed value
INTEGER  ::kp_follower(NSTRN)  !< follower condition

INTEGER::ne !< Number of nonzero coefficients
INTEGER,ALLOCATABLE::irn(:) !< line index of nonzero coefficients
INTEGER,ALLOCATABLE::jcn(:) !< column index of nonzero coefficients
REAL(DBL),ALLOCATABLE::coef(:)  !< value of nonzero coefficients

!=============================================
CONTAINS


!************************************************************
!*                                                          *                                      
!>  Assemble the coefficient matrix of the beam system      *
!*															*
!************************************************************ 
SUBROUTINE AssembleJacobian(ndof_el,niter,memb_info,v_root_a,omega_a,member,error,&
			& ncond_mb,mb_condition,distr_fun,dof_con,x,aero_flag,grav_flag,init_cond)

INTEGER,INTENT(IN)::ndof_el !<#ioaero::ndof_el
INTEGER,INTENT(IN)::niter   !<#ioaero::niter
INTEGER,INTENT(IN)::aero_flag   !<#ioaero::aero_flag
INTEGER,INTENT(IN)::grav_flag   !<#ioaero::grav_flag
REAL(DBL),INTENT(IN)::v_root_a(:)   !< linear velocity of frame a
REAL(DBL),INTENT(IN)::omega_a(:)    !< angular velocity of frame a
TYPE (MemberInf),INTENT(IN)::memb_info(:)   !< contains the member parameters of the whole structure

REAL(DBL),INTENT(IN) ::distr_fun(:,:)   !<#ioaero::distr_fun
INTEGER,INTENT(IN)   ::member(:,:)  !<#ioaero::member
INTEGER,INTENT(IN)   ::ncond_mb !<#ioaero::ncond_mb
INTEGER              ::dof_con(:) !< note dof_con is passed by value, hence what is changed in this subroutine will not affect the original value.

REAL(DBL),INTENT(IN) :: x(:)    !< solution vector
TYPE(PrescriInf),INTENT(IN)::mb_condition(:)    !<#ioaero::mb_condition
CHARACTER(*),INTENT(OUT)::error !<#ioaero::error

REAL(DBL),OPTIONAL,INTENT(IN)::init_cond(:,:)   !<#ioaero::init_cond

INTEGER::i,j

INTEGER::pre_dof !< prescibed dof of the connection point

REAL(DBL),ALLOCATABLE::coef_memb(:) !< coefficient matrix of a member

INTEGER:: nrow,ncol,nrow_end,ncol_end
INTEGER,ALLOCATABLE::irn_memb(:),jcn_memb(:)

INTEGER::nz_memb

REAL(DBL)::coef_tm(NDOF_ND,ndof_el)


ALLOCATE(coef(neMax),STAT=allo_stat)
IF(MemoryError('coef',error)) GOTO 9999


ALLOCATE(irn(neMax),STAT=allo_stat)
IF(MemoryError('irn',error)) GOTO 9999


ALLOCATE(jcn(neMax),STAT=allo_stat)
IF(MemoryError('jcn',error)) GOTO 9999


ne=0

DO i=1, SIZE(member,1)

   CALL ExtractMemberProperties(i,memb_info(i),member,ncond_mb,mb_condition,distr_fun,error,init_cond)      
      
   ALLOCATE(coef_memb(ndiv*nzElemMax),STAT=allo_stat)
   IF(MemoryError('coef_memb',error)) GOTO 9999
   ALLOCATE(irn_memb(ndiv*nzElemMax),STAT=allo_stat)
   IF(MemoryError('irn_memb',error)) GOTO 9999
   ALLOCATE(jcn_memb(ndiv*nzElemMax),STAT=allo_stat)
   IF(MemoryError('jcn_memb',error)) GOTO 9999


   ! Dealing with the starting point
   !-----------------------------------
   !Extract the solution from the previous step for the member, excluding the terminal points
   !If the starting point is not used yet, then the first NSTRN values belong to the point
   !If the starting point has been used already in the assembly, then x_memb starts with ncol+1
   !------------------------------------------------------------------------------------------
   nrow=index_mb(i,1)  ! the starting row of the member block excluding the first blocking connecting to the trailing point
   ncol=index_mb(i,2)  ! the starting column of the member block
   nrow_end=nrow+ncol_memb-NDOF_ND  ! the starting row of the ending block of the member
   ncol_end=ncol+ncol_memb-ndof_el  ! the starting column of the ending block of the member
  
!------------------------------------------------------------------------------------------
   CALL AssembleMemberJacobian(ndof_el,niter,memb_info(i),v_root_a,omega_a,x(ncol:ncol+ncol_memb-1),&
                           & nz_memb,irn_memb,jcn_memb,coef_memb,aero_flag,grav_flag) ! assemble the coefficient matrix within the member
  
   CALL  Extract2DElement(nz_memb,irn_memb,jcn_memb,coef_memb,coef_tm,0,0)
   CALL AssemblePointJ(-1,member(i,1),nrow-NSTRN,ncol) ! assemble the Jacobian related with the starting point

   
   !insert the member coef into the global matrices excluding the starting portion of the first element
   !and the ending portion of the last element. If there is only one element for static analysis, it is not necessary
   !-------------------------------------------------------------------------------------------------------------------
!	IF(ncol_memb>NDOF_ND) coef(nrow:nrow_end-1,ncol:ncol+ncol_memb-1)=coef_memb(NDOF_ND+1:ncol_memb,:) 

	IF(ncol_memb>NDOF_ND) THEN
	   DO j=1,nz_memb
	      IF(irn_memb(j)>NDOF_ND.AND.irn_memb(j)<=ncol_memb) THEN
		      ne=ne+1
		      irn(ne)=nrow-1+irn_memb(j)-NDOF_ND
			  jcn(ne)=ncol-1+jcn_memb(j)
              coef(ne)=coef_memb(j)
		  ENDIF
	   ENDDO
	ENDIF


   ! Dealing with the ending point
   !-----------------------------------
   CALL  Extract2DElement(nz_memb,irn_memb,jcn_memb,coef_memb,coef_tm,ncol_memb,ncol_memb-ndof_el)

   CALL AssemblePointJ(1,member(i,2),nrow_end,ncol_end) 


   DEALLOCATE(jcn_memb,irn_memb,coef_memb)
   IF(init_flag/=0) DEALLOCATE(init_memb)

ENDDO


! If some connection points have prescribed displacements, coef should be modified correspondingly
! know the corresponding equilibrium equation will be eliminated and replaced with an identity equation
! for the prescribed displacement
!-------------------------------------------------------------------------------------------------------------------
IF(assemble_flag==0) THEN
	DO i=1,SIZE(dof_con)
	  IF(dof_con(i)/=0) THEN ! a connection point
		kp_dof=dof_all(i,:)

		nrow=index_kp(i,1)-1
        ncol=index_kp(i,2)-1
        DO j=1,NSTRN
		    pre_dof=kp_dof(j)
			IF(pre_dof>0.AND.pre_dof<=NSTRN) THEN ! a displacement/rotation condition
			  !coef(nrow+pre_dof,:)=0.0D0
			  WHERE(irn==nrow+pre_dof) coef=0.0D0
			  !coef(nrow+pre_dof,ncol+pre_dof)=1.0D0
			  ne=ne+1
              irn(ne)=nrow+pre_dof
			  jcn(ne)=ncol+pre_dof
			  coef(ne)=1.0D0
			ENDIF
        ENDDO
	ENDIF

  ENDDO
ENDIF

IF(DEBUG)THEN
	WRITE(IOUT,*) 'COEFFICIENT MATRIX'
	WRITE(IOUT,*)ne
	DO j=1, ne
		WRITE(IOUT,'(1X,2i5,ES20.12)')irn(j),jcn(j),coef(j)
	ENDDO
ENDIF

9999 IF(error/='')THEN
    WRITE(0,*) 'System.f90 : ',error
    ERROR STOP 125
ENDIF

CONTAINS

!> Assemble the Jacobian related with trailing points       														
!-----------------------------------------------------------------
SUBROUTINE  AssemblePointJ(flag,kp,nrow_mb,ncol_mb) 
     
INTEGER,INTENT(IN)::flag !< indicating whether it is the starting point or the ending point: -1-starting, 1-ending
INTEGER,INTENT(IN)::kp !< the point number
INTEGER,INTENT(IN)::nrow_mb,ncol_mb

INTEGER::nrow_kp,ncol_kp

REAL(DBL)::ekttek(NDIM,NDIM,NDIM)   ! e_k.\theta^T+ \theta.e_k^T
REAL(DBL)::theta(NDIM)
REAL(DBL)::eCTtheta(NDIM,NDIM,NDIM)  !derivatives of C^T w.r.t theta, eCTtheta(k,:,:)= d_C^T/d_thetak
INTEGER:: i,j

!Extract point information
!----------------------------
nrow_kp=index_kp(kp,1)  ! starting row for the point
ncol_kp=index_kp(kp,2)  ! starting column for the point

x_pt=x(ncol_kp:ncol_kp+5) ! the point variables
kp_dof=dof_all(kp,:) ! the prescribed dof
kp_cond=cond_all(kp,:)  ! the prescribed values
kp_follower=follower_all(kp,:)  ! whether it is a follower

! Deal with follower conditions
!-----------------------------------------
IF(ANY(kp_follower==1)) THEN
     IF(kp_dof(1)==7)  THEN 
		theta=x_pt(4:6)
	 ELSE IF (kp_dof(1)==10) THEN
        theta=x_pt(1:3)
	 ENDIF
	 CALL CT_THETA(theta,DirCosineTRodrigues(theta),ekttek,eCTtheta)
ENDIF
   

IF(dof_con(kp)==0) THEN  !If the trailing point is a boundary point, 

!	coef(nrow_kp:nrow_kp+11,ncol_mb:ncol_mb+ndof_el-1)=coef_tm
	CALL Insert1DElement(ne,coef_tm,irn,jcn,coef,nrow_kp-1,ncol_mb-1)
	IF(assemble_flag==0)THEN
		DO j=1,NSTRN
			pre_dof=kp_dof(j)
			IF(pre_dof/=0) THEN
				IF(pre_dof>NSTRN) THEN
					!coef(nrow_kp+pre_dof-1,ncol_kp+j-1)=flag ! force type conditions
					ne=ne+1
					irn(ne)=nrow_kp+pre_dof-1
					jcn(ne)=ncol_kp+j-1
					coef(ne)=flag
				 ENDIF
				IF(pre_dof<=NSTRN) THEN
					!coef(nrow_kp+pre_dof-1,ncol_kp+j-1)=-flag ! displacement type conditions
					ne=ne+1
					irn(ne)=nrow_kp+pre_dof-1
					jcn(ne)=ncol_kp+j-1
					coef(ne)=-flag
				ENDIF
			ENDIF
		ENDDO

		IF(ANY(kp_follower==1))  CALL PointFollowerJ(flag,nrow_kp-1,ncol_kp-1,eCTtheta)  ! add the contribution due to follower forces/moments
	ENDIF


ELSE IF(dof_con(kp)==1) THEN ! the trailing point  is a connection point 
	
	    dof_con(kp)=2 ! indicated the connection point has been used
	
		!coef(nrow_kp:nrow_kp+11,ncol_mb:ncol_mb+ndof_el-1)=coef_tm
		CALL Insert1DElement(ne,coef_tm,irn,jcn,coef,nrow_kp-1,ncol_mb-1)

		IF(assemble_flag==0)THEN
			DO i=1,NSTRN
			    ne=ne+1
				!coef(nrow_kp+5+i,ncol_kp+i-1)=flag  
				irn(ne)=nrow_kp+5+i
				jcn(ne)=ncol_kp+i-1
				coef(ne)=flag  
			ENDDO

			IF(ANY(kp_dof>NSTRN).AND.ANY(kp_follower==1))&
			& CALL PointFollowerJ(1,nrow_kp-1,ncol_kp-1,eCTtheta)  ! add the contribution due to follower forces/moments
		ENDIF

ELSE IF(dof_con(kp)==2) THEN

	IF(assemble_flag==0)THEN
		DO i=1,NSTRN
			!coef(nrow_mb+i-1,ncol_kp+i-1)=flag
			ne=ne+1
			irn(ne)=nrow_mb+i-1
			jcn(ne)=ncol_kp+i-1
			coef(ne)=flag  
		ENDDO
	ENDIF
   	!coef(nrow_kp:nrow_kp+5,ncol_mb:ncol_mb+ndof_el-1)=coef_tm(:NSTRN,:)
	CALL Insert1DElement(ne,coef_tm(:NSTRN,:),irn,jcn,coef,nrow_kp-1,ncol_mb-1)
	!coef(nrow_mb:nrow_mb+5,ncol_mb:ncol_mb+ndof_el-1)=coef_tm(7:NDOF_ND,:)
	CALL Insert1DElement(ne,coef_tm(7:NDOF_ND,:),irn,jcn,coef,nrow_mb-1,ncol_mb-1)
ENDIF       


END SUBROUTINE  AssemblePointJ      
!----------------------------------------------------------

END SUBROUTINE AssembleJacobian
!***********************************************************





!************************************************************
!*                                                          *                                      
!>  Assemble the right hand side
!*															*
!************************************************************ 
SUBROUTINE AssembleRHS(ndof_el,memb_info,v_root_a,omega_a,member,error,&
			& ncond_mb,mb_condition,distr_fun,dof_con,x,rhs,aero_flag,grav_flag,init_cond)

INTEGER,INTENT(IN)::ndof_el !<#ioaero::ndof_el
INTEGER,INTENT(IN)::aero_flag    !<#ioaero::aero_flag
INTEGER,INTENT(IN)::grav_flag   !<#ioaero::grav_flag
TYPE (MemberInf),INTENT(IN)::memb_info(:)   !< contains the member parameters of the whole structure

REAL(DBL),INTENT(IN)::v_root_a(:)   !< linear velocity of frame a
REAL(DBL),INTENT(IN)::omega_a(:)    !< angular velocity of frame a
REAL(DBL),INTENT(IN) ::distr_fun(:,:)   !<#ioaero::distr_fun
INTEGER,INTENT(IN)   ::member(:,:)  !<#ioaero::member
INTEGER,INTENT(IN)   ::ncond_mb     !<#ioaero::ncond_mb
INTEGER              ::dof_con(:) 
REAL(DBL),INTENT(IN) :: x(:)    !<solution vector
TYPE(PrescriInf),INTENT(IN)::mb_condition(:)    !<#ioaero::mb_condition

REAL(DBL),INTENT(OUT)::rhs(:)
CHARACTER(*),INTENT(OUT)::error

REAL(DBL),OPTIONAL,INTENT(IN)::init_cond(:,:)   !<#ioaero::init_cond

REAL(DBL),ALLOCATABLE::rhs_memb(:)    !< rhs for a member


INTEGER::i,j
INTEGER:: nrow,ncol,nrow_end
INTEGER::pre_dof !< prescibed dof of the connection point

rhs=0.0D0

DO i=1, SIZE(member,1)
  
   CALL ExtractMemberProperties(i,memb_info(i),member,ncond_mb,mb_condition,distr_fun,error,init_cond)      
      
   ALLOCATE(rhs_memb(ncol_memb+NDOF_ND),STAT=allo_stat)
   IF(MemoryError('rhs_memb',error)) GOTO 9999
 

   nrow=index_mb(i,1)  ! the starting row of the member block
   ncol=index_mb(i,2)  ! the starting column of the member block
   nrow_end=nrow+ncol_memb-NDOF_ND  ! the starting row of the ending block of the member
      
   CALL AssembleMemberRHS(ndof_el,memb_info(i),v_root_a,omega_a,x(ncol:ncol+ncol_memb-1),rhs_memb,aero_flag,grav_flag,init_cond)
  
   CALL AssemblePointRHS(-1,member(i,1),rhs_memb(:NDOF_ND),nrow-NSTRN) ! Assemble the starting point
 
   !insert the member rhs into the global matrices excluding the starting portion of the first element
   !and the ending portion of the last element. If there is only one element within the member
   ! for static analysis, it is not necessary
   !-------------------------------------------------------------------------------------------------------------------
   IF(ncol_memb>NDOF_ND) rhs(nrow:nrow_end-1)=rhs_memb(NDOF_ND+1:ncol_memb)     
  
   CALL AssemblePointRHS(1,member(i,2),rhs_memb(ncol_memb+1:),nrow_end) ! assemble the ending point

   DEALLOCATE(rhs_memb)
   IF(init_flag/=0) DEALLOCATE(init_memb)


ENDDO

! If some connection points have prescribed displacements/rotations, rhs should be modified correspondingly
!-------------------------------------------------------------------------------------------------------------------
DO i=1,SIZE(dof_con)
	IF(dof_con(i)/=0) THEN ! a connection point
		kp_dof=dof_all(i,:)
        kp_cond=cond_all(i,:)
        
		nrow=index_kp(i,1)-1
        DO j=1,NSTRN
		    pre_dof=kp_dof(j)
			IF(pre_dof>0.AND.pre_dof<=NSTRN) rhs(nrow+pre_dof)=kp_cond(j)			     
        ENDDO
	ENDIF

ENDDO


IF(DEBUG)THEN
	WRITE(IOUT,*) 'RIGHT HAND SIDE'
	write(iout,*)nsize
	DO j=1, nsize
		WRITE(IOUT,'(1X,i5,1ES20.12)')j,	rhs(j)
	ENDDO	
ENDIF



9999 IF(error/='')THEN
    WRITE(0,*) 'System.f90 : ',error
    ERROR STOP 126
ENDIF


CONTAINS

!>  Assemble the Jacobian related with trailing points      
!----------------------------------------------------------
SUBROUTINE  AssemblePointRHS(flag,term_pt,rhs_tm,nrow_mb) 
     
INTEGER,INTENT(IN)::flag !< indicating whether it is the starting point or the ending point: -1-starting, 1-ending
INTEGER,INTENT(IN)::term_pt !< the point number
REAL(DBL),INTENT(IN)::rhs_tm(:) !< the right hand side of the member associated with the end point
INTEGER,INTENT(IN)::nrow_mb

INTEGER::nrow_kp,ncol_kp
INTEGER::free_dof,for_dof
INTEGER:: i,j

nrow_kp=index_kp(term_pt,1)  ! starting row for the point
ncol_kp=index_kp(term_pt,2)  ! starting column for the point

x_pt=x(ncol_kp:ncol_kp+5) ! the point variables
kp_dof=dof_all(term_pt,:) ! the prescribed dof
kp_cond=cond_all(term_pt,:)  ! the prescribed values

IF(ANY(follower_all(term_pt,:)==1)) &
& kp_cond=UpdateFollower(kp_dof,follower_all(term_pt,:),kp_cond,x_pt)
   
 
!If the trailing point is a boundary point, 
!insert the corresponding portion of the attaching  element 
!-------------------------------------------------------- 
IF(dof_con(term_pt)==0) THEN
	DO j=1,NSTRN
	
		pre_dof=kp_dof(j)
		IF(pre_dof>NSTRN) THEN ! force type conditions
	
			free_dof=pre_dof-NSTRN       
			rhs(nrow_kp+free_dof-1)=rhs(nrow_kp+free_dof-1)+kp_cond(j)  ! prescribed conditions are external force/moments at the boundary 
			rhs(nrow_kp+pre_dof-1)=rhs(nrow_kp+pre_dof-1)-flag*x_pt(j)

		ELSE IF(pre_dof<=NSTRN) THEN  ! displacement type conditions

			free_dof=pre_dof+NSTRN	
			rhs(nrow_kp+free_dof-1)=rhs(nrow_kp+free_dof-1)-kp_cond(j)*flag
 			rhs(nrow_kp+pre_dof-1)=rhs(nrow_kp+pre_dof-1)+flag*x_pt(j) ! solved values are interal forces/moments 

		ENDIF
	ENDDO

	rhs(nrow_kp:nrow_kp+11)=rhs(nrow_kp:nrow_kp+11)+rhs_tm

ELSE IF(dof_con(term_pt)==1) THEN ! it is a connection point not used yet

    dof_con(term_pt)=2 ! indicated the connection point has been used

	DO i=1,NSTRN
	! assemble force type conditions into rhs
	!----------------------------------------------
		for_dof=kp_dof(i)-NSTRN
		IF(for_dof>0) rhs(nrow_kp+for_dof-1)=rhs(nrow_kp+for_dof-1)+kp_cond(i)
	ENDDO

	rhs(nrow_kp+NSTRN:nrow_kp+11)=rhs(nrow_kp+NSTRN:nrow_kp+11)-flag*x_pt  ! correction from the calculated displacement/rotations at the rotation point
	rhs(nrow_kp:nrow_kp+11)=rhs(nrow_kp:nrow_kp+11)+rhs_tm

ELSE 
 
	rhs(nrow_kp:nrow_kp+5)=rhs(nrow_kp:nrow_kp+5)+rhs_tm(1:NSTRN)

	rhs(nrow_mb:nrow_mb+5)=rhs(nrow_mb:nrow_mb+5)+rhs_tm(NSTRN+1:NDOF_ND)-flag*x_pt

ENDIF       
  

END SUBROUTINE  AssemblePointRHS      
!---------------------------------------

END SUBROUTINE AssembleRHS
!***********************************************************




!************************************************************
!*                                                          *
!>  Add the contribution to Jacobian matrix due to follower
!! point force or moments
!*                                                          * 
!************************************************************
SUBROUTINE PointFollowerJ(flag,nrow,ncol,eCTtheta)

INTEGER,INTENT(IN)::flag
INTEGER,INTENT(IN)::nrow,ncol
REAL(DBL),INTENT(IN)::eCTtheta(:,:,:)


INTEGER:: str_row_fol !< starting row for inserting the jacobian due to follower conditions

IF(ANY(kp_follower(1:3)==1)) THEN
   IF(kp_dof(1)==7) str_row_fol=nrow
   IF(kp_dof(1)==10) str_row_fol=nrow+3 
!   coef(str_row_fol+1:str_row_fol+3,ncol+4:ncol+6)=&
!  & coef(str_row_fol+1:str_row_fol+3,ncol+4:ncol+6)-flag*FollowerJ(kp_follower(1:3),flag*kp_cond(1:3),eCTtheta)
  CALL Insert1DElement(ne,-flag*FollowerJ(kp_follower(1:3),flag*kp_cond(1:3),eCTtheta),&
                & irn,jcn,coef,str_row_fol,ncol+3)

 !& coef(str_row_fol+1:str_row_fol+3,ncol+4:ncol+6)-flag*FollowerJ(kp_follower(1:3),kp_cond(1:3),eCTtheta)
ENDIF

IF(ANY(kp_follower(4:6)==1)) THEN
   IF(kp_dof(4)==7) str_row_fol=nrow
   IF(kp_dof(4)==10) str_row_fol=nrow+3 
!   coef(str_row_fol+1:str_row_fol+3,ncol+4:ncol+6)=&
!   & coef(str_row_fol+1:str_row_fol+3,ncol+4:ncol+6)-flag*FollowerJ(kp_follower(4:6),flag*kp_cond(4:6),eCTtheta)
   !coef(str_row_fol+1:str_row_fol+3,ncol+4:ncol+6)=-flag*FollowerJ(kp_follower(4:6),flag*kp_cond(4:6),eCTtheta)
   CALL Insert1DElement(ne,-flag*FollowerJ(kp_follower(4:6),flag*kp_cond(4:6),eCTtheta),&
                     & irn,jcn,coef,str_row_fol,ncol+3)
! & coef(str_row_fol+1:str_row_fol+3,ncol+4:ncol+6)-flag*FollowerJ(kp_follower(4:6),kp_cond(4:6),eCTtheta)
ENDIF


END SUBROUTINE PointFollowerJ
!************************************************************




END MODULE System
