!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************

!***************************************************************
!*                                                             *
! This module contains information and calculation for an     *
!* element within a member 
!* Outputs: ElemEqn,ElemJacobian                               *
!***************************************************************
!>This module contains information and calculation for an element within a member 
MODULE Element

USE InternalData
USE PrescribedCondition
USE GlobalDataFun

IMPLICIT NONE

PRIVATE       ! So everything is private, except declared by PUBLIC

PUBLIC ElemEqn,ElemJacobian,ElemMass,ExtractElementProperties
PUBLIC load,exist_load,follower_load

! private data needed inside this module
!--------------------------------------------
REAL(DBL)::dL              !< length of the element
REAL(DBL)::Le              !< the ending arc length of the current element
REAL(DBL)::eCab(NDIM,NDIM) !< direction cosine matrix of the undeformed element
REAL(DBL)::eFlex(NSTRN,NSTRN) !<flexibility matrix of the elment
REAL(DBL)::eMass(NSTRN,NSTRN) !< inverse of the mass matrix of the element
TYPE(DistriLoad)::load        !< distributed load
LOGICAL::exist_load,follower_load !< flags to indicate whether distributed load exist and whether they are follower forces

REAL(DBL)::Ui(NDIM),theta(NDIM),Fi(NDIM),Mi(NDIM),e1GammaD(NDIM),kappa(NDIM),ePi(NDIM),Hi(NDIM),Vi(NDIM),OMEGAi(NDIM) ! e1GammaD=e1+GammaD
REAL(DBL)::ev_i(NDIM)    !< initial velocity of the mid point of the element
REAL(DBL)::eOmega_a(NDIM) !< initial angular velocity of the element

REAL(DBL)::eCT(NDIM,NDIM)         !< the transpose of the direction cosine matrix corresponding to elastic rotation
REAL(DBL)::eCTCAB(NDIM,NDIM)      !< eCT.Cab
REAL(DBL)::eCabhalfL(NDIM,NDIM)   !< eCab*dL/2
REAL(DBL)::eCTCabhalfL(NDIM,NDIM) !< eCTCab*dL/2

REAL(DBL)::UiDot(NDIM),ThetaDot(NDIM),CTCabPdot(NDIM),CTCabHdot(NDIM)

! data needed for aerodynamic forces
REAL(DBL)::Rho,Chord,X_cg,alpha,alphadot,hdot,aw,bw,U,hdotdot,alphadotdot,alpha_ac,beta_ac,beta
INTEGER  :: a_flag, g_flag
REAL(DBL)            :: lambda(NSTATES),lambdadot(NSTATES),P(NSTATES,2*NSTATES+2) ! induces-flow states
REAL(DBL)            :: Ident(NSTATES,NSTATES)
REAL(DBL)            :: lambda0 ! induces-flow states
REAL(DBL)            :: A(NSTATES,NSTATES),b(NSTATES),c(NSTATES)
REAL(DBL)            :: dir_moment(NDIM),dir_lift(NDIM),Jdir_moment(NDIM,NDIM),Jdir_lift(NDIM,NDIM)
REAL(DBL)            ::Jalpha_THETA(NDIM),Jalphadot_THETA(NDIM),Jhdot_THETA(NDIM),Jalphadot_PH(NDIM+NDIM),Jhdot_PH(NDIM+NDIM)
REAL(DBL)            ::Jalphadotdot_PHDOT(NDIM+NDIM),Jhdotdot_PHDOT(NDIM+NDIM)
REAL(DBL)            ::Jalphadotdot_PH(NDIM+NDIM),Jhdotdot_PH(NDIM+NDIM)
REAL(DBL)            ::Jlambda0_PHDOT(NDIM+NDIM),Jlambda0_LAMBDA(NSTATES)
REAL(DBL)            ::Jlift(3,18+NSTATES),Jmoment(3,18+NSTATES)
REAL(DBL)            ::eCaf(NDIM,NDIM),Wind(NDIM)

!=============================================

CONTAINS


!************************************************************
!*                                                          *
!* Evaluate nonlinear equations for each element, including *
!* the equation value due to x_0 and the term due to        *
!* distributed load. Note, it is put in the right hand side *
!*															*
!************************************************************ 
!> Compute the value of the Right Hand Side for a finite element
FUNCTION ElemEqn(ndof_el)

IMPLICIT NONE

INTEGER,INTENT(IN)  :: ndof_el !<#ioaero::ndof_el

REAL(DBL)::ElemEqn(ndof_el+NDOF_ND) ! the functional value for each element
REAL(DBL):: tmpR(NDIM),tmp,tmpN(NSTATES)

!initialisation
ElemEqn=0.D0

ElemEqn(1:3)=MATMUL(eCTCab,Fi)
ElemEqn(ndof_el+1:ndof_el+3)=-ElemEqn(1:3)

ElemEqn(4:6)= MATMUL(eCTCabhalfL,CrossProduct(e1GammaD,Fi))
tmpR=MATMUL(eCTCab,Mi)
ElemEqn(ndof_el+4:ndof_el+6)=-tmpR+ElemEqn(4:6)
ElemEqn(4:6)=tmpR+ElemEqn(4:6)

ElemEqn(7:9)=MATMUL(eCTCabhalfL,e1GammaD)-eCabhalfL(:,1)
ElemEqn(ndof_el+7:ndof_el+9)=Ui+ElemEqn(7:9)
ElemEqn(7:9)=-Ui+ElemEqn(7:9)

tmpR=MATMUL(eCabhalfL,kappa)
ElemEqn(10:12)=tmpR+CrossProduct(theta*0.5d0,tmpR)+MATMUL(OuterProduct(theta,theta),tmpR*0.25d0)
ElemEqn(ndof_el+10:ndof_el+12)=theta+ElemEqn(10:12)
ElemEqn(10:12)=-theta+ElemEqn(10:12)

IF(exist_load) THEN
	ElemEqn(1:6)=ElemEqn(1:6)+GetLoad(-1,dL,Le,eCT,load,follower_load)  ! add the contribution due to distributed load
	ElemEqn(ndof_el+1:ndof_el+6)=ElemEqn(ndof_el+1:ndof_el+6)+GetLoad(1,dL,Le,eCT,load,follower_load)  ! add the contribution due to distributed load
ENDIF

! for dynamic analysis
!------------------------------------------
IF(ndof_el>=18) THEN

	! modify the right hand side for dynamic analysis including 
	! both steady state and initial conditions and time marching
	!--------------------------------------------------
	
	tmpR=CrossProduct(eOmega_a,MATMUL(eCTCabhalfL,ePi))
	IF(init_flag==2) tmpR=tmpR+two_divide_dt*MATMUL(eCTCabhalfL,ePi)	
	ElemEqn(1:3)=ElemEqn(1:3)-tmpR
    ElemEqn(ndof_el+1:ndof_el+3)=ElemEqn(ndof_el+1:ndof_el+3)-tmpR
    
	tmpR=CrossProduct(eOmega_a,MATMUL(eCTCabhalfL,Hi))+MATMUL(eCTCabhalfL,CrossProduct(Vi,ePi))
    IF(init_flag==2) tmpR=tmpR+two_divide_dt*MATMUL(eCTCabhalfL,Hi)
	ElemEqn(4:6)=ElemEqn(4:6)-tmpR
    ElemEqn(ndof_el+4:ndof_el+6)=ElemEqn(ndof_el+4:ndof_el+6)-tmpR

	! the functions of fpi and fhi
	!------------------------------------------
	ElemEqn(NDOF_ND+1:NDOF_ND+3)=ev_i+CrossProduct(eOmega_a,Ui)-MATMUL(eCTCab,Vi)
    IF(init_flag==2) ElemEqn(NDOF_ND+1:NDOF_ND+3)=ElemEqn(NDOF_ND+1:NDOF_ND+3)+two_divide_dt*Ui
	ElemEqn(NDOF_ND+4:NDOF_ND+6)=MATMUL(TRANSPOSE(eCTCab),eOmega_a)-OMEGAi
ENDIF

IF(init_flag/=0) THEN !needed for the initial step and time marching
	ElemEqn(1:3)=ElemEqn(1:3)-dL*0.5D0*CTCabPdot
	ElemEqn(ndof_el+1:ndof_el+3)=ElemEqn(ndof_el+1:ndof_el+3)-dL*0.5D0*CTCabPdot
    
	ElemEqn(4:6)=ElemEqn(4:6)-dL*0.5D0*CTCabHdot
	ElemEqn(ndof_el+4:ndof_el+6)=ElemEqn(ndof_el+4:ndof_el+6)-dL*0.5D0*CTCabHdot

    ElemEqn(NDOF_ND+1:NDOF_ND+3)=ElemEqn(NDOF_ND+1:NDOF_ND+3)+UiDot
	ElemEqn(NDOF_ND+4:NDOF_ND+6)=ElemEqn(NDOF_ND+4:NDOF_ND+6)+ &
	&   MATMUL(TRANSPOSE(eCab),(ThetaDot-CrossProduct(theta*.5D0,ThetaDot))/(1+DOT_PRODUCT(theta,theta)*0.25D0) )
ENDIF

! quasi steady aerodynamics
IF (a_flag == 1) THEN
	!lift
	tmp = dL*PI*Rho*bw*U*(hdot+bw*(0.5-aw)*alphadot+U*alpha)
	tmpR = tmp*dir_lift
	ElemEqn(1:3)=ElemEqn(1:3) + tmpR
	ElemEqn(ndof_el+1:ndof_el+3)=ElemEqn(ndof_el+1:ndof_el+3)+ tmpR
	!moment
	tmp = bw*(0.5+aw)*tmp -0.5*pi*Rho*U*bw**3*alphadot*0.5*dL
	tmpR = tmp*dir_moment
	ElemEqn(4:6)=ElemEqn(4:6)+ tmpR
	ElemEqn(ndof_el+4:ndof_el+6)=ElemEqn(ndof_el+4:ndof_el+6)+ tmpR

! quasi steady aerodynamics with added mass
ELSEIF (a_flag >= 2) THEN
	!case a_flag = 2 : added mass only
	IF (a_flag == 2) lambda0 = 0
		
	!lift
	tmp = 0.5*dL*PI*Rho*bw**2*(hdotdot-bw*aw*alphadotdot+U*alphadot) + &
	& dL*PI*Rho*bw*U*(hdot+bw*(0.5-aw)*alphadot+U*alpha-lambda0)
	tmpR = tmp*dir_lift
	ElemEqn(1:3)=ElemEqn(1:3)+ tmpR
	ElemEqn(ndof_el+1:ndof_el+3)=ElemEqn(ndof_el+1:ndof_el+3)+ tmpR
	!moment
	tmp = bw*(0.5+aw)*tmp-0.5*dL*PI*Rho*bw**3*(0.5*hdotdot+U*alphadot+bw*(0.125-0.5*aw)*alphadotdot)
	tmpR = tmp*dir_moment
	ElemEqn(4:6)=ElemEqn(4:6)+ tmpR
	ElemEqn(ndof_el+4:ndof_el+6)=ElemEqn(ndof_el+4:ndof_el+6)+ tmpR
	
	IF (a_flag ==3) THEN
		! induced-flow equation
		tmpN = -bw*MATMUL(A,lambdadot)-U*lambda+bw*(hdotdot+U*alphadot+bw*(0.5-aw)*alphadotdot)*c
		ElemEqn(ndof_el-NSTATES+1:ndof_el) = ElemEqn(ndof_el-NSTATES+1:ndof_el) + tmpN
	ENDIF 
	
ENDIF

!with gravity
IF (g_flag /= 0) THEN
	tmpR = GRAV/eMass(3,3)*0.5*dL*(/0.D0,0.D0,1.D0/)
	ElemEqn(1:3)=ElemEqn(1:3)+ tmpR
	ElemEqn(ndof_el+1:ndof_el+3)=ElemEqn(ndof_el+1:ndof_el+3)+ tmpR

	tmpR = -GRAV/eMass(3,3)*X_cg*0.5*dL*(/0.D0,0.D0,1.D0/)
	ElemEqn(4:6)=ElemEqn(4:6)+ tmpR
	ElemEqn(ndof_el+4:ndof_el+6)=ElemEqn(ndof_el+4:ndof_el+6)+ tmpR
ENDIF

END FUNCTION ElemEqn
!***************************************************


!************************************************************
!*                                                          *
!>  Caculate the Jacobian matrix for each element
!*                                                          *
!*															*
!************************************************************ 
SUBROUTINE ElemJacobian(ndof_el,niter,elemJac)

IMPLICIT NONE

INTEGER,INTENT(IN)  :: ndof_el !< #ioaero::ndof_el
INTEGER,INTENT(IN)  :: niter !< #ioaero::niter

REAL(DBL),INTENT(OUT)::ElemJac(:,:) !< the coefficient matrix for each element

REAL(DBL)::ekttek(NDIM,NDIM,NDIM)   ! e_k.\theta^T+ \theta.e_k^T
REAL(DBL)::eCTtheta(NDIM,NDIM,NDIM),eCTtheta2(NDIM,NDIM,NDIM)  !derivatives of C^T w.r.t theta, eCTtheta(k,:,:)= d_C^T/d_thetak
REAL(DBL)::temp6(NDIM+NDIM),temp31(NDIM),temp32(NDIM),temp36(nDIM,NDIM+NDIM),temp333(NDIM,NDIM,NDIM) ! temporary arrays 
REAL(DBL)::tempN6(NSTATES,6),tempNN(NSTATES,NSTATES),temp33(NDIM,NDIM),tempN(NSTATES)

REAL(DBL)::tilde_omega(NDIM,NDIM),tt4

INTEGER:: i

ElemJac=0.0D0

CALL CT_THETA(theta,eCT,ekttek,eCTtheta)  ! obtain derivatives of C^T w.r.t theta

! d_fu/d_theta
!--------------------
temp31=MATMUL(eCab,Fi)
ElemJac(1:3,4:6)=-MATMUL3(eCTtheta,temp31)   ! coefficient of theta of Fu 
ElemJac(ndof_el+1:ndof_el+3,4:6)=-ElemJac(1:3,4:6)

! d_fu/d_F
!--------------------		
ElemJac(1:3,7:9)=-eCTCab  
ElemJac(ndof_el+1:ndof_el+3,7:9)=eCTCab

! d_fGammaD/d_theta
!--------------------
temp31=MATMUL(eCab,Mi)
temp32=-MATMUL(eCab,CrossProduct(e1GammaD,dL*0.5d0*Fi))
ElemJac(4:6,4:6)=MATMUL3(eCTtheta,temp32-temp31)  ! coefficient of theta of FGammaD1 
ElemJac(ndof_el+4:ndof_el+6,4:6)=MATMUL3(eCTtheta,temp32+temp31) 

! d_fGammaD/d_F
!--------------------
ElemJac(4:6,7:9)=-MATMUL(eCTCabhalfL,Tilde(e1GammaD)-MATMUL(Tilde(Fi),eFlex(1:3,1:3)))  
ElemJac(ndof_el+4:ndof_el+6,7:9)=ElemJac(4:6,7:9)

! d_fGammaD/d_M
!--------------------
temp33=MATMUL(eCTCabhalfL,MATMUL(Tilde(Fi),eFlex(1:3,4:6)))	
ElemJac(4:6,10:12)=temp33-eCTCab
ElemJac(ndof_el+4:ndof_el+6,10:12)=temp33+eCTCab

! d_fF/d_u
!--------------------
ElemJac(7:9,1:3)=I3 
ElemJac(ndof_el+7:ndof_el+9,1:3)=-I3 

! d_fF/d_theta
!--------------------
temp31=-MATMUL(eCabhalfL,e1GammaD)
ElemJac(7:9,4:6)=MATMUL3(eCTtheta,temp31)   
ElemJac(ndof_el+7:ndof_el+9,4:6)=ElemJac(7:9,4:6)


! d_fF/d_F,d_fF/d_M 
!--------------------
ElemJac(7:9, 7:12)=-MATMUL(eCTCabhalfL,eFlex(1:3,:)) 
ElemJac(ndof_el+7:ndof_el+9,7:12)=ElemJac(7:9, 7:12)

! d_fM/d_theta
!--------------------
DO i=1,NDIM
	temp333(i,:,:)=Tilde(I3(i,:))*0.5d0+ekttek(i,:,:)*0.25d0  ! \tilde{ek}/2+ekttek/4
ENDDO

temp33=MATMUL3(temp333,MATMUL(eCabhalfL,kappa))
ElemJac(10:12,4:6)=I3-temp33
ElemJac(ndof_el+10:ndof_el+12,4:6)=-I3-temp33


! d_fM/d_F,d_fM/d_M 
!--------------------
ElemJac(10:12,7:12)=-MATMUL(MATMUL(I3+Tilde(theta)*.5d0+OuterProduct(theta,theta)*.25d0,eCabhalfL),& 
                        & eFlex(4:6,:)) 
ElemJac(ndof_el+10:ndof_el+12,7:12)=ElemJac(10:12,7:12)

! add the contribution due to follower distributed load:\diff C^T/\theta.f
!----------------------------------------------------------
IF(follower_load)THEN
	ElemJac(1:6,4:6)=ElemJac(1:6,4:6)-GetLoadJ(-1,dL,Le,eCTtheta,load)
	ElemJac(ndof_el+1:ndof_el+6,4:6)=ElemJac(ndof_el+1:ndof_el+6,4:6)-GetLoadJ(1,dL,Le,eCTtheta,load)  
ENDIF

! for dynamic analysis
!------------------------------------------
IF(ndof_el>=18) THEN
	tilde_omega=Tilde(eOmega_a)
	IF(init_flag==2) tilde_omega=tilde_omega+two_divide_dt*I3

    ! additional terms for fui due to P
	!--------------------------------------------
	! coefficient of theta of Fu due to ePi
    temp31=MATMUL(eCabhalfL,ePi)
	temp33=MATMUL(tilde_omega,MATMUL3(eCTtheta,temp31))
	ElemJac(1:3,4:6)=ElemJac(1:3,4:6)+temp33   
	ElemJac(ndof_el+1:ndof_el+3,4:6)=ElemJac(ndof_el+1:ndof_el+3,4:6)+temp33   

	! coefficient of p of Fu
	ElemJac(1:3,13:15)=MATMUL(tilde_omega,eCTCabhalfL)    
	ElemJac(ndof_el+1:ndof_el+3,13:15)=ElemJac(1:3,13:15)

    ! additional terms for fGammaD due to P,H
	!--------------------------------------------
	! coefficient of theta of fGammaD due to P,H
    temp33=MATMUL3(eCTtheta,MATMUL(eCabHalfL,CrossProduct(Vi,ePi)))+ MATMUL(tilde_omega,MATMUL3(eCTtheta,MATMUL(eCabhalfL,Hi)))
	ElemJac(4:6,4:6)=ElemJac(4:6,4:6)+ temp33  
    ElemJac(ndof_el+4:ndof_el+6,4:6)=ElemJac(ndof_el+4:ndof_el+6,4:6)+ temp33  

	! coefficient of P,H of FGammaD
	ElemJac(4:6,13:15)=MATMUL(eCTCabhalfL,Tilde(Vi)-MATMUL(Tilde(ePi),eMass(1:3,1:3)))  
	ElemJac(4:6,16:18)=MATMUL(tilde_omega,eCTCabhalfL)-MATMUL(eCTCabhalfL,MATMUL(Tilde(ePi),eMass(1:3,4:6)))
    ElemJac(ndof_el+4:ndof_el+6,13:18)=ElemJac(4:6,13:18)
	
	! the Jacobian of fPi
	!------------------------------------------
	ElemJac(NDOF_ND+1:NDOF_ND+3,1:3)=-tilde_omega
    ElemJac(NDOF_ND+1:NDOF_ND+3,4:6)=MATMUL3(eCTtheta,MATMUL(eCab,Vi))
    ElemJac(NDOF_ND+1:NDOF_ND+3,13:18)=MATMUL(eCTCab,eMass(1:3,:))
    
	! the Jacobian of fHi
	!------------------------------------------
	DO i=1,NDIM
		eCTtheta2(i,:,:)=TRANSPOSE(eCTtheta(i,:,:))
	ENDDO
	ElemJac(NDOF_ND+4:NDOF_ND+6,4:6)=-MATMUL(TRANSPOSE(eCab),MATMUL3(eCTtheta2,eOmega_a))
    ElemJac(NDOF_ND+4:NDOF_ND+6,13:18)=eMass(4:6,:)    

	IF(init_flag==2) THEN
	   IF(niter/=1)THEN  ! nonlinear solution
			tt4=DOT_PRODUCT(theta,theta)*0.25D0
    	
			DO i=1,NDIM
				temp333(i,:,:)=(Tilde(I3(i,:)*0.5d0*(1+tt4))+MATMUL(I3-Tilde(theta*0.5d0),ekttek(i,:,:)*0.25d0))/(1+tt4)**2
			ENDDO
			temp33=MATMUL3(temp333,ThetaDot)
			temp33=temp33-(I3-Tilde(theta*.5D0))/(1+tt4)*two_divide_dt
        ELSE  ! linear solution
		    temp33=-I3*two_divide_dt
		ENDIF

		ElemJac(NDOF_ND+4:NDOF_ND+6,4:6)=	ElemJac(NDOF_ND+4:NDOF_ND+6,4:6)+MATMUL(TRANSPOSE(eCab),temp33)
	ENDIF

ENDIF




! quasi-steady aerodynamics
IF(a_flag == 1) THEN
	!the Jacobian of the lift
	! coefficient theta of lift
	temp33 = OuterProduct(dir_lift,dL*PI*Rho*bw*U*(Jhdot_THETA+bw*(0.5-aw)*Jalphadot_THETA+U*Jalpha_THETA))+ &
	& dL*PI*Rho*bw*U*(hdot+bw*(0.5-aw)*alphadot+U*alpha)*Jdir_lift
	ElemJac(1:3,4:6)=ElemJac(1:3,4:6)-temp33
	ElemJac(ndof_el+1:ndof_el+3,4:6)=ElemJac(ndof_el+1:ndof_el+3,4:6)-temp33

	! coefficient PH of lift	
	temp36 = OuterProduct(dir_lift,dL*PI*Rho*bw*U*(Jhdot_PH+bw*(0.5-aw)*Jalphadot_PH))
	ElemJac(1:3,13:18)=ElemJac(1:3,13:18)-temp36
	ElemJac(ndof_el+1:ndof_el+3,13:18)=ElemJac(ndof_el+1:ndof_el+3,13:18)-temp36
	
	!the Jacobian of the moment
	! coefficient theta of moment
	temp31 = PI*Rho*U*(aw + 0.5)*bw**2*dL*(U*Jalpha_THETA + Jhdot_THETA) + &
	& (PI*Rho*U*(aw + 0.5)*(-aw + 0.5)*bw**3*dL - 0.25*PI*Rho*U*bw**3*dL)*Jalphadot_THETA
	temp33 = OuterProduct(dir_moment,temp31)+(U*alpha + hdot)*PI*Rho*U*(aw + 0.5)*bw**2*dL + (PI*Rho*U*(aw + 0.5)*(-aw &
	& + 0.5)*bw**3*dL - 0.25*PI*Rho*U*bw**3*dL)*alphadot*Jdir_moment
	ElemJac(4:6,4:6)=ElemJac(4:6,4:6)-temp33
	ElemJac(ndof_el+4:ndof_el+6,4:6)=ElemJac(ndof_el+4:ndof_el+6,4:6)-temp33
	! coefficient PH of moment
	temp36 = OuterProduct(dir_moment,PI*Rho*U*(aw + 0.5)*bw**2*dL*Jhdot_PH + (PI*Rho*U*(aw + 0.5)*(-aw &
	& + 0.5)*bw**3*dL-0.25*PI*Rho*U*bw**3*dL)*Jalphadot_PH)
	ElemJac(4:6,13:18)=ElemJac(4:6,13:18)-temp36
	ElemJac(ndof_el+4:ndof_el+6,13:18)=ElemJac(ndof_el+4:ndof_el+6,13:18)-temp36	
ENDIF

!~ ! quasi steady aerodynamic with added mass and Peters

IF(a_flag >= 2) THEN
	!the Jacobian of lift
	!-------------------------
!~ 		lift = 0.5*dL*PI*Rho*bw**2*(hdotdot-bw*aw*alphadotdot+U*alphadot) + &
!~ 		& dL*PI*Rho*bw*U*(hdot+bw*(0.5-aw)*alphadot+U*alpha-lambda0)
	Jlift = 0.D0
	! coefficient theta
	temp31= dL*PI*Rho*bw*U*(U*Jalpha_THETA)
	Jlift(:,4:6) = OuterProduct(dir_lift,temp31)
	!coefficient PH
	temp6 = 0.5*dL*PI*Rho*bw**2*(Jhdotdot_PH-bw*aw*Jalphadotdot_PH+U*Jalphadot_PH) + &
		& dL*PI*Rho*bw*U*(Jhdot_PH+bw*(0.5-aw)*Jalphadot_PH)
	Jlift(:,13:18) = OuterProduct(dir_lift,temp6)
	IF (a_flag ==3) THEN
		!coefficient Lambda
		tempN = dL*PI*Rho*bw*U*(-Jlambda0_LAMBDA)
		Jlift(:,19:ndof_el) = OuterProduct(dir_lift,tempN)
	ENDIF

	!the Jacobian of moment
	!---------------------------------
!~ 		moment = bw*(0.5+aw)*(0.5*dL*PI*Rho*bw**2*(hdotdot-bw*aw*alphadotdot+U*alphadot) + &
!~ 		& dL*PI*Rho*bw*U*(hdot+bw*(0.5-aw)*alphadot+U*alpha-lambda0)) - &
!~ 		& 0.5*dL*PI*Rho*bw**3*(0.5*hdotdot+U*alphadot+bw*(0.125-0.5*aw)*alphadotdot)	
	Jmoment = 0.D0
	! coefficent theta
	temp31 = bw*(0.5+aw)*(dL*PI*Rho*bw*U*(U*Jalpha_THETA))
	Jmoment(:,4:6) = OuterProduct(dir_moment,temp31)
	!coefficient PH
	temp6 = bw*(0.5+aw)*(0.5*dL*PI*Rho*bw**2*(Jhdotdot_PH-bw*aw*Jalphadotdot_PH+U*Jalphadot_PH) + &
		& dL*PI*Rho*bw*U*(Jhdot_PH+bw*(0.5-aw)*Jalphadot_PH)) - &
		& 0.5*dL*PI*Rho*bw**3*(0.5*Jhdotdot_PH+U*Jalphadot_PH+bw*(0.125-0.5*aw)*Jalphadotdot_PH)	
	Jmoment(:,13:18) = OuterProduct(dir_moment,temp6)
	IF (a_flag ==3) THEN
		!coefficient Lambda
		tempN = bw*(0.5+aw)*(dL*PI*Rho*bw*U*(-Jlambda0_LAMBDA))
		Jmoment(:,19:ndof_el) = OuterProduct(dir_moment,tempN)
	ENDIF
	
	! Contribution of lift and moment to the Jacobian
	! sign - come from the -fi and -mi in the expression of fui and fPhii
	ElemJac(1:3,:) = ElemJac(1:3,:) - Jlift(:,1:ndof_el)
	ElemJac(ndof_el+1:ndof_el+3,:) = ElemJac(ndof_el+1:ndof_el+3,:) - Jlift(:,1:ndof_el)
	
	ElemJac(4:6,:) = ElemJac(4:6,:) - Jmoment(:,1:ndof_el)
	ElemJac(ndof_el+4:ndof_el+6,:) = ElemJac(ndof_el+4:ndof_el+6,:) - Jmoment(:,1:ndof_el)
	

	
	IF (a_flag == 3) THEN
		! The Jacobian of flambdai
		! coefficient PH
		temp6 = -bw*(Jhdotdot_PH+U*Jalphadot_PH+bw*(0.5-aw)*Jalphadotdot_PH)
		tempN6 = OuterProduct(c,temp6)
		ElemJac(ndof_el-NSTATES+1:ndof_el,13:18)=ElemJac(ndof_el-NSTATES+1:ndof_el,13:18)+tempN6	
		! coefficient lambda
!~ 		tempNN = U*Ident
        IF (U<TOLERANCE) THEN ! avoiding the zero speed singularity and force the decoupling between structural and induced flow modes at very low speed
		    tempNN = 1e12*Ident
        ELSE
            tempNN = max(1.D0/max(1e-12,U),U)*Ident
        ENDIF    



		IF(init_flag==2) 	tempNN = tempNN + two_divide_dt*bw*A
		ElemJac(ndof_el-NSTATES+1:ndof_el,ndof_el-NSTATES+1:ndof_el)=ElemJac(ndof_el-NSTATES+1:ndof_el,ndof_el-NSTATES+1:ndof_el)+tempNN	
	ENDIF
	
ENDIF

IF(init_flag==1) THEN !needed for the initial step calculating the initial conditions
	ElemJac(:,1:6)=0.D0
	ElemJac(1:3,1:3)=dL*.5d0*I3
    ElemJac(4:6,4:6)=ElemJac(1:3,1:3)
    ElemJac(ndof_el+1:ndof_el+3,1:3)=ElemJac(1:3,1:3)
    ElemJac(ndof_el+4:ndof_el+6,4:6)=ElemJac(1:3,1:3)
ENDIF

END SUBROUTINE ElemJacobian
!***************************************************



!************************************************************
!*                                                          *
!> Caculate the mass matrix for each element         
!*                                                          *
!*															*
!************************************************************ 
SUBROUTINE ElemMass(ndof_el,elemM)

IMPLICIT NONE

INTEGER,INTENT(IN)  :: ndof_el !< #ioaero::ndof_el

REAL(DBL),INTENT(OUT)::ElemM(:,:) !< the mass matrix for each element
REAL(DBL) tmp31(NDIM),tmp36(NDIM,NDIM+NDIM),tmp6(NDIM+NDIM)
REAL(DBL) tmpN6(NSTATES,6),tmpN3(NSTATES,3),tmpNN(NSTATES,NSTATES)

ElemM=0.0D0

! d_fu/d_theta_t
!--------------------
ElemM(1:3,4:6)=CT_THETA_T(theta,eCT,MATMUL(eCabhalfL,ePi))
ElemM(ndof_el+1:ndof_el+3,4:6)=ElemM(1:3,4:6)

! d_fu/d_P_t
!--------------------		
ElemM(1:3,13:15)=eCTCabhalfL  
ElemM(ndof_el+1:ndof_el+3,13:15)=eCTCabhalfL

! d_fGammaD/d_theta_t
!--------------------
ElemM(4:6,4:6)=CT_THETA_T(theta,eCT,MATMUL(eCabhalfL,Hi)) 
ElemM(ndof_el+4:ndof_el+6,4:6)=ElemM(4:6,4:6)

! d_fGammaD/d_H_t
!--------------------
ElemM(4:6,16:18)=eCTCabhalfL 
ElemM(ndof_el+4:ndof_el+6,16:18)=eCTCabhalfL 
	
! d_fPi/d_u_t
!------------------------------------------
ElemM(13:15,1:3)=-I3

    
! d_fHi/d_theta_t
!------------------------------------------
ElemM(16:18,4:6)=MATMUL(TRANSPOSE(eCab),(Tilde(theta*.5D0)-I3)/(1+DOT_PRODUCT(theta,theta)*0.25D0))

!~ IF(a_flag == 4) THEN
!~ 	! case a_flag = 2 : added mass only
!~ 	IF (a_flag >= 2) Jlambda0_PHDOT = 0.D0

!~ 	! d_fu/d_PHdot
!~ 	tmp6 = 0.5*dL*PI*Rho*bw**2*(Jhdotdot_PHDOT-bw*aw*Jalphadotdot_PHDOT)+&
!~ 	& dL*PI*Rho*U*bw*0*(Jhdotdot_PHDOT+bw*(0.5-aw)*Jalphadotdot_PHDOT-Jlambda0_PHDOT)
!~  	tmp36 = OuterProduct(dir_lift,tmp6)
!~ 	ElemM(1:3,13:18) = ElemM(1:3,13:18) - tmp36
!~ 	ElemM(ndof_el+1:ndof_el+3,13:18) = ElemM(1:3,13:18)
	
!~ 	! d_ftheta/d_PHdot	
!~ 	tmp6 = bw*(0.5+aw)*tmp6 - 0.5*dL*PI*Rho*bw**3*(0.5*Jhdotdot_PHDOT+bw*(0.125-0.5*aw)*Jalphadotdot_PHDOT)
!~ 	tmp36 = OuterProduct(dir_moment,tmp6)
!~ 	ElemM(4:6,13:18) = ElemM(4:6,13:18) - tmp36	
!~ 	ElemM(ndof_el+4:ndof_el+6,13:18) = ElemM(4:6,13:18)
	
!~ 	IF (a_flag == 3) THEN
!~ 		!d_flambda/d_THETAdot
!~ 		tmp31 = -bw*(U*dir_moment)
!~ 		tmpN3 = OuterProduct(c,tmp31)
!~ 		ElemM(ndof_el-NSTATES+1:ndof_el,4:6) = ElemM(ndof_el-NSTATES+1:ndof_el,4:6) + tmpN3
		
!~ 		!d_flambda/d_PHdot
!~ 		tmp6 = -bw*(Jhdotdot_PHDOT+bw*(0.5-aw)*Jalphadotdot_PHDOT)
!~ 		tmpN6 = OuterProduct(c,tmp6)
!~ 		ElemM(ndof_el-NSTATES+1:ndof_el,13:18) = ElemM(ndof_el-NSTATES+1:ndof_el,13:18) + tmpN6

!~ 		!d_flambda/d_lambdadot
!~ 		tmpNN = bw*A
!~ 		ElemM(ndof_el-NSTATES+1:ndof_el,ndof_el-NSTATES+1:ndof_el) = ElemM(ndof_el-NSTATES+1:ndof_el,ndof_el-NSTATES+1:ndof_el) + tmpNN
!~ 	ENDIF

!~ ENDIF

IF(a_flag >= 2) THEN
	!the Jacobian of lift
	!-------------------------
!~ 		lift = 0.5*dL*PI*Rho*bw**2*(hdotdot-bw*aw*alphadotdot+U*alphadot) + &
!~ 		& dL*PI*Rho*bw*U*(hdot+bw*(0.5-aw)*alphadot+U*alpha-lambda0)
	Jlift = 0.D0
	!coefficient PHDOT
	tmp6 = 0.5*dL*PI*Rho*bw**2*(Jhdotdot_PHDOT-bw*aw*Jalphadotdot_PHDOT)
	Jlift(:,13:18) = OuterProduct(dir_lift,tmp6)

	!the Jacobian of moment
	!---------------------------------
!~ 		moment = bw*(0.5+aw)*(0.5*dL*PI*Rho*bw**2*(hdotdot-bw*aw*alphadotdot+U*alphadot) + &
!~ 		& dL*PI*Rho*bw*U*(hdot+bw*(0.5-aw)*alphadot+U*alpha-lambda0)) - &
!~ 		& 0.5*dL*PI*Rho*bw**3*(0.5*hdotdot+U*alphadot+bw*(0.125-0.5*aw)*alphadotdot)	
	Jmoment = 0.D0
	! coefficient PHDOT
	tmp6= bw*(0.5+aw)*(0.5*dL*PI*Rho*bw**2*(Jhdotdot_PHDOT-bw*aw*Jalphadotdot_PHDOT)) - &
	& 0.5*dL*PI*Rho*bw**3*(0.5*Jhdotdot_PHDOT+bw*(0.125-0.5*aw)*Jalphadotdot_PHDOT)
	Jmoment(:,13:18) = OuterProduct(dir_moment,tmp6)
		
	! Contribution of lift and moment to the Jacobian
	! sign - come from the -fi and -mi in the expression of fui and fPhii
	ElemM(1:3,:) = ElemM(1:3,:) - Jlift(:,1:ndof_el)
	ElemM(ndof_el+1:ndof_el+3,:) = ElemM(ndof_el+1:ndof_el+3,:) - Jlift(:,1:ndof_el)
	
	ElemM(4:6,:) = ElemM(4:6,:) - Jmoment(:,1:ndof_el)
	ElemM(ndof_el+4:ndof_el+6,:) = ElemM(ndof_el+4:ndof_el+6,:) - Jmoment(:,1:ndof_el)

	IF (a_flag ==3) THEN
		! The Jacobian of flambdai
!~ 			flambdai = bw*MATMUL(A,lambdadot)+U*lambda-bw*(hdotdot+U*alphadot+bw*(0.5-aw)*alphadotdot)*c
		! coefficient PHDOT
		tmp6 = -bw*(Jhdotdot_PHDOT+bw*(0.5-aw)*Jalphadotdot_PHDOT)
		tmpN6 = OuterProduct(c,tmp6)
		ElemM(ndof_el-NSTATES+1:ndof_el,13:18) = ElemM(ndof_el-NSTATES+1:ndof_el,13:18) + tmpN6
		! coefficient LAMBDADOT
		tmpNN = bw*A
		ElemM(ndof_el-NSTATES+1:ndof_el,ndof_el-NSTATES+1:ndof_el) = ElemM(ndof_el-NSTATES+1:ndof_el,ndof_el-NSTATES+1:ndof_el) + tmpNN
	ENDIF
		
ENDIF


END SUBROUTINE ElemMass
!***************************************************



!************************************************************
!*                                                          *                                      
!>  Extract element properties needed for element assembly 
!************************************************************
SUBROUTINE  ExtractElementProperties(elem_no,memb_info_i,x_elem,v_root_a,omega_a,ndof_el,init_elem,aero_flag,grav_flag)  

INTEGER,INTENT(IN)::elem_no !<Element index
TYPE (MemberInf),INTENT(IN)::memb_info_i !< the paramater of the members (see memberinf Type)

REAL(DBL),INTENT(IN)::x_elem(:) !< the 12 variables of the element u_i, \theta_i, F_i, M_i, for dynamic analysis, 6 more Pi, Hi. for initial step, x_elem contains, CTCabPdot, CTCabHdot, F_i, M_i, P_i, H_i. Finally for Peters aero the Ns induced flow states
REAL(DBL),INTENT(IN)::v_root_a(:) !< linear velocity of the root point in inertial frame
REAL(DBL),INTENT(IN)::omega_a(:) !< angular velocity of the root point in inertial frame

REAL(DBL),INTENT(INOUT)::init_elem(:,:) !< initial step: the 12 initial values of the element u_i,\f$ \theta_i, \dot{u}_i, \dot{theta}_i.\f$ Time marching: \f$2/dt ui+\dot{u}_i, 2/dt thetai+\dot{theta}_i, 2/dt CTCabP+\dot, 2/dt CTCabP+dot \f$
INTEGER,INTENT(IN)::ndof_el !<#ioaero::nedof_el
INTEGER,INTENT(IN)::aero_flag !< #ioaero::aero_flag
INTEGER,INTENT(IN)::grav_flag !< #ioaero::grav_flag

! coordinate system parameter
REAL(DBL)            :: xflow(NDIM),xB(NDIM),yB(NDIM),zB(NDIM)
! additionnal variables used for the aerodynamic model
REAL(DBL)            :: PHDotB(2*NDIM),ViDot(NDIM),OMEGAiDot(NDIM),init_PH(2*NDIM)

!----------------------------------------------------
! Extract sectional properties for this element
!----------------------------------------------------
eFlex=memb_info_i%mate(elem_no,1:NSTRN,:) ! element flexibility matrix
IF(ndof_el>=18) eMass=memb_info_i%mate(elem_no,NSTRN+1:,:)! element inverse mass matrix for all the analyses which are not static

!------------------------------------------------------------------

!----------------------------------------------------------------------------
! Extract b frame for this element
!----------------------------------------------------------------------------
eCab=memb_info_i%triad(elem_no,:,:)
dL=memb_info_i%dL
Le=memb_info_i%Le(elem_no) ! the ending x_1 of the element

!Angular velocity and linear velocity
!---------------------------------------
IF(ndof_el>=18)	THEN
    eOmega_a=omega_a  ! angular velocity

	!linear velocity at the middle of the element
	!-----------------------------------------------
	ev_i=v_root_a+CrossProduct(eOmega_a,memb_info_i%coordinate(elem_no,:)-xyz_pt1)  
ENDIF

!----------------------------------------------------------------------------
! Extract internal variables for this element
!----------------------------------------------------------------------------
IF(init_flag==1) THEN  ! initial step
    Ui       = init_elem(elem_no,1:3)          ! given initial displacements
	theta    = init_elem(elem_no,4:6)          ! given initial rotations
	UiDot    = init_elem(elem_no,7:9)          ! given initial linear velocities
	ThetaDot = init_elem(elem_no,10:12)        ! given initial angular velocities
	CTCabPdot= x_elem(1:3)             ! calculated CTCabPdot from previous step
    CTCabHdot= x_elem(4:6)             ! calculated CTCabHdot from previous step
    IF (a_flag==3) THEN
		lambda = 0.D0
		lambdadot = 0.D0
	ENDIF
ELSE 
	Ui   =x_elem(1:3)
	theta=x_elem(4:6)
	IF(init_flag==2) THEN  ! time marching
		UiDot=-init_elem(elem_no,1:3)
		ThetaDot=two_divide_dt*theta-init_elem(elem_no,4:6)
		CTCabPdot=-init_elem(elem_no,7:9)
		CTCabHdot=-init_elem(elem_no,10:12)
	ENDIF
ENDIF

Fi   =x_elem(7:9)
Mi   =x_elem(10:12)

IF (ndof_el>=18) THEN
	ePi   =x_elem(13:15)
	Hi    =x_elem(16:18)	
ENDIF

e1GammaD=MATMUL(eFlex(1:3,:),x_elem(7:12))+e1
kappa=MATMUL(eFlex(4:6,:),x_elem(7:12))

IF(ndof_el>=18)THEN
	Vi=MATMUL(eMass(1:3,:),x_elem(13:18))
	OMEGAi=MATMUL(eMass(4:6,:),x_elem(13:18))
ENDIF

! Some numbers frequently needed for the element
!-----------------------------------------
eCabhalfL   =eCab*dL*0.5d0
eCT=DirCosineTRodrigues(theta)
eCTCab       =MATMUL(eCT,eCab)
eCTCabhalfL =eCTCab*dL*0.5d0

a_flag = aero_flag
g_flag = grav_flag

! Numbers related to aerodynamics (explanation in "modele aeroelastique couple")
IF (a_flag /= 0) THEN
	! for all model
	U =  memb_info_i%aerodyn_coef(elem_no,1)
	Rho = memb_info_i%aerodyn_coef(elem_no,2)
	Chord = memb_info_i%aerodyn_coef(elem_no,3)
	aw = memb_info_i%aerodyn_coef(elem_no,4)
	bw = 0.5*Chord
	
	! angles between flow and aircraft coordinate system 

	alpha_ac = memb_info_i%aerodyn_coef(elem_no,5)
	beta_ac = memb_info_i%aerodyn_coef(elem_no,6)
	
!~ 	! transfert matrix from flow to aircraft coordinates system eCaf
!~ 	eCaf(:,1) = (/COS(alpha_ac)*COS(beta_ac),-COS(alpha_ac)*SIN(beta_ac),-SIN(alpha_ac)/)
!~ 	eCaf(:,2) = (/SIN(beta_ac),COS(beta_ac),0.D0/)
!~ 	eCaf(:,3) = (/SIN(alpha_ac)*COS(beta_ac),-SIN(alpha_ac)*SIN(beta_ac),COS(alpha_ac)/)
!~ 	eCaf = TRANSPOSE(eCaf)
    Wind = (/-COS(alpha_ac)*COS(beta_ac),-SIN(beta_ac),-SIN(alpha_ac)*COS(beta_ac)/)
    xflow = -Wind
	! extraction of aerodynamic section parameter from the coordinate system
	! all vector are defined in frame a
	!--------------------------------------------------------------------------------------------------------
    xB = eCTCAB(:,1) 
    yB = eCTCAB(:,2) 
    zB = eCTCAB(:,3) 

    ! direction of the airfoil coordinate system (xflow,dir_moment,dir_lift) in the aircraft reference
    dir_lift = CrossProduct(xB,xflow)
    dir_lift = 1/Norm(dir_lift)*dir_lift
    dir_moment = CrossProduct(xflow,dir_lift)
    dir_moment = 1/Norm(dir_moment)*dir_moment

    ! airfoil parameter
    ! alpha ~ alpha_deformation + alpha_aircraft
    alpha = ASIN(DOT_PRODUCT(yB,dir_lift))

    alphadot = DOT_PRODUCT(MATMUL(eCTCAB,OMEGAi),dir_moment)
!~     alphadot = OMEGAi(1)
    hdot = -DOT_PRODUCT(MATMUL(eCTCAB,Vi),dir_lift)

!~     beta = ASIN(DOT_PRODUCT(xB,xflow))

	! flow velocity correction with beta :
!~ 	U = U*COS(beta)
    ! coefficient theta of Jacobian
    !---------------------------------------
    Jalpha_THETA = dir_moment
    Jhdot_THETA = 0.D0

	! coefficient PH of Jacobian
	Jalphadot_PH = MATMUL(dir_moment,MATMUL(eCTCab,eMass(4:6,:)))
!~ 	Jalphadot_PH = eMass(4,:)
	Jalphadotdot_PHDOT = Jalphadot_PH
	Jhdot_PH = -MATMUL(dir_lift,MATMUL(eCTCab,eMass(1:3,:)))
	Jhdotdot_PHDOT = Jhdot_PH

	! determination of alphadotdot and hdotdot for dynamic simulation
	IF (init_flag==2) THEN

		
		! transfer PHDot at t=t0 from frame a to B
		init_PH(1:3) = MATMUL(TRANSPOSE(eCTCAB),init_elem(elem_no,7:9))
		init_PH(4:6) = MATMUL(TRANSPOSE(eCTCAB),init_elem(elem_no,10:12))		
		! evaluate PHDot at t=t0+dt
		PHDotB = two_divide_dt*x_elem(13:18)-init_PH		
		! evaluate linear and angular acceleration in frame a using PHDot
		ViDot = MATMUL(eMass(1:3,:),PHDotB)
		OMEGAiDot = MATMUL(eMass(4:6,:),PHDotB)
		! alphadotdot = angular acceleration about dir_moment; hdotdot = -linear acceleration about dir_lift
		alphadotdot = DOT_PRODUCT(OMEGAiDot,MATMUL(TRANSPOSE(eCTCAB),dir_moment))
		hdotdot = -DOT_PRODUCT(ViDot,MATMUL(TRANSPOSE(eCTCAB),dir_lift))			
		Jalphadotdot_PH = two_divide_dt*MATMUL(eMass(:,4:6),MATMUL(TRANSPOSE(eCTCAB),dir_moment))
		Jhdotdot_PH = -two_divide_dt*MATMUL(eMass(:,1:3),MATMUL(TRANSPOSE(eCTCAB),dir_lift))
		
	ELSE ! for modal analysis alphadotdot and hdotdot are described in the mass matrix
		alphadotdot = 0.D0
		hdotdot = 0.D0
		Jalphadotdot_PH = 0.D0
		Jhdotdot_PH = 0.D0
	ENDIF

ENDIF

!Induced-flow states
IF (a_flag > 2) THEN
	
	P = Peters(NSTATES)
	A = P(:,1:NSTATES)
	b = P(:,NSTATES+1)
	c = P(:,NSTATES+2)
	Ident = P(:,NSTATES+3:)

	lambda = x_elem(ndof_el-NSTATES+1:ndof_el)
	IF (init_flag == 2) THEN
		lambdadot = two_divide_dt*lambda-init_elem(elem_no,12+1:12+NSTATES)
	ELSE
		lambdadot = 0.D0
	ENDIF
	lambda0 = 0.5*DOT_PRODUCT(b,lambda)
	Jlambda0_LAMBDA = 0.5*b

ENDIF

! Numbers related to gravity
IF (g_flag > 0) THEN
	X_cg = memb_info_i%aerodyn_coef(elem_no,8)
ENDIF

END SUBROUTINE ExtractElementProperties
!**********************************************************
   
END MODULE Element
