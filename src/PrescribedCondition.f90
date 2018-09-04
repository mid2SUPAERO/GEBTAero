!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************
!=============================================================
!
! A module for defining prescribed conditions including both
! concentrated information and distributed information 
!
!=============================================================
MODULE PrescribedCondition

USE GlobalDataFun
IMPLICIT NONE

PRIVATE       ! So everything is private, except declared by PUBLIC
PUBLIC DistriLoad,PrescriInf

PUBLIC ExistPI,InitPI,InputEchoPrescribedConditions,GetPrescribedDOF,GetPrescribedVal,UpdatePI,UpdateFollower,&
   & GetDistributedLoad,GetLoad,GetLoadJ,FollowerJ,InitPIAero

! Define the prescribed condition
!-----------------------------------
TYPE PrescriInf
     PRIVATE
	 INTEGER   ::id                 ! where it is applied, could be a node number or member number
	 INTEGER   ::dof(NSTRN)         ! maximum 6 degrees of freedom can be prescribed, 	                         
	                                ! for distributed loads, it is used to denote the distribution function no
     REAL(DBL) ::value(NSTRN)       ! the magnitude of the prescribed values
	 INTEGER   ::time_fun_no(NSTRN) ! which time function is used
	 INTEGER   ::follower(NSTRN)    ! whether the prescribed quantity is a follower or not: 1 is a follower; 0 is not
	 REAL(DBL) ::value_current(NSTRN) ! indicate the current functional value updated by time steps, calculated internally
END TYPE PrescriInf

! Define the distributed load condition 
!-------------------------------------------
TYPE DistriLoad
     PRIVATE
     REAL(DBL) ::value(NSTRN)           ! the current functional value of the load
	 REAL(DBL) ::distr_fun(NSTRN,NSTRN) ! the distribution function for each load
	 INTEGER   ::follower(NSTRN)        ! whether the force vector/moment vector are follower quantities
END TYPE DistriLoad

CONTAINS
!========================================================


!*****************************************************************
!*                                                               *
!*  Determine whether prescribed condition exist, and whether    *
!*  any of such conditions is a follower condition               *
!*                                                               *
!*****************************************************************
SUBROUTINE ExistPI(location,prescri_inf,exist_pi,follower_pi)
INTEGER,INTENT(IN)   ::location
TYPE(PrescriInf),INTENT(IN)::prescri_inf(:) 
LOGICAL,INTENT(OUT)::exist_pi,follower_pi
INTEGER::i

DO i=1,SIZE(prescri_inf)
   IF(location==prescri_inf(i)%id) THEN
     exist_pi=.TRUE.
     IF(ANY(prescri_inf(i)%follower==1))follower_pi=.TRUE.
	 RETURN
   ENDIF
ENDDO

END SUBROUTINE ExistPI
!*****************************************************************



!************************************************************
!*                                                          *
!*  Obtain the distributed load condition                   *
!*                                       					*
!************************************************************
FUNCTION GetDistributedLoad(memb_no,mb_condition,distr_fun) RESULT(res)

INTEGER,INTENT(IN)   ::memb_no
REAL(DBL),INTENT(IN) ::distr_fun(:,:)
TYPE(PrescriInf),INTENT(IN)::mb_condition(:) ! distributed load information

TYPE(DistriLoad)::res

REAL(DBL)::tmpR(NSTRN,NSTRN) ! hold the distribution functions for the six dofs
INTEGER:: i, j,distr_fun_no

DO i=1,SIZE(mb_condition)
   IF(memb_no==mb_condition(i)%id) THEN
	 tmpR=0.d0

     DO j=1,NSTRN
	    distr_fun_no=mb_condition(i)%DOF(j)
		IF(distr_fun_no/=0) tmpR(j,:)=distr_fun(distr_fun_no,:) !distribution function for this dof
     ENDDO
     res=DistriLoad(mb_condition(i)%value_current,tmpR,mb_condition(i)%follower) 
   ENDIF

ENDDO

END FUNCTION GetDistributedLoad
!************************************************************




!********************************************************
! Obtain the distributed load, transform if follower
!********************************************************
FUNCTION GetLoad(flag,dL,Le,eCT,load,follower_load) RESULT(res)

INTEGER,INTENT(IN)::flag ! if flag=-1, the starting portion, if flag=1, the ending portion
REAL(DBL),INTENT(IN)::dL,Le,eCT(:,:)
TYPE (DistriLoad),INTENT(IN)::load
LOGICAL,INTENT(IN)::follower_load

REAL(DBL)::res(NSTRN)
INTEGER:: i

res=load%value

! Obtain the load value
!---------------------------
DO i=1,NSTRN
	IF(ABS(res(i))>TOLERANCE) res(i)=res(i)*LoadIntegration(flag,dL,Le,load%distr_fun(i,:)) 
ENDDO

IF(follower_load) THEN
    res(1:3)=TransferFollower(load%follower(1:3),res(1:3),eCT)
    res(4:6)=TransferFollower(load%follower(4:6),res(4:6),eCT)
ENDIF

END FUNCTION GetLoad
!****************************************



!************************************************************
!*                                                          *
!*    Obtain the jacobian due to follower distributed load  *
!*															*
!************************************************************
FUNCTION GetLoadJ(flag,dL,Le,eCTtheta,load) RESULT(res)

INTEGER,INTENT(IN)::flag ! if flag=-1, the starting portion, if flag=1, the ending portion
REAL(DBL),INTENT(IN)::dL,Le,eCTtheta(:,:,:)
TYPE (DistriLoad),INTENT(IN)::load

REAL(DBL)::res(NSTRN,3)
REAL(DBL)::val(NSTRN)
INTEGER:: i

! Obtain the load value
!---------------------------
val=load%value

DO i=1,NSTRN
	IF(ABS(val(i))>TOLERANCE) val(i)=val(i)*LoadIntegration(flag,dL,Le,load%distr_fun(i,:)) 
ENDDO

res=0.0D0

IF(ANY(load%follower(1:3)==1)) res(1:3,:)=FollowerJ(load%follower(1:3),val(1:3),eCTtheta)
IF(ANY(load%follower(4:6)==1)) res(4:6,:)=FollowerJ(load%follower(4:6),val(4:6),eCTtheta)

END FUNCTION GetLoadJ
!************************************************************




!*****************************************************************
!*                                                               *
!*  Obtain Prescribed dof and follower condition                 *
!*															     *
!*****************************************************************
SUBROUTINE GetPrescribedDOF(nkp,pt_condition,kp_dof,kp_follower)

INTEGER,INTENT(IN)::nkp
TYPE(PrescriInf),INTENT(IN)::pt_condition(:) ! point condition

INTEGER,INTENT(OUT)::kp_dof(:,:)
INTEGER,INTENT(OUT)::kp_follower(:,:)

INTEGER:: i,j

DO j=1,nkp
	DO i=1,SIZE(pt_condition)
		IF(j==pt_condition(i)%id) THEN
			kp_dof(j,:)=pt_condition(i)%dof
			kp_follower(j,:)=pt_condition(i)%follower
		ENDIF
	ENDDO
ENDDO

END SUBROUTINE GetPrescribedDOF
!************************************************************



!*****************************************************************
!*                                                               *
!*  Obtain Prescribed value                                      *
!*															     *
!*****************************************************************
SUBROUTINE GetPrescribedVal(nkp,pt_condition,kp_cond)

INTEGER,INTENT(IN)::nkp
TYPE(PrescriInf),INTENT(IN)::pt_condition(:) ! point condition

REAL(DBL),INTENT(OUT)::kp_cond(:,:)

INTEGER:: i,j

DO j=1,nkp
	DO i=1,SIZE(pt_condition)
		IF(j==pt_condition(i)%id) kp_cond(j,:)=pt_condition(i)%value_current
	ENDDO
ENDDO

END SUBROUTINE GetPrescribedVal
!************************************************************



!************************************************************
!*                                                          *
!*  Initialize Prescribed Conditions                        *
!*															*
!************************************************************
ELEMENTAL FUNCTION InitPI() RESULT(res)

TYPE (PrescriInf)::res       

res%id=0
res%dof=0
res%value=0._DBL
res%time_fun_no=0
res%follower=0
res%value_current=0._DBL

END FUNCTION InitPI
!************************************************************



!************************************************************
!*                                                          *
!*  Input and echo Prescribed Conditions                    *
!*															*
!************************************************************
FUNCTION InputEchoPrescribedConditions(IN,EIN,error) RESULT(res)

INTEGER,INTENT(IN)::IN,EIN ! file units for input and echo files, respectively
CHARACTER(*),INTENT(OUT)::error
TYPE (PrescriInf)::res       
   
READ(IN,*,IOSTAT=in_stat) res%id   
IF(IOError('read location for prescribed condition',error)) GOTO 9999
    
READ(IN,*,IOSTAT=in_stat)res%dof      
IF(IOError('read dof or distribution function no for prescribed condition',error)) GOTO 9999

READ(IN,*,IOSTAT=in_stat)res%value 
IF(IOError('read values for each dof',error)) GOTO 9999

READ(IN,*,IOSTAT=in_stat)res%time_fun_no 
IF(IOError('read time function #',error)) GOTO 9999

READ(IN,*,IOSTAT=in_stat)res%follower 
IF(IOError('read follower conditions',error)) GOTO 9999

res%value_current=res%value ! the current value is initialized to be the original functional value, if not time varing
                            ! it will remain the same. If time varying, it will be updated.
WRITE(EIN,*) 
WRITE(EIN,*) 'Location=', res%id
WRITE(EIN,*) '----------------------------------'	
CALL WriteVec(EIN, res%dof)	
CALL WriteVec(EIN, res%value)
CALL WriteVec(EIN, res%time_fun_no)	
CALL WriteVec(EIN, res%follower)	

9999 IF (error/="") THEN
    WRITE(0,*) "PrescribedCondition.f90 : ",error
    ERROR STOP 110
ENDIF
	
END FUNCTION InputEchoPrescribedConditions
!************************************************************


!************************************************************
!*                                                          *
!*  Obtain Prescribed DOF and value needed for rhs          *
!*  assume only follower force/moments, and no displacements*
!*  or rotations can be prescribed for follower quantities. *
!*  And the first three prescribed dofs for the point with a* 
!*  follower component should be either 7 8 9 or 10 11 12.  *
!*  This assumption is made for the easiness to locate the  *
!*  rotation parameters                                     *
!*															*
!************************************************************
FUNCTION UpdateFollower(kp_dof,kp_follower,kp_cond,x_pt) RESULT(res)

INTEGER,INTENT(IN)::kp_dof(:),kp_follower(:)
REAL(DBL),INTENT(IN)::kp_cond(:),x_pt(:)
REAL(DBL)::res(NSTRN)

REAL(DBL)::CT(3,3) ! direction cosine matrix

res=kp_cond

IF(kp_dof(1)==7) THEN
	 CT=DirCosineTRodrigues(x_pt(4:6))  
ELSE 
     CT=DirCosineTRodrigues(x_pt(1:3))
ENDIF

IF(ANY(kp_follower(1:3)==1)) &
	& res(1:3)=TransferFollower(kp_follower(1:3),kp_cond(1:3),CT)
IF(ANY(kp_follower(4:6)==1)) &
	& res(4:6)=TransferFollower(kp_follower(4:6),kp_cond(4:6),CT)


END FUNCTION UpdateFollower
!************************************************************


!************************************************************
!*                                                          *
!*  Update the prescribed information based on the current  *
!*  time the value is stored in: value_current              *
!*															*
!************************************************************
SUBROUTINE UpdatePI(prescri_inf,time_fun,t) 

USE TimeFunctionModule

TYPE (PrescriInf),INTENT(INOUT)::prescri_inf(:) 
TYPE (TimeFunction),INTENT(IN)::time_fun(:)
REAL(DBL),INTENT(IN)::t  ! current time

INTEGER:: i,j,time_no
REAL(DBL)::val(NSTRN)

DO j=1,SIZE(prescri_inf)
	val=prescri_inf(j)%value
	DO i=1,NSTRN
	   time_no=prescri_inf(j)%time_fun_no(i) ! obtain time function #
	   IF(time_no/=0) THEN
	      IF(val(i)/=0.0D0)  prescri_inf(j)%value_current(i)=val(i)*GetTimeFunction(time_fun(time_no),t)
       ELSE
		  IF(val(i)/=0.0D0)  prescri_inf(j)%value_current(i)=val(i)
	   ENDIF
	ENDDO     

ENDDO

END SUBROUTINE UpdatePI
!************************************************************

!==================================================
!
!Subroutines or functions used internally
!
!==================================================




!************************************************************
!*                                                          *
!* Calculating the Jacobian due to follower conditions      *
!* J=\diff{C^T.vec}/\diff\theta, return a 3x3 matrix        *
!* with ith column corresponding to the derivative          *
!* withe respect to \theta_i                                *
!*                                                          *
!************************************************************
FUNCTION FollowerJ(follower,vec,eCTtheta) RESULT(res)

INTEGER,INTENT(IN)::follower(:)
REAL(DBL),INTENT(IN)::vec(:)
REAL(DBL),INTENT(IN)::eCTtheta(:,:,:)
REAL(DBL):: res(SIZE(vec),SIZE(vec))

INTEGER:: i,j

res=0.d0

DO i=1,3 ! loop through \theta_i
		DO j=1,3 ! loop through the components of the vector
			IF(follower(j)==1.AND.ABS(vec(j))>TOLERANCE) res(1:3,i)=res(1:3,i)+eCTtheta(i,:,j)*vec(j)
		ENDDO
ENDDO

END FUNCTION FollowerJ
!********************************************************



!************************************************************
!*                                                          *
!*  Caculate the load using Chebychev polynomials           *
!*                                                          *
!*															*
!************************************************************ 
FUNCTION LoadIntegration(flag,dL,Le,func) 

INTEGER,  INTENT(IN):: flag ! if flag=-1, it is for the starting point, if flag=1, it is for the ending point
REAL(DBL),INTENT(IN)::dL ! length of the element
REAL(DBL),INTENT(IN)::Le ! the ending point of the element at the arc length of the member
REAL(DBL),INTENT(IN)::func(NSTRN) ! distributed load function for this element

REAL(DBL)::LoadIntegration ! the load for each element

REAL(DBL)::dL2, dL3, dL4, dL5,dL6 ! the powers of dL
REAL(DBL)::Le2, Le3, Le4, Le5,Le6 ! the powers of Le

REAL(DBL)::fun_coef(6)  ! coefficient of Chebychev polynomials after integration

dL2=dL*dL; dL3=dL*dL2;dL4=dL2*dL2;dL5=dL2*dL3; dL6=dL3*dL3
Le2=Le*Le; Le3=Le*Le2;Le4=Le2*Le2;Le5=Le2*Le3; Le6=Le3*Le3


SELECT CASE (flag)     
	   CASE (-1)      ! starting point
            fun_coef=(/ 0.5d0*dL,&
			        &  -0.333333333333333d0*dL2 + 0.5d0*dL*Le, &
					&  -0.5d0*dL + 0.5d0*dL3 - 1.33333333333333d0*dL2*Le + dL*Le2,&
					&   dL2 - 0.8d0*dL4 - 1.5d0*dL*Le + 3*dL3*Le - 4*dL2*Le2 + 2*dL*Le3, &
                    &   0.5d0*dL - 2*dL3 + 1.33333333333333d0*dL5 + 5.33333333333333d0*dL2*Le &
					&   - 6.4d0*dL4*Le - 4.0d0*dL*Le2 + 12.0d0*dL3*Le2 - 10.6666666666667d0*dL2*Le3 + 4.0d0*dL*Le4, &
                    &  -1.66666666666667d0*dL2 + 4.0d0*dL4 - 2.28571428571429d0*dL6 + 2.5d0*dL*Le &
					&  - 15.0d0*dL3*Le + 13.3333333333333d0*dL5*Le + 20.0d0*dL2*Le2 - 32.0d0*dL4*Le2   &
					&  - 10.0d0*dL*Le3 + 40.0d0*dL3*Le3 - 26.6666666666667d0*dL2*Le4 + 8.0d0*dL*Le5 &
					&   /)



	   CASE(1)  ! ending point
            fun_coef=(/ 0.5d0*dL,&
			        &  -0.00476190476190476d0*dL*(35.0d0*dL - 105.0d0*Le),&
					&  -0.00476190476190476d0*dL*(105.0d0 - 35.0d0*dL2 + 140.0d0*dL*Le - 210.0d0*Le2),&
					&  -0.00476190476190476d0*dL*(-105.0d0*dL + 42.0d0*dL3 + 315.0d0*Le - 210.0d0*dL2*Le + 420.0d0*dL*Le2 - 420.0d0*Le3), &
                    &  -0.00476190476190476d0*dL*(-105.0d0 + 140.0d0*dL2 - 56.0d0*dL4 - 560.0d0*dL*Le + 336.0d0*dL3*Le &
					&    + 840.0d0*Le2 - 840.0d0*dL2*Le2 + 1120.0d0*dL*Le3 - 840.0d0*Le4),&
                    &  -0.00476190476190476d0*dL*(175.0d0*dL - 210.0d0*dL3 + 80.0d0*dL5 - 525.0d0*Le + 1050.0d0*dL2*Le &
					&   - 560.0d0*dL4*Le - 2100.0d0*dL*Le2 + 1680.0d0*dL3*Le2 + 2100.0d0*Le3 - 2800.0d0*dL2*Le3 &
					&    + 2800.0d0*dL*Le4 - 1680.0d0*Le5) &
					&   /)
END SELECT

LoadIntegration=DOT_PRODUCT(fun_coef,func)


END FUNCTION LoadIntegration
!************************************************************


!************************************************************
!*                                                          *
!* Transfer follower according to C^T.vec                   *
!* note vec is a 3x1 vector, and whether a component        *
!* is a follower or not is determined by follower           *
!************************************************************
FUNCTION TransferFollower(follower,vec,CT) RESULT(res)

INTEGER,INTENT(IN)::follower(:)
REAL(DBL),INTENT(IN)::vec(:)
REAL(DBL),INTENT(IN)::CT(:,:)
REAL(DBL):: res(SIZE(vec))
INTEGER:: i

res=0.d0
DO i=1,SIZE(vec)
   IF(ABS(vec(i))>TOLERANCE)THEN  ! the component is not zero
		IF(follower(i)==1) THEN      ! the  component is following
			res=res+CT(:,i)*vec(i) 
		ELSE
			res(i)=res(i)+vec(i) ! for dead force
		ENDIF
	ENDIF
ENDDO

END FUNCTION TransferFollower

!************************************************************
ELEMENTAL FUNCTION InitPIAero(i) RESULT(res)

TYPE (PrescriInf)::res       
INTEGER, INTENT(IN) ::i
res%id=i
res%dof=(/1,1,1,1,1,1/)
res%value=(/0,0,0,0,0,0/)
res%time_fun_no=(/0,0,0,0,0,0/)
res%follower=(/0,0,0,0,0,0/)
res%value_current=0._DBL

END FUNCTION InitPIAero
!************************************************************

END MODULE PrescribedCondition
