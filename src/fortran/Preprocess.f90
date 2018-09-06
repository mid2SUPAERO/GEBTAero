!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************

!***************************************************************
!*                                                             *
!* This module preprocess the finite element model including   *
!* connectivity and member information. This information are   *
!* time step indepedent                                        *
!***************************************************************
MODULE PreproModule
USE InternalData

IMPLICIT NONE

PRIVATE       ! So everything is private, except declared by PUBLIC
PUBLIC Preprocess 


INTEGER::ndiv,ncol_memb
!=============================================
CONTAINS


!***********************************************************
!*                                                         *
!* Obtaining the connection condition for each key point   *
!* if a point is connected to more than one member, it is  *
!* a connection point, otherwise it is a boundary point.   *
!* It also calculates the size of the problem              *
!*                                                         *
!***********************************************************
SUBROUTINE Preprocess(nkp,nelem,ndof_el,member,material,frame,coord,curvature,&
                    & dof_con,memb_info,error,aero_flag,grav_flag,aerodyn_coef)

IMPLICIT NONE

INTEGER,INTENT(IN)   ::nkp,nelem,aero_flag,grav_flag
INTEGER,INTENT(IN)   ::ndof_el
INTEGER,INTENT(IN)   ::member(:,:)
REAL(DBL),INTENT(INOUT)::material(:,:,:)
REAL(DBL),INTENT(IN)   ::frame(:,:,:),coord(:,:),curvature(:,:)
REAL(DBL),INTENT(IN),OPTIONAL :: aerodyn_coef(:,:)

INTEGER,INTENT(OUT)  ::dof_con(:)
TYPE (MemberInf),INTENT(OUT)::memb_info(:)

CHARACTER(*),INTENT(OUT)::error

REAL(DBL)::tempR(NSTRN,NSTRN)
INTEGER:: i           

INTEGER:: ntimes,nmemb  ! times a point appears in some array
INTEGER:: nrow,ncol

!--------------------------------------------------
! Obtain the inverse of the stiffness matrix, and
! store in the original position
!------------------------------------------
!DO i=1,SIZE(material,1)
!	CALL Invert(material(i,1:NSTRN,:),tempR,'Stiffness Matrix',error)
!	IF(error/='') GOTO 9999
!	material(i,1:NSTRN,:)=tempR
!ENDDO

!--------------------------------------------------
! Obtain the inverse of the mass matrix, and
! store in the original position
!------------------------------------------
IF(ndof_el>=18)THEN
	DO i=1,SIZE(material,1)
		CALL Invert(material(i,NSTRN+1:12,:),tempR,'MASS MATRIX',error)
		material(i,NSTRN+1:12,:)=tempR
	ENDDO
ENDIF

!----------------------------------
!Connectivity of the key points
!----------------------------------
dof_con=0 ! initialize

DO i=1,nkp
	ntimes=COUNT(i==member(:,1:2)) 

	IF(ntimes==0) THEN  
		WRITE(error,*)'The point ',i,' is not connected to any member!'
        WRITE(0,*) 'Preprocess.f90 : ',error
        ERROR STOP 112
	ELSE 
		IF(ntimes/=1) dof_con(i)=1  ! a connection point 
	ENDIF
ENDDO

nmemb=SIZE(member,1)

IF(ALLOCATED(index_mb))DEALLOCATE(index_mb)
IF(ALLOCATED(index_kp)) DEALLOCATE (index_kp)

ALLOCATE(index_kp(nkp,2),STAT=allo_stat)
IF(MemoryError('index_kp',error)) GOTO 9999
index_kp=0

ALLOCATE(index_mb(nmemb,2),STAT=allo_stat)
IF(MemoryError('index_mb',error)) GOTO 9999
index_mb=0

! Introduce the new assembly based on the starting row and column of the key points
! and first element of the member
!------------------------------------------------------------------------------------------
nrow=1; ncol=1

DO i=1,nmemb
   
   IF(index_kp(member(i,1),1)==0) THEN  ! The point is not used yet
		index_kp(member(i,1),1)=nrow
		index_kp(member(i,1),2)=ncol      
		nrow=nrow+NDOF_ND ! every point add NDOF_ND equations
	    ncol=ncol+NSTRN ! every point add NSTRN unknowns
   ELSE  ! The point has been used previously, it must connect to more than one members
        nrow=nrow+NSTRN ! every previously used connection point adds six more equations
   ENDIF

   index_mb(i,1)=nrow  ! member starting row excluding the starting NDOF_ND related with starting point
   index_mb(i,2)=ncol  ! member starting column 
   
   ndiv=member(i,6)
   ncol_memb=ndiv*ndof_el
   nrow=nrow+ncol_memb-NDOF_ND  ! every member add ndiv*ndof_el-2*NDOF_ND equations 
   ncol=ncol+ncol_memb          ! every member add ndiv*ndof_el unknowns
 
   IF(index_kp(member(i,2),1)==0) THEN ! The point is not used yet
		index_kp(member(i,2),1)=nrow
		index_kp(member(i,2),2)=ncol    
		nrow=nrow+NDOF_ND  ! every point add NDOF_ND equations
	    ncol=ncol+NSTRN ! every point add NSTRN unknowns
   ELSE
        nrow=nrow+NSTRN
   ENDIF

   CALL MemberProperties(i,ndof_el,member,material,frame,coord,curvature,&
                       & memb_info(i),error,aero_flag,grav_flag,aerodyn_coef)

ENDDO

neMax=nelem*nzElemMax+nkp*36 ! the maximum number of nonzero entries

IF(nrow/=ncol) THEN
	error='Something wrong with the mesh'
    WRITE(0,*) 'Preprocess.f90 : ',error
    ERROR STOP 113
ELSE
	nsize=nrow-1
ENDIF

9999 IF(error/='') THEN
    WRITE(0,*) 'Preprocess.f90 : ',error
    ERROR STOP 114
ENDIF

END SUBROUTINE Preprocess
!***********************************************************



!*************************************************************
!*                                                           *   
!* Extract member properties for each division               *
!*                                                           *
!*===========================================================*
!* Inputs:                                                   *
!* Output	                                                 *
!*************************************************************
	SUBROUTINE MemberProperties(memb_no,ndof_el,member,material,frame,coord,curvature,&
                          & memb_info_i,error,aero_flag,grav_flag,aerodyn_coef)

IMPLICIT NONE

INTEGER,INTENT(IN)  ::memb_no,ndof_el,member(:,:),aero_flag,grav_flag
REAL(DBL),INTENT(IN)::material(:,:,:),frame(:,:,:),coord(:,:),curvature(:,:)
REAL(DBL),INTENT(IN),OPTIONAL :: aerodyn_coef(:,:)

TYPE (MemberInf),INTENT(OUT)::memb_info_i
CHARACTER(*),INTENT(OUT)::error

INTEGER:: j,tmpN,mat1,mat2,frame_no,curve_no
REAL(DBL)::mCab(3,3),mCoord(2,3),mL,dL
REAL(DBL)::mCurv(3),kn
REAL(DBL)::k12,kn2,kn4,kkTkn2(3,3),tmp33_1(3,3),tmp33_2(3,3),knx1

memb_info_i%ndiv     =ndiv
memb_info_i%ncol_memb=ncol_memb

!-------------------------------
! Member sectional properties
!-------------------------------	
tmpN=NSTRN
IF(ndof_el>=18) tmpN=2*NSTRN

ALLOCATE(memb_info_i%mate(ndiv,tmpN,NSTRN),STAT=allo_stat)
IF(MemoryError('memb_info_i%mate',error)) GOTO 9999

mat1=member(memb_no,3);	mat2=member(memb_no,4)
DO j=1,ndiv
	IF(mat1==mat2)THEN
		memb_info_i%mate(j,:,:)=material(mat1,:,:)
	ELSE
		memb_info_i%mate(j,:,:)=material(mat1,:,:)+(j-0.5D0)/ndiv* &
		        &   (material(mat2,:,:)-material(mat1,:,:))/ndiv
	ENDIF
ENDDO
!---------------------------------------------------
IF (aero_flag>0 .OR. grav_flag>0) THEN
	ALLOCATE(memb_info_i%aerodyn_coef(ndiv,8),STAT=allo_stat)
	IF(MemoryError('memb_info_i%aerodyn_coef',error)) GOTO 9999
ENDIF
!aerodynamics coef
IF(aero_flag /=0 .OR. grav_flag/=0) THEN
	DO j=1,ndiv
		IF(mat1==mat2)THEN
			memb_info_i%aerodyn_coef(j,:)=aerodyn_coef(mat1,:)
		ELSE
			memb_info_i%aerodyn_coef(j,:)=aerodyn_coef(mat1,:)+(j-0.5D0)/ndiv* &
					&   (aerodyn_coef(mat2,:)-aerodyn_coef(mat1,:))/ndiv
		ENDIF
	ENDDO
ENDIF




!------------------------------------------------
! Length of the member, divisions, and ending
! arc length of the division
!------------------------------------------------
mCoord(1,:)=coord(member(memb_no,1),:) ! starting point of this member
mCoord(2,:)=coord(member(memb_no,2),:) ! ending point of this member
xyz_pt1=coord(member(1,1),:)       ! starting point of the first member

mL=Norm(mCoord(2,:)-mCoord(1,:))  !length of the member for prismatic beams, the distance between the end points: |r_a-r_0|

curve_no=member(memb_no,7)
IF(curve_no/=0) THEN
	mCurv=curvature(curve_no,:)
	kn=Norm(mCurv)
	IF(kn/=0.0d0) THEN
	   IF(mCurv(1)==0.0D0) THEN
			mL=2.0D0*ASIN(0.5D0*kn*mL)/kn
	   ELSE IF(mCurv(2)/=0.0D0.OR.mCurv(3)/=0.0D0) THEN
		   k12=mCurv(1)*mCurv(1)
		   kn2=kn*kn
		   kn4=kn2*kn2
		   mL=Rtbis(CurveBeamFun,kn,mL,kn2,k12,kn4,mL,mL*kn/mCurv(1),TOLERANCE*100,100,error)
	   ENDIF
     ENDIF
ENDIF

dL=mL/ndiv
memb_info_i%dL=dL

ALLOCATE(memb_info_i%Le(ndiv),STAT=allo_stat)
IF(MemoryError('memb_info_i%Le',error)) GOTO 9999

DO j=1,ndiv
	memb_info_i%Le(j)=j*dL
ENDDO
!-----------------------------------------------------


!------------------------------------------------
! Frame b at the middle point of each division
! and coordinate at the middle point of each division
!------------------------------------------------
frame_no=member(memb_no,5)
IF(frame_no/=0) THEN 
	mCab =frame(frame_no,:,:)
ELSE 
	mCab=I3 ! default global frame
ENDIF

ALLOCATE(memb_info_i%triad(ndiv,3,3),STAT=allo_stat)
IF(MemoryError('memb_info_i%triad',error)) GOTO 9999

ALLOCATE(memb_info_i%coordinate(ndiv,3),STAT=allo_stat)
IF(MemoryError('memb_info_i%coordinate',error)) GOTO 9999

IF(curve_no/=0) THEN
    kkTkn2=OuterProduct(mCurv,mCurv)/(kn*kn)
	tmp33_1=I3-kkTkn2
	tmp33_2=Tilde(mCurv)/kn
	DO j=1,ndiv
		knx1=kn*(memb_info_i%Le(j)-0.5d0*dL)
		memb_info_i%triad(j,:,:)=MATMUL(mCab,tmp33_1*Cos(knx1)+tmp33_2*Sin(knx1)+kkTkn2)
		memb_info_i%coordinate(j,:)=mCoord(1,:)+MATMUL(mCab, MATMUL(tmp33_1*Sin(knx1) & 
		                        &  +tmp33_2*(1.0D0-Cos(knx1))+kkTkn2*knx1,e1/kn))
	ENDDO
ELSE
	DO j=1,ndiv
		memb_info_i%triad(j,:,:)=mCab
		memb_info_i%coordinate(j,:)=mCoord(1,:)+(j-0.5D0)/ndiv*(mCoord(2,:)-mCoord(1,:) )
	ENDDO
ENDIF
!------------------------------------------------------

9999 IF(error/='') THEN
    WRITE(0,*) 'Preprocess.f90 : ',error
    ERROR STOP 115
ENDIF

END SUBROUTINE MemberProperties
!*************************************************************


!************************************************************
! Function for evaluating arc length of initially curved and twisted beams
!****************************************************************
FUNCTION  CurveBeamFun(kn,mL,kn2,k12,kn4,xvar)
IMPLICIT NONE
REAL(DBL),INTENT(IN)::kn,mL,kn2,k12,kn4
REAL(DBL),INTENT(IN)::xvar
REAL(DBL)			::CurveBeamFun
REAL(DBL)           ::knx
	
knx=kn*xvar
CurveBeamFun=mL*mL-( 2.0D0*(kn2-k12)*(1.0D0-COS(knx))+k12*knx*knx )/kn4

END FUNCTION  CurveBeamFun
!**************************************************************************





!*************************************************************
!*                                                           *   
!*  Use biosection to find root of a function                *
!* from the book of Numerical Recipes                        *
!*                                                           *
!*===========================================================*
!* Inputs:                                                   *
!*  func     --	  the given nonlinear equation func==0.0d0   *
!*  x1, x2   --   the root is known to lie between x1, x2    *
!*  xacc     --   a small number indicating the accuracy     *
!*  maxit    --   max number of iterations                   *
!* Output	                                                 *
!*  root     --  the root                                    *
!*************************************************************
FUNCTION Rtbis(func,kn,mL,kn2,k12,kn4,x1,x2,xacc,maxit,error) RESULT(root)

IMPLICIT NONE
REAL(DBL),INTENT(IN)::kn,mL,kn2,k12,kn4
REAL(DBL),INTENT(IN)	::x1,x2,xacc
INTEGER,INTENT(IN)		::maxit
CHARACTER(*),INTENT(OUT)::error
REAL(DBL)			    ::root

INTERFACE
    FUNCTION func(kn,mL,kn2,k12,kn4,x)
	    USE GlobalDataFun
		IMPLICIT NONE
		REAL(DBL),INTENT(IN)::kn,mL,kn2,k12,kn4
		REAL(DBL),INTENT(IN)::x
		REAL(DBL)           ::func
	END FUNCTION func
END INTERFACE

INTEGER  :: j
REAL(DBL):: dx,f,fmid,xmid

fmid=func(kn,mL,kn2,k12,kn4,x2)
f=func(kn,mL,kn2,k12,kn4,x1)
IF(f*fmid>=0.0D0) THEN
    error='rtbis: root is not bracketed.'
    WRITE(0,*) 'Preprocess.f90 : ',error
    ERROR STOP 116
ENDIF    
IF(f<0.0D0)THEN
	root=x1
	dx=x2-x1
ELSE
	root=x2
	dx=x1-x2
ENDIF
DO j=1,maxit
	dx=dx*0.5D0
	xmid=root+dx
	fmid=func(kn,mL,kn2,k12,kn4,xmid)
	IF(fmid<=0.0D0)root=xmid
	IF(ABS(dx)<xacc.OR.fmid==0.0D0)RETURN
ENDDO

error='Rtbis: the number of bisection exceeds the allowed.'
WRITE(0,*) 'Preprocess.f90 : ',error
ERROR STOP 117

END FUNCTION Rtbis
!**********************************************************


END MODULE PreproModule
!*********************************************
