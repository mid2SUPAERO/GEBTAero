!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************
! Programming Style
!------------------------
!1. Global constants are in CAPITALS
! 
!2. Fortran statements and keywords are in CAPITALS
!
!3. Intrinsic function names are in CAPITALS
! 
!4. User defined functions, subroutines, and modules, capitalize
!   the first letter of each word only, such as MySubroutine
!
!5. All other variables, including global and local variables, are 
!   in lower cases with underscore to separate words, eg. my_variable
!   
!6. Never use stop in a subroutine or function, but return an error message
!
!7. If at all possible, do not modify an input for functions, which means 
!   all the arguments for functions should have INTENT(IN)
!
!8. Allow only necessary information of a module accessible to an outsite procedure
!   use USE ModuleName, ONLY: a list of names
!
!9. Use PRIVATE to protect any data which will not pass information to outsite procedures.
!
!10.Use meaningful names whenever possible
!
!11.Include a brief data dictionary in the header of any program unit
!
!12.Echo any variable input from the end user so that the user can check the data
!
!13.Always initialize all variables in a program before using them through declaration (preferred), read, assignment
!
!14.Always use IMPLICIT NONE.
!
!15.Always include a default case in the case construction
!
!16.Never modify the value of a DO loop index variable while inside the loop. Never use the value of the
!   index variable after the DO loop completes normaly. 
!
!17.Use ES to get conventional scientific description.
!
!18.Begin every format associated with a WRITE statement with the 1X descriptor. Hence all the output will begin from a new line
!
!19.Break large programs into procedures whenever practical
!
!20.Declare the intent of every dummy argument in every procedure.
!
!21.Goto is only used to direct to return to the calling program due to critical errors 
!=========================================================================================

!---------------------------------------------------------------
!
! This module contains general-purpose global constants,
! I/O functions/subroutines and math functions/subroutines
!
! Gobal constants/variables
!----------------------------
! DBL,PI,DEG_2_RAD,RAD_2_DEG,TOLERANCE, I3,in_stat,allo_stat
!
! I/O functions/subroutines
!------------------------------
! FUNCTION FileOpen (file_unit,file_name,sta_type,rw_type,error)
! FUNCTION IOError(message,error)
! FUNCTION ItoChar(n)
! FUNCTION LowerCase(string) 
! FUNCTION  MemoryError(allo_stat,vari_name,error)
! SUBROUTINE TitlePrint(file_unit, title)
! FUNCTION UpperCase(string) 
! SUBROUTINE WriteError(EIN,error)
! SUBROUTINE WriteVec(file_unit,vec)
!
! math functions/subroutines
!------------------------------
! SUBROUTINE CT_THETA(theta,eCT,ekttek,eCTtheta)
! FUNCTION DirCosineTRodrigues(theta)
! SUBROUTINE Invert(matrix_in,matrix,vari_name,error)
! FUNCTION MATMUL3(mat,vec)
! FUNCTION Norm(vector)
! FUNCTION OuterProduct(vec1, vec2)
! FUNCTION Tilde(vect)
!---------------------------------------------------------------

!> This module contains general-purpose global constants,I/O functions/subroutines and math functions/subroutines

MODULE GlobalDataFun

IMPLICIT NONE

PRIVATE ! everything is private except declared by PUBLIC

PUBLIC DBL,PI,DEG_2_RAD,RAD_2_DEG,TOLERANCE,I3,FMT_INT,FMT_REAL,NDIM,NDOF_ND,NSTRN,e1
PUBLIC in_stat, allo_stat,MEMB_CONST,GRAV,NSTATES,Peters,RUNMOD,SOLVER,flutter_flag,MatA,ARPACK_MOD
PUBLIC EIGEN_OUTPUT,FLUTTER_LIMIT

PUBLIC FileOpen,IOError,MemoryError,TitlePrint,WriteError,WriteVec
PUBLIC CT_THETA,CT_THETA_T,DirCosineTRodrigues,Invert,MATMUL3,Norm,OuterProduct,Tilde,CrossProduct
PUBLIC Insert1DElement,Extract2DElement,MATMUL_sparse
INTERFACE WriteVec
 	  MODULE PROCEDURE WriteIntVector,WriteRealVector ! Write a vector
END INTERFACE


!=============================================================================
!
! Global constants/variables
!
!=============================================================================
INTEGER,PARAMETER:: NDIM=3         !< All the beams could behavior in the 3D space 
INTEGER,PARAMETER:: NDOF_ND=12     !< degrees of freedom per node/element is 12.
INTEGER,PARAMETER:: NSTRN=6        !< Number of strain measures/dofs in Timoshenko model 
INTEGER,PARAMETER:: MEMB_CONST=7   !< Number of labels needed for member properties
INTEGER,PARAMETER:: NSTATES= 6 !< Number of induces-flow states in Peters theory

INTEGER, PARAMETER:: DBL=SELECTED_REAL_KIND(15,307)

REAL(DBL),PARAMETER:: PI        =    3.1415926535897932D0    
REAL(DBL),PARAMETER:: DEG_2_RAD =    1.7453292519943296D-2 !< the ratio between radians and degrees
REAL(DBL),PARAMETER:: RAD_2_DEG =    5.7295779513082321D1  !< convert radian to degree
REAL(DBL),PARAMETER:: TOLERANCE = EPSILON(1.0_DBL)         !< a smart number of the double precision real number
REAL(DBL),PARAMETER:: GRAV      =    9.81                  !< gravity acceleration

REAL(DBL),PARAMETER::I3(3,3) = RESHAPE((/1.D0, 0.D0, 0.D0,& 
                                       & 0.D0, 1.D0, 0.D0,&
									   & 0.D0, 0.D0, 1.D0/),(/3,3/))  !< The 3x3 identity matrix
REAL(DBL),PARAMETER::e1(3)=(/1._DBL,0._DBL,0._DBL/)  !< The e1 unit vector

INTEGER:: in_stat      !< flag to indicate if the I/O process is successful: if positive, an error occured;
                       !< if negative, an end-of-file or end-of-record condition occurred; 
					   !< zero, no error, end-of-file, or end-of-record condition occurred.  
INTEGER:: allo_stat    !< flag to indicate status of allocating memory 


CHARACTER(*),PARAMETER :: FMT_REAL='ES15.7' !< format for output real numbers
CHARACTER(*),PARAMETER :: FMT_INT='I8'       !< format for output integer numbers

INTEGER:: RUNMOD=0  !< Define the output behavior of the program; 0: legacy mode of the computation code with output of a .out text file; 1: mode compatible with the python pre/postrpocessor (argument -p in the terminal), 2: silent mode (argument -s in the terminal)
INTEGER:: ARPACK_MOD=0  !< parameter WHICH of arpack solver (1:LI, 2:LM, 3:LR, 4:SR, default : LM in dnaupd and LI in dneupd) =>cf Arpack doc
INTEGER:: EIGEN_OUTPUT=0    !< define wich eigenvalue data to output (0: eigenvalues and eigenvectors, 1: eigenvalues only)
CHARACTER(10) :: SOLVER='MUMPS'       !< linear solver used (HSL : ddep.f, mc19.f + ma28 or MUMPS : linux library)
INTEGER:: flutter_flag=0        !< used in temporal simulation : 0= deformation are under a "flutter" state; 1= deformation are over a "flutter" state
REAL(DBL):: FLUTTER_LIMIT       !< the value of maximale angular deformaton use to trigger the flutter flag


!=========================================================================================


CONTAINS
!=================================================================
!
! The following are general purpose I/O functions/subroutines
!
!=================================================================

!************************************************************
!*                                                          *
!>    To open an old or new file for reading or writing     *
!*															*
!************************************************************
FUNCTION 	FileOpen (file_unit,file_name,sta_type,rw_type,error)

LOGICAL                   ::FileOpen
INTEGER,INTENT(IN)        ::file_unit !<File Unit (see fortran IO doc)
CHARACTER(*),INTENT(IN)   ::file_name
CHARACTER(*),INTENT(IN)   ::sta_type    !< status type
CHARACTER(*),INTENT(IN)   ::rw_type     !<rewrite configuration
CHARACTER(*),INTENT(OUT)  ::error   !<#ioaero::error

error=''
FileOpen=.FALSE.

OPEN (UNIT=file_unit, file=file_name,STATUS=sta_type,ACTION = rw_type,IOSTAT=in_stat)

IF (in_stat/=0) THEN
  IF(rw_type=='READ') error='Cannot open the file '//TRIM(file_name)//' for reading!'
  IF(rw_type=='WRITE')error='Cannot open the file '//TRIM(file_name)//' for writing!'
ENDIF

IF(error/='')FileOpen=.TRUE.

END FUNCTION FileOpen
!***********************************************************



!************************************************************
!*                                                          *
!>        Check the error of I/O processing                 *
!*															*
!************************************************************
FUNCTION  IOError(message,error)

LOGICAL                 ::IOError
CHARACTER(*),INTENT(IN) ::message        !< a character variable to hold error message
CHARACTER(*),INTENT(OUT)::error !<#ioaero::error

error=''
IOError=.FALSE.

IF(in_stat/=0) THEN 
    error='I/O error: '//TRIM(message)
    IOError=.TRUE.
ENDIF

END FUNCTION IOError
!***********************************************************



!************************************************************
!*                                                          *
!>         Convert an integer to character                  *
!*                                                          *
!************************************************************
FUNCTION ItoChar(n) RESULT(char)

	INTEGER,INTENT(IN):: n
	CHARACTER(20):: char

	 WRITE(char, *) n

END FUNCTION ItoChar
!***********************************************************



!~ !************************************************************
!~ !*                                                          *
!~ !*      Convert a string or character to lower case         *
!~ !*                                                          *
!~ !************************************************************
!~ FUNCTION LowerCase(string) RESULT(lc_string)

!~ CHARACTER(*),INTENT(IN):: string
!~ CHARACTER(LEN=LEN(string)):: lc_string

!~ CHARACTER(LEN=26),PARAMETER:: UPPER='ABCDEFGHIJKLMNOPQRSTUVWXYZ',&
!~ 	                          lower='abcdefghijklmnopqrstuvwxyz'

!~ INTEGER:: k  ! loop counter
!~ INTEGER::loc ! position in alphabet

!~ lc_string=string
    
!~ DO k=1,len(string)
!~ 	loc=INDEX(UPPER,string(k:k))
!~ 	IF(loc/=0)lc_string(k:k)=lower(loc:loc)
!~ ENDDO

!~ END FUNCTION LowerCase
!~ !************************************************************



!************************************************************
!*                                                          *
!>        Check the error of memory allocation              *
!*															*
!************************************************************
FUNCTION  MemoryError(vari_name,error)

LOGICAL                 ::MemoryError
CHARACTER(*),INTENT(IN) ::vari_name         !< a character variable to hold variable name
CHARACTER(*),INTENT(OUT)  ::error   !<#ioaero::error

error=''
MemoryError=.FALSE.

IF(allo_stat/=0) THEN
	error='Memory error: allocate '//TRIM(vari_name)
    MemoryError=.TRUE.
ENDIF


END FUNCTION MemoryError
!************************************************************



!************************************************************
!*                                                          *
!>        To print a title for a block of data              *
!*															*
!************************************************************
SUBROUTINE TitlePrint(file_unit, title)

INTEGER,INTENT(IN)        ::file_unit 
CHARACTER(*),INTENT(IN)   ::title

WRITE(file_unit,*) 
WRITE(file_unit,*) 
WRITE(file_unit,'(1x,100A)') title
WRITE(file_unit,*)'========================================================'

END SUBROUTINE TitlePrint
!***********************************************************



!~ !************************************************************
!~ !*                                                          *
!~ !*    Convert a string or character to upper case           *
!~ !*                                                          *
!~ !************************************************************
!~ FUNCTION UpperCase(string) RESULT(uc_string)

!~ CHARACTER(*),INTENT(IN):: string
!~ CHARACTER(LEN=LEN(string)):: uc_string

!~ CHARACTER(LEN=26),PARAMETER:: UPPER='ABCDEFGHIJKLMNOPQRSTUVWXYZ',&
!~ 		                      lower='abcdefghijklmnopqrstuvwxyz'

!~ INTEGER:: k  ! loop counter
!~ INTEGER::loc ! position in alphabet

!~ uc_string=string
    
!~ DO k=1,len(string)
!~ 	loc=INDEX(lower,string(k:k))
!~ 	IF(loc/=0)uc_string(k:k)=UPPER(loc:loc)
!~ ENDDO

!~ END FUNCTION UpperCase
!~ !************************************************************



!************************************************************
!*                                                          *
!>    Write error to the echo file                          *
!*															*
!************************************************************
SUBROUTINE WriteError(EIN,error)

INTEGER,INTENT(IN)::EIN !< file unit to write the error message
CHARACTER(*),INTENT(IN)  ::error

LOGICAL file_opened

INQUIRE (EIN,  OPENED = file_opened) ! Check whether the file is already opened, if yes, then dump the error message to this file

IF(file_opened)THEN

	WRITE(EIN,*) 
	IF(error/='')THEN

		CALL TitlePrint(EIN, 'Error Message')
		WRITE(EIN,'(1x, 300A)') error 

	ELSE

		WRITE(EIN,*) 'Congratulations! No errors!'
	ENDIF

	CLOSE(EIN)
ENDIF

END SUBROUTINE WriteError
!********************************************************



!************************************************************
!*                                                          *
!>  Write an integer vector to the file_unit                *
!*															*
!************************************************************
SUBROUTINE WriteIntVector(file_unit,vec)

INTEGER,INTENT(IN)::file_unit   !<File unit to write the vector
INTEGER,INTENT(IN)::vec(:)

WRITE(file_unit,'(1x,'//TRIM(ItoChar(SIZE(vec)))//FMT_INT//')')vec

END SUBROUTINE WriteIntVector
!********************************************************



!************************************************************
!*                                                          *
!>    Write a real vector to the file_unit                  *
!*															*
!************************************************************
SUBROUTINE WriteRealVector(file_unit,vec)

INTEGER,INTENT(IN)::file_unit   !<File unit to write the vector
REAL(DBL),INTENT(IN)::vec(:)
!REAL(DBL)::tmp(SIZE(vec)),vec_norm

!tmp=vec
!vec_norm=Norm(vec)

!IF(vec_norm>TOLERANCE) THEN
!	WHERE(ABS(tmp/vec_norm)<TOLERANCE)tmp=0.0D0  ! not output components which are negligible comparing to the biggest term in the vector
!ELSE
!	tmp=0.0D0
!ENDIF

WRITE(file_unit,'(1x,'//TRIM(ItoChar(SIZE(vec)))//FMT_REAL//')')vec

END SUBROUTINE WriteRealVector
!********************************************************



!===================================================================
!
! The following are general purpose math functions/subroutines
!
!==================================================================

!*************************************************************
!*                                                           *   
!>  Calculate eC^T derivative w.r.t theta                     *
!!  return derivative and ekttek                             *
!*                                                           *   
!*************************************************************
SUBROUTINE CT_THETA(theta,eCT,ekttek,eCTtheta)

REAL(DBL),INTENT(IN):: theta(:) !<Rodrigues rotation parameters
REAL(DBL),INTENT(IN):: eCT(:,:) !< Direction Cosine matrix between frame b and B
REAL(DBL),INTENT(OUT)::ekttek(:,:,:) !< =OuterProduct(ek,theta)+OuterProduct(theta,ek)
REAL(DBL),INTENT(OUT)::eCTtheta(:,:,:) !< Derivatives of eCT relative to theta

INTEGER:: k
REAL(DBL)::temp,ek(3)

temp=1.D0/(1.D0+DOT_PRODUCT(theta,theta)*0.25D0)

DO k=1,3  
    ek=I3(k,:)
	ekttek(k,:,:)=OuterProduct(ek,theta)+OuterProduct(theta,ek)
    eCTtheta(k,:,:)=( (ekttek(k,:,:)-theta(k)*(I3+eCT) )*0.5D0 +Tilde(ek) )*temp
ENDDO

END SUBROUTINE CT_Theta
!**********************************************************


!*************************************************************
!*                                                           *   
!>  Calculate \f$ \dot{eC^T}.x \f$ derivative w.r.t \f$ \dot{theta}  \f$    *
!*                                                           *   
!*************************************************************
FUNCTION CT_THETA_T(theta,eCT,x) RESULT(res)

REAL(DBL),INTENT(IN):: theta(:) !<Rodrigues rotation parameters
REAL(DBL),INTENT(IN):: eCT(:,:) !< Direction Cosine matrix between frame b and B
REAL(DBL),INTENT(IN):: x(:) !<  vector
REAL(DBL)            ::res(3,3)

res=( (DOT_PRODUCT(theta,x)*I3+OuterProduct(theta,x)-OuterProduct(MATMUL(I3+eCT,x),theta) )*.5d0 &
  &   -Tilde(x) )/(1.D0+DOT_PRODUCT(theta,theta)*0.25D0) 

END FUNCTION CT_Theta_T
!**********************************************************


!*************************************************************
!*                                                           *   
!>  Calculate the transpose of the direction cosine in       *
!!  terms of rodrigues parameters                            *
!*                                                           *
!*===========================================================*
!* Input:                                                    *
!*  theta   --	  rodrigues parameters                       *
!* Output	                                                 *
!*  DirCosine_Rodrigues --the transpose 3X3 direction        *
!*                         cosine matrix                     *
!*************************************************************
FUNCTION DirCosineTRodrigues(theta)

REAL(DBL),INTENT(IN):: theta(:)
REAL(DBL)::DirCosineTRodrigues(3,3) 
REAL(DBL):: t2_4

t2_4=DOT_PRODUCT(theta,theta)*0.25D0

DirCosineTRodrigues=((1.0D0-t2_4)*I3+Tilde(theta)+OuterProduct(theta,theta)*0.5D0)/(1.0d0+t2_4)

END FUNCTION DirCosineTRodrigues
!**********************************************************


!***********************************************
!*                                             *
!>       Invert a small square matrix          *
!*                                             * 
!***********************************************
SUBROUTINE Invert(matrix_in,matrix,vari_name,error)
 
REAL(DBL),INTENT(IN)::matrix_in(:,:) !< the matrix to be inverted
CHARACTER(*),INTENT(IN)::vari_name

REAL(DBL),INTENT(OUT)::matrix(:,:)    !< the inverse of the matrix

CHARACTER(*),INTENT(OUT)::error

INTEGER::i,k,n 
REAL(DBL)::con,diag_sum,zero

matrix=matrix_in
  
n= SIZE(matrix,1)

diag_sum=0.0d0

DO i=1,n
	diag_sum=diag_sum+matrix(i,i)
ENDDO

zero=TOLERANCE*diag_sum/n
DO k=1,n
   con=matrix(k,k);
   
    IF(ABS(CON)<zero) THEN
        WRITE(0,*) "GlobalDataFun.f90 : ",error
        ERROR STOP 111
    ENDIF

   matrix(k,k)=1.d0
   matrix(k,:)=matrix(k,:)/con
   
   DO i=1,n
      IF(i/=k) THEN
         con=matrix(i,k); matrix(i,k)=0.d0
         matrix(i,:)=matrix(i,:) - matrix(k,:)*con
      END IF
   ENDDO
   
ENDDO

9999 RETURN

END SUBROUTINE
!***********************************************************



!*************************************************************
!*                                                           *   
!>  Multiply a rank 3 matrix with a vector with every colum  *
!! of the resulting matrix is equal to the multiplication    *
!*                                                           *
!*************************************************************
FUNCTION MATMUL3(mat,vec)

REAL(DBL),INTENT(IN):: mat(:,:,:),vec(:)
REAL(DBL)::MATMUL3(SIZE(mat,1),SIZE(mat,2))

INTEGER::i

DO i=1,SIZE(mat,1) 

	MATMUL3(:,i)=MATMUL(mat(i,:,:),vec)

ENDDO


END FUNCTION MATMUL3
!**********************************************************



!*************************************************************
!*                                                           *   
!>  Calculate the L2 norm of a real vector                   *
!*                                                           *
!*************************************************************
FUNCTION Norm(vector)

REAL(DBL),INTENT(IN):: vector(:)
REAL(DBL)::Norm 

Norm=SQRT(DOT_PRODUCT(vector,vector))

END FUNCTION Norm
!**********************************************************



!*************************************************************
!*                                                           *   
!>  Calculate the outer product of two vectors               *
!*                                                           *
!*************************************************************
FUNCTION OuterProduct(vec1, vec2)

REAL(DBL),INTENT(IN):: vec1(:), vec2(:)
REAL(DBL)::OuterProduct(SIZE(vec1),SIZE(vec2)) 

INTEGER:: i,j,n1,n2

n1=SIZE(vec1); n2=SIZE(vec2)

DO i=1, n1
	DO j=1, n2
		OuterProduct(i,j)=vec1(i)*vec2(j)
	ENDDO
ENDDO

END FUNCTION OuterProduct
!**********************************************************



!*************************************************************
!*                                                           *   
!>  Carry out the tilde operation for a real vector          *
!*                                                           *
!*===========================================================*
!* Input:                                                    *
!*  vect     --	  a vector with three components             *
!* Output	                                                 *
!*  tilde    --  the 3X3 antisymmetric matrix                *
!*************************************************************
FUNCTION Tilde(vect)

REAL(DBL),INTENT(IN):: vect(3)
REAL(DBL):: Tilde(3,3) 

Tilde=0.0d0

Tilde(1,2)= -vect(3)
Tilde(1,3)=  vect(2)
Tilde(2,1)=  vect(3)
Tilde(2,3)= -vect(1)
Tilde(3,1)= -vect(2)
Tilde(3,2)=  vect(1)

END FUNCTION
!**********************************************************


!*************************************************************
!*                                                           *   
!>  Carry out cross product of two real vectors              *
!*                                                           *
!*===========================================================*
!* Input:                                                    *
!*  a,b     --	  two vectors with three components          *
!* Output	                                                 *
!*  CrossProduct  --  the result vector                      *
!*************************************************************
FUNCTION CrossProduct(a,b)

REAL(DBL),INTENT(IN):: a(3), b(3)
REAL(DBL):: CrossProduct(3) 


CrossProduct(1)=a(2)*b(3)-a(3)*b(2)
CrossProduct(2)=a(3)*b(1)-a(1)*b(3)
CrossProduct(3)=a(1)*b(2)-a(2)*b(1)

END FUNCTION
!**********************************************************



!*******************************************************************
!> Insert a real matrix into the 1D coefficient matrix
!*******************************************************************
SUBROUTINE Insert1DElement(nz,tmpR,irn,jcn,elemCoef1D,str_r1,str_c1,str_r2,str_c2,str_r3,str_c3,str_r4,str_c4)
   
IMPLICIT NONE

INTEGER,INTENT(INOUT)::nz,irn(:),jcn(:)
REAL(DBL),INTENT(INOUT)::elemCoef1D(:)
REAL(DBL),INTENT(IN)::tmpR(:,:) 
INTEGER,INTENT(IN)::str_r1,str_c1
INTEGER,OPTIONAL,INTENT(IN)::str_r2,str_c2 ! these two values should appear together, 
                                           ! implying the same matrix fill in 2nd block
INTEGER,OPTIONAL,INTENT(IN)::str_r3,str_c3 ! these two values should appear together, 
                                           ! implying the same matrix fill in the 3rd block

INTEGER,OPTIONAL,INTENT(IN)::str_r4,str_c4 ! these two values should appear together, 
                                           ! implying the same matrix fill in the 4th block

INTEGER:: i,j
DO i=1,SIZE(tmpR,1)
DO j=1,SIZE(tmpR,2)
	IF(ABS(tmpR(i,j))>TOLERANCE) THEN
	     nz=nz+1
         irn(nz)=str_r1+i
         jcn(nz)=str_c1+j
         elemCoef1D(nz)=tmpR(i,j)
         
		 IF(PRESENT(str_r2)) THEN
			nz=nz+1
			irn(nz)=str_r2+i
			jcn(nz)=str_c2+j
			elemCoef1D(nz)=tmpR(i,j)
	     ENDIF
		 IF(PRESENT(str_r3)) THEN
			nz=nz+1
			irn(nz)=str_r3+i
			jcn(nz)=str_c3+j
			elemCoef1D(nz)=tmpR(i,j)
	     ENDIF
		 IF(PRESENT(str_r4)) THEN
			nz=nz+1
			irn(nz)=str_r4+i
			jcn(nz)=str_c4+j
			elemCoef1D(nz)=tmpR(i,j)
	     ENDIF
      ENDIF
ENDDO
ENDDO

END SUBROUTINE Insert1DElement
!*************************************************************


!*******************************************************************
!> Back a 2D array from the 1D coefficient matrix
!*******************************************************************
SUBROUTINE Extract2DElement(nz,irn,jcn,elemCoef1D,tmpR,str_r1,str_c1)

IMPLICIT NONE

INTEGER,INTENT(IN)::nz,irn(:),jcn(:)
REAL(DBL),INTENT(IN)::elemCoef1D(:)
REAL(DBL),INTENT(OUT)::tmpR(:,:) 
INTEGER,INTENT(IN)::str_r1,str_c1

INTEGER:: i,j,ii, irow,jcol

tmpR=0.0D0
DO i=1,SIZE(tmpR,1)
	irow=i+str_r1

	DO j=1,SIZE(tmpR,2)
	    jcol=j+str_c1
		DO ii=1,nz
			IF(irn(ii)==irow.AND.jcn(ii)==jcol) THEN
				tmpR(i,j)=elemCoef1D(ii)
				EXIT
			ENDIF
		ENDDO
	ENDDO
ENDDO

END SUBROUTINE Extract2DElement
!****************************************

!**************************************************
!> Matmul(vector, matrix)
!! with matrix stored in a spare format
!*****************************************************
FUNCTION MATMUL_sparse(vector,nsize,ne,irn,jcn,matrix1D) RESULT (res)

IMPLICIT NONE

REAL(DBL),INTENT(IN)::vector(:),matrix1D(:)
INTEGER,INTENT(IN)::nsize,ne,irn(:),jcn(:)
REAL(DBL):: res(nsize)
INTEGER:: i,j,k

res=0.0D0

DO k = 1,ne
   i = irn(k)
   j = jcn(k)
   res(j) = res(j)+ vector(i)*matrix1D(k)
ENDDO
 
END FUNCTION MATMUL_sparse
!*****************************************************************

!> Functions used for the Finite State Unsteady Thin Airfoil Theory of Peters
!! see Introduction to Structural Dynamics and Aeroelasticity by Hodges p 139
FUNCTION Peters(n) RESULT(res)
INTEGER, INTENT(IN) ::n !<Number of induced flow states (Ns)
REAL(DBL) :: res(n,2*n+2)   !< Matrix containing ordered in line : the A matrix, the B vector the C vector, the Identity matrix
INTEGER   :: i
res =0.D0
res(:,1:n) = MatA(n)
res(:,n+1) = VecB(n)
res(:,n+2) = VecC(n)

DO i=1,n
	res(i,n+2+i) = 1
ENDDO

END FUNCTION Peters

!> Compute the Peters A matrix
FUNCTION MatA(n) RESULT(res)
IMPLICIT NONE

INTEGER, INTENT(IN) ::n
REAL(DBL) :: res(n,n)

res = 0.D0
res = MatD(n)+Prod(VecD(n),VecB(n),n)+Prod(VecC(n),VecD(n),n)+0.5D0*Prod(VecC(n),VecB(n),n)
END FUNCTION MatA

!> Compute the Peters D matrix
FUNCTION MatD(n) RESULT(res)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
INTEGER             :: i,j
REAL(DBL) :: res(n,n)

res = 0.D0
DO i=1,n
	DO j=1,n
	IF(i==j+1) res(i,j)=1.D0/(2.D0*i)
	IF(i==j-1) res(i,j)=-1.D0/(2.D0*i)
	ENDDO
ENDDO

END FUNCTION MatD

!> Compute the Peters B vector
FUNCTION VecB(n) RESULT(res)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
INTEGER             :: i,j
REAL(DBL) :: res(n)

res = 0.D0
DO i=1,n-1
	res(i) = (-1.0)**(i-1)
	DO j=n-i,n+i-1
		res(i) = res(i)*j
	ENDDO
	res(i)=res(i)/(Factoriel(i)**2)
ENDDO
res(n) = (-1.0)**(i-1)




END FUNCTION VecB

!> Compute Peters C vector
FUNCTION VecC(n) RESULT(res)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
INTEGER             :: i
REAL(DBL) :: res(n)

res = 0.D0
DO i=1,n
	res(i)=2.0/i
ENDDO

END FUNCTION VecC

!> Compute Peters D vector
FUNCTION VecD(n) RESULT(res)
IMPLICIT NONE

INTEGER, INTENT(IN) ::n
REAL(DBL) :: res(n)

res = 0.D0
res(1) = 0.5

END FUNCTION VecD

!> Compute the product of two square matrix
FUNCTION Prod(Vec1,Vec2,n) RESULT(res)
INTEGER,INTENT(IN)  :: n
REAL(DBL),INTENT(IN):: Vec1(n),Vec2(n)
REAL(DBL)           :: res(n,n)
INTEGER             :: i,j

res = 0.D0
DO i=1,n
	DO j=1,n
	res(i,j)=Vec1(i)*Vec2(j)
	ENDDO
ENDDO

END FUNCTION Prod

!> Compute the value of factoriel(n)
FUNCTION Factoriel(n) RESULT(res)
INTEGER,INTENT(IN)  :: n
INTEGER             :: i
REAL(DBL)           :: res

res = 1
DO i=2,n
res = res*i
ENDDO


END FUNCTION Factoriel

END MODULE GlobalDataFun
!=============================================================================






