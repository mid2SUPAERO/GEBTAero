!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************
!=================================================================
!
! A module for defining time functions needed for both prescribed  
! concentrated and distributed conditions 
!
!=================================================================
MODULE TimeFunctionModule

USE GlobalDataFun

IMPLICIT NONE

PRIVATE       ! So everything is private, except declared by PUBLIC
PUBLIC TimeFunction,GetTimeFunction,InitTF,InputEchoTimeFunctions,CurrentValues


TYPE TimeFunction
!~      PRIVATE
	 INTEGER   ::fun_type             !function type: 0, user defined; 1, harmonic
	 REAL(DBL) ::ts, te               ! starting and ending time
	 INTEGER   ::entries              ! number of entries
	 REAL(DBL),POINTER::time_val(:)   ! the ith time, in increasing order; the amplitude of the harmonic
	 REAL(DBL),POINTER::fun_val(:)    ! the ith functional value; the period of the harmonic
	 REAL(DBL),POINTER::phase_val(:)  ! the phase value of the harmonic
END TYPE TimeFunction

CONTAINS

!*********************************************************
!*                                                       *
!* get the time function value for any arbitrary time    *
!* from a piecewise linear function or harmonic function *
!* decide where t is located, then interpret the value   *
!*                                                       *
!*********************************************************
FUNCTION GetTimeFunction(tf,t) RESULT(res)

TYPE (TimeFunction),INTENT(IN)::tf  
REAL(DBL),INTENT(IN)::t 
REAL(DBL):: res
REAL(DBL)::t1, t2 ! t1<=t<=t2
REAL(DBL)::f1,f2   ! function values corresponding to t1, t2, respectively.
	
INTEGER:: n
INTEGER:: i
REAL(DBL):: tt ! a temp variable for the current time t, if t<ts, t=ts, if t>te, t=te

IF(tf%fun_type==0)THEN ! user defined piecewise linear time function
	n=COUNT(t>tf%time_val) 

	IF(n==0) THEN 		
		res=tf%fun_val(1)   ! t is smaller than or equal to the time of first entry, and the fun equal to the first entry
	ELSE IF(n==tf%entries) THEN 
		res=tf%fun_val(tf%entries) ! t is larger than the time of last entry, and the fun equal to the last entry
	ELSE 
		t1=tf%time_val(n);   f1=tf%fun_val(n)
		t2=tf%time_val(n+1); f2=tf%fun_val(n+1)
		res=(f2-f1)/(t2-t1)*(t-t1)+f1 ! linear interplation 	
	ENDIF    
ELSE IF(tf%fun_type==1) THEN ! harmonic functions
      res=0.0d0
   IF(t<tf%ts) THEN
		tt=tf%ts
   ELSE IF(t>tf%te) THEN
		tt=tf%te
   ELSE
        tt=t
   ENDIF
   DO i=1, tf%entries
      f1=tf%time_val(i)
	  IF(tf%fun_val(i)/=0.D0) THEN
		  IF(tf%fun_val(i)==-1.0D0) THEN
			f1=f1*SIN(2.0D0*PI*tf%phase_val(i))
		  ELSE
			f1=f1*SIN(2.0D0*PI*(tt/tf%fun_val(i)+tf%phase_val(i)))
		  ENDIF
	  ENDIF
      res=res+f1
   ENDDO 
ELSE ! Other types of time functions to be implemented


ENDIF
	
END FUNCTION GetTimeFunction
!*********************************************************



!*********************************************************
!*                                                       *
!*             Initialize the time function              *
!*                                                       *
!*********************************************************

ELEMENTAL FUNCTION InitTF() RESULT(res)

TYPE (TimeFunction)::res

res%fun_type=0  
res%ts=0.d0
res%te=0.d0 
res%entries=0
NULLIFY(res%time_val,res%fun_val, res%phase_val)	 
	 	 
END FUNCTION InitTF
!*********************************************************



!************************************************************
!*                                                          *
!*  Input and echo Time Functions                           *
!*															*
!************************************************************
FUNCTION InputEchoTimeFunctions(IN,EIN,error) RESULT(res)

INTEGER,INTENT(IN)::IN,EIN ! file units for input anf echo files, respectively
CHARACTER(*),INTENT(OUT)::error
TYPE (TimeFunction)::res       
INTEGER:: i,n

READ(IN,*,IOSTAT=in_stat) res%fun_type   
IF(IOError('read function type for a time function',error)) RETURN

READ(IN,*,IOSTAT=in_stat) res%ts, res%te   
IF(IOError('read starting and ending time for a time function',error)) RETURN

READ(IN,*,IOSTAT=in_stat) n   
IF(IOError('read number of entries for a time function',error)) RETURN

res%entries=n

ALLOCATE(res%time_val(n),STAT=allo_stat)
IF(MemoryError('res%time_val',error)) GOTO 9999

ALLOCATE(res%fun_val(n),STAT=allo_stat)
IF(MemoryError('res%fun_val',error)) GOTO 9999

IF(res%fun_type==1) THEN
	ALLOCATE(res%phase_val(n),STAT=allo_stat)
	IF(MemoryError('res%phase_val',error)) GOTO 9999
ENDIF

DO i=1,   n   
	IF(res%fun_type==0) THEN
		READ(IN,*,IOSTAT=in_stat)res%time_val(i),res%fun_val(i)
	ELSEIF(res%fun_type==1)THEN
		READ(IN,*,IOSTAT=in_stat)res%time_val(i),res%fun_val(i),res%phase_val(i)
	ELSE
	    error='Time function type is not defined'
        WRITE(0,*) 'TimeFunction.f90 : ',error
        ERROR STOP 127
    ENDIF
	
	IF(IOError('read values needed to define time functions',error)) RETURN
ENDDO


WRITE(EIN,*) 
WRITE(EIN,*) 'Function Type=',     res%fun_type
WRITE(EIN,*) 'Starting/Ending Time=', res%ts, res%te
WRITE(EIN,*) 'Number of Entries=', res%entries
WRITE(EIN,*) '----------------------------------'	
	
DO i=1,res%entries
	IF( res%fun_type==0)THEN
		CALL WriteVec(EIN,(/res%time_val(i),res%fun_val(i)/))
	ELSEIF(res%fun_type==1)THEN
		CALL WriteVec(EIN,(/res%time_val(i),res%fun_val(i),res%phase_val(i)/))
	ELSE
	ENDIF
ENDDO
	
9999 IF(error/='')THEN
    WRITE(0,*) 'TimeFunction.f90 : ',error
    ERROR STOP 128
ENDIF

END FUNCTION InputEchoTimeFunctions
!************************************************************



!************************************************************
!*                                                          *
!*  Evaluate current function based on magnitude and time   *
!* function and current time                                *
!* vec: magnitude                                           *
!* vec_tf: time function #   								*
!* time_fun: array holding all time functions               *
!* time_current:  current time for evaluation               *
!************************************************************
FUNCTION CurrentValues(vec,vec_tf,time_fun,time_current)

REAL(DBL),INTENT(IN)::vec(:)
INTEGER,  INTENT(IN)::vec_tf(:)
TYPE (TimeFunction),INTENT(IN)::time_fun(:)
REAL(DBL),INTENT(IN)::time_current
REAL(DBL)           ::CurrentValues(SIZE(vec))
INTEGER:: j

CurrentValues=vec  ! initialize to be the magnitude

DO j=1, SIZE(vec)
   IF(vec_tf(j)/=0.AND.vec(j)/=0.0D0) &  
     CurrentValues(j)  =vec(j) *GetTimeFunction(time_fun(vec_tf(j)),time_current)
ENDDO

END FUNCTION CurrentValues
!************************************************************



!================================================
! Subroutines/functions used internally
!================================================



END MODULE TimeFunctionModule
!===========================================================
