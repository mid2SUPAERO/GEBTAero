!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************

!===========================================================
! A module use to calculate the computation time
! USE:
!
! USE CPUTime
! 
! CALL TIC  ! start clock
!...........................
! The code needed to evaluate the computing time
!.....................
! computing_time=TOC() ! ending. 
!===========================================================


!>A module use to calculate the computation time
MODULE CPUTime

IMPLICIT NONE

PRIVATE                              ! So everything is private
PUBLIC TIC, TOC                      ! execpt ...

INTEGER(8)::start, rate, finish
!=====================================================

CONTAINS
!==================================================

!***************************************
!*                                     *
!*    Starting the system clock        *
!*                                     *
!***************************************
!> Start the timer
SUBROUTINE TIC

CALL SYSTEM_CLOCK(start,rate)

END SUBROUTINE TIC
!***************************************


!****************************************************
!*                                                  *
!*   Ending the system clock and calculate time     *
!*                                                  *
!****************************************************
!> Stop the timer
FUNCTION TOC() RESULT(sec)

REAL::sec !<computation time (seconds)

CALL SYSTEM_CLOCK(finish)
IF(finish>start) THEN 
	sec=REAL(finish-start)/REAL(rate)
ELSE 
    sec=0.0
ENDIF 

END FUNCTION TOC
!****************************************************


END MODULE CPUTime
!==========================================
