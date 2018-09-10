!***********************************************************
!*                          GEBT                           *
!*                                                         *
!*          --Geometrically Exact Beam Theory ---          *
!*                                                         *
!*=========================================================*
!*                                                         *
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!*=========================================================*
!*                                                         *
!*         Created on      : June 26, 2008                 *
!*         Modified on     : July 15, 2009                 *
!*         Modified on     : July 6,  2010                 *
!*=========================================================*
!*                                                         *
!*                OBJECTIVE OF GEBT                        *
!*            =============================                *
!*                                                         * 
!* This code implements the mixed-formulation of intrinsic *
!* formulation of the geometrically exact beam theory      *
!* developed by Prof. Hodges. It is fully compatible with  *
!* VABS and it is designed to work together with VABS      *
!*                                                         * 
!*                      Alpha Version                      *
!*            =============================                *
!*                                                         * 
!* This version has the capability for linear static       *
!* analysis, which was finished during a summer faculty    *
!* fellowship at WPAFB                                     *
!*                                                         *
!*                      Beta Version                       *
!*            =============================                *
!*                                                         * 
!* This version uses the DLL concept. Seeks a steady state *  
!* solution of the linear theory.                          *
!*                                                         *
!*                                                         *
!*                      1.0                                *
!*            =============================                *
!*                                                         * 
!* This version enables the nonlinearity capability of GEBT*  
!*                                                         *
!*                      1.5                                *
!*            =============================                *
!*                                                         * 
!* This version enables the follower force capability      *  
!*                                                         *
!*                                                         *
!*                      2.0                                *
!*            =============================                *
!*                                                         * 
!* This version enables AD capability using DNAD           *  
!*                                                         *
!*                                                         *
!*                      3.0                                *
!*            =============================                *
!*                                                         * 
!* This version enables dynamic response                   *  
!*                                                         *
!*                                                         *
!*                      4.0                                *
!*            =============================                *
!*                                                         * 
!* New assembly without storing the NxN matrix             *
!*                                                         *
!*                                                         *
!*                                                         *
!*                      4.2                                *
!*            =============================                *
!*                                                         * 
!* Use MA28 for linea solver and Arpack for eigen solver   *
!* Use long names of input file, including spaces in the   *
!* path and file names  (7/18/2011)                        *
!* input stiffness matrix (08/07/2011)                     *
!* fix a bug related with input stiffness matrix (10/25/2011)*
!*                                                         *
!*                                                         *
!*                                                         *
!***********************************************************

!> This program launch the aeroelastic simulation by reading the .dat file (and possibly the .ini file), execute the analysis subroutine and output the result in .out test file or/and vtk files (readable in paraview)

PROGRAM GEBT

!-----------------------------------------------------------------
! This serves as a general-purpose driving program
! Goto is used here for the purpose to terminate the execution
! due to critical errors.
!------------------------------------------------------------------
USE CPUTime
USE IOaero
USE PrescribedCondition
USE GlobalDataFun
 
IMPLICIT NONE ! Cancelling the naming convention of Fortran

CALL TIC ! start the clock

CALL Input ! Read inputs for the beam analysis: in IOaero

CALL Analysis(nkp,nelem,ndof_el,nmemb,ncond_pt,nmate, nframe,ndistrfun,ncurv,coord, &
		  &  member,pt_condition,material,aerodyn_coef,niter,nstep,sol_pt,sol_mb,error,ncond_mb, &
		  & ntimefun,frame,mb_condition,distr_fun,curvature,omega_a0,omega_a_tf, &
		  & v_root_a0,v_root_a_tf,simu_time,time_function,analysis_flag,init_cond, &
		  &  nev,eigen_val,eigen_vec_pt,eigen_vec_mb,aero_flag,grav_flag)

CALL Output ! Report results: in .out file

IF (nvtk>0) CALL OutputVtk ! Creating a vtk file for each istep containing sol_mb information

IF(ALLOCATED(eigen_vec_mb))  DEALLOCATE(eigen_vec_mb)
IF(ALLOCATED(eigen_vec_pt))  DEALLOCATE(eigen_vec_pt)		
IF(ALLOCATED(eigen_val)) 	 DEALLOCATE(eigen_val)
IF(ALLOCATED(sol_mb))		 DEALLOCATE(sol_mb)
IF(ALLOCATED(sol_pt))	     DEALLOCATE(sol_pt)
IF(ALLOCATED(init_cond))	 DEALLOCATE(init_cond)
IF(ALLOCATED(time_function)) DEALLOCATE(time_function)
IF(ALLOCATED(curvature))	 DEALLOCATE(curvature)
IF(ALLOCATED(distr_fun))     DEALLOCATE(distr_fun)
IF(ALLOCATED(mb_condition))  DEALLOCATE(mb_condition)
IF(ALLOCATED(frame)) 		 DEALLOCATE(frame)
IF(ALLOCATED(material))		 DEALLOCATE(material)
IF(ALLOCATED(aerodyn_coef))	 DEALLOCATE(aerodyn_coef)
IF(ALLOCATED(pt_condition))  DEALLOCATE(pt_condition)
IF(ALLOCATED(member))		 DEALLOCATE(member)
IF(ALLOCATED(coord)) 		 DEALLOCATE(coord)

IF (RUNMOD == 0) THEN
	WRITE(*,*)'GEBT finished successfully'

!~ 	! Determine the time spend by the system
!~ 	!===========================================
	WRITE(*,*)
	WRITE(*,*)'Program runs  ', TOC(), 'seconds'
ENDIF

IF (RUNMOD ==0) CALL WriteError(EIN,error) ! Write error to the echo file: in GEBTIO
CLOSE(EIN)

END PROGRAM GEBT
!**************************************************************


