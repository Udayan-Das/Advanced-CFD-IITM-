!***************************************************************************************
!****************** WRITTEN BY DR. NAGABHUSHANA RAO VADLAMANI **************************
!***************BASED ON THE HIGH ORDER COMPSQUARE SOLVER DEVELOPED BY DR. NRV**********
!***************DISTRIBUTED AS A PART OF AS6041 COURSE ON ADVANCED CFD *****************
!***************************************************************************************

      PROGRAM MISSION_AS6041
	  	  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
!@@@@@@@@@@@@@ PRE PROCESSING STEPS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	

            
!****** Declare Arrays, Integers, Real numbers *********************************	
      use declare_variables
      implicit none
	  
      integer step,var   
      character(30) :: monitorfile
	  
      print*, 'Declared Variables'

!******** READ INPUT ***********************************************************
      print*, 'Reading input file....'
      call READ_INPUT()
      print*, 'Reading input file Done'		  
!****** INITIALIZE ARRAYS ******************************************************
      call ALLOCATE_ROUTINE()
      res = 0.d0
      print*, 'Allocated Integers, Real numbers and Arrays Done'	
!****** GENERATE GRID ******************************************************
      call GENERATE_GRID()
      print*, 'Grid Generation Done'	  
!******** Initialize and Non-dimensionalize arrays *****************************
      print*, 'Initializing and non-dimensionalizing arrays....'
      call INITIALIZE_NON_DIMENSIONALIZE()
      print*, 'Initializing and non-dimensionalizing arrays done'
!******** Read Discretization coeffs for flow !*********************************
      print*, 'Obtaining Discretization Filter RK coefficients for flow....'
      call DISCRETIZATION_FILTER_RK_VALS()
      print*, 'Obtaining Discretization Filter RK coefficients for flow Done'           
!******** Compute Metric terms *************************************************
      print*, 'Computing Metrics....'
      call METRICS()
      print*, 'Computing Metrics Done'
!****************** Monitor residuals ******************************************

      

      if (testcase.eq.1) then
            if (dscheme.eq.1) monitorfile = 'Monitor_E2.out'
            if (dscheme.eq.2) monitorfile = 'Monitor_E4.out'
            if (dscheme.eq.3) monitorfile = 'Monitor_C4.out'
            if (dscheme.eq.4) monitorfile = 'Monitor_C6.out'
            if (dscheme.eq.5) monitorfile = 'Monitor_WENO_JS.out'
            if (dscheme.eq.6) monitorfile = 'Monitor_WENO_Z.out'
            if (dscheme.eq.7) monitorfile = 'Monitor_WENO_CU6.out'
            if (dscheme.eq.8) monitorfile = 'Monitor_WENO_CU6M.out'
      elseif (testcase.eq.2) then
            write(monitorfile,'(a,i3.3,a)') 'Monitor_COVO_N',NImax,'.out'

      elseif (testcase.eq.3) then
            write(monitorfile,'(a,i3.3,a)') 'Monitor_COVO_RAND_N',NImax,'.out'
      elseif (testcase.eq.4) then
            write(monitorfile,'(a,i3.3,a)') 'Monitor_COVO_SINE_N',NImax,'.out'

      endif

      if(restart.eq.1) OPEN(fresidual,file=monitorfile,form='formatted',access='append')
      if(restart.eq.0) OPEN(fresidual,file=monitorfile,form='formatted')



      ! if(restart.eq.1) OPEN(fresidual,file='Monitor.out',form='formatted',access='append')	  	  
      ! if(restart.eq.0) OPEN(fresidual,file='Monitor.out',form='formatted')	
      
      


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!@@@@@@@@@@@@@@@@@@@@@@@ MAIN TIME LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      time = 0.d0

      
      DO iter = 1,nsteps   !EXPLICIT TIME STEPPING
      
	      Qcini = Qc         !Conservative variables before entering RK time integration
	      Qcnew = Qc         !Conservative variable for the new time step, initialized to Old value
		
       Do step = 1,rk_steps
       
        !Estimate Inviscid and Viscous fluxes + Time integration
            ! print*, "line 58"
        call UNSTEADY(step)

      !   print*, "line 61"i

        
        !Filter the solution in all three directions for stability        
        if(step.eq.4) then
         call FILTERING_I(Qc,nconserv)
      !    print*, "line 66"		
         call FILTERING_J(Qc,nconserv)	
      !    print*, "line 68"	
         if(grid2d.ne.1) call FILTERING_K(Qc,nconserv)	
      !    print*, "line 70"			
        endif

      !   print*, "line 70"
        
        !Estimate Primitive variables from Conservative variables
        call SET_PRIMITIVES()

      !   print*, "line75"
        
       Enddo
		
        !Any additional post processing specific to test case
        if(taylor.eq.1) call volume_integral()


      !   print*, time, (res(var),var = 1,nconserv), tke, enstpt !tke - turbulent ke. enstpt - enstrophy 

      !   write(fresidual,*) time, tke, enstpt

      !   if (mod(iter,10).eq.0) print*, time, (res(var),var=1,nconserv), tke, enstpt
      !   if (mod(iter,10).eq.0) write(fresidual,*) time, tke, enstpt

        if (mod(iter,10).eq.0) print*, time, (res(var),var=1,nconserv), tke, enstpt
        if (mod(iter,10).eq.0 .and. testcase.eq.1) write(fresidual,*) time, tke, enstpt
      !   if (testcase.eq.2 .and. mod(iter,10).eq.0) write(fresidual,*) time
        if (covo.eq.1 .and. mod(iter,10).eq.0) write(fresidual,*) time

		res = 0.d0 !Reinitializing


      

        !Animation		
        if(mod(iter,animfreq).eq.0) call OUTPUT(1)		


            
        
        !Time increment
		time = time + time_step		
		
      ENDDO	

      ! if (testcase.eq.2) call covo_error()
      if (covo.eq.1) call covo_error()

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!@@@@@@@@@@@@@@@@@@@@@@@ END OF MAIN TIME LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      CLOSE(fresidual)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!@@@@@@@@@@@@@@@@@@@@@@@ POST PROCESSING STEPS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
	!****** Write the output ***************
101	  print*, 'Writing Output....'
      call SET_PRIMITIVES()
      call OUTPUT(0)
      print*, 'Writing output done'
 
	!****** Deallocate Arrays **************
      call DEALLOCATE_ROUTINE()
      print*, 'Deallocated Arrays'	  

      END PROGRAM

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!@@@@@@@@@@@@@@@@@@@@@@@ SUBROUTINES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!********* Subroutines in Pre-processing, Solver, Postprocessing.f90 files ****************