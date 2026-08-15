!********************* READ INPUT ***********************************************************
      SUBROUTINE READ_INPUT()
      use declare_variables
      implicit none	  
	  
	  OPEN (finput,file='input.dat',form='formatted')
	  
	  read(finput,*) restart, nblocks
	  read(finput,*) NImax,NJmax, NKmax
	  read(finput,*) testcase, viscous
	  read(finput,*) Re, Mach, gamma, prandtl, T_ref
	  read(finput,*) nprims, nconserv
	  read(finput,*) dscheme, fscheme, alpha_f
	  read(finput,*) rk_steps, nsteps, time_step, animfreq




	!   ----------------check if this is correct----------------
	!   spoiler alert: it isnt.

		


	!   if (NKmax .eq. 1) then
	! 	grid2d = 1
	  
	!   end if	

	!   if (testcase == 1) then
	! 	taylor = 1

	!   end if

	  grid2d = 0
	  taylor = 0
	  covo   = 0

	  if (NKmax .eq. 1) then
		grid2d = 1
	  end if	

	!   if (testcase == 1) then
	! 	taylor = 1
	!   elseif (testcase == 2) then
	! 	covo = 1
	!   end if

	  if (testcase == 1) then
		taylor = 1
	  elseif (testcase == 2 .or. testcase == 3 .or. testcase == 4) then
		covo = 1
	  end if

	  Ptsmax = max(NImax,NJmax,NKmax)
	  
	!   END

	  
	  close (finput)
      
      END
!********************************************************************************************

!********************* ALLOCATE_ROUTINE *****************************************************
      SUBROUTINE ALLOCATE_ROUTINE()
      use declare_variables
      implicit none	  
	  
	  allocate(NI(nblocks))
	  allocate(NJ(nblocks))
	  allocate(NK(nblocks))
	  
	  
	  do nbl = 1,nblocks
		NI(nbl) = NIMax
		NJ(nbl) = NJMax
		NK(nbl) = NKMax
	  
	  
	  
	  Enddo

	  allocate(xgrid(NIMax,NJMax,NKMax,nblocks))
	  allocate(ygrid(NIMax,NJMax,NKMax,nblocks))
	  allocate(zgrid(NIMax,NJMax,NKMax,nblocks))

	!   allocating the x, y, j derivs wrt i,j,k for the grid, which are used in the metrics and flux calculations.
		
	  allocate(xi(NIMax,NJMax,NKMax,nblocks))
	  allocate(yi(NIMax,NJMax,NKMax,nblocks))
	  allocate(zi(NIMax,NJMax,NKMax,nblocks))

	  allocate(xj(NIMax,NJMax,NKMax,nblocks))
	  allocate(yj(NIMax,NJMax,NKMax,nblocks))
	  allocate(zj(NIMax,NJMax,NKMax,nblocks))

	  allocate(xk(NIMax,NJMax,NKMax,nblocks))
	  allocate(yk(NIMax,NJMax,NKMax,nblocks))
	  allocate(zk(NIMax,NJMax,NKMax,nblocks))


	!   allocationg ix, iy, iz, jx, jy, jz, kx, ky, kz which are the metrics of the grid used in flux calculations.

	  allocate(ix(NIMax,NJMax,NKMax,nblocks))
	  allocate(iy(NIMax,NJMax,NKMax,nblocks))
	  allocate(iz(NIMax,NJMax,NKMax,nblocks))

	  allocate(jx(NIMax,NJMax,NKMax,nblocks))
	  allocate(jy(NIMax,NJMax,NKMax,nblocks))
	  allocate(jz(NIMax,NJMax,NKMax,nblocks))

	  allocate(kx(NIMax,NJMax,NKMax,nblocks))
	  allocate(ky(NIMax,NJMax,NKMax,nblocks))
	  allocate(kz(NIMax,NJMax,NKMax,nblocks))

	!   allocating the Jacobian 

	  allocate(Jac(NIMax,NJMax,NKMax,nblocks))

	  allocate(tked(NIMAx,NJMax,NKMax,nblocks))
	  allocate(enst(NIMAx,NJMax,NKMax,nblocks))


	!   seems like we might have to deallocate these bad bois later.


	  allocate(Qp(NIMax,NJMax,NKMax,nblocks,nprims))
      allocate(Qc(NIMax,NJMax,NKMax,nblocks,nconserv))

	  


	!   allocating the flux terms F G H


	  allocate(Fflux(NIMax,NJMax,NKMax,nblocks,nconserv))
	  allocate(Hflux(NIMax,NJMax,NKMax,nblocks,nconserv))
	  allocate(Gflux(NIMax,NJMax,NKMax,nblocks,nconserv))
	  allocate(fluxD(NIMax,NJMax,NKMax,nblocks,nconserv))
	  allocate(net_flux(NIMax,NJMax,NKMax,nblocks,nconserv))



	!   for the rk4 timestepping----------------------

	!   allocate(Qcini(NIMax,NJMax,NKMax,nblocks,nprims)) !another monument to my sins

	  allocate(Qcini(NIMax,NJMax,NKMax,nblocks,nconserv))

      allocate(Qcnew(NIMax,NJMax,NKMax,nblocks,nconserv))



	!   visous terms? aaaaaaaaaaaaaaaaaaaaaaaaaa

	  if (viscous.eq.1) then


		allocate(Qpi(NIMax,NJMax,NKMax,nblocks,nprims))
		allocate(Qpj(NIMax,NJMax,NKMax,nblocks,nprims))
		allocate(Qpk(NIMax,NJMax,NKMax,nblocks,nprims))

	  endif



	  allocate(AMD(Ptsmax)) ! lower diagonal for discretization
	  allocate(ACD(Ptsmax)) ! diagonal for discretization
	  allocate(APD(Ptsmax)) ! upper diagonal for discretization

	  allocate(AMF(Ptsmax)) !filtering majiggy
	  allocate(ACF(Ptsmax))
	  allocate(APF(Ptsmax))

	  allocate(fac_RK(rk_steps)) ! RK coefficients
	  allocate(fac_qini(rk_steps)) ! RK coefficients

	  allocate(res(nconserv))
	  allocate(fcoeff(6))




	!   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa



	!   allocate(PHIi(NImax,NJmax,NKmax,nblocks,nprims))
	!   allocate(PHIj(NImax,NJmax,NKmax,nblocks,nprims))
	!   allocate(PHIk(NImax,NJmax,NKmax,nblocks,nprims))


! 	  allocate(PHIi(NImax,NJmax,NKmax,nblocks,nprims))
	! allocate(PHIj(NImax,NJmax,NKmax,nblocks,nprims))
	! allocate(PHIk(NImax,NJmax,NKmax,nblocks,nprims))
		

		


      END 
!********************************************************************************************

!********************* GENERATE_GRID *****************************************************
      SUBROUTINE GENERATE_GRID()
      use declare_variables 
      implicit none	  
	  
	  !*****For TGV Case*****
	  if (testcase.eq.1) then
		! Lx = 2.d0*22.d0/7.d0
		! Ly = 2.d0*22.d0/7.d0
		! Lz = 2.d0*22.d0/7.d0

		Lx = 2.d0*acos(-1.d0)
		Ly = 2.d0*acos(-1.d0)
		Lz = 2.d0*acos(-1.d0)

		Do nbl = 1,nblocks
		Do k = 1,NK(nbl)
		Do j = 1,NJ(nbl)
		Do i = 1,NI(nbl)
			xgrid(i,j,k,nbl) = Lx*real(i-1.d0)/real(NI(nbl)-1.d0)
			ygrid(i,j,k,nbl) = Ly*real(j-1.d0)/real(NJ(nbl)-1.d0)
			zgrid(i,j,k,nbl) = Lz*real(k-1.d0)/real(NK(nbl)-1.d0)
		Enddo
		Enddo
		Enddo
		Enddo

	elseif (testcase.eq.2) then
		Lx = 16.d0
		Ly = 16.d0
		Lz = 1.d0

		Do nbl = 1,nblocks
		Do k = 1,NK(nbl)
		Do j = 1,NJ(nbl)
		Do i = 1,NI(nbl)
			xgrid(i,j,k,nbl) = -8.d0 + Lx*real(i-1)/real(NI(nbl)-1)
			ygrid(i,j,k,nbl) = -8.d0 + Ly*real(j-1)/real(NJ(nbl)-1)
			zgrid(i,j,k,nbl) = 0.d0
		Enddo
		Enddo
		Enddo
		Enddo

	

	! dcell   = (Lx/real(NImax-1)) * (Ly/real(NJmax-1)) * (Lz/real(NKmax-1))
	! Vdomain = Lx * Ly * Lz



	elseif (testcase.eq.3) then


	! covo on randomised mesh (visbal and gaitonde 2002, fig 5)
    ! uniform grid on [-8,8] with interior points randomly
    ! perturbed by 20% of spacing because life wasnt hard enough
    ! boundary points left alone so periodic bcs dont implode


      Lx = 16.d0
      Ly = 16.d0
      Lz = 1.d0



      Do nbl = 1,nblocks
      Do k = 1,NK(nbl)
      Do j = 1,NJ(nbl)
      Do i = 1,NI(nbl)
        


        xgrid(i,j,k,nbl) = Lx*real(i-1.d0)/real(NI(nbl)-1) - Lx/2.d0
        ygrid(i,j,k,nbl) = Ly*real(j-1.d0)/real(NJ(nbl)-1) - Ly/2.d0
        zgrid(i,j,k,nbl) = 0.d0

        
        if (xgrid(i,j,k,nbl).ge.-5.d0 .and. xgrid(i,j,k,nbl).le.5.d0 .and. &
            ygrid(i,j,k,nbl).ge.-5.d0 .and. ygrid(i,j,k,nbl).le.5.d0) then
          xgrid(i,j,k,nbl) = xgrid(i,j,k,nbl) + &
            0.2d0*(Lx/(NI(nbl)-1.d0))*(1.d0-rand())/2.d0
          ygrid(i,j,k,nbl) = ygrid(i,j,k,nbl) + &
            0.2d0*(Ly/(NJ(nbl)-1.d0))*(1.d0-rand())/2.d0
        endif

      Enddo
      Enddo
      Enddo
      Enddo

    elseif (testcase.eq.4) then


    ! covo on sinusoidal mesh (visbal and gaitonde 2002, eq 21, fig 6)
    ! analytically generated mesh that looks smooth but is horribly skewed
    ! frozen at max distortion because why would we make it easy
    !
    ! x(i,j) = xmin + dx_o * [(i-1) + Ax * sin(nx*pi*(j-1)*dy_o/Ly)]
    ! y(i,j) = ymin + dy_o * [(j-1) + Ay * sin(ny*pi*(i-1)*dx_o/Lx)]
    !
    ! paper uses Ax=1, Ay=2, nx=ny=6 on [-6,6]
    ! we stretch it to [-8,8] so the vortex has room to suffer
		! why


      Lx = 16.d0
      Ly = 16.d0
      Lz = 1.d0

      block
      real :: dx_o, dy_o, Ax_amp, Ay_amp
      integer :: nx_wave, ny_wave

      dx_o = Lx / real(NI(1)-1)
      dy_o = Ly / real(NJ(1)-1)

      
      Ax_amp  = 1.0d0      ! x-distortion amplitude
      Ay_amp  = 2.0d0      ! y-distortion amplitude
      nx_wave = 6           ! number of x-direction waves
      ny_wave = 6           ! number of y-direction waves

      Do nbl = 1,nblocks
      Do k = 1,NK(nbl)
      Do j = 1,NJ(nbl)
      Do i = 1,NI(nbl)
        xgrid(i,j,k,nbl) = -Lx/2.d0 + dx_o * ( real(i-1) + &
          Ax_amp * sin(nx_wave * acos(-1.d0) * real(j-1) * dy_o / Ly) )
        ygrid(i,j,k,nbl) = -Ly/2.d0 + dy_o * ( real(j-1) + &
          Ay_amp * sin(ny_wave * acos(-1.d0) * real(i-1) * dx_o / Lx) )
        zgrid(i,j,k,nbl) = 0.d0
      Enddo
      Enddo
      Enddo
      Enddo

      end block

	  
	endif

	if (grid2d.eq.1) then
		dcell   = (Lx/real(NImax-1)) * (Ly/real(NJmax-1))
		Vdomain = Lx * Ly
	else
		dcell   = (Lx/real(NImax-1)) * (Ly/real(NJmax-1)) * (Lz/real(NKmax-1))
		Vdomain = Lx * Ly * Lz
	endif
	  
	  
      END 
	
!********************************************************************************************

! !********************* INITIALIZE_NON_DIMENSIONALIZE ****************************************
!       SUBROUTINE INITIALIZE_NON_DIMENSIONALIZE()
!       use declare_variables 
!       implicit none	  
	  
!       real xl,yl,zl, Etotal, ul, vl, wl, rhl, pl, Tl
!       integer nvars, n

! 	! print*, 'lin85'
! 	  Do nbl = 1,nblocks
! 	  Do k = 1,NK(nbl)
! 	  Do j = 1,NJ(nbl)
! 	  Do i = 1,NI(nbl)

! 		xl = xgrid(i,j,k,nbl)
! 		yl = ygrid(i,j,k,nbl)
! 		zl = zgrid(i,j,k,nbl)
		
! 		! print*, 'lin95',i,j,k,nbl





! 		! for taylor green vortex, the initial conditions are given by the following equations. 

! 		! probably incorrect. ill comment them out again.

! 		! Qp(i,j,k,nbl,2) = sin(xl)*sin(yl)*sin(zl)	!u
! 		! Qp(i,j,k,nbl,3) = -cos(xl)*sin(yl)*cos(zl)	!

! 		! ! print*, 'lin100',i,j,k,nbl
! 		! Qp(i,j,k,nbl,4) = 0.d0	!w maybe?
! 		! Qp(i,j,k,nbl,5) = 1.d0/(gamma*Mach**2) + 1.d0/16.d0*(cos(2*xl)+cos(2*yl)+cos(2*zl)+2.d0) !pres

! 		! u-velocity (Fixed sin/cos combination)
! 		Qp(i,j,k,nbl,2) = sin(xl)*cos(yl)*cos(zl)   

! 		! v-velocity (Assuming yours was already correct, it should be this)
! 		Qp(i,j,k,nbl,3) = -cos(xl)*sin(yl)*cos(zl)  

! 		! w-velocity
! 		Qp(i,j,k,nbl,4) = 0.d0                      

! 		! Pressure (Fixed the multiplication of the spatial groupings)
! 		Qp(i,j,k,nbl,5) = 1.d0/(gamma*Mach**2) + 1.d0/16.d0 * (cos(2.d0*xl) + cos(2.d0*yl)) * (cos(2.d0*zl) + 2.d0)


! 		Qp(i,j,k,nbl,6) = 1.d0	!temp
! 		Qp(i,j,k,nbl,1) = gamma*Mach**2*Qp(i,j,k,nbl,5)/Qp(i,j,k,nbl,6)  
! 		! print*, 'lin105',i,j,k,nbl


		


! 		! local variables for readability

! 		pl = Qp(i,j,k,nbl,5)
! 		Tl = Qp(i,j,k,nbl,6)
! 		rhl = Qp(i,j,k,nbl,1)

! 		ul = Qp(i,j,k,nbl,2)
! 		vl = Qp(i,j,k,nbl,3)
! 		wl = Qp(i,j,k,nbl,4)

		

! 		Etotal = Qp(i,j,k,nbl,6)/(gamma*(gamma-1.d0)*Mach**2) + 0.5d0*(ul**2+vl**2+wl**2)

! 		Qc(i,j,k,nbl,1) = rhl	!rho
! 		Qc(i,j,k,nbl,2) = rhl*ul	!rhou
! 		Qc(i,j,k,nbl,3) = rhl*vl	!rhov
! 		Qc(i,j,k,nbl,4) = rhl*wl	!rhow
! 		Qc(i,j,k,nbl,5) = rhl*Etotal	!rhoE

! 		! wisdom - 0.5d0 is written for double precision.

! 		! ------------------------------------

! 		! write the restart bit here

! 		! ---------------------------------------


		

! 	  Enddo
! 	  Enddo
! 	  Enddo
! 	  Enddo





! 	  if (testcase.eq.2) then

!     real :: xc, yc, r2, C_vor, R_vor, pinf

!     xc    = 0.d0
!     yc    = 0.d0
!     C_vor = 0.02d0
!     R_vor = 1.d0
!     pinf  = 1.d0/(gamma*Mach**2)

!     Do nbl = 1,nblocks
!     Do k = 1,NK(nbl)
!     Do j = 1,NJ(nbl)
!     Do i = 1,NI(nbl)

!         xl = xgrid(i,j,k,nbl)
!         yl = ygrid(i,j,k,nbl)

!         r2 = ((xl-xc)**2 + (yl-yc)**2) / R_vor**2

!         ul  = 1.d0 - C_vor*(yl-yc)/R_vor**2 * exp(-r2/2.d0)
!         vl  = C_vor*(xl-xc)/R_vor**2 * exp(-r2/2.d0)
!         wl  = 0.d0
!         pl  = pinf - 1.d0*C_vor**2/(2.d0*R_vor**2) * exp(-r2)
!         rhl = 1.d0
!         Tl  = gamma*Mach**2*pl/rhl

!         Qp(i,j,k,nbl,1) = rhl
!         Qp(i,j,k,nbl,2) = ul
!         Qp(i,j,k,nbl,3) = vl
!         Qp(i,j,k,nbl,4) = wl
!         Qp(i,j,k,nbl,5) = pl
!         Qp(i,j,k,nbl,6) = Tl

!         Etotal = Tl/(gamma*(gamma-1.d0)*Mach**2) + 0.5d0*(ul**2+vl**2+wl**2)

!         Qc(i,j,k,nbl,1) = rhl
!         Qc(i,j,k,nbl,2) = rhl*ul
!         Qc(i,j,k,nbl,3) = rhl*vl
!         Qc(i,j,k,nbl,4) = rhl*wl
!         Qc(i,j,k,nbl,5) = rhl*Etotal

!     Enddo
!     Enddo
!     Enddo
!     Enddo

! 	endif
	  

	  
!       END 
! !********************************************************************************************


	!   cleaned up:

!********************* INITIALIZE_NON_DIMENSIONALIZE ****************************************
      SUBROUTINE INITIALIZE_NON_DIMENSIONALIZE()
      use declare_variables
      implicit none

      real :: xl, yl, zl, Etotal, ul, vl, wl, rhl, pl, Tl
      real :: xc, yc, r2, C_vor, R_vor, pinf
      integer :: nvars, n

      if (testcase.eq.1) then

         Do nbl = 1,nblocks
         Do k = 1,NK(nbl)
         Do j = 1,NJ(nbl)
         Do i = 1,NI(nbl)

            xl = xgrid(i,j,k,nbl)
            yl = ygrid(i,j,k,nbl)
            zl = zgrid(i,j,k,nbl)

            Qp(i,j,k,nbl,2) = sin(xl)*cos(yl)*cos(zl)
            Qp(i,j,k,nbl,3) = -cos(xl)*sin(yl)*cos(zl)
            Qp(i,j,k,nbl,4) = 0.d0
            Qp(i,j,k,nbl,5) = 1.d0/(gamma*Mach**2) + 1.d0/16.d0 * (cos(2.d0*xl) + cos(2.d0*yl)) * (cos(2.d0*zl) + 2.d0)
            Qp(i,j,k,nbl,6) = 1.d0
            Qp(i,j,k,nbl,1) = gamma*Mach**2*Qp(i,j,k,nbl,5)/Qp(i,j,k,nbl,6)

            pl  = Qp(i,j,k,nbl,5)
            Tl  = Qp(i,j,k,nbl,6)
            rhl = Qp(i,j,k,nbl,1)
            ul  = Qp(i,j,k,nbl,2)
            vl  = Qp(i,j,k,nbl,3)
            wl  = Qp(i,j,k,nbl,4)

            Etotal = Tl/(gamma*(gamma-1.d0)*Mach**2) + 0.5d0*(ul**2+vl**2+wl**2)

            Qc(i,j,k,nbl,1) = rhl
            Qc(i,j,k,nbl,2) = rhl*ul
            Qc(i,j,k,nbl,3) = rhl*vl
            Qc(i,j,k,nbl,4) = rhl*wl
            Qc(i,j,k,nbl,5) = rhl*Etotal

         Enddo
         Enddo
         Enddo
         Enddo

      
	! elseif (testcase.eq.2) then
	elseif (testcase.eq.2 .or. testcase.eq.3 .or. testcase.eq.4) then

         xc    = 0.d0
         yc    = 0.d0
         C_vor = 0.02d0
         R_vor = 1.d0
         pinf  = 1.d0/(gamma*Mach**2)

         Do nbl = 1,nblocks
         Do k = 1,NK(nbl)
         Do j = 1,NJ(nbl)
         Do i = 1,NI(nbl)

            xl = xgrid(i,j,k,nbl)
            yl = ygrid(i,j,k,nbl)

            r2  = ((xl-xc)**2 + (yl-yc)**2) / R_vor**2

            ul  = 1.d0 - C_vor*(yl-yc)/R_vor**2 * exp(-r2/2.d0)
            vl  = C_vor*(xl-xc)/R_vor**2 * exp(-r2/2.d0)
            wl  = 0.d0
            pl  = pinf - C_vor**2/(2.d0*R_vor**2) * exp(-r2)
            rhl = 1.d0
            Tl  = gamma*Mach**2*pl/rhl

            Qp(i,j,k,nbl,1) = rhl
            Qp(i,j,k,nbl,2) = ul
            Qp(i,j,k,nbl,3) = vl
            Qp(i,j,k,nbl,4) = wl
            Qp(i,j,k,nbl,5) = pl
            Qp(i,j,k,nbl,6) = Tl

            Etotal = Tl/(gamma*(gamma-1.d0)*Mach**2) + 0.5d0*(ul**2+vl**2+wl**2)

            Qc(i,j,k,nbl,1) = rhl
            Qc(i,j,k,nbl,2) = rhl*ul
            Qc(i,j,k,nbl,3) = rhl*vl
            Qc(i,j,k,nbl,4) = rhl*wl
            Qc(i,j,k,nbl,5) = rhl*Etotal

         Enddo
         Enddo
         Enddo
         Enddo

      endif

      END
!********************************************************************************************

!********************* DISCRETIZATION_VALS **************************************************
SUBROUTINE DISCRETIZATION_FILTER_RK_VALS()
	use declare_variables 
	implicit none	  

	if (dscheme.eq.1) then

		alpha = 0.d0
		adisc = 1.d0
		bdisc = 0.d0

	elseif (dscheme.eq.2) then

		alpha = 0.d0
		adisc = 4.d0/3.d0
		bdisc = -1.d0/3.d0

	elseif (dscheme.eq.3) then

		alpha = 1.d0/4.d0
		adisc = 3.d0/2.d0
		bdisc = 0.d0

		


	elseif (dscheme.eq.4) then

		alpha = 1.d0/3.d0
		adisc = 14.d0/9.d0
		bdisc = 1.d0/9.d0

	
	else   ! dscheme >= 5: coefficients for EXP grid/viscous differentiation
		alpha = 0.d0
		adisc = 4.d0/3.d0
		bdisc = -1.d0/3.d0


	endif

	AMD = alpha
	ACD = 1.d0
	APD = alpha

	! -------------- RK coeffs?-----------------

	fac_qini(1) = 0.d0
	fac_qini(2) = 0.5d0
	fac_qini(3) = 0.5d0	
	fac_qini(4) = 1.d0

	fac_RK(1) = 1.d0/6.d0
	fac_RK(2) = 2.d0/6.d0
	fac_RK(3) = 2.d0/6.d0
	fac_RK(4) = 1.d0/6.d0


	! -----------------filtering schemes!!!!!----------gaitonde 2002


	! fill em up!!!!
	! ----------
	! -----------
	! -----------

	if (fscheme == 2) then

		! aaaaaaaaaaa
		fcoeff(1) = 1.d0/2.d0 + alpha_f
		fcoeff(2) = 1.d0/2.d0 + alpha_f
		fcoeff(3) = 0.d0
		fcoeff(4) = 0.d0
		fcoeff(5) = 0.d0
		fcoeff(6) = 0.d0

	endif

	if (fscheme == 4) then

		fcoeff(1) = 5.d0/8.d0 + (3.d0*alpha_f)/4.d0
		fcoeff(2) = 1.d0/2.d0 + alpha_f
		fcoeff(3) = -1.d0/8.d0 + alpha_f/4.d0
		fcoeff(4) = 0.d0
		fcoeff(5) = 0.d0
		fcoeff(6) = 0.d0



		! aaaaaaaaaaa

	endif

	if (fscheme == 6) then

		! aaaaaaaaaaa

		fcoeff(1) = 11.d0/16.d0 + (5.d0*alpha_f)/8.d0
		fcoeff(2) = 15.d0/32.d0 + (17.d0*alpha_f)/16.d0
		fcoeff(3) = -3.d0/16.d0 + (3.d0*alpha_f)/8.d0
		fcoeff(4) = 1.d0/32.d0 - alpha_f/16.d0
		fcoeff(5) = 0.d0
		fcoeff(6) = 0.d0

	endif

	if (fscheme == 8) then

		! aaaaaaaaaaa

		fcoeff(1) = (93.d0 + 70.d0*alpha_f)/128.d0
		fcoeff(2) = (7.d0 + 18.d0*alpha_f)/16.d0
		fcoeff(3) = (-7.d0 + 14.d0*alpha_f)/32.d0
		fcoeff(4) = 1.d0/16.d0 - alpha_f/8.d0
		fcoeff(5) = -1.d0/128.d0 + alpha_f/64.d0
		fcoeff(6) = 0.d0

	endif

	if (fscheme == 10) then

		! aaaaaaaaaaa

		fcoeff(1) = (193.d0 + 126.d0*alpha_f)/256.d0
		fcoeff(2) = (105.d0 + 302.d0*alpha_f)/256.d0
		fcoeff(3) = 15.d0*(-1.d0 + 2.d0*alpha_f)/64.d0
		fcoeff(4) = 45.d0*(1.d0 - 2.d0*alpha_f)/512.d0
		fcoeff(5) = 5.d0*(-1.d0 + 2.d0*alpha_f)/256.d0
		fcoeff(6) = (1.d0 - 2.d0*alpha_f)/512.d0

	endif

	AMF = alpha_f
	ACF = 1.d0
	APF = alpha_f

	



	! -----------------------------------------------------------------



		
END SUBROUTINE
!********************************************************************************************

!********************* METRICS **************************************************************
	! what does this bad boi do? grid derivs

	SUBROUTINE METRICS()
	use declare_variables 
	implicit none	  
	
	real xil,yil,zil,xjl,yjl,zjl,xkl,ykl,zkl	  
	real vol

	! if (dscheme.eq.1 .or. dscheme.eq.2) then

		
	! 	call DISCRETIZATION_I_EXP_GRID(xgrid,xi)
	! 	call DISCRETIZATION_I_EXP_GRID(ygrid,yi)
	! 	call DISCRETIZATION_I_EXP_GRID(zgrid,zi)


	! 	! check these once:

	! 	call DISCRETIZATION_J_EXP_GRID(xgrid,xj)
	! 	call DISCRETIZATION_J_EXP_GRID(ygrid,yj)
	! 	call DISCRETIZATION_J_EXP_GRID(zgrid,zj)

	! 	call DISCRETIZATION_K_EXP_GRID(xgrid,xk)
	! 	call DISCRETIZATION_K_EXP_GRID(ygrid,yk)
	! 	call DISCRETIZATION_K_EXP_GRID(zgrid,zk)

	! elseif (dscheme.eq.3 .or. dscheme.eq.4) then

	! 	call DISCRETIZATION_I_COMP_GRID(xgrid,xi)
	! 	call DISCRETIZATION_I_COMP_GRID(ygrid,yi)
	! 	call DISCRETIZATION_I_COMP_GRID(zgrid,zi)

	! 	call DISCRETIZATION_J_COMP_GRID(xgrid,xj)
	! 	call DISCRETIZATION_J_COMP_GRID(ygrid,yj)
	! 	call DISCRETIZATION_J_COMP_GRID(zgrid,zj)

	! 	call DISCRETIZATION_K_COMP_GRID(xgrid,xk)
	! 	call DISCRETIZATION_K_COMP_GRID(ygrid,yk)
	! 	call DISCRETIZATION_K_COMP_GRID(zgrid,zk)

	! endif


	if (dscheme.eq.1 .or. dscheme.eq.2) then

		call DISCRETIZATION_I_EXP_GRID(xgrid,xi)
		call DISCRETIZATION_I_EXP_GRID(ygrid,yi)
		call DISCRETIZATION_I_EXP_GRID(zgrid,zi)

		call DISCRETIZATION_J_EXP_GRID(xgrid,xj)
		call DISCRETIZATION_J_EXP_GRID(ygrid,yj)
		call DISCRETIZATION_J_EXP_GRID(zgrid,zj)

		if (grid2d.ne.1) then
			call DISCRETIZATION_K_EXP_GRID(xgrid,xk)
			call DISCRETIZATION_K_EXP_GRID(ygrid,yk)
			call DISCRETIZATION_K_EXP_GRID(zgrid,zk)
		endif

	elseif (dscheme.eq.3 .or. dscheme.eq.4) then

		call DISCRETIZATION_I_COMP_GRID(xgrid,xi)
		call DISCRETIZATION_I_COMP_GRID(ygrid,yi)
		call DISCRETIZATION_I_COMP_GRID(zgrid,zi)

		call DISCRETIZATION_J_COMP_GRID(xgrid,xj)
		call DISCRETIZATION_J_COMP_GRID(ygrid,yj)
		call DISCRETIZATION_J_COMP_GRID(zgrid,zj)

		if (grid2d.ne.1) then
			call DISCRETIZATION_K_COMP_GRID(xgrid,xk)
			call DISCRETIZATION_K_COMP_GRID(ygrid,yk)
			call DISCRETIZATION_K_COMP_GRID(zgrid,zk)
		endif

	else   ! dscheme >= 5 (WENO): use explicit central for grid derivatives

		call DISCRETIZATION_I_EXP_GRID(xgrid,xi)
		call DISCRETIZATION_I_EXP_GRID(ygrid,yi)
		call DISCRETIZATION_I_EXP_GRID(zgrid,zi)

		call DISCRETIZATION_J_EXP_GRID(xgrid,xj)
		call DISCRETIZATION_J_EXP_GRID(ygrid,yj)
		call DISCRETIZATION_J_EXP_GRID(zgrid,zj)

		if (grid2d.ne.1) then
			call DISCRETIZATION_K_EXP_GRID(xgrid,xk)
			call DISCRETIZATION_K_EXP_GRID(ygrid,yk)
			call DISCRETIZATION_K_EXP_GRID(zgrid,zk)
		endif

	endif

	! 2D case, set K-direction grid derivatives manually
	if (grid2d.eq.1) then
		xk = 0.d0
		yk = 0.d0
		zk = 1.d0   ! dz/dk = 1 (unit spacing in computational space)
	endif



	Do nbl = 1,nblocks
		Do k = 1,NK(nbl)
			Do j = 1,NJ(nbl)
				Do i = 1,NI(nbl)

					xil = xi(i,j,k,nbl)
					yil = yi(i,j,k,nbl)	
					zil = zi(i,j,k,nbl)

					xjl = xj(i,j,k,nbl)
					yjl = yj(i,j,k,nbl)
					zjl = zj(i,j,k,nbl)

					xkl = xk(i,j,k,nbl)
					ykl = yk(i,j,k,nbl)	
					zkl = zk(i,j,k,nbl)

					Vol = xil*(yjl*zkl-zjl*ykl) - xjl*(yil*zkl-zil*ykl) + xkl*(yil*zjl-zil*yjl) 	!the inverse of the jacobian. refer to profs slides
					Jac(i,j,k,nbl) = 1.d0/Vol


					! -------------compute the metrics----------------


					! -------------------------------------------------

					! this is generated by copilot and the indices are fucked. refer to the slides and rectify them!!!!!!!!!!!!!!!

					! --------------------------------------------------

					! ! incorrect indices! refer to the slides and rectify them!!!!!!!!!!!!!!!

					! ix(i,j,k,nbl) = (yj(i,j,k,nbl)*zk(i,j,k,nbl) - zj(i,j,k,nbl)*yk(i,j,k,nbl))*Jac(i,j,k,nbl)
					! iy(i,j,k,nbl) = (zi(i,j,k,nbl)*xk(i,j,k,nbl) - xi(i,j,k,nbl)*zk(i,j,k,nbl))*Jac(i,j,k,nbl)
					! iz(i,j,k,nbl) = (xi(i,j,k,nbl)*yk(i,j,k,nbl) - yi(i,j,k,nbl)*xk(i,j,k,nbl))*Jac(i,j,k,nbl)


					! inside the i,j,k,nbl loops, after Jac(i,j,k,nbl) is set
					ix(i,j,k,nbl) = (yj(i,j,k,nbl)*zk(i,j,k,nbl) - yk(i,j,k,nbl)*zj(i,j,k,nbl)) * Jac(i,j,k,nbl)
					jx(i,j,k,nbl) = (yk(i,j,k,nbl)*zi(i,j,k,nbl) - yi(i,j,k,nbl)*zk(i,j,k,nbl)) * Jac(i,j,k,nbl)
					kx(i,j,k,nbl) = (yi(i,j,k,nbl)*zj(i,j,k,nbl) - yj(i,j,k,nbl)*zi(i,j,k,nbl)) * Jac(i,j,k,nbl)

					iy(i,j,k,nbl) = (xk(i,j,k,nbl)*zj(i,j,k,nbl) - xj(i,j,k,nbl)*zk(i,j,k,nbl)) * Jac(i,j,k,nbl)
					jy(i,j,k,nbl) = (xi(i,j,k,nbl)*zk(i,j,k,nbl) - xk(i,j,k,nbl)*zi(i,j,k,nbl)) * Jac(i,j,k,nbl)
					ky(i,j,k,nbl) = (xj(i,j,k,nbl)*zi(i,j,k,nbl) - xi(i,j,k,nbl)*zj(i,j,k,nbl)) * Jac(i,j,k,nbl)

					iz(i,j,k,nbl) = (xj(i,j,k,nbl)*yk(i,j,k,nbl) - xk(i,j,k,nbl)*yj(i,j,k,nbl)) * Jac(i,j,k,nbl)
					jz(i,j,k,nbl) = (xk(i,j,k,nbl)*yi(i,j,k,nbl) - xi(i,j,k,nbl)*yk(i,j,k,nbl)) * Jac(i,j,k,nbl)
					kz(i,j,k,nbl) = (xi(i,j,k,nbl)*yj(i,j,k,nbl) - xj(i,j,k,nbl)*yi(i,j,k,nbl)) * Jac(i,j,k,nbl)

					! print*, ix(i,j,k,nbl), iy(i,j,k,nbl), iz(i,j,k,nbl)

					


				enddo
			enddo

			
		enddo

		

	enddo

	! see if this can be done

	! deallocate (xi, yi, zi, xj, yj, zj, xk, yk, zk) ! deallocating the derivatives of the grid since we have calculated the metrics and we dont need them anymore. this is also more memory efficient.

	! print*, ix(i,j,k,nbl), iy(i,j,k,nbl), iz(i,j,k,nbl)



	


	
	  
	END 
!********************************************************************************************