!********************* OUTPUT GRID & FLOW *****************************************************
      SUBROUTINE OUTPUT(flag)	  
      use declare_variables
      implicit none	  
	  !**************************** Writing Grid file***********************************************
      integer n,nvars,flag
	  character(20)::filename

	  
	  if(flag.eq.1) write(filename,'(a,i5.5,a)') 'flow',iter,'.xyz' ! For animation flag = 1
	  if(flag.eq.0) filename = 'flow.xyz'
	  
	  print*, flag, filename
	  
	  open(fgrid, file= 'grid.xyz', form = 'unformatted')
	  
	  
	  write(fgrid) nblocks
	  print*, nblocks, (NI(nbl),NJ(nbl),NK(nbl), nbl= 1, nblocks)
	  write(fgrid) (NI(nbl),NJ(nbl),NK(nbl), nbl= 1, nblocks)
	  do nbl = 1, nblocks
			write(fgrid)(((xgrid(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))&
			,           (((ygrid(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))&
			,			(((zgrid(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))
	  enddo
	  close(fgrid)
	  !**************************** Writing Flow file************************************************
	  
	  !  open(fflow, file= 'filename.xyz', form = 'unformatted')
     open(fflow, file= filename, form = 'unformatted')
	   nvars = nprims + nconserv + 9
	  
	   write(fflow) nblocks
	   write(fflow) (NI(nbl),NJ(nbl),NK(nbl),nvars, nbl= 1, nblocks)
	   do nbl = 1, nblocks
			write(fflow)((((Qp(i,j,k,nbl,n), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl)),n =1, nprims) & 
			,((((Qc(i,j,k,nbl,n), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl)),n =1, nconserv)&
			, (((ix(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))&
			, (((iy(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))&
			, (((iz(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))&
			, (((jx(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))&
			, (((jy(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))&
			, (((jz(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))&
			, (((kx(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))&
			, (((ky(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))&
			, (((kz(i,j,k,nbl), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl))
			
	   
	              
	   enddo
	   close(fflow)
	  

!***************Restart File************************!
	   !open(fflow, file= 'flow.xyz', form = 'unformatted') check for this
	   !write(fflow) nblocks
	   
	   
	  
      END SUBROUTINE
	  
      SUBROUTINE DEALLOCATE_ROUTINE()
      use declare_variables
      implicit none	  
	  
	  
      END	

      ! SUBROUTINE volume_integral()
      ! use declare_variables
      ! implicit none	  
      ! real,dimension(NJmax,NKmax,nblocks,2) :: plane
      ! real,dimension(NKmax,nblocks,2) :: line	  
      ! real muavg, dcell

	  
      ! END SUBROUTINE




! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

! SUBROUTINE volume_integral()
! use declare_variables
! implicit none

! ! real,allocatable :: PHIi(:,:,:,:,:), PHIj(:,:,:,:,:), PHIk(:,:,:,:,:)
! real :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
! real :: omx, omy, omz
! real :: rhl, ul, vl, wl
! real :: ixl, iyl, izl, jxl, jyl, jzl, kxl, kyl, kzl
! real :: dcell, Vdomain
! real :: tke_sum, enst_sum
! integer :: nvp



! nvp = nprims  ! we differentiate all 6 primitive vars but only need 2,3,4 (u,v,w)



! ! allocate(PHIi(NImax,NJmax,NKmax,nblocks,nvp))
! ! allocate(PHIj(NImax,NJmax,NKmax,nblocks,nvp))
! ! allocate(PHIk(NImax,NJmax,NKmax,nblocks,nvp))




! ! --- get computational-space derivatives of primitive variables ---
! if (dscheme.eq.1 .or. dscheme.eq.2) then
!    call DISCRETIZATION_I_EXP(Qp, PHIi, nvp)
!    call DISCRETIZATION_J_EXP(Qp, PHIj, nvp)
!    call DISCRETIZATION_K_EXP(Qp, PHIk, nvp)
! else

!    call DISCRETIZATION_I_COMP(Qp, PHIi, nvp)
!    call DISCRETIZATION_J_COMP(Qp, PHIj, nvp)
!    call DISCRETIZATION_K_COMP(Qp, PHIk, nvp)
! endif




! ! cell volume for uniform Cartesian grid (same everywhere)
! dcell = (Lx/real(NImax-1)) * (Ly/real(NJmax-1)) * (Lz/real(NKmax-1))
! Vdomain = Lx * Ly * Lz

! tke_sum  = 0.d0
! enst_sum = 0.d0

! do nbl = 1,nblocks
!    do k = 1,NK(nbl)
!       do j = 1,NJ(nbl)
!          do i = 1,NI(nbl)

!             rhl = Qp(i,j,k,nbl,1)
!             ul  = Qp(i,j,k,nbl,2)
!             vl  = Qp(i,j,k,nbl,3)
!             wl  = Qp(i,j,k,nbl,4)

!             ixl = ix(i,j,k,nbl);  iyl = iy(i,j,k,nbl);  izl = iz(i,j,k,nbl)
!             jxl = jx(i,j,k,nbl);  jyl = jy(i,j,k,nbl);  jzl = jz(i,j,k,nbl)
!             kxl = kx(i,j,k,nbl);  kyl = ky(i,j,k,nbl);  kzl = kz(i,j,k,nbl)

!             ! Cartesian velocity gradients from computational derivatives + metrics
!             u_x = PHIi(i,j,k,nbl,2)*ixl + PHIj(i,j,k,nbl,2)*jxl + PHIk(i,j,k,nbl,2)*kxl
!             u_y = PHIi(i,j,k,nbl,2)*iyl + PHIj(i,j,k,nbl,2)*jyl + PHIk(i,j,k,nbl,2)*kyl
!             u_z = PHIi(i,j,k,nbl,2)*izl + PHIj(i,j,k,nbl,2)*jzl + PHIk(i,j,k,nbl,2)*kzl

!             v_x = PHIi(i,j,k,nbl,3)*ixl + PHIj(i,j,k,nbl,3)*jxl + PHIk(i,j,k,nbl,3)*kxl
!             v_y = PHIi(i,j,k,nbl,3)*iyl + PHIj(i,j,k,nbl,3)*jyl + PHIk(i,j,k,nbl,3)*kyl
!             v_z = PHIi(i,j,k,nbl,3)*izl + PHIj(i,j,k,nbl,3)*jzl + PHIk(i,j,k,nbl,3)*kzl

!             w_x = PHIi(i,j,k,nbl,4)*ixl + PHIj(i,j,k,nbl,4)*jxl + PHIk(i,j,k,nbl,4)*kxl
!             w_y = PHIi(i,j,k,nbl,4)*iyl + PHIj(i,j,k,nbl,4)*jyl + PHIk(i,j,k,nbl,4)*kyl
!             w_z = PHIi(i,j,k,nbl,4)*izl + PHIj(i,j,k,nbl,4)*jzl + PHIk(i,j,k,nbl,4)*kzl

!             ! vorticity components
!             omx = w_y - v_z
!             omy = u_z - w_x
!             omz = v_x - u_y

!             ! accumulate — multiply by dcell for the integral
!             tke_sum  = tke_sum  + 0.5d0*rhl*(ul**2 + vl**2 + wl**2) * dcell
!             enst_sum = enst_sum + 0.5d0*rhl*(omx**2 + omy**2 + omz**2) * dcell

!          enddo
!       enddo
!    enddo
! enddo

! ! normalize by domain volume (rho0 = 1 in non-dim)
! tke   = tke_sum  / Vdomain
! enstpt = enst_sum / Vdomain

! ! deallocate(PHIi, PHIj, PHIk)








! END SUBROUTINE




SUBROUTINE volume_integral()
use declare_variables
implicit none

real :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
real :: omx, omy, omz
real :: rhl, ul, vl, wl
real :: ixl, iyl, izl, jxl, jyl, jzl, kxl, kyl, kzl
real :: tke_sum, enst_sum

tke_sum  = 0.d0
enst_sum = 0.d0

if (viscous.eq.1) then

   do nbl = 1,nblocks
      do k = 1,NK(nbl)-1
         do j = 1,NJ(nbl)-1
            do i = 1,NI(nbl)-1

               rhl = Qp(i,j,k,nbl,1)
               ul  = Qp(i,j,k,nbl,2)
               vl  = Qp(i,j,k,nbl,3)
               wl  = Qp(i,j,k,nbl,4)

               ixl = ix(i,j,k,nbl);  iyl = iy(i,j,k,nbl);  izl = iz(i,j,k,nbl)
               jxl = jx(i,j,k,nbl);  jyl = jy(i,j,k,nbl);  jzl = jz(i,j,k,nbl)
               kxl = kx(i,j,k,nbl);  kyl = ky(i,j,k,nbl);  kzl = kz(i,j,k,nbl)

               u_x = Qpi(i,j,k,nbl,2)*ixl + Qpj(i,j,k,nbl,2)*jxl + Qpk(i,j,k,nbl,2)*kxl
               u_y = Qpi(i,j,k,nbl,2)*iyl + Qpj(i,j,k,nbl,2)*jyl + Qpk(i,j,k,nbl,2)*kyl
               u_z = Qpi(i,j,k,nbl,2)*izl + Qpj(i,j,k,nbl,2)*jzl + Qpk(i,j,k,nbl,2)*kzl

               v_x = Qpi(i,j,k,nbl,3)*ixl + Qpj(i,j,k,nbl,3)*jxl + Qpk(i,j,k,nbl,3)*kxl
               v_y = Qpi(i,j,k,nbl,3)*iyl + Qpj(i,j,k,nbl,3)*jyl + Qpk(i,j,k,nbl,3)*kyl
               v_z = Qpi(i,j,k,nbl,3)*izl + Qpj(i,j,k,nbl,3)*jzl + Qpk(i,j,k,nbl,3)*kzl

               w_x = Qpi(i,j,k,nbl,4)*ixl + Qpj(i,j,k,nbl,4)*jxl + Qpk(i,j,k,nbl,4)*kxl
               w_y = Qpi(i,j,k,nbl,4)*iyl + Qpj(i,j,k,nbl,4)*jyl + Qpk(i,j,k,nbl,4)*kyl
               w_z = Qpi(i,j,k,nbl,4)*izl + Qpj(i,j,k,nbl,4)*jzl + Qpk(i,j,k,nbl,4)*kzl

               omx = w_y - v_z
               omy = u_z - w_x
               omz = v_x - u_y

               tke_sum  = tke_sum  + 0.5d0*rhl*(ul**2 + vl**2 + wl**2) * dcell
               enst_sum = enst_sum + 0.5d0*rhl*(omx**2 + omy**2 + omz**2) * dcell

            enddo
         enddo
      enddo
   enddo

else

   do nbl = 1,nblocks
      do k = 1,NK(nbl)-1
         do j = 1,NJ(nbl)-1
            do i = 1,NI(nbl)-1



               ! THIS IS INCORRECT! DO NOT RUN INVISCID AND CALL THIS SUBROUT.






               rhl = Qp(i,j,k,nbl,1)
               ul  = Qp(i,j,k,nbl,2)
               vl  = Qp(i,j,k,nbl,3)
               wl  = Qp(i,j,k,nbl,4)

               ixl = ix(i,j,k,nbl);  iyl = iy(i,j,k,nbl);  izl = iz(i,j,k,nbl)
               jxl = jx(i,j,k,nbl);  jyl = jy(i,j,k,nbl);  jzl = jz(i,j,k,nbl)
               kxl = kx(i,j,k,nbl);  kyl = ky(i,j,k,nbl);  kzl = kz(i,j,k,nbl)

               u_x = (Qp(min(i+1,NI(nbl)),j,k,nbl,2) - Qp(max(i-1,1),j,k,nbl,2)) * ixl * 0.5d0
               u_y = (Qp(i,min(j+1,NJ(nbl)),k,nbl,2) - Qp(i,max(j-1,1),k,nbl,2)) * jyl * 0.5d0
               u_z = (Qp(i,j,min(k+1,NK(nbl)),nbl,2) - Qp(i,j,max(k-1,1),nbl,2)) * kzl * 0.5d0

               v_x = (Qp(min(i+1,NI(nbl)),j,k,nbl,3) - Qp(max(i-1,1),j,k,nbl,3)) * ixl * 0.5d0
               v_y = (Qp(i,min(j+1,NJ(nbl)),k,nbl,3) - Qp(i,max(j-1,1),k,nbl,3)) * jyl * 0.5d0
               v_z = (Qp(i,j,min(k+1,NK(nbl)),nbl,3) - Qp(i,j,max(k-1,1),nbl,3)) * kzl * 0.5d0

               w_x = (Qp(min(i+1,NI(nbl)),j,k,nbl,4) - Qp(max(i-1,1),j,k,nbl,4)) * ixl * 0.5d0
               w_y = (Qp(i,min(j+1,NJ(nbl)),k,nbl,4) - Qp(i,max(j-1,1),k,nbl,4)) * jyl * 0.5d0
               w_z = (Qp(i,j,min(k+1,NK(nbl)),nbl,4) - Qp(i,j,max(k-1,1),nbl,4)) * kzl * 0.5d0

               omx = w_y - v_z
               omy = u_z - w_x
               omz = v_x - u_y

               tke_sum  = tke_sum  + 0.5d0*rhl*(ul**2 + vl**2 + wl**2) * dcell
               enst_sum = enst_sum + 0.5d0*rhl*(omx**2 + omy**2 + omz**2) * dcell

            enddo
         enddo
      enddo
   enddo

endif

tke    = tke_sum  / Vdomain
enstpt = enst_sum / Vdomain

END SUBROUTINE






SUBROUTINE covo_error()
use declare_variables
implicit none

real :: xl, yl, r2, v_exact, v_num, error_max
real :: xc, yc, C_vor, R_vor
real :: dist, mindist
real :: dudx, dudy, dvdx, dvdy, omega_z, vort_mag
real :: dx_local, dy_local
integer :: jmid, jj, ip1, im1, jp1, jm1
character(60) :: cfile, vfile

! --- Vortex parameters ---
C_vor = 0.02d0
R_vor = 1.d0

! Vortex convects at u_inf = 1 in x-direction
! Wrap back into periodic domain [-8, 8]
xc = mod(1.d0 * time + 8.d0, 16.d0) - 8.d0
yc = 0.d0

! Find j-index closest to y = yc (general, works for any grid)
mindist = 1.d10
jmid = 1
do jj = 1, NJ(1)
   dist = abs(ygrid(1, jj, 1, 1) - yc)
   if (dist .lt. mindist) then
      mindist = dist
      jmid = jj
   endif
enddo


! compute and print linf error in v along centreline

error_max = 0.d0

do i = 1, NI(1)
   xl = xgrid(i, jmid, 1, 1)
   yl = ygrid(i, jmid, 1, 1)
   r2 = ((xl - xc)**2 + (yl - yc)**2) / R_vor**2
   v_exact = C_vor * (xl - xc) / R_vor**2 * exp(-r2 / 2.d0)
   v_num = Qp(i, jmid, 1, 1, 3)
   error_max = max(error_max, abs(v_num - v_exact))
enddo

print*, 'COVO L-inf error in v:', error_max
write(fresidual,*) 'COVO_ERROR:', error_max


! dump centreline v velocity profile 
!    x   v_numerical   v_exact

write(cfile,'(a,i3.3,a)') 'covo_centreline_N', NImax, '.dat'
open(unit=77, file=cfile, form='formatted')
write(77,*) '# x   v_numerical   v_exact'

do i = 1, NI(1)
   xl = xgrid(i, jmid, 1, 1)
   yl = ygrid(i, jmid, 1, 1)
   r2 = ((xl - xc)**2 + (yl - yc)**2) / R_vor**2
   v_exact = C_vor * (xl - xc) / R_vor**2 * exp(-r2 / 2.d0)
   v_num = Qp(i, jmid, 1, 1, 3)
   write(77,*) xl, v_num, v_exact
enddo

close(77)
print*, 'Wrote centreline data to ', trim(cfile)



! dump 2d vorticity magnitude field for contour plots
! uses 2nd order central diffs with periodic wrapping
! blank lines between j-rows so numpy doesnt have a seizure


write(vfile,'(a,i3.3,a)') 'covo_vorticity_N', NImax, '.dat'
open(unit=78, file=vfile, form='formatted')
write(78,*) '# x   y   vorticity_magnitude'

do j = 1, NJ(1)
   do i = 1, NI(1)
      ! Periodic index wrapping
      ip1 = i + 1; if (ip1 .gt. NI(1)) ip1 = 2
      im1 = i - 1; if (im1 .lt. 1)     im1 = NI(1) - 1
      jp1 = j + 1; if (jp1 .gt. NJ(1)) jp1 = 2
      jm1 = j - 1; if (jm1 .lt. 1)     jm1 = NJ(1) - 1

      dx_local = xgrid(ip1,j,1,1) - xgrid(im1,j,1,1)
      dy_local = ygrid(i,jp1,1,1) - ygrid(i,jm1,1,1)

      dvdx = (Qp(ip1,j,1,1,3) - Qp(im1,j,1,1,3)) / dx_local
      dudy = (Qp(i,jp1,1,1,2) - Qp(i,jm1,1,1,2)) / dy_local

      omega_z = dvdx - dudy
      vort_mag = abs(omega_z)

      write(78,*) xgrid(i,j,1,1), ygrid(i,j,1,1), vort_mag
   enddo
   write(78,*)  ! blank line between j-rows
enddo

close(78)
print*, 'Wrote vorticity field to ', trim(vfile)

END SUBROUTINE



	  
	  !  SUBROUTINE RESTART_call	
	  !  use declare_variables
	  !  implicit none
	  !  integer nvars,n
	  !  open(fflow, file= 'flow.xyz', form = 'unformatted')
	   
	  !  read(fflow) nblocks
	  !  read(fflow) (NI(nbl),NJ(nbl),NK(nbl),nvars, nbl= 1, nblocks)
	  !  do nbl = 1, nblocks
		! 	! write(fflow)((((Qp(i,j,k,nbl,n), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl)),n =1, nprims) & 
		! 	! ,((((Qc(i,j,k,nbl,n), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl)),n =1, nconserv)  
    !   read(fflow)((((Qp(i,j,k,nbl,n), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl)),n =1, nprims) & 
		! 	,((((Qc(i,j,k,nbl,n), i= 1, NI(NBL)),j = 1, NJ(nbl)),k = 1,NK(nbl)),n =1, nconserv)        
	  !  enddo
	  !  close(fflow)
	  ! END



SUBROUTINE RESTART_call	
  use declare_variables
  implicit none
  integer nvars,n

  ! Temporary arrays to absorb metric data we don't need
  ! (metrics get recomputed in METRICS() anyway)
  real,allocatable :: dummy3d(:,:,:)

  open(fflow, file= 'flow.xyz', form = 'unformatted')
  
  read(fflow) nblocks
  read(fflow) (NI(nbl),NJ(nbl),NK(nbl),nvars, nbl= 1, nblocks)

  allocate(dummy3d(NImax,NJmax,NKmax))

  do nbl = 1, nblocks
  read(fflow)((((Qp(i,j,k,nbl,n), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl)),n=1,nprims) & 
  ,((((Qc(i,j,k,nbl,n), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl)),n=1,nconserv)&
  , (((dummy3d(i,j,k), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))&
  , (((dummy3d(i,j,k), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))&
  , (((dummy3d(i,j,k), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))&
  , (((dummy3d(i,j,k), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))&
  , (((dummy3d(i,j,k), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))&
  , (((dummy3d(i,j,k), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))&
  , (((dummy3d(i,j,k), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))&
  , (((dummy3d(i,j,k), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))&
  , (((dummy3d(i,j,k), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))
  enddo

  deallocate(dummy3d)
  close(fflow)
END



    ! old code to check for errors


! !********************* OUTPUT GRID & FLOW *****************************************************
!       SUBROUTINE OUTPUT(flag)	  
!       use declare_variables
!       implicit none	  
	  
!       integer n,nvars,flag
	  
	  


!     ! writing grid files------
	  
! 	  open(fgrid,file='grid.xyz',form='unformatted')
	  
! 	  write(fgrid) nblocks
! 	  write(fgrid) (NI(nbl),NJ(nbl),NK(nbl), nbl=1,nblocks)
	  
! 	  do nbl=1,nblocks
! 		write(fgrid) (((xgrid(i,j,k,nbl),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))&
! 		,			   (((ygrid(i,j,k,nbl),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))&
! 		,			   (((zgrid(i,j,k,nbl),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))
! 	  enddo
	  
! 	  close(fgrid)


!     ! writing flow files------


!     nvars = nprims + nconserv + 6 !for xi, yi, zi, xj, yj, zj

!     open(fflow,file='flow.xyz',form='unformatted')
	  
! 	  write(fflow) nblocks
! 	  write(fflow) (NI(nbl),NJ(nbl),NK(nbl),nvars,nbl=1,nblocks)


!     do nbl=1,nblocks
! 		write(fflow) ((((Qp(i,j,k,nbl,n),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl)),n=1,nprims) &
!     ,			   ((((Qc(i,j,k,nbl,n),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl)),n=1,nconserv) &
!     ,			   (((xi(i,j,k,nbl),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl)) &
!     ,			   (((yi(i,j,k,nbl),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl)) &
!     ,			   (((zi(i,j,k,nbl),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl)) &
!     ,			   (((xj(i,j,k,nbl),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl)) &
!     ,			   (((yj(i,j,k,nbl),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl)) &
!     ,			   (((zj(i,j,k,nbl),i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl)) 

!     ! written xj, yj, zj here too.

!     ! can check if these work, and then replace them with the metric terms (ix, iy, iz etc)

    



!     ! ------------write a restart file-----------------

!     ! add shit in



! 	  enddo


	  
!       END
	  
!       SUBROUTINE DEALLOCATE_ROUTINE()
!       use declare_variables
!       implicit none	  
	  
	  
!       END	

!       SUBROUTINE volume_integral()
!       use declare_variables
!       implicit none	  
!       real,dimension(NJmax,NKmax,nblocks,2) :: plane
!       real,dimension(NKmax,nblocks,2) :: line	  
!       real muavg, dcell

	  
!       END	  