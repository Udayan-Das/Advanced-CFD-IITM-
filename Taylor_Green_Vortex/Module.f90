      MODULE declare_variables
      implicit none
	  
      integer, parameter :: fresidual = 1, finput = 2, fgrid = 3, fflow = 4
      integer,allocatable :: NI(:), NJ(:), NK(:)
      integer i,j,k,nbl, nblocks
      integer NImax,NJmax,NKmax,Ptsmax	  
      integer restart,iter, nsteps, nprims, nconserv, rk_steps, grid2d
      integer perI, perJ, perK	
      
      integer viscous
	  
	  !Discretization scheme coefficients
      integer dscheme, dschemek, bscheme1, bscheme2   !Discretization scheme in (i,j) and k, Boundary schemes
      integer fscheme             !Filtering scheme in (i,j) and k	  
      real adisc, bdisc, alpha, adisck, bdisck, alphak, alpha_f
      real alpha1, alpha2	  
	  
      real time, time_step
      real Lx, Ly, Lz
      real Re, Mach, gamma, prandtl, T_ref
	  
      real,allocatable :: xgrid(:,:,:,:), ygrid(:,:,:,:), zgrid(:,:,:,:)
      real,allocatable :: Qp(:,:,:,:,:) 		!Primitive variables: Rh,u,v,w,P,T
      real,allocatable :: Qc(:,:,:,:,:), Qcnew(:,:,:,:,:), Qcini(:,:,:,:,:) 		!Conservative variables in current and next timestep : RhU,RhV,RhW,RhE
      real,allocatable :: mu(:,:,:,:) 		    !Non dimensional Viscosity set based on Sutherland's law
      real,allocatable :: Qpi(:,:,:,:,:), Qpj(:,:,:,:,:),Qpk(:,:,:,:,:) 		!Derivatives of Primitive variables: Rh,u,v,w,P,T
      real,allocatable :: Fflux(:,:,:,:,:), Gflux(:,:,:,:,:),Hflux(:,:,:,:,:),net_flux(:,:,:,:,:),fluxD(:,:,:,:,:) 	!Flux in I,J,K directions, netflux & Derivative of flux in corresponding directions  
      real,allocatable :: xi(:,:,:,:), yi(:,:,:,:),zi(:,:,:,:) 		!Derivatives of grid
      real,allocatable :: xj(:,:,:,:), yj(:,:,:,:),zj(:,:,:,:) 		!Derivatives of grid
      real,allocatable :: xk(:,:,:,:), yk(:,:,:,:),zk(:,:,:,:) 		!Derivatives of grid

      real,allocatable :: ix(:,:,:,:), iy(:,:,:,:), iz(:,:,:,:) 		!Metrics of grid
      real,allocatable :: jx(:,:,:,:), jy(:,:,:,:), jz(:,:,:,:) 		!Metrics of grid
      real,allocatable :: kx(:,:,:,:), ky(:,:,:,:), kz(:,:,:,:) 		!Metrics of grid	  
      real,allocatable :: Jac(:,:,:,:), res(:)                        	!Jacobian, residual

      real,allocatable :: fac_qini(:),  fac_RK(:)                                                 !Coefficients for RK Time integration
      real,allocatable :: fcoeff(:), AMF(:), ACF(:), APF(:), AMD(:), ACD(:), APD(:)               !Coefficients for filtering, Discretization	  
      real,allocatable :: fdisc1(:), fdisc2(:)                                                    !Discretization coefficients pt1, pt2
	  
	  !******* For specific test cases*******************************************************
      integer testcase, taylor, covo, animfreq
      real,allocatable :: tked(:,:,:,:), enst(:,:,:,:)   		    !Turbulent Kinetic energy, Enstrophy
      real tke, enstpt


      ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

      ! VEEEERY TEMPRARY FIX RIGHT HERE. WOULD NEED TO ADDRESS SOME SHIT LATER
      

      ! real,allocatable :: PHIi(:,:,:,:,:), PHIj(:,:,:,:,:), PHIk(:,:,:,:,:)
      real :: dcell, Vdomain
      

      ! aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
	  
      END MODULE