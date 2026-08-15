!********************* SOLVER ROUTINES *****************************************************
      SUBROUTINE UNSTEADY(stepl)
      use declare_variables
      implicit none	  

      ! see slide 27 for the math
	  
      real Ucont,Vcont,Wcont
      real rhl, ul, vl, wl, pl, Tl, El, mul, vol, Elocal
      real ixl, jxl, kxl, iyl, jyl, kyl, izl, jzl, kzl	
      real uil, ujl, ukl, vil, vjl, vkl, wil, wjl, wkl, Til, Tjl, Tkl
      real bx, by, bz, u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, T_x, T_y, T_z
      real Txx, Txy, Txz, Tyy, Tyz, Tzz
      real div2b3, facprM, Sterm
      integer stepl, var	 	  
      real Pr
	  
	  !*********** Inviscid flux estimation ***********************


      Do nbl = 1,nblocks
         Do k = 1,NK(nbl)
            Do j = 1,NJ(nbl)
               Do i = 1,NI(nbl)


                  vol = 1.d0/Jac(i,j,k,nbl)

                  rhl = Qp(i,j,k,nbl,1)
                  ul = Qp(i,j,k,nbl,2)
                  vl = Qp(i,j,k,nbl,3)
                  wl = Qp(i,j,k,nbl,4)
                  pl = Qp(i,j,k,nbl,5)
                  Tl = Qp(i,j,k,nbl,6)
                  El = Tl/(gamma*(gamma-1.d0)*Mach**2) + 0.5d0*(ul**2 + vl**2 +wl**2)
                  Elocal = El

                  ! do not mind the two variabless for elocal. im too lazy to make a unif code
                  
                  

                  ixl = ix(i,j,k,nbl)
                  iyl = iy(i,j,k,nbl)
                  izl = iz(i,j,k,nbl)
                  jxl = jx(i,j,k,nbl)
                  jyl = jy(i,j,k,nbl)
                  jzl = jz(i,j,k,nbl)
                  kxl = kx(i,j,k,nbl)
                  kyl = ky(i,j,k,nbl)
                  kzl = kz(i,j,k,nbl)

                  ! the contravariant velocities

                  Ucont = ul*ixl + vl*iyl + wl*izl
                  Vcont = ul*jxl + vl*jyl + wl*jzl
                  Wcont = ul*kxl + vl*kyl + wl*kzl


                  ! the inviscid fluxes F G H

                  ! we do negative flux to send the vortex to the rhs. note this as we are doing the net flux update in the 
                  ! weno subroutine itself, so the updated bits would not need the minus sign!

                  ! if WENO == YES:

                  if (dscheme.le.4) then

                     Fflux(i,j,k,nbl,1) = -rhl*Ucont * vol
                     Fflux(i,j,k,nbl,2) = -((rhl*ul*Ucont + pl*ixl) * vol)
                     Fflux(i,j,k,nbl,3) = -((rhl*vl*Ucont + pl*iyl) * vol)
                     Fflux(i,j,k,nbl,4) = -((rhl*wl*Ucont + pl*izl) * vol)
                     Fflux(i,j,k,nbl,5) = -((rhl*Elocal*Ucont + pl*Ucont) * vol)


                     Gflux(i,j,k,nbl,1) = -rhl*Vcont * vol
                     Gflux(i,j,k,nbl,2) = -((rhl*ul*Vcont + pl*jxl) * vol)
                     Gflux(i,j,k,nbl,3) = -((rhl*vl*Vcont + pl*jyl) * vol)
                     Gflux(i,j,k,nbl,4) = -((rhl*wl*Vcont + pl*jzl) * vol)
                     Gflux(i,j,k,nbl,5) = -((rhl*Elocal*Vcont + pl*Vcont) * vol)


                     Hflux(i,j,k,nbl,1) = -rhl*Wcont * vol
                     Hflux(i,j,k,nbl,2) = -((rhl*ul*Wcont + pl*kxl) * vol)
                     Hflux(i,j,k,nbl,3) = -((rhl*vl*Wcont + pl*kyl) * vol)
                     Hflux(i,j,k,nbl,4) = -((rhl*wl*Wcont + pl*kzl) * vol)
                     Hflux(i,j,k,nbl,5) = -((rhl*Elocal*Wcont + pl*Wcont) * vol)

                  else

                     ! so flux splitting dies if we try to do zis.
                     ! does zis mean we are disregarding aristotle?
                     ! nein.

                     ! one must put le minus sign back in le weno disc subroutine end
                     ! ne pas le faire ferait de cette personne un idiot d'une ampleur inimaginable.

                     Fflux(i,j,k,nbl,1) = rhl*Ucont * vol
                     Fflux(i,j,k,nbl,2) = ((rhl*ul*Ucont + pl*ixl) * vol)
                     Fflux(i,j,k,nbl,3) = ((rhl*vl*Ucont + pl*iyl) * vol)
                     Fflux(i,j,k,nbl,4) = ((rhl*wl*Ucont + pl*izl) * vol)
                     Fflux(i,j,k,nbl,5) = ((rhl*Elocal*Ucont + pl*Ucont) * vol)


                     Gflux(i,j,k,nbl,1) = rhl*Vcont * vol
                     Gflux(i,j,k,nbl,2) = ((rhl*ul*Vcont + pl*jxl) * vol)
                     Gflux(i,j,k,nbl,3) = ((rhl*vl*Vcont + pl*jyl) * vol)
                     Gflux(i,j,k,nbl,4) = ((rhl*wl*Vcont + pl*jzl) * vol)
                     Gflux(i,j,k,nbl,5) = ((rhl*Elocal*Vcont + pl*Vcont) * vol)


                     Hflux(i,j,k,nbl,1) = rhl*Wcont * vol
                     Hflux(i,j,k,nbl,2) = ((rhl*ul*Wcont + pl*kxl) * vol)
                     Hflux(i,j,k,nbl,3) = ((rhl*vl*Wcont + pl*kyl) * vol)
                     Hflux(i,j,k,nbl,4) = ((rhl*wl*Wcont + pl*kzl) * vol)
                     Hflux(i,j,k,nbl,5) = ((rhl*Elocal*Wcont + pl*Wcont) * vol)

                  endif



                  ! else if WENO == NO:

                  ! the same fluxes without the ninus sign.
                  ! realisation: that aint gonna work chico.as soon as you turn on vicscous its gonna explode and crash and burn

                  ! endif



                  

               enddo
            enddo
         enddo
      enddo



      ! horrible mathematics here:
      ! can we even do this here? wakaranai. life is wakaranai.


      ! what is looove?
      ! baby dont hurt me
      ! dont hurt me no more

      net_flux = 0.d0	! initialize net flux to zero. we will add the inviscid and viscous fluxes to it later.


      if (dscheme.ge.5) then

         call WENO_I(Fflux, fluxD, nconserv)	
         net_flux = net_flux + fluxD
         call WENO_J(Gflux, fluxD, nconserv)
         net_flux = net_flux + fluxD
         if (grid2d.ne.1) then
            call WENO_K(Hflux, fluxD, nconserv)
            net_flux = net_flux + fluxD
         endif
      endif
	  



	  !*********** Viscous flux estimation ************************

      if (viscous.eq.1) then



      ! do we leave these as they are for WENO? so many questions
         

  

      if (dscheme.eq.1 .or. dscheme.eq.2) then
      
         call DISCRETIZATION_I_EXP(Qp, Qpi, nprims)	
         call DISCRETIZATION_J_EXP(Qp, Qpj, nprims)
         if (grid2d.ne.1) then
            call DISCRETIZATION_K_EXP(Qp, Qpk, nprims)
         else
            Qpk = 0.d0
         endif

      else if (dscheme.eq.3 .or. dscheme.eq.4) then

         call DISCRETIZATION_I_COMP(Qp,Qpi,nprims)
         call DISCRETIZATION_J_COMP(Qp,Qpj,nprims)
         if (grid2d.ne.1) then
            call DISCRETIZATION_K_COMP(Qp,Qpk,nprims)
         else
            Qpk = 0.d0
         endif


      ! should we do this?
         ! why not?
         ! yes well it would take 50000 years per timestep to compute
         ! but this is what aristotle intended.
         ! so we shall do this.


         ! we can go explicit schemes if we are a normie but we aint a normie mwaaaahahahahahaha



         

      ! else
      !    call DISCRETIZATION_I_COMP(Qp,Qpi,nprims)
      !    call DISCRETIZATION_J_COMP(Qp,Qpj,nprims)
      !    if (grid2d.ne.1) then
      !       call DISCRETIZATION_K_COMP(Qp,Qpk,nprims)
      !    else
      !       Qpk = 0.d0
      !    endif


         ! apparently one should do explicit here?
         ! apparently it doesnt blow up for non periodic cases?
         ! apparently it is recommended by professor claude?

         ! nani sore?

      else
         call DISCRETIZATION_I_EXP(Qp,Qpi,nprims)
         call DISCRETIZATION_J_EXP(Qp,Qpj,nprims)
         if (grid2d.ne.1) then
            call DISCRETIZATION_K_EXP(Qp,Qpk,nprims)
         else
            Qpk = 0.d0
         endif

      endif


      ! --------------------------------------------------------

      do nbl = 1,nblocks
         do k = 1,NK(nbl)
            do j = 1,NJ(nbl)
               do i = 1,NI(nbl)

                  rhl = Qp(i,j,k,nbl,1)
                  ul = Qp(i,j,k,nbl,2)
                  vl = Qp(i,j,k,nbl,3)
                  wl = Qp(i,j,k,nbl,4)
                  pl = Qp(i,j,k,nbl,5)
                  Tl = Qp(i,j,k,nbl,6)
                  El = Qp(i,j,k,nbl,6)/(gamma*(gamma-1.d0)*Mach**2) + 0.5d0*(ul**2+vl**2+wl**2)
                  


                  ixl = ix(i,j,k,nbl)
                  iyl = iy(i,j,k,nbl)	
                  izl = iz(i,j,k,nbl)

                  jxl = jx(i,j,k,nbl)
                  jyl = jy(i,j,k,nbl)
                  jzl = jz(i,j,k,nbl)

                  kxl = kx(i,j,k,nbl)
                  kyl = ky(i,j,k,nbl)	
                  kzl = kz(i,j,k,nbl)

                  ! populate the derivatives of primitive variables local:

                  uil = Qpi(i,j,k,nbl,2)
                  ujl = Qpj(i,j,k,nbl,2)
                  ukl = Qpk(i,j,k,nbl,2)

                  vil = Qpi(i,j,k,nbl,3)
                  vjl = Qpj(i,j,k,nbl,3)
                  vkl = Qpk(i,j,k,nbl,3)

                  wil = Qpi(i,j,k,nbl,4)
                  wjl = Qpj(i,j,k,nbl,4)
                  wkl = Qpk(i,j,k,nbl,4)

                  Til = Qpi(i,j,k,nbl,6)
                  Tjl = Qpj(i,j,k,nbl,6)
                  Tkl = Qpk(i,j,k,nbl,6)

                  ! cartesian derivatives of the primitive variables (local ofc):

                  ! god these are so wrong. i shall preserve them in comments as a monument to my sins

                  ! u_x = uil*ixl + ujl*iyl + ukl*izl
                  ! u_y = uil*jxl + ujl*jyl + ukl*jzl
                  ! u_z = uil*kxl + ujl*kyl + ukl*kzl

                  ! v_x = vil*ixl + vjl*iyl + vkl*izl
                  ! v_y = vil*jxl + vjl*jyl + vkl*jzl
                  ! v_z = vil*kxl + vjl*kyl + vkl*kzl

                  ! w_x = wil*ixl + wjl*iyl + wkl*izl
                  ! w_y = wil*jxl + wjl*jyl + wkl*jzl
                  ! w_z = wil*kxl + wjl*kyl + wkl*kzl

                  ! T_x = Til*ixl + Tjl*iyl + Tkl*izl
                  ! T_y = Til*jxl + Tjl*jyl + Tkl*jzl
                  ! T_z = Til*kxl + Tjl*kyl + Tkl*kzl

                  ! corrected cartesian derivatives of the primitive variables (local ofc):

                  u_x = uil*ixl + ujl*jxl + ukl*kxl
                  u_y = uil*iyl + ujl*jyl + ukl*kyl
                  u_z = uil*izl + ujl*jzl + ukl*kzl

                  v_x = vil*ixl + vjl*jxl + vkl*kxl
                  v_y = vil*iyl + vjl*jyl + vkl*kyl
                  v_z = vil*izl + vjl*jzl + vkl*kzl

                  w_x = wil*ixl + wjl*jxl + wkl*kxl
                  w_y = wil*iyl + wjl*jyl + wkl*kyl
                  w_z = wil*izl + wjl*jzl + wkl*kzl

                  T_x = Til*ixl + Tjl*jxl + Tkl*kxl
                  T_y = Til*iyl + Tjl*jyl + Tkl*kyl
                  T_z = Til*izl + Tjl*jzl + Tkl*kzl







                  ! now we can go ahead and calculate the stresses 

                  ! ---------------COMPUTE VISCOUS STRESSES-----------------

                  div2b3 = 2.d0/3.d0 * (u_x + v_y + w_z)

                  ! using sutheerland's law to estimate the non dim viscosity? it is non dim tho

                  mul = (Tl)**(3.d0/2.d0) * (1.d0 + 110.4d0/T_ref) / (Tl + 110.4d0/T_ref)	! non dimensional viscosity. check if this is correct. look at the ipad. also check if we need to add a reference viscosity here or not. i think we dont need to since we are non dimensionalizing with the reference viscosity. but check again.

                  ! lambda = - (2/3)*mu

                  Txx = (2.d0 * u_x - div2b3)*mul/Re
                  Tyy = (2.d0 * v_y - div2b3)*mul/Re
                  Tzz = (2.d0 * w_z - div2b3)*mul/Re

                  ! invoking symmetry straight off in the stress tensor. Txy = Tyx, Txz = Tzx, Tyz = Tzy

                  Txy = (u_y + v_x)*mul/Re
                  Txz = (u_z + w_x)*mul/Re
                  Tyz = (v_z + w_y)*mul/Re

                  Pr = prandtl !i think this might be needed

                  facprM = 1.d0 / ((gamma-1.d0)*Mach**2*Pr)


                  ! -------------------------------------------------

                  ! where did these come from?

                  ! probably incorrect with an extra facprm. commenting for now

                  ! bx =  facprM * (ul*Txx + vl*Txy + wl*Txz + mul/Re*facprM*T_x)
                  ! by =  facprM * (ul*Txy + vl*Tyy + wl*Tyz + mul/Re*facprM*T_y)
                  ! bz =  facprM * (ul*Txz + vl*Tyz + wl*Tzz + mul/Re*facprM*T_z)


                  ! probably correct ones:

                  ! Remove the leading facprM * from bx, by, and bz
                  bx = ul*Txx + vl*Txy + wl*Txz + (mul/Re)*facprM*T_x
                  by = ul*Txy + vl*Tyy + wl*Tyz + (mul/Re)*facprM*T_y
                  bz = ul*Txz + vl*Tyz + wl*Tzz + (mul/Re)*facprM*T_z





                  ! ---------------the rest of the fluxes-----------------

                  ! Fflux(i,j,k,nbl,1) = -rhl*Ucont * vol
                  ! Fflux(i,j,k,nbl,2) = -((rhl*ul*Ucont + pl*ixl) * vol)
                  ! Fflux(i,j,k,nbl,3) = -((rhl*vl*Ucont + pl*iyl) * vol)
                  ! Fflux(i,j,k,nbl,4) = -((rhl*wl*Ucont + pl*izl) * vol)
                  ! Fflux(i,j,k,nbl,5) = -((rhl*Elocal*Ucont + pl*Ucont) * vol)


                  ! Gflux(i,j,k,nbl,1) = -rhl*Vcont * vol
                  ! Gflux(i,j,k,nbl,2) = -((rhl*ul*Vcont + pl*jxl) * vol)
                  ! Gflux(i,j,k,nbl,3) = -((rhl*vl*Vcont + pl*jyl) * vol)
                  ! Gflux(i,j,k,nbl,4) = -((rhl*wl*Vcont + pl*jzl) * vol)
                  ! Gflux(i,j,k,nbl,5) = -((rhl*Elocal*Vcont + pl*Vcont) * vol)


                  ! Hflux(i,j,k,nbl,1) = -rhl*Wcont * vol
                  ! Hflux(i,j,k,nbl,2) = -((rhl*ul*Wcont + pl*kxl) * vol)
                  ! Hflux(i,j,k,nbl,3) = -((rhl*vl*Wcont + pl*kyl) * vol)
                  ! Hflux(i,j,k,nbl,4) = -((rhl*wl*Wcont + pl*kzl) * vol)
                  ! Hflux(i,j,k,nbl,5) = -((rhl*Elocal*Wcont + pl*Wcont) * vol)


                  ! structure here: 
                  ! total flux = -flux(invisc) + flux(viscous) = essentially the RHS of the eqn.


                  if (dscheme.le.4) then

                     ! we add the viscous fluxes to the inviscid fluxes here. we can also do this after discretization but it might be more expensive to do it there since we will have to do the disc of the viscous fluxes as well. we will see which one is better later.


                     Fflux(i,j,k,nbl,1) = Fflux(i,j,k,nbl,1) 
                     Fflux(i,j,k,nbl,2) = Fflux(i,j,k,nbl,2) + (Txx*ixl + Txy*iyl + Txz*izl)*vol
                     Fflux(i,j,k,nbl,3) = Fflux(i,j,k,nbl,3) + (Txy*ixl + Tyy*iyl + Tyz*izl)*vol
                     Fflux(i,j,k,nbl,4) = Fflux(i,j,k,nbl,4) + (Txz*ixl + Tyz*iyl + Tzz*izl)*vol
                     Fflux(i,j,k,nbl,5) = Fflux(i,j,k,nbl,5) + (ixl*bx + iyl*by + izl*bz)*vol !are we sure its correct?

                     Gflux(i,j,k,nbl,1) = Gflux(i,j,k,nbl,1)
                     Gflux(i,j,k,nbl,2) = Gflux(i,j,k,nbl,2) + (Txx*jxl + Txy*jyl + Txz*jzl)*vol
                     Gflux(i,j,k,nbl,3) = Gflux(i,j,k,nbl,3) + (Txy*jxl + Tyy*jyl + Tyz*jzl)*vol
                     Gflux(i,j,k,nbl,4) = Gflux(i,j,k,nbl,4) + (Txz*jxl + Tyz*jyl + Tzz*jzl)*vol
                     Gflux(i,j,k,nbl,5) = Gflux(i,j,k,nbl,5) + (jxl*bx + jyl*by + jzl*bz)*vol

                     Hflux(i,j,k,nbl,1) = Hflux(i,j,k,nbl,1)
                     Hflux(i,j,k,nbl,2) = Hflux(i,j,k,nbl,2) + (Txx*kxl + Txy*kyl + Txz*kzl)*vol
                     Hflux(i,j,k,nbl,3) = Hflux(i,j,k,nbl,3) + (Txy*kxl + Tyy*kyl + Tyz*kzl)*vol
                     Hflux(i,j,k,nbl,4) = Hflux(i,j,k,nbl,4) + (Txz*kxl + Tyz*kyl + Tzz*kzl)*vol
                     Hflux(i,j,k,nbl,5) = Hflux(i,j,k,nbl,5) + (kxl*bx + kyl*by + kzl*bz)*vol

                  else

                     Fflux(i,j,k,nbl,1) = 0.d0 
                     Fflux(i,j,k,nbl,2) = (Txx*ixl + Txy*iyl + Txz*izl)*vol
                     Fflux(i,j,k,nbl,3) = (Txy*ixl + Tyy*iyl + Tyz*izl)*vol
                     Fflux(i,j,k,nbl,4) = (Txz*ixl + Tyz*iyl + Tzz*izl)*vol
                     Fflux(i,j,k,nbl,5) = (ixl*bx + iyl*by + izl*bz)*vol !are we sure its correct?

                     Gflux(i,j,k,nbl,1) = 0.d0
                     Gflux(i,j,k,nbl,2) = (Txx*jxl + Txy*jyl + Txz*jzl)*vol
                     Gflux(i,j,k,nbl,3) = (Txy*jxl + Tyy*jyl + Tyz*jzl)*vol
                     Gflux(i,j,k,nbl,4) = (Txz*jxl + Tyz*jyl + Tzz*jzl)*vol
                     Gflux(i,j,k,nbl,5) = (jxl*bx + jyl*by + jzl*bz)*vol

                     Hflux(i,j,k,nbl,1) = 0.d0
                     Hflux(i,j,k,nbl,2) = (Txx*kxl + Txy*kyl + Txz*kzl)*vol
                     Hflux(i,j,k,nbl,3) = (Txy*kxl + Tyy*kyl + Tyz*kzl)*vol
                     Hflux(i,j,k,nbl,4) = (Txz*kxl + Tyz*kyl + Tzz*kzl)*vol
                     Hflux(i,j,k,nbl,5) = (kxl*bx + kyl*by + kzl*bz)*vol

                  endif






                  
               enddo
            enddo
         enddo
      enddo

      endif

      ! -------------end of viscous flux estimation------------------


      ! net_flux = 0.d0	! initialize net flux to zero. we will add the inviscid and viscous fluxes to it later.

      ! usually net_flux = visous - invisc

      ! in the eqn keep d/dt(U/J) is the only thing in the LHS of the eqn. all of the other terms are in the RHS. thats the residue. and we use the prev timestep value there?
      ! either way this needs some reading.

      ! this is just one of the ways to do it. 

      ! i think the shit below is just inviscid? life is wakaranai

      ! netflux is the rhs. 

      ! if (dscheme.eq.1 .or. dscheme.eq.2) then

      !    call DISCRETIZATION_I_EXP(Fflux, fluxD, nconserv)	! gives us the derivatives of the fluxes in the i direction. we will do this for all the directions and then add them up to get the net flux. we can also add the viscous fluxes to the inviscid fluxes before doing this. we will see which one is better later.
      !    net_flux = net_flux + fluxD
      !    call DISCRETIZATION_J_EXP(Gflux, fluxD, nconserv)
      !    net_flux = net_flux + fluxD
      !    call DISCRETIZATION_K_EXP(Hflux, fluxD, nconserv)
      !    net_flux = net_flux + fluxD

      ! else if (dscheme.eq.3 .or. dscheme.eq.4) then

      !    call DISCRETIZATION_I_COMP(Fflux, fluxD, nconserv)	! gives us the derivatives of the fluxes in the i direction. we will do this for all the directions and then add them up to get the net flux. we can also add the viscous fluxes to the inviscid fluxes before doing this. we will see which one is better later.
      !    net_flux = net_flux + fluxD
      !    call DISCRETIZATION_J_COMP(Gflux, fluxD, nconserv)
      !    net_flux = net_flux + fluxD
      !    call DISCRETIZATION_K_COMP(Hflux, fluxD, nconserv)
      !    net_flux = net_flux + fluxD

      ! endif





      ! NON NON NON NON NON NON
      ! ce n'est pas bon, c'est horrible, horrible ! aaaaa !


      ! THIS IS THE SUM OF THE SHIT STILL FLUXING IS IT? IS IT NOT? AAAAAAAAAAAA

      if (dscheme.eq.1 .or. dscheme.eq.2) then

         call DISCRETIZATION_I_EXP(Fflux, fluxD, nconserv)	
         net_flux = net_flux + fluxD
         call DISCRETIZATION_J_EXP(Gflux, fluxD, nconserv)
         net_flux = net_flux + fluxD
         if (grid2d.ne.1) then
            call DISCRETIZATION_K_EXP(Hflux, fluxD, nconserv)
            net_flux = net_flux + fluxD
         endif

      else if (dscheme.eq.3 .or. dscheme.eq.4) then

         call DISCRETIZATION_I_COMP(Fflux, fluxD, nconserv)	
         net_flux = net_flux + fluxD
         call DISCRETIZATION_J_COMP(Gflux, fluxD, nconserv)
         net_flux = net_flux + fluxD
         if (grid2d.ne.1) then
            call DISCRETIZATION_K_COMP(Hflux, fluxD, nconserv)
            net_flux = net_flux + fluxD
         endif

      else if (dscheme.ge.5 .and. viscous.eq.1) then
         call DISCRETIZATION_I_EXP(Fflux, fluxD, nconserv)
         net_flux = net_flux + fluxD
         call DISCRETIZATION_J_EXP(Gflux, fluxD, nconserv)
         net_flux = net_flux + fluxD
         if (grid2d.ne.1) then
            call DISCRETIZATION_K_EXP(Hflux, fluxD, nconserv)
            net_flux = net_flux + fluxD
         endif

      endif



      ! print*, fluxD

      ! pause

      

      ! time stepping rk4


      do var = 1,nconserv
         do nbl = 1,nblocks
            do k = 1,NK(nbl)
               do j = 1,NJ(nbl)
                  do i = 1,NI(nbl)

                     vol = 1.d0/Jac(i,j,k,nbl)
      
                     !********* Post processing steps (if any) ************************************

                     !*************** ESTIMATING PART OF NEW CONSERVATIVE VARIABLES ***********************************	
                        
                     ! 

                     Qcnew(i,j,k,nbl,var) = Qcnew(i,j,k,nbl,var) + (time_step * net_flux(i,j,k,nbl,var)) * fac_RK(stepl) / vol  !why not just *jac

                     ! fac_rk in the next stepl goes to fac_rk(stepl+1) 

                     !*************** ESTIMATING PART OF SUB-RK_STAGE CONSERVATIVE VARIABLE to estimate new net flux ***********************************

                     if (stepl <= rk_steps-1) then
                        Qc(i,j,k,nbl,var) = Qcini(i,j,k,nbl,var) + (time_step * net_flux(i,j,k,nbl,var)) * fac_qini(stepl+1) / vol 
                        ! but if stepl goes >=4, it blows up.

                     
                     

                     !*************** ESTIMATING MAXIMUM RESIDUAL ***********************************	 
                     else if (stepl == rk_steps) then
                        res(var) = max(res(var), abs(Qcnew(i,j,k,nbl,var) - Qcini(i,j,k,nbl,var)))


                        

                        ! check this stuff once my pc be dying
                     
                        

                     !*************** UPDATING CONSERVATIVE VARIABLES ***********************************

                     Qc(i,j,k,nbl,var) = Qcnew(i,j,k,nbl,var)

                     ! do we not update the conserv var?

                     endif

                     ! does the endif come after qc update?

                  enddo
               enddo
            enddo
         enddo
      enddo


      END



	  
      SUBROUTINE SET_PRIMITIVES()
      use declare_variables	 
      implicit none	 

      real rhl, ul, vl, wl, pl, Tl, El	  

      do nbl = 1,nblocks
         do k = 1,NK(nbl)
            do j = 1,NJ(nbl)
               do i = 1,NI(nbl)

                  rhl = Qc(i,j,k,nbl,1)
                  ul = Qc(i,j,k,nbl,2)/rhl
                  vl = Qc(i,j,k,nbl,3)/rhl
                  wl = Qc(i,j,k,nbl,4)/rhl
                  El = Qc(i,j,k,nbl,5)/rhl
                  Tl = (El - 0.5d0*(ul**2+vl**2+wl**2))*(gamma*(gamma-1.d0)*Mach**2)
					pl = rhl*Tl/(gamma*Mach**2)

                  Qp(i,j,k,nbl,1) = rhl
                  Qp(i,j,k,nbl,2) = ul
                  Qp(i,j,k,nbl,3) = vl
                  Qp(i,j,k,nbl,4) = wl
                  Qp(i,j,k,nbl,5) = pl
                  Qp(i,j,k,nbl,6) = Tl

               enddo
            enddo
         enddo
      enddo

      END	 
	  
!******************************************************************************************************************	  
!*********************************** PERIODIC DISCRETIZATION ROUTINES FOR VARIABLES *******************************
!******************************************************************************************************************	  

SUBROUTINE DISCRETIZATION_I_EXP(PHI,PHID,nvars)    !literally gives us the derivatives of whatever we throw into it PHI, in the i direc. returns PHID
use declare_variables	 
implicit none	  

integer var,ip1,ip2,im1,im2,nvars
real bb4, ab2	  
real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID



bb4 = bdisc/4.d0
ab2 = adisc/2.d0



Do var = 1,nvars   !iterating the variables



   Do nbl = 1,nblocks


      Do k = 1,NK(nbl)


         Do j = 1,NJ(nbl)



            Do i = 3,NI(nbl)-2

               ip2 = i+2
               ip1 = i+1  
               im1 = i-1
               im2 = i-2
               ! PHID(i,j,k,nbl,var) = (bb4*PHI(ip2,j,k,nbl,var) + ab2*PHI(ip1,j,k,nbl,var) - ab2*PHI(im1,j,k,nbl,var) - bb4*PHI(im2,j,k,nbl,var))
               PHID(i,j,k,nbl,var) = &
               bb4*(PHI(ip2,j,k,nbl,var) - PHI(im2,j,k,nbl,var)) &
               + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var))

            enddo

            ! Boundary points defined separately so that it dont fuck up.

            i=1
            ip2 = i+2
            ip1 = i+1  
            im1 = NI(nbl)-1
            im2 = NI(nbl)-2

            ! LLx = PHI(NI(nbl),j,k,nbl,var) - PHI(1,j,k,nbl,var)

            ! flow variables are periodic. no need for LLx

            ! PHID(i,j,k,nbl,var) = (bb4*(PHI(ip2,j,k,nbl,var) - PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var)))

            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(ip2,j,k,nbl,var) - PHI(im2,j,k,nbl,var)) &
            + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var))

            i=2
            ip2 = i+2  
            ip1 = i+1 
            im1 = i-1
            im2 = NI(nbl)-1

            ! check if the llx implementation is correct here

            ! PHID(i,j,k,nbl,var) = (bb4*(PHI(ip2,j,k,nbl,var) - PHI(im2,j,k,nbl,var))+ ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var)))

            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(ip2,j,k,nbl,var) - PHI(im2,j,k,nbl,var)) &
            + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var))

            i=NI(nbl)-1
            ip1 = i+1
            ip2 = 2
            im1 = i-1
            im2 = i-2

            ! PHID(i,j,k,nbl,var) = (bb4*(PHI(ip2,j,k,nbl,var) - PHI(im2,j,k,nbl,var))+ ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var)))

            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(ip2,j,k,nbl,var) - PHI(im2,j,k,nbl,var)) &
            + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var))
            
            i=NI(nbl)
            ip2 = 3
            ip1 = 2
            im1 = i-1
            im2 = i-2

            ! PHID(i,j,k,nbl,var) = (bb4*(PHI(ip2,j,k,nbl,var) - PHI(im2,j,k,nbl,var))+ ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var)))

            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(ip2,j,k,nbl,var) - PHI(im2,j,k,nbl,var)) &
            + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var))

            ! COMPLETE THEM



         enddo
      enddo
   enddo
enddo



END	






SUBROUTINE DISCRETIZATION_J_EXP(PHI,PHID,nvars)
use declare_variables	 
implicit none	 

integer var,jp1,jp2,jm1,jm2,nvars,coeff,jl
real bb4, ab2	  
real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID


bb4 = bdisc/4.d0
ab2 = adisc/2.d0



! this specific loop is for the j direction only.

Do var = 1,nvars   !iterating the variables


   Do nbl = 1,nblocks


      Do k = 1,NK(nbl)


         Do i = 1,NI(nbl)



            Do j = 3,NJ(nbl)-2

               jp2 = j+2
               jp1 = j+1  
               jm1 = j-1
               jm2 = j-2
               ! PHID(i,j,k,nbl,var) = (bb4*PHI(i,jp2,k,nbl,var) + ab2*PHI(i,jp1,k,nbl,var) - ab2*PHI(i,jm1,k,nbl,var) - bb4*PHI(i,jm2,k,nbl,var))
               PHID(i,j,k,nbl,var) = &
                                    bb4*(PHI(i,jp2,k,nbl,var) - PHI(i,jm2,k,nbl,var)) &
                                    + ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var))

            enddo

            ! Boundary points defined separately so that it dont fuck up.

            j=1
            jp2 = j+2
            jp1 = j+1  
            jm1 = NJ(nbl)-1
            jm2 = NJ(nbl)-2

            ! LLy = PHI(i,NJ(nbl),k,nbl) - PHI(i,1,k,nbl)  !not needed

            ! PHID(i,j,k,nbl,var) = (bb4*(PHI(i,jp2,k,nbl,var) - PHI(i,jm2,k,nbl,var))+ ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var)))
            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(i,jp2,k,nbl,var) - PHI(i,jm2,k,nbl,var)) &
            + ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var))


            j=2
            jp2 = j+2  
            jp1 = j+1 
            jm1 = j-1
            jm2 = NJ(nbl)-1

            ! check if the llx implementation is correct here

            ! PHID(i,j,k,nbl,var) = (bb4*(PHI(i,jp2,k,nbl,var) - PHI(i,jm2,k,nbl,var))+ ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var)))

            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(i,jp2,k,nbl,var) - PHI(i,jm2,k,nbl,var)) + &
            ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var))

            j=NJ(nbl)-1
            jp1 = j+1
            jp2 = 2
            jm1 = j-1
            jm2 = j-2

            ! PHID(i,j,k,nbl,var) = (bb4*(PHI(i,jp2,k,nbl,var) - PHI(i,jm2,k,nbl,var))+ ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var)))

            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(i,jp2,k,nbl,var) - PHI(i,jm2,k,nbl,var)) + &
            ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var))

            j=NJ(nbl)
            jp2 = 3
            jp1 = 2
            jm1 = j-1
            jm2 = j-2

            ! PHID(i,j,k,nbl,var) = (bb4*(PHI(i,jp2,k,nbl,var) - PHI(i,jm2,k,nbl,var))+ ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var)))
            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(i,jp2,k,nbl,var) - PHI(i,jm2,k,nbl,var)) + &
            ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var))
            ! COMPLETE THEM



         enddo
      enddo
   enddo


   
enddo 





END	

SUBROUTINE DISCRETIZATION_K_EXP(PHI,PHID,nvars)
use declare_variables	 
implicit none	  

integer var,kp1,kp2,km1,km2,nvars
real bb4, ab2	  
real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID




bb4 = bdisc/4.d0
ab2 = adisc/2.d0



! this specific loop is for the k direction only.

Do var = 1,nvars   !iterating the variables

   Do nbl = 1,nblocks


      Do i = 1,NI(nbl)


         Do j = 1,NJ(nbl)



            Do k = 3,NK(nbl)-2

               kp2 = k+2
               kp1 = k+1  
               km1 = k-1
               km2 = k-2
               ! PHID(i,j,k,nbl,var) = (bb4*PHI(i,j,kp2,nbl,var) + ab2*PHI(i,j,kp1,nbl,var) - ab2*PHI(i,j,km1,nbl,var) - bb4*PHI(i,j,km2,nbl,var))

               PHID(i,j,k,nbl,var) = &
                                    bb4*(PHI(i,j,kp2,nbl,var) - PHI(i,j,km2,nbl,var)) &
                                    + ab2*(PHI(i,j,kp1,nbl,var) - PHI(i,j,km1,nbl,var))

            enddo

            ! Boundary points defined separately so that it dont fuck up.

            k=1
            kp2 = k+2
            kp1 = k+1  
            km1 = NK(nbl)-1
            km2 = NK(nbl)-2

            ! LLz = PHI(i,NK(nbl),k,nbl) - PHI(i,1,k,nbl)

            PHID(i,j,k,nbl,var) = &
                                 bb4*(PHI(i,j,kp2,nbl,var) - PHI(i,j,km2,nbl,var)) &
                                 + ab2*(PHI(i,j,kp1,nbl,var) - PHI(i,j,km1,nbl,var))

            k=2
            kp2 = k+2  
            kp1 = k+1 
            km1 = k-1
            km2 = NK(nbl)-1

            ! check if the llx implementation is correct here

            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(i,j,kp2,nbl,var) - PHI(i,j,km2,nbl,var)) &
            + ab2*(PHI(i,j,kp1,nbl,var) - PHI(i,j,km1,nbl,var))

            k=NK(nbl)-1
            kp1 = k+1
            kp2 = 2
            km1 = k-1
            km2 = k-2

            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(i,j,kp2,nbl,var) - PHI(i,j,km2,nbl,var)) &
            + ab2*(PHI(i,j,kp1,nbl,var) - PHI(i,j,km1,nbl,var))

            k=NK(nbl)
            kp2 = 3
            kp1 = 2
            km1 = k-1
            km2 = k-2

            PHID(i,j,k,nbl,var) = &
            bb4*(PHI(i,j,kp2,nbl,var) - PHI(i,j,km2,nbl,var)) &
            + ab2*(PHI(i,j,kp1,nbl,var) - PHI(i,j,km1,nbl,var))
            
            ! COMPLETE THEM



         enddo
      enddo
   enddo

enddo

END	

SUBROUTINE DISCRETIZATION_K2D_EXP(PHI,PHID,nvars)
use declare_variables	 
implicit none	  

integer var,kp1,km1,nvars
real bb4, ab2	  
real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID

END		  

!****************************** COMPACT SCHEMES *******************************************************************

! SUBROUTINE DISCRETIZATION_I_COMP(PHI,PHID,nvars)
! use declare_variables	 
! implicit none	  

! integer var,ip1,ip2,im1,im2,nvars
! real bb4, ab2	  
! real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
! real,dimension(NImax) :: RHS	  

! END	

! SUBROUTINE DISCRETIZATION_J_COMP(PHI,PHID,nvars)
! use declare_variables	 
! implicit none	 

! integer var,jp1,jp2,jm1,jm2,nvars
! real bb4, ab2	  
! real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
! real,dimension(NJmax) :: RHS	  

! END	

! SUBROUTINE DISCRETIZATION_K_COMP(PHI,PHID,nvars)
! use declare_variables	 
! implicit none	  

! integer var,kp1,kp2,km1,km2,nvars
! real bb4, ab2	  
! real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
! real,dimension(NKmax) :: RHS	  


! END		  

SUBROUTINE DISCRETIZATION_I_COMP(PHI,PHID,nvars)
use declare_variables	 
implicit none	  

integer var,ip1,ip2,im1,im2,nvars, nm1
real bb4, ab2	  
real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
real,dimension(NImax) :: RHS	 

bb4 = bdisc/4.d0
ab2 = adisc/2.d0
do var = 1,nvars
do nbl = 1, nblocks
   do k=1, NK(nbl)
      do j=1, NJ(nbl)
         do i=3, NI(nbl)-2
            ip2 = i+2
            ip1 = i+1
            im1 = i-1
            im2 = i-2
            RHS(i) = bb4*(PHI(ip2,j,k,nbl,var) -PHI(im2, j, k,nbl,var))  + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var))
         enddo

         
         i = 1
         ip2 = i+2
         ip1 = i+1
         im1 = NI(nbl)-1
         im2 = NI(nbl)-2
         RHS(i) = bb4*(PHI(ip2,j,k,nbl,var) -PHI(im2, j, k,nbl,var))  + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var))

         i = 2
         ip2 = i+2
         ip1 = i+1
         im1 = i-1
         im2 = NI(nbl)-1
         RHS(i) = bb4*(PHI(ip2,j,k,nbl,var) -PHI(im2, j, k,nbl,var))  + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var))

         i = NI(nbl)-1
         ip2 = 2
         ip1 = i+1
         im1 = i-1
         im2 = i-2
         RHS(i)= bb4*(PHI(ip2,j,k,nbl,var) -PHI(im2, j, k,nbl,var))  + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var))
         
         i = NI(nbl)
         ip2 = 3
         ip1 = 2
         im1 = i-1
         im2 = i-2
         RHS(i) = bb4*(PHI(ip2,j,k,nbl,var) -PHI(im2, j, k,nbl,var))  + ab2*(PHI(ip1,j,k,nbl,var) - PHI(im1,j,k,nbl,var))
         nm1 = NI(nbl)-1			
         
         call TDMAP(1,NI(nbl)-1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)
         PHID(1:nm1, j,k,nbl,var)= RHS(1:nm1)
         PHID(NI(nbl),j,k,nbl,var) = PHID(1,j,k,nbl,var)
      
      enddo
   enddo
enddo	
enddo 	  

END SUBROUTINE

SUBROUTINE DISCRETIZATION_J_COMP(PHI,PHID,nvars)
use declare_variables	 
implicit none	 

integer var,jp1,jp2,jm1,jm2,nvars,nm1
real bb4, ab2	  
real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
real,dimension(NJmax) :: RHS
bb4 = bdisc/4.d0
ab2 = adisc/2.d0
do var = 1,nvars
do nbl = 1, nblocks
do k=1, NK(nbl)
   do i=1, NI(nbl)
      do j=3, NJ(nbl)-2
         jp2 = j+2
         jp1 = j+1
         jm1 = j-1
         jm2 = j-2
         RHS(j) = bb4*(PHI(i,jp2,k,nbl,var) -PHI(i, jm2, k,nbl,var))  + ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var))
      enddo
   
      
      j = 1
      jp2 = j+2
      jp1 = j+1
      jm1 = NJ(nbl)-1
      jm2 = NJ(nbl)-2
   
      RHS(j) = bb4*(PHI(i,jp2,k,nbl,var) -PHI(i, jm2, k,nbl,var))  + ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var))
      
      j= 2
      jp2 = j+2
      jp1 = j+1
      jm1 = j-1
      jm2 = NJ(nbl)-1
      RHS(j) = bb4*(PHI(i,jp2,k,nbl,var) -PHI(i, jm2, k,nbl,var))  + ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var))
      
      j = NJ(nbl)-1
      jp2 = 2
      jp1 =j+1
      jm1 =j-1
      jm2 =j-2
      RHS(j) = bb4*(PHI(i,jp2,k,nbl,var) -PHI(i, jm2, k,nbl,var))  + ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var))
      
      j = NJ(nbl)
      jp2 = 3
      jp1 = 2
      jm1 = j-1
      jm2 = j-2
      RHS(j) = bb4*(PHI(i,jp2,k,nbl,var) -PHI(i, jm2, k,nbl,var))  + ab2*(PHI(i,jp1,k,nbl,var) - PHI(i,jm1,k,nbl,var))   	  
      nm1 = NJ(nbl)-1			
      
      call TDMAP(1,NJ(nbl)-1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)
      PHID(i,1:nm1,k,nbl,var)= RHS(1:nm1)
      PHID(i,NJ(nbl),k,nbl,var) = PHID(i,1,k,nbl,var)	       
      enddo
enddo
enddo
enddo



END SUBROUTINE	

SUBROUTINE DISCRETIZATION_K_COMP(PHI,PHID,nvars)
use declare_variables	 
implicit none	  

integer var,kp1,kp2,km1,km2,nvars,nm1
real bb4, ab2	  
real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
real,dimension(NKmax) :: RHS
bb4 = bdisc/4.d0
ab2 = adisc/2.d0
do var = 1,nvars
do nbl = 1, nblocks
do i=1, NI(nbl)
   do j=1, NJ(nbl)
      do k=3, NK(nbl)-2

         kp2 = k+2
         kp1 = k+1
         km1 = k-1
         km2 = k-2
         RHS(k) = bb4*(PHI(i,j,kp2,nbl,var) -PHI(i,j,km2,nbl,var))  + ab2*(PHI(i,j,kp1,nbl,var) - PHI(i,j,km1,nbl,var))
      enddo
         
      
      k = 1
      kp2 = k+2
      kp1 = k+1
      km1 = NK(nbl)-1
      km2 = NK(nbl)-2
      
      RHS(k) = bb4*(PHI(i,j,kp2,nbl,var) - PHI(i,j,km2,nbl,var))  + ab2*(PHI(i,j,kp1,nbl,var) - PHI(i,j,km1,nbl,var))
      
      k= 2
      kp2 = k+2
      kp1 = k+1
      km1 = k-1
      km2 = NK(nbl)-1
      RHS(k) = bb4*(PHI(i,j,kp2,nbl,var) - PHI(i,j,km2,nbl,var))  + ab2*(PHI(i,j,kp1,nbl,var) - PHI(i,j,km1,nbl,var))
      
      k = NK(nbl)-1
      kp2 = 2
      kp1 =k+1
      km1 =k-1
      km2 =k-2
      RHS(k) = bb4*(PHI(i,j,kp2,nbl,var) - PHI(i,j,km2,nbl,var))  + ab2*(PHI(i,j,kp1,nbl,var) - PHI(i,j,km1,nbl,var))
      
      k = NK(nbl)
      kp2 = 3
      kp1 = 2
      km1 = k-1
      km2 = k-2
      RHS(k) = bb4*(PHI(i,j,kp2,nbl,var) - PHI(i,j,km2,nbl,var))  + ab2*(PHI(i,j,kp1,nbl,var) - PHI(i,j,km1,nbl,var))
         
      nm1 = NK(nbl)-1			
      
         call TDMAP(1,NK(nbl)-1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)
         PHID(i,j,1:nm1,nbl,var) = RHS(1:nm1)
         PHID(i,j,NK(nbl),nbl,var) = PHID(i,j,1,nbl,var)
   enddo
enddo
enddo
enddo	  
 

END SUBROUTINE






! -----------WENO SUBROUTINES-------------------

!=========================================================================
! WENO_I: WENO flux derivative in the I (xi) direction.
!
! Drop-in replacement for DISCRETIZATION_I_EXP when dscheme >= 5.
!
!   PHI  = inviscid flux at nodes (NO negative sign!)
!   PHID = output derivative. Has the negative sign baked in so the
!          caller just does:  net_flux = net_flux + PHID
!
! Internally reads Qc, Qp, Jac, ix, iy, iz from declare_variables.
!
! Calling code in UNSTEADY (replaces the dscheme 1-4 block):
!
!   if (dscheme.ge.5) then
!     call WENO_I(Fflux_raw, fluxD, nconserv)
!     net_flux = net_flux + fluxD
!     call WENO_J(Gflux_raw, fluxD, nconserv)
!     net_flux = net_flux + fluxD
!     if (grid2d.ne.1) then
!       call WENO_K(Hflux_raw, fluxD, nconserv)
!       net_flux = net_flux + fluxD
!     endif
!   endif
!
! dscheme selects variant:  5=JS, 6=Z, 7=CU6, 8=CU6-M
!=========================================================================

SUBROUTINE WENO_I(PHI, PHID, nvars)
use declare_variables
implicit none

integer, intent(in) :: nvars
real, dimension(NImax,NJmax,NKmax,nblocks,nvars), intent(in)  :: PHI
real, dimension(NImax,NJmax,NKmax,nblocks,nvars), intent(out) :: PHID

! --- 1D working arrays along one i-line ---
! real :: Fp(NImax), Fm(NImax)
! real :: Fhat_total(NImax)

real(8) :: Fp(NImax), Fm(NImax), Fhat_total(NImax)

! --- Scalars ---
real(8) :: alpha_w                        ! global max wave speed (I-dir)
real(8) :: cl, Ucont_l, grad_xi_l, vol_l  ! for alpha computation
real(8) :: fhat_plus, fhat_minus           ! reconstructed interface fluxes
real(8) :: f0, f1, f2, f3                 ! candidate fluxes (left-biased)
real(8) :: g0, g1, g2, g3                 ! candidate fluxes (right-biased)
real(8) :: b0, b1, b2, b6                 ! smoothness indicators
real(8) :: tau5, tau6                      ! global smoothness refs
real(8) :: a0, a1, a2, a3, asum           ! unnormalised weights

! Index helpers
integer :: var, im2, im1, ic, ip1, ip2, ip3, im1f
integer :: NI_l

! Scheme parameters
real(8) :: eps_js, eps_z
real(8) :: d0, d1, d2, d3                 ! optimal weights
real(8) :: C_cu, q_cu                     ! CU6 parameters
logical :: is_cu6                         ! true for dscheme 7 or 8

! ================================================================
! Set scheme parameters
! ================================================================
eps_js = 1.d-6
eps_z  = 1.d-40

if (dscheme.eq.5 .or. dscheme.eq.6) then
  ! 5th-order schemes (3 stencils)
  d0 = 1.d0/10.d0
  d1 = 6.d0/10.d0
  d2 = 3.d0/10.d0
  is_cu6 = .false.
else
  ! 6th-order schemes (4 stencils)
  d0 = 1.d0/20.d0
  d1 = 9.d0/20.d0
  d2 = 9.d0/20.d0
  d3 = 1.d0/20.d0
  is_cu6 = .true.
  if (dscheme.eq.7) then
    C_cu = 20.d0;   q_cu = 1.d0
  else
    C_cu = 1000.d0;  q_cu = 4.d0
  endif
endif

! ================================================================
! STEP 1: Global max wave speed in I-direction
! ================================================================
! alpha_w = max over ALL grid points of { |U_cont| + c * |grad_xi| }
! No /J because metrics already include Jacobian.

alpha_w = 0.d0
do nbl = 1, nblocks
  do k = 1, NK(nbl)
    do j = 1, NJ(nbl)
      do i = 1, NI(nbl)
        cl = dsqrt(gamma * Qp(i,j,k,nbl,5) / Qp(i,j,k,nbl,1))
        Ucont_l = Qp(i,j,k,nbl,2)*ix(i,j,k,nbl) &
                + Qp(i,j,k,nbl,3)*iy(i,j,k,nbl) &
                + Qp(i,j,k,nbl,4)*iz(i,j,k,nbl)
        grad_xi_l = dsqrt(ix(i,j,k,nbl)**2 &
                        + iy(i,j,k,nbl)**2 &
                        + iz(i,j,k,nbl)**2)
        alpha_w = dmax1(alpha_w, dabs(Ucont_l) + cl*grad_xi_l)
      enddo
    enddo
  enddo
enddo

! ================================================================
! STEP 2: For each variable, each j-k line: split → reconstruct → diff
! ================================================================

do var = 1, nvars
  do nbl = 1, nblocks
    NI_l = NI(nbl)
    do k = 1, NK(nbl)
      do j = 1, NJ(nbl)

        ! ---- LF splitting along this i-line ----
        do i = 1, NI_l
          vol_l = 1.d0 / Jac(i,j,k,nbl)
          Fp(i) = 0.5d0*(PHI(i,j,k,nbl,var) &
                       + alpha_w * Qc(i,j,k,nbl,var) * vol_l)
          Fm(i) = 0.5d0*(PHI(i,j,k,nbl,var) &
                       - alpha_w * Qc(i,j,k,nbl,var) * vol_l)
        enddo

        ! ---- Reconstruct at each interface i+1/2 ----
        do i = 1, NI_l

          ! Periodic wrapping (NI coincides with 1)
          ic  = i
          im2 = i-2;  if (im2.lt.1)    im2 = im2 + NI_l - 1
          im1 = i-1;  if (im1.lt.1)    im1 = im1 + NI_l - 1
          ip1 = i+1;  if (ip1.gt.NI_l) ip1 = ip1 - NI_l + 1
          ip2 = i+2;  if (ip2.gt.NI_l) ip2 = ip2 - NI_l + 1
          ip3 = i+3;  if (ip3.gt.NI_l) ip3 = ip3 - NI_l + 1

          ! ===========================================================
          ! LEFT-BIASED RECONSTRUCTION OF F+
          ! ===========================================================

          ! Candidate fluxes
          f0 = (2.d0*Fp(im2) - 7.d0*Fp(im1) + 11.d0*Fp(ic)) / 6.d0
          f1 = (-Fp(im1) + 5.d0*Fp(ic) + 2.d0*Fp(ip1)) / 6.d0
          f2 = (2.d0*Fp(ic) + 5.d0*Fp(ip1) - Fp(ip2)) / 6.d0

          ! Smoothness indicators (Jiang-Shu)
          b0 = (13.d0/12.d0)*(Fp(im2)-2.d0*Fp(im1)+Fp(ic))**2 &
             + 0.25d0*(Fp(im2)-4.d0*Fp(im1)+3.d0*Fp(ic))**2
          b1 = (13.d0/12.d0)*(Fp(im1)-2.d0*Fp(ic)+Fp(ip1))**2 &
             + 0.25d0*(Fp(im1)-Fp(ip1))**2
          b2 = (13.d0/12.d0)*(Fp(ic)-2.d0*Fp(ip1)+Fp(ip2))**2 &
             + 0.25d0*(3.d0*Fp(ic)-4.d0*Fp(ip1)+Fp(ip2))**2

          ! ---------- Weights ----------

          if (dscheme.eq.5) then                    ! WENO-JS
            a0 = d0 / (eps_js + b0)**2
            a1 = d1 / (eps_js + b1)**2
            a2 = d2 / (eps_js + b2)**2
            asum = a0 + a1 + a2
            fhat_plus = (a0*f0 + a1*f1 + a2*f2) / asum

          elseif (dscheme.eq.6) then                ! WENO-Z
            tau5 = dabs(b0 - b2)
            a0 = d0 * (1.d0 + tau5/(eps_z + b0))
            a1 = d1 * (1.d0 + tau5/(eps_z + b1))
            a2 = d2 * (1.d0 + tau5/(eps_z + b2))
            asum = a0 + a1 + a2
            fhat_plus = (a0*f0 + a1*f1 + a2*f2) / asum

          elseif (is_cu6) then                      ! CU6 or CU6-M
            ! Extra downwind candidate
            f3 = (11.d0*Fp(ip1)-7.d0*Fp(ip2)+2.d0*Fp(ip3)) / 6.d0
            ! Full 6-point smoothness indicator
            b6 = (1.d0/120960.d0) * ( &
                271779.d0*Fp(im2)**2 &
              + Fp(im2)*(-2380800.d0*Fp(im1) + 4086352.d0*Fp(ic) &
                       -  3462252.d0*Fp(ip1) + 1458762.d0*Fp(ip2) &
                       -   245620.d0*Fp(ip3)) &
              + Fp(im1)*( 5653317.d0*Fp(im1) - 20427884.d0*Fp(ic) &
                       + 17905032.d0*Fp(ip1) -  7727988.d0*Fp(ip2) &
                       +  1325006.d0*Fp(ip3)) &
              + Fp(ic) *(19510972.d0*Fp(ic) - 35817664.d0*Fp(ip1) &
                       + 15929912.d0*Fp(ip2) -  2792660.d0*Fp(ip3)) &
              + Fp(ip1)*(17195652.d0*Fp(ip1) - 15880404.d0*Fp(ip2) &
                       +  2863984.d0*Fp(ip3)) &
              + Fp(ip2)*( 3824847.d0*Fp(ip2) -  1429976.d0*Fp(ip3)) &
              + 139633.d0*Fp(ip3)**2 )
            tau6 = b6 - (b0 + b2 + 4.d0*b1) / 6.d0
            a0 = d0 * (C_cu + tau6/(b0 + eps_z))**q_cu
            a1 = d1 * (C_cu + tau6/(b1 + eps_z))**q_cu
            a2 = d2 * (C_cu + tau6/(b2 + eps_z))**q_cu
            a3 = d3 * (C_cu + tau6/(b6 + eps_z))**q_cu
            asum = a0 + a1 + a2 + a3
            fhat_plus = (a0*f0 + a1*f1 + a2*f2 + a3*f3) / asum
          endif

          ! ===========================================================
          ! RIGHT-BIASED RECONSTRUCTION OF F-
          ! Mirror: stencil flipped around interface i+1/2
          !   "im2" role → ip3,  "im1" role → ip2,  "ic" role → ip1
          !   "ip1" role → ic,   "ip2" role → im1,  "ip3" role → im2
          ! ===========================================================

          ! Candidate fluxes (mirrored)
          g0 = (2.d0*Fm(ip3) - 7.d0*Fm(ip2) + 11.d0*Fm(ip1)) / 6.d0
          g1 = (-Fm(ip2) + 5.d0*Fm(ip1) + 2.d0*Fm(ic)) / 6.d0
          g2 = (2.d0*Fm(ip1) + 5.d0*Fm(ic) - Fm(im1)) / 6.d0

          ! Smoothness indicators (mirrored)
          b0 = (13.d0/12.d0)*(Fm(ip3)-2.d0*Fm(ip2)+Fm(ip1))**2 &
             + 0.25d0*(Fm(ip3)-4.d0*Fm(ip2)+3.d0*Fm(ip1))**2
          b1 = (13.d0/12.d0)*(Fm(ip2)-2.d0*Fm(ip1)+Fm(ic))**2 &
             + 0.25d0*(Fm(ip2)-Fm(ic))**2
          b2 = (13.d0/12.d0)*(Fm(ip1)-2.d0*Fm(ic)+Fm(im1))**2 &
             + 0.25d0*(3.d0*Fm(ip1)-4.d0*Fm(ic)+Fm(im1))**2

          ! ---------- Weights (mirrored, same formulas) ----------

          if (dscheme.eq.5) then
            a0 = d0 / (eps_js + b0)**2
            a1 = d1 / (eps_js + b1)**2
            a2 = d2 / (eps_js + b2)**2
            asum = a0 + a1 + a2
            fhat_minus = (a0*g0 + a1*g1 + a2*g2) / asum

          elseif (dscheme.eq.6) then
            tau5 = dabs(b0 - b2)
            a0 = d0 * (1.d0 + tau5/(eps_z + b0))
            a1 = d1 * (1.d0 + tau5/(eps_z + b1))
            a2 = d2 * (1.d0 + tau5/(eps_z + b2))
            asum = a0 + a1 + a2
            fhat_minus = (a0*g0 + a1*g1 + a2*g2) / asum

          elseif (is_cu6) then
            g3 = (11.d0*Fm(ic)-7.d0*Fm(im1)+2.d0*Fm(im2)) / 6.d0
            b6 = (1.d0/120960.d0) * ( &
                271779.d0*Fm(ip3)**2 &
              + Fm(ip3)*(-2380800.d0*Fm(ip2) + 4086352.d0*Fm(ip1) &
                       -  3462252.d0*Fm(ic)  + 1458762.d0*Fm(im1) &
                       -   245620.d0*Fm(im2)) &
              + Fm(ip2)*( 5653317.d0*Fm(ip2) - 20427884.d0*Fm(ip1) &
                       + 17905032.d0*Fm(ic)  -  7727988.d0*Fm(im1) &
                       +  1325006.d0*Fm(im2)) &
              + Fm(ip1)*(19510972.d0*Fm(ip1) - 35817664.d0*Fm(ic) &
                       + 15929912.d0*Fm(im1) -  2792660.d0*Fm(im2)) &
              + Fm(ic)  *(17195652.d0*Fm(ic) - 15880404.d0*Fm(im1) &
                       +  2863984.d0*Fm(im2)) &
              + Fm(im1)*( 3824847.d0*Fm(im1) -  1429976.d0*Fm(im2)) &
              + 139633.d0*Fm(im2)**2 )
            tau6 = b6 - (b0 + b2 + 4.d0*b1) / 6.d0
            a0 = d0 * (C_cu + tau6/(b0 + eps_z))**q_cu
            a1 = d1 * (C_cu + tau6/(b1 + eps_z))**q_cu
            a2 = d2 * (C_cu + tau6/(b2 + eps_z))**q_cu
            a3 = d3 * (C_cu + tau6/(b6 + eps_z))**q_cu
            asum = a0 + a1 + a2 + a3
            fhat_minus = (a0*g0 + a1*g1 + a2*g2 + a3*g3) / asum
          endif

          ! Total interface flux at i+1/2
          Fhat_total(i) = fhat_plus + fhat_minus

        enddo  ! end of interface loop (i)

        ! ---- Flux derivative along this line ----
        ! PHID = -(Fhat(i+1/2) - Fhat(i-1/2))
        ! Negative sign so caller just does: net_flux = net_flux + PHID
        do i = 1, NI_l
          im1f = i - 1;  if (im1f.lt.1) im1f = im1f + NI_l - 1

         !  PHID(i,j,k,nbl,var) = (Fhat_total(i) - Fhat_total(im1f))

          PHID(i,j,k,nbl,var) = -(Fhat_total(i) - Fhat_total(im1f))

        enddo

      enddo  ! j
    enddo    ! k
  enddo      ! nbl
enddo        ! var

END SUBROUTINE WENO_I





!=========================================================================
! WENO_J: WENO flux derivative in the J (eta) direction.
! Same structure as WENO_I, sweeps along j instead of i.
! Uses jx/jy/jz metrics and Vcont for alpha computation.
!=========================================================================

SUBROUTINE WENO_J(PHI, PHID, nvars)
use declare_variables
implicit none

integer, intent(in) :: nvars
real, dimension(NImax,NJmax,NKmax,nblocks,nvars), intent(in)  :: PHI
real, dimension(NImax,NJmax,NKmax,nblocks,nvars), intent(out) :: PHID

! --- 1D working arrays along one j-line ---
! real :: Fp(NJmax), Fm(NJmax)
! real :: Fhat_total(NJmax)

real(8) :: Fp(NJmax), Fm(NJmax), Fhat_total(NJmax)

! --- Scalars ---
real(8) :: alpha_w
real(8) :: cl, Vcont_l, grad_eta_l, vol_l
real(8) :: fhat_plus, fhat_minus
real(8) :: f0, f1, f2, f3
real(8) :: g0, g1, g2, g3
real(8) :: b0, b1, b2, b6
real(8) :: tau5, tau6
real(8) :: a0, a1, a2, a3, asum

integer :: var, jm2, jm1, jc, jp1, jp2, jp3, jm1f
integer :: NJ_l

real(8) :: eps_js, eps_z
real(8) :: d0, d1, d2, d3
real(8) :: C_cu, q_cu
logical :: is_cu6

! ================================================================
! Set scheme parameters
! ================================================================
eps_js = 1.d-6
eps_z  = 1.d-40

if (dscheme.eq.5 .or. dscheme.eq.6) then
  d0 = 1.d0/10.d0;  d1 = 6.d0/10.d0;  d2 = 3.d0/10.d0
  is_cu6 = .false.
else
  d0 = 1.d0/20.d0;  d1 = 9.d0/20.d0;  d2 = 9.d0/20.d0;  d3 = 1.d0/20.d0
  is_cu6 = .true.
  if (dscheme.eq.7) then
    C_cu = 20.d0;   q_cu = 1.d0
  else
    C_cu = 1000.d0;  q_cu = 4.d0
  endif
endif

! ================================================================
! STEP 1: Global max wave speed in J-direction
! ================================================================
alpha_w = 0.d0
do nbl = 1, nblocks
  do k = 1, NK(nbl)
    do j = 1, NJ(nbl)
      do i = 1, NI(nbl)
        cl = dsqrt(gamma * Qp(i,j,k,nbl,5) / Qp(i,j,k,nbl,1))
        Vcont_l = Qp(i,j,k,nbl,2)*jx(i,j,k,nbl) &
                + Qp(i,j,k,nbl,3)*jy(i,j,k,nbl) &
                + Qp(i,j,k,nbl,4)*jz(i,j,k,nbl)
        grad_eta_l = dsqrt(jx(i,j,k,nbl)**2 &
                         + jy(i,j,k,nbl)**2 &
                         + jz(i,j,k,nbl)**2)
        alpha_w = dmax1(alpha_w, dabs(Vcont_l) + cl*grad_eta_l)
      enddo
    enddo
  enddo
enddo

! ================================================================
! STEP 2: For each variable, each i-k line: split → reconstruct → diff
! ================================================================

do var = 1, nvars
  do nbl = 1, nblocks
    NJ_l = NJ(nbl)
    do k = 1, NK(nbl)
      do i = 1, NI(nbl)

        ! ---- LF splitting along this j-line ----
        do j = 1, NJ_l
          vol_l = 1.d0 / Jac(i,j,k,nbl)
          Fp(j) = 0.5d0*(PHI(i,j,k,nbl,var) &
                       + alpha_w * Qc(i,j,k,nbl,var) * vol_l)
          Fm(j) = 0.5d0*(PHI(i,j,k,nbl,var) &
                       - alpha_w * Qc(i,j,k,nbl,var) * vol_l)
        enddo

        ! ---- Reconstruct at each interface j+1/2 ----
        do j = 1, NJ_l

          jc  = j
          jm2 = j-2;  if (jm2.lt.1)    jm2 = jm2 + NJ_l - 1
          jm1 = j-1;  if (jm1.lt.1)    jm1 = jm1 + NJ_l - 1
          jp1 = j+1;  if (jp1.gt.NJ_l) jp1 = jp1 - NJ_l + 1
          jp2 = j+2;  if (jp2.gt.NJ_l) jp2 = jp2 - NJ_l + 1
          jp3 = j+3;  if (jp3.gt.NJ_l) jp3 = jp3 - NJ_l + 1

          ! ---- LEFT-BIASED RECONSTRUCTION OF F+ ----
          f0 = (2.d0*Fp(jm2) - 7.d0*Fp(jm1) + 11.d0*Fp(jc)) / 6.d0
          f1 = (-Fp(jm1) + 5.d0*Fp(jc) + 2.d0*Fp(jp1)) / 6.d0
          f2 = (2.d0*Fp(jc) + 5.d0*Fp(jp1) - Fp(jp2)) / 6.d0

          b0 = (13.d0/12.d0)*(Fp(jm2)-2.d0*Fp(jm1)+Fp(jc))**2 &
             + 0.25d0*(Fp(jm2)-4.d0*Fp(jm1)+3.d0*Fp(jc))**2
          b1 = (13.d0/12.d0)*(Fp(jm1)-2.d0*Fp(jc)+Fp(jp1))**2 &
             + 0.25d0*(Fp(jm1)-Fp(jp1))**2
          b2 = (13.d0/12.d0)*(Fp(jc)-2.d0*Fp(jp1)+Fp(jp2))**2 &
             + 0.25d0*(3.d0*Fp(jc)-4.d0*Fp(jp1)+Fp(jp2))**2

          if (dscheme.eq.5) then
            a0 = d0 / (eps_js + b0)**2
            a1 = d1 / (eps_js + b1)**2
            a2 = d2 / (eps_js + b2)**2
            asum = a0 + a1 + a2
            fhat_plus = (a0*f0 + a1*f1 + a2*f2) / asum

          elseif (dscheme.eq.6) then
            tau5 = dabs(b0 - b2)
            a0 = d0 * (1.d0 + tau5/(eps_z + b0))
            a1 = d1 * (1.d0 + tau5/(eps_z + b1))
            a2 = d2 * (1.d0 + tau5/(eps_z + b2))
            asum = a0 + a1 + a2
            fhat_plus = (a0*f0 + a1*f1 + a2*f2) / asum

          elseif (is_cu6) then
            f3 = (11.d0*Fp(jp1)-7.d0*Fp(jp2)+2.d0*Fp(jp3)) / 6.d0
            b6 = (1.d0/120960.d0) * ( &
                271779.d0*Fp(jm2)**2 &
              + Fp(jm2)*(-2380800.d0*Fp(jm1) + 4086352.d0*Fp(jc) &
                       -  3462252.d0*Fp(jp1) + 1458762.d0*Fp(jp2) &
                       -   245620.d0*Fp(jp3)) &
              + Fp(jm1)*( 5653317.d0*Fp(jm1) - 20427884.d0*Fp(jc) &
                       + 17905032.d0*Fp(jp1) -  7727988.d0*Fp(jp2) &
                       +  1325006.d0*Fp(jp3)) &
              + Fp(jc) *(19510972.d0*Fp(jc) - 35817664.d0*Fp(jp1) &
                       + 15929912.d0*Fp(jp2) -  2792660.d0*Fp(jp3)) &
              + Fp(jp1)*(17195652.d0*Fp(jp1) - 15880404.d0*Fp(jp2) &
                       +  2863984.d0*Fp(jp3)) &
              + Fp(jp2)*( 3824847.d0*Fp(jp2) -  1429976.d0*Fp(jp3)) &
              + 139633.d0*Fp(jp3)**2 )
            tau6 = b6 - (b0 + b2 + 4.d0*b1) / 6.d0
            a0 = d0 * (C_cu + tau6/(b0 + eps_z))**q_cu
            a1 = d1 * (C_cu + tau6/(b1 + eps_z))**q_cu
            a2 = d2 * (C_cu + tau6/(b2 + eps_z))**q_cu
            a3 = d3 * (C_cu + tau6/(b6 + eps_z))**q_cu
            asum = a0 + a1 + a2 + a3
            fhat_plus = (a0*f0 + a1*f1 + a2*f2 + a3*f3) / asum
          endif

          ! ---- RIGHT-BIASED RECONSTRUCTION OF F- ----
          g0 = (2.d0*Fm(jp3) - 7.d0*Fm(jp2) + 11.d0*Fm(jp1)) / 6.d0
          g1 = (-Fm(jp2) + 5.d0*Fm(jp1) + 2.d0*Fm(jc)) / 6.d0
          g2 = (2.d0*Fm(jp1) + 5.d0*Fm(jc) - Fm(jm1)) / 6.d0

          b0 = (13.d0/12.d0)*(Fm(jp3)-2.d0*Fm(jp2)+Fm(jp1))**2 &
             + 0.25d0*(Fm(jp3)-4.d0*Fm(jp2)+3.d0*Fm(jp1))**2
          b1 = (13.d0/12.d0)*(Fm(jp2)-2.d0*Fm(jp1)+Fm(jc))**2 &
             + 0.25d0*(Fm(jp2)-Fm(jc))**2
          b2 = (13.d0/12.d0)*(Fm(jp1)-2.d0*Fm(jc)+Fm(jm1))**2 &
             + 0.25d0*(3.d0*Fm(jp1)-4.d0*Fm(jc)+Fm(jm1))**2

          if (dscheme.eq.5) then
            a0 = d0 / (eps_js + b0)**2
            a1 = d1 / (eps_js + b1)**2
            a2 = d2 / (eps_js + b2)**2
            asum = a0 + a1 + a2
            fhat_minus = (a0*g0 + a1*g1 + a2*g2) / asum

          elseif (dscheme.eq.6) then
            tau5 = dabs(b0 - b2)
            a0 = d0 * (1.d0 + tau5/(eps_z + b0))
            a1 = d1 * (1.d0 + tau5/(eps_z + b1))
            a2 = d2 * (1.d0 + tau5/(eps_z + b2))
            asum = a0 + a1 + a2
            fhat_minus = (a0*g0 + a1*g1 + a2*g2) / asum

          elseif (is_cu6) then
            g3 = (11.d0*Fm(jc)-7.d0*Fm(jm1)+2.d0*Fm(jm2)) / 6.d0
            b6 = (1.d0/120960.d0) * ( &
                271779.d0*Fm(jp3)**2 &
              + Fm(jp3)*(-2380800.d0*Fm(jp2) + 4086352.d0*Fm(jp1) &
                       -  3462252.d0*Fm(jc)  + 1458762.d0*Fm(jm1) &
                       -   245620.d0*Fm(jm2)) &
              + Fm(jp2)*( 5653317.d0*Fm(jp2) - 20427884.d0*Fm(jp1) &
                       + 17905032.d0*Fm(jc)  -  7727988.d0*Fm(jm1) &
                       +  1325006.d0*Fm(jm2)) &
              + Fm(jp1)*(19510972.d0*Fm(jp1) - 35817664.d0*Fm(jc) &
                       + 15929912.d0*Fm(jm1) -  2792660.d0*Fm(jm2)) &
              + Fm(jc)  *(17195652.d0*Fm(jc) - 15880404.d0*Fm(jm1) &
                       +  2863984.d0*Fm(jm2)) &
              + Fm(jm1)*( 3824847.d0*Fm(jm1) -  1429976.d0*Fm(jm2)) &
              + 139633.d0*Fm(jm2)**2 )
            tau6 = b6 - (b0 + b2 + 4.d0*b1) / 6.d0
            a0 = d0 * (C_cu + tau6/(b0 + eps_z))**q_cu
            a1 = d1 * (C_cu + tau6/(b1 + eps_z))**q_cu
            a2 = d2 * (C_cu + tau6/(b2 + eps_z))**q_cu
            a3 = d3 * (C_cu + tau6/(b6 + eps_z))**q_cu
            asum = a0 + a1 + a2 + a3
            fhat_minus = (a0*g0 + a1*g1 + a2*g2 + a3*g3) / asum
          endif

          Fhat_total(j) = fhat_plus + fhat_minus

        enddo  ! j reconstruction

        ! ---- Flux derivative ----
        do j = 1, NJ_l
          jm1f = j - 1;  if (jm1f.lt.1) jm1f = jm1f + NJ_l - 1
          PHID(i,j,k,nbl,var) = -(Fhat_total(j) - Fhat_total(jm1f))
        enddo

      enddo  ! i
    enddo    ! k
  enddo      ! nbl
enddo        ! var

END SUBROUTINE WENO_J





!=========================================================================
! WENO_K: WENO flux derivative in the K (zeta) direction.
! Same structure as WENO_I, sweeps along k instead of i.
! Uses kx/ky/kz metrics and Wcont for alpha computation.
!=========================================================================

SUBROUTINE WENO_K(PHI, PHID, nvars)
use declare_variables
implicit none

integer, intent(in) :: nvars
real, dimension(NImax,NJmax,NKmax,nblocks,nvars), intent(in)  :: PHI
real, dimension(NImax,NJmax,NKmax,nblocks,nvars), intent(out) :: PHID

! --- 1D working arrays along one k-line ---
! real :: Fp(NKmax), Fm(NKmax)
! real :: Fhat_total(NKmax)

real(8) :: Fp(NKmax), Fm(NKmax), Fhat_total(NKmax)

! --- Scalars ---
real(8) :: alpha_w
real(8) :: cl, Wcont_l, grad_zeta_l, vol_l
real(8) :: fhat_plus, fhat_minus
real(8) :: f0, f1, f2, f3
real(8) :: g0, g1, g2, g3
real(8) :: b0, b1, b2, b6
real(8) :: tau5, tau6
real(8) :: a0, a1, a2, a3, asum

integer :: var, km2, km1, kc, kp1, kp2, kp3, km1f
integer :: NK_l

real(8) :: eps_js, eps_z
real(8) :: d0, d1, d2, d3
real(8) :: C_cu, q_cu
logical :: is_cu6

! ================================================================
! Set scheme parameters
! ================================================================
eps_js = 1.d-6
eps_z  = 1.d-40

if (dscheme.eq.5 .or. dscheme.eq.6) then
  d0 = 1.d0/10.d0;  d1 = 6.d0/10.d0;  d2 = 3.d0/10.d0
  is_cu6 = .false.
else
  d0 = 1.d0/20.d0;  d1 = 9.d0/20.d0;  d2 = 9.d0/20.d0;  d3 = 1.d0/20.d0
  is_cu6 = .true.
  if (dscheme.eq.7) then
    C_cu = 20.d0;   q_cu = 1.d0
  else
    C_cu = 1000.d0;  q_cu = 4.d0
  endif
endif

! ================================================================
! STEP 1: Global max wave speed in K-direction
! ================================================================
alpha_w = 0.d0
do nbl = 1, nblocks
  do k = 1, NK(nbl)
    do j = 1, NJ(nbl)
      do i = 1, NI(nbl)
        cl = dsqrt(gamma * Qp(i,j,k,nbl,5) / Qp(i,j,k,nbl,1))
        Wcont_l = Qp(i,j,k,nbl,2)*kx(i,j,k,nbl) &
                + Qp(i,j,k,nbl,3)*ky(i,j,k,nbl) &
                + Qp(i,j,k,nbl,4)*kz(i,j,k,nbl)
        grad_zeta_l = dsqrt(kx(i,j,k,nbl)**2 &
                          + ky(i,j,k,nbl)**2 &
                          + kz(i,j,k,nbl)**2)
        alpha_w = dmax1(alpha_w, dabs(Wcont_l) + cl*grad_zeta_l)
      enddo
    enddo
  enddo
enddo

! ================================================================
! STEP 2: For each variable, each i-j line: split → reconstruct → diff
! ================================================================

do var = 1, nvars
  do nbl = 1, nblocks
    NK_l = NK(nbl)
    do j = 1, NJ(nbl)
      do i = 1, NI(nbl)

        ! ---- LF splitting along this k-line ----
        do k = 1, NK_l
          vol_l = 1.d0 / Jac(i,j,k,nbl)
          Fp(k) = 0.5d0*(PHI(i,j,k,nbl,var) &
                       + alpha_w * Qc(i,j,k,nbl,var) * vol_l)
          Fm(k) = 0.5d0*(PHI(i,j,k,nbl,var) &
                       - alpha_w * Qc(i,j,k,nbl,var) * vol_l)
        enddo

        ! ---- Reconstruct at each interface k+1/2 ----
        do k = 1, NK_l

          kc  = k
          km2 = k-2;  if (km2.lt.1)    km2 = km2 + NK_l - 1
          km1 = k-1;  if (km1.lt.1)    km1 = km1 + NK_l - 1
          kp1 = k+1;  if (kp1.gt.NK_l) kp1 = kp1 - NK_l + 1
          kp2 = k+2;  if (kp2.gt.NK_l) kp2 = kp2 - NK_l + 1
          kp3 = k+3;  if (kp3.gt.NK_l) kp3 = kp3 - NK_l + 1

          ! ---- LEFT-BIASED RECONSTRUCTION OF F+ ----
          f0 = (2.d0*Fp(km2) - 7.d0*Fp(km1) + 11.d0*Fp(kc)) / 6.d0
          f1 = (-Fp(km1) + 5.d0*Fp(kc) + 2.d0*Fp(kp1)) / 6.d0
          f2 = (2.d0*Fp(kc) + 5.d0*Fp(kp1) - Fp(kp2)) / 6.d0

          b0 = (13.d0/12.d0)*(Fp(km2)-2.d0*Fp(km1)+Fp(kc))**2 &
             + 0.25d0*(Fp(km2)-4.d0*Fp(km1)+3.d0*Fp(kc))**2
          b1 = (13.d0/12.d0)*(Fp(km1)-2.d0*Fp(kc)+Fp(kp1))**2 &
             + 0.25d0*(Fp(km1)-Fp(kp1))**2
          b2 = (13.d0/12.d0)*(Fp(kc)-2.d0*Fp(kp1)+Fp(kp2))**2 &
             + 0.25d0*(3.d0*Fp(kc)-4.d0*Fp(kp1)+Fp(kp2))**2

          if (dscheme.eq.5) then
            a0 = d0 / (eps_js + b0)**2
            a1 = d1 / (eps_js + b1)**2
            a2 = d2 / (eps_js + b2)**2
            asum = a0 + a1 + a2
            fhat_plus = (a0*f0 + a1*f1 + a2*f2) / asum

          elseif (dscheme.eq.6) then
            tau5 = dabs(b0 - b2)
            a0 = d0 * (1.d0 + tau5/(eps_z + b0))
            a1 = d1 * (1.d0 + tau5/(eps_z + b1))
            a2 = d2 * (1.d0 + tau5/(eps_z + b2))
            asum = a0 + a1 + a2
            fhat_plus = (a0*f0 + a1*f1 + a2*f2) / asum

          elseif (is_cu6) then
            f3 = (11.d0*Fp(kp1)-7.d0*Fp(kp2)+2.d0*Fp(kp3)) / 6.d0
            b6 = (1.d0/120960.d0) * ( &
                271779.d0*Fp(km2)**2 &
              + Fp(km2)*(-2380800.d0*Fp(km1) + 4086352.d0*Fp(kc) &
                       -  3462252.d0*Fp(kp1) + 1458762.d0*Fp(kp2) &
                       -   245620.d0*Fp(kp3)) &
              + Fp(km1)*( 5653317.d0*Fp(km1) - 20427884.d0*Fp(kc) &
                       + 17905032.d0*Fp(kp1) -  7727988.d0*Fp(kp2) &
                       +  1325006.d0*Fp(kp3)) &
              + Fp(kc) *(19510972.d0*Fp(kc) - 35817664.d0*Fp(kp1) &
                       + 15929912.d0*Fp(kp2) -  2792660.d0*Fp(kp3)) &
              + Fp(kp1)*(17195652.d0*Fp(kp1) - 15880404.d0*Fp(kp2) &
                       +  2863984.d0*Fp(kp3)) &
              + Fp(kp2)*( 3824847.d0*Fp(kp2) -  1429976.d0*Fp(kp3)) &
              + 139633.d0*Fp(kp3)**2 )
            tau6 = b6 - (b0 + b2 + 4.d0*b1) / 6.d0
            a0 = d0 * (C_cu + tau6/(b0 + eps_z))**q_cu
            a1 = d1 * (C_cu + tau6/(b1 + eps_z))**q_cu
            a2 = d2 * (C_cu + tau6/(b2 + eps_z))**q_cu
            a3 = d3 * (C_cu + tau6/(b6 + eps_z))**q_cu
            asum = a0 + a1 + a2 + a3
            fhat_plus = (a0*f0 + a1*f1 + a2*f2 + a3*f3) / asum
          endif

          ! ---- RIGHT-BIASED RECONSTRUCTION OF F- ----
          g0 = (2.d0*Fm(kp3) - 7.d0*Fm(kp2) + 11.d0*Fm(kp1)) / 6.d0
          g1 = (-Fm(kp2) + 5.d0*Fm(kp1) + 2.d0*Fm(kc)) / 6.d0
          g2 = (2.d0*Fm(kp1) + 5.d0*Fm(kc) - Fm(km1)) / 6.d0

          b0 = (13.d0/12.d0)*(Fm(kp3)-2.d0*Fm(kp2)+Fm(kp1))**2 &
             + 0.25d0*(Fm(kp3)-4.d0*Fm(kp2)+3.d0*Fm(kp1))**2
          b1 = (13.d0/12.d0)*(Fm(kp2)-2.d0*Fm(kp1)+Fm(kc))**2 &
             + 0.25d0*(Fm(kp2)-Fm(kc))**2
          b2 = (13.d0/12.d0)*(Fm(kp1)-2.d0*Fm(kc)+Fm(km1))**2 &
             + 0.25d0*(3.d0*Fm(kp1)-4.d0*Fm(kc)+Fm(km1))**2

          if (dscheme.eq.5) then
            a0 = d0 / (eps_js + b0)**2
            a1 = d1 / (eps_js + b1)**2
            a2 = d2 / (eps_js + b2)**2
            asum = a0 + a1 + a2
            fhat_minus = (a0*g0 + a1*g1 + a2*g2) / asum

          elseif (dscheme.eq.6) then
            tau5 = dabs(b0 - b2)
            a0 = d0 * (1.d0 + tau5/(eps_z + b0))
            a1 = d1 * (1.d0 + tau5/(eps_z + b1))
            a2 = d2 * (1.d0 + tau5/(eps_z + b2))
            asum = a0 + a1 + a2
            fhat_minus = (a0*g0 + a1*g1 + a2*g2) / asum

          elseif (is_cu6) then
            g3 = (11.d0*Fm(kc)-7.d0*Fm(km1)+2.d0*Fm(km2)) / 6.d0
            b6 = (1.d0/120960.d0) * ( &
                271779.d0*Fm(kp3)**2 &
              + Fm(kp3)*(-2380800.d0*Fm(kp2) + 4086352.d0*Fm(kp1) &
                       -  3462252.d0*Fm(kc)  + 1458762.d0*Fm(km1) &
                       -   245620.d0*Fm(km2)) &
              + Fm(kp2)*( 5653317.d0*Fm(kp2) - 20427884.d0*Fm(kp1) &
                       + 17905032.d0*Fm(kc)  -  7727988.d0*Fm(km1) &
                       +  1325006.d0*Fm(km2)) &
              + Fm(kp1)*(19510972.d0*Fm(kp1) - 35817664.d0*Fm(kc) &
                       + 15929912.d0*Fm(km1) -  2792660.d0*Fm(km2)) &
              + Fm(kc)  *(17195652.d0*Fm(kc) - 15880404.d0*Fm(km1) &
                       +  2863984.d0*Fm(km2)) &
              + Fm(km1)*( 3824847.d0*Fm(km1) -  1429976.d0*Fm(km2)) &
              + 139633.d0*Fm(km2)**2 )
            tau6 = b6 - (b0 + b2 + 4.d0*b1) / 6.d0
            a0 = d0 * (C_cu + tau6/(b0 + eps_z))**q_cu
            a1 = d1 * (C_cu + tau6/(b1 + eps_z))**q_cu
            a2 = d2 * (C_cu + tau6/(b2 + eps_z))**q_cu
            a3 = d3 * (C_cu + tau6/(b6 + eps_z))**q_cu
            asum = a0 + a1 + a2 + a3
            fhat_minus = (a0*g0 + a1*g1 + a2*g2 + a3*g3) / asum
          endif

          Fhat_total(k) = fhat_plus + fhat_minus

        enddo  ! k reconstruction

        ! ---- Flux derivative ----
        do k = 1, NK_l
          km1f = k - 1;  if (km1f.lt.1) km1f = km1f + NK_l - 1
          PHID(i,j,k,nbl,var) = -(Fhat_total(k) - Fhat_total(km1f))
        enddo

      enddo  ! i
    enddo    ! j
  enddo      ! nbl
enddo        ! var

END SUBROUTINE WENO_K



	  
!******************************************************************************************************************	  
!***********************************PERIODIC DISCRETIZATION ROUTINES FOR GRID COORDINATES *************************
!******************************************************************************************************************

      ! implicit

      ! ------------------------------------------------------------------------------------

      ! THE LLX AND LLY AND LLZ ADDITIONS HAVE MISTAKES HERE. REDO THE MATH ONCE TO VERIFY!!!!

      ! ------------------------------------------------------------------------------------

      SUBROUTINE DISCRETIZATION_I_EXP_GRID(PHI,PHID)
      use declare_variables	 
      implicit none	  
	  
      integer var,ip1,ip2,im1,im2
      real bb4, ab2, LLx  
      real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID

      bb4 = bdisc/4.d0
      ab2 = adisc/2.d0

      ! okay, the professor seems to do the order first with i as the inner loop (as we have here),
      ! then with j as the inner loop, then with k as the inner loop.
      ! why?
      
      ! where are we slapping in the boundary conds for the other directions?

      ! we can do that but we have to copy pasta 

      ! this specific loop is for the i direction only.

      Do nbl = 1,nblocks


         Do k = 1,NK(nbl)


            Do j = 1,NJ(nbl)



               Do i = 3,NI(nbl)-2

                  ip2 = i+2
                  ip1 = i+1  
                  im1 = i-1
                  im2 = i-2
                  ! PHID(i,j,k,nbl) = (bb4*PHI(ip2,j,k,nbl) + ab2*PHI(ip1,j,k,nbl) - ab2*PHI(im1,j,k,nbl) - bb4*PHI(im2,j,k,nbl))
                  PHID(i,j,k,nbl) = bb4*(PHI(ip2,j,k,nbl) - PHI(im2,j,k,nbl)) + ab2*(PHI(ip1,j,k,nbl) - PHI(im1,j,k,nbl))

               enddo

               ! Boundary points defined separately so that it dont fuck up.

               i=1
               ip2 = i+2
               ip1 = i+1  
               im1 = NI(nbl)-1
               im2 = NI(nbl)-2

               LLx = PHI(NI(nbl),j,k,nbl) - PHI(1,j,k,nbl)

               PHID(i,j,k,nbl) = bb4*(PHI(ip2,j,k,nbl) - PHI(im2,j,k,nbl) + LLx)+ ab2*(PHI(ip1,j,k,nbl) - PHI(im1,j,k,nbl) + LLx)

               i=2
               ip2 = i+2  
               ip1 = i+1 
               im1 = i-1
               im2 = NI(nbl)-1

               ! check if the llx implementation is correct here

               PHID(i,j,k,nbl) = bb4*(PHI(ip2,j,k,nbl) - PHI(im2,j,k,nbl) + LLx)+ ab2*(PHI(ip1,j,k,nbl) - PHI(im1,j,k,nbl))

               i=NI(nbl)-1
               ip1 = i+1
               ip2 = 2
               im1 = i-1
               im2 = i-2

               PHID(i,j,k,nbl) = bb4*(PHI(ip2,j,k,nbl) - PHI(im2,j,k,nbl) + LLx)+ ab2*(PHI(ip1,j,k,nbl) - PHI(im1,j,k,nbl))

               i=NI(nbl)
               ip2 = 3
               ip1 = 2
               im1 = i-1
               im2 = i-2

               PHID(i,j,k,nbl) = bb4*(PHI(ip2,j,k,nbl) - PHI(im2,j,k,nbl) + LLx)+ ab2*(PHI(ip1,j,k,nbl) - PHI(im1,j,k,nbl) + LLx)

               ! COMPLETE THEM



            enddo
         enddo
      enddo





      ! copypasta the same loops for the j and k directions. 
      ! make sure to edit the other codes that wew have to do before we do this.






      END	






      ! check the ones below!!!!!!




      SUBROUTINE DISCRETIZATION_J_EXP_GRID(PHI,PHID)
      use declare_variables	 
      implicit none	 

      integer var,jp1,jp2,jm1,jm2
      real bb4, ab2, LLy 
      real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID

      bb4 = bdisc/4.d0
      ab2 = adisc/2.d0

      

      ! this specific loop is for the j direction only.

      Do nbl = 1,nblocks


         Do k = 1,NK(nbl)


            Do i = 1,NI(nbl)



               Do j = 3,NJ(nbl)-2

                  jp2 = j+2
                  jp1 = j+1  
                  jm1 = j-1
                  jm2 = j-2
                  ! PHID(i,j,k,nbl) = (bb4*PHI(i,jp2,k,nbl) + ab2*PHI(i,jp1,k,nbl) - ab2*PHI(i,jm1,k,nbl) - bb4*PHI(i,jm2,k,nbl))
                  PHID(i,j,k,nbl) = bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl)) + ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl))
                  
               enddo

               ! Boundary points defined separately so that it dont fuck up.

               j=1
               jp2 = j+2
               jp1 = j+1  
               jm1 = NJ(nbl)-1
               jm2 = NJ(nbl)-2

               LLy = PHI(i,NJ(nbl),k,nbl) - PHI(i,1,k,nbl)

               PHID(i,j,k,nbl) = bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl) + LLy)+ ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl) + LLy)

               j=2
               jp2 = j+2  
               jp1 = j+1 
               jm1 = j-1
               jm2 = NJ(nbl)-1

               ! check if the llx implementation is correct here

               PHID(i,j,k,nbl) = bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl) + LLy)+ ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl))

               j=NJ(nbl)-1
               jp1 = j+1
               jp2 = 2
               jm1 = j-1
               jm2 = j-2

               PHID(i,j,k,nbl) = bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl) + LLy)+ ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl))

               j=NJ(nbl)
               jp2 = 3
               jp1 = 2
               jm1 = j-1
               jm2 = j-2

               PHID(i,j,k,nbl) = (bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl) + LLy)+ ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl) + LLy))

               ! COMPLETE THEM



            enddo
         enddo
      enddo


      END	



      SUBROUTINE DISCRETIZATION_K_EXP_GRID(PHI,PHID)
      use declare_variables	 
      implicit none	  
	  
      integer var,kp1,kp2,km1,km2
      real bb4, ab2, LLz
      real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID	  

      bb4 = bdisc/4.d0
      ab2 = adisc/2.d0

      

      ! this specific loop is for the k direction only.

      Do nbl = 1,nblocks


         Do i = 1,NI(nbl)


            Do j = 1,NJ(nbl)



               Do k = 3,NK(nbl)-2

                  kp2 = k+2
                  kp1 = k+1  
                  km1 = k-1
                  km2 = k-2
                  ! PHID(i,j,k,nbl) = (bb4*PHI(i,j,kp2,nbl) + ab2*PHI(i,j,kp1,nbl) - ab2*PHI(i,j,km1,nbl) - bb4*PHI(i,j,km2,nbl))
                  PHID(i,j,k,nbl) = bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl)) + ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl))

               enddo

               ! Boundary points defined separately so that it dont fuck up.

               k=1
               kp2 = k+2
               kp1 = k+1  
               km1 = NK(nbl)-1
               km2 = NK(nbl)-2

               LLz = PHI(i,j,NK(nbl),nbl) - PHI(i,j,1,nbl)

               PHID(i,j,k,nbl) = bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl) + LLz)+ ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl) + LLz)

               k=2
               kp2 = k+2  
               kp1 = k+1 
               km1 = k-1
               km2 = NK(nbl)-1

               ! check if the llx implementation is correct here

               PHID(i,j,k,nbl) = bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl) + LLz) + ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl))

               k=NK(nbl)-1
               kp1 = k+1
               kp2 = 2
               km1 = k-1
               km2 = k-2

               PHID(i,j,k,nbl) = bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl) + LLz)+ ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl))

               k=NK(nbl)
               kp2 = 3
               kp1 = 2
               km1 = k-1
               km2 = k-2

               PHID(i,j,k,nbl) = bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl) + LLz)+ ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl) + LLz)

               ! COMPLETE THEM



            enddo
         enddo
      enddo

      END	

      SUBROUTINE DISCRETIZATION_K2D_EXP_GRID(PHI,PHID)
      use declare_variables	 
      implicit none	  
	  
      integer var,kp1,km1
      real bb4, ab2	  
      real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
  

      END	






















!****************************** COMPACT SCHEMES FOR GRID *******************************************************************


      ! ------------------------------------------------------------------------------------

      ! THE LLX AND LLY AND LLZ ADDITIONS HAVE MISTAKES HERE. REDO THE MATH ONCE TO VERIFY!!!!

      ! ------------------------------------------------------------------------------------


      ! implicit

      SUBROUTINE DISCRETIZATION_I_COMP_GRID(PHI,PHID)
      use declare_variables	 
      implicit none	  
	  
      integer var,ip1,ip2,im1,im2, nm1
      real bb4, ab2, LLx	  
      real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
      real,dimension(NImax) :: RHS	  

      bb4 = bdisc/4.d0
      ab2 = adisc/2.d0

      ! copied it from the loops above; need to mod it for the compact scheme.

      ! this specific loop is for the i direction only.

      Do nbl = 1,nblocks


         Do k = 1,NK(nbl)


            Do j = 1,NJ(nbl)



               Do i = 3,NI(nbl)-2

                  ip2 = i+2
                  ip1 = i+1  
                  im1 = i-1
                  im2 = i-2
                  ! RHS(i) = (bb4*PHI(ip2,j,k,nbl) + ab2*PHI(ip1,j,k,nbl) - ab2*PHI(im1,j,k,nbl) - bb4*PHI(im2,j,k,nbl))
                  RHS(i) = bb4*(PHI(ip2,j,k,nbl) - PHI(im2,j,k,nbl)) + ab2*(PHI(ip1,j,k,nbl) - PHI(im1,j,k,nbl))

               enddo

               ! Boundary points defined separately so that it dont fuck up.

               i=1
               ip2 = i+2
               ip1 = i+1  
               im1 = NI(nbl)-1
               im2 = NI(nbl)-2

               LLx = PHI(NI(nbl),j,k,nbl) - PHI(1,j,k,nbl)

               RHS(i) = bb4*(PHI(ip2,j,k,nbl) - PHI(im2,j,k,nbl) + LLx)+ ab2*(PHI(ip1,j,k,nbl) - PHI(im1,j,k,nbl) + LLx)

               i=2
               ip2 = i+2  
               ip1 = i+1 
               im1 = i-1
               im2 = NI(nbl)-1

               ! check if the llx implementation is correct here

               RHS(i) = bb4*(PHI(ip2,j,k,nbl) - PHI(im2,j,k,nbl) + LLx) + ab2*(PHI(ip1,j,k,nbl) - PHI(im1,j,k,nbl))

               i=NI(nbl)-1
               ip1 = i+1
               ip2 = 2
               im1 = i-1
               im2 = i-2

               RHS(i) = bb4*(PHI(ip2,j,k,nbl) - PHI(im2,j,k,nbl) + LLx)+ ab2*(PHI(ip1,j,k,nbl) - PHI(im1,j,k,nbl))

               i=NI(nbl)
               ip2 = 3
               ip1 = 2
               im1 = i-1
               im2 = i-2

               RHS(i) = bb4*(PHI(ip2,j,k,nbl) - PHI(im2,j,k,nbl) + LLx)+ ab2*(PHI(ip1,j,k,nbl) - PHI(im1,j,k,nbl) + LLx)

               ! COMPLETE THEM

               ! call TDMAP(ji,jf,ap,ac,am,fi,NMAXL)


               nm1 = NI(nbl)-1
               ! print*, RHS(1:nm1)
               call TDMAP(1,nm1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)

               PHID(1:nm1,j,k,nbl) = RHS(1:nm1)

               PHID(NI(nbl),j,k,nbl) = PHID(1,j,k,nbl) ! periodic BC; first and last point are the same

               ! print*, RHS(1:nm1)

               !ap - super diagonal
               !ac - diagonal
               !am - sub diagonal
               !ji - rhs

               ! the fucntion overwrites the RHS with the derivatives

               ! pause


            enddo
         enddo
      enddo


      END	

      SUBROUTINE DISCRETIZATION_J_COMP_GRID(PHI,PHID)
      use declare_variables	 
      implicit none	 

      integer var,jp1,jp2,jm1,jm2, nm1
      real bb4, ab2, LLy	  
      real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
      real,dimension(NJmax) :: RHS	  

      bb4 = bdisc/4.d0
      ab2 = adisc/2.d0

      

      ! this specific loop is for the j direction only.

      Do nbl = 1,nblocks


         Do k = 1,NK(nbl)


            Do i = 1,NI(nbl)



               Do j = 3,NJ(nbl)-2

                  jp2 = j+2
                  jp1 = j+1  
                  jm1 = j-1
                  jm2 = j-2
                  ! PHID(i,j,k,nbl) = (bb4*PHI(i,jp2,k,nbl) + ab2*PHI(i,jp1,k,nbl) - ab2*PHI(i,jm1,k,nbl) - bb4*PHI(i,jm2,k,nbl))
                  RHS(j) = bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl)) + ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl))
               enddo

               ! Boundary points defined separately so that it dont fuck up.

               j=1
               jp2 = j+2
               jp1 = j+1  
               jm1 = NJ(nbl)-1
               jm2 = NJ(nbl)-2

               LLy = PHI(i,NJ(nbl),k,nbl) - PHI(i,1,k,nbl)

               ! PHID(i,j,k,nbl) = (bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl) - LLy)+ ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl) - LLy))
               RHS(j) = bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl)+ LLy) + ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl)+ LLy)


               j=2
               jp2 = j+2  
               jp1 = j+1 
               jm1 = j-1
               jm2 = NJ(nbl)-1

               ! check if the llx implementation is correct here

               ! PHID(i,j,k,nbl) = (bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl))+ ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl)))
               RHS(j) = bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl)+ LLy) + ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl))
	  


               j=NJ(nbl)-1
               jp1 = j+1
               jp2 = 2
               jm1 = j-1
               jm2 = j-2

               ! PHID(i,j,k,nbl) = (bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl) + LLy)+ ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl)))
               RHS(j) = bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl)+ LLy) + ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl))


               j=NJ(nbl)
               jp2 = 3
               jp1 = 2
               jm1 = j-1
               jm2 = j-2

               ! PHID(i,j,k,nbl) = (bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl) + LLy)+ ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl) + LLy))
               RHS(j) = bb4*(PHI(i,jp2,k,nbl) - PHI(i,jm2,k,nbl)+ LLy) + ab2*(PHI(i,jp1,k,nbl) - PHI(i,jm1,k,nbl)+ LLy)


               

               nm1 = NJ(nbl)-1
               ! print*, RHS(1:nm1)
               call TDMAP(1,nm1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)

               PHID(i,1:nm1,k,nbl) = RHS(1:nm1)

               PHID(i,NJ(nbl),k,nbl) = PHID(i,1,k,nbl) ! periodic BC; first and last point are the same

               ! print*, RHS(1:nm1)

               ! pause



            enddo
         enddo
      enddo
      




      END	

      SUBROUTINE DISCRETIZATION_K_COMP_GRID(PHI,PHID)
      use declare_variables	 
      implicit none	  
	  
      integer var,kp1,kp2,km1,km2,nvars, nm1
      real bb4, ab2, LLz
      real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
      real,dimension(NKmax) :: RHS	  


      bb4 = bdisc/4.d0
      ab2 = adisc/2.d0

      

      ! this specific loop is for the k direction only.

      Do nbl = 1,nblocks


         Do i = 1,NI(nbl)


            Do j = 1,NJ(nbl)



               Do k = 3,NK(nbl)-2

                  kp2 = k+2
                  kp1 = k+1  
                  km1 = k-1
                  km2 = k-2
                  ! PHID(i,j,k,nbl) = (bb4*PHI(i,j,kp2,nbl) + ab2*PHI(i,j,kp1,nbl) - ab2*PHI(i,j,km1,nbl) - bb4*PHI(i,j,km2,nbl))
                  RHS(k) = bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl)) + ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl))
	  
               enddo

               ! Boundary points defined separately so that it dont fuck up.

               k=1
               kp2 = k+2
               kp1 = k+1  
               km1 = NK(nbl)-1
               km2 = NK(nbl)-2

               LLz = PHI(i,j,NK(nbl),nbl) - PHI(i,j,1,nbl)

               ! PHID(i,j,k,nbl) = (bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl) - LLz)+ ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl) - LLz))
               RHS(k) = bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl)+ LLz) + ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl)+ LLz)


               k=2
               kp2 = k+2  
               kp1 = k+1 
               km1 = k-1
               km2 = NK(nbl)-1

               ! check if the llx implementation is correct here

               ! PHID(i,j,k,nbl) = (bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl))+ ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl)))
               RHS(k) = bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl)+ LLz) + ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl))
	  


               k=NK(nbl)-1
               kp1 = k+1
               kp2 = 2
               km1 = k-1
               km2 = k-2

               ! PHID(i,j,k,nbl) = (bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl) + LLz)+ ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl)))
               RHS(k) = bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl)+ LLz) + ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl))


               k=NK(nbl)
               kp2 = 3
               kp1 = 2
               km1 = k-1
               km2 = k-2

               ! PHID(i,j,k,nbl) = (bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl) + LLz)+ ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl) + LLz))
               RHS(k) = bb4*(PHI(i,j,kp2,nbl) - PHI(i,j,km2,nbl)+ LLz) + ab2*(PHI(i,j,kp1,nbl) - PHI(i,j,km1,nbl)+ LLz)


              

               nm1 = NK(nbl)-1
               ! print*, RHS(1:nm1)
               call TDMAP(1,nm1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)

               PHID(i,j,1:nm1,nbl) = RHS(1:nm1)

               PHID(i,j,NK(nbl),nbl) = PHID(i,j,1,nbl) ! periodic BC; first and last point are the same

               ! print*, RHS(1:nm1)

               !ap - super diagonal
               !ac - diagonal
               !am - sub diagonal
               !ji - rhs

               ! the fucntion overwrites the RHS with the derivatives

               ! pause



            enddo
         enddo
      enddo

      END		  
	  	  
!*************************************************************************************
!**************************  FILTERING SUBROUTINES ***********************************
!*************************************************************************************

      ! SUBROUTINE FILTERING_I(PHI,nvars)
      ! use declare_variables	 
      ! implicit none	  
	  
      ! integer var,ipc,imc,coeff,nvars  
      ! real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
      ! real,dimension(NImax)	:: RHS

      ! END		  


      ! SUBROUTINE FILTERING_J(PHI,nvars)
      ! use declare_variables	 
      ! implicit none	  
	  
      ! integer var,jpc,jmc,coeff,nvars  
      ! real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
      ! real,dimension(NJmax)	:: RHS


      ! END	


      ! SUBROUTINE FILTERING_K(PHI,nvars)
      ! use declare_variables	 
      ! implicit none	  
	  
      ! integer var,kpc,kmc,coeff,nvars  
      ! real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
      ! real,dimension(NKmax)	:: RHS

	  
      ! END	
      
      
SUBROUTINE FILTERING_I(PHI,nvars) !Qc = PHI
use declare_variables	 
implicit none	

! print*, "line1799, fil i"

integer var,ipc,imc,coeff,nvars,nm1
real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
real,dimension(NImax) :: RHS

! print*, "line1805, fil i"

coeff = fscheme/2 + 1

! print*, "line1809, fil i"
	  
do var = 1,nvars
   do nbl = 1,nblocks
      do k=1,NK(nbl)
         do j=1,NJ(nbl)
            do i=1,NI(nbl)
               RHS(i) = 0.d0

               ! print*, "line1818, fil i"

               RHS(i) = fcoeff(1) * PHI(i,j,k,nbl,var)   ! n=0 term: a0 * phi_i


               do coeff =1,fscheme/2
                  ipc = i+coeff
                  imc = i-coeff
                  
                  if(ipc.gt.NI(nbl)) ipc = ipc - NI(nbl) + 1
                  ! print*, "line1825, fil i"
                  if(imc.lt.1) 	   imc = imc + NI(nbl) - 1
                     RHS(i) = RHS(i) + 0.5d0*fcoeff(coeff+1)*(PHI(ipc,j,k,nbl,var)+PHI(imc,j,k,nbl,var))

                  ! print*, "line1829, fil i"
               enddo

               ! print*, "line 1829, fil i"

            enddo
            nm1 = NI(nbl)-1			
   
            call TDMAP(1,nm1,APF(1:nm1),ACF(1:nm1),AMF(1:nm1),RHS(1:nm1),nm1)
            !do i=1,nm1
            !print*,'filtering:', PHI(i,j,k,nbl,var), RHS(i) !RHS = filtered
            !enddo
            !pause

            ! print*, "line 1840, fil i"

            PHI(1:nm1, j,k,nbl,var)= RHS(1:nm1)
            PHI(NI(nbl),j,k,nbl,var) = PHI(1,j,k,nbl,var) 

            ! print*, "line 1845, fil i"
                  
         enddo
      enddo

   enddo
enddo

END SUBROUTINE
						
					


					
SUBROUTINE FILTERING_J(PHI,nvars)
use declare_variables
implicit none

integer var,jpc,jmc,coeff,nvars,nm1
real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
real,dimension(NJmax) :: RHS

do var = 1,nvars
   do nbl = 1,nblocks
      do k = 1,NK(nbl)
         do i = 1,NI(nbl)
            do j = 1,NJ(nbl)

               RHS(j) = fcoeff(1) * PHI(i,j,k,nbl,var)   ! n=0 term

               do coeff = 1,fscheme/2
                  jpc = j + coeff
                  jmc = j - coeff

                  if(jpc.gt.NJ(nbl)) jpc = jpc - NJ(nbl) + 1
                  if(jmc.lt.1)       jmc = jmc + NJ(nbl) - 1

                  RHS(j) = RHS(j) + 0.5d0*fcoeff(coeff+1)*(PHI(i,jpc,k,nbl,var)+PHI(i,jmc,k,nbl,var))
               enddo

            enddo
            nm1 = NJ(nbl)-1

            call TDMAP(1,nm1,APF(1:nm1),ACF(1:nm1),AMF(1:nm1),RHS(1:nm1),nm1)

            PHI(i,1:nm1,k,nbl,var) = RHS(1:nm1)
            PHI(i,NJ(nbl),k,nbl,var) = PHI(i,1,k,nbl,var)

         enddo
      enddo
   enddo
enddo

END SUBROUTINE










SUBROUTINE FILTERING_K(PHI,nvars)
use declare_variables
implicit none

integer var,kpc,kmc,coeff,nvars,nm1
real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
real,dimension(NKmax) :: RHS

do var = 1,nvars
   do nbl = 1,nblocks
      do j = 1,NJ(nbl)
         do i = 1,NI(nbl)
            do k = 1,NK(nbl)

               RHS(k) = fcoeff(1) * PHI(i,j,k,nbl,var)   ! n=0 term

               do coeff = 1,fscheme/2
                  kpc = k + coeff
                  kmc = k - coeff

                  if(kpc.gt.NK(nbl)) kpc = kpc - NK(nbl) + 1
                  if(kmc.lt.1)       kmc = kmc + NK(nbl) - 1

                  RHS(k) = RHS(k) + 0.5d0*fcoeff(coeff+1)*(PHI(i,j,kpc,nbl,var)+PHI(i,j,kmc,nbl,var))
               enddo

            enddo
            nm1 = NK(nbl)-1

            call TDMAP(1,nm1,APF(1:nm1),ACF(1:nm1),AMF(1:nm1),RHS(1:nm1),nm1)

            PHI(i,j,1:nm1,nbl,var) = RHS(1:nm1)
            PHI(i,j,NK(nbl),nbl,var) = PHI(i,j,1,nbl,var)

         enddo
      enddo
   enddo
enddo

END SUBROUTINE
	  
	  

   


   !    SUBROUTINE FILTERING_J(PHI,nvars)
   !    use declare_variables	 
   !    implicit none	  
	  
   !    integer var,jpc,jmc,coeff,nvars,nm1
   !    real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
   !    real,dimension(NJmax)	:: RHS
	  
	!   do var = 1,nvars
	! 	do nbl = 1,nblocks
	! 		do k=1,NK(nbl)
	! 		do i=1,NI(nbl)
	! 			do j=1,NJ(nbl)
					
	! 				RHS(j) = 0.d0
   !             RHS(j) = fcoeff(1) * PHI(i,j,k,nbl,var)   ! n=0 term: a0 * phi_j
   !             ! is this correct? life is wakaranai.

	! 					do coeff =1,fscheme/2
	! 					jpc = j+coeff
	! 					jmc = j-coeff
						
	! 					if(jpc.gt.NJ(nbl)) jpc = jpc - NJ(nbl) + 1
	! 					if(jmc.lt.1) 	   jmc = jmc + NJ(nbl) - 1
	! 						RHS(j) = RHS(j) + 0.5d0*fcoeff(coeff+1)*(PHI(i,jpc,k,nbl,var)+PHI(i,jmc,k,nbl,var))
	! 						enddo
	! 						enddo
	! 						nm1 = NJ(nbl)-1			
				
	! 						call TDMAP(1,nm1,APF(1:nm1),ACF(1:nm1),AMF(1:nm1),RHS(1:nm1),nm1)
	! 						!do j=1,nm1
	! 						!print*,'filtering:', PHI(i,j,k,nbl,var), RHS(j) !RHS = filtered
	! 						!enddo
	! 						!pause
	! 						PHI(i, 1:nm1,k,nbl,var)= RHS(1:nm1)
	! 						PHI(i,NJ(nbl),k,nbl,var) = PHI(i,1,k,nbl,var) 
							
	! 						enddo
	! 						enddo
	! 						enddo
	! 						enddo


   !    END SUBROUTINE


   !    SUBROUTINE FILTERING_K(PHI,nvars)
   !    use declare_variables	 
   !    implicit none	  
	  
   !    integer var,kpc,kmc,coeff,nvars,nm1
   !    real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
   !    real,dimension(NKmax)	:: RHS
	  
	!   do var = 1,nvars
	! 	do nbl = 1,nblocks
			
			
	! 			do j=1,NJ(nbl)
	! 			do i=1,NI(nbl)
	! 			do k=1,NK(nbl)
					
	! 				RHS(k) = 0.d0
	! 					do coeff =1,fscheme/2+1
	! 					kpc = k+coeff
	! 					kmc = k-coeff
						
	! 					if(kpc.gt.NK(nbl)) kpc = kpc - NK(nbl) + 1
	! 					if(kmc.lt.1) 	   kmc = kmc + NK(nbl) - 1
	! 						RHS(k) = RHS(k) + 0.5d0*fcoeff(coeff+1)*(PHI(i,j,kpc,nbl,var)+PHI(i,j,kmc,nbl,var))
	! 						enddo
	! 						enddo
	! 						nm1 = NK(nbl)-1			
				
	! 						call TDMAP(1,nm1,APF(1:nm1),ACF(1:nm1),AMF(1:nm1),RHS(1:nm1),nm1)
	! 						!do k=1,nm1
	! 						!print*,'filtering:', PHI(i,j,k,nbl,var), RHS(k) !RHS = filtered
	! 						!enddo
	! 						!pause
	! 						PHI(i, j,1:nm1,nbl,var)= RHS(1:nm1)
	! 						PHI(i,j,NK(nbl),nbl,var) = PHI(i,j,1,nbl,var) 
							
	! 						enddo
	! 						enddo
	! 						enddo
	! 						enddo

	  
   !    End SUBROUTINE














!************************** CYCLIC TDMA SOLVER *************************************	  

subroutine TDMAP(ji,jf,ap,ac,am,fi,NMAXL)

!ap - super diagonal
!ac - diagonal
!am - sub diagonal

implicit none

! -----------------------
!  Input/Output variables
! -----------------------
   integer:: ji, jf, NMAXL
   real:: ap(NMAXL), ac(NMAXL), am(NMAXL), fi(NMAXL)

! -------------------
!  Internal variables
! -------------------
   integer:: i, j, ja, jj
   real:: fnn, pp
   real:: qq(NMAXL), ss(NMAXL), fei(NMAXL)

   ja=ji+1
   jj=ji+jf

   qq(ji)=-ap(ji)/ac(ji)
   ss(ji)=-am(ji)/ac(ji)
   fnn=fi(jf)
   fi(ji)=fi(ji)/ac(ji)

!  forward elimination sweep
!----------------------------
   do j=ja,jf
      pp=1.0d0/(ac(j)+am(j)*qq(j-1))
          qq(j)=-ap(j)*pp
          ss(j)=-am(j)*ss(j-1)*pp
          fi(j)=(fi(j)-am(j)*fi(j-1))*pp
   enddo
   
!  backward pass
!----------------

   ss(jf)=1.0d0
   fei(jf)=0.0d0
 
   do i=ja,jf
      j=jj-i
          ss(j)=ss(j)+qq(j)*ss(j+1)
          fei(j)=fi(j)+qq(j)*fei(j+1)
   enddo
   
   fi(jf)=(fnn-ap(jf)*fei(ji)-am(jf)*fei(jf-1))/    &
   &      (ap(jf)*ss(ji)+am(jf)*ss(jf-1)+ac(jf)) 
   
!  backward substitution
!------------------------
   do i=ja,jf
      j=jj-i
          fi(j)=fi(jf)*ss(j)+fei(j)
   enddo

end subroutine TDMAP