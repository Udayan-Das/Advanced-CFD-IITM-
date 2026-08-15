

if (dscheme <= 4) then
    ! ==== EXISTING CENTRED SCHEMES (E2, E4, C4, C6) ==== 
    call DISCRETIZATION_I_EXP(Fflux, fluxD, nconserv) 
    net_flux = net_flux + fluxD
    ... (same as before)
else
! ==== WENO SCHEMES (dscheme = 5,6,7,8) ==== 
    call WENO_FLUX_I(Qc, Fflux_raw, net_flux, dscheme) 
    call WENO_FLUX_J(Qc, Gflux_raw, net_flux, dscheme)
    if (grid2d .ne. 1) call WENO_FLUX_K(Qc, Hflux_raw, net_flux, dscheme) 

endif


subroutine WENO_FLUX_I(Qc, Fflux_raw, net_flux, wscheme)



    !=================================================================
    ! Step 1: Find the maximum wave speed alpha across the whole grid
    !=================================================================
    alpha = 0.d0
    do nbl = 1, nblocks
      do k = 1, NK(nbl)
        do j = 1, NJ(nbl)
          do i = 1, NI(nbl)

            c = sqrt(gamma * Qp(i,j,k,nbl,5) / Qp(i,j,k,nbl,1))

            Ucont = math the contravariant velocity again. in i direction

            grad_xi = sqrt(ix(i,j,k,nbl)**2 + iy(i,j,k,nbl)**2 + iz(i,j,k,nbl)**2)

            alpha = max(alpha, abs(Ucont) + c*grad_xi)

          enddo
        enddo
      enddo
    enddo


    ! ==================================================================
    ! step 2: loopeth over each variable (5 components) and each (j,k) line.
    ! for each line, do LF splitting + WENO reconstruction + derivative.
    ! ==================================================================
    do var = 1, 5
        do nbl = 1, nblocks
            do k = 1, NK(nbl)
                do j = 1, NJ(nbl)

                    ! -------- step 2a: LF splitting along this line --------
                    do i = 1, NI(nbl)
                    vol = 1.d0 / Jac(i,j,k,nbl)
                    Fp(i) = 0.5d0 * (Fflux_raw(i,j,k,nbl,var) + alpha * Qc(i,j,k,nbl,var) * vol)
                    Fm(i) = 0.5d0 * (Fflux_raw(i,j,k,nbl,var) - alpha * Qc(i,j,k,nbl,var) * vol)
                    enddo

                    ! -------- step 2b: reconstruct F_hat at each interface i+1/2 --------
                    do i = 1, NI(nbl)

                        ! periodic wrapping (point NI is the same as point 1)
                        im2 = i-2;  if (im2 < 1)  im2 = im2 + NI(nbl) - 1
                        im1 = i-1;  if (im1 < 1)  im1 = im1 + NI(nbl) - 1
                        ic  = i
                        ip1 = i+1;  if (ip1 > NI(nbl)) ip1 = ip1 - NI(nbl) + 1
                        ip2 = i+2;  if (ip2 > NI(nbl)) ip2 = ip2 - NI(nbl) + 1
                        ip3 = i+3;  if (ip3 > NI(nbl)) ip3 = ip3 - NI(nbl) + 1

                        ! ===== GETTING FHAT FROM F+ =====

                        ! Candidate fluxes (3-point stencils)
                        fhat0 = (2*Fp(im2) - 7*Fp(im1) + 11*Fp(ic)) / 6.d0
                        fhat1 = (-Fp(im1) + 5*Fp(ic) + 2*Fp(ip1)) / 6.d0
                        fhat2 = (2*Fp(ic) + 5*Fp(ip1) - Fp(ip2)) / 6.d0

                        ! Extra downwind candidate (only for CU6 / CU6-M)
                        if (wscheme == 3 .or. wscheme == 4) then
                            fhat3 = (11*Fp(ip1) - 7*Fp(ip2) + 2*Fp(ip3)) / 6.d0
                        endif

                        ! Smoothness indicators for the 3-point stencils
                        b0 = (13.d0/12)*(Fp(im2) - 2*Fp(im1) + Fp(ic))**2 &
                            + 0.25d0*(Fp(im2) - 4*Fp(im1) + 3*Fp(ic))**2
                        b1 = (13.d0/12)*(Fp(im1) - 2*Fp(ic) + Fp(ip1))**2 &
                            + 0.25d0*(Fp(im1) - Fp(ip1))**2
                        b2 = (13.d0/12)*(Fp(ic) - 2*Fp(ip1) + Fp(ip2))**2 &
                            + 0.25d0*(3*Fp(ic) - 4*Fp(ip1) + Fp(ip2))**2

                        ! Full-stencil smoothness indicator (only for CU6 / CU6-M)
                        if (wscheme == 3 .or. wscheme == 4) then
                            b6 = (1.d0/120960.d0) * ( &
                                    271779.d0*Fp(im2)**2 &
                                + Fp(im2)*(-2380800*Fp(im1) + 4086352*Fp(ic) &
                                            - 3462252*Fp(ip1) + 1458762*Fp(ip2) &
                                            - 245620*Fp(ip3)) &
                                + Fp(im1)*(5653317*Fp(im1) - 20427884*Fp(ic) &
                                            + 17905032*Fp(ip1) - 7727988*Fp(ip2) &
                                            + 1325006*Fp(ip3)) &
                                + Fp(ic)*(19510972*Fp(ic) - 35817664*Fp(ip1) &
                                            + 15929912*Fp(ip2) - 2792660*Fp(ip3)) &
                                + Fp(ip1)*(17195652*Fp(ip1) - 15880404*Fp(ip2) &
                                            + 2863984*Fp(ip3)) &
                                + Fp(ip2)*(3824847*Fp(ip2) - 1429976*Fp(ip3)) &
                                + 139633.d0*Fp(ip3)**2 )
                        endif

                        ! ===== Compute nonlinear weights based on scheme =====

                        if (wscheme == 1) then              ! -------- WENO-JS --------
                            eps = 1.d-6
                            a0 = 0.1d0 / (eps + b0)**2
                            a1 = 0.6d0 / (eps + b1)**2
                            a2 = 0.3d0 / (eps + b2)**2
                            asum = a0 + a1 + a2
                            w0 = a0/asum;  w1 = a1/asum;  w2 = a2/asum
                            Fhat_plus = w0*fhat0 + w1*fhat1 + w2*fhat2

                        elseif (wscheme == 2) then          ! -------- WENO-Z --------
                            eps_z = 1.d-40
                            tau5 = abs(b0 - b2)
                            a0 = 0.1d0 * (1.d0 + tau5/(eps_z + b0))
                            a1 = 0.6d0 * (1.d0 + tau5/(eps_z + b1))
                            a2 = 0.3d0 * (1.d0 + tau5/(eps_z + b2))
                            asum = a0 + a1 + a2
                            w0 = a0/asum;  w1 = a1/asum;  w2 = a2/asum
                            Fhat_plus = w0*fhat0 + w1*fhat1 + w2*fhat2

                        elseif (wscheme == 3) then          ! -------- WENO-CU6 --------
                            eps = 1.d-8
                            C_cu = 20.d0
                            tau6 = b6 - (1.d0/6.d0)*(b0 + b2 + 4.d0*b1)
                            a0 = (1.d0/20.d0) * (C_cu + tau6/(b0 + eps))
                            a1 = (9.d0/20.d0) * (C_cu + tau6/(b1 + eps))
                            a2 = (9.d0/20.d0) * (C_cu + tau6/(b2 + eps))
                            a3 = (1.d0/20.d0) * (C_cu + tau6/(b6 + eps))
                            asum = a0 + a1 + a2 + a3
                            w0 = a0/asum;  w1 = a1/asum;  w2 = a2/asum;  w3 = a3/asum
                            Fhat_plus = w0*fhat0 + w1*fhat1 + w2*fhat2 + w3*fhat3

                        elseif (wscheme == 4) then          ! -------- WENO-CU6-M --------
                            eps = 1.d-8
                            C_cu = 1000.d0
                            tau6 = b6 - (1.d0/6.d0)*(b0 + b2 + 4.d0*b1)
                            a0 = (1.d0/20.d0) * (C_cu + tau6/(b0 + eps))**4
                            a1 = (9.d0/20.d0) * (C_cu + tau6/(b1 + eps))**4
                            a2 = (9.d0/20.d0) * (C_cu + tau6/(b2 + eps))**4
                            a3 = (1.d0/20.d0) * (C_cu + tau6/(b6 + eps))**4
                            asum = a0 + a1 + a2 + a3
                            w0 = a0/asum;  w1 = a1/asum;  w2 = a2/asum;  w3 = a3/asum
                            Fhat_plus = w0*fhat0 + w1*fhat1 + w2*fhat2 + w3*fhat3
                        endif

                        ! ===== RIGHT-BIASED RECONSTRUCTION OF F- =====
                        ! (mirror image of the above; see Section 5 for details)
                        ! Fhat_minus = ... (analogous code using Fm with reversed indexing)

                        ! Total interface flux
                        Fhat_total(i) = Fhat_plus + Fhat_minus

                    enddo    ! end of interface loop

                    ! -------- step 2c: math the flux derivative along this line --------
                    do i = 1, NI(nbl)
                    im1_iface = i-1;  if (im1_iface < 1) im1_iface = im1_iface + NI(nbl) - 1
                    net_flux(i,j,k,nbl,var) = net_flux(i,j,k,nbl,var) &
                        - (Fhat_total(i) - Fhat_total(im1_iface))
                    enddo

                enddo    ! end of le j loop
            enddo      ! end of le k loop
        enddo        ! end of le nbl loop
    enddo          ! end of le var loop

end subroutine WENO_FLUX_I




SUBROUTINE UNSTEADY(stepl)
  ...

  ! step 1: math fluxes at each node (same as before, 
  ! but WITHOUT the negative sign when using WENO)
  if (dscheme <= 4) then
    ! centred schemes: bake in the minus sign (existing behaviour)
    Fflux(i,j,k,nbl,1) = -rhl*Ucont * vol
    ...
  else
    ! WENO: no minus sign (we will apply it in WENO_FLUX_I)
    Fflux_raw(i,j,k,nbl,1) = rhl*Ucont * vol
    ...
  endif

  ! step 2: diff
  net_flux = 0.d0
  if (dscheme <= 4) then
    ! (existing)
    call DISCRETIZATION_I_EXP(Fflux, fluxD, nconserv)
    net_flux = net_flux + fluxD
    ... (J, K)
  else
    ! WENO path (new)
    call WENO_FLUX_I(Qc, Fflux_raw, net_flux, dscheme)
    call WENO_FLUX_J(Qc, Gflux_raw, net_flux, dscheme)
    if (grid2d .ne. 1) call WENO_FLUX_K(Qc, Hflux_raw, net_flux, dscheme)
  endif

  ! step 3: viscous terms (unchanged)
  if (viscous .eq. 1) then
    ... (existing viscous computation)
  endif

  ! step 4: RK4 time update (unchanged!)
  do var = 1, nconserv
    do nbl, k, j, i
      Qcnew(...) += time_step * net_flux(...) * fac_RK(stepl) / vol
      if (stepl <= rk_steps-1) then
        Qc(...) = Qcini(...) + time_step * net_flux(...) * fac_qini(stepl+1) / vol
      else
        Qc(...) = Qcnew(...)
      endif
    enddo
  enddo
END SUBROUTINE
