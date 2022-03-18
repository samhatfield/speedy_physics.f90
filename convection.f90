module convection

    implicit none

    private
    public get_convection_tendencies

    ! Constants for convection
    real, parameter :: psmin  = 0.8 !! Minimum (normalised) surface pressure for the occurrence of
                                    !! convection
    real, parameter :: trcnv  = 6.0 !! Time of relaxation (in hours) towards reference state
    real, parameter :: rhbl   = 0.9 !! Relative humidity threshold in the boundary layer
    real, parameter :: rhil   = 0.7 !! Relative humidity threshold in intermeduate layers for
                                    !! secondary mass flux
    real, parameter :: entmax = 0.5 !! Maximum entrainment as a fraction of cloud-base mass flux
    real, parameter :: smf    = 0.8 !! Ratio between secondary and primary mass flux at cloud-base
 
contains
    subroutine get_convection_tendencies(prog_sp, stat_en, prog_q, qsat, &
                                       & ngp, nlev, sig, &
                                       & itop, cbmf, precnv, dfse, dfqa)
        use physical_constants, only: p0, alhc, grav
        use utils, only: get_sigma_thicknesses, get_vertical_interpolation_weights

        real, intent(in) :: prog_sp(:)
        real, intent(in) :: stat_en(:,:)
        real, intent(in) :: prog_q(:,:)
        real, intent(in) :: qsat(:,:)
        integer, intent(in) :: ngp
        integer, intent(in) :: nlev
        real, intent(in) :: sig(nlev)
        integer, intent(out) :: itop(:)
        real, intent(out) :: cbmf(:)
        real, intent(out) :: precnv(:)
        real, intent(out) :: dfse(:,:)
        real, intent(out) :: dfqa(:,:)

        integer :: nlp, ktop1, ktop2, nl1, j, k, k1

        real :: mss(ngp,2:nlev), mse0, mse1, mss0, mss2, msthr
        real :: qdif(ngp)
        real :: entr(2:nlev-1), delq, enmass, fdq, fds, fm0, fmass, fpsa, fqmax
        real :: fsq, fuq, fus, qb, qmax, qsatb, qthr0, qthr1, rdps, rlhc, sb, sentr
        real :: wvi(ngp,2)
        real :: dsig(nlev)
        
        logical :: lqthr

        dsig = get_sigma_thicknesses(sig)
        wvi = get_vertical_interpolation_weights(sig)

        ! 1. Initialization of output and workspace arrays
        nl1 = nlev - 1
        nlp = nlev + 1
        fqmax = 5.0

        fm0 = p0*dsig(nlev)/(grav*trcnv*3600.0)
        rdps=2.0/(1.0 - psmin)
 
        dfse = 0.0
        dfqa = 0.0
        cbmf = 0.0
        precnv = 0.0
  
        ! Entrainment profile (up to sigma = 0.5)
        sentr = 0.0
        do k = 2, nl1
            entr(k) = (max(0.0, sig(k) - 0.5))**2.0
            sentr = sentr + entr(k)
        end do
   
        sentr = entmax/sentr
        entr(2:nl1) = entr(2:nl1)*sentr

        ! 2. Check of conditions for convection
        rlhc = 1.0/alhc

        do j = 1, ngp
            itop(j) = nlp
  
            if (prog_sp(j) > psmin) then
                ! Minimum of moist static energy in the lowest two levels
                mse0 = stat_en(j,nlev) + alhc*prog_q(j,nlev)
                mse1 = stat_en(j,nl1)  + alhc*prog_q(j,nl1)
                mse1 = min(mse0, mse1)
  
                ! Saturation (or super-saturated) moist static energy in PBL 
                mss0 = max(mse0, mss(j,nlev))
  
                ktop1 = nlev
                ktop2 = nlev
  
                do k = nlev - 3, 3, -1
                    mss2 = mss(j,k) + wvi(k,2)*(mss(j,k+1) - mss(j,k))
  
                    ! Check 1: conditional instability 
                    !          (MSS in PBL > MSS at top level)
                    if (mss0 > mss2) then
                        ktop1 = k
                    end if
  
                    ! Check 2: gradient of actual moist static energy 
                    !          between lower and upper troposphere                     
                    if (mse1 > mss2) then
                       ktop2 = k
                       msthr = mss2
                    end if
                end do
  
                if (ktop1 < nlev) then
                    ! Check 3: RH > RH_c at both k=NLEV and k=NL1
                    qthr0 = rhbl*qsat(j,nlev)
                    qthr1 = rhbl*qsat(j,nl1)
                    lqthr = (prog_q(j,nlev) > qthr0 .and. prog_q(j,nl1) > qthr1)
  
                    if (ktop2 < nlev) then
                        itop(j) = ktop1
                        qdif(j) = max(prog_q(j,nlev) - qthr0, (mse0 - msthr)*rlhc)
                    else if (lqthr) then
                        itop(j) = ktop1
                        qdif(j) = prog_q(j,nlev) - qthr0
                    endif
                end if
            end if
        end do

        ! 3. Convection over selected grid-points
        do j = 1, ngp
            if (itop(j) == nlp) cycle

            ! 3.1 Boundary layer (cloud base)
            k = nlev
            k1 = k-1

            ! Maximum specific humidity in the PBL
            qmax = max(1.01*prog_q(j,k), qsat(j,k))

            ! Dry static energy and moisture at upper boundary
            sb = stat_en(j,k1) + wvi(k1,2)*(stat_en(j,k) - stat_en(j,k1))
            qb = prog_q(j,k1) + wvi(k1,2)*(prog_q(j,k) - prog_q(j,k1))
            qb = min(qb, prog_q(j,k))

            ! Cloud-base mass flux, computed to satisfy:
            ! fmass*(qmax-qb)*(g/dp)=qdif/trcnv
            fpsa = prog_sp(j)*min(1.0, (prog_sp(j) - psmin)*rdps)
            fmass = fm0*fpsa*min(fqmax, qdif(j)/(qmax - qb))
            cbmf(j) = fmass

            ! Upward fluxes at upper boundary
            fus = fmass*stat_en(j,k)
            fuq = fmass*qmax

            ! Downward fluxes at upper boundary
            fds = fmass*sb
            fdq = fmass*qb

            ! Net flux of dry static energy and moisture
            dfse(j,k) = fds - fus
            dfqa(j,k) = fdq - fuq

            ! 3.2 Intermediate layers (entrainment)
            do k = nlev - 1, itop(j) + 1, -1
                k1 = k - 1

                ! Fluxes at lower boundary
                dfse(j,k) = fus - fds
                dfqa(j,k) = fuq - fdq

                ! Mass entrainment
                enmass = entr(k)*prog_sp(j)*cbmf(j)
                fmass = fmass + enmass

                ! Upward fluxes at upper boundary
                fus = fus + enmass*stat_en(j,k)
                fuq = fuq + enmass*prog_q(j,k)

                ! Downward fluxes at upper boundary
                sb = stat_en(j,k1) + wvi(k1,2)*(stat_en(j,k) - stat_en(j,k1))
                qb = prog_q(j,k1) + wvi(k1,2)*(prog_q(j,k) - prog_q(j,k1))
                fds = fmass*sb
                fdq = fmass*qb

                ! Net flux of dry static energy and moisture
                dfse(j,k) = dfse(j,k) + fds - fus
                dfqa(j,k) = dfqa(j,k) + fdq - fuq

                ! Secondary moisture flux
                delq = rhil*qsat(j,k) - prog_q(j,k)
                if (delq > 0.0) then
                    fsq = smf*cbmf(j)*delq
                    dfqa(j,k)    = dfqa(j,k)    + fsq 
                    dfqa(j,nlev) = dfqa(j,nlev) - fsq
                end if
            end do

            ! 3.3 Top layer (condensation and detrainment)
            k = itop(j)

            ! Flux of convective precipitation
            qsatb = qsat(j,k) + wvi(k,2)*(qsat(j,k+1) - qsat(j,k))

            precnv(j) = max(fuq - fmass*qsatb, 0.0)

            ! Net flux of dry static energy and moisture
            dfse(j,k) = fus - fds+alhc*precnv(j)
            dfqa(j,k) = fuq - fdq - precnv(j)
        end do
    end subroutine get_convection_tendencies
end module convection

