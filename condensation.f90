module condensation

    implicit none

    private
    public get_large_scale_condensation_tendencies

    ! Constants for large-scale condensation
    real, parameter :: trlsc  = 4.0  !! Relaxation time (in hours) for specific humidity
    real, parameter :: rhlsc  = 0.9  !! Maximum relative humidity threshold (at sigma=1)
    real, parameter :: drhlsc = 0.1  !! Vertical range of relative humidity threshold
    real, parameter :: rhblsc = 0.95 !! Relative humidity threshold for boundary layer

contains
    subroutine get_large_scale_condensation_tendencies(psa, qa, qsat, itop, &
                                       & ngp, nlev, sig, &
                                       & precls, dtlsc, dqlsc)
        use physical_constants, only: p0, cp, alhc, grav
        use utils, only: get_sigma_thicknesses

        real, intent(in) :: psa(ngp)
        real, intent(in) :: qa(ngp,nlev)
        real, intent(in) :: qsat(ngp,nlev)
        integer, intent(inout) :: itop(ngp)
        real, intent(out) :: precls(ngp)
        integer, intent(in) :: ngp
        integer, intent(in) :: nlev
        real, intent(in) :: sig(nlev)
        real, intent(out) :: dtlsc(ngp,nlev)
        real, intent(out) :: dqlsc(ngp,nlev)

        real :: psa2(ngp), dqa, tfact, sig2, rtlsc, rhref, qsmax, prg, pfact, dqmax
        real :: dsig(nlev)
        integer :: j, k

        ! 1. Initialization
        dsig = get_sigma_thicknesses(sig)
        qsmax = 10.0

        rtlsc = 1.0/(trlsc*3600.0)
        tfact = alhc/cp
        prg = p0/grav

        do j = 1, ngp
            dtlsc(j,1) = 0.0
            dqlsc(j,1) = 0.0
            precls(j)  = 0.0
            psa2(j)    = psa(j)*psa(j)
        end do

        ! 2. Tendencies of temperature and moisture
        !    NB. A maximum heating rate is imposed to avoid 
        !        grid-point-storm instability 
        do k = 2, nlev
            sig2 = sig(k)*sig(k)
            rhref = rhlsc + drhlsc*(sig2 - 1.0)
            if (k == nlev) rhref = max(rhref, rhblsc)
            dqmax = qsmax*sig2*rtlsc
  
            do j = 1, ngp
                dqa = rhref*qsat(j,k) - qa(j,k)
                if (dqa < 0.0) then
                    itop(j)    = min(k, itop(j))
                    dqlsc(j,k) = dqa*rtlsc
                    dtlsc(j,k) = tfact*min(-dqlsc(j,k), dqmax*psa2(j))
                else
                    dqlsc(j,k) = 0.0
                    dtlsc(j,k) = 0.0
                endif
            end do
        end do

        ! 3. Large-scale precipitation
        do k = 2, nlev
            pfact = dsig(k)*prg
            do j = 1, ngp
                precls(j) = precls(j) - pfact*dqlsc(j,k)
            end do
        end do

        do j = 1, ngp
            precls(j) = precls(j)*psa(j)
        end do
    end subroutine get_large_scale_condensation_tendencies
end module condensation
