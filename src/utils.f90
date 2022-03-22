module utils
    implicit none

    private
    public get_sigma_thicknesses
    public get_vertical_interpolation_weights
    public get_flux_conv_uvq

contains
    function get_half_sigma(sig) result(sigh)
        real, intent(in) :: sig(:)

        real :: sigh(size(sig)+1)

        integer :: k

        sigh(1) = 0.0
        do k = 2, size(sig)+1
            sigh(k) = 2*sig(k-1) - sigh(k-1)
        end do
    end function get_half_sigma

    function get_sigma_thicknesses(sig) result(dsig)
        real, intent(in) :: sig(:)

        real :: dsig(size(sig))

        real :: hsg(size(sig)+1)

        integer :: k

        hsg = get_half_sigma(sig)

        do k = 1, size(sig)
            dsig(k) = hsg(k+1) - hsg(k)
        end do
    end function get_sigma_thicknesses

    function get_flux_conv_uvq(sig) result(grdsig)
        use physical_constants, only: grav, p0

        real, intent(in) :: sig(:)

        real :: dsig(size(sig))

        real :: grdsig(size(sig))

        dsig = get_sigma_thicknesses(sig)

        grdsig = grav/(dsig*p0)
    end function get_flux_conv_uvq

    function get_vertical_interpolation_weights(sig) result(wvi)
        use physical_constants, only: grav, cp, p0

        real, intent(in) :: sig(:)

        real :: wvi(size(sig),2)

        integer :: k, nlev
        real :: sigh(size(sig)+1), sigl(size(sig))

        nlev = size(sig)

        ! 1.2 Functions of sigma and latitude
        sigh = get_half_sigma(sig)
        sigl = log(sig)

!        do k = 1, nlev
!            sigl(k) = log(fsg(k))
!            grdsig(k) = grav/(dhs(k)*p0)
!            grdscp(k) = grdsig(k)/cp
!        end do

        ! Weights for vertical interpolation at half-levels(1,nlev) and surface
        ! Note that for phys.par. half-lev(k) is between full-lev k and k+1
        ! Fhalf(k) = Ffull(k)+WVI(K,2)*(Ffull(k+1)-Ffull(k))
        ! Fsurf = Ffull(nlev)+WVI(nlev,2)*(Ffull(nlev)-Ffull(nlev-1))
        do k = 1, nlev-1
            wvi(k,1) = 1.0/(sigl(k+1) - sigl(k))
            wvi(k,2) = (log(sigh(k+1)) - sigl(k))*wvi(k,1)
        end do

        wvi(nlev,1) = 0.0
        wvi(nlev,2) = (log(0.99) - sigl(nlev))*wvi(nlev-1,1)
    end function get_vertical_interpolation_weights
end module utils

