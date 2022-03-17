module humidity
    use params, only: ngp

    implicit none

    private
    public spec_hum_to_rel_hum
    public rel_hum_to_spec_hum
    public get_qsat

contains
    subroutine spec_hum_to_rel_hum(prog_t, prog_sp, sig, prog_q, rh, qsat)
        real, intent(in) :: prog_t(ngp)
        real, intent(in) :: prog_sp(ngp)
        real, intent(in) :: sig
        real, intent(in) :: prog_q(ngp)
        real, intent(inout) :: rh(ngp)
        real, intent(inout) :: qsat(ngp)

        qsat = get_qsat(prog_t, prog_sp, sig)
        rh = prog_q/qsat
    end subroutine spec_hum_to_rel_hum

    subroutine rel_hum_to_spec_hum(prog_t, prog_sp, sig, rh, prog_q, qsat)
        real, intent(in) :: prog_t(ngp)
        real, intent(in) :: prog_sp(ngp)
        real, intent(in) :: sig
        real, intent(in) :: rh(ngp)
        real, intent(inout) :: prog_q(ngp)
        real, intent(inout) :: qsat(ngp)

        qsat = get_qsat(prog_t, prog_sp, sig)
        prog_q = rh*qsat
    end subroutine rel_hum_to_spec_hum

    function get_qsat(prog_t, prog_sp, sig) result(qsat)
        real, intent(in) :: prog_t(ngp) !! Absolute temperature
        real, intent(in) :: prog_sp(ngp) !! Normalized pressure (p/1000 hPa)
        real, intent(in) :: sig     !! Sigma level
        real :: qsat(ngp)           !! Saturation specific humidity in g/kg

        real :: e0, c1, c2, t0, t1, t2

        integer :: i

        ! 1. Compute Qsat (g/kg) from T (degK) and normalized pres. P (= p/1000_hPa)
        ! If sig > 0, P = Ps * sigma, otherwise P = Ps(1) = const.
        e0 = 6.108e-3
        c1 = 17.269
        c2 = 21.875
        t0 = 273.16
        t1 = 35.86
        t2 = 7.66

        do i = 1, ngp
            if (prog_t(i) >= t0) then
                qsat(i) = e0*exp(c1*(prog_t(i) - t0)/(prog_t(i) - t1))
            else
                qsat(i) = e0*exp(c2*(prog_t(i) - t0)/(prog_t(i) - t2))
            end if
        end do

        if (sig <= 0.0) then
            qsat = 622.0*qsat/(prog_sp(1) - 0.378*qsat)
        else
            qsat = 622.0*qsat/(sig*prog_sp - 0.378*qsat)
        end if
    end function get_qsat
end module humidity
