module radiation

    implicit none

    private
    public sol_oz, cloud, radsw

    ! Shortwave radiation and cloud constants
    real, parameter :: solc    = 342.0 !! Solar constant (area averaged) in W/m^2
    real, parameter :: rhcl1   = 0.30  !! Relative humidity threshold corresponding to
                                       !! cloud cover = 0
    real, parameter :: rhcl2   = 1.00  !! Relative humidity correponding to cloud cover = 1
    real, parameter :: qcl     = 0.20  !! Specific humidity threshold for cloud cover
    real, parameter :: wpcl    = 0.2   !! Cloud cover weight for the square-root of precipitation
                                       !! (for p = 1 mm/day)
    real, parameter :: pmaxcl  = 10.0  !! Maximum value of precipitation (mm/day) contributing to
                                       !! cloud cover
    real, parameter :: clsmax  = 0.60  !! Maximum stratiform cloud cover
    real, parameter :: clsminl = 0.15  !! Minimum stratiform cloud cover over land (for RH = 1)
    real, parameter :: gse_s0  = 0.25  !! Gradient of dry static energy corresponding to
                                       !! stratiform cloud cover = 0
    real, parameter :: gse_s1  = 0.40  !! Gradient of dry static energy corresponding to
                                       !! stratiform cloud cover = 1
    real, parameter :: albcl   = 0.43  !! Cloud albedo (for cloud cover = 1)
    real, parameter :: albcls  = 0.50  !! Stratiform cloud albedo (for st. cloud cover = 1)
    real, parameter :: epssw   = 0.020 !! Fraction of incoming solar radiation absorbed by ozone

    ! Shortwave absorptivities (for dp = 10^5 Pa)
    real, parameter :: absdry  = 0.033 !! Absorptivity of dry air (visible band)
    real, parameter :: absaer  = 0.033 !! Absorptivity of aerosols (visible band)
    real, parameter :: abswv1  = 0.022 !! Absorptivity of water vapour
                                       !! (visible band, for dq = 1 g/kg)
    real, parameter :: abswv2  = 15.0  !! Absorptivity of water vapour
                                       !! (near IR band, for dq = 1 g/kg)
    real, parameter :: abscl1  = 0.015 !! Absorptivity of clouds (visible band, maximum value)
    real, parameter :: abscl2  = 0.15  !! Absorptivity of clouds
                                       !! (visible band, for dq_base = 1 g/kg)

    ! Longwave absorptivities (per dp = 10^5 Pa)
    real, parameter :: ablwin =  0.3  !! Absorptivity of air in "window" band
    real            :: ablco2 =  6.0  !! Absorptivity of air in CO2 band
    real, parameter :: ablwv1 =  0.7  !! Absorptivity of water vapour in H2O band 1 (weak),
                                      !! (for dq = 1 g/kg)
    real, parameter :: ablwv2 = 50.0  !! Absorptivity of water vapour in H2O band 2 (strong),
                                      !! (for dq = 1 g/kg)
    real, parameter :: ablcl1 = 12.0  !! Absorptivity of "thick" clouds in window band
                                      !! (below cloud top)
    real, parameter :: ablcl2 =  0.6  !! Absorptivity of "thin" upper clouds in window and H2O
                                      !! bands

    real, parameter :: epslw = 0.05 !! Fraction of blackbody spectrum absorbed/emitted by PBL only

    real, parameter :: pi = 3.14159265359

contains

    subroutine sol_oz(tyear, nlon, lats, topsr, fsol, ozupp, ozone, zenit, stratz)
        real, intent(in) :: tyear
        integer, intent(in) :: nlon
        real, intent(in) :: lats(:)
        real, intent(out) :: topsr(size(lats))
        real, intent(out) :: fsol(nlon*size(lats))
        real, intent(out) :: ozupp(nlon*size(lats))
        real, intent(out) :: ozone(nlon*size(lats))
        real, intent(out) :: zenit(nlon*size(lats))
        real, intent(out) :: stratz(nlon*size(lats))

        integer :: ngp, nlat, nzen, i, j, j0

        real :: alpha, dalpha, coz1, coz2, azen, rzen, czen, szen, fs0, flat2
        real :: clat(size(lats)), slat(size(lats))

        nlat = size(lats)
        ngp = nlon*nlat

        ! Compute cosine and sine of latitude
        clat = cos(lats*pi/180.0)
        slat = sin(lats*pi/180.0)

        ! alpha = year phase ( 0 - 2pi, 0 = winter solstice = 22dec.h00 )
        alpha = 4.0*asin(1.0)*(tyear+10.0/365.0)
        dalpha = 0.0

        coz1 = 1.0*max(0.0, cos(alpha - dalpha))
        coz2 = 1.8
        azen = 1.0
        nzen = 2

        rzen = -cos(alpha)*23.45*asin(1.0)/90.0
        czen = cos(rzen)
        szen = sin(rzen)

        fs0 = 6.0

        ! Solar radiation at the top
        call solar(tyear, 4.0*solc, lats, topsr)

        do j = 1, nlat
            j0 = 1 + nlon*(j - 1)
            flat2 = 1.5*slat(j)**2 - 0.5

            ! Solar radiation at the top
            fsol(j0) = topsr(j)

            ! Ozone depth in upper and lower stratosphere
            ozupp(j0) = 0.5*epssw
            ozone(j0) = 0.4*epssw*(1.0 + coz1*slat(j) + coz2*flat2)

            ! Zenith angle correction to (downward) absorptivity
            zenit(j0) = 1.0 + azen*(1.0 - (clat(j)*czen + slat(j)*szen))**nzen

            ! Ozone absorption in upper and lower stratosphere
            ozupp(j0) = fsol(j0)*ozupp(j0)*zenit(j0)
            ozone(j0) = fsol(j0)*ozone(j0)*zenit(j0)

            ! Polar night cooling in the stratosphere
            stratz(j0) = max(fs0 - fsol(j0), 0.0)

            do i = 1, nlon - 1
                fsol  (i+j0) = fsol  (j0)
                ozone (i+j0) = ozone (j0)
                ozupp (i+j0) = ozupp (j0)
                zenit (i+j0) = zenit (j0)
                stratz(i+j0) = stratz(j0)
            end do
        end do
    end subroutine sol_oz

    subroutine solar(tyear, csol, lats, topsr)
        ! Average daily flux of solar radiation, from Hartmann (1994)

        real, intent(in) :: tyear
        real, intent(in) :: csol
        real, intent(in) :: lats(:)
        real, intent(out) :: topsr(size(lats))

        integer :: nlat, j
        real :: clat(size(lats)), slat(size(lats))
        real :: ca1, ca2, ca3, cdecl, ch0, csolp, decl, fdis, h0, alpha, pigr, sa1
        real :: sa2, sa3, sdecl, sh0, tdecl

        nlat = size(lats)

        ! Compute cosine and sine of latitude
        clat = cos(lats*pi/180.0)
        slat = sin(lats*pi/180.0)

        ! 1. Compute declination angle and Earth-Sun distance factor
        pigr  = 2.0*asin(1.0)
        alpha = 2.0*pigr*tyear

        ca1 = cos(alpha)
        sa1 = sin(alpha)
        ca2 = ca1*ca1-sa1*sa1
        sa2 = 2.*sa1*ca1
        ca3 = ca1*ca2-sa1*sa2
        sa3 = sa1*ca2+sa2*ca1

        decl = 0.006918 - 0.399912*ca1 + 0.070257*sa1 - 0.006758*ca2 + 0.000907*sa2&
            & - 0.002697*ca3 + 0.001480*sa3

        fdis = 1.000110 + 0.034221*ca1 + 0.001280*sa1 + 0.000719*ca2 + 0.000077*sa2

        cdecl = cos(decl)
        sdecl = sin(decl)
        tdecl = sdecl/cdecl

        ! 2. Compute daily-average insolation at the atm. top
        csolp = csol/pigr

        do j = 1, nlat
            ch0 = min(1.0, max(-1.0, -tdecl*slat(j)/clat(j)))
            h0  = acos(ch0)
            sh0 = sin(h0)

            topsr(j) = csolp*fdis*(h0*slat(j)*sdecl + sh0*clat(j)*cdecl)
        end do
    end subroutine solar

    subroutine solar_diurnal(tyear, rday, csol, lats, nlon, topsr)
        !-----------------------------------------------------------------------
        ! From ECBILT-CLIO model
        ! Calculates incoming solar radiation as a function of latitude
        ! for each day of the year, given the orbital parameters (see PMIP)
        ! One year has 360 days.  Reference: A. Berger, JAS, 35, 2362-2367,1978
        !-----------------------------------------------------------------------

        real, intent(in) :: tyear
        real, intent(in) :: rday
        real, intent(in) :: csol
        real, intent(in) :: lats(:)
        integer, intent(in) :: nlon
        real, intent(out) :: topsr(nlon,size(lats))

        integer :: i, j, l, NVE, nlat
        real :: clat(size(lats)), slat(size(lats))

        real :: beta, alam, alam0, ala, ala0
        real :: fac1, fac2, fac3, ro
        real :: deltal, sindl, cosdl, tandl
        real :: rkosz, rkosha1, ha1
        real :: deg2rad, day2rad
        real :: solard, solarf, solarcf, rlon, omweb, obl, hangle, ecc, cmu

        nlat = size(lats)

        ! Compute cosine and sine of latitude
        clat = cos(lats*pi/180.0)
        slat = sin(lats*pi/180.0)

        deg2rad = 2.0*asin(1.0)/180.0
        day2rad = 2.0*asin(1.0)/180.0

        ! Present-day orbital parameters: eccentricity ecc, obliquity obl and
        ! angle om between Vernal Equinox and Perihelion (angles all given
        ! in degrees and converted to radians). Solarc is the solar constant.
        ! NVE is day of the Vernal Equinox, set at 21 MARCH
        ! Implementatie van Nanne

        !mbp_s
        ecc = 0.016724
        !      ecc=0.016724*3.
        !mbp_e
        obl = 23.446*deg2rad
        !mbp_s
        omweb = (102.04 + 180.00)*deg2rad
        !      omweb=(102.04+180.00-180.0)*deg2rad
        !mbp_e
        !      csol=1365.
        NVE = 30 + 30 + 21

        ! In old SW-routine of ECBilt-model values were as follows:
        !
        !      ecc=0.0
        !      csol=1353.
        !      NVE=90
        !
        ! At 6000 years BP (Before Present) values were as follows:
        !
        !      ecc=0.018682
        !      obl=24.105*deg2rad
        !      omweb=(0.87+180.00)*deg2rad
        !
        ! :LGM:
        !     ecc=0.018994
        !     obl=22.949*deg2rad
        !     omweb=(114.42+180.00)*deg2rad

        ! First compute alam0 (the starting point). Then convert days to
        ! the true longitude ala, using Berger's formulae.
        ! Longitude ala loops from 0 (Vernal Equinox) to 359, ro is earth-sun
        ! distance relative to the major axis, del is declination.
        ala0 = 0.0
        beta = (1.0 - ecc**2.0)**0.5
        fac1 = (0.5*ecc + 0.125*ecc**3.0)*(1.0 + beta)*sin(ala0 - omweb)
        fac2 = 0.25*ecc**2.0*(0.5 + beta)*sin(2.0*(ala0 - omweb))
        fac3 = 0.125*ecc**3.0*(1.0/3.0 + beta)*sin(3.0*(ala0 -omweb))
        alam0 = ala0-2.0*(fac1 - fac2 + fac3)

        !      l=(imonth-1)*30+iday-NVE
        l = int(tyear*360) - NVE

        if (l < 0) l = l + 360
        alam = alam0 + l*1.0*2.0*asin(1.0)/180.0

        fac1 = (2.0*ecc - 0.25*ecc**3.0)*sin(alam - omweb)
        fac2 = 1.25*ecc**2.0*sin(2.0*(alam - omweb))
        fac3 = (13.0/12.0)*ecc**3.0*sin(3.0*(alam - omweb))
        ala = alam + fac1 + fac2 + fac3
        ro = (1.0 - ecc**2.0)/(1.0 + ecc*cos(ala - omweb))
        deltal = asin(sin(obl)*sin(ala))

        sindl = sin(deltal)
        cosdl = cos(deltal)
        tandl = tan(deltal)

        ! factor voor variable afstand Aarde-Zon (Berger, p.2367; Velds, p. 99)
        solard = 1.0/ro**2.0
        solarcf = csol
        ! beide effecten samen
        solarf = solarcf*solard

        do i = 1, nlon
            do j = 1, nlat
                rkosha1 = -(slat(j)/clat(j))*tandl
                rkosha1 = sign(min(abs(rkosha1), 1.0), rkosha1)
                ha1 = acos(rkosha1)
                rlon = float(i)/float(nlon)
                hangle = 2.0*acos(-1.0)*(tyear + rday + rlon)
                cmu = slat(j)*sindl - clat(j)*cosdl*cos(hangle)
                !fk
                !fk  NCAR Technical Note NCAR/TN-485+STR (CAM4.0), page 116
                !fk  note that only if cmu pos it is counted (Wiki pedia
                !fk  Solar Irradiance, https://en.wikipedia.org/wiki/Solar_irradiance)
                !fk
                !fk   daily mean formulae
                !fk
                !fk          rkosz=slat(j)*sindl*(ha1 - tan(ha1))/(2.*asin(1.))
                !fk   daily variing formulae
                !fk
                rkosz = max(cmu,0.0)
                topsr(i,j) = rkosz*csol*solard
            end do
        end do
    end subroutine solar_diurnal

    subroutine cloud(prog_q, rh, cnv_prec, precls, cnv_top, gse, fmask, &
                   & ngp, nlev, &
                   & icltop, cloudc, clstr, qcloud)
        real, intent(in) :: prog_q(ngp,nlev)
        real, intent(in) :: rh(ngp,nlev)
        real, intent(in) :: cnv_prec(ngp)
        real, intent(in) :: precls(ngp)
        integer, intent(in) :: cnv_top(ngp)
        real, intent(in) :: gse(ngp)
        real, intent(in) :: fmask(ngp)
        integer, intent(in) :: ngp
        integer, intent(in) :: nlev
        integer, intent(out) :: icltop(ngp)
        real, intent(out) :: cloudc(ngp)
        real, intent(out) :: clstr(ngp)
        real, intent(out) :: qcloud(ngp)

        integer :: nl1, nlp, j, k
        real :: clfact, clstrl, drh, fstab, pr1, rgse, rrcl

        nl1  = nlev - 1
        nlp  = nlev + 1
        rrcl = 1.0/(rhcl2 - rhcl1)

        ! 1. Cloud cover, defined as the sum of:
        !  - a term proportional to the square-root of precip. rate
        !  - a quadratic function of the max. relative humidity
        !    in tropospheric layers above pbl where q > prog_qcl :
        !     ( = 0 for rhmax < rhcl1, = 1 for rhmax > rhcl2 )
        !    Cloud-top level: defined as the highest (i.e. least sigma)
        !      between the top of convection/condensation and
        !      the level of maximum relative humidity.
        do j = 1, ngp
            if (rh(j,nl1) > rhcl1) then
                cloudc(j) = rh(j,nl1) - rhcl1
                icltop(j) = nl1
            else
                cloudc(j) = 0.0
                icltop(j) = nlp
            endif
        end do

        do k = 3, nlev - 2
            do j = 1, ngp
                drh = rh(j,k) - rhcl1
                if (drh > cloudc(j) .and. prog_q(j,k) > qcl) then
                    cloudc(j) = drh
                    icltop(j) = k
                endif
            end do
        end do

        do j = 1, ngp
            pr1 = min(pmaxcl, 86.4*(cnv_prec(j) + precls(j)))
            cloudc(j) = min(1.0, wpcl*sqrt(pr1) + min(1.0, cloudc(j)*rrcl)**2.0)
            icltop(j) = min(cnv_top(j), icltop(j))
        end do

        ! 2.  Equivalent specific humidity of clouds
        do j = 1, ngp
            qcloud(j) = prog_q(j,nl1)
        end do

        ! 3. Stratiform clouds at the top of pbl
        clfact = 1.2
        rgse   = 1.0/(gse_s1 - gse_s0)

        do j = 1, ngp
            ! Stratocumulus clouds over sea
            fstab = max(0.0, min(1.0, rgse*(gse(j) - gse_s0)))
            clstr(j) = fstab*max(clsmax - clfact*cloudc(j), 0.0)

            ! Stratocumulus clouds over land
            clstrl   = max(clstr(j), clsminl)*rh(j,nlev)
            clstr(j) = clstr(j) + fmask(j)*(clstrl - clstr(j))
        end do
    end subroutine cloud

    subroutine radsw(prog_sp, prog_q, icltop, cloudc, clstr, ozupp, ozone, zenit, stratz, fsol, qcloud, albsfc, &
                   & ngp, nlev, sig, dsig, &
                   & ssrd, ssr, tsr, tau2, tend_t_rsw)
        real, intent(in) :: prog_sp(ngp)
        real, intent(in) :: prog_q(ngp,nlev)
        integer, intent(in) :: icltop(ngp)
        real, intent(in) :: cloudc(ngp)
        real, intent(in) :: clstr(ngp)
        real, intent(in) :: ozupp(ngp)
        real, intent(in) :: ozone(ngp)
        real, intent(in) :: zenit(ngp)
        real, intent(in) :: stratz(ngp)
        real, intent(in) :: fsol(ngp)
        real, intent(in) :: qcloud(ngp)
        real, intent(in) :: albsfc(ngp)
        integer, intent(in) :: ngp
        integer, intent(in) :: nlev
        real, intent(in) :: sig(nlev)
        real, intent(in) :: dsig(nlev)
        real, intent(out) :: ssrd(ngp)
        real, intent(out) :: ssr(ngp)
        real, intent(out) :: tsr(ngp)
        real, intent(out) :: tau2(ngp,nlev,4)
        real, intent(out) :: tend_t_rsw(ngp,nlev)
!
!      equivalence (frefl(1,1),tau2(1,1,3))
        integer :: nl1, j, k
        real :: fband1, fband2
        real :: acloud(ngp), psaz(ngp), abs1, deltap, flux(ngp,2), eps1, stratc(ngp,2), acloud1

        nl1 = nlev - 1

        fband2 = 0.05
        fband1 = 1.0 - fband2

        ! 1. Initialization
        tau2 = 0.0

        do j = 1, ngp
            ! Change to ensure only ICLTOP <= NLEV used
            if (icltop(j) <= nlev) then
                tau2(j,icltop(j),3) = albcl*cloudc(j)
            end if
            ! end change
            tau2(j,nlev,3) = albcls*clstr(j)
        end do

        ! 2. Shortwave transmissivity:
        ! function of layer mass, ozone (in the statosphere),
        ! abs. humidity and cloud cover (in the troposphere)
        do j = 1, ngp
            psaz(j) = prog_sp(j)*zenit(j)
            acloud(j) = cloudc(j)*min(abscl1*qcloud(j), abscl2)
        end do

        do j = 1, ngp
            deltap = psaz(j)*dsig(1)
            tau2(j,1,1) = exp(-deltap*absdry)
        end do

        do k = 2, nl1
            abs1 = absdry + absaer*sig(k)**2
            do j = 1, ngp
                deltap = psaz(j)*dsig(k)
                if (k >= icltop(j)) then
                    tau2(j,k,1) = exp(-deltap*(abs1 + abswv1*prog_q(j,k) + acloud(j)))
                else
                    tau2(j,k,1) = exp(-deltap*(abs1 + abswv1*prog_q(j,k)))
                end if
            end do
        end do

        abs1 = absdry + absaer*sig(nlev)**2
        do j = 1, ngp
            deltap = psaz(j)*dsig(nlev)
            tau2(j,nlev,1) = exp(-deltap*(abs1 + abswv1*prog_q(j,nlev)))
        end do

        do k = 2, nlev
            do j = 1, ngp
                deltap = psaz(j)*dsig(k)
                tau2(j,k,2) = exp(-deltap*abswv2*prog_q(j,k))
            end do
        end do

        ! 3. Shortwave downward flux
        ! 3.1 Initialization of fluxes
        do j = 1, ngp
            tsr(j)    = fsol(j)
            flux(j,1) = fsol(j)*fband1
            flux(j,2) = fsol(j)*fband2
        end do

        ! 3.2 Ozone and dry-air absorption in the stratosphere
        k = 1
        do j = 1, ngp
            tend_t_rsw(j,k) = flux(j,1)
            flux (j,1) = tau2(j,k,1)*(flux(j,1) - ozupp(j)*prog_sp(j))
            tend_t_rsw(j,k) = tend_t_rsw(j,k) - flux(j,1)
        end do

        k = 2
        do j = 1, ngp
            tend_t_rsw(j,k) = flux(j,1)
            flux (j,1) = tau2(j,k,1)*(flux(j,1) - ozone(j)*prog_sp(j))
            tend_t_rsw(j,k) = tend_t_rsw(j,k) - flux(j,1)
        end do

        ! 3.3  Absorption and reflection in the troposphere
        do k = 3, nlev
            do j = 1, ngp
                tau2(j,k,3)      = flux(j,1)*tau2(j,k,3)
                flux (j,1)      = flux(j,1) - tau2(j,k,3)
                tend_t_rsw(j,k) = flux(j,1)
                flux (j,1)      = tau2(j,k,1)*flux(j,1)
                tend_t_rsw(j,k) = tend_t_rsw(j,k) - flux(j,1)
            end do
        end do

        do k = 2, nlev
            do j = 1, ngp
                tend_t_rsw(j,k) = tend_t_rsw(j,k) + flux(j,2)
                flux (j,2)      = tau2(j,k,2)*flux(j,2)
                tend_t_rsw(j,k) = tend_t_rsw(j,k) - flux(j,2)
            end do
        end do


        ! 4. Shortwave upward flux
        ! 4.1  Absorption and reflection at the surface
        do j = 1, ngp
            ssrd(j)   = flux(j,1)+flux(j,2)
            flux(j,1) = flux(j,1)*albsfc(j)
            ssr(j)    = ssrd(j) - flux(j,1)
        end do

        ! 4.2  Absorption of upward flux
        do k = nlev, 1, -1
            do j = 1, ngp
                tend_t_rsw(j,k) = tend_t_rsw(j,k) + flux(j,1)
                flux (j,1)      = tau2(j,k,1)*flux(j,1)
                tend_t_rsw(j,k) = tend_t_rsw(j,k) - flux(j,1)
                flux (j,1)      = flux(j,1) + tau2(j,k,3)
            end do
        end do

        ! 4.3  Net solar radiation = incoming - outgoing
        do j = 1, ngp
            tsr(j) = tsr(j) - flux(j,1)
        end do

        ! 5.  Initialization of longwave radiation model
        ! 5.1  Longwave transmissivity:
        !      function of layer mass, abs. humidity and cloud cover.

        !    Cloud-free levels (stratosphere + PBL)
        k = 1
        do j=1,ngp
            deltap = prog_sp(j)*dsig(k)
            tau2(j,k,1) = exp(-deltap*ablwin)
            tau2(j,k,2) = exp(-deltap*ablco2)
            tau2(j,k,3) = 1.0
            tau2(j,k,4) = 1.0
        end do

        do k = 2, nlev, nlev - 2
            do j = 1, ngp
                deltap      = prog_sp(j)*dsig(k)
                tau2(j,k,1) = exp(-deltap*ablwin)
                tau2(j,k,2) = exp(-deltap*ablco2)
                tau2(j,k,3) = exp(-deltap*ablwv1*prog_q(j,k))
                tau2(j,k,4) = exp(-deltap*ablwv2*prog_q(j,k))
            end do
        end do

        ! Cloudy layers (free troposphere)
        do j = 1, ngp
            acloud(j) = cloudc(j)*ablcl2
        end do

        do k = 3, nl1
            do j = 1, ngp
                deltap = prog_sp(j)*dsig(k)
                if (k < icltop(j)) then
                    acloud1 = acloud(j)
                else
                    acloud1 = ablcl1*cloudc(j)
                endif
                tau2(j,k,1) = exp(-deltap*(ablwin+acloud1))
                tau2(j,k,2) = exp(-deltap*ablco2)
                tau2(j,k,3) = exp(-deltap*max(ablwv1*prog_q(j,k),acloud(j)))
                tau2(j,k,4) = exp(-deltap*max(ablwv2*prog_q(j,k),acloud(j)))
            end do
        end do

        ! 5.2  Stratospheric correction terms
        eps1 = epslw/(dsig(1) + dsig(2))
        do j = 1, ngp
            stratc(j,1) = stratz(j)*prog_sp(j)
            stratc(j,2) = eps1*prog_sp(j)
        end do
    end subroutine radsw

!      SUBROUTINE RADLW (IMODE,TA,TS,
!     &                  FSFCD,FSFCU,
!     &                  FSFC,FTOP,DFABS)
!
!C--
!C--   SUBROUTINE RADLW (IMODE,TA,TS,
!C--  &                  FSFCD,FSFCU,
!C--  &                  FSFC,FTOP,DFABS)
!C--
!C--   Purpose: Compute the absorption of longwave radiation
!C--   Input:   IMODE  = index for operation mode 
!C--                     -1 : downward flux only
!C--                      0 : downward + upward flux 
!C--                     +1 : upward flux only
!C--            TA     = absolute temperature (3-dim)
!C--            TS     = surface temperature                    [if IMODE=0]
!C--            FSFCD  = downward flux of lw rad. at the sfc.   [if IMODE=1]
!C--            FSFCU  = surface blackbody emission (upward)    [if IMODE=1]
!C--            DFABS  = DFABS output from RADLW(-1,... )       [if IMODE=1]
!C--   Output:  FSFCD  = downward flux of lw rad. at the sfc.[if IMODE=-1,0]
!C--            FSFCU  = surface blackbody emission (upward)  [if IMODE=  0]
!C--            FSFC   = net upw. flux of lw rad. at the sfc. [if IMODE=0,1]
!C--            FTOP   = outgoing flux of lw rad. at the top  [if IMODE=0,1]
!C--            DFABS  = flux of lw rad. absorbed by each atm. layer (3-dim)
!C--
!C     Resolution parameters
!C
!      include "atparam.h"
!      include "atparam1.h"
!C
!      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )
!
!C     Number of radiation bands with tau < 1
!      PARAMETER ( NBAND=4 )
!
!C     Constants + functions of sigma and latitude
!
!      include "com_physcon.h"
!
!C     Radiation parameters
!
!      include "com_radcon.h"
!
!      REAL TA(NGP,NLEV), TS(NGP)
!
!      REAL FSFCD(NGP), FSFCU(NGP)
!
!      REAL FTOP(NGP), FSFC(NGP), DFABS(NGP,NLEV)
!
!C
!      NL1=NLEV-1
!
!      REFSFC=1.-EMISFC
!
!      IF (IMODE.EQ.1) GO TO 410
!
!C---  1. Blackbody emission from atmospheric levels.
!C        The linearized gradient of the blakbody emission is computed
!C        from temperatures at layer boundaries, which are interpolated 
!C        assuming a linear dependence of T on log_sigma.
!C        Above the first (top) level, the atmosphere is assumed isothermal.
!
!C     Temperature at level boundaries
!      DO K=1,NL1
!       DO J=1,NGP
!         ST4A(J,K,1)=TA(J,K)+WVI(K,2)*(TA(J,K+1)-TA(J,K))
!       ENDDO
!      ENDDO 
!
!C     Mean temperature in stratospheric layers 
!      DO J=1,NGP
!        ST4A(J,1,2)=0.75*TA(J,1)+0.25* ST4A(J,1,1)
!        ST4A(J,2,2)=0.50*TA(J,2)+0.25*(ST4A(J,1,1)+ST4A(J,2,1))
!      ENDDO
!
!C     Temperature gradient in tropospheric layers 
!
!      ANIS =1.0
!      ANISH=0.5*ANIS
!
!      DO K=3,NL1
!       DO J=1,NGP
!         ST4A(J,K,2)=ANISH*MAX(ST4A(J,K,1)-ST4A(J,K-1,1),0.)
!       ENDDO
!      ENDDO
!
!      DO J=1,NGP
!        ST4A(J,NLEV,2)=ANIS*MAX(TA(J,NLEV)-ST4A(J,NL1,1),0.)
!      ENDDO
!
!C     Blackbody emission in the stratosphere
!      DO K=1,2
!       DO J=1,NGP
!         ST4A(J,K,1)=SBC*ST4A(J,K,2)**4
!         ST4A(J,K,2)=0.
!       ENDDO
!      ENDDO
!
!C     Blackbody emission in the troposphere
!      DO K=3,NLEV
!       DO J=1,NGP
!         ST3A=SBC*TA(J,K)**3
!         ST4A(J,K,1)=ST3A*TA(J,K)
!         ST4A(J,K,2)=4.*ST3A*ST4A(J,K,2)
!       ENDDO
!      ENDDO
!
!
!C---  2. Initialization of fluxes
!
!      DO J=1,NGP
!        FSFCD(J)=0.
!      ENDDO
!
!      DO K=1,NLEV
!       DO J=1,NGP
!         DFABS(J,K)=0.
!       ENDDO
!      ENDDO
!
!C---  3. Emission ad absorption of longwave downward flux.
!C        For downward emission, a correction term depending on the      
!C        local temperature gradient and on the layer transmissivity is  
!C        added to the average (full-level) emission of each layer. 
!	
!C     3.1  Stratosphere
!
!      K=1
!      DO JB=1,2
!       DO J=1,NGP
!         EMIS=1.-TAU2(J,K,JB)
!         BRAD=FBAND(NINT(TA(J,K)),JB)*(ST4A(J,K,1)+EMIS*ST4A(J,K,2))
!         FLUX(J,JB)=EMIS*BRAD
!         DFABS(J,K)=DFABS(J,K)-FLUX(J,JB)
!       ENDDO
!      ENDDO
!
!      DO JB=3,NBAND
!       DO J=1,NGP
!         FLUX(J,JB)=0.
!       ENDDO
!      ENDDO
!	
!C     3.2  Troposphere
!
!      DO JB=1,NBAND
!       DO K=2,NLEV
!        DO J=1,NGP
!          EMIS=1.-TAU2(J,K,JB)
!          BRAD=FBAND(NINT(TA(J,K)),JB)*(ST4A(J,K,1)+EMIS*ST4A(J,K,2))
!          DFABS(J,K)=DFABS(J,K)+FLUX(J,JB)
!          FLUX(J,JB)=TAU2(J,K,JB)*FLUX(J,JB)+EMIS*BRAD
!          DFABS(J,K)=DFABS(J,K)-FLUX(J,JB)
!        ENDDO
!       ENDDO
!      ENDDO
!
!C     3.3 Surface downward flux
!
!      DO JB=1,NBAND
!       DO J=1,NGP
!         FSFCD(J)=FSFCD(J)+EMISFC*FLUX(J,JB)
!       ENDDO
!      ENDDO
!
!C     3.4 Correction for "black" band (incl. surface reflection)
!
!      EPS1=EPSLW*EMISFC
!      DO J=1,NGP
!        CORLW=EPS1*ST4A(J,NLEV,1)
!        DFABS(J,NLEV)=DFABS(J,NLEV)-CORLW
!        FSFCD(J)     =FSFCD(J)     +CORLW
!      ENDDO
!
!      IF (IMODE.EQ.-1) RETURN
!
!C---  4. Emission ad absorption of longwave upward flux. 
!C        For upward emission, a correction term depending on the      
!C        local temperature gradient and on the layer transmissivity is  
!C        subtracted from the average (full-level) emission of each layer. 
!	
!C     4.1  Surface
!
!C     Black-body (or grey-body) emission 
!      ESBC=EMISFC*SBC
!      DO J=1,NGP
!        TSQ=TS(J)*TS(J)
!        FSFCU(J)=ESBC*TSQ*TSQ
!      ENDDO
!
!C     Entry point for upward-only mode (IMODE=1)
! 410  CONTINUE
!
!      DO J=1,NGP
!        FSFC(J)=FSFCU(J)-FSFCD(J)
!      ENDDO
!
!      DO JB=1,NBAND
!       DO J=1,NGP
!         FLUX(J,JB)=FBAND(NINT(TS(J)),JB)*FSFCU(J)
!     &              +REFSFC*FLUX(J,JB)
!       ENDDO
!      ENDDO
!	
!C     4.2  Troposphere
!
!C     Correction for "black" band
!      DO J=1,NGP
!        DFABS(J,NLEV)=DFABS(J,NLEV)+EPSLW*FSFCU(J)
!      ENDDO
!
!      DO JB=1,NBAND
!       DO K=NLEV,2,-1
!        DO J=1,NGP
!          EMIS=1.-TAU2(J,K,JB)
!          BRAD=FBAND(NINT(TA(J,K)),JB)*(ST4A(J,K,1)-EMIS*ST4A(J,K,2))
!          DFABS(J,K)=DFABS(J,K)+FLUX(J,JB)
!          FLUX(J,JB)=TAU2(J,K,JB)*FLUX(J,JB)+EMIS*BRAD
!          DFABS(J,K)=DFABS(J,K)-FLUX(J,JB)
!        ENDDO
!       ENDDO
!      ENDDO
!	
!C     4.3  Stratosphere
!
!      K=1
!      DO JB=1,2
!       DO J=1,NGP
!         EMIS=1.-TAU2(J,K,JB)
!         BRAD=FBAND(NINT(TA(J,K)),JB)*(ST4A(J,K,1)-EMIS*ST4A(J,K,2))
!         DFABS(J,K)=DFABS(J,K)+FLUX(J,JB)
!         FLUX(J,JB)=TAU2(J,K,JB)*FLUX(J,JB)+EMIS*BRAD
!         DFABS(J,K)=DFABS(J,K)-FLUX(J,JB)
!       ENDDO
!      ENDDO
!
!C     Correction for "black" band and polar night cooling
!
!      DO J=1,NGP
!        CORLW1=DSIG(1)*STRATC(J,2)*ST4A(J,1,1)+STRATC(J,1)
!        CORLW2=DSIG(2)*STRATC(J,2)*ST4A(J,2,1)
!        DFABS(J,1)=DFABS(J,1)-CORLW1
!        DFABS(J,2)=DFABS(J,2)-CORLW2
!        FTOP(J)   =CORLW1+CORLW2
!      ENDDO
!
!C     4.4  Outgoing longwave radiation 
!
!      DO JB=1,NBAND
!       DO J=1,NGP
!         FTOP(J)=FTOP(J)+FLUX(J,JB)
!       ENDDO
!      ENDDO
!
!C---						
!      RETURN
!      END
!
!
!      SUBROUTINE RADSET
!C--
!C--   SUBROUTINE RADSET
!C--
!C--   Purpose: compute energy fractions in LW bands
!C--            as a function of temperature
!C--   Initialized common blocks: RADFIX
!
!C     Resolution parameters
!
!      include "atparam.h"
!      include "atparam1.h"
!
!      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )
!
!C     Radiation constants
!      include "com_radcon.h"
!
!      EPS1=1.-EPSLW
!
!      DO JTEMP=200,320
!        FBAND(JTEMP,2)=(0.148-3.0e-6*(JTEMP-247)**2)*EPS1
!        FBAND(JTEMP,3)=(0.356-5.2e-6*(JTEMP-282)**2)*EPS1
!        FBAND(JTEMP,4)=(0.314+1.0e-5*(JTEMP-315)**2)*EPS1
!        FBAND(JTEMP,1)=EPS1-(FBAND(JTEMP,2)+
!     &                       FBAND(JTEMP,3)+FBAND(JTEMP,4))
!      ENDDO
!
!      DO JB=1,4
!        DO JTEMP=100,199
!          FBAND(JTEMP,JB)=FBAND(200,JB)
!        ENDDO
!        DO JTEMP=321,400
!          FBAND(JTEMP,JB)=FBAND(320,JB)
!        ENDDO
!      ENDDO
!
!C--						
!      RETURN
!      END
end module radiation
