module physics

    implicit none

    private
    public get_physical_tendencies

contains

    !> Compute physical parametrization tendencies for u, v, t, q and add them
    !  to the dynamical grid-point tendencies
    subroutine get_physical_tendencies(prog_u, prog_v, prog_t, prog_q, prog_phi, prog_sp, &
                                     & ngp, nlev, nlat, sig, lats, tyear, lsm, &
                                     & tend_u, tend_v, tend_t, tend_q)
        ! Resolution parameters
        use physical_constants, only: cp
        use humidity, only: spec_hum_to_rel_hum
        use convection, only: get_convection_tendencies
        use condensation, only: get_large_scale_condensation_tendencies
        use radiation, only: sol_oz
        use utils, only: get_flux_conv_uvq

        ! Constants + functions of sigma and latitude
        !include "com_physcon.h"

        ! Model variables, tendencies and fluxes on gaussian grid
        !include "com_physvar2.h"

        ! Surface properties (time-inv.)
        !include "com_surfcon.h"

        ! Surface fields (daily averages)
        !include "com_cli_sea.h"
        !include "com_cli_land.h"
        !include "com_var_sea.h"
        !include "com_var_land.h"

        ! Logical and coupling flags
        !include "com_lflags.h"
        !include "com_cpl_flags.h"

        ! Date 
        !include "com_date.h"

        real, intent(in) :: prog_u(ngp,nlev)
        real, intent(in) :: prog_v(ngp,nlev)
        real, intent(in) :: prog_t(ngp,nlev)
        real, intent(in) :: prog_q(ngp,nlev)
        real, intent(in) :: prog_phi(ngp,nlev)
        real, intent(in) :: prog_sp(ngp)
        integer, intent(in) :: ngp
        integer, intent(in) :: nlev
        integer, intent(in) :: nlat
        real, intent(in) :: sig(nlev)
        real, intent(in) :: lats(nlat)
        real, intent(in) :: tyear
        real, intent(in) :: lsm(ngp)

        real, intent(out) :: tend_u(ngp,nlev)
        real, intent(out) :: tend_v(ngp,nlev)
        real, intent(out) :: tend_t(ngp,nlev)
        real, intent(out) :: tend_q(ngp,nlev)

        real :: prog_sp_inv(ngp)
        real :: stat_en(ngp,nlev)
        real :: rh(ngp,nlev)
        real :: qsat(ngp,nlev)

        integer :: nlon

        ! Intermediate variables for precipitation
        integer :: cnv_top(ngp)
        real :: cldbse_mss_flx(ngp)
        real :: cnv_prec(ngp)
        real :: tend_t_prc(ngp,nlev)
        real :: tend_q_prc(ngp,nlev)
        real :: precls(ngp)

        integer :: k
        integer :: icltop(ngp,2)
        integer :: icnv(ngp)
        real :: gse(ngp)

        ! Intermediate variables for radiation
        real :: topsr(nlat)

        ! Array for converting fluxes of prognostics into tendencies
        real :: grdsig(nlev)

        nlon = ngp/nlat

        grdsig = get_flux_conv_uvq(sig)

        ! Just make the tendencies equal to the prognostics divided by 10, as a test
        tend_u = prog_u/10.0
        tend_v = prog_v/10.0

        ! =========================================================================
        ! Compute thermodynamic variables
        ! =========================================================================

        ! Inverse of surface pressure
        prog_sp_inv = 1.0/prog_sp

        ! Dry static energy
        stat_en = cp*prog_t + prog_phi

        ! Calculate relative humidity from specific humidity (along with the
        ! saturation specific humidity)
        do k = 1, nlev
            call spec_hum_to_rel_hum(prog_t(:,k), prog_sp, sig(k), prog_q(:,k), &
                                   & rh(:,k), qsat(:,k))
        end do

        ! =========================================================================
        ! Precipitation
        ! =========================================================================

        ! Deep convection
        tend_t_prc = 0.0
        tend_q_prc = 0.0
        call get_convection_tendencies(prog_sp, stat_en, prog_q, qsat, &
                                     & ngp, nlev, sig, &
                                     & cnv_top, cldbse_mss_flx, cnv_prec, tend_t_prc, tend_q_prc)

        do k = 2, nlev
            tend_t(:,k) = tend_t_prc(:,k)*prog_sp_inv*grdsig(k)/cp
            tend_q(:,k) = tend_q_prc(:,k)*prog_sp_inv*grdsig(k)
        end do

        icnv(:) = nlev - cnv_top(:)

        ! Large-scale condensation
        tend_t_prc = 0.0
        tend_q_prc = 0.0
        call get_large_scale_condensation_tendencies(prog_sp, prog_q, qsat, cnv_top, &
                                                   & ngp, nlev, sig, &
                                                   & precls, tend_t_prc, tend_q_prc)

        do k = 2, nlev
            tend_t(:,k) = tend_t(:,k) + tend_t_prc(:,k)
            tend_q(:,k) = tend_q(:,k) + tend_q_prc(:,k)
        end do

        ! =========================================================================
        ! Radiation (shortwave and longwave) and surface fluxes
        ! =========================================================================

        ! Compute shortwave tendencies and initialize lw transmissivity
        ! The sw radiation may be called at selected time steps
!       if (lradsw) then
!            gse(:) = (stat_en(:,nlev-1) - stat_en(:,nlev))/(prog_phi(:,nlev-1) - prog_phi(:,nlev))
!
            call sol_oz(tyear, nlon, lats, topsr)
!
!            call cloud(qg1,rh,precnv,precls,iptop,gse,fmask1,icltop,cloudc,clstr)
!
!            cltop(:) = sigh(icltop(:,1)-1)*prog_sp(:)
!            prtop(:) = float(iptop(:))
!
!            call radsw(prog_sp,qg1,icltop,cloudc,clstr,ssrd,ssr,tsr,tt_rsw)
!
!            do k = 1, nlev
!                tt_rsw(:,k) = tt_rsw(:,k)*rps(:)*grdscp(k)
!            end do
!       end if

!C     3.2 Compute downward longwave fluxes 
!
!      CALL RADLW (-1,TG1,TS,
!     &            SLRD,SLRU(1,3),
!     &            SLR,OLR,TT_RLW)
!
!C     3.3. Compute surface fluxes and land skin temperature
!
!      if (iitest.eq.1) then 
!         print *, ' 3.3 in PHYPAR'
!         print *, 'mean(STL_AM) =', sum(STL_AM(:))/ngp
!         print *, 'mean(SST_AM) =', sum(SST_AM(:))/ngp
!      endif
!
!      CALL SUFLUX (PSG,UG1,VG1,TG1,QG1,RH,PHIG1,
!     &             PHIS0,FMASK1,STL_AM,SST_AM,SOILW_AM,SSRD,SLRD,
!     &             USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
!     &             TS,TSKIN,U0,V0,T0,Q0,.true.)
!     
!C--  3.3.1. Recompute sea fluxes in case of anomaly coupling
!
!      IF (ICSEA .GT. 0) THEN 
!
!         CALL SUFLUX (PSG,UG1,VG1,TG1,QG1,RH,PHIG1,
!     &             PHIS0,FMASK1,STL_AM,SSTI_OM,SOILW_AM,SSRD,SLRD,
!     &             USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
!     &             TS,TSKIN,U0,V0,T0,Q0,.false.)
!
!      ENDIF   
!
!C     3.4 Compute upward longwave fluxes, convert them to tendencies 
!C         and add shortwave tendencies
!
!      if (iitest.eq.1) print *, ' 3.4 in PHYPAR'
!
!      CALL RADLW (1,TG1,TS,
!     &            SLRD,SLRU(1,3),
!     &            SLR,OLR,TT_RLW)
!
!      DO K=1,NLEV
!        TT_RLW(:,K) = TT_RLW(:,K)*RPS(:)*GRDSCP(K)
!        TTEND (:,K) = TTEND(:,K)+TT_RSW(:,K)+TT_RLW(:,K)
!      ENDDO
!
!C--   4. PBL interactions with lower troposphere
!
!C     4.1 Vertical diffusion and shallow convection
!
!      CALL VDIFSC (UG1,VG1,SE,RH,QG1,QSAT,PHIG1,ICNV,
!     &             UT_PBL,VT_PBL,TT_PBL,QT_PBL)
!
!C     4.2 Add tendencies due to surface fluxes 
!
!      UT_PBL(:,NLEV) = UT_PBL(:,NLEV)+USTR(:,3)*RPS(:)*GRDSIG(NLEV)
!      VT_PBL(:,NLEV) = VT_PBL(:,NLEV)+VSTR(:,3)*RPS(:)*GRDSIG(NLEV)
!      TT_PBL(:,NLEV) = TT_PBL(:,NLEV)+ SHF(:,3)*RPS(:)*GRDSCP(NLEV)
!      QT_PBL(:,NLEV) = QT_PBL(:,NLEV)+EVAP(:,3)*RPS(:)*GRDSIG(NLEV)
!
!      DO K=1,NLEV
!        UTEND(:,K) = UTEND(:,K)+UT_PBL(:,K)
!        VTEND(:,K) = VTEND(:,K)+VT_PBL(:,K)
!        TTEND(:,K) = TTEND(:,K)+TT_PBL(:,K)
!        QTEND(:,K) = QTEND(:,K)+QT_PBL(:,K)
!      ENDDO
!
!C--   5. Store all fluxes for coupling and daily-mean output
!
!c      CALL DMFLUX (1)
    end subroutine get_physical_tendencies
end module physics

