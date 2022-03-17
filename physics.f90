module physics

    implicit none

    private
    public get_physical_tendencies

contains

    !> Compute physical parametrization tendencies for u, v, t, q and add them
    !  to the dynamical grid-point tendencies
    subroutine get_physical_tendencies(ug1, vg1, tg1, qg1, phig1, pslg1, utend, vtend, ttend, qtend)
        ! Resolution parameters
        use params, only: nlon, nlat, nlev, ngp

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

        real :: ug1(ngp,nlev)
        real :: vg1(ngp,nlev)
        real :: tg1(ngp,nlev)
        real :: qg1(ngp,nlev)
        real :: phig1(ngp,nlev)
        real :: pslg1(ngp)

        real :: utend(ngp,nlev)
        real :: vtend(ngp,nlev)
        real :: ttend(ngp,nlev)
        real :: qtend(ngp,nlev)

        integer :: iptop(ngp)
        integer :: icltop(ngp,2)
        integer :: icnv(ngp)
        real :: rps(ngp)
        real :: gse(ngp)

        integer :: iitest = 0

!C--   1. Compute sfc. pressure, dry static energy, relative humidity
!
!      if (iitest.eq.1) print *, ' 1. in PHYPAR'
!
!      PSG(:)=EXP(PSLG1(:))     
!      RPS(:)=1./PSG(:)
!
!      DO K=1,NLEV
!        SE(:,K)=CP*TG1(:,K)+PHIG1(:,K)
!      ENDDO
!
!      DO K=1,NLEV
!        CALL SHTORH (1,NGP,TG1(1,K),PSG,SIG(K),QG1(1,K),
!     &               RH(1,K),QSAT(1,K))
!      ENDDO
!
!C--   2. Precipitation 
!
!C     2.1 Deep convection
!
!      CALL CONVMF (PSG,SE,QG1,QSAT,
!     &             IPTOP,CBMF,PRECNV,TT_CNV,QT_CNV)
!
!      DO K=2,NLEV
!        TT_CNV(:,K) = TT_CNV(:,K)*RPS(:)*GRDSCP(K)
!        QT_CNV(:,K) = QT_CNV(:,K)*RPS(:)*GRDSIG(K)
!      ENDDO
!
!      ICNV(:) = NLEV-IPTOP(:)
!
!C     2.2 Large-scale condensation
!
!cfk#if !defined(KNMI)
!      CALL LSCOND (PSG,QG1,QSAT,
!     &             IPTOP,PRECLS,TT_LSC,QT_LSC)
!cfk#else
!cfk      CALL LSCOND (PSG,QG1,QSAT,TS,
!cfk     &             IPTOP,PRECLS,SNOWLS,TT_LSC,QT_LSC)
!cfk#endif
!
!      DO K=2,NLEV
!        TTEND(:,K) = TTEND(:,K)+TT_CNV(:,K)+TT_LSC(:,K)
!        QTEND(:,K) = QTEND(:,K)+QT_CNV(:,K)+QT_LSC(:,K)
!      ENDDO
!
!C--   3. Radiation (shortwave and longwave) and surface fluxes
!
!C     3.1 Compute shortwave tendencies and initialize lw transmissivity
!
!      if (iitest.eq.1) print *, ' 3.1 in PHYPAR'
!
!C     The sw radiation may be called at selected time steps
!
!      IF (LRADSW) THEN
!
!        GSE(:) = (SE(:,NLEV-1)-SE(:,NLEV))/
!     &           (PHIG1(:,NLEV-1)-PHIG1(:,NLEV))
!
!        CALL SOL_OZ (TYEAR,RDAY)
!
!        CALL CLOUD (QG1,RH,PRECNV,PRECLS,IPTOP,GSE,FMASK1,
!     &              ICLTOP,CLOUDC,CLSTR)
!
!        CLTOP(:) = SIGH(ICLTOP(:,1)-1)*PSG(:)
!        PRTOP(:) = float(IPTOP(:))
!
!        CALL RADSW (PSG,QG1,ICLTOP,CLOUDC,CLSTR,
!     &              SSRD,SSR,TSR,TT_RSW)
!
!        DO K=1,NLEV
!          TT_RSW(:,K)=TT_RSW(:,K)*RPS(:)*GRDSCP(K)
!        ENDDO
!
!      ENDIF
!
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
!
!C--   6. Random diabatic forcing 
!
!      IF (LRANDF) THEN
!
!C       6.1 Compute zonal-mean cross sections of diabatic forcing
!
!        IF (LRADSW) THEN
!          CALL XS_RDF (TT_LSC,TT_CNV,1)
!          CALL XS_RDF (TT_RSW,TT_RLW,2)
!        ENDIF
!
!C--     6.2 Compute and store 3-D pattern of random diabatic forcing
!        DO K=1,NLEV
!          TT_CNV(:,K) = TT_CNV(:,K)+TT_LSC(:,K)
!        ENDDO
!
!        CALL SETRDF (TT_LSC)
!
!        DO K=1,NLEV
!          TTEND(:,K) = TTEND(:,K)+TT_LSC(:,K)
!        ENDDO
!
!      ENDIF


!      include "phy_setrdf.f"
    end subroutine get_physical_tendencies
end module physics

