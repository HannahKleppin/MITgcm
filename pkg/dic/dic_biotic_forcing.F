#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: DIC_BIOTIC_FORCING

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_BIOTIC_FORCING(
     U                      PTR_DIC, PTR_ALK, PTR_PO4, PTR_DOP,
     U                      PTR_DET, PTR_PHY, PTR_ZOO
#ifdef ALLOW_O2
     U                      PTR_O2,
#endif
#ifdef ALLOW_FE
     U                      PTR_FE,
#endif
     I                      bi, bj, iMin, iMax, jMin, jMax,
     I                      myIter, myTime, myThid )

C !DESCRIPTION:
C updates all the tracers for the effects of air-sea exchange, biological
c activity and remineralization

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DIC_VARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"

C !INPUT/OUTPUT PARAMETERS: ===================================================
C  PTR_DIC              :: dissolced inorganic carbon
C  PTR_ALK              :: alkalinity
C  PTR_PO4              :: phosphate
c  PTR_DOP              :: dissolved organic phosphorous
C  PTR_DET              :: detritus
C  PTR_PHY              :: phytoplankton
C  PTR_ZOO              :: Zooplankton
c  PTR_O2               :: oxygen
C  PTR_FE               :: iron
c  bi, bj               :: current tile indices
C  myIter               :: current timestep
C  myTime               :: current time
C  myThid               :: thread number
      _RL  PTR_DIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_PO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_DOP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_PHY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_ZOO(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_DET(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#ifdef ALLOW_O2
      _RL  PTR_O2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
#ifdef ALLOW_FE
      _RL  PTR_FE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
      INTEGER bi, bj, iMin, iMax, jMin, jMax
      INTEGER myIter
      _RL myTime
      INTEGER myThid

#ifdef ALLOW_PTRACERS
#ifdef DIC_BIOTIC

C !LOCAL VARIABLES: ====================================================
C  i,j,k                  :: loop indices
C  G*                     :: tendency term for the tracers
C  SURA                   :: tendency of alkalinity due to freshwater
C  SURC                   :: tendency of DIC due to air-sea exchange
C                            and virtual flux
C  SURO                   :: tendency of O2 due to air-sea exchange
C  GPO4                   :: tendency of PO4 due to biological productivity,
C                            exchange with DOP pool/reminerization
C  GPHY                   :: tendency of Phytoplankton pool due to biological productivity
C                            (both new and regenerated), sinking out of the layer, 
C                            grazing of zooplankton, mortality and remineralization rate
C  GDET                   :: tendency of Detritus pool due to biological activity
C  GZOO                   :: tendency of zooplankton due to mortality, zooplankton
C                            grazing and loss rate to higher trophic levels 
C                            (is directly put back in to Detritus  pool)
C  GDOP                   :: tendency of DOP pool due to, detritus and phytoplankton degradation, 
C                            remineralization and regenerated biological productivity
C The 5 terms above also exist as _HAMOCC, with the corresponding formulations following Ilyina2013 (HAMOCC)
C  CAR                    :: carbonate changes due to biological
C                             productivity and remineralization
C  BIOac                  :: biological productivity
C  RDOP                   :: DOP sink due to remineralization
C  pflux                  :: changes to PO4 due to flux and remineralization
C  CAR_S                  :: carbonate sink
C  cflux                  :: carbonate changes due to flux and remineralization
C  freefe                 :: iron not bound to ligand
      _RL  GDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GPO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GDOP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GPHY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GDET(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GZOO(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GPO4_HAMOCC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GDOP_HAMOCC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GPHY_HAMOCC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GDET_HAMOCC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GZOO_HAMOCC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  SED_DET(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  SED_PHY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  sink_PHY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  sink_DET(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  SURA(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURO(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  CAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  BIOac(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  RDOP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  EX_PHY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  EX_ZOO(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  MORT_PHY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  MORT_ZOO(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  LTHL_ZOO(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GRAZ_ZOO(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GRAZ_ZOO_HAMOCC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  ASSIM_ZOO(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  ASSIM_ZOO_HAMOCC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  pflux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  exportflux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  CAR_S(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  cflux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#ifdef ALLOW_O2
      _RL  GO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
#ifdef ALLOW_FE
      _RL  GFE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  freefe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
      INTEGER i,j,k
#ifdef CAR_DISS
      INTEGER nCALCITEstep
#endif
#ifdef ALLOW_FE
# ifdef SEDFE
      INTEGER kBottom
# endif
#endif
CEOP

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('DIC_BIOTIC_FORCING',myThid)
#endif

      IF ( useThSIce .OR. useSEAICE .OR. useCoupler ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('DIC_FIELDS_UPDATE',myThid)
#endif
        CALL DIC_FIELDS_UPDATE(
     I                    bi, bj, myTime, myIter, myThid )
      ENDIF

       DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           RDOP(i,j,k) =0. _d 0
           GDIC(i,j,k) =0. _d 0
           GALK(i,j,k) =0. _d 0
           GPO4_HAMOCC(i,j,k) =0. _d 0
           GDOP_HAMOCC(i,j,k) =0. _d 0
           GPHY_HAMOCC(i,j,k) =0. _d 0
           GDET_HAMOCC(i,j,k) =0. _d 0
           GZOO_HAMOCC(i,j,k) =0. _d 0
           GPO4(i,j,k) =0. _d 0
           GDOP(i,j,k) =0. _d 0
           GPHY(i,j,k) =0. _d 0
           GDET(i,j,k) =0. _d 0
           GZOO(i,j,k) =0. _d 0
           CAR(i,j,k)  =0. _d 0
           BIOac(i,j,k)   =0. _d 0
           pflux(i,j,k)   =0. _d 0
           exportflux(i,j,k)=0. _d 0
           EXP_flux_PHY(i,j,k) = 0. _d 0
           EXP_fux_DET(i,j,k)= 0. _d 0
           cflux(i,j,k)   =0. _d 0
           CAR_S(i,j,k)   =0. _d 0
#ifdef ALLOW_O2
           GO2(i,j,k)     =0. _d 0
#endif
#ifdef ALLOW_FE
           GFE(i,j,k)     =0. _d 0
C  no longer needed after adding full initialisation of freefe in S/R FE_CHEM
c          freefe(i,j,k)  =0. _d 0
#endif
          ENDDO
         ENDDO
       ENDDO
       DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           SURA(i,j)   =0. _d 0
           SURC(i,j)   =0. _d 0
           SURO(i,j)   =0. _d 0
         ENDDO
       ENDDO

C carbon air-sea interaction
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('DIC_SURFFORCING',myThid)
#endif
       CALL DIC_SURFFORCING(
     I                    PTR_DIC, PTR_ALK, PTR_PO4,
     O                    SURC,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )

C alkalinity air-sea interaction
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('ALK_SURFFORCING',myThid)
#endif
       CALL ALK_SURFFORCING(
     I                    PTR_ALK,
     O                    SURA,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )

#ifdef ALLOW_O2
C oxygen air-sea interaction
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('O2_SURFFORCING',myThid)
#endif
       CALL O2_SURFFORCING(
     I                    PTR_O2,
     O                    SURO,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )
#endif

#ifdef ALLOW_FE
C find free iron
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('FE_CHEM',myThid)
#endif
       CALL FE_CHEM( bi, bj, iMin, iMax, jMin, jMax,
     U               PTR_FE,
     O               freefe,
     I               myIter, myThid )
#endif

C biological activity
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('BIO_EXPORT',myThid)
#endif
       CALL BIO_EXPORT(
     I               PTR_PO4, PTR_DOP
#ifdef ALLOW_FE
     I               PTR_FE,
#endif
     O               BIOac,
     I               bi, bj, iMin, iMax, jMin, jMax,
     I               myIter, myTime, myThid )

C flux of po4 from layers with biological activity
!C in our implementation particles are advected vertically, i.e.
!particles are sinking with explicit sinking velocity
!#ifdef ALLOW_DEBUG
!       IF (debugMode) CALL DEBUG_CALL('PHOS_FLUX',myThid)
!#endif
!       CALL PHOS_FLUX(
!     I               BIOac,
!     U               pflux, exportflux,
!     I               bi, bj, iMin, iMax, jMin, jMax,
!     I               myIter, myTime, myThid )

C- Carbonate sink
c- possibly sth needs still to be changed here?? Or only at carbonate flux? remove dop fraction
       DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
             CAR_S(i,j,k)=BIOac(i,j,k)*R_CP*rain_ratio(i,j,bi,bj)*
     &                    (1. _d 0-DOPfraction)
          ENDDO
         ENDDO
       ENDDO

C carbonate
#ifdef CAR_DISS
C dissolution only below saturation horizon
C code following method by Karsten Friis
        nCALCITEstep = 3600
        IF(myIter .lt. (nIter0+5) .or.
     &               mod(myIter,nCALCITEstep) .eq. 0)THEN
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL('CALCITE_SATURATION',myThid)
#endif
          CALL CALCITE_SATURATION(
     I                    PTR_DIC, PTR_ALK, PTR_PO4,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )
        ENDIF

#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('CAR_FLUX_OMEGA_TOP',myThid)
#endif
        CALL CAR_FLUX_OMEGA_TOP(
     I                    BIOac,
     O                    cflux,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )
#else
C old OCMIP way
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('CAR_FLUX',myThid)
#endif
        CALL CAR_FLUX(
     I                    CAR_S,
     U                    cflux,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )
#endif




!!!!!!!!!!!!!!!!!!!!! Active sinking of Phytoplankton and Detritus !!!!!!!!

C- Calculate sedimentation rate at base of each layer


       DO j=jmin,jmax
        DO i=imin,imax
         SED_DET(i,j,1)=0. _d 0
         SED_PHY(i,j,1)=0. _d 0
         SED_DET(i,j,nlev)=0. _d 0
         SED_PHY(i,j,nlev)=0. _d 0
         DO k=2,Nr-1
                  
             IF (hFacC(i,j,k,bi,bj).GT.0. _d 0) THEN
C--   If no layer below initial layer (because of bottom or
C--   topography), no flux
          IF (k.EQ.Nr) THEN
          SED_DET(i,j,k)=0. _d 0
           SED_PHY(i,j,k)=0. _d 0
          ELSEIF ( _hFacC(i,j,k+1,bi,bj).EQ.0. _d 0 ) THEN
           SED_DET(i,j,k)=0. _d 0
           SED_PHY(i,j,k)=0. _d 0
          ELSE
C- average to have concentration and sinking velocity at interface of layer        
           SED_PHY(i,j,k) = (PTR_PHY(i,j,k-1)*drF(k-1)*_hFacC(i,j,k-1,bi,bj)
    &       +PTR_PHY(i,j,k)*drF(k)* _hFacC(i,j,k,bi,bj))/(drF(k)*_hFacC(i,j,k,bi,bj)
    &         +drF(k-1)*_hFacC(i,j,k-1,bi,bj) )
    
           sink_PHY(i,j,k)=(w_sink(i,j,k-1)*drF(k-1)*_hFacC(i,j,k-1,bi,bj)
    &                      +w_sink(i,j,k)*drF(k)* _hFacC(i,j,k,bi,bj) )/
    &                  (drF(k)*_hFacC(i,j,k,bi,bj)+drF(k-1)*_hFacC(i,j,k-1,bi,bj) )
    
           SED_DET(i,j,k) = (PTR_DET(i,j,k-1)*drF(k-1)*_hFacC(i,j,k-1,bi,bj)
    &                    +PTR_DET(i,j,k)*drF(k)* _hFacC(i,j,k,bi,bj))/(drF(k)*
    &                    _hFacC(i,j,k,bi,bj)+drF(k-1)*_hFacC(i,j,k-1,bi,bj) )
    
           sink_DET(i,j,k)=(w_sink_d(i,j,k-1)*drF(k-1)*_hFacC(i,j,k-1,bi,bj)
    &                      +w_sink_d(i,j,k)*drF(k)* _hFacC(i,j,k,bi,bj) )/
    &                  (drF(k)*_hFacC(i,j,k,bi,bj)+drF(k-1)*_hFacC(i,j,k-1,bi,bj) )
    
    
          ENDIF
         ENDIF
     
        EXP_flux_PHY(i,j,k-1)=(SED_PHY(i,j,k-1)*sink_PHY(i,j,k-1)-SED_PHY(i,j,k)
    &                       *sink_PHY(i,j,k))/(drF(k-1)*_hFacC(i,j,k-1,bi,bj))
             
        EXP_flux_DET(i,j,k-1)=(SED_DET(i,j,k-1)*sink_DET(i,j,k-1)-SED_DET(i,j,k)
    &                       *sink_DET(i,j,k))/(drF(k-1)*_hFacC(i,j,k-1,bi,bj))
             
          ENDDO
        
        EXP_flux_PHY(i,j,Nr)=SED_PHY(i,j,Nr)*sink_PHY(i,j,Nr)
    &                        /(drF(Nr)*_hFacC(i,j,Nr,bi,bj))     
        EXP_flux_DET(i,j,Nr)=SED_DET(i,j,Nr)*sink_DET(i,j,Nr)
    &                        /(drF(Nr)*_hFacC(i,j,Nr,bi,bj))
        ENDDO
       ENDDO


C add all tendencies for PO4, DOP, ALK, DIC
       DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
#ifdef DIC_NO_NEG
           RDOP(i,j,k)= MAX(maskC(i,j,k,bi,bj)*KDOPRemin*PTR_DOP(i,j,k)
     &                     ,0. _d 0)
#else
           RDOP(i,j,k)= maskC(i,j,k,bi,bj)*KDOPRemin*PTR_DOP(i,j,k)
#endif
!C Source and sink terms for all biogeochemical tracers
           !Phytoplankton exudation to DOP 
           EX_PHY(i,j,k)= maskC(i,j,k,bi,bj)*r_PHY*PTR_PHY(i,j,k)
           !Zooplankton excretion to DOP (possibly a share could go to PO4)
           EX_ZOO(i,j,k)= maskC(i,j,k,bi,bj)*r_ZOO*PTR_ZOO(i,j,k)
           !remineralization of detritus to DOP
           RMIN_DET(i,j,k)= maskC(i,j,k,bi,bj)*K_remin_DET*PTR_DET(i,j,k)
           !phytoplankton mortality to Det
           MORT_PHY(i,j,k)= maskC(i,j,k,bi,bj)*mu_PHY*PTR_PHY(i,j,k)
           !Zooplankton mortality to Det
           MORT_ZOO(i,j,k)= maskC(i,j,k,bi,bj)*mu_ZOO*PTR_ZOO(i,j,k)
           !includes natural mortality and feeding of higher trophic levels
           MORT_ZOO_HAMOCC(i,j,k)= maskC(i,j,k,bi,bj)*mu_HAM_ZOO*
     &                            (PTR_ZOO(i,j,k)-ZOO_min)^2
           !Grazing of Zooplankton on PhytoP
           !Two different formulations for comparison, HAMOCC formulation has in all terms a minimum PHYtoplankton&ZOOplankton concentration that remains unaffected by grazing,
           !mortality, etc.
           GRAZ_ZOO(i,j,k)= maskC(i,j,k,bi,bj)*graz_ZOO_max*PTR_PHY(i,j,k)/
     &                      (PTR_PHY(i,j,k)+K_Zoo)*PTR_ZOO(i,j,k) 
           GRAZ_ZOO_HAMOCC(i,j,k)=maskC(i,j,k,bi,bj)*graz_ZOO_max*PTR_ZOO(i,j,k)*
     &                            (PTR_PHY(i,j,k)-PHY_MIN)/(K_Zoo+PTR_PHY(i,j,k))
           !Assimilation efficiency
           ASSIM_ZOO(i,j,k)= maskC(i,j,k,bi,bj)*GRAZ_ZOO*sigma_z
           ASSIM_ZOO_HAMOCC(i,j,k)= maskC(i,j,k,bi,bj)*GRAZ_ZOO_HAMOCC(i,j,k)
     &                             *sigma_assim*omega_assim
           !Loss of Zooplankton mass to higher trophic levels, returned to Det
           LTHL_ZOO(i,j,k)=maskC(i,j,k,bi,bj)*loss_HT_levels*PTR_ZOO(i,j,k)^2
!-----------------------------------------------------------------------------           
!C Tendency terms due to source and sinks of all biogeochemical tracers 
!Phosphate
           GPO4(i,j,k)=-BIOac(i,j,k)*PTR_PHY(i,j,k) + RDOP(i,j,k)+
     &                  RMIN_DET(i,j,k)
           GPO4_HAMOCC(i,j,k)=-BIOac(i,j,k)*PTR_PHY(i,j,k) + RDOP(i,j,k)+
     &                  RMIN_DET(i,j,k)+(GRAZ_ZOO_HAMOCC(i,j,k)*sigma_assim
     &                  -ASSIM_ZOO_HAMOCC(i,j,k))+MORT_ZOO_HAMOCC(i,j,k)*epsilon_z

!Dissolved organic P
           GDOP(i,j,k)=EX_PHY(i,j,k) + EX_ZOO(i,j,k)-RDOP(i,j,k)
     &                  +(GRAZ_ZOO(i,j,k)-ASSIM_ZOO(i,j,k))*omega_z
           GDOP_HAMOCC(i,j,k)= (EX_PHY(i,j,k)-r_PHY*PHY_MIN*maskC(i,j,k,bi,bj))
     &                        + (EX_ZOO(i,j,k)-r_ZOO*ZOO_min*maskC(i,j,k,bi,bj))
     &                        -RDOP(i,j,k)
!PHYTOPLANKTON
           GPHY(i,j,k) = (BIOac(i,j,k))*PTR_PHY(i,j,k)-EXP_flux_PHY(i,j,k)
     &                   -MORT_PHY(i,j,k)-EX_PHY(i,j,k)-GRAZ_ZOO(i,j,k)
           !HAMOCC PHYTO
           GPHY_HAMOCC(i,j,k) = (BIOac(i,j,k))*PTR_PHY(i,j,k)-EXP_flux_PHY(i,j,k)
     &                   -(MORT_PHY(i,j,k)-mu_PHY*PHY_MIN*maskC(i,j,k,bi,bj))
     &                   -(EX_PHY(i,j,k)-r_PHY*PHY_MIN*maskC(i,j,k,bi,bj))-GRAZ_ZOO_HAMOCC(i,j,k)
!DETRITUS
           GDET(i,j,k) = MORT_PHY(i,j,k)+ MORT_ZOO(i,j,k)-RMIN_DET(i,j,k)
     &                   +(GRAZ_ZOO(i,j,k)-ASSIM_ZOO(i,j,k))*(1. _d 0 -omega_z)+LTHL_ZOO(i,j,k)
     &                   - EXP_flux_DET(i,j,k)
           !HAMOCC DETRITUS
           GDET_HAMOCC(i,j,k)=MORT_PHY(i,j,k)+MORT_ZOO_HAMOCC(i,j,k)*(1. _d 0 -epsilon_z)
     &                        +(GRAZ_ZOO_HAMOCC(i,j,k)-ASSIM_ZOO_HAMOCC(i,j,k)/omega_assim)
     &                        -RMIN_DET(i,j,k)-EXP_flux_DET(i,j,k)
                             
!#           GDET_HAMOCC(i,j,k)=MORT_PHY(i,j,k)+MORT_ZOO(i,j,k)+GRAZ_ZOO_HAMOCC(i,j,k)
!     &                        -ASSIM_ZOO_HAMOCC(i,j,k)/sigma_z+LTHL_ZOO(i,j,k)*(1-0.95 _d 0)
!     &                         -RMIN_DET(i,j,k)-EXP_flux(i,j,k)
!Zooplankton
           GZOO_HAMOCC(i,j,k) = ASSIM_ZOO_HAMOCC(i,j,k)-MORT_ZOO_HAMOCC(i,j,k)
     &                          -(EX_ZOO(i,j,k)-r_ZOO*ZOO_min*maskC(i,j,k,bi,bj))
           GZOO(i,j,k) = ASSIM_ZOO(i,j,k)-MORT_ZOO(i,j,k)
     &                   -LTHL_ZOO(i,j,k)- EX_ZOO(i,j,k)
     
     
           car(i,j,k) = cflux(i,j,k) - CAR_S(i,j,k)
!Alkalinity
           GALK(i,j,k)=+2. _d 0 *car(i,j,k)-R_NP*GPO4(i,j,k)
!conversion to carbon
           GDIC(i,j,k)=car(i,j,k)+R_CP*GPO4(i,j,k)

#ifdef ALLOW_O2
           if (PTR_O2(i,j,k).GT.O2crit) then
             GO2(i,j,k)= R_OP*GPO4(i,j,k)
           else
             GO2(i,j,k)= 0. _d 0
           endif
#endif
#ifdef ALLOW_FE
           GFE(i,j,k) = R_FeP*GPO4(i,j,k)
     &                 -Kscav*freefe(i,j,k)
#endif
          ENDDO
         ENDDO
       ENDDO

       DO j=jMin,jMax
         DO i=iMin,iMax
               GALK(i,j,1)=GALK(i,j,1)+SURA(i,j)
               GDIC(i,j,1)=GDIC(i,j,1)+SURC(i,j)
#ifdef ALLOW_O2
               GO2(i,j,1) =GO2(i,j,1)+SURO(i,j)
#endif
#ifdef ALLOW_FE
               GFE(i,j,1)=GFE(i,j,1)+alpfe*
     &                    InputFe(i,j,bi,bj)*recip_drF(1)
     &                       *recip_hFacC(i,j,1,bi,bj)
# ifdef SEDFE
C include iron sediment source using the flux of po4 into bottom layer
               kBottom   = MAX(kLowC(i,j,bi,bj),1)
               GFE(i,j,kBottom)=GFE(i,j,kBottom)
     &              +( fesedflux_pcm*pflux(i,j,kBottom) + FeIntSec )
     &              *recip_drF(kBottom)*recip_hFacC(i,j,kBottom,bi,bj)
# endif
#endif
         ENDDO
       ENDDO

       IF ( useOBCS ) THEN
        DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
            GDIC(i,j,k) = GDIC(i,j,k)*maskInC(i,j,bi,bj)
            GALK(i,j,k) = GALK(i,j,k)*maskInC(i,j,bi,bj)
            GPO4(i,j,k) = GPO4(i,j,k)*maskInC(i,j,bi,bj)
            GDOP(i,j,k) = GDOP(i,j,k)*maskInC(i,j,bi,bj)
            GPHY(i,j,k) = GPHY(i,j,k)*maskInC(i,j,bi,bj)
            GDET(i,j,k) = GDET(i,j,k)*maskInC(i,j,bi,bj)
            GZOO(i,j,k) = GZOO(i,j,k)*maskInC(i,j,bi,bj)
#ifdef ALLOW_O2
            GO2(i,j,k)  =  GO2(i,j,k)*maskInC(i,j,bi,bj)
#endif
#ifdef ALLOW_FE
            GFE(i,j,k)  =  GFE(i,j,k)*maskInC(i,j,bi,bj)
#endif
          ENDDO
         ENDDO
        ENDDO
       ENDIF

C update
       DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
           PTR_DIC(i,j,k)=
     &      PTR_DIC(i,j,k)+GDIC(i,j,k)*PTRACERS_dTLev(k)
           PTR_ALK(i,j,k)=
     &      PTR_ALK(i,j,k)+GALK(i,j,k)*PTRACERS_dTLev(k)
           PTR_PO4(i,j,k)=
     &      PTR_PO4(i,j,k)+GPO4(i,j,k)*PTRACERS_dTLev(k)
           PTR_DOP(i,j,k)=
     &      PTR_DOP(i,j,k)+GDOP(i,j,k)*PTRACERS_dTLev(k)
           PTR_PHY(i,j,k)=
     &      PTR_PHY(i,j,k)+GPHY(i,j,k)*PTRACERS_dTLev(k)
           PTR_DET(i,j,k)=
     &      PTR_DET(i,j,k)+GDET(i,j,k)*PTRACERS_dTLev(k)
           PTR_ZOO(i,j,k)=
     &      PTR_ZOO(i,j,k)+GZOO(i,j,k)*PTRACERS_dTLev(k)
#ifdef ALLOW_O2
           PTR_O2(i,j,k)=
     &      PTR_O2(i,j,k)+GO2(i,j,k)*PTRACERS_dTLev(k)
#endif
#ifdef ALLOW_FE
           PTR_FE(i,j,k)=
     &      PTR_FE(i,j,k)+GFE(i,j,k)*PTRACERS_dTLev(k)
#endif
          ENDDO
         ENDDO
       ENDDO

       

#ifdef ALLOW_FE
#ifdef MINFE
c find free iron and get rid of insoluble part
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('FE_CHEM',myThid)
#endif
       CALL FE_CHEM( bi, bj, iMin, iMax, jMin, jMax,
     U               PTR_FE,
     O               freefe,
     I               myIter, myThid )
#endif
#endif

#ifdef ALLOW_TIMEAVE
C save averages
       IF ( PTRACERS_taveFreq.GT.0. ) THEN
        DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
            BIOave(i,j,k,bi,bj)   =BIOave(i,j,k,bi,bj)+
     &                             BIOac(i,j,k)*deltaTClock
            PHYave(i,j,k,bi,bj)   =PHYave(i,j,k,bi,bj)+
     &                             PTR_PHY(i,j,k)*deltaTClock
            ZOOave(i,j,k,bi,bj)   =ZOOave(i,j,k,bi,bj)+
     &                             PTR_ZOO(i,j,k)*deltaTClock
            DETave(i,j,k,bi,bj)   =DETave(i,j,k,bi,bj)+
     &                             PTR_DET(i,j,k)*deltaTClock
            DOPave(i,j,k,bi,bj)   =DOPave(i,j,k,bi,bj)+
     &                             PTR_DOP(i,j,k)*deltaTClock
            CARave(i,j,k,bi,bj)   =CARave(i,j,k,bi,bj)+
     &                             CAR(i,j,k)*deltaTClock
            OmegaCave(i,j,k,bi,bj)=OmegaCave(i,j,k,bi,bj)+
     &                             OmegaC(i,j,k,bi,bj)*deltaTClock
            pfluxave(i,j,k,bi,bj) =pfluxave(i,j,k,bi,bj) +
     &                             pflux(i,j,k)*deltaTClock
            epfluxave(i,j,k,bi,bj)=epfluxave(i,j,k,bi,bj) +
     &                             exportflux(i,j,k)*deltaTClock
            cfluxave(i,j,k,bi,bj) =cfluxave(i,j,k,bi,bj) +
     &                             cflux(i,j,k)*deltaTClock
          ENDDO
         ENDDO
        ENDDO
         DO j=jMin,jMax
          DO i=iMin,iMax
              SURave(i,j,bi,bj)    =SURave(i,j,bi,bj)+
     &                              SURC(i,j)*deltaTClock
#ifdef ALLOW_O2
              SUROave(i,j,bi,bj)   =SUROave(i,j,bi,bj)+
     &                              SURO(i,j)*deltaTClock
#endif
              pCO2ave(i,j,bi,bj)   =pCO2ave(i,j,bi,bj)+
     &                              pCO2(i,j,bi,bj)*deltaTClock
              pHave(i,j,bi,bj)     =pHave(i,j,bi,bj)+
     &                              pH(i,j,bi,bj)*deltaTClock
              fluxCO2ave(i,j,bi,bj)=fluxCO2ave(i,j,bi,bj)+
     &                           fluxCO2(i,j,bi,bj)*deltaTClock
          ENDDO
         ENDDO
         DIC_timeAve(bi,bj) = DIC_timeAve(bi,bj)+deltaTClock
       ENDIF
#endif /* ALLOW_TIMEAVE*/

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL DIAGNOSTICS_FILL(BIOac  ,'DICBIOA ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(CAR    ,'DICCARB ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pCO2   ,'DICPCO2 ',0,1 ,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(fluxCO2,'DICCFLX ',0,1 ,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pH     ,'DICPHAV ',0,1 ,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(SURC   ,'DICTFLX ',0,1 ,2,bi,bj,myThid)
#ifdef ALLOW_O2
        CALL DIAGNOSTICS_FILL(SURO   ,'DICOFLX ',0,1 ,2,bi,bj,myThid)
#endif
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('DIC_BIOTIC_FORCING',myThid)
#endif

#endif /* DIC_BIOTIC */
#endif /* ALLOW_PTRACERS */

       RETURN
       END
