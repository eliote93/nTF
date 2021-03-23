#include <defines.h>
SUBROUTINE GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, Fxr, ifxr, Tempref, ig1, ig2, Core, ipin, iz, GroupInfo, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type, XsMac_Type,  FxrInfo_Type,  GroupInfo_Type,  PE_Type
USE CNTL,          ONLY : nTracerCntl_Type
USE TH_Mod,         ONLY : GetPinFuelTemp
USE MacXsLib_mod
USE BasicOperation, ONLY : CP_VA
USE GC_mod, ONLY : isosize
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type),pointer :: Fxr(:,:)
TYPE(FxrInfo_Type),pointer :: myFxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: ig1, ig2, ipin, iz, ifxr
REAL :: TempRef
REAL :: IsoMacXsSm(isosize,GroupInfo%ng,GroupInfo%ng)
REAL :: IsoMacXsSmP1(isosize,GroupInfo%ng,GroupInfo%ng)
!REAL :: IsoMacXsSm1(600,GroupInfo%ng,GroupInfo%ng)
!REAL,POINTER :: IsoMacXsSm(:,:,:)

INTEGER :: ig, igg, iresoGrp1, iresoGrp2, i, ng, niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)!, EQXS(:, :), SubGrpLv(:, :)
REAL :: Temp! TempSubGrp 
REAL :: SigFOld(1000,500), IsoXsMacStr(1000,500), MacSM(500, 500)!SigNfOld(1000,500)
TYPE(XSMac_Type), SAVE :: XsMacDat
LOGICAL :: FlagF

myFxr => Fxr(ifxr,iz)

IsoMacXsSm=0;IsoMacXsSmP1=0;
!ALLOCATE(IsoMacXsSm(200,GroupInfo%ng,GroupInfo%ng)) ! 14/02/27
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
ng = GroupInfo%ng
pnum => myFxr%pnum; idiso => myFxr%idiso; niso = myFxr%niso; 

CALL MacXsBase(XsMac, myFxr, ig1, ig2, ng, 1.0, FALSE, .TRUE.)
!Kaffa Fission
IF(myFxr%lres) THEN
  temp = myFxr%temp
  !FlagF = .TRUE.
  !IF (myFxr%lCLD) FlagF=.FALSE.
  !CALL CP_VA(SigFOld(1:niso, ig1:ig2), XsMac%IsoXsMacF(1:niso, ig1:ig2), niso, ig2-ig1+1)
  DO ig = ig1, ig2
    IF(ig .LT. iresogrp1 .OR. ig .GT. iresogrp2) CYCLE
    !Fission Cross Section 
    !DO i = 1, niso
    !  IF(SigFOld(i, ig) * XsMac%IsoXsMacF(i, ig) .NE. 0._8) THEN 
    !     XsMac%IsoXsMacNf(i, ig) = XsMac%IsoXsMacNf(i, ig) * XsMac%IsoXsMacF(i, ig) / SigFOld(i, ig)
    !     XsMac%IsoXsMacKF(i, ig) = XsMac%IsoXsMacKF(i, ig) * XsMac%IsoXsMacF(i, ig) / SigFOld(i, ig)
    !  ENDIF
    !ENDDO
    DO i = 1, niso
      XsMac%IsoXsMacA(i, ig) = XsMac%IsoXsMacA(i, ig) * myFxr%fresoAIso(i,ig)
      IF (XsMac%IsoXsMacNF(i, ig).ne.0)THEN
        XsMac%IsoXsMacF(i, ig) = XsMac%IsoXsMacF(i, ig) * myFxr%fresoFIso(i,ig)  
        XsMac%IsoXsMacNF(i, ig) = XsMac%IsoXsMacNF(i, ig) * myFxr%fresoFIso(i,ig)  
        XsMac%IsoXsMacKF(i, ig) = XsMac%IsoXsMacKF(i, ig) * myFxr%fresoFIso(i,ig)  
      ENDIF
    ENDDO   
  ENDDO
ENDIF
!Transport Cross Section
DO ig = ig1, ig2
  IF(ig .GE. iresogrp1 .AND. ig .LE. iresogrp2) THEN
    DO i = 1, niso
      XsMac%IsoXsMacTr(i, ig) = XsMac%IsoXsMaca(i, ig) + (XsMac%IsoXsMacS0(i, ig) - XsMac%IsoXsMacS1(i, ig))
    ENDDO
  ENDIF
ENDDO

!--- S0 matrix
DO i = 1, niso
  CALL MacXsScatMatrix(XsMacDat, myFxr%temp, 1, idiso(i:i), pnum(i:i), 1, ng, ng, GroupInfo, &
                     .FALSE.)
  DO ig = 1, ng
      DO igg = 1, ng
          isoMacXsSm(i,ig,igg)=XsMacDat%XsMacSm(ig,igg)
      ENDDO
  ENDDO
  !CALL CP_VA(IsoMacXsSm(i, 1:ng, 1:ng), XsMacDat%XsMacSm(1:ng, 1:ng), ng, ng)
ENDDO                            
                            
                            
!--- S1 matrix  14/06/01
#define p1sm
#ifdef p1sm
DO i = 1, niso
  CALL MacP1XsScatMatrix(XsMacDat, myFxr%temp, 1, idiso(i:i), pnum(i:i), 1, ng, ng, GroupInfo)
  DO ig = 1, ng
      DO igg = 1, ng
          isoMacXsSmP1(i,ig,igg)=XsMacDat%XsMacP1Sm(ig,igg)
      ENDDO
  ENDDO
ENDDO
#endif

!IF(myFxr%lres) THEN
!  IF (nTracerCntl%lRST) THEN
!    DO ig = ig1, ig2
!      IF(ig .LT. iresogrp1 .OR. ig .GT. iresogrp2) CYCLE
!      DO i = 1, niso
!          XsMac%IsoXsMacS0(i, ig) = XsMac%IsoXsMacS0(i, ig) * myFxr%fresoSIso(i,ig)
!          XsMac%IsoXsMacS1(i, ig) = XsMac%IsoXsMacS1(i, ig) * myFxr%fresoS1Iso(i,ig)
!      ENDDO
!    ENDDO
!  ENDIF
!ENDIF
!Scattering Cross section
!CALL MacXsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo, .FALSE., .FALSE.)
END SUBROUTINE

SUBROUTINE GenFxrMacXS(XsMac, MacXsSm, MacXsSmP1, Fxr, ifxr, Tempref, ig1, ig2, Core, ipin, iz, GroupInfo, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type, XsMac_Type,  FxrInfo_Type,  GroupInfo_Type,  PE_Type
USE CNTL,          ONLY : nTracerCntl_Type
USE MacXsLib_mod
USE BasicOperation, ONLY : CP_VA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type),pointer :: Fxr(:,:)
TYPE(FxrInfo_Type),pointer :: myFxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: ig1, ig2, ipin, iz, ifxr
REAL :: TempRef
REAL :: MacXsSm(GroupInfo%ng,GroupInfo%ng)
REAL :: MacXsSmP1(GroupInfo%ng,GroupInfo%ng)

INTEGER :: ig, igg, iresoGrp1, iresoGrp2, i, ng, niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)!, EQXS(:, :), SubGrpLv(:, :)
REAL :: Temp, siglp_cf!, TempSubGrp 
REAL :: SigFOld(300,47), IsoXsMacStr(300,47) !(niso,ng)
TYPE(XSMac_Type), SAVE :: XsMacDat
LOGICAL :: FlagF
REAL, POINTER :: fresoSIso(:,:), fresoSSIso(:,:), fresoS1Iso(:,:)

! myFxr => Fxr(ifxr,iz)
! 
! MacXsSm=0;MacXsSmP1=0;
! iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
! ng = GroupInfo%ng
! pnum => myFxr%pnum; idiso => myFxr%idiso; niso = myFxr%niso; 
! 
! CALL MacXsBase(XsMac, myFxr, ig1, ig2, ng, 1.0, FALSE, .TRUE.)
! !Kaffa Fission
! IF(myFxr%lres) THEN
!     temp = myFxr%temp
!     FlagF = .TRUE.
!     IF (myFxr%lCLD) FlagF=.FALSE.
!     
!     CALL CP_VA(SigFOld(1:niso, ig1:ig2), XsMac%IsoXsMacF(1:niso, ig1:ig2), niso, ig2-ig1+1)
!     DO ig = ig1, ig2
!       IF(ig .LT. iresogrp1 .OR. ig .GT. iresogrp2) CYCLE
! 
!       !Sigf, Siga <-self sheilding
!       IF (.not. FlagF) CALL GetFuelSigLP(siglp_cf,Core,Fxr,ipin,iz,ig)
!       CALL EffMacXs(XsMac, myFxr, TempRef, niso, ig, ng, siglp_cf, FlagF, TRUE)
!       !Fission Cross Section 
!       DO i = 1, niso
!         IF(SigFOld(i, ig) * XsMac%IsoXsMacF(i, ig) .NE. 0._8) THEN 
!            XsMac%IsoXsMacNf(i, ig) = XsMac%IsoXsMacNf(i, ig) * XsMac%IsoXsMacF(i, ig) / SigFOld(i, ig)
!            XsMac%IsoXsMacKf(i, ig) = XsMac%IsoXsMacKf(i, ig) * XsMac%IsoXsMacF(i, ig) / SigFOld(i, ig)
!         ENDIF
!       ENDDO
!     ENDDO
! ENDIF
! !Transport Cross Section
! DO ig = ig1, ig2
!     IF(ig .GE. iresogrp1 .AND. ig .LE. iresogrp2) THEN
!       DO i = 1, niso
!         XsMac%IsoXsMacTr(i, ig) = XsMac%IsoXsMaca(i, ig) + (XsMac%IsoXsMacS0(i, ig) - XsMac%IsoXsMacS1(i, ig))
!       ENDDO
!     ENDIF
! ENDDO
! 
! !--- S0 matrix
! IF(ntracercntl%lXsLib) Then
! fresoSIso => myFxr%fresoSIso; fresoSSIso => myFxr%fresoSSIso; fresoS1Iso => myFxr%fresoS1Iso
! CALL MacXsScatMatrix(XsMacDat, myFxr%temp, 1, idiso(1:niso), pnum(1:niso), 1, ng, ng, GroupInfo, &
!                      .FALSE., myFxr%lres, fresoSIso, fresoSSIso, fresoS1Iso)
! MacXsSm=XsMacDat%XsMacSm
! !--- S1 matrix  14/06/01
! #define p1sm
! #ifdef p1sm
! IF(nTRACERCntl%ScatOd .GE. 1)THEN
! CALL MacP1XsScatMatrix(XsMacDat, myFxr%temp, niso, idiso(1:niso), pnum(1:niso), 1, ng, ng, GroupInfo)
! XsMac%XsMacP1Sm=XsMacDat%XsMacP1Sm
! DO i = 1, niso
!   CALL MacP1XsScatMatrix(XsMacDat, myFxr%temp, 1, idiso(i:i), pnum(i:i), 1, ng, ng, GroupInfo)
!   DO ig = 1, ng
!       DO igg = 1, ng
!           MacXsSmP1(ig,igg)=MacXsSmP1(ig,igg)+XsMacDat%XsMacP1Sm(ig,igg)
!       ENDDO
!   ENDDO
! ENDDO
! ENDIF
! #endif
! !--- S2 matrix  16/06/09
! IF(nTRACERCntl%ScatOd .GE. 2)THEN
!     CALL MacP2XsScatMatrix(XsMacDat, myFxr%temp, niso, idiso(1:niso), pnum(1:niso), 1, ng, ng, GroupInfo)
!     XsMac%XsMacP2Sm=XsMacDat%XsMacP2Sm
! ENDIF
! !--- S3 matrix  16/06/09
! IF(nTRACERCntl%ScatOd .GE. 3)THEN
!     CALL MacP3XsScatMatrix(XsMacDat, myFxr%temp, niso, idiso(1:niso), pnum(1:niso), 1, ng, ng, GroupInfo)
!     XsMac%XsMacP3Sm=XsMacDat%XsMacP3Sm
! ENDIF
! ENDIF
! !Scattering Cross section
! !CALL MacXsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo, .FALSE., .FALSE.)
 END SUBROUTINE
