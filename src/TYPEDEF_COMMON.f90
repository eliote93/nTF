MODULE TYPEDEF_COMMON

USE XSLIB_MOD, ONLY : riftype
IMPLICIT NONE

TYPE superPin_Type
  INTEGER :: nx, ny, nxy
  INTEGER :: ix, iy, iFuelPin
  INTEGER :: NeighIdx(15), NeighSurfIdx(15)
  INTEGER, POINTER :: pin(:), pin2D(:, :)
  REAL :: BdLength(6), Center2SurfaceL(6), Area
  LOGICAL, POINTER :: lFuel(:)

  ! HEX
  INTEGER :: nNgh = 4, nBdmPin(6)
  INTEGER :: BdMPidx(2, 6), BdMPsuf(2, 6)
  INTEGER :: NghBd(15)
  REAL    :: NghLgh(15)
END TYPE

TYPE CoreXsMac_Type
  REAL, POINTER :: XSt(:, :)
  REAL, POINTER :: XStr(:, :)
  REAL, POINTER :: XSa(:, :)
  REAL, POINTER :: XSnf(:, :)
  REAL, POINTER :: XSkf(:, :)
  REAL, POINTER :: XSsm(:, :), XSsmP1(:, :), XSsmP2(:, :), XSsmP3(:, :)
END TYPE

TYPE BlockFxrDepl_Type
  INTEGER :: nFxrD, nResD, nCldD, nGenD, ntfsr
  INTEGER, POINTER :: mapfxrid(:), mapglobalid(:)
  INTEGER, POINTER :: mapR2G(:), mapG2R(:), mapGR2G(:), mapG2GR(:), mapC2G(:), mapG2C(:)
  INTEGER, POINTER :: MapNucl(:,:), MapNuclp13(:,:), MapFis(:,:), MapNuclRes(:,:)
  INTEGER, POINTER :: FPIsoIdx(:)

  LOGICAL, DIMENSION(:), POINTER :: lfuel, lCLD, lres, lGD, lh2o, lvoid, l_pnum_all

  INTEGER, DIMENSION(:), POINTER :: niso, niso_depl, niso_past, ndim, ipin, fsrst, nfsr, mapfsr
  REAL, DIMENSION(:), POINTER :: Burnup, Burnup_Past, Hmkg0, temp, area, h2ofrac, fsrvol

  INTEGER, POINTER :: itTemp(:,:), itTempPx(:,:)
  REAL, POINTER :: wtTemp(:,:), wtTempPx(:,:)

  INTEGER, DIMENSION(:,:), POINTER :: idiso, idiso_past
  REAL, DIMENSION(:,:), POINTER :: pnum, ind, pnum_past, pnum_all, pnum_past_all, chi
  REAL(4), DIMENSION(:,:), POINTER :: xseq_f_1g, xseq_c_1g, FnAdj, &
    fresoa, fresof, fresos, fresostr, fresonf, fresokf
  REAL(4), POINTER :: xseq_f_mg(:,:,:), xseq_c_mg(:,:,:), Ftadj(:,:,:), xseq(:,:,:,:), NDAF(:,:,:,:), &
    fresoAIso(:,:,:), fresoFIso(:,:,:), fresoSIso(:,:,:), fresoSSIso(:,:,:), fresoS1Iso(:,:,:)
END TYPE

TYPE BlockFxrTri_Type
  INTEGER :: nFxrTri, nResTri, nAICTri, nCldTri, nGenTri, NtTri, NtRTri, ntfsr, nisomax
  INTEGER, POINTER :: AccNt(:), AccNtR(:)
  INTEGER, POINTER :: mapfxrid(:), mapglobalid(:)
  INTEGER, POINTER :: mapR2G(:), mapG2R(:), mapGR2G(:), mapG2GR(:), mapC2G(:), mapG2C(:), mapA2G(:), mapG2A(:)
  INTEGER, POINTER :: MapNucl(:), MapNuclp13(:), MapNuclRes(:)

  LOGICAL, DIMENSION(:), POINTER :: lAIC, lCLD, lres, lh2o, lvoid

  INTEGER, DIMENSION(:), POINTER :: niso, ndim, ipin, fsrst, nfsr, mapfsr
  REAL, DIMENSION(:), POINTER :: temp, area, h2ofrac, fsrvol

  INTEGER, POINTER :: itTemp(:,:), itTempPx(:,:)
  REAL, POINTER :: wtTemp(:,:), wtTempPx(:,:)

  INTEGER, DIMENSION(:), POINTER :: idiso
  REAL, DIMENSION(:), POINTER :: pnum
  REAL(4), DIMENSION(:,:), POINTER :: xseq_f_1g, xseq_c_1g, FnAdj, &
    fresoa, fresof, fresos, fresostr, fresonf, fresokf
  REAL(4), DIMENSION(:), POINTER :: fresoAIso, fresoFIso, fresoSIso, fresoSSIso, fresoS1Iso
  REAL(4), POINTER :: xseq_f_mg(:,:,:), xseq_c_mg(:,:,:), Ftadj(:,:,:), xseq(:,:,:,:), NDAF(:,:,:,:)
END TYPE

TYPE GlobalPinInfo_Type
  LOGICAL, POINTER :: lres(:), lAIC(:)
  INTEGER, POINTER :: FxrStD(:), FxrStT(:), GenStD(:), GenStT(:), AICStT(:)
END TYPE

TYPE BlockResVarPin_Type
  LOGICAL, POINTER :: lres(:), lresA(:), lresC(:)
  REAL(4), POINTER :: avgxseq_mg(:,:,:), avgxseq_1g(:,:)
  REAL(4), POINTER :: rifa(:,:,:), rifs(:,:,:), riff(:,:,:)
  REAL(4), POINTER :: FnAdj(:,:)
  INTEGER, POINTER :: niso(:), igt(:), idiso(:,:), mapnucl(:,:)
  REAL, POINTER :: pnum(:,:), temp(:), ind(:,:)
END TYPE

TYPE BlockIsodata_Type
  INTEGER :: ntiso, ng
  INTEGER :: nttemp, nttempPx, Np0, Np1, Np2, Np3
  REAL, POINTER :: ptsTemp(:), ptsTempPx(:)
  INTEGER, POINTER :: ptrTemp(:), ptrTempPx(:)

  REAL, DIMENSION(:,:), POINTER :: siga, sigf, sigkf, signf, sigs, sigstr
  INTEGER, DIMENSION(:,:), POINTER :: gidSmp0, gidSmp1, gidSmp2, gidSmp3
  INTEGER, DIMENSION(:,:), POINTER :: ptrSmp0, ptrSmp1, ptrSmp2, ptrSmp3
  REAL, POINTER :: smp0(:), smp1(:), smp2(:), smp3(:)
END TYPE

TYPE BlockResIsoData_Type
  ! ----------- SubGrpFspNM ----------------------
  ! [Should contain]
  ! - lamsigp, sigp, lreso, nsig0, rtempsq, ri_a
  ! - ptrRTemp
  ! ----------------------------------------------
  ! ----------- EFFXXXXXX ------------------------
  ! [Not contain]
  ! - sigp, sig0sq, ri_a, lvabs_log
  ! ----------------------------------------------
  INTEGER :: ntrtemp, ntsig0, ntlv, nresiso, ntrif, ntrifRat, ntempPot, ntempLv, ntRifRx
  LOGICAL, POINTER :: lclad(:)
  INTEGER, POINTER :: ityp(:), rifid(:,:), IsoG2R(:)
  INTEGER, POINTER :: ptrRTemp(:), ptrSig0(:), ptrNlv(:), ptrRifRat(:)
  INTEGER, POINTER :: ptrTempPot(:), ptrTempLv(:), ptrRifRx(:)
  REAL, POINTER :: sigp(:), lamsigp(:,:) ! (1) & (igresb:igrese)
  REAL, POINTER :: sig0sq(:), rtempsq(:) ! (nsig0) & (nrtemp)
  REAL(4), DIMENSION(:,:), POINTER :: xsalog, ri_a ! (nsig0, igresb:igrese, nrtemp)
  REAL, DIMENSION(:,:), POINTER :: lvabs, lvfis, lvabs_log ! (nlv, igresb:igrese)
  REAL, DIMENSION(:,:), POINTER ::  wgtabs, wgtfis ! (nlv, nrtemp, igresb:igrese)
  REAL, POINTER :: ratlog(:)
  REAL(4), POINTER :: rifabs(:,:), riffis(:,:)
END TYPE

TYPE BlockMLGLib_Type
  INTEGER :: f_nmaclv, f_nmaclv1G, c_nmaclv1G
  REAL, POINTER, DIMENSION(:) :: c_maclv1G, c_maclv1G_log
  REAL, POINTER, DIMENSION(:,:) :: f_maclv, f_maclv1G
  REAL, POINTER, DIMENSION(:,:) :: f_maclv_log, f_maclv1G_log
END TYPE
END MODULE
