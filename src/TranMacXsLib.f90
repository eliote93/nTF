MODULE TranMacXsLib_Mod
USE PARAM
USE TYPEDEF,         ONLY : XsMac_Type,     GroupInfo_Type,      Fxrinfo_type
USE XSLIB_MOD
USE XsUtil_mod,      ONLY : AllocXsMac,          AllocMacIsoXs,  &
                            LineIntPol,          XsTempInterpolation, P1XsTempInterpolation
USE BasicOperation,  ONLY : CP_CA, CP_VA,   MULTI_CA,            DotProduct 
USE IOUTIL,          ONLY : TERMINATE
!Default Kinetic Parameter for XS Library Calculation Mode
INTEGER, PARAMETER :: nprec_LibBase = 6
INTEGER, PARAMETER :: ng_LibBase = 47;

REAL :: DecayConst(6), DelaySpectrum(ng_LibBase)
REAL :: Alpha(6), Zeta(6) ! Coefficients for Decay Heat Calculation
DATA DecayConst / 0.0128_8, 0.0318_8, 0.119_8, 0.3181_8, 1.4027_8, 3.9286_8/
DATA DelaySpectrum / 4*0._8, 0.005_8, 0.021_8, 0.269_8, 0.247_8, 0.429_8, 0.029_8, 37*0._8/

DATA Alpha / 2.35402E-2, 1.89077E-2, 1.39236E-2, 6.90315E-3, 3.56888E-3, 3.31633E-3 /
DATA Zeta / 1.05345E-1, 8.37149E-3, 5.20337E-4, 4.73479E-5, 3.28153E-6, 1.17537E-11 /

CONTAINS

SUBROUTINE SetTranXsLibInfo(GroupInfo)
TYPE(GroupInfo_Type) :: GroupInfo
GroupInfo%nPrec = nprec_LibBase
IF(GroupInfo%ng .NE. ng_LibBase) THEN
  CALL TERMINATE('For Transient Calculation w/ XS Lib., the 47 Group Structure is allowed only')
ENDIF
END SUBROUTINE

SUBROUTINE InitChidLib(Core, FmInfo, GroupInfo, PE, nTracerCntl, Chid, ng)
USE TYPEDEF,          ONLY : CoreInfo_Type,   FmInfo_Type,    GroupInfo_Type,   PE_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE SUBGRP_MOD,       ONLY : FxrChiGen
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: Chid(ng)
INTEGER :: ng

! REAL :: chidg(8)
! INTEGER :: ig, gidx

IF(nTracerCntl%lchidgen) THEN 
  CALL GetAverageChid(Core, FmInfo, GroupInfo, PE, nTracerCntl, Chid, ng)

ELSE
  Chid(1:ng) = DelaySpectrum
END IF

! chidg = 0.
! DO ig = 1, ng
!   gidx = GroupInfo%INvGCStruct(ig)
!   chidg(gidx) = chidg(gidx) + chid(ig)
! END DO


! IF(PE%master) THEN 
!  PRINT '(10es14.6)', chidg
! END IF
END SUBROUTINE

SUBROUTINE InitChidklib(Core, FmInfo, GroupInfo, PE, nTracerCntl, chid, Chidk, ng, nprec)
USE TYPEDEF,          ONLY : CoreInfo_Type,   FmInfo_Type,    GroupInfo_Type,   PE_Type,    &
                             FxrInfo_Type,    Cell_Type,       Pin_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE MacXsLib_Mod,     ONLY : GetMacChidk
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: chid(ng)
REAL :: Chidk(ng, nprec)
INTEGER :: ng, nprec

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
REAL, POINTER :: phis(:,:,:)
REAL :: spectrum(ng), locchidk(ng, nprec)
REAL :: volsum, phisum(ng), vol, psisum
REAL :: chidksum(ng, nprec), locpsidsum
INTEGER :: myzb, myze, nxy
INTEGER :: nLocalFxr, nFsrInFxr, FxrIdxSt, FsrIdxSt
INTEGER :: iz, ipin, icel, ifxr, ifsr, ig, ireg, iprec
INTEGER :: i, j , ierr

IF(.NOT. nTracerCntl%lxslib) RETURN

Cell => Core%CellInfo
Pin => Core%Pin
phis => FmInfo%phis

myzb = PE%myzb
myze = PE%myze
ng = GroupInfo%ng
nxy = Core%nxy

chidksum = 0. 
DO iz = myzb, myze
  DO ipin = 1, nxy
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel = Pin(ipin)%Cell(iz)
    nLocalFxr = Cell(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      myFxr => FmInfo%Fxr(ifxr,iz)
      IF(.NOT. myFxr%lfuel) CYCLE
      nFsrInFxr = Cell(icel)%nFsrInFxr(j)
      spectrum = 0.
      phisum = 0.
      volsum = 0.
      locpsidsum = 0.
      DO i = 1, nFsrInFxr
        ireg = Cell(icel)%MapFxr2FsrIdx(i, j)
        ifsr = FsrIdxSt + ireg - 1
        vol = Cell(icel)%vol(ireg)
        volsum = volsum + vol
        DO ig = 1, ng
          phisum(ig) = phisum(ig) + vol * phis(ifsr, iz, ig)
        END DO 
        locpsidsum = locpsidsum + vol * FmInfo%psi(ifsr, iz)
      END DO
      DO ig = 1, ng
        spectrum(ig) = phisum(ig) / volsum
      END DO 
      CALL GetMacChidk(myFXR, Spectrum, locchidk, 1, ng, myFxr%niso, ng, nprec)
      DO iprec = 1, nprec
        DO ig = 1, ng
          chidksum(ig, iprec) = chidksum(ig, iprec) + locpsidsum * locchidk(ig, iprec) * myFxr%beta(iprec)          
        END DO
      END DO
    END DO
  END DO
END DO 

NULLIFY(Cell, Pin, phis, myFxr)

CALL MPI_ALLREDUCE(chidksum, chidk, ng*nprec, MPI_DOUBLE, MPI_SUM, PE%MPI_CMFD_COMM, ierr)

psisum = sum(chidk(:,:))
chid = 0.
DO iprec = 1, nprec
  DO ig = 1, ng
    chid(ig) = chid(ig) + chidk(ig, iprec)
  END DO
END DO
chid = chid / psisum
DO iprec = 1, nprec
  psisum = sum(chidk(1:ng, iprec))
  DO ig = 1, ng
    chidk(ig, iprec) = chidk(ig, iprec) / psisum
  END DO
END DO

!IF(PE%master) THEN 
!  DO iprec = 1, nprec
!    PRINT '(10es14.6)', chidk(:,iprec)
!  END DO 
!END IF


END SUBROUTINE

SUBROUTINE GetAverageChid(Core, FmInfo, GroupInfo, PE, nTracerCntl, avgchid, ng)
USE TYPEDEF,          ONLY : CoreInfo_Type,     FmInfo_Type,    GroupInfo_Type,   PE_Type,  &
                             FxrInfo_Type,      Pin_Type,       Cell_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE MacXsLib_Mod,     ONLY : GetMacChid
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: avgchid(ng)
INTEGER :: ng

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
REAL, POINTER :: phis(:,:,:)
REAL :: spectrum(ng), locchid(ng)
REAL :: volsum, phisum(ng), vol, psisum
REAL :: chidsum(ng), locpsidsum
INTEGER :: myzb, myze, nxy
INTEGER :: nLocalFxr, nFsrInFxr, FxrIdxSt, FsrIdxSt
INTEGER :: iz, ipin, icel, ifxr, ifsr, ig, ireg
INTEGER :: i, j , ierr

IF(.NOT. nTracerCntl%lxslib) RETURN

Cell => Core%CellInfo
Pin => Core%Pin
phis => FmInfo%phis

myzb = PE%myzb
myze = PE%myze
ng = GroupInfo%ng
nxy = Core%nxy

chidsum = 0. 
DO iz = myzb, myze
  DO ipin = 1, nxy
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel = Pin(ipin)%Cell(iz)
    nLocalFxr = Cell(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      myFxr => FmInfo%Fxr(ifxr,iz)
      IF(.NOT. myFxr%lfuel) CYCLE
      nFsrInFxr = Cell(icel)%nFsrInFxr(j)
      spectrum = 0.
      phisum = 0.
      volsum = 0.
      locpsidsum = 0.
      DO i = 1, nFsrInFxr
        ireg = Cell(icel)%MapFxr2FsrIdx(i, j)
        ifsr = FsrIdxSt + ireg - 1
        vol = Cell(icel)%vol(ireg)
        volsum = volsum + vol
        DO ig = 1, ng
          phisum(ig) = phisum(ig) + vol * phis(ifsr, iz, ig)
        END DO 
        locpsidsum = locpsidsum + vol * FmInfo%psi(ifsr, iz)
      END DO
      DO ig = 1, ng
        spectrum(ig) = phisum(ig) / volsum
      END DO 
      CALL GetMacChid(myFXR, Spectrum, locchid, 1, ng, myFxr%niso, ng)
      DO ig = 1, ng
        chidsum(ig) = chidsum(ig) + locpsidsum * locchid(ig) * myFxr%betat
      END DO 
    END DO
  END DO
END DO 

NULLIFY(Cell, Pin, phis, myFxr)

CALL MPI_ALLREDUCE(chidsum, avgchid, ng, MPI_DOUBLE, MPI_SUM, PE%MPI_CMFD_COMM, ierr)

psisum = sum(avgchid(1:ng))
DO ig = 1, ng
  avgchid(ig) = avgchid(ig) / psisum
END DO

END SUBROUTINE

SUBROUTINE GetAverageChi(Core, FmInfo, PE, avgchi, ng)
USE TYPEDEF,          ONLY : CoreInfo_Type,     FmInfo_Type,    PE_Type,    Pin_Type,  &
                             Cell_Type
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_Type) :: PE
REAL :: avgchi(ng)
INTEGER ::  ng

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL :: chisum(ng), locpsisum, psisum
INTEGER :: FxrIdxSt, FsrIdxSt
INTEGER :: myzb, myze, nxy, nLocalFxr, nFsrInFxr
INTEGER :: ipin, iz, icel, ifxr, ifsr, ig, ireg
INTEGER :: i, j, ierr

Pin => Core%Pin
Cell => Core%CellInfo
myzb = PE%myzb
myze = PE%myze
nxy = Core%nxy

chisum = 0. 
DO iz = myzb, myze
  DO ipin = 1, nxy
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel = Pin(ipin)%Cell(iz)
    nLocalFxr = Cell(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      IF(.NOT. FmInfo%Fxr(ifxr, iz)%lfuel) CYCLE
      nFsrInFxr = Cell(icel)%nFsrInFxr(j)
      locpsisum = 0.
      DO i = 1, nFsrInFxr
        ireg = Cell(icel)%MapFxr2FsrIdx(i,j)
        ifsr = FsrIdxSt + ireg - 1
        locpsisum = locpsisum + FmInfo%psi(ifsr, iz) * Cell(icel)%vol(ireg)
      END DO 
      DO ig = 1, ng
        chisum(ig) = chisum(ig) + locpsisum * FmInfo%Fxr(ifxr, iz)%Chi(ig)
      END DO 
    END DO 
  END DO 
END DO 
CALL MPI_ALLREDUCE(chisum, avgchi, ng, MPI_DOUBLE, MPI_SUM, PE%MPI_CMFD_COMM, ierr)

psisum = sum(avgchi(1:ng))
DO ig = 1, ng
  avgchi(ig) = avgchi(ig) / psisum
END DO
WRITE(408, '(10es14.6)') avgchi

END SUBROUTINE

SUBROUTINE GetAverageBeta(Core, FmInfo, PE, avgbeta, nprec)
USE TYPEDEF,          ONLY : CoreInfo_Type,     FmInfo_Type,    PE_Type,    Pin_Type,  &
                             Cell_Type
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_Type) :: PE
REAL :: avgbeta
INTEGER :: nprec

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL :: betasum, locpsisum, psisum, betat
REAL :: Buf0(2), Buf(2)
INTEGER :: FxrIdxSt, FsrIdxSt
INTEGER :: myzb, myze, nxy, nLocalFxr, nFsrInFxr
INTEGER :: ipin, iz, icel, ifxr, ifsr, ig, ireg
INTEGER :: i, j, ierr

Pin => Core%Pin
Cell => Core%CellInfo
myzb = PE%myzb
myze = PE%myze
nxy = Core%nxy

betasum = 0. 
psisum = 0.
DO iz = myzb, myze
  DO ipin = 1, nxy
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel = Pin(ipin)%Cell(iz)
    nLocalFxr = Cell(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      nFsrInFxr = Cell(icel)%nFsrInFxr(j)
      locpsisum = 0.
      DO i = 1, nFsrInFxr
        ireg = Cell(icel)%MapFxr2FsrIdx(i,j)
        ifsr = FsrIdxSt + ireg - 1
        locpsisum = locpsisum + FmInfo%psi(ifsr, iz) * Cell(icel)%vol(ireg)
      END DO 
      betat = sum(FmInfo%Fxr(ifxr, iz)%beta(1:nprec))
      psisum = psisum + locpsisum
      betasum = betasum + locpsisum * betat
    END DO 
  END DO 
END DO 

Buf0 = (/betasum, psisum/)
CALL MPI_ALLREDUCE(Buf0, Buf, 2, MPI_DOUBLE, MPI_SUM, PE%MPI_CMFD_COMM, ierr)
betasum = Buf(1) 
psisum = Buf(2)

avgbeta = betasum / psisum

END SUBROUTINE

SUBROUTINE InitLambdaLib(Lambda, refnid, llibdcy_del, lfitbeta)
REAL :: Lambda(nPrec_LibBase)
INTEGER :: refnid
LOGICAL :: llibdcy_del, lfitbeta

REAL :: A(3,4), y(3,6), fmult, temp
REAL :: tmpr(4)
INTEGER :: id, i, j, jl, ju
INTEGER :: ir, ic, k
INTEGER :: r1, r2, r3


IF(llibdcy_del) THEN
  id = MapNucl(refnid)
  Lambda = ldiso(id)%dcy_del
  IF(lfitbeta) THEN
    DO i = 1, nelthel
      IF(ldiso(i)%ifis .NE. 1) CYCLE
      IF(ldiso(i)%beta(1) .LE. 0.) CYCLE
      IF(i .EQ. id) CYCLE
      !PRINT*, 'NID: ', ldiso(i)%nid
      !print'(6ES14.6)', Lambda
      !print'(6ES14.6)', ldiso(i)%dcy_del
      !PRINT'(a20, 7ES14.6)', 'Before fitting: ', ldiso(i)%beta(0:6)
      DO j = 1, 6
        IF(j .EQ. 1) THEN
          jl = 6
        ELSE
          jl = j - 1
        END IF
        IF(j .EQ. 6) THEN
          ju = 1
        ELSE
          ju = j + 1
        END IF
        r1 = 2
        IF(lambda(jl) .GT. 1) THEN 
          r2 = 1
          r3 = 3
        ELSE
          r3 = 1
          r2 = 3
        END IF
        
        A(r1,1) = 1.
        A(r2,1) = Lambda(jl)
        A(r3,1) = 1./Lambda(jl)
        A(r1,2) = 1.
        A(r2,2) = Lambda(j)
        A(r3,2) = 1./Lambda(j)
        A(r1,3) = 1.
        A(r2,3) = Lambda(ju)
        A(r3,3) = 1./Lambda(ju)
        !- preserving total delayed neutron fraction
        A(r1,4) = ldiso(i)%beta(j)
        !- preserving average decay constant
        A(r2,4) = ldiso(i)%beta(j) * ldiso(i)%dcy_del(j)
        !A(r2,4) = ldiso(i)%beta(j) * Lambda(j)
        !- preserving mean decay time
        A(r3,4) = ldiso(i)%beta(j) / ldiso(i)%dcy_del(j)
        !A(r3,4) = ldiso(i)%beta(j) / Lambda(j)
        !-Gaussian elimination
        DO k = 1, 2
          DO ir = k+1, 3
            fmult = A(ir,k) / A(k,k)
            DO ic = 1, 4
              A(ir,ic) = A(ir, ic) - fmult*A(k,ic)
            END DO 
          END DO
          IF(k .EQ. 1) THEN
            IF(abs(A(2,2)) .LE. abs(A(3,2))) THEN
              tmpr(:) = A(3,:)
              A(3,:) = A(2,:)
              A(2,:) = tmpr(:)
            END IF
          END IF
        END DO 
        !-Backward substitution
        y(3,j) = A(3,4) / A(3,3)
        DO k = 2, 1, -1
          temp = 0.
          DO ic = k+1,3
            temp = temp + A(k,ic)*y(ic,j)
          END DO 
          y(k,j) = (A(k,4)-temp) / A(k,k)
        END DO 
      END DO
      !-Generate new beta
      DO j = 1, 6
        IF(j .EQ. 1) THEN
          jl = 6
        ELSE
          jl = j - 1
        END IF
        IF(j .EQ. 6) THEN
          ju = 1
        ELSE
          ju = j + 1
        END IF
        ldiso(i)%beta(j) = y(3,jl) + y(2,j) + y(1,ju)
      END DO
      !DO k = 1, 3
      !  print'(6ES14.6)', y(k,:)
      !END DO 
      ldiso(i)%beta(0) = sum(ldiso(i)%beta(1:6))
      !PRINT'(a20, 7ES14.6)', 'After fitting: ', ldiso(i)%beta(0:6)
    END DO
  END IF
ELSE
  Lambda = DecayConst
END IF

END SUBROUTINE

SUBROUTINE FxrBeta(XsMac, Fxr, FxrPhi, ng, iResBeg, iResEnd)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrINfo_Type) :: Fxr
REAL :: FxrPhi(ng)
INTEGER :: ng, iResBeg, iResEnd
!EffMacXs and FxrAvgPhi should be called prior to this routine
INTEGER, POINTER :: IdIso(:)
REAL, POINTER :: pnum(:)

REAL :: beta(nprec_LibBase)

REAL :: FisIso                !Fission Rate of One Isotope
REAL :: FisSum

INTEGER :: niso
INTEGER :: ig, iso, id, iprec, ifis

niso = Fxr%niso
IdIso => Fxr%IdIso; pnum => Fxr%pnum
FisSum = 0; Beta = 0
DO iso = 1, niso
  id = IdIso(iso);  ifis = MapFis(id)
  IF(ifis .EQ. 0) CYCLE
  id = MapNucl(idiso(iso)); 
  FisIso = 0;
  DO ig = 1, ng
    FisIso = FisIso + XsMac%IsoXsMacNf(iso, ig) * FxrPhi(ig)
    IF(Fxr%lres .AND. ig .GE. iResBeg .AND. ig .LE. iResEnd) THEN
      FisIso = FisIso + XsMac%IsoXsMacNf(iso, ig) * Fxr%fresoFIso(iso,ig) * FxrPhi(ig)
    ELSE
      FisIso = FisIso + XsMac%IsoXsMacNf(iso, ig) * FxrPhi(ig)
    END IF

  ENDDO
  FisSum = FisSum + FisIso
  
  DO iprec = 1, nprec_LibBase
    Beta(iprec) = Beta(iprec) + FisIso * ldiso(id)%Beta(iprec)
  ENDDO
ENDDO
IF(FisSum .GT. 0) THEN
  Fissum = 1._8 / FisSum
  DO iprec = 1, nprec_LibBase
    Beta(iprec) = Beta(iprec) * FisSum 
  ENDDO
ENDIF
Fxr%Beta = Beta
END SUBROUTINE

SUBROUTINE FxrVelo(Fxr, Temp, ng, lfixvel)
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
REAL :: Temp
INTEGER :: ng
LOGICAL :: lfixvel

REAL, POINTER :: pnum(:)
REAL, POINTER :: invvelsum(:)
REAL :: sigab10, velot, ubar, ebar, ndsum
REAL :: wt1, wt2
INTEGER :: it1, it2
INTEGER :: ig, id
INTEGER :: iso, niso

IF(lfixvel) THEN
  id = MapNucl(5010)
  CALL  XsTempInterpolation(id, ldiso(id), temp, wt1, wt2, it1, it2)

  DO ig = 1, ng
    IF(enbhel(ig) .LT. 0.1e+6) THEN  !less than 0.1 MeV
      sigab10 = (wt2 * ldiso(id)%siga(ig, it2) + wt1 * ldiso(id)%siga(ig, it1))
      velot = 3837._8 * 2.2e5_8/sigab10
    ELSE
      ubar=0.5_8*(uhel(ig)+uhel(ig+1))
      ebar=1.0e7_8/exp(ubar)
      velot=2.2e5_8*sqrt(ebar/0.0253_8)
    ENDIF
    Fxr%Velo(ig) = Velot
  ENDDO
ELSE
  niso = Fxr%niso
  pnum => Fxr%pnum
  ALLOCATE(invvelsum(ng))
  invvelsum = 0.
  ndsum = 0.
  DO iso = 1, niso
    id = MapNucl(Fxr%idiso(iso))
    CALL XsTempInterpolation(id, ldiso(id), temp, wt1, wt2, it1, it2)
    DO ig = 1, ng
      invvelsum(ig) = invvelsum(ig) + pnum(iso) * (wt1*ldiso(id)%invV(ig,it1) + wt2*ldiso(id)%invV(ig,it2))
    END DO
    ndsum = ndsum + pnum(iso)
  END DO
  DO ig = 1, ng
    Fxr%Velo(ig) = ndsum / invvelsum(ig) * 100
  END DO 
  NULLIFY(pnum)
  DEALLOCATE(invvelsum)
  !print*, 'njoyvel', it1, it2, wt1, wt2, temp
  !print'(47es14.6)', ldiso(id)%invvel(:, it1)
END IF
!write(211,'(47ES14.6)') Fxr%Velo
!STOP

END SUBROUTINE


END MODULE