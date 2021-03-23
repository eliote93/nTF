#include <defines.h>
SUBROUTINE OutputEdit
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,    PowerDist_Type,       DancoffDist_Type
USE GEOM,              ONLY : Core,             ng
USE Core_mod,          ONLY : FmInfo,           CmInfo,                     eigv,        &
                              ThInfo,           GroupInfo
USE TH_mod,            ONLY : THVar
USE CNTL,              ONLY : nTracerCntl
USE FILES,             ONLY : IO8,              io15
USE PE_mod,            ONLY : PE
USE ioutil,            ONLY : message,          PrintReal1DarrayTo2Darray
USE VTK_Mod,           ONLY : ProcessVTK
USE BinOutp_Mod,       ONLY : BinaryOutput
#ifdef MPI_ENV
USE MpiComm_mod, ONLY : BCAST, MPI_SYNC
#endif
IMPLICIT NONE
TYPE(PowerDist_Type) :: PowerDist
TYPE(DancoffDist_Type) :: DancoffDist
INTEGER :: io
LOGICAL :: Master

Master = TRUE
#ifdef MPI_ENV
Master = PE%Cmfdmaster
#endif

IF(Master) THEN
  WRITE(mesg, *) 'Writing Calculation Results in Output File...'
  CALL message(io8, TRUE, TRUE, mesg)
  WRITE(mesg, '(14x, a, I5, 4x a, I5)') 'Case No :', nTracerCntl%CaseNo, 'Calculation No :', nTracerCntl%CalNoId
  CALL message(io8, FALSE, TRUE, mesg)
ENDIF


IF(Master) THEN
  io = io8
  WRITE(io, *)
  WRITE(io, '(a)') '========================================================================'
  WRITE(io, '(5x, a)') 'Results'
  WRITE(io, '(5x, a, I5)') 'Case No        :', nTracerCntl%CaseNo
  WRITE(io, '(5x, a, I5)') 'Calculation No :', nTracerCntl%CalNoId
  WRITE(io, '(a)') '========================================================================'

  WRITE(io, *)
  WRITE(io, '(7x, A7, 3x, F8.5)') 'k-eff =', eigv
  WRITE(io, *)
ENDIF
!
CALL CorePowerCal(Core, CmInfo, PowerDIst, ng, nTracerCntl, PE)
IF (nTracerCntl%lED .AND. nTracerCntl%lrestrmt) CALL TreatDancoffData(Core,DancoffDist, PE)
#ifdef MPI_ENV
CALL BCAST(PowerDist%Fm3DNormalizer, PE%MPI_CMFD_COMM)
CALL BCAST(PowerDist%Pin3DNormalizer, PE%MPI_CMFD_COMM)
#endif
!
IF(Master) THEN
  CALL PrintRadialPower(io, Core, PowerDist)
!
  CALL PrintLocalPinPower(io, Core, PowerDist)
!
  CALL Asy3DPower(io, Core, PowerDist)
!
  CALL AxialAvgAsy2DPower(io, Core, PowerDist)

  CALL Axial1DPower(io, Core, PowerDist)
  IF (nTracerCntl%lED) CALL PrintDancoffFactor(io, Core, DancoffDist)
ENDIF

CALL AxAvg2DFlux(io, core, CmInfo, PowerDist, ng, PE)

CALL Asy3DFlux(io, core, CmInfo, PowerDist, ng, PE)

CALL RadAvg1DFlux(io, core, CmInfo, PowerDist, ng, PE)

CALL CoreSpectrum(io, core, CmInfo, PowerDist, ng, PE)

CALL Pin3DFlux(io, core, CmInfo, PowerDist, ng, nTracerCntl, PE)

CALL PrintDepletion(io, Core, FmInfo, CmInfo, nTracerCntl, PE)

IF(nTracerCntl%lfeedback .AND. MASTER) THEN
  CALL PrintFuelAvgTemp(io, Core, THInfo, THVar)

  CALL PrintFuelCentTemp(io, Core, THInfo, THVar)

  CALL PrintFuelSurfTemp(io, Core, THInfo, THVar)

  CALL PrintCoolantTemp(io, Core, THInfo, THVar)
  !CALL PrintCoolantNDAvgTemp(Core, FmInfo)
  !CALL PrintOthersTemp(Core, FmInfo)
ENDIF

!CALL PrintPinXS(io, Core, CmInfo, ng, PE)
IF(Master) WRITE(io, '(a)') '========================================================================'


IF(nTRACERCntl%OutpCntl%lBoutp) CALL BinaryOutput()

IF(nTRACERCntl%OutpCntl%VisMod .GT. 0) CALL ProcessVTK(Core, FMInfo, CmInfo, PowerDist, ThInfo, GroupInfo, nTracerCntl, PE)

IF(nTracerCntl%OutpCntl%lPinXsOut) CALL PrintPinXs(io15, Core, CmInfo, GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%lCspGenOut) CALL CspGenOut(io, Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)

IF(nTracerCntl%OutpCntl%nRegXsOut .GT. 0) CALL ProcessEffXs(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%nRegMGXsOut .GT. 0) CALL ProcessEffMGXs(Core, FmInfo, CmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%nRegMAT0Out .GT. 0) CALL ProcessEffMAT0(Core, FmInfo, GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%nRegPhimOut .GT. 0) CALL ProcessPhims(Core, FmInfo, GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%lSSPHout) CALL ProcessSSPHfactor(Core,nTracerCntl)

CALL FreePowerDist(PowerDist, PE)
END SUBROUTINE

SUBROUTINE PrintRadialPower(io, Core, PowerDist)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,             PowerDist_Type
USE ioutil,  ONLY : PrintReal1DarrayTo2Darray
USE BasicOperation, ONLY : CP_VA
IMPLICIT NONE
TYPE(PowerDist_Type) :: PowerDist
TYPE(CoreInfo_Type) :: Core
INTEGER, INTENT(IN) :: io

REAL, POINTER :: Dat(:)

INTEGER:: iasy, ixa, iya, itype
INTEGER :: nxya, nxy, nx, ny, ndat
INTEGER :: i, j, k, l

CHARACTER(120) :: formatdata

nDat = Core%nxyc
nxya = Core%nxya
ALLOCATE(DAT(nDat))
iasy = 0
301 FORMAT('Assembly', I5, x,'at' ,x, '(', I7, I7,' )')
WRITE(formatdata,'(a)') '(200F9.5)'
DO i = 1, nxya
  ixa = Core%Asy(i)%ixa; iya = Core%Asy(i)%iya
  itype = Core%Asy(i)%AsyType
  IF(.NOT. Core%AsyInfo(itype)%lFuel) CYCLE
  nx = Core%asyInfo(itype)%nx; ny = Core%asyInfo(itype)%ny
  nxy = Core%ASyInfo(itype)%nxy

  CALL EditOutputDat(Core, iType, PowerDist%PinPower2D(1:nxy, i), DAT(1:nDat), nxy, nDat)

  nxy = Core%nxyc0; nx = Core%nxc0; ny = Core%nyc0
  iasy = iasy + 1
  WRITE(io, 301) iasy, ixa, iya
  CALL PrintReal1DArrayTo2Darray(io, DAT(1:nxy), nxy, nx, ny, formatdata)
ENDDO
WRITE(io, *)
DEALLOCATE(DAT)
END SUBROUTINE

SUBROUTINE PrintDancoffFactor(io, Core, DancoffDist)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, DancoffDist_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
USE BasicOperation, ONLY : CP_VA
IMPLICIT NONE
INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(DancoffDist_Type) :: DancoffDist

REAL, POINTER :: DAT(:)

INTEGER :: nxy, nx, ny, nxya, nz
INTEGER ::  ndat
INTEGER :: iasy, ixa, iya, itype, iz
INTEGER :: i, j, k, l
CHARACTER(120) :: formatdata

nxya = Core%nxya; nz = Core%nz
nDat = Core%nxyc

WRITE(io, '(A)') '- Local Pin Dancoff -'
301 FORMAT('Assembly', I5, x,'at' ,x, '(', I7, I7,' )')
WRITE(formatdata,'(a)') '(200F9.5)'

ALLOCATE(DAT(nDat))

DO iz = 1, nz
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  WRITE(io, '(x, A5, I5)'), 'Plane', iz
  iasy = 0
  DO i = 1, nxya
    ixa = Core%Asy(i)%ixa; iya = Core%Asy(i)%iya
    itype = Core%Asy(i)%AsyType
    IF(.NOT. Core%AsyInfo(itype)%lFuel) CYCLE
    nx = Core%asyInfo(itype)%nx; ny = Core%asyInfo(itype)%ny
    nxy = Core%ASyInfo(itype)%nxy

    CALL EditOutputDat(Core, iType, DancoffDist%Dancoff3D(1:nxy, i, iz), DAT(1:nDat), nxy, nDat)

    nxy = Core%nxyc0; nx = Core%nxc0; ny = Core%nyc0
    iasy = iasy + 1
    WRITE(io, 301) iasy, ixa, iya
    CALL PrintReal1DarrayTo2Darray(io, DAT(1:nxy), nxy, nx, ny, formatdata)
  ENDDO
  WRITE(io, *)
ENDDO
WRITE(io, *)
DEALLOCATE(DAT)
END SUBROUTINE

SUBROUTINE PrintLocalPinPower(io, Core, PowerDist)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
USE BasicOperation, ONLY : CP_VA
IMPLICIT NONE
INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(PowerDist_Type) :: PowerDist

REAL, POINTER :: DAT(:)

INTEGER :: nxy, nx, ny, nxya, nz
INTEGER ::  ndat
INTEGER :: iasy, ixa, iya, itype, iz
INTEGER :: i, j, k, l
CHARACTER(120) :: formatdata

nxya = Core%nxya; nz = Core%nz
nDat = Core%nxyc

WRITE(io, '(A)') '- Local Pin Power -'
301 FORMAT('Assembly', I5, x,'at' ,x, '(', I7, I7,' )')
WRITE(formatdata,'(a)') '(200F9.5)'

ALLOCATE(DAT(nDat))

DO iz = 1, nz
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  WRITE(io, '(x, A5, I5)'), 'Plane', iz
  iasy = 0
  DO i = 1, nxya
    ixa = Core%Asy(i)%ixa; iya = Core%Asy(i)%iya
    itype = Core%Asy(i)%AsyType
    IF(.NOT. Core%AsyInfo(itype)%lFuel) CYCLE
    nx = Core%asyInfo(itype)%nx; ny = Core%asyInfo(itype)%ny
    nxy = Core%ASyInfo(itype)%nxy

    CALL EditOutputDat(Core, iType, PowerDist%PinPower3D(1:nxy, i, iz), DAT(1:nDat), nxy, nDat)

    nxy = Core%nxyc0; nx = Core%nxc0; ny = Core%nyc0
    iasy = iasy + 1
    WRITE(io, 301) iasy, ixa, iya
    CALL PrintReal1DarrayTo2Darray(io, DAT(1:nxy), nxy, nx, ny, formatdata)
  ENDDO
  WRITE(io, *)
ENDDO
WRITE(io, *)
DEALLOCATE(DAT)
END SUBROUTINE

SUBROUTINE Avg2DPower(io, Core, PowerDist)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
IMPLICIT NONE
INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(PowerDist_Type) :: PowerDist
REAL :: maxv
INTEGER :: nxy, nx, ny
INTEGER :: i, j
CHARACTER(120) :: formatdata

nxy = Core%nxya
nx = Core%nxa; ny = Core%nya
WRITE(io, '(A)') '- Axially Averaged 2D Assembly Power -'
WRITE(io, *)
WRITE(formatdata,'(a)') '(8x, 200F8.4)'
CALL PrintReal1DarrayTo2Darray(io,  PowerDist%AsyPower2D(1:nxy), nxy, nx, ny, formatdata)

maxv = maxval(PowerDist%AsyPower2D(1:nxy))
WRITE(io, *)
WRITE(io, '(5x, A7, 2x, F7.3)') 'Max. = ', maxv
WRITE(io, *)
END SUBROUTINE

SUBROUTINE AxialAvgAsy2DPower(io, Core, PowerDist)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(PowerDist_Type) :: PowerDist
REAL :: maxv
REAL, POINTER :: Dat(:)
INTEGER :: nxy, nx, ny, nxy0
INTEGER :: i, j,  ix, iy
CHARACTER(120) :: formatdata

nxy = Core%nxya
nx = Core%nxa; ny = Core%nya
nxy0 = nx * ny
ALLOCATE(Dat(nxy0))
CALL CP_CA(Dat, 0._8, nxy0)
DO i = 1, nxy
  ix = Core%Asy(i)%ixa; iy = Core%Asy(i)%iya
  j = nx * (iy - 1) + ix
  Dat(j) = PowerDist%AsyPower2D(i)
ENDDO
WRITE(io, '(A)') '- Axially Averaged 2D Assembly Power -'
WRITE(io, *)
WRITE(formatdata,'(a)') '(8x, 200F8.4)'
CALL PrintReal1DarrayTo2Darray(io,  Dat(1:nxy0), nxy0, nx, ny, formatdata)
maxv = maxval(PowerDist%AsyPower2D(1:nxy))
WRITE(io, *)
WRITE(io, '(5x, A7, 2x, F7.3)') 'Max. = ', maxv
WRITE(io, *)
DEALLOCATE(Dat)
END SUBROUTINE

SUBROUTINE AxialAvgA2DAsyPower(io, Core, PowerDist)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
IMPLICIT NONE
INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(PowerDist_Type) :: PowerDist
REAL :: maxv
INTEGER :: nxy, nx, ny
INTEGER :: i, j
CHARACTER(120) :: formatdata

nxy = Core%nxya
nx = Core%nxa; ny = Core%nya
WRITE(io, '(A)') '- Axially Averaged Assembly 2D Power -'
WRITE(formatdata,'(a)') '(8x, 200F8.4)'
CALL PrintReal1DarrayTo2Darray(io,  PowerDist%AsyPower2D(1:nxy), nxy, nx, ny, formatdata)

maxv = maxval(PowerDist%AsyPower2D(1:nxy))
WRITE(io, *)
WRITE(io, '(5x, A7, 2x, F7.3)') 'Max. = ', maxv
WRITE(io, *)
END SUBROUTINE


SUBROUTINE Asy3DPower(io, Core, PowerDist)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(PowerDist_Type) :: PowerDist
REAL :: maxv
REAL, POINTER :: Dat(:)
INTEGER :: nxy, nx, ny, nz, nxy0
INTEGER :: i, j, iz, ixa, iya
CHARACTER(120) :: formatdata

nxy = Core%nxya; nz = Core%nz
nx = Core%nxa; ny = Core%nya
nxy0 = nx * ny

ALLOCATE(Dat(nxy0))

WRITE(io, '(A)') ' - Assembly 3D Power -'
WRITE(formatdata,'(a)') '(8x, 200F8.4)'
DO iz = 1, nz
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  WRITE(io, '(5x, A5, I5)'), 'Plane', iz
  CALL CP_CA(Dat, 0._8, nxy0)
  DO i = 1, nxy
    ixa = Core%Asy(i)%ixa; iya = Core%Asy(i)%iya
    j = nx * (iya - 1) + ixa
    Dat(j) =  PowerDist%AsyPower3D(iz, i)
  ENDDO
  CALL PrintReal1DarrayTo2Darray(io, Dat(1:nxy0), nxy0, nx, ny, formatdata)
ENDDO
WRITE(io, *)
END SUBROUTINE

SUBROUTINE Axial1DPower(io, Core, PowerDist)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
IMPLICIT NONE
INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(PowerDist_Type) :: PowerDist
INTEGER :: nz, iz
CHARACTER(120) :: formatdata

nz = Core%nz
WRITE(io, '(A)') ' - Radially Averaged 1D Power -'
WRITE(formatdata,'(a)') '(8x, 200F8.4)'
WRITE(io, '(2x, A8, 3x, A5)') 'Plane #', 'Power'
DO iz = 1, nz
  WRITE(io, '(I10, F8.4)') iz, PowerDist%Axial1DPower(iz)
ENDDO
WRITE(io, '(A)')
END SUBROUTINE

SUBROUTINE AxAvg2DFlux(io, core, CmInfo, PowerDist, ng, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,   PowerDist_Type,  CMInfo_Type,     PE_TYPE,   &
                    AsyInfo_Type,    Asy_Type,        PinXs_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray

#ifdef MPI_ENV
USE BasicOperation, ONLY : CP_CA, CP_VA
USE MpiComm_mod, ONLY : Reduce
#endif
IMPLICIT NONE
INTEGER :: io, ng
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(PE_TYPE) :: PE

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(PinXS_Type), POINTER :: PinXs(:, :)
REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :), AsyVol(:, :)
REAL, POINTER :: hz(:)
INTEGER :: nz, nxy, nxya, nxa, nya, myzb, myze
INTEGER :: nAsyType, AsyType
INTEGER :: iz, ig, iasy, ixy, ixa, iya, i, j, k
LOGICAL :: master, slave
REAL,POINTER :: Avg2DFlx(:, :,  :), row(:)
#ifdef MPI_ENV
REAL, POINTER :: Buf(:, :)
#endif

Master = PE%Cmfdmaster; Slave = .NOT. Master

PhiC => CmInfo%PhiC; PinXs => CmInfo%PinXs
AsyInfo => Core%AsyInfo; Asy => Core%Asy
hz => Core%hz;
PinVol => Core%PinVol; AsyVol => Core%AsyVol
nz = Core%nz; myzb = PE%myzb; myze =PE%myze
nAsyType = Core%nAsyType; nxya = Core%nxya;
nxy = Core%nxy; nxa = Core%nxa; nya = Core%nya

ALLOCATE(Avg2DFlx(nxa,nya, ng))
#ifdef MPI_ENV
ALLOCATE(Buf(nxa, nya))
#endif
DO ig = 1, ng
  Avg2DFlx(:, :, ig) = 0
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
    nxy = AsyInfo(AsyType)%nxy
    Do iz = myzb, myze  !Axial Sweep
      DO i = 1, nxy  !Cell Sweep
        ixy = Asy(iasy)%GlobalPinIdx(i)
        Avg2DFlx(ixa, iya, ig) = Avg2DFlx(ixa, iya, ig) + Pinvol(ixy, iz) * Phic(ixy, iz, ig)
      ENDDO
    ENDDO
  ENDDO !
#ifdef MPI_ENV
  CALL CP_CA(Buf, zero, nxa, nya)
  CALL REDUCE(Avg2DFlx(:, :, ig), Buf(:, :), nxa, nya, PE%MPI_CMFD_COMM, .FALSE.)
  IF(Master) CALL CP_VA(Avg2DFlx(:, :, ig), Buf, nxa, nya)
#endif
ENDDO

#ifdef MPI_ENV
DEALLOCATE(Buf)
#endif

IF(Master) THEN
  WRITE(io, '(A)') ' -  Axially Averaged 2D Flux -'
  WRITE(io, *)
  DO iya = 1, nya
    DO ig = 1, ng
      WRITE(io, '(8x, 200(1pe15.4))') (Avg2DFlx(ixa, iya, ig), ixa = 1, nxa)
    ENDDO
    WRITE(io, *)
  ENDDO
ENDIF
DEALLOCATE(Avg2DFlx)

NULLIFY(PhiC); NULLIFY(PinXs)
NULLIFY(AsyInfo); NULLIFY(Asy)
NULLIFY(hz); NULLIFY(PinVol)
NULLIFY(AsyVol)
END SUBROUTINE

SUBROUTINE RadAvg1DFlux(io, core, CmInfo, PowerDist, ng, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,   PowerDist_Type,  CMInfo_Type,     PE_TYPE,   &
                    AsyInfo_Type,    Asy_Type,        PinXs_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
#ifdef MPI_ENV
USE MpiComm_mod, ONLY : SENDRECV
#endif
IMPLICIT NONE
INTEGER :: io, ng
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(PE_TYPE) :: PE

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(PinXS_Type), POINTER :: PinXs(:, :)
REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :), AsyVol(:, :)
REAL, POINTER :: hz(:)
INTEGER :: nz, nxy, nxya, nxa, nya, myzb, myze
INTEGER :: nAsyType, AsyType
INTEGER :: iz, ig, iasy, ixy, ixa, iya, i, j, k
LOGICAL :: master, slave
REAL,POINTER :: Avg1DFlx(:, :)
#ifdef MPI_ENV
INTEGER :: nproc, comm
INTEGER :: iz1, iz2
#endif

master = PE%CmfdMaster; slave = .NOT. master

PhiC => CmInfo%PhiC; PinXs => CmInfo%PinXs
AsyInfo => Core%AsyInfo; Asy => Core%Asy
hz => Core%hz;
PinVol => Core%PinVol; AsyVol => Core%AsyVol
nz = Core%nz
myzb = PE%myzb; myze = PE%myze
nAsyType = Core%nAsyType; nxya = Core%nxya;
nxy = Core%nxy; nxa = Core%nxa; nya = Core%nya

ALLOCATE(Avg1DFlx(ng, nz))

Avg1DFlx = 0
DO ig = 1, ng
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    nxy = AsyInfo(AsyType)%nxy
    Do iz = myzb, myze   !Axial Sweep
      DO i = 1, nxy  !Cell Sweep
        ixy = Asy(iasy)%GlobalPinIdx(i)
        Avg1DFlx(ig, iz) = Avg1DFlx(ig, iz) + Pinvol(ixy, iz) * Phic(ixy, iz, ig)
      ENDDO
    ENDDO
  ENDDO !
ENDDO

#ifdef MPI_ENV
nproc = PE%nCmfdProc; Comm = PE%MPI_CMFD_COMM
IF(Master) THEN
  DO i = 1, nproc - 1
    iz1 = PE%AxDomRange(1, i); iz2 = PE%AxDomRange(2, i)
    CALL SendRecv(Avg1DFlx(:, iz1:iz2), ng, iz2 - iz1 + 1, i, FALSE, COMM)
  ENDDO
ELSE
  CALL SendRecv(Avg1DFlx(:, myzb:myze), ng, myze - myzb + 1, 0, TRUE, COMM)
ENDIF
#endif

IF(master) THEN
  WRITE(io, '(A)') ' -  Radially Averaged 1D Flux -'
  WRITE(io, *)
  WRITE(io, '(A)') '                --- > Plane #'
  WRITE(io, '(2x,A,200(I10,5X))') 'Energy Group',(i, i=1,nz)
  DO ig = 1, ng
    WRITE(io, '(I10,4x, 200(1pe15.4))') ig, (Avg1DFlx(ig, iz), iz = 1, nz)
  ENDDO
  WRITE(io, *)
ENDIF

NULLIFY(PhiC); NULLIFY(PinXs)
NULLIFY(AsyInfo); NULLIFY(Asy)
NULLIFY(hz); NULLIFY(PinVol)
NULLIFY(AsyVol)

END SUBROUTINE

SUBROUTINE Asy3DFlux(io, core, CmInfo, PowerDist, ng, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,   PowerDist_Type,  CMInfo_Type,     PE_TYPE,   &
                    AsyInfo_Type,    Asy_Type,        PinXs_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
#ifdef MPI_ENV
USE MpiComm_mod, ONLY : SENDRECV
#endif
IMPLICIT NONE
INTEGER :: io, ng
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(PE_TYPE) :: PE

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(PinXS_Type), POINTER :: PinXs(:, :)
REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :), AsyVol(:, :)
REAL, POINTER :: hz(:)
INTEGER :: nz, nxy, nxya, nxa, nya, myzb, myze
INTEGER :: nAsyType, AsyType
INTEGER :: iz, ig, iasy, ixy, ixa, iya, i, j, k
LOGICAL :: master, slave
REAL,POINTER :: Asy3DFlx(:,:,:,:), row(:)
#ifdef MPI_ENV
INTEGER :: COMM, nproc
INTEGER :: iz1, iz2
#endif

Master = PE%CmfdMaster; Slave = .NOT. Master
PhiC => CmInfo%PhiC; PinXs => CmInfo%PinXs
AsyInfo => Core%AsyInfo; Asy => Core%Asy
hz => Core%hz;
PinVol => Core%PinVol; AsyVol => Core%AsyVol
nz = Core%nz
nAsyType = Core%nAsyType; nxya = Core%nxya;
nxy = Core%nxy; nxa = Core%nxa; nya = Core%nya
myzb = PE%myzb; myze = PE%myze

#ifndef MPI_ENV
ALLOCATE(Asy3DFlx(nxa,nya, ng, nz))
#else
IF(Master) THEN
  ALLOCATE(Asy3DFlx(nxa, nya, ng, nz))
ELSE
  ALLOCATE(Asy3DFlx(nxa, nya, ng, myzb:myze))
ENDIF
#endif
Do iz = myzb, myze   !Axial Sweep
  DO ig = 1, ng
    Asy3DFlx(:, :, ig, iz) = 0
    DO iasy = 1, nxya
      AsyType = Asy(iasy)%AsyType
      ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
      nxy = AsyInfo(AsyType)%nxy
      DO i = 1, nxy  !Cell Sweep
        ixy = Asy(iasy)%GlobalPinIdx(i)
        Asy3DFlx(ixa, iya, ig, iz) = Asy3DFlx(ixa, iya, ig, iz) + Pinvol(ixy, iz) * Phic(ixy, iz, ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO

#ifdef MPI_ENV
nproc = PE%nCmfdProc
COMM = PE%MPI_CMFD_COMM
IF(Master) THEN
  DO i = 1, nproc - 1
    iz1 = PE%AxDomRange(1, i); iz2 = PE%AxDomRange(2, i)
    CALL SendRecv(Asy3DFlx(:, :, :, iz1:iz2), nxa, nya, ng, iz2 - iz1 + 1,       &
                  i, FALSE, COMM)
  ENDDO
ELSE
  CALL SendRecv(Asy3DFlx(:, :, :, myzb:myze), nxa, nya, ng, myze - myzb + 1,       &
                0, TRUE, COMM)
  !ipartner, lSend, comm
ENDIF
#endif

IF(Master) THEN
  WRITE(io, '(A)') ' -  Planewise Radial Flux  -'
  WRITE(io, *)
  DO iz = 1, nz
    WRITE(io, '(5x, A5, I5)'), 'Plane', iz
    DO iya = 1, nya
      DO ig = 1, ng
        WRITE(io, '(8x, 200(1pe15.4))') (Asy3DFlx(ixa, iya, ig, iz), ixa = 1, nxa)
      ENDDO
      WRITE(io, *)
    ENDDO
  ENDDO
ENDIF
DEALLOCATE(Asy3DFlx)

NULLIFY(PhiC); NULLIFY(PinXs)
NULLIFY(AsyInfo); NULLIFY(Asy)
NULLIFY(hz); NULLIFY(PinVol)
NULLIFY(AsyVol)

END SUBROUTINE

SUBROUTINE Pin3DFlux(io, core, CmInfo, PowerDist, ng, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,   PowerDist_Type,    CmInfo_Type,  &
                                PE_Type,                                          &
                                Asy_Type,        AsyInfo_Type
USE CNTL,                ONLY : nTracerCntl_Type

USE BasicOperation,      ONLY : CP_CA,           CP_VA
#ifdef MPI_ENV
USE MPIComm_Mod,         ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: io, ng

TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
REAL, POINTER :: PhiC(:, :, :)

INTEGER, PARAMETER :: noutline = 10

INTEGER :: i, j, k, ibeg, iend
INTEGER :: ixy, ixya, ixa, iya, ix, iy, iz, ig
INTEGER :: iasytype

INTEGER :: nz, myzb, myze
INTEGER :: npinout, nstep, nout
INTEGER :: idx(4, noutline)

CHARACTER(5) :: GrpIdx
REAL, POINTER :: OutDat(:, :, :), buf(:, :, :)


Asy => Core%Asy; AsyInfo => Core%AsyInfo
PhiC => CmInfo%PhiC

nz = Core%nz; myzb = PE%myzb; myze = PE%myze
npinout = nTracerCntl%OutpCntl%FluxOutList(1, 0)

ALLOCATE(OutDat(ng, nz, noutline))
ALLOCATE(Buf(ng, nz, noutline))

nstep = INT(npinout / noutline)
IF(mod(npinout, noutline) .NE. 0) nstep = nstep + 1
IF(PE%MASTER) THEN
  WRITE(io, '(A)') ' -  Flux Distribution for selected Pins -  '
  WRITE(io, '(7x, A)') 'Pin Identifier'
  DO i = 1, npinout
    WRITE(io, '(7x, I5, 3x, A2, 3x, (A1, 4I4, A3))') i, '=>', '[', nTracerCntl%OutpCntl%FluxOutList(1:4, i), ']'
  ENDDO
  WRITE(io, *)
ENDIF
DO j = 1, nstep
  ibeg = (j - 1) * noutline + 1; iend = j * noutline
  iend = min(iend, npinout)
  nout = iend - ibeg +1
  CALL CP_CA(OutDat, 0._8, ng, nz, noutline)
  DO i = ibeg, iend
    k  = i - ibeg + 1
    ixa = nTracerCntl%OutpCntl%FluxOutList(1, i); iya = nTracerCntl%OutpCntl%FluxOutList(2, i)
    ix = nTracerCntl%OutpCntl%FluxOutList(3, i);  iy = nTracerCntl%OutpCntl%FluxOutList(4, i)
    ixya = Core%CoreIdx(ixa, iya); iasytype = Core%CoreMap(ixya)
    ixy = AsyInfo(iasytype)%Pin2DIdx(ix, iy)
    ixy = Asy(ixya)%GlobalPinIdx(ixy)
    idx(:, k) = nTracerCntl%OutpCntl%FluxOutList(2, i)
    DO ig = 1, ng
      DO iz = myzb, myze
        OutDat(ig, iz, k)= PhiC(ixy, iz, ig)
      ENDDO
    ENDDO
  ENDDO
#ifdef MPI_ENV
  CALL REDUCE(OutDat, Buf, ng, nz, noutline, PE%MPI_CMFD_COMM, .FALSE.)
  CALL CP_VA(OutDat, buf, ng, nz, noutline)
#endif
  IF(.NOT. PE%MASTER) CYCLE
  WRITE(io, '(7x, A)') 'Pin Index ->'
  WRITE(io, '(4x, A5, 2x, A5, 2x, 100(I10,2x))') 'GROUP', 'PLANE', (i, i = ibeg, iend)
  DO ig = 1, ng
    WRITE(GrpIdx, '(I5)') ig
    DO iz = 1, nz
      IF(iz .NE. 1) GrpIdx=''
      WRITE(io, '(4x, A5, 2x, I5, 2x, 100(e10.3,2x))') GrpIdx, iz, (OutDat(ig, iz, k), k = 1, nout)
    ENDDO
  ENDDO
  WRITE(io, *)
ENDDO


NULLIFY(Asy, ASyInfo)

END SUBROUTINE

SUBROUTINE CoreSpectrum(io, core, CmInfo, PowerDist, ng, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,   PowerDist_Type,  CMInfo_Type,   &
                    PE_TYPE
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
#ifdef MPI_ENV
USE MpiComm_mod, ONLY : REDUCE
#endif
IMPLICIT NONE
INTEGER :: io, ng
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(PE_TYPE) :: PE

REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :)
REAL, POINTER :: hz(:)
INTEGER :: nz, nxy, nxya, nxa, nya, myzb, myze
INTEGER :: nAsyType, AsyType
INTEGER :: iz, ig, iasy, ixy, ixa, iya, i, j, k
LOGICAL :: master, slave
REAL :: Spectrum(ng)
#ifdef MPI_ENV
REAL :: buf(ng)
#endif

Master = PE%Master; Slave = .NOT. Master
myzb = PE%myzb; myze = PE%myze
PhiC => CmInfo%PhiC;
hz => Core%hz;
PinVol => Core%PinVol
nxy = Core%nxy; nz = Core%nz

DO ig = 1, ng
  Spectrum(ig) = 0
  DO iz = myzb, myze
    DO ixy = 1, nxy
      Spectrum(ig) = Spectrum(ig) + Phic(ixy, iz, ig) * PinVol(ixy, iz)
    ENDDO
  ENDDO
ENDDO

#ifdef MPI_ENV
CALL REDUCE(Spectrum, buf, ng, PE%MPI_CMFD_COMM, FALSE)
Spectrum=buf
#endif

IF(Master) THEN
  WRITE(io, '(A)') ' - Core Spectrum -'
  WRITE(io, *)
  DO ig = 1, ng
    WRITE(io, '(5x,I5,1pe15.4)') ig, Spectrum(ig)
  ENDDO
  WRITE(io, *)
ENDIF
NULLIFY(PhiC)
NULLIFY(hz)
NULLIFY(PinVol)

END SUBROUTINE

SUBROUTINE FreePowerDist(PowerDist, PE)
USE PARAM
USE TypeDef,     ONLY : PowerDist_Type,      PE_TYPE
IMPLICIT NONE
TYPE(PowerDist_Type) :: PowerDist
TYPE(PE_TYPE) :: PE

DEALLOCATE(PowerDist%PinPower3D)
IF(PE%CMFDmaster) THEN
  DEALLOCATE(PowerDist%PinPower2D)
  DEALLOCATE(PowerDist%AsyPower2D)
  DEALLOCATE(PowerDist%AsyPower3D)
  DEALLOCATE(PowerDist%Axial1DPower)
ENDIF
END SUBROUTINE

SUBROUTINE FillCentDat(Dat, nxy0, nx0, ny0 , nxy, nx, ny)
USE BASICOPERATION,  ONLY : CP_CA
IMPLICIT NONE
REAL :: DAT(nxy)
INTEGER :: nx0, ny0, nxy0, nxy, nx, ny
REAL :: Buf(nx, ny)

INTEGER :: i, j, k
INTEGER :: ixy, ix, iy, ixy0, ix0, iy0

CALL CP_CA(Buf(1:nx, 1:ny), 0._8, nx, ny)
DO iy0 = 1, ny0
  DO ix0 = 1, nx0
    ixy0 = nx0 * (iy0 -1) + ix0
    ix = nx0 - ix0 + 1; iy = ny0 - iy0 + 1
    ix = nx - ix + 1; iy = ny - iy + 1;
    Buf(ix, iy) = DAT(ixy0)
  ENDDO
ENDDO
CALL CP_CA(Dat(1:nxy), 0._8, nxy)
ixy = 0
DO iy = 1, ny
  DO ix = 1, nx
    ixy = ixy + 1
    DAT(ixy) = Buf(ix, iy)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE EditOutputDat(Core, iAsyType, DAT0, DAT1, nDat0, nDat1)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type,     AsyInfo_Type,     PinInfo_type,         &
                       Cell_type
USE BasicOperation, ONLY : CP_VA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
INTEGER :: iAsyType
REAL :: Dat0(nDat0), Dat1(nDat1)
INTEGER:: nDat0, nDat1
INTEGER :: nxy, nx, ny, nxy0, nx0, ny0

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: ix, iy, ixy, ipin, icel

LOGICAL :: lCent
AsyInfo => Core%AsyInfo; PinInfo => Core%PinInfo
CellInfo => Core%CellInfo
nxy0 = AsyInfo(iAsyType)%nxy
IF(Core%lGap) THEN
  nxy = 0
  DO ixy = 1, nxy0
    ipin = AsyInfo(iAsyType)%Pin(ixy)
    icel = PinInfo(iPin)%Cell(1)
    IF(CellInfo(icel)%lGap) CYCLE
    nxy = nxy + 1
    Dat1(nxy) = Dat0(ixy)
  ENDDO
ELSE
  nxy = nxy0
  CALL CP_VA(Dat1(1:nxy), Dat0(1:nxy), nxy)
ENDIF
nx = AsyInfo(iAsyType)%nx; ny = AsyInfo(iAsyType)%ny
IF(Core%lGap) THEN
  IF(AsyInfo(iAsyType)%lCentXY) THEN
    nx= nx - 1; ny = ny - 1
  ELSEIF(AsyInfo(iAsyType)%lCentX) THEN
    nx = nx -2; ny = ny - 1
  ELSEIF(AsyInfo(iAsyType)%lCentY) THEN
    nx = nx -1; ny = ny -2
  ELSE
    nx = nx - 2; ny = ny - 2
  ENDIF
ENDIF

lCent = Core%AsyInfo(iAsyType)%lCentX
lCent = lCent .OR. Core%AsyInfo(iAsyType)%lCentY
lCent = lCent .OR. Core%AsyInfo(iAsyType)%lCentXY
IF(lCent) THEN
  nx0 = nx; ny0 = ny; nxy0 = nxy
  nxy = Core%nxyc0; nx = Core%nxc0; ny = Core%nyc0
  CALL FillCentDat(Dat1(1:nxy), nxy0, nx0, ny0 , nxy, nx, ny)
ENDIF
END SUBROUTINE







SUBROUTINE PrintFXRXs_CEA(Core, THInfo, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
USE param,   only : CKELVIN, PI
USE TYPEDEF,      ONLY : CoreInfo_Type,  THInfo_type, CmInfo_Type, FmInfo_type,    &
                         PE_TYPE,        GroupInfo_Type, PinXS_TYPE, &
                         AsyInfo_type, Asy_Type, Pin_Type, Cell_Type, FxrInfo_type , XsMac_type, PE_TYPE
USE CNTL,         ONLY : nTracerCntl_Type,           OutpCntl_Type
USE FILES,        ONLY : LOCALFN, caseid
USE IOUTIL,       ONLY : OpenFile
USE Depl_Mod,  ONLY : DeplCntl
USE BenchXs,        ONLY : XsBaseBen, xssm1ben, xssm2ben, xssm3ben
USE XsUtil_mod,   ONLY : AllocXsMac
USE TH_Mod,           ONLY : GetPinFuelTemp
IMPLICIT NONE
!
TYPE(CoreInfo_Type) :: Core
TYPE(THInfo_Type) :: THInfo
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(XsMac_Type) :: XsMac
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
REAL :: TempRef
INTEGER :: io

REAL, POINTER :: MacXsSm(:,:)  !IsoMacXsSm(isosize,ng,ng)
REAL, POINTER :: MacXsSmP1(:,:) !IsoMacXsSmP1(isosize,ng,ng)
!
CHARACTER(256) :: PinXsFn
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!
INTEGER :: nxy, nz, myzb, myze, ng, nxya, nlocalFxr, nIsoInFXR, nFSRinFXR, nCoreFXR, nchi, nx, ny
INTEGER, POINTER :: nPlaneFXR(:), lmix(:)
INTEGER :: ixy, ixyl, iz, ig, igg, igb, ige, ifxr, ifsr, icel, iasytype, iasy, fxridxst, fsridxst, ix, iy
INTEGER :: i, j, k, itype, idx, sidx, iscat
REAL, POINTER :: sigt(:), siga(:), sigf(:), signf(:), chi(:), sigs(:,:)
INTEGER, POINTER :: prof(:), prof0(:), prof1(:), prof2(:), prof3(:)
INTEGER :: nscat, scatod, gb, ge, fxrscatod
INTEGER :: imix
LOGICAL :: lfuel, lres
REAL :: dum, scale
!
CHARACTER(256) :: fn, matname, mnz, mnza, mnzap, mnzapf !materialname z-asy-pin-fxr
LOGICAL :: lxslib

INTERFACE

SUBROUTINE GenFxrMacXS(XsMac, MacXsSm, MacXsSmP1, Fxr, ifxr, Tempref, ig1, ig2, Core, ipin, iz, GroupInfo, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type, XsMac_Type,  FxrInfo_Type,  GroupInfo_Type,  PE_Type
USE CNTL,          ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type),pointer :: Fxr(:,:)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: ig1, ig2, ipin, iz, ifxr
REAL :: TempRef
REAL :: MacXsSm(GroupInfo%ng,GroupInfo%ng)
REAL :: MacXsSmP1(GroupInfo%ng,GroupInfo%ng)
END SUBROUTINE

END INTERFACE
write(*,*) '                  entering CEA_XS ... BYS edit 16/06/09'
io=43
lxslib=ntracerCntl%lxslib
AsyInfo => Core%AsyInfo; Asy => Core%Asy
nCorefxr= core%nCoreFxr
Pin => Core%Pin; Cell => Core%CellInfo
nxy = Core%nxy; nz = Core%nz; ng = GroupInfo%ng; nchi= GroupInfo%nChi
nxya= Core%nxya ! n-assemblies
ScatOd=nTracerCntl%ScatOd

!write(*,*) 'ncoreFXR', ncoreFXR

ALLOCATE(MacXsSm(ng,ng))
ALLOCATE(MacXsSmP1(ng,ng))
ALLOCATE(sigt(ng), siga(ng), sigf(ng), signf(ng), chi(ng), sigs(4,ng*ng) )
ALLOCATE(prof(2*ng), prof0(2*ng),prof1(2*ng),prof2(2*ng),prof3(2*ng))
ALLOCATE(nPlaneFXR(nz),lmix(100))

nCoreFXR=0;
myzb = PE%myzb; myze = PE%myze
DO iz = myzb, myze
    lmix=0;
    nPlaneFXR(iz)=0
    DO iasy = 1, nxya
        iasytype=Asy(iasy)%AsyType
        nxy=AsyInfo(iasytype)%nxy
        DO ixyl = 1, nxy
            ixy = Asy(iasy)%GlobalPinIdx(ixyl)
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            nCoreFXR=nCoreFXR+nlocalFXR
            !nPlaneFXR(iz)=nPlaneFXR(iz)+nlocalFXR
            DO j = 1, nLocalFxr
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                imix=myFXR%imix
                lfuel=myFXR%lfuel
                lres=myFXR%lres
                IF( lmix(imix).NE.0 .AND. .NOT.lRes)THEN
                    CONTINUE
                ELSE
                    lmix(imix)=lmix(imix)+1
                    nPlaneFXR(iz)=nPlaneFXR(iz)+1
                ENDIF
            ENDDO

        ENDDO
    ENDDO
ENDDO


idx=0
DO iz = myzb, myze
    lmix=0;
    WRITE(fn,'(A)') TRIM(caseid)
    IF( DeplCntl%nburnupstep .GT. 1 )THEN
        WRITE(fn,'(A,A)') TRIM(fn),'_d'
        IF( DeplCntl%nowstep .LT. 10 )THEN
            WRITE(fn,'(A,A,i1)') TRIM(fn),'0',DeplCntl%nowstep
        ELSE
            WRITE(fn,'(A,i2)') TRIM(fn),DeplCntl%nowstep
        ENDIF
    ENDIF
    WRITE(fn,'(A,A)') TRIM(fn),'_z'
    IF( iz .LT. 10 )THEN
        WRITE(fn,'(A,A,i1)') TRIM(fn),'0',iz
    ELSE
        WRITE(fn,'(A,i2)') TRIM(fn),iz
    ENDIF

    WRITE(fn,'(A,A)') TRIM(fn),'_CEA.xsl'
    OPEN(unit=io, file = fn, status = 'replace')

WRITE(io,'(a,i4,i8)') 'MacroscopicCrossSections[', ng, nPlaneFXR(iz)
    DO iasy = 1, nxya
        iasytype=Asy(iasy)%AsyType
        nxy=AsyInfo(iasytype)%nxy
        nx=AsyInfo(iasytype)%nx;ny=AsyInfo(iasytype)%ny
        DO iy = 1, ny
        DO ix = 1, nx
            ixyl = AsyInfo(iasytype)%Pin2dIdx(ix,iy)
            ixy = Asy(iasy)%GlobalPinIdx(ixyl)
            !ix = Pin(ixy)%ix
            !iy = Pin(ixy)%iy
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            IF( lxslib )    Tempref = GetPinFuelTemp(Core, FmInfo%FXR, iz, ixy)
            DO j = 1, nLocalFxr
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                imix=myFXR%imix
                lfuel=myFXR%lfuel
                lres=myFXR%lres
                IF( lmix(imix).NE.0 .AND. .NOT.lRes)THEN
                    CONTINUE
                ELSE
                    lmix(imix)=lmix(imix)+1

                nFsrInFxr = myFxr%nFsrInFxr
                nisoInFxr = myFxr%niso
                IF( lxslib )THEN
                    CALL GenFxrMacXS(XsMac, MacXsSm, MacXsSmP1, FmInfo%FXR, ifxr, Tempref, 1, ng, Core, ixy, iz, GroupInfo, nTRACERCntl, PE)
                ELSE
                    IF(.NOT. XsMac%lalloc) THEN
                        XsMac%ng = ng
                        CALL AllocXsMac(XsMac)
                    ENDIF
                    itype=myFxr%imix
                    CALL xsbaseBen(itype, 1, ng, 1, ng, nTracerCntl%lscat1, XsMac)
                ENDIF
                idx=idx+1
                prof0=0;prof1=0;prof2=0;prof3=0;

                !--- Transport correction rollback : tr > total  / r544 edit
                IF( nTracerCntl%scatod .GT. 0 )THEN
                    DO ig = 1, ng
                        !XsMac%xsmactr(ig)=XsMac%xsmaca(ig) + sum(XsMac%XsMacSm(ig,1:ng))
                    ENDDO
                    DO ig = 1, ng
                        !XsMac%XsMacSm(ig,ig)=XsMac%XsMacSm(ig,ig) + ( XsMac%XsMact(ig)-XsMac%XsMactr(ig) )
                    ENDDO
                    DO ig = 1, ng
                        XsMac%xsmact(ig)=XsMac%xsmaca(ig) + sum(XsMac%XsMacSm(ig,1:ng))
                    ENDDO
                    DO ig = 1, ng
                    sigt(ig)=XsMac%xsmact(ig)
                    ENDDO
                ELSE
                    DO ig = 1, ng
                    sigt(ig)=XsMac%xsmactr(ig)
                    ENDDO
                ENDIF


                DO ig = 1, ng
                    siga(ig)=XsMac%xsmaca(ig)
                    sigf(ig)=XsMac%xsmacf(ig)
                    signf(ig)=XsMac%xsmacnf(ig)
                    !IF( ig .LE. nchi .AND. XsMac%lfuel )THEN
                    !    lfuel=.TRUE.
                    IF( ig .LE. nchi .AND. lfuel )THEN
                        IF( lxslib )THEN
                            chi(ig)=myFxr%chi(ig)
                        ELSE
                            chi(ig)=XsMac%chi(ig)
                        ENDIF
                    ELSE
                        chi(ig)=0._8
                    ENDIF
                    DO igg = 1, ng
                        IF( XsMac%XsmacSm(ig,igg) .NE. 0 )THEN
                            IF( prof0(2*igg-1).EQ.0 )THEN
                                prof0(2*igg-1)=ig
                            ENDIF
                            prof0(2*igg)=ig
                        ENDIF
                    ENDDO
                    IF(ScatOd .GE. 1)THEN
                    DO igg = 1, ng
                        IF( XsMac%XsmacP1Sm(ig,igg) .NE. 0 )THEN
                            IF( prof1(2*igg-1).EQ.0 )THEN
                                prof1(2*igg-1)=ig
                            ENDIF
                            prof1(2*igg)=ig
                        ENDIF
                    ENDDO
                        IF(ScatOd .GE. 2)THEN
                        DO igg = 1, ng
                            IF( XsMac%XsmacP2Sm(ig,igg) .NE. 0 )THEN
                                IF( prof2(2*igg-1).EQ.0 )THEN
                                    prof2(2*igg-1)=ig
                                ENDIF
                                prof2(2*igg)=ig
                            ENDIF
                        ENDDO
                        IF(ScatOd .GE. 3)THEN
                            DO igg = 1, ng
                                IF( XsMac%XsmacP3Sm(ig,igg) .NE. 0 )THEN
                                    IF( prof3(2*igg-1).EQ.0 )THEN
                                        prof3(2*igg-1)=ig
                                    ENDIF
                                    prof3(2*igg)=ig
                                ENDIF
                            ENDDO
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO
                FXRScatOd=0;prof=prof0
                IF(sum(prof1).NE.0)THEN
                    FXRScatOd=1
                    DO ig=1,ng
                        prof(2*ig-1)=min(prof(2*ig-1),prof1(2*ig-1))
                        prof(2*ig)=max(prof(2*ig),prof1(2*ig))
                    ENDDO
                    IF(sum(prof2).NE.0)THEN
                        FXRScatOd=2
                        DO ig=1,ng
                            prof(2*ig-1)=min(prof(2*ig-1),prof2(2*ig-1))
                            prof(2*ig)=max(prof(2*ig),prof2(2*ig))
                        ENDDO
                        IF(sum(prof3).NE.0)THEN
                            FXRScatOd=3
                            DO ig=1,ng
                                prof(2*ig-1)=min(prof(2*ig-1),prof3(2*ig-1))
                                prof(2*ig)=max(prof(2*ig),prof3(2*ig))
                            ENDDO
                        ENDIF
                    ENDIF
                ENDIF
                sidx=0
                iscat=0;
                scale=(2._8*iscat+1._8)
                DO ig = 1, ng ! to
                    DO igg=prof(2*ig-1),prof(2*ig) !from
                        sidx=sidx+1
                        sigs(1,sidx)=XsMac%xsmacsm(igg,ig)*scale
                    ENDDO
                ENDDO
                IF(FXRScatOd.GE.1)THEN
                    sidx=0
                    iscat=1;
                    scale=(2._8*iscat+1._8)
                    DO ig = 1, ng ! to
                        DO igg=prof(2*ig-1),prof(2*ig) !from
                            sidx=sidx+1
                            sigs(2,sidx)=XsMac%xsmacp1sm(igg,ig)*scale
                        ENDDO
                    ENDDO
                    IF(FXRScatOd.GE.2)THEN
                        sidx=0
                        iscat=2;
                        scale=(2._8*iscat+1._8)
                        DO ig = 1, ng ! to
                            DO igg=prof(2*ig-1),prof(2*ig) !from
                                sidx=sidx+1
                                sigs(3,sidx)=XsMac%xsmacp2sm(igg,ig)*scale
                            ENDDO
                        ENDDO
                        IF(FXRScatOd.GE.3)THEN
                            sidx=0
                            iscat=3;
                            scale=(2._8*iscat+1._8)
                            DO ig = 1, ng ! to
                                DO igg=prof(2*ig-1),prof(2*ig) !from
                                    sidx=sidx+1
                                    sigs(4,sidx)=XsMac%xsmacp3sm(igg,ig)*scale
                                ENDDO
                            ENDDO
                        ENDIF
                    ENDIF
                ENDIF

                nscat=sidx

                IF( lres )THEN
                    WRITE(matname,'(a,i3,a,i3,a,i3,a,i3,a,i6,a,i2,a)') 'imix',imix, '_z', iz, '_asy', iasy, '_', ix, '_', iy, '_fxr', j, '_resT'
                ELSE
                    WRITE(matname,'(a,i3,a,i3,a,i3,a,i3,a,i6,a,i2,a)') 'imix',imix, '_z', iz, '_asy', iasy, '_', ix, '_', iy, '_fxr', j, '_resF'
                ENDIF
                CALL StripSpaces(matname)
                WRITE(io,'(a)') matname
                IF( lfuel )THEN
                    WRITE(io,'(a)') 'FissileIsotopes   1'
                ELSE
                    WRITE(io,'(a)') 'FissileIsotopes   0'
                ENDIF
                WRITE(io,'(a)') 'Total'
                WRITE(io,'(t2,1p,5e14.6)') sigt(:)
                WRITE(io,'(a)') 'Absorption'
                WRITE(io,'(t2,1p,5e14.6)') siga(:)
                IF( lfuel )THEN
                    WRITE(io,'(a)') 'Fission'
                    WRITE(io,'(t2,1p,5e14.6)') sigf(:)
                    WRITE(io,'(a)') 'NuFission'
                    WRITE(io,'(t2,1p,5e14.6)') signf(:)
                    WRITE(io,'(a)') 'FissionSpectrum'
                    WRITE(io,'(t2,1p,5e14.6)') chi(:)
                ENDIF
                !IF(
                WRITE(io,'(a,i2)') 'Transfer   ',FXRScatOd
                WRITE(io,'(a)') 'Profile'
                WRITE(io,'(12i6)') prof(:)
                DO iscat = 1, FXRScatOd+1
                    WRITE(io,'(t2,1p,5e14.6)') sigs(iscat,1:nscat)
                ENDDO

                ENDIF
            ENDDO ! FXR in Asy
        ENDDO ! Pinx in Asy
        ENDDO ! Piny in Asy
    ENDDO ! Global Asy
    WRITE(io,'(a)') ']'
CLOSE(io)
ENDDO !Z-axial plane

END SUBROUTINE



    subroutine StripSpaces(string)
    character(len=*) :: string

    integer :: stringLen
    integer :: last, actual

    stringLen = len (string)
    last = 1
    actual = 1

    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do

    end subroutine





SUBROUTINE PrintFXRXs(Core, THInfo, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
USE param,   only : CKELVIN
USE TYPEDEF,      ONLY : CoreInfo_Type,  THInfo_type, CmInfo_Type, FmInfo_type,    &
                         PE_TYPE,        GroupInfo_Type, PinXS_TYPE, &
                         AsyInfo_type, Asy_Type, Pin_Type, Cell_Type, FxrInfo_type , XsMac_type
USE CNTL,         ONLY : nTracerCntl_Type,           OutpCntl_Type
USE FILES,        ONLY : LOCALFN, caseid
USE IOUTIL,       ONLY : OpenFile
USE Depl_Mod,  ONLY : DeplCntl
USE BenchXs,        ONLY : XsBaseBen, xssm1ben, xssm2ben, xssm3ben
USE XsUtil_mod,   ONLY : AllocXsMac
USE TH_Mod,           ONLY : GetPinFuelTemp
IMPLICIT NONE
!
TYPE(CoreInfo_Type) :: Core
TYPE(THInfo_Type) :: THInfo
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(XsMac_Type) :: XsMac
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
REAL :: TempRef
INTEGER :: io, io2

REAL, POINTER :: MacXsSm(:,:)  !IsoMacXsSm(isosize,ng,ng)
REAL, POINTER :: MacXsSmP1(:,:) !IsoMacXsSmP1(isosize,ng,ng)
!
CHARACTER(256) :: PinXsFn
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!
INTEGER :: nxy, nz, myzb, myze, ng, nxya, nlocalFxr, nIsoInFXR, nFSRinFXR, nchi
INTEGER :: ixy, ixyl, iz, ig, igg, igb, ige, ifxr, ifsr, icel, iasytype, iasy, fxridxst, fsridxst, fidx
INTEGER :: i, j, k, idx, itype, imix
REAL :: chi
LOGICAL :: lxslib
!
CHARACTER(256) :: fn, fn2, matname
INTERFACE

SUBROUTINE GenFxrMacXS(XsMac, MacXsSm, MacXsSmP1, Fxr, ifxr, Tempref, ig1, ig2, Core, ipin, iz, GroupInfo, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type, XsMac_Type,  FxrInfo_Type,  GroupInfo_Type,  PE_Type
USE CNTL,          ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type),pointer :: Fxr(:,:)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: ig1, ig2, ipin, iz, ifxr
REAL :: TempRef
REAL :: MacXsSm(GroupInfo%ng,GroupInfo%ng)
REAL :: MacXsSmP1(GroupInfo%ng,GroupInfo%ng)
END SUBROUTINE

END INTERFACE
io=45
io2=46

AsyInfo => Core%AsyInfo; Asy => Core%Asy
Pin => Core%Pin; Cell => Core%CellInfo
nz = Core%nz; ng = GroupInfo%ng
nxya= Core%nxya ! n-assemblies
nchi= GroupInfo%nChi
ALLOCATE(MacXsSm(ng,ng))
ALLOCATE(MacXsSmP1(ng,ng))


lxslib=nTracerCntl%lxslib
idx=0

DO iz = PE%myzb, PE%myze
WRITE(fn,'(A)') TRIM(caseid)
IF( DeplCntl%nburnupstep .GT. 1 )THEN
    WRITE(fn,'(A,A)') TRIM(fn),'_d'
    IF( DeplCntl%nowstep .LT. 10 )THEN
        WRITE(fn,'(A,A,i1)') TRIM(fn),'0',DeplCntl%nowstep
    ELSE
        WRITE(fn,'(A,i2)') TRIM(fn),DeplCntl%nowstep
    ENDIF
ENDIF
    WRITE(fn,'(A,A)') TRIM(fn),'_z'
    IF( iz .LT. 10 )THEN
        WRITE(fn,'(A,A,i1)') TRIM(fn),'0',iz
    ELSE
        WRITE(fn,'(A,i2)') TRIM(fn),iz
    ENDIF
!WRITE(fn,'(A,i)') TRIM(fn), fidx

WRITE(fn,'(A,A)') TRIM(fn),'_nTFXR.xsl'
OPEN(unit=io, file = fn, status = 'replace')

WRITE(fn2,'(A)') TRIM(caseid)
IF( DeplCntl%nburnupstep .GT. 1 )THEN
    WRITE(fn2,'(A,A)') TRIM(fn),'_d'
    IF( DeplCntl%nowstep .LT. 10 )THEN
        WRITE(fn2,'(A,A,i1)') TRIM(fn2),'0',DeplCntl%nowstep
    ELSE
        WRITE(fn2,'(A,i2)') TRIM(fn2),DeplCntl%nowstep
    ENDIF
ENDIF
    WRITE(fn2,'(A,A)') TRIM(fn),'_z'
    IF( iz .LT. 10 )THEN
        WRITE(fn2,'(A,A,i1)') TRIM(fn),'0',iz
    ELSE
        WRITE(fn2,'(A,i2)') TRIM(fn),iz
    ENDIF

    DO iasy = 1, nxya
        iasytype=Asy(iasy)%AsyType
        nxy=AsyInfo(iasytype)%nxy
        DO ixyl = 1, nxy
            ixy = Asy(iasy)%GlobalPinIdx(ixyl)
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            IF( lxslib )    Tempref = GetPinFuelTemp(Core, FmInfo%FXR, iz, ixy)
            DO j = 1, nLocalFxr
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                imix=myFXR%imix
                nFsrInFxr = myFxr%nFsrInFxr
                !FsrIdxSt = myFxr%FsrIdxSt
                nisoInFxr = myFxr%niso
                IF( lxslib )THEN
                    CALL GenFxrMacXS(XsMac, MacXsSm, MacXsSmP1, FmInfo%FXR, ifxr, Tempref, 1, ng, Core, ixy, iz, GroupInfo, nTRACERCntl, PE)
                ELSE
                    IF(.NOT. XsMac%lalloc) THEN
                        XsMac%ng = ng
                        CALL AllocXsMac(XsMac)
                    ENDIF
                    itype=myFxr%imix
                    CALL xsbaseBen(itype, 1, ng, 1, ng, nTracerCntl%lscat1, XsMac)
                ENDIF
                idx=idx+1
                WRITE(io, '(a,i10)') ' base_micro ', idx
                !--- Transport correction rollback : tr > total  / r544 edit
                IF( nTracerCntl%scatod .GT. 0 )THEN
                    DO ig = 1, ng
                        XsMac%xsmactr(ig)=XsMac%xsmaca(ig) + sum(XsMac%XsMacSm(ig,1:ng))
                    ENDDO
                    DO ig = 1, ng
                        !XsMac%XsMacSm(ig,ig)=XsMac%XsMacSm(ig,ig) + ( XsMac%XsMact(ig)-XsMac%XsMactr(ig) )
                    ENDDO
                    DO ig = 1, ng
                        XsMac%xsmact(ig)=XsMac%xsmaca(ig) + sum(XsMac%XsMacSm(ig,1:ng))
                    ENDDO
                ENDIF

                IF( lxslib )THEN
                    DO ig = 1, ng
                        IF( ig .LE. nchi .AND. XsMac%xsmacnf(ig) .NE. 0 )THEN
                            chi=myFxr%chi(ig)
                        ELSE
                            chi=0._8
                        ENDIF
                        IF( nTracerCntl%scatod .GT. 0 )THEN
                        WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmact(ig), XsMac%xsmaca(ig), XsMac%xsmacnf(ig),XsMac%xsmacnf(ig), chi
                        ELSE
                        WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmactr(ig), XsMac%xsmaca(ig), XsMac%xsmacnf(ig),XsMac%xsmacnf(ig), chi
                        ENDIF

                    ENDDO ! Group sweep
                    DO ig = 1, ng
                        WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmacsm(ig,:)
                        !WRITE(io,'(t2,1p,100e14.6)') MacXsSm(ig,:)
                    ENDDO ! Group sweep
                    IF( nTracerCntl%scatod .GE. 1 )THEN
                        DO ig = 1, ng
                            WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmacp1sm(ig,:)
                        ENDDO ! Group sweep
                        IF( nTracerCntl%scatod .GE. 2 )THEN
                            DO ig = 1, ng
                                WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmacp2sm(ig,:)
                            ENDDO ! Group sweep
                            IF( nTracerCntl%scatod .GE. 3 )THEN
                                DO ig = 1, ng
                                    WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmacp3sm(ig,:)
                                ENDDO ! Group sweep
                            ENDIF
                        ENDIF
                    ENDIF
                ELSE
                    DO ig = 1, ng
                        IF( nTracerCntl%scatod .GT. 0 )THEN
                        WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmact(ig), XsMac%xsmaca(ig), XsMac%xsmacnf(ig),XsMac%xsmackf(ig),  XsMac%chi(ig)
                        ELSE
                        WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmactr(ig), XsMac%xsmaca(ig), XsMac%xsmacnf(ig),XsMac%xsmackf(ig),  XsMac%chi(ig)
                        ENDIF
                    ENDDO ! Group sweep
                    DO ig = 1, ng
                        WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmacsm(ig,:)
                    ENDDO ! Group sweep
                    IF( nTracerCntl%scatod .GE. 1 )THEN
                        CALL xssm1ben(itype, 1, ng, 1, ng, XsMac%xsmacp1sm)
                        DO ig = 1, ng
                            WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmacp1sm(ig,:)
                        ENDDO ! Group sweep
                        IF( nTracerCntl%scatod .GE. 2 )THEN
                            CALL xssm2ben(itype, 1, ng, 1, ng, XsMac%xsmacp2sm)
                            DO ig = 1, ng
                                WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmacp2sm(ig,:)
                            ENDDO ! Group sweep
                            IF( nTracerCntl%scatod .GE. 3 )THEN
                                CALL xssm3ben(itype, 1, ng, 1, ng, XsMac%xsmacp3sm)
                                DO ig = 1, ng
                                    WRITE(io,'(t2,1p,100e14.6)') XsMac%xsmacp3sm(ig,:)
                                ENDDO ! Group sweep
                            ENDIF
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO ! FXR in Asy
        ENDDO ! Pin in Asy
    ENDDO ! Global Asy
ENDDO !Z-axial plane
close(io)
close(io2)
END SUBROUTINE

SUBROUTINE PrintRayFile()
USE RAYS,   ONLY : RayInfo,            AziAngle,         PolarAngle,     &
                   ModRay,             AsyRay,           CoreRay,        &
                   RotRay,             CellRayBase
USE GEOM,   ONLY : AsyPitch,           AsyGeom,          lEdge,          &
                   CellInfo,           Core,             Asy,            &
                   nCellType,          nCellType0
USE FILES,        ONLY : caseid
IMPLICIT NONE
INTEGER :: i,j,k, iazi
INTEGER :: io
INTEGER :: nazi, npol, nmod, nasy, ncore, nrot, nphisv, nbd
LOGICAL :: lhex = .FALSE.
REAL :: del
CHARACTER*256 :: fn
REAL, POINTER :: nsegang(:)
CONTINUE
nAzi=RayInfo%nAziAngle
nPol=RayInfo%nPolarAngle
nmod=RayInfo%nModRay
nasy=RayInfo%nAsyRay
nCore=RayInfo%nCoreRay
nRot=RayInfo%nRotRay
nphiSv=RayInfo%nPhiangSv
Del=RayInfo%del

!--- File Open
io=47
WRITE(fn,'(A)') TRIM(caseid)
WRITE(fn,'(A,A)') TRIM(fn),'.rsf' !Ray Structure File
OPEN(unit=io, file = fn, status = 'replace')

!--- ANGLE
WRITE(io, '(a)') '! RayInfo- Del, nAzi, nPol / nMod, nAsy, nCore, nRot, nPhiSv'
WRITE(io, '(f13.4,2i13)') del,nazi,npol
WRITE(io, '(a)') '! Azimuthal Angle- '
DO i = 1, nAzi
    WRITE(io, '(i13,100E13.6)') i, AziAngle(i)%ang, AziAngle(i)%sinv, AziAngle(i)%cosv, AziAngle(i)%tanv, Aziangle(i)%Del, AziAngle(i)%weight
ENDDO
!--- Modular Ray
WRITE(io, '(a)') '! Modular Ray- NextRayIdx I/O, IO surf I/O, AngleIdx'
DO i = 1, nMod
    WRITE(io, '(i13,100i13)') i, ModRay(i)%NextRayIdx, ModRay(i)%InOutSurf, ModRay(i)%iAziAngIdx
ENDDO
!--- Assembly Ray
WRITE(io, '(a)') '! Assembly Ray- '
DO i = 1, nAsy
    WRITE(io, '(i13,100i13)') i, AsyRay(i)%AziAngIdx, AsyRay(i)%ModRayIdx, AsyRay(i)%NextRayIdx, AsyRay(i)%InOutSurf, AsyRay(i)%nCellRay, AsyRay(i)%PartialAsyRayFlag, AsyRay(i)%PartialAsyRayIdx
    DO j = 1, AsyRay(i)%nCellRay
        WRITE(io, '(a,i13,100i13)') '             ', j, AsyRay(i)%PinRaySurf(:,j), AsyRay(i)%PinIdx(j), AsyRay(i)%PinRayIdx(j)
    ENDDO
ENDDO
!--- Core Ray
WRITE(io, '(a)') '! Core Ray- '
DO i = 1, nCore
    WRITE(io, '(i13,100i13)') i, CoreRay(i)%nRay, CoreRay(i)%iAng, CoreRay(i)%InOutSurf
    DO j = 1, CoreRay(i)%nRay
        WRITE(io, '(a,i13,100i13)') '             ', j, CoreRay(i)%AsyRayIdx(j), CoreRay(i)%AsyIdx(j), CoreRay(i)%IXA(j), CoreRay(i)%IYA(j)
    ENDDO
ENDDO
!--- Rotational Ray
WRITE(io, '(a)') '! Rotational Ray'
DO i = 1, nRot
    WRITE(io, '(i13,100i13)') i, RotRay(i)%nRay, RotRay(i)%nSeg
    DO j = 1, RotRay(i)%nRay
        WRITE(io, '(a,i13,100i13)') '             ', j, RotRay(i)%RayIdx(j), RotRay(i)%Dir(j)
    ENDDO
ENDDO
CONTINUE
!--- Base Cell Ray
WRITE(io, '(a)') '! Base Cell Ray'
DO i = 1, nCellType
    IF( i .EQ. CellRayBase(i)%CellType )THEN
        WRITE(io, '(i13,100i13)') i, CellRayBase(i)%CellType, CellRayBase(i)%nLine
        DO j = 1, nAzi
            WRITE(io, '(a,i13,100i13)') '             ', j, CellRayBase(i)%AziAngIdx(j), CellRayBase(i)%nLines(j)
        ENDDO
        DO j = 1, CellRayBase(i)%nline
            WRITE(io, '(a,i13,100i13)') '             ', j, CellRayBase(i)%CellRay(j)%IdxSt, CellRayBase(i)%CellRay(j)%IdxEnd, CellRayBase(i)%CellRay(j)%nSeg
            DO  k = 1, CellRayBase(i)%CellRay(j)%nSeg
                WRITE(io, '(2a,i13,i13,e13.6)') '             ','             ', k, CellRayBase(i)%CellRay(j)%LocalFsrIdx(k), CellRayBase(i)%CellRay(j)%LenSeg(k)
            ENDDO
        ENDDO
    ELSE
        WRITE(io, '(i13,100i13)') i, CellRayBase(i)%CellType, 0
    ENDIF
ENDDO
WRITE(io, '(a)') '.'


ALLOCATE(nsegang(nazi))
nsegang=0
DO i = 1, nAsy
    iazi=AsyRay(i)%aziangidx
    DO j= 1, AsyRay(i)%ncellray
        k=AsyRay(i)%PinRayIdx(j)
        nsegang(iazi)=nsegang(iazi)+CellRayBase(1)%Cellray(k)%nseg
    enddo
enddo



CLOSE(io)



!--- GEOM output
!--- File Open
io=48
WRITE(fn,'(A)') TRIM(caseid)
WRITE(fn,'(A,A)') TRIM(fn),'.gsf' !Ray Structure File
OPEN(unit=io, file = fn, status = 'replace')

WRITE(io, '(a)') '! Base Geometry Information'
WRITE(io, '(100i13)') core%nxy, core%nxya, core%nCoreFxr, core%nCoreFsr, core%npintype, core%nCellType, core%nAsyType, core%nz
!--- Boundary Condition
WRITE(io, '(a)') '! Boundary Info'
DO i = 1, 6
    j=0
    IF(core%radsym(i)) j=1
    WRITE(io,'(2i13)') j, core%radbc(i)
ENDDO
!--- Axial Mesh height
DO i = 1, core%nz
    WRITE(io,'(i13,f13.6)') i, core%hz(i)
ENDDO

!------ PROBLEM DEFINED VARIABLES
WRITE(io, '(a)') '! Pin '
DO i = 1, core%nxy
    WRITE(io,'(i13,100i13)') i, core%pin(i)%ix, core%pin(i)%iy, core%pin(i)%iasy, core%pin(i)%ipin, core%pin(i)%pintype, core%pin(i)%asytype, core%pin(i)%ncell, core%pin(i)%nFsrMax, core%pin(i)%nFxrMax
    WRITE(io,'(i13,100i13)') core%pin(i)%nBD, core%pin(i)%FsrIdxSt, core%pin(i)%FxrIdxSt, core%pin(i)%lfuel, core%pin(i)%lgt, core%pin(i)%lradref, core%pin(i)%lBaffle, core%pin(i)%lBarrel, core%pin(i)%lGd
    DO j = 1, core%pin(i)%nBD
        WRITE(io,'(a,3i13)') '             ', j, core%Pin(i)%neighidx(j), core%Pin(i)%neighsurfidx(j)
    ENDDO
    DO j = 1, core%nz
        WRITE(io,'(a,2i13,f13.6,i13)') '             ', j, core%pin(i)%cell(j), core%pinvol(i,j), core%pin(i)%lmox(j)
    ENDDO
ENDDO
!--- Assembly
nbd=4
IF( lHex ) nbd=6
WRITE(io, '(a)') '! Assembly '
DO i = 1, Core%nxya
    WRITE(io,'(i13,100i13)') i, Core%Asy(i)%ixa, Core%Asy(i)%iya, Core%Asy(i)%asytype, Core%Asy(i)%partialAsyFlag
    DO j = 1, nbd
        WRITE(io,'(a,2i13)') '             ', core%asy(i)%neighIdx(j)
    ENDDO
    DO j = 1, Core%AsyInfo(Core%Asy(i)%asytype)%nxy
        WRITE(io,'(a,2i13)') '             ', j, core%asy(i)%GlobalPinIdx(j)
    ENDDO
ENDDO

!------ PRE DEFINED VALUES ---
!--- CellInfo
WRITE(io, '(a)') '! Cell Info'
DO i = 1,  core%ncellType
    WRITE(io,'(i13,100i13)') i, Core%CellInfo(i)%nDivAzi, Core%CellInfo(i)%nFsr, Core%CellInfo(i)%nFXR,             &
                                Core%CellInfo(i)%lEmpty, Core%CellInfo(i)%luse, Core%CellInfo(i)%lrect,             &
                                Core%CellInfo(i)%lres, Core%CellInfo(i)%lfuel, Core%CellInfo(i)%lgd
    WRITE(io,'(i13,100i13)') Core%CellInfo(i)%lmox, Core%CellInfo(i)%lccell, Core%CellInfo(i)%lcentx,               &
                             Core%CellInfo(i)%lcenty, Core%CellInfo(i)%lcentxy, Core%CellInfo(i)%lgap
    DO j = 1, Core%CellInfo(i)%nFxr
        WRITE(io,'(a,3i13)') '             ', j, Core%CellInfo(i)%FxrIdxSt(j), Core%CellInfo(i)%nFsrInFxr(j)
    ENDDO
    DO j = 1, Core%CellInfo(i)%nFsr
        WRITE(io,'(1p,a,2i13,e13.6)') '             ', j, Core%CellInfo(i)%ireg(j), core%CellInfo(i)%vol(j)
    ENDDO
ENDDO

!--- Pin Info
WRITE(io, '(a)') '! Pin Info'
DO i = 1,  core%npintype
    WRITE(io,'(i13,100i13)') i, Core%PinInfo(i)%ncell, Core%PinInfo(i)%nFsrMax, Core%PinInfo(i)%nFXRmax,            &
                                Core%PinInfo(i)%lEmpty, Core%PinInfo(i)%luse, Core%PinInfo(i)%lfuel,                &
                                Core%PinInfo(i)%lcentx, Core%PinInfo(i)%lcenty, Core%PinInfo(i)%lcentxy,            &
                                Core%PinInfo(i)%lgap
    DO j = 1, core%PinInfo(i)%ncell
        WRITE(io,'(a,2i13)') '             ', j, Core%PinInfo(i)%Cell(j)
    ENDDO
ENDDO



WRITE(io, '(a)') '.'
CLOSE(io)

ENDSUBROUTINE



SUBROUTINE ReadRayFile()
USE RAYS,   ONLY : RayInfo,            AziAngle,         PolarAngle,     &
                   ModRay,             AsyRay,           CoreRay,        &
                   RotRay,             CellRayBase
USE GEOM,   ONLY : AsyPitch,           AsyGeom,          lEdge,          &
                   CellInfo,           Core,             Asy,            &
                   nCellType,          nCellType0
USE FILES,        ONLY : caseid
IMPLICIT NONE
INTEGER :: i,j,k
INTEGER :: io
INTEGER :: nazi, npol, nmod, nasy, ncore, nrot, nphisv, nbd
LOGICAL :: lhex = .FALSE.
REAL :: del
CHARACTER*256 :: fn
CONTINUE
nAzi=RayInfo%nAziAngle
nPol=RayInfo%nPolarAngle
nmod=RayInfo%nModRay
nasy=RayInfo%nAsyRay
nCore=RayInfo%nCoreRay
nRot=RayInfo%nRotRay
nphiSv=RayInfo%nPhiangSv
Del=RayInfo%del

!--- File Open
io=47
WRITE(fn,'(A)') TRIM(caseid)
WRITE(fn,'(A,A)') TRIM(fn),'.rsf' !Ray Structure File
OPEN(unit=io, file = fn, status='replace')

!--- ANGLE
WRITE(io, '(a)') '! RayInfo- Del, nAzi, nPol / nMod, nAsy, nCore, nRot, nPhiSv'
WRITE(io, '(f13.4,2i13)') del,nazi,npol
WRITE(io, '(a)') '! Azimuthal Angle- '
DO i = 1, nAzi
    WRITE(io, '(i13,100E13.6)') i, AziAngle(i)%ang, AziAngle(i)%sinv, AziAngle(i)%cosv, AziAngle(i)%tanv, Aziangle(i)%Del, AziAngle(i)%weight
ENDDO
!--- Modular Ray
WRITE(io, '(a)') '! Modular Ray- NextRayIdx I/O, IO surf I/O, AngleIdx'
DO i = 1, nMod
    WRITE(io, '(i13,100i13)') i, ModRay(i)%NextRayIdx, ModRay(i)%InOutSurf, ModRay(i)%iAziAngIdx
ENDDO
!--- Assembly Ray
WRITE(io, '(a)') '! Assembly Ray- '
DO i = 1, nAsy
    WRITE(io, '(i13,100i13)') i, AsyRay(i)%AziAngIdx, AsyRay(i)%ModRayIdx, AsyRay(i)%NextRayIdx, AsyRay(i)%InOutSurf, AsyRay(i)%nCellRay, AsyRay(i)%PartialAsyRayFlag, AsyRay(i)%PartialAsyRayIdx
    DO j = 1, AsyRay(i)%nCellRay
        WRITE(io, '(a,i13,100i13)') '             ', j, AsyRay(i)%PinRaySurf(:,j), AsyRay(i)%PinIdx(j), AsyRay(i)%PinRayIdx(j)
    ENDDO
ENDDO
!--- Core Ray
WRITE(io, '(a)') '! Core Ray- '
DO i = 1, nCore
    WRITE(io, '(i13,100i13)') i, CoreRay(i)%nRay, CoreRay(i)%iAng, CoreRay(i)%InOutSurf
    DO j = 1, CoreRay(i)%nRay
        WRITE(io, '(a,i13,100i13)') '             ', j, CoreRay(i)%AsyRayIdx(j), CoreRay(i)%AsyIdx(j), CoreRay(i)%IXA(j), CoreRay(i)%IYA(j)
    ENDDO
ENDDO
!--- Rotational Ray
WRITE(io, '(a)') '! Rotational Ray'
DO i = 1, nRot
    WRITE(io, '(i13,100i13)') i, RotRay(i)%nRay, RotRay(i)%nSeg
    DO j = 1, RotRay(i)%nRay
        WRITE(io, '(a,i13,100i13)') '             ', j, RotRay(i)%RayIdx(j), RotRay(i)%Dir(j)
    ENDDO
ENDDO
CONTINUE
!--- Base Cell Ray
WRITE(io, '(a)') '! Base Cell Ray'
DO i = 1, nCellType
    IF( i .EQ. CellRayBase(i)%CellType )THEN
        WRITE(io, '(i13,100i13)') i, CellRayBase(i)%CellType, CellRayBase(i)%nLine
        DO j = 1, nAzi
            WRITE(io, '(a,i13,100i13)') '             ', j, CellRayBase(i)%AziAngIdx(j), CellRayBase(i)%nLines(j)
        ENDDO
        DO j = 1, CellRayBase(i)%nline
            WRITE(io, '(a,i13,100i13)') '             ', j, CellRayBase(i)%CellRay(j)%IdxSt, CellRayBase(i)%CellRay(j)%IdxEnd, CellRayBase(i)%CellRay(j)%nSeg
            DO  k = 1, CellRayBase(i)%CellRay(j)%nSeg
                WRITE(io, '(2a,i13,i13,e13.6)') '             ','             ', k, CellRayBase(i)%CellRay(j)%LocalFsrIdx(k), CellRayBase(i)%CellRay(j)%LenSeg(k)
            ENDDO
        ENDDO
    ELSE
        WRITE(io, '(i13,100i13)') i, CellRayBase(i)%CellType, 0
    ENDIF
ENDDO
WRITE(io, '(a)') '.'
CLOSE(io)

!--- GEOM output
!--- File Open
io=48
WRITE(fn,'(A)') TRIM(caseid)
WRITE(fn,'(A,A)') TRIM(fn),'.gsf' !Ray Structure File
OPEN(unit=io, file = fn, status = 'replace')

WRITE(io, '(a)') '! Base Geometry Information'
WRITE(io, '(100i13)') core%nxy, core%nxya, core%nCoreFxr, core%nCoreFsr, core%npintype, core%nCellType, core%nAsyType, core%nz
!--- Boundary Condition
WRITE(io, '(a)') '! Boundary Info'
DO i = 1, 6
    j=0
    IF(core%radsym(i)) j=1
    WRITE(io,'(2i13)') j, core%radbc(i)
ENDDO
!--- Axial Mesh height
DO i = 1, core%nz
    WRITE(io,'(i13,f13.6)') i, core%hz(i)
ENDDO

!------ PROBLEM DEFINED VARIABLES
WRITE(io, '(a)') '! Pin '
DO i = 1, core%nxy
    WRITE(io,'(i13,100i13)') i, core%pin(i)%ix, core%pin(i)%iy, core%pin(i)%iasy, core%pin(i)%ipin, core%pin(i)%pintype, core%pin(i)%asytype, core%pin(i)%ncell, core%pin(i)%nFsrMax, core%pin(i)%nFxrMax
    WRITE(io,'(i13,100i13)') core%pin(i)%nBD, core%pin(i)%FsrIdxSt, core%pin(i)%FxrIdxSt, core%pin(i)%lfuel, core%pin(i)%lgt, core%pin(i)%lradref, core%pin(i)%lBaffle, core%pin(i)%lBarrel, core%pin(i)%lGd
    DO j = 1, core%pin(i)%nBD
        WRITE(io,'(a,3i13)') '             ', j, core%Pin(i)%neighidx(j), core%Pin(i)%neighsurfidx(j)
    ENDDO
    DO j = 1, core%nz
        WRITE(io,'(a,2i13,f13.6,i13)') '             ', j, core%pin(i)%cell(j), core%pinvol(i,j), core%pin(i)%lmox(j)
    ENDDO
ENDDO
!--- Assembly
nbd=4
IF( lHex ) nbd=6
WRITE(io, '(a)') '! Assembly '
DO i = 1, Core%nxya
    WRITE(io,'(i13,100i13)') i, Core%Asy(i)%ixa, Core%Asy(i)%iya, Core%Asy(i)%asytype, Core%Asy(i)%partialAsyFlag
    DO j = 1, nbd
        WRITE(io,'(a,2i13)') '             ', core%asy(i)%neighIdx(j)
    ENDDO
    DO j = 1, Core%AsyInfo(Core%Asy(i)%asytype)%nxy
        WRITE(io,'(a,2i13)') '             ', j, core%asy(i)%GlobalPinIdx(j)
    ENDDO
ENDDO

!------ PRE DEFINED VALUES ---
!--- CellInfo
WRITE(io, '(a)') '! Cell Info'
DO i = 1,  core%ncellType
    WRITE(io,'(i13,100i13)') i, Core%CellInfo(i)%nDivAzi, Core%CellInfo(i)%nFsr, Core%CellInfo(i)%nFXR,             &
                                Core%CellInfo(i)%lEmpty, Core%CellInfo(i)%luse, Core%CellInfo(i)%lrect,             &
                                Core%CellInfo(i)%lres, Core%CellInfo(i)%lfuel, Core%CellInfo(i)%lgd
    WRITE(io,'(i13,100i13)') Core%CellInfo(i)%lmox, Core%CellInfo(i)%lccell, Core%CellInfo(i)%lcentx,               &
                             Core%CellInfo(i)%lcenty, Core%CellInfo(i)%lcentxy, Core%CellInfo(i)%lgap
    DO j = 1, Core%CellInfo(i)%nFxr
        WRITE(io,'(a,3i13)') '             ', j, Core%CellInfo(i)%FxrIdxSt(j), Core%CellInfo(i)%nFsrInFxr(j)
    ENDDO
    DO j = 1, Core%CellInfo(i)%nFsr
        WRITE(io,'(1p,a,2i13,e13.6)') '             ', j, Core%CellInfo(i)%ireg(j), core%CellInfo(i)%vol(j)
    ENDDO
ENDDO

!--- Pin Info
WRITE(io, '(a)') '! Pin Info'
DO i = 1,  core%npintype
    WRITE(io,'(i13,100i13)') i, Core%PinInfo(i)%ncell, Core%PinInfo(i)%nFsrMax, Core%PinInfo(i)%nFXRmax,            &
                                Core%PinInfo(i)%lEmpty, Core%PinInfo(i)%luse, Core%PinInfo(i)%lfuel,                &
                                Core%PinInfo(i)%lcentx, Core%PinInfo(i)%lcenty, Core%PinInfo(i)%lcentxy,            &
                                Core%PinInfo(i)%lgap
    DO j = 1, core%PinInfo(i)%ncell
        WRITE(io,'(a,2i13)') '             ', j, Core%PinInfo(i)%Cell(j)
    ENDDO
ENDDO



WRITE(io, '(a)') '.'
CLOSE(io)

ENDSUBROUTINE



SUBROUTINE PrintMOCUnderRelaxation(io, FmInfo, GroupInfo, PE)
USE PARAM
USE TYPEDEF,   ONLY : FmInfo_Type, GroupInfo_Type, PE_TYPE
IMPLICIT NONE
INTEGER :: io
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE

INTEGER :: i, j, ibeg, iend
WRITE(io, *)
WRITE(io, '(x, A)') 'MOC Under-Relaxation Factor'
i=0
DO
  i = i + 1
  ibeg = 10 * (i-1) + 1;iend = 10 * i
  iend = min(iend, GroupInfo%ng)
  WRITE(io, '(x,A5,x, I3,A3,x,I3, x, A, 20F7.3)') 'Group', ibeg,'-', iend, ':', (FmInfo%w(j), j=ibeg, iend)
  IF(iend .EQ. GroupInfo%ng) EXIT
ENDDO
WRITE(io,*)
END SUBROUTINE




SUBROUTINE PrintFSRPhi(Core, THInfo, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
USE param,   only : CKELVIN
USE TYPEDEF,      ONLY : CoreInfo_Type,  THInfo_type, CmInfo_Type, FmInfo_type,    &
                         PE_TYPE,        GroupInfo_Type, PinXS_TYPE, &
                         AsyInfo_type, Asy_Type, Pin_Type, Cell_Type, FxrInfo_type , XsMac_type
USE CNTL,         ONLY : nTracerCntl_Type,           OutpCntl_Type
USE FILES,        ONLY : LOCALFN, caseid
USE IOUTIL,       ONLY : OpenFile
USE Depl_Mod,  ONLY : DeplCntl
USE BenchXs,        ONLY : XsBaseBen, xssm1ben, xssm2ben, xssm3ben
USE XsUtil_mod,   ONLY : AllocXsMac
USE TH_Mod,           ONLY : GetPinFuelTemp
IMPLICIT NONE
!
TYPE(CoreInfo_Type) :: Core
TYPE(THInfo_Type) :: THInfo
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(XsMac_Type) :: XsMac
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
REAL :: TempRef
INTEGER :: io, io2

REAL, POINTER :: MacXsSm(:,:)  !IsoMacXsSm(isosize,ng,ng)
REAL, POINTER :: MacXsSmP1(:,:) !IsoMacXsSmP1(isosize,ng,ng)
!
CHARACTER(256) :: PinXsFn
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!
INTEGER :: nxy, nz, myzb, myze, ng, nxya, nlocalFxr, nIsoInFXR, nFSRinFXR, nchi
INTEGER :: ixy, ixyl, iz, ig, igg, igb, ige, ifxr, ifsr, icel, iasytype, iasy, fxridxst, fsridxst, fidx
INTEGER :: i, j, k, idx, itype, imix, ifsrl
REAL :: chi, absr, fisr
LOGICAL :: lxslib
!
CHARACTER(256) :: fn, fn2, matname
INTERFACE
    SUBROUTINE GenFxrMacXS(XsMac, MacXsSm, MacXsSmP1, Fxr, ifxr, Tempref, ig1, ig2, Core, ipin, iz, GroupInfo, nTRACERCntl, PE)
    USE PARAM
    USE TYPEDEF,      ONLY : CoreInfo_Type, XsMac_Type,  FxrInfo_Type,  GroupInfo_Type,  PE_Type
    USE CNTL,          ONLY : nTracerCntl_Type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(XsMac_Type) :: XsMac
    TYPE(FxrInfo_Type),pointer :: Fxr(:,:)
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    INTEGER :: ig1, ig2, ipin, iz, ifxr
    REAL :: TempRef
    REAL :: MacXsSm(GroupInfo%ng,GroupInfo%ng)
    REAL :: MacXsSmP1(GroupInfo%ng,GroupInfo%ng)
    END SUBROUTINE
END INTERFACE

io=45
io2=46

AsyInfo => Core%AsyInfo; Asy => Core%Asy
Pin => Core%Pin; Cell => Core%CellInfo
nz = Core%nz; ng = GroupInfo%ng
nxya= Core%nxya ! n-assemblies
nchi= GroupInfo%nChi
ALLOCATE(MacXsSm(ng,ng))
ALLOCATE(MacXsSmP1(ng,ng))


lxslib=nTracerCntl%lxslib
idx=0

DO iz = PE%myzb, PE%myze
    DO iasy = 1, nxya
        iasytype=Asy(iasy)%AsyType
        nxy=AsyInfo(iasytype)%nxy
        DO ixyl = 1, nxy
            WRITE(fn,'(A)') TRIM(caseid)
            IF( DeplCntl%nburnupstep .GT. 1 )THEN
                WRITE(fn,'(A,A)') TRIM(fn),'_d'
                IF( DeplCntl%nowstep .LT. 10 )THEN
                    WRITE(fn,'(A,A,i1)') TRIM(fn),'0',DeplCntl%nowstep
                ELSE
                    WRITE(fn,'(A,i2)') TRIM(fn),DeplCntl%nowstep
                ENDIF
            ENDIF
            WRITE(fn,'(A,A)') TRIM(fn),'_z'
            IF( iz .LT. 10 )THEN
                WRITE(fn,'(A,A,i1)') TRIM(fn),'0',iz
            ELSE
                WRITE(fn,'(A,i2)') TRIM(fn),iz
            ENDIF

            WRITE(fn,'(A,A)') TRIM(fn),'_a'
            IF( iasy .LT. 10 )THEN
                WRITE(fn,'(A,A,i1)') TRIM(fn),'0',iasy
            ELSEIF( iasy .LT. 100 )THEN
                WRITE(fn,'(A,i2)') TRIM(fn),iasy
            ELSE
                WRITE(fn,'(A,i3)') TRIM(fn),iasy
            ENDIF

            WRITE(fn,'(A,A)') TRIM(fn),'_p'
            IF( ixyl .LT. 10 )THEN
                WRITE(fn,'(A,A,i1)') TRIM(fn),'00',ixyl
            ELSEIF( ixyl .LT. 100 )THEN
                WRITE(fn,'(A,A,i2)') TRIM(fn),'0',ixyl
            ELSE
                WRITE(fn,'(A,i3)') TRIM(fn),ixyl
            ENDIF
            WRITE(fn,'(A,A)') TRIM(fn),'.spf'
            OPEN(unit=io, file = fn, status = 'replace')

            ixy = Asy(iasy)%GlobalPinIdx(ixyl)
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            WRITE(io,'(i3)') Cell(icel)%BaseCellStr
            IF( lxslib )    Tempref = GetPinFuelTemp(Core, FmInfo%FXR, iz, ixy)
            DO j = 1, nLocalFxr
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                imix=myFXR%imix
                nFsrInFxr = myFxr%nFsrInFxr
                !FsrIdxSt = myFxr%FsrIdxSt
                nisoInFxr = myFxr%niso
                IF( lxslib )THEN
                    CALL GenFxrMacXS(XsMac, MacXsSm, MacXsSmP1, FmInfo%FXR, ifxr, Tempref, 1, ng, Core, ixy, iz, GroupInfo, nTRACERCntl, PE)
                ELSE
                    IF(.NOT. XsMac%lalloc) THEN
                        XsMac%ng = ng
                        CALL AllocXsMac(XsMac)
                    ENDIF
                    itype=myFxr%imix
                    CALL xsbaseBen(itype, 1, ng, 1, ng, nTracerCntl%lscat1, XsMac)
                ENDIF
                fsridxst=myfxr%fsridxst
                nfsrinfxr=myfxr%nfsrinfxr
                ifsrl=0
                DO ifsr=fsridxst, fsridxst+nfsrinfxr-1
                    absr=0
                    fisr=0
                    ifsrl=ifsrl+1
                    DO ig = 1, ng
                        absr=absr+fminfo%phis(ifsr,iz,ig)*XsMac%xsmaca(ig)
                        fisr=fisr+fminfo%phis(ifsr,iz,ig)*XsMac%xsmacnf(ig)
                    ENDDO
                    WRITE(io,'(1p, 3i13, 2e14.7)') imix, j, ifsrl, fisr, absr
                ENDDO
            ENDDO ! FXR in Asy
            WRITE(io,'(a1)') '.'
            close(io)
        ENDDO ! Pin in Asy
    ENDDO ! Global Asy
ENDDO !Z-axial plane
END SUBROUTINE



SUBROUTINE PrintPinFisPhi(Core, THInfo, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
USE param,   only : CKELVIN
USE TYPEDEF,      ONLY : CoreInfo_Type,  THInfo_type, CmInfo_Type, FmInfo_type,    &
                         PE_TYPE,        GroupInfo_Type, PinXS_TYPE, &
                         AsyInfo_type, Asy_Type, Pin_Type, Cell_Type, FxrInfo_type , XsMac_type
USE CNTL,         ONLY : nTracerCntl_Type,           OutpCntl_Type
USE FILES,        ONLY : LOCALFN, caseid
USE IOUTIL,       ONLY : OpenFile
USE Depl_Mod,  ONLY : DeplCntl
USE BenchXs,        ONLY : XsBaseBen, xssm1ben, xssm2ben, xssm3ben
USE XsUtil_mod,   ONLY : AllocXsMac
USE TH_Mod,           ONLY : GetPinFuelTemp
IMPLICIT NONE
!
TYPE(CoreInfo_Type) :: Core
TYPE(THInfo_Type) :: THInfo
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(XsMac_Type) :: XsMac
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
REAL :: TempRef
INTEGER :: io, io2

REAL, POINTER :: MacXsSm(:,:)  !IsoMacXsSm(isosize,ng,ng)
REAL, POINTER :: MacXsSmP1(:,:) !IsoMacXsSmP1(isosize,ng,ng)

REAL, POINTER :: pin3dpwr(:,:,:) ! (x,y,z)
REAL, POINTER :: pin2dpwr(:,:) ! (x,y)
REAL, POINTER :: pin3dabs(:,:,:) ! (x,y,z)
REAL, POINTER :: pin2dabs(:,:) ! (x,y)
REAL, POINTER :: asy3dpwr(:,:,:) ! (x,y,z)
REAL, POINTER :: asy2dpwr(:,:) ! (x,y)
!
CHARACTER(256) :: PinXsFn
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!
INTEGER :: nxy, nz, myzb, myze, ng, nxya, nlocalFxr, nIsoInFXR, nFSRinFXR, nchi
INTEGER :: ixy, ixyl, iz, ig, igg, igb, ige, ifxr, ifsr, icel, iasytype, iasy, fxridxst, fsridxst, fidx
INTEGER :: i, j, k, idx, itype, imix, ifsrl, ix, iy, ixa, iya
INTEGER :: npinx, npiny, nasyx, nasyy
REAL :: chi, absr, fisr
LOGICAL :: lxslib, master
REAL :: hz
CHARACTER(256) :: fn, fn2, matname

INTERFACE
    SUBROUTINE GenFxrMacXS(XsMac, MacXsSm, MacXsSmP1, Fxr, ifxr, Tempref, ig1, ig2, Core, ipin, iz, GroupInfo, nTRACERCntl, PE)
    USE PARAM
    USE TYPEDEF,      ONLY : CoreInfo_Type, XsMac_Type,  FxrInfo_Type,  GroupInfo_Type,  PE_Type
    USE CNTL,          ONLY : nTracerCntl_Type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(XsMac_Type) :: XsMac
    TYPE(FxrInfo_Type),pointer :: Fxr(:,:)
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    INTEGER :: ig1, ig2, ipin, iz, ifxr
    REAL :: TempRef
    REAL :: MacXsSm(GroupInfo%ng,GroupInfo%ng)
    REAL :: MacXsSmP1(GroupInfo%ng,GroupInfo%ng)
    END SUBROUTINE
END INTERFACE

io=45
io2=46

AsyInfo => Core%AsyInfo; Asy => Core%Asy
Pin => Core%Pin; Cell => Core%CellInfo
nz = Core%nz; ng = GroupInfo%ng
npinx = Core%nx  ! npinx, npiny
npiny = Core%ny  ! npinx, npiny
nasyx = Core%nxa
nasyy = Core%nya
nxya= Core%nxya ! n-assemblies
nchi= GroupInfo%nChi
ALLOCATE(MacXsSm(ng,ng))
ALLOCATE(MacXsSmP1(ng,ng))

ALLOCATE(pin3dpwr(npinx,npiny,nz),pin2dpwr(npinx,npiny))
ALLOCATE(pin3dabs(npinx,npiny,nz),pin2dabs(npinx,npiny))
ALLOCATE(asy3dpwr(npinx,npiny,nz),asy2dpwr(npinx,npiny))
pin3dpwr=0;pin2dpwr=0;
pin3dabs=0;pin2dabs=0;
asy3dpwr=0;asy2dpwr=0;


lxslib=nTracerCntl%lxslib
idx=0

DO iz = PE%myzb, PE%myze
!DO iz = 1, Core%nz
    hz=Core%hz(iz)
    DO iasy = 1, nxya
        iasytype=Asy(iasy)%AsyType
        nxy=AsyInfo(iasytype)%nxy
        ixa=Asy(iasy)%ixa
        iya=Asy(iasy)%iya
        DO ixyl = 1, nxy
            ixy = Asy(iasy)%GlobalPinIdx(ixyl)
            ix = Pin(ixy)%ix
            iy = Pin(ixy)%iy
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            IF( lxslib )    Tempref = GetPinFuelTemp(Core, FmInfo%FXR, iz, ixy)
            DO j = 1, nLocalFxr
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                imix=myFXR%imix
                nFsrInFxr = myFxr%nFsrInFxr
                nisoInFxr = myFxr%niso
                IF( lxslib )THEN
                    CALL GenFxrMacXS(XsMac, MacXsSm, MacXsSmP1, FmInfo%FXR, ifxr, Tempref, 1, ng, Core, ixy, iz, GroupInfo, nTRACERCntl, PE)
                ELSE
                    IF(.NOT. XsMac%lalloc) THEN
                        XsMac%ng = ng
                        CALL AllocXsMac(XsMac)
                    ENDIF
                    itype=myFxr%imix
                    CALL xsbaseBen(itype, 1, ng, 1, ng, nTracerCntl%lscat1, XsMac)
                ENDIF
                fsridxst=myfxr%fsridxst
                nfsrinfxr=myfxr%nfsrinfxr
                ifsrl=0
                DO ifsr=fsridxst, fsridxst+nfsrinfxr-1
                    absr=0
                    fisr=0
                    ifsrl=ifsrl+1
                    DO ig = 1, ng
                        absr=absr+fminfo%phis(ifsr,iz,ig)*XsMac%xsmaca(ig)
                        fisr=fisr+fminfo%phis(ifsr,iz,ig)*XsMac%xsmacf(ig)
                    ENDDO
                    pin3dpwr(ix,iy,iz)=pin3dpwr(ix,iy,iz)+fisr*myFxr%area*hz
                    pin3dabs(ix,iy,iz)=pin3dabs(ix,iy,iz)+absr*myFxr%area*hz
                    pin2dpwr(ix,iy)=pin2dpwr(ix,iy)+fisr*myFxr%area*hz
                ENDDO
            ENDDO ! FXR in pin
            asy3dpwr(ixa,iya,iz)=asy3dpwr(ixa,iya,iz)+pin3dpwr(ix,iy,iz)
            asy2dpwr(ixa,iya)=asy2dpwr(ixa,iya)+pin3dpwr(ix,iy,iz)
        ENDDO ! Pin in Asy
    ENDDO ! Global Asy
ENDDO !Z-axial plane

!DO iz = 1, Core%nz
DO iz = PE%myzb, PE%myze

    if( iz .lt. 10 )then
        WRITE(fn,'(A,A,i1,A)') TRIM(caseid),'_3dpin_0',iz ,'.fph'
    else
        WRITE(fn,'(A,A,i2,A)') TRIM(caseid),'_3dpin_',iz ,'.fph'
    endif
    OPEN(unit=io, file = fn, status = 'replace')
    !WRITE(io,'(i)') iz
    DO iy = 1, npiny
        WRITE(io,'(i5, 1000e14.7)') iy, pin3dpwr(:,iy,iz)
    ENDDO
    WRITE(io,'(a)') ' '
    close(io)
ENDDO

END SUBROUTINE

SUBROUTINE TreatDancoffData(Core,DancoffDist, PE)
USE TYPEDEF,        ONLY : CoreInfo_Type, DancoffDist_Type, PE_Type, ResVarPin_Type, AsyInfo_Type, Asy_Type
#ifdef MPI_ENV
USE MPIComm_mod,    ONLY : SENDRECV
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(DancoffDist_Type) :: DancoffDist
TYPE(PE_Type) :: PE
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
LOGICAL :: master, slave
INTEGER :: nxy, nxymax,nz, nxya, myzb, myze
INTEGER :: nAsyType, AsyType
INTEGER :: ixy, iz, iasy
INTEGER :: i, iz1, iz2

master = PE%Cmfdmaster; slave = PE%Cmfdslave
AsyInfo => Core%AsyInfo; Asy => Core%Asy
ResVarPin => Core%ResVarPin

nAsyType = Core%nAsyType
nxya = Core%nxya; nz = Core%nz
myzb = PE%myzb; myze = PE%myze

nxy = 0
DO i = 1, nAsyType
  nxy = max(AsyInfo(i)%nxy, nxy)
ENDDO
nxymax = nxy
IF(Master) THEN
  ALLOCATE(DancoffDist%Dancoff3D(nxy, nxya, nz))
ELSE
  ALLOCATE(DancoffDist%Dancoff3D(nxy, nxya, myzb:myze))
ENDIF

!Communicate

DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxy = AsyInfo(AsyType)%nxy
  Do iz = myzb, myze   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    DO i = 1, nxy  !Cell Sweep
      ixy = Asy(iasy)%GlobalPinIdx(i)
      DancoffDist%Dancoff3D(i, iasy, iz) = ResVarPin(ixy,iz)%Dancoff
    ENDDO
  ENDDO
ENDDO !

#ifdef MPI_ENV 
IF(Master) THEN
  DO i = 1, PE%nCmfdProc - 1
    iz1 = PE%AxDomRange(1, i); iz2 = PE%AxDomRange(2, i)
    CALL SendRecv(DancoffDist%Dancoff3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya,    &
                  iz2 - iz1 + 1, i, .FALSE., PE%MPI_CMFD_COMM )
  ENDDO
ELSE
    iz1 = myzb; iz2 = myze
    CALL SendRecv(DancoffDist%Dancoff3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya,    &
                  iz2 - iz1 + 1, 0, .TRUE., PE%MPI_CMFD_COMM )
ENDIF

IF(slave) THEN
  NULLIFY(AsyInfo); NULLIFY(Asy)
  RETURN
ENDIF
#endif

END SUBROUTINE