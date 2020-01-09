SUBROUTINE PrintFuelAvgTemp(io, Core, THInfo, THVar)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,      ThInfo_Type ,      CoolantTh_Type,  &
                    FuelTH_Type,        THVar_Type,                         &
                    AsyInfo_Type,       Asy_Type
IMPLICIT NONE
INTEGER :: io

TYPE(CoreInfo_Type) :: Core
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)

TYPE(ThInfo_Type) :: ThInfo
TYPE(FuelTH_Type), POINTER :: FuelTH(:)
TYPE(THVar_Type) :: THVar

REAL, POINTER :: TavgFuel(:)
REAL, POINTER :: OutPut(:, :), OutPut1(:, :)
INTEGER :: iz, ixy, ixya, ixa, iya, iasy, AsyType
INTEGER :: nz, nxy, nxya, nxa, nya, ndat
INTEGER :: navg

nz = Core%nz; nxy = Core%nxy
nxya = Core%nxya;ndat = Core%nxyc
AsyInfo => Core%AsyInfo; Asy => Core%Asy

FuelTH => THInfo%FuelTH
navg = THVar%npr5

ALLOCATE(OutPut1(nz, nDat))

301 FORMAT('Assembly', I5, x,'at' ,x,'(', I7, I7,' )')

!ALLOCATE(hcool(nz)); hcool = 0
WRITE(io, '(A)') '- Fuel Pin Average Temperature -'
WRITE(io, *) 
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxya = AsyInfo(AsyType)%nxy
  nxa = AsyInfo(AsyType)%nx; nya = AsyInfo(AsyType)%ny
  ALLOCATE(OutPut(nz, nxya))
  DO ixya = 1, nxya
    ixy = Asy(iasy)%GlobalPinIdx(ixya)
    OutPut(1:nz, ixya) = 0.
    IF(.NOT. FuelTh(ixy)%lfuel) CYCLE
    TavgFuel => FuelTh(ixy)%tfuel(navg, :)
    OutPut(1:nz, ixya) =TavgFuel(1:nz)
    NULLIFY(TavgFuel)
  ENDDO
  DO iz = 1, nz
    CALL EditOutputDat(Core, AsyType, OutPut(iz, 1:nxya), OutPut1(iz, 1:ndat), nxya, ndat)
  ENDDO
  ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
  WRITE(io, 301) iasy, ixa, iya
  WRITE(io, *) 
  nxya = Core%nxyc0; nxa = Core%nxc0; nya = Core%nyc0
  CAll PrintTHResults(io, OutPut1(1:nz, 1:nxya), nxya, nxa, nya, nz)
  DEALLOCATE(OutPut)
ENDDO

NULLIFY(AsyInfo)
NULLIFY(Asy)
NULLIFY(FuelTH)
DEALLOCATE(OutPut1)
END SUBROUTINE

SUBROUTINE PrintFuelCentTemp(io, Core, THInfo, THVar)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,      ThInfo_Type ,      CoolantTh_Type,  &
                    FuelTH_Type,        THVar_Type,                         &
                    AsyInfo_Type,       Asy_Type
IMPLICIT NONE
INTEGER :: io

TYPE(CoreInfo_Type) :: Core
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)

TYPE(ThInfo_Type) :: ThInfo
!TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
TYPE(FuelTH_Type), POINTER :: FuelTH(:)
TYPE(THVar_Type) :: THVar

REAL, POINTER :: Tfuel(:)
REAL, POINTER :: OutPut(:, :), OutPut1(:, :)
INTEGER :: iz, ixy, ixya, ixa, iya, iasy, AsyType
INTEGER :: nz, nxy, nxya, nxa, nya, ndat
INTEGER :: j, jbeg, jend

nz = Core%nz; nxy = Core%nxy
nxya = Core%nxya;ndat = Core%nxyc
AsyInfo => Core%AsyInfo; Asy => Core%Asy

!CoolantTH => THInfo%CoolantTH; 
FuelTH => THInfo%FuelTH
!TCool => THInfo%Tcool;

301 FORMAT('Assembly', I5, x,'at' ,x,'(', I7, I7,' )')

!ALLOCATE(hcool(nz)); hcool = 0
ALLOCATE(OutPut1(nz, ndat))
WRITE(io, '(A)') '- Fuel Center-line Temperature -'
WRITE(io, *) 
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxya = AsyInfo(AsyType)%nxy
  nxa = AsyInfo(AsyType)%nx; nya = AsyInfo(AsyType)%ny
  ALLOCATE(OutPut(nz, nxya))
  DO ixya = 1, nxya
    ixy = Asy(iasy)%GlobalPinIdx(ixya)
    OutPut(1:nz, ixya) = 0.
    IF(.NOT. FuelTh(ixy)%lfuel) CYCLE
    Tfuel => FuelTh(ixy)%Tfuel(1, :)
    !TavgFuel => FuelTh(ixy)%tfuel(navg, :)
    OutPut(1:nz, ixya) =TFuel(1:nz)
    NULLIFY(TFuel)
  ENDDO
  DO iz = 1, nz
    CALL EditOutputDat(Core, AsyType, OutPut(iz, 1:nxya), OutPut1(iz, 1:ndat), nxya, ndat)
  ENDDO
  ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
  WRITE(io, 301) iasy, ixa, iya
  WRITE(io, *) 
  nxya = Core%nxyc0; nxa = Core%nxc0; nya = Core%nyc0
  CAll PrintTHResults(io, OutPut1, nxya, nxa, nya, nz)
  DEALLOCATE(OutPut)
ENDDO
!WRITE(*, '(A, F10.3)'), 'TEFF=', FuelTh(1)%teff(1)
NULLIFY(AsyInfo)
NULLIFY(Asy)
NULLIFY(FuelTH)
DEALLOCATE(Output1)
END SUBROUTINE

SUBROUTINE PrintFuelSurfTemp(io, Core, THInfo, THVar)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,      ThInfo_Type ,      CoolantTh_Type,  &
                    FuelTH_Type,        THVar_Type,                         &
                    AsyInfo_Type,       Asy_Type
IMPLICIT NONE
INTEGER :: io

TYPE(CoreInfo_Type) :: Core
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)

TYPE(ThInfo_Type) :: ThInfo
!TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
TYPE(FuelTH_Type), POINTER :: FuelTH(:)
TYPE(THVar_Type) :: THVar

REAL, POINTER :: Tfuel(:)
REAL, POINTER :: OutPut(:, :), OutPut1(:, :)
INTEGER :: iz, ixy, ixya, ixa, iya, iasy, AsyType
INTEGER :: nz, nxy, nxya, nxa, nya, ndat
INTEGER :: j, jbeg, jend

nz = Core%nz; nxy = Core%nxy
nxya = Core%nxya;ndat = Core%nxyc
AsyInfo => Core%AsyInfo; Asy => Core%Asy

!CoolantTH => THInfo%CoolantTH; 
FuelTH => THInfo%FuelTH
!TCool => THInfo%Tcool;

301 FORMAT('Assembly', I5, x,'at' ,x,'(', I7, I7,' )')

!ALLOCATE(hcool(nz)); hcool = 0
ALLOCATE(OutPut1(nz, ndat))
WRITE(io, '(A)') '- Fuel Surface Temperature -'
WRITE(io, *) 
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxya = AsyInfo(AsyType)%nxy
  nxa = AsyInfo(AsyType)%nx; nya = AsyInfo(AsyType)%ny
  ALLOCATE(OutPut(nz, nxya))
  DO ixya = 1, nxya
    ixy = Asy(iasy)%GlobalPinIdx(ixya)
    OutPut(1:nz, ixya) = 0.
    IF(.NOT. FuelTh(ixy)%lfuel) CYCLE
    Tfuel => FuelTh(ixy)%Tfuel(ThVar%npr1, :)
    OutPut(1:nz, ixya) =TFuel(1:nz)
    NULLIFY(TFuel)
  ENDDO
  DO iz = 1, nz
    CALL EditOutputDat(Core, AsyType, OutPut(iz, 1:nxya), OutPut1(iz, 1:ndat), nxya, ndat)
  ENDDO
  ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
  WRITE(io, 301) iasy, ixa, iya
  WRITE(io, *) 
  nxya = Core%nxyc0; nxa = Core%nxc0; nya = Core%nyc0
  CAll PrintTHResults(io, OutPut1, nxya, nxa, nya, nz)
  DEALLOCATE(OutPut)
ENDDO
!WRITE(*, '(A, F10.3)'), 'TEFF=', FuelTh(1)%teff(1)
NULLIFY(AsyInfo)
NULLIFY(Asy)
NULLIFY(FuelTH)
DEALLOCATE(Output1)
END SUBROUTINE

SUBROUTINE PrintCoolantTemp(io, Core, THInfo, THVar)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,      ThInfo_Type ,      CoolantTh_Type,  &
                    FuelTH_Type,        THVar_Type,                         &
                    AsyInfo_Type,       Asy_Type
IMPLICIT NONE
INTEGER :: io

TYPE(CoreInfo_Type) :: Core
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)

TYPE(ThInfo_Type) :: ThInfo
TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
!TYPE(FuelTH_Type), POINTER :: FuelTH(:)
TYPE(THVar_Type) :: THVar

REAL, POINTER :: Tcool(:)
REAL, POINTER :: OutPut(:, :), OutPut1(:, :)
INTEGER :: iz, ixy, ixya, ixa, iya, iasy, AsyType
INTEGER :: nz, nxy, nxya, nxa, nya, ndat
INTEGER :: j, jbeg, jend

nz = Core%nz; nxy = Core%nxy
nxya = Core%nxya; ndat = Core%nxyc
AsyInfo => Core%AsyInfo; Asy => Core%Asy

CoolantTH => THInfo%CoolantTH; 
!FuelTH => THInfo%FuelTH
!TCool => THInfo%Tcool;

301 FORMAT('Assembly', I5, x,'at' ,x,'(', I7, I7,' )')

ALLOCATE(OutPut1(nz, ndat))
WRITE(io, '(A)') '- Coolant Temperature -'
WRITE(io, *) 
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxya = AsyInfo(AsyType)%nxy
  nxa = AsyInfo(AsyType)%nx; nya = AsyInfo(AsyType)%ny
  ALLOCATE(OutPut(nz, nxya))
  DO ixya = 1, nxya
    ixy = Asy(iasy)%GlobalPinIdx(ixya)
    OutPut(1:nz, ixya) = 0.
    Tcool => THInfo%Tcool(:, ixy)
    OutPut(1:nz, ixya) =Tcool(1:nz)
    NULLIFY(Tcool)
  ENDDO
  DO iz = 1, nz
    CALL EditOutputDat(Core, AsyType, OutPut(iz, 1:nxya), OutPut1(iz, 1:ndat), nxya, ndat)
  ENDDO
  ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
  WRITE(io, 301) iasy, ixa, iya
  WRITE(io, *) 
  nxya = Core%nxyc0; nxa = Core%nxc0; nya = Core%nyc0
  CAll PrintTHResults(io, OutPut1, nxya, nxa, nya, nz)
  Deallocate(OutPut)
ENDDO
NULLIFY(AsyInfo)
NULLIFY(Asy)
NULLIFY(CoolantTH)
DEALLOCATE(Output1)
END SUBROUTINE


SUBROUTINE PrintTHResults(io, Output, nxya, nxa, nya, nz)
USE PARAM
IMPLICIT NONE
INTEGER :: io, nxya, nxa, nya, nz
REAL :: OutPut(nz, nxya)
INTEGER :: iya, ixa, iz
INTEGER :: j, jbeg, jend

WRITE(io, '(A7, A7, 4x, A10)') 'Y-Idx', 'Z-Idx','--> X-Idx'
WRITE(io, '(14x, 200I8)') (j, j = 1, nxa)
DO iya = 1, nya
  jbeg = nxa * (iya-1) + 1
  jend = nxa * iya
  DO iz = 1, nz
    IF(iz .EQ. 1) THEN
      WRITE(io, '(2I5, 4x, 200F8.2)') iya, iz, (OutPut(iz, j), j = jbeg, jend)
    ELSE
      WRITE(io, '(I10, 4x, 200F8.2)') iz, (OutPut(iz, j), j= jbeg, jend)
    ENDIF
  ENDDO
  WRITE(io, *) 
ENDDO 
WRITE(io, *)
END SUBROUTINE

SUBROUTINE PrintCoolantNDAvgTemp(Core, FmInfo)
USE TYPEDEF, ONLY : CoreInfo_Type, FmInfo_Type, AsyInfo_Type, Asy_Type, Cell_Type, Pin_Type, FxrInfo_Type
USE FILES,   ONLY : io16, caseid
USE ioutil,  ONLY : openfile
USE Depl_Mod,  ONLY : DeplCntl
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: PIN(:), mypin
TYPE(FxrInfo_Type), POINTER :: myFxr

INTEGER :: iz,nz,iasy,nxya,nx,ny,hnx,hny,ix,iy,ixy,nfxr,icel,ir,ifxr,iso,niso,iasytype
REAL :: nd(4),volsum,vol,totvolsum,tempsum,temp,m
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)
LOGICAL :: lhalf

CHARACTER(256) :: fn

WRITE(fn,'(A)') TRIM(caseid)
IF( DeplCntl%nburnupstep .GT. 1 )THEN
    WRITE(fn,'(A,A)') TRIM(fn),'_d'
    IF( DeplCntl%nowstep .LT. 10 )THEN
        WRITE(fn,'(A,A,i1)') TRIM(fn),'0',DeplCntl%nowstep
    ELSE
        WRITE(fn,'(A,i2)') TRIM(fn),DeplCntl%nowstep
    ENDIF
ENDIF       

WRITE(fn,'(A,A)') TRIM(fn),'.cool'

CALL openfile(io16,.FALSE.,.FALSE.,.FALSE., fn)

nz = Core%nz; nxya = Core%nxya
Pin => Core%Pin; CellInfo => Core%CellInfo
AsyInfo => Core%AsyInfo; Asy => Core%Asy

DO iz = 1, nz
  DO iasy = 1, nxya
    iasytype = Asy(iasy)%AsyType
    nx = AsyInfo(iasytype)%nx; ny = AsyInfo(iasytype)%ny
    if (mod(nx,2).eq.0) then
        hnx = nx/2
        hny = ny/2
        lhalf=.false.
    else
        hnx = nx/2 + 1
        hny = ny/2 + 1
        lhalf=.true.
    endif
    tempsum = 0._8; totvolsum = 0._8
    DO ix = hnx, nx
      DO iy = ix, ny
        ixy = Core%AsyInfo(iasyType)%pin2DIdx(ix, iy)  
        mypin => Pin(ixy); icel = mypin%Cell(iz)
        IF(CellInfo(icel)%lGap) CYCLE
        nfxr=CellInfo(icel)%nFXR
        nd = 0._8; volsum = 0._8
        m = 1._8
        if (lhalf) then
          if ((iy.eq.hny).and.(ix.eq.hnx)) then
            m = 0.25_8
          elseif ((ix.eq.iy).or.(ix.eq.hnx)) then
            m = 0.5_8
          endif
        else
          if (ix.eq.iy) m = 0.5_8
        endif
        DO ir = 1, nfxr
          ifxr = mypin%FxrIdxst + ir - 1; myFxr => FmInfo%Fxr(ifxr, iz)
          if (.not.myFxr%lh2o) CYCLE
          temp = myFxr%temp; 
          vol = myFxr%area * m
          niso = myFxr%niso; idiso => myFxr%idiso; pnum => myFxr%pnum
          do iso = 1, niso
            nd(iso) = nd(iso) + pnum(iso)*vol
          enddo
          volsum = volsum + vol
          totvolsum = totvolsum + vol
          tempsum = tempsum + temp * vol
        ENDDO
        nd = nd / volsum   
        write(io16,*)
        if (iy-hny.lt.hny.and.ix-hnx.lt.hnx) then
          write(io16,'(a3,a1,i1,a1,i1,5X,ES16.6)') '  W','0',iy-hny,'0',ix-hnx,sum(nd)
        elseif (iy-hny.lt.hny.and.ix-hnx.ge.hnx) then
          write(io16,'(a3,a1,i1,i2,5X,ES16.6)') '  W','0',iy-hny,ix-hnx,sum(nd)
        elseif (iy-hny.ge.hny.and.ix-hnx.lt.hnx) then
          write(io16,'(a3,i2,a1,i1,5X,ES16.6)') '  W',iy-hny,'0',ix-hnx,sum(nd)
        elseif (iy-hny.ge.hny.and.ix-hnx.ge.hnx) then
          write(io16,'(a3,i2,i2,5X,ES16.6)') '  W',iy-hny,ix-hnx,sum(nd)
        endif
        do iso = 1, niso
          if (idiso(iso).eq.1001) then
            write(io16,'(2i7,ES14.6,a)') 0,idiso(iso),nd(iso),' TEMPDEPT + hh2o'
          else
            write(io16,'(2i7,ES14.6,a)') 0,idiso(iso),nd(iso),' TEMPDEPT'  
          endif
        enddo
      ENDDO 
    ENDDO
    tempsum = tempsum / totvolsum
    write(io16,*)
    write(io16,'(2i3,ES16.6)') iz,iasy,tempsum
    write(io16,*)
  ENDDO  
ENDDO
close(io16)

NULLIFY(mypin)
NULLIFY(idiso)
NULLIFY(pnum)
NULLIFY(myFxr)
NULLIFY(Pin)
NULLIFY(CellInfo)
NULLIFY(AsyInfo)
NULLIFY(Asy)

END SUBROUTINE

SUBROUTINE PrintOthersTemp(Core, FmInfo)
USE TYPEDEF, ONLY : CoreInfo_Type, FmInfo_Type, AsyInfo_Type, Asy_Type, Cell_Type, Pin_Type, FxrInfo_Type
USE FILES,   ONLY : io16, caseid
USE ioutil,  ONLY : openfile
USE Depl_Mod,  ONLY : DeplCntl
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: PIN(:), mypin
TYPE(FxrInfo_Type), POINTER :: myFxr

INTEGER :: iz,nz,iasy,nxya,nx,ny,hnx,hny,ix,iy,ixy,nfxr,icel,ir,ifxr,iasytype,m,tempi
REAL :: temp
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)

CHARACTER(256) :: fn

WRITE(fn,'(A)') TRIM(caseid)
IF( DeplCntl%nburnupstep .GT. 1 )THEN
    WRITE(fn,'(A,A)') TRIM(fn),'_d'
    IF( DeplCntl%nowstep .LT. 10 )THEN
        WRITE(fn,'(A,A,i1)') TRIM(fn),'0',DeplCntl%nowstep
    ELSE
        WRITE(fn,'(A,i2)') TRIM(fn),DeplCntl%nowstep
    ENDIF
ENDIF       

WRITE(fn,'(A,A)') TRIM(fn),'.others'
CALL openfile(io16,.FALSE.,.FALSE.,.FALSE., fn)

nz = Core%nz; nxya = Core%nxya
Pin => Core%Pin; CellInfo => Core%CellInfo
AsyInfo => Core%AsyInfo; Asy => Core%Asy

DO iz = 1, nz
  DO iasy = 1, nxya
    iasytype = Asy(iasy)%AsyType
    nx = AsyInfo(iasytype)%nx; ny = AsyInfo(iasytype)%ny
    if (mod(nx,2).eq.0) then
        hnx = nx/2
        hny = ny/2
    else
        hnx = nx/2 + 1
        hny = ny/2 + 1
    endif
    DO ix = hnx, nx
      DO iy = ix, ny
        ixy = Core%AsyInfo(iasyType)%pin2DIdx(ix, iy)  
        mypin => Pin(ixy); icel = mypin%Cell(iz)
        if (.NOT.CellInfo(icel)%lfuel) CYCLE
        nfxr=CellInfo(icel)%nFXR
        DO ir = nfxr, 1, -1
          ifxr = mypin%FxrIdxst + ir - 1; myFxr => FmInfo%Fxr(ifxr, iz)
          if (myFxr%lh2o) CYCLE
          temp = myFxr%temp
          m = mod(nint(temp),5)
          if (m.gt.2) then
            tempi = nint(temp)+(5-m)  
          else
            tempi = nint(temp)-m
          endif
          if (iy-hny.lt.hny.and.ix-hnx.lt.hnx) then
            write(io16,'(a1,i2,a1,a1,i1,a1,i1,i6)') 'T',nfxr-ir+1,'_','0',iy-hny,'0',ix-hnx,tempi
          elseif (iy-hny.lt.hny.and.ix-hnx.ge.hnx) then
            write(io16,'(a1,i2,a1,a1,i1,i2,i6)') 'T',nfxr-ir+1,'_','0',iy-hny,ix-hnx,tempi
          elseif (iy-hny.ge.hny.and.ix-hnx.lt.hnx) then
            write(io16,'(a1,i2,a1,i2,a1,i1,i6)') 'T',nfxr-ir+1,'_',iy-hny,'0',ix-hnx,tempi
          elseif (iy-hny.ge.hny.and.ix-hnx.ge.hnx) then
            write(io16,'(a1,i2,a1,i2,i2,i6)') 'T',nfxr-ir+1,'_',iy-hny,ix-hnx,tempi
          endif          
        ENDDO
      ENDDO 
    ENDDO
  ENDDO  
ENDDO
close(io16)

NULLIFY(mypin)
NULLIFY(myFxr)
NULLIFY(Pin)
NULLIFY(CellInfo)
NULLIFY(AsyInfo)
NULLIFY(Asy)


END SUBROUTINE