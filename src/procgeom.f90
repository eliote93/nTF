SUBROUTINE ReadNInit_Cell(Dataline0, lCAD, lGAP, lCellInit)
use param
use typedef
use allocs
use cntl,         only : nTracerCntl
use geom,         only :  CellInfo,  CellPitch,    nCellType
use ioutil,       only :  toupper,   IFnumeric,    nfields,   fndchara,  fndchara512, &
                          nfieldto, nfieldto512
use GeomTreatment, ONLY : AnnularFsrArea, AnnularRegionDivision
use BenchXs,      only :  MacXsBen, DynMacXsBen
use Material_Mod, ONLY : Mixture
use SPH_mod,      ONLY : calcCellSSPH,calcAICCellSSPH
use XSLIB_MOD,    ONLY : igresb, igrese,nofghel
implicit none
character(512),intent(in) :: dataline0
LOGICAL :: LCAD, lGAP, lCellInit

INTEGER :: ndiv(0:300) = 0, ndivx(0:100) = 0, ndivy(0:100) = 0, ipos(100) = 0, iReg(0:300) = 0
INTEGER :: icel,item,imat(100),nmat,nreg,idiv(100)
INTEGER :: ndatafield,ndata,ndatax,ndatay
INTEGER :: nspt, nfueldiv
INTEGER :: i,j, k, m, n, iRegFr,iRegTo
INTEGER :: irfr, ir, iv,ist, ix, iy
INTEGER :: nFsrInFxr
INTEGER :: CCentIdx(2, 4)
INTEGER :: ibFuel, ieFuel
character(512) :: dataline,chline, CentPos

REAL :: rr(300), subrr(100), subvol(100), irad(0:100)
REAL :: delx(100), dely(100)
REAL :: Pitch,HalfPitch
REAL :: tvol,dvol,vol,theta,vol_

DATA CCentIdx / 1, 1,  1, 2,  2, 2,  2, 1 /

IF(.NOT. lCellInit) then
  allocate(CellInfo(nCellType))
  DO i = 1, nCellType
    CellInfo(i)%lempty = TRUE
    CellInfo(i)%luse = FALSE
  ENDDO
  lCellInit = TRUE
ENDIF
IF(LCAD) THEN
  CALL ReadNInit_CadCell(Dataline0)
  RETURN
ENDIF
IF(LGAP) THEN
  CALL ReadNInit_UsrGapCell(Dataline0)
  RETURN
ENDIF


Pitch = CellPitch;  HalfPitch = 0.5_8*CellPitch

dataline = dataline0
READ(dataline,*) icel
! icel = CellTypeMap(icel)   !--- CNJ Edit : Flexible Cell & Pin Numbering
CellInfo(icel)%icel0 = icel
CellInfo(icel)%lempty = FALSE
CellInfo(icel)%lrect = FALSE
CellInfo(icel)%lCCell = .FALSE.
CellInfo(icel)%geom%lCCent  = FALSE
CellInfo(icel)%lFuel = FALSE
CellInfo(icel)%lRes = FALSE
CellInfo(icel)%lgap = FALSE
CellInfo(icel)%lGd = FALSE
CellInfo(icel)%lMox = .FALSE.
CellInfo(icel)%lAIC = .FALSE.
CellInfo(icel)%lCentX = FALSE; CellInfo(icel)%lCentY = FALSE; CellInfo(icel)%lCentXY = FALSE; 
nDataField = len_trim(dataline)


CALL fndchara512(dataline,ipos,nspt,SLASH)

IF(nspt .eq. 2 .or. nspt .eq. 3) THEN
  CellInfo(icel)%nDivAzi = 8
  ndata = nfieldto512(dataline,SLASH)-1
  
  CellInfo(icel)%Geom%lCircle = TRUE
  CellInfo(icel)%Geom%lRect = FALSE
  CellInfo(icel)%Geom%nCircle = 0
  CellInfo(icel)%Geom%nx = 3
  CellInfo(icel)%Geom%ny = 3
  
  read(dataline(ipos(2)+1:nDataField),*) (nDiv(i),i=1,ndata)
  !Circle Division Read
  DO i =1, ndata
    CellInfo(icel)%geom%nCircle = CellInfo(icel)%geom%nCircle + nDiv(i)
  ENDDO
  CellInfo(icel)%nFXR = CellInfo(icel)%geom%nCircle + 1
  !Off-Center Pincell
  CellInfo(icel)%geom%cx = 0 
  CellInfo(icel)%geom%cy = 0
  IF(nspt .eq. 3) then
    CellInfo(icel)%geom%lCCent  = TRUE
    CellInfo(icel)%lCCell = TRUE
    CellInfo(icel)%nDivAzi = 2
    read(dataline(ipos(3)+1:nDataField),*) CentPos
    IF(CentPos .EQ. 'SW')  CellInfo(icel)%geom%CCentType = 1
    IF(CentPos .EQ. 'NW')  CellInfo(icel)%geom%CCentType = 2
    IF(CentPos .EQ. 'NE')  CellInfo(icel)%geom%CCentType = 3
    IF(CentPos .EQ. 'SE')  CellInfo(icel)%geom%CCentType = 4
    !IF( CellInfo(icel)%geom%CCentType) CALL TERMINATE('Cell Input Error')
  ENDIF
  CellInfo(icel)%nFSR = (CellInfo(icel)%geom%nCircle + 1)*CellInfo(icel)%nDivAzi
  !Composition Read
  iRegFr=0
  iRegTo=nData
  read(dataline(ipos(1)+1:ipos(2)-1),*) (iReg(i),i=iRegTo,iRegFr,-1)
  !read(dataline(ipos(3)+1:nDataField),*) (iReg(i),i=1,ndata)
ELSE
  CellInfo(icel)%lRect = TRUE
  CellInfo(icel)%geom%lrect = TRUE
  CellInfo(icel)%geom%lcircle = FALSE
  CellInfo(icel)%Geom%nCircle =0
  CellInfo(icel)%geom%cx = 0 
  CellInfo(icel)%geom%cy = 0
  ndatax = nfieldto(dataline,SLASH)
  read(dataline, *) item, (delx(i), i = 1, ndatax-1)
  delx(ndatax) = CellPitch - sum(delx(1:ndatax-1))
  
  !---BYS edit for GC GeomReStart
  CellInfo(icel)%geom%inpnx=ndatax
  CellInfo(icel)%geom%inpdelx(1:ndatax-1)=delx(1:ndatax-1)
  !---BYS edit for GC GeomReStart END
  read(dataline(ipos(3)+1:ipos(4)-1),*) (ndivx(i), i=1,ndatax)
  !---BYS edit for GC GeomReStart
  CellInfo(icel)%geom%inpdivx(1:ndatax)=ndivx(1:ndatax)
  !---BYS edit for GC GeomReStart END
  CellInfo(icel)%geom%nx =1
  DO i = 1, ndatax
    CellInfo(icel)%geom%nx = CellInfo(icel)%geom%nx + ndivx(i)
    delx(i) = delx(i) / ndivx(i)
  ENDDO
  chline = dataline(ipos(1)+1:256)
  ndatay = nfieldto(chline,SLASH) + 1
  read(dataline(ipos(1)+1:ipos(2)-1), *) (dely(i), i = 1, ndatay-1)
  !IF(ndatay .EQ. 1) dely(1) = CellPitch
  dely(ndatay) = CellPitch - sum(dely(1:ndatay-1))
  !---BYS edit for GC GeomReStart
  CellInfo(icel)%geom%inpny=ndatay
  CellInfo(icel)%geom%inpdely(1:ndatay-1)=dely(1:ndatay-1)
  !---BYS edit for GC GeomReStart END
  read(dataline(ipos(4)+1:256),*) (ndivy(i), i=1,ndatay)
    !---BYS edit for GC GeomReStart
  CellInfo(icel)%geom%inpdivy(1:ndatay)=ndivy(1:ndatay)
  !---BYS edit for GC GeomReStart END
  CellInfo(icel)%geom%ny =1
  DO i = 1, ndatay
    CellInfo(icel)%geom%ny = CellInfo(icel)%geom%ny + ndivy(i)
    dely(i) = dely(i) / ndivy(i)
  ENDDO
  iRegFr = 1
  iRegTo = nDataX*nDataY
  read(dataline(ipos(2)+1:ipos(3)-1),*) (ireg(i),i=iRegFr,iRegTo)
  
  CellInfo(icel)%nFSR = (CellInfo(icel)%geom%nx-1)*(CellInfo(icel)%geom%ny-1)
  !CellInfo(icel)%nFXR = 1
  CellInfo(icel)%nFXR = nDataX * nDataY
  CellInfo(icel)%geom%nx=CellInfo(icel)%geom%nx-1
  CellInfo(icel)%geom%ny=CellInfo(icel)%geom%ny-1
ENDIF
!Composition Check
CellInfo(icel)%Geom%lx = CellPitch; CellInfo(icel)%geom%ly = CellPitch
CellInfo(icel)%Geom%x(1) = -CellPitch*Half; CellInfo(icel)%Geom%x(2) = CellPitch*Half
CellInfo(icel)%Geom%y(1) = -CellPitch*Half; CellInfo(icel)%Geom%y(2) = CellPitch*Half
!CALL dmalloc(CellInfo(icell)%FxrIdxSt,CellInfo(icell)%nFxr)
!CALL dmalloc(CellInfo(icell)%FxrIdxSt,CellInfo(icell)%nFxr)
!
IF(.not. cellInfo(icel)%lRect) then
  IF(CellInfo(icel)%Geom%lCCent) THEN
    HalfPitch = HalfPitch  * 2._8
    !Cellinfo(icel)%geom%cx = CellInfo(icel)%Geom%x(CCentIdx(1, CellInfo(icel)%Geom%CCentType))
    !Cellinfo(icel)%geom%cy = CellInfo(icel)%Geom%y(CCentIdx(2, CellInfo(icel)%Geom%CCentType))    
  ENDIF
  read(dataline,*) item, (rr(j),j=ndata,1,-1)    !Read Annunal ring Info
  read(dataline(ipos(2)+1:ipos(3)-1),*) (ndiv(j),j=ndata,1,-1)
  
  !lfuel and lRes Setting
  DO j=0, ndata
    IF(nTracerCntl%lXsLib) THEN
      CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. Mixture(ireg(j))%lfuel
      CellInfo(icel)%lres = CellInfo(icel)%lres .or. Mixture(ireg(j))%lres
      CellInfo(icel)%lGd = CellInfo(icel)%lGd .or. Mixture(ireg(j))%lGd
      CellInfo(icel)%lMOX = CellInfo(icel)%lMOX .or. Mixture(ireg(j))%lMOX
      CellInfo(icel)%lAIC = CellInfo(icel)%lAIC .or. Mixture(ireg(j))%lAIC
      CellInfo(icel)%lsSPH = CellInfo(icel)%lfuel .or. CellInfo(icel)%lAIC
    ELSE
      IF(nTracerCntl%lDynamicBen) THEN
        CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. DynMacXsBen(ireg(j))%lfuel
        CellInfo(icel)%lCR = CellInfo(icel)%lCR .or. DynMacXsBen(ireg(j))%lCR
      ELSE
        CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. MacXsBen(ireg(j))%lfuel
        CellInfo(icel)%lCR = CellInfo(icel)%lCR .or. MacXsBen(ireg(j))%lCR
      END IF
    ENDIF
  ENDDO
  
  imat=0; irad=0._8; idiv=0
  imat(1)=iReg(ndata); nmat=1
  irad(1)=rr(ndata)
  idiv(1)=ndiv(ndata)
  do i=iRegTo-1,iRegFr,-1
    if (iReg(i).eq.iReg(i+1)) then
        if (i.eq.0) then
            idiv(nmat)=idiv(nmat)+1
            irad(nmat)=pitch/dsqrt(PI)
        else
            idiv(nmat)=idiv(nmat)+ndiv(i)
        endif
        cycle
    endif
    nmat=nmat+1
    imat(nmat)=iReg(i)
    if (i.gt.0) then
        irad(nmat)=rr(i)
        idiv(nmat)=ndiv(i)
    else
        irad(nmat)=pitch/dsqrt(PI)
        idiv(nmat)=1
    endif
  enddo
  allocate(CellInfo(icel)%matidx(nmat))
  CellInfo(icel)%nmat=nmat
  CellInfo(icel)%matidx(1:nmat)=imat(1:nmat)
  allocate(CellInfo(icel)%matrad(nmat))
  CellInfo(icel)%matrad(1:nmat)=irad(1:nmat)
  
  nreg = sum(idiv(1:nmat))
  CellInfo(icel)%nreg_cp=nreg
  allocate(CellInfo(icel)%rad_cp(nreg),CellInfo(icel)%fxrvol(nreg))
  CellInfo(icel)%rad_cp=0._8
  vol=0._8
  k=1
  do i = 1, nmat
    vol_ = PI*(irad(i)*irad(i)-irad(i-1)*irad(i-1))/idiv(i)  
    do j = 1, idiv(i)    
      vol = vol + vol_
      CellInfo(icel)%rad_cp(k) = dsqrt(vol*INVPI)  
      CellInfo(icel)%fxrvol(k) = vol_
      k=k+1
    enddo
  enddo
  
  if (CellInfo(icel)%lfuel) then
    do i = nmat, 1, -1
      IF(nTracerCntl%lXsLib) THEN
        IF (Mixture(imat(i))%lfuel) EXIT 
      ENDIF
    enddo
    do k = i+1, nmat
      IF(nTracerCntl%lXsLib) THEN
        IF (Mixture(imat(k))%lh2o) EXIT 
      ENDIF
    enddo
    k=sum(idiv(1:k-1));
    CellInfo(icel)%fuelgapcldvol=sum(CellInfo(icel)%fxrvol(1:k))
    CellInfo(icel)%cldfxridx=k
    CellInfo(icel)%nmodfxr=nreg-k
    CellInfo(icel)%invnmodfxr=1._8/real(CellInfo(icel)%nmodfxr,8)
  elseif (CellInfo(icel)%lAIC) then
    do i = nmat, 1, -1
      IF(nTracerCntl%lXsLib) THEN
        IF (Mixture(imat(i))%lAIC) EXIT 
      ENDIF
    enddo    
    do k = i+1, nmat
      IF(nTracerCntl%lXsLib) THEN
        IF (Mixture(imat(k))%lh2o) EXIT 
      ENDIF
    enddo
    k=sum(idiv(1:k-1));
    CellInfo(icel)%fuelgapcldvol=sum(CellInfo(icel)%fxrvol(1:k))
    CellInfo(icel)%cldfxridx=k
    CellInfo(icel)%nmodfxr=nreg-k
    CellInfo(icel)%invnmodfxr=1._8/real(CellInfo(icel)%nmodfxr,8)   
  endif
  
  ! Calculate volume of each MATERIAL region (not each fxr)
  ALLOCATE(CellInfo(icel)%vol(0:ndata))
  vol=0._8
  DO j=ndata,1,-1
      IF (rr(j).gt.HalfPitch) THEN
          theta=acos(HalfPitch/rr(j))
          CellInfo(icel)%vol(j)=4._8*rr(j)*(HalfPitch*sin(theta)+rr(j)*(PI*0.25_8-theta))-vol
      ELSE
          CellInfo(icel)%vol(j)=PI*rr(j)*rr(j)-vol
      ENDIF
      vol=vol+CellInfo(icel)%vol(j)
  ENDDO
  CellInfo(icel)%vol(0)=Pitch*Pitch-vol  
  IF (CellInfo(icel)%lfuel) THEN
    IF (nTracerCntl%lXsLib) THEN
      CellInfo(icel)%lhole = .NOT. Mixture(ireg(ndata))%lfuel
      DO j = ndata, 0, -1
        IF (Mixture(ireg(j))%lfuel) EXIT
      END DO
    ELSEIF(nTracerCntl%lDynamicBen) THEN 
      CellInfo(icel)%lhole = .NOT. DynMacXsBen(ireg(ndata))%lfuel
      DO j = ndata, 0, -1
        IF (DynMacXsBen(ireg(j))%lfuel) EXIT
      END DO
    ELSE
      CellInfo(icel)%lhole = .NOT. MacXsBen(ireg(ndata))%lfuel
      DO j = ndata, 0, -1
        IF (MacXsBen(ireg(j))%lfuel) EXIT
      END DO
    END IF
    CellInfo(icel)%ibFuel = j
    DO j=ndata-1,0,-1
      IF(nTracerCntl%lXsLib) THEN
        IF (Mixture(ireg(j+1))%lfuel.and..not.Mixture(ireg(j))%lfuel) EXIT
      ELSE
        IF(nTracerCntl%lDynamicBen) THEN
          IF (DynMacXsBen(ireg(j+1))%lfuel .AND. .not. DynMacXsBen(ireg(j))%lfuel) EXIT
        ELSE
          IF (MacXsBen(ireg(j+1))%lfuel .AND. .not.MacXsBen(ireg(j))%lfuel) EXIT
        END IF
      ENDIF
    ENDDO
    j=j+1
    CellInfo(icel)%ieFuel = j
    CellInfo(icel)%FuelRad0 = rr(CellInfo(icel)%ieFuel)
  ELSEIF (CellInfo(icel)%lAIC) THEN
    DO j=ndata,0,-1
        IF (.not.Mixture(ireg(j))%lAIC) EXIT ! Based on the assumption that there's no hole in an AIC pin.
    ENDDO
    j=j+1
    CellInfo(icel)%ieFuel = j
    CellInfo(icel)%ibFuel = ndata
    CellInfo(icel)%FuelRad0 = rr(CellInfo(icel)%ieFuel)
  ENDIF
  ieFuel = CellInfo(icel)%ieFuel
  ibFuel = CellInfo(icel)%ibFuel
  
  nfueldiv=0
  DO k=ibFuel,ieFuel,-1
      nfueldiv=nfueldiv+ndiv(k)
  ENDDO
  CellInfo(icel)%nfueldiv=nfueldiv
  
  !! SPH factor
  IF (nTracerCntl%lSSPH) THEN
      IF (CellInfo(icel)%lAIC) THEN
          call calcAICCellSSPH(CellInfo,icel,ndata,iReg,rr,ndiv,j,nfueldiv)
      ELSE
          call calcCellSSPH(CellInfo,icel,ndata,iReg,rr,ndiv,ibFuel,ieFuel,nfueldiv,CellInfo(icel)%lhole)
      ENDIF
  ENDIF
  
  !Allocation
  CALL dmalloc(CellInfo(icel)%iReg, CellInfo(icel)%nFSR)
  CALL dmalloc(CellInfo(icel)%FxrIdxSt,CellInfo(icel)%nFxr)
  CALL dmalloc(CellInfo(icel)%nFsrInFxr,CellInfo(icel)%nFxr)
  CALL dmalloc(CellInfo(icel)%MapFxr2FsrIdx, CellInfo(icel)%nFSR, CellInfo(icel)%nFxr)
  nFsrInFxr = CellInfo(icel)%nDivAzi
  ndiv(0) =1 
  
  rr(ndata+1) = 0
  tvol = PI * rr(1) * rr(1)
  CellInfo(icel)%FxrIdxSt(1) = 1
  CellInfo(icel)%nFsrInFxr(1) = nFsrInFxr
  DO j = 1, nFsrInFxr
    CellInfo(icel)%ireg(j) = ireg(0)
  ENDDO
  irfr = 1
  CALL dmalloc(CellInfo(icel)%geom%Circle,3,CellInfo(icel)%geom%nCircle)
  DO j = 1, ndata
    dvol = PI*(rr(j)*rr(j)-rr(j+1)*rr(j+1))/ndiv(j)
    tvol = PI*rr(j)*rr(j)
    IF(HalfPitch .LE. rr(j)) THEN
      subrr(1:ndiv(j)) = AnnularRegionDivision(HalfPitch, rr(j:j+1), ndiv(j))
    ENDIF
    DO ir = 1, ndiv(j)
      Cellinfo(icel)%FxrIdxSt(irfr+ir) = (irfr+ir-1)*nFsrInFxr + 1   !Staring Index of local FSR index in a FXR
      Cellinfo(icel)%nFsrInFxr(irfr+ir) = nFsrinFxr                  !Number of FSRs which are beleong to a certain FXR
      IF(HalfPitch .GT. rr(j)) THEN
        Cellinfo(icel)%geom%circle(3,irfr+ir-1) = sqrt(tvol*INVPI)        !Annularing 
        tvol = tvol - dvol
      ELSE
        Cellinfo(icel)%geom%circle(3,irfr+ir-1) = Subrr(ir)
      ENDIF
      IF(CellInfo(icel)%Geom%lCCent) THEN
        Cellinfo(icel)%geom%circle(1,irfr+ir-1) = CellInfo(icel)%Geom%x(CCentIdx(1, CellInfo(icel)%Geom%CCentType))
        Cellinfo(icel)%geom%circle(2,irfr+ir-1) = CellInfo(icel)%Geom%y(CCentIdx(2, CellInfo(icel)%Geom%CCentType))
      ENDIF
      !Assign the Mixture or composition number in a FSR
      ist = (irfr+ir-1)*nFsrInFxr
      DO iv = 1, nFsrInFxr
        CellInfo(icel)%iReg(ist+iv) = ireg(j)
      ENDDO
    ENDDO
    irfr = irfr + ndiv(j)
  ENDDO
  !FXS region Index Starting Point-> FSR Region Idx  
  DO ir = 1, Cellinfo(icel)%nFxr
    ist = CellInfo(icel)%FxrIdxSt(ir)
    DO j = 1, Cellinfo(icel)%nFsrInFxr(ir)
      Cellinfo(icel)%MapFxr2FsrIdx(j,ir) = ist+j-1
    ENDDO
  ENDDO
ELSE  !Rectangular Geomemtry
  !Allocation
  CALL dmalloc(CellInfo(icel)%iReg, CellInfo(icel)%nFSR)
  CALL dmalloc(CellInfo(icel)%FxrIdxSt,CellInfo(icel)%nFxr)
  CALL dmalloc(CellInfo(icel)%nFsrInFxr,CellInfo(icel)%nFxr)
  CALL dmalloc(CellInfo(icel)%MapFxr2FsrIdx, CellInfo(icel)%nFSR, CellInfo(icel)%nFxr)
  DO i = iregfr,iregto
    IF(nTracerCntl%lXsLib) THEN
      CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. Mixture(ireg(i))%lfuel
    ELSE
      IF(nTracerCntl%lDynamicBen) THEN
        CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. DynMacXsBen(ireg(i))%lfuel
      ELSE
        CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. MacXsBen(ireg(i))%lfuel
      END IF
    ENDIF
  ENDDO
  i = 0 
  DO iy = 1, nDatay
    DO ix = 1, nDatax
      i = i + 1
      CellInfo(icel)%nFsrInFxr(i) = ndivx(ix) * ndivy(iy)
    ENDDO
  ENDDO
  
  CellInfo(icel)%FxrIdxSt(1) = 1
  
  DO i = 2, nDatay * nDatax
    CellInfo(icel)%FxrIdxSt(i) = CellInfo(icel)%FxrIdxSt(i-1) + CellInfo(icel)%nFsrInFxr(i-1)
  ENDDO
  
  DO i = 2, 100
    ndivy(i) = ndivy(i-1) + ndivy(i)
    ndivx(i) = ndivx(i-1) + ndivx(i)
  ENDDO

  CALL Dmalloc(CellInfo(icel)%Geom%delx, CellInfo(icel)%Geom%nx)
  CALL Dmalloc(CellInfo(icel)%Geom%dely, CellInfo(icel)%Geom%ny)

  DO ix = 1, ndatax
    DO j= ndivx(ix-1) + 1, ndivx(ix)
      CellInfo(icel)%Geom%DelX(j) = DelX(ix) 
    ENDDO
  ENDDO 

  DO iy = 1, ndatay
    DO j= ndivy(iy-1) + 1, ndivy(iy)
      CellInfo(icel)%Geom%Dely(j) = Dely(iy) 
    ENDDO
  ENDDO

  m = 0; n = 0
  DO iy = 1, ndatay
    DO ix = 1, ndatax
      m = m + 1; n = 0
      DO j = ndivy(iy-1) + 1, ndivy(iy)
        DO i = ndivx(ix-1) + 1, ndivx(ix)
          n = n + 1
          k = CellInfo(icel)%Geom%nx * (j - 1) + i
          CellInfo(icel)%MapFxr2FsrIdx(n, m) = k
          CellInfo(icel)%iReg(k) = ireg(m)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  !CellInfo(icel)%FxrIdxSt(1) = 1
  !CellInfo(icel)%nFsrInFxr(1) = CellInfo(icel)%nFSR
  
  !DO i = 1, CellInfo(icel)%nFSR
  !  CellInfo(icel)%MapFxr2FsrIdx(i,1) = i
  !  CellInfo(icel)%iReg(i) = iReg(1)
  !ENDDO
  !TH setting
  
ENDIF
CellInfo(icel)%CellInput = Dataline
END SUBROUTINE

SUBROUTINE Init_GapCell(DataLineIn, lCellInit)
!Initialize Gap Cell
use param
use typedef
use allocs
use cntl
use geom,   only :  lGap, CellInfo,  CellPitch,     AsyPitch,  nCellType0,  nCellType,   &
                    nCellX0
!use ioutil, only :  toupper,   IFnumeric,    nfields,   fndchara,    nfieldto
use BenchXs,only :  MacXsBen
implicit none
character(256),intent(in) :: datalineIn
LOGICAL :: lCellInit

INTEGER :: icel, icel1, icel2, icel3
INTEGER :: i,j
INTEGER :: GapComp
REAL :: Pitch,HalfPitch,GapPitch
REAL  :: tvol,dvol
logical,save :: first
data first /.TRUE./

if(.NOT. lgap) return

IF(.NOT. lCellInit) then
  allocate(CellInfo(nCellType))
  DO i = 1, nCellType
    CellInfo(i)%lempty = TRUE
    CellInfo(i)%luse = FALSE
  ENDDO
  lCellInit = TRUE
ENDIF

read(DataLineIn, *) lGap, GapComp

GapPitch = (AsyPitch - CellPitch*nCellX0)*HALF

icel1 = nCellType0 - 2;  icel2 = nCellType0 - 1;  icel3 = nCellType0
DO icel = icel1, icel3
  Cellinfo(icel)%lfuel = FALSE; CellInfo(icel)%lempty = FALSE;
  CellInfo(icel)%lCentX = FALSE; CellInfo(icel)%lCentY = FALSE;
  CellInfo(icel)%lCentXY = FALSE; CellInfo(icel)%lgap = TRUE; 
  CellInfo(icel)%Geom%Lcircle = FALSE
  CellInfo(icel)%Geom%nCircle =0
  CellInfo(icel)%Geom%lCCent = FALSE
  CellInfo(icel)%lRect = TRUE
  CellInfo(icel)%geom%cx = 0 
  CellInfo(icel)%geom%cy = 0
ENDDO

DO icel = icel1, icel3
  CALL Dmalloc(CellInfo(icel)%FxrIdxSt, 1) 
  CALL Dmalloc(CellInfo(icel)%nFsrInFxr, 1)
  CALL Dmalloc(CellInfo(icel)%MapFxr2FsrIdx,8,1)
  CALL Dmalloc(CellInfo(icel)%iReg, 8)       !Mixture or Composition Idx
ENDDO
!Up and Bottom Gap Cell
icel = icel1
CellInfo(icel)%geom%lx = CellPitch 
CellInfo(icel)%Geom%ly = GapPitch
CellInfo(icel)%Geom%nx = 4                     !Region 
CellInfo(icel)%Geom%ny = 2
CellInfo(icel)%GapType = 1

!Left and Right Gap Cell
icel = icel2
CellInfo(icel)%geom%ly = CellPitch
CellInfo(icel)%Geom%lx = GapPitch
CellInfo(icel)%Geom%nx = 2
CellInfo(icel)%Geom%ny = 4
CellInfo(icel)%GapType = 2
!Coner Gap Cell
icel = icel3
CellInfo(icel)%Geom%lx = GapPitch
CellInfo(icel)%Geom%ly = GapPitch
CellInfo(icel)%Geom%nx = 2
CellInfo(icel)%Geom%ny = 2
CellInfo(icel)%GapType = 3
DO icel = icel1, icel3
  
  CALL Dmalloc(CellInfo(icel)%Geom%delx, CellInfo(icel)%Geom%nx)
  CALL Dmalloc(CellInfo(icel)%Geom%dely, CellInfo(icel)%Geom%ny)
  CellInfo(icel)%Geom%lcircle = FALSE
  CellInfo(icel)%Geom%nCircle = 0
  CellInfo(icel)%Geom%lrect = TRUE
  CellInfo(icel)%lfuel = .FALSE.
  !CellInfo(icel)%ireg(1) = GapComp
  CellInfo(icel)%nFXR = 1
  CellInfo(icel)%nFSR = CellInfo(icel)%geom%nx*CellInfo(icel)%geom%ny
  CellInfo(icel)%nFsrInFxr(1) = CellInfo(icel)%nFSR
  CellInfo(icel)%FxrIdxSt(1) = 1
  DO i = 1, CellInfo(icel)%nFSR
    CellInfo(icel)%ireg(i) = GapComp
    CellInfo(icel)%MapFxr2FsrIdx(i, 1) = i
  ENDDO
  CellInfo(icel)%Geom%x(1) = -CellInfo(icel)%Geom%lx * Half
  CellInfo(icel)%Geom%x(2) = CellInfo(icel)%Geom%lx * Half
  CellInfo(icel)%Geom%y(1) = -CellInfo(icel)%Geom%ly*Half
  CellInfo(icel)%Geom%y(2) = CellInfo(icel)%Geom%ly*Half
  
  CALL Dmalloc(CellInfo(icel)%Geom%delx, CellInfo(icel)%Geom%nx)
  CALL Dmalloc(CellInfo(icel)%Geom%dely, CellInfo(icel)%Geom%ny)
  
  CellInfo(icel)%Geom%DelX(:) = CellInfo(icel)%Geom%lx / CellInfo(icel)%Geom%nx
  CellInfo(icel)%Geom%DelY(:) = CellInfo(icel)%Geom%ly / CellInfo(icel)%Geom%ny
  !CellInfo(icel)%Geom%DelX()
  CONTINUE
  CellInfo(icel)%ireg(:) = GapComp
  continue
ENDDO

END SUBROUTINE

SUBROUTINE ReadNInit_UsrGapCell(Dataline0)
use param
use typedef
use allocs
use cntl,         only : nTracerCntl
use geom,         only :  CellInfo,  CellPitch,    AsyPitch,  nCellX0,     nCellType,lGap,          &
                          !--- CNJ Edit : Flexible Cell & Pin Numbering
                          CellTypeMap,  PinTypeMap
use ioutil,       only :  toupper,   IFnumeric,    nfields,   fndchara,    nfieldto, terminate
use GeomTreatment, ONLY : AnnularFsrArea, AnnularRegionDivision
use BenchXs,      only :  MacXsBen
use Material_Mod, ONLY : Mixture
implicit none
character(256),intent(in) :: dataline0

INTEGER :: ndiv(0:100), ndivx(0:100), ndivy(0:100), ipos(100),iReg(0:100)
INTEGER :: icel,item, igaptype
INTEGER :: ndatafield,ndata,ndatax,ndatay
INTEGER :: nspt
INTEGER :: i,j, k, m, n, iRegFr,iRegTo
INTEGER :: irfr, ir, iv,ist, ix, iy
INTEGER :: nFsrInFxr
character(256) :: dataline,chline, GapType

REAL :: GapPitch0(2,3)
REAL :: lx, ly, delx(100), dely(100)
REAL :: Pitch,HalfPitch, GapPitch
REAL :: tvol,dvol

logical,save :: first
data first /.TRUE./


IF(.NOT. lGap) RETURN

GapPitch = (AsyPitch - CellPitch*nCellX0)*HALF
Pitch = CellPitch;  HalfPitch = 0.5_8*CellPitch
GapPitch0 = RESHAPE((/CellPitch, GapPitch,  GapPitch, CellPitch, GapPitch, GapPitch/), (/2,3/))

dataline = dataline0
READ(dataline,*) icel
! icel = CellTypeMap(icel)   !--- CNJ Edit : Flexible Cell & Pin Numbering
CellInfo(icel)%lempty = FALSE
CellInfo(icel)%lrect = FALSE
CellInfo(icel)%lCCell = .FALSE.
CellInfo(icel)%geom%lCCent  = FALSE
CellInfo(icel)%lFuel = FALSE
CellInfo(icel)%lRes = FALSE
CellInfo(icel)%lgap = TRUE
CellInfo(icel)%lGd = FALSE
CellInfo(icel)%lMox = .FALSE.
CellInfo(icel)%lCentX = FALSE; CellInfo(icel)%lCentY = FALSE; CellInfo(icel)%lCentXY = FALSE; 
nDataField = len_trim(dataline)


CALL fndchara(dataline,ipos,nspt,SLASH)

read(dataline(ipos(5)+1:256),*) GapType
CALL TOUPPER(GapType); iGapType=0
IF(GapType .EQ. 'TB') iGapType=1
IF(GapType .EQ. 'LR') iGapType=2
IF(GapType .EQ. 'CN') iGapType=3
IF(iGapType .EQ. 0) THEN
  CALL terminate('Error in GapCell Card: Not Proper Gap Cell Type')
ENDIF

lx = GapPitch0(1, iGapType); ly = GapPitch0(2, iGapType)

CellInfo(icel)%lRect = TRUE
CellInfo(icel)%geom%lrect = TRUE
CellInfo(icel)%geom%lcircle = FALSE
CellInfo(icel)%geom%nCircle = 0
CellInfo(icel)%geom%cx = 0 
CellInfo(icel)%geom%cy = 0
ndatax = nfieldto(dataline,SLASH)
read(dataline, *) item, (delx(i), i = 1, ndatax-1)
delx(ndatax) = lx - sum(delx(1:ndatax-1))
  
ndivx = 0; ndivy = 0
read(dataline(ipos(3)+1:ipos(4)-1),*) (ndivx(i), i=1,ndatax)
CellInfo(icel)%geom%nx =1
DO i = 1, ndatax
  CellInfo(icel)%geom%nx = CellInfo(icel)%geom%nx + ndivx(i)
  delx(i) = delx(i) / ndivx(i)
ENDDO
chline = dataline(ipos(1)+1:256)
ndatay = nfieldto(chline,SLASH) + 1
read(dataline(ipos(1)+1:ipos(2)-1), *) (dely(i), i = 1, ndatay-1)
!IF(ndatay .EQ. 1) dely(1) = CellPitch
dely(ndatay) = ly - sum(dely(1:ndatay-1))
read(dataline(ipos(4)+1:256),*) (ndivy(i), i=1,ndatay)
CellInfo(icel)%geom%ny =1
DO i = 1, ndatay
  CellInfo(icel)%geom%ny = CellInfo(icel)%geom%ny + ndivy(i)
  dely(i) = dely(i) / ndivy(i)
ENDDO
iRegFr = 1
iRegTo = nDataX*nDataY
read(dataline(ipos(2)+1:ipos(3)-1),*) (ireg(i),i=iRegFr,iRegTo)
 

CellInfo(icel)%GapType = iGapType
CellInfo(icel)%nFSR = (CellInfo(icel)%geom%nx-1)*(CellInfo(icel)%geom%ny-1)
!CellInfo(icel)%nFXR = 1
CellInfo(icel)%nFXR = nDataX * nDataY
CellInfo(icel)%geom%nx=CellInfo(icel)%geom%nx-1
CellInfo(icel)%geom%ny=CellInfo(icel)%geom%ny-1
!Composition Check
CellInfo(icel)%Geom%lx = lx; CellInfo(icel)%geom%ly = ly
CellInfo(icel)%Geom%x(1) = -lx*Half; CellInfo(icel)%Geom%x(2) = lx*Half
CellInfo(icel)%Geom%y(1) = -ly*Half; CellInfo(icel)%Geom%y(2) = ly*Half
!Allocation
CALL dmalloc(CellInfo(icel)%iReg, CellInfo(icel)%nFSR)
CALL dmalloc(CellInfo(icel)%FxrIdxSt,CellInfo(icel)%nFxr)
CALL dmalloc(CellInfo(icel)%nFsrInFxr,CellInfo(icel)%nFxr)
CALL dmalloc(CellInfo(icel)%MapFxr2FsrIdx, CellInfo(icel)%nFSR, CellInfo(icel)%nFxr)

DO i = iregfr,iregto
  IF(nTracerCntl%lXsLib) THEN
    CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. Mixture(ireg(i))%lfuel
  ELSE
    CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. MacXsBen(ireg(i))%lfuel
  ENDIF  
ENDDO
i = 0 
DO iy = 1, nDatay
  DO ix = 1, nDatax
    i = i + 1
    CellInfo(icel)%nFsrInFxr(i) = ndivx(ix) * ndivy(iy)
  ENDDO
ENDDO
  
CellInfo(icel)%FxrIdxSt(1) = 1
  
DO i = 2, nDatay * nDatax
  CellInfo(icel)%FxrIdxSt(i) = CellInfo(icel)%FxrIdxSt(i-1) + CellInfo(icel)%nFsrInFxr(i-1)
ENDDO

DO i = 2, 100
  ndivy(i) = ndivy(i-1) + ndivy(i)
  ndivx(i) = ndivx(i-1) + ndivx(i)
ENDDO

CALL Dmalloc(CellInfo(icel)%Geom%delx, CellInfo(icel)%Geom%nx)
CALL Dmalloc(CellInfo(icel)%Geom%dely, CellInfo(icel)%Geom%ny)

DO ix = 1, ndatax
  DO j= ndivx(ix-1) + 1, ndivx(ix)
    CellInfo(icel)%Geom%DelX(j) = DelX(ix) 
  ENDDO
ENDDO 

DO iy = 1,ndatay
  DO j= ndivy(iy-1) + 1, ndivy(iy)
    CellInfo(icel)%Geom%Dely(j) = Dely(iy) 
  ENDDO
ENDDO

m = 0; n = 0
DO iy = 1, ndatay
  DO ix = 1, ndatax
    m = m + 1; n = 0
    DO j = ndivy(iy-1) + 1, ndivy(iy)
      DO i = ndivx(ix-1) + 1, ndivx(ix)
        n = n + 1
        k = CellInfo(icel)%Geom%nx * (j - 1) + i
        CellInfo(icel)%MapFxr2FsrIdx(n, m) = k
        CellInfo(icel)%iReg(k) = ireg(m)
      ENDDO
    ENDDO
  ENDDO
ENDDO
CellInfo(icel)%CellInput = Dataline
END SUBROUTINE

SUBROUTINE ReadNinit_Pin(datalineIn, lGAP)
!Read Pin CARD Input line
use param
use typedef
use allocs
use cntl
use ioutil,         only : terminate
use geom,           only : CellInfo,      BaseCellInfo,    PinInfo,        nz,       nPinType,     nPinType0,   &
                           nCellType0,     nCellType,                                          &
                           !--- CNJ Edit : Flexible Cell & Pin Numbering
                           CellTypeMap,  PinTypeMap
IMPLICIT NONE

character(256),intent(in) :: datalineIn
LOGICAL :: lGAP
character(256) :: dataline
INTEGER :: ipin, icel, ibcel
INTEGER :: i,j,k
INTEGER :: nFSR, nFXR    !Number of Flat Source and XS regsion
logical,save :: first
data first /.TRUE./

IF(first) then
  allocate(PinInfo(nPinType))
  DO i=1, nPinType
    PinInfo(i)%lempty = .TRUE.
  ENDDO
  first = FALSE
ENDIF

IF(lGAP) THEN
  CALL ReadNinit_UsrGapPin(datalineIn)
  RETURN
ENDIF

dataline = datalineIn
read(dataline,*) ipin                          !Read the Pin Number
! ipin = PinTypeMap(ipin)   !--- CNJ Edit : Flexible Cell & Pin Numbering
CALL dmalloc(PinInfo(ipin)%cell,nz)

read(dataline,*) i , (PinInfo(ipin)%cell(k), k = 1, nz)  !Read the Pin Configuration

! !--- CNJ Edit : Flexible Cell & Pin Numbering
! DO k = 1, nz
!   PinInfo(ipin)%cell(k) = CellTypeMap(PinInfo(ipin)%cell(k))
! ENDDO

PinInfo(ipin)%ncell = nz; PinInfo(ipin)%lGap = .FALSE.
PinInfo(ipin)%lEmpty= FALSE;  PinInfo(ipin)%lUse = FALSE; PinInfo(ipin)%lFuel = FALSE
PinInfo(ipin)%lCentX = FALSE; PinInfo(ipin)%lCentY = FALSE; PinInfo(ipin)%lCentXY = FALSE 
PinInfo(ipin)%nFsrMax = 0; PinInfo(ipin)%nFxrMax = 0

DO k= 1, nz
  icel = PinInfo(ipin)%cell(k)
  !--- JSR Edit : nTIG Restart
  IF (nTracerCntl%lnTIGRst) THEN
    ibcel = CellInfo(icel)%basecellstr
    nFSR = BaseCellInfo(ibcel)%nFSR
    nFXR = BaseCellInfo(ibcel)%nFXR
  ELSE
    nFSR = CellInfo(icel)%nFSR
    nFXR = CellInfo(icel)%nFXR
  ENDIF
  PinInfo(ipin)%nFsrMax = max(nFsr, PinInfo(ipin)%nFsrMax)
  PinInfo(ipin)%nFxrMax = max(nFxr, PinInfo(ipin)%nFxrMax)
  PinInfo(ipin)%lfuel = PinInfo(ipin)%lFuel .or. CellInfo(icel)%lFuel
ENDDO
END SUBROUTINE

SUBROUTINE ReadNinit_UsrGapPin(datalineIn)
!Read Pin CARD Input line
use param
use typedef
use allocs
use cntl
use ioutil,         only : terminate
use geom,           only : CellInfo,      PinInfo,        nz,       nPinType,     nPinType0,   &
                           nCellType0,     nCellType,                                          &
                           !--- CNJ Edit : Flexible Cell & Pin Numbering
                           CellTypeMap,  PinTypeMap
IMPLICIT NONE

character(256),intent(in) :: datalineIn
character(256) :: dataline
INTEGER :: ipin, icel
INTEGER :: i,j,k
INTEGER :: nFSR, nFXR    !Number of Flat Source and XS regsion
INTEGER :: GapType
dataline = datalineIn
read(dataline,*) ipin                          !Read the Pin Number
! ipin = PinTypeMap(ipin)   !--- CNJ Edit : Flexible Cell & Pin Numbering
CALL dmalloc(PinInfo(ipin)%cell,nz)

read(dataline,*) i , (PinInfo(ipin)%cell(k), k = 1, nz)  !Read the Pin Configuration

! !--- CNJ Edit : Flexible Cell & Pin Numbering
! DO k = 1, nz
!   PinInfo(ipin)%cell(k) = CellTypeMap(PinInfo(ipin)%cell(k))
! ENDDO

PinInfo(ipin)%ncell = nz; PinInfo(ipin)%lGap = .TRUE.
PinInfo(ipin)%lEmpty= FALSE;  PinInfo(ipin)%lUse = FALSE; PinInfo(ipin)%lFuel = FALSE
PinInfo(ipin)%lCentX = FALSE; PinInfo(ipin)%lCentY = FALSE; PinInfo(ipin)%lCentXY = FALSE 
PinInfo(ipin)%nFsrMax = 0; PinInfo(ipin)%nFxrMax = 0

GapType = CellInfo(PinInfo(ipin)%cell(1))%GapType

DO k= 1, nz
  icel = PinInfo(ipin)%cell(k)
  IF(k .EQ. 1) GapType = CellInfo(icel)%GapType
  IF(icel .gt. nCellType0 .or. CellInfo(icel)%lEmpty) then
    write(mesg,'(2(a,i3),a)') 'CELL (',icel,') in Pin ',ipin,' is not defined'
    CALL terminate(mesg) 
  ENDIF

  IF(.NOT. CellInfo(icel)%lGap) then
    write(mesg,'(2(a,i3),a)') 'GAPCELL (',icel,') in Pin ',ipin,' is not defined as gap cell'
    CALL terminate(mesg) 
  ENDIF  
  
  IF(GapType .NE. CellInfo(icel)%GapType) then
    write(mesg,'(2(a,i3),a)') 'GAPCELL (',icel,') in Pin ',ipin,' is not proper gap type'
    CALL terminate(mesg) 
  ENDIF
  
  nFSR = CellInfo(icel)%nFSR
  nFXR = CellInfo(icel)%nFXR
  PinInfo(ipin)%nFsrMax = max(nFsr, PinInfo(ipin)%nFsrMax)
  PinInfo(ipin)%nFxrMax = max(nFxr, PinInfo(ipin)%nFxrMax)
  PinInfo(ipin)%lfuel = PinInfo(ipin)%lFuel .or. CellInfo(icel)%lFuel
  PinInfo(ipin)%GapType = GapType
ENDDO
END SUBROUTINE

SUBROUTINE Init_GapPin
!Initiation 
USE param
USE typedef
USE allocs
USE cntl
USE ioutil,      ONLY : terminate
USE Geom,        ONLY : nCellType0,   nCellType,     nPinType0,      nPinType,     nz,     &
                        CellInfo,     PinInfo
IMPLICIT NONE
INTEGER :: icel, ipin
INTEGER :: icel1, icel2, icel3, ipin1, ipin2, ipin3

icel1 = nCellType0 - 2;  icel2 = nCellType0 - 1;  icel3 = nCellType0
ipin1 = nPinType0 -2;    ipin2 = nPinType0 - 1;   ipin3 = nPinType0  

DO ipin = ipin1, ipin3
  CALL dmalloc(PinInfo(ipin)%cell,nz)
  PinInfo(ipin)%ncell = nz
  PinInfo(ipin)%lEmpty= FALSE
  PinInfo(ipin)%lUse = FALSE
  PinInfo(ipin)%lFuel = FALSE
ENDDO



PinInfo(ipin1)%Cell(1:nz) = icel1
PinInfo(ipin1)%lGap = TRUE; PinInfo(ipin1)%GapType = 1
PinInfo(ipin2)%Cell(1:nz) = icel2
PinInfo(ipin2)%lGap = TRUE; PinInfo(ipin2)%GapType = 2
PinInfo(ipin3)%Cell(1:nz) = icel3
PinInfo(ipin3)%lGap = TRUE; PinInfo(ipin2)%GapType = 3
DO ipin = ipin1, ipin3
  icel = PinInfo(ipin)%Cell(1)
  PinInfo(ipin)%nFsrMax = CellInfo(icel)%nFsr
  PinInfo(ipin)%nFxrMax = CellInfo(icel)%nFxr
  PinInfo(ipin)%lFuel = CellInfo(icel)%lFuel
ENDDO
END SUBROUTINE

SUBROUTINE ReadNInit_ASSEMBLY(indev,datalineIn)
!Read ASSEMBLY CARD Input line
use param
use typedef
use allocs
use cntl
USE PE_MOD,         ONLY : PE
use files,          only : io5,         io8    
use inputcards ,    only : oneline,     probe
use ioutil,         only : terminate,   toupper,       IFnumeric,   nfields,     message
use inputcards ,    only : oneline,     probe
use geom,           only : AsyInfo,     nAsyType0,     nAsyType,     nCellX0,    PinInfo,      &
                           nPinType,    nPinType0,     nCellType0,   nCellx,     lEdge,        &
                           lGap
                           
implicit none
INTEGER :: AsyConfig0(nCellMax),AsyConfig(ncellxmax,ncellxmax)
INTEGER,intent(iN) :: indev 
character(256),intent(in) :: datalineIn

character(512) :: dataline
INTEGER :: ipin, icel, igap, iii
INTEGER :: i,j,k,n,ix,iy,jx,jy
INTEGER :: iasy, iAngAsy         !Assembly Number and Symmetry Anlge
INTEGER :: nread, nxread, nline, nfld
INTEGER :: jfr, jto,jfr0,jto0
LOGICAL :: master
logical,save :: first
data first /.TRUE./

IF(first) then
  Allocate(AsyInfo(nAsyType))
  DO i = 1, nAsyType
    AsyInfo(i)%lempty = .true. 
  ENDDO
  first = .false.
ENDIF



dataline = datalinein
nfld = nfields(dataline)
read(dataline, *) iasy

AsyInfo(iasy)%lEmpty = FALSE;  AsyInfo(iasy)%lFuel = FALSE
AsyInfo(iasy)%luse = FALSE
AsyInfo(iasy)%lCentX = FALSE; AsyInfo(iasy)%lCentY = FALSE; 
AsyInfo(iasy)%lCentXY = FALSE; 
AsyInfo(iasy)%GapType = 0

read(dataline, *) iii , iAngAsy
IF(nfld .GE. 3) THEN
  read(dataline, *) iii , iAngAsy, igap
  AsyInfo(iasy)%GapType = iGap
ENDIF




nline = nCellX0
nxread = nCellX0/2 + mod(nCellX0,2)
nread = nline*nline 
IF(iAngAsy .eq. 90) then
  nread = nxread*nxread
elseIF(iAngAsy .eq. 45) then
  nread = nxread*(1+nxread)/2
ENDIF

jfr = 1
jto = 0
!Read Assembly Configuration 
master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(Master) CALL message(io8,FALSE, FALSE,oneline)
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  n = nFields(oneline)
  jto = jfr + n -1
  read(oneline, *) (AsyConfig0(j), j = jfr, jto)
  !Pin Check
  DO j = jfr, jto
    ipin = AsyCOnfig0(j)
    IF(ipin .gt. nPinType0 .or. PinInfo(ipin)%lempty) then
        write(mesg,'(2(a,i3),a)') 'Pin (',ipin,') in Assembly ',iasy,' not defined'
        CALL terminate(mesg)       
    ENDIF
  ENDDO
  jfr = jfr + n 
  IF(jto .eq. nread) exit
ENDDO
jfr = 1;
jto = nCellX;
AsyConfig = 0 
IF(lGap) then    !Gap Cell Assignment 
  AsyConfig(2:nCellX-1, 1) = nPinType0-2
  AsyConfig(2:nCellX-1, nCellX) = nPinType0-2
  
  AsyConfig(1, 2:nCellX-1) = nPinType0-1; AsyConfig(nCellX, 2:nCellX-1) = nPinType0-1
  !Conner Gap
  AsyConfig(1, 1) = nPinType0;AsyConfig(nCellX, 1) = nPinType0;
  AsyConfig(1, nCellX) = nPinType0;AsyConfig(nCellX, nCellX) = nPinType0;
  jfr = 2; jto = nCellX-1
ENDIF

!One Dimensional Array -> Two Dimensional Array
SELECT CASE(iAngAsy)
  case(360)
    k=0
    DO i = jfr, jto
      DO j = jfr, jto
        k = k + 1
        AsyConfig(j, i) = AsyConfig0(k)  
      ENDDO
    ENDDO
  case(90)
    k=0
    jto0 = jto; jfr0 = jto - nxread +1
    DO i = jfr0, jto0
      DO j = jfr0, jto0
        k = k + 1
        AsyConfig(j, i) = AsyConfig0(k)
      ENDDO
    ENDDO
  CASE(45)
    k=0
    jto0 = jto; jfr0 = jto - nxread +1
    DO i = jfr0, jto0
      DO j = jfr0, i !jfr0+(i-jfr0+1)
        k = k + 1
        AsyConfig(j, i) = AsyConfig0(k)
      ENDDO
    ENDDO
    !45 -> 90 Symmetry
    k= nread+1 
    DO i = jto0, jfr0, -1
      DO j = i, jfr0,-1
      !DO j = jto0, i-1,-1
         k=k-1
         AsyConfig(i, j) = AsyConfig0(k)
      ENDDO
    ENDDO         
END SELECT
IF( iAngAsy .NE. 360) then
  !90 Symmetry -> 360 
  DO i = jfr0, jto0
    iy = nCellX - i + 1 
    DO j = jfr0, jto0
      jx = nCellX - j + 1 
      AsyConfig(jx, i) = AsyConfig(j, i)
      AsyConfig(j, iy) = AsyConfig(j, i)
      AsyConfig(jx, iy) = AsyConfig(j, i)
    ENDDO 
  ENDDO
ENDIF
AsyInfo(iasy)%nx = nCellX
AsyInfo(iasy)%ny = nCellX
AsyInfo(iasy)%nxy = nCellX*nCellX
CALL Dmalloc(AsyInfo(iasy)%Pin,AsyInfo(iasy)%nxy)
CALL Dmalloc(AsyInfo(iasy)%Pin2DIdx,nCellx,nCellx)
!Temporal Variable -> Assembly Type  
!AsyInfo
k=0 
DO i = 1, nCellX
  DO j = 1, nCellx
    k = k + 1 
    AsyInfo(iasy)%pin(k) = AsyConfig(j, i)
    AsyInfo(iasy)%pin2DIdx(j, i) = k
  ENDDO
ENDDO

AsyInfo(iasy)%nFuelPin = 0
DO k = 1, AsyInfo(iasy)%nxy
  ipin = AsyInfo(iasy)%pin(k);
  AsyInfo(iasy)%lfuel = AsyInfo(iasy)%lfuel .OR. PinInfo(ipin)%lfuel
  IF(PinInfo(ipin)%lfuel) AsyInfo(iasy)%nFuelPin = AsyInfo(iasy)%nFuelPin + 1
ENDDO
END SUBROUTINE

SUBROUTINE ReadNInit_AsyGapConf(indev, datalineIn)
USE PARAM
USE TYPEDEF
USE ALLOCS 
USE CNTL
USE PE_MOD,         ONLY : PE
USE FILES,              ONLY : io5,          io8
USE inputcards ,        ONLY : oneline,      probe
USE GEOM,               ONLY : AsyGap,       nAsyGapType,              &
                               AsyInfo,      PinInfo,                  &
                               nPinType,     nPinType0,    nCellX0,       nCellX      
USE ioutil,             ONLY : terminate,    toupper,       nfields,   &
                               message
USE BasicOperation,     ONLY : CP_CA
IMPLICIT NONE
INTEGER, INTENT(IN) :: indev
CHARACTER(256), INTENT(IN) :: DataLineIn

LOGICAL,SAVE :: First
LOGICAL :: MASTER
LOGICAL :: lSimInp
INTEGER :: UsrGapConfSim(3,3)
INTEGER :: UsrGapConf(nCellXMax,nCellXMax)
INTEGER :: ChkRange(4, 4, 3)
INTEGER :: ix, iy, i, k, l, igap
INTEGER :: n

character(512) :: dataline
CHARACTER(10) :: InpType


DATA FIRST /.TRUE./

IF(.NOT. nTracerCntl%lUsrGap .AND. .NOT. nTracerCntl%lGap) THEN
  IF(Master) CALL message(io8,FALSE, FALSE, 'Warning : ASYGAP Card is not neccesary in this input deck')
  RETURN
ENDIF

MASTER = PE%MASTER

IF(First) THEN
  ALLOCATE(AsyGap(nAsyGapType))
  first = .false.
  DO i = 1, nAsyGapType
    AsyGap(i)%lEmpty = .FALSE.
  ENDDO
  First = .FALSE.
ENDIF

Dataline = DatalineIn
n = nfields(Dataline)
lSimInp = .FALSE.
IF(n .EQ. 2) THEN
  READ(DataLine, *) igap, InpType
  CALL toupper(InpType)
  IF(InpType .EQ. 'S') lSimInp = .TRUE.
  IF(InpType .EQ. 'D') lSimInp = .FALSE.
ELSE
   READ(DataLine, *) igap
ENDIF

CALL CP_CA(UsrGapConf(1:nCellX, 1:nCellX), 0, nCellX, nCellX)

IF(lSimInp) THEN
  i=0
  DO WHILE(TRUE)
    READ(indev, '(a256)') oneline
    IF(Master) CALL message(io8,FALSE, FALSE, oneline)
    IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
    i=i+1;
    READ(oneline, *) UsrGapConfSim(1:3, i)
    IF(i .EQ. 3) EXIT
  ENDDO
ELSE
  i=0
  DO WHILE(TRUE)
    READ(indev, '(a256)') oneline
    IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
    i=i+1;
    READ(oneline, *) UsrGapConf(1:nCellX, i)
    IF(i .EQ. nCellX) EXIT
  ENDDO
ENDIF

IF(lSimInp) THEN
  UsrGapConf(1, 1) = UsrGapConfSim(1, 1);      UsrGapConf(nCellX, 1) = UsrGapConfSim(3, 1)
  UsrGapConf(1, nCellX) = UsrGapConfSim(1, 3); UsrGapConf(nCellX, nCellX) = UsrGapConfSim(3, 3)
  
  UsrGapConf(2:nCellX-1, 1) = UsrGapConfSim(2, 1)  ! TOP
  UsrGapConf(2:nCellX-1, nCellX) = UsrGapConfSim(2,3) ! Bottom
  
  UsrGapConf(1, 2:nCellX-1) = UsrGapConfSim(1, 2)  ! LEFT
  UsrGapConf(nCellX, 2:nCellX-1) = UsrGapConfSim(3, 2)  ! Right
ENDIF

!Gap Type Check
ChkRange = RESHAPE((/      2, nCellX-1,        1,        1,  &
                           2, nCellX-1,   nCellX,   nCellX,  &
                           0,       -1,        0,       -1,  &
                           0,       -1,        0,       -1,  &
                           1,        1,        2, nCellX-1,  &
                      nCellX,   nCellX,        2, nCellX-1,  &
                           0,       -1,        0,       -1,  &
                           0,       -1,        0,       -1,  &
                           1,        1,        1,        1,  &
                      nCellX,   nCellX,        1,        1,  &
                           1,        1,   nCellX,   nCellX,  &
                      nCellX,   nCellX,   nCellX,   nCellX/),&
                      SHAPE(ChkRange))
DO l = 1, 3
  DO k = 1, 4
     DO iy = ChkRange(3, k, l), ChkRange(4, k, l)
       DO ix = ChkRange(1, k, l), ChkRange(2, k, l)
          i = UsrGapConf(ix, iy)
          IF(i .GT. nPinType0 - 3) THEN
            write(mesg,'(2(a,i3),a)') 'Pin (',i,') in ASYGAP (',igap,') is not defined'
            CALL terminate(mesg)             
          ENDIF
          IF(.NOT. PinInfo(i)%lGap) THEN
            write(mesg,'(2(a,i3),a)') 'Pin (',i,') in ASYGAP (',igap,') is not defined as GAPPIN'
            CALL terminate(mesg)                  
          ENDIF
          IF(PinInfo(i)%GapType .NE. l) THEN
            write(mesg,'((a,i3),a, i3)') 'ASYGAP (',igap,') is not constructed with proper gap type : Gap Pin', UsrGapConf(ix, iy)
            CALL terminate(mesg)             
          ENDIF
       ENDDO
     ENDDO
  ENDDO
ENDDO

AsyGap(igap)%nx = nCellX; AsyGap(igap)%ny = nCellX
CALL Dmalloc(AsyGap(igap)%GapPin, nCellX * nCellX)
CALL Dmalloc(AsyGap(igap)%GapPin2D, nCellX, nCellX)

i=0
DO iy = 1, nCellX
  DO ix = 1, nCellX
    i = i + 1
    AsyGap(igap)%GapPin(i) = UsrGapConf(ix, iy)
    AsyGap(igap)%GapPin2D(ix, iy) = UsrGapConf(ix, iy)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE ReadNInit_Core(indev)
use param
use typedef
use allocs
use cntl
USE PE_MOD,         ONLY : PE
use files,          only : io5,         io8    
use inputcards ,    only : oneline,     probe
use ioutil,         only : toupper,     IFnumeric,    nfields,    message,  terminate
use geom,           only : CORE,        AsyInfo,      albedo,     lRot,     lCbd,    lCoreAng
implicit none
INTEGER,intent(in) :: indev
INTEGER :: ia, ja
INTEGER :: j, jfr, jto
INTEGER :: nFieldsLine
LOGICAL :: master
ja = 0
master = PE%master
CALL dmalloc(Core%CoreMap, Core%nxya)
CALL dmalloc(Core%CoreIdx, Core%nxa,Core%nya)
jfr = 1

DO WHILE(TRUE)
  read(indev,'(a256)') oneline
  IF(Master )CALL message(io8,FALSE,FALSE,oneline)
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  ja = ja + 1
  nFieldsLine = nfields(oneline)
  jto = jfr + nFieldsLine -1
  read(oneline, *) (Core%CoreMap(j), j = jfr, jto)
  ia=0
  IF(lCoreAng .EQ. 360 .AND. nFieldsLine .NE. Core%nxa) THEN
    ia = (Core%nxa - nFieldsLine) / 2
  ENDIF
  DO j = jfr ,jto
    ia = ia + 1
    Core%CoreIdx(ia, ja) = j
  ENDDO
  jfr = jto + 1 
  IF(ja .eq. core%nya) exit
ENDDO



Core%RadSym = FALSE
Core%RadBC = VoidCell
Core%AxBC = RefCell
DO j = 1, 4
  IF(Albedo(j) .EQ. ZERO) then
    Core%RadSym(j) = TRUE
    Core%RadBC(j) = RefCell
    IF(lROT) THEN  !Rotational Boundary
      IF(j .EQ. NORTH) Core%RadBC(j) = RotCell
      IF(j .EQ. WEST) Core%RadBC(j) = RotCell
    ELSEIF(lCBD) THEN
      Core%RadBC(j) = CbdCell
    ENDIF
  ENDIF
  IF(lCbd .AND. Albedo(j) .NE. ZERO) THEN
    CALL TERMINATE('All Radial B.C should be reflective in the Checkerboard Configuration')
  ENDIF
ENDDO

IF(Albedo(5) .NE. ZERO) Core%AxBc(1) = VoidCell
IF(Albedo(6) .NE. ZERO) Core%AxBc(2) = VoidCell

END SUBROUTINE

SUBROUTINE ReadNInit_Barrel(dataline)
USE PARAM
USE TYPEDEF,      ONLY : MISCSTRUCT_TYPE
USE GEOM,         ONLY : MISCSTRUCT
USE ALLOCS
USE CNTL
USE PE_MOD,         ONLY : PE
USE FILES,          ONLY : IO5,         IO8    
USE INPUTCARDS ,    ONLY : ONELINE,     PROBE
USE IOUTIL,         ONLY : TOUPPER,     IFNUMERIC,    NFIELDS,    MESSAGE,   &
                           FNDCHARA,    TERMINATE

character(256),intent(in) :: dataline

INTEGER :: nspt, ipos(100)

CALL FndChara(dataline,ipos,nspt,SLASH)

IF(NSPT .LT. 2) CALL TERMINATE('INPUT ERROR IN BARREL CARD')
!Read Barrel Radius
READ(dataline(1:ipos(1)-1), *) MiscStruct%rad_ring(1:2, 0)

READ(dataline(ipos(1)+1:ipos(2)-1), *) MiscStruct%ring_plnbeg(0), MiscStruct%ring_plnend(0)

!Read Mixture ID
Read(dataline(ipos(2)+1:256), *) MiscStruct%mix_ring(0)             

IF(NSPT .GE. 3) Read(dataline(ipos(3)+1:256), *) MiscStruct%lVolCor_ring(0)             
MiscStruct%lBarrel = .TRUE.
nTracerCntl%lAutoBarrel = .TRUE.
nTracerCntl%lMiscStruct = .TRUE.
END SUBROUTINE

SUBROUTINE ReadNInit_RingStruct(dataline)
USE PARAM
USE TYPEDEF,      ONLY : MISCSTRUCT_TYPE
USE GEOM,         ONLY : MISCSTRUCT
USE ALLOCS
USE CNTL
USE PE_MOD,         ONLY : PE
USE FILES,          ONLY : IO5,         IO8    
USE INPUTCARDS ,    ONLY : ONELINE,     PROBE
USE IOUTIL,         ONLY : TOUPPER,     IFNUMERIC,    NFIELDS,    MESSAGE,   &
                           FNDCHARA,    TERMINATE

character(256),intent(in) :: dataline

INTEGER :: nspt, ipos(100), id, id0

CALL FndChara(dataline,ipos,nspt,SLASH)

IF(NSPT .LT. 1) CALL TERMINATE('INPUT ERROR IN BARREL CARD')
READ(dataline, *) id
!Read Barrel Radius
READ(dataline(1:ipos(1)-1), *) id0, MiscStruct%rad_ring(1:2, id)

READ(dataline(ipos(1)+1:ipos(2)-1), *) MiscStruct%ring_plnbeg(id), MiscStruct%ring_plnend(id)

!Read Mixture ID
Read(dataline(ipos(2)+1:256), *) MiscStruct%mix_ring(id)             
IF(NSPT .GE. 2) Read(dataline(ipos(3)+1:256), *) MiscStruct%lVolCor_ring(id)             
MiscStruct%nring = max(MiscStruct%nring, id)

MiscStruct%lRing = .TRUE.
nTracerCntl%lAutoRingStruct = .TRUE.
nTracerCntl%lMiscStruct = .TRUE.
END SUBROUTINE

!--- JSR Edit : nTIG Restart

SUBROUTINE ReadNInit_Cell_nTIG(Dataline0, lCAD, lGAP, lCellInit)
use param
use typedef
use allocs
use cntl,         only : nTracerCntl
use geom,         only :  CellInfo,  CellPitch,    nCellType, nBasecell, BasecellInfo
use ioutil,       only :  toupper,   IFnumeric,    nfields,   fndchara,  fndchara512, &
                          nfieldto, nfieldto512
use GeomTreatment, ONLY : AnnularFsrArea, AnnularRegionDivision
use BenchXs,      only :  MacXsBen
use Material_Mod, ONLY : Mixture
use SPH_mod,      ONLY : calcCellSSPH
use XSLIB_MOD,    ONLY : igresb, igrese,nofghel
implicit none
character(512),intent(in) :: dataline0
LOGICAL :: LCAD, lGAP, lCellInit

INTEGER :: ndiv(0:300) = 0, ndivx(0:100) = 0, ndivy(0:100) = 0, ipos(100) = 0, iReg(0:300) = 0
INTEGER :: icel,item,imat(100),nmat,nreg,idiv(100)
INTEGER :: ndatafield,ndata,ndatax,ndatay
INTEGER :: nspt
INTEGER :: i,j, k, m, n, iRegFr,iRegTo
INTEGER :: irfr, ir, iv,ist, ix, iy
INTEGER :: nFsrInFxr
INTEGER :: CCentIdx(2, 4)
INTEGER :: ibasecell ! --- 180723
character(512) :: dataline,chline, CentPos

REAL :: rr(300), subrr(100), subvol(100), irad(0:100)
REAL :: delx(100), dely(100)
REAL :: tvol,dvol,vol,theta,vol_

DATA CCentIdx / 1, 1,  1, 2,  2, 2,  2, 1 /

IF(.NOT. lCellInit) then
  allocate(CellInfo(nCellType))
  DO i = 1, nCellType
    CellInfo(i)%lempty = TRUE
    CellInfo(i)%luse = FALSE
  ENDDO
  lCellInit = TRUE
ENDIF
IF(LCAD) THEN
  CALL ReadNInit_CadCell(Dataline0)
  RETURN
ENDIF
IF(LGAP) THEN
  CALL ReadNInit_UsrGapCell(Dataline0)
  RETURN
ENDIF

dataline = dataline0
READ(dataline,*) icel, ibasecell
CellInfo(icel)%icel0 = icel
CellInfo(icel)%basecellstr = ibasecell
CellInfo(icel)%lempty = FALSE
CellInfo(icel)%lFuel = FALSE
CellInfo(icel)%lRes = FALSE
CellInfo(icel)%lgap = FALSE
CellInfo(icel)%lGd = FALSE
CellInfo(icel)%lMox = .FALSE.
CellInfo(icel)%lAIC = .FALSE.
nDataField = len_trim(dataline)

ndata = BaseCellInfo(ibasecell)%geom%inpn
rr = BaseCellInfo(ibasecell)%rr
ndiv(1:ndata) = BaseCellInfo(ibasecell)%geom%inpdiv(1:ndata)

CALL fndchara512(dataline,ipos,nspt,SLASH)

IF(.not. BaseCellInfo(ibasecell)%lRect) then
  iRegFr=0
  iRegTo=nData
  read(dataline(ipos(1)+1:nDataField),*) (iReg(i),i=iRegTo,iRegFr,-1)
  CellInfo(icel)%iReg0(iRegFr:iRegTo) = iReg(iRegFr:iRegTo)
  
  imat=0; irad=0._8; idiv=0
  imat(1)=iReg(ndata); nmat=1
  irad(1)=rr(ndata)
  idiv(1)=ndiv(ndata)
  do i=iRegTo-1,iRegFr,-1
    if (iReg(i).eq.iReg(i+1)) then
        if (i.eq.0) then
            idiv(nmat)=idiv(nmat)+1
            irad(nmat)=CellPitch/dsqrt(PI)
        else
            idiv(nmat)=idiv(nmat)+ndiv(i)
        endif
        cycle
    endif
    nmat=nmat+1
    imat(nmat)=iReg(i)
    if (i.gt.0) then
        irad(nmat)=rr(i)
        idiv(nmat)=ndiv(i)
    else
        irad(nmat)=CellPitch/dsqrt(PI)
        idiv(nmat)=1
    endif
  enddo
  allocate(CellInfo(icel)%matidx(nmat))
  CellInfo(icel)%nmat=nmat
  CellInfo(icel)%matidx(1:nmat)=imat(1:nmat)
  allocate(CellInfo(icel)%matrad(nmat))
  CellInfo(icel)%matrad(1:nmat)=irad(1:nmat)
  
  nreg = sum(idiv(1:nmat))
  CellInfo(icel)%nreg_cp=nreg
  allocate(CellInfo(icel)%rad_cp(nreg),CellInfo(icel)%fxrvol(nreg))
  CellInfo(icel)%rad_cp=0._8
  vol=0._8
  k=1
  do i = 1, nmat
    if (i.eq.nmat) then
        CellInfo(icel)%fuelgapcldvol=vol
        CellInfo(icel)%cldfxridx=k-1
        CellInfo(icel)%nmodfxr=nreg-k+1
        CellInfo(icel)%invnmodfxr=1._8/real(CellInfo(icel)%nmodfxr,8)
    endif
    vol_ = PI*(irad(i)*irad(i)-irad(i-1)*irad(i-1))/idiv(i)  
    do j = 1, idiv(i)    
      vol = vol + vol_
      CellInfo(icel)%rad_cp(k) = dsqrt(vol*INVPI)  
      CellInfo(icel)%fxrvol(k) = vol_
      k=k+1
    enddo
  enddo
  
 DO j=0, ndata
   IF(nTracerCntl%lXsLib) THEN
     CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. Mixture(ireg(j))%lfuel
     CellInfo(icel)%lres = CellInfo(icel)%lres .or. Mixture(ireg(j))%lres
     CellInfo(icel)%lGd = CellInfo(icel)%lGd .or. Mixture(ireg(j))%lGd
     CellInfo(icel)%lMOX = CellInfo(icel)%lMOX .or. Mixture(ireg(j))%lMOX
     CellInfo(icel)%lAIC = CellInfo(icel)%lAIC .or. Mixture(ireg(j))%lAIC
     CellInfo(icel)%lsSPH = CellInfo(icel)%lfuel .or. CellInfo(icel)%lAIC
   ELSE
     CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. MacXsBen(ireg(j))%lfuel
     CellInfo(icel)%lCR = CellInfo(icel)%lCR .or. MacXsBen(ireg(j))%lCR
   ENDIF
 ENDDO
 

 !! SPH factor
 ndiv(1:ndata) = BaseCellInfo(ibasecell)%geom%inpdiv(1:ndata)
 DO j = 1, ndata
   ndiv(j) = BaseCellInfo(ibasecell)%geom%inpdiv(ndata-j+1)
 END DO
 
 !Allocation
 CALL dmalloc(CellInfo(icel)%iReg, BaseCellInfo(ibasecell)%nFSR)
 nFsrInFxr = BaseCellInfo(ibasecell)%nDivAzi
 ndiv(0) =1 
 
 DO j = 1, nFsrInFxr
   CellInfo(icel)%ireg(j) = ireg(0)
 ENDDO
 irfr = 1
 DO j = 1, ndata
  DO ir = 1, ndiv(j)
     ist = (irfr+ir-1)*nFsrInFxr
     DO iv = 1, nFsrInFxr
       CellInfo(icel)%iReg(ist+iv) = ireg(j)
     ENDDO
   ENDDO
   irfr = irfr + ndiv(j)
 ENDDO
ELSE  !Rectangular Geomemtry
 ndatax = BaseCellInfo(ibasecell)%geom%inpnx; ndatay = BaseCellInfo(ibasecell)%geom%inpny
 iRegFr=1
 iRegTo=nDatax*nDatay
 read(dataline(ipos(1)+1:nDataField),*) (iReg(i),i=iRegfr, iRegto)
 CALL dmalloc(CellInfo(icel)%iReg, BaseCellInfo(ibasecell)%nFSR)
 !Allocation
 DO i = iregfr,iregto
   IF(nTracerCntl%lXsLib) THEN
     CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. Mixture(ireg(i))%lfuel
   ELSE
     CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. MacXsBen(ireg(i))%lfuel
   ENDIF  
 ENDDO
 ndivx(1:ndatax) = BaseCellInfo(ibasecell)%geom%inpdivx(1:ndatax)
 ndivy(1:ndatay) = BaseCellInfo(ibasecell)%geom%inpdivy(1:ndatay)
 m = 0; n = 0
 DO iy = 1, ndatay
   DO ix = 1, ndatax
     m = m + 1; n = 0
     DO j = ndivy(iy-1) + 1, ndivy(iy)
       DO i = ndivx(ix-1) + 1, ndivx(ix)
         n = n + 1
         k = BaseCellInfo(ibasecell)%Geom%nx * (j - 1) + i
         CellInfo(icel)%iReg(k) = ireg(m)
       ENDDO
     ENDDO
   ENDDO
 ENDDO
ENDIF
CellInfo(icel)%CellInput = Dataline ! **************
END SUBROUTINE

SUBROUTINE ReadNInit_BaseCell_nTIG(Dataline0, lCAD, lGAP, lBaseCellInit) ! --- 180722
use param
use typedef
use allocs
use cntl,         only : nTracerCntl
use geom,         only :  CellInfo,  CellPitch,    nCellType, nBasecell, BaseCellInfo
use ioutil,       only :  toupper,   IFnumeric,    nfields,   fndchara,  fndchara512, &
                          nfieldto, nfieldto512
use GeomTreatment, ONLY : AnnularFsrArea, AnnularRegionDivision
use BenchXs,      only :  MacXsBen
use Material_Mod, ONLY : Mixture
use SPH_mod,      ONLY : calcCellSSPH
use XSLIB_MOD,    ONLY : igresb, igrese,nofghel
implicit none
character(512),intent(in) :: dataline0
LOGICAL :: LCAD, lGAP, lBaseCellInit

INTEGER :: ndiv(0:300) = 0, ndivx(0:100) = 0, ndivy(0:100) = 0, ipos(100) = 0, iReg(0:300) = 0
INTEGER :: icel,item,item2 ! --- 180626 JSR
INTEGER :: ndatafield,ndata,ndatax,ndatay
INTEGER :: nspt
INTEGER :: i,j, k, m, n
INTEGER :: irfr, ir, iv,ist, ix, iy
INTEGER :: nFsrInFxr
INTEGER :: CCentIdx(2, 4)
character(512) :: dataline,chline, CentPos

REAL :: rr(300), subrr(100), subvol(100)
REAL :: delx(100), dely(100)
REAL :: Pitch,HalfPitch
REAL :: tvol,dvol,vol,theta

INTEGER :: ibasecell ! --- 180622 JSR

DATA CCentIdx / 1, 1,  1, 2,  2, 2,  2, 1 /

IF(.NOT. lBaseCellInit) then
  allocate(BaseCellInfo(nBasecell))
  DO i = 1, nBasecell
    BaseCellInfo(i)%lempty = TRUE
    BaseCellInfo(i)%luse = FALSE
  ENDDO
  lBaseCellInit = TRUE
ENDIF

Pitch = CellPitch;  HalfPitch = 0.5_8*CellPitch
dataline = dataline0
READ(dataline,*) ibasecell !, ibasecell ! --- 180622 JSR (+ ibasecell)
! icel = CellTypeMap(icel)   !--- CNJ Edit : Flexible Cell & Pin Numbering
icel = ibasecell
BaseCellInfo(ibasecell)%lempty = FALSE
BaseCellInfo(ibasecell)%lrect = FALSE
BaseCellInfo(ibasecell)%lCCell = .FALSE.
BaseCellInfo(ibasecell)%geom%lCCent  = FALSE
BaseCellInfo(ibasecell)%lFuel = FALSE
BaseCellInfo(ibasecell)%lRes = FALSE
BaseCellInfo(ibasecell)%lgap = FALSE
BaseCellInfo(ibasecell)%lGd = FALSE
BaseCellInfo(ibasecell)%lMox = .FALSE.
BaseCellInfo(ibasecell)%lCentX = FALSE; BaseCellInfo(icel)%lCentY = FALSE; BaseCellInfo(icel)%lCentXY = FALSE; 
nDataField = len_trim(dataline)


CALL fndchara512(dataline,ipos,nspt,SLASH)

IF(nspt .eq. 1 .or. nspt .eq. 2) THEN
  BaseCellInfo(icel)%nDivAzi = 8
  ndata = nfieldto512(dataline,SLASH)-1 ! --- 180622 JSR (icell -1 -> icell, ibasecell -2)
  
  BaseCellInfo(icel)%Geom%lCircle = TRUE
  BaseCellInfo(icel)%Geom%lRect = FALSE
  BaseCellInfo(icel)%Geom%nCircle = 0
  BaseCellInfo(icel)%Geom%nx = 3
  BaseCellInfo(icel)%Geom%ny = 3
  
  BaseCellInfo(icel)%geom%inpn = ndata
  read(dataline(ipos(1)+1:nDataField),*) (nDiv(i),i=1,ndata)
  BaseCellInfo(icel)%geom%inpdiv(1:ndata) = nDiv(1:ndata) ! --- 180723
  
  !Circle Division Read
  DO i =1, ndata
    BaseCellInfo(icel)%geom%nCircle = BaseCellInfo(icel)%geom%nCircle + nDiv(i)
  ENDDO
  BaseCellInfo(icel)%nFXR = BaseCellInfo(icel)%geom%nCircle + 1
  !Off-Center Pincell
  BaseCellInfo(icel)%geom%cx = 0 
  BaseCellInfo(icel)%geom%cy = 0
  IF(nspt .eq. 2) then
    BaseCellInfo(icel)%geom%lCCent  = TRUE
    BaseCellInfo(icel)%lCCell = TRUE
    BaseCellInfo(icel)%nDivAzi = 2
    read(dataline(ipos(2)+1:nDataField),*) CentPos
    IF(CentPos .EQ. 'SW')  BaseCellInfo(icel)%geom%CCentType = 1
    IF(CentPos .EQ. 'NW')  BaseCellInfo(icel)%geom%CCentType = 2
    IF(CentPos .EQ. 'NE')  BaseCellInfo(icel)%geom%CCentType = 3
    IF(CentPos .EQ. 'SE')  BaseCellInfo(icel)%geom%CCentType = 4
    !IF( BaseCellInfo(icel)%geom%CCentType) CALL TERMINATE('Cell Input Error')
  ENDIF
  BaseCellInfo(icel)%nFSR = (BaseCellInfo(icel)%geom%nCircle + 1)*BaseCellInfo(icel)%nDivAzi
ELSE
  BaseCellInfo(icel)%lRect = TRUE
  BaseCellInfo(icel)%geom%lrect = TRUE
  BaseCellInfo(icel)%geom%lcircle = FALSE
  BaseCellInfo(icel)%Geom%nCircle =0
  BaseCellInfo(icel)%geom%cx = 0 
  BaseCellInfo(icel)%geom%cy = 0
  ndatax = nfieldto(dataline,SLASH) ! --- 180628 JSR (exclude for when it accounts basecell idx either)
  read(dataline, *) item, (delx(i), i = 1, ndatax-1)
  delx(ndatax) = CellPitch - sum(delx(1:ndatax-1))
  
  !---BYS edit for GC GeomReStart
  BaseCellInfo(icel)%geom%inpnx=ndatax
  BaseCellInfo(icel)%geom%inpdelx(1:ndatax-1)=delx(1:ndatax-1)
  !---BYS edit for GC GeomReStart END
  
  read(dataline(ipos(2)+1:ipos(3)-1),*) (ndivx(i), i=1,ndatax)
  !---BYS edit for GC GeomReStart
  BaseCellInfo(icel)%geom%inpdivx(1:ndatax)=ndivx(1:ndatax)
  !---BYS edit for GC GeomReStart END
    
  BaseCellInfo(icel)%geom%nx =1
  DO i = 1, ndatax
    BaseCellInfo(icel)%geom%nx = BaseCellInfo(icel)%geom%nx + ndivx(i)
    delx(i) = delx(i) / ndivx(i)
  ENDDO
  chline = dataline(ipos(1)+1:256)
  ndatay = nfieldto(chline,SLASH) + 1
  read(dataline(ipos(1)+1:ipos(2)-1), *) (dely(i), i = 1, ndatay-1)
  !IF(ndatay .EQ. 1) dely(1) = CellPitch
  dely(ndatay) = CellPitch - sum(dely(1:ndatay-1))  
  
  !---BYS edit for GC GeomReStart
  BaseCellInfo(icel)%geom%inpny=ndatay
  BaseCellInfo(icel)%geom%inpdely(1:ndatay-1)=dely(1:ndatay-1)
  !---BYS edit for GC GeomReStart END
  
  read(dataline(ipos(3)+1:256),*) (ndivy(i), i=1,ndatay)
    !---BYS edit for GC GeomReStart
  BaseCellInfo(icel)%geom%inpdivy(1:ndatay)=ndivy(1:ndatay)
  !---BYS edit for GC GeomReStart END
  
  BaseCellInfo(icel)%geom%ny =1
  DO i = 1, ndatay
    BaseCellInfo(icel)%geom%ny = BaseCellInfo(icel)%geom%ny + ndivy(i)
    dely(i) = dely(i) / ndivy(i)
  ENDDO
  
  BaseCellInfo(icel)%nFSR = (BaseCellInfo(icel)%geom%nx-1)*(BaseCellInfo(icel)%geom%ny-1)
  !BaseCellInfo(icel)%nFXR = 1
  BaseCellInfo(icel)%nFXR = nDataX * nDataY
  BaseCellInfo(icel)%geom%nx=BaseCellInfo(icel)%geom%nx-1
  BaseCellInfo(icel)%geom%ny=BaseCellInfo(icel)%geom%ny-1
ENDIF
!Composition Check
BaseCellInfo(icel)%Geom%lx = CellPitch; BaseCellInfo(icel)%geom%ly = CellPitch
BaseCellInfo(icel)%Geom%x(1) = -CellPitch*Half; BaseCellInfo(icel)%Geom%x(2) = CellPitch*Half
BaseCellInfo(icel)%Geom%y(1) = -CellPitch*Half; BaseCellInfo(icel)%Geom%y(2) = CellPitch*Half
!CALL dmalloc(BaseCellInfo(icell)%FxrIdxSt,BaseCellInfo(icell)%nFxr)
!CALL dmalloc(BaseCellInfo(icell)%FxrIdxSt,BaseCellInfo(icell)%nFxr)
!
IF(.not. BaseCellInfo(icel)%lRect) then
  IF(BaseCellInfo(icel)%Geom%lCCent) THEN
    HalfPitch = HalfPitch  * 2._8
    !BaseCellInfo(icel)%geom%cx = BaseCellInfo(icel)%Geom%x(CCentIdx(1, BaseCellInfo(icel)%Geom%CCentType))
    !BaseCellInfo(icel)%geom%cy = BaseCellInfo(icel)%Geom%y(CCentIdx(2, BaseCellInfo(icel)%Geom%CCentType))    
  ENDIF
  read(dataline,*) item, (rr(j),j=ndata,1,-1)    !Read Annunal ring Info
  BasecellInfo(icel)%rr(1:ndata) = rr(1:ndata) ! --- 180914
  read(dataline(ipos(1)+1:ipos(2)-1),*) (ndiv(j),j=ndata,1,-1)
  !lfuel and lRes Setting

  ! Calculate volume of each MATERIAL region (not each fxr)
  ALLOCATE(BaseCellInfo(icel)%vol(0:ndata))
  vol=0._8
  DO j=ndata,1,-1
      IF (rr(j).gt.HalfPitch) THEN
          theta=acos(HalfPitch/rr(j))
          !BaseCellInfo(icel)%vol(j-1)=Pitch*Pitch-rr(j)*(PI*rr(j)-4._8*(rr(j)*theta-HalfPitch*sin(theta)))-vol
          BaseCellInfo(icel)%vol(j)=4._8*rr(j)*(HalfPitch*sin(theta)+rr(j)*(PI*0.25_8-theta))-vol
      ELSE
          !BaseCellInfo(icel)%vol(j-1)=Pitch*Pitch-PI*rr(j)*rr(j)-vol
          BaseCellInfo(icel)%vol(j)=PI*rr(j)*rr(j)-vol
      ENDIF
      !vol=vol+BaseCellInfo(icel)%vol(j-1)
      vol=vol+BaseCellInfo(icel)%vol(j)
  ENDDO
  BaseCellInfo(icel)%vol(0)=Pitch*Pitch-vol  
  !! SPH factor
!  IF (nTracerCntl%lSSPH) call calcCellSSPH(CellInfo,icel,ndata,iReg,rr,ndiv)          !!!!!!!!!!!!!!!!!!!
  
  !Allocation
  CALL dmalloc(BaseCellInfo(icel)%iReg, BaseCellInfo(icel)%nFSR)                                                                                       !!!!!!!!!!!!!!!!!!
  CALL dmalloc(BaseCellInfo(icel)%FxrIdxSt,BaseCellInfo(icel)%nFxr)                                                                                    !!!!!!!!!!!!!!!!!!
  CALL dmalloc(BaseCellInfo(icel)%nFsrInFxr,BaseCellInfo(icel)%nFxr)                                                                                   !!!!!!!!!!!!!!!!!!
  CALL dmalloc(BaseCellInfo(icel)%MapFxr2FsrIdx, BaseCellInfo(icel)%nFSR, BaseCellInfo(icel)%nFxr)                                                         !!!!!!!!!!!!!!!!!!
  nFsrInFxr = BaseCellInfo(icel)%nDivAzi                                                                                                           !!!!!!!!!!!!!!!!!!
  ndiv(0) =1                                                                                                                                   !!!!!!!!!!!!!!!!!!
                                                                                                                                               !!!!!!!!!!!!!!!!!!
  rr(ndata+1) = 0                                                                                                                              !!!!!!!!!!!!!!!!!!
  tvol = PI * rr(1) * rr(1)                                                                                                                    !!!!!!!!!!!!!!!!!!
  BaseCellInfo(icel)%FxrIdxSt(1) = 1                                                                                                               !!!!!!!!!!!!!!!!!!
  BaseCellInfo(icel)%nFsrInFxr(1) = nFsrInFxr                                                                                                      !!!!!!!!!!!!!!!!!!
!  DO j = 1, nFsrInFxr                                                                                                                          !!!!!!!!!!!!!!!!!!
!    BaseCellInfo(icel)%ireg(j) = ireg(0)                                                                                                           !!!!!!!!!!!!!!!!!!
!  ENDDO                                                                                                                                        !!!!!!!!!!!!!!!!!!
  irfr = 1                                                                                                                                     !!!!!!!!!!!!!!!!!!
  CALL dmalloc(BaseCellInfo(icel)%geom%Circle,3,BaseCellInfo(icel)%geom%nCircle)                                                                       !!!!!!!!!!!!!!!!!!
  ! print *, ndata                                                                                                                             !!!!!!!!!!!!!!!!!!
  DO j = 1, ndata                                                                                                                              !!!!!!!!!!!!!!!!!!
  !  print *, j                                                                                                                                !!!!!!!!!!!!!!!!!!
    dvol = PI*(rr(j)*rr(j)-rr(j+1)*rr(j+1))/ndiv(j)                                                                                            !!!!!!!!!!!!!!!!!!
    tvol = PI*rr(j)*rr(j)                                                                                                                      !!!!!!!!!!!!!!!!!!
    IF(HalfPitch .LE. rr(j)) THEN                                                                                                              !!!!!!!!!!!!!!!!!!
      subrr(1:ndiv(j)) = AnnularRegionDivision(HalfPitch, rr(j:j+1), ndiv(j))                                                                  !!!!!!!!!!!!!!!!!!
    ENDIF                                                                                                                                      !!!!!!!!!!!!!!!!!!
    DO ir = 1, ndiv(j)                                                                                                                         !!!!!!!!!!!!!!!!!!
      BaseCellInfo(icel)%FxrIdxSt(irfr+ir) = (irfr+ir-1)*nFsrInFxr + 1   !Staring Index of local FSR index in a FXR                                !!!!!!!!!!!!!!!!!!
      BaseCellInfo(icel)%nFsrInFxr(irfr+ir) = nFsrinFxr                  !Number of FSRs which are beleong to a certain FXR                        !!!!!!!!!!!!!!!!!!
      IF(HalfPitch .GT. rr(j)) THEN                                                                                                            !!!!!!!!!!!!!!!!!!
        BaseCellInfo(icel)%geom%circle(3,irfr+ir-1) = sqrt(tvol/PI)        !Annularing                                                             !!!!!!!!!!!!!!!!!!
        tvol = tvol - dvol                                                                                                                     !!!!!!!!!!!!!!!!!!
      ELSE                                                                                                                                     !!!!!!!!!!!!!!!!!!
        BaseCellInfo(icel)%geom%circle(3,irfr+ir-1) = Subrr(ir)                                                                                    !!!!!!!!!!!!!!!!!!
      ENDIF                                                                                                                                    !!!!!!!!!!!!!!!!!!
      IF(BaseCellInfo(icel)%Geom%lCCent) THEN                                                                                                      !!!!!!!!!!!!!!!!!!
        BaseCellInfo(icel)%geom%circle(1,irfr+ir-1) = BaseCellInfo(icel)%Geom%x(CCentIdx(1, BaseCellInfo(icel)%Geom%CCentType))                            !!!!!!!!!!!!!!!!!!
        BaseCellInfo(icel)%geom%circle(2,irfr+ir-1) = BaseCellInfo(icel)%Geom%y(CCentIdx(2, BaseCellInfo(icel)%Geom%CCentType))                            !!!!!!!!!!!!!!!!!!
      ENDIF                                                                                                                                    !!!!!!!!!!!!!!!!!!
      !Assign the Mixture or composition number in a FSR                                                                                       !!!!!!!!!!!!!!!!!!
      ist = (irfr+ir-1)*nFsrInFxr                                                                                                              !!!!!!!!!!!!!!!!!!
!      DO iv = 1, nFsrInFxr                                                                                                                     !!!!!!!!!!!!!!!!!!
!        BaseCellInfo(icel)%iReg(ist+iv) = ireg(j)                                                                                                  !!!!!!!!!!!!!!!!!!
!      ENDDO                                                                                                                                    !!!!!!!!!!!!!!!!!!
    ENDDO                                                                                                                                      !!!!!!!!!!!!!!!!!!
    irfr = irfr + ndiv(j)                                                                                                                      !!!!!!!!!!!!!!!!!!
  ENDDO                                                                                                                                        !!!!!!!!!!!!!!!!!!
  !FXS region Index Starting Point-> FSR Region Idx                                                                                            !!!!!!!!!!!!!!!!!!
  DO ir = 1, BaseCellInfo(icel)%nFxr                                                                                                               !!!!!!!!!!!!!!!!!!
    ist = BaseCellInfo(icel)%FxrIdxSt(ir)                                                                                                          !!!!!!!!!!!!!!!!!!
    DO j = 1, BaseCellInfo(icel)%nFsrInFxr(ir)                                                                                                     !!!!!!!!!!!!!!!!!!
      BaseCellInfo(icel)%MapFxr2FsrIdx(j,ir) = ist+j-1                                                                                             !!!!!!!!!!!!!!!!!!
    ENDDO                                                                                                                                      !!!!!!!!!!!!!!!!!!
  ENDDO                                                                                                                                        !!!!!!!!!!!!!!!!!!
ELSE  !Rectangular Geomemtry
  !Allocation
  CALL dmalloc(BaseCellInfo(icel)%iReg, BaseCellInfo(icel)%nFSR)
  CALL dmalloc(BaseCellInfo(icel)%FxrIdxSt,BaseCellInfo(icel)%nFxr)
  CALL dmalloc(BaseCellInfo(icel)%nFsrInFxr,BaseCellInfo(icel)%nFxr)
  CALL dmalloc(BaseCellInfo(icel)%MapFxr2FsrIdx, BaseCellInfo(icel)%nFSR, BaseCellInfo(icel)%nFxr)                                                                                                                                !!!!!!!!!!!!!!!!!!
  i = 0    
  DO iy = 1, nDatay
    DO ix = 1, nDatax
      i = i + 1
      BaseCellInfo(icel)%nFsrInFxr(i) = ndivx(ix) * ndivy(iy)
    ENDDO
  ENDDO
  
  BaseCellInfo(icel)%FxrIdxSt(1) = 1
  
  DO i = 2, nDatay * nDatax
    BaseCellInfo(icel)%FxrIdxSt(i) = BaseCellInfo(icel)%FxrIdxSt(i-1) + BaseCellInfo(icel)%nFsrInFxr(i-1)
  ENDDO
  
  DO i = 2, 100
    ndivy(i) = ndivy(i-1) + ndivy(i)
    ndivx(i) = ndivx(i-1) + ndivx(i)
  ENDDO

  CALL Dmalloc(BaseCellInfo(icel)%Geom%delx, BaseCellInfo(icel)%Geom%nx)
  CALL Dmalloc(BaseCellInfo(icel)%Geom%dely, BaseCellInfo(icel)%Geom%ny)

  DO ix = 1, ndatax
    DO j= ndivx(ix-1) + 1, ndivx(ix)
      BaseCellInfo(icel)%Geom%DelX(j) = DelX(ix) 
    ENDDO
  ENDDO 

  DO iy = 1, ndatay
    DO j= ndivy(iy-1) + 1, ndivy(iy)
      BaseCellInfo(icel)%Geom%Dely(j) = Dely(iy) 
    ENDDO
  ENDDO

  m = 0; n = 0
  DO iy = 1, ndatay
    DO ix = 1, ndatax
      m = m + 1; n = 0
      DO j = ndivy(iy-1) + 1, ndivy(iy)
        DO i = ndivx(ix-1) + 1, ndivx(ix)
          n = n + 1
          k = BaseCellInfo(icel)%Geom%nx * (j - 1) + i
          BaseCellInfo(icel)%MapFxr2FsrIdx(n, m) = k                                            !!!!!!!!!!!!!!!!!!!!!!!!!!
!          BaseCellInfo(icel)%iReg(k) = ireg(m)                                                  !!!!!!!!!!!!!!!!!!!!!!!!!!
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  !BaseCellInfo(icel)%FxrIdxSt(1) = 1
  !BaseCellInfo(icel)%nFsrInFxr(1) = BaseCellInfo(icel)%nFSR
  
  !DO i = 1, BaseCellInfo(icel)%nFSR
  !  BaseCellInfo(icel)%M/.apFxr2FsrIdx(i,1) = i
  !  BaseCellInfo(icel)%iReg(i) = iReg(1)
  !ENDDO
  !TH setting
  
ENDIF
BaseCellInfo(icel)%CellInput = Dataline
END SUBROUTINE

SUBROUTINE Init_GapCell_nTIG(DataLineIn, lBaseCellInit)
!Initialize Gap Cell
use param
use typedef
use allocs
use cntl
use geom,   only :  lGap, CellInfo,  CellPitch,     AsyPitch,  nCellType0,  nCellType,   &
                    nCellX0, BasecellInfo, nBaseCell, nBaseCell0
!use ioutil, only :  toupper,   IFnumeric,    nfields,   fndchara,    nfieldto
use BenchXs,only :  MacXsBen
implicit none
character(256),intent(in) :: datalineIn
LOGICAL :: lBaseCellInit

INTEGER :: icel, icel1, icel2, icel3
INTEGER :: i,j
INTEGER :: GapComp
REAL :: Pitch,HalfPitch,GapPitch
REAL  :: tvol,dvol
logical,save :: first
data first /.TRUE./

if(.NOT. lgap) return

IF(.NOT. lBaseCellInit) then
  allocate(BaseCellInfo(nBasecell))
  DO i = 1, nBasecell
    BaseCellInfo(i)%lempty = TRUE
    BaseCellInfo(i)%luse = FALSE
  ENDDO
  lBaseCellInit = TRUE
ENDIF

read(DataLineIn, *) lGap, GapComp

GapPitch = (AsyPitch - CellPitch*nCellX0)*HALF

icel1 = nBaseCell0 - 2;  icel2 = nBaseCell0 - 1;  icel3 = nBaseCell0
DO icel = icel1, icel3
  BaseCellInfo(icel)%lfuel = FALSE; BaseCellInfo(icel)%lempty = FALSE;
  BaseCellInfo(icel)%lCentX = FALSE; BaseCellInfo(icel)%lCentY = FALSE;
  BaseCellInfo(icel)%lCentXY = FALSE; BaseCellInfo(icel)%lgap = TRUE; 
  BaseCellInfo(icel)%Geom%Lcircle = FALSE
  BaseCellInfo(icel)%Geom%nCircle =0
  BaseCellInfo(icel)%Geom%lCCent = FALSE
  BaseCellInfo(icel)%lRect = TRUE
  BaseCellInfo(icel)%geom%cx = 0 
  BaseCellInfo(icel)%geom%cy = 0
ENDDO

DO icel = icel1, icel3
  CALL Dmalloc(BaseCellInfo(icel)%FxrIdxSt, 1) 
  CALL Dmalloc(BaseCellInfo(icel)%nFsrInFxr, 1)
  CALL Dmalloc(BaseCellInfo(icel)%MapFxr2FsrIdx,8,1)
  CALL Dmalloc(BaseCellInfo(icel)%iReg, 8)
ENDDO

!Up and Bottom Gap Cell
icel = icel1
BaseCellInfo(icel)%geom%lx = CellPitch 
BaseCellInfo(icel)%Geom%ly = GapPitch
BaseCellInfo(icel)%Geom%nx = 4                     !Region 
BaseCellInfo(icel)%Geom%ny = 2
BaseCellInfo(icel)%GapType = 1

!Left and Right Gap Cell
icel = icel2
BaseCellInfo(icel)%geom%ly = CellPitch
BaseCellInfo(icel)%Geom%lx = GapPitch
BaseCellInfo(icel)%Geom%nx = 2
BaseCellInfo(icel)%Geom%ny = 4
BaseCellInfo(icel)%GapType = 2
!Coner Gap Cell
icel = icel3
BaseCellInfo(icel)%Geom%lx = GapPitch
BaseCellInfo(icel)%Geom%ly = GapPitch
BaseCellInfo(icel)%Geom%nx = 2
BaseCellInfo(icel)%Geom%ny = 2
BaseCellInfo(icel)%GapType = 3
DO icel = icel1, icel3
  
  CALL Dmalloc(BaseCellInfo(icel)%Geom%delx, BaseCellInfo(icel)%Geom%nx)
  CALL Dmalloc(BaseCellInfo(icel)%Geom%dely, BaseCellInfo(icel)%Geom%ny)
  BaseCellInfo(icel)%Geom%lcircle = FALSE
  BaseCellInfo(icel)%Geom%nCircle = 0
  BaseCellInfo(icel)%Geom%lrect = TRUE
  BaseCellInfo(icel)%lfuel = .FALSE.
  !BaseCellInfo(icel)%ireg(1) = GapComp
  BaseCellInfo(icel)%nFXR = 1
  BaseCellInfo(icel)%nFSR = BaseCellInfo(icel)%geom%nx*BaseCellInfo(icel)%geom%ny
  BaseCellInfo(icel)%nFsrInFxr(1) = BaseCellInfo(icel)%nFSR
  BaseCellInfo(icel)%FxrIdxSt(1) = 1
  DO i = 1, BaseCellInfo(icel)%nFSR
    BaseCellInfo(icel)%ireg(i) = GapComp
    BaseCellInfo(icel)%MapFxr2FsrIdx(i, 1) = i
  ENDDO
  BaseCellInfo(icel)%Geom%x(1) = -BaseCellInfo(icel)%Geom%lx * Half
  BaseCellInfo(icel)%Geom%x(2) = BaseCellInfo(icel)%Geom%lx * Half
  BaseCellInfo(icel)%Geom%y(1) = -BaseCellInfo(icel)%Geom%ly*Half
  BaseCellInfo(icel)%Geom%y(2) = BaseCellInfo(icel)%Geom%ly*Half
  
  CALL Dmalloc(BaseCellInfo(icel)%Geom%delx, BaseCellInfo(icel)%Geom%nx)
  CALL Dmalloc(BaseCellInfo(icel)%Geom%dely, BaseCellInfo(icel)%Geom%ny)
  
  BaseCellInfo(icel)%Geom%DelX(:) = BaseCellInfo(icel)%Geom%lx / BaseCellInfo(icel)%Geom%nx
  BaseCellInfo(icel)%Geom%DelY(:) = BaseCellInfo(icel)%Geom%ly / BaseCellInfo(icel)%Geom%ny
  !BaseCellInfo(icel)%Geom%DelX()
  continue
ENDDO

icel1 = nCellType0 - 2;  icel2 = nCellType0 - 1;  icel3 = nCellType0
DO icel = icel1, icel3
  CALL Dmalloc(CellInfo(icel)%iReg, 8)
ENDDO

!Up and Bottom Gap Cell
icel = icel1
CellInfo(icel)%nFXR = 1
CellInfo(icel)%nFSR = 8

!Left and Right Gap Cell
icel = icel2
CellInfo(icel)%nFXR = 1
CellInfo(icel)%nFSR = 8

!Coner Gap Cell
icel = icel3
CellInfo(icel)%nFXR = 1
CellInfo(icel)%nFSR = 4

DO icel = icel1, icel3
  CellInfo(icel)%ireg(:) = GapComp
ENDDO

END SUBROUTINE

SUBROUTINE CopyBaseCell()
USE geom,           ONLY : nCellType, CellInfo, BaseCellInfo, nCellType0, nBaseCell0, lgap
USE CNTL,           ONLY : nTracerCntl
use BenchXs,        only : MacXsBen
use Material_Mod,   ONLY : Mixture
USE PE_Mod,         ONLY : PE
use SPH_mod,        ONLY : calcCellSSPH,calcAICCellSSPH
USE OMP_LIB
IMPLICIT NONE

INTEGER :: i, j, icel, ibasecell
INTEGER :: ndata, nfueldiv
INTEGER :: iReg(0:300) = 0, ndiv(0:300) = 0
REAL :: rr(300) = 0.

IF (lgap) THEN
  DO icel = nCellType0, nCellType0-2, -1 ! ---- 180622 JSR (If lGap. Though, unless there are three gap cells, this should be revised)
    CellInfo(icel)%basecellstr = nBaseCell0+(icel-nCellType0)
  END DO
END IF

IF (nCellType0 .NE. nCellType) THEN
  DO icel = nCellType0 + 1, nCellType
    CellInfo(icel)%basecellstr = CellInfo(icel - nCellType0)%basecellstr + nBaseCell0
  ENDDO
ENDIF

CALL omp_set_num_threads(PE%nThread)

!$OMP PARALLEL PRIVATE(i, j, ibasecell, ireg, rr, ndata, ndiv, nfueldiv)
!$OMP DO SCHEDULE(GUIDED)
DO icel = 1, nCellType

  ibasecell = CellInfo(icel)%basecellstr
  
  CellInfo(icel)%lEmpty = .FALSE.
  CellInfo(icel)%lUse = .TRUE.
  
  CellInfo(icel)%vol                => BasecellInfo(ibasecell)%vol     
  
  CellInfo(icel)%nDivAzi            = BasecellInfo(ibasecell)%nDivAzi     
  CellInfo(icel)%nBd                = BasecellInfo(ibasecell)%nBd         
  CellInfo(icel)%nFSR               = BasecellInfo(ibasecell)%nFSR        
  CellInfo(icel)%nFXR               = BasecellInfo(ibasecell)%nFXR        
  CellInfo(icel)%nCellRay           = BasecellInfo(ibasecell)%nCellRay    
  CellInfo(icel)%EdgeCellIdx        = BasecellInfo(ibasecell)%EdgeCellIdx 
  CellInfo(icel)%GapType            = BasecellInfo(ibasecell)%GapType       
  CellInfo(icel)%FuelRad0           = BasecellInfo(ibasecell)%FuelRad0    
                                    
  CellInfo(icel)%lempty             = BasecellInfo(ibasecell)%lempty 
  CellInfo(icel)%luse               = BasecellInfo(ibasecell)%luse   
  CellInfo(icel)%lrect              = BasecellInfo(ibasecell)%lrect   
  CellInfo(icel)%lCCell             = BasecellInfo(ibasecell)%lCCell 
  CellInfo(icel)%lcentX             = BasecellInfo(ibasecell)%lcentX 
  CellInfo(icel)%lCentY             = BasecellInfo(ibasecell)%lCentY 
  CellInfo(icel)%lCentXY            = BasecellInfo(ibasecell)%lCentXY
  CellInfo(icel)%lgap               = BasecellInfo(ibasecell)%lgap   
  CellInfo(icel)%lcad               = BasecellInfo(ibasecell)%lcad   
  CellInfo(icel)%lCrCell            = BasecellInfo(ibasecell)%lCrCell
    
  CellInfo(icel)%FxrIdxSt           => BasecellInfo(ibasecell)%FxrIdxSt     
  CellInfo(icel)%nFsrInFxr          => BasecellInfo(ibasecell)%nFsrInFxr    
  CellInfo(icel)%MapFxr2FsrIdx      => BasecellInfo(ibasecell)%MapFxr2FsrIdx
  CellInfo(icel)%THCell             => BasecellInfo(ibasecell)%THCell       
  CellInfo(icel)%CadGeom            => BasecellInfo(ibasecell)%CadGeom      
  CellInfo(icel)%CrCell             => BasecellInfo(ibasecell)%CrCell       

  CellInfo(icel)%geom%nx            = BasecellInfo(ibasecell)%geom%nx     
  CellInfo(icel)%geom%ny            = BasecellInfo(ibasecell)%geom%ny     
  CellInfo(icel)%geom%cx            = BasecellInfo(ibasecell)%geom%cx     
  CellInfo(icel)%geom%cy            = BasecellInfo(ibasecell)%geom%cy     
  CellInfo(icel)%geom%lx            = BasecellInfo(ibasecell)%geom%lx     
  CellInfo(icel)%geom%ly            = BasecellInfo(ibasecell)%geom%ly     
  CellInfo(icel)%geom%x             = BasecellInfo(ibasecell)%geom%x      
  CellInfo(icel)%geom%y             = BasecellInfo(ibasecell)%geom%y      
  CellInfo(icel)%geom%inpnx         = BasecellInfo(ibasecell)%geom%inpnx  
  CellInfo(icel)%geom%inpny         = BasecellInfo(ibasecell)%geom%inpny  
  CellInfo(icel)%geom%inpdivx       = BasecellInfo(ibasecell)%geom%inpdivx
  CellInfo(icel)%geom%inpdivy       = BasecellInfo(ibasecell)%geom%inpdivy
  CellInfo(icel)%geom%inpdelx       = BasecellInfo(ibasecell)%geom%inpdelx
  CellInfo(icel)%geom%inpdely       = BasecellInfo(ibasecell)%geom%inpdely
                                    
  CellInfo(icel)%geom%lcircle       = BasecellInfo(ibasecell)%geom%lcircle     
  CellInfo(icel)%geom%lrect         = BasecellInfo(ibasecell)%geom%lrect       
  CellInfo(icel)%geom%lCCent        = BasecellInfo(ibasecell)%geom%lCCent
  
  CellInfo(icel)%geom%CCentType     = BasecellInfo(ibasecell)%geom%CCentType    
  CellInfo(icel)%geom%ncircle       = BasecellInfo(ibasecell)%geom%ncircle  
  
  CellInfo(icel)%geom%delx          => BasecellInfo(ibasecell)%geom%delx        
  CellInfo(icel)%geom%dely          => BasecellInfo(ibasecell)%geom%dely        
  CellInfo(icel)%geom%circle        => BasecellInfo(ibasecell)%geom%circle      

  IF (icel .GT. nCellType0) CYCLE
  
  IF (nTracerCntl%lSSPH) THEN
    IF (.NOT. BasecellInfo(ibasecell)%lRect) THEN
      ndata = BasecellInfo(ibasecell)%geom%inpn
      rr = BasecellInfo(ibasecell)%rr
      ndiv(1:ndata) = BasecellInfo(ibasecell)%geom%inpdiv(1:ndata)
      iReg = CellInfo(icel)%iReg0
      j = 0; nfueldiv=0
      IF (CellInfo(icel)%lfuel) THEN
        DO j=ndata-1,0,-1
          IF(nTracerCntl%lXsLib) THEN
            IF (Mixture(ireg(j+1))%lfuel.and..not.Mixture(ireg(j))%lfuel) EXIT 
          ELSE
            IF (MacXsBen(ireg(j+1))%lfuel.and..not.MacXsBen(ireg(j))%lfuel) EXIT
          ENDIF
        ENDDO
        j=j+1
        CellInfo(icel)%FuelRad0=rr(j)
        CellInfo(icel)%iefuel = j
        CellInfo(icel)%lhole = .NOT. Mixture(ireg(ndata))%lfuel
        DO j = ndata, 0, -1
          IF (Mixture(ireg(j))%lfuel) EXIT
        END DO
        CellInfo(icel)%ibfuel = j
      ELSEIF (CellInfo(icel)%lAIC) THEN
        DO j=ndata,0,-1
            IF (.not.Mixture(ireg(j))%lAIC) EXIT ! Based on the assumption that there's no hole in an AIC pin.
        ENDDO
        j=j+1
        CellInfo(icel)%FuelRad0=rr(j)  
        CellInfo(icel)%iefuel = j
        CellInfo(icel)%ibfuel = ndata
      ENDIF
      
      DO i=CellInfo(icel)%ibFuel,CellInfo(icel)%ieFuel,-1
        nfueldiv=nfueldiv+ndiv(i)
      ENDDO
      CellInfo(icel)%nfueldiv=nfueldiv
      IF (CellInfo(icel)%lAIC) THEN
        CALL calcAICCellSSPH(CellInfo,icel,ndata,iReg,rr,ndiv,j,nfueldiv)
      ELSE
        CALL calcCellSSPH(CellInfo,icel,ndata,iReg,rr,ndiv,CellInfo(icel)%ibFuel,CellInfo(icel)%ieFuel,nfueldiv, CellInfo(icel)%lhole)
      ENDIF
    ENDIF
  ENDIF
  
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE
  
SUBROUTINE CopyBaseCellGeom()
USE geom,     ONLY : nCellType, CellInfo, BaseCellInfo, nCellType0, nBaseCell0
USE PE_Mod,   ONLY : PE
USE OMP_LIB
IMPLICIT NONE

INTEGER :: icel, ibasecell

CALL omp_set_num_threads(PE%nThread)

!$OMP PARALLEL PRIVATE(ibasecell)
!$OMP DO SCHEDULE(GUIDED)
DO icel = 1, nCellType

  ibasecell = CellInfo(icel)%basecellstr  
  
  CellInfo(icel)%vol                => BasecellInfo(ibasecell)%vol
  
  CellInfo(icel)%geom%bdline        => BasecellInfo(ibasecell)%geom%bdline      
  CellInfo(icel)%geom%bdlineRange   => BasecellInfo(ibasecell)%geom%bdlineRange 
  CellInfo(icel)%geom%line          => BasecellInfo(ibasecell)%geom%line        
  
  CellInfo(icel)%geom%nbd          =  BasecellInfo(ibasecell)%geom%nbd        
  CellInfo(icel)%geom%nline        =  BasecellInfo(ibasecell)%geom%nline             
    
END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE
