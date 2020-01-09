SUBROUTINE WriteFxr(Core, FmInfo, PE)
USE PARAM
USE typedef, ONLY : CoreInfo_Type, FmInfo_Type, FxrInfo_Type, PE_TYPE
USE files,     ONLY : io9, io10, io11, io13, caseid, localfn
USE Depl_Mod,  ONLY : DeplCntl

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_TYPE) :: PE

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: niso, nFxr, nchar
INTEGER :: ifxr, iz

CHARACTER(256) :: fn, fn0,caseid0

nFxr = Core%nCoreFxr
Fxr => FmInfo%Fxr

DO i = 1, 80
  if(caseid(i:i) .EQ. '') EXIT
ENDDO
i = i-1
caseid0(1:i) = caseid(1:i)
IF(DeplCntl%NowStep .LT. 10) THEN
  WRITE(caseid0(i+1:256),'(A, A, I1)') '_cycle', '00', DeplCntl%NowStep
ELSEIF(DeplCntl%NowStep .LT. 100) THEN
  WRITE(caseid0(i+1:256),'(A, A, I2)') '_cycle', '0', DeplCntl%NowStep
ELSE
  WRITE(caseid0(i+1:256),'(A, I3)') '_cycle', DeplCntl%NowStep
ENDIF
i=i+9; nchar = i
DO iz = PE%myzb, PE%myze
  IF(iz .LT. 10) THEN
    WRITE(fn,'(A,A,A,I1,A)') caseid0(1:nchar), '_pln','00', iz, '.iso'
    !WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i), '_pln','00', iz
  ELSEIF(iz .LT.100) THEN
    WRITE(fn,'(A,A,A,I2,A)') caseid0(1:nchar), '_pln','0', iz, '.iso'
    !WRITE(fn_base,'(A,A,A,I2,A)') caseid0(1:i), '_pln','0', iz
  ELSE
    WRITE(fn,'(A,A,I3,A)') caseid0(1:nchar), '_pln', iz, '.iso'
    !WRITE(fn_base,'(A,A,I3,A)') caseid0(1:i), '_pln', iz
  ENDIF
  !WRITE(fn,'(A,A,A,I1,A)') caseid0(1:i), '.vtk'
  !WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i)
  OPEN(unit=io13, file = fn, status = 'replace')
  DO i = 1, nFxr
    WRITE(io13, '(I15, I15, F20.7)') i, Fxr(i, iz)%niso, Fxr(i, iz)%burnup
    DO j = 1, Fxr(i, iz)%niso
      WRITE(io13, '(I10, e40.25)') Fxr(i, iz)%idiso(j),Fxr(i, iz)%pnum(j)
    ENDDO
  ENDDO
  CLOSE(io13)
ENDDO

NULLIFY(FXR)
END SUBROUTINE


SUBROUTINE ReadFxr(Core, FmInfo, PE)
USE PARAM
USE typedef, ONLY : CoreInfo_Type, FmInfo_Type, FxrInfo_Type, PE_TYPE
USE files,     ONLY : io9, io10, io11, io13, caseid, localfn
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_TYPE) :: PE

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: niso, nFxr, nchar
INTEGER :: ifxr, iz
INTEGER :: i, j, k
REAL :: burnup


CHARACTER(256) :: fn, fn0,caseid0
nFxr = Core%nCoreFxr
Fxr => FmInfo%Fxr
caseid0 = '../TDSG/burnup/pc2q_47g_1_cycle030'

DO i = 1, 80
  if(caseid0(i:i) .EQ. '') EXIT
ENDDO
nchar = i - 1
DO iz = PE%myzb, PE%myze
  IF(iz .LT. 10) THEN
    WRITE(fn,'(A,A,A,I1,A)') caseid0(1:nchar), '_pln','00', iz, '.iso'
    !WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i), '_pln','00', iz
  ELSEIF(iz .LT.100) THEN
    WRITE(fn,'(A,A,A,I2,A)') caseid0(1:nchar), '_pln','0', iz, '.iso'
    !WRITE(fn_base,'(A,A,A,I2,A)') caseid0(1:i), '_pln','0', iz
  ELSE
    WRITE(fn,'(A,A,I3,A)') caseid0(1:nchar), '_pln', iz, '.iso'
    !WRITE(fn_base,'(A,A,I3,A)') caseid0(1:i), '_pln', iz
  ENDIF

  OPEN(unit=io13, file = fn, status = 'old')
  DO i = 1, nFxr
    READ(io13, *) ifxr, niso, burnup
    IF(niso == 0) CYCLE
    DO j = 1, niso
      READ(io13, *) Fxr(i, iz)%idiso(j), Fxr(i, iz)%pnum(j)
      IF(Fxr(i, iz)%idiso(j) .EQ. 94240) THEN
        Fxr(i, iz)%pnum(j) = 1.5_8 * Fxr(i, iz)%pnum(j)
      ENDIF
    ENDDO
    Fxr(i,iz)%niso = niso; Fxr(i, iz)%burnup = burnup
  ENDDO
  CLOSE(io13)
ENDDO
NULLIFY(Fxr)
END SUBROUTINE

SUBROUTINE calnum(Core, FmInfo, PE)
USE PARAM
USE typedef, ONLY : CoreInfo_Type, FmInfo_Type, FxrInfo_Type, PE_TYPE
USE files,     ONLY : io9, io10, io11, io13, caseid, localfn
USE Depl_Mod,  ONLY : DeplCntl

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_TYPE) :: PE

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: niso, nFxr, nchar
INTEGER :: ifxr, iz
real :: psum

CHARACTER(256) :: fn, fn0,caseid0

nFxr = Core%nCoreFxr
Fxr => FmInfo%Fxr

DO i = 1, 80
  if(caseid(i:i) .EQ. '') EXIT
ENDDO
i = i-1
caseid0(1:i) = caseid(1:i)
DO istep = 1, 50
  DO i = 1, 80
    if(caseid(i:i) .EQ. '') EXIT
  ENDDO
  i = i-1
  caseid0(1:i) = caseid(1:i)
  IF(istep .LT. 10) THEN
    WRITE(caseid0(i+1:256),'(A, A, I1)') '_cycle', '00',istep
  ELSEIF(istep .LT. 100) THEN
    WRITE(caseid0(i+1:256),'(A, A, I2)') '_cycle', '0',istep
  ELSE
    WRITE(caseid0(i+1:256),'(A, I3)') '_cycle',istep
  ENDIF
  i=i+9; nchar = i
  DO iz = PE%myzb, PE%myze
    IF(iz .LT. 10) THEN
      WRITE(fn,'(A,A,A,I1,A)') caseid0(1:nchar), '_pln','00', iz, '.iso'
      !WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i), '_pln','00', iz
    ELSEIF(iz .LT.100) THEN
      WRITE(fn,'(A,A,A,I2,A)') caseid0(1:nchar), '_pln','0', iz, '.iso'
      !WRITE(fn_base,'(A,A,A,I2,A)') caseid0(1:i), '_pln','0', iz
    ELSE
      WRITE(fn,'(A,A,I3,A)') caseid0(1:nchar), '_pln', iz, '.iso'
      !WRITE(fn_base,'(A,A,I3,A)') caseid0(1:i), '_pln', iz
    ENDIF
    !WRITE(fn,'(A,A,A,I1,A)') caseid0(1:i), '.vtk'
    !WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i)
    OPEN(unit=io13, file = fn, status = 'old')
    psum = 0
    DO k = 1, nFxr
      READ(io13, *) ifxr, niso    
      IF(niso == 0) CYCLE
      DO j = 1, niso
        READ(io13, *) Fxr(k, iz)%idiso(j), Fxr(k, iz)%pnum(j)
        IF(Fxr(k, iz)%idiso(j) .EQ. 94238) THEN
          psum = psum +  Fxr(k, iz)%pnum(j)/9._8
        ENDIF
      ENDDO
      Fxr(k,iz)%niso = niso   
    ENDDO
    CLOSE(io13)
    WRITE(101, *) istep, psum
  ENDDO
ENDDO
NULLIFY(FXR)
stop
END SUBROUTINE

