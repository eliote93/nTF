MODULE MCP_Restart_Main
  CONTAINS
  SUBROUTINE MCP_Initialize()
    USE param
    USE files,                  ONLY : io8,             Caseid
    USE ioutil,                 ONLY : message,         IFnumeric,          Terminate
    USE PE_mod,                 ONLY : PE
    USE nuclidmap_mod,          ONLY : lExistIsotope
    USE MCP_Util
    USE Material_mod,           ONLY : Mixture

    IMPLICIT NONE
    LOGICAL :: master
    INTEGER :: iMCP
    INTEGER :: i, j, k, l, m
    INTEGER, PARAMETER :: InputMCP = 85
    CHARACTER(512) :: Oneline, String
    INTEGER :: nFXR
    INTEGER :: Macro_region, Micro_region, nMix, nMix_nT
    CHARACTER(10) :: nMix_Char, Temperature_Char
    REAL :: Temperature, vol
    CHARACTER(10) :: Dummy1, Dummy2, Dummy3
    CHARACTER(4) :: Dummy4
    INTEGER :: isotope_number, idx, nMixInput, iMixSearch, nIsoInput, Macro_Index, Macro_Index2, Micro_Index
    REAL :: isotope_density
    LOGICAL :: lExist_isotope
    INTEGER :: Ring_Start, Ring_End, nRing, Ring_Start_Multi
    LOGICAL :: lOctant, lQuarter
    INTEGER :: index_file, MCP_no, Mat_no
    CHARACTER(512) :: fn
    INTEGER :: nIsotope, Idx_Iso, Idx_Cell, Idx_Cell2, Idx_Pin, Cell_line, MCP_CELL_MAT_INDEX, MCP_ASEM_MAT_INDEX, Num_Start
    LOGICAL :: lDepl, lFuel, lres
    CHARACTER(1) :: MixLogic
    CHARACTER(1) :: c1, c2, c3, c4, c5
    REAL :: r1, r2
    INTEGER :: i_sol, i_mod
    LOGICAL :: l_1000
    INTEGER :: nPlane, index_MCP_start, index_MCP_end, iPlane, iPin, nCell

    TYPE(TYPE_MCP_Info), ALLOCATABLE :: MCP_Info(:)

    master = PE%master
    nMixInput = SIZE(Mixture)

    nPlane = MOD(nMCP, nMCP_Plane)
    IF(nPlane .NE. 0) THEN
      mesg = 'nMCP is not matched with nMCP_Plane'
      IF(MASTER) CALL Terminate(mesg)
    END IF
    nPlane = nMCP/nMCP_Plane
    IF (MASTER) THEN
      mesg = 'MCP Restart Initilization...'
      CALL message(io8, TRUE, TRUE, mesg)
    ENDIF

    ALLOCATE(MCP_Info(1:nMCP))
    !Scannint MCP file real quick
    DO iMCP = 1, nMCP
      j = 0
      OPEN(UNIT = InputMCP, File = MCP_Path(iMCP), STATUS='OLD', ACTION = 'READ', POSITION = 'REWIND')
      DO
        READ(InputMCP, *, END = 1000) Oneline
        String = 'Macro-region:'
        IF(Oneline .EQ. String) THEN
          j = j + 1
        END IF
      END DO
      1000 CONTINUE

      MCP_Info(iMCP)%nFXR = j
      IF(j .GT. 0) ALLOCATE(MCP_Info(iMCP)%IntInfo(1:j))
      IF(j .GT. 0) ALLOCATE(MCP_Info(iMCP)%Mixture(1:j))
      CLOSE(InputMCP)
    END DO

    DO iMCP = 1, nMCP
      OPEN(UNIT = InputMCP, File = MCP_Path(iMCP), STATUS='OLD', ACTION = 'READ', POSITION = 'REWIND')
      DO
        READ(InputMCP, *, END = 1002) Oneline
        String = 'PWR'
        IF(Oneline .EQ. String) THEN
          BACKSPACE(InputMCP)
          READ(InputMCP, *) Oneline, MCP_Info(iMCP)%nCell, r1, r2, c1, c2, c3, c4, MCP_Info(iMCP)%Assem_Type
          EXIT
        END IF
      END DO
      1002 CONTINUE
      CLOSE(InputMCP)

      IF(MOD(MCP_Info(iMCP)%nCell, 2) .EQ. 0) THEN
        MCP_Info(iMCP)%nCell45 = MCP_Info(iMCP)%nCell/2
      ELSE IF (MOD(MCP_Info(iMCP)%nCell, 2) .EQ. 1) THEN
        MCP_Info(iMCP)%nCell45 = MCP_Info(iMCP)%nCell/2 + 1
      END IF

      IF(MCP_Info(iMCP)%Assem_Type .EQ. 8) THEN
        ALLOCATE(MCP_info(iMCP)%Micro_Number_Start(1:MCP_Info(iMCP)%nCell45))
        ALLOCATE(MCP_info(iMCP)%Micro_Number_End(1:MCP_Info(iMCP)%nCell45))
      ELSEIF(MCP_Info(iMCP)%Assem_Type .EQ. 1) THEN
        ALLOCATE(MCP_info(iMCP)%Micro_Number_Start(1:MCP_Info(iMCP)%nCell))
        ALLOCATE(MCP_info(iMCP)%Micro_Number_End(1:MCP_Info(iMCP)%nCell))
      ELSE
        WRITE(mesg, *) 'This MCP assembly type is not for PWR'
        CALL Terminate(mesg)
      END IF
    END DO

    DO iMCP = 1, nMCP
      OPEN(UNIT = InputMCP, File = MCP_Path(iMCP), STATUS='OLD', ACTION = 'READ', POSITION = 'REWIND')
      DO
        READ(InputMCP, *, END = 1001) Oneline
        String = 'Macro-region'
        IF(Oneline .EQ. String) THEN
          READ(InputMCP, *) Oneline
          READ(Oneline,'(i6)') MCP_Info(iMCP)%nMaxMacro
          IF (MCP_Info(iMCP)%Assem_Type .EQ. 1) MCP_Info(iMCP)%nMaxMacro = MCP_Info(iMCP)%nMaxMacro + MCP_Info(iMCP)%nCell - 1
          EXIT
        END IF
      END DO
      1001 CONTINUE
      CLOSE(InputMCP)

      Macro_Index = MCP_Info(iMCP)%nMaxMacro
      Macro_Index2 = MCP_Info(iMCP)%nMaxMacro
      IF(MCP_Info(iMCP)%Assem_Type .EQ. 8) THEN
        DO k = 1, MCP_Info(iMCP)%nCell45
          MCP_Info(iMCP)%Micro_Number_Start(k) = Macro_Index
          MCP_Info(iMCP)%Micro_Number_End(k) = Macro_Index2
          Macro_Index = Macro_Index - k - 1
          Macro_Index2 = Macro_Index2 - k
        END DO
        IF(MCP_Info(iMCP)%Micro_Number_Start(MCP_Info(iMCP)%nCell45) .NE. 1) PRINT*, 'Error at Micro Region Index'

        Macro_Index = MCP_Info(iMCP)%nMaxMacro
        ALLOCATE(MCP_info(iMCP)%Micro_Char(1:2*MCP_Info(iMCP)%nCell45))
        ALLOCATE(MCP_info(iMCP)%Micro_Start(1:Macro_Index))
        ALLOCATE(MCP_info(iMCP)%Micro_End(1:Macro_Index))
      ELSEIF(MCP_Info(iMCP)%Assem_Type .EQ. 1) THEN
        DO k = 1, MCP_Info(iMCP)%nCell
          MCP_Info(iMCP)%Micro_Number_Start(k) = Macro_Index - MCP_Info(iMCP)%nCell + 1
          MCP_Info(iMCP)%Micro_Number_End(k) = Macro_Index
          Macro_Index = Macro_Index - MCP_Info(iMCP)%nCell
        END DO
        IF(MCP_Info(iMCP)%Micro_Number_Start(MCP_Info(iMCP)%nCell) .NE. 1) PRINT*, 'Error at Micro Region Index'
        Macro_Index = MCP_Info(iMCP)%nMaxMacro
        ALLOCATE(MCP_info(iMCP)%Micro_Char(1:3*8))
        ALLOCATE(MCP_info(iMCP)%Micro_Start(1:Macro_Index))
        ALLOCATE(MCP_info(iMCP)%Micro_End(1:Macro_Index))
      END IF
    END DO

    DO iMCP = 1, nMCP
      OPEN(UNIT = InputMCP, File = MCP_Path(iMCP), STATUS='OLD', ACTION = 'READ', POSITION = 'REWIND')
      IF(MCP_Info(iMCP)%Assem_Type .EQ. 8) THEN
        DO
          READ(InputMCP, *, END = 1003) Oneline
          String = 'Micro-region'
          IF(Oneline .EQ. String) THEN
            DO k = 1, MCP_Info(iMCP)%nCell45
              READ(InputMCP, '(a512)') Oneline
              READ(Oneline, *) (MCP_Info(iMCP)%Micro_Char(l), l=1,2*k)
              Macro_Index = MCP_Info(iMCP)%Micro_Number_Start(k)
              m = 0
              DO l = 1, 2*k
                m = m + 1
                IF(MOD(m, 2) .EQ. 1) THEN
                  MCP_Info(iMCP)%Micro_Char(m)(LEN_TRIM(MCP_Info(iMCP)%Micro_Char(m)):LEN_TRIM(MCP_Info(iMCP)%Micro_Char(m))) = ''
                  READ(MCP_Info(iMCP)%Micro_Char(m),*) MCP_Info(iMCP)%Micro_Start(Macro_Index)
                ELSE
                  READ(MCP_Info(iMCP)%Micro_Char(m),*) MCP_Info(iMCP)%Micro_End(Macro_Index)
                  Macro_Index = Macro_Index + 1
                END IF
              END DO
            END DO
            EXIT
          END IF
        END DO
      ELSEIF(MCP_Info(iMCP)%Assem_Type .EQ. 1) THEN
        DO
          READ(InputMCP, *, END = 1003) Oneline
          String = 'Macro-micro'
          i_sol = MCP_Info(iMCP)%nMaxMacro / 8
          i_sol = i_sol + 1
          i_mod = mod(MCP_Info(iMCP)%nMaxMacro, 8)
          l_1000 = .FALSE.
          Macro_Index = 1
          IF(Oneline .EQ. String) THEN
            READ(InputMCP, '(a512)') Oneline
            DO k = 1, i_sol - 1
              IF(.NOT. l_1000) THEN
                READ(InputMCP, '(a512)') Oneline
                READ(InputMCP, '(a512)') Oneline
                READ(Oneline, *) c1, c2, (MCP_Info(iMCP)%Micro_Char(l), l=1,3*8)
                DO m = 1, 8
                  READ(MCP_Info(iMCP)%Micro_Char(3*m-2), *) MCP_Info(iMCP)%Micro_Start(Macro_Index)
                  READ(MCP_Info(iMCP)%Micro_Char(3*m), *) MCP_Info(iMCP)%Micro_End(Macro_Index)
                  IF(MCP_Info(iMCP)%Micro_End(Macro_Index) .GE. 999) l_1000 = .TRUE.
                  Macro_Index = Macro_Index + 1
                END DO
              ELSEIF(l_1000) THEN
                READ(InputMCP, '(a512)') Oneline
                READ(InputMCP, '(a512)') Oneline
                READ(Oneline, *) c1, (MCP_Info(iMCP)%Micro_Char(l), l = 1,3*8)
                DO m = 1, 8
                  IF(m .EQ. 1) MCP_Info(iMCP)%Micro_Char(m)(1:1) = ''
                  READ(MCP_Info(iMCP)%Micro_Char(3*m-2), *) MCP_Info(iMCP)%Micro_Start(Macro_Index)
                  READ(MCP_Info(iMCP)%Micro_Char(3*m), *) MCP_Info(iMCP)%Micro_End(Macro_Index)
                  Macro_Index = Macro_Index + 1
                END DO
              END IF
            END DO
            IF(i_mod .NE. 0) THEN
              READ(InputMCP, '(a512)') Oneline
              READ(InputMCP, '(a512)') Oneline
              IF(l_1000) THEN
                READ(Oneline, *) c1, (MCP_Info(iMCP)%Micro_Char(l), l = 1,3*i_mod)
                DO m = 1, i_mod
                  IF(m .EQ. 1) MCP_Info(iMCP)%Micro_Char(m)(1:1) = ''
                  READ(MCP_Info(iMCP)%Micro_Char(3*m-2), *) MCP_Info(iMCP)%Micro_Start(Macro_Index)
                  READ(MCP_Info(iMCP)%Micro_Char(3*m), *) MCP_Info(iMCP)%Micro_End(Macro_Index)
                  Macro_Index = Macro_Index + 1
                END DO
              ELSEIF(.NOT. l_1000) THEN
                READ(Oneline, *) c1, c2, (MCP_Info(iMCP)%Micro_Char(l), l = 1,3*i_mod)
                DO m = 1, i_mod
                  READ(MCP_Info(iMCP)%Micro_Char(3*m-2), *) MCP_Info(iMCP)%Micro_Start(Macro_Index)
                  READ(MCP_Info(iMCP)%Micro_Char(3*m), *) MCP_Info(iMCP)%Micro_End(Macro_Index)
                  Macro_Index = Macro_Index + 1
                END DO
              END IF
            END IF
          END IF
        END DO
      END IF
      1003 CONTINUE
      CLOSE(InputMCP)
      MCP_Info(iMCP)%nMaxMicro = MCP_Info(iMCP)%Micro_End(MCP_Info(iMCP)%nMaxMacro)
    END DO

    !Save the Isotope Number Densities in MCP File
    DO iMCP = 1, nMCP
      j = 0
      OPEN(UNIT = InputMCP, File = MCP_Path(iMCP), STATUS='OLD', ACTION = 'READ', POSITION = 'REWIND')
      DO
        READ(InputMCP, *, END = 10001) Oneline
        String = 'Macro-region:'
        IF(Oneline .EQ. String) THEN
          j = j + 1
          BACKSPACE(InputMCP)
          READ(InputMCP, *) Dummy1, Macro_region, Dummy2, Micro_region, Dummy3, Vol, Dummy4, Temperature_Char, nMix_Char
          MCP_Info(iMCP)%IntInfo(j)%Macro_region = Macro_region; MCP_Info(iMCP)%IntInfo(j)%Micro_region = Micro_region
          MCP_Info(iMCP)%IntInfo(j)%vol = vol
          READ(Temperature_Char, '(F7.3)') Temperature
          MCP_Info(iMCP)%Mixture(j)%lempty = .FALSE.
          MCP_Info(iMCP)%Mixture(j)%lMox = .FALSE.
          MCP_Info(iMCP)%Mixture(j)%name = dummy4
          MCP_Info(iMCP)%Mixture(j)%temp = Temperature ! Unit = Kelvin
          idx = 0
          IF(IFnumeric(nMix_Char)) THEN
            READ(nMix_Char, '(i4)') nMix
            IF(nMix .GE. 308) THEN
              IF(nMix .EQ. 308) nMix_nT = nMix - 146            !
              IF(nMix .EQ. 328) nMix_nT = nMix - 146 - 4        !
              MCP_Info(iMCP)%Mixture(j)%niso = nMix_nT
              MCP_Info(iMCP)%Mixture(j)%lRes = .TRUE.
              MCP_Info(iMCP)%Mixture(j)%lFuel = .TRUE.
              MCP_Info(iMCP)%Mixture(j)%lDepl = .TRUE.
              ALLOCATE(MCP_Info(iMCP)%Mixture(j)%idiso(1:nMix_nT))
              ALLOCATE(MCP_Info(iMCP)%Mixture(j)%fweig(1:nMix_nT))
              ALLOCATE(MCP_Info(iMCP)%Mixture(j)%pnum(1:nMix_nT))
              DO k = 1, nMix
                READ(InputMCP, *) isotope_number, Dummy1, Dummy2, isotope_density
                lExist_isotope = .FALSE.
                CALL isotope_number_chk(isotope_number, lExist_isotope, iMCP, j)
                IF(lExist_isotope) THEN
                  !isotope nubmering check for nTRACER library
                  IF(.NOT. lExistIsotope(isotope_number)) THEN
                    WRITE(mesg, *) 'This MCP converted nulide is not included in the library', iMCP, isotope_number
                    CALL Terminate(mesg)
                  END IF
                  idx = idx + 1
                  MCP_Info(iMCP)%Mixture(j)%idiso(idx) = isotope_number
                  MCP_Info(iMCP)%Mixture(j)%fweig(idx) = isotope_density
                  MCP_Info(iMCP)%Mixture(j)%pnum(idx) = isotope_density
                END IF
              END DO
            ELSE
              SELECT CASE(Dummy4)
                CASE('HEL ')
                  Dummy4 = 'AIR '
                CASE('CAN ')
                  Dummy4 = 'CLD '
                CASE('BOX ')
                  Dummy4 = 'CLD '
                CASE('COO ')
                  Dummy4 = 'MOD '
                CASE('AIR ')
                  Dummy4 = 'AIR '
              END SELECT
              MCP_Info(iMCP)%Mixture(j)%name = Dummy4
              DO k = 1, nMixInput
                IF(Mixture(k)%name .EQ. Dummy4) THEN
                  iMixSearch = k
                END IF
              END DO
              nIsoInput = Mixture(iMixSearch)%niso
              MCP_Info(iMCP)%Mixture(j)%niso = nIsoInput
              MCP_Info(iMCP)%Mixture(j)%lRes = Mixture(iMixSearch)%lRes
              MCP_Info(iMCP)%Mixture(j)%lFuel = Mixture(iMixSearch)%lFuel
              IF(MCP_Info(iMCP)%Mixture(j)%lFuel) PRINT*, 'MCP Isotope Match Fail : Non-Fuel Region'
              MCP_Info(iMCP)%Mixture(j)%lDepl = Mixture(iMixSearch)%lDepl
              ALLOCATE(MCP_Info(iMCP)%Mixture(j)%idiso(1:nIsoInput))
              ALLOCATE(MCP_Info(iMCP)%Mixture(j)%fweig(1:nIsoInput))
              ALLOCATE(MCP_Info(iMCP)%Mixture(j)%pnum(1:nIsoInput))
              DO k = 1, nIsoInput
                MCP_Info(iMCP)%Mixture(j)%idiso(k) = Mixture(iMixSearch)%idiso(k)
                MCP_Info(iMCP)%Mixture(j)%fweig(k) = Mixture(iMixSearch)%fweig(k)
                MCP_Info(iMCP)%Mixture(j)%pnum(k) = Mixture(iMixSearch)%pnum(k)
              END DO
            END IF
          ELSE
            SELECT CASE(Dummy4)
              CASE('HEL ')
                Dummy4 = 'AIR '
              CASE('CAN ')
                Dummy4 = 'CLD '
              CASE('BOX ')
                Dummy4 = 'CLD '
              CASE('COO ')
                Dummy4 = 'MOD '
              CASE('AIR ')
                Dummy4 = 'AIR '
            END SELECT
            MCP_Info(iMCP)%Mixture(j)%name = Dummy4
            DO k = 1, nMixInput
              IF(Mixture(k)%name .EQ. Dummy4) THEN
                iMixSearch = k
              END IF
            END DO
            nIsoInput = Mixture(iMixSearch)%niso
            MCP_Info(iMCP)%Mixture(j)%niso = nIsoInput
            MCP_Info(iMCP)%Mixture(j)%lRes = Mixture(iMixSearch)%lRes
            MCP_Info(iMCP)%Mixture(j)%lFuel = Mixture(iMixSearch)%lFuel
            IF(MCP_Info(iMCP)%Mixture(j)%lFuel) PRINT*, 'MCP Isotope Match Fail : Non-Fuel Region'
            MCP_Info(iMCP)%Mixture(j)%lDepl = Mixture(iMixSearch)%lDepl
            ALLOCATE(MCP_Info(iMCP)%Mixture(j)%idiso(1:nIsoInput))
            ALLOCATE(MCP_Info(iMCP)%Mixture(j)%fweig(1:nIsoInput))
            ALLOCATE(MCP_Info(iMCP)%Mixture(j)%pnum(1:nIsoInput))
            DO k = 1, nIsoInput
              MCP_Info(iMCP)%Mixture(j)%idiso(k) = Mixture(iMixSearch)%idiso(k)
              MCP_Info(iMCP)%Mixture(j)%fweig(k) = Mixture(iMixSearch)%fweig(k)
              MCP_Info(iMCP)%Mixture(j)%pnum(k) = Mixture(iMixSearch)%pnum(k)
            END DO
          ENDIF
        END IF
      END DO
      10001 CONTINUE
      CLOSE(InputMCP)
    END DO

    !Volume Calculation -> Radius and Number of Cell Data Collection
    DO iMCP = 1, nMCP
      idx = 0
      Macro_Index = MCP_Info(iMCP)%nMaxMacro
      Micro_Index = MCP_Info(iMCP)%nMaxMicro
      ALLOCATE(MCP_Info(iMCP)%Radius(1:Macro_Index))
      DO i = 1, Macro_Index
        nRing = MCP_Info(iMCP)%Micro_End(i) - MCP_Info(iMCP)%Micro_Start(i) + 1
        MCP_Info(iMCP)%Radius(i)%nRing = nRing
        ALLOCATE(MCP_Info(iMCP)%Radius(i)%Ring_Volume(1:nRing))
        ALLOCATE(MCP_Info(iMCP)%Radius(i)%Ring_Radius(1:nRing-1))
        Ring_Start = MCP_Info(iMCP)%Micro_Start(i)
        Ring_End = MCP_Info(iMCP)%Micro_End(i)

        IF(MCP_Info(iMCP)%Assem_Type .EQ. 8) THEN
          !Scanning
          lOctant = .FALSE.
          lQuarter = .FALSE.
          DO k = 1, MCP_Info(iMCP)%nCell45
            IF(MCP_Info(iMCP)%Micro_Number_Start(k) .EQ. i) lQuarter = .TRUE.
            IF(MCP_Info(iMCP)%Micro_Number_End(k) .EQ. i) lQuarter = .TRUE.
          END DO
          IF(MCP_Info(iMCP)%nMaxMacro .EQ. i) lOctant = .TRUE.

          IF(lQuarter .AND. lOctant) THEN
            DO j = 1, nRing
              idx = idx + 1
              MCP_Info(iMCP)%Radius(i)%Ring_Volume(j) = MCP_Info(iMCP)%intInfo(idx)%vol * 1.25_8 * 8._8 / 10._8
            END DO
          ELSE IF(lQuarter .AND. .NOT. lOctant) THEN
            DO j = 1, nRing
              idx = idx + 1
              MCP_Info(iMCP)%Radius(i)%Ring_Volume(j) = MCP_Info(iMCP)%intInfo(idx)%vol * 1.25_8 * 2._8 / 10._8
            END DO
          ELSE
            DO j = 1, nRing
              idx = idx + 1
              MCP_Info(iMCP)%Radius(i)%Ring_volume(j) = MCP_Info(iMCP)%intInfo(idx)%vol * 1.25_8 / 10._8
            END DO
          END IF
          vol = 0
          DO j = 1, nRing - 1
            vol = vol + MCP_Info(iMCP)%Radius(i)%Ring_Volume(j)
            MCP_Info(iMCP)%Radius(i)%Ring_Radius(j) = SQRT(vol/PI)
          END DO
        ELSEIF(MCP_Info(iMCP)%Assem_Type .EQ. 1) THEN
          DO j = 1, nRing
            idx = idx + 1
            MCP_Info(iMCP)%Radius(i)%Ring_Volume(j) = MCP_Info(iMCP)%intInfo(idx)%vol
          END DO
          vol = 0
          DO j = 1, nRing - 1
            vol = vol + MCP_Info(iMCP)%Radius(i)%Ring_Volume(j)
            MCP_Info(iMCP)%Radius(i)%Ring_Radius(j) = SQRT(vol/PI)
          END DO
        END IF
      END DO
    END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Write Material BLOCK for Restart
    index_file = 100
    MCP_no = 0
    Idx_Iso = 0
    Idx_Cell = 0
    MCP_CELL_MAT_INDEX = 0
    DO iMCP = 1, nMCP
      !MATERIAL INPUT WRITE
      IF(iMCP .LT. 10) THEN
        WRITE(fn, '(A,A,A,i1,A)') 'MCP_',TRIM(CASEID),'_MATERIAL_000',iMCP,'.include'
      ELSE IF(iMCP .GE. 10 .AND. iMCP .LT. 100) THEN
        WRITE(fn, '(A,A,A,i2,A)') 'MCP_',TRIM(CASEID),'_MATERIAL_00',iMCP,'.include'
      ELSE IF(iMCP .GE. 100 .AND. iMCP .LT. 1000) THEN
        WRITE(fn, '(A,A,A,i3,A)') 'MCP_',TRIM(CASEID),'_MATERIAL_0',iMCP,'.include'
      ELSE IF(iMCP .GE. 1000) THEN
        WRITE(fn, '(A,A,A,i4,A)') 'MCP_',TRIM(CASEID),'_MATERIAL_',iMCP,'.include'
      END IF
      OPEN(UNIT = index_file, File = fn, STATUS = 'UNKNOWN', ACTION = 'WRITE', POSITION = 'REWIND')
      WRITE(Oneline, '(A,i4)') '!This Material Include file is from MCP Number_', iMCP
      WRITE(index_file,'(A)') TRIM(Oneline)
      WRITE(Oneline, '(A)') '!Density of Material will be automatically calculated durint nTRACER cal'
      WRITE(index_file,'(A)') TRIM(Oneline)
      Micro_index = MCP_Info(iMCP)%nMaxMicro
      DO i = 1, Micro_index
        nIsotope = MCP_Info(iMCP)%Mixture(i)%niso
        lDepl = .FALSE.; lFuel = .FALSE.; lres = .FALSE.
        Idx_Iso = Idx_Iso + 1
        MixLogic = '0'
        DO j = 1, nIsotope
          IF(j .EQ. 1) THEN
            lDepl = MCP_Info(iMCP)%Mixture(i)%lDepl;
            lFuel = MCP_Info(iMCP)%Mixture(i)%lFuel;
            lRes = MCP_Info(iMCP)%Mixture(i)%lRes;
            IF(lRes) MixLogic = '1'
            IF(lDepl) MixLogic = '3'
            IF(lFuel) MixLogic = '2'
            WRITE(Oneline, '(A,2X,I5,7X,A,1X,A,1X,F10.3,1X,F10.2,3X,A,3X,I5,3X,1PE15.6)') &
              '  mixture',Idx_Iso,MCP_Info(iMCP)%Mixture(i)%name,MixLogic,10.000,MCP_Info(iMCP)%Mixture(i)%Temp-ckelvin,'/',MCP_Info(iMCP)%Mixture(i)%idiso(j), MCP_Info(iMCP)%Mixture(i)%pnum(j)
            WRITE(index_file, '(A)') TRIM(Oneline)
          ELSE
            WRITE(Oneline, '(58X,I5,3X,1PE15.6)') MCP_Info(iMCP)%Mixture(i)%idiso(j), MCP_Info(iMCP)%Mixture(i)%pnum(j)
            WRITE(index_file, '(A)') TRIM(Oneline)
          END IF
          Oneline = ''
        END DO
      END DO
      CLOSE(index_file)
      !MATERIAL WRITE DONE
    END DO

    DO iMCP = 1, nMCP
      !CELL INPUT WRITE
      IF(iMCP .LT. 10) THEN
        WRITE(fn, '(A,A,A,i1,A)') 'MCP_',TRIM(CASEID),'_CELL_000',iMCP,'.include'
      ELSE IF(iMCP .GE. 10 .AND. iMCP .LT. 100) THEN
        WRITE(fn, '(A,A,A,i2,A)') 'MCP_',TRIM(CASEID),'_CELL_00',iMCP,'.include'
      ELSE IF(iMCP .GE. 100 .AND. iMCP .LT. 1000) THEN
        WRITE(fn, '(A,A,A,i3,A)') 'MCP_',TRIM(CASEID),'_CELL_0',iMCP,'.include'
      ELSE IF(iMCP .GE. 1000) THEN
        WRITE(fn, '(A,A,A,i4,A)') 'MCP_',TRIM(CASEID),'_CELL_',iMCP,'.include'
      END IF
      OPEN(UNIT = index_file, File = fn, STATUS = 'UNKNOWN', ACTION = 'WRITE', POSITION = 'REWIND')
      WRITE(Oneline, '(A,i4)') '!This Cell Include file is from MCP Number_', iMCP
      WRITE(index_file,'(A)') TRIM(Oneline)
      Macro_index = MCP_Info(iMCP)%nMaxMacro
      DO i = 1, Macro_index
        Oneline = ''
        Idx_Cell = Idx_Cell + 1
        WRITE(Oneline(3:6),'(A)') 'cell'
        WRITE(Oneline(7:11),'(i4)') Idx_Cell
        Cell_Line = 12
        nFXR = MCP_Info(iMCP)%Radius(i)%nRing
        DO j = 1, nFXR - 1
          !Cell_Line = Cell_Line + 1
          WRITE(Oneline(Cell_Line:Cell_Line+5), '(F6.4)') MCP_info(iMCP)%Radius(i)%Ring_Radius(j)
          Cell_Line = Cell_Line + 7
        END DO
        !Cell_Line = Cell_Line + 1
        WRITE(Oneline(Cell_Line:Cell_Line), '(A)') '/'
        Cell_Line = Cell_Line + 2

        Ring_Start = MCP_Info(iMCP)%Micro_Start(i) + MCP_CELL_MAT_INDEX

        DO j = 1, nFXR
          !Cell_Line = Cell_Line + 1
          WRITE(Oneline(Cell_Line:Cell_Line+4), '(i5)') Ring_Start
          Ring_Start = Ring_Start + 1
          Cell_Line = Cell_Line + 7
        END DO

        !Cell_Line = Cell_Line + 1
        WRITE(Oneline(Cell_Line:Cell_Line), '(A)') '/'
        Cell_Line = Cell_Line + 2
        WRITE(Oneline(Cell_line:Cell_Line+1), '(i2)') nFXR-1
        Cell_Line = Cell_Line + 2
        WRITE(Oneline(Cell_line:Cell_Line),'(A)') '*'
        Cell_Line = Cell_Line + 1
        WRite(Oneline(Cell_line:Cell_Line),'(i1)') 1
        !DO j = 1, nFXR - 1
        !  WRITE(Oneline(Cell_Line:Cell_Line), '(i1)') 1
        !  Cell_Line = Cell_Line + 2
        !END DO

        WRITE(index_file, '(A)') TRIM(Oneline)
      END DO
      CLOSE(index_file)
      MCP_CELL_MAT_INDEX = MCP_CELL_MAT_INDEX + MCP_Info(iMCP)%nMaxMicro
    END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Idx_Cell = 0
    Idx_Cell2 = 0
    Idx_Pin = 0
    MCP_ASEM_MAT_INDEX = 0
    index_MCP_Start = 1
    index_MCP_end = 0

    DO iPlane = 1, nPlane
      !Pin Input WRITE
      IF(iPlane .LT. 10) THEN
        WRITE(fn, '(A,A,A,i1,A)') 'MCP_',TRIM(CASEID),'_PIN_000',iPlane,'.include'
      ELSE IF(iPlane .GE. 10 .AND. iPlane .LT. 100) THEN
        WRITE(fn, '(A,A,A,i2,A)') 'MCP_',TRIM(CASEID),'_PIN_00',iPlane,'.include'
      ELSE IF(iPlane .GE. 100 .AND. iPlane .LT. 1000) THEN
        WRITE(fn, '(A,A,A,i3,A)') 'MCP_',TRIM(CASEID),'_PIN_0',iPlane,'.include'
      ELSE IF(iPlane .GE. 1000) THEN
        WRITE(fn, '(A,A,A,i4,A)') 'MCP_',TRIM(CASEID),'_PIN_',iPlane,'.include'
      END IF
      OPEN(UNIT = index_file, File = fn, STATUS = 'UNKNOWN', ACTION = 'WRITE', POSITION = 'REWIND')

      Macro_index = MCP_Info(index_MCP_start)%nMaxMacro

      index_MCP_end = index_MCP_end + nMCP_Plane
      DO i = 1, Macro_index
        Oneline = ''
        Idx_Cell = Idx_Cell + 1
        Idx_Pin = Idx_Pin + 1
        Idx_Cell = Idx_Pin
        Idx_Cell2 = Idx_Cell2 + 1
        Cell_Line = 3
        WRITE(Oneline(Cell_Line:Cell_Line+2), '(A)') 'pin'
        Cell_Line = 6
        WRITE(Oneline(Cell_Line:Cell_Line+4), '(i5)') Idx_Cell2
        Cell_Line = 13
        DO iPin = 1, nMCP_Plane
          Cell_Line = Cell_Line + 1
          WRITE(Oneline(Cell_Line:Cell_Line+1), '(i1)') 1
          Cell_Line = Cell_Line + 1
          WRITE(Oneline(Cell_Line:Cell_Line+1), '(A)') '*'
          Cell_Line = Cell_Line + 1

          IF(Idx_Cell .LT. 10) THEN
            WRITE(Oneline(Cell_Line:Cell_line), '(i1)') Idx_Cell
            Cell_Line = Cell_Line + 5
          ELSEIF(Idx_Cell .GE. 10 .AND. Idx_Cell .LT. 100) THEN
            WRITE(Oneline(Cell_Line:Cell_Line+1), '(i2)') Idx_Cell
            Cell_Line = Cell_Line + 5
          ELSEIF(Idx_Cell .GE. 100 .AND. Idx_Cell .LT. 1000) THEN
            WRITE(Oneline(Cell_Line:Cell_Line+2), '(i3)') Idx_Cell
            Cell_Line = Cell_Line + 5
          ELSEIF(Idx_Cell .GE. 1000 .AND. Idx_Cell .LT. 10000) THEN
            WRITE(Oneline(Cell_Line:Cell_Line+3), '(i4)') Idx_Cell
            Cell_Line = Cell_Line + 5
          END IF
          IF (iPin .NE. nMCP_Plane) Idx_Cell = Idx_Cell + Macro_index
        END DO
        WRITE(index_file, '(A)') TRIM(Oneline)
      END DO
      Idx_Pin = Idx_Pin + Macro_index * (nMCP_Plane - 1)
      index_MCP_start = index_MCP_start + nMCP_Plane
      CLOSE(index_file)
    END DO
    nCell = Idx_Cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO iMCP = 1, nPlane
      !Assembly Input WRITE
      IF(iMCP .LT. 10) THEN
        WRITE(fn, '(A,A,A,i1,A)') 'MCP_',TRIM(CASEID),'_ASSEMBLY_000',iMCP,'.include'
      ELSE IF(iMCP .GE. 10 .AND. iMCP .LT. 100) THEN
        WRITE(fn, '(A,A,A,i2,A)') 'MCP_',TRIM(CASEID),'_ASSEMBLY_00',iMCP,'.include'
      ELSE IF(iMCP .GE. 100 .AND. iMCP .LT. 1000) THEN
        WRITE(fn, '(A,A,A,i3,A)') 'MCP_',TRIM(CASEID),'_ASSEMBLY_0',iMCP,'.include'
      ELSE IF(iMCP .GE. 1000) THEN
        WRITE(fn, '(A,A,A,i4,A)') 'MCP_',TRIM(CASEID),'_ASSEMBLY_',iMCP,'.include'
      END IF
      OPEN(UNIT = index_file, File = fn, STATUS = 'UNKNOWN', ACTION = 'WRITE', POSITION = 'REWIND')
      WRITE(Oneline, '(A,i4)') '!This Assembly Include file is from MCP Number_', iMCP
      WRITE(index_file,'(A)') TRIM(Oneline)
      Macro_index = MCP_Info(iMCP)%nMaxMacro
      Oneline = ''
      IF(MCP_Info(iMCP)%Assem_Type .EQ. 8) THEN
        WRITE(Oneline(3:10),'(A)') 'assembly'
        WRITE(Oneline(12:16),'(i5)') iMCP
        WRITE(Oneline(20:21),'(i2)') 45
        WRITE(Oneline(23:24),'(i1)') 1
        WRITE(index_file, '(A)') TRIM(Oneline)

        DO i = 1, MCP_Info(iMCP)%nCell45
          Oneline = ''
          Cell_Line = 3
          Num_Start = MCP_Info(iMCP)%Micro_Number_Start(i) + MCP_ASEM_MAT_INDEX
          DO j = 1, i
            WRITE(Oneline(Cell_Line:Cell_Line+4),'(i5)') Num_Start
            Cell_Line = Cell_Line + 7
            Num_Start = Num_Start + 1
          END DO
          WRITE(index_file, '(A)') TRIM(Oneline)
        END DO
      ELSEIF(MCP_Info(iMCP)%Assem_Type .EQ. 1) THEN
        WRITE(Oneline(3:10),'(A)') 'assembly'
        WRITE(Oneline(12:16),'(i5)') iMCP
        WRITE(Oneline(19:21),'(i3)') 360
        WRITE(Oneline(23:24),'(i1)') 1
        WRITE(index_file, '(A)') TRIM(Oneline)

        DO i = 1, MCP_Info(iMCP)%nCell
          Oneline = ''
          Cell_Line = 3
          Num_Start = MCP_Info(iMCP)%Micro_Number_Start(i) + MCP_ASEM_MAT_INDEX
          DO j = 1, MCP_Info(iMCP)%nCell
            WRITE(Oneline(Cell_Line:Cell_Line+5),'(i5)') Num_Start
            Cell_Line = Cell_Line + 7
            Num_Start = Num_Start + 1
          END DO
          WRITE(index_file, '(A)') TRIM(Oneline)
        END DO
      END IF
      CLOSE(index_file)
      MCP_ASEM_MAT_INDEX = MCP_ASEM_MAT_INDEX + Macro_index
    END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (MASTER) THEN
      mesg = 'MCP File Generation Done...'
      CALL message(io8, TRUE, TRUE, mesg)
      PRINT '(A,2X,i5)', '  Total Number of Mixture      = ', MCP_CELL_MAT_INDEX
      PRINT '(A,2X,i5)', '  Total Number of Cell         = ', nCell
      PRINT '(A,2X,i5)', '  Total Number of Pin          = ', MCP_ASEM_MAT_INDEX
      PRINT '(A,2X,i5)', '  Total Number of Assembly     = ', nPlane
      PRINT '(A,2X,i5)', '  Total Number of Include File = ', nMCP * 2 + nPlane * 2
    ENDIF
    STOP
  END SUBROUTINE
END MODULE
