#include <defines.h>
SUBROUTINE MakeRays()
USE PARAM
USE GEOM,   ONLY : AsyPitch,           AsyGeom,          lEdge,          &
                   CellInfo,           Core,             Asy,            &
                   nCellType,          nCellType0,       nBaseCell,      &
                   BaseCellInfo
USE RAYS,   ONLY : RayInfo,            AziAngle,         PolarAngle,     &
                   ModRay,             AsyRay,           CoreRay,        &
                   RotRay,             CellRayBase,      RayInfo4Cmfd,   &
                   CellRay1D,          FastCoreRayDat,                   &
                   !--- CNJ Edit : Domain Decomposition
                   DcmpAsyRay
USE SetRay, ONLY : SetModularRayAngle, SetPolarRayAngle, SetModRay,      &
                   SetAsyRay,          CellRayGen,       CoreRayGen,     &
                   RotationRayGen,     RayInfo4CmfdGen,                  &
                   !--- CNJ Edit : Domain Decomposition
                   DcmpRayGen
USE PE_MOD, ONLY : PE
#ifdef MPI_ENV
USE MPICOMM_MOD,    ONLY : MPI_SYNC
USE MPIConfig_Mod,  ONLY :SetRayPEVariables
#endif
USE IOUTIL, ONLY : message, terminate    
USE FILES,  ONLY : io8
USE ALLOCS
USE CNTL
IMPLICIT NONE
INTEGER :: iz, icel, jcel
LOGICAL :: lsame
LOGICAL :: master

master = PE%master
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
mesg = 'Setting up Ray Tracing Data...'
IF(master) CALL message(io8,TRUE, TRUE, mesg)

IF (nTracerCntl%lRayGen) THEN

#ifdef DetailedIO
mesg = '  Generating Modular and Assembly Rays...'
IF(master) CALL message(io8,TRUE, TRUE, mesg)
#endif
CALL SetModularRayAngle(AziAngle, RayInfo, AsyPitch)
CALL SetPolarRayAngle(PolarAngle, RayInfo)
CALL SetModRay(ModRay, RayInfo, AsyGeom, ASyPitch, lEdge)
CALL SetAsyRay(AsyRay, RayInfo, AsyGeom, lEdge)

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

IF (.NOT. nTracerCntl%lnTIGRst) THEN

  mesg = '  Finding Base Cell Structures...'
  IF(master) CALL message(io8,TRUE, TRUE, mesg)
  CellInfo(1)%basecellstr=1
  DO icel = 2, nCellType0
    CellInfo(icel)%basecellstr=icel
    IF (nTracerCntl%lBCRopt .AND. .NOT. CellInfo(icel)%lGap) THEN
      !--- CNJ Edit : User Defined Gap Test
      !--start searching
      DO jcel = 1, icel-1
        IF (CellInfo(jcel)%basecellstr .EQ. jcel) THEN !is basecell
          CALL ChkSameGeomStr(CellInfo(icel), CellInfo(jcel), lsame)
          IF (lSame) THEN
            CellInfo(iCel)%basecellstr=jcel
            CONTINUE
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (nCellType0 .NE. nCellType) THEN
    DO icel = nCellType0 + 1, nCellType
      CellInfo(icel)%basecellstr = CellInfo(icel - nCellType0)%basecellstr + nCellType0
    ENDDO
  ENDIF

#ifdef DetailedIO
  mesg = '  Generating Base Cell Rays...'
  IF(master) CALL message(io8,TRUE, TRUE, mesg)
#endif

  ALLOCATE(CellRayBase(nCellType))

  CALL omp_set_num_threads(PE%nThread)

  !$OMP PARALLEL
  !$OMP DO SCHEDULE(DYNAMIC)
  DO icel = 1, nCellType
    IF (CellInfo(iCel)%BaseCellStr .EQ. icel) THEN
      IF (CellInfo(icel)%luse .OR. nTracerCntl%lBCRopt .AND. .NOT. (CellInfo(icel)%lCCell .AND. (icel .GT. nCellType0))) THEN
        CALL SetCellRayBase(CellInfo(icel), RayInfo, CellRayBase(icel))
      ENDIF
    ENDIF
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  
ELSE

#ifdef DetailedIO
  mesg = '  Generating Base Cell Rays...'
  IF(master) CALL message(io8,TRUE, TRUE, mesg)
#endif

  ALLOCATE(CellRayBase(nBaseCell))

  CALL omp_set_num_threads(PE%nThread)

  !$OMP PARALLEL
  !$OMP DO SCHEDULE(DYNAMIC)
  DO icel = 1, nBaseCell
    IF (BaseCellInfo(icel)%luse .OR. nTracerCntl%lBCRopt .AND. .NOT. (BaseCellInfo(icel)%lCCell)) THEN
      CALL SetCellRayBase(BaseCellInfo(icel), RayInfo, CellRayBase(icel))
      CellInfo(icel)%CellRay => BaseCellInfo(icel)%CellRay
    ENDIF
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  
ENDIF

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

#ifdef DetailedIO
mesg = '  Generating Actual Cell Rays...'
IF(master) CALL message(io8,TRUE, TRUE, mesg)
#endif
CALL CellRayGen(RayInfo, CellRayBase, AsyGeom, lEdge)

IF(nTracerCntl%lFastMOC) CALL MakeCellRay1DArray(RayInfo, CellRayBase, CellInfo, nCellType)
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

#ifdef DetailedIO
mesg = '  Generating Core Rays...'
IF(master) CALL message(io8,TRUE, TRUE, mesg)
#endif
CALL CoreRayGen(RayInfo, CoreRay, Core, Asy, AsyGeom, lEdge)

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

#ifdef DetailedIO
mesg = '  Generating Cycle Rays...'
IF(master) CALL message(io8,TRUE, TRUE, mesg)
#endif
CALL RotationRayGen(RayInfo, RotRay, Core, Asy, lEdge)

IF (nTracerCntl%RayGenOpt .GE. 1) THEN !Ray Structure printout
  !Print RayFile
  CALL PrintRayFile()
  IF(nTracerCntl%RayGenOpt .EQ. 2) THEN !Ray Structure printout +Terminate
    !Terminate
    CALL Terminate('Terminated after generating the Ray Structure File')
  ENDIF
ENDIF

ELSE !lRayGen = .FALSE.
  CALL ReadRayFile()
  !READ RayFile
ENDIF

CALL RayInfoMaxSize(Core, RayInfo, PE%myzb, PE%myze)
IF(nTracerCntl%lFastMOC) CALL SetFastMOC(RayInfo, Core, nTracerCntl, PE)
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

#ifdef MPI_ENV
CALL SetRayPEVariables(PE, Core, RayInfo)
#endif

!--- CNJ Edit : Domain Decomposition
IF (nTracerCntl%lDomainDcmp) THEN
#ifdef DetailedIO
  mesg = 'Generating Decomposed Assembly Rays...'
  IF(master) CALL message(io8,TRUE, TRUE, mesg)
#endif
  CALL DcmpRayGen(Core, RayInfo, DcmpAsyRay)
ENDIF

!--- CNJ Edit : CASMO Linear Source
IF (nTracerCntl%lLinSrcCASMO) THEN
#ifdef DetailedIO
  mesg = '  Calculating FSR Centroids and Moments...'
  IF(master) CALL message(io8,TRUE, TRUE, mesg)
#endif
  DO iz = PE%myzb, PE%myze
    CALL GetFSRCentroid(RayInfo, Core, iz, nTracerCntl%lHybrid)
  ENDDO
#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
ENDIF

CALL RayInfo4CmfdGen(RayInfo, RayInfo4Cmfd, Core)

#ifdef DetailedIO
mesg = '  Finished Generating Ray Data...'
IF(master) CALL message(io8,TRUE, TRUE, mesg)
#endif

END SUBROUTINE