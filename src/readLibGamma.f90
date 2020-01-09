#include <defines.h>
#ifdef __GAMMA_TRANSPORT
SUBROUTINE readlibGamma(indev,filenm,ng,nofg,norg,ntiso)
!--- JSU EDIT 20170721

USE PARAM
USE ALLOCS
USE IOUTIL,                 ONLY : OpenFile, ToUpper
USE GammaLibdata

IMPLICIT NONE

CHARACTER :: probe, probel
CHARACTER(80) :: filenm
CHARACTER(1000) :: onelinel
CHARACTER(256) :: oneline
CHARACTER(20) :: blockname
INTEGER :: indev,ng,nofg,norg,ntiso


INTEGER :: nid,ntemp
INTEGER :: ityp,ifis,iradcap,iinel
CHARACTER*8  :: aid
REAL :: aw

INTEGER :: i,ig,ix,it,id,ic
INTEGER :: igx,itx,nx1,nx2,jg,ip,ib,iy

REAL :: xtemp(1000)

! VARIABLES TO HANDLE GAMMA TRANSPORT VARIABLES
INTEGER :: ngg, igg, imt

EQUIVALENCE(probe, oneline)
EQUIVALENCE(probel, onelinel)

!!    1. open file
CLOSE(indev)
CALL openfile(indev,TRUE,FALSE, FALSE, filenm)

!!    2. Read basic information
Do WHILE(.TRUE.)
  READ(indev,'(a1000)', END=201) onelinel
  onelinel=trim(adjustl(onelinel))
  IF(probel.eq.'/') EXIT
  IF(probel.eq.' ' .or. probel.eq.'' .or. onelinel.eq.' '.or. onelinel.eq.'\t') CYCLE
  READ(onelinel,*) blockname
  blockname=trim(adjustl(blockname))
  CALL toupper(blockname)

  SELECT CASE(blockname)

! Dimension
  CASE('!DIMENSION')
    READ(indev,*) nogGAM,noggGAM,nofgGAM,norgGAM,notgGAM,neltGAM
    ng=nogGAM
    ntiso=neltGAM
    ngg =noggGAM

! Enery Group Boundary       [upper bound, e(ng+1)=1.e-4 eV]
  CASE('!NEUTRON')
    CALL dmalloc(NeuENB,nogGAM)
    CALL dmalloc(NeutronU,nogGAM)
    ! read(indev, *) ! "Group Boundary"
    READ(indev, '(10E14.7)') (NeutronU(ig), ig=1, nogGAM)
    DO i=1,nogGAM
      NeuENB(i)=1.0e7_8/exp(NeutronU(i))
    END DO
! Photon Energy Group Boundary,   [ upper bound, e(ng+1) = 1.e+3 eV]
  CASE('!GAMMA')  
    CALL dmalloc(GamENB, NOGGGAM)
    CALL dmalloc(GammaU, NOGGGAM)
    READ(indev, '(10E14.7)') (GammaU(igg), igg = 1, noggGAM)
    Do igg = 1, noggGAM
        GamENB(igg) = 1.0E7_8/EXP(GammaU(igg))
    End Do
! Nuclide information
  CASE('!DIR')
    ALLOCATE(Gldiso(neltGAM))
    ! read(indev,*)  !"%DIR:"
    DO i=1,neltGAM
      READ(indev,*) ix,nid,aw,ityp,ntemp,aid
      
      IF(ix.ne.i) THEN
        mesg='readlibGamma.F - Check index of !DIR'
        WRITE(*,'(a)') mesg
        STOP
      END IF
      Gldiso(i)%nid=nid
      Gldiso(i)%aw=aw
      Gldiso(i)%ityp=ityp
      !Gldiso(i)%ifis=ifis
      !Gldiso(i)%iradcap=iradcap
      !Gldiso(i)%iinel=iinel
      Gldiso(i)%ntemp=ntemp
      Gldiso(i)%aid=aid
    END DO
  CASE('!NUC:')
    BACKSPACE(indev)
    EXIT
  CASE DEFAULT
    mesg='readlibGamma - BLOCK NAME '//blockname//' Not Allowed...'
    WRITE(*,'(a)') mesg
    STOP 'readlibGamma.f90'
  END SELECT
END DO

!!    3. Check nog and nogg (neutron and gamma energy group structure) 
IF(nogGAM.gt.1000) THEN
  mesg='readlibPLC.F - Check xtemp size '
  WRITE(*,'(a)') mesg
  STOP
END IF


!!    4. Memory Allocate
DO i=1, neltGAM
  nid = Gldiso(i)%nid
  ityp = Gldiso(i)%ityp
  IF (ityp .GT. 0) THEN !
    ntemp = Gldiso(i)%ntemp
  ! PHOTON PRODUCTION DATA (INDUCED BY NEUTRON)
  ! Temperature
    CALL dmalloc(Gldiso(i)%temp,ntemp)
  ! Production Matrix
    ALLOCATE(Gldiso(i)%GPM(NGG,NTEMP,0:3))
    ALLOCATE(Gldiso(i)%GPMOUTB(NG,NTEMP,0:3))
    ALLOCATE(Gldiso(i)%GPMOUTE(NG,NTEMP,0:3))
  END IF
  IF (ityp .NE. 1) THEN
  ! PHOTON REACTION DATA
  ! XSP+ (photon cross section), siga, sigtr, sigtrs, sigs0
    CALL dmalloc(Gldiso(i)%GSIGA,NGG)
    CALL dmalloc(Gldiso(i)%GSIGTR,NGG)
    CALL dmalloc(Gldiso(i)%GSIGS,NGG)
    CALL dmalloc(Gldiso(i)%GSIGSTR,NGG)
    CALL dmalloc(Gldiso(i)%KERMA,NGG)
  ! P0~P3 Scattering for photon
    CALL dmalloc(Gldiso(i)%GSIGSP1, NGG)
    CALL dmalloc(Gldiso(i)%GSIGSP2, NGG)
    CALL dmalloc(Gldiso(i)%GSIGSP3, NGG)
    ALLOCATE(Gldiso(i)%GSM(ngg))
    ALLOCATE(Gldiso(i)%GSM1(ngg))
    ALLOCATE(Gldiso(i)%GSM2(ngg))
    ALLOCATE(Gldiso(i)%GSM3(ngg))
  END IF
END DO

!!    5. read nuclide information
DO WHILE(TRUE)
  READ(indev,'(a1000)', END=202) onelinel
  onelinel=trim(adjustl(onelinel))
  IF(probel.eq.'/') EXIT
  IF(probel.eq.' ' .or. probel.eq.'' .or. onelinel.eq.' '.or. onelinel.eq.'\t') CYCLE
  READ(onelinel,*) blockname
  CALL toupper(blockname)

  ! check line
  READ(indev,'(a256)') oneline
  oneline=trim(adjustl(oneline))
  BACKSPACE(indev);
  IF((oneline.eq.' ' .or. probe.eq.'!') .and. blockname.ne.'!NUC:') CYCLE

  SELECT CASE(blockname)
    CASE('!NUC:')
      READ(onelinel, *) mesg,                                                                       &
            id,nid,aw,ityp,ntemp,aid      
    CASE('!XSP+')
      IF (Gldiso(id)%ityp .EQ. 1) STOP 'SUBROUTINE ReadLibGamma - ityp = 1 and reaction exists'
      DO igg=1,ngg
        READ(indev,*) igx,Gldiso(id)%KERMA(igg),Gldiso(id)%GSIGA(igg),                             &
                Gldiso(id)%GSIGTR(igg),Gldiso(id)%GSIGS(igg),                                      &
                nx1,nx2,(xtemp(jg),jg=nx1,nx2)
!        Gldiso(id)%GSIGSTR = Gldiso(id)%GSIGTR(igg) - Gldiso(id)%GSIGA(igg)
        Gldiso(id)%GSM(igg)%ib=nx1
        Gldiso(id)%GSM(igg)%ie=nx2
        ALLOCATE(Gldiso(id)%GSM(igg)%from(nx1:nx2))
        DO jg=nx1,nx2
           Gldiso(id)%GSM(igg)%from(jg)=xtemp(jg)
        END DO
      END DO
    CASE('!PA1+')
      DO igg=1,ngg
         READ(indev,*) igx,Gldiso(id)%GSIGSP1(igg),nx1,nx2,(xtemp(jg),jg=nx1,nx2)
         Gldiso(id)%GSM1(igg)%ib=nx1
         Gldiso(id)%GSM1(igg)%ie=nx2
         ALLOCATE(Gldiso(id)%GSM1(igg)%from(nx1:nx2))
         DO jg=nx1,nx2
             Gldiso(id)%GSM1(igg)%from(jg)=xtemp(jg)
         END DO
      END DO
    CASE('!PA2+');
      DO igg=1,ngg
         READ(indev,*) igx,Gldiso(id)%GSIGSP2(igg),nx1,nx2,(xtemp(jg),jg=nx1,nx2)
         Gldiso(id)%GSM2(igg)%ib=nx1
         Gldiso(id)%GSM2(igg)%ie=nx2
         ALLOCATE(Gldiso(id)%GSM2(igg)%from(nx1:nx2))
         DO jg=nx1,nx2
             Gldiso(id)%GSM2(igg)%from(jg)=xtemp(jg)
         END DO
      END DO
    CASE('!PA3+');
      DO igg=1,ngg
         READ(indev,*) igx,Gldiso(id)%GSIGSP3(igg),nx1,nx2,(xtemp(jg),jg=nx1,nx2)
         Gldiso(id)%GSM3(igg)%ib=nx1
         Gldiso(id)%GSM3(igg)%ie=nx2
         ALLOCATE(Gldiso(id)%GSM3(igg)%from(nx1:nx2))
         DO jg=nx1,nx2
           Gldiso(id)%GSM3(igg)%from(jg)=xtemp(jg)
         END DO
      END DO
    CASE('!TP1+')
      IF (Gldiso(id)%ityp .EQ. 0) STOP 'SUBROUTINE ReadLibGamma - ityp = 0 and production exists'
      READ(indev, *) (Gldiso(id)%temp(it), it=1, ntemp)
    CASE('!GAM+')
      IF (Gldiso(id)%ityp .EQ. 0) STOP 'SUBROUTINE ReadLibGamma - ityp = 0 and production exists'
      READ(indev,*) (Gldiso(id)%lphoton(imt),imt=1,3)
      READ(indev,*) Gldiso(id)%ppgrLOW, Gldiso(id)%ppgrUP
      Gldiso(id)%lphoton(imtTOT) = any(Gldiso(id)%lphoton(1:3))
      IF (Gldiso(id)%lphoton(imtFIS)) Gldiso(id)%ifis = 1
      IF (Gldiso(id)%lphoton(imtRAD)) Gldiso(id)%iradcap = 1
      IF (Gldiso(id)%lphoton(imtINEL)) Gldiso(id)%iinel = 1
      DO imt=0,3
        IF (.not. Gldiso(id)%lphoton(imt)) CYCLE
        READ(indev,*)
        DO igg=Gldiso(id)%ppgrlow(imt),Gldiso(id)%ppgrup(imt)
          DO it=1,gldiso(id)%ntemp
            READ(indev,*) igx,itx,nx1,nx2,(xtemp(jg),jg=nx1,nx2)
            Gldiso(id)%GPM(igg,it,imt)%ib=nx1
            Gldiso(id)%GPM(igg,it,imt)%ie=nx2
            ALLOCATE(Gldiso(id)%GPM(igg,it,imt)%from(nx1:nx2))
            DO jg=nx1,nx2
              Gldiso(id)%GPM(igg,it,imt)%from(jg)=xtemp(jg)
            END DO
          END DO
        END DO
      END DO
    ! EXPLICIT KAPPA  -- JSU EDIT 2017/09/12
    CASE('!KAW+') ! FISSION FRAGMENTWISE KAPPA
      READ(indev, *) (Gldiso(id)%exkappa(ic), ic = 1, 6) 
      Gldiso(id)%lfis = .TRUE.
    CASE('!KA0+') ! INCIDENT ENERGY INDEPENDENT KAPPA
      READ(indev, *) (Gldiso(id)%exkappa0(ic), ic = 1, 6)
    CASE('!KCA+') ! RADIATIVE CAPTURE KAPPA
      READ(indev, *) Gldiso(id)%capkappa
    CASE('!QCA+')
      READ(indev, *) Gldiso(id)%CapQ, Gldiso(id)%n2nQ, Gldiso(id)%n3nQ
    CASE DEFAULT
      mesg='BLOCK NAME '//blockname//' Not Allowed...'
      WRITE(*,'(a)') mesg
      STOP 'readlibGamma'
  END SELECT
END DO
202 CONTINUE

!!    6. close file
CLOSE(indev)

!!    7. mapping
CALL mappingGAM

!!    8. ordering for scattering      
CALL orderGAMlib

!!    9. Group Information for Gamma
CALL SetGamGroupInfo  
RETURN

201   CONTINUE
STOP 'readlibGamma - stop 1'

RETURN
END SUBROUTINE


SUBROUTINE mappingGAM
!  -------------------------------------------------------------------------------------------------
!  Mapping for DeCART VHTR to Helios Libaray related Variable
!  -------------------------------------------------------------------------------------------------
USE PARAM
USE GammaLibdata
USE ALLOCS

IMPLICIT NONE

INTEGER :: i,j,k
INTEGER :: ifis, irad, iinel
INTEGER :: idelem, ielm
INTEGER :: it1, it2, nid, m
INTEGER :: itemp2, itemp1

nfisGAM  = 0
nradGAM  = 0
ninelGAM = 0

! Nuclide List
CALL dmalloc(nuclidGAM,neltGAM)
CALL dmalloc(elementGAM,neltGAM)
CALL dmalloc(iso2elm,neltGAM)
DO i=1,neltGAM
  nuclidGAM(i)=Gldiso(i)%nid
END DO

!  for ityp > 0 (photon source isotopes)
DO i=1,neltGAM
  IF(Gldiso(i)%lphoton(imtFIS)) nfisGAM = nfisGAM + 1 
  IF(Gldiso(i)%lphoton(imtRAD)) nradGAM = nradGAM + 1
  IF(Gldiso(i)%lphoton(imtINEL)) ninelGAM = ninelGAM + 1
END DO

CALL dmalloc(idfisGAM,nfisGAM)   !  list of fission isotope (imt=1)
CALL dmalloc(idradGAM,nradGAM)   !  list of radioactive capture isotope (imt=2)
CALL dmalloc(idinelGAM,ninelGAM) !  list of inelastic scattering isotope (imt=3)
ifis = 0
irad = 0
iinel = 0
DO i = 1, neltGAM
  IF(Gldiso(i)%lphoton(imtFIS)) THEN
    ifis = ifis + 1
    idfisGAM(ifis)=i
  END IF
  IF(Gldiso(i)%lphoton(imtRAD)) THEN
    irad = irad + 1
    idradGAM(irad)=i
  END IF
  IF(Gldiso(i)%lphoton(imtINEL)) THEN
    iinel = iinel + 1
    idinelGAM(iinel)=i
  END IF
END DO

! for ityp = 0 (Photon Reaction Data)
!     matching isotopes to corresponding element
nelmGAM = 0
DO i = 1,neltGAM
!  idelem = (nuclidGAM(i) / 1000) * 1000
!  if (nuclidGAM(i) .eq. idelem) nelmGAM = nelmGAM + 1
  IF (Gldiso(i)%ityp .EQ. 1) CYCLE ! (ityp=1) is pure production
  nelmGAM = nelmGAM + 1
END DO
ALLOCATE(idelmGAM(nelmGAM)) ! library data index
ielm = 0
Do i = 1,neltGAM
  idelem = (nuclidGAM(i) / 1000) * 1000
  IF (nuclidGAM(i) .NE. idelem) CYCLE
  ielm = ielm + 1
  idelmGAM(ielm) = i
End Do
! isotope -> element index matching
Do i = 1, neltGAM
  idelem = (nuclidGAM(i) / 1000) * 1000
  Do j = 1, nelmGAM
    IF (nuclidGAM(idelmGAM(j)).eq.idelem) EXIT
  End Do
  elementGAM(i) = idelem
  iso2elm(i) = idelmGAM(j)
End Do

mapnuclGAM = 0
DO k=1,neltGAM
  mapnuclGAM(nuclidGAM(k)) = k
END DO

mapnuclELM = 0
DO k = 1, neltGAM
  mapnuclELM(nuclidGAM(k)) = iso2elm(k)
END DO

CALL dmalloc(itempmapGAM,nxtempmap,neltGAM)
itempmapGAM = 0
tempmax=0
DO nid=1,neltGAM
  IF (Gldiso(nid)%ntemp .lt. 1) CYCLE
  it1 = Gldiso(nid)%ntemp
  IF(Gldiso(nid)%temp(it1).gt.tempmax) tempmax = Gldiso(nid)%temp(it1)
END DO

DO nid = 1, neltGAM
  IF (Gldiso(nid)%ityp .lt. 1) CYCLE
  it1=Gldiso(nid)%ntemp
  it2=0
  IF (it1.gt.1) THEN
    itemp2 = NINT(Gldiso(nid)%temp(1))
    DO i=1,itemp2
      itempmapGAM(i,nid) = 1
    END DO
    DO m=2,it1
      itemp1=itemp2+1
      itemp2=NINT(Gldiso(nid)%temp(m))
      DO i=itemp1,itemp2
        itempmapGAM(i,nid)=m-1
      END DO
    END DO
    itemp1=itemp2+1
    DO i=itemp1,nxtempmap
      itempmapGAM(i,nid)=it1-1
    END DO
  ELSE
    DO i=1,nxtempmap
      itempmapGAM(i,nid)=1
    END DO
  END IF
END DO

RETURN
END SUBROUTINE
      
      
SUBROUTINE orderGAMlib

USE PARAM
USE GammaLibdata,    ONLY :  Gldiso, nelmGAM, noggGAM, idelmGAM, nogGAM, neltGAM

IMPLICIT NONE

INTEGER :: i,igg,it,igx,ntemp, ig, imt
INTEGER :: ielm

!--Out Scattering Range
DO ielm = 1, nelmGAM
  i = idelmGAM(ielm)   ! library data index
  DO igg = 1, noggGAM
    DO igx = 1, noggGAM         ! igg -> igx, beginning
      IF(igg.ge.Gldiso(i)%GSM(igx)%ib .and. igg.le.Gldiso(i)%GSM(igx)%ie) THEN
        Gldiso(i)%GSM(igg)%ioutsb=igx  ! upper(energy) group boundary from igg
        EXIT
      END IF
    END DO
    DO igx=noggGAM,1,-1         ! igg -> igx, ending
      IF(igg.ge.Gldiso(i)%GSM(igx)%ib .and. igg.le.Gldiso(i)%GSM(igx)%ie) THEN
        Gldiso(i)%GSM(igg)%ioutse=igx  ! upper(energy) group boundary from igg
        EXIT
      END IF
    END DO
  END DO ! group loop
END DO ! element loop

!- PRODUCTION RANGE
!      RANGE OF PHORON GROUP CORRESPONDING NEUTRON GROUP
DO i = 1,neltGAM
  IF (Gldiso(i)%ityp .EQ. 0) CYCLE
  Gldiso(i)%GPMOUTB = 47
  Gldiso(i)%GPMOUTE = 0
  DO imt = 1, 3
    IF (.NOT. Gldiso(i)%LPHOTON(imt)) CYCLE
    DO it = 1, Gldiso(i)%ntemp
      DO ig = 1, nogGAM
        DO igg = 1, noggGAM         ! ig -> igg, beginning
          IF(ig .GE. Gldiso(i)%GPM(igg,it,imt)%ib .AND. ig .LE. Gldiso(i)%GPM(igg,it,imt)%ie) THEN
            Gldiso(i)%GPMOUTB(ig,it,imt) = igg  ! upper(energy) group boundary from igg
            EXIT
          END IF
        END DO
        DO igg=noggGAM,1,-1         ! ig -> igg, ending
          IF(ig .GE. Gldiso(i)%GPM(igg,it,imt)%ib .AND. ig .LE. Gldiso(i)%GPM(igg,it,imt)%ie) THEN
            Gldiso(i)%GPMOUTE(ig,it,imt) = igg  ! upper(energy) group boundary from igg
            EXIT
          END IF
        END DO
        Gldiso(i)%GPMOUTB(ig,it,0) = MIN(Gldiso(i)%GPMOUTB(ig,it,0),Gldiso(i)%GPMOUTB(ig,it,imt))
        Gldiso(i)%GPMOUTE(ig,it,0) = MAX(Gldiso(i)%GPMOUTE(ig,it,0),Gldiso(i)%GPMOUTE(ig,it,imt))
      END DO ! NEUTRON GROUP LOOP (ig)
    END DO ! TEMPERATURE LOOP (it)
  END DO ! REACTION TYPE LOOP (imt)
END DO

DO ielm = 1,nelmGAM
  i = idelmGAM(ielm)   ! library data index
  Gldiso(i)%GSIGSTR=zero ! igg -> igx, FROM GROUP XS SUMMATION
  ! Diagonal term of scattering matrix is already transport corrected (outflow for gamma)
  DO igg=1,noggGAM
    DO igx= Gldiso(i)%GSM(igg)%ib, Gldiso(i)%GSM(igg)%ie
      Gldiso(i)%GSIGSTR(igx)=Gldiso(i)%GSIGSTR(igx) + Gldiso(i)%GSM(igg)%from(igx)
    END DO
  END DO
END DO! element loop

RETURN
END SUBROUTINE
      
!Set Group Sweeping and Scattering Range Setting
SUBROUTINE SetGamGroupInfo
USE PARAM
USE ALLOCS
USE GammaLibdata,      ONLY : neltGAM,        noggGAM,        nelmGAM,       nogGAM,          &
                               Gldiso,        idelmGAM,       NeuENB,        GamENB
USE GammaCore_mod,     ONLY : GamGroupInfo
IMPLICIT NONE

INTEGER :: ng, ngg
INTEGER :: i, j, k, igg, igg2, iso, ig
INTEGER :: nid, ntemp

ng = nogGAM
ngg = noggGAM

GamGroupInfo%ng = ng
GamGroupInfo%ngg = ngg

ALLOCATE(GamGroupInfo%InScatRange(2, ngg))
ALLOCATE(GamGroupInfo%OutScatRange(2, ngg))
GamGroupInfo%InScatRange(1, :) = ngg;  GamGroupInfo%InScatRange(2, :) = 1;
GamGroupInfo%OutScatRange(1, :) = ngg; GamGroupInfo%OutScatRange(2, :) = 1;

!InScattering Info
GamGroupInfo%UpScatRange = ngg
DO iso = 1, nelmGAM
  nid = idelmGAM(iso)
  DO igg = 1, ngg
    !Inscattering Range
    GamGroupInfo%InScatRange(1, igg) = min(GamGroupInfo%InScatRange(1, igg), Gldiso(nid)%GSM(igg)%ib)
    GamGroupInfo%InScatRange(2, igg) = max(GamGroupInfo%InScatRange(2, igg), Gldiso(nid)%GSM(igg)%ie)          
    !Out Scattering Range
    DO igg2 = Gldiso(nid)%GSM(igg)%ib, Gldiso(nid)%GSM(igg)%ie
      GamGroupInfo%OutScatRange(1, igg2) = min(Gldiso(nid)%GSM(igg)%ib, GamGroupInfo%OutScatRange(1, igg2))
      GamGroupInfo%OutScatRange(2, igg2) = max(Gldiso(nid)%GSM(igg)%ie, GamGroupInfo%OutScatRange(2, igg2))
    ENDDO
    !Upscatering Range
    IF(Gldiso(nid)%GSM(igg)%ie .GT. igg) THEN
      GamGroupInfo%UpScatRange(1) = min(igg, GamGroupInfo%UpScatRange(1))
      CONTINUE
    ENDIF
  ENDDO  !Group Sweep
ENDDO  !Isotope Sweep

GamGroupInfo%UpScatRange(2) = ngg
GamGroupInfo%lUpScat = FALSE
DO igg = 1, ngg
  IF(GamGroupInfo%InScatRange(2, igg) .GT. igg) THEN
    GamGroupInfo%lUpScat = TRUE
    GamGroupInfo%UpScatRange(1) = igg
    EXIT
  ENDIF
ENDDO

! GROUP REPRESENTATIVE ENERGY (GEOMETRIC AVERAGE)                         |-- JSU EDIT 2017.09.14. |
ALLOCATE(GamGroupInfo%GamAvgE(ngg), GamGroupInfo%NeuAvgE(ng))
GamGroupInfo%GamAvgE = 0.
DO igg = 1, ngg
  IF (igg .EQ. ngg) THEN
    GamGroupInfo%GamAvgE(igg) = SQRT(GamENB(ngg) * 1e+3)
  ELSE
    GamGroupInfo%GamAvgE(igg) = SQRT(GamENB(igg) * GamENB(igg + 1))
  END IF
END DO
GamGroupInfo%NeuAvgE = 0.
DO ig = 1, ng
  IF (ig .EQ. ng) THEN
    GamGroupInfo%NeuAvgE(ig) = SQRT(NeuENB(ig) * 1e-4)
  ELSE
    GamGroupInfo%NeuAvgE(ig) = SQRT(NeuENB(ig) * NeuENB(ig + 1))
  END IF
END DO

GamGroupInfo%GamAvgE = GamGroupInfo%GamAvgE * eVtoJ  ! eV to Joule
GamGroupInfo%NeuAvgE = GamGroupInfo%NeuAvgE * eVtoJ  ! eV to Joule
END SUBROUTINE
#endif