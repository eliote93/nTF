#include <defines.h>

SUBROUTINE ProcessPhims(Core, FmInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF
USE CNTL,             ONLY : nTracerCntl_Type,    OutpCntl_Type
USE FILES,            ONLY : io16,                caseid
USE MOC_Mod,          ONLY : FxrAvgPhi
USE XSLIB_MOD
USE ioutil
USE BasicOperation
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE


INTEGER :: ixa, iya, ixc, iyc, nxc, nyc, ixy, ixya, ixyc
INTEGER :: iz, iasytype, icel
INTEGER :: i, j, k
INTEGER :: io, ilpho, impho, nmpho, nphiout
INTEGER :: nfxr, ifxr, ir, ifsr, ifsrinfxr, nFsrInFxr, ilocalfsr
INTEGER :: ig, ng

!Output Data
REAL, POINTER :: Area(:)
REAL :: vol
CHARACTER(256) :: fn

#ifdef __GFORTRAN__
CHARACTER(10) :: str2num
#endif

TYPE phim_type
    INTEGER :: nazimom
    REAL,POINTER :: aziphiout(:,:)               ! ng,# of azimuthal moments
END TYPE

TYPE regphimout_type                             ! nTracerCntl%OutpCntl%PhimOrdOutList(0, i)
    INTEGER :: nphimout
    TYPE(phim_type),POINTER :: phimout(:)
ENDTYPE

TYPE(regphimout_type), POINTER :: regphimout(:)  ! nTracerCntl%OutpCntl%nRegPhimOut

INTEGER :: ordermap(4,3),nordermap(0:3)
DATA nordermap / 1, 2, 3, 4 /
DATA ordermap  / 1, 2, 0, 0, &
                 3, 4, 5, 0, &
                 6, 7, 8, 9/

ng = GroupInfo%ng
ALLOCATE(Area(nTracerCntl%OutpCntl%nRegPhimOut))
ALLOCATE(regphimout(nTracerCntl%OutpCntl%nRegPhimOut))
DO i = 1, nTracerCntl%OutpCntl%nRegPhimOut
  iz = nTracerCntl%OutpCntl%RegPhimOutList(1, i)
  ixa = nTracerCntl%OutpCntl%RegPhimOutList(2, i); iya = nTracerCntl%OutpCntl%RegPhimOutList(3, i)
  ixya = Core%CoreIdx(ixa, iya); iasytype = Core%CoreMap(ixya)
  nphiout = nTracerCntl%OutpCntl%PhimOrdOutList(0, i)
  regphimout(i)%nphimout=nphiout
  ALLOCATE(regphimout(i)%phimout(nphiout))
  IF (nTracerCntl%OutpCntl%RegPhimOutASM(i)) THEN
      DO io=1,nphiout
          ilpho=nTracerCntl%OutpCntl%PhimOrdOutList(io, i)
          nmpho=nordermap(ilpho); regphimout(i)%phimout(io)%nazimom=nmpho
          ALLOCATE(regphimout(i)%phimout(io)%aziphiout(ng,nmpho))
          regphimout(i)%phimout(io)%aziphiout=0
          nxc = Core%AsyInfo(iasytype)%nx; nyc = Core%AsyInfo(iasytype)%ny;
          AREA(i) = 0
          DO iyc=1,nyc
              DO ixc=1,nxc
                  ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
                  ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Core%Pin(ixy)%Cell(iz)
                  nfxr = Core%CellInfo(icel)%nFXR
                  DO ifxr= 1, nfxr
                      nFsrInFxr = Core%CellInfo(icel)%nFsrInFxr(ifxr)
                      DO ifsrinfxr = 1, nFsrInFxr
                          ilocalfsr=Core%CellInfo(icel)%MapFxr2FsrIdx(ifsrinfxr,ifxr)
                          vol = Core%CellInfo(icel)%vol(ilocalfsr); AREA(i) = Area(i) + vol;
                          ifsr=Core%Pin(ixy)%FsrIdxst +ilocalfsr-1
                          !write(101,'(5i4,ES16.6)') ixy,icel,ifxr,ilocalfsr,ifsr,vol
                          DO ig=1,ng
                              DO impho=1,nmpho
                                  if (ilpho.eq.0) then
                                      regphimout(i)%phimout(io)%aziphiout(ig,impho) = regphimout(i)%phimout(io)%aziphiout(ig,impho) + vol * FMInfo%Phis(ifsr,iz,ig)
                                  else
                                      regphimout(i)%phimout(io)%aziphiout(ig,impho) = regphimout(i)%phimout(io)%aziphiout(ig,impho) + vol * FMInfo%Phim(ordermap(impho,ilpho),ifsr,iz,ig)
                                  endif
                              ENDDO ! DO impho=1,nmpho
                          ENDDO ! DO ig=1,ng
                      ENDDO ! DO ifsrinfxr = 1, nFsrInFxr
                  ENDDO ! DO ifxr= 1, nfxr
              ENDDO ! DO ixc=1,nxc
          ENDDO ! DO iyc=1,nyc
      ENDDO ! DO io=1,nphiout
      DO io=1,nphiout
          DO impho=1,regphimout(i)%phimout(io)%nazimom
              DO ig=1,ng
                 regphimout(i)%phimout(io)%aziphiout(ig,impho)=regphimout(i)%phimout(io)%aziphiout(ig,impho)/AREA(i)
              ENDDO
          ENDDO
      ENDDO
  ELSEIF (nTracerCntl%OutpCntl%RegPhimOutPIN(i)) THEN
      DO io=1,nphiout
          ilpho=nTracerCntl%OutpCntl%PhimOrdOutList(io, i)
          nmpho=nordermap(ilpho); regphimout(i)%phimout(io)%nazimom=nmpho
          ALLOCATE(regphimout(i)%phimout(io)%aziphiout(ng,nmpho))
          regphimout(i)%phimout(io)%aziphiout=0
          ixc = nTracerCntl%OutpCntl%RegPhimOutList(4, i); iyc = nTracerCntl%OutpCntl%RegPhimOutList(5, i)
          AREA(i) = 0
          if (ixc.ne.0.and.iyc.ne.0) then
              ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc); ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc)
              icel = Core%Pin(ixy)%Cell(iz); nfxr = Core%CellInfo(icel)%nFXR
              DO ifxr= 1, nfxr
                  nFsrInFxr = Core%CellInfo(icel)%nFsrInFxr(ifxr)
                  DO ifsrinfxr = 1, nFsrInFxr
                      ilocalfsr=Core%CellInfo(icel)%MapFxr2FsrIdx(ifsrinfxr,ifxr)
                      vol = Core%CellInfo(icel)%vol(ilocalfsr); AREA(i) = Area(i) + vol;
                      ifsr=Core%Pin(ixy)%FsrIdxst+ilocalfsr-1
                      !write(102,'(5i4,ES16.6)') ixy,icel,ifxr,ilocalfsr,ifsr,vol
                      DO ig=1,ng
                          DO impho=1,nmpho
                              if (ilpho.eq.0) then
                                  regphimout(i)%phimout(io)%aziphiout(ig,impho) = regphimout(i)%phimout(io)%aziphiout(ig,impho) + vol * FMInfo%Phis(ifsr,iz,ig)
                              else
                                  regphimout(i)%phimout(io)%aziphiout(ig,impho) = regphimout(i)%phimout(io)%aziphiout(ig,impho) + vol * FMInfo%Phim(ordermap(impho,ilpho),ifsr,iz,ig)
                              endif
                          ENDDO ! DO impho=1,nmpho
                      ENDDO ! DO ig=1,ng
                  ENDDO ! DO ifsrinfxr = 1, nFsrInFxr
              ENDDO ! DO ifxr= 1, nfxr
          elseif (ixc.eq.0) then
              nxc=Core%AsyInfo(iasyType)%nx
              do ixc=1,nxc
                  ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc); ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc)
                  icel = Core%Pin(ixy)%Cell(iz); nfxr = Core%CellInfo(icel)%nFXR
                  DO ifxr= 1, nfxr
                      nFsrInFxr = Core%CellInfo(icel)%nFsrInFxr(ifxr)
                      DO ifsrinfxr = 1, nFsrInFxr
                          ilocalfsr=Core%CellInfo(icel)%MapFxr2FsrIdx(ifsrinfxr,ifxr)
                          vol = Core%CellInfo(icel)%vol(ilocalfsr); AREA(i) = Area(i) + vol;
                          ifsr=Core%Pin(ixy)%FsrIdxst+ilocalfsr-1
                          !write(102,'(5i4,ES16.6)') ixy,icel,ifxr,ilocalfsr,ifsr,vol
                          DO ig=1,ng
                              DO impho=1,nmpho
                                  if (ilpho.eq.0) then
                                      regphimout(i)%phimout(io)%aziphiout(ig,impho) = regphimout(i)%phimout(io)%aziphiout(ig,impho) + vol * FMInfo%Phis(ifsr,iz,ig)
                                  else
                                      regphimout(i)%phimout(io)%aziphiout(ig,impho) = regphimout(i)%phimout(io)%aziphiout(ig,impho) + vol * FMInfo%Phim(ordermap(impho,ilpho),ifsr,iz,ig)
                                  endif
                              ENDDO ! DO impho=1,nmpho
                          ENDDO ! DO ig=1,ng
                      ENDDO ! DO ifsrinfxr = 1, nFsrInFxr
                  ENDDO ! DO ifxr= 1, nfxr
              enddo
          else
              nyc=Core%AsyInfo(iasyType)%ny
              do iyc=1,nyc
                  ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc); ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc)
                  icel = Core%Pin(ixy)%Cell(iz); nfxr = Core%CellInfo(icel)%nFXR
                  DO ifxr= 1, nfxr
                      nFsrInFxr = Core%CellInfo(icel)%nFsrInFxr(ifxr)
                      DO ifsrinfxr = 1, nFsrInFxr
                          ilocalfsr=Core%CellInfo(icel)%MapFxr2FsrIdx(ifsrinfxr,ifxr)
                          vol = Core%CellInfo(icel)%vol(ilocalfsr); AREA(i) = Area(i) + vol;
                          ifsr=Core%Pin(ixy)%FsrIdxst+ilocalfsr-1
                          !write(102,'(5i4,ES16.6)') ixy,icel,ifxr,ilocalfsr,ifsr,vol
                          DO ig=1,ng
                              DO impho=1,nmpho
                                  if (ilpho.eq.0) then
                                      regphimout(i)%phimout(io)%aziphiout(ig,impho) = regphimout(i)%phimout(io)%aziphiout(ig,impho) + vol * FMInfo%Phis(ifsr,iz,ig)
                                  else
                                      regphimout(i)%phimout(io)%aziphiout(ig,impho) = regphimout(i)%phimout(io)%aziphiout(ig,impho) + vol * FMInfo%Phim(ordermap(impho,ilpho),ifsr,iz,ig)
                                  endif
                              ENDDO ! DO impho=1,nmpho
                          ENDDO ! DO ig=1,ng
                      ENDDO ! DO ifsrinfxr = 1, nFsrInFxr
                  ENDDO ! DO ifxr= 1, nfxr
              enddo
          endif
      ENDDO ! DO io=1,nphiout
      DO io=1,nphiout
          DO impho=1,regphimout(i)%phimout(io)%nazimom
              DO ig=1,ng
                 regphimout(i)%phimout(io)%aziphiout(ig,impho)=regphimout(i)%phimout(io)%aziphiout(ig,impho)/AREA(i)
              ENDDO
          ENDDO
      ENDDO
  ELSEIF (nTracerCntl%OutpCntl%RegPhimOutFXR(i)) THEN
      ixc = nTracerCntl%OutpCntl%RegPhimOutList(4, i); iyc = nTracerCntl%OutpCntl%RegPhimOutList(5, i)
      ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc); ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc)
      icel = Core%Pin(ixy)%Cell(iz); nfxr = Core%CellInfo(icel)%nFXR;
      DO io=1,nphiout
          ilpho=nTracerCntl%OutpCntl%PhimOrdOutList(io, i)
          nmpho=nordermap(ilpho); regphimout(i)%phimout(io)%nazimom=nmpho
          ALLOCATE(regphimout(i)%phimout(io)%aziphiout(ng,nmpho))
          regphimout(i)%phimout(io)%aziphiout=0
          AREA(i) = 0
          DO ir= nTracerCntl%OutpCntl%RegPhimOutList(6, i), nTracerCntl%OutpCntl%RegPhimOutList(7, i)
              nFsrInFxr = Core%CellInfo(icel)%nFsrInFxr(ir)
              DO ifsrinfxr = 1, nFsrInFxr
                  ilocalfsr = Core%CellInfo(icel)%MapFxr2FsrIdx(ifsrinfxr,ir)
                  vol = Core%CellInfo(icel)%vol(ilocalfsr); AREA(i) = Area(i) + vol;
                  ifsr=Core%Pin(ixy)%FsrIdxst+ilocalfsr-1
                  !write(103,'(5i4,ES16.6)') ixy,icel,ir,ilocalfsr,ifsr,vol
                  DO ig=1,ng
                      DO impho=1,nmpho
                          if (ilpho.eq.0) then
                              regphimout(i)%phimout(io)%aziphiout(ig,impho) = regphimout(i)%phimout(io)%aziphiout(ig,impho) + vol * FMInfo%Phis(ifsr,iz,ig)
                          else
                              regphimout(i)%phimout(io)%aziphiout(ig,impho) = regphimout(i)%phimout(io)%aziphiout(ig,impho) + vol * FMInfo%Phim(ordermap(impho,ilpho),ifsr,iz,ig)
                          endif
                      ENDDO ! DO impho=1,nmpho
                  ENDDO ! DO ig=1,ng
              ENDDO ! DO ifsrinfxr = 1, nFsrInFxr
          ENDDO ! DO ir= nTracerCntl%OutpCntl%RegXsOutList(6, i), nTracerCntl%OutpCntl%RegXsOutList(7, i)
      ENDDO ! DO io=1,nphiout
      DO io=1,nphiout
          DO impho=1,regphimout(i)%phimout(io)%nazimom
              DO ig=1,ng
                 regphimout(i)%phimout(io)%aziphiout(ig,impho)=regphimout(i)%phimout(io)%aziphiout(ig,impho)/AREA(i)
              ENDDO
          ENDDO
      ENDDO
  ELSE
      RETURN
  ENDIF
ENDDO

fn=trim(caseid)//'.phim'

CALL openfile(io16,FALSE,FALSE,FALSE, fn)
DO i = 1, nTracerCntl%OutpCntl%nRegPhimOut
    nphiout = nTracerCntl%OutpCntl%PhimOrdOutList(0, i)
#ifndef __GFORTRAN__
    IF (nTracerCntl%OutpCntl%RegPhimOutASM(i)) THEN
        WRITE(io16, '(A10, I5, A3, 2I5, A3, <nphiout>I5)') 'Region :',  nTracerCntl%OutpCntl%RegPhimOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(2:3, i), '/',  (nTracerCntl%OutpCntl%PhimOrdOutList(k, i),k=1,nphiout)
    ELSEIF (nTracerCntl%OutpCntl%RegPhimOutPIN(i)) THEN
        WRITE(io16, '(A10, I5, A3, 2(2I5, A3), <nphiout>I5)') 'Region :',  nTracerCntl%OutpCntl%RegPhimOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(4:5, i), '/', (nTracerCntl%OutpCntl%PhimOrdOutList(k, i),k=1,nphiout)
    ELSEIF (nTracerCntl%OutpCntl%RegPhimOutFXR(i)) THEN
        WRITE(io16, '(A10, I5, A3, 3(2I5, A3), <nphiout>I5)') 'Region :',  nTracerCntl%OutpCntl%RegPhimOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(4:5, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(6:7, i), '/', (nTracerCntl%OutpCntl%PhimOrdOutList(k, i),k=1,nphiout)
    ENDIF
#else
    READ(str2num, *) nphiout
    IF (nTracerCntl%OutpCntl%RegPhimOutASM(i)) THEN
        WRITE(io16, '(A10, I5, A3, 2I5, A3, ' // TRIM(str2num) // 'I5)') 'Region :',  nTracerCntl%OutpCntl%RegPhimOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(2:3, i), '/',  (nTracerCntl%OutpCntl%PhimOrdOutList(k, i),k=1,nphiout)
    ELSEIF (nTracerCntl%OutpCntl%RegPhimOutPIN(i)) THEN
        WRITE(io16, '(A10, I5, A3, 2(2I5, A3), ' // TRIM(str2num) // 'I5)') 'Region :',  nTracerCntl%OutpCntl%RegPhimOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(4:5, i), '/', (nTracerCntl%OutpCntl%PhimOrdOutList(k, i),k=1,nphiout)
    ELSEIF (nTracerCntl%OutpCntl%RegPhimOutFXR(i)) THEN
        WRITE(io16, '(A10, I5, A3, 3(2I5, A3), ' // TRIM(str2num) // 'I5)') 'Region :',  nTracerCntl%OutpCntl%RegPhimOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(4:5, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhimOutList(6:7, i), '/', (nTracerCntl%OutpCntl%PhimOrdOutList(k, i),k=1,nphiout)
    ENDIF
#endif
    WRITE(io16, '(A10, ES16.6)') 'AREA :',Area(i)
    WRITE(io16, '(A5)', ADVANCE='NO') 'GRP'
    DO io=1,nTracerCntl%OutpCntl%PhimOrdOutList(0, i)
        k=nTracerCntl%OutpCntl%PhimOrdOutList(io, i)
        WRITE(io16, '(I16)', ADVANCE='NO') k
        DO j = 1, 16*(nordermap(k)-1)
           WRITE(io16, '(A)', ADVANCE='NO') ' '
        ENDDO
    ENDDO
    WRITE(io16,*)
    DO ig=1,ng
        WRITE(io16, '(I5)', ADVANCE='NO') ig
        DO io=1,nTracerCntl%OutpCntl%PhimOrdOutList(0, i)
            k=nTracerCntl%OutpCntl%PhimOrdOutList(io, i)
            DO impho=1,nordermap(k)
                WRITE(io16, '(ES16.6)', ADVANCE='NO') regphimout(i)%phimout(io)%aziphiout(ig,impho)
            ENDDO
        ENDDO
        WRITE(io16,*)
    ENDDO
ENDDO

close(io16)


END SUBROUTINE
