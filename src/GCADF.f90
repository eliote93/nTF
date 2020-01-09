
SUBROUTINE GCSurf(Core, FmInfo, ng)
    USE TYPEDEF,         ONLY : CoreInfo_Type, FmInfo_Type
    USE GC_mod,          ONLY : GrpIdx, nGrp, iasy, Cell, Asy, AsyInfo, &
                                nx, ny, nxy, lZdir, iasytype, iz, SurfJG, SurfJ, SurfPhiG, hz, kinf
    USE GCpin_mod
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FmInfo_Type) :: FmInfo
    INTEGER :: ng
    
    REAL, POINTER :: CellSurfJ(:,:,:,:,:), CellSurfPhi(:,:,:,:) !(innout2, nbd4, nx, ny, ngrp), (nbd4, nx, ny, ngrp)
    REAL, POINTER :: CellSurfJG(:,:,:,:,:), CellSurfPhiG(:,:,:,:) !CellSurfJG(innout2, nbd4, nx, ny, ngrp)
    INTEGER :: ninnout, nbd
    INTEGER :: innout, i, ig, gidx, gidx2, ix, iy, ixy, ixyl, icell
    REAL :: wsum, w
    INTEGER :: ig2
    
    ninnout=2; nbd=4  
    !CellSurfJ=>Cminfo%radJout
    CellSurfJ=>Fminfo%radJout(1:2,:,:,:,:) !ninnout/ 1:in 2:out 3:surfphi
    CellSurfPhi=>FmInfo%radJout(3,:,:,:,:)
    
    ALLOCATE(CellSurfJG(0:ninnout, nbd, nx, ny, ngrp)) ! 0 for out
    ALLOCATE(CellSurfPhiG(nbd, nx, ny, ngrp)) ! 0 for out
    CellSurfJG=0;CellSurfPhiG=0;
    
    !----------Assembly surfJ in 2G    
    DO ig = 1, ng
        gidx=GrpIdx(ig)
        DO ix = 1, nx
            DO iy = 1, ny
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
                ixy = Asy(iasy)%GlobalPinIdx(ixyl)
                !---
                pinR(ixy,ig)=pinFis(ixy,ig)/kinf
                DO ig2 = 1, ng
                    gidx2=GrpIdx(ig2)
                    pinR(ixy,ig)=pinR(ixy,ig)+pinSS(ixy,ig2,ig)
                    IF( gidx2 .NE. gidx )THEN
                        pinSoutG(ixy,ig)=pinSoutG(ixy,ig)+pinSS(ixy,ig,ig2)
                        pinSinG(ixy,ig)=pinSinG(ixy,ig)+pinSS(ixy,ig2,ig)
                    ENDIF
                    pinSout(ixy,ig)=pinSout(ixy,ig)+pinSS(ixy,ig,ig2)
                    pinSin(ixy,ig)=pinSin(ixy,ig)+pinSS(ixy,ig2,ig)
                    
                    IF( ig2 .NE. ig )THEN
                    !    pinR(ixyl,ig)=pinR(ixyl,ig)+pinSS(ixyl,ig2,ig)
                    ENDIF
                ENDDO
                !===
                !srcr(ixy,ig)=pinR(ixy,ig)/pinfisr(ixy,ig)/hz
                !pinR(ixy,ig)=pinR(ixy,ig)-pinRmv(ixy,ig)
                !DO innout = 1, ninnout
                !    DO i = 1, nbd
                !        CellSurfJG(innout, i, ix, iy, gidx)=CellSurfJG(innout, i, ix, iy, gidx)+CellSurfJ(innout, i, ixy, iz, ig)
                !        !---
                !        if( innout .EQ. 1 )THEN
                !            pinJ(ixy,ig)=pinJ(ixy,ig)-CellSurfJ(innout, i, ixy, iz, ig)
                !            pinJdir(ixy,ig,i)=pinJdir(ixy,ig,i)-CellSurfJ(innout, i, ixy, iz, ig)
                !        else
                !            pinJ(ixy,ig)=pinJ(ixy,ig)+CellSurfJ(innout, i, ixy, iz, ig)
                !            pinJdir(ixy,ig,i)=pinJdir(ixy,ig,i)+CellSurfJ(innout, i, ixy, iz, ig)
                !        endif
                !        !===
                !    ENDDO
                !ENDDO
                DO i = 1, nbd
                    CellSurfPhiG(i, ix, iy, gidx)=CellSurfPhiG(i, ix, iy, gidx)+CellSurfPhi(i, ixy, iz, ig)
                    DO innout = 1, ninnout
                        CellSurfJG(innout, i, ix, iy, gidx)=CellSurfJG(innout, i, ix, iy, gidx)+CellSurfJ(innout, i, ixy, iz, ig)
                        !---
                        if( innout .EQ. 1 )THEN
                            pinJ(ixy,ig)=pinJ(ixy,ig)-CellSurfJ(innout, i, ixy, iz, ig)
                            pinJdir(ixy,ig,i)=pinJdir(ixy,ig,i)-CellSurfJ(innout, i, ixy, iz, ig)
                        else
                            pinJ(ixy,ig)=pinJ(ixy,ig)+CellSurfJ(innout, i, ixy, iz, ig)
                            pinJdir(ixy,ig,i)=pinJdir(ixy,ig,i)+CellSurfJ(innout, i, ixy, iz, ig)
                        endif
                        !===
                    ENDDO
                ENDDO
                !---
                !pinR(ixyl,ig)=pinR(ixyl,ig)/pinJ(ixyl,ig)/hz
                !pinR(ixy,ig)=pinR(ixy,ig)/hz-pinJ(ixy,ig)
                !lossr(ixy,ig)=pinrmv(ixy,ig)/pinrmvr(ixy,ig)/hz
                !PinJR(ixy,ig)=PinJ(ixy,ig)/pinJR(ixy,ig)
                !===
            ENDDO
        ENDDO
    ENDDO
    
    SurfJ=0;
    DO ig = 1, ng
        wsum=0
        DO i = 1, 4
            wsum=0
            SELECT CASE(i)
            CASE(1) !---SOUTH
                iy=ny
                DO ix = 1, nx
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LX
                    !w=1                    
                    DO innout = 1, ninnout
                        SurfJ(innout, i, ig)=SurfJ(innout, i, ig)+CellSurfJ(innout, i, ixy, iz, ig)!*w
                    ENDDO
                    wsum=wsum+w
                ENDDO
            CASE(2) !---WEST
                ix=1
                DO iy = 1, ny
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LY
                    !w=1
                    DO innout = 1, ninnout
                        SurfJ(innout, i, ig)=SurfJ(innout, i, ig)+CellSurfJ(innout, i, ixy, iz, ig)!*w
                    ENDDO
                    wsum=wsum+w
                ENDDO
            CASE(3) !---NORTH
                iy=1
                DO ix = 1, nx
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LX
                    !w=1
                    DO innout = 1, ninnout
                        SurfJ(innout, i, ig)=SurfJ(innout, i, ig)+CellSurfJ(innout, i, ixy, iz, ig)!*w
                    ENDDO
                    wsum=wsum+w
                ENDDO
            CASE(4) !---EAST
                ix=nx
                DO iy = 1, ny
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LY
                    !w=1
                    DO innout = 1, ninnout
                        SurfJ(innout, i, ig)=SurfJ(innout, i, ig)+CellSurfJ(innout, i, ixy, iz, ig)!*w
                    ENDDO
                    wsum=wsum+w
                ENDDO                
            ENDSELECT
            SurfJ(0, i, ig)=SurfJ(2, i, ig)-SurfJ(1, i, ig)
        ENDDO !--- Dir sweep    
        DO i=1,4
            DO innout = 0, ninnout
                SurfJ(innout, i, ig)=SurfJ(innout, i, ig)*hz
            ENDDO
        ENDDO    
    ENDDO !--- group sweep
    
    
    SurfJG=0;SurfPhiG=0;
    DO gidx = 1, ngrp
        wsum=0
        DO i = 1, 4
            wsum=0
            SELECT CASE(i)
            CASE(1) !---SOUTH
                iy=ny
                DO ix = 1, nx
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LX
                    !w=1                    
                    DO innout = 1, ninnout
                        SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)+CellSurfJG(innout, i, ix, iy, gidx)!*w
                    ENDDO
                    SurfPhiG(i, gidx)=SurfPhiG(i, gidx)+CellSurfPhiG(i, ix, iy, gidx)
                    wsum=wsum+w
                ENDDO
            CASE(2) !---WEST
                ix=1
                DO iy = 1, ny
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LY
                    !w=1
                    DO innout = 1, ninnout
                        SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)+CellSurfJG(innout, i, ix, iy, gidx)!*w
                    ENDDO
                    SurfPhiG(i, gidx)=SurfPhiG(i, gidx)+CellSurfPhiG(i, ix, iy, gidx)
                    wsum=wsum+w
                ENDDO
            CASE(3) !---NORTH
                iy=1
                DO ix = 1, nx
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LX
                    !w=1
                    DO innout = 1, ninnout
                        SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)+CellSurfJG(innout, i, ix, iy, gidx)!*w
                    ENDDO
                    SurfPhiG(i, gidx)=SurfPhiG(i, gidx)+CellSurfPhiG(i, ix, iy, gidx)
                    wsum=wsum+w
                ENDDO
            CASE(4) !---EAST
                ix=nx
                DO iy = 1, ny
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LY
                    !w=1
                    DO innout = 1, ninnout
                        SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)+CellSurfJG(innout, i, ix, iy, gidx)!*w
                    ENDDO
                    SurfPhiG(i, gidx)=SurfPhiG(i, gidx)+CellSurfPhiG(i, ix, iy, gidx)
                    wsum=wsum+w
                ENDDO                
            ENDSELECT
            DO innout = 1, ninnout
                !SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)/wsum
            ENDDO
            SurfPhiG(i, gidx)=SurfPhiG(i, gidx)/wsum
            SurfJG(0, i, gidx)=SurfJG(2, i, gidx)-SurfJG(1, i, gidx)
        ENDDO !--- Dir sweep    
        DO i=1,4
            DO innout = 0, ninnout
                !SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)/wsum
                SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)*hz
            ENDDO
        ENDDO    
    ENDDO !--- gidx sweep
    
    CONTINUE
    DEALLOCATE(CellSurfJG, CellSurfPhiG)    
    
    
ENDSUBROUTINE

SUBROUTINE GCSurf_Pin(Core, FmInfo, ng)
    USE TYPEDEF,         ONLY : CoreInfo_Type, FmInfo_Type
    USE GC_mod,          ONLY : GrpIdx, nGrp, iasy, Cell, Asy, AsyInfo, &
                                nx, ny, nxy, lZdir, iasytype, iz, SurfJG, SurfJ, SurfPhiG, hz, kinf,&
                                xbg, xed, yst, yed
    USE GCpin_mod
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FmInfo_Type) :: FmInfo
    INTEGER :: ng
    
    REAL, POINTER :: CellSurfJ(:,:,:,:,:), CellSurfPhi(:,:,:,:) !(innout2, nbd4, nx, ny, ngrp), (nbd4, nx, ny, ngrp)
    REAL, POINTER :: CellSurfJG(:,:,:,:,:), CellSurfPhiG(:,:,:,:) !CellSurfJG(innout2, nbd4, nx, ny, ngrp)
    INTEGER :: ninnout, nbd
    INTEGER :: innout, i, ig, gidx, gidx2, ix, iy, ixy, ixyl, icell
    REAL :: wsum, w
    INTEGER :: ig2
    
    ninnout=2; nbd=4  
    !CellSurfJ=>Cminfo%radJout
    CellSurfJ=>Fminfo%radJout(1:2,:,:,:,:) !ninnout/ 1:in 2:out 3:surfphi
    CellSurfPhi=>FmInfo%radJout(3,:,:,:,:)
    
    ALLOCATE(CellSurfJG(0:ninnout, nbd, nx, ny, ngrp)) ! 0 for out
    ALLOCATE(CellSurfPhiG(nbd, nx, ny, ngrp)) ! 0 for out
    CellSurfJG=0;CellSurfPhiG=0;
    
    !----------Assembly surfJ in 2G    
    DO ig = 1, ng
        gidx=GrpIdx(ig)
        DO ix = xbg, xed
            DO iy = yst, yed
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
                ixy = Asy(iasy)%GlobalPinIdx(ixyl)
                !---
                pinR(ixy,ig)=pinFis(ixy,ig)/kinf
                DO ig2 = 1, ng
                    gidx2=GrpIdx(ig2)
                    pinR(ixy,ig)=pinR(ixy,ig)+pinSS(ixy,ig2,ig)
                    IF( gidx2 .NE. gidx )THEN
                        pinSoutG(ixy,ig)=pinSoutG(ixy,ig)+pinSS(ixy,ig,ig2)
                        pinSinG(ixy,ig)=pinSinG(ixy,ig)+pinSS(ixy,ig2,ig)
                    ENDIF
                    pinSout(ixy,ig)=pinSout(ixy,ig)+pinSS(ixy,ig,ig2)
                    pinSin(ixy,ig)=pinSin(ixy,ig)+pinSS(ixy,ig2,ig)
                    
                    IF( ig2 .NE. ig )THEN
                    !    pinR(ixyl,ig)=pinR(ixyl,ig)+pinSS(ixyl,ig2,ig)
                    ENDIF
                ENDDO
                !===
                !srcr(ixy,ig)=pinR(ixy,ig)/pinfisr(ixy,ig)/hz
                !pinR(ixy,ig)=pinR(ixy,ig)-pinRmv(ixy,ig)
                !DO innout = 1, ninnout
                !    DO i = 1, nbd
                !        CellSurfJG(innout, i, ix, iy, gidx)=CellSurfJG(innout, i, ix, iy, gidx)+CellSurfJ(innout, i, ixy, iz, ig)
                !        !---
                !        if( innout .EQ. 1 )THEN
                !            pinJ(ixy,ig)=pinJ(ixy,ig)-CellSurfJ(innout, i, ixy, iz, ig)
                !            pinJdir(ixy,ig,i)=pinJdir(ixy,ig,i)-CellSurfJ(innout, i, ixy, iz, ig)
                !        else
                !            pinJ(ixy,ig)=pinJ(ixy,ig)+CellSurfJ(innout, i, ixy, iz, ig)
                !            pinJdir(ixy,ig,i)=pinJdir(ixy,ig,i)+CellSurfJ(innout, i, ixy, iz, ig)
                !        endif
                !        !===
                !    ENDDO
                !ENDDO
                DO i = 1, nbd
                    CellSurfPhiG(i, ix, iy, gidx)=CellSurfPhiG(i, ix, iy, gidx)+CellSurfPhi(i, ixy, iz, ig)
                    DO innout = 1, ninnout
                        CellSurfJG(innout, i, ix, iy, gidx)=CellSurfJG(innout, i, ix, iy, gidx)+CellSurfJ(innout, i, ixy, iz, ig)
                        !---
                        if( innout .EQ. 1 )THEN
                            pinJ(ixy,ig)=pinJ(ixy,ig)-CellSurfJ(innout, i, ixy, iz, ig)
                            pinJdir(ixy,ig,i)=pinJdir(ixy,ig,i)-CellSurfJ(innout, i, ixy, iz, ig)
                        else
                            pinJ(ixy,ig)=pinJ(ixy,ig)+CellSurfJ(innout, i, ixy, iz, ig)
                            pinJdir(ixy,ig,i)=pinJdir(ixy,ig,i)+CellSurfJ(innout, i, ixy, iz, ig)
                        endif
                        !===
                    ENDDO
                ENDDO
                !---
                !pinR(ixyl,ig)=pinR(ixyl,ig)/pinJ(ixyl,ig)/hz
                !pinR(ixy,ig)=pinR(ixy,ig)/hz-pinJ(ixy,ig)
                !lossr(ixy,ig)=pinrmv(ixy,ig)/pinrmvr(ixy,ig)/hz
                !PinJR(ixy,ig)=PinJ(ixy,ig)/pinJR(ixy,ig)
                !===
            ENDDO
        ENDDO
    ENDDO
    
    SurfJ=0;
    DO ig = 1, ng
        wsum=0
        DO i = 1, 4
            wsum=0
            SELECT CASE(i)
            CASE(1) !---SOUTH
                iy=yed
                DO ix = xbg, xed
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LX
                    !w=1                    
                    DO innout = 1, ninnout
                        SurfJ(innout, i, ig)=SurfJ(innout, i, ig)+CellSurfJ(innout, i, ixy, iz, ig)!*w
                    ENDDO
                    wsum=wsum+w
                ENDDO
            CASE(2) !---WEST
                ix=xbg
                DO iy = yst, yed
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LY
                    !w=1
                    DO innout = 1, ninnout
                        SurfJ(innout, i, ig)=SurfJ(innout, i, ig)+CellSurfJ(innout, i, ixy, iz, ig)!*w
                    ENDDO
                    wsum=wsum+w
                ENDDO
            CASE(3) !---NORTH
                iy=yst
                DO ix = xbg, xed
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LX
                    !w=1
                    DO innout = 1, ninnout
                        SurfJ(innout, i, ig)=SurfJ(innout, i, ig)+CellSurfJ(innout, i, ixy, iz, ig)!*w
                    ENDDO
                    wsum=wsum+w
                ENDDO
            CASE(4) !---EAST
                ix=xed
                DO iy = yst, yed
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LY
                    !w=1
                    DO innout = 1, ninnout
                        SurfJ(innout, i, ig)=SurfJ(innout, i, ig)+CellSurfJ(innout, i, ixy, iz, ig)!*w
                    ENDDO
                    wsum=wsum+w
                ENDDO                
            ENDSELECT
            SurfJ(0, i, ig)=SurfJ(2, i, ig)-SurfJ(1, i, ig)
        ENDDO !--- Dir sweep    
        DO i=1,4
            DO innout = 0, ninnout
                SurfJ(innout, i, ig)=SurfJ(innout, i, ig)*hz
            ENDDO
        ENDDO    
    ENDDO !--- group sweep
    
    
    SurfJG=0;SurfPhiG=0;
    DO gidx = 1, ngrp
        wsum=0
        DO i = 1, 4
            wsum=0
            SELECT CASE(i)
            CASE(1) !---SOUTH
                iy=yed
                DO ix = xbg, xed
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LX
                    !w=1                    
                    DO innout = 1, ninnout
                        SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)+CellSurfJG(innout, i, ix, iy, gidx)!*w
                    ENDDO
                    SurfPhiG(i, gidx)=SurfPhiG(i, gidx)+CellSurfPhiG(i, ix, iy, gidx)
                    wsum=wsum+w
                ENDDO
            CASE(2) !---WEST
                ix=xbg
                DO iy = yst, yed
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LY
                    !w=1
                    DO innout = 1, ninnout
                        SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)+CellSurfJG(innout, i, ix, iy, gidx)!*w
                    ENDDO
                    SurfPhiG(i, gidx)=SurfPhiG(i, gidx)+CellSurfPhiG(i, ix, iy, gidx)
                    wsum=wsum+w
                ENDDO
            CASE(3) !---NORTH
                iy=yst
                DO ix = xbg, xed
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LX
                    !w=1
                    DO innout = 1, ninnout
                        SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)+CellSurfJG(innout, i, ix, iy, gidx)!*w
                    ENDDO
                    SurfPhiG(i, gidx)=SurfPhiG(i, gidx)+CellSurfPhiG(i, ix, iy, gidx)
                    wsum=wsum+w
                ENDDO
            CASE(4) !---EAST
                ix=xed
                DO iy = yst, yed
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    icell=Core%Pin(ixy)%Cell(iz)
                    w=Cell(icell)%Geom%LY
                    !w=1
                    DO innout = 1, ninnout
                        SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)+CellSurfJG(innout, i, ix, iy, gidx)!*w
                    ENDDO
                    SurfPhiG(i, gidx)=SurfPhiG(i, gidx)+CellSurfPhiG(i, ix, iy, gidx)
                    wsum=wsum+w
                ENDDO                
            ENDSELECT
            DO innout = 1, ninnout
                !SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)/wsum
            ENDDO
            SurfPhiG(i, gidx)=SurfPhiG(i, gidx)/wsum
            SurfJG(0, i, gidx)=SurfJG(2, i, gidx)-SurfJG(1, i, gidx)
        ENDDO !--- Dir sweep    
        DO i=1,4
            DO innout = 0, ninnout
                !SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)/wsum
                SurfJG(innout, i, gidx)=SurfJG(innout, i, gidx)*hz
            ENDDO
        ENDDO    
    ENDDO !--- gidx sweep
    
    CONTINUE
    DEALLOCATE(CellSurfJG, CellSurfPhiG)    
    
    
ENDSUBROUTINE

SUBROUTINE GCADF(Core, ng)
    USE TYPEDEF,         ONLY : CoreInfo_Type
    USE GC_mod,          ONLY : GrpIdx, iasy, Cell, Asy, AsyInfo, PinPhiVolG, ADF2G, AvgADF2G, ADF2GSurf, &
                                nx, ny, nxy, PinPhiVol, nGrp, AvgBPhiG, PhiG, SurfPhiG, lZdir, iasytype, iz, nExtPol, lpin
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    INTEGER :: ng
    
    INTEGER :: i, j, k
    INTEGER :: ig
    INTEGER :: ix, iy, ixy, ixyl, icell
    INTEGER, POINTER :: pidx(:,:)   ! 2 = two point exptrapolation !!!WARNING : MIGHT NOT VALID FOR WATER GAP ASY > MODIFIED in r532
    INTEGER ::  gidx, gidx2
    REAL, POINTER :: tx(:), w(:), xx(:), yy(:)
    REAL :: wsum
    REAL :: PhiVolSum, VolSum
    ADF2G=1
    avgADF2g=1
    IF( .NOT. lPin )THEN
    !!---ADF for single assembly-----------------------------------------------------
    ALLOCATE(PinPhiVolG(ngrp,nx,ny))
    ALLOCATE(pidx(nExtPol,nxy))
    ALLOCATE(tx(nxy), w(nxy))
    ALLOCATE(xx(nExtPol), yy(nExtPol))
    
    !----------Assembly Discontinuity Factor in 2-Groups
    PinPhiVolg=0
    DO ig = 1, ng
        gidx=GrpIdx(ig)
        DO ix = 1, nx
            DO iy = 1, ny
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
                ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
                PinPhiVolg(gidx,ix,iy) = PinPhiVolg(gidx,ix,iy) + PinPhiVol(ig,ixy)
            ENDDO
        ENDDO
    ENDDO
    VolSum=0
    DO ix = 1, nx
        DO iy = 1, ny
            ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
            VolSum = VolSum + Core%PinVol(ixy,iz)
        ENDDO
    ENDDO
    
    DO gidx = 1, ngrp
        PhiVolSum=0
        DO ix = 1, nx
            DO iy = 1, ny
                PhiVolSum=PhiVolSum + PinPhiVolg(gidx,ix,iy)
            ENDDO
        ENDDO
        PhiVolSum=VolSum/PhiVolSum
        DO ix = 1, nx
            DO iy = 1, ny
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
                ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
                PinPhiVolg(gidx,ix,iy)=PinPhiVolg(gidx,ix,iy)*PhiVolSum/Core%PinVol(ixy,iz)
            ENDDO
        ENDDO
        DO i = 1, 4
            ADF2GSurf(i, gidx)=SurfPhiG(i, gidx)*PhiVolSum
        ENDDO        
    ENDDO
    
    if( nx .EQ. 1 )then
        nextpol=1;
    endif
    if( ny .EQ. 1 )then
        nextpol=1;
    endif
    
    if( nx .LT. nExtPol )THEN
        nextpol=2
    endif
    
    
    DO gidx = 1, ngrp
    DO i = 1, 4
        SELECT CASE(i)
        CASE(1) !---SOUTH
            ix=1
            DO j = 1, nExtPol
                iy=ny+1-j
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                icell=Core%Pin(ixy)%Cell(iz)
                tx(j)=Cell(icell)%Geom%LY
            ENDDO
            iy=ny
            wsum=0;
            DO ix = 1, nx
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                icell=Core%Pin(ixy)%Cell(iz)
                w(ix)=Cell(icell)%Geom%LX
                wsum=wsum+w(ix)
            ENDDO
            j=1
            xx(j) = tx(j)/2
            DO j = 2, nExtPol
                xx(j) = xx(j-1) + (tx(j-1)+tx(j)) /2
            ENDDO
            
            yy=0
            DO j = 1, nExtPol
                iy=ny+1-j
                DO ix = 1, nx
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    yy(j)=yy(j)+PinPhiVolg(gidx, ix, iy)*w(ix)
                ENDDO
                yy(j)=yy(j)/wsum
            ENDDO
        CASE(2) !---WEST
            iy=1
            DO j = 1, nExtPol
                ix=j
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                icell=Core%Pin(ixy)%Cell(iz)
                tx(j)=Cell(icell)%Geom%LX
            ENDDO
            ix=1
            wsum=0;
            DO iy = 1, ny
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                icell=Core%Pin(ixy)%Cell(iz)
                w(iy)=Cell(icell)%Geom%LY
                wsum=wsum+w(iy)
            ENDDO
            j=1
            xx(j) = tx(j)/2
            DO j = 2, nExtPol
                xx(j) = xx(j-1) + (tx(j-1)+tx(j)) /2
            ENDDO
            
            yy=0
            DO j = 1, nExtPol
                ix=j
                DO iy = 1, ny
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    yy(j)=yy(j)+PinPhiVolg(gidx, ix, iy)*w(iy)
                ENDDO
                yy(j)=yy(j)/wsum
            ENDDO
        CASE(3) !---NORTH
            ix=1
            DO j = 1, nExtPol
                iy=j
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                icell=Core%Pin(ixy)%Cell(iz)
                tx(j)=Cell(icell)%Geom%LY
            ENDDO
            iy=1
            wsum=0;
            DO ix = 1, nx
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                icell=Core%Pin(ixy)%Cell(iz)
                w(ix)=Cell(icell)%Geom%LX
                wsum=wsum+w(ix)
            ENDDO
            j=1
            xx(j) = tx(j)/2
            DO j = 2, nExtPol
                xx(j) = xx(j-1) + (tx(j-1)+tx(j)) /2
            ENDDO
            
            yy=0
            DO j = 1, nExtPol
                iy=j
                DO ix = 1, nx
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    yy(j)=yy(j)+PinPhiVolg(gidx, ix, iy)*w(ix)
                ENDDO
                yy(j)=yy(j)/wsum
            ENDDO
        CASE(4) !---EAST
            iy=1
            DO j = 1, nExtPol
                ix=nx+1-j
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                icell=Core%Pin(ixy)%Cell(iz)
                tx(j)=Cell(icell)%Geom%LX
            ENDDO
            ix=nx
            wsum=0;
            DO iy = 1, ny
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                icell=Core%Pin(ixy)%Cell(iz)
                w(iy)=Cell(icell)%Geom%LY
                wsum=wsum+w(iy)
            ENDDO
            j=1
            xx(j) = tx(j)/2
            DO j = 2, nExtPol
                xx(j) = xx(j-1) + (tx(j-1)+tx(j)) /2
            ENDDO
            
            yy=0
            DO j = 1, nExtPol
                ix=nx+1-j
                DO iy = 1, ny
                    ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy) ! Local Pin Index
                    ixy = Asy(iasy)%GlobalPinIdx(ixyl)       ! Global Pin Index
                    yy(j)=yy(j)+PinPhiVolg(gidx, ix, iy)*w(iy)
                ENDDO
                yy(j)=yy(j)/wsum
            ENDDO
        ENDSELECT       
        
        SELECT CASE(nExtPol)
        CASE(2)
            ADF2G(i,gidx)=xx(2)*yy(1)-xx(1)*yy(2)
            ADF2G(i,gidx)=ADF2G(i,gidx)/(xx(2)-xx(1))
        case(3)
            ADF2G(i,gidx)=+(yy(1)/xx(1)/xx(1)-yy(2)/xx(2)/xx(2))*xx(1)*xx(2)/(xx(2)-xx(1))
            ADF2G(i,gidx)=-(yy(2)/xx(2)/xx(2)-yy(3)/xx(3)/xx(3))*xx(2)*xx(3)/(xx(3)-xx(2))+ADF2G(i,gidx)
            ADF2G(i,gidx)=ADF2G(i,gidx)*xx(1)*xx(2)*xx(3)/(xx(1)*xx(3)+xx(2)*xx(3)-xx(1)*xx(2)-xx(1)*xx(3))
        ENDSELECT
    ENDDO !--- Dir sweep
    
    ENDDO !--- gidx sweep
    
    ENDIF
    CONTINUE
    
    
    
ENDSUBROUTINE

SUBROUTINE GCADF_Pin(Core, ng)
    USE TYPEDEF,         ONLY : CoreInfo_Type
    USE GC_mod,          ONLY : GrpIdx, iasy, Cell, Asy, AsyInfo, PinPhiVolG, ADF2G, AvgADF2G, ADF2GSurf, &
                                nx, ny, nxy, PinPhiVol, nGrp, AvgBPhiG, PhiG, SurfPhiG, lZdir, iasytype, iz, nExtPol, lpin, &
                                xbg, xed, yst, yed
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    INTEGER :: ng
    
    INTEGER :: i, j, k
    INTEGER :: ig
    INTEGER :: ix, iy, ixy, ixyl, icell
    INTEGER, POINTER :: pidx(:,:)   ! 2 = two point exptrapolation !!!WARNING : MIGHT NOT VALID FOR WATER GAP ASY > MODIFIED in r532
    INTEGER ::  gidx, gidx2
    REAL, POINTER :: tx(:), w(:), xx(:), yy(:)
    REAL :: wsum
    REAL :: PhiVolSum, VolSum
    ADF2G=1
    avgADF2g=1
    !!---ADF for single assembly-----------------------------------------------------
    ALLOCATE(PinPhiVolG(ngrp,nx,ny))
    ALLOCATE(pidx(nExtPol,nxy))
    ALLOCATE(tx(nxy), w(nxy))
    ALLOCATE(xx(nExtPol), yy(nExtPol))
    
    !----------Assembly Discontinuity Factor in 2-Groups
    PinPhiVolg=0
    DO ig = 1, ng
        gidx=GrpIdx(ig)
        DO ix = xbg, xed
            DO iy = yst, yed
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
                ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
                PinPhiVolg(gidx,ix,iy) = PinPhiVolg(gidx,ix,iy) + PinPhiVol(ig,ixy)
            ENDDO
        ENDDO
    ENDDO
    VolSum=0
    DO ix = xbg, xed
        DO iy = yst, yed
            ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
            VolSum = VolSum + Core%PinVol(ixy,iz)
        ENDDO
    ENDDO
    
    DO gidx = 1, ngrp
        PhiVolSum=0
        DO ix = xbg, xed
            DO iy = yst, yed
                PhiVolSum=PhiVolSum + PinPhiVolg(gidx,ix,iy)
            ENDDO
        ENDDO
        PhiVolSum=VolSum/PhiVolSum
        DO ix = xbg, xed
            DO iy = yst, yed
                ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
                ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
                PinPhiVolg(gidx,ix,iy)=PinPhiVolg(gidx,ix,iy)*PhiVolSum/Core%PinVol(ixy,iz)
            ENDDO
        ENDDO
        DO i = 1, 4
            ADF2GSurf(i, gidx)=SurfPhiG(i, gidx)*PhiVolSum
            ADF2G(i,gidx)=1
        ENDDO        
    ENDDO
    
ENDSUBROUTINE





SUBROUTINE GCFluxNorm(ng)
    !--- Flux power normalization
    USE GC_mod,       ONLY : ngrp, PhiG, avgBPhiG, conPhiG, Phi, iasy, NormFac, Asy, SurfPhiG, ADF2GSurf, lpin, NormFac, PFac
    IMPLICIT NONE
    INTEGER :: ng
    
    INTEGER :: ig, gidx, i
    REAL ::  p_ratio, phisumG, phisum
    
    
    IF( .NOT. lPin )THEN
        phisumG=0;phisum=0;
        p_ratio=NormFac/Asy(iasy)%wt
        !--- Few G spectrum
        DO gidx = 1, ngrp
            phisumG=phisumG+PhiG(gidx)
            !PhiG(gidx)=PhiG(gidx)* p_ratio
            avgbPhiG(gidx)=avgbPhiG(gidx)*p_ratio
            DO i = 1, 3
                conPhiG(i,gidx)=conPhiG(i,gidx)*p_ratio
            ENDDO
            DO i = 1, 4
                SurfPhiG(i, gidx)=ADF2GSurf(i, gidx)*PhiG(gidx)
            ENDDO        
        ENDDO
        
        DO ig = 1, ng
            phisum=phisum+Phi(ig)
        ENDDO
        p_ratio=p_ratio*phisumG/phisum
        DO ig = 1, ng
            Phi(ig)=Phi(ig)*p_ratio 
        ENDDO
        PFac=p_ratio
    ELSE
        p_ratio=PFac
        DO ig = 1, ng
            Phi(ig)=Phi(ig)*p_ratio 
        ENDDO
    ENDIF    
    
ENDSUBROUTINE


SUBROUTINE GCJNorm(ng)
    !--- Current power normalization
    USE GC_mod,       ONLY : ngrp, PhiG, SurfJG, isoMacXs2g, ChiG, kinf !MAC : D, a, r, f, nf, kf
    IMPLICIT NONE
    INTEGER :: ng
    
    INTEGER :: ig, gidx, gidx2, i, inout
    REAL :: FS
    REAL :: netJ, refJ, Jratio
    
    FS=0
    DO gidx = 1, ngrp
        FS=FS+isoMacXs2G(0, 5, gidx)*PhiG(gidx)
    ENDDO
    
    DO gidx = 1, ngrp
        netJ=0
        DO i = 1, 4
            netJ=netJ+SurfJG(0, i, gidx)
        ENDDO
        refJ=ChiG(gidx)*FS/kinf
        IF( gidx.EQ.1 )THEN
            gidx2=2
        ELSE
            gidx2=1
        ENDIF
        refJ=refJ+(isoMacXs2G(0, 3, gidx2)-isoMacXs2G(0, 2, gidx2))*PhiG(gidx2)
        refJ=refJ-isoMacXs2G(0, 3, gidx)*PhiG(gidx)
        Jratio=refJ/netJ
        DO i = 1, 4
            DO inout = 0, 2
                !SurfJG(inout, i, gidx)=SurfJG(inout, i, gidx)*Jratio ! JA
            ENDDO
        ENDDO
    ENDDO
    
    
ENDSUBROUTINE



