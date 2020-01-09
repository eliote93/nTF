module Xsec4MC
    real(8), private, pointer, dimension(:,:,:) :: SctCDF, xssm2g
    real(8), private, pointer, dimension(:,:) :: xsaa
    integer, pointer, dimension(:) :: g2gc
contains

subroutine InitXsec4MC(ngMC, gThrm)
    implicit none
    integer, intent(in) :: ngMC, gThrm
    ! required to be modified
    allocate(g2gc(ngMC))
    g2gc(1:gThrm)=1
    g2gc(gThrm+1:ngMC)=2
    ! MJ Fixing negatvie scattering cross section
!    call ConvNSct2Absr()
!    call adjXsec4MC()
!    call TuneXsec()
    call CalSctCDF() ! this first and
    call CalSct2G()  ! then this later
end subroutine

subroutine TuneXsec()
    use BenchXs
    implicit none
    integer :: i, g, gp
    real(8) :: sct
    do i=1, nxsltype
        do g=1, ngben
            do gp=1, ngben
                if (MacXsBen(i)%xss(g,gp)<0) MacXsBen(i)%xss(g,gp)=0
            enddo
            sct=sum(MacXsBen(i)%xss(g,:))
            MacXsBen(i)%xst(g)=MacXsBen(i)%xsa(g)+sct
        enddo
    enddo
end subroutine

subroutine SetG2GC(g2gc_) 
    implicit none
    integer, pointer, dimension(:) :: g2gc_
    g2gc_=g2gc
end subroutine

function GetG2GC(g) result(gc)
    implicit none
    integer :: gc
    integer, intent(in) :: g
    gc=g2gc(g)
end function

subroutine AdrSctCDF(xs, txs)
    use MCDefine
    implicit none
    type(XsecSet), intent(inout) :: xs
    integer, intent(in) :: txs
    xs%SctCDF=>SctCDF(:,:,txs)
end subroutine

subroutine AdrSct2G(xs, txs)
    use MCDefine
    implicit none
    type(XsecSet), intent(inout) :: xs
    integer, intent(in) :: txs
    xs%xssm2g=>xssm2g(:,:,txs)
end subroutine

subroutine AdrAbsrA(xs, txs)
    use MCDefine
    implicit none
    type(XsecSet), intent(inout) :: xs
    integer, intent(in) :: txs
    xs%xsaa=>xsaa(:,txs)
end subroutine

subroutine ConvNSct2Absr()
    use BenchXs
    implicit none
    integer :: i, g, gp
    real(8) :: xsct
    allocate(xsaa(ngben, nxsltype))
    do i=1, nxsltype
        do g=1, ngben
            do gp=1, ngben
                if (MacXsBen(i)%xss(g,gp)<0) then
                    if (g==gp) xsaa(g,i)=-MacXsBen(i)%xss(g,gp)
                    MacXsBen(i)%xss(g,gp)=0
                endif
            enddo
        enddo
    enddo
end subroutine

subroutine CalSctCDF()
    use BenchXs
    USE HighOrderSC, ONLY : SetAnIsoGaussian_MAC
    implicit none
    integer :: i, g, gp, ng, igs, ig
    real(8) :: xsct
    allocate(SctCDF(ngben,ngben,nxsltype))
    SctCDF=0.
    ng=ngben
    do i=1, nxsltype
        IF( benchxstype.EQ.4 )THEN  
            ALLOCATE(MacXsBen(i)%xss0(ng,ng))
            DO ig = 1, ng
                DO igs = 1, ng
                    MacXsBen(i)%xss0(igs,ig)=0
                ENDDO
                DO igs = MacXsBen(i)%PnSm(1,ig)%ib, MacXsBen(i)%PnSm(1,ig)%ie
                    MacXsBen(i)%xss0(igs,ig)=MacXsBen(i)%PnSm(1,ig)%from(igs)
                ENDDO
            ENDDO
            MacXsBen(i)%xss=>MacXsBen(i)%xss0
            IF( scatod .GE. 1 )THEN
                ALLOCATE(MacXsBen(i)%xss1(ng,ng))
                DO ig = 1, ng
                    DO igs = 1, ng
                        MacXsBen(i)%xss1(igs,ig)=0
                    ENDDO
                    DO igs = MacXsBen(i)%PnSm(2,ig)%ib, MacXsBen(i)%PnSm(2,ig)%ie
                        MacXsBen(i)%xss1(igs,ig)=MacXsBen(i)%PnSm(2,ig)%from(igs)
                    ENDDO
                    !MacXsBen(i)%xss0(ig,ig)=MacXsBen(i)%xss0(ig,ig)+MacXsBen(i)%xs0sum(ig)-MacXsBen(i)%xs1sum(ig) !regain total scattering
                ENDDO
                IF( scatod .GE. 2 )THEN
                    ALLOCATE(MacXsBen(i)%xss2(ng,ng))
                    DO ig = 1, ng
                        DO igs = 1, ng
                            MacXsBen(i)%xss2(igs,ig)=0
                        ENDDO
                        DO igs = MacXsBen(i)%PnSm(3,ig)%ib, MacXsBen(i)%PnSm(3,ig)%ie
                            MacXsBen(i)%xss2(igs,ig)=MacXsBen(i)%PnSm(3,ig)%from(igs)
                        ENDDO
                    ENDDO
                    IF( scatod .GE. 3 )THEN
                        ALLOCATE(MacXsBen(i)%xss3(ng,ng))
                        DO ig = 1, ng
                            DO igs = 1, ng
                                MacXsBen(i)%xss3(igs,ig)=0
                            ENDDO
                            DO igs = MacXsBen(i)%PnSm(4,ig)%ib, MacXsBen(i)%PnSm(4,ig)%ie
                                MacXsBen(i)%xss3(igs,ig)=MacXsBen(i)%PnSm(4,ig)%from(igs)
                            ENDDO
                        ENDDO
                    ENDIF
                ENDIF
            ENDIF           
        ENDIF
        CALL SetAnIsoGaussian_MAC(MacXsBen(i), scatod)
        
        do g=1, ngben
            xsct=MacXsBen(i)%xss0(g,1)
            if (xsct<0) then
                xsct=0.
            endif
            
            SctCDF(1,g,i)=xsct
            do gp=2, ngben
                xsct=MacXsBen(i)%xss0(g,gp)
                if (xsct<0) then
                    xsct=0.
                endif
                
                SctCDF(gp,g,i)=SctCDF(gp-1,g,i)+xsct
            enddo
            SctCDF(:,g,i)=SctCDF(:,g,i)/SctCDF(ngben,g,i)
        enddo
        !DEALLOCATE(MacXsBen(i)%xss0)
    enddo
end subroutine



subroutine CalSct2G() ! for 2G CMFD
    use BenchXs
    implicit none
    integer :: i, g, gp
    
    allocate(xssm2g(ngben,2,nxsltype))
    xssm2g=0.
    do i=1, nxsltype
        do g=1, ngben
            do gp=1, ngben
                xssm2g(g,g2gc(gp),i)=xssm2g(g,g2gc(gp),i)+MacXsBen(i)%xss(g,gp)
            enddo
        enddo
    enddo
end subroutine

subroutine adjXsec4MC
    use BenchXs
    implicit none
    integer :: i, g, gp
    
1010 format ("Negative Self Scattering Fixed for iso", 1x, i2, 2x, "group", 1x, i2, ".")
1020 format ("Negative Scattering Fixed for iso", 1x, i2, "g(", i2, ") -> g(",i2,")")
    do i=1, nxsltype
        do g=1, ngben
            do gp=1, ngben
                if (MacXsBen(i)%xss(g,gp)<0) then
                    if (g == gp) then
                        MacXsBen(i)%xst(g)=MacXsBen(i)%xst(g)-MacXsBen(i)%xss(g,gp)
!                        MacXsBen(i)%xsa(g)=MacXsBen(i)%xsa(g)-MacXsBen(i)%xss(g,gp)
                        MacXsBen(i)%xss(g,gp)=0
                        write(*, 1010) i, g
                    else
                        MacXsBen(i)%xss(g,gp)=0
                        write(*, 1020) i, g, gp
                    endif
                endif
            enddo
        enddo
    enddo
end subroutine
    
end module
