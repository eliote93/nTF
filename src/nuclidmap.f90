#include <defines.h>
module nuclidmap_mod
USE XSLIB_MOD
USE ALLOCS
IMPLICIT NONE

CONTAINS
SUBROUTINE nuclidmap
!
! setup mapping vectors for nuclide IDs and temperatures
!        itempmap(i, j)    : mapping vector
!                           i = integer temperature
!                           j = isotope
!     ****************************************************************************
IMPLICIT NONE
logical lexiso
INTEGER :: I,  K,  K2,  KIR,  M,  INUCL,  IPOS,  IPOSISO,  IDISO, IG
INTEGER :: IT1,  IT2,  ITEMP1,  ITEMP2,  IDDIF,  IP1
INTEGER :: NID
REAL :: TEMPMAX,  YDI,  YDXE,  YDPM,  VAL
INTEGER :: imt

! All Nuclide List
call dmalloc(nuclidhel,nelthel)
do i=1,nelthel
    nuclidhel(i)=ldiso(i)%nid
enddo
do k=1,maxnid
  mapnucl(k) = 0
ENDDO
DO k = 1, nelthel
  mapnucl(nuclidhel(k)) = k
ENDDO
! List of nuclides that contained P1 scattering CX
np1hel=0
do i=1,nelthel
    if (ldiso(i)%np1temp>0) np1hel=np1hel+1
enddo
k=0
call dmalloc(idnp1hel,np1hel)
do i=1,nelthel
    if(ldiso(i)%np1temp>0) then
        k=k+1
        idnp1hel(k)=i
    endif
enddo
do k=1,maxnid
  mapnuclp13(k) = 0
ENDDO
DO k = 1, np1hel
  k2 = idnp1hel(k)
  mapnuclp13(nuclidhel(k2)) = k
ENDDO
! List of fissionable material(Depletable Material)
nfishel=0
do i=1,nelthel
    if(ldiso(i)%ifis.gt.0) nfishel=nfishel+1
enddo
call dmalloc(idfishel,nfishel)
k=0
do i=1,nelthel
    if(ldiso(i)%ifis.gt.0) then
        k=k+1
        idfishel(k)=i
    endif
enddo
do k=1,maxnid
  mapfis(k) = 0
ENDDO
DO k = 1, nfishel
  mapfis(nuclidhel(idfishel(k))) = k
ENDDO
! temperature mapping
tempmax = 0
call dmalloc(itempmap,nxtempmap,nelthel)
DO nid = 1, nelthel
    it1=ldiso(nid)%ntemp
    if(ldiso(nid)%temp(it1).gt.tempmax) tempmax=ldiso(nid)%temp(it1)
ENDDO
DO nid = 1, nelthel
    it1=ldiso(nid)%ntemp
    it2=0
  IF (it1.gt.1) then
        itemp2 = NINT(ldiso(nid)%temp(1))
    DO i = 1, itemp2
      itempmap(i, nid) = 1
    ENDDO
    DO m = 2, it1
      itemp1 = itemp2+1
            itemp2=NINT(ldiso(nid)%temp(m))
      DO i = itemp1, itemp2
        itempmap(i, nid) = m-1
      ENDDO
    ENDDO
    itemp1 = itemp2+1
    DO i = itemp1, nxtempmap
      itempmap(i, nid) = it1-1
    ENDDO
  else
    DO i = 1, nxtempmap
      itempmap(i, nid) = 1
    ENDDO
  endIF
ENDDO
call dmalloc(itempmapp1,nxtempmap,np1hel)
DO nid = 1, np1hel
    it1=ldiso(nid)%np1temp
    it2=0
  IF (it1.gt.1) then
        itemp2=nint(ldiso(nid)%p1temp(1))
    DO i = 1, itemp2
      itempmapp1(i, nid) = 1
    ENDDO
    DO m = 2, it1
      itemp1 = itemp2+1
            itemp2=nint(ldiso(nid)%p1temp(m))
      DO i = itemp1, itemp2
        itempmapp1(i, nid) = m-1
      ENDDO
    ENDDO
    itemp1 = itemp2+1
    DO i = itemp1, nxtempmap
      itempmapp1(i, nid) = it1-1
    ENDDO
  else
    DO i = 1, nxtempmap
      itempmapp1(i, nid) = 1
    ENDDO
  endIF
ENDDO
! N2N
n2nhel=0
do i=1,nelthel
    if(ldiso(i)%inmn .ge. 2) n2nhel=n2nhel+1
enddo
call dmalloc(idn2nhel,n2nhel)
allocate(xsn2nhel(n2nhel,noghel))
k=0
do i=1,nelthel
    if(ldiso(i)%inmn .ge. 2) then
        k=k+1
        idn2nhel(k)=i
        do ig=1,noghel
            xsn2nhel(k,ig)=ldiso(i)%sign2n(ig)
        enddo
    endif
enddo
call dmalloc(mapn2n,nelthel)
DO i = 1, n2nhel
  mapn2n(idn2nhel(i)) = i
ENDDO
! N3N
n3nhel=0
do i=1,nelthel
    if(ldiso(i)%inmn .eq. 3) n3nhel=n3nhel+1
enddo
call dmalloc(idn3nhel,n3nhel)
allocate(xsn3nhel(n3nhel,noghel))
k=0
do i=1,nelthel
    if(ldiso(i)%inmn .eq. 3) then
        k=k+1
        idn3nhel(k)=i
        do ig=1,noghel
            xsn3nhel(k,ig)=ldiso(i)%sign3n(ig)
        end do
    endif
enddo
call dmalloc(mapn3n,nelthel)
DO i = 1, n3nhel
  mapn3n(idn3nhel(i)) = i
ENDDO
! Photon Production Mapping...
! photon production Fission(imt=1)/Radioactive(imt=2)/Inelastic Scattering(imt=3)/non-elastic(imt=4)/
imt = 1
nphprod=0
DO imt = 1, 4
  DO i=1,nelthel
    if(ldiso(i)%lphoton(imt)) nphprod(imt)=nphprod(imt)+1
  END DO
END DO
call dmalloc(idphprod,MAXVAL(nphprod),4)
DO imt = 1, 4
  k=0
  DO i=1,nelthel
    IF(ldiso(i)%lphoton(imt)) THEN
        k=k+1
        idphprod(k,imt)=i
    endif
  END DO
END DO
mapnuclpp = 0
DO imt = 1, 4
  DO k = 1, nphprod(imt)
    mapnuclpp(nuclidhel(idphprod(k,imt)),imt) = k
  END DO
END DO
! burnable
nburhel=0
do i=1,nelthel
    if(ldiso(i)%ibur.gt.0) nburhel=nburhel+1
enddo
call dmalloc(idburhel,nburhel)
k=0
do i=1,nelthel
    if(ldiso(i)%ibur.gt.0) then
        k=k+1
        idburhel(k)=i
    endif
enddo
!About Decayhel
allocate(decayhel(nburhel)); decayhel=0
k=0
do i=1, nelthel
    if (ldiso(i)%ibur.gt.0) then
        k=k+1
        decayhel(k)=ldiso(i)%dcy
    endif
end do
!About Yieldhel
allocate(yieldhel(nfishel,nburhel)); yieldhel=0
k=0
do i=1, nelthel
    if (ldiso(i)%ifis.gt.0) then
        k=k+1
        ! for possible future use
        call dmalloc(ldiso(i)%yield,nburhel)
        ldiso(i)%yield=0.
        yieldhel(k,:)=ldiso(i)%yield(:)
    endif
end do
DO i = 1, nburhel
  IF(nuclidhel(idburhel(i)).eq.53635) iphel_i = i    ! search Iodine position
  IF(nuclidhel(idburhel(i)).eq.54635) iphel_xe = i   ! search Xenon position
  IF(nuclidhel(idburhel(i)).eq.61649) iphel_pm = i   ! search Samarium position
ENDDO
!
RETURN
END SUBROUTINE
!
SUBROUTINE yieldxesm(idiso, ydi, ydxe, ydpm)
INTEGER :: idiso, inucl
REAL :: ydi, ydxe, ydpm
inucl = mapfis(idiso)
ydi = yieldhel(inucl, iphel_i)
ydxe = yieldhel(inucl, iphel_xe)
ydpm = yieldhel(inucl, iphel_pm)
return
END SUBROUTINE
!
SUBROUTINE decayconst(idiso, val, ipos, lexiso)
INTEGER :: idiso, ipos, I
REAL :: VAL
LOGICAL :: lexiso
lexiso = .false.
DO i = 1, nburhel
  IF(nuclidhel(idburhel(i)).eq.idiso) then
    lexiso = .true.
    ipos = i
    val = decayhel(i)
    exit
  endIF
ENDDO
return
END SUBROUTINE
!
function iposiso(idiso)
IMPLICIT NONE
INTEGER :: iposiso
INTEGER :: idiso

iposiso = mapnucl(idiso)
end function
!
function nuclide(iso)
IMPLICIT NONE
INTEGER :: nuclide,  iso
nuclide = nuclidhel(iso)
return
end function
!
function iposnucl(idiso)
IMPLICIT NONE
INTEGER :: iposnucl
INTEGER :: idiso
iposnucl = mapnucl(idiso)
return
end function
!
function lExistIsotope(idiso)
IMPLICIT NONE
logical lExistIsotope
INTEGER :: idiso
lExistIsotope = .true.
IF(mapnucl(idiso).eq.0) lExistIsotope = .false.
return
end function
!
function lfissile(idiso)
IMPLICIT NONE
logical lfissile
INTEGER :: idiso
lfissile = .true.
IF(mapfis(idiso).eq.0) lfissile = .false.
return
end function
!
function AtomicWeight(idiso)
IMPLICIT NONE
REAL :: AtomicWeight
INTEGER :: idiso
!AtomicWeight = awhel(mapnucl(idiso))
AtomicWeight = ldiso(mapnucl(idiso))%aw
return
end function
!
function pndcrit(iso)
IMPLICIT NONE
real :: pndcrit
integer :: iso
pndcrit = ldiso(iso)%crit_nd
return
end function
!
function decayxe135()
IMPLICIT NONE
REAL(4) :: decayxe135
decayxe135 = decayhel(iphel_xe)
return
end function

function decayI135()
IMPLICIT NONE
REAL(4) :: decayI135
decayI135 = decayhel(iphel_i)
return
end function

subroutine InitResoCat()
USE XSLIB_MOD, ONLY : libdata, ldiso, CoreResIsoUpdate, nmaxz
IMPLICIT NONE
INTEGER :: niso,sum,iiso,i,ic,nid,idiso(100),xxx
TYPE(libdata),POINTER :: lib

#ifndef newcat
nCat=4
allocate(ResoCat(nCat))
allocate(ResoCatUse(nCat,nmaxz))
ResoCatUse=.false.

ResoCat(1)%repid=92238
ResoCat(2)%repid=92235
ResoCat(3)%repid=40000
ResoCat(4)%repid=72177

sum=0
niso=1; ResoCat(1)%niso=niso; sum=sum+niso
allocate(ResoCat(1)%idiso(niso))
niso=18; ResoCat(2)%niso=niso; sum=sum+niso
allocate(ResoCat(2)%idiso(niso))
niso=6; ResoCat(3)%niso=niso; sum=sum+niso
allocate(ResoCat(3)%idiso(niso))
niso=43; ResoCat(4)%niso=niso; sum=sum+niso
allocate(ResoCat(4)%idiso(niso))

iiso=1
ResoCat(1)%idiso(iiso)=92238

iiso=1
ResoCat(2)%idiso(iiso)=92235; iiso=iiso+1
ResoCat(2)%idiso(iiso)=90232; iiso=iiso+1
ResoCat(2)%idiso(iiso)=92233; iiso=iiso+1
ResoCat(2)%idiso(iiso)=92236; iiso=iiso+1
ResoCat(2)%idiso(iiso)=93237; iiso=iiso+1
ResoCat(2)%idiso(iiso)=94238; iiso=iiso+1
ResoCat(2)%idiso(iiso)=94239; iiso=iiso+1
ResoCat(2)%idiso(iiso)=94240; iiso=iiso+1
ResoCat(2)%idiso(iiso)=94241; iiso=iiso+1
ResoCat(2)%idiso(iiso)=94242; iiso=iiso+1
ResoCat(2)%idiso(iiso)=95241; iiso=iiso+1
ResoCat(2)%idiso(iiso)=95243; iiso=iiso+1
ResoCat(2)%idiso(iiso)=42095; iiso=iiso+1
ResoCat(2)%idiso(iiso)=43099; iiso=iiso+1
ResoCat(2)%idiso(iiso)=46108; iiso=iiso+1
ResoCat(2)%idiso(iiso)=54131; iiso=iiso+1
ResoCat(2)%idiso(iiso)=54135; iiso=iiso+1
ResoCat(2)%idiso(iiso)=55133

iiso=1
ResoCat(3)%idiso(iiso)=40000; iiso=iiso+1
ResoCat(3)%idiso(iiso)=40090; iiso=iiso+1
ResoCat(3)%idiso(iiso)=40091; iiso=iiso+1
ResoCat(3)%idiso(iiso)=40092; iiso=iiso+1
ResoCat(3)%idiso(iiso)=40094; iiso=iiso+1
ResoCat(3)%idiso(iiso)=40096

iiso=1
ResoCat(4)%idiso(iiso)=45103; iiso=iiso+1
ResoCat(4)%idiso(iiso)=47107; iiso=iiso+1
ResoCat(4)%idiso(iiso)=47109; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48000; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48106; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48108; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48110; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48111; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48112; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48113; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48114; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48116; iiso=iiso+1
ResoCat(4)%idiso(iiso)=49113; iiso=iiso+1
ResoCat(4)%idiso(iiso)=49115; iiso=iiso+1
ResoCat(4)%idiso(iiso)=60143; iiso=iiso+1
ResoCat(4)%idiso(iiso)=61147; iiso=iiso+1
ResoCat(4)%idiso(iiso)=61148; iiso=iiso+1
ResoCat(4)%idiso(iiso)=62149; iiso=iiso+1
ResoCat(4)%idiso(iiso)=62151; iiso=iiso+1
ResoCat(4)%idiso(iiso)=62152; iiso=iiso+1
ResoCat(4)%idiso(iiso)=63151; iiso=iiso+1
ResoCat(4)%idiso(iiso)=63152; iiso=iiso+1
ResoCat(4)%idiso(iiso)=63153; iiso=iiso+1
ResoCat(4)%idiso(iiso)=63154; iiso=iiso+1
ResoCat(4)%idiso(iiso)=63155; iiso=iiso+1
ResoCat(4)%idiso(iiso)=64155; iiso=iiso+1
ResoCat(4)%idiso(iiso)=64156; iiso=iiso+1
ResoCat(4)%idiso(iiso)=64157; iiso=iiso+1
ResoCat(4)%idiso(iiso)=64158; iiso=iiso+1
ResoCat(4)%idiso(iiso)=66160; iiso=iiso+1
ResoCat(4)%idiso(iiso)=66161; iiso=iiso+1
ResoCat(4)%idiso(iiso)=66162; iiso=iiso+1
ResoCat(4)%idiso(iiso)=66163; iiso=iiso+1
ResoCat(4)%idiso(iiso)=66164; iiso=iiso+1
ResoCat(4)%idiso(iiso)=68166; iiso=iiso+1
ResoCat(4)%idiso(iiso)=68167; iiso=iiso+1
ResoCat(4)%idiso(iiso)=72000; iiso=iiso+1
ResoCat(4)%idiso(iiso)=72174; iiso=iiso+1
ResoCat(4)%idiso(iiso)=72176; iiso=iiso+1
ResoCat(4)%idiso(iiso)=72177; iiso=iiso+1
ResoCat(4)%idiso(iiso)=72178; iiso=iiso+1
ResoCat(4)%idiso(iiso)=72179; iiso=iiso+1
ResoCat(4)%idiso(iiso)=72108

#else

nCat=5
allocate(ResoCat(nCat))
allocate(ResoCatUse(nCat,nmaxz))
ResoCatUse=.false.

ResoCat(1)%repid=92238
ResoCat(2)%repid=92235
ResoCat(3)%repid=40000
ResoCat(4)%repid=47107
ResoCat(5)%repid=94239

sum=0
niso=1; ResoCat(1)%niso=niso; sum=sum+niso
allocate(ResoCat(1)%idiso(niso))
niso=27; ResoCat(2)%niso=niso; sum=sum+niso
allocate(ResoCat(2)%idiso(niso))
niso=13; ResoCat(3)%niso=niso; sum=sum+niso
allocate(ResoCat(3)%idiso(niso))
niso=13; ResoCat(4)%niso=niso; sum=sum+niso
allocate(ResoCat(4)%idiso(niso))
!niso=nreshel-sum; ResoCat(5)%niso=niso
niso=26; ResoCat(5)%niso=niso
allocate(ResoCat(5)%idiso(niso))

iiso=1
ResoCat(1)%idiso(iiso)=92238
iiso=1
ResoCat(2)%idiso(iiso)=92235; iiso=iiso+1
ResoCat(2)%idiso(iiso)=90232; iiso=iiso+1
ResoCat(2)%idiso(iiso)=92233; iiso=iiso+1
ResoCat(2)%idiso(iiso)=92236; iiso=iiso+1
ResoCat(2)%idiso(iiso)=93237; iiso=iiso+1
ResoCat(2)%idiso(iiso)=40100; iiso=iiso+1
ResoCat(2)%idiso(iiso)=40190; iiso=iiso+1
ResoCat(2)%idiso(iiso)=40191; iiso=iiso+1
ResoCat(2)%idiso(iiso)=40192; iiso=iiso+1
ResoCat(2)%idiso(iiso)=40194; iiso=iiso+1
ResoCat(2)%idiso(iiso)=40196; iiso=iiso+1
ResoCat(2)%idiso(iiso)=42092; iiso=iiso+1
ResoCat(2)%idiso(iiso)=42094; iiso=iiso+1
ResoCat(2)%idiso(iiso)=42095; iiso=iiso+1
ResoCat(2)%idiso(iiso)=42096; iiso=iiso+1
ResoCat(2)%idiso(iiso)=42097; iiso=iiso+1
ResoCat(2)%idiso(iiso)=42098; iiso=iiso+1
ResoCat(2)%idiso(iiso)=42100; iiso=iiso+1
ResoCat(2)%idiso(iiso)=43099; iiso=iiso+1
ResoCat(2)%idiso(iiso)=45103; iiso=iiso+1
ResoCat(2)%idiso(iiso)=46108; iiso=iiso+1
ResoCat(2)%idiso(iiso)=54131; iiso=iiso+1
ResoCat(2)%idiso(iiso)=54135; iiso=iiso+1
ResoCat(2)%idiso(iiso)=55133; iiso=iiso+1
ResoCat(2)%idiso(iiso)=60143; iiso=iiso+1
ResoCat(2)%idiso(iiso)=61147; iiso=iiso+1
ResoCat(2)%idiso(iiso)=61148
iiso=1
ResoCat(3)%idiso(iiso)=40000; iiso=iiso+1
ResoCat(3)%idiso(iiso)=40090; iiso=iiso+1
ResoCat(3)%idiso(iiso)=40091; iiso=iiso+1
ResoCat(3)%idiso(iiso)=40092; iiso=iiso+1
ResoCat(3)%idiso(iiso)=40094; iiso=iiso+1
ResoCat(3)%idiso(iiso)=40096; iiso=iiso+1
ResoCat(3)%idiso(iiso)=72000; iiso=iiso+1
ResoCat(3)%idiso(iiso)=72174; iiso=iiso+1
ResoCat(3)%idiso(iiso)=72176; iiso=iiso+1
ResoCat(3)%idiso(iiso)=72177; iiso=iiso+1
ResoCat(3)%idiso(iiso)=72178; iiso=iiso+1
ResoCat(3)%idiso(iiso)=72179; iiso=iiso+1
ResoCat(3)%idiso(iiso)=72180
iiso=1
ResoCat(4)%idiso(iiso)=47107; iiso=iiso+1
ResoCat(4)%idiso(iiso)=47109; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48000; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48106; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48108; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48110; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48111; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48112; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48113; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48114; iiso=iiso+1
ResoCat(4)%idiso(iiso)=48116; iiso=iiso+1
ResoCat(4)%idiso(iiso)=49113; iiso=iiso+1
ResoCat(4)%idiso(iiso)=49115
iiso=1
ResoCat(5)%idiso(iiso)=94238; iiso=iiso+1
ResoCat(5)%idiso(iiso)=94239; iiso=iiso+1
ResoCat(5)%idiso(iiso)=94240; iiso=iiso+1
ResoCat(5)%idiso(iiso)=94241; iiso=iiso+1
ResoCat(5)%idiso(iiso)=94242; iiso=iiso+1
ResoCat(5)%idiso(iiso)=95241; iiso=iiso+1
ResoCat(5)%idiso(iiso)=95243; iiso=iiso+1
ResoCat(5)%idiso(iiso)=62149; iiso=iiso+1
ResoCat(5)%idiso(iiso)=62151; iiso=iiso+1
ResoCat(5)%idiso(iiso)=62152; iiso=iiso+1
ResoCat(5)%idiso(iiso)=63151; iiso=iiso+1
ResoCat(5)%idiso(iiso)=63152; iiso=iiso+1
ResoCat(5)%idiso(iiso)=63153; iiso=iiso+1
ResoCat(5)%idiso(iiso)=63154; iiso=iiso+1
ResoCat(5)%idiso(iiso)=63155; iiso=iiso+1
ResoCat(5)%idiso(iiso)=64155; iiso=iiso+1
ResoCat(5)%idiso(iiso)=64156; iiso=iiso+1
ResoCat(5)%idiso(iiso)=64157; iiso=iiso+1
ResoCat(5)%idiso(iiso)=64158; iiso=iiso+1
ResoCat(5)%idiso(iiso)=66160; iiso=iiso+1
ResoCat(5)%idiso(iiso)=66161; iiso=iiso+1
ResoCat(5)%idiso(iiso)=66162; iiso=iiso+1
ResoCat(5)%idiso(iiso)=66163; iiso=iiso+1
ResoCat(5)%idiso(iiso)=66164; iiso=iiso+1
ResoCat(5)%idiso(iiso)=68166; iiso=iiso+1
ResoCat(5)%idiso(iiso)=68167

#endif

do i = 1, nelthel
  lib => ldiso(i)
  lib%icat=0
  if (.not.lib%lreso) cycle
  nid = lib%nid
  xxx = mod(nid,1000)
  if (xxx.gt.500) then
    nid = nid - 500
  endif
  do ic = 1, nCat
    do iiso = 1, ResoCat(ic)%niso
      if (nid.eq.ResoCat(ic)%idiso(iiso)) lib%icat=ic
    enddo
  enddo
enddo         
end subroutine
end module