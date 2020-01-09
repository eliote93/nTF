MODULE SP3DEF
USE TYPEDEF, ONLY : PinXS_Type, AxFLX_TYPE
!TYPE PinXS_Type
!  SEQUENCE
!!  REAL,POINTER :: TRSPT(:),AB(:),RMV(:),NUFS(:),CHI(:),SC(:,:)
!!  REAL,POINTER :: D0(:),D2(:)
!  TYPE(SCATMAT),POINTER,DIMENSION(:) :: XSS
!  REAL,POINTER,DIMENSION(:,:) :: DTIL, DHAT         !Radial Dhat and Dtil
!  REAL,POINTER,DIMENSION(:,:) :: AxDtil, AxDhat     !
!  REAL,POINTER,DIMENSION(:) :: XSD, XSD2, XST, XSTR, XSR, XSNF, XSKF, CHI !,PHI,FAC,XSRD
!  REAL, POINTER, DIMENSION(:) :: PHI, FAC, XSRD
!!REAL :: BETAT,OMEGA,BETAP        
!  !INTEGER :: NG
!END TYPE

TYPE SCATMAT
  SEQUENCE
  INTEGER IB,IE
  REAL,POINTER :: FROM(:), DUMMY
  REAL :: WITHINGROUPSCAT
END TYPE

!TYPE AXFLX_TYPE
!  SEQUENCE
!  REAL,POINTER :: PHI(:,:,:) !FLUX, Modlar Space Flux
!  REAL,POINTER :: PSI(:)  !Fission Source
!  REAL,POINTER :: TLKG(:,:,:) !Source
!  REAL,POINTER :: JOUT(:,:,:), DHAT(:,:), DTIL(:, :), PDHAT(:, :) !CURRENT
!  REAL,POINTER :: S(:,:), KSQ(:,:), QT(:,:) !EIGENVECTOR,KAFFA SQUARE,QHAT TRANS
!END TYPE

TYPE MODALFLX_TYPE
  SEQUENCE
  REAL,POINTER :: PHI(:,:,:) !FLUX, Modlar Space Flux
  REAL,POINTER :: QH(:,:) !Source, MODAL SPACE SOURCE
  REAL,POINTER :: A(:,:),B(:,:),PSOL(:,:,:)
END TYPE

TYPE AxGEOM_TYPE
  SEQUENCE
!  INTEGER :: NG
!  INTEGER :: NMESH,NCOMP
!  INTEGER :: BC(2)
!  REAL,POINTER :: H(:),HINV(:) 
!  INTEGER,POINTER :: COMP(:)  
  
  INTEGER :: IX, IY, NG, idum
  INTEGER :: myzb, myze, myzbf, myzef
  INTEGER :: NMESH, NCOMP
  INTEGER :: BC(2)
  REAL :: Area
  REAL, POINTER :: H(:), HINV(:)
  INTEGER, POINTER :: COMP(:), idum2    
  LOGICAL :: lTransient = .FALSE.
END TYPE     

TYPE AxSolItrCntl_TYPE
  SEQUENCE
  INTEGER :: GroupItr = 10
  INTEGER :: NodalItr = 5
  INTEGER :: SourceItr = 1
  
!  INTEGER :: GroupItr = 5
!  INTEGER :: NodalItr = 10 
!  INTEGER :: SourceItr = 1  
END TYPE

REAL, PARAMETER :: MaxK = 400

TYPE SENMMAT_TYPE
  SEQUENCE
  !REAL,POINTER ::DIAG(:, :), OFFDIAG(:, :, :)
  !INTEGER, POINTER :: nOffDiag(:), OffDiagLoc(:, :)
  !INTEGER, POINTER :: L(:, :, :), U(:, :, :)
  INTEGER, POINTER :: nElmt(:), ElmtLoc(:, :), DiagIdx(:)
  REAL, POINTER :: Elmt(:, :, :)
  INTEGER, POINTER :: LowerElmtIdx(:, :), nLowerElmt(:)
  REAL, POINTER :: LU(:, :, :) !, U(:, :, :)
ENDTYPE

END MODULE