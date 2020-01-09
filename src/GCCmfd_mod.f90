MODULE GCCmfd_mod
USE PARAM
USE TYPEDEF,         ONLY : GroupInfo_Type,       PinXS_TYPE,       CMFDLS_TYPE
IMPLICIT NONE
INTEGER, PRIVATE :: nglist(4)
REAL(4), PRIVATE :: GrpLBList(8,4)

TYPE(PinXS_Type), SAVE, POINTER :: GcPinXS(:, :)
TYPE(CMFDLS_Type), SAVE, POINTER :: GcCMFDLS(:)
INTEGER :: NGC

DATA nglist    /2, 4, 8, 1/
DATA GrpLBList / 6.2506E-01,         0._8,       0._8,         0._8,       0._8,        0._8,         0._8,  0._8,     &
                 9.1188E+3_8, 3.9279E+0_8, 6.2506E-1_8,        0._8,       0._8,        0._8,         0._8,  0._8,     &
                 2.2313E+6_8, 8.2085E+5_8, 9.1188E+3_8, 1.3007E+2_8, 3.9279E+0_8, 6.2506E-1_8, 1.4572E-1_8,  0._8,     &
                        0._8,        0._8,       0._8,         0._8,       0._8,        0._8,         0._8,  0._8/
INTERFACE

SUBROUTINE MakeGcPinXS(PinXS, GcPinXS, phi, ng, ngc, GroupInfo, GcGroupInfo)
USE PARAM
USE TYPEDEF,    ONLY : PinXS_Type,     GroupInfo_Type
IMPLICIT NONE
TYPE(PinXS_Type) :: PinXS, GcPinXS
TYPE(GroupInfo_Type) :: GroupInfo,     GcGroupInfo
REAL :: Phi(ng)
INTEGER :: ng, ngc

END SUBROUTINE

SUBROUTINE GenGcHomoXs(Core, PinXS, GcPinXS, phi, GcPhi, ng, ngc, GroupInfo, GcGroupInfo, PE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type,         PinXS_Type,         GroupInfo_Type,       &
                       PE_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PinXS_Type), POINTER :: PinXS(:, :), GcPinXS(:, :)
TYPE(GroupInfo_Type) :: GroupInfo,     GcGroupInfo
TYPE(PE_TYPE) :: PE
REAL, POINTER :: Phi(:, :, :), GcPhi(:, :, :)
INTEGER :: ng, ngc

REAL :: PHI0(ng)


END SUBROUTINE


SUBROUTINE GcRadCouplingCoeff(Core, PinXS, GcPinXS, Phi, ng, ngc, GroupInfo, PE)
USE PARAM
USE TYPEDEF,   ONLY : CoreInfo_Type,          PinXS_Type,          GroupInfo_Type,    &
                      Pin_Type,               PE_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PinXS_Type), POINTER :: PinXs(:, :), GcPinXs(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
REAL, POINTER :: phi(:, :, :)
INTEGER :: ng, ngc

END SUBROUTINE

SUBROUTINE GcAxCouplingCoeff(Core, GcPinXS, PHI, GcPhi, AxDtil, AxDhat, ng, ngc, GcGroupInfo, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,          AxFlx_Type,          PinXS_Type,       &
                           GroupInfo_Type,         PE_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PinXs_Type), POINTER :: GcPinXs(:, :)
TYPE(GroupInfo_Type) :: GcGroupInfo
TYPE(PE_Type) :: PE
REAL, POINTER :: PHI(:, :, :), GcPhi(:, :, :)
REAL, POINTER :: AxDtil(:, :, :, :), AxDhat(:, :, :, :)
INTEGER :: ng, ngc

END SUBROUTINE

SUBROUTINE SetGcCmfdLinearSystem(lcmfd, l3dim)
USE PARAM
IMPLICIT NONE
LOGICAL :: lcmfd, l3dim

END SUBROUTINE

SUBROUTINE GcCmfdPsiUpdt(PhiC, PsiC)
USE PARAM
IMPLICIT NONE
REAL, POINTER :: PhiC(:, :, :), PsiC(:, :)

END SUBROUTINE


SUBROUTINE GcCmfdSrcUpdt(SRC, psi, phi, eigv, ig)
USE PARAM
IMPLICIT NONE
REAL, POINTER :: SRC(:, :), psi(:, :), phi(:, :, :)
REAL :: Eigv
INTEGER :: ig


END SUBROUTINE

SUBROUTINE GcCmfdEigUpdt(psi, psid, eigv, psierr, PE)
USE PARAM
USE TYPEDEF,  ONLY : PE_TYPE
IMPLICIT NONE
TYPE(PE_TYPE), OPTIONAL :: PE
REAL, POINTER :: psi(:, :), psid(:, :)
REAL :: psierr, eigv
END SUBROUTINE

FUNCTION GcResidualError(phi, psi, eigv,  PE, constsrc)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE
IMPLICIT NONE
REAL :: GcResidualError
REAL, POINTER :: phi(:, :, :), psi(:, :)
REAL :: eigv
TYPE(PE_TYPE) :: PE
REAL, OPTIONAL :: ConstSrc

END FUNCTION

SUBROUTINE MgCmfdSolUpdt(CmInfo, GroupInfo, w)
USE PARAM
USE TYPEDEF,   ONLY : CmInfo_Type, GroupInfo_Type
IMPLICIT NONE
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
REAL :: w
END SUBROUTINE

END INTERFACE 

CONTAINS

SUBROUTINE UserdefineGrpStruct1(GroupInfo, nGC, ELBList)
USE XSLIB_MOD,  ONLY : enbhel
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: nGC
REAL :: ELBList(ngc-1)

INTEGER :: ng

INTEGER :: ig, ig0
INTEGER :: i, j
REAL :: EStruct(NGMAX),Eavg(NGMAX), temp, ELB
INTEGER :: EIdx(NGMAX)
ng = GroupInfo%ng
EStruct(1:ng) = enbhel(1:ng)
EStruct(ng+1) = 1.0E-4_8
DO i = 1, ng
  Eavg(i) = (EStruct(i) + EStruct(i+1))/2._8
ENDDO

DO i = 1, nGC-1
  ELB = ELBList(i)
  DO j = 1, ng-1
    temp = (ELB-Eavg(j))*(ELB-Eavg(j+1))
    IF(temp .LE. 0._8) THEN
      EIdx(i)=j
      EXIT
    ENDIF
  ENDDO
ENDDO

GroupInfo%GCStruct(1:2, 1) = (/1, EIdx(1)/)
DO ig = 2, nGC-1
  GroupInfo%GCStruct(1:2, ig) = (/EIdx(ig-1)+1, EIdx(ig)/)
ENDDO
GroupInfo%GCStruct(1:2, nGC) = (/EIdx(nGC-1)+1, ng/)
DO ig0 = 1, nGC
  DO ig = GroupInfo%GCStruct(1, ig0), GroupInfo%GCStruct(2, ig0)
    GroupInfo%InvGCStruct(ig) = ig0
  ENDDO
ENDDO
END SUBROUTINE


SUBROUTINE UserDefineGrpStruct2(GroupInfo, nGC, GroupLBIdx)
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: nGC
INTEGER :: GroupLBIdx(nGC-1)
INTEGER :: EIdx(ngmax)
INTEGER :: ig, ig0
INTEGER :: ng
ng = GroupInfo%ng

EIdx(1) = 1; EIdx(ngc+1) = ng
EIdx(2:ngc) = GroupLBIdx(1:ngc-1)
GroupInfo%GCStruct(1:2, 1) = (/1, EIdx(2)/)
DO ig = 2, nGC
  GroupInfo%GCStruct(1:2, ig) = (/EIdx(ig)+1, EIdx(ig+1)/)
ENDDO
!GroupInfo%GCStruct(1:2, nGC) = (/EIdx(nGC-1)+1, ng/) 
DO ig0 = 1, nGC
  DO ig = GroupInfo%GCStruct(1, ig0), GroupInfo%GCStruct(2, ig0)
    GroupInfo%InvGCStruct(ig) = ig0
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE PredefineGrpStruct(GroupInfo, nGC)
USE XSLIB_MOD,  ONLY : enbhel
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: nGC
INTEGER :: ng

INTEGER :: ig, ig0
INTEGER :: i, j, m
REAL :: EStruct(ngmax),Eavg(ngmax), temp, ELB
INTEGER :: EIdx(ngmax)
ng = GroupInfo%ng
EStruct(1:ng) = enbhel(1:ng)
EStruct(ng+1) = 1.0E-4_8
DO i = 1, ng
  Eavg(i) = (EStruct(i) + EStruct(i+1))/2._8
ENDDO
DO i = 1, 3
  IF(nGC .EQ. nGList(i)) EXIT
ENDDO
m = i

DO i = 1, nGC-1
  ELB = GrpLBList(i, m)
  DO j = 1, ng-1
    temp = (ELB-Eavg(j))*(ELB-Eavg(j+1))
    IF(temp .LE. 0._8) THEN
      EIdx(i)=j
      EXIT
    ENDIF
  ENDDO
ENDDO

IF( nGC .NE. 1 )THEN
    GroupInfo%GCStruct(1:2, 1) = (/1, EIdx(1)/)
    DO ig = 2, nGC-1
    GroupInfo%GCStruct(1:2, ig) = (/EIdx(ig-1)+1, EIdx(ig)/)
    ENDDO
    GroupInfo%GCStruct(1:2, nGC) = (/EIdx(nGC-1)+1, ng/)
ELSE
    GroupInfo%GCStruct(1:2, 1)=(/1, ng/)
ENDIF

DO ig0 = 1, nGC
    DO ig = GroupInfo%GCStruct(1, ig0), GroupInfo%GCStruct(2, ig0)
        GroupInfo%InvGCStruct(ig) = ig0
    ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SetGcCmfdEnv(CmInfo, ngc0)
USE PARAM
USE TYPEDEF,     ONLY : CmInfo_Type
IMPLICIT NONE
TYPE(CmInfo_Type) :: CmInfo
INTEGER :: ngc0

ngc = ngc0
GcPinXs => CmInfo%GcPinXs
GcCmfdLs => CmInfo%GcCmfdLs
END SUBROUTINE



END MODULE