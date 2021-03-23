module files
character(256), save :: caseid,localfn,args
character(256) :: filename(40)
character(1024) :: workingdir
character(1024) :: outputdir
INTEGER, PARAMETER :: XsFileIdx = 3
INTEGER, PARAMETER :: SteamTblIdx= 4
INTEGER, PARAMETER :: InputFileIdx = 5
INTEGER, PARAMETER :: TempFileIdx = 6
INTEGER, PARAMETER :: DeplFileIdx = 7
INTEGER, PARAMETER :: XeEqFileIdx = 8
INTEGER, PARAMETER :: RstFileIdx = 9
INTEGER, PARAMETER :: DeplRstFileIdx = 10
INTEGER, PARAMETER :: ModTFileIdx = 11
INTEGER, PARAMETER :: FuelTFileIdx = 12
INTEGER, PARAMETER :: io_quick = 84
INTEGER, PARAMETER :: TempInclIdx = 13
INTEGER, PARAMETER :: rlbIdx = 14
INTEGER, PARAMETER :: slbIdx = 15
INTEGER, PARAMETER :: phlIdx = 16
INTEGER, PARAMETER :: pmatrxIdx = 18
INTEGER, PARAMETER :: pxsIdx = 17
!85 will use for MCP
INTEGER, save :: ibase,ixslib,ioutp,irstin,irstout,isum,icrho,IFbrho
INTEGER, save :: io1,io2,io3,io4,io5,io6,io7,io8,io9,io10
INTEGER, save :: io11,io12,io13,io14,io15,io16,io17,io18,io19,io20
INTEGER, save :: io21,io22,io23,io24,io25,io26,io27,io28,io29,io30
LOGICAL :: lIncfile=.FALSE.

!io4 : steamtable
!io5 : inputfile
!io6 : input file copy
!io8 : Output file
!io9 ~ io12 : VTK file
!io13 : isotope data out
!io14 : restart file
!io15 : Pin XS Printing & Cusping
!io16 : Effective Cross Section
!io17 : Binary Output File
!io18 : EFT something
!io19 : pin-wise power dist. converging histories
!io_quick : XS file read
data io1 ,io2 ,io3 ,io4 ,io5 ,io6 ,io7 ,io8 ,io9 ,io10,            &
     io11,io12,io13,io14,io15,io16,io17,io18,io19,io20,            &
     io21,io22,io23,io24,io25,io26,io27,io28,io29,io30             &
     / 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20, &
      21,22,23,24,25,26,27,28,29,30/
data caseid /'ntacer'/
data filename(5) /'ntracer.inp'/
END module
