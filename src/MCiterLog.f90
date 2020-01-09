module MCIterLog
    USE files, ONLY : io8
    public
    private :: WriteLogCaption, WriteLog
    integer, private :: icyc, nsrc 
    real(8), private :: keff, kavg, kstd, shnMC, tcyc, tacc
    real(8), private :: ksum, ksqsum
    integer, private :: nacc
    integer, private :: istepcmfd
    real(8), private :: tcmfd, kcmfd, shncmfd
    logical, private :: lfdb
    logical, private :: lcmfd
    character *1025 logstr
contains

subroutine InitMCIterLog(lcmfd_)
    implicit none
    logical, intent(in) :: lcmfd_
    nacc=0
    tacc=0
    ksum=0
    ksqsum=0
    lcmfd=lcmfd_
    call WriteLogCaption()
end subroutine

subroutine WriteLogMC()
    implicit none
    if (lcmfd) then
        call WriteLogMCCMFD()
    else
        call WriteLogMCwoCMFD
    endif
end subroutine

subroutine ResetLogK()
    implicit none
    nacc=0
    ksum=0
    ksqsum=0
end subroutine

subroutine ReSetLogTime()
    implicit none
    tacc=0
end subroutine

subroutine SetLogCycle(icyc_)
    implicit none
    integer :: icyc_
    icyc=icyc_
end subroutine

subroutine SetLogNSSRC(nsrc_)
    implicit none
    integer :: nsrc_
    nsrc=nsrc_
end subroutine

subroutine SetLogTCYC(tcyc_)
    implicit none
    real(8) :: tcyc_
    tcyc=tcyc_
    tacc=tacc+tcyc
end subroutine

subroutine SetLogTCMFD(tcmfd_)
    implicit none
    real(8) :: tcmfd_
    tcmfd=tcmfd_
end subroutine

subroutine SetLogKeff(keff_)
    implicit none
    real(8) :: kefF_
    nacc=nacc+1
    keff=keff_
    ksum=ksum+keff
    ksqsum=ksqsum+keff**2
    kavg=ksum/nacc
    if (nacc .ne. 1) then
        kstd=sqrt((ksqsum-ksum**2/nacc)/nacc/(nacc-1))
    else
        kstd=0.
    endif
end subroutine

subroutine SetLogShn(shn)
    implicit none
    real(8), intent(in) :: shn
    shnMC=shn
end subroutine

subroutine SetLogCMFD(istep, fdb, time, k, shn)
    implicit none
    integer, intent(in) :: istep
    logical, intent(in) :: fdb
    real(8), intent(in) :: time, k, shn

    istepcmfd=istep
    lfdb=fdb
    tcmfd=time
    kcmfd=k
    shncmfd=shn
end subroutine

subroutine WriteLogCaption()
    implicit none
    if (lcmfd) then
        write(logstr, '(a)') " CYCLE      NHT    T[sec]   TA[sec]    k      AVG_k    STDEV_k    SHN_F  CMFD [STEP/FDB/Time/k-eff/SHN]"
    else                                                                                      
        write(logstr, '(a)') " CYCLE      NHT    T[sec]   TA[sec]    k      AVG_k    STDEV_k    SHN_F"
    endif
    call WriteLog()
end subroutine

SUBROUTINE OutputEditMC_2
  IMPLICIT NONE
END SUBROUTINE

subroutine WriteLogMCwoCMFD()
    implicit none

101 format(i5, 1x,i10, 1x,f8.3, 1x,f8.3, 1x,f8.5, 1x,f9.6, 1x,f9.6, 1x, f8.4)    
    write(logstr, 101), icyc, nsrc, tcyc, tacc, keff, kavg, kstd, shnMC
    CALL OutputEditMC_2()
    call WriteLog()
end subroutine
subroutine WriteLogMCCMFD()
    implicit none

201 format(1x,i5, 1x,i10, 1x,f8.3, 1x,f8.3, 1x,f8.5, 1x,f9.6, 1x,f9.6, 1x, f8.4, 2x,i5,1x, l, 1x, f8.3, 1x, f8.5,1x, f8.4,1x)    
    write(logstr, 201), icyc, nsrc, tcyc, tacc, keff, kavg, kstd, shnMC, istepcmfd, lfdb, tcmfd, kcmfd, shncmfd
    call WriteLog()
end subroutine

subroutine WriteLog()
    implicit none
    write(*,'(a)') trim(logstr)
    WRITE(io8, '(a)') trim(logstr)
end subroutine
    
end module
