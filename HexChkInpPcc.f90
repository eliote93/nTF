!SUBROUTINE HexChkInpPcc
!! ASSUME : defined value must be REAL
!! ASSUME : only two operations = 1. A == B / 2. defined(A)
!
!USE PARAM,      ONLY : TRUE, FALSE, POUND, ZERO
!USE files,      ONLY : io5, InputFIleIdx, InputFIleIdx
!USE inputcards, ONLY : oneline, probe
!USE ioutil,     ONLY : nfields
!
!IMPLICIT NONE
!
!INTEGER :: indev, nppc, ndef, idef
!CHARACTER*12 :: tmp, ASTRING
!
!CHARACTER*12, POINTER, DIMENSION(10) :: cdef ! Fixed
!REAL, POINTER, DIMENSION(10) :: rdef
!! ----------------------------------------------------
!
!! CHK : PPC exists
!indev = io5
!
!CALL OpenFile(io5, TRUE, FALSE, FALSE, filename(InputFIleIdx))
!
!DO
!  READ (indev, '(A512)', END = 1000) oneline
!  
!  IF (probe .NE. POUND) CYCLE
!  
!  GO TO 2000
!END DO
!
!1000 CONTINUE
!
!CLOSE (indev)
!
!RETURN
!
!2000 CONTINUE
!
!! FND : def. 1
!CALL OpenFile(io5, TRUE, FALSE, FALSE, filename(InputFIleIdx))
!
!ndef = 0
!
!DO
!  READ (indev, '(A512)', END = 3000) oneline
!  
!  IF (probe .NE. POUND) CYCLE
!  
!  READ (oneline, *) tmp
!  
!  IF (tmp .NE. '#define') CYCLE
!  
!  ndef = ndef + 1
!  
!  SELECT CASE (nfields(oneline))
!    CASE (2); READ (oneline, *) ASTRING, cdef(ndef)
!    CASE (3); READ (oneline, *) ASTRING, cdef(ndef), rdef(ndef)
!    CASE DEFAULT; CALL terminate("PPC ERR - DEF")
!  END SELECT
!END DO
!
!3000 CONTINUE
!
!REWIND (indev)
!
!! FND : def. 2
!DO
!  READ (indev, '(A512)', END = 4000) oneline
!  
!  IF (probe .NE. POUND) CYCLE
!  
!  READ (oneline, *) tmp
!  
!  IF (tmp .EQ. '#ifdef') THEN
!    READ (oneline, *) ASTRING, tmp
!    
!    IF (chkdef(tmp)) THEN
!      DO
!        READ (indev,   *) oneline
!        READ (oneline, *) tmp
!        
!        IF (tmp .EQ. '#define') THEN
!          ndef = ndef + 1
!          
!          SELECT CASE (nfields(oneline))
!            CASE (2); READ (oneline, *) ASTRING, cdef(ndef)
!            CASE (3); READ (oneline, *) ASTRING, cdef(ndef), rdef(ndef)
!            CASE DEFAULT; CALL terminate("PPC ERR - DEF")
!          END SELECT
!        ELSE IF (tmp .EQ. '#elif') THEN
!          EXIT
!        ELSE IF (tmp .EQ. '#endif') THEN
!          EXIT
!        END IF
!      END DO
!      
!      DO
!        
!      END DO
!    ELSE
!      DO
!        
!      END DO
!    END IF
!  END IF
!  
!  IF (tmp .EQ. '#if defined') THEN
!    
!  END IF
!  
!  ndef = ndef + 1
!  
!  SELECT CASE (nfields(oneline))
!    CASE (2); READ (oneline, *) ASTRING, cdef(ndef)
!    CASE (3); READ (oneline, *) ASTRING, cdef(ndef), rdef(ndef)
!    CASE DEFAULT; CALL terminate("PPC ERR - DEF")
!  END SELECT
!END DO
!
!3000 CONTINUE
!
!REWIND (indev)
!! ----------------------------------------------------
!CONTAINS
!
!FUNCTION chkdef(cc)
!
!IMPLICIT NONE
!
!CHARACTER*15 :: cc
!
!LOGICAL :: chkdef
!
!DO idef = 1, ndef
!  IF (cc .EQ. cdef(idef)) EXIT
!END DO
!
!chkdef = idef .LE. ndef
!
!END FUNCTION chkdef
!! ----------------------------------------------------
!
!END SUBROUTINE HexChkInpPcc