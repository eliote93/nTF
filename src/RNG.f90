module rng
    integer(8), private     :: RN_MULT, RN_ADD, RN_SEED0, RN_BITS, RN_MASK, SEED, RN_STRIDE, RN_STRIDE2
    real(8), private        :: RN_NORM
    
contains    
subroutine long2short(long, short)
    integer(8), intent(in) :: long
    integer, intent(out) :: short(2)
    short(1)=rshift(long,32)
    short(2)=rshift(lshift(long,32),32)
end subroutine
subroutine short2long(short, long)
    integer, intent(in) :: short(2)
    integer(8), intent(out) :: long
    integer(8) :: temp
    temp=short(1)
    long=lshift(temp,32)
    long=long+short(2)
end subroutine

    subroutine initRN
      RN_MULT    = 5_8**19
      RN_ADD     = 0_8
      RN_STRIDE  = 152917_8
      RN_STRIDE2 = 315417_8
      RN_SEED0   = 5_8**19
      RN_BITS    = 48
      RN_MOD     = ishft( 1_8,       RN_BITS )
      RN_MASK    = ishft( not(0_8),  RN_BITS-64 )
      RN_NORM    = 1.0_8/2._8**48
      SEED       = RN_SEED0
    end subroutine
    
    function getRN() result(rn)
        real(8) :: rn
        integer(8) :: ts
        SEED  = iand( iand( RN_MULT*SEED, RN_MASK) + RN_ADD,  RN_MASK)
        rn    = SEED*RN_NORM
!        print *, "Global RN"
    end function

    function getRNS(tseed) result(rn)
        real(8) :: rn
        integer(8) :: tseed
        integer(8) :: ab
        integer :: a, b
!       rn=getRN()
        tseed  = iand( RN_MULT*tseed, RN_MASK)
        tseed  = iand(tseed + RN_ADD,  RN_MASK)
        rn  = tseed*RN_NORM
      
!        ! step 1
!        call long2short(tseed,a,b)
!        call short2long(a,b,ab)
        
    end function

!	subroutine changeseed(lseed)
!		integer(8) :: lseed
!		SEED = lseed
!	end subroutine
	
	function getSEED() 
	    integer(8) :: getSEED
	    getSEED = SEED
	end function
	
	subroutine setSEED(tseed)
	    integer(8) :: tseed
	    SEED = tseed
	end subroutine
	
	function shift_seed(tseed)
	    integer(8) :: tseed, shift_seed
	    shift_seed = tseed + RN_STRIDE-1
	end function
	
	function stride(tseed)
	    integer(8) :: tseed, stride
	    stride=tseed
!	    stride = iand(RN_MULT**RN_STRIDE*tseed + RN_ADD*(RN_MULT**RN_STRIDE-1)/(RN_MULT-1), RN_MASK)
	    stride = iand(RN_MULT**RN_STRIDE*tseed , RN_MASK)
	end function
	

	function strideN(tseed, N)
	    integer(8) :: tseed, strideN
	    integer :: N
	    integer(8) :: skip
	    skip=RN_STRIDE*N
	    strideN = iand(RN_MULT**skip*tseed , RN_MASK)
    endfunction	    

	subroutine changeseed(lseed)
		integer(8) :: lseed
		SEED = lseed
	end subroutine

end module
