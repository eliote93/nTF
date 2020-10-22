module allocs

    INTEGER(8) :: nbytesf = 0
    INTEGER(8) :: memprofile(1000)
    interface dmalloc
    module procedure mallocf1
    module procedure mallocd1
    module procedure malloci1
    module procedure mallocl1
    module procedure mallocf2
    module procedure mallocd2
    module procedure malloci2
    module procedure mallocl2
    module procedure mallocf3
    module procedure mallocd3
    module procedure malloci3
    module procedure mallocl3
    module procedure mallocf4
    module procedure mallocd4
    module procedure malloci4
    module procedure mallocf5
    module procedure mallocd5
    module procedure malloci5
    END interface
    interface dmalloc0
    module procedure mallocf01
    module procedure mallocd01
    module procedure malloci01
    module procedure mallocl01
    module procedure mallocf02
    module procedure mallocd02
    module procedure malloci02
    module procedure mallocl02
    module procedure mallocf03
    module procedure mallocd03
    module procedure malloci03
    module procedure mallocl03
    module procedure mallocf04
    module procedure mallocd04
    module procedure malloci04
    module procedure mallocl04
    module procedure mallocf05
    module procedure mallocd05
    module procedure malloci05
    module procedure mallocl05
    END interface
    !
    contains

    SUBROUTINE mallocf1(a,n1)

    REAL(4),pointer :: a(:)
    allocate(a(n1))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocd1(a,n1)

    REAL(8),pointer :: a(:)
    allocate(a(n1))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE malloci1(a,n1)

    INTEGER,pointer :: a(:)
    allocate(a(n1))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocl1(a,n1)

    logical,pointer :: a(:)
    allocate(a(n1))
    a=.FALSE.
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocf2(a,n1,n2)

    REAL(4),pointer :: a(:,:)
    allocate(a(n1,n2))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocd2(a,n1,n2)

    REAL(8),pointer :: a(:,:)
    allocate(a(n1,n2))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE malloci2(a,n1,n2)

    INTEGER,pointer :: a(:,:)
    allocate(a(n1,n2))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocl2(a,n1,n2)

    logical,pointer :: a(:,:)
    allocate(a(n1,n2))
    a=.FALSE.
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocf3(a,n1,n2,n3)

    REAL(4),pointer :: a(:,:,:)
    allocate(a(n1,n2,n3))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocd3(a,n1,n2,n3)

    REAL(8),pointer :: a(:,:,:)
    allocate(a(n1,n2,n3))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE malloci3(a,n1,n2,n3)

    INTEGER,pointer :: a(:,:,:)
    allocate(a(n1,n2,n3))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocl3(a,n1,n2,n3)

    LOGICAL,pointer :: a(:,:,:)
    allocate(a(n1,n2,n3))
    a=.FALSE.
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocf4(a,n1,n2,n3,n4)

    REAL(4),pointer :: a(:,:,:,:)
    allocate(a(n1,n2,n3,n4))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocd4(a,n1,n2,n3,n4)

    REAL(8),pointer :: a(:,:,:,:)
    allocate(a(n1,n2,n3,n4))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE malloci4(a,n1,n2,n3,n4)

    INTEGER,pointer :: a(:,:,:,:)
    allocate(a(n1,n2,n3,n4))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE

    SUBROUTINE mallocf5(a,n1,n2,n3,n4,n5)

    REAL(4),pointer :: a(:,:,:,:,:)
    allocate(a(n1,n2,n3,n4,n5))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocd5(a,n1,n2,n3,n4,n5)

    REAL(8),pointer :: a(:,:,:,:,:)
    allocate(a(n1,n2,n3,n4,n5))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE malloci5(a,n1,n2,n3,n4,n5)

    INTEGER,pointer :: a(:,:,:,:,:)
    allocate(a(n1,n2,n3,n4,n5))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE

    !
    SUBROUTINE mallocf01(a,nb1,ne1)

    REAL(4),pointer :: a(:)
    allocate(a(nb1:ne1))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocd01(a,nb1,ne1)

    REAL(8),pointer :: a(:)
    allocate(a(nb1:ne1))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE malloci01(a,nb1,ne1)

    INTEGER,pointer :: a(:)
    allocate(a(nb1:ne1))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocl01(a,nb1,ne1)

    LOGICAL,pointer :: a(:)
    allocate(a(nb1:ne1))
    a=.FALSE.
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocf02(a,nb1,ne1,nb2,ne2)

    REAL(4),pointer :: a(:,:)
    allocate(a(nb1:ne1,nb2:ne2))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocd02(a,nb1,ne1,nb2,ne2)

    REAL(8),pointer :: a(:,:)
    allocate(a(nb1:ne1,nb2:ne2))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE malloci02(a,nb1,ne1,nb2,ne2)

    INTEGER,pointer :: a(:,:)
    allocate(a(nb1:ne1,nb2:ne2))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocl02(a,nb1,ne1,nb2,ne2)

    LOGICAL,pointer :: a(:,:)
    allocate(a(nb1:ne1,nb2:ne2))
    a=.FALSE.
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocf03(a,nb1,ne1,nb2,ne2,nb3,ne3)

    REAL(4),pointer :: a(:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocd03(a,nb1,ne1,nb2,ne2,nb3,ne3)

    REAL(8),pointer :: a(:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE malloci03(a,nb1,ne1,nb2,ne2,nb3,ne3)

    INTEGER,pointer :: a(:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocl03(a,nb1,ne1,nb2,ne2,nb3,ne3)

    LOGICAL,pointer :: a(:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3))
    a=.FALSE.
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocf04(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4)

    REAL(4),pointer :: a(:,:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocd04(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4)

    REAL(8),pointer :: a(:,:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE malloci04(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4)

    INTEGER,pointer :: a(:,:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocl04(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4)

    LOGICAL,pointer :: a(:,:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4))
    a=.FALSE.
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocf05(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5)

    REAL(4),pointer :: a(:,:,:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocd05(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5)

    REAL(8),pointer :: a(:,:,:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE malloci05(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5)

    INTEGER,pointer :: a(:,:,:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5))
    a=0
    nbytesf=nbytesf+size(a)
    END SUBROUTINE
    !
    SUBROUTINE mallocl05(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5)

    LOGICAL,pointer :: a(:,:,:,:,:)
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5))
    a=.FALSE.
    nbytesf=nbytesf+size(a)
    END SUBROUTINE

END module