SUBROUTINE QuadCoefGen(xs,ys,cs)
  IMPLICIT NONE
  REAL(8) :: xs(1:), ys(1:), cs(1:)
  REAL(8) :: det, a, b, c

  a = xs(3)-xs(2); b = xs(1) - xs(3); c = xs(2)-xs(1);
  det = 1./(a*xs(2)*xs(3)+b*xs(3)*xs(1)+c*xs(2)*xs(1))

  a = a*det*ys(1); cs(:) = (/xs(3)*xs(2), -xs(3)-xs(2), 1./)*a;
  b = b*det*ys(2); cs(:) = cs(:)+(/xs(3)*xs(1), -xs(3)-xs(1), 1./)*b;
  c = c*det*ys(3); cs(:) = cs(:)+(/xs(2)*xs(1), -xs(2)-xs(1), 1./)*c;
END SUBROUTINE QuadCoefGen

REAL(8) FUNCTION QuadFunc(cs,x)
  IMPLICIT NONE
  REAL(8) :: cs(:), x
  QuadFunc = cs(1)+cs(2)*x+cs(3)*x*x
END FUNCTION

SUBROUTINE PostCorrection(DeplFxrBundle, lCorrector, lfirst, Nths, pndcrit)
  USE HPDeplType
  USE OMP_LIB
  IMPLICIT NONE
  INTERFACE
    INTEGER FUNCTION IsoidSrch(IsoArr, NofIso, Idiso)
    USE HPDeplType
    INTEGER, INTENT(IN) :: NofIso, Idiso
    INTEGER, INTENT(IN) :: IsoArr(*)
    END FUNCTION
  END INTERFACE
  TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
  LOGICAL :: lCorrector, lfirst
  INTEGER :: Nths
  REAL(8), OPTIONAL :: pndcrit(:)

  INTEGER :: nTrueGd, nIsoGd
  INTEGER, POINTER :: IdTrueGd(:), IdIsoGd(:)
  TYPE(DeplFxr_Type), POINTER :: aFxr
  TYPE(GdFxrInfo_Type), POINTER :: GdFxrs(:)
  REAL(8), POINTER :: pnum(:)

  INTEGER, POINTER :: IdIsoEig(:)
  INTEGER :: i,j,k,l, iGd

  REAL(8) :: f, f1, f2, corND, a, b, c, d
  REAL(8),POINTER :: Gdpnum(:,:,:), f_pc(:)

  LOGICAL :: passcrit

  nIsoGd = DeplFxrBundle%nIsoGd; nTrueGd = DeplFxrBUndle%nTrueGd
  IdTrueGd => DeplFxrBundle%IdTrueGd; IdIsoGd => DeplFxrBundle%IdIsoGd

  GdFxrs => DeplFxrBundle%GdFxrBundle

  call OMP_SET_NUM_THREADS(Nths);
  !$OMP PARALLEL DO PRIVATE(Gdpnum, aFxr, pnum, j, iGd, passcrit)
  DO i = 1, nTrueGd
    aFxr => GdFxrs(i)%aFxr; Gdpnum => GdFxrs(i)%Gdpnum
    IF (lCorrector) THEN
      pnum=>aFxr%pnum_depl
    ELSE
      pnum=>aFxr%pnum_pre
    END IF
    IF (lCorrector) THEN
      !Gdpnum(2,-1,:) = Gdpnum(2,0,:)
      !Gdpnum(2, 0,:) = Gdpnum(2,1,:)
      Gdpnum(2,1,:) =  1.e-40;
      DO iGd = 1, nIsoGd
        passcrit = .FALSE.
        j = IsoIdSrch(aFxr%IdIsoEig, aFxr%NIsoEig, IdIsoGd(iGd))
        IF (PRESENT(pndcrit)) passcrit = (pnum(iGd).LT.pndcrit(j))
        IF (passcrit) CYCLE
        Gdpnum(2, 1,iGd) = pnum(IdIsoGd(iGd))
      END DO
      !Gdpnum(2, 0,:) = pnum(IdIsoGd(:))
    ELSE
      !Gdpnum(1,-1,:) = Gdpnum(1,0,:)
      !Gdpnum(1, 0,:) = Gdpnum(1,1,:)
      Gdpnum(1,1,:) =  1.e-40;
      DO iGd = 1, nIsoGd
        passcrit = .FALSE.
        j = IsoIdSrch(aFxr%IdIsoEig, aFxr%NIsoEig, IdIsoGd(iGd))
        IF (PRESENT(pndcrit)) passcrit = (pnum(iGd).LT.pndcrit(j))
        IF (passcrit) CYCLE
        Gdpnum(1, 1,iGd) = pnum(IdIsoGd(iGd))
      END DO
      IF (lfirst) THEN
        Gdpnum(2,0,:) = aFxr%pnum_depl(IdIsoGd(:))
      END IF
    END IF
  END DO
  !$OMP END PARALLEL DO
  IF (.NOT. lCorrector) THEN
    !$OMP PARALLEL DO PRIVATE(passcrit, Gdpnum, aFxr, f_pc, IdIsoEig, iGd, f, l, f1, f2, a, b, c, d, corND, pnum)
    DO i = 1, nTrueGd
      aFxr => GdFxrs(i)%aFxr;
      Gdpnum=>GdFxrs(i)%Gdpnum; f_pc=>GdFxrs(i)%f_pc
      IdIsoEig => aFxr%IdIsoEig
      DO l = 1, nIsoGd
        iGd = IsoidSrch(IdIsoEig, aFxr%NisoEig, IdIsoGd(l))
        IF (iGD .EQ. 0) CYCLE
        !f = Gdpnum(1,0,l)*Gdpnum(2,0,l)*Gdpnum(1,-1,l)*Gdpnum(2,-1,l)
        f = Gdpnum(1,0,l)*Gdpnum(2,0,l)*Gdpnum(2,-1,l)*Gdpnum(1,-1,l)
        IF (abs(f).LT.1.e-20) THEN
          f_pc(l) = 1.; CYCLE
        END IF
        f = Gdpnum(2,0,l)/Gdpnum(2,-1,l)
        IF (abs(f-1.).LT.1.e-20) THEN
          f_pc(l) = 1.; CYCLE
        END IF
        f1 = log(f)
        !IF (.NOT.(f.GE.0.5 .AND. f.LE.1.5)) CYCLE ! To clear NaN also
        f = Gdpnum(1,0,l)/Gdpnum(2,-1,l)
        IF (abs(f-1.).LT.1.e-20) THEN
          f_pc(l) = 1.; CYCLE
        END IF
        f2 = log(f)

        f = f1/f2

        f_pc(l) = f

        f = Gdpnum(1,1,l)*Gdpnum(2,0,l); IF(abs(f).LT.1.e-20) CYCLE
        f = Gdpnum(1,1,l)/Gdpnum(2,0,l); IF(abs(f-1.).LT.1.e-20) CYCLE
        a = Gdpnum(1,0,l); b = Gdpnum(1,1,l); c = Gdpnum(2,-1,l); d = Gdpnum(2,0,l);
        f1 = log(f)*f_pc(l)
        corND = Gdpnum(2,0,l)*exp(f1); IF(corND.LT.1.e-30) CYCLE

        passcrit = .FALSE.
        IF (PRESENT(pndcrit)) passcrit = (aFxr%pnum_sseig(iGd).LT.pndcrit(iGd))
        IF (passcrit) CYCLE

        f = corND/aFxr%pnum_sseig(iGd)
        IF (f.GT.1.5) f=1.5
        IF (f.LT.0.5) f=0.5
        corND = f*aFxr%pnum_sseig(iGd)
        IF (.NOT. (corND.GT.0 .OR. corND.LE.0)) THEN
          WRITE(*, '(A,I)') 'Nan during PostCorr of', l
          WRITE(*, '(A,ES12.5,A,ES12.5)') 'Pre. Nd1', a, 'Pre. Nd2', b
          WRITE(*, '(A,ES12.5,A,ES12.5)') 'Corr. Nd1', c, 'Corr. Nd2', d
          STOP
        END IF
        aFxr%pnum_sseig(iGd) = corND
        !!Gdpnum(1, 1, l) = corND
        !iGd = IdIsoGd(l)
        !aFxr%pnum_pre(iGd) = corND
      END DO
    END DO
    !$OMP END PARALLEL DO
  END IF
  !$OMP PARALLEL DO PRIVATE(Gdpnum, aFxr)
  DO i = 1, nTrueGd
    aFxr => GdFxrs(i)%aFxr; Gdpnum => GdFxrs(i)%Gdpnum
    IF (lCorrector) THEN
      Gdpnum(:,-1,:) = Gdpnum(:,0,:)
      Gdpnum(:, 0,:) = Gdpnum(:,1,:)
    END IF
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE

SUBROUTINE SetGdVars(DeplFxrBundle, lCorrector, Nths, pndcrit)
  USE HPDeplType
  USE OMP_LIB
  IMPLICIT NONE
  INTERFACE
    SUBROUTINE QuadCoefGen(xs,ys,cs)
    REAL(8) :: xs(1:), ys(1:), cs(1:)
    END SUBROUTINE QuadCoefGen

    INTEGER FUNCTION IsoidSrch(IsoArr, NofIso, Idiso)
    USE HPDeplType
    INTEGER, INTENT(IN) :: NofIso, Idiso
    INTEGER, INTENT(IN) :: IsoArr(*)
    END FUNCTION
  END INTERFACE
  TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
  LOGICAL :: lCorrector
  INTEGER :: Nths
  REAL(8), OPTIONAL :: pndcrit(:)

  INTEGER :: nTrueGd, nIsoGd
  INTEGER, POINTER :: IdTrueGd(:), IdIsoGd(:)
  TYPE(DeplFxr_Type), POINTER :: aFxr
  TYPE(GdFxrInfo_Type), POINTER :: GdFxrs(:)
  REAL(8), POINTER :: pnum(:)

  INTEGER :: i,j,k,l, iGd,i155

  REAL(8), POINTER :: GdRR(:,:), c_qd(:,:), Gd155(:)
  REAL(8) :: GdXs(-1:1)
  REAL(8) :: phi1g, ND155

  LOGICAL :: passcrit;

  IdTrueGd => DeplFxrBundle%IdTrueGd; IdIsoGd => DeplFxrBundle%IdIsoGd;
  GdFxrs => DeplFxrBundle%GdFxrBundle

  nTrueGd = DeplFxrBundle%nTrueGd; nIsoGd = DeplFxrBundle%nIsoGd; i155 = DeplFxrBundle%i155

  call OMP_SET_NUM_THREADS(Nths);

  IF (lCorrector) THEN
    !DEALLOCATE(GdXs);RETURN
    !$OMP PARALLEL DO PRIVATE(GdXs, aFxr, j, l, Gd155, GdRR, c_qd, pnum, ND155, phi1g, iGd)
    DO i = 1, nTrueGd
      aFxr => GdFxrs(i)%aFxr; Gd155 => GdFxrs(i)%Gd155
      GdRR => GdFxrs(i)%GdRR; c_qd => GdFxrs(i)%c_qd

      c_qd = 0.;
      pnum => aFxr%pnum_sseig
      j = IsoIdSrch(aFxr%IdIsoEig, aFxr%NIsoEig, i155)
      ND155 = pnum(j); phi1g = aFxr%NormFlux1g
      passcrit = .TRUE.
      IF (PRESENT(pndcrit)) passcrit = (ND155.GE.pndcrit(j))
      IF (passcrit.AND.(ND155.GE.1.E-20)) Gd155(1) = ND155
      !IF (ND155 .LT. 1.e-7) CYCLE
      DO l = 1, nIsoGd
        iGd = IdIsoGd(l); GdXs = 0.

        passcrit = .FALSE.
        j = IsoIdSrch(aFxr%IdIsoEig, aFxr%NIsoEig, iGd)
        IF (PRESENT(pndcrit)) passcrit = (pnum(iGd).LT.pndcrit(j))
        IF ((pnum(iGd) .LT. 1.e-20).OR.passcrit) CYCLE
        GdRR(1,l) = aFxr%xs1g(RctIdCAP,iGd)+aFxr%xs1g(RctIdNA,iGd)+aFxr%xs1g(RctIdNP,iGd)+aFxr%xs1g(RctIdCAPm,iGd)
        GdRR(1,l) = GdRR(1,l)*phi1g

        Gdxs(-1:1) = GdRR(-1:1,l)/phi1g

        IF (ANY(Gdxs(-1:1).LT.1.e-30)) CYCLE

        CALL QuadCoefGen(Gd155(:), Gdxs(:), c_qd(:,l))
      END DO
    END DO
    !$OMP END PARALLEL DO
  ELSE
    !i155 = IsoidSrch(IdIsoGd,nIsoGd,i155)
    !$OMP PARALLEL DO PRIVATE(GdXs, aFxr, j, l, Gd155, GdRR, c_qd, pnum, ND155, phi1g, iGd)
    DO i = 1, nTrueGd
      aFxr => GdFxrs(i)%aFxr; Gd155 => GdFxrs(i)%Gd155
      GdRR => GdFxrs(i)%GdRR; c_qd => GdFxrs(i)%c_qd

      GdRR(-1,:) = GdRR(0,:);  Gd155(-1) = Gd155(0);
      GdRR(0,:) = 0.0; Gd155(0) = 0.0;
      c_qd = 0.;
      pnum => aFxr%pnum_depl
      !ND155 = GdFxrs(i)%Gdpnum(1,1,i155);
      j = IsoIdSrch(aFxr%IdIsoEig, aFxr%NIsoEig, i155)
      ND155 = pnum(i155);
      IF (PRESENT(pndcrit)) passcrit = (ND155.GE.pndcrit(j))
      IF (passcrit.AND.(ND155.GE.1.E-20)) Gd155(0) = ND155
      phi1g = aFxr%NormFlux1g
      !Gd155(1) = ND155
      DO l = 1, nIsoGd
        iGd = IdIsoGd(l); GdXs = 0.

        passcrit = .FALSE.
        j = IsoIdSrch(aFxr%IdIsoEig, aFxr%NIsoEig, iGd)
        IF (PRESENT(pndcrit)) passcrit = (pnum(iGd).LT.pndcrit(j))
        IF ((pnum(iGd) .LT. 1.e-20).OR.passcrit) CYCLE
        GdRR(0,l) = aFxr%xs1g(RctIdCAP,iGd)+aFxr%xs1g(RctIdNA,iGd)+aFxr%xs1g(RctIdNP,iGd)+aFxr%xs1g(RctIdCAPm,iGd)
        GdRR(0,l) = GdRR(0,l)*phi1g

        !GdRR(1,l) = aFxr%xs1g(RctIdCAP,iGd)+aFxr%xs1g(RctIdNA,iGd)+aFxr%xs1g(RctIdNP,iGd)+aFxr%xs1g(RctIdCAPm,iGd)
        !GdRR(1,l) = GdRR(1,l)*phi1g

        !IF (GdRR(-1,l) .LT. 1.e-30) CYCLE
        !IF (GdRR(0,l) .LT. 1.e-30) CYCLE
        !Gdxs(-1:1) = GdRR(-1:1,l)/phi1g
        !
        !CALL QuadCoefGen(Gd155(:), Gdxs(:), c_qd(:,l))
      END DO
    END DO
    !$OMP END PARALLEL DO
  END IF
END SUBROUTINE
