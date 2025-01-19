SUBROUTINE DPNM(X,N,P)
    IMPLICIT REAL*8(D-H,O-Z)
    DIMENSION P(*)
    STH=DSQRT(1.D00-X*X)
    PI4=16.D00*DATAN(1.D00)
    P(1)=1.D00/DSQRT(PI4)
    IF(N.LE.0) RETURN
    SOU=1.D00
    ZN=-1.
    NP1=N+1
    IFLAG=0
    DO 3 MS=1,NP1
    M=MS-1
    F1=DFLOAT(M+M)
    I=(M+1)*(M+2)/2
    IF(M.EQ.0) GOTO 1
    SOU=SOU*(F1+1.D00)/F1
    IF(IFLAG-1.LE.0) THEN
    STHM=STH**M
    IF(STHM.LT.1.D-55) IFLAG=1
    GOTO 25
    ELSE
    STHM=0.D00
    ENDIF
 25 CONTINUE
    P(I)=ZN*DSQRT(SOU/PI4)*STHM
    IF(M.EQ.N) GOTO 3
    ZN=-ZN
  1 I1=I+M+1
    P(I1)=DSQRT(F1+3.D00)*X*P(I)
    MP2=M+2
    IF(MP2.GT.N) GOTO 3
    DO 2 J=MP2,N
    F1=DFLOAT((J-M)*(J+M))
    F2=DFLOAT((J+J-1)*(J+J+1))
    F3=DFLOAT((J+J+1)*(J-M-1)*(J+M-1))/DFLOAT(J+J-3)
    F2=DSQRT(F2/F1)*X
    F3=DSQRT(F3/F1)
    I=J*(J+1)/2+M+1
    I1=I-J
    I2=I1-J+1
    P(I)=F2*P(I1)-F3*P(I2)
  2 CONTINUE
  3 CONTINUE
    RETURN
    END

SUBROUTINE bandec(a,n,m1,m2,np,mp,al,mpl,indx,d)
    INTEGER m1,m2,mp,mpl,n,np,indx(*)
    REAL*8 d,a(np,mp),al(np,mpl),TINY
    PARAMETER (TINY=1.d-20)
    INTEGER i,j,k,l,mm
    REAL*8 dum
    mm=m1+m2+1
    if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in bandec'
    l=m1
    do 13 i=1,m1
        do 11 j=m1+2-i,mm
        a(i,j-l)=a(i,j)
    11      continue
        l=l-1
        do 12 j=mm-l,mm
        a(i,j)=0d0
    12      continue
    13    continue
    d=1d0
    l=m1
    do 18 k=1,n
        dum=a(k,1)
        i=k
        if(l.lt.n)l=l+1
        do 14 j=k+1,l
        if(abs(a(j,1)).gt.abs(dum))then
            dum=a(j,1)
            i=j
        endif
    14      continue
        indx(k)=i
        if(dum.eq.0d0) a(k,1)=TINY
        if(i.ne.k)then
        d=-d
        do 15 j=1,mm
            dum=a(k,j)
            a(k,j)=a(i,j)
            a(i,j)=dum
    15        continue
        endif
        do 17 i=k+1,l
        dum=a(i,1)/a(k,1)
        al(k,i-k)=dum
        do 16 j=2,mm
            a(i,j-1)=a(i,j)-dum*a(k,j)
    16        continue
        a(i,mm)=0d0
    17      continue
    18    continue
    return
    END subroutine bandec
        
SUBROUTINE banbks(a,n,m1,m2,np,mp,al,mpl,indx,b)
    INTEGER m1,m2,mp,mpl,n,np,indx(*)
    REAL*8 a(np,mp),al(np,mpl),b(*)
    INTEGER i,k,l,mm
    REAL*8 dum
    mm=m1+m2+1
    if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in banbks'
    l=m1
    do 12 k=1,n
        i=indx(k)
        if(i.ne.k)then
        dum=b(k)
        b(k)=b(i)
        b(i)=dum
        endif
        if(l.lt.n)l=l+1
        do 11 i=k+1,l
        b(i)=b(i)-al(k,i-k)*b(k)
    11      continue
    12    continue
    l=1
    do 14 i=n,1,-1
        dum=b(i)
        do 13 k=2,l
        dum=dum-a(i,k)*b(k+i-1)
    13      continue
        b(i)=dum/a(i,1)
        if(l.lt.mm) l=l+1
    14    continue
    return
    END subroutine banbks

SUBROUTINE test_equation_of_continuity(number_of_layers, jmax, j, m, radius, delta_r, displacement)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: radius, delta_r
    COMPLEX(KIND=8), INTENT(IN) :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)

    ! Declare local variables
    REAL(KIND=8) :: j1, factor_radius
    REAL(KIND=8) :: a1, a2, b1, b2
    COMPLEX(KIND=16) :: error, highest_error
    INTEGER :: i, index_highest_error
    INTEGER :: disp_idx

    ! Precompute degree-based constants
    j1 = DBLE(j)  ! Double precision conversion of j
    a1 = SQRT(j1 / (2.0D0 * j1 + 1.0D0))
    a2 = - (j1 - 1.0D0) * a1
    b1 = - SQRT((j1 + 1.0D0) / (2.0D0 * j1 + 1.0D0))
    b2 = (j1 + 2.0D0) * b1

    ! Initialize values
    highest_error = CMPLX(0.0D0, 0.0D0, KIND=16)  ! Complex zero
    index_highest_error = 0

    ! Loop through layers
    DO i = 1, number_of_layers
        ! Precompute terms reused in calculations
        factor_radius = radius + (i - 1) * delta_r
        disp_idx = 3 * (j * (j + 1) / 2 + m)

        ! Compute error terms
        error = ((a2 / (2.0D0 * factor_radius)) - (a1 / delta_r)) * displacement(i, disp_idx - 1)
        error = error + ((b2 / (2.0D0 * factor_radius)) - (b1 / delta_r)) * displacement(i, disp_idx + 1)
        error = error + ((a2 / (2.0D0 * factor_radius)) + (a1 / delta_r)) * displacement(i + 1, disp_idx - 1)
        error = error + ((b2 / (2.0D0 * factor_radius)) + (b1 / delta_r)) * displacement(i + 1, disp_idx + 1)

        ! Update largest error if this error is greater
        IF (ABS(error) > ABS(highest_error)) THEN
            highest_error = error
            index_highest_error = i
        END IF
    END DO

    ! Print results
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the equation of continuity for degree j =", j, " and order m =", m
    PRINT '(A, I4)', "Largest error is at layer:", index_highest_error
    PRINT '(A, 2ES25.16)', "Value of the largest error:", highest_error

END SUBROUTINE test_equation_of_continuity

SUBROUTINE test_equation_of_motion(number_of_layers, jmax, j, m, radius, delta_r, cauchy, cauchy_isotropic, volume_force)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: radius, delta_r
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_isotropic(number_of_layers, jmax + 1, jmax + 1)
    COMPLEX(KIND=8), INTENT(IN) :: volume_force(number_of_layers - 1, 2)

    ! Declare local variables
    REAL(KIND=8) :: j1, factor_radius, factor_delta_r
    REAL(KIND=8) :: c1, c2, d1, d2, e1, e2, f1, f2, g1, g2, h1, h2
    COMPLEX(KIND=8) :: error, highest_error_1, highest_error_2
    INTEGER :: i, index_highest_error_1, index_highest_error_2
    INTEGER :: cauchy_idx

    ! Precompute constants
    j1 = DBLE(j)
    factor_delta_r = 1.0D0 / delta_r

    c1 = -SQRT(j1 / (3.0D0 * (2.0D0 * j1 + 1.0D0)))
    c2 = (j1 + 1.0D0) * c1
    d1 = SQRT((j1 - 1.0D0) / (2.0D0 * j1 - 1.0D0))
    d2 = -(j1 - 2.0D0) * d1
    e1 = -SQRT((j1 + 1.0D0) * (2.0D0 * j1 + 3.0D0) / (6.0D0 * (2.0D0 * j1 - 1.0D0) * (2.0D0 * j1 + 1.0D0)))
    e2 = (j1 + 1.0D0) * e1
    f1 = SQRT((j1 + 1.0D0) / (3.0D0 * (2.0D0 * j1 + 1.0D0)))
    f2 = -j1 * f1
    g1 = SQRT(j1 * (2.0D0 * j1 - 1.0D0) / (6.0D0 * (2.0D0 * j1 + 1.0D0) * (2.0D0 * j1 + 3.0D0)))
    g2 = -j1 * g1
    h1 = -SQRT((j1 + 2.0D0) / (2.0D0 * j1 + 3.0D0))
    h2 = (j1 + 3.0D0) * h1

    ! Initialize error tracking
    highest_error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
    highest_error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)
    index_highest_error_1 = 0
    index_highest_error_2 = 0


    ! Loop through layers for both equations
    DO i = 1, number_of_layers - 1
        factor_radius = 1.0D0 / (2.0D0 * (radius + (i - 0.5D0) * delta_r))
        cauchy_idx = 5 * (j * (j + 1) / 2 + m)

        ! First equation test
        error = (c2 * factor_radius - c1 * factor_delta_r) * cauchy_isotropic(i, j + 1, m + 1) + &
                (d2 * factor_radius - d1 * factor_delta_r) * cauchy(i, cauchy_idx - 7) + &
                (e2 * factor_radius - e1 * factor_delta_r) * cauchy(i, cauchy_idx - 5) + &
                (c2 * factor_radius + c1 * factor_delta_r) * cauchy_isotropic(i + 1, j + 1, m + 1) + &
                (d2 * factor_radius + d1 * factor_delta_r) * cauchy(i + 1, cauchy_idx - 7) + &
                (e2 * factor_radius + e1 * factor_delta_r) * cauchy(i + 1, cauchy_idx - 5)

        IF (j == 2) THEN
            SELECT CASE (m)
                CASE (0)
                    error = error + volume_force(i, 1)
                CASE (2)
                    error = error + volume_force(i, 2)
            END SELECT
        END IF

        IF (ABS(error) > ABS(highest_error_1)) THEN
            highest_error_1 = error
            index_highest_error_1 = i
        END IF

        ! Second equation test
        error = (f2 * factor_radius - f1 * factor_delta_r) * cauchy_isotropic(i, j + 1, m + 1) + &
                (g2 * factor_radius - g1 * factor_delta_r) * cauchy(i, cauchy_idx - 5) + &
                (h2 * factor_radius - h1 * factor_delta_r) * cauchy(i, cauchy_idx - 3) + &
                (f2 * factor_radius + f1 * factor_delta_r) * cauchy_isotropic(i + 1, j + 1, m + 1) + &
                (g2 * factor_radius + g1 * factor_delta_r) * cauchy(i + 1, cauchy_idx - 5) + &
                (h2 * factor_radius + h1 * factor_delta_r) * cauchy(i + 1, cauchy_idx - 3)

        IF (ABS(error) > ABS(highest_error_2)) THEN
            highest_error_2 = error
            index_highest_error_2 = i
        END IF
    END DO

    ! Print results
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the equation of motion for degree j =", j, " and order m =", m
    PRINT '(A, I4)', "Largest error in first equation is at layer:", index_highest_error_1
    PRINT '(A, 2ES25.16)', "Value of the largest error in first equation:", highest_error_1
    PRINT '(A, I4)', "Largest error in second equation is at layer:", index_highest_error_2
    PRINT '(A, 2ES25.16)', "Value of the largest error in second equation:", highest_error_2

END SUBROUTINE

SUBROUTINE test_rheology(number_of_layers, jmax, j, m, radius, delta_r, mu, cauchy, cauchy_integral, displacement)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: radius, delta_r, mu
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_integral(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)

    ! Declare local variables
    REAL(KIND=8) :: j1, factor_radius, factor_delta_r
    REAL(KIND=8) :: p1, p2, q1, q2, r1, r2, s1, s2
    COMPLEX(KIND=8) :: error, highest_error_1, highest_error_2, highest_error_3
    INTEGER :: i, index_highest_error_1, index_highest_error_2, index_highest_error_3
    INTEGER :: cauchy_idx, disp_idx

    ! Precompute constants
    j1 = DBLE(j)
    factor_delta_r = 1.0D0 / delta_r

    p1 = -SQRT((j1 - 1.0D0) / (2.0D0 * j1 - 1.0D0))
    p2 = j1 * p1
    q1 = SQRT((j1 + 1.0D0) * (2.0D0 * j1 + 3.0D0) / (6.0D0 * (2.0D0 * j1 - 1.0D0) * (2.0D0 * j1 + 1.0D0)))
    q2 = -(j1 - 1.0D0) * q1
    r1 = -SQRT((j1 * (2.0D0 * j1 - 1.0D0)) / (6.0D0 * (2.0D0 * j1 + 1.0D0) * (2.0D0 * j1 + 3.0D0)))
    r2 = (j1 + 2.0D0) * r1
    s1 = SQRT((j1 + 2.0D0) / (2.0D0 * j1 + 3.0D0))
    s2 = -(j1 + 1.0D0) * s1

    ! Initialize error tracking
    highest_error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
    highest_error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)
    highest_error_3 = CMPLX(0.0D0, 0.0D0, KIND=8)
    index_highest_error_1 = 0
    index_highest_error_2 = 0
    index_highest_error_3 = 0

    cauchy_idx = 5 * (j * (j + 1) / 2 + m)
    disp_idx = 3 * (j * (j + 1) / 2 + m)

    ! Precompute factors and process each equation test
    DO i = 1, number_of_layers
        factor_radius = 1.0D0 / (2.0D0 * (radius + (i - 1.0D0) * delta_r))

        ! ---------------------------------
        ! First equation test
        ! ---------------------------------
        error = cauchy(i, cauchy_idx - 5) / (2.0D0 * mu) &
              + (q2 * factor_radius - q1 * factor_delta_r) * displacement(i, disp_idx - 1) &
              + (r2 * factor_radius - r1 * factor_delta_r) * displacement(i, disp_idx + 1) &
              + (q2 * factor_radius + q1 * factor_delta_r) * displacement(i + 1, disp_idx - 1) &
              + (r2 * factor_radius + r1 * factor_delta_r) * displacement(i + 1, disp_idx + 1) &
              - cauchy_integral(i, cauchy_idx - 5)

        IF (ABS(error) > ABS(highest_error_1)) THEN
            highest_error_1 = error
            index_highest_error_1 = i
        END IF

        ! ---------------------------------
        ! Second equation test
        ! ---------------------------------
        error = cauchy(i, cauchy_idx - 7) / (2.0D0 * mu) &
              + (p2 * factor_radius - p1 * factor_delta_r) * displacement(i, disp_idx - 1) &
              + (p2 * factor_radius + p1 * factor_delta_r) * displacement(i + 1, disp_idx - 1) &
              - cauchy_integral(i, cauchy_idx - 7)

        IF (ABS(error) > ABS(highest_error_2)) THEN
            highest_error_2 = error
            index_highest_error_2 = i
        END IF

        ! ---------------------------------
        ! Third equation test
        ! ---------------------------------
        error = cauchy(i, cauchy_idx - 3) / (2.0D0 * mu) &
              + (s2 * factor_radius - s1 * factor_delta_r) * displacement(i, disp_idx + 1) &
              + (s2 * factor_radius + s1 * factor_delta_r) * displacement(i + 1, disp_idx + 1) &
              - cauchy_integral(i, cauchy_idx - 3)

        IF (ABS(error) > ABS(highest_error_3)) THEN
            highest_error_3 = error
            index_highest_error_3 = i
        END IF
    END DO

    ! Print the index and value of the largest error
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the rheology equations for degree j =", j, " and order m =", m
    PRINT '(A, I4)', "Largest error in first equation is at layer:", index_highest_error_1
    PRINT '(A, 2ES25.16)', "Value of the largest error in first equation:", highest_error_1
    PRINT '(A, I4)', "Largest error in second equation is at layer:", index_highest_error_2
    PRINT '(A, 2ES25.16)', "Value of the largest error in second equation:", highest_error_2
    PRINT '(A, I4)', "Largest error in third equation is at layer:", index_highest_error_3
    PRINT '(A, 2ES25.16)', "Value of the largest error in third equation:", highest_error_3

END SUBROUTINE

SUBROUTINE test_upper_boundary_condition(number_of_layers, jmax, j, m, ice_density, surface_g, cauchy, cauchy_isotropic, displacement)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: ice_density, surface_g
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_isotropic(number_of_layers, jmax + 1, jmax + 1)
    COMPLEX(KIND=8), INTENT(IN) :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)

    ! Declare local variables
    REAL(KIND=8) :: j1, factor_density, factor_displacement, factor_displacement_root
    REAL(KIND=8) :: a, b, c, d, e, f
    COMPLEX(KIND=8) :: error_1, error_2
    INTEGER :: cauchy_idx, disp_idx

    ! Compute constants
    j1 = DBLE(j)
    factor_density = (ice_density * surface_g) / (2.0D0 * (2.0D0 * j1 + 1.0D0))
    factor_displacement = factor_density * j1
    factor_displacement_root = factor_density * SQRT(j1 * (j1 + 1.0D0))

    cauchy_idx = 5 * (j * (j + 1) / 2 + m)
    disp_idx = 3 * (j * (j + 1) / 2 + m)

    a = SQRT(j1 / (3.0D0 * (2.0D0 * j1 + 1.0D0)))
    b = SQRT((j1 + 1.0D0) * (2.0D0 * j1 + 3.0D0) / (6.0D0 * (2.0D0 * j1 + 1.0D0) * (2.0D0 * j1 - 1.0D0)))
    c = SQRT((j1 - 1.0D0) / (2.0D0 * j1 - 1.0D0))
    d = SQRT((j1 + 1.0D0) / (3.0D0 * (2.0D0 * j1 + 1.0D0)))
    e = SQRT(j1 * (2.0D0 * j1 - 1.0D0) / (6.0D0 * (2.0D0 * j1 + 1.0D0) * (2.0D0 * j1 + 3.0D0)))
    f = SQRT((j1 + 2.0D0) / (2.0D0 * j1 + 3.0D0))

    ! Initialize errors
    error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
    error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)


    ! ---------------------------------
    ! First equation test
    ! ---------------------------------
    error_1 = -a * cauchy_isotropic(number_of_layers, j + 1, m + 1) &
              -b * cauchy(number_of_layers, cauchy_idx - 5) &
              +c * cauchy(number_of_layers, cauchy_idx - 7)

    error_1 = error_1 + factor_displacement * displacement(number_of_layers, disp_idx - 1) &
              +factor_displacement * displacement(number_of_layers + 1, disp_idx - 1)

    error_1 = error_1 - factor_displacement_root * displacement(number_of_layers, disp_idx + 1) &
              -factor_displacement_root * displacement(number_of_layers + 1, disp_idx + 1)

    ! ---------------------------------
    ! Second equation test
    ! ---------------------------------
    error_2 = d * cauchy_isotropic(number_of_layers, j + 1, m + 1) &
              +e * cauchy(number_of_layers, cauchy_idx - 5) &
              -f * cauchy(number_of_layers, cauchy_idx - 3)

    error_2 = error_2 - factor_displacement_root * displacement(number_of_layers, disp_idx - 1) &
              -factor_displacement_root * displacement(number_of_layers + 1, disp_idx - 1)

    error_2 = error_2 + factor_density * (j1 + 1.0D0) * displacement(number_of_layers, disp_idx + 1) &
              +factor_density * (j1 + 1.0D0) * displacement(number_of_layers + 1, disp_idx + 1)

    ! Print results
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the Upper boundary condition for degree j =", j, " and order m =", m
    PRINT '(A, 2ES25.16)', "Value of the error in first equation:", error_1
    PRINT '(A, 2ES25.16)', "Value of the error in second equation:", error_2

END SUBROUTINE

SUBROUTINE test_lower_boundary_condition(number_of_layers, jmax, j, m, bottom_g, delta_rho, cauchy, cauchy_isotropic, displacement, bottom_force)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: bottom_g, delta_rho
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_isotropic(number_of_layers, jmax + 1, jmax + 1)
    COMPLEX(KIND=8), INTENT(IN) :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)
    COMPLEX(KIND=8), INTENT(IN) :: bottom_force(4)

    ! Declare local variables
    REAL(KIND=8) :: j1, factor_density, factor_displacement, factor_displacement_root
    REAL(KIND=8) :: a, b, c, d, e, f
    COMPLEX(KIND=8) :: error_1, error_2
    INTEGER :: cauchy_idx, disp_idx

    ! Precompute constants
    j1 = DBLE(j)
    factor_density = (delta_rho * bottom_g) / (2.0D0 * (2.0D0 * j1 + 1.0D0))
    factor_displacement = factor_density * j1
    factor_displacement_root = factor_density * SQRT(j1 * (j1 + 1.0D0))

    cauchy_idx = 5 * (j * (j + 1) / 2 + m)
    disp_idx = 3 * (j * (j + 1) / 2 + m)

    a = SQRT(j1 / (3.0D0 * (2.0D0 * j1 + 1.0D0)))
    b = SQRT((j1 + 1.0D0) * (2.0D0 * j1 + 3.0D0) / (6.0D0 * (2.0D0 * j1 + 1.0D0) * (2.0D0 * j1 - 1.0D0)))
    c = SQRT((j1 - 1.0D0) / (2.0D0 * j1 - 1.0D0))
    d = SQRT((j1 + 1.0D0) / (3.0D0 * (2.0D0 * j1 + 1.0D0)))
    e = SQRT(j1 * (2.0D0 * j1 - 1.0D0) / (6.0D0 * (2.0D0 * j1 + 1.0D0) * (2.0D0 * j1 + 3.0D0)))
    f = SQRT((j1 + 2.0D0) / (2.0D0 * j1 + 3.0D0))

    ! Initialize errors
    error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
    error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)


    ! ---------------------------------
    ! First equation test
    ! ---------------------------------
    error_1 = a * cauchy_isotropic(1, j + 1, m + 1) &
              + b * cauchy(1, cauchy_idx - 5) &
              - c * cauchy(1, cauchy_idx - 7)

    error_1 = error_1 + factor_displacement * (displacement(1, disp_idx - 1) + displacement(2, disp_idx - 1))

    error_1 = error_1 - factor_displacement_root * (displacement(1, disp_idx + 1) + displacement(2, disp_idx + 1))

    IF (j == 2) THEN
        SELECT CASE (m)
            CASE (0)
                error_1 = error_1 - bottom_force(1)
            CASE (2)
                error_1 = error_1 - bottom_force(3)
        END SELECT
    END IF

    ! ---------------------------------
    ! Second equation test
    ! ---------------------------------
    error_2 = -d * cauchy_isotropic(1, j + 1, m + 1) &
              - e * cauchy(1, cauchy_idx - 5) &
              + f * cauchy(1, cauchy_idx - 3)

    error_2 = error_2 - factor_displacement_root * (displacement(1, disp_idx - 1) + displacement(2, disp_idx - 1))

    error_2 = error_2 + factor_density * (j1 + 1.0D0) * (displacement(1, disp_idx + 1) + displacement(2, disp_idx + 1))

    IF (j == 2) THEN
        SELECT CASE (m)
            CASE (0)
                error_2 = error_2 - bottom_force(2)
            CASE (2)
                error_2 = error_2 - bottom_force(4)
        END SELECT
    END IF

    ! Print results
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the Lower boundary condition for degree j =", j, " and order m =", m
    PRINT '(A, 2ES25.16)', "Value of the error in first equation:", error_1
    PRINT '(A, 2ES25.16)', "Value of the error in the second equation:", error_2

END SUBROUTINE

SUBROUTINE test_toroidial_equation_of_motion(number_of_layers, jmax, j, m, radius, delta_r, cauchy)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: radius, delta_r
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)

    ! Declare variables
    REAL(KIND=8) :: alpha1, alpha2, beta1, beta2, j1, current_radius, delta_r_inv
    COMPLEX(KIND=8) :: error, highest_error
    INTEGER :: i, index_highest_error, cauchy_index

    ! Precompute constants
    j1 = DBLE(j)

    alpha1 = SQRT((j1 - 1.0D0) / (2.0D0 * (2.0D0 * j1 + 1.0D0)))
    alpha2 = -(j1 - 1.0D0) * alpha1
    beta1 = -SQRT((j1 + 2.0D0) / (2.0D0 * (2.0D0 * j1 + 1.0D0)))
    beta2 = (j1 + 2.0D0) * beta1
    delta_r_inv = 1.0D0 / delta_r

    ! Initialize error tracking
    highest_error = CMPLX(0.0D0, 0.0D0, KIND=8)
    index_highest_error = 0

    ! Compute the base index for cauchy array once
    cauchy_index = 5 * (j * (j + 1) / 2 + m) - 5

    ! Loop through layers
    DO i = 1, number_of_layers - 1

        current_radius = radius + (i - 0.5D0) * delta_r

        error = ((alpha2 / (2.0D0 * current_radius)) - (alpha1 * delta_r_inv)) * cauchy(i, cauchy_index - 1)
        error = error + ((beta2 / (2.0D0 * current_radius)) - (beta1 * delta_r_inv)) * cauchy(i, cauchy_index + 1)
        error = error + ((alpha2 / (2.0D0 * current_radius)) + (alpha1 * delta_r_inv)) * cauchy(i + 1, cauchy_index - 1)
        error = error + ((beta2 / (2.0D0 * current_radius)) + (beta1 * delta_r_inv)) * cauchy(i + 1, cauchy_index + 1)

        ! Track largest error
        IF (ABS(error) > ABS(highest_error)) THEN
            highest_error = error
            index_highest_error = i
        END IF

    END DO

    ! Print Results
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the toroidial equation of motion for degree j =", j, " and order m =", m
    PRINT '(A, I4)', "Largest error is at layer:", index_highest_error
    PRINT '(A, 2ES25.16)', "Value of the largest error:", highest_error

END SUBROUTINE

SUBROUTINE test_toroidial_rheology(number_of_layers, jmax, j, m, radius, delta_r, mu, cauchy, cauchy_integral, displacement)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: radius, delta_r, mu
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_integral(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)

    ! Declare variables
    REAL(KIND=8) :: gamma1, gamma2, delta1, delta2, j1, current_radius, delta_r_inv
    COMPLEX(KIND=8) :: error, highest_error_1, highest_error_2
    INTEGER :: i, index_highest_error_1, index_highest_error_2, cauchy_index, disp_index

    ! Compute constants
    j1 = DBLE(j)

    gamma1 = -SQRT((j1 - 1.0D0) / (2.0D0 * (2.0D0 * j1 + 1.0D0)))
    gamma2 = (j1 + 1.0D0) * gamma1
    delta1 = SQRT((j1 + 2.0D0) / (2.0D0 * (2.0D0 * j1 + 1.0D0)))
    delta2 = -j1 * delta1
    delta_r_inv = 1.0D0 / delta_r

    ! Initialize error tracking
    highest_error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
    highest_error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)
    index_highest_error_1 = 0
    index_highest_error_2 = 0

    ! Compute base indices for cauchy and displacement arrays
    cauchy_index = 5 * (j * (j + 1) / 2 + m) - 5
    disp_index = 3 * (j * (j + 1) / 2 + m)

    ! ---------------------------------
    ! First equation test
    ! ---------------------------------
    DO i = 1, number_of_layers

        current_radius = radius + (i - 1.0D0) * delta_r

        error = cauchy(i, cauchy_index - 1) / (2.0D0 * mu)
        error = error + (((gamma2 / (2.0D0 * current_radius)) - (gamma1 * delta_r_inv)) * displacement(i, disp_index))
        error = error + (((gamma2 / (2.0D0 * current_radius)) + (gamma1 * delta_r_inv)) * displacement(i + 1, disp_index))
        error = error - cauchy_integral(i, cauchy_index - 1)

        ! Track largest error
        IF (ABS(error) > ABS(highest_error_1)) THEN
            highest_error_1 = error
            index_highest_error_1 = i
        END IF

    END DO

    ! ---------------------------------
    ! Second equation test
    ! ---------------------------------
    DO i = 1, number_of_layers

        current_radius = radius + (i - 1.0D0) * delta_r

        error = cauchy(i, cauchy_index + 1) / (2.0D0 * mu)
        error = error + (((delta2 / (2.0D0 * current_radius)) - (delta1 * delta_r_inv)) * displacement(i, disp_index))
        error = error + (((delta2 / (2.0D0 * current_radius)) + (delta1 * delta_r_inv)) * displacement(i + 1, disp_index))
        error = error - cauchy_integral(i, cauchy_index + 1)

        ! Track largest error
        IF (ABS(error) > ABS(highest_error_2)) THEN
            highest_error_2 = error
            index_highest_error_2 = i
        END IF

    END DO

    ! Print the index and value of the largest error
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the toroidial rheology equations for degree j =", j, " and order m =", m
    PRINT '(A, I4)', "Largest error in first equation is at layer:", index_highest_error_1
    PRINT '(A, 2ES25.16)', "Value of the largest error in first equation:", highest_error_1
    PRINT '(A, I4)', "Largest error in second equation is at layer:", index_highest_error_2
    PRINT '(A, 2ES25.16)', "Value of the largest error in second equation:", highest_error_2

END SUBROUTINE

SUBROUTINE test_toroidial_boundary_conditions(number_of_layers, jmax, j, m, cauchy)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)

    ! Declare variables
    REAL(KIND=8) :: a, b, j1
    COMPLEX(KIND=8) :: error_1, error_2
    INTEGER :: cauchy_idx

    ! Compute constants
    j1 = DBLE(j)

    a = SQRT((j1 - 1.0D0) / (2.0D0 * (2.0D0 * j1 + 1.0D0)))
    b = SQRT((j1 + 2.0D0) / (2.0D0 * (2.0D0 * j1 + 1.0D0)))

    ! Compute the base index for the cauchy array
    cauchy_idx = 5 * (j * (j + 1) / 2 + m) - 5

    ! Initialize errors
    error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
    error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)

    ! ---------------------------------
    ! Lower boundary condition
    ! ---------------------------------
    error_1 = a * cauchy(1, cauchy_idx - 1) - b * cauchy(1, cauchy_idx + 1)

    ! ---------------------------------
    ! Upper boundary condition
    ! ---------------------------------
    error_2 = a * cauchy(number_of_layers, cauchy_idx - 1) - b * cauchy(number_of_layers, cauchy_idx + 1)

    ! Print the results
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the toroidial boundary conditions for degree j =", j, " and order m =", m
    PRINT '(A, 2ES25.16)', "Value of the error in the lower boundary condition:", error_1
    PRINT '(A, 2ES25.16)', "Value of the error in the upper boundary condition:", error_2

END SUBROUTINE

SUBROUTINE test_equation_of_continuity_for_j1(number_of_layers, jmax, j, m, radius, delta_r, displacement)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: radius, delta_r
    COMPLEX(KIND=8), INTENT(IN) :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)

    ! Declare variables
    REAL(KIND=8) :: a1, b1, b2, delta_r_inv, current_radius
    COMPLEX(KIND=8) :: error, highest_error
    INTEGER :: i, index_highest_error, base_index

    ! Compute constants
    a1 = SQRT(1.0D0 / 3.0D0)
    b1 = -SQRT(2.0D0 / 3.0D0)
    b2 = 3.0D0 * b1
    delta_r_inv = 1.0D0 / delta_r

    ! Precompute the base index for displacement array
    base_index = 3 * (j * (j + 1) / 2 + m)

    ! Initialize values
    highest_error = CMPLX(0.0D0, 0.0D0, KIND=8)
    index_highest_error = 0

    ! Loop through each layer to calculate errors
    DO i = 1, number_of_layers

        current_radius = radius + (i - 1.0D0) * delta_r

        ! Initialize error for the current layer
        error = CMPLX(0.0D0, 0.0D0, KIND=8)

        ! Accumulate error terms
        error = error - (a1 * delta_r_inv) * displacement(i, base_index - 1)
        error = error + ((b2 / (2.0D0 * current_radius)) - (b1 * delta_r_inv)) * displacement(i, base_index + 1)
        error = error + (a1 * delta_r_inv) * displacement(i + 1, base_index - 1)
        error = error + ((b2 / (2.0D0 * current_radius)) + (b1 * delta_r_inv)) * displacement(i + 1, base_index + 1)

        ! Check if this error is the largest so far
        IF (ABS(error) > ABS(highest_error)) THEN
            highest_error = error
            index_highest_error = i
        END IF

    END DO

    ! Print the index and value of the largest error
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the equation of continuity for degree j =", j, " and order m =", m
    PRINT '(A, I4)', "Largest error is at layer:", index_highest_error
    PRINT '(A, 2ES25.16)', "Value of the largest error:", highest_error

END SUBROUTINE

    
SUBROUTINE test_equation_of_motion_for_j1(number_of_layers, jmax, j, m, radius, delta_r, cauchy, cauchy_isotropic)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: radius, delta_r
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_isotropic(number_of_layers, jmax + 1, jmax + 1)

    ! Declare variables
    REAL(KIND=8) :: c1, c2, e1, e2, f1, f2, g1, g2, h1, h2, delta_r_inv, current_radius
    COMPLEX(KIND=8) :: error, highest_error_1, highest_error_2
    INTEGER :: i, index_highest_error_1, index_highest_error_2, cauchy_idx

    ! Compute constants
    c1 = -(1.0D0) / (3.0D0)
    c2 = 2.0D0 * c1
    e1 = -SQRT(5.0D0) / 3.0D0
    e2 = 2.0D0 * e1
    f1 = SQRT(2.0D0) / 3.0D0
    f2 = -f1
    g1 = SQRT(1.0D0 / 90.0D0)
    g2 = -g1
    h1 = -SQRT(3.0D0 / 5.0D0)
    h2 = 4.0D0 * h1
    delta_r_inv = 1.0D0 / delta_r

    ! Precompute the base index for cauchy array
    cauchy_idx = 3 * ((j * (j + 1)) / 2 + m)

    ! Initialize error tracking
    highest_error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
    highest_error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)
    index_highest_error_1 = 0
    index_highest_error_2 = 0

    ! ---------------------------------
    ! First equation test
    ! ---------------------------------
    DO i = 1, number_of_layers - 1

        current_radius = radius + (i - 0.5D0) * delta_r

        error = CMPLX(0.0D0, 0.0D0, KIND=8)
        error = error + (((c2 / (2.0D0 * current_radius)) - (c1 * delta_r_inv)) * cauchy_isotropic(i, j + 1, m + 1))
        error = error + (((e2 / (2.0D0 * current_radius)) - (e1 * delta_r_inv)) * cauchy(i, cauchy_idx - 1))
        error = error + (((c2 / (2.0D0 * current_radius)) + (c1 * delta_r_inv)) * cauchy_isotropic(i + 1, j + 1, m + 1))
        error = error + (((e2 / (2.0D0 * current_radius)) + (e1 * delta_r_inv)) * cauchy(i + 1, cauchy_idx - 1))

        ! Track largest error
        IF (ABS(error) > ABS(highest_error_1)) THEN
            highest_error_1 = error
            index_highest_error_1 = i
        END IF

    END DO

    ! ---------------------------------
    ! Second equation test
    ! ---------------------------------
    DO i = 1, number_of_layers - 1

        current_radius = radius + (i - 0.5D0) * delta_r

        error = CMPLX(0.0D0, 0.0D0, KIND=8)
        error = error + (((f2 / (2.0D0 * current_radius)) - (f1 * delta_r_inv)) * cauchy_isotropic(i, j + 1, m + 1))
        error = error + (((g2 / (2.0D0 * current_radius)) - (g1 * delta_r_inv)) * cauchy(i, cauchy_idx - 1))
        error = error + (((h2 / (2.0D0 * current_radius)) - (h1 * delta_r_inv)) * cauchy(i, cauchy_idx + 1))
        error = error + (((f2 / (2.0D0 * current_radius)) + (f1 * delta_r_inv)) * cauchy_isotropic(i + 1, j + 1, m + 1))
        error = error + (((g2 / (2.0D0 * current_radius)) + (g1 * delta_r_inv)) * cauchy(i + 1, cauchy_idx - 1))
        error = error + (((h2 / (2.0D0 * current_radius)) + (h1 * delta_r_inv)) * cauchy(i + 1, cauchy_idx + 1))

        ! Track largest error
        IF (ABS(error) > ABS(highest_error_2)) THEN
            highest_error_2 = error
            index_highest_error_2 = i
        END IF

    END DO

    ! ---------------------------------
    ! Print Results
    ! ---------------------------------
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the equation of motion for degree j =", j, " and order m =", m
    PRINT '(A, I4)', "Largest error in first equation is at layer:", index_highest_error_1
    PRINT '(A, 2ES25.16)', "Value of the largest error in first equation:", highest_error_1
    PRINT '(A, I4)', "Largest error in second equation is at layer:", index_highest_error_2
    PRINT '(A, 2ES25.16)', "Value of the largest error in second equation:", highest_error_2

END SUBROUTINE

SUBROUTINE test_rheology_for_j1(number_of_layers, jmax, j, m, radius, delta_r, mu, cauchy, cauchy_integral, displacement)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: radius, delta_r, mu
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_integral(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)

    ! Declare variables
    REAL(KIND=8) :: q1, r1, r2, s1, s2, delta_r_inv, current_radius
    COMPLEX(KIND=8) :: error_1, error_2, highest_error_1, highest_error_2
    INTEGER :: i, index_highest_error_1, index_highest_error_2, base_index

    ! Compute constants
    q1 = SQRT(5.0D0) / 3.0D0
    r1 = -SQRT(1.0D0 / 90.0D0)
    r2 = 3.0D0 * r1
    s1 = SQRT(3.0D0 / 5.0D0)
    s2 = -2.0D0 * s1
    delta_r_inv = 1.0D0 / delta_r

    ! Precompute the base index for displacement array
    base_index = 3 * ((j * (j + 1)) / 2 + m)

    ! Initialize error tracking
    highest_error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
    highest_error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)
    index_highest_error_1 = 0
    index_highest_error_2 = 0

    ! ---------------------------------
    ! Combined loop for both equations
    ! ---------------------------------
    DO i = 1, number_of_layers

        current_radius = radius + (i - 1.0D0) * delta_r

        ! Initialize errors for the current layer
        error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
        error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)

        ! First equation error accumulation
        error_1 = error_1 + cauchy(i, base_index - 1) / (2.0D0 * mu)
        error_1 = error_1 + -(q1 * delta_r_inv) * displacement(i, base_index - 1)
        error_1 = error_1 + ((r2 / (2.0D0 * current_radius)) - (r1 * delta_r_inv)) * displacement(i, base_index + 1)
        error_1 = error_1 + (q1 * delta_r_inv) * displacement(i + 1, base_index - 1)
        error_1 = error_1 + ((r2 / (2.0D0 * current_radius)) + (r1 * delta_r_inv)) * displacement(i + 1, base_index + 1)
        error_1 = error_1 - cauchy_integral(i, base_index - 1)

        ! Track largest error for the first equation
        IF (ABS(error_1) > ABS(highest_error_1)) THEN
            highest_error_1 = error_1
            index_highest_error_1 = i
        END IF

        ! Second equation error accumulation
        error_2 = error_2 + cauchy(i, base_index + 1) / (2.0D0 * mu)
        error_2 = error_2 + ((s2 / (2.0D0 * current_radius)) - (s1 * delta_r_inv)) * displacement(i, base_index + 1)
        error_2 = error_2 + ((s2 / (2.0D0 * current_radius)) + (s1 * delta_r_inv)) * displacement(i + 1, base_index + 1)
        error_2 = error_2 - cauchy_integral(i, base_index + 1)

        ! Track largest error for the second equation
        IF (ABS(error_2) > ABS(highest_error_2)) THEN
            highest_error_2 = error_2
            index_highest_error_2 = i
        END IF

    END DO

    ! Print the index and value of the largest errors
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the rheology equations for degree j =", j, " and order m =", m
    PRINT '(A, I4)', "Largest error in first equation is at layer:", index_highest_error_1
    PRINT '(A, 2ES25.16)', "Value of the largest error in first equation:", highest_error_1
    PRINT '(A, I4)', "Largest error in second equation is at layer:", index_highest_error_2
    PRINT '(A, 2ES25.16)', "Value of the largest error in second equation:", highest_error_2

END SUBROUTINE

SUBROUTINE test_upper_boundary_condition_for_j1(number_of_layers, jmax, j, m, ice_density, surface_g, cauchy, cauchy_isotropic, displacement)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: ice_density, surface_g
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_isotropic(number_of_layers, jmax + 1, jmax + 1)
    COMPLEX(KIND=8), INTENT(IN) :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)

    ! Declare variables
    REAL(KIND=8) :: a, b, d, e, f, coeff_1, coeff_2
    COMPLEX(KIND=8) :: error_1, error_2
    INTEGER :: base_index

    ! Compute constants
    a = SQRT(1.0D0 / 9.0D0)
    b = SQRT(5.0D0 / 9.0D0)
    d = SQRT(2.0D0 / 9.0D0)
    e = SQRT(1.0D0 / 90.0D0)
    f = SQRT(3.0D0 / 5.0D0)

    coeff_1 = (ice_density * surface_g) / 6.0D0
    coeff_2 = coeff_1 * SQRT(2.0D0)

    ! Precompute the base index for the displacement array
    base_index = 3 * (j * (j + 1) / 2 + m)

    ! Initialize errors
    error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
    error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)

    ! ---------------------------------
    ! First equation test
    ! ---------------------------------
    error_1 = error_1 - a * cauchy_isotropic(number_of_layers, j + 1, m + 1)
    error_1 = error_1 - b * cauchy(number_of_layers, base_index - 1)
    error_1 = error_1 + coeff_1 * (displacement(number_of_layers, base_index - 1) + displacement(number_of_layers + 1, base_index - 1))
    error_1 = error_1 - coeff_2 * (displacement(number_of_layers, base_index + 1) + displacement(number_of_layers + 1, base_index + 1))

    ! ---------------------------------
    ! Second equation test
    ! ---------------------------------
    error_2 = error_2 + d * cauchy_isotropic(number_of_layers, j + 1, m + 1)
    error_2 = error_2 + e * cauchy(number_of_layers, base_index - 1)
    error_2 = error_2 - f * cauchy(number_of_layers, base_index + 1)
    error_2 = error_2 - coeff_2 * (displacement(number_of_layers, base_index - 1) + displacement(number_of_layers + 1, base_index - 1))
    error_2 = error_2 + 2.0D0 * coeff_1 * (displacement(number_of_layers, base_index + 1) + displacement(number_of_layers + 1, base_index + 1))

    ! Print the index and value of the largest error
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the Upper boundary condition for degree j =", j, " and order m =", m
    PRINT '(A, 2ES25.16)', "Value of the error in first equation:", error_1
    PRINT '(A, 2ES25.16)', "Value of the error in the second equation:", error_2

END SUBROUTINE

SUBROUTINE test_lower_boundary_condition_for_j1(number_of_layers, jmax, j, m, bottom_g, delta_rho, cauchy, cauchy_isotropic, displacement)
    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax, j, m
    REAL(KIND=8), INTENT(IN) :: bottom_g, delta_rho
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_isotropic(number_of_layers, jmax + 1, jmax + 1)
    COMPLEX(KIND=8), INTENT(IN) :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)

    ! Declare variables
    REAL(KIND=8) :: a, b, d, e, f, coeff_1, coeff_2
    COMPLEX(KIND=8) :: error_1, error_2
    INTEGER :: base_index

    ! Compute constants
    a = SQRT(1.0D0 / 9.0D0)
    b = SQRT(5.0D0 / 9.0D0)
    d = SQRT(2.0D0 / 9.0D0)
    e = SQRT(1.0D0 / 90.0D0)
    f = SQRT(3.0D0 / 5.0D0)

    coeff_1 = (delta_rho * bottom_g) / 6.0D0
    coeff_2 = coeff_1 * SQRT(2.0D0)

    ! Precompute the base index for the displacement array
    base_index = 3 * (j * (j + 1) / 2 + m)

    ! Initialize errors
    error_1 = CMPLX(0.0D0, 0.0D0, KIND=8)
    error_2 = CMPLX(0.0D0, 0.0D0, KIND=8)

    ! ---------------------------------
    ! First equation test
    ! ---------------------------------
    error_1 = error_1 + a * cauchy_isotropic(1, j + 1, m + 1)
    error_1 = error_1 + b * cauchy(1, base_index - 1)
    error_1 = error_1 + coeff_1 * (displacement(1, base_index - 1) + displacement(2, base_index - 1))
    error_1 = error_1 - coeff_2 * (displacement(1, base_index + 1) + displacement(2, base_index + 1))

    ! ---------------------------------
    ! Second equation test
    ! ---------------------------------
    error_2 = error_2 - d * cauchy_isotropic(1, j + 1, m + 1)
    error_2 = error_2 - e * cauchy(1, base_index - 1)
    error_2 = error_2 + f * cauchy(1, base_index + 1)
    error_2 = error_2 - coeff_2 * (displacement(1, base_index - 1) + displacement(2, base_index - 1))
    error_2 = error_2 + 2.0D0 * coeff_1 * (displacement(1, base_index + 1) + displacement(2, base_index + 1))

    ! Print the index and value of the largest error
    PRINT *
    PRINT '(A, I2, A, I2)', "Test of the Lower boundary condition for degree j =", j, " and order m =", m
    PRINT '(A, 2ES25.16)', "Value of the error in first equation:", error_1
    PRINT '(A, 2ES25.16)', "Value of the error in the second equation:", error_2

END SUBROUTINE

SUBROUTINE test_equations_solution(number_of_layers, jmax, radius, delta_r, mu, ice_density, surface_g, bottom_g, delta_rho, volume_force, bottom_force, cauchy_integral, cauchy, cauchy_isotropic, displacement)

    IMPLICIT NONE

    ! Declare input arguments
    INTEGER, INTENT(IN) :: number_of_layers, jmax
    REAL(KIND=8), INTENT(IN) :: radius, delta_r, mu, ice_density, surface_g, bottom_g, delta_rho
    COMPLEX(KIND=8), INTENT(IN) :: volume_force(number_of_layers - 1, 2)
    COMPLEX(KIND=8), INTENT(IN) :: bottom_force(4)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_integral(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    COMPLEX(KIND=8), INTENT(IN) :: cauchy_isotropic(number_of_layers, jmax + 1, jmax + 1)
    COMPLEX(KIND=8), INTENT(IN) :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)

    ! Local variables
    INTEGER :: j, m

    ! Test for j = 1
    DO m = 0, 1
        CALL test_equation_of_continuity_for_j1(number_of_layers, jmax, 1, m, radius, delta_r, displacement)
        CALL test_equation_of_motion_for_j1(number_of_layers, jmax, 1, m, radius, delta_r, cauchy, cauchy_isotropic)
        CALL test_rheology_for_j1(number_of_layers, jmax, 1, m, radius, delta_r, mu, cauchy, cauchy_integral, displacement)
        CALL test_upper_boundary_condition_for_j1(number_of_layers, jmax, 1, m, ice_density, surface_g, cauchy, cauchy_isotropic, displacement)
        CALL test_lower_boundary_condition_for_j1(number_of_layers, jmax, 1, m, bottom_g, delta_rho, cauchy, cauchy_isotropic, displacement)
    END DO

    ! Test for j >= 2
    DO j = 2, jmax
        DO m = 0, j
            CALL test_equation_of_continuity(number_of_layers, jmax, j, m, radius, delta_r, displacement)
            CALL test_equation_of_motion(number_of_layers, jmax, j, m, radius, delta_r, cauchy, cauchy_isotropic, volume_force)
            CALL test_rheology(number_of_layers, jmax, j, m, radius, delta_r, mu, cauchy, cauchy_integral, displacement)
            CALL test_upper_boundary_condition(number_of_layers, jmax, j, m, ice_density, surface_g, cauchy, cauchy_isotropic, displacement)
            CALL test_lower_boundary_condition(number_of_layers, jmax, j, m, bottom_g, delta_rho, cauchy, cauchy_isotropic, displacement, bottom_force)
            CALL test_toroidial_equation_of_motion(number_of_layers, jmax, j, m, radius, delta_r, cauchy)
            CALL test_toroidial_rheology(number_of_layers, jmax, j, m, radius, delta_r, mu, cauchy, cauchy_integral, displacement)
            CALL test_toroidial_boundary_conditions(number_of_layers, jmax, j, m, cauchy)
        END DO
    END DO

END SUBROUTINE

SUBROUTINE calculate_radial_displacement_at_layer(i, jmax, number_of_layers, displacement, radial_displacement)

    integer :: i, jmax, number_of_layers
    complex*16 :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)
    complex*16 :: radial_displacement(jmax+1, jmax+1)

    real*8 j1
    integer :: j, m

    j1 = real(j)

    do j=1, jmax
    
        do m=0, j

            radial_displacement(j+1,m+1) = 1/(sqrt(2*j1+1))*(sqrt(j1))*(displacement(i, 3 * (j * (j + 1) / 2 + m) - 1)+displacement(i+1, 3 * (j * (j + 1) / 2 + m) - 1))/2.0 - 1/(sqrt(2*j1+1))*(sqrt(j1+1))*(displacement(i, 3 * (j * (j + 1) / 2 + m) + 1)+displacement(i+1, 3 * (j * (j + 1) / 2 + m) + 1))/2.0

        end do
        
    end do

    END SUBROUTINE

SUBROUTINE fill_matrix(number_of_layers, matrix, j, mu, ice_density, radius, delta_r, delta_rho, surface_g, bottom_g)

    INTEGER number_of_layers
    REAL*8 matrix(6*number_of_layers+2,6*number_of_layers+2), j, mu, ice_density, radius, delta_r, delta_rho, surface_g, bottom_g
    REAL*8 :: a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2, g1, g2, h1, h2, q1, q2, r1, r2, p1, p2, s1, s2

    a1=sqrt((j)/(2*j+1))
    a2=-(j-1)*a1
    b1=-sqrt((j+1)/(2*j+1))
    b2=(j+2)*b1
    c1=-sqrt((j)/(3*(2*j+1)))
    c2=(j+1)*c1
    d1=sqrt((j-1)/(2*j-1))
    d2=-(j-2)*d1
    e1=-sqrt(((j+1)*(2*j+3))/(6*(2*j-1)*(2*j+1)))
    e2=(j+1)*e1
    f1=sqrt((j+1)/(3*(2*j+1)))
    f2=-j*f1
    g1=sqrt((j*(2*j-1))/(6*(2*j+1)*(2*j+3)))
    g2=-j*g1
    h1=-sqrt((j+2)/(2*j+3))
    h2=(j+3)*h1
    q1=sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))
    q2=-(j-1)*q1
    r1=-sqrt((j*(2*j-1))/(6*(2*j+1)*(2*j+3)))
    r2=(j+2)*r1
    p1=-sqrt((j-1)/(2*j-1))
    p2=j*p1
    s1=sqrt((j+2)/(2*j+3))
    s2=-(j+1)*s1

    matrix(1,1)=((a2)/(2*radius)-(a1)/(delta_r))
    matrix(1,2)=((b2)/(2*radius)-(b1)/(delta_r))
    matrix(1,7)=((a2)/(2*radius)+(a1)/(delta_r))
    matrix(1,8)=((b2)/(2*radius)+(b1)/(delta_r))

    matrix(2,3)=sqrt((j)/(3*(2*j+1)))
    matrix(2,4)=-sqrt((j-1)/(2*j-1))
    matrix(2,5)=sqrt(((j+1)*(2*j+3))/(6*(2*j-1)*(2*j+1)))
    matrix(2,1)=(delta_rho*bottom_g*j)/((2*j+1)*2)
    matrix(2,7)=(delta_rho*bottom_g*j)/((2*j+1)*2)
    matrix(2,2)=-(delta_rho*bottom_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix(2,8)=-(delta_rho*bottom_g*sqrt(j*(j+1)))/((2*j+1)*2)

    matrix(3,3)=-sqrt((j+1)/(3*(2*j+1)))
    matrix(3,5)=-sqrt((j*(2*j-1))/(6*(2*j+1)*(2*j+3)))
    matrix(3,6)=sqrt((j+2)/(2*j+3))
    matrix(3,1)=-(delta_rho*bottom_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix(3,7)=-(delta_rho*bottom_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix(3,2)=(delta_rho*bottom_g*(j+1))/((2*j+1)*2)
    matrix(3,8)=(delta_rho*bottom_g*(j+1))/((2*j+1)*2)

    matrix(4,1)=((p2)/(2*radius)-(p1)/(delta_r))
    matrix(4,7)=((p2)/(2*radius)+(p1)/(delta_r))
    matrix(4,4)=1.0/(2*mu)

    matrix(5,1)=((q2)/(2*radius)-(q1)/(delta_r))
    matrix(5,2)=((r2)/(2*radius)-(r1)/(delta_r))
    matrix(5,5)=1.0/(2*mu)
    matrix(5,7)=((q2)/(2*radius)+(q1)/(delta_r))
    matrix(5,8)=((r2)/(2*radius)+(r1)/(delta_r))

    matrix(6,2)=((s2)/(2*radius)-(s1)/(delta_r))
    matrix(6,8)=((s2)/(2*radius)+(s1)/(delta_r))
    matrix(6,6)=1.0/(2*mu)


    do i=2,number_of_layers

        matrix(6*(i-1)+1,6*(i-1)+1)=((a2)/(2*(radius+(i-1.0)*delta_r))-(a1)/(delta_r))
        matrix(6*(i-1)+1,6*(i-1)+2)=((b2)/(2*(radius+(i-1.0)*delta_r))-(b1)/(delta_r))
        matrix(6*(i-1)+1,6*i+1)=((a2)/(2*(radius+(i-1.0)*delta_r))+(a1)/(delta_r))
        matrix(6*(i-1)+1,6*i+2)=((b2)/(2*(radius+(i-1.0)*delta_r))+(b1)/(delta_r))

        matrix(6*(i-1)+2,6*(i-1)-3)=((c2)/(2*(radius+(i-1.5)*delta_r))-(c1)/(delta_r))
        matrix(6*(i-1)+2,6*(i-1)-2)=((d2)/(2*(radius+(i-1.5)*delta_r))-(d1)/(delta_r))
        matrix(6*(i-1)+2,6*(i-1)-1)=((e2)/(2*(radius+(i-1.5)*delta_r))-(e1)/(delta_r))
        matrix(6*(i-1)+2,6*(i-1)+3)=((c2)/(2*(radius+(i-1.5)*delta_r))+(c1)/(delta_r))
        matrix(6*(i-1)+2,6*(i-1)+4)=((d2)/(2*(radius+(i-1.5)*delta_r))+(d1)/(delta_r))
        matrix(6*(i-1)+2,6*(i-1)+5)=((e2)/(2*(radius+(i-1.5)*delta_r))+(e1)/(delta_r))

        matrix(6*(i-1)+3,6*(i-1)-3)=((f2)/(2*(radius+(i-1.5)*delta_r))-(f1)/(delta_r))
        matrix(6*(i-1)+3,6*(i-1)-1)=((g2)/(2*(radius+(i-1.5)*delta_r))-(g1)/(delta_r))
        matrix(6*(i-1)+3,6*(i-1))=((h2)/(2*(radius+(i-1.5)*delta_r))-(h1)/(delta_r))
        matrix(6*(i-1)+3,6*(i-1)+3)=((f2)/(2*(radius+(i-1.5)*delta_r))+(f1)/(delta_r))
        matrix(6*(i-1)+3,6*(i-1)+5)=((g2)/(2*(radius+(i-1.5)*delta_r))+(g1)/(delta_r))
        matrix(6*(i-1)+3,6*i)=((h2)/(2*(radius+(i-1.5)*delta_r))+(h1)/(delta_r))

        matrix(6*(i-1)+4,6*(i-1)+1)=((p2)/(2*(radius+(i-1.0)*delta_r))-(p1)/(delta_r))
        matrix(6*(i-1)+4,6*i+1)=((p2)/(2*(radius+(i-1.0)*delta_r))+(p1)/(delta_r))
        matrix(6*(i-1)+4,6*(i-1)+4)=1.0/(2*mu)

        matrix(6*(i-1)+5,6*(i-1)+1)=((q2)/(2*(radius+(i-1.0)*delta_r))-(q1)/(delta_r))
        matrix(6*(i-1)+5,6*(i-1)+2)=((r2)/(2*(radius+(i-1.0)*delta_r))-(r1)/(delta_r))
        matrix(6*(i-1)+5,6*i+1)=((q2)/(2*(radius+(i-1.0)*delta_r))+(q1)/(delta_r))
        matrix(6*(i-1)+5,6*i+2)=((r2)/(2*(radius+(i-1.0)*delta_r))+(r1)/(delta_r))
        matrix(6*(i-1)+5,6*(i-1)+5)=1.0/(2*mu)

        matrix(6*(i-1)+6,6*(i-1)+2)=((s2)/(2*(radius+(i-1.0)*delta_r))-(s1)/(delta_r))
        matrix(6*(i-1)+6,6*i+2)=((s2)/(2*(radius+(i-1.0)*delta_r))+(s1)/(delta_r))
        matrix(6*(i-1)+6,6*i)=1.0/(2*mu)

    end do

    matrix(6*number_of_layers+1,6*number_of_layers-3)=-sqrt((j)/(3*(2*j+1)))
    matrix(6*number_of_layers+1,6*number_of_layers-2)=sqrt((j-1)/(2*j-1))
    matrix(6*number_of_layers+1,6*number_of_layers-1)=-sqrt(((j+1)*(2*j+3))/(6*(2*j-1)*(2*j+1)))
    matrix(6*number_of_layers+1,6*number_of_layers-5)=(ice_density*surface_g*j)/((2*j+1)*2)
    matrix(6*number_of_layers+1,6*number_of_layers+1)=(ice_density*surface_g*j)/((2*j+1)*2)
    matrix(6*number_of_layers+1,6*number_of_layers-4)=-(ice_density*surface_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix(6*number_of_layers+1,6*number_of_layers+2)=-(ice_density*surface_g*sqrt(j*(j+1)))/((2*j+1)*2)

    matrix(6*number_of_layers+2,6*number_of_layers-3)=sqrt((j+1)/(3*(2*j+1)))
    matrix(6*number_of_layers+2,6*number_of_layers-1)=sqrt((j*(2*j-1))/(6*(2*j+1)*(2*j+3)))
    matrix(6*number_of_layers+2,6*number_of_layers)=-sqrt((j+2)/(2*j+3))
    matrix(6*number_of_layers+2,6*number_of_layers-5)=-(ice_density*surface_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix(6*number_of_layers+2,6*number_of_layers+1)=-(ice_density*surface_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix(6*number_of_layers+2,6*number_of_layers-4)=(ice_density*surface_g*(j+1))/((2*j+1)*2)
    matrix(6*number_of_layers+2,6*number_of_layers+2)=(ice_density*surface_g*(j+1))/((2*j+1)*2)

    END SUBROUTINE

SUBROUTINE fill_toroidial_matrix(number_of_layers, toroidial_matrix, j, mu, radius, delta_r)

    INTEGER number_of_layers
    REAL*8 toroidial_matrix(3*number_of_layers+1, 3*number_of_layers+1), j, mu, radius, delta_r
    REAL*8 :: alpha1, alpha2, beta1, beta2, gamma1, gamma2, delta1, delta2

    alpha1 = sqrt((j-1)/(2*(2*j+1)))
    alpha2 = -(j-1)*alpha1
    beta1 = -sqrt((j+2)/(2*(2*j+1)))
    beta2 = (j+2)*beta1
    gamma1 = -sqrt((j-1)/(2*(2*j+1)))
    gamma2 = (j+1)*gamma1
    delta1 = sqrt((j+2)/(2*(2*j+1)))
    delta2 = -j*delta1

    toroidial_matrix(1,2) = sqrt((j-1)/(2*(2*j+1)))
    toroidial_matrix(1,3) = -sqrt((j+2)/(2*(2*j+1)))

    toroidial_matrix(2,1) = ((gamma2)/(2*radius)-(gamma1)/(delta_r))
    toroidial_matrix(2,4) = ((gamma2)/(2*radius)+(gamma1)/(delta_r))
    toroidial_matrix(2,2) = 1.0/(2*mu)

    toroidial_matrix(3,1) = ((delta2)/(2*radius)-(delta1)/(delta_r))
    toroidial_matrix(3,4) = ((delta2)/(2*radius)+(delta1)/(delta_r))
    toroidial_matrix(3,3) = 1.0/(2*mu)

    do i=2,number_of_layers

        toroidial_matrix(3*(i-1)+1,3*(i-1)-1) = ((alpha2)/(2*(radius+(i-1.5)*delta_r))-(alpha1)/(delta_r))
        toroidial_matrix(3*(i-1)+1,3*(i-1)) = ((beta2)/(2*(radius+(i-1.5)*delta_r))-(beta1)/(delta_r))
        toroidial_matrix(3*(i-1)+1,3*(i-1)+2) = ((alpha2)/(2*(radius+(i-1.5)*delta_r))+(alpha1)/(delta_r))
        toroidial_matrix(3*(i-1)+1,3*(i-1)+3) = ((beta2)/(2*(radius+(i-1.5)*delta_r))+(beta1)/(delta_r))

        toroidial_matrix(3*(i-1)+2,3*(i-1)+1) = ((gamma2)/(2*(radius+(i-1.0)*delta_r))-(gamma1)/(delta_r))
        toroidial_matrix(3*(i-1)+2,3*(i-1)+4) = ((gamma2)/(2*(radius+(i-1.0)*delta_r))+(gamma1)/(delta_r))
        toroidial_matrix(3*(i-1)+2,3*(i-1)+2) = 1.0/(2*mu)

        toroidial_matrix(3*(i-1)+3,3*(i-1)+1) = ((delta2)/(2*(radius+(i-1.0)*delta_r))-(delta1)/(delta_r))
        toroidial_matrix(3*(i-1)+3,3*(i-1)+4) = ((delta2)/(2*(radius+(i-1.0)*delta_r))+(delta1)/(delta_r))
        toroidial_matrix(3*(i-1)+3,3*(i-1)+3) = 1.0/(2*mu)

    end do

    toroidial_matrix(3*number_of_layers + 1, 3*number_of_layers-1) = sqrt((j-1)/(2*(2*j+1)))
    toroidial_matrix(3*number_of_layers + 1, 3*number_of_layers) = -sqrt((j+2)/(2*(2*j+1)))

    END SUBROUTINE

SUBROUTINE fill_matrix_for_j1(number_of_layers, matrix_for_j1, mu, ice_density, radius, delta_r, delta_rho, surface_g, bottom_g)

    INTEGER number_of_layers
    REAL*8 matrix_for_j1(5*number_of_layers+2,5*number_of_layers+2), mu, ice_density, radius, delta_r, delta_rho, surface_g, bottom_g
    REAL*8 :: a1, b1, b2, c1, c2, e1, e2, f1, f2, g1, g2, h1, h2, q1, q2, r1, r2, s1, s2
    REAL*8 :: j

    j=1.0

    a1=sqrt((j)/(2*j+1))

    b1=-sqrt((j+1)/(2*j+1))
    b2=(j+2)*b1
    c1=-sqrt((j)/(3*(2*j+1)))
    c2=(j+1)*c1

    e1=-sqrt(((j+1)*(2*j+3))/(6*(2*j-1)*(2*j+1)))
    e2=(j+1)*e1
    f1=sqrt((j+1)/(3*(2*j+1)))
    f2=-j*f1
    g1=sqrt((j*(2*j-1))/(6*(2*j+1)*(2*j+3)))
    g2=-j*g1
    h1=-sqrt((j+2)/(2*j+3))
    h2=(j+3)*h1
    q1=sqrt((j+1)*(2*j+3)/(6*(2*j-1)*(2*j+1)))
    q2=-(j-1)*q1
    r1=-sqrt((j*(2*j-1))/(6*(2*j+1)*(2*j+3)))
    r2=(j+2)*r1

    s1=sqrt((j+2)/(2*j+3))
    s2=-(j+1)*s1

    matrix_for_j1(1,1)=-(a1)/(delta_r)
    matrix_for_j1(1,2)=((b2)/(2*radius)-(b1)/(delta_r))
    matrix_for_j1(1,6)=(a1)/(delta_r)
    matrix_for_j1(1,7)=((b2)/(2*radius)+(b1)/(delta_r))

    matrix_for_j1(2,3)=sqrt((j)/(3*(2*j+1)))
    matrix_for_j1(2,4)=sqrt(((j+1)*(2*j+3))/(6*(2*j-1)*(2*j+1)))
    matrix_for_j1(2,1)=(delta_rho*bottom_g*j)/((2*j+1)*2)
    matrix_for_j1(2,6)=(delta_rho*bottom_g*j)/((2*j+1)*2)
    matrix_for_j1(2,2)=-(delta_rho*bottom_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix_for_j1(2,7)=-(delta_rho*bottom_g*sqrt(j*(j+1)))/((2*j+1)*2)

    matrix_for_j1(3,3)=-sqrt((j+1)/(3*(2*j+1)))
    matrix_for_j1(3,4)=-sqrt((j*(2*j-1))/(6*(2*j+1)*(2*j+3)))
    matrix_for_j1(3,5)=sqrt((j+2)/(2*j+3))
    matrix_for_j1(3,1)=-(delta_rho*bottom_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix_for_j1(3,6)=-(delta_rho*bottom_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix_for_j1(3,2)=(delta_rho*bottom_g*(j+1))/((2*j+1)*2)
    matrix_for_j1(3,7)=(delta_rho*bottom_g*(j+1))/((2*j+1)*2)

    matrix_for_j1(4,1)=((q2)/(2*radius)-(q1)/(delta_r))
    matrix_for_j1(4,2)=((r2)/(2*radius)-(r1)/(delta_r))
    matrix_for_j1(4,4)=1.0/(2*mu)
    matrix_for_j1(4,6)=((q2)/(2*radius)+(q1)/(delta_r))
    matrix_for_j1(4,7)=((r2)/(2*radius)+(r1)/(delta_r))

    matrix_for_j1(5,2)=((s2)/(2*radius)-(s1)/(delta_r))
    matrix_for_j1(5,7)=((s2)/(2*radius)+(s1)/(delta_r))
    matrix_for_j1(5,5)=1.0/(2*mu)


    do i=2,number_of_layers

        matrix_for_j1(5*(i-1)+1,5*(i-1)+1)=-(a1)/(delta_r)
        matrix_for_j1(5*(i-1)+1,5*(i-1)+2)=((b2)/(2*(radius+(i-1.0)*delta_r))-(b1)/(delta_r))
        matrix_for_j1(5*(i-1)+1,5*i+1)=(a1)/(delta_r)
        matrix_for_j1(5*(i-1)+1,5*i+2)=((b2)/(2*(radius+(i-1.0)*delta_r))+(b1)/(delta_r))

        matrix_for_j1(5*(i-1)+2,5*(i-1)-2)=((c2)/(2*(radius+(i-1.5)*delta_r))-(c1)/(delta_r))
        matrix_for_j1(5*(i-1)+2,5*(i-1)-1)=((e2)/(2*(radius+(i-1.5)*delta_r))-(e1)/(delta_r))
        matrix_for_j1(5*(i-1)+2,5*(i-1)+3)=((c2)/(2*(radius+(i-1.5)*delta_r))+(c1)/(delta_r))
        matrix_for_j1(5*(i-1)+2,5*(i-1)+4)=((e2)/(2*(radius+(i-1.5)*delta_r))+(e1)/(delta_r))

        matrix_for_j1(5*(i-1)+3,5*(i-1)-2)=((f2)/(2*(radius+(i-1.5)*delta_r))-(f1)/(delta_r))
        matrix_for_j1(5*(i-1)+3,5*(i-1)-1)=((g2)/(2*(radius+(i-1.5)*delta_r))-(g1)/(delta_r))
        matrix_for_j1(5*(i-1)+3,5*(i-1))=((h2)/(2*(radius+(i-1.5)*delta_r))-(h1)/(delta_r))
        matrix_for_j1(5*(i-1)+3,5*(i-1)+3)=((f2)/(2*(radius+(i-1.5)*delta_r))+(f1)/(delta_r))
        matrix_for_j1(5*(i-1)+3,5*(i-1)+4)=((g2)/(2*(radius+(i-1.5)*delta_r))+(g1)/(delta_r))
        matrix_for_j1(5*(i-1)+3,5*i)=((h2)/(2*(radius+(i-1.5)*delta_r))+(h1)/(delta_r))

        matrix_for_j1(5*(i-1)+4,5*(i-1)+1)=((q2)/(2*(radius+(i-1.0)*delta_r))-(q1)/(delta_r))
        matrix_for_j1(5*(i-1)+4,5*(i-1)+2)=((r2)/(2*(radius+(i-1.0)*delta_r))-(r1)/(delta_r))
        matrix_for_j1(5*(i-1)+4,5*i+1)=((q2)/(2*(radius+(i-1.0)*delta_r))+(q1)/(delta_r))
        matrix_for_j1(5*(i-1)+4,5*i+2)=((r2)/(2*(radius+(i-1.0)*delta_r))+(r1)/(delta_r))
        matrix_for_j1(5*(i-1)+4,5*(i-1)+4)=1.0/(2*mu)

        matrix_for_j1(5*(i-1)+5,5*(i-1)+2)=((s2)/(2*(radius+(i-1.0)*delta_r))-(s1)/(delta_r))
        matrix_for_j1(5*(i-1)+5,5*i+2)=((s2)/(2*(radius+(i-1.0)*delta_r))+(s1)/(delta_r))
        matrix_for_j1(5*(i-1)+5,5*i)=1.0/(2*mu)

    end do

    matrix_for_j1(5*number_of_layers+1,5*number_of_layers-2)=-sqrt((j)/(3*(2*j+1)))
    matrix_for_j1(5*number_of_layers+1,5*number_of_layers-1)=-sqrt(((j+1)*(2*j+3))/(6*(2*j-1)*(2*j+1)))
    matrix_for_j1(5*number_of_layers+1,5*number_of_layers-4)=(ice_density*surface_g*j)/((2*j+1)*2)
    matrix_for_j1(5*number_of_layers+1,5*number_of_layers+1)=(ice_density*surface_g*j)/((2*j+1)*2)
    matrix_for_j1(5*number_of_layers+1,5*number_of_layers-3)=-(ice_density*surface_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix_for_j1(5*number_of_layers+1,5*number_of_layers+2)=-(ice_density*surface_g*sqrt(j*(j+1)))/((2*j+1)*2)

    matrix_for_j1(5*number_of_layers+2,5*number_of_layers-2)=sqrt((j+1)/(3*(2*j+1)))
    matrix_for_j1(5*number_of_layers+2,5*number_of_layers-1)=sqrt((j*(2*j-1))/(6*(2*j+1)*(2*j+3)))
    matrix_for_j1(5*number_of_layers+2,5*number_of_layers)=-sqrt((j+2)/(2*j+3))
    matrix_for_j1(5*number_of_layers+2,5*number_of_layers-4)=-(ice_density*surface_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix_for_j1(5*number_of_layers+2,5*number_of_layers+1)=-(ice_density*surface_g*sqrt(j*(j+1)))/((2*j+1)*2)
    matrix_for_j1(5*number_of_layers+2,5*number_of_layers-3)=(ice_density*surface_g*(j+1))/((2*j+1)*2)
    matrix_for_j1(5*number_of_layers+2,5*number_of_layers+2)=(ice_density*surface_g*(j+1))/((2*j+1)*2)

    END SUBROUTINE

SUBROUTINE calculate_forces(number_of_layers, t, volume_force, bottom_force, &
              radius, delta_r, angular_speed, delta_t, ice_density, delta_rho, excentricity)
    ! subroutine calculate_forces computes the volume and bottom forces
    ! applied on Europa's ice shell due to tidal deformations.
    !
    ! Arguments:
    !   number_of_layers - the number of layers in the simulation
    !   t - the current time step index
    !   jmax - maximum degree of spherical harmonics
    !   volume_force - array to hold the volume forces computed
    !   bottom_force - array to hold the bottom forces computed
    !   radius - the radius of Europa
    !   delta_r - the change in radius per layer
    !   angular_speed - angular speed of Europa's orbit
    !   delta_t - the time increment per step
    !   ice_density - the density of Europa's ice layer
    !   excentricity - the orbital excentricity of Europa

    ! Declare the types of the arguments
    integer, intent(in) :: number_of_layers, t
    complex*16, intent(out) :: volume_force(number_of_layers-1, 2), bottom_force(4)
    real*8, intent(in) :: radius, delta_r, angular_speed, delta_t, ice_density, delta_rho, excentricity

    ! Declare local variables
    real*8 :: konst, real_, imag_, t_, ksi
    integer :: i

    t_ = dble(t)

    ! Compute a constant factor used in the force calculations
    konst=ice_density*angular_speed*angular_speed*excentricity

    ! Loop over each layer to compute the volume forces
    do i=1,number_of_layers-1
        ksi=konst*(radius+(i-1/2.0)*delta_r)
        volume_force(i,1) = dcmplx(-ksi*sqrt(18*4.D0*datan(1.D0))*dcos(angular_speed*t_*delta_t), 0.D0)

        real_=ksi*sqrt(27*4.D0*datan(1.D0))*dcos(angular_speed*t_*delta_t)
        imag_=-ksi*sqrt(48*4.D0*datan(1.D0))*dsin(angular_speed*t_*delta_t)!<
        volume_force(i,2)=dcmplx(real_,imag_)
    end do

    ! Compute a modified constant factor for bottom force calculation
    ksi = (ice_density + delta_rho) * angular_speed**2 * excentricity * radius**2

    ! Compute the bottom forces at specific indices
    bottom_force(1)=dcmplx(-ksi*sqrt((18*4.D0*datan(1.D0))/(25))*dcos(angular_speed*t_*delta_t),0)
    bottom_force(2)=dcmplx(ksi*sqrt((27*4.D0*datan(1.D0))/(25))*dcos(angular_speed*t_*delta_t),0)

    real_=ksi*sqrt((27*4.D0*datan(1.D0))/(25))*dcos(angular_speed*t_*delta_t)
    imag_=-ksi*sqrt((48*4.D0*datan(1.D0))/(25))*dsin(angular_speed*t_*delta_t)
    bottom_force(3)=dcmplx(real_,imag_)

    real_=-ksi*sqrt((81*4.D0*datan(1.D0))/(50))*dcos(angular_speed*t_*delta_t)
    imag_=ksi*sqrt((72*4.D0*datan(1.D0))/(25))*dsin(angular_speed*t_*delta_t)
    bottom_force(4)=dcmplx(real_,imag_)

    end subroutine

subroutine update_cauchy_integral(jmax, number_of_layers, cauchy_integral, cauchy, delta_t, eta, mu)

    integer :: jmax, number_of_layers
    complex*16 :: cauchy_integral(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    complex*16 :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    real*8 :: delta_t, eta, mu

    integer :: j, m, k

    do j=1, jmax

        do m=0, j
        
            do k=1, number_of_layers
                
                cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 3) = cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 3) -&
                (1.0/eta) * cauchy(k, 5 * (j * (j + 1) / 2 + m) - 3) * delta_t
                cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 4) = cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 4) -&
                (1.0/eta) * cauchy(k, 5 * (j * (j + 1) / 2 + m) - 4) * delta_t
                cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 5) = cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 5)-&
                (1.0/eta) * cauchy(k, 5 * (j * (j + 1) / 2 + m) - 5) * delta_t
                cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 6) = cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 6)-&
                (1.0/eta) * cauchy(k, 5 * (j * (j + 1) / 2 + m) - 6) * delta_t
                cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 7) = cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 7)-&
                (1.0/eta) * cauchy(k, 5 * (j * (j + 1) / 2 + m) - 7) * delta_t

            end do

        end do

    end do

    end subroutine

subroutine solve_linear_system(number_of_layers, jmax, j, m, matrix,&
         volume_force, bottom_force, cauchy_integral, cauchy, cauchy_isotropic, displacement)
    implicit none

    integer, intent(in) :: number_of_layers, jmax, j, m
    real*8, intent(in) :: matrix(6*number_of_layers+2,6*number_of_layers+2)
    complex*16, intent(in) :: volume_force(number_of_layers-1,3*(jmax*(jmax+1)/2+jmax)+1), &
                        bottom_force(3*(jmax*(jmax+1)/2+jmax)+1), &
                        cauchy_integral(number_of_layers,5*(jmax*(jmax+1)/2+jmax)-3)
    complex*16, intent(inout) :: cauchy(number_of_layers,5*(jmax*(jmax+1)/2+jmax)-3), cauchy_isotropic(number_of_layers, jmax + 1, jmax +1), &
    displacement(number_of_layers+1,3*(jmax*(jmax+1)/2+jmax)+1)

    ! Declare workspace matrices and vectors with fixed bounds.
    real*8 :: A(1000,1000), AL(1000,1000)
    real*8 :: yr(6*number_of_layers+2), yi(6*number_of_layers+2), d

    ! Declare integer variables for indexing and dimensioning purposes.
    integer :: i, k, n, m1, m2, np, mp, mpl

    ! Allocate the index array for the LU decomposition.
    integer, allocatable :: INDX(:)

    ! Initialize the real and imaginary parts of the solution vectors.
    yr = 0.0d0
    yi=0.0d0

    ! Set the dimensions and bandwidths for the matrix A.
    n = 6*number_of_layers+2  ! Dimension of the matrix B.
    m1 = 7  ! Number of sub-diagonals in matrix B.
    m2 = 7  ! Number of super-diagonals in matrix B.
    np = 1000  ! Number of rows in the matrix A, which holds the band of matrix B.
    mp = 1000  ! Number of columns in the matrix A.
    mpl = 1000  ! Number of columns in the matrix AL (auxiliary matrix).

    ! Initialize the band matrix A to zero before populating with values from 'matrix'.
    do i = 1, n
        do k = 1, m1 + m2 + 1
            A(i, k) = 0d0
        end do
    end do


    ! Populate the band matrix A from 'matrix'.
    do i = 1, n
        do k = 1, m1 + m2 + 1
            if (i + k - m1 - 1 > 0 .and. i + k - m1 - 1 <= n) then
                A(i, k) = matrix(i, i + k - m1 - 1)
            endif
        end do
    end do

    ! Allocate the index array for LU decomposition.
    ALLOCATE(INDX(n))

    ! Perform LU decomposition on the band matrix A.
    call bandec (A, n, m1, m2, np, mp, AL, mpl, INDX, d)

    ! Prepare the right-hand side vector from 'volume_force' and 'bottom_force'.
    if (j == 2) then
        select case (m)
            case (0)
                do i = 1, number_of_layers - 1
                    yr(6 * i + 2) = -dreal(volume_force(i, 1))
                end do
            case (2)
                do i = 1, number_of_layers - 1
                    yr(6 * i + 2) = -dreal(volume_force(i, 2))
                    yi(6 * i + 2) = -dimag(volume_force(i, 2))
                end do
        end select
    end if

    ! Set the first two elements from 'bottom_force'.
    if (j == 2) then
        select case (m)
            case (0)
                yr(2) = dreal(bottom_force(1))
                yr(3) = dreal(bottom_force(2))
            case (2)
                yr(2) = dreal(bottom_force(3))
                yi(2) = dimag(bottom_force(3))
                yr(3) = dreal(bottom_force(4))
                yi(3) = dimag(bottom_force(4))
        end select
    end if

    do i = 1, number_of_layers
        yr(i*6-2) = yr(i*6-2) + dreal(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 7))
        yr(i*6-1) = yr(i*6-1) + dreal(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 5))
        yr(i*6) = yr(i*6) + dreal(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 3))
        yi(i*6-2) = yi(i*6-2) + dimag(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 7))
        yi(i*6-1) = yi(i*6-1) + dimag(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 5))
        yi(i*6) = yi(i*6) +dimag(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 3))
    end do


    ! Solve the system for the real and imaginary parts.
    call banbks(A, n, m1, m2, np, mp, AL, mpl, INDX, yr)  ! Solve for real part.
    call banbks(A, n, m1, m2, np, mp, AL, mpl, INDX, yi)  ! Solve for imaginary part.

    ! Deallocate the index array.
    DEALLOCATE(INDX)

    ! Store the solution into 'displacement' and 'cauchy' arrays.
    do i = 1, number_of_layers
    displacement(i, 3 * ((j * (j + 1)) / 2 + m) - 1) = dcmplx(yr(6 * (i - 1) + 1), yi(6 * (i - 1) + 1))
    displacement(i, 3 * ((j * (j + 1)) / 2 + m) + 1) = dcmplx(yr(6 * (i - 1) + 2), yi(6 * (i - 1) + 2))
    cauchy_isotropic(i, j+1, m+1) = dcmplx(yr(6 * (i - 1) + 3), yi(6 * (i - 1) + 3))
    cauchy(i, 5 * ((j * (j + 1)) / 2 + m) - 7) = dcmplx(yr(6 * (i - 1) + 4), yi(6 * (i - 1) + 4))
    cauchy(i, 5 * ((j * (j + 1)) / 2 + m) - 5) = dcmplx(yr(6 * (i - 1) + 5), yi(6 * (i - 1) + 5))
    cauchy(i, 5 * ((j * (j + 1)) / 2 + m) - 3) = dcmplx(yr(6 * (i - 1) + 6), yi(6 * (i - 1) + 6))
    end do

    ! Store the solution for the last layer into 'displacement'.
    displacement(number_of_layers + 1, 3 * ((j * (j + 1)) / 2 + m) - 1) = dcmplx(yr(6 * number_of_layers + 1),&
                                                                        yi(6 * number_of_layers + 1))
    displacement(number_of_layers + 1, 3 * ((j * (j + 1)) / 2 + m) + 1) = dcmplx(yr(6 * number_of_layers + 2),&
                                                                        yi(6 * number_of_layers + 2))


end subroutine

subroutine solve_toroidial_system(number_of_layers, jmax, j, m, toroidial_matrix, cauchy_integral, cauchy, displacement)
    implicit none

    integer, intent(in) :: number_of_layers, jmax, j, m
    real*8, intent(in) :: toroidial_matrix(3*number_of_layers+1,3*number_of_layers+1)
    complex*16, intent(in) :: cauchy_integral(number_of_layers,5*(jmax*(jmax+1)/2+jmax)-3)
    complex*16, intent(inout) :: cauchy(number_of_layers,5*(jmax*(jmax+1)/2+jmax)-3), displacement(number_of_layers+1,3*(jmax*(jmax+1)/2+jmax)+1)

    ! Declare workspace matrices and vectors with fixed bounds.
    real*8 :: A(1000,1000), AL(1000,1000)
    real*8 :: yr(3*number_of_layers+1), yi(3*number_of_layers+1), d

    ! Declare integer variables for indexing and dimensioning purposes.
    integer :: i, k, n, m1, m2, np, mp, mpl

    ! Allocate the index array for the LU decomposition.
    integer, allocatable :: INDX(:)

    ! Initialize the real and imaginary parts of the solution vectors.
    yr = 0.0d0
    yi=0.0d0

    ! Set the dimensions and bandwidths for the matrix A.
    n = 3*number_of_layers+1  ! Dimension of the matrix B.
    m1 = 2  ! Number of sub-diagonals in matrix B.
    m2 = 2  ! Number of super-diagonals in matrix B.
    np = 1000  ! Number of rows in the matrix A, which holds the band of matrix B.
    mp = 1000  ! Number of columns in the matrix A.
    mpl = 1000  ! Number of columns in the matrix AL (auxiliary matrix).

    ! Initialize the band matrix A to zero before populating with values from 'matrix'.
    do i = 1, n
        do k = 1, m1 + m2 + 1
            A(i, k) = 0d0
        end do
    end do


    ! Populate the band matrix A from 'matrix'.
    do i = 1, n
        do k = 1, m1 + m2 + 1
            if (i + k - m1 - 1 > 0 .and. i + k - m1 - 1 <= n) then
                A(i, k) = toroidial_matrix(i, i + k - m1 - 1)
            endif
        end do
    end do

    ! Allocate the index array for LU decomposition.
    ALLOCATE(INDX(n))

    ! Perform LU decomposition on the band matrix A.
    call bandec (A, n, m1, m2, np, mp, AL, mpl, INDX, d)

    do i = 1, number_of_layers
        yr(i*3-1) = yr(i*3-1) + dreal(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 6))
        yr(i*3) = yr(i*3) + dreal(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 4))
        yi(i*3-1) = yi(i*3-1) + dimag(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 6))
        yi(i*3) = yi(i*3) + dimag(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 4))
    end do

    ! Solve the system for the real and imaginary parts.
    call banbks(A, n, m1, m2, np, mp, AL, mpl, INDX, yr)  ! Solve for real part.
    call banbks(A, n, m1, m2, np, mp, AL, mpl, INDX, yi)  ! Solve for imaginary part.

    ! Deallocate the index array.
    DEALLOCATE(INDX)

    ! Store the solution into 'displacement' and 'cauchy' arrays.
    do i = 1, number_of_layers
        displacement(i, 3 * ((j * (j + 1)) / 2 + m)) = dcmplx(yr(3 * (i - 1) + 1), yi(3 * (i - 1) + 1))
        cauchy(i, 5 * ((j * (j + 1)) / 2 + m) - 6) = dcmplx(yr(3 * (i - 1) + 2), yi(3 * (i - 1) + 2))
        cauchy(i, 5 * ((j * (j + 1)) / 2 + m) - 4) = dcmplx(yr(3 * (i - 1) + 3), yi(3 * (i - 1) + 3))
    end do

    ! Store the solution for the last layer into 'displacement'.
    displacement(number_of_layers + 1, 3 * ((j * (j + 1)) / 2 + m)) = dcmplx(yr(3 * number_of_layers + 1), yi(3 * number_of_layers + 1))

end subroutine

subroutine solve_system_for_j1(number_of_layers, jmax, j, m, matrix_for_j1, cauchy_integral, cauchy, cauchy_isotropic, displacement)
    implicit none

    integer, intent(in) :: number_of_layers, jmax, j, m
    real*8, intent(in) :: matrix_for_j1(5*number_of_layers+2,5*number_of_layers+2)
    complex*16, intent(in) :: cauchy_integral(number_of_layers,5*(jmax*(jmax+1)/2+jmax)-3)
    complex*16, intent(inout) :: cauchy(number_of_layers,5*(jmax*(jmax+1)/2+jmax)-3), cauchy_isotropic(number_of_layers, jmax + 1, jmax + 1), &
    displacement(number_of_layers+1,3*(jmax*(jmax+1)/2+jmax)+1)

    ! Declare workspace matrices and vectors with fixed bounds.
    real*8 :: A(1000,1000), AL(1000,1000)
    real*8 :: yr(5*number_of_layers+2), yi(5*number_of_layers+2), d

    ! Declare integer variables for indexing and dimensioning purposes.
    integer :: i, k, n, m1, m2, np, mp, mpl

    ! Allocate the index array for the LU decomposition.
    integer, allocatable :: INDX(:)

    ! Initialize the real and imaginary parts of the solution vectors.
    yr = 0.0d0
    yi=0.0d0

    ! Set the dimensions and bandwidths for the matrix A.
    n = 5*number_of_layers+2  ! Dimension of the matrix B.
    m1 = 6  ! Number of sub-diagonals in matrix B.
    m2 = 6  ! Number of super-diagonals in matrix B.
    np = 1000  ! Number of rows in the matrix A, which holds the band of matrix B.
    mp = 1000  ! Number of columns in the matrix A.
    mpl = 1000  ! Number of columns in the matrix AL (auxiliary matrix).

    ! Initialize the band matrix A to zero before populating with values from 'matrix'.
    do i = 1, n
    do k = 1, m1 + m2 + 1
        A(i, k) = 0d0
    end do
    end do


    ! Populate the band matrix A from 'matrix'.
    do i = 1, n
    do k = 1, m1 + m2 + 1
        if (i + k - m1 - 1 > 0 .and. i + k - m1 - 1 <= n) then
            A(i, k) = matrix_for_j1(i, i + k - m1 - 1)
        endif
    end do
    end do

    ! Allocate the index array for LU decomposition.
    ALLOCATE(INDX(n))

    ! Perform LU decomposition on the band matrix A.
    call bandec (A, n, m1, m2, np, mp, AL, mpl, INDX, d)

    do i = 1, number_of_layers
        yr(i*5-1) = yr(i*5-1) + dreal(cauchy_integral(i, 3 * (j * (j + 1) / 2 + m) - 1))
        yr(i*5) = yr(i*5) + dreal(cauchy_integral(i, 3 * (j * (j + 1) / 2 + m) + 1))
        yi(i*5-1) = yi(i*5-1) + dimag(cauchy_integral(i, 3 * (j * (j + 1) / 2 + m) - 1))
        yi(i*5) = yi(i*5) +dimag(cauchy_integral(i, 3 * (j * (j + 1) / 2 + m) + 1))
    end do

    ! Solve the system for the real and imaginary parts.
    call banbks(A, n, m1, m2, np, mp, AL, mpl, INDX, yr)  ! Solve for real part.
    call banbks(A, n, m1, m2, np, mp, AL, mpl, INDX, yi)  ! Solve for imaginary part.

    ! Deallocate the index array.
    DEALLOCATE(INDX)

    ! Store the solution into 'displacement' and 'cauchy' arrays.
    do i = 1, number_of_layers
        displacement(i, 3 * ((j * (j + 1)) / 2 + m) - 1) = dcmplx(yr(5 * (i - 1) + 1), yi(5 * (i - 1) + 1))
        displacement(i, 3 * ((j * (j + 1)) / 2 + m) + 1) = dcmplx(yr(5 * (i - 1) + 2), yi(5 * (i - 1) + 2))
        cauchy_isotropic(i, j+1, m+1) = dcmplx(yr(5 * (i - 1) + 3), yi(5 * (i - 1) + 3))
        cauchy(i, 3 * ((j * (j + 1)) / 2 + m) - 1) = dcmplx(yr(5 * (i - 1) + 4), yi(5 * (i - 1) + 4))
        cauchy(i, 3 * ((j * (j + 1)) / 2 + m) + 1) = dcmplx(yr(5 * (i - 1) + 5), yi(5 * (i - 1) + 5))
    end do

    ! Store the solution for the last layer into 'displacement'.
    displacement(number_of_layers + 1, 3 * ((j * (j + 1)) / 2 + m) - 1) = dcmplx(yr(5 * number_of_layers + 1), yi(5 * number_of_layers + 1))
    displacement(number_of_layers + 1, 3 * ((j * (j + 1)) / 2 + m) + 1) = dcmplx(yr(5 * number_of_layers + 2), yi(5 * number_of_layers + 2))

end subroutine

! This program simulates the deformation of Jupiter's moon Europa
program Europa_simulation
    implicit none   ! Require all variables to be explicitly declared

    ! Constants for the simulation
    integer, parameter :: number_of_layers=101
    integer, parameter :: jmax=2
    integer, parameter :: number_of_time_steps=1000
    real*8, parameter :: delta_t=450

    ! Physical parameters for Europa's ice layer
    real*8, parameter :: radius = 1510000.0d0        ! Radius of Europa in meters
    real*8, parameter :: thickness = 30000.0d0       ! Thickness of ice in meters
    real*8, parameter :: mu = 3.3E9                  ! Shear modulus in Pascal
    real*8, parameter :: eta = 1.0E13                ! Vizkozity in Pascal*s
    real*8, parameter :: surface_g = 1.314           ! Surface gravity in m/s^2
    real*8, parameter :: bottom_g = 1.314            ! Bottom gravity in m/s^2 (same as surface gravity)
    real*8, parameter :: ice_density = 970.0d0       ! Density of ice in kg/m^3
    real*8, parameter :: delta_rho = 30.0d0          ! Density difference in kg/m^3
    real*8, parameter :: angular_speed = 2.0478E-5   ! Angular speed in rad/s
    real*8, parameter :: excentricity = 0.009        ! Eccentricity of Europa's orbit

    real*8, parameter :: delta_r = thickness / (number_of_layers - 1)  ! Radial distance between layers


    ! Allocate fields and matrices for storing calculations
    complex*16 :: displacement(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)
    complex*16 :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    complex*16 :: cauchy_isotropic(number_of_layers, jmax+1, jmax+1)
    complex*16 :: cauchy_integral(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)

    complex*16 :: volume_force(number_of_layers - 1, 2)   ! The volume force calculated as the gradient of the tidal potential, it has only two non-zero components at index 8 and 14
    complex*16 :: bottom_force(4)    ! The force caused by the pressure at the lower boundary, calculated from the tidal potential, it has 4 nontrivial components: 8, 10, 14, 16
    
    real*8 :: matrix(6 * number_of_layers + 2, 6 * number_of_layers + 2)
    real*8 :: matrix_for_j1(5 * number_of_layers + 2, 5 * number_of_layers + 2)
    real*8 :: toroidial_matrix(3 * number_of_layers + 1, 3 * number_of_layers + 1)


    ! Loop variables and arrays for results
    integer ::t,j,i,m,k
    real*8 :: j1

    complex*16 :: radial_displacement(jmax+1, jmax+1) 
    real*8 :: Q_in_time(number_of_time_steps, number_of_layers)
    real*8 :: Qcum(number_of_layers)
    real*8 :: Q, fac
    real*8 :: Q_average(number_of_time_steps, number_of_layers)
    real*8 :: Q_total(number_of_time_steps)
    complex*16 :: t_comp

    logical :: test
    logical :: write_deformation_data

    character(len=50) :: folder_name
    character(len=8)  :: date
    character(len=10) :: time
    integer, dimension(8) :: values
    integer :: status

    test = .TRUE.
    write_deformation_data = .FALSE.

    if (write_deformation_data) then

        ! Get the current date and time
        call date_and_time(date, time, values=values)

        ! Create timestamp in the format YYYYMMDD_HHMMSS
        write(folder_name, '(I4.4,I2.2,I2.2,"_",I2.2,I2.2,I2.2)') &
            values(1), values(2), values(3), values(5), values(6), values(7)

        ! Append folder description
        folder_name = trim(folder_name) // "_simulation_data"

        ! Create the directory
        call execute_command_line("mkdir -p " // trim(folder_name), exitstat=status)
        if (status /= 0) then
            print *, "Error: Unable to create directory!"
            stop
        end if

        print *, "Folder created: ", trim(folder_name)

        open(unit=10, file=trim(folder_name)//"/displacement101.txt",status="replace", action="write")
        close(10)

    end if

    ! Initialize arrays with zero values
    displacement = 0.0d0
    cauchy = 0.0d0
    cauchy_isotropic = 0.0d0
    cauchy_integral = 0.0d0

    bottom_force = 0.0d0
    volume_force = 0.0d0

    matrix = 0.0d0
    toroidial_matrix = 0.0d0
    matrix_for_j1 = 0.0d0

    !!!!!!!
    radial_displacement = 0.0d0
    Q_in_time = 0.0d0
    Qcum = 0.0d0
    Q = 0.0d0
    t_comp = 0
    Q_average = 0.0d0
    Q_total = 0.0d0


    ! Main simulation loop over time steps, careful to number_of_time_steps - 1
    do t=0,number_of_time_steps-1

        ! Calculate the forces at this time step
        call calculate_forces(number_of_layers, t, volume_force, bottom_force, radius, delta_r, angular_speed, delta_t, ice_density, delta_rho, excentricity)

        call update_cauchy_integral(jmax, number_of_layers, cauchy_integral, cauchy, delta_t, eta, mu)

        j=1

        do m=0,1

            call fill_matrix_for_j1(number_of_layers, matrix_for_j1, mu, ice_density, radius, delta_r, delta_rho, surface_g, bottom_g)
            
            call solve_system_for_j1(number_of_layers, jmax, j, m, matrix_for_j1, cauchy_integral, cauchy, cauchy_isotropic, displacement)

        end do

        ! Loop over each harmonic degree 'j'
        do j=2,jmax

            j1 = real(j)

            call fill_matrix(number_of_layers, matrix, j1, mu, ice_density, radius, delta_r, delta_rho, surface_g, bottom_g)

            call fill_toroidial_matrix(number_of_layers, toroidial_matrix, j1, mu, radius, delta_r)

            ! Loop over each order 'm' for the current degree 'j'
            do m=0,j

                ! Solve the linear system for the current harmonic order
                call solve_linear_system(number_of_layers, jmax, j, m, matrix, volume_force, bottom_force, cauchy_integral, cauchy, cauchy_isotropic, displacement)

                call solve_toroidial_system(number_of_layers, jmax, j, m, toroidial_matrix, cauchy_integral, cauchy, displacement)

            end do
        end do

        if (test) then
            
            print *
            print '(A, I0)', "Running tests at time step ", t

            call test_equations_solution(number_of_layers, jmax, radius, delta_r, mu, ice_density, surface_g, bottom_g, delta_rho, volume_force, bottom_force, cauchy_integral, cauchy, cauchy_isotropic, displacement)

        end if

        if (write_deformation_data) then

            call calculate_radial_displacement_at_layer(100, jmax, number_of_layers, displacement, radial_displacement)

            open(unit=10, file=trim(folder_name)//"/displacement101.txt", status="old", action="write", position="append")
            ! Write time step
            write(10, *) t * delta_t
            ! Write harmonic coefficients (radial displacement)
            do j = 0, jmax
                do m = 0, j
                    write(10,*) j, m, real(radial_displacement(j+1, m+1)), aimag(radial_displacement(j+1, m+1))
                end do
            end do
            close(10)
        
            print *, "Radial displacement data written to: ", trim(folder_name)//"/displacement101.txt"
        
            
        end if

        !!!!!!!!!!!!!
        ! Pocitani zahrivani
        do i=1, number_of_layers

            Qcum=0
            do j=2, jmax

                do m=0,j
                    fac=2d0
                    if(m.eq.0) fac=1d0
                    t_comp = cauchy(i, 5 * (j * (j + 1) / 2 + m) - 3)
                    Qcum(i) = Qcum(i) + fac*(real(t_comp)**2+aimag(t_comp)**2)
                    t_comp = cauchy(i, 5 * (j * (j + 1) / 2 + m) - 5)
                    Qcum(i) = Qcum(i) + fac*(real(t_comp)**2+aimag(t_comp)**2)
                    t_comp = cauchy(i, 5 * (j * (j + 1) / 2 + m) - 7)
                    Qcum(i) = Qcum(i) + fac*(real(t_comp)**2+aimag(t_comp)**2)
                end do

            end do

            Q_in_time(t+1, i) = Qcum(i)/(2*eta)
        end do
        if (t>20) then
            do i=21,t
                do m=1, number_of_layers
                    Q_average(t, m) = Q_average(t, m) + Q_in_time(i, m)
                end do
            end do
        do m=1, number_of_layers
            Q_average(t, m) = Q_average(t, m)/(t-20)
            !print*, Q_average(t, m)
        end do
        do m=1, number_of_layers -1
            Q_total(t) = Q_total(t) + ((Q_average(t,m)+Q_average(t,m+1))/2.0)*(radius+m*delta_r/2.0)**2*delta_r
        end do
        print*, t, Q_total(t)


        endif

    end do

    open(4, file = 'Q.dat')
        ! Write the time and displacement data to file
        do i=21, number_of_time_steps-1
            write(4,*) (i-1)*(delta_t)*angular_speed/(2.0*acos(-1.0d0)), Q_total(i)/1e9
        end do
    close(4)


end program Europa_simulation
