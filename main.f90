
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

  !matrix(2,1)=1/2.0
  !matrix(2,7)=1/2.0

  !matrix(3,2)=1/2.0
  !matrix(3,8)=1/2.0

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

subroutine calculate_forces(number_of_layers, t, jmax, volume_force, bottom_force, &
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
  integer, intent(in) :: number_of_layers, t, jmax
  complex*16, intent(out) :: volume_force(number_of_layers-1, 3*(jmax*(jmax+1)/2+jmax)+1),&
                             bottom_force(3*(jmax*(jmax+1)/2+jmax)+1)
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
      volume_force(i,8)=cmplx(-ksi*sqrt(18*4.D0*datan(1.D0))*dcos(angular_speed*t_*delta_t),0)

      real_=ksi*sqrt(27*4.D0*datan(1.D0))*dcos(angular_speed*t_*delta_t)
      imag_=-ksi*sqrt(48*4.D0*datan(1.D0))*dsin(angular_speed*t_*delta_t)!<
      volume_force(i,14)=cmplx(real_,imag_)
  end do

  ! Compute a modified constant factor for bottom force calculation
  ksi = (ice_density + delta_rho) * angular_speed**2 * excentricity * radius**2

  ! Compute the bottom forces at specific indices
  bottom_force(8)=cmplx(-ksi*sqrt((18*4.D0*datan(1.D0))/(25))*dcos(angular_speed*t_*delta_t),0)
  bottom_force(10)=cmplx(ksi*sqrt((27*4.D0*datan(1.D0))/(25))*dcos(angular_speed*t_*delta_t),0)

  real_=ksi*sqrt((27*4.D0*datan(1.D0))/(25))*dcos(angular_speed*t_*delta_t)
  imag_=-ksi*sqrt((48*4.D0*datan(1.D0))/(25))*dsin(angular_speed*t_*delta_t)
  bottom_force(14)=cmplx(real_,imag_)

  real_=-ksi*sqrt((81*4.D0*datan(1.D0))/(50))*dcos(angular_speed*t_*delta_t)
  imag_=ksi*sqrt((72*4.D0*datan(1.D0))/(25))*dsin(angular_speed*t_*delta_t)
  bottom_force(16)=cmplx(real_,imag_)

end subroutine

subroutine solve_linear_system(number_of_layers, jmax, j, m, matrix,&
         volume_force, bottom_force, cauchy_integral, cauchy, posunuti)
  implicit none

  integer, intent(in) :: number_of_layers, jmax, j, m
  real*8, intent(in) :: matrix(6*number_of_layers+2,6*number_of_layers+2)
  complex*16, intent(in) :: volume_force(number_of_layers-1,3*(jmax*(jmax+1)/2+jmax)+1), &
                        bottom_force(3*(jmax*(jmax+1)/2+jmax)+1), &
                        cauchy_integral(number_of_layers,5*(jmax*(jmax+1)/2+jmax)-3)
  complex*16, intent(inout) :: cauchy(number_of_layers,5*(jmax*(jmax+1)/2+jmax)-3), &
                         posunuti(number_of_layers+1,3*(jmax*(jmax+1)/2+jmax)+1)

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
  do i = 1, number_of_layers - 1
      yr(6 * i + 2) = -dreal(volume_force(i, 3 * ((j * (j + 1)) / 2 + m) - 1))
      yi(6 * i + 2) = -dimag(volume_force(i, 3 * ((j * (j + 1)) / 2 + m) - 1))
      yr(6 * i + 3) = -dreal(volume_force(i, 3 * ((j * (j + 1)) / 2 + m) + 1))
      yi(6 * i + 3) = -dimag(volume_force(i, 3 * ((j * (j + 1)) / 2 + m) + 1))
  end do

  ! Set the first two elements from 'bottom_force'.
  yr(2) = dreal(bottom_force(3 * ((j * (j + 1)) / 2 + m) - 1))
  yi(2) = dimag(bottom_force(3 * ((j * (j + 1)) / 2 + m) - 1))
  yr(3) = dreal(bottom_force(3 * ((j * (j + 1)) / 2 + m) + 1))
  yi(3) = dimag(bottom_force(3 * ((j * (j + 1)) / 2 + m) + 1))

  do i = 1, number_of_layers
      yr(i*6-2) = yr(i*6-2) + dreal(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 7))
      yr(i*6-1) = yr(i*6-1) + dreal(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 5))
      yr(i*6) = yr(i*6) + dreal(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 3))
      yi(i*6-2) = yi(i*6-2) + dimag(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 7))
      yi(i*6-1) = yi(i*6-1) + dimag(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 5))
      yi(i*6) = yi(i*6) +dimag(cauchy_integral(i, 5 * (j * (j + 1) / 2 + m) - 3))
  end do

      !if (j==2) then
       !   if (m==0) then
        !      open(6, file = 'rhs.dat')
        !!      do i=1, 6*number_of_layers+2
        !          write(6,*) yr(i)
         !     end do
         !     close(6)

         !     open(7, file= 'matrix.dat')
         !         do i = 1, 6*number_of_layers+2
        !          write(7, *) (matrix(i, k), k=1, 6*number_of_layers+2)
        !          end do
        !      close(7)
        !  end if
     ! end if



  ! Solve the system for the real and imaginary parts.
  call banbks(A, n, m1, m2, np, mp, AL, mpl, INDX, yr)  ! Solve for real part.
  call banbks(A, n, m1, m2, np, mp, AL, mpl, INDX, yi)  ! Solve for imaginary part.

  ! Deallocate the index array.
  DEALLOCATE(INDX)

  ! Store the solution into 'posunuti' and 'cauchy' arrays.
  do i = 1, number_of_layers
      posunuti(i, 3 * ((j * (j + 1)) / 2 + m) - 1) = cmplx(yr(6 * (i - 1) + 1), yi(6 * (i - 1) + 1))
      posunuti(i, 3 * ((j * (j + 1)) / 2 + m) + 1) = cmplx(yr(6 * (i - 1) + 2), yi(6 * (i - 1) + 2))
      cauchy(i, 5 * ((j * (j + 1)) / 2 + m) - 7) = cmplx(yr(6 * (i - 1) + 4), yi(6 * (i - 1) + 4))
      cauchy(i, 5 * ((j * (j + 1)) / 2 + m) - 5) = cmplx(yr(6 * (i - 1) + 5), yi(6 * (i - 1) + 5))
      cauchy(i, 5 * ((j * (j + 1)) / 2 + m) - 3) = cmplx(yr(6 * (i - 1) + 6), yi(6 * (i - 1) + 6))
  end do

  ! Store the solution for the last layer into 'posunuti'.
  posunuti(number_of_layers + 1, 3 * ((j * (j + 1)) / 2 + m) - 1) = cmplx(yr(6 * number_of_layers + 1),&
                                                                      yi(6 * number_of_layers + 1))
  posunuti(number_of_layers + 1, 3 * ((j * (j + 1)) / 2 + m) + 1) = cmplx(yr(6 * number_of_layers + 2),&
                                                                      yi(6 * number_of_layers + 2))


!!!!!!
! Tady ukladam vysledky do souboru
  if (j==2) then
      if (m==0) then
      open(4, file = 'posunuti.dat')
          do i=1, number_of_layers+1
             write(4,*) 76000 + 800*(i-3.0/2.0), real(posunuti(i, 8)), real(posunuti(i, 10))
          end do
         close(4)
          open(5, file = 'napeti.dat')
          do i=1, number_of_layers
              write(5,*) 76000 + 800*(i-1), yr(6 * (i - 1) + 3), yr(6 * (i - 1) + 4), yr(6 * (i - 1) + 5), yr(6 * (i - 1) + 6)
          end do
          close(5)
     end if
  end if

  !!!!!!!!

end subroutine

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

subroutine calculate_and_print_surface_radial_displacement(number_of_layers, jmax, t, number_of_time_steps, uu2, uu0, delta_t, posunuti)
      integer number_of_layers, t, jmax, number_of_time_steps
      real*8 delta_t
      complex*16 posunuti(number_of_layers+1,3*(jmax*(jmax+1)/2+jmax)+1), uu2(number_of_time_steps), uu0(number_of_time_steps),&
       koef1, koef2
      integer j, m, th, ph, i
      real*8 fac, x, thr, phr, j1, p(1000), xm
      complex*8 :: koeficient(30)


      th=90 ! zemepisne souradnice bodu, kde chceme pocitat hodnotu
      ph=0
      thr=th*datan(1d0)/45d0 ! prepocet na radiany
      phr=ph*datan(1d0)/45d0

      x=dcos(thr) ! k vypoctu PLF potrebujeme cos(th)

  !    vypocet PLF v bode cos(th); hodnoty ulozeny v jednorozmernem poli p
      call DPNM(x,jmax,p)

  !   vypocet SH funkci a soucet rady v bode th,ph
      do j=2,jmax
          do m=0,j
          j1=real(j)
          xm=dble(m)
          fac=2d0
          if(m.eq.0) fac=1d0
          i=j*(j+1)/2+m+1 ! sdruzeny index
          koef1=sqrt((j1)/(2*j1+1))*(1/2.0)*(posunuti(number_of_layers,3*(j*(j+1)/2+m)-1)+ &
                      posunuti(number_of_layers+1,3*(j*(j+1)/2+m)-1))
          koef2=sqrt((j1+1)/(2*j1+1))*(1/2.0)*(posunuti(number_of_layers,3*(j*(j+1)/2+m)+1)+ &
                      posunuti(number_of_layers+1,3*(j*(j+1)/2+m)+1))
          !koef1=sqrt((j1)/(2*j1+1))*(1/2.0)*(posunuti(1,3*(j*(j+1)/2+m)-1)+ &
           !           posunuti(2,3*(j*(j+1)/2+m)-1))
          !koef2=sqrt((j1+1)/(2*j1+1))*(1/2.0)*(posunuti(1,3*(j*(j+1)/2+m)+1)+ &
           !           posunuti(2,3*(j*(j+1)/2+m)+1))
          !koeficient(i)=koef1-koef2
          !uu(t+1)=uu(t+1)+dreal(dcmplx(p(i),0)*koeficient(i))*fac
          if (m.eq.0) uu0(t+1) = koef1-koef2
          if (m.eq.2) uu2(t+1) = koef1-koef2
          enddo

      end do

      print*, t*(delta_t)/3600, uu0(t+1), uu2(t+1)
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
  real*8, parameter :: radius = 1501000.0d0        ! Radius of Europa in meters
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



  ! Allocate matrices for storing calculations
  real*8 :: matrix(6 * number_of_layers + 2, 6 * number_of_layers + 2)
  complex*16 :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
  complex*16 :: cauchy_integral(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
  complex*16 :: posunuti(number_of_layers + 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)
  complex*16 :: bottom_force(3 * (jmax * (jmax + 1) / 2 + jmax) + 1)
  complex*16 :: volume_force(number_of_layers - 1, 3 * (jmax * (jmax + 1) / 2 + jmax) + 1)

  ! Loop variables and arrays for results
  integer ::t,j,i,m,k
  real*8 :: j1
  complex*16 :: uu2(number_of_time_steps)
  complex*16 :: uu0(number_of_time_steps)
  real*8 :: Q_in_time(number_of_time_steps, number_of_layers)
  real*8 :: Qcum(number_of_layers)
  real*8 :: Q, fac
  real*8 :: Q_average(number_of_time_steps, number_of_layers)
  real*8 :: Q_total(number_of_time_steps)
  complex*16 :: t_comp

  !!!!!!!!!!!!!
  real*8 :: koef1, koef2, ksi
  complex*16 :: volume_force_8(number_of_time_steps)
  complex*16 :: volume_force_14(number_of_time_steps)
  complex*16 :: bottom_force_0(number_of_time_steps)
  complex*16 :: real_bottom_force(number_of_time_steps)
  !!!!!!!!!!!!!

  ! Initialize arrays with zero values
  matrix = 0.0d0
  cauchy = 0.0d0
  cauchy_integral = 0.0d0
  posunuti = 0.0d0
  bottom_force = 0.0d0
  volume_force = 0.0d0

  !!!!!!!
  uu2 = 0.0d0
  uu0 = 0.0d0
  Q_in_time = 0.0d0
  Qcum = 0.0d0
  Q = 0.0d0
  t_comp = 0
  Q_average = 0.0d0
  Q_total = 0.0d0

  !!!!!!!!
  ! Pro testovani
  volume_force_14 = 0.0d0
  volume_force_8 = 0.0d0
  bottom_force_0 = 0.0d0
  real_bottom_force = 0.0d0
  !!!!!!!!

  ! Main simulation loop over time steps, careful to number_of_time_steps - 1
  do t=0,number_of_time_steps-1
      ! Calculate the forces at this time step
      call calculate_forces(number_of_layers, t, jmax, volume_force, bottom_force, radius, delta_r, angular_speed, delta_t, ice_density, delta_rho, excentricity)


      !!!!!!!!!!!!!
      ! Testovaci kod
      volume_force_8(t+1) = (1/2.0)*(volume_force(73,8)+volume_force(74,8))
      volume_force_14(t+1) = (1/2.0)*(volume_force(73,14)+volume_force(74,14))

      koef1=sqrt((2.0)/(5.0))*bottom_force(8)
      koef2=sqrt((3.0)/(5.0))*bottom_force(10)
      bottom_force_0(t+1)= koef1-koef2
      ksi = (ice_density + delta_rho) * angular_speed**2 * excentricity * radius**2
      real_bottom_force(t+1) = cmplx(-ksi*sqrt((18.D0*4*datan(1.D0))/(10))*dcos(angular_speed*t*delta_t),0)
      !!!!!!!!!!!!!

      !print*, bottom_force(8), bottom_force(10)



      ! Loop over each harmonic degree 'j'
      do j=2,jmax
          j1 = real(j)
          call fill_matrix(number_of_layers, matrix, j1, mu, ice_density, radius, delta_r, delta_rho, surface_g, bottom_g)

          ! Loop over each order 'm' for the current degree 'j'
          do m=0,j

              do k=1, number_of_layers
                  cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 3) = cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 3) -&
                   (1.0/eta) * cauchy(k, 5 * (j * (j + 1) / 2 + m) - 3) * delta_t
                  cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 5) = cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 5)-&
                   (1.0/eta) * cauchy(k, 5 * (j * (j + 1) / 2 + m) - 5) * delta_t
                  cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 7) = cauchy_integral(k, 5 * (j * (j + 1) / 2 + m) - 7)-&
                   (1.0/eta) * cauchy(k, 5 * (j * (j + 1) / 2 + m) - 7) * delta_t
              end do

              ! Prozatim bez objemove sily a ciste elasticke teleso
              ! volume_force = 0.0d0
              ! cauchy_integral = 0.0d0

              ! Solve the linear system for the current harmonic order
              call solve_linear_system(number_of_layers, jmax, j, m, matrix, volume_force,&
                   bottom_force, cauchy_integral, cauchy, posunuti)
          end do
      end do

      ! Calculate and print the radial displacement at the surface
      call calculate_and_print_surface_radial_displacement(number_of_layers, jmax, t, number_of_time_steps, uu2, uu0, delta_t, posunuti)

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



  open(1, file = 'file.dat')
      ! Write the time and displacement data to file
      do i=1, number_of_time_steps
          write(1,*) (i-1)*(delta_t)*angular_speed/(2.0*acos(-1.0d0)), real(uu0(i)), real(uu2(i)), aimag(uu2(i))
      end do
  close(1)

  !!!!!!!!!!!!!
  open(2, file = 'volume_force.dat')
  ! Write the time and displacement data to file
  do i=1, number_of_time_steps
      write(2,*) volume_force_8(i), volume_force_14(i)
  end do
  close(2)
  open(3, file = 'bottom_force.dat')
  ! Write the time and displacement data to file
  do i=1, number_of_time_steps
      write(3,*) bottom_force_0(i)
  end do
  close(3)
  open(3, file = 'real_bottom_force.dat')
  ! Write the time and displacement data to file
  do i=1, number_of_time_steps
      write(3,*) real_bottom_force(i)
  end do
  close(3)
  !!!!!!!
      open(4, file = 'Q.dat')
      ! Write the time and displacement data to file
      do i=21, number_of_time_steps-1
          write(4,*) (i-1)*(delta_t)*angular_speed/(2.0*acos(-1.0d0)), Q_total(i)/1e9
      end do
  close(4)


end program Europa_simulation
