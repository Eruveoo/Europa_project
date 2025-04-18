!C ----------------------------------------------------------------
SUBROUTINE VCSUM(NP,AJM,BJM,NR,CJM)
!C ----------------------------------------------------------------
    IMPLICIT REAL*8(D-H,O-Z),COMPLEX*16(A-C)
    DIMENSION AJM(*),BJM(*),CJM(*)
    DIMENSION CSUMA(4,301),COMP(4)
    DIMENSION RE(2049),QIM(2049)
    COMMON /DD0/ NGQ,NFOURC,ROOTS(200),WGHTS(200),PNMARR(200,7000),CEMARR(10192,200)

    NFOUR=NFOURC
    NROOT=NGQ/2
    NP1=NP+1
    NR1=NR+1
    DNFR=DFLOAT(NFOUR)
    PI2=8.D00*DATAN(1.D00)
    PI4=PI2+PI2

    DO 20 IR=1,NROOT
    POMWT=PI4*WGHTS(IR)
    DO 6 MS=1,NP1
    DO 4 K=1,4
4 CSUMA(K,MS)=(0.D00,0.D00)
    DO 6 JS=MS,NP1
    JM=(JS-1)*JS/2+MS
    DO 5 K=1,4
    CPOM=AJM(JM)
    IF(K.GT.2) CPOM=BJM(JM)
    ZK=(-1.D00)**(JS+MS)
    ZK=ZK**(K+1)
5 CSUMA(K,MS)=CSUMA(K,MS)+CPOM*PNMARR(IR,JM)*ZK
6 CONTINUE
!C     EVALUATE THE PRODUCTS AROUND THE LATITUDE CIRCLE
    DO 11 JFI=1,NFOUR
    DO 7 K=1,4
7 COMP(K)=(0.D00,0.D00)
    DO 9 M=1,NP
    MS=M+1
    DO 8 K=1,4
8 COMP(K)=COMP(K)+CSUMA(K,MS)*CEMARR(JFI,M)
9 CONTINUE
    DO 10 K=1,4
10 COMP(K)=COMP(K)+DCONJG(COMP(K))+CSUMA(K,1)
    FABP=DBLE(COMP(1)*COMP(3))
    FABN=DBLE(COMP(2)*COMP(4))
    RE(JFI)=FABP
11 QIM(JFI)=FABN
!C     PERFORM THE FFT OF TWO REAL SERIES
    CALL DFFTAR(RE,QIM,NFOUR,0)
    RE(NFOUR+1)=RE(1)
    QIM(NFOUR+1)=QIM(1)
    DO 12 MS=1,NR1
    NFMS=NFOUR+2-MS
    POM1=(RE(MS)+RE(NFMS))/2.D00/DNFR
    POM2=(QIM(MS)-QIM(NFMS))/2.D00/DNFR
    CAP=DCMPLX(POM1,POM2)
    POM1=(QIM(MS)+QIM(NFMS))/2.D00/DNFR
    POM2=(-RE(MS)+RE(NFMS))/2.D00/DNFR
    CAN=DCMPLX(POM1,POM2)
!C     FINALLY, EVALUATE THE IR-TH TERM OF GAUSS-LEGENDRE QUADRATURE
    DO 12 JS=MS,NR1
    JM=(JS-1)*JS/2+MS
    ZN=(-1.D00)**(JS+MS)
    IF(IR.EQ.1) CJM(JM)=(0.D00,0.D00)
12 CJM(JM)=CJM(JM)+POMWT*PNMARR(IR,JM)*(CAP+CAN*ZN)
20 CONTINUE
    RETURN
    END

!........................................................................
subroutine CLEB1(j1,m1,j2,m2,j,m,cc)

!     Clebch-Gordanovy koeficienty pro j2=1	C_(j1,m1,j2,m2)^(j,m)
!........................................................................
    implicit real*8 (a-h,o-z)
    cc=0.d00
    if(j2.ne.1) return
    if(iabs(j1-j).gt.1.or.(j1+j).eq.0) return
    if(iabs(m2).gt.1.or.iabs(m1).gt.j1) return
    if(m.ne.(m1+m2)) return
    c=dble(j)
    g=dble(m)
    zn=1d00
    if(m2) 10,20,30
10 if(j1-j) 12,14,16
12 cc=(c-g-1d00)*(c-g)/(c+c-1d00)/(c+c)
    goto 100
14 cc=(c+g+1d00)*(c-g)/(c+1d00)/(c+c)
    goto 100
16 cc=(c+g+2d00)*(c+g+1d00)/(c+c+2d00)/(c+c+3d00)
    goto 100
20 if(j1-j) 22,24,26
22 cc=(c+g)*(c-g)/(c+c-1d00)/c
    goto 100
24 cc=g/dsqrt(c*(c+1d00))
    return
26 cc=(c+g+1d00)*(c-g+1d00)/(c+1d00)/(c+c+3d00)
    zn=-1d00
    goto 100
30 if(j1-j) 32,34,36
32 cc=(c+g-1d00)*(c+g)/(c+c-1d00)/(c+c)
    goto 100
34 cc=(c+g)*(c-g+1d00)/(c+1d00)/(c+c)
    zn=-1d00
    goto 100
36 cc=(c-g+1d00)*(c-g+2d00)/(c+c+2d00)/(c+c+3d00)
100 cc=zn*dsqrt(cc)
    return
    end

!.....................................................................
    subroutine SIXJ1(j1,j2,j,l2,l,c6)
!
!     6-j symboly:  {j1,j2,j
!                     1,l2,l}
!.....................................................................
    implicit real*8 (a-h,o-z)
    c6=0.d00
    if(l2.gt.(j+1).or.iabs(j-1).gt.l2) return
    if(j.gt.(j1+j2).or.iabs(j1-j2).gt.j) return
    if(l.gt.(j2+1).or.iabs(j2-1).gt.l) return
    a=dfloat(j1)
    b=dfloat(j2)
    c=dfloat(j)
    s=a+b+c
    zn=(-1.d00)**(j1+j2+j)
    if(j-l2) 10,20,30
10 if(j2-l) 12,14,16
12 hor=(s+2)*(s+3)*(s-a-a+1)*(s-a-a+2)
    dol=(b+b+1)*(b+1)*(b+b+3)*(c+c+1)*(c+1)*(c+c+3)
    goto 100
14 hor=(s+2)*(s-c-c)*(s-b-b+1)*(s-a-a+1)
    dol=b*(b+b+1)*(b+1)*(c+c+1)*(c+1)*(c+c+3)
    zn=-zn
    goto 100
16 hor=(s-c-c-1)*(s-c-c)*(s-b-b+1)*(s-b-b+2)
    dol=(b+b-1)*b*(b+b+1)*(c+c+1)*(c+1)*(c+c+3)
    goto 100
20 if(j2-l) 22,24,26
22 hor=(s+2)*(s-c-c+1)*(s-b-b)*(s-a-a+1)
    dol=(b+b+1)*(b+1)*(b+b+3)*c*(c+c+1)*(c+1)
    zn=-zn
    goto 100
24 hor=-a*(a+1)+b*(b+1)+c*(c+1)
    dol=b*(b+b+1)*(b+1)*c*(c+c+1)*(c+1)
    c6=-zn*hor/dsqrt(dol)/2.d00
    return
26 hor=(s+1)*(s-c-c)*(s-b-b+1)*(s-a-a)
    dol=(b+b-1)*b*(b+b+1)*c*(c+c+1)*(c+1)
    goto 100
30 if(j2-l) 32,34,36
32 hor=(s-c-c+1)*(s-c-c+2)*(s-b-b-1)*(s-b-b)
    dol=(b+b+1)*(b+1)*(b+b+3)*(c+c-1)*c*(c+c+1)
    goto 100
34 hor=(s+1)*(s-c-c+1)*(s-b-b)*(s-a-a)
    dol=b*(b+b+1)*(b+1)*(c+c-1)*c*(c+c+1)
    goto 100
36 hor=s*(s+1)*(s-a-a-1)*(s-a-a)
    dol=(b+b-1)*b*(b+b+1)*(c+c-1)*c*(c+c+1)
100 c6=zn*dsqrt(hor/dol)/2.d00 
    return
end

SUBROUTINE shz_tens(jmax, tens, xx, xy, xz, yx, yy, yz, zx, zy, zz)        ! pocita kartezske zlozky tenzoru 3x3 zadaneho v sfer. harmonikach
    implicit none
    INTEGER :: jmax
    COMPLEX(8) :: tens((jmax+3)*(jmax-2)*9/2+17+jmax*9+4+2+2+1), xx((jmax + 3) * (jmax + 4) / 2), xy((jmax + 3) * (jmax + 4) / 2), xz((jmax + 3) * (jmax + 4) / 2), yx((jmax + 3) * (jmax + 4) / 2), yy((jmax + 3) * (jmax + 4) / 2), yz((jmax + 3) * (jmax + 4) / 2), zx((jmax + 3) * (jmax + 4) / 2), zy((jmax + 3) * (jmax + 4) / 2), zz((jmax + 3) * (jmax + 4) / 2), mkonst
    COMPLEX(8), ALLOCATABLE :: xvec(:), yvec(:), zvec(:)
    REAL(8) :: cg, nj
    INTEGER :: indmax1,i, j1, m1, l1, k1, j, m, l, j1max, jnowmin, jnowmax, indvec, INDV, INDS, INDT
    REAL(8), PARAMETER :: PI=3.141592653589793238462643383279502884197

    xx=0    
    xy=0
    xz=0
    yx=0
    yy=0
    yz=0
    zx=0
    zy=0
    zz=0

    indmax1=size(tens)
    j=-1
    do while(.true.)                        !urcenie hodnoty jmax
            j=j+1
            i=INDT(j,j,j+2,2)
            if (i==indmax1) exit
    enddo

    j1max=j

    if (.not.(INDS(j1max+2,j1max+2)==size(xx).and.size(xx)==size(xy).and.size(xy)==size(xz).and.&                !kontrola velkosti vstupnych poli
            INDS(j1max+2,j1max+2)==size(yx).and.size(yx)==size(yy).and.size(yy)==size(yz).and.&
            INDS(j1max+2,j1max+2)==size(zx).and.size(zx)==size(zy).and.size(zy)==size(zz))) then
            write(*,*)'Problem s dimenziami vstupnych poli. - tens'
            return
    endif

    ALLOCATE (xvec(INDV(j1max+1,j1max+1,j1max+2)),yvec(INDV(j1max+1,j1max+1,j1max+2)),zvec(INDV(j1max+1,j1max+1,j1max+2))) 

    !nasobenie jednotk. vektormi 1.-krat

    xvec=0
    yvec=0
    zvec=0

    do j1=0,j1max

            do m1=0,j1
                    do k1=0,2
                            do l1=abs(j1-k1),j1+k1
                                    if (tens(INDT(j1,m1,l1,k1))==0) cycle
                                    l=l1                                                !z vlastnosti c-g koef.
                                    
                                    jnowmin=max(abs(j1-1),abs(l-1))
                                    jnowmax=min(j1+1,l+1)

                                    do j=jnowmin,jnowmax
                                            call ninej0m(j1,j,1,k1,1,l1,l,nj)        !podla Matasa a po 2 permutaciach
                                            mkonst=(-1)**(k1+1)*dSQRT((2*j1+1.0)*3.D0*(2*l1+1)*(2*k1+1.0))/dSQRT(4*pi)
                                            mkonst=mkonst*tens(INDT(j1,m1,l1,k1))*nj
                                                    indvec=INDV(j,0,l)
                                                    call CLEB1(j1,m1,1,1,j,0,cg)
                                                    xvec(indvec)=xvec(INDV(j,0,l))-mkonst*cg
                                                    yvec(indvec)=yvec(INDV(j,0,l))+(0.D0,1.D0)*mkonst*cg
                                                    call CLEB1(j1,m1,1,-1,j,0,cg)
                                                    xvec(indvec)=xvec(INDV(j,0,l))+mkonst*cg
                                                    yvec(indvec)=yvec(INDV(j,0,l))+(0.D0,1.D0)*mkonst*cg
                                                    call CLEB1(j1,m1,1,0,j,0,cg)
                                                    zvec(indvec)=zvec(INDV(j,0,l))+mkonst*cg
                                            do m=1,j
                                                    indvec=INDV(j,m,l)
                                                    call CLEB1(j1,m1,1,1,j,m,cg)
                                                    xvec(indvec)=xvec(INDV(j,m,l))-mkonst*cg
                                                    yvec(indvec)=yvec(INDV(j,m,l))+(0.D0,1.D0)*mkonst*cg
                                                    call CLEB1(j1,m1,1,1,j,-m,cg)
                                                    xvec(indvec)=xvec(INDV(j,m,l))-(-1)**(j+m+l+1)*CONJG(mkonst*cg)
                                                    yvec(indvec)=yvec(INDV(j,m,l))+(-1)**(j+m+l+1)*CONJG((0.D0,1.D0)*mkonst*cg)
                                                    call CLEB1(j1,m1,1,-1,j,m,cg)
                                                    xvec(indvec)=xvec(INDV(j,m,l))+mkonst*cg
                                                    yvec(indvec)=yvec(INDV(j,m,l))+(0.D0,1.D0)*mkonst*cg
                                                    call CLEB1(j1,m1,1,-1,j,-m,cg)
                                                    xvec(indvec)=xvec(INDV(j,m,l))+(-1)**(j+m+l+1)*CONJG(mkonst*cg)
                                                    yvec(indvec)=yvec(INDV(j,m,l))+(-1)**(j+m+l+1)*CONJG((0.D0,1.D0)*mkonst*cg)
                                                    call CLEB1(j1,m1,1,0,j,m,cg)
                                                    zvec(indvec)=zvec(INDV(j,m,l))+mkonst*cg
                                                    call CLEB1(j1,m1,1,0,j,-m,cg)
                                                    zvec(indvec)=zvec(INDV(j,m,l))+(-1)**(j+m+l+1)*CONJG((0.D0,1.D0)*mkonst*cg)
                                            enddo
                                    enddo
                            enddo
                    enddo
            enddo
    enddo

    xvec=xvec*dSQRT(2*pi)
    yvec=yvec*dSQRT(2*pi)
    zvec=zvec*2*dSQRT(pi)

    !nasobenie jednotkovymi vektormi druhykrat resp. vypocet zloziek vektorov

    call shz_vec(j1max, xvec,xx,xy,xz)        
    call shz_vec(j1max, yvec,yx,yy,yz)
    call shz_vec(j1max, zvec,zx,zy,zz)

END SUBROUTINE shz_tens


SUBROUTINE shz_vec(jmax, vec,x,y,z)                !vypocita kartezske zlozky vektoru zadaneho pomocou vek. sfer. harm.
    implicit none

    INTEGER :: jmax
    COMPLEX(8) :: vec(3 * ((jmax + 1) * (jmax + 2) / 2 + jmax+1) + 1), x((jmax + 3) * (jmax + 4) / 2), y((jmax + 3) * (jmax + 4) / 2), z((jmax + 3) * (jmax + 4) / 2), mkonst
    REAL(8) :: cg, sj
    INTEGER :: ind1, indmax1, i, j1, m1, l1, j, m, j1max, indscal, INDV, INDS
    REAL(8), PARAMETER :: PI=3.141592653589793238462643383279502884197


    indmax1=size(vec)
    j=-1
    do while(.true.)                                !urcenie hodnoty jmax
            j=j+1
            i=(j*(j+3)/2)*3+1
            if (i==indmax1) exit
    enddo

    j1max=j
    if (.not.((j1max+3)*(j1max+2)/2==(size(x)).and.(size(x)==size(y)).and.(size(y)==size(z)))) then         !kontrola velkosti vstupnych poli
            write(*,*)'Problem s dimenziami vstupnych poli. - vec'
            return
    endif


    x=0
    y=0
    z=0
    do j1=0,j1max

            do m1=0,j1
                    do l1=abs(j1-1),j1+1
                            
                            ind1 = INDV(j1,m1,l1)
                            if (vec(ind1)==0) cycle

                            j=l1                        !j musi byt rovne l1 z vlastnosti c-g koef.

                            !spocitam konstantu mkonst, ktora nezavisi na m a len nou nasobim vysledok
                            sj=(-1)**(j+j1+1)/dsqrt((2.D0*j+1)*3)               !6-j symbol s nulou sa da vyjadrit takto jednoducho
                            mkonst=-dSQRT(2.D0*j1+1)*dSQRT(3.D0)*SQRT(2.D0*l1+1)/dSQRT(4*pi)
                            mkonst=mkonst*vec(ind1)*sj/dSQRT(2.D0*j+1)


                            !zlozky x,y,z
                                    indscal=INDS(j,0)
                                    call CLEB1(j1,m1,1,1,j,0,cg)
                                    x(indscal)=x(INDS(j,0))-mkonst*cg
                                    y(indscal)=y(INDS(j,0))+(0.D0,1.D0)*mkonst*cg
                                    call CLEB1(j1,m1,1,-1,j,0,cg)
                                    x(indscal)=x(INDS(j,0))+mkonst*cg
                                    y(indscal)=y(INDS(j,0))+(0.D0,1.D0)*mkonst*cg
                                    call CLEB1(j1,m1,1,0,j,0,cg)
                                    z(indscal)=z(INDS(j,0))+mkonst*cg        
                            do m=1,j
                                    indscal=INDS(j,m)
                                    call CLEB1(j1,m1,1,1,j,m,cg)
                                    x(indscal)=x(INDS(j,m))-mkonst*cg
                                    y(indscal)=y(INDS(j,m))+(0.D0,1.D0)*mkonst*cg
                                    call CLEB1(j1,m1,1,1,j,-m,cg)
                                    x(indscal)=x(INDS(j,m))-(-1)**m*CONJG(mkonst*cg)
                                    y(indscal)=y(INDS(j,m))+(-1)**m*CONJG((0.D0,1.D0)*mkonst*cg)
                                    call CLEB1(j1,m1,1,-1,j,m,cg)
                                    x(indscal)=x(INDS(j,m))+mkonst*cg
                                    y(indscal)=y(INDS(j,m))+(0.D0,1.D0)*mkonst*cg
                                    call CLEB1(j1,m1,1,-1,j,-m,cg)
                                    x(indscal)=x(INDS(j,m))+(-1)**m*CONJG(mkonst*cg)
                                    y(indscal)=y(INDS(j,m))+(-1)**m*CONJG((0.D0,1.D0)*mkonst*cg)
                                    call CLEB1(j1,m1,1,0,j,m,cg)        
                                    z(indscal)=z(INDS(j,m))+mkonst*cg
                                    call CLEB1(j1,m1,1,0,j,-m,cg)
                                    z(indscal)=z(INDS(j,m))+(-1)**m*CONJG(mkonst*cg)        
                            enddo
                    enddo
            enddo
    enddo

    x=x*dSQRT(2*pi)
    y=y*dSQRT(2*pi)
    z=z*2*dSQRT(pi)

END SUBROUTINE shz_vec


! ----------------------------------------------------------------------
FUNCTION INDV(j,m,l)
INTEGER :: INDV, j, m, l
! v e k t o r o v y   i n d e x pre m=0,j
INDV=(j*(j+1)/2+m)*3+l-j
END FUNCTION INDV
! ----------------------------------------------------------------------
FUNCTION INDSmin(j,m)
INTEGER :: INDSmin, j, m
! s k a l a r n y  i n d e x  pre m=-j,j
INDSmin=j*(j+1)+m+1
END FUNCTION INDSmin
! ----------------------------------------------------------------------
FUNCTION INDS(j,m)
INTEGER :: INDS, j, m
! s k a l a r n y  i n d e x  pre m=0,j
INDS=j*(j+1)/2+m+1
END FUNCTION INDS
! ----------------------------------------------------------------------
FUNCTION INDT(j,m,l,k)
INTEGER :: INDT, j, m, l, k
! t e n z o r o v y i n d e x pre m=0,j
INDT=(j+3)*(j-2)*9/2+17+m*9+k*k+l-j+k+1
if (j==0.and.m==0.and.l==0.and.k==0) INDT=1
if (j==0.and.m==0.and.l==1.and.k==1) INDT=2
if (j==0.and.m==0.and.l==2.and.k==2) INDT=3
if (j==1.and.m==0.and.l==1.and.k==0) INDT=4
if (j==1.and.m==0.and.l==0.and.k==1) INDT=5
if (j==1.and.m==0.and.l==1.and.k==1) INDT=6
if (j==1.and.m==0.and.l==2.and.k==1) INDT=7
if (j==1.and.m==0.and.l==1.and.k==2) INDT=8
if (j==1.and.m==0.and.l==2.and.k==2) INDT=9
if (j==1.and.m==0.and.l==3.and.k==2) INDT=10
if (j==1.and.m==1.and.l==1.and.k==0) INDT=11
if (j==1.and.m==1.and.l==0.and.k==1) INDT=12
if (j==1.and.m==1.and.l==1.and.k==1) INDT=13
if (j==1.and.m==1.and.l==2.and.k==1) INDT=14
if (j==1.and.m==1.and.l==1.and.k==2) INDT=15
if (j==1.and.m==1.and.l==2.and.k==2) INDT=16
if (j==1.and.m==1.and.l==3.and.k==2) INDT=17
END FUNCTION INDT
! ----------------------------------------------------------------------
FUNCTION INDVmin(j,m,l)
INTEGER :: INDVmin, j, m, l
! v e k t o r o v y   i n d e x pre m=-j,j
INDVmin=(j*(j+1)+m)*3+l-j
END FUNCTION INDVmin
! ----------------------------------------------------------------------


SUBROUTINE ninej0m(j1,j2,j3,j4,j6,j7,j8,nj)
!spocita 9-j Wignerov symbol, pre j9=0 a j5=1 
!        {j1,j2,j3}
!        {j4,j5,j6}
!        {j7,j8,j9}

INTEGER, INTENT(IN) :: j1,j2,j3,j4,j6,j7,j8
REAL(8) :: nj, sj

    call SIXJ1(j1,j2,j3,j4,j7,sj)

nj=(-1)**(j2+j3+j4+j7)/SQRT((2.D0*j3+1.D0)*(2.D0*j7+1.D0))*sj

if (j3.ne.j6.or.j7.ne.j8) nj=0


END SUBROUTINE ninej0m

SUBROUTINE Genersg (krok,sg)                      !nageneruje velkosti integracnych uhlov sg

INTEGER, INTENT(IN) :: krok
REAL (8) :: sg(0:180/krok), theta
INTEGER :: it

sg(0)=(1-cos(krok*pi/180/2))*2*pi
do it=1,180/krok-1
        theta=(it+1/2.D0)*krok*pi/180
        sg(it)= (1-cos(theta))*2*pi
enddo

do it=180/krok-1,1,-1
        sg(it)=(sg(it)-sg(it-1))/(360/krok)
enddo
sg(0)=sg(0)/(360/krok)
sg(180/krok)=sg(0)

END SUBROUTINE Genersg

!..................................................................
SUBROUTINE HARMAN(NTH,JMAX,DATA,CRHS)
!..................................................................
!     THIS SUBROUTINE CREATES THE ARRAY OF THE RIGHT-HAND SIDES
!     OF NORMAL EQUATIONS OCCURING IN THE LEAST SQUARES ADJUSTMENT
!
!     DESCRIPTION OF PARAMETERS:
!     NTH   THE NUMBER OF THE LATITUDE CIRCLES (NYQUIST FREQUENCY)
!           OF AN EQUAL ANGULAR GRID. THE STEP OF THE GRID IS THEN
!           STEP=180/NTH. IN ORDER TO ENSURE THE SYMMETRY OF THE
!           GRID WITH RESPECT TO THE EQUATOR, NTH NUST BE AN EVEN
!           NUMBER.
!     JMAX  THE DESIRED CUT-OFF DEGREE OF THE TRUNCATED SPHERICAL
!           HARMONIC SERIES FOR THE SPECTRAL REPRESENTATION OF A
!           FUNCTION. TO AVOID THE DISTORTION OF THE SPHERICAL
!           HARMONIC COEFFICIENTS BY THE ALIAS EFFECT, IT IS NECE-
!           SSARY TO PUT JMAX.LT.NTH.
!     DATA  TWO DIMENSIONAL ARRAY OF INPUT DATA VALUES. EACH COLUMN
!           OF THIS ARRAY CONTAINS DATA POINTS ALONG ONE PARALLEL.
!           THE PARALLELS ARE STORED FROM SOUTH TO NORTH POLE IN
!           SUCH A WAY THAT THE FIRST COLUMN IS THE CIRCLE WITH THE
!           180-STEP/2 CO-LATITUDE, AND THE LAST COLUMN OF THE
!           DATA ARRAY CONTAINS THE NORTHEST CIRCLE WITH THE
!           CO-LATITUDE STEP/2. THE FIRST ROW OF THE DATA ARRAY
!           CORRESPONDS TO THE ZERO-MERIDIAN. THE LONGITUDE
!           OF THE OTHER MERIDIANS INCREASES EASTWARDS WITH THE 
!           STEP-INCREAMENT.
!     CRHS  THE RIGHT-HAND SIDES OF NORMAL EQUATIONS. DUE
!           TO EVEN-ODD SYMMETRY FOR THE NEGATIVE ORDERS M, ONLY
!           THE COEFFICIENTS WITH NON-NEGATIVE ORDERS ARE COMPUTED.
!           THE COUPLE OF INDICES (J,M), J=0,1,...,JMAX; M=0,1,
!           ...,J, IS REDUCED TO THE ONE-DIMENSIONAL INDEX JM AS:
!
!    (1)    JM=J*(J+1)/2+M+1
!
!           THEREFORE, THE COEFFICIENTS ARE ARRANGED FIRST BY THE
!           DEGREE J, AND THEN BY THE ORDER M: 00,10,11,20,21,22,
!           ...,(JMAX,JMAX).
!     DIMENSIONS:
!           DATA(2*NTH,NTH),CRHS((JMAX+1)*(JMAX+2)/2)
!     DIMENS. OF AUXILIARY ARRAYS:
!           PNM(JMAX+1),RE(2*NTH+1),QIM(2*NTH+1)
!     NOTE:
!           THE PRESENT VALUES OF THE DIMENSIONS OF THE AUXILIARY
!           ARRAYS ARE SET-UP FOR NTH.LE.180 AND JMAX.LE.180.
!           IF THESE BOUNDS ARE OVERCOME, THE DIMENSIONS MUST BE
!           INCREASED.
!     THE USED SUBROUTINES:
!           C06ECF - CALCULATES FAST FOURIER TRANSFORM OF MIX-RADIX.
!                    IT WAS TAKEN FROM THE DOUBLE PRECISION NAG
!                    SUBROUTINE LIBRARY.
!           DPMM   - COMPUTES THE FULLY NORMALIZED ASSOCIATED LEGENDRE
!                    FUNCTIONS.
!........................................................
        IMPLICIT REAL*8(A,B,D-H,O-Z), COMPLEX*16(C)
        DIMENSION DATA(360,180),CRHS(*)
        DIMENSION PNM(181),RE(361),QIM(361)
    200 FORMAT(1X,'IFAIL=',I2)
        STEP=180.D00/DFLOAT(NTH)
        NPH=2*NTH
        NTH2=NTH/2
        FAC=DATAN(1.D00)/45.D00
        TH0=180.D00-STEP/2.D00
        JMAX1=JMAX+1
        DNFR=DFLOAT(NPH)
        SNFR=DSQRT(DNFR)
        DO 5 ITH2=1,NTH2
        THS=TH0-(ITH2-1)*STEP
        TH=THS*FAC
        XTH=DCOS(TH)
        ITHF=ITH2
        ITHS=NTH-ITH2+1
            DO 1 IPH=1,NPH
            RE(IPH)=DATA(IPH,ITHF)
    1    QIM(IPH)=DATA(IPH,ITHS)
        IFAIL=0
        CALL C06ECF(RE,QIM,NPH,IFAIL)
        IF(IFAIL.EQ.0) GOTO 2
        STOP
    2 CONTINUE
        RE(NPH+1)=RE(1)
        QIM(NPH+1)=QIM(1)
            DO 3 MS=1,JMAX1
            M=MS-1
            CALL DPMM(XTH,JMAX,M,PNM)
            NFMS=NPH+2-MS
            POM1=(RE(MS)+RE(NFMS))/2.D00
            POM2=(QIM(MS)-QIM(NFMS))/2.D00
            CAF=DCMPLX(POM1,POM2)*SNFR
            POM1=(QIM(MS)+QIM(NFMS))/2.D00
            POM2=(-RE(MS)+RE(NFMS))/2.D00
            CAS=DCMPLX(POM1,POM2)*SNFR
            ZN=-1.D00
            DO 3 JS=MS,JMAX1
            ZN=-ZN
            JM=(JS-1)*JS/2+MS
            JMP=JS-MS+1
            IF(ITH2.EQ.1) CRHS(JM)=(0.D00,0.D00)
    3       CRHS(JM)=CRHS(JM)+PNM(JMP)*(CAF+CAS*ZN)
    5 CONTINUE
        RETURN
        END
!..............................................................
        SUBROUTINE HARMLS(NTH,JMAX,CRHS,COEF)
!..............................................................
!     THE SUBROUTINE COMPUTES THE OPTIMAL LEAST SQUARES ESTIMATES
!     OF THE SPHERICAL HARMONIC COEFFICIENTS. THE ALGORITHM EXPLOITS
!     THE FACT THAT THE NORMAL MATRIX IS SPARSE FOR REGULARLY SPACED
!     GRID. THE NORMAL MATRIX IS REORDERED INTO THE BLOCK-DIAGONAL
!     FORM AND THEN EACH OF THE SUB-SYSTEM IS SOLVED BY THE NAG
!     SUBROUTINE F04ABF PERFORMING CHOLESKY'S FACTORIZATION.
!
!     DESCRIPTION OF PARAMETERS:
!     NTH    THE SAME AS IN HARMAN
!     JMAX   THE SAME AS IN HARMAN
!     CRHS   THE SAME AS IN HARMAN
!     COEF   THE LEAST SQUARES ESTIMATES OF THE SPHERICAL HARMONIC
!            COEFFICIENTS. THE ONE-DIMENSIONAL ARRAY COEF IS ARRANGED
!            ACCORDING TO THE INDEX (1)-VIZ SUBROUTINE HARMAN.
!     DIMENSIONS:
!            CRHS((JMAX+1)*(JMAX+2)/2),COEF((JMAX+1)*(JMAX+2)/2)
!     DIMENS. OF AUXILIARY ARRAYS:
!            PNM(JMAX+1),AMTRX(JMAX/2+1,JMAX/2+1),AUX(JMAX/2+1,2)
!            PARAM(JMAX/2+1,2),BB(JMAX/2+1,2),WORK(JMAX/2+1)
!     NOTE:
!            THE PRESENT VALUES OF THE DIMENSIONS OF THE AUXILIARY
!            ARRAYS ARE SET-UP FOR NTH.LE.180 AND JMAX.LE.180.
!            IF THESE BOUNDS ARE OVERCOME, THE DIMENSIONS MUST BE
!            INCREASED.
!     THE USED SUBROUTINES:
!            FO4ABF - CALCULATES THE ACCURATE SOLUTION OF A SET OF
!                     REAL SYMMETRIC POSITIVE DEFINITE LINEAR
!                     EQUATIONS WITH MULTIPLE RIGHT HAND SIDES BY
!                     CHOLESKY'S DECOMPOSITION METHOD.
!                     IT WAS TAKEN FROM THE DOUBLE PRECISION NAG 
!                     SUBROUTINE LIBRARY.
!            DPMM   - COMPUTES THE FULLY NORMALIZED ASSOCIATED LEGENDRE
!                     FUNCTIONS.
!..............................................................
        IMPLICIT REAL*8(A,B,D-H,O-Z), COMPLEX*16(C)
        DIMENSION CRHS(*),COEF(*)
        DIMENSION PNM(181),AMTRX(91,91),PARAM(91,2)
    200 FORMAT(1X,'IFAIL=',I2)
        STEP=180.D00/DFLOAT(NTH)
        NTH2=NTH/2
        FAC=DATAN(1.D00)/45.D00
        JMAX1=JMAX+1
        TH0=180.D00-STEP/2.D00
        DO 6 M1S=1,JMAX1
        M1=M1S-1
        M1S1=M1S+1
        DO 6 MP=M1S,M1S1
        NPM=(JMAX1-MP+2)/2
        IF(MP.GT.JMAX1) GOTO 6
            DO 2 ITH=1,NTH2
            THS=TH0-(ITH-1)*STEP
            TH=THS*FAC
            XTH=DCOS(TH)
            CALL DPMM(XTH,JMAX,M1,PNM)
            IJ1=0
            DO 1 J1S=MP,JMAX1,2
            IJ1=IJ1+1
            J1M1=J1S-M1S+1
            POM1=4.D00*NTH*PNM(J1M1)
            IF(DABS(POM1).LT.1.D-30) POM1=0.D00
            IJ2=0
            DO 1 J2S=MP,JMAX1,2
            IJ2=IJ2+1
            J2M1=J2S-M1S+1
            POM2=PNM(J2M1)
            IF(DABS(POM2).LT.1.D-30) POM2=0.D00
            IF(ITH.EQ.1) AMTRX(IJ1,IJ2)=0.D00
            AMTRX(IJ1,IJ2)=AMTRX(IJ1,IJ2)+POM1*POM2
    1       CONTINUE
    2    CONTINUE
            IJ1=0
            DO 3 J1S=MP,JMAX1,2
            IJ1=IJ1+1
            J1M1=(J1S-1)*J1S/2+M1S
            AMTRX(IJ1,NPM+1)=DBLE(CRHS(J1M1))
    3       AMTRX(IJ1,NPM+2)=DIMAG(CRHS(J1M1))
        CALL SIMZM2(NPM,AMTRX,2,PARAM)
            IJ1=0
            DO 5 J1S=MP,JMAX1,2
            IJ1=IJ1+1
            J1M1=(J1S-1)*J1S/2+M1S
    5    COEF(J1M1)=DCMPLX(PARAM(IJ1,1),PARAM(IJ1,2))
    6 CONTINUE
        RETURN
        END
!..................................................................
        SUBROUTINE HARMSY(NTH,JMAX,COEF,DATA)
!..................................................................
!     THIS SUBROUTINE EVALUATES THE TRUNCATED SUM OF THE SPHERICAL
!     HARMONICS COMPLETE TO DEGREE AND ORDER JMAX AT EACH ONE OF
!     THE 2*NTH**2 POINTS IN THE GRID.
!
!     DESCRIPTION OF PARAMETERS:
!          NTH   THE SAME AS IN HARMAN
!          JMAX  THE SAME AS IN HARMAN
!          COEF  THE SAME AS IN HARMLS
!          DATA  TWO-DIMENSIONAL ARRAY OF THE DATA VALUES COMPUTED
!                IN THE NODES OF THE REGULAR GRID. THIS ARRAY IS
!                ORGANIZED IN THE SAME WAY AS THE ARRAY DATA IN THE 
!                SUBROUTINE HARMAN.
!     DIMENSIONS:
!           DATA(2*NTH,NTH),COEF((JMAX+1)*(JMAX+2)/2)
!     DIMENS. OF AUXILIARY ARRAYS:
!           PNM(JMAX+1),RE1(2*NTH+1),QIM1(2*NTH+1),RE2(2*NTH+1),
!           QIM2(2*NTH+1)
!     NOTE:
!           THE PRESENT VALUES OF THE DIMENSIONS OF THE AUXILIARY
!           ARRAYS ARE SET-UP FOR NTH.LE.180 AND JMAX.LE.180.
!           IF THESE BOUNDS ARE OVERCOME, THE DIMENSIONS MUST BE
!           INCREASED.
!     THE USED SUBROUTINES:
!           C06ECF - CALCULATES FAST FOURIER TRANSFORM OF MIX-RADIX.
!                    IT WAS TAKEN FROM THE DOUBLE PRECISION NAG
!                    SUBROUTINE LIBRARY.
!           DPMM   - COMPUTES THE FULLY NORMALIZED ASSOCIATED LEGENDRE
!                    FUNCTIONS.
!..................................................................
        IMPLICIT REAL*8(A,B,D-H,O-Z), COMPLEX*16(C)
        DIMENSION DATA(360,180),COEF(*)
        DIMENSION PNM(181),RE1(361),QIM1(361),RE2(361),QIM2(361)
    200 FORMAT(1X,'IFAIL=',I2)
        STEP=180.D00/DFLOAT(NTH)
        NPH=2*NTH
        NTH2=NTH/2
        FAC=DATAN(1.D00)/45.D00
        TH0=180.D00-STEP/2.D00
        JMAX1=JMAX+1
        JMAX2=JMAX+2
        DNFR=DFLOAT(NPH)
        SNFR=DSQRT(DNFR)
        DO 5 ITH2=1,NTH2
        THS=TH0-(ITH2-1)*STEP
        TH=THS*FAC
        XTH=DCOS(TH)
            DO 2 MS=1,JMAX1
            M=MS-1
            CALL DPMM(XTH,JMAX,M,PNM)
            CSUM1=(0.D00,0.D00)
            CSUM2=(0.D00,0.D00)
            ZN=-1.D00
            DO 1 JS=MS,JMAX1
            ZN=-ZN
            JM=(JS-1)*JS/2+MS
            JMP=JS-MS+1
            CSUM1=CSUM1+COEF(JM)*PNM(JMP)
    1       CSUM2=CSUM2+COEF(JM)*PNM(JMP)*ZN
            RE1(MS)=DBLE(CSUM1)
            QIM1(MS)=-DIMAG(CSUM1)
            RE2(MS)=DBLE(CSUM2)
    2    QIM2(MS)=-DIMAG(CSUM2)
        RE1(1)=RE1(1)/2.D00
        RE2(1)=RE2(1)/2.D00
            DO 3 MS=JMAX2,NPH
            RE1(MS)=0.D00
            QIM1(MS)=0.D00
            RE2(MS)=0.D00
    3    QIM2(MS)=0.D00
        IFAIL=0
        CALL C06ECF(RE1,QIM1,NPH,IFAIL)
        CALL C06ECF(RE2,QIM2,NPH,IFAIL)
        IF(IFAIL.EQ.0) GOTO 4
        STOP
    4 CONTINUE
            DO 5 IPH=1,NPH
            DATA(IPH,ITH2)=2.D00*RE1(IPH)*SNFR
    5    DATA(IPH,NTH-ITH2+1)=2.D00*RE2(IPH)*SNFR
        RETURN
        END
!..................................................................
        SUBROUTINE DPMM(X,NN,M,P)
!..................................................................
        IMPLICIT REAL*8(D-H,O-Z)
        DIMENSION P(*)

        STH=DSQRT(1.D00-X*X)
        PI4=16.D00*DATAN(1.D00)
        P(1)=1.D00/DSQRT(PI4)
        IF(NN.LE.0) RETURN
        F1=0.D00
        IF(M.EQ.0) GOTO 1
            SOU=1.D00
            IFLAG=0
            DO 3 MP=1,M
            F1=DFLOAT(MP+MP)
            SOU=SOU*(F1+1.D00)/F1
            IF(IFLAG-1.LT.0) THEN
    10    STHM=STH**MP
            IF(STHM.LT.1.D-55) IFLAG=1
            GOTO 25
            ELSE
    20    STHM=0.D00
            ENDIF
    25    CONTINUE
    3    CONTINUE
        P(1)=(-1.D00)**M*DSQRT(SOU/PI4)*STHM
        IF(M.EQ.NN) GOTO 33
    1 P(2)=DSQRT(F1+3.D00)*X*P(1)
        MP2=M+2
        IF(MP2.GT.NN) GOTO 33
            I=1
            DO 2 J=MP2,NN
            F1=DFLOAT((J-M)*(J+M))
            F2=DFLOAT((J+J-1)*(J+J+1))
            F3=DFLOAT((J+J+1)*(J-M-1)*(J+M-1))/DFLOAT(J+J-3)
            F2=DSQRT(F2/F1)*X
            F3=DSQRT(F3/F1)
            I=I+1
            P(I+1)=F2*P(I)-F3*P(I-1)
    2    CONTINUE
    33 CONTINUE
        RETURN
        END
!....................................................................
        SUBROUTINE SIMZM2(N,A,NRHS,X)
!     Gaussova eliminace - bez pivotace a urcovani determinantu
!....................................................................
        REAL*8 A,X,PIVOT,AIK
        DIMENSION A(91,91),X(91,2)
        MAX=N+NRHS
        DO 18 K=1,N
        PIVOT=A(K,K)
        DO 14 J=1,MAX
    14 A(K,J)=A(K,J)/PIVOT
        A(K,K)=(1.D00,0.D00)/PIVOT
        DO 18 I=1,N
        AIK=A(I,K)
        IF(I.EQ.K) GO TO 18
        A(I,K)=-AIK/PIVOT
        DO 17 J=1,MAX
    17 IF(J.NE.K) A(I,J)=A(I,J)-AIK*A(K,J)
    18 CONTINUE
        do 20 IRHS=1,NRHS
        DO 20 I=1,N
    20 X(I,IRHS)=A(I,N+IRHS)
        RETURN
        END
        SUBROUTINE C06ECF(X, Y, PTS, IFAIL)
        INTEGER IFAIL, PTS
        REAL*8 X(*), Y(*)
!     DOUBLE PRECISION SRNAME
        REAL*8 SQPTS
        INTEGER IERROR, IPTS, PMAX, PSYM, TWOGRP
        INTEGER FACTOR(21), SYM(21), UNSYM(21)
        REAL*8 DSQRT
!     DATA SRNAME /'  C06ECE'/
        DATA PMAX /19/
        DATA TWOGRP /8/
        IF (PTS.LE.1) GO TO 40
        IERROR = 0
        CALL EAZC06(PTS,PMAX,TWOGRP,FACTOR,SYM,PSYM,UNSYM,IERROR)
        IF (IERROR.NE.0) GO TO 60
        CALL ECWC06(X, Y, PTS, FACTOR)
        CALL ECYC06(X, Y, PTS, SYM, PSYM, UNSYM)
        SQPTS = DSQRT(DFLOAT(PTS))
        DO 20 IPTS=1,PTS
            X(IPTS) = X(IPTS)/SQPTS
            Y(IPTS) = Y(IPTS)/SQPTS
    20 CONTINUE
        IFAIL = 0
        GO TO 80
    40 IERROR = 3
    60 CONTINUE
    80 RETURN
        END
        SUBROUTINE ECWC06(X, Y, PTS, FACTOR)
        INTEGER PTS
        REAL*8 X(PTS), Y(PTS)
        INTEGER FACTOR(21)
        INTEGER F, M1, M2, M3, M4, M5, M6, M7, M, P
        F = 0
        M = PTS
    20 CONTINUE
        F = F + 1
        P = FACTOR(F)
        IF (P.EQ.0) RETURN
        IF (P.EQ.1) GO TO 20
        M = M/P
        M1 = PTS - M
        M2 = M1 - M
        M3 = M2 - M
        M4 = M3 - M
        M5 = M4 - M
        M6 = M5 - M
        M7 = M6 - M
        IF (P.EQ.2) GO TO 40
        IF (P.EQ.3) GO TO 60
        IF (P.EQ.4) GO TO 80
        IF (P.EQ.5) GO TO 100
        IF (P.EQ.8) GO TO 120
        GO TO 140
    40 CONTINUE
        CALL ECVC06(X(1), Y(1), PTS, X(M+1), Y(M+1), M1, M)
        GO TO 20
    60 CONTINUE
        CALL ECUC06(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1),M2,M)
        GO TO 20
    80 CONTINUE
        CALL ECTC06(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1),M2,X(3*M+1),Y(3*M+1),M3,M)
        GO TO 20
    100 CONTINUE
        CALL ECSC06(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1),M2,X(3*M+1),Y(3*M+1),M3,X(4*M+1),Y(4*M+1),M4,M)
        GO TO 20
    120 CONTINUE
        CALL ECRC06(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1),M2,X(3*M+1),Y(3*M+1),M3,&
        X(4*M+1),Y(4*M+1),M4,X(5*M+1),Y(5*M+1),M5,X(6*M+1),Y(6*M+1),M6,X(7*M+1),Y(7*M+1),M7,M)
        GO TO 20
    140 CONTINUE
        CALL ECQC06(X, Y, PTS, M, P)
        GO TO 20
        END
        SUBROUTINE ECYC06(X, Y, PTS, SYM, PSYM, UNSYM)
        INTEGER PSYM, PTS
        REAL*8 X(PTS), Y(PTS)
        INTEGER SYM(21), UNSYM(21)
        REAL*8 T
        INTEGER DK,I,II,IL,J,JJ,JL,K,KK,KS,LK,MODS,MULT,NEST,PUNSYM,TEST
        LOGICAL ONEMOD
        INTEGER MODULO(20)
        DATA NEST /20/
        CALL ECXC06(X, Y, PTS, SYM)
        IF (UNSYM(1).EQ.0) GO TO 280
        PUNSYM = PTS/PSYM**2
        MULT = PUNSYM/UNSYM(1)
        TEST = (UNSYM(1)*UNSYM(2)-1)*MULT*PSYM
        LK = MULT
        DK = MULT
        DO 20 K=2,NEST
            IF (UNSYM(K).EQ.0) GO TO 40
            LK = LK*UNSYM(K-1)
            DK = DK/UNSYM(K)
            MODULO(K) = (LK-DK)*PSYM
            MODS = K
    20 CONTINUE
    40 CONTINUE
        ONEMOD = MODS.LT.3
        IF (ONEMOD) GO TO 80
        K = (MODS+3)/2
        DO 60 J=3,K
            JJ = MODS + 3 - J
            KK = MODULO(J)
            MODULO(J) = MODULO(JJ)
            MODULO(JJ) = KK
    60 CONTINUE
    80 CONTINUE
        JL = (PUNSYM-3)*PSYM
        KS = PUNSYM*PSYM
        DO 260 J=PSYM,JL,PSYM
            JJ = J
    100    CONTINUE
            JJ = JJ*MULT
            IF (ONEMOD) GO TO 140
            DO 120 I=3,MODS
            JJ = JJ - (JJ/MODULO(I))*MODULO(I)
    120    CONTINUE
    140    CONTINUE
            IF (JJ.GE.TEST) GO TO 160
            JJ = JJ - (JJ/MODULO(2))*MODULO(2)
            GO TO 180
    160    CONTINUE
            JJ = JJ - (JJ/MODULO(2))*MODULO(2) + MODULO(2)
    180    CONTINUE
            IF (JJ.LT.J) GO TO 100
            IF (JJ.EQ.J) GO TO 240
            LK = JJ - J
            II = J + 1
            IL = J + PSYM
            DO 220 I=II,IL
            DO 200 K=I,PTS,KS
                KK = K + LK
                T = X(K)
                X(K) = X(KK)
                X(KK) = T
                T = Y(K)
                Y(K) = Y(KK)
                Y(KK) = T
    200       CONTINUE
    220    CONTINUE
    240    CONTINUE
    260 CONTINUE
    280 CONTINUE
        RETURN
        END
        SUBROUTINE ECQC06(X, Y, PTS, M, P)
        INTEGER M, P, PTS
        REAL*8 X(PTS), Y(PTS)
        REAL*8 ANGLE, IS, IU, RS, RU, T, TWOPI, XT, YT
        INTEGER J, JJ, K0, K, KS1, KS2, MOVER2, MP, PM, PP, U, V
        LOGICAL FOLD, ZERO
        REAL*8 A(18),AA(9,9),B(18),BB(9,9),C(18),IA(9),IB(9),RA(9)
        REAL*8 RB(9),S(18)
        REAL*8 DCOS, DSIN
        twopi=8.0*datan(1.d00)
        MOVER2 = M/2 + 1
        MP = M*P
        PP = P/2
        PM = P - 1
        DO 20 U=1,PP
            JJ = P - U
            ANGLE = TWOPI*FLOAT(U)/FLOAT(P)
            A(U) = DCOS(ANGLE)
            B(U) = DSIN(ANGLE)
            A(JJ) = A(U)
            B(JJ) = -B(U)
    20 CONTINUE
        DO 60 U=1,PP
            DO 40 V=1,PP
            JJ = U*V - ((U*V)/P)*P
            AA(V,U) = A(JJ)
            BB(V,U) = B(JJ)
    40    CONTINUE
    60 CONTINUE
        DO 300 J=1,MOVER2
            FOLD = J.GT.1 .AND. 2*J.LT.M + 2
            K0 = J
            ANGLE = TWOPI*FLOAT(J-1)/FLOAT(MP)
            ZERO = ANGLE.EQ.0.0
            C(1) = DCOS(ANGLE)
            S(1) = DSIN(ANGLE)
            DO 80 U=2,PM
            C(U) = C(U-1)*C(1) - S(U-1)*S(1)
            S(U) = S(U-1)*C(1) + C(U-1)*S(1)
    80    CONTINUE
            GO TO 140
    100    CONTINUE
            FOLD = .FALSE.
            K0 = M + 2 - J
            DO 120 U=1,PM
            T = C(U)*A(U) + S(U)*B(U)
            S(U) = -S(U)*A(U) + C(U)*B(U)
            C(U) = T
    120    CONTINUE
    140    CONTINUE
            DO 280 K=K0,PTS,MP
            XT = X(K)
            YT = Y(K)
            KS1 = M + K
            KS2 = (P-1)*M + K
            RS = X(KS1) + X(KS2)
            IS = Y(KS1) + Y(KS2)
            RU = X(KS1) - X(KS2)
            IU = Y(KS1) - Y(KS2)
            DO 160 U=1,PP
                RA(U) = XT + RS*AA(U,1)
                IA(U) = YT + IS*AA(U,1)
                RB(U) = RU*BB(U,1)
                IB(U) = IU*BB(U,1)
    160       CONTINUE
            XT = XT + RS
            YT = YT + IS
            DO 200 U=2,PP
                JJ = P - U
                KS1 = U*M + K
                KS2 = JJ*M + K
                RS = X(KS1) + X(KS2)
                IS = Y(KS1) + Y(KS2)
                RU = X(KS1) - X(KS2)
                IU = Y(KS1) - Y(KS2)
                XT = XT + RS
                YT = YT + IS
                DO 180 V=1,PP
                    RA(V) = RA(V) + RS*AA(V,U)
                    IA(V) = IA(V) + IS*AA(V,U)
                    RB(V) = RB(V) + RU*BB(V,U)
                    IB(V) = IB(V) + IU*BB(V,U)
    180          CONTINUE
    200       CONTINUE
            X(K) = XT
            Y(K) = YT
            DO 260 U=1,PP
                JJ = P - U
                IF (ZERO) GO TO 220
                XT = RA(U) + IB(U)
                YT = IA(U) - RB(U)
                KS1 = U*M + K
                X(KS1) = XT*C(U) + YT*S(U)
                Y(KS1) = YT*C(U) - XT*S(U)
                XT = RA(U) - IB(U)
                YT = IA(U) + RB(U)
                KS1 = JJ*M + K
                X(KS1) = XT*C(JJ) + YT*S(JJ)
                Y(KS1) = YT*C(JJ) - XT*S(JJ)
                GO TO 240
    220          CONTINUE
                KS1 = U*M + K
                X(KS1) = RA(U) + IB(U)
                Y(KS1) = IA(U) - RB(U)
                KS1 = JJ*M + K
                X(KS1) = RA(U) - IB(U)
                Y(KS1) = IA(U) + RB(U)
    240          CONTINUE
    260       CONTINUE
    280    CONTINUE
            IF (FOLD) GO TO 100
    300 CONTINUE
        RETURN
        END
        SUBROUTINE ECRC06(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,X4,Y4,M4,X5,Y5,M5,X6,Y6,M6,X7,Y7,M7,M)
        INTEGER M1, M2, M3, M4, M5, M6, M7, M, PTS
        REAL*8 X0(PTS),X1(M1),X2(M2),X3(M3),X4(M4),X5(M5),X6(M6)
        REAL*8 X7(M7),Y0(PTS),Y1(M1),Y2(M2),Y3(M3),Y4(M4),Y5(M5)
        REAL*8 Y6(M6),Y7(M7)
        REAL*8 ANGLE, C1, C2, C3, C4, C5, C6, C7, E, I1, I2, I3, I4
        REAL*8 I5, I6, I7, IS0, IS1, IS2, IS3, ISS0, ISS1, ISU0, ISU1, IU0
        REAL*8 IU1, IU2, IU3, IUS0, IUS1, IUU0, IUU1, R1, R2, R3, R4, R5
        REAL*8 R6, R7, RS0, RS1, RS2, RS3, RSS0, RSS1, RSU0, RSU1, RU0
        REAL*8 RU1, RU2, RU3, RUS0, RUS1, RUU0, RUU1, S1, S2, S3, S4, S5
        REAL*8 S6, S7, T, TWOPI
        INTEGER J, K0, K, M8, MOVER2
        LOGICAL FOLD, ZERO
        REAL*8 DCOS, DSIN
        M8 = M*8
        MOVER2 = M/2 + 1
        twopi=8.0*datan(1.d00)
        E = DCOS(TWOPI/8.0D00)
        DO 120 J=1,MOVER2
            FOLD = J.GT.1 .AND. 2*J.LT.M + 2
            K0 = J
            ANGLE = TWOPI*FLOAT(J-1)/FLOAT(M8)
            ZERO = ANGLE.EQ.0.0
            C1 = DCOS(ANGLE)
            S1 = DSIN(ANGLE)
            C2 = C1*C1 - S1*S1
            S2 = S1*C1 + C1*S1
            C3 = C2*C1 - S2*S1
            S3 = S2*C1 + C2*S1
            C4 = C2*C2 - S2*S2
            S4 = S2*C2 + C2*S2
            C5 = C4*C1 - S4*S1
            S5 = S4*C1 + C4*S1
            C6 = C4*C2 - S4*S2
            S6 = S4*C2 + C4*S2
            C7 = C4*C3 - S4*S3
            S7 = S4*C3 + C4*S3
            GO TO 40
    20    CONTINUE
            FOLD = .FALSE.
            K0 = M + 2 - J
            T = (C1+S1)*E
            S1 = (C1-S1)*E
            C1 = T
            T = S2
            S2 = C2
            C2 = T
            T = (-C3+S3)*E
            S3 = (C3+S3)*E
            C3 = T
            C4 = -C4
            T = -(C5+S5)*E
            S5 = (-C5+S5)*E
            C5 = T
            T = -S6
            S6 = -C6
            C6 = T
            T = (C7-S7)*E
            S7 = -(C7+S7)*E
            C7 = T
    40    CONTINUE
            DO 100 K=K0,PTS,M8
            RS0 = X0(K) + X4(K)
            IS0 = Y0(K) + Y4(K)
            RU0 = X0(K) - X4(K)
            IU0 = Y0(K) - Y4(K)
            RS1 = X1(K) + X5(K)
            IS1 = Y1(K) + Y5(K)
            RU1 = X1(K) - X5(K)
            IU1 = Y1(K) - Y5(K)
            RS2 = X2(K) + X6(K)
            IS2 = Y2(K) + Y6(K)
            RU2 = X2(K) - X6(K)
            IU2 = Y2(K) - Y6(K)
            RS3 = X3(K) + X7(K)
            IS3 = Y3(K) + Y7(K)
            RU3 = X3(K) - X7(K)
            IU3 = Y3(K) - Y7(K)
            RSS0 = RS0 + RS2
            ISS0 = IS0 + IS2
            RSU0 = RS0 - RS2
            ISU0 = IS0 - IS2
            RSS1 = RS1 + RS3
            ISS1 = IS1 + IS3
            RSU1 = RS1 - RS3
            ISU1 = IS1 - IS3
            RUS0 = RU0 - IU2
            IUS0 = IU0 + RU2
            RUU0 = RU0 + IU2
            IUU0 = IU0 - RU2
            RUS1 = RU1 - IU3
            IUS1 = IU1 + RU3
            RUU1 = RU1 + IU3
            IUU1 = IU1 - RU3
            T = (RUS1+IUS1)*E
            IUS1 = (IUS1-RUS1)*E
            RUS1 = T
            T = (RUU1+IUU1)*E
            IUU1 = (IUU1-RUU1)*E
            RUU1 = T
            X0(K) = RSS0 + RSS1
            Y0(K) = ISS0 + ISS1
            IF (ZERO) GO TO 60
            R1 = RUU0 + RUU1
            I1 = IUU0 + IUU1
            R2 = RSU0 + ISU1
            I2 = ISU0 - RSU1
            R3 = RUS0 + IUS1
            I3 = IUS0 - RUS1
            R4 = RSS0 - RSS1
            I4 = ISS0 - ISS1
            R5 = RUU0 - RUU1
            I5 = IUU0 - IUU1
            R6 = RSU0 - ISU1
            I6 = ISU0 + RSU1
            R7 = RUS0 - IUS1
            I7 = IUS0 + RUS1
            X4(K) = R1*C1 + I1*S1
            Y4(K) = I1*C1 - R1*S1
            X2(K) = R2*C2 + I2*S2
            Y2(K) = I2*C2 - R2*S2
            X6(K) = R3*C3 + I3*S3
            Y6(K) = I3*C3 - R3*S3
            X1(K) = R4*C4 + I4*S4
            Y1(K) = I4*C4 - R4*S4
            X5(K) = R5*C5 + I5*S5
            Y5(K) = I5*C5 - R5*S5
            X3(K) = R6*C6 + I6*S6
            Y3(K) = I6*C6 - R6*S6
            X7(K) = R7*C7 + I7*S7
            Y7(K) = I7*C7 - R7*S7
            GO TO 80
    60       CONTINUE
            X4(K) = RUU0 + RUU1
            Y4(K) = IUU0 + IUU1
            X2(K) = RSU0 + ISU1
            Y2(K) = ISU0 - RSU1
            X6(K) = RUS0 + IUS1
            Y6(K) = IUS0 - RUS1
            X1(K) = RSS0 - RSS1
            Y1(K) = ISS0 - ISS1
            X5(K) = RUU0 - RUU1
            Y5(K) = IUU0 - IUU1
            X3(K) = RSU0 - ISU1
            Y3(K) = ISU0 + RSU1
            X7(K) = RUS0 - IUS1
            Y7(K) = IUS0 + RUS1
    80       CONTINUE
    100    CONTINUE
            IF (FOLD) GO TO 20
    120 CONTINUE
        RETURN
        END
        SUBROUTINE ECSC06(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,X4,Y4,M4,M)
        INTEGER M1, M2, M3, M4, M, PTS
        REAL*8 X0(PTS), X1(M1), X2(M2), X3(M3), X4(M4), Y0(PTS)
        REAL*8 Y1(M1), Y2(M2), Y3(M3), Y4(M4)
        REAL*8 A1, A2, ANGLE, AS, AU, B1, B2, C1, C2, C3, C4, I0, I1
        REAL*8 I2, I3, I4, IA1, IA2, IAS, IAU, IB1, IB2, IS1, IS2, ISS
        REAL*8 IU1, IU2, R0, R1, R2, R3, R4, RA1, RA2, RAS, RAU, RB1, RB2
        REAL*8 RS1, RS2, RSS, RU1, RU2, S1, S2, S3, S4, T, TWOPI
        INTEGER J, K0, K, M5, MOVER2
        LOGICAL FOLD, ZERO
        REAL*8 DCOS, DSIN, DSQRT
        M5 = M*5
        MOVER2 = M/2 + 1
        twopi=8.0*datan(1.d00)
        A1 = DCOS(TWOPI/5.0D00)
        B1 = DSIN(TWOPI/5.0D00)
        A2 = DCOS(2.0*TWOPI/5.0D00)
        B2 = DSIN(2.0*TWOPI/5.0D00)
        AS = -1.0/4.0
        AU = DSQRT(5.0D00)/4.0
        DO 120 J=1,MOVER2
            FOLD = J.GT.1 .AND. 2*J.LT.M + 2
            K0 = J
            ANGLE = TWOPI*FLOAT(J-1)/FLOAT(M5)
            ZERO = ANGLE.EQ.0.0
            C1 = DCOS(ANGLE)
            S1 = DSIN(ANGLE)
            C2 = C1*C1 - S1*S1
            S2 = S1*C1 + C1*S1
            C3 = C2*C1 - S2*S1
            S3 = S2*C1 + C2*S1
            C4 = C2*C2 - S2*S2
            S4 = S2*C2 + C2*S2
            GO TO 40
    20    CONTINUE
            FOLD = .FALSE.
            K0 = M + 2 - J
            T = C1*A1 + S1*B1
            S1 = C1*B1 - S1*A1
            C1 = T
            T = C2*A2 + S2*B2
            S2 = C2*B2 - S2*A2
            C2 = T
            T = C3*A2 - S3*B2
            S3 = -C3*B2 - S3*A2
            C3 = T
            T = C4*A1 - S4*B1
            S4 = -C4*B1 - S4*A1
            C4 = T
    40    CONTINUE
            DO 100 K=K0,PTS,M5
            R0 = X0(K)
            I0 = Y0(K)
            RS1 = X1(K) + X4(K)
            IS1 = Y1(K) + Y4(K)
            RU1 = X1(K) - X4(K)
            IU1 = Y1(K) - Y4(K)
            RS2 = X2(K) + X3(K)
            IS2 = Y2(K) + Y3(K)
            RU2 = X2(K) - X3(K)
            IU2 = Y2(K) - Y3(K)
            RSS = RS1 + RS2
            ISS = IS1 + IS2
            RAS = R0 + RSS*AS
            IAS = I0 + ISS*AS
            RAU = (RS1-RS2)*AU
            IAU = (IS1-IS2)*AU
            RA1 = RAS + RAU
            IA1 = IAS + IAU
            RA2 = RAS - RAU
            IA2 = IAS - IAU
            RB1 = RU1*B1 + RU2*B2
            IB1 = IU1*B1 + IU2*B2
            RB2 = RU1*B2 - RU2*B1
            IB2 = IU1*B2 - IU2*B1
            X0(K) = R0 + RSS
            Y0(K) = I0 + ISS
            IF (ZERO) GO TO 60
            R1 = RA1 + IB1
            I1 = IA1 - RB1
            R2 = RA2 + IB2
            I2 = IA2 - RB2
            R3 = RA2 - IB2
            I3 = IA2 + RB2
            R4 = RA1 - IB1
            I4 = IA1 + RB1
            X1(K) = R1*C1 + I1*S1
            Y1(K) = I1*C1 - R1*S1
            X2(K) = R2*C2 + I2*S2
            Y2(K) = I2*C2 - R2*S2
            X3(K) = R3*C3 + I3*S3
            Y3(K) = I3*C3 - R3*S3
            X4(K) = R4*C4 + I4*S4
            Y4(K) = I4*C4 - R4*S4
            GO TO 80
    60       CONTINUE
            X1(K) = RA1 + IB1
            Y1(K) = IA1 - RB1
            X2(K) = RA2 + IB2
            Y2(K) = IA2 - RB2
            X3(K) = RA2 - IB2
            Y3(K) = IA2 + RB2
            X4(K) = RA1 - IB1
            Y4(K) = IA1 + RB1
    80       CONTINUE
    100    CONTINUE
            IF (FOLD) GO TO 20
    120 CONTINUE
        RETURN
        END
        SUBROUTINE ECTC06(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,M)
        INTEGER M1, M2, M3, M, PTS
        REAL*8 X0(PTS),X1(M1),X2(M2),X3(M3),Y0(PTS),Y1(M1),Y2(M2),Y3(M3)
        REAL*8 ANGLE, C1, C2, C3, I1, I2, I3, IS0, IS1, IU0, IU1, R1
        REAL*8 R2, R3, RS0, RS1, RU0, RU1, S1, S2, S3, T, TWOPI
        INTEGER J, K0, K, M4, MOVER2
        LOGICAL FOLD, ZERO
        REAL*8 DCOS, DSIN
        M4 = M*4
        MOVER2 = M/2 + 1
        twopi=8.0*datan(1.d00)
        DO 120 J=1,MOVER2
            FOLD = J.GT.1 .AND. 2*J.LT.M + 2
            K0 = J
            ANGLE = TWOPI*FLOAT(J-1)/FLOAT(M4)
            ZERO = ANGLE.EQ.0.0
            C1 = DCOS(ANGLE)
            S1 = DSIN(ANGLE)
            C2 = C1*C1 - S1*S1
            S2 = S1*C1 + C1*S1
            C3 = C2*C1 - S2*S1
            S3 = S2*C1 + C2*S1
            GO TO 40
    20    CONTINUE
            FOLD = .FALSE.
            K0 = M + 2 - J
            T = C1
            C1 = S1
            S1 = T
            C2 = -C2
            T = C3
            C3 = -S3
            S3 = -T
    40    CONTINUE
            DO 100 K=K0,PTS,M4
            RS0 = X0(K) + X2(K)
            IS0 = Y0(K) + Y2(K)
            RU0 = X0(K) - X2(K)
            IU0 = Y0(K) - Y2(K)
            RS1 = X1(K) + X3(K)
            IS1 = Y1(K) + Y3(K)
            RU1 = X1(K) - X3(K)
            IU1 = Y1(K) - Y3(K)
            X0(K) = RS0 + RS1
            Y0(K) = IS0 + IS1
            IF (ZERO) GO TO 60
            R1 = RU0 + IU1
            I1 = IU0 - RU1
            R2 = RS0 - RS1
            I2 = IS0 - IS1
            R3 = RU0 - IU1
            I3 = IU0 + RU1
            X2(K) = R1*C1 + I1*S1
            Y2(K) = I1*C1 - R1*S1
            X1(K) = R2*C2 + I2*S2
            Y1(K) = I2*C2 - R2*S2
            X3(K) = R3*C3 + I3*S3
            Y3(K) = I3*C3 - R3*S3
            GO TO 80
    60       CONTINUE
            X2(K) = RU0 + IU1
            Y2(K) = IU0 - RU1
            X1(K) = RS0 - RS1
            Y1(K) = IS0 - IS1
            X3(K) = RU0 - IU1
            Y3(K) = IU0 + RU1
    80       CONTINUE
    100    CONTINUE
            IF (FOLD) GO TO 20
    120 CONTINUE
        RETURN
        END
        SUBROUTINE ECUC06(X0, Y0, PTS, X1, Y1, M1, X2, Y2, M2, M)
        INTEGER M1, M2, M, PTS
        REAL*8 X0(PTS), X1(M1), X2(M2), Y0(PTS), Y1(M1), Y2(M2)
        REAL*8 A, ANGLE, B, C1, C2, I0, I1, I2, IA, IB, IS, R0, R1, R2
        REAL*8 RA, RB, RS, S1, S2, T, TWOPI
        INTEGER J, K0, K, M3, MOVER2
        LOGICAL FOLD, ZERO
        REAL*8 DCOS, DSIN, DSQRT
        M3 = M*3
        MOVER2 = M/2 + 1
        twopi=8.0*datan(1.d00)
        A = -0.5
        B = DSQRT(0.75D00)
        DO 120 J=1,MOVER2
            FOLD = J.GT.1 .AND. 2*J.LT.M + 2
            K0 = J
            ANGLE = TWOPI*FLOAT(J-1)/FLOAT(M3)
            ZERO = ANGLE.EQ.0.0
            C1 = DCOS(ANGLE)
            S1 = DSIN(ANGLE)
            C2 = C1*C1 - S1*S1
            S2 = S1*C1 + C1*S1
            GO TO 40
    20    CONTINUE
            FOLD = .FALSE.
            K0 = M + 2 - J
            T = C1*A + S1*B
            S1 = C1*B - S1*A
            C1 = T
            T = C2*A - S2*B
            S2 = -C2*B - S2*A
            C2 = T
    40    CONTINUE
            DO 100 K=K0,PTS,M3
            R0 = X0(K)
            I0 = Y0(K)
            RS = X1(K) + X2(K)
            IS = Y1(K) + Y2(K)
            X0(K) = R0 + RS
            Y0(K) = I0 + IS
            RA = R0 + RS*A
            IA = I0 + IS*A
            RB = (X1(K)-X2(K))*B
            IB = (Y1(K)-Y2(K))*B
            IF (ZERO) GO TO 60
            R1 = RA + IB
            I1 = IA - RB
            R2 = RA - IB
            I2 = IA + RB
            X1(K) = R1*C1 + I1*S1
            Y1(K) = I1*C1 - R1*S1
            X2(K) = R2*C2 + I2*S2
            Y2(K) = I2*C2 - R2*S2
            GO TO 80
    60       CONTINUE
            X1(K) = RA + IB
            Y1(K) = IA - RB
            X2(K) = RA - IB
            Y2(K) = IA + RB
    80       CONTINUE
    100    CONTINUE
            IF (FOLD) GO TO 20
    120 CONTINUE
        RETURN
        END
        SUBROUTINE ECVC06(X0, Y0, PTS, X1, Y1, M1, M)
        INTEGER M1, M, PTS
        REAL*8 X0(PTS), X1(M1), Y0(PTS), Y1(M1)
        REAL*8 ANGLE, C, IS, IU, RS, RU, S, TWOPI
        INTEGER J, K0, K, M2, MOVER2
        LOGICAL FOLD, ZERO
        REAL*8 DCOS, DSIN
        M2 = M*2
        MOVER2 = M/2 + 1
        twopi=8.0*datan(1.d00)
        DO 120 J=1,MOVER2
            FOLD = J.GT.1 .AND. 2*J.LT.M + 2
            K0 = J
            ANGLE = TWOPI*FLOAT(J-1)/FLOAT(M2)
            ZERO = ANGLE.EQ.0.0
            C = DCOS(ANGLE)
            S = DSIN(ANGLE)
            GO TO 40
    20    CONTINUE
            FOLD = .FALSE.
            K0 = M + 2 - J
            C = -C
    40    CONTINUE
            DO 100 K=K0,PTS,M2
            RS = X0(K) + X1(K)
            IS = Y0(K) + Y1(K)
            RU = X0(K) - X1(K)
            IU = Y0(K) - Y1(K)
            X0(K) = RS
            Y0(K) = IS
            IF (ZERO) GO TO 60
            X1(K) = RU*C + IU*S
            Y1(K) = IU*C - RU*S
            GO TO 80
    60       CONTINUE
            X1(K) = RU
            Y1(K) = IU
    80       CONTINUE
    100    CONTINUE
            IF (FOLD) GO TO 20
    120 CONTINUE
        RETURN
        END
        SUBROUTINE ECXC06(X, Y, PTS, SYM)
        INTEGER PTS
        REAL*8 X(PTS), Y(PTS)
        INTEGER SYM(21)
        REAL*8 T
        INTEGER I10, I1, I2, I3, I4, I5, I6, I7, I8, I9, J, JJ, K10
        INTEGER K11, K1, K2, K3, K4, K5, K6, K7, K8, K9, KK, L10, L1, L2
        INTEGER L3, L4, L5, L6, L7, L8, L9, LEVEL, LOOP, NEST
        INTEGER I(20), K(20), L(20)
        DATA NEST /20/
        DATA LOOP /10/
        IF (SYM(1).EQ.0) GO TO 360
        DO 20 J=1,NEST
            L(J) = 1
            I(J) = 1
    20 CONTINUE
        KK = PTS
        DO 40 J=1,NEST
            IF (SYM(J).EQ.0) GO TO 60
            L(J) = KK
            I(J) = KK/SYM(J)
            KK = KK/SYM(J)
    40 CONTINUE
    60 CONTINUE
        L1 = L(1)
        L2 = L(2)
        L3 = L(3)
        L4 = L(4)
        L5 = L(5)
        L6 = L(6)
        L7 = L(7)
        L8 = L(8)
        L9 = L(9)
        L10 = L(10)
        I1 = I(1)
        I2 = I(2)
        I3 = I(3)
        I4 = I(4)
        I5 = I(5)
        I6 = I(6)
        I7 = I(7)
        I8 = I(8)
        I9 = I(9)
        I10 = I(10)
        KK = 0
        LEVEL = NEST
        K(LEVEL) = 1
        GO TO 100
    80 CONTINUE
        IF (LEVEL.GE.NEST) GO TO 360
        LEVEL = LEVEL + 1
        K(LEVEL) = K(LEVEL) + I(LEVEL)
        IF (K(LEVEL).GT.L(LEVEL)) GO TO 80
    100 CONTINUE
        LEVEL = LEVEL - 1
        DO 120 J=LOOP,LEVEL
            JJ = LEVEL + LOOP - J
            K(JJ) = K(JJ+1)
    120 CONTINUE
        K11 = K(11)
        DO 340 K10=K11,L10,I10
            K(10) = K10
            DO 320 K9=K10,L9,I9
            K(9) = K9
            DO 300 K8=K9,L8,I8
                K(8) = K8
                DO 280 K7=K8,L7,I7
                    K(7) = K7
                    DO 260 K6=K7,L6,I6
                        K(6) = K6
                        DO 240 K5=K6,L5,I5
                        K(5) = K5
                        DO 220 K4=K5,L4,I4
                            K(4) = K4
                            DO 200 K3=K4,L3,I3
                                K(3) = K3
                                DO 180 K2=K3,L2,I2
                                    K(2) = K2
                                    DO 160 K1=K2,L1,I1
                                    K(1) = K1
                                    KK = KK + 1
                                    IF (KK.GE.K1) GO TO 140
                                    T = X(KK)
                                    X(KK) = X(K1)
                                    X(K1) = T
                                    T = Y(KK)
                                    Y(KK) = Y(K1)
                                    Y(K1) = T
    140                               CONTINUE
    160                            CONTINUE
    180                         CONTINUE
    200                      CONTINUE
    220                   CONTINUE
    240                CONTINUE
    260             CONTINUE
    280          CONTINUE
    300       CONTINUE
    320    CONTINUE
    340 CONTINUE
        LEVEL = LOOP
        GO TO 80
    360 CONTINUE
        RETURN
        END
        SUBROUTINE EAZC06(PTS,PMAX,TWOGRP,FACTOR,SYM,PSYM,UNSYM,IERROR)
        INTEGER IERROR, PMAX, PSYM, PTS, TWOGRP
        INTEGER FACTOR(21), SYM(21), UNSYM(21)
        INTEGER F, J, JJ, N, NEST, P, PTWO, Q, R
        INTEGER PP(10), QQ(20)
        DATA NEST /20/
        N = PTS
        PSYM = 1
        F = 2
        P = 0
        Q = 0
    20 CONTINUE
        IF (N.LE.1) GO TO 100
        DO 40 J=F,PMAX
            IF (N.EQ.(N/J)*J) GO TO 60
    40 CONTINUE
        GO TO 280
    60 CONTINUE
        IF (2*P+Q.GE.NEST) GO TO 300
        F = J
        N = N/F
        IF (N.EQ.(N/F)*F) GO TO 80
        Q = Q + 1
        QQ(Q) = F
        GO TO 20
    80 CONTINUE
        N = N/F
        P = P + 1
        PP(P) = F
        PSYM = PSYM*F
        GO TO 20
    100 CONTINUE
        R = 1
        IF (Q.EQ.0) R = 0
        IF (P.LT.1) GO TO 140
        DO 120 J=1,P
            JJ = P + 1 - J
            SYM(J) = PP(JJ)
            FACTOR(J) = PP(JJ)
            JJ = P + Q + J
            FACTOR(JJ) = PP(J)
            JJ = P + R + J
            SYM(JJ) = PP(J)
    120 CONTINUE
    140 CONTINUE
        IF (Q.LT.1) GO TO 180
        DO 160 J=1,Q
            JJ = P + J
            UNSYM(J) = QQ(J)
            FACTOR(JJ) = QQ(J)
    160 CONTINUE
        SYM(P+1) = PTS/PSYM**2
    180 CONTINUE
        JJ = 2*P + Q
        FACTOR(JJ+1) = 0
        PTWO = 1
        J = 0
    200 CONTINUE
        J = J + 1
        IF (FACTOR(J).EQ.0) GO TO 240
        IF (FACTOR(J).NE.2) GO TO 200
        PTWO = PTWO*2
        FACTOR(J) = 1
        IF (PTWO.GE.TWOGRP) GO TO 220
        IF (FACTOR(J+1).EQ.2) GO TO 200
    220 CONTINUE
        FACTOR(J) = PTWO
        PTWO = 1
        GO TO 200
    240 CONTINUE
        IF (P.EQ.0) R = 0
        JJ = 2*P + R
        SYM(JJ+1) = 0
        IF (Q.LE.1) Q = 0
        UNSYM(Q+1) = 0
        IERROR = 0
    260 CONTINUE
        RETURN
    280 CONTINUE
        IERROR = 1
        GO TO 260
    300 CONTINUE
        IERROR = 2
        GO TO 260
        END
!     SUBROUTINE C06GCE(Y, PTS, IFAIL)
!     INTEGER IFAIL, PTS
!     REAL*8 Y(PTS)
!     DOUBLE PRECISION SRNAME
!     INTEGER IERROR, J
!     DATA SRNAME /'  C06GCE'/
!     IF (PTS.LE.0) GO TO 40
!     IERROR = 0
!     DO 20 J=1,PTS
!        Y(J) = -Y(J)
!  20 CONTINUE
!     GO TO 60
!  40 IERROR = 1
!  60 CONTINUE
!     RETURN
!     END

!c....................................................................
subroutine VCST(np,cajm,cbjml,cjml)
!c
!           vector-coupled sum of the product of
!c
!                 SCALAR  and  TENSOR
!c
!      vstup:
!               np       ..  maximalni stupen obou vstupnich rozvoju 
!               cajm     ..  koeficienty skalaru
!               cbjml    ..  koeficienty tensoru
!      vystup:
!               cjml     ..  koeficienty soucinu skalaru s tensorem
!c....................................................................
    implicit real*8(a,b,d-h,o-z), complex*16(c)
    dimension cajm(*),cbjml(*),cjml(*)
    complex*16, allocatable :: caux(:), cu(:,:), cv(:,:)

        ! Allocate dynamic arrays based on input parameters
    allocate(caux(5))          ! Temporary array for intermediate calculations
    allocate(cu(5, 7000))      ! Dynamic tensor data
    allocate(cv(5, 7000))      ! Resultant tensor data

!c
    cunit=dcmplx(0.d00,1.d00)
    np2=np+2
!c
    do 1 j=np+1,np+2
    do 1 ms=1,j+1
    m=ms-1
    jm=j*(j+1)/2+m+1
1 cajm(jm)=(0.d00,0.d00)
!c
    call CONVCU(np,cbjml,cu)
    call VCSUM5(np2,cajm,cu,cv)
!c
    do 10 js=1,np2+1
    j=js-1
    do 10 ms=1,js
    m=ms-1
    do 8 l=iabs(j-2),j+2
    call indt0(j,m,l,jml)
    csum=(0.d00,0.d00)
    do 6 iss=1,5
    is=iss-3
    mi=m-is
    call cleb2(l,mi,2,is,j,m,qq)
    mia=iabs(mi)
    lmi=l*(l+1)/2+mia+1
    do 7 i=1,5
    caux(i)=cv(i,lmi)
7 if(mi.lt.0) caux(i)=(-1.d00)**mia*dconjg(caux(i))
    if(is.eq.(-2)) cgama=(caux(1)+cunit*caux(5))/2.d00
    if(is.eq.(-1)) cgama=(caux(2)-cunit*caux(4))/2.d00
    if(is.eq.0) cgama=caux(3)
    if(is.eq.1) cgama=-(caux(2)+cunit*caux(4))/2.d00
    if(is.eq.2) cgama=(caux(1)-cunit*caux(5))/2.d00
    csum=csum+cgama*qq
6 continue
    cjml(jml)=csum  
8 continue
10 continue
    return

    ! Deallocate dynamically allocated arrays
    deallocate(caux, cu, cv)
    end
!c....................................................................
subroutine CONVCU(np,ctjml,cu)
!c....................................................................
    implicit real*8(a,b,d-h,o-z), complex*16(c)
    dimension ctjml(*),cu(5,7000)
!c
    cunit=dcmplx(0.d00,1.d00)
!c
    do 5 iss=1,5
    is=iss-3
    do 3 ls=1,np+3
    l=ls-1
    do 3 mis=1,ls
    mi=mis-1
    m=mi+is
    ma=iabs(m)
    lmi=l*(l+1)/2+mi+1
    csum=(0.d00,0.d00)
!     if(ma.gt.np) goto 3
    ipom=iabs(l-2)
    jdol=max0(ma,ipom)+1
    jhor=min0(np,l+2)+1
    do 2 js=jdol,jhor
    j=js-1
    call indt0(j,m,l,jml)
    cpom=ctjml(jml)
    if(m.lt.0) cpom=(-1.d00)**(j+m+l)*dconjg(cpom)
    call CLEB2(l,mi,2,is,j,m,qq)
    csum=csum+cpom*qq
2 continue
4 if(is) 11,12,13
11 if(is.eq.(-2)) cu(1,lmi)=csum
    if(is.eq.(-1)) cu(2,lmi)=csum
    goto 3
12 cu(3,lmi)=csum
    goto 3
13 if(is.eq.1) goto 14 
    cpom1=cu(1,lmi)
    cu(1,lmi)=cpom1+csum
    cu(5,lmi)=cunit*(-cpom1+csum)
    goto 3
14 cpom1=cu(2,lmi)
    cu(2,lmi)=cpom1-csum
    cu(4,lmi)=cunit*(cpom1+csum)
3 continue
5 continue 
!c
    return
end
!C.................................................................
SUBROUTINE VCSUM5(NP,CAJM,CU,CV)
!C.................................................................
    IMPLICIT REAL*8(D-H,O-Z),COMPLEX*16(A-C)
    DIMENSION CAJM(*),CU(5,7000),CV(5,7000)
    DIMENSION CSUMA(4,200),COMP(4)
    !DIMENSION RE(10192),QIM(10192)
    real*8, allocatable :: RE(:), QIM(:)

    COMMON /DD0/ NGQ,NFOURC,ROOTS(200),WGHTS(200),PNMARR(200,7000),CEMARR(10192,200)

    ! Allocate large arrays dynamically
    allocate(RE(10192))
    allocate(QIM(10192))

!C
    nfour=nfourc
    NROOT=NGQ/2
    NP1=NP+1
    DNFR=DFLOAT(NFOUR)
    PI2=8.D00*DATAN(1.D00)
    PI4=PI2+PI2
!C
    DO 20 IR=1,NROOT
!     ROOT=ROOTS(IR)
    POMWT=PI4*WGHTS(IR)
    do 30 ip=1,5
    DO 6 MS=1,NP1
    DO 4 K=1,4
4  CSUMA(K,MS)=(0.D00,0.D00)
    DO 6 JS=MS,NP1
    JM=(JS-1)*JS/2+MS
    CPOM=CU(ip,JM)
    ZK=(-1.D00)**(JS+MS)
    CSUMA(1,MS)=CSUMA(1,MS)+CPOM*PNMARR(IR,JM)
    CSUMA(2,MS)=CSUMA(2,MS)+CPOM*PNMARR(IR,JM)*ZK
    CSUMA(3,MS)=CSUMA(3,MS)+CAJM(JM)*PNMARR(IR,JM)
    CSUMA(4,MS)=CSUMA(4,MS)+CAJM(JM)*PNMARR(IR,JM)*ZK
6 CONTINUE
!    EVALUATE THE PRODUCTS AROUND THE LATITUDE CIRCLE
    DO 11 JFI=1,NFOUR
    DO 7 K=1,4
7 COMP(K)=(0.D00,0.D00)
    DO 9 M=1,NP
    MS=M+1
    DO 8 K=1,4
8 COMP(K)=COMP(K)+CSUMA(K,MS)*CEMARR(JFI,M)
9 CONTINUE
    DO 10 K=1,4
10 COMP(K)=COMP(K)+DCONJG(COMP(K))+CSUMA(K,1)
    FABP=DBLE(COMP(1)*COMP(3))
    FABN=DBLE(COMP(2)*COMP(4))
    RE(JFI)=FABP
11 QIM(JFI)=FABN
!    PERFORM THE FFT OF TWO REAL SERIES
    CALL DFFTAR(RE,QIM,NFOUR,0)
    RE(NFOUR+1)=RE(1)
    QIM(NFOUR+1)=QIM(1)
    DO 12 MS=1,NP1
    NFMS=NFOUR+2-MS
    POM1=(RE(MS)+RE(NFMS))/2.D00/DNFR
    POM2=(QIM(MS)-QIM(NFMS))/2.D00/DNFR
    CAP=DCMPLX(POM1,POM2)
    POM1=(QIM(MS)+QIM(NFMS))/2.D00/DNFR
    POM2=(-RE(MS)+RE(NFMS))/2.D00/DNFR
    CAN=DCMPLX(POM1,POM2)
!    FINALLY, EVALUATE THE IR-TH TERM OF GAUSS-LEGENDRE QUADRATURE
    DO 12 JS=MS,NP1
    JM=(JS-1)*JS/2+MS
    ZN=(-1.D00)**(JS+MS)
    IF(IR.EQ.1) CV(ip,JM)=(0.D00,0.D00)
12 CV(ip,JM)=CV(ip,JM)+POMWT*PNMARR(IR,JM)*(CAP+CAN*ZN)
30 continue
20 CONTINUE
    RETURN
    deallocate(RE, QIM)
    END
!....................................................................
subroutine indt0(j,m,l,ix)
!
!    index tensoru 2.radu
!....................................................................
    if(j.gt.1) goto 1
    ix=3*(j*(j+1)/2+iabs(m))+l-j-1
    return
1 ix=5*(j*(j+1)/2+iabs(m))+l-j-5
    return
end
!......................................................................     
SUBROUTINE CLEB2(J1,M1,J2,M2,J,M,C2)
!.....................................................................
!      UCEL:   PODPROGRAM SLOUZI K VYPOCTU C.-G. KOEFICIENTU S INDEXEM
!              J2=2 (VIZ VARSALOVIC ET AL., TABULKA 8.4, STR.228-229).
!       VSTUP:  INDEXY J,M (HORNI DVOJICE INDEXU),
!               INDEXY J1,M1 (INDEXY VLEVO DOLE),
!               INDEX M2 (ZAVISLY INDEX VPRAVO DOLE, ABS(M2).LE.2).
!      VYSTUP: C2 - PRISLUSNY C.-G. KOEFICIENT
!
!                                    J,M
!              C2(J,M,J1,M1,2,M2) = C
!                                    J1,M1 2,M2
!......................................................................
        REAL*8 C2
!
    C2=0.D00
!
    if(j2.ne.2) return
    IF((M1+M2).NE.M) return
    IF(IABS(M).GT.J) return
    IF(IABS(M1).GT.J1) return
    IF(iabs(j-j1).gt.2) return
    IF((j+j1).LT.2) return
!
    IF(M2) 1,2,3
!
1    I=-M
    GOTO 5
!
2   	K=J+M
L=J-M
    N=J+J
        IF(J.EQ.(J1+2)) THEN
    C2=DFLOAT(3*(K-1)*K*(L-1)*L)/DFLOAT((N-3)*(N-2)*(N-1)*J)
        C2=DSQRT(C2)
        ELSEIF(J.EQ.(J1+1)) THEN
        C2=DFLOAT(3*K*L)/DFLOAT((J-1)*(N-1)*J*(J+1))
        C2=DFLOAT(M)*DSQRT(C2)
    ELSEIF(J.EQ.J1) THEN
    C2=DFLOAT(3*M*M-J*(J+1))/DSQRT(DFLOAT((N-1)*J*(J+1)*(N+3)))
        ELSEIF(J.EQ.(J1-1)) THEN
        C2=DFLOAT(3*(K+1)*(L+1))/DFLOAT(J*(J+1)*(N+3)*(J+2))
    C2=-DFLOAT(M)*DSQRT(C2)
    ELSE
    C2=DSQRT(DFLOAT(3*(K+1)*(K+2)*(L+1)*(L+2))/DFLOAT((J+1)*(N+3)*(N+4)*(N+5)))
    ENDIF
GOTO 99
!
3    I=M
5    K=J+I
L=J-I
N=J+J
IF(IABS(M2).EQ.2) GOTO 4
        IF(J.EQ.(J1+2)) THEN
    C2=DFLOAT((K-2)*(K-1)*K*L)/DFLOAT((N-3)*(J-1)*(N-1)*J)
    C2=DSQRT(C2)
    ELSEIF(J.EQ.(J1+1)) THEN
    C2=DFLOAT((K-1)*K)/DFLOAT((N-2)*(N-1)*J*(J+1))
    C2=DFLOAT(-M2*(J-I-I+1))*DSQRT(C2)
    ELSEIF(J.EQ.J1) THEN
    C2=DFLOAT(3*K*(L+1))/DFLOAT(2*(N-1)*J*(J+1)*(N+3))
    C2=DFLOAT(1-I-I)*DSQRT(C2)
    ELSEIF(J.EQ.(J1-1)) THEN
    C2=DFLOAT((L+1)*(L+2))/DFLOAT(J*(J+1)*(N+3)*(N+4))
    C2=DFLOAT(M2*(J+I+I))*DSQRT(C2)
    ELSE
    C2=-DSQRT(DFLOAT((K+1)*(L+1)*(L+2)*(L+3))/DFLOAT((J+1)*(N+3)*(J+2)*(N+5)))
    ENDIF
    GOTO 99
!
4    IS=-M2/IABS(M2)
    IF(J.EQ.(J1+2)) THEN
    C2=DFLOAT((K-3)*(K-2)*(K-1)*K)/DFLOAT((N-3)*(N-2)*(N-1)*N)
    C2=DSQRT(C2)
    ELSEIF(J.EQ.(J1+1)) THEN
    C2=DFLOAT((K-2)*(K-1)*K*(L+1))/DFLOAT((N-2)*(N-1)*J*(J+1))
    C2=DFLOAT(IS)*DSQRT(C2)
    ELSEIF(J.EQ.J1) THEN
    C2=DFLOAT(3*(K-1)*K*(L+1)*(L+2))/DFLOAT(2*(N-1)*J*(J+1)*(N+3))
    C2=DSQRT(C2)
    ELSEIF(J.EQ.(J1-1)) THEN
    C2=DFLOAT(K*(L+1)*(L+2)*(L+3))/DFLOAT(J*(J+1)*(N+3)*(N+4))
    C2=DFLOAT(IS)*DSQRT(C2)
    ELSE
C2=DSQRT(DFLOAT((L+1)*(L+2)*(L+3)*(L+4))/DFLOAT((N+2)*(N+3)*(N+4)*(N+5)))
ENDIF
!
99   RETURN
    END
!.................................................................
SUBROUTINE GNDD0(NP)
!
!   Generovani pomocnych poli pro subr. VCDD0 a VCDERD.
!            Dimenze nastaveny pro NP=10.
!.................................................................
    IMPLICIT REAL*8(D-H,O-Z),COMPLEX*16(A-C)
    REAL*8 LEGE
    EXTERNAL LEGE
    LOGICAL IFLAG
    DIMENSION PNM(7000)
    COMMON /CORDER/ DUMMY(2),NGQ1
    COMMON /DD0/ NGQ,NFOUR,ROOTS(200),WGHTS(200),PNMARR(200,7000),CEMARR(10192,200)
!
!     FIND THE DEGREE OF GAUSS-LEGENDRE QUADRATURE FORMULA
    NGQ=3*(NP+2)/2+1
    NGQ=(NGQ+1)/2*2
    NGQ=MAX0(2*(NP+1),NGQ)
    NGQ1=NGQ
!
    NFOUR=MAX0(4*NP+1,3*(NP+2)+1)
    DO 1 K=1,20
    NAUX=2**K
1 IF(NAUX.GE.NFOUR) GOTO 2
2 NFOUR=NAUX
!
!      IRS=.TRUE.
    CALL NODES(LEGE,ROOTS,IFLAG)
    IF(IFLAG) GOTO 3
!     IF NOT SUCCESSFUL, SET IRS=.FALSE. AND RETURN
!      IRS=.FALSE.
    RETURN
3 CALL WGHLE(WGHTS,ROOTS)
!    DEFINE THE AUXILIARY PARAMETERS
!
    NROOT=NGQ/2
    NP3=NP+3
    DNFR=DFLOAT(NFOUR)
    PI2=8.D00*DATAN(1.D00)
    npall=np3*(np3+1)/2
!
    DO 4 IR=1,NROOT
    ROOT=ROOTS(IR)
    CALL DPNM(ROOT,NP+2,PNM)
    do 4 i=1,npall
4 pnmarr(ir,i)=pnm(i)
!
    DO 5 JFI=1,NFOUR
    PH=DFLOAT(JFI-1)*PI2/DNFR
    DO 5 M=1,NP+2
    REM=PH*DFLOAT(M)
    CEM=DCMPLX(0.D00,REM)
5 CEMARR(JFI,M)=CDEXP(CEM)
!
    RETURN
END
!.................................................................
!
!     SUBROUTINE DFFTAR
!
!     THIS PROGRAM IMPLEMENTS THE FFT ALGORITHM TO COMPUTE
!     THE DISCRETE FOURIER COEFFICIENTS OF A DATA SEQUENCE
!     OF N POINTS
!
!     CALLING SEQUENCE FROM THE MAIN PROGRAM:
!     CALL DFFTAR(A,B,N,INV)
!     N: NUMBER OF DATA POINTS
!     A,B: REAL AND IMAGINARY PARTS OF THE DATA SEQUENCE.
!          IN THE END DFT COEFFS. ARE RETURNED IN THE ARRAY.
!     INV: FLAG FOR INVERSE
!          INV=0  FOR FORWARD TRANSFORM
!          INV=1  FOR INVERSE TRANSFORM
!
!     FOR REFERENCE SEE: N.AHMED, K.R.RAO: ORTHOGONAL TRANSFORMS
!     FOR DIGITAL SIGNAL PROCESSING. SPRINGER VERLAG, BERLIN, 1975
!    (ALSO RUSSIAN TRANSLATION, MOSKVA, 1980)
!
!...................................................................
    SUBROUTINE DFFTAR(A,B,N,INV)
    REAL*8 A(*),B(*),WPWR,ARG,C,S,R,Q,T,PI
    PI=4.D00*DATAN(1.D00)
!     COMPUTE THE NUMBER OF ITERATIONS (LOG.N TO THE BASE 2)
    ITER=0
    IREM=N
10 IREM=IREM/2
    IF(IREM.EQ.0) GOTO 20
    ITER=ITER+1
    GOTO 10
20 CONTINUE
    ASIGN=-1.
    IF(INV.EQ.1) ASIGN=1.
    NXP2=N
    DO 50 IT=1,ITER
!   COMPUTATION FOR EACH ITERATION
!     NXP:   NUMBER OF POINTS IN A PARTITION
!     NXP2:  NXP/2
    NXP=NXP2
    NXP2=NXP/2
    WPWR=PI/DFLOAT(NXP2)
    DO 40 M=1,NXP2
!     CALCULATE THE MULTIPLIER
    ARG=DFLOAT(M-1)*WPWR
    C=DCOS(ARG)
    S=ASIGN*DSIN(ARG)
    DO 40 MXP=NXP,N,NXP
!     COMPUTATION FOR EACH PARTITION
    J1=MXP-NXP+M
    J2=J1+NXP2
    R=A(J1)-A(J2)
    Q=B(J1)-B(J2)
    A(J1)=A(J1)+A(J2)
    B(J1)=B(J1)+B(J2)
    A(J2)=R*C-Q*S
40 B(J2)=R*S+Q*C
50 CONTINUE
!     UNSCRAMBLE THE BIT-REVERSED DFT COEFFS.
    N2=N/2
    N1=N-1
    J=1
    DO 65 I=1,N1
    IF(I.GE.J) GOTO 55
    T=A(J)
    A(J)=A(I)
    A(I)=T
    T=B(J)
    B(J)=B(I)
    B(I)=T
55 K=N2
60 IF(K.GE.J) GOTO 65
    J=J-K
    K=K/2
    GOTO 60
65 J=J+K
    IF(INV.EQ.0) GOTO 75
    DO 70 I=1,N
    A(I)=A(I)/DFLOAT(N)
70 B(I)=B(I)/DFLOAT(N)
75 CONTINUE
    RETURN
END
!--------------------------------------------------------------
SUBROUTINE NWLEGP(NP,NR)
!--------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    REAL*8 LEGE
    EXTERNAL LEGE
    LOGICAL IFLAG,IRS
    COMMON / GQCOM / WGHTS(200),ROOTS(200),NGQ
    COMMON / CORDER / DUMMY(2),NLEG
!
    NGQ=NP+NP+1
    IF(NR.LE.NP) NGQ=3*NP/2+1
    NGQ=(NGQ+1)/2*2
    NLEG=NGQ
    IRS=.TRUE.
! --- DETERMINE NODES OF THE LEGENDRE POLYNOMIAL
    CALL NODES(LEGE,ROOTS,IFLAG)
    IF (IFLAG) GO TO 1
! --- IF NOT SUCCESSFULL, SET IRS=.FALSE. AND STOP
    IRS=.FALSE.
    STOP
! --- DETERMINE WEIGHTS
1 CALL WGHLE(WGHTS,ROOTS)
    RETURN
    END

SUBROUTINE NODES(P,T,IFLAG)                    
!     ===========================                                       AAPE0321
! --- THIS SUBR. DETERMINES ALL NODES OF THE FUNCTION 'P'.              AAPE0322
! --- THE PRESENT ENTRY IS USED FOR 'P' BEING THE LEGENDRE              AAPE0323
! --- POLYNOMIAL. THE NODES ARE STORED INTO THE ARRAY 'T' IN            AAPE0324
! --- THE ORDER 0<T(1)<T(2)<...<T(N/2) WHERE N IS THE DEGREE            AAPE0325
! --- OF THE POLYNOMIAL. T(N/2+1)=0 FOR N ODD.                          AAPE0326
! --- THE NODES ARE DETERMINED WITH THE RELATIVE ACCURACY OF 1D-15      AAPE0327
!                                                                       AAPE0328
    IMPLICIT REAL*8 (A-H,O-Z)                                         
    LOGICAL IFLAG                                                     
    COMMON / CORDER / ALPHA,BETA,N                                    
    DIMENSION T(*)                                                    
    EXTERNAL P                                                        

    MM=N/2                                                            
    XMIN=0D0                                                          
    IF (2*MM.NE.N) XMIN=1D-9                                          
    XINCR=(1D0-XMIN)/DFLOAT(N)                                        
    T(MM+1)=0D0                                                       
    NN=N                                                              
    J=0                                                               
    IFLAG=.TRUE.                                                      
10 NCNT=0                                                            
    X=XMIN                                                            
    FX=P(X)                                                           
    DO 1 I=1,NN                                                       
    Y=X+XINCR                                                         
    FY=P(Y)                                                           
    IF (FX*FY.LT.0D0) NCNT=NCNT+1                                     
    X=Y                                                               
1 FX=FY                                                            
    IF (NCNT.EQ.MM) GO TO 2                                           
    J=J+1                                                             
    NN=NN*2                                                           
    XINCR=XINCR/2D0                                                   
    IF (J.LT.10) GO TO 10                                             
    IFLAG=.FALSE.                                                     
    RETURN                                                            
2 X=XMIN                                                            
    I=1                                                               
    FX=P(X)                                                           
4 Y=X+XINCR                                                         
    FY=P(Y)                                                           
    IF (FX*FY.LT.0D0) GO TO 3                                         
8 X=Y                                                               
    FX=FY                                                             
    GO TO 4                                                           
3 T(I)=XNODE(P,X,Y,FX,FY)                                           
    I=I+1                                                             
    IF (I.GT.MM) RETURN                                               
    GO TO 8                                                           
    END                                                               
    FUNCTION LEGE(X)                                                  
    IMPLICIT REAL*8(L,P,X)                                            
    COMMON / CORDER / PDUMMY(2),N                                     
    IF (N-1) 1,2,3                                                    
1 LEGE=1D0                                                          
    RETURN                                                            
2 LEGE=X                                                            
    RETURN                                                            
3 P1=1D0                                                            
    P2=X                                                              
    DO 4 I=2,N                                                        
    XX=P2*X                                                           
    PP=XX-P1+XX-(XX-P1)/DFLOAT(I)                                     
    P1=P2                                                             
4 P2=PP                                                             
    LEGE=P2                                                           
    RETURN                                                            
    END                                                               
    SUBROUTINE WGHLE(A,T)                                             
    IMPLICIT REAL*8 (A-H,L,O-Z)                                       
    COMMON / CORDER / DUMMY(2),N                                      
    DIMENSION A(*),T(*)                                               
    M=N/2                                                             
    IF (2*M.NE.N) M=M+1                                               
    ORD=DFLOAT(N)**2                                                  
    N=N-1                                                             
    DO 1 I=1,M                                                        
    A(I)=(1D0-T(I)**2)/(ORD*LEGE(T(I))**2)                            
1 CONTINUE                                                          
    N=N+1                                                             
    RETURN                                                            
    END                                                               
    FUNCTION XNODE(P,X,Y,FX,FY)                                       
    IMPLICIT REAL*8 (A-H,O-Z)                                         
    X1=X                                                              
    X2=Y                                                              
    FX1=FX                                                            
    FX2=FY                                                            
7 TT=(X1+X2)/2D0                                                    
    FFX=P(TT)                                                         
    IF (DABS(FFX).LT.1D-15) GO TO 29                                
    IF (FX1*FFX.LT.0D0) GO TO 21
    FX1=FFX                                                           
    X1=TT                                                             
    GO TO 22                                                          
21 FX2=FFX                                                           
    X2=TT                                                             
22 IF ((X2-X1)/DABS(X1).GT.1D-15) GO TO 7                            
    XNODE=(X1+X2)/2D0                                                 
    RETURN                                                            
29 XNODE=TT                                                          
    RETURN                                                            
END           

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

    do j=1, jmax

        j1 = real(j)
    
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
        imag_=-ksi*sqrt(48*4.D0*datan(1.D0))*dsin(angular_speed*t_*delta_t)
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

subroutine update_cauchy_integral(jmax, number_of_layers, delta_t, eta, mu, cauchy_integral, cauchy, fluidity_2, cauchy_times_fluidity_2)
    implicit none

    ! Inputs
    integer, intent(in) :: jmax, number_of_layers
    real*8, intent(in) :: delta_t, eta, mu
    complex*16, intent(inout) :: cauchy_integral(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    complex*16, intent(in) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    complex*16, intent(in) :: fluidity_2(number_of_layers, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: cauchy_times_fluidity_2(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)

    ! Local variables
    integer :: k, max_tensor_index, max_scalar_index, i
    complex*16, allocatable :: cajm(:), ctjml(:), cjml(:)
    real :: start_time, end_time, elapsed_time

    ! Constants
    max_tensor_index = 5 * (jmax * (jmax + 1) / 2 + jmax) - 3
    max_scalar_index = (jmax + 1) * (jmax + 2) / 2

    ! Allocate arrays dynamically
    allocate(cjml(2*max_tensor_index))
    allocate(cajm(2*max_scalar_index))
    allocate(ctjml(2*max_tensor_index))

    ! Initialize arrays to debug-friendly values
    cajm = (0.0d0, 0.0d0)  ! Non-zero imaginary part for debug visibility
    ctjml = (0.0d0, 0.0d0)  ! Non-zero real part for debug visibility
    cjml = (0.0d0, 0.0d0)   ! Non-zero values for both parts

    call CPU_TIME(start_time)

    !Main loop
    do k = 1, number_of_layers

        ! Fill cajm with fluidity values
        cajm(:max_scalar_index) = fluidity_2(k, :max_scalar_index)

        ! Copy cauchy values for the current layer to ctjml
        ctjml(:max_tensor_index) = cauchy(k, :max_tensor_index)

        ! Call VCST with the relevant slice directly
        call VCST(jmax, cajm, ctjml, cjml)

        ! Write calculated data to the cauchy_times_fluidity_2 tensor

        cauchy_times_fluidity_2(k, :) = cjml

        ! Update cauchy_integral using array operations
        cauchy_integral(k, :) = cauchy_integral(k, :) - delta_t * cjml

    end do

    call CPU_TIME(end_time)

    elapsed_time = end_time - start_time

    !print*, "Elapsed time inside the funcion :", elapsed_time, "seconds"

    deallocate(cajm, ctjml, cjml)

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

subroutine calculate_global_dissipation(number_of_layers, jmax, number_of_time_steps_for_deformaion, t, radius, delta_r, cauchy_times_fluidity_2, cauchy, Q_in_time)
    implicit none

    integer, intent(in) :: number_of_layers, jmax, number_of_time_steps_for_deformaion, t
    real*8 :: radius, delta_r
    complex*16, intent(in) :: cauchy_times_fluidity_2(number_of_layers,5*(jmax*(jmax+1)/2+jmax)-3), cauchy(number_of_layers,5*(jmax*(jmax+1)/2+jmax)-3)
    real*8, intent(inout) :: Q_in_time(number_of_time_steps_for_deformaion)

    integer :: i, j, m
    real*8 :: Qcum, fac
    complex*16 :: t_comp, t_times_f_comp
    real*8 :: Q(number_of_layers)

    do i = 1, number_of_layers
        Qcum = 0d0

        j=1
        do m=0, 1
            fac = 2d0
            if (m .eq. 0) fac = 1d0  ! Set fac=1d0 only for m=0

            t_comp = cauchy(i, 3 * ((j * (j + 1)) / 2 + m)-1)
            t_times_f_comp = cauchy_times_fluidity_2(i, 3 * ((j * (j + 1)) / 2 + m)-1)
            Qcum = Qcum + fac * (dreal(t_comp)*dreal(t_times_f_comp) + dimag(t_comp)*dimag(t_times_f_comp))

            t_comp = cauchy(i, 3 * ((j * (j + 1)) / 2 + m)+1)
            t_times_f_comp = cauchy_times_fluidity_2(i, 3 * ((j * (j + 1)) / 2 + m)+1)

            Qcum = Qcum + fac * (dreal(t_comp)*dreal(t_times_f_comp) + dimag(t_comp)*dimag(t_times_f_comp))

        end do
            
    
        do j = 2, jmax
            do m = 0, j
                fac = 2d0
                if (m .eq. 0) fac = 1d0  ! Set fac=1d0 only for m=0
                
                ! Accumulate energy components from cauchy tensor
                t_comp = cauchy(i, 5 * (j * (j + 1) / 2 + m) - 3)
                t_times_f_comp = cauchy_times_fluidity_2(i, 5 * (j * (j + 1) / 2 + m) - 3)

                Qcum = Qcum + fac * (dreal(t_comp)*dreal(t_times_f_comp) + dimag(t_comp)*dimag(t_times_f_comp))

                t_comp = cauchy(i, 5 * (j * (j + 1) / 2 + m) - 4)
                t_times_f_comp = cauchy_times_fluidity_2(i, 5 * (j * (j + 1) / 2 + m) - 4)

                Qcum = Qcum + fac * (dreal(t_comp)*dreal(t_times_f_comp) + dimag(t_comp)*dimag(t_times_f_comp))

                t_comp = cauchy(i, 5 * (j * (j + 1) / 2 + m) - 5)
                t_times_f_comp = cauchy_times_fluidity_2(i, 5 * (j * (j + 1) / 2 + m) - 5)

                Qcum = Qcum + fac * (dreal(t_comp)*dreal(t_times_f_comp) + dimag(t_comp)*dimag(t_times_f_comp))

                t_comp = cauchy(i, 5 * (j * (j + 1) / 2 + m) - 6)
                t_times_f_comp = cauchy_times_fluidity_2(i, 5 * (j * (j + 1) / 2 + m) - 6)

                Qcum = Qcum + fac * (dreal(t_comp)*dreal(t_times_f_comp) + dimag(t_comp)*dimag(t_times_f_comp))

                t_comp = cauchy(i, 5 * (j * (j + 1) / 2 + m) - 7)
                t_times_f_comp = cauchy_times_fluidity_2(i, 5 * (j * (j + 1) / 2 + m) - 7)

                Qcum = Qcum + fac * (dreal(t_comp)*dreal(t_times_f_comp) + dimag(t_comp)*dimag(t_times_f_comp))

            end do
        end do
    
        Q(i) = Qcum
    end do
    
    ! Apply the trapezoidal rule to approximate the integral
    do m = 1, number_of_layers - 1
        Q_in_time(t) = Q_in_time(t) + ((Q(m) * (radius + (m - 1) * delta_r)**2 + Q(m + 1) * (radius + m * delta_r)**2) / 2.0d0) * delta_r
    end do
    
    ! Output the result
    print*, t, Q_in_time(t)

end subroutine

subroutine update_grid_dissipation(number_of_layers, jmax, number_of_time_steps_for_deformaion, cauchy, fluidity_2, averaged_dissipation_on_grid)
    implicit none

    integer, intent(in) :: number_of_layers, jmax, number_of_time_steps_for_deformaion
    complex*16, intent(in) :: cauchy(number_of_layers, 5*(jmax*(jmax+1)/2+jmax)-3)
    complex*16, intent(in) :: fluidity_2(number_of_layers, (jmax + 1) * (jmax + 2) / 2)
    real*8, intent(inout) :: averaged_dissipation_on_grid(number_of_layers, 360, 180)

    integer :: max_scalar_index_2, max_scalar_index, jmax12, jmax3     ! The maximal indexes for harmonic series of a scalar, vector and a deviatoric tensor 
    integer :: i, k, j, m, base_idx, INDT, jm
    real*8 :: multiplicator

    ! Declare arrays as allocatable
    complex*16, allocatable :: cauchy_auxiliary(:)
    complex*16 :: cummulated_series((jmax + 1) * (jmax + 2) / 2)

    complex*16, allocatable :: xx(:), xy(:), xz(:)
    complex*16, allocatable :: yx(:), yy(:), yz(:)
    complex*16, allocatable :: zx(:), zy(:), zz(:)

    real*8, allocatable :: data1(:,:), cummulated_data(:,:)

    complex*16 coef(12000), crhs(12000)

    complex*16 sh1(100000),sh2(100000),sh3(100000)

    allocate(data1(360, 180), cummulated_data(360, 180))

    ! Compute max indices
    max_scalar_index_2 = (jmax + 3) * (jmax + 4) / 2

    max_scalar_index = (jmax + 1) * (jmax + 2) / 2

    ! Now allocate the arrays dynamically
    allocate(cauchy_auxiliary(INDT(jmax, jmax, jmax + 2, 2)))

    allocate(xx(max_scalar_index_2), xy(max_scalar_index_2), xz(max_scalar_index_2))
    allocate(yx(max_scalar_index_2), yy(max_scalar_index_2), yz(max_scalar_index_2))
    allocate(zx(max_scalar_index_2), zy(max_scalar_index_2), zz(max_scalar_index_2))

    jmax12 = jmax
    jmax3 = jmax

    do i=1, number_of_layers

        cauchy_auxiliary = 0.D0
        data1 = 0.D0
        cummulated_data = 0.D0
        coef = 0
        crhs = 0
        cummulated_series = 0

        j=1

        do m=0, j

            if (m==0) then
                multiplicator = 1.0
            else
                multiplicator = 2.0
            end if

            base_idx = 3 * ((j * (j + 1)) / 2 + m)

            cauchy_auxiliary(INDT(j,m,j,2)) = multiplicator*cauchy(i, base_idx-1)
            cauchy_auxiliary(INDT(j,m,j+2,2)) = multiplicator*cauchy(i, base_idx+1)

        end do

        do j=2, jmax

            do m=0, j

                if (m==0) then
                    multiplicator = 1.0
                else
                    multiplicator = 2.0
                end if

                base_idx = 5 * ((j * (j + 1)) / 2 + m)

                cauchy_auxiliary(INDT(j,m,j-2,2)) = multiplicator*cauchy(i, base_idx-7)
                cauchy_auxiliary(INDT(j,m,j-1,2)) = multiplicator*cauchy(i, base_idx-6)
                cauchy_auxiliary(INDT(j,m,j,2)) = multiplicator*cauchy(i, base_idx-5)
                cauchy_auxiliary(INDT(j,m,j+1,2)) = multiplicator*cauchy(i, base_idx-4)
                cauchy_auxiliary(INDT(j,m,j+2,2)) = multiplicator*cauchy(i, base_idx-3)

            end do

        end do

        call shz_tens(jmax, cauchy_auxiliary, xx, xy, xz, yx, yy, yz, zx, zy, zz)

        do j=0, jmax

            do m=1, j

                k=j*(j+1)/2+m+1

                xx(k) = xx(k)/2.0
                yy(k) = yy(k)/2.0
                zz(k) = zz(k)/2.0
                xy(k) = xy(k)/2.0
                xz(k) = xz(k)/2.0
                yz(k) = yz(k)/2.0

            end do

        end do

        call HARMSY(180, jmax, xx, data1)

        do j=1, 180

            do k=1, 360

                cummulated_data(k, j) = cummulated_data(k, j) + data1(k, j)*data1(k, j)

            end do

        end do

        call HARMSY(180, jmax, yy, data1)

        do j=1, 180

            do k=1, 360

                cummulated_data(k, j) = cummulated_data(k, j) + data1(k, j)*data1(k, j)

            end do

        end do

        call HARMSY(180, jmax, zz, data1)

        do j=1, 180

            do k=1, 360

                cummulated_data(k, j) = cummulated_data(k, j) + data1(k, j)*data1(k, j)

            end do

        end do

        call HARMSY(180, jmax, xy, data1)

        do j=1, 180

            do k=1, 360

                cummulated_data(k, j) = cummulated_data(k, j) + 2.0*data1(k, j)*data1(k, j)

            end do

        end do

        call HARMSY(180, jmax, xz, data1)

        do j=1, 180

            do k=1, 360

                cummulated_data(k, j) = cummulated_data(k, j) + 2.0*data1(k, j)*data1(k, j)

            end do

        end do

        call HARMSY(180, jmax, yz, data1)

        do j=1, 180

            do k=1, 360

                cummulated_data(k, j) = cummulated_data(k, j) + 2.0*data1(k, j)*data1(k, j)

            end do

        end do
        

        call HARMSY(180, jmax, fluidity_2(i,:), data1)

        do j=1, 180

            do k=1, 360

                averaged_dissipation_on_grid(i, k, j) = averaged_dissipation_on_grid(i, k, j) + cummulated_data(k, j)*data1(k,j)*(1.0/(number_of_time_steps_for_deformaion-100))

            end do

        end do

    end do


end subroutine


subroutine initial_explicit_euler_update(number_of_layers, jmax, heat_equation_delta_t, radius, delta_r, k_0, ice_density, c_p, temperature, log_temperature, dissipation, upper_dirichlet_temperature, lower_dirichlet_temperature)
    implicit none

    integer, intent(in) :: number_of_layers, jmax
    real*8, intent(in) :: heat_equation_delta_t, radius, delta_r, k_0, ice_density, c_p
    complex*16, intent(inout) :: temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: log_temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(in) :: dissipation(number_of_layers, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(in) :: upper_dirichlet_temperature((jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(in) :: lower_dirichlet_temperature((jmax + 1) * (jmax + 2) / 2)

    integer :: i, j, m, base_idx
    complex*16 :: second_derivative, first_derivative, third_term, j1


    do i=2, number_of_layers

        do j=0, jmax

            j1 = real(j)

            do m=0, j

                base_idx = j*(j+1)/2+m+1

                second_derivative = (1.0/(delta_r*delta_r))*(log_temperature(i+1, base_idx)-2*log_temperature(i, base_idx)+log_temperature(i-1, base_idx))
                first_derivative = (2.0/((radius+(i-(3.0)/(2.0))*delta_r)*2.0*delta_r))*(log_temperature(i+1, base_idx) - log_temperature(i-1, base_idx))
                third_term = - (j1*(j1+1))/((radius+(i-(3.0)/(2.0))*delta_r)*(radius+(i-(3.0)/(2.0))*delta_r))*log_temperature(i, base_idx)
                temperature(i, base_idx) = temperature(i, base_idx) + heat_equation_delta_t*((k_0)/(ice_density*c_p))*(second_derivative+first_derivative+third_term) + heat_equation_delta_t*(1.0/(ice_density*c_p))*(dissipation(i, base_idx)+dissipation(i-1, base_idx))/(2.0)

            end do

        end do

    end do

    do j =0, jmax

        do m=0,j

            base_idx = j*(j+1)/2+m+1

            temperature(1,i) = 2.0*lower_dirichlet_temperature(i)-temperature(2,i)
            temperature(number_of_layers+1,i) = 2.0*upper_dirichlet_temperature(i)-temperature(number_of_layers,i)

        end do

    end do

end subroutine

subroutine update_temperature(number_of_layers, jmax, heat_equation_delta_t, radius, delta_r, k_0, ice_density, c_p, temperature, log_temperature, log_temperature_previous, dissipation, upper_dirichlet_temperature, lower_dirichlet_temperature)
    implicit none

    integer, intent(in) :: number_of_layers, jmax
    real*8, intent(in) :: heat_equation_delta_t, radius, delta_r, k_0, ice_density, c_p
    complex*16, intent(inout) :: temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: log_temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: log_temperature_previous(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(in) :: upper_dirichlet_temperature((jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(in) :: lower_dirichlet_temperature((jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(in) :: dissipation(number_of_layers, (jmax + 1) * (jmax + 2) / 2)

    integer :: i, j, m, base_idx
    complex*16 :: second_derivative, first_derivative, third_term, second_derivative_, first_derivative_, third_term_, first_part, second_part
    real*8 :: j1

    do i=2, number_of_layers

        do j=0, jmax

            j1 = real(j)

            do m=0, j

                base_idx = j*(j+1)/2+m+1

                second_derivative = (1.0/(delta_r*delta_r))*(log_temperature(i+1, base_idx)-2*log_temperature(i, base_idx)+log_temperature(i-1, base_idx))
                first_derivative = (2.0/((radius+(i-(3.0)/(2.0))*delta_r)*2.0*delta_r))*(log_temperature(i+1, base_idx) - log_temperature(i-1, base_idx))
                third_term = - (j1*(j1+1))/((radius+(i-(3.0)/(2.0))*delta_r)*(radius+(i-(3.0)/(2.0))*delta_r))*log_temperature(i, base_idx)
                first_part = ((k_0)/(ice_density*c_p))*(second_derivative+first_derivative+third_term)
                
                second_derivative_ = (1.0/(delta_r*delta_r))*(log_temperature_previous(i+1, base_idx)-2*log_temperature_previous(i, base_idx)+log_temperature_previous(i-1, base_idx))
                first_derivative_ = (2.0/((radius+(i-(3.0)/(2.0))*delta_r)*2.0*delta_r))*(log_temperature_previous(i+1, base_idx) - log_temperature_previous(i-1, base_idx))
                third_term_ = - (j1*(j1+1))/((radius+(i-(3.0)/(2.0))*delta_r)*(radius+(i-(3.0)/(2.0))*delta_r))*log_temperature_previous(i, base_idx)
                second_part = ((k_0)/(ice_density*c_p))*(second_derivative+first_derivative+third_term)

                temperature(i, base_idx) = temperature(i, base_idx) + heat_equation_delta_t*((3.0/2.0)*first_part - (1.0/2.0)*second_part) + heat_equation_delta_t*(1.0/(ice_density*c_p))*(dissipation(i, base_idx)+dissipation(i-1, base_idx))/(2.0)

            end do

        end do

    end do

    do j =0, jmax

        do m=0,j

            base_idx = j*(j+1)/2+m+1

            temperature(1,i) = 2.0*lower_dirichlet_temperature(i)-temperature(2,i)
            temperature(number_of_layers+1,i) = 2.0*upper_dirichlet_temperature(i)-temperature(number_of_layers,i)

        end do

    end do

end subroutine

subroutine update_log_temperature(number_of_layers, jmax, temperature, log_temperature, log_temperature_previous)
    implicit none

    integer, intent(in) :: number_of_layers, jmax
    complex*16, intent(inout) :: temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: log_temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: log_temperature_previous(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)

    integer :: i, k, l, m, max_scalar_index
    
    real*8, allocatable :: temperature_data(:,:), log_temperature_data(:,:)
    complex*16 coef(12000),crhs(12000),coef1(12000)
    allocate(temperature_data(360, 180), log_temperature_data(360, 180))

    max_scalar_index = (jmax + 1) * (jmax + 2) / 2


    do i=1, number_of_layers+1

        do k=1, max_scalar_index

            log_temperature_previous(i,k) = log_temperature(i,k)

        end do


        call HARMSY(180, jmax, temperature(i,:), temperature_data)

        do l=1, 180

            do m=1, 360

                log_temperature_data(m, l) = DLOG(temperature_data(m, l))
            end do

        end do

        call HARMAN(180,jmax,log_temperature_data,crhs)
        call HARMLS(180,jmax,crhs,coef)

        do k=1, max_scalar_index

            log_temperature(i,k) = coef(k)

        end do

    end do

end subroutine

subroutine write_temperature_data(number_of_layers, jmax, t, k, number_of_steps_for_heat_equation, heat_equation_delta_t, delta_r, thickness, &
    temperature, folder_name, num_angles, angles)
    implicit none

    ! Inputs
    integer, intent(in) :: number_of_layers, jmax, t, k, number_of_steps_for_heat_equation
    real*8, intent(in) :: heat_equation_delta_t, delta_r, thickness
    complex*16, intent(inout) :: temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    character(len=*), intent(in) :: folder_name  ! Folder path where files will be saved
    integer, intent(in) :: num_angles  ! Number of predefined angles
    integer, dimension(num_angles,2), intent(in) :: angles  ! Array of (latitude, longitude) pairs

    ! Local variables
    real*8, allocatable :: data(:,:)
    integer :: i, layer_index, file_unit, l, m
    character(len=100) :: file_path

    ! Allocate memory for data
    allocate(data(360, 180))

    ! ===== 1. Write temperature at layers (5, 10, ..., 45) =====
    do i = 1, 9
        layer_index = i * 5  ! 5, 10, 15, ..., 45
        write(file_path, '(A,"/temperature",I1,".txt")') trim(folder_name), i

        ! Process temperature data
        call HARMSY(180, jmax, temperature(layer_index,:), data)

        ! Assign a unique file unit (11-19)
        file_unit = 10 + i  

        ! Open file in append mode
        open(unit=file_unit, file=trim(file_path), status="old", action="write", position="append")

        ! Write time information
        write(file_unit, *) (k-1) * number_of_steps_for_heat_equation * heat_equation_delta_t / (31536E9)

        ! Write spatial data
        do l = 1, 180
            do m = 1, 360
                write(file_unit, *) 180-l + 0.5, m-1, data(m,l)
            end do
        end do

        ! Close file
        close(file_unit)
    end do

        ! ===== 2. Write temperature at predefined angles =====
    do i = 1, num_angles
        ! Generate filename
        write(file_path, '(A,"/temperature_",I3.3,"_",I3.3,".txt")') &
        trim(folder_name), angles(i,1), angles(i,2)

        ! Assign a unique file unit (21-26, or more depending on num_angles)
        file_unit = 20 + i  

        ! Open file in append mode
        open(unit=file_unit, file=trim(file_path), status="old", action="write", position="append")

        ! Write time information
        write(file_unit, *) (k-1) * number_of_steps_for_heat_equation * heat_equation_delta_t / (31536E9)

        call HARMSY(180, jmax, (temperature(1,:)+temperature(2,:))/(2.0), data)
        write(file_unit, *) 0, &
            data(angles(i,2)+1, 180-angles(i,1))

        ! Loop over layers to write temperature at specific lat/lon
        do layer_index = 2, number_of_layers
            call HARMSY(180, jmax, temperature(layer_index,:), data)

            ! Select temperature at the corresponding (lat, lon)
            write(file_unit, *) (layer_index-(3.0)/(2.0))*delta_r/(1000.0), &
            data(angles(i,2)+1, 180-angles(i,1))
        end do

        call HARMSY(180, jmax, (temperature(number_of_layers,:)+temperature(number_of_layers+1,:))/(2.0), data)
        write(file_unit, *) (thickness)/(1000.0), &
            data(angles(i,2)+1, 180-angles(i,1))

        ! Close file
        close(file_unit)
    end do

    ! Deallocate memory
    deallocate(data)

end subroutine write_temperature_data


subroutine write_dissipation(number_of_layers, jmax, k, number_of_steps_for_heat_equation, heat_equation_delta_t, dissipation, Total_averaged_dissipation, folder_name)
    implicit none

    ! Inputs
    integer, intent(in) :: number_of_layers, jmax, k, number_of_steps_for_heat_equation
    real*8, intent(in) :: heat_equation_delta_t
    complex*16, intent(inout) :: dissipation(number_of_layers, (jmax + 1) * (jmax + 2) / 2)
    real*8, intent(in) :: Total_averaged_dissipation
    character(len=*), intent(in) :: folder_name  ! Folder path where files will be saved

    ! Local variables
    real*8, allocatable :: data(:,:)
    integer :: i, layer_index, file_unit, l, m
    character(len=100) :: file_path

    ! Allocate memory for data
    allocate(data(360, 180))

    ! ===== 1. Write dissipation at layers (5, 10, ..., 45) =====
    do i = 1, 9
        layer_index = i * 5  ! 5, 10, 15, ..., 45
        write(file_path, '(A,"/dissipation",I1,".txt")') trim(folder_name), i

        ! Process dissipation data
        call HARMSY(180, jmax, dissipation(layer_index,:), data)

        ! Assign a unique file unit (14-22)
        file_unit = 13 + i  

        ! Open file in append mode
        open(unit=file_unit, file=trim(file_path), status="old", action="write", position="append")

        ! Write time information
        write(file_unit, *) (k-1) * number_of_steps_for_heat_equation * heat_equation_delta_t / (31536E9)

        ! Write spatial data
        do l = 1, 180
            do m = 1, 360
                write(file_unit, *) 180-l + 0.5, m-1, data(m,l)
            end do
        end do

        ! Close file
        close(file_unit)
    end do

    ! ===== 2. Write total dissipation value =====
    write(file_path, '(A,"/total_dissipation.txt")') trim(folder_name)

    ! Open total dissipation file in append mode
    open(unit=23, file=trim(file_path), status="old", action="write", position="append")

    ! Write total dissipation value
    write(23, *) (k-1) * number_of_steps_for_heat_equation * heat_equation_delta_t / (31536E9), &
                 Total_averaged_dissipation / (1.0E9)

    ! Close total dissipation file
    close(23)

    ! Deallocate memory
    deallocate(data)

end subroutine write_dissipation


subroutine calculate_fluidity_2(number_of_layers, jmax, eta_0, temperature, fluidity_2)
    implicit none

    integer, intent(in) :: number_of_layers, jmax
    real*8, intent(in) :: eta_0
    complex*16, intent(in) :: temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: fluidity_2(number_of_layers, (jmax + 1) * (jmax + 2) / 2)

    real*8, allocatable :: temperature_data(:,:), fluidity_2_data(:,:)
    integer :: l, m, i, max_scalar_index
    complex*16 coef(12000),crhs(12000),coef1(12000)

    max_scalar_index = (jmax + 1) * (jmax + 2) / 2

    allocate(temperature_data(360, 180), fluidity_2_data(360, 180))


    do i=1, number_of_layers

        call HARMSY(180, jmax, (temperature(i,:)+temperature(i+1,:))/2.0, temperature_data)

        do l=1, 180

            do m=1, 360
    
                fluidity_2_data(m,l) = 1.0/(2.0*eta_0*dexp(((59000.0)/(8.314*273.0))*(((273.0)/(temperature_data(m,l)))-1.0)))
    
            end do
    
        end do

        call HARMAN(180,jmax,fluidity_2_data,crhs)
        call HARMLS(180,jmax,crhs,coef)

        do l=1, max_scalar_index

            fluidity_2(i,l) = coef(l)

        end do

    end do

end subroutine


subroutine write_fluidity_2_data(number_of_layers, jmax, t, k, number_of_steps_for_heat_equation, heat_equation_delta_t, fluidity_2, folder_name)
    implicit none

    ! Inputs
    integer, intent(in) :: number_of_layers, jmax, t, k, number_of_steps_for_heat_equation
    real*8, intent(in) :: heat_equation_delta_t
    complex*16, intent(inout) :: fluidity_2(number_of_layers, (jmax + 1) * (jmax + 2) / 2)
    character(len=*), intent(in) :: folder_name  ! Folder path where files will be saved

    ! Local variables
    real*8, allocatable :: data(:,:)
    integer :: i, layer_index, file_unit, l, m
    character(len=100) :: file_path

    ! Allocate memory for data
    allocate(data(360, 180))

    ! Loop over indices (5, 10, 15, ..., 45) and write corresponding files
    do i = 1, 9
        layer_index = i * 5  ! 5, 10, 15, ..., 45
        write(file_path, '(A,I1,A)') trim(folder_name)//"/fluidity_2", i, ".txt"

        ! Process the data for the given layer index
        call HARMSY(180, jmax, fluidity_2(layer_index,:), data)

        ! Assign a unique unit number for each file
        file_unit = 30 + i  ! Ensures units 31-39 for fluidity files

        ! Open the file in append mode and write data
        open(unit=file_unit, file=trim(file_path), status="old", action="write", position="append")
        
        ! Write time information
        write(file_unit, *) (k-1) * number_of_steps_for_heat_equation * heat_equation_delta_t / (31536E9)

        ! Write spatial data
        do l = 1, 180
            do m = 1, 360
                write(file_unit, *) 180-l + 0.5, m-1, data(m,l)
            end do
        end do

        ! Close file
        close(file_unit)
    end do


    ! Deallocate memory
    deallocate(data)

end subroutine write_fluidity_2_data


subroutine set_up_initial_values_for_temperature_and_log_temperature(number_of_layers, jmax, temperature, log_temperature, log_temperature_previous, upper_dirichlet_temperature, lower_dirichlet_temperature)
    implicit none

    ! Inputs
    integer, intent(in) :: number_of_layers, jmax
    complex*16, intent(inout) :: temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: log_temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: log_temperature_previous(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: upper_dirichlet_temperature((jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: lower_dirichlet_temperature((jmax + 1) * (jmax + 2) / 2)

    ! Local variables
    real*8, allocatable :: temperature_data(:,:), log_temperature_data(:,:)
    integer :: j, m, i, k, max_scalar_index, l
    integer :: q, w
    complex*16 :: coef_, delta_coef
    character(len=100) :: file_path_1, file_path_2, file_path_3
    complex*16 coef(12000),crhs(12000),coef1(12000)

    max_scalar_index = (jmax + 1) * (jmax + 2) / 2

    ! Allocate memory for data
    allocate(temperature_data(360, 180), log_temperature_data(360, 180))

    open(1,file='outer_dirichlet.dat') ! cteni SH koeficientu ze souboru

    do j=0,30
        do m=0,j
            read(1,*) q, w, coef_
            i=j*(j+1)/2+m+1
            upper_dirichlet_temperature(i) = coef_
        end do
    end do

    close(1)
    
    lower_dirichlet_temperature(1) = 967.7

    do j=0,30 ! cyklus pres stupen
        do m=0,j    ! cyklus pres rad
            i=j*(j+1)/2+m+1

            coef_ = upper_dirichlet_temperature(i)

            if (j==0) then
                delta_coef = (coef_ - lower_dirichlet_temperature(1))/(number_of_layers-1)
            else
                delta_coef = (coef_)/(number_of_layers-1)
            end if

            do l=1, number_of_layers+1

                temperature(l,i) = lower_dirichlet_temperature(i) + (l-(3.0)/(2.0))*delta_coef

            end do
        enddo
    enddo

    do i=1, number_of_layers+1

        call HARMSY(180, jmax, temperature(i,:), temperature_data)

        do j=1, 180

            do m=1, 360

                log_temperature_data(m, j) = DLOG(temperature_data(m, j))

            end do

        end do
        
        call HARMAN(180,jmax,log_temperature_data,crhs)
        call HARMLS(180,jmax,crhs,coef)

        do k=1, max_scalar_index

            log_temperature(i,k) = coef(k)
            log_temperature_previous(i,k) = coef(k)


        end do

    end do

end subroutine set_up_initial_values_for_temperature_and_log_temperature

subroutine update_cauchy_times_fluidity_2(jmax, number_of_layers, delta_t, eta, mu, cauchy_integral, cauchy, fluidity_2, cauchy_times_fluidity_2)
    implicit none

    ! Inputs
    integer, intent(in) :: jmax, number_of_layers
    real*8, intent(in) :: delta_t, eta, mu
    complex*16, intent(inout) :: cauchy_integral(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    complex*16, intent(in) :: cauchy(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)
    complex*16, intent(in) :: fluidity_2(number_of_layers, (jmax + 1) * (jmax + 2) / 2)
    complex*16, intent(inout) :: cauchy_times_fluidity_2(number_of_layers, 5 * (jmax * (jmax + 1) / 2 + jmax) - 3)

    ! Local variables
    integer :: k, max_tensor_index, max_scalar_index, i
    complex*16, allocatable :: cajm(:), ctjml(:), cjml(:)
    real :: start_time, end_time, elapsed_time

    ! Constants
    max_tensor_index = 5 * (jmax * (jmax + 1) / 2 + jmax) - 3
    max_scalar_index = (jmax + 1) * (jmax + 2) / 2

    ! Allocate arrays dynamically
    allocate(cjml(2*max_tensor_index))
    allocate(cajm(2*max_scalar_index))
    allocate(ctjml(2*max_tensor_index))

    ! Initialize arrays to debug-friendly values
    cajm = (0.0d0, 0.0d0)  ! Non-zero imaginary part for debug visibility
    ctjml = (0.0d0, 0.0d0)  ! Non-zero real part for debug visibility
    cjml = (0.0d0, 0.0d0)   ! Non-zero values for both parts

    call CPU_TIME(start_time)

    !Main loop
    do k = 1, number_of_layers

        ! Fill cajm with fluidity values
        cajm(:max_scalar_index) = fluidity_2(k, :max_scalar_index)

        ! Copy cauchy values for the current layer to ctjml
        ctjml(:max_tensor_index) = cauchy(k, :max_tensor_index)

        ! Call VCST with the relevant slice directly
        call VCST(jmax, cajm, ctjml, cjml)

        ! Write calculated data to the cauchy_times_fluidity_2 tensor

        cauchy_times_fluidity_2(k, :) = cjml

    end do

    call CPU_TIME(end_time)

    elapsed_time = end_time - start_time

    !print*, "Elapsed time inside the funcion :", elapsed_time, "seconds"

    deallocate(cajm, ctjml, cjml)

end subroutine

subroutine dissipation_from_grid_and_total(number_of_layers, jmax, radius, delta_r, averaged_dissipation_on_grid, dissipation, Total_averaged_dissipation)
    implicit none

    integer, intent(in) :: number_of_layers, jmax
    real*8 :: radius, delta_r
    real*8 averaged_dissipation_on_grid(number_of_layers,360, 180)
    complex*16 :: dissipation(number_of_layers, (jmax + 1) * (jmax + 2) / 2)
    real*8, intent(inout) :: Total_averaged_dissipation

    complex*16 :: sqrt_dissipation(number_of_layers, (jmax + 1) * (jmax + 2) / 2)
    real*8 :: sqrt_on_grid(360, 180), data1(360,180)

    complex*16 coef(12000),crhs(12000),coef1(12000)

    integer :: i, j, m, k, o, l
    real*8 :: Qcum, fac
    real*8 :: Q(number_of_layers)
    real*8 :: Q_total
    real*8 :: theta, weight, pi

    do i=1, number_of_layers

        call HARMAN(180,jmax,averaged_dissipation_on_grid(i,:,:),crhs)
        call HARMLS(180,jmax,crhs,coef)

        do j=0,jmax
            do m=0,j
            o=j*(j+1)/2+m+1
            dissipation(i,o) = coef(o)
            enddo
        enddo

    end do

    Q_total=0

    pi = 4.0d0 * datan(1.0d0)

    ! Loop over radial layers
    do i = 1, number_of_layers
        Qcum = 0.0d0

        call HARMSY(180, jmax, dissipation(i,:), data1)

        ! Loop over latitude and longitude
        do k = 1, 180
            theta = (k - 0.5) * (pi / 180.0)  ! Midpoint latitude in radians
            weight = dsin(theta)  ! Proper weight factor for latitude integration

            do l = 1, 360
                Qcum = Qcum + averaged_dissipation_on_grid(i, l, k) * &
                    (pi / 180.0) * (2.0d0 * pi / 360.0) * weight
            end do
        end do

        Q(i) = Qcum
    end do

    ! Apply the trapezoidal rule to approximate the integral
    do m = 1, number_of_layers - 1
        Q_total = Q_total + ((Q(m) * (radius + (m - 1) * delta_r)**2 + Q(m + 1) * (radius + m * delta_r)**2) / 2.0d0) * delta_r
    end do

    Total_averaged_dissipation = Q_total

    ! Output the result
    print*, 'The total dissipation from integral of averaged dissipation', Q_total/(1.0E9)

end subroutine

subroutine calculate_and_write_heat_flux(jmax, number_of_layers, k, number_of_steps_for_heat_equation, heat_equation_delta_t, radius, delta_r, thickness, k_0, temperature,  folder_name)
    implicit none

    ! Inputs
    integer, intent(in) :: jmax, number_of_layers, k, number_of_steps_for_heat_equation
    real*8, intent(in) :: heat_equation_delta_t,radius, delta_r, thickness, k_0
    complex*16, intent(in) :: temperature(number_of_layers+1, (jmax + 1) * (jmax + 2) / 2)
    character(len=*), intent(in) :: folder_name  ! Folder path where files will be saved

    complex*16 :: lower_temperature_gradient(3 * (jmax * (jmax + 1) / 2 + jmax) + 1), upper_temperature_gradient(3 * (jmax * (jmax + 1) / 2 + jmax) + 1)
    complex*16 :: lower_radial_gradiant((jmax + 1) * (jmax + 2) / 2), upper_radial_gradiant((jmax + 1) * (jmax + 2) / 2)
    integer :: vector_index, scalar_index, j, m, l
    real*8 :: j1
    real*8, allocatable :: temperature_data(:,:), readial_gradient_data(:,:)
    real*8 :: Qcum, pi, weight, theta

    pi = 4.0d0 * datan(1.0d0)

    allocate(temperature_data(360, 180), readial_gradient_data(360, 180))

    lower_temperature_gradient=0D0
    upper_temperature_gradient=0D0

    do j=0,jmax

        j1 = real(j)

        do m=0, j

            scalar_index = j * (j + 1) / 2 + m + 1
            vector_index = 3 * (j * (j + 1) / 2 + m)

            if (j==0) then

                lower_temperature_gradient(vector_index+1)= - (dsqrt(j1+1))/(dsqrt(2*j1+1))*((Temperature(2,scalar_index)-Temperature(1,scalar_index))/(delta_r)-(j1)/(radius)*(Temperature(2,scalar_index)+Temperature(1,scalar_index))/(2.0))

                upper_temperature_gradient(vector_index+1)= - (dsqrt(j1+1))/(dsqrt(2*j1+1))*((Temperature(number_of_layers+1,scalar_index)-Temperature(number_of_layers,scalar_index))/(delta_r)-(j1)/(radius+thickness)*(Temperature(number_of_layers+1,scalar_index)+Temperature(number_of_layers,scalar_index))/(2.0))

            else

                lower_temperature_gradient(vector_index-1)= (dsqrt(j1))/(dsqrt(2*j1+1))*((Temperature(2,scalar_index)-Temperature(1,scalar_index))/(delta_r)+(j1+1)/(radius)*(Temperature(2,scalar_index)+Temperature(1,scalar_index))/(2.0))

                lower_temperature_gradient(vector_index+1)= - (dsqrt(j1+1))/(dsqrt(2*j1+1))*((Temperature(2,scalar_index)-Temperature(1,scalar_index))/(delta_r)-(j1)/(radius)*(Temperature(2,scalar_index)+Temperature(1,scalar_index))/(2.0))

                upper_temperature_gradient(vector_index-1)= (dsqrt(j1))/(dsqrt(2*j1+1))*((Temperature(number_of_layers+1,scalar_index)-Temperature(number_of_layers,scalar_index))/(delta_r)+(j1+1)/(radius+thickness)*(Temperature(number_of_layers+1,scalar_index)+Temperature(number_of_layers,scalar_index))/(2.0))

                upper_temperature_gradient(vector_index+1)= - (dsqrt(j1+1))/(dsqrt(2*j1+1))*((Temperature(number_of_layers+1,scalar_index)-Temperature(number_of_layers,scalar_index))/(delta_r)-(j1)/(radius+thickness)*(Temperature(number_of_layers+1,scalar_index)+Temperature(number_of_layers,scalar_index))/(2.0))

            end if

        end do

    end do

    do j=0,jmax

        j1 = real(j)

        do m=0, j

            scalar_index = j * (j + 1) / 2 + m + 1
            vector_index = 3 * (j * (j + 1) / 2 + m)

            if (j==0) then

                upper_radial_gradiant(scalar_index) = -dsqrt((j1+1)/(2*j1+1))*upper_temperature_gradient(vector_index+1)
                lower_radial_gradiant(scalar_index) = -dsqrt((j1+1)/(2*j1+1))*lower_temperature_gradient(vector_index+1)

            else

                upper_radial_gradiant(scalar_index) = dsqrt((j1)/(2*j1+1))*upper_temperature_gradient(vector_index-1)-dsqrt((j1+1)/(2*j1+1))*upper_temperature_gradient(vector_index+1)
                lower_radial_gradiant(scalar_index) = dsqrt((j1)/(2*j1+1))*lower_temperature_gradient(vector_index-1)-dsqrt((j1+1)/(2*j1+1))*lower_temperature_gradient(vector_index+1)

            end if

        end do

    end do

    open(unit=60, file=trim(trim(folder_name) // "/lower_heat_flux.txt"), status="old", action="write", position="append")
        
    ! Write time information
    write(60, *) (k-1) * number_of_steps_for_heat_equation * heat_equation_delta_t / (31536E9)

    call HARMSY(180, jmax, lower_radial_gradiant(:), readial_gradient_data)
    call HARMSY(180, jmax, (temperature(1,:)+temperature(2,:))/2.0, temperature_data)


    ! Write spatial data
    do l = 1, 180
        do m = 1, 360
            write(60, *) 180-l + 0.5, m-1, -(k_0*readial_gradient_data(m,l))/(temperature_data(m,l))
        end do
    end do

        ! Close file
    close(60)

    Qcum = 0.0d0
    ! Loop over latitude and longitude
    do l = 1, 180
        theta = (l - 0.5) * (pi / 180.0)  ! Midpoint latitude in radians
        weight = dsin(theta)  ! Proper weight factor for latitude integration

        do m = 1, 360
            Qcum = Qcum + -(k_0*readial_gradient_data(m,l))/(temperature_data(m,l)) * &
                (pi / 180.0) * (2.0d0 * pi / 360.0) * weight
        end do
    end do

    Qcum = Qcum*(radius)*(radius)

    open(unit=62, file=trim(trim(folder_name) // "/total_lower_heat_flux.txt"), status="old", action="write", position="append")

    ! Write total dissipation value
    write(62, *) (k-1) * number_of_steps_for_heat_equation * heat_equation_delta_t / (31536E9), &
    Qcum / (1.0E9)
    close(62)

    print*, 'Total lower heat flux:', Qcum / (1.0E9)

    open(unit=61, file=trim(trim(folder_name) // "/upper_heat_flux.txt"), status="old", action="write", position="append")
        
    ! Write time information
    write(61, *) (k-1) * number_of_steps_for_heat_equation * heat_equation_delta_t / (31536E9)

    call HARMSY(180, jmax, upper_radial_gradiant(:), readial_gradient_data)
    call HARMSY(180, jmax, (temperature(number_of_layers+1,:)+temperature(number_of_layers,:))/2.0, temperature_data)


    ! Write spatial data
    do l = 1, 180
        do m = 1, 360
            write(61, *) 180-l + 0.5, m-1, -(k_0*readial_gradient_data(m,l))/(temperature_data(m,l))
        end do
    end do

    ! Close file
    close(61)


    Qcum = 0.0d0
    ! Loop over latitude and longitude
    do l = 1, 180
        theta = (l - 0.5) * (pi / 180.0)  ! Midpoint latitude in radians
        weight = dsin(theta)  ! Proper weight factor for latitude integration

        do m = 1, 360
            Qcum = Qcum + -(k_0*readial_gradient_data(m,l))/(temperature_data(m,l)) * &
                (pi / 180.0) * (2.0d0 * pi / 360.0) * weight
        end do
    end do

    Qcum = Qcum*(radius+thickness)*(radius+thickness)

    open(unit=63, file=trim(trim(folder_name) // "/total_upper_heat_flux.txt"), status="old", action="write", position="append")

    ! Write total dissipation value
    write(63, *) (k-1) * number_of_steps_for_heat_equation * heat_equation_delta_t / (31536E9), &
    Qcum / (1.0E9)
    close(63)

    print*, 'Total upper heat flux:', Qcum / (1.0E9)

end subroutine calculate_and_write_heat_flux

! This program simulates the deformation of Jupiter's moon Europa
program Europa_simulation
    implicit none   ! Require all variables to be explicitly declared

    ! The resolution parameters of the simulation 
    integer, parameter :: number_of_layers=50                           ! number_of_layers determines the radial resolution
    integer, parameter :: jmax=30                                       ! jmax determines the spherical resolution, that is the maximal degree j used in spherical harmonic series

    ! The length of the discrete time step
    real*8, parameter :: delta_t=3.551181*24*36                         ! This is 0.01 of the period of the Europa's revolution around Jupiter, the period is 3.551181 days
    real*8, parameter :: heat_equation_delta_t=2E10

    integer, parameter :: number_or_total_rounds=200
    ! Number of time steps, determines the length of the simulation
    integer, parameter :: number_of_time_steps_for_deformaion=200
    integer, parameter :: number_of_steps_for_heat_equation=100
    
    ! These are the physical paramteres for Europa's icy layer, the main parameters of the simulation
    real*8, parameter :: radius = 1531000.0d0                           ! The inner radius of Europa in meters
    real*8, parameter :: thickness = 30000.0d0                          ! Thickness of the ice layer in meters
    real*8, parameter :: mu = 3.3E9                                     ! Shear modulus of ice in Pascal
    real*8, parameter :: eta = 1E15                                     ! Viscosity of ice in Pascal*s
    real*8, parameter :: surface_g = 1.314                              ! Surface gravity in m/s^2
    real*8, parameter :: bottom_g = 1.314                               ! Bottom gravity in m/s^2 (same as surface gravity)
    real*8, parameter :: ice_density = 920.0d0                          ! Density of ice in kg/m^3
    real*8, parameter :: delta_rho = 80.0d0                             ! Difference between density of water and ice in kg/m^3
    real*8, parameter :: angular_speed = 2.047827249E-5                 ! Angular speed of the revolution around Jupiter in rad/s
    real*8, parameter :: excentricity = 0.009                           ! Eccentricity of Europa's orbit
    real*8, parameter :: k_0 = 651
    real*8, parameter :: c_p = 2108
    real*8, parameter :: eta_0 = 1E13

    real*8, parameter :: delta_r = thickness / (number_of_layers - 1)   ! Radial distance between layers

    ! Allocatable arrays
    complex*16, allocatable :: displacement(:,:)                        ! The harmonic series of the displacement vector field, the first index is layer index, displacement is calculated between layers, so it has number_of_layers + 1 dimension, second is the vector multiindex
    complex*16, allocatable :: cauchy(:,:)                              ! The harmonic series of the deviatoric part of the cauchy stress tensor field, the first index is layer index, second is the tensor multiindex
    complex*16, allocatable :: cauchy_isotropic(:,:,:)                  ! The harmonic series of the isotropic part of the cauchy stress tensor field, the first index is layer index, second the degree index(!starting at 1), and third the order index(!Also starting at 1), this will be probably changed to scalar multiindex
                                                                        ! The isotropic part of the cauchy stress tensor, it is not really that important in the simulation. It is used mainly in tests of the satisfaction of deformation equations
    complex*16, allocatable :: cauchy_integral(:,:)                     ! This is the viscous time integral of cauchy stress tensor divided by 2 times viscosity, which appears at the right hand side of the constitutive relation, it is updated every time step
    complex*16, allocatable :: fluidity_2(:,:)                          ! The harmonic series of the fluidity scalar field, fluidity is the inverted value of viscosity divided by two, the first index is layer index and the second is the scalar multiindex
    complex*16, allocatable :: cauchy_times_fluidity_2(:,:)             ! This is a deviatoric tensor, which is computed as a product of scalar fludity_2 and the deviatoric cauchy stress tensor, it is needed in the right hand side of the constitutive relation and in the dissipation calculation

    complex*16, allocatable :: volume_force(:,:)                        ! Thetidal deformation force, the first index is the layer index, volume force is calculated in the middle of each layer so it has number_of_layers - 1 dimensions, and the second index has only two dimensions, corresponding to the vector indeces 8 and 14 
    complex*16, allocatable :: bottom_force(:)                          ! The tidal deformation force acting on the bottom of the ice layer, it has four dimensions, corresponding to the vector indeces 8, 10, 14, 16

    real*8, allocatable :: matrix(:,:)                                  ! The matrix of the linear system which is the discretization of the spherical deformation equations for degree j>1
    real*8, allocatable :: matrix_for_j1(:,:)                           ! The matrix of the linear system which is the discretization of the spherical deformation equations for degree j=1
    real*8, allocatable :: toroidial_matrix(:,:)                        ! The matrix of the linear system which is the discretization of the toroidial deformation equations for degree j>1

    complex*16, allocatable :: radial_displacement(:,:)                 ! The harmonic series of the radial displacement at the outer boundary, mainly used for benchmark and debugging purposes, two indeces first order(!starting at 1) and second degree(!starting at one)

    real*8, allocatable :: Q_in_time(:)                                 ! The total dissipation (heat production) in the whole ice crust as a function of time, each dimension is a time step
    real*8 :: Total_averaged_dissipation                                ! This is the total averaged dissipation that is averaged over few periods of revolution of Europa around Jupiter
    complex*16, allocatable :: dissipation(:,:)                         ! The harmonic series of the dissipation/heating, at each layer we have a harmonic series

    complex*16, allocatable :: temperature(:,:)
    complex*16, allocatable :: log_temperature(:,:)
    complex*16, allocatable :: log_temperature_previous(:,:)

    complex*16, allocatable :: upper_dirichlet_temperature(:)
    complex*16, allocatable :: lower_dirichlet_temperature(:)

    complex*16, allocatable :: lower_temperature_gradient(:)
    complex*16, allocatable :: upper_temperature_gradient(:)

    real*8, allocatable :: averaged_dissipation_on_grid(:,:,:)

    ! Loop variables and other results
    integer :: t, j, i, m, k                                            ! Integer variables used in do loops
    real*8 :: j1                                                        ! the real value of the integer degree j used in filling the matrices of linear systems

    ! Used while reading fluidity values from a file
    integer :: q_, w_, e_, o, base_idx                                               ! just some randomm variables used only while loading the coefficinets for the fluidity from a file
    complex*16 :: coef_                                                  ! An auxiliary variable used while loading fluidity coefficients from a file

    integer :: max_scalar_index, max_vector_index, max_tensor_index     ! The maximal indexes for harmonic series of a scalar, vector and a deviatoric tensor 

    logical :: test                                                     ! A boolean variable, whether to test or not test the solution of equations, if true then tests of equation satisfaction are run
    logical :: write_deformation_data                                   ! A boolean variable, whether to write or not some simulated data to a file, if true then data of radial displacement at boundary are written to a file in a folder
    logical :: write_simulation_data

    ! This variables are for writing of data to a file
    integer, parameter :: num_angles = 7
    integer, dimension(num_angles, 2) :: angles
    character(len=80) :: folder_name, file_name
    character(len=8)  :: date
    character(len=10) :: time
    integer, dimension(8) :: values

    complex*16 coef(12000),crhs(12000),coef1(12000)

    real :: start_time, end_time, elapsed_time

    angles(1,1) = 0   ; angles(1,2) = 0
    angles(2,1) = 90  ; angles(2,2) = 0
    angles(3,1) = 90  ; angles(3,2) = 90
    angles(4,1) = 90  ; angles(4,2) = 180
    angles(5,1) = 90  ; angles(5,2) = 270
    angles(6,1) = 179 ; angles(6,2) = 0
    angles(7,1) = 30   ; angles(7,2) = 90

    max_scalar_index = (jmax + 1) * (jmax + 2) / 2
    max_vector_index = 3 * (jmax * (jmax + 1) / 2 + jmax) + 1
    max_tensor_index = 5 * (jmax * (jmax + 1) / 2 + jmax) - 3
    
    ! Allocate arrays
    allocate(displacement(number_of_layers + 1, max_vector_index))
    allocate(cauchy(number_of_layers, max_tensor_index))
    allocate(cauchy_isotropic(number_of_layers, jmax+1, jmax+1))
    allocate(cauchy_integral(number_of_layers, max_tensor_index))
    allocate(fluidity_2(number_of_layers, max_scalar_index))
    allocate(cauchy_times_fluidity_2(number_of_layers, max_tensor_index))

    allocate(volume_force(number_of_layers - 1, 2))
    allocate(bottom_force(4))

    allocate(matrix(6 * number_of_layers + 2, 6 * number_of_layers + 2))
    allocate(matrix_for_j1(5 * number_of_layers + 2, 5 * number_of_layers + 2))
    allocate(toroidial_matrix(3 * number_of_layers + 1, 3 * number_of_layers + 1))

    allocate(radial_displacement(jmax + 1, jmax + 1))

    allocate(Q_in_time(number_of_time_steps_for_deformaion))
    allocate(dissipation(number_of_layers, max_scalar_index))

    allocate(temperature(number_of_layers+1, max_scalar_index))
    allocate(log_temperature(number_of_layers+1, max_scalar_index))
    allocate(log_temperature_previous(number_of_layers+1, max_scalar_index))

    allocate(upper_dirichlet_temperature(max_scalar_index))
    allocate(lower_dirichlet_temperature(max_scalar_index))

    allocate(lower_temperature_gradient(max_vector_index))
    allocate(upper_temperature_gradient(max_vector_index))

    allocate(averaged_dissipation_on_grid(number_of_layers, 360, 180))


    test = .FALSE.
    write_simulation_data = .TRUE.

    if (write_simulation_data) then
        ! Get the current date and time
        call date_and_time(date, time, values=values)
    
        ! Create timestamp in the format YYYYMMDD_HHMMSS
        write(folder_name, '(I4.4,I2.2,I2.2,"_",I2.2,I2.2,I2.2)') &
            values(1), values(2), values(3), values(5), values(6), values(7)
    
        ! Append folder description
        folder_name = trim(folder_name) // "_simulation_data"
    
        ! Create the directory
        call system("mkdir -p " // trim(folder_name))
    
        print *, "Folder created: ", trim(folder_name)
    
        ! Open and close required files using loops
    
        ! Temperature files
        do i = 1, 9
            write(file_name, '(A,I1,A)') trim(folder_name)//"/temperature", i, ".txt"
            open(unit=10+i, file=trim(file_name), status="replace", action="write")
            close(10+i)
        end do
    
        ! Dissipation files
        do i = 1, 9
            write(file_name, '(A,I1,A)') trim(folder_name)//"/dissipation", i, ".txt"
            open(unit=20+i, file=trim(file_name), status="replace", action="write")
            close(20+i)
        end do
    
        ! Fluidity files
        do i = 21, 29  ! Keeping the numbering consistent with the previous convention
            write(file_name, '(A,I2,A)') trim(folder_name)//"/fluidity_", i, ".txt"
            open(unit=30+i-21, file=trim(file_name), status="replace", action="write")
            close(30+i-21)
        end do

        ! Special file: total_dissipation.txt
        file_name = trim(folder_name) // "/total_dissipation.txt"
        open(unit=50, file=file_name, status="replace", action="write")
        close(50)
    
        do i = 1, num_angles
            write(file_name, '(A,"/temperature_",I3.3,"_",I3.3,".txt")') &
                trim(folder_name), angles(i,1), angles(i,2)
            open(unit=40+i, file=trim(file_name), status="replace", action="write")
            close(40+i)
        end do

        !heat_flux_at_the_bottom
        file_name = trim(folder_name) // "/lower_heat_flux.txt"
        open(unit=60, file=file_name, status="replace", action="write")
        close(60)

        file_name = trim(folder_name) // "/upper_heat_flux.txt"
        open(unit=61, file=file_name, status="replace", action="write")
        close(61)

        !heat_flux_at_the_bottom
        file_name = trim(folder_name) // "/total_lower_heat_flux.txt"
        open(unit=62, file=file_name, status="replace", action="write")
        close(62)

        file_name = trim(folder_name) // "/total_upper_heat_flux.txt"
        open(unit=63, file=file_name, status="replace", action="write")
        close(63)
    
        print *, "Files successfully created in folder: ", trim(folder_name)
    else
        print *, "File writing is disabled."
    end if
    

    ! Initialize arrays with zero values
    displacement = 0.0d0
    cauchy = 0.0d0
    cauchy_isotropic = 0.0d0
    cauchy_integral = 0.0d0
    fluidity_2 = 0.0d0
    cauchy_times_fluidity_2 = 0.d0

    bottom_force = 0.0d0
    volume_force = 0.0d0

    matrix = 0.0d0
    matrix_for_j1 = 0.0d0
    toroidial_matrix = 0.0d0
    
    Q_in_time = 0.0d0
    Total_averaged_dissipation = 0.0D0
    dissipation = 0.D0
    temperature = 0.D0
    log_temperature = 0.D0
    log_temperature_previous = 0.D0
    upper_dirichlet_temperature = 0.D0
    lower_dirichlet_temperature = 0.D0


    call set_up_initial_values_for_temperature_and_log_temperature(number_of_layers, jmax, temperature, log_temperature, log_temperature_previous, upper_dirichlet_temperature, lower_dirichlet_temperature)

    call GNDD0(jmax)

    do k=1, number_or_total_rounds

        print*, "Running total round:", k

        dissipation = 0
        Q_in_time = 0
        Total_averaged_dissipation = 0
        averaged_dissipation_on_grid = 0

        if (write_simulation_data) then
            call calculate_and_write_heat_flux(jmax, number_of_layers, k, number_of_steps_for_heat_equation, heat_equation_delta_t, radius, delta_r, thickness, k_0, temperature,  folder_name)
        end if

        if (write_simulation_data) then
            call write_temperature_data(number_of_layers, jmax, t, k, number_of_steps_for_heat_equation, heat_equation_delta_t, delta_r, thickness, temperature, folder_name, num_angles, angles)
        end if

        call calculate_fluidity_2(number_of_layers, jmax, eta_0, temperature, fluidity_2)

        if (write_simulation_data) then
            call write_fluidity_2_data(number_of_layers, jmax, t, k, number_of_steps_for_heat_equation, heat_equation_delta_t, fluidity_2, folder_name)
        end if

        ! Main simulation loop over time steps, careful to number_of_time_steps - 1
        do t=0,number_of_time_steps_for_deformaion-1

            !call CPU_TIME(end_time)

            !elapsed_time = end_time - start_time

            !print*, "Elapsed time in one cycle :", elapsed_time, "seconds"

            call CPU_TIME(start_time)

            call calculate_forces(number_of_layers, t, volume_force, bottom_force, radius, delta_r, angular_speed, delta_t, ice_density, delta_rho, excentricity)

            call update_cauchy_integral(jmax, number_of_layers, delta_t, eta, mu, cauchy_integral, cauchy, fluidity_2, cauchy_times_fluidity_2)

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
            
            if (t>99) then

                call update_grid_dissipation(number_of_layers, jmax, number_of_time_steps_for_deformaion, cauchy, fluidity_2, averaged_dissipation_on_grid)

            end if

        end do

        call dissipation_from_grid_and_total(number_of_layers, jmax, radius, delta_r, averaged_dissipation_on_grid, dissipation, Total_averaged_dissipation)

        ! Check if Total_averaged_dissipation is NaN or negative
        if (Total_averaged_dissipation /= Total_averaged_dissipation .or. Total_averaged_dissipation < 0.0) then
            print*, "Error: Total averaged dissipation is NaN or negative. Terminating program."
            stop
        end if

        !!!!!!
        !dissipation = 0

        if (write_simulation_data) then
            call write_dissipation(number_of_layers, jmax, k, number_of_steps_for_heat_equation, heat_equation_delta_t, dissipation, Total_averaged_dissipation, folder_name)
        end if

        call initial_explicit_euler_update(number_of_layers, jmax, heat_equation_delta_t, radius, delta_r, k_0, ice_density, c_p, temperature, log_temperature, dissipation, upper_dirichlet_temperature, lower_dirichlet_temperature)

        call update_log_temperature(number_of_layers, jmax, temperature, log_temperature, log_temperature_previous)

        do t=1, number_of_steps_for_heat_equation - 1
            !print*, 'Heat equation step:', t

            call update_temperature(number_of_layers, jmax, heat_equation_delta_t, radius, delta_r, k_0, ice_density, c_p, temperature, log_temperature, log_temperature_previous, dissipation, upper_dirichlet_temperature, lower_dirichlet_temperature)

            call update_log_temperature(number_of_layers, jmax, temperature, log_temperature, log_temperature_previous)

        end do

    end do

    print*, "Final Total averaged dissipation:", Total_averaged_dissipation

    open(4, file = 'Q.dat')
        ! Write the time and displacement data to file
        do i=100, number_of_time_steps_for_deformaion-1
            write(4,*) i*(delta_t)*angular_speed/(2.0*acos(-1.0d0)), Q_in_time(i)/1e9
        end do
    close(4)

    open(5, file = 'dissipation.dat')
        ! Write the time and displacement data to file
        do j=0, jmax

            do m=0,j

                i=j*(j+1)/2+m+1

                write(5,*) j, m, dissipation(25, i)

            end do

        end do
    close(5)


end program Europa_simulation

! gfortran -ffree-line-length-none main.f90 -o main && ./main