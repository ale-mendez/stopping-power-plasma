** WE FIND THE PROTON STOPPING POWER (Lindhard 1954)
	PROGRAM SLPA
      IMPLICIT REAL*8 (A-H,O-Z)
      character(len=30) input,filename
      character(len=:),allocatable:: str
      real*8,allocatable :: xarray(:)
      real*8,allocatable :: carray(:)
      real*8,allocatable :: karray(:)
      integer m1,m2,l,col
      parameter (NVMX=28)

      common /plas/ xkdau,vthau,xkfc,xdeg,xneta,enl,tev
      COMMON /CINDIV/ VV
	common/hart/xnocc,m1,m2,l
      common/norma/xnormval
      common/rbck/rmax
      common/tolbck/TOLK,TOLW,TOLR

      EXTERNAL FUN1,GABQ,ga,wavefun,fint,fnelect,fnorm

      dimension vp_au(NVMX)

c.... Default velocity grid in a.u.
      data vp_au/   0.1D0,   0.2D0,    0.3D0,    0.4D0,    0.6D0,
     |              0.8D0,   1.0D0,    1.2D0,    1.4D0,    1.6D0,
     |              1.8D0,   2.0D0,    2.4D0,    2.8D0,    3.0D0,
     |              4.0D0,   4.5D0,    5.0D0,   6.32D0,   8.94D0,
     |           10.954D0, 12.65D0, 14.142D0, 16.733D0,  20.00D0,
     |            28.28D0, 47.00D0,  63.24D0/

      PI = 3.1415926536
      tev = 0.d0
      conv_au_to_eVcm2 = 0.76196d0
      xnormval = 1.d0
      rmax = 10.0d0
      TOLK = 1.0e-2
      TOLW = 1.0e-3
      TOLR = 1.0e-3
************      READ   DATA   ***************************
      open(10,file='slpa_ll.inp',status='old')
      READ (10,*)vmin
      READ (10,*)vmax
      CLOSE(10)
*************     READ THE HARTREE COEF      **********************
      write(*,*) 'Enter the input file'
      read(*,*) input 
      str = trim(input)
	open(11,file='Data/'//str//'.inp',status='old')
	READ(11,*) m1
	READ(11,*) m2
	allocate (xarray(m1))
	allocate (carray(m1))
	allocate (karray(m2))
      READ(11,*) (xarray(col),col=1,m1)
      READ(11,*) (carray(col),col=1,m1)  
      READ(11,*) (karray(col),col=1,m2)
      READ(11,*) l
      READ(11,*) enl
      READ(11,*) xnocc
      close(11)
c      WRITE(*,*) xarray(:)
c      WRITE(*,*) carray(:)
c      WRITE(*,*) enl
***********************************************************************
C... find maximum radial value for wavefunction integration
      call find_rmax(rmax)

c...  wavefunction normalization
      xnel = fint(fnelect)
      xnormval = fint(fnorm)
      npts = 2000
      h = rmax/npts
      do i=1,npts
        r=(i-1)*h
        write(444,*) r, (wavefun(r)*r) ** 2.0 * 4.0 * pi
      enddo
      write(*,1000) xnel,xnormval
      write(*,*) xnel,xnormval
      write(*,*) "**** RENORMALIZATION ****"
      xnel = fint(fnelect)
      write(*,1000) xnel,fint(fnorm)
      do i=1,npts
        r=(i-1)*h
        write(555,*) r, (wavefun(r)*r) ** 2.0 * 4.0 * pi
      enddo
      write(*,*) 
1000  format(1x,'Ne=',f6.2,'   Norm = ',f6.2)

      pause
***********************************************************************
      OPEN(7,FILE='SLPA-LL_'//str//'_T0.dat',status='unknown')
      WRITE(7,*)'# PROTON STOPPING POWER USING SLPA METHOD'
      write(7,*)'# WITH LEVINE-LINDHARD DIELECTRIC FUNCTION'
      write(7,*)'# TARGET '//str
      write(7,*)'# T(eV): ',tev
      WRITE(7,1002)
 1002 FORMAT(1X,'# V(au)',3X,'SP(au/atom)',4X,'SP(10^-15eVcm2/atom)')
      close(7)

* COMPUTE THE STOPPING POWER OVER A GIVEN DISTRIBUTION OF VELOCITIES
      do 200 iv=1,NVMX
         v=vp_au(iv)
         VV=v
c....    add the minimum and maximum velocity given
         if ((vp_au(iv+1).gt.vmin).and.(v.lt.vmin)) v = vmin
         if ((v.gt.vmax).and.(vp_au(iv-1).lt.vmax)) v = vmax
         if ((v.lt.vmin).or.(v.gt.vmax)) go to 200         
* INTEGRATION LIMITS FOR ALL EXCITATIONS:
* THE K-INTEGRATION GOES TO INFINITY, BUT IT IS SUFFICIENT TO
* TAKE EIK2 TWICE THE VALUE FOR LINDHARD EIK2, I.E. THE MAXIMUM
* TRANSFERABLE MOMENTUM IN ABSENCE OF DAMPING
         EIK1 = 1.D-5
         EIK2 = 4.d0*(v+xkfc)
* STOPPING
         CALL GABQ(FUN1,EIK1,EIK2,SUM,TOLK,IER)
         spp = 2.d0/PI/V**2 * SUM
         write(*,*)'v=',v,'  sp',spp
         OPEN(7,FILE='SLPA-LL_'//str//'_T0.dat',status='old',
     |            access='append')
         WRITE(7,1001) v,spp,spp*conv_au_to_eVcm2
         close(7)

 200  CONTINUE
 1001 FORMAT(2X,F6.3,3X,1pE12.6,3x,e12.6,3x,e12.6)
      stop
	deallocate(xarray)
	deallocate(carray)
	deallocate(karray)
      END

*********  FUNCTIONS  AND  SUBROUTINES    **************************
      subroutine find_rmax(rmax)
      IMPLICIT REAL*8 (A-H,O-Z)
      external wavefun

      eps = 1.0e-6
      do i=1,100
        r=rmax+i*0.5
        wavemax = wavefun(r)
        if (wavemax.lt.eps) goto 100
      enddo
 100  rmax = r

      return
      end
C  **************************************************************
      FUNCTION FUN1(K)
*     ----------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 K,KK,wp(20),gamma(20),alpha(20),vf(20)
c
      COMMON /CINDIV/ V
      COMMON /CK/ KK,W
      EXTERNAL GABQ1,FW
      common/tolbck/TOLK,TOLW,TOLR
c
      KK  = K
* INTEGRATE OVER OMEGA EVERYTHING FROM ZERO TO KV
      W1  = 1.E-5
      W2  = K*V
      CALL GABQ1(FW,W1,W2,SUM1,TOLW,IER)
      FUN1 = 1.d0/K* SUM1
      RETURN
      END
C  **************************************************************
      FUNCTION FW(W)
      IMPLICIT REAL*8(A-H,K,O-Z)
	REAL*8 K,W,WW
      COMMON/CK/K,WW
      common/rbck/rmax
      common/tolbck/TOLK,TOLW,TOLR
	EXTERNAL GABQ2,FR
      PI=3.1415926536
	WW=W
*  Integrate over r from zero to 10
	R1 = 1.e-5
	R2 = rmax
      CALL GABQ2(FR,R1,R2,RES,TOLR,IER)
      RLPA = 4.d0*pi*RES
      FW=W*RLPA
      RETURN
      END
C  **************************************************************
      FUNCTION FR(R)
      IMPLICIT REAL*8(A-B,D-H,K,O-Z), COMPLEX*16(C)
	REAL*8 K,W
	COMMON/CK/K,W
      DATA C1/(1.D0,0.D0)/, C0/(0.D0,0.D0)/, CI/(0.D0,1.D0)/

      call elflev(w,k,r,cepsilon)
      if (cepsilon.eq.0) then
         feps=0.d0
      else
         FEPS=DIMAG(-C1/CEPSILON)
      endif
      FR=(r**2.d0)*FEPS

      RETURN
      END
C  **************************************************************
      subroutine elflev(w,xk,r,celf)
      IMPLICIT REAL*8(A-B,D-H,K,O-Z), COMPLEX*16(C)
      common/plas/xkdau,vthau,xkfc,xdeg,xneta,enl,tev
      DATA C1/(1.D0,0.D0)/, C0/(0.D0,0.D0)/, CI/(0.D0,1.D0)/
      external density_tot, lindhard_full, lindhard

      pi = 3.1415926536

      if(dabs(w).ge.enl) then
        wg = sqrt(w**2.d0-enl**2.d0)
        if (tev.eq.0.d0) then
            vf = (3*pi**2.d0*density_tot(r))**(1.d0/3.d0)
            call lindhard_full(wg,xk,vf,celf)
c            call lindhard(wg,xk,vf,celf)
        end if
      else
        celf=C0
      end if
      return
      end
C  **************************************************************
      subroutine lindhard_full(w,xk,vf,CEPSILON)
      IMPLICIT REAL*8(A-B,D-H,K,O-Z), COMPLEX*16(C)
      DATA C1/(1.D0,0.D0)/, C0/(0.D0,0.D0)/, CI/(0.D0,1.D0)/

      zero = 1.0d-4
      pi = 3.1415926536
      cxi2 = C1/(pi*vf)
      cu = (W+CI*ZERO)/(xk*vf)
      zz = xk/(2.d0*vf)

      CEPSILON=C1+cxi2/(ZZ**2)*(0.5D0 + 1.D0/(8.D0*ZZ)*(
     |    (C1-(ZZ+CU)**2)*CDLOG((ZZ+CU+C1)/(ZZ+CU-C1))+
     |    (C1-(ZZ-CU)**2)*CDLOG((ZZ-CU+C1)/(ZZ-CU-C1))))

      RETURN
      END
C  **************************************************************
      subroutine lindhard(w,xk,v0,celf)
      IMPLICIT REAL*8(A-B,D-H,K,O-Z), COMPLEX*16(C)
      DATA C1/(1.D0,0.D0)/, C0/(0.D0,0.D0)/, CI/(0.D0,1.D0)/
      external gfun

      PI=3.1415926536
      xi2=1.d0/(pi*v0)
      ukw=w/(xk*v0)
      zkw=xk/(2.d0*v0)

      f1=0.5d0+1.d0/(8.d0*zkw)*(gfun(zkw-ukw)+gfun(zkw+ukw))

      if((zkw+ukw).lt.1.d0) then
         f2=pi*ukw/2.d0
      else if(dabs(zkw-ukw).lt.1.d0.and.(zkw+ukw).gt.1.d0) then
         f2=pi*(1-(zkw-ukw)**2)/(8*zkw)
      else if(dabs(zkw-ukw).gt.1.d0) then 
         f2=0.d0
      endif

      celf=1.d0+xi2/(zkw**2.d0)*dcmplx(f1,f2)

      RETURN
      END
C  **************************************************************
      function gfun(x)
      IMPLICIT REAL*8(A-H,K,O-Z)
      gfun=(1.d0-x**2.d0)*dlog(dabs((x+1.d0)/(x-1.d0)))
      return
      end
C  **************************************************************
      function gdeg(x)
      implicit double precision (a-h,o-z)
      gdeg=x+(1-x**2.d0)*dlog(dabs((1+x)/(x-1)))/2
      return
      end
C  **************************************************************
C                  radial functions
C  **************************************************************
      function wavefun(r)
C     This function takes the input values saved in karray, xarray and
C     carray and constructs the Harte-Fock wave function for the
C     selected shell. For more information about this inputs see C1s.inp
      implicit double precision (a-h,o-z)
      real*8,allocatable :: xarray(:)
      real*8,allocatable :: carray(:)
      real*8,allocatable :: karray(:)
      integer m1,m2,l,n,ksize,s,i,j
      common/hart/xnocc,m1,m2,l
      common/norma/xnormval
      external xi

      j = 1
      ksize = size(karray)/2
      s = 0
      orb = 0.d0
      DO WHILE (j.le.ksize)
        k = karray(2*j)
        n = karray(2*j-1)
        DO i = 1, k
c            if (r.eq.0) write(*,*) n, l, 'C=', carray(i+s)
            orb = orb + carray(i+s)*xi(r,n,l,i+s)
        END DO
        s = s + k
        j = j + 1
      END DO
      wavefun = orb / sqrt(xnormval)
      return
      end
C  **************************************************************
      function density_tot(r)
C     This function takes the input values saved in karray, xarray and
C     carray and constructs the Harte-Fock wave function for the
C     selected shell. For more information about this inputs see C1s.inp
      implicit double precision (a-h,o-z)
      real*8,allocatable :: xarray(:)
      real*8,allocatable :: carray(:)
      real*8,allocatable :: karray(:)
      common/hart/xnocc,m1,m2,l
      external wavefun
      density_tot = xnocc*(wavefun(r)**2.d0)
      return
      end
C  **************************************************************
      function xkf(r)
      implicit double precision (a-h,o-z)
      pi=3.1415926536
      xkf = (3*density_tot(r)*pi**2.d0)**(1.d0/3.d0)
      return
      end
C  **************************************************************
	function xi(r,n,l,i)
	implicit double precision (a-h,o-z)
      real*8,allocatable :: xarray(:)
      real*8,allocatable :: carray(:)
      real*8,allocatable :: karray(:)
	integer n,i
	real*8 N10,xa,fact
	pi=3.1415926536
C	This function constructs the X functions for each electron
C	in a given shell. r is the current radius, n is the shell 
C	level, l is the orbital type and i is the index of the 
C	coefficient for each electron.
      fact = 1.0
      do j = 2, 2*n
        fact = fact*j
      end do
      xa = xarray(i)
      Yl = sqrt((2.d0*l+1)/4.d0/pi)
      N10 = 1.0/sqrt(fact)*(2.d0*xa)**(n+0.5)
      xi = N10*(r**(n-1))*exp(-xa*r)*Yl
c      if (r.eq.0) write(*,*) n,l,fact,'xa=',xa,' N=',N10
      return
      end
C  **************************************************************
	real function fact(n)
	integer i
	fact = 1
	do i = 2, n
	   fact = fact*i
	end do
	return
	end
C  **************************************************************
      function fint(func)
      implicit double precision (a-h,o-z)
      external gabq,func
      R1=0.005
      R2=10.d0
      tol=1.d-4
      CALL GABQ(func,R1,R2,sume,tol,ie)
      fint=sume
      return
      end
C  **************************************************************
      function fnelect(r)
      implicit double precision (a-h,o-z)
      common/hart/xnocc,m1,m2,l
      external wavefun
      pi=3.1415926536
      fnelect=4*pi*xnocc*(r*wavefun(r))**2.d0
      return
      end
C  **************************************************************
      function fnorm(r)
      implicit double precision (a-h,o-z)
      common/hart/xnocc,m1,m2,l
      external density_tot
      pi=3.1415926536
      fnorm=4*pi*(r*wavefun(r))**2.d0
      return
      end
C  **************************************************************
C                      SUBROUTINE GABQ
C  **************************************************************
      SUBROUTINE GABQ(FCT,XL,XU,SUM,TOL,IER)
C
C      THIS INTEGRATION SUBROUTINE APPLIES THE GAUSS METHOD WITH
C   AN ADAPTIVE-BIPARTITION SCHEME.
C      FCT IS THE (EXTERNAL) FUNCTION BEING INTEGRATED OVER THE
C   INTERVAL (XL,XU). SUM IS THE RESULTANT VALUE OF THE INTEGRAL.
C      TOL IS THE TOLERANCE, I.E. MAXIMUM RELATIVE ERROR REQUIRED
C   ON THE COMPUTED VALUE (SUM). TOL SHOULD NOT EXCEED 1.0D-13.
C      IER IS AN ERROR CONTROL PARAMETER; ITS OUTPUT VALUE IS
C   IER=0 IF THE INTEGRATION ALGORITHM HAS BEEN ABLE TO GET THE
C   REQUIRED ACCURACY AND IER=1 OTHERWISE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(10),W(10),S(128),SN(128),L(128),LN(128)
C  ****  PRINTED OUTPUT OF PARTIAL RESULTS: SET IWR=1.
      DATA IWR/0/
C  ****  COEFFICIENTS FOR GAUSS 20-POINT INTEGRATION.
      DATA NP,NP2,NP4/10,20,40/
C  ****  ABSCISAS.
      DATA X/7.6526521133497334D-02,2.2778585114164508D-01,
     1       3.7370608871541956D-01,5.1086700195082710D-01,
     2       6.3605368072651503D-01,7.4633190646015079D-01,
     3       8.3911697182221882D-01,9.1223442825132591D-01,
     4       9.6397192727791379D-01,9.9312859918509492D-01/
C  ****  WEIGHTS.
      DATA W/1.5275338713072585D-01,1.4917298647260375D-01,
     1       1.4209610931838205D-01,1.3168863844917663D-01,
     2       1.1819453196151842D-01,1.0193011981724044D-01,
     3       8.3276741576704749D-02,6.2672048334109064D-02,
     4       4.0601429800386941D-02,1.7614007139152118D-02/
C  ****  CORRECTED TOLERANCE.
      CTOL=DMAX1(TOL,1.0D-13)
      PTOL=0.01D0*CTOL
      H=XU-XL
C
      IF(IWR.EQ.1) THEN
      WRITE(6,10)
   10 FORMAT(///5X,'GAUSS ADAPTIVE-BIPARTITION QUADRATURE')
      WRITE(6,11) XL,XU,TOL
   11 FORMAT(/5X,'XL = ',1PD15.8,', XU = ',D15.8,', TOL =',
     1 D8.1)
      ENDIF
      IER=0
C  ****  GAUSS INTEGRATION FROM XL TO XU.
      A=0.5D0*(XU-XL)
      B=0.5D0*(XL+XU)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 1 I1=2,NP
      C=A*X(I1)
    1 D=D+W(I1)*(FCT(B+C)+FCT(B-C))
      SUM=D*A
C  ****  ADAPTIVE BIPARTITION SCHEME.
      ICALL=NP2
      LH=1
      S(1)=SUM
      L(1)=1
    2 HO=H
      H=0.5D0*H
      ASUM=SUM
      LHN=0
      DO 5 I=1,LH
      K=L(I)
      SI=S(I)
      XA=XL+(K-1)*HO
      XB=XA+H
      XC=XA+HO
      A=0.5D0*(XB-XA)
      B=0.5D0*(XB+XA)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 3 I2=2,NP
      C=A*X(I2)
    3 D=D+W(I2)*(FCT(B+C)+FCT(B-C))
      S1=D*A
      A=0.5D0*(XC-XB)
      B=0.5D0*(XC+XB)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 4 I3=2,NP
      C=A*X(I3)
    4 D=D+W(I3)*(FCT(B+C)+FCT(B-C))
      S2=D*A
      ICALL=ICALL+NP4
      S12=S1+S2
      SUM=SUM+S12-SI
      IF(DABS(S12-SI).LT.DMAX1(PTOL*DABS(S12),1.0D-35)) GO TO 5
      LHN=LHN+2
      IF(LHN.GT.128.OR.ICALL.GT.9999) GO TO 8
      SN(LHN)=S2
      LN(LHN)=K+K
      SN(LHN-1)=S1
      LN(LHN-1)=LN(LHN)-1
   5  CONTINUE
      ERR=DABS(SUM-ASUM)/DMAX1(DABS(SUM),1.0D-35)
      IF(IWR.EQ.1) WRITE(6,12) ICALL,SUM,ERR,LHN
   12 FORMAT(5X,'N =',I5,', SUM =',1PD19.12,', ERR =',D8.1,
     1 ', LH =',I3)
      IF(ERR.GT.CTOL.AND.LHN.GT.0) GO TO 6
      IF(IWR.EQ.1) WRITE(6,13)
   13 FORMAT(5X,'END OF GAUSS-BIPARTITION PROCEDURE'///)
      RETURN
    6 LH=LHN
      DO 7 I=1,LH
      S(I)=SN(I)
    7 L(I)=LN(I)
      GO TO 2
C  ****  WARNING (LOW ACCURACY) MESSAGE.
    8 WRITE(6,14)
   14 FORMAT(/5X,'*** LOW ACCURACY IN SUBROUTINE GABQ.')
      WRITE(6,11) XL,XU,TOL
      WRITE(6,15) SUM,ERR
   15 FORMAT(5X,'SUM =',1PD19.12,', ERR =',D8.1//)
      IER=1
      RETURN
      END
C  **************************************************************
C                      SUBROUTINE GABQ1
C  **************************************************************
      SUBROUTINE GABQ1(FCT,XL,XU,SUM,TOL,IER)
C
C      THIS INTEGRATION SUBROUTINE APPLIES THE GAUSS METHOD WITH
C   AN ADAPTIVE-BIPARTITION SCHEME.
C      FCT IS THE (EXTERNAL) FUNCTION BEING INTEGRATED OVER THE
C   INTERVAL (XL,XU). SUM IS THE RESULTANT VALUE OF THE INTEGRAL.
C      TOL IS THE TOLERANCE, I.E. MAXIMUM RELATIVE ERROR REQUIRED
C   ON THE COMPUTED VALUE (SUM). TOL SHOULD NOT EXCEED 1.0D-13.
C      IER IS AN ERROR CONTROL PARAMETER; ITS OUTPUT VALUE IS
C   IER=0 IF THE INTEGRATION ALGORITHM HAS BEEN ABLE TO GET THE
C   REQUIRED ACCURACY AND IER=1 OTHERWISE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(10),W(10),S(128),SN(128),L(128),LN(128)
C  ****  PRINTED OUTPUT OF PARTIAL RESULTS: SET IWR=1.
      DATA IWR/0/
C  ****  COEFFICIENTS FOR GAUSS 20-POINT INTEGRATION.
      DATA NP,NP2,NP4/10,20,40/
C  ****  ABSCISAS.
      DATA X/7.6526521133497334D-02,2.2778585114164508D-01,
     1       3.7370608871541956D-01,5.1086700195082710D-01,
     2       6.3605368072651503D-01,7.4633190646015079D-01,
     3       8.3911697182221882D-01,9.1223442825132591D-01,
     4       9.6397192727791379D-01,9.9312859918509492D-01/
C  ****  WEIGHTS.
      DATA W/1.5275338713072585D-01,1.4917298647260375D-01,
     1       1.4209610931838205D-01,1.3168863844917663D-01,
     2       1.1819453196151842D-01,1.0193011981724044D-01,
     3       8.3276741576704749D-02,6.2672048334109064D-02,
     4       4.0601429800386941D-02,1.7614007139152118D-02/
C  ****  CORRECTED TOLERANCE.
      CTOL=DMAX1(TOL,1.0D-13)
      PTOL=0.01D0*CTOL
      H=XU-XL
C
      IF(IWR.EQ.1) THEN
      WRITE(6,10)
   10 FORMAT(///5X,'GAUSS ADAPTIVE-BIPARTITION QUADRATURE')
      WRITE(6,11) XL,XU,TOL
   11 FORMAT(/5X,'XL = ',1PD15.8,', XU = ',D15.8,', TOL =',
     1 D8.1)
      ENDIF
      IER=0
C  ****  GAUSS INTEGRATION FROM XL TO XU.
      A=0.5D0*(XU-XL)
      B=0.5D0*(XL+XU)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 1 I1=2,NP
      C=A*X(I1)
    1 D=D+W(I1)*(FCT(B+C)+FCT(B-C))
      SUM=D*A
C  ****  ADAPTIVE BIPARTITION SCHEME.
      ICALL=NP2
      LH=1
      S(1)=SUM
      L(1)=1
    2 HO=H
      H=0.5D0*H
      ASUM=SUM
      LHN=0
      DO 5 I=1,LH
      K=L(I)
      SI=S(I)
      XA=XL+(K-1)*HO
      XB=XA+H
      XC=XA+HO
      A=0.5D0*(XB-XA)
      B=0.5D0*(XB+XA)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 3 I2=2,NP
      C=A*X(I2)
    3 D=D+W(I2)*(FCT(B+C)+FCT(B-C))
      S1=D*A
      A=0.5D0*(XC-XB)
      B=0.5D0*(XC+XB)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 4 I3=2,NP
      C=A*X(I3)
    4 D=D+W(I3)*(FCT(B+C)+FCT(B-C))
      S2=D*A
      ICALL=ICALL+NP4
      S12=S1+S2
      SUM=SUM+S12-SI
      IF(DABS(S12-SI).LT.DMAX1(PTOL*DABS(S12),1.0D-35)) GO TO 5
      LHN=LHN+2
      IF(LHN.GT.128.OR.ICALL.GT.9999) GO TO 8
      SN(LHN)=S2
      LN(LHN)=K+K
      SN(LHN-1)=S1
      LN(LHN-1)=LN(LHN)-1
    5 CONTINUE
      ERR=DABS(SUM-ASUM)/DMAX1(DABS(SUM),1.0D-35)
      IF(IWR.EQ.1) WRITE(6,12) ICALL,SUM,ERR,LHN
   12 FORMAT(5X,'N =',I5,', SUM =',1PD19.12,', ERR =',D8.1,
     1 ', LH =',I3)
      IF(ERR.GT.CTOL.AND.LHN.GT.0) GO TO 6
      IF(IWR.EQ.1) WRITE(6,13)
   13 FORMAT(5X,'END OF GAUSS-BIPARTITION PROCEDURE'///)
      RETURN
    6 LH=LHN
      DO 7 I=1,LH
      S(I)=SN(I)
    7 L(I)=LN(I)
      GO TO 2
C  ****  WARNING (LOW ACCURACY) MESSAGE.
    8 WRITE(6,14)
   14 FORMAT(/5X,'*** LOW ACCURACY IN SUBROUTINE GABQ1.')
      WRITE(6,11) XL,XU,TOL
      WRITE(6,15) SUM,ERR
   15 FORMAT(5X,'SUM =',1PD19.12,', ERR =',D8.1//)
      IER=1
      RETURN
      END
C  **************************************************************
C                      SUBROUTINE GABQ2
C  **************************************************************
      SUBROUTINE GABQ2(FCT,XL,XU,SUM,TOL,IER)
C
C      THIS INTEGRATION SUBROUTINE APPLIES THE GAUSS METHOD WITH
C   AN ADAPTIVE-BIPARTITION SCHEME.
C      FCT IS THE (EXTERNAL) FUNCTION BEING INTEGRATED OVER THE
C   INTERVAL (XL,XU). SUM IS THE RESULTANT VALUE OF THE INTEGRAL.
C      TOL IS THE TOLERANCE, I.E. MAXIMUM RELATIVE ERROR REQUIRED
C   ON THE COMPUTED VALUE (SUM). TOL SHOULD NOT EXCEED 1.0D-13.
C      IER IS AN ERROR CONTROL PARAMETER; ITS OUTPUT VALUE IS
C   IER=0 IF THE INTEGRATION ALGORITHM HAS BEEN ABLE TO GET THE
C   REQUIRED ACCURACY AND IER=1 OTHERWISE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(10),W(10),S(128),SN(128),L(128),LN(128)
C  ****  PRINTED OUTPUT OF PARTIAL RESULTS: SET IWR=1.
      DATA IWR/0/
C  ****  COEFFICIENTS FOR GAUSS 20-POINT INTEGRATION.
      DATA NP,NP2,NP4/10,20,40/
C  ****  ABSCISAS.
      DATA X/7.6526521133497334D-02,2.2778585114164508D-01,
     1       3.7370608871541956D-01,5.1086700195082710D-01,
     2       6.3605368072651503D-01,7.4633190646015079D-01,
     3       8.3911697182221882D-01,9.1223442825132591D-01,
     4       9.6397192727791379D-01,9.9312859918509492D-01/
C  ****  WEIGHTS.
      DATA W/1.5275338713072585D-01,1.4917298647260375D-01,
     1       1.4209610931838205D-01,1.3168863844917663D-01,
     2       1.1819453196151842D-01,1.0193011981724044D-01,
     3       8.3276741576704749D-02,6.2672048334109064D-02,
     4       4.0601429800386941D-02,1.7614007139152118D-02/
C  ****  CORRECTED TOLERANCE.
      CTOL=DMAX1(TOL,1.0D-13)
      PTOL=0.01D0*CTOL
      H=XU-XL
C
      IF(IWR.EQ.1) THEN
      WRITE(6,10)
   10 FORMAT(///5X,'GAUSS ADAPTIVE-BIPARTITION QUADRATURE')
      WRITE(6,11) XL,XU,TOL
   11 FORMAT(/5X,'XL = ',1PD15.8,', XU = ',D15.8,', TOL =',
     1 D8.1)
      ENDIF
      IER=0
C  ****  GAUSS INTEGRATION FROM XL TO XU.
      A=0.5D0*(XU-XL)
      B=0.5D0*(XL+XU)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 1 I1=2,NP
      C=A*X(I1)
    1 D=D+W(I1)*(FCT(B+C)+FCT(B-C))
      SUM=D*A
C  ****  ADAPTIVE BIPARTITION SCHEME.
      ICALL=NP2
      LH=1
      S(1)=SUM
      L(1)=1
    2 HO=H
      H=0.5D0*H
      ASUM=SUM
      LHN=0
      DO 5 I=1,LH
      K=L(I)
      SI=S(I)
      XA=XL+(K-1)*HO
      XB=XA+H
      XC=XA+HO
      A=0.5D0*(XB-XA)
      B=0.5D0*(XB+XA)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 3 I2=2,NP
      C=A*X(I2)
    3 D=D+W(I2)*(FCT(B+C)+FCT(B-C))
      S1=D*A
      A=0.5D0*(XC-XB)
      B=0.5D0*(XC+XB)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 4 I3=2,NP
      C=A*X(I3)
    4 D=D+W(I3)*(FCT(B+C)+FCT(B-C))
      S2=D*A
      ICALL=ICALL+NP4
      S12=S1+S2
      SUM=SUM+S12-SI
      IF(DABS(S12-SI).LT.DMAX1(PTOL*DABS(S12),1.0D-35)) GO TO 5
      LHN=LHN+2
      IF(LHN.GT.128.OR.ICALL.GT.9999) GO TO 8
      SN(LHN)=S2
      LN(LHN)=K+K
      SN(LHN-1)=S1
      LN(LHN-1)=LN(LHN)-1
    5 CONTINUE
      ERR=DABS(SUM-ASUM)/DMAX1(DABS(SUM),1.0D-35)
      IF(IWR.EQ.1) WRITE(6,12) ICALL,SUM,ERR,LHN
   12 FORMAT(5X,'N =',I5,', SUM =',1PD19.12,', ERR =',D8.1,
     1 ', LH =',I3)
      IF(ERR.GT.CTOL.AND.LHN.GT.0) GO TO 6
      IF(IWR.EQ.1) WRITE(6,13)
   13 FORMAT(5X,'END OF GAUSS-BIPARTITION PROCEDURE'///)
      RETURN
    6 LH=LHN
      DO 7 I=1,LH
      S(I)=SN(I)
    7 L(I)=LN(I)
      GO TO 2
C  ****  WARNING (LOW ACCURACY) MESSAGE.
    8 WRITE(6,14)
   14 FORMAT(/5X,'*** LOW ACCURACY IN SUBROUTINE GABQ2.')
      WRITE(6,11) XL,XU,TOL
      WRITE(6,15) SUM,ERR
   15 FORMAT(5X,'SUM =',1PD19.12,', ERR =',D8.1//)
      IER=1
      RETURN
      END
