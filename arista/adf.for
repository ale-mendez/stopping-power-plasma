c
c  This program computes the PROTON STOPPING POWER for all plasma 
c  degeneracy using Arista's dielectric function (PRA84)
c
      program ABtheory
      IMPLICIT REAL*8 (A-H,O-Z)
      character (24) filename
      character(len=:),allocatable:: str
      parameter (NVMX=28)

      COMMON/CINDIV/VV
      common/tolbck/TOLK,TOLW,TOLX
      common/plas_atom/xkf
      common/plas_temp/tev,xdeg,xneta

      EXTERNAL FUN1,GABQ
c      external ga

      dimension vp_au(NVMX)

c.... Default velocity grid in a.u.
      data vp_au/   0.1D0,   0.2D0,    0.3D0,    0.4D0,    0.6D0,
     |              0.8D0,   1.0D0,    1.2D0,    1.4D0,    1.6D0,
     |              1.8D0,   2.0D0,    2.4D0,    2.8D0,    3.0D0,
     |              4.0D0,   4.5D0,    5.0D0,   6.32D0,   8.94D0,
     |           10.954D0, 12.65D0, 14.142D0, 16.733D0,  20.00D0,
     |            28.28D0, 47.00D0,  63.24D0/

      pi = 3.1415926536
      a0 = 0.529177210903E-8
      xkB = 8.617333262E-5    ! eV/K 
      conv_cm3_to_au = 0.08916
      conv_au_to_eVcm2 = 0.76196d0
      conv_au_to_eV = 27.2114
      conv_eV_to_au = 1.d0/conv_au_to_eV
      conv_eV_to_K = 11604.525
      TOLK = 1.0e-3
      TOLW = 1.0e-3
      TOLX = 1.0e-3
      theta = 0.d0
      xdeg = 0.d0
      xneta = 0.d0

c... Read input data 
      open(10,file='adf.inp',status='old')
      READ(10,*) vmin
      READ(10,*) vmax
      read(10,*) tev
      read(10,*) iopt
      ! choose if input is ratio atomic_density/mass or electron density (cm^-3) 
      if (iopt.eq.1) then
         read(10,*) atomic_density_ratio
      else if(iopt.eq.2) then
         read(10,*) electron_density !cm^-3
      else
         print*, 'iopt must be'
         print*, '  = 1 : ratio =  density [g/cm^3] / atomic_mass [g]'
         print*, '  = 2 : electron density [cm^-3]'
         stop
      endif
      read(10,*) xnumber_electrons
      close(10)

c... Compute parameters and conversions
      if (iopt.eq.1) then
         atomic_density_au = conv_cm3_to_au * atomic_density_ratio
         atomic_density = atomic_density_au / a0 ** 3
         electron_density_au = xnumber_electrons * atomic_density_au
         electron_density = electron_density_au / a0 ** 3
      else if (iopt.eq.2) then
         electron_density_au = electron_density * a0 ** 3
         atomic_density_au = electron_density_au / xnumber_electrons
         atomic_density = atomic_density_au / a0 ** 3
      endif
      xkf = (3*pi**2*electron_density_au)**(1.d0/3.d0)
      wp = sqrt(4.0*pi*electron_density_au)
      ef_au = 0.5d0*xkf**2.d0

      tK = tev*conv_eV_to_K
      xkBT_au = xkB * tK * conv_eV_to_au
      if (tev.eq.0) goto 100

      theta = xkBT_au/ef_au
      xdeg = 1.d0/theta
      xneta = finv12_mathcad(2.d0/3.d0*xdeg**1.5d0)

100   continue
      WRITE(*,1000) atomic_density,atomic_density_au,
     |          xnumber_electrons,electron_density,electron_density_au,
     |          tev,xdeg,xneta

c... Write output file
      write(*,*) 'Enter the filename'
      read(*,*) filename

      str=trim(filename)//'_ADF.dat'
      OPEN(7,FILE=str,status='unknown')
      WRITE(7,1002) trim(filename)
      write(7,1000) atomic_density,atomic_density_au,
     |          xnumber_electrons,electron_density,electron_density_au,
     |          tev
      WRITE(7,1003)
      close(7)

1002  format('# PROTON STOPPING POWER ON PLASMA'/,
     |       '# USING THE ARISTA DIELECTRIC FUNCTION'/,
     |       '# Target:',1x,a)
1000  format('# Atomic density:',1x,1pe12.4,1x,'cm^-3,',1x,
     |                           0pf8.4,1x,'a.u.'/,
     |       '# Number of electrons: ',f4.2/,
     |       '# Electron density:',1x,1pe12.4,1x,'cm^-3,',1x,
     |                              0pf8.4,1x,'a.u.'/,
     |       '# Temperature: ',f8.4,1x,'eV'/,
     |       '# Degeneracy: ',f8.4/,
     |       '# Eta: ',f8.4/)
1003  FORMAT(/'# V(au)',3X,'SP(au)',9X,'SPCS(10^-15eVcm2/atom)')

c... Compute stopping power over a given distribution of velocities
      do 200 iv=1,NVMX
        v = vp_au(iv)
        VV = v

c.... add the minimum and maximum velocity given
        if ((vp_au(iv+1).gt.vmin).and.(v.lt.vmin)) v = vmin
        if ((v.gt.vmax).and.(vp_au(iv-1).lt.vmax)) v = vmax
        if ((v.lt.vmin).or.(v.gt.vmax)) go to 200

c.... define integration limits for k variable
        EIK1 = 1.D-5
        if (v.lt.1) epsmax = 50.d0
        if (v.ge.1) epsmax = 5.d0
        EIK2 = epsmax*v

c.... compute integral
        CALL GABQ(FUN1,EIK1,EIK2,SUM,TOLK,IER)
        spp = 2./PI/V**2 * SUM
        spp_conv = spp*conv_au_to_eVcm2/atomic_density_au

c.... print results in terminal and output file
        write(*,*) 'v=',v,'  sp=',spp_conv
        open(7,file=str,status='old',access='append')
        write(7,1004) v,spp,spp_conv
        close(7)
200   continue

1004  format(2X,F6.3,3X,1pE12.6,3x,e12.6,3x,e12.6)
      stop
      deallocate(str)
      end program
c
c-----------------------------------------------------------------------
c
c Energy loss integrand over k
c
c-----------------------------------------------------------------------
c
      FUNCTION FUN1(K)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 K,KK,wp(20),gamma(20),alpha(20),vf(20)

      COMMON/CINDIV/V
      common/CK/KK
      common/tolbck/TOLK,TOLW,TOLX
      EXTERNAL GABQ1,fun2

      KK  = K
      ! INTEGRATE OVER OMEGA EVERYTHING FROM ZERO TO KV
      W1  = 1.E-5
      W2  = K*V
      CALL GABQ1(fun2,W1,W2,SUM1,TOLW,IER)
      FUN1 = 1./K* SUM1
      RETURN
      END
c
c-----------------------------------------------------------------------
c
c Energy loss integrand over omega
c
c-----------------------------------------------------------------------
c
      FUNCTION Fun2(W)
      IMPLICIT REAL*8 (A-Z)
      COMMON/CK/K
      common/plas_temp/tev,xdeg,xneta
      external elfplasdeg,elfplas
      external lindhard_full,lindhard_gamma0

      if (tev.eq.0) then
         call lindhard_full(w,k,inv_elf)
         !call elfplasdeg(w,k,ximelf)
      else
         call elfplas(w,k,inv_elf)
      end if
      Fun2=inv_elf*W

      RETURN
      END
c
c-----------------------------------------------------------------------
c
c Full dielectric constant function - Lindhard 1954
c
c-----------------------------------------------------------------------
c
      subroutine lindhard_full(w,xk,FEPS)
      IMPLICIT REAL*8(A-B,D-H,K,O-Z), COMPLEX*16(C)
      common/plas_atom/xkf
      DATA C1/(1.D0,0.D0)/, C0/(0.D0,0.D0)/, CI/(0.D0,1.D0)/

      zero = 1.0d-4
      pi = 3.1415926536
      cxi2 = C1/(pi*xkf)
      cu = (W+CI*ZERO)/(xk*xkf)
      zz = xk/(2.d0*xkf)

      CEPSILON=C1+cxi2/(ZZ**2)*(0.5D0 + 1.D0/(8.D0*ZZ)*(
     |    (C1-(ZZ+CU)**2)*CDLOG((ZZ+CU+C1)/(ZZ+CU-C1))+
     |    (C1-(ZZ-CU)**2)*CDLOG((ZZ-CU+C1)/(ZZ-CU-C1))))
	
      FEPS=DIMAG(-C1/CEPSILON)
      RETURN
      END
c
c-----------------------------------------------------------------------
c
c Approximated dielectric constant function - Lindhard 1954
c
c-----------------------------------------------------------------------
c
      subroutine lindhard_gamma0(w,xk,FEPS)
      IMPLICIT REAL*8(A-B,D-H,K,O-Z), COMPLEX*16(C)
      common/plas_atom/xkf
      external g_lindhard
      DATA C1/(1.D0,0.D0)/, C0/(0.D0,0.D0)/, CI/(0.D0,1.D0)/

      pi = 3.1415926536
      xi2 = 1.d0/(pi*xkf)
      uu = w/(xk*xkf)
      zz = xk/(2.d0*xkf)

      feps1 = 0.5+1.0/(8.d0*zz)*(g_lindhard(zz-uu)+g_lindhard(zz+uu))
      if ((zz+uu).le.1.d0) then
         feps2 = pi*uu/2
      elseif (dabs(zz-uu).lt.1.d0.and.(zz+uu).gt.1.d0) then
         feps2 = pi/(8.d0*zz)*(1.d0-(zz-uu)**2.d0)
      elseif (dabs(zz-uu).gt.1.d0) then
         feps2 = 0.d0
      endif

      ceL = C1+xi2/(zz**2)*dcmplx(feps1,feps2)
      FEPS=DIMAG(-C1/ceL)
      
      RETURN
      END
c
c-----------------------------------------------------------------------
c
c Function G in Lindhard's dielectric approximation - Lindhard 1964
c
c-----------------------------------------------------------------------
c
      function g_lindhard(x)
      implicit double precision (a-h,o-z)
      g_lindhard=(1-x**2.d0)*dlog(dabs((x+1)/(x-1)))
      return
      end
c
c-----------------------------------------------------------------------
c
c 
c
c-----------------------------------------------------------------------
c
      subroutine elfplasdeg(w,xk,elf)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 df
      common/plas_atom/xkf

      PI=3.1415926536
	ukw=w/(xk*xkf)
	zkw=xk/(2.d0*xkf)
      wgranre=1+(gdeg(ukw+zkw)-gdeg(ukw-zkw))/(4.d0*pi*xkf*zkw**3.d0)

	if((zkw+ukw).lt.1.d0) f2=pi*ukw/2.d0
	if(dabs(zkw-ukw).lt.1.d0.and.(zkw+ukw).gt.1.d0) 
	1  f2=pi*(1.d0-(zkw-ukw)**2)/(8.d0*zkw)
	if(dabs(zkw-ukw).gt.1.d0) f2=0.d0

      wgranim=f2/(pi*xkf*zkw**2.d0)
	df=dcmplx(wgranre,wgranim)
      elf=dimag(-1.0d0/df)

      return
      end
c
c-----------------------------------------------------------------------
c
c 
c
c-----------------------------------------------------------------------
c
      function gdeg(x)
      implicit double precision (a-h,o-z)
      gdeg=x+0.5d0*(1-x**2.d0)*dlog(dabs((1+x)/(1-x)))
      return
      end
c
c-----------------------------------------------------------------------
c
c Arista dielectric function - Arista 1984
c
c-----------------------------------------------------------------------
c
      subroutine elfplas(w,xk,elf)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 df
      common/plas_atom/xkf
      common/plas_temp/tev,xdeg,xneta
      external ga

      PI=3.1415926536
	ukw=w/(xk*xkf)
	zkw=xk/(2.d0*xkf)
      eps1=1.d0+(ga(ukw+zkw)-ga(ukw-zkw))/(4.d0*pi*xkf*zkw**3.d0)

	xdenime=1.d0+dexp(xneta-xdeg*dabs((ukw-zkw)**2.d0))
	xdenima=1.d0+dexp(xneta-xdeg*dabs((ukw+zkw)**2.d0))
      eps2=dlog(xdenime/xdenima)/(8.d0*xdeg*xkf*zkw**3.d0)

	df=dcmplx(eps1,eps2)
      elf=dimag(-1.0d0/df)

      return 
      end
c
c-----------------------------------------------------------------------
c
c Function g in Arista dielectric function - Arista 1984
c
c-----------------------------------------------------------------------
c
      function ga(x)
      implicit double precision (a-h,o-z)
	common/gai/xx
      common/plas_temp/tev,xdeg,xneta
      common/tolbck/TOLK,TOLW,TOLX
      external gaint,gabq2

	xx=x
      eiy1=1.d-5
      eiy2=dsqrt((709.d0+xneta)/xdeg)
      if (x.eq.0) print*,xneta,xdeg,'xmax=',eiy2
      call gabq2(gaint,eiy1,eiy2,xsum,tolx,ier)
	ga=xsum

      return
      end
c
c-----------------------------------------------------------------------
c
c Integrand defined in g(x) function in Arista dielectric function - Arista 1984
c
c-----------------------------------------------------------------------
c
      function gaint(y)
      implicit double precision (a-h,o-z)
      common/plas_temp/tev,xdeg,xneta
	common/gai/x

	xln=dlog(dabs((x+y)/(x-y)))
	gaint=y*xln/(dexp(xdeg*y**2.d0-xneta)+1.d0)

      return
      end
c
c-----------------------------------------------------------------------
c
c Inverse approximation of Fermi integral of order 1/2
c
c-----------------------------------------------------------------------
c
      function finv12_mathcad(y)
      implicit double precision (a-h,o-z)
      dimension a(3), b(3), c(3), d(3)
      data a/44.593646d0, 11.288764d0,1.d0/
      data b/39.519346d0,-5.7517464d0,0.26594291d0/
      data c/34.873722d0,-26.922515d0,1.d0/
      data d/26.612832d0,-20.452930d0,11.808945d0/

      if (y.lt.4.d0) then
        R1 = (a(1) + a(2)*y + a(3)*y*y)/(b(1) + b(2)*y + b(3)*y*y)
        x = dlog(y*R1)
      else if (y.ge.4.d0) then
        yp = y ** (-1.d0/(1.d0+0.5d0))
        R2 = (c(1) + c(2)*yp + c(3)*yp*yp)/(d(1) + d(2)*yp + d(3)*yp*yp)
        x = y ** (1.d0/(1.d0+0.5d0)) * R2
      endif
      finv12_mathcad = x

      return
      end
c
c-----------------------------------------------------------------------
c
c Interpolation of inverse Fermi function results given in Arista 1984
c
c-----------------------------------------------------------------------
c
      FUNCTION FINV12_cspline(Y)
      implicit double precision (a-h,o-z)
      DIMENSION E0(10),V(10)
      DATA A0/0.497563805159625769d0/
      DATA A1/-0.929879930744899986d0/
      DATA A2/0.34452875966595822d0/
      DATA A3/0.125896539214666517d0/
      DATA A4/0.131353080600043063d-1/
      DATA A5/0.510979233863366652d-3/
      DATA B0/0.440954441108243294d0/
      DATA B1/-1.0d0/
      DATA B2/0.671962074844631728d0/
      DATA B3/-0.868779148174670455d-1/
      DATA B4/0.505599207758294305d-2/
      DATA B5/-0.12821991486196731d-3/
      DATA C0/-0.314130764307281443d0/
      DATA C1/-0.652493120734329158d0/
      DATA C2/1.0d0/
      DATA C3/0.87131946995868494d0/
      DATA C4/0.117446015507808365d0/
      DATA C5/0.192964073354827288d-2/
      DATA D0/0.923099668858667349d-1/
      DATA D1/0.771768060514091258d0/
      DATA D2/0.906449500868349234d0/
      DATA D3/0.215506116242806155d0/
      DATA D4/0.775738101060550066d-2/
      DATA D5/0.102785488357280288d-4/
      DATA YMIN/0.129301317461299139d0/
      DATA YMAX/0.416279240188605876d0/
      DATA E0(1)/0.593713769772101558d-1/
      DATA E0(2)/0.564617640522927885d-1/
      DATA E0(3)/0.865275738864030861d-2/
      DATA E0(4)/0.442063532429317682d-3/
      DATA E0(5)/0.657336866246965319d-4/
      DATA E0(6)/0.599193349974194576d-5/
      DATA E0(7)/0.293776657109598689d-6/
      DATA E0(8)/0.980093595837500779d-8/
      DATA E0(9)/0.847860996944584164d-8/
      DATA E0(10)/0.764722847007293529d-9/
      DATA G0/0.303029119211352746d0/
      DATA G1/-0.794138569610706479d0/
      DATA G2/0.246820656511986589d-2/
      DATA G3/-0.157979823631583854d-1/
      DATA G4/-0.356952563200875895d0/
      DATA G5/1.0d0/
      DATA H0/0.231254499113070739d0/
      DATA H1/-0.606041154122076666d0/
      DATA H2/0.188359278609885606d-2/
      DATA H3/-0.120560117511733673d-1/
      DATA H4/-0.161642516107504751d0/
      DATA H5/0.472938502365927207d0/
      IF(Y .LT. 1.39637528063739999d0) THEN
      Y1 = Y
      Y2 = Y * Y
      Y3 = Y2 * Y
      Y4 = Y3 * Y
      Y5 = Y4 * Y
      Y6 = Y5 * Y
      ZNUM=A0*Y1+A1*Y2+A2*Y3+A3*Y4+A4*Y5+A5*Y6
      ZDEN=B0+B1*Y1+B2*Y2+B3*Y3+B4*Y4+B5*Y5
      FINV12_cspline = LOG(ZNUM/ZDEN)
      RETURN
      ELSE IF(Y.LT.5.77072652741599997d0) THEN
      Y1=Y
      Y2 = Y1*Y
      Y3 = Y2*Y
      Y4 = Y3*Y
      Y5 = Y4*Y
      ZNUM=C0+C1*Y1+C2*Y2+C3*Y3+C4*Y4+C5*Y5
      ZDEN=D0+D1*Y1+D2*Y2+D3*Y3+D4*Y4+D5*Y5
      FINV12_cspline = ZNUM / ZDEN
      RETURN
      ELSE IF(Y .LT. 0.598127953660604876d2) THEN
      !TRANSFORMATION DE LA VARIABLE
      Z=1.d0/SQRT(Y)
      ZZ=((Z-YMIN)-(YMAX-Z))/(YMAX-YMIN)
      V(1)=1.d0
      V(2)=ZZ
      DO 1515 KK=3,10
      V(KK)=2*ZZ*V(KK-1)-V(KK-2)
1515  continue
      P=0.D0
      DO 1516 KK=10,1,-1
      P=P+E0(KK)*V(KK)
1516  continue
      FINV12_cspline=(1.d0/P**2)**(1./3.)
      RETURN
      ELSE
      Y1=Y**(2./3.)
      Y2=(Y)**(1./3.)
      Y4=1.0d0/Y2
      Y5=1.0d0/Y1
      Y6=1.0d0/Y
      Y7=Y6*(Y6)**(1./3.)
      Y8=Y6*Y6**(2./3.)
      ZNUM=G0*Y1+G1*Y2+G2+G3*Y4+G4*Y5+G5*Y6
      ZDEN=H0+H1*Y4+H2*Y5+H3*Y6+H4*Y7+H5*Y8
      FINV12_cspline = ZNUM / ZDEN

      RETURN
      END IF
      END
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
     1D8.1)
      ENDIF
      IER=0
C  ****  GAUSS INTEGRATION FROM XL TO XU.
      A=0.5D0*(XU-XL)
      B=0.5D0*(XL+XU)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 1 I1=2,NP
      C=A*X(I1)
      D=D+W(I1)*(FCT(B+C)+FCT(B-C))
    1 continue
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
      D=D+W(I2)*(FCT(B+C)+FCT(B-C))
    3 continue 
      S1=D*A
      A=0.5D0*(XC-XB)
      B=0.5D0*(XC+XB)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 4 I3=2,NP
      C=A*X(I3)
      D=D+W(I3)*(FCT(B+C)+FCT(B-C))
    4 continue 
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
     1', LH =',I3)
      IF(ERR.GT.CTOL.AND.LHN.GT.0) GO TO 6
      IF(IWR.EQ.1) WRITE(6,13)
   13 FORMAT(5X,'END OF GAUSS-BIPARTITION PROCEDURE'///)
      RETURN
    6 LH=LHN
      DO 7 I=1,LH
      S(I)=SN(I)
      L(I)=LN(I)
    7 continue
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
C  ****  PRINTED OUTPUT OF PARTIAL RESULTS: SET IWR=1.
      DATA IWR/0/
c      DIMENSION X(10),W(10),S(128),SN(128),L(128),LN(128)
c C  ****  COEFFICIENTS FOR GAUSS 20-POINT INTEGRATION.
c       DATA NP,NP2,NP4/10,20,40/
c C  ****  ABSCISAS.
c       DATA X/7.6526521133497334D-02,2.2778585114164508D-01,
c      1       3.7370608871541956D-01,5.1086700195082710D-01,
c      2       6.3605368072651503D-01,7.4633190646015079D-01,
c      3       8.3911697182221882D-01,9.1223442825132591D-01,
c      4       9.6397192727791379D-01,9.9312859918509492D-01/
c C  ****  WEIGHTS.
c       DATA W/1.5275338713072585D-01,1.4917298647260375D-01,
c      1       1.4209610931838205D-01,1.3168863844917663D-01,
c      2       1.1819453196151842D-01,1.0193011981724044D-01,
c      3       8.3276741576704749D-02,6.2672048334109064D-02,
c      4       4.0601429800386941D-02,1.7614007139152118D-02/

c C  ****  COEFFICIENTS FOR GAUSS 40-POINT INTEGRATION.
c       DIMENSION X(20),W(20),S(128),SN(128),L(128),LN(128)
c       DATA NP,NP2,NP4/20,40,80/
c C  ****  ABSCISAS.
c       DATA X/0.0387724175060508,0.1160840706752552,
c      1       0.1926975807013711,0.2681521850072537,
c      2       0.3419940908257585,0.4137792043716050,
c      3       0.4830758016861787,0.5494671250951282,
c      4       0.6125538896679802,0.6719566846141796,
c      5       0.7273182551899271,0.7783056514265194,
c      6       0.8246122308333117,0.8659595032122595,
c      7       0.9020988069688743,0.9328128082786765,
c      8       0.9579168192137917,0.9772599499837743,
c      9       0.9907262386994570,0.9982377097105593/
c C  ****  WEIGHTS.
c       DATA W/0.0775059479784248,0.0770398181642480,
c      1       0.0761103619006262,0.0747231690579683,
c      2       0.0728865823958041,0.0706116473912868,
c      3       0.0679120458152339,0.0648040134566010,
c      4       0.0613062424929289,0.0574397690993916,
c      5       0.0532278469839368,0.0486958076350722,
c      6       0.0438709081856733,0.0387821679744720,
c      7       0.0334601952825478,0.0279370069800234,
c      8       0.0222458491941670,0.0164210583819079,
c      9       0.0104982845311528,0.0045212770985332/

C  ****  COEFFICIENTS FOR GAUSS 64-POINT INTEGRATION.
      DIMENSION X(32),W(32),S(128),SN(128),L(128),LN(128)
      DATA NP,NP2,NP4/32,64,128/
C  ****  ABSCISAS.
      DATA X/0.0243502926634244,0.072993121787799,
     1       0.121462819296121,0.169644420423993,
     2       0.217423643740007,0.264687162208767,
     3       0.311322871990211,0.357220158337668,
     4       0.402270157963992,0.446366017253464,
     5       0.489403145707053,0.531279464019895,
     6       0.571895646202634,0.611155355172393,
     7       0.648965471254657,0.685236313054233,
     8       0.719881850171611,0.752819907260532,
     9       0.783972358943341,0.813265315122798,
     1       0.84062929625258,0.865999398154093,
     2       0.889315445995114,0.910522137078503,
     3       0.92956917213194,0.946411374858403,
     4       0.961008799652054,0.973326827789911,
     5       0.983336253884626,0.991013371476744,
     6       0.996340116771955,0.999305041735772/
C  ****  WEIGHTS.
      DATA W/0.0486909570091397,0.0485754674415034,
     1       0.048344762234803,0.0479993885964583,
     2       0.0475401657148303,0.04696818281621,
     3       0.0462847965813144,0.0454916279274181,
     4       0.0445905581637566,0.0435837245293235,
     5       0.0424735151236536,0.0412625632426235,
     6       0.0399537411327203,0.0385501531786156,
     7       0.03705512854024,0.0354722132568824,
     8       0.0338051618371416,0.0320579283548516,
     9       0.0302346570724025,0.0283396726142595,
     1       0.0263774697150547,0.0243527025687109,
     2       0.0222701738083833,0.0201348231535302,
     3       0.0179517157756973,0.0157260304760247,
     4       0.0134630478967186,0.0111681394601311,
     5       0.0088467598263639,0.0065044579689784,
     6       0.0041470332605625,0.0017832807216964/

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
     1D8.1)
      ENDIF
      IER=0
C  ****  GAUSS INTEGRATION FROM XL TO XU.
      A=0.5D0*(XU-XL)
      B=0.5D0*(XL+XU)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 1 I1=2,NP
      C=A*X(I1)
      D=D+W(I1)*(FCT(B+C)+FCT(B-C))
    1 continue
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
          D=D+W(I2)*(FCT(B+C)+FCT(B-C))
    3   continue  
        S1=D*A
        A=0.5D0*(XC-XB)
        B=0.5D0*(XC+XB)
        C=A*X(1)
        D=W(1)*(FCT(B+C)+FCT(B-C))
        DO 4 I3=2,NP
          C=A*X(I3)
          D=D+W(I3)*(FCT(B+C)+FCT(B-C))
    4   continue 
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
     1', LH =',I3)
      IF(ERR.GT.CTOL.AND.LHN.GT.0) GO TO 6
      IF(IWR.EQ.1) WRITE(6,13)
   13 FORMAT(5X,'END OF GAUSS-BIPARTITION PROCEDURE'///)
      RETURN
    6 LH=LHN
      DO 7 I=1,LH
      S(I)=SN(I)
      L(I)=LN(I)
    7 continue
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
C  ****  PRINTED OUTPUT OF PARTIAL RESULTS: SET IWR=1.
      DATA IWR/0/
c      DIMENSION X(10),W(10),S(128),SN(128),L(128),LN(128)
c C  ****  COEFFICIENTS FOR GAUSS 20-POINT INTEGRATION.
c       DATA NP,NP2,NP4/10,20,40/
c C  ****  ABSCISAS.
c       DATA X/7.6526521133497334D-02,2.2778585114164508D-01,
c      1       3.7370608871541956D-01,5.1086700195082710D-01,
c      2       6.3605368072651503D-01,7.4633190646015079D-01,
c      3       8.3911697182221882D-01,9.1223442825132591D-01,
c      4       9.6397192727791379D-01,9.9312859918509492D-01/
c C  ****  WEIGHTS.
c       DATA W/1.5275338713072585D-01,1.4917298647260375D-01,
c      1       1.4209610931838205D-01,1.3168863844917663D-01,
c      2       1.1819453196151842D-01,1.0193011981724044D-01,
c      3       8.3276741576704749D-02,6.2672048334109064D-02,
c      4       4.0601429800386941D-02,1.7614007139152118D-02/

c C  ****  COEFFICIENTS FOR GAUSS 40-POINT INTEGRATION.
c       DIMENSION X(20),W(20),S(128),SN(128),L(128),LN(128)
c       DATA NP,NP2,NP4/20,40,80/
c C  ****  ABSCISAS.
c       DATA X/0.0387724175060508,0.1160840706752552,
c      1       0.1926975807013711,0.2681521850072537,
c      2       0.3419940908257585,0.4137792043716050,
c      3       0.4830758016861787,0.5494671250951282,
c      4       0.6125538896679802,0.6719566846141796,
c      5       0.7273182551899271,0.7783056514265194,
c      6       0.8246122308333117,0.8659595032122595,
c      7       0.9020988069688743,0.9328128082786765,
c      8       0.9579168192137917,0.9772599499837743,
c      9       0.9907262386994570,0.9982377097105593/
c C  ****  WEIGHTS.
c       DATA W/0.0775059479784248,0.0770398181642480,
c      1       0.0761103619006262,0.0747231690579683,
c      2       0.0728865823958041,0.0706116473912868,
c      3       0.0679120458152339,0.0648040134566010,
c      4       0.0613062424929289,0.0574397690993916,
c      5       0.0532278469839368,0.0486958076350722,
c      6       0.0438709081856733,0.0387821679744720,
c      7       0.0334601952825478,0.0279370069800234,
c      8       0.0222458491941670,0.0164210583819079,
c      9       0.0104982845311528,0.0045212770985332/

C  ****  COEFFICIENTS FOR GAUSS 64-POINT INTEGRATION.
      DIMENSION X(32),W(32),S(128),SN(128),L(128),LN(128)
      DATA NP,NP2,NP4/32,64,128/
C  ****  ABSCISAS.
      DATA X/0.0243502926634244,0.072993121787799,
     1       0.121462819296121,0.169644420423993,
     2       0.217423643740007,0.264687162208767,
     3       0.311322871990211,0.357220158337668,
     4       0.402270157963992,0.446366017253464,
     5       0.489403145707053,0.531279464019895,
     6       0.571895646202634,0.611155355172393,
     7       0.648965471254657,0.685236313054233,
     8       0.719881850171611,0.752819907260532,
     9       0.783972358943341,0.813265315122798,
     1       0.84062929625258,0.865999398154093,
     2       0.889315445995114,0.910522137078503,
     3       0.92956917213194,0.946411374858403,
     4       0.961008799652054,0.973326827789911,
     5       0.983336253884626,0.991013371476744,
     6       0.996340116771955,0.999305041735772/
C  ****  WEIGHTS.
      DATA W/0.0486909570091397,0.0485754674415034,
     1       0.048344762234803,0.0479993885964583,
     2       0.0475401657148303,0.04696818281621,
     3       0.0462847965813144,0.0454916279274181,
     4       0.0445905581637566,0.0435837245293235,
     5       0.0424735151236536,0.0412625632426235,
     6       0.0399537411327203,0.0385501531786156,
     7       0.03705512854024,0.0354722132568824,
     8       0.0338051618371416,0.0320579283548516,
     9       0.0302346570724025,0.0283396726142595,
     1       0.0263774697150547,0.0243527025687109,
     2       0.0222701738083833,0.0201348231535302,
     3       0.0179517157756973,0.0157260304760247,
     4       0.0134630478967186,0.0111681394601311,
     5       0.0088467598263639,0.0065044579689784,
     6       0.0041470332605625,0.0017832807216964/

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
     1D8.1)
      ENDIF
      IER=0
C  ****  GAUSS INTEGRATION FROM XL TO XU.
      A=0.5D0*(XU-XL)
      B=0.5D0*(XL+XU)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 1 I1=2,NP
      C=A*X(I1)
      D=D+W(I1)*(FCT(B+C)+FCT(B-C))
    1 continue  
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
      D=D+W(I2)*(FCT(B+C)+FCT(B-C))
    3 continue  
      S1=D*A
      A=0.5D0*(XC-XB)
      B=0.5D0*(XC+XB)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO 4 I3=2,NP
      C=A*X(I3)
      D=D+W(I3)*(FCT(B+C)+FCT(B-C))
    4 continue  
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
     1', LH =',I3)
      IF(ERR.GT.CTOL.AND.LHN.GT.0) GO TO 6
      IF(IWR.EQ.1) WRITE(6,13)
   13 FORMAT(5X,'END OF GAUSS-BIPARTITION PROCEDURE'///)
      RETURN
    6 LH=LHN
      DO 7 I=1,LH
      S(I)=SN(I)
      L(I)=LN(I)
    7 continue
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
