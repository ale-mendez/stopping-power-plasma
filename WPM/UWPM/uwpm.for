c
c   This program computes stopping cross sections using the Unified Wave 
c   Packet Model (UWPM) proposed by Archubi & Arista (2020)
c
      program UWPM
      IMPLICIT REAL*8 (A-H,O-Z) 
      character (24) systime
      character(len=30) filename
      character(len=:),allocatable:: str
      parameter (NVMX=28)

      common/CINDIV/VV
      common/tolbck/TOLK,TOLW,TOLX
      common/plas/enl,chi2,qbar

      EXTERNAL fun1,GABQ

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
      conv_au_to_A = 51.4194586
      conv_au_to_eV = 27.2114
      conv_eV_to_au = 1.d0/conv_au_to_eV
      conv_eV_to_K = 11604.525
      TOLK = 1.0e-3
      TOLW = 1.0e-3
      TOLX = 1.0e-3

      enl = 0.d0
c**********      READ DATA      ***********************
      open(10,file='uwpm.inp',status='old')
      read(10,*) vmin
      read(10,*) vmax
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
c      read(10,*) number_electrons
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
      rs_au = (3.0/(4.0*pi*electron_density_au))**(1.0/3.0)

      tK = tev*conv_eV_to_K
      xkBT_au = xkB * tK * conv_eV_to_au
      xlambda = 0.4d0
      xkBTeff = (xkBT_au**2.d0 + xlambda * ef_au**2.d0)**(1.d0/2.d0)
      theta_eff = xkBTeff / ef_au
      qbar = xkf * theta_eff**(1.d0/2.d0)
      chi2 = 4.d0/(3.d0*pi**(3.d0/2.d0)*xkf*theta_eff**2.d0)

      write(*,*)'Compute FEG stopping with UWPM for plasmas'
      write(*,*)'Temperature (eV)',tev,xkBT_au,xlambda*ef_au**2.d0
      write(*,*)'kf=',xkf,'ef=',ef_au
      write(*,*)'theta=',theta_eff,xkBTeff
      write(*,*)'chi^2=',chi2
      write(*,*)'qbar=',qbar

      WRITE(*,1000) atomic_density,atomic_density_au,
     |          xnumber_electrons,electron_density,electron_density_au,
     |          tev

c**********      OUTPUT FILE     ***************************
      WRITE(*,*)'Enter the output file name'
      READ(*,*) filename

      str = trim(filename)//'_UWPM.dat'
      open(7,file=str,status='unknown')
      WRITE(7,1002) trim(filename)
      write(7,1000) atomic_density,atomic_density_au,
     |          xnumber_electrons,electron_density,electron_density_au,
     |          tev
      WRITE(7,1003)
      CLOSE(7)

1002  format('# PROTON STOPPING POWER ON PLASMA'/,
     |       '# USING THE ARISTA DIELECTRIC FUNCTION'/,
     |       '# Target:',1x,a)
1000  format('# Atomic density:',1x,1pe12.4,1x,'cm^-3,',1x,
     |                           0pf8.4,1x,'a.u.'/,
     |       '# Number of electrons: ',f6.2/,
     |       '# Electron density:',1x,1pe12.4,1x,'cm^-3,',1x,
     |                              0pf8.4,1x,'a.u.'/,
     |       '# Temperature: ',f8.4,1x,'eV')
1001  format('# Degeneracy: ',f8.4/,
     |       '# Eta: ',f8.4/)
1003  FORMAT(/'# V(au)',4X,'SP(au)',9X,'SPCS(10^-15eVcm2/atom)')

***********************************************************************
* COMPUTE THE STOPPING POWER OVER A GIVEN DISTRIBUTION OF VELOCITIES
      do 200 iv=1,NVMX
        v = vp_au(iv)
        VV = v
c.... add the minimum and maximum velocity given
        if ((vp_au(iv+1).gt.vmin).and.(v.lt.vmin)) v = vmin
        if ((v.gt.vmax).and.(vp_au(iv-1).lt.vmax)) v = vmax
        if ((v.lt.vmin).or.(v.gt.vmax)) go to 200

c.... define integration limits for k variable
        EIK1 = 1.D-5
c        EIK2 = 4.d0*(v+xkf)
        if (v.lt.1) epsmax = 50.d0
        if (v.ge.1) epsmax = 5.d0
        EIK2 = epsmax*v

c.... compute integral
        call gabq(fun1,EIK1,EIK2,sum,TOLK,IER)
        spp = 2/(pi*v**2) * sum
        spp_conv = spp*conv_au_to_eVcm2/atomic_density_au

c.... print results in terminal and output file
        write(*,*) 'v=',v,'sp=', spp_conv
        open(7,file=str,status='old',access='append')
        write(7,1004) v, spp, spp_conv
        close(7)
200   continue

1004  format(2X,F6.3,3X,1pE12.6,3x,e12.6,3x,e12.6)
      stop
      deallocate(str)
      end program

c********  FUNCTIONS AND SUBROUTINES **********************
      FUNCTION FUN1(K)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 K,KK
c
        COMMON /CINDIV/ V
        COMMON /CK/ KK
        EXTERNAL GABQ,fun2
        KK  = K
        W1  = 1.E-5
        W2  = K*V
        TOL = 1.E-6
        CALL GABQ(fun2,W1,W2,SUM1,TOL,IER)
        FUN1 = 1.d0/K* SUM1
        RETURN
       END
*********************************************************
      FUNCTION fun2(W)
       IMPLICIT REAL*8 (A-H,O-Z)
       REAL*8 K,W
       COMMON/CK/K
       EXTERNAL elfuwpm

       CALL elfuwpm(w,k,ximelf)
       fun2 = ximelf*W

       RETURN
       END
*********************************************************
      SUBROUTINE elfuwpm(w,xk,elf)
       IMPLICIT REAL*8 (A-H,O-Z)
       COMPLEX*16 df
       common/plas/enl,chi2,xkf

       pi = 3.1415926536
       if (enl.eq.0) then
         wg = w
       elseif (enl.gt.0) then
         print*, 'not implemented. only FEG'
         stop
       endif
       ukw = wg / (xk*xkf)
       zkw = xk / (2.d0*xkf)
       e1 = 1 + chi2*(G(ukw+zkw) - G(ukw-zkw)) / (8.d0*zkw**3.d0)
       e2 = pi*chi2*(exp(-(ukw-zkw)**2.d0) - exp(-(ukw+zkw)**2.d0))
       e2 = e2 / (8.d0*zkw**3.d0)
       df = dcmplx(e1,e2)
       elf = dimag(-1.0d0/df)

       return
       END
*********************************************************
      function G(x)
       implicit double precision (a-h,o-z)
       external GABQ, gaint
       common/tolbck/TOLK,TOLW,TOLX
       common/gai/xx

       pi = 3.1415926536

       xx = x
       eit1 = 1.d-5
       eit2 = 1.0d0
       call GABQ(gaint,eit1,eit2,sum2,tolx,ier)
       G = 2.d0*sqrt(pi)*x*sum2

       return
       end
*********************************************************
      function gaint(t)
       implicit double precision (a-h,o-z)
       common/gai/x
       xind = (t**2.d0 - 1)*x**2.d0
       gaint=dexp(xind)
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
