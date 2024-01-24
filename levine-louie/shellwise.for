*
*  This program computes stopping cross sections using the Levine and 
*  Louie dielectric function (1982)
*
	program shellwise
      implicit real*8 (a-h,o-z)
      complex*8 ceps
      character(len=30) filename
      character(len=:),allocatable:: str
      parameter (NVMX=28)

      common/cindiv/vv
      common/tolbck/TOLK,TOLW
      common/plas_atom/xkf,enl

      external fun1,gabq1

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
      conv_cm3_to_au = 0.08916
      conv_au_to_eVcm2 = 0.76196d0
      TOLK = 1.0e-3
      TOLW = 1.0e-3

************      READ   DATA   ***************************
      open(10,file='shellwise.inp',status='old')
      READ(10,*) vmin
      READ(10,*) vmax
      read(10,*) atomic_density_ratio
      read(10,*) number_electrons
      read(10,*) enl
      read(10,*) filename
      close(10)

      atomic_density_au = conv_cm3_to_au * atomic_density_ratio
      atomic_density = atomic_density_au / a0 ** 3
      electron_density_au = number_electrons * atomic_density_au
      electron_density = electron_density_au / a0 ** 3
      xkf = (3*pi**2*electron_density_au)**(1.d0/3.d0)
      wp=sqrt(4.0*pi*electron_density_au)
      
      WRITE(*,1000) atomic_density,atomic_density_au,
     |           number_electrons,electron_density,electron_density_au

************    OUTPUT FILE   ***************************
      str = trim(filename)//'_shellwise'//'.dat'
      OPEN(7,FILE=str,status='unknown')
      WRITE(7,1001) trim(filename)
      WRITE(7,1000) atomic_density,atomic_density_au,
     |           number_electrons,electron_density,electron_density_au
      WRITE(7,1002)
      close(7)
1001  format('# PROTON STOPPING POWER USING SHELLWISE METHOD'/,
     |       '# WITH THE LEVINE-LOUIE DIELECTRIC FUNCTION'/,
     |       '# target:',1x,a)

1000  format('# Atomic density:',1x,1pe12.4,1x,'cm^-3,',1x,
     |                           0pf8.4,1x,'a.u.'/,
     |       '# Number of electrons: ',i4/,
     |       '# Electron density:'1x,1pe10.4,1x,'cm^-3,',
     |                              0pf8.4,1x,'a.u.')
1002  format(/'# V(au)',3X,'SP(au/atom)',4X,'SP(10^-15eVcm2/atom)')

***********************************************************************
* COMPUTE THE STOPPING POWER OVER A GIVEN DISTRIBUTION OF VELOCITIES
      do 200 iv = 1,NVMX
         v = vp_au(iv)
         VV = v

c.... add the minimum and maximum velocity given
         if ((vp_au(iv+1).gt.vmin).and.(v.lt.vmin)) v = vmin
         if ((v.gt.vmax).and.(vp_au(iv-1).lt.vmax)) v = vmax
         if ((v.lt.vmin).or.(v.gt.vmax)) go to 200         

c.... define integration limits for k variable
         EIK1 = 1.D-5
         if (v.lt.1) epsmax = 50.d0
         if (v.ge.1) epsmax = 10.d0
         EIK2 = epsmax*v

c.... compute integral
         call gabq(fun1,EIK1,EIK2,sum,TOLK,IER)
         spp = 2.d0/(pi*v**2) * sum
         spp_conv = spp*conv_au_to_eVcm2/atomic_density_au

c.... print results in terminal and output file
         write(*,*)'v=',v,'  sp=',spp_conv         
         open(7,file=str,status='old',access='append')
         write(7,1003) v,spp,spp_conv
         close(7)
 200  continue

 1003 format(2X,F6.3,3X,1pE12.6,3x,e12.6,3x,e12.6)
      stop
      deallocate(str)
      end program

C  **************************************************************
      FUNCTION FUN1(K)
       IMPLICIT REAL*8 (A-H,O-Z)
       REAL*8 K,KK,wp(20),gamma(20),alpha(20),vf(20)

       COMMON /CINDIV/ V
       COMMON /CK/ KK,W
       EXTERNAL GABQ1,FW
       common/tolbck/TOLK,TOLW

       KK  = K
       W1  = 1.E-5
       W2  = K*V
       CALL GABQ1(FW,W1,W2,SUM1,TOLW,IER)
       FUN1 = 1.d0/K* SUM1
       RETURN
       END
c
C  **************************************************************
      FUNCTION FW(W)
       IMPLICIT REAL*8(A-H,K,O-Z)
       complex*8 ceps,C1
       REAL*8 K,W,WW
       COMMON/CK/K,WW
       common/plas_atom/xkf,enl
       EXTERNAL levinelouie
       DATA C1/(1.D0,0.D0)/

       WW=W
       call levinelouie(w,k,ceps)
       feps = imag(-C1/ceps)
       if (feps.ne.feps) feps = 0.d0
       FW = W*feps

       RETURN
       END
c
C  **************************************************************
      subroutine levinelouie(w,xk,ceps)
       implicit real*8(A-B,D-H,K,O-Z), complex*8(C)
       common/plas_atom/xkf,enl
       external fM, lindhard_full
       DATA C1/(1.D0,0.D0)/, C0/(0.D0,0.D0)/, CI/(0.D0,1.D0)/

       wg = sqrt(w**2.d0-enl**2.d0)
       if(dabs(w).ge.enl) then
         call lindhard_full(wg,xk,cepsilon)
         ceps = cepsilon
       else
         ceps = fM(xk,w)
       endif

       return
       end
c
      subroutine lindhard_full(w,xk,cepsilon)
       implicit real*8(A-B,D-H,K,O-Z), complex*8(C)
       common/plas_atom/xkf,enl
       DATA C1/(1.D0,0.D0)/, C0/(0.D0,0.D0)/, CI/(0.D0,1.D0)/
       external cgfun

       zero = 1.0d-4
       pi = 3.1415926536
       cxi2 = C1/(pi*xkf)
       cu = (W+CI*ZERO)/(xk*xkf)
       cz = C1*xk/(2.d0*xkf)

       cterm = 1.D0/(8.D0*cz)*(cgfun(cz+cu) + cgfun(cz-cu))
       cepsilon = C1 + cxi2/(cz**2)*(0.5D0 + cterm)

       RETURN
       END
c
      function cgfun(cx)
       IMPLICIT COMPLEX*8(c)
       DATA C1/(1.D0,0.D0)/, C0/(0.D0,0.D0)/, CI/(0.D0,1.D0)/
       cgfun=(c1-cx**2.d0)*clog((cx+c1)/(cx-c1))
       return
       end
c
      function fM(xk,w)
        implicit real*8 (a-h,o-z)
        common/plas_atom/xkf,enl

        pi=3.1415926536
        ef=xkf**2/2
        xlambda=enl/ef

        Qk=xk/xkf
        Dw=sqrt(xlambda**2-(w/ef)**2)
        xi2=1.d0/(pi*xkf)

        fatan=atan((2*Qk+Qk**2)/Dw)+atan((2*Qk-Qk**2)/Dw)
        famp=Dw**2/(8*Qk**5)+1.d0/(2*Qk**3)-1.d0/(8*Qk)
        flog=dlog((Dw**2+(2*Qk+Qk**2)**2)/(Dw**2+(2*Qk-Qk**2)**2))
        fM=1.d0+2.d0*xi2*(1/Qk**2-Dw/(2.d0*Qk**3)*fatan+famp*flog)
        return
        end
c
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
   14 FORMAT(/5X,'*** LOW ACCURACY IN SUBROUTINE GABQ1.')
      WRITE(6,11) XL,XU,TOL
      WRITE(6,15) SUM,ERR
   15 FORMAT(5X,'SUM =',1PD19.12,', ERR =',D8.1//)
      IER=1
      RETURN
      END