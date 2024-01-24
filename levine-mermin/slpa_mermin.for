** WE FIND THE PROTON STOPPING POWER AND THE STRAGGLING FUNCTION
** para toda degeneración del plasma (AristaPRA84)
	PROGRAM SLPA
C     use PORTLIB
	use commonarrays
      IMPLICIT REAL*8 (A-H,O-Z)
c   integer*4 i1,i2,i3
      character (24) systime
      common /plas/ xkdau,vthau,xkfc,xnu,xdegc,xnetac,enl
      COMMON /CINDIV/ VV
      EXTERNAL FUN1,GABQ,FFUN1,ga
      integer m1,m2,l,col
      character(len=30) input,filename
      character(len=:),allocatable:: str
	common/hart/m1,m2,l


      PI = 3.1415926536
************      READ   DATA   ***************************
      open(10,file='slpa_mermin.inp',status='old')
* VELOCITY PARAMETERS: standard 1-32
      READ (10,*)VI1
      READ (10,*)VI2
      READ (10,*)VI3
      msp=nint((vi2-vi1)/vi3)+2
      read(10,*)cne
      read(10,*)tev
      read(10,*)xnu
      CLOSE(10)
c
      cld=743.0d0*dsqrt(tev)/dsqrt(cne)
      xkdau=0.529d-8/cld
      vthau=4.19d5*dsqrt(tev)/2.18769d6
      xneau=cne*(0.529d-8)**3.d0
      wpau=dsqrt(4.d0*pi*xneau)
      xkfc=(3.d0*pi**2.d0*xneau)**(1.d0/3.d0)
      xdegc=(xkfc**2.d0)/(2.d0*vthau**2.d0)
      xnetac=finv12((2.d0*xdegc**1.5d0)/3.d0)
      xcoupc=xkfc/(xkfc**2.d0+2.d0*vthau**2.d0)/pi
      write(*,*)'Temperature (eV)',tev
      write(*,*)'Electronic density (cm-3)',cne
      write(*,*)'Degeneracy parameter',xdegc
      write(*,*)'Reduced temperature',1.d0/xdegc
      write(*,*)'Coupling parameter',xcoupc

c   pause
*************     READ THE HARTREE COEF      **********************
      write(*,*) 'Enter the input file'
      read(*,*) input 
      str = trim(input)
	open(11,file='./Data/'//str//'.inp',status='old')
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
      close(11)
      WRITE(*,*) xarray(:)
      WRITE(*,*) carray(:)
      WRITE(*,*) enl
      WRITE(*,*) 'Enter the output file name'
      READ(*,*) filename
      filename = trim(filename)
      pause
***********************************************************************
      OPEN(7,FILE=filename//'_mermin.dat')
      WRITE(7,*)'PROTON STOPPING POWER AND STRAGGLING'
      write(7,*)'FOR PLASMA DIELECTRIC FUNCTION'
      write(7,*)'TARGET PlasDeut '
      write(7,*)'tev',tev
      write(7,*)'cne',cne
      WRITE(7,*)'Kmax = 2*V'
      WRITE(7,1002)
 1002 FORMAT('V(au)',5X,'SP (au)',5X,'SP Bethe')
      close(7)
* THE FOLLOWING IS A PIECEWISE EQUIDISTANT DISTRIBUTION OF VELOCITIES
      DO 200 VR=VI1,VI2,VI3
      v=vr
      VV=Vr
      write(*,*)'v=',v
      systime = CTIME (TIME())
      write(*,*) 'Date and time: ',systime
* INTEGRATION LIMITS FOR ALL EXCITATIONS:
* THE K-INTEGRATION GOES TO INFINITY, BUT IT IS SUFFICIENT TO
* TAKE EIK2 TWICE THE VALUE FOR LINDHARD EIK2, I.E. THE MAXIMUM
* TRANSFERABLE MOMENTUM IN ABSENCE OF DAMPING
      EIK1 = 1.D-5
c      EIK2 = 4.*(V+VFer)
      eik2=4.d0*(v+xkfc)
c      eik2=2.5d2*v
      TOL = 1.E-2
* STOPPING
      CALL GABQ(FUN1,EIK1,EIK2,SUM,TOL,IER)
c      SP(IVR) = 2./PI/V**2 * SUM
      SPp = 2./PI/V**2 * SUM
	write(*,*)'sp',spp
* STRAGGLING
C      CALL GABQ(FFUN1,EIK1,EIK2,SSUM,TOL,IER)
c      STRAG(IVR) = 2./PI/V**2 * SSUM
C      STRAGg = 2./PI/V**2 * SSUM
c      WRITE(7,1001)VR11(IVR),Sp(IVR),strag(ivr)
      spbethe=(wpau/v)**2*dlog(2.d0*v**2.d0/wpau)
      OPEN(7,FILE=filename//'_mermin.dat',status='old')
      do i=1,msp+8
        read(7,*,end=75)
      end do
75    continue
      backspace 7
      WRITE(7,1001)VR,Spp,spbethe
      close(7)
 200  CONTINUE
 1001 FORMAT(2X,F6.3,3X,E12.6,3x,e12.6,3x,e12.6)
      stop
	deallocate(xarray)
	deallocate(carray)
	deallocate(karray)

      END

*********  FUNCTIONS  AND  SUBROUTINES    **************************
* stopping power
      FUNCTION FUN1(K)
*     ----------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 K,KK,wp(20),gamma(20),alpha(20),vf(20)
c
      COMMON /CINDIV/ V
      COMMON /CK/ KK,W
      EXTERNAL GABQ1,fun2
c
      KK  = K
* INTEGRATE OVER OMEGA EVERYTHING FROM ZERO TO KV
      W1  = 1.E-5
      W2  = K*V
      TOL = 1.E-2
      CALL GABQ1(fun2,W1,W2,SUM1,TOL,IER)
      FUN1 = 1.d0/K* SUM1
      RETURN
      END
********************************************************************
* straggling
      FUNCTION FFUN1(K)
*     ----------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 K,KK,wp(20),gamma(20),alpha(20),vf(20)
c
      COMMON /CINDIV/ V
      COMMON/CK/KK,W
      EXTERNAL GABQ1,FFun2
c
      KK  = K
* INTEGRATE OVER OMEGA EVERYTHING FROM ZERO TO KV
      W1  = 1.E-5
      W2  = K*V
      TOL = 1.E-2
      CALL GABQ1(FFun2,W1,W2,SSUM1,TOL,IER)
      FFUN1 = 1./K* SSUM1
      RETURN
      END
C  **************************************************************
      FUNCTION Fun2(W)
      IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 K,W,WW
      COMMON/CK/K,WW
	EXTERNAL GABQ2,FUN3

c
cc  open(1,file='elfarista.dat')
cc  do k=0.01d0,1.5d0,0.01d0
c   k=1.d-3*xkdau
cc  do w=0.01d0,1.5d0,0.01d0
cc      CALL elfplaspardeg(w,k,ximelf)
cc  write(*,*)k,w,ximelf
cc  write(1,*)k,w,ximelf
cc  enddo
cc  enddo
cc  close(1)
cc  pause
c
	WW=W
*  Integrate over r from zero to 10
	R1 = 1.e-5
	R2 = 10.e0
	TOL = 1.E-2
      CALL GABQ2(FUN3,R1,R2,SUM2,TOL,IER)
c      if (k.gt.7.and.k.lt.10.and.w.gt.10.and.w.lt.15) then
c            write(*,*) 'k',k, 'w', w, 'sum2', sum2
c      end if
      Fun2=sum2*W
      RETURN
      END
C  **************************************************************
      FUNCTION FUN3(R)
      IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 K,W
	COMMON/CK/K,W
      PI=3.1415926536

      CALL elflev(w,k,r,ximelf)
      fun3=ximelf*4.d0*pi*(r**2.d0)
      RETURN
      END
C  **************************************************************
      FUNCTION FFun2(W)
      IMPLICIT REAL*8 (A-Z)
      COMMON/CK/K,WW
      common/plas/xkdau,vthau,xkfc,xnu,xdegc,xnetac,enl
      if(w/(vthau**2.d0).gt.400) then
        ngran=1.0d0
        goto 1
      endif
      invn=dexp(w/(vthau**2.d0))-1.0d0
      ngran=2.0d0/invn+1.0d0
1     CALL elflev(w,k,r,ximelf)
      FFun2=ximelf*W*W*ngran
      RETURN
      END
C  **************************************************************
      subroutine elflev(w,xk,r,elf)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/plas/xkdau,vthau,xkfc,xnu,xdegc,xnetac,enl

      if(dabs(w).ge.enl) then
        wg = sqrt(w**2.d0-enl**2.d0)
        call elfpmer(wg,xk,r,elf)
      else
        elf=0
      end if
      return
      end
C  ****************************************************************
      subroutine elfpmer(w,xk,r,ximelf)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 elfw0,elf1,elf,cimg
	common/plas/xkdau,vthau,xkfc,xnu,xdegc,xnetac,enl
      cimg=dcmplx(0.d0,1.0d0)
	call elfplaslin(0.d0+cimg*0.d0,xk,r,elfw0)
	call elfplaslin(w+cimg*xnu,xk,r,elf1)
	elf=1.d0+(w+cimg*xnu)*(elf1-1.d0)*(elfw0-1.d0)/(w*(elfw0-1.d0)+
	1    cimg*xnu*(elf1-1.d0))																			   
      ximelf=DIMAG(-1.0D0/(elf))
      RETURN
      END
C  **************************************************************
      subroutine elfplaslin(w,xk,r,df)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 w,ukw,xm1,xm2,wgranim,wgranre,cimg
      COMPLEX*16 df
	common/plas/xkdau,vthau,xkfc,xnu,xdegc,xnetac,enl
      PI=3.1415926536
      vf = xkf(r)
	ukw=w/(xk*vf)
	zkw=xk/(2.d0*vf)
      xm1=(1.d0-(zkw-ukw)**2)*cdlog((zkw-ukw+1.d0)/(zkw-ukw-1.d0)) 
	xm2=(1.d0-(zkw+ukw)**2)*cdlog((zkw+ukw+1.d0)/(zkw+ukw-1.d0)) 
      wgranre=1.0d0+(4.d0*zkw+xm1+xm2)/(8.d0*pi*vf*zkw**3.d0)
      wgranim=0.d0
      cimg=dcmplx(0.d0,1.0d0)
	df=wgranre+cimg*wgranim
      elf=dimag(-1.0d0/(df))
      RETURN
      END
C  ********************************************************************
      SUBROUTINE LOS_MER(w,xk,r,VF0,LOSMER)
*     ----------------------------
C THIS SUBROUTINE FINDS   IM{-1/EM},
C WHERE EM IS THE MERMIN DIELECTRIC FUNCTION
C SEE THE PAPER FROM AGM ET AL. PHYS.STAT.SOL.
C
C Z     IS K/(2*KF) OR K/(2*VF)
C GAM   IS GAMMA0/EF
C GAMMA0 IS THE DAMPING OF THE ELECTRON GAS
C X     IS H*OMEGA/EF OR 2*OMEGA/VF/VF
C CHISQ IS 1/(PI*VF) IN ATOMIC UNITS

      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 LOSMER,K,GAMMA0
      COMPLEX*16 RSLT1,RSLT2,W,E,EPSM,CIMG,E0Z,PT2,PT3,PT4,E0ZCN
      

      PI=3.1415926536
      GAMMA0 = 0.001
C      WRITE(7,*)' DATOS EN LA SUBRUTINA LOS_MER'
C      WRITE(7,*)'K=',K
C      WRITE(7,*)'OMEGA=',OMEGA
C      WRITE(7,*)'VF0=',VF0
C      WRITE(7,*)'GAMMA0=',GAMMA0

C CIMG IS I
      CIMG=DCMPLX(0.0E0,1.0E0)

      vf = xkf(r)	
      Z     = K/(vf+vf)
      X     = (w+w)/vf/vf
      CHISQ = 1.E0/PI/vf
      GAM   = (GAMMA0+GAMMA0)/vf/vf

C  PT4 IS (Z+1)/(Z-1)
      PT4 = DCMPLX((Z+1.0D0)/(Z-1.0D0),0.0D0)

C E0Z IS EL(Z,0)-1, WHERE EL IS THE LINDHARD DIELECTRIC FUNCTION
      E0Z   = CHISQ*(4.0*Z+2.*(1.-Z*Z)*CDLOG(PT4))/(8.0*Z**3)
      E0ZCN = DCONJG(E0Z)
      E0Z   = .5D0*(E0Z+E0ZCN)
C W IS X+I*GAM
      omega  = DCMPLX(X,GAM)
      PT2 = omega/(4.*Z)

C IN SMALZ SUBROUTINE: Z AND W ARE INPUT, RSLT IS OUTPUT.
C RSLT IS THE LOGARITHM TERM IN THE LINDHARD DIELECTRIC FUNCTION
      CALL SMALZ(Z,RSLT1,-W)
      CALL SMALZ(Z,RSLT2,W)

C E IS EL(Z,W)-1, WHERE EL(Z,W) IS THE LINDHARD DIELECTRIC FUNCTION
C IN THIS CASE, IN THE LINDHARD FUNCTION WE MADE THE SUBSTITUTION
C OF OMEGA BY OMEGA+I*GAMMA0
C WITH THIS SUBSTITUTION THE USUAL PARAMETETRIZATION IN EL(Z,U)
C IN THIS CASE U=W/(4Z)
      E   = CHISQ*(4.*Z+(1.-(Z-PT2)**2)*RSLT1+
     .                  (1.-(Z+PT2)**2)*RSLT2)/(8.0*Z**3)
c      xloslind = DIMAG(-1.0D0/(E+1.0d0))
      PT3 = GAM*E/E0Z

C EPSM IS EM, THE MERMIN DIELECTRIC FUNCTION
      EPSM = 1.0 + W*E/(X+CIMG*PT3)
C LOSFUN1 IS IM(-1/EM)

      LOSMER = DIMAG(-1.0D0/(EPSM))
C
      RETURN
      END
***********************************************************************

      SUBROUTINE SMALZ(Z,RSLT,W)
*     --------------------------

      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 RSLT,CON2,W,ARG,PT2

*-------------------

      PT2 = W/(4.*Z)
      X   = DREAL(W)
      GAM = DIMAG(W)

      FAC   = (4.*Z*Z+X-4.*Z)
      DENOM = FAC*FAC+GAM*GAM
      EPS1  =  8.*Z*FAC/DENOM
      EPS2  = -8.*Z*GAM/DENOM
      EPS3  = 2.0*EPS1+EPS1*EPS1+EPS2*EPS2

      IF(DABS(EPS3)-0.2) 10, 10, 12

  10  CONTINUE
c      PRINT*,'SMALZ: SUM. Z, W =',Z,W
      SUM  = EPS3
      SMLT = -1.0
      FAC  = EPS3

      DO 11 I=2,60
      FAC  = FAC*EPS3
      SUM  = SUM+SMLT*FAC/DBLE(I)
      SMLT = -SMLT
  11  CONTINUE

      RLLPT = .5D0*SUM
      ARG   = DCMPLX(1.0+EPS1,EPS2)
      CON2  = CDLOG(ARG)
      CIMG  = DIMAG(CON2)
      RSLT  = DCMPLX(RLLPT,CIMG)
      RETURN

  12  CONTINUE
c      PRINT*,'SMALZ: EVALUATING LOGARITHM DIRECTLY',z,w
      RSLT = CDLOG((Z+PT2+1.0)/(Z+PT2-1.0))

      RETURN
      END
c  ************************************************************
      subroutine elfplasdeg(w,xk,r,elf)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 df
      common/plas/xkdau,vthau,xkfc,xnu,xdegc,xnetac,enl

      PI=3.1415926536
      vf = xkf(r)	
      ukw=w/(xk*vf)
      zkw=xk/(2.d0*vf)
      wgranre=1+(gdeg(ukw+zkw)-gdeg(ukw-zkw))/(4.d0*pi*vf*(zkw**3.d0))
c  **************************************************************
c       Aproximacion degenerada D>>1, e2D   , AristaPRA98
c  **************************************************************
      if((zkw+ukw).lt.1.d0) f2=2.d0*pi*ukw/(vf**2.d0)
      if(dabs(zkw-ukw).lt.1.d0.and.(zkw+ukw).gt.1.d0)
     1  f2=pi*(1.d0-((ukw-zkw)**2.d0))
      if(dabs(zkw-ukw).gt.1.d0) f2=0.d0
      wgranim=f2/(8.d0*pi*vf*zkw**3.d0)
c  **************************************************************
      df=dcmplx(wgranre,wgranim)
      elf=dimag(-1.0d0/df)
c      elf=dabs(dimag(-1.0d0/df)) !!VALOR ABSOLUTO!!
	
      RETURN
      END
C  **************************************************************
C                Functions from AristaPRA84
C  **************************************************************
      function ga(x)
      implicit double precision (a-h,o-z)
      common/gai/xx
      common/plas/xkdau,vthau,xkfc,xnu,xdegc,xnetac,enl
      external gaint,gdeg,gabq3
      xx=x
      eiy1=1.d-5
      eiy2=dsqrt((709.d0+xnetac)/xdegc)
      tol=1.d-4
      call gabq3(gaint,eiy1,eiy2,sum,tol,ier)
      ga=sum
      return
      end
C  **************************************************************
      function gaint(y)
      implicit double precision (a-h,o-z)
      common/plas/xkdau,vthau,xkfc,xnu,xdegc,xnetac,enl
      common/gai/x
      xdegc=(xkfc**2.d0)/(2.d0*vthau**2.d0)
      xnuc=finv12((2.d0*xdegc**1.5d0)/3.d0)
      xln=dlog(dabs((x+y)/(x-y)))
      gaint=y*xln/(dexp(xdegc*y**2.d0-xnuc)+1.d0)
      return
      end
C  **************************************************************
      function gdeg(x)
      implicit double precision (a-h,o-z)
      gdeg=x+(1-x**2.d0)*dlog(dabs((1+x)/(x-1)))/2
      return
      end
C  **************************************************************
C                      SLPA functions
C  **************************************************************
      function wavefun(r)
	  use commonarrays
C     This function takes the input values saved in karray, xarray and
C     carray and constructs the Harte-Fock wave function for the
C     selected shell
      implicit double precision (a-h,o-z)
	  integer m1,m2,l,n,ksize,s,i,j
	  common/hart/m1,m2,l
      j = 1
      ksize = size(karray)/2
      s = 0
      orb = 0.d0
      DO WHILE (j.le.ksize)
        k = karray(2*j)
        n = karray(2*j-1)
        DO i = 1, k	
c            write(*,*) 'n', n, 'l', l
            orb = orb + carray(i+s)*xi(r,n,l,i+s)
        END DO
        s = s + k
        j = j + 1
      END DO
      wavefun = 2.d0*(orb**2.d0)
      return
      end
C  **************************************************************
      function xkf(r)
      implicit double precision (a-h,o-z)
      pi=3.1415926536
        xkf = (3*wavefun(r)*pi**2.d0)**(1.d0/3.d0)
      return
      end
C  **************************************************************
	function xi(r,n,l,i)
	use commonarrays
	implicit double precision (a-h,o-z)
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
c      write(*,*) 'xa', xa
      Yl = sqrt((2.d0*l+1)/4.d0/pi)
      N10 = 1/sqrt(fact)*(2.d0*xa)**(n+0.5)
      xi = N10*(r**(n-1))*exp(-xa*r)*Yl
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
C                Inverse de la function de Fermi 1/2
C  **************************************************************
      FUNCTION FINV12(Y)
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
      FINV12 = LOG(ZNUM/ZDEN)
      RETURN
      ELSE IF(Y.LT.5.77072652741599997d0) THEN
      Y1=Y
      Y2 = Y1*Y
      Y3 = Y2*Y
      Y4 = Y3*Y
      Y5 = Y4*Y
      ZNUM=C0+C1*Y1+C2*Y2+C3*Y3+C4*Y4+C5*Y5
      ZDEN=D0+D1*Y1+D2*Y2+D3*Y3+D4*Y4+D5*Y5
      FINV12 = ZNUM / ZDEN
       RETURN
      ELSE IF(Y .LT. 0.598127953660604876d2) THEN
C TRANSFORMATION DE LA VARIABLE
      Z=1.d0/SQRT(Y)
      ZZ=((Z-YMIN)-(YMAX-Z))/(YMAX-YMIN)
      V(1)=1.d0
      V(2)=ZZ
      DO 1515 KK=3,10
1515  V(KK)=2*ZZ*V(KK-1)-V(KK-2)
      P=0.D0
      DO 1516 KK=10,1,-1
1516  P=P+E0(KK)*V(KK)
      FINV12=(1.d0/P**2)**(1./3.)
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
      FINV12 = ZNUM / ZDEN
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
   14 FORMAT(/5X,'*** LOW ACCURACY IN SUBROUTINE GABQ.')
      WRITE(6,11) XL,XU,TOL
      WRITE(6,15) SUM,ERR
   15 FORMAT(5X,'SUM =',1PD19.12,', ERR =',D8.1//)
      IER=1
      RETURN
      END
C  **************************************************************
C                      SUBROUTINE GABQ3
C  **************************************************************
      SUBROUTINE GABQ3(FCT,XL,XU,SUM,TOL,IER)
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

