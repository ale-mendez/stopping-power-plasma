** WE FIND THE PROTON STOPPING POWER (Lindhard 1954)
	PROGRAM SLPA
      
      IMPLICIT REAL*8 (A-H,O-Z)
      character(len=30) input,filename
      character(len=:),allocatable:: str
      character*3 sorb,ssorb,target
      character*20 enerfile,wfolder
      parameter(NVMX=28)
      parameter(NOR=30)

      common/plas/e_thr,xnocc,ja
      COMMON/CINDIV/VV
      common/norma/xnormval
      common/rbck/rmax
      common/tolbck/TOLK,TOLW,TOLR

      common/target/noc(NOR),gtype
      common/gridtype/igrid,maxpts
      common/orbitalbck/nmax,sorb(NOR)
      common/namesbck/enerfile,wfolder
      common/bindener/enl(NOR)
      common/meanval/rnl(-1:3,NOR),rmaximo(NOR)

      EXTERNAL FUN1,GABQ

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
      TOLW = 1.0e-2
      TOLR = 1.0e-2

      gtype = 0
      igrid = 1

c.... Read input file
      open(10,file='slpa_ll_v2.inp',status='old')
      read(10,*) vmin
      read(10,*) vmax
      read(10,*) target
      read(10,*) enerfile
      read(10,*) wfolder
      do 100 ja = 1,NOR
        read(10,*,end=101) ssorb,nnocc
        sorb(ja) = ssorb
        noc(ja) = nnocc
100   continue

c.... Define number of orbitals to compute stopping
101   nmax = ja - 1
      close(10)

c.... Read target input data and compute wave mean values
      call radial_grid()
      call input_target()
      call rasymptotic()

      ja = 1
      e_thr = enl(ja)
      xnocc = noc(ja)
      rmax = rmaximo(ja)
      print*,sorb(ja),xnocc,e_thr,rmax,rnl(0,ja)

***********************************************************************
      str = 'SLPA-LL_'//trim(target)//trim(sorb(ja))//'.dat'
      open(7,file=str,status='unknown')
      write(7,1000) trim(target),trim(sorb(ja))
      write(7,1002)
      close(7)

1000  format('# Proton stopping power using the SLPA method'/,
     |       '# with the Levine-Lindhard dielectric function'/,
     |       '# Target: ',a,"(",a,')'/)
1002  format('# V(au)',3X,'SP(au/atom)',4X,'SP(10^-15eVcm2/atom)')

c.... compute stopping over a given distribution of velocities
      do 200 iv=1,NVMX
         v=vp_au(iv)
         VV=v

c....    add the minimum and maximum velocity given
         if ((vp_au(iv+1).gt.vmin).and.(v.lt.vmin)) v = vmin
c         if (iv.gt.1) then
           if ((v.gt.vmax).and.(vp_au(iv-1).lt.vmax)) v = vmax
c         endifZ
         if ((v.lt.vmin).or.(v.gt.vmax)) go to 200

c....    define integration limits for k variable
         EIK1 = 1.D-5
         if (v.lt.1) epsmax = 20.d0
         if (v.ge.1) epsmax = 5.d0
         EIK2 = epsmax*v

* STOPPING
         CALL GABQ(FUN1,EIK1,EIK2,SUM,TOLK,IER)
         spp = 2.d0/PI/V**2 * SUM
         write(*,*)'v=',v,'  sp',spp
         OPEN(7,FILE=str,status='old',access='append')
         WRITE(7,1001) v,spp,spp*conv_au_to_eVcm2
         close(7)

 200  CONTINUE
 1001 FORMAT(2X,F6.3,3X,1pE12.6,3x,e12.6,3x,e12.6)
      stop
      END

C  **************************************************************
      FUNCTION FUN1(K)
*     ----------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 K,KK
c
      COMMON /CINDIV/ V
      COMMON /CK/ KK,W
      EXTERNAL GABQ1,FW
      common/tolbck/TOLK,TOLW,TOLR
c
      KK  = K
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
	EXTERNAL FR

      PI=3.1415926536
	WW=W
*  Integrate over r from zero to rmax
	R1 = 1.e-5
	R2 = 2.9
      CALL GABQ2(FR,R1,R2,RES,TOLR,IER)
      RLPA = 4.d0*pi*RES
c      CALL INTEG4(R1,R2,CFR,CRES,ERR,0,0)
c      RLPA = 4.d0*pi*real(CRES)
      FW = W*RLPA

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
      subroutine elflev(w,xk,rr,celf)
      IMPLICIT REAL*8(A-B,D-H,K,O-Z), COMPLEX*16(C)

      common/plas/enl,xnocc,ja
      DATA C1/(1.D0,0.D0)/, C0/(0.D0,0.D0)/, CI/(0.D0,1.D0)/
      external lindhard_full, density_spline

      pi = 3.1415926536

      if(dabs(w).ge.enl) then
         wg = sqrt(w**2.d0-enl**2.d0)
         density_tot = xnocc * density_spline(ja,rr)
c         write(555,*) w,xk,rr,xnocc,enl,density_tot
         vf = (3*pi**2.d0*density_tot)**(1.d0/3.d0)
         call lindhard_full(wg,xk,vf,celf)
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

c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================

c
c***********************************************************************
      subroutine radial_grid()
c***********************************************************************
c
      implicit real*8 (a-b,d-h,o-z)
      parameter(NGP=500)

      common/radialgrid/rs(NGP),rsp(NGP),h
      common/gridtype/igrid,maxpts

      data rzero,half,one,two/0.0d0,0.5d0,1.0d0,2.0d0/

      rs = rzero

c.... Radial grid used in radial integrals
c     r(i) can not be 0 at the origin because log(r) is computed 
      if (igrid.eq.1) then
         rmax = 100.0d0
         r0 = 0.01d0
         h = dlog(rmax/r0+1)/(NGP-1)
         maxpts = NGP
         do 100 i = 1,NGP
            t = h*i
            rsp(i) = r0*dexp(t) 
            rs(i) = rsp(i) - r0
100      continue
      else
         print*,'not implemented igrid .ne. 1'
         stop
c      elseif (igrid.eq.0) then
c         maxpts = NRMAX_UNIVERSAL
c         do 200 i = 1,maxpts
c            rs(i) = RRR0(i)
c 200      continue
      endif
      
      return
      end
c
c***********************************************************************
      subroutine rasymptotic()
c***********************************************************************
c
c     Compute integration radii for bound electrons
c
      implicit real*8 (a-b,d-h,o-z)
      character*3 sorb
      parameter(NOR=30,NGP=500)

      common/orbitalbck/nmax,sorb(NOR)
      common/waves/r(NGP,NOR),p(NGP,NOR),npts(NOR)
      common/meanval/rnl(-1:3,NOR),rmaximo(NOR)

      data eps/1E-4/

      do 100 ja = 1,nmax
c.... use 1st criteria: p(rmax) <= eps = 1E-4 (DEFAULT)
         do 200 i=1,npts(ja)
            pp = dabs(p(i,ja))
            rr = r(i,ja)
            if (rr.ge.rnl(1,ja).and.pp.le.eps) goto 201
200      continue
201   rasympt = rr
      rmaximo(ja) = rasympt
c      print*,'rasympt=',rasympt
100   continue

9     return
      end
c
c***********************************************************************
      subroutine input_target()
c***********************************************************************
c

c.... Read binding energy values
      call binding_energies()

c.... Read orbital wavefunctions
      call wavefunctions()

c.... Define weight for Auger 
c      call define_weight()

c.... Define energy gap
c      call define_wgap()

      return
      end
c
c***********************************************************************
      subroutine binding_energies()
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      character*3 header
      character*3 sorb,ssorb
      character*20 enerfile,wfolder
      character(len=:),allocatable:: str
      parameter(NOR=30)

      common/bindener/enl(NOR)
      common/orbitalbck/nmax,sorb(NOR)
      common/target/noc(NOR),gtype
      common/namesbck/enerfile,wfolder

      data rzero,half,one,two/0.d0,0.5d0,1.0d0,2.0d0/

      enl = rzero
      str = trim(enerfile)
      open(unit=30,file=str,status='unknown')

c...  Check headers
      read(30,*) header
      if (header(1:3).eq.'Orb'
     |      .or.header(1:2).eq.'nl'
     |      .or.header(1:1).eq.'#') then
         continue
      else
         rewind(30)
      endif

c...  Hartree-Fock (Fischer) - energies in Hartree
      if (gtype.eq.0) then
         do 100 i=1,NOR
            read(30,*,end=300) ssorb,enerEh
c            print*,ssorb,enerEh
c     match read values with input orbitals
            ii=0
            do 110 ja=1,nmax
               if (ssorb.eq.sorb(ja)) enl(ja) = dabs(enerEh)
110         continue
100      continue
c...  Autostructure (Model potential + CI)
      else if (gtype.eq.1) then
         print*, 'not implemented yet'
         stop
c...  Hullac (perturbative Dirac) - energies in Hartree
      else if (gtype.eq.2) then
         do 200 i=1,NOR
            read(30,*,end=300) ssorb,enerEh
c     match read values with input orbitals
            do 210 ja=1,NOR
               if (ssorb.eq.sorb(ja)) enl(ja) = -enerEh
210         continue
200      continue
      endif

300   nmax = ja - 1
      close(unit=30)

      return 
      end
c
c***********************************************************************
      subroutine wavefunctions()
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      character*3,sorb
      parameter (NOR=30)

      common/orbitalbck/nmax,sorb(NOR)
      common/target/noc(NOR),gtype
      common/meanval/rnl(-1:3,NOR),rmaximo(NOR)

      call readwavefun()
      call wavecubicspline()
      call checkwavespline()

c       ELECT = 0.D0
c       RR2 = 0.D0
c       DO 100 J = 1,NMAX
c          RR2 = RR2 + rnl(2,J)*noc(J)
c          ELECT = ELECT + noc(J)
c 100   CONTINUE
c       RR2 = RR2/ELECT

      return 
      end
c
c***********************************************************************
      subroutine readwavefun()
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      character*3 header
      character*22 fname
      character*3 sorb,ssorb
      character*20 enerfile,wfolder
      character(len=:),allocatable:: str

      parameter(NGP=500,NOR=30)

      common/target/noc(NOR),gtype
      common/waves/r(NGP,NOR),p(NGP,NOR),npts(NOR)
      common/bindener/enl(NOR)
      common/orbitalbck/nmax,sorb(NOR)
      common/namesbck/enerfile,wfolder

      data rzero,half,one,two/0.0d0,0.5d0,1.0d0,2.0d0/

      do 100 ja = 1,nmax
c.... Initiate wavefunction matrix
         do 110 i = 1,NGP
            p(i,ja) = rzero
110      continue
c.... Open file and read wavefunction
         ssorb = sorb(ja)
         str = trim(wfolder)//'/'//'wave'//trim(ssorb)//'.dat'          ! cambiar en WINDOWS
         open(unit=20,file=str)
c...  Check headers
         read(20,*) header
         if (header(1:3).eq.'Orb'
     |      .or.header(1:1).eq.'r'
     |      .or.header(1:2).eq.'nl'
     |      .or.header(1:1).eq.'#') then
            continue
         else
            rewind(20)
         endif
c.... Include zero point in grid and wavefunction
         r(1,ja) = rzero
         p(1,ja) = rzero
c.... Read wavefunctions
         do 120 i = 2,NGP
c...  Hullac (perturbative Dirac)
            if (gtype.eq.2) then
               read(20,*,end=130) x,pp,qq
               y = dsqrt(pp*pp+qq*qq)
c...  Hartree-Fock (Fischer) / Autostructure (Model potential + CI)
            else
               read(20,*,end=130) x,y
               if (abs(y).lt.1E-16) y = rzero
c               write(40+ja,*) x,y
            endif
            r(i,ja) = x
            p(i,ja) = y
120      continue
130      npts(ja) = i-1
         close(unit=20)
100   continue

      return
      end
c
c***********************************************************************
      subroutine wavecubicspline()
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      character*3,sorb

      parameter(NGP=500,NOR=30)

      common/orbitalbck/nmax,sorb(NOR)
      common/waves/r(NGP,NOR),p(NGP,NOR),npts(NOR)
      common/cspline/A0(NGP,NOR),B0(NGP,NOR),D0(NGP,NOR),E0(NGP,NOR)

      dimension rr(NGP),wave(NGP)
      dimension Ai(NGP),Bi(NGP),Di(NGP),Ei(NGP)

      data rzero,half,one,two/0.0d0,0.5d0,1.0d0,2.0d0/

      A0 = rzero
      B0 = rzero
      D0 = rzero
      E0 = rzero
      rr = rzero
      wave = rzero
      do 100 ja = 1,nmax
         ipts = npts(ja)
c.... Initiate cubic spline matrix and radii/wave vectors 
         do 110 i = 1,ipts
            rr(i) = r(i,ja)
            wave(i) = p(i,ja)
110      continue
c.... Make cubic spline of radii/wave vectors
         S1 = rzero
         SN = rzero
         call spline(rr,wave,Ai,Bi,Di,Ei,S1,SN,ipts)
         do 120 k = 1,ipts
           A0(k,ja) = Ai(k)
           B0(k,ja) = Bi(k)
           D0(k,ja) = Di(k)
           E0(k,ja) = Ei(k)
c          write(50+ja,*) r(k,ja), A0(k,ja)
120      continue
100   continue

      return
      end
c
c***********************************************************************
      subroutine checkwavespline()
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      character*3,sorb

      parameter(NGP=500,NOR=30)

      common/orbitalbck/nmax,sorb(NOR)
      common/waves/r(NGP,NOR),p(NGP,NOR),npts(NOR)
      common/cspline/A0(NGP,NOR),B0(NGP,NOR),D0(NGP,NOR),E0(NGP,NOR)
      common/radialgrid/rs(NGP),rsp(NGP),h
      common/gridtype/igrid,maxpts
      common/meanval/rnl(-1:3,NOR),rmaximo(NOR)
c      common/meanval/rnl(-1:3,NOR),irmax(NOR),rscale,RR2,maxir

      dimension rr(NGP),ws(NGP)

      data rzero,half,one,two/0.0d0,0.5d0,1.0d0,2.0d0/

      do 200 ja = 1,nmax
         do 210 i = 1,NGP
            rr(i) = r(i,ja)
            ws(i) = rzero
210      continue
         ipts = npts(ja)
         do 220 i = 1,maxpts
            ri = rs(i)
            call findi(rr,ri,ipts,k)
c.... If upper bound, use asymptotic condition (A,B,C,D=0)
            if (k.eq.ipts-1) k = ipts
            ws(i) = A0(k,ja)+ri*(B0(k,ja)+ri*(D0(k,ja)+ri*E0(k,ja)))
c            write(60+ja,*) ri,ws(i) ! check cubicsplined wavefunctions
220      continue

c         if (igrid.eq.0) call meanvalues_trap(ja,h,rs,ws)
         if (igrid.eq.1) call meanvalues_quad(ja,h,rs,rsp,ws)
         write(111,*) ja,sorb(ja),rnl(0,ja)
200   continue

      return
      end

c
c***********************************************************************
      subroutine meanvalues_quad(ja,h,r,rp,ws)
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
      parameter(NGP=500,NOR=30)

      common/meanval/rnl(-1:3,NOR),rmaximo(NOR)

      dimension r(NGP),rp(NGP),ws(NGP)
      dimension u(NGP)

      data rzero,half,one,two/0.0d0,0.5d0,1.0d0,2.0d0/

c.... Compute r^k integrals
      do k = -1,3
         u(1) = rzero
         do i = 2,NGP
            if (k.ge.0) rk = r(i)**k
            if (k.lt.0) rk = one/r(i)
            u(i) = ws(i)*ws(i)*rp(i)*rk
         enddo
         rkint = rint(u,1,NGP,7,h)
         rnl(k,ja) = rkint
      enddo

      return
      end
c
c***********************************************************************
      function density_spline(ja, rr)
c***********************************************************************
c
      implicit real*8(a-b,d-h,o-z)
      implicit double complex (c)
      parameter(NGP=500,NOR=30)

      common/waves/r(NGP,NOR),p(NGP,NOR),npts(NOR)
      common/cspline/A0(NGP,NOR),B0(NGP,NOR),D0(NGP,NOR),E0(NGP,NOR)

      dimension rj(NGP)

      DATA PI/3.14159265358979324/
      data rzero,half,one,two,four/0.0d0,0.5d0,1.0d0,2.0d0,4.0d0/

      if (rr.eq.0) then
         den = rzero
         goto 200
      endif
      ipts = npts(ja)
      do 110 i = 1,ipts
         rj(i) = r(i,ja)
110   continue
      call findi(rj,rr,ipts,k) 
      u = A0(k,ja)+rr*(B0(k,ja)+rr*(D0(k,ja)+rr*E0(k,ja)))
      wave = u/rr
      density_spline = wave*wave/(four*PI)            ! => (R Y00)**2 = (P/r Y00)**2

200   return
      end
c
c***********************************************************************
      function rint (f,na,nb,nq,h)
c***********************************************************************
c
c  this program calculates the integral of the function f from point na
c  to point nb using a nq points quadrature ( nq is any integer between
c  1 and 14 ).  h is the grid size.
c                                      written by c. c. j. roothaan
c
      implicit real*8(a-h,o-z)
      dimension c(105),c1(78),c2(27),d(14),f(nb)
      equivalence (c1(1),c(1)) ,(c2(1),c(79))
      data c1/
     a 1.,
     b 2.,1.,
     c 23.,28.,9.,
     d 25.,20.,31.,8.,
     e 1413.,1586.,1104.,1902.,475.,
     f 1456.,1333.,1746.,944.,1982.,459.,
     g 119585.,130936.,89437.,177984.,54851.,176648.,36799.,
     h 122175.,111080.,156451.,46912.,220509.,29336.,185153.,35584.,
     i 7200319.,7783754.,5095890.,12489922.,-1020160.,16263486.,261166.,
     i 11532470.,2082753.,
     j 7305728 ,6767167.,9516362.,1053138.,18554050.,-7084288.,
     j 20306238.,-1471442.,11965622.,2034625.,
     k 952327935.,1021256716.,636547389.,1942518504.,-1065220914.,
     k 3897945600.,-2145575886.,3373884696.,-454944189.,1637546484.,
     k 262747265.,
     l 963053825.,896771060.,1299041091.,-196805736.,3609224754.,
     l-3398609664.,6231334350.,-3812282136.,4207237821.,-732728564.,
     l 1693103359., 257696640./
      data c2 / 5206230892907.,5551687979302.,3283609164916.,
     m 12465244770050.,-13155015007785.,39022895874876.,
     m-41078125154304.,53315213499588.,-32865015189975.,28323664941310.,
     m-5605325192308.,9535909891802.,1382741929621.,
     n 5252701747968.,4920175305323.,7268021504806.,-3009613761932.,
     n 28198302087170.,-41474518178601.,76782233435964.,
     n-78837462715392.,81634716670404.,-48598072507095.,
     n 34616887868158.,-7321658717812.,9821965479386.,1360737653653./
      data d/2.,2.,24.,24.,1440.,1440.,120960.,120960.,7257600.,
     a  7257600.,958003200.,958003200.,5230697472000.,5230697472000./

      a = 0.0
      l = na
      m = nb
      i = nq*(nq+1)/2
      do 100 j = 1,nq
         a = a + c(i)*( f(l) + f(m) )
         l = l + 1
         m = m - 1
         i = i - 1
 100  continue
      a = a/d(nq)
      do 200 n = l,m
        a = a + f(n)
 200  continue
      rint = a*h
      return
      end
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SUBROUTINE SPLINE(X,Y,A,B,C,D,S1,SN,N)                            
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C                                                                       
C      CUBIC SPLINE INTERPOLATION BETWEEN TABULATED DATA.               
C   INPUT:                                                              
C     X(I) (I=1, ...,N) ........ GRID POINTS.                           
C                    (THE X VALUES MUST BE IN INCREASING ORDER).        
C     Y(I) (I=1, ...,N) ........ CORRESPONDING FUNCTION VALUES.         
C     S1,SN ..... SECOND DERIVATIVES AT X(1) AND X(N).                  
C            (THE NATURAL SPLINE CORRESPONDS TO TAKING S1=SN=0).        
C     N ........................ NUMBER OF GRID POINTS.                 
C      THE INTERPOLATING POLYNOMIAL IN THE I-TH INTERVAL, FROM          
C   X(I) TO X(I+1), IS                                                  
C            PI(X) = A(I)+X*(B(I)+X*(C(I)+X*D(I)))                      
C   OUTPUT:                                                             
C     A(I),B(I),C(I),D(I) ...... SPLINE COEFFICIENTS.                   
C                                                                       
C      REF.: M.J. MARON, 'NUMERICAL ANALYSIS: A PRACTICAL               
C            APPROACH', MACMILLAN PUBL. CO., NEW YORK 1982.             
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      PARAMETER(NMX=500)                               
      DIMENSION X(NMX),Y(NMX),A(NMX),B(NMX),C(NMX),D(NMX)                           
      IF(N.LT.4) THEN                                                   
      WRITE(6,10) N                                                     
   10 FORMAT(5X,'SPLINE INTERPOLATION CANNOT BE PERFORMED WITH',        
     1I4,' POINTS. STOP.')                                              
      STOP                                                              
      ENDIF                                                             
      N1=N-1                                                            
      N2=N-2                                                            
C  ****  AUXILIARY ARRAYS H(=A) AND DELTA(=D).                          
      DO 1 I=1,N1                                                       
      IF(X(I+1)-X(I).LT.1.0D-10) THEN                                   
      WRITE(6,11)                                                       
   11 FORMAT(5X,'SPLINE X VALUES NOT IN INCREASING ORDER. STOP.')       
      STOP                                                              
      ENDIF                                                             
      A(I)=X(I+1)-X(I)                                                  
      D(I)=(Y(I+1)-Y(I))/A(I)
    1 CONTINUE
C  ****  SYMMETRIC COEFFICIENT MATRIX (AUGMENTED).                      
      DO 2 I=1,N2                                                       
      B(I)=2.0D0*(A(I)+A(I+1))                                          
      K=N1-I+1                                                          
      D(K)=6.0D0*(D(K)-D(K-1))
    2 CONTINUE
      D(2)=D(2)-A(1)*S1                                                 
      D(N1)=D(N1)-A(N1)*SN                                              
C  ****  GAUSS SOLUTION OF THE TRIDIAGONAL SYSTEM.                      
      DO 3 I=2,N2                                                       
      R=A(I)/B(I-1)                                                     
      B(I)=B(I)-R*A(I)                                                  
      D(I+1)=D(I+1)-R*D(I)
    3 CONTINUE
C  ****  THE SIGMA COEFFICIENTS ARE STORED IN ARRAY D.                  
      D(N1)=D(N1)/B(N2)                                                 
      DO 4 I=2,N2                                                       
      K=N1-I+1                                                          
      D(K)=(D(K)-A(K)*D(K+1))/B(K-1)
    4 CONTINUE
      D(N)=SN                                                           
C  ****  SPLINE COEFFICIENTS.                                           
      SI1=S1                                                            
      DO 5 I=1,N1                                                       
      SI=SI1                                                            
      SI1=D(I+1)                                                        
      H=A(I)                                                            
      HI=1.0D0/H                                                        
      A(I)=(HI/6.0D0)*(SI*X(I+1)**3-SI1*X(I)**3)                        
     1    +HI*(Y(I)*X(I+1)-Y(I+1)*X(I))                                 
     2    +(H/6.0D0)*(SI1*X(I)-SI*X(I+1))                               
      B(I)=(HI/2.0D0)*(SI1*X(I)**2-SI*X(I+1)**2)                        
     1    +HI*(Y(I+1)-Y(I))+(H/6.0D0)*(SI-SI1)                          
      C(I)=(HI/2.0D0)*(SI*X(I+1)-SI1*X(I))                              
      D(I)=(HI/6.0D0)*(SI1-SI)
    5 CONTINUE
      RETURN                                                            
      END                                                               

C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SUBROUTINE FINDI(X,XC,N,I)                                        
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C                                                                       
C      FINDS THE INTERVAL (X(I),X(I+1)) CONTAINING THE VALUE XC.        
C   INPUT:                                                              
C     X(I) (I=1, ...,N) ........ GRID POINTS.                           
C                    (THE X VALUES MUST BE IN INCREASING ORDER).        
C     XC ....................... POINT TO BE LOCATED.                   
C     N ........................ NUMBER OF GRID POINTS.                 
C   OUTPUT:                                                             
C     I ........................ INTERVAL INDEX.                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      PARAMETER(NMX=500)                               
      DIMENSION X(NMX)
                                                          
      IF(XC.GT.X(N)) THEN                                               
      I=N-1                                                             
      RETURN                                                            
      ENDIF                                                             
      IF(XC.LT.X(1)) THEN                                               
      I=1                                                               
      RETURN                                                            
      ENDIF                                                             
      I=1                                                               
      I1=N                                                              
    1 IT=(I+I1)/2                                                       
      IF(XC.GT.X(IT)) I=IT                                              
      IF(XC.LE.X(IT)) I1=IT                                             
      IF(I1-I.GT.1) GO TO 1                                             
      RETURN                                                            
      END                                                               

c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================



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
      HO=H
    2 continue
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
   14 FORMAT(/5X,'*** LOW ACCURACY IN SUBROUTINE GABQ (k integral).')
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
      D=D+W(I1)*(FCT(B+C)+FCT(B-C))
    1 continue
      SUM=D*A
C  ****  ADAPTIVE BIPARTITION SCHEME.
      ICALL=NP2
      LH=1
      S(1)=SUM
      L(1)=1
      HO=H
    2 continue
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
   14 FORMAT(/5X,'*** LOW ACCURACY IN SUBROUTINE GABQ1 (w integral).')
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
      D=D+W(I1)*(FCT(B+C)+FCT(B-C))
    1 continue
      SUM=D*A
C  ****  ADAPTIVE BIPARTITION SCHEME.
      ICALL=NP2
      LH=1
      S(1)=SUM
      L(1)=1
      HO=H
    2 continue
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
   14 FORMAT(/5X,'*** LOW ACCURACY IN SUBROUTINE GABQ2 (R integral).')
      WRITE(6,11) XL,XU,TOL
      WRITE(6,15) SUM,ERR
   15 FORMAT(5X,'SUM =',1PD19.12,', ERR =',D8.1//)
      IER=1
      RETURN
      END

