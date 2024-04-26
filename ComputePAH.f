	module optECdata
	IMPLICIT NONE
	integer nEg_optEC,nr_optEC,nl_optEC
	parameter(nEg_optEC=14,nr_optEC=7,nl_optEC=1000)
	real*8 l_optEC(nl_optEC),e1_optEC(nl_optEC,nEg_optEC,nr_optEC),e2_optEC(1000,nEg_optEC,nr_optEC)
	real*8 Eg_optEC(nEg_optEC),r_optEC(nr_optEC)
	end module optECdata

	subroutine Compute_optEC(Ca_cont,Cs_cont,computelamcloud)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	logical computelamcloud(nlam),ionized
	real*8 Ca_cont(nlam),Cs_cont(nlam),Ca(nlam),Cs(nlam),Mc,Nc,Mpart
	parameter(Mc=2d-23) ! in grams

	if(mixrat_optEC.le.0d0) return

	Mpart=1.5d0*(4d0*pi*(rad_optEC*1d-4)**3)/3d0
	Nc=Mpart/Mc
	
	Ca=0d0
	Cs=0d0
	call Make_optEC(lam*1d4,Ca,Cs,rad_optEC,Eg_optEC,nlam,computelamcloud)

	Ca_cont=Ca_cont+Ca*mixrat_optEC/Nc
	Cs_cont=Cs_cont+Cs*mixrat_optEC/Nc

	return
	end


	subroutine Init_optEC()
	use optECdata
	IMPLICIT NONE
	character*100 filebase,filename,homedir
	character*4 ext(7)
	integer i,j,k
	real*8 dummy
	external Tholin
	logical truefalse

	call getenv('HOME',homedir)
	filebase=trim(homedir) // '/ARCiS/Data/optEC/'
	r_optEC(1:nr_optEC) = (/ 0.33,0.5,1.0,3.0,10.0,30.0,100.0 /)
	r_optEC=r_optEC*1d-3
	ext(1:nr_optEC) = (/ "0_33", "0_5 ", "1   ", "3   ", "10  ", "30  ", "100 " /)
	Eg_optEC(1:nEg_optEC) = (/ -0.1,0.0,0.1,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.67 /)

	do i=1,nr_optEC
		filename=trim(filebase) // "n" // trim(ext(i)) // "nm.dat"
		inquire(file=filename,exist=truefalse)
		if(.not.truefalse) then
			write(*,'("WARNING: optEC data not found")')
			write(*,'("         optEC opacities set to tholin")')
			write(*,'("==================================================================")')
			flush(9)
			do j=1,nl_optEC
				l_optEC(j)=10d0**(log10(0.2)+log10(1000.0/0.2)*real(j-1)/real(nl_optEC-1))
			enddo
			call RegridDataLNK(Tholin,l_optEC,e1_optEC(1:nl_optEC,1,1),e2_optEC(1:nl_optEC,1,1),nl_optEC,.true.)
			do j=1,nl_optEC
				e1_optEC(j,1:nEg_optEC,1:nr_optEC)=e1_optEC(j,1,1)
				e2_optEC(j,1:nEg_optEC,1:nr_optEC)=e2_optEC(j,1,1)
			enddo
			return
		endif
		open(unit=25,file=filename,FORM="FORMATTED")
		filename=trim(filebase) // "k" // trim(ext(i)) // "nm.dat"
		open(unit=26,file=filename,FORM="FORMATTED")
		do j=1,nl_optEC
			read(25,*) l_optEC(j),dummy,e1_optEC(j,1:nEg_optEC,i)
			read(26,*) l_optEC(j),dummy,e2_optEC(j,1:nEg_optEC,i)
		enddo
		close(unit=25)
		close(unit=26)
	enddo

	return
	end

	subroutine Make_optEC(lam,Cabs,Csca,rad,Eg,nlam,computelam)
	use optECdata
	IMPLICIT NONE
	integer nlam,ir,iEg,ilam
	real*8 Cabs(nlam),Csca(nlam),rad,Eg,QEX,QSC,QAB,G,lam(nlam),pi
	parameter(pi=3.14159265358979323846264338328d0)
	real*8 wr1,wr2,wEg1,wEg2,e1x(nlam),e2x(nlam),e1(nl_optEC),e2(nl_optEC)
	logical computelam(nlam)

	call RefInd_optEC(lam,e1x,e2x,rad,Eg,nlam)

	do ilam=1,nlam
		if(computelam(ilam)) then
			call Q_MIE(e1x(ilam),e2x(ilam),lam(ilam),rad,QEX,QSC,QAB,G)
			Cabs(ilam)=1d-8*(qab*pi*rad**2)
			Csca(ilam)=1d-8*(qsc*pi*rad**2)
		endif
	enddo
	
	return
	end
	


	subroutine RefInd_optEC(lam,e1x,e2x,rad,Eg,nlam)
	use optECdata
	IMPLICIT NONE
	integer nlam,ir,iEg,ilam
	real*8 rad,Eg,QEX,QSC,QAB,G,lam(nlam),pi
	parameter(pi=3.14159265358979323846264338328d0)
	real*8 wr1,wr2,wEg1,wEg2,e1x(nlam),e2x(nlam),e1(nl_optEC),e2(nl_optEC)

	if(rad.le.r_optEC(1)) then
		ir=1
		wr1=1d0
		wr2=0d0
	else if(rad.ge.r_optEC(nr_optEC)) then
		ir=nr_optEC-1
		wr1=0d0
		wr2=1d0
	else
		do ir=1,nr_optEC-1
			if(rad.ge.r_optEC(ir).and.rad.lt.r_optEC(ir+1)) exit
		enddo
		wr1=log(r_optEC(ir+1)/rad)/log(r_optEC(ir+1)/r_optEC(ir))
		wr2=1d0-wr1
	endif

	if(Eg.le.Eg_optEC(1)) then
		iEg=1
		wEg1=1d0
		wEg2=0d0
	else if(Eg.ge.Eg_optEC(nEg_optEC)) then
		iEg=nEg_optEC-1
		wEg1=0d0
		wEg2=1d0
	else
		do iEg=1,nEg_optEC-1
			if(Eg.ge.Eg_optEC(iEg).and.Eg.lt.Eg_optEC(iEg+1)) exit
		enddo
		wEg1=(Eg_optEC(iEg+1)-Eg)/(Eg_optEC(iEg+1)-Eg_optEC(iEg))
		wEg2=1d0-wEg1
	endif

	e1(1:nl_optEC)=wr1*wEg1*e1_optEC(1:nl_optEC,iEg,ir)+
     &			wr1*wEg2*e1_optEC(1:nl_optEC,iEg+1,ir)+
     &			wr2*wEg1*e1_optEC(1:nl_optEC,iEg,ir+1)+
     &			wr2*wEg2*e1_optEC(1:nl_optEC,iEg+1,ir+1)
	e2(1:nl_optEC)=wr1*wEg1*e2_optEC(1:nl_optEC,iEg,ir)+
     &			wr1*wEg2*e2_optEC(1:nl_optEC,iEg+1,ir)+
     &			wr2*wEg1*e2_optEC(1:nl_optEC,iEg,ir+1)+
     &			wr2*wEg2*e2_optEC(1:nl_optEC,iEg+1,ir+1)
	call regridarray(l_optEC,e1,nl_optEC,lam,e1x,nlam)
	call regridarray(l_optEC,e2,nl_optEC,lam,e2x,nlam)

	return
	end



	subroutine ComputePAH(Ca_cont,Cs_cont,computelamcloud)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	logical computelamcloud(nlam),ionized
	real*8 HC,Ca_cont(nlam),Cs_cont(nlam),Ca(nlam),Cs(nlam),Mc
	parameter(Mc=2d-23) ! in grams
	
	if(mixrat_PAH.le.0d0) return
	
	if(nC_PAH.lt.25) then
		HC=0.5d0
	else if(nC_PAH.lt.100) then
		HC=0.5/sqrt(nC_PAH/25d0)
	else
		HC=0.25d0
	endif
	ionized=.false.

	Ca=0d0
	Cs=0d0
	call MakePAH(lam*1d4,Ca,Cs,nC_PAH,HC,nlam,ionized,computelamcloud)

	Ca_cont=Ca_cont+Ca*Mc*nC_PAH*mixrat_PAH
	Cs_cont=Cs_cont+Cs*Mc*nC_PAH*mixrat_PAH

	return
	end
	

	
	subroutine MakePAH(lam,Cabs,Csca,Nc,HC,nlam,ionized,computelam)
	IMPLICIT NONE
	integer nlam,i,j,ilam,nj
	parameter(nj=30)
	real*8 lam(nlam),Cabs(nlam),Csca(nlam),Nc,HC,x,cutoffPAH,Mc
	real*8 e1x(nlam),e2x(nlam),e1y(nlam),e2y(nlam),CabsGra
	real*8 a,fpah,qgra
	real*8 r,QEX,QSC,QAB,G
	character*100 filename
	logical ionized,computelam(nlam)
	real*8 lj(nj),gj(nj),sjn(nj),sji(nj),S(nj),pi
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(Mc=2d-23) ! in grams
	data (lj(j),j=1,30) / 0.0722d0, 0.2175d0,1.050d0,1.260d0,1.905d0,3.300d0,5.270d0,5.700d0,
     &	6.220d0,6.690d0,7.417d0,7.598d0,7.850d0,8.330d0,8.610d0,10.68d0,11.23d0,11.33d0,11.99d0,12.62d0,
     &	12.69d0,13.48d0,14.19d0,15.90d0,16.45d0,17.04d0,17.375d0,17.87d0,18.92d0,15.0d0 /
	data (gj(j),j=1,30) / 0.195d0,0.217d0,0.055d0,0.11d0,0.09d0,0.012d0,0.034d0,0.035d0,0.030d0,
     &	0.070d0,0.126d0,0.044d0,0.053d0,0.052d0,0.039d0,0.020d0,0.012d0,0.032d0,0.045d0,0.042d0,0.013d0,
     &	0.040d0,0.025d0,0.020d0,0.014d0,0.065d0,0.012d0,0.016d0,0.10d0,0.8d0 / 
	data (sjn(j),j=1,30) / 7.97d7,1.23d7,0.0d0,0.0d0,0.0d0,394.0d0,2.5d0,4.0d0,29.4d0,7.35d0,
     &	20.8d0,18.1d0,21.9d0,6.94d0,27.8d0,0.3d0,18.9d0,52.0d0,24.2d0,35.0d0,1.3d0,8.0d0,0.45d0,0.04d0,0.5d0,
     &	2.22d0,0.11d0,0.067d0,0.10d0,50.0d0 /
	data (sji(j),j=1,30) / 7.97d7,1.23d7,2.0d4,0.078d0,-146.5d0,89.4d0,20.0d0,32.0d0,
     &	235.0d0,59.0d0,181.0d0,163.0d0,197.0d0,48.0d0,194.0d0,0.3d0,17.7d0,49.0d0,20.5d0,31.0d0,1.3d0,8.0d0,
     &	0.45d0,0.04d0,0.5d0,2.22d0,0.11d0,0.067d0,0.17d0,50.0d0 /
	external Graphite_x,Graphite_y

	r=1d-3*(Nc/468d0)**(1d0/3d0)

	call RegridDataLNK(Graphite_x,lam,e1x,e2x,nlam,.false.)
	call RegridDataLNK(Graphite_y,lam,e1y,e2y,nlam,.false.)

	do ilam=1,nlam
		if(computelam(ilam)) then
		call Q_MIE(e1x(ilam),e2x(ilam),lam(ilam),r,QEX,QSC,QAB,G)
		CabsGra=1d-8*(qab*pi*r**2)/3d0
		Csca(ilam)=1d-8*(qsc*pi*r**2)/3d0
		call Q_MIE(e1y(ilam),e2y(ilam),lam(ilam),r,QEX,QSC,QAB,G)
		CabsGra=CabsGra+2d0*1d-8*(qab*pi*r**2)/3d0
		Csca(ilam)=Csca(ilam)+2d0*1d-8*(qsc*pi*r**2)/3d0
		
		CabsGra=CabsGra/Nc
		Csca(ilam)=Csca(ilam)/NC
		
		x=1d0/lam(ilam)
		do j=1,nj
			S(j)=2d0*gj(j)*lj(j)*1d-4/(pi*((lam(ilam)/lj(j)-lj(j)/lam(ilam))**2+gj(j)**2))
			if(ionized) then
				if(j.eq.6) then
c use the 3.3 micron feature from Visser 2007
					S(j)=S(j)*sjn(j)*1d-20/(1d0+41d0/(Nc-14d0))
				else
					S(j)=S(j)*sji(j)*1d-20
				endif
			else
				S(j)=S(j)*sjn(j)*1d-20
			endif
		enddo
		S(6)=S(6)*HC
		do j=14,22
			S(j)=S(j)*HC
		enddo
		if(x.gt.17.25d0) then
			Cabs(ilam)=CabsGra
		else if(x.gt.15d0) then
			Cabs(ilam)=(126.0-6.4943*x)*1e-18
		else if(x.gt.10d0) then
			Cabs(ilam)=S(1)+(-3.0+1.35*x)*1e-18
		else if(x.gt.7.7d0) then
			Cabs(ilam)=(66.302-24.367*x+2.950*x**2-0.1057*x**3)*1e-18
		else if(x.gt.5.9d0) then
			Cabs(ilam)=S(2)+(1.8687+0.1905*x+0.4175*(x-5.9)**2+0.04370*(x-5.9)**3)*1e-18
		else if(x.gt.3.3d0) then
			Cabs(ilam)=S(2)+(1.8687+0.1905*x)*1e-18
		else
			Cabs(ilam)=34.58*10d0**(-18d0-3.431/x)*cutoffPAH(lam(ilam),Nc,ionized)
			do j=3,nj
				Cabs(ilam)=Cabs(ilam)+S(j)
			enddo
		endif
		a=(50d-4/r)**3
		if(a.gt.1d0) a=1d0
		qgra=0.01d0
		fpah=(1d0-qgra)*a
		Cabs(ilam)=(Cabs(ilam)*fpah+CabsGra*(1d0-fpah))/Mc
		if(Cabs(ilam).lt.0d0) Cabs(ilam)=0d0
		Csca(ilam)=Csca(ilam)/Mc
		endif
	enddo
			
	return
	end
	



	real*8 function cutoffPAH(lam,Nc,ionized)
	IMPLICIT NONE
	real*8 lam,Nc,y,M
	logical ionized
			
	if(Nc.gt.40d0) then
		M=0.4*Nc
	else
		M=0.3*Nc
	endif
	
	if(ionized) then
		y=1d0/(2.282*M**(-0.5)+0.889)
	else
		y=1d0/(3.804*M**(-0.5)+1.052)
	endif		
	y=y/lam
	
	cutoffPAH=atan(10d3*(y-1d0)**3/y)/3.1415926536+0.5d0
	
	return
	end
	

	
***********************************************************************
*       New Mie subroutine that approximates for big grains           *
*                                                                     *
***********************************************************************
      SUBROUTINE Q_MIE(E1,E2,LAM,RAD,QEX,QSC,QAB,G)
      IMPLICIT real*8 (A-H,O-Z)
      real*8 LAM,RAD,T,QEX,QSC,QAB,E1,E2,G,EV


C
C     MIE THEORY EFFICIENCY FACTORS FOR SPHERICAL PARTICLES OF
C     RADIUS 'RAD' AT WAVELENGTH 'LAM'.
C     E=E1 + I*E2 IS THE SQUARE OF THE COMPLEX REFRACTIVE INDEX.
C     THE REFRACTIVE INDEX IS GIVEN BY SUBROUTINE 'EPS'
C                         
      COMPLEX*16 E,RM,Y,ZN,ZN1,ZN2,C,A,B,AO,RRAT,A1,ANM1,BNM1
	complex*16,allocatable :: AN(:)

	T=0.0

      E=DCMPLX(E1,-E2)
      E=E**2.


      X=6.2831853*RAD/LAM
      IF(X.LT.0.001)THEN
C
C        USE SMALL PARTICLE FORMULAE.
C        CHANGED CRITERIION FROM X < 0.01 TO 0.001 BECAUSE SILICATE
C        SCATTERING WAS NOT CORRECT.
C	15-8-2001: Changed scattering formula from QSC=(X**4/.375)*DBLE(C**2)
C				into the correct formula QSC=(X**4/.375)*DABS(C)**2
C	Michiel Min
C
         C=(E-1.)/(E+2.)
         QSC=(X**4/.375)*cDABS(C)**2
         A=DIMAG(-4.*C)
         B=DIMAG(-C*(E*E+27.*E+38.)/(2.*E+3.)/3.75)
         QAB=X*(A+X*X*B)
         QEX=QAB+QSC
C
C        G THE ASYMMETRY PARAMETER IS ALWAYS NEGLIGIBLE FOR SMALL PARTICLES.
C
         G=0.0
         RETURN
      END IF
C
C     FULL MIE THEORY CALCULATION.
C     RM - COMPLEX REFRACTIVE INDEX
C
      RM=CDSQRT(E)
      EN1=DBLE(RM)
      EN2=DIMAG(RM)
      Y=X*RM
      ZN2=DCMPLX(DCOS(X),-DSIN(X))
      ZN1=DCMPLX(DSIN(X),DCOS(X))
      RIND=EN1**2+EN2**2     ! Rind = |rm|²
      NTIL=1.5*SQRT(RIND)*X+1

c	Number of iterations changed to improve for small |m| (Michiel Min)
	if(real(ntil).lt.(1.5*x))ntil=1.5*x

      NTOT=MAX0(20,NTIL)
c
      if (ntot.le.70000) then    ! go ahead with full Mie theory
	allocate(AN(NTOT))
c
      AN(NTOT)=DCMPLX(0,0)
      SUME=0.
      SUMS=0.
      SUMG1=0.
      SUMG2=0.
      PSG1=0.
      PSG2=0.
      NTOTA=NTOT
  100 P=DFLOAT(NTOTA)
      AN(NTOTA-1)=P/Y-(1./(P/Y+AN(NTOTA)))
      NTOTA=NTOTA-1
      IF(NTOTA.EQ.1) GOTO 101
      GOTO 100
  101 AO1=DSIN(EN1*X)*DCOS(EN1*X)
      EN2P=-EN2
c      IF(EN2P*X.GE.44.)WRITE(6,*)'EN2P,X,LAM,RAD,E1,E2',EN2P,X,LAM,
c     >RAD,E1,E2
      if(EN2P*X.GE.350.) then
         AO=dcmplx(0.0,1.0)
      else
        AO2=DSINH(EN2P*X)*DCOSH(EN2P*X)
        AO3=(DSIN(EN1*X))**2+(DSINH(EN2P*X))**2
        AO=DCMPLX(AO1,AO2)
        AO=AO/AO3
      endif
      A1=-1./Y+(1./(1./Y-AO))
      RRAT=A1/AN(1)
      f=2.0/(x*x)
      DO 4 N=1,NTOT
         AN(N)=AN(N)*RRAT
    4 CONTINUE 
      DO 2 N=1,NTOT
         P=DFLOAT(N)
         ZN=DFLOAT(2*N-1)*ZN1/X-ZN2
         C=AN(N)/RM+P/X
         A=C*DBLE(ZN)-DBLE(ZN1)
         A=A/(C*ZN-ZN1)
         C=RM*AN(N)+P/X
         B=C*DBLE(ZN)-DBLE(ZN1)
         B=B/(C*ZN-ZN1)
C
C        PP, PPG1, PPG2 ARE CONSTANTS CONTAINING THE N TERMS IN THE 
C        SUMMATIONS.
C
         PP=DFLOAT(2*N+1)
C
         PSS=PP*(A*dCONJG(A)+B*dCONJG(B))
         PSE=PP*DBLE(A+B)
         IF(N.GT.1)THEN
C
C           CALCULATE G USING FORMULA ON P.128 OF VAN DE HULST'S BOOK.
C           HAVE REPLACED N BY (N-1) IN THE FORMULA SO THAT WE CAN USE
C           PREVIOUS A(N) AND B(N) INSTEAD OF A(N+1) AND B(N+1)
C
            REN=DFLOAT(N)
            PPG1=(REN-1.)*(REN+1.)/REN
            PPG2=(2.*REN-1.)/((REN-1.)*REN)
            PSG1=PPG1*DBLE(ANM1*dCONJG(A)+BNM1*dCONJG(B))
            PSG2=PPG2*DBLE(ANM1*dCONJG(BNM1))
         END IF
         SUME=SUME+PSE
         SUMS=SUMS+PSS
         SUMG1=SUMG1+PSG1
         SUMG2=SUMG2+PSG2
         D1=ABS(PSE/SUME)
         D2=ABS(PSS/SUMS)
C        IF(D1.LT.1.E-7.AND.D2.LT.1.E-7) GO TO 5
         PT=ABS(PSS/PP)
         IF(PT.LE.1.E-20) GOTO 5
C
C        SAVE PREVIOUS A AND B FOR CALCULATION OF G THE ASYMMETRY PARAMETER
C
         ANM1=A
         BNM1=B
         ZN2=ZN1
         ZN1=ZN
    2 CONTINUE
    5 F=2.0/(X*X)
      QEX=F*SUME
      QSC=F*SUMS
      QAB=F*(SUME-SUMS)
      G=2.0*F*(SUMG1+SUMG2)/QSC
	deallocate(AN)
      RETURN
      else
c               Geometrical optics for big spheres
      call geopt(rm,ans)
      qex =2.0d0
      g=9.23d-01   !approx true for D&L silicate.......
      qsc=ans
      end if
      return
      END          
c******************************************************************************
      subroutine geopt(m,ans)
c      intgrates the reflection coefficient
c      trapezium rule integration from 0 to pi/2
      implicit real*8 (a-h,o-z)
      complex*16 m
      a=0.0d0
      b=1.570796327d0
      nstrip = 5000
      tot=0
      h=(b-a)/dfloat(nstrip)   !strip width
      tot=tot+0.5*ref(m,a)*h   !1st term
      do i=1,nstrip-1
       x=a+h*dfloat(i)
       tot=tot+ref(m,x)*h      !middle terms
      end do
      tot=tot+0.5*ref(m,b)*h   !last term
      ans=1.+2.*tot    !ans is Qsca
      return
      end
      
c******************************************************************************
                                                                               
      function ref(m,thetai)
c         Calculates Reflection coeffs
      implicit real*8 (a-h,o-z)
      complex*16 sinTHETAt,cosTHETAt ,m,rpll,rper          
      sinTHETAt=sin(THETAi)/m
      cosTHETAt=cdsqrt(1-(sinTHETAt*sinTHETAt))
c       r for E parallel to plane
      rpll = (cosTHETAt-m*cos(THETAi)) / (cosTHETAt+m*cos(THETAi))
c       r for E perp. to plane
      rper = (cos(THETAi)-m*cosTHETAt) / (cos(THETAi)+m*cosTHETAt)
C       R = ½(|rpll|²+|rper|²)
      R= (abs(rpll)*abs(rpll) + abs(rper)*abs(rper))/2.0
      ref=r*sin(THETAi)*cos(THETAi)
      return                                            
      end                 

C
C
      SUBROUTINE INTERP(X,Y,NPTS,NTERMS,XIN,YOUT)
      REAL*8 DELTAX,PROD,SUM,X(3000),Y(3000),XIN,YOUT
      REAL*8 DELTA(10),A(10)
      REAL*8 DENOM
**************************************************   
*     SEARCH FOR AN APPROPRIATE VALUE OF X(1)    *
**************************************************
   11 DO 19 I=1,NPTS
      IF (XIN-X(I)) 13,17,19
   13 I1=I-NTERMS/2
      IF(I1) 15,15,21
   15 I1=1
      GOTO 21
   17 YOUT=Y(I)
   18 GOTO 61
   19 CONTINUE
      I1=NPTS-NTERMS+1
   21 I2=I1+NTERMS-1
      IF (NPTS-I2) 23,31,31
   23 I2=NPTS
      I1=I2-NTERMS+1
   25 IF (I1) 26,26,31
   26 I1=1
   27 NTERMS=I2-I1+1
C
C  EVALUATE DEVIATIONS DELTA
C
   31 DENOM=X(I1+1)-X(I1)
      DELTAX=(XIN-X(I1))/DENOM
      DO 35 I=1,NTERMS
         IX=I1+I-1
   35 DELTA(I)=(X(IX)-X(I1)) / DENOM
**********************************************
*           ACCUMULATE COEFFICIENTS A        *
**********************************************
   40 A(1)=Y(I1)
   41 DO 50 K=2,NTERMS
         PROD=1.
         SUM=0.
         IMAX=K-1
         IXMAX=I1+IMAX
         DO 49 I=1,IMAX
            J=K-I
            PROD=PROD*(DELTA(K)-DELTA(J))
   49    SUM=SUM-A(J)/PROD
   50 A(K)=SUM+Y(IXMAX)/PROD
***********************************************
*         ACCUMULATE SUM OF EXPANSION         *
***********************************************
   51 SUM=A(1)
      DO 57 J=2,NTERMS
         PROD=1.
         IMAX=J-1
         DO 56 I=1,IMAX
   56    PROD=PROD*(DELTAX-DELTA(I))
   57 SUM=SUM+A(J)* PROD
   60 YOUT=SUM
   61 RETURN
      END                

