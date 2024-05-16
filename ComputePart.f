	subroutine ComputePart(C,ii,isize,computelamcloud)
	use GlobalSetup
	use Constants
	use AImod
	IMPLICIT NONE

	type(CloudType) C,Cdust
	integer MAXMAT
	parameter(MAXMAT=40)
	integer ii,isize

	real cext0,csca0,maxf
	real minlog,maxlog,pow,cabs0,totA
	real e1av,e2av,rad,r1,r2,tot,lmax,lmin,Mass,tot2,Ntot
	real lambda,Vol,rho_av
	real,allocatable :: r0(:),nr0(:,:),f(:),wf(:),rho(:)
	real,allocatable :: e1(:,:),e2(:,:)
	complex*16,allocatable :: m12(:)
	real*8 e1blend,e2blend
	real*8,allocatable :: frac(:)
	real,allocatable :: e1d(:,:),e2d(:,:)
	integer i,j,k,l,na,nf,ns,nm,ilam,Err,spheres,toolarge
	complex m,min,mav,alpha
	real QEXT, QSCA, QBS, GQSC,wvno,scale
	character*3 meth
	character*500 filename(MAXMAT),grid,tmp,tmp2,partfile,lnkfile

	real*8 rmie,lmie,e1mie,e2mie,csmie,cemie,KR,theta,dummy,amin,amax,rcore,camie,fmie
	logical truefalse,checkparticlefile,lnkloglog
	integer abun_in_name,LL,LLmax
	parameter(abun_in_name=2)
	real*8 Kabs(nlam+1),Ksca(nlam+1),Kext(nlam+1),lgrid(nlam+1)
	logical fcomputed,computelamcloud(nlam)
	real*8 csmie_fcomp,cemie_fcomp,gasdev

	write(meth,100)
100	format('DHS')

	na=180

	allocate(e1(MAXMAT,C%nlam))
	allocate(e2(MAXMAT,C%nlam))

	allocate(frac(MAXMAT))
	allocate(rho(MAXMAT))

	amin=C%rv(isize)
	amax=C%rv(isize)
	if((Cloud(ii)%type.eq."DIFFUSE".or.Cloud(ii)%type.eq."WATER").and.cloud_dens(isize,ii).lt.1d-40) then
		C%M(isize)=(3d0*4d0*pi*(amin*1d-4)**3)/3d0
		C%rho=3d0
		do ilam=1,C%nlam
			C%Kabs(isize,ilam)=1d0
			C%Ksca(isize,ilam)=1d0
			C%Kext(isize,ilam)=1d0
		enddo
		goto 301
	endif

	minlog=log10(amin)
	maxlog=log10(amax)
	pow=-3.5
	maxf=C%fmax
	
	lgrid(1:nlam)=lam(1:nlam)*1d4
	lgrid(nlam+1)=C%lref

	nm=0
	if(C%blend) then
		nm=maxval(C%nax(1:C%nmat))
		allocate(e1d(C%nmat+1,C%nlam))
		allocate(e2d(C%nmat+1,C%nlam))
		frac(1:C%nmat)=C%frac(isize,1:C%nmat)
		tot=0d0
		do i=1,C%nmat
			tot=tot+frac(i)
		enddo
		if(tot.gt.0d0) then
			frac=frac/tot
		else
			frac=1d0/real(C%nmat)
		endif
		frac(C%nmat+1)=C%porosity
		frac(1:C%nmat)=frac(1:C%nmat)*(1d0-C%porosity)
		e1d(C%nmat+1,1:C%nlam)=1d0
		e2d(C%nmat+1,1:C%nlam)=0d0
		rho(1:C%nmat)=C%rho_mat(1:C%nmat)
		rho(C%nmat+1)=0d0
		do i=1,nm
			do j=1,C%nmat
				e1d(j,1:C%nlam)=C%e1(j,minval((/i,C%nax(j)/)),1:C%nlam)
				e2d(j,1:C%nlam)=C%e2(j,minval((/i,C%nax(j)/)),1:C%nlam)
			enddo
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(j,e1blend,e2blend,m12)
!$OMP& SHARED(C,i,e1,e2,frac,e1d,e2d)
			allocate(m12(C%nmat+1))
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC,1)
			do j=1,C%nlam
				m12(1:C%nmat+1)=dcmplx(e1d(1:C%nmat+1,j),e2d(1:C%nmat+1,j))
				call Blender(m12(1:C%nmat+1),frac(1:C%nmat+1),C%nmat+1,e1blend,e2blend)
				e1(i,j)=e1blend
				e2(i,j)=e2blend
			enddo
!$OMP END DO
			deallocate(m12)
!$OMP FLUSH
!$OMP END PARALLEL
			rho_av=0d0
			do j=1,C%nmat+1
				rho_av=rho_av+frac(j)*rho(j)
			enddo
			rho(i)=rho_av
		enddo
		frac(1:nm)=1d0/real(nm)
		deallocate(e1d)
		deallocate(e2d)
	else
		do j=1,C%nmat
			do i=1,C%nax(j)
				if(C%frac(isize,j).gt.0d0) then
					nm=nm+1
					e1(nm,1:C%nlam)=C%e1(j,i,1:C%nlam)
					e2(nm,1:C%nlam)=C%e2(j,i,1:C%nlam)
					rho(nm)=C%rho_mat(j)
					frac(nm)=C%frac(isize,j)/(real(C%nax(j))*rho(nm))
				endif
			enddo
		enddo
		tot=0d0
		do i=1,nm
			tot=tot+frac(i)
		enddo
		if(tot.gt.0d0) then
			frac=frac/tot
		else
			frac=1d0/real(nm)
		endif
	endif
	min=dcmplx(1d0,0d0)

	if(maxf.eq.0d0) then
		nf=1
	else
		nf=20
	endif
	allocate(f(nf),wf(nf))
	if(nf.gt.1.and.maxf.gt.0.01e0) then
		call gauleg2(0.01e0,maxf,f(1:nf),wf(1:nf),nf)
	else if(maxf.eq.0e0) then
		f(1:nf)=0d0
		wf(1:nf)=1d0/real(nf)
	else
		f(1)=maxf
		wf(1)=1d0
	endif

	ns=1
	allocate(r0(ns),nr0(nm,ns))
	do l=1,nm
c		if(.false.) then
c			j=0
c			if(C%sigma(isize).le.1d-3.or.ns.eq.1) then
c				ns=1
c				r0(1)=C%rv(isize)
c				nr0(l,1)=1d0
c				tot=r0(1)**3
c			else
c101			tot=0d0
c			do k=1,ns
c				r0(k)=10d0**(minlog+(maxlog-minlog)*real(k-1)/real(ns-1))
c				nr0(l,k)=exp(-((r0(k)-C%rv(isize))/C%sigma(isize))**2)/(C%sigma(isize))
c				nr0(l,k)=nr0(l,k)*r0(k)
c				tot=tot+nr0(l,k)*r0(k)**3
c			enddo
c			if(.not.tot.gt.0d0) then
c				ns=1
c				r0(1)=C%rv(isize)
c				if(r0(1).lt.amin) r0(1)=amin
c				if(r0(1).gt.amax) r0(1)=amax
c				nr0(l,1)=1d0
c				tot=r0(1)**3
c			endif
c			endif
c			do k=1,ns
c				nr0(l,k)=frac(l)*nr0(l,k)/tot
c			enddo
c		else if(ns.eq.1) then
			r0(1)=10d0**((minlog+maxlog)/2d0)
			nr0(l,1)=frac(l)
			C%rv(isize)=r0(1)
c		else
c			tot=0d0
c			C%rv(isize)=0d0
c			do k=1,ns
c				r0(k)=10d0**(minlog
c     &				+(maxlog-minlog)*real(k-1)/real(ns-1))
c				nr0(l,k)=r0(k)**(pow+1d0)
c				tot=tot+nr0(l,k)*r0(k)**3
c			enddo
c			do k=1,ns
c				nr0(l,k)=frac(l)*nr0(l,k)/tot
c				C%rv(isize)=C%rv(isize)+nr0(l,k)*r0(k)**2
c			enddo
c			C%rv(isize)=sqrt(C%rv(isize))
c		endif
	enddo

	Mass=0d0
	Vol=0d0
	Ntot=0d0
	do l=1,nm
		if(frac(l).gt.0d0) then
			do k=1,ns
				r1=r0(k)
				do i=1,nf
					Mass=Mass+wf(i)*nr0(l,k)*rho(l)*4d0*pi*r1**3/3d0
					Vol=Vol+wf(i)*nr0(l,k)*4d0*pi*r1**3/3d0
					Ntot=Ntot+wf(i)*nr0(l,k)
				enddo
			enddo
		endif
	enddo
	C%M(isize)=Mass*(1d-4)**3/Ntot
	C%rho=Mass/Vol
	rho_av=Mass/Vol

	fmie=maxf
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,csca0,cabs0,cext0,theta,i,l,tot,k,Err,spheres,toolarge,
!$OMP&         rad,wvno,m,r1,rcore,qext,qsca,qbs,gqsc,rmie,lmie,e1mie,e2mie,
!$OMP&         csmie,cemie,camie,tot2,j,fcomputed)
!$OMP& SHARED(C,na,nm,ns,frac,minlog,maxlog,f,e1,e2,wf,isize,computelamcloud,Mass,
!$OMP&        pow,lgrid,rho,nf,r0,nr0,Kabs,Ksca,Kext,nlam,fmie,useDLMie,
!$OMP&		  DLMie_e1max,DLMie_e2max)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
	do ilam=1,C%nlam
	
	if(ilam.le.nlam) then
		if(.not.computelamcloud(ilam)) then
			Kabs(ilam)=1d-10
			Ksca(ilam)=1d-10
			Kext(ilam)=1d-10
			goto 11
		endif
	endif
	csca0=0d0
	cabs0=0d0
	cext0=0d0

	do l=1,nm
	if(frac(l).eq.0d0) goto 10
	do k=1,ns
	r1=r0(k)
	if(.not.r1.gt.0.001d0) r1=0.001d0
	if(.not.e1(l,ilam).gt.1d-6) e1(l,ilam)=1.0d-6
	if(e1(l,ilam).eq.1d0) e1(l,ilam)=1.001d0
	if(.not.e2(l,ilam).gt.1d-6) e2(l,ilam)=1.0d-6
	Err=0
	spheres=0
	toolarge=0
	fcomputed=.false.
	e1mie=e1(l,ilam)
	e2mie=e2(l,ilam)
	if(useDLMie.and.e1mie.lt.DLmie_e1max.and.e2mie.lt.DLmie_e2max) then
		rmie=r1
		lmie=lgrid(ilam)
		call DLMie(rmie,lmie,e1mie,e2mie,fmie,camie,csmie)
		cext0=cext0+nr0(l,k)*(camie+csmie)
		csca0=csca0+nr0(l,k)*csmie
	   	cabs0=cabs0+nr0(l,k)*camie
	else
	do i=1,nf
		rad=r1/(1d0-f(i))**(1d0/3d0)
		m=dcmplx(e1(l,ilam),-e2(l,ilam))
		wvno=2d0*3.1415926536/(lgrid(ilam))
		if(f(i).eq.0d0) then
			spheres=1
			goto 20
		endif
		if(r1*wvno.gt.1000d0) then
			toolarge=1
			goto 20
		endif
		rcore=rad*f(i)**(1d0/3d0)
		rmie=rad
		lmie=lgrid(ilam)
		e1mie=e1(l,ilam)
		e2mie=e2(l,ilam)
		if(rmie/lmie.lt.10d0) then
			call callBHCOAT(rmie,rcore,lmie,e1mie,e2mie,csmie,cemie,Err)
		else
			lmie=rmie/10d0
			call callBHCOAT(rmie,rcore,lmie,e1mie,e2mie,csmie,cemie,Err)
		endif
		if(.not.csmie.gt.0d0) then
			Err=1
		endif
20		if(Err.eq.1.or.spheres.eq.1.or.toolarge.eq.1) then
			rad=r1
			rcore=rad
			rmie=rad
			lmie=lgrid(ilam)
			e1mie=e1(l,ilam)
			e2mie=e2(l,ilam)
			if(Err.eq.1.or.i.eq.1) then
				if(rmie/lmie.lt.10d0) then
					call callBHMIE(rmie,lmie,e1mie,e2mie,csmie,cemie)
				else
					lmie=rmie/10d0
					call callBHMIE(rmie,lmie,e1mie,e2mie,csmie,cemie)
				endif
			endif
		endif
		cext0=cext0+wf(i)*nr0(l,k)*cemie
		csca0=csca0+wf(i)*nr0(l,k)*csmie
	   	cabs0=cabs0+wf(i)*nr0(l,k)*(cemie-csmie)
	enddo
	endif
	enddo
10	continue
	enddo

	Kabs(ilam)=1d4*cabs0/Mass
	Ksca(ilam)=1d4*csca0/Mass
	Kext(ilam)=1d4*cext0/Mass

	Kext(ilam)=Kabs(ilam)+Ksca(ilam)
	if(Kabs(ilam)/Kext(ilam).lt.1d-4) then
		Kabs(ilam)=Kext(ilam)*1d-4
		Kext(ilam)=Kabs(ilam)+Ksca(ilam)
	endif

11	continue

	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	C%Kabs(isize,1:C%nlam)=Kabs(1:C%nlam)
	C%Kext(isize,1:C%nlam)=Kext(1:C%nlam)
	C%Ksca(isize,1:C%nlam)=Ksca(1:C%nlam)
	
300	continue	

	deallocate(r0)
	deallocate(nr0)
	deallocate(f)
	deallocate(wf)

301	continue

	deallocate(e1)
	deallocate(e2)

	deallocate(frac)
	deallocate(rho)


	return
	end


	
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine readrefindCP(input,grid,e1,e2,n,loglog)
	IMPLICIT NONE
	real*8 grid(n)
	real*8 e1(n),e2(n),x0,y01,y02,x1,y11,y12,wp,gamma
	complex*16 m0,m1,m
	integer i,j,n
	character*500 input
	logical loglog

	open(unit=20,file=input,FORM="FORMATTED")
	i=1
1	read(20,*,end=102,err=1) x0,y01,y02
	if(y02.lt.1d-8) y02=1d-8
	wp=(1d0-y01)/x0**2
	gamma=y02/x0**3
103	if(x0.ge.grid(i)) then
		e1(i)=1d0-wp*grid(i)**2
		e2(i)=gamma*grid(i)**3
		e1(i)=y01
		e2(i)=y02
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) x1,y11,y12
	if(y12.lt.1d-8) y12=1d-8
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		e1(i)=10d0**(log10(y11)+(log10(grid(i)/x1))*(log10(y01/y11))/(log10(x0/x1)))
		e2(i)=10d0**(log10(y12)+(log10(grid(i)/x1))*(log10(y02/y12))/(log10(x0/x1)))
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y01=y11
	y02=y12
	goto 100
102	continue

	if(loglog) then
		m0=dcmplx(e1(i-1),e2(i-1))
		if(abs(m0).gt.2d0.and..false.) then
c don't use the conducting extrapolation since it is not very accurate
			do j=i,n
				m=m0*sqrt(grid(j)/grid(i-1))
				e1(j)=real(m)
				e2(j)=dimag(m)
			enddo
		else
c use loglog extrapolation
			m0=dcmplx(e1(i-2),e2(i-2))
			m1=dcmplx(e1(i-1),e2(i-1))
			do j=i,n
c				m=10d0**(log10(m0)+log10(m1/m0)*log10(grid(i-2)/grid(j))/log10(grid(i-2)/grid(i-1)))
c				e1(j)=real(m)
c				e2(j)=dimag(m)
				e1(j)=10d0**(log10(e1(i-2))+log10(e1(i-1)/e1(i-2))*log10(grid(i-2)/grid(j))/log10(grid(i-2)/grid(i-1)))
				e2(j)=10d0**(log10(e2(i-2))+log10(e2(i-1)/e2(i-2))*log10(grid(i-2)/grid(j))/log10(grid(i-2)/grid(i-1)))
			enddo
		endif
	else
c use the dielectric extrapolation, this is the default
		do j=i,n
			e1(j)=e1(i-1)
			e2(j)=e2(i-1)*grid(i-1)/grid(j)
		enddo
	endif

	close(unit=20)
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine Blender(m,abun,nm,e1out,e2out)
	IMPLICIT NONE
	integer,intent(in) :: nm
	real*8,intent(in) :: abun(nm)
	complex*16,intent(in) :: m(nm)
	real*8,intent(out) :: e1out,e2out
	integer j,iter
	complex*16,save :: mm,me,sum,m2(100)
	logical,save :: doit(100)
!$OMP THREADPRIVATE(mm,me,sum,m2,doit)


c LLL mixing rule (not preferred)
	mm=0d0
	do j=1,nm
		doit(j)=(abun(j).gt.0d0)
		if(doit(j)) then
			m2(j)=m(j)**2
			mm=mm+m(j)**(2d0/3d0)*abun(j)
		endif
	enddo
	mm=mm**3d0

	do iter=1,100
		sum=0d0
		do j=1,nm
			if(doit(j)) sum=sum+((m2(j)-mm)/(m2(j)+2d0*mm))*abun(j)
		enddo
		me=(2d0*sum+1d0)/(1d0-sum)
		me=mm*me
		mm=me
		if(cdabs(sum).lt.1d-8) exit
	enddo

	e1out=dreal(cdsqrt(me))
	e2out=dimag(cdsqrt(me))

	return
	end
	


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE gauleg2(x1,x2,x,w,n)
      INTEGER n
      REAL x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END


	subroutine callBHMIE(rmie,lmie,e1mie,e2mie,csmie,cemie)
	IMPLICIT NONE
	real*8 rmie,lmie,e1mie,e2mie,csmie,cemie
	real*8 pi
	parameter(pi=3.14159265358979323846264338328d0)
      INTEGER NANG,Err
      REAL GSCA,QBACK,QEXT,QSCA,X
      COMPLEX REFREL
      COMPLEX S1(2001),S2(2001)
	NANG=2
	X=2d0*pi*rmie/lmie
	REFREL=cmplx(e1mie,e2mie)
	call BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA,Err)
	if(Err.gt.0) then
1		x=x/2d0
		call BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA,Err)
		if(Err.gt.0d0) goto 1
	endif
	csmie=pi*rmie**2*QSCA
	cemie=pi*rmie**2*QEXT
	
	return
	end

	subroutine callBHCOAT(rmie,rcore,lmie,e1mie,e2mie,csmie,cemie,Err)
	IMPLICIT NONE
	real*8 rmie,lmie,e1mie,e2mie,csmie,cemie,rcore
	real*8 pi,theta,test
	integer Err,i
	parameter(pi=3.14159265358979323846264338328d0)
      REAL GSCA,QBACK,QEXT,QSCA,X,Y, WVNO
      COMPLEX REFRL1,REFRL2

      INTEGER   MAXANG, NUMANG
      INTEGER	MAXANG1,NUMANG1
      parameter(MAXANG1=20,NUMANG1=2)
      REAL      MU( NUMANG1 ), D21( MAXANG1, 2 ), M1( MAXANG1, 2 ),
     &          M2( MAXANG1, 2 ), S21( MAXANG1, 2 )


	X=2d0*pi*rcore/lmie
	Y=2d0*pi*rmie/lmie
	test=abs(cmplx(e1mie,e2mie))
	if(test.lt.3d0.and.test.gt.0.5d0.and.e2mie.lt.2d0) then
		REFRL1=cmplx(1.001,1e-8)
		REFRL2=cmplx(e1mie,e2mie)
		Err=0
		call BHCOAT(X,Y,REFRL1,REFRL2,QEXT,QSCA,QBACK,GSCA)
		if(.not.QEXT.gt.0e0) Err=1
	else
		Err=1
	endif
	if(Err.eq.1) then
		Err=0
		wvno=1d0
		REFRL1=cmplx(1.001,-1e-8)
		REFRL2=cmplx(e1mie,-e2mie)
		NUMANG=NUMANG1
		MAXANG=MAXANG1
		do i=1,NUMANG
			theta=(real(i)-0.5)/real(NUMANG)*3.1415926536/2d0
			mu(i)=cos(theta)
		enddo

		call DMiLay( X, Y, WVNO, REFRL2, REFRL1, MU,
     &                   NUMANG, QEXT, QSCA, QBACK, GSCA, 
     &                   M1, M2, S21, D21, MAXANG, Err)
	endif

	csmie=pi*rmie**2*QSCA
	cemie=pi*rmie**2*QEXT
     	
	return
	end


      SUBROUTINE BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA,Err)
      IMPLICIT NONE

C Declare parameters:
C Note: important that MXNANG be consistent with dimension of S1 and S2
C       in calling routine!

      INTEGER MXNANG,NMXX
C      PARAMETER(MXNANG=1000,NMXX=15000)
      PARAMETER(MXNANG=1000,NMXX=15000)

C Arguments:

      INTEGER NANG,Err
      REAL GSCA,QBACK,QEXT,QSCA,X
      COMPLEX REFREL
      COMPLEX S1(2*MXNANG-1),S2(2*MXNANG-1)

C Local variables:

      LOGICAL SINGLE
      INTEGER J,JJ,N,NSTOP,NMX,NN
      DOUBLE PRECISION CHI,CHI0,CHI1,DANG,DX,EN,FN,P,PII,PSI,PSI0,PSI1,
     &                 THETA,XSTOP,YMOD
      DOUBLE PRECISION
     &   AMU(MXNANG),
     &   PI(MXNANG),
     &   PI0(MXNANG),
     &   PI1(MXNANG),
     &   TAU(MXNANG)
      DOUBLE COMPLEX
     &   DCXS1(2*MXNANG-1),
     &   DCXS2(2*MXNANG-1)

C***********************************************************************
C
C Subroutine BHMIE is derived from the Bohren-Huffman Mie scattering
C     subroutine to calculate scattering and absorption by a homogenous
C     isotropic sphere.
C Given:
C    X = 2*pi*a/lambda
C    REFREL = (complex refr. index of sphere)/(real index of medium)
C    NANG = number of angles between 0 and 90 degrees
C           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
C           if called with NANG<2, will set NANG=2 and will compute
C           scattering for theta=0,90,180.
C Returns:
C    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
C                                scatt. E perp. to scatt. plane)
C    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
C                                scatt. E parr. to scatt. plane)
C    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
C    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
C    QBACK = 4.*pi*(dC_sca/domega)/pi*a**2
C          = backscattering efficiency
C    GSCA = <cos(theta)> for scattering
C
C S1 and S2 are the diagonal elements of the "amplitude scattering matrix"
C (see eq. 3.12 of Bohren & Huffman 1983) -- the off-diagonal elements
C vanish for a spherical target.
C For unpolarized incident light, the intensity of scattered light a
C distance r from the sphere is just
C          1
C  I_s = ------ * I_in * S_11
C        (kr)^2
C
C where k=2*pi/lambda 
C and the "Muller matrix element" S_11 = 0.5*( |S_1|^2 + |S_2|^2 )
C
C for incident light polarized perp to the scattering plane,
C the scattered light is polarized perp to the scattering plane
C with intensity I_s = I_in * |S_1|^2 / (kr)^2
C
C for incident light polarized parallel to the scattering plane,
C the scattered light is polarized parallel to the scattering plane
C with intensity I_s = I_in * |S_2|^2 / (kr)^2
C
C History:
C Original program taken from Bohren and Huffman (1983), Appendix A
C Modified by B.T.Draine, Princeton Univ. Obs., 90.10.26
C in order to compute <cos(theta)>
C 91.05.07 (BTD): Modified to allow NANG=1
C 91.08.15 (BTD): Corrected error (failure to initialize P)
C 91.08.15 (BTD): Modified to enhance vectorizability.
C 91.08.15 (BTD): Modified to make NANG=2 if called with NANG=1
C 91.08.15 (BTD): Changed definition of QBACK.
C 92.01.08 (BTD): Converted to full double precision and double complex
C                 eliminated 2 unneed lines of code
C                 eliminated redundant variables (e.g. APSI,APSI0)
C                 renamed RN -> EN = double precision N
C                 Note that DOUBLE COMPLEX and DCMPLX are not part
C                 of f77 standard, so this version may not be fully
C                 portable.  In event that portable version is
C                 needed, use src/bhmie_f77.f
C 93.06.01 (BTD): Changed AMAX1 to generic function MAX
C 98.09.17 (BTD): Added variable "SINGLE" and warning in event that
C                 code is used with single-precision arithmetic (i.e.,
C                 compiler does not support DOUBLE COMPLEX)
C 99.02.17 (BTD): Replaced calls to REAL() and IMAG() by
C                 REALPART() and IMAGPART() for compatibility with g77
C                 Note that when code is used with standard f77 
C                 compilers, it is now necessary to enable two lines
C                 defining functions REALPART(X) and IMAGPART(X)
C 99.02.19 (BTD): added lines to be enabled to properly define
C                 REALPART() and IMAGPART() if NOT using g77
C                 ***see below!!***
C 01.02.16 (BTD): added IMPLICIT NONE
C 01.02.27 (BTD): changed definition of QBACK back to convention of
C                 Bohren & Huffman and others:
C                 Q_back = 4.*pi*(dC_sca/dOmega)/(pi*a^2) in backward
C                          direction
c 02.03.09 (BTD): defined statement function REALPART_SP to
c                 avoid warning regarding type conversion when taking
c                 real part of S1(1) to evaluate QEXT
c                 some cleanup regarding type conversion
c 02.05.30 (BTD): introduced internal double complex arrays DCXS1,DCXS2
c                 to possibly increase accuracy during summations.
c                 After summations, output scattering amplitudes
c                 via single complex arrays S1,S2 as before.
c                 Usage of this routine is unaffected by change.
c                 Note: no longer need statement function REALPART_SP
c 02.09.18 (BTD): Error in evaluation of QBACK reported by Okada Yasuhiko
c                 Was calculating QBACK using S1 rather than DCXS1
c                 Corrected.
c 02.10.16 (BTD): Added comments explaining definition of S_1 and S_2 .
C end history
C
C***********************************************************************
C 
C This module is dependent on whether compiler supports double precision
C complex variables:
C
C If your compiler does NOT support double complex, comment out following
C three lines, and uncomment corresponding 3 lines further below
C
      DOUBLE COMPLEX AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
      DOUBLE COMPLEX D(NMXX)
      PARAMETER(SINGLE=.FALSE.)

C      COMPLEX AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
C      COMPLEX D(NMXX)
C      PARAMETER(SINGLE=.TRUE.)

C**********************************************************************

C Following five statements should be enabled if NOT using g77.
C They assume that the compiler supports double complex, since the
C statements DBLE and DIMAG are used.  If double complex is not available
C (see above) you will need to change DIMAG to AIMAG
C
C If using g77, following statements could be commented out, as 
C REALPART and IMAGPART are g77 intrinsic functions
C However, they do not need to be commented out.

      DOUBLE COMPLEX DPCX
      DOUBLE PRECISION REALPART
      DOUBLE PRECISION IMAGPART
      REALPART(DPCX)=(DBLE(DPCX))
      IMAGPART(DPCX)=(DIMAG(DPCX))

	Err=0      
C***********************************************************************
C*** Safety checks

      IF(SINGLE)WRITE(0,*)'Warning: this version of bhmie uses only ',
     &          'single precision complex numbers!'
      IF(NANG.GT.MXNANG)NANG=MXNANG
      IF(NANG.LT.2)NANG=2

C*** Obtain pi:

      PII=4.D0*ATAN(1.D0)
      DX=X
      DREFRL=REFREL
      Y=X*DREFRL
      YMOD=ABS(Y)

C*** Series expansion terminated after NSTOP terms
C    Logarithmic derivatives calculated from NMX on down
      XSTOP=X+4.*X**0.3333+2.
      NMX=NINT(MAX(XSTOP,YMOD))+15
C BTD experiment 91.1.15: add one more term to series and compare results
C      NMX=MAX(XSTOP,YMOD)+16
C test: compute 7001 wavelengths between .0001 and 1000 micron
C for a=1.0micron SiC grain.  When NMX increased by 1, only a single
C computed number changed (out of 4*7001) and it only changed by 1/8387
C conclusion: we are indeed retaining enough terms in series!

      NSTOP=NINT(XSTOP)

      IF(NMX.GT.NMXX.or..not.NMX.gt.0)THEN
c         WRITE(0,*)'Error: NMX > NMXX=',NMXX,' for |m|x=',YMOD
         Err=1
         return
      ENDIF

C*** Require NANG.GE.1 in order to calculate scattering intensities

      DANG=0.
      IF(NANG.GT.1)DANG=.5*PII/DBLE(NANG-1)
      DO J=1,NANG
         THETA=DBLE(J-1)*DANG
         AMU(J)=COS(THETA)
      ENDDO
      DO J=1,NANG
         PI0(J)=0.
         PI1(J)=1.
      ENDDO
      NN=2*NANG-1
      DO J=1,NN
         DCXS1(J)=(0.D0,0.D0)
         DCXS2(J)=(0.D0,0.D0)
      ENDDO

C*** Logarithmic derivative D(J) calculated by downward recurrence
C    beginning with initial value (0.,0.) at J=NMX

      D(NMX)=(0.,0.)
      NN=NMX-1
      DO N=1,NN
         EN=NMX-N+1
         D(NMX-N)=(EN/Y)-(1./(D(NMX-N+1)+EN/Y))
      ENDDO

C*** Riccati-Bessel functions with real argument X
C    calculated by upward recurrence

      PSI0=COS(DX)
      PSI1=SIN(DX)
      CHI0=-SIN(DX)
      CHI1=COS(DX)
      XI1=DCMPLX(PSI1,-CHI1)
      QSCA=0.E0
      GSCA=0.E0
      P=-1.
      DO N=1,NSTOP
         EN=N
         FN=(2.E0*EN+1.)/(EN*(EN+1.))

C for given N, PSI  = psi_n        CHI  = chi_n
C              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
C              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
C Calculate psi_n and chi_n

         PSI=(2.E0*EN-1.)*PSI1/DX-PSI0
         CHI=(2.E0*EN-1.)*CHI1/DX-CHI0
         XI=DCMPLX(PSI,-CHI)

C*** Store previous values of AN and BN for use
C    in computation of g=<cos(theta)>

         IF(N.GT.1)THEN
            AN1=AN
            BN1=BN
         ENDIF

C*** Compute AN and BN:

         AN=(D(N)/DREFRL+EN/DX)*PSI-PSI1
         AN=AN/((D(N)/DREFRL+EN/DX)*XI-XI1)
         BN=(DREFRL*D(N)+EN/DX)*PSI-PSI1
         BN=BN/((DREFRL*D(N)+EN/DX)*XI-XI1)

C*** Augment sums for Qsca and g=<cos(theta)>

         QSCA=QSCA+REAL((2.*EN+1.)*(ABS(AN)**2+ABS(BN)**2))
         GSCA=GSCA+REAL(((2.*EN+1.)/(EN*(EN+1.)))*
     &        (REALPART(AN)*REALPART(BN)+IMAGPART(AN)*IMAGPART(BN)))
         IF(N.GT.1)THEN
            GSCA=GSCA+REAL(((EN-1.)*(EN+1.)/EN)*
     &      (REALPART(AN1)*REALPART(AN)+IMAGPART(AN1)*IMAGPART(AN)+
     &      REALPART(BN1)*REALPART(BN)+IMAGPART(BN1)*IMAGPART(BN)))
         ENDIF

C*** Now calculate scattering intensity pattern
C    First do angles from 0 to 90

         DO J=1,NANG
            JJ=2*NANG-J
            PI(J)=PI1(J)
            TAU(J)=EN*AMU(J)*PI(J)-(EN+1.)*PI0(J)
            DCXS1(J)=DCXS1(J)+FN*(AN*PI(J)+BN*TAU(J))
            DCXS2(J)=DCXS2(J)+FN*(AN*TAU(J)+BN*PI(J))
         ENDDO

C*** Now do angles greater than 90 using PI and TAU from
C    angles less than 90.
C    P=1 for N=1,3,...; P=-1 for N=2,4,...

         P=-P
         DO J=1,NANG-1
            JJ=2*NANG-J
            DCXS1(JJ)=DCXS1(JJ)+FN*P*(AN*PI(J)-BN*TAU(J))
            DCXS2(JJ)=DCXS2(JJ)+FN*P*(BN*PI(J)-AN*TAU(J))
         ENDDO
         PSI0=PSI1
         PSI1=PSI
         CHI0=CHI1
         CHI1=CHI
         XI1=DCMPLX(PSI1,-CHI1)

C*** Compute pi_n for next value of n
C    For each angle J, compute pi_n+1
C    from PI = pi_n , PI0 = pi_n-1

         DO J=1,NANG
            PI1(J)=((2.*EN+1.)*AMU(J)*PI(J)-(EN+1.)*PI0(J))/EN
            PI0(J)=PI(J)
         ENDDO
      ENDDO

C*** Have summed sufficient terms.
C    Now compute QSCA,QEXT,QBACK,and GSCA

      GSCA=REAL(2.D0*GSCA/QSCA)
      QSCA=REAL((2.D0/(DX*DX))*QSCA)
      QEXT=REAL((4.D0/(DX*DX))*REALPART(DCXS1(1)))
      QBACK=REAL(4.D0*(ABS(DCXS1(2*NANG-1))/DX)**2)

C prepare single precision complex scattering amplitude for output

      DO J=1,2*NANG-1
         S1(J)=CMPLX(DCXS1(J))
         S2(J)=CMPLX(DCXS2(J))
      ENDDO

      RETURN
      END



      SUBROUTINE BHCOAT(XX,YY,RRFRL1,RRFRL2,QQEXT,QQSCA,QBACK,GSCA)
      IMPLICIT NONE

! Arguments:

      REAL GSCA,QBACK,QQEXT,QQSCA,XX,YY
      COMPLEX RRFRL1,RRFRL2

! Local variables:

      DOUBLE COMPLEX II
      PARAMETER(II=(0.0d0,1.D0))
      DOUBLE PRECISION DEL
      PARAMETER(DEL=1.D-8)

      INTEGER IFLAG,N,NSTOP
      DOUBLE PRECISION
     &   CHI0Y,CHI1Y,CHIY,EN,PSI0Y,PSI1Y,PSIY,QEXT,QSCA,RN,X,Y,YSTOP
      DOUBLE COMPLEX
     &   AMESS1,AMESS2,AMESS3,AMESS4,AN,AN1,ANCAP,
     &   BN,BN1,BNCAP,BRACK,
     &   CHI0X2,CHI0Y2,CHI1X2,CHI1Y2,CHIX2,CHIPX2,CHIPY2,CHIY2,CRACK,
     &   D0X1,D0X2,D0Y2,D1X1,D1X2,D1Y2,DNBAR,GNBAR,
     &   REFREL,RFREL1,RFREL2,
     &   XBACK,XI0Y,XI1Y,XIY,
     &   X1,X2,Y2

!***********************************************************************
!
! Subroutine BHCOAT calculates Q_ext, Q_sca, Q_back, g=<cos> 
! for coated sphere.
! All bessel functions computed by upward recurrence.
! Input:
!        X = 2*PI*RCORE*REFMED/WAVEL
!        Y = 2*PI*RMANT*REFMED/WAVEL
!        RFREL1 = REFCOR/REFMED
!        RFREL2 = REFMAN/REFMED 
! where  REFCOR = complex refr.index of core)
!        REFMAN = complex refr.index of mantle)
!        REFMED = real refr.index of medium)
!        RCORE = radius of core
!        RMANT = radius of mantle
!        WAVEL = wavelength of light in ambient medium

! returns:
!        QQEXT = C_ext/pi*rmant^2
!        QQSCA = C_sca/pi*rmant^2
!        QBACK = 4*pi*(dQ_sca/dOmega)
!              = "radar backscattering efficiency factor"
!        GSCA  = <cos(theta)> for scattered power
!
! Routine BHCOAT is taken from Bohren & Huffman (1983)
! extended by Prof. Francis S. Binkowski of The University of North
! Carolina at Chapel Hill to evaluate GSCA=<cos(theta)>
! History:
! 92.11.24 (BTD) Explicit declaration of all variables
! 00.05.16 (BTD) Added IMPLICIT NONE
! 12.04.10 (FSB) Modified by Prof. Francis S. Binkowski of
!                The University of North Carolina at Chapel Hill
!                to evaluate GSCA=<cos(theta)>
! 12.06.15 (BTD) Cosmetic changes
!***********************************************************************

      X=XX
      Y=YY
      RFREL1=RRFRL1
      RFREL2=RRFRL2
!         -----------------------------------------------------------
!              del is the inner sphere convergence criterion
!         -----------------------------------------------------------
      X1=RFREL1*X
      X2=RFREL2*X
      Y2=RFREL2*Y
      YSTOP=Y+4.*Y**0.3333+2.0
      REFREL=RFREL2/RFREL1
      NSTOP=YSTOP
!         -----------------------------------------------------------
!              series terminated after nstop terms
!         -----------------------------------------------------------
      D0X1=COS(X1)/SIN(X1)
      D0X2=COS(X2)/SIN(X2)
      D0Y2=COS(Y2)/SIN(Y2)
      PSI0Y=COS(Y)
      PSI1Y=SIN(Y)
      CHI0Y=-SIN(Y)
      CHI1Y=COS(Y)
      XI0Y=PSI0Y-II*CHI0Y
      XI1Y=PSI1Y-II*CHI1Y
      CHI0Y2=-SIN(Y2)
      CHI1Y2=COS(Y2)
      CHI0X2=-SIN(X2)
      CHI1X2=COS(X2)
      QSCA=0.0
      QEXT=0.0
      XBACK=(0.0,0.0)
      IFLAG=0
      DO N=1,NSTOP
         RN=N
         EN=RN
         PSIY=(2.0*RN-1.)*PSI1Y/Y-PSI0Y
         CHIY=(2.0*RN-1.)*CHI1Y/Y-CHI0Y
         XIY=PSIY-II*CHIY
         D1Y2=1.0/(RN/Y2-D0Y2)-RN/Y2
         IF(IFLAG.EQ.0)THEN

! calculate inner sphere ancap, bncap
!           and brack and crack

            D1X1=1.0/(RN/X1-D0X1)-RN/X1
            D1X2=1.0/(RN/X2-D0X2)-RN/X2
            CHIX2=(2.0*RN-1.0)*CHI1X2/X2-CHI0X2
            CHIY2=(2.0*RN-1.0)*CHI1Y2/Y2-CHI0Y2
            CHIPX2=CHI1X2-RN*CHIX2/X2
            CHIPY2=CHI1Y2-RN*CHIY2/Y2
            ANCAP=REFREL*D1X1-D1X2
            ANCAP=ANCAP/(REFREL*D1X1*CHIX2-CHIPX2)
            ANCAP=ANCAP/(CHIX2*D1X2-CHIPX2)
            BRACK=ANCAP*(CHIY2*D1Y2-CHIPY2)
            BNCAP=REFREL*D1X2-D1X1
            BNCAP=BNCAP/(REFREL*CHIPX2-D1X1*CHIX2)
            BNCAP=BNCAP/(CHIX2*D1X2-CHIPX2)
            CRACK=BNCAP*(CHIY2*D1Y2-CHIPY2)

! calculate convergence test expressions for inner sphere
! see pp 483-485 of Bohren & Huffman for definitions

            AMESS1=BRACK*CHIPY2
            AMESS2=BRACK*CHIY2
            AMESS3=CRACK*CHIPY2
            AMESS4=CRACK*CHIY2

         ENDIF ! test on iflag.eq.0

! now test for convergence for inner sphere
! all four criteria must be satisfied
! see p 484 of Bohren & Huffman

         IF(ABS(AMESS1).LT.DEL*ABS(D1Y2).AND.
     &      ABS(AMESS2).LT.DEL.AND.
     &      ABS(AMESS3).LT.DEL*ABS(D1Y2).AND.
     &      ABS(AMESS4).LT.DEL)THEN

! convergence for inner sphere

            BRACK=(0.,0.)
            CRACK=(0.,0.)
            IFLAG=1
         ELSE

! no convergence yet

            IFLAG=0

         ENDIF
         DNBAR=D1Y2-BRACK*CHIPY2
         DNBAR=DNBAR/(1.0-BRACK*CHIY2)
         GNBAR=D1Y2-CRACK*CHIPY2
         GNBAR=GNBAR/(1.0-CRACK*CHIY2)

! store previous values of an and bn for use in computation of 
! g=<cos(theta)>

         IF(N.GT.1)THEN
            AN1=AN
            BN1=BN
         ENDIF

! update an and bn

         AN=(DNBAR/RFREL2+RN/Y)*PSIY-PSI1Y
         AN=AN/((DNBAR/RFREL2+RN/Y)*XIY-XI1Y)
         BN=(RFREL2*GNBAR+RN/Y)*PSIY-PSI1Y
         BN=BN/((RFREL2*GNBAR+RN/Y)*XIY-XI1Y)

! calculate sums for qsca,qext,xback

         QSCA=QSCA+(2.0*RN+1.0)*(ABS(AN)*ABS(AN)+ABS(BN)*ABS(BN))
         XBACK=XBACK+(2.0*RN+1.0)*(-1.)**N*(AN-BN)
         QEXT=QEXT+(2.0*RN+1.0)*(DBLE(AN)+DBLE(BN))

! (FSB) calculate the sum for the asymmetry factor

         GSCA=GSCA+((2.0*EN+1.)/(EN*(EN+1.0)))*
     &        (REAL(AN)*REAL(BN)+IMAG(AN)*IMAG(BN))
         IF(N.GT.1)THEN
            GSCA=GSCA+((EN-1.)*(EN+1.)/EN)*
     &           (REAL(AN1)*REAL(AN)+IMAG(AN1)*IMAG(AN)+
     &            REAL(BN1)*REAL(BN)+IMAG(BN1)*IMAG(BN))
         ENDIF

! continue update for next iteration

         PSI0Y=PSI1Y
         PSI1Y=PSIY
         CHI0Y=CHI1Y
         CHI1Y=CHIY
         XI1Y=PSI1Y-II*CHI1Y
         CHI0X2=CHI1X2
         CHI1X2=CHIX2
         CHI0Y2=CHI1Y2
         CHI1Y2=CHIY2
         D0X1=D1X1
         D0X2=D1X2
         D0Y2=D1Y2
      ENDDO

! have summed sufficient terms
! now compute QQSCA,QQEXT,QBACK, and GSCA
      
      QQSCA=(2.0/(Y*Y))*QSCA
      QQEXT=(2.0/(Y*Y))*QEXT
      QBACK=(ABS(XBACK))**2
      QBACK=(1.0/(Y*Y))*QBACK
      GSCA=2.0*GSCA/QSCA
      RETURN
      END



c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


      SUBROUTINE DMiLay( RCORE, RSHELL, WVNO, RINDSH, RINDCO, MU,
     &                   NUMANG, QEXT, QSCA, QBS, GQSC, 
     &                   M1, M2, S21, D21, MAXANG, Err)

c **********************************************************************
c    DOUBLE PRECISION version of MieLay, which computes electromagnetic 
c    scattering by a stratified sphere, i.e. a particle consisting of a 
c    spherical core surrounded by a spherical shell.  The surrounding 
c    medium is assumed to have refractive index unity.  The formulas, 
c    manipulated to avoid the ill-conditioning that plagued earlier 
c    formulations, were published in:

c        Toon, O. and T. Ackerman, Applied Optics 20, 3657 (1981)

c    Further documentation, including definitons of input and output
c    arguments, is inside the single precision version of this program
c    (SUBROUTINE MieLay, available by anonymous ftp from 
c    climate.gsfc.nasa.gov in directory pub/wiscombe).

c    It is recommended to use this DOUBLE PRECISION version for IEEE 
c    arithmetic (32-bit floating point) computers, just to be safe.
c    If computer time is critical, back-to-back tests with the single
c    precision version should be done for ranges of radii and refractive
c    index relevant to your particular problem, before adopting the
c    single precision version.  This version is also recommended for
c    cases of large size parameter (bigger than 10 or so) and/or large 
c    imaginary refractive index (bigger than 1 or so) and also whenever 
c    overflows or strange behavior are encountered in running the
c    single precision version.  Sometimes the bigger exponent range in
c    DOUBLE PRECISION is as important as the added precision.

c    This version is designed to be interchangeable with the single
c    precision version:  all the floating-point input arguments are
c    still single precision.  Only the name of the routine has been
c    changed to prevent confusion (and it is strongly urged not to
c    change it to MieLay for the same reason).

c **********************************************************************

c     .. Parameters ..

      INTEGER   MXANG, LL, Err
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER ( MXANG = 200, LL = 200000, ZERO = 0.0D0, ONE = 1.0D0,
     &            TWO = 2.0D0 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   MAXANG, NUMANG
      REAL      GQSC, QBS, QEXT, QSCA, RCORE, RSHELL, WVNO
      COMPLEX   RINDCO, RINDSH
c     ..
c     .. Array Arguments ..

      REAL      MU( NUMANG ), D21( MAXANG, 2 ), M1( MAXANG, 2 ),
     &          M2( MAXANG, 2 ), S21( MAXANG, 2 )
c     ..
c     .. Local Scalars ..

      LOGICAL   INPERR, PASS1
      INTEGER   J, K, M, N, NMX1, NMX2, NN

      DOUBLE PRECISION  AA, AIM, AM1IM, AM1RE, ARE, BB, BIM, BM1IM,
     &                  BM1RE, BRE, CC, COSX1, COSX4, DD, DENOM, 
     &                  DGQSC, DQEXT, DQSCA, E2Y1,
     &                  EY1, EY1MY4, EY1PY4, EY4, FOURPI, PINUM,
     &                  RMM, RX, SINX1, SINX4, TOLER, X1, X4,
     &                  XCORE, XSHELL, Y1, Y4

      DOUBLE COMPLEX  AC, ACOE, ACOEM1, BC, BCOE, BCOEM1, CI, CZERO,
     &                DH1, DH2, DH4, DUMMY, DUMSQ, K1, K2, K3, 
     &                P24H21, P24H24, RRFX, SBACK, WM1
c     ..
c     .. Local Arrays ..

      DOUBLE PRECISION  PI( MXANG, 3 ), SI2THT( MXANG ), T( 5 ),
     &                  TA( 4 ), TAU( MXANG, 3 )

      DOUBLE COMPLEX  S1( MXANG, 2 ), S2( MXANG, 2 ),
     &                U( 8 ), WFN( 2 ), Z( 4 )
	double complex,allocatable :: W(:,:),acap(:)
c     ..
c     .. External Functions ..

      LOGICAL   WRTBAD, WRTDIM
      EXTERNAL  WRTBAD, WRTDIM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, DIMAG, ASIN, DCMPLX, COS, EXP, MOD, DBLE, SIN
c     ..
c     .. Save statement ..

      SAVE  PINUM, PASS1
c     ..
c     .. Data statements ..

      DATA      PASS1 / .True. / , TOLER / 1.D-6 / ,
     &          CZERO / ( 0.D0, 0.D0 ) / , CI / ( 0.D0, 1.D0 ) /
c     ..

c	allocate(w(3,LL))
c	allocate(acap(LL))

c      IF( PASS1 ) THEN

         PINUM  = TWO*ASIN( ONE )
         PASS1  = .False.

c      END IF

      XSHELL = RSHELL*WVNO
      XCORE  = RCORE*WVNO
      T( 1 ) = XSHELL*ABS( RINDSH )
      NMX1   = 1.1D0*T( 1 )
      NMX2   = T( 1 )

      IF( NMX1.LE.150 ) THEN

         NMX1   = 150
         NMX2   = 135

      END IF

c                        ** Check input arguments for gross errors
      INPERR = .False.

      IF( WVNO.LE.0.0 ) INPERR = WRTBAD( 'WVNO' )

      IF( RSHELL.LE.0.0 ) INPERR = WRTBAD( 'Rshell' )

      IF( RCORE.LE.0.0 .OR. RCORE.GT.RSHELL ) 
     &    INPERR = WRTBAD( 'Rcore' )

      IF( REAL(RINDSH).LE.0.0 .OR. AIMAG(RINDSH).GT.0.0 )
     &    INPERR = WRTBAD( 'RindSh' )

      IF( REAL(RINDCO).LE.0.0 .OR. AIMAG(RINDCO).GT.0.0 ) 
     &    INPERR = WRTBAD( 'RindCo' )

      IF( NUMANG.LT.0 ) INPERR = WRTBAD( 'NumAng' )

      IF( NUMANG.GT.MXANG ) INPERR = WRTDIM( 'MxAng', NUMANG )

      IF( NUMANG.GT.MAXANG ) INPERR = WRTDIM( 'MaxAng', NUMANG )

      IF( NMX1 + 1 .GT. LL ) INPERR = WRTDIM( 'LL', NMX1 + 1 )

      DO 10 J = 1, NUMANG
         IF( MU(J).LT.- TOLER .OR. MU(J).GT. 1.0+TOLER )
     &        INPERR = WRTBAD( 'MU' )
   10 CONTINUE

      IF( INPERR ) then
		Err=1
		return
		CALL ERRMSG(
     &    'MIELAY--Input argument errors.  Aborting...', .True. )
	endif
	allocate(w(3,NMX1+1))
	allocate(acap(NMX1+1))


      K1     = RINDCO*WVNO
      K2     = RINDSH*WVNO
      K3     = DCMPLX( WVNO )
      Z( 1 ) = RINDSH*XSHELL
      Z( 2 ) = XSHELL
      Z( 3 ) = RINDCO*XCORE
      Z( 4 ) = RINDSH*XCORE
      X1     =  DBLE( Z(1) )
      Y1     = DIMAG( Z(1) )
      X4     =  DBLE( Z(4) )
      Y4     = DIMAG( Z(4) )
      RX     = ONE / XSHELL

c                                ** Down-recurrence for A function
      ACAP( NMX1 + 1 ) = CZERO
      DO 20 M = 1, 3
         W( M, NMX1 + 1 ) = CZERO
   20 CONTINUE

      RRFX  = ONE / ( RINDSH*XSHELL)
      DO 40 NN = NMX1, 1, - 1

         ACAP( NN ) = ( ( NN + 1)*RRFX ) -
     &                ONE / ( ( (NN + 1)*RRFX) + ACAP( NN + 1) )

         DO 30 M = 1, 3

            W( M, NN ) = ( ( NN + 1) / Z( M + 1) ) - 
     &                   ONE / ( ( (NN + 1)/Z(M + 1)) + W( M, NN + 1) )

   30    CONTINUE

   40 CONTINUE


      DO 50 J = 1, NUMANG

         SI2THT( J ) = ONE - MU( J )**2
         PI( J, 1 ) = ZERO
         PI( J, 2 ) = ONE
         TAU( J, 1 ) = ZERO
         TAU( J, 2 ) = MU( J )

   50 CONTINUE

c                          ** Initialization of homogeneous sphere

      T( 1 ) = COS( XSHELL )
      T( 2 ) = SIN( XSHELL )
      WM1      = DCMPLX( T(1), - T(2) )
      WFN( 1 ) = DCMPLX( T(2), T(1) )
      TA( 1 ) = T( 2 )
      TA( 2 ) = T( 1 )
      WFN( 2 ) = RX*WFN( 1 ) - WM1
      TA( 3 ) =  DBLE( WFN(2) )
      TA( 4 ) = DIMAG( WFN(2) )

c                      ** Initialization procedure for stratified sphere
      N      = 1
      SINX1  = SIN( X1 )
      SINX4  = SIN( X4 )
      COSX1  = COS( X1 )
      COSX4  = COS( X4 )
      EY1    = EXP( Y1 )
      E2Y1   = EY1**2
      EY4    = EXP( Y4 )
      EY1MY4 = EXP( Y1 - Y4 )
      EY1PY4 = EY1*EY4
      AA     = SINX4*( EY1PY4 + EY1MY4 )
      BB     = COSX4*( EY1PY4 - EY1MY4 )
      CC     = SINX1*( E2Y1 + ONE )
      DD     = COSX1*( E2Y1 - ONE )
      DENOM  = ONE + E2Y1*( 4.0D0*SINX1**2 - TWO + E2Y1 )
      DUMMY  = DCMPLX( ( AA*CC + BB*DD) / DENOM,
     &                 ( BB*CC - AA*DD) / DENOM )
      DUMMY  = DUMMY*( ACAP(N) + N / Z(1) ) / ( W(3, N) + N / Z(4) )
      DUMSQ  = DUMMY**2

      P24H24 = 0.5D0 + DCMPLX( SINX4**2 - 0.5D0, COSX4*SINX4 )*EY4**2
      P24H21 = 0.5D0*DCMPLX( SINX1*SINX4 - COSX1*COSX4,
     &                       SINX1*COSX4 + COSX1*SINX4 )*EY1PY4
     &       + 0.5D0*DCMPLX( SINX1*SINX4 + COSX1*COSX4,
     &                     - SINX1*COSX4 + COSX1*SINX4 )*EY1MY4
      DH1    = Z( 1 ) / ( ONE + CI*Z( 1) ) - ONE / Z( 1 )
      DH2    = Z( 2 ) / ( ONE + CI*Z( 2) ) - ONE / Z( 2 )
      DH4    = Z( 4 ) / ( ONE + CI*Z( 4) ) - ONE / Z( 4 )
      P24H24 = P24H24 / ( ( DH4 + N/Z(4))*( W(3, N) + N/Z(4)) )
      P24H21 = P24H21 / ( ( DH1 + N/Z(1))*( W(3, N) + N/Z(4)) )

      U( 1 ) = K3*ACAP( N ) - K2*W( 1, N )
      U( 2 ) = K3*ACAP( N ) - K2*DH2
      U( 3 ) = K2*ACAP( N ) - K3*W( 1, N )
      U( 4 ) = K2*ACAP( N ) - K3*DH2
      U( 5 ) = K1*W( 3, N ) - K2*W( 2, N )
      U( 6 ) = K2*W( 3, N ) - K1*W( 2, N )
      U( 7 ) = - CI*( DUMMY*P24H21 - P24H24 )
      U( 8 ) = TA( 3 ) / WFN( 2 )

      ACOE  = U( 8 )*( U(1)*U(5)*U(7) + K1*U(1) - DUMSQ*K3*U(5) ) /
     &               ( U(2)*U(5)*U(7) + K1*U(2) - DUMSQ*K3*U(5) )

      BCOE  = U( 8 )*( U(3)*U(6)*U(7) + K2*U(3) - DUMSQ*K2*U(6) ) /
     &               ( U(4)*U(6)*U(7) + K2*U(4) - DUMSQ*K2*U(6) )

      ACOEM1 = ACOE
      BCOEM1 = BCOE
      ARE    =  DBLE( ACOE )
      AIM    = DIMAG( ACOE )
      BRE    =  DBLE( BCOE )
      BIM    = DIMAG( BCOE )

      DQEXT  = 3.D0*( ARE + BRE )
      DQSCA  = 3.D0*( ARE**2 + AIM**2 + BRE**2 + BIM**2 )
      DGQSC  = ZERO
      SBACK  = 3.D0*( ACOE - BCOE )
      RMM    = ONE

      AC  = 1.5D0*ACOE
      BC  = 1.5D0*BCOE
      DO 60 J = 1, NUMANG

         S1( J, 1 ) = AC*PI( J, 2 ) + BC*TAU( J, 2 )
         S1( J, 2 ) = AC*PI( J, 2 ) - BC*TAU( J, 2 )
         S2( J, 1 ) = BC*PI( J, 2 ) + AC*TAU( J, 2 )
         S2( J, 2 ) = BC*PI( J, 2 ) - AC*TAU( J, 2 )

   60 CONTINUE

c ***************** Start of Mie summing loop ******************

      N  = 2
   70 CONTINUE
c                              ** Recurrences for functions little-pi,
c                                 little-tau of Mie theory
      T( 1 ) = 2*N - 1
      T( 2 ) = N - 1
      DO 80 J = 1, NUMANG

         PI( J, 3 ) = ( T( 1)*PI( J, 2)*MU( J) - N*PI( J, 1) ) / T( 2 )

         TAU( J, 3 ) = MU( J )*( PI( J, 3) - PI( J, 1) ) -
     &                 T( 1 )*SI2THT( J )*PI( J, 2 ) + TAU( J, 1 )

   80 CONTINUE

c                                 ** Here set up homogeneous sphere
      WM1    = WFN( 1 )
      WFN( 1 ) = WFN( 2 )
      WFN( 2 ) = T( 1 )*RX*WFN( 1 ) - WM1
      TA( 1 ) =  DBLE( WFN( 1) )
      TA( 2 ) = DIMAG( WFN( 1) )
      TA( 3 ) =  DBLE( WFN( 2) )
      TA( 4 ) = DIMAG( WFN( 2) )

c                                 ** Here set up stratified sphere

      DH1    = - N / Z( 1 ) + ONE / ( N / Z( 1) - DH1 )
      DH2    = - N / Z( 2 ) + ONE / ( N / Z( 2) - DH2 )
      DH4    = - N / Z( 4 ) + ONE / ( N / Z( 4) - DH4 )
      P24H24 = P24H24 / ( ( DH4 + N/Z(4))*( W(3, N) + N/Z(4)) )
      P24H21 = P24H21 / ( ( DH1 + N/Z(1))*( W(3, N) + N/Z(4)) )
      DUMMY  = DUMMY*( ACAP(N) + N / Z(1) ) / ( W(3, N) + N / Z(4) )
      DUMSQ  = DUMMY**2

      U( 1 ) = K3*ACAP( N ) - K2*W( 1, N )
      U( 2 ) = K3*ACAP( N ) - K2*DH2
      U( 3 ) = K2*ACAP( N ) - K3*W( 1, N )
      U( 4 ) = K2*ACAP( N ) - K3*DH2
      U( 5 ) = K1*W( 3, N ) - K2*W( 2, N )
      U( 6 ) = K2*W( 3, N ) - K1*W( 2, N )
      U( 7 ) = - CI*( DUMMY*P24H21 - P24H24 )
      U( 8 ) = TA( 3 ) / WFN( 2 )

      ACOE  = U( 8 )*( U(1)*U(5)*U(7) + K1*U(1) - DUMSQ*K3*U(5) ) /
     &               ( U(2)*U(5)*U(7) + K1*U(2) - DUMSQ*K3*U(5) )

      BCOE  = U( 8 )*( U(3)*U(6)*U(7) + K2*U(3) - DUMSQ*K2*U(6) ) /
     &               ( U(4)*U(6)*U(7) + K2*U(4) - DUMSQ*K2*U(6) )
      ARE  =  DBLE( ACOE )
      AIM  = DIMAG( ACOE )
      BRE  =  DBLE( BCOE )
      BIM  = DIMAG( BCOE )

c                           ** Increment sums for efficiency factors

      AM1RE  =  DBLE( ACOEM1 )
      AM1IM  = DIMAG( ACOEM1 )
      BM1RE  =  DBLE( BCOEM1 )
      BM1IM  = DIMAG( BCOEM1 )
      T( 4 ) = (2*N - ONE) / ( N*(N - ONE) )
      T( 2 ) = (N - ONE)*(N + ONE) / N
      DGQSC  = DGQSC + T( 2 )*( AM1RE*ARE + AM1IM*AIM +
     &                          BM1RE*BRE + BM1IM*BIM ) +
     &                 T( 4 )*( AM1RE*BM1RE + AM1IM*BM1IM )

      T( 3 )  = 2*N + 1
      DQEXT   = DQEXT + T( 3 )*( ARE + BRE )
      T( 4 )  = ARE**2 + AIM**2 + BRE**2 + BIM**2
      DQSCA   = DQSCA + T( 3 )*T( 4 )
      RMM     = - RMM
      SBACK  = SBACK + T( 3 ) * RMM *( ACOE - BCOE )

      T( 2 ) = N*( N + 1 )
      T( 1 ) = T( 3 ) / T( 2 )

      AC  = T( 1 )*ACOE
      BC  = T( 1 )*BCOE
      DO 90 J = 1, NUMANG
         S1( J, 1 ) = S1( J, 1 ) + AC*PI( J, 3 ) + BC*TAU( J, 3 )
         S2( J, 1 ) = S2( J, 1 ) + BC*PI( J, 3 ) + AC*TAU( J, 3 )
   90 CONTINUE

c                               ** Scattering matrix elements for
c                                  supplements of 0-90 degree scattering
c                                  angles submitted by user
      IF( MOD(N, 2).EQ.0 ) THEN

         DO 100 J = 1, NUMANG
            S1( J, 2 ) = S1( J, 2 ) - AC*PI( J, 3 ) + BC*TAU( J, 3 )
            S2( J, 2 ) = S2( J, 2 ) - BC*PI( J, 3 ) + AC*TAU( J, 3 )
  100    CONTINUE

      ELSE

         DO 110 J = 1, NUMANG
            S1( J, 2 ) = S1( J, 2 ) + AC*PI( J, 3 ) - BC*TAU( J, 3 )
            S2( J, 2 ) = S2( J, 2 ) + BC*PI( J, 3 ) - AC*TAU( J, 3 )
  110    CONTINUE

      END IF

c                                      ** Test for convergence of sums
      IF( T(4).GE.1.0D-14 ) THEN

         N  = N + 1

         IF( N.GT.NMX2 .or..not.N.gt.0) then
		Err=1
		return
	   	CALL ERRMSG(
     &       'MIELAY--Dimensions for W,ACAP not enough. Suggest'//
     &       ' get detailed output, modify routine', .True. )
	   endif

         DO 120 J = 1, NUMANG

            PI( J, 1 ) = PI( J, 2 )
            PI( J, 2 ) = PI( J, 3 )
            TAU( J, 1 ) = TAU( J, 2 )
            TAU( J, 2 ) = TAU( J, 3 )

  120    CONTINUE

         ACOEM1 = ACOE
         BCOEM1 = BCOE

         GO TO 70

      END IF

c ***************** End of summing loop ******************

c                            ** Transform complex scattering amplitudes
c                               into elements of real scattering matrix

      DO 140 J = 1, NUMANG

         DO 130 K = 1, 2

            M1( J, K ) = DBLE( S1(J, K) )**2 + DIMAG( S1(J, K) )**2
            M2( J, K ) = DBLE( S2(J, K) )**2 + DIMAG( S2(J, K) )**2
            S21( J, K ) = DBLE(  S1(J, K) )*DBLE(  S2(J, K) ) +
     &                    DIMAG( S1(J, K) )*DIMAG( S2(J, K) )
            D21( J, K ) = DIMAG( S1(J, K) )*DBLE( S2(J, K) ) -
     &                    DIMAG( S2(J, K) )*DBLE( S1(J, K) )

  130    CONTINUE

  140 CONTINUE


      T( 1 ) = TWO*RX**2
      QEXT   = T( 1 )*DQEXT
      QSCA   = T( 1 )*DQSCA
      GQSC   = TWO*T( 1 )*DGQSC
      SBACK  = 0.5*SBACK
      QBS    = ( DBLE(SBACK)**2 + DIMAG(SBACK)**2 ) / (PINUM*XSHELL**2)

	return
      END




      SUBROUTINE ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error
c        after making symbolic dump (machine-specific)

c       Provenance:  the 3 error-handlers ErrMsg, WrtBad, WrtDim are
c                    borrowed from MIEV, the Wiscombe Mie program

c     .. Scalar Arguments ..

      CHARACTER MESSAG*( * )
      LOGICAL   FATAL
c     ..
c     .. Local Scalars ..

      LOGICAL   MSGLIM
      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

c cccc EXTERNAL  SYMDUMP
c     ..
c     .. Save statement ..

      SAVE      MAXMSG, NUMMSG, MSGLIM
c     ..
c     .. Data statements ..

      DATA      NUMMSG / 0 / , MAXMSG / 100 / , MSGLIM / .FALSE. /
c     ..

      IF( FATAL ) THEN

c         WRITE( *, '(//,2A,//)' ) ' ****** ERROR *****  ', MESSAG

c                                 ** Example symbolic dump call for Cray
c cccc    CALL SYMDUMP( '-B -c3' )

c         write(*,*) 'I should actually stop, but... whatever'!STOP

      END IF


      NUMMSG = NUMMSG + 1

      IF( MSGLIM ) RETURN

      IF( NUMMSG.LE.MAXMSG ) THEN

c         WRITE( *, '(/,2A,/)' ) ' ****** WARNING *****  ', MESSAG

      ELSE

c         WRITE( *, '(//,A,//)' )
c     &      ' ****** TOO MANY WARNING MESSAGES --  ' //
c     &      'They will no longer be printed *******'

         MSGLIM = .True.

      END IF

      END

      LOGICAL FUNCTION WrtBad( VARNAM )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

c     .. Scalar Arguments ..

      CHARACTER VARNAM*( * )
c     ..
c     .. Local Scalars ..

      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Save statement ..

      SAVE      NUMMSG, MAXMSG
c     ..
c     .. Data statements ..

      DATA      NUMMSG / 0 / , MAXMSG / 50 /
c     ..

      WRTBAD = .TRUE.
      NUMMSG = NUMMSG + 1
c      WRITE( *, '(3A)' ) ' ****  Input variable  ', VARNAM,
c     &   '  in error  ****'

      IF( NUMMSG.EQ.MAXMSG ) CALL ERRMSG(
     &    'Too many input errors.  Aborting...', .TRUE. )

      END

      LOGICAL FUNCTION WrtDim( DIMNAM, MINVAL )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

c     .. Scalar Arguments ..

      CHARACTER DIMNAM*( * )
      INTEGER   MINVAL
c     ..

c      WRITE( *, '(3A,I7)' ) ' ****  Symbolic dimension  ',
c     &   DIMNAM, '  should be increased to at least ', MINVAL

      WRTDIM = .TRUE.

      END

	
