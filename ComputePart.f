	subroutine ComputePart(C,ii,isize,computelamcloud)
	use GlobalSetup
	use Constants
	use AImod
	IMPLICIT NONE

	type(CloudType) C,Cdust
	integer MAXMAT
	parameter(MAXMAT=60)
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
	real QEXT, QSCA, QBS, GQSC,wvno,scale,gsca,rcore4
	character*3 meth
	character*500 filename(MAXMAT),grid,tmp,tmp2,partfile,lnkfile

	real*8 rmie,lmie,e1mie,e2mie,csmie,cemie,KR,theta,dummy,amin,amax,camie,fmie,gmie,rcore
	logical truefalse,checkparticlefile,lnkloglog
	integer abun_in_name,LL,LLmax
	parameter(abun_in_name=2)
	real*8 Kabs(nlam+1),Ksca(nlam+1),Kext(nlam+1),lgrid(nlam+1),G(nlam+1)
	logical computelamcloud(nlam)
	real*8 csmie_fcomp,cemie_fcomp,gasdev

	real*8,allocatable :: Mief11(:),f11(:,:)
	real,allocatable :: mu(:),M1(:,:),M2(:,:),S21(:,:),D21(:,:)

	write(meth,100)
100	format('DHS')

	na=180

	allocate(e1(MAXMAT,C%nlam))
	allocate(e2(MAXMAT,C%nlam))

	allocate(frac(MAXMAT))
	allocate(rho(MAXMAT))
	allocate(f11(C%nlam,na))

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
		frac(C%nmat+1)=C%porosity(isize)
		frac(1:C%nmat)=frac(1:C%nmat)*(1d0-C%porosity(isize))
		e1d(C%nmat+1,1:C%nlam)=1d0
		e2d(C%nmat+1,1:C%nlam)=0d0
		rho(1:C%nmat)=C%rho_mat(1:C%nmat)
		rho(C%nmat+1)=0d0
		do i=1,nm
			do j=1,C%nmat
				e1d(j,1:C%nlam)=C%e1(j,minval((/i,C%nax(j)/)),1:C%nlam)
				e2d(j,1:C%nlam)=C%e2(j,minval((/i,C%nax(j)/)),1:C%nlam)
			enddo
!$OMP PARALLEL IF(useomp)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(j,e1blend,e2blend)
!$OMP& SHARED(C,i,e1,e2,frac,e1d,e2d)
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC,1)
			do j=1,C%nlam
				call Blender(e1d(1:C%nmat+1,j),e2d(1:C%nmat+1,j),frac(1:C%nmat+1),C%nmat+1,e1blend,e2blend)
				e1(i,j)=e1blend
				e2(i,j)=e2blend
			enddo
!$OMP END DO
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

	do i=1,nlam
		do j=1,na
			f11(i,j)=0d0
		enddo
	enddo

	fmie=maxf
!$OMP PARALLEL IF(useomp)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,csca0,cabs0,cext0,theta,i,l,tot,k,Err,spheres,toolarge,
!$OMP&         rad,wvno,m,r1,rcore,qext,qsca,qbs,gqsc,rmie,lmie,e1mie,e2mie,
!$OMP&         csmie,cemie,camie,tot2,j,gsca,gmie,MieF11,mu,M1,M2,S21,D21,rcore4)
!$OMP& SHARED(C,na,nm,ns,frac,minlog,maxlog,f,e1,e2,wf,isize,computelamcloud,Mass,
!$OMP&        pow,lgrid,rho,nf,r0,nr0,Kabs,Ksca,Kext,nlam,fmie,useDLMie,
!$OMP&		  DLMie_e1max,DLMie_e2max,G,f11,min,anisoscattstar)
	allocate(Mief11(na))
	allocate(mu(na))
	allocate(M1(na,2))
	allocate(M2(na,2))
	allocate(S21(na,2))
	allocate(D21(na,2))
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
	gsca=0d0

	do i=1,na/2
		theta=(real(i)-0.5)/real(na/2)*3.1415926536/2d0
		mu(i)=cos(theta)
	enddo

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
	e1mie=e1(l,ilam)
	e2mie=e2(l,ilam)
	if(useDLMie.and.e1mie.lt.DLmie_e1max.and.e2mie.lt.DLmie_e2max) then
		rmie=r1
		lmie=lgrid(ilam)
		call DLMie(rmie,lmie,e1mie,e2mie,fmie,camie,csmie)
		cext0=cext0+nr0(l,k)*(camie+csmie)
		csca0=csca0+nr0(l,k)*csmie
	   	cabs0=cabs0+nr0(l,k)*camie
	   	gsca=0d0
	else
	do i=1,nf
		rad=r1/(1d0-f(i))**(1d0/3d0)
		m=dcmplx(e1(l,ilam),-e2(l,ilam))
		wvno=2d0*3.1415926536/lgrid(ilam)

		if(f(i).eq.0d0) then
			spheres=1
			goto 20
		endif
		if(r1*wvno.gt.1000d0.and.anisoscattstar) then
			toolarge=1
			goto 20
		endif
		rcore=rad*f(i)**(1d0/3d0)
		rmie=rad
		lmie=lgrid(ilam)
		e1mie=e1(l,ilam)
		e2mie=e2(l,ilam)

		if(anisoscattstar) then
			if(rmie/lmie.gt.1000d0) lmie=rmie/1000d0
			wvno=2d0*3.1415926536/lmie
			rcore4=rcore
			call DMiLay(rcore4, rad, WVNO, m, min, MU,
     &                   NA/2, QEXT, QSCA, QBS, GQSC, 
     &                   M1, M2, S21, D21, NA ,Err)
			cemie=qext*pi*rad**2
			csmie=qsca*pi*rad**2
			do j=1,na/2
				Mief11(j)=(M2(j,1)+M1(j,1))/csmie/wvno**2*2d0*pi
				Mief11(na-j+1)=(M2(j,2)+M1(j,2))/csmie/wvno**2*2d0*pi
			enddo
		else
			if(rmie/lmie.gt.10d0) lmie=rmie/10d0
			wvno=2d0*3.1415926536/lmie
			call callBHCOAT(rmie,rcore,lmie,e1mie,e2mie,csmie,cemie,gmie,Err)
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
				if(anisoscattstar) then
					if(rmie/lmie.gt.1000d0) lmie=rmie/1000d0
					call MeerhoffMie(rmie,lmie,e1mie,e2mie,csmie,cemie
     &								,Mief11,na)
				else
					if(rmie/lmie.gt.10d0) lmie=rmie/10d0
					call callBHMIE(rmie,lmie,e1mie,e2mie,csmie,cemie,gmie)
				endif
			endif
		endif

		do j=1,na
			f11(ilam,j)=f11(ilam,j)+wf(i)*nr0(l,k)*Mief11(j)*csmie
		enddo

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
	if(anisoscattstar) then
		do ilam=1,C%nlam
			C%F11(isize,ilam,1:180)=f11(ilam,1:180)
			tot=0d0
			tot2=0d0
			do j=1,180
				tot=tot+C%F11(isize,ilam,j)*sin(pi*(real(j)-0.5)/real(na))
				tot2=tot2+sin(pi*(real(j)-0.5)/real(na))
			enddo
			if(tot.gt.0d0) then
				C%F11(isize,ilam,1:180)=C%F11(isize,ilam,1:180)*tot2/tot
				tot2=0d0
				tot=0d0
				do j=1,180
					tot=tot+C%F11(isize,ilam,j)*sin(pi*(real(j)-0.5)/real(na))*cos(pi*(real(j)-0.5)/real(na))
					tot2=tot2+C%F11(isize,ilam,j)*sin(pi*(real(j)-0.5)/real(na))
				enddo
				C%G(isize,ilam)=tot/tot2
			else
				C%F11(isize,ilam,1:180)=1d0
				C%G(isize,ilam)=0d0
			endif
		enddo
	else
		C%F11(isize,1:C%nlam,1:180)=1d0
		C%G(isize,1:C%nlam)=0d0
	endif
	if(outputopacity) call WriteOpacity(isize,"wolk",1d4/lgrid(1:C%nlam),Kabs(1:C%nlam),C%nlam,1)

	
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
		if(i.gt.n) return
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

	subroutine BlenderOld(m,abun,nm,e1out,e2out)
	IMPLICIT NONE
	integer,intent(in) :: nm
	real*8,intent(in) :: abun(nm)
	complex*16,intent(in) :: m(nm)
	real*8,intent(out) :: e1out,e2out
	integer j,iter
	complex*16,save :: mm,me,sum,m2(100)
	logical,save :: doit(100)
!$OMP THREADPRIVATE(mm,me,sum,m2,doit)

c	call BlenderCDE(m,abun,nm,e1out,e2out)
c	return

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
	

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

	subroutine Blender(e1in,e2in,abun,nm,e1out,e2out)
	IMPLICIT NONE
	integer nm,j,iter,niter,iguess,nguess,i1,i2
	real e1in(nm),e2in(nm)
	real*8 e1out,e2out,abun(nm)
	complex*16 mm,m(nm),me,fx,dfx,dme,m2
	complex*16 a1,a2,a3,a4,a5,mguess(10)
	logical CDE
	real*8 fmax,eps,f,amax,e1max,e2max,minfx

	do j=1,nm
		m(j)=dcmplx(e1in(j),e2in(j))
	enddo
	call BlenderOld(m,abun,nm,e1out,e2out)
	return

	eps=1d-5
	niter=100
	f=1d0

	me=0d0
	e1max=0d0
	e2max=0d0
	amax=0d0
	do j=1,nm
		me=me+m(j)**(2d0/3d0)*abun(j)
		if(real(m(j)).gt.e1max) e1max=real(m(j))
		if(imag(m(j)).gt.e2max) e2max=imag(m(j))
		if(abun(j).gt.amax) then
			mguess(1)=m(j)
			amax=abun(j)
		endif
	enddo
	me=me**(3d0/2d0)
	mguess(2)=me

	mguess(4)=cmplx(1d0,0d0)
	mguess(5)=cmplx(e1max,e2max)
	mguess(6)=cmplx((e1max-1d0)/2d0+1d0,e2max/2d0)
	nguess=6
	
	fmax=0.0001d0
	CDE=.false.

	iguess=1

1	continue

	if(iguess.eq.3) then
		mm=mguess(2)
		do iter=1,niter
			fx=0d0
			do j=1,nm
				fx=fx+((m(j)**2-mm)/(m(j)**2+2d0*mm))*abun(j)
			enddo
			me=(2d0*fx+1d0)/(1d0-fx)
			me=mm*me
			mm=me
			if(cdabs(fx).lt.eps) exit
		enddo
		mguess(3)=cdsqrt(me)
	endif
	me=mguess(iguess)
	do iter=1,niter
		fx=0d0
		dfx=0d0
		do j=1,nm
			m2=(m(j)/me)**2
			if(CDE) then
!	CDE
				if(abs(m2-1d0).lt.1d-2) then
					fx=fx+abun(j)*(-1d0+(m2-1d0)/2d0-(m2-1d0)**2/6d0)
					dfx=dfx-abun(j)*(m(j)*(4d0*me-m(j))/(3d0*me**3))
				else
					fx=fx+abun(j)*(2d0*m2*log(m2)/(m2-1d0)-2d0)
					dfx=dfx+abun(j)*(4d0*m(j)**2*((log(m2)+1d0)*me**2-m(j)**2)/(me*(me**2-m(j)**2)**2))
				endif
!	DHS
			else
				m2=(m(j)/me)**2
				a1=m2+2d0
				a2=m2-1d0
				a3=2d0*m2+1d0
				a4=6d0*fmax*a2**2*(m2+1d0)
				a5=a1*a3-2d0*fmax*a2**2
				if(abs(a2).lt.1d-2) then
					fx=fx+abun(j)*3d0*a2/a1
					dfx=dfx+abun(j)*((18d0*m2*(a1-3d0*(m2+1d0)))/(a1**2*a3)+(18d0*a2**2*fmax**2*m2*(a1-6d0*(m2+1d0)))/(a1**3*a3**2))
				else
					fx=fx+abun(j)*((3d0*a3)/(2d0*a2*fmax))*cdlog((a1*a3)/a5)
					dfx=dfx+abun(j)*9d0*m2*(a1*a5*cdlog(a1*a3/a5)-a4)/(fmax*a1*a5*a2**2)
				endif
			endif
		enddo
		dme=fx/dfx
		if(abs(dme/me).lt.eps) exit
2		me=me-dme*f
		e1out=dreal(me)
		e2out=dimag(me)
		if(e1out.lt.0d0.or.e1out.gt.e1max.or.e2out.lt.0d0.or.e2out.gt.e2max) then
			me=me+dme*f
			f=f/2d0
			goto 2
		endif
	enddo
	e1out=dreal(me)
	e2out=dimag(me)
	if(e1out.lt.0d0.or.e2out.lt.0d0.or.abs(dme/me).gt.1d-2) then
		if(f.gt.0.1) then
			f=f/2d0
			goto 1
		else
			f=1d0
			if(iguess.lt.nguess) then
				iguess=iguess+1
				goto 1
			else
				print*,'failed!'
				me=mguess(2)
				e1out=dreal(me)
				e2out=dimag(me)
			endif		
		endif
	endif

	return
	end
	


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

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


	subroutine callBHMIE(rmie,lmie,e1mie,e2mie,csmie,cemie,gmie)
	IMPLICIT NONE
	real*8 rmie,lmie,e1mie,e2mie,csmie,cemie,gmie
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
	gmie=GSCA
	
	return
	end

	subroutine callBHCOAT(rmie,rcore,lmie,e1mie,e2mie,csmie,cemie,gmie,Err)
	IMPLICIT NONE
	real*8 rmie,lmie,e1mie,e2mie,csmie,cemie,rcore,gmie
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
	gmie=GSCA
     	
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
	double complex,allocatable,save :: W(:,:),acap(:)
!$OMP THREADPRIVATE(w,acap)
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
!$OMP THREADPRIVATE(PINUM,PASS1)

c     ..
c     .. Data statements ..

      DATA      PASS1 / .True. / , TOLER / 1.D-6 / ,
     &          CZERO / ( 0.D0, 0.D0 ) / , CI / ( 0.D0, 1.D0 ) /
c     ..


      IF( PASS1 ) THEN
		allocate(w(3,LL))
		allocate(acap(LL))

         PINUM  = TWO*ASIN( ONE )
         PASS1  = .False.

      END IF

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
!$OMP THREADPRIVATE(MAXMSG, NUMMSG, MSGLIM)
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
!$OMP THREADPRIVATE(NUMMSG, MAXMSG)
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

	








































	subroutine MeerhoffMie(rad,lam,e1,e2,Csca,Cext,F11,na)
	IMPLICIT NONE
	integer na,nangle
	real*8 rad,lam,e1,e2,csca,cext,F11(na)

	integer nparts,develop,nsubr(1),ngaur(1),idis(1),ndis,i
	real*8 delta,cutoff,thmin,thmax,step,wavel(1),Rem(1),fImm(1)
	real*8 par1(1),par2(1),par3(1),rdis(1,300),nwrdis(1,300)
	real*8 F(4,6000),theta(6000)
	
	nangle=na
	nparts=1
	develop=0
	delta=1d-8
	cutoff=1d-8
	thmin=180d0*(real(1)-0.5)/real(na)
	thmax=180d0*(real(na)-0.5)/real(na)
	step=(thmax-thmin)/real(na-1)
	wavel=lam
	Rem=e1
	fImm=e2
	nsubr=1
	ngaur=1
	idis=0
	par1=rad
	par2=0d0
	par3=0d0
	ndis=1
	rdis(1,1)=rad
	nwrdis(1,1)=1d0

      call mie(nparts, develop, delta, cutoff, thmin, thmax, step
     +               , wavel, Rem, fImm, nsubr , ngaur, idis
     +               , par1, par2, par3, ndis, rdis, nwrdis
     +               , nangle, theta, F ,Cext, Csca)

	do i=1,nangle
		F11(i)=F(1,nangle-i+1)
	enddo

	return
	end
	



      subroutine mie(nparts, develop, delta, cutoff, thmin, thmax, step
     +               , wavel, Rem, fImm, nsubr , ngaur, idis
     +               , par1, par2, par3, ndis, rdis, nwrdis
     +               , nangle, theta, F ,outCext, outCsca)


      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter ( NDang=6000, NDcoef=6000, NDpart = 4, nrunit=88 )
      parameter ( NDdis=300 )
      double precision lambda, nr, ni, miec, nwrdis,outCext,outCsca
      integer develop
      double complex m
      character*60 s1, s2, s3, s4, s5, s6
      character*10 scfile
      dimension F(4,NDang)
     +        , miec(13),u(NDang),w(NDang),coefs(4,4,0:NDcoef)
     +        , wavel(*), Rem(*), fImm(*), nsubr(*)
     +        , ngaur(*), idis(*)
     +        , par1(*), par2(*), par3(*)
     +        , rdis(1,*), nwrdis(1,*)
     +        , scfile(NDpart), theta(NDang)
      pi     = dacos(-1.D0)
      radtod = 180.D0/pi
      scfile(1) = 'mie1_sc  '
      scfile(2) = 'mie2_sc  '
      scfile(3) = 'mie3_sc  '
      scfile(4) = 'mie4_sc  '
************************************************************************
*  Start a loop over different 'particles' referring to the different  *
*  cases specified in the input                                        *
************************************************************************
      do 50 iparts=1,nparts
c      write(7,920) iparts
************************************************************************
*  In case of development in generalized spherical functions a file    *
*  is opened on which the expansion coefficients will be written       *
************************************************************************
c      if (develop .ne. 0)
c     +    call opensc(scfile(iparts),nrunit,'new')
************************************************************************
*  Determine the integration bounds for the radius integration         *
************************************************************************
      call rminmax( idis(iparts), nsubr(iparts), ngaur(iparts)
     +   , par1(iparts), par2(iparts), par3(iparts), cutoff
     +   , rmin, rmax
     +   , iparts, ndis, rdis, nwrdis )
*
      lambda = wavel(iparts)
      nr     = Rem(iparts)
      ni     = fImm(iparts)
      m      = dcmplx(nr,-ni)
************************************************************************
*  Calculate the scattering matrix : this is most of the work !        *
************************************************************************
      call scatmat( u, w, m, wavel(iparts), idis(iparts)
     +            , thmin, thmax, step, develop
     +            , nsubr(iparts), ngaur(iparts), rmin, rmax
     +            , par1(iparts), par2(iparts), par3(iparts)
     +            , delta, F, miec, nangle 
     +            , ndis, rdis, nwrdis, iparts )
************************************************************************
*  Test if the number of coefficients needed fits in the array coefs   *
*  The number of coefficients is equal to the number of integration    *
*  points in scattering angle !! (think about this !)                  *
************************************************************************
      ncoef = nangle
      if ((ncoef .gt. NDcoef) .and. (develop .eq. 1)) then
         write(*,*) ' main: too many coefficients needed :',ncoef
         write(*,*) '       maximum array size NDcoef = ',NDcoef
         write(*,*) '       Therefore I cannot do the expansion.'
         develop=0
         call clossc(nrunit)
      endif
************************************************************************
*  If required, calculate the expansion coefficients                   *
************************************************************************
      if (develop.eq.1) then
         call devel(ncoef,nangle,u,w,F,coefs)
         npunt = idnint((thmax-thmin)/step)+1
         if (npunt .gt. NDang) then
            write(*,*) ' main: too many angles for table npunt=',npunt
            write(*,*) '       maximum array size NDang = ',NDang
            npunt=61
            write(*,*) '      I have set the number of angles to ',npunt
         endif
      endif
************************************************************************
*  Print to standard output the following data :                       *
*                                                                      *
*  miec(1) = Csca : the average scattering cross section               *
*  miec(2) = Cext : the average extinction cross section               *
*  miec(3) = Qsca : the scattering efficiency factor                   *
*  miec(4) = Qext : the extinction efficiency factor                   *
*  miec(5) = a    : the single scattering albedo                       *
*  miec(6) = G    : the average geometrical cross section              *
*  miec(7) = reff : the average effective radius                       *
*  miec(8) = xeff : the average effective size parameter               *
*  miec(9) = num  : the integrated number of particles                 *
*  miec(10)= V    : the average volume                                 *
************************************************************************
c      write(7,922) miec(1)
c      write(7,923) miec(2)
c      write(7,924) miec(3)
c      write(7,925) miec(4)
c      write(7,926) miec(5)
c      write(7,927) miec(6)
c      write(7,928) miec(7)
c      write(7,929) miec(8)
c      write(7,930) miec(9)
c      write(7,931) miec(10)

	outCsca=miec(1)
	outCext=miec(2)

c      if (develop .eq. 0) goto 1020

************************************************************************
*  In case a development in GSF is required, prepare the heading of    *
*  coefficient file.                                                   *
************************************************************************
* ruler: '.........1.........2.........3.........4.........5.........6'
c      s1='PRODUCED BY THE MEERHOFF MIE PROGRAM VERSION 3.1            '
c      write(s2,1001) wavel(iparts), nr, ni
c 1001 format('lambda= ',f11.7,' Re(m) = ',f11.7,' Im(m) = ',f11.7)
************************************************************************
c      if (idis(iparts) .ne. 0) then
c          write(s5,1007) rmin, rmax
c          write(s6,1009) nsubr(iparts), ngaur(iparts)
c      endif
c 1007 format('rmin  = ',f11.7,' rmax  = ',f11.7,'                     ')
c 1009 format('r integration with ',i3,' interval(s) of '
c     +                                            ,i3,' Gauss points  ')
************************************************************************
c      if (idis(iparts) .eq. 0) then
c      s5='                                                            '
c      s6='                                                            '
c      s3='No size distribution (single particle)                      '
c      write(s4,1003) par1(iparts)
c 1003 format('Radius = ',1f11.7
c     +                       ,'                                       ')
************************************************************************
c      else if (idis(iparts) .eq. 1) then
c      s3='Size distribution : two parameter gamma                     '
c      write(s4,1011) par1(iparts), par2(iparts)
c 1011 format('alpha = ',f11.7,' b    = ',f11.7,'                    ')
************************************************************************
c      else if (idis(iparts) .eq. 2) then
c      s3='Size distribution : two parameter gamma                     '
c      write(s4,1012) par1(iparts), par2(iparts)
c 1012 format('reff  = ',f11.7,' veff  = ',f11.7,'                   ')
************************************************************************
c      else if (idis(iparts) .eq. 3) then
c      s3='Size distribution : bimodal gamma with equal mode weights   '
c      write(s4,1013) par1(iparts), par2(iparts), par3(iparts)
c 1013 format('reff1 = ',f11.7,' reff2 = ',f11.7,' veff = ',f11.7)
************************************************************************
c      else if (idis(iparts) .eq. 4) then
c      s3='Size distribution : log normal                              '
c      write(s4,1014) par1(iparts), par2(iparts)
c 1014 format('rg     = ',f11.7,' sigma = ',f11.7,'                    ')
************************************************************************
c      else if (idis(iparts) .eq. 5) then
c      s3='Size distribution : log normal                              '
c      write(s4,1015) par1(iparts), par2(iparts)
c 1015 format('reff   = ',f11.7,' veff  = ',f11.7,'                    ')
************************************************************************
c      else if (idis(iparts) .eq. 6) then
c      s3='Size distribution : power law                               '
c      write(s4,1016) par1(iparts), par2(iparts), par3(iparts)
c 1016 format('alpha = ',f11.7,' rmin  = ',f11.7,' rmax  = ',f11.7)
************************************************************************
c      else if (idis(iparts) .eq. 7) then
c      s3='Size distribution : modified gamma                          '
c      write(s4,1017) par1(iparts), par2(iparts), par3(iparts)
c 1017 format('alpha = ',f11.7,' rc    = ',f11.7,' gamma = ',f11.7)
************************************************************************
c      else if (idis(iparts) .eq. 8) then
c      s3='Size distribution : modified gamma                          '
c      write(s4,1018) par1(iparts), par2(iparts), par3(iparts)
c 1018 format('alpha = ',f11.7,' b     = ',f11.7,' gamma = ',f11.7)
************************************************************************
c      else if (idis(iparts) .eq. 9) then
c      s3='Size distribution : discrete values                         '
c      write(s4,1019) par1(iparts), par2(iparts), par3(iparts)
c 1019 format('min radius = ',f11.7,' max = ',f11.7,
c     +       ' ngaus = ',f5.1)
************************************************************************
c      else
c          write(7,*) ' MAIN: illegal size distribution nr:',idis(iparts)
c      endif
      if (nangle .ge. 1) then
          cosbar = coefs(1,1,1)/3.D0
      else
          cosbar = 0.D0
      endif
************************************************************************
*  Write the coefficients out onto file                                *
************************************************************************
c      call writsc(nrunit,coefs,NDcoef,ncoef,miec(5),cosbar
c     +        ,miec(4),miec(3),miec(6),miec(10),s1, s2, s3, s4, s5, s6 )
************************************************************************
*  Print a table of the scattering matrix and the linear polarization  *
************************************************************************
c 1020 write (7,911)
      do 30 j=1,nangle
          i     = nangle-j+1
          xp    = -100.d0*F(2,i)/F(1,i)
          theta1 = radtod*dacos(u(i))
c          write(7,913) j,theta1,F(1,i),F(2,i),F(3,i),F(4,i),xp
   30 continue
************************************************************************
*  Print the coefficients on output.                                   *
*  Evaluate the expansion in GSF at an equidistant set of scattering   *
*  angles and print the resulting scattering matrix.                   *
************************************************************************
      if (develop.eq.1) then
c          write (7,914)
          do 40 i=0,ncoef
c              write(7,915) i,coefs(1,1,i),coefs(2,2,i),coefs(3,3,i),
c     +                       coefs(4,4,i),coefs(1,2,i),coefs(3,4,i)
   40    continue
         call expand( ncoef, npunt, coefs, u, F, thmin, thmax, step )
c         write (7,911)
         do 45 i=1,npunt
             xp    = -100.d0*F(2,i)/F(1,i)
             theta1 = radtod*dacos(u(i))
c             write(7,913) i,theta1,F(1,i),F(2,i),F(3,i),F(4,i),xp
   45    continue
      endif
   50 continue
************************************************************************
*  End of loop over 'particles'.                                       *
************************************************************************
c      if (develop .ne. 0) call clossc(nrunit)

c  800 format(1x )
c  801 format(51x,i2)
c  802 format(51x,e8.2)
c  803 format(51x,i1)
c  804 format(51x,f5.2)
c  805 format(8x,3(2x,f7.4),4x,i1,4x,i2,4x,i4)
c  806 format(8x,5(2x,f7.2))

c  911 format(/,'    num   scat.angle     F11              F21   ',
c     +           '           F33             F34         pol',
c     +'arization',/,t90,'[proc]')
c  913 format(4x,i3,4x,f6.2,4(3x,e14.8),3x,f7.3)
c  914 format(1h1,13x,'alpha1',11x,'alpha2',11x,'alpha3',11x,
c     +'alpha4',12x,'beta1',12x,'beta2',/)
c  915 format(1h ,'l= ',i3,6f17.12)
c  920 format('1PARTICLE NUMBER ',i3,/)
c  921 format(1h ,'the number of computed scattering angles was :'
c     +,i7)
c  922 format(1h ,'Csca : the average scattering cross section  :'
c     +,e24.14)
c  923 format(1h ,'Cext : the average extinction cross section  :'
c     +,e24.14)
c  924 format(1h ,'Qsca : the scattering efficiency factor      :'
c     +,e24.14)
c  925 format(1h ,'Qext : the extinction efficiency factor      :'
c     +,e24.14)
c  926 format(1h ,'a    : the single scattering albedo          :'
c     +,e24.14)
c  927 format(1h ,'G    : the average geometrical cross section :'
c     +,e24.14)
c  928 format(1h ,'reff : the average effective radius          :'
c     +,e24.14)
c  929 format(1h ,'xeff : the average effective size parameter  :'
c     +,e24.14)
c  930 format(1h ,'num  : the integrated number of particles    :'
c     +,e24.14)
c  931 format(1h ,'V    : the average volume                    :'
c     +,e24.14)

c      stop 'program mie terminated'
	return
      end

      subroutine anbn( m, x, nmax, psi, chi, d, an, bn )
************************************************************************
*  Calculate the Mie coefficients an and bn.                           *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( NDn=10000000 )
      double complex m, zn, znm1, save, perm
      double complex an(nmax), bn(nmax), D(nmax)
      dimension psi(0:nmax), chi(0:nmax)
      perm = 1.D0/m
      perx = 1.D0/x
      xn   = 0.D0
* DEBUG
*     write(7,*) ' anbn:'
*     write(7,*) '     Re(an)           Im(an)       Re(bn)      Im(bn)'
* END DEBUG
      do 100 n=1, nmax
          zn   = dcmplx( psi(n),   chi(n))
          znm1 = dcmplx( psi(n-1), chi(n-1))
          xn   = dble(n)*perx
          save = D(n)*perm+xn
          an(n)= (save*psi(n)-psi(n-1)) / (save*zn-znm1)
          save = m*D(n)+xn
          bn(n)= (save*psi(n)-psi(n-1)) / (save*zn-znm1)
* DEBUG
*          write(72,*) n,an(n),bn(n),psi(n),chi(n),D(n)
* END DEBUG
  100 continue
      return
      end
      subroutine clossc(iunit)
************************************************************************
*  On entry :                                                          *
*      iunit     number of the unit to be closed                       *
*  On exit :                                                           *
*      The file is closed.                                             *
*                                                 VLD January 9, 1989  *
************************************************************************
c      close(unit=iunit,err=999)
      return
c  999 write(7,*) ' clossc: error in closing file with unit number',iunit
c      stop 'in clossc error in closing file'
      end
      subroutine devel(ncoef,nangle,u,w,F,coefs)
************************************************************************
*  Calculate the expansion coefficients of the scattering matrix in    *
*  generalized spherical functions by numerical integration over the   *
*  scattering angle.                                                   *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( NDcoef=3000, NDang = 6000 )
      dimension u(NDang),w(NDang),F(4,NDang)
      dimension coefs(4,4,0:NDcoef), P00(NDang,2), P02(NDang,2)
     +                             , P22(NDang,2), P2m2(NDang,2)
************************************************************************
*  Initialization                                                      *
************************************************************************
      qroot6 = -0.25D0*dsqrt(6.D0)
      do 20 j=0,NDcoef
          do 20 i=1,4
              do 20 ii=1,4
                  coefs(ii,i,j)=0.D0
   20         continue
   30     continue
   40 continue
************************************************************************
*  Multiply the scattering matrix F with the weights w for all angles  *
*  We do this here because otherwise it should be done for each l      *
************************************************************************
      do 60 k=1,4
          do 50 i=1, nangle
              F(k,i) = w(i)*F(k,i)
   50     continue
   60 continue
************************************************************************
*  Start loop over the coefficient index l                             *
*  first update generalized spherical functions, then calculate coefs. *
*  lold and lnew are pointer-like indices used in recurrence           *
************************************************************************
      lnew = 1
      lold = 2

      do 70 l=0, ncoef
               if (l .eq. 0) then
************************************************************************
*             Adding paper Eq. (77) with m=0                           *
************************************************************************
              do 80 i=1, nangle
                  P00(i,lold) = 1.D0
                  P00(i,lnew) = 0.D0
                  P02(i,lold) = 0.D0
                  P22(i,lold) = 0.D0
                  P2m2(i,lold)= 0.D0
                  P02(i,lnew) = 0.D0
                  P22(i,lnew) = 0.D0
                  P2m2(i,lnew)= 0.D0
   80         continue
          else
              fac1 = (2.D0*l-1.d0)/dble(l)
              fac2 = dble(l-1.d0)/dble(l)
************************************************************************
*             Adding paper Eq. (81) with m=0                           *
************************************************************************
              do 90 i=1, nangle
                  P00(i,lold) = fac1*u(i)*P00(i,lnew) - fac2*P00(i,lold)
   90         continue
          endif
          if (l .eq. 2) then
************************************************************************
*             Adding paper Eqs. (78) and (80)  with m=2                *
*             sql4 contains the factor dsqrt(l*l-4) needed in          *
*             the recurrence Eqs. (81) and (82)                        *
************************************************************************
              do 100 i=1, nangle
                  P02(i,lold) = qroot6*(1.D0-u(i)*u(i))
                  P22(i,lold) = 0.25D0*(1.D0+u(i))*(1.D0+u(i))
                  P2m2(i,lold)= 0.25D0*(1.D0-u(i))*(1.D0-u(i))
                  P02(i,lnew) = 0.D0
                  P22(i,lnew) = 0.D0
                  P2m2(i,lnew)= 0.D0
  100         continue
              sql41 = 0.D0
          else if (l .gt. 2) then
************************************************************************
*             Adding paper Eq. (82) with m=0 and m=2                   *
************************************************************************
              sql4  = sql41
              sql41 = dsqrt(dble(l*l)-4.d0)
              twol1 = 2.D0*dble(l)-1.d0
              tmp1  = twol1/sql41
              tmp2  = sql4/sql41
              denom = (dble(l)-1.d0)*(dble(l*l)-4.d0)
              fac1  = twol1*(dble(l)-1.d0)*dble(l)/denom
              fac2  = 4.D0*twol1/denom
              fac3  = dble(l)*((dble(l)-1.d0)*(dble(l)-1.d0)-4.d0)/denom
              do 110 i=1, nangle
                  P02(i,lold) = tmp1*u(i)*P02(i,lnew) - tmp2*P02(i,lold)
                  P22(i,lold) = (fac1*u(i)-fac2)*P22(i,lnew)
     +                                                - fac3*P22(i,lold)
                  P2m2(i,lold)= (fac1*u(i)+fac2)*P2m2(i,lnew)
     +                                               - fac3*P2m2(i,lold)
  110         continue
          endif
************************************************************************
*         Switch indices so that lnew indicates the function with      *
*         the present index value l, this mechanism prevents swapping  *
*         of entire arrays.                                            *
************************************************************************
          itmp = lnew
          lnew = lold
          lold = itmp
************************************************************************
*         Now calculate the coefficients by integration over angle     *
*         See de Haan et al. (1987) Eqs. (68)-(73).                    *
*         Remember for Mie scattering : F11 = F22 and F33 = F44        *
************************************************************************
          alfap = 0.D0
          alfam = 0.D0
          do 120 i=1, nangle
              coefs(1,1,l) = coefs(1,1,l) + P00(i,lnew)*F(1,i)
              alfap = alfap + P22(i,lnew)*(F(1,i)+F(3,i))
              alfam = alfam + P2m2(i,lnew)*(F(1,i)-F(3,i))
              coefs(4,4,l) = coefs(4,4,l) + P00(i,lnew)*F(3,i)
              coefs(1,2,l) = coefs(1,2,l) + P02(i,lnew)*F(2,i)
              coefs(3,4,l) = coefs(3,4,l) + P02(i,lnew)*F(4,i)
  120     continue
************************************************************************
*         Multiply with trivial factors like 0.5D0*(2*l+1)             *
************************************************************************
          fl = dble(l)+0.5D0
          coefs(1,1,l) =  fl*coefs(1,1,l)
          coefs(2,2,l) =  fl*0.5D0*(alfap+alfam)
          coefs(3,3,l) =  fl*0.5D0*(alfap-alfam)
          coefs(4,4,l) =  fl*coefs(4,4,l)
          coefs(1,2,l) =  fl*coefs(1,2,l)
          coefs(3,4,l) =  fl*coefs(3,4,l)
          coefs(2,1,l) =     coefs(1,2,l)
          coefs(4,3,l) =    -coefs(3,4,l)
   70 continue
************************************************************************
*     End of loop over index l                                         *
************************************************************************

************************************************************************
*     Remove the weight factor from the scattering matrix              *
************************************************************************
      do 140 k=1, 4
          do 130 i=1, nangle
              F(k,i) = F(k,i)/w(i)
  130     continue
  140 continue
      return
      end
      subroutine expand(ncoef,nangle,coefs,u,F,thmin,thmax,step)
c
c     ***********************************************************
c
      implicit double precision (a-h,o-z)
      parameter( pi=3.141592653589793238462643D0, radfac= pi/180.D0 )
      parameter( NDcoef=3000, NDang=6000 )
      dimension u(NDang),F(4,NDang)
      dimension coefs(4,4,0:NDcoef), P00(NDang,2), P02(NDang,2)
c
c     **********************************************************
c
c                           initialize

      do 10 i=1, nangle
          u(i) = dcos(radfac*(thmin+dble(i-1)*step))
   10 continue
      qroot6 = -0.25D0*dsqrt(6.D0)
************************************************************************
*  Set scattering matrix F to zero                                     *
************************************************************************
      do 30 k=1,4
          do 20 i=1, nangle
              F(k,i) = 0.D0
   20     continue
   30 continue
************************************************************************
*  Start loop over the coefficient index l                             *
*  first update generalized spherical functions, then calculate coefs. *
*  lold and lnew are pointer-like indices used in recurrence           *
************************************************************************
      lnew = 1
      lold = 2
      do 90 l=0, ncoef
               if (l .eq. 0) then
************************************************************************
*             Adding paper Eqs. (76) and (77) with m=0                 *
************************************************************************
              do 40 i=1, nangle
                  P00(i,lold) = 1.D0
                  P00(i,lnew) = 0.D0
                  P02(i,lold) = 0.D0
                  P02(i,lnew) = 0.D0
   40         continue
          else
              fac1 = (2.D0*dble(l)-1.d0)/dble(l)
              fac2 = (dble(l)-1.d0)/dble(l)
************************************************************************
*             Adding paper Eq. (81) with m=0                           *
************************************************************************
              do 50 i=1, nangle
                  P00(i,lold) = fac1*u(i)*P00(i,lnew) - fac2*P00(i,lold)
   50         continue
          endif
          if (l .eq. 2) then
************************************************************************
*             Adding paper Eq. (78)                                    *
*             sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in  *
*             the recurrence Eqs. (81) and (82)                        *
************************************************************************
              do 60 i=1, nangle
                  P02(i,lold) = qroot6*(1.D0-u(i)*u(i))
                  P02(i,lnew) = 0.D0
   60         continue
              sql41 = 0.D0
          else if (l .gt. 2) then
************************************************************************
*             Adding paper Eq. (82) with m=0                           *
************************************************************************
              sql4  = sql41
              sql41 = dsqrt(dble(l*l)-4.d0)
              tmp1  = (2.D0*dble(l)-1.d0)/sql41
              tmp2  = sql4/sql41
              do 70 i=1, nangle
                  P02(i,lold) = tmp1*u(i)*P02(i,lnew) - tmp2*P02(i,lold)
   70         continue
          endif
************************************************************************
*         Switch indices so that lnew indicates the function with      *
*         the present index value l, this mechanism prevents swapping  *
*         of entire arrays.                                            *
************************************************************************
          itmp = lnew
          lnew = lold
          lold = itmp
************************************************************************
*         Now add the l-th term to the scattering matrix.              *
*         See de Haan et al. (1987) Eqs. (68)-(73).                    *
*         Remember for Mie scattering : F11 = F22 and F33 = F44        *
************************************************************************
          do 80 i=1, nangle
              F(1,i) = F(1,i) + coefs(1,1,l)*P00(i,lnew)
              F(2,i) = F(2,i) + coefs(1,2,l)*P02(i,lnew)
              F(3,i) = F(3,i) + coefs(4,4,l)*P00(i,lnew)
              F(4,i) = F(4,i) + coefs(3,4,l)*P02(i,lnew)
   80     continue
   90 continue
************************************************************************
*     End of loop over index l                                         *
************************************************************************
      return
      end
      subroutine fichid( m, x, nchi, nmax, nD, psi, chi, D )
************************************************************************
*  Calculate functions psi(x)  chi(x) and D(z) where z = mx.           *
*  On entry, the following should be supplied :                        *
*      m      : complex index of refraction                            *
*      x      : sizeparameter                                          *
*      nchi   : starting index for backwards recurrence of chi         *
*      nmax   : number of chi, psi and D that must be available        *
*      nD     : starting index for backwards recurrence of D           *
*  On exit, the desired functions are returned through psi, chi and D  *
************************************************************************
      implicit double precision (a-h,o-z)
      double complex D,m,z, perz, zn1
      dimension psi(0:nchi), chi(0:nmax+1), D(nd)
*
      z = m*x
      perz = 1.D0/z
      perx = 1.D0/x
      sinx = dsin(x)
      cosx = dcos(x)
      psi(nchi)=0d0
************************************************************************
*  (mis-) use the array psi to calculate the functions rn(x)
*  De Rooij and van der Stap Eq. (A6)
************************************************************************
      do 10 n=nchi-1, 0, -1
          psi(n) = 1.D0 / (dble(2*n+1)/x - psi(n+1))
   10 continue
************************************************************************
*  Calculate functions D(z) by backward recurrence
*  De Rooij and van der Stap Eq. (A11)
************************************************************************
      D(nD) = dcmplx(0.D0,0.D0)
      do 20 n=nD - 1, 1, -1
          zn1 = dble(n+1)*perz
          D(n) = zn1 - 1.D0/(D(n+1)+zn1)
   20 continue
************************************************************************
*  De Rooij and van der Stap Eqs. (A3) and (A1)
*  De Rooij and van der Stap Eq. (A5) where psi(n) was equal to r(n)
*  and Eq. (A2)
************************************************************************
      psi(0) = sinx
      psi1   = psi(0)*perx - cosx
      if (dabs(psi1) .gt. 1.d-4) then
          psi(1) = psi1
          do 30 n=2,nmax
              psi(n) = psi(n)*psi(n-1)
   30     continue
      else
          do 35 n=1,nmax
              psi(n) = psi(n)*psi(n-1)
   35     continue
      endif
************************************************************************
*  De Rooij and van der Stap Eqs. (A4) and (A2)
************************************************************************
      chi(0) = cosx
      chi(1) = chi(0)*perx + sinx
      do 40 n=1, nmax-1
          chi(n+1) = dble(2*n+1)*chi(n)*perx - chi(n-1)
   40 continue
* DEBUG
*     write(7,*) ' fichid: x = ',x
*     write(7,12)
*     do 26 n=0, nchi
*         write(7,11) n,psi(n),chi(n),D(n)
*  26 continue
*  11 format(i4,4e24.14)
*  12 format('   n',t20,'psi(n)',t44,'chi(n)',t68
*    +      ,'Re(D(n))',t92,'Im(D(n))',/
*    +      ,' ----------------------------------------------------'
*    +      ,'-----------------------------------------------------')
* END DEBUG
      return
      end
      subroutine gaulegCP(ndim,ngauss,a,b,x,w)
************************************************************************
*   Given the lower and upper limits of integration a and b, and given *
*   the number of Gauss-Legendre points ngauss, this routine returns   *
*   through array x the abscissas and through array w the weights of   *
*   the Gauss-Legendre quadrature formula. Eps is the desired accuracy *
*   of the abscissas. This routine is documented further in :          *
*   W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)   *
*   page 125 ISBN 0-521-30811-9                                        *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( eps = 1.d-14 )
      double precision x(ndim), w(ndim)
      double precision a,b,xm,xl,z,p1,p2,p3,pp,z1,pi
      pi=4.D0*datan(1.D0)
      m=(ngauss+1)/2
      xm=0.5D0*(a+b)
      xl=0.5D0*(b-a)
      do 12 i=1,m
*         THIS IS A REALLY CLEVER ESTIMATE :
          z= dcos(pi*(dble(i)-0.25D0)/(dble(ngauss)+0.5D0))
    1     continue
              p1=1.D0
              p2=0.D0
              do 11 j=1,ngauss
                  p3= p2
                  p2= p1
                  p1=((dble(2*j)-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
   11         continue
              pp=ngauss*(z*p1-p2)/(z*z-1.d0)
              z1= z
              z= z1-p1/pp
          if (dabs(z-z1) .gt. eps) goto 1
          x(i)= xm-xl*z
          x(ngauss+1-i)= xm+xl*z
          w(i)=2.D0*xl/((1.D0-z*z)*pp*pp)
******          write(7,*) ' gaulegCP: ',i,'  x=',x(i),' w=',w(i)
          w(ngauss+1-i)= w(i)
   12 continue
      return
      end
      subroutine gausspt(ndim,ngauss,a,b,x,w)
c
c     ***********************************************************
c
c     put the gauss-legendre points of order ndim in array x,
c     the weights in array w. the points are transformed to the
c     interval [a,b]
c
c     ***********************************************************
c
      implicit double precision (a-h,o-z)
      dimension x(ndim),w(ndim)
c
c     ***********************************************************
c
c                     find starting values
c
      gn=0.5D0/dble(ngauss)
      extra=1.0D0/(.4D0*dble(ngauss)*dble(ngauss)+5.0D0)
      xz=-gn
      nt=0
      nteken=0
    5 pnm2=1.0D0
      pnm1= xz
      do 10 i=2,ngauss
          pnm1xz= pnm1*xz
          pn=2.0D0*pnm1xz-pnm2-(pnm1xz-pnm2)/dble(i)
          pnm2= pnm1
          pnm1= pn
   10 continue
      mteken=1
      if (pn .le. 0.0D0) mteken=-1
      if ((mteken+nteken) .eq. 0) then
          nt=nt+1
          x(nt)= xz
      endif
      nteken=mteken
      if ((1.0D0-xz) .le. extra) go to 30
      xz= xz+(1.D0-xz*xz)*gn+extra
      go to 5
   30 continue
c
c     ***********************************************************
c
c                determine zero's and weights
c
      do 60 i=1,nt
          xz= x(i)
          delta2=1.D0
   35         pnm2=1.0D0
              pnm1= xz
              pnm1af=1.0D0
              z=.5D0+1.5D0*xz*xz
              do 40 k=2,ngauss
                  pnm1xz= pnm1*xz
                  pn=2.0D0*pnm1xz-pnm2-(pnm1xz-pnm2)/dble(k)
                  pnaf= xz*pnm1af+k*pnm1
                  z= z+(dble(k)+0.5D0)*pn*pn
                  pnm2= pnm1
                  pnm1= pn
                  pnm1af= pnaf
   40         continue
              delta1= pn/pnaf
              xz= xz-delta1
              if(delta1.lt.0.0D0) delta1=-delta1
              if((delta1.ge.delta2).and.(delta2.lt.1.d-14)) go to 50
              delta2= delta1
          go to 35
   50     x(i)= xz
          w(i)=1.0D0/z
   60 continue
c
c     ***********************************************************
c
c                 transform to the interval [a,b]
c
      nghalf=ngauss/2
      ngp1=ngauss+1
      ntp1=nt+1
      apb= a+b
      bmag2=(b-a)/2.0D0
      do 90 i=1,nghalf
          x(ngp1-i)= b-bmag2*(1.0D0-x(ntp1-i))
   90     w(ngp1-i)= bmag2*w(ntp1-i)
      if (nghalf .ne. nt) then
          x(nt)= apb/2.0D0
          w(nt)= w(1)*bmag2
      endif
      do 120 i=1,nghalf
          x(i)= apb-x(ngp1-i)
          w(i)= w(ngp1-i)
  120 continue
      return
      end


      subroutine pitau(u,nmax,pi,tau)
c     ***********************************************************
c     calculates pi,n(u) and tau,n(u) with upward recursion
c
c     ***********************************************************
      implicit double precision (a-h,o-z)
      dimension pi(nmax),tau(nmax)
c
c       starting values:
      pi(1) = 1.D0
      pi(2) = 3.D0*u
      delta = 3.D0*u*u-1.d0
      tau(1)= u
      tau(2)= 2.D0*delta-1.d0

c       upward recursion:
      do 10 n=2, nmax-1
          pi(n+1) = dble(n+1)/dble(n)*delta + u*pi(n)
          delta   = u*pi(n+1) - pi(n)
          tau(n+1)= dble(n+1)*delta - pi(n)
   10 continue
      return
      end
      subroutine rminmax( idis, nsubr, ngaur, par1, par2, par3, cutoff
     +                  , rmin, rmax 
     +                        , iparts, ndis, rdis, nwrdis )
************************************************************************
*  Find the integration bounds rmin and rmax for the integration over  *
*  a size distribution. These bounds are chosen such that the size     *
*  distribution falls below the user specified cutoff. It is essential *
*  that the size distribution is normalized such that the integral     *
*  over all r is equal to one !                                        *
*  This is programmed rather clumsy and will in the future be changed  *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( eps = 1.d-10, NDpart = 4, NDdis=300 )
      double precision nwithr, rdis(NDpart,NDdis), nwrdis(NDpart,NDdis)
      dimension r(1), nwithr(1)
*
      if (idis.eq.0) then
         rmin= par1
         rmax= par1
      else
          goto (10,20,30,40,50,60,70,80,90) idis
          write(*,*) ' rminmax: illegal size distribution index :',idis
          stop 'in rminmax illegal size distribution index'
************************************************************************
   10     sef = 1.D0/dsqrt(par2+3.D0)
          ref = 1.D0/(sef*sef*par2)
          rref= ref
          goto 100

   20     ref = par1
          sef = dsqrt(par2)
          rref= ref
          goto 100

   30     sef = dsqrt(par3)
          ref = dmax1(par1,par2)+sef
          rref= dmin1(par1,par2)
          goto 100

   40     sef = dsqrt(dexp(dlog(par2)**2)-1.d0)
          ref = par1*(1.D0+sef*sef)**0.4D0
          rref= ref
          goto 100

   50     ref = par1
          sef = dsqrt(ref)
          rref= ref
          goto 100

   60     rmin= par2
          rmax= par3
          goto 999

   70     ref = par2
          sef = 2.D0*ref
          rref=0.5D0*ref
          goto 100

   80     ref = (par1/(par2*par3))**par3
          sef = 2.D0*ref
          rref= 0.5D0*ref
          goto 100

   90     rmin = par1
          rmax = par2
          goto 999

************************************************************************
*  search for a value of r such that the size distribution
*  is less than the cutoff. Start the search at ref+sef which          *
*  guarantees that such a value will be found on the TAIL of the       *
*  distribution.                                                       *
************************************************************************
  100     r(1) = ref+sef
          r0   = ref
  200          call sizedis( idis, par1, par2, par3, r, 1, nwithr 
     +                                , iparts, ndis, rdis, nwrdis )
          if (nwithr(1) .gt. cutoff) then
              r0   = r(1)
              r(1) = 2.D0*r(1)
              goto 200
          endif
          r1 = r(1)
************************************************************************
*  Now the size distribution assumes the cutoff value somewhere        *
*  between r0 and r1  Use bisection to find the corresponding r        *
************************************************************************
  300     r(1) = 0.5D0*(r0+r1)
          call sizedis( idis, par1, par2, par3, r, 1, nwithr 
     +                                , iparts, ndis, rdis, nwrdis )
          if (nwithr(1) .gt. cutoff) then
              r0 = r(1)
          else
              r1 = r(1)
          endif
          if ((r1-r0) .gt. eps) goto 300
          rmax = 0.5D0*(r0+r1)
************************************************************************
*  Search for a value of r on the low end of the size distribution     *
*  such that the distribution falls below the cutoff. There is no      *
*  guarantee that such a value exists, so use an extra test to see if  *
*  the search comes very near to r = 0                                 *
************************************************************************
          r1 = rref
          r0 = 0.D0
  400     r(1) = 0.5D0*r1
          call sizedis( idis, par1, par2, par3, r, 1, nwithr 
     +                                , iparts, ndis, rdis, nwrdis )
          if (nwithr(1) .gt. cutoff) then
              r1 = r(1)
              if (r1 .gt. eps) goto 400
          else
              r0 = r(1)
          endif
************************************************************************
*  Possibly the size distribution goes through cutoff between r0       *
*  and r1 try to find the exact value of r where this happens by       *
*  bisection.                                                          *
*  In case there is no solution, the algorithm will terminate soon.    *
************************************************************************
  500     r(1) = 0.5D0*(r0+r1)
          call sizedis( idis, par1, par2, par3, r, 1, nwithr 
     +                                , iparts, ndis, rdis, nwrdis )
          if (nwithr(1) .gt. cutoff) then
              r1 = r(1)
          else
              r0 = r(1)
          endif
          if ((r1-r0) .gt. eps) goto 500
          if(r1 .le. eps) then
              rmin = 0.D0
          else
              rmin = 0.5D0*(r0+r1)
          endif
      endif

  999 continue!write(7,*) ' rminmax: found rmin = ',rmin,' rmax = ',rmax
      return
      end

      subroutine scatmat( u, wth, m, lambda, idis
     +                  , thmin, thmax, step, develop
     +                  , nsub, ngauss, rmin, rmax
     +                  , par1, par2, par3
     +                  , delta, F, miec, nangle 
     +                  , ndis, rdis, nwrdis, iparts )
************************************************************************
*  Calculate the scattering matrix of an ensemble of homogenous        *
*  spheres. On entry, the following must be supplied :                 *
*     m            : complex index of refraction                       *
*     lambda       : wavelength                                        *
*     idis         : index of the size distribution                    *
*     nsub         : number of subintervals for integration over r     *
*     ngauss       : number of Gauss points used per subinterval       *
*     rmin         : lower bound for integration over r                *
*     rmax         : upper bound for integration over r                *
*     par1,2,3     : parameters of the size distribution               *
*     delta        : cutoff used in truncation of the Mie sum          *
*     thmin        : minimum scattering angle in degrees               *
*     thmax        : maximum scattering angle in degrees               *
*     step         : step in scattering angle in degrees               *
*     develop      : expansion in GSF (1) or not (0)                   *
*  On exit, the following results are returned :                       *
*     u            : cosines of scattering angles                      *
*     wth          : Gaussian weights associated with u                *
*     F            : scattering matrix for all cosines in u            *
*     miec         : array containing cross sections etc.              *
*     nangle       : the number of scattering angles                   *
*  When develop=1 u contains Gauss points, when develop=0 u contains   *
*  cosines of equidistant angles between thmin and thmax with step.    *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( NDn=10000000, NDr=1000, NDang=6000, NDdis=300, NDpart = 4)
      double complex   m, ci, Splusf, Sminf, cSplusf
      double complex   cSminf, Splusb, Sminb, cSplusb, cSminb

      double complex, allocatable :: an(:), bn(:),D(:)
      double precision,allocatable :: pi(:), tau(:), fi(:), chi(:)
      double precision,allocatable :: facf(:), facb(:)

      double precision lambda, nwithr, miec, numpar, thmin, thmax, step
     +               , nwrdis
      integer     develop
      dimension   u(NDang), wth(NDang), F(4,NDang),rdis(NDpart, NDdis)
      dimension   miec(13), nwrdis(NDpart, NDdis),r(NDr), w(NDr), nwithr(NDr)
      logical     symth
************************************************************************
*  Initialization                                                      *
************************************************************************
      do 10 j=1,NDang
          do 10 k=1,4
   10         F(k,j)=0.D0
      Csca  = 0.D0
      Cext  = 0.D0
      numpar= 0.D0
      G     = 0.D0
      reff  = 0.D0
      nfou  = 0
      fac90 = 1.D0
      ci    = dcmplx(0.D0,1.D0)
      call tstsym( develop, thmin, thmax, step, symth )
************************************************************************
*  Constants                                                           *
************************************************************************
      pie   = dacos(-1.d0)
      radfac= pie/180.D0
      rtox  = 2.D0*pie/lambda
      fakt  = lambda*lambda/(2.D0*pie)
* nfac is the number of precomputed factors (2n+1)/(n*(n+1))
      nfac  = 0
************************************************************************
*  distinguish between distribution or not                             *
************************************************************************
      if (idis.eq.0) then
         w(1)     = 1.D0
         r(1)     = rmin
         nwithr(1)= 1.D0
         nsub     = 1
         ngauss   = 1
         dr       = 0.D0
      else
         dr = (rmax-rmin)/dble(nsub)
*        call gausspt( ngauss, ngauss, (rmax-dr) , rmax ,r ,w )
         call gaulegCP(  ngauss, ngauss, (rmax-dr) , rmax ,r ,w )
         call sizedis( idis, par1, par2, par3, r, ngauss, nwithr 
     +               , iparts, ndis, rdis, nwrdis )
      endif
************************************************************************
*  Start integration over radius r with largest radius !               *
************************************************************************
      do 60 l=nsub,1,-1
      do 50 i=ngauss,1,-1
*
      sw   = nwithr(i)*w(i)
      x    = rtox*r(i)
      nmax = x + 4.05D0*x**(1.D0/3.D0) + 20
      nfi  = nmax+60
      zabs = x*cdabs(m)
      nD   = zabs + 4.05D0*zabs**(1.D0/3.D0) + 70

	allocate(an(max(nD,nfi,nmax)))
	allocate(bn(max(nD,nfi,nmax)))
	allocate(pi(max(nD,nfi,nmax)))
	allocate(tau(max(nD,nfi,nmax)))
	allocate(fi(0:max(nD,nfi,nmax)))
	allocate(chi(0:max(nD,nfi,nmax)))
	allocate(D(max(nD,nfi,nmax)))
	allocate(facf(max(nD,nfi,nmax)))
	allocate(facb(max(nD,nfi,nmax)))

      if ((nD.gt.NDn) .or. (nfi.gt.NDn)) then
          write(*,*) ' scatmat: estimated number of Mie-terms:',nD
          write(*,*) '          for particle sizeparameter   :',x
          write(*,*) '          maximum NDn is only          : ',NDn
          stop 'in scatmat no room for all Mie terms'
      endif
      call fichid( m, x, nfi, nmax, nD, fi, chi, D )
      call anbn( m, x, nmax, fi, chi, D, an, bn )
************************************************************************
*  Precompute the factor (2n+1)/(n*(n+1)) needed in Mie sum over n     *
************************************************************************
      if (nmax .gt. nfac) then
          do 26 n=nfac+1, nmax
              facf(n) = dble(2*n+1)/dble(n*(n+1))
              facb(n) = facf(n)
              if (mod(n,2) .eq. 1) facb(n) = -facb(n)
   26     continue
          nfac = nmax
      endif
************************************************************************
*  Calculate extinction and scattering cross section                   *
*  Use the convergence criterion to determine the number of terms that *
*  will later be used in the mie sum for the scattering matrix itself  *
************************************************************************
      Cextsum = 0.D0
      Cscasum = 0.D0
      nstop = nmax
      do 52 n=1, nmax
          aux = (2.D0*dble(n)+1.D0) *
     +             dabs(dble(an(n)*conjg(an(n)) + bn(n)*conjg(bn(n))))
          Cscasum = Cscasum + aux
          Cextsum = Cextsum + (2.D0*n+1.D0)*dble(an(n)+bn(n))
          if (aux .lt. delta) then
              nstop = n
              goto 53
          endif
   52 continue
   53 nfou = nstop
      if (nfou .ge. nmax) then
          write(*,*) ' WARNING from scatmat : Mie sum not converged for'
     +           ,' scattering cross section'
          write(*,*) '   radius r = ',r(i),' sizeparameter x = ',x
     +           ,' sizedistribution nr. ',idis
          write(*,*) '   Re(m) = ',dble(m),' Im(m) = ',dimag(m)
          write(*,*) '   a priori estimate of number of Mie terms:',nmax
          write(*,*) '   term ',nmax,' for Csca was ',aux
          write(*,*) '   should have been less than ',delta
          write(*,*) '   the a priori estimate will be used'
      endif
************************************************************************
*  Only for the first run through the loop set points in u= dcos(th)   *
************************************************************************
      if ((l.eq.nsub) .and. (i.eq.ngauss)) then
************************************************************************
*  In case of expansion in GSF : set Gauss points for dcos(th)         *
************************************************************************
          if (develop .eq. 1) then
************************************************************************
*  Ensure accurate integrations: add two terms: nangle = 2*nfou+2      *
*  One should be sufficient, but total should be even !                *
************************************************************************
              nangle = 2*nfou+2
              if (nangle .gt. NDang) then
                 write(*,*) ' scatmat: need too many integration angles'
     +                  ,' nangle=',nangle
                 write(*,*) '       maximum array size NDang = ',NDang
                 stop 'too many integration angles'
              endif
*             call gausspt(nangle,nangle,-1.d0,1.D0,u,wth)
              call gaulegCP(nangle,nangle,-1.d0,1.D0,u,wth)
          else
************************************************************************
*  In case no expansion in GSF is desired : set u= dcos(th) for        *
*  for equidistant angles between thmin and thmax.                     *
************************************************************************
              nangle = idnint((thmax-thmin)/step) + 1
              if (nangle .gt. NDang) then
                 write(*,*) ' scatmat : need too many scattering angles'
     +                  ,' nangle=',nangle
                 write(*,*) '       maximum array size NDang = ',NDang
                 stop 'too many scattering angles'
              endif
              wfac = 2.d0/dble(nangle)
              do 260 iang=1, nangle
                  th = thmin + dble(iang-1)*step
                  u(nangle+1-iang) = dcos( radfac*th )
                  wth(iang) = wfac
  260         continue
          endif
      endif
************************************************************************
*  Integration for normalization of size distibution, geometrical      *
*  cross section and effective radius                                  *
************************************************************************
      numpar = numpar+sw
      G      = G     +sw*r(i)*r(i)
      reff   = reff  +sw*r(i)*r(i)*r(i)
      if (symth) then
************************************************************************
*  Start loop over scattering angles, do only half and use symmetry    *
*  between forward and backward scattering angles                      *
*  The factor fac90 will later be used to correct for the fact that    *
*  for a symmetrical set of scattering angles with an odd number of    *
*  angles the scattering matrix is a factor 2 too big at 90 degrees    *
*  due to the way we programmed the symmetry relations                 *
************************************************************************
        if (mod(nangle,2) .eq. 1) then
            nhalf = (nangle+1)/2
            fac90 = 0.5D0
        else
            nhalf = nangle/2
        endif
*
        do 40 j=1, nhalf
            call pitau( u(j), nmax, pi, tau )
            Splusf = dcmplx(0.D0,0.D0)
            Sminf  = dcmplx(0.D0,0.D0)
            Splusb = dcmplx(0.D0,0.D0)
            Sminb  = dcmplx(0.D0,0.D0)
*  THIS IS THE INNERMOST LOOP !! (Mie sum)
*  can be programmed more efficiently by taking the facf multiplication
*  outside the angle loop over index j 
            do 20 n=1,nfou
                Splusf = Splusf + facf(n)*(an(n)+bn(n)) * (pi(n)+tau(n))
                Sminf  = Sminf  + facf(n)*(an(n)-bn(n)) * (pi(n)-tau(n))
                Splusb = Splusb + facb(n)*(an(n)+bn(n)) * (pi(n)-tau(n))
                Sminb  = Sminb  + facb(n)*(an(n)-bn(n)) * (pi(n)+tau(n))
   20       continue
            cSplusf = conjg(Splusf)
            cSminf  = conjg(Sminf )
            cSplusb = conjg(Splusb)
            cSminb  = conjg(Sminb )
            k = nangle-j+1
*  the forward scattering elements
            F(1,j) = F(1,j) +    sw*(Splusf*cSplusf + Sminf *cSminf)
            F(2,j) = F(2,j) -    sw*(Sminf *cSplusf + Splusf*cSminf)
            F(3,j) = F(3,j) +    sw*(Splusf*cSplusf - Sminf *cSminf)
            F(4,j) = F(4,j) + ci*sw*(Sminf *cSplusf - Splusf*cSminf)
*  the backward scattering elements
            F(1,k) = F(1,k) +    sw*(Splusb*cSplusb + Sminb *cSminb)
            F(2,k) = F(2,k) -    sw*(Sminb *cSplusb + Splusb*cSminb)
            F(3,k) = F(3,k) +    sw*(Splusb*cSplusb - Sminb *cSminb)
            F(4,k) = F(4,k) + ci*sw*(Sminb *cSplusb - Splusb*cSminb)
   40   continue
      else
************************************************************************
*  Start loop over scattering angles, do all angles                    *
************************************************************************
        do 400 j=1, nangle
            call pitau( u(j), nmax, pi, tau )
            Splusf = dcmplx(0.D0,0.D0)
            Sminf  = dcmplx(0.D0,0.D0)
*  THIS IS THE INNERMOST LOOP !! (Mie sum)
            do 200 n=1,nfou
                Splusf = Splusf + facf(n)*(an(n)+bn(n)) * (pi(n)+tau(n))
                Sminf  = Sminf  + facf(n)*(an(n)-bn(n)) * (pi(n)-tau(n))
  200       continue
            cSplusf = conjg(Splusf)
            cSminf  = conjg(Sminf )
            k = nangle-j+1
*  the forward scattering elements
            F(1,j) = F(1,j) +      sw*(Splusf*cSplusf + Sminf *cSminf)
            F(2,j) = F(2,j) -      sw*(Sminf *cSplusf + Splusf*cSminf)
            F(3,j) = F(3,j) +      sw*(Splusf*cSplusf - Sminf *cSminf)
            F(4,j) = F(4,j) + ci*sw*(Sminf *cSplusf - Splusf*cSminf)
  400   continue
      endif
************************************************************************
*  Integration for cross sections, shift radius to next subinterval    *
************************************************************************
      Csca = Csca + sw*Cscasum
      Cext = Cext + sw*Cextsum
      r(i) = r(i) - dr
   50 continue
      if (l .ne. 1)
     +         call sizedis( idis, par1, par2, par3, r, ngauss, nwithr 
     +                     , iparts, ndis, rdis, nwrdis )
   60 continue
************************************************************************
*  End of integration over size distribution                           *
*  Some final corrections :                                            *
************************************************************************
      do 70 j=1,nangle
          do 70 k=1,4
   70         F(k,j)= F(k,j)/(2.D0*Csca)
      if (symth) then
          do 80 k=1,4
   80         F(k,nhalf) = fac90*F(k,nhalf)
      endif
      G     = pie*G
      Csca  = Csca*fakt
      Cext  = Cext*fakt
      Qsca  = Csca/G
      Qext  = Cext/G
      albedo= Csca/Cext
      volume= (4.d0/3.d0)*pie*reff
      reff  = pie*reff/G
      xeff  = rtox*reff
*
      miec(1) = Csca
      miec(2) = Cext
      miec(3) = Qsca
      miec(4) = Qext
      miec(5) = albedo
      miec(6) = G
      miec(7) = reff
      miec(8) = xeff
      miec(9) = numpar
      miec(10)= volume
*

	deallocate(an)
	deallocate(bn)
	deallocate(pi)
	deallocate(tau)
	deallocate(fi)
	deallocate(chi)
	deallocate(D)
	deallocate(facf)
	deallocate(facb)

      return
      end
      subroutine sizedis( idis, par1, par2, par3, r, numr, nwithr
     +                  , iparts, ndis, rdis, nwrdis )
************************************************************************
*  Calculate the size distribution n(r) for the numr radius values     *
*  contained in array r and return the results through the array nwithr*
*  The size distributions are normalized such that the integral over   *
*  all r is equal to one.                                              *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter ( NDdis=300, NDpart = 4 )
      double precision nwithr, logC, logC1, logC2, nwrdis, nwrdisp
      dimension r(numr),nwithr(numr),nwrdis(NDpart,NDdis)
     +        , nwrdisp(NDdis), rdis(NDpart,NDdis), rdisp(NDdis)
     +        , y2(NDdis)
*
      pi     = dacos(-1.d0)
      root2p = dsqrt(pi+pi)
      if (idis .eq. 0) return
      goto (10, 20, 30, 40, 50, 60, 70, 80, 90 ) idis
      write(*,*) ' sizedis: illegal index : ',idis
      stop 'illegal index in sizedis'
************************************************************************
*  1 TWO PARAMETER GAMMA with alpha and b given                        *
************************************************************************
   10 alpha = par1
      b     = par2
      alpha1= alpha+1.D0
      logC  = alpha1*dlog(b)-gammlnCP(alpha1)
      do 11 i=1, numr
          nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
   11 continue
      goto 999
************************************************************************
*  2 TWO PARAMETER GAMMA with par1= reff and par2= veff given          *
************************************************************************
   20 alpha = 1.D0/par2 - 3.D0
      b     = 1.D0/(par1*par2)
      alpha1= alpha+1.D0
      logC  = alpha1*dlog(b)-gammlnCP(alpha1)
      do 21 i=1, numr
          nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
   21 continue
      goto 999
************************************************************************
*  3 BIMODAL GAMMA with equal mode weights                             *
************************************************************************
   30 alpha = 1.D0/par3 - 3.D0
      b1    = 1.D0/(par1*par3)
      b2    = 1.D0/(par2*par3)
      gamlna= gammlnCP(alpha+1.D0)
      logC1 = (alpha+1.D0)*dlog(b1)-gamlna
      logC2 = (alpha+1.D0)*dlog(b2)-gamlna
      do 31 i=1, numr
          nwithr(i) = 0.5D0*( dexp(logC1+alpha*dlog(r(i))-b1*r(i))
     +                      + dexp(logC2+alpha*dlog(r(i))-b2*r(i)) )
   31 continue
      goto 999
************************************************************************
*  4 LOG NORMAL with rg and sigma given                                *
************************************************************************
   40 flogrg = dlog(par1)
      flogsi = dabs(dlog(par2))
      C      = 1.D0/(root2p*flogsi)
      fac    = -0.5D0/(flogsi*flogsi)
      do 41 i=1, numr
          nwithr(i) = C * dexp( fac*(dlog(r(i))-flogrg)**2 ) / r(i)
   41 continue
      goto 999
************************************************************************
*  5 LOG NORMAL with reff and veff given                               *
************************************************************************
   50 rg     = par1/(1.D0+par2)**2.5D0
      flogrg = dlog(rg)
      flogsi = dsqrt(dlog(1.D0+par2))
      C      = 1.D0/(root2p*flogsi)
      fac    = -0.5D0/(flogsi*flogsi)
      do 51 i=1, numr
          nwithr(i) = C * dexp( fac*(dlog(r(i))-flogrg)**2 ) / r(i)
   51 continue
      goto 999
************************************************************************
*  6 POWER LAW                                                         *
************************************************************************
   60 alpha = par1
      rmin  = par2
      rmax  = par3
* DEBUG
*     write(7,*) ' sizedis : power law with alpha = ',alpha
*     write(7,*) '                          rmin  = ',rmin
*     write(7,*) '                          rmax  = ',rmax
* END DEBUG
      if (dabs(alpha+1.D0) .lt. 1.d-10) then
          C = 1.D0/dlog(rmax/rmin)
      else
          alpha1 = alpha-1.d0
          C = alpha1 * rmax**alpha1 / ((rmax/rmin)**alpha1-1.d0)
      endif
      do 61 i=1, numr
          if ((r(i) .lt. rmax) .and. (r(i) .gt. rmin)) then
              nwithr(i) = C*r(i)**(-alpha)
          else
              nwithr(i) = 0.D0
          endif
   61 continue
      goto 999
************************************************************************
*  7 MODIFIED GAMMA with alpha, rc and gamma given                     *
************************************************************************
   70 alpha = par1
      rc    = par2
      gamma = par3
      b     = alpha / (gamma*rc**gamma)
      aperg = (alpha+1.D0)/gamma
      logC  = dlog(gamma) + aperg*dlog(b) - gammlnCP(aperg)
      do 71 i=1, numr
          nwithr(i) = dexp( logC + alpha*dlog(r(i)) - b*r(i)**gamma )
   71 continue
      goto 999
************************************************************************
*  8 MODIFIED GAMMA with alpha, b and gamma given                      *
************************************************************************
   80 alpha = par1
      b     = par2
      gamma = par3
      aperg = (alpha+1.D0)/gamma
      logC  = dlog(gamma) + aperg*dlog(b) - gammlnCP(aperg)
      do 81 i=1, numr
          nwithr(i) = dexp( logC + alpha*dlog(r(i)) - b*r(i)**gamma )
   81 continue
      goto 999
************************************************************************
*  9 DISCRETE VALUED DISTRIBUTION                                      *
************************************************************************
   90 do 95 k=1, ndis
          rdisp(k)   = rdis(iparts,k)
          nwrdisp(k) = nwrdis(iparts,k)
   95 continue
      call splineCP(rdisp,nwrdisp,ndis,1.d30,1.d30,y2,NDdis,NDdis)
      do 91 j=1,numr
          call splintCP(rdisp,nwrdisp,y2,ndis,r(j),nwithr(j),NDdis,NDdis)
   91 continue
      goto 999
************************************************************************
  999 if (numr .le. 1) return
* DEBUG
*     write(7,*) ' sizedis:'
*     write(7,*) '   i             r(i)               n(r(i))'
*     write(7,*) ' ----------------------------------------------------'
*     do 1000 i=1, numr
*         write(7,1001) i,r(i),nwithr(i)
*1000 continue
*1001 format(i4,1pe24.14,1pe22.14)
* END DEBUG
      return
      end
 
 
      function gammlnCP(xarg)
************************************************************************
*  Return the value of the natural logarithm of the gamma function.    *
*  The argument xarg must be real and positive.                        *
*  This function is documented in :                                    *
*                                                                      *
*  W.H. Press et al. 1986, 'Numerical Recipes' Cambridge Univ. Pr.     *
*  page 157 (ISBN 0-521-30811)                                         *
*                                                                      *
*  When the argument xarg is between zero and one, the relation (6.1.4)*
*  on page 156 of the book by Press is used.                           *
*                                         V.L. Dolman April 18 1989    *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( eps = 1.d-7, one = 1.D0, two = 2.D0, half = 0.5D0
     +         , fpf = 5.5D0 )
      dimension cof(6)
      data cof,stp/76.18009173D0,-86.50532033D0, 24.01409822D0
     +   ,-1.231739516D0, 0.120858003D-2, -0.536382D-5, 2.50662827465D0/
      pi = 4.D0*datan(1.D0)
      if (xarg .le. 0.D0) then
        write(*,*) ' gammlnCP: called with negative argument xarg = ',xarg
        stop 'function gammlnCP called with negative value'
      endif
      if (dabs(xarg-one) .lt. eps) then
          write(*,*) ' gammlnCP: argument too close to one for algorithm'
          stop ' in function gammlnCP argument too close to one'
      endif
      if (xarg .ge. one) then
          xx = xarg
      else
          xx = xarg+two
      endif
      x = xx-one
      tmp = x+fpf
      tmp = (x+half)*dlog(tmp)-tmp
      ser = one
      do 11 j=1, 6
          x = x+one
          ser = ser+cof(j)/x
   11 continue
      gtmp = tmp+dlog(stp*ser)
      if  (xarg .gt. one) then
          gammlnCP = gtmp
      else
          pix = pi*(one-xarg)
          gammlnCP = dlog(pix/dsin(pix))-gtmp
      endif
      return
      end
      

      subroutine splineCP(x, y, n, yp1, ypn, y2,NDmui,nn)
***********************************************************************
** Spline interpolation routine from Press et al. (1986, p.88).      **
**                                                                   **
** Given arrays x and y of length n containing a tabulated function, **
** i.e. y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1  **
** and ypn for the first derivative of the interpolating function at **
** points 1 and n respectively, this routine returns an array y2 of  **
** length n which contains the second derivatives of the interpola-  **
** ting function at the tabulated points x(i).                       **
** If yp1 and/or yp2 are equal to 1x10^30 or larger, the routine is  **
** signalled to set the corresponding boundary condition for a natu- **
** ral spline, with zero second derivative on that boundary.         **
***********************************************************************
      implicit double precision (a-h,o-z)
      parameter (nmax=1500)
      dimension x(NDmui),y(nn),y2(nn),u(nmax)

      if (nn .ne. NDmui) then
          write(*,*) ' splineCP : uncomfortable error occurred, '
          write(*,*) ' nn .ne. NDmui '
          stop ' dimension error in splineCP '
      endif
      if (nmax .lt. nn) then
          write(*,*) ' splineCP : nmax may not be smaller than nn '
          write(*,*) ' nmax = ',nmax,' nn = ',nn
          stop ' dimension error in splineCP '
      endif
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     +       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
   11 continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
   12 continue
      return
      end
      subroutine splintCP(xa, ya, y2a, n, x, y,NDmui,nn)
***********************************************************************
** Spline interpolation routine from Press et al. (1986, p.88).      **
**                                                                   **
** Given the arrays xa and ya of length n, which tabulate a function **
** (with the xa(i)'s in order), and given the array y2a, which is    **
** the output from splineCP above, and given a value of x, this        **
** routine returns a cubic-spline interpolated value y.              **
***********************************************************************
      implicit double precision (a-h,o-z)

      dimension xa(NDmui),ya(nn),y2a(nn)

      if (nn .ne. NDmui) then
          write(*,*) ' a very stupid an unnecessary error occurred '
          write(*,*) ' in splintCP: nn .ne. NDmui '
          stop ' dimension error in splintCP '
      endif
      klo=1
      khi=n

    1 if (khi-klo.gt.1) then
        k=(klo+khi)/2.d0  
        if (xa(k).gt.x) then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (dabs(h).lt.1.d-10) write (7,10) 
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     +  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
   10 format('ERROR in routine splintCP: Bad XA input.')
      end
      subroutine tstsym( develop, thmin, thmax, step, symth )
************************************************************************
*  Test if the set of theta points is symmetrical around 90 degrees    *
*  and return this information through logical symth                   *
*  In case of development in GSF we have a symmetrical Gauss set !!    *
************************************************************************
      implicit double precision(a-h,o-z)
      parameter( eps=1.d-6, heps=0.5d0*eps )
      integer develop
      logical symth
      symth = .false.
      if (develop .eq. 1) then
          symth = .true.
      else
          if ( (dabs( 180.d0 - thmin - thmax ) .lt. eps)  .and.
     +         (dmod( thmax-thmin+heps, step ) .lt. eps) )
     +              symth = .true.
      endif
* DEBUG
c      if (symth) then
c          write(7,*) ' tstsym: theta points symmetrical'
c      else
c          write(7,*) ' tstsym: theta points NOT symmetrical !!'
c      endif
* END DEBUG
      return
      end



