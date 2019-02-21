	subroutine MCDoComputeT(converged,f)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iphase
	real*8 z,dz,E,tau,Planck,random,v,dx,dy,f
	real*8,allocatable :: Ca(:,:,:),Cs(:,:),Ce(:,:,:),g(:,:),spec(:,:,:),EJv_omp(:),EJv2_omp(:)
	real*8,allocatable :: specsource_star(:),specsource_planet(:),EJv_phot(:)
	real*8 tau_v,x,y,EJv(nr),E0_star,E0_planet,EJv2(nr),dEJv(nr)
	real*8 tot,Lum,tot2,Cplanck(nr),dTdP_ad,dTdP,must,Eprev(nr),rr,cwg(ng),dT0(nr)
	real*8 CabsL(nlam),T0(nr),E0,chi2,nabla,scale,Crw(nr),Cpl(nr),dlnT,dlnP,kapmax
	integer iphot,ir,Nphot,ilam,ig,nscat,jrnext,NphotStar,NphotPlanet,jr,jr0
	integer iT1,iT2,iT,i,icount,iopenmp0,omp_get_thread_num,iter,niter
	logical docloud0(max(nclouds,1)),goingup,onedge,converged,dorw(nr),RandomWalkT
	type(Mueller),allocatable :: M(:,:)
	
	allocate(specsource_star(nlam-1))
	allocate(specsource_planet(nlam-1))
	allocate(spec(nlam-1,ng,nr))
	allocate(Ca(nr,nlam,ng))
	allocate(Ce(nr,nlam,ng))
	allocate(Cs(nr,nlam))
	allocate(g(nr,nlam))
	allocate(M(nr,nlam))

	niter=10

	NphotPlanet=max(Nphot0/10/niter,10)
	NphotStar=max(Nphot0/niter,100)

	T0=T

	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo

	call InitRandomWalk()

	must=betaT
	if(i2d.ne.0) then
		if(i2d.eq.1) then
			call ComputeBeta(90d0,twind,must)
		else if(i2d.eq.2) then
			call ComputeBeta(270d0,twind,must)
		else if(i2d.eq.3) then
			call ComputeBeta(0d0,twind,must)
		else if(i2d.eq.4) then
			call ComputeBeta(180d0,twind,must)
		endif
		must=must*betaT
	endif

	do ilam=1,nlam
		do ig=1,ng
			do ir=1,nr
				call Crossections(ir,ilam,ig,Ca(ir,ilam,ig),Cs(ir,ilam),docloud0)
				Ce(ir,ilam,ig)=Ca(ir,ilam,ig)+Cs(ir,ilam)
			enddo
		enddo
	enddo
	do ir=1,nr
		do ilam=1,nlam-1
			call GetMatrix(ir,ilam,M(ir,ilam),docloud0)
			g(ir,ilam)=0d0
			tot=0d0
			do iphase=1,180
				g(ir,ilam)=g(ir,ilam)+M(ir,ilam)%F11(iphase)*costheta(iphase)*sintheta(iphase)
				tot=tot+M(ir,ilam)%F11(iphase)*sintheta(iphase)
			enddo
			g(ir,ilam)=g(ir,ilam)/tot
		enddo
	enddo

	EJv=0d0
	EJv2=0d0

	do iter=1,niter

	tau=0d0
	do ir=nr,1,-1
		iT=T0(ir)+1
		iT=min(max(iT,1),nBB)
		tot=0d0
		tot2=0d0
		Crw(ir)=0d0
		Cpl(ir)=0d0
		do ilam=1,nlam-1
			if(.not.Cs(ir,ilam).gt.0d0) then
				Cs(ir,ilam)=0d0
			endif
			if(.not.g(ir,ilam).gt.0d0) then
				g(ir,ilam)=0d0
			endif
			do ig=1,ng
				if(.not.Ca(ir,ilam,ig).gt.0d0) then
					if(ig.eq.1) then
						Ca(ir,ilam,ig)=1d-30
					else
						Ca(ir,ilam,ig)=Ca(ir,ilam,ig-1)
					endif
				endif
				Cpl(ir)=Cpl(ir)+wgg(ig)*dfreq(ilam)*BB(iT,ilam)*Ca(ir,ilam,ig)
				Crw(ir)=Crw(ir)+wgg(ig)*dfreq(ilam)*BB(iT,ilam)/(Ca(ir,ilam,ig)+Cs(ir,ilam)*(1d0-g(ir,ilam)))
			enddo
			tot=tot+dfreq(ilam)*BB(iT,ilam)
			tot2=tot2+dfreq(ilam)*BB(iT,ilam)
		enddo
		if(iter.eq.1) Eprev(ir)=Cpl(ir)*abs(R(ir+1)-R(ir))
		Cpl(ir)=Cpl(ir)/tot
		Crw(ir)=tot2/Crw(ir)

		kapmax=(2d0/7d0)*((T(ir)/TeffP)**4)*16d0/(3d0*Hp(ir))
		if(Crw(ir).gt.kapmax) then
			scale=kapmax/Crw(ir)
			Ca(ir,1:nlam,1:ng)=Ca(ir,1:nlam,1:ng)*scale
			Cs(ir,1:nlam)=Cs(ir,1:nlam)*scale
			Ce(ir,1:nlam,1:ng)=Ce(ir,1:nlam,1:ng)*scale
			Crw(ir)=Crw(ir)*scale
			Cpl(ir)=Cpl(ir)*scale
		endif
		tau=tau+Crw(ir)*abs(R(ir+1)-R(ir))
		dorw(ir)=(tau.gt.factRW)
	enddo

	do ir=1,nr
		iT=max(1,min(int(T0(ir)),nBB))
		tot=0d0
		do ilam=1,nlam-1
			do ig=1,ng
				spec(ilam,ig,ir)=Ca(ir,ilam,ig)*BB(iT,ilam)
				tot=tot+wgg(ig)*dfreq(ilam)*spec(ilam,ig,ir)
			enddo
		enddo
		spec(1:nlam-1,1:ng,ir)=spec(1:nlam-1,1:ng,ir)/tot
	enddo


	Nphot=0
	if(TeffP.gt.0d0) then
		Nphot=Nphot+NphotPlanet
		tot=0d0
		do ilam=1,nlam-1
			specsource_planet(ilam)=Planck(TeffP,freq(ilam))
			tot=tot+dfreq(ilam)*specsource_planet(ilam)
		enddo
		specsource_planet=specsource_planet/tot
		Lum=((2d0*(pi*kb*TeffP)**4)/(15d0*hplanck**3*clight**3))
		E0_planet=Lum/real(NphotPlanet)
	endif
	if(must.gt.0d0) then
		Nphot=Nphot+NphotStar
		Lum=0d0
		tot=0d0
		do ilam=1,nlam-1
			specsource_star(ilam)=Fstar(ilam)/(4d0*pi*Dplanet**2)
			Lum=Lum+dfreq(ilam)*specsource_star(ilam)
			tot=tot+dfreq(ilam)*Planck(Tstar,freq(ilam))
		enddo
		specsource_star=specsource_star/Lum
		scale=((2d0*(pi*kb*Tstar)**4)/(15d0*hplanck**3*clight**3))/tot
		Lum=Lum*scale
		E0_star=Lum*max(must,0d0)/real(NphotStar)
	endif

	if(Nphot.eq.0) then
		T=2.7d0
		return
	endif

	do ir=1,nr
		call cummulR(spec(1:nlam-1,1:ng,ir),(nlam-1)*ng)
	enddo
	call cummul(specsource_star,nlam-1)
	call cummul(specsource_planet,nlam-1)

	cwg(1)=wgg(1)
	do ig=2,ng
		cwg(ig)=wgg(ig)+cwg(ig-1)
	enddo
	cwg(1:ng)=cwg(1:ng)/cwg(ng)

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(iphot,x,y,z,dx,dy,dz,goingup,E0,jr,ilam,ig,onedge,EJv_omp,icount,rr,jr0,EJv_phot,EJv2_omp)
!$OMP& SHARED(Nphot,NphotPlanet,TeffP,specsource_planet,specsource_star,must,ng,nlam,nr,Crw,Cpl,Ce,Ca,Cs,g,M,EJv,spec,
!$OMP& 	dorw,E0_planet,E0_star,R,iopenmp0,niter,iter,cwg,EJv2)
	allocate(EJv_omp(nr))
	EJv_omp=0d0
	allocate(EJv2_omp(nr))
	EJv2_omp=0d0
	allocate(EJv_phot(nr))
!$OMP DO
	do iphot=1,Nphot
		EJv_phot=0d0
		if(iphot.le.NphotPlanet.and.TeffP.gt.0d0) then
			call emit(specsource_planet,ilam)
			rr=random(idum)
			call hunt(cwg,ng,rr,ig)
			ig=ig+1
			ig=min(max(1,ig),ng)
			x=0d0
			y=0d0
			z=R(1)
			call randomdirection(dx,dy,dz)
			dz=sign(dz,1d0)
			goingup=.true.
			onedge=.true.
			jr=1
			E0=E0_planet
		else
			call emit(specsource_star,ilam)
			rr=random(idum)
			call hunt(cwg,ng,rr,ig)
			ig=ig+1
			ig=min(max(1,ig),ng)
			x=0d0
			y=0d0
			z=R(nr+1)
			dz=-max(min(must,1d0),0d0)
			dx=sqrt(1d0-dz**2)
			dy=0d0
			goingup=.false.
			onedge=.true.
			jr=nr
			E0=E0_star
		endif
		do while((jr.le.nr.and.jr.ge.1).and.random(idum).gt.1d-8)
1			continue
			if((dorw(jr).and.jr.lt.nr).and.random(idum).gt.1d-8) then
				if(RandomWalkT(x,y,z,dx,dy,dz,E0,Crw,Cpl,jr,EJv_phot,dorw)) goto 1
			endif
			call travelcomputeT(x,y,z,dx,dy,dz,E0,jr,onedge,goingup,
     &			Ce,Ca,Cs,g,M,EJv_phot,ilam,ig)
			do jr=1,nr
				if(z.ge.R(jr).and.z.lt.R(jr+1)) exit
			enddo
			if(jr.le.nr.and.jr.ge.1) then
				call emitR(spec,ilam,ig,jr)
				call randomdirection(dx,dy,dz)
				goingup=(dz.gt.0d0)
				onedge=.false.
			endif
		enddo
		EJv_omp=EJv_omp+EJv_phot
		EJv2_omp=EJv2_omp+(EJv_phot)**2
	enddo
!$OMP END DO
!$OMP CRITICAL
	EJv(1:nr)=EJv(1:nr)+EJv_omp(1:nr)
	EJv2(1:nr)=EJv2(1:nr)+EJv2_omp(1:nr)
!$OMP END CRITICAL
	deallocate(EJv_omp)
!$OMP FLUSH
!$OMP END PARALLEL

	chi2=0d0
	do ir=nr,1,-1
		CabsL=0d0
		do ilam=1,nlam-1
			do ig=1,ng
				CabsL(ilam)=CabsL(ilam)+Ca(ir,ilam,ig)*wgg(ig)
			enddo
		enddo
		E0=EJv(ir)/real(niter)+Eprev(ir)*(1d0-real(iter)/real(niter))
		call increaseT(T(ir),T0(ir),E0,CabsL,ir)
	enddo

	do ir=nr-1,1,-1
		dlnP=log(P(ir+1)/P(ir))
		dlnT=log(T0(ir+1)/T0(ir))
		if((dlnT/dlnP).gt.nabla_ad(ir)) then
			dlnT=(nabla_ad(ir))*dlnP
			T0(ir)=T0(ir+1)/exp(dlnT)
		endif
	enddo

	call output_erase(trim(dbl2string(1d0*real(100*iter)/real(niter),'(f5.1)')) // " %")

	enddo

	call output("")

	do ir=nr,1,-1
		dEJv(ir)=sqrt(EJv2(ir))
		if((EJv(ir)**2)/EJv2(ir).gt.2d0) then
			dT0(ir)=0.25d0*T0(ir)*dEJv(ir)/EJv(ir)
		else
			dT0(ir)=1d5
			if(ir.lt.nr) then
				T0(ir)=T0(ir+1)
			endif
		endif
	enddo
	call Smooth(T0,dT0,T0,nr,idum)

	converged=.true.
	do ir=1,nr
		if(abs(T(ir)-T0(ir))/(T(ir)+T0(ir)).gt.epsiter) converged=.false.
	enddo

	chi2=0d0
	do ir=1,nr
		chi2=chi2+((min(T(ir),2900d0)-min(T0(ir),2900d0))/((min(T0(ir),2900d0)+min(T(ir),2900d0))*epsiter))**2
c		T(ir)=T(ir)**(1d0-f)*T0(ir)**f
		T(ir)=T(ir)*(1d0-f)+T0(ir)*f
	enddo
	chi2=chi2/real(nr)

	call output("Chi2: " // trim(dbl2string(chi2,'(f7.2)')))

	deallocate(specsource_star)
	deallocate(specsource_planet)
	deallocate(spec)
	deallocate(Ca)
	deallocate(Ce)
	deallocate(Cs)
	deallocate(g)
	deallocate(M)

	call WriteStructure

c	do ir=1,nr
c		if(T(ir).gt.2900d0) T(ir)=2900d0
c	enddo
	
	return
	end

	subroutine ConvectionTransport(x,y,z,ir)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 lr,x,y,z,rho1,rho2,random
	integer ir

	rho1=sqrt(x**2+y**2+z**2)
	lr=-log(random(idum))*Hp(ir)
	lr=random(idum)*Hp(ir)
	lr=lr/rho1
	x=x*(1d0+lr)
	y=y*(1d0+lr)
	z=z*(1d0+lr)
	rho2=sqrt(x**2+y**2+z**2)
	do ir=1,nr
		if(rho2.lt.R(ir+1).and.rho2.gt.R(ir)) exit
	enddo
	
	return
	end	


	subroutine increaseT(T1,T2,E,CabsL,ir)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ilam,ir,iTmax,iTmin,iT
	real*8 T1,T2,E,E0,CabsL(nlam),Ca0,Cs0,V

	V=(R(ir+1)-R(ir))

	iTmin=1

	iTmax=nBB

	iT=0.5d0*(iTmax+iTmin)

	do while(abs(iTmax-iTmin).gt.1)
		E0=0d0
		do ilam=1,nlam-1
			E0=E0+dfreq(ilam)*CabsL(ilam)*BB(iT,ilam)
		enddo
		E0=E0*V
		if(E0.gt.E) then
			iTmax=iT
		else
			iTmin=iT
		endif
		iT=0.5d0*(iTmax+iTmin)
		if(iT.lt.1) iT=1
		if(iT.gt.nBB) iT=nBB
	enddo
	T2=real(iT)
	
	return
	end
	

	subroutine cummul(spec,n)
	use GlobalSetup
	IMPLICIT NONE
	integer i,n,ig,ilam
	real*8 spec(n)

	spec(1)=spec(1)*dfreq(1)
	do i=2,n
		ilam=i
		spec(i)=spec(i-1)+dfreq(ilam)*spec(i)
	enddo
	spec=spec/spec(n)

	return
	end
	
	subroutine cummulR(spec,n)
	use GlobalSetup
	IMPLICIT NONE
	integer i,n,ig,ilam
	real*8 spec(n)

	spec(1)=spec(1)*dfreq(1)
	do i=2,n
		ig=(i-1)/(nlam-1)+1
		ilam=i-(nlam-1)*(ig-1)
		spec(i)=spec(i-1)+dfreq(ilam)*wgg(ig)*spec(i)
	enddo
	spec=spec/spec(n)

	return
	end
	

	subroutine emit(spec,ilam)
	use GlobalSetup
	IMPLICIT NONE
	integer ilam
	real*8 spec(nlam),Ltot,Lr
	real*8 random

	Lr=random(idum)

	call hunt(spec,nlam-1,Lr,ilam)
	ilam=max(1,ilam)
	
	return
	end
	


	subroutine emitR(spec,ilam,ig,jr)
	use GlobalSetup
	IMPLICIT NONE
	integer ilam,ig,jr,i,n
	real*8 spec(nlam-1,ng,nr),Ltot,Lr
	real*8 random

	Lr=random(idum)

	n=(nlam-1)*ng
	call hunt(spec(1:nlam-1,1:ng,jr),n,Lr,i)
	i=max(1,i)

	ig=(i-1)/(nlam-1)+1
	ilam=i-(nlam-1)*(ig-1)

	return
	end
	




	subroutine travelcomputeT(x,y,z,dx,dy,dz,E0,jr,onedge,goingup,Ce,Ca,Cs,g,M,EJv,ilam,ig)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,y,z,dx,dy,dz,Ce(nr,nlam,ng),Ca(nr,nlam,ng),Cs(nr,nlam),g(nr,nlam),EJv(nr),E0
	integer jr,nscat,jrnext,nscat0,ilam,ig
	logical onedge,goingup
	real*8 tau,v,tau_v,albedo,random,znext
	type(Mueller) M(nr,nlam)

1	continue
	tau=-log(random(idum))
2	continue
	goingup=(dz.gt.0d0)
	if(goingup) then
		v=abs((R(jr+1)-z)/dz)
		jrnext=jr+1
		znext=R(jrnext)
	else
		v=abs((R(jr)-z)/dz)
		jrnext=jr-1
		znext=R(jr)
	endif

	tau_v=v*Ce(jr,ilam,ig)
	albedo=Cs(jr,ilam)/Ce(jr,ilam,ig)
	if(tau_v.lt.tau) then
		x=x+v*dx
		y=y+v*dy
		z=z+v*dz
		tau=tau-tau_v
		EJv(jr)=EJv(jr)+E0*tau_v*(1d0-albedo)
		jr=jrnext
		z=znext
		if(jr.gt.nr) return
		if(jr.lt.1) then
			jr=1
			dz=sign(dz,1d0)
			z=R(1)
			goingup=.true.
			nscat=nscat+1
		endif
		onedge=.true.
		goto 2
	endif
	v=tau/Ce(jr,ilam,ig)
	x=x+v*dx
	y=y+v*dy
	z=z+v*dz
	EJv(jr)=EJv(jr)+E0*tau*(1d0-albedo)

	if(random(idum).gt.albedo) then
		return
	else
		call scattangle(M(jr,ilam),dx,dy,dz)
	endif
	onedge=.false.
	goto 1

	return
	end
	



	logical function RandomWalkT(x,y,z,dx,dy,dz,E0,Crw,Cpl,jr,EJv,dorw)
	use GlobalSetup
	use Constants
	use RandomWalkModule
	IMPLICIT NONE
	real*8 dmin,v,ry,random,lr,x,y,z,dx,dy,dz,Crw(nr),Cpl(nr),EJv(nr),E0,d1,d2,d,tau,znext,EE(nr)
	real*8 tauprev
	integer jr,iy,jr1,jr2,jrnext,djr
	logical dorw(nr)

	RandomWalkT=.false.

	d1=abs(z-R(jr))*Crw(jr)
	djr=-1
	jr1=jr+djr
	if(jr1.lt.1) then
		jr1=1
		djr=1
	endif
	do while(dorw(jr1).and.jr1.lt.nr)
		d1=d1+abs(R(jr1+1)-R(jr1))*Crw(jr1)
		jr1=jr1+djr
		if(jr1.lt.1) then
			jr1=1
			djr=1
		endif
	enddo
	d2=abs(z-R(jr+1))*Crw(jr)
	djr=1
	jr2=jr+djr
	do while(dorw(jr2).and.jr2.lt.nr)
		d2=d2+abs(R(jr2+1)-R(jr2))*Crw(jr2)
		jr2=jr2+djr
	enddo

	dmin=min(d1,d2)

	if(dmin.le.factRW) return

	RandomWalkT=.true.

	ry=random(idum)
	iy=1
	call hunt(phi,NY,ry,iy)

	v=-3d0*log(yy(iy))*dmin**2/(pi**2)

	call randomdirection(dx,dy,dz)

	EE=0d0
	tauprev=0d0

1	continue
	if(dz.gt.0d0) then
		d=(R(jr+1)-z)/dz
		tau=Crw(jr)*d
		jrnext=jr+1
		znext=R(jrnext)
	else
		d=(R(jr)-z)/dz
		tau=Crw(jr)*d
		jrnext=jr-1
		znext=R(jr)
	endif
	if(tau.lt.dmin) then
		x=x+d*dx
		y=y+d*dy
		z=z+d*dz
		EE(jr)=EE(jr)+((tauprev+tau)**2-tauprev**2)
		tauprev=tauprev+tau
		dmin=dmin-tau
		jr=jrnext
		z=znext
		if(jr.lt.1) then
			jr=1
			dz=abs(dz)
			z=R(1)
		endif
		if(jr.gt.nr) goto 2
		goto 1
	endif
	d=dmin/Crw(jr)
	x=x+d*dx
	y=y+d*dy
	z=z+d*dz
	EE(jr)=EE(jr)+((tauprev+dmin)**2-tauprev**2)

2	continue

	EE=EE/sum(EE)
	EJv=EJv+EE*E0*v*Cpl/Crw
	
	return
	end

	
