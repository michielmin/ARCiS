	module modComputeT
	IMPLICIT NONE
	real*8,allocatable :: CrV_prev(:),CrT_prev(:)
	real*8,allocatable :: Taverage(:)
	integer iaverage
	end module modComputeT
	
	subroutine DoComputeT(converged,f)
	use GlobalSetup
	use Constants
	use modComputeT
	use CloudModule
	IMPLICIT NONE
	integer iphase,iter
	real*8 tau_V,tau_T,Planck,CrV(nr),CrT(nr),CrVe(nr),f
	real*8 Cs,g,dlnT,dlnP
	real*8,allocatable :: Ce(:,:,:),Ca(:,:,:)
	real*8 tot,tot2,tot3,tot4,chi2,must,gamma,dP,Tirr,T0(nr),must_i
	integer ir,ilam,ig,i,iT,niter
	logical docloud0(max(nclouds,1)),converged
	type(Mueller) M	

c	call MCDoComputeT(converged)
c	converged=.false.
c	return
	
	if(.not.allocated(CrV_prev)) allocate(CrV_prev(nr),CrT_prev(nr),Taverage(nr))

	allocate(Ce(nr,nlam,ng))
	allocate(Ca(nr,nlam,ng))

	if(nTiter.eq.maxiter.and.iaverage.gt.0) then
		T=Taverage/real(iaverage)
		call WriteStructure
		return
	endif

	if(nTiter.eq.maxiter/2) then
		iaverage=0
		Taverage=0d0
	endif

	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo

	T0(1:nr)=T(1:nr)

	do ir=1,nr
		do ilam=1,nlam-1
			call GetMatrix(ir,ilam,M,docloud0)
			g=0d0
			tot=0d0
			do iphase=1,180
				g=g+M%F11(iphase)*costheta(iphase)*sintheta(iphase)
				tot=tot+M%F11(iphase)*sintheta(iphase)
			enddo
			g=g/tot
			do ig=1,ng
				call Crossections(ir,ilam,ig,Ca(ir,ilam,ig),Cs,docloud0)
				Ca(ir,ilam,ig)=Ca(ir,ilam,ig)/dens(ir)
				Ce(ir,ilam,ig)=(Ca(ir,ilam,ig)+Cs*(1d0-g))/dens(ir)
			enddo
		enddo
	enddo

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

	do iter=1,5

	do ir=1,nr
		iT=T0(ir)+1
		if(iT.gt.nBB-1) iT=nBB-1
		CrT(ir)=0d0
		tot2=0d0
		CrV(ir)=0d0
		tot3=0d0
		CrVe(ir)=0d0
		tot4=0d0
		do ilam=1,nlam-2
			do ig=1,ng
				CrT(ir)=CrT(ir)+dfreq(ilam)*BB(iT,ilam)/Ca(ir,ilam,ig)
				tot2=tot2+dfreq(ilam)*BB(iT,ilam)

				CrV(ir)=CrV(ir)+dfreq(ilam)*Fstar(ilam)/Ca(ir,ilam,ig)
				tot3=tot3+dfreq(ilam)*Fstar(ilam)

				CrVe(ir)=CrVe(ir)+dfreq(ilam)*Fstar(ilam)/Ce(ir,ilam,ig)
				tot4=tot4+dfreq(ilam)*Fstar(ilam)
			enddo
		enddo
		CrT(ir)=tot2/CrT(ir)
		CrV(ir)=tot3/CrV(ir)
		CrVe(ir)=tot4/CrVe(ir)
	enddo

	if(nTiter.ne.0) then
		CrV=CrV**f*CrV_prev**(1d0-f)
		CrT=CrT**f*CrT_prev**(1d0-f)
	endif

	chi2=0d0

	Tirr=sqrt(Rstar/Dplanet)*Tstar
	tau_V=0d0
	tau_T=0d0
	do ir=nr,1,-1
		if(ir.eq.nr) then
			dP=P(ir)
		else
			dP=abs(P(ir+1)-P(ir))
		endif
		dP=dP/100d0
		T0(ir)=0d0
		gamma=CrV(ir)/CrT(ir)
		do i=1,100
			tau_V=tau_V+CrV(ir)*1d6*dP/grav(ir)
			tau_T=tau_T+CrT(ir)*1d6*dP/grav(ir)
			gamma=tau_V/tau_T
			T0(ir)=T0(ir)+3d0*TeffP**4*(2d0/3d0+tau_T)/4d0+
     &	3d0*Tirr**4*(must*2d0/3d0+must**2/gamma+(gamma/3d0-must**2/gamma)*exp(-tau_V/must))/4d0
		enddo
		T0(ir)=T0(ir)/100d0
		T0(ir)=T0(ir)**0.25

		if(.not.T0(ir).gt.3d0) T0(ir)=3d0

		if(ir.lt.nr) then
			dlnP=log(P(ir+1)/P(ir))
			dlnT=log(T0(ir+1)/T0(ir))
			if((dlnT/dlnP).gt.nabla_ad(ir)) then
				dlnT=(nabla_ad(ir))*dlnP
				T0(ir)=T0(ir+1)/exp(dlnT)
			endif
		endif
	enddo
	
	enddo

	open(unit=35,file=trim(outputdir) // "/RosselandOpacity.dat",RECL=1000)
	tau_V=0d0
	tau_T=0d0
	do ir=nr,1,-1
		if(ir.eq.nr) then
			dP=P(ir)
		else
			dP=abs(P(ir+1)-P(ir))
		endif
		tau_V=tau_V+CrV(ir)*1d6*dP/grav(ir)
		tau_T=tau_T+CrT(ir)*1d6*dP/grav(ir)
		write(35,*) P(ir),CrV(ir),CrT(ir),tau_V,tau_T
	enddo
	close(unit=35)

	converged=.true.
	do ir=1,nr
		if(abs(T(ir)-T0(ir))/(T(ir)+T0(ir)).gt.epsiter) converged=.false.
	enddo

	chi2=0d0
	do ir=1,nr
		chi2=chi2+((min(T(ir),2900d0)-min(T0(ir),2900d0))/((min(T0(ir),2900d0)+min(T(ir),2900d0))*epsiter))**2
		T(ir)=T(ir)**(1d0-f)*T0(ir)**f
	enddo
	chi2=chi2/real(nr)

	CrV_prev=CrV
	CrT_prev=CrT

	if(nTiter.ge.maxiter/2) then
		iaverage=iaverage+1
		Taverage=Taverage+T0
	endif

	call WriteStructure

	deallocate(Ce)
	deallocate(Ca)

	return
	end




























	subroutine DoComputeTfluxes(converged,niter,f)
	use GlobalSetup
	use Constants
	use CloudModule
	IMPLICIT NONE
	integer iphase,iter
	real*8 Planck,Pb(nr+1),f,dtau,expint,cp
	real*8,allocatable :: Ca(:,:,:),Cs(:,:),g(:,:),Tr(:,:,:)
	real*8,allocatable :: Fl(:,:),Ff(:,:),Pl(:,:),T0(:),PlStar(:),PlEffP(:)
	real*8 tot,tot2,tot3,must,tstep,deltaF,err,chi2,dlnT,dlnP
	real*8 ComputeSolidAbun,Tmaxi,Tmini,f1,f0,deltaT,xx
	integer ir,ilam,ig,i,niter
	logical docloud0(max(nclouds,1)),converged
	type(Mueller) M
	
	allocate(Tr(nr,nlam,ng))
	allocate(T0(nr))
	allocate(Fl(2,nr+1))
	allocate(Ca(nr,nlam,ng))
	allocate(Cs(nr,nlam))
	allocate(g(nr,nlam))
	allocate(Ff(2,nr+1))
	allocate(Pl(nr,nlam))
	allocate(PlStar(nlam))
	allocate(PlEffP(nlam))

	if(nTiter.eq.0) call MCDoComputeT(converged)

	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo

	T0=T

	do ilam=1,nlam-1
		do ig=1,ng
			do ir=1,nr
				call Crossections(ir,ilam,ig,Ca(ir,ilam,ig),Cs(ir,ilam),docloud0)
				Ca(ir,ilam,ig)=Ca(ir,ilam,ig)/dens(ir)
				Cs(ir,ilam)=Cs(ir,ilam)/dens(ir)
			enddo
		enddo
	enddo

	Pb(1)=P(1)
	do i=2,nr
		Pb(i)=sqrt(P(i-1)*P(i))
	enddo
	Pb(nr+1)=P(nr)

	do ir=1,nr
		do ilam=1,nlam-1
			call GetMatrix(ir,ilam,M,docloud0)
			g(ir,ilam)=0d0
			tot=0d0
			do iphase=1,180
				g(ir,ilam)=g(ir,ilam)+M%F11(iphase)*costheta(iphase)*sintheta(iphase)
				tot=tot+M%F11(iphase)*sintheta(iphase)
			enddo
			g(ir,ilam)=g(ir,ilam)/tot
		enddo
	enddo

!$OMP PARALLEL IF(nlam.gt.200)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,ir,ig,dtau,xx,tot)
!$OMP& SHARED(nlam,nr,ng,Pb,Ca,Cs,g,grav,Tr)
!$OMP DO
	do ilam=1,nlam-1
		do ig=1,ng
			do ir=1,nr
				dtau=abs(Pb(ir+1)-Pb(ir))*1d6*(Ca(ir,ilam,ig)+Cs(ir,ilam)*(1d0-g(ir,ilam)))/grav(ir)
				if(dtau.gt.1d-2) then
					Tr(ir,ilam,ig)=(1d0-dtau)*exp(-dtau)+dtau**2*expint(1,dtau)
				else
					Tr(ir,ilam,ig)=0d0
					do i=1,100
						xx=(real(i)-0.5)/real(100)
						Tr(ir,ilam,ig)=Tr(ir,ilam,ig)+xx*exp(-dtau/xx)
					enddo
					Tr(ir,ilam,ig)=Tr(ir,ilam,ig)*2d0/real(100)
				endif
			enddo
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	
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

	tot=0d0
	do ilam=1,nlam-1
		PlStar(ilam)=Planck(Tstar,freq(ilam))
		tot=tot+PlStar(ilam)*dfreq(ilam)
	enddo
	tot=tot*pi*clight
	if(tot.eq.0d0) then
		PlStar=0d0
	else
		PlStar(1:nlam-1)=PlStar(1:nlam-1)*sigma*Tstar**4/tot
		PlStar(1:nlam-1)=PlStar(1:nlam-1)*must*(Rstar/Dplanet)**2
	endif
	tot=0d0
	do ilam=1,nlam-1
		PlEffP(ilam)=Planck(max(TeffP,3d0),freq(ilam))
		tot=tot+PlEffP(ilam)*dfreq(ilam)
	enddo
	tot=tot*pi*clight
	if(tot.eq.0d0) then
		PlEffP=0d0
	else
		PlEffP(1:nlam-1)=PlEffP(1:nlam-1)*sigma*max(TeffP,3d0)**4/tot
	endif
	
	do iter=1,niter

	if(niter.gt.500) call tellertje(iter,niter)
	Fl=0d0

	do ir=1,nr
		tot=0d0
		do ilam=1,nlam-1
			Pl(ir,ilam)=Planck(T(ir),freq(ilam))
			tot=tot+Pl(ir,ilam)*dfreq(ilam)
		enddo
		tot=tot*pi*clight
		if(tot.eq.0d0) then
 			Pl(ir,1:nlam-1)=0d0
 		else
 			Pl(ir,1:nlam-1)=Pl(ir,1:nlam-1)*sigma*T(ir)**4/tot
		endif
	enddo

	Ff=0d0
	do ilam=1,nlam-1
		do ig=1,ng
			Fl(1,nr+1)=pi*PlStar(ilam)
			do ir=nr,1,-1
				Fl(1,ir)=Tr(ir,ilam,ig)*Fl(1,ir+1)+pi*Pl(ir,ilam)*(1d0-Tr(ir,ilam,ig))
			enddo
			Fl(2,1)=Fl(1,1)+pi*PlEffP(ilam)
			do ir=2,nr+1
				Fl(2,ir)=Tr(ir-1,ilam,ig)*Fl(2,ir-1)+pi*Pl(ir-1,ilam)*(1d0-Tr(ir-1,ilam,ig))
			enddo

			do ir=1,nr+1
				Ff(1,ir)=Ff(1,ir)+dfreq(ilam)*Fl(1,ir)
				Ff(2,ir)=Ff(2,ir)+dfreq(ilam)*Fl(2,ir)
			enddo
		enddo
	enddo
	Ff=Ff/real(ng)
	
	converged=.true.
	do ir=nr,1,-1
		cp=3.5*kb/(mp*MMW(ir))
		tstep=cp*P(ir)/(sigma*grav(ir)*T(ir)**3)
		deltaF=((Ff(2,ir+1)-Ff(1,ir+1))-(Ff(2,ir)-Ff(1,ir)))/abs(R(ir+1)-R(ir))
		tstep=tstep*1d5/(abs(deltaF)**0.9)
		deltaT=deltaF*tstep/(dens(ir)*cp)
c		if(abs(deltaT).gt.(T(ir)*0.1d0)) deltaT=min(deltaT,deltaT*min(100d0,T(ir))*0.1d0/abs(deltaT))
		T(ir)=T(ir)-deltaT
		if(.not.T(ir).gt.3d0) T(ir)=3d0

		if(ir.lt.nr) then
			dlnP=log(P(ir+1)/P(ir))
			dlnT=log(T(ir+1)/T(ir))
			if((dlnT/dlnP).gt.nabla_ad(ir)) then
				dlnT=(nabla_ad(ir))*dlnP
				T(ir)=T(ir+1)/exp(dlnT)
			endif
		endif
	enddo

	do ir=1,nr
		err=clight*pi*((Ff(2,ir+1)-Ff(1,ir+1))-(Ff(2,ir)-Ff(1,ir)))/(sigma*T(ir)**4)
		if(abs(err).gt.1d-5) converged=.false.
	enddo
	if(converged) exit

	enddo

	converged=.true.
	do ir=1,nr
		if(abs(T(ir)-T0(ir))/(T(ir)+T0(ir)).gt.epsiter) converged=.false.
	enddo

	if(.not.allocated(ATP).or.nTiter.eq.0) then
		do ir=1,nr
			T(ir)=T0(ir)**(1d0-f)*T(ir)**f
		enddo
	else
		do ir=1,nr
			f0=ComputeSolidAbun(T0(ir),ir)
			f1=ComputeSolidAbun(T(ir),ir)
			if(f0.lt.0d0.and.f1.lt.0d0) then
				T(ir)=T0(ir)**(1d0-f)*T(ir)**f
			else
				f0=(f0+f1)/2d0
				Tmaxi=max(T(ir),T0(ir))
				Tmini=min(T(ir),T0(ir))
				T(ir)=(Tmaxi+Tmini)/2d0
				do iter=1,1000
					f1=ComputeSolidAbun(T(ir),ir)
					if(f1.lt.f0) then
						Tmaxi=T(ir)
						T(ir)=(T(ir)+Tmini)/2d0
					else
						Tmini=T(ir)
						T(ir)=(T(ir)+Tmaxi)/2d0
					endif
				enddo
			endif
		enddo
	endif

	deallocate(Tr)
	deallocate(T0)
	deallocate(Fl)
	deallocate(Ca)
	deallocate(Cs)
	deallocate(g)
	deallocate(Ff)
	deallocate(Pl)
	deallocate(PlStar)
	deallocate(PlEffP)

	call WriteStructure

	return
	end


	real*8 function ComputeSolidAbun(T0,ir)
	use GlobalSetup
	use Constants
	use CloudModule
	use AtomsModule
	IMPLICIT NONE
	real*8 T0,densv(nCS),fmax,f
	integer ir,iCS,iCS0
	
	densv=10d0**(BTP-ATP/T0-log10(T0))
	fmax=-1d200
	iCS0=1
	do iCS=1,nCS
		if(xv_bot(iCS).gt.0d0) then
			if(T0.gt.maxT(iCS)) densv(iCS)=densv(iCS)+10d0**(BTP(iCS)-ATP(iCS)/(T0*10d0)-log10(T0*10d0))
			f=xv_bot(iCS)*(1d0-densv(iCS)/(xv_bot(iCS)*dens(ir)))
			if(f.gt.fmax) then
				iCS0=iCS
				fmax=f
			endif
		endif
	enddo

	ComputeSolidAbun=fmax

	return
	end




	subroutine ComputeTevap()
	use GlobalSetup
	use Constants
	use CloudModule
	use AtomsModule
	IMPLICIT NONE
	real*8 dens0(nCS),densv(nCS),Tmini(nCS),Tmaxi(nCS),T0(nCS)
	integer ir,iCS,iter
	character*500 form
	
	open(unit=36,file=trim(outputdir) // "/Tevap.dat",RECL=1000)
	form='("#",a18,' // trim(int2string(nCS,'(i4)')) // 'a23,a19)'
	write(36,form) "P[bar]",CSname(1:nCS),"T[K]"
	form='(es19.7E3,' // trim(int2string(nCS,'(i4)')) // 'es23.7E3,es19.7E3)'
	
	do ir=1,nnr
		Tmini=3d0
		Tmaxi=3d5
		dens0=xv_bot*Clouddens(ir)
		do iter=1,1000
			T0=(Tmini+Tmaxi)/2d0
			densv=10d0**(BTP-ATP/T0-log10(T0))
			do iCS=1,nCS
				if(densv(iCS).gt.dens0(iCS)) then
					Tmaxi(iCS)=T0(iCS)
				else
					Tmini(iCS)=T0(iCS)
				endif
			enddo
		enddo
		write(36,form) CloudP(ir),T0(1:nCS)
	enddo

	close(unit=36)

	return
	end






































	subroutine MCDoComputeT(converged)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iphase
	real*8 z,dz,E,tau,Planck,random,v,dx,dy
	real*8,allocatable :: Ca(:,:,:),Cs(:,:),Ce(:,:,:),g(:,:),spec(:,:,:)
	real*8,allocatable :: specsource(:,:),Eplanck(:,:)
	real*8 vR1,vR2,b,rr,R1,R2,tau_v,x,y,theta,ct1,ct2,Tc(nr),Ec(nr),EJv(nr),EJvTot(nr)
	real*8 tot,Lum,tot2,mu0,gamma,Cplanck(nr),Cstar(nr),dTdP_ad,dTdP,albedo,Tmini,Tmaxi
	real*8 CabsL(nlam),Ca0,Cs0,T0,E0,Eabs,chi2,CabsLG(nlam,ng),Crw(nr),nabla,rho,scale,Pl
	integer iphot,ir,Nphot,ilam,ig,nscat,jrnext,NphotStar,NphotPlanet,jr,ir0,jr0
	integer iT1,iT2,iT,i
	logical docloud0(max(nclouds,1)),goingup,onedge,dorw(nr),converged,RandomWalkT
	external RandomWalkT
	type(Mueller),allocatable :: M(:,:)
	
	allocate(specsource(nlam,1))
	allocate(spec(nlam,ng,nr))
	allocate(Ca(nr,nlam,ng))
	allocate(Ce(nr,nlam,ng))
	allocate(Cs(nr,nlam))
	allocate(g(nr,nlam))
	allocate(M(nr,nlam))
	allocate(Eplanck(nr,nBB))
	
	Eplanck=-1d0
	
	call output("Temperature computation (in beta phase!!)")

	NphotPlanet=2500
	NphotStar=10000

	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo

	call InitRandomWalk()

	ir0=nr
	do ilam=1,nlam
		do ig=1,ng
			do ir=1,nr
				call Crossections(ir,ilam,ig,Ca(ir,ilam,ig),Cs(ir,ilam),docloud0)
				Ce(ir,ilam,ig)=Ca(ir,ilam,ig)+Cs(ir,ilam)
			enddo
		enddo
	enddo

	Tc=2.7d0
	Ec=0d0
	EJvTot=0d0

	jr0=0

	do ir=1,nr
		Crw(ir)=0d0
		tot=0d0
		do ilam=1,nlam-2
			Pl=Planck(T(ir),freq(ilam))
			do ig=1,ng
				Crw(ir)=Crw(ir)+dfreq(ilam)*Pl/Ce(ir,ilam,ig)
				tot=tot+dfreq(ilam)*Pl
			enddo
		enddo
		Crw(ir)=tot/Crw(ir)
	enddo

	do ir=1,nr
		dorw(ir)=((R(ir+1)-R(ir))*Crw(ir).gt.factRW)		
	enddo

	do ir=1,nr
		do ilam=1,nlam-1
			tot=Planck(T(ir),freq(ilam))
			do ig=1,ng
				spec(ilam,ig,ir)=Ca(ir,ilam,ig)*tot
			enddo
		enddo
	enddo

	if(TeffP.gt.0d0) then
	call output("Internal heating")
	EJv=0d0
	Lum=0d0
	do ilam=1,nlam-1
		specsource(ilam,1)=Planck(max(TeffP,3d0),freq(ilam))
		Lum=Lum+dfreq(ilam)*specsource(ilam,1)
	enddo
	Lum=Lum*pi*clight
	if(Lum.eq.0d0) then
		specsource(1:nlam-1,1)=0d0
	else
		specsource(1:nlam-1,1)=specsource(1:nlam-1,1)*sigma*max(TeffP,3d0)**4/Lum
	endif
	Lum=0d0
	do ilam=1,nlam-1
		Lum=Lum+specsource(ilam,1)*dfreq(ilam)
	enddo
	Nphot=NphotPlanet
	E0=Lum/real(Nphot)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(iphot,ilam,ig,x,y,z,jr,onedge,dx,dy,dz,goingup,v,jrnext,tau,tau_v,albedo)
!$OMP& SHARED(Nphot,R,nr,Ce,Cs,EJv,spec,ng,nlam,betaT,specsource,Crw,dorw)
!$OMP DO SCHEDULE (STATIC,1)
	do iphot=1,Nphot
c		call tellertje(iphot,Nphot)
		call emit(specsource,1,ilam,ig)
		ig=random(idum)*real(ng)+1
		call randomdirection(x,y,z)
		x=0d0
		y=0d0
		z=R(1)
		jr=1
		onedge=.true.
1		call randomdirection(dx,dy,dz)
		goingup=(dz.gt.0d0)
		if(.not.goingup) goto 1
		tau=-log(random(idum))
		do while(jr.le.nr)
			if(goingup) then
				v=(R(jr+1)-z)/dz
				jrnext=jr+1
			else
				v=(R(jr)-z)/dz
				jrnext=jr-1
			endif

			tau_v=v*Ce(jr,ilam,ig)
			albedo=Cs(jr,ilam)/Ce(jr,ilam,ig)
			if(tau_v.lt.tau) then
				x=x+v*dx
				y=y+v*dy
				z=z+v*dz
				tau=tau-tau_v
				EJv(jr)=EJv(jr)+tau_v*(1d0-albedo)
				jr=jrnext
				if(jr.lt.1) then
					jr=1
					dz=-dz
					goingup=.true.
				endif
				onedge=.true.
			else
				v=tau/Ce(jr,ilam,ig)
				x=x+v*dx
				y=y+v*dy
				z=z+v*dz
				EJv(jr)=EJv(jr)+tau*(1d0-albedo)
				call randomdirection(dx,dy,dz)
				if(random(idum).gt.albedo) then
					call emit(spec(1:nlam,1:ng,jr),ng,ilam,ig)
				endif
				if(dorw(jr)) then
					do while(RandomWalkT(x,y,z,dx,dy,dz,Crw(jr),jr,EJv(jr)))
					enddo
				endif
				tau=-log(random(idum))
				goingup=(dz.gt.0d0)
				onedge=.false.
			endif
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	EJvTot=EJvTot+EJv*E0
	endif

	call output("Stellar heating")
	EJv=0d0
	Lum=0d0
	do ilam=1,nlam-1
		specsource(ilam,1)=Planck(Tstar,freq(ilam))
		Lum=Lum+specsource(ilam,1)*dfreq(ilam)
	enddo
	Lum=Lum*pi*clight
	if(Lum.eq.0d0) then
		specsource(1:nlam-1,1)=0d0
	else
		specsource(1:nlam-1,1)=specsource(1:nlam-1,1)*sigma*Tstar**4/Lum
		specsource(1:nlam-1,1)=specsource(1:nlam-1,1)*betaT*(Rstar/Dplanet)**2
	endif
	Lum=0d0
	do ilam=1,nlam-1
		Lum=Lum+specsource(ilam,1)*dfreq(ilam)
	enddo
	Nphot=NphotStar
	E0=Lum/real(Nphot)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(iphot,ilam,ig,x,y,z,jr,onedge,dx,dy,dz,goingup,v,jrnext,tau,tau_v,albedo)
!$OMP& SHARED(Nphot,R,nr,Ce,Cs,EJv,spec,ng,nlam,betaT,specsource,Crw,dorw)
!$OMP DO SCHEDULE (STATIC,1)
	do iphot=1,Nphot
c		call tellertje(iphot,Nphot)
		call emit(specsource,1,ilam,ig)
		ig=random(idum)*real(ng)+1
		x=0d0
		y=0d0
		z=R(nr+1)
		dz=-1d0*betaT
		dx=sqrt(1d0-dz**2)
		dy=0d0
		goingup=.false.
		onedge=.true.
		jr=nr
		tau=-log(random(idum))
		do while(jr.le.nr)
			if(goingup) then
				v=(R(jr+1)-z)/dz
				jrnext=jr+1
			else
				v=(R(jr)-z)/dz
				jrnext=jr-1
			endif

			tau_v=v*Ce(jr,ilam,ig)
			albedo=Cs(jr,ilam)/Ce(jr,ilam,ig)
			if(tau_v.lt.tau) then
				x=x+v*dx
				y=y+v*dy
				z=z+v*dz
				tau=tau-tau_v
				EJv(jr)=EJv(jr)+tau_v*(1d0-albedo)
				jr=jrnext
				if(jr.lt.1) then
					jr=1
					dz=-dz
					goingup=.true.
				endif
				onedge=.true.
			else
				v=tau/Ce(jr,ilam,ig)
				x=x+v*dx
				y=y+v*dy
				z=z+v*dz
				EJv(jr)=EJv(jr)+tau*(1d0-albedo)
				call randomdirection(dx,dy,dz)
				if(random(idum).gt.albedo) then
					call emit(spec(1:nlam,1:ng,jr),ng,ilam,ig)
				endif
				if(dorw(jr)) then
					do while(RandomWalkT(x,y,z,dx,dy,dz,Crw(jr),jr,EJv(jr)))
					enddo
				endif
				tau=-log(random(idum))
				goingup=(dz.gt.0d0)
				onedge=.false.
			endif
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	EJvTot=EJvTot+EJv*E0


	chi2=0d0
	converged=.true.
	do ir=nr,1,-1
		CabsL=0d0
		do ilam=1,nlam-1
			do ig=1,ng
				CabsL(ilam)=CabsL(ilam)+Ca(ir,ilam,ig)/real(ng)
			enddo
		enddo

		Tmaxi=20000d0
		Tmini=3d0
		T0=T(ir)
		do i=1,20
			Lum=0d0
			tot=0d0
			do ilam=1,nlam-1
				specsource(ilam,1)=Planck(T0,freq(ilam))
				tot=tot+specsource(ilam,1)*dfreq(ilam)
				do ig=1,ng
					Lum=Lum+specsource(ilam,1)*dfreq(ilam)*Ca(ir,ilam,ig)/real(ng)
				enddo
			enddo
			tot=tot*pi*clight
			if(tot.eq.0d0) then
 				Lum=0d0
 			else
	 			Lum=Lum*sigma*T0**4/tot
			endif
			Lum=Lum*(R(ir+1)-R(ir))
			if(Lum.lt.EjvTot(ir)) then
				Tmini=T0
				T0=(T0+Tmaxi)/2d0
			else
				Tmaxi=T0
				T0=(T0+Tmini)/2d0
			endif
		enddo

		if(abs(T(ir)-T0)/(T(ir)+T0).gt.epsiter) converged=.false.
		chi2=chi2+((T(ir)-min(T0,2900d0))/((min(T0,2900d0)+T(ir))*0.1))**2
		T(ir)=sqrt(T(ir)*T0)
	enddo
	chi2=chi2/real(nr)
	call output("Chi2: " // trim(dbl2string(chi2,'(f7.2)')))

	deallocate(specsource)
	deallocate(spec)
	deallocate(Ca)
	deallocate(Ce)
	deallocate(Cs)
	deallocate(g)
	deallocate(M)
	deallocate(Eplanck)

c	do ir=1,nr
c		if(dochemistry) then
c			do i=1,nmol
c				call MorleyChemistry(mixrat_r(ir,i),T(ir),P(ir),molname(i),metallicity)
c			enddo
c		endif
c	enddo

	call WriteStructure

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


	subroutine increaseT(T1,T2,E,Eabs,CabsL,ir,Eplanck)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ilam,ir,iTmax,iTmin,iT
	real*8 T1,T2,E,E0,Eabs,CabsL(nlam),Ca0,Cs0,V,Eplanck(nr,nBB)

	V=4d0*pi*(R(ir+1)**3-R(ir)**3)/3d0

	E=E+Eabs

	iTmin=1

	iTmax=nBB

	iT=0.5d0*(iTmax+iTmin)

	do while(abs(iTmax-iTmin).gt.1)
		if(Eplanck(ir,iT).lt.0d0) then
			E0=0d0
			do ilam=1,nlam-1
				E0=E0+dfreq(ilam)*CabsL(ilam)*BB(iT,ilam)
			enddo
			E0=E0*V
			Eplanck(ir,iT)=E0
		else
			E0=Eplanck(ir,iT)
		endif
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
	

	subroutine emit(spec,ng0,ilam,ig)
	use GlobalSetup
	IMPLICIT NONE
	integer ilam,ng0,ig
	real*8 spec(nlam,ng0),Ltot,Lr
	real*8 random

	Ltot=0d0
	do ilam=1,nlam-1
		do ig=1,ng0
			Ltot=Ltot+dfreq(ilam)*spec(ilam,ig)
		enddo
	enddo

	Lr=random(idum)*Ltot

	Ltot=0d0
	do ilam=1,nlam-1
		do ig=1,ng0
			Ltot=Ltot+dfreq(ilam)*spec(ilam,ig)
			if(Ltot.gt.Lr) return
		enddo
	enddo

	return
	end
	
		
