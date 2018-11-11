	subroutine DoComputeT(converged,nTiter)
	use GlobalSetup
	use Constants
	use CloudModule
	IMPLICIT NONE
	integer iphase,nTiter,iter
	real*8 Planck,Pb(nr+1),f,dtau,expint,cp
	real*8,allocatable :: Ca(:,:,:),Cs(:,:),Ce(:,:,:),g(:,:),Tr(:,:,:)
	real*8,allocatable :: Fl(:,:,:,:),Ff(:,:),Pl(:,:),T0(:)
	real*8 tot,tot2,tot3,must,tstep,deltaF,err
	integer ir,ilam,ig,i
	logical docloud0(max(nclouds,1)),converged
	type(Mueller),allocatable :: M(:,:)
	
	allocate(Tr(nr,nlam,ng))
	allocate(T0(nr))
	allocate(Fl(2,nr+1,nlam,ng))
	allocate(Ca(nr,nlam,ng))
	allocate(Ce(nr,nlam,ng))
	allocate(Cs(nr,nlam))
	allocate(g(nr,nlam))
	allocate(M(nr,nlam))
	allocate(Ff(2,nr+1))
	allocate(Pl(nr,nlam))
	
	call output("Temperature computation (in beta phase!!)")

	docloud0=.false.
	if(nTiter.ne.0) then
		do i=1,nclouds
			if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
		enddo
	endif

	T0=T

	do ilam=1,nlam
		do ig=1,ng
			do ir=1,nr
				call Crossections(ir,ilam,ig,Ca(ir,ilam,ig),Cs(ir,ilam),docloud0)
				Ca(ir,ilam,ig)=Ca(ir,ilam,ig)/dens(ir)
				Cs(ir,ilam)=Cs(ir,ilam)/dens(ir)
				Ce(ir,ilam,ig)=Ca(ir,ilam,ig)+Cs(ir,ilam)
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

	do ir=1,nr
		do ilam=1,nlam-1
			do ig=1,ng
				dtau=abs(Pb(ir+1)-Pb(ir))*1d6*(Ca(ir,ilam,ig)+Cs(ir,ilam)*(1d0-g(ir,ilam)))/grav(ir)
				if(dtau.gt.0d0) then
					Tr(ir,ilam,ig)=(1d0-dtau)*exp(-dtau)+dtau**2*expint(1,dtau)
				else
					Tr(ir,ilam,ig)=1d0
				endif
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

	do iter=1,100

	Fl=0d0

	do ir=1,nr
		tot=0d0
		do ilam=1,nlam-1
			Pl(ir,ilam)=Planck(T(ir),freq(ilam))
			tot=tot+Pl(ir,ilam)*dfreq(ilam)
		enddo
		tot=tot*pi*clight
 		Pl(ir,1:nlam)=Pl(ir,1:nlam)*sigma*T(ir)**4/tot
	enddo

	do ilam=1,nlam-1
		Fl(1,nr+1,ilam,1:ng)=must*Fstar(ilam)/Dplanet**2
		do ir=nr,1,-1
			do ig=1,ng
				Fl(1,ir,ilam,ig)=Tr(ir,ilam,ig)*Fl(1,ir+1,ilam,ig)+pi*Pl(ir,ilam)*(1d0-Tr(ir,ilam,ig))
			enddo
		enddo
		Fl(2,1,ilam,1:ng)=Fl(1,1,ilam,1:ng)+pi*Planck(max(3d0,TeffP),freq(ilam))
		do ir=2,nr+1
			do ig=1,ng
				Fl(2,ir,ilam,ig)=Tr(ir-1,ilam,ig)*Fl(2,ir-1,ilam,ig)+pi*Pl(ir-1,ilam)*(1d0-Tr(ir-1,ilam,ig))
			enddo
		enddo
	enddo
	do ir=1,nr+1
		Ff(1,ir)=0d0
		Ff(2,ir)=0d0
		do ilam=1,nlam-1
			do ig=1,ng
				Ff(1,ir)=Ff(1,ir)+dfreq(ilam)*Fl(1,ir,ilam,ig)/real(ng)
				Ff(2,ir)=Ff(2,ir)+dfreq(ilam)*Fl(2,ir,ilam,ig)/real(ng)
			enddo
		enddo
	enddo
	
	do ir=1,nr
		cp=3.5*kb/(mp*MMW(ir))
		tstep=cp*P(ir)/(sigma*grav(ir)*T(ir)**3)
		deltaF=((Ff(2,ir+1)-Ff(1,ir+1))-(Ff(2,ir)-Ff(1,ir)))/abs(R(ir+1)-R(ir))
		tstep=tstep*1d5/(abs(deltaF)**0.9)
		T0(ir)=T0(ir)-deltaF*tstep/(dens(ir)*cp)
		if(.not.T0(ir).gt.3d0) T0(ir)=3d0
		err=((Ff(2,ir+1)-Ff(1,ir+1))-(Ff(2,ir)-Ff(1,ir)))/(sigma*T(ir)**4)
	enddo

	enddo

	chi2=0d0
	do ir=1,nr
		chi2=chi2+((min(T(ir),2900d0)-min(T0(ir),2900d0))/((min(T0(ir),2900d0)+min(T(ir),2900d0))*epsiter))**2
		f=max(0.5e0,real(nTiter)/real(maxiter))
		if(nTiter.ne.0) then
			T(ir)=T(ir)**(1d0-f)*T0(ir)**f
		else
			T(ir)=T0(ir)
		endif
	enddo
	chi2=chi2/real(nr)

	converged=.false.

	deallocate(Tr)
	deallocate(Fl)
	deallocate(Ca)
	deallocate(Ce)
	deallocate(Cs)
	deallocate(g)
	deallocate(M)
	deallocate(Ff)
	deallocate(Pl)

	call WriteStructure

	return
	end

