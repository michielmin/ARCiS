	module modComputeT
	IMPLICIT NONE
	real*8,allocatable :: CrV_prev(:),CrT_prev(:)
	end module modComputeT
	
	subroutine DoComputeT(converged,nTiter)
	use GlobalSetup
	use Constants
	use modComputeT
	IMPLICIT NONE
	integer iphase,nTiter,iter
	real*8 tau_V,tau_T,Planck,CrV(nr),CrT(nr),f
	real*8,allocatable :: Ca(:,:,:),Cs(:,:),Ce(:,:,:),g(:,:)
	real*8 tot,tot2,tot3,chi2,must,gamma,dP,Tirr,T0(nr),must_i
	integer ir,ilam,ig,i,iT
	logical docloud0(max(nclouds,1)),converged
	type(Mueller),allocatable :: M(:,:)
	
	allocate(Ca(nr,nlam,ng))
	allocate(Ce(nr,nlam,ng))
	allocate(Cs(nr,nlam))
	allocate(g(nr,nlam))
	allocate(M(nr,nlam))
	if(.not.allocated(CrV_prev)) allocate(CrV_prev(nr),CrT_prev(nr))
	
	call output("Temperature computation (in beta phase!!)")

	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo

	do ilam=1,nlam
		do ig=1,ng
			do ir=1,nr
				call Crossections(ir,ilam,ig,Ca(ir,ilam,ig),Cs(ir,ilam),docloud0)
				Ce(ir,ilam,ig)=Ca(ir,ilam,ig)+Cs(ir,ilam)
			enddo
		enddo
	enddo

	T0(1:nr)=T(1:nr)
	
	do iter=1,5

	do ir=1,nr
		iT=T0(ir)+1
		if(iT.gt.nBB-1) iT=nBB-1
		CrT(ir)=0d0
		tot2=0d0
		CrV(ir)=0d0
		tot3=0d0
		do ilam=1,nlam-1
			call GetMatrix(ir,ilam,M(ir,ilam),docloud0)
			g(ir,ilam)=0d0
			tot=0d0
			do iphase=1,180
				g(ir,ilam)=g(ir,ilam)+M(ir,ilam)%F11(iphase)*costheta(iphase)*sintheta(iphase)
				tot=tot+M(ir,ilam)%F11(iphase)*sintheta(iphase)
			enddo
			g(ir,ilam)=g(ir,ilam)/tot
			do ig=1,ng
				CrT(ir)=CrT(ir)+dfreq(ilam)*BB(iT,ilam)
     &				/(Ca(ir,ilam,ig)+Cs(ir,ilam)*(1d0-g(ir,ilam)))
				tot2=tot2+dfreq(ilam)*BB(iT,ilam)

				CrV(ir)=CrV(ir)+dfreq(ilam)*Fstar(ilam)
     &				/(Ca(ir,ilam,ig)+Cs(ir,ilam)*(1d0-g(ir,ilam)))
				tot3=tot3+dfreq(ilam)*Fstar(ilam)
			enddo
		enddo
		CrT(ir)=tot2/CrT(ir)
		CrV(ir)=tot3/CrV(ir)
	enddo
	CrT=CrT/dens
	CrV=CrV/dens

	if(nTiter.ne.0) then
		CrV=sqrt(CrV*CrV_prev)
		CrT=sqrt(CrT*CrT_prev)
	endif

	chi2=0d0
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
     &		3d0*Tirr**4*must*(2d0/3d0+must/gamma+(gamma/(3d0*must)-must/gamma)*exp(-tau_V/must))/4d0
		enddo
		T0(ir)=T0(ir)/100d0
		T0(ir)=T0(ir)**0.25
	enddo
	
	enddo

	chi2=0d0
	do ir=1,nr
		chi2=chi2+((min(T(ir),2900d0)-min(T0(ir),2900d0))/((min(T0(ir),2900d0)+min(T(ir),2900d0))*epsiter))**2
		f=real(nTiter)/real(maxiter)
		if(nTiter.ne.0) then
			T(ir)=T(ir)**(1d0-f)*T0(ir)**f
		else
			T(ir)=T0(ir)
		endif
	enddo
	chi2=chi2/real(nr)

	converged=.false.
	if(chi2.lt.1d0) converged=.true.
	call output("Chi2: " // trim(dbl2string(chi2,'(f7.2)')))

	CrV_prev=CrV
	CrT_prev=CrT

	deallocate(Ca)
	deallocate(Ce)
	deallocate(Cs)
	deallocate(g)
	deallocate(M)

	call WriteStructure

	return
	end

