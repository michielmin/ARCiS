	subroutine DoComputeT(converged,f)
	use GlobalSetup
	use Constants
	use CloudModule
	use TimingModule
	IMPLICIT NONE
	integer iphase,Titer
	real*8 tau_V,tau_T,Planck,CrV(nr),CrT(nr),CrVe(nr),f
	real*8 Cs,g,dlnT,dlnP
	real*8,allocatable :: Ce(:,:,:),Ca(:,:,:)
	real*8 tot,tot2,tot3,tot4,chi2,must,gamma,dP,Tirr,T0(nr),must_i
	integer ir,ilam,ig,i,iT,niter
	logical docloud0(max(nclouds,1)),converged
	real*8 time
	integer itime

	call cpu_time(time)
	timetemp=timetemp-time
	call system_clock(itime)
	itimetemp=itimetemp-itime
	ctimetemp=ctimetemp+1

	call DoComputeTeddington(converged,f)

	call cpu_time(time)
	timetemp=timetemp+time
	call system_clock(itime)
	itimetemp=itimetemp+itime

	return
	end


	subroutine DoComputeTeddington(converged,f)
	use GlobalSetup
	use Constants
	use CloudModule
	use Struct3D
	IMPLICIT NONE
	integer iphase,iter,iter2
	real*8 tau_V,tau_T,Planck,f
	real*8 g,dlnT,dlnP,d,tau,tautot,fact,contr,tau_a,exp_tau
	real*8,allocatable :: Ce(:,:,:),Ca(:,:,:),Cs(:,:,:),taustar(:,:),tauR_nu(:,:,:)
	real*8,allocatable :: wscat(:,:,:),wabs(:,:,:)
	real*8 tot,tot2,tot3,tot4,chi2,must,gamma,dP,Tirr,T0(nr),T1(nr),must_i,E,E0,Tinp(nr)
	real*8 z,Hstar(nr),Jtot,Htot,Ktot
	real*8 fedd(nr),Hedd(nr),lH1,lH2,P1,P2
	real*8,allocatable :: Jnu(:,:,:),Hnu(:,:,:),Knu(:,:,:)
	integer ir,ilam,ig,i,iT,niter,inu,nnu,jr,iTmin,iTmax
	logical docloud0(max(nclouds,1)),converged,stopscat
	real*8 tauf(nr),Si(0:nr),B1,B2,x1,x2,dx1,dx2,ax,bx,ff,TT,maxErr,Fl0(nr),IntH0(nr,nr)
	integer info,IWORK(10*(nr+1)*(nr+1)),NRHS,ii(3),iscat,nscat
	real*8 tau1,tau2,ee0,ee1,ee2,tauR(nr),Ij(nr),Ih(nr),scale,dtauR(nr),EabDirect(nr)
	integer nlam_LR
	real*8 IntH(nr,nr),Fl(nr),Ts(nr),minFl(nr),maxFl(nr),maxfact,err,tot1,Pb(nr+1)
	real*8,allocatable :: lam_LR(:),dfreq_LR(:),freq_LR(:),BB_LR(:,:),IntHnu(:,:,:),dtauR_nu(:,:,:)
	integer i1,i2,ngF,j
	real*8 ww,w1,w2,FstarBottom,tauRoss
	real*8,allocatable :: Si_omp(:,:),Ih_omp(:),Ij_omp(:),tauR_omp(:),Hsurf(:,:),Hstar_omp(:),EabDirect_omp(:)
	real*8,allocatable :: Fstar_LR(:),SurfEmis_LR(:),IntEab(:,:,:),IntHnuSurf(:,:),IntEabSurf(:,:)
	integer,allocatable :: IP(:)
	real*8,allocatable :: WS(:),nu(:),wnu(:),Hstar_lam(:),Hsurf_lam(:)
	real*8,allocatable :: IhN(:,:,:),IjN(:,:,:),HBottom(:)
	character*500 file
	integer icloud
	logical Convec(0:nr)
	real*8,allocatable,save :: Hstar0(:),EabDirect0(:)

	if(deepredist.and.deepredisttype.eq.'fixflux') then
		if(.not.allocated(Hstar0)) allocate(Hstar0(nr))
		if(.not.allocated(EabDirect0)) allocate(EabDirect0(nr))
	endif

	maxfact=4d0

	allocate(IP(nr*4+2))
	IP(1)=(2*(nr)+nr*3+(nr*2+2)*(nr+7))
	IP(2)=nr*4+2
	allocate(WS(IP(1)))
	
	nlam_LR=0
	do ilam=1,nlam
		if(RTgridpoint(ilam)) then
			nlam_LR=nlam_LR+1
		endif
	enddo

	allocate(lam_LR(nlam_LR))
	allocate(freq_LR(nlam_LR))
	allocate(dfreq_LR(nlam_LR))
	nlam_LR=0
	do ilam=1,nlam
		if(RTgridpoint(ilam)) then
			nlam_LR=nlam_LR+1
			lam_LR(nlam_LR)=lam(ilam)
			freq_LR(nlam_LR)=freq(ilam)
			dfreq_LR(nlam_LR)=dfreq(ilam)
		endif
	enddo

	allocate(BB_LR(nlam_LR,nBB))

	BB_LR=0d0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(TT,j,i,tot,scale)
!$OMP& SHARED(nlam_LR,BB_LR,freq_LR,dfreq_LR)
!$OMP DO
	do j=nBB,1,-1
		TT=real(j)
		tot=0d0
		do i=1,nlam_LR
			BB_LR(i,j)=Planck(TT,freq_LR(i))
			tot=tot+dfreq_LR(i)*BB_LR(i,j)
		enddo
		scale=((2d0*(pi*kb*TT)**4)/(15d0*hplanck**3*clight**3))/tot
		if(scale.gt.0d0.and.tot.gt.0d0) then
			BB_LR(1:nlam_LR,j)=BB_LR(1:nlam_LR,j)*scale
		else if(j.lt.nBB) then
			scale=((2d0*(pi*kb*TT)**4)/(15d0*hplanck**3*clight**3))/((2d0*(pi*kb*real(j+1))**4)/(15d0*hplanck**3*clight**3))
			BB_LR(1:nlam_LR,j)=BB_LR(1:nlam_LR,j+1)*scale
		endif
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	allocate(taustar(nlam_LR,ng))
	allocate(Ce(nr,nlam_LR,ng))
	allocate(Ca(nr,nlam_LR,ng))
	allocate(Cs(nr,nlam_LR,ng))
	allocate(Fstar_LR(nlam_LR))
	allocate(tauR_nu(0:nr,nlam_LR,ng))
	allocate(dtauR_nu(0:nr,nlam_LR,ng))
	allocate(IntHnu(nlam_LR,nr,nr))
	allocate(IntHnuSurf(nlam_LR,nr))
	allocate(IntEab(nlam_LR,nr,nr))
	allocate(IntEabSurf(nlam_LR,nr))
	allocate(SurfEmis_LR(nlam_LR))

	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo
	
	niter=50
	nscat=10

	must=betaT

	do ir=nr,1,-1
		if(T(ir).gt.0d0) then
			Tinp(ir)=T(ir)
		else if(ir.lt.nr) then
			Tinp(ir)=T(ir+1)
		else
			Tinp(ir)=abs(must)*sqrt(Rstar/(2d0*Dplanet))*Tstar
		endif
		T(ir)=Tinp(ir)
	enddo

	i=0
	do ilam=1,nlam
		if(RTgridpoint(ilam)) then
			i=i+1
			SurfEmis_LR(i)=surface_emis(ilam)
			Fstar_LR(i)=Fstar(ilam)
			do ir=1,nr
				do ig=1,ng
					call Crossections(ir,ilam,ig,Ca(ir,i,ig),Cs(ir,i,ig),docloud0)
					Ca(ir,i,ig)=Ca(ir,i,ig)/dens(ir)
					Cs(ir,i,ig)=Cs(ir,i,ig)/dens(ir)
				enddo
			enddo
		endif
	enddo

	Ce=Ca+Cs
	allocate(wscat(nr,nlam_LR,ng),wabs(nr,nlam_LR,ng))
	do ir=1,nr
		do ilam=1,nlam_LR
			do ig=1,ng
				wabs(ir,ilam,ig)=Ca(ir,ilam,ig)/Ce(ir,ilam,ig)
				if(.not.wabs(ir,ilam,ig).gt.1d-4) then
					wabs(ir,ilam,ig)=1d-4
					Ca(ir,ilam,ig)=Cs(ir,ilam,ig)/(1d0/wabs(ir,ilam,ig)-1d0)
					Ce(ir,ilam,ig)=Ca(ir,ilam,ig)+Cs(ir,ilam,ig)
				endif
				wscat(ir,ilam,ig)=Cs(ir,ilam,ig)/Ce(ir,ilam,ig)
				if(.not.wscat(ir,ilam,ig).gt.0d0) then
					wscat(ir,ilam,ig)=0d0
					Cs(ir,ilam,ig)=0d0
					Ce(ir,ilam,ig)=Ca(ir,ilam,ig)+Cs(ir,ilam,ig)
				endif
			enddo
		enddo
	enddo

	converged=.true.

	Pb(1)=P(1)
	do i=2,nr
		Pb(i)=sqrt(P(i-1)*P(i))
	enddo
	Pb(nr+1)=0d0
	do ilam=1,nlam_LR
		do ig=1,ng
			do jr=nr,1,-1
				if(jr.eq.nr) then
					d=(P(jr)-Pb(jr+1))*1d6/grav(jr)
					tau=d*Ce(jr,ilam,ig)
				else
					d=(Pb(jr+1)-P(jr+1))*1d6/grav(jr+1)
					tau=d*Ce(jr+1,ilam,ig)
					d=(P(jr)-Pb(jr+1))*1d6/grav(jr)
					tau=tau+d*Ce(jr,ilam,ig)
				endif
				if(.not.tau.gt.0d0) then
					tau=0d0
				endif
				tau=tau+1d-10
				tau=1d0/(1d0/tau+1d-10)
				dtauR_nu(jr,ilam,ig)=tau
			enddo
			dtauR_nu(nr,ilam,ig)=dtauR_nu(nr-1,ilam,ig)**2/dtauR_nu(nr-2,ilam,ig)
			do jr=nr,1,-1
				if(jr.eq.nr) then
					tauR(jr)=dtauR_nu(nr,ilam,ig)
				else
					tauR(jr)=tauR(jr+1)+dtauR_nu(jr,ilam,ig)
				endif
			enddo
			tauR_nu(1:nr,ilam,ig)=tauR(1:nr)
		enddo
	enddo

	nnu=5
	allocate(nu(nnu),wnu(nnu))
	call gauleg(0d0,1d0,nu,wnu,nnu)

c===============================================================
c quick thing to read in a file!
c	file='sedlhs.txt'
c	call regridlog(file,1d4*lam_LR,Fstar_LR,nlam_LR)
c	Fstar_LR=Fstar_LR*distance**2/1e23
c===============================================================
	E=0d0
	do ilam=1,nlam_LR
c		Fstar_LR(ilam)=Planck(Tstar,freq_LR(ilam))*pi*Rstar**2
		E=E+dfreq_LR(ilam)*Fstar_LR(ilam)
	enddo
	scale=pi*Rstar**2*((2d0*(pi*kb*Tstar)**4)/(15d0*hplanck**3*clight**3))/E
	Fstar_LR=Fstar_LR*scale
	
	ff=1d0
	
	Hstar(1:nr)=0d0
	IntHnu(1:nlam_LR,1:nr,1:nr)=0d0
	IntEab(1:nlam_LR,1:nr,1:nr)=0d0
	IntHnuSurf(1:nlam_LR,1:nr)=0d0
	EabDirect(1:nr)=0d0
	IntEabSurf(1:nlam_LR,1:nr)=0d0

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(Si_omp,tauR_omp,Ih_omp,Ij_omp,ilam,ig,ir,inu,jr,EabDirect_omp,HBottom,
!$OMP&			Hstar_omp,contr,FstarBottom,Hstar_lam,Hsurf_lam,tot,IhN,IjN)
!$OMP& SHARED(nlam_LR,ng,nr,nnu,tauR_nu,nu,wnu,dfreq_LR,wgg,IntHnu,SurfEmis_LR,dtauR_nu,Ca,Ce,Cs,Hsurf,
!$OMP&			Hstar,Dplanet,Fstar_LR,must,wabs,wscat,EabDirect,IntEab,IntHnuSurf,IntEabSurf,betaF,isoFstar)
	allocate(Si_omp(nr,0:nr+1),tauR_omp(nr),Ih_omp(nr),Ij_omp(nr))
	allocate(IhN(nr,0:nr+1,nnu),IjN(nr,0:nr+1,nnu))
	allocate(Hstar_omp(nr),Hstar_lam(nr),Hsurf_lam(nr),EabDirect_omp(nr),HBottom(nr))
	Hstar_omp=0d0
	EabDirect_omp=0d0
!$OMP DO
	do ilam=1,nlam_LR
		call tellertje(ilam,nlam_LR)
		do ig=1,ng
			Hstar_lam(1:nr)=0d0
			Hsurf_lam(1:nr)=0d0
			contr=(Fstar_LR(ilam)/(pi*Dplanet**2))
			if(isoFstar) then
				FstarBottom=0d0
				do inu=1,nnu
					tauR_omp(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(nu(inu))
					Ij_omp(1:nr)=contr*exp(-tauR_omp(1:nr))
					Ih_omp(1:nr)=-betaF*Ij_omp(1:nr)*nu(inu)
					Hstar_lam(1:nr)=Hstar_lam(1:nr)+2d0*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
					EabDirect_omp(1:nr)=EabDirect_omp(1:nr)+wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)*Ca(1:nr,ilam,ig)
					FstarBottom=FstarBottom+2d0*wnu(inu)*abs(Ih_omp(1))
				enddo
			else
				tauR_omp(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(max(must,1d-5))
				Ij_omp(1:nr)=contr*exp(-tauR_omp(1:nr))
				Ih_omp(1:nr)=-betaF*Ij_omp(1:nr)
				Hstar_lam(1:nr)=Hstar_lam(1:nr)+dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
				EabDirect_omp(1:nr)=EabDirect_omp(1:nr)+dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)*Ca(1:nr,ilam,ig)
				FstarBottom=abs(Ih_omp(1))
			endif

c Si_omp(0:nr,0) is the direct stellar contribution
			Si_omp(1:nr,0)=Ij_omp(1:nr)*wscat(1:nr,ilam,ig)/2d0

c Si_omp(0:nr,1:nr) is the direct contribution from the atmosphere
			do ir=1,nr
				Si_omp(1:nr,ir)=0d0
				Si_omp(ir,ir)=wabs(ir,ilam,ig)
			enddo

c Si_omp(0:nr,nr+1) is the direct contribution from the surface
			Si_omp(1:nr,nr+1)=0d0
			do inu=1,nnu
				tauR_omp(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(nu(inu))
				Ij_omp(1:nr)=exp(-abs(tauR_omp(1:nr)-tauR_omp(1)))
				Ih_omp(1:nr)=Ij_omp(1:nr)
				Si_omp(1:nr,nr+1)=Si_omp(1:nr,nr+1)+wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)/2d0
				Hsurf_lam(1:nr)=Hsurf_lam(1:nr)+2d0*FstarBottom*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*
     & 								Ih_omp(1:nr)*(1d0-SurfEmis_LR(ilam))
				EabDirect_omp(1:nr)=EabDirect_omp(1:nr)+2d0*FstarBottom*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*
     & 								Ij_omp(1:nr)*(1d0-SurfEmis_LR(ilam))*Ca(1:nr,ilam,ig)
				IntHnuSurf(ilam,1:nr)=IntHnuSurf(ilam,1:nr)+2d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
				IntEabSurf(ilam,1:nr)=IntEabSurf(ilam,1:nr)+2d0*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)*Ca(1:nr,ilam,ig)
			enddo
			Si_omp(1:nr,nr+1)=Si_omp(1:nr,nr+1)*wscat(1:nr,ilam,ig)

			call AddScatter(Si_omp(1:nr,0:nr+1),tauR_nu(1:nr,ilam,ig),
     &					Ca(1:nr,ilam,ig),Cs(1:nr,ilam,ig),Ce(1:nr,ilam,ig),(1d0-SurfEmis_LR(ilam)),nr,nu,wnu,nnu,nr+2)

			do inu=1,nnu
				tauR_omp(1:nr)=dtauR_nu(1:nr,ilam,ig)/abs(nu(inu))
				call SolveIjhExpN(tauR_omp(1:nr),Si_omp(1:nr,0:nr+1),IhN(1:nr,0:nr+1,inu),IjN(1:nr,0:nr+1,inu),nr,nr+2)
			enddo

			do inu=1,nnu
				Ih_omp(1:nr)=IhN(1:nr,0,inu)
				Ij_omp(1:nr)=IjN(1:nr,0,inu)
				Hstar_lam(1:nr)=Hstar_lam(1:nr)+2d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
				EabDirect_omp(1:nr)=EabDirect_omp(1:nr)+2d0*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)*Ca(1:nr,ilam,ig)
			enddo
			HBottom=0d0
			do ir=1,nr
				do inu=1,nnu
					Ih_omp(1:nr)=IhN(1:nr,ir,inu)
					Ij_omp(1:nr)=IjN(1:nr,ir,inu)
					IntHnu(ilam,1:nr,ir)=IntHnu(ilam,1:nr,ir)+4d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
					IntEab(ilam,1:nr,ir)=IntEab(ilam,1:nr,ir)+4d0*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)*Ca(1:nr,ilam,ig)
					HBottom(ir)=HBottom(ir)+4d0*nu(inu)*wnu(inu)*Ij_omp(1)
				enddo
			enddo

			do inu=1,nnu
				tauR_omp(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(nu(inu))
				Ij_omp(1:nr)=exp(-abs(tauR_omp(1:nr)-tauR_omp(1)))
				Ih_omp(1:nr)=Ij_omp(1:nr)
				do ir=1,nr
					IntHnu(ilam,1:nr,ir)=IntHnu(ilam,1:nr,ir)+2d0*HBottom(ir)*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*
     &						Ih_omp(1:nr)*(1d0-SurfEmis_LR(ilam))
					IntEab(ilam,1:nr,ir)=IntEab(ilam,1:nr,ir)+2d0*HBottom(ir)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*
     &						Ij_omp(1:nr)*Ca(1:nr,ilam,ig)*(1d0-SurfEmis_LR(ilam))
				enddo
			enddo
			do inu=1,nnu
				Ih_omp(1:nr)=IhN(1:nr,nr+1,inu)
				Ij_omp(1:nr)=IjN(1:nr,nr+1,inu)
				Hsurf_lam(1:nr)=Hsurf_lam(1:nr)+FstarBottom*4d0*nu(inu)*wnu(inu)*Ih_omp(1:nr)*(1d0-SurfEmis_LR(ilam))
				IntHnuSurf(ilam,1:nr)=IntHnuSurf(ilam,1:nr)+4d0*nu(inu)*wnu(inu)*Ih_omp(1:nr)
				IntEabSurf(ilam,1:nr)=IntEabSurf(ilam,1:nr)+4d0*wnu(inu)*Ij_omp(1:nr)*Ca(1:nr,ilam,ig)
				EabDirect_omp(1:nr)=EabDirect_omp(1:nr)+FstarBottom*4d0*wnu(inu)*
     &						Ij_omp(1:nr)*(1d0-SurfEmis_LR(ilam))*Ca(1:nr,ilam,ig)
			enddo

			do ir=1,nr
				Hstar_omp(ir)=Hstar_omp(ir)+min(Hstar_lam(ir),0d0)+max(Hsurf_lam(ir),0d0)
			enddo
		enddo
	enddo
!$OMP END DO
!$OMP CRITICAL
	Hstar=Hstar+Hstar_omp
	EabDirect=EabDirect+EabDirect_omp
!$OMP END CRITICAL
	deallocate(Si_omp,tauR_omp,Ih_omp,Ij_omp,IhN,IjN)
	deallocate(Hstar_omp,Hstar_lam,Hsurf_lam,EabDirect_omp,HBottom)
!$OMP FLUSH
!$OMP END PARALLEL

	if(do3D.and.deepredist.and.deepredisttype.eq.'fixflux') then
		if(i3D.eq.n3D) then
			Hstar0(1:nr)=Hstar(1:nr)/betaF
			EabDirect0(1:nr)=EabDirect(1:nr)
		else
			Hstar(1:nr)=Hstar0(1:nr)*betaF
			EabDirect(1:nr)=EabDirect0(1:nr)
		endif
	endif

	iter=1
	iter2=1
	Tsurface=T(1)
	do while(iter.le.niter.and.iter2.le.niter*5)
	iter=iter+1

	Ts(1:nr)=T(1:nr)

	Hedd=0d0
	E0=2d0*(((pi*kb*TeffP)**4)/(15d0*hplanck**3*clight**3))
	do ir=1,nr
		Hedd(ir)=Hedd(ir)+E0-Hstar(ir)
	enddo

	IntH=0d0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(tot,iT,scale,jr,ir)
!$OMP& SHARED(nlam_LR,BB_LR,IntH,IntHnu,nr,T,IntHnuSurf,SurfEmis_LR,dfreq_LR)
!$OMP DO
	do jr=1,nr
		iT=T(jr)+1
		if(iT.gt.nBB-1) iT=nBB-1
		if(iT.lt.1) iT=1
		scale=(T(jr)/real(iT))**4
		do ir=1,nr
			tot=0d0
			do ilam=1,nlam_LR
				tot=tot+scale*BB_LR(ilam,iT)*IntHnu(ilam,ir,jr)
			enddo
			IntH(ir,jr)=tot
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	iT=Tsurface+1
	if(iT.gt.nBB-1) iT=nBB-1
	if(iT.lt.1) iT=1
	scale=(Tsurface/real(iT))**4
	do ir=1,nr
		do ilam=1,nlam_LR
			IntH(ir,1)=IntH(ir,1)+scale*BB_LR(ilam,iT)*IntHnuSurf(ilam,ir)*SurfEmis_LR(ilam)
		enddo
	enddo

	Fl=Hedd
	
	do ir=1,nr-1
		dlnP=log(P(ir+1)/P(ir))
		dlnT=(nabla_ad(ir))*dlnP
		maxFl(ir)=(exp(dlnT)*(T(ir)/T(ir+1)))**4d0	!Fl(ir+1)/Fl(ir)
	enddo
	
	minFl=0.5d0
	maxFl=2.0

	IP(1)=(2*(nr)+nr*3+(nr*2+2)*(nr+7))
	IP(2)=nr*4+2
	NRHS=1

	Fl0=Fl
	IntH0=IntH
	Convec(0:nr)=.false.

1	continue
	iter2=iter2+1

	call DGESV( nr, NRHS, IntH, nr, IWORK, Fl, nr, info )
c	call PosSolve(IntH,Fl,minFl,maxFl,nr,IP,WS)

	do ir=1,nr
		Fl(ir)=min(max(0.5d0,Fl(ir)),2d0)
	enddo

	maxErr=0d0
	j=1
	do ir=1,nr-2
		dlnP=log(P(ir+1)/P(ir))
		dlnT=log(T(ir+1)/T(ir))+0.25*log(Fl(ir+1)/Fl(ir))
		err=(dlnT/dlnP)/nabla_ad(ir)
		if(abs(err).gt.abs(maxErr).and..not.Convec(ir)) then
c		if(err.gt.maxErr.and..not.Convec(ir)) then
			maxErr=err
			j=ir
		endif
	enddo
	if(abs(maxErr).gt.1d0) then
		dlnP=log(P(j+1)/P(j))
		maxErr=nabla_ad(j)*abs(maxErr)/maxErr
		maxFl(j)=exp(4d0*(dlnP*maxErr-log(T(j+1)/T(j))))	!Fl(ir+1)/Fl(ir)
		
		IntH0(j,1:nr)=0d0
		IntH0(j,j)=maxFl(j)
		IntH0(j,j+1)=-1d0
		Fl0(j)=0d0
		Fl=Fl0
		IntH=IntH0
		Convec(j)=.true.
		goto 1
	endif		
				
	Ts=T*Fl**0.25

	maxErr=0d0
	do ir=1,nr-1
		if(abs(T(ir)-Ts(ir))/(T(ir)+Ts(ir)).gt.maxErr) maxErr=abs(T(ir)-Ts(ir))/(T(ir)+Ts(ir))
	enddo

	do ir=1,nr
		if(.not.Ts(ir).gt.3d0) Ts(ir)=3d0
		T(ir)=sqrt(Ts(ir)*T(ir))
	enddo
	Tsurface=T(1)

	if(maxErr.lt.(epsiter/5d0).and.iter.gt.5) exit
	enddo

	maxErr=0d0
	do ir=1,nr-1
		if(abs(T(ir)-Tinp(ir))/(T(ir)+Tinp(ir)).gt.maxErr) maxErr=abs(T(ir)-Tinp(ir))/(T(ir)+Tinp(ir))
		if(maxErr.gt.epsiter) converged=.false.
	enddo
	if(.not.allocated(Tdist)) allocate(Tdist(nr,maxiter))
	Tdist(1:nr,nTiter)=T(1:nr)
	call output("Maximum error on T-struct: " // dbl2string(maxErr*100d0,'(f5.1)') // "%")
	if(do3D.and..not.retrieval) print*,"Maximum error on T-struct: " // dbl2string(maxErr*100d0,'(f5.1)') // "%"
	T0(1:nr)=Tinp(1:nr)
	T1(1:nr)=T(1:nr)
	do ir=1,nr
		call computeav50(Tdist(ir,1:nTiter),nTiter,T1(ir))
	enddo
	do ir=1,nr
		T(ir)=f*T1(ir)+(1d0-f)*T0(ir)
	enddo

	Tsurface=T(1)
	call output("Surface temperature: " // dbl2string(Tsurface,'(f8.2)') // " K")
	if(do3D.and..not.retrieval) print*,"Surface temperature: " // dbl2string(Tsurface,'(f8.2)') // " K"

	if(.not.retrieval) then
		open(unit=26,file=trim(outputdir) // 'convection.dat',FORM="FORMATTED",ACCESS="STREAM")
		do ir=1,nr
			if(Convec(ir)) then
				write(26,*) T(ir),P(ir)
			else if(ir.ne.1) then
				if(Convec(ir-1)) then
					write(26,*) T(ir),P(ir)
				else
					write(26,*)
				endif
			else
				write(26,*)
			endif
		enddo
		close(unit=26)
	
		call WriteStructure
	endif
	call tellertje(niter,niter)

	deallocate(lam_LR)
	deallocate(freq_LR)
	deallocate(dfreq_LR)
	deallocate(BB_LR)
	deallocate(taustar)
	deallocate(Ce)
	deallocate(Ca)
	deallocate(Cs)
	deallocate(Fstar_LR)
	deallocate(tauR_nu)
	deallocate(IntHnu)
	deallocate(IntEab)
	deallocate(IntEabSurf)
	deallocate(IntHnuSurf)
	deallocate(SurfEmis_LR)

	return
	end


	subroutine ComputeDeriv(x,y,dy,n,yp1,ypn)
	IMPLICIT NONE
	integer n,i
	real*8 x(n),y(n),dy(n),d2y(n),yp1,ypn
	real*8 A,B

	do i=1,n-1
		dy(i)=(y(i+1)-y(i))/(x(i+1)-x(i))
	enddo
	dy(n)=ypn
	return

	call spline(x,y,n,yp1,ypn,d2y)
	dy(1)=yp1
	dy(n)=ypn
	do i=2,n-1
		dy(i)=(y(i+1)-y(i))/(x(i+1)-x(i))+(x(i+1)-x(i))*(d2y(i)-d2y(i+1))/24d0
	enddo

	return
	end
	

	subroutine ComputeI12(tau1,tau2,S1,S2,I12)
	IMPLICIT NONE
	real*8 tau1,tau2,S1,S2,I12,dtau,exptau
	
	dtau=abs(tau1-tau2)

	if(dtau.lt.1d-4) then
		I12=0.5*(S1+S2)*dtau
	else if(dtau.gt.1d4) then
		I12=S2+(S1-S2)/dtau
	else
		I12=(-(S1*dtau-S2+S1)*exp(-dtau)+S2*dtau-S2+S1)/dtau
	endif

	return
	end


	subroutine SolveIjhN(dtau,Si,Ih,Ij,nr,nrhs)
	IMPLICIT NONE
	integer ir,nr,jr,nrhs,i,info
	real*8 Ih(nr,nrhs),Ij(nr,nrhs),Si(nr,nrhs),dtau(nr)
	real*8 y(nr),Ma(nr),Mb(nr),Mc(nr),x(nr),fact,yp1,ypn
	real*8 Ij1(nrhs),IjN(nrhs),Ih1(nrhs),IhN(nrhs)

	call SolveIjhExpN(dtau,Si,Ih,Ij,nr,nrhs)
	Ij1=Ij(1,1:nrhs)
	Ih1=Ih(1,1:nrhs)
	IjN=Ij(nr,1:nrhs)
	IhN=Ih(nr,1:nrhs)

	Ma=0d0
	Mb=0d0
	Mc=0d0
	ir=1
	do ir=2,nr-1
		fact=1d0/(0.5d0*(dtau(ir)+dtau(ir-1)))
		Mb(ir)=1d0+fact*(1d0/(dtau(ir))+1d0/(dtau(ir-1)))
		Ma(ir)=-fact*1d0/(dtau(ir-1))
		Mc(ir)=-fact*1d0/(dtau(ir))
	enddo
	Mb(1)=1d0
	Mb(nr)=1d0

	do i=1,nrhs
		x(2:nr-1)=Si(2:nr-1,i)
		x(1)=Ij1(i)
		x(nr)=IjN(i)
		info=0
		call tridag(Ma,Mb,Mc,x,y,nr,info)
		Ij(1:nr,i)=y(1:nr)
	enddo

	do ir=nr-1,2,-1
		Ih(ir,1:nrhs)=(Ij(ir-1,1:nrhs)-Ij(ir,1:nrhs))/dtau(ir-1)
	enddo

	x(nr)=0d0
	do i=nr-1,1,-1
		x(i)=x(i+1)+dtau(i)
	enddo
	do i=1,nrhs
		x(2:nr-1)=-(Si(2:nr-1,i)-Si(1:nr-2,i))/dtau(2:nr-1)
		x(1)=Ih1(i)
		x(nr)=IhN(i)
		info=0
		call tridag(Ma,Mb,Mc,x,y,nr,info)
		Ih(1:nr,i)=y(1:nr)
	enddo

	return
	end
	


	subroutine SolveIjhExpN(dtau,Si,Ih,Ij,nr,nrhs)
	IMPLICIT NONE
	integer ir,nr,jr,nrhs
	real*8 Ih(nr,nrhs),Ij(nr,nrhs),Si(nr,nrhs),dtau(nr)
	real*8 exptau(nr),s0(nrhs),s1(nrhs),x1,yj(nrhs),yh(nrhs),x2
	real*8 Ijm(nr,nrhs),Ijp(nr,nrhs),Ihm(nr,nrhs),Ihp(nr,nrhs)

	do ir=1,nr
		exptau(ir)=exp(-dtau(ir))
	enddo

	Ijm=0d0
	Ihm=0d0
	do ir=nr,1,-1
		if(ir.eq.nr) then
			s0=Si(nr,1:nrhs)
		else
			s0=Si(ir+1,1:nrhs)
		endif
		s1=Si(ir,1:nrhs)
		x1=dtau(ir)
		if(x1.lt.1d-4) then
			yh=0.5d0*x1*s0
			yj=0.5d0*x1*(s0+s1)
		else if(x1.gt.1d4) then
			yh=s0/x1
			yj=s1+(s0-s1)/x1
		else
			yh=(-(s0*x1+s0)*exptau(ir)+s0)/x1
			yj=(-(s0*x1-s1+s0)*exptau(ir)+s1*x1-s1+s0)/x1
		endif
		if(ir.eq.nr) then
			Ijm(ir,1:nrhs)=yj(1:nrhs)
			Ihm(ir,1:nrhs)=yh(1:nrhs)
		else
			Ijm(ir,1:nrhs)=Ijm(ir+1,1:nrhs)*exptau(ir)+yj(1:nrhs)
			Ihm(ir,1:nrhs)=Ijm(ir+1,1:nrhs)*exptau(ir)+yh(1:nrhs)
		endif
	enddo

	Ijp=0d0
	Ihp=0d0
	do ir=1,nr-1
		s0=Si(ir,1:nrhs)
		s1=Si(ir+1,1:nrhs)
		x1=dtau(ir)
		if(x1.lt.1d-4) then
			yh=0.5d0*x1*s0
			yj=0.5d0*x1*(s0+s1)
		else if(x1.gt.1d4) then
			yh=s0/x1
			yj=s1+(s0-s1)/x1
		else
			yh=(-(s0*x1+s0)*exptau(ir)+s0)/x1
			yj=(-(s0*x1-s1+s0)*exptau(ir)+s1*x1-s1+s0)/x1
		endif
		Ijp(ir+1,1:nrhs)=Ijp(ir,1:nrhs)*exptau(ir)+yj(1:nrhs)
		Ihp(ir+1,1:nrhs)=Ijp(ir,1:nrhs)*exptau(ir)+yh(1:nrhs)
	enddo

	Ij=(Ijp+Ijm)/2d0

	Ih=(Ihp-Ihm)/2d0

	do ir=1,nr
		s1=Si(ir,1:nrhs)
		x1=dtau(ir)
		if(ir.gt.1) then
			x2=dtau(ir-1)
			if(x1.lt.1d-4.and.x2.lt.1d-4) then
				yh=0.5d0*(x2-x1)*s1
			else
				yh=s1*((exptau(ir-1)-1d0)*x1-(exptau(ir)-1d0)*x2)
				yh=yh/(x1*x2)
			endif
		else
			if(x1.lt.1d-4) then
				yh=-0.5d0*x1*s1
			else if(x1.gt.1d4) then
				yh=-s1+s1/x1
			else
				yh=-s1*(exptau(ir)+x1-1d0)/x1
			endif
		endif
		Ih(ir,1:nrhs)=Ih(ir,1:nrhs)+yh(1:nrhs)/2d0
	enddo

c	x1=dtau(nr)
c	do ir=nr-1,2,-1
c		x1=x1+dtau(ir)
c		if(x1.gt.1d-6) then
c			s0=Ij(ir,1:nrhs)
c			s1=Ij(ir+1,1:nrhs)
c			Ih(ir,1:nrhs)=(s0-s1)/dtau(ir)
c		endif
c	enddo

	return
	end
	



	subroutine SolveIjExp(tau,Si,Ij,Ih,nr)
	IMPLICIT NONE
	integer ir,nr,jr
	real*8 Ij(nr),Si(nr),tau(nr),Ip(nr),Im(nr),Ih(nr)
	real*8 exptau(nr),dtau(nr),s0,s1,x1,y

	dtau(nr)=tau(nr)
	do ir=nr-1,1,-1
		dtau(ir)=abs(tau(ir+1)-tau(ir))
	enddo
	do ir=1,nr
		exptau(ir)=exp(-dtau(ir))
	enddo

	Im=0d0
	do ir=nr,1,-1
		if(ir.eq.nr) then
			s0=0d0
		else
			s0=Si(ir+1)
		endif
		s1=Si(ir)
		x1=dtau(ir)
		if(x1.lt.1d-4) then
			y=0.5d0*x1*(s0+s1)
		else if(x1.gt.1d4) then
			y=s1+(s0-s1)/x1
		else
			y=(-(s0*x1-s1+s0)*exptau(ir)+s1*x1-s1+s0)/x1
		endif
		if(ir.eq.nr) then
			Im(ir)=y
		else
			Im(ir)=Im(ir+1)*exptau(ir)+y
		endif
	enddo

	Ip=0d0
	do ir=1,nr-1
		s0=Si(ir)
		s1=Si(ir+1)
		x1=dtau(ir)
		if(x1.lt.1d-4) then
			y=0.5d0*x1*(s0+s1)
		else if(x1.gt.1d4) then
			y=s1+(s0-s1)/x1
		else
			y=(-(s0*x1-s1+s0)*exptau(ir)+s1*x1-s1+s0)/x1
		endif
		Ip(ir+1)=Ip(ir)*exptau(ir)+y
	enddo

	Ij=(Ip+Im)/2d0
	Ih=(Ip-Im)/2d0

	return
	end
	


	subroutine InvertIjExp(tau,Linv,nr)
	IMPLICIT NONE
	integer ir,nr,jr
	real*8 Linv(nr,nr)
	real*8 Si(nr),tau(nr),Ip(nr),Im(nr)
	real*8 exptau(nr),dtau(nr),s0,s1,x1,y

	dtau(nr)=tau(nr)
	do ir=nr-1,1,-1
		dtau(ir)=abs(tau(ir+1)-tau(ir))
	enddo
	do ir=1,nr
		exptau(ir)=exp(-dtau(ir))
	enddo

	do jr=1,nr
		Im=0d0

		ir=jr
		s0=1d0
		s1=0d0
		x1=dtau(ir)
		if(x1.lt.1d-4) then
			y=0.5d0*x1*(s0+s1)
		else if(x1.gt.1d4) then
			y=s0+(s1-s0)/x1
		else
			y=(s0*x1+s1-s0-exptau(ir)*(s1*x1+s1-s0))/x1
		endif
		Im(jr)=y

		if(jr.gt.1) then
			ir=jr-1
			s0=0d0
			s1=1d0
			x1=dtau(ir)
			if(x1.lt.1d-4) then
				y=0.5d0*x1*(s0+s1)
			else if(x1.gt.1d4) then
				y=s0+(s1-s0)/x1
			else
				y=(s0*x1+s1-s0-exptau(ir)*(s1*x1+s1-s0))/x1
			endif
			Im(ir)=Im(jr)*exptau(ir)+y
			if(jr.gt.2) then
				do ir=jr-2,1,-1
					Im(ir)=Im(ir+1)*exptau(ir)
					if(Im(ir).le.0d0) exit
				enddo
			endif
		endif
		
		Ip=0d0
		if(jr.gt.1) then
			ir=jr-1
			s0=1d0
			s1=0d0
			x1=dtau(ir)
			if(x1.lt.1d-4) then
				y=0.5d0*x1*(s0+s1)
			else if(x1.gt.1d4) then
				y=s0+(s1-s0)/x1
			else
				y=(s0*x1+s1-s0-exptau(ir)*(s1*x1+s1-s0))/x1
			endif
			Ip(jr)=y
			if(jr.lt.nr) then
				ir=jr
				s0=0d0
				s1=1d0
				x1=dtau(ir)
				if(x1.lt.1d-4) then
					y=0.5d0*x1*(s0+s1)
				else if(x1.gt.1d4) then
					y=s0+(s1-s0)/x1
				else
					y=(s0*x1+s1-s0-exptau(ir)*(s1*x1+s1-s0))/x1
				endif
				Ip(jr+1)=Ip(jr)*exptau(ir)+y
				if(jr.lt.nr-1) then
					do ir=jr+2,nr
						Ip(ir)=Ip(ir-1)*exptau(ir-1)
						if(Ip(ir).le.0d0) exit
					enddo
				endif
			endif
		endif
		Linv(1:nr,jr)=(Ip(1:nr)+Im(1:nr))/2d0
	enddo

	return
	end


	
      SUBROUTINE tridag(a,b,c,r,u,n,info)
      INTEGER n
      REAL*8 a(n),b(n),c(n),r(n),u(n)
      INTEGER j,info
      REAL*8 bet,gam(n)
      info=0
      if(b(1).eq.0d0) then
      	info=1
      	return
      endif
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0d0) then
	      	info=1
   	   		return
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END


	subroutine spldiff(x,y,n,u)
	IMPLICIT NONE
c Calculates the numerical derivative using spline interpolation
c Calls the Numerical Recipes routine TRIDIAG (double precision)
c REMARK
c End point behaviour is not so good
c INPUTS
c x(n),y(n) - arrays comprising x & y data (equispaced x)
c h - spacing in x values
c n - array size
c OUTPUT
c u(n) - first order derivative (dy/dx)
c
c P. Arumugam, IIT Roorkee - 07-June-2010
c
	integer n,i,info
	real*8 x(n),y(n)
	real*8 a(n),b(n),c(n),r(n),u(n)
	b(1)=1.d0
	c(1)=1.d0
	r(1)=2.d0*(y(2)-y(1))/(x(2)-x(1))
	a(n)=1.d0
	b(n)=1.d0
	r(n)=2.d0*(y(n)-y(n-1))/(x(n)-x(n-1))

	a(n)=0.d0
	b(n)=1.d0
	r(n)=y(n)
	do i=2,n-1
		a(i)=1.d0
		b(i)=4.d0
		c(i)=1.d0
		r(i)=6.d0*(y(i+1)-y(i-1))/(x(i+1)-x(i-1))
	end do
	call tridag(a,b,c,r,u,n,info)
	return
	end
	
	
	subroutine splintegral(x,dy,n,y)
	IMPLICIT NONE
	integer n,i,info
	real*8 x(n),y(n),dy(n)
	real*8 a(n),b(n),c(n),dx(n),r(n)
	
	do i=1,n-1
		dx(i)=x(i+1)-x(i)
	enddo
	dx(n)=dx(n-1)
	r=dy
	do i=2,n-1
		a(i)=-dx(i)/(dx(i-1)*(dx(i-1)+dx(i)))
		b(i)=dx(i)/(dx(i-1)*(dx(i-1)+dx(i)))-dx(i-1)/(dx(i)*(dx(i-1)+dx(i)))
		c(i)=dx(i-1)/(dx(i)*(dx(i-1)+dx(i)))
	enddo
	a(1)=0d0
	b(1)=-1d0/dx(1)
	c(1)=1d0/dx(1)
	a(n)=0d0
	b(n)=1d0
	c(n)=0d0
	r(n)=0d0
	call tridag(a,b,c,r,y,n,info)
	return
	end
		
	
      SUBROUTINE banmul(a,n,m1,m2,np,mp,x,b)
      INTEGER m1,m2,mp,n,np
      REAL*8 a(np,mp),b(n),x(n)
      INTEGER i,j,k
      do 12 i=1,n
        b(i)=0d0
        k=i-m1-1
        do 11 j=max(1,1-k),min(m1+m2+1,n-k)
          b(i)=b(i)+a(i,j)*x(j+k)
11      continue
12    continue
      return
      END


      SUBROUTINE banbks(a,n,m1,m2,np,mp,al,mpl,indx,b)
      INTEGER m1,m2,mp,mpl,n,np,indx(n)
      REAL*8 a(np,mp),al(np,mpl),b(n)
      INTEGER i,k,l,mm
      REAL*8 dum
      mm=m1+m2+1
      if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in banbks'
      l=m1
      do 12 k=1,n
        i=indx(k)
        if(i.ne.k)then
          dum=b(k)
          b(k)=b(i)
          b(i)=dum
        endif
        if(l.lt.n)l=l+1
        do 11 i=k+1,l
          b(i)=b(i)-al(k,i-k)*b(k)
11      continue
12    continue
      l=1
      do 14 i=n,1,-1
        dum=b(i)
        do 13 k=2,l
          dum=dum-a(i,k)*b(k+i-1)
13      continue
        b(i)=dum/a(i,1)
        if(l.lt.mm) l=l+1
14    continue
      return
      END

      SUBROUTINE bandec(a,n,m1,m2,np,mp,al,mpl,indx,d)
      INTEGER m1,m2,mp,mpl,n,np,indx(n)
      REAL*8 d,a(np,mp),al(np,mpl),TINY
      PARAMETER (TINY=1d-30)
      INTEGER i,j,k,l,mm
      REAL*8 dum
      mm=m1+m2+1
      if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in bandec'
      l=m1
      do 13 i=1,m1
        do 11 j=m1+2-i,mm
          a(i,j-l)=a(i,j)
11      continue
        l=l-1
        do 12 j=mm-l,mm
          a(i,j)=0d0
12      continue
13    continue
      d=1.
      l=m1
      do 18 k=1,n
        dum=a(k,1)
        i=k
        if(l.lt.n)l=l+1
        do 14 j=k+1,l
          if(abs(a(j,1)).gt.abs(dum))then
            dum=a(j,1)
            i=j
          endif
14      continue
        indx(k)=i
        if(dum.eq.0.) a(k,1)=TINY
        if(i.ne.k)then
          d=-d
          do 15 j=1,mm
            dum=a(k,j)
            a(k,j)=a(i,j)
            a(i,j)=dum
15        continue
        endif
        do 17 i=k+1,l
          dum=a(i,1)/a(k,1)
          al(k,i-k)=dum
          do 16 j=2,mm
            a(i,j-1)=a(i,j)-dum*a(k,j)
16        continue
          a(i,mm)=0d0
17      continue
18    continue
      return
      END



	real*8 function ComputeSolidAbun(T0,ir)
	use GlobalSetup
	use Constants
	use CloudModule
	use AtomsModule
	IMPLICIT NONE
	real*8 T0,densv(nCS),fmax,f
	integer ir,iCS,iCS0
	
	densv=(mu*mp/(kb*T0))*exp(BTP-ATP/T0)
	fmax=-1d200
	iCS0=1
	do iCS=1,nCS
		if(xv_bot(iCS).gt.0d0) then
			if(T0.gt.maxT(iCS)) densv(iCS)=densv(iCS)+(mu(iCS)*mp/(kb*T0*10d0))*exp(BTP(iCS)-ATP(iCS)/(T0*10d0))
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





	subroutine PosSolve(A,x,minx,maxx,N,IP,WS)
	IMPLICIT NONE
	integer N,NN,i
	real*8 A(N,N),x(N),W(3*N,N+1),minx(N),maxx(N)

	integer ME,MA,MG,MODE,MDW
	real*8 PRGOPT(10),RNORME,RNORML
	integer IP(*)
	real*8 WS(*)

	NN=5*N
	
	ME=0
	MA=N
	MG=2*N
	MDW=3*N

	PRGOPT(1)=4
	PRGOPT(2)=1
	PRGOPT(3)=1
	PRGOPT(4)=7
	PRGOPT(5)=10
	PRGOPT(6)=1
	PRGOPT(7)=1

	do i=1,N
		W(1:N,i)=A(1:N,i)/x(i)
	enddo
	W(1:N,N+1)=1d0

	W(N+1:2*N,1:N)=0d0
	W(N+1:2*N,N+1)=minx(1:N)
	do i=1,N
		W(i+N,i)=1d0
	enddo

	W(2*N+1:3*N-1,1:N+1)=0d0
	W(2*N+1:3*N,N+1)=-maxx(1:N)
	do i=1,N-1
		W(i+2*N,i)=-1d0
	enddo


	call dlsei(W, MDW, ME, MA, MG, N, PRGOPT, x, RNORME,
     +   RNORML, MODE, WS, IP)

	return
	end
	



	subroutine AddScatter(Si_in,tauR_in,Ca,Cs,Ce,SurfAlb,nr,nu,wnu,nnu,NRHS)
	use Constants
	IMPLICIT NONE
	integer inu,nnu,ilam,ir,info,NRHS,nr,i,j,jr
	real*8 tauR(nr),tauR_in(nr),Si_in(nr,NRHS)
	real*8 Si(nr,NRHS),Ca(nr),Cs(nr),Ce(nr)
	real*8 nu(nnu),wnu(nnu),albedo(1:nr),SurfAlb
	real*8 Ij(nr),Itot(nr,NRHS),Linv(nr,nr),Lmat(nr,nr),Hsurf(nr)
	integer IWORKomp(nr)

	do ir=1,nr
		albedo(ir)=Cs(ir)/Ce(ir)
		if(.not.albedo(ir).lt.1d0/(1d0+1d-4)) albedo(ir)=1d0/(1d0+1d-4)
	enddo

	Linv=0d0
	Hsurf=0d0
	do inu=1,nnu
		tauR(1:nr)=tauR_in(1:nr)/abs(nu(inu))
		call InvertIjExp(tauR,Lmat,nr)
		do ir=1,nr
			Linv(ir,1:nr)=Linv(ir,1:nr)+wnu(inu)*Lmat(ir,1:nr)*albedo(ir)
		enddo
		Hsurf(1:nr)=Hsurf(1:nr)+2d0*nu(inu)*wnu(inu)*Lmat(1,1:nr)
	enddo
	do inu=1,nnu
		tauR(1:nr)=abs((tauR_in(1:nr)-tauR_in(1))/nu(inu))
		do ir=1,nr
			Linv(ir,1:nr)=Linv(ir,1:nr)+wnu(inu)*SurfAlb*Hsurf(1:nr)*albedo(ir)*exp(-tauR(ir))
		enddo
	enddo

	if(.false.) then
		Itot=Si_in
		i=0
1		continue
		Si=0d0
		do j=1,NRHS
			Si(1:nr,j)=matmul(Linv,Itot(1:nr,j))
		enddo
		Si_in=Si_in+Si
		Itot=Si
		i=i+1
		if(maxval(Si/Si_in).gt.1d-3) goto 1
		return
	endif
	
	Linv=-Linv
	do ir=1,nr
		Linv(ir,ir)=1d0+Linv(ir,ir)
	enddo
	Itot=0d0
	Itot(1:nr,1:NRHS)=Si_in(1:nr,1:NRHS)
	info=0
	call DGESV( nr, NRHS, Linv, nr, IWORKomp, Itot, nr, info)
	if(info.ne.0) then
		print*,"problem in scattering",info
		Si=Si_in
		return
	endif

	Si(1:nr,1:NRHS)=Itot(1:nr,1:NRHS)
	do i=1,NRHS
	do ir=1,nr
		if(.not.Si(ir,i).ge.Si_in(ir,i)) then
			Si(ir,i)=Si_in(ir,i)
		endif
	enddo
	enddo
	Si_in=Si

	return
	end


