	subroutine DoComputeT(converged,f)
	use GlobalSetup
	use Constants
	use CloudModule
	IMPLICIT NONE
	integer iphase,Titer
	real*8 tau_V,tau_T,Planck,CrV(nr),CrT(nr),CrVe(nr),f
	real*8 Cs,g,dlnT,dlnP
	real*8,allocatable :: Ce(:,:,:),Ca(:,:,:)
	real*8 tot,tot2,tot3,tot4,chi2,must,gamma,dP,Tirr,T0(nr),must_i
	integer ir,ilam,ig,i,iT,niter
	logical docloud0(max(nclouds,1)),converged

	if(doMCcompute.and.nTiter.gt.0) then
		call MCDoComputeT(converged,f)
		converged=.false.
	else
		call DoComputeTeddington(converged,f)
		return
		call DoComputeTfluxes(converged,f)
	endif

	return
	end


	subroutine DoComputeTeddington(converged,f)
	use GlobalSetup
	use Constants
	use CloudModule
	IMPLICIT NONE
	integer iphase,iter
	real*8 tau_V,tau_T,Planck,f
	real*8 g,dlnT,dlnP,d,tau,tautot,fact,contr,tau_a,exp_tau
	real*8,allocatable :: Ce(:,:,:),Ca(:,:,:),Cs(:,:,:),taustar(:,:),tauR_nu(:,:,:)
	real*8 tot,tot2,tot3,tot4,chi2,must,gamma,dP,Tirr,T0(nr),T1(nr),must_i,E,E0,Tinp(nr)
	real*8 z,Hstar(0:nr),Jtot,Htot,Ktot
	real*8 fedd(nr),Hedd(0:nr),lH1,lH2,P1,P2,directHstar(0:nr)
	real*8,allocatable :: Jnu(:,:,:),Hnu(:,:,:),Knu(:,:,:)
	integer ir,ilam,ig,i,iT,niter,inu,nnu,jr,iTmin,iTmax
	logical docloud0(max(nclouds,1)),converged,stopscat
	real*8 tauf(nr),Si(0:nr),B1,B2,x1,x2,dx1,dx2,ax,bx,ff,TT,maxErr
	integer info,IWORK(10*(nr+1)*(nr+1)),NRHS,ii(3),iscat,nscat
	real*8 tau1,tau2,ee0,ee1,ee2,tauR(0:nr),Ij(0:nr),Ih(0:nr),scale,dtauR(0:nr)
	integer nlam_LR
	real*8 specres_LR,IntH(nr,nr),Fl(nr),Ts(nr),minFl(nr),maxFl(nr),maxfact
	real*8,allocatable :: lam_LR(:),dfreq_LR(:),freq_LR(:),BB_LR(:,:),IntHnu(:,:,:),dtauR_nu(:,:,:)
	integer i1,i2,ngF,j
	real*8 ww,w1,w2,SurfStar,SurfStar_omp,FstarBottom,tauRoss
	real*8,allocatable :: Si_omp(:,:),Ih_omp(:),Ij_omp(:),tauR_omp(:),Hsurf(:,:),Hstar_omp(:)
	real*8,allocatable :: temp_a(:),wtemp(:),Ca_HR(:,:),Cs_HR(:,:),Fstar_LR(:),SurfEmis_LR(:)
	integer,allocatable :: IP(:)
	real*8,allocatable :: WS(:),nu(:),wnu(:),Hstar_lam(:),Hsurf_lam(:),directHstar_omp(:)
	real*8,allocatable :: x_SIj(:),y_SIj(:),tauR_SIj(:),Ma_SIj(:),Mb_SIj(:),Mc_SIj(:)
	character*500 file
	integer icloud

	maxfact=4d0

	allocate(IP(nr*4+2))
	IP(1)=(2*(nr)+nr*3+(nr*2+2)*(nr+7))
	IP(2)=nr*4+2
	allocate(WS(IP(1)))
	
	specres_LR=min(specres/1.5,11d0)
	tot=lam(1)
	nlam_LR=1
	do while(tot.lt.lam(nlam))
		tot=tot*(1d0+1d0/specres_LR)
		nlam_LR=nlam_LR+1
	enddo
	allocate(lam_LR(nlam_LR))
	allocate(freq_LR(nlam_LR))
	allocate(dfreq_LR(nlam_LR))
	lam_LR(1)=lam(1)
	do i=2,nlam_LR-1
		lam_LR(i)=lam_LR(i-1)*(1d0+1d0/specres_LR)
	enddo
	lam_LR(nlam_LR)=lam(nlam)
	do i=1,nlam_LR
		freq_LR(i)=1d0/lam_LR(i)
	enddo
	do i=1,nlam_LR-1
		dfreq_LR(i)=abs(freq_LR(i+1)-freq_LR(i))
	enddo

	allocate(BB_LR(nBB,nlam_LR))

	BB_LR=0d0
	do j=nBB,1,-1
		TT=real(j)
		tot=0d0
		do i=1,nlam_LR-1
			BB_LR(j,i)=Planck(TT,freq_LR(i))
			tot=tot+dfreq_LR(i)*BB_LR(j,i)
		enddo
		scale=((2d0*(pi*kb*TT)**4)/(15d0*hplanck**3*clight**3))/tot
		if(scale.gt.0d0.and.tot.gt.0d0) then
			BB_LR(j,1:nlam_LR)=BB_LR(j,1:nlam_LR)*scale
		else if(j.lt.nBB) then
			scale=((2d0*(pi*kb*TT)**4)/(15d0*hplanck**3*clight**3))/((2d0*(pi*kb*real(j+1))**4)/(15d0*hplanck**3*clight**3))
			BB_LR(j,1:nlam_LR)=BB_LR(j+1,1:nlam_LR)*scale
		endif
	enddo


	allocate(taustar(nlam_LR,ng))
	allocate(Ca_HR(nlam,ng))
	allocate(Cs_HR(nlam,ng))
	allocate(Ce(nr,nlam_LR,ng))
	allocate(Ca(nr,nlam_LR,ng))
	allocate(Cs(nr,nlam_LR,ng))
	allocate(Fstar_LR(nlam_LR))
	allocate(tauR_nu(0:nr,nlam_LR,ng))
	allocate(dtauR_nu(0:nr,nlam_LR,ng))
	allocate(IntHnu(nlam_LR,0:nr,0:nr))
	allocate(SurfEmis_LR(nlam_LR))

	call ComputeSurface()
	do ilam=1,nlam_LR-1
		i1=0
		i2=0
		do i=1,nlam-1
			if(lam_LR(ilam).ge.lam(i).and.lam_LR(ilam).lt.lam(i+1)) i1=i
			if(lam_LR(ilam+1).ge.lam(i).and.lam_LR(ilam+1).lt.lam(i+1)) i2=i+1
		enddo
		if(ilam.eq.1) i1=1
		if(ilam.eq.nlam_LR-1) i2=nlam
		if(i1.gt.0.and.i2.gt.0) then
			tot=0d0
			SurfEmis_LR(ilam)=0d0
			do i=i1,i2
				if(i1.eq.i2) then
					ww=1d0
				else if(i.eq.i1) then
					ww=abs(lam_LR(ilam)-lam(i+1))
				else if(i.eq.i2) then
					ww=abs(lam_LR(ilam+1)-lam(i))
				else
					ww=abs(lam(i)-lam(i+1))
				endif
				SurfEmis_LR(ilam)=SurfEmis_LR(ilam)+ww*surface_emis(i)
				tot=tot+ww
			enddo
			SurfEmis_LR(ilam)=SurfEmis_LR(ilam)/tot
		else
			SurfEmis_LR(ilam)=1d0
		endif
	enddo


	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo
	
	niter=50
	nscat=10
	epsiter=3d-2

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

	do ir=nr,1,-1
		if(T(ir).gt.0d0) then
			Tinp(ir)=T(ir)
		else if(ir.lt.nr) then
			Tinp(ir)=T(ir+1)
		else
			Tinp(ir)=must*sqrt(Rstar/(2d0*Dplanet))*Tstar
		endif
		T(ir)=Tinp(ir)
	enddo

	allocate(temp_a(ng*nlam))
	allocate(wtemp(ng*nlam))
	Fstar_LR=0d0
	Cs=0d0
	Ca=0d0
	do ir=1,nr
		do ilam=1,nlam-1
			g=0d0
			do ig=1,ng
				call Crossections(ir,ilam,ig,Ca_HR(ilam,ig),Cs_HR(ilam,ig),docloud0)
				Ca_HR(ilam,ig)=Ca_HR(ilam,ig)/dens(ir)
				Cs_HR(ilam,ig)=Cs_HR(ilam,ig)*(1d0-g)/dens(ir)
			enddo
		enddo

		do ilam=1,nlam_LR-1
			i1=0
			i2=0
			do i=1,nlam-1
				if(lam_LR(ilam).ge.lam(i).and.lam_LR(ilam).lt.lam(i+1)) i1=i
				if(lam_LR(ilam+1).ge.lam(i).and.lam_LR(ilam+1).lt.lam(i+1)) i2=i+1
			enddo
			if(ilam.eq.1) i1=1
			if(ilam.eq.nlam_LR-1) i2=nlam
			if(i1.gt.0.and.i2.gt.0) then
				ngF=0
				tot=0d0
				Fstar_LR(ilam)=0d0
				Cs(ir,ilam,1:ng)=0d0
				do i=i1,i2
					if(i1.eq.i2) then
						ww=1d0
					else if(i.eq.i1) then
						ww=abs(lam_LR(ilam)-lam(i+1))
					else if(i.eq.i2) then
						ww=abs(lam_LR(ilam+1)-lam(i))
					else
						ww=abs(lam(i)-lam(i+1))
					endif
					do ig=1,ng
						ngF=ngF+1
						temp_a(ngF)=Ca_HR(i,ig)
						wtemp(ngF)=ww*wgg(ig)
					enddo
					Fstar_LR(ilam)=Fstar_LR(ilam)+ww*Fstar(i)
					Cs(ir,ilam,1:ng)=Cs(ir,ilam,1:ng)+ww*Cs_HR(i,1:ng)
					tot=tot+ww
				enddo
				Fstar_LR(ilam)=Fstar_LR(ilam)/tot
				Cs(ir,ilam,1:ng)=Cs(ir,ilam,1:ng)/tot
				tot=0d0
				do ig=1,ngF
					tot=tot+temp_a(ig)*wtemp(ig)
				enddo
				tot=tot/sum(wtemp(1:ngF))

				call sortw(temp_a,wtemp,ngF)
				if(ng.eq.1) then
					tot=0d0
					do ig=1,ngF
						tot=tot+temp_a(ig)*wtemp(ig)
					enddo
					tot=tot/sum(wtemp(1:ngF))
					Ca(ir,ilam,1)=tot
				else
					do ig=2,ngF
						wtemp(ig)=wtemp(ig)+wtemp(ig-1)
					enddo
					wtemp(1:ngF)=wtemp(1:ngF)/wtemp(ngF)
					do ig=1,ng
						call hunt(wtemp,ngF,gg(ig),j)
						if(j.eq.0) then
							Ca(ir,ilam,ig)=temp_a(1)
						else
							w1=(gg(ig)-wtemp(j+1))/(wtemp(j)-wtemp(j+1))
							Ca(ir,ilam,ig)=temp_a(j)*w1+temp_a(j+1)*(1d0-w1)
						endif
					enddo
					tot2=0d0
					do ig=1,ng
						tot2=tot2+wgg(ig)*Ca(ir,ilam,ig)
					enddo
					if(tot2.gt.0d0) then
						Ca(ir,ilam,1:ng)=Ca(ir,ilam,1:ng)*tot/tot2
					else
						Ca(ir,ilam,1:ng)=tot
					endif
				endif
			else
				Ca(ir,ilam,1:ng)=0d0
				Cs(ir,ilam,1:ng)=0d0
				Fstar_LR(ilam)=0d0
			endif
		enddo
	enddo
	deallocate(temp_a)
	deallocate(wtemp)

	Ce=Ca+Cs
	do ir=1,nr
		do ilam=1,nlam_LR-1
			do ig=1,ng
				if(Ca(ir,ilam,ig)/Ce(ir,ilam,ig).lt.1d-4) then
					Ca(ir,ilam,ig)=Cs(ir,ilam,ig)/(1d4-1d0)
					Ce(ir,ilam,ig)=Ca(ir,ilam,ig)+Cs(ir,ilam,ig)
				endif
			enddo
		enddo
	enddo

	converged=.true.

	do ilam=1,nlam_LR-1
		do ig=1,ng
			do jr=nr,0,-1
				if(jr.eq.nr) then
					ir=jr
					d=0.5*P(ir)*1d6/grav(ir)
					tau=d*Ce(ir,ilam,ig)
				else if(jr.eq.nr-1) then
					ir=jr
					d=abs(0.5*P(ir+1)-P(ir+1))*1d6/grav(ir)
					tau=d*Ce(ir+1,ilam,ig)
					d=abs(sqrt(P(ir+1)*P(ir))-P(ir+1))*1d6/grav(ir)
					tau=tau+d*Ce(ir,ilam,ig)
				else if(jr.eq.0) then
					ir=1
					d=abs(sqrt(P(ir+1)*P(ir))-P(ir))*1d6/grav(ir)
					tau=d*Ce(ir,ilam,ig)
					d=abs(P(ir)*sqrt(P(ir)/P(ir+1))-P(ir))*1d6/grav(ir)
					tau=tau+d*Ce(ir,ilam,ig)
				else
					ir=jr
					d=abs(sqrt(P(ir+2)*P(ir+1))-P(ir+1))*1d6/grav(ir)
					tau=d*Ce(ir+1,ilam,ig)
					d=abs(sqrt(P(ir+1)*P(ir))-P(ir+1))*1d6/grav(ir)
					tau=tau+d*Ce(ir,ilam,ig)
				endif
				if(P(ir).gt.Psimplecloud) then
					tau=tau+1d4
				endif
				if(.not.tau.gt.1d-8) then
					tau=1d-8
				endif
				if(tau.gt.1d8) then
					tau=1d8
				endif
				if(ir.lt.nr) then
					tauR(jr)=tauR(jr+1)+tau
				else
					tauR(jr)=tau
				endif
				dtauR_nu(jr,ilam,ig)=tau
			enddo
			tauR_nu(0:nr,ilam,ig)=tauR(0:nr)
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
	do ilam=1,nlam_LR-1
c		Fstar_LR(ilam)=Planck(Tstar,freq_LR(ilam))*pi*Rstar**2
		E=E+dfreq_LR(ilam)*Fstar_LR(ilam)
	enddo
	scale=pi*Rstar**2*((2d0*(pi*kb*Tstar)**4)/(15d0*hplanck**3*clight**3))/E
	Fstar_LR=Fstar_LR*scale
	
	ff=1d0
	
	Hstar(0:nr)=0d0
	directHstar(0:nr)=0d0
	SurfStar=0d0
	IntHnu(1:nlam_LR,0:nr,0:nr)=0d0
	
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(x_SIj,y_SIj,tauR_SIj,Ma_SIj,Mb_SIj,Mc_SIj,Si_omp,tauR_omp,Ih_omp,Ij_omp,ilam,ig,ir,inu,jr,
!$OMP&			Hstar_omp,SurfStar_omp,contr,FstarBottom,Hstar_lam,Hsurf_lam,directHstar_omp)
!$OMP& SHARED(nlam_LR,ng,nr,nnu,tauR_nu,nu,wnu,dfreq_LR,wgg,IntHnu,SurfEmis_LR,dtauR_nu,Ca,Ce,Cs,Hsurf,
!$OMP&			Hstar,SurfStar,Dplanet,Fstar_LR,must,directHstar)
	allocate(x_SIj(nr+2),y_SIj(nr+2),tauR_SIj(0:nr+1),Ma_SIj(nr+2),Mb_SIj(nr+2),Mc_SIj(nr+2))
	allocate(Si_omp(0:nr,0:nr+1),tauR_omp(0:nr),Ih_omp(0:nr),Ij_omp(0:nr))
	allocate(Hstar_omp(0:nr),Hstar_lam(0:nr),Hsurf_lam(0:nr),directHstar_omp(0:nr))
	Hstar_omp=0d0
	SurfStar_omp=0d0
	directHstar_omp=0d0
!$OMP DO
	do ilam=1,nlam_LR-1
		call tellertje(ilam,nlam_LR-1)
		do ig=1,ng
			Hstar_lam(0:nr)=0d0
			Hsurf_lam(0:nr)=0d0
			contr=(Fstar_LR(ilam)/(pi*Dplanet**2))
			tauR_omp(0:nr)=tauR_nu(0:nr,ilam,ig)/abs(must)
			Ij_omp(0:nr)=contr*exp(-tauR_omp(0:nr))/(4d0*pi)
			Ih_omp(0:nr-1)=-must*contr*exp(-(tauR_omp(0:nr-1)+tauR_omp(1:nr))/2d0)
			Ih_omp(nr)=-must*contr*exp(-tauR_omp(nr))

			Hstar_lam(0:nr)=Hstar_lam(0:nr)+dfreq_LR(ilam)*wgg(ig)*Ih_omp(0:nr)
			directHstar_omp(0:nr)=directHstar_omp(0:nr)+dfreq_LR(ilam)*wgg(ig)*Ih_omp(0:nr)
			SurfStar_omp=SurfStar_omp+dfreq_LR(ilam)*wgg(ig)*abs(Ih_omp(0))*SurfEmis_LR(ilam)
			FstarBottom=abs(Ih_omp(0))

c Si_omp(0:nr,0) is the direct stellar contribution
			Si_omp(1:nr,0)=(Ij_omp(1:nr)*Cs(1:nr,ilam,ig)/Ce(1:nr,ilam,ig))
			Si_omp(0,0)=Ij_omp(0)*(1d0-SurfEmis_LR(ilam))

c Si_omp(0:nr,1:nr) is the direct contribution from the atmosphere
			do ir=1,nr
				Si_omp(0:nr,ir)=0d0
				Si_omp(ir,ir)=Ca(ir,ilam,ig)/Ce(ir,ilam,ig)
			enddo

c Si_omp(0:nr,nr+1) is the direct contribution from the surface
			Si_omp(0:nr,nr+1)=0d0
			do inu=1,nnu
				tauR_omp(0:nr)=tauR_nu(0:nr,ilam,ig)/abs(nu(inu))
				Ij_omp(0:nr)=exp(-abs(tauR_omp(0:nr)-tauR_omp(0)))
				call ComputeDeriv(tauR_omp(0:nr),Ij_omp(0:nr),Ih_omp(0:nr),nr+1,Ij_omp(0),Ij_omp(nr))
				Ih_omp(0:nr)=Ij_omp(0:nr)
				Si_omp(0:nr,nr+1)=Si_omp(0:nr,nr+1)+wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ij_omp(0:nr)/(4d0*pi)
				IntHnu(ilam,0:nr,0)=IntHnu(ilam,0:nr,0)+2d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(0:nr)
				Hsurf_lam(0:nr)=Hsurf_lam(0:nr)+FstarBottom*2d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*
     & 								Ih_omp(0:nr)*(1d0-SurfEmis_LR(ilam))
			enddo
			Si_omp(0,nr+1)=0d0
			Si_omp(1:nr,nr+1)=Si_omp(1:nr,nr+1)*Cs(1:nr,ilam,ig)/Ce(1:nr,ilam,ig)

			call AddScatter(Si_omp(1:nr,0:nr+1),tauR_nu(1:nr,ilam,ig),
     &					Ca(1:nr,ilam,ig),Cs(1:nr,ilam,ig),Ce(1:nr,ilam,ig),nr,nu,wnu,nnu,nr+2)

			do inu=1,nnu
				tauR_omp(0:nr)=tauR_nu(0:nr,ilam,ig)/abs(nu(inu))
				call SolveIj(tauR_omp(0:nr),Si_omp(0:nr,0),Ij_omp(0:nr),nr,x_SIj,y_SIj,tauR_SIj(0:nr+1),Ma_SIj,Mb_SIj,Mc_SIj)
				call ComputeDeriv(tauR_omp(0:nr),Ij_omp(0:nr),Ih_omp(0:nr),nr+1,-Ij_omp(0),Ij_omp(nr))
				Hstar_lam(0:nr)=Hstar_lam(0:nr)+8d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(0:nr)
				SurfStar_omp=SurfStar_omp+8d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ij_omp(1)*SurfEmis_LR(ilam)
			enddo

			do ir=1,nr
				do inu=1,nnu
					tauR_omp(0:nr)=tauR_nu(0:nr,ilam,ig)/abs(nu(inu))
					call SolveIj(tauR_omp(0:nr),Si_omp(0:nr,ir),Ij_omp(0:nr),nr,x_SIj,y_SIj,tauR_SIj(0:nr+1),Ma_SIj,Mb_SIj,Mc_SIj)
					call ComputeDeriv(tauR_omp(0:nr),Ij_omp(0:nr),Ih_omp(0:nr),nr+1,-Ij_omp(0),Ij_omp(nr))
					IntHnu(ilam,1:nr,ir)=IntHnu(ilam,1:nr,ir)+8d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
				enddo
			enddo

			do inu=1,nnu
				tauR_omp(0:nr)=tauR_nu(0:nr,ilam,ig)/abs(nu(inu))
				call SolveIj(tauR_omp(0:nr),Si_omp(0:nr,nr+1),Ij_omp(0:nr),nr,x_SIj,y_SIj,tauR_SIj(0:nr+1),Ma_SIj,Mb_SIj,Mc_SIj)
				call ComputeDeriv(tauR_omp(0:nr),Ij_omp(0:nr),Ih_omp(0:nr),nr+1,-Ij_omp(0),Ij_omp(nr))
				IntHnu(ilam,0:nr,0)=IntHnu(ilam,0:nr,0)+8d0*nu(inu)*wnu(inu)*Ih_omp(0:nr)
				Hsurf_lam(0:nr)=Hsurf_lam(0:nr)+FstarBottom*8d0*nu(inu)*wnu(inu)*Ih_omp(0:nr)*(1d0-SurfEmis_LR(ilam))
			enddo

			do ir=0,nr
				Hstar_omp(ir)=Hstar_omp(ir)+min(min(Hstar_lam(ir),0d0)+max(Hsurf_lam(ir),0d0),0d0)
			enddo
		enddo
	enddo
!$OMP END DO
!$OMP CRITICAL
	Hstar(0:nr)=Hstar(0:nr)+Hstar_omp(0:nr)
	directHstar(0:nr)=directHstar(0:nr)+directHstar_omp(0:nr)
	SurfStar=SurfStar+SurfStar_omp
!$OMP END CRITICAL
	deallocate(x_SIj,y_SIj,tauR_SIj,Ma_SIj,Mb_SIj,Mc_SIj)
	deallocate(Si_omp,tauR_omp,Ih_omp,Ij_omp)
	deallocate(Hstar_omp,Hstar_lam,Hsurf_lam,directHstar_omp)
!$OMP FLUSH
!$OMP END PARALLEL


	do iter=1,niter

	Ts(1:nr)=T(1:nr)

	Hedd=0d0
	E0=4d0*(((pi*kb*TeffP)**4)/(15d0*hplanck**3*clight**3))
	do ir=0,nr
		Hedd(ir)=Hedd(ir)+E0+abs(Hstar(ir))
	enddo

	IntH=0d0
	do ir=1,nr
		do ilam=1,nlam_LR-1
			do jr=1,nr
				iT=T(jr)+1
				if(iT.gt.nBB-1) iT=nBB-1
				if(iT.lt.1) iT=1
				scale=(T(jr)/real(iT))**4
				IntH(ir,jr)=IntH(ir,jr)+scale*BB_LR(iT,ilam)*IntHnu(ilam,ir,jr)
			enddo
		enddo
	enddo

c=========== begin experimental redistribution ===========================================
	if(deepredist) then

	E0=(Rstar/Dplanet)**2*((2d0*(pi*kb*Tstar)**4)/(15d0*hplanck**3*clight**3))*(f_deepredist-must)
	do ir=1,nr
		iT=T(ir)+1
		if(iT.gt.nBB-1) iT=nBB-1
		if(iT.lt.1) iT=1
		tauRoss=0d0
		E=0d0
		do ilam=1,nlam_LR-1
			do ig=1,ng
				tauRoss=tauRoss+dfreq_LR(ilam)*wgg(ig)*BB_LR(iT,ilam)/max(tauR_nu(ir,ilam,ig),1d-6)
				E=E+dfreq_LR(ilam)*wgg(ig)*BB_LR(iT,ilam)
			enddo
		enddo
		tauRoss=E/tauRoss
		Hedd(ir)=Hedd(ir)+max(-abs(Hstar(ir)),E0*exp(-tauRoss))
	enddo

	endif
c=========== end experimental redistribution =============================================


	Fl=0d0
	iT=Tsurface+1
	if(iT.gt.nBB-1) iT=nBB-1
	if(iT.lt.1) iT=1
	scale=(Tsurface/real(iT))**4
	do ir=1,nr
		do ilam=1,nlam_LR-1
			Fl(ir)=Fl(ir)-abs(IntHnu(ilam,ir,0))*scale*BB_LR(iT,ilam)*SurfEmis_LR(ilam)
		enddo
	enddo
	do ir=1,nr
		Fl(ir)=Fl(ir)+Hedd(ir)
	enddo
	
	minFl=0.5
	maxFl=2.0

	IP(1)=(2*(nr)+nr*3+(nr*2+2)*(nr+7))
	IP(2)=nr*4+2
	NRHS=1

	call DGESV( nr, NRHS, IntH, nr, IWORK, Fl, nr, info )
c	call PosSolve(IntH,Fl,minFl,maxFl,nr,IP,WS)

	do ir=1,nr
		Fl(ir)=min(max(minFl(ir),Fl(ir)),maxFl(ir))
	enddo

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ir,E,E0,ilam,iTmin,iTmax,ig,iT)
!$OMP& SHARED(nr,nlam_LR,dfreq_LR,wgg,BB_LR,Ca,Fl,ng,T,Ts)
!$OMP DO
	do ir=1,nr
		iT=T(ir)+1
		if(iT.gt.nBB-1) iT=nBB-1
		if(iT.lt.1) iT=1
		E=0d0
		do ilam=1,nlam_LR-1
			do ig=1,ng
				E=E+dfreq_LR(ilam)*wgg(ig)*BB_LR(iT,ilam)*Ca(ir,ilam,ig)
			enddo
		enddo
		E=E*Fl(ir)
		iTmin=1
		iTmax=nBB
		iT=(iTmax+iTmin)/2
		do while(abs(iTmax-iTmin).gt.1)
			E0=0d0
			do ilam=1,nlam_LR-1
				do ig=1,ng
					E0=E0+dfreq_LR(ilam)*wgg(ig)*BB_LR(iT,ilam)*Ca(ir,ilam,ig)
				enddo
			enddo
			if(E0.gt.E) then
				iTmax=iT
			else
				iTmin=iT
			endif
			iT=0.5d0*(iTmax+iTmin)
			if(iT.lt.1) iT=1
			if(iT.gt.nBB) iT=nBB
		enddo
		Ts(ir)=real(iT)*(E/E0)**0.25
		if(.not.Ts(ir).gt.3d0) Ts(ir)=3d0
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	do ir=nr-1,1,-1
		if(P(ir).gt.Pplanet) exit
	enddo
	j=min(max(ir,2),nr-1)
	do ir=j-1,1,-1
		if(ir.lt.nr) then
			dlnP=log(P(ir+1)/P(ir))
			dlnT=log(Ts(ir+1)/Ts(ir))
			if((dlnT/dlnP).gt.nabla_ad(ir)) then
				dlnT=(nabla_ad(ir))*dlnP
				Ts(ir)=Ts(ir+1)/exp(dlnT)
			endif
			if((dlnT/dlnP).lt.-nabla_ad(ir)) then
				dlnT=(-nabla_ad(ir))*dlnP
				Ts(ir)=Ts(ir+1)/exp(dlnT)
			endif
		endif
	enddo

	do ir=j,nr
		dlnP=log(P(ir)/P(ir-1))
		dlnT=log(Ts(ir)/Ts(ir-1))
		if((dlnT/dlnP).gt.nabla_ad(ir)) then
			dlnT=(nabla_ad(ir))*dlnP
			Ts(ir)=Ts(ir-1)*exp(dlnT)
		endif
		if((dlnT/dlnP).lt.-nabla_ad(ir)) then
			dlnT=(-nabla_ad(ir))*dlnP
			Ts(ir)=Ts(ir-1)*exp(dlnT)
		endif
	enddo

	maxErr=0d0
	do ir=1,nr-1
		if(abs(T(ir)-Ts(ir))/(T(ir)+Ts(ir)).gt.maxErr) maxErr=abs(T(ir)-Ts(ir))/(T(ir)+Ts(ir))
	enddo

	do ir=1,nr
		T(ir)=(Ts(ir)+T(ir))/2d0
	enddo

	E0=(((pi*kb*TeffP)**4)/(15d0*hplanck**3*clight**3))
	E=SurfStar+E0
	do ilam=1,nlam_LR-1
		do ig=1,ng
			do inu=1,nnu
				tauR(0:nr)=(tauR_nu(1,ilam,ig)-tauR_nu(0:nr,ilam,ig))/abs(nu(inu))
				dtauR(0:nr)=dtauR_nu(0:nr,ilam,ig)/abs(nu(inu))
				do ir=1,nr
					iT=T(ir)+1
					if(iT.gt.nBB-1) iT=nBB-1
					if(iT.lt.1) iT=1
					scale=(T(ir)/real(iT))**4
					contr=nu(inu)*2d0*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*exp(-tauR(ir))
					E=E+contr*scale*BB_LR(iT,ilam)*(1d0-exp(-dtauR(ir)))*SurfEmis_LR(ilam)*Ca(ir,ilam,ig)/Ce(ir,ilam,ig)
					if(tauR(ir).gt.10d0) exit
				enddo
			enddo
		enddo
	enddo

	iTmin=1
	iTmax=nBB
	iT=(iTmax+iTmin)/2
	do while(abs(iTmax-iTmin).gt.1)
		E0=0d0
		do ilam=1,nlam_LR-1
			E0=E0+dfreq_LR(ilam)*BB_LR(iT,ilam)*SurfEmis_LR(ilam)
		enddo
		if(E0.gt.E) then
			iTmax=iT
		else
			iTmin=iT
		endif
		iT=0.5d0*(iTmax+iTmin)
		if(iT.lt.1) iT=1
		if(iT.gt.nBB) iT=nBB
	enddo
	Tsurface=real(iT)*(E/E0)**0.25

	if(.not.Tsurface.gt.3d0.or.E0.eq.0d0) Tsurface=3d0

	if(maxErr.lt.(epsiter/5d0).and.iter.gt.5) exit
	enddo

	call output("Surface temperature: " // dbl2string(Tsurface,'(f8.2)') // " K")

	maxErr=0d0
	do ir=1,nr-1
		if(abs(T(ir)-Tinp(ir))/(T(ir)+Tinp(ir)).gt.maxErr) maxErr=abs(T(ir)-Tinp(ir))/(T(ir)+Tinp(ir))
		if(maxErr.gt.epsiter) converged=.false.
	enddo
	call output("Maximum error on T-struct: " // dbl2string(maxErr*100d0,'(f5.1)') // "%")
	T0(1:nr)=Tinp(1:nr)
	T1(1:nr)=T(1:nr)
	do ir=1,nr
		T(ir)=f*T1(ir)+(1d0-f)*T0(ir)
	enddo

	call tellertje(niter,niter)
	call WriteStructure

	deallocate(lam_LR)
	deallocate(freq_LR)
	deallocate(dfreq_LR)
	deallocate(BB_LR)
	deallocate(taustar)
	deallocate(Ca_HR)
	deallocate(Cs_HR)
	deallocate(Ce)
	deallocate(Ca)
	deallocate(Cs)
	deallocate(Fstar_LR)
	deallocate(tauR_nu)
	deallocate(IntHnu)

	return
	end


	subroutine ComputeDeriv(x,y,dy,n,yp1,ypn)
	IMPLICIT NONE
	integer n,i
	real*8 x(n),y(n),dy(n),d2y(n),yp1,ypn
	real*8 A,B

	call spline(x,y,n,yp1,ypn,d2y)
	dy(1)=yp1
	dy(n)=ypn
	do i=2,n-1
		dy(i)=(y(i+1)-y(i))/(x(i+1)-x(i))+(x(i+1)-x(i))*(d2y(i)-d2y(i+1))/24d0
	enddo

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

	Ip=0d0
	do ir=1,nr-1
		s0=Si(ir)
		s1=Si(ir+1)
		x1=dtau(ir)
		y=(-(s0*x1-s1+s0)*exptau(ir)+s1*x1-s1+s0)/x1
		Ip(ir+1)=Ip(ir)*exptau(ir)+y
	enddo

	Im=0d0
	do ir=nr-1,1,-1
		s0=Si(ir+1)
		s1=Si(ir)
		x1=dtau(ir)
		y=(-(s0*x1-s1+s0)*exptau(ir)+s1*x1-s1+s0)/x1
		Im(ir)=Im(ir+1)*exptau(ir)+y
	enddo

	Ij=(Ip+Im)/2d0
	Ih=(Ip-Im)/2d0

	return
	end
	


	subroutine SolveIjExp2(tau,Si,Ij,Ih,nr)
	IMPLICIT NONE
	integer ir,nr,jr
	real*8 Ij(nr),Si(nr),tau(nr),Ip(nr),Im(nr),Ih(nr)
	real*8 exptau(nr),fact,exptau2,Q,dtau(nr)

	dtau(nr)=tau(nr)
	do ir=nr-1,1,-1
		dtau(ir)=abs(tau(ir+1)-tau(ir))
	enddo
	do ir=1,nr
		exptau(ir)=exp(-dtau(ir))
	enddo

	Ip=0d0
	do ir=1,nr-1
		Q=(Si(ir)*(1d0-(1d0+dtau(ir))*exptau(ir))+Si(ir+1)*(dtau(ir)-1d0+exptau(ir)))/dtau(ir)
		Ip(ir+1)=Ip(ir)*exptau(ir)+Q
	enddo

	Im=0d0
	do ir=nr,2,-1
		Q=(Si(ir)*(1d0-(1d0+dtau(ir))*exptau(ir))+Si(ir-1)*(dtau(ir)-1d0+exptau(ir)))/dtau(ir)
		Im(ir-1)=Im(ir)*exptau(ir)+Q
	enddo

	Ij=(Ip+Im)/2d0
	Ih=(Ip-Im)/2d0

	return
	end
	


	subroutine SolveIj(tauR_in,Si,Ij,nr,x,y,tauR,Ma,Mb,Mc)
	IMPLICIT NONE
	integer ir,nr,info
	real*8 tauR_in(0:nr),Ij(0:nr),Si(0:nr),fact,d
	real*8,allocatable :: MM(:,:),MMal(:,:)
	integer,allocatable :: indx(:)
	real*8 x(nr+2),y(nr+2),tauR(0:nr+1),Ma(nr+2),Mb(nr+2),Mc(nr+2)

	tauR(0:nr)=tauR_in(0:nr)
	tauR(nr+1)=0d0

	Ma=0d0
	Mb=0d0
	Mc=0d0
	ir=1
	do ir=1,nr
		fact=1d0/(0.5d0*(tauR(ir+1)+tauR(ir))-0.5d0*(tauR(ir)+tauR(ir-1)))
		Mb(ir+1)=1d0+fact*(1d0/(tauR(ir+1)-tauR(ir))+1d0/(tauR(ir)-tauR(ir-1)))
		Ma(ir+1)=-fact*1d0/(tauR(ir)-tauR(ir-1))
		Mc(ir+1)=-fact*1d0/(tauR(ir+1)-tauR(ir))
	enddo
	Mb(1)=1d0+1d0/(tauR(0)-tauR(1))
	Mc(1)=-1d0/(tauR(0)-tauR(1))
	Ma(nr+2)=1d0/(tauR(nr))
	Mb(nr+2)=-1d0-1d0/(tauR(nr))

	x=0d0
	x(2:nr+1)=Si(1:nr)
	info=0
	call tridag(Ma,Mb,Mc,x,y,nr+2,info)
	Ij(0:nr)=y(1:nr+1)
	do ir=0,nr
		if(Ij(ir).lt.0d0) info=1
	enddo
	if(info.ne.0) then
		allocate(MM(nr+2,3),MMal(nr+2,1))
		allocate(indx(nr+2))
		do ir=1,nr+2
			MM(ir,1)=Ma(ir)
			MM(ir,2)=Mb(ir)
			MM(ir,3)=Mc(ir)
		enddo
		call bandec(MM,nr+2,1,1,nr+2,3,MMal,1,indx,d)
		x=0d0
		x(2:nr+1)=Si(1:nr)
		call banbks(MM,nr+2,1,1,nr+2,3,MMal,1,indx,x)
		Ij(0:nr)=y(1:nr+1)
		deallocate(MM,MMal,indx)
	endif
	do ir=0,nr
		if(Ij(ir).lt.0d0) then
			Ij(ir)=0d0
		endif
	enddo

	return
	end
	

	subroutine SolveIjSurface(tauR_in,Ij,nr)
	IMPLICIT NONE
	integer ir,nr
	real*8 tauR_in(0:nr),Ij(0:nr),x(nr+2),y(nr+2),fact,d,tauR(0:nr+1)
	real*8 MM(nr+2,3),MMal(nr+2,1),Ma(nr+2),Mb(nr+2),Mc(nr+2)
	integer indx(nr+2),info

	tauR(0:nr)=tauR_in(0:nr)
	tauR(nr+1)=0d0
	
	Ma=0d0
	Mb=0d0
	Mc=0d0
	do ir=1,nr
		fact=1d0/(0.5d0*(tauR(ir+1)+tauR(ir))-0.5d0*(tauR(ir)+tauR(ir-1)))
		Mb(ir+1)=1d0+fact*(1d0/(tauR(ir+1)-tauR(ir))+1d0/(tauR(ir)-tauR(ir-1)))
		Ma(ir+1)=-fact*1d0/(tauR(ir)-tauR(ir-1))
		Mc(ir+1)=-fact*1d0/(tauR(ir+1)-tauR(ir))
	enddo
	Mb(1)=1d0/(tauR(0)-tauR(1))
	Mc(1)=-1d0/(tauR(0)-tauR(1))
	Ma(nr+2)=1d0/(tauR(nr))
	Mb(nr+2)=-1d0-1d0/(tauR(nr))
	x=0d0
	x(1)=1d0
	info=0
	call tridag(Ma,Mb,Mc,x,y,nr+2,info)
	Ij(0:nr)=y(1:nr+1)
	do ir=0,nr
		if(Ij(ir).lt.0d0) info=1
	enddo
	if(info.ne.0) then
		do ir=1,nr+2
			MM(ir,1)=Ma(ir)
			MM(ir,2)=Mb(ir)
			MM(ir,3)=Mc(ir)
		enddo
		call bandec(MM,nr+2,1,1,nr+2,3,MMal,1,indx,d)
		x=0d0
		x(1)=1d0
		call banbks(MM,nr+2,1,1,nr+2,3,MMal,1,indx,x)
		Ij(0:nr)=x(1:nr+1)
	endif
	do ir=0,nr
		if(Ij(ir).lt.0d0) then
			Ij(ir)=0d0
		endif
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

	W(1:N,1:N)=A(1:N,1:N)
	W(1:N,N+1)=x(1:N)

	W(N+1:2*N,1:N)=0d0
	W(N+1:2*N,N+1)=minx(1:N)
	do i=1,N
		W(i+N,i)=1d0
	enddo

	W(2*N+1:3*N,1:N)=0d0
	W(2*N+1:3*N,N+1)=-maxx(1:N)
	do i=1,N
		W(i+2*N,i)=-1d0
	enddo

	call dlsei(W, MDW, ME, MA, MG, N, PRGOPT, x, RNORME,
     +   RNORML, MODE, WS, IP)

	return
	end
	





	subroutine InvertIj(tauR_in,Linv,nr)
	IMPLICIT NONE
	integer ir,nr,iir,nrr
	real*8 tauR_in(nr),Linv(nr,nr),Lmat(0:nr+1,0:nr+1)
	real*8 tauR(0:nr+1),fact
	real*8 Ma(nr+2),Mb(nr+2),Mc(nr+2)
	integer info

	tauR(1:nr)=tauR_in(1:nr)
	tauR(nr+1)=0d0
	tauR(0)=tauR(1)+1d0

	Ma=0d0
	Mb=0d0
	Mc=0d0
	ir=1
	do ir=1,nr
		fact=1d0/(0.5d0*(tauR(ir+1)+tauR(ir))-0.5d0*(tauR(ir)+tauR(ir-1)))
		Mb(ir+1)=1d0+fact*(1d0/(tauR(ir+1)-tauR(ir))+1d0/(tauR(ir)-tauR(ir-1)))
		Ma(ir+1)=-fact*1d0/(tauR(ir)-tauR(ir-1))
		Mc(ir+1)=-fact*1d0/(tauR(ir+1)-tauR(ir))
	enddo
	Mb(1)=1d0+1d0/(tauR(0)-tauR(1))
	Mc(1)=-1d0/(tauR(0)-tauR(1))
	Ma(nr+2)=1d0/(tauR(nr))
	Mb(nr+2)=-1d0-1d0/(tauR(nr))
	Linv=0d0

	Lmat(0:nr+1,0:nr+1)=0d0
	do iir=0,nr+1
		Lmat(iir,iir)=1d0
	enddo
	nrr=nr+2
	call dgtsv(nrr,nrr,Ma(2:nr+2),Mb(1:nr+2),Mc(1:nr+1),Lmat(0:nr+1,0:nr+1),nrr,info)
	do ir=1,nr
		do iir=1,nr
			Linv(ir,iir)=Lmat(ir,iir)
		enddo
	enddo

	return
	end


	subroutine AddScatter(Si_in,tauR_in,Ca,Cs,Ce,nr,nu,wnu,nnu,NRHS)
	use Constants
	IMPLICIT NONE
	integer inu,nnu,ilam,ir,info,NRHS,nr,i,j
	real*8 tauR(nr),tauR_in(nr),Si_in(nr,NRHS)
	real*8 Si(nr,NRHS),Ca(nr),Cs(nr),Ce(nr)
	real*8 nu(nnu),wnu(nnu),albedo(1:nr)
	real*8 Ij(nr),Itot(nr,NRHS),Linv(nr,nr),Lmat(nr,nr)
	integer,allocatable :: IWORKomp(:)

	allocate(IWORKomp(10*nr*nr))

	do ir=1,nr
		albedo(ir)=min(1d0/(1d0+1d-4),Cs(ir)/Ce(ir))
	enddo

	Linv=0d0
	do inu=1,nnu
		tauR(1:nr)=tauR_in(1:nr)/abs(nu(inu))
		call InvertIj(tauR,Lmat,nr)
		do ir=1,nr
			Linv(ir,1:nr)=Linv(ir,1:nr)+wnu(inu)*Lmat(ir,1:nr)*albedo(1:nr)
		enddo
	enddo

	do ir=1,nr
		Linv(ir,ir)=1d0-Linv(ir,ir)
	enddo
	Itot=0d0
	Itot(1:nr,1:NRHS)=Si_in(1:nr,1:NRHS)
	info=0
	call DGESV( nr, NRHS, Linv, nr, IWORKomp, Itot, nr, info)

	Si(1:nr,1:NRHS)=Itot(1:nr,1:NRHS)
	do i=1,NRHS
	do ir=1,nr
		if(.not.Si(ir,i).gt.Si_in(ir,i)) then
			Si(ir,i)=Si_in(ir,i)
		endif
	enddo
	enddo
	Si_in=Si

	deallocate(IWORKomp)

	return
	end



	subroutine AddScatter_old(Si,Ca,Cs)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer inu,nnu,ir,info,NRHS
	parameter(nnu=5)
	real*8 tau,d,tauR_nu(nr),contr
	real*8 Si(nr),Si_in(nr),Ca(nr),Cs(nr),Ce(nr)
	real*8 nu(nnu),wnu(nnu),must
	real*8,allocatable :: tauR(:),Ij(:),Itot(:),Linv(:,:),Lmat(:,:),Iprev(:)
	integer,allocatable :: IWORKomp(:)

	Si_in=Si
	if(.not.scattering) return
	Ce=Ca+Cs
	tauR_nu=0d0
			do ir=nr,1,-1
				d=R(ir+1)-R(ir)
				tau=d*Ce(ir)
				if(P(ir).gt.Psimplecloud) then
					tau=tau+1d4
				endif
				if(.not.tau.gt.1d-10) then
					tau=1d-10
				endif
				if(tau.gt.1d10) then
					tau=1d10
				endif
				if(ir.eq.nr) then
					tauR_nu(ir)=tau
				else
					tauR_nu(ir)=tauR_nu(ir+1)+tau
				endif
			enddo

	call gauleg(0d0,1d0,nu,wnu,nnu)
	
		NRHS=1
		allocate(tauR(nr))
		allocate(Itot(nr))
		allocate(Linv(nr,nr))
		allocate(Lmat(nr,nr))
		allocate(IWORKomp(10*nr*nr))
					Linv=0d0
					do inu=1,nnu
						tauR(1:nr)=tauR_nu(1:nr)/abs(nu(inu))
						call InvertIj(tauR,Lmat,nr)
						do ir=1,nr
							Linv(ir,1:nr)=Linv(ir,1:nr)+wnu(inu)*Lmat(ir,1:nr)*Cs(1:nr)/(Ce(1:nr))
						enddo
					enddo
					Linv=-Linv
					do ir=1,nr
						Linv(ir,ir)=1d0+Linv(ir,ir)
					enddo
					Itot(1:nr)=Si_in(1:nr)
					call DGESV( nr, NRHS, Linv, nr, IWORKomp, Itot, nr, info )
					Si(1:nr)=Itot(1:nr)
					do ir=1,nr
						if(.not.Si(ir).gt.0d0) then
							Si(ir)=Si_in(ir)
						endif
					enddo
		deallocate(tauR)
		deallocate(Itot)
		deallocate(Linv)
		deallocate(Lmat)
		deallocate(IWORKomp)
	
	return
	end
	






