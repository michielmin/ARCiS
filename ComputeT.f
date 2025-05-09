	subroutine DoComputeT(converged,f)
	use GlobalSetup
	use Constants
	use CloudModule
	use TimingModule
	IMPLICIT NONE
	real*8 f
	logical converged
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
	real*8,allocatable,save :: Ce(:,:,:,:),Ca(:,:,:,:),Cs(:,:,:,:),tauR_nu(:,:,:,:)
	real*8,allocatable,save :: wscat(:,:,:,:),wabs(:,:,:,:),dtauR_nu(:,:,:,:)
	real*8 tot,tot2,tot3,tot4,chi2,must,gamma,dP,Tirr,T0(nr),T1(nr),must_i,E,E0,Tinp(nr)
	real*8 z,Hstar(nr),Jtot,Htot,Ktot
	real*8 fedd(nr),Hedd(nr),lH1,lH2,P1,P2
	real*8,allocatable,save :: Jnu(:,:,:),Hnu(:,:,:),Knu(:,:,:)
	integer ir,ilam,ig,i,iT,niter,inu,nnu,jr,iTmin,iTmax
	logical docloud0(max(nclouds,1)),converged,stopscat
	real*8 tauf(nr),Si(0:nr),B1,B2,x1,x2,dx1,dx2,ax,bx,ff,TT,maxErr,Fl0(nr),IntH0(nr,nr),error(nr)
	integer info,IWORK(10*(nr+1)*(nr+1)),NRHS,ii(3),iscat,nscat
	real*8 tau1,tau2,ee0,ee1,ee2,tauR(nr),Ij(nr),Ih(nr),scale,dtauR(nr),F11(180),Gsca
	integer nlam_LR,icc
	real*8 IntH(nr,nr),Fl(nr),Ts(nr),minFl(nr),maxFl(nr),maxfact,err,tot1,Pb(nr+1)
	real*8,allocatable,save :: lam_LR(:),dfreq_LR(:),freq_LR(:),BB_LR(:,:),IntHnu(:,:,:)
	integer i1,i2,ngF,j
	real*8 ww,w1,w2,FstarBottom,tauRoss,HUVstar,HUVstar_omp,SUMC
	real*8,allocatable,save :: Si_omp(:,:),Ih_omp(:),Ij_omp(:),tauR_omp(:),Hsurf(:,:)
	real*8,allocatable,save :: Fstar_LR(:),SurfEmis_LR(:),IntHnuSurf(:,:),tobesummed(:)
	integer,allocatable,save :: IP(:)
	real*8,allocatable,save :: WS(:),nu(:),wnu(:),Hstar_lam(:,:),Hsurf_lam(:)
	real*8,allocatable,save :: IhN(:,:,:),IjN(:,:,:),UVstar(:),UVstar_omp(:),HBottom(:)
!$OMP THREADPRIVATE(Si_omp,tauR_omp,Ih_omp,Ij_omp,IhN,IjN,Hsurf_lam,HBottom,UVstar_omp)
	logical,save :: first_entry=.true.
	character*500 file
	integer icloud
	logical Convec(0:nr),fixed(nr)
	real*8,allocatable,save :: Hstar0(:),deltaT(:,:),prevT(:,:),AT(:,:),al(:,:)
	real*8,save :: maxerr_prev=1d0

	if(deepredist.and.deepredisttype.eq.'fixflux') then
		if(.not.allocated(Hstar0)) allocate(Hstar0(nr))
	endif

	maxfact=4d0

	if(first_entry) then
		allocate(IP(nr*4+2))
		IP(1)=(2*(nr+1)+nr*3+(nr*2+2)*(nr*2+7))
		IP(2)=nr*4+2
		allocate(WS(IP(1)))
		allocate(deltaT(nr,maxiter),prevT(nr,maxiter),AT(nr,maxiter),al(0:maxiter,1))
	else
		IP(1)=(2*(nr)+nr*3+(nr*2+2)*(nr+7))
		IP(2)=nr*4+2
	endif
		
	nlam_LR=0
	do ilam=1,nlam
		if(RTgridpoint(ilam)) then
			nlam_LR=nlam_LR+1
		endif
	enddo
	nnu=5

	if(first_entry) then
		allocate(lam_LR(nlam_LR))
		allocate(freq_LR(nlam_LR))
		allocate(dfreq_LR(nlam_LR))
		allocate(BB_LR(nlam_LR,nBB))
		nlam_LR=0
		do ilam=1,nlam
			if(RTgridpoint(ilam)) then
				nlam_LR=nlam_LR+1
				lam_LR(nlam_LR)=lam(ilam)
				freq_LR(nlam_LR)=freq(ilam)
				dfreq_LR(nlam_LR)=dfreq(ilam)
			endif
		enddo
		BB_LR=0d0
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

		allocate(Ce(nr,nlam_LR,ng,ncc))
		allocate(Ca(nr,nlam_LR,ng,ncc))
		allocate(Cs(nr,nlam_LR,ng,ncc))
		allocate(Fstar_LR(nlam_LR))
		allocate(tauR_nu(0:nr,nlam_LR,ng,ncc))
		allocate(dtauR_nu(0:nr,nlam_LR,ng,ncc))
		allocate(IntHnu(nlam_LR,nr,nr))
		allocate(IntHnuSurf(nlam_LR,nr))
		allocate(SurfEmis_LR(nlam_LR))
		allocate(wscat(nr,nlam_LR,ng,ncc),wabs(nr,nlam_LR,ng,ncc))
		allocate(tobesummed(nlam_LR))
		allocate(UVstar(nr))
		allocate(Hstar_lam(nlam,nr),nu(nnu),wnu(nnu))
		call gauleg(0d0,1d0,nu,wnu,nnu)
!$OMP PARALLEL IF(useomp)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(nr,nnu)
		allocate(Si_omp(nr,0:nr+1),tauR_omp(nr),Ih_omp(nr),Ij_omp(nr))
		allocate(IhN(nr,0:nr+1,nnu),IjN(nr,0:nr+1,nnu))
		allocate(Hsurf_lam(nr),HBottom(nr))
		allocate(UVstar_omp(nr))
!$OMP FLUSH
!$OMP END PARALLEL
		first_entry=.false.
	endif


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
					do icc=1,ncc
						call Crossections(ir,ilam,ig,Ca(ir,i,ig,icc),Cs(ir,i,ig,icc),docloud(icc,1:nclouds),0,F11,Gsca,anisoscattstar)
						Ca(ir,i,ig,icc)=Ca(ir,i,ig,icc)/dens(ir)
						Cs(ir,i,ig,icc)=Cs(ir,i,ig,icc)*(1d0-Gsca)/dens(ir)
					enddo
				enddo
			enddo
		endif
	enddo

	Ce=Ca+Cs
	do ir=1,nr
		do ilam=1,nlam_LR
			do ig=1,ng
				do icc=1,ncc
					wabs(ir,ilam,ig,icc)=Ca(ir,ilam,ig,icc)/Ce(ir,ilam,ig,icc)
					if(.not.wabs(ir,ilam,ig,icc).gt.1d-4) then
						wabs(ir,ilam,ig,icc)=1d-4
						Ca(ir,ilam,ig,icc)=Cs(ir,ilam,ig,icc)/(1d0/wabs(ir,ilam,ig,icc)-1d0)
						Ce(ir,ilam,ig,icc)=Ca(ir,ilam,ig,icc)+Cs(ir,ilam,ig,icc)
					endif
					wscat(ir,ilam,ig,icc)=Cs(ir,ilam,ig,icc)/Ce(ir,ilam,ig,icc)
					if(.not.wscat(ir,ilam,ig,icc).gt.0d0) then
						wscat(ir,ilam,ig,icc)=0d0
						Cs(ir,ilam,ig,icc)=0d0
						Ce(ir,ilam,ig,icc)=Ca(ir,ilam,ig,icc)+Cs(ir,ilam,ig,icc)
					endif
				enddo
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
			do icc=1,ncc
				do jr=nr,1,-1
					if(jr.eq.nr) then
						d=(P(jr)-Pb(jr+1))*1d6/grav(jr)
						tau=d*Ce(jr,ilam,ig,icc)
					else
						d=(Pb(jr+1)-P(jr+1))*1d6/grav(jr+1)
						tau=d*Ce(jr+1,ilam,ig,icc)
						d=(P(jr)-Pb(jr+1))*1d6/grav(jr)
						tau=tau+d*Ce(jr,ilam,ig,icc)
					endif
					if(.not.tau.gt.0d0) then
						tau=0d0
					endif
					tau=tau+1d-10
					tau=1d0/(1d0/tau+1d-10)
					dtauR_nu(jr,ilam,ig,icc)=tau
				enddo
				dtauR_nu(nr,ilam,ig,icc)=dtauR_nu(nr-1,ilam,ig,icc)**2/dtauR_nu(nr-2,ilam,ig,icc)
				do jr=nr,1,-1
					if(jr.eq.nr) then
						tauR(jr)=dtauR_nu(nr,ilam,ig,icc)
					else
						tauR(jr)=tauR(jr+1)+dtauR_nu(jr,ilam,ig,icc)
					endif
				enddo
				tauR_nu(1:nr,ilam,ig,icc)=tauR(1:nr)
			enddo
		enddo
	enddo

c===============================================================
c quick thing to read in a file!
c	file='sedlhs.txt'
c	call regridlog(file,1d4*lam_LR,Fstar_LR,nlam_LR)
c	Fstar_LR=Fstar_LR*distance**2/1e23
c===============================================================
c	E=0d0
c	do ilam=1,nlam_LR
c		Fstar_LR(ilam)=Planck(Tstar,freq_LR(ilam))*pi*Rstar**2
c		E=E+dfreq_LR(ilam)*Fstar_LR(ilam)
c	enddo
c	scale=pi*Rstar**2*((2d0*(pi*kb*Tstar)**4)/(15d0*hplanck**3*clight**3))/E
c	Fstar_LR=Fstar_LR*scale
	
	ff=1d0
	
	Hstar(1:nr)=0d0
	do ir=1,nr
		do jr=1,nr
			IntHnu(:,ir,jr)=0d0
		enddo
		IntHnuSurf(:,ir)=0d0
	enddo
	UVstar=0d0
	HUVstar=0d0

!$OMP PARALLEL IF(useomp)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,ig,ir,inu,jr,contr,FstarBottom,tot,HUVstar_omp,icc)
!$OMP& SHARED(nlam_LR,ng,nr,nnu,tauR_nu,nu,wnu,dfreq_LR,wgg,IntHnu,SurfEmis_LR,dtauR_nu,Ca,Ce,Cs,Hsurf,night2day,deepredist,
!$OMP&			Hstar,Dplanet,Fstar_LR,must,wabs,wscat,IntHnuSurf,betaF,isoFstar,do3D,UVstar,HUVstar,
!$OMP&			lam_LR,deepredisttype,init3D,distrUV,Hstar_lam,ncc,cloudfrac)
	UVstar_omp=0d0
	HUVstar_omp=0d0
!$OMP DO
	do ilam=1,nlam_LR
		call tellertje(ilam,nlam_LR)
		Hstar_lam(ilam,1:nr)=0d0
		do icc=1,ncc
		do ig=1,ng
			Hsurf_lam(1:nr)=0d0
			contr=cloudfrac(icc)*(Fstar_LR(ilam)/(pi*Dplanet**2))
c Si_omp(0:nr,0) is the direct stellar contribution
			Si_omp(1:nr,0)=0d0
			FstarBottom=0d0
			if(isoFstar) then
				do inu=1,nnu
					tauR_omp(1:nr)=tauR_nu(1:nr,ilam,ig,icc)/abs(nu(inu))
					Ij_omp(1:nr)=contr*exp(-tauR_omp(1:nr))*betaF*2d0
					Ih_omp(1:nr)=-Ij_omp(1:nr)*nu(inu)
					Hstar_lam(ilam,1:nr)=Hstar_lam(ilam,1:nr)+wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
					Si_omp(1:nr,0)=Si_omp(1:nr,0)+wnu(inu)*Ij_omp(1:nr)*wscat(1:nr,ilam,ig,icc)/4d0
					FstarBottom=FstarBottom+wnu(inu)*abs(Ih_omp(1))
					if((.not.do3D.and..not.init3D).or.distrUV) then
						if(lam_LR(ilam).lt.0.4e-4) then
							UVstar_omp(1:nr)=UVstar_omp(1:nr)+2d0*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)
							HUVstar_omp=HUVstar_omp+wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(nr)
						endif
					endif
				enddo
				if((do3D.or.init3D).and..not.distrUV) then
					tauR_omp(1:nr)=tauR_nu(1:nr,ilam,ig,icc)/abs(max(must,1d-5))
					if(must.eq.0d0) then
						Ij_omp(1:nr)=0d0
					else
						Ij_omp(1:nr)=contr*exp(-tauR_omp(1:nr))
					endif
					Ih_omp(1:nr)=-betaF*Ij_omp(1:nr)
					if(lam_LR(ilam).lt.0.4e-4) then
						UVstar_omp(1:nr)=UVstar_omp(1:nr)+dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)
						HUVstar_omp=HUVstar_omp+dfreq_LR(ilam)*wgg(ig)*Ih_omp(nr)
					endif
				endif
			else
				tauR_omp(1:nr)=tauR_nu(1:nr,ilam,ig,icc)/abs(max(must,1d-5))
				Ij_omp(1:nr)=contr*exp(-tauR_omp(1:nr))
				Ih_omp(1:nr)=-betaF*Ij_omp(1:nr)
				Hstar_lam(ilam,1:nr)=Hstar_lam(ilam,1:nr)+dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
				FstarBottom=abs(Ih_omp(1))
				Si_omp(1:nr,0)=Si_omp(1:nr,0)+Ij_omp(1:nr)*wscat(1:nr,ilam,ig,icc)/4d0
				if(lam_LR(ilam).lt.0.4e-4) then
					UVstar_omp(1:nr)=UVstar_omp(1:nr)+dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)
					HUVstar_omp=HUVstar_omp+dfreq_LR(ilam)*wgg(ig)*Ih_omp(nr)
				endif
			endif
			if(deepredist.and.deepredisttype.eq.'deep'.and.betaF.gt.must) then
				Hstar_lam(ilam,1:nr)=Hstar_lam(ilam,1:nr)*must/betaF-(betaF-must)*dfreq_LR(ilam)*wgg(ig)*contr
			endif
			
c Si_omp(0:nr,1:nr) is the direct contribution from the atmosphere
			do ir=1,nr
				Si_omp(1:nr,ir)=0d0
				Si_omp(ir,ir)=wabs(ir,ilam,ig,icc)*cloudfrac(icc)
			enddo
c Si_omp(0:nr,nr+1) is the direct contribution from the surface
			Si_omp(1:nr,nr+1)=0d0
			do inu=1,nnu
				tauR_omp(1:nr)=tauR_nu(1:nr,ilam,ig,icc)/abs(nu(inu))
				Ij_omp(1:nr)=exp(-abs(tauR_omp(1:nr)-tauR_omp(1)))*cloudfrac(icc)
				Ih_omp(1:nr)=Ij_omp(1:nr)
				Si_omp(1:nr,nr+1)=Si_omp(1:nr,nr+1)+wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)/2d0
				Hsurf_lam(1:nr)=Hsurf_lam(1:nr)+2d0*FstarBottom*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*
     & 								Ih_omp(1:nr)*(1d0-SurfEmis_LR(ilam))
				IntHnuSurf(ilam,1:nr)=IntHnuSurf(ilam,1:nr)+2d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
			enddo
			Si_omp(1:nr,nr+1)=Si_omp(1:nr,nr+1)*wscat(1:nr,ilam,ig,icc)

			call AddScatter_lin(Si_omp(1:nr,0:nr+1),tauR_nu(1:nr,ilam,ig,icc),
     &					Ca(1:nr,ilam,ig,icc),Cs(1:nr,ilam,ig,icc),Ce(1:nr,ilam,ig,icc),(1d0-SurfEmis_LR(ilam)),nr,nu,wnu,nnu,nr+2)

			do inu=1,nnu
				tauR_omp(1:nr)=dtauR_nu(1:nr,ilam,ig,icc)/abs(nu(inu))
				call SolveIjhExpN(tauR_omp(1:nr),Si_omp(1:nr,0:nr+1),IhN(1:nr,0:nr+1,inu),IjN(1:nr,0:nr+1,inu),nr,nr+2)
			enddo

			do inu=1,nnu
				Ih_omp(1:nr)=IhN(1:nr,0,inu)
				Ij_omp(1:nr)=IjN(1:nr,0,inu)
				Hstar_lam(ilam,1:nr)=Hstar_lam(ilam,1:nr)+4d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
				if(lam_LR(ilam).lt.0.4e-4) then
					UVstar_omp(1:nr)=UVstar_omp(1:nr)+4d0*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ij_omp(1:nr)
					HUVstar_omp=HUVstar_omp+2d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(nr)
				endif
			enddo
			HBottom=0d0
			do ir=1,nr
				do inu=1,nnu
					Ih_omp(1:nr)=IhN(1:nr,ir,inu)
					Ij_omp(1:nr)=IjN(1:nr,ir,inu)
					IntHnu(ilam,1:nr,ir)=IntHnu(ilam,1:nr,ir)+4d0*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih_omp(1:nr)
					HBottom(ir)=HBottom(ir)+4d0*nu(inu)*wnu(inu)*Ij_omp(1)
				enddo
			enddo

			do inu=1,nnu
				tauR_omp(1:nr)=tauR_nu(1:nr,ilam,ig,icc)/abs(nu(inu))
				Ij_omp(1:nr)=exp(-abs(tauR_omp(1:nr)-tauR_omp(1)))
				Ih_omp(1:nr)=Ij_omp(1:nr)
				do ir=1,nr
					IntHnu(ilam,1:nr,ir)=IntHnu(ilam,1:nr,ir)+2d0*HBottom(ir)*nu(inu)*wnu(inu)*dfreq_LR(ilam)*wgg(ig)*
     &						Ih_omp(1:nr)*(1d0-SurfEmis_LR(ilam))
				enddo
			enddo
			do inu=1,nnu
				Ih_omp(1:nr)=IhN(1:nr,nr+1,inu)
				Ij_omp(1:nr)=IjN(1:nr,nr+1,inu)
				Hsurf_lam(1:nr)=Hsurf_lam(1:nr)+FstarBottom*4d0*nu(inu)*wnu(inu)*Ih_omp(1:nr)*(1d0-SurfEmis_LR(ilam))
				IntHnuSurf(ilam,1:nr)=IntHnuSurf(ilam,1:nr)+4d0*nu(inu)*wnu(inu)*Ih_omp(1:nr)
			enddo
			if(lam_LR(ilam).lt.0.4e-4) then
				HUVstar_omp=HUVstar_omp+dfreq_LR(ilam)*wgg(ig)*FstarBottom*SurfEmis_LR(ilam)
				HUVstar_omp=HUVstar_omp+Hsurf_lam(nr)
			endif

			do ir=1,nr
				Hstar_lam(ilam,ir)=Hstar_lam(ilam,ir)+Hsurf_lam(ir)
			enddo
		enddo
		enddo
	enddo
!$OMP END DO
!$OMP CRITICAL
	UVstar=UVstar+UVstar_omp
	HUVstar=HUVstar+HUVstar_omp
!$OMP END CRITICAL
!$OMP FLUSH
!$OMP END PARALLEL
	do ir=1,nr
		Hstar(ir)=SUMC(Hstar_lam(1:nlam_LR,ir),nlam_LR)
	enddo
	scaleUV=max(0d0,-(pi*HUVstar/3.4947466112306125E-009))
	UVstar(1:nr)=UVstar(1:nr)/UVstar(nr)
	if(nTiter.eq.1) then
		tauUV=UVstar
	else
		tauUV=(tauUV+UVstar)/2d0
	endif
	if(do3D) then
		tot=0d0
		do ilam=1,nlam_LR
			tot=tot+dfreq_LR(ilam)*(Fstar_LR(ilam)/(pi*Dplanet**2))
		enddo
		local_albedo(i3D)=1d0+Hstar(nr)/(tot*betaF)
	endif
	do ir=nr,1,-1
		tot=0d0
		kappaUV(ir)=0d0
		do ilam=1,nlam_LR
			if(lam_LR(ilam).lt.0.4e-4) then
				do ig=1,ng
				do icc=1,ncc
					kappaUV(ir)=kappaUV(ir)+Ca(ir,ilam,ig,icc)*Fstar_LR(ilam)*dfreq_LR(ilam)*wgg(ig)*cloudfrac(icc)
					tot=tot+Fstar_LR(ilam)*dfreq_LR(ilam)*wgg(ig)*cloudfrac(icc)
				enddo
				enddo
			endif
		enddo
		kappaUV(ir)=kappaUV(ir)/tot
	enddo
c	tot=0d0
c	do ilam=1,nlam_LR
c		if(lam_LR(ilam).lt.0.4e-4) then
c			tot=tot+dfreq_LR(ilam)*Fstar_LR(ilam)/(Dplanet**2)
c		endif
c	enddo
c	print*,'EUV: ',tot

	if(init3D.and.deepredist.and.deepredisttype.eq.'fixflux') then
		Hstar0(1:nr)=Hstar(1:nr)/betaF
	endif
	if(do3D.and.deepredist.and.deepredisttype.eq.'fixflux') then
		if(betaF.gt.0d0) then
			Hstar(1:nr)=Hstar(1:nr)*betaT/betaF
		endif
		do ir=1,nr
			Hstar(ir)=min(Hstar(ir)+Hstar0(ir)*(betaF-betaT),0d0)
		enddo
	endif

	iter=1
	iter2=1
	Tsurface=T(1)

c=========================================================================================
c=== start of loop =======================================================================
c=========================================================================================
	do while(iter.le.niter.and.iter2.le.niter*5)
	iter=iter+1

	Ts(1:nr)=T(1:nr)

	Hedd=0d0
	E0=2d0*(((pi*kb*TeffP)**4)/(15d0*hplanck**3*clight**3))
	do ir=1,nr
		Hedd(ir)=Hedd(ir)+E0-Hstar(ir)
	enddo

	do ir=1,nr
		IntH(:,ir)=0d0
	enddo
	do jr=1,nr
		iT=T(jr)+1
		if(iT.gt.nBB-1) iT=nBB-1
		if(iT.lt.1) iT=1
		scale=(T(jr)/real(iT))**4
		do ir=1,nr
			do ilam=1,nlam_LR
				tobesummed(ilam)=scale*BB_LR(ilam,iT)*IntHnu(ilam,ir,jr)
			enddo
			IntH(ir,jr)=SUMC(tobesummed,nlam_LR)
		enddo
	enddo
	iT=Tsurface+1
	if(iT.gt.nBB-1) iT=nBB-1
	if(iT.lt.1) iT=1
	scale=(Tsurface/real(iT))**4
	do ir=1,nr
		do ilam=1,nlam_LR
			tobesummed(ilam)=scale*BB_LR(ilam,iT)*IntHnuSurf(ilam,ir)*SurfEmis_LR(ilam)
		enddo
		IntH(ir,1)=IntH(ir,1)+SUMC(tobesummed,nlam_LR)
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
	do ir=1,nr
		IntH0(:,ir)=IntH(:,ir)
	enddo
	Convec(0:nr)=.false.
	fixed=.false.

1	continue
	iter2=iter2+1

	call DGESV( nr, NRHS, IntH, nr, IWORK, Fl, nr, info )
c	call PosSolve(IntH,Fl,minFl,maxFl,nr,IP,WS)

	do ir=1,nr
		Fl(ir)=min(max(1d-5,Fl(ir)),1d5)
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

	do ir=1,nr
		Fl(ir)=min(max(0.5d0,Fl(ir)),2d0)
	enddo
				
	Ts=T*Fl**0.25

c	if(i3D.ne.n3D.and..not.init3D) then
c		do ir=1,nr
c			if(Ts(ir).gt.Tprev3D(ir).and..not.fixed(ir)) then
c				IntH0(ir,1:nr)=0d0
c				IntH0(ir,ir)=1d0
c				Fl0(ir)=(Tprev3D(ir)/T(ir))**4
c				Fl=Fl0
c				IntH=IntH0
c				fixed(ir)=.true.
c				goto 1
c			endif
c		enddo
c	endif

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
c=========================================================================================
c=== end of loop =========================================================================
c=========================================================================================

c=========================================================================================
c=== computation of convective flux ======================================================
c=========================================================================================
	if(convectKzz) then
	Hedd=0d0
	E0=2d0*(((pi*kb*TeffP)**4)/(15d0*hplanck**3*clight**3))
	do ir=1,nr
		Hedd(ir)=Hedd(ir)+E0-Hstar(ir)
		IntH(:,ir)=0d0
	enddo

	do jr=1,nr
		iT=T(jr)+1
		if(iT.gt.nBB-1) iT=nBB-1
		if(iT.lt.1) iT=1
		scale=(T(jr)/real(iT))**4
		do ir=1,nr
			do ilam=1,nlam_LR
				tobesummed(ilam)=scale*BB_LR(ilam,iT)*IntHnu(ilam,ir,jr)
			enddo
			IntH(ir,jr)=SUMC(tobesummed,nlam_LR)
		enddo
	enddo
	iT=Tsurface+1
	if(iT.gt.nBB-1) iT=nBB-1
	if(iT.lt.1) iT=1
	scale=(Tsurface/real(iT))**4
	do ir=1,nr
		do ilam=1,nlam_LR
			tobesummed(ilam)=scale*BB_LR(ilam,iT)*IntHnuSurf(ilam,ir)*SurfEmis_LR(ilam)
		enddo
		IntH(ir,1)=IntH(ir,1)+SUMC(tobesummed,nlam_LR)
	enddo

	Fl=Hedd
	do ir=1,nr
		if(Convec(ir)) then
			Fl(ir)=Fl(ir)-SUMC(IntH(ir,1:nr),nr)
			Fl(ir)=Fl(ir)*sigma/(2d0*(((pi*kb)**4)/(15d0*hplanck**3*clight**3)))
			Kzz_convect(ir)=(1d0/3d0)*Hp(ir)*((exp_ad-1d0)*abs(Fl(ir))/(exp_ad*MMW(ir)*dens(ir)))**(1d0/3d0)
		else
			Kzz_convect(ir)=0d0
		endif
	enddo
	T1=Kzz_convect
	do ir=2,nr-1
		Kzz_convect(ir)=(T1(ir-1)*T1(ir+1)*T1(ir)**2)**0.25
	enddo
	
	endif
c=========================================================================================

	T1=T
	do ir=2,nr-1
		T(ir)=(T1(ir-1)*T1(ir+1)*T1(ir)**2)**0.25
	enddo

	if(fixnight2day) then
		tot=0d0
		tauLW=0d0
		iT=Tsurface+1
		if(iT.gt.nBB-1) iT=nBB-1
		if(iT.lt.1) iT=1
		do ilam=1,nlam_LR
			do ig=1,ng
			do icc=1,ncc
				tauLW=tauLW+BB_LR(ilam,iT)*dfreq_LR(ilam)*wgg(ig)*cloudfrac(icc)/tauR_nu(1,ilam,ig,icc)
				tot=tot+BB_LR(ilam,iT)*dfreq_LR(ilam)*wgg(ig)*cloudfrac(icc)
			enddo
			enddo
		enddo
		tauLW=tot/tauLW
	endif
	

	maxErr=0d0
	do ir=1,nr-1
		error(ir)=abs(T(ir)-Tinp(ir))/(T(ir)+Tinp(ir))
		if(error(ir).gt.maxErr) maxErr=error(ir)
	enddo
	call sort(error,nr)
	converged=.false.
	j=9*nr/10
	if(maxErr.lt.epsiter*3d0.and.error(j).lt.epsiter) converged=.true.

	if(.not.allocated(Tdist)) allocate(Tdist(nr,maxiter))
	Tdist(1:nr,nTiter)=T(1:nr)
	call output("Maximum error on T-struct: " // dbl2string(maxErr*100d0,'(f6.2)') // "%")
	if(converged.and.maxErr.gt.epsiter) call output("      90% of values below: " // dbl2string(error(j)*100d0,'(f6.2)') // "%")
	if(do3D.and..not.retrieval.and..not.dopostequalweights) then
		print*,"Maximum error on T-struct: " // dbl2string(maxErr*100d0,'(f6.2)') // "%"
		if(converged.and.maxErr.gt.epsiter) print*,"      90% of values below: " // dbl2string(error(j)*100d0,'(f6.2)') // "%"
	endif

	deltaT(1:nr,nTiter)=T(1:nr)-Tinp(1:nr)
	prevT(1:nr,nTiter)=T(1:nr)
	
	T0(1:nr)=Tinp(1:nr)
	T1(1:nr)=T(1:nr)
	if(nTiter.gt.2) then
		j=max(2,min(min(nTiter/2-1,nr-2),8))
		call FindNext(deltaT,prevT,T,nr,nTiter,j,IP,WS)
		T=f*T1+(1d0-f)*T
		tot=max(2d0*epsiter,maxErr)
		if(nTiter.gt.10) then
			tot=tot*(1d0/real(nTiter-10))
		endif
		do ir=1,nr
			if(T(ir).gt.Tinp(ir)*(1d0+tot)**1.5) T(ir)=Tinp(ir)*(1d0+tot)**1.5
			if(T(ir).lt.Tinp(ir)*(1d0-tot)**1.5) T(ir)=Tinp(ir)*(1d0-tot)**1.5
		enddo
	else
		do ir=1,nr
			T(ir)=f*T1(ir)+(1d0-f)*T0(ir)
		enddo
	endif
	f=max(f/1.5d0,0.05d0)

	maxerr_prev=maxErr

	Tsurface=T(1)
	call output("Surface temperature: " // dbl2string(Tsurface,'(f8.2)') // " K")
	if(do3D.and..not.retrieval.and..not.dopostequalweights) print*,"Surface temperature: " // dbl2string(Tsurface,'(f8.2)') // " K"

	if(writefiles) then
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

	return
	end


	subroutine MidpointDeriv(dx,y,dy,n)
	IMPLICIT NONE
	integer n,i
	real*8 dx(n),y(n),dy(n),dx1,dx2

	dy(1)=(y(2)-y(1))/dx(1)
	do i=2,n-1
		dx1=dx(i)
		dx2=dx(i+1)
		dy(i)=(y(i)-y(i-1))/dx1+dx1*((y(i+1)-y(i))/dx2-(y(i)-y(i-1))/dx1)/(dx1+dx2)
	enddo
	dy(n)=(y(n)-y(n-1))/dx(n)

	return
	end
	

	subroutine ComputeDeriv(x,y,dy,n,yp1,ypn)
	IMPLICIT NONE
	integer n,i
	real*8 x(n),y(n),dy(n),d2y(n),yp1,ypn,dx1,dx2
	real*8 A,B

	dy(1)=(y(2)-y(1))/(x(2)-x(1))
	dy(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
	do i=2,n-1
		dx1=(x(i)-x(i-1))
		dx2=(x(i+1)-x(i))
		dy(i)=(y(i)-y(i-1))/dx1+dx1*((y(i+1)-y(i))/dx2-(y(i)-y(i-1))/dx1)/(dx1+dx2)
	enddo
	return

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
	real*8 tau1,tau2,S1,S2,I12,dtau
	
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
	return

	x1=dtau(nr)
	do ir=nr-1,2,-1
		x1=x1+dtau(ir)
		if(x1.gt.1d-6) then
			s0=Ij(ir,1:nrhs)
			s1=Ij(ir+1,1:nrhs)
			Ih(ir,1:nrhs)=(s0-s1)/dtau(ir)
		endif
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
		endif
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
	


	subroutine FindNext(A,P,x,M,N,Nmax,IP,WS)
	IMPLICIT NONE
	integer N,NN,i,M,Nmax
	real*8 A(M,N),x(M),W(M+Nmax+1,Nmax+1),y(M),P(M,N),ymin

	integer ME,MA,MG,MODE,MDW
	real*8 PRGOPT(10),RNORME,RNORML
	integer IP(*)
	real*8 WS(*)

	ymin=0d0

	ME=1
	MA=M
	MG=Nmax
	MDW=M+Nmax+1

	PRGOPT(1)=1!4
	PRGOPT(2)=1
	PRGOPT(3)=1
	PRGOPT(4)=7
	PRGOPT(5)=10
	PRGOPT(6)=1
	PRGOPT(7)=1

	W(1,1:Nmax+1)=1d0
	do i=1,M
		W(1+i,1:Nmax)=A(i,N-Nmax+1:N)/P(i,N)
		W(1+i,Nmax+1)=0d0
	enddo

	W(M+2:M+1+Nmax,1:Nmax+1)=0d0
	do i=1,Nmax
		W(M+1+i,i)=1d0
	enddo
	W(M+2:M+1+Nmax,Nmax+1)=ymin

	call dlsei(W, MDW, ME, MA, MG, Nmax, PRGOPT, y, RNORME,
     +   RNORML, MODE, WS, IP)

	do i=1,Nmax
		if(y(i).lt.ymin) y(i)=ymin
	enddo
	y=y/sum(y(1:Nmax))
c	print*,y(1:Nmax)
	x=0d0
	do i=1,Nmax
		x(1:M)=x(1:M)+y(i)*P(1:M,N-Nmax+i)
	enddo

	return
	end
	



	subroutine AddScatter_lin(Si_in,tauR_in,Ca,Cs,Ce,SurfAlb,nr,nu,wnu,nnu,NRHS)
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
		Si(1:nr,1:NRHS)=matmul(Linv,Itot(1:nr,1:NRHS))
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



	subroutine AddScatter(Si_in,tauR_in,Ca,Cs,Ce,SurfAlb,nr,nu,wnu,nnu,NRHS)
	use Constants
	IMPLICIT NONE
	integer nextra
	parameter(nextra=50)
	integer inu,nnu,ilam,ir,info,NRHS,nr,i,j,jr,nre
	real*8 tauR(nr+nextra),tauR_in(nr),Si_in(nr,NRHS),tauR_e(nr+nextra)
	real*8 Ca(nr),Cs(nr),Ce(nr),Si_e(nr+nextra,NRHS)
	real*8 nu(nnu),wnu(nnu),albedo(nr+nextra),SurfAlb,temp(nr)
	real*8 Linv(nr+nextra,nr+nextra),Lmat(nr+nextra,nr+nextra),Hsurf(nr+nextra)
	integer IWORKomp(nr+nextra),ii(nr+nextra)
	logical doit

	Si_e=0d0
	ii=0
	do ir=1,nr
		albedo(ir)=Cs(ir)/Ce(ir)
		if(.not.albedo(ir).lt.1d0/(1d0+1d-4)) albedo(ir)=1d0/(1d0+1d-4)
		Si_e(ir,1:NRHS)=Si_in(ir,1:NRHS)
		ii(ir)=ir
	enddo

	tauR_e(1:nr)=tauR_in(1:nr)

	doit=.true.
	nre=nr
	tauR(1:nre)=tauR_e(1:nre)
	do while(doit)
		doit=.false.
		do ir=nre,2,-1
			if(tauR_e(ir-1).gt.1d-2.and.tauR_e(ir-1).lt.5d0) then
				if(tauR_e(ir-1)-tauR_e(ir).gt.tauR_e(ir)*1.25.and.nre.lt.nr+nextra) then
					ii(ir+1:nre+1)=ii(ir:nre)
					ii(ir)=0
					tauR(ir+1:nre+1)=tauR(ir:nre)
					tauR(ir)=(tauR(ir-1)+tauR(ir+1))/2d0
					albedo(ir+1:nre+1)=albedo(ir:nre)
					albedo(ir)=(albedo(ir-1)+albedo(ir+1))/2d0
					do i=1,NRHS
						Si_e(ir+1:nre+1,i)=Si_e(ir:nre,i)
						if(Si_e(ir+1,i).le.0d0.or.Si_e(ir-1,i).le.0d0) then
							Si_e(ir,i)=(Si_e(ir-1,i)+Si_e(ir+1,i))/2d0
						else
							Si_e(ir,i)=sqrt(Si_e(ir-1,i)*Si_e(ir+1,i))
						endif
					enddo
					doit=.true.
					nre=nre+1
				endif
			endif
		enddo
		tauR_e(1:nre)=tauR(1:nre)
	enddo

	Linv=0d0
	Hsurf=0d0
	do inu=1,nnu
		tauR(1:nre)=tauR_e(1:nre)/abs(nu(inu))
		call InvertIjExp(tauR,Lmat(1:nre,1:nre),nre)
		do ir=1,nre
			Linv(ir,1:nre)=Linv(ir,1:nre)+wnu(inu)*Lmat(ir,1:nre)*albedo(ir)
		enddo
		Hsurf(1:nre)=Hsurf(1:nre)+2d0*nu(inu)*wnu(inu)*Lmat(1,1:nre)
	enddo
	do inu=1,nnu
		tauR(1:nre)=abs((tauR_e(1:nre)-tauR_e(1))/nu(inu))
		do ir=1,nre
			Linv(ir,1:nre)=Linv(ir,1:nre)+wnu(inu)*SurfAlb*Hsurf(1:nre)*albedo(ir)*exp(-tauR(ir))
		enddo
	enddo

	Linv=-Linv
	do ir=1,nre
		Linv(ir,ir)=1d0+Linv(ir,ir)
	enddo
	info=0
	call DGESV( nre, NRHS, Linv, nr+nextra, IWORKomp, Si_e, nr+nextra, info)
	if(info.ne.0) then
		print*,"problem in scattering",info
		return
	endif

	do i=1,NRHS
		do ir=1,nre
			if(ii(ir).gt.0) then
				if(Si_e(ir,i).ge.Si_in(ii(ir),i)) then
					Si_in(ii(ir),i)=Si_e(ir,i)
				endif
			endif
		enddo
	enddo

	return
	end


