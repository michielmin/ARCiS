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

	if(doMCcompute.and.nTiter.gt.0) then
		call MCDoComputeT(converged,f)
		converged=.false.
		return
	else
		call DoComputeTfluxes(converged,f)
		return
	endif

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
				CrT(ir)=CrT(ir)+wgg(ig)*dfreq(ilam)*BB(iT,ilam)/Ca(ir,ilam,ig)
				tot2=tot2+wgg(ig)*dfreq(ilam)*BB(iT,ilam)

				CrV(ir)=CrV(ir)+wgg(ig)*dfreq(ilam)*Fstar(ilam)/Ca(ir,ilam,ig)
				tot3=tot3+wgg(ig)*dfreq(ilam)*Fstar(ilam)

				CrVe(ir)=CrVe(ir)+wgg(ig)*dfreq(ilam)*Fstar(ilam)/Ce(ir,ilam,ig)
				tot4=tot4+wgg(ig)*dfreq(ilam)*Fstar(ilam)
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
     &	3d0*Tirr**4*must*(2d0/3d0+1d0/(sqrt(3d0)*gamma)+(gamma/sqrt(3d0)-1d0/(sqrt(3d0)*gamma))*exp(-tau_T*gamma*sqrt(3d0)))/4d0
		enddo
		T0(ir)=T0(ir)/100d0
		T0(ir)=T0(ir)**0.25

		if(.not.T0(ir).gt.3d0) T0(ir)=3d0

c		if(ir.lt.nr) then
c			dlnP=log(P(ir+1)/P(ir))
c			dlnT=log(T0(ir+1)/T0(ir))
c			if((dlnT/dlnP).gt.nabla_ad(ir)) then
c				dlnT=(nabla_ad(ir))*dlnP
c				T0(ir)=T0(ir+1)/exp(dlnT)
c			endif
c		endif
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




























	subroutine DoComputeTfluxes(converged,f)
	use GlobalSetup
	use Constants
	use modComputeT
	use CloudModule
	IMPLICIT NONE
	integer iphase,iter
	real*8 tau_V,tau_T,Planck,Cp(nr),f
	real*8 g,dlnT,dlnP,d,tau,tautot,fact,contr,tau_a,exp_tau,nu
	real*8,allocatable :: Ce(:,:,:),Ca(:,:,:),Cs(:,:,:),taustar(:,:)
	real*8 tot,tot2,tot3,tot4,chi2,must,gamma,dP,Tirr,T0(nr),must_i,E,E0
	real*8 Cjstar(nr),Jedd(nr),Cj(nr),Ch(nr),z,Hstar(nr),Jtot,Htot,Ktot,AStar(nr)
	real*8 Jstar(nr),fedd(nr),Hedd(nr),lH1,lH2,P1,P2,scale,IntH(nr,nr),Fl2(nr)
	real*8 B1,B2,Si(0:nr+1),Ts(nr),tauR(nr+1),tottauR(nr+1),Fl(nr),minFl(nr),maxFl(nr)
	real*8,allocatable :: Jnu(:,:,:),Hnu(:,:,:),Knu(:,:,:),IntHnu(:,:,:),TempH(:)
	integer ir,ilam,ig,i,iT,niter,inu,nnu,jr,iTmin,iTmax,info,NRHS,IWORK(10*nr*nr)
	logical docloud0(max(nclouds,1)),converged
	type(Mueller) M	
	integer,allocatable :: IP(:)
	real*8,allocatable :: WS(:)

	allocate(IP(nr*4+2))
	IP(1)=(2*(nr)+nr*3+(nr*2+2)*(nr+7))
	IP(2)=nr*4+2
	allocate(WS(IP(1)))

	allocate(IntHnu(nlam,nr,nr))
	allocate(taustar(nlam,ng))
	allocate(Ce(nr,nlam,ng))
	allocate(Ca(nr,nlam,ng))
	allocate(Cs(nr,nlam,ng))
	allocate(Jnu(nr,nlam,ng))
	allocate(Hnu(nr,nlam,ng))
	allocate(Knu(nr,nlam,ng))

	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo

	T0(1:nr)=T(1:nr)

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
				call Crossections(ir,ilam,ig,Ca(ir,ilam,ig),Cs(ir,ilam,ig),docloud0)
				Ce(ir,ilam,ig)=Ca(ir,ilam,ig)+Cs(ir,ilam,ig)
			enddo
		enddo
	enddo
	Ce=Ca
	Cs=0d0

	nnu=5
	IntHnu=0d0
	Hstar=0d0
	Cjstar=0d0
	tauR=0d0
	tottauR=0d0
	tauR(nr+1)=0d0
	tottauR(nr+1)=1d0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ir,ilam,ig,inu,jr,tautot,TempH,fact,nu,d,tau_a,tau,exp_tau,contr)
!$OMP& SHARED(nr,nlam,ng,nnu,R,P,Ca,Cs,Psimplecloud,dfreq,wgg,Cjstar,Hstar,Fstar,
!$OMP&	Dplanet,must,IntHnu)
	allocate(TempH(nr))
!$OMP DO
!$OMP& SCHEDULE(STATIC)
	do ir=nr,1,-1
		do ilam=1,nlam-1
			TempH=0d0
			do ig=1,ng
				do inu=1,nnu
					tautot=0d0
					fact=1d0
					nu=-1d0*(real(inu)-0.5)/real(nnu)
					do jr=ir,nr
						d=abs(R(jr+1)-R(jr))/abs(nu)
						tau_a=d*Ca(jr,ilam,ig)
						tau=tau_a+d*Cs(jr,ilam,ig)
						if(P(jr).gt.Psimplecloud) tau=1d4
						exp_tau=exp(-tau)
						tautot=tautot+tau
						contr=((tau-1d0+exp_tau)/tau)*fact*tau_a/tau
						TempH(jr)=TempH(jr)+dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						if(jr.lt.nr) then
							contr=((1d0-(tau+1d0)*exp_tau)/tau)*fact*tau_a/tau
							TempH(jr+1)=TempH(jr+1)+dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						endif
						fact=fact*exp_tau
						if(tautot.gt.20d0) exit
					enddo
					tautot=0d0
					fact=1d0
					nu=1d0*(real(inu)-0.5)/real(nnu)
					do jr=ir,1,-1
						if(jr.gt.1) then
							d=abs(R(jr-1)-R(jr))/abs(nu)
							tau_a=d*Ca(jr-1,ilam,ig)
							tau=tau_a+d*Cs(jr-1,ilam,ig)
						else
							d=abs(R(jr+1)-R(jr))/abs(nu)
							tau_a=d*Ca(jr,ilam,ig)
							tau=tau_a+d*Cs(jr,ilam,ig)
						endif
						if(P(jr).gt.Psimplecloud) tau=1d4
						exp_tau=exp(-tau)
						tautot=tautot+tau
						contr=((tau-1d0+exp_tau)/tau)*fact*tau_a/tau
						TempH(jr)=TempH(jr)+dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						contr=((1d0-(tau+1d0)*exp_tau)/tau)*fact*tau_a/tau
						if(jr.gt.1) then
							TempH(jr-1)=TempH(jr-1)+dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						else
							TempH(jr)=TempH(jr)+dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						endif
						fact=fact*exp_tau
						if(tautot.gt.20d0) exit
					enddo
				enddo
				tautot=0d0
				fact=1d0
				do jr=ir,nr
					d=abs(R(jr+1)-R(jr))/abs(must)
					tau_a=d*Ca(jr,ilam,ig)
					tau=tau_a+d*Cs(jr,ilam,ig)
					if(P(jr).gt.Psimplecloud) tau=1d4
					tautot=tautot+tau
					exp_tau=exp(-tau)
					fact=fact*exp_tau
					if(tautot.gt.20d0) exit
				enddo
				contr=(Fstar(ilam)/(4d0*pi*Dplanet**2))*fact
				Hstar(ir)=Hstar(ir)-dfreq(ilam)*wgg(ig)*contr*must
				Cjstar(ir)=Cjstar(ir)+dfreq(ilam)*wgg(ig)*contr*Ca(ir,ilam,ig)
			enddo
			IntHnu(ilam,ir,1:nr)=TempH(1:nr)
		enddo
	enddo
!$OMP END DO
	deallocate(TempH)
!$OMP FLUSH
!$OMP END PARALLEL
	
	do iter=1,50
	
	Ts(1:nr)=T(1:nr)

	tauR=0d0
	tottauR=0d0
	tauR(nr+1)=0d0
	tottauR(nr+1)=1d0
	IntH=0d0
	Si(nr+1)=0d0
	do ilam=1,nlam-1
		do jr=nr,1,-1
			iT=T(jr)+1
			if(iT.gt.nBB-1) iT=nBB-1
			if(iT.lt.1) iT=1
			scale=(T(jr)/real(iT))**4
			Si(jr)=scale*BB(iT,ilam)
		enddo
		Si(0)=Si(1)
		do ir=nr,1,-1
			do jr=1,nr
				IntH(ir,jr)=IntH(ir,jr)+Si(jr)*IntHnu(ilam,ir,jr)
			enddo
			do ig=1,ng
				do jr=nr,1,-1
					if(jr.eq.ir) then
						tauR(ir)=tauR(ir)+dfreq(ilam)*wgg(ig)*abs(R(ir+1)-R(ir))*Ce(ir,ilam,ig)*Si(ir)/(must*0.5d0+0.5d0)
						tottauR(ir)=tottauR(ir)+dfreq(ilam)*wgg(ig)*Si(ir)
					endif	
				enddo
			enddo
		enddo
	enddo
	tauR=tauR/tottauR
	do ir=nr,1,-1
		tauR(ir)=tauR(ir+1)+tauR(ir)
	enddo

	E0=0d0
	iT=TeffP+1
	if(iT.gt.nBB-1) iT=nBB-1
	if(iT.lt.1) iT=1
	scale=(TeffP/real(iT))**4
	do ilam=1,nlam-1
		E0=E0+scale*dfreq(ilam)*BB(iT,ilam)
	enddo
	E0=E0/4d0

	do ir=1,nr
		iT=T(ir)+1
		if(iT.gt.nBB-1) iT=nBB-1
		if(iT.lt.1) iT=1
		E=0d0
		do ilam=1,nlam-1
			do ig=1,ng
				E=E+wgg(ig)*dfreq(ilam)*BB(iT,ilam)*Ca(ir,ilam,ig)
			enddo
		enddo
		minFl(ir)=Cjstar(ir)/E
		maxFl(ir)=max(minFl(ir)*10d0,100d0)
	enddo

	do ir=1,nr
		Fl(ir)=E0-Hstar(ir)
		Fl(ir)=Fl(ir)*P(ir)**(-0.25)
		IntH(ir,1:nr)=IntH(ir,1:nr)*P(ir)**(-0.25)
	enddo
	IP(1)=(2*(nr)+nr*3+(nr*2+2)*(nr+7))
	IP(2)=nr*4+2
	call PosSolve(IntH,Fl,minFl,maxFl,nr,IP,WS)

	Fl2=Fl
	do ir=1,nr
		Fl(ir)=0d0
		tot2=0d0
		do jr=1,nr
			tot=exp(-(real(ir-jr)/2d0)**2)
			tot2=tot2+tot
			Fl(ir)=Fl(ir)+Fl2(jr)*tot
		enddo
		Fl(ir)=Fl(ir)/tot2
	enddo

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ir,iT,E,iTmin,iTmax,E0,ilam,ig)
!$OMP& SHARED(nr,wgg,BB,Ca,nlam,ng,dfreq,Fl,Ts,T)
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC,1)
	do ir=1,nr
		iT=T(ir)+1
		if(iT.gt.nBB-1) iT=nBB-1
		if(iT.lt.1) iT=1
		E=0d0
		do ilam=1,nlam-1
			do ig=1,ng
				E=E+wgg(ig)*dfreq(ilam)*BB(iT,ilam)*Ca(ir,ilam,ig)
			enddo
		enddo
		E=E*min(max(1d-2,Fl(ir)),100d0)
		iTmin=1
		iTmax=nBB
		iT=(iTmax+iTmin)/2
		do while(abs(iTmax-iTmin).gt.1)
			E0=0d0
			do ilam=1,nlam-1
				do ig=1,ng
					E0=E0+wgg(ig)*dfreq(ilam)*BB(iT,ilam)*Ca(ir,ilam,ig)
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
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	do ir=nr,1,-1
		if(.not.Ts(ir).gt.3d0) then
			if(ir.eq.nr) then
				Ts(ir)=T(ir)
			else
				Ts(ir)=Ts(ir+1)
			endif
		endif
	enddo
	converged=.true.
	do ir=1,nr
		if(abs(T(ir)-Ts(ir))/(T(ir)+Ts(ir)).gt.epsiter) converged=.false.
	enddo

	do ir=1,nr
		T(ir)=sqrt(T(ir)*Ts(ir))
	enddo

	if(converged.and.iter.gt.5) exit
	enddo

	do ir=nr-1,1,-1
		dlnP=log(P(ir+1)/P(ir))
		dlnT=log(T(ir+1)/T(ir))
		if((dlnT/dlnP).gt.nabla_ad(ir)) then
			dlnT=(nabla_ad(ir))*dlnP
			T(ir)=T(ir+1)/exp(dlnT)
		endif
	enddo

	do ir=1,nr
		T(ir)=T0(ir)**(1d0-f)*T(ir)**f
	enddo

	converged=.true.
	do ir=1,nr
		if(abs(T(ir)-T0(ir))/(T(ir)+T0(ir)).gt.epsiter) converged=.false.
	enddo

	call WriteStructure

	deallocate(IntHnu)
	deallocate(taustar)
	deallocate(Ce)
	deallocate(Ca)
	deallocate(Cs)
	deallocate(Jnu)
	deallocate(Hnu)
	deallocate(Knu)

	deallocate(IP,WS)

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
	





















	subroutine DoComputeTeddington(converged,f)
	use GlobalSetup
	use Constants
	use modComputeT
	use CloudModule
	IMPLICIT NONE
	integer iphase,iter
	real*8 tau_V,tau_T,Planck,Cp(nr),f
	real*8 g,dlnT,dlnP,d,tau,tautot,fact,contr,tau_a,exp_tau,nu
	real*8,allocatable :: Ce(:,:,:),Ca(:,:,:),Cs(:,:,:),taustar(:,:)
	real*8 tot,tot2,tot3,tot4,chi2,must,gamma,dP,Tirr,T0(nr),must_i,E,E0
	real*8 Cjstar(nr),Jedd(nr),Cj(nr),Ch(nr),z,Hstar(nr),Jtot,Htot,Ktot
	real*8 Jstar(nr),fedd(nr),Hedd(nr),lH1,lH2,P1,P2,scale
	real*8,allocatable :: Jnu(:,:,:),Hnu(:,:,:),Knu(:,:,:)
	integer ir,ilam,ig,i,iT,niter,inu,nnu,jr,iTmin,iTmax
	logical docloud0(max(nclouds,1)),converged
	type(Mueller) M	
	real*8 tauf(nr),Si(0:nr+1),B1,B2,x1,x2,dx1,dx2,ax,bx
	integer info,IWORK(10*nr*nr),NRHS

	allocate(taustar(nlam,ng))
	allocate(Ce(nr,nlam,ng))
	allocate(Ca(nr,nlam,ng))
	allocate(Cs(nr,nlam,ng))
	allocate(Jnu(nr,nlam,ng))
	allocate(Hnu(nr,nlam,ng))
	allocate(Knu(nr,nlam,ng))

	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo
	
	do iter=1,50
	
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
				call Crossections(ir,ilam,ig,Ca(ir,ilam,ig),Cs(ir,ilam,ig),docloud0)
				Ca(ir,ilam,ig)=Ca(ir,ilam,ig)/dens(ir)
				Cs(ir,ilam,ig)=Cs(ir,ilam,ig)/dens(ir)
				Ce(ir,ilam,ig)=(Ca(ir,ilam,ig)+Cs(ir,ilam,ig)*(1d0-g))
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

	do ir=1,nr
		Hstar(ir)=0d0
		Jstar(ir)=0d0
		Cjstar(ir)=0d0
		do ilam=1,nlam-1
			do ig=1,ng
				tautot=0d0
				fact=1d0
				do jr=ir+1,nr
					d=abs(R(jr+1)-R(jr))/abs(must)
					tau_a=d*Ca(jr,ilam,ig)*dens(jr)
					tau=tau_a+d*Cs(jr,ilam,ig)*dens(jr)
					if(P(jr).gt.Psimplecloud) tau=1d4
					tautot=tautot+tau
					exp_tau=exp(-tau)
					fact=fact*exp_tau
					if(tautot.gt.20d0) exit
				enddo
				contr=(Fstar(ilam)/(pi*Dplanet**2))*fact
				Hstar(ir)=Hstar(ir)-must*dfreq(ilam)*wgg(ig)*contr
				Jstar(ir)=Jstar(ir)+dfreq(ilam)*wgg(ig)*contr
				Cjstar(ir)=Cjstar(ir)+dfreq(ilam)*wgg(ig)*contr*Ca(ir,ilam,ig)
			enddo
		enddo
	enddo

	E0=0d0
	iT=TeffP+1
	if(iT.gt.nBB-1) iT=nBB-1
	if(iT.lt.1) iT=1
	scale=(TeffP/real(iT))**4
	do ilam=1,nlam-1
		E0=E0+scale*dfreq(ilam)*BB(iT,ilam)
	enddo
	E0=E0*4d0*pi
	do ir=1,nr
		Hedd(ir)=E0-Hstar(ir)
	enddo

	Jnu=0d0
	Hnu=0d0
	Knu=0d0

	nnu=5
	do ilam=1,nlam-1
		iT=T(nr)+1
		if(iT.gt.nBB-1) iT=nBB-1
		if(iT.lt.1) iT=1
		scale=(T(nr)/real(iT))**4
		Si(nr+1)=0d0
		do jr=nr,1,-1
			iT=T(jr)+1
			if(iT.gt.nBB-1) iT=nBB-1
			if(iT.lt.1) iT=1
			scale=(T(jr)/real(iT))**4
			Si(jr)=scale*BB(iT,ilam)
		enddo
		Si(0)=Si(1)
		do ig=1,ng
			do ir=1,nr
				do inu=1,nnu
					tautot=0d0
					fact=1d0
					nu=-1d0*(real(inu)-0.5)/real(nnu)
					do jr=ir,nr
						B1=Si(jr)
						B2=Si(jr+1)
						d=abs(R(jr+1)-R(jr))/abs(nu)
						tau_a=d*Ca(jr,ilam,ig)
						tau=tau_a+d*Cs(jr,ilam,ig)
						if(P(jr).gt.Psimplecloud) tau=1d4
						exp_tau=exp(-tau)
						tautot=tautot+tau
						contr=(B1*(tau-1d0+exp_tau)/tau)*fact*tau_a/tau
						Jnu(ir,ilam,ig)=Jnu(ir,ilam,ig)+dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						Hnu(ir,ilam,ig)=Hnu(ir,ilam,ig)+nu*dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						Knu(ir,ilam,ig)=Knu(ir,ilam,ig)+nu*nu*dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						if(jr.lt.nr) then
							contr=(B2*(1d0-(tau+1d0)*exp_tau)/tau)*fact*tau_a/tau
							Jnu(ir,ilam,ig)=Jnu(ir,ilam,ig)+dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
							Hnu(ir,ilam,ig)=Hnu(ir,ilam,ig)+nu*dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
							Knu(ir,ilam,ig)=Knu(ir,ilam,ig)+nu*nu*dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						endif
						fact=fact*exp_tau
						if(tautot.gt.20d0) exit
					enddo
				enddo
				do inu=1,nnu
					tautot=0d0
					fact=1d0
					nu=1d0*(real(inu)-0.5)/real(nnu)
					do jr=ir,1,-1
						B1=Si(jr)
						B2=Si(jr-1)
						if(jr.gt.1) then
							d=abs(R(jr-1)-R(jr))/abs(nu)
							tau_a=d*Ca(jr-1,ilam,ig)
							tau=tau_a+d*Cs(jr-1,ilam,ig)
						else
							d=abs(R(jr+1)-R(jr))/abs(nu)
							tau_a=d*Ca(jr,ilam,ig)
							tau=tau_a+d*Cs(jr,ilam,ig)
						endif
						if(P(jr).gt.Psimplecloud) tau=1d4
						exp_tau=exp(-tau)
						tautot=tautot+tau
						contr=(B1*(tau-1d0+exp_tau)/tau)*fact*tau_a/tau
						Jnu(ir,ilam,ig)=Jnu(ir,ilam,ig)+dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						Hnu(ir,ilam,ig)=Hnu(ir,ilam,ig)+nu*dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						Knu(ir,ilam,ig)=Knu(ir,ilam,ig)+nu*nu*dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						contr=(B2*(1d0-(tau+1d0)*exp_tau)/tau)*fact*tau_a/tau
						Jnu(ir,ilam,ig)=Jnu(ir,ilam,ig)+dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						Hnu(ir,ilam,ig)=Hnu(ir,ilam,ig)+nu*dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						Knu(ir,ilam,ig)=Knu(ir,ilam,ig)+nu*nu*dfreq(ilam)*wgg(ig)*nu*contr/real(nnu*2)
						fact=fact*exp_tau
						if(tautot.gt.20d0) exit
					enddo
				enddo
			enddo
		enddo
	enddo

	do ir=1,nr
		Jtot=0d0
		Htot=0d0
		Ktot=0d0
		Cj(ir)=0d0
		Ch(ir)=0d0
		do ilam=1,nlam-1
			do ig=1,ng
				Jtot=Jtot+dfreq(ilam)*wgg(ig)*Jnu(ir,ilam,ig)
				Htot=Htot+dfreq(ilam)*wgg(ig)*Hnu(ir,ilam,ig)
				Ktot=Ktot+dfreq(ilam)*wgg(ig)*Knu(ir,ilam,ig)
				Cj(ir)=Cj(ir)+wgg(ig)*dfreq(ilam)*Jnu(ir,ilam,ig)*Ca(ir,ilam,ig)
				Ch(ir)=Ch(ir)+abs(wgg(ig)*dfreq(ilam)*Hnu(ir,ilam,ig)*Ce(ir,ilam,ig))
			enddo
		enddo
		fedd(ir)=Ktot/Jtot
		Cj(ir)=Cj(ir)/Jtot
		Ch(ir)=Ch(ir)/Htot
		print*,ir,fedd(ir)
	enddo

	Jedd(nr)=Hedd(nr)*Jtot/Htot
	do ir=nr-1,1,-1
		dx1=-dens(ir+1)*(Ch(ir+1)*Hedd(ir+1)/4d0)
		dx2=-dens(ir)*(Ch(ir)*Hedd(ir)/4d0)
		x1=R(ir+1)
		x2=R(ir)
		ax=0.5d0*(dx1-dx2)/(x1-x2)
		bx=dx1-2d0*ax*x1
		Jedd(ir)=fedd(ir+1)*Jedd(ir+1)-ax*x1**2-bx*x1+ax*x2**2+bx*x2
	print*,ir,fedd(ir+1)*Jedd(ir+1),(R(ir)-R(ir+1))*dens(ir+1)*(Ch(ir+1)*Hedd(ir+1)/4d0)
		Jedd(ir)=max(0d0,Jedd(ir)/fedd(ir))
	enddo
	
	do ir=nr,1,-1
		E=(Cjstar(ir)+Cj(ir)*Jedd(ir))
		iTmin=1
		iTmax=nBB
		iT=0.5d0*(iTmax+iTmin)
		do while(abs(iTmax-iTmin).gt.1)
			E0=0d0
			do ilam=1,nlam-1
				do ig=1,ng
					E0=E0+wgg(ig)*dfreq(ilam)*BB(iT,ilam)*Ca(ir,ilam,ig)
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
		T0(ir)=real(iT)*(E/E0)**0.25

		if(.not.T0(ir).gt.3d0) then
			print*,T0(ir),iT,ir,fedd(ir),Ch(ir)
			T0(ir)=3d0
		endif
c		if(ir.lt.nr) then
c			dlnP=log(P(ir+1)/P(ir))
c			dlnT=log(T0(ir+1)/T0(ir))
c			if((dlnT/dlnP).gt.nabla_ad(ir)) then
c				dlnT=(nabla_ad(ir))*dlnP
c				T0(ir)=T0(ir+1)/exp(dlnT)
c			endif
c		endif
	enddo

	converged=.true.
	do ir=1,nr
		if(abs(T(ir)-T0(ir))/(T(ir)+T0(ir)).gt.epsiter) converged=.false.
	enddo

	chi2=0d0
	do ir=1,nr
		chi2=chi2+((min(T(ir),2900d0)-min(T0(ir),2900d0))/((min(T0(ir),2900d0)+min(T(ir),2900d0))*epsiter))**2
	f=0.1
		T(ir)=T(ir)**(1d0-f)*T0(ir)**f
	enddo
	chi2=chi2/real(nr)

	call WriteStructure

	enddo
	converged=.false.

	deallocate(taustar)
	deallocate(Ce)
	deallocate(Ca)
	deallocate(Cs)
	deallocate(Jnu)
	deallocate(Hnu)
	deallocate(Knu)


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
			densv=(mu*mp/(kb*T0))*exp(BTP-ATP/T0)
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





