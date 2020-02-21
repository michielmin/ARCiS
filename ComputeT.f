	subroutine DoComputeT(converged,f)
	use GlobalSetup
	use Constants
	use modComputeT
	use CloudModule
	IMPLICIT NONE
	integer iphase,Titer
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
	use modComputeT
	IMPLICIT NONE
	integer iphase,iter
	real*8 tau_V,tau_T,Planck,Cp(nr),f
	real*8 g,dlnT,dlnP,d,tau,tautot,fact,contr,tau_a,exp_tau
	real*8,allocatable :: Ce(:,:,:),Ca(:,:,:),Cs(:,:,:),taustar(:,:),tauR_nu(:,:,:)
	real*8 tot,tot2,tot3,tot4,chi2,must,gamma,dP,Tirr,T0(nr),must_i,E,E0,Tinp(nr)
	real*8 Cjstar(nr),Jedd(nr),Cj(nr),Ch(nr),z,Hstar(nr),Jtot,Htot,Ktot
	real*8 fedd(nr),Hedd(nr),lH1,lH2,P1,P2
	real*8,allocatable :: Jnu(:,:,:),Hnu(:,:,:),Knu(:,:,:)
	integer ir,ilam,ig,i,iT,niter,inu,nnu,jr,iTmin,iTmax
	logical docloud0(max(nclouds,1)),converged,stopscat
	type(Mueller) M	
	real*8 tauf(nr),Si(1:nr),B1,B2,x1,x2,dx1,dx2,ax,bx,ff,TT
	integer info,IWORK(10*nr*nr),NRHS,ii(3),iscat,nscat
	real*8 tau1,tau2,ee0,ee1,ee2,tauR(nr),Ij(nr),Ih(nr),Itot(nr),scale
	integer nlam_LR
	real*8 specres_LR,IntH(nr,nr),Fl(nr),Ts(nr),Si_p(nr),minFl(nr),maxFl(nr)
	real*8,allocatable :: lam_LR(:),dfreq_LR(:),freq_LR(:),BB_LR(:,:),IntHnu(:,:,:)
	integer i1,i2,ngF,j
	real*8 ww,w1,w2
	real*8,allocatable :: temp_a(:),wtemp(:),Ca_HR(:,:),Cs_HR(:,:),Fstar_LR(:)
	integer,allocatable :: IP(:)
	real*8,allocatable :: WS(:),nu(:),wnu(:)
	
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
		if(scale.gt.2d0.and.j.lt.nBB) then
			do iT=1,j
				BB_LR(iT,1:nlam_LR)=BB_LR(j+1,1:nlam_LR)*(real(iT)/real(j+1))**4
			enddo
			exit
		endif
		BB_LR(j,i)=BB_LR(j,i)*scale
	enddo


	allocate(taustar(nlam_LR,ng))
	allocate(Ca_HR(nlam,ng))
	allocate(Cs_HR(nlam,ng))
	allocate(Ce(nr,nlam_LR,ng))
	allocate(Ca(nr,nlam_LR,ng))
	allocate(Cs(nr,nlam_LR,ng))
	allocate(Fstar_LR(nlam_LR))
	allocate(tauR_nu(nr,nlam_LR,ng))
	allocate(IntHnu(nlam_LR,nr,nr))

	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo
	
	niter=50
	nscat=10
	epsiter=1d-3

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
			call GetMatrix(ir,ilam,M,docloud0)
			g=0d0
			tot=0d0
			do iphase=1,180
				g=g+M%F11(iphase)*costheta(iphase)*sintheta(iphase)
				tot=tot+M%F11(iphase)*sintheta(iphase)
			enddo
			g=g/tot
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
					if(tot2.ne.0d0) then
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

	do ilam=1,nlam_LR-1
		do ig=1,ng
			do ir=nr,1,-1
				d=abs(P(ir+1)-P(ir))*1d6/grav(ir)
				tau=d*Ce(ir,ilam,ig)
				if(P(ir).gt.Psimplecloud) then
					tau=tau+1d4
					Ce(ir,ilam,ig)=tau/d
					Ca(ir,ilam,ig)=Ce(ir,ilam,ig)-Cs(ir,ilam,ig)
				endif
				if(tau.lt.1d-5) then
					tau=1d-5
					Ce(ir,ilam,ig)=tau/d
					Ca(ir,ilam,ig)=Ce(ir,ilam,ig)-Cs(ir,ilam,ig)
				endif
				if(tau.gt.1d10) then
					tau=1d10
					Ce(ir,ilam,ig)=tau/d
					Ca(ir,ilam,ig)=Ce(ir,ilam,ig)
					Cs(ir,ilam,ig)=0d0
				endif
				if(ir.lt.nr) then
					tauR(ir)=tauR(ir+1)+tau
				else
					tauR(ir)=tau
				endif
			enddo
			tauR_nu(1:nr,ilam,ig)=tauR(1:nr)
		enddo
	enddo

	nnu=5
	allocate(nu(nnu),wnu(nnu))
	call gauleg(0d0,1d0,nu,wnu,nnu)


	Hstar=0d0
	Cjstar=0d0
	do ilam=1,nlam_LR-1
		do ig=1,ng
			do ir=1,nr
				fact=exp(-tauR_nu(ir,ilam,ig)/abs(must))
				contr=(Fstar_LR(ilam)/(4d0*pi*Dplanet**2))*fact
				Hstar(ir)=Hstar(ir)-must*dfreq_LR(ilam)*wgg(ig)*contr
				Cjstar(ir)=Cjstar(ir)+dfreq_LR(ilam)*wgg(ig)*contr*Ca(ir,ilam,ig)
			enddo
		enddo
	enddo


	E0=((2d0*(pi*kb*TeffP)**4)/(15d0*hplanck**3*clight**3))
	do ir=1,nr
		Hedd(ir)=E0/2d0-Hstar(ir)
	enddo

	ff=0.5d0
	
	T0(1:nr)=T(1:nr)

	IntHnu=0d0
	do ilam=1,nlam_LR-1
		call tellertje(ilam,nlam_LR-1)
		do ig=1,ng
			do ir=1,nr
				tautot=0d0
				do jr=ir,nr
					tautot=tautot+tauR_nu(jr,ilam,ig)/abs(must)
				enddo
				fact=exp(-tautot)
				contr=must*(Fstar_LR(ilam)/(4d0*pi*Dplanet**2))*fact
				iT=T(ir)+1
				if(iT.gt.nBB-1) iT=nBB-1
				if(iT.lt.1) iT=1
				scale=(T(ir)/real(iT))**4
				Si(ir)=1d0!scale*BB_LR(iT,ilam)*Ca(ir,ilam,ig)/Ce(ir,ilam,ig)+contr*Cs(ir,ilam,ig)/(Ce(ir,ilam,ig)*4d0*pi)
			enddo
			do ir=1,nr
				Si_p=0d0
				Si_p(ir)=1d0
c				tauR(1:nr)=tauR_nu(1:nr,ilam,ig)
c				call AddScatter(Si_p,tauR,Ca(1:nr,ilam,ig),Cs(1:nr,ilam,ig),Ce(1:nr,ilam,ig),nr,nu,wnu,nnu)
				do inu=1,nnu
					tauR(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(nu(inu))
					call SolveIj(tauR,Si_p,Ij,nr)
c					call spldiff(tauR,Ij,nr,Ih)
					do jr=1,nr-1
						Ih(jr)=(Ij(jr)-Ij(jr+1))/(tauR(jr)-tauR(jr+1))
					enddo
					Ih(nr)=Ij(nr)
					IntHnu(ilam,1:nr,ir)=IntHnu(ilam,1:nr,ir)+wnu(inu)*dfreq_LR(ilam)*wgg(ig)*Ih(1:nr)
				enddo
			enddo
		enddo
	enddo
	

	do iter=1,50

	Ts(1:nr)=T(1:nr)

	IntH=0d0
	do ir=1,nr
		do ilam=1,nlam_LR-1
			do jr=nr,1,-1
				iT=T(jr)+1
				if(iT.gt.nBB-1) iT=nBB-1
				if(iT.lt.1) iT=1
				scale=(T(jr)/real(iT))**4
				IntH(ir,jr)=IntH(ir,jr)+scale*BB_LR(iT,ilam)*IntHnu(ilam,ir,jr)
			enddo
		enddo
		Fl(ir)=Hedd(ir)
	enddo
	minFl=0.1d0
	maxFl=10d0

	IP(1)=(2*(nr)+nr*3+(nr*2+2)*(nr+7))
	IP(2)=nr*4+2
	NRHS=1

c	call DGESV( nr, NRHS, IntH, nr, IWORK, Fl, nr, info )
	call PosSolve(IntH,Fl,minFl,maxFl,nr,IP,WS)

	do ir=1,nr
		print*,Fl(ir)
		Fl(ir)=min(max(0.1d0,Fl(ir)),10d0)
	enddo

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
	enddo

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
		T(ir)=0.5d0*(T(ir)+Ts(ir))
	enddo

	if(converged.and.iter.gt.5) exit
	enddo

	do ir=nr-1,1,-1
		if(ir.lt.nr) then
			dlnP=log(P(ir+1)/P(ir))
			dlnT=log(T(ir+1)/T(ir))
			if((dlnT/dlnP).gt.nabla_ad(ir)) then
				dlnT=(nabla_ad(ir))*dlnP
				T(ir)=T(ir+1)/exp(dlnT)
			endif
		endif
	enddo

	converged=.true.
	do ir=1,nr
		if(abs(T(ir)-Tinp(ir))/(T(ir)+Tinp(ir)).gt.epsiter) converged=.false.
		T(ir)=Tinp(ir)*(1d0-f)+T(ir)*f
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

	subroutine SolveIj(tauR_in,Si,Ij,nr)
	IMPLICIT NONE
	integer ir,nr
	real*8 tauR_in(nr),Ij(nr),Si(nr),x(nr+2),y(nr+2),fact,d,tauR(0:nr+1)
	real*8 MM(nr+2,3),MMal(nr+2,1),Ma(nr+2),Mb(nr+2),Mc(nr+2)
	integer indx(nr+2),info

	tauR(1:nr)=tauR_in(1:nr)
	tauR(0)=tauR(1)*2d0
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
	Mb(1)=-1d0
	Mc(1)=1d0
	Ma(nr+2)=1d0/(tauR(nr))
	Mb(nr+2)=-1d0-1d0/(tauR(nr))
	x=0d0
	x(2:nr+1)=Si(1:nr)
	info=0
	call tridag(Ma,Mb,Mc,x,y,nr+2,info)
	Ij(1:nr)=y(2:nr+1)
	do ir=1,nr
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
		x(2:nr+1)=Si(1:nr)
		call banbks(MM,nr+2,1,1,nr+2,3,MMal,1,indx,x)
		Ij(1:nr)=x(2:nr+1)
	endif
	do ir=1,nr
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
	





	subroutine InvertIj(tauR,Linv,nr)
	IMPLICIT NONE
	integer ir,nr,iir
	real*8 tauR(nr),Ij(nr),Si(nr),Linv(nr,nr),Lmat(nr,nr)
	real*8 x(nr),y(nr),fact,d
	real*8 MM(nr,3),MMal(nr,1),Ma(nr),Mb(nr),Mc(nr)
	integer indx(nr),info

	Ma=0d0
	Mb=0d0
	Mc=0d0
	do ir=2,nr-1
		fact=1d0/(0.5d0*(tauR(ir+1)+tauR(ir))-0.5d0*(tauR(ir)+tauR(ir-1)))
		Mb(ir)=1d0+fact*(1d0/(tauR(ir+1)-tauR(ir))+1d0/(tauR(ir)-tauR(ir-1)))
		Ma(ir)=-fact*1d0/(tauR(ir)-tauR(ir-1))
		Mc(ir)=-fact*1d0/(tauR(ir+1)-tauR(ir))
	enddo
	Mb(1)=1d0/(tauR(1)-tauR(2))
	Mc(1)=-1d0/(tauR(1)-tauR(2))
	Ma(nr)=1d0/(tauR(nr-1)-tauR(nr))
	Mb(nr)=-1d0-1d0/(tauR(nr-1)-tauR(nr))
	Linv=0d0

	Lmat=0d0
	do iir=1,nr
		Lmat(iir,iir)=1d0
	enddo
	call dgtsv(nr,nr,Ma(2:nr),Mb(1:nr),Mc(1:nr-1),Lmat,nr,info)
	do ir=1,nr
		do iir=1,nr
			Linv(ir,iir)=Lmat(ir,iir)
		enddo
	enddo

	return

	do iir=2,nr-1
		x=0d0
		x(iir)=1d0
		info=0
		call tridag(Ma,Mb,Mc,x,y,nr,info)
		Ij(1:nr)=y(1:nr)
		do ir=1,nr
			if(Ij(ir).lt.0d0) info=1
		enddo
		if(info.ne.0) then
			do ir=1,nr
				MM(ir,1)=Ma(ir)
				MM(ir,2)=Mb(ir)
				MM(ir,3)=Mc(ir)
			enddo
			call bandec(MM,nr,1,1,nr,3,MMal,1,indx,d)
			x=0d0
			x(iir)=1d0
			call banbks(MM,nr,1,1,nr,3,MMal,1,indx,x)
			Ij(1:nr)=x(1:nr)
		endif
		do ir=1,nr
			if(Ij(ir).lt.0d0) then
				Ij(ir)=0d0
			endif
		enddo
		Linv(iir,1:nr)=Ij(1:nr)
	enddo


	return
	end


	subroutine SolveIjStar(tauR,I0,Ij,nr)
	IMPLICIT NONE
	integer ir,nr
	real*8 tauR(nr),Ij(nr),I0
	real*8 x(nr),y(nr),fact,d
	real*8 MM(nr,3),MMal(nr,1),Ma(nr),Mb(nr),Mc(nr)
	integer indx(nr),info

	Ma=0d0
	Mb=0d0
	Mc=0d0
	do ir=2,nr-1
		fact=1d0/(0.5d0*(tauR(ir+1)+tauR(ir))-0.5d0*(tauR(ir)+tauR(ir-1)))
		Mb(ir)=1d0+fact*(1d0/(tauR(ir+1)-tauR(ir))+1d0/(tauR(ir)-tauR(ir-1)))
		Ma(ir)=-fact*1d0/(tauR(ir)-tauR(ir-1))
		Mc(ir)=-fact*1d0/(tauR(ir+1)-tauR(ir))
	enddo
	Mb(1)=1d0/(tauR(1)-tauR(2))
	Mc(1)=-1d0/(tauR(1)-tauR(2))
	Ma(nr)=0d0
	Mb(nr)=1d0
	x=0d0
	x(nr)=I0
	info=0
	call tridag(Ma,Mb,Mc,x,y,nr,info)
	Ij(1:nr)=y(1:nr)
	do ir=1,nr
		if(Ij(ir).lt.0d0) info=1
	enddo
	if(info.ne.0) then
		do ir=1,nr
			MM(ir,1)=Ma(ir)
			MM(ir,2)=Mb(ir)
			MM(ir,3)=Mc(ir)
		enddo
		call bandec(MM,nr,1,1,nr,3,MMal,1,indx,d)
		x=0d0
		x(nr)=I0
		call banbks(MM,nr,1,1,nr,3,MMal,1,indx,x)
		Ij(1:nr)=x(1:nr)
	endif
	do ir=1,nr
		if(Ij(ir).lt.0d0) then
			Ij(ir)=0d0
		endif
	enddo

	return
	end
	


	subroutine AddScatter(Si,tauR,Ca,Cs,Ce,nr,nu,wnu,nnu)
	use Constants
	IMPLICIT NONE
	integer inu,nnu,ilam,ir,info,NRHS,nr
	real*8 tauR(nr)
	real*8 Si(nr),Ca(nr),Cs(nr),Ce(nr)
	real*8 nu(nnu),wnu(nnu)
	real*8 Ij(nr),Itot(nr),Linv(nr,nr),Lmat(nr,nr)
	integer,allocatable :: IWORKomp(:)

	NRHS=1
	allocate(IWORKomp(10*nr*nr))
	Linv=0d0
	do inu=1,nnu
		call InvertIj(tauR/nu(inu),Lmat,nr)
		do ir=1,nr
			Linv(ir,1:nr)=Linv(ir,1:nr)+wnu(inu)*Lmat(ir,1:nr)*Cs(1:nr)/(Ce(1:nr))
		enddo
	enddo
	Linv=-Linv
	do ir=1,nr
		Linv(ir,ir)=1d0+Linv(ir,ir)
	enddo
	Itot(1:nr)=Si(1:nr)
	call DGESV( nr, NRHS, Linv, nr, IWORKomp, Itot, nr, info )
	Si(1:nr)=Itot(1:nr)
	deallocate(IWORKomp)

	return
	end










