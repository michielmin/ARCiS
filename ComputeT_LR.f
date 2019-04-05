	subroutine DoComputeT(converged,f,Titer)
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
	IMPLICIT NONE
	integer iphase,iter
	real*8 tau_V,tau_T,Planck,Cp(nr),f
	real*8 g,dlnT,dlnP,d,tau,tautot,fact,contr,tau_a,exp_tau,nu
	real*8,allocatable :: Ce(:,:,:),Ca(:,:,:),Cs(:,:,:),taustar(:,:),tauR_nu(:,:,:)
	real*8 tot,tot2,tot3,tot4,chi2,must,gamma,dP,Tirr,T0(nr),must_i,E,E0,Tinp(nr)
	real*8 Cjstar(nr),Jedd(nr),Cj(nr),Ch(nr),z,Hstar(nr),Jtot,Htot,Ktot
	real*8 fedd(nr),Hedd(nr),lH1,lH2,P1,P2
	real*8,allocatable :: Jnu(:,:,:),Hnu(:,:,:),Knu(:,:,:),Hstar_nu(:,:,:),Jstar_nu(:,:,:)
	integer ir,ilam,ig,i,iT,niter,inu,nnu,jr,iTmin,iTmax
	logical docloud0(max(nclouds,1)),converged,stopscat
	type(Mueller) M	
	real*8 tauf(nr),Si(1:nr),B1,B2,x1,x2,dx1,dx2,ax,bx,ff
	integer info,IWORK(10*nr*nr),NRHS,ii(3),iscat,nscat
	real*8 tau1,tau2,ee0,ee1,ee2,tauR(0:nr+1),Ij(nr),Ih(nr),Itot(nr)
	integer nlam_LR
	real*8 specres_LR
	real*8,allocatable :: lam_LR(:),dfreq_LR(:),freq_LR(:)
	integer i1,i2,ngF,j
	real*8 ww,w1,w2
	real*8,allocatable :: temp_a(:),wtemp(:),Ca_HR(:,:),Cs_HR(:,:),Fstar_LR(:)
	
	specres_LR=10d0
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

	allocate(taustar(nlam_LR,ng))
	allocate(Ca_HR(nlam,ng))
	allocate(Cs_HR(nlam,ng))
	allocate(Ce(nr,nlam_LR,ng))
	allocate(Ca(nr,nlam_LR,ng))
	allocate(Cs(nr,nlam_LR,ng))
	allocate(Fstar_LR(nlam_LR))
	allocate(tauR_nu(0:nr+1,nlam_LR,ng))
	allocate(Jnu(nr,nlam_LR,ng))
	allocate(Hnu(nr,nlam_LR,ng))
	allocate(Knu(nr,nlam_LR,ng))

	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
	enddo

	Tinp(1:nr)=T(1:nr)
	
	niter=50
	nscat=5
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

	open(unit=20,file='HRopac.dat',RECL=6000)

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
				Fstar_LR(ilam)=0d0
			endif
		enddo
		if(ir.eq.nr/2) then
			do ilam=1,nlam-1
				do ig=1,ng
					write(20,*) lam(ilam),Ca_HR(ilam,ig)
				enddo
			enddo
		endif
	enddo
	deallocate(temp_a)
	deallocate(wtemp)
	Ce=Ca+Cs

	close(unit=20)

	open(unit=20,file='LRopac.dat',RECL=6000)
	do ilam=1,nlam_LR-1
		do ig=1,ng
			write(20,*) lam_LR(ilam),Ca(1:nr,ilam,ig)
		enddo
	enddo
	close(unit=20)

	do ir=1,nr
		Hstar(ir)=0d0
		Cjstar(ir)=0d0
		do ilam=1,nlam_LR-1
			do ig=1,ng
				tautot=0d0
				fact=1d0
				do jr=ir,nr
					d=abs(R(jr+1)-R(jr))*dens(jr)/abs(must)
					tau_a=d*Ca(jr,ilam,ig)
					tau=tau_a+d*Cs(jr,ilam,ig)
					if(P(jr).gt.Psimplecloud) tau=1d4
					tautot=tautot+tau
					exp_tau=exp(-tau)
					fact=fact*exp_tau
					if(tautot.gt.20d0) exit
				enddo
				contr=(Fstar_LR(ilam)/(4d0*pi*Dplanet**2))*fact/2d0
				Hstar(ir)=Hstar(ir)-must*dfreq_LR(ilam)*wgg(ig)*contr
				Cjstar(ir)=Cjstar(ir)+dfreq_LR(ilam)*wgg(ig)*contr*Ca(ir,ilam,ig)
			enddo
		enddo
	enddo

	E0=0d0
	do ilam=1,nlam_LR-1
		E0=E0+dfreq_LR(ilam)*Planck(TeffP,freq_LR(ilam))
	enddo
	E0=E0/(2d0*pi)
	do ir=1,nr
		Hedd(ir)=E0-Hstar(ir)
	enddo

	do ilam=1,nlam_LR-1
		do ig=1,ng
			tauR(nr+1)=0d0
			do ir=nr,1,-1
				d=abs(R(ir+1)-R(ir))*dens(ir)
				tau=d*Ce(ir,ilam,ig)
				if(P(ir).gt.Psimplecloud) tau=1d4
				if(tau.lt.1d-6) tau=1d-6
				tauR(ir)=tauR(ir+1)+tau
			enddo
			tauR(0)=2d0*tauR(1)
			tauR_nu(0:nr+1,ilam,ig)=tauR(0:nr+1)
		enddo
	enddo


	ff=0.5d0
	do iter=1,niter
	call tellertje(iter,niter+1)
	
	T0(1:nr)=T(1:nr)

	Jnu=0d0
	Hnu=0d0
	Knu=0d0

	nnu=5
	do ilam=1,nlam_LR-1
		do ig=1,ng
			do ir=1,nr
				fact=1d0
				tautot=0d0
				do jr=ir,nr
					tautot=tautot+tauR_nu(jr,ilam,ig)/abs(must)
				enddo
				fact=exp(-tautot)
				contr=(Fstar_LR(ilam)/(4d0*pi*Dplanet**2))*fact/2d0
				Si(ir)=Planck(T(ir),freq_LR(ilam))*Ca(ir,ilam,ig)/Ce(ir,ilam,ig)+contr*Cs(ir,ilam,ig)/(Ce(ir,ilam,ig)*4d0*pi)
			enddo
			do iscat=1,nscat
				Itot=0d0
				do inu=1,nnu
					nu=(real(inu)-0.5)/real(nnu)
					tauR(0:nr+1)=tauR_nu(0:nr+1,ilam,ig)/abs(nu)
					call SolveIj(tauR,Si,Ij,nr)
					Itot=Itot+Ij/real(nnu)
				enddo
				Itot(1:nr)=Si(1:nr)+Itot(1:nr)*Cs(1:nr,ilam,ig)/(Ce(1:nr,ilam,ig)*4d0*pi)
				stopscat=.true.
				do ir=1,nr
					if(abs((Itot(ir)-Si(ir))/(Itot(ir)+Si(ir))).gt.1d-1) stopscat=.false.
				enddo
				Si=Itot
				if(stopscat) exit
			enddo
			do inu=1,nnu
				nu=(real(inu)-0.5)/real(nnu)
				tauR(0:nr+1)=tauR_nu(0:nr+1,ilam,ig)/abs(nu)
				call SolveIj(tauR,Si,Ij,nr)
				do ir=1,nr-1
					Ih(ir)=-(Ij(ir)-Ij(ir+1))/(tauR(ir)-tauR(ir+1))
				enddo
				Ih(nr)=-(Ij(nr-1)-Ij(nr))/(tauR(nr-1)-tauR(nr))
				do ir=1,nr
					Jnu(ir,ilam,ig)=Jnu(ir,ilam,ig)+Ij(ir)/real(nnu)
					Hnu(ir,ilam,ig)=Hnu(ir,ilam,ig)-nu*Ih(ir)/real(nnu)
					Knu(ir,ilam,ig)=Knu(ir,ilam,ig)+nu*nu*Ij(ir)/real(nnu)
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
		do ilam=1,nlam_LR-1
			do ig=1,ng
				Jtot=Jtot+dfreq_LR(ilam)*wgg(ig)*Jnu(ir,ilam,ig)
				Htot=Htot+dfreq_LR(ilam)*wgg(ig)*Hnu(ir,ilam,ig)
				Ktot=Ktot+dfreq_LR(ilam)*wgg(ig)*Knu(ir,ilam,ig)
				Cj(ir)=Cj(ir)+wgg(ig)*dfreq_LR(ilam)*Jnu(ir,ilam,ig)*Ca(ir,ilam,ig)
				Ch(ir)=Ch(ir)+wgg(ig)*dfreq_LR(ilam)*Hnu(ir,ilam,ig)*Ce(ir,ilam,ig)
			enddo
		enddo
		fedd(ir)=Ktot/Jtot
		Cj(ir)=Cj(ir)/Jtot
		Ch(ir)=abs(Ch(ir)/Htot)
	enddo

	ir=nr
	Jtot=0d0
	Htot=0d0
	Ktot=0d0
	do ilam=1,nlam_LR-1
		do ig=1,ng
			Jtot=Jtot+dfreq_LR(ilam)*wgg(ig)*Jnu(ir,ilam,ig)
			Htot=Htot+dfreq_LR(ilam)*wgg(ig)*Hnu(ir,ilam,ig)
			Ktot=Ktot+dfreq_LR(ilam)*wgg(ig)*Knu(ir,ilam,ig)
		enddo
	enddo

	Jedd(nr)=max(Jtot,Hedd(nr)*Jtot/Htot)
	do ir=nr-1,1,-1
		dx1=-Ch(ir)*Hedd(ir)*dens(ir)
		x1=R(ir+1)
		x2=R(ir)
		Jedd(ir)=fedd(ir+1)*Jedd(ir+1)+dx1*(x2-x1)
		Jedd(ir)=max(0d0,Jedd(ir)/fedd(ir))
	enddo
	
	do ir=nr,1,-1
		E=Cjstar(ir)+Cj(ir)*Jedd(ir)
		E0=0d0
		do ilam=1,nlam_LR-1
			do ig=1,ng
				E0=E0+wgg(ig)*dfreq_LR(ilam)*Planck(T(ir),freq_LR(ilam))*Ca(ir,ilam,ig)
			enddo
		enddo

		T0(ir)=T(ir)*(E/E0)**0.25d0

		iTmin=1
		iTmax=nBB
		iT=T0(ir)+1
		if(iT.gt.nBB-1) iT=nBB-1
		if(iT.lt.1) iT=1
		do while(abs(iTmax-iTmin).gt.1)
			E0=0d0
			do ilam=1,nlam_LR-1
				do ig=1,ng
					E0=E0+wgg(ig)*dfreq_LR(ilam)*Planck(T0(ir),freq_LR(ilam))*Ca(ir,ilam,ig)
				enddo
			enddo
			if(E0.gt.E) then
				iTmax=iT
			else
				iTmin=iT
			endif
			if(abs(iTmax-iTmin).gt.1) then
				iT=0.5d0*(iTmax+iTmin)
				if(iT.lt.1) iT=1
				if(iT.gt.nBB) iT=nBB
				T0(ir)=real(iT)
			endif
		enddo

		if(.not.T0(ir).gt.3d0) then
			print*,T0(ir),iT,ir,fedd(ir),Ch(ir)
			T0(ir)=T(ir)
		endif
	enddo

	converged=.true.
	do ir=1,nr
		if(abs(T(ir)-T0(ir))/(T(ir)+T0(ir)).gt.epsiter) converged=.false.
	enddo

	do ir=1,nr
		T(ir)=(1d0-ff)*T(ir)+ff*T0(ir)
	enddo

	if(converged) exit
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

	deallocate(taustar)
	deallocate(Ce)
	deallocate(Ca)
	deallocate(Cs)
	deallocate(Hnu)
	deallocate(Knu)
	deallocate(tauR_nu)


	return
	end


	subroutine SolveIj(tauR,Si,Ij,nr)
	IMPLICIT NONE
	integer ir,nr
	real*8 tauR(0:nr+1),Ij(nr),Si(nr),x(nr+2),y(nr+2),fact,d
	real*8 MM(nr+2,3),MMal(nr+2,1),Ma(nr+2),Mb(nr+2),Mc(nr+2)
	integer indx(nr+2),info

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
      if(b(1).eq.0.) then
      	info=1
      	return
      endif
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
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


      SUBROUTINE banmul(a,n,m1,m2,np,mp,x,b)
      INTEGER m1,m2,mp,n,np
      REAL*8 a(np,mp),b(n),x(n)
      INTEGER i,j,k
      do 12 i=1,n
        b(i)=0.
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
      PARAMETER (TINY=1.e-20)
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
          a(i,j)=0.
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
          a(i,mm)=0.
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








