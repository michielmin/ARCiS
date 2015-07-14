	subroutine DoComputeT(converged)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nphase,iphase
	real*8 z,dz,E,tau,Planck,random,v,dx,dy
	real*8,allocatable :: Ca(:,:,:),Cs(:,:),Ce(:,:,:),g(:,:),spec(:,:)
	real*8 vR1,vR2,b,rr,R1,R2,tau_v,x,y,theta,ct1,ct2,Tc(nr),Ec(nr),EJv(nr)
	real*8 tot,specsurface(nlam,1),Lplanet,tot2
	real*8 CabsL(nlam),Ca0,Cs0,T0,E0,Eabs,chi2,CabsLG(nlam,ng),Crw(nr),nabla,rho
	integer iphot,ir,Nphot,ilam,ig,nscat,jrnext,NphotStar,NphotPlanet,jr,ir0,jr0
	integer iT1,iT2,iT
	logical docloud(nclouds),goingup,onedge,dorw(nr),converged
	type(Mueller),allocatable :: M(:,:)
	
	allocate(spec(nlam,ng))
	allocate(Ca(nr,nlam,ng))
	allocate(Ce(nr,nlam,ng))
	allocate(Cs(nr,nlam))
	allocate(g(nr,nlam))
	allocate(M(nr,nlam))
	
	call output("Temperature computation (in beta phase!!)")

	NphotPlanet=100
	NphotStar=1000

	call InitRandomWalk()

	Lplanet=0d0
	do ilam=1,nlam-1
		specsurface(ilam,1)=Planck(TeffP,freq(ilam))*4d0*pi*R(1)**2
		Lplanet=Lplanet+ABS(freq(ilam)-freq(ilam+1))*specsurface(ilam,1)
	enddo
	ir0=nr
	do ilam=1,nlam
		do ig=1,ng
			tau=0d0
			do ir=1,nr
				call Crossections(ir,ilam,ig,Ca(ir,ilam,ig),Cs(ir,ilam),docloud)
				Ce(ir,ilam,ig)=Ca(ir,ilam,ig)+Cs(ir,ilam)
			enddo
		enddo
	enddo
	do ilam=1,nlam
		do ig=1,ng
			tau=0d0
			do ir=nr,2,-1
				tau=tau+Ce(ir,ilam,ig)*(R(ir+1)-R(ir))
				if(tau.gt.maxtau) exit
			enddo
			if(ir.lt.ir0) ir0=ir
		enddo
	enddo
	ir0=1

	do ir=1,nr
		iT=T(ir)+1
		if(iT.gt.nBB-1) iT=nBB-1
		Crw(ir)=0d0
		tot2=0d0
		do ilam=1,nlam-1
			call GetMatrix(ir,ilam,M(ir,ilam),docloud)
			g(ir,ilam)=0d0
			tot=0d0
			do iphase=1,180
				g(ir,ilam)=g(ir,ilam)+M(ir,ilam)%F11(iphase)*costheta(iphase)*sintheta(iphase)
				tot=tot+M(ir,ilam)%F11(iphase)*sintheta(iphase)
			enddo
			g(ir,ilam)=g(ir,ilam)/tot
			do ig=1,ng
				Crw(ir)=Crw(ir)+ABS(freq(ilam)-freq(ilam+1))*(BB(iT+1,ilam)-BB(iT,ilam))
     &				/(Ca(ir,ilam,ig)+Cs(ir,ilam)*(1d0-g(ir,ilam)))
				tot2=tot2+ABS(freq(ilam)-freq(ilam+1))*(BB(iT+1,ilam)-BB(iT,ilam))
			enddo
		enddo
		Crw(ir)=tot2/Crw(ir)
		dorw(ir)=.false.
		if((R(ir+1)-R(ir))*Crw(ir).gt.factRW) dorw(ir)=.true.
	enddo

	Tc=2.7d0
	Ec=0d0
	EJv=0d0

	jr0=0

	docloud=.false.
	E0=Lplanet/real(NphotPlanet)
	do iphot=1,NphotPlanet
		call tellertje(iphot,NphotPlanet)
		call emit(specsurface,1,ilam,ig)
		call randomdirection(x,y,z)
		rr=R(ir0)
		x=x*rr
		y=y*rr
		z=z*rr
		jr=ir0
		onedge=.true.
1		call randomdirection(dx,dy,dz)
		goingup=((x*dx+y*dy+z*dz).gt.0d0)
		if(.not.goingup) goto 1
		do while(jr.le.nr)
			call travel(x,y,z,dx,dy,dz,jr,onedge,goingup,E,nscat,
     &			Ce(1:nr,ilam,ig),Ca(1:nr,ilam,ig),Cs(1:nr,ilam),Crw,g(1:nr,ilam),M(1:nr,ilam),dorw,1d0,Eabs,EJv)
			if(jr.le.nr) then
			if(jr.ne.jr0) then
				CabsL=0d0
				do ilam=1,nlam-1
					do ig=1,ng
						CabsL(ilam)=CabsL(ilam)+Ca(jr,ilam,ig)/real(ng)
					enddo
				enddo
				jr0=jr
			endif
			call increaseT(Tc(jr),T0,Ec(jr),E0*Eabs,CabsL,jr)
			iT1=Tc(jr)+1
			if(iT1.gt.nBB-1) iT1=nBB-1
			iT2=T0+1
			if(iT1.gt.nBB) iT1=nBB
			if(iT1.eq.iT2) iT2=iT1+1
			do ilam=1,nlam-1
				do ig=1,ng
					spec(ilam,ig)=Ca(jr,ilam,ig)*abs(BB(iT2,ilam)-BB(iT1,ilam))
				enddo
			enddo
			call emit(spec,ng,ilam,ig)
			call randomdirection(dx,dy,dz)
			goingup=((x*dx+y*dy+z*dz).gt.0d0)
			onedge=.false.
			Tc(jr)=T0
			endif
		enddo
	enddo

	EJv=EJv*E0
	chi2=0d0
	do ir=nr,1,-1
		CabsL=0d0
		do ilam=1,nlam-1
			do ig=1,ng
				CabsL(ilam)=CabsL(ilam)+Ca(ir,ilam,ig)/real(ng)
			enddo
		enddo
		call increaseT(T(ir),T0,EJv(ir),0d0,CabsL,ir)
c		if(ir.lt.nr) then
c			nabla=log(T(ir+1)/T0)/log(P(ir+1)/P(ir))
c			if(nabla.gt.(2d0/7d0)) then
c				nabla=2d0/7d0
c				T0=T(ir+1)/exp(nabla*log(P(ir+1)/P(ir)))
c			endif
c		endif
		chi2=chi2+((T(ir)-T0)/((T0+T(ir))*0.1))**2
		T(ir)=sqrt(T(ir)*T0)
	enddo
	chi2=chi2/real(nr)
	converged=.false.
	if(chi2.lt.10d0) converged=.true.
	call output("Chi2: " // trim(dbl2string(chi2,'(f7.2)')))

	deallocate(spec)
	deallocate(Ca)
	deallocate(Ce)
	deallocate(Cs)
	deallocate(g)
	deallocate(M)
	
	return
	end

	subroutine increaseT(T1,T2,E,Eabs,CabsL,ir)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ilam,ir,iTmax,iTmin,iT
	real*8 T1,T2,E,E0,Eabs,CabsL(nlam),Ca0,Cs0,V

	V=4d0*pi*(R(ir+1)**3-R(ir)**3)/3d0

	E=E+Eabs
	
	iTmin=1

	iTmax=nBB

	iT=0.5d0*(iTmax+iTmin)

	do while(abs(iTmax-iTmin).gt.1)
		E0=0d0
		do ilam=1,nlam-1
			E0=E0+ABS(freq(ilam)-freq(ilam+1))*CabsL(ilam)*BB(iT,ilam)
		enddo
		E0=E0*V
		if(E0.gt.E) then
			iTmax=iT
		else
			iTmin=iT
		endif
		iT=0.5d0*(iTmax+iTmin)
	enddo
	T2=real(iT)
	
	return
	end
	

	subroutine emit(spec,ng0,inu,ig)
	use GlobalSetup
	IMPLICIT NONE
	integer inu,ng0,j,ig,i
	real*8 spec(nlam,ng0),Ltot,Lr,Lt
	real*8 random

	Ltot=0d0
	do i=1,nlam-1
		do j=1,ng0
			Ltot=Ltot+ABS(freq(i)-freq(i+1))*spec(i,j)
		enddo
	enddo

	Lr=random(idum)*Ltot

	Lt=0d0
	do i=1,nlam-1
		do j=1,ng0
			Lt=Lt+ABS(freq(i)-freq(i+1))*spec(i,j)
			if(Lt.gt.Lr) then
				inu=i
				ig=j
				return
			endif
		enddo
	enddo

	return
	end
	
