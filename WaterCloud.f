	subroutine WaterCloud(ii)
	use GlobalSetup
	use Constants
	use AtomsModule
	use CloudModule
	use TimingModule
	IMPLICIT NONE
	real*8,allocatable :: x(:),vsed(:),dx(:),vth(:)
	real*8,allocatable :: Sn(:),mpart(:)
	real*8,allocatable :: An(:,:),y(:,:),Ma(:),Mb(:),Mc(:),CloudHp(:)
	real*8,allocatable :: at_ab(:,:)
	real*8,allocatable,save :: Sc(:),vthv(:),Aomp(:,:),xomp(:),IWORKomp(:),AB(:,:)
	real*8,allocatable :: tcinv(:,:),rho_av(:),densv(:,:),Kd(:),Kg(:),Km(:)
	real*8 dz,z12,z13,z12_2,z13_2,g,rr,mutot,npart,tot,lambda,densv_t,tot1,tot2,tot3
	integer info,i,j,iter,NN,NRHS,niter,ii,k,ihaze,kl,ku
	real*8 cs,eps,frac_nuc,m_nuc,tcoaginv,Dp,vmol,f,mm,err,maxerr
	real*8 Pv,molfracs_atoms0(N_atoms),NKn,Kzz_r(nr),vBM,scale,mixrat_bot
	integer,allocatable :: IWORK(:),ixv(:,:),ixc(:,:)
	real*8 sigmastar,Sigmadot,Pstar,gz,sigmamol,COabun,lmfp,fstick,kappa_cloud,fmin,rho_nuc
	logical ini,Tconverged
	character*500 cloudspecies(max(nclouds,1)),form
	real*8,allocatable :: CrV_prev0(:),CrT_prev0(:),CloudMR(:),xn_iter(:,:),xm_iter(:,:),xc_iter(:,:,:),xv_iter(:,:,:)
	logical,allocatable :: empty(:),docondense(:),liq(:)
	integer iCS,ir,nrdo,nconv,iconv
	real*8 logP(nr),logx(nr),dlogx(nr),SiSil,OSil,St,dztot
	real*8,allocatable :: logCloudP(:),CloudtauUV(:),CloudkappaUV(:)
	integer INCFD,IERR
	logical SKIP,freeflow
	real*8 time,tcrystinv,nucryst,Tcryst,fsed
	integer itime
!$OMP THREADPRIVATE(Sc,vthv,Aomp,xomp,IWORKomp,AB)

	logical dochemR(nr)

	call cpu_time(time)
	timecloud=timecloud-time
	call system_clock(itime)
	itimecloud=itimecloud-itime
	ctimecloud=ctimecloud+1

	nnr=(nr-1)*nr_cloud+1
	nCS=1
	if(.not.allocated(CloudP)) then
		allocate(CloudP(nnr))
		allocate(CloudT(nnr))
		allocate(CloudR(nnr))
		allocate(Clouddens(nnr))
		allocate(xv(nCS,nnr))
		allocate(xc(nCS,nnr))
		allocate(xn(nnr))
		allocate(xm(nnr))
		allocate(xa(nnr))
		allocate(rpart(nnr))
	endif
	allocate(Kd(nnr),Kg(nnr),Km(nnr))
	allocate(logCloudP(nnr))
	allocate(liq(nnr),CloudMR(nnr))
	allocate(CloudHp(nnr))
	
	niter=500
	nconv=20
	if(computeT) then
		if(nTiter.eq.1) then
			niter=20
			nconv=5
		else if(nTiter.le.4) then
			niter=100
			nconv=10
		endif
	endif

	allocate(densv(nnr,nCS),docondense(nCS))

	if(.not.allocated(ATP)) then
		allocate(ATP(nCS))
		allocate(BTP(nCS))
		allocate(rhodust(nCS))
		allocate(atoms_cloud(nCS,N_atoms))
		allocate(maxT(nCS))
		allocate(xv_bot(nCS))
		allocate(mu(nCS))
		allocate(CSname(nCS))
		allocate(CSnmol(nCS))
		allocate(ice(nCS))
		allocate(Tevap(nr,nCS))
	endif

	maxT=1d200

	atoms_cloud=0
	i=0
c H2O: 7
	i=i+1
	CSname(i)='H2O'
	ATP(i)=6511.
	BTP(i)=33.1
	atoms_cloud(i,1)=2
	atoms_cloud(i,5)=1
	rhodust(i)=1.0d0
	CSnmol(i)=1d0
	ice(i)=.true.

	nCS=i

	do i=1,nCS
		mu(i)=sum(mass_atoms(1:N_atoms)*atoms_cloud(i,1:N_atoms))/CSnmol(i)
	enddo

	mutot=0d0
	xv_bot=mixrat_r(1,1)
	mixrat_bot=mixrat_r(1,1)
	do i=1,nmol
		if(includemol(i)) mutot=mutot+mixrat_r(1,i)*Mmol(i)
	enddo
	xv_bot=xv_bot*mu*CSnmol/mutot

	Cloud(ii)%frac=0d0

	sigmastar=0.1
	Pstar=60d-6

	Pstar=Cloud(ii)%P
	sigmastar=log(Cloud(ii)%dP)

	fstick=1d0
	
	sigmamol=8d-15

	eps=1d-3

	Sigmadot=Cloud(ii)%Sigmadot
	
	select case(Cloud(ii)%hazetype)
		case("SOOT","soot","Soot")
			rho_nuc=1.00
			ihaze=10
		case("THOLIN","tholin","Tholin")
			rho_nuc=1.00
			ihaze=10
		case("optEC")
			rho_nuc=1.50
			ihaze=10
		case("SiC")
			rho_nuc=3.22
			ihaze=11
		case("CARBON","Carbon","carbon")
			rho_nuc=1.80
			ihaze=10
		case("CORRUNDUM","Corrundum","corrundum","Al2O3")
			rho_nuc=3.97
			ihaze=3
		case("IRON","Iron","iron","Fe")
			rho_nuc=7.87
			ihaze=8
		case("SiO")
			rho_nuc=2.18
			ihaze=4
		case("TiO2")
			rho_nuc=4.23
			ihaze=1
		case("Enstatite","enstatite","ENSTATITE")
			rho_nuc=3.20
			ihaze=5
		case default
			call output("hazetype unknown")
			stop
	end select

	m_nuc=4d0*pi*Cloud(ii)%rnuc**3*rho_nuc/3d0

	allocate(mpart(nnr))
	allocate(rho_av(nnr))
	allocate(y(nnr,5))
	allocate(Sn(nnr))
	allocate(vth(nnr))
	allocate(tcinv(niter,nnr))
	allocate(xn_iter(niter,nnr),xm_iter(niter,nnr))
	allocate(xc_iter(niter,nCS,nnr),xv_iter(niter,nCS,nnr))
	allocate(vsed(nnr))

	allocate(ixv(nCS,nnr))
	allocate(ixc(nCS,nnr))

	do iCS=1,nCS
		j=0
		do i=1,nnr
			j=j+1
			ixc(iCS,i)=j
			j=j+1
			ixv(iCS,i)=j
		enddo
	enddo

	rho_av=sum(rhodust)/real(nCS)

	do i=1,nnr
		CloudP(i)=10d0**(log10(P(1))+log10(P(nr)/P(1))*real(i-1)/real(nnr-1))
	enddo

	logP(1:nr)=-log(P(1:nr))
	logCloudP(1:nnr)=-log(CloudP(1:nnr))

	SKIP=.false.
	INCFD=1
	logx(1:nr)=log(T(1:nr))
	call DPCHIM(nr,logP,logx,dlogx,INCFD)
	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudT,IERR)
	CloudT(1:nnr)=exp(CloudT(1:nnr))

	SKIP=.false.
	INCFD=1
	logx(1:nr)=log(dens(1:nr))
	call DPCHIM(nr,logP,logx,dlogx,INCFD)
	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,Clouddens,IERR)
	Clouddens(1:nnr)=exp(Clouddens(1:nnr))

	SKIP=.false.
	INCFD=1
	logx(1:nr)=log(R(1:nr))
	call DPCHIM(nr,logP,logx,dlogx,INCFD)
	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudR,IERR)
	CloudR(1:nnr)=exp(CloudR(1:nnr))

	SKIP=.false.
	INCFD=1
	logx(1:nr)=log(mixrat_r(1:nr,1))
	call DPCHIM(nr,logP,logx,dlogx,INCFD)
	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudMR,IERR)
	CloudMR(1:nnr)=exp(CloudMR(1:nnr))

	if(complexKzz) then
		SKIP=.false.
		INCFD=1
		logx(1:nr)=log(Hp(1:nr))
		call DPCHIM(nr,logP,logx,dlogx,INCFD)
		call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudHp,IERR)
		CloudHp(1:nnr)=exp(CloudHp(1:nnr))
	endif
	if(Cloud(ii)%globalKzz.or.Cloud(ii)%Kzz.le.0d0) then
		SKIP=.false.
		INCFD=1
		logx(1:nr)=log10(Kzz_b(1:nr))
		call DPCHIM(nr,logP,logx,dlogx,INCFD)
		call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,Km,IERR)
		Km=10d0**Km
		do i=1,nnr
			Kd(i)=Km(i)
			if(complexKzz) then
				SKIP=.false.
				INCFD=1
				logx(1:nr)=log10(Kzz_g(1:nr))
				call DPCHIM(nr,logP,logx,dlogx,INCFD)
				call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,Kg,IERR)
				Kg=10d0**Kg
			else
				Kg(i)=Kd(i)
			endif
		enddo
	else
		Km=Cloud(ii)%Kzz
		Kd=Km
		Kg=Kd
	endif

	f=0.1d0

	allocate(empty(nnr))
	NN=2*nnr
	allocate(x(NN))
	allocate(IWORK(10*NN*NN))

	allocate(An(nnr,nnr))
	allocate(Ma(nnr))
	allocate(Mb(nnr))
	allocate(Mc(nnr))

	rpart=Cloud(ii)%rnuc
	xn=0d0
	xm=0d0
	xv=0d0
	xc=0d0
	do j=1,nCS
		do i=1,nnr
			xv(j,i)=xv_bot(j)
		enddo
	enddo
	tcinv=0d0
	
	docondense=.false.
	if(Cloud(ii)%hazetype.eq.'optEC') then
		allocate(CloudtauUV(nnr),CloudkappaUV(nnr))
		do ir=1,nr
			if(kappaUV0.gt.0d0) then
				tauUV(ir)=exp(-kappaUV0*1d6*P(ir)/grav(ir))
			else if(tauUV(ir).lt.0d0) then
				tauUV(ir)=exp(-1d6*P(ir)/grav(ir))
			endif
		enddo
		SKIP=.false.
		INCFD=1
		logx(1:nr)=tauUV(1:nr)
		call DPCHIM(nr,logP,logx,dlogx,INCFD)
		call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudtauUV,IERR)
		logx(1:nr)=kappaUV(1:nr)
		call DPCHIM(nr,logP,logx,dlogx,INCFD)
		call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudkappaUV,IERR)
		tot=0d0
	endif
	do i=1,nnr
		call PvapH2O(CloudT(i),densv(i,1),liq(i))
		densv(i,1)=1d6*densv(i,1)*(mu(1)*mp/(kb*CloudT(i)))
		if(densv(i,1).lt.Clouddens(i)*xv_bot(1)*CloudMR(i)/mixrat_bot) docondense(1)=.true.
		gz=Ggrav*Mplanet/CloudR(i)**2
		if(Cloud(ii)%hazetype.eq.'optEC') then
			Sn(i)=Clouddens(i)*CloudP(i)*CloudtauUV(i)
			if(i.eq.nnr) then
				tot=tot+abs(CloudR(i-1)-CloudR(i))*Sn(i)
			else if(i.eq.1) then
				tot=tot+abs(CloudR(i+1)-CloudR(i))*Sn(i)
			else
				tot=tot+abs(CloudR(i-1)-CloudR(i+1))*0.5*Sn(i)
			endif
		else
			Sn(i)=(Clouddens(i)*gz*Sigmadot/(sigmastar*CloudP(i)*1d6*sqrt(2d0*pi)))*exp(-log(CloudP(i)/Pstar)**2/(2d0*sigmastar**2))
		endif
	enddo
	if(Cloud(ii)%hazetype.eq.'optEC') then
		if(tot.gt.0d0) then
			Sn=Sn*scaleUV*Sigmadot/tot
		else
			Sn=0d0
		endif
	endif

!$OMP PARALLEL IF(useomp)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(nnr,NN)
	allocate(vthv(nnr))
	allocate(Sc(nnr))
	allocate(xomp(NN))
	allocate(IWORKomp(NN))
	allocate(Aomp(NN,NN))
	allocate(AB(7,NN))
!$OMP END PARALLEL
	
	iconv=0
c start the loop
	do iter=1,niter
	call tellertje(iter,niter)

	vsed=0d0
	do i=1,nnr
		cs=sqrt(kb*CloudT(i)/(2.3d0*mp))
		vth(i)=sqrt(8d0*kb*CloudT(i)/(pi*2.3d0*mp))
		if(Cloud(ii)%usefsed.and.iter.ne.1) then
			fsed=Cloud(ii)%fsed_alpha*exp((CloudR(i)-Rplanet)/(6d0*Cloud(ii)%fsed_beta*CloudHp(i)))+1d-2
			vsed(i)=-fsed*Kd(i)/CloudHp(i)
			rpart(i)=vsed(i)/(-sqrt(pi)*rho_av(i)*Ggrav*Mplanet/(2d0*Clouddens(i)*vth(i)*CloudR(i)**2))
			if(rpart(i).lt.Cloud(ii)%rnuc) then
				rpart(i)=Cloud(ii)%rnuc
				vsed(i)=-sqrt(pi)*rpart(i)*rho_av(i)*Ggrav*Mplanet/(2d0*Clouddens(i)*vth(i)*CloudR(i)**2)
			endif
		else
			vsed(i)=-sqrt(pi)*rpart(i)*rho_av(i)*Ggrav*Mplanet/(2d0*Clouddens(i)*vth(i)*CloudR(i)**2)
		endif
		mpart(i)=rho_av(i)*4d0*pi*rpart(i)**3/3d0
	enddo
	if(complexKzz) then
		do i=1,nnr
			St=rpart(i)*rho_av(i)*Km(i)/(vth(i)*Clouddens(i)*CloudHp(i)**2)
			Kd(i)=Km(i)/(1d0+St)
		enddo
	endif

	empty=.false.
	freeflow=Cloud(ii)%freeflow_nuc

	if(Cloud(ii)%usefsed.and.iter.ne.1) then
		do i=1,nnr
			tot=sum(xc(1:nCS,i))+xm(i)
			xn(i)=tot/mpart(i)
			if(Cloud(ii)%haze) then
				xm(i)=xn(i)*m_nuc
			else
				xm=0d0
			endif
		enddo
	else
	if(Cloud(ii)%hazetype.eq.'optEC') then
c equations for mass in vapor creating nuclii
		Ma=0d0
		Mb=0d0
		Mc=0d0
		x=0d0
		Mb(1)=1d0
		x(1)=mixrat_r(1,6)
		do i=2,nnr-1
			dztot=(CloudR(i+1)-CloudR(i-1))/2d0
			dz=(CloudR(i)-CloudR(i+1))
			Mc(i)=Mc(i)-(0.5d0*(Kg(i+1)*Clouddens(i+1)+Kg(i)*Clouddens(i))/dz)/dztot
			Mb(i)=Mb(i)+(0.5d0*(Kg(i+1)*Clouddens(i+1)+Kg(i)*Clouddens(i))/dz)/dztot
			dz=(CloudR(i-1)-CloudR(i))
			Ma(i-1)=Ma(i-1)-0.5d0*(Kg(i-1)*Clouddens(i-1)+Kg(i)*Clouddens(i))/dz/dztot
			Mb(i)=Mb(i)+(0.5d0*(Kg(i-1)*Clouddens(i-1)+Kg(i)*Clouddens(i))/dz)/dztot

			Mb(i)=Mb(i)-Sn(i)*mutot/Mmol(6)
			x(i)=0d0
		enddo
		i=nnr
		j=j+1
		dz=CloudR(i)-CloudR(i-1)
		Mb(nnr)=Kg(i)/dz
		Ma(nnr-1)=-Kg(i)/dz
		x(nnr)=0d0

		NRHS=1
c		call DGESV( nnr, NRHS, An, nnr, IWORK, x, nnr, info )
		call dgtsv(nnr,NRHS,Ma,Mb,Mc,x,nnr,info)
	
		do i=1,nnr
			if(x(i).lt.0d0) x(i)=0d0
		enddo
		xa(1:nnr)=x(1:nnr)
	else
		xa=1d0
	endif

c equations for number of Nuclii
	Ma=0d0
	Mb=0d0
	Mc=0d0
	x=0d0
	if(freeflow) then
		i=1
		dztot=(CloudR(i+1)-CloudR(i))
		dz=(CloudR(i)-CloudR(i+1))
		Mc(i)=Mc(i)-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		Mb(i)=Mb(i)+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		Mb(i)=Mb(i)+Clouddens(i)*vsed(i)/dztot
		x(i)=-Sn(i)*xa(i)
c coagulation
		if(Cloud(ii)%coagulation) then
			npart=xn(i)*Clouddens(i)
			lmfp=2.3d0*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
			vmol=0.5d0*lmfp*vth(i)
			Dp=kb*CloudT(i)/(6d0*pi*rpart(i)*vmol*Clouddens(i))
			tcoaginv=sqrt(pi)*3d0*(sum(xc(1:nCS,i))+xm(i))*Ggrav*Mplanet/(8d0*vth(i)*CloudR(i)**2)
			vBM=sqrt(16d0*kb*CloudT(i)/(pi*mpart(i)))
			if(Dp/rpart(i).lt.vBM) vBM=Dp/rpart(i)
			tcoaginv=tcoaginv+2d0*pi*rpart(i)**2*npart*vBM
			if(.not.tcoaginv.gt.0d0) tcoaginv=0d0
			tcinv(iter,i)=tcoaginv
			tcoaginv=sum(tcinv(1:iter,i))/real(iter)
			Mb(i)=Mb(i)-Clouddens(i)*tcoaginv
		endif
	else
		Mb(1)=1d0
		x(1)=Cloud(ii)%xm_bot
	endif
	do i=2,nnr-1
		dztot=(CloudR(i+1)-CloudR(i-1))/2d0
		dz=(CloudR(i)-CloudR(i+1))
		Mc(i)=Mc(i)-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		Mb(i)=Mb(i)+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		dz=(CloudR(i-1)-CloudR(i))
		Ma(i-1)=Ma(i-1)-0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz/dztot
		Mb(i)=Mb(i)+(0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz)/dztot
		Mb(i)=Mb(i)+Clouddens(i)*vsed(i)/dztot

		x(i)=-Sn(i)*xa(i)

c coagulation
		if(Cloud(ii)%coagulation) then
			npart=xn(i)*Clouddens(i)
			lmfp=2.3d0*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
			vmol=0.5d0*lmfp*vth(i)
			Dp=kb*CloudT(i)/(6d0*pi*rpart(i)*vmol*Clouddens(i))

c			tcoaginv=2d0*npart*pi*rpart(i)**2*abs(vsed(i))
c rewritten for better convergence
			tcoaginv=sqrt(pi)*3d0*(sum(xc(1:nCS,i))+xm(i))*Ggrav*Mplanet/(8d0*vth(i)*CloudR(i)**2)

			vBM=sqrt(16d0*kb*CloudT(i)/(pi*mpart(i)))
			if(Dp/rpart(i).lt.vBM) vBM=Dp/rpart(i)
			tcoaginv=tcoaginv+2d0*pi*rpart(i)**2*npart*vBM

			if(.not.tcoaginv.gt.0d0) tcoaginv=0d0

			tcinv(iter,i)=tcoaginv
			
c			call computemedian(tcinv(1:iter,i),iter,tcoaginv)
			tcoaginv=sum(tcinv(1:iter,i))/real(iter)

			Mb(i)=Mb(i)-Clouddens(i)*tcoaginv
		endif
	enddo
	dz=CloudR(i)-CloudR(i-1)
	Mb(nnr)=Kd(i)/dz-vsed(i)
	Ma(nnr-1)=-Kd(i)/dz
	x(nnr)=0d0

	NRHS=1
c	call DGESV( nnr, NRHS, An, nnr, IWORK, x, nnr, info )
	call dgtsv(nnr,NRHS,Ma,Mb,Mc,x,nnr,info)

	do i=1,nnr
		if(x(i).lt.0d0) x(i)=0d0
	enddo
	xn(1:nnr)=x(1:nnr)

	if(Cloud(ii)%coagulation.and.Cloud(ii)%haze) then
c equations for mass in Nuclii
		Ma=0d0
		Mb=0d0
		Mc=0d0
		x=0d0
		if(freeflow) then
			i=1
			dztot=(CloudR(i+1)-CloudR(i))
			dz=(CloudR(i)-CloudR(i+1))
			Mc(i)=Mc(i)-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
			Mb(i)=Mb(i)+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
			Mb(i)=Mb(i)+Clouddens(i)*vsed(i)/dztot
			x(i)=-Sn(i)
		else
			Mb(1)=1d0
			x(1)=Cloud(ii)%xm_bot
		endif
		do i=2,nnr-1
			dztot=(CloudR(i+1)-CloudR(i-1))/2d0
			dz=(CloudR(i)-CloudR(i+1))
			Mc(i)=Mc(i)-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
			Mb(i)=Mb(i)+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
			dz=(CloudR(i-1)-CloudR(i))
			Ma(i-1)=Ma(i-1)-0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz/dztot
			Mb(i)=Mb(i)+(0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz)/dztot
			Mb(i)=Mb(i)+Clouddens(i)*vsed(i)/dztot

			x(i)=-Sn(i)
		enddo
		i=nnr
		j=j+1
		dz=CloudR(i)-CloudR(i-1)
		Mb(nnr)=Kd(i)/dz-vsed(i)
		Ma(nnr-1)=-Kd(i)/dz
		x(nnr)=0d0

		NRHS=1
c		call DGESV( nnr, NRHS, An, nnr, IWORK, x, nnr, info )
		call dgtsv(nnr,NRHS,Ma,Mb,Mc,x,nnr,info)
	
		do i=1,nnr
			if(x(i).lt.0d0) x(i)=0d0
		enddo
		xm(1:nnr)=x(1:nnr)
	else
		xm=xn
	endif

	do i=1,nnr
		if(.not.xm(i).gt.0d0.or..not.xn(i).gt.0d0) then
			xn(i)=0d0
			xm(i)=0d0
			empty(i)=.true.
		endif
		if(xn(i).gt.xm(i)) then
			xn(i)=xm(i)
		endif
	enddo

	if(.not.Cloud(ii)%haze) then
		xm=0d0
	endif
	xn(1:nnr)=xn(1:nnr)/m_nuc

	endif

	freeflow=Cloud(ii)%freeflow_con
	iCS=1

	do i=1,nnr
		cs=sqrt(kb*CloudT(i)/(2.3d0*mp))
		vthv(i)=sqrt(8d0*kb*CloudT(i)/(pi*mu(iCS)*mp))
		Sc(i)=min(vthv(i)*rpart(i),kb*CloudT(i)*sqrt(8d0*kb*CloudT(i)/(pi*2.3*mp))/(3d0*CloudP(i)*1d6*8e-15))
     &				*4d0*pi*rpart(i)*Clouddens(i)
		Sc(i)=fstick*Sc(i)
	enddo

c equations for material
	Aomp=0d0
	xomp=0d0
	j=0
	i=1
	j=j+1
	if(freeflow) then
c assume continuous flux at the bottom (dF/dz=Sc=0)
		xomp(j)=0d0

		dztot=(CloudR(i+1)-CloudR(i))
		dz=(CloudR(i)-CloudR(i+1))
		Aomp(j,ixc(iCS,i+1))=Aomp(j,ixc(iCS,i+1))-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+Clouddens(i)*vsed(i)/dztot

		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))+Sc(i)*xn(i)*Clouddens(i)!*CloudMR(i)/mixrat_bot
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-Sc(i)*densv(i,iCS)/mpart(i)
	else
		Aomp(j,ixc(iCS,i))=1d0
		xomp(j)=0d0
	endif
	j=j+1
	Aomp(j,ixv(iCS,i))=1d0
	xomp(j)=xv_bot(iCS)
	do i=2,nnr-1
		j=j+1

		if(.not.empty(i)) then

		dztot=(CloudR(i+1)-CloudR(i-1))/2d0
		dz=(CloudR(i)-CloudR(i+1))
		Aomp(j,ixc(iCS,i+1))=Aomp(j,ixc(iCS,i+1))-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		dz=(CloudR(i-1)-CloudR(i))
		Aomp(j,ixc(iCS,i-1))=Aomp(j,ixc(iCS,i-1))-0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz/dztot
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+(0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz)/dztot
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+Clouddens(i)*vsed(i)/dztot

		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))+Sc(i)*xn(i)*Clouddens(i)!*CloudMR(i)/mixrat_bot
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-Sc(i)*densv(i,iCS)/mpart(i)
		else
		Aomp(j,ixc(iCS,i))=1d0
		endif

		j=j+1

		dztot=(CloudR(i+1)-CloudR(i-1))/2d0
		dz=(CloudR(i)-CloudR(i+1))
		Aomp(j,ixv(iCS,i+1))=Aomp(j,ixv(iCS,i+1))-(0.5d0*(Kg(i+1)*Clouddens(i+1)+Kg(i)*Clouddens(i))/dz)/dztot
		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))+(0.5d0*(Kg(i+1)*Clouddens(i+1)+Kg(i)*Clouddens(i))/dz)/dztot
		dz=(CloudR(i-1)-CloudR(i))
		Aomp(j,ixv(iCS,i-1))=Aomp(j,ixv(iCS,i-1))-0.5d0*(Kg(i-1)*Clouddens(i-1)+Kg(i)*Clouddens(i))/dz/dztot
		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))+(0.5d0*(Kg(i-1)*Clouddens(i-1)+Kg(i)*Clouddens(i))/dz)/dztot

		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))-Sc(i)*xn(i)*Clouddens(i)!*CloudMR(i)/mixrat_bot
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+Sc(i)*densv(i,iCS)/mpart(i)
	enddo
	i=nnr
	dz=CloudR(i)-CloudR(i-1)
	j=j+1
	Aomp(j,ixc(iCS,i))=Kd(i)/dz-vsed(i)
	Aomp(j,ixc(iCS,i-1))=-Kd(i)/dz
	xomp(j)=0d0!-Mc_top/Clouddens(i)

	j=j+1
	Aomp(j,ixv(iCS,i))=Kg(i)/dz
	Aomp(j,ixv(iCS,i-1))=-Kg(i)/dz
	xomp(j)=0d0!Mc_top/Clouddens(i)

10	continue
	NRHS=1
	info=0
c Use Band matrix algorithm
	KL=2
	KU=2
	AB=0d0
	do j=1,NN
		do i=max(1,j-KU),min(j+KL,NN)
			AB(KL+KU+1+i-j,j) = Aomp(i,j)
		enddo
	enddo
	j=7
	call DGBSV(NN,KL,KU,NRHS,AB,j,IWORKomp,xomp,NN,INFO)	

	do i=1,NN
		if(.not.xomp(i).gt.0d0) xomp(i)=0d0
	enddo

	do i=1,nnr
		if(empty(i)) xomp(ixc(iCS,i))=0d0
		xc(iCS,i)=xomp(ixc(iCS,i))
		xv(iCS,i)=xomp(ixv(iCS,i))
	enddo
	do i=1,nnr
		if(xc(iCS,i).lt.0d0) xc(iCS,i)=0d0
		if(xv(iCS,i).lt.0d0) xv(iCS,i)=0d0
	enddo
	xc_iter(iter,1:nCS,1:nnr)=xc(1:nCS,1:nnr)
	xv_iter(iter,1:nCS,1:nnr)=xv(1:nCS,1:nnr)
	xn_iter(iter,1:nnr)=xn(1:nnr)
	xm_iter(iter,1:nnr)=xm(1:nnr)

	maxerr=0d0
	do i=1,nnr
		tot=xm(i)/rho_nuc
		do iCS=1,nCS
			tot=tot+xc(iCS,i)/rhodust(iCS)
		enddo
		if(tot.gt.0d0) then
			rho_av(i)=(sum(xc(1:nCS,i))+xm(i))/tot
		else
			rho_av(i)=sum(rhodust(1:nCS))/real(nCS)
		endif
		tot=sum(xc(1:nCS,i))+xm(i)
		if(xn(i).gt.0d0) then
			rr=(3d0*(tot/xn(i))/(4d0*pi*rho_av(i)))**(1d0/3d0)
			if(.not.rr.ge.Cloud(ii)%rnuc) then
				rr=Cloud(ii)%rnuc
				xn(i)=(3d0*(tot/(rr**3))/(4d0*pi*rho_av(i)))
			endif
		else
			rr=Cloud(ii)%rnuc
			xn(i)=(3d0*(tot/(rr**3))/(4d0*pi*rho_av(i)))
		endif
		err=abs(rr-rpart(i))/(rr+rpart(i))
		if(err.gt.maxerr) maxerr=err
		rpart(i)=sqrt(rr*rpart(i))
	enddo
	if(maxerr.lt.eps) then
		iconv=iconv+1
		if(iconv.gt.nconv) exit
	else
		iconv=0
	endif
	enddo
c end the loop
	if(iter.gt.niter) then
c		if(iconv.eq.0) print*,'Cloud formation not converged: ',maxerr
		iter=niter
	endif
	xn=0d0
	xm=0d0
	xc=0d0
	xv=0d0
	do i=iter-nconv+1,iter
		xn(1:nnr)=xn(1:nnr)+xn_iter(i,1:nnr)/real(nconv)
		xm(1:nnr)=xm(1:nnr)+xm_iter(i,1:nnr)/real(nconv)
		xc(1:nCS,1:nnr)=xc(1:nCS,1:nnr)+xc_iter(i,1:nCS,1:nnr)/real(nconv)
		xv(1:nCS,1:nnr)=xv(1:nCS,1:nnr)+xv_iter(i,1:nCS,1:nnr)/real(nconv)
	enddo
	do i=1,nnr
		tot=xm(i)/rho_nuc
		do iCS=1,nCS
			tot=tot+xc(iCS,i)/rhodust(iCS)
		enddo
		if(tot.gt.0d0) then
			rho_av(i)=(sum(xc(1:nCS,i))+xm(i))/tot
		else
			rho_av(i)=sum(rhodust(1:nCS))/real(nCS)
		endif
		tot=sum(xc(1:nCS,i))+xm(i)
		if(xn(i).gt.0d0) then
			rr=(3d0*(tot/xn(i))/(4d0*pi*rho_av(i)))**(1d0/3d0)
			if(.not.rr.ge.Cloud(ii)%rnuc) then
				rr=Cloud(ii)%rnuc
				xn(i)=(3d0*(tot/(rr**3))/(4d0*pi*rho_av(i)))
			endif
		else
			rr=Cloud(ii)%rnuc
			xn(i)=(3d0*(tot/(rr**3))/(4d0*pi*rho_av(i)))
		endif
		rpart(i)=rr
	enddo


!$OMP PARALLEL IF(useomp)
!$OMP& DEFAULT(NONE)
	deallocate(vthv)
	deallocate(Sc)
	deallocate(xomp)
	deallocate(IWORKomp)
	deallocate(Aomp)
	deallocate(AB)
!$OMP END PARALLEL

	allocate(dx(nnr))
	logP(1:nr)=-log(P(1:nr))
	logCloudP(1:nnr)=-log(CloudP(1:nnr))

	x(1:nnr)=xm(1:nnr)*Clouddens(1:nnr)
	do iCS=1,nCS
		x(1:nnr)=x(1:nnr)+xc(iCS,1:nnr)*Clouddens(1:nnr)
	enddo
	call regridarray(logCloudP,x,nnr,logP,cloud_dens(1:nr,ii),nr)

	x(1:nnr)=rpart(1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%rv(1:nr),nr)

	Cloud(ii)%frac(1:nr,1:20)=0d0

c H2O
	x(1:nnr)=xc(1,1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,12),nr)
c Seed particles
	x(1:nnr)=xm(1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,15),nr)

	if(writefiles) then
		if(do3D) then
			open(unit=20,file=trim(outputdir) // '/cloudstructure' // trim(int2string(i3D,'(i0.4)')) // '.dat',
     &             FORM="FORMATTED",ACCESS="STREAM")
		else
			open(unit=20,file=trim(outputdir) // '/cloudstructure.dat',FORM="FORMATTED",ACCESS="STREAM")
		endif
		form='("#",a18,a19,a19,' // trim(int2string(nCS,'(i4)')) // 'a23,a23,a19,a19,a19)'
		write(20,form) "P[bar]","dens[g/cm^3]","xn",(trim(CSname(i)),i=1,nCS),"H2O[v]","r[micron]","T[K]","Jstar"
		form='(es19.7E3,es19.7E3,es19.7E3,' // trim(int2string(nCS,'(i4)')) // 'es23.7E3,es23.7E3,es19.7E3,es19.7E3,es19.7E3)'
		do i=1,nnr
			write(20,form) CloudP(i),Clouddens(i),xm(i),xc(1:nCS,i),xv(1,i),rpart(i),CloudT(i),
     &						Sn(i)/m_nuc
		enddo
		close(unit=20)
	endif

	do i=1,nr
		tot=sum(Cloud(ii)%frac(i,1:20))
		if(tot.gt.0d0) then
			Cloud(ii)%frac(i,1:20)=Cloud(ii)%frac(i,1:20)/tot
		else
			Cloud(ii)%frac(i,1:20)=1d0/20d0
			cloud_dens(i,ii)=0d0
		endif
	enddo


	x(1:nnr)=xv(1,1:nnr)
	call regridarray(logCloudP,x,nnr,logP,logx,nr)
	mixrat_r(2:nr,1)=logx(2:nr)*mutot/(mu(1)*CSnmol(1))

	if(Cloud(ii)%hazetype.eq.'optEC') then
c Seed vapor
		x(1:nnr)=xa(1:nnr)
		call regridarray(logCloudP,x,nnr,logP,logx,nr)
		mixrat_r(2:nr,6)=logx(2:nr)
	endif
	
	do i=1,nr
		tot=0d0
		do j=1,nmol
			if(mixrat_r(i,j).gt.0d0) tot=tot+mixrat_r(i,j)
		enddo
		if(tot.gt.0d0) mixrat_r(i,1:nmol)=mixrat_r(i,1:nmol)/tot
	enddo

	call cpu_time(time)
	timecloud=timecloud+time
	call system_clock(itime)
	itimecloud=itimecloud+itime

	deallocate(densv,docondense)
	deallocate(mpart)
	deallocate(rho_av)
	deallocate(y)
	deallocate(Sn)
	deallocate(vth)
	deallocate(tcinv,xn_iter,xm_iter,xc_iter,xv_iter)
	deallocate(vsed)
	deallocate(ixv)
	deallocate(ixc)
	deallocate(x,dx)
	deallocate(empty)
	deallocate(IWORK)
	deallocate(An)
	deallocate(logCloudP)
	deallocate(Kd,Kg,Km)
	deallocate(Ma,Mb,Mc)
	if(Cloud(ii)%hazetype.eq.'optEC') deallocate(CloudtauUV,CloudkappaUV)

	return
	end
	
	

	subroutine PvapH2O(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0_i,c1_i,c2_i,c3_i
	real*8 c0_l,c1_l,c2_l,c3_l,c4_l
	parameter(	
     &	c0_i=6111.5d0,
     &	c1_i=23.036d0,
     &	c2_i=-333.7d0,
     &	c3_i=279.82d0,
     &	c0_l=2.98605d+01,
     &	c1_l=-3.15220d+03,
     &	c2_l=-7.30370d+00,
     &	c3_l=2.42470d-09,
     &	c4_l=1.80900d-06)

	if(T.gt.747d0) then
		Pi=1d200
	else
		Pi=c0_i*exp((c1_i*(T-273.15) + (T-273.15)**2/c2_i)/((T-273.15) + c3_i))
		Pi=1d-6*Pi
	endif
		
	if(T.gt.747d0) then
		Pl=1d200
	else if(T.lt.198.0) then
		Pl=1d200
	else
		Pl=c0_l+c1_l/T+c2_l*log10(T)+c3_l*T+c4_l*T**2
		Pl=0.00133322368*10d0**Pl
	endif
	if(Pi.lt.Pl) then
		Pvap=Pi
		liquid=.false.
	else
		Pvap=Pl
		liquid=.true.
	endif
	
	return
	end
	