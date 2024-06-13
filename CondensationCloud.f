	subroutine CondensationCloud(ii)
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
	real*8,allocatable,save :: ScN(:),ScS(:),vthv(:),Aomp(:,:),xomp(:),IWORKomp(:),AB(:,:)
	real*8,allocatable :: tcinv(:,:),rho_av(:),densv(:,:),Kd(:),Kg(:),Km(:)
	real*8 dz,z12,z13,z12_2,z13_2,g,rr,mutot,npart,tot,lambda,densv_t,tot1,tot2,tot3
	integer info,i,j,iter,NN,NRHS,niter,ii,k,ihaze,kl,ku
	real*8 cs,eps,frac_nuc,m_nuc,tcoaginv,Dp,vmol,f,mm,err,maxerr,dztot
	real*8 Pv,molfracs_atoms0(N_atoms),NKn,Kzz_r(nr),vBM,scale
	integer,allocatable :: IWORK(:),ixv(:,:),ixc(:,:)
	real*8 sigmastar,Sigmadot,Pstar,gz,sigmamol,COabun,lmfp,fstick,kappa_cloud,fmin,rho_nuc
	logical ini,Tconverged
	character*500 cloudspecies(max(nclouds,1)),form
	real*8,allocatable :: CrV_prev0(:),CrT_prev0(:),xn_iter(:,:),xm_iter(:,:),xc_iter(:,:,:),xv_iter(:,:,:)
	logical,allocatable :: empty(:),docondense(:)
	integer iCS,ir,nrdo,iconv,nconv
	real*8 logP(nr),logx(nr),dlogx(nr),SiSil,OSil,St,fsed
	real*8,allocatable :: logCloudP(:),CloudtauUV(:),CloudkappaUV(:),ScTot(:),xs(:)
	integer INCFD,IERR
	logical SKIP,freeflow,liq
	real*8 time
	integer itime
!$OMP THREADPRIVATE(ScN,ScS,vthv,Aomp,xomp,IWORKomp,AB)

	logical dochemR(nr)

	call cpu_time(time)
	timecloud=timecloud-time
	call system_clock(itime)
	itimecloud=itimecloud-itime
	ctimecloud=ctimecloud+1

	nnr=(nr-1)*nr_cloud+1
	nCS=Cloud(ii)%nmat
	if(Cloud(ii)%haze) nCS=nCS-1
	if(.not.allocated(CloudP)) then
		allocate(CloudP(nnr))
		allocate(CloudT(nnr))
		allocate(CloudR(nnr))
		allocate(Clouddens(nnr))
		allocate(xv(nCS,nnr))
		allocate(xc(nCS,nnr))
		allocate(xn(nnr))
		allocate(xm(nnr))
		allocate(xnv(nnr))
		allocate(rpart(nnr))
	endif
	allocate(Kd(nnr),Kg(nnr),Km(nnr))
	allocate(logCloudP(nnr))
	allocate(CloudHp(nnr))
	allocate(xs(nnr),ScTot(nnr))
	ScTot=0d0
	
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

	call SetAbun

	COabun=min(molfracs_atoms(3),molfracs_atoms(5))
	molfracs_atoms(3)=molfracs_atoms(3)-COabun
	molfracs_atoms(5)=molfracs_atoms(5)-COabun

	ATP=28030d0		!6.89d4
	BTP=12.471d0		!37.8

	maxT=1d200

	atoms_cloud=0
	
	do i=1,nCS
		select case(Cloud(ii)%condensate(i))
			case('TiO2')
c TiO2: 1
				CSname(i)='TiO2'
				ATP(i)=74734.7	! Value from Woitke & Helling 2004 (differences are so small that the old version is kept for backward compatibility)
				BTP(i)=35.8027
				atoms_cloud(i,15)=1
				atoms_cloud(i,5)=2
				rhodust(i)=7.0d0	! moet ik nog checken
				CSnmol(i)=1d0
				ice(i)=.false.
			case('VO')
c VO: 2
				CSname(i)='VO'
				ATP(i)=74734.7	! Value from Woitke & Helling 2004 for TiO2
				BTP(i)=35.8027
				atoms_cloud(i,16)=1
				atoms_cloud(i,5)=1
				rhodust(i)=7.0d0	! moet ik nog checken
				CSnmol(i)=1d0
				ice(i)=.false.
			case('Al2O3','CORRUNDUM')
c Al2O3: 3
				CSname(i)='Al2O3'
				ATP(i)=77365.
				BTP(i)=39.3
				atoms_cloud(i,8)=2
				atoms_cloud(i,5)=3
				rhodust(i)=4.0d0
				CSnmol(i)=2d0
				ice(i)=.false.
			case('SiO2','QUARTZ')
c SiO2: 4
				CSname(i)='SiO2'
				ATP(i)=69444.
				BTP(i)=33.1
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=2
				rhodust(i)=2.2d0
				CSnmol(i)=1d0
				ice(i)=.false.
			case('SILICATE')
c Mg/Ca Silicates: 5
				ATP(i)=68908.	! MgSiO3 for now
				BTP(i)=38.1
c ===========================================================
c = New version, using Ca and Mg to condense =========
c ===========================================================
				atoms_cloud(i,1:N_atoms)=0
				atoms_cloud(i,7)=molfracs_atoms(7)
				atoms_cloud(i,14)=molfracs_atoms(14)
				tot=sum(atoms_cloud(i,1:N_atoms))
				atoms_cloud(i,1:N_atoms)=atoms_cloud(i,1:N_atoms)/tot
				atoms_cloud(i,5)=1

				write(CSname(i),'("Mg",f3.1,"Ca",f3.1,"SiO",f3.1)') atoms_cloud(i,7),atoms_cloud(i,14),atoms_cloud(i,5)+2
c ===========================================================
				rhodust(i)=2.0d0
				CSnmol(i)=2d0
				ice(i)=.false.
			case('H2O','WATER')
c H2O: 7
				CSname(i)='H2O'
				ATP(i)=6511.
				BTP(i)=33.1
				atoms_cloud(i,1)=2
				atoms_cloud(i,5)=1
				rhodust(i)=1.0d0
				CSnmol(i)=1d0
				ice(i)=.true.
			case('Fe','IRON')
c Fe: 8
				CSname(i)='Fe'
				ATP(i)=48354.
				BTP(i)=29.2
				atoms_cloud(i,17)=1
				rhodust(i)=7.8d0
				CSnmol(i)=1d0
				ice(i)=.false.
			case('FeS','TROILITE')
c FeS: 9
				CSname(i)='FeS'
				ATP(i)=48354.
				BTP(i)=29.2
				atoms_cloud(i,17)=1
				atoms_cloud(i,11)=1
				rhodust(i)=4.8d0
				maxT(i)=680d0
				CSnmol(i)=1d0
				ice(i)=.false.
			case('C','CARBON')
c C: 10
				CSname(i)='C'
				ATP(i)=93646.
				BTP(i)=36.7
				atoms_cloud(i,3)=1
				rhodust(i)=1.8d0
				CSnmol(i)=1d0
				ice(i)=.false.
			case('SiC')
c SiC: 11
				CSname(i)='SiC'
				ATP(i)=78462.
				BTP(i)=37.8
				atoms_cloud(i,9)=1
				atoms_cloud(i,3)=1
				rhodust(i)=2.5d0
				CSnmol(i)=1d0
				ice(i)=.false.
			case("NH3","AMONIA")
c NH3
				CSname(i)='NH3'
				atoms_cloud(i,1)=3
				atoms_cloud(i,4)=1
				rhodust(i)=0.87d0
				CSnmol(i)=1d0
				ice(i)=.true.
			case("NaCl")
c NaCl
				CSname(i)='NaCl'
				atoms_cloud(i,6)=1
				atoms_cloud(i,12)=1
				rhodust(i)=2.17d0
				CSnmol(i)=1d0
				ice(i)=.true.
			case("KCl")
c KCl
				CSname(i)='KCl'
				atoms_cloud(i,13)=1
				atoms_cloud(i,12)=1
				rhodust(i)=1.99d0
				CSnmol(i)=1d0
				ice(i)=.true.
			case default
				call output("Unknown condensate")
				stop
		end select
	enddo

	docondense=.true.

	do i=1,nCS
		mu(i)=sum(mass_atoms(1:N_atoms)*atoms_cloud(i,1:N_atoms))/CSnmol(i)
	enddo

	mutot=0d0
	xv_bot=1d200
	do i=1,N_atoms
		mutot=mutot+mass_atoms(i)*molfracs_atoms(i)
	enddo
	if(dochemistry) then
		do iCS=1,nCS
			if(docondense(iCS)) then
				do k=1,N_atoms
					if(atoms_cloud(iCS,k).gt.0d0) then
						f=molfracs_atoms(k)/atoms_cloud(iCS,k)
						if(f.lt.xv_bot(iCS)) then
							xv_bot(iCS)=f
						endif
					endif
				enddo
				do k=1,N_atoms
					molfracs_atoms(k)=molfracs_atoms(k)-xv_bot(iCS)*atoms_cloud(iCS,k)
					if(molfracs_atoms(k).lt.0d0) molfracs_atoms(k)=0d0
				enddo
			else
				xv_bot(iCS)=0d0
			endif
		enddo
	else
		do iCS=1,nCS
			xv_bot(iCS)=Cloud(ii)%xv_bot(iCS)
		enddo
		mutot=0d0
		do i=1,nmol
			if(includemol(i)) mutot=mutot+mixrat_r(1,i)*Mmol(i)
		enddo
		do iCS=1,nCS
			if(CSname(iCS).eq."H2O") then
				xv_bot(iCS)=mixrat_r(1,1)*mu(iCS)*CSnmol(iCS)/mutot
			endif
		enddo
	endif
	
	molfracs_atoms0=molfracs_atoms
	xv_bot=xv_bot*mu*CSnmol/mutot

	if(Cloud(ii)%rainout) then
		densv(1,1:nCS)=(mu*mp/(kb*T(1)))*exp(BTP(1:nCS)-ATP(1:nCS)/T(1))
		do iCS=1,nCS
			if(CSname(iCS).eq."H2O") then
				call PvapH2O(T(1),densv(1,iCS),liq)
				densv(1,iCS)=1d6*densv(1,iCS)*(mu(iCS)*mp/(kb*T(1)))
			else if(CSname(iCS).eq."NH3") then
				call PvapNH3(T(1),densv(1,iCS),liq)
				densv(1,iCS)=1d6*densv(1,iCS)*(mu(iCS)*mp/(kb*T(1)))
			else if(CSname(iCS).eq."KCl") then
				call PvapKCl(T(1),densv(1,iCS),liq)
				densv(1,iCS)=1d6*densv(1,iCS)*(mu(iCS)*mp/(kb*T(1)))
			else if(CSname(iCS).eq."NaCl") then
				call PvapNaCl(T(1),densv(1,iCS),liq)
				densv(1,iCS)=1d6*densv(1,iCS)*(mu(iCS)*mp/(kb*T(1)))
			endif
			if(T(1).gt.maxT(iCS)) densv(1,iCS)=densv(1,iCS)+(mu(iCS)*mp/(kb*T(1)*10d0))*exp(BTP(iCS)-ATP(iCS)/(T(1)*10d0))
		enddo
		do iCS=1,nCS
			if(densv(1,iCS).lt.dens(1)*xv_bot(iCS)) then
				call output("Rainout of: " // trim(CSname(iCS)) // "(" //
     &		trim(dbl2string(100d0*(1d0-densv(1,iCS)/(dens(1)*xv_bot(iCS))),'(f5.1)')) // " %)")
				xv_bot(iCS)=densv(1,iCS)/dens(1)
			endif
		enddo
	endif

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
				tauUV(ir)=kappaUV0*1d6*P(ir)/grav(ir)
			else if(tauUV(ir).lt.0d0) then
				tauUV(ir)=1d6*P(ir)/grav(ir)
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
		densv(i,1:nCS)=(mu*mp/(kb*CloudT(i)))*exp(BTP(1:nCS)-ATP(1:nCS)/CloudT(i))
		do iCS=1,nCS
			if(CSname(iCS).eq."H2O") then
				call PvapH2O(CloudT(i),densv(i,iCS),liq)
				densv(i,iCS)=1d6*densv(i,iCS)*(mu(iCS)*mp/(kb*CloudT(i)))
			else if(CSname(iCS).eq."NH3") then
				call PvapNH3(CloudT(i),densv(i,iCS),liq)
				densv(i,iCS)=1d6*densv(i,iCS)*(mu(iCS)*mp/(kb*CloudT(i)))
			else if(CSname(iCS).eq."NaCl") then
				call PvapNaCl(CloudT(i),densv(i,iCS),liq)
				densv(i,iCS)=1d6*densv(i,iCS)*(mu(iCS)*mp/(kb*CloudT(i)))
			else if(CSname(iCS).eq."KCl") then
				call PvapKCl(CloudT(i),densv(i,iCS),liq)
				densv(i,iCS)=1d6*densv(i,iCS)*(mu(iCS)*mp/(kb*CloudT(i)))
			endif
			if(cloudT(i).gt.maxT(iCS)) densv(i,iCS)=densv(i,iCS)+
     &                    (mu(iCS)*mp/(kb*CloudT(i)*10d0))*exp(BTP(iCS)-ATP(iCS)/(CloudT(i)*10d0))
		enddo
		do iCS=1,nCS
			if(densv(i,iCS).lt.Clouddens(i)*xv_bot(iCS)) docondense(iCS)=.true.
		enddo
		gz=Ggrav*Mplanet/CloudR(i)**2
		if(Cloud(ii)%hazetype.eq.'optEC') then
			Sn(i)=Clouddens(i)*CloudP(i)*exp(-CloudtauUV(i))
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
	if(Cloud(ii)%hazetype.eq.'optEC') Sn=Sn*scaleUV*Sigmadot/tot

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(nnr,NN)
	allocate(vthv(nnr))
	allocate(ScS(nnr),ScN(nnr))
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
				vsed(i)=-rpart(i)*rho_av(i)*Ggrav*Mplanet/(Clouddens(i)*vth(i)*CloudR(i)**2)
				lmfp=2.3d0*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
				vsed(i)=vsed(i)*max(1d0,4d0*rpart(i)/(9d0*lmfp))
			endif
		else
			vsed(i)=-Cloud(ii)%rnuc*rho_nuc*Ggrav*Mplanet/(Clouddens(i)*vth(i)*CloudR(i)**2)
			lmfp=2.3d0*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
			vsed(i)=vsed(i)*max(1d0,4d0*Cloud(ii)%rnuc/(9d0*lmfp))
		endif
		mpart(i)=rho_nuc*4d0*pi*Cloud(ii)%rnuc**3/3d0
	enddo
	if(complexKzz) then
		do i=1,nnr
			St=Cloud(ii)%rnuc*rho_nuc*Km(i)/(vth(i)*Clouddens(i)*CloudHp(i)**2)
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

c==================================================
c============ compute seed particles ==============
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
c loss term to rain particles
		Mb(i)=Mb(i)-ScTot(i)
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

c loss term to rain particles
		Mb(i)=Mb(i)-ScTot(i)

		x(i)=-Sn(i)
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
	xs(1:nnr)=x(1:nnr)

c============ end seed particles ==================
c==================================================

	if(iter.ne.1) then

	vsed=0d0
	do i=1,nnr
		cs=sqrt(kb*CloudT(i)/(2.3d0*mp))
		vth(i)=sqrt(8d0*kb*CloudT(i)/(pi*2.3d0*mp))
		vsed(i)=-rpart(i)*rho_av(i)*Ggrav*Mplanet/(Clouddens(i)*vth(i)*CloudR(i)**2)
		lmfp=2.3d0*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
		vsed(i)=vsed(i)*max(1d0,4d0*rpart(i)/(9d0*lmfp))
		mpart(i)=rho_av(i)*4d0*pi*rpart(i)**3/3d0
	enddo
	if(complexKzz) then
		do i=1,nnr
			St=rpart(i)*rho_av(i)*Km(i)/(vth(i)*Clouddens(i)*CloudHp(i)**2)
			Kd(i)=Km(i)/(1d0+St)
		enddo
	endif

c equations for number of Nuclii
	Ma=0d0
	Mb=0d0
	Mc=0d0
	x=0d0
	i=1
	dztot=(CloudR(i+1)-CloudR(i))
	dz=(CloudR(i)-CloudR(i+1))
	Mc(i)=Mc(i)-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
	Mb(i)=Mb(i)+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
	Mb(i)=Mb(i)+Clouddens(i)*vsed(i)/dztot
	x(i)=-ScTot(i)*xs(i)
c coagulation
	if(Cloud(ii)%coagulation) then
		npart=xn(i)*Clouddens(i)
		lmfp=2.3d0*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
		vmol=0.5d0*lmfp*vth(i)
		Dp=kb*CloudT(i)/(6d0*pi*rpart(i)*vmol*Clouddens(i))
		tcoaginv=sqrt(pi)*3d0*(sum(xc(1:nCS,i))+xm(i))*Ggrav*Mplanet/(8d0*vth(i)*CloudR(i)**2)
c		tcoaginv=tcoaginv*exp(-(vsed(i)/vfrag)**2)
		vBM=sqrt(16d0*kb*CloudT(i)/(pi*mpart(i)))
		if(Dp/rpart(i).lt.vBM) vBM=Dp/rpart(i)
		tcoaginv=tcoaginv+2d0*pi*rpart(i)**2*npart*vBM!*exp(-(vBM/vfrag)**2)
		if(.not.tcoaginv.gt.0d0) tcoaginv=0d0
		tcinv(iter,i)=tcoaginv
		tcoaginv=sum(tcinv(1:iter,i))/real(iter)
		Mb(i)=Mb(i)-Clouddens(i)*tcoaginv
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

		x(i)=-ScTot(i)*xs(i)

c coagulation
		if(Cloud(ii)%coagulation) then
			npart=xn(i)*Clouddens(i)
			lmfp=2.3d0*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
			vmol=0.5d0*lmfp*vth(i)
			Dp=kb*CloudT(i)/(6d0*pi*rpart(i)*vmol*Clouddens(i))

c			tcoaginv=2d0*npart*pi*rpart(i)**2*abs(vsed(i))
c rewritten for better convergence
			tcoaginv=sqrt(pi)*3d0*(sum(xc(1:nCS,i))+xm(i))*Ggrav*Mplanet/(8d0*vth(i)*CloudR(i)**2)
c			tcoaginv=tcoaginv*exp(-(vsed(i)/vfrag)**2)
			vBM=sqrt(16d0*kb*CloudT(i)/(pi*mpart(i)))
			if(Dp/rpart(i).lt.vBM) vBM=Dp/rpart(i)
			tcoaginv=tcoaginv+2d0*pi*rpart(i)**2*npart!*vBM*exp(-(vBM/vfrag)**2)

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
		i=1
		dztot=(CloudR(i+1)-CloudR(i))
		dz=(CloudR(i)-CloudR(i+1))
		Mc(i)=Mc(i)-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		Mb(i)=Mb(i)+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		Mb(i)=Mb(i)+Clouddens(i)*vsed(i)/dztot
		x(i)=-ScTot(i)*xs(i)
		do i=2,nnr-1
			dztot=(CloudR(i+1)-CloudR(i-1))/2d0
			dz=(CloudR(i)-CloudR(i+1))
			Mc(i)=Mc(i)-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
			Mb(i)=Mb(i)+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
			dz=(CloudR(i-1)-CloudR(i))
			Ma(i-1)=Ma(i-1)-0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz/dztot
			Mb(i)=Mb(i)+(0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz)/dztot
			Mb(i)=Mb(i)+Clouddens(i)*vsed(i)/dztot

			x(i)=-ScTot(i)*xs(i)
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

	else
		xn=xs
		xm=xs
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
	xs(1:nnr)=xs(1:nnr)/m_nuc

	endif

	if(Cloud(ii)%condensates) then
	freeflow=Cloud(ii)%freeflow_con

	ScTot=0d0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(cs,iCS,i,j,dz,NRHS,INFO,kl,ku,dztot)
!$OMP& SHARED(nCS,nnr,CloudT,Clouddens,CloudP,mu,fstick,CloudR,densv,Kd,Kg,xn,empty,docondense,
!$OMP&		NN,rpart,ixc,vsed,ixv,m_nuc,mpart,xv_bot,xc,xv,iter,nTiter,i3D,xm,freeflow,ScTot,Cloud,ii,xs)
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC, 1)
	do iCS=1,nCS
	if(docondense(iCS)) then
	do i=1,nnr
		cs=sqrt(kb*CloudT(i)/(2.3d0*mp))
		vthv(i)=sqrt(8d0*kb*CloudT(i)/(pi*mu(iCS)*mp))
		ScS(i)=min(vthv(i)*Cloud(ii)%rnuc,kb*CloudT(i)*sqrt(8d0*kb*CloudT(i)/(pi*2.3*mp))/(3d0*CloudP(i)*1d6*8e-15))
     &				*4d0*pi*Cloud(ii)%rnuc*Clouddens(i)
		ScS(i)=fstick*ScS(i)
		ScN(i)=min(vthv(i)*rpart(i),kb*CloudT(i)*sqrt(8d0*kb*CloudT(i)/(pi*2.3*mp))/(3d0*CloudP(i)*1d6*8e-15))
     &				*4d0*pi*rpart(i)*Clouddens(i)
		ScN(i)=fstick*ScN(i)
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

		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))+(ScS(i)*xs(i)+ScN(i)*xn(i))*Clouddens(i)
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-ScN(i)*densv(i,iCS)/mpart(i)
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

		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))+(ScS(i)*xs(i)+ScN(i)*xn(i))*Clouddens(i)
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-ScN(i)*densv(i,iCS)/mpart(i)
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

		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))-(ScS(i)*xs(i)+ScN(i)*xn(i))*Clouddens(i)
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+ScN(i)*densv(i,iCS)/mpart(i)
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

!$OMP CRITICAL
	do i=1,nnr
		ScTot(i)=ScTot(i)+(ScN(i)*xn(i)+ScS(i)*xs(i))*Clouddens(i)*xv(iCS,i)-ScN(i)*densv(i,iCS)*xc(iCS,i)/mpart(i)
	enddo
!$OMP END CRITICAL

	else

	do i=1,nnr
		xc(iCS,i)=0d0
		xv(iCS,i)=xv_bot(iCS)
	enddo
	
	endif

	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	else
	do iCS=1,nCS
		do i=1,nnr
			xc(iCS,i)=0d0
			xv(iCS,i)=xv_bot(iCS)
		enddo
	enddo
	endif

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
		if(err.gt.maxerr.and..not.empty(i)) maxerr=err
		rpart(i)=sqrt(rr*rpart(i))
	enddo
	do i=1,nnr
		tot=0d0
		do iCS=1,nCS
			tot=tot+xc(iCS,i)
		enddo
		if(tot.gt.0d0) then
			ScTot(i)=max(0d0,ScTot(i)/tot)
		else
			ScTot(i)=0d0
		endif
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

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
	deallocate(vthv)
	deallocate(ScN,ScS)
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

	Cloud(ii)%frac(1:nr,1:Cloud(ii)%nmat)=0d0

	do iCS=1,nCS
		x(1:nnr)=xc(iCS,1:nnr)
		call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,iCS),nr)
	enddo
c Seed particles
	if(Cloud(ii)%haze) then
		x(1:nnr)=xm(1:nnr)
		call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,Cloud(ii)%nmat),nr)
	endif

	if(.not.retrieval) then
		if(do3D) then
			open(unit=20,file=trim(outputdir) // '/cloudstructure' // trim(int2string(i3D,'(i0.4)')) // '.dat',
     &             FORM="FORMATTED",ACCESS="STREAM")
		else
			open(unit=20,file=trim(outputdir) // '/cloudstructure.dat',FORM="FORMATTED",ACCESS="STREAM")
		endif
		form='("#",a18,a19,a19,a19,' // trim(int2string(nCS+1,'(i4)')) // 'a23,a19,a19,a19)'
		write(20,form) "P[bar]","dens[g/cm^3]","xs","xn",(trim(CSname(i)),i=1,nCS),"r[cm]","T[K]","Jstar"
		form='(es19.7E3,es19.7E3,es19.7E3,es19.7E3,' // trim(int2string(nCS,'(i4)')) // 'es23.7E3,es19.7E3,es19.7E3,es19.7E3)'
		do i=1,nnr
			write(20,form) CloudP(i),Clouddens(i),xs(i)*m_nuc,xn(i)*m_nuc,xc(1:nCS,i),rpart(i),CloudT(i),Sn(i)/m_nuc
		enddo
		close(unit=20)
	endif

	do i=1,nr
		tot=sum(Cloud(ii)%frac(i,1:Cloud(ii)%nmat))
		if(tot.gt.0d0) then
			Cloud(ii)%frac(i,1:Cloud(ii)%nmat)=Cloud(ii)%frac(i,1:Cloud(ii)%nmat)/tot
		else
			Cloud(ii)%frac(i,1:Cloud(ii)%nmat)=1d0/real(Cloud(ii)%nmat)
			cloud_dens(i,ii)=0d0
		endif
	enddo


c Elemental abundances
	allocate(at_ab(nr,N_atoms))
	at_ab=0d0
	do iCS=1,nCS
		x(1:nnr)=xv(iCS,1:nnr)
		call regridarray(logCloudP,x,nnr,logP,logx,nr)
		do j=1,N_atoms
			at_ab(1:nr,j)=at_ab(1:nr,j)+logx(1:nr)*mutot*atoms_cloud(iCS,j)/(mu(iCS)*CSnmol(iCS))
		enddo
	enddo

	call cpu_time(time)
	timecloud=timecloud+time
	call system_clock(itime)
	itimecloud=itimecloud+itime

c	open(unit=20,file=trim(outputdir) // '/atoms.dat',FORM="FORMATTED",ACCESS="STREAM")
	if(dochemistry) then
		dochemR=.false.
		dochemR(1)=.true.
		dochemR(nr)=.true.
		do i=1,nr,nrstepchem
			dochemR(i)=.true.
		enddo
		ini=.true.
		do i=1,nr
			call tellertje(i,nr)
			if(dochemR(i)) then
			molfracs_atoms(1:N_atoms)=at_ab(i,1:N_atoms)
			molfracs_atoms=molfracs_atoms+molfracs_atoms0
			molfracs_atoms(3)=molfracs_atoms(3)+COabun
			molfracs_atoms(5)=molfracs_atoms(5)+COabun
			do j=1,N_atoms
				if(.not.molfracs_atoms(j).gt.0d0) then
					molfracs_atoms(j)=0d0
				endif
			enddo
			tot=sum(molfracs_atoms(1:N_atoms))
			molfracs_atoms=molfracs_atoms/tot
			do j=1,N_atoms
				if(.not.molfracs_atoms(j).gt.1d-50) then
					molfracs_atoms(j)=1d-50
				endif
			enddo
			if(nPhotoReacts.gt.0) call doPhotoChemAtom(i)
			if((P(i).ge.mixP.or.i.eq.1).and.dochemistry) then
				call call_chemistry(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,.false.,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol,.false.)
			else
				mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
				XeqCloud(i,1:nclouds)=XeqCloud(i-1,1:nclouds)
				nabla_ad(i)=nabla_ad(i-1)
				MMW(i)=MMW(i-1)
				didcondens(i)=didcondens(i-1)
			endif
c			write(20,*) P(i),molfracs_atoms(1:N_atoms)
			endif
		enddo
		if(nrstepchem.ne.1) then
			do i=1,nmol
				if(includemol(i).or.diseqmol(i)) then
					call fillblanks(P,mixrat_r(1:nr,i),nr,dochemR,.true.)
				endif
			enddo
			call fillblanks(P,MMW,nr,dochemR,.true.)
			call fillblanks(P,nabla_ad,nr,dochemR,.true.)
		endif
		if(fixMMW) MMW=MMW0
c		close(unit=20)
		if(disequilibrium) then
c			call disequilibrium code
c			input: 	R(1:nr+1) : These are the radial boundaries of the layers (bottom to top)
c			P(1:nr),T(1:nr) : These are the pressure and temperature inside the layers
c			molname(1:nmol) : names of the molecules included
c			Kzz_r(1:nr) : Diffusion coefficient
c			input/output:	mixrat_r(1:nr,1:nmol) : number densities inside each layer. Now set to equilibrium abundances.
			call output("==================================================================")
			call output("Computing disequilibrium chemistry")
			call diseq_calc(nr,R(1:nr+1),P(1:nr),T(1:nr),nmol,molname(1:nmol),mixrat_r(1:nr, 1:nmol),COratio,Kzz_g(1:nr))
		endif
	else
		do iCS=1,nCS
			if(CSname(iCS).eq."H2O") then
				x(1:nnr)=xv(iCS,1:nnr)
				call regridarray(logCloudP,x,nnr,logP,logx,nr)
				mixrat_r(2:nr,1)=logx(2:nr)*mutot/(mu(iCS)*CSnmol(iCS))
			endif
		enddo
		do i=1,nr
			tot=0d0
			do j=1,nmol
				if(mixrat_r(i,j).gt.0d0) tot=tot+mixrat_r(i,j)
			enddo
			if(tot.gt.0d0) mixrat_r(i,1:nmol)=mixrat_r(i,1:nmol)/tot
		enddo
	endif
	deallocate(at_ab)

	do i=1,nr
		do j=1,nmol
			if(.not.mixrat_r(i,j).gt.0d0) mixrat_r(i,j)=0d0
		enddo
	enddo
	if(nfixmol.gt.0) then
		do i=1,nfixmol
			mixrat_r(1:nr,ifixmol(i))=fixmol_abun(i)*exp(-P(1:nr)/fixmol_P(i))
		enddo
	endif
	call doPhotoChemMol()
	do j=1,nmol
		if(isotope(j).gt.0) then
			do i=1,nr
				mixrat_r(i,isotope(j))=mixrat_r(i,j)/f_isotope(j)
				mixrat_r(i,j)=mixrat_r(i,j)*(1d0-1d0/f_isotope(j))
			enddo
		endif
	enddo
	
	if(.not.retrieval) then
		call ComputeTevap()
		if(complexKzz) then
			open(unit=50,file=trim(outputdir) // 'cloudKzz.dat',FORM="FORMATTED",ACCESS="STREAM")
			form='("#",a12,a13,a13,a13)'
			write(50,trim(form)) "Kzz [cm^2/s]","P [bar]","Kpart","Kgas"
			form='(es13.3E3,es13.3E3,es13.3E3,es13.3E3)'
			do i=1,nnr
				write(50,trim(form)) Km(i),CloudP(i),Kd(i),Kg(i)
			enddo
			close(unit=50)
		endif
	endif


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
	
	


	subroutine PvapNH3(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0,c1,c2
	parameter(c0=10.53,c1=-2161.0,c2=-86596.0)

	Pvap=exp(c0+c1/T+c2/T**2)
	liquid=.true.

	return
	end
	
	

	subroutine PvapNaCl(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0_i,c1_i,c2_i,c3_i,c4_i
	real*8 c0_l,c1_l,c2_l,c3_l,c4_l
	parameter(	
     &	c0_i=-2.79146E+04,
     &	c1_i=3.46023E+01,
     &	c2_i=-3.11287E-03,
     &	c3_i=5.30965E-07,
     &	c4_i=-2.59584E-12,
     &	c0_l=-2.48880E+04,
     &	c1_l=3.18494E+01,
     &	c2_l=-3.08748E-03,
     &	c3_l=4.84990E-07,
     &	c4_l=-2.60359E-11)

	Pi=c0_i/T+c1_i+c2_i*T+c3_i*T**2+c4_i*T**3
	Pl=c0_l/T+c1_l+c2_l*T+c3_l*T**2+c4_l*T**3

	if(Pi.lt.Pl) then
		Pvap=Pi
		liquid=.false.
	else
		Pvap=Pl
		liquid=.true.
	endif
	Pvap=1d-6*exp(Pvap)
	
	return
	end
		
	
	subroutine PvapKCl(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0_i,c1_i,c2_i,c3_i,c4_i
	real*8 c0_l,c1_l,c2_l,c3_l,c4_l
	parameter(	
     &	c0_i=-2.69250E+04,
     &	c1_i=3.39574E+01,
     &	c2_i=-2.04903E-03,
     &	c3_i=-2.83957E-07,
     &	c4_i=1.82974E-10,
     &	c0_l=-2.50293E+04,
     &	c1_l=3.39453E+01,
     &	c2_l=-4.61815E-03,
     &	c3_l=7.36857E-07,
     &	c4_l=0.00000E+00)

	Pi=c0_i/T+c1_i+c2_i*T+c3_i*T**2+c4_i*T**3
	Pl=c0_l/T+c1_l+c2_l*T+c3_l*T**2+c4_l*T**3

	if(Pi.lt.Pl) then
		Pvap=Pi
		liquid=.false.
	else
		Pvap=Pl
		liquid=.true.
	endif
	Pvap=1d-6*exp(Pvap)
		
	return
	end
	