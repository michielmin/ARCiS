	subroutine CondensationCloud(ii)
	use GlobalSetup
	use Constants
	use AtomsModule
	use CloudModule
	use TimingModule
	IMPLICIT NONE
	real*8,allocatable :: x(:),vsed(:),dx(:),vth(:)
	real*8,allocatable :: Sn(:),mpart(:)
	real*8,allocatable :: y(:,:),CloudHp(:),CloudMMW(:)
	real*8,allocatable :: at_ab(:,:)
	real*8,allocatable,save :: Sc(:,:),vthv(:),IWORKomp(:),AB(:,:)
	real*8,allocatable :: tcinv(:,:),rho_av(:),Kd(:),Kg(:),Km(:),tcrystinv(:)
	real*8 dz,g,rr,mutot,npart,tot,lambda,tot1,tot2,tot3,nucryst,Tcryst
	integer info,i,j,iter,NN,NRHS,niter,ii,k,ihaze,kl,ku,jCS
	real*8 cs,eps,frac_nuc,m_nuc,tcoaginv,Dp,vmol,f,mm,err,maxerr,dztot
	real*8 Pv,molfracs_atoms0(N_atoms),NKn,Kzz_r(nr),vBM,scale
	integer,allocatable :: IWORK(:),ixv(:,:),ixc(:,:),iVL(:,:),ixn(:),icryst(:),iamorph(:),inuc(:)
	real*8 sigmastar,Sigmadot,Pstar,sigmamol,COabun,lmfp,fstick,kappa_cloud,fmin,rho_nuc,Gibbs
	logical ini,Tconverged,haze_nuc
	character*500 cloudspecies(max(nclouds,1)),form
	real*8,allocatable :: xn_iter(:,:),xc_iter(:,:,:),xv_iter(:,:,:)
	logical,allocatable :: docondense(:)
	integer iCS,ir,nrdo,iconv,nconv,iVS,nVS,NStot,ik
	real*8 logP(nr),logx(nr),dlogx(nr),St,fsed
	real*8,allocatable :: logCloudP(:),CloudtauUV(:),CloudkappaUV(:),CloudG(:)
	character*10,allocatable :: v_names(:),v_names_out(:)
	logical,allocatable :: v_include(:),c_rainout(:),do_nuc(:)
	integer INCFD,IERR
	logical SKIP,liq
	real*8 time,kp,Otot(nr),Ctot(nr),Ntot(nr)
	integer itime
	real*8,allocatable :: v_atoms(:,:),muC(:),muV(:),v_cloud(:,:),Sat(:,:),Sat0(:,:),fSat(:,:),v_H2(:)
	real*8,allocatable :: xv_out(:),Jn_xv(:,:),A_J(:),B_j(:),sigma_nuc(:),r0_nuc(:),Nf_nuc(:),Nc_nuc(:,:),Jn_out(:)

	logical dochemR(nr)

	call cpu_time(time)
	timecloud=timecloud-time
	call system_clock(itime)
	itimecloud=itimecloud-itime
	ctimecloud=ctimecloud+1

	nVS=15
	allocate(v_names(nVS),v_atoms(nVS,N_atoms),v_include(nVS))
	
	v_atoms=0d0

	v_names(1)="SiO"
	v_atoms(1,9)=1
	v_atoms(1,5)=1

	v_names(2)="TiO"
	v_atoms(2,5)=1
	v_atoms(2,15)=1

	v_names(3)="Mg"
	v_atoms(3,7)=1

	v_names(4)="H2O"
	v_atoms(4,1)=2
	v_atoms(4,5)=1

	v_names(5)="H2S"
	v_atoms(5,1)=2
	v_atoms(5,11)=1

	v_names(6)="Fe"
	v_atoms(6,17)=1

	v_names(7)="Al"
	v_atoms(7,8)=1

	v_names(8)="Na"
	v_atoms(8,6)=1

	v_names(9)="K"
	v_atoms(9,13)=1

	v_names(10)="HCl"
	v_atoms(10,1)=1
	v_atoms(10,12)=1

	v_names(11)="NH3"
	v_atoms(11,1)=3
	v_atoms(11,4)=1

	v_names(12)="Zn"
	v_atoms(12,30)=1

	v_names(13)="Mn"
	v_atoms(13,27)=1

	v_names(14)="Cr"
	v_atoms(14,26)=1

	v_names(15)="CH4"
	v_atoms(15,1)=4
	v_atoms(15,3)=1

	v_include=.false.

	nnr=(nr-1)*nr_cloud+1
	nCS=Cloud(ii)%nmat
	if(Cloud(ii)%haze) nCS=nCS-1
	allocate(Kd(nnr),Kg(nnr),Km(nnr))
	allocate(logCloudP(nnr))
	allocate(CloudMMW(nnr),CloudHp(nnr),Sat(nnr,nCS),Sat0(nnr,nCS),fSat(nnr,nCS),CloudG(nnr))
	
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
	
	allocate(docondense(nCS))

	if(.not.allocated(CloudP)) then
		allocate(CloudP(nnr))
		allocate(CloudT(nnr))
		allocate(CloudR(nnr))
		allocate(Clouddens(nnr))
		allocate(xv(nVS,nnr))
		allocate(xc(nCS,nnr))
		allocate(xn(nnr))
		allocate(xnv(nnr))
		allocate(rpart(nnr))
	endif
	if(.not.allocated(ATP)) then
		allocate(ATP(nCS),BTP(nCS))
		allocate(rhodust(nCS))
		allocate(atoms_cloud(nCS,N_atoms))
		allocate(xv_bot(nVS))
		allocate(mu(nCS))
		allocate(CSname(nCS),maxT(nCS))
	endif
	allocate(muC(nCS))
	allocate(muV(nVS))
	allocate(v_cloud(nCS,nVS),iVL(nnr,nCS),v_H2(nCS))
	allocate(icryst(nCS),iamorph(nCS),tcrystinv(nnr))
	allocate(A_J(nCS),B_J(nCS),do_nuc(nCS),Jn_xv(nnr,nCS),sigma_nuc(nCS),r0_nuc(nCS),Nf_nuc(nCS),Nc_nuc(nnr,nCS))
	allocate(inuc(nCS))

	call SetAbun

	COabun=min(molfracs_atoms(3),molfracs_atoms(5))
	molfracs_atoms(3)=molfracs_atoms(3)-COabun
	molfracs_atoms(5)=molfracs_atoms(5)-COabun

	atoms_cloud=0
	v_cloud=0d0
	v_H2=0d0
	icryst=0
	iamorph=0
	do_nuc=.false.
	do i=1,nCS
		select case(Cloud(ii)%condensate(i))
			case('SiO2','QUARTZ')
				CSname(i)='SiO2'
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=2
				v_cloud(i,1)=1
				v_cloud(i,4)=1
				v_H2(i)=-1
				v_include(1)=.true.
				v_include(4)=.true.
				rhodust(i)=2.65
			case('MgSiO3','ENSTATITE')
				CSname(i)=trim(Cloud(ii)%condensate(i))
				atoms_cloud(i,7)=1
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=3
				v_cloud(i,1)=1
				v_cloud(i,3)=1
				v_cloud(i,4)=2
				v_H2(i)=-2
				v_include(1)=.true.
				v_include(3)=.true.
				v_include(4)=.true.
				rhodust(i)=3.19
			case('Mg2SiO4','FORSTERITE')
				CSname(i)=trim(Cloud(ii)%condensate(i))
				atoms_cloud(i,7)=2
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=4
				v_cloud(i,1)=1
				v_cloud(i,3)=2
				v_cloud(i,4)=3
				v_H2(i)=-3
				v_include(1)=.true.
				v_include(3)=.true.
				v_include(4)=.true.
				rhodust(i)=3.21
			case('MgO')
				CSname(i)='MgO'
				atoms_cloud(i,7)=1
				atoms_cloud(i,5)=1
				v_cloud(i,3)=1
				v_cloud(i,4)=1
				v_H2(i)=-1
				v_include(3)=.true.
				v_include(4)=.true.
				rhodust(i)=3.58
			case('H2O','WATER')
				CSname(i)='H2O'
				atoms_cloud(i,1)=2
				atoms_cloud(i,5)=1
				v_cloud(i,4)=1
				v_include(4)=.true.
				rhodust(i)=0.93
				maxT(i)=747d0
				do_nuc(i)=Cloud(ii)%ComputeJn
				inuc(i)=4
				sigma_nuc(i)=109.0
				r0_nuc(i)=1.973e-8
				Nf_nuc(i)=1d0
			case('Fe','IRON')
				CSname(i)='Fe'
				atoms_cloud(i,17)=1
				v_cloud(i,6)=1
				v_include(6)=.true.
				rhodust(i)=7.87
			case('FeS','TROILITE')
				CSname(i)='FeS'
				atoms_cloud(i,17)=1
				atoms_cloud(i,11)=1
				v_cloud(i,5)=1
				v_cloud(i,6)=1
				v_H2(i)=-1
				v_include(5)=.true.
				v_include(6)=.true.
				maxT(i)=680d0
				rhodust(i)=4.83
			case('FeO')
				CSname(i)='FeO'
				atoms_cloud(i,17)=1
				atoms_cloud(i,5)=1
				v_cloud(i,6)=1
				v_cloud(i,4)=1
				v_H2(i)=-1
				v_include(6)=.true.
				v_include(4)=.true.
				rhodust(i)=5.9
			case('Fe2O3')
				CSname(i)='Fe2O3'
				atoms_cloud(i,17)=2
				atoms_cloud(i,5)=3
				v_cloud(i,6)=2
				v_cloud(i,4)=3
				v_H2(i)=-3
				v_include(6)=.true.
				v_include(4)=.true.
				rhodust(i)=5.24
			case('Al2O3','CORRUNDUM')
				CSname(i)='Al2O3'
				atoms_cloud(i,5)=3
				atoms_cloud(i,8)=2
				v_cloud(i,4)=3
				v_cloud(i,7)=2
				v_H2(i)=-3
				v_include(4)=.true.
				v_include(7)=.true.
				rhodust(i)=3.97
			case("NaCl")
				CSname(i)='NaCl'
				atoms_cloud(i,6)=1
				atoms_cloud(i,12)=1
				v_cloud(i,8)=1
				v_cloud(i,10)=1
				v_H2(i)=-0.5
				v_include(8)=.true.
				v_include(10)=.true.
				rhodust(i)=2.17
				do_nuc(i)=Cloud(ii)%ComputeJn
				inuc(i)=10
				sigma_nuc(i)=113.3
				r0_nuc(i)=2.205e-8
				Nf_nuc(i)=1d0
			case("KCl")
				CSname(i)='KCl'
				atoms_cloud(i,13)=1
				atoms_cloud(i,12)=1
				v_cloud(i,9)=1
				v_cloud(i,10)=1
				v_H2(i)=-0.5
				v_include(9)=.true.
				v_include(10)=.true.
				rhodust(i)=1.99
				do_nuc(i)=Cloud(ii)%ComputeJn
				inuc(i)=10
				sigma_nuc(i)=100.3
				r0_nuc(i)=2.462e-8
				Nf_nuc(i)=1d0
			case("Na2S")
				CSname(i)='Na2S'
				atoms_cloud(i,6)=1
				atoms_cloud(i,11)=2
				v_cloud(i,8)=1
				v_cloud(i,5)=1
				v_H2(i)=-1
				v_include(8)=.true.
				v_include(5)=.true.
				rhodust(i)=1.86
			case("NH3","AMONIA")
				CSname(i)='NH3'
				atoms_cloud(i,1)=3
				atoms_cloud(i,4)=1
				v_cloud(i,11)=1
				v_include(11)=.true.
				rhodust(i)=0.87
				maxT(i)=220d0
				do_nuc(i)=Cloud(ii)%ComputeJn
				inuc(i)=11
				sigma_nuc(i)=23.4
				r0_nuc(i)=1.980e-8
				Nf_nuc(i)=1d0
			case('TiO2')
				CSname(i)='TiO2'
				atoms_cloud(i,5)=1
				atoms_cloud(i,15)=2
				v_cloud(i,2)=1
				v_cloud(i,4)=1
				v_H2(i)=-1
				v_include(2)=.true.
				v_include(4)=.true.
				rhodust(i)=4.23
				do_nuc(i)=Cloud(ii)%ComputeJn
				A_J(i)=1.112e12
				B_J(i)=-24.0
				inuc(i)=2
				sigma_nuc(i)=480.6
				r0_nuc(i)=1.956e-8
				Nf_nuc(i)=0d0
			case('H2SO4')
				CSname(i)='H2SO4'
				atoms_cloud(i,1)=2
				atoms_cloud(i,5)=4
				atoms_cloud(i,11)=1
				v_cloud(i,4)=4
				v_cloud(i,5)=1
				v_H2(i)=-4
				v_include(4)=.true.
				v_include(5)=.true.
				rhodust(i)=1.84d0
			case('ZnS')
				CSname(i)='ZnS'
				atoms_cloud(i,30)=1
				atoms_cloud(i,11)=1
				v_cloud(i,12)=1
				v_cloud(i,5)=1
				v_H2(i)=-1
				v_include(12)=.true.
				v_include(5)=.true.
				rhodust(i)=4.09
			case('MnS')
				CSname(i)='MnS'
				atoms_cloud(i,27)=1
				atoms_cloud(i,11)=1
				v_cloud(i,13)=1
				v_cloud(i,5)=1
				v_H2(i)=-1
				v_include(13)=.true.
				v_include(5)=.true.
				rhodust(i)=4.08
			case('Zn')
				CSname(i)='Zn'
				atoms_cloud(i,30)=1
				v_cloud(i,12)=1
				v_include(12)=.true.
				rhodust(i)=7.14d0
			case('Mn')
				CSname(i)='Mn'
				atoms_cloud(i,27)=1
				v_cloud(i,13)=1
				v_include(13)=.true.
				rhodust(i)=7.43d0
			case('Cr')
				CSname(i)='Cr'
				atoms_cloud(i,26)=1
				v_cloud(i,14)=1
				v_include(14)=.true.
				rhodust(i)=7.19d0
				do_nuc(i)=Cloud(ii)%ComputeJn
				inuc(i)=14
				sigma_nuc(i)=3330.0
				r0_nuc(i)=1.421
				Nf_nuc(i)=1d0
			case('NH4Cl')
				CSname(i)='NH4Cl'
				atoms_cloud(i,1)=4
				atoms_cloud(i,4)=1
				atoms_cloud(i,12)=1
				v_cloud(i,10)=1
				v_cloud(i,11)=1
				v_include(10)=.true.
				v_include(11)=.true.
				rhodust(i)=1.53
			case('SiO')
				CSname(i)='SiO'
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=1
				v_cloud(i,1)=1
				v_include(1)=.true.
				rhodust(i)=2.18
				maxT(i)=5000d0
				do_nuc(i)=Cloud(ii)%computeJn
				A_J(i)=4.4e12
				B_J(i)=1.33
				inuc(i)=1
				sigma_nuc(i)=849.4
				r0_nuc(i)=2.0e-8
				Nf_nuc(i)=1d0
			case default
				call output("Unknown condensate")
				stop
		end select
	enddo
	
	do iCS=1,nCS
		do jCS=1,nCS
			if(Cloud(ii)%condensate(iCS).eq.'FORSTERITE'.and.Cloud(ii)%condensate(jCS).eq.'Mg2SiO4') then
				iamorph(iCS)=jCS
				icryst(jCS)=iCS
			endif
			if(Cloud(ii)%condensate(iCS).eq.'ENSTATITE'.and.Cloud(ii)%condensate(jCS).eq.'MgSiO3') then
				iamorph(iCS)=jCS
				icryst(jCS)=iCS
			endif
		enddo
	enddo

	docondense=.true.

	haze_nuc=.false.
	select case(Cloud(ii)%hazetype)
		case("SOOT","soot","Soot")
			rho_nuc=1.00
		case("THOLIN","tholin","Tholin")
			rho_nuc=1.00
		case("optEC")
			rho_nuc=1.50
			ihaze=15
			v_include(15)=(.not.dochemistry)
			haze_nuc=(.not.dochemistry)
		case("SiC")
			rho_nuc=3.22
		case("CARBON","Carbon","carbon")
			rho_nuc=1.80
		case("CORRUNDUM","Corrundum","corrundum","Al2O3")
			rho_nuc=3.97
		case("IRON","Iron","iron","Fe")
			rho_nuc=7.87
		case("SiO")
			rho_nuc=2.18
		case("TiO2")
			rho_nuc=4.23
		case("Enstatite","enstatite","ENSTATITE")
			rho_nuc=3.20
		case default
			call output("hazetype unknown")
			stop
	end select

	do i=1,nCS
		muC(i)=sum(mass_atoms(1:N_atoms)*atoms_cloud(i,1:N_atoms))
	enddo
	do i=1,nVS
		muV(i)=sum(mass_atoms(1:N_atoms)*v_atoms(i,1:N_atoms))
	enddo

	if(dochemistry) then
		mutot=COabun*(mass_atoms(3)+mass_atoms(5))
		xv_bot=1d200
		do i=1,N_atoms
			mutot=mutot+mass_atoms(i)*molfracs_atoms(i)
		enddo
	else
		mutot=0d0
		do i=1,nmol
			if(includemol(i)) then
				mutot=mutot+Mmol(i)*mixrat_r(1,i)
			endif
		enddo
	endif
	do iVS=1,nVS
		if(v_include(iVS)) then
			if(dochemistry) then
				do k=1,N_atoms
					if(v_atoms(iVS,k).gt.0d0) then
						f=molfracs_atoms(k)/v_atoms(iVS,k)
						if(f.lt.xv_bot(iVS)) then
							xv_bot(iVS)=f
						endif
					endif
				enddo
				do k=1,N_atoms
					molfracs_atoms(k)=molfracs_atoms(k)-xv_bot(iVS)*v_atoms(iVS,k)
					if(molfracs_atoms(k).lt.0d0) molfracs_atoms(k)=0d0
				enddo
			else
				do k=1,nmol
					if(molname(k).eq.v_names(iVS)) then
						xv_bot(iVS)=mixrat_r(1,k)
					endif
				enddo
			endif
		endif
	enddo
	
	molfracs_atoms0=molfracs_atoms
	xv_bot=xv_bot*muV/mutot
	do iVS=1,nVS
		if(.not.v_include(iVS)) xv_bot(iVS)=0d0
	enddo

c	print*,xv_bot(1:7)
c	xv_bot(1:7) = 10.0**metallicity*(/ 6.1e-4, 2.4e-6, 4.1e-4, 1.4e-3, 1.9e-4, 7.6e-4, 3.2e-5 /)
c	print*,xv_bot(1:7)

	do iCS=1,nCS
		if(do_nuc(iCS).and.CSname(iCS).ne.'SiO'.and.CSname(iCS).eq.'TiO2') then
			A_J(iCS)=(16.*pi*sigma_nuc(iCS)**3*(muC(iCS)*mp/rhodust(iCS))**2/(3.*kb**3))
			B_J(iCS)=log(sqrt(2.*sigma_nuc(iCS)/(pi*muC(iCS)*mp))*(muC(iCS)*mp/rhodust(iCS)))
		endif
	enddo

	Cloud(ii)%frac=0d0

	sigmastar=0.1
	Pstar=60d-6

	Pstar=Cloud(ii)%P
	sigmastar=log(Cloud(ii)%dP)

	fstick=1d0
	
	sigmamol=8d-15

	eps=1d-3

	Sigmadot=Cloud(ii)%Sigmadot
	

	m_nuc=4d0*pi*Cloud(ii)%rnuc**3*rho_nuc/3d0

	allocate(mpart(nnr))
	allocate(rho_av(nnr))
	allocate(y(nnr,5))
	allocate(Sn(nnr))
	allocate(vth(nnr))
	allocate(tcinv(niter,nnr))
	allocate(xn_iter(niter,nnr))
	allocate(xc_iter(niter,nCS,nnr),xv_iter(niter,nVS,nnr))
	allocate(vsed(nnr))

	allocate(ixv(nVS,nnr))
	allocate(ixc(nCS,nnr))
	allocate(ixn(nnr))

	j=0
	do i=1,nnr
		if(.not.Cloud(ii)%usefsed) then
			j=j+1
			ixn(i)=j
		endif
		do iCS=1,nCS
			j=j+1
			ixc(iCS,i)=j
		enddo
		do iVS=1,nVS
			if(v_include(iVS)) then
				j=j+1
				ixv(iVS,i)=j
			endif
		enddo
		if(i.eq.1) NStot=j
	enddo
	NN=j

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
	logx(1:nr)=MMW(1:nr)
	call DPCHIM(nr,logP,logx,dlogx,INCFD)
	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudMMW,IERR)

	if(complexKzz.or.Cloud(ii)%usefsed) then
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

	allocate(x(NN))
	allocate(IWORK(10*NN*NN))

c	if(.not.computeT.or.nTiter.eq.1) then
		rpart=Cloud(ii)%rnuc
		xn=0d0
		xv=0d0
		xc=0d0
		do j=1,nVS
			do i=1,nnr
				xv(j,i)=xv_bot(j)
			enddo
		enddo
c	endif
	tcinv=0d0
	Jn_xv=0d0
	Nc_nuc=0d0
	
	docondense=.true.
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

	if(.not.retrieval) then
		open(unit=20,file='Tevap.dat',RECL=6000)
		form='("#",a18,' // trim(int2string(nCS,'(i4)')) // 'a23)'
		write(20,form) "T[K]",(trim(CSname(i)),i=1,nCS)
		form='(es19.7E3,' // trim(int2string(nCS,'(i4)')) // 'es23.7E3)'
		do i=10,5000,10
			tot1=real(i)
			do iCS=1,nCS
				select case(CSname(iCS))
					case("H2O")
						call PvapH2O(tot1,Sat(1,iCS),liq)
						Sat(1,iCS)=1d0/Sat(1,iCS)
					case("NH3")
c	Gibbs energy as derived from Eq from GGChem paper does not work at high pressures. Use adjusted Pvap.
						call PvapNH3(tot1,Sat(1,iCS),liq)
						Sat(1,iCS)=1d0/Sat(1,iCS)
					case("SiO")
c	Gibbs energy as derived from Eq from GGChem paper does not work at high pressures. Use adjusted Pvap.
						call PvapSiO(tot1,Sat(1,iCS),liq)
						Sat(1,iCS)=1d0/Sat(1,iCS)
					case default
						select case(CSname(iCS))
							case("SiO2")
								call Gibbs_SiO2_s(tot1,Gibbs,maxT(iCS))
							case("MgSiO3","ENSTATITE")
								call Gibbs_MgSiO3_s(tot1,Gibbs,maxT(iCS))
							case("Mg2SiO4","FORSTERITE")
								call Gibbs_Mg2SiO4_s(tot1,Gibbs,maxT(iCS))
							case("MgO")
								call Gibbs_MgO_s(tot1,Gibbs,maxT(iCS))
							case("H2O")
								call Gibbs_H2O_s(tot1,Gibbs,maxT(iCS))
							case("Fe")
								call Gibbs_Fe_s(tot1,Gibbs,maxT(iCS))
							case("FeO")
								call Gibbs_FeO_s(tot1,Gibbs,maxT(iCS))
							case("Fe2O3")
								call Gibbs_Fe2O3_s(tot1,Gibbs,maxT(iCS))
							case("FeS")
								call Gibbs_FeS_s(tot1,Gibbs,maxT(iCS))
							case("Al2O3")
								call Gibbs_Al2O3_s(tot1,Gibbs,maxT(iCS))
							case("NaCl")
								call Gibbs_NaCl_s(tot1,Gibbs,maxT(iCS))
							case("KCl")
								call Gibbs_KCl_s(tot1,Gibbs,maxT(iCS))
							case("Na2S")
								call Gibbs_Na2S_s(tot1,Gibbs,maxT(iCS))
							case("TiO2")
								call Gibbs_TiO2_s(tot1,Gibbs,maxT(iCS))
							case("H2SO4")
								call Gibbs_H2SO4_s(tot1,Gibbs,maxT(iCS))
							case("NH3")
								call Gibbs_NH3_s(tot1,Gibbs,maxT(iCS))
							case("ZnS")
								call Gibbs_ZnS_s(tot1,Gibbs,maxT(iCS))
							case("Zn")
								call Gibbs_Zn_s(tot1,Gibbs,maxT(iCS))
							case("MnS")
								Gibbs=4.184*(1.12482E+05/tot1-1.81938E+05
     &			+5.87107E+01*tot1+8.89360E-05*tot1**2-4.20876E-09*tot1**3)
								maxT(iCS)=5000d0
							case("Mn")
								call Gibbs_Mn_s(tot1,Gibbs,maxT(iCS))
							case("Cr")
								call Gibbs_Cr_s(tot1,Gibbs,maxT(iCS))
							case default
								print*,'Unknown condensate'
								stop
						end select
						Sat(1,iCS)=Gibbs
						do iVS=1,nVS
							if(v_cloud(iCS,iVS).gt.0d0.and.v_include(iVS)) then
								select case(v_names(iVS))
									case("SiO")
										call Gibbs_SiO_g(tot1,Gibbs)
									case("H2O")
										call Gibbs_H2O_g(tot1,Gibbs)
									case("Mg")
										call Gibbs_Mg_g(tot1,Gibbs)
									case("Fe")
										call Gibbs_Fe_g(tot1,Gibbs)
									case("H2S")
										call Gibbs_H2S_g(tot1,Gibbs)
									case("Al")
										call Gibbs_Al_g(tot1,Gibbs)
									case("K")
										call Gibbs_K_g(tot1,Gibbs)
									case("Na")
										call Gibbs_Na_g(tot1,Gibbs)
									case("HCl")
										call Gibbs_HCl_g(tot1,Gibbs)
									case("NH3")
										call Gibbs_NH3_g(tot1,Gibbs)
									case("TiO")
										call Gibbs_TiO_g(tot1,Gibbs)
									case("Zn")
										call Gibbs_Zn_g(tot1,Gibbs)
									case("Mn")
										call Gibbs_Mn_g(tot1,Gibbs)
									case("Cr")
										call Gibbs_Cr_g(tot1,Gibbs)
									case default
										print*,'Unknown gas phase'
										stop
								end select
								Sat(1,iCS)=Sat(1,iCS)-Gibbs*v_cloud(iCS,iVS)
							endif
						enddo
						Sat(1,iCS)=exp(-Sat(1,iCS)/(8.314e-3*tot1))
				end select
				if(tot1.gt.maxT(iCS)) Sat(1,iCS)=0d0
				if(Sat(1,iCS).lt.1d-40) Sat(1,iCS)=1d-40
				do iVS=1,nVS
					if(v_include(iVS)) then
						Sat(1,iCS)=Sat(1,iCS)*(xv_bot(iVS)*CloudMMW(1)/muV(iVS))**v_cloud(iCS,iVS)
					endif
				enddo
				Sat(1,iCS)=(1d0/Sat(1,iCS))**(1d0/(v_H2(iCS)+sum(v_cloud(iCS,1:nVS))))			
			enddo
			write(20,form) tot1,Sat(1,1:nCS)
		enddo
		close(unit=20)
	endif


c Compute crystallinity
c values from Fabian et al. 2000
	nucryst=2d13
	Tcryst=41000d0

	do i=1,nnr
		do iCS=1,nCS
			select case(CSname(iCS))
				case("H2O")
					call PvapH2O(CloudT(i),Sat(i,iCS),liq)
					Sat(i,iCS)=CloudP(i)/Sat(i,iCS)
				case("NH3")
c	Gibbs energy as derived from Eq from GGChem paper does not work at high pressures. Use adjusted Pvap.
					call PvapNH3(CloudT(i),Sat(i,iCS),liq)
					Sat(i,iCS)=CloudP(i)/Sat(i,iCS)
				case("SiO")
c	Gibbs energy as derived from Eq from GGChem paper does not work at high pressures. Use adjusted Pvap.
					call PvapSiO(CloudT(i),Sat(i,iCS),liq)
					Sat(i,iCS)=CloudP(i)/Sat(i,iCS)
				case default
					select case(CSname(iCS))
						case("SiO2")
							call Gibbs_SiO2_s(CloudT(i),Gibbs,maxT(iCS))
						case("MgSiO3","ENSTATITE")
							call Gibbs_MgSiO3_s(CloudT(i),Gibbs,maxT(iCS))
						case("Mg2SiO4","FORSTERITE")
							call Gibbs_Mg2SiO4_s(CloudT(i),Gibbs,maxT(iCS))
						case("MgO")
							call Gibbs_MgO_s(CloudT(i),Gibbs,maxT(iCS))
						case("H2O")
							call Gibbs_H2O_s(CloudT(i),Gibbs,maxT(iCS))
						case("Fe")
							call Gibbs_Fe_s(CloudT(i),Gibbs,maxT(iCS))
						case("FeO")
							call Gibbs_FeO_s(CloudT(i),Gibbs,maxT(iCS))
						case("Fe2O3")
							call Gibbs_Fe2O3_s(CloudT(i),Gibbs,maxT(iCS))
						case("FeS")
							call Gibbs_FeS_s(CloudT(i),Gibbs,maxT(iCS))
						case("Al2O3")
							call Gibbs_Al2O3_s(CloudT(i),Gibbs,maxT(iCS))
						case("NaCl")
							call Gibbs_NaCl_s(CloudT(i),Gibbs,maxT(iCS))
						case("KCl")
							call Gibbs_KCl_s(CloudT(i),Gibbs,maxT(iCS))
						case("Na2S")
							call Gibbs_Na2S_s(CloudT(i),Gibbs,maxT(iCS))
						case("TiO2")
							call Gibbs_TiO2_s(CloudT(i),Gibbs,maxT(iCS))
						case("H2SO4")
							call Gibbs_H2SO4_s(CloudT(i),Gibbs,maxT(iCS))
						case("NH3")
							call Gibbs_NH3_s(CloudT(i),Gibbs,maxT(iCS))
						case("ZnS")
							call Gibbs_ZnS_s(CloudT(i),Gibbs,maxT(iCS))
						case("Zn")
							call Gibbs_Zn_s(CloudT(i),Gibbs,maxT(iCS))
						case("MnS")
							Gibbs=4.184*(1.12482E+05/CloudT(i)-1.81938E+05
     &			+5.87107E+01*CloudT(i)+8.89360E-05*CloudT(i)**2-4.20876E-09*CloudT(i)**3)
							maxT(iCS)=5000d0
						case("Mn")
							call Gibbs_Mn_s(CloudT(i),Gibbs,maxT(iCS))
						case("Cr")
							call Gibbs_Cr_s(CloudT(i),Gibbs,maxT(iCS))
						case default
							print*,'Unknown condensate'
							stop
					end select
					Sat(i,iCS)=Gibbs
					do iVS=1,nVS
						if(v_cloud(iCS,iVS).gt.0d0.and.v_include(iVS)) then
							select case(v_names(iVS))
								case("SiO")
									call Gibbs_SiO_g(CloudT(i),Gibbs)
								case("H2O")
									call Gibbs_H2O_g(CloudT(i),Gibbs)
								case("Mg")
									call Gibbs_Mg_g(CloudT(i),Gibbs)
								case("Fe")
									call Gibbs_Fe_g(CloudT(i),Gibbs)
								case("H2S")
									call Gibbs_H2S_g(CloudT(i),Gibbs)
								case("Al")
									call Gibbs_Al_g(CloudT(i),Gibbs)
								case("K")
									call Gibbs_K_g(CloudT(i),Gibbs)
								case("Na")
									call Gibbs_Na_g(CloudT(i),Gibbs)
								case("HCl")
									call Gibbs_HCl_g(CloudT(i),Gibbs)
								case("NH3")
									call Gibbs_NH3_g(CloudT(i),Gibbs)
								case("TiO")
									call Gibbs_TiO_g(CloudT(i),Gibbs)
								case("Zn")
									call Gibbs_Zn_g(CloudT(i),Gibbs)
								case("Mn")
									call Gibbs_Mn_g(CloudT(i),Gibbs)
								case("Cr")
									call Gibbs_Cr_g(CloudT(i),Gibbs)
								case default
									print*,'Unknown gas phase'
									stop
							end select
							Sat(i,iCS)=Sat(i,iCS)-Gibbs*v_cloud(iCS,iVS)
						endif
					enddo
					Sat(i,iCS)=exp(-Sat(i,iCS)/(8.314e-3*CloudT(i)))*CloudP(i)
			end select
			if(CloudT(i).gt.maxT(iCS)) Sat(i,iCS)=0d0
			if(Sat(i,iCS).lt.1d-40) Sat(i,iCS)=1d-40
		enddo
		CloudG(i)=Ggrav*Mplanet/CloudR(i)**2
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
			Sn(i)=(Clouddens(i)*CloudG(i)*Sigmadot/(sigmastar*CloudP(i)*1d6*sqrt(2d0*pi)))*exp(-log(CloudP(i)/Pstar)**2/(2d0*sigmastar**2))
		endif
		tcrystinv(i)=nucryst*exp(-Tcryst/CloudT(i))
	enddo
	if(Cloud(ii)%hazetype.eq.'optEC') Sn=Sn*scaleUV*Sigmadot/tot

	Sat0=Sat

	if(Cloud(ii)%rainout) then
		allocate(c_rainout(nCS))
		c_rainout=.false.
		xv(1,1:nVS)=0d0
		f=0.5
		i=1
		xv(1:nVS,i)=xv_bot(1:nVS)
		do iter=1,100
1			continue
			do iCS=1,nCS
				tot1=1d200
				iVL(i,iCS)=1
				do iVS=1,nVS
					if(v_cloud(iCS,iVS).gt.0d0.and.v_include(iVS)) then
						tot2=xv_bot(iVS)/v_cloud(iCS,iVS)
						if(tot2.lt.tot1) then
							tot1=tot2
							iVL(i,iCS)=iVS
						endif
					endif
				enddo
				fSat(i,iCS)=CloudP(i)**v_H2(iCS)
				do iVS=1,nVS
					if(v_include(iVS).and.v_cloud(iCS,iVS).ne.0d0) then
						fSat(i,iCS)=fSat(i,iCS)*(CloudP(i)*xv_bot(iVS)*CloudMMW(i)/muV(iVS))**v_cloud(iCS,iVS)
					endif
				enddo
				tot1=Sat(i,iCS)*fSat(i,iCS)/Cloud(ii)%Srainout
				if(tot1.gt.1d0) then
					tot1=xv_bot(iVL(i,iCS))*(1d0-1d0/tot1**f)
					c_rainout(iCS)=.true.
					do iVS=1,nVS
						if(v_cloud(iCS,iVS).gt.0d0.and.v_include(iVS)) then
							xv(iVS,i)=xv(iVS,i)-tot1*v_cloud(iCS,iVS)/v_cloud(iCS,iVL(i,iCS))
							if(xv(iVS,i).lt.0d0) then
								f=f/2d0
								xv(1:nVS,i)=xv_bot(1:nVS)
								goto 1
							endif
						endif
					enddo
				endif
			enddo
			maxerr=0d0
			do iVS=1,nVS
				if((xv_bot(iVS)+xv(iVS,i)).gt.0d0) then
					tot1=abs(xv_bot(iVS)-xv(iVS,i))/(xv_bot(iVS)+xv(iVS,i))
					if(tot1.gt.maxerr) maxerr=tot1
				endif
			enddo
			xv_bot(1:nVS)=xv(1:nVS,i)
			f=(0.5d0+f)/2d0
			if(maxerr.lt.1d-2) exit
		enddo
		do iCS=1,nCS
			if(c_rainout(iCS)) then
				print*,'Partial rainout of ',CSname(iCS)
			endif
		enddo
		deallocate(c_rainout)
	endif
	f=0.1

	allocate(vthv(nnr))
	allocate(Sc(nnr,nCS))
	allocate(IWORKomp(NN))
	allocate(AB(2*NStot+NStot+1,NN))

	iconv=0
c start the loop
	do iter=1,niter
	call tellertje(iter,niter)

	vsed=0d0
	do i=1,nnr
		do iCS=1,nCS
			tot1=1d200
			iVL(i,iCS)=1
			do iVS=1,nVS
				if(v_cloud(iCS,iVS).gt.0d0.and.v_include(iVS)) then
					tot2=xv(iVS,i)*(sqrt(8d0*kb*CloudT(i)/(pi*muV(iVS)*mp)))/(muV(iVS)*v_cloud(iCS,iVS))
					if(tot2.lt.tot1) then
						tot1=tot2
						iVL(i,iCS)=iVS
					endif
				endif
			enddo
			fSat(i,iCS)=CloudP(i)**v_H2(iCS)*(CloudMMW(i)/muV(iVL(i,iCS)))*
     &				(CloudP(i)*xv(iVL(i,iCS),i)*CloudMMW(i)/muV(iVL(i,iCS)))**(v_cloud(iCS,iVL(i,iCS))-1d0)
			do iVS=1,nVS
				if(v_include(iVS).and.iVS.ne.iVL(i,iCS)) then
					fSat(i,iCS)=fSat(i,iCS)*(CloudP(i)*xv(iVS,i)*CloudMMW(i)/muV(iVS))**v_cloud(iCS,iVS)
				endif
			enddo
			if(fSat(i,iCS).lt.1d-40) fSat(i,iCS)=1d-40
		enddo

		cs=sqrt(kb*CloudT(i)/(CloudMMW(i)*mp))
		vth(i)=sqrt(8d0*kb*CloudT(i)/(pi*CloudMMW(i)*mp))
		if(Cloud(ii)%usefsed) then
			fsed=Cloud(ii)%fsed_alpha*exp((CloudR(i)-Rplanet)/(6d0*Cloud(ii)%fsed_beta*CloudHp(i)))+1d-2
			vsed(i)=-fsed*Kd(i)/CloudHp(i)
			rpart(i)=vsed(i)/(-sqrt(pi)*rho_av(i)*CloudG(i)/(2d0*Clouddens(i)*vth(i)))
			if(rpart(i).lt.Cloud(ii)%rnuc) then
				rpart(i)=Cloud(ii)%rnuc
				vsed(i)=-rpart(i)*rho_av(i)*CloudG(i)/(Clouddens(i)*vth(i))
				lmfp=2.3d0*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
				vsed(i)=vsed(i)*max(1d0,4d0*rpart(i)/(9d0*lmfp))
			endif
			mpart(i)=rho_av(i)*4d0*pi*rpart(i)**3/3d0
			if(iter.eq.1) then
				tot=0d0
				do iCS=1,nCS
					do iVS=1,nVS
						tot=tot+xv_bot(iVL(i,iCS))*v_cloud(iCS,iVS)*muV(iVS)/muV(iVL(i,iCS))
					enddo
				enddo
			else
				tot=sum(xc(1:nCS,i))
			endif
			xn(i)=(xn(i)+tot/mpart(i))/2d0
		else
			vsed(i)=-rpart(i)*rho_av(i)*CloudG(i)/(Clouddens(i)*vth(i))
			lmfp=CloudMMW(i)*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
			vsed(i)=vsed(i)*max(1d0,4d0*rpart(i)/(9d0*lmfp))
			mpart(i)=rho_av(i)*4d0*pi*rpart(i)**3/3d0
		endif
	enddo
	Sat=Sat0*fSat
	if(complexKzz) then
		do i=1,nnr
			St=rpart(i)*rho_av(i)*Km(i)/(vth(i)*Clouddens(i)*CloudHp(i)**2)
			Kd(i)=Km(i)/(1d0+St)
		enddo
	endif

C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================
c equations for material
C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================
	do i=1,nnr
		do iCS=1,nCS
			vthv(i)=sqrt(8d0*kb*CloudT(i)/(pi*muV(iVL(i,iCS))*mp))

			Dp=kb*CloudT(i)*vthv(i)/(3d0*CloudP(i)*1d6*sigmamol)
			Sc(i,iCS)=fstick*Clouddens(i)**2*(muC(iCS)/(muV(iVL(i,iCS))*v_cloud(iCS,iVL(i,iCS))))*
     &						pi*rpart(i)*min(rpart(i)*vthv(i),4d0*Dp)
			if(do_nuc(iCS)) then
				tot1=Sat0(i,iCS)*fSat(i,iCS)*xv(iVL(i,iCS),i)
				tot2=CloudP(i)*CloudMMW(i)/(muV(inuc(iCS))*kb*CloudT(i))
				select case(CSname(iCS))
					case('SiO','TiO2')
						call ComputeJ_xv(xv(inuc(iCS),i),tot2,CloudT(i),tot1,Jn_xv(i,iCS),A_J(iCS),B_J(iCS))
						Nc_nuc(i,iCS)=m_nuc/(muC(iCS)*mp)
					case default
						tot2=tot2*xv(inuc(iCS),i)
						call ComputeJ(CloudT(i),tot1,tot2,vthv(i),sigma_nuc(iCS),r0_nuc(iCS),Nf_nuc(iCS),Jn_xv(i,iCS),Nc_nuc(i,iCS))
				end select
			else
				Jn_xv(i,iCS)=0d0
			endif
		enddo
	enddo
	Nc_nuc=Nc_nuc*mp

	x=0d0

	AB=0d0

	KL=NStot
	KU=NStot

	j=0
	i=1

	if(.not.Cloud(ii)%usefsed) then
		j=j+1
		if(Cloud(ii)%freeflow_nuc) then
			x(j)=0d0
	
			dztot=(CloudR(i+1)-CloudR(i))
			dz=(CloudR(i)-CloudR(i+1))
	
			ik=KL+KU+1+j-ixn(i+1)
			AB(ik,ixn(i+1))=AB(ik,ixn(i+1))-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
	
			ik=KL+KU+1+j-ixn(i)
			AB(ik,ixn(i))=AB(ik,ixn(i))+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
			AB(ik,ixn(i))=AB(ik,ixn(i))+Clouddens(i)*vsed(i)/dztot
		else
			ik=KL+KU+1+j-ixn(i)
			AB(ik,ixn(i))=1d0
			x(j)=Cloud(ii)%xm_bot
		endif
	endif
	do iCS=1,nCS
		j=j+1
		if(Cloud(ii)%freeflow_con) then
c assume continuous flux at the bottom (dF/dz=Sc=0)
			x(j)=0d0

			dztot=(CloudR(i+1)-CloudR(i))
			dz=(CloudR(i)-CloudR(i+1))

			ik=KL+KU+1+j-ixc(iCS,i+1)
			AB(ik,ixc(iCS,i+1))=AB(ik,ixc(iCS,i+1))-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot

			ik=KL+KU+1+j-ixc(iCS,i)
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+Clouddens(i)*vsed(i)/dztot

c			if(iamorph(iCS).ne.0) then
c				ik=KL+KU+1+j-ixc(iamorph(iCS),i)
c				AB(ik,ixc(iamorph(iCS),i))=AB(ik,ixc(iamorph(iCS),i))+Clouddens(i)*tcrystinv(i)
c			else
c				ik=KL+KU+1+j-ixv(iVL(i,iCS),i)
c				AB(ik,ixv(iVL(i,iCS),i))=AB(ik,ixv(iVL(i,iCS),i))+Sc(i,iCS)*xn(i)
c			endif
c			if(icryst(iCS).ne.0) then
c				ik=KL+KU+1+j-ixc(iCS,i)
c				AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))-Clouddens(i)*tcrystinv(i)
c			endif
c
c			ik=KL+KU+1+j-ixc(iCS,i)
c			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))-Sc(i,iCS)/(Sat(i,iCS)*mpart(i))
		else
			ik=KL+KU+1+j-ixc(iCS,i)
			AB(ik,ixc(iCS,i))=1d0
			x(j)=0d0
		endif
	enddo
	if(Cloud(ii)%freeflow_con) then
		do iVS=1,nVS
			if(v_include(iVS)) then
				j=j+1
				ik=KL+KU+1+j-ixv(iVS,i)
				AB(ik,ixv(iVS,i))=1d0
				x(j)=xv_bot(iVS)

				do iCS=1,nCS
					ik=KL+KU+1+j-ixc(iCS,i)
					AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+v_cloud(iCS,iVS)*muV(iVS)/muC(iCS)
				enddo
			endif
		enddo
	else
		do iVS=1,nVS
			if(v_include(iVS)) then
				j=j+1
				ik=KL+KU+1+j-ixv(iVS,i)
				AB(ik,ixv(iVS,i))=1d0
				x(j)=xv_bot(iVS)
			endif
		enddo
	endif
	do i=2,nnr-1
		if(.not.Cloud(ii)%usefsed) then
			j=j+1
	
			dztot=(CloudR(i+1)-CloudR(i-1))/2d0
			dz=(CloudR(i)-CloudR(i+1))
			ik=KL+KU+1+j-ixn(i+1)
			AB(ik,ixn(i+1))=AB(ik,ixn(i+1))-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
	
			ik=KL+KU+1+j-ixn(i)
			AB(ik,ixn(i))=AB(ik,ixn(i))+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
	
			dz=(CloudR(i-1)-CloudR(i))
			ik=KL+KU+1+j-ixn(i-1)
			AB(ik,ixn(i-1))=AB(ik,ixn(i-1))-0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz/dztot
	
			ik=KL+KU+1+j-ixn(i)
			AB(ik,ixn(i))=AB(ik,ixn(i))+(0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz)/dztot
			AB(ik,ixn(i))=AB(ik,ixn(i))+Clouddens(i)*vsed(i)/dztot
	
			if(haze_nuc) then
				ik=KL+KU+1+j-ixv(ihaze,i)
				AB(ik,ixv(ihaze,i))=AB(ik,ixv(ihaze,i))+(CloudMMW(i)/muV(ihaze))*Sn(i)/m_nuc
			else
				x(j)=-Sn(i)/m_nuc
			endif
			
			do iCS=1,nCS
				if(do_nuc(iCS)) then
					ik=KL+KU+1+j-ixv(inuc(iCS),i)
					AB(ik,ixv(inuc(iCS),i))=AB(ik,ixv(inuc(iCS),i))+Jn_xv(i,iCS)
				endif
			enddo
c coagulation
			if(Cloud(ii)%coagulation) then
				npart=xn(i)*Clouddens(i)
				lmfp=CloudMMW(i)*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
				vmol=0.5d0*lmfp*vth(i)
				Dp=kb*CloudT(i)/(6d0*pi*rpart(i)*vmol*Clouddens(i))
				vBM=sqrt(16d0*kb*CloudT(i)/(pi*mpart(i)))
				if(Dp/rpart(i).lt.vBM) vBM=Dp/rpart(i)
				tcoaginv=npart*pi*rpart(i)**2*(abs(vsed(i))*exp(-(vsed(i)/vfrag)**2)+2d0*vBM*exp(-(vBM/vfrag)**2))
	
				if(.not.tcoaginv.gt.0d0) tcoaginv=0d0
	
				tcinv(iter,i)=tcoaginv
				tcoaginv=sum(tcinv(max(1,iter-10):iter,i))/real(iter-max(1,iter-10)+1)
	
				ik=KL+KU+1+j-ixn(i)
				AB(ik,ixn(i))=AB(ik,ixn(i))-Clouddens(i)*tcoaginv
			endif
		endif
		do iCS=1,nCS
			j=j+1
	
			dztot=(CloudR(i+1)-CloudR(i-1))/2d0
			dz=(CloudR(i)-CloudR(i+1))
			ik=KL+KU+1+j-ixc(iCS,i+1)
			AB(ik,ixc(iCS,i+1))=AB(ik,ixc(iCS,i+1))-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot

			ik=KL+KU+1+j-ixc(iCS,i)
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot

			dz=(CloudR(i-1)-CloudR(i))
			ik=KL+KU+1+j-ixc(iCS,i-1)
			AB(ik,ixc(iCS,i-1))=AB(ik,ixc(iCS,i-1))-0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz/dztot

			ik=KL+KU+1+j-ixc(iCS,i)
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+(0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz)/dztot
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+Clouddens(i)*vsed(i)/dztot
	
			if(iamorph(iCS).ne.0) then
				ik=KL+KU+1+j-ixc(iamorph(iCS),i)
				AB(ik,ixc(iamorph(iCS),i))=AB(ik,ixc(iamorph(iCS),i))+Clouddens(i)*tcrystinv(i)
			else
				ik=KL+KU+1+j-ixv(iVL(i,iCS),i)
				AB(ik,ixv(iVL(i,iCS),i))=AB(ik,ixv(iVL(i,iCS),i))+Sc(i,iCS)*xn(i)
			endif
			if(icryst(iCS).ne.0) then
				ik=KL+KU+1+j-ixc(iCS,i)
				AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))-Clouddens(i)*tcrystinv(i)
			endif

			if(do_nuc(iCS)) then
				ik=KL+KU+1+j-ixv(inuc(iCS),i)
				AB(ik,ixv(inuc(iCS),i))=AB(ik,ixv(inuc(iCS),i))+Jn_xv(i,iCS)*Nc_nuc(i,iCS)*muC(iCS)
			endif

			ik=KL+KU+1+j-ixc(iCS,i)
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))-Sc(i,iCS)/(Sat(i,iCS)*mpart(i))
		enddo

		do iVS=1,nVS
			if(v_include(iVS)) then
				j=j+1
	
				dztot=(CloudR(i+1)-CloudR(i-1))/2d0
				dz=(CloudR(i)-CloudR(i+1))
				ik=KL+KU+1+j-ixv(iVS,i+1)
				AB(ik,ixv(iVS,i+1))=AB(ik,ixv(iVS,i+1))-(0.5d0*(Kg(i+1)*Clouddens(i+1)+Kg(i)*Clouddens(i))/dz)/dztot

				ik=KL+KU+1+j-ixv(iVS,i)
				AB(ik,ixv(iVS,i))=AB(ik,ixv(iVS,i))+(0.5d0*(Kg(i+1)*Clouddens(i+1)+Kg(i)*Clouddens(i))/dz)/dztot

				dz=(CloudR(i-1)-CloudR(i))
				ik=KL+KU+1+j-ixv(iVS,i-1)
				AB(ik,ixv(iVS,i-1))=AB(ik,ixv(iVS,i-1))-0.5d0*(Kg(i-1)*Clouddens(i-1)+Kg(i)*Clouddens(i))/dz/dztot

				ik=KL+KU+1+j-ixv(iVS,i)
				AB(ik,ixv(iVS,i))=AB(ik,ixv(iVS,i))+(0.5d0*(Kg(i-1)*Clouddens(i-1)+Kg(i)*Clouddens(i))/dz)/dztot
	
				do iCS=1,nCS
					if(iamorph(iCS).eq.0) then
						ik=KL+KU+1+j-ixv(iVL(i,iCS),i)
						AB(ik,ixv(iVL(i,iCS),i))=AB(ik,ixv(iVL(i,iCS),i))-Sc(i,iCS)*xn(i)*v_cloud(iCS,iVS)*(muV(iVS)/muC(iCS))
					endif
					ik=KL+KU+1+j-ixc(iCS,i)
					AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+Sc(i,iCS)*v_cloud(iCS,iVS)*(muV(iVS)/muC(iCS))/(Sat(i,iCS)*mpart(i))

					if(do_nuc(iCS)) then
						ik=KL+KU+1+j-ixv(inuc(iCS),i)
						AB(ik,ixv(inuc(iCS),i))=AB(ik,ixv(inuc(iCS),i))-Jn_xv(i,iCS)*Nc_nuc(i,iCS)*v_cloud(iCS,iVS)*muV(iVS)
					endif
				enddo

				if(haze_nuc.and.iVS.eq.ihaze) then
					ik=KL+KU+1+j-ixv(ihaze,i)
					AB(ik,ixv(ihaze,i))=AB(ik,ixv(ihaze,i))-(CloudMMW(i)/muV(ihaze))*Sn(i)
				endif

			endif
		enddo
	enddo
	i=nnr
	dz=CloudR(i)-CloudR(i-1)

	if(.not.Cloud(ii)%usefsed) then
		j=j+1
		ik=KL+KU+1+j-ixn(i)
		AB(ik,ixn(i))=Kd(i)/dz-vsed(i)
	
		ik=KL+KU+1+j-ixn(i-1)
		AB(ik,ixn(i-1))=-Kd(i)/dz
	
		x(j)=0d0
c		x(j)=F_IDP/((4d0*pi*(0.1e-4)**3*rho_av(i)/3d0)*Clouddens(i))
	endif
	do iCS=1,nCS
		j=j+1
		ik=KL+KU+1+j-ixc(iCS,i)
		AB(ik,ixc(iCS,i))=Kd(i)/dz-vsed(i)

		ik=KL+KU+1+j-ixc(iCS,i-1)
		AB(ik,ixc(iCS,i-1))=-Kd(i)/dz

		x(j)=0d0
c		if(CSname(iCS).eq.'MgSiO3') then
c			x(j)=F_IDP/(Clouddens(i))
c		endif
	enddo
	do iVS=1,nVS
		if(v_include(iVS)) then
			j=j+1
			ik=KL+KU+1+j-ixv(iVS,i)
			AB(ik,ixv(iVS,i))=Kg(i)/dz

			ik=KL+KU+1+j-ixv(iVS,i-1)
			AB(ik,ixv(iVS,i-1))=-Kg(i)/dz

			x(j)=0d0
		endif
	enddo
	
10	continue
	NRHS=1
	info=0
		
	j=2*KL+KU+1

	call DGBSV(NN,KL,KU,NRHS,AB,j,IWORKomp,x,NN,INFO)	

	do i=1,nnr
		if(.not.x(ixn(i)).gt.0d0) x(ixn(i))=0d0
		do iCS=1,nCS
			if(.not.x(ixc(iCS,i)).gt.0d0) x(ixc(iCS,i))=0d0
		enddo
		do iVS=1,nVS
			if(v_include(iVS)) then
				if(.not.x(ixv(iVS,i)).gt.0d0) x(ixv(iVS,i))=0d0
			endif
		enddo
	enddo

	f=0.5
	do i=1,nnr
		if(.not.Cloud(ii)%usefsed) then
			if(iter.eq.1) then
				xn(i)=x(ixn(i))
			else
				xn(i)=(xn(i)*f+x(ixn(i))*(1d0-f))
			endif
		endif
		do iCS=1,nCS
			if(iter.eq.1) then
				xc(iCS,i)=x(ixc(iCS,i))
			else
				xc(iCS,i)=(xc(iCS,i)*f+x(ixc(iCS,i))*(1d0-f))
			endif
		enddo
		do iVS=1,nVS
			if(v_include(iVS)) then
				if(iter.eq.1) then
					xv(iVS,i)=x(ixv(iVS,i))
				else
					xv(iVS,i)=(xv(iVS,i)*f+x(ixv(iVS,i))*(1d0-f))
				endif
			endif
		enddo
	enddo

C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================


	xc_iter(iter,1:nCS,1:nnr)=xc(1:nCS,1:nnr)
	xv_iter(iter,1:nVS,1:nnr)=xv(1:nVS,1:nnr)
	xn_iter(iter,1:nnr)=xn(1:nnr)
	maxerr=0d0
	do i=1,nnr
		tot=0d0
		do iCS=1,nCS
			tot=tot+xc(iCS,i)/rhodust(iCS)
		enddo
		if(tot.gt.0d0) then
			rho_av(i)=sum(xc(1:nCS,i))/tot
		else
			rho_av(i)=sum(rhodust(1:nCS))/real(nCS)
		endif
		tot=sum(xc(1:nCS,i))
		if(xn(i).gt.0d0) then
			if(Cloud(ii)%computeJn) then
				rr=(3d0*(tot/xn(i))/(4d0*pi*rho_av(i)))**(1d0/3d0)
			else
				rr=(3d0*(tot/xn(i))/(4d0*pi*rho_av(i))+Cloud(ii)%rnuc**3)**(1d0/3d0)
			endif
c			if(.not.rr.ge.Cloud(ii)%rnuc) rr=Cloud(ii)%rnuc
		else
			rr=Cloud(ii)%rnuc
		endif
		err=abs(rr-rpart(i))/(rr+rpart(i))
		if(err.gt.maxerr.and.tot.gt.1d-20) then
			maxerr=err
			j=i
		endif
		rpart(i)=sqrt(rr*rpart(i))
		if(Cloud(ii)%computeJn) xn(i)=tot/(4d0*rpart(i)**3*pi*rho_av(i)/3d0)
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
		if(iconv.eq.0.and.nTiter.gt.4) print*,'Cloud formation not converged: ',maxerr
		iter=niter
	endif
	xn=0d0
	xc=0d0
	xv=0d0
	do i=iter-nconv+1,iter
		xn(1:nnr)=xn(1:nnr)+xn_iter(i,1:nnr)/real(nconv)
		xc(1:nCS,1:nnr)=xc(1:nCS,1:nnr)+xc_iter(i,1:nCS,1:nnr)/real(nconv)
		xv(1:nVS,1:nnr)=xv(1:nVS,1:nnr)+xv_iter(i,1:nVS,1:nnr)/real(nconv)
	enddo
	do i=1,nnr
		tot=0d0
		do iCS=1,nCS
			tot=tot+xc(iCS,i)/rhodust(iCS)
		enddo
		if(xn(i).gt.0d0) then
			if(Cloud(ii)%haze) then
				rr=(3d0*(tot/xn(i))/(4d0*pi))**(1d0/3d0)
			else
				rr=(3d0*(tot/xn(i))/(4d0*pi)+Cloud(ii)%rnuc**3)**(1d0/3d0)
			endif
			if(.not.rr.ge.Cloud(ii)%rnuc) then
				rr=Cloud(ii)%rnuc
				xn(i)=3d0*(tot/(rr**3))/(4d0*pi)
			endif
		else
			rr=Cloud(ii)%rnuc
			xn(i)=3d0*(tot/(rr**3))/(4d0*pi)
		endif
		if(tot.gt.0d0) then
			rho_av(i)=sum(xc(1:nCS,i))/tot
		else
			rho_av(i)=sum(rhodust(1:nCS))/real(nCS)
		endif
		rpart(i)=rr
	enddo

	deallocate(vthv)
	deallocate(Sc)
	deallocate(IWORKomp)
	deallocate(AB)

	allocate(dx(nnr))
	logP(1:nr)=-log(P(1:nr))
	logCloudP(1:nnr)=-log(CloudP(1:nnr))

	x(1:nnr)=0d0
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

	if(.not.retrieval) then
		allocate(Jn_out(nCS))
		do i=1,nnr
			do iCS=1,nCS
				if(do_nuc(iCS)) then
					Jn_xv(i,iCS)=Jn_xv(i,iCS)*xv(inuc(iCS),i)*Nc_nuc(i,iCS)*muC(iCS)/m_nuc
					Sn(i)=Sn(i)+Jn_xv(i,iCS)
				endif
			enddo
		enddo
		if(do3D) then
			open(unit=20,file=trim(outputdir) // '/cloudstructure' // trim(int2string(i3D,'(i0.4)')) // '.dat',
     &             FORM="FORMATTED",ACCESS="STREAM")
		else
			open(unit=20,file=trim(outputdir) // '/cloudstructure.dat',FORM="FORMATTED",ACCESS="STREAM")
		endif
		allocate(v_names_out(nVS),xv_out(nVS))
		j=0
		do iVS=1,nVS
			if(v_include(iVS)) then
				j=j+1
				v_names_out(j)=v_names(iVS)
			endif
		enddo
		form='("#",a18,a19,a19,' // trim(int2string(nCS+j+1,'(i4)')) // 'a23,a19,a19,a19,a19)'
		write(20,form) "P[bar]","dens[g/cm^3]","xn",(trim(CSname(i))//"[s]",i=1,nCS),
     &				(trim(v_names_out(i))//"[v]",i=1,j),"r[cm]","T[K]","Jstar","fsed"
		form='(es19.7E3,es19.7E3,es19.7E3,' // trim(int2string(nCS+j,'(i4)')) // 'es23.7E3,es19.7E3,es19.7E3,es19.7E3,es19.7E3)'
		do i=1,nnr
			j=0
			do iVS=1,nVS
				if(v_include(iVS)) then
					j=j+1
					xv_out(j)=xv(iVS,i)
				endif
			enddo
			write(20,form) CloudP(i),Clouddens(i),xn(i)*m_nuc,xc(1:nCS,i),xv_out(1:j),rpart(i),
     &							CloudT(i),Sn(i),-vsed(i)*CloudHp(i)/Kd(i)
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
	do iVS=1,nVS
		if(v_include(iVS)) then
			x(1:nnr)=xv(iVS,1:nnr)
			call regridarray(logCloudP,x,nnr,logP,logx,nr)
			do j=1,N_atoms
				at_ab(1:nr,j)=at_ab(1:nr,j)+logx(1:nr)*mutot*v_atoms(iVS,j)/muV(iVS)
			enddo
		endif
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
			call AddDiseqAtoms(Otot,Ctot,Ntot)
c			call disequilibrium code
c			input: 	R(1:nr+1) : These are the radial boundaries of the layers (bottom to top)
c			P(1:nr),T(1:nr) : These are the pressure and temperature inside the layers
c			molname(1:nmol) : names of the molecules included
c			Kzz_r(1:nr) : Diffusion coefficient
c			input/output:	mixrat_r(1:nr,1:nmol) : number densities inside each layer. Now set to equilibrium abundances.
			call output("==================================================================")
			call output("Computing disequilibrium chemistry")
			call diseq_calc(nr,R(1:nr+1),P(1:nr),T(1:nr),nmol,molname(1:nmol),mixrat_r(1:nr, 1:nmol),COratio,Kzz_g(1:nr))
			call CorrectDiseqAtoms(Otot,Ctot,Ntot)
		endif
	else
		do iVS=1,nVS
			if(v_include(iVS)) then
				do k=1,nmol
					if(includemol(k)) then
						if(v_names(iVS).eq.molname(k)) then
							x(1:nnr)=xv(iVS,1:nnr)
							call regridarray(logCloudP,x,nnr,logP,logx,nr)
							mixrat_r(1:nr,k)=logx(1:nr)*mutot/(muV(iVS))
						endif
					endif
				enddo
			endif
		enddo
	endif
	do i=1,nr
		tot=0d0
		do j=1,nmol
			if(mixrat_r(i,j).gt.0d0) tot=tot+mixrat_r(i,j)
		enddo
		if(tot.gt.0d0) mixrat_r(i,1:nmol)=mixrat_r(i,1:nmol)/tot
	enddo
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


	deallocate(Sat,docondense)
	deallocate(mpart)
	deallocate(rho_av)
	deallocate(y)
	deallocate(Sn)
	deallocate(vth)
	deallocate(tcinv,xn_iter,xc_iter,xv_iter)
	deallocate(vsed)
	deallocate(ixv)
	deallocate(ixc)
	deallocate(x,dx)
	deallocate(IWORK)
	deallocate(logCloudP)
	deallocate(Kd,Kg,Km)
	if(Cloud(ii)%hazetype.eq.'optEC') deallocate(CloudtauUV,CloudkappaUV)

	return
	end

	subroutine ComputeJ_xv(xv,scale,T,Sat,J,A,B)
	IMPLICIT NONE
	real*8 xv,scale,T,Sat,J,A,B
	
	if(Sat.lt.1d0) then
		J=0d0
	else
		J=xv*scale**2*exp(B-A/(T**3*(log(Sat))**2))
	endif
	if(.not.J.gt.0d0) J=0d0
	
	return
	end


	subroutine ComputeJ(T,Sat,nx,vth,sigma,r0,Nf,J,Nc)
	IMPLICIT NONE
	real*8 T,nx,sigma,Sat,vth,J
	real*8 Nstar1,theta,r0,A0,logSat
	real*8 Tmax,mu,tau,xv,MMW
	real*8 pi,Z,dGRT,Nstar_inf,Nf,kb,mp,Nc
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(kb=1.3806503d-16)
	parameter(mp=1.660539040d-24)	!atomic mass unit
	
	logSat=log(Sat)

	A0=4d0*pi*r0**2

	theta=A0*sigma/(kb*T)

	Nstar_inf=(2d0*theta/(3d0*logSat))**3
	Nstar1=(Nstar_inf/8d0)*(1d0+sqrt(1d0+2d0*(Nf/Nstar_inf)**(1d0/3d0))-2d0*(Nf/Nstar_inf)**(1d0/3d0))**3

	Z=sqrt(theta*(2d0*Nstar1**(1d0/3d0)+3d0*Nf**(1d0/3d0))/(6d0*pi))/(Nstar1**(1d0/3d0)+Nf**(1d0/3d0))

	dGRT=theta*Nstar1/(Nstar1**(1d0/3d0)+Nf**(1d0/3d0))

	tau=nx*vth*(Nstar1+1d0)**(2./3.)*A0

	J=(nx*tau)*Z*exp(Nstar1*logSat-dGRT)
	Nc=Nstar1+1d0

	if(.not.J.gt.0d0) then
		J=0d0
		Nc=1d0
	endif

	return	
	end


	
	subroutine AddDiseqAtoms(Otot,Ctot,Ntot)	
	use GlobalSetup
	IMPLICIT NONE
	real*8 Otot(nr),Ctot(nr),Ntot(nr)
	integer iCO,iCO2,iH2O,iCH4,iN2,iNH3
	
	iH2O=1
	iCO2=2
	iCO=5
	iCH4=6
	iNH3=11
	iN2=22
	
	Otot(1:nr)=mixrat_r(1:nr,iH2O)+2.0*mixrat_r(1:nr,iCO2)+mixrat_r(1:nr,iCO)
	Ctot(1:nr)=mixrat_r(1:nr,iCO2)+mixrat_r(1:nr,iCO)+mixrat_r(1:nr,iCH4)
	Ntot(1:nr)=mixrat_r(1:nr,iNH3)+2.0*mixrat_r(1:nr,iN2)
	
	return
	end

	subroutine CorrectDiseqAtoms(Otot,Ctot,Ntot)	
	use GlobalSetup
	IMPLICIT NONE
	real*8 Otot(nr),Ctot(nr),Ntot(nr),f
	integer iCO,iCO2,iH2O,iCH4,iN2,iNH3
	integer i
	
	iH2O=1
	iCO2=2
	iCO=5
	iCH4=6
	iNH3=11
	iN2=22
	
	do i=1,nr
		f=(Otot(i)-mixrat_r(i,iCO))/(2.0*mixrat_r(i,iCO2)+mixrat_r(i,iH2O))
		if(f.lt.0d0) f=0d0
		mixrat_r(i,iCO2)=f*mixrat_r(i,iCO2)
		mixrat_r(i,iH2O)=f*mixrat_r(i,iH2O)
		f=(Ctot(i)-mixrat_r(i,iCO2))/(mixrat_r(i,iCH4)+mixrat_r(i,iCO))
		if(f.lt.0d0) f=0d0
		mixrat_r(i,iCO)=f*mixrat_r(i,iCO)
		mixrat_r(i,iCH4)=f*mixrat_r(i,iCH4)
		f=Ntot(i)/(mixrat_r(i,iNH3)+2.0*mixrat_r(i,iN2))
		mixrat_r(i,iNH3)=f*mixrat_r(i,iNH3)
		mixrat_r(i,iN2)=f*mixrat_r(i,iN2)
	enddo
	
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
	

	subroutine PvapSiO(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0,c1
	parameter(c0=-4.95200E+04,c1=3.25200E+01)

	Pvap=1d-6*exp(c0/T+c1)
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



	subroutine Gibbs_FeS_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          39),Ggibbs(          39)
	parameter(nG=39)
	data (Tgibbs(j),j=1,          39) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 
     &  0.25000E+04, 0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 
     &  0.30000E+04, 0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 
     &  0.35000E+04, 0.36000E+04, 0.37000E+04, 0.38000E+04 /
	data (Ggibbs(j),j=1,          39) /
     & -0.10196E+03,-0.10190E+03,-0.10195E+03,-0.10195E+03,-0.10209E+03, 
     & -0.10249E+03,-0.10303E+03,-0.10364E+03,-0.10415E+03,-0.10344E+03, 
     & -0.97670E+02,-0.91766E+02,-0.85786E+02,-0.79802E+02,-0.73918E+02, 
     & -0.68943E+02,-0.65382E+02,-0.61837E+02,-0.58273E+02,-0.54008E+02, 
     & -0.49656E+02,-0.45287E+02,-0.40900E+02,-0.36496E+02,-0.32076E+02, 
     & -0.27639E+02,-0.23187E+02,-0.18718E+02,-0.14235E+02,-0.97360E+01, 
     & -0.52230E+01,-0.69500E+00, 0.11270E+02, 0.26912E+02, 0.42510E+02, 
     &  0.58067E+02, 0.73585E+02, 0.89066E+02, 0.10451E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Fe_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          64),Ggibbs(          64)
	parameter(nG=64)
	data (Tgibbs(j),j=1,          64) /
     &  0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 0.41000E+04, 
     &  0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 0.46000E+04, 
     &  0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 0.51000E+04, 
     &  0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 0.56000E+04, 
     &  0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          64) /
     &  0.39993E+03, 0.38481E+03, 0.37718E+03, 0.36980E+03, 0.36952E+03, 
     &  0.36186E+03, 0.35421E+03, 0.34656E+03, 0.33894E+03, 0.32376E+03, 
     &  0.30872E+03, 0.29383E+03, 0.27913E+03, 0.26465E+03, 0.25050E+03, 
     &  0.23659E+03, 0.22288E+03, 0.20926E+03, 0.19574E+03, 0.18232E+03, 
     &  0.16901E+03, 0.15584E+03, 0.14348E+03, 0.13131E+03, 0.11925E+03, 
     &  0.10730E+03, 0.95456E+02, 0.83705E+02, 0.72043E+02, 0.60467E+02, 
     &  0.48972E+02, 0.37554E+02, 0.26209E+02, 0.14933E+02, 0.37240E+01, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Fe_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          32),Ggibbs(          32)
	parameter(nG=          32)
	data (Tgibbs(j),j=1,          32) /
     &  0.29815E+03, 0.30000E+03, 0.35000E+03, 0.40000E+03, 0.45000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04, 0.21000E+04, 0.22000E+04, 0.32000E+04, 0.33000E+04, 
     &  0.34000E+04, 0.35000E+04, 0.36000E+04, 0.37000E+04, 0.38000E+04, 
     &  0.39000E+04, 0.40000E+04 /
	data (Ggibbs(j),j=1,          32) /
     & -0.53360E+01,-0.53210E+01,-0.49050E+01,-0.44840E+01,-0.40640E+01, 
     & -0.36480E+01,-0.28420E+01,-0.20860E+01,-0.13980E+01,-0.80200E+00, 
     & -0.33800E+00,-0.91000E-01,-0.10000E-01, 0.52000E-01, 0.65000E-01, 
     &  0.56000E-01, 0.29000E-01, 0.19000E-01, 0.85000E-01, 0.87300E+00, 
     &  0.17610E+01, 0.26760E+01, 0.36120E+01, 0.74240E+01, 0.18511E+02, 
     &  0.29541E+02, 0.40516E+02, 0.51439E+02, 0.62313E+02, 0.73140E+02, 
     &  0.83921E+02, 0.94660E+02 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_H2O_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          61),Ggibbs(          61)
	parameter(nG=61)
	data (Tgibbs(j),j=1,          61) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 
     &  0.25000E+04, 0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 
     &  0.30000E+04, 0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 
     &  0.35000E+04, 0.36000E+04, 0.37000E+04, 0.38000E+04, 0.39000E+04, 
     &  0.40000E+04, 0.41000E+04, 0.42000E+04, 0.43000E+04, 0.44000E+04, 
     &  0.45000E+04, 0.46000E+04, 0.47000E+04, 0.48000E+04, 0.49000E+04, 
     &  0.50000E+04, 0.51000E+04, 0.52000E+04, 0.53000E+04, 0.54000E+04, 
     &  0.55000E+04, 0.56000E+04, 0.57000E+04, 0.58000E+04, 0.59000E+04, 
     &  0.60000E+04 /
	data (Ggibbs(j),j=1,          61) /
     & -0.23658E+03,-0.23277E+03,-0.22858E+03,-0.22850E+03,-0.22390E+03, 
     & -0.21905E+03,-0.21401E+03,-0.20881E+03,-0.20350E+03,-0.19808E+03, 
     & -0.19259E+03,-0.18703E+03,-0.18143E+03,-0.17577E+03,-0.17009E+03, 
     & -0.16438E+03,-0.15864E+03,-0.15288E+03,-0.14711E+03,-0.14132E+03, 
     & -0.13553E+03,-0.12972E+03,-0.12391E+03,-0.11808E+03,-0.11225E+03, 
     & -0.10642E+03,-0.10058E+03,-0.94729E+02,-0.88878E+02,-0.83023E+02, 
     & -0.77163E+02,-0.71298E+02,-0.65430E+02,-0.59558E+02,-0.53681E+02, 
     & -0.47801E+02,-0.41916E+02,-0.36027E+02,-0.30133E+02,-0.24236E+02, 
     & -0.18334E+02,-0.12427E+02,-0.65160E+01,-0.60000E+00, 0.53200E+01, 
     &  0.11245E+02, 0.17175E+02, 0.23111E+02, 0.29052E+02, 0.34998E+02, 
     &  0.40949E+02, 0.46906E+02, 0.52869E+02, 0.58838E+02, 0.64811E+02, 
     &  0.70791E+02, 0.76777E+02, 0.82769E+02, 0.88767E+02, 0.94770E+02, 
     &  0.10078E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_H2O_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          12),Ggibbs(          12)
	parameter(nG=12)
	data (Tgibbs(j),j=1,          12) /
     &   0.29815E+03, 0.30000E+03, 0.32000E+03, 0.34000E+03, 
     &  0.36000E+03, 0.38000E+03, 0.40000E+03, 0.42000E+03, 0.44000E+03, 
     &  0.46000E+03, 0.48000E+03, 0.50000E+03 /
	data (Ggibbs(j),j=1,          12) /
     &  -0.23714E+03,-0.23684E+03,-0.23360E+03,-0.23040E+03, 
     & -0.22723E+03,-0.22410E+03,-0.22101E+03,-0.21794E+03,-0.21491E+03, 
     & -0.21191E+03,-0.20894E+03,-0.20600E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_H2S_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          61),Ggibbs(          61)
	parameter(nG=61)
	data (Tgibbs(j),j=1,          61) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 
     &  0.25000E+04, 0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 
     &  0.30000E+04, 0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 
     &  0.35000E+04, 0.36000E+04, 0.37000E+04, 0.38000E+04, 0.39000E+04, 
     &  0.40000E+04, 0.41000E+04, 0.42000E+04, 0.43000E+04, 0.44000E+04, 
     &  0.45000E+04, 0.46000E+04, 0.47000E+04, 0.48000E+04, 0.49000E+04, 
     &  0.50000E+04, 0.51000E+04, 0.52000E+04, 0.53000E+04, 0.54000E+04, 
     &  0.55000E+04, 0.56000E+04, 0.57000E+04, 0.58000E+04, 0.59000E+04, 
     &  0.60000E+04 /
	data (Ggibbs(j),j=1,          61) /
     & -0.23519E+02,-0.28754E+02,-0.33329E+02,-0.33408E+02,-0.37343E+02, 
     & -0.40179E+02,-0.42399E+02,-0.44201E+02,-0.45694E+02,-0.45870E+02, 
     & -0.40984E+02,-0.36066E+02,-0.31129E+02,-0.26179E+02,-0.21222E+02, 
     & -0.16261E+02,-0.11298E+02,-0.63360E+01,-0.13740E+01, 0.35870E+01, 
     &  0.85470E+01, 0.13505E+02, 0.18461E+02, 0.23416E+02, 0.28370E+02, 
     &  0.33323E+02, 0.38275E+02, 0.43227E+02, 0.48178E+02, 0.53129E+02, 
     &  0.58080E+02, 0.63031E+02, 0.67982E+02, 0.72935E+02, 0.77887E+02, 
     &  0.82841E+02, 0.87796E+02, 0.92752E+02, 0.97710E+02, 0.10267E+03, 
     &  0.10763E+03, 0.11259E+03, 0.11756E+03, 0.12253E+03, 0.12750E+03, 
     &  0.13247E+03, 0.13745E+03, 0.14243E+03, 0.14742E+03, 0.15240E+03, 
     &  0.15739E+03, 0.16239E+03, 0.16739E+03, 0.17239E+03, 0.17740E+03, 
     &  0.18241E+03, 0.18743E+03, 0.19245E+03, 0.19748E+03, 0.20251E+03, 
     &  0.20754E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_MgSiO3_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          31),Ggibbs(          31)
	parameter(nG=31)
	data (Tgibbs(j),j=1,          31) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 
     &  0.25000E+04, 0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 
     &  0.30000E+04 /
	data (Ggibbs(j),j=1,          31) /
     & -0.15181E+04,-0.14904E+04,-0.14620E+04,-0.14615E+04,-0.14323E+04, 
     & -0.14031E+04,-0.13740E+04,-0.13449E+04,-0.13160E+04,-0.12873E+04, 
     & -0.12580E+04,-0.12285E+04,-0.11991E+04,-0.11698E+04,-0.11375E+04, 
     & -0.10991E+04,-0.10608E+04,-0.10222E+04,-0.98118E+03,-0.94229E+03, 
     & -0.90565E+03,-0.86922E+03,-0.83298E+03,-0.79693E+03,-0.76105E+03, 
     & -0.72535E+03,-0.68980E+03,-0.65441E+03,-0.61917E+03,-0.58407E+03, 
     & -0.54909E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Mg_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          64),Ggibbs(          64)
	parameter(nG=64)
	data (Tgibbs(j),j=1,          64) /
     &  0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 0.41000E+04, 
     &  0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 0.46000E+04, 
     &  0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 0.51000E+04, 
     &  0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 0.56000E+04, 
     &  0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          64) /
     &  0.13569E+03, 0.12397E+03, 0.11812E+03, 0.11252E+03, 0.11231E+03, 
     &  0.10653E+03, 0.10078E+03, 0.95067E+02, 0.89387E+02, 0.78122E+02, 
     &  0.66981E+02, 0.55961E+02, 0.45062E+02, 0.35000E+02, 0.25283E+02, 
     &  0.15690E+02, 0.62090E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_SiO2_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          46),Ggibbs(          46)
	parameter(nG=46)
	data (Tgibbs(j),j=1,          46) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 
     &  0.25000E+04, 0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 
     &  0.30000E+04, 0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 
     &  0.35000E+04, 0.36000E+04, 0.37000E+04, 0.38000E+04, 0.39000E+04, 
     &  0.40000E+04, 0.41000E+04, 0.42000E+04, 0.43000E+04, 0.44000E+04, 
     &  0.45000E+04 /
	data (Ggibbs(j),j=1,          46) /
     & -0.89166E+03,-0.87422E+03,-0.85644E+03,-0.85611E+03,-0.83781E+03, 
     & -0.81954E+03,-0.80137E+03,-0.78333E+03,-0.76545E+03,-0.74778E+03, 
     & -0.73026E+03,-0.71280E+03,-0.69543E+03,-0.67811E+03,-0.66087E+03, 
     & -0.64368E+03,-0.62655E+03,-0.60906E+03,-0.58956E+03,-0.57018E+03, 
     & -0.55091E+03,-0.53174E+03,-0.51267E+03,-0.49370E+03,-0.47481E+03, 
     & -0.45600E+03,-0.43728E+03,-0.41863E+03,-0.40005E+03,-0.38154E+03, 
     & -0.36309E+03,-0.34471E+03,-0.32638E+03,-0.30812E+03,-0.28991E+03, 
     & -0.27175E+03,-0.24318E+03,-0.21417E+03,-0.18521E+03,-0.15631E+03, 
     & -0.12747E+03,-0.98680E+02,-0.69941E+02,-0.41252E+02,-0.12612E+02, 
     &  0.15981E+02 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_SiO_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          64),Ggibbs(          64)
	parameter(nG=64)
	data (Tgibbs(j),j=1,          64) /
     &  0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 0.41000E+04, 
     &  0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 0.46000E+04, 
     &  0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 0.51000E+04, 
     &  0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 0.56000E+04, 
     &  0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          64) /
     & -0.10930E+03,-0.11839E+03,-0.12294E+03,-0.12731E+03,-0.12747E+03, 
     & -0.13196E+03,-0.13641E+03,-0.14082E+03,-0.14519E+03,-0.15384E+03, 
     & -0.16238E+03,-0.17082E+03,-0.17917E+03,-0.18743E+03,-0.19561E+03, 
     & -0.20371E+03,-0.21174E+03,-0.21971E+03,-0.22760E+03,-0.23542E+03, 
     & -0.24274E+03,-0.24746E+03,-0.25214E+03,-0.25677E+03,-0.26136E+03, 
     & -0.26591E+03,-0.27041E+03,-0.27488E+03,-0.27931E+03,-0.28370E+03, 
     & -0.28806E+03,-0.29238E+03,-0.29668E+03,-0.30094E+03,-0.30516E+03, 
     & -0.30936E+03,-0.31353E+03,-0.31767E+03,-0.32178E+03,-0.31540E+03, 
     & -0.30850E+03,-0.30159E+03,-0.29466E+03,-0.28772E+03,-0.28076E+03, 
     & -0.27379E+03,-0.26681E+03,-0.25981E+03,-0.25280E+03,-0.24578E+03, 
     & -0.23874E+03,-0.23170E+03,-0.22464E+03,-0.21757E+03,-0.21049E+03, 
     & -0.20340E+03,-0.19630E+03,-0.18919E+03,-0.18207E+03,-0.17494E+03, 
     & -0.16780E+03,-0.16065E+03,-0.15350E+03,-0.14633E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Al2O3_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          41),Ggibbs(          41)
	parameter(nG=          41)
	data (Tgibbs(j),j=1,          41) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 
     &  0.25000E+04, 0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 
     &  0.30000E+04, 0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 
     &  0.35000E+04, 0.36000E+04, 0.37000E+04, 0.38000E+04, 0.39000E+04, 
     &  0.40000E+04 /
	data (Ggibbs(j),j=1,          41) /
     & -0.16416E+04,-0.16127E+04,-0.15823E+04,-0.15817E+04,-0.15502E+04, 
     & -0.15187E+04,-0.14873E+04,-0.14561E+04,-0.14249E+04,-0.13939E+04, 
     & -0.13614E+04,-0.13283E+04,-0.12952E+04,-0.12623E+04,-0.12294E+04, 
     & -0.11966E+04,-0.11639E+04,-0.11313E+04,-0.10988E+04,-0.10664E+04, 
     & -0.10341E+04,-0.10018E+04,-0.96968E+03,-0.93759E+03,-0.90913E+03, 
     & -0.88224E+03,-0.85563E+03,-0.82929E+03,-0.80128E+03,-0.75443E+03, 
     & -0.70789E+03,-0.66166E+03,-0.61572E+03,-0.57006E+03,-0.52468E+03, 
     & -0.47956E+03,-0.43470E+03,-0.39009E+03,-0.34571E+03,-0.30157E+03, 
     & -0.25766E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Al_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          64),Ggibbs(          64)
	parameter(nG=          64)
	data (Tgibbs(j),j=1,          64) /
     &  0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 0.41000E+04, 
     &  0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 0.46000E+04, 
     &  0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 0.51000E+04, 
     &  0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 0.56000E+04, 
     &  0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          64) /
     &  0.31603E+03, 0.30248E+03, 0.29564E+03, 0.28907E+03, 0.28882E+03, 
     &  0.28202E+03, 0.27524E+03, 0.26850E+03, 0.26179E+03, 0.24844E+03, 
     &  0.23522E+03, 0.22211E+03, 0.20913E+03, 0.19703E+03, 0.18542E+03, 
     &  0.17392E+03, 0.16250E+03, 0.15117E+03, 0.13992E+03, 0.12874E+03, 
     &  0.11763E+03, 0.10658E+03, 0.95600E+02, 0.84673E+02, 0.73800E+02, 
     &  0.62979E+02, 0.52209E+02, 0.41486E+02, 0.30808E+02, 0.20175E+02, 
     &  0.95830E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_HCl_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          64),Ggibbs(          64)
	parameter(nG=          64)
	data (Tgibbs(j),j=1,          64) /
     &  0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 0.41000E+04, 
     &  0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 0.46000E+04, 
     &  0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 0.51000E+04, 
     &  0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 0.56000E+04, 
     &  0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          64) /
     & -0.93236E+02,-0.94284E+02,-0.94809E+02,-0.95300E+02,-0.95318E+02, 
     & -0.95809E+02,-0.96280E+02,-0.96732E+02,-0.97166E+02,-0.97985E+02, 
     & -0.98747E+02,-0.99465E+02,-0.10015E+03,-0.10080E+03,-0.10143E+03, 
     & -0.10204E+03,-0.10264E+03,-0.10323E+03,-0.10381E+03,-0.10439E+03, 
     & -0.10496E+03,-0.10552E+03,-0.10608E+03,-0.10663E+03,-0.10718E+03, 
     & -0.10773E+03,-0.10827E+03,-0.10881E+03,-0.10935E+03,-0.10988E+03, 
     & -0.11041E+03,-0.11093E+03,-0.11145E+03,-0.11197E+03,-0.11248E+03, 
     & -0.11299E+03,-0.11349E+03,-0.11399E+03,-0.11448E+03,-0.11497E+03, 
     & -0.11545E+03,-0.11593E+03,-0.11640E+03,-0.11686E+03,-0.11733E+03, 
     & -0.11778E+03,-0.11823E+03,-0.11868E+03,-0.11912E+03,-0.11956E+03, 
     & -0.11999E+03,-0.12042E+03,-0.12084E+03,-0.12126E+03,-0.12168E+03, 
     & -0.12209E+03,-0.12250E+03,-0.12291E+03,-0.12332E+03,-0.12372E+03, 
     & -0.12412E+03,-0.12452E+03,-0.12492E+03,-0.12532E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_KCl_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          21),Ggibbs(          21)
	parameter(nG=          21)
	data (Tgibbs(j),j=1,          21) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04 /
	data (Ggibbs(j),j=1,          21) /
     & -0.42778E+03,-0.41805E+03,-0.40876E+03,-0.40859E+03,-0.39884E+03, 
     & -0.38894E+03,-0.37915E+03,-0.36950E+03,-0.35999E+03,-0.35062E+03, 
     & -0.34142E+03,-0.32920E+03,-0.31536E+03,-0.30180E+03,-0.28850E+03, 
     & -0.27545E+03,-0.26262E+03,-0.25000E+03,-0.23759E+03,-0.22535E+03, 
     & -0.21330E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_K_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          64),Ggibbs(          64)
	parameter(nG=          64)
	data (Tgibbs(j),j=1,          64) /
     &  0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 0.41000E+04, 
     &  0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 0.46000E+04, 
     &  0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 0.51000E+04, 
     &  0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 0.56000E+04, 
     &  0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          64) /
     &  0.80028E+02, 0.70006E+02, 0.65116E+02, 0.60476E+02, 0.60299E+02, 
     &  0.55652E+02, 0.51335E+02, 0.47086E+02, 0.42894E+02, 0.34652E+02, 
     &  0.26568E+02, 0.18614E+02, 0.10773E+02, 0.30340E+01, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Na2S_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          29),Ggibbs(          29)
	parameter(nG=          29)
	data (Tgibbs(j),j=1,          29) /
     &  0.00000E+00, 0.30000E+03, 0.40000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04 /
	data (Ggibbs(j),j=1,          29) /
     &  0.96232E+02,-0.35448E+03,-0.35013E+03,-0.34409E+03,-0.33776E+03, 
     & -0.33132E+03,-0.32484E+03,-0.31728E+03,-0.30488E+03,-0.29265E+03, 
     & -0.27586E+03,-0.24818E+03,-0.22143E+03,-0.19594E+03,-0.17130E+03, 
     & -0.14685E+03,-0.12258E+03,-0.98480E+02,-0.74541E+02,-0.50753E+02, 
     & -0.27110E+02,-0.36040E+01, 0.19772E+02, 0.43024E+02, 0.66157E+02, 
     &  0.89176E+02, 0.11209E+03, 0.13489E+03, 0.15760E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_NaCl_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          26),Ggibbs(          26)
	parameter(nG=          26)
	data (Tgibbs(j),j=1,          26) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 
     &  0.25000E+04 /
	data (Ggibbs(j),j=1,          26) /
     & -0.40246E+03,-0.39305E+03,-0.38402E+03,-0.38386E+03,-0.37464E+03, 
     & -0.36501E+03,-0.35549E+03,-0.34610E+03,-0.33683E+03,-0.32772E+03, 
     & -0.31877E+03,-0.31067E+03,-0.30226E+03,-0.28823E+03,-0.27441E+03, 
     & -0.26079E+03,-0.24735E+03,-0.23408E+03,-0.22096E+03,-0.20800E+03, 
     & -0.19518E+03,-0.18250E+03,-0.16994E+03,-0.15750E+03,-0.14518E+03, 
     & -0.13297E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Na_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          64),Ggibbs(          64)
	parameter(nG=          64)
	data (Tgibbs(j),j=1,          64) /
     &  0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 0.41000E+04, 
     &  0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 0.46000E+04, 
     &  0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 0.51000E+04, 
     &  0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 0.56000E+04, 
     &  0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          64) /
     &  0.97574E+02, 0.86977E+02, 0.81776E+02, 0.76825E+02, 0.76636E+02, 
     &  0.71560E+02, 0.66753E+02, 0.62161E+02, 0.57626E+02, 0.48697E+02, 
     &  0.39919E+02, 0.31263E+02, 0.22710E+02, 0.14246E+02, 0.58650E+01, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Mg2SiO4_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          41),Ggibbs(          41)
	parameter(nG=          41)
	data (Tgibbs(j),j=1,          41) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 
     &  0.25000E+04, 0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 
     &  0.30000E+04, 0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 
     &  0.35000E+04, 0.36000E+04, 0.37000E+04, 0.38000E+04, 0.39000E+04, 
     &  0.40000E+04 /
	data (Ggibbs(j),j=1,          41) /
     & -0.21345E+04,-0.20968E+04,-0.20579E+04,-0.20571E+04,-0.20171E+04, 
     & -0.19772E+04,-0.19374E+04,-0.18978E+04,-0.18584E+04,-0.18191E+04, 
     & -0.17786E+04,-0.17377E+04,-0.16970E+04,-0.16564E+04,-0.16095E+04, 
     & -0.15506E+04,-0.14919E+04,-0.14331E+04,-0.13720E+04,-0.13113E+04, 
     & -0.12508E+04,-0.11906E+04,-0.11316E+04,-0.10752E+04,-0.10190E+04, 
     & -0.96315E+03,-0.90750E+03,-0.85206E+03,-0.79683E+03,-0.74182E+03, 
     & -0.68699E+03,-0.63235E+03,-0.57790E+03,-0.52361E+03,-0.46949E+03, 
     & -0.41553E+03,-0.35125E+03,-0.28664E+03,-0.22217E+03,-0.15784E+03, 
     & -0.93665E+02 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_TiO_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          64),Ggibbs(          64)
	parameter(nG=          64)
	data (Tgibbs(j),j=1,          64) /
     &  0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 0.41000E+04, 
     &  0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 0.46000E+04, 
     &  0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 0.51000E+04, 
     &  0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 0.56000E+04, 
     &  0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          64) /
     &  0.44842E+02, 0.34491E+02, 0.29387E+02, 0.24535E+02, 0.24350E+02, 
     &  0.19375E+02, 0.14456E+02, 0.95860E+01, 0.47600E+01,-0.47710E+01, 
     & -0.14161E+02,-0.23425E+02,-0.32573E+02,-0.41609E+02,-0.50527E+02, 
     & -0.59198E+02,-0.67545E+02,-0.75814E+02,-0.84005E+02,-0.92116E+02, 
     & -0.10014E+03,-0.10809E+03,-0.11594E+03,-0.12326E+03,-0.13015E+03, 
     & -0.13691E+03,-0.14355E+03,-0.15007E+03,-0.15647E+03,-0.16276E+03, 
     & -0.16895E+03,-0.17504E+03,-0.18103E+03,-0.18693E+03,-0.19274E+03, 
     & -0.19846E+03,-0.20410E+03,-0.20966E+03,-0.21515E+03,-0.22056E+03, 
     & -0.21811E+03,-0.21213E+03,-0.20611E+03,-0.20006E+03,-0.19398E+03, 
     & -0.18787E+03,-0.18172E+03,-0.17554E+03,-0.16932E+03,-0.16307E+03, 
     & -0.15679E+03,-0.15048E+03,-0.14414E+03,-0.13776E+03,-0.13136E+03, 
     & -0.12492E+03,-0.11845E+03,-0.11196E+03,-0.10543E+03,-0.98873E+02, 
     & -0.92288E+02,-0.85674E+02,-0.79032E+02,-0.72361E+02 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_TiO2_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          39),Ggibbs(          39)
	parameter(nG=          39)
	data (Tgibbs(j),j=1,          39) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04 /
	data (Ggibbs(j),j=1,          39) /
     & -0.92551E+03,-0.90758E+03,-0.88941E+03,-0.85216E+03,-0.83394E+03, 
     & -0.81587E+03,-0.79794E+03,-0.78013E+03,-0.76243E+03,-0.74480E+03, 
     & -0.72711E+03,-0.70926E+03,-0.69150E+03,-0.67380E+03,-0.65616E+03, 
     & -0.63857E+03,-0.62101E+03,-0.60349E+03,-0.58553E+03,-0.56726E+03, 
     & -0.55119E+03,-0.53612E+03,-0.52112E+03,-0.50618E+03,-0.49129E+03, 
     & -0.47646E+03,-0.46168E+03,-0.44695E+03,-0.43227E+03,-0.41763E+03, 
     & -0.40303E+03,-0.38848E+03,-0.37396E+03,-0.35948E+03,-0.34504E+03, 
     & -0.32284E+03,-0.29721E+03,-0.27164E+03,-0.24614E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_H2SO4_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          11),Ggibbs(          11)
	parameter(nG=          11)
	data (Tgibbs(j),j=1,          11) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04 /
	data (Ggibbs(j),j=1,          11) /
     & -0.77730E+03,-0.73377E+03,-0.68992E+03,-0.68915E+03,-0.64797E+03, 
     & -0.60725E+03,-0.56748E+03,-0.52885E+03,-0.49149E+03,-0.45438E+03, 
     & -0.41379E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Zn_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          64),Ggibbs(          64)
	parameter(nG=          64)
	data (Tgibbs(j),j=1,          64) /
     &  0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 0.41000E+04, 
     &  0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 0.46000E+04, 
     &  0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 0.51000E+04, 
     &  0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 0.56000E+04, 
     &  0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          64) /
     &  0.11880E+03, 0.10664E+03, 0.10062E+03, 0.94859E+02, 0.94638E+02, 
     &  0.88695E+02, 0.82788E+02, 0.76916E+02, 0.71078E+02, 0.59500E+02, 
     &  0.48130E+02, 0.37878E+02, 0.27759E+02, 0.17758E+02, 0.78630E+01, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_ZnS_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          11),Ggibbs(          11)
	parameter(nG=          11)
	data (Tgibbs(j),j=1,          11) /
     &  0.30000E+03, 0.40000E+03, 0.50000E+03, 0.60000E+03, 0.70000E+03, 
     &  0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 0.12000E+04, 
     &  0.13000E+04 /
	data (Ggibbs(j),j=1,          11) /
     & -0.23930E+03,-0.22970E+03,-0.22020E+03,-0.21070E+03,-0.20130E+03, 
     & -0.19100E+03,-0.18070E+03,-0.17040E+03,-0.16020E+03,-0.14810E+03, 
     & -0.12820E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_MgO_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          51),Ggibbs(          51)
	parameter(nG=          51)
	data (Tgibbs(j),j=1,          51) /
     &  0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 0.40000E+03, 
     &  0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 
     &  0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 
     &  0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 
     &  0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 
     &  0.25000E+04, 0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 
     &  0.30000E+04, 0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 
     &  0.35000E+04, 0.36000E+04, 0.37000E+04, 0.38000E+04, 0.39000E+04, 
     &  0.40000E+04, 0.41000E+04, 0.42000E+04, 0.43000E+04, 0.44000E+04, 
     &  0.45000E+04, 0.46000E+04, 0.47000E+04, 0.48000E+04, 0.49000E+04, 
     &  0.50000E+04 /
	data (Ggibbs(j),j=1,          51) /
     & -0.58960E+03,-0.57949E+03,-0.56895E+03,-0.56875E+03,-0.55790E+03, 
     & -0.54708E+03,-0.53631E+03,-0.52560E+03,-0.51493E+03,-0.50429E+03, 
     & -0.49295E+03,-0.48140E+03,-0.46984E+03,-0.45829E+03,-0.44357E+03, 
     & -0.42275E+03,-0.40203E+03,-0.38139E+03,-0.36085E+03,-0.34039E+03, 
     & -0.32002E+03,-0.29973E+03,-0.27951E+03,-0.25937E+03,-0.23930E+03, 
     & -0.21930E+03,-0.19937E+03,-0.17951E+03,-0.15971E+03,-0.13998E+03, 
     & -0.12031E+03,-0.10070E+03,-0.83538E+02,-0.66584E+02,-0.49708E+02, 
     & -0.32905E+02,-0.16174E+02, 0.48700E+00, 0.17083E+02, 0.33615E+02, 
     &  0.50084E+02, 0.66495E+02, 0.82849E+02, 0.99149E+02, 0.11539E+03, 
     &  0.13159E+03, 0.14774E+03, 0.16384E+03, 0.17990E+03, 0.19591E+03, 
     &  0.21189E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_NH3_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          60),Ggibbs(          60)
	parameter(nG=          60)
	data (Tgibbs(j),j=1,          60) /
     &  0.10000E+03, 0.20000E+03, 0.30000E+03, 0.40000E+03, 0.50000E+03, 
     &  0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 
     &  0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 
     &  0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 
     &  0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 
     &  0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 
     &  0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 
     &  0.36000E+04, 0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 
     &  0.41000E+04, 0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 
     &  0.46000E+04, 0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 
     &  0.51000E+04, 0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 
     &  0.56000E+04, 0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          60) /
     & -0.50676E+02,-0.29793E+02,-0.10202E+02, 0.95258E+01, 0.29510E+02, 
     &  0.49710E+02, 0.70073E+02, 0.90553E+02, 0.11112E+03, 0.13174E+03, 
     &  0.15240E+03, 0.17308E+03, 0.19377E+03, 0.21447E+03, 0.23516E+03, 
     &  0.25584E+03, 0.27651E+03, 0.29716E+03, 0.31780E+03, 0.33842E+03, 
     &  0.35902E+03, 0.37960E+03, 0.40016E+03, 0.42070E+03, 0.44123E+03, 
     &  0.46173E+03, 0.48221E+03, 0.50268E+03, 0.52313E+03, 0.54356E+03, 
     &  0.56398E+03, 0.58438E+03, 0.60476E+03, 0.62513E+03, 0.64549E+03, 
     &  0.66583E+03, 0.68616E+03, 0.70649E+03, 0.72680E+03, 0.74710E+03, 
     &  0.76739E+03, 0.78767E+03, 0.80795E+03, 0.82822E+03, 0.84849E+03, 
     &  0.86875E+03, 0.88900E+03, 0.90925E+03, 0.92950E+03, 0.94975E+03, 
     &  0.96999E+03, 0.99024E+03, 0.10105E+04, 0.10307E+04, 0.10510E+04, 
     &  0.10712E+04, 0.10915E+04, 0.11117E+04, 0.11320E+04, 0.11522E+04 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_NH3_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          62),Ggibbs(          62)
	parameter(nG=          62)
	data (Tgibbs(j),j=1,          62) /
     &  0.00000E+00, 0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.40000E+03, 0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 
     &  0.90000E+03, 0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 
     &  0.14000E+04, 0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 
     &  0.19000E+04, 0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 
     &  0.24000E+04, 0.25000E+04, 0.26000E+04, 0.27000E+04, 0.28000E+04, 
     &  0.29000E+04, 0.30000E+04, 0.31000E+04, 0.32000E+04, 0.33000E+04, 
     &  0.34000E+04, 0.35000E+04, 0.36000E+04, 0.37000E+04, 0.38000E+04, 
     &  0.39000E+04, 0.40000E+04, 0.41000E+04, 0.42000E+04, 0.43000E+04, 
     &  0.44000E+04, 0.45000E+04, 0.46000E+04, 0.47000E+04, 0.48000E+04, 
     &  0.49000E+04, 0.50000E+04, 0.51000E+04, 0.52000E+04, 0.53000E+04, 
     &  0.54000E+04, 0.55000E+04, 0.56000E+04, 0.57000E+04, 0.58000E+04, 
     &  0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          62) /
     & -0.38907E+02,-0.34034E+02,-0.25679E+02,-0.16367E+02,-0.16183E+02, 
     & -0.59410E+01, 0.48000E+01, 0.15879E+02, 0.27190E+02, 0.38662E+02, 
     &  0.50247E+02, 0.61910E+02, 0.73625E+02, 0.85373E+02, 0.97141E+02, 
     &  0.10892E+03, 0.12070E+03, 0.13247E+03, 0.14423E+03, 0.15599E+03, 
     &  0.16772E+03, 0.17945E+03, 0.19115E+03, 0.20284E+03, 0.21451E+03, 
     &  0.22616E+03, 0.23779E+03, 0.24941E+03, 0.26100E+03, 0.27258E+03, 
     &  0.28414E+03, 0.29569E+03, 0.30722E+03, 0.31873E+03, 0.33023E+03, 
     &  0.34172E+03, 0.35319E+03, 0.36465E+03, 0.37610E+03, 0.38754E+03, 
     &  0.39897E+03, 0.41038E+03, 0.42180E+03, 0.43320E+03, 0.44459E+03, 
     &  0.45598E+03, 0.46736E+03, 0.47874E+03, 0.49012E+03, 0.50149E+03, 
     &  0.51286E+03, 0.52422E+03, 0.53559E+03, 0.54695E+03, 0.55832E+03, 
     &  0.56968E+03, 0.58104E+03, 0.59241E+03, 0.60378E+03, 0.61515E+03, 
     &  0.62652E+03, 0.63789E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_NH4Cl_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          17),Ggibbs(          17)
	parameter(nG=          17)
	data (Tgibbs(j),j=1,          17) /
     &  0.00000E+00, 0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.40000E+03, 0.45770E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 
     &  0.90000E+03, 0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 
     &  0.14000E+04, 0.15000E+04 /
	data (Ggibbs(j),j=1,          17) /
     & -0.31139E+03,-0.27639E+03,-0.23983E+03,-0.20309E+03,-0.20240E+03, 
     & -0.16505E+03, 0.50000E+03,-0.92110E+02,-0.56079E+02,-0.20332E+02, 
     &  0.15037E+02, 0.49958E+02, 0.84389E+02, 0.11831E+03, 0.15172E+03, 
     &  0.18463E+03, 0.21703E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Zn_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          25),Ggibbs(          25)
	parameter(nG=          25)
	data (Tgibbs(j),j=1,          25) /
     &  0.00000E+00, 0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 
     &  0.30000E+03, 0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 
     &  0.60000E+03, 0.69273E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 
     &  0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 
     &  0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04 /
	data (Ggibbs(j),j=1,          25) /
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.70000E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.19360E+01, 0.11646E+02, 0.21274E+02, 0.30827E+02, 
     &  0.40309E+02, 0.49725E+02, 0.59078E+02, 0.68373E+02, 0.77612E+02 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_FeO_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          49),Ggibbs(          49)
	parameter(nG=          49)
	data (Tgibbs(j),j=1,          49) /
     &  0.00000E+00, 0.30000E+03, 0.40000E+03, 0.50000E+03, 0.60000E+03, 
     &  0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 0.11000E+04, 
     &  0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 0.16000E+04, 
     &  0.16500E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 0.21000E+04, 
     &  0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 0.26000E+04, 
     &  0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 0.31000E+04, 
     &  0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 0.36000E+04, 
     &  0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 0.41000E+04, 
     &  0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 0.46000E+04, 
     &  0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04 /
	data (Ggibbs(j),j=1,          49) /
     &  0.60752E+02,-0.25130E+03,-0.24454E+03,-0.23802E+03,-0.23167E+03, 
     & -0.22543E+03,-0.21926E+03,-0.21312E+03,-0.20696E+03,-0.20067E+03, 
     & -0.19432E+03,-0.18795E+03,-0.18165E+03,-0.17541E+03,-0.16924E+03, 
     &  0.17000E+04,-0.15916E+03,-0.15383E+03,-0.14845E+03,-0.14309E+03, 
     & -0.13774E+03,-0.13241E+03,-0.12709E+03,-0.12178E+03,-0.11648E+03, 
     & -0.11119E+03,-0.10592E+03,-0.10065E+03,-0.95385E+02,-0.90131E+02, 
     & -0.77460E+02,-0.61133E+02,-0.44869E+02,-0.28665E+02,-0.12518E+02, 
     &  0.35750E+01, 0.19616E+02, 0.35608E+02, 0.51552E+02, 0.67453E+02, 
     &  0.83310E+02, 0.99128E+02, 0.11491E+03, 0.13065E+03, 0.14636E+03, 
     &  0.16203E+03, 0.17768E+03, 0.19330E+03, 0.20888E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Fe2O3_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          27),Ggibbs(          27)
	parameter(nG=          27)
	data (Tgibbs(j),j=1,          27) /
     &  0.00000E+00, 0.10000E+03, 0.20000E+03, 0.29815E+03, 0.30000E+03, 
     &  0.40000E+03, 0.50000E+03, 0.60000E+03, 0.70000E+03, 0.80000E+03, 
     &  0.90000E+03, 0.10000E+04, 0.11000E+04, 0.12000E+04, 0.13000E+04, 
     &  0.14000E+04, 0.15000E+04, 0.16000E+04, 0.17000E+04, 0.18000E+04, 
     &  0.19000E+04, 0.20000E+04, 0.21000E+04, 0.22000E+04, 0.23000E+04, 
     &  0.24000E+04, 0.25000E+04 /
	data (Ggibbs(j),j=1,          27) /
     & -0.81902E+03,-0.79744E+03,-0.77057E+03,-0.74352E+03,-0.74301E+03, 
     & -0.71573E+03,-0.68893E+03,-0.66264E+03,-0.63684E+03,-0.61148E+03, 
     & -0.58651E+03,-0.56187E+03,-0.53717E+03,-0.51237E+03,-0.48756E+03, 
     & -0.46289E+03,-0.43835E+03,-0.41391E+03,-0.38952E+03,-0.36512E+03, 
     & -0.33934E+03,-0.31341E+03,-0.28747E+03,-0.26153E+03,-0.23558E+03, 
     & -0.20962E+03,-0.18365E+03 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Mn_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          65),Ggibbs(          65)
	parameter(nG=          65)
	data (Tgibbs(j),j=1,          65) /
     &  0.00000E+00, 0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 
     &  0.30000E+03, 0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 
     &  0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 
     &  0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 
     &  0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 
     &  0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 
     &  0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 
     &  0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 
     &  0.36000E+04, 0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 
     &  0.41000E+04, 0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 
     &  0.46000E+04, 0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 
     &  0.51000E+04, 0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 
     &  0.56000E+04, 0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          65) /
     &  0.28205E+03, 0.26932E+03, 0.25500E+03, 0.24785E+03, 0.24101E+03, 
     &  0.24075E+03, 0.23369E+03, 0.22667E+03, 0.21971E+03, 0.21280E+03, 
     &  0.19912E+03, 0.18562E+03, 0.17230E+03, 0.15916E+03, 0.14624E+03, 
     &  0.13367E+03, 0.12125E+03, 0.10899E+03, 0.96927E+02, 0.85238E+02, 
     &  0.74375E+02, 0.63820E+02, 0.53414E+02, 0.43148E+02, 0.33015E+02, 
     &  0.23008E+02, 0.13121E+02, 0.33480E+01, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Mn_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          35),Ggibbs(          35)
	parameter(nG=          35)
	data (Tgibbs(j),j=1,          35) /
     &  0.00000E+00, 0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 
     &  0.30000E+03, 0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 
     &  0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 
     &  0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 
     &  0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 
     &  0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 
     &  0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04 /
	data (Ggibbs(j),j=1,          35) /
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.63150E+01, 0.15874E+02, 
     &  0.25334E+02, 0.34698E+02, 0.43970E+02, 0.53155E+02, 0.62256E+02 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Cr_s(T,G,Tmax)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          40),Ggibbs(          40)
	parameter(nG=          40)
	data (Tgibbs(j),j=1,          40) /
     &  0.00000E+00, 0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 
     &  0.30000E+03, 0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 
     &  0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 
     &  0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 
     &  0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 
     &  0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 
     &  0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 
     &  0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04 /
	data (Ggibbs(j),j=1,          40) /
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.55080E+01, 
     &  0.16982E+02, 0.28430E+02, 0.39855E+02, 0.51258E+02, 0.62641E+02 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
	subroutine Gibbs_Cr_g(T,G)
	IMPLICIT NONE
	integer nG,j,i
	real*8 T,G,Tmax
	real*8 Tgibbs(          65),Ggibbs(          65)
	parameter(nG=          65)
	data (Tgibbs(j),j=1,          65) /
     &  0.00000E+00, 0.10000E+03, 0.20000E+03, 0.25000E+03, 0.29815E+03, 
     &  0.30000E+03, 0.35000E+03, 0.40000E+03, 0.45000E+03, 0.50000E+03, 
     &  0.60000E+03, 0.70000E+03, 0.80000E+03, 0.90000E+03, 0.10000E+04, 
     &  0.11000E+04, 0.12000E+04, 0.13000E+04, 0.14000E+04, 0.15000E+04, 
     &  0.16000E+04, 0.17000E+04, 0.18000E+04, 0.19000E+04, 0.20000E+04, 
     &  0.21000E+04, 0.22000E+04, 0.23000E+04, 0.24000E+04, 0.25000E+04, 
     &  0.26000E+04, 0.27000E+04, 0.28000E+04, 0.29000E+04, 0.30000E+04, 
     &  0.31000E+04, 0.32000E+04, 0.33000E+04, 0.34000E+04, 0.35000E+04, 
     &  0.36000E+04, 0.37000E+04, 0.38000E+04, 0.39000E+04, 0.40000E+04, 
     &  0.41000E+04, 0.42000E+04, 0.43000E+04, 0.44000E+04, 0.45000E+04, 
     &  0.46000E+04, 0.47000E+04, 0.48000E+04, 0.49000E+04, 0.50000E+04, 
     &  0.51000E+04, 0.52000E+04, 0.53000E+04, 0.54000E+04, 0.55000E+04, 
     &  0.56000E+04, 0.57000E+04, 0.58000E+04, 0.59000E+04, 0.60000E+04 /
	data (Ggibbs(j),j=1,          65) /
     &  0.39534E+03, 0.38238E+03, 0.36737E+03, 0.35982E+03, 0.35255E+03, 
     &  0.35227E+03, 0.34475E+03, 0.33725E+03, 0.32979E+03, 0.32235E+03, 
     &  0.30755E+03, 0.29288E+03, 0.27831E+03, 0.26385E+03, 0.24950E+03, 
     &  0.23526E+03, 0.22114E+03, 0.20713E+03, 0.19325E+03, 0.17949E+03, 
     &  0.16586E+03, 0.15237E+03, 0.13901E+03, 0.12578E+03, 0.11269E+03, 
     &  0.99743E+02, 0.87589E+02, 0.75788E+02, 0.64044E+02, 0.52350E+02, 
     &  0.40702E+02, 0.29096E+02, 0.17528E+02, 0.59940E+01, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 
     &  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 /
	Tmax=Tgibbs(nG)
	if(T.lt.Tgibbs(1)) then
		G=Ggibbs(1)
	else if(T.gt.Tgibbs(nG)) then
		G=Ggibbs(nG)+(T-Tgibbs(nG))*(Ggibbs(nG-1)-Ggibbs(nG))/(Tgibbs(nG-1)-Tgibbs(nG))
	else
		do i=1,nG-1
			if(T.ge.Tgibbs(i).and.T.le.Tgibbs(i+1)) then
				G=Ggibbs(i)+(T-Tgibbs(i))*(Ggibbs(i+1)-Ggibbs(i))/(Tgibbs(i+1)-Tgibbs(i))
				return
			endif
		enddo
	endif
	return
	end
