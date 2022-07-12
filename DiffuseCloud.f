	subroutine DiffuseCloud(ii)
	use GlobalSetup
	use Constants
	use AtomsModule
	use CloudModule
	use TimingModule
	IMPLICIT NONE
	real*8,allocatable :: x(:),vsed(:),dx(:),vth(:),vthv(:)
	real*8,allocatable :: Sc(:),Sn(:),mpart(:),xMgO(:)
	real*8,allocatable :: An(:,:),y(:,:)
	real*8,allocatable :: Aomp(:,:),xomp(:),at_ab(:,:)
	real*8,allocatable :: drhoKd(:),drhovsed(:),tcinv(:,:),rho_av(:),densv(:,:),Kd(:)
	real*8 dz,z12,z13,z12_2,z13_2,g,rr,mutot,npart,tot,lambda,densv_t,tot1,tot2,tot3
	integer info,i,j,iter,NN,NRHS,niter,ii,k,ihaze
	real*8 cs,eps,frac_nuc,m_nuc,tcoaginv,Dp,vmol,f,mm,ComputeKzz,err,maxerr
	real*8 Pv,w_atoms(N_atoms),molfracs_atoms0(N_atoms),NKn,Kzz_r(nr)
	integer,allocatable :: IWORK(:),ixv(:,:),ixc(:,:),IWORKomp(:)
	real*8 sigmastar,Sigmadot,Pstar,gz,sigmamol,COabun,lmfp,fstick,kappa_cloud,fmin,rho_nuc
	logical ini,Tconverged
	character*500 cloudspecies(max(nclouds,1)),form
	real*8,allocatable :: CrV_prev0(:),CrT_prev0(:)
	logical,allocatable :: empty(:)
	integer iCS,ir,nrdo
	real*8 logP(nr),logx(nr),dlogx(nr)
	real*8,allocatable :: logCloudP(:)
	integer INCFD,IERR
	logical SKIP
	real*8 time
	integer itime

	call cpu_time(time)
	timecloud=timecloud-time
	call system_clock(itime)
	itimecloud=itimecloud-itime
	ctimecloud=ctimecloud+1

	nnr=(nr-1)*nr_cloud+1
	nCS=10
	if(.not.allocated(CloudP)) then
		allocate(CloudP(nnr))
		allocate(CloudT(nnr))
		allocate(CloudR(nnr))
		allocate(Clouddens(nnr))
		allocate(xv(nCS,nnr))
		allocate(xc(nCS,nnr))
		allocate(xn(nnr))
		allocate(xm(nnr))
		allocate(rpart(nnr))
	endif
	allocate(Kd(nnr))
	allocate(logCloudP(nnr))
	
	niter=500
	if(computeT) then
		if(nTiter.eq.1) then
			niter=20
		else if(nTiter.le.4) then
			niter=100
		endif
	endif

	w_atoms(1) = 1.00794		!'H'
	w_atoms(2) = 4.002602		!'He'
	w_atoms(3) = 12.011			!'C'
	w_atoms(4) = 14.00674		!'N'
	w_atoms(5) = 15.9994		!'O'
	w_atoms(6) = 22.989768		!'Na'
	w_atoms(7) = 24.3050		!'Mg'
	w_atoms(8) = 26.981539		!'Al'
	w_atoms(9) = 28.0855		!'Si'
	w_atoms(10) = 30.973762 	!'P'
	w_atoms(11) = 32.066 		!'S'
	w_atoms(12) = 35.4527 		!'Cl'
	w_atoms(13) = 39.0983 		!'K'
	w_atoms(14) = 40.078 		!'Ca'
	w_atoms(15) = 47.867 		!'Ti'
	w_atoms(16) = 50.9415 		!'V'
	w_atoms(17) = 55.845 		!'Fe'
	w_atoms(18) = 58.6934 		!'Ni'
	
	allocate(densv(nnr,nCS))

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
	i=0
c TiO2: 1
	i=i+1
	CSname(i)='TiO2'
	ATP(i)=77365.	! Al2O3 for now
	BTP(i)=39.3
	atoms_cloud(i,15)=1
	atoms_cloud(i,5)=2
	rhodust(i)=7.0d0	! moet ik nog checken
	CSnmol(i)=1d0
	ice(i)=.false.
c VO: 2
	i=i+1
	CSname(i)='VO'
	ATP(i)=77365.	! Al2O3 for now
	BTP(i)=39.3
	atoms_cloud(i,16)=1
	atoms_cloud(i,5)=1
	rhodust(i)=7.0d0	! moet ik nog checken
	CSnmol(i)=1d0
	ice(i)=.false.
c Al2O3: 3
	i=i+1
	CSname(i)='Al2O3'
	ATP(i)=77365.
	BTP(i)=39.3
	atoms_cloud(i,8)=2
	atoms_cloud(i,5)=3
	rhodust(i)=4.0d0
	CSnmol(i)=2d0
	ice(i)=.false.
c SiO2: 4
	i=i+1
	CSname(i)='SiO2'
	ATP(i)=69444.
	BTP(i)=33.1
	atoms_cloud(i,9)=1
	atoms_cloud(i,5)=2
	rhodust(i)=2.2d0
	CSnmol(i)=1d0
	ice(i)=.false.
c Silicates: 5
	i=i+1
	ATP(i)=68908.	! MgSiO3 for now
	BTP(i)=38.1

c ===========================================================
c = Old version, only using Mg to condense ==================
c ===========================================================
	atoms_cloud(i,1:N_atoms)=0
	atoms_cloud(i,7)=1
	atoms_cloud(i,5)=1
c ===========================================================
c = New version, using Na, K, Ca and Mg to condense =========
c ===========================================================
	atoms_cloud(i,1:N_atoms)=0
	atoms_cloud(i,6)=molfracs_atoms(6)
	atoms_cloud(i,7)=molfracs_atoms(7)
	atoms_cloud(i,13)=molfracs_atoms(13)
	atoms_cloud(i,14)=molfracs_atoms(14)
	tot=sum(atoms_cloud(i,1:N_atoms))
	atoms_cloud(i,1:N_atoms)=atoms_cloud(i,1:N_atoms)/tot
	atoms_cloud(i,5)=1
	write(CSname(i),'("Na",f3.1,"Mg",f3.1,"K",f3.1,"Ca",f3.1,"SiO",f3.1)') atoms_cloud(i,6),atoms_cloud(i,7),atoms_cloud(i,13)
     &				,atoms_cloud(i,14),atoms_cloud(i,5)+2
c ===========================================================
	rhodust(i)=2.0d0
	CSnmol(i)=3d0
	ice(i)=.false.
c H2O: 6
	i=i+1
	CSname(i)='H2O'
	ATP(i)=6511.
	BTP(i)=33.1
	atoms_cloud(i,1)=2
	atoms_cloud(i,5)=1
	rhodust(i)=1.0d0
	CSnmol(i)=1d0
	ice(i)=.true.
c Fe: 7
	i=i+1
	CSname(i)='Fe'
	ATP(i)=48354.
	BTP(i)=29.2
	atoms_cloud(i,17)=1
	rhodust(i)=7.8d0
	CSnmol(i)=1d0
	ice(i)=.false.
c FeS: 8
	i=i+1
	CSname(i)='FeS'
	ATP(i)=48354.
	BTP(i)=29.2
c	atoms_cloud(i,17)=1
	atoms_cloud(i,11)=1
	rhodust(i)=4.8d0
	maxT(i)=680d0
	CSnmol(i)=1d0
	ice(i)=.false.
c C: 9
	i=i+1
	CSname(i)='C'
	ATP(i)=93646.
	BTP(i)=36.7
	atoms_cloud(i,3)=1
	rhodust(i)=1.8d0
	CSnmol(i)=1d0
	ice(i)=.false.
c SiC: 10
	i=i+1
	CSname(i)='SiC'
	ATP(i)=78462.
	BTP(i)=37.8
	atoms_cloud(i,9)=1
c	atoms_cloud(i,3)=1
	rhodust(i)=2.5d0
	CSnmol(i)=1d0
	ice(i)=.false.

	nCS=i

	do i=1,nCS
		mu(i)=sum(w_atoms(1:N_atoms)*atoms_cloud(i,1:N_atoms))/CSnmol(i)
	enddo

	mutot=0d0
	xv_bot=1d200
	do i=1,N_atoms
		mutot=mutot+w_atoms(i)*molfracs_atoms(i)
	enddo
	do iCS=1,nCS
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
	enddo
	if(xv_bot(10).gt.xv_bot(9)) xv_bot(10)=xv_bot(9)
	if(xv_bot(8).gt.xv_bot(7)) xv_bot(8)=xv_bot(7)
	
	molfracs_atoms0=molfracs_atoms
	xv_bot=xv_bot*mu*CSnmol/mutot

	if(rainout) then
		densv(1,1:nCS)=(mu*mp/(kb*T(1)))*exp(BTP(1:nCS)-ATP(1:nCS)/T(1))
		do iCS=1,nCS
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
			ihaze=9
		case("THOLIN","tholin","Tholin")
			rho_nuc=1.00
			ihaze=9
		case("SiC")
			rho_nuc=3.22
			ihaze=10
		case("CARBON","Carbon","carbon")
			rho_nuc=1.80
			ihaze=9
		case("CORRUNDUM","Corrundum","corrundum","Al2O3")
			rho_nuc=3.97
			ihaze=3
		case("IRON","Iron","iron","Fe")
			rho_nuc=7.87
			ihaze=7
		case("SiO")
			rho_nuc=2.18
			ihaze=4
		case("TiO2")
			rho_nuc=4.23
			ihaze=1
		case("Enstatite","enstatite","ENSTATITE")
			rho_nuc=3.20
			ihaze=5
		case("MIX")
			tot=Cloud(ii)%fHazeSiO+Cloud(ii)%fHazeAl2O3+Cloud(ii)%fHazeTiO2+Cloud(ii)%fHazeTholin+Cloud(ii)%fHazeFe
			Cloud(ii)%fHazeSiO=Cloud(ii)%fHazeSiO/tot
			Cloud(ii)%fHazeAl2O3=Cloud(ii)%fHazeAl2O3/tot
			Cloud(ii)%fHazeTiO2=Cloud(ii)%fHazeTiO2/tot
			Cloud(ii)%fHazeTholin=Cloud(ii)%fHazeTholin/tot
			Cloud(ii)%fHazeFe=Cloud(ii)%fHazeFe/tot
			rho_nuc=1d0/(Cloud(ii)%fHazeSiO/2.18+Cloud(ii)%fHazeAl2O3/3.97+Cloud(ii)%fHazeTiO2/4.23+
     &					Cloud(ii)%fHazeTholin/1.00+Cloud(ii)%fHazeFe/7.87)
			ihaze=1
		case default
			call output("hazetype unknown")
			stop
	end select

	m_nuc=4d0*pi*r_nuc**3*rho_nuc/3d0

	allocate(mpart(nnr))
	allocate(rho_av(nnr))
	allocate(y(nnr,5))
	allocate(Sn(nnr))
	allocate(vth(nnr))
	allocate(drhoKd(nnr))
	allocate(drhovsed(nnr))
	allocate(tcinv(niter,nnr))
	allocate(vsed(nnr))

	allocate(ixv(nCS,nnr))
	allocate(ixc(nCS,nnr))

	allocate(xMgO(nnr))
	xMgO=0d0
	
	do iCS=1,nCS
		j=0
		do i=1,nnr
			j=j+1
			ixc(iCS,i)=j
		enddo
		do i=1,nnr
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
	call regridarray(logP,logx,nr,logCloudP,CloudT,nnr)
c	call DPCHIM(nr,logP,logx,dlogx,INCFD)
c	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudT,IERR)
	CloudT(1:nnr)=exp(CloudT(1:nnr))

	SKIP=.false.
	INCFD=1
	logx(1:nr)=log(dens(1:nr))
	call regridarray(logP,logx,nr,logCloudP,Clouddens,nnr)
c	call DPCHIM(nr,logP,logx,dlogx,INCFD)
c	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,Clouddens,IERR)
	Clouddens(1:nnr)=exp(Clouddens(1:nnr))

	SKIP=.false.
	INCFD=1
	logx(1:nr)=log(R(1:nr))
	call regridarray(logP,logx,nr,logCloudP,CloudR,nnr)
c	call DPCHIM(nr,logP,logx,dlogx,INCFD)
c	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudR,IERR)
	CloudR(1:nnr)=exp(CloudR(1:nnr))

	if((Kzz_deep.gt.0d0.and.Kzz_1bar.gt.0d0).or.Cloud(ii)%Kzz.le.0d0) then
		do i=1,nnr
			Kd(i)=ComputeKzz(CloudP(i))
		enddo
	else
		Kd=Cloud(ii)%Kzz
	endif

	f=0.1d0

	allocate(empty(nnr))
	NN=2*nnr
	allocate(x(NN))
	allocate(IWORK(10*NN*NN))

	allocate(An(nnr,nnr))

	rpart=r_nuc
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
	
	do i=1,nnr
		densv(i,1:nCS)=(mu*mp/(kb*CloudT(i)))*exp(BTP(1:nCS)-ATP(1:nCS)/CloudT(i))
		do iCS=1,nCS
			if(cloudT(i).gt.maxT(iCS)) densv(i,iCS)=densv(i,iCS)+
     &                    (mu(iCS)*mp/(kb*CloudT(i)*10d0))*exp(BTP(iCS)-ATP(iCS)/(CloudT(i)*10d0))
		enddo
		if(i.eq.1) then
			drhoKd(i)=(Kd(i+1)*Clouddens(i+1)-Kd(i)*Clouddens(i))/(CloudR(i+1)-CloudR(i))
		else if(i.eq.nnr) then
			drhoKd(i)=(Kd(i)*Clouddens(i)-Kd(i-1)*Clouddens(i-1))/(CloudR(i)-CloudR(i-1))
		else
			drhoKd(i)=(Kd(i+1)*Clouddens(i+1)-Kd(i-1)*Clouddens(i-1))/(CloudR(i+1)-CloudR(i-1))
		endif

		gz=Ggrav*Mplanet/CloudR(i)**2

		Sn(i)=(Clouddens(i)*gz*Sigmadot/(sigmastar*CloudP(i)*1d6*sqrt(2d0*pi)))*exp(-log(CloudP(i)/Pstar)**2/(2d0*sigmastar**2))
	enddo

c start the loop
	do iter=1,niter
	call tellertje(iter,niter)

	vsed=0d0
	do i=1,nnr
		cs=sqrt(kb*CloudT(i)/(2.3d0*mp))
		vth(i)=sqrt(8d0*kb*CloudT(i)/(pi*2.3d0*mp))
		vsed(i)=-sqrt(pi)*rpart(i)*rho_av(i)*Ggrav*Mplanet/(2d0*Clouddens(i)*vth(i)*CloudR(i)**2)
		mpart(i)=rho_av(i)*4d0*pi*rpart(i)**3/3d0
	enddo

	do i=1,nnr
		if(i.eq.nnr) then
			drhovsed(i)=(vsed(i)*Clouddens(i)-vsed(i-1)*Clouddens(i-1))/(CloudR(i)-CloudR(i-1))
		else
			drhovsed(i)=(vsed(i+1)*Clouddens(i+1)-vsed(i)*Clouddens(i))/(CloudR(i+1)-CloudR(i))
		endif
	enddo
	empty=.false.

c equations for number of Nuclii
	An=0d0
	x=0d0
	j=0
	do i=2,nnr-1
		j=j+1

		An(j,i)=An(j,i)-drhovsed(i)

		dz=CloudR(i+1)-CloudR(i)
		An(j,i+1)=An(j,i+1)+(-Clouddens(i)*vsed(i)+drhoKd(i))/dz
		An(j,i)=An(j,i)-(-Clouddens(i)*vsed(i)+drhoKd(i))/dz

		dz=CloudR(i+1)-CloudR(i-1)
		An(j,i+1)=An(j,i+1)+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i+1)-CloudR(i)))
		An(j,i-1)=An(j,i-1)+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i)-CloudR(i-1)))
		An(j,i)=An(j,i)-2d0*Clouddens(i)*Kd(i)*(1d0/(dz*(CloudR(i+1)-CloudR(i)))+1d0/(dz*(CloudR(i)-CloudR(i-1))))

		x(j)=-Sn(i)

c coagulation
		if(coagulation) then
			npart=xn(i)*Clouddens(i)
			lmfp=2.3d0*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
			vmol=0.5d0*lmfp*vth(i)
			Dp=kb*CloudT(i)/(6d0*pi*rpart(i)*vmol*Clouddens(i))

c			tcoaginv=npart*pi*(2d0*rpart(i))**2*abs(vsed(i))
c rewritten for better convergence
			tcoaginv=sqrt(pi)*3d0*(sum(xc(1:nCS,i))+xm(i))*Ggrav*Mplanet/(8d0*vth(i)*CloudR(i)**2)

			if(sqrt(16d0*kb*CloudT(i)/(pi*mpart(i)))*rpart(i).lt.Dp) then
				tcoaginv=tcoaginv+2d0*pi*rpart(i)**2*npart*sqrt(16d0*kb*CloudT(i)/(pi*mpart(i)))
			else
				tcoaginv=tcoaginv+2d0*pi*rpart(i)*npart*Dp
			endif
			tcoaginv=tcoaginv*exp(-(vsed(i)/vfrag)**2)

			if(.not.tcoaginv.gt.0d0) tcoaginv=0d0

			tcinv(iter,i)=tcoaginv
			
c			call computemedian(tcinv(1:iter,i),iter,tcoaginv)
			tcoaginv=sum(tcinv(1:iter,i))/real(iter)
c	call printstats(tcinv(1:iter,i),iter)

			An(j,i)=An(j,i)-Clouddens(i)*tcoaginv
		endif
	enddo
	i=nnr
	j=j+1
	dz=CloudR(i)-CloudR(i-1)
	An(j,i)=Kd(i)/dz-vsed(i)
	An(j,i-1)=-Kd(i)/dz
	x(j)=0d0
	i=1
	j=j+1
	An(j,i)=1d0
	x(j)=0d0

	NRHS=1
	call DGESV( nnr, NRHS, An, nnr, IWORK, x, nnr, info )
	
	xn(1:nnr)=x(1:nnr)

	if(coagulation.and.Cloud(ii)%haze) then
c equations for mass in Nuclii
		An=0d0
		x=0d0
		j=0
		do i=2,nnr-1
			j=j+1

			An(j,i)=An(j,i)-drhovsed(i)

			dz=CloudR(i+1)-CloudR(i)
			An(j,i+1)=An(j,i+1)+(-Clouddens(i)*vsed(i)+drhoKd(i))/dz
			An(j,i)=An(j,i)-(-Clouddens(i)*vsed(i)+drhoKd(i))/dz

			dz=CloudR(i+1)-CloudR(i-1)
			An(j,i+1)=An(j,i+1)+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i+1)-CloudR(i)))
			An(j,i-1)=An(j,i-1)+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i)-CloudR(i-1)))
			An(j,i)=An(j,i)-2d0*Clouddens(i)*Kd(i)*(1d0/(dz*(CloudR(i+1)-CloudR(i)))+1d0/(dz*(CloudR(i)-CloudR(i-1))))

			x(j)=-Sn(i)
		enddo
		i=nnr
		j=j+1
		dz=CloudR(i)-CloudR(i-1)
		An(j,i)=Kd(i)/dz-vsed(i)
		An(j,i-1)=-Kd(i)/dz
		x(j)=0d0
		i=1
		j=j+1
		An(j,i)=1d0
		x(j)=0d0

		NRHS=1
		call DGESV( nnr, NRHS, An, nnr, IWORK, x, nnr, info )
	
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

	if(Cloud(ii)%haze) then
		do i=1,nnr
			xm(i)=xm(i)*max(0d0,1d0-densv(i,ihaze)/(Clouddens(i)*xv_bot(ihaze)))
			if(xn(i).gt.xm(i)) then
				xn(i)=xm(i)
			endif
		enddo
	else
		xm=0d0
	endif
	xn(1:nnr)=xn(1:nnr)/m_nuc

	if(Cloud(ii)%condensates) then

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(Sc,vthv,cs,Aomp,xomp,IWORKomp,iCS,i,j,dz,NRHS,INFO)
!$OMP& SHARED(nCS,nnr,CloudT,Clouddens,CloudP,mu,fstick,CloudR,densv,drhovsed,Kd,xn,empty,
!$OMP&		NN,rpart,ixc,vsed,drhoKd,ixv,m_nuc,mpart,xv_bot,xc,xv,iter,nTiter,i3D)
	allocate(vthv(nnr))
	allocate(Sc(nnr))
	allocate(xomp(NN))
	allocate(IWORKomp(10*NN*NN))
	allocate(Aomp(NN,NN))
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC, 1)
	do iCS=1,nCS

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
	do i=2,nnr-1
		j=j+1

		if(.not.empty(i)) then
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-drhovsed(i)

		dz=CloudR(i+1)-CloudR(i)
		Aomp(j,ixc(iCS,i+1))=Aomp(j,ixc(iCS,i+1))+(-Clouddens(i)*vsed(i)+drhoKd(i))/dz
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-(-Clouddens(i)*vsed(i)+drhoKd(i))/dz

		dz=CloudR(i+1)-CloudR(i-1)
		Aomp(j,ixc(iCS,i+1))=Aomp(j,ixc(iCS,i+1))+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i+1)-CloudR(i)))
		Aomp(j,ixc(iCS,i-1))=Aomp(j,ixc(iCS,i-1))+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i)-CloudR(i-1)))
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-2d0*Clouddens(i)*Kd(i)*(1d0/(dz*(CloudR(i+1)-CloudR(i)))
     &					+1d0/(dz*(CloudR(i)-CloudR(i-1))))

		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))+Sc(i)*xn(i)*Clouddens(i)
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-Sc(i)*densv(i,iCS)/mpart(i)
		else
		Aomp(j,ixc(iCS,i))=1d0
		endif
		xomp(j)=0d0

		j=j+1

		dz=CloudR(i+1)-CloudR(i)
		Aomp(j,ixv(iCS,i+1))=Aomp(j,ixv(iCS,i+1))+(drhoKd(i))/dz
		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))-(drhoKd(i))/dz

		dz=CloudR(i+1)-CloudR(i-1)
		Aomp(j,ixv(iCS,i+1))=Aomp(j,ixv(iCS,i+1))+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i+1)-CloudR(i)))
		Aomp(j,ixv(iCS,i-1))=Aomp(j,ixv(iCS,i-1))+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i)-CloudR(i-1)))
		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))-2d0*Clouddens(i)*Kd(i)*(1d0/(dz*(CloudR(i+1)-CloudR(i)))
     &					+1d0/(dz*(CloudR(i)-CloudR(i-1))))

		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))-Sc(i)*xn(i)*Clouddens(i)
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+Sc(i)*densv(i,iCS)/mpart(i)
		xomp(j)=0d0
	enddo
	i=1
	j=j+1
	Aomp(j,ixc(iCS,i))=1d0
	xomp(j)=0d0

	j=j+1
	Aomp(j,ixv(iCS,i))=1d0
	xomp(j)=xv_bot(iCS)

	i=nnr
	dz=CloudR(i)-CloudR(i-1)
	j=j+1
	Aomp(j,ixc(iCS,i))=Kd(i)/dz-vsed(i)
	Aomp(j,ixc(iCS,i-1))=-Kd(i)/dz
	xomp(j)=0d0!-Mc_top/Clouddens(i)

	j=j+1
	Aomp(j,ixv(iCS,i))=Kd(i)/dz
	Aomp(j,ixv(iCS,i-1))=-Kd(i)/dz
	xomp(j)=0d0!Mc_top/Clouddens(i)

	NRHS=1
	info=0
	call DGESV( NN, NRHS, Aomp, NN, IWORKomp, xomp, NN, info )

	do i=1,NN
		if(.not.xomp(i).gt.0d0) xomp(i)=0d0
	enddo

	do i=1,nnr
		xc(iCS,i)=xomp(ixc(iCS,i))
		xv(iCS,i)=xomp(ixv(iCS,i))
		if(empty(i)) xc(iCS,i)=0d0
	enddo
	do i=1,nnr
		if(xc(iCS,i).lt.0d0) xc(iCS,i)=0d0
		if(xv(iCS,i).lt.0d0) xv(iCS,i)=0d0
	enddo

	enddo
!$OMP END DO
!$OMP FLUSH
	deallocate(vthv)
	deallocate(Sc)
	deallocate(xomp)
	deallocate(IWORKomp)
	deallocate(Aomp)
!$OMP END PARALLEL

	else
	do iCS=1,nCS
		do i=1,nnr
			xc(iCS,i)=0d0
			xv(iCS,i)=xv_bot(iCS)
		enddo
	enddo
	endif

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
			if(.not.rr.ge.r_nuc) then
				rr=r_nuc
				xn(i)=(3d0*(tot/(rr**3))/(4d0*pi*rho_av(i)))
			endif
			if(.not.rr.le.1d0) then
				rr=1d0
				xn(i)=(3d0*(tot/(rr**3))/(4d0*pi*rho_av(i)))
			endif
		else
			rr=r_nuc
			xn(i)=(3d0*(tot/(rr**3))/(4d0*pi*rho_av(i)))
		endif
		err=abs(rr-rpart(i))/(rr+rpart(i))
		if(err.gt.maxerr) maxerr=err
		rpart(i)=sqrt(rr*rpart(i))
	enddo
	if(maxerr.lt.1d-3) exit
	enddo
c end the loop


	do k=1,nnr
c correction for silicates
		f=mu(4)*CSnmol(4)+mu(5)*CSnmol(5)
		mm=(f/(mu(5)*CSnmol(5))-1d0)*xc(5,k)
		if(xc(4,k).lt.mm) then
			f=xc(4,k)/mm
			xc(4,k)=0d0
			xMgO(k)=xc(5,k)*(1d0-f)
			xc(5,k)=xc(5,k)*f*(mu(5)*CSnmol(5)+w_atoms(9)+2d0*w_atoms(5))/(mu(5)*CSnmol(5))
		else
			xc(4,k)=xc(4,k)-mm
			xc(5,k)=xc(5,k)+mm
			xMgO(k)=0d0
		endif
c correction for FeS
		f=mu(7)*CSnmol(7)+mu(8)*CSnmol(8)
		mm=(f/(mu(8)*CSnmol(8))-1d0)*xc(8,k)
		if(xc(7,k).lt.mm) then
			f=xc(7,k)/mm
			xc(7,k)=0d0
			xv(8,k)=xv(8,k)+xc(8,k)*(1d0-f)
			xc(8,k)=xc(8,k)*f*(mu(7)*CSnmol(7)+mu(8)*CSnmol(8))/(mu(8)*CSnmol(8))
		else
			xc(7,k)=xc(7,k)-mm
			xc(8,k)=xc(8,k)+mm
		endif
c correction for SiC
		f=mu(9)*CSnmol(9)+mu(10)*CSnmol(10)
		mm=(f/(mu(10)*CSnmol(10))-1d0)*xc(10,k)
		if(xc(9,k).lt.mm) then
			f=xc(9,k)/mm
			xc(9,k)=0d0
			xv(10,k)=xv(10,k)+xc(10,k)*(1d0-f)
			xc(10,k)=xc(10,k)*f*(mu(9)*CSnmol(9)+mu(10)*CSnmol(10))/(mu(10)*CSnmol(10))
		else
			xc(9,k)=xc(9,k)-mm
			xc(10,k)=xc(10,k)+mm
		endif
	enddo

	allocate(dx(nnr))
	logP(1:nr)=-log(P(1:nr))
	logCloudP(1:nnr)=-log(CloudP(1:nnr))

	SKIP=.false.
	INCFD=1
	x(1:nnr)=xm(1:nnr)*Clouddens(1:nnr)
	do iCS=1,nCS
		x(1:nnr)=x(1:nnr)+xc(iCS,1:nnr)*Clouddens(1:nnr)
	enddo
	call regridarray(logCloudP,x,nnr,logP,cloud_dens(1:nr,ii),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,cloud_dens(1:nr,ii),IERR)

	SKIP=.false.
	INCFD=1
	x(1:nnr)=rpart(1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%rv(1:nr),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%rv(1:nr),IERR)

	Cloud(ii)%frac(1:nr,1:20)=0d0
c TiO
	SKIP=.false.
	INCFD=1
	x(1:nnr)=xc(1,1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,1),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%frac(1:nr,1),IERR)
	Cloud(ii)%frac(1:nr,1)=Cloud(ii)%frac(1:nr,1)/3d0
	Cloud(ii)%frac(1:nr,2)=Cloud(ii)%frac(1:nr,1)
	Cloud(ii)%frac(1:nr,3)=Cloud(ii)%frac(1:nr,1)
c Al2O3
	SKIP=.false.
	INCFD=1
	x(1:nnr)=xc(3,1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,10),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%frac(1:nr,10),IERR)
c Silicates
	SKIP=.false.
	INCFD=1
	x(1:nnr)=xc(5,1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,13),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%frac(1:nr,13),IERR)
	Cloud(ii)%frac(1:nr,13)=Cloud(ii)%frac(1:nr,13)/3d0
	Cloud(ii)%frac(1:nr,14)=Cloud(ii)%frac(1:nr,13)
	Cloud(ii)%frac(1:nr,15)=Cloud(ii)%frac(1:nr,13)
c SiO2
	SKIP=.false.
	INCFD=1
	x(1:nnr)=xc(4,1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,8),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%frac(1:nr,8),IERR)
c H2O
	SKIP=.false.
	INCFD=1
	x(1:nnr)=xc(6,1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,18),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%frac(1:nr,18),IERR)
c Fe+FeS
	SKIP=.false.
	INCFD=1
	x(1:nnr)=(xc(7,1:nnr)+xc(8,1:nnr))
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,9),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%frac(1:nr,9),IERR)
c SiC
	SKIP=.false.
	INCFD=1
	x(1:nnr)=xc(10,1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,17),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%frac(1:nr,17),IERR)
c C
	SKIP=.false.
	INCFD=1
	x(1:nnr)=xc(9,1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,16),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%frac(1:nr,16),IERR)
c MgO
	SKIP=.false.
	INCFD=1
	x(1:nnr)=xMgO(1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,12),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%frac(1:nr,12),IERR)
c Seed particles
	SKIP=.false.
	INCFD=1
	x(1:nnr)=xm(1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,19),nr)
c	call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c	call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,Cloud(ii)%frac(1:nr,19),IERR)

	do i=1,nr
		tot=sum(Cloud(ii)%frac(i,1:19))
		if(tot.gt.0d0) then
			Cloud(ii)%frac(i,1:19)=Cloud(ii)%frac(i,1:19)/tot
		else
			Cloud(ii)%frac(i,1:19)=1d0/19d0
			cloud_dens(i,ii)=0d0
		endif
	enddo

	if(.not.retrieval) then
		open(unit=20,file=trim(outputdir) // '/cloudstructure.dat',RECL=1000)
		form='("#",a18,a19,a19,' // trim(int2string(nCS+1,'(i4)')) // 'a23,a19,a19,a19)'
		write(20,form) "P[bar]","dens[g/cm^3]","xn",(trim(CSname(i)),i=1,nCS),"MgO","r[micron]","T[K]","Jstar"
		form='(es19.7E3,es19.7E3,es19.7E3,' // trim(int2string(nCS+1,'(i4)')) // 'es23.7E3,es19.7E3,es19.7E3,es19.7E3)'
		do i=1,nnr
			write(20,form) CloudP(i),Clouddens(i),xn(i)*m_nuc,xc(1:nCS,i),xMgO(i),rpart(i),CloudT(i),Sn(i)/m_nuc
		enddo
		close(unit=20)
	endif

c Elemental abundances
	allocate(at_ab(nr,N_atoms))
	at_ab=0d0
	do iCS=1,nCS
		SKIP=.false.
		INCFD=1
		x(1:nnr)=xv(iCS,1:nnr)
		call regridarray(logCloudP,x,nnr,logP,logx,nr)
c		call DPCHIM(nnr,logCloudP,x,dx,INCFD)
c		call DPCHFE(nnr,logCloudP,x,dx,INCFD,SKIP,nr,logP,logx,IERR)
		do j=1,N_atoms
			at_ab(1:nr,j)=at_ab(1:nr,j)+logx(1:nr)*mutot*atoms_cloud(iCS,j)/(mu(iCS)*CSnmol(iCS))
		enddo
	enddo

c	open(unit=20,file=trim(outputdir) // '/atoms.dat',RECL=6000)
	ini=.true.
	do i=1,nr
		call tellertje(i,nr)
		molfracs_atoms(1:N_atoms)=at_ab(i,1:N_atoms)
		molfracs_atoms=molfracs_atoms+molfracs_atoms0
		molfracs_atoms(3)=molfracs_atoms(3)+COabun
		molfracs_atoms(5)=molfracs_atoms(5)+COabun
		do j=1,N_atoms
			if(.not.molfracs_atoms(j).gt.1d-50) then
				molfracs_atoms(j)=1d-50
			endif
		enddo
		tot=sum(molfracs_atoms(1:N_atoms))
		molfracs_atoms=molfracs_atoms/tot
		if(P(i).ge.mixP.or.i.eq.1) then
			call call_chemistry(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,.false.,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol)
   		else
   			mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
   			XeqCloud(i,1:nclouds)=XeqCloud(i-1,1:nclouds)
   			nabla_ad(i)=nabla_ad(i-1)
   			MMW(i)=MMW(i-1)
   			didcondens(i)=didcondens(i-1)
   		endif
c		write(20,*) P(i),molfracs_atoms(1:N_atoms)
	enddo
c	close(unit=20)
	deallocate(at_ab)

	if(disequilibrium) then
c       call disequilibrium code
c       input: 	R(1:nr+1) : These are the radial boundaries of the layers (bottom to top)
c       P(1:nr),T(1:nr) : These are the pressure and temperature inside the layers
c       molname(1:nmol) : names of the molecules included
c       Kzz_r(1:nr) : Diffusion coefficient
c       input/output:	mixrat_r(1:nr,1:nmol) : number densities inside each layer. Now set to equilibrium abundances.
	   call output("==================================================================")
	   call output("Computing disequilibrium chemistry")
		do i=1,nr
			Kzz_r(i)=ComputeKzz(P(i))
		enddo
	   call diseq_calc(nr,R(1:nr+1),P(1:nr),T(1:nr),nmol,molname(1:nmol),mixrat_r(1:nr, 1:nmol),COratio,Kzz_r(1:nr))
	endif
	do i=1,nr
		do j=1,nmol
			if(.not.mixrat_r(i,j).gt.0d0) mixrat_r(i,j)=0d0
		enddo
	enddo

	deallocate(densv)
	deallocate(mpart)
	deallocate(rho_av)
	deallocate(y)
	deallocate(Sn)
	deallocate(vth)
	deallocate(drhoKd)
	deallocate(drhovsed)
	deallocate(tcinv)
	deallocate(vsed)
	deallocate(ixv)
	deallocate(ixc)
	deallocate(x,dx)
	deallocate(empty)
	deallocate(IWORK)
	deallocate(An)
	deallocate(logCloudP)
	deallocate(Kd)

	call cpu_time(time)
	timecloud=timecloud+time
	call system_clock(itime)
	itimecloud=itimecloud+itime

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
	
	do ir=1,nr
		Tmini=3d0
		Tmaxi=3d5
		dens0=xv_bot*dens(ir)
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
		Tevap(ir,1:nCS)=T0(1:nCS)
		write(36,form) P(ir),T0(1:nCS)
	enddo

	close(unit=36)

	return
	end






