	subroutine DiffuseCloud(ii)
	use GlobalSetup
	use Constants
	use AtomsModule
	use CloudModule
	use modComputeT
	IMPLICIT NONE
	real*8,allocatable :: x(:),vsed(:),xtot(:),vth(:),vthv(:)
	real*8,allocatable :: Sc(:),Sn(:),rpart(:),mpart(:),xMgO(:)
	real*8,allocatable :: An(:,:),y(:,:),xv(:,:),xn(:),xc(:,:)
	real*8,allocatable :: Aomp(:,:),xomp(:)
	real*8,allocatable :: drho(:),drhovsed(:),tcinv(:),rho_av(:),densv(:),Kd(:)
	real*8 dz,z12,z13,z12_2,z13_2,g,rr,mutot,npart,tot,lambda
	integer info,i,j,iter,NN,NRHS,niter,ii,k
	real*8 cs,err,maxerr,eps,frac_nuc,m_nuc,tcoaginv,Dp,vmol,f,T0(nr),mm
	real*8 af,bf,f1,f2,Pv,w_atoms(N_atoms),molfracs_atoms0(N_atoms),NKn
	integer,allocatable :: IWORK(:),ixv(:,:),ixc(:,:),IWORKomp(:)
	real*8 sigmastar,Sigmadot,Pstar,gz,sigmamol,COabun,lmfp,fstick,kappa_cloud,fmin
	logical quadratic,ini,Tconverged,cloudsform
	character*500 cloudspecies(max(nclouds,1)),form
	real*8,allocatable :: CrV_prev0(:),CrT_prev0(:)

	real*8 Mc_top,Mn_top,IDP_dens,IDP_rad,fact,NOratio,SiOratio,Zout
	integer iCS

	nnr=(nr-1)*nr_cloud+1
	if(.not.allocated(CloudP)) then
		allocate(CloudP(nnr))
		allocate(CloudT(nnr))
		allocate(CloudR(nnr))
		allocate(Clouddens(nnr))
	endif
	allocate(Kd(nnr))

	T0=T

	niter=20

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
	
	nCS=10
	allocate(densv(nCS))

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
	endif

	quadratic=.false.

	call SetAbun

	COabun=min(molfracs_atoms(3),molfracs_atoms(5))
	molfracs_atoms(3)=molfracs_atoms(3)-COabun
	molfracs_atoms(5)=molfracs_atoms(5)-COabun

	ATP=28030d0		!6.89d4
	BTP=12.471d0		!37.8

	maxT=1d200

	atoms_cloud=0
	i=0
c TiO2
	i=i+1
	CSname(i)='TiO2'
	ATP(i)=77365.	! Al2O3 for now
	BTP(i)=39.3
	atoms_cloud(i,15)=1
	atoms_cloud(i,5)=2
	rhodust(i)=7.0d0	! moet ik nog checken
	CSnmol(i)=1d0
	ice(i)=.false.
c VO
	i=i+1
	CSname(i)='VO'
	ATP(i)=77365.	! Al2O3 for now
	BTP(i)=39.3
	atoms_cloud(i,16)=1
	atoms_cloud(i,5)=1
	rhodust(i)=7.0d0	! moet ik nog checken
	CSnmol(i)=1d0
	ice(i)=.false.
c Al2O3
	i=i+1
	CSname(i)='Al2O3'
	ATP(i)=77365.
	BTP(i)=39.3
	atoms_cloud(i,8)=2
	atoms_cloud(i,5)=3
	rhodust(i)=4.0d0
	CSnmol(i)=2d0
	ice(i)=.false.
c SiO2
	i=i+1
	CSname(i)='SiO2'
	ATP(i)=69444.
	BTP(i)=33.1
	atoms_cloud(i,9)=1
	atoms_cloud(i,5)=2
	rhodust(i)=2.2d0
	CSnmol(i)=1d0
	ice(i)=.false.
c Silicates
	i=i+1
	ATP(i)=68908.	! MgSiO3 for now
	BTP(i)=38.1
	atoms_cloud(i,9)=1
c	atoms_cloud(i,6)=min(molfracs_atoms(6),molfracs_atoms(8))/molfracs_atoms(9)
c	atoms_cloud(i,8)=min(molfracs_atoms(6),molfracs_atoms(8))/molfracs_atoms(9)
	atoms_cloud(i,7)=molfracs_atoms(7)/molfracs_atoms(9)
c	atoms_cloud(i,13)=molfracs_atoms(13)/molfracs_atoms(9)
c	atoms_cloud(i,14)=molfracs_atoms(14)/molfracs_atoms(9)
	atoms_cloud(i,5)=atoms_cloud(i,6)+atoms_cloud(i,7)+atoms_cloud(i,8)+atoms_cloud(i,14)+atoms_cloud(i,13)+2d0
c	write(CSname(i),'("Al",f3.1,"Na",f3.1,"Mg",f3.1,"SiO",f3.1)') atoms_cloud(i,8),atoms_cloud(i,6),atoms_cloud(i,7)
c     &				,atoms_cloud(i,5)
	write(CSname(i),'("Mg",f3.1,"SiO",f3.1)') atoms_cloud(i,7),atoms_cloud(i,5)
	atoms_cloud(i,9)=0
	atoms_cloud(i,5)=atoms_cloud(i,5)-2
	rhodust(i)=2.0d0
	CSnmol(i)=3d0
	ice(i)=.false.
c H2O
	i=i+1
	CSname(i)='H2O'
	ATP(i)=6511.
	BTP(i)=33.1
	atoms_cloud(i,1)=2
	atoms_cloud(i,5)=1
	rhodust(i)=1.0d0
	CSnmol(i)=1d0
	ice(i)=.true.
c Fe
	i=i+1
	CSname(i)='Fe'
	ATP(i)=48354.
	BTP(i)=29.2
	atoms_cloud(i,17)=1
	rhodust(i)=7.8d0
	CSnmol(i)=1d0
	ice(i)=.false.
c FeS
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
c C
	i=i+1
	CSname(i)='C'
	ATP(i)=93646.
	BTP(i)=36.7
	atoms_cloud(i,3)=1
	rhodust(i)=1.8d0
	CSnmol(i)=1d0
	ice(i)=.false.
c SiC
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
			if(atoms_cloud(iCS,k).gt.0) then
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
	
	molfracs_atoms0=molfracs_atoms
	xv_bot=xv_bot*mu*CSnmol/mutot

	if(Tform.le.0d0.and..not.trend_compute) then
		do iCS=1,nCS
			if(ice(iCS)) then
				xv_bot(iCS)=xv_bot(iCS)*f_wet
			else
				xv_bot(iCS)=xv_bot(iCS)*f_dry
			endif
		enddo
	endif

	if(rainout) then
		densv=(mu*mp/(kb*T(1)))*exp(BTP-ATP/T(1))
		do iCS=1,nCS
			if(T(1).gt.maxT(iCS)) densv(iCS)=densv(iCS)+(mu(iCS)*mp/(kb*T(1)*10d0))*exp(BTP(iCS)-ATP(iCS)/(T(1)*10d0))
		enddo
		do iCS=1,nCS
			if(densv(iCS).lt.dens(1)*xv_bot(iCS)) then
				call output("Rainout of: " // trim(CSname(iCS)) // "(" //
     &		trim(dbl2string(100d0*(1d0-densv(iCS)/(dens(1)*xv_bot(iCS))),'(f5.1)')) // " %)")
				xv_bot(iCS)=densv(iCS)/dens(1)
			endif
		enddo
	endif


	cloudsform=.false.
	do i=1,nr
		densv=(mu*mp/(kb*T(i)))*exp(BTP-ATP/T(i))
		do iCS=1,nCS
			if(densv(iCS).lt.dens(i)*xv_bot(iCS)) cloudsform=.true.
		enddo
	enddo
	if(.not.cloudsform) then!.or.(computeT.and.nTiter.lt.min(2,maxiter/2))) then
		cloud_dens=0d0
		do i=1,nr
			Cloud(ii)%rv(i)=r_nuc
			Cloud(ii)%frac(i,1:18)=1d0/18d0
		enddo
		return
	endif

	IDP_dens=0d0
	IDP_rad=0.01d0*micron

	Mc_top=-IDP_dens*sqrt(Ggrav*Mplanet/Dplanet)
	Mn_top=-Mc_top*(r_nuc/IDP_rad)**3

	sigmastar=0.1
	Pstar=60d-6

	Pstar=Cloud(ii)%P
	sigmastar=log(Cloud(ii)%dP)

	fstick=1d0
	
	sigmamol=8d-15

	eps=1d-3

	Kd=Cloud(ii)%Kzz*Clouddens**Cloud(ii)%Kzz_pow
	Sigmadot=Cloud(ii)%Sigmadot
	
	m_nuc=4d0*pi*r_nuc**3*6d0/3d0

	allocate(rpart(nnr))
	allocate(mpart(nnr))
	allocate(rho_av(nnr))
	allocate(y(nnr,5))
	allocate(Sn(nnr))
	allocate(vth(nnr))
	allocate(drho(nnr))
	allocate(drhovsed(nnr))
	allocate(xv(nCS,nnr))
	allocate(xc(nCS,nnr))
	allocate(xn(nnr))
	allocate(tcinv(nnr))
	allocate(vsed(nnr))

	allocate(ixv(nCS,nnr))
	allocate(ixc(nCS,nnr))

	allocate(xMgO(nnr))
	
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

	k=1
	do i=1,nr-1
		do j=1,nr_cloud
			CloudP(k)=10d0**(log10(P(i))+log10(P(i+1)/P(i))*real(j-1)/real(nr_cloud))
			CloudT(k)=10d0**(log10(T(i))+log10(T(i+1)/T(i))*real(j-1)/real(nr_cloud))
			CloudR(k)=10d0**(log10(R(i))+log10(R(i+1)/R(i))*real(j-1)/real(nr_cloud))
			Clouddens(k)=10d0**(log10(dens(i))+log10(dens(i+1)/dens(i))*real(j-1)/real(nr_cloud))
			k=k+1
		enddo
	enddo
	CloudP(k)=P(nr)
	CloudT(k)=T(nr)
	CloudR(k)=R(nr)
	Clouddens(k)=dens(nr)

	f=0.1d0
	rpart=r_nuc

	NN=2*nnr
	allocate(x(NN))
	allocate(IWORK(10*NN*NN))

	allocate(An(nnr,nnr))

	xn=0d0
	
	xv=0d0
	xc=0d0
	do j=1,nCS
		do i=1,nnr
			xv(j,i)=xv_bot(j)
		enddo
	enddo
	tcinv=0d0

c start the loop
	do iter=1,niter
	call tellertje(iter,niter)

	vsed=0d0
	do i=1,nnr
		cs=sqrt(kb*CloudT(i)/(2.3*mp))
		vth(i)=sqrt(8d0*kb*CloudT(i)/(pi*2.3*mp))
		vsed(i)=-sqrt(pi)*rpart(i)*rho_av(i)*Ggrav*Mplanet/(2d0*Clouddens(i)*vth(i)*CloudR(i)**2)
		mpart(i)=rho_av(i)*4d0*pi*rpart(i)**3/3d0
	enddo

	do i=1,nnr
		if(i.eq.1) then
			drho(i)=(Clouddens(i+1)-Clouddens(i))/(CloudR(i+1)-CloudR(i))
			drhovsed(i)=(vsed(i+1)*Clouddens(i+1)-vsed(i)*Clouddens(i))/(CloudR(i+1)-CloudR(i))
		else if(i.eq.nnr) then
			drho(i)=(Clouddens(i)-Clouddens(i-1))/(CloudR(i)-CloudR(i-1))
			drhovsed(i)=(vsed(i)*Clouddens(i)-vsed(i-1)*Clouddens(i-1))/(CloudR(i)-CloudR(i-1))
		else
			if(quadratic) then
			f1=(CloudR(i-1)**2-CloudR(i)**2)/(CloudR(i+1)**2-CloudR(i)**2)
			f2=(CloudR(i-1)-CloudR(i))/(CloudR(i+1)-CloudR(i))
			af=((Clouddens(i-1)-Clouddens(i))-f2*(Clouddens(i+1)-Clouddens(i)))/((CloudR(i-1)**2-CloudR(i)**2)-
     &					f2*(CloudR(i+1)**2-CloudR(i)**2))
			bf=((Clouddens(i-1)-Clouddens(i))-f1*(Clouddens(i+1)-Clouddens(i)))/((CloudR(i-1)-CloudR(i))-
     &					f1*(CloudR(i+1)-CloudR(i)))
			drho(i)=2d0*af*CloudR(i)+bf
			af=((vsed(i-1)*Clouddens(i-1)-vsed(i)*Clouddens(i))-f2*(vsed(i+1)*Clouddens(i+1)-vsed(i)*Clouddens(i)))/
     &				((CloudR(i-1)**2-CloudR(i)**2)-f2*(CloudR(i+1)**2-CloudR(i)**2))
			bf=((vsed(i-1)*Clouddens(i-1)-vsed(i)*Clouddens(i))-f1*(vsed(i+1)*Clouddens(i+1)-vsed(i)*Clouddens(i)))/
     &				((CloudR(i-1)-CloudR(i))-f1*(CloudR(i+1)-CloudR(i)))
			drhovsed(i)=2d0*af*CloudR(i)+bf
			else
			drho(i)=(Clouddens(i+1)-Clouddens(i-1))/(CloudR(i+1)-CloudR(i-1))
			drhovsed(i)=(vsed(i+1)*Clouddens(i+1)-vsed(i-1)*Clouddens(i-1))/(CloudR(i+1)-CloudR(i-1))
			endif
		endif

		densv=(mu*mp/(kb*CloudT(i)))*exp(BTP-ATP/CloudT(i))
		gz=Ggrav*Mplanet/CloudR(i)**2

		Sn(i)=(Clouddens(i)*gz*Sigmadot/(sigmastar*CloudP(i)*1d6*sqrt(2d0*pi)))*exp(-log(CloudP(i)/Pstar)**2/(2d0*sigmastar**2))
	enddo
	

c equations for Nuclii
	An=0d0
	x=0d0
	do i=2,nnr-1
		j=i+1

		dz=CloudR(i+1)-CloudR(i-1)

		An(j,i)=An(j,i)-drhovsed(i)

		if(quadratic) then
		f1=(CloudR(i-1)**2-CloudR(i)**2)/(CloudR(i+1)**2-CloudR(i)**2)
		f2=(CloudR(i-1)-CloudR(i))/(CloudR(i+1)-CloudR(i))
		af=1d0/((CloudR(i-1)**2-CloudR(i)**2)-f2*(CloudR(i+1)**2-CloudR(i)**2))
		bf=1d0/((CloudR(i-1)-CloudR(i))-f1*(CloudR(i+1)-CloudR(i)))

		An(j,i-1)=An(j,i-1)+(2d0*af*CloudR(i)+bf)*(-Clouddens(i)*vsed(i)+Kd(i)*drho(i))
		An(j,i+1)=An(j,i+1)-(2d0*af*f2*CloudR(i)+bf*f1)*(-Clouddens(i)*vsed(i)+Kd(i)*drho(i))
		An(j,i)=An(j,i)+(2d0*af*(f2-1d0)*CloudR(i)+bf*(f1-1d0))*(-Clouddens(i)*vsed(i)+Kd(i)*drho(i))

		An(j,i-1)=An(j,i-1)+2d0*af*Clouddens(i)*Kd(i)
		An(j,i+1)=An(j,i+1)-2d0*f2*af*Clouddens(i)*Kd(i)
		An(j,i)=An(j,i)+2d0*(f2-1d0)*af*Clouddens(i)*Kd(i)
		else
		An(j,i+1)=An(j,i+1)+(-Clouddens(i)*vsed(i)+Kd(i)*drho(i))/dz
		An(j,i-1)=An(j,i-1)-(-Clouddens(i)*vsed(i)+Kd(i)*drho(i))/dz

		An(j,i+1)=An(j,i+1)+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i+1)-CloudR(i)))
		An(j,i-1)=An(j,i-1)+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i)-CloudR(i-1)))
		An(j,i)=An(j,i)-2d0*Clouddens(i)*Kd(i)*(1d0/(dz*(CloudR(i+1)-CloudR(i)))+1d0/(dz*(CloudR(i)-CloudR(i-1))))
		endif

		x(j)=-Sn(i)

c coagulation
		if(coagulation) then
			npart=xn(i)*Clouddens(i)/m_nuc
			lmfp=2.3*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
			vmol=0.5d0*lmfp*vth(i)
			Dp=kb*CloudT(i)/(6d0*pi*rpart(i)*vmol*Clouddens(i))

c			tcoaginv=npart*pi*(2d0*rpart(i))**2*abs(vsed(i))
c rewritten for better convergence
			tcoaginv=sqrt(pi)*3d0*sum(xc(1:nCS,i))*Ggrav*Mplanet/(8d0*vth(i)*CloudR(i)**2)

			if(sqrt(16d0*kb*CloudT(i)/(pi*mpart(i)))*rpart(i).lt.Dp) then
				tcoaginv=tcoaginv+2d0*pi*rpart(i)**2*npart*sqrt(16d0*kb*CloudT(i)/(pi*mpart(i)))
			else
				tcoaginv=tcoaginv+2d0*pi*rpart(i)*npart*Dp
			endif

			if(.not.tcoaginv.gt.0d0) tcoaginv=0d0

			tcinv(i)=(tcinv(i)+tcoaginv)/2d0

			An(j,i)=An(j,i)-Clouddens(i)*tcinv(i)
		endif
	enddo
	i=nnr
	dz=CloudR(i)-CloudR(i-1)
	An(1,i)=Kd(i)/dz-vsed(i)
	An(1,i-1)=-Kd(i)/dz
	x(1)=0d0!-Mn_top/Clouddens(i)
	i=1
	An(2,i)=1d0
	x(2)=0d0

	NRHS=1
	call DGESV( nnr, NRHS, An, nnr, IWORK, x, nnr, info )
	
	do i=1,nnr
		if(.not.x(i).gt.0d0) x(i)=0d0
	enddo
	xn(1:nnr)=x(1:nnr)

	do i=1,nnr
		if(xn(i).lt.0d0) xn(i)=0d0
	enddo




!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(Sc,vthv,cs,Aomp,xomp,IWORKomp,iCS,i,j,dz,f1,f2,af,bf,NRHS,INFO)
!$OMP& SHARED(nCS,nnr,CloudT,Clouddens,CloudP,mu,fstick,CloudR,densv,ATP,BTP,maxT,drhovsed,Kd,xn,
!$OMP&		NN,rpart,ixc,quadratic,vsed,drho,ixv,m_nuc,mpart,xv_bot,xc,xv)
	allocate(vthv(nnr))
	allocate(Sc(nnr))
	allocate(xomp(NN))
	allocate(IWORKomp(10*NN*NN))
	allocate(Aomp(NN,NN))
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC, 1)
	do iCS=1,nCS

	do i=1,nnr
		cs=sqrt(kb*CloudT(i)/(2.3*mp))
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
		dz=CloudR(i+1)-CloudR(i-1)
		densv(iCS)=(mu(iCS)*mp/(kb*CloudT(i)))*exp(BTP(iCS)-ATP(iCS)/CloudT(i))
		if(cloudT(i).gt.maxT(iCS)) densv(iCS)=densv(iCS)+(mu(iCS)*mp/(kb*CloudT(i)*10d0))*exp(BTP(iCS)-ATP(iCS)/(CloudT(i)*10d0))

		f1=(CloudR(i-1)**2-CloudR(i)**2)/(CloudR(i+1)**2-CloudR(i)**2)
		f2=(CloudR(i-1)-CloudR(i))/(CloudR(i+1)-CloudR(i))
		af=1d0/((CloudR(i-1)**2-CloudR(i)**2)-f2*(CloudR(i+1)**2-CloudR(i)**2))
		bf=1d0/((CloudR(i-1)-CloudR(i))-f1*(CloudR(i+1)-CloudR(i)))

		j=j+1

		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-drhovsed(i)

		if(quadratic) then
		Aomp(j,ixc(iCS,i-1))=Aomp(j,ixc(iCS,i-1))+(2d0*af*CloudR(i)+bf)*(-Clouddens(i)*vsed(i)+Kd(i)*drho(i))
		Aomp(j,ixc(iCS,i+1))=Aomp(j,ixc(iCS,i+1))-(2d0*af*f2*CloudR(i)+bf*f1)*(-Clouddens(i)*vsed(i)+Kd(i)*drho(i))
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+(2d0*af*(f2-1d0)*CloudR(i)+bf*(f1-1d0))*(-Clouddens(i)*vsed(i)+Kd(i)*drho(i))

		Aomp(j,ixc(iCS,i-1))=Aomp(j,ixc(iCS,i-1))+2d0*af*Clouddens(i)*Kd(i)
		Aomp(j,ixc(iCS,i+1))=Aomp(j,ixc(iCS,i+1))-2d0*f2*af*Clouddens(i)*Kd(i)
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+2d0*(f2-1d0)*af*Clouddens(i)*Kd(i)
		else
		Aomp(j,ixc(iCS,i+1))=Aomp(j,ixc(iCS,i+1))+(-Clouddens(i)*vsed(i)+Kd(i)*drho(i))/dz
		Aomp(j,ixc(iCS,i-1))=Aomp(j,ixc(iCS,i-1))-(-Clouddens(i)*vsed(i)+Kd(i)*drho(i))/dz

		Aomp(j,ixc(iCS,i+1))=Aomp(j,ixc(iCS,i+1))+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i+1)-CloudR(i)))
		Aomp(j,ixc(iCS,i-1))=Aomp(j,ixc(iCS,i-1))+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i)-CloudR(i-1)))
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-2d0*Clouddens(i)*Kd(i)*(1d0/(dz*(CloudR(i+1)-CloudR(i)))
     &					+1d0/(dz*(CloudR(i)-CloudR(i-1))))
		endif

		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))+Sc(i)*xn(i)*Clouddens(i)/m_nuc
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))-Sc(i)*densv(iCS)/mpart(i)
		xomp(j)=0d0

		j=j+1

		if(quadratic) then
		Aomp(j,ixv(iCS,i-1))=Aomp(j,ixv(iCS,i-1))+(2d0*af*CloudR(i)+bf)*(Kd(i)*drho(i))
		Aomp(j,ixv(iCS,i+1))=Aomp(j,ixv(iCS,i+1))-(2d0*af*f2*CloudR(i)+bf*f1)*(Kd(i)*drho(i))
		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))+(2d0*af*(f2-1d0)*CloudR(i)+bf*(f1-1d0))*(Kd(i)*drho(i))

		Aomp(j,ixv(iCS,i-1))=Aomp(j,ixv(iCS,i-1))+2d0*af*Clouddens(i)*Kd(i)
		Aomp(j,ixv(iCS,i+1))=Aomp(j,ixv(iCS,i+1))-2d0*f2*af*Clouddens(i)*Kd(i)
		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))+2d0*(f2-1d0)*af*Clouddens(i)*Kd(i)
		else
		Aomp(j,ixv(iCS,i+1))=Aomp(j,ixv(iCS,i+1))+(Kd(i)*drho(i))/dz
		Aomp(j,ixv(iCS,i-1))=Aomp(j,ixv(iCS,i-1))-(Kd(i)*drho(i))/dz

		Aomp(j,ixv(iCS,i+1))=Aomp(j,ixv(iCS,i+1))+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i+1)-CloudR(i)))
		Aomp(j,ixv(iCS,i-1))=Aomp(j,ixv(iCS,i-1))+2d0*Clouddens(i)*Kd(i)/(dz*(CloudR(i)-CloudR(i-1)))
		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))-2d0*Clouddens(i)*Kd(i)*(1d0/(dz*(CloudR(i+1)-CloudR(i)))
     &					+1d0/(dz*(CloudR(i)-CloudR(i-1))))
		endif

		Aomp(j,ixv(iCS,i))=Aomp(j,ixv(iCS,i))-Sc(i)*xn(i)*Clouddens(i)/m_nuc
		Aomp(j,ixc(iCS,i))=Aomp(j,ixc(iCS,i))+Sc(i)*densv(iCS)/mpart(i)
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
	call DGESV( NN, NRHS, Aomp, NN, IWORKomp, xomp, NN, info )

	do i=1,NN
		if(.not.xomp(i).gt.1d-200) xomp(i)=1d-200
	enddo

	do i=1,nnr
		xc(iCS,i)=xomp(ixc(iCS,i))
		xv(iCS,i)=xomp(ixv(iCS,i))
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
		if(xn(i).gt.0d0) then
			tot=sum(xc(1:nCS,i))
			rr=(3d0*(tot*m_nuc/xn(i))/(4d0*pi*rho_av(i)))**(1d0/3d0)
			if(.not.rr.gt.r_nuc) rr=r_nuc
		else
			rr=r_nuc
		endif
		rpart(i)=sqrt(rr*rpart(i))
	enddo

	if(computeT.and.nTiter.lt.maxiter.and.nTiter.gt.0.and..not.doMCcompute.and..false.) then
		k=1
		do i=1,nr
			cloud_dens(i,ii)=0d0
			do j=1,nr_cloud
				do iCS=1,nCS
					cloud_dens(i,ii)=cloud_dens(i,ii)+xc(iCS,k)*Clouddens(k)/real(nr_cloud)
				enddo
				k=k+1
				if(k.gt.nnr) k=nnr
			enddo
		enddo
		call DoComputeT(Tconverged,0.1d0)
		
		k=1
		do i=1,nr-1
			do j=1,nr_cloud
				CloudT(k)=10d0**(log10(T(i))+log10(T(i+1)/T(i))*real(j-1)/real(nr_cloud))
				k=k+1
				if(k.gt.nnr) k=nnr
			enddo
		enddo
		CloudT(k)=T(nr)
	endif

	enddo
c end the loop

	k=1
	do i=1,nr
		cloud_dens(i,ii)=0d0
		Cloud(ii)%rv(i)=0d0
		Cloud(ii)%frac(i,1:18)=1d-200
		do j=1,nr_cloud
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

			do iCS=1,nCS
				cloud_dens(i,ii)=cloud_dens(i,ii)+xc(iCS,k)*Clouddens(k)/real(nr_cloud)
			enddo
			Cloud(ii)%rv(i)=Cloud(ii)%rv(i)+rpart(k)/real(nr_cloud)
			Cloud(ii)%frac(i,1:3)=Cloud(ii)%frac(i,1:3)+xc(1,k)/3d0		! TiO
			Cloud(ii)%frac(i,10)=Cloud(ii)%frac(i,10)+xc(3,k)			! Al2O3
			Cloud(ii)%frac(i,13:15)=Cloud(ii)%frac(i,13:15)+xc(5,k)/3d0	! Silicates
			Cloud(ii)%frac(i,8)=Cloud(ii)%frac(i,8)+xc(4,k)				! SiO2
			Cloud(ii)%frac(i,18)=Cloud(ii)%frac(i,18)+xc(6,k)			! H2O			
			Cloud(ii)%frac(i,9)=Cloud(ii)%frac(i,9)+xc(7,k)+xc(8,k)		! Fe + FeS
			Cloud(ii)%frac(i,17)=Cloud(ii)%frac(i,17)+xc(10,k)			! SiC			
			Cloud(ii)%frac(i,16)=Cloud(ii)%frac(i,16)+xc(9,k)			! C
			Cloud(ii)%frac(i,12)=Cloud(ii)%frac(i,12)+xMgO(k)			! MgO
			k=k+1
			if(k.gt.nnr) k=nnr
		enddo
		tot=sum(Cloud(ii)%frac(i,1:18))
		Cloud(ii)%frac(i,1:18)=Cloud(ii)%frac(i,1:18)/tot
	enddo

	open(unit=20,file=trim(outputdir) // '/cloudstructure.dat',RECL=1000)
	form='("#",a18,a19,a19,' // trim(int2string(nCS+1,'(i4)')) // 'a23,a19,a19,a19)'
	write(20,form) "P[bar]","dens[g/cm^3]","xn",(trim(CSname(i)),i=1,nCS),"MgO","r[micron]","T[K]","Jstar"
	form='(es19.7E3,es19.7E3,es19.7E3,' // trim(int2string(nCS+1,'(i4)')) // 'es23.7E3,es19.7E3,es19.7E3,es19.7E3)'
	do i=1,nnr
		densv=(mu*mp/(kb*CloudT(i)))*exp(BTP-ATP/CloudT(i))
		do iCS=1,nCS
			if(cloudT(i).gt.maxT(iCS)) densv(iCS)=densv(iCS)+(mu(iCS)*mp/(kb*CloudT(i)*10d0))
     &												*exp(BTP(iCS)-ATP(iCS)/(CloudT(i)*10d0))
		enddo
		write(20,form) CloudP(i),Clouddens(i),xn(i),xc(1:nCS,i),xMgO(i),rpart(i),CloudT(i),Sn(i)/m_nuc
	enddo
	close(unit=20)

	call ComputeTevap

	T=T0

	open(unit=20,file=trim(outputdir) // '/atoms.dat',RECL=6000)
	ini=.true.
	do i=1,nr
		call tellertje(i,nr)
		molfracs_atoms=molfracs_atoms0
		do iCS=1,nCS
			do j=1,N_atoms
				molfracs_atoms(j)=molfracs_atoms(j)+xv(iCS,(i-1)*nr_cloud+1)*mutot*atoms_cloud(iCS,j)/(mu(iCS)*CSnmol(iCS))
			enddo
		enddo
		molfracs_atoms(3)=molfracs_atoms(3)+COabun
		molfracs_atoms(5)=molfracs_atoms(5)+COabun
		do j=1,N_atoms
			if(.not.molfracs_atoms(j).gt.1d-10) molfracs_atoms(j)=1d-10
		enddo
		tot=sum(molfracs_atoms(1:N_atoms))
		molfracs_atoms=molfracs_atoms/tot
		call call_chemistry(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,.false.,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol)
		write(20,*) P(i),molfracs_atoms(1:N_atoms)
	enddo
	close(unit=20)

	deallocate(densv)
	deallocate(rpart)
	deallocate(mpart)
	deallocate(rho_av)
	deallocate(y)
	deallocate(Sn)
	deallocate(vth)
	deallocate(drho)
	deallocate(drhovsed)
	deallocate(xv)
	deallocate(xc)
	deallocate(xn)
	deallocate(tcinv)
	deallocate(vsed)
	deallocate(ixv)
	deallocate(ixc)
	deallocate(x)
	deallocate(IWORK)
	deallocate(An)

	return
	end
	
	
