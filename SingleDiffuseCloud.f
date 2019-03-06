	subroutine DiffuseCloud(ii)
	use GlobalSetup
	use Constants
	use AtomsModule
	IMPLICIT NONE
	real*8,allocatable :: x(:),vsed(:),xtot(:),vth(:),vthv(:)
	real*8,allocatable :: Sc(:),Sn(:),rpart(:),mpart(:)
	real*8,allocatable :: An(:,:),y(:,:),xv(:),xn(:),xc(:),A(:,:)
	real*8,allocatable :: drho(:),drhovsed(:),tcinv(:)
	real*8 K,dz,z12,z13,z12_2,z13_2,g,mu,rhodust,rr,Kc,Kp,xv_bot,densv,mutot,npart
	integer info,i,j,iter,NN,NRHS,niter,ii
	real*8 cs,err,maxerr,eps,SupSat,frac_nuc,m_nuc,tcoaginv,Dp,vmol,f
	real*8 af,bf,f1,f2,Pv,w_atoms(N_atoms),molfracs_atoms0(N_atoms),NKn
	integer,allocatable :: IWORK(:)
	real*8 sigmastar,Sigmadot,Pstar,gz,sigmamol,COabun,lmfp,fstick,kappa_cloud
	logical quadratic,ini
	integer atoms_cloud(N_atoms)
	character*500 cloudspecies(max(nclouds,1))
	integer nr_cloud,nnr
	parameter(nr_cloud=10)
	real*8 CloudP((nr-1)*nr_cloud+1),CloudT((nr-1)*nr_cloud+1),CloudR((nr-1)*nr_cloud+1),Clouddens((nr-1)*nr_cloud+1)

	real*8 Mc_top,Mn_top,IDP_dens,IDP_rad

	IDP_dens=0d0
	IDP_rad=0.01d0*micron

	Mc_top=-IDP_dens*sqrt(Ggrav*Mplanet/Dplanet)
	Mn_top=-Mc_top*(r_nuc/IDP_rad)**3

	quadratic=.true.

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

	if(Tform.gt.0d0) then
		call FormAbun(Tform,f_dry,f_wet,COratio,metallicity0,metallicity)
	else
		call easy_chem_set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
   	endif

	COabun=min(molfracs_atoms(3),molfracs_atoms(5))
	molfracs_atoms(3)=molfracs_atoms(3)-COabun
	molfracs_atoms(5)=molfracs_atoms(5)-COabun
c MgSiO3
	atoms_cloud=0
	atoms_cloud(7)=1
	atoms_cloud(9)=1
	atoms_cloud(5)=3
	mutot=0d0
	xv_bot=1d200
	do i=1,N_atoms
		mutot=mutot+w_atoms(i)*molfracs_atoms(i)
		if(atoms_cloud(i).gt.0) then
			f=molfracs_atoms(i)/real(atoms_cloud(i))
			if(f.lt.xv_bot) xv_bot=f
		endif
	enddo
	mu=0d0
	do i=1,N_atoms
		mu=mu+w_atoms(i)*real(atoms_cloud(i))
		molfracs_atoms(i)=molfracs_atoms(i)-xv_bot*real(atoms_cloud(i))
	enddo
	molfracs_atoms0=molfracs_atoms
	rhodust=2.8d0
	xv_bot=xv_bot*mu/mutot

	sigmastar=0.1
	Pstar=60d-6

	Pstar=Cloud(ii)%P
	sigmastar=log(Cloud(ii)%dP)

	fstick=1d0
	
	sigmamol=8d-15

	eps=1d-3

	K=Cloud(ii)%Kzz
	Sigmadot=Cloud(ii)%Sigmadot

	niter=40
	
	Kc=K
	Kp=K

c	AA=28030d0		!6.89d4
c	BB=12.471d0		!37.8
	
	m_nuc=4d0*pi*r_nuc**3*rhodust/3d0

	nnr=(nr-1)*nr_cloud+1

	allocate(rpart(nnr))
	allocate(mpart(nnr))
	allocate(y(nnr,5))
	allocate(Sc(nnr))
	allocate(Sn(nnr))
	allocate(vth(nnr))
	allocate(vthv(nnr))
	allocate(drho(nnr))
	allocate(drhovsed(nnr))
	allocate(xv(nnr))
	allocate(xc(nnr))
	allocate(xn(nnr))
	allocate(tcinv(nnr))
	allocate(vsed(nnr))

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
	allocate(A(NN,NN))

	allocate(An(nnr,nnr))

	xn=0d0
	xc=0d0
	xv=xv_bot
	tcinv=0d0

c start the loop
	do iter=1,niter
	if(iter.le.2) xv=xv_bot
	call tellertje(iter,niter)

	vsed=0d0
	do i=1,nnr
		cs=sqrt(kb*CloudT(i)/(2.3*mp))
		vth(i)=sqrt(8d0*kb*CloudT(i)/(pi*2.3*mp))
		vthv(i)=sqrt(8d0*kb*CloudT(i)/(pi*mu*mp))
		vsed(i)=-sqrt(pi)*rpart(i)*rhodust*Ggrav*Mplanet/(2d0*Clouddens(i)*vth(i)*CloudR(i)**2)
		mpart(i)=rhodust*4d0*pi*rpart(i)**3/3d0

		Sc(i)=min(vthv(i)*rpart(i),kb*CloudT(i)*sqrt(8d0*kb*CloudT(i)/(pi*2.3*mp))/(3d0*CloudP(i)*1d6*8e-15))
     &					*4d0*pi*rpart(i)*Clouddens(i)
		Sc(i)=fstick*Sc(i)
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
			af=((Clouddens(i-1)-Clouddens(i))-f2*(Clouddens(i+1)-Clouddens(i)))/((CloudR(i-1)**2-CloudR(i)**2)
     &				-f2*(CloudR(i+1)**2-CloudR(i)**2))
			bf=((Clouddens(i-1)-Clouddens(i))-f1*(Clouddens(i+1)-Clouddens(i)))/((CloudR(i-1)-CloudR(i))
     &				-f1*(CloudR(i+1)-CloudR(i)))
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

		Pv=1.04e17*exp(-58663d0/CloudT(i))
c		densv=10d0**(BB-AA/CloudT(i)-log10(CloudT(i)))
		densv=Pv/vthv(i)**2
		SupSat=Clouddens(i)*xv(i)/densv
		gz=Ggrav*Mplanet/CloudR(i)**2

		Sn(i)=(Clouddens(i)*gz*Sigmadot/(sigmastar*CloudP(i)*1d6*sqrt(2d0*pi)))*exp(-log(CloudP(i)/Pstar)**2/(2d0*sigmastar**2))

c		kappa_cloud=1d-6
c		Sn(i)=Clouddens(i)*gz*Sigmadot*(gz)*(kappa_cloud)*exp(-gz*kappa_cloud*CloudP(i)*1d6)
	enddo
	

c equations for Nuclii
	An=0d0
	x=0d0
	do i=2,nnr-1

		Pv=1.04e17*exp(-58663d0/CloudT(i))
c		densv=10d0**(BB-AA/CloudT(i)-log10(CloudT(i)))
		densv=Pv/vthv(i)**2

		j=i+1

		dz=CloudR(i+1)-CloudR(i-1)

		An(j,i)=An(j,i)-drhovsed(i)

		if(quadratic) then
		f1=(CloudR(i-1)**2-CloudR(i)**2)/(CloudR(i+1)**2-CloudR(i)**2)
		f2=(CloudR(i-1)-CloudR(i))/(CloudR(i+1)-CloudR(i))
		af=1d0/((CloudR(i-1)**2-CloudR(i)**2)-f2*(CloudR(i+1)**2-CloudR(i)**2))
		bf=1d0/((CloudR(i-1)-CloudR(i))-f1*(CloudR(i+1)-CloudR(i)))

		An(j,i-1)=An(j,i-1)+(2d0*af*CloudR(i)+bf)*(-Clouddens(i)*vsed(i)+Kp*drho(i))
		An(j,i+1)=An(j,i+1)-(2d0*af*f2*CloudR(i)+bf*f1)*(-Clouddens(i)*vsed(i)+Kp*drho(i))
		An(j,i)=An(j,i)+(2d0*af*(f2-1d0)*CloudR(i)+bf*(f1-1d0))*(-Clouddens(i)*vsed(i)+Kp*drho(i))

		An(j,i-1)=An(j,i-1)+2d0*af*Clouddens(i)*Kp
		An(j,i+1)=An(j,i+1)-2d0*f2*af*Clouddens(i)*Kp
		An(j,i)=An(j,i)+2d0*(f2-1d0)*af*Clouddens(i)*Kp
		else
		An(j,i+1)=An(j,i+1)+(-Clouddens(i)*vsed(i)+Kp*drho(i))/dz
		An(j,i-1)=An(j,i-1)-(-Clouddens(i)*vsed(i)+Kp*drho(i))/dz

		An(j,i+1)=An(j,i+1)+2d0*Clouddens(i)*Kp/(dz*(CloudR(i+1)-CloudR(i)))
		An(j,i-1)=An(j,i-1)+2d0*Clouddens(i)*Kp/(dz*(CloudR(i)-CloudR(i-1)))
		An(j,i)=An(j,i)-2d0*Clouddens(i)*Kp*(1d0/(dz*(CloudR(i+1)-CloudR(i)))+1d0/(dz*(CloudR(i)-CloudR(i-1))))
		endif

		x(j)=-Sn(i)

c coagulation
		if(coagulation) then
			npart=xn(i)*Clouddens(i)/m_nuc
			lmfp=2.3*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
			vmol=0.5d0*lmfp*vth(i)
			Dp=kb*CloudT(i)/(6d0*pi*rpart(i)*vmol*Clouddens(i))

c			tcoaginv=npart*pi*rpart(i)**2*abs(vsed(i))
c rewritten for better convergence
			tcoaginv=sqrt(pi)*3d0*xc(i)*Ggrav*Mplanet/(8d0*vth(i)*CloudR(i)**2)

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
	An(1,i)=Kp/dz-vsed(i)
	An(1,i-1)=-Kp/dz
	x(1)=-Mn_top/Clouddens(i)
	i=1
	An(2,i)=1d0
	x(2)=0d0

	NRHS=1
	call DGESV( nnr, NRHS, An, nnr, IWORK, x, nnr, info )

	do i=1,nnr
		if(.not.x(i).gt.0d0) x(i)=0d0
	enddo
	xn(1:nnr)=x(1:nnr)

c equations for material
	A=0d0
	x=0d0
	do i=2,nnr-1
		dz=CloudR(i+1)-CloudR(i-1)
		Pv=1.04e17*exp(-58663d0/CloudT(i))
c		densv=10d0**(BB-AA/CloudT(i)-log10(CloudT(i)))
		densv=Pv/vthv(i)**2

		f1=(CloudR(i-1)**2-CloudR(i)**2)/(CloudR(i+1)**2-CloudR(i)**2)
		f2=(CloudR(i-1)-CloudR(i))/(CloudR(i+1)-CloudR(i))
		af=1d0/((CloudR(i-1)**2-CloudR(i)**2)-f2*(CloudR(i+1)**2-CloudR(i)**2))
		bf=1d0/((CloudR(i-1)-CloudR(i))-f1*(CloudR(i+1)-CloudR(i)))

		j=i+1

		A(j,i)=A(j,i)-drhovsed(i)

		if(quadratic) then
		A(j,i-1)=A(j,i-1)+(2d0*af*CloudR(i)+bf)*(-Clouddens(i)*vsed(i)+Kp*drho(i))
		A(j,i+1)=A(j,i+1)-(2d0*af*f2*CloudR(i)+bf*f1)*(-Clouddens(i)*vsed(i)+Kp*drho(i))
		A(j,i)=A(j,i)+(2d0*af*(f2-1d0)*CloudR(i)+bf*(f1-1d0))*(-Clouddens(i)*vsed(i)+Kp*drho(i))

		A(j,i-1)=A(j,i-1)+2d0*af*Clouddens(i)*Kp
		A(j,i+1)=A(j,i+1)-2d0*f2*af*Clouddens(i)*Kp
		A(j,i)=A(j,i)+2d0*(f2-1d0)*af*Clouddens(i)*Kp
		else
		A(j,i+1)=A(j,i+1)+(-Clouddens(i)*vsed(i)+Kc*drho(i))/dz
		A(j,i-1)=A(j,i-1)-(-Clouddens(i)*vsed(i)+Kc*drho(i))/dz

		A(j,i+1)=A(j,i+1)+2d0*Clouddens(i)*Kc/(dz*(CloudR(i+1)-CloudR(i)))
		A(j,i-1)=A(j,i-1)+2d0*Clouddens(i)*Kp/(dz*(CloudR(i)-CloudR(i-1)))
		A(j,i)=A(j,i)-2d0*Clouddens(i)*Kp*(1d0/(dz*(CloudR(i+1)-CloudR(i)))+1d0/(dz*(CloudR(i)-CloudR(i-1))))
		endif

		A(j,nnr+i)=A(j,nnr+i)+Sc(i)*xn(i)*Clouddens(i)/m_nuc
		A(j,i)=A(j,i)-Sc(i)*(densv)/mpart(i)
		x(j)=0d0

		j=i+nnr+1

		if(quadratic) then
		A(j,i+nnr-1)=A(j,i+nnr-1)+(2d0*af*CloudR(i)+bf)*(Kc*drho(i))
		A(j,i+nnr+1)=A(j,i+nnr+1)-(2d0*af*f2*CloudR(i)+bf*f1)*(Kc*drho(i))
		A(j,i+nnr)=A(j,i+nnr)+(2d0*af*(f2-1d0)*CloudR(i)+bf*(f1-1d0))*(Kc*drho(i))

		A(j,i+nnr-1)=A(j,i+nnr-1)+2d0*af*Clouddens(i)*Kc
		A(j,i+nnr+1)=A(j,i+nnr+1)-2d0*f2*af*Clouddens(i)*Kc
		A(j,i+nnr)=A(j,i+nnr)+2d0*(f2-1d0)*af*Clouddens(i)*Kc
		else
		A(j,i+nnr+1)=A(j,i+nnr+1)+(Kc*drho(i))/dz
		A(j,i+nnr-1)=A(j,i+nnr-1)-(Kc*drho(i))/dz

		A(j,i+nnr+1)=A(j,i+nnr+1)+2d0*Clouddens(i)*Kc/(dz*(CloudR(i+1)-CloudR(i)))
		A(j,i+nnr-1)=A(j,i+nnr-1)+2d0*Clouddens(i)*Kc/(dz*(CloudR(i)-CloudR(i-1)))
		A(j,i+nnr)=A(j,i+nnr)-2d0*Clouddens(i)*Kc*(1d0/(dz*(CloudR(i+1)-CloudR(i)))+1d0/(dz*(CloudR(i)-CloudR(i-1))))		
		endif

		A(j,nnr+i)=A(j,nnr+i)-Sc(i)*xn(i)*Clouddens(i)/m_nuc
		A(j,i)=A(j,i)+Sc(i)*(densv)/mpart(i)
		x(j)=0d0
	enddo
	i=nnr
	dz=CloudR(i)-CloudR(i-1)
	A(nnr+1,i)=Kp/dz-vsed(i)
	A(nnr+1,i-1)=-Kp/dz
	x(nnr+1)=-Mc_top/Clouddens(i)
	A(nnr+2,nnr+i)=Kc/dz
	A(nnr+2,nnr+i-1)=-Kc/dz
	x(nnr+2)=0d0!Mc_top/Clouddens(i)
	i=1
	A(1,i)=1d0
	x(1)=0d0
	A(2,nnr+i)=1d0
	x(2)=xv_bot

	NRHS=1
	call DGESV( NN, NRHS, A(1:NN,1:NN), NN, IWORK, x, NN, info )

	do i=1,NN
		if(.not.x(i).gt.1d-100) x(i)=1d-100
	enddo

	xc(1:nnr)=x(1:nnr)
	xv(1:nnr)=x(nnr+1:nnr*2)
	
	do i=1,nnr
		if(xc(i).lt.0d0) xc(i)=0d0
		if(xv(i).lt.0d0) xv(i)=0d0
		if(xn(i).lt.0d0) xn(i)=0d0
	enddo

	open(unit=20,file='output.dat',RECL=1000)
	do i=1,nnr
		Pv=1.04e17*exp(-58663d0/CloudT(i))
c		densv=10d0**(BB-AA/CloudT(i)-log10(CloudT(i)))
		densv=Pv/vthv(i)**2
		write(20,*) CloudP(i),Clouddens(i),xn(i),xc(i),xv(i),rpart(i),Clouddens(i)*xv(i)/densv,CloudT(i),Sn(i)
	enddo
	close(unit=20)

	do i=1,nnr
		if(xn(i).gt.0d0) then
			rr=(3d0*(xc(i)*m_nuc/xn(i))/(4d0*pi*rhodust))**(1d0/3d0)
			if(rr.lt.r_nuc) rr=r_nuc
		else
			rr=r_nuc
		endif
		rpart(i)=sqrt(rr*rpart(i))
	enddo

	enddo
c end the loop

	open(unit=20,file='output.dat',RECL=1000)
	do i=1,nnr
		Pv=1.04e17*exp(-58663d0/CloudT(i))
c		densv=10d0**(BB-AA/CloudT(i)-log10(CloudT(i)))
		densv=Pv/vthv(i)**2
		write(20,*) CloudP(i),Clouddens(i),xn(i),xc(i),xv(i),rpart(i),Clouddens(i)*xv(i)/densv,CloudT(i),Sn(i)
	enddo
	close(unit=20)

	k=1
	do i=1,nr-1
		cloud_dens(i,ii)=0d0
		Cloud(ii)%rv(i)=0d0
		do j=1,nr_cloud
			cloud_dens(i,ii)=cloud_dens(i,ii)+xc(k)*Clouddens(k)/real(nr_cloud)
			Cloud(ii)%rv(i)=Cloud(ii)%rv(i)+rpart(k)/real(nr_cloud)
			k=k+1
		enddo
	enddo
	cloud_dens(nr,ii)=xc(k)*Clouddens(k)
	Cloud(ii)%rv(nr)=rpart(k)

	do i=1,nr
		call tellertje(i,nr)
		ini=.true.
		do j=1,N_atoms
			molfracs_atoms(j)=molfracs_atoms0(j)+xv((i-1)*nr_cloud+1)*mutot*real(atoms_cloud(j))/mu
		enddo
		molfracs_atoms(3)=molfracs_atoms(3)+COabun
		molfracs_atoms(5)=molfracs_atoms(5)+COabun
		call call_easy_chem(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,.false.,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),.false.,f_dry,MMW(i),didcondens(i),fast_chem,includemol)
	enddo

	deallocate(rpart)
	deallocate(mpart)
	deallocate(y)
	deallocate(Sc)
	deallocate(Sn)
	deallocate(vth)
	deallocate(vthv)
	deallocate(drho)
	deallocate(drhovsed)
	deallocate(xv)
	deallocate(xc)
	deallocate(xn)
	deallocate(tcinv)
	deallocate(vsed)
	deallocate(x)
	deallocate(IWORK)
	deallocate(A)
	deallocate(An)

	return
	end
	
	
