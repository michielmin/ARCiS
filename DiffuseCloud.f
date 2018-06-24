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
	real*8 cs,err,maxerr,eps,SupSat,frac_nuc,r_nuc,m_nuc,tcoaginv,Dp,vmol,f
	real*8 af,bf,f1,f2,Pv,w_atoms(N_atoms),molfracs_atoms0(N_atoms),NKn
	integer,allocatable :: IWORK(:)
	real*8 sigmastar,Sigmadot,Pstar,gz,sigmamol,COabun,lmfp,fstick
	logical quadratic,ini
	integer atoms_cloud(N_atoms)
	character*500 cloudspecies(max(nclouds,1))

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
		call FormAbun(Tform,f_enrich,COratio,metallicity0,metallicity)
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

	niter=100
	
	Kc=K
	Kp=K

c	AA=28030d0		!6.89d4
c	BB=12.471d0		!37.8
	
	r_nuc=0.001*1d-4
	m_nuc=4d0*pi*r_nuc**3*rhodust/3d0

	allocate(rpart(nr))
	allocate(mpart(nr))
	allocate(y(nr,5))
	allocate(Sc(nr))
	allocate(Sn(nr))
	allocate(vth(nr))
	allocate(vthv(nr))
	allocate(drho(nr))
	allocate(drhovsed(nr))
	allocate(xv(nr))
	allocate(xc(nr))
	allocate(xn(nr))
	allocate(tcinv(nr))
	allocate(vsed(nr))

	f=0.1d0
	rpart=r_nuc

	NN=2*nr
	allocate(x(NN))
	allocate(IWORK(10*NN*NN))
	allocate(A(NN,NN))

	allocate(An(nr,nr))

	xn=0d0
	xc=0d0
	xv=xv_bot
	tcinv=0d0

c start the loop
	do iter=1,niter
	if(iter.le.2) xv=xv_bot
	call tellertje(iter,niter)

	vsed=0d0
	do i=1,nr
		cs=sqrt(kb*T(i)/(2.3*mp))
		vth(i)=sqrt(8d0*kb*T(i)/(pi*2.3*mp))
		vthv(i)=sqrt(8d0*kb*T(i)/(pi*mu*mp))
		vsed(i)=-sqrt(pi)*rpart(i)*rhodust*Ggrav*Mplanet/(2d0*dens(i)*vth(i)*R(i)**2)
		mpart(i)=rhodust*4d0*pi*rpart(i)**3/3d0

		Sc(i)=min(vthv(i)*rpart(i),kb*T(i)*sqrt(8d0*kb*T(i)/(pi*2.3*mp))/(3d0*P(i)*1d6*8e-15))*4d0*pi*rpart(i)*dens(i)
		Sc(i)=fstick*Sc(i)
	enddo

	do i=1,nr
		if(i.eq.1) then
			drho(i)=(dens(i+1)-dens(i))/(R(i+1)-R(i))
			drhovsed(i)=(vsed(i+1)*dens(i+1)-vsed(i)*dens(i))/(R(i+1)-R(i))
		else if(i.eq.nr) then
			drho(i)=(dens(i)-dens(i-1))/(R(i)-R(i-1))
			drhovsed(i)=(vsed(i)*dens(i)-vsed(i-1)*dens(i-1))/(R(i)-R(i-1))
		else
			if(quadratic) then
			f1=(R(i-1)**2-R(i)**2)/(R(i+1)**2-R(i)**2)
			f2=(R(i-1)-R(i))/(R(i+1)-R(i))
			af=((dens(i-1)-dens(i))-f2*(dens(i+1)-dens(i)))/((R(i-1)**2-R(i)**2)-f2*(R(i+1)**2-R(i)**2))
			bf=((dens(i-1)-dens(i))-f1*(dens(i+1)-dens(i)))/((R(i-1)-R(i))-f1*(R(i+1)-R(i)))
			drho(i)=2d0*af*R(i)+bf
			af=((vsed(i-1)*dens(i-1)-vsed(i)*dens(i))-f2*(vsed(i+1)*dens(i+1)-vsed(i)*dens(i)))/
     &				((R(i-1)**2-R(i)**2)-f2*(R(i+1)**2-R(i)**2))
			bf=((vsed(i-1)*dens(i-1)-vsed(i)*dens(i))-f1*(vsed(i+1)*dens(i+1)-vsed(i)*dens(i)))/
     &				((R(i-1)-R(i))-f1*(R(i+1)-R(i)))
			drhovsed(i)=2d0*af*R(i)+bf
			else
			drho(i)=(dens(i+1)-dens(i-1))/(R(i+1)-R(i-1))
			drhovsed(i)=(vsed(i+1)*dens(i+1)-vsed(i-1)*dens(i-1))/(R(i+1)-R(i-1))
			endif
		endif

		Pv=1.04e17*exp(-58663d0/T(i))
c		densv=10d0**(BB-AA/T(i)-log10(T(i)))
		densv=Pv/vthv(i)**2
		SupSat=dens(i)*xv(i)/densv
		gz=Ggrav*Mplanet/R(i)**2

		Sn(i)=(dens(i)*gz*Sigmadot/(sigmastar*P(i)*1d6*sqrt(2d0*pi)))*exp(-log(P(i)/Pstar)**2/(2d0*sigmastar**2))
	enddo
	

c equations for Nuclii
	An=0d0
	x=0d0
	do i=2,nr-1

		Pv=1.04e17*exp(-58663d0/T(i))
c		densv=10d0**(BB-AA/T(i)-log10(T(i)))
		densv=Pv/vthv(i)**2

		j=i+1

		dz=R(i+1)-R(i-1)

		An(j,i)=An(j,i)-drhovsed(i)

		if(quadratic) then
		f1=(R(i-1)**2-R(i)**2)/(R(i+1)**2-R(i)**2)
		f2=(R(i-1)-R(i))/(R(i+1)-R(i))
		af=1d0/((R(i-1)**2-R(i)**2)-f2*(R(i+1)**2-R(i)**2))
		bf=1d0/((R(i-1)-R(i))-f1*(R(i+1)-R(i)))

		An(j,i-1)=An(j,i-1)+(2d0*af*R(i)+bf)*(-dens(i)*vsed(i)+Kp*drho(i))
		An(j,i+1)=An(j,i+1)-(2d0*af*f2*R(i)+bf*f1)*(-dens(i)*vsed(i)+Kp*drho(i))
		An(j,i)=An(j,i)+(2d0*af*(f2-1d0)*R(i)+bf*(f1-1d0))*(-dens(i)*vsed(i)+Kp*drho(i))

		An(j,i-1)=An(j,i-1)+2d0*af*dens(i)*Kp
		An(j,i+1)=An(j,i+1)-2d0*f2*af*dens(i)*Kp
		An(j,i)=An(j,i)+2d0*(f2-1d0)*af*dens(i)*Kp
		else
		An(j,i+1)=An(j,i+1)+(-dens(i)*vsed(i)+Kp*drho(i))/dz
		An(j,i-1)=An(j,i-1)-(-dens(i)*vsed(i)+Kp*drho(i))/dz

		An(j,i+1)=An(j,i+1)+2d0*dens(i)*Kp/(dz*(R(i+1)-R(i)))
		An(j,i-1)=An(j,i-1)+2d0*dens(i)*Kp/(dz*(R(i)-R(i-1)))
		An(j,i)=An(j,i)-2d0*dens(i)*Kp*(1d0/(dz*(R(i+1)-R(i)))+1d0/(dz*(R(i)-R(i-1))))
		endif

		x(j)=-Sn(i)

c coagulation
		if(coagulation) then
			npart=xn(i)*dens(i)/m_nuc
			lmfp=2.3*mp/(sqrt(2d0)*dens(i)*sigmamol)
			vmol=0.5d0*lmfp*vth(i)
			Dp=kb*T(i)/(6d0*pi*rpart(i)*vmol*dens(i))

c			tcoaginv=npart*pi*rpart(i)**2*abs(vsed(i))
c rewritten for better convergence
			tcoaginv=sqrt(pi)*3d0*xc(i)*Ggrav*Mplanet/(8d0*vth(i)*R(i)**2)

			if(sqrt(16d0*kb*T(i)/(pi*mpart(i)))*rpart(i).lt.Dp) then
				tcoaginv=tcoaginv+2d0*pi*rpart(i)**2*npart*sqrt(16d0*kb*T(i)/(pi*mpart(i)))
			else
				tcoaginv=tcoaginv+2d0*pi*rpart(i)*npart*Dp
			endif

			if(.not.tcoaginv.gt.0d0) tcoaginv=0d0

			tcinv(i)=(tcinv(i)+tcoaginv)/2d0

			An(j,i)=An(j,i)-dens(i)*tcinv(i)
		endif
	enddo
	i=nr
	dz=R(i)-R(i-1)
	An(1,i)=Kp/dz+vsed(i)
	An(1,i-1)=-Kp/dz
	x(1)=0d0
	i=1
	An(2,i)=1d0
	x(2)=0d0

	NRHS=1
	call DGESV( nr, NRHS, An, nr, IWORK, x, nr, info )

	do i=1,nr
		if(.not.x(i).gt.0d0) x(i)=0d0
	enddo
	xn(1:nr)=x(1:nr)

c equations for material
	A=0d0
	x=0d0
	do i=2,nr-1
		dz=R(i+1)-R(i-1)
		Pv=1.04e17*exp(-58663d0/T(i))
c		densv=10d0**(BB-AA/T(i)-log10(T(i)))
		densv=Pv/vthv(i)**2

		f1=(R(i-1)**2-R(i)**2)/(R(i+1)**2-R(i)**2)
		f2=(R(i-1)-R(i))/(R(i+1)-R(i))
		af=1d0/((R(i-1)**2-R(i)**2)-f2*(R(i+1)**2-R(i)**2))
		bf=1d0/((R(i-1)-R(i))-f1*(R(i+1)-R(i)))

		j=i+1

		A(j,i)=A(j,i)-drhovsed(i)

		if(quadratic) then
		A(j,i-1)=A(j,i-1)+(2d0*af*R(i)+bf)*(-dens(i)*vsed(i)+Kp*drho(i))
		A(j,i+1)=A(j,i+1)-(2d0*af*f2*R(i)+bf*f1)*(-dens(i)*vsed(i)+Kp*drho(i))
		A(j,i)=A(j,i)+(2d0*af*(f2-1d0)*R(i)+bf*(f1-1d0))*(-dens(i)*vsed(i)+Kp*drho(i))

		A(j,i-1)=A(j,i-1)+2d0*af*dens(i)*Kp
		A(j,i+1)=A(j,i+1)-2d0*f2*af*dens(i)*Kp
		A(j,i)=A(j,i)+2d0*(f2-1d0)*af*dens(i)*Kp
		else
		A(j,i+1)=A(j,i+1)+(-dens(i)*vsed(i)+Kc*drho(i))/dz
		A(j,i-1)=A(j,i-1)-(-dens(i)*vsed(i)+Kc*drho(i))/dz

		A(j,i+1)=A(j,i+1)+2d0*dens(i)*Kc/(dz*(R(i+1)-R(i)))
		A(j,i-1)=A(j,i-1)+2d0*dens(i)*Kp/(dz*(R(i)-R(i-1)))
		A(j,i)=A(j,i)-2d0*dens(i)*Kp*(1d0/(dz*(R(i+1)-R(i)))+1d0/(dz*(R(i)-R(i-1))))
		endif

		A(j,nr+i)=A(j,nr+i)+Sc(i)*xn(i)*dens(i)/m_nuc
		A(j,i)=A(j,i)-Sc(i)*(densv)/mpart(i)
		x(j)=0d0

		j=i+nr+1

		if(quadratic) then
		A(j,i+nr-1)=A(j,i+nr-1)+(2d0*af*R(i)+bf)*(Kc*drho(i))
		A(j,i+nr+1)=A(j,i+nr+1)-(2d0*af*f2*R(i)+bf*f1)*(Kc*drho(i))
		A(j,i+nr)=A(j,i+nr)+(2d0*af*(f2-1d0)*R(i)+bf*(f1-1d0))*(Kc*drho(i))

		A(j,i+nr-1)=A(j,i+nr-1)+2d0*af*dens(i)*Kc
		A(j,i+nr+1)=A(j,i+nr+1)-2d0*f2*af*dens(i)*Kc
		A(j,i+nr)=A(j,i+nr)+2d0*(f2-1d0)*af*dens(i)*Kc
		else
		A(j,i+nr+1)=A(j,i+nr+1)+(Kc*drho(i))/dz
		A(j,i+nr-1)=A(j,i+nr-1)-(Kc*drho(i))/dz

		A(j,i+nr+1)=A(j,i+nr+1)+2d0*dens(i)*Kc/(dz*(R(i+1)-R(i)))
		A(j,i+nr-1)=A(j,i+nr-1)+2d0*dens(i)*Kc/(dz*(R(i)-R(i-1)))
		A(j,i+nr)=A(j,i+nr)-2d0*dens(i)*Kc*(1d0/(dz*(R(i+1)-R(i)))+1d0/(dz*(R(i)-R(i-1))))		
		endif

		A(j,nr+i)=A(j,nr+i)-Sc(i)*xn(i)*dens(i)/m_nuc
		A(j,i)=A(j,i)+Sc(i)*(densv)/mpart(i)
		x(j)=0d0
	enddo
	i=nr
	dz=R(i)-R(i-1)
	A(nr+1,i)=Kp/dz+vsed(i)
	A(nr+1,i-1)=-Kp/dz
	x(nr+1)=0d0
	A(nr+2,nr+i)=Kc/dz
	A(nr+2,nr+i-1)=-Kc/dz
	x(nr+2)=0d0
	i=1
	A(1,i)=1d0
	x(1)=0d0
	A(2,nr+i)=1d0
	x(2)=xv_bot

	NRHS=1
	call DGESV( NN, NRHS, A(1:NN,1:NN), NN, IWORK, x, NN, info )

	do i=1,NN
		if(.not.x(i).gt.1d-100) x(i)=1d-100
	enddo

	xc(1:nr)=x(1:nr)
	xv(1:nr)=x(nr+1:nr*2)
	
	do i=1,nr
		if(xc(i).lt.0d0) xc(i)=0d0
		if(xv(i).lt.0d0) xv(i)=0d0
		if(xn(i).lt.0d0) xn(i)=0d0
	enddo

	open(unit=20,file='output.dat',RECL=1000)
	do i=1,nr
		Pv=1.04e17*exp(-58663d0/T(i))
c		densv=10d0**(BB-AA/T(i)-log10(T(i)))
		densv=Pv/vthv(i)**2
		write(20,*) P(i),dens(i),xn(i),xc(i),xv(i),rpart(i),dens(i)*xv(i)/densv,T(i)
	enddo
	close(unit=20)

	do i=1,nr
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
	do i=1,nr
		Pv=1.04e17*exp(-58663d0/T(i))
c		densv=10d0**(BB-AA/T(i)-log10(T(i)))
		densv=Pv/vthv(i)**2
		write(20,*) P(i),dens(i),xn(i),xc(i),xv(i),rpart(i),dens(i)*xv(i)/densv,T(i)
	enddo
	close(unit=20)

	cloud_dens(1:nr,ii)=xc(1:nr)*dens(1:nr)
	Cloud(ii)%rv(1:nr)=rpart(1:nr)

	do i=1,nr
		call tellertje(i,nr)
		ini=.true.
		do j=1,N_atoms
			molfracs_atoms(j)=molfracs_atoms0(j)+xv(i)*mutot*real(atoms_cloud(j))/mu
		enddo
		molfracs_atoms(3)=molfracs_atoms(3)+COabun
		molfracs_atoms(5)=molfracs_atoms(5)+COabun
		call call_easy_chem(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,.false.,cloudspecies,
     &						XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),.false.,f_enrich,MMW(i),didcondens(i))
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
	
	
