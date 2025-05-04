	subroutine SetupCloud(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 pp,tot,column,tt,z,densdust,CloudMass,x
	integer ii,i,j,nsubr,isize,ilam
	real*8 Xc,Xc1,lambdaC,Ca,Cs,tau,P_SI1,P_SI2,veff,frac(nr,10)
	real*8 CPmin,CPmax,CPtau
	logical cl

	if(.not.allocated(Cloud(ii)%rv)) then
		allocate(Cloud(ii)%rv(nr))
		allocate(Cloud(ii)%sigma(nr))
	endif
	
	if(Cloud(ii)%x_slider.ge.1d0.and.Cloud(ii)%x_slider.lt.60d0) then
		i=Cloud(ii)%x_slider
		Cloud(ii)%abun(1:60)=0d0
		Cloud(ii)%abun(i)=1d0-(Cloud(ii)%x_slider-real(i))
		Cloud(ii)%abun(i+1)=1d0-Cloud(ii)%abun(i)
	endif
	select case(Cloud(ii)%type)
		case("DIFFUSE")
			Cloud(ii)%nlam=nlam
			call DiffuseCloud(ii)
			Cloud(ii)%rv=Cloud(ii)%rv*1d4
			Cloud(ii)%sigma=1d-10
			Cloud(ii)%onepart=.false.
			call SetupPartCloud(ii)
		case("WATER")
			Cloud(ii)%nlam=nlam
			call WaterCloud(ii)
			Cloud(ii)%rv=Cloud(ii)%rv*1d4
			Cloud(ii)%sigma=1d-10
			Cloud(ii)%onepart=.false.
			call SetupPartCloud(ii)
		case("CONDENSATION")
			Cloud(ii)%nlam=nlam
			call CondensationCloud(ii)
			Cloud(ii)%rv=Cloud(ii)%rv*1d4
			Cloud(ii)%sigma=1d-10
			Cloud(ii)%onepart=.false.
			call SetupPartCloud(ii)
		case("FILE")
			Cloud(ii)%nlam=nlam
			call regridN(Cloud(ii)%file,P,cloud_dens(1:nr,ii),nr,2,6,1,1,.false.,.false.)
			call regridN(Cloud(ii)%file,P,Cloud(ii)%rv(1:nr),nr,2,5,1,1,.false.,.true.)
c 10% iron
			Cloud(ii)%frac(1:nr,9)=0.1d0
c 90% MgSiO3
			Cloud(ii)%frac(1:nr,13:15)=0.9d0/3d0
			densdust=3.7
			Cloud(ii)%hazetype='THOLIN'
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*(4d0*pi*densdust*Cloud(ii)%rv(1:nr)**3)/3d0
			Cloud(ii)%rv=Cloud(ii)%rv*1d4
			Cloud(ii)%sigma=1d-10
			Cloud(ii)%onepart=.false.
			call SetupPartCloud(ii)
		case("FILEDRIFT")
			Cloud(ii)%nlam=nlam
			call regridN(Cloud(ii)%file,P*1d6,frac,nr,2,9,10,4,.true.,.true.)
			Cloud(ii)%frac(1:nr,1)=frac(1:nr,1)/3d0
			Cloud(ii)%frac(1:nr,2)=frac(1:nr,1)/3d0
			Cloud(ii)%frac(1:nr,3)=frac(1:nr,1)/3d0

			Cloud(ii)%frac(1:nr,4)=frac(1:nr,2)/3d0
			Cloud(ii)%frac(1:nr,5)=frac(1:nr,2)/3d0
			Cloud(ii)%frac(1:nr,6)=frac(1:nr,2)/3d0

			Cloud(ii)%frac(1:nr,7)=frac(1:nr,3)
			Cloud(ii)%frac(1:nr,8)=frac(1:nr,4)
			Cloud(ii)%frac(1:nr,9)=frac(1:nr,5)
			Cloud(ii)%frac(1:nr,10)=frac(1:nr,6)
			Cloud(ii)%frac(1:nr,11)=frac(1:nr,7)
			Cloud(ii)%frac(1:nr,12)=frac(1:nr,8)

			Cloud(ii)%frac(1:nr,13)=frac(1:nr,9)/3d0
			Cloud(ii)%frac(1:nr,14)=frac(1:nr,9)/3d0
			Cloud(ii)%frac(1:nr,15)=frac(1:nr,9)/3d0

			Cloud(ii)%frac(1:nr,16)=frac(1:nr,10)

			Cloud(ii)%frac(1:nr,17)=0d0
			Cloud(ii)%frac(1:nr,18)=0d0

			call regridN(Cloud(ii)%file,P*1d6,cloud_dens(1:nr,ii),nr,2,43,1,4,.true.,.false.)
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)

			call regridN(Cloud(ii)%file,P*1d6,Cloud(ii)%rv(1:nr),nr,2,7,1,4,.true.,.false.)
			Cloud(ii)%rv=Cloud(ii)%rv*1d4
			Cloud(ii)%sigma=1d-10
			Cloud(ii)%onepart=.false.
			call SetupPartCloud(ii)
		case("LAYER")
			Cloud(ii)%nlam=nlam+1
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			Cloud(ii)%onepart=(Cloud(ii)%rpow.eq.0d0)
			do i=1,nr
				Cloud(ii)%frac(i,1:60)=Cloud(ii)%abun(1:60)
			enddo
			call SetupPartCloud(ii)
			CPmin=Cloud(ii)%Pmin
			CPmax=Cloud(ii)%Pmax
			CPtau=Cloud(ii)%Ptau
			if(CPmin.gt.CPmax) then
				x=CPmin
				CPmin=CPmax
				CPmax=x
			endif
			if(CPtau.lt.0d0) then
				CPtau=CPmax
			endif
			if(CPtau.le.CPmin) then
				CPtau=CPmin*1.0001
			endif
			do i=1,nr
				Cloud(ii)%Kref=Cloud(ii)%Kext(i,nlam+1)
				if(P(i).gt.CPmin.and.P(i).lt.CPmax) then
					if(abs(Cloud(ii)%xi).gt.1d-8) then
						cloud_dens(i,ii)=(grav(i)*Cloud(ii)%xi*Cloud(ii)%tau*P(i)**(Cloud(ii)%xi-1d0))/
     &					(Cloud(ii)%Kref*1d6*(CPtau**Cloud(ii)%xi-CPmin**Cloud(ii)%xi))
					else
						cloud_dens(i,ii)=(grav(i)*Cloud(ii)%tau*P(i)**(Cloud(ii)%xi-1d0))/
     &					(Cloud(ii)%Kref*1d6*log(CPtau/CPmin))
					endif
				else
					cloud_dens(i,ii)=0d0
				endif
			enddo
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)
			if(Cloud(ii)%Ptau.ge.Cloud(ii)%Pmax) then
				tot=0d0
				do i=1,nr
					tot=tot+cloud_dens(i,ii)*Cloud(ii)%Kext(i,nlam+1)*(R(i+1)-R(i))
				enddo
				if(tot.gt.0d0) then
					cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*Cloud(ii)%tau/tot
				endif
			endif
		case("SLAB")
			Cloud(ii)%nlam=nlam+1
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			Cloud(ii)%onepart=(Cloud(ii)%rpow.eq.0d0)
			do i=1,nr
				Cloud(ii)%frac(i,1:60)=Cloud(ii)%abun(1:60)
			enddo
			call SetupPartCloud(ii)
			CPmin=Cloud(ii)%Pmin
			CPmax=Cloud(ii)%Pmax
			CPtau=Cloud(ii)%Ptau
			if(CPmin.gt.CPmax) then
				x=CPmin
				CPmin=CPmax
				CPmax=x
			endif
			do i=1,nr
				Cloud(ii)%Kref=Cloud(ii)%Kext(i,nlam+1)
				if(P(i).gt.CPmin.and.P(i).lt.CPmax) then
					cloud_dens(i,ii)=(grav(i)*2d0*Cloud(ii)%tau*P(i))/
     &					(Cloud(ii)%Kref*1d6*(CPmax*CPmax-CPmin*CPmin))
				else
					cloud_dens(i,ii)=0d0
				endif
			enddo
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)
		case("DECK")
			Cloud(ii)%nlam=nlam+1
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			Cloud(ii)%onepart=(Cloud(ii)%rpow.eq.0d0)
			do i=1,nr
				Cloud(ii)%frac(i,1:60)=Cloud(ii)%abun(1:60)
			enddo
			call SetupPartCloud(ii)
			do i=1,nr
				Cloud(ii)%Kref=Cloud(ii)%Kext(i,nlam+1)
				cloud_dens(i,ii)=(grav(i)*exp((P(i)-Cloud(ii)%Ptau)/Cloud(ii)%Phi))/
     &				(Cloud(ii)%Kref*1d6*Cloud(ii)%Phi*(1d0-exp(-Cloud(ii)%Ptau/Cloud(ii)%Phi)))
			enddo
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)
		case("HOMOGENEOUS")
			Cloud(ii)%nlam=nlam
			do i=1,nr
				if(P(i).lt.Cloud(ii)%Pmin) then
					cloud_dens(i,ii)=0d0
				else
					cloud_dens(i,ii)=dens(i)*Cloud(ii)%mixrat
				endif
			enddo
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			Cloud(ii)%onepart=(Cloud(ii)%rpow.eq.0d0)
			do i=1,nr
				Cloud(ii)%frac(i,1:60)=Cloud(ii)%abun(1:60)
			enddo
			call SetupPartCloud(ii)
		case("GAUSS","HALFGAUSS")
			Cloud(ii)%nlam=nlam+1
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			Cloud(ii)%onepart=(Cloud(ii)%rpow.eq.0d0)
			do i=1,nr
				Cloud(ii)%frac(i,1:60)=Cloud(ii)%abun(1:60)
			enddo
			call SetupPartCloud(ii)
			do i=1,nr
				Cloud(ii)%Kref=Cloud(ii)%Kext(i,nlam+1)
				if(Cloud(ii)%type.eq."HALFGAUSS".and.P(i).gt.Cloud(ii)%P) then
				cloud_dens(i,ii)=(grav(i)*Cloud(ii)%tau*(Cloud(ii)%P**(Cloud(ii)%xi-2d0))/
     &					(Cloud(ii)%Kref*1d6*(P(i)**(Cloud(ii)%xi-1d0))*Cloud(ii)%dP*sqrt(2d0*pi)))*
     &					exp(-0.5d0*(log(P(i)/Cloud(ii)%P)/Cloud(ii)%dP)**8-0.5d0*Cloud(ii)%dP*(Cloud(ii)%xi-2d0)**2)
				else
				cloud_dens(i,ii)=(grav(i)*Cloud(ii)%tau*(Cloud(ii)%P**(Cloud(ii)%xi-2d0))/
     &					(Cloud(ii)%Kref*1d6*(P(i)**(Cloud(ii)%xi-1d0))*Cloud(ii)%dP*sqrt(2d0*pi)))*
     &					exp(-0.5d0*(log(P(i)/Cloud(ii)%P)/Cloud(ii)%dP)**2-0.5d0*Cloud(ii)%dP*(Cloud(ii)%xi-2d0)**2)
				endif
			enddo
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)
			tot=0d0
			do i=1,nr
				tot=tot+cloud_dens(i,ii)*Cloud(ii)%Kext(i,nlam+1)*(R(i+1)-R(i))
			enddo
			if(tot.gt.0d0) then
				cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*Cloud(ii)%tau/tot
			endif
		case("RING")
			Cloud(ii)%nlam=nlam+1
			Cloud(ii)%rv(1:nr)=Cloud(ii)%reff
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			Cloud(ii)%onepart=.true.
			do i=1,nr
				Cloud(ii)%frac(i,1:60)=Cloud(ii)%abun(1:60)
			enddo
			call SetupPartCloud(ii)
			Cloud(ii)%Kref=Cloud(ii)%Kext(1,nlam+1)
			cloud_dens(1:nr,ii)=0d0
			doRingCloud=.true.
		case default
			call output("Cloud type unknown: " // trim(Cloud(ii)%type))
			stop
	end select

	CloudMass=0d0
	do i=1,nr
		CloudMass=CloudMass+cloud_dens(i,ii)*4d0*pi*abs(R(i+1)**3-R(i)**3)/3d0
	enddo
	call output("Cloud mass: " // dbl2string(CloudMass/Mearth,'(es13.4E3)') // " Mearth")

	return
	end	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine SetupPartCloud(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,is,ilam,j
	real*8 phi,thet,tot,tot2,fact,tautot(nlam),HG,kabs,ksca
	logical computelamcloud(nlam),restrictcomputecloud
	real*8 dens1bar,haze_scale
	character*100 form

	if(.not.allocated(Cloud(ii)%rv)) allocate(Cloud(ii)%rv(nr))
	if(.not.allocated(Cloud(ii)%M)) allocate(Cloud(ii)%M(nr))
	if(.not.allocated(Cloud(ii)%Kabs)) then
		allocate(Cloud(ii)%Kabs(nr,nlam+1))
		allocate(Cloud(ii)%Ksca(nr,nlam+1))
		allocate(Cloud(ii)%Kext(nr,nlam+1))
		allocate(Cloud(ii)%g(nr,nlam+1))
		allocate(Cloud(ii)%F11(nr,nlam+1,180))
	endif

	select case(Cloud(ii)%opacitytype)
		case("OPACITY")
			call output("Using precomputed cloud particles")
			Cloud(ii)%Kext=0d0
			Cloud(ii)%Kabs=0d0
			Cloud(ii)%Ksca=0d0
			do j=1,Cloud(ii)%nmat
				do ilam=1,nlam+1
					do is=1,nr
						Cloud(ii)%Kext(is,ilam)=Cloud(ii)%Kext(is,ilam)+Cloud(ii)%frac(is,j)*Cloud(ii)%KeFile(j,ilam)
						Cloud(ii)%Kabs(is,ilam)=Cloud(ii)%Kabs(is,ilam)+Cloud(ii)%frac(is,j)*Cloud(ii)%KaFile(j,ilam)
						Cloud(ii)%Ksca(is,ilam)=Cloud(ii)%Ksca(is,ilam)+Cloud(ii)%frac(is,j)*Cloud(ii)%KsFile(j,ilam)
					enddo
				enddo
			enddo
			Cloud(ii)%g=Cloud(ii)%g0
			do j=1,180
c Henyey greenstein phase function
				Cloud(ii)%F11(1:nr,1:nlam+1,j)=(1d0-Cloud(ii)%g**2)/
     &				((1d0+Cloud(ii)%g**2-2d0*Cloud(ii)%g*cos(pi*(real(j)-0.5)/180d0))**(2d0/3d0))
			enddo
		case("FILE")
			call output("Using precomputed cloud particles")
			Cloud(ii)%Kext=0d0
			Cloud(ii)%Kabs=0d0
			Cloud(ii)%Ksca=0d0
			do j=1,Cloud(ii)%nmat
				do ilam=1,nlam+1
					do is=1,nr
						Cloud(ii)%Kext(is,ilam)=Cloud(ii)%Kext(is,ilam)+Cloud(ii)%frac(is,j)*Cloud(ii)%KeFile(j,ilam)
						Cloud(ii)%Kabs(is,ilam)=Cloud(ii)%Kabs(is,ilam)+Cloud(ii)%frac(is,j)*Cloud(ii)%KaFile(j,ilam)
						Cloud(ii)%Ksca(is,ilam)=Cloud(ii)%Ksca(is,ilam)+Cloud(ii)%frac(is,j)*Cloud(ii)%KsFile(j,ilam)
						Cloud(ii)%g(is,ilam)=Cloud(ii)%gFile(j,ilam)+Cloud(ii)%gFile(j,ilam)*Cloud(ii)%Ksca(is,ilam)
					enddo
				enddo
			enddo
			do ilam=1,nlam+1
				do is=1,nr
					if(Cloud(ii)%Ksca(is,ilam).gt.0d0) then
						Cloud(ii)%g(is,ilam)=Cloud(ii)%g(is,ilam)/Cloud(ii)%Ksca(is,ilam)
					else
						Cloud(ii)%g(is,ilam)=0d0
					endif
				enddo
			enddo
			do j=1,180
c Henyey greenstein phase function
				Cloud(ii)%F11(1:nr,1:nlam+1,j)=(1d0-Cloud(ii)%g**2)/
     &				((1d0+Cloud(ii)%g**2-2d0*Cloud(ii)%g*cos(pi*(real(j)-0.5)/180d0))**(2d0/3d0))
			enddo
		case("PARAMETERISED")
			call output("Computing parameterised cloud particles")
			do ilam=1,nlam
				Cloud(ii)%Kext(1:nr,ilam)=Cloud(ii)%kappa/(1d0+(lam(ilam)*1d4/Cloud(ii)%klam)**Cloud(ii)%kpow)
			enddo
			Cloud(ii)%Kext(1:nr,nlam+1)=Cloud(ii)%kappa/(1d0+(Cloud(ii)%lref/Cloud(ii)%klam)**Cloud(ii)%kpow)
			Cloud(ii)%Ksca(1:nr,1:nlam)=Cloud(ii)%Kext(1:nr,1:nlam)*Cloud(ii)%albedo
			Cloud(ii)%Kabs(1:nr,1:nlam)=Cloud(ii)%Kext(1:nr,1:nlam)*(1d0-Cloud(ii)%albedo)
			Cloud(ii)%g=Cloud(ii)%g0
			do j=1,180
c Henyey greenstein phase function
				Cloud(ii)%F11(1:nr,1:nlam+1,j)=(1d0-Cloud(ii)%g**2)/
     &				((1d0+Cloud(ii)%g**2-2d0*Cloud(ii)%g*cos(pi*(real(j)-0.5)/180d0))**(2d0/3d0))
c Rayleigh scattering phase function
c				Cloud(ii)%F11(1:nr,1:nlam+1,j)=3d0*(1d0+cos(pi*(real(j)-0.5)/180d0)**2)/4d0
			enddo
		case("MATERIAL","REFIND")
			call output("Computing inhomogeneous cloud particles")
			computelamcloud(1:nlam)=computelam(1:nlam)
			tautot=0d0
			restrictcomputecloud=(.not.computeT.and..not.scattering)
			if(Cloud(ii)%onepart) then
				is=1
				call ComputePart(Cloud(ii),ii,is,computelamcloud)
				do is=2,nr
					Cloud(ii)%Kext(is,1:nlam+1)=Cloud(ii)%Kext(1,1:nlam+1)
					Cloud(ii)%Kabs(is,1:nlam+1)=Cloud(ii)%Kabs(1,1:nlam+1)
					Cloud(ii)%Ksca(is,1:nlam+1)=Cloud(ii)%Ksca(1,1:nlam+1)
					Cloud(ii)%g(is,1:nlam+1)=Cloud(ii)%g(1,1:nlam+1)
					Cloud(ii)%F11(is,1:nlam+1,1:180)=Cloud(ii)%F11(1,1:nlam+1,1:180)
				enddo
			else
				do is=nr,1,-1
					call tellertje(nr-is+1,nr)
					call ComputePart(Cloud(ii),ii,is,computelamcloud)
					if(restrictcomputecloud) then
						do ilam=1,nlam
							if(computelamcloud(ilam)) then
								tautot(ilam)=tautot(ilam)+Cloud(ii)%Kext(is,ilam)*cloud_dens(is,ii)*(R(is+1)-R(is))
								if(tautot(ilam).gt.maxtau) computelamcloud(ilam)=.false.
							endif
						enddo
					endif
				enddo
			endif
			if(writefiles) then
				open(unit=28,file=trim(outputdir) // 'cloudcomp.dat',RECL=1000)				
				form='("#",a18,a23,' // trim(int2string(Cloud(ii)%nmat,'(i4)')) // 'a23)'
				write(28,form) "P[bar]","dens",(trim(Cloud(ii)%condensate(j)) // '[s]',j=1,Cloud(ii)%nmat)
				form='(es19.5,es23.5,' // trim(int2string(Cloud(ii)%nmat,'(i4)')) // 'f23.5)'
				do is=nr,1,-1
					write(28,form) P(is),cloud_dens(is,ii),Cloud(ii)%frac(is,1:Cloud(ii)%nmat)
				enddo
				close(unit=28)
			endif
		case default
			call output("Opacity type unknown: " // trim(Cloud(ii)%opacitytype))
			stop
	end select
		
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------



	subroutine SetupMaterialCloud()
	use GlobalSetup
	IMPLICIT NONE
	integer ii,i,j,iref,ngrid
	real*8 lgrid(nlam+1),e1d(nlam+1),e2d(nlam+1),kap
	logical lnkloglog
	external Carbon_BE_Zubko1996,Mg07Fe03SiO3_Dorschner1995,AstroSilicate
	external Enstatite_X,Enstatite_Y,Enstatite_Z,checkparticlefile
	external Forsterite_X,Forsterite_Y,Forsterite_Z
	external Rutile_xy,Rutile_z,Water,OrganicsHenning,Soot,Tholin
	external SiO,SiO2,Corrundum,Iron,FeO,Mg06Fe04O,MgO,SiC,H2SO4,AmorphSiO2

	lnkloglog=.true.
	ngrid=nlam+1

	do ii=1,nclouds
		lgrid(1:nlam)=lam(1:nlam)*1d4
		lgrid(nlam+1)=Cloud(ii)%lref
		iref=nlam+1
		do i=1,nlam
			if(lam(i)*1d4.gt.Cloud(ii)%lref) then
				lgrid(i+1:nlam+1)=lam(i:nlam)*1d4
				lgrid(i)=Cloud(ii)%lref
				iref=i
				exit
			endif
		enddo
		if(Cloud(ii)%opacitytype.eq."MATERIAL".or.Cloud(ii)%opacitytype.eq."REFIND") then
		allocate(Cloud(ii)%e1(Cloud(ii)%nmat,3,nlam+1),Cloud(ii)%e2(Cloud(ii)%nmat,3,nlam+1))
		if(Cloud(ii)%opacitytype.eq.'REFIND') then
			Cloud(ii)%e1(1,1,1:ngrid)=Cloud(ii)%e1_par
			Cloud(ii)%e2(1,1,1:ngrid)=Cloud(ii)%e2_par
			Cloud(ii)%nmat=1
			Cloud(ii)%nax(1)=1
		else
			do i=1,Cloud(ii)%nmat
				select case(Cloud(ii)%material(i))
					case("FILE")
						if(Cloud(ii)%lnkfile(i,2).eq.' ') then
							Cloud(ii)%nax(i)=1
						else
							Cloud(ii)%nax(i)=3
						endif
						do j=1,Cloud(ii)%nax(i)
							call readrefindCP(Cloud(ii)%lnkfile(i,j),lgrid,e1d,e2d,ngrid,lnkloglog)
							Cloud(ii)%e1(i,j,1:ngrid)=e1d(1:ngrid)
							Cloud(ii)%e2(i,j,1:ngrid)=e2d(1:ngrid)
						enddo
					case('ENSTATITE')
						Cloud(ii)%nax(i)=3
						call RegridDataLNK(Enstatite_X,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						call RegridDataLNK(Enstatite_Y,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,2,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,2,1:ngrid)=e2d(1:ngrid)
						call RegridDataLNK(Enstatite_Z,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,3,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,3,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=2.8
					case('FORSTERITE')
						Cloud(ii)%nax(i)=3
						call RegridDataLNK(Forsterite_X,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						call RegridDataLNK(Forsterite_Y,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,2,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,2,1:ngrid)=e2d(1:ngrid)
						call RegridDataLNK(Forsterite_Z,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,3,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,3,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=3.33
					case('BROOKITE','RUTILE','TiO2') 
						Cloud(ii)%nax(i)=3
						call RegridDataLNK(Rutile_xy,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						call RegridDataLNK(Rutile_xy,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,2,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,2,1:ngrid)=e2d(1:ngrid)
						call RegridDataLNK(Rutile_z,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,3,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,3,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=2.80
					case('WATER')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(Water,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=1.00
					case('CARBON')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(Carbon_BE_Zubko1996,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=1.80
					case('QUARTZ','SiO2')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(SiO2,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=2.648
					case('A-SiO2')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(AmorphSiO2,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=2.648
					case('SiO')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(SiO,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=2.18
					case('SiC')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(SiC,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=3.22
					case('IRON')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(Iron,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=7.87
					case('CORRUNDUM')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(Corrundum,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=3.97
					case('ORGANICS')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(OrganicsHenning,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=1.80
					case('THOLIN')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(Tholin,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=1.00
					case('ASTROSIL',"OLIVINE")
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(AstroSilicate,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=3.00
					case("PYROXENE")
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(Mg07Fe03SiO3_Dorschner1995,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=3.01
					case('H2SO4')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(H2SO4,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=1.00
					case('FeO')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(FeO,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=5.70
					case('MgO')
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(MgO,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=3.58
					case("SOOT")
						Cloud(ii)%nax(i)=1
						call RegridDataLNK(Soot,lgrid,e1d,e2d,ngrid,.true.)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=1.00
					case("optEC")
						Cloud(ii)%nax(i)=1
						call RefInd_optEC(lgrid,e1d,e2d,rad_optEC,Eg_optEC,ngrid)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=1.50
					case default
						call output("Material unknown: " // trim(Cloud(ii)%material(i)))
					stop
				end select
				do j=1,Cloud(ii)%nax(i)
					e1d(1:iref-1)=Cloud(ii)%e1(i,j,1:iref-1)
					e1d(iref:nlam)=Cloud(ii)%e1(i,j,iref+1:nlam+1)
					e1d(nlam+1)=Cloud(ii)%e1(i,j,iref)
					Cloud(ii)%e1(i,j,1:nlam+1)=e1d(1:nlam+1)
					e2d(1:iref-1)=Cloud(ii)%e2(i,j,1:iref-1)
					e2d(iref:nlam)=Cloud(ii)%e2(i,j,iref+1:nlam+1)
					e2d(nlam+1)=Cloud(ii)%e2(i,j,iref)
					Cloud(ii)%e2(i,j,1:nlam+1)=e2d(1:nlam+1)
				enddo
			enddo
		endif
		else if(Cloud(ii)%opacitytype.eq."OPACITY".or.Cloud(ii)%opacitytype.eq."FILE") then
			if(Cloud(ii)%opacitytype.eq."OPACITY") then
				j=0
			else
				j=1
			endif
			allocate(Cloud(ii)%KeFile(Cloud(ii)%nmat,ngrid))
			allocate(Cloud(ii)%KaFile(Cloud(ii)%nmat,ngrid))
			allocate(Cloud(ii)%KsFile(Cloud(ii)%nmat,ngrid))
			allocate(Cloud(ii)%gFile(Cloud(ii)%nmat,ngrid))
			do i=1,Cloud(ii)%nmat
				call regridOpacity(Cloud(ii)%material(i),lgrid,
     &	Cloud(ii)%KeFile(i,1:ngrid),Cloud(ii)%KaFile(i,1:ngrid),Cloud(ii)%KsFile(i,1:ngrid),
     &  Cloud(ii)%gFile(i,1:ngrid),ngrid,j)
				kap=Cloud(ii)%KeFile(i,iref)
				Cloud(ii)%KeFile(i,iref:nlam)=Cloud(ii)%KeFile(i,iref+1:nlam+1)
				Cloud(ii)%KeFile(i,nlam+1)=kap
				kap=Cloud(ii)%KaFile(i,iref)
				Cloud(ii)%KaFile(i,iref:nlam)=Cloud(ii)%KaFile(i,iref+1:nlam+1)
				Cloud(ii)%KaFile(i,nlam+1)=kap
				kap=Cloud(ii)%KsFile(i,iref)
				Cloud(ii)%KsFile(i,iref:nlam)=Cloud(ii)%KsFile(i,iref+1:nlam+1)
				Cloud(ii)%KsFile(i,nlam+1)=kap
				kap=Cloud(ii)%gFile(i,iref)
				Cloud(ii)%gFile(i,iref:nlam)=Cloud(ii)%gFile(i,iref+1:nlam+1)
				Cloud(ii)%gFile(i,nlam+1)=kap
			enddo
		endif
	enddo


	return
	end
	


	subroutine RefreshMaterialCloud()
	use GlobalSetup
	IMPLICIT NONE
	integer ii,i,j,iref,ngrid
	real*8 lgrid(nlam+1),e1d(nlam+1),e2d(nlam+1),kap
	logical lnkloglog
	external Carbon_BE_Zubko1996,Mg07Fe03SiO3_Dorschner1995,AstroSilicate
	external Enstatite_X,Enstatite_Y,Enstatite_Z,checkparticlefile
	external Forsterite_X,Forsterite_Y,Forsterite_Z
	external Rutile_xy,Rutile_z,Water,OrganicsHenning,Soot,Tholin
	external SiO,SiO2,Corrundum,Iron,FeO,Mg06Fe04O,MgO,SiC,H2SO4,AmorphSiO2

	lnkloglog=.true.
	ngrid=nlam+1

	do ii=1,nclouds
		lgrid(1:nlam)=lam(1:nlam)*1d4
		lgrid(nlam+1)=Cloud(ii)%lref
		iref=nlam+1
		do i=1,nlam
			if(lam(i)*1d4.gt.Cloud(ii)%lref) then
				lgrid(i+1:nlam+1)=lam(i:nlam)*1d4
				lgrid(i)=Cloud(ii)%lref
				iref=i
				exit
			endif
		enddo
		if(Cloud(ii)%opacitytype.eq."MATERIAL".or.Cloud(ii)%opacitytype.eq."REFIND") then
		if(Cloud(ii)%opacitytype.eq.'REFIND') then
			Cloud(ii)%e1(1,1,1:ngrid)=Cloud(ii)%e1_par
			Cloud(ii)%e2(1,1,1:ngrid)=Cloud(ii)%e2_par
			Cloud(ii)%nmat=1
			Cloud(ii)%nax(1)=1
		else
			do i=1,Cloud(ii)%nmat
				select case(Cloud(ii)%material(i))
					case("optEC")
						Cloud(ii)%nax(i)=1
						call RefInd_optEC(lgrid,e1d,e2d,rad_optEC,Eg_optEC,ngrid)
						Cloud(ii)%e1(i,1,1:ngrid)=e1d(1:ngrid)
						Cloud(ii)%e2(i,1,1:ngrid)=e2d(1:ngrid)
						Cloud(ii)%rho_mat(i)=1.50
						do j=1,Cloud(ii)%nax(i)
							e1d(1:iref-1)=Cloud(ii)%e1(i,j,1:iref-1)
							e1d(iref:nlam)=Cloud(ii)%e1(i,j,iref+1:nlam+1)
							e1d(nlam+1)=Cloud(ii)%e1(i,j,iref)
							Cloud(ii)%e1(i,j,1:nlam+1)=e1d(1:nlam+1)
							e2d(1:iref-1)=Cloud(ii)%e2(i,j,1:iref-1)
							e2d(iref:nlam)=Cloud(ii)%e2(i,j,iref+1:nlam+1)
							e2d(nlam+1)=Cloud(ii)%e2(i,j,iref)
							Cloud(ii)%e2(i,j,1:nlam+1)=e2d(1:nlam+1)
						enddo
				end select
			enddo
		endif
		endif
	enddo


	return
	end
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine regridOpacity(input,grid,Ke,Ka,Ks,g,n,filetype)
	IMPLICIT NONE
	integer i,j,n,filetype
	real*8 grid(n),Ke(n),Ka(n),Ks(n),x0,Ke0,Ka0,Ks0,x1,Ke1,Ka1,Ks1
	real*8 g(n),g0,g1
	character*500 input
	logical truefalse
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,FORM="FORMATTED",ACCESS="STREAM")
	i=1
1	if(filetype.eq.0) then
		read(20,*,end=102,err=1) x0,Ke0,Ka0,Ks0
		g0=0d0
	else
		read(20,*,end=102,err=1) x0,Ke0,Ks0,g0
		Ke0=Ke0*10d0
		Ks0=Ke0*Ks0
		Ka0=Ke0-Ks0
	endif
103	if(x0.ge.grid(i)) then
		Ke(i)=Ke0
		Ka(i)=Ka0
		Ks(i)=Ks0
		g(i)=g0
		i=i+1
		goto 103
	endif
100	if(filetype.eq.0) then
		read(20,*,end=102,err=100) x1,Ke1,Ka1,Ks1
		g1=0d0
	else
		read(20,*,end=102,err=100) x1,Ke1,Ks1,g1
		Ke1=Ke1*10d0
		Ks1=Ke1*Ks1
		Ka1=Ke1-Ks1
	endif
101	if(grid(i).le.x1.and.grid(i).ge.x0) then
		Ke(i)=Ke1+(grid(i)-x1)*(Ke0-Ke1)/(x0-x1)
		Ka(i)=Ka1+(grid(i)-x1)*(Ka0-Ka1)/(x0-x1)
		Ks(i)=Ks1+(grid(i)-x1)*(Ks0-Ks1)/(x0-x1)
		g(i)=g1+(grid(i)-x1)*(g0-g1)/(x0-x1)
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	Ke0=Ke1
	Ka0=Ka1
	Ks0=Ks1
	g0=g1
	goto 100
102	continue
	do j=i,n
		Ka(j)=Ka(i-1)*(grid(i-1)/grid(j))**2
		Ks(j)=Ks(i-1)*(grid(i-1)/grid(j))**4
		Ke(j)=Ka(j)+Ks(j)
		g(j)=g(i-1)
	enddo
	close(unit=20)
	return
	end


	subroutine AddRingCloud()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,i
	logical doR(nclouds)
	real*8 tau

	doR=.false.
	do ii=1,nclouds
		if(Cloud(ii)%type.eq."RING") doR(ii)=.true.
	enddo
	
	do i=1,nlam
		tau=0d0
		do ii=1,nclouds
			if(doR(ii)) tau=tau+Cloud(ii)%tau*Cloud(ii)%Kext(1,i)/Cloud(ii)%Kref
		enddo
		obsA(0,i)=obsA(0,i)+pi*((Rring+dRring)**2-Rring**2)*Rplanet**2*(1d0-exp(-tau))
	enddo
	
	return
	end
	
