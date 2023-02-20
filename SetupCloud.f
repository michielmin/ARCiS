	subroutine SetupCloud(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 pp,tot,column,tt,z,densdust,CloudMass
	integer ii,i,j,nsubr,isize,ilam
	real*8 Xc,Xc1,lambdaC,Ca,Cs,tau,P_SI1,P_SI2,veff,frac(nr,10)
	logical cl

	if(.not.allocated(Cloud(ii)%rv)) then
		allocate(Cloud(ii)%rv(nr))
		allocate(Cloud(ii)%sigma(nr))
	endif
	
	select case(Cloud(ii)%type)
		case("DIFFUSE")
			Cloud(ii)%nlam=nlam
			call DiffuseCloud(ii)
			Cloud(ii)%rv=Cloud(ii)%rv*1d4
			Cloud(ii)%sigma=1d-10
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
			call SetupPartCloud(ii)
		case("LAYER")
			Cloud(ii)%nlam=nlam+1
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			Cloud(ii)%onepart=(Cloud(ii)%rpow.eq.0d0)
			do i=1,nr
				Cloud(i)%frac(i,1:40)=Cloud(i)%abun(1:40)
			enddo
			call SetupPartCloud(ii)
			do i=1,nr
				Cloud(ii)%Kref=Cloud(ii)%Kext(i,nlam+1)
				if(P(i).gt.Cloud(ii)%Pmin.and.P(i).lt.Cloud(ii)%Pmax) then
					cloud_dens(i,ii)=(grav(i)*Cloud(ii)%xi*Cloud(ii)%tau*P(i)**(Cloud(ii)%xi-1d0))/
     &					(Cloud(ii)%Kref*1d6*(Cloud(ii)%Ptau**Cloud(ii)%xi-Cloud(ii)%Pmin**Cloud(ii)%xi))
				else
					cloud_dens(i,ii)=0d0
				endif
			enddo
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)
		case("SLAB")
			Cloud(ii)%nlam=nlam+1
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			Cloud(ii)%onepart=(Cloud(ii)%rpow.eq.0d0)
			do i=1,nr
				Cloud(i)%frac(i,1:40)=Cloud(i)%abun(1:40)
			enddo
			call SetupPartCloud(ii)
			do i=1,nr
				Cloud(ii)%Kref=Cloud(ii)%Kext(i,nlam+1)
				if(P(i).gt.Cloud(ii)%Pmin.and.P(i).lt.Cloud(ii)%Pmax) then
					cloud_dens(i,ii)=(grav(i)*2d0*Cloud(ii)%tau*P(i))/
     &					(Cloud(ii)%Kref*1d6*(Cloud(ii)%Pmax*Cloud(ii)%Pmax-Cloud(ii)%Pmin*Cloud(ii)%Pmin))
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
				Cloud(i)%frac(i,1:40)=Cloud(i)%abun(1:40)
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
			cloud_dens(1:nr,ii)=dens(1:nr)*Cloud(ii)%mixrat
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			Cloud(ii)%onepart=(Cloud(ii)%rpow.eq.0d0)
			do i=1,nr
				Cloud(i)%frac(i,1:40)=Cloud(i)%abun(1:40)
			enddo
			call SetupPartCloud(ii)
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

	if(.not.allocated(Cloud(ii)%rv)) allocate(Cloud(ii)%rv(nr))
	if(.not.allocated(Cloud(ii)%M)) allocate(Cloud(ii)%M(nr))
	if(.not.allocated(Cloud(ii)%Kabs)) then
		allocate(Cloud(ii)%Kabs(nr,nlam+1))
		allocate(Cloud(ii)%Ksca(nr,nlam+1))
		allocate(Cloud(ii)%Kext(nr,nlam+1))
	endif

	call output("Computing inhomogeneous cloud particles")

	select case(Cloud(ii)%opacitytype)
		case("PARAMETERISED")
			do ilam=1,nlam
				Cloud(ii)%Kext(1:nr,ilam)=Cloud(ii)%kappa/(1d0+(lam(ilam)*1d4/Cloud(ii)%klam)**Cloud(ii)%kpow)
			enddo
			Cloud(ii)%Kext(1:nr,nlam+1)=Cloud(ii)%kappa/(1d0+(Cloud(ii)%lref/Cloud(ii)%klam)**Cloud(ii)%kpow)
			Cloud(ii)%Ksca(1:nr,1:nlam)=Cloud(ii)%Kext(1:nr,1:nlam)*Cloud(ii)%albedo
			Cloud(ii)%Kabs(1:nr,1:nlam)=Cloud(ii)%Kext(1:nr,1:nlam)*(1d0-Cloud(ii)%albedo)
		case("MATERIAL","REFIND")
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
		case default
			call output("Opacity type unknown: " // trim(Cloud(ii)%opacitytype))
			stop
	end select
		
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------



	subroutine SetupRefIndCloud()
	use GlobalSetup
	IMPLICIT NONE
	integer ii,i,j,iref,ngrid
	real*8 lgrid(nlam+1),e1d(nlam+1),e2d(nlam+1)
	logical lnkloglog
	external Carbon_BE_Zubko1996,Mg07Fe03SiO3_Dorschner1995,AstroSilicate
	external Enstatite_X,Enstatite_Y,Enstatite_Z,checkparticlefile
	external Forsterite_X,Forsterite_Y,Forsterite_Z
	external Rutile_xy,Rutile_z,Water,OrganicsHenning,Soot,Tholin
	external SiO,SiO2,Corrundum,Iron,FeO,Mg06Fe04O,MgO,SiC,H2SO4

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
	enddo


	return
	end
	


