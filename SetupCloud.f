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
			call DiffuseCloud(ii)
			Cloud(ii)%rv=Cloud(ii)%rv*1d4
			Cloud(ii)%sigma=1d-10
			call SetupPartCloud(ii)
		case("FILE")
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
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			call SetupPartCloud(ii)
			do i=1,nr
				if(P(i).gt.Cloud(ii)%Pmin.and.P(i).lt.Cloud(ii)%Pmax) then
					cloud_dens(i,ii)=(grav(i)*Cloud(ii)%xi*Cloud(ii)%tau*P(i)**(Cloud(ii)%xi-1d0))/
     &					(Cloud(ii)%Kref*(Cloud(ii)%Ptau**Cloud(ii)%xi-Cloud(ii)%Pmin**Cloud(ii)%xi))
				else
					cloud_dens(i,ii)=0d0
				endif
			enddo
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)
		case("SLAB")
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			call SetupPartCloud(ii)
			do i=1,nr
				if(P(i).gt.Cloud(ii)%Pmin.and.P(i).lt.Cloud(ii)%Pmax) then
					cloud_dens(i,ii)=(grav(i)*2d0*Cloud(ii)%tau*P(i))/
     &					(Cloud(ii)%Kref*(Cloud(ii)%Pmax*Cloud(ii)%Pmax-Cloud(ii)%Pmin*Cloud(ii)%Pmin))
				else
					cloud_dens(i,ii)=0d0
				endif
			enddo
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)
		case("DECK")
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
			call SetupPartCloud(ii)
			do i=1,nr
				cloud_dens(i,ii)=(grav(i)*exp((P(i)-Cloud(ii)%Ptau)/Cloud(ii)%Phi))/
     &				(Cloud(ii)%Kref*Cloud(ii)%Phi*(1d0-exp(-Cloud(ii)%Ptau/Cloud(ii)%Phi)))
			enddo
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)
		case("HOMOGENEOUS")
			cloud_dens(1:nr,ii)=dens(1:nr)*Cloud(ii)%mixrat
			Cloud(ii)%rv(1:nr)=Cloud(ii)%rnuc+(Cloud(ii)%reff-Cloud(ii)%rnuc)*(P(1:nr)/Cloud(ii)%Pref)**Cloud(ii)%rpow
			Cloud(ii)%sigma(1:nr)=Cloud(ii)%veff
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

	if(.not.allocated(Cloud(ii)%rv)) allocate(Cloud(ii)%rv(Cloud(ii)%nr))
	if(.not.allocated(Cloud(ii)%M)) allocate(Cloud(ii)%M(Cloud(ii)%nr))
	if(.not.allocated(Cloud(ii)%Kabs)) then
		allocate(Cloud(ii)%Kabs(Cloud(ii)%nr,nlam))
		allocate(Cloud(ii)%Ksca(Cloud(ii)%nr,nlam))
		allocate(Cloud(ii)%Kext(Cloud(ii)%nr,nlam))
	endif

	call output("Computing inhomogeneous cloud particles")

	j=0
	veff=Cloud(ii)%veff
1	tot=0d0
	do i=1,Cloud(ii)%nr
		Cloud(ii)%w(i)=(1d4*Cloud(ii)%rv(i))**(1d0+(1d0-3d0*veff)/veff)*
     &					exp(-1d0*Cloud(ii)%rv(i)/(Cloud(ii)%reff*veff))
		Cloud(ii)%w(i)=Cloud(ii)%w(i)*Cloud(ii)%rv(i)**3
		tot=tot+Cloud(ii)%w(i)
	enddo
	if(.not.tot.gt.0d0) then
		if(j.lt.10) then
			veff=veff*2d0
			call output("increasing veff")
			j=j+1
			goto 1
		else
			tot=1d0
			Cloud(ii)%w(1:Cloud(ii)%nr)=1d0
		endif
	endif
	Cloud(ii)%w=Cloud(ii)%w/tot

	select case(Cloud(ii)%opacitytype)
		case("PARAMETERISED")
			do ilam=1,nlam
				Cloud(ii)%Kext(1:nr,ilam)=Cloud(ii)%kappa/(1d0+(lam(ilam)/Cloud(ii)%klam)**Cloud(ii)%kpow)
			enddo
			Cloud(ii)%Kref=Cloud(ii)%kappa/(1d0+(Cloud(ii)%lref/Cloud(ii)%klam)**Cloud(ii)%kpow)
			Cloud(ii)%Ksca(1:nr,1:nlam)=Cloud(ii)%Kext(1:nr,1:nlam)*Cloud(ii)%albedo
			Cloud(ii)%Kabs(1:nr,1:nlam)=Cloud(ii)%Kext(1:nr,1:nlam)*(1d0-Cloud(ii)%albedo)
		case("MATERIAL","REFIND")
			computelamcloud(1:nlam)=computelam(1:nlam)
			tautot=0d0
			restrictcomputecloud=(.not.computeT.and..not.scattering)
			do is=Cloud(ii)%nr,1,-1
				call tellertje(Cloud(ii)%nr-is+1,Cloud(ii)%nr)
				call ComputePart(Cloud(ii),ii,is,computelamcloud)
				if(restrictcomputecloud.and.(useDRIFT.or.cloudcompute)) then
					do ilam=1,nlam
						if(computelamcloud(ilam)) then
							tautot(ilam)=tautot(ilam)+Cloud(ii)%Kext(is,ilam)*cloud_dens(is,ii)*(R(is+1)-R(is))
							if(tautot(ilam).gt.maxtau) computelamcloud(ilam)=.false.
						endif
					enddo
				endif
			enddo
		case default
			call output("Opacity type unknown: " // trim(Cloud(ii)%opacitytype))
			stop
	end select
		
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


