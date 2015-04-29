c=========================================================================================
c module containing the physical constants in cgs (moet nog naar SI)
c=========================================================================================
	module Constants
	IMPLICIT NONE
	real*8 pi,G,Msun,AU,clight,Rsun,mp,kb,hplanck,parsec,Lsun,sigma
	real*8 Mearth,Rearth,Mjup,Rjup,year
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(clight=2.9979245800d10) !cm/s
	parameter(AU=1.49598d13)
	parameter(parsec=3.08568025d18)
	parameter(Rsun=6.955d10)
	parameter(Msun=1.98892d33)
	parameter(Lsun=3.827d33)
	parameter(kb=1.3806503d-16)
	parameter(sigma=5.6704d-5)
	parameter(mp=1.67262178d-24)	!proton mass
	parameter(G=6.67300d-8) ! in cm^3/g/s^2
	parameter(hplanck=6.626068d-27) ! cm^2 g/s
	parameter(Mearth=5.97219d27)
	parameter(Mjup=1.89813d30)
	parameter(Rearth=6.3781d8)
	parameter(Rjup=7.1492d9)
	parameter(year=24d0*60d0*60d0*265.25d0)
	
	end module Constants

c=========================================================================================
c global setup for the code
c=========================================================================================
	module GlobalSetup
	IMPLICIT NONE
	real*8 Mplanet,Rplanet									! mass and radius of the planet
	real*8,allocatable :: dens(:),T(:),P(:)					! radius
	real*8,allocatable :: gas_dens(:,:),dust_dens(:,:)		! radius, component
	real*8,allocatable :: R(:)								! radius
	real*8,allocatable :: opac(:,:,:,:)						! component,wav,T,P
	integer nT,np,nrad,ncomp,nlam,nobs		! #T, #P, #radial points, #components, #wavelength bins, #obs
	character*500 outputdir

c string converting functions
	character*20 int2string,dbl2string
	external int2string,dbl2string

	logical retrieval

	type SettingKey
		character*100 key1,key2,value
		integer nr1,nr2
		logical last
		type(SettingKey),pointer :: next
	end type SettingKey

	type Observation
		character*500 filename
		character*10 type
		real*8,allocatable :: lam(:),flux(:)
	end type Observation

	type(Observation),allocatable :: obs(:)

	end module GlobalSetup
	
