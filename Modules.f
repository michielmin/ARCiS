c=========================================================================================
c module containing the physical constants in cgs (moet nog naar SI)
c=========================================================================================
	module Constants
	IMPLICIT NONE
	real*8 pi,Ggrav,Msun,AU,clight,Rsun,mp,kb,hplanck,parsec,Lsun,sigma
	real*8 Mearth,Rearth,Mjup,Rjup,year,micron,Rgas,Avogadro,atm
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(clight=2.9979245800d10) !cm/s
	parameter(AU=1.49598d13)
	parameter(parsec=3.08568025d18)
	parameter(Rsun=6.955d10)
	parameter(Msun=1.98892d33)
	parameter(Lsun=3.827d33)
	parameter(kb=1.3806503d-16)
	parameter(sigma=5.6704d-5)
c	parameter(mp=1.67262178d-24)	!proton mass
	parameter(mp=1.660540210d-24)	!atomic mass unit
	parameter(Ggrav=6.6725985d-8) ! in cm^3/g/s^2
	parameter(hplanck=6.626068d-27) ! cm^2 g/s
	parameter(Mearth=5.97219d27)
	parameter(Mjup=1.898d30)
	parameter(Rearth=6.3781d8)
	parameter(Rjup=6.9911e9)		!7.1492d9)
	parameter(year=24d0*60d0*60d0*265.25d0)
	parameter(micron=1d-4)
	parameter(Rgas=8.3144621e7)
	parameter(Avogadro=6.022136736e23)
	parameter(atm=1.01325)	! bar
	
	end module Constants

c=========================================================================================
c global setup for the code
c=========================================================================================
	module GlobalSetup
	IMPLICIT NONE
	real*8 Mplanet,Rplanet									! mass and radius of the planet
	real*8 Tstar,Rstar,logg,Dplanet
	real*8,allocatable :: dens(:),T(:),P(:),Ndens(:),Tin(:)	! radius
	real*8,allocatable :: dust_dens(:,:)					! radius, component
	real*8,allocatable :: R(:),Hp(:)						! radius
	real*8,allocatable :: mixrat(:)							! component
	real*8,allocatable :: mixrat_r(:,:),mixrat_old_r(:,:)	! radius,component
	real*8,allocatable :: Cabs(:,:,:),Csca(:,:)				! radius,wav,g
	real*8,allocatable :: cloud_dens(:,:)					! radius, cloud
	real*8,allocatable :: Fstar(:)							! wavelength
	real*8,allocatable :: tau1depth(:,:),cloudtau(:,:)		! ncc,wavelength
	integer nangle_Jscat
	parameter(nangle_Jscat=60)
	real*8,allocatable :: Jscat(:,:)						! radius, angle
	integer nT,np,nr,nmol,nlam		! #T, #P, #radial points, #molecules, #wavelength bins, #obs
	integer nlines,ng,ncia,nclouds
	character*500 outputdir,HITRANdir,HITEMPdir
	integer idum,maxiter,Nphot0
!$OMP THREADPRIVATE(idum)
	logical retrieval,outputopacity,do_cia,gridTPfile,scattering,scattstar,computeT
	logical dochemistry,retrieve_profile,condensates
	logical,allocatable :: includemol(:)
	real*8 lam1,lam2,specres,Pmin,Pmax,epsCk,distance,TP0,dTP,TeffP
	real*8 gammaT1,gammaT2,kappaT,betaT,alphaT
	logical mixratfile,par_tprofile
	character*500 TPfile
	real*8 metallicity,COratio,PQ,mixP,PRplanet,mixratHaze
	real*8 cutoff_abs,cutoff_lor,eps_lines,maxtau,factRW
	real*8,allocatable :: lam(:),freq(:),dfreq(:)
	real*8,allocatable :: ZZ(:,:,:),TZ(:)	! partition function
	integer nTZ,nspike
	integer,allocatable :: niso(:)
	real*8 Mmol(59),mu
	integer nBB
	parameter(nBB=10000)
	real*8,allocatable :: BB(:,:)						! nBB,nlam
	character*10 molname(59)
	parameter(molname = (/'H2O   ','CO2   ','O3    ','N2O   ','CO    ','CH4   ',
     &	'O2    ','NO    ','SO2   ','NO2   ','NH3   ','HNO3  ','OH    ','HF    ',
     &	'HCl   ','HBr   ','HI    ','ClO   ','OCS   ','H2CO  ','HOCl  ','N2    ',
     &	'HCN   ','CH3Cl ','H2O2  ','C2H2  ','C2H6  ','PH3   ','COF2  ','SF6   ',
     &	'H2S   ','HCOOH ','HO2   ','O     ','ClONO2','NO+   ','HOBr  ','C2H4  ',
     &	'CH3OH ','CH3Br ','CH3CN ','CF4   ','C4H2  ','HC3N  ','H2    ','CS    ',
     &	'SO3   ','He    ','X     ','X     ','X     ','X     ','X     ','X     ',
     &  'X     ','Na    ','K     ','TiO   ','VO    ' /))
	parameter(Mmol = (/     17.8851,  43.6918,  47.6511,  43.6947,  27.8081,  15.9272,  
     &	31.7674,  29.7889,  63.5840,  45.6607,  16.9072,  62.5442,  16.8841,  19.8619,  
     &	36.1973,  80.3271,   1.0070,  51.0760,  59.6379,  29.8088,  52.0765,  27.8112,  
     &	26.8304,  50.1116,  33.7598,  25.8499,  29.8518,  33.7516,  65.5261, 144.9081,  
     &	33.8332,  45.6731,  32.7593,  15.8794,  96.7366,  29.7813,  96.2063,  27.8507,  
     &	31.7949,  94.2413,  40.7302,  87.3580,  49.6543,  50.6424,   2.0014,  43.7539,  
     &	79.3792,   4.0030,   0.0000,   0.0000,   0.0000,   0.0000,   0.0000,   0.0000,
     &   0.0000,  22.9900,  39.0980,  63.8660,  66.9410 /))
	real*8,allocatable :: a_therm(:),a_press(:)
	integer n_voigt
	logical HITEMP,opacitymode,compute_opac
	integer nPom,nTom
	character*500 opacitydir,specresfile
	real*8 Tmin,Tmax,minTprofile,maxTprofile
	real*8 sintheta(360),costheta(360)
	logical,allocatable :: do_dB(:)
	
	logical sinkZ
	real*8 alphaZ

	real*8,allocatable :: flux(:,:),obsA(:,:),phase(:,:,:)
	integer ncc,nphase
	logical cloudcompute
	logical,allocatable :: docloud(:,:)
	real*8,allocatable :: cloudfrac(:),XCloud(:,:),XeqCloud(:,:),XeqCloud_old(:,:)

	integer,allocatable,dimension(:) :: L_imol,L_iiso,L_nclose,L_ilam
	real*8,allocatable,dimension(:) :: L_Aul,L_freq,L_Elow,L_lam,L_S0,L_S
	real*8,allocatable,dimension(:) :: L_gamma_air,L_gamma_self
	real*8,allocatable,dimension(:) :: L_a_therm,L_a_press,L_n
	real*8,allocatable,dimension(:) :: L_gu,L_gl,L_Saver
	logical,allocatable,dimension(:) :: L_do
	
	integer,allocatable :: ig_comp(:,:,:)

	type CIA_pair
		character*20 name
		character*500 filename
		real*8,allocatable :: T(:),Cabs(:,:)
		integer nT,imol1,imol2
	end type CIA_pair
	
	type(CIA_pair),allocatable :: CIA(:)
	real*8 cia_mixrat(59)

	type Mueller
		real*8 F11(180),F12(180),F22(180)
		real*8 F33(180),F44(180),F34(180)
		real*8 IF11(180),IF12(180)
	end type Mueller
	
	type(Mueller) Rayleigh

	type CloudType
		real*8 P,dP,s,column
		real*8 coverage,frain
		real*8,allocatable :: rv(:),w(:)							! dimension nsize
		real*8 rho,amin,amax,fmax,porosity,fcarbon,reff,veff
		logical blend,haze
		real*8 fcond,mixrat
		real*8,allocatable :: Kabs(:,:),Ksca(:,:),Kext(:,:)			! dimension nsize,nlam
		type(Mueller),allocatable :: F(:,:)							! dimension nsize,nlam
		character*500 file
		character*20 standard,ptype
		character*500 species
		integer nsize,nsubgrains
	end type CloudType

	type(CloudType),allocatable :: Cloud(:) 

cPoints for the temperature structure
	real*8,allocatable :: P_point(:),T_point(:)
	integer n_points

	type RetrievalPar
		character*500 keyword
		real*8 xmin,xmax,x0,dx,value,error1,error2
		logical logscale,squarescale,opacitycomp
		integer n
	end type RetrievalPar

	type(RetrievalPar),allocatable :: RetPar(:)
	integer n_ret
	
	type ObservedSpec
		character*500 file
		character*10 type
		real*8,allocatable :: lam(:),y(:),dy(:),R(:),Rexp(:),model(:)
		real*8 beta,scale
		integer nlam
		logical spec
	end type ObservedSpec
	type(ObservedSpec),allocatable :: ObsSpec(:)

	integer npop,ngen,nobs
	logical gene_cross
	

c========================================================
c Interfaces for input/output subroutines
c========================================================

	interface
		subroutine outputform(string,form)
			character :: string*(*)
			character,intent(in),optional :: form*(*)
		end subroutine outputform
	end interface

	interface
		function int2string(i,form)
			character*20 :: int2string
			integer :: i
			character,intent(in),optional :: form*(*)
		end function int2string
	end interface

	interface
		function dbl2string(x,form)
			character*20 :: dbl2string
			real*8 :: x
			character,intent(in),optional :: form*(*)
		end function dbl2string
	end interface


	end module GlobalSetup
	



