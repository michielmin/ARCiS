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
	parameter(mp=1.660539040d-24)	!atomic mass unit
c	parameter(Ggrav=6.6725985d-8) ! in cm^3/g/s^2
	parameter(Ggrav=6.6740831d-8) ! in cm^3/g/s^2
	parameter(hplanck=6.626068d-27) ! cm^2 g/s
	parameter(Mearth=5.97219d27)
	parameter(Mjup=1.898d30)
	parameter(Rearth=6.3781d8)
	parameter(Rjup=6.9911e9)
c	parameter(Rjup=7.1492d9)
	parameter(year=24d0*60d0*60d0*365.25d0)
	parameter(micron=1d-4)
	parameter(Rgas=8.3144621e7)
	parameter(Avogadro=6.022136736e23)
	parameter(atm=1.01325)	! bar
	
	end module Constants



c========================================================
c Interfaces for input/output subroutines
c========================================================
	module OutputModule
	IMPLICIT NONE

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


	end module OutputModule
	


c=========================================================================================
c global setup for the code
c=========================================================================================
	module GlobalSetup
	use OutputModule
	IMPLICIT NONE
	real*8 Mplanet,Rplanet,Pplanet,loggPlanet				! mass and radius of the planet at pressure Pplanet
	real*8 Tstar,Rstar,Mstar,logg,Dplanet
	real*8 orbit_P,orbit_e,orbit_omega,orbit_inc
	real*8,allocatable :: dens(:),T(:),P(:),Ndens(:),Tin(:)	! radius
	real*8,allocatable :: dust_dens(:,:)					! radius, component
	real*8,allocatable :: R(:),Hp(:)						! radius
	real*8,allocatable :: mixrat(:)							! component
	real*8,allocatable :: mixrat_r(:,:),mixrat_old_r(:,:)	! radius,component
	real*8,allocatable :: Cabs(:,:,:),Csca(:,:)				! radius,wav,g
	real*8,allocatable :: cloud_dens(:,:)					! radius, cloud
	real*8,allocatable :: Fstar(:)							! wavelength
	real*8,allocatable :: tau1depth(:,:),cloudtau(:,:)		! ncc,wavelength
	real*8,allocatable :: Cabs_mol(:,:,:,:),Cext_cont(:,:)
	integer nangle_Jscat
	parameter(nangle_Jscat=60)
	real*8,allocatable :: Jscat(:,:)						! radius, angle
	integer nT,np,nr,nmol,nlam		! #T, #P, #radial points, #molecules, #wavelength bins, #obs
	integer nlines,ng,ncia,nclouds,nTiter
	character*500 outputdir,HITRANdir,HITEMPdir
	integer idum,maxiter,Nphot0,idum0,nlamdust,iWolk
!$OMP THREADPRIVATE(idum)
	logical retrieval,outputopacity,do_cia,gridTPfile,scattering,scattstar,computeT,computecontrib
	logical dochemistry,retrieve_profile,condensates,faircoverage,speclimits,mapCOratio,randomseed
	logical,allocatable :: includemol(:),didcondens(:),lamemis(:),lamtrans(:),opacitymol(:)
	real*8 lam1,lam2,specres,Pmin,Pmax,epsCk,distance,TP0,dTP,TeffP,specresdust,twind,epsiter
	real*8 gammaT1,gammaT2,kappaT,betaT,alphaT,Tchem,Pchem,Psimplecloud,metallicity0
	logical mixratfile,par_tprofile,adiabatic_tprofile,domakeai,modelsucces,PTchemAbun,useobsgrid
	logical didcondens_chem,coagulation,resume_multinest,doMCcompute,disequilibrium,const_eff_multinest
	character*500 TPfile,particledir,retrievaltype,planetparameterfile,planetname,element_abun_file
	real*8 metallicity,COratio,PQ,mixP,PRplanet,maxchemtime,TiScale,f_multinest,tol_multinest
	real*8 mixratHaze,PHaze,dPHaze,kappaHaze,Kzz,SiOratio,NOratio,fDay,betapow,Kxx,Kyy,vxx,powvxx,night2day
	real*8 Kzz_deep,Kzz_1bar,Kzz_P,SOratio,Tsurface
	logical enhancecarbon,fast_chem,gamma_equal,dopostequalweights
	logical transspec,emisspec,rainout,computeLC,doscaleR
	real*8 cutoff_abs,cutoff_lor,eps_lines,maxtau,factRW,Tform,Pform,f_dry,f_wet,scale_fe
	real*8,allocatable :: lam(:),freq(:),dfreq(:),lamdust(:),dlam(:),surface_emis(:)
	real*8,allocatable :: gg(:),wgg(:),obsA_contr(:,:),flux_contr(:,:),obsA_LC(:,:)
	real*8,allocatable :: ZZ(:,:,:),TZ(:)	! partition function
	real*8 planetform_fdust,planetform_fplan,planetform_Rstart,planetform_Mstart
	logical planetform,massprior
	real*8 Mp_prior,dMp_prior
	character*20 surfacetype
	integer nTZ,nspike,nai,nboot,npew
	integer,allocatable :: niso(:),instr_nobs(:)
	real*8,allocatable :: MMW(:)
	real*8,allocatable :: PTaverage3D(:,:),mixrat_average3D(:,:,:)
	logical fulloutput3D
	integer nBB
	parameter(nBB=10000)
	real*8,allocatable :: BB(:,:)						! nBB,nlam

	logical do3D
	real*8 Kzz3D_1,Kzz3D_2,beta3D_1,beta3D_2,Sdot3D_1,Sdot3D_2,long_shift
	real*8 kappa3D_1,kappa3D_2,gamma3D_1,gamma3D_2

	integer nmol_data
	parameter(nmol_data=112)
	real*8 Mmol(nmol_data)
	character*10 molname(nmol_data)
	parameter(molname = (/'H2O   ','CO2   ','O3    ','N2O   ','CO    ','CH4   ',
     &	'O2    ','NO    ','SO2   ','NO2   ','NH3   ','HNO3  ','OH    ','HF    ',
     &	'HCl   ','HBr   ','HI    ','ClO   ','OCS   ','H2CO  ','HOCl  ','N2    ',
     &	'HCN   ','CH3Cl ','H2O2  ','C2H2  ','C2H6  ','PH3   ','COF2  ','SF6   ',
     &	'H2S   ','HCOOH ','HO2   ','O     ','ClONO2','NO+   ','HOBr  ','C2H4  ',
     &	'CH3OH ','CH3Br ','CH3CN ','CF4   ','C4H2  ','HC3N  ','H2    ','CS    ',
     &	'SO3   ','He    ','AlCl  ','AlF   ','AlH   ','AlO   ','SiS   ','SiO   ',
     &  'C2    ','Na    ','K     ','TiO   ','VO    ','FeH   ','C     ','CH2OH ',
     &  'CH3   ','H     ','N     ','NH    ','NH2   ','N2H3  ','AsH3  ','BeH   ',
     &  'CaF   ','CaH   ','CaO   ','CH    ','CH3F  ','CN    ','CP    ','CrH   ',
     &  'H2+   ','H3+   ','HeH+  ','KCl   ','KF    ','LiCl  ','LiF   ','LiH   ',
     &  'LiH+  ','MgF   ','MgH   ','MgO   ','NaCl  ','NaF   ','NH    ','NS    ',
     &  'OH+   ','P2H2_c','P2H2_t','PN    ','PO    ','PS    ','ScH   ','SH    ',
     &  'SiH   ','SiH4  ','TiH   ','H2Oem ','Na+   ','K+    ','H-    ','CH2   ',
     &  'H3O+  ','SiO2  '/))
	parameter(Mmol = (/     17.8851,  43.6918,  47.6511,  43.6947,  27.8081,  15.9272,  
     &	31.7674,  29.7889,  63.5840,  45.6607,  16.9072,  62.5442,  16.8841,  19.8619,  
     &	36.1973,  80.3271,   1.0070,  51.0760,  59.6379,  29.8088,  52.0765,  27.8112,  
     &	26.8304,  50.1116,  33.7598,  25.8499,  29.8518,  33.7516,  65.5261, 144.9081,  
     &	33.8332,  45.6731,  32.7593,  15.8794,  96.7366,  29.7813,  96.2063,  27.8507,  
     &	31.7949,  94.2413,  40.7302,  87.3580,  49.6543,  50.6424,   2.0014,  43.7539,  
     &	79.3792,   4.0030,  51.4530,  45.9799,  27.9895,  42.9800,  60.1500,  44.0845,
     &  24.0214,  22.9900,  39.0980,  63.8660,  66.9410,   56.853,  12.0107,  31.0339,
     &	15.0345,   1.0079,  14.0067,  15.0146,  16.0225,   31.037,  77.9500,  10.0200,
     &  59.0764,  41.0859,  56.0774,  13.0179,  34.0300,  26.0174,  42.9800,  53.0040,
     &   2.0151,   3.0238,   5.0105,  74.5513,  58.0967,  42.3940,  25.9390,   7.9500,
     &   7.9500,  43.3030,  25.3130,  40.3040,  58.4400,  41.9882,  15.0146,  46.0717,
     &  17.0070,  63.9600,  63.9600,  44.9800,  46.9700,  31.9800,  45.9600,  33.0729,
     &  29.0934,  32.1173,  48.8700,  17.8851,  22.9900,  39.0980,   1.0079,  14.0266,
     &  19.0232,  60.0800/))
	integer Catoms(nmol_data),Oatoms(nmol_data),Hatoms(nmol_data)
	parameter(Catoms = (/0,1,0,0,1,1,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,0,0,1,1,0,0,
     &				 1,1,0,2,2,0,1,0,
     &				 0,1,0,0,0,0,0,2,
     &				 1,1,2,1,4,3,0,1,
     &				 0,0,0,0,0,0,0,0,
     &				 2,0,0,0,0,0,1,1,
     &				 1,0,0,0,0,0,0,0,
     &                           0,0,0,1,1,1,1,0,
     &                           0,0,0,0,0,0,0,0,
     &                           0,0,0,0,0,0,0,0,
     &                           0,0,0,0,0,0,0,0,
     &                           0,0,0,0,0,0,0,1,
     &				 0,0 /))
	parameter(Oatoms = (/1,2,3,1,1,0,
     &				 2,1,2,2,0,3,1,0,
     &				 0,0,0,1,1,1,1,0,
     &				 0,0,2,0,0,0,1,0,
     &				 0,2,2,1,3,1,1,0,
     &				 1,0,0,0,0,0,0,0,
     &				 3,0,0,0,0,1,0,1,
     &				 0,0,0,1,1,0,0,1,
     &				 0,0,0,0,0,0,0,0,
     &                           0,0,1,0,0,0,0,0,
     &                           0,0,0,0,0,0,0,0,
     &                           0,0,0,1,0,0,0,0,
     &                           1,0,0,0,1,0,0,0,
     &                           0,0,0,1,0,0,0,0,
     &				 1,2  /))
	parameter(Hatoms = (/2,0,0,0,0,4,
     &				 0,0,0,0,3,1,1,1,
     &				 1,1,1,0,0,2,1,0,
     &				 1,3,2,2,6,3,0,0,
     &				 2,2,0,0,0,0,1,4,
     &				 4,3,3,0,2,1,2,0,
     &                           0,0,0,0,1,0,0,0,
     &				 0,0,0,0,0,1,0,3,
     &				 3,1,0,1,2,3,3,1,
     &                           0,1,0,1,3,0,0,1,
     &                           2,3,1,0,0,0,0,1,
     &                           1,0,1,0,0,0,1,0,
     &                           1,2,2,0,0,0,1,1,
     &                           1,4,1,2,0,0,1,2,
     &				 3,0 /))
	real*8,allocatable :: a_therm(:),a_press(:)
	integer n_voigt,n_instr
	logical HITEMP,opacitymode,compute_opac,Mp_from_logg,trend_compute
	integer nPom,nTom
	character*500 opacitydir,specresfile,starfile
	character*500,allocatable :: instrument(:)
	real*8,allocatable :: instr_ntrans(:)
	real*8 Tmin,Tmax,minTprofile,maxTprofile,fcloud_default,chimax,r_nuc
	real*8 sintheta(360),costheta(360)
	logical,allocatable :: do_dB(:)
	real*8 COret,COerr(2)
	
	logical sinkZ
	real*8 alphaZ

	real*8,allocatable :: flux(:,:),obsA(:,:),phase(:,:,:),obsLightCurve(:,:)
	real*8,allocatable :: timeLightCurve(:),theta_phase(:)
	integer ncc,nphase,n2d,i2d,nLightCurve
	logical cloudcompute,useDRIFT,singlecloud,makeimage
	logical,allocatable :: docloud(:,:)
	real*8,allocatable :: cloudfrac(:),XCloud(:,:),XeqCloud(:,:),XeqCloud_old(:,:)
	real*8,allocatable :: nabla_ad(:),grav(:)

	integer,allocatable,dimension(:) :: L_imol,L_iiso,L_nclose,L_ilam
	real*8,allocatable,dimension(:) :: L_Aul,L_freq,L_Elow,L_lam,L_S0,L_S
	real*8,allocatable,dimension(:) :: L_gamma_air,L_gamma_self
	real*8,allocatable,dimension(:) :: L_a_therm,L_a_press,L_n
	real*8,allocatable,dimension(:) :: L_gu,L_gl,L_Saver
	logical,allocatable,dimension(:) :: L_do
	
	integer,allocatable :: ig_comp(:)
	integer ng_comp

	type CIA_pair
		character*20 name
		character*500 filename
		real*8,allocatable :: T(:),Cabs(:,:)
		integer nT,imol1,imol2
	end type CIA_pair
	
	type(CIA_pair),allocatable :: CIA(:)
	real*8 cia_mixrat(nmol_data)

	type CloudType
		real*8 P,dP,s,column
		real*8 coverage,frain
		real*8,allocatable :: rv(:),w(:),M(:)					! dimension nsize
		real*8,allocatable :: frac(:,:),sigma(:)
		real*8 rho,amin,amax,fmax,porosity,fcarbon,reff,veff
		logical blend,haze,condensates
		real*8 fcond,mixrat,tau,lam
		real*8,allocatable :: Kabs(:,:),Ksca(:,:),Kext(:,:)			! dimension nsize,nlam
		character*500 file,Kzzfile
		character*20 standard,ptype,type
		character*500 species,hazetype
		real*8 fHazeSiO,fHazeTiO2,fHazeTholin,fHazeAl2O3,fHazeFe
		integer nr,nsubgrains
		real*8 tmix,betamix,Kzz,Kscale,Sigmadot
c simple cloud from Jo
		real*8 ff,g1,g2,kappa,albedo,kappa_haze,albedo_haze
		logical simplecloud
	end type CloudType

	type(CloudType),allocatable :: Cloud(:) 

cPoints for the temperature structure
	real*8,allocatable :: P_point(:),T_point(:)
	integer n_points

	type RetrievalPar
		character*500 keyword
		real*8 xmin,xmax,x0,dx,value,error1,error2
		logical logscale,squarescale,opacitycomp,increase
		integer n
	end type RetrievalPar

	type(RetrievalPar),allocatable :: RetPar(:)
	integer n_ret

	type Parameter3D
		character*500 keyword
		real*8 xmin,xmax,x,pow
		logical logscale
	end type Parameter3D

	type(Parameter3D),allocatable :: Par3D(:)
	integer n_Par3D
	
	type ObservedSpec
		character*500 file
		character*10 type
		real*8,allocatable :: lam(:),y(:),dy(:),R(:),Rexp(:),model(:)
		integer,allocatable :: ilam(:)
		real*8 beta,scale
		integer ndata,i2d,iphase
		logical spec,scaling
		integer nlam,nt
		real*8,allocatable :: LC(:,:),dLC(:,:),t(:),dt(:)
	end type ObservedSpec
	type(ObservedSpec),allocatable :: ObsSpec(:)

	integer npop,ngen,nobs
	logical gene_cross
	
	end module GlobalSetup
	

	module modComputeT
	IMPLICIT NONE
	real*8,allocatable :: Si_prev(:,:,:),Cp_prev(:)
	integer iaverage
	end module modComputeT
	


	module CloudModule
	IMPLICIT NONE
	integer nr_cloud
	real*8,allocatable :: CloudP(:),CloudT(:),CloudR(:),Clouddens(:),CSnmol(:)
	real*8,allocatable :: ATP(:),BTP(:),rhodust(:),atoms_cloud(:,:),maxT(:),mu(:),xv_bot(:)
	character*25,allocatable :: CSname(:)
	integer nCS,nnr
	logical,allocatable :: ice(:)

	end module



	module RetrievalMod
	implicit none
	integer imodel
	real*8 bestlike,chi2_0
	real*8,allocatable :: dvarq(:),bestvar(:)
	real*8,allocatable :: obsA0(:),obsA1(:),obsA2(:),dobsA(:,:)
	real*8,allocatable :: emis0(:),emis1(:),emis2(:),demis(:,:)
	real*8,allocatable :: emisR0(:),emisR1(:),emisR2(:),demisR(:,:)
	end module RetrievalMod




	module RandomWalkModule
	IMPLICIT NONE
	integer NY
	parameter(NY=1000)
	real*8 phi(NY),yy(NY)
	end module RandomWalkModule

	
	module Struct3D
	implicit none
	integer nlong,nlatt
	integer n3D,nnu0
c	parameter(n3D=10,nnu0=10)
c	parameter(nlong=36,nlatt=18)
c	real*8 long(nlong),latt(nlatt)	!(Lambda, Phi)
c	real*8 tanx(nlong),tany(nlong)
c	real*8 cost2(nlatt),beta3D_eq(nlong),ibeta3D_eq(nlong)
c	integer ibeta(nlong,nlatt),inu3D(nlong,nlatt)

	real*8,allocatable :: long(:),latt(:)	!(Lambda, Phi)
	real*8,allocatable :: tanx(:),tany(:)
	real*8,allocatable :: cost2(:),beta3D_eq(:),ibeta3D_eq(:)
	integer,allocatable :: ibeta(:,:),inu3D(:,:)

	real*8,allocatable :: R3D(:,:),R3D2(:,:)
	end module Struct3D



