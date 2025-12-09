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
	parameter(Rjup=6.9911d9)
c	parameter(Rjup=7.1492d9)
	parameter(year=24d0*60d0*60d0*365.25d0)
	parameter(micron=1d-4)
	parameter(Rgas=8.3144621d7)
	parameter(Avogadro=6.022136736d23)
	parameter(atm=1.01325d0)	! bar
	
	end module Constants


	module dbl2str_mod
	implicit none
	private
	public :: dbl2string

	interface dbl2string
		procedure dbl2string_real, dbl2string_dbl, dbl2string_qdr, dbl2string_int
	end interface dbl2string
	contains
	function dbl2string_real(x,form) result(res)
		real x
		character(len=20) res
		character,intent(in),optional :: form*(*)

		if(present(form)) then
			write(res,form) x
		else
			write(res,*) x
		endif
	end function dbl2string_real

	function dbl2string_dbl(x,form) result(res)
		real*8 x
		character(len=20) res
		character,intent(in),optional :: form*(*)

		if(present(form)) then
			write(res,form) x
		else
			write(res,*) x
		endif
	end function dbl2string_dbl

	function dbl2string_qdr(x,form) result(res)
		real*16 x
		character(len=20) res
		character,intent(in),optional :: form*(*)

		if(present(form)) then
			write(res,form) x
		else
			write(res,*) x
		endif
	end function dbl2string_qdr

	function dbl2string_int(x,form) result(res)
		integer x
		character(len=20) res
		character,intent(in),optional :: form*(*)
		res=' '
		return
	end function dbl2string_int

	end module dbl2str_mod


	module int2str_mod
	implicit none
	private
	public :: int2string

	interface int2string
		procedure int2string_real, int2string_dbl, int2string_qdr, int2string_int
	end interface int2string
	contains
	function int2string_real(x,form) result(res)
		real x
		character(len=20) res
		character,intent(in),optional :: form*(*)
		res=' '
		return
	end function int2string_real

	function int2string_dbl(x,form) result(res)
		real*8 x
		character(len=20) res
		character,intent(in),optional :: form*(*)
		res=' '
		return
	end function int2string_dbl

	function int2string_qdr(x,form) result(res)
		real*16 x
		character(len=20) res
		character,intent(in),optional :: form*(*)
		res=' '
		return
	end function int2string_qdr

	function int2string_int(x,form) result(res)
		integer x
		character(len=20) res
		character,intent(in),optional :: form*(*)

		if(present(form)) then
			write(res,form) x
		else
			write(res,*) x
		endif
	end function int2string_int

	end module int2str_mod


c========================================================
c Interfaces for input/output subroutines
c========================================================
	module OutputModule
	use dbl2str_mod
	use int2str_mod
	IMPLICIT NONE

	interface
		subroutine outputform(string,form)
			character :: string*(*)
			character,intent(in),optional :: form*(*)
		end subroutine outputform
	end interface


	end module OutputModule
	


c=========================================================================================
c global setup for the code
c=========================================================================================
	module GlobalSetup
	use OutputModule
	IMPLICIT NONE
	real*8 Mplanet,Rplanet,Pplanet,loggPlanet				! mass and radius of the planet at pressure Pplanet
	real*8 Tstar,Rstar,Mstar,logg,Dplanet,Dplanet0,Lplanet
	real*8 orbit_P,orbit_e,orbit_omega,orbit_inc
	real*8,allocatable :: dens(:),T(:),P(:),Ndens(:),Tin(:)	! radius
	real*8,allocatable :: dust_dens(:,:)					! radius, component
	real*8,allocatable :: R(:),Hp(:)						! radius
	real*8,allocatable :: mixrat(:)							! component
	real*8,allocatable :: mixrat_r(:,:),mixrat_old_r(:,:)	! radius,component
	real*8,allocatable :: Cabs(:,:,:,:),Csca(:,:)			! radius,wav,g,velocity
	real*8,allocatable :: cloud_dens(:,:)					! radius, cloud
	real*8,allocatable :: Fstar(:)							! wavelength
	real*8,allocatable :: tau1depth(:,:),cloudtau(:,:)		! ncc,wavelength
	real*8,allocatable :: Cabs_mol(:,:,:,:),Cext_cont(:,:)
	real*8,allocatable :: Pswitch_mol(:),abun_switch_mol(:)
	real*8,allocatable :: x_el(:)
	integer nangle_Jscat,nvel
	parameter(nangle_Jscat=60)
	real*8,allocatable :: Jscat(:,:)						! radius, angle
	integer nT,np,nr,nmol,nlam		! #T, #P, #radial points, #molecules, #wavelength bins, #obs
	integer nlines,ng,ncia,nclouds,nTiter,i3D,i_alb,nest_update
	character*1000 outputdir
	character*1000,allocatable :: commandargs(:),line_add_ret(:)
	integer ncommandargs,n_add_ret
	integer idum,maxiter,miniter,Nphot0,idum0,iWolk
!$OMP THREADPRIVATE(idum)
	logical retrieval,outputopacity,do_cia,gridTPfile,scattering,scattstar,anisoscattstar,computeT,computecontrib,do_rayleigh,isoFstar,writefiles
	logical dochemistry,free_tprofile,condensates,faircoverage,speclimits,mapCOratio,randomseed,useXS,modelfail,projectedD
	logical,allocatable :: includemol(:),diseqmol(:),didcondens(:),lamemis(:),lamtrans(:),opacitymol(:)
	logical,allocatable :: includemol_raytrace(:),includemol_default(:)
	real*8 lam1,lam2,specres,Pmin,Pmax,epsCk,distance,TP0,dTP,TeffP,twind,epsiter,specres_LR
	real*8 gammaT1,gammaT2,kappaT,betaT,alphaT,metallicity0,vfrag,betaF,kappaUV0,scaleUV,gammaUV,eps_dup
	logical mixratfile,par_tprofile,adiabatic_tprofile,domakeai,modelsucces,useobsgrid,blackbodystar,grey_isoT
	logical didcondens_chem,resume_multinest,disequilibrium,const_eff_multinest,outflow,useDLMie,importance_nested_sampling
	character*500 TPfile,particledir,retrievaltype,planetparameterfile,planetname,element_abun_file,pargridfile,deepredisttype
	real*8 metallicity,COratio,PQ,mixP,PRplanet,maxchemtime,TiScale,f_multinest,tol_multinest,dmetallicity
	real*8 Kzz,Kzz_offset,Kzz_max,SiOratio,NOratio,fDay,betapow,Kxx,Kyy,vxx,powvxx,night2day,pole2eq,Rp_range,tauLW
	real*8 Kzz_deep,Kzz_1bar,Kzz_P,Kzz_contrast,SOratio,Tsurface,hotspotshift0,exp_ad,Tsurface0
	logical gamma_equal,dopostequalweights,inverseCOratio,setsurfpressure,fixnight2day,tidallock,distrUV
	logical transspec,emisspec,dosimplerainout,computeLC,doscaleR,complexKzz,convectKzz,SCKzz,writeWolk,dotranshide,ComputeTeff
	real*8 cutoff_abs,cutoff_lor,eps_lines,maxtau,factRW,Tform,Pform,f_dry,f_wet,scale_fe
	real*8,allocatable :: lam(:),freq(:),dfreq(:),dlam(:),blam(:,:),surface_emis(:),surface_props(:,:)
	real*8 f_ice,f_grass,f_snow,f_water,f_surface(5)
	real*8 bdrf_args(4,5)
	integer bdrf_type(5),n_surface,nRTatm
	character*500 bdrf_file(5)
	real*8,allocatable :: gg(:),wgg(:),obsA_contr(:,:),flux_contr(:,:),obsA_LC(:,:)
	real*8,allocatable :: ZZ(:,:,:),TZ(:)	! partition function
	real*8 planetform_fdust,planetform_fplan,planetform_Mstart,planetform_SolidC,vrot0,vrot_max
	real*8,allocatable :: velocity(:),Kzz_convect(:),Kzz_g(:),Kzz_b(:)
	real*8 planetform_Macc,planetform_Dmigrate,planetform_Rend,TeffPoutput,Hydrogenloss,taurexsmooth
	logical planetform,massprior,retrievestar,simAb_converge,log_emis,randomstart,logTprofile,taurexprofile
c for exchange when computing secondary atmosphere
	real*8 Toutgas,Poutgas
	real*8 molfracs_atoms_outgas(41),MMW0,fH2O	! fH2O is the total mass fraction in Water compared to the mass of the planet
	logical secondary_atmosphere,constant_g,forceEbalance,fixMMW,WaterWorld,useomp,fullcovmat,cloudoverlap
	logical,allocatable :: RTgridpoint(:),computelam(:)
	character*20 aim3D
	character*2 standardstarname
	
	integer nrsurf,nrstepchem
	real*8 Psurf
	
	real*8,allocatable :: model_err_rel(:),model_err_abs(:),model_err_lam(:)
	integer nmodel_err
	
	real*8 nC_PAH,mixrat_PAH,rad_optEC,Eg_optEC,mixrat_optEC,mixrat_optEC0
	real*8,allocatable :: mixrat_optEC_r(:)

	real*8 Mp_prior,dMp_prior,surfacealbedo,MSimAb
	character*20 surfacetype
	character*500 surfacefile,Full3Ddir
	integer nTZ,nspike,nai,nboot,npew,nscaleR
	integer,allocatable :: niso(:),instr_nobs(:)
	real*8,allocatable :: MMW(:),tauUV(:),kappaUV(:),Tprev3D(:)
	real*8,allocatable :: PTaverage3D(:,:),mixrat_average3D(:,:,:)
	logical fulloutput3D,deepredist,readFull3D
	integer nBB
	parameter(nBB=10000)
	character*500 formationcommand

	real*8,allocatable :: tau_Vpoint(:),tau_IRpoint(:),dT_Vpoint(:),dT_IRpoint(:)
	real*8,allocatable :: Ppoint(:),dTpoint(:)
	logical pos_dT_lowest,pos_dT,freePT_fitT,freePT_fitP,useEOS
	real*8 PrefTpoint,wiggle_err
	integer nVpoints,nIRpoints,nTpoints

	real*8 Rring,dRring,tauRing
	real*8,allocatable :: FRing(:)
	logical doRing,doRingCloud

	logical do3D,init3D
	real*8 par3Dsteepness
	
	real*8 global_like,global_chi2
	
	real*8,allocatable :: Cov_L_loc(:),Cov_a_loc(:),Cov_lam_loc(:)
	integer Cov_n_loc

	integer n_ff
	parameter(n_ff=1)
	integer molnr_ff(n_ff)
	parameter(molnr_ff = (/ 45 /))
	logical use_ff
	integer nmol_data
	parameter(nmol_data=126)
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
     &  'H3O+  ','SiO2  ','S2    ','SO    ','15NH3 ','SiH2  ','NaH   ','Mg    ',
     &  'Fe    ','Al    ','H2SO4 ','Ca    ','Ti    ','Si    ','S     ','P     '/))
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
     &  19.0232,  60.0800,  64.1300,  48.0644,  16.9072,  30.1014,  23.9977,  24.3050,
     &  55.8450,  26.9815,  98.0790,  40.0784,  47.8671,  27.9769,  32.0660,  30.9738/))
	integer Catoms(nmol_data),Oatoms(nmol_data),Natoms(nmol_data),Hatoms(nmol_data),tot_atoms(nmol_data)
	parameter(Natoms = (/0,0,0,1,0,0,
     &				 0,1,0,1,1,1,0,0,
     &				 0,0,0,0,0,0,0,2,
     &				 1,0,0,0,0,0,0,0,
     &				 0,0,0,0,1,1,0,0,
     &				 0,0,1,0,0,1,0,0,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,1,1,1,2,0,0,
     &				 0,0,0,0,0,1,0,0,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,0,0,0,0,1,1,
     &				 0,0,0,1,0,0,0,0,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,0,0,1,0,0,0,
     &				 0,0,0,0,0,0,0,0  /))
	parameter(Catoms = (/0,1,0,0,1,1,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,0,0,1,1,0,0,
     &				 1,1,0,2,2,0,1,0,
     &				 0,1,0,0,0,0,0,2,
     &				 1,1,2,1,4,3,0,1,
     &				 0,0,0,0,0,0,0,0,
     &				 2,0,0,0,0,0,1,1,
     &				 1,0,0,0,0,0,0,0,
     &				 0,0,0,1,1,1,1,0,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,0,0,0,0,0,1,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,0,0,0,0,0,0 /))
	parameter(Oatoms = (/1,2,3,1,1,0,
     &				 2,1,2,2,0,3,1,0,
     &				 0,0,0,1,1,1,1,0,
     &				 0,0,2,0,0,0,1,0,
     &				 0,2,2,1,3,1,1,0,
     &				 1,0,0,0,0,0,0,0,
     &				 3,0,0,0,0,1,0,1,
     &				 0,0,0,1,1,0,0,1,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,1,0,0,0,0,0,
     &				 0,0,0,0,0,0,0,0,
     &				 0,0,0,1,0,0,0,0,
     &				 1,0,0,0,1,0,0,0,
     &				 0,0,0,1,0,0,0,0,
     &				 1,2,0,1,0,0,0,0,
     &				 0,0,4,0,0,0,0,0  /))
	parameter(Hatoms = (/2,0,0,0,0,4,
     &				 0,0,0,0,3,1,1,1,
     &				 1,1,1,0,0,2,1,0,
     &				 1,3,2,2,6,3,0,0,
     &				 2,2,0,0,0,0,1,4,
     &				 4,3,3,0,2,1,2,0,
     &				 0,0,0,0,1,0,0,0,
     &				 0,0,0,0,0,1,0,3,
     &				 3,1,0,1,2,3,3,1,
     &				 0,1,0,1,3,0,0,1,
     &				 2,3,1,0,0,0,0,1,
     &				 1,0,1,0,0,0,1,0,
     &				 1,2,2,0,0,0,1,1,
     &				 1,4,1,2,0,0,1,2,
     &				 3,0,0,0,3,2,1,0,
     &				 0,0,2,0,0,0,0,0 /))
	parameter(tot_atoms = (/3,3,3,3,2,5,
     &				 2,2,3,3,4,5,2,2,
     &				 2,2,1,2,3,4,3,2,
     &				 3,5,4,4,8,4,4,7,
     &				 3,5,3,1,5,2,3,6,
     &				 6,5,6,5,6,5,2,2,
     &				 4,1,2,2,2,2,2,2,
     &				 2,1,1,2,2,2,1,5,
     &				 4,1,1,2,3,5,4,2,
     &				 2,2,2,2,5,2,2,2,
     &				 2,3,2,2,2,2,2,2,
     &				 2,2,2,2,2,2,2,2,
     &				 2,4,4,2,2,2,2,2,
     &				 2,5,2,3,1,1,1,3,
     &				 4,3,2,2,4,3,2,1,
     &				 1,1,7,1,1,1,1,1 /))
	real*8,allocatable :: a_therm(:),a_press(:)
	integer n_voigt,n_instr
	logical opacitymode,Mp_from_logg,trend_compute,backgroundgas(nmol_data),dobackgroundgas
	integer nPom,nTom
	character*500 opacitydir,specresfile,starfile
	character*500,allocatable :: instrument(:)
	real*8,allocatable :: instr_ntrans(:)
	real*8 Tmin,Tmax,minTprofile,maxTprofile,chimax,PphotMol(nmol_data)
	real*8 sintheta(360),costheta(360)
	logical,allocatable :: do_dB(:)
	real*8 COret,COerr(2),WWInit_Pmax,WWInit_Pmin,WWInit_Pplanet,WWInit_mixrat(nmol_data)
	
	character*10 fixmol_name(nmol_data)
	real*8 fixmol_abun(nmol_data),fixmol_P(nmol_data)
	integer nfixmol,ifixmol(nmol_data)

	integer isotope(nmol_data)
	real*8 f_isotope(nmol_data)

	logical sinkZ
	real*8 alphaZ

	real*8,allocatable :: flux(:,:),obsA(:,:),phase(:,:,:),obsLightCurve(:,:)
	real*8,allocatable :: timeLightCurve(:),theta_phase(:),obsA_split(:,:)
	integer ncc,nphase,n2d,i2d,nLightCurve
	logical makeimage,makemovie
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
		character*20 opacitytype,type
		real*8 P,dP,xi,Pmax,Pmin,Ptau,Phi,coverage
		real*8,allocatable :: rv(:),M(:)					! dimension nsize
		real*8,allocatable :: frac(:,:),sigma(:),cryst(:,:),abun(:),xv_bot(:),porosity(:)
		real*8 rho,fmax,porosity0,reff,veff,rpow,Pref,rnuc,fractalDim,rnuc_phot
		logical blend,haze,condensates,rainout,globalKzz,computecryst,coagulation
		logical onepart,freeflow_nuc,freeflow_con,condenseNaK
		real*8 mixrat,tau,lref,cryst0,e1_par,e2_par,Kref
		real*8,allocatable :: Kabs(:,:),Ksca(:,:),Kext(:,:),g(:,:)		! dimension nsize,nlam
		real*8,allocatable :: F11(:,:,:) ! dimension nsize,nlam,nangle
		character*500 species,hazetype,file,composition
		integer nmat,nlam,fixcloud
c			the parameter fixcloud fixes the cloud after 'fixcloud' PT iterations. 
c			Setting it to 1 means the cloud is not iterated with the PT structure.
c			Setting it to 0 means the cloud is always recomputed every iteration.
		character*500,allocatable :: lnkfile(:,:),material(:),condensate(:)
		real*8 Kzz,Sigmadot,xm_bot,Sigmadot_phot
		real*8 kappa,albedo,kpow,klam,g0
		real*8 kappa_Gauss,lam_Gauss,dlam_Gauss
		real*8,allocatable :: e1(:,:,:),e2(:,:,:),rho_mat(:),KeFile(:,:),KaFile(:,:),KsFile(:,:),gFile(:,:)
		integer,allocatable :: nax(:)
		logical usefsed,computeJn,EqChemBoundary,globalGasMixing
		real*8 fsed_alpha,fsed_beta,Srainout,fstick,x_slider
	end type CloudType

	type(CloudType),allocatable :: Cloud(:) 

	type photochem
		real*8,allocatable :: react(:),product(:),abun(:,:)
		integer nreact
		logical atomic
		real*8 f_eff,haze,scaleKappa
	end type photochem
	type(photochem),allocatable :: PhotoReacts(:)
	integer nPhotoReacts

	type RetrievalPar
		character*500 keyword
		real*8 xmin,xmax,x0,dx,value,error1,error2
		logical logscale,squarescale,opacitycomp,increase
	end type RetrievalPar

	type(RetrievalPar),allocatable :: RetPar(:)
	integer n_ret

	type Parameter3D
		character*500 keyword
		real*8 xmin,xmax,x
		logical logscale,multiply
	end type Parameter3D

	type(Parameter3D),allocatable :: Par3D(:)
	integer n_Par3D
	
	type ObservedSpec
		character*500 file,filter
		character*10 type
		real*8,allocatable :: lam(:),y(:),dy(:),R(:),Rexp(:),model(:)
		integer,allocatable :: ilam(:)
		real*8 beta,scale,slope,adderr,dscale,fscale,offset
		integer ndata,i2d,iphase
		logical spec,scaling
		integer nlam,nt
		real*8,allocatable :: LC(:,:),dLC(:,:),t(:),dt(:)
		real*8,allocatable :: f(:,:)
		real*8 Cov_a,Cov_L,Cov_offset
		real*8,allocatable :: Cov_L_loc(:),Cov_a_loc(:),Cov_lam_loc(:)
		integer Cov_n_loc
	end type ObservedSpec
	type(ObservedSpec),allocatable :: ObsSpec(:)

	integer npop,ngen,nobs,npost
	real*8 epsinit_MCMC
	logical gene_cross
	
	end module GlobalSetup
	

	module CloudModule
	IMPLICIT NONE
	integer nr_cloud
	real*8,allocatable :: CloudP(:),CloudT(:),CloudR(:),Clouddens(:),CSnmol(:),SatRat(:)
	real*8,allocatable :: xv(:,:),xn(:),xc(:,:),xm(:),rpart(:),xa(:),xs(:)
	real*8,allocatable :: ATP(:),BTP(:),rhodust(:),atoms_cloud(:,:),maxT(:),mu(:),xv_bot(:)
	character*25,allocatable :: CSname(:)
	integer nCS,nnr
	logical,allocatable :: ice(:)
	real*8,allocatable :: Tevap(:,:),Tdist(:,:),molfracs_atoms_cloud(:,:)

	end module

	module TimingModule
	IMPLICIT NONE
	real*8 timechem,timecloud,timetemp
	integer itimechem,itimecloud,itimetemp
	integer ctimechem,ctimecloud,ctimetemp,rate
	end module

	module RetrievalMod
	implicit none
	integer imodel
	real*8 bestlike,chi2_0,bestchi2
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
	integer nlong,nlatt,nalbedo_iter
	integer n3D,nnu0

	real*8,allocatable :: long(:),latt(:)	!(Lambda, Phi)
	real*8,allocatable :: tanx(:),tany(:)
	real*8,allocatable :: cost2(:),beta3D_eq(:),x3D_eq(:)
	real*8,allocatable :: beta3D(:),x3D(:)
	integer,allocatable :: ibeta(:,:)

	real*8,allocatable :: R3D(:,:),R3D2(:,:),local_albedo(:)
	real*8 hotspotshift
	end module Struct3D


c===================================================================================
c Module containing the setup for the formation subroutine
c===================================================================================
	module FormationModule
	IMPLICIT NONE
c Disk setup
	integer Nr_disk
c Static settings for the disk temperature and atomic abundances
	real*8,allocatable :: R_disk(:),T_disk(:)
c The abun_???(:,:) arrays contain the elemental abundances in each of the disk components: 
c			gas, dust, planetessimals
c Units are number of atoms per unit gram. This way we can easily accrete a given mass and 
c compute what this means for the total number of atoms for every element
	real*8,allocatable :: abun_gas(:,:),abun_dust(:,:),abun_planet(:,:)
c Disk density setup where the dust and planet densities are the max values possible
	real*8,allocatable :: d2g_disk(:),p2g_disk(:)
	real*8 Mstar,Lstar,mu,kappa_r,alpha_disk,M_acc,d2g_T
	
	end module


c===================================================================================
c Module containing the abundances of the atoms and their names
c===================================================================================
	module AtomsModule
	IMPLICIT NONE

	INTEGER,parameter            :: N_atoms = 41
	CHARACTER*40                 :: names_atoms(N_atoms)
	DOUBLE PRECISION             :: molfracs_atoms(N_atoms)
	DOUBLE PRECISION             :: mass_atoms(N_atoms)

	end module AtomsModule
c===================================================================================
c===================================================================================
