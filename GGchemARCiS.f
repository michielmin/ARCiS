**********************************************************************
      MODULE ARCiS_GGCHEM
**********************************************************************
      implicit none
      integer,allocatable :: linkmol(:),linkele(:)
      logical GGCHEM_P_iter
      end MODULE ARCiS_GGCHEM


***********************************************************************
	subroutine init_GGchem(mol_names_in,n_mol_in,condensates)
***********************************************************************
      use PARAMETERS,ONLY: elements,abund_pick,model_dim,model_pconst,
     >     model_struc,model_eqcond,Npoints,useDatabase,verbose,
     >     Tfast,Tmin,Tmax,pmin,pmax,nHmin,nHmax,pick_mfrac,
     >     abund_file,struc_file,remove_condensates,phyllosilicates,
     >     initchem_info,auto_atmos,stop_after_init,Mpl,Rpl,gamma,
     >     adapt_cond,method_eqcond,Nseq,Tseq,Tmin_atmos
      use EXCHANGE,ONLY: chemcall,chemiter,ieqcond,ieqconditer,
     >                   itransform,
     >                   H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,
     >                   Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,
     >                   As,Se,Br,Kr,Rb,Sr,Y,Zr,W
      use DATABASE,ONLY: NLAST

      use CHEMISTRY,ONLY: NewBackIt,NewFullIt,NewBackFac,NewPreMethod,
     >                    NewFastLevel,dispol_file,NMOLE,NELM,cmol,Natmax
      use DUST_DATA,ONLY: DustChem_file,NELEM,eps=>eps0,mass,muH,elnam,amu,bar,MEarth,REarth,dust_nam
      use ARCiS_GGCHEM
      implicit none
      integer,parameter :: qp=selected_real_kind(33,4931)
	integer n_mol_in,ii
	character*10 :: mol_names_in(n_mol_in),uppername
      integer :: i,j,nr
      real(kind=qp) :: m,val,abund(74,4),eps0(NELEM),epsH,mfrac(NELEM)
       real(kind=qp) :: addH2O
	  character(len=2) :: el
      character(len=20) :: elname
      character(len=10) :: source(4)
      character(len=200) :: line
      logical :: found,condensates
	character*100 homedir
	call getenv('HOME',homedir) 

      !-------------------------
      ! ***  default values  ***
      !-------------------------
      dispol_file(1) = trim(homedir) // '/ARCiS/Data/GGchem/dispol_BarklemCollet.dat'
      dispol_file(2) = trim(homedir) // '/ARCiS/Data/GGchem/dispol_StockKitzmann_withoutTsuji.dat'
      dispol_file(3) = trim(homedir) // '/ARCiS/Data/GGchem/dispol_WoitkeRefit.dat'
      dispol_file(4) = trim(homedir) // '/ARCiS/src/dispol_Burcat.dat'
      DustChem_file  = trim(homedir) // '/ARCiS/Data/GGchem/DustChem.dat'

      elements     = 'H He C N O Na Mg Si Fe Al Ca Ti S Cl K Li P V el'

      auto_atmos         = .false.
      adapt_cond         = .false.
      stop_after_init    = .false.
      abund_pick   = 3
      model_dim    = 1
      model_struc  = 0
      verbose      = 0
      Mpl          = MEarth
      Rpl          = REarth
      gamma        = 7.0/5.0
      Tmin_atmos   = 0.0
      method_eqcond= 2
      Natmax       = 16
      Nseq         = 1
      pick_mfrac         = .false.
      abund_pick   = 3
      model_eqcond = condensates
      remove_condensates = condensates
      phyllosilicates    = condensates
      model_dim    = 0
      model_pconst = GGCHEM_P_iter
      model_struc  = 0
      initchem_info      = .false.
      Npoints      = 100
      Tfast        = 1000.d0
      Tmin         = 100.d0
      Tmax         = 6000.d0
      pmin         = 1.d0*bar
      pmax         = 1.d0*bar
      nHmin        = 4.d+19
      nHmax        = 4.d+19
      UseDataBase  = .true.
      NewFullIt    = .true.
      NewBackIt    = 10
      NewBackFac   = 1.E+10
      NewFastLevel = 1
      NewPreMethod = 2


      write(*,*) 
      write(*,*) "elemental abundances and masses ..."
      write(*,*) "==================================="
      elnam(1)  = 'H '
      elnam(2)  = 'He'
      elnam(3)  = 'Li'
      elnam(4)  = 'Be'
      elnam(5)  = 'B '
      elnam(6)  = 'C '
      elnam(7)  = 'N '
      elnam(8)  = 'O '
      elnam(9)  = 'F '
      elnam(10) = 'Ne'
      elnam(11) = 'Na'
      elnam(12) = 'Mg'
      elnam(13) = 'Al'
      elnam(14) = 'Si'
      elnam(15) = 'P '
      elnam(16) = 'S '
      elnam(17) = 'Cl'
      elnam(18) = 'Ar'
      elnam(19) = 'K '
      elnam(20) = 'Ca'
      elnam(21) = 'Sc'
      elnam(22) = 'Ti'
      elnam(23) = 'V '
      elnam(24) = 'Cr'
      elnam(25) = 'Mn'
      elnam(26) = 'Fe'
      elnam(27) = 'Co'
      elnam(28) = 'Ni'
      elnam(29) = 'Cu'
      elnam(30) = 'Zn'
      elnam(31) = 'Ga'
      elnam(32) = 'Ge'
      elnam(33) = 'As'
      elnam(34) = 'Se'
      elnam(35) = 'Br'
      elnam(36) = 'Kr'
      elnam(37) = 'Rb'
      elnam(38) = 'Sr'
      elnam(39) = 'Y '
      elnam(40) = 'Zr'
      elnam(41) = 'W '

*     --------------------
*     ***  Atommassen  ***
*     --------------------
      mass(H)  = 1.008  * amu  
      mass(He) = 4.0026 * amu
      mass(Li) = 6.94   * amu  
      mass(Be) = 9.0122 * amu
      mass(B)  = 10.81  * amu  
      mass(C)  = 12.011 * amu  
      mass(N)  = 14.007 * amu  
      mass(O)  = 15.999 * amu  
      mass(F)  = 18.998 * amu
      mass(Ne) = 20.180 * amu 
      mass(Na) = 22.990 * amu
      mass(Mg) = 24.305 * amu 
      mass(Al) = 26.982 * amu
      mass(Si) = 28.085 * amu  
      mass(P)  = 30.974 * amu
      mass(S)  = 32.06  * amu  
      mass(Cl) = 35.45  * amu  
      mass(Ar) = 39.948 * amu  
      mass(K)  = 39.098 * amu 
      mass(Ca) = 40.078 * amu  
      mass(Sc) = 44.956 * amu
      mass(Ti) = 47.867 * amu  
      mass(V)  = 50.942 * amu 
      mass(Cr) = 51.996 * amu 
      mass(Mn) = 54.938 * amu
      mass(Fe) = 55.845 * amu  
      mass(Co) = 58.933 * amu
      mass(Ni) = 58.693 * amu 
      mass(Cu) = 63.546 * amu  
      mass(Zn) = 65.38  * amu  
      mass(Ga) = 69.723 * amu  
      mass(Ge) = 72.63  * amu  
      mass(As) = 74.922 * amu
      mass(Se) = 78.96  * amu  
      mass(Br) = 79.904 * amu  
      mass(Kr) = 83.798 * amu  
      mass(Rb) = 85.468 * amu 
      mass(Sr) = 87.62  * amu  
      mass(Y ) = 88.906 * amu
      mass(Zr) = 91.224 * amu  
      mass(W ) = 183.84 * amu       

*     ---------------------------------------
*     ***      element abundancies        ***
*     ---------------------------------------
*     Grevesse + Noels (1996, "photosphere"):
*     ---------------------------------------
      eps(H)  = 12.00 D0
      eps(He) = 10.99 D0
      eps(Li) =  1.16 D0
      eps(C)  =  8.55 D0
      eps(N)  =  7.97 D0
      eps(O)  =  8.87 D0
      eps(Ne) =  8.08 D0
      eps(Na) =  6.33 D0
      eps(Mg) =  7.58 D0
      eps(Al) =  6.47 D0
      eps(Si) =  7.55 D0
      eps(S)  =  7.33 D0
      eps(Cl) =  5.50 D0
      eps(K)  =  5.12 D0
      eps(Ca) =  6.36 D0
      eps(Ti) =  5.02 D0
      eps(Cr) =  5.67 D0
      eps(Mn) =  5.39 D0
      eps(Fe) =  7.50 D0
      eps(Ni) =  6.25 D0

*     ---------------------------------
*     Grevesse, Asplund, Sauval (2007):
*     ---------------------------------
      eps(H)  = 12.00 D0
      eps(He) = 10.93 D0
      eps(Li) =  1.10 D0  ! Lodders, Palme Gail 2009
      eps(C)  =  8.39 D0  
      eps(N)  =  7.78 D0  
      eps(O)  =  8.66 D0  
      eps(F)  =  4.56 D0
      eps(Ne) =  7.84 D0
      eps(Na) =  6.17 D0
      eps(Mg) =  7.53 D0
      eps(Al) =  6.37 D0
      eps(Si) =  7.51 D0
      eps(S)  =  7.14 D0
      eps(Cl) =  5.50 D0
      eps(K)  =  5.08 D0
      eps(Ca) =  6.31 D0
      eps(Ti) =  4.90 D0
      eps(Cr) =  5.64 D0
      eps(Mn) =  5.39 D0
      eps(Fe) =  7.45 D0
      eps(Ni) =  6.23 D0

*     ----------------------------------------------------
*     Asplund et al. (2009):
*     http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A
*     ----------------------------------------------------
      eps(H)  = 12.00    
      eps(He) = 10.93  
      eps(Li) = 1.05     
      eps(Be) = 1.38   
      eps(B)  = 2.70     
      eps(C)  = 8.43     
      eps(N)  = 7.83     
      eps(O)  = 8.69     
      eps(F)  = 4.56  
      eps(Ne) = 7.93    
      eps(Na) = 6.24  
      eps(Mg) = 7.60    
      eps(Al) = 6.45  
      eps(Si) = 7.51     
      eps(P)  = 5.41  
      eps(S)  = 7.12     
      eps(Cl) = 5.50     
      eps(Ar) = 6.40     
      eps(K)  = 5.03    
      eps(Ca) = 6.34     
      eps(Sc) = 3.15  
      eps(Ti) = 4.95     
      eps(V)  = 3.93    
      eps(Cr) = 5.64    
      eps(Mn) = 5.43  
      eps(Fe) = 7.50     
      eps(Co) = 4.99  
      eps(Ni) = 6.22    
      eps(Cu) = 4.19     
      eps(Zn) = 4.56     
      eps(Ga) = 3.04     
      eps(Ge) = 3.65     
      eps(As) = -40.   
      eps(Se) = -40.     
      eps(Br) = -40.     
      eps(Kr) = 3.25     
      eps(Rb) = 2.52    
      eps(Sr) = 2.87     
      eps(Y ) = 2.21  
      eps(Zr) = 2.58  
      eps(W ) = 0.85   

      do i=1,NELEM
        eps(i) = 10.Q0 ** (eps(i)-12.Q0)
      enddo
      eps0 = eps
  

      call INIT_CHEMISTRY
      if(condensates) call INIT_DUSTCHEM

	if(.not.allocated(linkmol)) then
		allocate(linkmol(n_mol_in))
		allocate(linkele(n_mol_in))
		linkmol=0
		linkele=0
		do j=1,n_mol_in
			call To_upper_ARCiS(mol_names_in(j),uppername)
			do i=1,NMOLE
				if(uppername.eq.cmol(i)) then
					linkmol(j)=i
				endif
			enddo
			do i=1,NELEM
				if(mol_names_in(j).eq.elnam(i)) then
					linkele(j)=i
				endif
			enddo
		enddo
	endif

	return
	end



***********************************************************************
	subroutine call_GGchem(Tin,Pin,atom_names_in,atom_abuns_in,n_atom_in,mol_names_in,mol_abuns_in,n_mol_in,
     >							MMW,condensates,atom_abuns_out)
***********************************************************************
      use PARAMETERS,ONLY: elements,abund_pick,model_dim,model_pconst,
     >                     model_struc,model_eqcond,Npoints,useDatabase,
     >                     Tfast,Tmin,Tmax,pmin,pmax,nHmin,nHmax,
     >                     abund_file,struc_file,remove_condensates
      use EXCHANGE,ONLY: chemcall,chemiter,ieqcond,ieqconditer,
     >                   itransform,nel,nat,nion,nmol,mmol,H,C,N,O,W
      use DATABASE,ONLY: NLAST

      use CHEMISTRY,ONLY: NewBackIt,NewFullIt,NewBackFac,NewPreMethod,
     >                    NewFastLevel,dispol_file,
     >                    NMOLE,NELM,m_kind,elnum,cmol,el,charge
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,muH,amu,mass,
     >                    dust_nel,dust_el,dust_nu,dust_nam,dust_mass,
     >                    dust_Vol,mass,mel
      use ARCiS_GGCHEM
      implicit none
	integer :: n_atom_in,n_mol_in,verbose,i,j
	real*8 :: Tin,Pin,atom_abuns_in(n_atom_in),mol_abuns_in(n_mol_in),MMW,atom_abuns_out(n_atom_in)
	character*40 :: atom_names_in(n_atom_in)
	character*10 :: mol_names_in(n_mol_in),uppername,elnam_UPPER

      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real :: p,rhog,rhod,dustV,nges,mges,kT,pges,mu
      real*8 :: Tg,nHges,dfdmu,dmu,ff,fold,muold,ngas,pgas
      integer :: it

      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST),tot,maxabun
	logical condensates,merk

	model_eqcond=condensates
	eps=1d-50
	do i=1,n_atom_in
		do j=1,NELEM
			if(trim(elnam(j)).eq.trim(atom_names_in(i))) then
				eps(j)=atom_abuns_in(i)
			endif
		enddo
	enddo
	eps0(1:NELEM)=eps(1:NELEM)

      muH = 0.d0
      do i=1,NELEM
        if (index(trim(elements)," "//trim(elnam(i))//" ")>0) then 
c          write(*,'(1x,a2,1x,0pF8.3,1pE12.4,0pF8.3)') 
c     &            elnam(i),12.d0+LOG10(eps(i)),eps(i),mass(i)/amu
          muH = muH + mass(i)*eps(i)
        endif  
      enddo
c      write(*,'("rho = n<H> *",1pE12.4," amu")') muH/amu
c      write(*,'("C/O =",0pF6.3)') eps(C)/eps(O)
        eldust = 0.Q0
	mu=muH

	Tg=Tin
	p=Pin*bar
	nHges = p*mu/(bk*Tg)/muH

	verbose=0
	merk=.false.

      eldust = 0.Q0

        do it=1,999
          if (model_pconst) nHges = p*mu/(bk*Tg)/muH
          if (condensates) then
            call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
          endif
          call GGCHEM(nHges,Tg,eps,.false.,verbose)

          ngas = nel
          do j=1,NELEM
            ngas = ngas + nat(j)
          enddo
          do j=1,NMOLE
            ngas = ngas + nmol(j)
          enddo
          pgas  = ngas*bk*Tg
          ff    = p-pgas
          if (it==1) then
            muold = mu
            mu = nHges/pgas*(bk*Tg)*muH
            dmu = mu-muold
            if (.not.model_pconst) exit
          else
            dfdmu = (ff-fold)/(mu-muold)
            dmu   = -ff/dfdmu
            !write(98,'(I3,99(1pE14.7))') it,muold,mu,fold,ff,dfdmu,dmu/mu
            muold = mu
            if ((dmu>0.0).or.ABS(dmu/mu)<0.7) then
              mu = muold+dmu
            else
              mu = nHges/pgas*(bk*Tg)*muH
            endif  
          endif
          fold = ff
c          print '("p-it=",i3,"  mu=",2(1pE20.12))',it,mu/amu,dmu/mu
		if (.not.ABS(dmu/mu).gt.1.E-10) exit
        enddo  


	if(condensates) then
		atom_abuns_out=1d-50
		do i=1,n_atom_in
			do j=1,NELEM
				if(trim(elnam(j)).eq.trim(atom_names_in(i))) then
					atom_abuns_out(i)=eps(j)
				endif
			enddo
		enddo
	endif

	tot=sum(nmol(1:NMOLE))+sum(nat(1:NELEM))

c		print*,Pin,Tin
c		do i=1,NMOLE
c			if(nmol(i)/tot.gt.1d-8) then
c				print*,nmol(i)/tot,trim(cmol(i))
c			endif
c		enddo


	mol_abuns_in=0d0
	do j=1,n_mol_in
c		call To_upper_ARCiS(mol_names_in(j),uppername)
c		do i=1,NMOLE
c			if(uppername.eq.cmol(i)) then
c				mol_abuns_in(j)=mol_abuns_in(j)+nmol(i)/tot
c			endif
c		enddo
c		do i=1,NELEM
c			if(mol_names_in(j).eq.elnam(i)) then
c				mol_abuns_in(j)=mol_abuns_in(j)+nat(i)/tot
c			endif
c		enddo
		i=linkmol(j)
		if(i.ne.0) then
			mol_abuns_in(j)=nmol(i)/tot
		endif
		i=linkele(j)
		if(i.ne.0) then
			mol_abuns_in(j)=nat(i)/tot
		endif
	enddo

	MMW=0d0
	do i=1,NMOLE
		MMW=MMW+nmol(i)*mmol(i)/tot
	enddo
	do i=1,NELEM
		MMW=MMW+nat(i)*mass(i)/tot
	enddo
	MMW=MMW/amu
	
      end

	subroutine To_upper_ARCiS(strin,strout)
	character(*) :: strin,strout
	integer :: i

	strout = strin
	do i = 1, len_trim(strout)
		select case(strout(i:i))
			case("a":"z")
				strout(i:i) = achar(iachar(strout(i:i))-32)
		end select
	end do
	return
	end


	subroutine call_SuperSat_ARCiS(Tg)
	use CloudModule
	use AtomsModule
	use EXCHANGE,ONLY: nat,nmol
	use DUST_DATA,ONLY: NDUST,dust_nam
	implicit none
	integer,parameter :: qp = selected_real_kind ( 33, 4931 )
	real(kind=qp) :: Sat(NDUST)
	integer i,j,l
	real*8 :: Tg
	
	nCS=10
	if(.not.allocated(SatRat)) then
		allocate(SatRat(nCS))
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

		i=0
c TiO2
		i=i+1
		CSname(i)='TiO2'
c VO
		i=i+1
		CSname(i)='VO'
c Al2O3
		i=i+1
		CSname(i)='Al2O3'
c SiO2
		i=i+1
		CSname(i)='SiO2'
c Silicates
		i=i+1
		CSname(i)='MgSiO3'
c H2O
		i=i+1
		CSname(i)='H2O'
c Fe
		i=i+1
		CSname(i)='Fe'
c FeS
		i=i+1
		CSname(i)='FeS'
c C
		i=i+1
		CSname(i)='C'
c SiC
		i=i+1
		CSname(i)='SiC'

		nCS=i
	endif


	call SUPERSAT(Tg,nat,nmol,Sat)
	SatRat(1:nCS)=1d-100
	do j=1,nCS
		l=len_trim(CSname(j))
		do i=1,NDUST
			if(CSname(j)(1:l).eq.dust_nam(i)(1:l)) then
				if(Sat(i).gt.SatRat(j)) SatRat(j)=Sat(i)
			endif
		enddo
		write(*,'(a10,se20.4)') CSname(j),SatRat(j)
	enddo

	return
	end

