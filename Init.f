c==============================================================================	
c Module for reading keywords
c==============================================================================	
	module ReadKeywords
	IMPLICIT NONE

	type SettingKey
		character*500 key1,key2,value,key
		integer nr1,nr2,key2d
		logical last
		type(SettingKey),pointer :: next
	end type SettingKey

	contains
	
	subroutine GetKeywords(firstkey)
	use GlobalSetup
	IMPLICIT NONE
	type(SettingKey),target :: firstkey
	type(SettingKey),pointer :: key
	integer ncla	! number of command line arguments
	character*1000 readline,inputfile,command
	logical readfile

	call getarg(1,inputfile)
	if(inputfile.eq.' ') inputfile='input.dat'

	open(unit=20,file=inputfile,RECL=1000)

	call output("input file: " // trim(inputfile))

	call system("cp " // trim(inputfile) // " " // trim(outputdir) // "input.dat")
	open(unit=21,file=trim(outputdir) // "input.dat",RECL=1000,access='APPEND')
	write(21,'("*** command line keywords ***")')
	
	ncla=0

	key => firstkey

	readfile=.true.
	
	goto 10
20	readfile=.false.
	close(unit=20)
10	continue
	if(readfile) then
		call ignorestar(20)
		read(20,'(a1000)',end=20,err=20) readline
	else
		call getarg(2+ncla,readline)
		if(readline(1:2).eq.'-s') then
			ncla=ncla+1
			call getarg(2+ncla,readline)
			call output("Command line argument: " // trim(readline))
			write(21,'(a)') trim(readline)
			ncla=ncla+1
		else if(readline(1:2).eq.'-o') then
			ncla=ncla+2
			goto 10
		else
			if(readline.ne.' ') then
c				try to read another command line argument
				open(unit=20,file=readline,RECL=1000)
				readfile=.true.
				ncla=ncla+1
				goto 10
			else
c				all arguments are read
				goto 30
			endif
		endif
	endif

	if(readline.eq.' ') goto 10

	allocate(key%next)
	key => key%next
	key%last=.false.
	call get_key_value(readline,key%key,key%key1,key%key2,key%value,key%nr1,key%nr2,key%key2d)

c read another command, so go back
	goto 10

30	continue
	close(unit=21)
	allocate(key%next)
	key => key%next
	key%last=.true.

	return
	end subroutine GetKeywords
	

c=========================================================================================
c This subroutine just seperates the key and value component of a string given
c key=value syntax. Key is transformed to lowercase.
c=========================================================================================
	subroutine get_key_value(line_in,key,key1,key2,value,nr1,nr2,key2d)
	use GlobalSetup
	IMPLICIT NONE
	character*1000 line_in,line
	character*500 key,key1,key2,value
	integer i,nr1,nr2,ikey1,ikey2,key2d
	
	line=line_in
	key2d=0
	if(line(1:2).eq.'2d'.and.line(4:4).eq.':') then
		read(line(3:3),*) key2d
		if(key2d.gt.n2d) n2d=key2d
		write(line,'(a)') line_in(5:1000)
	endif
	
	ikey1=index(line,'=')
	ikey2=index(line,':')
	if(ikey1.eq.0) ikey1=len_trim(line)+1

	nr1=1
	nr2=1
	if(ikey2.gt.0) then
		key1=line(1:ikey2-1)
		key=key1
		key2=line(ikey2+1:ikey1-1)
		call checknr(key1,nr1)
		call checknr(key2,nr2)
	else
		key1=line(1:ikey1-1)
		key=key1
		key2=' '
		call checknr(key1,nr1)
	endif

	value=line(index(line,'=')+1:len_trim(line))
	if(value(1:1).eq.'"'.or.value(1:1).eq."'") then
		value=value(2:len_trim(value)-1)
	endif

	do i=1,len_trim(key1)
		if(iachar(key1(i:i)).ge.65.and.iachar(key1(i:i)).le.90) then
			key1(i:i)=achar(iachar(key1(i:i))+32)
		endif
	enddo
	do i=1,len_trim(key2)
		if(iachar(key2(i:i)).ge.65.and.iachar(key2(i:i)).le.90) then
			key2(i:i)=achar(iachar(key2(i:i))+32)
		endif
	enddo

	return
	end subroutine get_key_value
	
	
	subroutine checknr(key,nr)
	IMPLICIT NONE
	character*500 key
	integer nr,i,n
	
	n=len_trim(key)
	i=n
1	read(key(i:n),*,err=2) nr
	i=i-1
	if(i.eq.0) goto 2
	goto 1
2	continue
	if(i.eq.n) then
		nr=1
	else
		read(key(i+1:n),*,err=3) nr	
	endif
3	continue
	key=key(1:i)
	return
	end subroutine checknr
	



	subroutine CountStuff(firstkey)
	use GlobalSetup
	IMPLICIT NONE
	type(SettingKey),target :: firstkey
	type(SettingKey),pointer :: key
	integer i,j,ncia0,n
	character*500 homedir,h2h2file,h2hefile,h2ch4file
	character*10 names(59)
	logical existh2h2,existh2he,existh2ch4


	nmol=1
	nobs=0
	ncia=0
	nclouds=0
	n_points=0
	n_ret=0
	n_instr=0
	j=0
	mixratfile=.false.
	fcloud_default=1d0

	i2d=0

	nr=20
	key => firstkey
	do while(.not.key%last)
		select case(key%key1)
			case("nr")
				read(key%value,*) nr
		end select
		key=>key%next
	enddo

	key => firstkey
	do while(.not.key%last)
		select case(key%key1)
			case("i2d")
				read(key%value,*) i2d
			case("obs")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nobs) nobs=key%nr1
			case("instrument")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.n_instr) n_instr=key%nr1
			case("cia")
				if(key%key2.eq.' ') then
					read(key%value,*) do_cia
				else
					if(key%nr1.eq.0) key%nr1=1
					if(key%nr2.eq.0) key%nr2=1
					if(key%nr1.gt.ncia) ncia=key%nr1
				endif
			case("mixratfile")
				read(key%value,*) mixratfile
			case("tpfile")
				read(key%value,'(a)') TPfile
			case("cloud")
c				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nclouds) nclouds=key%nr1
			case("fcloud")
				read(key%value,*) fcloud_default
			case("point")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.n_points) n_points=key%nr1
			case("retpar","fitpar")
				if(key%key2.eq.'keyword') then
					if(key%value.eq.'tprofile') then
						n_ret=n_ret+nr+2
					else
						n_ret=n_ret+1
					endif
				endif
			case default
				do i=1,59
					if(key%key.eq.molname(i)) then
						if(i.gt.nmol) nmol=i
						j=j+1
					endif
				enddo
		end select
		key=>key%next
	enddo

	if(mixratfile) then
		open(unit=30,file=TPfile,RECL=6000)
		read(30,*) n
		read(30,*) names(1:n)
		do j=1,n
			do i=1,59
				if(names(j).eq.molname(i)) then
					if(i.gt.nmol) nmol=i
				endif
			enddo
		enddo
	endif

	if(do_cia) nmol=max(nmol,48)

	allocate(mixrat(nmol))
	allocate(includemol(nmol))
	allocate(Cloud(max(nclouds,1)))
	allocate(XeqCloud(nr,max(nclouds,1)))
	allocate(XeqCloud_old(nr,max(nclouds,1)))
	allocate(XCloud(nr,max(nclouds,1)))
	allocate(P_point(max(n_points,1)))
	allocate(T_point(max(n_points,1)))
	allocate(RetPar(max(n_ret,1)))
	allocate(ObsSpec(max(nobs,1)))
	allocate(Tin(nr))
	allocate(instrument(max(n_instr,1)))
	allocate(instr_ntrans(max(n_instr,1)))

	ncia0=0
	existh2h2=.false.
	existh2he=.false.
	existh2ch4=.false.
	if(do_cia) then
		call getenv('HOME',homedir)
c find H2-H2 cia file
		h2h2file=trim(homedir) // '/HITRAN/H2-H2_2011.cia'
		inquire(file=h2h2file,exist=existh2h2)
		if(existh2h2) then
			ncia0=ncia0+1
		else
			h2h2file=trim(homedir) // '/HITRAN/CIA/H2-H2_2011.cia'
			inquire(file=h2h2file,exist=existh2h2)
			if(existh2h2) then
				ncia0=ncia0+1
			else
				h2h2file=trim(homedir) // '/ARCiS/Data/CIA/H2-H2_2011.cia'
				inquire(file=h2h2file,exist=existh2h2)
				if(existh2h2) then
					ncia0=ncia0+1
				endif
			endif
		endif
c find H2-He cia file
		h2hefile=trim(homedir) // '/HITRAN/H2-He_2011.cia'
		inquire(file=h2hefile,exist=existh2he)
		if(existh2he) then
			ncia0=ncia0+1
		else
			h2hefile=trim(homedir) // '/HITRAN/CIA/H2-He_2011.cia'
			inquire(file=h2hefile,exist=existh2he)
			if(existh2he) then
				ncia0=ncia0+1
			else
				h2hefile=trim(homedir) // '/ARCiS/Data/CIA/H2-He_2011.cia'
				inquire(file=h2hefile,exist=existh2he)
				if(existh2he) then
					ncia0=ncia0+1
				endif
			endif
		endif
c find H2-CH4 cia file
c		h2ch4file=trim(homedir) // '/HITRAN/H2-CH4_eq_2011.cia'
c		inquire(file=h2ch4file,exist=existh2ch4)
c		if(existh2ch4) then
c			ncia0=ncia0+1
c		else
c			h2ch4file=trim(homedir) // '/HITRAN/CIA/H2-CH4_eq_2011.cia'
c			inquire(file=h2ch4file,exist=existh2ch4)
c			if(existh2ch4) then
c				ncia0=ncia0+1
c			else
c				h2ch4file=trim(homedir) // '/CIA/H2-CH4_eq_2011.cia'
c				inquire(file=h2ch4file,exist=existh2ch4)
c				if(existh2ch4) then
c					ncia0=ncia0+1
c				endif
c			endif
c		endif
	endif

	allocate(CIA(max(ncia+ncia0,1)))

	if(existh2h2) then
		ncia=ncia+1
		CIA(ncia)%filename=h2h2file
	endif
	if(existh2he) then
		ncia=ncia+1
		CIA(ncia)%filename=h2hefile
	endif
	if(existh2ch4) then
		ncia=ncia+1
		CIA(ncia)%filename=h2ch4file
	endif

	call output('Number of molecules:       ' // int2string(j,'(i4)'))
	call output('Number of clouds:          ' // int2string(nclouds,'(i4)'))
	call output('Number of collision pairs: ' // int2string(ncia,'(i4)'))
	call output('Number of observations:    ' // int2string(nobs,'(i4)'))
	

	return
	end subroutine CountStuff
	
	end module ReadKeywords
	
c==============================================================================	
c==============================================================================	


	subroutine Init()
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey),pointer :: key,first
	type(SettingKey) keyret
	integer i,j,omp_get_max_threads,omp_get_thread_num
	real*8 tot,tot2,theta,Planck
	real*8,allocatable :: var(:),dvar(:)
	character*1000 line

	allocate(key)
	first => key

	n2d=0

	call GetKeywords(key)
c Count the number of zones, particles, and stars
c allocate the arrays
	key => first%next

	do_cia=.true.
	call CountStuff(key)

	call SetDefaults

	n_ret=0

	key => first%next
	do while(.not.key%last)
		call ReadAndSetKey(key)
		key => key%next
	enddo

	if(randomseed) call system_clock(idum0)
	idum=-idum0
	
#ifdef USE_OPENMP
	j=omp_get_max_threads()
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(j,idum0)
!$OMP DO SCHEDULE(STATIC,j)
	do i=1,j
		idum=-idum0-omp_get_thread_num()
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
#endif

	if(n_ret.gt.0) n_ret=n_ret+RetPar(n_ret)%n-1

	if(retrieval) then

	key => first%next

	do while(.not.key%last)

		do i=1,n_ret
			line=trim(RetPar(i)%keyword) // "=0d0"
			call get_key_value(line,keyret%key,keyret%key1,keyret%key2,keyret%value,keyret%nr1,keyret%nr2,keyret%key2d)
			if(trim(keyret%key1).eq.trim(key%key1).and.trim(keyret%key2).eq.trim(key%key2).and.
     &		   keyret%nr1.eq.key%nr1.and.keyret%nr2.eq.key%nr2.and.keyret%key2d.eq.key%key2d) then
				read(key%value,*) RetPar(i)%x0
			endif
		enddo

		key => key%next
	
	enddo
	
	endif

	call ConvertUnits()

c	condensates=(condensates.or.cloudcompute)

	if(specresdust.gt.specres) specresdust=specres

	allocate(gg(ng),wgg(ng))
	call gauleg(0d0,1d0,gg,wgg,ng)
	
	call InitFreq()

	if(nphase.le.0) then
		nphase=2d0*pi/asin(Rstar/Dplanet)
		if(nphase.gt.90) then
			call output("Adjusting number of phase angles to 90")
			nphase=90
		endif
	endif

	if(particledir.eq.' ') particledir=trim(outputdir)

	if(opacitydir(len_trim(opacitydir)-1:len_trim(opacitydir)).ne.'/') then
		opacitydir=trim(opacitydir) // '/'
	endif
	if(opacitymode) then
		call InitOpacityMode()
	else
		allocate(dens(nr))
		allocate(Ndens(nr))
		allocate(R(nr+1))
		allocate(T(nr))
		allocate(P(nr))
		allocate(Hp(nr))
		allocate(nabla_ad(nr))
		allocate(grav(nr))
		allocate(MMW(nr))
		allocate(didcondens(nr))
		allocate(mixrat_r(nr,nmol))
		allocate(mixrat_old_r(nr,nmol))
		allocate(cloud_dens(nr,max(nclouds,1)))
		call InitDens()
		call InitObs()
		do i=1,nclouds
			call output("==================================================================")
			call output("Setting up cloud: " // trim(int2string(i,'(i4)')))
			if(.not.useDRIFT) then
				call SetupPartCloud(i)
				allocate(Cloud(i)%w(Cloud(i)%nr))
				if(Cloud(i)%species.eq.' ') call NameCloudSpecies(Cloud(i)%standard,Cloud(i)%species)
			else
				Cloud(i)%nr=nr
				allocate(Cloud(i)%w(Cloud(i)%nr))
				allocate(Cloud(i)%frac(nr,20))
				Cloud(i)%species="DRIFT"
				Cloud(i)%standard="DRIFT"
			endif
		enddo
	endif

	metallicity0=metallicity

	allocate(Cabs(nr,nlam,ng))
	allocate(Csca(nr,nlam))
	do i=1,360
		theta=(real(i)-0.5d0)*pi/180d0
		sintheta(i)=sin(theta)
		costheta(i)=cos(theta)
	enddo
	do i=1,180
		Rayleigh%F11(i)=(1d0+costheta(i)**2)/2d0
	enddo
	tot=0d0
	tot2=0d0
	do i=1,180
		tot=tot+Rayleigh%F11(i)*sintheta(i)
		tot2=tot2+sintheta(i)
	enddo
	Rayleigh%F11=Rayleigh%F11*tot2/tot
	Rayleigh%IF11=0d0
	j=1
	Rayleigh%IF11(j)=sintheta(j)*Rayleigh%F11(j)
	do j=2,180
		Rayleigh%IF11(j)=Rayleigh%IF11(j-1)+sintheta(j)*Rayleigh%F11(j)
	enddo

	allocate(Fstar(nlam))
	call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam,starfile)
	Fstar=Fstar*pi*Rstar**2

	call output("==================================================================")

	if(compute_opac) then
		call ReadData()
	else
		do i=1,nmol
			if(includemol(i)) call InitReadOpacityFITS(i)
		enddo
	endif

	call ReadDataCIA()

	if(Pform.lt.0d0) Pform=-Pform*82.05736*1.01325*Tform/(mp*2.2*Avogadro)

	if(retrieval) then
		do i=1,n_ret
			if(RetPar(i)%x0.lt.-1d150) then
				if(RetPar(i)%keyword(1:6).eq.'tvalue') then
					RetPar(i)%x0=0d0
				else
					if(RetPar(i)%logscale) then
						RetPar(i)%x0=sqrt(RetPar(i)%xmax*RetPar(i)%xmin)
					else
						RetPar(i)%x0=0.5d0*(RetPar(i)%xmax+RetPar(i)%xmin)
					endif
				endif
			endif
			if(RetPar(i)%dx.lt.0d0) then
				if(RetPar(i)%logscale) then
					RetPar(i)%dx=5d0*(RetPar(i)%xmax/RetPar(i)%xmin)
				else
					RetPar(i)%dx=5d0*(RetPar(i)%xmax-RetPar(i)%xmin)
				endif
			endif
		enddo			
	endif

	if(dochemistry) call init_GGchem()

	if(iWolk.gt.0) then
		open(unit=50,file=trim(outputdir) // "/Wolk.dat",RECL=6000)
		allocate(var(n_ret),dvar(n_ret))
		do i=1,iWolk
			read(50,*) j,tot,var(1:n_ret)
		enddo
		call MapRetrieval(var,dvar)
		close(unit=50)
	endif
	call output("==================================================================")
	
	return
	end


	subroutine ReadAndSetKey(key)
	use GlobalSetup
	use ReadKeywords
	use CloudModule
	IMPLICIT NONE
	type(SettingKey) key
	integer i

	if(key%key2d.eq.i2d.or.key%key2d.eq.0) then

	select case(key%key1)
		case("nr")
c is already set in CountStuff
c			read(key%value,*) nr
		case("nrcloud","nrhr")
			read(key%value,*) nr_cloud
			nr_cloud=max(1,nr_cloud/nr)
		case("mp")
			Mp_from_logg=.false.
			read(key%value,*) Mplanet
		case("rp")
			read(key%value,*) Rplanet
		case("pp")
			read(key%value,*) Pplanet
		case("loggp")
			Mp_from_logg=.true.
			read(key%value,*) loggPlanet
		case("rstar")
			read(key%value,*) Rstar
		case("mstar")
			read(key%value,*) Mstar
		case("tstar")
			read(key%value,*) Tstar
		case("starfile")
c starfile should be in W/(m^2 Hz) at the stellar surface
			starfile=key%value
		case("logg")
			read(key%value,*) logg
		case("dp","dplanet")
			read(key%value,*) Dplanet
		case("retrieval")
			read(key%value,*) retrieval
		case("contrib","computecontrib")
			read(key%value,*) computecontrib
		case("outputopacity","writeopacity")
			read(key%value,*) outputopacity
		case("nphase")
			read(key%value,*) nphase
		case("cia")
			call ReadCIA(key)
		case("cutoff","cutoff_lor")
			read(key%value,*) cutoff_lor
		case("cutoff_abs")
			read(key%value,*) cutoff_abs
		case("lam")
			if(key%nr1.eq.1) read(key%value,*) lam1
			if(key%nr1.eq.2) read(key%value,*) lam2
		case("lmin")
			read(key%value,*) lam1
		case("lmax")
			read(key%value,*) lam2
		case("pmin")
			read(key%value,*) pmin
		case("pmax")
			read(key%value,*) pmax
		case("pcloud","psimplecloud")
			read(key%value,*) psimplecloud
		case("tmin")
			read(key%value,*) Tmin
		case("tmax")
			read(key%value,*) Tmax
		case("eps","epsck")
			read(key%value,*) epsCk
		case("epslines","eps_lines")
			read(key%value,*) eps_lines
		case("factrw")
			read(key%value,*) factRW
		case("maxtau","max_tau")
			read(key%value,*) maxtau
		case("specres")
			read(key%value,*) specres
		case("specresdust")
			read(key%value,*) specresdust
		case("specresfile")
			specresfile=trim(key%value)
		case("particledir","dirparticle")
			particledir=trim(key%value)
		case("tpfile")
			TPfile=key%value
		case("mixratfile")
			read(key%value,*) mixratfile
		case("gridtpfile")
			read(key%value,*) gridTPfile
		case("maxt","maxtprofile")
			read(key%value,*) maxTprofile
		case("tp")
			read(key%value,*) TP0
		case("dtp")
			read(key%value,*) dTP
		case("gamma","gammat")
			if(key%nr1.eq.1) read(key%value,*) gammaT1
			if(key%nr1.eq.2) then
				gamma_equal=.false.
				read(key%value,*) gammaT2
			endif
		case("kappa","kappat")
			read(key%value,*) kappaT
		case("alpha","alphat")
			read(key%value,*) alphaT
		case("beta","betat")
			read(key%value,*) betaT
		case("tmix","twind")
			read(key%value,*) twind
			if(twind.gt.0d0) n2d=max(2,n2d)
		case("partprofile","par_tprofile")
			read(key%value,*) par_tprofile
		case("ng")
			read(key%value,*) ng
		case("distance")
			read(key%value,*) distance
		case("hitemp")
			read(key%value,*) HITEMP
		case("cloud")
			call ReadCloud(key)
		case("mixrathaze","haze")
			read(key%value,*) mixratHaze
		case("phaze")
			read(key%value,*) PHaze
		case("dphaze")
			read(key%value,*) dPHaze
		case("kappahaze")
			read(key%value,*) kappaHaze
		case("scattering")
			read(key%value,*) scattering
		case("scattstar","starscatt")
			read(key%value,*) scattstar
		case("compute","compute_opac","computeopac")
			read(key%value,*) compute_opac
		case("opacitymode")
			read(key%value,*) opacitymode
		case("np")
			read(key%value,*) nPom
		case("nt")
			read(key%value,*) nTom
		case("opacitydir")
			opacitydir=trim(key%value)
		case("computet")
			read(key%value,*) computeT
		case("teffp","tplanet")
			read(key%value,*) TeffP
		case("maxiter")
			read(key%value,*) maxiter
		case("epsiter")
			read(key%value,*) epsiter
		case("chemistry")
			read(key%value,*) dochemistry
		case("fastchem","fast_chem")
			read(key%value,*) fast_chem
		case("metallicity")
			read(key%value,*) metallicity
		case("coratio")
			read(key%value,*) COratio
		case("tiscale")
			read(key%value,*) TiScale
		case("enhancecarbon")
			read(key%value,*) enhancecarbon
		case("condensates")
			read(key%value,*) condensates
		case("cloudcompute")
			read(key%value,*) cloudcompute
		case("usedrift")
			read(key%value,*) useDRIFT
		case("mixp")
			read(key%value,*) mixP
		case("sinkz")
			read(key%value,*) sinkZ
		case("alphaz")
			read(key%value,*) alphaZ
		case("nspike")
			read(key%value,*) nspike
		case("nphot")
			read(key%value,*) Nphot0
		case("domccompute","mccompute","mccomputet")
			read(key%value,*) doMCcompute
		case("point")
			call ReadPoint(key)
		case("retpar","fitpar")
			call ReadRetrieval(key)
		case("obs")
			call ReadObsSpec(key)
		case("instrument")
			call ReadInstrument(key)
		case("chimax","chi2max")
			read(key%value,*) chimax
		case("npop")
			read(key%value,*) npop
		case("ngen")
			read(key%value,*) ngen
		case("gene_cross")
			read(key%value,*) gene_cross
		case("resume_nest")
			read(key%value,*) resume_multinest
		case("fmultinest","efr")
			read(key%value,*) f_multinest
		case("tolmultinest","ftol")
			read(key%value,*) tol_multinest
		case("retrievaltype")
			read(key%value,*) retrievaltype
		case("tvalue")
			read(key%value,*) Tin(key%nr1)
		case("retrieve_profile")
			read(key%value,*) retrieve_profile
		case("faircoverage")
			read(key%value,*) faircoverage
		case("speclimits")
			read(key%value,*) speclimits
		case("adiabatic","adiabatic_tprofile")
			read(key%value,*) adiabatic_tprofile
		case("tform")
			read(key%value,*) Tform
		case("pform")
			read(key%value,*) Pform
		case("fenrich","f_enrich","enrich","fdry","f_dry")
			read(key%value,*) f_dry
		case("fwet","f_wet")
			read(key%value,*) f_wet
		case("scalefe","scale_fe","f_fe")
			read(key%value,*) scale_fe
		case("tchem")
			read(key%value,*) Tchem
		case("pchem")
			read(key%value,*) Pchem
		case("chemabun","ptchemabun")
			read(key%value,*) PTchemAbun
		case("rnuc","r_nuc")
			read(key%value,*) r_nuc
		case("rhoform","densform")
			read(key%value,*) Pform
			Pform=-Pform
		case("makeai")
			read(key%value,*) domakeai
		case("pew","postequalweights")
			read(key%value,*) dopostequalweights
		case("mapcoratio")
			read(key%value,*) mapCOratio
		case("nai")
			read(key%value,*) nai
		case("maxchemtime")
			read(key%value,*) maxchemtime
		case("idum","seed")
			read(key%value,*) idum0
		case("randomseed")
			read(key%value,*) randomseed
		case("fcloud")
			read(key%value,*) fcloud_default
		case("singlecloud","exactcoverage")
			read(key%value,*) singlecloud
		case("coagulation")
			read(key%value,*) coagulation
		case("parameterfile","planetparameterfile")
			planetparameterfile=trim(key%value)
		case("planetname")
			read(key%value,*) planetname
			call ReadPlanetName
		case("trend_compute","dotrend")
			read(key%value,*) trend_compute
		case("i2d")
			read(key%value,*) i2d
		case("n2d")
			read(key%value,*) i
			n2d=max(i,n2d)
		case("iwolk")
			read(key%value,*) iWolk
		case default
			do i=1,59
				if(key%key.eq.molname(i)) then
					read(key%value,*) mixrat(i)
					includemol(i)=.true.
					goto 1
				endif
			enddo
			call output("Keyword not recognised: " // trim(key%key1))
			stop
1			continue
	end select

	endif
	
	return
	end

	subroutine ConvertUnits()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 tot
	integer i
	
	Rplanet=Rplanet*Rjup
	Mplanet=Mplanet*Mjup

	if(Mp_from_logg) then
		Mplanet=(Rplanet**2)*(10d0**(loggPlanet))/Ggrav
		call output("Planet mass: " // dbl2string(Mplanet/Mjup,'(f8.3)') // " Mjup")
	else
		loggPlanet=log10(Ggrav*Mplanet/(Rplanet**2))
		call output("Planet logg: " // dbl2string(loggPlanet,'(f8.3)'))
	endif	
	if(gamma_equal.or.gammaT2.lt.0d0) then
		gammaT2=gammaT1
		gamma_equal=.true.
	endif

	Rstar=Rstar*Rsun
	Mstar=Mstar*Msun
	Dplanet=Dplanet*AU
	
	lam1=lam1*micron
	lam2=lam2*micron
	
	r_nuc=r_nuc*micron

	distance=distance*parsec

c	tot=0d0
c	do i=1,nmol
c		if(mixrat(i).gt.0d0) tot=tot+mixrat(i)
c	enddo
c	if(tot.gt.1d0) then
c		call output("Summed mixing ratio above 1. Renormalizing.")
c		mixrat=mixrat/tot
c	endif

	return
	end


	subroutine InitDens()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,n
	real*8 g,dp,dz,P0(nr),T0(nr),pp,tt,mr0(nr,nmol),mm(nmol),yp1,ypn
	real*8,allocatable :: y2(:)
	character*10 names(nmol)
	integer imol(nmol)
	
	do i=1,nr
		mixrat_r(i,1:nmol)=mixrat(1:nmol)
	enddo
	
	if(mixratfile) then
		open(unit=20,file=TPfile,RECL=1000)
		read(20,*) n
		read(20,*) names(1:n)
		do j=1,n
			do i=1,59
				if(names(j).eq.molname(i)) then
					imol(j)=i
					includemol(imol(j))=.true.
				endif
			enddo
		enddo
		i=1
1		read(20,*,err=1,end=2) pp,tt,mm(1:n)
		if(i.gt.nr) then
			call output("Increase number of radial points for this file")
			stop
		endif
		P0(i)=pp
		T0(i)=tt
		mr0(i,1:n)=mm(1:n)
		i=i+1
		goto 1
2		close(unit=20)
		nr=i-1
	else if(gridTPfile) then
		open(unit=20,file=TPfile)
		i=1
3		read(20,*,err=3,end=4) pp,tt
		if(i.gt.nr) then
			call output("Increase number of radial points for this file")
			stop
		endif
		P0(i)=pp
		T0(i)=tt
		i=i+1
		goto 3
4		close(unit=20)
		nr=i-1
	else
		do i=1,nr
			P0(i)=exp(log(Pmin)+log(Pmax/Pmin)*real(i-1)/real(nr-1))
		enddo
		if(TPfile.ne.' ') then
			call regridlog(TPfile,P0,T0,nr)
		else if(n_points.gt.1) then
			P_point=log10(P_point)
			T_point=log10(T_point)
			allocate(y2(n_points))
			yp1=0d0
			ypn=0d0
			call spline(P_point,T_point,n_points,yp1,ypn,y2)
			do i=1,nr
				call splint(P_point,T_point,y2,n_points,log10(P0(i)),T0(i))
				T0(i)=10d0**T0(i)
			enddo
			deallocate(y2)
			P_point=10d0**P_point
			T_point=10d0**T_point
		else
			do i=1,nr
				T0(i)=10d0**(log10(TP0)+dTP*log10(P0(i)))
			enddo
		endif
	endif

c	if(retrieval) then
c		do i=1,nr
c			if(T0(i).gt.Tmax) T0(i)=Tmax
c			if(T0(i).lt.Tmin) T0(i)=Tmin
c		enddo
c	endif

	do i=1,nr
		T(i)=T0(nr+1-i)
		P(i)=P0(nr+1-i)
	enddo
c	if(par_tprofile) call ComputeParamT(T)
	do i=1,nr
		if(retrieve_profile) then
			T(i)=T(i)+Tin(i)
			if(T(i).gt.maxTprofile) T(i)=maxTprofile
			if(T(i).lt.minTprofile) T(i)=minTprofile
		endif
		if(mixratfile) then
			do j=1,n
				mixrat_r(i,imol(j))=mr0(nr+1-i,j)
			enddo
		endif
		if(IsNaN(T(i))) T(i)=sqrt(minTprofile*maxTprofile)
	enddo

	return
	end	


	subroutine InitOpacityMode()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,n
	real*8 g,dp,dz,P0(nr),T0(nr),pp,tt,mr0(nr,nmol),mm(nmol)
	character*10 names(nmol)
	integer imol(nmol)

	nr=nPom*nTom
	maxtau=-1d0

	allocate(dens(nr))
	allocate(Ndens(nr))
	allocate(R(nr+1))
	allocate(T(nr))
	allocate(P(nr))
	allocate(Hp(nr))
	allocate(nabla_ad(nr))
	allocate(grav(nr))
	allocate(mixrat_r(nr,nmol))
	allocate(mixrat_old_r(nr,nmol))
	allocate(cloud_dens(nr,max(nclouds,1)))
	do i=1,nclouds
		call output("==================================================================")
		call output("Setting up cloud: " // trim(int2string(i,'(i4)')))
		call SetupPartCloud(i)
		allocate(Cloud(i)%w(Cloud(i)%nr))
		if(Cloud(i)%species.eq.' ') call NameCloudSpecies(Cloud(i)%standard,Cloud(i)%species)
	enddo

	do i=1,nr
		mixrat_r(i,1:nmol)=mixrat(1:nmol)
	enddo

	n=0
	do i=1,nPom
		do j=1,nTom	
			n=n+1
			P(n)=10d0**(log10(Pmin)+log10(Pmax/Pmin)*real(i-1)/real(nPom-1))
			T(n)=10d0**(log10(Tmin)+log10(Tmax/Tmin)*real(j-1)/real(nTom-1))
		enddo
	enddo

	return
	end	


	subroutine SetDefaults()
	use GlobalSetup
	use Constants
	use CloudModule
	IMPLICIT NONE
	integer i
	character*100 homedir
	
	idum0=42
	randomseed=.true.

	iWolk=0
	
	particledir='./Particles/'

	Mplanet=1d0
	Rplanet=1d0
	Pplanet=10d0
	loggPlanet=2.5d0
	Mp_from_logg=.false.
	
	Tstar=5777d0
	Rstar=1d0
	Mstar=1d0
	Dplanet=1d0
	logg=4.5d0
	
	lam1=1d0
	lam2=15d0
	specres=10d0
	specresdust=10d0
	
	trend_compute=.false.

	mixrat=0d0
	includemol=.false.
	par_tprofile=.false.
	
	nphase=0
	Nphot0=2500

	dochemistry=.false.
	fast_chem=.false.
	metallicity=0d0
	condensates=.true.
	COratio=0.55
	TiScale=1d0
	enhancecarbon=.false.
	mixP=0d0
	mixratHaze=0d0
	PHaze=1d0
	dPHaze=1d6
	kappaHaze=0d0
	Psimplecloud=1d9
	coagulation=.true.
	singlecloud=.false.
	
	twind=-1d0
	
	PRplanet=10d0

	TPfile=' '
	mixratfile=.false.
	Tin=0d0

	starfile=' '
	
	specresfile=' '
	
	speclimits=.true.
	
	sinkZ=.false.
	alphaZ=1d0
	
	faircoverage=.false.

	Tform=-10d0
	Pform=1d0
	f_dry=1d0
	f_wet=1d0
	scale_fe=1d0
	
	r_nuc=1d-3
	
	PTchemAbun=.false.
	Tchem=500d0
	Pchem=1d0
	
	adiabatic_tprofile=.false.

	domakeai=.false.
	nai=1000
	dopostequalweights=.false.

	computecontrib=.false.
	
	instrument="ARIEL"
	instr_ntrans=1d0

	nr_cloud=10

	do i=1,nclouds
		Cloud(i)%P=1d-4
		Cloud(i)%dP=10d0
		Cloud(i)%s=2d0
		Cloud(i)%column=0d0
		Cloud(i)%tau=-1d0
		Cloud(i)%lam=0.55d0
		Cloud(i)%file=' '
		Cloud(i)%Kzzfile=' '
		Cloud(i)%standard='ASTROSIL'
		Cloud(i)%nr=50
		Cloud(i)%nsubgrains=1
		Cloud(i)%amin=0.1
		Cloud(i)%amax=100d0
		Cloud(i)%fmax=0d0
		Cloud(i)%porosity=0d0
		Cloud(i)%fcarbon=0d0
		Cloud(i)%blend=.true.
		Cloud(i)%reff=1d0
		Cloud(i)%veff=0.1
		Cloud(i)%coverage=fcloud_default
		Cloud(i)%ptype='COMPUTE'
		Cloud(i)%frain=1d0
		Cloud(i)%species=''
		Cloud(i)%haze=.false.
		Cloud(i)%tmix=300d0
		Cloud(i)%betamix=2.2
		Cloud(i)%Kscale=1d0
		Cloud(i)%Kzz=1d8
		Cloud(i)%Kzz_pow=0d0
		Cloud(i)%Sigmadot=1d-17
	enddo
	cloudcompute=.false.
	useDRIFT=.false.
	nspike=0

	retrieval=.false.
	do i=1,n_ret
		RetPar(i)%x0=-1d200
		RetPar(i)%dx=-1d0
		RetPar(i)%logscale=.false.
		RetPar(i)%squarescale=.false.
		RetPar(i)%opacitycomp=.true.
		RetPar(i)%increase=.false.
	enddo
	npop=30
	ngen=0
	gene_cross=.false.
	resume_multinest=.false.
	f_multinest=0.3d0
	tol_multinest=0.5d0
	retrievaltype='MN'

	do i=1,nobs
		ObsSpec(i)%beta=-1d0
		ObsSpec(i)%scale=-1d-1
		ObsSpec(i)%spec=.true.
		ObsSpec(i)%i2d=0
	enddo

	computeT=.false.
	doMCcompute=.true.
	TeffP=600d0
	outputopacity=.false.

	Pmin=1d-6
	Pmax=1d+3

	call getenv('HOME',homedir) 

	HITRANdir=trim(homedir) // '/HITRAN/'
	HITEMPdir=trim(homedir) // '/HITEMP/'

	HITEMP=.false.

	planetparameterfile=trim(homedir) // '/ARCiS/Data/allplanets-ascii.txt'
	planetname=' '

	ng=25
	epsCk=0.25d0
	
	eps_lines=0d0
	maxtau=50d0
	factRW=10d0
	
	distance=10d0
	
	cutoff_abs=1d200
	cutoff_lor=1d200
	
	gridTPfile=.false.
	
	scattering=.true.
	scattstar=.false.
	
	compute_opac=.true.
	
	opacitymode=.false.
	opacitydir=trim(homedir) // '/ARCiS/Data/Opacities'

	nTom=100
	nPom=50
	Tmin=71d0
	Tmax=2900d0
	
	TP0=600d0
	dTP=0.1
	
	gammaT1=1.58e-1
	gammaT2=1.58e-1
	gamma_equal=.true.
	kappaT=3d-4
	betaT=1d0
	alphaT=1d0
	
	retrieve_profile=.false.
	chimax=1d0
	
	maxTprofile=1d6
	
	maxiter=6
	epsiter=1d-2

	maxchemtime=1d200
	
	return
	end

	
	subroutine InitObs()
	use GlobalSetup
	use Constants
	IMPLICIT NONE

c number of cloud/nocloud combinations
	ncc=2**nclouds
	if(singlecloud) ncc=1+nclouds
	allocate(docloud(ncc,nclouds))
	allocate(cloudfrac(ncc))
	allocate(flux(0:ncc,nlam))
	allocate(obsA(0:ncc,nlam))
	allocate(tau1depth(ncc,nlam))
	allocate(cloudtau(ncc,nlam))
	allocate(phase(nphase,0:ncc,nlam))

	return
	end
	
	
	
	subroutine ReadPoint(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	i=key%nr1
	
	select case(key%key2)
		case("p")
			read(key%value,*) P_point(i)
		case("t")
			read(key%value,*) T_point(i)
		case default
			call output("Keyword not recognised: " // trim(key%key2))
	end select
	
	return
	end

	subroutine ReadRetrieval(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i,j
	i=n_ret
	
	select case(key%key2)
		case("keyword","parameter")
			if(i.eq.0) then
				n_ret=n_ret+1
			else
				n_ret=n_ret+RetPar(i)%n
			endif
			i=n_ret
			read(key%value,*) RetPar(i)%keyword
			if(RetPar(i)%keyword.eq.'tprofile') then
				RetPar(i)%n=nr+2
				retrieve_profile=.true.
			else
				RetPar(i)%n=1
			endif
			if(RetPar(i)%keyword.eq.'tprofile') then
				do j=1,nr
					RetPar(i+j-1)%keyword='tvalue' // trim(int2string(j,'(i0.3)'))
				enddo
				RetPar(i+nr)%keyword='tp'
				RetPar(i+nr+1)%keyword='dtp'
			endif
		case("min","xmin")
			read(key%value,*) RetPar(i)%xmin
			do j=1,RetPar(i)%n
				RetPar(i+j-1)%xmin=RetPar(i)%xmin
			enddo
			if(RetPar(i)%keyword.eq.'tvalue001') then
				minTprofile=RetPar(i)%xmin
				RetPar(i+RetPar(i)%n-1)%xmin=-0.2
			endif
		case("max","xmax")
			read(key%value,*) RetPar(i)%xmax
			do j=1,RetPar(i)%n
				RetPar(i+j-1)%xmax=RetPar(i)%xmax
			enddo
			if(RetPar(i)%keyword.eq.'tvalue001') then
c				do j=1,nr
c					RetPar(i+j-1)%xmin=-RetPar(i)%xmax
c				enddo
				maxTprofile=RetPar(i)%xmax
				RetPar(i+RetPar(i)%n-1)%xmax=0.2
			endif
		case("init","x")
			read(key%value,*) RetPar(i)%x0
			do j=1,RetPar(i)%n
				RetPar(i+j-1)%x0=RetPar(i)%x0
			enddo
			if(RetPar(i)%keyword.eq.'tvalue001') then
				RetPar(i+RetPar(i)%n-1)%x0=0d0
			endif
		case("spread","width","dx")
			read(key%value,*) RetPar(i)%dx
			do j=1,RetPar(i)%n
				RetPar(i+j-1)%dx=RetPar(i)%dx
			enddo
			if(RetPar(i)%keyword.eq.'tvalue001') then
				RetPar(i+RetPar(i)%n-1)%dx=1d0
			endif
		case("log","logscale")
			read(key%value,*) RetPar(i)%logscale
			do j=1,RetPar(i)%n
				RetPar(i+j-1)%logscale=RetPar(i)%logscale
			enddo
			if(RetPar(i)%keyword.eq.'tvalue001') then
				do j=1,RetPar(i)%n
					RetPar(i+j-1)%logscale=.false.
				enddo
				read(key%value,*) RetPar(i+RetPar(i)%n-2)%logscale
			endif
		case("square","squarescale")
			read(key%value,*) RetPar(i)%squarescale
			do j=1,RetPar(i)%n
				RetPar(i+j-1)%squarescale=RetPar(i)%squarescale
			enddo
		case("opacity","opacitycomp")
			read(key%value,*) RetPar(i)%opacitycomp
			do j=1,RetPar(i)%n
				RetPar(i+j-1)%opacitycomp=RetPar(i)%opacitycomp
			enddo
		case("increase")
			read(key%value,*) RetPar(i)%increase
			do j=1,RetPar(i)%n
				RetPar(i+j-1)%increase=RetPar(i)%increase
			enddo
		case default
			call output("Keyword not recognised: " // trim(key%key2))
	end select
	
	return
	end

	subroutine ReadObsSpec(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	i=key%nr1
	
	select case(key%key2)
		case("type")
			read(key%value,*) ObsSpec(i)%type
		case("file")
			ObsSpec(i)%file=key%value
		case("beta","weight")
			read(key%value,*) ObsSpec(i)%beta
		case("scale")
			read(key%value,*) ObsSpec(i)%scale
		case("i2d")
			read(key%value,*) ObsSpec(i)%i2d
			if(ObsSpec(i)%i2d.gt.n2d) n2d=ObsSpec(i)%i2d
		case default
			call output("Keyword not recognised: " // trim(key%key2))
	end select
	
	return
	end

	subroutine ReadInstrument(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	i=key%nr1
	
	select case(key%key2)
		case("name")
			instrument(i)=key%value
		case("ntrans")
			read(key%value,*) instr_ntrans(i)
		case default
			call output("Keyword not recognised: " // trim(key%key2))
	end select
	
	return
	end

	subroutine ReadCIA(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	i=key%nr1
	
	if(key%key2.eq.' ') then
		read(key%value,*) do_cia
	else
		select case(key%key2)
			case("file")
				CIA(i)%filename=key%value
			case default
				call output("Keyword not recognised: " // trim(key%key2))
		end select
	endif
		
	return
	end

		
	subroutine GetOutputDir
	use GlobalSetup
	IMPLICIT NONE
	integer ncla
	character*500 readline,command

	outputdir='./outputARCiS/'

	ncla=2
1	continue
	call getarg(ncla,readline)
	if(readline(1:2).eq.'-o') then
		call getarg(1+ncla,outputdir)
		if(outputdir(len_trim(outputdir)-1:len_trim(outputdir)).ne.'/') then
			outputdir=trim(outputdir) // '/'
		endif
		ncla=ncla+1
		goto 1
	endif
	ncla=ncla+1
	if(readline.ne.' ') goto 1

	write(command,'("mkdir -p ",a)') trim(outputdir)
	call system(command)
	
	return
	end


	subroutine InitFreq()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 lam0,T0,Planck,tot
	integer i,j
	
	lam0=lam1
	nlam=1
	do while(lam0.le.lam2)
		lam0=lam0+lam0/specres
		nlam=nlam+1
	enddo
	
	allocate(lam(nlam))
	allocate(freq(nlam))
	allocate(dfreq(nlam))
	
	i=1
	lam(i)=lam1
	do while(lam(i).le.lam2)
		i=i+1
		lam(i)=lam(i-1)+lam(i-1)/specres
	enddo
	lam(nlam)=lam2
	
	do i=1,nlam
		freq(i)=1d0/lam(i)
	enddo

	do i=1,nlam-1
		dfreq(i)=abs(freq(i+1)-freq(i))
	enddo

	allocate(BB(nBB,nlam))

	do j=1,nBB
		T0=real(j)
		tot=0d0
		do i=1,nlam-1
			BB(j,i)=Planck(T0,freq(i))
			tot=tot+dfreq(i)*BB(j,i)
		enddo
		BB(j,i)=BB(j,i)*((2d0*(pi*kb*T0)**4)/(15d0*hplanck**3*clight**3))/tot
	enddo


	lam0=lam1
	nlamdust=1
	do while(lam0.le.lam2)
		lam0=lam0+lam0/specresdust
		nlamdust=nlamdust+1
	enddo
	
	allocate(lamdust(nlamdust))
	
	i=1
	lamdust(i)=lam1
	do while(lamdust(i).le.lam2)
		i=i+1
		lamdust(i)=lamdust(i-1)+lamdust(i-1)/specresdust
	enddo
	lamdust(nlamdust)=lam2
	

	return
	end


	subroutine ReadCloud(key)
	use GlobalSetup
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer j,j1,j2

	Cloud(key%nr1)%ptype="COMPUTE"

	if(key%nr1.eq.0) then
		j1=1
		j2=nclouds
	else
		j1=key%nr1
		j2=key%nr1
	endif

	do j=j1,j2
	select case(key%key2)
		case("file")
			Cloud(j)%file=trim(key%value)
		case("ngrains","nsize","nr")
			read(key%value,*) Cloud(j)%nr
		case("nsubgrains")
			read(key%value,*) Cloud(j)%nsubgrains
		case("amin")
			read(key%value,*) Cloud(j)%amin
		case("amax")
			read(key%value,*) Cloud(j)%amax
		case("fmax")
			read(key%value,*) Cloud(j)%fmax
		case("blend")
			read(key%value,*) Cloud(j)%blend
		case("porosity")
			read(key%value,*) Cloud(j)%porosity
		case("standard")
			Cloud(j)%standard=trim(key%value)
		case("species")
			Cloud(j)%species=trim(key%value)
		case("fcarbon")
			read(key%value,*) Cloud(j)%fcarbon
		case("pressure","p")
			read(key%value,*) Cloud(j)%P
		case("dp")
			read(key%value,*) Cloud(j)%dP
		case("s","sharpness")
			read(key%value,*) Cloud(j)%s
		case("column")
			read(key%value,*) Cloud(j)%column
		case("tau")
			read(key%value,*) Cloud(j)%tau
		case("lam")
			read(key%value,*) Cloud(j)%lam
		case("reff")
			read(key%value,*) Cloud(j)%reff
		case("veff")
			read(key%value,*) Cloud(j)%veff
		case("coverage")
			read(key%value,*) Cloud(j)%coverage
		case("frain")
			read(key%value,*) Cloud(j)%frain
		case("haze")
			read(key%value,*) Cloud(j)%haze
		case("mixrat")
			read(key%value,*) Cloud(j)%mixrat
		case("fcond")
			read(key%value,*) Cloud(j)%fcond
		case("tmix")
			read(key%value,*) Cloud(j)%tmix
		case("betamix")
			read(key%value,*) Cloud(j)%betamix
		case("kzzfile")
			Cloud(j)%Kzzfile=trim(key%value)
		case("kzz","k")
			read(key%value,*) Cloud(j)%Kzz
		case("kzz_pow","k_pow")
			read(key%value,*) Cloud(j)%Kzz_pow
		case("kscale")
			read(key%value,*) Cloud(j)%Kscale
		case("sigmadot","nucleation")
			read(key%value,*) Cloud(j)%Sigmadot
		case default
			call output("Unknown cloud keyword: " // trim(key%key2))
			stop
	end select
	enddo

	return
	end
	


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine SetupPartCloud(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,is,ilam,j
	real*8 phi,thet,tot,tot2,fact
	
	if(.not.allocated(Cloud(ii)%rv)) allocate(Cloud(ii)%rv(Cloud(ii)%nr))
	if(.not.allocated(Cloud(ii)%M)) allocate(Cloud(ii)%M(Cloud(ii)%nr))
	if(.not.allocated(Cloud(ii)%Kabs)) then
		allocate(Cloud(ii)%Kabs(Cloud(ii)%nr,nlam))
		allocate(Cloud(ii)%Ksca(Cloud(ii)%nr,nlam))
		allocate(Cloud(ii)%Kext(Cloud(ii)%nr,nlam))
		allocate(Cloud(ii)%F(Cloud(ii)%nr,nlam))
	endif

	select case(Cloud(ii)%ptype)
		case("COMPUTE")
			do is=1,Cloud(ii)%nr
				call tellertje(is,Cloud(ii)%nr)
				call ComputePart(Cloud(ii),ii,is)
			enddo
c		case("PARTFILE")
c			call ReadParticle(Cloud(ii),ii)
		case default
			call output("I did not understand what I was trying to do. Sorry!")
	end select
	
	if(nspike.gt.0.and.nspike.le.180) call output("making the first " // trim(int2string(nspike,'(i3)')) // " degrees isotropic")
	do ilam=1,nlam
		do is=1,Cloud(ii)%nr
			tot=0d0
			tot2=0d0
			do j=1,180
				tot=tot+Cloud(ii)%F(is,ilam)%F11(j)*sin(pi*(real(j)-0.5)/180d0)
				tot2=tot2+sin(pi*(real(j)-0.5)/180d0)
			enddo
			do j=1,180
				Cloud(ii)%F(is,ilam)%F11(j)=tot2*Cloud(ii)%F(is,ilam)%F11(j)/tot
				Cloud(ii)%F(is,ilam)%F12(j)=tot2*Cloud(ii)%F(is,ilam)%F12(j)/tot
				Cloud(ii)%F(is,ilam)%F22(j)=tot2*Cloud(ii)%F(is,ilam)%F22(j)/tot
				Cloud(ii)%F(is,ilam)%F33(j)=tot2*Cloud(ii)%F(is,ilam)%F33(j)/tot
				Cloud(ii)%F(is,ilam)%F34(j)=tot2*Cloud(ii)%F(is,ilam)%F34(j)/tot
				Cloud(ii)%F(is,ilam)%F44(j)=tot2*Cloud(ii)%F(is,ilam)%F44(j)/tot
			enddo

			if(nspike.gt.0.and.nspike.le.180) then
c the nspike parameter removes the n degree spike in the forward direction.
				do j=1,nspike
					fact=Cloud(ii)%F(is,ilam)%F11(nspike+1)/Cloud(ii)%F(is,ilam)%F11(j)
					Cloud(ii)%F(is,ilam)%F12(j)=Cloud(ii)%F(is,ilam)%F12(j)*fact
					Cloud(ii)%F(is,ilam)%F22(j)=Cloud(ii)%F(is,ilam)%F22(j)*fact
					Cloud(ii)%F(is,ilam)%F33(j)=Cloud(ii)%F(is,ilam)%F33(j)*fact
					Cloud(ii)%F(is,ilam)%F34(j)=Cloud(ii)%F(is,ilam)%F34(j)*fact
					Cloud(ii)%F(is,ilam)%F44(j)=Cloud(ii)%F(is,ilam)%F44(j)*fact
					Cloud(ii)%F(is,ilam)%F11(j)=Cloud(ii)%F(is,ilam)%F11(nspike+1)
				enddo

				tot=0d0
				tot2=0d0
				do j=1,180
					tot=tot+Cloud(ii)%F(is,ilam)%F11(j)*sin(pi*(real(j)-0.5)/180d0)
					tot2=tot2+sin(pi*(real(j)-0.5)/180d0)
				enddo
				Cloud(ii)%Ksca(is,ilam)=Cloud(ii)%Ksca(is,ilam)*tot/tot2
				Cloud(ii)%Kext(is,ilam)=Cloud(ii)%Kabs(is,ilam)+Cloud(ii)%Ksca(is,ilam)
				do j=1,180
					Cloud(ii)%F(is,ilam)%F11(j)=tot2*Cloud(ii)%F(is,ilam)%F11(j)/tot
					Cloud(ii)%F(is,ilam)%F12(j)=tot2*Cloud(ii)%F(is,ilam)%F12(j)/tot
					Cloud(ii)%F(is,ilam)%F22(j)=tot2*Cloud(ii)%F(is,ilam)%F22(j)/tot
					Cloud(ii)%F(is,ilam)%F33(j)=tot2*Cloud(ii)%F(is,ilam)%F33(j)/tot
					Cloud(ii)%F(is,ilam)%F34(j)=tot2*Cloud(ii)%F(is,ilam)%F34(j)/tot
					Cloud(ii)%F(is,ilam)%F44(j)=tot2*Cloud(ii)%F(is,ilam)%F44(j)/tot
				enddo
			endif
		enddo
	enddo
	do ilam=1,nlam
		do is=1,Cloud(ii)%nr
			Cloud(ii)%F(is,ilam)%IF11=0d0
			j=1
			thet=pi*(real(j)-0.5d0)/180d0
			Cloud(ii)%F(is,ilam)%IF11(j)=sin(thet)*Cloud(ii)%F(is,ilam)%F11(j)
			do j=2,180
				thet=pi*(real(j)-0.5d0)/180d0
				Cloud(ii)%F(is,ilam)%IF11(j)=Cloud(ii)%F(is,ilam)%IF11(j-1)+sin(thet)
     &			*Cloud(ii)%F(is,ilam)%F11(j)
			enddo
		enddo
	enddo
	
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------



	subroutine NameCloudSpecies(standard,species)
	IMPLICIT NONE
	character*20 standard
	character*500 species
	
	select case(standard)
		case("ENSTATITE","DIANA")
			species='MgSiO3(c)'
		case("ASTROSIL","FORSTERITE")
			species='Mg2SiO4(c)'
		case("VO")
			species='VO(c)'
		case("SiC")
			species='SiC(c)'
		case("Al2O3","CORRUNDUM")
			species='Al2O3(c)'
		case("Fe","IRON")
			species='Fe(c)'
		case("WATER")
			species='H2O(L)'
		case("BROOKITE","TiO2")	
c not entirely correct...
			species='TiO(c)'
		case("ICE")
			species='H2O(c)'
	end select

	return
	end
	
	
	
	subroutine ReadPlanetNameCSV()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,dR1,dR2,dM1,dM2
	character*100 name,namestar
	integer i,n
	character*10 Zc
	character*1000 line,key(100),value(100)

	open(unit=72,file=planetparameterfile,RECL=6000)
	read(72,'(a1000)') line
	call getkeys(line,key,n)

1	read(72,'(a1000)',end=2) line
	call getkeys(line,value,n)
	do i=1,n
		if(value(i).ne.' ') then
			select case(key(i))
				case('Planet_Name')
					name=trim(value(i))
				case('Star_Temperature_[K]')
					read(value(i),*) Tstar
				case('Star_Metallicity')
					read(value(i),*) Zc
				case('Star_Mass_[Ms]')
					read(value(i),*) Mstar
				case('Star_Radius_[Rs]')
					read(value(i),*) Rstar
				case('Star_Distance_[pc]')
					read(value(i),*) Distance
				case('Planet_Semi-major_Axis_[m]')
					read(value(i),*) Dplanet
					Dplanet=Dplanet*1d2/AU
				case('Planet_Mass_[Me]')
					read(value(i),*) Mplanet
					Mplanet=Mplanet*Mearth/MJup
				case('Planet_Radius_[Re]')
					read(value(i),*) Rplanet
					Rplanet=Rplanet*Rearth/RJup
			end select
		endif
		logg=log10(Ggrav*(Mstar*Msun)/((Rstar*Rsun)**2))
	enddo
	if(trim(name).eq.trim(planetname)) then
		close(unit=72)
		call output("Stellar mass:     " // dbl2string(Mstar,'(f7.2)') // "Msun")
		call output("Stellar T:        " // dbl2string(Tstar,'(f7.2)') // "K")
		call output("Stellar radius:   " // dbl2string(Rstar,'(f7.2)') // "Rsun")
		call output("Stellar logg:     " // dbl2string(logg,'(f7.2)'))
		call output("Stellar distance: " // dbl2string(Distance,'(f7.2)') // "pc")
		call output("Metallicity:      " // dbl2string(metallicity,'(f7.2)'))
		call output("Planet orbit:     " // dbl2string(Dplanet,'(f7.4)') // "AU")
		call output("Planet radius:    " // dbl2string(Rplanet,'(f7.2)') // "Rjup")
		call output("Planet mass:      " // dbl2string(Mplanet,'(f7.2)') // "Mjup")
		call output("Blackbody T:      " // dbl2string(sqrt(Rstar*Rsun/(2d0*Dplanet*AU))*Tstar,'(f7.2)') // "K")
		dR1=Rplanet*0.1
		dR2=Rplanet*0.1
		dM1=Mplanet*0.1
		dM2=Mplanet*0.2
		do i=1,n_ret
			select case(RetPar(i)%keyword)
				case("Rp","rp","RP")
					RetPar(i)%x0=Rplanet
					RetPar(i)%xmin=max(0d0,Rplanet-5d0*dR1)
					RetPar(i)%xmax=Rplanet+5d0*dR2
					if(RetPar(i)%xmin*Rjup.lt.0.1d0*Rearth) RetPar(i)%xmin=0.1d0*Rearth/Rjup
		call output("Minimum radius: " // dbl2string(RetPar(i)%xmin,'(f7.2)') // "Rjup")
		call output("Maximum radius: " // dbl2string(RetPar(i)%xmax,'(f7.2)') // "Rjup")
				case("Mp","mp","MP")
					RetPar(i)%x0=Mplanet
					RetPar(i)%xmin=max(0d0,Mplanet-5d0*dM1)
					RetPar(i)%xmax=Mplanet+5d0*dM2
					if(RetPar(i)%xmin*Mjup.lt.0.1d0*Mearth) RetPar(i)%xmin=0.1d0*Mearth/Mjup
		call output("Minimum mass:   " // dbl2string(RetPar(i)%xmin,'(f7.2)') // "Mjup")
		call output("Maximum mass:   " // dbl2string(RetPar(i)%xmax,'(f7.2)') // "Mjup")
				case("loggP","loggp")
					RetPar(i)%x0=log10(Ggrav*(Mplanet*Mjup)/((Rplanet*Rjup)**2))
					RetPar(i)%xmin=max(0.1,log10(Ggrav*(max(Mjup*(Mplanet-5d0*dM1),0.1d0*Mearth))/(((Rplanet+5d0*dR2)*Rjup)**2)))
					RetPar(i)%xmax=log10(Ggrav*((Mplanet+5d0*dM2)*Mjup)/((max((Rplanet-5d0*dR1)*Rjup,0.1d0*Rearth))**2))
		call output("Minimum logg:   " // dbl2string(RetPar(i)%xmin,'(f7.2)'))
		call output("Maximum logg:   " // dbl2string(RetPar(i)%xmax,'(f7.2)'))
			end select
		enddo
		return
	endif
	goto 1
2	close(unit=72)
	call output("Planet not found: " // trim(planetname))
	stop
	
	return
	end
	
	
	
	subroutine getkeys(line,key,n)
	character*1000 line
	character*1000 key(100)
	integer icomma,i,l,n,j

	icomma=1
	l=0
	i=1
	do while(icomma.ne.0)
		icomma=index(line(l+1:1000),',')
		if(icomma.ne.0) then
			key(i)=line(l+1:l+icomma-1)
			do j=1,icomma-1
				if(key(i)(j:j).eq.' ') key(i)(j:j)='_'
			enddo
			i=i+1
			l=l+icomma
		endif
	enddo
	n=i-1
	
	return
	end

	
	
	subroutine ReadPlanetName()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,dR1,dR2,dM1,dM2
	character*100 name,namestar
	integer i
	character*10 Zc

	print*,trim(planetname)

	i=len_trim(planetparameterfile)
	if(planetparameterfile(i-3:i).eq.'.csv') then
		call ReadPlanetNameCSV()
		return
	endif

	i=len_trim(planetname)
	if(planetname(i:i).eq.'b') then
		planetname(i:i)=' '
	endif
	
	open(unit=72,file=planetparameterfile,RECL=6000)
	read(72,*)
1	read(72,*,end=2) name,Tstar,x,x,Zc,x,x,Mstar,x,x,Rstar,x,x,logg,x,x,x,x,x,x,x,x,x,Dplanet,x,x,Mplanet,dM1,dM2,Rplanet,dR1,dR2
	if(Zc.eq.'-1') then
		metallicity=0d0
	else
		read(Zc,*) metallicity
	endif
	i=len_trim(name)
	if(name(i:i).eq.'b') then
		name(i:i)=' '
	endif
	if(trim(planetname).eq.trim(name)) then	
		close(unit=72)
		call output("Stellar mass:   " // dbl2string(Mstar,'(f7.2)') // "Msun")
		call output("Stellar T:      " // dbl2string(Tstar,'(f7.2)') // "K")
		call output("Stellar radius: " // dbl2string(Rstar,'(f7.2)') // "Rsun")
		call output("Stellar logg:   " // dbl2string(logg,'(f7.2)'))
		call output("Metallicity:    " // dbl2string(metallicity,'(f7.2)'))
		call output("Planet orbit:   " // dbl2string(Dplanet,'(f7.4)') // "AU")
		call output("Planet radius:  " // dbl2string(Rplanet,'(f7.2)') // "Rjup")
		call output("Planet mass:    " // dbl2string(Mplanet,'(f7.2)') // "Mjup")
		call output("Blackbody T:    " // dbl2string(sqrt(Rstar*Rsun/(2d0*Dplanet*AU))*Tstar,'(f7.2)') // "K")
		do i=1,n_ret
			select case(RetPar(i)%keyword)
				case("Rp","rp","RP")
					RetPar(i)%x0=Rplanet
					RetPar(i)%xmin=max(0d0,Rplanet-5d0*dR1)
					RetPar(i)%xmax=Rplanet+5d0*dR2
					if(RetPar(i)%xmin*Rjup.lt.0.1d0*Rearth) RetPar(i)%xmin=0.1d0*Rearth/Rjup
		call output("Minimum radius: " // dbl2string(RetPar(i)%xmin,'(f7.2)') // "Rjup")
		call output("Maximum radius: " // dbl2string(RetPar(i)%xmax,'(f7.2)') // "Rjup")
				case("Mp","mp","MP")
					RetPar(i)%x0=Mplanet
					RetPar(i)%xmin=max(0d0,Mplanet-5d0*dM1)
					RetPar(i)%xmax=Mplanet+5d0*dM2
					if(RetPar(i)%xmin*Mjup.lt.0.1d0*Mearth) RetPar(i)%xmin=0.1d0*Mearth/Mjup
		call output("Minimum mass:   " // dbl2string(RetPar(i)%xmin,'(f7.2)') // "Mjup")
		call output("Maximum mass:   " // dbl2string(RetPar(i)%xmax,'(f7.2)') // "Mjup")
				case("loggP","loggp")
					RetPar(i)%x0=log10(Ggrav*(Mplanet*Mjup)/((Rplanet*Rjup)**2))
					RetPar(i)%xmin=max(0.1,log10(Ggrav*(max(Mjup*(Mplanet-5d0*dM1),0.1d0*Mearth))/(((Rplanet+5d0*dR2)*Rjup)**2)))
					RetPar(i)%xmax=log10(Ggrav*((Mplanet+5d0*dM2)*Mjup)/((max((Rplanet-5d0*dR1)*Rjup,0.1d0*Rearth))**2))
		call output("Minimum logg:   " // dbl2string(RetPar(i)%xmin,'(f7.2)'))
		call output("Maximum logg:   " // dbl2string(RetPar(i)%xmax,'(f7.2)'))
			end select
		enddo
		return
	endif
	goto 1
2	close(unit=72)
	call output("Planet not found: " // trim(planetname))
	stop
	
	return
	end
	