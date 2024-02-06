c==============================================================================	
c Module for reading keywords
c==============================================================================	
	module ReadKeywords
	IMPLICIT NONE

	type SettingKey
		character*500 key1,key2,value,key
		character*500 orkey1,orkey2
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
	logical readfile,exist

	if(ncommandargs.gt.0) then
		inputfile=commandargs(1)
	else
		inputfile='input.dat'
	endif

	inquire(file=inputfile,exist=exist)
	if(.not.exist) then
		call output("Input file does not exist!")
		call output("filename: " // trim(inputfile))
		stop
	endif

	open(unit=20,file=inputfile,FORM="FORMATTED")

	call output("input file: " // trim(inputfile))

	call system("cp " // trim(inputfile) // " " // trim(outputdir) // "input.dat")
	open(unit=21,file=trim(outputdir) // "input.dat",FORM="FORMATTED",ACCESS="APPEND")
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
		readline=commandargs(2+ncla)
		if(readline(1:2).eq.'-s') then
			ncla=ncla+1
			readline=commandargs(2+ncla)
			call output("Command line argument: " // trim(readline))
			write(21,'(a)') trim(readline)
			ncla=ncla+1
		else if(readline(1:2).eq.'-o') then
			ncla=ncla+2
			goto 10
		else
			if(readline(1:1).ne.' ') then
c				try to read another command line argument
				close(unit=21)
				call system("cat " // trim(readline) // " >> " // trim(outputdir) // "input.dat")
				open(unit=21,file=trim(outputdir) // "input.dat",FORM="FORMATTED",ACCESS="APPEND")
				open(unit=20,file=readline,FORM="FORMATTED")
				readfile=.true.
				ncla=ncla+1
				goto 10
			else
c				all arguments are read
				goto 30
			endif
		endif
	endif

	if(readline(1:1).eq.' ') goto 10

	allocate(key%next)
	key => key%next
	key%last=.false.
	call get_key_value(readline,key%key,key%key1,key%key2,key%orkey1,key%orkey2,key%value,key%nr1,key%nr2,key%key2d)

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
	subroutine get_key_value(line_in,key,key1,key2,orkey1,orkey2,value,nr1,nr2,key2d)
	use GlobalSetup
	IMPLICIT NONE
	character*1000 line_in,line
	character*500 key,key1,key2,value,orkey1,orkey2
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
		orkey1=key1
		orkey2=key2
		call checknr(key1,nr1)
		call checknr(key2,nr2)
	else
		key1=line(1:ikey1-1)
		key=key1
		key2=' '
		orkey1=key1
		orkey2=key2
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
	use AtomsModule
	IMPLICIT NONE
	type(SettingKey),target :: firstkey
	type(SettingKey),pointer :: key
	integer i,j,ncia0,n
	character*500 homedir,h2h2file,h2hefile,h2ch4file
	integer ndiseq_list
	parameter(ndiseq_list=20)
	character*10 names(nmol_data),diseq_list(20)
	parameter(diseq_list = (/ "CH4       ","CO        ","CO2       ","H2        ","H2O       ",
     &						  "NH3       ","N2        ","C         ","CH2OH     ","CH3       ",
     &						  "CH3OH     ","C2H2      ","H         ","O         ","OH        ",
     &						  "N         ","NH        ","NH2       ","NO        ","N2H3      " /))
	logical existh2h2,existh2he,existh2ch4


	nmol=1
	nobs=0
	ncia=0
	nclouds=0
	n_ret=0
	n_Par3D=0
	n_instr=0
	nphase=0
	nVpoints=0
	nIRpoints=0
	nPhotoReacts=0
	nmodel_err=1
	j=0
	mixratfile=.false.
	Pmin=1d-6
	Pmax=1d+3
	freePT_fitT=.false.
	freePT_fitP=.true.
	Rp_range=20d0

	i2d=0

	nr=20
	nrsurf=0
	key => firstkey
	do while(.not.key%last)
		select case(key%key1)
			case("nr")
				read(key%value,*) nr
			case("nrsurf")
				read(key%value,*) nrsurf
		end select
		key=>key%next
	enddo
	nr=nr+nrsurf

	nTpoints=0
	key => firstkey
	do while(.not.key%last)
		select case(key%key1)
			case("ntpoints","nfreet")
				read(key%value,*) nTpoints
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
			case("phase")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nphase) nphase=key%nr1
			case("tauvpoint")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nVpoints) nVpoints=key%nr1
			case("tauirpoint")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nIRpoints) nIRpoints=key%nr1
			case("par3d")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.n_Par3D) n_Par3D=key%nr1
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
			case("diseq")
				read(key%value,*) disequilibrium
			case("rp_range")
				read(key%value,*) Rp_range
			case("tpfile")
				read(key%value,'(a)') TPfile
			case("cloud")
c				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nclouds) nclouds=key%nr1
			case("retpar","fitpar")
				if(key%key2.eq.'keyword') then
					if(key%value.eq.'tprofile') then
						free_tprofile=.true.
						n_ret=n_ret+nTpoints
						if(freePT_fitP) n_ret=n_ret+nTpoints-2
					else
						n_ret=n_ret+1
					endif
				endif
			case("free_fitt","freetp_fitt","freept_fitt")
				read(key%value,*) freePT_fitT
			case("free_fitp","freetp_fitp","freept_fitp")
				read(key%value,*) freePT_fitP
			case("free_fitdt","freetp_fitdt","freept_fitdt")
				read(key%value,*) freePT_fitT
				freePT_fitT=.not.freePT_fitT
			case("pmin")
				read(key%value,*) pmin
			case("pmax")
				read(key%value,*) pmax
			case("model_err_abs","model_err_rel")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nmodel_err) nmodel_err=key%nr1
			case("photoreac","photoreactant","photoprod","photoproduct","photoeff")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr1.gt.nPhotoReacts) nPhotoReacts=key%nr1
			case default
				do i=1,nmol_data
					if(key%key.eq.molname(i)) then
						if(i.gt.nmol) nmol=i
						j=j+1
					endif
				enddo
		end select
		key=>key%next
	enddo

	if(mixratfile) then
		open(unit=30,file=TPfile,FORM="FORMATTED")
		read(30,*) n
		read(30,*) names(1:n)
		close(unit=30)
		do j=1,n
			do i=1,nmol_data
				if(names(j).eq.molname(i)) then
					if(i.gt.nmol) nmol=i
				endif
			enddo
		enddo
	endif

	if(disequilibrium) then
c select at least the species relevant for disequilibrium chemistry
		do j=1,ndiseq_list
			do i=1,nmol_data
				if(diseq_list(j).eq.molname(i)) then
					if(i.gt.nmol) nmol=i
				endif
			enddo
		enddo
	endif
	
	if(do_cia) nmol=max(nmol,48)

	allocate(mixrat(nmol))
	allocate(Pswitch_mol(nmol),abun_switch_mol(nmol))
	allocate(includemol(nmol))
	allocate(diseqmol(nmol))
	allocate(opacitymol(nmol))
	allocate(Cloud(max(nclouds,1)))
	do i=1,nclouds
		allocate(Cloud(i)%abun(40))
		allocate(Cloud(i)%nax(40))
		allocate(Cloud(i)%rho_mat(40))
		allocate(Cloud(i)%lnkfile(40,3))
		allocate(Cloud(i)%material(40))
	enddo
	allocate(XeqCloud(nr,max(nclouds,1)))
	allocate(XeqCloud_old(nr,max(nclouds,1)))
	allocate(XCloud(nr,max(nclouds,1)))
	allocate(RetPar(max(n_ret,1)))
	allocate(tau_Vpoint(max(nVpoints,1)),dT_Vpoint(max(nVpoints,1)))
	allocate(tau_IRpoint(max(nIRpoints,1)),dT_IRpoint(max(nIRpoints,1)))
	allocate(dTpoint(max(nTpoints,1)))
	allocate(Ppoint(max(nTpoints,1)))
	allocate(ObsSpec(max(nobs,1)))
	allocate(Tin(nr))
	allocate(instrument(max(n_instr,1)))
	allocate(instr_ntrans(max(n_instr,1)))
	allocate(instr_nobs(max(n_instr,1)))
	allocate(Par3D(max(n_Par3D,1)))
	allocate(theta_phase(max(nphase,1)))
	allocate(model_err_abs(max(nmodel_err,1)))
	allocate(model_err_rel(max(nmodel_err,1)))
	allocate(model_err_lam(max(nmodel_err,1)))
	allocate(PhotoReacts(max(nPhotoReacts,1)))
	do i=1,nPhotoReacts
		allocate(PhotoReacts(i)%react(max(nmol_data,N_atoms)))
		allocate(PhotoReacts(i)%product(nmol_data))
		allocate(PhotoReacts(i)%abun(nr,nmol_data))
	enddo

	do i=1,nTpoints
		Ppoint(i)=exp(log(Pmin)+log(Pmax/Pmin)*real(i-1)/real(nTpoints-1))
	enddo
	PrefTpoint=Pmin

	ncia0=0
	existh2h2=.false.
	existh2he=.false.
	existh2ch4=.false.
	if(do_cia) then
		call getenv('HOME',homedir)
c find H2-H2 cia file
		h2h2file=trim(homedir) // '/HITRAN/H2-H2_combined.cia'
		inquire(file=h2h2file,exist=existh2h2)
		if(existh2h2) then
			ncia0=ncia0+1
		else
			h2h2file=trim(homedir) // '/HITRAN/CIA/H2-H2_combined.cia'
			inquire(file=h2h2file,exist=existh2h2)
			if(existh2h2) then
				ncia0=ncia0+1
			else
				h2h2file=trim(homedir) // '/ARCiS/Data/CIA/H2-H2_combined.cia'
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
		if(ncia0.eq.0) then
			call output("NO CIA FILES FOUND: rerun with cia=.false.")
			stop
		endif
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


	subroutine getcommandline()
	use GlobalSetup
	IMPLICIT NONE
	integer i
	
	ncommandargs=iargc()
	allocate(commandargs(ncommandargs+1))
	do i=1,ncommandargs
		call getarg(i,commandargs(i))
	enddo
	commandargs(ncommandargs+1)=' '
	
	return
	end


	subroutine Init()
	use GlobalSetup
	use Constants
	use ReadKeywords
	use Struct3D
	use TimingModule
	use ARCiS_GGCHEM
	IMPLICIT NONE
	type(SettingKey),pointer :: key,first
	type(SettingKey) keyret
	integer i,j,omp_get_max_threads,omp_get_thread_num
	real*8 tot,tot2,theta,Planck
	real*8,allocatable :: var(:),dvar(:)
	character*1000 line
	character*500 file,homedir
	integer ndiseq_list
	parameter(ndiseq_list=20)
	character*10 diseq_list(20)
	parameter(diseq_list = (/ "CH4       ","CO        ","CO2       ","H2        ","H2O       ",
     &						  "NH3       ","N2        ","C         ","CH2OH     ","CH3       ",
     &						  "CH3OH     ","C2H2      ","H         ","O         ","OH        ",
     &						  "N         ","NH        ","NH2       ","NO        ","N2H3      " /))

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

	timechem=0d0
	timecloud=0d0
	timetemp=0d0
	itimechem=0
	itimecloud=0
	itimetemp=0
	ctimechem=0
	ctimecloud=0
	ctimetemp=0
	call system_clock(count_rate=rate)

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

	if(useobsgrid) then
		if(nobs.le.0) useobsgrid=.false.
	endif

	if(doscaleR.and.log_emis) then
		call output("ScaleR and log_emis settings both to true does not work toegther")
		call output("switching to log_emis=.false.")
		log_emis=.false.
	endif

	if(retrieval) then

	key => first%next

	do while(.not.key%last)

		do i=1,n_ret
			line=trim(RetPar(i)%keyword) // "=0d0"
			call get_key_value(line,keyret%key,keyret%key1,keyret%key2,keyret%orkey1,keyret%orkey2,
     &						keyret%value,keyret%nr1,keyret%nr2,keyret%key2d)
			if(trim(keyret%key1).eq.trim(key%key1).and.trim(keyret%key2).eq.trim(key%key2).and.
     &		   keyret%nr1.eq.key%nr1.and.keyret%nr2.eq.key%nr2.and.keyret%key2d.eq.key%key2d) then
				read(key%value,*) RetPar(i)%x0
			endif
		enddo

		key => key%next
	
	enddo
	
	endif

	if(orbit_P.lt.0d0) orbit_P=sqrt(Dplanet**3/Mstar)*365.25*86400d0

	allocate(tauUV(nr),kappaUV(nr))

	call ConvertUnits()

c	condensates=(condensates.or.cloudcompute)

	allocate(gg(ng),wgg(ng))
	call gauleg(0d0,1d0,gg,wgg,ng)
	
	call InitFreq()
	if(nphase.le.0) then
		nphase=1
		theta_phase(1)=180d0
	endif

	if(particledir.eq.' ') particledir=trim(outputdir)

	if(opacitydir(len_trim(opacitydir)-1:len_trim(opacitydir)).ne.'/') then
		opacitydir=trim(opacitydir) // '/'
	endif
	call Init_optEC
	if(opacitymode) then
		call InitOpacityMode()
	else
		allocate(dens(nr))
		allocate(Ndens(nr))
		allocate(R(nr+1))
		allocate(T(nr))
		allocate(P(nr+1))
		allocate(Hp(nr))
		allocate(nabla_ad(nr))
		allocate(grav(nr))
		allocate(MMW(nr))
		allocate(didcondens(nr))
		allocate(mixrat_r(nr,nmol))
		allocate(mixrat_old_r(nr,nmol))
		allocate(mixrat_optEC_r(nr))
		allocate(cloud_dens(nr,max(nclouds,1)))
		cloud_dens=0d0
		call InitDens()
		call InitObs()
		do i=1,nclouds
			call output("==================================================================")
			call output("Setting up cloud: " // trim(int2string(i,'(i4)')))
			allocate(Cloud(i)%frac(nr,40))
			do j=1,nr
				Cloud(i)%frac(j,1:40)=Cloud(i)%abun(1:40)
			enddo
			allocate(Cloud(i)%cryst(nr,40))
			Cloud(i)%cryst=Cloud(i)%cryst0
			if(Cloud(i)%type.eq.'DIFFUSE'.or.Cloud(i)%type.eq.'WATER') then
				Cloud(i)%opacitytype="MATERIAL"
				Cloud(i)%nmat=14
				Cloud(i)%material(1)='RUTILE'
				Cloud(i)%material(2)='FORSTERITE'
				Cloud(i)%material(3)='SiO'
				Cloud(i)%material(4)='SiO2'
				Cloud(i)%material(5)='IRON'
				Cloud(i)%material(6)='CORRUNDUM'
				Cloud(i)%material(7)='FeO'
				Cloud(i)%material(8)='MgO'
				Cloud(i)%material(9)='ENSTATITE'
				Cloud(i)%material(10)='CARBON'
				Cloud(i)%material(11)='SiC'
				Cloud(i)%material(12)='WATER'
				Cloud(i)%material(13)='PYROXENE'
				Cloud(i)%material(14)='A-SiO2'
				if(Cloud(i)%haze) then
					Cloud(i)%nmat=Cloud(i)%nmat+1
					Cloud(i)%material(15)=Cloud(i)%hazetype
				endif
			endif
		enddo
		call SetupMaterialCloud
	endif

	metallicity0=metallicity

c	allocate(Cabs_mol(nr,ng,nmol,nlam)) ! efficient, though unlogical storage
	allocate(Cabs_mol(ng,nlam,nmol,nr))
	allocate(Cext_cont(nr,nlam))
	allocate(Cabs(nr,nlam,ng))
	allocate(Csca(nr,nlam))
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(Cabs_mol,Cext_cont,Cabs,Csca,nlam,ng,nr,nmol)
!$OMP DO SCHEDULE(STATIC)
		do i=1,nlam
			Cabs_mol(1:ng,i,1:nmol,1:nr)=0d0
			Cext_cont(1:nr,i)=0d0
			Cabs(1:nr,i,1:ng)=0d0
			Csca(1:nr,i)=0d0
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	do i=1,360
		theta=(real(i)-0.5d0)*pi/180d0
		sintheta(i)=sin(theta)
		costheta(i)=cos(theta)
	enddo

	allocate(Fstar(nlam))
	call StarSpecSetup(Tstar,logg,1d4*lam,1d4*blam,Fstar,nlam,starfile,blackbodystar)
	Fstar=Fstar*pi*Rstar**2

	call output("==================================================================")

	diseqmol=.false.
	if(disequilibrium) then
c select at least the species relevant for disequilibrium chemistry
		do j=1,ndiseq_list
			do i=1,nmol_data
				if(diseq_list(j).eq.molname(i)) then
					diseqmol(i)=.true.
				endif
			enddo
		enddo
	endif
	do i=1,nmol
		opacitymol(i)=.false.
		if(includemol(i)) call InitReadOpacityFITS(i)
	enddo

	call ReadDataCIA()

	if(retrieval) then
		do i=1,n_ret
			if(RetPar(i)%x0.lt.-1d150) then
				if(RetPar(i)%logscale) then
					RetPar(i)%x0=sqrt(RetPar(i)%xmax*RetPar(i)%xmin)
				else
					RetPar(i)%x0=0.5d0*(RetPar(i)%xmax+RetPar(i)%xmin)
				endif
			endif
			if(RetPar(i)%dx.lt.0d0) then
				if(RetPar(i)%logscale) then
					RetPar(i)%dx=5d0*(RetPar(i)%xmax/RetPar(i)%xmin)
				else
					RetPar(i)%dx=5d0*(RetPar(i)%xmax-RetPar(i)%xmin)
				endif
			endif
			if(RetPar(i)%keyword.eq.'Tstar'.or.RetPar(i)%keyword.eq.'tstar'.or.RetPar(i)%keyword.eq.'Rstar'
     &	.or.RetPar(i)%keyword.eq.'rstar'.or.RetPar(i)%keyword.eq.'logg'.or.RetPar(i)%keyword.eq.'TStar'
     &	.or.RetPar(i)%keyword.eq.'TSTAR'.or.RetPar(i)%keyword.eq.'RStar'.or.RetPar(i)%keyword.eq.'RSTAR'
     &	.or.RetPar(i)%keyword.eq.'LOGG') retrievestar=.true.
		enddo
	endif

	dobackgroundgas=.false.
	do i=1,nmol
		if(backgroundgas(i)) dobackgroundgas=.true.
	enddo
	if(dochemistry.or.secondary_atmosphere) then
		outputdirGGchem=outputdir
		allocate(usemolGGchem(nmol))
		usemolGGchem(1:nmol)=includemol(1:nmol)
		if(secondary_atmosphere) then
			call init_GGchem(molname,nmol,.true.)
		else
			call init_GGchem(molname,nmol,condensates)
		endif
		dobackgroundgas=.false.
	endif

	if(pos_dT_lowest) then
		do i=1,n_ret
			if(RetPar(i)%keyword(1:7).eq."dTpoint") then
				read(RetPar(i)%keyword(8:len(RetPar(i)%keyword)),*) j
				if(j.eq.nTpoints) RetPar(i)%xmin=0d0
			endif
		enddo
	endif

	do i=1,n_ret
		if(RetPar(i)%keyword(1:6).eq."Tpoint") then
			RetPar(i)%xmin=Tmin
			RetPar(i)%xmax=Tmax
			RetPar(i)%logscale=logTprofile
		endif
	enddo

	if(pos_dT) then
		do i=1,n_ret
			if(RetPar(i)%keyword(1:7).eq."dTpoint") then
				RetPar(i)%xmin=0d0
			endif
			if(RetPar(i)%keyword(1:6).eq."Tpoint") then
				read(RetPar(i)%keyword(7:len_trim(RetPar(i)%keyword)),*) j
				if(j.gt.1) RetPar(i)%increase=.true.
			endif
		enddo
	endif

	if(iWolk.gt.0) then
		open(unit=50,file=trim(outputdir) // "/Wolk.dat",FORM="FORMATTED",ACCESS="STREAM")
		allocate(var(n_ret),dvar(n_ret))
		do i=1,iWolk
			read(50,*) j,tot,var(1:n_ret)
		enddo
		call MapRetrieval(var,dvar)
		do i=1,n_ret
			call output(RetPar(i)%keyword(1:15) // trim(dbl2string(RetPar(i)%value,'(es15.3)')))
		enddo
		close(unit=50)
	endif
	call output("==================================================================")
	
	allocate(lamemis(nlam),lamtrans(nlam))
	lamemis=emisspec
	lamtrans=transspec

	allocate(surface_emis(nlam))
	surface_emis=1d0
	select case(surfacetype)
		case("FILE","file")
			call regridSimple(surfacefile,lam*1d4,surface_emis,nlam)
			surface_emis=1d0-surface_emis/100d0
		case("Earth","EARTH","earth")
			allocate(surface_emis_ice(nlam))
			allocate(surface_emis_snow(nlam))
			allocate(surface_emis_grass(nlam))
			allocate(surface_emis_sand(nlam))
			allocate(surface_emis_water(nlam))
			call getenv('HOME',homedir)
			file=trim(homedir) // '/ARCiS/Data/Surface/Ice.dat'
			call regridSimple(file,lam*1d4,surface_emis_ice,nlam)
			surface_emis_ice=1d0-surface_emis_ice/100d0
			file=trim(homedir) // '/ARCiS/Data/Surface/Snow.dat'
			call regridSimple(file,lam*1d4,surface_emis_snow,nlam)
			surface_emis_snow=1d0-surface_emis_snow/100d0
			file=trim(homedir) // '/ARCiS/Data/Surface/Grass.dat'
			call regridSimple(file,lam*1d4,surface_emis_grass,nlam)
			surface_emis_grass=1d0-surface_emis_grass/100d0
			file=trim(homedir) // '/ARCiS/Data/Surface/brown-darkbrown-sand.dat'
			call regridSimple(file,lam*1d4,surface_emis_sand,nlam)
			surface_emis_sand=1d0-surface_emis_sand/100d0
			file=trim(homedir) // '/ARCiS/Data/Surface/Water.dat'
			call regridSimple(file,lam*1d4,surface_emis_water,nlam)
			surface_emis_water=1d0-surface_emis_water/100d0
	end select
	
	if(makemovie) makeimage=.true.

c If reading in a full 3D model (from e.g. a GCM model) the number of 3D models needs to be equal to nlatt*nlong
	if(readFull3D) then
		n3D=(nlong-1)*(nlatt-1)
		call output("Full 3D mode: Number of 1D models:  " // int2string(n3D,'(i4)'))
		call output("              Number of longitudes: " // int2string(nlong,'(i4)'))
		call output("              Number of latitudes:  " // int2string(nlatt,'(i4)'))
c In this case the beta map should be the static one. Make sure this is set properly.
		night2day=0d0
		pole2eq=1d0
		tidallock=.true.
		vxx=0d0
		fDay=1d0
		Kxx=1d0
		Kyy=1d0
		powvxx=0d0
		hotspotshift0=-1d5
	endif
	
	if(deepredist.and..not.do3D) deepredist=.false.
	if(deepredist) then
		do i=1,len_trim(deepredisttype)
			if(iachar(deepredisttype(i:i)).ge.65.and.iachar(deepredisttype(i:i)).le.90) then
				deepredisttype(i:i)=achar(iachar(deepredisttype(i:i))+32)
			endif
		enddo
	endif
	if(do3D) allocate(Tprev3D(nr))
		
	allocate(long(nlong),latt(nlatt))
	allocate(tanx(nlong),tany(nlong))
	allocate(cost2(nlatt),beta3D_eq(nlong),x3D_eq(nlong))
	allocate(ibeta(nlong,nlatt))
	allocate(beta3D(n3D),x3D(n3D))

	if(fulloutput3D) then
		allocate(PTaverage3D(0:nphase,nr))
		allocate(mixrat_average3D(0:nphase,nr,nmol))
	endif
	
	if(doRing) allocate(FRing(nlam))

	if(planetform) call InitFormation(Mstar,Tstar,Rstar,planetform_SolidC,planetform_Macc)
	
	if(useDLMie) call InitDLMie()
	
	do i=1,nmodel_err-1
		model_err_lam(i)=10d0**(log10(lam(1))+log10(lam(nlam)/lam(1))*(real(i)/real(nmodel_err)))
	enddo
	
	return
	end


	subroutine ReadAndSetKey(key)
	use GlobalSetup
	use ReadKeywords
	use CloudModule
	use Struct3D
	use ARCiS_GGCHEM
	IMPLICIT NONE
	type(SettingKey) key
	integer i

	if(key%key2d.eq.i2d.or.key%key2d.eq.0) then

	select case(key%key1)
		case("compute")
			call output("keyword 'compute' no longer supported")
		case("nr")
c is already set in CountStuff
c			read(key%value,*) nr
		case("nrsurf")
c is already set in CountStuff
c			read(key%value,*) nrsurf
		case("psurf")
			read(key%value,*) Psurf
		case("nrcloud","nrhr")
			read(key%value,*) nr_cloud
			nr_cloud=max(1,nr_cloud/nr)
		case("nrstepchem")
			read(key%value,*) nrstepchem
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
		case("constant_g")
			read(key%value,*) constant_g
		case("fixmmw")
			read(key%value,*) fixMMW
		case("mmw")
			read(key%value,*) MMW0
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
		case("bbstar")
			read(key%value,*) blackbodystar
		case("dp","dplanet")
			read(key%value,*) Dplanet
		case("retrieval")
			read(key%value,*) retrieval
		case("doscaler","scaler")
			read(key%value,*) doscaleR
		case("nscaler")
			read(key%value,*) nscaleR
		case("randomstart")
			read(key%value,*) randomstart
		case("contrib","computecontrib")
			read(key%value,*) computecontrib
		case("outputopacity","writeopacity")
			read(key%value,*) outputopacity
		case("phase")
			read(key%value,*) theta_phase(key%nr1)
		case("nphase")
			print*,'WARNING: no longer supported keyword nphase'
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
		case("tmin")
			read(key%value,*) Tmin
		case("tmax")
			read(key%value,*) Tmax
		case("logtprofile")
			read(key%value,*) logTprofile
		case("taurexprofile")
			read(key%value,*) taurexprofile
		case("taurexsmooth")
			read(key%value,*) taurexsmooth
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
		case("specreslr","specres_lr","specresrt","specres_rt")
			read(key%value,*) specres_LR
		case("specresdust")
			call output("Parameter specresdust is not in use anymore!")
		case("usexs")
			read(key%value,*) useXS
		case("specresfile")
			specresfile=trim(key%value)
		case("particledir","dirparticle")
			particledir=trim(key%value)
		case("rayleigh")
			read(key%value,*) do_rayleigh
		case("tpfile")
			TPfile=key%value
		case("mixratfile")
			read(key%value,*) mixratfile
		case("gridtpfile")
			read(key%value,*) gridTPfile
		case("maxt","maxtprofile")
			read(key%value,*) maxTprofile
		case("pswitch")
			call ReadPswitch(key)
		case("abunswitch","abun_switch")
			call ReadAbunSwitch(key)
		case("background","backgroundgas")
			call ReadBackgroundgas(key)
		case("isotope","f_isotope")
			call ReadIsotope(key)
		case("setsurfp","setsurfpressure")
			read(key%value,*) setsurfpressure
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
		case("cloud")
			call ReadCloud(key)
		case("complexkzz")
			read(key%value,*) complexKzz
		case("scattering")
			read(key%value,*) scattering
		case("scattstar","starscatt")
			read(key%value,*) scattstar
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
		case("exp_ad")
			read(key%value,*) exp_ad
		case("isofstar")
			read(key%value,*) isoFstar
		case("forceebalance")
			read(key%value,*) forceEbalance
		case("teffp","tplanet")
			read(key%value,*) TeffP
		case("maxiter")
			read(key%value,*) maxiter
		case("miniter")
			read(key%value,*) miniter
		case("epsiter")
			read(key%value,*) epsiter
		case("chemistry")
			read(key%value,*) dochemistry
		case("diseq")
			read(key%value,*) disequilibrium
		case("ggchem_piter")
			read(key%value,*) GGCHEM_P_iter
		case("ggchem_tfast")
			read(key%value,*) GGCHEM_Tfast
		case("kzz")
			read(key%value,*) Kzz
		case("kzz_deep")
			read(key%value,*) Kzz_deep
		case("kzz_1bar","kzz_up")
			read(key%value,*) Kzz_1bar
		case("kzz_p","kzz_pow")
			read(key%value,*) Kzz_P
		case("kzz_contrast")
			read(key%value,*) Kzz_contrast
		case("metallicity")
			read(key%value,*) metallicity
		case("coratio")
			read(key%value,*) COratio
		case("sioratio")
			read(key%value,*) SiOratio
		case("noratio")
			read(key%value,*) NOratio
		case("soratio")
			read(key%value,*) SOratio
		case("tiscale")
			read(key%value,*) TiScale
		case("inversecoratio")
			read(key%value,*) inverseCOratio
		case("elementfile")
			element_abun_file=key%value
		case("condensates")
			read(key%value,*) condensates
		case("elementlist")
			elements_ARCiS=key%value
		case("secondary_atmosphere")
			read(key%value,*) secondary_atmosphere
		case("waterworld")
			read(key%value,*) WaterWorld
		case("simplerainout")
			read(key%value,*) dosimplerainout
		case("mixp")
			read(key%value,*) mixP
		case("sinkz")
			read(key%value,*) sinkZ
		case("alphaz")
			read(key%value,*) alphaZ
		case("nspike")
			read(key%value,*) nspike
		case("dlmie")
			read(key%value,*) useDLMie
		case("nphot")
			read(key%value,*) Nphot0
		case("retpar","fitpar")
			call ReadRetrieval(key)
		case("obs")
			call ReadObsSpec(key)
		case("model_err_abs")
			read(key%value,*) model_err_abs(key%nr1)
		case("model_err_rel")
			read(key%value,*) model_err_rel(key%nr1)
		case("par3d")
			call ReadPar3D(key)
		case("useobsgrid")
			read(key%value,*) useobsgrid
		case("logemis","log_emis")
			read(key%value,*) log_emis
		case("massprior")
			read(key%value,*) massprior
		case("rp_range")
			read(key%value,*) Rp_range
		case("mp_prior")
			read(key%value,*) Mp_prior
		case("dmp_prior")
			read(key%value,*) dMp_prior
		case("nboot")
			read(key%value,*) nboot
		case("npew")
			read(key%value,*) npew
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
		case("consteff")
			read(key%value,*) const_eff_multinest
		case("fmultinest","efr")
			read(key%value,*) f_multinest
		case("tolmultinest","ftol")
			read(key%value,*) tol_multinest
		case("nestupdate")
			read(key%value,*) nest_update
		case("retrievaltype")
			read(key%value,*) retrievaltype
		case("writewolk")
			read(key%value,*) writeWolk
		case("tauvpoint")
			read(key%value,*) tau_Vpoint(key%nr1)
		case("tauirpoint")
			read(key%value,*) tau_IRpoint(key%nr1)
		case("dtvpoint")
			read(key%value,*) dT_Vpoint(key%nr1)
		case("dtirpoint")
			read(key%value,*) dT_IRpoint(key%nr1)
		case("free_tprofile")
			read(key%value,*) free_tprofile
		case("free_fitt","freetp_fitt","freept_fitt")
*			read(key%value,*) freePT_fitT
		case("free_fitp","freetp_fitp","freept_fitp")
*			read(key%value,*) freePT_fitP
		case("free_fitdt","freetp_fitdt","freept_fitdt")
*			read(key%value,*) freePT_fitT
*			freePT_fitT=.not.freePT_fitT
		case("wiggle_err")
			read(key%value,*) wiggle_err
		case("tpoint","dtpoint")
			read(key%value,*) dTpoint(key%nr1)
		case("ppoint")
			read(key%value,*) Ppoint(key%nr1)
		case("ntpoints","nfreet")
c is already set in CountStuff
c			read(key%value,*) nTpoints
		case("preftpoint","pref")
			read(key%value,*) PrefTpoint
		case("pos_dt_lowest")
			read(key%value,*) pos_dT_lowest
		case("pos_dt")
			read(key%value,*) pos_dT
		case("faircoverage")
			read(key%value,*) faircoverage
		case("speclimits")
			read(key%value,*) speclimits
		case("adiabatic","adiabatic_tprofile")
			read(key%value,*) adiabatic_tprofile
		case("outflow")
			read(key%value,*) outflow
		case("fday")
			read(key%value,*) fDay
		case("kxx")
			read(key%value,*) Kxx
		case("kyy")
			read(key%value,*) Kyy
		case("vxx")
			read(key%value,*) vxx
		case("powvxx","powv","vxx_pow")
			read(key%value,*) powvxx
		case("night2day")
			read(key%value,*) night2day
		case("fixnight2day","fixn2d")
			read(key%value,*) fixnight2day
		case("pole2eq")
			read(key%value,*) pole2eq
		case("tidallock")
			read(key%value,*) tidallock
		case("hotspotshift")
			read(key%value,*) hotspotshift0
		case("n3d")
			read(key%value,*) n3D
		case("par3dsteepness","steepness3d")
			read(key%value,*) par3Dsteepness
		case("nnu","nnustar")
			read(key%value,*) nnu0
		case("nlong")
			read(key%value,*) nlong
		case("nlatt")
			read(key%value,*) nlatt
		case("betapow")
			read(key%value,*) betapow
		case("makeai")
			read(key%value,*) domakeai
		case("parametergridfile")
			pargridfile=key%value
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
		case("vfrag")
			read(key%value,*) vfrag
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
		case("do3d","run3d")
			read(key%value,*) do3D
		case("output3d")
			read(key%value,*) fulloutput3D
		case("deepredist")
			read(key%value,*) deepredist
		case("deepredisttype")
			read(key%value,*) deepredisttype
		case("tsurface")
			read(key%value,*) Tsurface0
		case("readfull3d")
			read(key%value,*) readFull3D
		case("computealbedo","planetalbedo")
			read(key%value,*) computealbedo
		case("iwolk")
			read(key%value,*) iWolk
		case("emisspec")
			read(key%value,*) emisspec
		case("transspec")
			read(key%value,*) transspec
		case("dotranshide")
			read(key%value,*) dotranshide
		case("makeimage")
			read(key%value,*) makeimage
		case("makemovie")
			read(key%value,*) makemovie
		case("orbit")
			select case(key%key2)
				case("p")
					read(key%value,*) orbit_P
				case("e")
					read(key%value,*) orbit_e
				case("omega")
					read(key%value,*) orbit_omega
				case("inc")
					read(key%value,*) orbit_inc
			end select
		case("computelc")
			read(key%value,*) computeLC
		case("vrot")
			read(key%value,*) vrot0
		case("planetform")
			read(key%value,*) planetform
		case("fdust")
			read(key%value,*) planetform_fdust
		case("fplanet")
			read(key%value,*) planetform_fplan
		case("rstart","rcore")
			call output("Switching to depreciated mode of SimAb!!!")
			read(key%value,*) planetform_Dmigrate
			planetform_Dmigrate=planetform_Dmigrate-Dplanet
			planetform_Rend=Dplanet
			call output("Using dMigrate    = " // dbl2string(planetform_Dmigrate,"(f6.3)") // "AU")
			call output("Using RendMigrate = " // dbl2string(planetform_Rend,"(f6.3)") // "AU")
		case("dmigrate")
			read(key%value,*) planetform_Dmigrate
		case("rendmigrate")
			read(key%value,*) planetform_Rend
		case("mstart","mcore")
			read(key%value,*) planetform_Mstart
		case("sootline","solidc","fsolidc")
			read(key%value,*) planetform_SolidC
		case("mdotdisk")
			read(key%value,*) planetform_Macc
		case("surfacetype")
			read(key%value,*) surfacetype
		case("surfacefile")
			surfacefile=trim(key%value)
			call checkfile(surfacefile)
		case("surfacealbedo")
			read(key%value,*) surfacealbedo
		case("fwater","focean")
			read(key%value,*) f_surface_water
		case("ncpah","nc_pah")
			read(key%value,*) nC_PAH
		case("pah")
			read(key%value,*) mixrat_PAH
		case("optec")
			read(key%value,*) mixrat_optEC0
		case("rad_optec")
			read(key%value,*) rad_optEC
		case("eg_optec")
			read(key%value,*) Eg_optEC
		case("fixmol")
			call ReadFixMol(key)
		case("photoreac","photoreactant","photoprod","photoproduct")
			call ReadPhotoChem(key)
		case("photoeff")
			read(key%value,*) PhotoReacts(key%nr1)%f_eff
		case("photohaze")
			read(key%value,*) PhotoReacts(key%nr1)%haze
		case("kappauv")
			read(key%value,*) kappaUV0
		case("hydrogenloss")
			read(key%value,*) Hydrogenloss
		case("rring")
			read(key%value,*) Rring
		case("drring")
			read(key%value,*) dRring
		case("tauring")
			read(key%value,*) tauRing
		case("doring")
			read(key%value,*) doRing
		case("adderr")
			do i=key%nr1,key%nr2
				read(key%value,*) ObsSpec(i)%adderr
				print*,i,ObsSpec(i)%adderr
			enddo
		case default
			do i=1,nmol_data
				if(key%key.eq.molname(i)) then
					read(key%value,*) mixrat(i)
					includemol(i)=(mixrat(i).ge.0d0)
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
	real*8 tot,k,kap,Teq,f
	integer i
	
	Rplanet=Rplanet*Rjup
	Mplanet=Mplanet*Mjup

	if(Mp_from_logg) then
		Mplanet=(Rplanet**2)*(10d0**(loggPlanet))/Ggrav
		call output("Planet mass: " // dbl2string(Mplanet/Mjup,'(f8.3)') // " Mjup")
		call output("Planet mass: " // dbl2string(Mplanet/Mearth,'(f8.3)') // " Mearth")
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
	
	do i=1,nclouds
		Cloud(i)%rnuc=Cloud(i)%rnuc*micron
	enddo

	distance=distance*parsec
	
	orbit_inc=orbit_inc*pi/180d0

c	tot=0d0
c	do i=1,nmol
c		if(mixrat(i).gt.0d0) tot=tot+mixrat(i)
c	enddo
c	if(tot.gt.1d0) then
c		call output("Summed mixing ratio above 1. Renormalizing.")
c		mixrat=mixrat/tot
c	endif

	if(setsurfpressure) then
		Pmax=0d0
		do i=1,nmol
			Pmax=Pmax+mixrat(i)
		enddo
	endif

	if(fixnight2day) call ComputeNight2Day(.true.)
	tauUV=-1d0
	scaleUV=1d0

	return
	end


	subroutine ComputeNight2Day(usekap1)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 k,Teq,f,kap
	logical usekap1
	
c Use formalism from Koll (2022)
	k=1.2	! value for TRAPPIST-1b
c k is proportional to Rplanet and grav
	k=k*(Ggrav*Mplanet/(Rplanet)**2)/795d0
	k=k*Rplanet/(1.12d0*Rearth)
	if(usekap1) then
		kap=1d0
	else
		kap=tauLW/Pmax
	endif
	Teq=sqrt(Rstar/(2d0*Dplanet))*Tstar
	f=(kap**(1./3.)*Pmax*(Teq/600d0)**(-4./3.))
	night2day=f/(2d0*k+f)
	betaT=2d0/3d0-(5d0/12d0)*f/(k+f)
	if(.not.retrieval) then
		call output("night2day contrast: " // dbl2string(night2day,'(f7.4)'))
		call output("1D beta value:      " // dbl2string(betaT,'(f7.4)'))
	endif

	return
	end
	


	subroutine InitDens()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,k,n,nrsurf_tot
	real*8 g,dp,dz,P0(nr),T0(nr),pp,tt,mr0(nr,nmol),mm(nmol),yp1,ypn
	real*8,allocatable :: y2(:)
	character*10 names(nmol)
	integer imol(nmol)
		
	do i=1,nr
		mixrat_r(i,1:nmol)=mixrat(1:nmol)
	enddo
	
	if(mixratfile) then
		open(unit=20,file=TPfile,FORM="FORMATTED")
		read(20,*) n
		read(20,*) names(1:n)
		do j=1,n
			imol(j)=0
			do i=1,nmol_data
				if(names(j).eq.molname(i)) then
					imol(j)=i
					includemol(imol(j))=.true.
				endif
			enddo
			if(imol(j).eq.0) then
				call output("molecule " // trim(names(j)) // " not included yet in ARCiS")
				stop
			endif
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
		open(unit=20,file=TPfile,FORM="FORMATTED")
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
		if(Pmin.gt.Pmax/Psurf.or.nrsurf.le.1) then
			do i=1,nr
				P0(i)=exp(log(Pmin)+log(Pmax/Pmin)*real(i-1)/real(nr-1))
			enddo
		else
			nrsurf_tot=nrsurf+nr*(log(Psurf)/log(Pmax/Pmin))
			if(nrsurf_tot.gt.nr-5) nrsurf_tot=nr-5
			do i=1,nr-nrsurf_tot
				P0(i)=exp(log(Pmin)+log((Pmax/Psurf)/Pmin)*real(i-1)/real(nr-nrsurf_tot))
			enddo
			do i=nr-nrsurf_tot+1,nr
				P0(i)=exp(log(Pmax/Psurf)+log(Psurf)*real(i-nr+nrsurf_tot-1)/real(nrsurf_tot-1))
			enddo
			call sort(P0,nr)
		endif
		if(TPfile.ne.' ') then
			call regridlog(TPfile,P0,T0,nr)
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

	P(nr+1)=0d0
	do i=1,nr
		T(i)=T0(nr+1-i)
		P(i)=P0(nr+1-i)
	enddo

c	if(par_tprofile) call ComputeParamT(T)
	do i=1,nr
		if(free_tprofile) then
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
	Tsurface=T(1)

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
	allocate(P(nr+1))
	allocate(Hp(nr))
	allocate(nabla_ad(nr))
	allocate(grav(nr))
	allocate(mixrat_r(nr,nmol))
	allocate(mixrat_old_r(nr,nmol))
	allocate(cloud_dens(nr,max(nclouds,1)))

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
	use Struct3D
	use ARCiS_GGCHEM
	IMPLICIT NONE
	integer i
	character*100 homedir

	call SetupAtoms
	
	pargridfile=" "
	idum0=42
	randomseed=.true.

	iWolk=0
	
	particledir='./Particles/'

	transspec=.true.
	emisspec=.true.
	dotranshide=.false.

	Mplanet=1d0
	Rplanet=1d0
	Pplanet=10d0
	loggPlanet=2.5d0
	Mp_from_logg=.false.
	constant_g=.false.
	fixMMW=.false.
	MMW0=2.2
	Psurf=10d0
	
	Tstar=5777d0
	Rstar=1d0
	Mstar=1d0
	Dplanet=1d0
	logg=4.5d0
	blackbodystar=.false.
	retrievestar=.false.
	
	fDay=0.5d0
	betapow=1d0
	Kxx=1d0
	Kyy=1d0
	vxx=0d0
	powvxx=0d0
	night2day=0.5d0
	hotspotshift0=-1d5
	pole2eq=1d0
	fixnight2day=.false.
	tidallock=.true.
	
	orbit_P=-1d0
	orbit_e=0d0
	orbit_omega=0d0
	orbit_inc=90d0
	
	computeLC=.false.
	
	lam1=1d0
	lam2=15d0
	specres=10d0
	specres_LR=10d0
	
	useXS=.false.
	
	trend_compute=.false.

	mixrat=0d0
	includemol=.false.
	par_tprofile=.false.
	do_rayleigh=.true.
	Pswitch_mol=0d0
	abun_switch_mol=0d0
	setsurfpressure=.false.
	backgroundgas=.false.
	
	Nphot0=2500
	
	makeimage=.false.
	makemovie=.false.
	
	surfacetype='BLACK'
	surfacealbedo=0.5d0
	f_surface_water=0.6

	dochemistry=.false.
	elements_ARCiS= 'H He C N O Na Mg Si Fe Al Ca Ti S Cl K Li P V el'
	disequilibrium=.false.
	nfixmol=0
	fixmol_P=1d20
	Kzz=1d8
	metallicity=0d0
	condensates=.false.
	dosimplerainout=.false.
	COratio=0.5495407855762011
	SiOratio=0.06606931168616334
	NOratio=0.13803841820153123
	SOratio=0.026915346428939845
	TiScale=1d0
	inverseCOratio=.false.
	element_abun_file=' '
	complexKzz=.false.
	mixP=0d0
	vfrag=100d0	!cm/s
	
	vrot0=0d0
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
	
	secondary_atmosphere=.false.
	WaterWorld=.false.
	Poutgas=0d0
	Toutgas=0d0

c This parameter should be true! Only set this to false to reproduce earlier computations when
c  GGchem was still implemented slightly wrong.
	GGCHEM_P_iter=.true.

	GGCHEM_Tfast=700d0
	
	adiabatic_tprofile=.false.
	outflow=.false.

	Rring=1.5d0
	dRring=2.0d0
	tauRing=0.5d0
	doRing=.false.

	domakeai=.false.
	nai=1000
	dopostequalweights=.false.
	npew=-1
	nboot=1

	do3D=.false.
	init3D=.false.

	n3D=10
	nnu0=10
	nlong=36
	nlatt=19
	i3D=1
	
	deepredist=.false.
	deepredisttype='fixbeta'

	readFull3D=.false.
	computealbedo=.false.

	Kzz_deep=1d2
	Kzz_1bar=-1d0
	Kzz_P=0.5d0
	Kzz_contrast=-1d0

	computecontrib=.false.

	nC_PAH=25
	mixrat_PAH=0d0

	mixrat_optEC0=0d0
	rad_optEC=0.01
	Eg_optEC=1.0

	kappaUV0=-1d0

	instrument="ARIEL"
	instr_ntrans=1d0
	do i=1,n_instr
		instr_nobs(i)=i
	enddo

	nrstepchem=1
	nr_cloud=10

	do i=1,nclouds
		Cloud(i)%opacitytype=' '
		Cloud(i)%type=' '
		Cloud(i)%P=1d-4
		Cloud(i)%dP=10d0
		Cloud(i)%xi=2d0
		Cloud(i)%Pmax=1d10
		Cloud(i)%Pmin=0d0
		Cloud(i)%Ptau=1d0
		Cloud(i)%Phi=2d0
		Cloud(i)%coverage=1d0
		Cloud(i)%abun=1d0
		Cloud(i)%fmax=0d0
		Cloud(i)%porosity=0d0
		Cloud(i)%reff=1d0
		Cloud(i)%veff=0.1
		Cloud(i)%rpow=0d0
		Cloud(i)%Pref=1d0
		Cloud(i)%rnuc=0.001
		Cloud(i)%xm_bot=0d0
		Cloud(i)%blend=.true.
		Cloud(i)%haze=.false.
		Cloud(i)%condensates=.true.
		Cloud(i)%rainout=.false.
		Cloud(i)%coagulation=.true.
		Cloud(i)%computecryst=.false.
		Cloud(i)%mixrat=0d0
		Cloud(i)%tau=1d0
		Cloud(i)%lref=1d0
		Cloud(i)%cryst0=1d0
		Cloud(i)%e1_par=1.5
		Cloud(i)%e2_par=0.01
		Cloud(i)%hazetype='SOOT'
		Cloud(i)%nmat=1
		Cloud(i)%material=' '
		Cloud(i)%lnkfile=' '
		Cloud(i)%Kzz=-1d0
		Cloud(i)%Sigmadot=1d-17
		Cloud(i)%kappa=1d-2
		Cloud(i)%albedo=0.99d0
		Cloud(i)%kpow=4d0
		Cloud(i)%klam=1d0
		Cloud(i)%rho_mat=3.0
		Cloud(i)%nax=1
		Cloud(i)%globalKzz=.false.
	enddo
	nspike=0
	useDLMie=.false.

	isotope=0
	f_isotope=0d0

	retrieval=.false.
	doscaleR=.false.
	nscaleR=-1
	massprior=.false.
	useobsgrid=.false.
	log_emis=.true.
	model_err_abs=0d0
	model_err_rel=0d0
	do i=1,n_ret
		RetPar(i)%x0=-1d200
		RetPar(i)%dx=-1d0
		RetPar(i)%logscale=.false.
		RetPar(i)%squarescale=.false.
		RetPar(i)%opacitycomp=.true.
		RetPar(i)%increase=.false.
		RetPar(i)%xmin=0d0
		RetPar(i)%xmax=1d0		
	enddo
	do i=1,n_Par3D
		Par3D(i)%logscale=.false.
		Par3D(i)%multiply=.false.
	enddo
	par3Dsteepness=1d-4
	npop=30
	ngen=0
	gene_cross=.false.
	randomstart=.false.
	resume_multinest=.false.
	f_multinest=0.3d0
	tol_multinest=0.5d0
	nest_update=100
	const_eff_multinest=.false.
	retrievaltype='MN'
	writeWolk=.true.

	do i=1,nobs
		ObsSpec(i)%beta=1d0
		ObsSpec(i)%scale=1d0
		ObsSpec(i)%scaling=.false.
		ObsSpec(i)%dscale=-1d0
		ObsSpec(i)%spec=.true.
		ObsSpec(i)%i2d=0
		ObsSpec(i)%iphase=1
		ObsSpec(i)%slope=0d0
		ObsSpec(i)%adderr=0d0
		ObsSpec(i)%filter=' '
	enddo
	
	computeT=.false.
	isoFstar=.false.
	TeffP=600d0
	outputopacity=.false.
	forceEbalance=.false.
	exp_ad=7d0/5d0		! adiabatic exponent
	Tsurface0=-1d0

	call getenv('HOME',homedir) 

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
	
	scattering=.false.
	scattstar=.false.
	
	opacitymode=.false.
	opacitydir=trim(homedir) // '/ARCiS/Data/Opacities'

	nTom=100
	nPom=50
	Tmin=2.7d0
	Tmax=10000d0
	
	logTprofile=.true.
	taurexprofile=.false.
	taurexsmooth=1.2
	
	TP0=600d0
	dTP=0.1
	
	gammaT1=1.58e-1
	gammaT2=1.58e-1
	gamma_equal=.true.
	kappaT=3d-4
	betaT=1d0
	alphaT=1d0
	
	free_tprofile=.false.
	tau_Vpoint=1d0
	tau_IRpoint=1d0
	dT_Vpoint=2d0/7d0
	dT_IRpoint=2d0/7d0
	dTpoint=2d0/7d0
	chimax=1d0
	pos_dT_lowest=.false.	
	pos_dT=.false.	
	wiggle_err=-1d0

	maxTprofile=1d6
	
	maxiter=6
	miniter=5
	epsiter=3d-2

	maxchemtime=1d200
	
	fulloutput3D=.false.

	planetform=.false.
	planetform_fdust=0.1
	planetform_fplan=0.1
	planetform_Dmigrate=15d0
	planetform_Rend=-1d0
	planetform_Mstart=10d0
	planetform_SolidC=0d0
	planetform_Macc=1d-7
	Hydrogenloss=1d0

	do i=1,nPhotoReacts
		PhotoReacts(i)%react=0d0
		PhotoReacts(i)%product=0d0
		PhotoReacts(i)%f_eff=1d0
		PhotoReacts(i)%haze=0d0
		PhotoReacts(i)%nreact=0
		PhotoReacts(i)%atomic=.false.
	enddo

	return
	end

	
	subroutine InitObs()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j

c number of cloud/nocloud combinations
	ncc=1
	do i=1,nclouds
		if(Cloud(i)%coverage.ne.1d0) ncc=2**nclouds
	enddo
	do i=1,n_ret
		j=len_trim(RetPar(i)%keyword)
		if(j.gt.7) then
			if(RetPar(i)%keyword(j-7:j).eq.'coverage') ncc=2**nclouds
		endif
	enddo
	allocate(docloud(ncc,nclouds))
	allocate(cloudfrac(ncc))
	allocate(flux(0:ncc,nlam))
	allocate(obsA(0:ncc,nlam))
	allocate(obsA_split(nlam,2))
	allocate(tau1depth(ncc,nlam))
	allocate(cloudtau(ncc,nlam))
	allocate(phase(nphase,0:ncc,nlam))
	if(computealbedo) allocate(planet_albedo(nphase,nlam))

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
			n_ret=n_ret+1
			i=n_ret
			read(key%value,*) RetPar(i)%keyword
			if(RetPar(i)%keyword.eq.'tprofile') then
 				free_tprofile=.true.
				n_ret=n_ret+nTpoints-1
				if(freePT_fitT) then
 				do j=1,nTpoints
					RetPar(i+j-1)%keyword='Tpoint' // trim(int2string(j,'(i0.3)'))
					RetPar(i+j-1)%xmin=10d0
					RetPar(i+j-1)%xmax=10000d0
					RetPar(i+j-1)%logscale=.true.
				enddo
				else
 				do j=1,nTpoints
					RetPar(i+j-1)%keyword='dTpoint' // trim(int2string(j,'(i0.3)'))
					RetPar(i+j-1)%x0=2d0/7d0
					RetPar(i+j-1)%xmin=-4d0/7d0
					RetPar(i+j-1)%xmax=4d0/7d0
					RetPar(i+j-1)%logscale=.false.
				enddo
				endif
				if(freePT_fitP) then
				n_ret=n_ret+nTpoints-2
 				do j=2,nTpoints-1
					RetPar(i+nTpoints+j-2)%keyword='Ppoint' // trim(int2string(j,'(i0.3)'))
					RetPar(i+nTpoints+j-2)%xmin=pmin
					RetPar(i+nTpoints+j-2)%xmax=pmax
					RetPar(i+nTpoints+j-2)%logscale=.true.
					if(j.eq.2) then
						RetPar(i+nTpoints+j-2)%increase=.false.
					else
						RetPar(i+nTpoints+j-2)%increase=.true.
					endif
				enddo
				endif
			endif
		case("min","xmin")
			read(key%value,*) RetPar(i)%xmin
		case("max","xmax")
			read(key%value,*) RetPar(i)%xmax
		case("init","x")
			read(key%value,*) RetPar(i)%x0
		case("spread","width","dx")
			read(key%value,*) RetPar(i)%dx
		case("log","logscale")
			read(key%value,*) RetPar(i)%logscale
		case("square","squarescale")
			read(key%value,*) RetPar(i)%squarescale
		case("opacity","opacitycomp")
			read(key%value,*) RetPar(i)%opacitycomp
		case("increase")
			read(key%value,*) RetPar(i)%increase
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
		case("filter")
			ObsSpec(i)%filter=key%value
		case("beta","weight")
			read(key%value,*) ObsSpec(i)%beta
		case("scaling","scale")
			read(key%value,*) ObsSpec(i)%scaling
		case("dscale","sigscale")
			read(key%value,*) ObsSpec(i)%dscale
		case("i2d")
			read(key%value,*) ObsSpec(i)%i2d
			if(ObsSpec(i)%i2d.gt.n2d) n2d=ObsSpec(i)%i2d
		case("iphase")
			read(key%value,*) ObsSpec(i)%iphase
		case("slope")
			read(key%value,*) ObsSpec(i)%slope
		case("adderr")
			read(key%value,*) ObsSpec(i)%adderr
		case default
			call output("Keyword not recognised: " // trim(key%key2))
	end select
	
	return
	end

	subroutine ReadPar3D(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	i=key%nr1
	
	select case(key%key2)
		case("keyword","parameter")
			read(key%value,*) Par3D(i)%keyword
		case("min","xmin")
			read(key%value,*) Par3D(i)%xmin
		case("max","xmax")
			read(key%value,*) Par3D(i)%xmax
		case("log","logscale")
			read(key%value,*) Par3D(i)%logscale
		case("relative","multiply")
			read(key%value,*) Par3D(i)%multiply
		case default
			call output("Keyword not recognised: " // trim(key%key2))
	end select
	
	return
	end
	
	subroutine ReadPswitch(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	
	do i=1,nmol_data
		if(key%orkey2.eq.molname(i)) then
			if(i.le.nmol) then
				read(key%value,*) Pswitch_mol(i)
			endif
			return
		endif
	enddo
	call output("Molecule not recognised in Pswitch")
	stop
	
	return
	end
	
	subroutine ReadAbunSwitch(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	
	do i=1,nmol_data
		if(key%orkey2.eq.molname(i)) then
			if(i.le.nmol) then
				read(key%value,*) abun_switch_mol(i)
			endif
			return
		endif
	enddo
	call output("Molecule not recognised in abun_switch")
	stop
	
	return
	end
	
	subroutine ReadBackgroundgas(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	
	do i=1,nmol_data
		if(key%orkey2.eq.molname(i)) then
			if(i.le.nmol) then
				read(key%value,*) backgroundgas(i)
			endif
			return
		endif
	enddo
	call output("Molecule not recognised as background gas")
	stop
	
	return
	end
	
	subroutine ReadIsotope(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer imol,i
	
	do i=1,nmol_data
		if(key%orkey2.eq.molname(i)) exit
	enddo
	if(i.gt.nmol_data) then
		call output("Molecule not recognised in isotope")
		stop
	endif
	imol=i

	select case(key%key1)
		case("isotope")
			do i=1,nmol_data
				if(key%value.eq.molname(i)) exit
			enddo
			if(i.gt.nmol_data) then
				call output("Molecule not recognised in isotope")
				stop
			endif
			isotope(imol)=i
		case("f_isotope")
			read(key%value,*) f_isotope(imol)
		case default
			call output("Keyword not recognised: " // trim(key%key2))
	end select
	
	return
	end
	
	subroutine ReadFixMol(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i,j
	i=key%nr1
	if(i.gt.nfixmol) nfixmol=i
	
	select case(key%key2)
		case("name")
			fixmol_name(i)=trim(key%value)
			do j=1,nmol_data
				if(fixmol_name(i).eq.molname(j)) then
					ifixmol(i)=j
				endif
			enddo
		case("abun")
			read(key%value,*) fixmol_abun(i)
		case("p")
			read(key%value,*) fixmol_P(i)
		case default
			call output("Keyword not recognised: " // trim(key%key2))
	end select
	
	return
	end

	subroutine ReadPhotoChem(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	use AtomsModule
	IMPLICIT NONE
	type(SettingKey) key
	integer i,j
	j=key%nr1
	if(j.gt.nPhotoReacts) nPhotoReacts=j

	do i=1,N_atoms
		if(key%orkey2.eq.names_atoms(i)) then
			if(i.le.nmol) then
				select case(key%key1)
					case("photoreac","photoreactant")
						read(key%value,*) PhotoReacts(j)%react(i)
						if(.not.PhotoReacts(j)%atomic.and.PhotoReacts(j)%nreact.gt.0) then
							call output("Mixed atomic/molecular photoreaction not supported")
							stop
						endif
						PhotoReacts(j)%nreact=PhotoReacts(j)%nreact+1
						PhotoReacts(j)%atomic=.true.
				end select
			endif
			return
		endif
	enddo
	do i=1,nmol_data
		if(key%orkey2.eq.molname(i)) then
			if(i.le.nmol) then
				select case(key%key1)
					case("photoreac","photoreactant")
						read(key%value,*) PhotoReacts(j)%react(i)
						if(PhotoReacts(j)%atomic.and.PhotoReacts(j)%nreact.gt.0) then
							call output("Mixed atomic/molecular photoreaction not supported")
							stop
						endif
						PhotoReacts(j)%nreact=PhotoReacts(j)%nreact+1
						PhotoReacts(j)%atomic=.false.
					case("photoprod","photoproduct")
						read(key%value,*) PhotoReacts(j)%product(i)
				end select
			endif
			return
		endif
	enddo
	call output("Molecule not recognised in abun_switch")
	stop
	
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
		case("nobs","nr")
			read(key%value,*) instr_nobs(i)
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
				call checkfile(CIA(i)%filename)
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
	character*1000 readline,command

	outputdir='./outputARCiS/'

	ncla=2
1	continue
	readline=commandargs(ncla)
	if(readline(1:2).eq.'-o') then
		outputdir=commandargs(1+ncla)
		if(outputdir(len_trim(outputdir)-1:len_trim(outputdir)).ne.'/') then
			outputdir=trim(outputdir) // '/'
		endif
		ncla=ncla+1
		goto 1
	endif
	ncla=ncla+1
	if(readline(1:1).ne.' ') goto 1

	write(command,'("mkdir -p ",a)') trim(outputdir)
	call system(command)
	
	return
	end


	subroutine InitFreq()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 lam0,T0,Planck,tot,x,y,dy,dx,lminRT,lmaxRT
	integer i,j,ilam,nj,jlam
	logical truefalse
	
	lminRT=0.11d0*micron
	lmaxRT=47d0*micron

	nlam=0
	if(useobsgrid) then
		do i=1,nobs
			select case(ObsSpec(i)%type)
				case('tprofile','logtp','priors','prior')
					continue
				case('lightcurve')
					inquire(file=ObsSpec(i)%file,exist=truefalse)
					if(.not.truefalse) then
						call output("File does not exist" // trim(ObsSpec(i)%file))
						stop
					endif
					open(unit=30,file=ObsSpec(i)%file,FORM="FORMATTED")
11					read(30,*,err=11) nj
					do j=1,nj
						read(30,*)
					enddo
12					read(30,*,end=13,err=12) x,dx
					read(30,*)
					read(30,*)
					nlam=nlam+1
					goto 12
13					close(unit=30)					
				case default
					if(ObsSpec(i)%filter.ne.' ') then
						inquire(file=ObsSpec(i)%filter,exist=truefalse)
						if(.not.truefalse) then
							call output("File does not exist" // trim(ObsSpec(i)%filter))
							stop
						endif
						open(unit=30,file=ObsSpec(i)%filter,FORM="FORMATTED")
5						read(30,*,end=6,err=5) x,y
						nlam=nlam+1
						goto 5
6						close(unit=30)
					else
						inquire(file=ObsSpec(i)%file,exist=truefalse)
						if(.not.truefalse) then
							call output("File does not exist" // trim(ObsSpec(i)%file))
							stop
						endif
						open(unit=30,file=ObsSpec(i)%file,FORM="FORMATTED")
1						read(30,*,end=2,err=1) x,y,dy,dx
						nlam=nlam+1
						goto 1
2						close(unit=30)
					endif
			end select
		enddo
	else
		lam0=lam1
		nlam=nlam+1
		do while((lam0+lam0/specres).le.lam2)
			lam0=lam0+lam0/specres
			nlam=nlam+1
		enddo
	endif
	if(computeT.or.doRing) then
		lam0=lminRT
		nlam=nlam+1
		do while((lam0+lam0/specres_LR).le.lmaxRT)
			lam0=lam0+lam0/specres_LR
			nlam=nlam+1
		enddo
	endif
	allocate(lam(nlam))
	allocate(blam(2,nlam))
	allocate(freq(nlam))
	allocate(dfreq(nlam))
	allocate(dlam(nlam))
	allocate(RTgridpoint(nlam),computelam(nlam))
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(computelam,RTgridpoint,lam,nlam)
!$OMP DO SCHEDULE(STATIC)
		do i=1,nlam
			computelam(i)=.true.
			RTgridpoint(i)=.false.
			lam(i)=0d0
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	
	ilam=0
	if(useobsgrid) then
		do i=1,nobs
			select case(ObsSpec(i)%type)
				case('tprofile','logtp','priors','prior')
					continue
				case('lightcurve')
					open(unit=30,file=ObsSpec(i)%file,FORM="FORMATTED")
14					read(30,*,err=14) nj
					do j=1,nj
						read(30,*)
					enddo
15					read(30,*,end=16,err=15) x,dx
					read(30,*)
					read(30,*)
					ilam=ilam+1
					lam(ilam)=x*micron
	dx=100d0
					dx=x/dx
					dlam(ilam)=dx*micron
					goto 15
16					close(unit=30)					
				case default
					if(ObsSpec(i)%filter.ne.' ') then
						open(unit=30,file=ObsSpec(i)%filter,FORM="FORMATTED")
7						read(30,*,end=8,err=7) x,dx
						ilam=ilam+1
						lam(ilam)=x*micron
						dlam(ilam)=dx*micron
						goto 7
8						close(unit=30)
					else
						open(unit=30,file=ObsSpec(i)%file,FORM="FORMATTED")
3						read(30,*,end=4,err=3) x,y,dy,dx
						ilam=ilam+1
						lam(ilam)=x*micron
						dx=x/dx
						dlam(ilam)=dx*micron
						do jlam=1,ilam-1
							if(abs(lam(ilam)-lam(jlam)).lt.(dlam(jlam)*0.1d0).and.
     &							abs(dlam(jlam)-dlam(ilam))/(dlam(jlam)+dlam(ilam)).lt.0.1d0) then
     							ilam=ilam-1
     							goto 3
     						endif
     					enddo
						goto 3
4						close(unit=30)
					endif
			end select
		enddo
	else
		ilam=ilam+1
		lam(ilam)=lam1
		dlam(ilam)=lam(ilam)/specres
		do while((lam(ilam)+lam(ilam)/specres).le.lam2)
			ilam=ilam+1
			lam(ilam)=lam(ilam-1)+lam(ilam-1)/specres
			dlam(ilam)=lam(ilam)/specres
		enddo
	endif
	if(computeT.or.doRing) then
		ilam=ilam+1
		lam(ilam)=lminRT
		dlam(ilam)=-lam(ilam)/specres_LR
		do while((lam(ilam)+lam(ilam)/specres_LR).le.lmaxRT)
			ilam=ilam+1
			lam(ilam)=lam(ilam-1)+lam(ilam-1)/specres_LR
			dlam(ilam)=-lam(ilam)/specres_LR
		enddo
	endif
	nlam=ilam
	call sortw(lam,dlam,nlam)
	do ilam=1,nlam-1
		if(lam(ilam).eq.lam(ilam+1)) then
			if(ilam.eq.1) then
				lam(ilam+1)=(1d0-1d-3)*lam(ilam+1)+1d-3*lam(ilam+2)
			else
				lam(ilam)=(1d0-1d-3)*lam(ilam)+1d-3*lam(ilam-1)
			endif
		endif
	enddo
	RTgridpoint=.false.
	if(computeT.or.doRing) then
		do ilam=1,nlam
			if(dlam(ilam).lt.0d0) then
				RTgridpoint(ilam)=.true.
				dlam(ilam)=-dlam(ilam)
			endif
		enddo
	endif
	do i=1,nlam
		blam(1,i)=lam(i)/sqrt(1d0+dlam(i)/lam(i))
		blam(2,i)=lam(i)*sqrt(1d0+dlam(i)/lam(i))
		freq(i)=1d0/lam(i)
		dfreq(i)=abs(1d0/blam(1,i)-1d0/blam(2,i))
	enddo
	do i=1,nlam
		if(lam(i).lt.lam1) lam1=lam(i)
		if(lam(i).gt.lam2) lam2=lam(i)
	enddo

	return
	end


	subroutine ReadCloud(key)
	use GlobalSetup
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer j,j1,j2,i


	if(key%nr1.eq.0) then
		j1=1
		j2=nclouds
	else
		j1=key%nr1
		j2=key%nr1
	endif

	do j=j1,j2
	select case(key%key2)
		case("type")
			Cloud(j)%type=trim(key%value)
		case("opacitytype")
			Cloud(j)%opacitytype=trim(key%value)
		case("file")
			Cloud(j)%file=trim(key%value)
			call checkfile(Cloud(j)%file)
		case("fmax")
			read(key%value,*) Cloud(j)%fmax
		case("blend")
			read(key%value,*) Cloud(j)%blend
		case("porosity")
			read(key%value,*) Cloud(j)%porosity
		case("hazetype")
			Cloud(j)%hazetype=trim(key%value)
		case("pressure","p")
			read(key%value,*) Cloud(j)%P
		case("dp")
			read(key%value,*) Cloud(j)%dP
		case("ptop")
			read(key%value,*) Cloud(j)%Pmin
		case("pbottom")
			read(key%value,*) Cloud(j)%Pmax
		case("ptau")
			read(key%value,*) Cloud(j)%Ptau
		case("xi")
			read(key%value,*) Cloud(j)%xi
		case("dlogp")
			read(key%value,*) Cloud(j)%Phi
		case("tau")
			read(key%value,*) Cloud(j)%tau
		case("lam_ref")
			read(key%value,*) Cloud(j)%lref
		case("lam_kappa")
			read(key%value,*) Cloud(j)%klam
		case("pow_kappa")
			read(key%value,*) Cloud(j)%kpow
		case("albedo")
			read(key%value,*) Cloud(j)%albedo
		case("rnuc")
			read(key%value,*) Cloud(j)%rnuc
		case("pref")
			read(key%value,*) Cloud(j)%Pref
		case("pow_rad")
			read(key%value,*) Cloud(j)%rpow
		case("n")
			read(key%value,*) Cloud(j)%e1_par
		case("k")
			read(key%value,*) Cloud(j)%e2_par
		case("coverage")
			read(key%value,*) Cloud(j)%coverage
		case("haze")
			read(key%value,*) Cloud(j)%haze
		case("condensates")
			read(key%value,*) Cloud(j)%condensates
		case("mixrat")
			read(key%value,*) Cloud(j)%mixrat
		case("kzz")
			read(key%value,*) Cloud(j)%Kzz
		case("globalkzz")
			read(key%value,*) Cloud(j)%GlobalKzz
		case("sigmadot","nucleation")
			read(key%value,*) Cloud(j)%Sigmadot
		case("xm_bot")
			read(key%value,*) Cloud(j)%xm_bot
		case("reff")
			read(key%value,*) Cloud(j)%reff
		case("veff")
			read(key%value,*) Cloud(j)%veff
		case("kappa")
			read(key%value,*) Cloud(j)%kappa
		case("rainout")
			read(key%value,*) Cloud(j)%rainout
		case("computecryst")
			read(key%value,*) Cloud(j)%computecryst
		case("coagulation")
			read(key%value,*) Cloud(j)%coagulation
		case("cryst")
			read(key%value,*) Cloud(j)%cryst0
		case("abun")
			i=key%nr2
			if(i.gt.Cloud(j)%nmat) Cloud(j)%nmat=i
			read(key%value,*) Cloud(j)%abun(i)
		case("rho_mat")
			i=key%nr2
			if(i.gt.Cloud(j)%nmat) Cloud(j)%nmat=i
			read(key%value,*) Cloud(j)%rho_mat(i)
		case("material")
			i=key%nr2
			if(i.gt.Cloud(j)%nmat) Cloud(j)%nmat=i
			Cloud(j)%material(i)=trim(key%value)
		case("lnkfile")
			i=key%nr2
			if(i.gt.Cloud(j)%nmat) Cloud(j)%nmat=i
			Cloud(j)%lnkfile(i,1)=trim(key%value)
			call checkfile(Cloud(j)%lnkfile(i,1))
			Cloud(j)%nax(i)=1
		case("lnkfilex")
			i=key%nr2
			if(i.gt.Cloud(j)%nmat) Cloud(j)%nmat=i
			Cloud(j)%lnkfile(i,1)=trim(key%value)
			call checkfile(Cloud(j)%lnkfile(i,1))
			Cloud(j)%nax(i)=3
		case("lnkfiley")
			i=key%nr2
			if(i.gt.Cloud(j)%nmat) Cloud(j)%nmat=i
			Cloud(j)%lnkfile(i,2)=trim(key%value)
			call checkfile(Cloud(j)%lnkfile(i,2))
			Cloud(j)%nax(i)=3
		case("lnkfilez")
			i=key%nr2
			if(i.gt.Cloud(j)%nmat) Cloud(j)%nmat=i
			Cloud(j)%lnkfile(i,3)=trim(key%value)
			call checkfile(Cloud(j)%lnkfile(i,3))
			Cloud(j)%nax(i)=3
		case default
			call output("Unknown cloud keyword: " // trim(key%key2))
			stop
	end select
	enddo

	return
	end
	


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
	
	subroutine compactname(name,cname)
	IMPLICIT NONE
	character*500 name,cname
	integer i,n,j
	n=len_trim(name)
	j=1
	cname=' '
	do i=1,n
		if(name(i:i).ne.' '.and.name(i:i).ne.'_'.and.name(i:i).ne.'-') then
			cname(j:j)=name(i:i)
			j=j+1
		endif
	enddo
	
	return
	end

	
	subroutine ReadPlanetNameCSV()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,dR1,dR2,dM1,dM2,dsig
	character*500 name,namestar,cname,cplanetname
	integer i,n
	character*10 Zc
	character*1000 key(1000),value(1000)
	character*6000 line

	open(unit=72,file=planetparameterfile,FORM="FORMATTED")
	read(72,'(a6000)') line
	call getkeys(line,key,n)

	call compactname(planetname,cplanetname)

	dMp_prior=0d0
1	read(72,'(a6000)',end=2) line
	call getkeys(line,value,n)
	do i=1,n
		if(value(i).ne.' ') then
			select case(key(i))
				case('Planet_Name','name')
					name=trim(value(i))
				case('Star_Temperature_[K]','star_teff')
					read(value(i),*) Tstar
				case('Star_Metallicity','star_metallicity')
					read(value(i),*) Zc
				case('Star_Mass_[Ms]','star_mass')
					read(value(i),*) Mstar
				case('Star_Radius_[Rs]','star_radius')
					read(value(i),*) Rstar
				case('Star_Distance_[pc]','star_distance')
					read(value(i),*) Distance
				case('Planet_Semi-major_Axis_[m]')
					read(value(i),*) Dplanet
					Dplanet=Dplanet*1d2/AU
				case('semi_major_axis')
					read(value(i),*) Dplanet
				case('Planet_Mass_[Me]')
					read(value(i),*) Mp_prior
					Mp_prior=Mp_prior*Mearth/MJup
				case('mass')
					read(value(i),*) Mp_prior
				case('mass_error_min','mass_error_max')
					read(value(i),*) x
					dMp_prior=sqrt((dMp_prior**2+x**2)/2d0)
				case('Planet_Radius_[Re]')
					read(value(i),*) Rplanet
					Rplanet=Rplanet*Rearth/RJup
				case('radius')
					read(value(i),*) Rplanet
				case('Planet_Period_[days]','orbital_period')
					read(value(i),*) orbit_P
					orbit_P=orbit_P*86400d0
				case('Eccentricity','eccentricity')
					read(value(i),*) orbit_e
				case('Planet_Temperature_[K]','temp_calculated')
					read(value(i),*) TP0
			end select
		endif
		logg=log10(Ggrav*(Mstar*Msun)/((Rstar*Rsun)**2))
	enddo
	call compactname(name,cname)

	if(Mp_prior.le.0d0) then
		massprior=.false.
	else
		Mplanet=Mp_prior
	endif
	
	if(trim(cname).eq.trim(cplanetname)) then
		close(unit=72)
		call output("Stellar mass:     " // dbl2string(Mstar,'(f9.4)') // "Msun")
		call output("Stellar T:        " // dbl2string(Tstar,'(f9.4)') // "K")
		call output("Stellar radius:   " // dbl2string(Rstar,'(f9.4)') // "Rsun")
		call output("Stellar logg:     " // dbl2string(logg,'(f9.4)'))
		call output("Stellar distance: " // dbl2string(Distance,'(f9.4)') // "pc")
		call output("Metallicity:      " // dbl2string(metallicity,'(f9.4)'))
		call output("Planet orbit:     " // dbl2string(Dplanet,'(f7.4)') // "AU")
		call output("Planet radius:    " // dbl2string(Rplanet,'(f9.4)') // "Rjup")
		call output("Planet mass:      " // dbl2string(Mplanet,'(f9.4)') // "Mjup")
		call output("Planet T:         " // dbl2string(TP0,'(f9.4)') // "K")
		call output("Blackbody T:      " // dbl2string(sqrt(Rstar*Rsun/(2d0*Dplanet*AU))*Tstar,'(f9.4)') // "K")
		call output("Orbital period:   " // dbl2string(orbit_P/86400d0,'(f9.4)') // "days")
			
		open(unit=73,file=trim(outputdir) // "ListParameters",FORM="FORMATTED")
		write(73,'(a)') "Stellar mass:     " // dbl2string(Mstar,'(f9.4)') // "Msun"
		write(73,'(a)') "Stellar T:        " // dbl2string(Tstar,'(f9.4)') // "K"
		write(73,'(a)') "Stellar radius:   " // dbl2string(Rstar,'(f9.4)') // "Rsun"
		write(73,'(a)') "Stellar logg:     " // dbl2string(logg,'(f9.4)')
		write(73,'(a)') "Stellar distance: " // dbl2string(Distance,'(f9.4)') // "pc"
		write(73,'(a)') "Metallicity:      " // dbl2string(metallicity,'(f9.4)')
		write(73,'(a)') "Planet orbit:     " // dbl2string(Dplanet,'(f9.4)') // "AU"
		write(73,'(a)') "Planet radius:    " // dbl2string(Rplanet,'(f9.4)') // "Rjup"
		write(73,'(a)') "Planet mass:      " // dbl2string(Mplanet,'(f9.4)') // "Mjup"
		write(73,'(a)') "Planet T:         " // dbl2string(TP0,'(f9.4)') // "K"
		write(73,'(a)') "Blackbody T:      " // dbl2string(sqrt(Rstar*Rsun/(2d0*Dplanet*AU))*Tstar,'(f9.4)') // "K"
		write(73,'(a)') "Orbital period:   " // dbl2string(orbit_P/86400d0,'(f9.4)') // "days"
		close(unit=73)

		dR1=Rplanet*0.1
		dR2=Rplanet*0.1
		dM1=Mplanet*0.1
		dM2=Mplanet*0.2

		Mp_prior=Mplanet
		if(dMp_prior.le.0d0) then
			dMp_prior=sqrt(dM1**2+dM2**2)
		endif
		dsig=20d0
		do i=1,n_ret
			select case(RetPar(i)%keyword)
				case("Rp","rp","RP")
					RetPar(i)%x0=Rplanet
					RetPar(i)%xmin=max(0d0,Rplanet-Rp_range*dR1)
					RetPar(i)%xmax=Rplanet+Rp_range*dR2
					if(RetPar(i)%xmin*Rjup.lt.0.1d0*Rearth) RetPar(i)%xmin=0.1d0*Rearth/Rjup
		call output("Minimum radius:      " // dbl2string(RetPar(i)%xmin,'(f7.2)') // "Rjup")
		call output("Maximum radius:      " // dbl2string(RetPar(i)%xmax,'(f7.2)') // "Rjup")
				case("Mp","mp","MP")
					RetPar(i)%x0=Mplanet
					RetPar(i)%xmin=max(0d0,Mplanet-dsig*dM1)
					RetPar(i)%xmax=Mplanet+dsig*dM2
					if(RetPar(i)%xmin*Mjup.lt.0.1d0*Mearth) RetPar(i)%xmin=0.1d0*Mearth/Mjup
		call output("Minimum mass:        " // dbl2string(RetPar(i)%xmin,'(f7.2)') // "Mjup")
		call output("Maximum mass:        " // dbl2string(RetPar(i)%xmax,'(f7.2)') // "Mjup")
				case("loggP","loggp")
					RetPar(i)%x0=log10(Ggrav*(Mplanet*Mjup)/((Rplanet*Rjup)**2))
					RetPar(i)%xmin=max(0.1,log10(Ggrav*(max(Mjup*(Mplanet-dsig*dM1),0.1d0*Mearth))/(((Rplanet+dsig*dR2)*Rjup)**2)))
					RetPar(i)%xmax=log10(Ggrav*((Mplanet+dsig*dM2)*Mjup)/((max((Rplanet-dsig*dR1)*Rjup,0.1d0*Rearth))**2))
		call output("Minimum logg:        " // dbl2string(RetPar(i)%xmin,'(f7.2)'))
		call output("Maximum logg:        " // dbl2string(RetPar(i)%xmax,'(f7.2)'))
				case("Rstart","rstart")
					RetPar(i)%x0=Dplanet
					RetPar(i)%xmin=max(RetPar(i)%xmin,Dplanet)
		call output("Minimum Rstart:      " // dbl2string(RetPar(i)%xmin,'(f7.2)'))
				case("RendMigrate","rendmigrate","Rendmigrate","rendMigrate")
					RetPar(i)%x0=Dplanet
					RetPar(i)%xmin=max(RetPar(i)%xmin,Dplanet)
		call output("Minimum RendMigrate: " // dbl2string(RetPar(i)%xmin,'(f7.2)'))
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
	character*6000 line
	character*1000 key(1000)
	integer icomma,i,l,n,j

	icomma=1
	l=0
	i=1
	do while(icomma.ne.0)
		icomma=index(line(l+1:6000),',')
		if(icomma.ne.0) then
			key(i)=line(l+1:l+icomma-1)
			do j=1,icomma-1
				if(key(i)(j:j).eq.' ') key(i)(j:j)='_'
			enddo
			i=i+1
			l=l+icomma
		else
			key(i)=trim(line(l+1:l+1000))
			do j=1,len_trim(key(i))
				if(key(i)(j:j).eq.' ') key(i)(j:j)='_'
			enddo
			i=i+1
		endif
	enddo
	n=i-1
	
	return
	end

	
	
	subroutine ReadPlanetName()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,dR1,dR2,dM1,dM2,dsig
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
	
	open(unit=72,file=planetparameterfile,FORM="FORMATTED")
	read(72,*)
1	read(72,*,end=2) name,Tstar,x,x,Zc,x,x,Mstar,x,x,Rstar,x,x,logg,x,x,x,x,x,orbit_P,orbit_e,x,x,Dplanet,x,x,
     &					Mp_prior,dM1,dM2,Rplanet,dR1,dR2
	if(Zc.eq.'-1') then
		metallicity=0d0
	else
		read(Zc,*) metallicity
	endif
	i=len_trim(name)
	if(name(i:i).eq.'b') then
		name(i:i)=' '
	endif
	orbit_P=orbit_P*86400d0
	if(Mp_prior.lt.0d0) then
		massprior=.false.
	else if(Mp_prior.eq.0d0) then
		if(dM2.eq.0d0) then
			Mplanet=dM1
			dM2=dM1
		else
			massprior=.false.
		endif
	else
		Mplanet=Mp_prior
	endif
	Mp_prior=Mplanet
	dMp_prior=sqrt(dM1**2+dM2**2)
	if(trim(planetname).eq.trim(name)) then	
		close(unit=72)
		call output("Stellar mass:   " // dbl2string(Mstar,'(f9.4)') // "Msun")
		call output("Stellar T:      " // dbl2string(Tstar,'(f9.4)') // "K")
		call output("Stellar radius: " // dbl2string(Rstar,'(f9.4)') // "Rsun")
		call output("Stellar logg:   " // dbl2string(logg,'(f9.4)'))
		call output("Metallicity:    " // dbl2string(metallicity,'(f9.4)'))
		call output("Planet orbit:   " // dbl2string(Dplanet,'(f7.4)') // "AU")
		call output("Planet radius:  " // dbl2string(Rplanet,'(f9.4)') // "Rjup")
		call output("Planet mass:    " // dbl2string(Mplanet,'(f9.4)') // "Mjup")
		call output("Mass uncert.:   " // dbl2string(dMp_prior,'(f9.4)') // "Mjup")
		call output("Blackbody T:    " // dbl2string(sqrt(Rstar*Rsun/(2d0*Dplanet*AU))*Tstar,'(f9.4)') // "K")
		call output("Orbital period: " // dbl2string(orbit_P/86400d0,'(f9.4)') // "days")

		open(unit=73,file=trim(outputdir) // "ListParameters",FORM="FORMATTED")
		write(73,'(a)') "Stellar mass:     " // dbl2string(Mstar,'(f9.4)') // "Msun"
		write(73,'(a)') "Stellar T:        " // dbl2string(Tstar,'(f9.4)') // "K"
		write(73,'(a)') "Stellar radius:   " // dbl2string(Rstar,'(f9.4)') // "Rsun"
		write(73,'(a)') "Stellar logg:     " // dbl2string(logg,'(f9.4)')
		write(73,'(a)') "Stellar distance: " // dbl2string(Distance,'(f9.4)') // "pc"
		write(73,'(a)') "Metallicity:      " // dbl2string(metallicity,'(f9.4)')
		write(73,'(a)') "Planet orbit:     " // dbl2string(Dplanet,'(f9.4)') // "AU"
		write(73,'(a)') "Planet radius:    " // dbl2string(Rplanet,'(f9.4)') // "Rjup"
		write(73,'(a)') "Planet mass:      " // dbl2string(Mplanet,'(f9.4)') // "Mjup"
		write(73,'(a)') "Blackbody T:      " // dbl2string(sqrt(Rstar*Rsun/(2d0*Dplanet*AU))*Tstar,'(f9.4)') // "K"
		write(73,'(a)') "Orbital period:   " // dbl2string(orbit_P/86400d0,'(f9.4)') // "days"
		close(unit=73)

		dsig=20d0
		do i=1,n_ret
			select case(RetPar(i)%keyword)
				case("Rp","rp","RP")
					RetPar(i)%x0=Rplanet
					RetPar(i)%xmin=max(0d0,Rplanet-Rp_range*dR1)
					RetPar(i)%xmax=Rplanet+Rp_range*dR2
					if(RetPar(i)%xmin*Rjup.lt.0.1d0*Rearth) RetPar(i)%xmin=0.1d0*Rearth/Rjup
		call output("Minimum radius:      " // dbl2string(RetPar(i)%xmin,'(f7.2)') // "Rjup")
		call output("Maximum radius:      " // dbl2string(RetPar(i)%xmax,'(f7.2)') // "Rjup")
				case("Mp","mp","MP")
					RetPar(i)%x0=Mplanet
					RetPar(i)%xmin=max(0d0,Mplanet-dsig*dM1/4d0)
					RetPar(i)%xmax=Mplanet+dsig*dM2/4d0
					if(RetPar(i)%xmin*Mjup.lt.0.1d0*Mearth) RetPar(i)%xmin=0.1d0*Mearth/Mjup
		call output("Minimum mass:        " // dbl2string(RetPar(i)%xmin,'(f7.2)') // "Mjup")
		call output("Maximum mass:        " // dbl2string(RetPar(i)%xmax,'(f7.2)') // "Mjup")
				case("loggP","loggp")
					RetPar(i)%x0=log10(Ggrav*(Mplanet*Mjup)/((Rplanet*Rjup)**2))
					RetPar(i)%xmin=max(0.1,log10(Ggrav*(max(Mjup*(Mplanet-dsig*dM1/4d0),0.1d0*Mearth))/(((Rplanet+dsig*dR2)*Rjup)**2)))
					RetPar(i)%xmax=log10(Ggrav*((Mplanet+dsig*dM2/4d0)*Mjup)/((max((Rplanet-dsig*dR1)*Rjup,0.1d0*Rearth))**2))
		call output("Minimum logg:        " // dbl2string(RetPar(i)%xmin,'(f7.2)'))
		call output("Maximum logg:        " // dbl2string(RetPar(i)%xmax,'(f7.2)'))
				case("Rstart","rstart")
					RetPar(i)%x0=Dplanet
					RetPar(i)%xmin=max(RetPar(i)%xmin,Dplanet)
		call output("Minimum Rstart:      " // dbl2string(RetPar(i)%xmin,'(f7.2)'))
				case("RendMigrate","rendmigrate","Rendmigrate","rendMigrate")
					RetPar(i)%x0=Dplanet
					RetPar(i)%xmin=max(RetPar(i)%xmin,Dplanet)
		call output("Minimum RendMigrate: " // dbl2string(RetPar(i)%xmin,'(f7.2)'))
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
	

	
	real*8 function HG(g,theta)
	IMPLICIT NONE
	real*8 g,theta,t,tot,f
	integer i,n

c	HG=(1d0-g**2)/((1d0-2d0*g*cos(theta)+g**2)**(3.0/2.0))/2d0

	n=1000
	tot=0d0
	HG=0d0
	do i=1,n
		t=theta+(real(i-1)/real(n-1)-0.5)*3.1415926536/180d0
		f=(1d0-g**2)/((1d0-2d0*g*cos(t)+g**2)**(3.0/2.0))/2d0
		HG=HG+f*sin(t)
		tot=tot+sin(t)
	enddo
	HG=HG/tot

	return
	end


	subroutine InitFormation(Ms,Ts,Rs,frac_SolidC,Macc_in)
	use FormationModule
	use Constants
	IMPLICIT NONE
	real*8 Ms,frac_SolidC,Macc_in,Ts,Rs

	Mstar=Ms
	Lstar=Lsun*(Rs/Rsun)**2*(Ts/5777d0)**4
	call SetupAtoms

	call SetupPPdisk(frac_SolidC,Macc_in)

	return
	end
	
