	subroutine Init()
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey),pointer :: key,first
	integer i,j,omp_get_max_threads,omp_get_thread_num
	character*500 TPfile

	TPfile=' '

	idum=-42
#ifdef USE_OPENMP
	j=omp_get_max_threads()
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(j)
!$OMP DO SCHEDULE(STATIC,j)
	do i=1,j
		idum=-42-omp_get_thread_num()
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
#endif

	allocate(key)
	first => key

	call GetKeywords(key)
c Count the number of zones, particles, and stars
c allocate the arrays
	key => first%next

	do_cia=.true.
	call CountStuff(key)

	call SetDefaults

	key => first%next

	do while(.not.key%last)

	select case(key%key1)
		case("part")

		case("nr")
			read(key%value,*) nr
		case("mp")
			read(key%value,*) Mplanet
		case("rp")
			read(key%value,*) Rplanet
		case("retrieval")
			read(key%value,*) retrieval
		case("outputopacity","writeopacity")
			read(key%value,*) outputopacity
		case("obs")
			call ReadObs(key)
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
		case("eps","epsck")
			read(key%value,*) epsCk
		case("specres")
			read(key%value,*) specres
		case("tpfile")
			read(key%value,'(a)') TPfile
		case("ng")
			read(key%value,*) ng
		case("distance")
			read(key%value,*) distance
		case default
			do i=1,48
				if(key%key.eq.molname(i)) then
					read(key%value,*) mixrat(i)
					goto 1
				endif
			enddo
			call output("Keyword not recognised: " // trim(key%key1))
			stop
1			continue
	end select

	key => key%next
	
	enddo

	call ConvertUnits()

	call InitFreq()
	call InitDens(TPfile)
	call InitObs()

	allocate(opac(nr,nlam,ng))

	call output("==================================================================")

	call ReadHITRAN()

	call output("==================================================================")
	
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
	lam1=lam1*micron
	lam2=lam2*micron
	
	distance=distance*parsec
	
	tot=0d0
	do i=1,nmol
		if(mixrat(i).gt.0d0) tot=tot+mixrat(i)
	enddo
	if(tot.gt.1d0) then
		call output("Summed mixing ratio above 1. Renormalizing.")
		mixrat=mixrat/tot
	endif
	
	return
	end


	subroutine InitDens(TPfile)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i
	real*8 g,dp,dz,P0(nr),T0(nr)
	character*500 TPfile
	
	allocate(dens(nr))
	allocate(Ndens(nr))
	allocate(R(nr+1))
	allocate(T(nr))
	allocate(P(nr))
	
	do i=1,nr
		P0(i)=exp(log(Pmin)+log(Pmax/Pmin)*real(i-1)/real(nr-1))
	enddo
	if(TPfile.ne.' ') then
		call regridlog(TPfile,P0,T0,nr)
	else
		do i=1,nr
			T0(i)=exp(log(270d0)+log(10d0)*(real(i-1)/real(nr-1)))
		enddo
	endif

	do i=1,nr
		T(i)=T0(nr+1-i)
		P(i)=P0(nr+1-i)
	enddo

	return
	end	

	subroutine SetDefaults()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i
	character*100 homedir
	
	Mplanet=1d0
	Rplanet=1d0
	
	lam1=1d0
	lam2=15d0
	specres=10d0

	mixrat=-1d0
	
	do i=1,nobs
		obs(i)%type='EMIS'
		obs(i)%filename=' '
	enddo

	retrieval=.false.
	outputopacity=.false.
	nr=20
	
	Pmin=1d-5
	Pmax=10d0

	call getenv('HOME',homedir) 

	HITRANfile=trim(homedir) // '/HITRAN/HITRAN2012.par'
	
	ng=100
	epsCk=0.25d0
	
	distance=10d0
	
	cutoff_abs=1d200
	cutoff_lor=1d200
	
	return
	end

	
	subroutine InitObs()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs
	
	do iobs=1,nobs
		allocate(obs(iobs)%lam(nlam))
		allocate(obs(iobs)%flux(nlam))
	enddo
	
	return
	end
	
	
	subroutine ReadObs(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	i=key%nr1
	
	select case(key%key2)
		case("type")
			read(key%value,'(a)') obs(i)%type
		case("file")
			read(key%value,'(a)') obs(i)%filename
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
				read(key%value,'(a)') CIA(i)%filename
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

	outputdir='./outputSPARC/'

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
	real*8 lam0
	integer i
	
	lam0=lam1
	nlam=1
	do while(lam0.le.lam2)
		lam0=lam0+lam0/specres
		nlam=nlam+1
	enddo
	
	allocate(lam(nlam))
	allocate(freq(nlam))
	
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

	return
	end

