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
	call CountStuff(key)

	call SetDefaults

	key => first%next

	do while(.not.key%last)

	select case(key%key1)
		case("molecule","mol")
			call ReadMol(key)
		case("part")

		case("nr")
			read(key%value,*) nrad
		case("mp")
			read(key%value,*) Mplanet
		case("rp")
			read(key%value,*) Rplanet
		case("retrieval")
			read(key%value,*) Rplanet
		case("obs")
			call ReadObs(key)
		case("lam")
			if(key%nr1.eq.1) read(key%value,*) lam1
			if(key%nr1.eq.2) read(key%value,*) lam2
		case("specres")
			read(key%value,*) specres
		case("tpfile")
			read(key%value,'(a)') TPfile
		case default
			call output("Keyword not recognised: " // trim(key%key1))
			stop
	end select

	key => key%next
	
	enddo

	call ConvertUnits()

	call InitFreq()
	call InitDens(TPfile)

	call output("==================================================================")
	
	return
	end


	subroutine ConvertUnits()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	
	Rplanet=Rplanet*Rjup
	Mplanet=Mplanet*Mjup
	lam1=lam1*micron
	lam2=lam2*micron
	
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
			T0(i)=exp(log(4000d0)-log(20d0)*(real(i-1)/real(nr-1)))
		enddo
	endif

	do i=1,nr
		T(i)=T0(nr+1-i)
		P(i)=P0(nr+1-i)
	enddo

	g=Ggrav*Mplanet/Rplanet**2
	R(1)=Rplanet
	do i=1,nr
		if(i.eq.1) then
			dp=P(1)-P(2)
		else if(i.eq.nr) then
			dp=P(nr-1)-P(nr)
		else
			dp=(P(i-1)-P(i+1))/2d0
		endif
		dens(i)=P(i)/(kb*T(i))
		dz=dp/(dens(i)*mp*mu*g)
		R(i+1)=R(i)+dz
	enddo

	open(unit=50,file=trim(outputdir) // 'densityprofile.dat',RECL=1000)
	do i=1,nr
		write(50,*) sqrt(R(i)*R(i+1)),dens(i),T(i),P(i)
	enddo
	close(unit=20)

	return
	end	

	subroutine SetDefaults()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i
	
	Mplanet=1d0
	Rplanet=1d0
	
	lam1=1d0
	lam2=15d0
	specres=10d0
	
	do i=1,nobs
		obs(i)%type='EMIS'
		obs(i)%filename=' '
	enddo

	retrieval=.false.
	nr=100
	
	Pmin=1d-7
	Pmax=1d0

	HITRANfile='~/HITRAN/HITRAN2012.par'
		
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

	
	subroutine ReadMol(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	i=key%nr1
	
	select case(key%key2)
		case("abun")
			read(key%value,*) abun(i)
		case default
			call output("Keyword not recognised: " // trim(key%key2))
	end select
	
	return
	end

	
	subroutine GetOutputDir
	use GlobalSetup
	IMPLICIT NONE
	integer ncla
	character*500 readline,command

	outputdir='./outputELMO/'

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
		freq(i)=clight/lam(i)
	enddo

	return
	end
	