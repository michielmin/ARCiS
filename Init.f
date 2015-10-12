c==============================================================================	
c Module for reading keywords
c==============================================================================	
	module ReadKeywords
	IMPLICIT NONE

	type SettingKey
		character*100 key1,key2,value,key
		integer nr1,nr2
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
	call get_key_value(readline,key%key,key%key1,key%key2,key%value,key%nr1,key%nr2)

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
	subroutine get_key_value(line,key,key1,key2,value,nr1,nr2)
	IMPLICIT NONE
	character*1000 line
	character*100 key,key1,key2,value
	integer i,nr1,nr2,ikey1,ikey2
	
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
	character*100 key
	integer nr,i,n
	
	n=len_trim(key)
	i=n
1	read(key(i:n),*,err=2) nr
	i=i-1
	if(i.eq.0) goto 2
	goto 1
2	continue
	if(i.eq.n) then
		nr=0
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
	character*500 homedir,h2h2file,h2hefile,h2ch4file,TPfile
	character*10 names(59)
	logical existh2h2,existh2he,existh2ch4,mixratfile

	key => firstkey

	nmol=1
	nobs=0
	ncia=0
	nclouds=0
	nretr=0
	j=0
	mixratfile=.false.
	do while(.not.key%last)
		select case(key%key1)
			case("obs")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nobs) nobs=key%nr1
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
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nclouds) nclouds=key%nr1
			case("retrieval")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nretr) nretr=key%nr1
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

	allocate(mixrat(nmol))
	allocate(includemol(nmol))
	allocate(Cloud(max(nclouds,1)))

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
				h2h2file=trim(homedir) // '/CIA/H2-H2_2011.cia'
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
				h2hefile=trim(homedir) // '/CIA/H2-He_2011.cia'
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
	integer i,j,omp_get_max_threads,omp_get_thread_num
	logical mixratfile
	character*500 TPfile
	real*8 tot,tot2,theta,Planck

	TPfile=' '
	mixratfile=.false.

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
		case("nr")
			read(key%value,*) nr
		case("mp")
			read(key%value,*) Mplanet
		case("rp")
			read(key%value,*) Rplanet
		case("rstar")
			read(key%value,*) Rstar
		case("tstar")
			read(key%value,*) Tstar
		case("logg")
			read(key%value,*) logg
		case("dp","dplanet")
			read(key%value,*) Dplanet
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
		case("tpfile")
			read(key%value,'(a)') TPfile
		case("mixratfile")
			read(key%value,*) mixratfile
		case("gridtpfile")
			read(key%value,*) gridTPfile
		case("tp")
			read(key%value,*) TP0
		case("dtp")
			read(key%value,*) dTP
		case("ng")
			read(key%value,*) ng
		case("distance")
			read(key%value,*) distance
		case("hitemp")
			read(key%value,*) HITEMP
		case("cloud")
			call ReadCloud(key)
		case("scattering")
			read(key%value,*) scattering
		case("scattstar")
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
		case("chemistry")
			read(key%value,*) dochemistry
		case("metallicity")
			read(key%value,*) metallicity
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

	key => key%next
	
	enddo

	call ConvertUnits()

	call InitFreq()

	if(opacitydir(len_trim(opacitydir)-1:len_trim(opacitydir)).ne.'/') then
		opacitydir=trim(opacitydir) // '/'
	endif
	if(opacitymode) then
		call InitOpacityMode()
	else
		call InitDens(TPfile,mixratfile)
		call InitObs()
	endif

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
	do j=2,180
		Rayleigh%IF11(j)=Rayleigh%IF11(j-1)+sintheta(j)*Rayleigh%F11(j)
	enddo

	allocate(Fstar(nlam))
	call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam)
	Fstar=Fstar*pi*Rstar**2*pi/3.336e11

	call output("==================================================================")

	if(compute_opac) then
		call ReadData()
	else
		do i=1,nmol
			if(includemol(i)) call InitReadOpacityFITS(i)
		enddo
	endif

	call ReadDataCIA()

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

	Rstar=Rstar*Rsun
	Dplanet=Dplanet*AU
	
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


	subroutine InitDens(TPfile,mixratfile)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,n
	real*8 g,dp,dz,P0(nr),T0(nr),pp,tt,mr0(nr,nmol),mm(nmol)
	character*500 TPfile
	character*10 names(nmol)
	integer imol(nmol)
	logical mixratfile
	
	allocate(dens(nr))
	allocate(Ndens(nr))
	allocate(R(nr+1))
	allocate(T(nr))
	allocate(P(nr))
	allocate(Hp(nr))
	allocate(mixrat_r(nr,nmol))
	allocate(cloud_dens(nr,max(nclouds,1)))

	do i=1,nr
		mixrat_r(i,1:nmol)=mixrat(1:nmol)
	enddo
	
	if(mixratfile) then
		open(unit=20,file=TPfile)
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
		else
			do i=1,nr
				T0(i)=10d0**(log10(TP0)+dTP*log10(P0(i)))
			enddo
		endif
	endif

	do i=1,nr
		T(i)=T0(nr+1-i)
		P(i)=P0(nr+1-i)
		if(mixratfile) then
			do j=1,n
				mixrat_r(i,imol(j))=mr0(nr+1-i,j)
			enddo
		endif
		if(dochemistry) then
			do j=1,nmol
				call MorleyChemistry(mixrat_r(i,j),T(i),P(i),molname(j),metallicity)
			enddo
		endif
	enddo

	do i=1,nclouds
		call output("==================================================================")
		call output("Setting up cloud: " // trim(int2string(i,'(i4)')))
		call SetupPartCloud(i)
		allocate(Cloud(i)%w(Cloud(i)%nsize))
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

	nclouds=0
	ncia=0
	nr=nPom*nTom
	maxtau=-1d0

	compute_opac=.true.
	
	allocate(dens(nr))
	allocate(Ndens(nr))
	allocate(R(nr+1))
	allocate(T(nr))
	allocate(P(nr))
	allocate(mixrat_r(nr,nmol))
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
	IMPLICIT NONE
	integer i
	character*100 homedir
	
	Mplanet=1d0
	Rplanet=1d0
	
	Tstar=5777d0
	Rstar=1d0
	Dplanet=1d0
	logg=4.5d0
	
	lam1=1d0
	lam2=15d0
	specres=10d0

	mixrat=0d0
	includemol=.false.
	
	obs%nphase=45

	dochemistry=.false.
	metallicity=1d0

	do i=1,nclouds
		Cloud(i)%P=1d-4
		Cloud(i)%dP=10d0
		Cloud(i)%s=2d0
		Cloud(i)%column=0d0
		Cloud(i)%file=' '
		Cloud(i)%standard='ASTROSIL'
		Cloud(i)%nsize=50
		Cloud(i)%nsubgrains=1
		Cloud(i)%amin=0.1
		Cloud(i)%amax=100d0
		Cloud(i)%fmax=0d0
		Cloud(i)%porosity=0d0
		Cloud(i)%fcarbon=0d0
		Cloud(i)%blend=.true.
		Cloud(i)%reff=1d0
		Cloud(i)%veff=0.1
		Cloud(i)%coverage=1.0
		Cloud(i)%ptype='COMPUTE'
	enddo
	nspike=0

	retrieval=.false.
	computeT=.false.
	TeffP=600d0
	outputopacity=.false.
	nr=20
	
	Pmin=1d-6
	Pmax=1d+3

	call getenv('HOME',homedir) 

	HITRANdir=trim(homedir) // '/HITRAN/'
	HITEMPdir=trim(homedir) // '/HITEMP/'

	HITEMP=.false.

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
	opacitydir='Opacities'

	nTom=100
	nPom=50
	Tmin=71d0
	Tmax=2900d0
	
	maxiter=6
	
	return
	end

	
	subroutine InitObs()
	use GlobalSetup
	use Constants
	IMPLICIT NONE

c number of cloud/nocloud combinations
	obs%ncc=2**nclouds
	allocate(obs%docloud(obs%ncc,nclouds))
	allocate(obs%cloudfrac(obs%ncc))
	allocate(obs%lam(nlam))
	allocate(obs%flux(0:obs%ncc,nlam))
	allocate(obs%A(0:obs%ncc,nlam))
	allocate(obs%phase(obs%nphase,0:obs%ncc,nlam))
	allocate(obs%dflux(nretr,nlam))
	allocate(obs%retr_par(nretr))
	allocate(obs%retr_dpar(nretr))
	allocate(obs%retr_ir(nretr))

	return
	end
	
	
	subroutine ReadObs(key)
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	
	select case(key%key2)
		case("nphase")
			read(key%value,*) obs%nphase
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
	real*8 lam0,T0,Planck
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
		do i=1,nlam-1
			BB(j,i)=Planck(T0,freq(i))
		enddo
	enddo

	return
	end


	subroutine ReadCloud(key)
	use GlobalSetup
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey) key
	integer j

	Cloud(key%nr1)%ptype="COMPUTE"

	select case(key%key2)
		case("file")
			Cloud(key%nr1)%file=trim(key%value)
		case("ngrains","nsize")
			read(key%value,*) Cloud(key%nr1)%nsize
		case("nsubgrains")
			read(key%value,*) Cloud(key%nr1)%nsubgrains
		case("amin")
			read(key%value,*) Cloud(key%nr1)%amin
		case("amax")
			read(key%value,*) Cloud(key%nr1)%amax
		case("fmax")
			read(key%value,*) Cloud(key%nr1)%fmax
		case("blend")
			read(key%value,*) Cloud(key%nr1)%blend
		case("porosity")
			read(key%value,*) Cloud(key%nr1)%porosity
		case("standard")
			Cloud(key%nr1)%standard=trim(key%value)
		case("fcarbon")
			read(key%value,*) Cloud(key%nr1)%fcarbon
		case("pressure","p")
			read(key%value,*) Cloud(key%nr1)%P
		case("dp")
			read(key%value,*) Cloud(key%nr1)%dP
		case("s","sharpness")
			read(key%value,*) Cloud(key%nr1)%s
		case("column")
			read(key%value,*) Cloud(key%nr1)%column
		case("reff")
			read(key%value,*) Cloud(key%nr1)%reff
		case("veff")
			read(key%value,*) Cloud(key%nr1)%veff
		case("coverage")
			read(key%value,*) Cloud(key%nr1)%coverage
		case default
			call output("Unknown cloud keyword: " // trim(key%key2))
			stop
	end select

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
	
	allocate(Cloud(ii)%rv(Cloud(ii)%nsize))
	allocate(Cloud(ii)%Kabs(Cloud(ii)%nsize,nlam))
	allocate(Cloud(ii)%Ksca(Cloud(ii)%nsize,nlam))
	allocate(Cloud(ii)%Kext(Cloud(ii)%nsize,nlam))
	allocate(Cloud(ii)%F(Cloud(ii)%nsize,nlam))

	select case(Cloud(ii)%ptype)
		case("COMPUTE")
			do is=1,Cloud(ii)%nsize
				call tellertje(is,Cloud(ii)%nsize)
				call ComputePart(Cloud(ii),ii,is)
			enddo
c		case("PARTFILE")
c			call ReadParticle(Cloud(ii),ii)
		case default
			call output("I did not understand what I was trying to do. Sorry!")
	end select
	
	if(nspike.gt.0.and.nspike.le.180) call output("making the first " // trim(int2string(nspike,'(i3)')) // " degrees isotropic")
	do ilam=1,nlam
		do is=1,Cloud(ii)%nsize
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
		do is=1,Cloud(ii)%nsize
			Cloud(ii)%F(is,ilam)%IF11=0d0
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
