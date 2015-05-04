c=========================================================================================
c module containing the physical constants in cgs (moet nog naar SI)
c=========================================================================================
	module Constants
	IMPLICIT NONE
	real*8 pi,Ggrav,Msun,AU,clight,Rsun,mp,kb,hplanck,parsec,Lsun,sigma
	real*8 Mearth,Rearth,Mjup,Rjup,year,micron,mu
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
	parameter(mu=2.3) 		!mean molecular weight
	parameter(Ggrav=6.67300d-8) ! in cm^3/g/s^2
	parameter(hplanck=6.626068d-27) ! cm^2 g/s
	parameter(Mearth=5.97219d27)
	parameter(Mjup=1.89813d30)
	parameter(Rearth=6.3781d8)
	parameter(Rjup=7.1492d9)
	parameter(year=24d0*60d0*60d0*265.25d0)
	parameter(micron=1d-4)
	
	end module Constants

c=========================================================================================
c global setup for the code
c=========================================================================================
	module GlobalSetup
	IMPLICIT NONE
	real*8 Mplanet,Rplanet									! mass and radius of the planet
	real*8,allocatable :: dens(:),T(:),P(:),Ndens(:)		! radius
	real*8,allocatable :: dust_dens(:,:)					! radius, component
	real*8,allocatable :: R(:)								! radius
	real*8,allocatable :: mixrat(:)							! component
	real*8,allocatable :: opac(:,:,:,:)						! component,wav,T,P
	integer nT,np,nrad,nmol,nlam,nobs		! #T, #P, #radial points, #molecules, #wavelength bins, #obs
	integer nlines
	character*500 outputdir,HITRANfile
	integer nr,idum
!$OMP THREADPRIVATE(idum)
	logical retrieval
	real*8 lam1,lam2,specres,Pmin,Pmax
	real*8,allocatable :: lam(:),freq(:)
	real*8,allocatable :: ZZ(:,:,:),TZ(:)	! partition function
	integer nTZ
	integer,allocatable :: niso(:)
	real*8,allocatable :: Mmol(:)
	character*10 molname(47)
	parameter(molname = (/ 'H2O','CO2','O3','N2O','CO','CH4','O2','CO2','SO2','NO2',
     &	'NH3','HNO3','OH','HF','HCl','HBr','HI','ClO','OCS','H2CO','HOCl','N2',
     &	'HCN','CH3Cl','H2O2','C2H2','C2H6','PH3','COF2','SF6','H2S','HCOOH','HO2',
     &	'O','ClONO2','NO+','HOBr','C2H4','CH3OH','CH3Br','CH3CN','CF4','C4H2',
     &	'HC3N','H2','CS','SO3' /))


	type Observation
		character*500 filename
		character*10 type
		real*8,allocatable :: lam(:),flux(:)
	end type Observation

	type(Observation),allocatable :: obs(:)

	type Line
		integer imol,iiso
		real*8 Aul,freq,Eup,lam,S0,S
		real*8 gamma_air,gamma_self
		real*8 a_therm,a_press
		real*8 gu,gl
	end type Line

	type(Line),allocatable,target :: Lines(:)


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
	integer i,j

	key => firstkey

	nmol=1
	nobs=1
	j=0
	do while(.not.key%last)
		select case(key%key1)
			case("obs")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nobs) nobs=key%nr1
			case default
				do i=1,47
					if(key%key.eq.molname(i)) then
						if(i.gt.nmol) nmol=i
						j=j+1
					endif
				enddo
		end select
		key=>key%next
	enddo

	call output('Number of molecules:    ' // int2string(j,'(i4)'))
	call output('Number of observations: ' // int2string(nobs,'(i4)'))

	allocate(obs(nobs))
	allocate(Mmol(nmol))
	allocate(mixrat(nmol))

	return
	end subroutine CountStuff
	
	end module ReadKeywords
	
c==============================================================================	
c==============================================================================	
