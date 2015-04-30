	subroutine Init()
	use GlobalSetup
	use Constants
	use ReadKeywords
	IMPLICIT NONE
	type(SettingKey),pointer :: key,first

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
		case("comp","part","molecule")
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
		case default
			call output("Keyword not recognised: " // trim(key%key1))
			stop
	end select

	key => key%next
	
	enddo

	call output("==================================================================")

	
	return
	end
	

	subroutine SetDefaults()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i
	
	Mplanet=1d0
	Rplanet=1d0
	
	do i=1,nobs
		obs(i)%type='EMIS'
		obs(i)%filename=' '
	enddo

	retrieval=.false.
	nr=100
	
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

