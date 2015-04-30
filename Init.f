	subroutine Init()
	use GlobalSetup
	use Constants
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
	
	
	
	return
	end
	
	subroutine ReadObs(key)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(SettingKey) key
	integer i
	i=key%nr1
	
	select case(key%key2)
		case("type")
			read(key%value,'(a)') obs(i)%type
		case default
			call output("Keyword not recognised: " // trim(key%key2))
	end select
	
	return
	end
	

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
	call get_key_value(readline,key%key1,key%key2,key%value,key%nr1,key%nr2)

c read another command, so go back
	goto 10

30	continue
	close(unit=21)
	allocate(key%next)
	key => key%next
	key%last=.true.

	return
	end
	

c=========================================================================================
c This subroutine just seperates the key and value component of a string given
c key=value syntax. Key is transformed to lowercase.
c=========================================================================================
	subroutine get_key_value(line,key1,key2,value,nr1,nr2)
	IMPLICIT NONE
	character*1000 line
	character*100 key1,key2,value
	integer i,nr1,nr2,ikey1,ikey2
	
	ikey1=index(line,'=')
	ikey2=index(line,':')
	if(ikey1.eq.0) ikey1=len_trim(line)+1

	nr1=1
	nr2=1
	if(ikey2.gt.0) then
		key1=line(1:ikey2-1)
		key2=line(ikey2+1:ikey1-1)
		call checknr(key1,nr1)
		call checknr(key2,nr2)
	else
		key1=line(1:ikey1-1)
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
	end
	
	
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
	end
	



	subroutine CountStuff(firstkey)
	use GlobalSetup
	IMPLICIT NONE
	type(SettingKey),target :: firstkey
	type(SettingKey),pointer :: key
	integer i

	key => firstkey

	ncomp=1
	do while(.not.key%last)
		select case(key%key1)
			case("comp","part","molecule")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.ncomp) ncomp=key%nr1
		end select
		key=>key%next
	enddo

	call output('Number of components:   ' // int2string(ncomp,'(i4)'))

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


