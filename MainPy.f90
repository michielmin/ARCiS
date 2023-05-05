	module pyEx
	IMPLICIT NONE
	real*8,allocatable :: pyP(:),pyT(:),pytrans(:),pyemis(:,:),pystar(:),pylam(:)
	real*8,allocatable :: pyobs(:,:)
	integer pynvars
	end module pyEx
	

	subroutine pysetoutputmode(tf)
	use GlobalSetup
	IMPLICIT NONE
	logical tf

	call SetOutputMode(tf)
	
	return
	end
	

	subroutine pyInit(pyinputfile,pyoutputdir,cline)
	use GlobalSetup
	use RetrievalMod
	use pyEx
	IMPLICIT NONE
	character*500 VersionGIT
	character(*),optional :: cline
	character(*) pyinputfile,pyoutputdir
	character*1000 dummy
	integer i,j,k

	call SetOutputMode(.true.)
	
	j=1
	i=0
	if(present(cline)) then
		do while(j.lt.len_trim(cline))
			read(cline(j:len_trim(cline)),*) dummy
			j=j+len_trim(dummy)+1
			i=i+1
		enddo
	endif

	ncommandargs=3+i
	
	allocate(commandargs(ncommandargs+1))
	write(commandargs(1),'(a)') pyinputfile
	write(commandargs(2),'(a)') '-o'
	write(commandargs(3),'(a)') pyoutputdir
	if(present(cline)) then
		j=1
		i=0
		do while(j.lt.len_trim(cline))
			read(cline(j:len_trim(cline)),*) commandargs(i+4)
			j=j+len_trim(commandargs(i+4))+1
			i=i+1
		enddo
	endif
	write(commandargs(ncommandargs+1),'(a)') ' '

	call GetOutputDir

	open(unit=9,file=trim(outputdir) // "log.dat",FORM="FORMATTED",ACCESS="STREAM")
	call output("Output dir: " // trim(outputdir))

	call output("==================================================================")
	call output("         ARtful modelling code for exoplanet Science - ARCiS")
	call output("==================================================================")
c terms of use
	call output("By using ARCiS you agree to the terms of use.")
	call output("It basically means you offer us co-author rights on any paper")
	call output("that uses results computed with ARCiS.")

	call output("==================================================================")
	call output("Let's get the show on the road!!")
	call output("ARCiS version "//trim(VersionGIT()))
	call output("==================================================================")

	call Init()

	if(nobs.ne.0) call ReadObs()

	if(n_ret.ne.0) then
		allocate(obsA0(nlam))
		allocate(obsA1(nlam))
		allocate(obsA2(nlam))
		allocate(dobsA(n_ret,nlam))
		allocate(emis0(nlam))
		allocate(emis1(nlam))
		allocate(emis2(nlam))
		allocate(demis(n_ret,nlam))
		allocate(emisR0(nlam))
		allocate(emisR1(nlam))
		allocate(emisR2(nlam))
		allocate(demisR(n_ret,nlam))
		allocate(bestvar(n_ret))
	endif
	imodel=0
	bestlike=-1d200

	allocate(pyP(nr),pyT(nr),pylam(nlam),pytrans(nlam),pyemis(nphase,nlam),pystar(1:nlam))
	if(nobs.gt.0) then
		j=0
		do i=1,nobs
			j=j+ObsSpec(i)%ndata
		enddo
		allocate(pyobs(4,j))
		k=0
		do i=1,nobs
			do j=1,ObsSpec(i)%ndata
				k=k+1
				pyobs(1,k)=ObsSpec(i)%lam(j)*1d4
				pyobs(2,k)=ObsSpec(i)%y(j)
				pyobs(3,k)=ObsSpec(i)%dy(j)
			enddo
		enddo
	endif

	call output("==================================================================")

	call SetOutputMode(.false.)

	return
	end

	subroutine pyRunARCiS()
	use GlobalSetup
	IMPLICIT NONE

	if(dopostequalweights) then
		call PostEqualWeights()
	else if(domakeai) then
		call MakeAI()
	else if(retrieval) then
		call DoRetrieval()
	else
		call ComputeModel(.true.)
		call WriteStructure()
		call WriteOutput()
	endif

	call output("==================================================================")
	call output("All done!")
	call output("==================================================================")

	return
	end
	
	subroutine pyComputeModel()
	use GlobalSetup
	use pyEx
	use Constants
	IMPLICIT NONE
	real*8,allocatable :: spec(:),spectemp(:)
	integer i,j,k

	call InitDens()
	call ComputeModel(.true.)
	
	pylam(1:nlam)=lam(1:nlam)*1d4
	pyP(1:nr)=P(1:nr)
	pyT(1:nr)=pyT(1:nr)
	pytrans(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
	pyemis(1:nphase,1:nlam)=phase(1,1:nphase,1:nlam)+flux(1:nphase,1:nlam)
	pystar(1:nlam)=Fstar(1:nlam)

	if(nobs.gt.0) then
		k=0
		do i=1,nobs
			if(ObsSpec(i)%ndata.gt.k) k=ObsSpec(i)%ndata
		enddo
		allocate(spec(k),spectemp(nlam))
		k=0
		do i=1,nobs
			call RemapObs(i,spec,spectemp)
			do j=1,ObsSpec(i)%ndata
				k=k+1
				pyobs(4,k)=spec(j)
			enddo
		enddo
		deallocate(spec)
	endif

	
	return
	end
	
	subroutine pyWriteFiles()
	use GlobalSetup
	IMPLICIT NONE

	call WriteStructure()
	call WriteOutput()

	return
	end
	
	
	subroutine pySetKeyword(pykey,pyval)
	use GlobalSetup
	use ReadKeywords
	use Constants
	IMPLICIT NONE
	type(SettingKey) key
	character(*) pykey,pyval
	character*1000 readline
	integer i

	Rplanet=Rplanet/Rjup
	Mplanet=Mplanet/Mjup
	Rstar=Rstar/Rsun
	Mstar=Mstar/Msun
	Dplanet=Dplanet/AU
	lam1=lam1/micron
	lam2=lam2/micron
	distance=distance/parsec
	do i=1,nclouds
		Cloud(i)%rnuc=Cloud(i)%rnuc/micron
	enddo
	orbit_inc=orbit_inc*180d0/pi

	metallicity=metallicity0

	readline=trim(pykey) // "=" // trim(pyval)
	call get_key_value(readline,key%key,key%key1,key%key2,key%orkey1,key%orkey2,key%value,key%nr1,key%nr2,key%key2d)
	call ReadAndSetKey(key)

	call ConvertUnits()
	metallicity0=metallicity

	return
	end

	
	subroutine pySetValue(pykey,pyval)
	use GlobalSetup
	use ReadKeywords
	use Constants
	IMPLICIT NONE
	type(SettingKey) key
	character(*) pykey
	real*8 pyval
	character*1000 readline
	integer i

	Rplanet=Rplanet/Rjup
	Mplanet=Mplanet/Mjup
	Rstar=Rstar/Rsun
	Mstar=Mstar/Msun
	Dplanet=Dplanet/AU
	lam1=lam1/micron
	lam2=lam2/micron
	distance=distance/parsec
	do i=1,nclouds
		Cloud(i)%rnuc=Cloud(i)%rnuc/micron
	enddo
	orbit_inc=orbit_inc*180d0/pi

	metallicity=metallicity0

	readline=trim(pykey) // "=" // trim(dbl2string(pyval,'(es14.7)'))
	call get_key_value(readline,key%key,key%key1,key%key2,key%orkey1,key%orkey2,key%value,key%nr1,key%nr2,key%key2d)
	call ReadAndSetKey(key)

	call ConvertUnits()
	metallicity0=metallicity

	return
	end

	function pyGetLike(var) result(lnew)
	use GlobalSetup
	IMPLICIT NONE
	real*8 var(:)
Cf2py indent(inout) var
	real*8 lnew

	call slikelihood(var,n_ret,lnew)
	if(.not.lnew.gt.-1d100) lnew=-1d100

	return
	end


	function pyGetRetrievalNames(n) result(names)
	use GlobalSetup
	IMPLICIT NONE
	integer,intent(in) :: n
	character*500 :: names
	
	names=RetPar(n)%keyword

	return
	end

