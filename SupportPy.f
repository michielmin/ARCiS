	subroutine pyInitSupport(pyinputfile,pyoutputdir,cline,py_nlam,py_nr,py_n_ret,py_nobs)
	use GlobalSetup
	use RetrievalMod
	IMPLICIT NONE
	character*500 VersionGIT
	character(*),optional :: cline
	character(*) pyinputfile,pyoutputdir
	character*1000 dummy
	integer i,j,k,py_nlam,py_nr,py_n_ret,py_nobs

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
! terms of use
	call output("By using ARCiS you agree to the terms of use.")
	call output("It basically means you consult with any of the authors before")
	call output("publication if there is any doubt on the outcome of the results.")

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

	call output("==================================================================")

	call SetOutputMode(.false.)

	py_nlam=0
	do i=1,nlam
		if(.not.RTgridpoint(i)) py_nlam=py_nlam+1
	enddo
	py_nr=nr
	py_n_ret=n_ret
	py_nobs=nobs

	return
	end


	subroutine pyRunARCiSSupport()
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
	


	subroutine GetAndSetKeyPy(readline)
	use GlobalSetup
	use ReadKeywords
	use Constants
	IMPLICIT NONE
	character*1000 readline
	logical hasnr1,hasnr2
	integer i
	type(SettingKey) key

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
	
	call get_key_value(readline,key%key,key%key1,key%key2,key%orkey1,key%orkey2,key%value,key%nr1,key%nr2,hasnr1,hasnr2,key%key2d)
	call ReadAndSetKey(key)

	call ConvertUnits()
	metallicity0=metallicity
	
	return
	end


	subroutine pyWriteBestfitSupport()
	use GlobalSetup
	IMPLICIT NONE
	integer i,status

	status=system("cp " // trim(outputdir) // "input.dat " // trim(outputdir) // "bestfit.dat")
	open(unit=21,file=trim(outputdir) // "bestfit.dat",FORM="FORMATTED",access='APPEND')
	write(21,'("*** retrieval keywords ***")')
	write(21,'("retrieval=.false.")')
	do i=1,n_ret
		write(21,'(a," = ",es14.7)') trim(RetPar(i)%keyword),RetPar(i)%value
	enddo
	close(unit=21)	

	return
	end
	


	function pyTransformRetParSupport(var,i) result(x)
	use GlobalSetup
	IMPLICIT NONE
	real*8 var,x
	integer i,j

	if(i.gt.1.and.RetPar(i)%increase) then
		RetPar(i)%xmin=RetPar(i-1)%value
	endif
	if(RetPar(i)%logscale) then
!	log
		x=var
		if(RetPar(i)%keyword(1:6).eq."Ppoint") then
			read(RetPar(i)%keyword(7:len(RetPar(i)%keyword)),*) j
			x=1d0-x**(1d0/(real(nTpoints-j+1)))
		endif
		x=(log10(RetPar(i)%xmin)+log10(RetPar(i)%xmax/RetPar(i)%xmin)*x)
	else if(RetPar(i)%squarescale) then
!	square
		x=var
		x=sqrt(RetPar(i)%xmin**2+(RetPar(i)%xmax**2-RetPar(i)%xmin**2)*x)
	else
!	linear
		x=var
		x=RetPar(i)%xmin+(RetPar(i)%xmax-RetPar(i)%xmin)*x
	endif

	return
	end

	subroutine pyGetObsSupport(i,l,obs,dobs,model)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 l(:),obs(:),model(:),dobs(:)
	integer i,n
	n=ObsSpec(i)%ndata
	
	l(1:n)=ObsSpec(i)%lam(1:n)
	obs(1:n)=ObsSpec(i)%y(1:n)
	dobs(1:n)=ObsSpec(i)%dy(1:n)
	model(1:n)=ObsSpec(i)%model(1:n)

	return
	end


	function pyGetRetrievalNamesSupport(n) result(names)
	use GlobalSetup
	IMPLICIT NONE
	integer,intent(in) :: n
	character*500 :: names
	
	names=RetPar(n)%keyword

	return
	end

	function pyGetObsNSupport(i) result(n)
	use GlobalSetup
	integer n
	n=ObsSpec(i)%ndata
	return
	end



	subroutine pyGetLamSupport(l,py_nlam)
	use GlobalSetup
	IMPLICIT NONE
	integer py_nlam
	real*8 l(py_nlam)
	integer i,j
	
	j=0
	do i=1,nlam
		if(.not.RTgridpoint(i)) then
			j=j+1
			l(j)=lam(i)
		endif
	enddo

	return
	end

	subroutine pyGetTransSupport(l,trans,py_nlam)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer py_nlam
	real*8 l(py_nlam),trans(py_nlam)
	integer i,j
	
	j=0
	do i=1,nlam
		if(.not.RTgridpoint(i)) then
			j=j+1
			l(j)=lam(i)
			trans(j)=obsA(0,i)/(pi*Rstar**2)
		endif
	enddo

	return
	end

	subroutine pyGetEmisSupport(l,emis,ip,py_nlam)
	use GlobalSetup
	IMPLICIT NONE
	integer py_nlam
	real*8 l(py_nlam),emis(py_nlam)
	integer i,ip,j
	
	j=0
	do i=1,nlam
		if(.not.RTgridpoint(i)) then
			j=j+1
			l(j)=lam(i)
			emis(j)=phase(1,ip,i)+flux(ip,i)
		endif
	enddo

	return
	end

	subroutine pyGetStarSupport(l,star,py_nlam)
	use GlobalSetup
	IMPLICIT NONE
	integer py_nlam
	real*8 l(py_nlam),star(py_nlam)
	integer i,j
	
	j=0
	do i=1,nlam
		if(.not.RTgridpoint(i)) then
			j=j+1
			l(j)=lam(i)
			star(j)=Fstar(i)
		endif
	enddo

	return
	end


	subroutine pyGetPTSupport(Ppy,Tpy)
	use GlobalSetup
	IMPLICIT NONE
	real*8 Ppy(nr),Tpy(nr)

	Ppy(1:nr)=P(1:nr)
	Tpy(1:nr)=T(1:nr)

	return
	end

	subroutine pyVerboseSupport(switch)
	use GlobalSetup
	IMPLICIT NONE
	logical switch

	writefiles=switch
	call SetOutputMode(.false.)

	return
	end
	
	
	
