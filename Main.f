	program main
	use GlobalSetup
	IMPLICIT NONE
	character*500 VersionGIT
	logical converged
	integer i

	call GetOutputDir
	open(unit=9,file=trim(outputdir) // "log.dat",RECL=6000)
	call output("Output dir: " // trim(outputdir))

c terms of use
	call output("==================================================================")
	call output("By using ELMO you agree to the terms of use.")
	call output("It basically means you offer us co-author rights on any paper.")
	call output("that uses results computed with ELMO.")

	call output("==================================================================")
	call output("Let's get the show on the road!!")
	call output("ELMO "//trim(VersionGIT()))
	call output("==================================================================")

	call Init()

	if(retrieval) call ReadObs()

	converged=.false.
	do while(.not.converged)
		call SetupStructure()
		call SetupOpacities()
		do i=1,nobs
			call Raytrace(i)
		enddo
		if(retrieval) then
			call AdjustParameters(converged)
		else
			converged=.true.
		endif
	enddo

	call WriteOutput()

	call output("==================================================================")
	call output("All done!")
	call output("==================================================================")

	end
	
