	program main
	use GlobalSetup
	IMPLICIT NONE
	character*500 VersionGIT
	integer i
	real*8 starttime,stoptime

	call SetOutputMode(.true.)
	
	call cpu_time(starttime)

	call GetOutputDir
	open(unit=9,file=trim(outputdir) // "log.dat",RECL=6000)
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

	call cpu_time(stoptime)
	call output("Initialisation time: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	call output("==================================================================")

	if(opacitymode) then
		if(dochemistry) then
			call OnlyChemCompute
		else
			do i=1,nmol
				mixrat_r(:,i)=mixrat(i)
			enddo
		endif
		call SetupOpacities()
		call WriteOpacityFITS()
	else if(dopostequalweights) then
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

	call cpu_time(stoptime)
	call output("==================================================================")
	call output("Total runtime:       " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	call output("All done!")
	call output("==================================================================")

	end
	
	subroutine ComputeModel(recomputeopacities)
	use GlobalSetup
	IMPLICIT NONE
	logical recomputeopacities
	
	if(do3D) then
		call Run3D(recomputeopacities)
	else
		call ComputeModel1D(recomputeopacities)
	endif
	
	return
	end


	subroutine ComputeModel1D(recomputeopacities)
	use GlobalSetup
	use TimingModule
	IMPLICIT NONE
	logical Tconverged
	real*8 starttime,stoptime,starttime_w,stoptime_w,omp_get_wtime,f
	logical recomputeopacities
	logical computeopac,temp
	integer i
	
	computeopac=recomputeopacities

	call cpu_time(starttime)
	Tconverged=.false.
	nTiter=0
	if(computeT.and.computeopac) then
		temp=par_tprofile
		par_tprofile=.false.
		f=1d0
		computelam=RTgridpoint
		do nTiter=1,maxiter
			call output("Temperature computation (" // trim(int2string(nTiter,'(i3)')) // " of " 
     &					// trim(int2string(maxiter,'(i3)')) // ")")
			call SetupStructure(.true.)
			call SetupOpacities()
			if(nTiter.eq.1) then
				f=1d0
			else
				f=0.5d0
				if(forceEbalance) f=f+0.5d0*exp(-real(maxiter-nTiter)/5d0)
			endif
			if(f.gt.1d0) f=1d0
			call DoComputeT(Tconverged,f)
			if(Tconverged.and.nTiter.gt.4) exit
c			call SetoutputMode(.true.)
c			call output("Chemistry cpu time: " // trim(dbl2string(timechem,'(f10.4)')) // " s")
c			call output("Chemistry walltime: " // trim(dbl2string(dble(itimechem)/dble(rate),'(f10.4)')) // " s")
c			call output("Number of chemistry calls: " // trim(int2string(ctimechem,'(i5)')))
c			call output("Cloud cpu time:    " // trim(dbl2string(timecloud,'(f10.4)')) // " s")
c			call output("Cloud walltime:    " // trim(dbl2string(dble(itimecloud)/dble(rate),'(f10.4)')) // " s")
c			call output("Number of cloud calls:     " // trim(int2string(ctimecloud,'(i5)')))
c			call output("PTstruct cpu time: " // trim(dbl2string(timetemp,'(f10.4)')) // " s")
c			call output("PTstruct walltime: " // trim(dbl2string(dble(itimetemp)/dble(rate),'(f10.4)')) // " s")
c			call output("Number of PTstruct calls:  " // trim(int2string(ctimetemp,'(i5)')))
c			call SetoutputMode(.false.)
		enddo
		computelam=.not.RTgridpoint
		if(forceEbalance) computelam=.true.
		call SetupStructure(.true.)
		call SetupOpacities()
		if(forceEbalance) then
			f=1d0
			nTiter=1
			call DoComputeT(Tconverged,f)
			computelam=.not.RTgridpoint
		endif
		par_tprofile=temp
	else
		call SetupStructure(computeopac)
		if(domakeai.and..not.modelsucces) return
		if(computeopac) call SetupOpacities()
	endif
	call cpu_time(stoptime)
	call output("Opacity computation: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	if(.not.do3D) call Raytrace()
	call cpu_time(stoptime)
	call output("Model runtime:       " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")

	return
	end
