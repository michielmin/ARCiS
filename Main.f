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
	IMPLICIT NONE
	logical Tconverged
	real*8 starttime,stoptime,starttime_w,stoptime_w,omp_get_wtime,f
	logical recomputeopacities
	logical computeopac,temp
	real*8 srdtemp,srd,ldtemp(nlamdust),lam0
	integer nldtemp,nld,i
	
	computeopac=recomputeopacities

	call cpu_time(starttime)
	Tconverged=.false.
	nTiter=0
	fiter=1d0
	if(computeT.and.computeopac) then
		nldtemp=nlamdust
		srdtemp=specresdust
		ldtemp(1:nlamdust)=lamdust(1:nlamdust)
		specresdust=10d0
		if(useobsgrid) then
			nlamdust=nlam
			lamdust=lam
		else
			lam0=lam1
			nld=1
			do while(lam0.le.lam2)
				lam0=lam0+lam0/specresdust
				nld=nld+1
			enddo
			if(nld.lt.nlamdust) then
				nlamdust=nld
				i=1
				lamdust(i)=lam1
				do while(lamdust(i).le.lam2)
					i=i+1
					lamdust(i)=lamdust(i-1)+lamdust(i-1)/specresdust
				enddo
				lamdust(nlamdust)=lam2
			else
				nlamdust=nldtemp
				specresdust=srdtemp
				lamdust(1:nlamdust)=ldtemp(1:nlamdust)
			endif
		endif
	endif
	call SetupStructure(computeopac)
	if(domakeai.and..not.modelsucces) return
	if(computeopac) call SetupOpacities()
	if(computeT.and.computeopac) then
		temp=par_tprofile
		do while(.not.Tconverged.and.nTiter.le.maxiter)
			call output("Temperature computation (" // trim(int2string(nTiter,'(i3)')) // " of " 
     &					// trim(int2string(maxiter,'(i3)')) // ")")
			f=max(min(1d0/real(nTiter+1)+0.2d0,1d0),1d-2)
c			f=0.9
c			if(nTiter.ge.4) f=1d0/real(nTiter-2)**0.5+0.1
			fiter=f
			call DoComputeT(Tconverged,f)
			par_tprofile=.false.
			nTiter=nTiter+1
			if(Tconverged.or.nTiter.gt.maxiter) then
				nlamdust=nldtemp
				specresdust=srdtemp
				lamdust(1:nlamdust)=ldtemp(1:nlamdust)
			endif
			call SetupStructure(.true.)
			call SetupOpacities()
		enddo
		par_tprofile=temp
	endif
	call cpu_time(stoptime)
	call output("Opacity computation: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	if(.not.do3D) call Raytrace()
	call cpu_time(stoptime)
	call output("Model runtime:       " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")

	return
	end
