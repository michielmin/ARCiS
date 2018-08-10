	subroutine getloglike(Cube,ndims,nPar,lnew) 					!subroutine which gives lnew=loglike(Cube(ndims))
	IMPLICIT NONE								!distribution, parameter constraints, max loglike & log evidence values
	integer ndims,nPar
	real*8 Cube(ndims),lnew

	call slikelihood(Cube,ndims,lnew)

	return
	end

	subroutine dumper(nSamples,nlive,nPar,physLive,posterior, 
     &     		paramConstr,maxloglike,logZ,INSlogZ,logZerr,context)		!subroutine called after every updInt*10 iterations with the posterior 
	use GlobalSetup
	IMPLICIT NONE								!distribution, parameter constraints, max loglike & log evidence values
	integer nSamples,nlive,nPar
	real*8 physLive(nlive, nPar+1) 		     		!= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
	real*8 posterior(nSamples, nPar+2)		     		!= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
	real*8 paramConstr(1, 4*nPar)
	real*8 maxLogLike					    	!= maximum loglikelihood value
	real*8 logZ						     	!= log evidence value from the default (non-INS) mode
	real*8 INSlogZ						     	!= log evidence value from the INS mode
	real*8 logZerr						     	!= error on log evidence value
	integer context							!not required by MultiNest, any additional information user wants to pass

	return
	end
	
	
	subroutine slikelihood(var,nvars,lnew)
	use GlobalSetup
	use Constants
	use params_multinest
	IMPLICIT NONE
	integer nvars,i,j,nlamtot
	real*8 var(nvars),chi2obs(nobs),error(2,nvars),lnew
	real*8,allocatable :: spec(:)
	logical recomputeopac

	recomputeopac=.true.
	imodel=imodel+1
	call output("model number: " // int2string(imodel,'(i7)'))

	do i=1,nvars
		if(var(i).gt.1d0) var(i)=1d0
		if(var(i).lt.0d0) var(i)=0d0
	enddo
	error=0d0
	call MapRetrieval(var,error)

	do i=1,nvars
		var(i)=RetPar(i)%value
		if(RetPar(i)%logscale) var(i)=log10(var(i))
	enddo

	call InitDens()
	call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam,starfile)
	Fstar=Fstar*pi*Rstar**2
	call SetOutputMode(.false.)
	call ComputeModel(recomputeopac)
	call SetOutputMode(.true.)
	
	lnew=0d0
	nlamtot=0
	do i=1,nobs
		allocate(spec(ObsSpec(i)%nlam))
		call RemapObs(i,spec)
		chi2obs(i)=0d0
		do j=1,ObsSpec(i)%nlam
			chi2obs(i)=chi2obs(i)+((spec(j)-ObsSpec(i)%y(j))/ObsSpec(i)%dy(j))**2
		enddo
		nlamtot=nlamtot+ObsSpec(i)%nlam
		lnew=lnew+chi2obs(i)
		chi2obs(i)=chi2obs(i)/real(ObsSpec(i)%nlam)
		deallocate(spec)
	enddo
	lnew=lnew/real(nlamtot)

	write(72,*) imodel,lnew,var(1:nvars),COratio,metallicity
	call flush(72)

	if(lnew.lt.bestlike) then
		call WriteStructure()
		call WriteOutput()

		do i=1,nobs
			select case(ObsSpec(i)%type)
				case("trans","transmission","emisr","emisR","emisa","emis","emission","transC")
					open(unit=20,file=trim(outputdir) // "obs" // trim(int2string(i,'(i0.3)')),RECL=1000)
					do j=1,ObsSpec(i)%nlam
						write(20,*) ObsSpec(i)%lam(j)*1d4,ObsSpec(i)%model(j),ObsSpec(i)%y(j),ObsSpec(i)%dy(j)
					enddo
					close(unit=20)
			end select
		enddo

		call system("cp " // trim(outputdir) // "input.dat " // trim(outputdir) // "bestfit.dat")
		open(unit=21,file=trim(outputdir) // "bestfit.dat",RECL=1000,access='APPEND')
		write(21,'("*** retrieval keywords ***")')
		write(21,'("retrieval=.false.")')
		do i=1,n_ret
			write(21,'(a," = ",es14.7)') trim(RetPar(i)%keyword),RetPar(i)%value
		enddo
		close(unit=21)	

		bestlike=lnew
	endif
	lnew=-0.5d0*lnew
	
	return
	end
	
	
	
	subroutine doMultiNest
	use Nested
	use params_multinest
	use GlobalSetup
	IMPLICIT NONE
	external getloglike,dumper
	integer i,context
	integer nest_pWrap(n_ret)

	imodel=0
	bestlike=1d200
	write(nest_root,'(a,"/")') trim(outputdir)		
	nest_pWrap=0
	sdim=n_ret
	nest_nClsPar=n_ret
	nest_nlive=npop
	nest_resume=resume_multinest
	nest_efr=f_multinest
	nest_tol=tol_multinest

	if(nest_resume) then
		open(unit=72,file=trim(outputdir) // '/Wolk.dat',RECL=6000,ACCESS='APPEND')
	else
		open(unit=72,file=trim(outputdir) // '/Wolk.dat',RECL=6000)
	endif

	call nestRun(nest_IS,nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_efr,sdim,sdim, 
     & nest_nClsPar,nest_maxModes,nest_updInt,nest_Ztol,nest_root,nest_rseed,nest_pWrap, 
     & nest_fb,nest_resume,nest_outfile,nest_initMPI,nest_logZero,nest_maxIter,getloglike,dumper,context)

	return
	end

