#ifdef USE_MULTINEST
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
	integer i1,i2

	return

	if(n2d.eq.0) then
		i1=0
		i2=0
	else
		i1=1
		i2=n2d
	endif
	do i2d=i1,i2
		call WritePTlimitsMN
	enddo

	return
	end
	
	
	subroutine doMultiNest
	use Nested
	use params_multinest
	use GlobalSetup
	use RetrievalMod
	IMPLICIT NONE
	external getloglike,dumper
	integer i,context
	integer nest_pWrap(n_ret)

	imodel=0
	bestlike=-1d200
	write(nest_root,'(a,"/")') trim(outputdir)		
	nest_pWrap=0
	sdim=n_ret
	nest_nClsPar=n_ret
	nest_nlive=npop
	nest_resume=resume_multinest
	nest_efr=f_multinest
	nest_tol=tol_multinest
	nest_ceff=const_eff_multinest

	if(nest_resume) then
		open(unit=31,file=trim(outputdir) // '/Wolk.dat',RECL=6000,ACCESS='APPEND')
	else
		open(unit=31,file=trim(outputdir) // '/Wolk.dat',RECL=6000)
	endif

	call nestRun(nest_IS,nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_efr,sdim,sdim, 
     & nest_nClsPar,nest_maxModes,nest_updInt,nest_Ztol,nest_root,nest_rseed,nest_pWrap, 
     & nest_fb,nest_resume,nest_outfile,nest_initMPI,nest_logZero,nest_maxIter,getloglike,dumper,context)

	return
	end

#else
	subroutine doMultiNest
	use GlobalSetup
	IMPLICIT NONE
	call output("ARCiS was compiled without MultiNest")
	stop
	return
	end
#endif


	
	subroutine slikelihood(var,nvars,lnew)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nvars,i,j,nlamtot,k,ny,ii
	real*8 var(nvars),chi2obs(nobs),error(2,nvars),lnew,scale
	real*8,allocatable :: spec(:)
	logical recomputeopac
	real*8 tot
	character*500 keyword

	k=0
	do i=1,nobs
		do j=1,ObsSpec(i)%ndata
			k=k+1
		enddo
	enddo
	ny=k
	allocate(spec(ny))
	call mrqcomputeY(var,spec,nvars,ny,lnew,scale)

	do i=1,nvars
		var(i)=RetPar(i)%value
		if(RetPar(i)%logscale) var(i)=log10(var(i))
	enddo
	deallocate(spec)

	return
	end
	
	
