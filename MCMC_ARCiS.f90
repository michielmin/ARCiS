#ifdef USE_MCMCF90
subroutine doMCMCF90(var0,nvar)
	use GlobalSetup
	IMPLICIT NONE
	integer nvar,i,j
	real*8 var0(nvar),cov(nvar,nvar),w,error(2,nvar)

	open(unit=20,file=trim(outputdir) // 'mcmcpar.dat',RECL=6000)
	write(20,*) var0(1:nvar)
	close(unit=20)
	
	open(unit=20,file=trim(outputdir) // 'mcmccov.dat',RECL=6000)
	cov=0d0
	do i=1,nvar
		cov(i,i)=1d-2
	enddo
	do i=1,nvar
		write(20,*) cov(1:nvar,i)
	enddo
	close(unit=20)

	open(unit=20,file='mcmcinit.nml',RECL=1000)
	write(20,'("!!")')
	write(20,'("!! Run time parameters for the mcmc run")')
	write(20,'("!!")')
	write(20,'("&mcmc")')
	write(20,'("method = ",a)') "'dram'"
	write(20,'(" ",a12,a2,i9)')" nsimu","=",npop
	write(20,'(" ",a12,a2,i9)')" verbosity","=",0
	write(20,'(" ",a12,a2,i9)')" doadapt","=",0
	write(20,'(" ",a12,a2,i9)')" adaptint","=",100
	write(20,'(" ",a12,a2,i9)')" burnintime","=",min(1000,npop/10)
	write(20,'(" ",a12,a2,i9)')" doburnin","=",1
	write(20,'(" ",a12,a2,i9)')" drscale","=",0
	write(20,'(" ",a12,a2,i9)')" printint","=",100
	write(20,'(" ",a12,a2,i9)')" updatesigma","=",0
	write(20,'(" ",a12,a2,i9)')" N0","=",0
	write(20,'(" ",a12,a2,i9)')" S02","=",0
	write(20,'(" ",a12,a3,a,a)')" chainfile","= '",trim(outputdir) // 'chain.dat',"'"
	write(20,'(" ",a12,a3,a,a)')" ssfile","= '",trim(outputdir) // 'sschain.dat',"'"
	write(20,'(" ",a12,a3,a,a)')" parfile","= '",trim(outputdir) // 'mcmcpar.dat',"'"
	write(20,'(" ",a12,a3,a,a)')" cov0file","= '",trim(outputdir) // 'mcmccov.dat',"'"
	write(20,'("/")')
	close(unit=20)

	call mcmc_main()
	
	call SetOutputMode(.false.)
	open(unit=20,file=trim(outputdir) // 'chain.dat',RECL=6000)
	open(unit=21,file=trim(outputdir) // 'posterior.dat',RECL=6000)
1	read(20,*,end=2) var0(1:nvar),w
	call MapRetrieval(var0,error)
	do i=1,nvar
		var0(i)=RetPar(i)%value
		if(RetPar(i)%logscale) var0(i)=log10(var0(i))
	enddo
	write(21,*) var0(1:nvar),w
	goto 1
2	close(unit=20)
	close(unit=21)
	call SetOutputMode(.true.)

	return
end subroutine doMCMCF90


function ssfunction(theta,npar,ny) result(ss)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
  integer :: npar, ny
  real*8 :: theta(npar)
  real*8 :: ss(ny)
	integer k,i
	real*8 amoebafunk,chi2
	
	k=0
	do i=1,nobs
		k=k+ObsSpec(i)%ndata
	enddo
	if(massprior) k=k+1

	chi2=amoebafunk(theta,k)

	ss(1)=global_chi2*max(1,k-npar)

end function ssfunction

!!!
!!! this function returns false if any theta(i) is out of bounds
!!!
function checkbounds(theta)
  implicit none
  real*8 theta(:)
  logical checkbounds
  
!! example: all thetas must be positive
  checkbounds = .true.
  if (any(theta<=0.0).or.any(theta>=1.0)) checkbounds = .false.
  return

end function checkbounds
#else
	subroutine doMCMCF90(var0,nvar)
	use GlobalSetup
	IMPLICIT NONE
	integer nvar
	real*8 var0(nvar)
	call output("ARCiS was compiled without mcmcf90")
	stop
	return
	end
#endif

