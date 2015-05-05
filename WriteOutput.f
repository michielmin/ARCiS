	subroutine WriteOutput()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs,ilam
	
	do iobs=1,nobs
		open(unit=20,file=trim(outputdir) // "obs" // trim(int2string(iobs,'(i0.2)')),RECL=1000)
		do ilam=1,nlam-1
			write(20,*) obs(iobs)%lam(ilam)/micron,obs(iobs)%flux(ilam)
		enddo
		close(unit=20)
	enddo
	
	return
	end
	
	
