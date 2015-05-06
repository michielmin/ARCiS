	subroutine WriteOutput()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs,i
	character*500 filename
	
	do iobs=1,nobs
		call output("==================================================================")
		filename=trim(outputdir) // "obs" // trim(int2string(iobs,'(i0.2)'))
		call output("Writing spectrum to: " // trim(filename))
		open(unit=30,file=filename,RECL=1000)
		write(30,'("#",a13,a19)') "lambda [mu]","flux [Jy]"
		do i=1,nlam-1
			write(30,'(f12.6,e19.7)') sqrt(lam(i)*lam(i+1))/micron,obs(iobs)%flux(i)
		enddo
		close(unit=30)
	enddo
	
	return
	end
	
	
