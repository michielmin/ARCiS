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
		write(30,'("#",a13,a19,a19)') "lambda [mu]","flux [Jy]","R^2/Rp^2"
		do i=1,nlam-1
			write(30,'(f12.6,es19.7,f19.9)') sqrt(lam(i)*lam(i+1))/micron,obs(iobs)%flux(i),
     &										obs(iobs)%A(i)/(pi*Rstar**2)
		enddo
		close(unit=30)
	enddo
	
	return
	end
	
	


	subroutine WriteOpacity(ir,flag,nu0,kappa0,nnu0,ng0)
	use GlobalSetup
	IMPLICIT NONE
	character*500 file
	character*4 flag
	integer nnu0,i,ir,ng0,j
	real*8 nu0(nnu0),kappa0(nnu0,ng0)
	
	file=trim(outputdir) // "opacity_" // trim(flag) // "_" // trim(int2string(ir,'(i0.4)')) // ".dat"
	open(unit=30,file=file,RECL=100)
	write(30,'("# Pressure:    ",es10.3," Ba")') P(ir)
	write(30,'("# Temperature: ",f10.3," K")') T(ir)
	write(30,'("#",a13,a19)') "lambda [mu]","kappa [cm^2/mol]"
	do i=1,nnu0
		do j=1,ng0
			write(30,'(f12.6,es19.7)') 1d4/nu0(i),kappa0(i,j)
		enddo
	enddo
	close(unit=30)
	
	return
	end
