	subroutine WriteOutput()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs,i
	character*500 filename,form
	logical,allocatable :: docloud(:,:)

	do iobs=1,nobs
		allocate(docloud(nclouds,obs(iobs)%ncc))
		do i=1,nclouds
			docloud(i,:)=obs(iobs)%docloud(:,i)
		enddo
		call output("==================================================================")
		select case(obs(iobs)%type)
			case("EMIS","emis","emission","EMISSION")
				filename=trim(outputdir) // "emis" // trim(int2string(iobs,'(i0.2)'))
				call output("Writing spectrum to: " // trim(filename))
				open(unit=30,file=filename,RECL=1000)
				form='("#",a13,a19,' // trim(int2string(obs(iobs)%ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
				write(30,form) "lambda [mu]","flux [Jy]",docloud(1:nclouds,1:obs(iobs)%ncc)
				form='(f14.6,' // int2string(obs(iobs)%ncc+1,'(i3)') // 'es19.7)'
				do i=1,nlam-1
					write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					obs(iobs)%flux(0:obs(iobs)%ncc,i)
				enddo
				close(unit=30)
			case("TRANS","trans","transit","TRANSIT")
				filename=trim(outputdir) // "trans" // trim(int2string(iobs,'(i0.2)'))
				call output("Writing spectrum to: " // trim(filename))
				open(unit=30,file=filename,RECL=1000)
				form='("#",a13,a19,' // trim(int2string(obs(iobs)%ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
				write(30,form) "lambda [mu]","Rp^2/Rstar^2",docloud(1:nclouds,1:obs(iobs)%ncc)
				form='(f14.6,' // int2string(obs(iobs)%ncc+1,'(i3)') // 'es19.7)'
				do i=1,nlam-1
					write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					obs(iobs)%A(0:obs(iobs)%ncc,i)/(pi*Rstar**2)
				enddo
				close(unit=30)
		end select
		deallocate(docloud)
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
