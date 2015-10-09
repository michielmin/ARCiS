	subroutine WriteOutput()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i
	character*500 filename
	character*6000 form
	real*8,allocatable :: theta(:)
	logical,allocatable :: docloud(:,:)

	allocate(docloud(nclouds,obs%ncc))
	allocate(theta(obs%nphase))
	do i=1,obs%nphase
		theta(i)=acos(1d0-2d0*(real(i)-0.5d0)/real(obs%nphase))*180d0/pi
	enddo
	do i=1,nclouds
		docloud(i,:)=obs%docloud(:,i)
	enddo
	call output("==================================================================")
	filename=trim(outputdir) // "emis"
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,RECL=1000)
	if(nclouds.gt.0) then
		form='("#",a13,a19,' // trim(int2string(obs%ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]","flux [Jy]",docloud(1:nclouds,1:obs%ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","flux [Jy]"
	endif
	form='(f14.6,' // int2string(obs%ncc+1,'(i3)') // 'es19.7)'
	do i=1,nlam-1
		write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					obs%flux(0:obs%ncc,i)
	enddo
	close(unit=30)

	filename=trim(outputdir) // "emisR"
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,RECL=1000)
	if(nclouds.gt.0) then
		form='("#",a13,a19,' // trim(int2string(obs%ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]","flux [Jy]",docloud(1:nclouds,1:obs%ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","flux [Jy]"
	endif
	form='(f14.6,' // int2string(obs%ncc+1,'(i3)') // 'es19.7)'
	do i=1,nlam-1
		write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					obs%flux(0:obs%ncc,i)/(Fstar(i)*1d23/distance**2)
	enddo
	close(unit=30)


	filename=trim(outputdir) // "trans"
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,RECL=1000)
	if(nclouds.gt.0) then
		form='("#",a13,a19,' // trim(int2string(obs%ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]","Rp^2/Rstar^2",docloud(1:nclouds,1:obs%ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","Rp^2/Rstar^2"
	endif
	form='(f14.6,' // int2string(obs%ncc+1,'(i3)') // 'es19.7)'
	do i=1,nlam-1
		write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					obs%A(0:obs%ncc,i)/(pi*Rstar**2)
	enddo
	close(unit=30)



	filename=trim(outputdir) // "phase"
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,RECL=6000)
	form='("#",a13,' // trim(int2string(obs%nphase,'(i3)')) // 
     &				 '("   flux(",f5.1,") [Jy]"),"         fstar [Jy]")'
	write(30,form) "lambda [mu]",theta(1:obs%nphase)
	form='(f14.6,' // int2string(obs%nphase+1,'(i3)') // 'es19.7)'
	do i=1,nlam-1
		write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					obs%phase(1:obs%nphase,0,i)+obs%flux(0,i),
     &					Fstar(i)*1d23/distance**2
	enddo
	close(unit=30)
	deallocate(docloud)
	
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
