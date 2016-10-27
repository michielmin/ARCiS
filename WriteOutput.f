	subroutine WriteOutput()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i
	character*500 filename
	character*6000 form
	real*8,allocatable :: theta(:)
	real*8 Fp1,Fp2,ApAs
	logical,allocatable :: docloud0(:,:)
	real*8,allocatable :: spec(:,:),specR(:),lamR(:),specRexp(:),specSNR(:)
	real*8 x,specres_obs,expspecres_obs,SNR,gasdev
	integer ilam,j,nj,nlamR
	character*1000 line

	allocate(docloud0(max(nclouds,1),ncc))
	allocate(theta(nphase))
	do i=1,nphase
		theta(i)=acos(1d0-2d0*(real(i)-0.5d0)/real(nphase))*180d0/pi
	enddo
	do i=1,nclouds
		docloud0(i,:)=docloud(:,i)
	enddo
	call output("==================================================================")
	filename=trim(outputdir) // "emis"
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,RECL=1000)
	if(nclouds.gt.0) then
		form='("#",a13,a19,' // trim(int2string(ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]","flux [Jy]",docloud0(1:nclouds,1:ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","flux [Jy]"
	endif
	form='(f14.6,' // int2string(ncc+1,'(i3)') // 'es19.7)'
	do i=1,nlam-1
		write(30,form) sqrt(lam(i)*lam(i+1))/micron,
c     &					flux(0:ncc,i)
     &					(phase(1,j,i)+flux(j,i),j=0,ncc)
	enddo
	close(unit=30)

	filename=trim(outputdir) // "emisR"
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,RECL=1000)
	if(nclouds.gt.0) then
		form='("#",a13,a19,' // trim(int2string(ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]","flux [Jy]",docloud0(1:nclouds,1:ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","flux [Jy]"
	endif
	form='(f14.6,' // int2string(ncc+1,'(i3)') // 'es19.7)'
	do i=1,nlam-1
		write(30,form) sqrt(lam(i)*lam(i+1))/micron,
c     &					flux(0:ncc,i)/(Fstar(i)*1d23/distance**2)
     &					((phase(1,j,i)+flux(j,i))/(Fstar(i)*1d23/distance**2),j=0,ncc)
	enddo
	close(unit=30)


	filename=trim(outputdir) // "trans"
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,RECL=1000)
	if(nclouds.gt.0) then
		form='("#",a13,a19,' // trim(int2string(ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]","Rp^2/Rstar^2",docloud0(1:nclouds,1:ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","Rp^2/Rstar^2"
	endif
	form='(f14.6,' // int2string(ncc+1,'(i3)') // 'es19.7)'
	do i=1,nlam-1
		write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					obsA(0:ncc,i)/(pi*Rstar**2)
	enddo
	close(unit=30)


	if(.not.domakeai) then
		filename=trim(outputdir) // "phase"
		call output("Writing spectrum to: " // trim(filename))
		open(unit=30,file=filename,RECL=6000)
		form='("#",a13,' // trim(int2string(nphase,'(i4)')) // 
     &				 '("   flux(",f5.1,") [Jy]"),"         fstar [Jy]")'
		write(30,form) "lambda [mu]",theta(1:nphase)
		form='(f14.6,' // int2string(nphase+1,'(i3)') // 'es19.7)'
		do i=1,nlam-1
			write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					phase(1:nphase,0,i)+flux(0,i),
     &					Fstar(i)*1d23/distance**2
		enddo
		close(unit=30)
	endif

	filename=trim(outputdir) // "transC"
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,RECL=1000)
	write(30,'("#",a13,a19,a26)') "lambda [mu]","Rp^2/Rstar^2","Rp^2/Rstar^2 using eclipse"
	form='(f14.6,es19.7,es19.7)'
	do i=1,nlam-1
		Fp1=(phase(nphase-1,0,i)+flux(0,i))/(Fstar(i)*1d23/distance**2)
		Fp2=(phase(nphase,0,i)+flux(0,i))/(Fstar(i)*1d23/distance**2)
		ApAs=obsA(0,i)/(pi*Rstar**2)
		write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					(Fp1-Fp2+ApAs)/(Fp1+1d0),ApAs-Fp2
	enddo
	close(unit=30)



	filename=trim(outputdir) // "tau1depth"
	call output("Writing tau1depth to: " // trim(filename))
	open(unit=30,file=filename,RECL=1000)
	if(nclouds.gt.0) then
		form='("#",a13,' // trim(int2string(ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]",docloud0(1:nclouds,1:ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","P [bar]"
	endif
	form='(f14.6,' // int2string(ncc,'(i3)') // 'es19.7)'
	do i=1,nlam-1
		write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					tau1depth(1:ncc,i)
	enddo
	close(unit=30)


	filename=trim(outputdir) // "cloudtau"
	call output("Writing cloud optical depth to: " // trim(filename))
	open(unit=30,file=filename,RECL=1000)
	if(nclouds.gt.0) then
		form='("#",a13,' // trim(int2string(ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]",docloud0(1:nclouds,1:ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","optical depth"
	endif
	form='(f14.6,' // int2string(ncc,'(i3)') // 'es19.7)'
	do i=1,nlam-1
		write(30,form) sqrt(lam(i)*lam(i+1))/micron,
     &					cloudtau(1:ncc,i)
	enddo
	close(unit=30)


	if(specresfile.ne.' ') then

	open(unit=20,file=specresfile)
	j=1
	ilam=1
1	read(20,*,end=2,err=1) x,specres_obs
	if(x.gt.(lam(1)/micron).and.x.lt.(lam(nlam)/micron)) ilam=ilam+1
	j=j+1
	goto 1
2	nlamR=ilam-1
	nj=j-1
	rewind(20)
	allocate(lamR(nlamR))
	allocate(specR(nlamR))
	allocate(specRexp(nlamR))
	allocate(specSNR(nlamR))
	allocate(spec(nphase,nlamR))
	ilam=1
	do j=1,nj
3		read(20,'(a1000)',err=3) line
		specres_obs=specres
		expspecres_obs=2d0
		SNR=1d8
		read(line,*,err=3,end=4) x,specres_obs,expspecres_obs,SNR
4		if(x.gt.(lam(1)/micron).and.x.lt.(lam(nlam)/micron)) then
			lamR(ilam)=x*micron
			specR(ilam)=specres_obs
			specRexp(ilam)=expspecres_obs
			specSNR(ilam)=SNR
			ilam=ilam+1
		endif
	enddo
	close(unit=20)

	filename=trim(outputdir) // "obs_emisR"
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,RECL=6000)
	form='("#",a13,' // trim(int2string(nphase,'(i4)')) // 
     &				 '("   flux(",f5.1,") [Jy]"),"         fstar [Jy]")'
	write(30,form) "lambda [mu]",theta(1:nphase)
	form='(f14.6,' // int2string(nphase+1,'(i3)') // 'es19.7)'
	do i=1,nphase
		call regridspecres(lam,phase(i,0,1:nlam-1)+flux(0,1:nlam-1),nlam-1,
     &						lamR,spec(i,1:nlamR),specR,specRexp,nlamR)
	enddo
	spec=spec/(Fstar(i)*1d23/distance**2)
	do j=1,nphase
		do i=1,nlamR
			spec(j,i)=spec(j,i)*(1d0+gasdev(idum)/specSNR(i))
		enddo
	enddo
	do i=1,nlamR
		write(30,form) lamR(i)/micron,spec(1:nphase,i)
	enddo
	close(unit=30)

	filename=trim(outputdir) // "obs_trans"
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,RECL=1000)
	write(30,'("#",a13,3a19)') "lambda [mu]","Rp^2/Rstar^2","dlam [mu]","error"
	form='(f14.6,3es19.7)'
	call regridspecres(lam,obsA(0,1:nlam-1),nlam-1,
     &					lamR,spec(1,1:nlamR),specR,specRexp,nlamR)
	spec=spec/(pi*Rstar**2)
	do i=1,nlamR
		spec(1,i)=spec(1,i)*(1d0+gasdev(idum)/specSNR(i))
	enddo
	do i=1,nlamR
		write(30,form) lamR(i)/micron,spec(1,i),lamR(i)/micron/specR(i),spec(1,i)/specSNR(i)
	enddo
	close(unit=30)

	deallocate(lamR)
	deallocate(specR)
	deallocate(specRexp)
	deallocate(specSNR)
	deallocate(spec)
	endif
	



	deallocate(docloud0)
	deallocate(theta)
	call output("==================================================================")
	
	
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
	write(30,'("# Pressure:    ",es10.3," bar")') P(ir)
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
