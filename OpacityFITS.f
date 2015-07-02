	subroutine WriteOpacityFITS()
	use GlobalSetup
	IMPLICIT NONE
c================================
c for fitsio
	integer status,unit,blocksize,bitpix,naxis,naxes(4)
	integer group,fpixel,nelements
	logical simple,extend,truefalse
c================================
	character*500 filename,command
	integer imol,iT,iP,ig,ilam,ir,i
	real*8,allocatable :: array(:,:,:,:)

	write(command,'("mkdir -p ",a)') trim(opacitydir)
	call system(command)

	filename=trim(opacitydir) // "opacity"
	do i=1,nmol
		if(mixrat_r(1,i).gt.0d0) then
			filename=trim(filename) // "_" // trim(molname(i))
		endif
	enddo
	filename=trim(filename) // ".fits.gz"

	inquire(file=filename,exist=truefalse)
	if(truefalse) then
		call output("FITS file already exists, overwriting")
		open(unit=90,file=filename)
		close(unit=90,status='delete')
	endif

	
	status=0
C	 Get an unused Logical Unit Number to use to create the FITS file
	call ftgiou(unit,status)
C	 create the new empty FITS file
	blocksize=1
	call ftinit(unit,filename,blocksize,status)

	simple=.true.
	extend=.true.
	group=1
	fpixel=1

	bitpix=-64
	naxis=4
	naxes(1)=nlam
	naxes(2)=ng
	naxes(3)=nTom
	naxes(4)=nPom
	nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
	allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))

	! Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Write optional keywords to the header

	call ftpkyd(unit,'Tmin',Tmin,8,'[K]',status)
	call ftpkyd(unit,'Tmax',Tmax,8,'[K]',status)
	call ftpkyd(unit,'Pmin',Pmin,8,'[atm]',status)
	call ftpkyd(unit,'Pmax',Pmax,8,'[atm]',status)
	call ftpkyd(unit,'l_min',lam(1),8,'[micron]',status)
	call ftpkyd(unit,'l_max',lam(nlam),8,'[micron]',status)

	call ftpkyj(unit,'nT',nTom,' ',status)
	call ftpkyj(unit,'nP',nPom,' ',status)
	call ftpkyj(unit,'nlam',nlam,' ',status)
	call ftpkyj(unit,'ng',ng,' ',status)


	!  Write the array to the FITS file.

	!------------------------------------------------------------------------------
	! HDU 0: opacities 
	!------------------------------------------------------------------------------

	ir=0
	do iP=1,nPom
		do iT=1,nTom
			ir=ir+1
			do ilam=1,nlam
				do ig=1,ng
					array(ilam,ig,iT,iP)=Cabs(ir,ilam,ig)
				enddo
			enddo
		enddo
	enddo

	call ftpprd(unit,group,fpixel,nelements,array,status)

	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 1: temperature grid
	!------------------------------------------------------------------------------

	bitpix=-64
	naxis=1
	naxes(1)=nTom
	nelements=naxes(1)

	allocate(array(nTom,1,1,1))

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do iT=1,nTom
		array(iT,1,1,1)=10d0**(log10(Tmin)+log10(Tmax/Tmin)*real(iT-1)/real(nTom-1))
	enddo

	!  Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,array,status)

	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 2: pressure grid
	!------------------------------------------------------------------------------

	bitpix=-64
	naxis=1
	naxes(1)=nTom
	nelements=naxes(1)

	allocate(array(nPom,1,1,1))

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do iP=1,nPom
		array(iP,1,1,1)=10d0**(log10(Pmin)+log10(Pmax/Pmin)*real(iP-1)/real(nPom-1))
	enddo

	!  Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,array,status)

	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 3: wavelength grid
	!------------------------------------------------------------------------------

	bitpix=-64
	naxis=1
	naxes(1)=nlam
	nelements=naxes(1)

	allocate(array(nlam,1,1,1))

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do ilam=1,nlam
		array(ilam,1,1,1)=lam(ilam)
	enddo

	!  Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,array,status)

	deallocate(array)


	!------------------------------------------------------------------------------
	!  Close the file and free the unit number.
	!------------------------------------------------------------------------------

	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	if (status.gt.0) then
	   call output('error in export to fits file' // int2string(status,'(i6)'))
	end if


	return
	end


	module OpacityFITSdata
	IMPLICIT NONE
	
	type databaseKtable
		real*8,allocatable :: ktable(:,:,:,:)
		integer nlam,nP,nT,ng
		real*8 lam1,lam2,P1,P2,T1,T2
		logical available
	end type databaseKtable

	type(databaseKtable),allocatable :: Ktable(:)
	
	end module OpacityFITSdata


	subroutine InitReadOpacityFITS(imol)
	use GlobalSetup
	use OpacityFITSdata
	implicit none
	character*80 comment,errmessage
	character*30 errtext
	integer status,stat2,stat3,readwrite,unit,blocksize,nfound,group
	integer firstpix,nbuffer,npixels
	integer istat,stat4,tmp_int,stat5,stat6
	real*8  nullval
	logical anynull
	integer naxes(4)
	character*500 filename
	integer ig,ilam,iT,iP,imol,i,j

	if(.not.allocated(Ktable)) allocate(Ktable(nmol))

	! Get an unused Logical Unit Number to use to open the FITS file.
	status=0

	call ftgiou (unit,status)
	! Open file
	readwrite=0
	filename=trim(opacitydir) // "opacity"
	filename=trim(filename) // "_" // trim(molname(imol))
	filename=trim(filename) // ".fits.gz"
	call ftopen(unit,filename,readwrite,blocksize,status)
	if (status /= 0) then
		call output("Opacity file not available: " // trim(filename))
		call output("setting opacities to 0")
		Ktable(imol)%available=.false.
		return
	endif
	Ktable(imol)%available=.true.
	group=1
	nullval=-999

	call ftgkyd(unit,'Pmin',Ktable(imol)%P1,comment,status)
	call ftgkyd(unit,'Pmax',Ktable(imol)%P2,comment,status)
	call ftgkyd(unit,'Tmin',Ktable(imol)%T1,comment,status)
	call ftgkyd(unit,'Tmax',Ktable(imol)%T2,comment,status)
	call ftgkyd(unit,'l_min',Ktable(imol)%lam1,comment,status)
	call ftgkyd(unit,'l_max',Ktable(imol)%lam2,comment,status)

	call ftgkyj(unit,'nP',Ktable(imol)%nP,comment,status)
	call ftgkyj(unit,'nT',Ktable(imol)%nT,comment,status)
	call ftgkyj(unit,'nlam',Ktable(imol)%nlam,comment,status)
	call ftgkyj(unit,'ng',Ktable(imol)%ng,comment,status)

	call output("Reading in correlated k-tables for " // trim(molname(imol)))

	firstpix=1

	!------------------------------------------------------------------------
	! HDU0 : opacities
	!------------------------------------------------------------------------
	! Check dimensions
	call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

	! read_image

	naxes(1)=Ktable(imol)%nlam
	naxes(2)=Ktable(imol)%ng
	naxes(3)=Ktable(imol)%nT
	naxes(4)=Ktable(imol)%nP
	npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)
	if(.not.allocated(Ktable(imol)%ktable)) allocate(Ktable(imol)%ktable(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpvd(unit,group,firstpix,npixels,nullval,Ktable(imol)%ktable,anynull,status)
   				 
	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	!  Get the text string which describes the error
	if (status > 0) then
	   call ftgerr(status,errtext)
	   call output('error in reading fits file' // int2string(status,'(i6)'))

	   !  Read and print out all the error messages on the FITSIO stack
	   call ftgmsg(errmessage)
	   do while (errmessage .ne. ' ')
		  print *,errmessage
		  call ftgmsg(errmessage)
	   end do
	endif


	return

	end



	subroutine ReadOpacityFITS(kappa_mol,imol,ir)
	use GlobalSetup
	use OpacityFITSdata
	implicit none
	character*80 comment,errmessage
	character*30 errtext
	integer ig,ilam,iT,iP,imol,i,j,ir,ngF,i1,i2
	real*8 kappa_mol(nmol,nlam,ng),wP1,wP2,wT1,wT2,x1,x2
	real*8,allocatable :: temp(:),lamF(:)

	if(.not.Ktable(imol)%available) then
		kappa_mol(imol,1:nlam,1:ng)=0d0
		return
	endif

	allocate(temp(Ktable(imol)%ng*Ktable(imol)%nlam))
	allocate(lamF(Ktable(imol)%nlam))

	do ilam=1,Ktable(imol)%nlam
		lamF(ilam)=10d0**(log10(Ktable(imol)%lam1)+
     &			log10(Ktable(imol)%lam2/Ktable(imol)%lam1)*real(ilam-1)/real(Ktable(imol)%nlam-1))
	enddo

	iP=(log10(P(ir)/Ktable(imol)%P1)/log10(Ktable(imol)%P2/Ktable(imol)%P1))*real(Ktable(imol)%nP-1)+1d0
	if(iP.lt.1) then
		iP=1
		wP1=1d0
		wP2=0d0
	else if(iP.ge.Ktable(imol)%nP) then
		iP=Ktable(imol)%nP-1
		wP1=0d0
		wP2=1d0
	else
		x1=10d0**(log10(Ktable(imol)%P1)+log10(Ktable(imol)%P2/Ktable(imol)%P1)*real(iP-1)/real(Ktable(imol)%nP-1))
		x2=10d0**(log10(Ktable(imol)%P1)+log10(Ktable(imol)%P2/Ktable(imol)%P1)*real(iP)/real(Ktable(imol)%nP-1))
		wP1=1d0-log10(P(ir)/x1)/log10(x2/x1)
		wP2=1d0-wP1
	endif
	iT=(log10(T(ir)/Ktable(imol)%T1)/log10(Ktable(imol)%T2/Ktable(imol)%T1))*real(Ktable(imol)%nT-1)+1d0
	if(iT.lt.1) then
		iT=1
		wT1=1d0
		wT2=0d0
	else if(iT.ge.Ktable(imol)%nT) then
		iT=Ktable(imol)%nT-1
		wT1=0d0
		wT2=1d0
	else
		x1=10d0**(log10(Ktable(imol)%T1)+log10(Ktable(imol)%T2/Ktable(imol)%T1)*real(iT-1)/real(Ktable(imol)%nT-1))
		x2=10d0**(log10(Ktable(imol)%T1)+log10(Ktable(imol)%T2/Ktable(imol)%T1)*real(iT)/real(Ktable(imol)%nT-1))
		wT1=1d0-log10(T(ir)/x1)/log10(x2/x1)
		wT2=1d0-wT1
	endif
 
	do ilam=1,nlam-1
		i1=0
		i2=0
		do i=2,Ktable(imol)%nlam-1
			if(lam(ilam).le.lamF(i).and.lam(ilam).gt.lamF(i-1)) i1=i
			if(lam(ilam+1).ge.lamF(i).and.lam(ilam+1).lt.lamF(i+1)) i2=i
		enddo
		if(i1.eq.i2) i2=i1+1
		if(i1.gt.0) then
			ngF=0
			do i=i1,i2-1
				do ig=1,Ktable(imol)%ng
					ngF=ngF+1
					temp(ngF)=Ktable(imol)%ktable(i,ig,iT,iP)*wT1*wP1+
     &						  Ktable(imol)%ktable(i,ig,iT+1,iP)*wT2*wP1+
     &						  Ktable(imol)%ktable(i,ig,iT,iP+1)*wT1*wP2+
     &						  Ktable(imol)%ktable(i,ig,iT+1,iP+1)*wT2*wP2
				enddo
			enddo
			call sort(temp,ngF)
			do ig=1,ng
				j=1+real(ngF-1)*real(ig-1)/real(ng-1)
				kappa_mol(imol,ilam,ig)=temp(j)
			enddo
		else
			kappa_mol(imol,ilam,1:ng)=0d0
		endif
	enddo

	return
	end



