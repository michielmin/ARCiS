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
		real*8,allocatable :: ktable(:,:,:,:),P(:),T(:)
		real*8,allocatable :: g(:),wg(:),Cp(:,:)
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
	real*8  nullval,tot2,w1,ww,Pl,Planck
	real*8,allocatable :: lamF(:),Ktemp(:,:,:,:),temp(:,:,:),wtemp(:),tot(:,:)
	logical anynull,truefalse
	integer naxes(4)
	character*500 filename
	integer ig,ilam,iT,iP,imol,i,j,i1,i2,ngF
	integer*4 hdutype

	if(.not.allocated(Ktable)) allocate(Ktable(nmol))

	! Get an unused Logical Unit Number to use to open the FITS file.
	status=0

	call ftgiou (unit,status)
	! Open file
	readwrite=0
	status=0
	blocksize=0
	filename=trim(opacitydir) // "opacity"
	filename=trim(filename) // "_" // trim(molname(imol))
	filename=trim(filename) // ".fits"
	inquire(file=trim(filename),exist=truefalse)
	if(truefalse) then
		call ftopen(unit,trim(filename),readwrite,blocksize,status)
	endif
	if (.not.truefalse.or.status /= 0) then
		readwrite=0
		status=0
		blocksize=0
		print*,unit
		filename=trim(opacitydir) // "opacity"
		filename=trim(filename) // "_" // trim(molname(imol))
		filename=trim(filename) // ".fits.gz"
		inquire(file=trim(filename),exist=truefalse)
		if(truefalse) then
			call ftopen(unit,trim(filename),readwrite,blocksize,status)
		endif
		if (.not.truefalse.or.status /= 0) then
			call output("Opacity file not available: " // trim(filename))
			call output("setting opacities to 0")
			Ktable(imol)%available=.false.
			return
		endif
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
	call output("   wavelength range: " // trim(dbl2string(Ktable(imol)%lam1*1d4,'(es8.2)')) // " - " 
     &		// trim(dbl2string(Ktable(imol)%lam2*1d4,'(es8.2)')) // " micron")

	firstpix=1

	if(.not.allocated(Ktable(imol)%g)) then
		allocate(Ktable(imol)%g(Ktable(imol)%ng))
		allocate(Ktable(imol)%wg(Ktable(imol)%ng))
		call gauleg(0d0,1d0,Ktable(imol)%g,Ktable(imol)%wg,Ktable(imol)%ng)
	endif

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
	if(.not.allocated(Ktable(imol)%ktable)) allocate(Ktable(imol)%ktable(nlam,ng,naxes(3),naxes(4)))
	allocate(Ktemp(naxes(1),naxes(2),naxes(3),naxes(4)))

	call ftgpvd(unit,group,firstpix,npixels,nullval,Ktemp,anynull,status)

	allocate(lamF(Ktable(imol)%nlam+1))

	do ilam=1,Ktable(imol)%nlam+1
		lamF(ilam)=10d0**(log10(Ktable(imol)%lam1)+
     &			log10(Ktable(imol)%lam2/Ktable(imol)%lam1)*real(ilam-1)/real(Ktable(imol)%nlam))
	enddo

!$OMP PARALLEL IF(nlam.gt.200)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,i1,i2,i,ngF,ig,temp,j,tot,tot2,wtemp,ww,w1,iT,iP)
!$OMP& SHARED(nlam,Ktable,lam,lamF,imol,ng,gg,wgg,Ktemp)
	allocate(temp(Ktable(imol)%ng*Ktable(imol)%nlam,Ktable(imol)%nT,Ktable(imol)%nP))
	allocate(wtemp(Ktable(imol)%ng*Ktable(imol)%nlam))
	allocate(tot(Ktable(imol)%nT,Ktable(imol)%nP))
!$OMP DO
	do ilam=1,nlam-1
		i1=0
		i2=0
		do i=1,Ktable(imol)%nlam
			if(lam(ilam).ge.lamF(i).and.lam(ilam).lt.lamF(i+1)) i1=i
			if(lam(ilam+1).ge.lamF(i).and.lam(ilam+1).lt.lamF(i+1)) i2=i
		enddo
		if(i1.gt.0.and.i2.gt.0) then
			ngF=0
			do i=i1,i2
				if(i1.eq.i2) then
					ww=1d0
				else if(i.eq.i1) then
					ww=abs(lam(ilam)-lamF(i+1))
				else if(i.eq.i2) then
					ww=abs(lam(ilam)-lamF(i))
				else
					ww=abs(lamF(i)-lamF(i+1))
				endif
				do ig=1,Ktable(imol)%ng
					ngF=ngF+1
					temp(ngF,1:Ktable(imol)%nT,1:Ktable(imol)%nP)=Ktemp(i,ig,1:Ktable(imol)%nT,1:Ktable(imol)%nP)
					wtemp(ngF)=ww*Ktable(imol)%wg(ig)
				enddo
			enddo
			tot=0d0
			do ig=1,ngF
				tot(1:Ktable(imol)%nT,1:Ktable(imol)%nP)=tot(1:Ktable(imol)%nT,1:Ktable(imol)%nP)
     &					+temp(ig,1:Ktable(imol)%nT,1:Ktable(imol)%nP)*wtemp(ig)
			enddo
			tot=tot/sum(wtemp(1:ngF))
			call sortw(temp,wtemp,ngF)
			if(ng.eq.1) then
				Ktable(imol)%ktable(ilam,1,1:Ktable(imol)%nT,1:Ktable(imol)%nP)=tot(1:Ktable(imol)%nT,1:Ktable(imol)%nP)
			else
				do ig=2,ngF
					wtemp(ig)=wtemp(ig)+wtemp(ig-1)
				enddo
				wtemp(1:ngF)=wtemp(1:ngF)/wtemp(ngF)
				do ig=1,ng
					call hunt(wtemp,ngF,gg(ig),j)
					if(j.eq.0) then
						Ktable(imol)%ktable(ilam,ig,1:Ktable(imol)%nT,1:Ktable(imol)%nP)=
     &			temp(1,1:Ktable(imol)%nT,1:Ktable(imol)%nP)
					else
						w1=(gg(ig)-wtemp(j+1))/(wtemp(j)-wtemp(j+1))
						Ktable(imol)%ktable(ilam,ig,1:Ktable(imol)%nT,1:Ktable(imol)%nP)=
     &			temp(j,1:Ktable(imol)%nT,1:Ktable(imol)%nP)*w1+temp(j+1,1:Ktable(imol)%nT,1:Ktable(imol)%nP)*(1d0-w1)
					endif
				enddo
				do iT=1,Ktable(imol)%nT
				do iP=1,Ktable(imol)%nP
					tot2=0d0
					do ig=1,ng
						tot2=tot2+wgg(ig)*Ktable(imol)%ktable(ilam,ig,iT,iP)
					enddo
					if(tot2.ne.0d0) then
						Ktable(imol)%ktable(ilam,1:ng,iT,iP)=Ktable(imol)%ktable(ilam,1:ng,iT,iP)*tot(iT,iP)/tot2
					else
						Ktable(imol)%ktable(ilam,1:ng,iT,iP)=tot(iT,iP)
					endif
				enddo
				enddo
			endif
		else
			Ktable(imol)%ktable(ilam,1:ng,1:Ktable(imol)%nT,1:Ktable(imol)%nP)=0d0
		endif
	enddo
!$OMP END DO
	deallocate(temp)
	deallocate(wtemp)
	deallocate(tot)
!$OMP FLUSH
!$OMP END PARALLEL
	deallocate(Ktemp,lamF)

	!------------------------------------------------------------------------
	! HDU 1: temperature
	!------------------------------------------------------------------------
	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)	

	! Check dimensions
	call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)

	npixels=naxes(1)

	! read_image

	if(.not.allocated(Ktable(imol)%T)) allocate(Ktable(imol)%T(Ktable(imol)%nT))
	call ftgpvd(unit,group,firstpix,npixels,nullval,Ktable(imol)%T,anynull,status)

	!------------------------------------------------------------------------
	! HDU 2: pressure
	!------------------------------------------------------------------------
	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)	

	! Check dimensions
	call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)

	npixels=naxes(1)

	! read_image

	if(.not.allocated(Ktable(imol)%P)) allocate(Ktable(imol)%P(Ktable(imol)%nP))
	call ftgpvd(unit,group,firstpix,npixels,nullval,Ktable(imol)%P,anynull,status)

   				 
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
	
	allocate(Ktable(imol)%Cp(Ktable(imol)%nT,Ktable(imol)%nP))
	do iT=1,Ktable(imol)%nT
		Ktable(imol)%Cp(iT,1:Ktable(imol)%nP)=0d0
		do ilam=1,nlam-1
			Pl=Planck(Ktable(imol)%T(iT),freq(ilam))
			do iP=1,Ktable(imol)%nP
				do ig=1,ng
					Ktable(imol)%Cp(iT,iP)=Ktable(imol)%Cp(iT,iP)+wgg(ig)*dfreq(ilam)*Pl*Ktable(imol)%ktable(ilam,ig,iT,iP)
				enddo
			enddo
			tot2=tot2+dfreq(ilam)*Pl
		enddo
		Ktable(imol)%Cp(iT,1:Ktable(imol)%nP)=Ktable(imol)%Cp(iT,1:Ktable(imol)%nP)/tot2
	enddo

	return

	end

	subroutine GetCplanck(Cplanck,ir,T0)
	use GlobalSetup
	use OpacityFITSdata
	implicit none
	integer iT,iP,imol,ir
	real*8 Cplanck,wP1,wP2,wT1,wT2,T0

 	Cplanck=0d0
	do imol=1,nmol
		if(mixrat_r(ir,imol).gt.0d0.and.Ktable(imol)%available) then

		if(P(ir).lt.Ktable(imol)%P1) then
			iP=1
			wP1=1d0
			wP2=0d0
		else if(P(ir).gt.Ktable(imol)%P2) then
			iP=Ktable(imol)%nP-1
			wP1=0d0
			wP2=1d0
		else
			do iP=1,Ktable(imol)%nP
				if(P(ir).ge.Ktable(imol)%P(iP).and.P(ir).lt.Ktable(imol)%P(iP+1)) exit
			enddo
			wP1=1d0-log10(P(ir)/Ktable(imol)%P(iP))/log10(Ktable(imol)%P(iP+1)/Ktable(imol)%P(iP))
			wP2=1d0-wP1
		endif

		if(T0.lt.Ktable(imol)%T1) then
			iT=1
			wT1=1d0
			wT2=0d0
		else if(T0.gt.Ktable(imol)%T2) then
			iT=Ktable(imol)%nT-1
			wT1=0d0
			wT2=1d0
		else
			do iT=1,Ktable(imol)%nT
				if(T0.ge.Ktable(imol)%T(iT).and.T0.lt.Ktable(imol)%T(iT+1)) exit
			enddo
			wT1=1d0-log10(T0/Ktable(imol)%T(iT))/log10(Ktable(imol)%T(iT+1)/Ktable(imol)%T(iT))
			wT2=1d0-wT1
		endif

		Cplanck=Cplanck+mixrat_r(ir,imol)*(Ktable(imol)%Cp(iT,iP)*wT1*wP1+
     &						  Ktable(imol)%Cp(iT+1,iP)*wT2*wP1+
     &						  Ktable(imol)%Cp(iT,iP+1)*wT1*wP2+
     &						  Ktable(imol)%Cp(iT+1,iP+1)*wT2*wP2)

		endif
	enddo
	Cplanck=Cplanck*Ndens(ir)

	return
	end

	subroutine ReadOpacityFITS(kappa_mol,imol,ir)
	use GlobalSetup
	use OpacityFITSdata
	implicit none
	character*80 comment,errmessage
	character*30 errtext
	integer ig,ilam,iT,iP,imol,i,j,ir,ngF,i1,i2
	real*8 kappa_mol(ng,nmol,nlam),wP1,wP2,wT1,wT2,x1,x2,tot,tot2,random,w1,ww
	real*8,allocatable :: temp(:),wtemp(:)

	if(.not.Ktable(imol)%available.or..not.includemol(imol)) then
		kappa_mol(1:ng,imol,1:nlam)=0d0
		return
	endif

	if(P(ir).lt.Ktable(imol)%P1) then
		iP=1
		wP1=1d0
		wP2=0d0
	else if(P(ir).gt.Ktable(imol)%P2) then
		iP=Ktable(imol)%nP-1
		wP1=0d0
		wP2=1d0
	else
		do iP=1,Ktable(imol)%nP
			if(P(ir).ge.Ktable(imol)%P(iP).and.P(ir).lt.Ktable(imol)%P(iP+1)) exit
		enddo
		wP1=1d0-log10(P(ir)/Ktable(imol)%P(iP))/log10(Ktable(imol)%P(iP+1)/Ktable(imol)%P(iP))
		wP2=1d0-wP1
	endif

	if(T(ir).lt.Ktable(imol)%T1) then
		iT=1
		wT1=1d0
		wT2=0d0
	else if(T(ir).gt.Ktable(imol)%T2) then
		iT=Ktable(imol)%nT-1
		wT1=0d0
		wT2=1d0
	else
		do iT=1,Ktable(imol)%nT
			if(T(ir).ge.Ktable(imol)%T(iT).and.T(ir).lt.Ktable(imol)%T(iT+1)) exit
		enddo
		wT1=1d0-log10(T(ir)/Ktable(imol)%T(iT))/log10(Ktable(imol)%T(iT+1)/Ktable(imol)%T(iT))
		wT2=1d0-wT1
	endif
 
	do ilam=1,nlam-1
		do ig=1,ng
			kappa_mol(ig,imol,ilam)=Ktable(imol)%ktable(ilam,ig,iT,iP)*wT1*wP1+
     &						  Ktable(imol)%ktable(ilam,ig,iT+1,iP)*wT2*wP1+
     &						  Ktable(imol)%ktable(ilam,ig,iT,iP+1)*wT1*wP2+
     &						  Ktable(imol)%ktable(ilam,ig,iT+1,iP+1)*wT2*wP2
		enddo
	enddo

	return
	end



