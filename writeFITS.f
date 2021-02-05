	subroutine writefitsfile(filename,im,nlam,n)
	IMPLICIT NONE
	character*500 filename
	integer n,nlam
	real*8 im(n,n,nlam)

      integer status,unit,blocksize,bitpix,naxis,naxes(3)
      integer i,j,group,fpixel,nelements
      logical simple,extend,truefalse

	inquire(file=filename,exist=truefalse)
	if(truefalse) then
		call output("FITS file already exists, overwriting")
		open(unit=90,file=filename)
		close(unit=90,status='delete')
	endif

      status=0
C     Get an unused Logical Unit Number to use to create the FITS file
      call ftgiou(unit,status)
C     create the new empty FITS file
      blocksize=1
      call ftinit(unit,filename,blocksize,status)

C     initialize parameters about the FITS image (IMDIM x IMDIM 64-bit reals)
      simple=.true.
      bitpix=-64
	naxes(1)=n
	naxes(2)=n
	if(nlam.gt.1) then
		naxis=3
		naxes(3)=nlam
	else
		naxis=2
		naxes(3)=1
	endif
      extend=.true.

C     write the required header keywords
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

C     write the array to the FITS file
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)*naxes(3)
      call ftpprd(unit,group,fpixel,nelements,im,status)

C     close the file and free the unit number
      call ftclos(unit, status)
      call ftfiou(unit, status)


	return
	end


	subroutine writeContribution(filename,ptrace,lam,obsA,flux,ntrace,nlam)
	IMPLICIT NONE
	integer ntrace,nlam
	real*8 lam(nlam),ptrace(ntrace),obsA(ntrace,nlam),flux(ntrace,nlam)
	character*500 filename
	real*8,allocatable :: array(:,:)
	integer status,unit,blocksize,bitpix,naxis,naxes(2)
	integer group,fpixel,nelements,i,j
	logical simple,extend,truefalse

	if(filename(len_trim(filename)-4:len_trim(filename)).eq.'.fits') then
		filename=trim(filename)//'.gz'
	endif

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
	naxis=1
	naxes(1)=nlam
	naxes(2)=1
	nelements=naxes(1)*naxes(2)

	! Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Write optional keywords to the header

c	call ftpkyd(unit,'Rin',ZZ%Rin/AU,8,'[AU]',status)
c	call ftpkyd(unit,'Rout',ZZ%Rout/AU,8,'[AU]',status)

c	call ftpkyj(unit,'nR',ZZ%nr,' ',status)
c	call ftpkyj(unit,'nTheta',ZZ%nt,' ',status)
c	call ftpkyj(unit,'nPhi',ZZ%np,' ',status)
c	call ftpkyj(unit,'npart',npart,' ',status)
	call ftpkyj(unit,'nlam',nlam,' ',status)
	call ftpkyj(unit,'ntrace',ntrace,' ',status)

c	call ftpkyj(unit,'nHDU',nvars,' ',status)	
c	do i=1,nvars
c		write(hdu,'("HDU",i3)') i
c		call ftpkys(unit,hdu,trim(vars(i)),'',status)
c	enddo

	!  Write the array to the FITS file.

	!------------------------------------------------------------------------------
	! HDU 0: wavelength grid
	!------------------------------------------------------------------------------

	allocate(array(naxes(1),naxes(2)))

	do i=1,nlam
		array(i,1)=lam(i)
	enddo

	call ftpprd(unit,group,fpixel,nelements,array,status)
	
	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 1: ptrace grid
	!------------------------------------------------------------------------------
	! create new hdu
	call ftcrhd(unit, status)
	bitpix=-64
	naxes(1)=ntrace
	naxes(2)=1
	naxis=1
	nelements=naxes(1)*naxes(2)
	allocate(array(naxes(1),naxes(2)))
	do i=1,ntrace
		array(i,1)=ptrace(i)
	enddo

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	!  Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,array,status)

	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 2: optical depth
	!------------------------------------------------------------------------------
	! create new hdu
	call ftcrhd(unit, status)
	bitpix=-64
	naxes(1)=nlam
	naxes(2)=ntrace
	naxis=2
	nelements=naxes(1)*naxes(2)
	allocate(array(naxes(1),naxes(2)))
	do i=1,nlam
		do j=1,ntrace
			array(i,j)=obsA(j,i)
		enddo
	enddo

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	!  Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,array,status)

	deallocate(array)


	!------------------------------------------------------------------------------
	! HDU 3: flux contribution
	!------------------------------------------------------------------------------
	! create new hdu
	call ftcrhd(unit, status)
	bitpix=-64
	naxes(1)=nlam
	naxes(2)=ntrace
	naxis=2
	nelements=naxes(1)*naxes(2)
	allocate(array(naxes(1),naxes(2)))
	do i=1,nlam
		do j=1,ntrace
			array(i,j)=flux(j,i)
		enddo
	enddo

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	!  Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,array,status)

	deallocate(array)


	
	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	if (status.gt.0) then
	   print*,'error in export to fits file',status
	end if


	return
	end
	

