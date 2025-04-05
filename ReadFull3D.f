	module Read3DModule
	IMPLICIT NONE
	real*8,allocatable :: longb1(:),longb2(:),latb1(:),latb2(:)
	character*500,allocatable :: TPfile3D(:)
	integer nfiles3D
	end module

	subroutine GetIndxReadFull3D(i,ilong,ilatt)
	use Read3DModule
	use Struct3D
	use Constants
	IMPLICIT NONE
	integer i,ilong,ilatt
	real*8 lo,la
	
	lo=(long(ilong)+long(ilong+1))/2d0
	la=(latt(ilatt)+latt(ilatt+1))/2d0
	do i=1,nfiles3D
		if((lo.ge.longb1(i).and.lo.le.longb2(i)).or.
     &	  ((lo+2d0*pi).ge.longb1(i).and.(lo+2d0*pi).le.longb2(i))) then
			if(((la-pi/2d0).gt.latb1(i).and.(la-pi/2d0).le.latb2(i)).or. 
     &		   ((pi/2d0-la).gt.latb1(i).and.(pi/2d0-la).le.latb2(i))) return
		endif
	enddo

	return
	end
	
	subroutine SetTPfileReadFull3D(i)
	use GlobalSetup
	use Read3DModule
	use Struct3D
	IMPLICIT NONE
	integer i

	TPfile=TPfile3D(i)

	return
	end
	
	subroutine InitReadFull3D(dir)
	use Read3DModule
	use GlobalSetup
	use Struct3D
	IMPLICIT NONE
	character*500 dir,filename,line
	integer i
	
	filename=trim(dir) // "/overview.dat"
	open(unit=45,file=filename)
	read(45,*)
	read(45,*)
	read(45,*)
	read(45,*)
	read(45,*)
	read(45,*) nfiles3D
	n3D=nfiles3D
	allocate(longb1(nfiles3D),longb2(nfiles3D),latb1(nfiles3D),latb2(nfiles3D))
	allocate(TPfile3D(nfiles3D))
	do i=1,nfiles3D
		read(45,'(a500)') line(1:500)
		TPfile3D(i)=line
		read(45,*) longb1(i),longb2(i)
		read(45,*) latb1(i),latb2(i)
	enddo
	close(unit=45)
	
	return
	end
