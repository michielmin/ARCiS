	subroutine ReadHITRAN()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	character*12 imol,iiso,nu,S,A,gamma_air,gamma_self,E,n,delta,gu,gl
	character*100 dummy
	integer i,j,k,maxiiso,it
	type(Line),pointer :: L
	real*8 scale

	open(unit=30,file=HITRANfile,RECL=500)

	nlines=0
1	read(30,*,end=2)
	nlines=nlines+1
	goto 1
2	close(unit=30)

	allocate(Lines(nlines))

	call output("Reading HITRAN database")
	call output("number of lines: " // trim(int2string(nlines)))

	nmol=0
	maxiiso=0
	open(unit=30,file=HITRANfile,RECL=500)
	do i=1,nlines
		call tellertje(i,nlines)
		L => Lines(i)
		read(30,'(a2,a1,a12,a10,a10,a5,a5,a10,a4,a8,a79,a7,a7)') 
     &			imol,iiso,nu,S,A,gamma_air,gamma_self,E,n,delta,dummy,gu,gl
		read(imol,*) L%imol
		if(L%imol.gt.nmol) nmol=L%imol
		read(iiso,*) L%iiso
		if(L%iiso.gt.maxiiso) maxiiso=L%iiso
		read(nu,*) L%freq
		read(S,*) L%S
		read(A,*) L%Aul
		read(gamma_air,*) L%gamma_air
		read(gamma_self,*) L%gamma_self
		read(E,*) L%Eup
c		read(n,*) L%imol
c		read(delta,*) L%imol
		read(gu,*) L%gu
		read(gl,*) L%gl
	enddo

	print*,nmol,maxiiso
	allocate(niso(nmol))
	niso=0
	do i=1,nlines
		if(Lines(i)%iiso.gt.niso(Lines(i)%imol)) niso(Lines(i)%imol)=Lines(i)%iiso
	enddo

	nTZ=500
	allocate(ZZ(nmol,maxiiso,nTZ))
	allocate(TZ(nTZ))
	ZZ=0d0
	do it=1,nTZ
		TZ(it)=exp(log(71d0)+log(2900d0/71d0)*real(it-1)/real(nTZ-1))
		print*,TZ(it)
	enddo
	do j=1,nmol
		do k=1,niso(j)
			call TIPS_2011(j,k,295d0,scale)
			do it=1,nTZ
				print*,TZ(it)
				call TIPS_2011(j,k,TZ(it),ZZ(j,k,it))
				ZZ(j,k,it)=ZZ(j,k,it)/scale
			enddo
		enddo
	enddo
	
	return
	end
	