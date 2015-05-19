	subroutine ReadHITRAN()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	character*12 imol,iiso,nu,S,A,gamma_air,gamma_self,E,n,delta,gu,gl
	character*100 dummy
	integer i,j,k,maxiiso,it
	type(Line),pointer :: L
	real*8 scale

	call output("Reading HITRAN database")
	call output("file: " // trim(HITRANfile))

	open(unit=30,file=HITRANfile,RECL=500)

	nlines=0
1	read(30,'(a2)',end=2) imol
	read(imol,*) j
	if(j.le.nmol) then
		if(mixrat(j).gt.0d0) nlines=nlines+1
	endif
	goto 1
2	close(unit=30)

	allocate(Lines(nlines))

	call output("number of lines: " // trim(int2string(nlines)))

c	nmol=0
	maxiiso=0
	open(unit=30,file=HITRANfile,RECL=500)
	i=0
3	read(30,'(a2,a1,a12,a10,a10,a5,a5,a10,a4,a8,a79,a7,a7)',end=4) 
     &			imol,iiso,nu,S,A,gamma_air,gamma_self,E,n,delta,dummy,gu,gl
	read(imol,*) j
	if(j.le.nmol) then
	if(mixrat(j).gt.0d0) then
		i=i+1
		call tellertje(i,nlines)
		L => Lines(i)
		read(imol,*) L%imol
		if(L%imol.gt.nmol) nmol=L%imol
		read(iiso,*) L%iiso
		if(L%iiso.gt.maxiiso) maxiiso=L%iiso
		read(nu,*) L%freq
		read(S,*) L%S0
		read(A,*) L%Aul
		read(gamma_air,*) L%gamma_air
		read(gamma_self,*) L%gamma_self
		read(E,*) L%Elow
		read(n,*) L%n
c		read(delta,*) L%delta
		read(gu,*) L%gu
		read(gl,*) L%gl
	endif
	endif
	goto 3
4	close(unit=30)

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
	enddo
	do j=1,nmol
		do k=1,niso(j)
			call TIPS_2011(j,k,296d0,scale)
			do it=1,nTZ
				call TIPS_2011(j,k,TZ(it),ZZ(j,k,it))
				ZZ(j,k,it)=ZZ(j,k,it)/scale
			enddo
		enddo
	enddo

	if(ncia.gt.0) call output("Reading CIA opacities")
	do i=1,ncia
		call InitCIA(i)
		call output("CIA: " // trim(molname(CIA(i)%imol1)) // "-" // trim(molname(CIA(i)%imol1)))
	enddo

	cia_mixrat=-1d0
	do i=1,nmol
		cia_mixrat(i)=mixrat(i)
	enddo
c add Helium (arbitrary value for now...)
	cia_mixrat(48)=0.1
c set default for H2 to 1.0
	if(cia_mixrat(45).lt.0d0) cia_mixrat(45)=1d0
	
	return
	end
	