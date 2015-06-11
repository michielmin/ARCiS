	subroutine InitCIA(icia)
	use GlobalSetup
	IMPLICIT NONE
	integer i,j,iT,n,icia
	character*20 name,cmin,cmax,cn,cT
	real*8 x0,y0,x1,y1
	
	call output("Reading file:" // trim(CIA(icia)%filename))
	iT=0
	open(unit=20,file=CIA(icia)%filename,RECL=100)
1	read(20,'(a20,a10,a10,a7,a7)',err=1,end=2) name,cmin,cmax,cn,cT
	read(cn,*) n
	do i=1,n
		read(20,*)
	enddo
	iT=iT+1
	goto 1
2	close(unit=20)

	CIA(icia)%nT=iT
	allocate(CIA(icia)%T(CIA(icia)%nT))
	allocate(CIA(icia)%Cabs(CIA(icia)%nT,nlam))

	open(unit=20,file=CIA(icia)%filename,RECL=100)
	do iT=1,CIA(icia)%nT
3		read(20,'(a20,a10,a10,a7,a7)',err=3) name,cmin,cmax,cn,cT
		if(iT.eq.1) then
			CIA(icia)%name=name
		else
			if(CIA(icia)%name.ne.name) then
				print*,'error!'
				stop
			endif
		endif
		read(cT,*) CIA(icia)%T(iT)
		read(cn,*) n
		i=nlam
		read(20,*,end=102,err=1) x0,y0
103		if(x0.ge.freq(i)) then
			CIA(icia)%Cabs(iT,i)=y0
			i=i-1
			goto 103
		endif
		do j=2,n
			read(20,*) x1,y1
101			if(i.gt.0) then
				if(freq(i).le.x1.and.freq(i).ge.x0) then
					CIA(icia)%Cabs(iT,i)=y1+(freq(i)-x1)*(y0-y1)/(x0-x1)
					i=i-1
					goto 101
				endif
			endif
			x0=x1
			y0=y1
		enddo
102		continue
		do j=i,1,-1
			CIA(icia)%Cabs(iT,i)=CIA(icia)%Cabs(iT,i-1)*freq(i-1)/freq(j)
		enddo
	enddo
	close(unit=20)

	CIA(icia)%imol1=45
	CIA(icia)%imol2=45

	i=index(CIA(icia)%name,'-')
	read(CIA(icia)%name(1:i-1),*) name
	do i=1,48
		if(name.eq.molname(i)) then
			CIA(icia)%imol1=i
		endif
	enddo
	i=index(CIA(icia)%name,'-')
	read(CIA(icia)%name(i+1:20),*) name
	do i=1,48
		if(name.eq.molname(i)) then
			CIA(icia)%imol2=i
		endif
	enddo

	end
	
	
	