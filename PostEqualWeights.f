	subroutine PostEqualWeights()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 var(n_ret),error(n_ret)
	real*8,allocatable :: spectrans(:,:),specemis(:,:),specemisR(:,:),sorted(:)
	integer i,nmodels,ilam,im3,im1,ime,ip1,ip3
	
	open(unit=35,file=trim(outputdir) // "/post_equal_weights.dat",RECL=6000)
	
	i=1
1	read(35,*,end=2) var(1:n_ret)
	i=i+1
	goto 1
2	close(unit=35)
	nmodels=i-1
	print*,nmodels

	allocate(spectrans(nmodels,nlam))
	allocate(specemis(nmodels,nlam))
	allocate(specemisR(nmodels,nlam))
	allocate(sorted(nmodels))

	spectrans=0d0
	specemis=0d0
	specemisR=0d0

	open(unit=35,file=trim(outputdir) // "/post_equal_weights.dat",RECL=6000)

	do i=1,nmodels

	read(35,*) var(1:n_ret)

	call output("model number: " // int2string(i,'(i7)'))

	if(n2d.eq.0) then
		i2d=0
	else
		i2d=1
	endif

3	continue
	error=0d0
	call SetOutputMode(.false.)
	call MapRetrievalMN(var,error)
	call SetOutputMode(.true.)

	call InitDens()
	call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam,starfile)
	Fstar=Fstar*pi*Rstar**2
	call SetOutputMode(.false.)
	call ComputeModel(.true.)
	call SetOutputMode(.true.)
	
	if(i2d.eq.0) then
		spectrans(i,1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
	else if(i2d.le.2) then
		spectrans(i,1:nlam)=spectrans(i,1:nlam)+obsA(0,1:nlam)/(pi*Rstar**2)/2d0
	endif
	if(i2d.eq.0.or.i2d.eq.3) then
		specemisR(i,1:nlam)=(phase(1,0,1:nlam-1)+flux(0,1:nlam-1))/(Fstar(1:nlam-1)*1d23/distance**2)
		specemis(i,1:nlam)=phase(1,0,1:nlam-1)+flux(0,1:nlam-1)
	endif

	i2d=i2d+1
	if(i2d.le.n2d) goto 3

	
	if(i.gt.2) then
		open(unit=26,file=trim(outputdir) // "trans_limits",RECL=1000)
		open(unit=27,file=trim(outputdir) // "emis_limits",RECL=1000)
		open(unit=28,file=trim(outputdir) // "emisR_limits",RECL=1000)
		im1=real(i)/2d0-real(i)*0.682689492137086/2d0+0.5d0
		im3=real(i)/2d0-real(i)*0.997300203936740/2d0+0.5d0
		ip1=real(i)/2d0+real(i)*0.682689492137086/2d0+0.5d0
		ip3=real(i)/2d0+real(i)*0.997300203936740/2d0+0.5d0
		ime=real(i)/2d0+0.5d0
		if(im1.lt.1) im1=1
		if(im3.lt.1) im3=1
		if(ime.lt.1) ime=1
		if(ip1.gt.i) ip1=i
		if(ip3.gt.i) ip3=i
		if(ime.gt.i) ime=i
		do ilam=1,nlam-1
			sorted(1:i)=spectrans(1:i,ilam)
			call sort(sorted,i)
			write(26,*) lam(ilam)*1d4,sorted(im3),sorted(im1),sorted(ime),sorted(ip1),sorted(ip3)
			sorted(1:i)=specemis(1:i,ilam)
			call sort(sorted,i)
			write(27,*) lam(ilam)*1d4,sorted(im3),sorted(im1),sorted(ime),sorted(ip1),sorted(ip3)
			sorted(1:i)=specemisR(1:i,ilam)
			call sort(sorted,i)
			write(28,*) lam(ilam)*1d4,sorted(im3),sorted(im1),sorted(ime),sorted(ip1),sorted(ip3)
		enddo
		close(unit=26)
		close(unit=27)
		close(unit=28)		
	endif

	enddo
	close(unit=35)

	deallocate(spectrans)
	deallocate(specemis)
	deallocate(specemisR)
	deallocate(sorted)
	
	return
	end
	