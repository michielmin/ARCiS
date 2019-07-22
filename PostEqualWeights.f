	subroutine PostEqualWeights()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 error(n_ret),random,starttime,stoptime,remaining,omp_get_wtime,sig,aver
	real*8,allocatable :: spectrans(:,:),specemis(:,:),specemisR(:,:),sorted(:)
	real*8,allocatable :: PTstruct(:,:),var(:,:)
	integer i,nmodels,ilam,im3,im1,ime,ip1,ip3,im2,ip2,ir,imodel,iobs,donmodels
	logical,allocatable :: done(:)
	
	open(unit=35,file=trim(outputdir) // "/post_equal_weights.dat",RECL=6000)
	
	i=1
1	read(35,*,end=2) error(1:n_ret)
	i=i+1
	goto 1
2	close(unit=35)
	nmodels=i-1
	print*,nmodels

	allocate(spectrans(0:nmodels,nlam))
	allocate(specemis(0:nmodels,nlam))
	allocate(specemisR(0:nmodels,nlam))
	allocate(PTstruct(0:nmodels,nr))
	allocate(sorted(nmodels))
	allocate(done(nmodels))
	allocate(var(nmodels,n_ret))

	spectrans=0d0
	specemis=0d0
	specemisR=0d0
	PTstruct=0d0

	open(unit=35,file=trim(outputdir) // "/post_equal_weights.dat",RECL=6000)

	do i=1,nmodels
		read(35,*) var(i,1:n_ret)
	enddo

	done=.false.
	
	call cpu_time(starttime)
#if USE_OPENMP
	starttime=omp_get_wtime()
#endif

	if(npew.lt.0) then
		donmodels=nmodels
	else
		donmodels=min(npew,nmodels)
	endif

	do i=0,donmodels

	if(i.eq.0) then
		call output("best fit model")
	else
5		imodel=random(idum)*nmodels+1
		if(done(imodel)) goto 5
		done(imodel)=.true.

		call cpu_time(stoptime)
#if USE_OPENMP
		stoptime=omp_get_wtime()
#endif

		call output("model number: " // int2string(imodel,'(i7)') 
     &			// "(" // trim(dbl2string(real(i)*100d0/real(nmodels),'(f5.1)')) // "%)")
		remaining=(stoptime-starttime)*real(donmodels-i)/real(i)
		if(remaining.gt.3600d0*24d0) then
			call output("time remaining: " // trim(dbl2string(remaining/3600d0/24d0,'(f6.1)')) // " days")
		else if(remaining.gt.3600d0) then
			call output("time remaining: " // trim(dbl2string(remaining/3600d0,'(f6.1)')) // " hours")
		else if(remaining.gt.60d0) then
			call output("time remaining: " // trim(dbl2string(remaining/60d0,'(f6.1)')) // " minutes")
		else
			call output("time remaining: " // trim(dbl2string(remaining,'(f6.1)')) // " seconds")
		endif
	endif

	if(n2d.eq.0) then
		i2d=0
	else
		i2d=1
	endif

3	continue
	error=0d0
	call SetOutputMode(.false.)
	if(i.ne.0) call MapRetrievalMN(var(imodel,1:n_ret),error)
	call SetOutputMode(.true.)

	call InitDens()
	call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam,starfile)
	Fstar=Fstar*pi*Rstar**2
	call SetOutputMode(.false.)
	call ComputeModel(.true.)
	call SetOutputMode(.true.)
	
	if(i2d.ne.0.and.nobs.ne.0) then
		do iobs=1,nobs
			if(ObsSpec(iobs)%i2d.eq.i2d) then
				select case(ObsSpec(iobs)%type)
					case("trans","transmission","transC")
						spectrans(i,1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
					case("emisr","emisR","emisa","emis","emission")
						specemisR(i,1:nlam)=(phase(1,0,1:nlam-1)+flux(0,1:nlam-1))/(Fstar(1:nlam-1)*1d23/distance**2)
						specemis(i,1:nlam)=phase(1,0,1:nlam-1)+flux(0,1:nlam-1)
				end select
			endif
		enddo
	else
		if(i2d.eq.0) then
			spectrans(i,1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
		else if(i2d.le.2) then
			spectrans(i,1:nlam)=spectrans(i,1:nlam)+obsA(0,1:nlam)/(pi*Rstar**2)/2d0
		endif
		if(i2d.eq.0.or.i2d.eq.3) then
			specemisR(i,1:nlam)=(phase(1,0,1:nlam-1)+flux(0,1:nlam-1))/(Fstar(1:nlam-1)*1d23/distance**2)
			specemis(i,1:nlam)=phase(1,0,1:nlam-1)+flux(0,1:nlam-1)
		endif
	endif

	i2d=i2d+1
	if(i2d.le.n2d) goto 3
	PTstruct(i,1:nr)=T(1:nr)
	
	if(i.gt.2.and.(100*(i/100).eq.i.or.i.eq.donmodels.or.i.le.10)) then
		im1=real(i)/2d0-real(i)*0.682689492137086/2d0+0.5d0
		im2=real(i)/2d0-real(i)*0.954499736103642/2d0+0.5d0
		im3=real(i)/2d0-real(i)*0.997300203936740/2d0+0.5d0
		ip1=real(i)/2d0+real(i)*0.682689492137086/2d0+0.5d0
		ip2=real(i)/2d0+real(i)*0.954499736103642/2d0+0.5d0
		ip3=real(i)/2d0+real(i)*0.997300203936740/2d0+0.5d0
		ime=real(i)/2d0+0.5d0
		if(im1.lt.1) im1=1
		if(im3.lt.1) im3=1
		if(ime.lt.1) ime=1
		if(ip1.gt.i) ip1=i
		if(ip3.gt.i) ip3=i
		if(ime.gt.i) ime=i
		open(unit=26,file=trim(outputdir) // "trans_limits",RECL=1000)
		open(unit=27,file=trim(outputdir) // "emis_limits",RECL=1000)
		open(unit=28,file=trim(outputdir) // "emisR_limits",RECL=1000)
		do ilam=1,nlam-1
			sorted(1:i)=spectrans(1:i,ilam)
			call sort(sorted,i)
			write(26,*) lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
			sorted(1:i)=specemis(1:i,ilam)
			call sort(sorted,i)
			write(27,*) lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
			sorted(1:i)=specemisR(1:i,ilam)
			call sort(sorted,i)
			write(28,*) lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
		enddo
		close(unit=26)
		close(unit=27)
		close(unit=28)		

		open(unit=26,file=trim(outputdir) // "PT_limits",RECL=1000)
		do ir=1,nr
			sorted(1:i)=PTstruct(1:i,ir)
			call sort(sorted,i)
			write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
		enddo
		close(unit=26)

		open(unit=26,file=trim(outputdir) // "trans_sigma",RECL=1000)
		open(unit=27,file=trim(outputdir) // "emis_sigma",RECL=1000)
		open(unit=28,file=trim(outputdir) // "emisR_sigma",RECL=1000)
		do ilam=1,nlam-1
			sorted(1:i)=spectrans(1:i,ilam)
			call stats(sorted,i,aver,sig)
			write(26,*) lam(ilam)*1d4,spectrans(0,ilam),sig,aver
			sorted(1:i)=specemis(1:i,ilam)
			call stats(sorted,i,aver,sig)
			write(27,*) lam(ilam)*1d4,specemis(0,ilam),sig,aver
			sorted(1:i)=specemisR(1:i,ilam)
			call stats(sorted,i,aver,sig)
			write(28,*) lam(ilam)*1d4,specemisR(0,ilam),sig,aver
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
	deallocate(PTstruct)
	deallocate(sorted)
	
	return
	end
	