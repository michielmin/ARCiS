	subroutine PostEqualWeights()
	use GlobalSetup
	use Constants
	use Struct3D
	IMPLICIT NONE
	real*8 error(n_ret),random,starttime,stoptime,remaining,omp_get_wtime,sig,aver,xmin,xmax
	real*8,allocatable :: spectrans(:,:),specemis(:,:),specemisR(:,:),sorted(:),hotspotshift_der(:)
	real*8,allocatable :: PTstruct(:,:),var(:,:),values(:,:),COratio_der(:),Z_der(:),cloudstruct(:,:),MMW_der(:)
	integer i,nmodels,ilam,im3,im1,ime,ip1,ip3,im2,ip2,ir,imodel,iobs,donmodels,j,iphase,imol,k
	logical,allocatable :: done(:)
	real*8,allocatable :: PTstruct3D(:,:,:),mixrat3D(:,:,:,:),phase3D(:,:,:),phase3DR(:,:,:),var3D(:,:,:)
	real*8,allocatable :: dPTstruct3D(:,:,:),dPTstruct(:,:),Kzz_struct(:,:),mol_struct(:,:,:),like(:),Tplanet(:)
	real*8,allocatable :: speccloudtau(:,:),speccloudK(:,:,:),specobs(:,:,:),Cov(:,:)
	real*8 lbest,x1(n_ret),x2(n_ret),ctrans,cmax,lm1,cmin,spectemp(nlam),gasdev
	integer i1,i2,ibest,info
	character*6000 line
	integer*4 counts, count_rate, count_max
	logical variablePgrid,multinestpost
	real*8 Pgrid(nr),Pg1(nr),Pg2(nr),yy(nr),Pmin0,Pmax0,xx,xy
	character*500 lowkey
	integer ipmin,ipmax
	
	writefiles=.false.
	
	if(retrievaltype.eq.'MC'.or.retrievaltype.eq.'MCMC') then
		open(unit=35,file=trim(outputdir) // "/posteriorMCMC.dat",FORM="FORMATTED",ACCESS="STREAM")
		read(35,*)
		i=1
11		read(35,*,end=12) error(1:n_ret)
		i=i+1
		goto 11
12		close(unit=35)
		nmodels=i-1
	else
		inquire(file=trim(outputdir) // "/post_equal_weights.dat",exist=multinestpost)
		if(multinestpost) then
			open(unit=35,file=trim(outputdir) // "/post_equal_weights.dat",FORM="FORMATTED",ACCESS="STREAM")
		else
			open(unit=35,file=trim(outputdir) // "/post_equal_weights.txt",FORM="FORMATTED",ACCESS="STREAM")
		endif
		i=1
1		read(35,*,end=2) error(1:n_ret)
		i=i+1
		goto 1
2		close(unit=35)
		nmodels=i-1
	endif

	print*,nmodels

	allocate(spectrans(0:nmodels,nlam))
	allocate(specemis(0:nmodels,nlam))
	allocate(specemisR(0:nmodels,nlam))
	allocate(PTstruct(0:nmodels,nr))
	allocate(dPTstruct(0:nmodels,nr))
	allocate(Kzz_struct(0:nmodels,nr))
	if(dochemistry.or..true.) allocate(mol_struct(0:nmodels,nr,nmol))
	allocate(cloudstruct(0:nmodels,nr))
	allocate(speccloudtau(0:nmodels,nlam))
	allocate(speccloudK(0:nmodels,max(nclouds,1),nlam))
	allocate(values(0:nmodels,n_ret))
	allocate(COratio_der(0:nmodels))
	allocate(Z_der(0:nmodels))
	allocate(MMW_der(0:nmodels))
	allocate(Tplanet(0:nmodels))
	if(do3D) allocate(hotspotshift_der(0:nmodels))
	allocate(sorted(nmodels))
	allocate(done(nmodels))
	allocate(like(nmodels))
	allocate(var(nmodels,n_ret))
	if(do3D.and.fulloutput3D) then
		allocate(PTstruct3D(0:nmodels,0:nphase,nr))
		allocate(dPTstruct3D(0:nmodels,0:nphase,nr))
		allocate(mixrat3D(0:nmodels,0:nphase,nr,nmol))
		allocate(phase3D(0:nmodels,1:nphase,nlam))
		allocate(phase3DR(0:nmodels,1:nphase,nlam))
		allocate(var3D(0:nmodels,1:nlong,0:n_Par3D))
	endif
	if(useobsgrid) allocate(specobs(0:nmodels,nobs,nlam))

	spectrans=0d0
	specemis=0d0
	specemisR=0d0
	PTstruct=0d0
	dPTstruct=0d0
	cloudstruct=0d0
	speccloudtau=0d0
	speccloudK=0d0
	cmax=0d0

	if(retrievaltype.eq.'MC'.or.retrievaltype.eq.'MCMC') then
		open(unit=35,file=trim(outputdir) // "/posteriorMCMC.dat",FORM="FORMATTED",ACCESS="STREAM")
		read(35,*)
		do i=1,nmodels
			read(35,*) var(i,1:n_ret)
			like(i)=1d0
		enddo
	else
		inquire(file=trim(outputdir) // "/post_equal_weights.dat",exist=multinestpost)
		if(multinestpost) then
			open(unit=35,file=trim(outputdir) // "/post_equal_weights.dat",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nmodels
				read(35,*) var(i,1:n_ret),like(i)
			enddo
		else
			open(unit=35,file=trim(outputdir) // "/post_equal_weights.txt",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nmodels
				read(35,*) var(i,1:n_ret)
				like(i)=1d0
			enddo
		endif
	endif
	ipmin=0
	ipmax=0
	variablePgrid=(WaterWorld.or.setsurfpressure)
	do i=1,n_ret
		lowkey=RetPar(i)%keyword
		do j=1,len_trim(lowkey)
			if(iachar(lowkey(j:j)).ge.65.and.iachar(lowkey(j:j)).le.90) then
				lowkey(j:j)=achar(iachar(lowkey(j:j))+32)
			endif
		enddo
		if(lowkey.eq.'pmin') then
			variablePgrid=.true.
			ipmin=i
		endif
		if(lowkey.eq.'pmax') then
			variablePgrid=.true.
			ipmax=i
		endif
	enddo
	if(variablePgrid) then
		Pmin0=Pmin
		Pmax0=Pmax
		do i=1,nmodels
			call MapRetrievalMN(var(i,1:n_ret),error)
			if(ipmin.gt.0) then
				if(RetPar(ipmin)%value.lt.Pmin0) Pmin0=RetPar(ipmin)%value
			endif
			if(ipmax.gt.0) then
				if(RetPar(ipmax)%value.gt.Pmax0) Pmax0=RetPar(ipmax)%value
			endif
		enddo
		do i=1,nr
			Pgrid(i)=exp(log(Pmax0)+log(Pmin0/Pmax0)*real(i-1)/real(nr-1))
		enddo
	endif
	i1=nmodels/4
	i2=(3*nmodels)/4
	k=0
	do j=1,n_ret
		values(1:nmodels,j)=var(1:nmodels,j)
		call sort(values(1:nmodels,j),nmodels)
	enddo
15	continue
	do j=1,n_ret
		x1(j)=values(i1,j)
		x2(j)=values(i2,j)
	enddo
	lbest=0d0
	ibest=0
	j=0
	sorted(1:nmodels)=like(1:nmodels)
	call sort(sorted,nmodels)
	im1=real(nmodels)-real(nmodels)*0.682689492137086+0.5d0
	im2=real(nmodels)-real(nmodels)*0.954499736103642+0.5d0
	im3=real(nmodels)-real(nmodels)*0.997300203936740+0.5d0
	lm1=sorted(im1)

	do i=1,nmodels
		if(.not.(any(var(i,1:n_ret).lt.x1(1:n_ret)).or.any(var(i,1:n_ret).gt.x2(1:n_ret)))) then
			j=j+1
			if(like(i).gt.lbest) then
				ibest=i
				lbest=like(i)
			endif
		endif
	enddo
	if(j.ne.1) then
		if(j.gt.1.and.k.eq.0) then
			i1=i1+1
			i2=i2-1
			if(i1.lt.i2) goto 15
		else if(j.lt.1) then
			k=k+1
			i1=i1-1
			i2=i2+1
			if(i1.ge.1.and.i2.le.nmodels) goto 15
		endif
	endif
	call output("MPM model output written for all variable P > " // dbl2string(2d0*real(i1)/real(nmodels),'(f4.2)'))
	call MapRetrievalMN(var(ibest,1:n_ret),error)
	call system("cp " // trim(outputdir) // "input.dat " // trim(outputdir) // "mpm.dat")
	open(unit=21,file=trim(outputdir) // "mpm.dat",FORM="FORMATTED",access='APPEND')
	write(21,'("*** Median Probability Model parameters ***")')
	write(21,'("retrieval=.false.")')
	write(21,'("pew=.false.")')
	do i=1,n_ret
		write(21,'(a," = ",es14.7)') trim(RetPar(i)%keyword),RetPar(i)%value
	enddo
	if(free_tprofile.and..not.freePT_fitP) then
		do i=1,nTpoints
			write(21,'("Ppoint",i0.3," = ",es14.7)') i,Ppoint(i)
		enddo
	endif
	do i=1,n_add_ret
		write(21,'(a)') trim(line_add_ret(i))
	enddo
	close(unit=21)	

	done=.false.
	
c	call cpu_time(starttime)
	call SYSTEM_CLOCK(counts, count_rate, count_max)
	starttime = DBLE(counts)/DBLE(count_rate)

	if(npew.lt.0) then
		donmodels=nmodels
	else
		donmodels=min(npew,nmodels)
	endif

	open(unit=83,file=trim(outputdir) // "pew_output.dat",FORM="FORMATTED",ACCESS="STREAM")
	line="# "
	do i=1,n_ret
		if(i.eq.1) then
			write(line(3:12),'(a10)') RetPar(i)%keyword
		else
			write(line(i*12-11:i*12+1),'(" ",a11)') RetPar(i)%keyword
		endif
	enddo
	write(line(n_ret*12+1:6000),'(4a12)') "C/O","[Z]","Teff","MMW"
	write(83,'(a)') trim(line)

	do i=0,donmodels

	if(i.eq.0) then
		call output("best fit model")
	else
5		imodel=random(idum)*nmodels+1
		if(done(imodel)) goto 5
		done(imodel)=.true.

c		call cpu_time(stoptime)
		call SYSTEM_CLOCK(counts, count_rate, count_max)
		stoptime = DBLE(counts)/DBLE(count_rate)

		call output("model number: " // int2string(imodel,'(i7)') 
     &			// "(" // trim(dbl2string(real(i)*100d0/real(donmodels),'(f5.1)')) // "%)")
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
	do j=1,n_ret
		values(i,j)=RetPar(j)%value
	enddo
	COratio_der(i)=COratio
	Z_der(i)=metallicity
	MMW_der(i)=MMW(1)
	if(do3D) hotspotshift_der(i)=hotspotshift
	call SetOutputMode(.true.)

	call InitDens()
	if(retrievestar) then
		if(standardstarname.ne.' ') call StandardStar(standardstarname,Tstar,Rstar,Mstar)
		call StarSpecSetup(Tstar,logg,1d4*lam,1d4*blam,Fstar,nlam,starfile,blackbodystar)
		Fstar=Fstar*pi*Rstar**2
	endif
	call SetOutputMode(.false.)
	call ComputeModel(.true.)
	call SetOutputMode(.true.)
	
	if(i2d.ne.0.and.nobs.ne.0) then
		do iobs=1,nobs
			if(ObsSpec(iobs)%i2d.eq.i2d) then
				select case(ObsSpec(iobs)%type)
					case("trans","transmission","transC")
						spectrans(i,1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
					case("transM")
						if(do3D) then
							spectrans(i,1:nlam)=2d0*obsA(1:nlam,1)/(pi*Rstar**2)
						else
							spectrans(i,1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
						endif
					case("transE")
						if(do3D) then
							spectrans(i,1:nlam)=2d0*obsA(1:nlam,2)/(pi*Rstar**2)
						else
							spectrans(i,1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
						endif
					case("emisr","emisR","emisa","emis","emission")
						specemisR(i,1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar(1:nlam)*1d23/distance**2)
						specemis(i,1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
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
			specemisR(i,1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar(1:nlam)*1d23/distance**2)
			specemis(i,1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
		endif
	endif

	if(useobsgrid) then
		do iobs=1,nobs
			select case(ObsSpec(iobs)%type)
				case("trans","transmission","transC","transM","transE","emisr","emisR","emisa","emis","emission")
					call RemapObs(iobs,specobs(i,iobs,1:ObsSpec(iobs)%ndata),spectemp)
					ObsSpec(iobs)%scale=1d0
					if(ObsSpec(iobs)%scaling) then
						xy=0d0
						xx=0d0
						do ilam=1,ObsSpec(iobs)%ndata
							xy=xy+specobs(i,iobs,ilam)*(ObsSpec(iobs)%y(ilam)+ObsSpec(iobs)%offset)/ObsSpec(iobs)%dy(ilam)**2
							xx=xx+specobs(i,iobs,ilam)*specobs(i,iobs,ilam)/ObsSpec(iobs)%dy(ilam)**2
						enddo
						if(ObsSpec(iobs)%dscale.gt.0d0) then
							xx=xx+1d0/ObsSpec(iobs)%dscale**2
							xy=xy+1d0/ObsSpec(iobs)%dscale**2
						endif
						ObsSpec(iobs)%scale=1d0
						if(xx.gt.0d0) ObsSpec(iobs)%scale=xx/xy
					endif
					specobs(i,iobs,1:ObsSpec(iobs)%ndata)=specobs(i,iobs,1:ObsSpec(iobs)%ndata)/ObsSpec(iobs)%scale

					allocate(Cov(ObsSpec(iobs)%ndata,ObsSpec(iobs)%ndata))
					Cov(1:ObsSpec(iobs)%ndata,1:ObsSpec(iobs)%ndata)=0d0
					if(Cov_n_loc .gt. 0) then
						do j = 1, Cov_n_loc
							do ilam = 1, ObsSpec(iobs)%ndata
								xx = (ObsSpec(iobs)%lam(ilam)*1d4 - Cov_lam_loc(j)) / Cov_L_loc(j)
								spectemp(ilam) = Cov_a_loc(j)*exp(-0.5d0*xx**2)
							enddo
							call dger(ObsSpec(iobs)%ndata, ObsSpec(iobs)%ndata, 1d0, spectemp, 1, spectemp, 1, Cov, ObsSpec(iobs)%ndata)
						enddo
					endif
					do j=1,ObsSpec(iobs)%ndata
						Cov(j,j)=Cov(j,j)+(ObsSpec(iobs)%dy(j))**2
						do ilam=1,ObsSpec(iobs)%ndata
							xx=(ObsSpec(iobs)%lam(j)-ObsSpec(iobs)%lam(ilam))*1d4
							Cov(j,ilam)=Cov(j,ilam)+ObsSpec(iobs)%Cov_a**2*exp(-0.5*(xx/ObsSpec(iobs)%Cov_L)**2)
						enddo
					enddo
					call dpotrf('L', ObsSpec(iobs)%ndata, Cov, ObsSpec(iobs)%ndata, info)
      ! Fill spectemp with independent N(0,1) samples
					do ilam = 1, ObsSpec(iobs)%ndata
						spectemp(ilam) = gasdev(idum)
					enddo
      ! Multiply by Cholesky factor (L from DPOTRF in Cov)
					call dtrmv('L','N','N', ObsSpec(iobs)%ndata, Cov, ObsSpec(iobs)%ndata, spectemp, 1)
      ! Add to the deterministic model
					specobs(i,iobs,1:ObsSpec(iobs)%ndata)=specobs(i,iobs,1:ObsSpec(iobs)%ndata)
     &								+spectemp(1:ObsSpec(iobs)%ndata)
					deallocate(Cov)
			end select
		enddo
	endif

	if(imodel.gt.0) then
	if(like(imodel).gt.lm1.or.i.eq.1) then
		ctrans=maxval(spectrans(i,1:nlam))-minval(spectrans(i,1:nlam))
		if(ctrans.gt.cmax.or.i.eq.1) then
			cmax=ctrans
			call system("cp " // trim(outputdir) // "input.dat " // trim(outputdir) // "maxcontrast.dat")
			open(unit=21,file=trim(outputdir) // "maxcontrast.dat",FORM="FORMATTED",access='APPEND')
			write(21,'("*** Maximum contrast Model parameters ***")')
			write(21,'("retrieval=.false.")')
			write(21,'("pew=.false.")')
			do j=1,n_ret
				write(21,'(a," = ",es14.7)') trim(RetPar(j)%keyword),RetPar(j)%value
			enddo
			close(unit=21)	
		endif
		if(ctrans.lt.cmin.or.i.eq.1) then
			cmin=ctrans
			call system("cp " // trim(outputdir) // "input.dat " // trim(outputdir) // "mincontrast.dat")
			open(unit=21,file=trim(outputdir) // "mincontrast.dat",FORM="FORMATTED",access='APPEND')
			write(21,'("*** Minimum contrast Model parameters ***")')
			write(21,'("retrieval=.false.")')
			write(21,'("pew=.false.")')
			do j=1,n_ret
				write(21,'(a," = ",es14.7)') trim(RetPar(j)%keyword),RetPar(j)%value
			enddo
			close(unit=21)	
		endif
	endif	
	endif

	i2d=i2d+1
	if(i2d.le.n2d) goto 3

	if(variablePgrid) then
		Pg1(1:nr)=-log(P(1:nr))
		Pg2(1:nr)=-log(Pgrid(1:nr))

		yy(1:nr)=T(1:nr)
		call regridarray(Pg1,yy,nr,Pg2,T,nr)

		yy(1:nr)=log(Kzz_g(1:nr))
		call regridarray(Pg1,yy,nr,Pg2,Kzz_g,nr)
		Kzz_g(1:nr)=exp(Kzz_g(1:nr))

		do imol=1,nmol
			if(includemol(imol)) then
				yy(1:nr)=mixrat_r(1:nr,imol)
				call regridarray(Pg1,yy,nr,Pg2,mixrat_r(1:nr,imol),nr)
			endif
		enddo
		yy(1:nr)=cloud_dens(1:nr,1)
		call regridarray(Pg1,yy,nr,Pg2,cloud_dens(1:nr,1),nr)

		yy(1:nr)=log(dens(1:nr))
		call regridarray(Pg1,yy,nr,Pg2,dens,nr)
		dens(1:nr)=exp(dens(1:nr))

		if(do3D.and.fulloutput3D) then
			do j=0,nphase
				yy(1:nr)=PTaverage3D(j,1:nr)
				call regridarray(Pg1,yy,nr,Pg2,PTaverage3D(j,1:nr),nr)
				do imol=1,nmol
					if(includemol(imol)) then
						yy(1:nr)=mixrat_average3D(j,1:nr,imol)
						call regridarray(Pg1,yy,nr,Pg2,mixrat_average3D(j,1:nr,imol),nr)
					endif
				enddo
			enddo
		endif
		do j=1,nr
			if((Pgrid(j).lt.P(nr).and.ipmin.ne.0).or.(Pgrid(j).gt.P(1).and.(ipmax.ne.0.or.WaterWorld))) then
				mixrat_r(j,1:nmol)=1d-50
				cloud_dens(j,1)=0d0
				dens(j)=1d-50
				if(do3D.and.fulloutput3D) mixrat_average3D(0:nphase,j,1:nmol)=1d-50
			endif
		enddo
		P(1:nr)=Pgrid(1:nr)
	endif
	PTstruct(i,1:nr)=T(1:nr)
	dPTstruct(i,1)=log(T(2)/T(1))/log(P(2)/P(1))
	do j=1,nr
		Kzz_struct(i,j)=Kzz_g(j)
		if(dochemistry.or..true.) mol_struct(i,j,1:nmol)=mixrat_r(j,1:nmol)
	enddo

	do j=2,nr-1
		dPTstruct(i,j)=log(T(j+1)/T(j-1))/log(P(j+1)/P(j-1))
	enddo
	dPTstruct(i,nr)=log(T(nr)/T(nr-1))/log(P(nr)/P(nr-1))
	cloudstruct(i,1:nr)=cloud_dens(1:nr,1)
	cloudstruct(i,1:nr)=cloudstruct(i,1:nr)/dens(1:nr)
	do j=1,ncc
		speccloudtau(i,1:nlam)=speccloudtau(i,1:nlam)+cloudtau(j,1:nlam)
	enddo
	do j=1,nclouds
		if(Cloud(j)%onepart) then
			speccloudK(i,j,1:nlam)=Cloud(j)%Kext(nr/2,1:nlam)
		endif
	enddo
	Tplanet(i)=TeffPoutput
	if(do3D.and.fulloutput3D) then
		PTstruct3D(i,0:nphase,1:nr)=PTaverage3D(0:nphase,1:nr)
		dPTstruct3D(i,0:nphase,1)=log(PTaverage3D(0:nphase,2)/PTaverage3D(0:nphase,1))/log(P(2)/P(1))
		do j=2,nr-1
			dPTstruct3D(i,0:nphase,j)=log(PTaverage3D(0:nphase,j+1)/PTaverage3D(0:nphase,j-1))/log(P(j+1)/P(j-1))
		enddo
		dPTstruct3D(i,0:nphase,nr)=log(PTaverage3D(0:nphase,nr)/PTaverage3D(0:nphase,nr-1))/log(P(nr)/P(nr-1))
		mixrat3D(i,0:nphase,1:nr,1:nmol)=mixrat_average3D(0:nphase,1:nr,1:nmol)
		do iphase=1,nphase
			phase3D(i,iphase,1:nlam)=phase(iphase,0,1:nlam)+flux(0,1:nlam)
			phase3DR(i,iphase,1:nlam)=(phase(iphase,0,1:nlam)+flux(0,1:nlam))/(Fstar(1:nlam)*1d23/distance**2)
		enddo
		var3D(i,1:nlong,0)=beta3D_eq(1:nlong)
		do j=1,n_Par3D
			xmin=Par3D(j)%xmin
			if(Par3D(j)%multiply) then
				xmax=xmin*Par3D(j)%xmax
			else
				xmax=Par3D(j)%xmax
			endif			
			do ir=1,nlong-1
				if(Par3D(j)%logscale) then
					var3D(i,ir,j)=10d0**(log10(xmin)+
     &							  log10(xmax/xmin)*x3D_eq(ir))
				else
					var3D(i,ir,j)=xmin+(xmax-xmin)*x3D_eq(ir)
				endif
			enddo
		enddo
	endif
	if(i.ne.0) then
		write(line,'("(",i0.4,"es12.4)")') n_ret+4
		write(83,line) var(imodel,1:n_ret),COratio_der(i),Z_der(i),Tplanet(i),MMW_der(i)
	endif
	if(i.gt.2.and.(100*(i/100).eq.i.or.i.eq.donmodels.or.i.le.10)) then
		im1=real(i)/2d0-real(i)*0.682689492137086/2d0+0.5d0
		im2=real(i)/2d0-real(i)*0.954499736103642/2d0+0.5d0
		im3=real(i)/2d0-real(i)*0.997300203936740/2d0+0.5d0
		ip1=real(i)/2d0+real(i)*0.682689492137086/2d0+0.5d0
		ip2=real(i)/2d0+real(i)*0.954499736103642/2d0+0.5d0
		ip3=real(i)/2d0+real(i)*0.997300203936740/2d0+0.5d0
		ime=real(i)/2d0+0.5d0
		if(im1.lt.1) im1=1
		if(im2.lt.1) im2=1
		if(im3.lt.1) im3=1
		if(ime.lt.1) ime=1
		if(ip1.gt.i) ip1=i
		if(ip2.gt.i) ip2=i
		if(ip3.gt.i) ip3=i
		if(ime.gt.i) ime=i
		open(unit=26,file=trim(outputdir) // "trans_limits",FORM="FORMATTED",ACCESS="STREAM")
		open(unit=27,file=trim(outputdir) // "emis_limits",FORM="FORMATTED",ACCESS="STREAM")
		open(unit=28,file=trim(outputdir) // "emisR_limits",FORM="FORMATTED",ACCESS="STREAM")
		do ilam=1,nlam
		if(computelam(ilam)) then
			sorted(1:i)=spectrans(1:i,ilam)
			call sort(sorted,i)
			write(26,*) lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
			sorted(1:i)=specemis(1:i,ilam)
			call sort(sorted,i)
			write(27,*) lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
			sorted(1:i)=specemisR(1:i,ilam)
			call sort(sorted,i)
			write(28,*) lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
		endif
		enddo
		close(unit=26)
		close(unit=27)
		close(unit=28)
		
		if(useobsgrid) then
			do iobs=1,nobs
				open(unit=26,file=trim(outputdir) // "obs" // trim(int2string(iobs,'(i0.3)')) // "_limits",
     &						FORM="FORMATTED",ACCESS="STREAM")
				do ilam=1,ObsSpec(iobs)%ndata
					sorted(1:i)=specobs(1:i,iobs,ilam)
					call sort(sorted,i)
					write(26,*) ObsSpec(iobs)%lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
				enddo
			enddo
		endif
		
		open(unit=26,file=trim(outputdir) // "PT_limits",FORM="FORMATTED",ACCESS="STREAM")
		do ir=1,nr
			sorted(1:i)=PTstruct(1:i,ir)
			call sort(sorted,i)
			write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
		enddo
		close(unit=26)
		open(unit=26,file=trim(outputdir) // "dPT_limits",FORM="FORMATTED",ACCESS="STREAM")
		do ir=1,nr
			sorted(1:i)=dPTstruct(1:i,ir)
			call sort(sorted,i)
			write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
		enddo
		close(unit=26)
		open(unit=26,file=trim(outputdir) // "Kzz_limits",FORM="FORMATTED",ACCESS="STREAM")
		do ir=1,nr
			sorted(1:i)=Kzz_struct(1:i,ir)
			call sort(sorted,i)
			write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
		enddo
		close(unit=26)

		if(dochemistry.or..true.) then
		do imol=1,nmol
			if(includemol(imol)) then
				open(unit=26,file=trim(outputdir) // trim(molname(imol)) // "_limits",FORM="FORMATTED",ACCESS="STREAM")
				do ir=1,nr
					sorted(1:i)=mol_struct(1:i,ir,imol)
					call sort(sorted,i)
					write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
				enddo
				close(unit=26)
			endif
		enddo
		endif

		open(unit=26,file=trim(outputdir) // "cloud_limits",FORM="FORMATTED",ACCESS="STREAM")
		do ir=1,nr
			sorted(1:i)=cloudstruct(1:i,ir)
			call sort(sorted,i)
			write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
		enddo
		close(unit=26)

		open(unit=26,file=trim(outputdir) // "cloudtau_limits",FORM="FORMATTED",ACCESS="STREAM")
		do ilam=1,nlam
			if(computelam(ilam)) then
			sorted(1:i)=speccloudtau(1:i,ilam)
			call sort(sorted,i)
			write(26,*) lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
			endif
		enddo
		close(unit=26)

		do j=1,nclouds
			if(Cloud(j)%onepart) then
			open(unit=26,file=trim(outputdir) // "cloudK" // trim(int2string(j,'(i0.3)')) //"_limits",FORM="FORMATTED",ACCESS="STREAM")
			do ilam=1,nlam
				if(computelam(ilam)) then
				sorted(1:i)=speccloudK(1:i,j,ilam)
				call sort(sorted,i)
				write(26,*) lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
				endif
			enddo
			close(unit=26)
			endif
		enddo

		if(do3D.and.fulloutput3D) then
			open(unit=26,file=trim(outputdir) // "beta3D_eq_limits",FORM="FORMATTED",ACCESS="STREAM")
			do ir=1,nlong-1
				sorted(1:i)=var3D(1:i,ir,0)
				call sort(sorted,i)
				write(26,*) 90d0*(long(ir)+long(ir+1))/pi,sorted(im3),sorted(im2),sorted(im1),
     &								sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
			enddo
			close(unit=26)
			do j=1,n_Par3D
				open(unit=26,file=trim(outputdir) // trim(Par3D(j)%keyword) // "_eq_limits",FORM="FORMATTED",ACCESS="STREAM")
				do ir=1,nlong-1
					sorted(1:i)=var3D(1:i,ir,j)
					call sort(sorted,i)
					write(26,*) 90d0*(long(ir)+long(ir+1))/pi,sorted(im3),sorted(im2),sorted(im1),
     &								sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
				enddo
				close(unit=26)
			enddo
			open(unit=26,file=trim(outputdir) // "PT_trans_limits",FORM="FORMATTED",ACCESS="STREAM")
			do ir=1,nr
				sorted(1:i)=PTstruct3D(1:i,0,ir)
				call sort(sorted,i)
				write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
			enddo
			close(unit=26)
			open(unit=26,file=trim(outputdir) // "dPT_trans_limits",FORM="FORMATTED",ACCESS="STREAM")
			do ir=1,nr
				sorted(1:i)=dPTstruct3D(1:i,0,ir)
				call sort(sorted,i)
				write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
			enddo
			close(unit=26)
			do imol=1,nmol
				if(includemol(imol)) then
					open(unit=26,file=trim(outputdir) // trim(molname(imol)) // "_trans_limits",FORM="FORMATTED",ACCESS="STREAM")
					do ir=1,nr
						sorted(1:i)=mixrat3D(1:i,0,ir,imol)
						call sort(sorted,i)
						write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
					enddo
					close(unit=26)
				endif
			enddo
			do iphase=1,nphase
				open(unit=26,file=trim(outputdir) // "PT" // trim(int2string(int(theta_phase(iphase)),'(i0.3)')) // "_limits",FORM="FORMATTED",ACCESS="STREAM")
				do ir=1,nr
					sorted(1:i)=PTstruct3D(1:i,iphase,ir)
					call sort(sorted,i)
					write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
				enddo
				close(unit=26)
				open(unit=26,file=trim(outputdir) // "dPT" // trim(int2string(int(theta_phase(iphase)),'(i0.3)')) // "_limits",FORM="FORMATTED",ACCESS="STREAM")
				do ir=1,nr
					sorted(1:i)=dPTstruct3D(1:i,iphase,ir)
					call sort(sorted,i)
					write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
				enddo
				close(unit=26)
				do imol=1,nmol
					if(includemol(imol)) then
						open(unit=26,file=trim(outputdir) // trim(molname(imol)) // "_" // 
     &						trim(int2string(int(theta_phase(iphase)),'(i0.3)')) // "_limits",FORM="FORMATTED",ACCESS="STREAM")
						do ir=1,nr
							sorted(1:i)=mixrat3D(1:i,iphase,ir,imol)
							call sort(sorted,i)
							write(26,*) P(ir),sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
						enddo
						close(unit=26)
					endif
				enddo
				open(unit=26,file=trim(outputdir) // "phase" // trim(int2string(int(theta_phase(iphase)),'(i0.3)')) 
     &			// "_limits",FORM="FORMATTED",ACCESS="STREAM")
				open(unit=27,file=trim(outputdir) // "phaseR" // trim(int2string(int(theta_phase(iphase)),'(i0.3)')) 
     &			// "_limits",FORM="FORMATTED",ACCESS="STREAM")
				do ilam=1,nlam
				if(computelam(ilam)) then
					sorted(1:i)=phase3D(1:i,iphase,ilam)
					call sort(sorted,i)
					write(26,*) lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
					sorted(1:i)=phase3DR(1:i,iphase,ilam)
					call sort(sorted,i)
					write(27,*) lam(ilam)*1d4,sorted(im3),sorted(im2),sorted(im1),sorted(ime),sorted(ip1),sorted(ip2),sorted(ip3)
				endif
				enddo
				close(unit=26)
				close(unit=27)
			enddo
		endif

		open(unit=26,file=trim(outputdir) // "retrieval",FORM="FORMATTED",ACCESS="STREAM")
		write(26,'("# ",a8,3a12,"  (",2a12")")') "par","median","1sig lower","1sig upper","3sig lower","3sig upper"
		do j=1,n_ret
			sorted(1:i)=values(1:i,j)
			call sort(sorted,i)
			write(26,'(a10,3es12.4,"  (",2es12.4")")') trim(RetPar(j)%keyword),sorted(ime),sorted(im1),sorted(ip1),sorted(im3),sorted(ip3)
		enddo
		sorted(1:i)=COratio_der(1:i)
		call sort(sorted,i)
		write(26,'(a10,3es12.4,"  (",2es12.4")")') "C/O",sorted(ime),sorted(im1),sorted(ip1),sorted(im3),sorted(ip3)
		sorted(1:i)=Z_der(1:i)
		call sort(sorted,i)
		write(26,'(a10,3es12.4,"  (",2es12.4")")') "[Z]",sorted(ime),sorted(im1),sorted(ip1),sorted(im3),sorted(ip3)
		if(do3D) then
			sorted(1:i)=hotspotshift_der(1:i)
			call sort(sorted,i)
			write(26,'(a10,3es12.4,"  (",2es12.4")")') "hotspot",sorted(ime),sorted(im1),sorted(ip1),sorted(im3),sorted(ip3)
		endif
		if(emisspec.and..not.useobsgrid) then
			sorted(1:i)=Tplanet(1:i)
			call sort(sorted,i)
			write(26,'(a10,3es12.4,"  (",2es12.4")")') "Tplanet",sorted(ime),sorted(im1),sorted(ip1),sorted(im3),sorted(ip3)
		endif
		sorted(1:i)=MMW_der(1:i)
		call sort(sorted,i)
		write(26,'(a10,3es12.4,"  (",2es12.4")")') "MMW",sorted(ime),sorted(im1),sorted(ip1),sorted(im3),sorted(ip3)

		close(unit=26)

		open(unit=26,file=trim(outputdir) // "trans_sigma",FORM="FORMATTED",ACCESS="STREAM")
		open(unit=27,file=trim(outputdir) // "emis_sigma",FORM="FORMATTED",ACCESS="STREAM")
		open(unit=28,file=trim(outputdir) // "emisR_sigma",FORM="FORMATTED",ACCESS="STREAM")
		do ilam=1,nlam
		if(computelam(ilam)) then
			sorted(1:i)=spectrans(1:i,ilam)
			call stats(sorted,i,aver,sig)
			write(26,*) lam(ilam)*1d4,spectrans(0,ilam),sig,aver
			sorted(1:i)=specemis(1:i,ilam)
			call stats(sorted,i,aver,sig)
			write(27,*) lam(ilam)*1d4,specemis(0,ilam),sig,aver
			sorted(1:i)=specemisR(1:i,ilam)
			call stats(sorted,i,aver,sig)
			write(28,*) lam(ilam)*1d4,specemisR(0,ilam),sig,aver
		endif
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
	deallocate(speccloudtau)
	deallocate(PTstruct)
	deallocate(Kzz_struct)
	deallocate(like)
	deallocate(sorted)
	deallocate(values)
	deallocate(COratio_der)
	deallocate(Z_der)
	deallocate(MMW_der)
	if(do3D) deallocate(hotspotshift_der)
	deallocate(done)
	deallocate(var)
	if(do3D.and.fulloutput3D) deallocate(PTstruct3D,mixrat3D,phase3D,phase3DR)

	writefiles=.true.
		
	return
	end
	