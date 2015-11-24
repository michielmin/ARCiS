	subroutine ReadObs()
	use GlobalSetup
	IMPLICIT NONE
	integer i,j,ilam,nj
	real*8 x,y,dy
	
	do i=1,nobs
		if(ObsSpec(i)%type.eq.'tprofile'.or.ObsSpec(i)%type.eq.'logtp') then
			ObsSpec(i)%nlam=nr
			allocate(ObsSpec(i)%lam(ObsSpec(i)%nlam))
			allocate(ObsSpec(i)%y(ObsSpec(i)%nlam))
			allocate(ObsSpec(i)%dy(ObsSpec(i)%nlam))
			ObsSpec(i)%y=0d0
			ObsSpec(i)%dy=1d0
			ObsSpec(i)%spec=.false.
		else
			open(unit=20,file=ObsSpec(i)%file,RECL=1000)
			j=1
			ilam=1
1			read(20,*,end=2,err=1) x,y,dy
			if(x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)) ilam=ilam+1
			j=j+1
			goto 1
2			ObsSpec(i)%nlam=ilam-1
			nj=j-1
			rewind(20)
			allocate(ObsSpec(i)%lam(ObsSpec(i)%nlam))
			allocate(ObsSpec(i)%y(ObsSpec(i)%nlam))
			allocate(ObsSpec(i)%dy(ObsSpec(i)%nlam))
			ilam=1
			do j=1,nj
3				read(20,*,err=3) x,y,dy
				if(x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)) then
					ObsSpec(i)%lam(ilam)=x*1d-4
					ObsSpec(i)%y(ilam)=y
					ObsSpec(i)%dy(ilam)=dy
					ilam=ilam+1
				endif
			enddo
		endif
	enddo
	
	return
	end



	subroutine DoRetrieval()
	use GlobalSetup
	IMPLICIT NONE
	external ComputeChi2
	real*8 var0(n_ret),dvar0(2,n_ret),ComputeChi2,dvar(n_ret),x(n_ret),chi2min
	real*8,allocatable :: y(:),y1(:),y2(:),dy(:,:),yobs(:),dyobs(:),ybest(:)
	real*8 chi2obs(nobs),var(n_ret),chi2,gasdev,maxd,error(2,n_ret),var_best(n_ret)
	real*8 x1,x2,minT(nr),maxT(nr),ran1,tot,lambda,chi2_1,chi2_2,dchi2(n_ret),chi2prev
	integer imodel,ny,i,j,iter1,iter2,iy,k
	
	real*8 dvar_av(n_ret),Cov(n_ret,n_ret),b(n_ret),W(n_ret,n_ret),Winv(n_ret,n_ret),dmax,scale
	real*8 chi2_spec,chi2_prof
	integer na,map(n_ret),info
	logical dofit(n_ret),initfit
	
	do i=1,n_ret
10		var0(i)=gasdev(idum)*0.1+0.5
		if(var0(i).gt.1d0) goto 10
		if(var0(i).lt.0d0) goto 10
	enddo
	call MapTprofile(var0)

	dvar0=10d0
	open(unit=31,file=trim(outputdir) // "Wolk.dat",RECL=6000)

	if(ngen.gt.0) then
c first genetic algoritm to make the first estimate
		call Genetic(ComputeChi2,var0,dvar0,n_ret,nobs,npop,ngen,idum,gene_cross,.true.)
	endif

	imodel=npop*ngen
	ny=0
	do i=1,nobs
		ny=ny+ObsSpec(i)%nlam
	enddo
	allocate(y(ny))
	allocate(ybest(ny))
	allocate(y1(ny))
	allocate(y2(ny))
	allocate(dy(n_ret,ny))
	allocate(yobs(ny))
	allocate(dyobs(ny))
	iy=1
	do i=1,nobs
		yobs(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%y(1:ObsSpec(i)%nlam)
		dyobs(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%dy(1:ObsSpec(i)%nlam)
		iy=iy+ObsSpec(i)%nlam
	enddo
	dvar=0.1d0
	var=var0
	var_best=var0
	chi2min=1d200
	chi2prev=1d200
	error=0.1d0
	lambda=0.1

	do iter1=1,10
	lambda=0.1
	var=var_best
	var0=var
	do iter2=1,10
		j=0
		imodel=imodel+1
		call output("Iteration " // trim(int2string(iter2+10*(iter1-1),'(i4)')) // 
     &				" - model" // trim(int2string(imodel,'(i5)')))
		chi2=1d0/ComputeChi2(imodel,n_ret,var,nobs,chi2obs,.true.,error)
		iy=1
		do i=1,nobs
	 		call RemapObs(i,y(iy:iy+ObsSpec(i)%nlam-1))
			iy=iy+ObsSpec(i)%nlam
		enddo
		if(chi2.le.chi2prev) then
			lambda=lambda/2d0
			var0=var
			if(chi2.le.chi2min) then
				var_best=var
				ybest=y
				chi2_spec=0d0
				chi2_prof=0d0
				do j=1,nobs
					if(ObsSpec(j)%spec) then
						chi2_spec=chi2_spec+1d0/chi2obs(j)
					else
						chi2_prof=chi2_prof+1d0/chi2obs(j)
					endif
				enddo
				chi2min=chi2
				call output("Iteration improved" // trim(dbl2string(chi2,'(f8.3)')))
				call SetOutputMode(.false.)
				call WriteStructure()
				call WriteOutput()
				call WriteRetrieval(imodel,chi2_spec,var(1:n_ret))
				call WritePTlimits(var0,Cov,error)
				call SetOutputMode(.true.)
			endif
		else if(iter2.eq.1) then
			var0=var
		else
			lambda=lambda*2d0
			var=var0
			goto 1
		endif
		do i=1,n_ret
			dvar(i)=max(error(1,i),error(2,i))/2d0
     	enddo
		chi2prev=chi2
		do i=1,n_ret
			if(dvar(i).lt.1d-3) dvar(i)=1d-3
			if(dvar(i).gt.1d-1) dvar(i)=1d-1
20			var=var0
			var(i)=var(i)+dvar(i)
			if(var(i).gt.1d0) var(i)=1d0
			if(var(i).lt.0d0) var(i)=0d0
			x1=var(i)
			imodel=imodel+1
			chi2=1d0/ComputeChi2(imodel,n_ret,var,nobs,chi2obs,RetPar(i)%opacitycomp,error)
			iy=1
			do j=1,nobs
		 		call RemapObs(j,y1(iy:iy+ObsSpec(j)%nlam-1))
				iy=iy+ObsSpec(j)%nlam
			enddo
			chi2_1=chi2
			if(chi2.lt.chi2min.and.RetPar(i)%opacitycomp) then
				chi2min=chi2
				var_best=var
				ybest=y
				chi2_spec=0d0
				chi2_prof=0d0
				do j=1,nobs
					if(ObsSpec(j)%spec) then
						chi2_spec=chi2_spec+1d0/chi2obs(j)
					else
						chi2_prof=chi2_prof+1d0/chi2obs(j)
					endif
				enddo
				call output("Updating best fit " // trim(dbl2string(chi2,'(f8.3)')))
				call SetOutputMode(.false.)
				call WriteStructure()
				call WriteOutput()
				call WriteRetrieval(imodel,chi2_spec,var(1:n_ret))
				call WritePTlimits(var,Cov,error)
				call SetOutputMode(.true.)
			endif
			var=var0
			var(i)=var(i)-dvar(i)
			if(var(i).gt.1d0) var(i)=1d0
			if(var(i).lt.0d0) var(i)=0d0
			x2=var(i)
			imodel=imodel+1
			chi2=1d0/ComputeChi2(imodel,n_ret,var,nobs,chi2obs,RetPar(i)%opacitycomp,error)
			iy=1
			do j=1,nobs
		 		call RemapObs(j,y2(iy:iy+ObsSpec(j)%nlam-1))
				iy=iy+ObsSpec(j)%nlam
			enddo
			chi2_2=chi2
			if(chi2.lt.chi2min.and.RetPar(i)%opacitycomp) then
				chi2min=chi2
				var_best=var
				ybest=y
				chi2_spec=0d0
				chi2_prof=0d0
				do j=1,nobs
					if(ObsSpec(j)%spec) then
						chi2_spec=chi2_spec+1d0/chi2obs(j)
					else
						chi2_prof=chi2_prof+1d0/chi2obs(j)
					endif
				enddo
				call output("Updating best fit " // trim(dbl2string(chi2,'(f8.3)')))
				call SetOutputMode(.false.)
				call WriteStructure()
				call WriteOutput()
				call WriteRetrieval(imodel,chi2_spec,var(1:n_ret))
				call WritePTlimits(var,Cov,error)
				call SetOutputMode(.true.)
			endif
			dy(i,1:ny)=(y1(1:ny)-y2(1:ny))/abs(x1-x2)
			dchi2(i)=(chi2_1-chi2_2)/abs(x1-x2)
			if(abs(chi2_1-chi2_2)/(chi2_1+chi2_2).lt.1d-4) then
				if(dvar(i).lt.0.99d0) then
					dvar(i)=dvar(i)*5d0
					goto 20
				else if(dvar(i).lt.0.999d0) then
					dvar(i)=1d0
					goto 20
				endif
			endif
		enddo


1		continue

		dofit=.true.
		initfit=.true.

2		continue
		na=0
		do i=1,n_ret
			if(dofit(i)) then
				na=na+1
				map(na)=i
			endif
		enddo
		if(na.gt.1) then
		
		W=0d0
		b=0d0
		do k=1,ny
			do i=1,na
				do j=1,na
					W(i,j)=W(i,j)+dy(map(i),k)*dy(map(j),k)/dyobs(k)**2
				enddo
				b(i)=b(i)+(yobs(k)-y(k))*dy(map(i),k)/dyobs(k)**2
			enddo
		enddo

		call MatrixInvert(W(1:na,1:na),Winv(1:na,1:na),na,info)
		if(INFO.ne.0) then
			call output('Error in matrix inversion.')
			call output('Retrieving parameters with no influence on the observations?')
			stop
		endif
		if(initfit) then
			Cov=Winv
			do i=1,n_ret
				do j=1,2
					error(j,i)=sqrt(abs(Winv(i,i)))
				enddo
			enddo
		endif

		do i=1,na
			W(i,i)=W(i,i)*(1d0+lambda)
		enddo
		call MatrixInvert(W(1:na,1:na),Winv(1:na,1:na),na,info)
		if(INFO.ne.0) then
			call output('Error in matrix inversion.')
			call output('Retrieving parameters with no influence on the observations?')
			stop
		endif

		dvar=0d0
		do i=1,na
			do j=1,na
				dvar(map(i))=dvar(map(i))+b(j)*Winv(i,j)
			enddo
			var(map(i))=var0(map(i))+dvar(map(i))
		enddo

		dmax=0d0
		j=0
		do i=1,n_ret
			if(dofit(i)) then
			if(var(i).gt.1d0.or.var(i).lt.0d0) then
				if(abs(dvar(i)/error(2,i)).gt.dmax) then
					dmax=abs(dvar(i)/error(2,i))
					j=i
				endif
				if(abs(dvar(i)/error(1,i)).gt.dmax) then
					dmax=abs(dvar(i)/error(1,i))
					j=i
				endif
			endif
			endif
		enddo
		if(j.ne.0) then
			dofit(j)=.false.
			if(var(j).gt.1d0) then
				var(j)=var0(j)
			endif
			if(var(j).lt.0d0) then
				var(j)=var0(j)
			endif
			initfit=.false.
			goto 2
		endif
		endif
	enddo
	do i=1,nobs
		if(.not.ObsSpec(i)%spec.and.ObsSpec(i)%scale.gt.0d0) then
			print*,chi2_spec,chi2_prof,chi2min
			ObsSpec(i)%beta=log10(ObsSpec(i)%scale*ObsSpec(i)%beta*chi2_spec/chi2_prof)+log10(ObsSpec(i)%beta)
			ObsSpec(i)%beta=10d0**(ObsSpec(i)%beta/2d0)
			call output("Adjusting beta to " // trim(dbl2string(ObsSpec(i)%beta,'(es10.4)')))
			chi2min=1d200
		endif
	enddo
	enddo

	close(unit=31)
	
	return
	end



	subroutine WritePTlimits(var0,W,error)
	use GlobalSetup
	IMPLICIT NONE
	real*8 minT(nr),maxT(nr),var(n_ret),error(2,n_ret),tot,W(n_ret,n_ret)
	real*8 var0(n_ret),random,vec(n_ret),Tbest(nr),gasdev
	integer i,j,k,minC(nr),maxC(nr)

	maxT=0d0
	minT=0d0
	maxC=0
	minC=0
	call SetOutputMode(.false.)
	var=var0
	call MapRetrieval(var,error)
	call InitDens()
	call SetupStructure()
	Tbest=T
	do k=1,1000
		do i=1,n_ret
			if(random(idum).gt.0.5d0) then
				var(i)=var0(i)-abs(gasdev(idum))*error(1,i)
			else
				var(i)=var0(i)+abs(gasdev(idum))*error(2,i)
			endif
			if(var(i).lt.0d0) var(i)=0d0
			if(var(i).gt.1d0) var(i)=1d0
		enddo
		call MapRetrieval(var,error)
		call InitDens()
		call SetupStructure()
		do i=1,nr
			if(T(i).gt.Tbest(i)) then
				maxT(i)=maxT(i)+(T(i)-Tbest(i))**2
				maxC(i)=maxC(i)+1
			else
				minT(i)=minT(i)+(T(i)-Tbest(i))**2
				minC(i)=minC(i)+1
			endif
		enddo				
	enddo
	maxT=sqrt(maxT/real(maxC))
	minT=sqrt(minT/real(minC))

	call SetOutputMode(.true.)
	open(unit=45,file=trim(outputdir) // "limits.dat")
	do i=1,nr
		write(45,*) P(i),Tbest(i)-minT(i),Tbest(i)+maxT(i)
	enddo
	close(unit=45)

	return
	end



	subroutine RemapObs(i,spec)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,k
	real*8 spec(*),x
	real*8 lamobs(nlam-1)
	real*8 eta,Tirr,tau,expint

	do j=1,nlam-1
		lamobs(j)=sqrt(lam(j)*lam(j+1))
	enddo
	select case(ObsSpec(i)%type)
		case("trans","transmission")
			call regridarray(lamobs,obsA(0,1:nlam-1)/(pi*Rstar**2),nlam-1,
     &					ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
		case("emisr","emisR")
			call regridarray(lamobs,flux(0,1:nlam-1)/(Fstar(1:nlam-1)*1d23/distance**2),
     &					nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
		case("emisa","emis","emission")
			call regridarray(lamobs,flux(0,1:nlam-1),
     &					nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
		case("tprofile")
			tau=0d0
			Tirr=betaT*sqrt(Rstar/(2d0*Dplanet))*Tstar
			do j=nr,1,-1
				eta=2d0/3d0+(2d0/(3d0*gammaT))*(1d0+(gammaT*tau/2d0-1)*exp(-gammaT*tau))
     &					+(2d0*gammaT/3d0)*(1d0-tau**2/2d0)*expint(2,gammaT*tau)
				spec(j)=(3d0*TeffP**4/4d0)*(2d0/3d0+tau)+(3d0*Tirr**4/4d0)*eta
				spec(j)=spec(j)**0.25d0
				spec(j)=(spec(j)-T(j))/spec(j)
				tau=tau+kappaT*dens(j)*(R(j+1)-R(j))
			enddo
			spec(1:ObsSpec(i)%nlam)=spec(1:ObsSpec(i)%nlam)*ObsSpec(i)%beta
		case("logtp")
			spec(1)=(log10(T(1)/T(2))/log10(P(1)/P(2)))
			do j=2,nr-1
				spec(j)=(log10(T(j-1)/T(j))/log10(P(j-1)/P(j)))-(log10(T(j)/T(j+1))/log10(P(j)/P(j+1)))
				spec(j)=spec(j)+log10(T(j-1)/T(j+1))/log10(P(j-1)/P(j+1))
			enddo
			spec(nr)=(log10(T(nr-1)/T(nr))/log10(P(nr-1)/P(nr)))
			spec(1:ObsSpec(i)%nlam)=spec(1:ObsSpec(i)%nlam)/real(nr)
			spec(1:ObsSpec(i)%nlam)=spec(1:ObsSpec(i)%nlam)*ObsSpec(i)%beta
	end select

	return
	end
		
	
	real*8 function ComputeChi2(imodel,nvars,var,nobs0,chi2obs,recomputeopac,error0)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nvars,nobs0,i,j,imodel
	real*8 var(nvars),chi2obs(nobs0),error(2,n_ret)
	real*8,allocatable :: spec(:)
	real*8,intent(in),optional :: error0(2,n_ret)
	logical recomputeopac
	error=0d0
	if(present(error0)) error=error0

	call MapRetrieval(var,error)
	
	call InitDens()
	call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam)
	Fstar=Fstar*pi*Rstar**2*pi/3.336e11
	call SetOutputMode(.false.)
	call ComputeModel(recomputeopac)
	call SetOutputMode(.true.)
	
	ComputeChi2=0d0
	do i=1,nobs
		allocate(spec(ObsSpec(i)%nlam))
		call RemapObs(i,spec)
		chi2obs(i)=0d0
		do j=1,ObsSpec(i)%nlam
			chi2obs(i)=chi2obs(i)+((spec(j)-ObsSpec(i)%y(j))/ObsSpec(i)%dy(j))**2
		enddo
		chi2obs(i)=chi2obs(i)/real(ObsSpec(i)%nlam)
		ComputeChi2=ComputeChi2+chi2obs(i)
		deallocate(spec)
	enddo

	write(31,*) imodel,ComputeChi2,(trim(dbl2string(RetPar(i)%value,'(es18.7)')),i=1,n_ret)
	flush(31)

	chi2obs=1d0/chi2obs
	ComputeChi2=1d0/ComputeChi2
	
	return
	end
	
	
	subroutine MapRetrieval(var,dvar)
	use GlobalSetup
	use ReadKeywords
	use Constants
	IMPLICIT NONE
	type(SettingKey) key
	character*1000 readline
	integer i
	real*8 var(n_ret),dvar(2,n_ret),x

	do i=1,n_ret
		if(RetPar(i)%logscale) then
c	log
			x=var(i)
			RetPar(i)%value=10d0**(log10(RetPar(i)%xmin)+log10(RetPar(i)%xmax/RetPar(i)%xmin)*x)
			x=var(i)+dvar(1,i)
			RetPar(i)%error1=10d0**(log10(RetPar(i)%xmin)+log10(RetPar(i)%xmax/RetPar(i)%xmin)*x)
			RetPar(i)%error1=RetPar(i)%error1/RetPar(i)%value
			x=var(i)+dvar(2,i)
			RetPar(i)%error2=10d0**(log10(RetPar(i)%xmin)+log10(RetPar(i)%xmax/RetPar(i)%xmin)*x)
			RetPar(i)%error2=RetPar(i)%error2/RetPar(i)%value
		else if(RetPar(i)%squarescale) then
c	square
			x=var(i)
			RetPar(i)%value=sqrt(RetPar(i)%xmin**2+(RetPar(i)%xmax**2-RetPar(i)%xmin**2)*x)
			x=var(i)+dvar(1,i)
			RetPar(i)%error1=sqrt(RetPar(i)%xmin**2+(RetPar(i)%xmax**2-RetPar(i)%xmin**2)*x)
			RetPar(i)%error1=RetPar(i)%error1-RetPar(i)%value
			x=var(i)+dvar(2,i)
			RetPar(i)%error2=sqrt(RetPar(i)%xmin**2+(RetPar(i)%xmax**2-RetPar(i)%xmin**2)*x)
			RetPar(i)%error2=RetPar(i)%error2-RetPar(i)%value
		else
c	linear
			x=var(i)
			RetPar(i)%value=RetPar(i)%xmin+(RetPar(i)%xmax-RetPar(i)%xmin)*x
			x=var(i)+dvar(1,i)
			RetPar(i)%error1=RetPar(i)%xmin+(RetPar(i)%xmax-RetPar(i)%xmin)*x
			RetPar(i)%error1=RetPar(i)%error1-RetPar(i)%value
			x=var(i)+dvar(2,i)
			RetPar(i)%error2=RetPar(i)%xmin+(RetPar(i)%xmax-RetPar(i)%xmin)*x
			RetPar(i)%error2=RetPar(i)%error2-RetPar(i)%value
		endif
	enddo

	Rplanet=Rplanet/Rjup
	Mplanet=Mplanet/Mjup
	Rstar=Rstar/Rsun
	Dplanet=Dplanet/AU
	lam1=lam1/micron
	lam2=lam2/micron
	distance=distance/parsec
	do i=1,n_ret
		readline=trim(RetPar(i)%keyword) // "=" // trim(dbl2string(RetPar(i)%value,'(es14.7)'))
		call get_key_value(readline,key%key,key%key1,key%key2,key%value,key%nr1,key%nr2)
		call ReadAndSetKey(key)
	enddo
	call ConvertUnits()

	return
	end
	

	subroutine MapTprofile(var)
	use GlobalSetup
	use ReadKeywords
	use Constants
	IMPLICIT NONE
	integer i,j
	real*8 var(n_ret),x,Tmap(nr),tau,eta,Tirr,expint

	call InitDens()
	call SetupStructure()

	tau=0d0
	do i=nr,1,-1
		Tirr=betaT*sqrt(Rstar/(2d0*Dplanet))*Tstar
		eta=2d0/3d0+(2d0/(3d0*gammaT))*(1d0+(gammaT*tau/2d0-1)*exp(-gammaT*tau))
     &					+(2d0*gammaT/3d0)*(1d0-tau**2/2d0)*expint(2,gammaT*tau)
		Tmap(i)=(3d0*TeffP**4/4d0)*(2d0/3d0+tau)+(3d0*Tirr**4/4d0)*eta
		Tmap(i)=Tmap(i)**0.25d0
		tau=tau+kappaT*dens(i)*(R(i+1)-R(i))
	enddo
	do i=1,n_ret
		if(RetPar(i)%keyword(1:6).eq.'tvalue') then
			read(RetPar(i)%keyword(7:len_trim(RetPar(i)%keyword)),*) j
			if(RetPar(i)%logscale) then
c	log
				var(i)=(log10(Tmap(j))-log10(RetPar(i)%xmin))/log10(RetPar(i)%xmax/RetPar(i)%xmin)
			else if(RetPar(i)%squarescale) then
c	square
				var(i)=(Tmap(j)**2-RetPar(i)%xmin**2)/(RetPar(i)%xmax**2-RetPar(i)%xmin**2)
			else
c	linear
				var(i)=(Tmap(j)-RetPar(i)%xmin)/(RetPar(i)%xmax-RetPar(i)%xmin)
			endif
			if(var(i).gt.1d0) var(i)=1d0
			if(var(i).lt.0d0) var(i)=0d0
		endif
	enddo

	return
	end
	

	
	subroutine regridarray(x0,y0,n0,x1,y1,n1)
	IMPLICIT NONE
	integer i1,i0,n0,n1
	real*8 x0(n0),y0(n0),x1(n1),y1(n1)

	do i1=1,n1
		if(x1(i1).lt.x0(1)) then
			y1(i1)=y0(1)
		else if(x1(i1).gt.x0(n0)) then
			y1(i1)=y0(n0)
		else
			do i0=1,n0-1
				if(x1(i1).ge.x0(i0).and.x1(i1).le.x0(i0+1)) then
					y1(i1)=y0(i0)+(y0(i0+1)-y0(i0))*(x0(i0+1)-x1(i1))/(x0(i0+1)-x0(i0))
				endif
			enddo
		endif
	enddo

	return
	end
	
	subroutine WriteRetrieval(imodel,chi2,var)
	use GlobalSetup
	IMPLICIT NONE
	real*8 var(n_ret),chi2
	integer i,imodel

	open(unit=20,file=trim(outputdir) // "retrieval",RECL=1000)
	write(20,'("Model ",i)') imodel
	write(20,'("chi2=",f14.6)') chi2
	do i=1,n_ret
		if(RetPar(i)%logscale) then
c	log
			write(20,'(a15," = ",es14.7," x/: ",es11.4,es11.4)') trim(RetPar(i)%keyword),RetPar(i)%value,
     &					RetPar(i)%error2,RetPar(i)%error1
 		else
c	linear/squared
			write(20,'(a15," = ",es14.7," +/- ",es11.4,es11.4)') trim(RetPar(i)%keyword),RetPar(i)%value,
     &					RetPar(i)%error2,RetPar(i)%error1
		endif
	enddo
	close(unit=20)
		
		
	return
	end


	subroutine MatrixInvert(A,Ainv,N,INFO)
	IMPLICIT NONE
	integer n,ierr
	real*8 A(n,n),Ainv(n,n)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK
      INTEGER :: i,j, LWORK

      INTEGER	     INFO, LDA,	M
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV

      INTEGER DeAllocateStatus

      external DGETRF
      external DGETRI

	INFO=0
      LDA = N
      LWORK = N*N
      ALLOCATE (WORK(LWORK))
      ALLOCATE (IPIV(N))
C
C     DGETRF computes an LU factorization of a general M-by-N matrix A
C     using partial pivoting with row interchanges.

      M=N
      LDA=N

C  Store A in Ainv to prevent it from being overwritten by LAPACK

      Ainv = A

      CALL DGETRF( M, N, Ainv, LDA, IPIV, INFO )

C  DGETRI computes the inverse of a matrix using the LU factorization
C  computed by DGETRF.

      CALL DGETRI(N, Ainv, N, IPIV, WORK, LWORK, INFO)

      DEALLOCATE (IPIV, STAT = DeAllocateStatus)
      DEALLOCATE (WORK, STAT = DeAllocateStatus)

	return
	end
	