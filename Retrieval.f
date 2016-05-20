	subroutine ReadObs()
	use GlobalSetup
	IMPLICIT NONE
	integer i,j,ilam,nj
	real*8 x,y,dy,specres_obs,expspecres_obs
	character*1000 line
	
	do i=1,nobs
		select case(ObsSpec(i)%type)
			case('tprofile','logtp')
				ObsSpec(i)%nlam=nr
				allocate(ObsSpec(i)%lam(ObsSpec(i)%nlam))
				allocate(ObsSpec(i)%y(ObsSpec(i)%nlam))
				allocate(ObsSpec(i)%dy(ObsSpec(i)%nlam))
				ObsSpec(i)%y=0d0
				ObsSpec(i)%dy=1d0
				ObsSpec(i)%spec=.false.
			case('priors','prior')
				ObsSpec(i)%nlam=n_ret
				allocate(ObsSpec(i)%lam(ObsSpec(i)%nlam))
				allocate(ObsSpec(i)%y(ObsSpec(i)%nlam))
				allocate(ObsSpec(i)%dy(ObsSpec(i)%nlam))
				do j=1,n_ret
					if(RetPar(j)%logscale) then
						ObsSpec(i)%y(j)=log10(RetPar(j)%x0)
						ObsSpec(i)%dy(j)=log10(RetPar(j)%dx)
					else
						ObsSpec(i)%y(j)=RetPar(j)%x0
						ObsSpec(i)%dy(j)=RetPar(j)%dx/ObsSpec(i)%beta
					endif
				enddo
				ObsSpec(i)%spec=.false.
			case default
				open(unit=20,file=ObsSpec(i)%file,RECL=1000)
				j=1
				ilam=1
1				read(20,*,end=2,err=1) x,y,dy
				if(x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)) ilam=ilam+1
				j=j+1
				goto 1
2				ObsSpec(i)%nlam=ilam-1
				nj=j-1
				rewind(20)
				allocate(ObsSpec(i)%lam(ObsSpec(i)%nlam))
				allocate(ObsSpec(i)%y(ObsSpec(i)%nlam))
				allocate(ObsSpec(i)%dy(ObsSpec(i)%nlam))
				allocate(ObsSpec(i)%R(ObsSpec(i)%nlam))
				allocate(ObsSpec(i)%Rexp(ObsSpec(i)%nlam))
				allocate(ObsSpec(i)%model(ObsSpec(i)%nlam))
				ilam=1
				do j=1,nj
3					read(20,'(a1000)',err=3) line
					specres_obs=specres
					expspecres_obs=2d0
					read(line,*,err=3,end=4) x,y,dy,specres_obs,expspecres_obs
4					if(x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)) then
						ObsSpec(i)%lam(ilam)=x*1d-4
						ObsSpec(i)%y(ilam)=y
						ObsSpec(i)%dy(ilam)=dy/ObsSpec(i)%beta
						ObsSpec(i)%R(ilam)=specres_obs
						ObsSpec(i)%Rexp(ilam)=expspecres_obs
						ilam=ilam+1
					endif
				enddo
		end select
		if(ObsSpec(i)%beta.lt.0d0) ObsSpec(i)%beta=ObsSpec(i)%nlam
	enddo
	
	return
	end



	subroutine DoRetrieval()
	use GlobalSetup
	IMPLICIT NONE
	external ComputeChi2
	real*8 var0(n_ret),dvar0(2,n_ret),ComputeChi2,dvar(n_ret),x(n_ret),chi2min,random
	real*8,allocatable :: y(:),y1(:),y2(:),dy(:,:),yobs(:),dyobs(:),ybest(:),y0(:)
	real*8 chi2obs(nobs),var(n_ret),chi2,gasdev,maxd,error(2,n_ret),var_best(n_ret)
	real*8 x1,x2,minT(nr),maxT(nr),ran1,tot,lambda,chi2_1,chi2_2,dchi2(n_ret),chi2prev
	integer imodel,ny,i,j,iter1,iter2,iy,k
	
	real*8 dvar_av(n_ret),Cov(n_ret,n_ret),b(n_ret*3),W(n_ret*3,n_ret+1),Winv(n_ret,n_ret),dmax,scale
	real*8 chi2_spec,chi2_prof,maxsig,WLU(n_ret,n_ret),ErrVec(n_ret,n_ret)
	real*8 backup_beta(nobs)
	integer na,map(n_ret),info,iboot,nboot,niter1,niter2,n_not_improved
	logical dofit(n_ret),dofit_prev(n_ret),new_best
	logical,allocatable :: specornot(:)

	integer ME,MA,MG,MODE,MDW,ii
	integer,allocatable :: IP(:)
	real*8 PRGOPT(10),RNORME,RNORML
	real*8,allocatable :: WS(:)

	real*8 XeqCloud_best(nr,max(nclouds,1)),mixrat_best_r(nr,nmol)

	allocate(WS(2*(n_ret)+n_ret*3+(n_ret*2+2)*(n_ret+7)))
	allocate(IP(n_ret*4+2))
	IP(1)=2*(n_ret)+n_ret*3+(n_ret*2+2)*(n_ret+7)
	IP(2)=n_ret*4+2
	
	do i=1,n_ret
10		var0(i)=gasdev(idum)*0.1+0.5
		var0(i)=0.5d0
		if(var0(i).gt.1d0) goto 10
		if(var0(i).lt.0d0) goto 10
	enddo
c	call MapTprofile(var0)

	dvar0=10d0
	open(unit=31,file=trim(outputdir) // "Wolk.dat",RECL=6000)

	if(ngen.gt.0) then
c first genetic algoritm to make the first estimate
c		call Genetic(ComputeChi2,var0,dvar0,n_ret,nobs,npop,ngen,idum,gene_cross,.true.)
	endif

	imodel=npop*ngen
	ny=0
	do i=1,nobs
		ny=ny+ObsSpec(i)%nlam
	enddo
	allocate(y(ny))
	allocate(y0(ny))
	allocate(ybest(ny))
	allocate(y1(ny))
	allocate(y2(ny))
	allocate(dy(n_ret,ny))
	allocate(yobs(ny))
	allocate(dyobs(ny))
	allocate(specornot(ny))
	iy=1
	do i=1,nobs
		yobs(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%y(1:ObsSpec(i)%nlam)
		dyobs(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%dy(1:ObsSpec(i)%nlam)
		specornot(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%spec
		iy=iy+ObsSpec(i)%nlam		
	enddo
	dvar=0.1d0
	var=var0
	var_best=var0
	chi2min=1d200
	chi2prev=1d200
	error=0.1d0
	lambda=0.1
	Cov=0d0

	niter1=7
	niter2=15
	n_not_improved=0
	new_best=.false.

	if(ngen.gt.0) then
		do i=1,nobs
			backup_beta(i)=ObsSpec(i)%beta
			if(.not.ObsSpec(i)%spec) then
				ObsSpec(i)%beta=0d0
			endif
		enddo
		var0=var
		dvar0=10d0
		call Genetic(ComputeChi2,var0,dvar0,n_ret,nobs,npop,ngen,idum,gene_cross,.true.)
		chi2prev=1d200
		imodel=npop*ngen
		do i=1,nobs
			ObsSpec(i)%beta=backup_beta(i)
		enddo
		var=var0
	endif

	do iter1=1,niter1

	do iter2=1,niter2
		j=0
		imodel=imodel+1
		call output("Iteration " // trim(int2string(iter2+niter2*(iter1-1),'(i4)')) // 
     &				" - model" // trim(int2string(imodel,'(i5)')))
		chi2=1d0/ComputeChi2(imodel,n_ret,var,nobs,chi2obs,.true.,error)
		mixrat_old_r=mixrat_r
		XeqCloud_old=XeqCloud
		iy=1
		do i=1,nobs
	 		call RemapObs(i,y(iy:iy+ObsSpec(i)%nlam-1))
			iy=iy+ObsSpec(i)%nlam
		enddo
		if(chi2.le.chi2prev) then
			var0=var
			n_not_improved=0
			call output("Iteration improved" // trim(dbl2string(chi2,'(f8.3)')))
			lambda=lambda/5d0
			chi2_spec=0d0
			chi2_prof=0d0
			do j=1,nobs
				if(ObsSpec(j)%spec) then
					chi2_spec=chi2_spec+1d0/chi2obs(j)
				else
					chi2_prof=chi2_prof+1d0/chi2obs(j)
				endif
			enddo
			call output("   " // trim(dbl2string(chi2_spec,'(f8.3)')) // trim(dbl2string(chi2_prof,'(f8.3)')))
			call SetOutputMode(.false.)
			call WriteStructure()
			call WriteOutput()
			call WriteRetrieval(imodel,chi2_spec,var(1:n_ret))
			call WritePTlimits(var0,Cov,ErrVec,error,chi2*real(ny)/real(max(1,ny-n_ret)))
			call SetOutputMode(.true.)
			if(chi2.le.chi2min) then
				new_best=.false.
				var_best=var
				ybest=y
				mixrat_best_r=mixrat_r
				XeqCloud_best=XeqCloud
				chi2min=chi2
			else
				new_best=.false.
				var=var_best
				var0=var_best
				y=ybest
				mixrat_r=mixrat_best_r
				XeqCloud=XeqCloud_best
				chi2=chi2min
				n_not_improved=0
			endif
			chi2prev=chi2
		else if(iter2.eq.1) then
			dofit=.true.
			var=var0
		else if(n_not_improved.lt.10.or..not.new_best) then
			n_not_improved=n_not_improved+1
			lambda=lambda*5d0
			var=var0
			y=y0
			dofit=dofit_prev
			goto 1
		else
			new_best=.false.
			var=var_best
			var0=var_best
			y=ybest
			mixrat_r=mixrat_best_r
			XeqCloud=XeqCloud_best
			chi2=chi2prev
			n_not_improved=0
		endif
		do i=1,n_ret
			dvar(i)=max(error(1,i),error(2,i))
			if(dvar(i).eq.0d0) dvar(i)=0.1d0
			if(dvar(i).gt.0.1d0) dvar(i)=0.1d0
     	enddo
c		if(iter1*iter2.eq.1) dvar=1d0
c		chi2prev=chi2
		y0=y
		dofit=.true.
		dofit_prev=.true.
		do i=1,n_ret
			call output("varying " // trim(RetPar(i)%keyword))
			if(dvar(i).lt.1d-3) dvar(i)=1d-3
c			if(dvar(i).gt.1d-1) dvar(i)=1d-1
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
			if(chi2.lt.chi2min) then	!.and.RetPar(i)%opacitycomp) then
				new_best=.true.
				chi2min=chi2
				var_best=var
				ybest=y1
				mixrat_best_r=mixrat_r
				XeqCloud_best=XeqCloud
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
				call output("   " // trim(dbl2string(chi2_spec,'(f8.3)')) // trim(dbl2string(chi2_prof,'(f8.3)')))
c				call SetOutputMode(.false.)
c				call WriteStructure()
c				call WriteOutput()
c				call WriteRetrieval(imodel,chi2_spec,var(1:n_ret))
c				call WritePTlimits(var,Cov,ErrVec,error,chi2*real(ny)/real(max(1,ny-n_ret)))
c				call SetOutputMode(.true.)
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
			if(chi2.lt.chi2min) then	!.and.RetPar(i)%opacitycomp) then
				new_best=.true.
				chi2min=chi2
				var_best=var
				ybest=y2
				mixrat_best_r=mixrat_r
				XeqCloud_best=XeqCloud
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
				call output("   " // trim(dbl2string(chi2_spec,'(f8.3)')) // trim(dbl2string(chi2_prof,'(f8.3)')))
c				call SetOutputMode(.false.)
c				call WriteStructure()
c				call WriteOutput()
c				call WriteRetrieval(imodel,chi2_spec,var(1:n_ret))
c				call WritePTlimits(var,Cov,ErrVec,error,chi2*real(ny)/real(max(1,ny-n_ret)))
c				call SetOutputMode(.true.)
			endif
			dy(i,1:ny)=(y1(1:ny)-y2(1:ny))/(x1-x2)
			dchi2(i)=(chi2_1-chi2_2)/(x1-x2)
			if(abs(chi2_1-chi2_2)/(chi2_1+chi2_2).lt.1d-4) then
				if(dvar(i).lt.0.99d0) then
					dvar(i)=dvar(i)*5d0
					goto 20
				else if(dvar(i).lt.0.999d0) then
					dvar(i)=1d0
					goto 20
				endif
c				dofit(i)=.false.
c				var(i)=var0(i)
c				dvar(i)=1d0
c				call output("derivative too small")
c				call output("removing " // trim(RetPar(i)%keyword))
			endif
		enddo

		dofit_prev=dofit
1		continue

		do ii=1,n_ret
			W=0d0
			b=0d0
			do k=1,ny
				do i=1,n_ret
					do j=1,n_ret
						W(i,j)=W(i,j)+dy(i,k)*dy(j,k)/dyobs(k)**2
					enddo
				enddo
			enddo

			ME=0
			MA=n_ret
			MG=0
			MDW=n_ret*3

			PRGOPT(1)=1

			IP(1)=2*(n_ret)+n_ret*3+(n_ret*2+2)*(n_ret+7)
			IP(2)=n_ret*4+2

			W(1:n_ret,n_ret+1)=0d0
			W(ii,n_ret+1)=1d0
			call dlsei(w, MDW, ME, MA, MG, n_ret, PRGOPT, b, RNORME,
     +   RNORML, MODE, WS, IP)
   
			do j=1,n_ret
				Cov(ii,j)=b(j)
			enddo
		enddo
		Winv=Cov

		W=0d0
		b=0d0
		do k=1,ny
			do i=1,n_ret
				do j=1,n_ret
					W(i,j)=W(i,j)+dy(i,k)*dy(j,k)/dyobs(k)**2
				enddo
				b(i)=b(i)+(yobs(k)-y(k))*dy(i,k)/dyobs(k)**2
			enddo
		enddo
		do i=1,n_ret
			W(i,i)=W(i,i)*(1d0+lambda)
		enddo
		do j=1,n_ret
			W(j+n_ret,1:n_ret)=0d0
			b(j+n_ret)=-var0(j)
			W(j+n_ret,j)=1d0

			W(j+n_ret*2,1:n_ret)=0d0
			b(j+n_ret*2)=var0(j)-1d0
			W(j+n_ret*2,j)=-1d0
		enddo

		ME=0
		MA=n_ret
		MG=n_ret*2
		MDW=n_ret*3

		PRGOPT(1)=1

		IP(1)=2*(n_ret)+n_ret*3+(n_ret*2+2)*(n_ret+7)
		IP(2)=n_ret*4+2

		W(1:MDW,n_ret+1)=b(1:MDW)
		call dlsei(w, MDW, ME, MA, MG, n_ret, PRGOPT, b, RNORME,
     +   RNORML, MODE, WS, IP)

		dvar=b
		var=var0+dvar
		call WritePTlimits(var,Cov,ErrVec,error,chi2min)
		do i=1,n_ret
			if(var(i).gt.1d0) var(i)=1d0
			if(var(i).lt.0d0) var(i)=0d0
			dvar(i)=sqrt(error(1,i)**2+error(2,i)**2)*3d0
		enddo
	enddo
	do i=1,nobs
		if(.not.ObsSpec(i)%spec.and.ObsSpec(i)%scale.gt.0d0) then
			print*,chi2_spec,chi2_prof,chi2min
			ObsSpec(i)%beta=log10(ObsSpec(i)%scale*ObsSpec(i)%beta*sqrt(chi2_spec/chi2_prof))+log10(ObsSpec(i)%beta)
			ObsSpec(i)%beta=10d0**(ObsSpec(i)%beta/2d0)
			call output("Adjusting beta to " // trim(dbl2string(ObsSpec(i)%beta,'(es10.4)')))
			chi2min=1d200
		endif
	enddo
	enddo

	close(unit=31)
	
	return
	end



	subroutine WritePTlimits(var0,Cov,ErrVec,error,chi2)
	use GlobalSetup
	IMPLICIT NONE
	real*8 minT(nr),maxT(nr),var(n_ret),error(2,n_ret),tot,w(n_ret),chi2,ErrVec(n_ret,n_ret)
	real*8 var0(n_ret),random,vec(n_ret),Tbest(nr),gasdev,Cov(n_ret,n_ret)
	real*8 dvar(n_ret)
	integer i,j,k,info,nk,ierr
	character*6000 form

	call Eigenvalues(Cov,ErrVec,w,n_ret,INFO)

	call SetOutputMode(.false.)
	var=var0
	call MapRetrieval(var,error)
	call InitDens()
	call SetupStructure(.false.)
	Tbest=T
	maxT=0d0
	minT=1d200
	do i=1,n_ret
		error(1:2,i)=var0(i)
	enddo

	nk=n_ret*100
	open(unit=35,file=trim(outputdir) // "error_cloud.dat",RECL=6000)
	form='("#"' // trim(int2string(n_ret,'(i4)')) // 'a19)'
	write(35,form) (trim(int2string(i,'(i4)')) // " " // RetPar(i)%keyword,i=1,n_ret)
	form='(' // int2string(n_ret,'(i4)') // 'es19.7)'
	
	do k=1,nk
		vec=0d0
		tot=0d0
		do i=1,n_ret
			vec(i)=2d0*(random(idum)-0.5d0)
			tot=tot+vec(i)**2
		enddo
		vec=vec*(random(idum)**(1d0/real(n_ret)))
		vec=vec/sqrt(tot)
		vec=vec*sqrt(chi2)
		vec=vec*sqrt(w)
		do i=1,n_ret
			dvar(i)=0d0
			do j=1,n_ret
				dvar(i)=dvar(i)+vec(j)*ErrVec(i,j)
			enddo
			var(i)=var0(i)+dvar(i)
			if(var(i).gt.1d0) var(i)=1d0
			if(var(i).lt.0d0) var(i)=0d0
		enddo	

		call MapRetrieval(var,error)
		write(35,form) (RetPar(i)%value,i=1,n_ret)
		call InitDens()
		call SetupStructure(.false.)
		do i=1,nr
			if(T(i).gt.maxT(i)) maxT(i)=T(i)
			if(T(i).lt.minT(i)) minT(i)=T(i)
		enddo
		do i=1,n_ret
			if(var(i).lt.error(1,i)) error(1,i)=var(i)
			if(var(i).gt.error(2,i)) error(2,i)=var(i)
		enddo
	enddo
	close(unit=35)
	error(1,1:n_ret)=abs(error(1,1:n_ret)-var0(1:n_ret))
	error(2,1:n_ret)=abs(error(2,1:n_ret)-var0(1:n_ret))

	call MapRetrieval(var0,error)

	call SetOutputMode(.true.)
	open(unit=45,file=trim(outputdir) // "limits.dat")
	do i=1,nr
		write(45,*) P(i),minT(i),maxT(i)
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
c			call regridarray(lamobs,obsA(0,1:nlam-1)/(pi*Rstar**2),nlam-1,
c     &					ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
			call regridspecres(lamobs,obsA(0,1:nlam-1)/(pi*Rstar**2),nlam-1,
     &					ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%nlam)
     		ObsSpec(i)%model(1:ObsSpec(i)%nlam)=spec(1:ObsSpec(i)%nlam)
		case("emisr","emisR")
c			call regridarray(lamobs,flux(0,1:nlam-1)/(Fstar(1:nlam-1)*1d23/distance**2),
c     &					nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
c			call regridspecres(lamobs,flux(0,1:nlam-1)/(Fstar(1:nlam-1)*1d23/distance**2),
c     &					nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%nlam)
			call regridspecres(lamobs,(phase(1,0,1:nlam-1)+flux(0,1:nlam-1))/(Fstar(1:nlam-1)*1d23/distance**2),
     &					nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%nlam)
     		ObsSpec(i)%model(1:ObsSpec(i)%nlam)=spec(1:ObsSpec(i)%nlam)
		case("emisa","emis","emission")
c			call regridarray(lamobs,flux(0,1:nlam-1),
c     &					nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
			call regridspecres(lamobs,phase(1,0,1:nlam-1)+flux(0,1:nlam-1),
     &					nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%nlam)
     		ObsSpec(i)%model(1:ObsSpec(i)%nlam)=spec(1:ObsSpec(i)%nlam)
		case("tprofile")
			call ComputeParamT(spec(1:nr))
			spec(1:ObsSpec(i)%nlam)=spec(1:ObsSpec(i)%nlam)*ObsSpec(i)%beta
		case("logtp")
			spec(1)=(log10(abs(T(1)/T(2))/log10(P(1)/P(2))))
			do j=2,nr-1
				spec(j)=(log10(abs(T(j-1)/T(j)))/log10(P(j-1)/P(j)))-(log10(abs(T(j)/T(j+1)))/log10(P(j)/P(j+1)))
c				spec(j)=spec(j)+log10(abs(T(j-1)/T(j+1)))/log10(P(j-1)/P(j+1))
c	print*,j,T(j-1),T(j),T(j+1)
			enddo
			spec(nr)=(log10(abs(T(nr-1)/T(nr)))/log10(P(nr-1)/P(nr)))
			spec(1:ObsSpec(i)%nlam)=spec(1:ObsSpec(i)%nlam)/real(nr)
			spec(1:ObsSpec(i)%nlam)=spec(1:ObsSpec(i)%nlam)*ObsSpec(i)%beta
		case("prior","priors")
			do j=1,n_ret
				if(RetPar(j)%logscale) then
					spec(j)=log10(RetPar(j)%value)
				else
					spec(j)=RetPar(j)%value
				endif
			enddo
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

c	print*,imodel,(trim(dbl2string(RetPar(i)%value,'(es18.7)')),i=1,n_ret)

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
	real*8 var(n_ret),x,Tmap(nr)

	call InitDens()
	call SetupStructure(.false.)

	call ComputeParamT(Tmap)
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
	

	subroutine regridspecres(x0,y0,n0,x1,y1,R1,expR1,n1)
	IMPLICIT NONE
	integer i1,i0,n0,n1
	real*8 x0(n0),y0(n0),x1(n1),y1(n1),R1(n1),expR1(n1),w,tot

	do i1=1,n1
		if(x1(i1).lt.x0(1)) then
			y1(i1)=y0(1)
		else if(x1(i1).gt.x0(n0)) then
			y1(i1)=y0(n0)
		else
			y1(i1)=0d0
			tot=0d0
			do i0=1,n0
				w=exp(-abs((x0(i0)-x1(i1))*R1(i1)*2d0/x1(i1))**expR1(i1))
				if(i0.eq.1) then
					w=w*abs(x0(2)-x0(1))/2d0
				else if(i0.eq.n0) then
					w=w*abs(x0(n0)-x0(n0-1))/2d0
				else
					w=w*abs(x0(i0-1)-x0(i0+1))/2d0
				endif
				tot=tot+w
				y1(i1)=y1(i1)+w*y0(i0)
			enddo
			y1(i1)=y1(i1)/tot
		endif
	enddo

	return
	end
	

	subroutine WriteRetrieval(imodel,chi2,var)
	use GlobalSetup
	IMPLICIT NONE
	real*8 var(n_ret),chi2
	integer i,imodel,j

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
	
	do i=1,nobs
		select case(ObsSpec(i)%type)
			case("trans","transmission","emisr","emisR","emisa","emis","emission")
				open(unit=20,file=trim(outputdir) // "obs" // trim(int2string(i,'(i0.3)')),RECL=1000)
				do j=1,ObsSpec(i)%nlam
					write(20,*) ObsSpec(i)%lam(j)*1d4,ObsSpec(i)%model(j),ObsSpec(i)%y(j),ObsSpec(i)%dy(j)
				enddo
				close(unit=20)
		end select
	enddo

	call system("cp " // trim(outputdir) // "input.dat " // trim(outputdir) // "bestfit.dat")
	open(unit=21,file=trim(outputdir) // "bestfit.dat",RECL=1000,access='APPEND')
	write(21,'("*** retrieval keywords ***")')
	write(21,'("retrieval=.false.")')
	do i=1,n_ret
		write(21,'(a," = ",es14.7)') trim(RetPar(i)%keyword),RetPar(i)%value
	enddo
	close(unit=21)	

	return
	end


	subroutine MatrixInvert(A,Ainv,ALU,N,INFO)
	IMPLICIT NONE
	integer n,ierr
	real*8 A(n,n),Ainv(n,n),ALU(n,n)
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
      ALU=Ainv

C  DGETRI computes the inverse of a matrix using the LU factorization
C  computed by DGETRF.

      CALL DGETRI(N, Ainv, N, IPIV, WORK, LWORK, INFO)

      DEALLOCATE (IPIV, STAT = DeAllocateStatus)
      DEALLOCATE (WORK, STAT = DeAllocateStatus)

	return
	end
	
	
	
	subroutine Eigenvalues(A,E,v,N,INFO)
	IMPLICIT NONE
	integer N,INFO
	real*8 A(N,N),E(N,N),v(N),w(1000)
	real*8, ALLOCATABLE :: WORK(:)
	integer LWORK

	LWORK=-1
	call DSYEV ('V', 'U', N, E, N, v, w, LWORK, INFO)
	LWORK=max(1d0,w(1))
	allocate(WORK(LWORK))
	E=A
	call DSYEV ('V', 'U', N, E, N, v, WORK, LWORK, INFO)
	deallocate(WORK)
	
	return
	end

	
	
	