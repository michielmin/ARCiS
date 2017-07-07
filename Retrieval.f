	subroutine ReadObs()
	use GlobalSetup
	IMPLICIT NONE
	integer i,j,ilam,nj
	real*8 x,y,dy,specres_obs,expspecres_obs
	character*1000 line
	real*8 scale,scale_av,d,dmin
	integer nscale,i2,j2
	
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
						ObsSpec(i)%dy(j)=log10(RetPar(j)%dx)/ObsSpec(i)%beta
					else
						ObsSpec(i)%y(j)=RetPar(j)%x0
						ObsSpec(i)%dy(j)=RetPar(j)%dx/ObsSpec(i)%beta
					endif
				enddo
				ObsSpec(i)%spec=.false.
			case default
				ObsSpec(i)%spec=.true.
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
				allocate(ObsSpec(i)%model0(ObsSpec(i)%nlam))
				allocate(ObsSpec(i)%modelbest(ObsSpec(i)%nlam))
				ilam=1
				if(ObsSpec(i)%beta.lt.0d0) ObsSpec(i)%beta=ObsSpec(i)%nlam
				do j=1,nj
3					read(20,'(a1000)',err=3) line
					specres_obs=specres
					expspecres_obs=20d0
					read(line,*,err=3,end=4) x,y,dy,specres_obs,expspecres_obs
4					continue
					if(dy.lt.1d-2*y) dy=1d-2*y
					if(x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)) then
						ObsSpec(i)%lam(ilam)=x*1d-4
						if(ObsSpec(i)%type.eq."emisa".or.ObsSpec(i)%type.eq."emis".or.ObsSpec(i)%type.eq."emission") then
							dy=dy/y
							ObsSpec(i)%y(ilam)=log(y)
						else
							ObsSpec(i)%y(ilam)=y
						endif
						ObsSpec(i)%dy(ilam)=dy/ObsSpec(i)%beta
						ObsSpec(i)%R(ilam)=specres_obs
						ObsSpec(i)%Rexp(ilam)=expspecres_obs
						ilam=ilam+1
					endif
				enddo
				close(unit=20)
		end select
	enddo

	if(faircoverage) then
		scale_av=0d0
		nscale=0
		do i=1,nobs
			if(ObsSpec(i)%spec) then
				do j=1,ObsSpec(i)%nlam
					scale=0d0
					do i2=1,nobs
						if(ObsSpec(i2)%spec) then
							do j2=1,ObsSpec(i2)%nlam
								d=ObsSpec(i)%lam(j)/ObsSpec(i2)%lam(j2)
								d=log10(d)/log10(2d0)
								scale=scale+exp(-d**2)
							enddo
						endif
					enddo
					ObsSpec(i)%dy(j)=ObsSpec(i)%dy(j)*scale
					scale_av=scale_av+scale
					nscale=nscale+1
				enddo
			endif
		enddo
		scale_av=scale_av/real(nscale)
		do i=1,nobs
			if(ObsSpec(i)%spec) then
				do j=1,ObsSpec(i)%nlam
					ObsSpec(i)%dy(j)=ObsSpec(i)%dy(j)/scale_av
				enddo
			endif
		enddo
	endif
	
	return
	end



	subroutine DoRetrieval()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	external ComputeChi2
	real*8 var0(n_ret),dvar0(2,n_ret),ComputeChi2,dvar(n_ret),x(n_ret),chi2min,random
	real*8,allocatable :: y(:),y1(:),y2(:),dy(:,:),yobs(:),dyobs(:),ybest(:),y0(:)
	real*8 chi2obs(nobs),var(n_ret),chi2,gasdev,maxd,error(2,n_ret),var_best(n_ret),chi2_0
	real*8 x1,x2,minT(nr),maxT(nr),ran1,tot,lambda,chi2_1,chi2_2,dchi2(n_ret),chi2prev
	integer imodel,ny,i,j,iter1,iter2,iy,k
	
	real*8 dvar_av(n_ret),Cov(n_ret,n_ret),b(n_ret*3),W(n_ret*3,n_ret+1),Winv(n_ret,n_ret),dmax
	real*8 chi2_spec,chi2_prof,maxsig,WLU(n_ret,n_ret),ErrVec(n_ret,n_ret),dvar_prev(n_ret)
	real*8 backup_xmin(n_ret),backup_xmax(n_ret)
	real*8 obsA1(nlam),obsA2(nlam),dobsA(n_ret,nlam)
	real*8 emis1(nlam),emis2(nlam),demis(n_ret,nlam)
	real*8 emisR1(nlam),emisR2(nlam),demisR(n_ret,nlam)
	real*8 phase0(nlam),flux0(nlam),obsA0(nlam),scale,scalemin
	integer na,map(n_ret),info,iboot,nboot,niter1,niter2,n_not_improved
	logical dofit(n_ret),dofit_prev(n_ret),new_best,improved_iter
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
		RetPar(i)%value=RetPar(i)%x0
	enddo
	call MapRetrievalInverse(var0)
	do i=1,n_ret
		if(var0(i).gt.1d0) var0(i)=1d0
		if(var0(i).lt.0d0) var0(i)=0d0
	enddo

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
		dyobs(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%dy(1:ObsSpec(i)%nlam)*sqrt(real(ObsSpec(i)%nlam))
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
		do i=1,n_ret
			if(RetPar(i)%keyword(1:6).eq."tvalue") then
				backup_xmin(i)=RetPar(i)%xmin
				backup_xmax(i)=RetPar(i)%xmax
				RetPar(i)%xmin=-1d0
				RetPar(i)%xmax=1d0
			endif
		enddo
		var0=var
		dvar0=10d0
		call Genetic(ComputeChi2,var0,dvar0,n_ret,nobs,npop,ngen,idum,gene_cross,.true.)
		chi2prev=1d200
		imodel=npop*ngen
		do i=1,n_ret
			if(RetPar(i)%keyword(1:6).eq."tvalue") then
				RetPar(i)%xmin=backup_xmin(i)
				RetPar(i)%xmax=backup_xmax(i)
				var0(i)=-RetPar(i)%xmin/(RetPar(i)%xmax-RetPar(i)%xmin)
				if(var0(i).lt.0d0) var0(i)=0d0
			endif
		enddo
		var=var0
	endif


	do i=1,n_ret
		ErrVec(:,i)=0d0
		ErrVec(i,i)=dvar(i)
	enddo

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
		improved_iter=.false.

		if(chi2.le.chi2prev) then
			call output("Iteration improved" // trim(dbl2string(chi2,'(f13.3)')))
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
			call output("   " // trim(dbl2string(chi2_spec,'(f13.3)')) // trim(dbl2string(chi2_prof,'(f13.3)')))
		else
			lambda=lambda*5d0
		endif
		chi2_spec=0d0
		chi2_prof=0d0
		do j=1,nobs
			if(ObsSpec(j)%spec) then
				chi2_spec=chi2_spec+1d0/chi2obs(j)
			else
				chi2_prof=chi2_prof+1d0/chi2obs(j)
			endif
		enddo
		call output("   " // trim(dbl2string(chi2_spec,'(f13.3)')) // trim(dbl2string(chi2_prof,'(f13.3)')))

		if(chi2.le.chi2min) then
			new_best=.false.
			var_best=var
			var0=var
			ybest=y
			mixrat_best_r=mixrat_r
			XeqCloud_best=XeqCloud
			chi2min=chi2
			do i=1,nobs
				if(ObsSpec(i)%spec) then
					ObsSpec(i)%modelbest=ObsSpec(i)%model
				endif
			enddo
		else
			new_best=.false.
			var=var_best
			var0=var_best
			y=ybest
			mixrat_r=mixrat_best_r
			XeqCloud=XeqCloud_best
			chi2=chi2min
			n_not_improved=0
			do i=1,nobs
				if(ObsSpec(i)%spec) then
					ObsSpec(i)%model=ObsSpec(i)%modelbest
				endif
			enddo
		endif
		chi2prev=chi2

		improved_iter=.true.

C	if(chi2.le.chi2prev) then
C		var0=var
C		n_not_improved=0
C		improved_iter=.true.
C		call output("Iteration improved" // trim(dbl2string(chi2,'(f13.3)')))
C		lambda=lambda/5d0
C		chi2_spec=0d0
C		chi2_prof=0d0
C		do j=1,nobs
C			if(ObsSpec(j)%spec) then
C				chi2_spec=chi2_spec+1d0/chi2obs(j)
C			else
C				chi2_prof=chi2_prof+1d0/chi2obs(j)
C			endif
C		enddo
C		call output("   " // trim(dbl2string(chi2_spec,'(f13.3)')) // trim(dbl2string(chi2_prof,'(f13.3)')))
C		if(chi2.le.chi2min) then
C			new_best=.false.
C			var_best=var
C			ybest=y
C			mixrat_best_r=mixrat_r
C			XeqCloud_best=XeqCloud
C			chi2min=chi2
C			do i=1,nobs
C				if(ObsSpec(i)%spec) then
C					ObsSpec(i)%modelbest=ObsSpec(i)%model
C				endif
C			enddo
C		else
C			new_best=.false.
C			var=var_best
C			var0=var_best
C			y=ybest
C			mixrat_r=mixrat_best_r
C			XeqCloud=XeqCloud_best
C			chi2=chi2min
C			n_not_improved=0
C			do i=1,nobs
C				if(ObsSpec(i)%spec) then
C					ObsSpec(i)%model=ObsSpec(i)%modelbest
C				endif
C			enddo
C		endif
C		chi2prev=chi2
C	else if(iter2.eq.1) then
C		dofit=.true.
C		var=var0
C	else if(n_not_improved.lt.10.or..not.new_best) then
C		n_not_improved=n_not_improved+1
C		lambda=lambda*5d0
C		improved_iter=.false.
C		new_best=.false.
C		var=var_best
C		var0=var_best
C		y=ybest
C		mixrat_r=mixrat_best_r
C		XeqCloud=XeqCloud_best
C		chi2=chi2min
C		n_not_improved=0
C		do i=1,nobs
C			if(ObsSpec(i)%spec) then
C				ObsSpec(i)%model=ObsSpec(i)%modelbest
C			endif
C		enddo
C	else
C		improved_iter=.false.
C		lambda=0.01d0
C		new_best=.false.
C		var=var_best
C		var0=var_best
C		y=ybest
C		mixrat_r=mixrat_best_r
C		XeqCloud=XeqCloud_best
C		chi2=chi2prev
C		n_not_improved=0
C		do i=1,nobs
C			if(ObsSpec(i)%spec) then
C				ObsSpec(i)%model=ObsSpec(i)%modelbest
C			endif
C		enddo
C	endif

		do i=1,n_ret
			dvar(i)=max(error(1,i),error(2,i))
			if(dvar(i).eq.0d0) dvar(i)=0.1d0
			dvar(i)=(dvar(i)+dvar_prev(i))/2d0
     	enddo
		var=var0
		chi2=1d0/ComputeChi2(imodel,n_ret,var,nobs,chi2obs,.true.,error)
		mixrat_old_r=mixrat_r
		XeqCloud_old=XeqCloud
		iy=1
		do i=1,nobs
	 		call RemapObs(i,y(iy:iy+ObsSpec(i)%nlam-1))
			iy=iy+ObsSpec(i)%nlam
		enddo
		dvar_prev=dvar
		y0=y
		obsA0(1:nlam)=obsA(0,1:nlam)
		phase0(1:nlam)=phase(1,0,1:nlam)
		flux0(1:nlam)=flux(0,1:nlam)
		chi2_0=chi2
		do j=1,nobs
			if(ObsSpec(j)%spec) then
				ObsSpec(j)%model0=ObsSpec(j)%model
			endif
		enddo
		dofit=.true.
		dofit_prev=.true.
		do i=1,n_ret
			call output("varying " // trim(RetPar(i)%keyword))
			if(dvar(i).lt.1d-3) dvar(i)=1d-3
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
			chi2_1=0d0
			do j=1,nobs
				if(ObsSpec(j)%spec) then
					chi2_1=chi2_1+1d0/chi2obs(j)
				endif
			enddo
			obsA1(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
			emis1(1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
			emisR1(1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar*1d23/distance**2)
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
						ObsSpec(j)%modelbest=ObsSpec(j)%model
					else
						chi2_prof=chi2_prof+1d0/chi2obs(j)
					endif
				enddo
				call output("Updating best fit " // trim(dbl2string(chi2,'(f13.3)')))
				call output("   " // trim(dbl2string(chi2_spec,'(f13.3)')) // trim(dbl2string(chi2_prof,'(f13.3)')))
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
			chi2_2=0d0
			do j=1,nobs
				if(ObsSpec(j)%spec) then
					chi2_2=chi2_2+1d0/chi2obs(j)
				endif
			enddo
			obsA2(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
			emis2(1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
			emisR2(1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar*1d23/distance**2)
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
						ObsSpec(j)%modelbest=ObsSpec(j)%model
					else
						chi2_prof=chi2_prof+1d0/chi2obs(j)
					endif
				enddo
				call output("Updating best fit " // trim(dbl2string(chi2,'(f13.3)')))
				call output("   " // trim(dbl2string(chi2_spec,'(f13.3)')) // trim(dbl2string(chi2_prof,'(f13.3)')))
			endif
			dy(i,1:ny)=(y1(1:ny)-y2(1:ny))/(x1-x2)
			dchi2(i)=(chi2_1-chi2_2)/(x1-x2)
			dobsA(i,1:nlam)=(obsA1(1:nlam)-obsA2(1:nlam))/(x1-x2)
			demis(i,1:nlam)=(emis1(1:nlam)-emis2(1:nlam))/(x1-x2)
			demisR(i,1:nlam)=(emisR1(1:nlam)-emisR2(1:nlam))/(x1-x2)
		enddo

		dofit_prev=dofit
1		continue

c================================================
		W=0d0
		b=0d0
		do k=1,ny
			do i=1,n_ret
				do j=1,n_ret
					W(i,j)=W(i,j)+dy(i,k)*dy(j,k)/dyobs(k)**2
				enddo
			enddo
		enddo
		call MatrixInvert(W(1:n_ret,1:n_ret),Winv(1:n_ret,1:n_ret),WLU(1:n_ret,1:n_ret),n_ret,info)
		if(INFO.eq.0) then
			do i=1,n_ret
				do j=1,n_ret
					Cov(i,j)=Winv(i,j)
				enddo
			enddo
c================================================
		else
			call output('Error in matrix inversion.')
			call output('Retrieving parameters with no influence on the observations?')
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
		endif

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

		PRGOPT(1)=1

		IP(1)=2*(n_ret)+n_ret*3+(n_ret*2+2)*(n_ret+7)
		IP(2)=n_ret*4+2

		W(1:MDW,n_ret+1)=b(1:MDW)
		call dlsei(w, MDW, ME, MA, MG, n_ret, PRGOPT, b, RNORME,
     +   RNORML, MODE, WS, IP)

		dvar(1:n_ret)=b(1:n_ret)
		scalemin=1d0
		var=var0+dvar
		do i=1,n_ret
			if(var(i).gt.1d0) then
				var(i)=1d0
			endif
			if(var(i).lt.0d0) then
				var(i)=0d0
			endif
		enddo
		var=var0+0.75d0*scalemin*dvar
		do i=1,n_ret
			if(var(i).gt.1d0) var(i)=1d0
			if(var(i).lt.0d0) var(i)=0d0
		enddo
		if(improved_iter) then
			y=y0
			obsA(0,1:nlam)=obsA0(1:nlam)
			phase(1,0,1:nlam)=phase0(1:nlam)
			flux(0,1:nlam)=flux0(1:nlam)
			dofit=dofit_prev
			do j=1,nobs
				if(ObsSpec(j)%spec) then
					ObsSpec(j)%model=ObsSpec(j)%model0
				endif
			enddo
			iy=1
			do j=1,nobs
		 		call RemapObs(j,y(iy:iy+ObsSpec(j)%nlam-1))
				iy=iy+ObsSpec(j)%nlam
			enddo
			call SetOutputMode(.false.)
			call WritePTlimits(var0,Cov,ErrVec,error,chi2_0,dobsA,demis,demisR,.true.)
			call WriteRetrieval(imodel,chi2_0,var0,error)
			call WriteStructure()
			call WriteOutput()
			call SetOutputMode(.true.)
		endif
		call WritePTlimits(var,Cov,ErrVec,error,chi2_0,dobsA,demis,demisR,.false.)
		do i=1,n_ret
			if(var(i).gt.1d0) var(i)=1d0
			if(var(i).lt.0d0) var(i)=0d0
			dvar(i)=sqrt(error(1,i)**2+error(2,i)**2)*max(lambda*10d0,0.1d0)
		enddo
	enddo
	do i=1,nobs
		if(.not.ObsSpec(i)%spec.and.ObsSpec(i)%scale.gt.0d0) then
			print*,chi2_spec,chi2_prof,chi2min
			ObsSpec(i)%beta=log10(ObsSpec(i)%scale*ObsSpec(i)%beta*sqrt(chi2_spec/chi2_prof))+log10(ObsSpec(i)%beta)
			ObsSpec(i)%beta=10d0**(ObsSpec(i)%beta/2d0)
			call output("Adjusting beta to " // trim(dbl2string(ObsSpec(i)%beta,'(es10.4)')))
			chi2min=1d200
			chi2prev=1d200
		endif
	enddo
	enddo

	close(unit=31)
	
	return
	end



	subroutine WritePTlimits(var0,Cov,ErrVec,error,chi2,dobsA,demis,demisR,ioflag)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 minT(nr),maxT(nr),var(n_ret),error(2,n_ret),tot,w(n_ret),chi2,ErrVec(n_ret,n_ret)
	real*8 var0(n_ret),random,vec(n_ret),Tbest(nr),gasdev,Cov(n_ret,n_ret)
	real*8 dvar(n_ret),dobsA(n_ret,nlam),demis(n_ret,nlam),demisR(n_ret,nlam)
	real*8 obsAmin(nlam),obsAmax(nlam),emismin(nlam),emismax(nlam),emisRmin(nlam),emisRmax(nlam)
	real*8 emis0(nlam),obsA0(nlam),emisR0(nlam)
	integer i,j,k,info,nk,ierr,iter
	character*6000 form
	logical ioflag

	call Eigenvalues(Cov,ErrVec,w,n_ret,INFO)
	do i=1,n_ret
		tot=0d0
		do j=1,n_ret
			tot=tot+ErrVec(j,i)**2
		enddo
		ErrVec(1:n_ret,i)=ErrVec(1:n_ret,i)*sqrt(w(i))*chi2/sqrt(tot)
	enddo
	
	call SetOutputMode(.false.)
	var=var0
	call MapRetrieval(var,error)
	call InitDens()
	call SetupStructure(.false.)
	Tbest=T
	maxT=0d0
	minT=1d200
	obsAmin=1d200
	emismin=1d200
	emisRmin=1d200
	obsAmax=0d0
	emismax=0d0
	emisRmax=0d0
	COerr(1)=1d200
	COerr(2)=0d0

	do i=1,n_ret
		error(1:2,i)=var0(i)
	enddo

	nk=n_ret*100
	if(ioflag) then
		open(unit=35,file=trim(outputdir) // "error_cloud.dat",RECL=6000)
		form='("#"' // trim(int2string(n_ret,'(i4)')) // 'a19)'
		write(35,form) (trim(int2string(i,'(i4)')) // " " // RetPar(i)%keyword,i=1,n_ret)
		form='(' // int2string(n_ret,'(i4)') // 'es19.7)'
	endif
	
	do k=1,nk
		iter=0
1		continue
c		if(k.le.n_ret) then
c			vec=0d0
c			vec(k)=sqrt(chi2)
c			iter=1000
c		else if(k.le.n_ret*2) then
c			vec=0d0
c			vec(k-n_ret)=-sqrt(chi2)
c			iter=1000
c		else
			vec=0d0
			tot=0d0
			do i=1,n_ret
				vec(i)=2d0*(random(idum)-0.5d0)
				tot=tot+vec(i)**2
			enddo
			vec=vec*(random(idum)**(1d0/real(n_ret)))
			vec=vec/sqrt(tot)
c		endif
		do i=1,n_ret
			dvar(i)=0d0
			do j=1,n_ret
				dvar(i)=dvar(i)+vec(j)*ErrVec(i,j)
			enddo
		enddo
2		continue
		iter=iter+1
		if(speclimits) then
			obsA0=0d0
			emis0=0d0
			emisR0=0d0
		endif
		do i=1,n_ret
			var(i)=var0(i)+dvar(i)
			if(var(i).gt.1d0) then
				if(iter.lt.100) goto 1
				dvar=dvar*0.999d0*(1d0-var0(i))/(var(i)-var0(i))
				goto 2
			endif
			if(var(i).lt.0d0) then
				if(iter.lt.100) goto 1
				dvar=dvar*0.999d0*(-var0(i))/(var(i)-var0(i))
				goto 2
			endif
			if(speclimits) then
				obsA0(1:nlam)=obsA0(1:nlam)+dobsA(i,1:nlam)*dvar(i)
				emis0(1:nlam)=emis0(1:nlam)+demis(i,1:nlam)*dvar(i)
				emisR0(1:nlam)=emisR0(1:nlam)+demisR(i,1:nlam)*dvar(i)
			endif
		enddo
		if(speclimits) then
			obsA0=obsA0+obsA(0,1:nlam)/(pi*Rstar**2)
			emis0=emis0+phase(1,0,1:nlam)+flux(0,1:nlam)
			emisR0=emisR0+(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar*1d23/distance**2)
			do i=1,nlam
				if(obsA0(i).lt.0d0) obsA0(i)=0d0
				if(emis0(i).lt.0d0) emis0(i)=0d0
				if(emisR0(i).lt.0d0) emisR0(i)=0d0
			enddo
		endif

		call MapRetrieval(var,error)
		
		if(ioflag) write(35,form) (RetPar(i)%value,i=1,n_ret)
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
		if(COret.lt.COerr(1)) COerr(1)=COret
		if(COret.gt.COerr(2)) COerr(2)=COret
		if(speclimits) then
			do i=1,nlam
				if(obsA0(i).gt.obsAmax(i)) obsAmax(i)=obsA0(i)
				if(emis0(i).gt.emismax(i)) emismax(i)=emis0(i)
				if(emisR0(i).gt.emisRmax(i)) emisRmax(i)=emisR0(i)
				if(obsA0(i).lt.obsAmin(i)) obsAmin(i)=obsA0(i)
				if(emis0(i).lt.emismin(i)) emismin(i)=emis0(i)
				if(emisR0(i).lt.emisRmin(i)) emisRmin(i)=emisR0(i)
			enddo
		endif
	enddo
	error(1,1:n_ret)=abs(error(1,1:n_ret)-var0(1:n_ret))
	error(2,1:n_ret)=abs(error(2,1:n_ret)-var0(1:n_ret))

	call MapRetrieval(var0,error)
	call InitDens()
	call SetupStructure(.false.)

	COerr(1)=abs(COerr(1)-COret)
	COerr(2)=abs(COerr(2)-COret)

	call SetOutputMode(.true.)
	if(ioflag) then
		close(unit=35)

		open(unit=45,file=trim(outputdir) // "limits.dat")
		do i=1,nr
			write(45,*) P(i),minT(i),maxT(i)
		enddo
		close(unit=45)

		if(speclimits) then
			open(unit=45,file=trim(outputdir) // "emis_limits.dat")
			do i=1,nlam-1
				write(45,*) sqrt(lam(i)*lam(i+1))*1d4,emismin(i),emismax(i)
			enddo
			close(unit=45)

			open(unit=45,file=trim(outputdir) // "emisR_limits.dat")
			do i=1,nlam-1
				write(45,*) sqrt(lam(i)*lam(i+1))*1d4,emisRmin(i),emisRmax(i)
			enddo
			close(unit=45)

			open(unit=45,file=trim(outputdir) // "trans_limits.dat")
			do i=1,nlam-1
				write(45,*) sqrt(lam(i)*lam(i+1))*1d4,obsAmin(i),obsAmax(i)
			enddo
			close(unit=45)
		endif
	endif

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
		case("trans","transmission","transC")
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
			spec(1:ObsSpec(i)%nlam)=log(spec(1:ObsSpec(i)%nlam))
		case("tprofile")
			call ComputeParamT(spec(1:nr))
			spec(1:ObsSpec(i)%nlam)=spec(1:ObsSpec(i)%nlam)*ObsSpec(i)%beta
		case("logtp")
			spec(1)=0d0!(log10(abs(T(1)/T(2))/log10(P(1)/P(2))))
			do j=2,nr-1
				spec(j)=(log10(abs(T(j-1)/T(j)))/log10(P(j-1)/P(j)))-(log10(abs(T(j)/T(j+1)))/log10(P(j)/P(j+1)))
			enddo
			spec(nr)=0d0!(log10(abs(T(nr-1)/T(nr)))/log10(P(nr-1)/P(nr)))
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
	call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam,starfile)
	Fstar=Fstar*pi*Rstar**2
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
	integer i,j,k
	real*8 var(n_ret),dvar(2,n_ret),x,xx

	do k=1,2

	do i=1,n_ret
		if(RetPar(i)%keyword(1:6).eq.'tvalue') then
			read(RetPar(i)%keyword(7:len_trim(RetPar(i)%keyword)),*) j
			xx=10d0**(log10(TP0)+dTP*log10(P(j)))
			x=var(i)
			if(x.gt.0.5) then
				x=(x-0.5)*2d0
				RetPar(i)%value=(RetPar(i)%xmax-xx)*x
			else
				x=(0.5-x)*2d0
				RetPar(i)%value=(RetPar(i)%xmin-xx)*x
			endif
			x=var(i)+dvar(1,i)
			if(x.gt.0.5) then
				x=(x-0.5)*2d0
				RetPar(i)%error1=(RetPar(i)%xmax-xx)*x
				RetPar(i)%error1=RetPar(i)%error1-RetPar(i)%value
			else
				x=(0.5-x)*2d0
				RetPar(i)%error1=(RetPar(i)%xmin-xx)*x
				RetPar(i)%error1=RetPar(i)%error1-RetPar(i)%value
			endif
			x=var(i)+dvar(2,i)
			if(x.gt.0.5) then
				x=(x-0.5)*2d0
				RetPar(i)%error2=(RetPar(i)%xmax-xx)*x
				RetPar(i)%error2=RetPar(i)%error2-RetPar(i)%value
			else
				x=(0.5-x)*2d0
				RetPar(i)%error2=(RetPar(i)%xmin-xx)*x
				RetPar(i)%error2=RetPar(i)%error2-RetPar(i)%value
			endif
		else if(RetPar(i)%logscale) then
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
	Mstar=Mstar/Msun
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

	enddo

	return
	end
	


	
	subroutine MapRetrievalInverse(var)
	use GlobalSetup
	use ReadKeywords
	use Constants
	IMPLICIT NONE
	integer i,j
	real*8 var(n_ret),x,xx

	do i=1,n_ret
		if(RetPar(i)%keyword(1:6).eq.'tvalue') then
			read(RetPar(i)%keyword(7:len_trim(RetPar(i)%keyword)),*) j
			xx=10d0**(log10(TP0)+dTP*log10(P(j)))
			if(RetPar(i)%value.lt.0d0) then
				if(xx.lt.RetPar(i)%xmin) then
					var(i)=0d0
				else
					x=RetPar(i)%value/(RetPar(i)%xmin-xx)
					var(i)=0.5d0-x/2d0
				endif
			else
				if(xx.gt.RetPar(i)%xmax) then
					var(i)=1d0
				else
					x=RetPar(i)%value/(RetPar(i)%xmax-xx)
					var(i)=0.5d0+x/2d0
				endif
			endif
		else if(RetPar(i)%logscale) then
c	log
			var(i)=log10(RetPar(i)%value/RetPar(i)%xmin)/log10(RetPar(i)%xmax/RetPar(i)%xmin)
		else if(RetPar(i)%squarescale) then
c	square
			var(i)=(RetPar(i)%value**2-RetPar(i)%xmin**2)/(RetPar(i)%xmax**2-RetPar(i)%xmin**2)
		else
c	linear
			var(i)=(RetPar(i)%value-RetPar(i)%xmin)/(RetPar(i)%xmax-RetPar(i)%xmin)
		endif
	enddo

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
	

	subroutine WriteRetrieval(imodel,chi2,var,error)
	use GlobalSetup
	IMPLICIT NONE
	real*8 var(n_ret),chi2,error(2,n_ret),sig
	integer i,imodel,j

	call MapRetrieval(var,error)

	open(unit=20,file=trim(outputdir) // "retrieval",RECL=1000)
	write(20,'("Model ",i)') imodel
	write(20,'("chi2=",f14.6)') chi2
	do i=1,n_ret
		if(RetPar(i)%logscale) then
c	log
			sig=log(RetPar(i)%value/RetPar(i)%x0)/RetPar(i)%dx
			write(20,'(a15," = ",es14.7," x/: ",es11.4,es11.4,f9.2)') trim(RetPar(i)%keyword),RetPar(i)%value,
     &					RetPar(i)%error2,RetPar(i)%error1,sig
 		else
c	linear/squared
			sig=(RetPar(i)%value-RetPar(i)%x0)/RetPar(i)%dx
			write(20,'(a15," = ",es14.7," +/- ",es11.4,es11.4,f9.2)') trim(RetPar(i)%keyword),RetPar(i)%value,
     &					RetPar(i)%error2,RetPar(i)%error1,sig
		endif
	enddo
	if(.not.dochemistry) then
		write(20,'(a15," = ",es14.7," +/- ",es11.4,es11.4,f9.2)') 'COratio',COret,COerr(2),COerr(1)
	endif

	close(unit=20)
	
	do i=1,nobs
		select case(ObsSpec(i)%type)
			case("trans","transmission","emisr","emisR","emisa","emis","emission","transC")
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

	
	
	