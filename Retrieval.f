	module RetrievalMod
	implicit none
	integer imodel
	real*8 bestlike
	real*8,allocatable :: dvarq(:),bestvar(:)
	real*8,allocatable :: obsA0(:),obsA1(:),obsA2(:),dobsA(:,:)
	real*8,allocatable :: emis0(:),emis1(:),emis2(:),demis(:,:)
	real*8,allocatable :: emisR0(:),emisR1(:),emisR2(:),demisR(:,:)
	end module RetrievalMod



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
c					if(dy.lt.1d-2*y) dy=1d-2*y
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
	use RetrievalMod
	IMPLICIT NONE
	external ComputeChi2,mrqcomputemodel
	real*8 var0(n_ret),dvar0(2,n_ret),ComputeChi2,dvar(n_ret),x(n_ret),chi2min,random
	real*8,allocatable :: y(:),y1(:),y2(:),dy(:,:),yobs(:),dyobs(:),ybest(:),y0(:)
	real*8 chi2obs(nobs),var(n_ret),chi2,gasdev,maxd,error(2,n_ret),var_best(n_ret),chi2_0
	real*8 x1,x2,minT(nr),maxT(nr),ran1,tot,lambda,chi2_1,chi2_2,dchi2(n_ret),chi2prev
	integer ny,i,j,iter,itermax,iy,k
	
	real*8 maxsig,WLU(n_ret,n_ret),ErrVec(n_ret,n_ret),dvar_prev(n_ret)
	real*8 backup_xmin(n_ret),backup_xmax(n_ret),alphaW(3*n_ret,3*n_ret),Cov(3*n_ret,3*n_ret)
	real*8 beta(n_ret),da(n_ret)
	real*8 phase0(nlam),flux0(nlam),scale,scalemin,Covar(n_ret,n_ret)
	integer na,map(n_ret),info,iboot,nboot,niter1,niter2,n_not_improved,ia(n_ret)
	logical dofit(n_ret),dofit_prev(n_ret),new_best,improved_iter
	logical,allocatable :: specornot(:)

	integer ME,MA,MG,MODE,MDW,ii,nca
	integer,allocatable :: IP(:)
	real*8 PRGOPT(10),RNORME,RNORML
	real*8,allocatable :: WS(:)

	real*8 XeqCloud_best(nr,max(nclouds,1)),mixrat_best_r(nr,nmol)

	external amoebafunk
	real*8 pamoeba(n_ret+1,n_ret),yamoeba(n_ret+1),ftol,amoebafunk

	allocate(obsA0(nlam))
	allocate(obsA1(nlam))
	allocate(obsA2(nlam))
	allocate(dobsA(n_ret,nlam))
	allocate(emis0(nlam))
	allocate(emis1(nlam))
	allocate(emis2(nlam))
	allocate(demis(n_ret,nlam))
	allocate(emisR0(nlam))
	allocate(emisR1(nlam))
	allocate(emisR2(nlam))
	allocate(demisR(n_ret,nlam))

	if(retrievaltype.eq.'MN'.or.retrievaltype.eq.'MultiNest') then
		call doMultiNest
		return
	endif
	
	imodel=0
	bestlike=1d200

	do i=1,n_ret
		RetPar(i)%value=RetPar(i)%x0
	enddo
	call MapRetrievalInverse(var0)
	do i=1,n_ret
		if(var0(i).gt.1d0) var0(i)=1d0
		if(var0(i).lt.0d0) var0(i)=0d0
	enddo

	open(unit=31,file=trim(outputdir) // "Wolk.dat",RECL=6000)

	ny=0
	do i=1,nobs
		ny=ny+ObsSpec(i)%nlam
	enddo
	ny=ny+n_ret
	allocate(y(ny))
	allocate(y0(ny))
	allocate(ybest(ny))
	allocate(y1(ny))
	allocate(y2(ny))
	allocate(dy(n_ret,ny))
	allocate(yobs(ny))
	allocate(dyobs(ny))
	allocate(specornot(ny))
	allocate(dvarq(n_ret))
	allocate(bestvar(n_ret))
	iy=1
	do i=1,nobs
		yobs(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%y(1:ObsSpec(i)%nlam)
		dyobs(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%dy(1:ObsSpec(i)%nlam)*sqrt(real(ObsSpec(i)%nlam))
		specornot(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%spec
		iy=iy+ObsSpec(i)%nlam
	enddo
	do i=1,n_ret
		yobs(iy)=var0(i)
		dyobs(iy)=10d0
		iy=iy+1
	enddo
	
	Cov=0d0

	var=var0
	dvarq=0.02d0

c	goto 2	! skip the amoeba step

	do i=1,n_ret+1
		if(i.eq.1) then
			pamoeba(i,1:n_ret)=var0(1:n_ret)
		else
			do j=1,n_ret
				pamoeba(i,j)=random(idum)
			enddo
		endif
		do j=1,n_ret
			if(pamoeba(i,j).ge.1d0) then
				pamoeba(i,j)=25d0
			else if(pamoeba(i,j).le.0d0) then
				pamoeba(i,j)=-25d0
			else
				pamoeba(i,j)=-log(1d0/pamoeba(i,j)-1d0)
			endif
			if(pamoeba(i,j).gt.25d0) pamoeba(i,j)=25d0
			if(pamoeba(i,j).lt.-25d0) pamoeba(i,j)=-25d0
		enddo
		yamoeba(i)=amoebafunk(pamoeba(i,1:n_ret),ny)
	enddo
	ftol=0.2d0
	i=1
	nca=n_ret+1
	itermax=1000
	call amoeba(pamoeba,yamoeba,nca,n_ret,n_ret,ftol,amoebafunk,iter,ny,itermax)
	var=0d0
	do i=1,n_ret+1
		pamoeba(i,1:n_ret)=1d0/(1d0+exp(-pamoeba(i,1:n_ret)))
		var(1:n_ret)=var(1:n_ret)+pamoeba(i,1:n_ret)/real(n_ret+1)
	enddo
	do i=1,n_ret+1
		dvarq(1:n_ret)=dvarq(1:n_ret)+(var(1:n_ret)+pamoeba(i,1:n_ret))**2
	enddo
	dvarq=sqrt(dvarq/real(n_ret))

2	continue

	ia=1
	nca=3*n_ret
	lambda=-1d0
	n_not_improved=0
	chi2_0=1d200
	do i=1,100
		print*,"Iteration: ",i,chi2
		var=bestvar
		call mrqmin(yobs,dyobs,ny,var,ia,n_ret,Cov,alphaW,nca,chi2,mrqcomputemodel,lambda,beta)
		do j=1,n_ret
			if(var(j).gt.1d0) var(j)=1d0
			if(var(j).lt.0d0) var(j)=0d0
		enddo
		if((chi2_0-chi2).gt.0d0) then
			if((chi2_0-chi2).lt.0.01d0) then
				n_not_improved=n_not_improved+1
			else if((chi2_0-chi2).gt.0.5d0) then
				n_not_improved=0
			endif
			if(n_not_improved.gt.2) exit
		endif

		chi2_0=chi2
	enddo
	lambda=0d0
	call mrqmin(yobs,dyobs,ny,var,ia,n_ret,Cov,alphaW,nca,chi2,mrqcomputemodel,lambda,beta)
	do j=1,n_ret
		if(var(j).gt.1d0) var(j)=1d0
		if(var(j).lt.0d0) var(j)=0d0
	enddo

	call WritePTlimits(var,Cov(1:n_ret,1:n_ret),ErrVec,error,chi2,.true.)
	call WriteRetrieval(imodel,chi2,var,bestvar,error)

	close(unit=31)
	
	return
	end



	
	subroutine mrqcomputemodel(var0,ymod,dyda,nvars,ny)
	use GlobalSetup
	use Constants
	use RetrievalMod
	IMPLICIT NONE
	integer nvars,i,j,nlamtot,ny
	real*8 var(nvars),ymod(ny),dyda(ny,nvars),error(2,nvars),var0(nvars),lnew
	real*8 y1(ny),y2(ny),var1(nvars),var2(nvars),chi2_0,chi2_1,chi2_2,random
	real*8 aq,bq,cq,gasdev,dd
	real*8,allocatable :: spec(:)
	logical recomputeopac

	recomputeopac=.true.

	var=var0
	do i=1,nvars
		if(var(i).gt.1d0) var(i)=1d0
		if(var(i).lt.0d0) var(i)=0d0
	enddo
	call mrqcomputeY(var,ymod,nvars,ny,chi2_0)
	obsA0(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
	emis0(1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
	emisR0(1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar*1d23/distance**2)

c	goto 2	
	do i=1,nvars
3		var1=var
		if(var(i).lt.0.5d0) then
			var1(i)=var1(i)+dvarq(i)*(0.5d0+0.5d0*random(idum))
		else
			var1(i)=var1(i)-dvarq(i)*(0.5d0+0.5d0*random(idum))
		endif
		dd=gasdev(idum)
		var1(i)=var(i)+dd*dvarq(i)
		if(var1(i).lt.0d0) var1(i)=0d0
		if(var1(i).gt.1d0) var1(i)=1d0
		if(abs(var(i)-var1(i)).lt.1d-5) goto 3
		call mrqcomputeY(var1,y1,nvars,ny,chi2_1)
		obsA1(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
		emis1(1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
		emisR1(1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar*1d23/distance**2)
		dyda(1:ny,i)=(y1(1:ny)-ymod(1:ny))/(var1(i)-var(i))
		do j=1,nlam
			dobsA(i,j)=(obsA1(j)-obsA0(j))/(var1(i)-var(i))
			demis(i,j)=(emis1(j)-emis0(j))/(var1(i)-var(i))
			demisR(i,j)=(emisR1(j)-emisR0(j))/(var1(i)-var(i))
		enddo
		dvarq(i)=sqrt(dvarq(i)*(0.2d0*abs(var1(i)-var(i))*chi2_0/abs(dd*(chi2_0-chi2_1))))
		if(dvarq(i).gt.0.1d0) dvarq(i)=0.1d0
		if(dvarq(i).lt.1d-5) dvarq(i)=1d-5
	enddo
	return

2	continue
	do i=1,nvars
1		var1=var
		var2=var
		if(var(i).eq.0d0) then
			var1(i)=var1(i)+dvarq(i)*(0.5d0+0.5d0*random(idum))
			var2(i)=var2(i)+dvarq(i)*(0.5d0+0.5d0*random(idum))
		else if(var(i).eq.1d0) then
			var1(i)=var1(i)-dvarq(i)*(0.5d0+0.5d0*random(idum))
			var2(i)=var2(i)-dvarq(i)*(0.5d0+0.5d0*random(idum))
		else
			var1(i)=var1(i)-dvarq(i)*(0.5d0+0.5d0*random(idum))
			var2(i)=var2(i)+dvarq(i)*(0.5d0+0.5d0*random(idum))
		endif
		if(var1(i).lt.0d0) var1(i)=0d0
		if(var1(i).gt.1d0) var1(i)=1d0
		if(var2(i).lt.0d0) var2(i)=0d0
		if(var2(i).gt.1d0) var2(i)=1d0
		if(var1(i).eq.var(i).or.var2(i).eq.var(i).or.var1(i).eq.var2(i)) goto 1
		call mrqcomputeY(var1,y1,nvars,ny,chi2_1)
		obsA1(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
		emis1(1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
		emisR1(1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar*1d23/distance**2)
		call mrqcomputeY(var2,y2,nvars,ny,chi2_2)
		obsA2(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
		emis2(1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
		emisR2(1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar*1d23/distance**2)
		do j=1,ny
			call quadint(var1(i),var(i),var2(i),y1(j),ymod(j),y2(j),aq,bq,cq)
			dyda(j,i)=2d0*aq*var(i)+bq
		enddo
		do j=1,nlam
			call quadint(var1(i),var(i),var2(i),obsA1(j),obsA0(j),obsA2(j),aq,bq,cq)
			dobsA(i,j)=2d0*aq*var(i)+bq
			call quadint(var1(i),var(i),var2(i),emis1(j),emis0(j),emis2(j),aq,bq,cq)
			demis(i,j)=2d0*aq*var(i)+bq
			call quadint(var1(i),var(i),var2(i),emisR1(j),emisR0(j),emisR2(j),aq,bq,cq)
			demisR(i,j)=2d0*aq*var(i)+bq
		enddo
		dvarq(i)=sqrt(dvarq(i)*(0.2d0*abs(var1(i)-var2(i))*chi2_0/(abs(chi2_0-chi2_1)+abs(chi2_0-chi2_1))))
		if(dvarq(i).gt.0.1d0) dvarq(i)=0.1d0
		if(dvarq(i).lt.1d-5) dvarq(i)=1d-5
	enddo
		
	return
	end

	real*8 function amoebafunk(var_in,ny)
	use GlobalSetup
	use Constants
	use RetrievalMod
	IMPLICIT NONE
	integer ny,i
	real*8 ymod(ny),var(n_ret),chi2,var_in(n_ret)
	real*16 x

	do i=1,n_ret
		x=var_in(i)
		if(x.gt.25.0) then
			var(i)=1d0
		else if(x.lt.-25.0) then
			var(i)=0d0
		else
			var(i)=1d0/(1d0+qexp(-x))
		endif
	enddo
	call mrqcomputeY(var,ymod,n_ret,ny,chi2)
	amoebafunk=chi2
	
	return
	end


	subroutine mrqcomputeY(var,ymod,nvars,ny,lnew)
	use GlobalSetup
	use Constants
	use RetrievalMod
	IMPLICIT NONE
	integer nvars,i,j,nlamtot,ny,k
	real*8 var(nvars),ymod(ny),error(2,nvars),lnew
	real*8,allocatable :: spec(:)
	logical recomputeopac

	recomputeopac=.true.
	imodel=imodel+1
	call output("model number: " // int2string(imodel,'(i7)'))

	error=0d0
	call MapRetrieval(var,error)

	call InitDens()
	call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam,starfile)
	Fstar=Fstar*pi*Rstar**2
	call SetOutputMode(.false.)
	call ComputeModel(recomputeopac)
	call SetOutputMode(.true.)
	
	k=0
	lnew=0d0
	do i=1,nobs
		allocate(spec(ObsSpec(i)%nlam))
		call RemapObs(i,spec)
		do j=1,ObsSpec(i)%nlam
			k=k+1
			ymod(k)=spec(j)
			lnew=lnew+((spec(j)-ObsSpec(i)%y(j))/ObsSpec(i)%dy(j))**2
		enddo
		do j=1,n_ret
			k=k+1
			ymod(k)=var(j)
		enddo
		deallocate(spec)
	enddo
	lnew=lnew/real(k-1)

	write(31,*) imodel,lnew,var(1:nvars),COratio,metallicity
	call flush(31)

	if(lnew.lt.bestlike) then
		call WriteStructure()
		call WriteOutput()

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

		bestlike=lnew
		bestvar=var
	endif
	
	return
	end
	



	subroutine WritePTlimits(var0,Cov,ErrVec,error,chi2,ioflag)
	use GlobalSetup
	use Constants
	use RetrievalMod
	IMPLICIT NONE
	real*8 minT(nr),maxT(nr),var(n_ret),error(2,n_ret),tot,w(n_ret),chi2,ErrVec(n_ret,n_ret)
	real*8 var0(n_ret),random,vec(n_ret),Tbest(nr),gasdev,Cov(n_ret,n_ret)
	real*8 dvar(n_ret),value(n_ret)
	real*8 obsAerr(2,nlam),emiserr(2,nlam),emisRerr(2,nlam)
	integer iobsA(2,nlam),iemis(2,nlam),iemisR(2,nlam)
	real*8 Cinv(n_ret,n_ret),ALU(n_ret,n_ret),max(n_ret)
	integer i,j,k,info,nk,ierr,iter,iCO(2),ivar(2,n_ret),imaxT(nr),iminT(nr),seed
	character*6000 form
	logical ioflag
	seed=42

	call Eigenvalues(Cov,ErrVec,w,n_ret,INFO)

	do i=1,n_ret
		tot=0d0
		do j=1,n_ret
			tot=tot+ErrVec(j,i)**2
		enddo
		if(tot.le.0d0) tot=1d0
		ErrVec(1:n_ret,i)=ErrVec(1:n_ret,i)*sqrt(w(i))/sqrt(tot)
	enddo
	do i=1,n_ret
		max(i)=1d4
		do j=1,n_ret
			tot=1d0/abs(ErrVec(j,i))
			if(tot.lt.max(i)) max(i)=tot
		enddo
	enddo
	
	call SetOutputMode(.false.)
	var=var0
	call MapRetrieval(var,error)
	call InitDens()
	call SetupStructure(.false.)
	Tbest=T
	maxT=0d0
	minT=0d0
	COerr=0d0
	error=0d0
	ivar=0
	imaxT=0
	iminT=0
	iCO=0
	obsAerr=0d0
	iobsA=0d0
	emiserr=0d0
	iemis=0d0
	emisRerr=0d0
	iemisR=0d0

	nk=n_ret*2500
	if(ioflag) then
		open(unit=35,file=trim(outputdir) // "error_cloud.dat",RECL=6000)
		form='("#"' // trim(int2string(n_ret+2,'(i4)')) // 'a19)'
		write(35,form) (trim(int2string(i,'(i4)')) // " " // RetPar(i)%keyword,i=1,n_ret)
     &	,(trim(int2string(n_ret+1,'(i4)')) // " COratio"),(trim(int2string(n_ret+2,'(i4)')) // " metallicity")
		form='(' // int2string(n_ret+2,'(i4)') // 'es19.7)'
	endif
	
	do k=1,nk
		iter=0
1		continue
		vec=0d0
		tot=0d0
		do i=1,n_ret
			call truncated_normal_ab_sample ( 0d0, 1d0, -max(i), max(i), seed, vec(i) )
		enddo
c		vec=vec/sqrt(tot)
c		vec=vec*gasdev(idum)**(1d0/real(n_ret))
		dvar=0d0
		do j=1,n_ret
			do i=1,n_ret
				dvar(i)=dvar(i)+vec(j)*ErrVec(i,j)
			enddo
		enddo
		if(k.eq.1) dvar=0d0
2		continue
		iter=iter+1
		var=var0+dvar
		do i=1,n_ret
			if(var(i).gt.1d0) goto 1
			if(var(i).lt.0d0) goto 1
		enddo
		if(speclimits) then
			obsA1=obsA0
			emis1=emis0
			emisR1=emisR0
			do i=1,n_ret
				obsA1(1:nlam)=obsA1(1:nlam)+dobsA(i,1:nlam)*dvar(i)
				emis1(1:nlam)=emis1(1:nlam)+demis(i,1:nlam)*dvar(i)
				emisR1(1:nlam)=emisR1(1:nlam)+demisR(i,1:nlam)*dvar(i)
			enddo
			do i=1,nlam
				if(obsA1(i).lt.0d0) obsA1(i)=0d0
				if(emis1(i).lt.0d0) emis1(i)=0d0
				if(emisR1(i).lt.0d0) emisR1(i)=0d0
			enddo
		endif

		call MapRetrieval(var,error)
		
		do i=1,n_ret
			value(i)=RetPar(i)%value
			if(RetPar(i)%logscale) value(i)=log10(value(i))
		enddo
		if(ioflag) write(35,form) (value(i),i=1,n_ret),COratio,metallicity
		call InitDens()
		call SetupStructure(.false.)

		if(k.eq.1) then
			Tbest(1:nr)=T(1:nr)
			COratio=COret
		endif
		do i=1,nr
			if(T(i).gt.Tbest(i)) then
				maxT(i)=maxT(i)+(T(i)-Tbest(i))**2
				imaxT(i)=imaxT(i)+1
			endif
			if(T(i).lt.Tbest(i)) then
				minT(i)=minT(i)+(T(i)-Tbest(i))**2
				iminT(i)=iminT(i)+1
			endif
		enddo
		if(speclimits) then
			do i=1,nlam
				if(obsA1(i).lt.obsA0(i)) then
					obsAerr(1,i)=obsAerr(1,i)+(obsA1(i)-obsA0(i))**2
					iobsA(1,i)=iobsA(1,i)+1
				else
					obsAerr(2,i)=obsAerr(2,i)+(obsA1(i)-obsA0(i))**2
					iobsA(2,i)=iobsA(2,i)+1
				endif
				if(emis1(i).lt.emis0(i)) then
					emiserr(1,i)=emiserr(1,i)+(emis1(i)-emis0(i))**2
					iemis(1,i)=iemis(1,i)+1
				else
					emiserr(2,i)=emiserr(2,i)+(emis1(i)-emis0(i))**2
					iemis(2,i)=iemis(2,i)+1
				endif
				if(emisR1(i).lt.emisR0(i)) then
					emisRerr(1,i)=emisRerr(1,i)+(emisR1(i)-emisR0(i))**2
					iemisR(1,i)=iemisR(1,i)+1
				else
					emisRerr(2,i)=emisRerr(2,i)+(emisR1(i)-emisR0(i))**2
					iemisR(2,i)=iemisR(2,i)+1
				endif
			enddo
		endif
		do i=1,n_ret
			if(var(i).lt.var0(i)) then
				error(1,i)=error(1,i)+(var(i)-var0(i))**2
				ivar(1,i)=ivar(1,i)+1
			endif
			if(var(i).gt.var0(i)) then
				error(2,i)=error(2,i)+(var(i)-var0(i))**2
				ivar(2,i)=ivar(2,i)+1
			endif
		enddo
		if(COret.lt.COratio) then
			COerr(1)=COerr(1)+(COret-COratio)**2
			iCO(1)=iCO(1)+1
		endif
		if(COret.gt.COratio) then
			COerr(2)=COerr(2)+(COret-COratio)**2
			iCO(2)=iCO(2)+1
		endif
	enddo
	error(1,1:n_ret)=sqrt(error(1,1:n_ret)/real(ivar(1,1:n_ret)-1))
	error(2,1:n_ret)=sqrt(error(2,1:n_ret)/real(ivar(2,1:n_ret)-1))
	COerr(1:2)=sqrt(COerr(1:2)/real(iCO(1:2)-1))
	minT(1:nr)=sqrt(minT(1:nr)/real(iminT(1:nr)-1))
	maxT(1:nr)=sqrt(maxT(1:nr)/real(imaxT(1:nr)-1))

	call MapRetrieval(var0,error)
	call InitDens()
	call SetupStructure(.false.)

	call SetOutputMode(.true.)
	if(ioflag) then
		close(unit=35)

		open(unit=45,file=trim(outputdir) // "limits.dat",RECL=1000)
		do i=1,nr
			write(45,*) P(i),Tbest(i),minT(i),maxT(i)
		enddo
		close(unit=45)

		if(speclimits) then
			obsAerr=sqrt(obsAerr/real(iobsA-1))
			emiserr=sqrt(emiserr/real(iemis-1))
			emisRerr=sqrt(emisRerr/real(iemisR-1))

			open(unit=45,file=trim(outputdir) // "emis_limits.dat",RECL=1000)
			do i=1,nlam-1
				write(45,*) sqrt(lam(i)*lam(i+1))*1d4,emis0(i),emiserr(1,i),emiserr(2,i)
			enddo
			close(unit=45)

			open(unit=45,file=trim(outputdir) // "emisR_limits.dat",RECL=1000)
			do i=1,nlam-1
				write(45,*) sqrt(lam(i)*lam(i+1))*1d4,emisR0(i),emisRerr(1,i),emisRerr(2,i)
			enddo
			close(unit=45)

			open(unit=45,file=trim(outputdir) // "trans_limits.dat",RECL=1000)
			do i=1,nlam-1
				write(45,*) sqrt(lam(i)*lam(i+1))*1d4,obsA0(i),obsAerr(1,i),obsAerr(2,i)
			enddo
			close(unit=45)
		endif
	endif

	return
	end


	subroutine WritePTlimits_slow(var0,Cov,ErrVec,error,chi2,dobsA,demis,demisR,ioflag)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 minT(nr),maxT(nr),var(n_ret),error(2,n_ret),tot,w(n_ret),chi2,ErrVec(n_ret,n_ret)
	real*8 var0(n_ret),random,vec(n_ret),Tbest(nr),gasdev,Cov(n_ret,n_ret)
	real*8 dvar(n_ret),dobsA(n_ret,nlam),demis(n_ret,nlam),demisR(n_ret,nlam),value(n_ret)
	real*8 obsAmin(nlam),obsAmax(nlam),emismin(nlam),emismax(nlam),emisRmin(nlam),emisRmax(nlam)
	real*8 emis0(nlam),obsA0(nlam),emisR0(nlam),Cinv(n_ret,n_ret),ALU(n_ret,n_ret)
	integer i,j,k,info,nk,ierr,iter,iCO(2),ivar(2,n_ret),imaxT(nr),iminT(nr)
	character*6000 form
	logical ioflag


	call MatrixInvert(Cov,Cinv,ALU,n_ret,INFO)
	
	call SetOutputMode(.false.)
	var=var0
	call MapRetrieval(var,error)
	call InitDens()
	call SetupStructure(.false.)
	Tbest=T
	maxT=0d0
	minT=0d0
	obsAmin=1d200
	emismin=1d200
	emisRmin=1d200
	obsAmax=0d0
	emismax=0d0
	emisRmax=0d0
	COerr=0d0
	error=0d0
	ivar=0
	imaxT=0
	iminT=0
	iCO=0

	nk=n_ret*1000
	if(ioflag) then
		open(unit=35,file=trim(outputdir) // "error_cloud.dat",RECL=6000)
		form='("#"' // trim(int2string(n_ret,'(i4)')) // 'a19)'
		write(35,form) (trim(int2string(i,'(i4)')) // " " // RetPar(i)%keyword,i=1,n_ret)
		form='(' // int2string(n_ret,'(i4)') // 'es19.7)'
	endif
	
	do k=1,nk
		iter=0
1		continue
		dvar=0d0
		do i=1,n_ret
			var(i)=random(idum)
		enddo
		if(k.eq.1) var=var0
2		continue
		iter=iter+1
		dvar=var0-var
		vec=0d0
		do i=1,n_ret
			do j=1,n_ret
				vec(j)=vec(j)+dvar(i)*Cinv(i,j)
			enddo
		enddo
		tot=0d0
		do i=1,n_ret
			tot=tot+vec(i)*dvar(i)
		enddo
		if(tot.gt.abs(gasdev(idum))) goto 1

		if(speclimits) then
			obsA0=0d0
			emis0=0d0
			emisR0=0d0
		endif
		do i=1,n_ret
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
		
		do i=1,n_ret
			value(i)=RetPar(i)%value
			if(RetPar(i)%logscale) value(i)=log10(value(i))
		enddo
		if(ioflag) write(35,form) (value(i),i=1,n_ret)
		call InitDens()
		call SetupStructure(.false.)

		if(k.eq.1) then
			Tbest(1:nr)=T(1:nr)
			COratio=COret
		endif
		do i=1,nr
			if(T(i).gt.Tbest(i)) then
				maxT(i)=maxT(i)+(T(i)-Tbest(i))**2
				imaxT(i)=imaxT(i)+1
			endif
			if(T(i).lt.Tbest(i)) then
				minT(i)=minT(i)+(T(i)-Tbest(i))**2
				iminT(i)=iminT(i)+1
			endif
		enddo
		do i=1,n_ret
			if(var(i).lt.var0(i)) then
				error(1,i)=error(1,i)+(var(i)-var0(i))**2
				ivar(1,i)=ivar(1,i)+1
			endif
			if(var(i).gt.var0(i)) then
				error(2,i)=error(2,i)+(var(i)-var0(i))**2
				ivar(2,i)=ivar(2,i)+1
			endif
		enddo
		if(COret.lt.COratio) then
			COerr(1)=COerr(1)+(COret-COratio)**2
			iCO(1)=iCO(1)+1
		endif
		if(COret.gt.COratio) then
			COerr(2)=COerr(2)+(COret-COratio)**2
			iCO(2)=iCO(2)+1
		endif
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
	error(1,1:n_ret)=sqrt(error(1,1:n_ret)/real(ivar(1,1:n_ret)-1))
	error(2,1:n_ret)=sqrt(error(2,1:n_ret)/real(ivar(2,1:n_ret)-1))
	COerr(1:2)=sqrt(COerr(1:2)/real(iCO(1:2)-1))
	minT(1:nr)=sqrt(minT(1:nr)/real(iminT(1:nr)-1))
	maxT(1:nr)=sqrt(maxT(1:nr)/real(imaxT(1:nr)-1))

	call MapRetrieval(var0,error)
	call InitDens()
	call SetupStructure(.false.)

	call SetOutputMode(.true.)
	if(ioflag) then
		close(unit=35)

		open(unit=45,file=trim(outputdir) // "limits.dat")
		do i=1,nr
			write(45,*) P(i),Tbest(i),minT(i),maxT(i)
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
	if(gamma_equal) gammaT2=gammaT1
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
	integer i1,i0,n0,n1,i,nn0
	real*8 x0(n0),y0(n0),x1(n1),y1(n1),R1(n1),expR1(n1),w,tot
	real*8 xx0(n0+n1),yy0(n0+n1)

	do i=1,n1
		xx0(n0+i)=x1(i)
	enddo
	nn0=n0+n1
	call sort(xx0,nn0)

	call regridarray(x0,y0,n0,xx0,yy0,nn0)

	do i1=1,n1
		if(x1(i1).lt.xx0(1)) then
			y1(i1)=yy0(1)
		else if(x1(i1).gt.xx0(nn0)) then
			y1(i1)=yy0(nn0)
		else
			y1(i1)=0d0
			tot=0d0
			do i0=1,nn0
				w=exp(-abs((xx0(i0)-x1(i1))*R1(i1)*2d0/x1(i1))**expR1(i1))
				if(i0.eq.1) then
					w=w*abs(xx0(2)-xx0(1))/2d0
				else if(i0.eq.nn0) then
					w=w*abs(xx0(nn0)-xx0(nn0-1))/2d0
				else
					w=w*abs(xx0(i0-1)-xx0(i0+1))/2d0
				endif
				tot=tot+w
				y1(i1)=y1(i1)+w*yy0(i0)
			enddo
			y1(i1)=y1(i1)/tot
		endif
	enddo

	return
	end
	

	subroutine WriteRetrieval(imodel,chi2,var,bestvar,error)
	use GlobalSetup
	IMPLICIT NONE
	real*8 var(n_ret),chi2,error(2,n_ret),sig,bestvar(n_ret)
	integer i,imodel,j
	character*500 command

	call MapRetrieval(var,error)

98	open(unit=20,file=trim(outputdir) // "retrieval",RECL=1000,ERR=99)
	goto 100
99	write(command,'("mkdir -p ",a)') trim(outputdir)
	call system(command)
	call sleep(10)
	goto 98
100	continue
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
		write(20,'(a15," = ",es14.7," +/- ",es11.4,es11.4)') 'COratio',COret,COerr(2),COerr(1)
	endif
	if(Tform.gt.0d0) then
		if(domakeai.and..not.retrieval) then
			write(20,'(a15," = ",es14.7," +/- ",es11.4,es11.4)') 'Tform',Tform,0d0,0d0
			write(20,'(a15," = ",es14.7," +/- ",es11.4,es11.4)') 'enrich',f_enrich,0d0,0d0
		endif
		write(20,'(a15," = ",es14.7," +/- ",es11.4,es11.4)') 'COratio',COratio,0d0,0d0
		write(20,'(a15," = ",es14.7," +/- ",es11.4,es11.4)') 'metallicity',metallicity,0d0,0d0
	endif
	if(mapCOratio) then
		do i=1,nmol
			if(mixrat_r(1,i).gt.0d0) write(20,'(a15," = ",es14.7)') trim(molname(i)),mixrat_r(1,i)
		enddo
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

	call MapRetrieval(bestvar,error)

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

	
	
	