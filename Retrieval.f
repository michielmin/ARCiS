	subroutine ReadObs()
	use GlobalSetup
	IMPLICIT NONE
	integer i,j,ilam,nj,k,it
	real*8 x,y,dy,specres_obs,expspecres_obs
	character*6000 line
	real*8 scale,scale_av,d,dmin,maxsig,minsig
	integer nscale,i2,j2
	logical truefalse

	if(useobsgrid) then
		lamemis=.false.
		lamtrans=.false.
	endif
	if((computeT.or.doRing).and.useobsgrid) lamemis=RTgridpoint

	do i=1,nobs
		select case(ObsSpec(i)%type)
			case('tprofile','logtp')
				ObsSpec(i)%ndata=nr
				ObsSpec(i)%nlam=0
				allocate(ObsSpec(i)%lam(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%y(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%dy(ObsSpec(i)%ndata))
				ObsSpec(i)%y=0d0
				ObsSpec(i)%dy=1d0
				ObsSpec(i)%spec=.false.
			case('priors','prior')
				ObsSpec(i)%ndata=n_ret
				ObsSpec(i)%nlam=0
				allocate(ObsSpec(i)%lam(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%y(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%dy(ObsSpec(i)%ndata))
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
			case('lightcurve')
				ObsSpec(i)%spec=.true.
				inquire(file=ObsSpec(i)%file,exist=truefalse)
				if(.not.truefalse) then
					call output("File does not exist" // trim(ObsSpec(i)%file))
					stop
				endif
				open(unit=20,file=ObsSpec(i)%file,FORM="FORMATTED",ACCESS="STREAM")
10				read(20,*,err=10) ObsSpec(i)%nt
				allocate(ObsSpec(i)%t(ObsSpec(i)%nt))
				allocate(ObsSpec(i)%dt(ObsSpec(i)%nt))
				do j=1,ObsSpec(i)%nt
					read(20,*) ObsSpec(i)%t(j)!,ObsSpec(i)%dt(j)
				enddo
				j=1
				ilam=1
11				read(20,*,end=12) x
				read(20,*)
				read(20,*)
				if((x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)).or.useobsgrid) ilam=ilam+1
				j=j+1
				goto 11
12				ObsSpec(i)%nlam=ilam-1
				ObsSpec(i)%ndata=ObsSpec(i)%nlam*ObsSpec(i)%nt
				nj=j-1
				rewind(20)
20				read(20,*,err=20) ObsSpec(i)%nt
				do j=1,ObsSpec(i)%nt
					read(20,*) ObsSpec(i)%t(j)!,ObsSpec(i)%dt(j)
				enddo
				allocate(ObsSpec(i)%lam(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%y(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%dy(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%R(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%Rexp(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%model(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%LC(ObsSpec(i)%nlam,ObsSpec(i)%nt))
				allocate(ObsSpec(i)%dLC(ObsSpec(i)%nlam,ObsSpec(i)%nt))
				ilam=1
				if(ObsSpec(i)%beta.lt.0d0) ObsSpec(i)%beta=ObsSpec(i)%ndata
				k=0
				print*,nj
				do j=1,nj
13					read(20,'(a6000)') line
					specres_obs=specres
					expspecres_obs=20d0
					read(line,*,err=13,end=14) x,specres_obs,expspecres_obs
14					continue
					if((x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)).or.useobsgrid) then
						ObsSpec(i)%lam(ilam)=x*1d-4
						ObsSpec(i)%R(ilam)=specres_obs
						ObsSpec(i)%Rexp(ilam)=expspecres_obs
						read(20,*) ObsSpec(i)%LC(ilam,1:ObsSpec(i)%nt)
						read(20,*) ObsSpec(i)%dLC(ilam,1:ObsSpec(i)%nt)
						ilam=ilam+1
					else
						read(20,*)
						read(20,*)
					endif
				enddo
				close(unit=20)
				k=0
				do it=1,ObsSpec(i)%nt
					do ilam=1,ObsSpec(i)%nlam
						k=k+1
						ObsSpec(i)%y(k)=ObsSpec(i)%LC(ilam,it)
						if(ObsSpec(i)%LC(ilam,it).eq.1d0) then
							ObsSpec(i)%dy(k)=ObsSpec(i)%dLC(ilam,it)/1000d0
						else
							ObsSpec(i)%dy(k)=ObsSpec(i)%dLC(ilam,it)/100d0
						endif
					enddo
				enddo
			case default
				ObsSpec(i)%spec=.true.
				inquire(file=ObsSpec(i)%file,exist=truefalse)
				if(.not.truefalse) then
					call output("File does not exist" // trim(ObsSpec(i)%file))
					stop
				endif
				open(unit=20,file=ObsSpec(i)%file)
				j=1
				ilam=1
1				read(20,*,end=2,err=1) x,y,dy
				if((x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)).or.useobsgrid) ilam=ilam+1
				j=j+1
				goto 1
2				ObsSpec(i)%ndata=ilam-1
				ObsSpec(i)%nlam=ObsSpec(i)%ndata
				nj=j-1
				close(unit=20)
				open(unit=20,file=ObsSpec(i)%file,FORM="FORMATTED",ACCESS="STREAM")
				allocate(ObsSpec(i)%lam(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%y(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%dy(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%R(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%Rexp(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%model(ObsSpec(i)%ndata))
				allocate(ObsSpec(i)%ilam(ObsSpec(i)%ndata))
				ilam=1
				if(ObsSpec(i)%beta.lt.0d0) ObsSpec(i)%beta=ObsSpec(i)%ndata
				do j=1,nj
3					read(20,'(a1000)',err=3) line
					specres_obs=specres
					expspecres_obs=20d0
					read(line,*,err=3,end=4) x,y,dy,specres_obs,expspecres_obs
4					continue
					if((x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)).or.useobsgrid) then
						ObsSpec(i)%lam(ilam)=x*1d-4
						if(ObsSpec(i)%type.eq."emisa".or.ObsSpec(i)%type.eq."emis".or.ObsSpec(i)%type.eq."emission"
     &								.or.ObsSpec(i)%type.eq."phase") then
							if(log_emis) then
								dy=dy/y
								if(y.le.0d0) then
									call output("negative flux values while using logemis=.true.")
									call output("please use: logemis=.false.")
									call output("stopping...")
									stop
								endif
								ObsSpec(i)%y(ilam)=log(y)
							else
								ObsSpec(i)%y(ilam)=y
							endif
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
				do j=1,ObsSpec(i)%ndata
					dmin=1d200
					do ilam=1,nlam
						d=abs(lam(ilam)-ObsSpec(i)%lam(j))/dlam(ilam)
						if(abs(ObsSpec(i)%R(j)-lam(ilam)/dlam(ilam))/ObsSpec(i)%R(j).lt.0.1.and.d.lt.dmin.and..not.RTgridpoint(ilam)) then
							ObsSpec(i)%ilam(j)=ilam
							dmin=d
						endif
					enddo
					if(useobsgrid) then
						if(ObsSpec(i)%type(1:4).eq.'emis'.or.ObsSpec(i)%type(1:5).eq."phase") lamemis(ObsSpec(i)%ilam(j))=.true.
						if(ObsSpec(i)%type(1:5).eq.'trans') lamtrans(ObsSpec(i)%ilam(j))=.true.
					endif
				enddo
		end select
		if(ObsSpec(i)%filter.ne.' ') then
			if(ObsSpec(i)%type(1:4).eq.'emis'.or.ObsSpec(i)%type(1:5).eq."phase") lamemis=.true.
			if(ObsSpec(i)%type(1:5).eq.'trans') lamtrans=.true.
			allocate(ObsSpec(i)%f(nlam,ObsSpec(i)%ndata))
			ObsSpec(i)%f=0d0
			call regridNsimple(ObsSpec(i)%filter,lam*1d4,ObsSpec(i)%f,nlam,1,3,ObsSpec(i)%ndata)
		endif
	enddo

	if(faircoverage) then
		scale_av=0d0
		nscale=0
		do i=1,nobs
			if(ObsSpec(i)%spec) then
				do j=1,ObsSpec(i)%ndata
					scale=0d0
					do i2=1,nobs
						if(ObsSpec(i2)%spec) then
							do j2=1,ObsSpec(i2)%ndata
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
				do j=1,ObsSpec(i)%ndata
					ObsSpec(i)%dy(j)=ObsSpec(i)%dy(j)/scale_av
				enddo
			endif
		enddo
	endif

	if(retrieval) then
		do i=1,n_ret
			if(RetPar(i)%keyword(1:9).eq.'model_err') then
				read(RetPar(i)%keyword(14:len(trim(RetPar(i)%keyword))),*) it
				if(it.eq.1) then
					x=0d0
				else
					x=model_err_lam(it-1)
				endif
				if(it.eq.nmodel_err) then
					y=1d200
				else
					y=model_err_lam(it)
				endif
				maxsig=0d0
				minsig=1d200
				do k=1,nobs
					do j=1,ObsSpec(k)%ndata
						if(ObsSpec(k)%lam(j).ge.x.and.ObsSpec(k)%lam(j).le.y) then
							if(ObsSpec(k)%dy(j).gt.maxsig) maxsig=ObsSpec(k)%dy(j)
							if(ObsSpec(k)%dy(j).lt.minsig) minsig=ObsSpec(k)%dy(j)
						endif
					enddo
				enddo
				RetPar(i)%xmin=0.1*minsig
				RetPar(i)%xmax=10d0*maxsig
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
	external mrqcomputemodel
	real*8 var0(n_ret),dvar0(2,n_ret),dvar(n_ret),x(n_ret),chi2min,random
	real*8,allocatable :: y(:),y1(:),y2(:),dy(:,:),yobs(:),dyobs(:),ybest(:),y0(:)
	real*8 chi2obs(nobs),var(n_ret),chi2,gasdev,maxd,error(2,n_ret),var_best(n_ret)
	real*8 x1,x2,minT(nr),maxT(nr),ran1,tot,lambda,chi2_1,chi2_2,dchi2(n_ret),chi2prev
	integer ny,i,j,iter,itermax,iy,k
	
	real*8 maxsig,WLU(n_ret,n_ret),ErrVec(n_ret,n_ret),dvar_prev(n_ret)
	real*8 backup_xmin(n_ret),backup_xmax(n_ret),alphaW(3*n_ret,3*n_ret),Cov(3*n_ret,3*n_ret)
	real*8 beta(n_ret),da(n_ret)
	real*8 phase0(nlam),flux0(nlam),scale,scalemin,Covar(n_ret,n_ret)
	integer na,map(n_ret),info,niter1,niter2,n_not_improved,ia(n_ret)
	logical dofit(n_ret),dofit_prev(n_ret),new_best,improved_iter
	logical,allocatable :: specornot(:)

	integer ME,MA,MG,MODE,MDW,ii,nca,jj
	integer,allocatable :: IP(:)
	real*8 PRGOPT(10),RNORME,RNORML
	real*8,allocatable :: WS(:)

	real*8 XeqCloud_best(nr,max(nclouds,1)),mixrat_best_r(nr,nmol)

	external amoebafunk,lmcompute,MCMCfunc
	real*8 pamoeba(n_ret+1,n_ret),yamoeba(n_ret+1),ftol,amoebafunk,MCMCfunc

	integer iboot
	real*8,allocatable :: chi2_boot(:),count(:)
	real*8 chi2_boot_av,chi2_boot_sig1,chi2_boot_sig2
	
	if(.not.allocated(obsA0)) then
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
		allocate(bestvar(n_ret))
	endif

	if(retrievaltype.eq.'MN'.or.retrievaltype.eq.'MultiNest') then
		call doMultiNest
		return
	endif
	
	imodel=0
	bestlike=-1d200

	if(randomstart) then
		do i=1,n_ret
			var0(i)=random(idum)
		enddo
	else
		do i=1,n_ret
			RetPar(i)%value=RetPar(i)%x0
		enddo
		call MapRetrievalInverse(var0)
	endif
	call fold(var0,var0,n_ret)
	do i=1,n_ret
		if(var0(i).gt.1d0) var0(i)=1d0
		if(var0(i).lt.0d0) var0(i)=0d0
	enddo
c	do i=1,n_ret
c		var0(i)=-log(1d0/var0(i)-1d0)
c		if(var0(i).gt.25d0) var0(i)=25d0
c		if(var0(i).lt.-25d0) var0(i)=-25d0
c	enddo

	
	if(writeWolk) open(unit=31,file=trim(outputdir) // "Wolk.dat",FORM="FORMATTED",ACCESS="STREAM")

	ny=0
	do i=1,nobs
		ny=ny+ObsSpec(i)%ndata
	enddo
	allocate(y(ny))
	allocate(y0(ny))
	allocate(ybest(ny))
	allocate(y1(ny))
	allocate(y2(ny))
	allocate(count(ny))
	allocate(dy(n_ret,ny))
	allocate(yobs(ny))
	allocate(dyobs(ny))
	allocate(specornot(ny))
	allocate(dvarq(n_ret))

	if(retrievaltype.eq.'MC'.or.retrievaltype.eq.'MCMC') then
		call doMCMCF90(var0,n_ret)
c		call MCMC(MCMCfunc,var0,n_ret,npop,npop*100,ny)
		return
	endif
	
	allocate(chi2_boot(nboot))
	open(unit=90,file='chi2boot',FORM="FORMATTED",ACCESS="STREAM")
	do iboot=1,nboot
	imodel=0
	bestlike=-1d200

	iy=1
	do i=1,nobs
		yobs(iy:iy+ObsSpec(i)%ndata-1)=ObsSpec(i)%y(1:ObsSpec(i)%ndata)
		dyobs(iy:iy+ObsSpec(i)%ndata-1)=ObsSpec(i)%dy(1:ObsSpec(i)%ndata)
		specornot(iy:iy+ObsSpec(i)%ndata-1)=ObsSpec(i)%spec
		iy=iy+ObsSpec(i)%ndata
	enddo
	if(iboot.ne.nboot) then
c		count=1d-6
c		do i=1,ny
c			j=random(idum)*real(ny)+1
c			count(j)=count(j)+1d0
c		enddo
c		dyobs=dyobs/sqrt(count)
		do i=1,ny
			yobs(i)=yobs(i)+gasdev(idum)*dyobs(i)
			dyobs(i)=dyobs(i)*sqrt(2d0)
		enddo
	endif
	
	Cov=0d0

	var=var0
	dvarq=0.02d0

	jj=0
	if(retrievaltype.eq.'OE'.or.retrievaltype.eq.'oe') goto 2	! skip the amoeba step
1	continue
	jj=jj+1

	do i=1,n_ret+1
		if(i.eq.1.and.jj.eq.1) then
			pamoeba(i,1:n_ret)=var(1:n_ret)
		else
			pamoeba(i,1:n_ret)=var(1:n_ret)
			do j=1,n_ret
				pamoeba(i,j)=random(idum)
			enddo
		endif
		yamoeba(i)=amoebafunk(pamoeba(i,1:n_ret),ny)
	enddo
	ftol=0.2d0
	i=1
	nca=n_ret+1
	itermax=1000
	call amoeba(pamoeba,yamoeba,nca,n_ret,n_ret,ftol,amoebafunk,iter,ny,itermax)
	var=0d0
	do i=1,n_ret+1
		var(1:n_ret)=var(1:n_ret)+pamoeba(i,1:n_ret)/real(n_ret+1)
	enddo
	do i=1,n_ret+1
		dvarq(1:n_ret)=dvarq(1:n_ret)+(var(1:n_ret)+pamoeba(i,1:n_ret))**2
	enddo
	dvarq=sqrt(dvarq/real(n_ret))

2	continue

	do ii=1,3

	if(jj.eq.1.or.ii.gt.1) then
		if(jj.eq.1.and.ii.eq.1) then
			var=bestvar
		else
			if(chi2.lt.1d4) then
				var=bestvar
			else
				do i=1,n_ret
					var(i)=random(idum)
				enddo
			endif
		endif
	endif

	ia=1
	nca=3*n_ret
	lambda=-1d0
	n_not_improved=0
	chi2prev=1d200

	do i=1,100
c		print*,"Iteration: ",iboot,ii,i,chi2
		call mrqmin(yobs,dyobs,ny,var,ia,n_ret,Cov,alphaW,nca,chi2,mrqcomputemodel,lambda,beta)
		chi2=chi2/real(max(1,ny-n_ret))
		if((chi2prev-chi2).gt.0d0) then
			if((chi2prev-chi2).lt.0.01d0) then
				n_not_improved=n_not_improved+1
			else if((chi2prev-chi2).gt.0.5d0) then
				n_not_improved=0
			endif
			if(n_not_improved.gt.2) exit
		endif
		chi2prev=chi2
	enddo

	enddo

	lambda=-1d0
	call mrqmin(yobs,dyobs,ny,var,ia,n_ret,Cov,alphaW,nca,chi2,mrqcomputemodel,lambda,beta)
	lambda=0d0
	call mrqmin(yobs,dyobs,ny,var,ia,n_ret,Cov,alphaW,nca,chi2,mrqcomputemodel,lambda,beta)
	chi2=chi2/real(max(1,ny-n_ret))
	call fold(var,var,n_ret)
	do j=1,n_ret
		if(var(j).gt.1d0) var(j)=1d0
		if(var(j).lt.0d0) var(j)=0d0
	enddo

	chi2_boot(iboot)=chi2

	if(iboot.gt.1.and.iboot.lt.nboot) then
		call stati(chi2_boot(1:iboot),iboot,chi2_boot_av,chi2_boot_sig1,chi2_boot_sig2)
		write(90,*) chi2_boot_av,chi2_boot_sig1,chi2_boot_sig2
		flush(90)
	endif

	enddo

	call output("chi2 best:    " // trim(dbl2string(bestchi2,'(f10.3)')))
	if(nboot.gt.1) then
		call stati(chi2_boot(1:nboot-1),nboot-1,chi2_boot_av,chi2_boot_sig1,chi2_boot_sig2)
		call output("chi2 average: " // trim(dbl2string(chi2_boot_av,'(f10.3)')))
		call output("chi2 sigma-:  " // trim(dbl2string(chi2_boot_sig1,'(f10.3)')))
		call output("chi2 sigma+:  " // trim(dbl2string(chi2_boot_sig2,'(f10.3)')))
	endif
	call WritePTlimits(var,Cov(1:n_ret,1:n_ret),ErrVec,error,bestchi2,.true.)
	call WriteRetrieval(imodel,chi2,var,bestvar,error)

	if(writeWolk) close(unit=31)
	
	return
	end

	subroutine fold(var0,var1,nvar)
	IMPLICIT NONE
	integer nvar,i
	real*8 var0(nvar),var1(nvar)
	
	var1=var0
	do i=1,nvar
		if(.not.var1(i).lt.1d0) var1(i)=1d0
		if(.not.var1(i).gt.0d0) var1(i)=0d0
	enddo
	return
	
	do i=1,nvar
1		continue
		if(.not.var1(i).gt.0.5d0.and..not.var1(i).lt.0.5d0) var1(i)=0.5d0
		if(var1(i).gt.1d0) then
			var1(i)=2d0-var1(i)
			goto 1
		endif
		if(var1(i).lt.0d0) then
			var1(i)=-var1(i)
			goto 1
		endif
	enddo

	return
	end


	
	subroutine mrqcomputemodel(var0,ymod,dyda,nvars,ny,what)
	use GlobalSetup
	use Constants
	use RetrievalMod
	IMPLICIT NONE
	integer nvars,i,j,nlamtot,ny,what,ii
	real*8 var(nvars),ymod(ny),dyda(ny,nvars),error(2,nvars),var0(nvars),lnew
	real*8 y1(ny),y2(ny),var1(nvars),var2(nvars),chi2_1,chi2_2,random
	real*8 aq,bq,cq,gasdev,dd,scale
	real*8,allocatable :: spec(:)
	logical recomputeopac,recompute

	recomputeopac=.true.

	do i=1,nvars
		if(var0(i).gt.1d0) var0(i)=1d0
		if(var0(i).lt.0d0) var0(i)=0d0
	enddo

	var=var0
	if(what.eq.1) then
		call mrqcomputeY(var,ymod,nvars,ny,chi2_0,scale)
		chi2_0=global_chi2
		obsA0(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
		emis0(1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
		emisR0(1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar*1d23/distance**2)
	endif

	if(what.eq.2) then
	do i=1,nvars
		ii=0
3		ii=ii+1
		if(dvarq(i).gt.0.1d0) dvarq(i)=0.1d0
		if(dvarq(i).lt.1d-5) dvarq(i)=1d-5
4		var1=var
		dd=gasdev(idum)
		var1(i)=var(i)+dd*dvarq(i)
		if(var1(i).gt.1d0) goto 4
		if(var1(i).lt.0d0) goto 4
		if(abs(var(i)-var1(i)).lt.1d-5) goto 3
		call mrqcomputeY(var1,y1,nvars,ny,chi2_1,scale)
		chi2_1=global_chi2
		obsA1(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
		emis1(1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
		emisR1(1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar*1d23/distance**2)
		dyda(1:ny,i)=(y1(1:ny)-ymod(1:ny))/(var1(i)-var(i))
		do j=1,nlam
			dobsA(i,j)=(obsA1(j)-obsA0(j))/(var1(i)-var(i))
			demis(i,j)=(emis1(j)-emis0(j))/(var1(i)-var(i))
			demisR(i,j)=(emisR1(j)-emisR0(j))/(var1(i)-var(i))
		enddo
		if(abs(chi2_0-chi2_1).lt.1d-4.and.ii.lt.3.and.dd*dvarq(i).lt.0.5d0) then
			dvarq(i)=dvarq(i)*5d0
			if(dvarq(i).gt.0.5d0) dvarq(i)=0.5d0
			goto 3
		endif
		dvarq(i)=sqrt(dvarq(i)*(0.2d0*abs(var1(i)-var(i))*chi2_0/abs(chi2_0-chi2_1)))
		if(dvarq(i).gt.0.5d0) dvarq(i)=0.5d0
		if(dvarq(i).lt.1d-5) dvarq(i)=1d-5
		if(abs(chi2_0-chi2_1).gt.(chi2_0/2d0).and.ii.lt.3) goto 3
	enddo
	endif	
		
	return
	end

	real*8 function amoebafunk(var,ny)
	use GlobalSetup
	use Constants
	use RetrievalMod
	IMPLICIT NONE
	integer ny,i
	real*8 ymod(ny),var(n_ret),chi2,scale

	call mrqcomputeY(var,ymod,n_ret,ny,chi2,scale)
	amoebafunk=-chi2
	
	return
	end

	subroutine mrqcomputeY(var_in,ymod,nvars,ny,lnew,scale)
	use GlobalSetup
	use Constants
	use RetrievalMod
	use Struct3D
	IMPLICIT NONE
	integer nvars,i,j,nlamtot,ny,k,maxspec,im,ilam,status,system,ii
	real*8 var(nvars),ymod(ny),error(2,nvars),lnew,var_in(nvars),spectemp(nlam),specsave(nobs,nlam)
	real*8,allocatable :: spec(:),allspec(:,:)
	logical recomputeopac,truefalse,doscaleR2
	real*8 xx,xy,scale,dy(ny),tot
	character*100 command

	doscaleR2=doscaleR
	recomputeopac=.true.
	if(dopostequalweights) doscaleR2=.false.
	if(imodel.gt.nscaleR.and.nscaleR.gt.0) doscaleR2=.false.
	imodel=imodel+1
	if(.not.useobsgrid.or.100*(imodel/100).eq.imodel.or.(do3D.and.night2day.ne.1d0.and.n3D.gt.2)) call output("model number: " 
     &				// int2string(imodel,'(i7)') // dbl2string(bestchi2,'(f10.2)'))

	var=var_in
	call fold(var_in,var,n_ret)
2	do i=1,n_ret
		if(var(i).gt.1d0) var(i)=1d0
		if(var(i).lt.0d0) var(i)=0d0
	enddo

	maxspec=0
	do i=1,nobs
		if(ObsSpec(i)%ndata.gt.maxspec) maxspec=ObsSpec(i)%ndata
	enddo
	allocate(allspec(nobs,maxspec))

	allspec=0d0
	if(n2d.eq.0) then
		i2d=0
	else
		i2d=1
	endif

1	continue
	error=0d0
	call SetOutputMode(.false.)
	call MapRetrieval(var,error)
	call SetOutputMode(.true.)

	call InitDens()
	if(retrievestar) then
		if(standardstarname.ne.' ') call StandardStar(standardstarname,Tstar,Rstar,Mstar)
		call StarSpecSetup(Tstar,logg,1d4*lam,1d4*blam,Fstar,nlam,starfile,blackbodystar)
		Fstar=Fstar*pi*Rstar**2
	endif
	call SetOutputMode(.false.)
	call ComputeModel(recomputeopac)
	call SetOutputMode(.true.)
	
	do i=1,nobs
		allocate(spec(ObsSpec(i)%ndata))
		call RemapObs(i,spec,spectemp)
		if(ObsSpec(i)%i2d.eq.i2d) then
			allspec(i,1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)
			specsave(i,1:nlam)=spectemp(1:nlam)
		endif
		deallocate(spec)
	enddo
	i2d=i2d+1
	if(i2d.le.n2d) goto 1

	k=0
	do i=1,nobs
		do j=1,ObsSpec(i)%ndata
			k=k+1
			dy(k)=sqrt(ObsSpec(i)%dy(j)**2+ObsSpec(i)%adderr**2)
			do ii=1,nmodel_err-1
				if(ObsSpec(i)%lam(j).lt.model_err_lam(ii)) exit
			enddo
			select case(ObsSpec(i)%type)
				case("emisa","emis","emission","phase")
					dy(k)=sqrt(dy(k)**2+model_err_abs(ii)**2)
				case("trans","transmission","emisr","emisR","transC","phaser","phaseR","transM","transE")
					dy(k)=sqrt(dy(k)**2+model_err_rel(ii)**2)
			end select
		enddo
	enddo
	if(doscaleR2) then
		xy=0d0
		xx=0d0
		k=0
		do i=1,nobs
			do j=1,ObsSpec(i)%ndata
				k=k+1
				ymod(k)=allspec(i,j)
				if(.not.ObsSpec(i)%scaling) then
					xy=xy+ymod(k)*ObsSpec(i)%y(j)/dy(k)**2
					xx=xx+ymod(k)*ymod(k)/dy(k)**2
				endif
			enddo
		enddo
		scale=1d0
		if(xx.gt.0d0) scale=xy/xx
		do i=1,n_ret
			if(RetPar(i)%keyword.eq.'Rp'.or.RetPar(i)%keyword.eq.'rp') then
				RetPar(i)%value=RetPar(i)%value*sqrt(scale)
				if(.not.RetPar(i)%value.gt.RetPar(i)%xmin) RetPar(i)%value=RetPar(i)%xmin
				if(RetPar(i)%value.gt.RetPar(i)%xmax) RetPar(i)%value=RetPar(i)%xmax
				if(RetPar(i)%logscale) then
c	log
					var(i)=log10(RetPar(i)%value/RetPar(i)%xmin)/log10(RetPar(i)%xmax/RetPar(i)%xmin)
				else if(RetPar(i)%squarescale) then
c	square
					var(i)=(RetPar(i)%value**2-RetPar(i)%xmin**2)/(RetPar(i)%xmax**2-RetPar(i)%xmin**2)
				else
c	linear
					var(i)=(RetPar(i)%value-RetPar(i)%xmin)/(RetPar(i)%xmax-RetPar(i)%xmin)
				endif
				if(.not.var(i).gt.0d0) var(i)=0d0
				if(var(i).gt.1d0) var(i)=1d0
			endif
		enddo
		doscaleR2=.false.
		recomputeopac=.true.
		deallocate(allspec)
		goto 2
	else
		scale=1d0
	endif

	k=0
	do i=1,nobs
		if(ObsSpec(i)%scaling) then
			xy=0d0
			xx=0d0
			do j=1,ObsSpec(i)%ndata
				k=k+1
				ymod(k)=allspec(i,j)
				xy=xy+ymod(k)*ObsSpec(i)%y(j)/dy(k)**2
				xx=xx+ymod(k)*ymod(k)/dy(k)**2
			enddo
			if(ObsSpec(i)%dscale.gt.0d0) then
				xx=xx+1d0/ObsSpec(i)%dscale**2
				xy=xy+1d0/ObsSpec(i)%dscale**2
			endif
			ObsSpec(i)%scale=1d0
			if(xx.gt.0d0) ObsSpec(i)%scale=xx/xy
		else
			ObsSpec(i)%scale=1d0
			do j=1,ObsSpec(i)%ndata
				k=k+1
				ymod(k)=allspec(i,j)
			enddo
		endif
	enddo

	global_chi2=0d0
	tot=0d0
	k=0
	do i=1,nobs
		do j=1,ObsSpec(i)%ndata
			k=k+1
			ymod(k)=allspec(i,j)
			global_chi2=global_chi2+((ymod(k)-ObsSpec(i)%scale*ObsSpec(i)%y(j))/(ObsSpec(i)%scale*dy(k)))**2
			tot=tot-log(sqrt(2d0*pi)*dy(k))
			ObsSpec(i)%model(j)=allspec(i,j)
		enddo
		if(ObsSpec(i)%scaling.and.ObsSpec(i)%dscale.gt.0d0) then
			global_chi2=global_chi2+((ObsSpec(i)%scale-1d0)/ObsSpec(i)%dscale)**2
			tot=tot-log(sqrt(2d0*pi)*ObsSpec(i)%dscale)
		endif
	enddo
	if(planetform.and..not.simAb_converge) then
		global_chi2=global_chi2+((Mplanet-MSimAb)/(Mplanet*1d-3))**2
	endif
	if(massprior) then
		global_chi2=global_chi2+((Mplanet/Mjup-Mp_prior)/dMp_prior)**2
		tot=tot-log(sqrt(2d0*pi)*dMp_prior)
		k=k+1
	endif
	if(modelfail) global_chi2=global_chi2*1d50
	if(IsNaN(global_chi2)) global_chi2=1d50
	lnew=-global_chi2/2d0+tot
	global_chi2=global_chi2/real(max(1,k-n_ret))
	if(free_tprofile.and.wiggle_err.gt.0d0) then
		tot=-(real(nr-2)/2d0)*log(2d0*pi*wiggle_err)
		do i=2,nr-1
			tot=tot-log(T(i+1)*T(i-1)/T(i)**2)**2/(0.5d0*wiggle_err*log(P(i+1)/P(i-1))**2)
		enddo
		lnew=lnew+tot
	endif
	global_like=lnew

	if(writeWolk) write(31,*) imodel,global_chi2,var(1:nvars),COratio,metallicity
	if((.not.useobsgrid.or.dochemistry.or.do3D.or.computeT).and.writeWolk) call flush(31)

	if(lnew.gt.bestlike) then
		inquire(file="improve.sh",exist=truefalse)
		if(truefalse) then
			write(command,'("./improve.sh ",f10.3)') global_chi2
			status=system(command)
		endif

		i2d=0
		call WriteStructure()
		call WriteOutput()

		do i=1,nobs
			select case(ObsSpec(i)%type)
				case("trans","transmission","emisr","emisR","emisa","emis","emission","transC","phase","phaser","phaseR","transM","transE")
					open(unit=20,file=trim(outputdir) // "obs" // trim(int2string(i,'(i0.3)')),FORM="FORMATTED",ACCESS="STREAM")
					do j=1,ObsSpec(i)%ndata
						write(20,*) ObsSpec(i)%lam(j)*1d4,ObsSpec(i)%model(j)/ObsSpec(i)%scale,ObsSpec(i)%scale*ObsSpec(i)%y(j),ObsSpec(i)%dy(j)
					enddo
					close(unit=20)
					open(unit=20,file=trim(outputdir) // "fullobs" // trim(int2string(i,'(i0.3)')),FORM="FORMATTED",ACCESS="STREAM")
					do j=1,nlam
						write(20,*) lam(j)*1d4,specsave(i,j)/ObsSpec(i)%scale
					enddo
					close(unit=20)
				case("lightcurve")
					im=0
					do ilam=1,ObsSpec(i)%nlam
						open(unit=20,file=trim(outputdir) // "lightcurve" // trim(int2string(i,'(i0.4)')) // "_"  
     &								// trim(int2string(ilam,'(i0.4)')) // ".dat",FORM="FORMATTED",ACCESS="STREAM")
						do j=1,ObsSpec(i)%nt
							im=(j-1)*ObsSpec(i)%nlam+ilam
							write(20,*) ObsSpec(i)%t(j),ObsSpec(i)%model(im),ObsSpec(i)%y(im)
						enddo
					enddo
					close(unit=20)
			end select
		enddo

		status=system("cp " // trim(outputdir) // "input.dat " // trim(outputdir) // "bestfit.dat")
		open(unit=21,file=trim(outputdir) // "bestfit.dat",FORM="FORMATTED",access='APPEND')
		write(21,'("*** retrieval keywords ***")')
		write(21,'("retrieval=.false.")')
		do i=1,n_ret
			write(21,'(a," = ",es14.7)') trim(RetPar(i)%keyword),RetPar(i)%value
		enddo
		close(unit=21)	

		bestlike=lnew
		bestvar=var
		bestchi2=global_chi2
	endif
	
	deallocate(allspec)
	
	return
	end
	



	subroutine WritePTlimits_nolimits(var0,Cov,ErrVec,error,chi2,ioflag)
	use GlobalSetup
	use Constants
	use RetrievalMod
	IMPLICIT NONE
	real*8 minT(nr),maxT(nr),var(n_ret),error(2,n_ret),tot,w(n_ret),chi2,ErrVec(n_ret,n_ret)
	real*8 var0(n_ret),random,vec(n_ret),Tbest(nr),gasdev,Cov(n_ret,n_ret)
	real*8 dvar(n_ret),value(n_ret)
	real*8 obsAerr(2,nlam),emiserr(2,nlam),emisRerr(2,nlam),errdummy(2,n_ret)
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
		ErrVec(1:n_ret,i)=ErrVec(1:n_ret,i)*sqrt(chi2*w(i))/sqrt(tot)
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
	errdummy=0d0
	call MapRetrieval(var,errdummy)
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
		open(unit=35,file=trim(outputdir) // "error_cloud.dat",FORM="FORMATTED",ACCESS="STREAM")
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
			vec(i)=gasdev(idum)
		enddo
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

		call MapRetrieval(var,errdummy)
		
		do i=1,n_ret
			value(i)=RetPar(i)%value
			if(RetPar(i)%logscale) value(i)=log10(value(i))
		enddo
		if(ioflag) write(35,form) (value(i),i=1,n_ret),COratio,metallicity

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

		open(unit=45,file=trim(outputdir) // "limits.dat",FORM="FORMATTED",ACCESS="STREAM")
		do i=1,nr
			write(45,*) P(i),Tbest(i),minT(i),maxT(i)
		enddo
		close(unit=45)

		if(speclimits) then
			obsAerr=sqrt(obsAerr/real(iobsA-1))
			emiserr=sqrt(emiserr/real(iemis-1))
			emisRerr=sqrt(emisRerr/real(iemisR-1))

			open(unit=45,file=trim(outputdir) // "emis_limits.dat",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nlam
				write(45,*) lam(i)*1d4,emis0(i),emiserr(1,i),emiserr(2,i)
			enddo
			close(unit=45)

			open(unit=45,file=trim(outputdir) // "emisR_limits.dat",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nlam
				write(45,*) lam(i)*1d4,emisR0(i),emisRerr(1,i),emisRerr(2,i)
			enddo
			close(unit=45)

			open(unit=45,file=trim(outputdir) // "trans_limits.dat",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nlam
				write(45,*) lam(i)*1d4,obsA0(i),obsAerr(1,i),obsAerr(2,i)
			enddo
			close(unit=45)
		endif
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
	real*8 dvar(n_ret),value(n_ret),COret0
	real*8 obsAerr(2,nlam),emiserr(2,nlam),emisRerr(2,nlam),errdummy(2,n_ret)
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
	errdummy=0d0
	call MapRetrieval(var,errdummy)
	call InitDens()
	call SetupStructure(.false.)
	Tbest=T
	maxT=0d0
	minT=0d0
	COerr=0d0
	COret0=COret
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
		open(unit=35,file=trim(outputdir) // "error_cloud.dat",FORM="FORMATTED",ACCESS="STREAM")
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
c			vec(i)=gasdev(idum)
		enddo
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

		call MapRetrieval(var,errdummy)
		
		do i=1,n_ret
			value(i)=RetPar(i)%value
			if(RetPar(i)%logscale) value(i)=log10(value(i))
		enddo
		call InitDens()
		call SetupStructure(.false.)
		if(ioflag) write(35,form) (value(i),i=1,n_ret),COratio,metallicity

		if(k.eq.1) then
			Tbest(1:nr)=T(1:nr)
			COret0=COret
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
		if(COret.lt.COret0) then
			COerr(1)=COerr(1)+(COret-COret0)**2
			iCO(1)=iCO(1)+1
		endif
		if(COret.gt.COret0) then
			COerr(2)=COerr(2)+(COret-COret0)**2
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

		open(unit=45,file=trim(outputdir) // "limits.dat",FORM="FORMATTED",ACCESS="STREAM")
		do i=1,nr
			write(45,*) P(i),Tbest(i),minT(i),maxT(i)
		enddo
		close(unit=45)

		if(speclimits) then
			obsAerr=sqrt(obsAerr/real(iobsA-1))
			emiserr=sqrt(emiserr/real(iemis-1))
			emisRerr=sqrt(emisRerr/real(iemisR-1))

			open(unit=45,file=trim(outputdir) // "emis_limits.dat",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nlam
				write(45,*) lam(i)*1d4,emis0(i),emiserr(1,i),emiserr(2,i)
			enddo
			close(unit=45)

			open(unit=45,file=trim(outputdir) // "emisR_limits.dat",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nlam
				write(45,*) lam(i)*1d4,emisR0(i),emisRerr(1,i),emisRerr(2,i)
			enddo
			close(unit=45)

			open(unit=45,file=trim(outputdir) // "trans_limits.dat",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nlam
				write(45,*) lam(i)*1d4,obsA0(i),obsAerr(1,i),obsAerr(2,i)
			enddo
			close(unit=45)
		endif
	endif

	return
	end




	subroutine WritePTlimitsMN()
	use GlobalSetup
	use Constants
	use RetrievalMod
	IMPLICIT NONE
	real*8 minT(nr),maxT(nr),var(n_ret),error(2,n_ret),tot,chi2
	real*8 var0(n_ret),random,vec(n_ret),Tbest(nr),gasdev,Cov(n_ret,n_ret)
	real*8 dvar(n_ret),value(n_ret),COret0,errdummy(2,n_ret)
	real*8 iCO(2),ivar(2,n_ret),imaxT(nr),iminT(nr),w,wmax,dummy
	integer i,j,k,info,nk,ierr,iter,seed
	character*6000 form
	character*10 side
	logical ioflag,exist
	
	inquire(file=trim(outputdir)//".txt",exist=exist)
	if(.not.exist) return

	if(i2d.eq.0) then
		side=" "
	else
		write(side,'("_",i0.2)') i2d
	endif
	
	ioflag=.true.
	
	open(unit=40,file=trim(outputdir)//".txt",FORM="FORMATTED",ACCESS="STREAM")
	wmax=-1d200
3	read(40,*,end=4) w,dummy,var(1:n_ret)
	if(w.gt.wmax) then
		wmax=w
		var0=var
	endif
	goto 3
4	close(unit=40)
	
	call SetOutputMode(.false.)
	var=var0
	errdummy=0d0
	call MapRetrievalMN(var,errdummy)
	call InitDens()
	call SetupStructure(.false.)
	Tbest=T
	maxT=0d0
	minT=0d0
	COerr=0d0
	COret0=COret
	error=0d0
	ivar=0
	imaxT=0
	iminT=0
	iCO=0

	if(ioflag) then
		open(unit=35,file=trim(outputdir) // "error_cloud.dat",FORM="FORMATTED",ACCESS="STREAM")
		form='("#"' // trim(int2string(n_ret+2,'(i4)')) // 'a19)'
		write(35,form) (trim(int2string(i,'(i4)')) // " " // RetPar(i)%keyword,i=1,n_ret)
     &	,(trim(int2string(n_ret+1,'(i4)')) // " COratio"),(trim(int2string(n_ret+2,'(i4)')) // " metallicity")
		form='(' // int2string(n_ret+3,'(i4)') // 'es19.7)'
	endif

	open(unit=40,file=trim(outputdir)//".txt",FORM="FORMATTED",ACCESS="STREAM")
1	continue
		read(40,*,end=2) w,dummy,var(1:n_ret)
		call MapRetrievalMN(var,errdummy)
		
		do i=1,n_ret
			value(i)=RetPar(i)%value
			if(RetPar(i)%logscale) value(i)=log10(value(i))
		enddo
		call InitDens()
		call SetupStructure(.false.)
		if(ioflag) write(35,form) (value(i),i=1,n_ret),COratio,metallicity,w

		do i=1,nr
			if(T(i).gt.Tbest(i)) then
				maxT(i)=maxT(i)+w*(T(i)-Tbest(i))**2
				imaxT(i)=imaxT(i)+w
			endif
			if(T(i).lt.Tbest(i)) then
				minT(i)=minT(i)+w*(T(i)-Tbest(i))**2
				iminT(i)=iminT(i)+w
			endif
		enddo
		do i=1,n_ret
			if(var(i).lt.var0(i)) then
				error(1,i)=error(1,i)+w*(var(i)-var0(i))**2
				ivar(1,i)=ivar(1,i)+w
			endif
			if(var(i).gt.var0(i)) then
				error(2,i)=error(2,i)+w*(var(i)-var0(i))**2
				ivar(2,i)=ivar(2,i)+w
			endif
		enddo
		if(COret.lt.COret0) then
			COerr(1)=COerr(1)+w*(COret-COret0)**2
			iCO(1)=iCO(1)+w
		endif
		if(COret.gt.COret0) then
			COerr(2)=COerr(2)+w*(COret-COret0)**2
			iCO(2)=iCO(2)+w
		endif
	goto 1
2	continue
	close(unit=40)

	error(1,1:n_ret)=sqrt(error(1,1:n_ret)/real(ivar(1,1:n_ret)))
	error(2,1:n_ret)=sqrt(error(2,1:n_ret)/real(ivar(2,1:n_ret)))
	COerr(1:2)=sqrt(COerr(1:2)/real(iCO(1:2)))
	minT(1:nr)=sqrt(minT(1:nr)/real(iminT(1:nr)))
	maxT(1:nr)=sqrt(maxT(1:nr)/real(imaxT(1:nr)))

	call MapRetrievalMN(var0,error)
	call InitDens()
	call SetupStructure(.false.)

	call SetOutputMode(.true.)
	if(ioflag) then
		close(unit=35)

		open(unit=45,file=trim(outputdir) // "limits" // trim(side) //".dat",FORM="FORMATTED",ACCESS="STREAM")
		do i=1,nr
			write(45,*) P(i),Tbest(i),minT(i),maxT(i)
		enddo
		close(unit=45)
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
		open(unit=35,file=trim(outputdir) // "error_cloud.dat",FORM="FORMATTED",ACCESS="STREAM")
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

		open(unit=45,file=trim(outputdir) // "limits.dat",FORM="FORMATTED",ACCESS="STREAM")
		do i=1,nr
			write(45,*) P(i),Tbest(i),minT(i),maxT(i)
		enddo
		close(unit=45)

		if(speclimits) then
			open(unit=45,file=trim(outputdir) // "emis_limits.dat",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nlam
				write(45,*) lam(i)*1d4,emismin(i),emismax(i)
			enddo
			close(unit=45)

			open(unit=45,file=trim(outputdir) // "emisR_limits.dat",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nlam
				write(45,*) lam(i)*1d4,emisRmin(i),emisRmax(i)
			enddo
			close(unit=45)

			open(unit=45,file=trim(outputdir) // "trans_limits.dat",FORM="FORMATTED",ACCESS="STREAM")
			do i=1,nlam
				write(45,*) lam(i)*1d4,obsAmin(i),obsAmax(i)
			enddo
			close(unit=45)
		endif
	endif

	return
	end



	subroutine RemapObs(i,spec,specsave)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,k,ilam,nuse
	real*8 spec(*),x,specsave(*)
	real*8 lamobs(nlam),lamuse(nlam),specuse(nlam)
	real*8 eta,Tirr,tau,expint,starspec(nlam),starspecregrid(nlam)

	do j=1,nlam
		lamobs(j)=lam(j)
	enddo
	select case(ObsSpec(i)%type)
		case("trans","transmission","transC")
			specsave(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
			if(ObsSpec(i)%filter.ne.' ') then
				call regridfilter(lam,specsave,nlam,spec,ObsSpec(i)%ndata,ObsSpec(i)%f,computelam)
			else if(useobsgrid) then
				spec(1:ObsSpec(i)%ndata)=specsave(ObsSpec(i)%ilam)
			else
				nuse=0
				do j=1,nlam
					if(computelam(j)) then
						nuse=nuse+1
						lamuse(nuse)=lamobs(j)
						specuse(nuse)=specsave(j)
					endif
				enddo
				call regridspecres(lamuse,specuse,nuse,
     &					ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%ndata)
    		endif
    		do j=1,ObsSpec(i)%ndata
    			spec(j)=spec(j)+(ObsSpec(i)%lam(j)-ObsSpec(i)%lam(1))*1d4*ObsSpec(i)%slope
    		enddo
     		ObsSpec(i)%model(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)
		case("transM")
			if(do3D) then
				specsave(1:nlam)=2d0*obsA_split(1:nlam,1)/(pi*Rstar**2)
			else
				specsave(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
			endif
			if(ObsSpec(i)%filter.ne.' ') then
				call regridfilter(lam,specsave,nlam,spec,ObsSpec(i)%ndata,ObsSpec(i)%f,computelam)
			else if(useobsgrid) then
				spec(1:ObsSpec(i)%ndata)=specsave(ObsSpec(i)%ilam)
			else
				do j=1,nlam
					if(computelam(j)) then
						nuse=nuse+1
						lamuse(nuse)=lamobs(j)
						specuse(nuse)=specsave(j)
					endif
				enddo
				call regridspecres(lamuse,specuse,nuse,
     &					ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%ndata)
    		endif
     		ObsSpec(i)%model(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)
		case("transE")
			if(do3D) then
				specsave(1:nlam)=2d0*obsA_split(1:nlam,2)/(pi*Rstar**2)
			else
				specsave(1:nlam)=obsA(0,1:nlam)/(pi*Rstar**2)
			endif
			if(ObsSpec(i)%filter.ne.' ') then
				call regridfilter(lam,specsave,nlam,spec,ObsSpec(i)%ndata,ObsSpec(i)%f,computelam)
			else if(useobsgrid) then
				spec(1:ObsSpec(i)%ndata)=specsave(ObsSpec(i)%ilam)
			else
				do j=1,nlam
					if(computelam(j)) then
						nuse=nuse+1
						lamuse(nuse)=lamobs(j)
						specuse(nuse)=specsave(j)
					endif
				enddo
				call regridspecres(lamuse,specuse,nuse,
     &					ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%ndata)
    		endif
     		ObsSpec(i)%model(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)
		case("emisr","emisR")
			specsave(1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))/(Fstar(1:nlam)*1d23/distance**2)
			if(ObsSpec(i)%filter.ne.' ') then
				specuse(1:nlam)=(phase(1,0,1:nlam)+flux(0,1:nlam))
				starspec(1:nlam)=(Fstar(1:nlam)*1d23/distance**2)
				call regridfilter(lam,specuse,nlam,spec,ObsSpec(i)%ndata,ObsSpec(i)%f,computelam)
				call regridfilter(lam,starspec,nlam,starspecregrid,ObsSpec(i)%ndata,ObsSpec(i)%f,computelam)
				spec(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)/starspecregrid(1:ObsSpec(i)%ndata)
			else if(useobsgrid) then
				spec(1:ObsSpec(i)%ndata)=specsave(ObsSpec(i)%ilam)
			else
				do j=1,nlam
					if(computelam(j)) then
						nuse=nuse+1
						lamuse(nuse)=lamobs(j)
						specuse(nuse)=(phase(1,0,j)+flux(0,j))
						starspec(nuse)=(Fstar(j)*1d23/distance**2)
					endif
				enddo
				call regridspecres(lamuse,specuse,nuse,
     &					ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%ndata)
				call regridspecres(lamuse,starspec,nuse,
     &					ObsSpec(i)%lam,starspecregrid,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%ndata)
				spec(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)/starspecregrid(1:ObsSpec(i)%ndata)
    		endif
     		ObsSpec(i)%model(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)
		case("phaser","phaseR")
			specsave(1:nlam)=(phase(ObsSpec(i)%iphase,0,1:nlam)+flux(0,1:nlam))/(Fstar(1:nlam)*1d23/distance**2)
			if(ObsSpec(i)%filter.ne.' ') then
				specuse(1:nlam)=(phase(ObsSpec(i)%iphase,0,1:nlam)+flux(0,1:nlam))
				starspec(1:nlam)=(Fstar(1:nlam)*1d23/distance**2)
				call regridfilter(lam,specuse,nlam,spec,ObsSpec(i)%ndata,ObsSpec(i)%f,computelam)
				call regridfilter(lam,starspec,nlam,starspecregrid,ObsSpec(i)%ndata,ObsSpec(i)%f,computelam)
				spec(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)/starspecregrid(1:ObsSpec(i)%ndata)
			else if(useobsgrid) then
				spec(1:ObsSpec(i)%ndata)=specsave(ObsSpec(i)%ilam)
			else
				do j=1,nlam
					if(computelam(j)) then
						nuse=nuse+1
						lamuse(nuse)=lamobs(j)
						specuse(nuse)=(phase(ObsSpec(i)%iphase,0,j)+flux(0,j))
						starspec(nuse)=(Fstar(j)*1d23/distance**2)
					endif
				enddo
				call regridspecres(lamuse,specuse,nuse,
     &					ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%ndata)
				call regridspecres(lamuse,starspec,nuse,
     &					ObsSpec(i)%lam,starspecregrid,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%ndata)
				spec(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)/starspecregrid(1:ObsSpec(i)%ndata)
    		endif
     		ObsSpec(i)%model(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)
		case("phase")
			specsave(1:nlam)=phase(ObsSpec(i)%iphase,0,1:nlam)+flux(0,1:nlam)
			if(ObsSpec(i)%filter.ne.' ') then
				call regridfilter(lam,specsave,nlam,spec,ObsSpec(i)%ndata,ObsSpec(i)%f,computelam)
			else if(useobsgrid) then
				spec(1:ObsSpec(i)%ndata)=specsave(ObsSpec(i)%ilam)
			else
				do j=1,nlam
					if(computelam(j)) then
						nuse=nuse+1
						lamuse(nuse)=lamobs(j)
						specuse(nuse)=specsave(j)
					endif
				enddo
				call regridspecres(lamuse,specuse,nuse,
     &					ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%ndata)
    		endif
     		ObsSpec(i)%model(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)
		case("emisa","emis","emission")
			specsave(1:nlam)=phase(1,0,1:nlam)+flux(0,1:nlam)
			if(ObsSpec(i)%filter.ne.' ') then
				call regridfilter(lam,specsave,nlam,spec,ObsSpec(i)%ndata,ObsSpec(i)%f,computelam)
			else if(useobsgrid) then
				spec(1:ObsSpec(i)%ndata)=specsave(ObsSpec(i)%ilam)
			else
				do j=1,nlam
					if(computelam(j)) then
						nuse=nuse+1
						lamuse(nuse)=lamobs(j)
						specuse(nuse)=specsave(j)
					endif
				enddo
				call regridspecres(lamuse,specuse,nuse,
     &					ObsSpec(i)%lam,spec,ObsSpec(i)%R,ObsSpec(i)%Rexp,ObsSpec(i)%ndata)
     		endif
     		ObsSpec(i)%model(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)
			if(log_emis) spec(1:ObsSpec(i)%ndata)=log(spec(1:ObsSpec(i)%ndata))
		case("tprofile")
			call ComputeParamT(spec(1:nr))
			spec(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)*ObsSpec(i)%beta
		case("logtp")
			spec(1)=0d0!(log10(abs(T(1)/T(2))/log10(P(1)/P(2))))
			do j=2,nr-1
				spec(j)=(log10(abs(T(j-1)/T(j)))/log10(P(j-1)/P(j)))-(log10(abs(T(j)/T(j+1)))/log10(P(j)/P(j+1)))
			enddo
			spec(nr)=0d0!(log10(abs(T(nr-1)/T(nr)))/log10(P(nr-1)/P(nr)))
			spec(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)/real(nr)
			spec(1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)*ObsSpec(i)%beta
		case("prior","priors")
			do j=1,n_ret
				if(RetPar(j)%logscale) then
					spec(j)=log10(RetPar(j)%value)
				else
					spec(j)=RetPar(j)%value
				endif
			enddo
		case("lightcurve")
			do j=1,ObsSpec(i)%ndata
				spec(j)=ObsSpec(i)%model(j)
			enddo
	end select

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
	real*8 var(n_ret),dvar(2,n_ret),x,xx,erfinv

	do k=1,2

	do i=1,n_ret
		if(i.gt.1.and.RetPar(i)%increase) then
			RetPar(i)%xmin=RetPar(i-1)%value
		endif
		if(RetPar(i)%logscale) then
c	log
			x=var(i)
			if(RetPar(i)%keyword(1:6).eq."Ppoint") then
				read(RetPar(i)%keyword(7:len(RetPar(i)%keyword)),*) j
				x=1d0-x**(1d0/(real(nTpoints-j+1)))
			endif
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
	do i=1,nclouds
		Cloud(i)%rnuc=Cloud(i)%rnuc/micron
	enddo
	orbit_inc=orbit_inc*180d0/pi

	metallicity=metallicity0
	if(WaterWorld) then
		Pmax=WWInit_Pmax
		Pmin=WWInit_Pmin
		Pplanet=WWInit_Pplanet
		mixrat(1:nmol)=WWInit_mixrat(1:nmol)
	endif
	Kzz_convect=0d0
	do i=1,n_ret
		readline=trim(RetPar(i)%keyword) // "=" // trim(dbl2string(RetPar(i)%value,'(es14.7)'))
		call get_key_value(readline,key%key,key%key1,key%key2,key%orkey1,key%orkey2,key%value,
     &					key%nr1,key%nr2,key%hasnr1,key%hasnr2,key%key2d)
		call ReadAndSetKey(key)
	enddo
	call ConvertUnits()
	metallicity0=metallicity
	enddo
	call RefreshMaterialCloud()

	return
	end
	

	
	subroutine MapRetrievalMN(var,dvar)
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
		if(RetPar(i)%logscale) then
c	log
			x=var(i)
			RetPar(i)%value=10d0**x
		else
c	linear, square
			x=var(i)
			RetPar(i)%value=x
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
	do i=1,nclouds
		Cloud(i)%rnuc=Cloud(i)%rnuc/micron
	enddo
	orbit_inc=orbit_inc*180d0/pi

	if(WaterWorld) then
		Pmax=WWInit_Pmax
		Pmin=WWInit_Pmin
		Pplanet=WWInit_Pplanet
		mixrat(1:nmol)=WWInit_mixrat(1:nmol)
	endif
	Kzz_convect=0d0
	do i=1,n_ret
		readline=trim(RetPar(i)%keyword) // "=" // trim(dbl2string(RetPar(i)%value,'(es14.7)'))
		call get_key_value(readline,key%key,key%key1,key%key2,key%orkey1,key%orkey2,key%value,
     &					key%nr1,key%nr2,key%hasnr1,key%hasnr2,key%key2d)
		call ReadAndSetKey(key)
	enddo
	call ConvertUnits()
	metallicity0=metallicity
	enddo
	call RefreshMaterialCloud()

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
		if(i.gt.1.and.RetPar(i)%increase) then
			RetPar(i)%xmin=RetPar(i-1)%value
		endif
		if(RetPar(i)%logscale) then
c	log
			var(i)=log10(RetPar(i)%value/RetPar(i)%xmin)/log10(RetPar(i)%xmax/RetPar(i)%xmin)
			if(RetPar(i)%keyword(1:6).eq."Ppoint") then
				read(RetPar(i)%keyword(7:len(RetPar(i)%keyword)),*) j
				var(i)=(1d0-var(i))**(real(nTpoints-j+1))
			endif
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
					y1(i1)=y0(i0)+(y0(i0+1)-y0(i0))*(x1(i1)-x0(i0))/(x0(i0+1)-x0(i0))
				endif
			enddo
		endif
	enddo

	i1=n1-1
	do while(IsNaN(y1(n1)).and.i1.gt.0) 
		y1(n1)=y1(i1)
		i1=i1-1
	enddo

1	continue
	i0=0
	do i1=1,n1-1
		if(IsNaN(y1(i1))) then
			y1(i1)=y1(i1+1)
			i0=i0+1
		endif
	enddo
	if(i0.eq.n1-1) then
		y1=1d0
	else if(i0.gt.0) then
		goto 1
	endif

	return
	end
	

	subroutine regridspecres(x0,y0,n0,x1,y1,R1,expR1,n1)
	IMPLICIT NONE
	integer i1,i0,n0,n1,i,nn0
	real*8 x0(n0),y0(n0),x1(n1),y1(n1),R1(n1),expR1(n1),w,tot
	real*8 xx0(n0+n1),yy0(n0+n1)

	xx0(1:n0)=x0(1:n0)
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
				if(abs((xx0(i0)-x1(i1))*2d0*R1(i1)/(x1(i1)))**expR1(i1).lt.50d0) then
				w=exp(-abs((xx0(i0)-x1(i1))*2d0*R1(i1)/(x1(i1)))**expR1(i1))
				if(i0.eq.1) then
					w=w*abs(1d0/xx0(2)-1d0/xx0(1))/2d0
				else if(i0.eq.nn0) then
					w=w*abs(1d0/xx0(nn0)-1d0/xx0(nn0-1))/2d0
				else
					w=w*abs(1d0/xx0(i0-1)-1d0/xx0(i0+1))/2d0
				endif
				tot=tot+w
				y1(i1)=y1(i1)+w*yy0(i0)
				endif
			enddo
			y1(i1)=y1(i1)/tot
		endif
	enddo

	return
	end
	

	subroutine regridfilter(x0,y0,n0,y1,n1,filter,computelam)
	IMPLICIT NONE
	integer n0,n1,i0,i1
	real*8 x0(n0),y0(n0),y1(n1),filter(n0,n1),w,tot
	logical computelam(n0)

	do i1=1,n1
		y1(i1)=0d0
		tot=0d0
		do i0=1,n0
			if(computelam(i0)) then
				w=filter(i0,i1)
				tot=tot+w
				y1(i1)=y1(i1)+w*y0(i0)
			endif
		enddo
		y1(i1)=y1(i1)/tot
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

98	open(unit=20,file=trim(outputdir) // "retrieval",FORM="FORMATTED",ACCESS="STREAM",ERR=99)
	goto 100
99	write(command,'("mkdir -p ",a)') trim(outputdir)
	call system(command)
	call sleep(10)
	goto 98
100	continue
	write(20,'("Model ",i7)') imodel
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
			write(20,'(a15," = ",es14.7," +/- ",es11.4,es11.4)') 'f_dry',f_dry,0d0,0d0
			write(20,'(a15," = ",es14.7," +/- ",es11.4,es11.4)') 'f_wet',f_wet,0d0,0d0
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
	
c	do i=1,nobs
c		select case(ObsSpec(i)%type)
c			case("trans","transmission","emisr","emisR","emisa","emis","emission","transC")
c				open(unit=20,file=trim(outputdir) // "obs" // trim(int2string(i,'(i0.3)')),FORM="FORMATTED",ACCESS="STREAM")
c				do j=1,ObsSpec(i)%ndata
c					write(20,*) ObsSpec(i)%lam(j)*1d4,ObsSpec(i)%model(j),ObsSpec(i)%y(j),ObsSpec(i)%dy(j)
c				enddo
c				close(unit=20)
c		end select
c	enddo

	call MapRetrieval(bestvar,error)

	call system("cp " // trim(outputdir) // "input.dat " // trim(outputdir) // "bestfit.dat")
	open(unit=21,file=trim(outputdir) // "bestfit.dat",FORM="FORMATTED",access='APPEND')
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

	

	subroutine ComputeLike(lnew)
	use GlobalSetup
	use Constants
	use RetrievalMod
	use Struct3D
	IMPLICIT NONE
	integer nvars,i,j,nlamtot,ny,k,maxspec,im,ilam,status,system,ii
	real*8 ymod(nlam*nobs),lnew,spectemp(nlam*nobs),specsave(nobs,nlam*nobs)
	real*8,allocatable :: spec(:),allspec(:,:)
	logical recomputeopac,truefalse,doscaleR2
	real*8 xx,xy,scale,dy(nlam*nobs),tot
	character*100 command
	
	maxspec=0
	do i=1,nobs
		if(ObsSpec(i)%ndata.gt.maxspec) maxspec=ObsSpec(i)%ndata
	enddo
	allocate(allspec(nobs,maxspec))
	do i=1,nobs
		allocate(spec(ObsSpec(i)%ndata))
		call RemapObs(i,spec,spectemp)
		if(ObsSpec(i)%i2d.eq.i2d) then
			allspec(i,1:ObsSpec(i)%ndata)=spec(1:ObsSpec(i)%ndata)
			specsave(i,1:nlam)=spectemp(1:nlam)
		endif
		deallocate(spec)
	enddo

	k=0
	do i=1,nobs
		do j=1,ObsSpec(i)%ndata
			k=k+1
			dy(k)=sqrt(ObsSpec(i)%dy(j)**2+ObsSpec(i)%adderr**2)
			do ii=1,nmodel_err-1
				if(ObsSpec(i)%lam(j).lt.model_err_lam(ii)) exit
			enddo
			select case(ObsSpec(i)%type)
				case("emisa","emis","emission","phase")
					dy(k)=sqrt(dy(k)**2+model_err_abs(ii)**2)
				case("trans","transmission","emisr","emisR","transC","phaser","phaseR","transM","transE")
					dy(k)=sqrt(dy(k)**2+model_err_rel(ii)**2)
			end select
		enddo
	enddo

	k=0
	do i=1,nobs
		if(ObsSpec(i)%scaling) then
			xy=0d0
			xx=0d0
			do j=1,ObsSpec(i)%ndata
				k=k+1
				ymod(k)=allspec(i,j)
				xy=xy+ymod(k)*ObsSpec(i)%y(j)/dy(k)**2
				xx=xx+ymod(k)*ymod(k)/dy(k)**2
			enddo
			if(ObsSpec(i)%dscale.gt.0d0) then
				xx=xx+1d0/ObsSpec(i)%dscale**2
				xy=xy+1d0/ObsSpec(i)%dscale**2
			endif
			ObsSpec(i)%scale=1d0
			if(xx.gt.0d0) ObsSpec(i)%scale=xx/xy
		else
			ObsSpec(i)%scale=1d0
			do j=1,ObsSpec(i)%ndata
				k=k+1
				ymod(k)=allspec(i,j)
			enddo
		endif
	enddo

	global_chi2=0d0
	tot=0d0
	k=0
	do i=1,nobs
		do j=1,ObsSpec(i)%ndata
			k=k+1
			ymod(k)=allspec(i,j)
			global_chi2=global_chi2+((ymod(k)-ObsSpec(i)%scale*ObsSpec(i)%y(j))/(ObsSpec(i)%scale*dy(k)))**2
			tot=tot-log(sqrt(2d0*pi)*dy(k))
			ObsSpec(i)%model(j)=allspec(i,j)
		enddo
		if(ObsSpec(i)%scaling.and.ObsSpec(i)%dscale.gt.0d0) then
			global_chi2=global_chi2+((ObsSpec(i)%scale-1d0)/ObsSpec(i)%dscale)**2
			tot=tot-log(sqrt(2d0*pi)*ObsSpec(i)%dscale)
		endif
	enddo
	if(planetform.and..not.simAb_converge) then
		global_chi2=global_chi2+((Mplanet-MSimAb)/(Mplanet*1d-3))**2
	endif
	if(massprior) then
		global_chi2=global_chi2+((Mplanet/Mjup-Mp_prior)/dMp_prior)**2
		tot=tot-log(sqrt(2d0*pi)*dMp_prior)
		k=k+1
	endif
	lnew=-global_chi2/2d0+tot
	global_chi2=global_chi2/real(max(1,k-n_ret))
	if(free_tprofile.and.wiggle_err.gt.0d0) then
		tot=-(real(nr-2)/2d0)*log(2d0*pi*wiggle_err)
		do i=2,nr-1
			tot=tot-log(T(i+1)*T(i-1)/T(i)**2)**2/(0.5d0*wiggle_err*log(P(i+1)/P(i-1))**2)
		enddo
		lnew=lnew+tot
	endif
	global_like=lnew
	
	deallocate(allspec)
	
	return
	end
	


	subroutine MCMC(func,var0,nvar,nburn,nstep,ny)
	IMPLICIT NONE
	integer nvar,nburn,nstep,j,i,idum,iv,ny,ii,idstep
	real*8 var0(nvar),var(nvar,nstep),func,stepsize,varp(nvar),random
	real*8 l1,l2,f,r,dx(nvar),ll(nstep),accepted(nstep+nburn),accr,fs
	real*8 dstep(nvar),dav(nvar),gasdev
	external func

	idum=-42
	stepsize=0.1d0
	dav=0.1d0
	l1=1d-20
	varp=var0
	j=0
	fs=0.9
	idstep=0
	do i=1,nburn+nstep
		l2=func(varp,ny)
		f=exp(l2-l1)
		r=random(idum)
		if(f.gt.r.or.i.eq.1) then
			accepted(i)=1.0
c			if(i.gt.nburn) then
				j=j+1
				var(1:nvar,j)=varp(1:nvar)
				ll(j)=l2
c			endif
			var0(1:nvar)=varp(1:nvar)
			l1=l2
			if(i.gt.nburn) then
				do iv=1,nvar
					dav(iv)=(dav(iv)*real(idstep)+min(0.25,dstep(iv)))/real(idstep+1)
				enddo
				idstep=idstep+1
			endif
		else
			accepted(i)=0.0
		endif
		if(i.eq.nburn) j=0
		do iv=1,nvar
			dstep(iv)=gasdev(idum)*dav(iv)*stepsize
			varp(iv)=var0(iv)+dstep(iv)
			if(varp(iv).gt.1d0) varp(iv)=2d0-varp(iv)
			if(varp(iv).lt.0d0) varp(iv)=-varp(iv)
		enddo
	enddo

	print*,j
	open(unit=20,file='output.dat',RECL=6000)
	do i=1,j
		write(20,*) ll(i),var(1:nvar,i)
	enddo
	close(unit=20)
	
	return
	end

	subroutine randomdirectionN(dx,N,idum)
	IMPLICIT NONE
	integer N,i,idum
	real*8 dx(N),random,r

1	continue
	r=0d0
	do i=1,N
		dx(i)=2d0*(random(idum)-0.5d0)
		r=r+dx(i)**2
		if(r.gt.1d0) goto 1
	enddo
	dx=dx/sqrt(r)
	
	return
	end
	

	real*8 function MCMCfunc(var,ny)
	use GlobalSetup
	use Constants
	use RetrievalMod
	IMPLICIT NONE
	integer ny,i
	real*8 ymod(ny),var(n_ret),lnew,scale

	call mrqcomputeY(var,ymod,n_ret,ny,lnew,scale)
	MCMCfunc=lnew
	
	return
	end


	
	
	
	


	