	subroutine Run3D(recomputeopac)
	use GlobalSetup
	use Constants
	use Struct3D
	IMPLICIT NONE
	integer i,j,icloud,k,irtrace,iptrace,inu,imol
	real*8 beta(nlong,nlatt),Planck,phi,la,lo,A,rr
	real*8 b1,b2,betamin,betamax,freq0,Rmax,theta
	real*8,allocatable :: Ca(:,:,:,:,:),Cs(:,:,:),BBr(:,:),Si(:,:,:,:,:),Ca_mol(:,:,:,:,:),Ce_cont(:,:,:)
	integer ir,ilam,ig,isize,iRmax,ndisk,nsub,nrtrace,nptrace,ipc,npc,nmol_count
	logical recomputeopac
	logical,allocatable :: hit(:,:)
	real*8,allocatable :: rtrace(:),wrtrace(:),ftot(:),rphi_image(:,:,:),xy_image(:,:,:),cloud3D(:,:)
	real*8 x,y,z,vx,vy,vz,v,w1,w2
	integer edgeNR,i1,i2,i3,i1next,i2next,i3next,edgenext
	real*8,allocatable :: fluxp(:),tau(:,:),fact(:,:),tautot(:,:),exp_tau(:,:),obsA_split_omp(:,:),dtauR_nu(:,:,:,:,:)
	real*8,allocatable :: tauc(:),Afact(:),vv(:,:),obsA_omp(:),mixrat3D(:,:,:),T3D(:,:),fluxp_omp(:)
	real*8 g,tot,contr,tmp(nmol),Rmin_im,Rmax_im,random,xmin,xmax,Pb(nr+1)
	integer nx_im,ix,iy,ni,ilatt,ilong,imustar,ivel,miniter0
	character*500 file
	real*8 tau1,fact1,exp_tau1,maximage,beta_c,NormSig,Fstar_temp(nlam),SiR1,SiRalb1,tau0
	real*8,allocatable :: maxdet(:,:),SiSc(:,:,:,:,:),alb_omp(:),SiR0(:,:),SiRalb0(:,:),R3DC(:,:),lgrid(:)
	logical iterateshift,actually1D,do_ibeta(n3D)
	real*8 vxxmin,vxxmax,ComputeKzz,betamin_term,tot1,tot2,vrot,dlam_rot,albedo_day,scale,scale_prev
	logical docloud0(max(nclouds,1)),do_rot
	
	allocate(Ca(nlam,ng,nr,n3D,-nvel:nvel),Cs(nlam,nr,n3D),BBr(nlam,0:nr),Si(nlam,ng,0:nr,nnu0,n3D))
	if(computealbedo) allocate(SiSc(nlam,ng,0:nr,nnu0,n3D))
	allocate(Ca_mol(nlam,ng,nmol,nr,n3D),Ce_cont(nlam,nr,n3D))
	allocate(dtauR_nu(nlam,ng,n3D,nr,-nvel:nvel))
	allocate(R3D(n3D,nr+2))
	allocate(R3D2(n3D,nr+2))
	allocate(R3DC(n3D,nr+2))
	allocate(T3D(n3D,0:nr))
	allocate(local_albedo(n3D))
	if(.not.retrieval.and.fulloutput3D) allocate(cloud3D(n3D,0:nr))
	allocate(mixrat3D(n3D,nr,nmol))
	nx_im=200

	if(retrieval.or.dopostequalweights) call SetOutputMode(.false.)

c	recomputeopac=.true.
	cloudfrac=1d0
	docloud0=.false.
	do i=1,nclouds
		if(Cloud(i)%coverage.gt.0.5) docloud0(i)=.true.
		if(Cloud(i)%coverage.ne.1d0) then
			call output("WARNING! Cloud coverage fraction does not work in 3D structures!")
		endif
	enddo

	if(fixnight2day.or.(WaterWorld.and..not.setsurfpressure).or.(deepredist.and.deepredisttype.eq.'fixflux')) then
		do3D=.false.
		init3D=.true.
		if(.not.fixnight2day) betaT=2d0/3d0-(5d0/12d0)*night2day
		call SetOutputMode(.false.)
		call InitDens()
		call ComputeModel1D(recomputeopac)
		if(.not.retrieval.and..not.dopostequalweights) then
			call SetOutputMode(.true.)
			call output("night2day contrast: " // dbl2string(night2day,'(f7.4)'))
		endif
		init3D=.false.
		do3D=.true.
	endif

	iterateshift=.false.
	if(hotspotshift0.ge.-180d0.and.hotspotshift0.le.180d0.and.tidallock) then
		iterateshift=.true.
		vxxmin=-4d0
		vxxmax=10d0
		vxx=10d0**((vxxmin+vxxmax)/2d0)
		if(hotspotshift0.lt.0d0) vxx=-vxx
		if(hotspotshift0.eq.0d0) then
			vxx=0d0
			iterateshift=.false.
		endif
	endif

	i=0
	do while((iterateshift.and.
     &		abs(hotspotshift-hotspotshift0).gt.5d-3.and.
     &		abs((vxxmax-vxxmin)/(vxxmax+vxxmin)).gt.1d-4.and.
     &		i.lt.999).
     &		or.i.eq.0)

	call Setup3D(beta,long,latt,nlong,nlatt,Kxx,Kyy,vxx,powvxx,night2day,pole2eq,fDay,betamin,betamax,tidallock)

	if(n3D.eq.2) then
		do ilong=1,nlong-1
			do ilatt=1,nlatt-1
				lo=-pi+2d0*pi*(real(ilong)-0.5d0)/real(nlong-1)
				la=-pi/2d0+pi*(real(ilatt)-0.5d0)/real(nlatt-1)
				if(lo.le.0d0) then
					beta(ilong,ilatt)=0d0
				else
					beta(ilong,ilatt)=1d0
				endif
			enddo
		enddo
		betamin=0d0
		betamax=1d0
	endif

	if(vxx.eq.0d0.or..not.tidallock) then
		hotspotshift=0d0
	else
		call DetermineShift(long,beta,nlong,nlatt,hotspotshift)
		hotspotshift=hotspotshift-180d0
	endif

	if(iterateshift) then
		if(abs(hotspotshift).gt.abs(hotspotshift0)) then
			vxxmax=log10(abs(vxx))
			vxx=10d0**((vxxmin+vxxmax)/2d0)
		else
			vxxmin=log10(abs(vxx))
			if(i.eq.0) then
				vxx=10d0**vxxmax
			else
				vxx=10d0**((vxxmin+vxxmax)/2d0)
			endif
		endif
		if(hotspotshift0.lt.0d0) vxx=-vxx
	endif
	i=i+1

	enddo

	actually1D=.true.
	if(((vxx.ne.0d0.or.night2day.ne.1d0.or.pole2eq.ne.1d0).and.betamax.ne.betamin).or.(n3D.eq.2.and.n_Par3D.ne.0)) 
     &		actually1D=.false.

	call output("hotspot shift: " // dbl2string(hotspotshift,'(f6.2)') // " degrees")

	if(iterateshift.and.abs(hotspotshift-hotspotshift0).gt.5d-3) call output("Desired hotspot shift could not be obtained!!")

	if(.not.retrieval.and..not.dopostequalweights) then
		open(unit=20,file=trim(outputdir) // "structure3D.dat",FORM="FORMATTED",ACCESS="STREAM")
		do j=1,nlatt-1
			write(20,*) beta(1:nlong-1,j)
		enddo
		close(unit=20)
	endif

	do i=1,nlong
		if(long(i).lt.(pi/4d0)) then
			tanx(i)=sin(long(i))/cos(long(i))
			tany(i)=1d0
		else if(long(i).lt.(3d0*pi/4d0)) then
			tanx(i)=-1d0
			tany(i)=-cos(long(i))/sin(long(i))
		else if(long(i).lt.(5d0*pi/4d0)) then
			tanx(i)=sin(long(i))/cos(long(i))
			tany(i)=1d0
		else if(long(i).lt.(7d0*pi/4d0)) then
			tanx(i)=-1d0
			tany(i)=-cos(long(i))/sin(long(i))
		else
			tanx(i)=sin(long(i))/cos(long(i))
			tany(i)=1d0
		endif
	enddo
	do i=1,nlatt
		cost2(i)=cos(latt(i))**2
	enddo

	ibeta=1
	do_ibeta=.false.
	do i=1,nlong-1
		do j=1,nlatt-1
			if(readFull3D) then
				ibeta(i,j)=(i-1)*(nlatt-1)+j
			else
				if(betamax.eq.betamin.or.(night2day.eq.1d0.and.vxx.eq.0d0.and.pole2eq.eq.1d0.and.n3D.ne.2)) then
					ibeta(i,j)=1
				else
					ibeta(i,j)=(real(n3D)*((beta(i,j)-betamin)/(betamax-betamin)))+1
					if(ibeta(i,j).gt.n3D) ibeta(i,j)=n3D
					if(.not.ibeta(i,j).gt.1) ibeta(i,j)=1
				endif
			endif
			do_ibeta(ibeta(i,j))=.true.
		enddo
	enddo
	Rmax=0d0
	iRmax=1


	beta_c=0d0
	betamin_term=1d0
	i=max(1,nlong/4)
	do j=1,nlatt-1
		beta_c=beta_c+beta(i,j)
		if(beta(i,j).lt.betamin_term) betamin_term=beta(i,j)
	enddo
	i=max(1,3*nlong/4)
	do j=1,nlatt-1
		beta_c=beta_c+beta(i,j)
		if(beta(i,j).lt.betamin_term) betamin_term=beta(i,j)
	enddo
	beta_c=beta_c/real(2*(nlatt-1))	

	if(.not.retrieval.and..not.domakeai.and..not.dopostequalweights) then
		open(unit=20,file=trim(outputdir) // "parameter3D.dat",FORM="FORMATTED",ACCESS="STREAM")
		do j=1,nlatt-1
			do i=1,nlong-1
				beta3D_eq(i)=beta(i,j)
				x3D_eq(i)=NormSig(beta(i,j),par3Dsteepness,beta_c,betamin,betamax)
			enddo
			write(20,*) x3D_eq(1:nlong-1)
		enddo
		close(unit=20)
	endif

	if(fulloutput3D) then
		j=nlatt/2
		do i=1,nlong-1
			beta3D_eq(i)=beta(i,j)
			x3D_eq(i)=NormSig(beta3D_eq(i),par3Dsteepness,beta_c,betamin,betamax)
		enddo
	endif

	call output("Computing multiple 1D structures")

	albedo_day=0.5d0
	local_albedo=0.5d0
	scale=1d0

	miniter0=miniter
	do i_alb=1,nalbedo_iter

	call tellertje_perc(0,n3D)
	do i=n3D,1,-1
		call SetOutputMode(.false.)
		i3D=i
		beta3D(i)=betamin+(betamax-betamin)*(real(i3D)-0.5)/real(n3D)
		if(n3D.eq.2) then
			if(beta3D(i).lt.0.5d0) then
				beta3D(i)=0d0
			else
				beta3D(i)=1d0
			endif
			x3D(i)=beta3D(i)
		else
			x3D(i)=NormSig(beta3D(i),par3Dsteepness,beta_c,betamin,betamax)
		endif
		do j=1,n_Par3D
			xmin=Par3D(j)%xmin
			if(Par3D(j)%multiply) then
				xmax=xmin*Par3D(j)%xmax
			else
				xmax=Par3D(j)%xmax
			endif
			if(Par3D(j)%logscale) then
				Par3D(j)%x=10d0**(log10(xmin)+log10(xmax/xmin)*x3D(i))
			else
				Par3D(j)%x=xmin+(xmax-xmin)*x3D(i)
			endif
		enddo

		call MapPar3D()
		if(n3D.eq.2) beta3D(i)=betaT

		betaF=beta3D(i)*scale
		if(.not.deepredist.or.deepredisttype.ne.'fixbeta') then
			betaT=beta3D(i)
			tot=0d0
			betaT=0d0
			do ilong=1,nlong-1
				do ilatt=1,nlatt-1
					if(ibeta(ilong,ilatt).eq.i) then
						la=cos((latt(ilatt)+latt(ilatt+1))/2d0-pi/2d0)
						lo=(long(ilong)+long(ilong+1))/2d0-pi
						tot=tot+abs(la)
						if(abs(lo).le.pi/2d0) betaT=betaT+abs(la*la*cos(lo))
					endif
				enddo
			enddo
			betaT=betaT/tot
		endif
		if(.not.betaT.ge.0d0.or..not.betaT.le.1d0) betaT=betaF

c Now call the setup for the readFull3D part
		if(readFull3D) then
			ilong=i/(nlatt-1)+1
			ilatt=i-(ilong-1)*(nlatt-1)
			call DoReadFull3D(i,ilong,ilatt)
		endif

		if((.not.actually1D.and.do_ibeta(i)).or.i.eq.n3D) then
			call InitDens()
			if(i_alb.ne.1) then
				T(1:nr)=T3D(i,1:nr)
				Tsurface=T3D(i,0)
				miniter=min(miniter,2)
			else
				miniter=miniter0
			endif
			call ComputeModel1D(recomputeopac)
			if(i3D.eq.n3D.and.deepredist.and.deepredisttype.eq.'fixflux') then
				i3D=1
				call ComputeModel1D(recomputeopac)
				i3D=i
			endif

			if(R(nr+1).gt.Rmax) then
				Rmax=R(nr+1)
				iRmax=i
			endif
			do ir=1,nr+1
				R3D(i,ir)=R(ir)
			enddo
			do ir=1,nr
				T3D(i,ir)=T(ir)
				mixrat3D(i,ir,1:nmol)=mixrat_r(ir,1:nmol)
				if(.not.retrieval.and.fulloutput3D) cloud3D(i,ir)=sum(cloud_dens(ir,1:nclouds))
				do ilam=1,nlam
					freq0=freq(ilam)
					BBr(ilam,ir)=Planck(T(ir),freq0)
					do ig=1,ng
						do ivel=-nvel,nvel
							call Crossections(ir,ilam,ig,Ca(ilam,ig,ir,i,ivel),Cs(ilam,ir,i),docloud0,ivel)
						enddo
					enddo
				enddo
			enddo
			if(computeT) then
				T3D(i,0)=Tsurface
			else
				T3D(i,0)=T3D(i,1)
			endif
			do ilam=1,nlam
				freq0=freq(ilam)
				BBr(ilam,0)=Planck(T3D(i,0),freq0)
			enddo
			do ir=1,nr
				k=0
				do ilam=1,nlam
					Ce_cont(ilam,ir,i)=0d0
					do ig=1,ng
						Ce_cont(ilam,ir,i)=Ce_cont(ilam,ir,i)+wgg(ig)*Ca(ilam,ig,ir,i,0)
					enddo
				enddo
				do imol=1,nmol
					if(includemol(imol)) then
						k=k+1
						do ig=1,ng
							do ilam=1,nlam
								Ca_mol(ilam,ig,k,ir,i)=Cabs_mol(ig,ilam,imol,ir)
								Ce_cont(ilam,ir,i)=Ce_cont(ilam,ir,i)-wgg(ig)*Ca_mol(ilam,ig,k,ir,i)
							enddo
						enddo
					endif
				enddo
				nmol_count=k
			enddo
			if(emisspec) call ComputeScatter(BBr(1:nlam,0:nr),Si(1:nlam,1:ng,0:nr,1:nnu0,i),Ca(1:nlam,1:ng,1:nr,i,0),Cs(1:nlam,1:nr,i))
			if(computealbedo) then
				BBr(1:nlam,0:nr)=0d0
				Fstar_temp(1:nlam)=Fstar(1:nlam)
				Fstar=1d0
				call ComputeScatter(BBr(1:nlam,0:nr),SiSc(1:nlam,1:ng,0:nr,1:nnu0,i),Ca(1:nlam,1:ng,1:nr,i,0),Cs(1:nlam,1:nr,i))
				Fstar(1:nlam)=Fstar_temp(1:nlam)
			endif

			Pb(1)=P(1)
			do ir=2,nr
				Pb(ir)=sqrt(P(ir-1)*P(ir))
			enddo
			Pb(nr+1)=0d0
			do ilam=1,nlam
				do ig=1,ng
					do ir=nr,1,-1
						do ivel=-nvel,nvel
							if(ir.eq.nr) then
								v=(P(ir)-Pb(ir+1))*1d6/grav(ir)
								tau1=v*(Ca(ilam,ig,ir,i,ivel)+Cs(ilam,ir,i))/dens(ir)
							else
								v=(Pb(ir+1)-P(ir+1))*1d6/grav(ir+1)
								tau1=v*(Ca(ilam,ig,ir+1,i,ivel)+Cs(ilam,ir+1,i))/dens(ir+1)
								v=(P(ir)-Pb(ir+1))*1d6/grav(ir)
								tau1=tau1+v*(Ca(ilam,ig,ir,i,ivel)+Cs(ilam,ir,i))/dens(ir)
							endif
							if(.not.tau1.gt.0d0) then
								tau1=0d0
							endif
							tau1=tau1+1d-10
							tau1=1d0/(1d0/tau1+1d-10)
							dtauR_nu(ilam,ig,i,ir,ivel)=tau1
						enddo
					enddo
				enddo
			enddo
		else
			R3D(i,1:nr+1)=R3D(n3D,1:nr+1)
			T3D(i,0:nr)=T3D(n3D,0:nr)
			mixrat3D(i,1:nr,1:nmol)=mixrat3D(n3D,1:nr,1:nmol)
			Ca(1:nlam,1:ng,1:nr,i,-nvel:nvel)=Ca(1:nlam,1:ng,1:nr,n3D,-nvel:nvel)
			Cs(1:nlam,1:nr,i)=Cs(1:nlam,1:nr,n3D)
			Ce_cont(1:nlam,1:nr,i)=Ce_cont(1:nlam,1:nr,n3D)
			Ca_mol(1:nlam,1:ng,1:nmol_count,1:nr,i)=Ca_mol(1:nlam,1:ng,1:nmol_count,1:nr,n3D)
			dtauR_nu(1:nlam,1:ng,i,1:nr,-nvel:nvel)=dtauR_nu(1:nlam,1:ng,n3D,1:nr,-nvel:nvel)
			Si(1:nlam,1:ng,0:nr,1:nnu0,i)=Si(1:nlam,1:ng,0:nr,1:nnu0,n3D)
			if(computealbedo) SiSc(1:nlam,1:ng,0:nr,1:nnu0,i)=SiSc(1:nlam,1:ng,0:nr,1:nnu0,n3D)
			local_albedo(i)=local_albedo(n3D)
		endif
		if(.not.retrieval.and..not.dopostequalweights) then
			call SetOutputMode(.true.)
			open(unit=20,file=trim(outputdir) // "mixrat" // trim(int2string(i,'(i0.3)')),FORM="FORMATTED",ACCESS="STREAM")
			write(20,'("#",a9,a13,a13)') "T[K]","P[bar]","Kzz[cm^2/s]"
			do j=1,nr
				write(20,'(f10.3,2es13.3E3)') T(j),P(j),ComputeKzz(P(j),T(j),dens(j),Hp(j),complexKzz)
			enddo
			close(unit=20)
		endif
		do ir=1,nr
			Tprev3D(ir)=maxval(T3D(i3D:n3D,ir))
		enddo
		call tellertje_perc(n3D-i+1,n3D)
	enddo

	scale_prev=scale
	albedo_day=0d0
	tot1=0d0
	tot2=0d0
	scale=0d0
	do i=1,nlong-1
		lo=(long(i)+long(i+1))/2d0-pi
		if(abs(lo).le.pi/2d0) then
			lo=cos(lo)
			do j=1,nlatt-1
				la=cos((latt(j)+latt(j+1))/2d0-pi/2d0)
				albedo_day=albedo_day+local_albedo(ibeta(i,j))*la*lo
				tot1=tot1+la*lo
			enddo
		endif
		do j=1,nlatt-1
			la=cos((latt(j)+latt(j+1))/2d0-pi/2d0)
			scale=scale+(1d0-local_albedo(ibeta(i,j)))*la*beta3D(ibeta(i,j))
			tot2=tot2+la
		enddo
	enddo
	albedo_day=albedo_day/tot1
	scale=4d0*scale/tot2
	scale=(1d0-albedo_day)/scale
	call output("Day-side albedo: " // trim(dbl2string(albedo_day,'(f5.3)')))
	if(nalbedo_iter.gt.1) call output("Scaling emission by: " // trim(dbl2string(scale,'(f5.3)')))
	if(abs(scale-scale_prev)/(scale+scale_prev).lt.0.5d-2) exit
	enddo

	Rmax=Rmax*1.001
	R3D(1:n3D,nr+2)=Rmax
	R3DC(1:n3D,1:nr+1)=sqrt(R3D(1:n3D,1:nr+1)*R3D(1:n3D,2:nr+2))
	R3DC(1:n3D,nr+2)=Rmax
	R3D2(1:n3D,1:nr+2)=R3DC(1:n3D,1:nr+2)**2
	do ilam=1,nlam
		do ig=1,ng
			do ivel=-nvel,nvel
				dtauR_nu(ilam,ig,1:n3D,1:nr,ivel)=dtauR_nu(ilam,ig,1:n3D,1:nr,ivel)/(R3DC(1:n3D,2:nr+1)-R3DC(1:n3D,1:nr))
			enddo
		enddo
	enddo	

	if(.not.retrieval.and..not.dopostequalweights) then
		open(unit=20,file=trim(outputdir) // "surfacetemp3D.dat",FORM="FORMATTED",ACCESS="STREAM")
		do j=1,nlatt-1
			write(20,*) T3D(ibeta(1:nlong-1,j),0)
		enddo
		close(unit=20)
		if(fulloutput3D) then
			do ir=1,nr
				open(unit=20,file=trim(outputdir) // "temp3D_P" // trim(dbl2string(P(ir),'(es8.2)')) // ".dat",FORM="FORMATTED",ACCESS="STREAM")
				do j=1,nlatt-1
					write(20,*) T3D(ibeta(1:nlong-1,j),ir)
				enddo
				close(unit=20)
				open(unit=20,file=trim(outputdir) // "cloud3D_P" // trim(dbl2string(P(ir),'(es8.2)')) // ".dat",FORM="FORMATTED",ACCESS="STREAM")
				do j=1,nlatt-1
					write(20,*) cloud3D(ibeta(1:nlong-1,j),ir)
				enddo
				close(unit=20)
			enddo
		endif
	endif

	call output("Raytracing over the planet disk in 3D")

	do ilam=1,nlam
		cloudtau(1,ilam)=0d0
		do ir=nr,1,-1
			do icloud=1,nclouds
				cloudtau(1,ilam)=cloudtau(1,ilam)+(R(ir+1)-R(ir))*Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
				cloudtau(1,ilam)=cloudtau(1,ilam)+(R(ir+1)-R(ir))*Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
			enddo
		enddo
	enddo


	if(emisspec) then

	if(fulloutput3D) then
		PTaverage3D=0d0
		mixrat_average3D=0d0
	endif
	ndisk=20
	nsub=0

	nrtrace=(nr-1)*nsub+ndisk
	nptrace=(nlatt-1)*2
	if(actually1D.and.nphase.eq.1.and.theta_phase(1).eq.180d0) nptrace=1
	if(vrot0.ne.0d0) then
		do_rot=.true.
		nptrace=max(10,nptrace)
	else
		do_rot=.false.
	endif
	
	if(makeimage) then
		nrtrace=nrtrace*4
		nptrace=(nlatt-1)*8
		allocate(rphi_image(nlam,nrtrace,nptrace))
		allocate(xy_image(nx_im,nx_im,nlam))
	endif
	allocate(rtrace(nrtrace),wrtrace(nrtrace))

	call gauleg(0d0,Rmax,rtrace,wrtrace,nrtrace)

	allocate(fluxp(nlam))
	npc=nphase
	call tellertje_perc(0,npc)
	do ipc=1,npc
	if(fulloutput3D) PTaverage3D(ipc,1:nr)=0d0
	theta=2d0*pi*theta_phase(ipc)/360d0
	if(theta.gt.2d0*pi) theta=theta-2d0*pi
	fluxp=0d0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(irtrace,iptrace,A,phi,rr,y,z,x,vx,vy,vz,la,lo,i1,i2,i3,edgeNR,j,i,inu,fluxp_omp,w1,w2,SiR0,SiR1,tau0,
!$OMP&			i1next,i2next,i3next,edgenext,freq0,tot,v,ig,ilam,tau1,fact,exp_tau1,contr,ftot,alb_omp,SiRalb0,SiRalb1,
!$OMP&			vrot,dlam_rot,ivel)
!$OMP& SHARED(theta,fluxp,nrtrace,rtrace,wrtrace,nptrace,Rmax,nr,freq,ibeta,fulloutput3D,Rplanet,computeT,dtauR_nu,vrot0,lam,do_rot,
!$OMP&			rphi_image,makeimage,nnu0,nlong,nlatt,R3D,planet_albedo,SiSc,computealbedo,orbit_inc,maxtau,R3DC,computelam,nvel,
!$OMP&			Ca,Cs,wgg,Si,R3D2,latt,long,T,ng,nlam,ipc,PTaverage3D,mixrat_average3D,T3D,mixrat3D,nmol,surface_emis,lamemis)
	allocate(fact(nlam,ng))
	allocate(fluxp_omp(nlam))
	allocate(ftot(nlam),SiR0(nlam,ng))
	if(computealbedo) allocate(alb_omp(nlam),SiRalb0(nlam,ng))
	fluxp_omp=0d0
	if(computealbedo) alb_omp=0d0
!$OMP DO SCHEDULE(DYNAMIC,1)
	do irtrace=1,nrtrace
		A=2d0*pi*rtrace(irtrace)*wrtrace(irtrace)/real(nptrace)
		do iptrace=1,nptrace
			ftot=0d0
			fact=1d0
c Note we are here using the symmetry between North and South
			if(nptrace.eq.1) then
				if(2*(nlatt/2).eq.nlatt) then
					phi=pi*(real(iptrace)-0.5)/real((nlatt-1)*2)
				else
					phi=pi*(real(iptrace)-0.5)/real(nlatt-1)
				endif
			else
				phi=2d0*pi*(real(iptrace)-0.5)/real(nptrace)
			endif
			rr=rtrace(irtrace)
			y=rr*cos(phi)
			z=rr*sin(phi)
			x=sqrt(Rmax**2-y**2-z**2)
			if(do_rot) then
				vrot=vrot0*y/Rmax
				vx=vrot0*y/Rmax
				vy=vrot0*x/Rmax
				vz=0d0
				rr=sqrt(vx**2+vy**2+vz**2)
				call rotateY3D(vx,vy,vz,pi/2d0-orbit_inc)
				vrot=vx*rr/sqrt(vx**2+vy**2+vz**2)
				ivel=real(nvel)*vrot/vrot0
				ivel=(vrot/vrot0)*(real(nvel)+0.5)+0.5*sign(1d0,vrot/vrot0)
				if(ivel.lt.-nvel) ivel=-nvel
				if(ivel.gt.nvel) ivel=nvel
			else
				ivel=0
			endif
			vx=-1d0
			vy=0d0
			vz=0d0
			call rotateY3D(x,y,z,pi/2d0-orbit_inc)
			rr=sqrt(x**2+y**2+z**2)
			x=x*Rmax/rr
			y=y*Rmax/rr
			z=z*Rmax/rr
			call rotateZ3D(x,y,z,theta)
			rr=sqrt(x**2+y**2+z**2)
			x=x*Rmax/rr
			y=y*Rmax/rr
			z=z*Rmax/rr
			call rotateY3D(vx,vy,vz,pi/2d0-orbit_inc)
			rr=sqrt(vx**2+vy**2+vz**2)
			vx=vx/rr
			vy=vy/rr
			vz=vz/rr
			call rotateZ3D(vx,vy,vz,theta)
			rr=sqrt(vx**2+vy**2+vz**2)
			vx=vx/rr
			vy=vy/rr
			vz=vz/rr
			la=acos(z/sqrt(x**2+y**2+z**2))
			lo=acos(x/sqrt(x**2+y**2))
			if(y.gt.0d0) lo=2d0*pi-lo
			do i2=1,nlong-1
				if(lo.ge.long(i2).and.lo.lt.long(i2+1)) exit
			enddo
			do i3=1,nlatt-1
				if(la.ge.latt(i3).and.la.lt.latt(i3+1)) exit
			enddo
			i1=nr+1
			edgeNR=2
			i=ibeta(i2,i3)
			la=(-x/sqrt(x**2+y**2+z**2))
			if(la.lt.0d0) then
				inu=nnu0
			else
				inu=(real(nnu0-2)*la+0.5d0)+1
				if(inu.lt.1) inu=1
				if(.not.inu.lt.nnu0-1) inu=nnu0-1
			endif
			if(fulloutput3D) then
				PTaverage3D(ipc,1:nr)=PTaverage3D(ipc,1:nr)+T3D(i,1:nr)*A
				mixrat_average3D(ipc,1:nr,1:nmol)=mixrat_average3D(ipc,1:nr,1:nmol)+mixrat3D(i,1:nr,1:nmol)*A
			endif
			j=0
			tau0=0d0
			SiR0=0d0
			if(computealbedo) SiRalb0=0d0
1			continue
			call TravelSph(x,y,z,vx,vy,vz,edgeNR,i1,i2,i3,v,i1next,i2next,i3next,edgenext,nr,nlong,nlatt)
			if(i1.le.nr) then
				i=ibeta(i2,i3)
				la=(-x/sqrt(x**2+y**2+z**2))
				if(la.lt.0d0) then
					inu=nnu0
				else
					inu=(real(nnu0-2)*la+0.5d0)+1
					if(inu.lt.1) inu=1
					if(.not.inu.lt.nnu0-1) inu=nnu0-1
				endif
				do ilam=1,nlam
					if(lamemis(ilam).and.computelam(ilam)) then
					do ig=1,ng
						tau0=0d0
						tau1=v*dtauR_nu(ilam,ig,i,i1,ivel)
						exp_tau1=exp(-tau1)
						rr=sqrt((x+v*vx)**2+(y+v*vy)**2+(z+v*vz)**2)
						if(i1.lt.nr) then
							w1=(R3DC(i,i1+1)-rr)/(R3DC(i,i1+1)-R3DC(i,i1))
							w2=1d0-w1
							SiR1=w1*Si(ilam,ig,i1,inu,i)+w2*Si(ilam,ig,i1+1,inu,i)
							if(computealbedo) SiRalb1=w1*SiSc(ilam,ig,i1,inu,i)+w2*SiSc(ilam,ig,i1+1,inu,i)
						else
							SiR1=Si(ilam,ig,nr,inu,i)
							if(computealbedo) SiRalb1=SiSc(ilam,ig,nr,inu,i)
						endif
						call ComputeI12(tau1,tau0,SiR1,SiR0(ilam,ig),contr)
						ftot(ilam)=ftot(ilam)+A*wgg(ig)*contr*fact(ilam,ig)
						if(computealbedo) then
							call ComputeI12(tau1,tau0,SiRalb1,SiRalb0(ilam,ig),contr)
							alb_omp(ilam)=alb_omp(ilam)+A*wgg(ig)*contr*fact(ilam,ig)
							SiRalb0(ilam,ig)=SiRalb1
						endif
						fact(ilam,ig)=fact(ilam,ig)*exp_tau1
						SiR0(ilam,ig)=SiR1
					enddo
					endif
				enddo
			endif
			if(i1next.le.0) then
				i=ibeta(i2,i3)
				la=(-x/sqrt(x**2+y**2+z**2))
				if(la.lt.0d0) then
					inu=nnu0
				else
					inu=(real(nnu0-2)*la+0.5d0)+1
					if(inu.lt.1) inu=1
					if(.not.inu.lt.nnu0-1) inu=nnu0-1
				endif
				do ilam=1,nlam
					if(lamemis(ilam).and.computelam(ilam)) then
					do ig=1,ng
						contr=Si(ilam,ig,0,inu,i)
						ftot(ilam)=ftot(ilam)+A*wgg(ig)*contr*fact(ilam,ig)
						if(computealbedo) alb_omp(ilam)=alb_omp(ilam)+A*wgg(ig)*SiSc(ilam,ig,0,inu,i)*fact(ilam,ig)
					enddo
					endif
				enddo
				goto 2
			endif
			if(i1next.ge.nr+2) goto 2
			x=x+v*vx
			y=y+v*vy
			z=z+v*vz
			if(i2next.gt.nlong-1.or.i3next.gt.nlatt-1.or.i2next.lt.1.or.i3next.lt.1) goto 2
			if(ibeta(i2,i3).ne.ibeta(i2next,i3next)) then
				rr=x**2+y**2+z**2
				i=ibeta(i2next,i3next)
				if((.not.rr.ge.R3D2(i,1)).or.(.not.rr.le.R3D2(i,nr+2))) goto 2
				do i1=1,nr+1
					if(rr.gt.R3D2(i,i1).and.rr.le.R3D2(i,i1+1)) exit
				enddo
				if(i1.ge.nr+2) goto 2
			else
				i1=i1next
			endif
			i2=i2next
			i3=i3next
			edgeNR=edgenext
			j=j+1
			if(j.lt.nr*2*max(nlatt,nlong)) goto 1
2			continue
			fluxp_omp(1:nlam)=fluxp_omp(1:nlam)+ftot(1:nlam)
			if(makeimage) rphi_image(1:nlam,irtrace,iptrace)=ftot(1:nlam)
		enddo
	enddo
!$OMP END DO
!$OMP CRITICAL
	fluxp(1:nlam)=fluxp(1:nlam)+fluxp_omp(1:nlam)
	if(computealbedo) then
		planet_albedo(ipc,1:nlam)=planet_albedo(ipc,1:nlam)+alb_omp(1:nlam)
		deallocate(alb_omp,SiRalb0)
	endif
	deallocate(fluxp_omp)
!$OMP END CRITICAL
	deallocate(fact)
	deallocate(ftot)
	deallocate(SiR0)
!$OMP FLUSH
!$OMP END PARALLEL
	fluxp=fluxp*1d23/distance**2
	phase(ipc,0,1:nlam)=fluxp(1:nlam)
	flux(0,1:nlam)=0d0

	call tellertje_perc(ipc,npc)
	if(fulloutput3D) then
		PTaverage3D(ipc,1:nr)=PTaverage3D(ipc,1:nr)/(pi*Rmax**2)
		mixrat_average3D(ipc,1:nr,1:nmol)=mixrat_average3D(ipc,1:nr,1:nmol)/(pi*Rmax**2)
		open(unit=25,file=trim(outputdir) // "mixrat_phase" // trim(int2string(int(theta_phase(ipc)),'(i0.3)')),FORM="FORMATTED",ACCESS="STREAM")
		do ir=1,nr
			j=0
			do i=1,nmol
				if(includemol(i)) then
					j=j+1
					tmp(j)=mixrat_average3D(ipc,ir,i)
				endif
			enddo
			write(25,*) PTaverage3D(ipc,ir),P(ir),tmp(1:j)
		enddo
	endif
c	tot=0d0
c	do ilam=1,nlam
c		if(lamemis(ilam).and.computelam(ilam)) tot=tot+phase(ipc,0,ilam)*dfreq(ilam)
c	enddo
c	tot=tot*distance**2/1d23
c	print*,theta_phase(ipc),tot/(2d0*pi*(((pi*kb*Tstar)**4)/(15d0*hplanck**3*clight**3))*Rstar**2*Rplanet**2/(Dplanet**2))
c	print*,Tstar*(Rstar/Dplanet)**0.5

	if(makeimage) then
		file=trim(outputdir) // "image" //  trim(int2string(int(theta_phase(ipc)),'(i0.3)')) // ".fits"
		call output("Creating image: " // trim(file))
		Rmin_im=0d0
		xy_image=0d0
		do irtrace=1,nrtrace
			call tellertje(irtrace,nrtrace)
			A=2d0*pi*rtrace(irtrace)*wrtrace(irtrace)
			Rmax_im=sqrt((A+pi*Rmin_im**2)/pi)
			do iptrace=1,nptrace
				ni=125d0*real(nx_im*nx_im)/real(nrtrace*nptrace)
				do i=1,ni
					rr=sqrt(Rmin_im**2+random(idum)*(Rmax_im**2-Rmin_im**2))
					phi=2d0*pi*(real(iptrace)-random(idum))/real(nptrace)
					x=rr*cos(phi)
					y=rr*sin(phi)
					ix=(x+Rmax)*real(nx_im)/(2d0*Rmax)+1
					iy=(y+Rmax)*real(nx_im)/(2d0*Rmax)+1
					if(ix.lt.1) ix=1
					if(iy.lt.1) iy=1
					if(ix.gt.nx_im) ix=nx_im
					if(iy.gt.nx_im) iy=nx_im
					xy_image(ix,iy,1:nlam)=xy_image(ix,iy,1:nlam)+rphi_image(1:nlam,irtrace,iptrace)/real(ni*2)
c					ix=(x+Rmax)*real(nx_im)/(2d0*Rmax)+1
c					iy=(-y+Rmax)*real(nx_im)/(2d0*Rmax)+1
c					if(ix.lt.1) ix=1
c					if(iy.lt.1) iy=1
c					if(ix.gt.nx_im) ix=nx_im
c					if(iy.gt.nx_im) iy=nx_im
c					xy_image(ix,iy,1:nlam)=xy_image(ix,iy,1:nlam)+rphi_image(1:nlam,irtrace,iptrace)/real(ni*2)
				enddo
			enddo
			Rmin_im=Rmax_im
		enddo
		ni=nlam
		allocate(lgrid(nlam))
		lgrid(1:nlam)=lam(1:nlam)
		ni=0
		do i=1,nlam
			if(computelam(i)) then
				ni=ni+1
				if(ni.ne.i) then
					xy_image(1:nx_im,1:nx_im,ni)=xy_image(1:nx_im,1:nx_im,i)
					lgrid(ni)=lgrid(i)
				endif
			endif
		enddo
		
		call writefitsfile(file,xy_image(1:nx_im,1:nx_im,1:ni),ni,nx_im)

		file=trim(outputdir) // "imageRGB" //  trim(int2string(int(theta_phase(ipc)),'(i0.3)')) // ".fits"
		call output("Creating image: " // trim(file))
		i=0
		allocate(maxdet(nx_im*nx_im,3))
		do ix=1,nx_im
			do iy=1,nx_im
				call SPECTOXYZ(xy_image(ix,iy,1:ni),lgrid(1:ni),ni,x,y,z)
				xy_image(ix,iy,1)=x
				xy_image(ix,iy,2)=y
				xy_image(ix,iy,3)=z
				i=i+1
				if(x.ge.y.and.x.ge.z) then
					x=0d0
				else if(y.ge.x.and.y.ge.z) then
					y=0d0
				else
					z=0d0
				endif
				maxdet(i,1)=x
				maxdet(i,2)=y
				maxdet(i,3)=z
			enddo
		enddo
		deallocate(lgrid)
		if(ipc.eq.1.or..not.makemovie) then
		ni=nx_im*nx_im
		do i=1,ni/20
			maximage=0d0
			do j=1,ni
				if(maxdet(j,1).gt.maximage) then
					k=j
					maximage=maxdet(j,1)
				endif
				if(maxdet(j,2).gt.maximage) then
					k=j
					maximage=maxdet(j,2)
				endif
				if(maxdet(j,3).gt.maximage) then
					k=j
					maximage=maxdet(j,3)
				endif
			enddo
			if(k.gt.ni) k=ni
			maxdet(k,1:3)=0d0
		enddo
		endif
		deallocate(maxdet)
		xy_image(1:nx_im,1:nx_im,1:3)=xy_image(1:nx_im,1:nx_im,1:3)/maximage
		do ix=1,nx_im
			do iy=1,nx_im
				call XYZTORGB(xy_image(ix,iy,1),xy_image(ix,iy,2),xy_image(ix,iy,3),x,y,z)
				xy_image(ix,iy,1)=x
				xy_image(ix,iy,2)=y
				xy_image(ix,iy,3)=z
			enddo
		enddo
		ni=3
		call writefitsfile(file,xy_image,ni,nx_im)
		file=trim(outputdir) // "imageRGB" //  trim(int2string(int(theta_phase(ipc)),'(i0.3)')) // ".ppm"
		call output("Creating image: " // trim(file))
		call writeppmfile(file,xy_image,ni,nx_im)
	endif

	enddo
	
	if(.not.retrieval.or.doRing) then
		Lplanet=0d0
		do ilam=1,nlam
			if(doRing.and..not.computeT) then
				if(RTgridpoint(ilam)) Lplanet=Lplanet+dfreq(ilam)*phase(1,0,ilam)
			else if(computelam(ilam).and.lamemis(ilam)) then
				Lplanet=Lplanet+dfreq(ilam)*phase(1,0,ilam)
			endif
		enddo
		TeffPoutput=(Lplanet*distance**2*1e-23/(pi*Rplanet**2*((2d0*(pi*kb)**4)/(15d0*hplanck**3*clight**3))))**0.25d0
		call output("Teff: " // dbl2string(TeffPoutput,'(f10.2)') // "K" )
	endif

	if(makeimage) then
		deallocate(rphi_image)
		deallocate(xy_image)
	endif
	deallocate(rtrace,wrtrace)
	deallocate(fluxp)
	endif

	if(fulloutput3D) then
		PTaverage3D(0,1:nr)=0d0
		do j=1,nlatt
			i=nlong/4+1
			PTaverage3D(0,1:nr)=PTaverage3D(0,1:nr)+T3D(ibeta(i,j),1:nr)/real(nlatt*2)
			mixrat_average3D(0,1:nr,1:nmol)=mixrat_average3D(0,1:nr,1:nmol)+mixrat3D(ibeta(i,j),1:nr,1:nmol)/real(nlatt*2)
			i=3*nlong/4+1
			PTaverage3D(0,1:nr)=PTaverage3D(0,1:nr)+T3D(ibeta(i,j),1:nr)/real(nlatt*2)
			mixrat_average3D(0,1:nr,1:nmol)=mixrat_average3D(0,1:nr,1:nmol)+mixrat3D(ibeta(i,j),1:nr,1:nmol)/real(nlatt*2)
		enddo
		open(unit=25,file=trim(outputdir) // "mixrat_transit",FORM="FORMATTED",ACCESS="STREAM")
		do ir=1,nr
			j=0
			do i=1,nmol
				if(includemol(i)) then
					j=j+1
					tmp(j)=mixrat_average3D(0,ir,i)
				endif
			enddo
			write(25,*) PTaverage3D(0,ir),P(ir),tmp(1:j)
		enddo
	endif

	if(transspec) then

	R3D2(1:n3D,1:nr+2)=R3D(1:n3D,1:nr+2)**2

	ndisk=4
	nsub=3

	nrtrace=(nr-1)*nsub+ndisk
	if(2*(nlatt/2).eq.nlatt) then
		nptrace=(nlatt-1)*2
	else
		nptrace=nlatt-1
	endif
	if(actually1D) nptrace=1
	allocate(rtrace(nrtrace))
	if(n3D.eq.2) nptrace=2

	k=0
	do i=1,ndisk
		k=k+1
		rtrace(k)=R3D(iRmax,1)*real(i-1)/real(ndisk)
	enddo
	do i=1,nr-1
		do j=1,nsub
			k=k+1
			rtrace(k)=R3D(iRmax,i)+(R3D(iRmax,i+1)-R3D(iRmax,i))*real(j-1)/real(nsub)
		enddo
	enddo
	obsA=0d0
	obsA_split=0d0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(vv,tau,tauc,Afact,hit,irtrace,A,iptrace,phi,rr,x,y,z,vx,vy,vz,la,lo,i1,i2,i3,
!$OMP&			edgeNR,v,i1next,i2next,i3next,edgenext,i,imol,ir,ig,obsA_omp,obsA_split_omp)
!$OMP& SHARED(nrtrace,nptrace,nmol_count,nr,nlam,ng,rtrace,ibeta,R3D2,Ce_cont,Ca_mol,wgg,obsA,latt,long,Rmax,
!$OMP&          nnu0,n3D,nlong,nlatt,obsA_split)
	allocate(obsA_omp(nlam))
	allocate(obsA_split_omp(nlam,2))
	allocate(vv(nr,n3D))
	allocate(tau(nlam,ng))
	allocate(tauc(nlam),Afact(nlam))
	allocate(hit(nr+2,n3D))
	obsA_omp=0d0
	obsA_split_omp=0d0
!$OMP DO
	do irtrace=1,nrtrace-1
		A=pi*(rtrace(irtrace+1)**2-rtrace(irtrace)**2)/real(nptrace)
		do iptrace=1,nptrace
c Note we use the symmetry of the North and South here!
			if(nptrace.gt.2) then
				if(2*(nlatt/2).eq.nlatt) then
					phi=pi*(real(iptrace)-0.5)/real((nlatt-1)*2)
				else
					phi=pi*(real(iptrace)-0.5)/real(nlatt-1)
				endif
			else if(nptrace.eq.2) then
				phi=pi*(real(iptrace-1))/real(nptrace-1)
			else
				phi=2d0*pi*(real(iptrace)-0.5)/real(nptrace)
			endif
			rr=0.5d0*(rtrace(irtrace)+rtrace(irtrace+1))
			y=rr*cos(phi)
			z=rr*sin(phi)
			x=sqrt(Rmax**2-y**2-z**2)
			vx=-1d0
			vy=0d0
			vz=0d0
			la=acos(z/sqrt(x**2+y**2+z**2))
			lo=acos(x/sqrt(x**2+y**2))
			if(y.gt.0d0) lo=2d0*pi-lo
			do i2=1,nlong-1
				if(lo.ge.long(i2).and.lo.le.long(i2+1)) exit
			enddo
			do i3=1,nlatt-1
				if(la.ge.latt(i3).and.la.le.latt(i3+1)) exit
			enddo
			i1=nr+1
			edgeNR=2
			tauc=0d0
			vv=0d0
			hit=.false.
3			continue
			call TravelSph(x,y,z,vx,vy,vz,edgeNR,i1,i2,i3,v,i1next,i2next,i3next,edgenext,nr,nlong,nlatt)
			if(i1next.le.0) then
				tauc(1:nlam)=1d5
				goto 4
			endif
			if(i1next.ge.nr+2) goto 4
			if(i1.le.nr) then
				i=ibeta(i2,i3)
				vv(i1,i)=vv(i1,i)+v
				hit(i1,i)=.true.
				tauc(1:nlam)=tauc(1:nlam)+v*Ce_cont(1:nlam,i1,i)
			endif
			x=x+v*vx
			y=y+v*vy
			z=z+v*vz
			if(i2next.gt.nlong.or.i3next.gt.nlatt.or.i2next.lt.1.or.i3next.lt.1) goto 4
			if(ibeta(i2,i3).ne.ibeta(i2next,i3next)) then
				rr=x**2+y**2+z**2
				i=ibeta(i2next,i3next)
				if((.not.rr.ge.R3D2(i,1)).or.(.not.rr.le.R3D2(i,nr+2))) goto 4
				do i1=1,nr+1
					if(rr.gt.R3D2(i,i1).and.rr.le.R3D2(i,i1+1)) exit
				enddo
			else
				i1=i1next
			endif
			i2=i2next
			i3=i3next
			edgeNR=edgenext
			goto 3
4			continue
			Afact(1:nlam)=exp(-tauc(1:nlam))
			do imol=1,nmol_count
				tau=0d0
				do i=1,n3D
					do ir=1,nr
						if(hit(ir,i)) then
							tau(1:nlam,1:ng)=tau(1:nlam,1:ng)+vv(ir,i)*Ca_mol(1:nlam,1:ng,imol,ir,i)
						endif
					enddo
				enddo
				tauc=0d0
				do ig=1,ng
					tauc(1:nlam)=tauc(1:nlam)+wgg(ig)*exp(-tau(1:nlam,ig))
				enddo
				Afact=Afact*tauc
			enddo
			do imol=1,nmol_count
				tauc=0d0
			enddo
			obsA_omp(1:nlam)=obsA_omp(1:nlam)+A*(1d0-Afact(1:nlam))
			if(y.lt.0d0) then
				obsA_split_omp(1:nlam,1)=obsA_split_omp(1:nlam,1)+A*(1d0-Afact(1:nlam))
			else
				obsA_split_omp(1:nlam,2)=obsA_split_omp(1:nlam,2)+A*(1d0-Afact(1:nlam))
			endif
		enddo
	enddo
!$OMP END DO
!$OMP CRITICAL
	obsA(0,1:nlam)=obsA(0,1:nlam)+obsA_omp(1:nlam)
	obsA_split(1:nlam,1:2)=obsA_split(1:nlam,1:2)+obsA_split_omp(1:nlam,1:2)
!$OMP END CRITICAL
	deallocate(obsA_omp,obsA_split_omp)
	deallocate(vv)
	deallocate(tau)
	deallocate(tauc,Afact)
	deallocate(hit)
!$OMP FLUSH
!$OMP END PARALLEL
	deallocate(rtrace)

	endif
	
	deallocate(Ca,Cs,BBr,Si)
	if(computealbedo) deallocate(SiSc)
	deallocate(Ca_mol,Ce_cont)
	deallocate(R3D)
	deallocate(R3D2)
	deallocate(R3DC)
	deallocate(T3D,mixrat3D)
	deallocate(local_albedo)
	if(.not.retrieval.and.fulloutput3D) deallocate(cloud3D)

	if(retrieval.or.dopostequalweights) call SetOutputMode(.true.)
	
	return
	end
	

	subroutine SolveIj3D(tauR_in,Si,ftot,nr)
	IMPLICIT NONE
	integer ir,nr,info,nnr
	real*8 tauR_in(nr),Si(nr),fact,d,ftot
	real*8 x(nr+2),y(nr+2),tauR(0:nr+1),Ma(nr+2),Mb(nr+2),Mc(nr+2)

	nnr=1
	x=0d0
	do ir=1,nr
		tauR(nr-ir+1)=tauR_in(ir)
		x(nr+2-ir)=Si(ir)
	enddo
	tauR(0)=tauR(1)**2/tauR(2)
	tauR(nr+1)=0d0

	Ma=0d0
	Mb=0d0
	Mc=0d0
	ir=1
	do ir=1,nr
		fact=1d0/(0.5d0*(tauR(ir+1)+tauR(ir))-0.5d0*(tauR(ir)+tauR(ir-1)))
		Mb(ir+1)=1d0+fact*(1d0/(tauR(ir+1)-tauR(ir))+1d0/(tauR(ir)-tauR(ir-1)))
		Ma(ir+1)=-fact*1d0/(tauR(ir)-tauR(ir-1))
		Mc(ir+1)=-fact*1d0/(tauR(ir+1)-tauR(ir))
	enddo
	Mb(1)=1d0+1d0/(tauR(0)-tauR(1))
	Mc(1)=-1d0/(tauR(0)-tauR(1))
	Ma(nr+2)=1d0/(tauR(nr))
	Mb(nr+2)=-1d0-1d0/(tauR(nr))

	info=0

	call dgtsv(nr+2,nnr,Ma(2:nr+2),Mb(1:nr+2),Mc(1:nr+1),x,nr+2,info)
	ftot=x(nr+2)

	return
	end
	

	subroutine MapPar3D()
	use GlobalSetup
	use ReadKeywords
	use Constants
	IMPLICIT NONE
	type(SettingKey) key
	character*1000 readline
	integer i
	real*8 Z0
	
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

	if(dochemistry) metallicity=metallicity0
	do i=1,n_Par3D
		readline=trim(Par3D(i)%keyword) // "=" // trim(dbl2string(Par3D(i)%x,'(es14.7)'))
		call get_key_value(readline,key%key,key%key1,key%key2,key%orkey1,key%orkey2,key%value,key%nr1,key%nr2,key%key2d)
		call ReadAndSetKey(key)
	enddo
	call ConvertUnits()
	if(dochemistry) metallicity0=metallicity

	return
	end
	

	subroutine Setup3D(beta,long,latt,nlong,nlatt,Kxx,Kyy,vxx,powvxx,night2day,pole2eq,f,betamin,betamax,tidallock)
	IMPLICIT NONE
	integer i,j,nlong,nlatt
	real*8 pi
	parameter(pi=3.14159265358979323846264338328d0)
	real*8 long(nlong),latt(nlatt)	!(Lambda, Phi)
	real*8 beta(nlong,nlatt),la,lo,albedo,p,f,tot1,tot2,b1,b2
	real*8 Kxx,vxx,betamin,betamax,beta1(nlong),night2day,Kyy,powvxx,pole2eq
	logical tidallock

	do i=1,nlong
		long(i)=-pi+2d0*pi*real(i-1)/real(nlong-1)
	enddo
	do j=1,nlatt
		latt(j)=-pi/2d0+pi*real(j-1)/real(nlatt-1)
	enddo
	
	betamin=1d0
	betamax=0d0

	if(tidallock) then
		call DiffuseBeta(beta,long,latt,nlong,nlatt,Kxx,Kyy,vxx,powvxx,night2day)
	else
		call DiffuseBetaRotate(beta,long,latt,nlong,nlatt,Kxx,Kyy,vxx,powvxx,pole2eq)
	endif
	beta=beta*f

	do i=1,nlong-1
		do j=1,nlatt-1
			if(beta(i,j).gt.betamax) betamax=beta(i,j)
			if(beta(i,j).lt.betamin) betamin=beta(i,j)
		enddo
	enddo
	
	latt=latt+pi/2d0
	long=long+pi

	return
	end

	subroutine DetermineShift(long,beta,nlong,nlatt,shift)
	IMPLICIT NONE
	integer nlong,nlatt,i,j,k
	real*8 shift,beta(nlong,nlatt),max,x1,x2,x3,y1,y2,y3,pi,long(nlong)
	parameter(pi=3.14159265358979323846264338328d0)
	
	j=nlatt/2
	max=0d0
	k=nlong/2
	do i=1,nlong-1
		if(beta(i,j).gt.max) then
			max=beta(i,j)
			k=i
		endif
	enddo
	if(k.eq.1) then
		y1=beta(nlong-1,j)
		y2=beta(k,j)
		y3=beta(k+1,j)
		x1=(long(nlong-1)+long(nlong))/2d0
		x2=(long(k)+long(k+1))/2d0
		x3=(long(k+1)+long(k+2))/2d0
	else if(k.eq.nlong-1) then
		y1=beta(k-1,j)
		y2=beta(k,j)
		y3=beta(1,j)
		x1=(long(k-1)+long(k))/2d0
		x2=(long(k)+long(k+1))/2d0
		x3=(long(1)+long(2))/2d0
	else
		y1=beta(k-1,j)
		y2=beta(k,j)
		y3=beta(k+1,j)
		x1=(long(k-1)+long(k))/2d0
		x2=(long(k)+long(k+1))/2d0
		x3=(long(k+1)+long(k+2))/2d0
	endif
	shift=-(((y1-y2)/(x1-x2))*((x1*x1-x2*x2)/(x1-x2)-(x2*x2-x3*x3)/(x2-x3))/((y1-y2)/(x1-x2)-(y2-y3)/(x2-x3))
     &		-(x1*x1-x2*x2)/(x1-x2))/2d0
	
	shift=shift*180d0/pi
	
	return
	end


	subroutine DiffuseBeta(beta,long,latt,nlong,nlatt,Kxx,Kyy,vxx,powvxx,contrast)
	IMPLICIT NONE
	integer nlatt,nlong,NN
	integer,allocatable :: IWORK(:)
	integer i,info,j,NRHS,ii(nlong-1,nlatt-1)
	real*8 beta(nlong,nlatt),Kxx,vxx,long(nlong),latt(nlatt),S(nlong-1,nlatt-1)
	real*8 pi,tot,x((nlatt-1)*(nlong-1)),contrast,eps,Kyy,powvxx
	parameter(eps=1d-4)
	real*8 la(nlatt-1),lo(nlong-1),tp,tm,tot1,tot2,cmax,cmin
	real*8,allocatable :: A(:,:)
	integer j0,jm,jp,k,niter,maxiter
	parameter(pi=3.14159265358979323846264338328d0)
	real*8 smax,smin,scale,betamin,betamax,contr

	do i=1,nlong-1
		lo(i)=(long(i)+long(i+1))/2d0
	enddo
	do j=1,nlatt-1
		la(j)=(latt(j)+latt(j+1))/2d0
	enddo
	k=0
	do i=1,nlong-1
		do j=1,nlatt-1
			k=k+1
			ii(i,j)=k
		enddo
	enddo
	if(contrast.le.0d0) then
		beta=0d0
		do i=1,nlong-1
			do j=1,nlatt-1
				if(abs(lo(i)).le.pi/2d0) beta(i,j)=abs(cos(lo(i))*cos(la(j)))
			enddo
		enddo
		return
	else if(contrast.ge.1d0.and.vxx.eq.0d0) then
		do i=1,nlong-1
			do j=1,nlatt-1
				beta(i,j)=0.25
			enddo
		enddo
		return
	endif

	allocate(IWORK(nlatt*nlong))
	allocate(A((nlatt-1)*(nlong-1),(nlatt-1)*(nlong-1)))

	smax=-sqrt(2.4)
	smin=-sqrt(3.2)
	scale=1d0
	
	contr=contrast*10d0
	maxiter=500
	niter=0
	do while(abs((contrast-contr)/(contrast+contr)).gt.eps.and.abs((smax-smin)/(smax+smin)).gt.eps.and.niter.lt.maxiter)
	if(nlong*nlatt.gt.1600) print*,((contrast-contr)/(contrast+contr))/eps,contr,scale
	
	niter=niter+1
	
	S=0d0
	A=0d0
	tot=0d0
	do i=1,nlong-1
		do j=1,nlatt-1
			if(abs(lo(i)).le.pi/2d0) S(i,j)=-cos(lo(i))*cos(la(j))
		enddo
	enddo
	do i=1,nlong-1
		do j=1,nlatt-1
			j0=ii(i,j)
			x(j0)=S(i,j)
			A(j0,j0)=A(j0,j0)-1d0
			if(i.gt.1) then
				jm=ii(i-1,j)
			else
				jm=ii(nlong-1,j)
			endif
			if(i.lt.nlong-1) then
				jp=ii(i+1,j)
			else
				jp=ii(1,j)
			endif
			A(j0,jm)=A(j0,jm)+0.5d0*vxx*scale*real(nlong)*cos(la(j))**(powvxx-1d0)
			A(j0,jp)=A(j0,jp)-0.5d0*vxx*scale*real(nlong)*cos(la(j))**(powvxx-1d0)
			A(j0,jm)=A(j0,jm)+Kxx*scale*(real(nlong)/cos(la(j)))**2
			A(j0,j0)=A(j0,j0)-2d0*Kxx*scale*(real(nlong)/cos(la(j)))**2
			A(j0,jp)=A(j0,jp)+Kxx*scale*(real(nlong)/cos(la(j)))**2

			if(j.gt.1) then
				jm=ii(i,j-1)
				tm=latt(j)
			else
				k=i+nlong/2
				if(k.gt.nlong-1) k=k-nlong+1
				jm=ii(k,1)
				tm=latt(1)
			endif
			if(j.lt.nlatt-1) then
				jp=ii(i,j+1)
				tp=latt(j+1)
			else
				k=i+nlong/2
				if(k.gt.nlong-1) k=k-nlong+1
				jp=ii(k,nlatt-1)
				tp=latt(nlatt)
			endif
			A(j0,jp)=A(j0,jp)+Kxx*Kyy*scale*real(nlatt*2)**2*cos(tp)/cos(la(j))
			A(j0,jm)=A(j0,jm)+Kxx*Kyy*scale*real(nlatt*2)**2*cos(tm)/cos(la(j))
			A(j0,j0)=A(j0,j0)-Kxx*Kyy*scale*real(nlatt*2)**2*(cos(tp)+cos(tm))/cos(la(j))
		enddo
	enddo

	NRHS=1
	NN=(nlong-1)*(nlatt-1)
	call DGESV( NN, NRHS, A, NN, IWORK, x, NN, info )

	betamin=0d0
	betamax=0d0
	tot1=0d0
	tot2=0d0
	do i=1,nlong-1
		do j=1,nlatt-1
			j0=ii(i,j)
			beta(i,j)=x(j0)
			tot1=tot1+beta(i,j)*cos(la(j))
			tot2=tot2+cos(la(j))
			if(abs(lo(i)).le.pi/2d0) then
				betamax=betamax+beta(i,j)*abs(cos(la(j))*cos(lo(i)))
			else
				betamin=betamin+beta(i,j)*abs(cos(la(j))*cos(lo(i)))
			endif
		enddo
	enddo
	beta=beta*0.25*tot2/tot1
	contr=betamin/betamax
	if(contr.lt.contrast) then
		smin=scale
		cmin=contr
		if(smax.lt.0d0) then
			scale=scale*10d0
		else
			scale=scale+(smax-scale)*(contrast-contr)/(cmax-contr)
		endif
	else
		smax=scale
		cmax=contr
		if(smin.lt.0d0) then
			scale=scale/10d0
		else
			scale=scale+(smin-scale)*(contrast-contr)/(cmin-contr)
		endif
	endif
	enddo

	deallocate(IWORK)
	deallocate(A)
	
	return
	end
	


	subroutine DiffuseBetaRotate(beta,long,latt,nlong,nlatt,Kxx,Kyy,vxx,powvxx,contrast)
	IMPLICIT NONE
	integer nlatt,nlong,NN
	integer,allocatable :: IWORK(:)
	integer i,info,j,NRHS,ii(nlong-1,nlatt-1)
	real*8 beta(nlong,nlatt),Kxx,vxx,long(nlong),latt(nlatt),S(nlong-1,nlatt-1)
	real*8 pi,tot,x((nlatt-1)*(nlong-1)),contrast,eps,Kyy,powvxx
	parameter(eps=1d-4)
	real*8 la(nlatt-1),lo(nlong-1),tp,tm,tot1,tot2,cmax,cmin
	real*8,allocatable :: A(:,:)
	integer j0,jm,jp,k,niter,maxiter
	parameter(pi=3.14159265358979323846264338328d0)
	real*8 smax,smin,scale,betamin,betamax,contr

	do i=1,nlong-1
		lo(i)=(long(i)+long(i+1))/2d0
	enddo
	do j=1,nlatt-1
		la(j)=(latt(j)+latt(j+1))/2d0
	enddo
	k=0
	do i=1,nlong-1
		do j=1,nlatt-1
			k=k+1
			ii(i,j)=k
		enddo
	enddo
	if(contrast.le.0d0) then
		beta=0d0
		do i=1,nlong-1
			do j=1,nlatt-1
				beta(i,j)=abs(cos(la(j)))/pi
			enddo
		enddo
		return
	else if(contrast.ge.1d0.and.vxx.eq.0d0) then
		do i=1,nlong-1
			do j=1,nlatt-1
				beta(i,j)=0.25
			enddo
		enddo
		return
	endif

	allocate(IWORK(nlatt*nlong))
	allocate(A((nlatt-1)*(nlong-1),(nlatt-1)*(nlong-1)))

	smax=-sqrt(2.4)
	smin=-sqrt(3.2)
	scale=1d0
	
	contr=contrast*10d0
	maxiter=500
	niter=0
	do while(abs((contrast-contr)/(contrast+contr)).gt.eps.and.abs((smax-smin)/(smax+smin)).gt.eps.and.niter.lt.maxiter)
	if(nlong*nlatt.gt.1600) print*,((contrast-contr)/(contrast+contr))/eps,contr,scale
	
	niter=niter+1
	
	S=0d0
	A=0d0
	tot=0d0
	do i=1,nlong-1
		do j=1,nlatt-1
			S(i,j)=-cos(la(j))/pi
		enddo
	enddo
	do i=1,nlong-1
		do j=1,nlatt-1
			j0=ii(i,j)
			x(j0)=S(i,j)
			A(j0,j0)=A(j0,j0)-1d0
			if(i.gt.1) then
				jm=ii(i-1,j)
			else
				jm=ii(nlong-1,j)
			endif
			if(i.lt.nlong-1) then
				jp=ii(i+1,j)
			else
				jp=ii(1,j)
			endif
			A(j0,jm)=A(j0,jm)+0.5d0*vxx*scale*real(nlong)*cos(la(j))**(powvxx-1d0)
			A(j0,jp)=A(j0,jp)-0.5d0*vxx*scale*real(nlong)*cos(la(j))**(powvxx-1d0)
			A(j0,jm)=A(j0,jm)+Kxx*scale*(real(nlong)/cos(la(j)))**2
			A(j0,j0)=A(j0,j0)-2d0*Kxx*scale*(real(nlong)/cos(la(j)))**2
			A(j0,jp)=A(j0,jp)+Kxx*scale*(real(nlong)/cos(la(j)))**2

			if(j.gt.1) then
				jm=ii(i,j-1)
				tm=latt(j)
			else
				k=i+nlong/2
				if(k.gt.nlong-1) k=k-nlong+1
				jm=ii(k,1)
				tm=latt(1)
			endif
			if(j.lt.nlatt-1) then
				jp=ii(i,j+1)
				tp=latt(j+1)
			else
				k=i+nlong/2
				if(k.gt.nlong-1) k=k-nlong+1
				jp=ii(k,nlatt-1)
				tp=latt(nlatt)
			endif
			A(j0,jp)=A(j0,jp)+Kxx*Kyy*scale*real(nlatt*2)**2*cos(tp)/cos(la(j))
			A(j0,jm)=A(j0,jm)+Kxx*Kyy*scale*real(nlatt*2)**2*cos(tm)/cos(la(j))
			A(j0,j0)=A(j0,j0)-Kxx*Kyy*scale*real(nlatt*2)**2*(cos(tp)+cos(tm))/cos(la(j))
		enddo
	enddo

	NRHS=1
	NN=(nlong-1)*(nlatt-1)
	call DGESV( NN, NRHS, A, NN, IWORK, x, NN, info )

	betamin=1d0
	betamax=0d0
	tot1=0d0
	tot2=0d0
	do i=1,nlong-1
		do j=1,nlatt-1
			j0=ii(i,j)
			beta(i,j)=x(j0)
			tot1=tot1+beta(i,j)*cos(la(j))
			tot2=tot2+cos(la(j))
			if(beta(i,j).lt.betamin) betamin=beta(i,j)
			if(beta(i,j).gt.betamax) betamax=beta(i,j)
		enddo
	enddo
	beta=beta*0.25*tot2/tot1
	contr=betamin/betamax
	if(contr.lt.contrast) then
		smin=scale
		cmin=contr
		if(smax.lt.0d0) then
			scale=scale*10d0
		else
			scale=scale+(smax-scale)*(contrast-contr)/(cmax-contr)
		endif
	else
		smax=scale
		cmax=contr
		if(smin.lt.0d0) then
			scale=scale/10d0
		else
			scale=scale+(smin-scale)*(contrast-contr)/(cmin-contr)
		endif
	endif
	enddo

	deallocate(IWORK)
	deallocate(A)
	
	return
	end
	


	subroutine DiffuseBeta_new(beta,long,latt,nlong,nlatt,Kxx,Kyy,vxx,contrast)
	IMPLICIT NONE
	integer nlatt,nlong,NN
	integer,allocatable :: IWORK(:)
	integer i,info,j,NRHS,ii(nlong-1,nlatt-1)
	real*8 beta(nlong,nlatt),Kxx,vxx,long(nlong),latt(nlatt),S(nlong-1,nlatt-1)
	real*8 pi,tot,x((nlatt-1)*(nlong-1)),contrast,eps,Kyy
	parameter(eps=1d-4)
	real*8 la(nlatt-1),lo(nlong-1),tp,tm,tot1,tot2,cmax,cmin
	real*8,allocatable :: A(:,:)
	integer j0,jm,jp,k,niter,maxiter
	parameter(pi=3.14159265358979323846264338328d0)
	real*8 smax,smin,scale,betamin,betamax,contr
	real*8 f1m,f1p,f2m,f2p,am,ap,a0,bm,bp,b0

	allocate(IWORK(nlatt*nlong))
	allocate(A((nlatt-1)*(nlong-1),(nlatt-1)*(nlong-1)))

	smax=-sqrt(2.4)
	smin=-sqrt(3.2)
	scale=1d0
	
	contr=contrast*10d0
	maxiter=500
	niter=0
	do while(abs((contrast-contr)/(contrast+contr)).gt.eps.and.abs((smax-smin)/(smax+smin)).gt.eps.and.niter.lt.maxiter)
	if(nlong*nlatt.gt.1600) print*,((contrast-contr)/(contrast+contr))/eps,contr,scale
	
	niter=niter+1
	
	do i=1,nlong-1
		lo(i)=(long(i)+long(i+1))/2d0
	enddo
	do j=1,nlatt-1
		la(j)=(latt(j)+latt(j+1))/2d0
	enddo
	k=0
	do i=1,nlong-1
		do j=1,nlatt-1
			k=k+1
			ii(i,j)=k
		enddo
	enddo
	S=0d0
	A=0d0
	tot=0d0
	do i=1,nlong-1
		do j=1,nlatt-1
			if(abs(lo(i)).le.pi/2d0) S(i,j)=-cos(lo(i))*cos(la(j))
		enddo
	enddo
	do i=1,nlong-1
		do j=1,nlatt-1
			j0=ii(i,j)
			x(j0)=S(i,j)
			A(j0,j0)=A(j0,j0)-1d0
			if(i.gt.1) then
				jm=ii(i-1,j)
			else
				jm=ii(nlong-1,j)
			endif
			if(i.lt.nlong-2) then
				jp=ii(i+1,j)
			else
				jp=ii(1,j)
			endif
			A(j0,jm)=A(j0,jm)+0.5d0*vxx*scale*real(nlong)/cos(la(j))
			A(j0,jp)=A(j0,jp)-0.5d0*vxx*scale*real(nlong)/cos(la(j))
			A(j0,jm)=A(j0,jm)+Kxx*scale*(real(nlong)/cos(la(j)))**2
			A(j0,j0)=A(j0,j0)-2d0*Kxx*scale*(real(nlong)/cos(la(j)))**2
			A(j0,jp)=A(j0,jp)+Kxx*scale*(real(nlong)/cos(la(j)))**2

			if(j.gt.1) then
				jm=ii(i,j-1)
				tm=la(j-1)
			else
				k=i+nlong/2
				if(k.gt.nlong-1) k=k-nlong+1
				jm=ii(k,1)
				tm=-la(1)
			endif
			if(j.lt.nlatt-1) then
				jp=ii(i,j+1)
				tp=la(j+1)
			else
				k=i+nlong/2
				if(k.gt.nlong-1) k=k-nlong+1
				jp=ii(k,nlatt-1)
				tp=-la(nlatt-1)
			endif
		f1m=tm-la(j)
		f1p=tp-la(j)
		f2m=tm**2-la(j)**2
		f2p=tp**2-la(j)**2

		am=1d0/((f2m/f1m-f2p/f1p)*f1m)
		ap=1d0/((f2p/f1p-f2m/f1m)*f1p)
		a0=(1d0/(f2m/f1m-f2p/f1p))*(1d0/f1p-1d0/f1m)

		bm=1d0/((f1m/f2m-f1p/f2p)*f2m)
		bp=1d0/((f1p/f2p-f1m/f2m)*f2p)
		b0=(1d0/(f1m/f2m-f1p/f2p))*(1d0/f2p-1d0/f2m)

		A(j0,jm)=A(j0,jm)-Kxx*(2d0*am-(2d0*am*la(j)+bm)*tan(la(j)))/(2d0*pi)
		A(j0,jp)=A(j0,jp)-Kxx*(2d0*ap-(2d0*ap*la(j)+bp)*tan(la(j)))/(2d0*pi)
		A(j0,j0)=A(j0,j0)-Kxx*(2d0*a0-(2d0*a0*la(j)+b0)*tan(la(j)))/(2d0*pi)

		enddo
	enddo

	NRHS=1
	NN=(nlong-1)*(nlatt-1)
	call DGESV( NN, NRHS, A, NN, IWORK, x, NN, info )

	betamin=0d0
	betamax=0d0
	tot1=0d0
	tot2=0d0
	do i=1,nlong-1
		do j=1,nlatt-1
			j0=ii(i,j)
			beta(i,j)=x(j0)
			tot1=tot1+beta(i,j)*cos(la(j))
			tot2=tot2+cos(la(j))
			if(abs(lo(i)).le.pi/2d0) then
				betamax=betamax+beta(i,j)*abs(cos(la(j))*cos(lo(i)))
			else
				betamin=betamin+beta(i,j)*abs(cos(la(j))*cos(lo(i)))
			endif
		enddo
	enddo
	beta=beta*0.25*tot2/tot1
	contr=betamin/betamax
	if(contr.lt.contrast) then
		smin=scale
		cmin=contr
		if(smax.lt.0d0) then
			scale=scale*10d0
		else
			scale=scale+(smax-scale)*(contrast-contr)/(cmax-contr)
		endif
	else
		smax=scale
		cmax=contr
		if(smin.lt.0d0) then
			scale=scale/10d0
		else
			scale=scale+(smin-scale)*(contrast-contr)/(cmin-contr)
		endif
	endif
	enddo

	deallocate(IWORK)
	deallocate(A)
	
	return
	end
	



	subroutine DiffuseBetaFix(beta,long,latt,nlong,nlatt,Kxx,vxx)
	IMPLICIT NONE
	integer nlatt,nlong,NN
	integer,allocatable :: IWORK(:)
	integer i,info,j,NRHS,ii(nlong-1,nlatt-1)
	real*8 beta(nlong,nlatt),Kxx,vxx,long(nlong),latt(nlatt),S(nlong-1,nlatt-1)
	real*8 pi,tot,x((nlatt-1)*(nlong-1))
	real*8 la(nlatt-1),lo(nlong-1),tp,tm,tot1,tot2
	real*8,allocatable :: A(:,:)
	integer j0,jm,jp,k
	parameter(pi=3.14159265358979323846264338328d0)
	allocate(IWORK(50*nlatt*nlong*nlatt*nlong))
	allocate(A((nlatt-1)*(nlong-1),(nlatt-1)*(nlong-1)))

	do i=1,nlong-1
		lo(i)=(long(i)+long(i+1))/2d0
	enddo
	do j=1,nlatt-1
		la(j)=(latt(j)+latt(j+1))/2d0
	enddo
	k=0
	do i=1,nlong-1
		do j=1,nlatt-1
			k=k+1
			ii(i,j)=k
		enddo
	enddo
	S=0d0
	A=0d0
	tot=0d0
	do i=1,nlong-1
		do j=1,nlatt-1
			if(abs(lo(i)).le.pi/2d0) S(i,j)=-cos(lo(i))*cos(la(j))
		enddo
	enddo
	do i=1,nlong-1
		do j=1,nlatt-1
			j0=ii(i,j)
			x(j0)=S(i,j)
			A(j0,j0)=A(j0,j0)-1d0
			if(i.gt.1) then
				jm=ii(i-1,j)
			else
				jm=ii(nlong-1,j)
			endif
			if(i.lt.nlong-2) then
				jp=ii(i+1,j)
			else
				jp=ii(1,j)
			endif
			A(j0,jm)=A(j0,jm)+0.5d0*vxx*real(nlong)/cos(la(j))
			A(j0,jp)=A(j0,jp)-0.5d0*vxx*real(nlong)/cos(la(j))
			A(j0,jm)=A(j0,jm)+Kxx*(real(nlong)/cos(la(j)))**2
			A(j0,j0)=A(j0,j0)-2d0*Kxx*(real(nlong)/cos(la(j)))**2
			A(j0,jp)=A(j0,jp)+Kxx*(real(nlong)/cos(la(j)))**2

			if(j.gt.1) then
				jm=ii(i,j-1)
				tm=latt(j-1)
			else
				k=i+nlong/2
				if(k.gt.nlong-1) k=k-nlong+1
				jm=ii(k,1)
				tm=latt(1)
			endif
			if(j.lt.nlatt-2) then
				jp=ii(i,j+1)
				tp=latt(j+2)
			else
				k=i+nlong/2
				if(k.gt.nlong-1) k=k-nlong+1
				jp=ii(k,nlatt-1)
				tp=latt(nlatt)
			endif
			A(j0,jp)=A(j0,jp)+Kxx*real(nlatt*2)**2*cos(tp)/cos(la(j))
			A(j0,jm)=A(j0,jm)+Kxx*real(nlatt*2)**2*cos(tm)/cos(la(j))
			A(j0,j0)=A(j0,j0)-Kxx*real(nlatt*2)**2*(cos(tp)+cos(tm))/cos(la(j))
		enddo
	enddo

	NRHS=1
	NN=(nlong-1)*(nlatt-1)
	call DGESV( NN, NRHS, A, NN, IWORK, x, NN, info )

	do i=1,nlong-1
		do j=1,nlatt-1
			j0=ii(i,j)
			beta(i,j)=x(j0)
			tot1=tot1+beta(i,j)*cos(la(j))
			tot2=tot2+cos(la(j))
		enddo
	enddo
	beta=beta*0.25*tot2/tot1

	deallocate(IWORK)
	deallocate(A)
	
	return
	end
	



	subroutine DiffuseBeta1D(beta,theta,nbeta,K,v)
	IMPLICIT NONE
	integer nbeta
	integer i,info,IWORK(50*nbeta*nbeta),j,NRHS
	real*8 beta(nbeta),K,v,theta(nbeta),S(nbeta)
	real*8 pi,A(nbeta,nbeta),tot
	parameter(pi=3.14159265358979323846264338328d0)
	
	S=0d0
	A=0d0
	do i=1,nbeta
		if(abs(theta(i)).le.pi/2d0) S(i)=-cos(theta(i))
	enddo
	tot=-sum(S)
	do i=2,nbeta-1
		j=i
		A(j,i)=-1d0
		A(j,i-1)=A(j,i-1)+0.5d0*v*real(nbeta)
		A(j,i+1)=A(j,i+1)-0.5d0*v*real(nbeta)
		A(j,i-1)=A(j,i-1)+K*real(nbeta)**2
		A(j,i  )=A(j,i  )-2d0*K*real(nbeta)**2
		A(j,i+1)=A(j,i+1)+K*real(nbeta)**2
	enddo
	A(1,1)=-1d0
	A(1,nbeta)=A(1,nbeta)+0.5d0*v*real(nbeta)
	A(1,2)=A(1,2)-0.5d0*v*real(nbeta)
	A(1,nbeta)=A(1,nbeta)+K*real(nbeta)**2
	A(1,1)=A(1,1)-2d0*K*real(nbeta)**2
	A(1,2)=A(1,2)+K*real(nbeta)**2

	A(nbeta,1:nbeta)=1d0
	S(nbeta)=tot

	beta=S

	NRHS=1
	call DGESV( nbeta, NRHS, A, nbeta, IWORK, beta, nbeta, info )
	
	return
	end
	



	
	logical function HitR3D(Rad,r,a,b,v)
	IMPLICIT NONE
	real*8 Rad,r,a,b,cc,discr,vr1,vr2,v,q
	
	HitR3D=.false.
	v=1d200

	cc=r-Rad
	discr=(b**2-4d0*cc*a)
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(b.gt.0d0) then
			q=-0.5d0*(b+discr)
		else
			q=-0.5d0*(b-discr)
		endif
		vr1=q/a
		vr2=cc/q
		if(vr1.gt.0d0) then
			v=vr1
			HitR3D=.true.
		endif
		if(vr2.gt.0d0.and.vr2.lt.v) then
			v=vr2
			HitR3D=.true.
		endif
	endif
	return
	end


	
	logical function hitRin(Rad,r,a,b,v)
	IMPLICIT NONE
	real*8 Rad,r,a,b,cc,discr,vr1,vr2,v,q
	
	hitRin=.false.
	v=1d200

	cc=r-Rad
	discr=(b**2-4d0*cc*a)
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(b.gt.0d0) then
			q=-0.5d0*(b+discr)
		else
			q=-0.5d0*(b-discr)
		endif
		vr1=q/a
		vr2=cc/q
		v=vr1
		hitRin=.true.
		if(vr2.gt.v) then
			v=vr2
			hitRin=.true.
		endif
	endif
	if(.not.v.gt.0d0) hitRin=.false.

	return
	end


	
	logical function hitRout(Rad,r,a,b,v)
	IMPLICIT NONE
	real*8 Rad,r,a,b,cc,discr,vr1,vr2,v,q
	
	hitRout=.false.
	v=1d200

	cc=r-Rad
	discr=(b**2-4d0*cc*a)
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(b.gt.0d0) then
			q=-0.5d0*(b+discr)
		else
			q=-0.5d0*(b-discr)
		endif
		vr1=q/a
		vr2=cc/q
		v=vr1
		hitRout=.true.
		if(vr2.lt.v) then
			v=vr2
			hitRout=.true.
		endif
	endif
	if(.not.v.gt.0d0) hitRout=.false.
	return
	end


	logical function hitT(z,vz,Thet,r,a,b,v,midplane)
	IMPLICIT NONE
	real*8 Thet,r,a,b,at,bt,ct,discr,vt1,vt2,v,q,z,vz
	logical midplane

	hitT=.false.
	v=1d200

	if(midplane) then
		v=-z/vz
		if(v.gt.0d0) hitT=.true.
		return
	endif

	at=Thet*a-vz*vz
	bt=Thet*b-2d0*z*vz
	ct=Thet*r-z*z
	discr=bt*bt-4d0*at*ct
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(bt.gt.0d0) then
			q=-0.5d0*(bt+discr)
		else
			q=-0.5d0*(bt-discr)
		endif
		vt1=q/at
		vt2=ct/q
		if(vt1.gt.0d0) then
			v=vt1
			hitT=.true.
		endif
		if(vt2.gt.0d0.and.vt2.lt.v) then
			v=vt2
			hitT=.true.
		endif
	endif
	
	if(z*(z+v*vz).lt.0d0.and..not.midplane) hitT=.false.

	return
	end

	logical function hitTsame(z,vz,Thet,r,a,b,v)
	IMPLICIT NONE
	real*8 Thet,r,a,b,at,bt,v,z,vz

	hitTsame=.true.
	v=1d200

	bt=Thet*b-2d0*z*vz
	at=Thet*a-vz**2
	v=-bt/at
	if(.not.v.gt.0d0) hitTsame=.false.
	if(z*(z+v*vz).lt.0d0) hitTsame=.false.

	return
	end

	logical function hitP(tanx,tany,x0,vx,y0,vy,v)
	IMPLICIT NONE
	real*8 tanx,tany,x0,vx,y0,vy,v
	
	hitP=.true.
	v=(-tany*y0-tanx*x0)/(tanx*vx+tany*vy)
	
	if(.not.v.gt.0d0) hitP=.false.
	
	return
	end

	



	subroutine TravelSph(x,y,z,vx,vy,vz,edgeNR,i1,i2,i3,v,i1next,i2next,i3next,edgenext,n1,n2,n3)
	use Struct3D
	use Constants
	IMPLICIT NONE
	real*8 x,y,z,vx,vy,vz,v
	integer edgeNR,i1,i2,i3,i1next,i2next,i3next,edgenext,n1,n2,n3

	real*8 a,b,r,R1,R2,T1,T2,vR1,vR2,vT1,vT2,vP1,vP2
	logical hitR1,hitR2,hitR3D,hitT1,hitT2,hitT,hitTsame
	logical hitP1,hitP2,hitP,i1midplane,i2midplane
	real*8 tanx1,tanx2,tany1,tany2

	r=x**2+y**2+z**2
	R1=R3D2(ibeta(i2,i3),i1)
	R2=R3D2(ibeta(i2,i3),i1+1)
	T1=cost2(i3)
	T2=cost2(i3+1)
	tanx1=tanx(i2)
	tanx2=tanx(i2+1)
	tany1=tany(i2)
	tany2=tany(i2+1)

	b=2d0*(x*vx+y*vy+z*vz)
	a=vx**2+vy**2+vz**2

	i1midplane=(i3.eq.(nlatt+1)/2)
	i2midplane=((i3+1).eq.(nlatt+1)/2)

	hitT1=.false.
	hitT2=.false.
	select case(edgeNr)
		case(1)
			hitR1=.false.
			vR1=1d200
			hitR2=HitR3D(R2,r,a,b,vR2)
			if(i3.gt.1) hitT1=hitT(z,vz,T1,r,a,b,vT1,i1midplane)
			if(i3.lt.n3-1) hitT2=hitT(z,vz,T2,r,a,b,vT2,i2midplane)
			hitP1=hitP(tanx1,tany1,x,vx,y,vy,vP1)
			hitP2=hitP(tanx2,tany2,x,vx,y,vy,vP2)
		case(2)
			hitR1=HitR3D(R1,r,a,b,vR1)
			hitR2=.true.
			vR2=-b/a
			if(i3.gt.1) hitT1=hitT(z,vz,T1,r,a,b,vT1,i1midplane)
			if(i3.lt.n3-1) hitT2=hitT(z,vz,T2,r,a,b,vT2,i2midplane)
			hitP1=hitP(tanx1,tany1,x,vx,y,vy,vP1)
			hitP2=hitP(tanx2,tany2,x,vx,y,vy,vP2)
		case(3)
			hitR1=HitR3D(R1,r,a,b,vR1)
			hitR2=HitR3D(R2,r,a,b,vR2)
			if(latt(i3).le.pi/2d0.or.i1midplane) then
				hitT1=.false.
				vT1=1d200
			else if(i3.gt.1) then
				hitT1=hitTsame(z,vz,T1,r,a,b,vT1)
			endif
			if(i3.lt.n3-1) hitT2=hitT(z,vz,T2,r,a,b,vT2,i2midplane)
			hitP1=hitP(tanx1,tany1,x,vx,y,vy,vP1)
			hitP2=hitP(tanx2,tany2,x,vx,y,vy,vP2)
		case(4)
			hitR1=HitR3D(R1,r,a,b,vR1)
			hitR2=HitR3D(R2,r,a,b,vR2)
			if(i3.gt.1) hitT1=hitT(z,vz,T1,r,a,b,vT1,i1midplane)
			if(latt(i3+1).ge.pi/2d0.or.i2midplane) then
				hitT2=.false.
				vT2=1d200
			else if(i3.lt.n3-1) then
				hitT2=hitTsame(z,vz,T2,r,a,b,vT2)
			endif
			hitP1=hitP(tanx1,tany1,x,vx,y,vy,vP1)
			hitP2=hitP(tanx2,tany2,x,vx,y,vy,vP2)
		case(5)
			hitR1=HitR3D(R1,r,a,b,vR1)
			hitR2=HitR3D(R2,r,a,b,vR2)
			if(i3.gt.1) hitT1=hitT(z,vz,T1,r,a,b,vT1,i1midplane)
			if(i3.lt.n3-1) hitT2=hitT(z,vz,T2,r,a,b,vT2,i2midplane)
			hitP1=.false.
			vP1=1d200
			hitP2=hitP(tanx2,tany2,x,vx,y,vy,vP2)
		case(6)
			hitR1=HitR3D(R1,r,a,b,vR1)
			hitR2=HitR3D(R2,r,a,b,vR2)
			if(i3.gt.1) hitT1=hitT(z,vz,T1,r,a,b,vT1,i1midplane)
			if(i3.lt.n3-1) hitT2=hitT(z,vz,T2,r,a,b,vT2,i2midplane)
			hitP1=hitP(tanx1,tany1,x,vx,y,vy,vP1)
			hitP2=.false.
			vP2=1d200
		case default
			hitR1=HitR3D(R1,r,a,b,vR1)
			hitR2=HitR3D(R2,r,a,b,vR2)
			if(i3.gt.1) hitT1=hitT(z,vz,T1,r,a,b,vT1,i1midplane)
			if(i3.lt.n3-1) hitT2=hitT(z,vz,T2,r,a,b,vT2,i2midplane)
			hitP1=hitP(tanx1,tany1,x,vx,y,vy,vP1)
			hitP2=hitP(tanx2,tany2,x,vx,y,vy,vP2)
	end select

c	hitT1=.false.
c	hitT2=.false.
c	hitP1=.false.
c	hitP2=.false.

	v=1d200
	edgenext=1
	i1next=0
	i2next=0
	i3next=0
	if(hitR1.and.vR1.lt.v.and.vR1.gt.0d0) then
		v=vR1
		i1next=i1-1
		i2next=i2
		i3next=i3
		edgenext=2
	endif
	if(hitR2.and.vR2.lt.v.and.vR2.gt.0d0) then
		v=vR2
		i1next=i1+1
		i2next=i2
		i3next=i3
		edgenext=1
	endif
	if(hitT1.and.vT1.lt.v.and.vT1.gt.0d0) then
		v=vT1
		i1next=i1
		i2next=i2
		i3next=i3-1
		edgenext=4
		if(i3next.lt.1) then
			i3next=1
			edgenext=3
		endif
	endif
	if(hitT2.and.vT2.lt.v.and.vT2.gt.0d0) then
		v=vT2
		i1next=i1
		i2next=i2
		i3next=i3+1
		edgenext=3
		if(i3next.ge.n3) then
			i3next=n3-1
			edgenext=4
		endif
	endif
	if(hitP1.and.vP1.lt.v.and.vP1.gt.0d0) then
		v=vP1
		i1next=i1
		i2next=i2-1
		i3next=i3
		if(i2next.lt.1) i2next=n2-1
		edgenext=6
	endif
	if(hitP2.and.vP2.lt.v.and.vP2.gt.0d0) then
		v=vP2
		i1next=i1
		i2next=i2+1
		i3next=i3
		if(i2next.ge.n2) i2next=1
		edgenext=5
	endif

	return
	end
	


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotateZ3D(x,y,z,theta)
	IMPLICIT NONE
	real*8 x,y,z,xx,yy,r,theta
	real*8 cost,sint
	cost=cos(theta)
	sint=sin(theta)

	xx=x*cost-y*sint
	yy=y*cost+x*sint
	x=xx
	y=yy
	

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotateY3D(x,y,z,theta)
	IMPLICIT NONE
	real*8 x,y,z,xx,zz,r,theta
	real*8 cost,sint
	cost=cos(theta)
	sint=sin(theta)

	xx=x*cost-z*sint
	zz=z*cost+x*sint
	x=xx
	z=zz

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------



	subroutine ComputeScatter(BBr,Si,Ca,Cs)
	use GlobalSetup
	use Constants
	use Struct3D
	IMPLICIT NONE
	integer inu,nnu,ilam,ir,ig,inu0,iter,niter,info,NRHS
	parameter(nnu=10,niter=500)
	real*8 tau,d,tauR_nu(0:nr,nlam,ng),contr,Jstar_nu(nr,nlam,ng),Pb(nr+1)
	real*8 Si(nlam,ng,0:nr,nnu0),BBr(nlam,0:nr),Ca(nlam,ng,nr),Cs(nlam,nr),Ce(nlam,ng,nr)
	real*8 nu(nnu),wnu(nnu),must,tauRs(nr),Ijs(nr),eps,Planck,tot,wabs(nlam,ng,nr),wscat(nlam,ng,nr)
	logical err
	parameter(eps=1d-20)
	real*8,allocatable :: tauR(:),Ij(:),Ih(:)
	
	do ir=1,nr
		do ig=1,ng
			do ilam=1,nlam
				Ce(ilam,ig,ir)=Ca(ilam,ig,ir)+Cs(ilam,ir)
				wabs(ilam,ig,ir)=Ca(ilam,ig,ir)/Ce(ilam,ig,ir)
				if(.not.wabs(ilam,ig,ir).gt.1d-4) then
					wabs(ilam,ig,ir)=1d-4
					Ca(ilam,ig,ir)=Cs(ilam,ir)/(1d0/wabs(ilam,ig,ir)-1d0)
					Ce(ilam,ig,ir)=Ca(ilam,ig,ir)+Cs(ilam,ir)
				endif
				wscat(ilam,ig,ir)=Cs(ilam,ir)/Ce(ilam,ig,ir)
				if(.not.wscat(ilam,ig,ir).gt.0d0) then
					wscat(ilam,ig,ir)=0d0
					Cs(ilam,ir)=0d0
					Ce(ilam,ig,ir)=Ca(ilam,ig,ir)+Cs(ilam,ir)
				endif
			enddo
		enddo
	enddo

	do inu0=1,nnu0
		do ig=1,ng
			Si(1:nlam,ig,1:nr,inu0)=BBr(1:nlam,1:nr)*wabs(1:nlam,ig,1:nr)
			Si(1:nlam,ig,0,inu0)=BBr(1:nlam,0)*surface_emis(1:nlam)
		enddo
	enddo
	if(.not.scattering) return

	Pb(1)=P(1)
	do ir=2,nr
		Pb(ir)=sqrt(P(ir-1)*P(ir))
	enddo
	Pb(nr+1)=0d0
	tauR_nu=0d0
	do ilam=1,nlam
		do ig=1,ng
			do ir=nr,1,-1
				if(ir.eq.nr) then
					d=(P(ir)-Pb(ir+1))*1d6/grav(ir)
					tau=d*Ce(ilam,ig,ir)/dens(ir)
				else
					d=(Pb(ir+1)-P(ir+1))*1d6/grav(ir+1)
					tau=d*Ce(ilam,ig,ir+1)/dens(ir+1)
					d=(P(ir)-Pb(ir+1))*1d6/grav(ir)
					tau=tau+d*Ce(ilam,ig,ir)/dens(ir)
				endif
				if(.not.tau.gt.0d0) then
					tau=0d0
				endif
				tau=tau+1d-10
				tau=1d0/(1d0/tau+1d-10)
				if(ir.eq.nr) then
					tauR_nu(ir,ilam,ig)=tau
				else
					tauR_nu(ir,ilam,ig)=tauR_nu(ir+1,ilam,ig)+tau
				endif
			enddo
		enddo
	enddo

	call gauleg(0d0,1d0,nu,wnu,nnu)


!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(tauR,Ij,ilam,ig,inu0,inu,contr,must,Ih)
!$OMP& SHARED(nr,ng,nlam,Fstar,Dplanet,Si,Ca,Ce,Cs,nu,wnu,surface_emis,tauR_nu,BBr,scattstar,lamemis,
!$OMP&        nnu0,wscat)
	allocate(tauR(nr),Ij(nr),Ih(nr))
!$OMP DO
	do ilam=1,nlam
		if(lamemis(ilam)) then
		do ig=1,ng
			do inu0=1,nnu0
				if(inu0.ne.nnu0.and.scattstar) then
					must=(real(inu0)-0.5)/real(nnu0-1)
					contr=(Fstar(ilam)/(pi*Dplanet**2))
					tauR(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(must)
					Si(ilam,ig,1:nr,inu0)=Si(ilam,ig,1:nr,inu0)+contr*exp(-tauR(1:nr))*wscat(ilam,ig,1:nr)/4d0
					contr=contr*exp(-tauR(1))*must
				else
					contr=0d0
				endif
				Si(ilam,ig,0,inu0)=Si(ilam,ig,0,inu0)+contr*(1d0-surface_emis(ilam))
				do inu=1,nnu
					tauR(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(nu(inu))
					tauR(1:nr)=abs(tauR(1:nr)-tauR(1))
					Ij(1:nr)=(BBr(ilam,0)*surface_emis(ilam)+contr*(1d0-surface_emis(ilam)))*exp(-tauR(1:nr))
					Si(ilam,ig,1:nr,inu0)=Si(ilam,ig,1:nr,inu0)+wnu(inu)*Ij(1:nr)*wscat(ilam,ig,1:nr)
				enddo
			enddo
			call AddScatter(Si(ilam,ig,1:nr,1:nnu0),tauR_nu(1:nr,ilam,ig),
     &	Ca(ilam,ig,1:nr),Cs(ilam,1:nr),Ce(ilam,ig,1:nr),(1d0-surface_emis(ilam)),nr,nu,wnu,nnu,nnu0)
			do inu0=1,nnu0
				do inu=1,nnu
					tauR(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(nu(inu))
					call SolveIjExp(tauR,Si(ilam,ig,1:nr,inu0),Ij,Ih,nr)
					Si(ilam,ig,0,inu0)=Si(ilam,ig,0,inu0)-nu(inu)*Ih(1)*wnu(inu)*(1d0-surface_emis(ilam))
				enddo
			enddo
		enddo
		endif
	enddo
!$OMP END DO
	deallocate(tauR,Ij,Ih)
!$OMP FLUSH
!$OMP END PARALLEL
	
	return
	end


	subroutine ComputeHapke(emis,e1,e2,nlam,computelam)
	IMPLICIT NONE
	integer i,nlam,nt,j,k,nmu
	parameter(nmu=10)
	real*8 e1(nlam),e2(nlam),emis(nlam),Rs,Rp,theta,tot,pi,w,g
	real*8 mu(nmu),wmu(nmu),refl,H1,H2,lmie,rmie,csmie,cemie
	logical computelam(nlam)

	call gauleg(0d0,1d0,mu,wmu,nmu)
	pi=3.1415926536
	nt=10
	tot=0d0
	do j=1,nt
		theta=real(j-1)*pi/(2d0*real(nt-1))
		tot=tot+sin(theta)
	enddo
	do i=1,nlam
		if(computelam(i)) then
		rmie=1d0
		lmie=rmie/100d0
		call callBHMIE(rmie,lmie,e1(i),e2(i),csmie,cemie)
		w=csmie/cemie
		if(.not.w.lt.(1d0-1d-4)) w=1d0-1d-4
		if(.not.w.gt.0d0) w=0d0

		g=sqrt(1d0-w)
		refl=0d0
		do j=1,nmu
			do k=1,nmu
				H1=(1d0+2d0*mu(j))/(1d0+2d0*g*mu(j))
				H2=(1d0+2d0*mu(k))/(1d0+2d0*g*mu(k))
				refl=refl+wmu(j)*wmu(k)*mu(j)*H1*H2/(mu(j)+mu(k))
			enddo
		enddo
		emis(i)=1d0-w*refl/2d0
c		emis(i)=1d0+log(2d0*g+1d0)*(g-1)/(2d0*g)
		else
		emis(i)=1d0
		endif
	enddo

	return
	end	
	


	subroutine Reflect(n,k,theta,Rs,Rp)
	IMPLICIT NONE
	real*8 theta,n,k,Rs,Rp
	complex*16 m
	
	m=dcmplx(n,k)
	
	Rs=cdabs((cos(theta)-sqrt(m**2-sin(theta)**2))/
     &		(cos(theta)+sqrt(m**2-sin(theta)**2)))**2
	Rp=cdabs((sqrt(1d0-(sin(theta)/m)**2)-m*cos(theta))/
     &		(sqrt(1d0-(sin(theta)/m)**2)+m*cos(theta)))**2

	return
	end
	



	subroutine ComputeSurface()
	use GlobalSetup
	IMPLICIT NONE
	real*8 SurfEmis(nlam),e1(nlam),e2(nlam),tot,f_ice,f_grass,f_snow,f_sand,f_water
	external Enstatite_X,Enstatite_Y,Enstatite_Z
	external Forsterite_X,Forsterite_Y,Forsterite_Z
	external Labradorite_X,Labradorite_Y,Labradorite_Z,SiO2,FeO
	integer i

c=========================================
c compute emissivity of the surface
c=========================================

	select case(surfacetype)
		case("MIXED","MIX","mixed","mix")
			surface_emis=0d0
			tot=0.6d0
			call RegridDataLNK(Enstatite_X,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)/3d0
			call RegridDataLNK(Enstatite_Y,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)/3d0
			call RegridDataLNK(Enstatite_Z,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)/3d0

			tot=0.3d0
			call RegridDataLNK(Forsterite_X,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)/3d0
			call RegridDataLNK(Forsterite_Y,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)/3d0
			call RegridDataLNK(Forsterite_Z,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)/3d0

			tot=0.1d0
			call RegridDataLNK(SiO2,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)
		case("QUARTZ","quartz","SiO2")
			surface_emis=0d0
			tot=1d0
			call RegridDataLNK(SiO2,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)
		case("FeO")
			surface_emis=0d0
			tot=1d0
			call RegridDataLNK(FeO,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)
		case("labradorite","LABRADORITE")
			surface_emis=0d0
			tot=1.0d0
			call RegridDataLNK(Labradorite_X,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)/3d0
			call RegridDataLNK(Labradorite_Y,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)/3d0
			call RegridDataLNK(Labradorite_Z,lam(1:nlam)*1d4,e1(1:nlam),e2(1:nlam),nlam,.true.)
			call computeHapke(SurfEmis(1:nlam),e1(1:nlam),e2(1:nlam),nlam,computelam(1:nlam))
			surface_emis(1:nlam)=surface_emis(1:nlam)+tot*SurfEmis(1:nlam)/3d0
		case("BLACK","black")
			surface_emis(1:nlam)=1.0
		case("GREY","grey")
			surface_emis(1:nlam)=1d0-surfacealbedo
		case("WHITE","white")
			surface_emis(1:nlam)=1d-4
		case("Earth","EARTH","earth")
c ice fraction according to Ramirez 2023
			if(Tsurface.lt.239d0) then
				f_ice=1d0
			else if(Tsurface.lt.273.15) then
				f_ice=1d0-exp((Tsurface-273.15)/12.5)
			else
				f_ice=0d0
			endif
			if(Tsurface.gt.308.15) then
				f_sand=1d0
			else if(Tsurface.gt.273.15) then
				f_sand=((Tsurface-273.15)/(308.15-273.15))**4
			else
				f_sand=0d0
			endif
			f_grass=(1d0-f_ice)*(1d0-f_surface_water)*(1d0-f_sand)
			f_sand=(1d0-f_ice)*(1d0-f_surface_water)*f_sand
			f_water=(1d0-f_ice)*f_surface_water
			f_snow=f_ice*(1d0-f_surface_water/2d0)
			f_ice=f_ice*f_surface_water/2d0
			surface_emis=f_grass*surface_emis_grass+f_sand*surface_emis_sand+
     &					 f_water*surface_emis_water+f_snow*surface_emis_snow+f_ice*surface_emis_ice
			print*,'grass: ',f_grass
			print*,'sand:  ',f_sand
			print*,'water: ',f_water
			print*,'snow:  ',f_snow
			print*,'ice:   ',f_ice
		case("FILE","file")
		case default
			call output("Surface type not known!")
			stop
	end select

	if(.not.retrieval.and..not.dopostequalweights) then
		open(unit=93,file=trim(outputdir) // 'surfemis.dat',FORM="FORMATTED",ACCESS="STREAM")
		do i=1,nlam
			write(93,*) lam(i),surface_emis(i)
		enddo
		close(unit=93)
	endif

c=========================================
c=========================================


	return
	end




	real*8 function NormSig(x,a,c,x0,x1)
	IMPLICIT NONE
	real*8 x,a,y1,y2,xx,c,x0,x1,yy

	if(x0.eq.x1) then
		NormSig=x0
		return
	endif

	xx=(0d0-c)*a
	if(xx.lt.-40d0) then
		y1=exp(xx)
	else if(xx.gt.40d0) then
		y1=1d0
	else if(abs(xx).lt.0.1) then
		y1=0.5d0+xx/4d0-xx**3/48d0
	else
		y1=1d0/(exp(-xx)+1d0)
	endif
	xx=(1d0-c)*a
	if(xx.lt.-40d0) then
		y2=exp(xx)
	else if(xx.gt.40d0) then
		y2=1d0
	else if(abs(xx).lt.0.1) then
		y2=0.5d0+xx/4d0-xx**3/48d0
	else
		y2=1d0/(exp(-xx)+1d0)
	endif

	if(x1.eq.x0) then
		yy=1d0
	else
		xx=(x-x0)/(x1-x0)
		xx=(xx-c)*a
		if(xx.lt.-40d0) then
			yy=exp(xx)
		else if(xx.gt.40d0) then
			yy=1d0
		else if(abs(xx).lt.0.1) then
			yy=0.5d0+xx/4d0-xx**3/48d0
		else
			yy=1d0/(exp(-xx)+1d0)
		endif
	endif
	if(y2.eq.y1) then
		NormSig=1d0
	else
		NormSig=(yy-y1)/(y2-y1)
	endif
	
	return
	end

