	subroutine Raytrace()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 rr,xx1,xx2,si,exp_tau,A,d,s,fluxg,Planck,fact,tau,freq0,tau_a,tautot,Ag
	real*8 Ca,Cs,tot,contr,dP,tot2,TT
	integer icloud,isize
	real*8,allocatable :: rtrace(:),phase0(:)
	real*8,allocatable :: fluxg_contr(:),fact_contr(:),Ag_contr(:)
	integer irc,imolhide,jmolhide
	integer nrtrace,ndisk,i,ir,ir_next,ilam,ig,nsub,j,k,nk
	logical in,dohide
	integer icc,imol
	character*500 filename
	real*8,allocatable :: dtrace(:,:),CaCont(:,:),Ca_cloud(:,:),Cs_cloud(:,:),BBr(:)
	integer,allocatable :: irtrace(:,:),nirtrace(:)

	docloud=.false.
	do icc=2,ncc
		docloud(icc,1:nclouds)=docloud(icc-1,1:nclouds)
		i=0
10		i=i+1
		docloud(icc,i)=.not.docloud(icc,i)
		if(.not.docloud(icc,i)) goto 10
	enddo
	do icc=1,ncc
		cloudfrac(icc)=1d0
		do i=1,nclouds
			if(docloud(icc,i)) then
				cloudfrac(icc)=cloudfrac(icc)*Cloud(i)%coverage
			else
				cloudfrac(icc)=cloudfrac(icc)*(1d0-Cloud(i)%coverage)
			endif
		enddo
	enddo

	if(ncc.eq.1) then
		docloud=.true.
		cloudfrac=1d0
	endif

	if(.not.retrieval) then
	do icc=1,ncc
		do ilam=1,nlam
		tau=0d0
		do ir=nr,1,-1
			Ca=0d0
			Cs=0d0
			do icloud=1,nclouds
				if(docloud(icc,icloud)) then
					Ca=Ca+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
					Cs=Cs+Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
				endif
			enddo
			tau=tau+(R(ir+1)-R(ir))*(Ca+Cs)
		enddo
		cloudtau(icc,ilam)=tau
		enddo
	enddo
	endif

	flux=0d0

	call output("==================================================================")

	if(scattering) then
	call InitRandomWalk()
	call output("Scattered light contributions")

	call tellertje(1,nlam)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,fluxg,icc,phase0)
!$OMP& SHARED(nlam,nclouds,flux,phase,ncc,docloud,cloudfrac,nphase)
	allocate(phase0(nphase))
!$OMP DO SCHEDULE(DYNAMIC)
	do ilam=1,nlam
		call tellertje(ilam+1,nlam+1)
		flux(:,ilam)=0d0
		phase(:,:,ilam)=0d0
		do icc=1,ncc
			if(cloudfrac(icc).gt.0d0) then
				call MCRad(ilam,fluxg,phase0,docloud(icc,1:nclouds))
			endif
			flux(0,ilam)=flux(0,ilam)+cloudfrac(icc)*fluxg
			flux(icc,ilam)=flux(icc,ilam)+fluxg
			phase(1:nphase,icc,ilam)=phase(1:nphase,icc,ilam)+
     &				phase0(1:nphase)
			phase(1:nphase,0,ilam)=phase(1:nphase,0,ilam)+
     &				cloudfrac(icc)*phase0(1:nphase)
		enddo
	enddo
!$OMP END DO
	deallocate(phase0)
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(nlam,nlam)
	endif

	call output("Raytracing over the planet disk")

	ndisk=15
	nsub=3

	if(.not.transspec) then
		ndisk=15
		nsub=0
	endif
	if(.not.emisspec) then
		ndisk=2
		nsub=3
	endif
	
	nrtrace=(nr-1)*nsub+ndisk
	allocate(rtrace(nrtrace))

	k=0
	if(nsub.eq.0) then
		do i=1,ndisk
			k=k+1
			rtrace(k)=R(nr)*real(i-1)/real(ndisk-1)
		enddo
	else
		do i=1,ndisk
			k=k+1
			rtrace(k)=R(1)*real(i-1)/real(ndisk-1)
		enddo
	endif
	do i=1,nr-1
		do j=1,nsub
			k=k+1
			rtrace(k)=R(i)+(R(i+1)-R(i))*real(j-1)/real(nsub)
		enddo
	enddo

	allocate(dtrace(nr*2,nrtrace))
	allocate(irtrace(nr*2+1,nrtrace))
	allocate(nirtrace(nrtrace))
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,rr,ir,k,si,xx1,in,xx2,d,ir_next)
!$OMP& SHARED(nrtrace,rtrace,R,dtrace,irtrace,nirtrace,nr)
!$OMP DO SCHEDULE(DYNAMIC)
	do i=1,nrtrace-1
		rr=sqrt(rtrace(i)*rtrace(i+1))
		ir=nr
		k=0
		si=-1d0
		xx1=si*sqrt(abs(R(ir+1)**2-rr**2))
		in=.true.
1		continue
		k=k+1
		if(in) then
			xx2=(R(ir)**2-rr**2)
			if(xx2.gt.0d0) then
				xx2=si*sqrt(abs(xx2))
				d=abs(xx1-xx2)
				ir_next=ir-1
				goto 2
			else
				si=-si
				xx2=si*sqrt(abs(R(ir+1)**2-rr**2))
				d=abs(xx1-xx2)
				ir_next=ir+1
				in=.false.
				goto 2
			endif
		else
			xx2=si*sqrt(abs(R(ir+1)**2-rr**2))
			d=abs(xx1-xx2)
			ir_next=ir+1
			goto 2
		endif
2		continue
		irtrace(k,i)=ir
		dtrace(k,i)=d
		if(ir_next.gt.0.and.ir_next.le.nr) then
			ir=ir_next
			xx1=xx2
			goto 1
		endif
		if(ir_next.le.0) then
			k=k+1
			irtrace(k,i)=ir_next
			dtrace(k,i)=0d0
		endif
		nirtrace(i)=k
		irtrace(k+1,i)=nr+1
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	if(emisspec.or.computecontrib) then

	if(.not.allocated(flux_contr).and.computecontrib) then
		allocate(flux_contr(nr,nlam))
		allocate(obsA_contr(nr,nlam))
	endif
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(ilam,freq0,ig,i,fluxg,fact,A,rr,ir,si,xx1,in,xx2,d,ir_next,tau,exp_tau,tau_a,tautot,Ag,Ca_cloud,Cs_cloud,
!$OMP&         Ca,Cs,icloud,isize,BBr,imol,irc,contr,fact_contr,fluxg_contr,Ag_contr,nk)
!$OMP& SHARED(nlam,freq,obsA,flux,cloudfrac,ncc,docloud,nrtrace,ng,rtrace,nr,R,Ndens,Cabs,Csca,T,lam,maxtau,nclouds,Cloud,
!$OMP&			cloud_dens,P,flux_contr,obsA_contr,irtrace,dtrace,nirtrace,
!$OMP&			nmol,mixrat_r,includemol,computecontrib,wgg)
	allocate(fact_contr(nr))
	allocate(fluxg_contr(nr))
	allocate(Ag_contr(nr))
	allocate(Ca_cloud(ncc,nr),Cs_cloud(ncc,nr),BBr(0:nr))
!$OMP DO SCHEDULE(DYNAMIC)
	do ilam=1,nlam
		Ca_cloud(1:ncc,1:nr)=0d0
		Cs_cloud(1:ncc,1:nr)=0d0
		if(nclouds.gt.0) then
		do icc=1,ncc
			if(cloudfrac(icc).gt.0d0) then
				do ir=1,nr
					do icloud=1,nclouds
						if(docloud(icc,icloud)) then
							Ca_cloud(icc,ir)=Ca_cloud(icc,ir)+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
							Cs_cloud(icc,ir)=Cs_cloud(icc,ir)+Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
						endif
					enddo
				enddo
			endif
		enddo
		endif
		freq0=freq(ilam)
		obsA(:,ilam)=0d0
		if(computecontrib) then
			obsA_contr(1:nr,ilam)=0d0
			flux_contr(1:nr,ilam)=flux(0,ilam)
		endif
		BBr(0)=Planck(Tsurface,freq0)
		do ir=1,nr
			BBr(ir)=Planck(T(ir),freq0)
		enddo
		do ig=1,ng
			do icc=1,ncc
			if(cloudfrac(icc).gt.0d0) then
				fluxg=0d0
				Ag=0d0
				fluxg_contr=0d0
				Ag_contr=0d0
				ir_next=0
				do i=1,nrtrace-1
					fact=1d0
					fact_contr=1d0
					tautot=0d0
					A=pi*(rtrace(i+1)**2-rtrace(i)**2)
					nk=nirtrace(i)
					do k=1,nk

					ir=irtrace(k,i)
					d=dtrace(k,i)

					Ca=Cabs(ir,ilam,ig,0)*Ndens(ir)+Ca_cloud(icc,ir)
					Cs=Csca(ir,ilam)*Ndens(ir)+Cs_cloud(icc,ir)
					tau_a=d*Ca
					tau=tau_a+d*Cs
					exp_tau=exp(-tau)
					tautot=tautot+tau
					contr=A*BBr(ir)*(1d0-exp_tau)*fact*tau_a/tau
					fluxg=fluxg+contr
					fact=fact*exp_tau
					if(computecontrib) then
						do irc=1,nr
							if(ir.ne.irc) then
								contr=A*BBr(ir)*(1d0-exp_tau)*fact_contr(irc)*tau_a/tau
								fluxg_contr(irc)=fluxg_contr(irc)+contr
								fact_contr(irc)=fact_contr(irc)*exp_tau
							endif
						enddo
					endif
					if(k.lt.nk) ir_next=irtrace(k+1,i)
					if(ir_next.le.0.or.tautot.ge.maxtau) exit

					enddo

					if(ir_next.le.0) then
						fluxg=fluxg+A*BBr(0)*fact
						if(computecontrib) then
							do irc=1,nr
								if(ir.ne.irc) then
									contr=A*BBr(0)*fact_contr(irc)
									fluxg_contr(irc)=fluxg_contr(irc)+contr
								endif
							enddo
						endif
						fact_contr=0d0
					endif
					if(ir_next.le.0.or.tautot.ge.maxtau) fact=0d0
					Ag=Ag+A*(1d0-fact)
					if(computecontrib) then
						Ag_contr=Ag_contr+A*(1d0-fact_contr)
					endif
				enddo
				flux(0,ilam)=flux(0,ilam)+cloudfrac(icc)*fluxg*wgg(ig)
				obsA(0,ilam)=obsA(0,ilam)+cloudfrac(icc)*Ag*wgg(ig)
				if(computecontrib) then
					do irc=1,nr
						flux_contr(irc,ilam)=flux_contr(irc,ilam)+cloudfrac(icc)*fluxg_contr(irc)*wgg(ig)
						obsA_contr(irc,ilam)=obsA_contr(irc,ilam)+cloudfrac(icc)*Ag_contr(irc)*wgg(ig)
					enddo
				endif
				flux(icc,ilam)=flux(icc,ilam)+fluxg*wgg(ig)
				obsA(icc,ilam)=obsA(icc,ilam)+Ag*wgg(ig)
			endif
			enddo
		enddo
	enddo
!$OMP END DO
	deallocate(fact_contr)
	deallocate(fluxg_contr)
	deallocate(Ag_contr)
	deallocate(Ca_cloud,Cs_cloud,BBr)
!$OMP FLUSH
!$OMP END PARALLEL
	endif

	if(transspec) then
	if(.not.allocated(obsA_LC)) allocate(obsA_LC(nrtrace,nlam))
	allocate(CaCont(nr,nlam))


	if(dotranshide) then
	do jmolhide=1,2
	do imolhide=0,nmol
	obsA=0d0
	obsA_LC=0d0
	dohide=.false.
	if(imolhide.eq.0) then
		if(nclouds.gt.0) dohide=.true.
	else if(opacitymol(imolhide)) then
		dohide=.true.
	endif
	if(dohide) then
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,i,icc,imol,ig,tautot,Ag,A,ir,d,tau,k,ir_next,Ca,icloud,nk)
!$OMP& SHARED(nlam,ncc,nrtrace,nmol,ng,opacitymol,irtrace,dtrace,Cabs_mol,mixrat_r,
!$OMP&		wgg,P,Cloud,nclouds,cloud_dens,obsA,rtrace,docloud,Cext_cont,Rstar,
!$OMP&		cloudfrac,nirtrace,maxtau,ndisk,nr,obsA_LC,CaCont,imolhide,jmolhide)
!$OMP DO SCHEDULE(DYNAMIC)
	do ilam=1,nlam
		do icc=1,ncc
		if(cloudfrac(icc).gt.0d0) then
			if((imolhide.eq.0.and.jmolhide.eq.1).or.(imolhide.ne.0.and.jmolhide.eq.2)) then
				do ir=1,nr
					CaCont(ir,ilam)=Cext_cont(ir,ilam)
				enddo
			else
				do ir=1,nr
					Ca=Cext_cont(ir,ilam)
					do icloud=1,nclouds
						if(docloud(icc,icloud)) then
							Ca=Ca+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
							Ca=Ca+Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
						endif
					enddo
					CaCont(ir,ilam)=Ca
				enddo
			endif
			do i=1,nrtrace-1
				if(i.lt.ndisk) then !.or.jmolhide.eq.2.and.rtrace(i).lt.sqrt(0.0195)*Rstar) then
					A=0d0
				else
					A=1d0
					nk=nirtrace(i)
					tautot=0d0
					do k=1,nk
						ir=irtrace(k,i)
						d=dtrace(k,i)
						tau=d*CaCont(ir,ilam)
						tautot=tautot+tau
						ir_next=irtrace(k+1,i)
						if(tautot.gt.maxtau) then
							tautot=1d4
							A=0d0
							goto 19
						endif
						if(ir_next.lt.1) then
							tautot=1d4
							A=0d0
							goto 19
						endif
					enddo
					A=A*exp(-tautot)
					do imol=1,nmol
						if((opacitymol(imol).and.imol.ne.imolhide.and.jmolhide.eq.1).or.
     &					   (opacitymol(imol).and.imol.eq.imolhide.and.jmolhide.eq.2)) then
						Ag=0d0
						do ig=1,ng
							tautot=0d0
							do k=1,nk
								ir=irtrace(k,i)
								d=dtrace(k,i)
								tau=d*Cabs_mol(ig,ilam,imol,ir)
								tautot=tautot+tau
								ir_next=irtrace(k+1,i)
								if(tautot.gt.maxtau) goto 18
								if(ir_next.lt.1) goto 18
							enddo
							Ag=Ag+exp(-tautot)*wgg(ig)
18							continue
						enddo
						A=A*min(Ag,1d0)
						endif
					enddo
19					continue
				endif
				if(.not.A.gt.0d0) A=0d0
				obsA_LC(i,ilam)=obsA_LC(i,ilam)+cloudfrac(icc)*A
				A=pi*(rtrace(i+1)**2-rtrace(i)**2)*(1d0-A)
				obsA(icc,ilam)=obsA(icc,ilam)+A
				obsA(0,ilam)=obsA(0,ilam)+cloudfrac(icc)*A
			enddo
		endif
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	if(jmolhide.eq.1) then
		if(imolhide.eq.0) then
			open(unit=83,file=trim(outputdir) // "/trans_hide_clouds",FORM="FORMATTED",ACCESS="STREAM")
		else
			open(unit=83,file=trim(outputdir) // "/trans_hide_" // trim(molname(imolhide)),FORM="FORMATTED",ACCESS="STREAM")
		endif
		do ilam=1,nlam
			if(computelam(ilam)) write(83,*) lam(ilam)*1d4,obsA(0,ilam)/(pi*Rstar**2)
		enddo
	else
		if(imolhide.eq.0) then
			open(unit=83,file=trim(outputdir) // "/trans_only_clouds",FORM="FORMATTED",ACCESS="STREAM")
		else
			open(unit=83,file=trim(outputdir) // "/trans_only_" // trim(molname(imolhide)),FORM="FORMATTED",ACCESS="STREAM")
		endif
		do ilam=1,nlam
			if(computelam(ilam)) write(83,*) lam(ilam)*1d4,obsA(0,ilam)/(pi*Rstar**2)
		enddo
	endif
	close(unit=83)
	endif
	enddo
	enddo
	endif


	if(computecontrib) then
	obsA_contr=0d0

	do j=1,nr
		call tellertje(j,nr)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,i,icc,imol,ig,tautot,Ag,A,ir,d,tau,k,ir_next,Ca,icloud,nk)
!$OMP& SHARED(nlam,ncc,nrtrace,nmol,ng,opacitymol,irtrace,dtrace,Cabs_mol,mixrat_r,
!$OMP&		wgg,P,Cloud,nclouds,cloud_dens,obsA,rtrace,docloud,Cext_cont,
!$OMP&		cloudfrac,nirtrace,maxtau,ndisk,nr,obsA_contr,CaCont,j)
!$OMP DO SCHEDULE(DYNAMIC)
	do ilam=1,nlam
		do icc=1,ncc
		if(cloudfrac(icc).gt.0d0) then
			do ir=1,nr
				Ca=Cext_cont(ir,ilam)
				do icloud=1,nclouds
					if(docloud(icc,icloud)) then
						Ca=Ca+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
						Ca=Ca+Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
					endif
				enddo
				CaCont(ir,ilam)=Ca
			enddo
			do i=1,nrtrace-1
				if(i.lt.ndisk) then
					A=0d0
				else
					A=1d0
					nk=nirtrace(i)
					tautot=0d0
					do k=1,nk
						ir=irtrace(k,i)
						d=dtrace(k,i)
						if(ir.ne.j) then
							tau=d*CaCont(ir,ilam)
						else
							tau=0d0
						endif
						tautot=tautot+tau
						ir_next=irtrace(k+1,i)
						if(tautot.gt.maxtau) then
							tautot=1d4
							A=0d0
							goto 29
						endif
						if(ir_next.lt.1) then
							tautot=1d4
							A=0d0
							goto 29
						endif
					enddo
					A=A*exp(-tautot)
					do imol=1,nmol
						if(opacitymol(imol)) then
						Ag=0d0
						do ig=1,ng
							tautot=0d0
							do k=1,nk
								ir=irtrace(k,i)
								d=dtrace(k,i)
								if(ir.ne.j) then
									tau=d*Cabs_mol(ig,ilam,imol,ir)
								else
									tau=0d0
								endif
								tautot=tautot+tau
								ir_next=irtrace(k+1,i)
								if(tautot.gt.maxtau) goto 28
								if(ir_next.lt.1) goto 28
							enddo
							Ag=Ag+exp(-tautot)*wgg(ig)
28							continue
						enddo
						A=A*min(Ag,1d0)
						endif
					enddo
29					continue
				endif
				if(.not.A.gt.0d0) A=0d0
				A=pi*(rtrace(i+1)**2-rtrace(i)**2)*(1d0-A)
				obsA_contr(j,ilam)=obsA_contr(j,ilam)+cloudfrac(icc)*A
			enddo
		endif
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	enddo
	endif



	obsA=0d0
	obsA_LC=0d0

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,i,icc,imol,ig,tautot,Ag,A,ir,d,tau,k,ir_next,Ca,icloud,nk)
!$OMP& SHARED(nlam,ncc,nrtrace,nmol,ng,opacitymol,irtrace,dtrace,Cabs_mol,mixrat_r,
!$OMP&		wgg,P,Cloud,nclouds,cloud_dens,obsA,rtrace,docloud,Cext_cont,
!$OMP&		cloudfrac,nirtrace,maxtau,ndisk,nr,obsA_LC,CaCont)
!$OMP DO SCHEDULE(DYNAMIC)
	do ilam=1,nlam
		do icc=1,ncc
		if(cloudfrac(icc).gt.0d0) then
			do ir=1,nr
				Ca=Cext_cont(ir,ilam)
				do icloud=1,nclouds
					if(docloud(icc,icloud)) then
						Ca=Ca+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
						Ca=Ca+Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
					endif
				enddo
				CaCont(ir,ilam)=Ca
			enddo
			do i=1,nrtrace-1
				if(i.lt.ndisk) then
					A=0d0
				else
					A=1d0
					nk=nirtrace(i)
					tautot=0d0
					do k=1,nk
						ir=irtrace(k,i)
						d=dtrace(k,i)
						tau=d*CaCont(ir,ilam)
						tautot=tautot+tau
						ir_next=irtrace(k+1,i)
						if(tautot.gt.maxtau) then
							tautot=1d4
							A=0d0
							goto 9
						endif
						if(ir_next.lt.1) then
							tautot=1d4
							A=0d0
							goto 9
						endif
					enddo
					A=A*exp(-tautot)
					do imol=1,nmol
						if(opacitymol(imol)) then
						Ag=0d0
						do ig=1,ng
							tautot=0d0
							do k=1,nk
								ir=irtrace(k,i)
								d=dtrace(k,i)
								tau=d*Cabs_mol(ig,ilam,imol,ir)
								tautot=tautot+tau
								ir_next=irtrace(k+1,i)
								if(tautot.gt.maxtau) goto 8
								if(ir_next.lt.1) goto 8
							enddo
							Ag=Ag+exp(-tautot)*wgg(ig)
8							continue
						enddo
						A=A*min(Ag,1d0)
						endif
					enddo
9					continue
				endif
				if(.not.A.gt.0d0) A=0d0
				obsA_LC(i,ilam)=obsA_LC(i,ilam)+cloudfrac(icc)*A
				A=pi*(rtrace(i+1)**2-rtrace(i)**2)*(1d0-A)
				obsA(icc,ilam)=obsA(icc,ilam)+A
				obsA(0,ilam)=obsA(0,ilam)+cloudfrac(icc)*A
			enddo
		endif
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	deallocate(CaCont)
	obsA_LC=1d0-obsA_LC
	if(computeLC) then
c		call LightCurve(rtrace,nrtrace)
c		if(retrieval) call LightCurveRetrieval_Fit(rtrace,nrtrace)
		call LightCurveRetrieval_Fit(rtrace,nrtrace)
c		call LightCurveRetrieval(rtrace,nrtrace)
	endif
	endif

	deallocate(dtrace)
	deallocate(irtrace)
	deallocate(nirtrace)
	
	if(computecontrib) then
		filename=trim(outputdir) // "contribution.fits.gz"
		do irc=1,nr
			if(irc.eq.1) then
				dP=log10(P(1)/P(2))
			else if(irc.eq.nr) then
				dP=log10(P(nr-1)/P(nr))
			else
				dP=log10(P(irc-1)/P(irc+1))
			endif
			flux_contr(irc,1:nlam)=abs(flux(0,1:nlam)-flux_contr(irc,1:nlam))/dP
			obsA_contr(irc,1:nlam)=abs(obsA(0,1:nlam)-obsA_contr(irc,1:nlam))/dP
		enddo
		do ilam=1,nlam
			tot=sum(flux_contr(1:nr,ilam))
			flux_contr(1:nr,ilam)=flux_contr(1:nr,ilam)/tot
			tot=sum(obsA_contr(1:nr,ilam))
			obsA_contr(1:nr,ilam)=obsA_contr(1:nr,ilam)/tot
		enddo
		call writeContribution(filename,P,lam,obsA_contr,flux_contr,nr,nlam)
	endif

	flux=flux*1d23/distance**2
	phase=phase*1d23/distance**2

	deallocate(rtrace)

	return
	end
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	real*8 function Planck(T,nu)
	use Constants
	IMPLICIT NONE
	real*8 T,nu,x

	x=hplanck*nu*clight/(kb*T)
	if (x.gt.40d0) then
		Planck=(2d0*hplanck*nu**3*clight)*exp(-x)
	else if (x.lt.0.1) then
		Planck=(2d0*hplanck*nu**3*clight)*(-0.5d0+1.d0/x+x/12.d0-x**3/720.d0)
	else
		Planck=(2d0*hplanck*nu**3*clight)/(exp(x)-1d0)
	endif

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	real*8 function dPlanck(T,nu)
	use Constants
	IMPLICIT NONE
	real*8 T,nu,x,expx

	x=hplanck*nu*clight/(kb*T)
	expx=exp(x)
	dPlanck=(2d0*nu**2*x**2*kb*expx)/((expx-1d0)**2)

	return
	end

