	subroutine Raytrace()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 rr,xx1,xx2,si,exp_tau,A,d,s,fluxg,Planck,fact,tau,freq0,tau_a,tautot,Ag
	real*8 Ca,Cs,BBr(nr),tot,contr
	integer icloud,isize
	real*8,allocatable :: rtrace(:),phase0(:),ptrace(:)
	real*8,allocatable :: fluxg_contr(:),fact_contr(:),Ag_contr(:)
	integer irc
	integer nrtrace,ndisk,i,ir,ir_next,ilam,ig,nsub,j,k
	logical in
	integer icc,imol
	real*8 Ocolumn(2,nlam,ncc),Ccolumn(2,nlam,ncc),Hcolumn(2,nlam,ncc),Otot,Ctot,Htot,dP
	character*500 filename
	real*8,allocatable :: dtrace(:,:)
	integer,allocatable :: irtrace(:,:),nirtrace(:)

	docloud=.false.
	if(singlecloud) then
	cloudfrac(1)=0d0
	do i=1,nclouds
		cloudfrac(i+1)=Cloud(i)%coverage
		docloud(i+1,i)=.true.
	enddo
	else
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
	endif

	do icc=1,ncc
		do ilam=1,nlam
		tau=0d0
		Otot=0d0
		Ctot=0d0
		Htot=0d0
		do ir=nr,2,-1
			Ca=sum(wgg(1:ng)*Cabs(ir,ilam,1:ng))*Ndens(ir)
			Cs=Csca(ir,ilam)*Ndens(ir)
			do icloud=1,nclouds
				if(docloud(icc,icloud)) then
					if(useDRIFT) then
						Ca=Ca+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
						Cs=Cs+Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
					else
						do isize=1,Cloud(icloud)%nr
							Ca=Ca+
     &		Cloud(icloud)%Kabs(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
							Cs=Cs+
     &		Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
						enddo
					endif
				endif
			enddo
			tau_a=(R(ir)-R(ir-1))*(Ca+Cs)
			if((tau+tau_a).gt.1d0) then
				d=(1d0-tau)/tau_a
				tau1depth(icc,ilam)=10d0**(log10(P(ir))+log10(P(ir-1)/P(ir))*d)
				d=d*(R(ir)-R(ir-1))
				do imol=1,nmol
					if(includemol(imol)) then
						Otot=Otot+d*Ndens(ir)*mixrat_r(ir,imol)*real(Oatoms(imol))
						Ctot=Ctot+d*Ndens(ir)*mixrat_r(ir,imol)*real(Catoms(imol))
						Htot=Htot+d*Ndens(ir)*mixrat_r(ir,imol)*real(Hatoms(imol))
					endif
				enddo
				goto 3
			endif
			tau=tau+tau_a
			d=(R(ir)-R(ir-1))
			do imol=1,nmol
				if(includemol(imol)) then
					Otot=Otot+d*Ndens(ir)*mixrat_r(ir,imol)*real(Oatoms(imol))
					Ctot=Ctot+d*Ndens(ir)*mixrat_r(ir,imol)*real(Catoms(imol))
					Htot=Htot+d*Ndens(ir)*mixrat_r(ir,imol)*real(Hatoms(imol))
				endif
			enddo
		enddo
		tau1depth(icc,ilam)=P(1)
3		continue
		Ocolumn(2,ilam,icc)=Otot
		Ccolumn(2,ilam,icc)=Ctot
		Hcolumn(2,ilam,icc)=Htot
		enddo

		do ilam=1,nlam
		tau=0d0
		do ir=nr,2,-1
			Ca=0d0
			Cs=0d0
			do icloud=1,nclouds
				if(docloud(icc,icloud)) then
					if(useDRIFT) then
						Ca=Ca+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
						Cs=Cs+Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
					else
						do isize=1,Cloud(icloud)%nr
							Ca=Ca+
     &		Cloud(icloud)%Kabs(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
							Cs=Cs+
     &		Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
						enddo
					endif
				endif
			enddo
			tau=tau+(R(ir)-R(ir-1))*(Ca+Cs)
		enddo
		cloudtau(icc,ilam)=tau
		enddo
	enddo
		

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
!$OMP DO SCHEDULE(STATIC,1)
	do ilam=1,nlam-1
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

	ndisk=10
	nsub=3
	
	nrtrace=(nr-1)*nsub+ndisk
	allocate(rtrace(nrtrace))
	allocate(ptrace(nrtrace))

	k=0
	do i=1,ndisk
		k=k+1
		rtrace(k)=Rplanet*real(i-1)/real(ndisk)
		ptrace(k)=P(1)*(1d0+real(ndisk-i)/real(ndisk))
	enddo
	do i=1,nr-1
		do j=1,nsub
			k=k+1
			rtrace(k)=R(i)+(R(i+1)-R(i))*real(j-1)/real(nsub)
			ptrace(k)=P(i)+(P(i+1)-P(i))*real(j-1)/real(nsub)
		enddo
	enddo

	allocate(dtrace(nrtrace,nr*2))
	allocate(irtrace(nrtrace,nr*2))
	allocate(nirtrace(nrtrace))
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,rr,ir,k,si,xx1,in,xx2,d,ir_next)
!$OMP& SHARED(nrtrace,rtrace,R,dtrace,irtrace,nirtrace,nr)
!$OMP DO SCHEDULE(STATIC,1)
	do i=1,nrtrace-1
		rr=sqrt(rtrace(i)*rtrace(i+1))
		ir=nr
		k=0
		si=-1d0
		xx1=si*sqrt(R(ir+1)**2-rr**2)
		in=.true.
1		continue
		k=k+1
		if(in) then
			xx2=(R(ir)**2-rr**2)
			if(xx2.gt.0d0) then
				xx2=si*sqrt(xx2)
				d=abs(xx1-xx2)
				ir_next=ir-1
				goto 2
			else
				si=-si
				xx2=si*sqrt(R(ir+1)**2-rr**2)
				d=abs(xx1-xx2)
				ir_next=ir+1
				in=.false.
				goto 2
			endif
		else
			xx2=si*sqrt(R(ir+1)**2-rr**2)
			d=abs(xx1-xx2)
			ir_next=ir+1
			goto 2
		endif
2		continue
		irtrace(i,k)=ir
		dtrace(i,k)=d
		if(ir_next.gt.0.and.ir_next.le.nr) then
			ir=ir_next
			xx1=xx2
			goto 1
		endif
		if(ir_next.le.0) then
			k=k+1
			irtrace(i,k)=ir_next
			dtrace(i,k)=0d0
		endif
		nirtrace(i)=k
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL


	if(emisspec.or.computecontrib) then

	if(.not.allocated(flux_contr)) then
		allocate(flux_contr(nr,nlam))
		allocate(obsA_contr(nr,nlam))
	endif
	Ocolumn=0d0
	Ccolumn=0d0
	Hcolumn=0d0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(ilam,freq0,ig,i,fluxg,fact,A,rr,ir,si,xx1,in,xx2,d,ir_next,tau,exp_tau,tau_a,tautot,Ag,
!$OMP&         Ca,Cs,icloud,isize,BBr,Otot,Ctot,Htot,imol,irc,contr,fact_contr,fluxg_contr,Ag_contr)
!$OMP& SHARED(nlam,freq,obsA,flux,cloudfrac,ncc,docloud,nrtrace,ng,rtrace,nr,R,Ndens,Cabs,Csca,T,lam,maxtau,nclouds,Cloud,
!$OMP&			cloud_dens,useDRIFT,Psimplecloud,P,flux_contr,obsA_contr,irtrace,dtrace,nirtrace,
!$OMP&			Ocolumn,Ccolumn,Hcolumn,nmol,mixrat_r,includemol,computecontrib,wgg)
	allocate(fact_contr(nr))
	allocate(fluxg_contr(nr))
	allocate(Ag_contr(nr))
!$OMP DO SCHEDULE(STATIC,1)
	do ilam=1,nlam-1
		freq0=sqrt(freq(ilam)*freq(ilam+1))
		obsA(:,ilam)=0d0
		obsA_contr(1:nr,ilam)=0d0
		flux_contr(1:nr,ilam)=flux(0,ilam)
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
				do i=1,nrtrace-1
					fact=1d0
					fact_contr=1d0
					tautot=0d0
					A=pi*(rtrace(i+1)**2-rtrace(i)**2)
					do k=1,nirtrace(i)

					ir=irtrace(i,k)
					d=dtrace(i,k)

					Ca=Cabs(ir,ilam,ig)*Ndens(ir)
					Cs=Csca(ir,ilam)*Ndens(ir)
					do icloud=1,nclouds
						if(docloud(icc,icloud)) then
							if(useDRIFT) then
								Ca=Ca+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
								Cs=Cs+Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
							else
								do isize=1,Cloud(icloud)%nr
									Ca=Ca+
     &		Cloud(icloud)%Kabs(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
									Cs=Cs+
     &		Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
								enddo
							endif
						endif
					enddo
					tau_a=d*Ca
					tau=tau_a+d*Cs
					if(P(ir).gt.Psimplecloud) tau=1d4
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
					if(k.lt.nirtrace(i)) ir_next=irtrace(i,k+1)
					if(ir_next.le.0.or.tautot.ge.maxtau) exit

					enddo

					if(ir_next.le.0) then
						fluxg=fluxg+A*BBr(ir)*fact
						if(computecontrib) then
							do irc=1,nr
								if(ir.ne.irc) then
									contr=A*BBr(ir)*fact_contr(irc)
									fluxg_contr(irc)=fluxg_contr(irc)+contr
								endif
							enddo
						endif
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
!$OMP FLUSH
!$OMP END PARALLEL

	endif

	if(transspec) then
	obsA=0d0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,i,icc,imol,ig,tautot,Ag,A,ir,d,tau,k,ir_next,Ca,icloud)
!$OMP& SHARED(nlam,ncc,nrtrace,nmol,ng,includemol,irtrace,dtrace,Cabs_mol,mixrat_r,
!$OMP&		wgg,P,Cloud,nclouds,cloud_dens,obsA,rtrace,docloud,useDRIFT,Cext_cont,
!$OMP&		cloudfrac,nirtrace,Psimplecloud,maxtau)
!$OMP DO SCHEDULE(STATIC,1)
	do ilam=1,nlam-1
		do icc=1,ncc
		if(cloudfrac(icc).gt.0d0) then
			do i=1,nrtrace-1
				A=1d0
				do imol=1,nmol
					if(includemol(imol)) then
					Ag=0d0
					do ig=1,ng
						tautot=0d0
						do k=1,nirtrace(i)
							ir=irtrace(i,k)
							d=dtrace(i,k)
							tau=d*Cabs_mol(imol,ir,ilam,ig)*mixrat_r(ir,imol)
							if(P(ir).gt.Psimplecloud) tau=1d4
							tautot=tautot+tau
							if(k.lt.nirtrace(i)) ir_next=irtrace(i,k+1)
							if(ir_next.le.0.or.tautot.ge.maxtau) then
								tautot=1d4
								exit
							endif
						enddo
						Ag=Ag+exp(-tautot)*wgg(ig)
					enddo
					A=A*Ag
					endif
				enddo
				tautot=0d0
				do k=1,nirtrace(i)
					ir=irtrace(i,k)
					d=dtrace(i,k)
					Ca=Cext_cont(ir,ilam)
					do icloud=1,nclouds
						if(docloud(icc,icloud)) then
							if(useDRIFT) then
								Ca=Ca+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
								Ca=Ca+Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
							else
								do isize=1,Cloud(icloud)%nr
									Ca=Ca+
     &		Cloud(icloud)%Kabs(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
									Ca=Ca+
     &		Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
								enddo
							endif
						endif
					enddo
					tau=d*Ca
					if(P(ir).gt.Psimplecloud) tau=1d4
					tautot=tautot+tau
					if(k.lt.nirtrace(i)) ir_next=irtrace(i,k+1)
					if(ir_next.le.0.or.tautot.ge.maxtau) then
						tautot=1d4
						exit
					endif
				enddo
				A=A*exp(-tautot)
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
	deallocate(ptrace)

	return
	end
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	real*8 function Planck(T,nu)
	use Constants
	IMPLICIT NONE
	real*8 T,nu,x

	x=hplanck*nu*clight/(kb*T)
	Planck=(2d0*hplanck*nu**3*clight)/(exp(x)-1d0)

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

