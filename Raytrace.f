	subroutine Raytrace(iobs)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs
	real*8 rr,xx1,xx2,si,exp_tau,A,d,s,fluxg,Planck,fact,tau,freq0,tau_a,tautot,Ag
	real*8 Ca,Cs,BB(nr)
	integer icloud,isize
	real*8,allocatable :: rtrace(:),phase(:)
	integer nrtrace,ndisk,i,ir,ir_next,ilam,ig,nsub,j,k
	logical in
	integer icc

	allocate(phase(obs(iobs)%nphase))

	obs(iobs)%docloud=.false.
	do icc=2,obs(iobs)%ncc
		obs(iobs)%docloud(icc,1:nclouds)=obs(iobs)%docloud(icc-1,1:nclouds)
		i=0
10		i=i+1
		obs(iobs)%docloud(icc,i)=.not.obs(iobs)%docloud(icc,i)
		if(.not.obs(iobs)%docloud(icc,i)) goto 10
	enddo
	do icc=1,obs(iobs)%ncc
		obs(iobs)%cloudfrac(icc)=1d0
		do i=1,nclouds
			if(obs(iobs)%docloud(icc,i)) then
				obs(iobs)%cloudfrac(icc)=obs(iobs)%cloudfrac(icc)*Cloud(i)%coverage
			else
				obs(iobs)%cloudfrac(icc)=obs(iobs)%cloudfrac(icc)*(1d0-Cloud(i)%coverage)
			endif
		enddo
	enddo

	obs(iobs)%flux=0d0

	call output("==================================================================")

	if(scattering) then
	call output("Scattered light contributions")

	call tellertje(1,nlam)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,fluxg,icc,phase)
!$OMP& SHARED(nlam,nclouds,obs,iobs)
!$OMP DO
	do ilam=1,nlam-1
		call tellertje(ilam+1,nlam+1)
		obs(iobs)%flux(:,ilam)=0d0
		obs(iobs)%phase(:,:,ilam)=0d0
		do icc=1,obs(iobs)%ncc
			call MCRad(ilam,fluxg,phase,obs(iobs)%nphase,obs(iobs)%docloud(icc,1:nclouds))
			obs(iobs)%flux(0,ilam)=obs(iobs)%flux(0,ilam)+obs(iobs)%cloudfrac(icc)*fluxg
			obs(iobs)%flux(icc,ilam)=obs(iobs)%flux(icc,ilam)+fluxg
			obs(iobs)%phase(1:obs(iobs)%nphase,icc,ilam)=obs(iobs)%phase(1:obs(iobs)%nphase,icc,ilam)+
     &				phase(1:obs(iobs)%nphase)
			obs(iobs)%phase(1:obs(iobs)%nphase,0,ilam)=obs(iobs)%phase(1:obs(iobs)%nphase,0,ilam)+
     &				obs(iobs)%cloudfrac(icc)*phase(1:obs(iobs)%nphase)
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(nlam,nlam)
	endif

	call output("Raytracing over the planet disk")

	ndisk=10
	nsub=3
	
	nrtrace=nr*nsub+ndisk
	allocate(rtrace(nrtrace))

	k=0
	do i=1,ndisk
		k=k+1
		rtrace(k)=Rplanet*real(i-1)/real(ndisk-1)
	enddo
	do i=1,nr
		do j=1,nsub
			k=k+1
			rtrace(k)=R(i)+(R(i+1)-R(i))*real(j)/real(nsub)
		enddo
	enddo

	call tellertje(1,nlam)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,freq0,ig,i,fluxg,fact,A,rr,ir,si,xx1,in,xx2,d,ir_next,tau,exp_tau,tau_a,tautot,Ag,
!$OMP&         Ca,Cs,icloud,isize,BB)
!$OMP& SHARED(nlam,freq,obs,iobs,nrtrace,ng,rtrace,nr,R,Ndens,Cabs,Csca,T,lam,maxtau,nclouds,Cloud,
!$OMP&			cloud_dens)
!$OMP DO SCHEDULE(STATIC,1)
	do ilam=1,nlam-1
		call tellertje(ilam+1,nlam+1)
		freq0=sqrt(freq(ilam)*freq(ilam+1))
		obs(iobs)%lam(ilam)=sqrt(lam(ilam)*lam(ilam+1))
		obs(iobs)%A(:,ilam)=0d0
		do ir=1,nr
			BB(ir)=Planck(T(ir),freq0)
		enddo
		do ig=1,ng
			do icc=1,obs(iobs)%ncc
			if(obs(iobs)%cloudfrac(icc).gt.0d0) then
				fluxg=0d0
				Ag=0d0
				do i=1,nrtrace-1
					fact=1d0
					tautot=0d0
					A=pi*(rtrace(i+1)**2-rtrace(i)**2)
					rr=sqrt(rtrace(i)*rtrace(i+1))
					ir=nr
					si=-1d0
					xx1=si*sqrt(R(ir+1)**2-rr**2)
					in=.true.
1					continue
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
2					continue
					Ca=Cabs(ir,ilam,ig)*Ndens(ir)
					Cs=Csca(ir,ilam)*Ndens(ir)
					do icloud=1,nclouds
						if(obs(iobs)%docloud(icc,icloud)) then
							do isize=1,Cloud(icloud)%nsize
								Ca=Ca+
     &		Cloud(icloud)%Kabs(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
								Cs=Cs+
     &		Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
							enddo
						endif
					enddo
					tau_a=d*Ca
					tau=tau_a+d*Cs
					exp_tau=exp(-tau)
					tautot=tautot+tau
					fluxg=fluxg+A*BB(ir)*(1d0-exp_tau)*fact*tau_a/tau
					fact=fact*exp_tau
					if(ir_next.gt.0.and.ir_next.le.nr.and.tautot.lt.maxtau) then
						ir=ir_next
						xx1=xx2
						goto 1
					endif
					if(ir_next.le.0.or.tautot.ge.maxtau) fact=0d0
					Ag=Ag+A*(1d0-fact)
				enddo
				obs(iobs)%flux(0,ilam)=obs(iobs)%flux(0,ilam)+obs(iobs)%cloudfrac(icc)*fluxg/real(ng)
				obs(iobs)%A(0,ilam)=obs(iobs)%A(0,ilam)+obs(iobs)%cloudfrac(icc)*Ag/real(ng)
				obs(iobs)%flux(icc,ilam)=obs(iobs)%flux(icc,ilam)+fluxg/real(ng)
				obs(iobs)%A(icc,ilam)=obs(iobs)%A(icc,ilam)+Ag/real(ng)
			endif
			enddo
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(nlam,nlam)
	
	do ilam=1,nlam
		do icc=0,obs(iobs)%ncc
			obs(iobs)%A(icc,ilam)=obs(iobs)%A(icc,ilam)-
     &			obs(iobs)%phase(obs(iobs)%nphase,icc,ilam)/(Fstar(ilam)/(pi*Rstar**2))
		enddo
	enddo

	obs(iobs)%flux=obs(iobs)%flux*1d23/distance**2
	obs(iobs)%phase=obs(iobs)%phase*1d23/distance**2

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
	Planck=(2d0*hplanck*nu**3*clight)/(exp(x)-1d0)

	return
	end

