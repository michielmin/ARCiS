	module Struct3D
	implicit none
	integer nlong,nlatt
	integer n3D,nnu0
	parameter(n3D=10,nnu0=10)
	parameter(nlong=36,nlatt=18)
	real*8 long(nlong),latt(nlatt)	!(Lambda, Phi)
	real*8 tanx(nlong),tany(nlong)
	real*8 cost2(nlatt)
	real*8,allocatable :: R3D(:,:),R3D2(:,:)
	integer ibeta(nlong,nlatt),inu3D(nlong,nlatt)
	end module Struct3D

	subroutine Run3D(recomputeopac)
	use GlobalSetup
	use Constants
	use Struct3D
	IMPLICIT NONE
	real*8 beta3D(n3D),Kzz3D(n3D),Sdot3D(n3D)
	integer i,j,icloud,k,irtrace,iptrace,inu,imol
	real*8 beta(nlong,nlatt),Planck,phi,la,lo,A,rr
	real*8 long0,b1,b2,betamin,betamax,freq0,Rmax,theta
	real*8,allocatable :: Ca(:,:,:,:),Cs(:,:,:),BBr(:,:,:),Si(:,:,:,:,:),Ca_mol(:,:,:,:,:),Ce(:,:,:)
	integer ir,ilam,ig,isize,iRmax,ndisk,nsub,nrtrace,nptrace,ipc,npc,nmol_count
	logical recomputeopac
	logical,allocatable :: hit(:,:)
	real*8,allocatable :: rtrace(:),wrtrace(:)
	real*8 x,y,z,vx,vy,vz,v
	integer edgeNR,i1,i2,i3,i1next,i2next,i3next,edgenext
	real*8,allocatable :: fluxp(:),tau(:,:),fact(:,:),tautot(:,:),exp_tau(:,:),tauR(:,:,:),SiR(:,:,:),Ij(:),tauR1(:),SiR1(:)
	real*8,allocatable :: tauc(:),Afact(:),vv(:,:,:),obsA_omp(:),mixrat3D(:,:,:),T3D(:,:),fluxp_omp(:)
	type(Mueller) M
	real*8 g,tot,contr,tmp(nmol)
	character*500 file
	real*8 tau1,fact1,exp_tau1

	allocate(Ca(nlam,ng,nr,n3D),Cs(nlam,nr,n3D),BBr(nlam,nr,n3D),Si(nlam,ng,nr,nnu0,n3D))
	allocate(Ca_mol(nlam,ng,nmol,nr,n3D),Ce(nlam,nr,n3D))
	allocate(R3D(n3D,nr+2))
	allocate(R3D2(n3D,nr+2))
	allocate(T3D(n3D,nr))
	allocate(mixrat3D(n3D,nr,nmol))

	if(retrieval) call SetOutputMode(.false.)

c	recomputeopac=.true.
	docloud=.true.
	cloudfrac=1d0

	if(beta3D_1.ge.0d0) then
		b1=beta3D_1	! evening limb, limited to 0-0.25
		b2=beta3D_2	! morning limb, limited to 0-0.25
	else
		b1=betaT
		b2=betaT
	endif
	long0=long_shift*pi/180d0

c	call Setup3D_old(beta,long,latt,nlong,nlatt,long0,b1,b2,betapow,fDay,betamin,betamax)
	call Setup3D(beta,long,latt,nlong,nlatt,Kxx,vxx,night2day,fDay,betamin,betamax)
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
	do i=1,nlong-1
		do j=1,nlatt-1
			if(betamax.eq.betamin) then
				ibeta(i,j)=1
			else
				ibeta(i,j)=(real(n3D-1)*((beta(i,j)-betamin)/(betamax-betamin))+0.5d0)+1
				if(ibeta(i,j).gt.n3D) ibeta(i,j)=n3D
				if(.not.ibeta(i,j).gt.1) ibeta(i,j)=1
			endif
			la=0.5d0*(latt(j)+latt(j+1))-pi/2d0
			lo=0.5d0*(long(i)+long(i+1))-pi
			x=cos(lo)*cos(la)
			if(x.lt.0d0) then
				inu3D(i,j)=nnu0
			else
				inu3D(i,j)=(real(nnu0-2)*x+0.5d0)+1
				if(inu3D(i,j).lt.1) inu3D(i,j)=1
				if(.not.inu3D(i,j).lt.nnu0) inu3D(i,j)=nnu0
			endif				
		enddo
	enddo
	Rmax=0d0
	call output("Computing multiple 1D structures")

	call tellertje_perc(0,n3D)
	do i=1,n3D
		call SetOutputMode(.false.)
		beta3D(i)=betamin+(betamax-betamin)*real(i-1)/real(n3D-1)
		do j=1,n_Par3D
			if(Par3D(j)%logscale) then
c				Par3D(j)%x=10d0**(log10(Par3D(j)%xmin)+log10(Par3D(j)%xmax/Par3D(j)%xmin)*beta3D(i))
				Par3D(j)%x=10d0**(log10(Par3D(j)%xmin)+
     &							  log10(Par3D(j)%xmax/Par3D(j)%xmin)*(real(i-1)/real(n3D-1))**Par3D(j)%pow)
			else
c				Par3D(j)%x=Par3D(j)%xmin+(Par3D(j)%xmax-Par3D(j)%xmin)*beta3D(i)
				Par3D(j)%x=Par3D(j)%xmin+(Par3D(j)%xmax-Par3D(j)%xmin)*(real(i-1)/real(n3D-1))**Par3D(j)%pow
			endif
		enddo
		call MapPar3D()

		betaT=beta3D(i)

		call InitDens()
		call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam,starfile)
		Fstar=Fstar*pi*Rstar**2
c===============================================================
c quick thing to read in a file!
c	file='houghtonsolarwl.dat'
c	call regridlog(file,1d4*lam,Fstar,nlam)
c	Fstar=Fstar*lam**2/4d0
c===============================================================
		call ComputeModel1D(recomputeopac)

		if(R(nr+1).gt.Rmax) then
			Rmax=R(nr+1)
			iRmax=i
		endif
		do ir=1,nr+1
			R3D(i,ir)=R(ir)
			R3D2(i,ir)=R(ir)**2
		enddo
		do ir=1,nr
			T3D(i,ir)=T(ir)
			mixrat3D(i,ir,1:nmol)=mixrat_r(ir,1:nmol)
			do ilam=1,nlam-1
				if(useobsgrid) then
					freq0=freq(ilam)
				else
					freq0=sqrt(freq(ilam)*freq(ilam+1))
				endif
				BBr(ilam,ir,i)=Planck(T(ir),freq0)
				Ca(ilam,1:ng,ir,i)=0d0
				Cs(ilam,ir,i)=0d0
				do icloud=1,nclouds
					if(Cloud(icloud)%standard.eq.'MIX') then
						Ca(ilam,1:ng,ir,i)=Ca(ilam,1:ng,ir,i)+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
						Cs(ilam,ir,i)=Cs(ilam,ir,i)+Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
					else
						do isize=1,Cloud(icloud)%nr
							Ca(ilam,1:ng,ir,i)=Ca(ilam,1:ng,ir,i)+
     &		Cloud(icloud)%Kabs(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
							Cs(ilam,ir,i)=Cs(ilam,ir,i)+
     &		Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
						enddo
					endif
				enddo
				Cs(ilam,ir,i)=Cs(ilam,ir,i)+Csca(ir,ilam)*Ndens(ir)

c ===================================================================
c correction for anisotropic scattering	
c ===================================================================
 			if(scattering.and..false.) then
 			M%F11=Rayleigh%F11*Csca(ir,ilam)*Ndens(ir)
 			do icloud=1,nclouds
 				if(Cloud(icloud)%standard.eq.'MIX') then
 					M%F11=M%F11+Cloud(icloud)%F(ir,ilam)%F11*Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
 				else
 					do isize=1,Cloud(icloud)%nr
 						M%F11=M%F11+Cloud(icloud)%F(isize,ilam)%F11*
     &								Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
 					enddo
 				endif
 			enddo
 			g=0d0
 			tot=0d0
 			do j=1,180
 				g=g+M%F11(j)*costheta(j)*sintheta(j)
 				tot=tot+M%F11(j)*sintheta(j)
 			enddo
 			g=g/tot
 			if(.not.g.ge.-1d0) then
 				g=0.999d0
 			endif
 			if(.not.g.le.1d0) then
 				g=0.999d0
 			endif
 			Cs(ilam,ir,i)=Cs(ilam,ir,i)*(1d0-g)
 		endif
c ===================================================================
c ===================================================================

				Ce(ilam,ir,i)=Ca(ilam,1,ir,i)+Cs(ilam,ir,i)
				do ig=1,ng
					Ca(ilam,ig,ir,i)=Ca(ilam,ig,ir,i)+Cabs(ir,ilam,ig)*Ndens(ir)
				enddo
			enddo
		enddo
		do ir=1,nr
			k=0
			do imol=1,nmol
				if(includemol(imol)) then
					k=k+1
					do ig=1,ng
						do ilam=1,nlam
							Ca_mol(ilam,ig,k,ir,i)=Cabs_mol(ir,ig,imol,ilam)
						enddo
					enddo
				endif
			enddo
			nmol_count=k
		enddo
		if(i.eq.1.or.betamax.ne.betamin) then
			if(emisspec) call ComputeScatter(BBr(1:nlam,1:nr,i),Si(1:nlam,1:ng,1:nr,1:nnu0,i),Ca(1:nlam,1:ng,1:nr,i),Cs(1:nlam,1:nr,i))
		else
			Si(1:nlam,1:ng,1:nr,1:nnu0,i)=Si(1:nlam,1:ng,1:nr,1:nnu0,1)
		endif
		if(.not.retrieval) call SetOutputMode(.true.)
		call tellertje_perc(i,n3D)
	enddo
	Rmax=Rmax*1.001
	R3D(1:n3D,nr+2)=Rmax
	R3D2(1:n3D,nr+2)=Rmax**2

	call output("Raytracing over the planet disk in 3D")

	do ilam=1,nlam
		cloudtau(1,ilam)=0d0
		do ir=nr,1,-1
			do icloud=1,nclouds
				if(Cloud(icloud)%standard.eq.'MIX') then
					cloudtau(1,ilam)=cloudtau(1,ilam)+(R(ir+1)-R(ir))*Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
					cloudtau(1,ilam)=cloudtau(1,ilam)+(R(ir+1)-R(ir))*Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
				else
					do isize=1,Cloud(icloud)%nr
						cloudtau(1,ilam)=cloudtau(1,ilam)+
     &		(R(ir+1)-R(ir))*Cloud(icloud)%Kabs(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
						cloudtau(1,ilam)=cloudtau(1,ilam)+
     &		(R(ir+1)-R(ir))*Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
					enddo
				endif
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
	nptrace=nlatt
	allocate(rtrace(nrtrace),wrtrace(nrtrace))

	call gauleg(0d0,Rmax,rtrace,wrtrace,ndisk)

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
!$OMP& PRIVATE(irtrace,iptrace,A,phi,rr,y,z,x,vx,vy,vz,la,lo,i1,i2,i3,edgeNR,j,i,inu,fluxp_omp,
!$OMP&			i1next,i2next,i3next,edgenext,freq0,tot,v,ig,ilam,Ij,tau1,fact,exp_tau1,contr)
!$OMP& SHARED(theta,fluxp,nrtrace,rtrace,wrtrace,nptrace,Rmax,nr,useobsgrid,freq,ibeta,inu3D,fulloutput3D,
!$OMP&			Ca,Cs,wgg,Si,R3D2,latt,long,T,ng,nlam,ipc,PTaverage3D,mixrat_average3D,T3D,mixrat3D,nmol)
	allocate(fact(nlam,ng))
	allocate(fluxp_omp(nlam))
	fluxp_omp=0d0
!$OMP DO
	do irtrace=1,nrtrace
		A=2d0*pi*rtrace(irtrace)*wrtrace(irtrace)/real(nptrace)
		do iptrace=1,nptrace
			fact=1d0
c Note we are here using the symmetry between North and South
			phi=pi*(real(iptrace)-0.5)/real(nptrace)
			rr=rtrace(irtrace)
			y=rr*cos(phi)
			z=rr*sin(phi)
			x=sqrt(Rmax**2-y**2-z**2)
			vx=-1d0
			vy=0d0
			vz=0d0
			call rotateZ3D(x,y,z,theta)
			rr=sqrt(x**2+y**2+z**2)
			x=x*Rmax/rr
			y=y*Rmax/rr
			z=z*Rmax/rr
			call rotateZ3D(vx,vy,vz,theta)
			rr=sqrt(vx**2+vy**2+vz**2)
			vx=vx/rr
			vy=vy/rr
			vz=vz/rr
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
			i=ibeta(i2,i3)
			inu=inu3D(i2,i3)
			if(fulloutput3D) then
				PTaverage3D(ipc,1:nr)=PTaverage3D(ipc,1:nr)+T3D(i,1:nr)*A
				mixrat_average3D(ipc,1:nr,1:nmol)=mixrat_average3D(ipc,1:nr,1:nmol)+mixrat3D(i,1:nr,1:nmol)*A
			endif
1			continue
			call TravelSph(x,y,z,vx,vy,vz,edgeNR,i1,i2,i3,v,i1next,i2next,i3next,edgenext,nr,nlong,nlatt)
			if(i1next.le.0) then
				j=j+1
				do ilam=1,nlam-1
					if(useobsgrid) then
						freq0=freq(ilam)
					else
						freq0=sqrt(freq(ilam)*freq(ilam+1))
					endif
					contr=Planck(T(1),freq0)
					do ig=1,ng
						fluxp_omp(ilam)=fluxp_omp(ilam)+A*wgg(ig)*contr*fact(ilam,ig)
					enddo
				enddo
				goto 2
			endif
			if(i1next.ge.nr+2) goto 2
			if(i1.le.nr) then
				i=ibeta(i2,i3)
				inu=inu3D(i2,i3)
				do ig=1,ng
					do ilam=1,nlam
						tau1=v*(Ca(ilam,ig,i1,i)+Cs(ilam,i1,i))
						exp_tau1=exp(-tau1)
						fluxp_omp(ilam)=fluxp_omp(ilam)+A*wgg(ig)*Si(ilam,ig,i1,inu,i)*(1d0-exp_tau1)*fact(ilam,ig)
						fact(ilam,ig)=fact(ilam,ig)*exp_tau1
					enddo
				enddo
			endif
			x=x+v*vx
			y=y+v*vy
			z=z+v*vz
			if(ibeta(i2,i3).ne.ibeta(i2next,i3next)) then
				rr=x**2+y**2+z**2
				i=ibeta(i2next,i3next)
				if(rr.lt.R3D2(i,1).or.rr.gt.R3D2(i,nr+2)) goto 2
				do i1=1,nr+1
					if(rr.gt.R3D2(i,i1).and.rr.le.R3D2(i,i1+1)) exit
				enddo
			else
				i1=i1next
			endif
			i2=i2next
			i3=i3next
			edgeNR=edgenext
			goto 1
2			continue
		enddo
	enddo
!$OMP END DO
!$OMP CRITICAL
	fluxp(1:nlam)=fluxp(1:nlam)+fluxp_omp(1:nlam)
	deallocate(fluxp_omp)
!$OMP END CRITICAL
	deallocate(fact)
!$OMP FLUSH
!$OMP END PARALLEL
	fluxp=fluxp*1d23/distance**2
	phase(ipc,0,1:nlam)=fluxp(1:nlam)
	flux(0,1:nlam)=0d0
	call tellertje_perc(ipc,npc)
	if(fulloutput3D) then
		PTaverage3D(ipc,1:nr)=PTaverage3D(ipc,1:nr)/(pi*Rmax**2)
		mixrat_average3D(ipc,1:nr,1:nmol)=mixrat_average3D(ipc,1:nr,1:nmol)/(pi*Rmax**2)
		open(unit=25,file=trim(outputdir) // "mixrat_phase" // trim(int2string(int(theta_phase(ipc)),'(i0.3)')),RECL=6000)
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
	enddo
	
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
		open(unit=25,file=trim(outputdir) // "mixrat_transit",RECL=6000)
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

	ndisk=4
	nsub=3

	nrtrace=(nr-1)*nsub+ndisk
	nptrace=nlatt
	allocate(rtrace(nrtrace))

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
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(vv,tau,tauc,Afact,hit,irtrace,A,iptrace,phi,rr,x,y,z,vx,vy,vz,la,lo,i1,i2,i3,
!$OMP&			edgeNR,v,i1next,i2next,i3next,edgenext,i,imol,ir,ig,obsA_omp)
!$OMP& SHARED(nrtrace,nptrace,nmol_count,nr,nlam,ng,rtrace,ibeta,R3D2,Ce,Ca_mol,wgg,obsA,latt,long,Rmax)
	allocate(obsA_omp(nlam))
	allocate(vv(nmol_count,nr,n3D))
	allocate(tau(nlam,ng))
	allocate(tauc(nlam),Afact(nlam))
	allocate(hit(nr+2,n3D))
	obsA_omp=0d0
!$OMP DO
	do irtrace=1,nrtrace-1
		A=pi*(rtrace(irtrace+1)**2-rtrace(irtrace)**2)/real(nptrace)
		do iptrace=1,nptrace
c Note we use the symmetry of the North and South here!
			phi=pi*(real(iptrace)-0.5)/real(nptrace)
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
			if(i1next.le.0.or.i1next.ge.nr+2) goto 4
			if(i1.le.nr) then
				i=ibeta(i2,i3)
				vv(1:nmol_count,i1,i)=vv(1:nmol_count,i1,i)+v
				hit(i1,i)=.true.
				tauc(1:nlam)=tauc(1:nlam)+v*Ce(1:nlam,i1,i)
			endif
			x=x+v*vx
			y=y+v*vy
			z=z+v*vz
			if(ibeta(i2,i3).ne.ibeta(i2next,i3next)) then
				rr=x**2+y**2+z**2
				i=ibeta(i2next,i3next)
				if(rr.lt.R3D2(i,1).or.rr.gt.R3D2(i,nr+2)) goto 4
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
							tau(1:nlam,1:ng)=tau(1:nlam,1:ng)+vv(imol,ir,i)*Ca_mol(1:nlam,1:ng,imol,ir,i)
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
		enddo
	enddo
!$OMP END DO
!$OMP CRITICAL
	obsA(0,1:nlam)=obsA(0,1:nlam)+obsA_omp(1:nlam)
!$OMP END CRITICAL
	deallocate(obsA_omp)
	deallocate(vv)
	deallocate(tau)
	deallocate(tauc,Afact)
	deallocate(hit)
!$OMP FLUSH
!$OMP END PARALLEL
	deallocate(rtrace)

	endif
	
	deallocate(Ca,Cs,BBr,Si)
	deallocate(Ca_mol,Ce)
	deallocate(R3D)
	deallocate(R3D2)
	deallocate(T3D,mixrat3D)

	if(retrieval) call SetOutputMode(.true.)
	
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
	
	Rplanet=Rplanet/Rjup
	Mplanet=Mplanet/Mjup
	Rstar=Rstar/Rsun
	Mstar=Mstar/Msun
	Dplanet=Dplanet/AU
	lam1=lam1/micron
	lam2=lam2/micron
	distance=distance/parsec
	r_nuc=r_nuc/micron
	orbit_inc=orbit_inc*180d0/pi

	metallicity=metallicity0
	do i=1,n_Par3D
		readline=trim(Par3D(i)%keyword) // "=" // trim(dbl2string(Par3D(i)%x,'(es14.7)'))
		call get_key_value(readline,key%key,key%key1,key%key2,key%value,key%nr1,key%nr2,key%key2d)
		call ReadAndSetKey(key)
	enddo
	call ConvertUnits()
	metallicity0=metallicity

	return
	end
	

	subroutine Setup3D_old(beta,long,latt,nlong,nlatt,long0,b1,b2,p,f,betamin,betamax)
	IMPLICIT NONE
	integer i,j,nlong,nlatt
	real*8 pi
	parameter(pi=3.1415926536)
	real*8 long(nlong),latt(nlatt)	!(Lambda, Phi)
	real*8 beta(nlong,nlatt),la,lo,albedo,p
	real*8 long0,b1,b2,f,firr1,firr2,betamin,betamax

	do i=1,nlong
		long(i)=-pi+2d0*pi*real(i-1)/real(nlong-1)
	enddo
	do j=1,nlatt
		latt(j)=-pi/2d0+pi*real(j-1)/real(nlatt-1)
	enddo
	
	betamin=1d0
	betamax=0d0
	do i=1,nlong-1
		do j=1,nlatt-1
			lo=(long(i)+long(i+1))/2d0
			la=(latt(j)+latt(j+1))/2d0
			if(abs(lo).lt.pi/2d0) then	! dayside
				firr1=b1+(b2-b1)*((lo+pi/2d0)/pi)**p
			else if(lo.lt.(-pi/2d0)) then
				firr1=b1+(b2-b1)*((-lo-pi/2d0)/pi)**p
			else
				firr1=b1+(b2-b1)*(1d0-(lo-pi/2d0)/pi)**p
			endif
			firr2=0d0
			if(lo.gt.(-pi/2d0).and.lo.lt.(long0)) then
				firr2=f*cos((lo-long0)*(pi/2d0)/(pi/2d0+long0))
			else if(lo.lt.(pi/2d0).and.lo.gt.(long0)) then
				firr2=f*cos((lo-long0)*(pi/2d0)/(pi/2d0-long0))
			endif
			beta(i,j)=abs((b1+(b2-b1)/(p+1d0))
     &		+(firr1+firr2-(b1+(b2-b1)/(p+1d0)))*cos(la))
			if(beta(i,j).gt.betamax) betamax=beta(i,j)
			if(beta(i,j).lt.betamin) betamin=beta(i,j)
		enddo
	enddo
	
	latt=latt+pi/2d0
	long=long+pi

	return
	end
	

	subroutine Setup3D(beta,long,latt,nlong,nlatt,Kxx,vxx,night2day,f,betamin,betamax)
	IMPLICIT NONE
	integer i,j,nlong,nlatt
	real*8 pi
	parameter(pi=3.1415926536)
	real*8 long(nlong),latt(nlatt)	!(Lambda, Phi)
	real*8 beta(nlong,nlatt),la,lo,albedo,p,f,tot1,tot2,b1,b2
	real*8 Kxx,vxx,betamin,betamax,beta1(nlong),night2day

	do i=1,nlong
		long(i)=-pi+2d0*pi*real(i-1)/real(nlong-1)
	enddo
	do j=1,nlatt
		latt(j)=-pi/2d0+pi*real(j-1)/real(nlatt-1)
	enddo
	
	betamin=1d0
	betamax=0d0

	call DiffuseBeta(beta,long,latt,nlong,nlatt,Kxx,vxx,night2day)
	beta=beta*f

	do i=1,nlong-1
		do j=1,nlatt-1
			if(beta(i,j).gt.betamax) betamax=beta(i,j)
			if(beta(i,j).lt.betamin) betamin=beta(i,j)
		enddo
	enddo
	
	latt=latt+pi/2d0
	long=long+pi

	open(unit=20,file='output.dat',RECL=6000)
	do j=1,nlatt-1
		write(20,*) beta(1:nlong-1,j)
	enddo
	close(unit=20)

	return
	end
	

	subroutine DiffuseBeta(beta,long,latt,nlong,nlatt,Kxx,vxx,contrast)
	IMPLICIT NONE
	integer nlatt,nlong,NN
	integer,allocatable :: IWORK(:)
	integer i,info,j,NRHS,ii(nlong-1,nlatt-1)
	real*8 beta(nlong,nlatt),Kxx,vxx,long(nlong),latt(nlatt),S(nlong-1,nlatt-1)
	real*8 pi,tot,x((nlatt-1)*(nlong-1)),contrast,eps
	parameter(eps=1d-3)
	real*8 la(nlatt-1),lo(nlong-1),tp,tm,tot1,tot2
	real*8,allocatable :: A(:,:)
	integer j0,jm,jp,k,niter,maxiter
	parameter(pi=3.1415926536)
	real*8 smax,smin,scale,betamin,betamax,contr

	allocate(IWORK(50*nlatt*nlong*nlatt*nlong))
	allocate(A((nlatt-1)*(nlong-1),(nlatt-1)*(nlong-1)))

	smax=-1d0
	smin=-1d0
	scale=1d0
	
	contr=contrast*10d0
	maxiter=10000
	niter=0
	do while(abs((contrast-contr)/(contrast+contr)).gt.eps.or.niter.gt.maxiter)
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
			A(j0,jp)=A(j0,jp)+Kxx*scale*real(nlatt*2)**2*cos(tp)/cos(la(j))
			A(j0,jm)=A(j0,jm)+Kxx*scale*real(nlatt*2)**2*cos(tm)/cos(la(j))
			A(j0,j0)=A(j0,j0)-Kxx*scale*real(nlatt*2)**2*(cos(tp)+cos(tm))/cos(la(j))
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
		if(smax.lt.0d0) then
			scale=scale*2d0
		else
			scale=(scale+smax)/2d0
		endif
	else
		smax=scale
		if(smin.lt.0d0) then
			scale=scale/2d0
		else
			scale=(scale+smin)/2d0
		endif
	endif
	enddo
	
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
	parameter(pi=3.1415926536)
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
	
	return
	end
	



	subroutine DiffuseBeta1D(beta,theta,nbeta,K,v)
	IMPLICIT NONE
	integer nbeta
	integer i,info,IWORK(50*nbeta*nbeta),j,NRHS
	real*8 beta(nbeta),K,v,theta(nbeta),S(nbeta)
	real*8 pi,A(nbeta,nbeta),tot
	parameter(pi=3.1415926536)
	
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
	if(v.lt.0d0) hitRin=.false.

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
	if(v.lt.0d0) hitRout=.false.
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
	if(v.le.0d0) hitTsame=.false.

	return
	end

	logical function hitP(tanx,tany,x0,vx,y0,vy,v)
	IMPLICIT NONE
	real*8 tanx,tany,x0,vx,y0,vy,v
	
	hitP=.true.
	v=(tany*y0-tanx*x0)/(tanx*vx-tany*vy)
	
	if(v.lt.0d0) hitP=.false.
	
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

	i1midplane=0!(i3.eq.(nlatt+1)/2)
	i2midplane=0!((i3+1).eq.(nlatt+1)/2)

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

	v=1d200
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
		if(i3next.lt.1) i3next=1
		edgenext=4
	endif
	if(hitT2.and.vT2.lt.v.and.vT2.gt.0d0) then
		v=vT2
		i1next=i1
		i2next=i2
		i3next=i3+1
		if(i3next.gt.n3) i3next=n3
		edgenext=3
	endif
	if(hitP1.and.vP1.lt.v.and.vP1.gt.0d0) then
		v=vP1
		i1next=i1
		i2next=i2-1
		i3next=i3
		if(i2next.lt.1) i2next=nlong-1
		edgenext=6
	endif
	if(hitP2.and.vP2.lt.v.and.vP2.gt.0d0) then
		v=vP2
		i1next=i1
		i2next=i2+1
		i3next=i3
		if(i2next.ge.nlong) i2next=1
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
	parameter(nnu=5,niter=500)
	real*8 tau,d,tauR_nu(nr,nlam,ng),contr,Jstar_nu(nr,nlam,ng)
	real*8 Si(nlam,ng,nr,nnu0),BBr(nlam,nr),Ca(nlam,ng,nr),Cs(nlam,nr),Ce(nlam,ng,nr)
	real*8 nu(nnu),wnu(nnu),must,tauRs(nr),Ijs(nr),eps,Planck,tot
	logical err
	parameter(eps=1d-20)
	real*8,allocatable :: tauR(:),Ij(:),Itot(:),Linv(:,:),Lmat(:,:),Iprev(:)
	integer,allocatable :: IWORKomp(:)

	if(.not.scattering) then
		do inu0=1,nnu0
			do ig=1,ng
				Si(1:nlam,ig,1:nr,inu0)=BBr(1:nlam,1:nr)*Ca(1:nlam,ig,1:nr)/(Ca(1:nlam,ig,1:nr)+Cs(1:nlam,1:nr))
			enddo
		enddo
		return
	endif

	do ir=1,nr
		do ig=1,ng
			do ilam=1,nlam
				Ce(ilam,ig,ir)=Ca(ilam,ig,ir)+Cs(ilam,ir)
			enddo
		enddo
	enddo

	tauR_nu=0d0
	do ilam=1,nlam-1
		do ig=1,ng
			do ir=nr,1,-1
c				d=abs(P(ir+1)-P(ir))*1d6/grav(ir)
c				tau=d*Ce(ilam,ig,ir)/dens(ir)
				d=R(ir+1)-R(ir)
				tau=d*Ce(ilam,ig,ir)
				if(P(ir).gt.Psimplecloud) then
					tau=tau+1d4
				endif
				if(.not.tau.gt.1d-10) then
					tau=1d-10
				endif
				if(tau.gt.1d10) then
					tau=1d10
				endif
				if(ir.eq.nr) then
					tauR_nu(ir,ilam,ig)=tau
				else
					tauR_nu(ir,ilam,ig)=tauR_nu(ir+1,ilam,ig)+tau
				endif
			enddo
		enddo
	enddo

	call gauleg(0d0,1d0,nu,wnu,nnu)

	do inu0=1,nnu0
		if(inu0.eq.nnu0.or..not.scattstar) then
			Jstar_nu=0d0
			if(inu0.ne.1.and..not.scattstar) then
				Si(1:nlam,1:ng,1:nr,inu0)=Si(1:nlam,1:ng,1:nr,1)
				goto 1
			endif
		else
			must=(real(inu0)-0.5)/real(nnu0-1)
			do ilam=1,nlam-1
				do ig=1,ng
					tauRs(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(must)
					contr=(Fstar(ilam)/(4d0*pi*Dplanet**2))
					call SolveIjStar(tauRs,contr,Ijs,nr)
					Jstar_nu(1:nr,ilam,ig)=Ijs(1:nr)*2d0
				enddo
			enddo
		endif
	
		NRHS=1
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,ig,Linv,Lmat,tauR,Itot,IWORKomp,inu,ir,info)
!$OMP& SHARED(nlam,ng,lamemis,tauR_nu,nr,wnu,nu,BBr,Ca,Cs,Ce,Si,inu0,NRHS,Jstar_nu)
		allocate(tauR(nr))
		allocate(Itot(nr))
		allocate(Linv(nr,nr))
		allocate(Lmat(nr,nr))
		allocate(IWORKomp(10*nr*nr))
!$OMP DO
		do ilam=1,nlam-1
			if(lamemis(ilam)) then
				do ig=1,ng
					Linv=0d0
					do inu=1,nnu
						tauR(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(nu(inu))
						call InvertIj(tauR,Lmat,nr)
						do ir=1,nr
							Linv(ir,1:nr)=Linv(ir,1:nr)+wnu(inu)*Lmat(ir,1:nr)*Cs(ilam,1:nr)/(Ce(ilam,ig,1:nr))
						enddo
					enddo
					Linv=-Linv
					do ir=1,nr
						Linv(ir,ir)=1d0+Linv(ir,ir)
					enddo
					Itot(1:nr)=(BBr(ilam,1:nr)*Ca(ilam,ig,1:nr)+Jstar_nu(1:nr,ilam,ig)*Cs(ilam,1:nr))/Ce(ilam,ig,1:nr)
					call DGESV( nr, NRHS, Linv, nr, IWORKomp, Itot, nr, info )
					Si(ilam,ig,1:nr,inu0)=Itot(1:nr)
					do ir=1,nr
						if(.not.Si(ilam,ig,ir,inu0).gt.0d0) then
							Si(ilam,ig,ir,inu0)=BBr(ilam,ir)*Ca(ilam,ig,ir)/Ce(ilam,ig,ir)
						endif
					enddo
				enddo
			else
				do ig=1,ng
					do ir=1,nr
						Si(ilam,ig,ir,inu0)=BBr(ilam,ir)*Ca(ilam,ig,ir)/Ce(ilam,ig,ir)
					enddo
				enddo
			endif
		enddo
!$OMP END DO
		deallocate(tauR)
		deallocate(Itot)
		deallocate(Linv)
		deallocate(Lmat)
		deallocate(IWORKomp)
!$OMP FLUSH
!$OMP END PARALLEL

1		continue
	enddo	
	
	return
	end


	subroutine InvertIjExp(tau,Linv,nr)
	IMPLICIT NONE
	integer ir,nr,jr
	real*8 Linv(nr,nr),tau(nr),Lp(nr,nr),Lm(nr,nr)
	real*8 exptau(nr),exptau2,fact,Q

	do ir=1,nr
		exptau(ir)=exp(-tau(ir))
	enddo
	
	Lp=0d0
	do ir=1,nr-1
		Lp(ir+1,1:nr)=Lp(ir,1:nr)*exptau(ir)
		Q=(1d0-(1d0+tau(ir))*exptau(ir))/tau(ir)
		Lp(ir+1,ir)=Lp(ir+1,ir)+Q
		Q=(tau(ir)-1d0+exptau(ir))/tau(ir)
		Lp(ir+1,ir+1)=Lp(ir+1,ir+1)+Q
	enddo

	Lm=0d0
	do ir=nr,2,-1
		Lm(ir-1,1:nr)=Lm(ir,1:nr)*exptau(ir)
		Q=(1d0-(1d0+tau(ir))*exptau(ir))/tau(ir)
		Lm(ir-1,ir)=Lm(ir-1,ir)+Q
		Q=(tau(ir)-1d0+exptau(ir))/tau(ir)
		Lm(ir-1,ir-1)=Lm(ir-1,ir-1)+Q
	enddo

	Linv=(Lp+Lm)/2d0

	return
	end




	subroutine SolveIjStarExp(tau,I0,Ij,nr)
	IMPLICIT NONE
	integer ir,nr
	real*8 Ij(nr),I0,tau(nr)
	real*8 fact,d,exptau

	fact=1d0
	do ir=nr,1,-1
		exptau=exp(-tau(ir))
		Ij(ir)=(1d0-exptau)*fact*I0
		fact=fact*exptau
	enddo
	
	return
	end
	

	subroutine SolveIjExp(tau,Si,Ij,nr)
	IMPLICIT NONE
	integer ir,nr,jr
	real*8 Ij(nr),Si(nr),tau(nr),Ip(nr),Im(nr)
	real*8 exptau(nr),fact,exptau2,Q

	do ir=1,nr
		exptau(ir)=exp(-tau(ir))
	enddo
	
	Ip=0d0
	do ir=1,nr-1
		Q=(Si(ir)*(1d0-(1d0+tau(ir))*exptau(ir))+Si(ir+1)*(tau(ir)-1d0+exptau(ir)))/tau(ir)
		Ip(ir+1)=Ip(ir)*exptau(ir)+Q
	enddo

	Im=0d0
	do ir=nr,2,-1
		Q=(Si(ir)*(1d0-(1d0+tau(ir))*exptau(ir))+Si(ir-1)*(tau(ir)-1d0+exptau(ir)))/tau(ir)
		Im(ir-1)=Im(ir)*exptau(ir)+Q
	enddo

	Ij=(Ip+Im)/2d0

	return
	end
	

