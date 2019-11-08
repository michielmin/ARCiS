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

	subroutine Run3D
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
	real*8,allocatable :: rtrace(:)
	real*8 x,y,z,vx,vy,vz,v
	integer edgeNR,i1,i2,i3,i1next,i2next,i3next,edgenext
	real*8,allocatable :: fluxp(:),tau(:,:),fact(:,:),tautot(:,:),exp_tau(:,:),tau_a(:,:)
	real*8,allocatable :: tauc(:),Afact(:),vv(:,:,:),obsA_omp(:)
	type(Mueller) M
	real*8 g,tot

	allocate(Ca(nlam,ng,nr,n3D),Cs(nlam,nr,n3D),BBr(nlam,nr,n3D),Si(nlam,ng,nr,nnu0,n3D))
	allocate(Ca_mol(nlam,ng,nmol,nr,n3D),Ce(nlam,nr,n3D))
	allocate(R3D(n3D,nr+2))
	allocate(R3D2(n3D,nr+2))

	recomputeopac=.true.

	if(beta3D_1.gt.0d0) then
		b1=beta3D_1	! evening limb, limited to 0-0.25
		b2=beta3D_2	! morning limb, limited to 0-0.25
	else
		b1=betaT
		b2=betaT
	endif
	long0=long_shift*pi/180d0

	call Setup3D(beta,long,latt,nlong,nlatt,long0,b1,b2,Palbedo,betamin,betamax)
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

	do i=1,n3D
c		call tellertje_perc(i,n3D)
		call SetOutputMode(.false.)
		beta3D(i)=betamin+(betamax-betamin)*real(i-1)/real(n3D-1)
		do j=1,n_Par3D
			if(Par3D(j)%logscale) then
				Par3D(j)%x=10d0**(log10(Par3D(j)%xmin)+log10(Par3D(j)%xmax/Par3D(j)%xmin)*beta3D(i))
			else
				Par3D(j)%x=Par3D(j)%xmin+(Par3D(j)%xmax-Par3D(j)%xmin)*beta3D(i)
			endif
		enddo
		call MapPar3D()

		betaT=beta3D(i)

		call InitDens()
		call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam,starfile)
		Fstar=Fstar*pi*Rstar**2
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
C			M%F11=Rayleigh%F11*Csca(ir,ilam)*Ndens(ir)
C			do icloud=1,nclouds
C				if(Cloud(icloud)%standard.eq.'MIX') then
C					M%F11=M%F11+Cloud(icloud)%F(ir,ilam)%F11*Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
C				else
C					do isize=1,Cloud(icloud)%nr
C						M%F11=M%F11+Cloud(icloud)%F(isize,ilam)%F11*
C    &								Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
C					enddo
C				endif
C			enddo
C			g=0d0
C			tot=0d0
C			do j=1,180
C				g=g+M%F11(j)*costheta(j)*sintheta(j)
C				tot=tot+M%F11(j)*sintheta(j)
C			enddo
C			g=g/tot
C			if(.not.g.ge.-1d0) then
C				g=0.999d0
C			endif
C			if(.not.g.le.1d0) then
C				g=0.999d0
C			endif
C			Cs(ilam,ir,i)=Cs(ilam,ir,i)*(1d0-g)
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
		if(emisspec) call ComputeScatter(BBr(1:nlam,1:nr,i),Si(1:nlam,1:ng,1:nr,1:nnu0,i),Ca(1:nlam,1:ng,1:nr,i),Cs(1:nlam,1:nr,i))
		call SetOutputMode(.true.)
	enddo
	Rmax=Rmax*1.1
	R3D(1:n3D,nr+2)=Rmax
	R3D2(1:n3D,nr+2)=Rmax**2

	call output("Raytracing over the planet disk in 3D")


	if(emisspec) then
	ndisk=20
	nsub=0

	nrtrace=(nr-1)*nsub+ndisk
	nptrace=nlatt
	allocate(rtrace(nrtrace))

	k=0
	do i=1,ndisk
		k=k+1
		rtrace(k)=Rmax*real(i-1)/real(ndisk)
	enddo

	allocate(fluxp(nlam),tau(nlam,ng),fact(nlam,ng),tautot(nlam,ng),exp_tau(nlam,ng),tau_a(nlam,ng))
c	npc=8
c	do ipc=1,npc
c	call tellertje_perc(ipc,npc)
c	theta=2d0*pi*real(ipc-1)/real(npc)
	theta=pi
	fluxp=0d0
	do irtrace=1,nrtrace-1
		A=pi*(rtrace(irtrace+1)**2-rtrace(irtrace)**2)/real(nptrace)
		do iptrace=1,nptrace
			phi=2d0*pi*(real(iptrace)-0.5)/real(nptrace)
			rr=0.5d0*(rtrace(irtrace)+rtrace(irtrace+1))
			y=rr*sin(phi)
			z=rr*cos(phi)
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
			tautot=0d0
			fact=1d0
1			continue
			call TravelSph(x,y,z,vx,vy,vz,edgeNR,i1,i2,i3,v,i1next,i2next,i3next,edgenext,nlatt)
			if(i1next.le.0.or.i1next.ge.nr+2) goto 2
			if(i1.le.nr) then
				i=ibeta(i2,i3)
				inu=inu3D(i2,i3)
				tau_a(1:nlam,1:ng)=v*Ca(1:nlam,1:ng,i1,i)
				do ig=1,ng
					tau(1:nlam,ig)=tau_a(1:nlam,ig)+v*Cs(1:nlam,i1,i)
				enddo
				exp_tau=exp(-tau)
				tautot=tautot+tau
				do ig=1,ng
					fluxp(1:nlam)=fluxp(1:nlam)+
     &	wgg(ig)*A*Si(1:nlam,ig,i1,inu,i)*(1d0-exp_tau(1:nlam,ig))*fact(1:nlam,ig)*tau_a(1:nlam,ig)/tau(1:nlam,ig)
				enddo
				fact=fact*exp_tau
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
	fluxp=fluxp*1d23/distance**2
	phase(1,0,1:nlam)=fluxp(1:nlam)
	flux(0,1:nlam)=0d0
c	enddo
	
	deallocate(rtrace)
	deallocate(fluxp,tau,fact,tautot,exp_tau,tau_a)
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
			phi=2d0*pi*(real(iptrace)-0.5)/real(nptrace)
			rr=0.5d0*(rtrace(irtrace)+rtrace(irtrace+1))
			y=rr*sin(phi)
			z=rr*cos(phi)
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
			call TravelSph(x,y,z,vx,vy,vz,edgeNR,i1,i2,i3,v,i1next,i2next,i3next,edgenext,nlatt)
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
	

	subroutine Setup3D(beta,long,latt,nlong,nlatt,long0,b1,b2,albedo,betamin,betamax)
	IMPLICIT NONE
	integer i,j,nlong,nlatt
	real*8 pi
	parameter(pi=3.1415926536)
	real*8 long(nlong),latt(nlatt)	!(Lambda, Phi)
	real*8 beta(nlong,nlatt),la,lo,albedo
	real*8 long0,b1,b2,f,fact,betamin,betamax

	f=max(0d0,((1d0-albedo)-2d0*(b1+b2)))

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
				beta(i,j)=b1+(b2-b1)*(lo+pi/2d0)/pi
			else
				if(lo.lt.0d0) then
					beta(i,j)=b1-(b2-b1)*(lo+pi/2d0)/pi
				else
					beta(i,j)=b2+(b1-b2)*(lo-pi/2d0)/pi
				endif
			endif
			if((lo+long0).gt.0d0) then
				fact=cos((lo+long0)*(pi/2d0)/(long0+pi/2d0))
			else
				fact=cos(-(lo+long0)*(pi/2d0)/(long0-pi/2d0))
			endif
			if(abs(lo).lt.pi/2d0) then
				beta(i,j)=beta(i,j)+f*fact
			endif
			beta(i,j)=abs((b1+b2)/2d0+(beta(i,j)-(b1+b2)/2d0)*cos(la))
			if(beta(i,j).gt.betamax) betamax=beta(i,j)
			if(beta(i,j).lt.betamin) betamin=beta(i,j)
		enddo
	enddo
	
	latt=latt+pi/2d0
	long=long+pi

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

	



	subroutine TravelSph(x,y,z,vx,vy,vz,edgeNR,i1,i2,i3,v,i1next,i2next,i3next,edgenext,n3)
	use Struct3D
	use Constants
	IMPLICIT NONE
	real*8 x,y,z,vx,vy,vz,v
	integer edgeNR,i1,i2,i3,i1next,i2next,i3next,edgenext,n3

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
		edgenext=4
	endif
	if(hitT2.and.vT2.lt.v.and.vT2.gt.0d0) then
		v=vT2
		i1next=i1
		i2next=i2
		i3next=i3+1
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
	real*8 nu(nnu),wnu(nnu),must,tauRs(nr),Ijs(nr),eps
	logical err
	parameter(eps=1d-20)
	real*8,allocatable :: tauR(:),Ij(:),Itot(:),Linv(:,:),Lmat(:,:)
	integer,allocatable :: IWORKomp(:)

	if(.not.scattering) then
		do inu0=1,nnu0
			do ig=1,ng
				Si(1:nlam,ig,1:nr,inu0)=BBr(1:nlam,1:nr)
			enddo
		enddo
		return
	endif

	do ig=1,ng
		Ce(1:nlam,ig,1:nr)=Ca(1:nlam,ig,1:nr)+Cs(1:nlam,1:nr)
	enddo

	tauR_nu=0d0
	do ilam=1,nlam-1
		do ig=1,ng
			do ir=nr,1,-1
				d=abs(P(ir+1)-P(ir))*1d6/grav(ir)
				tau=d*Ce(ilam,ig,ir)/dens(ir)
				if(P(ir).gt.Psimplecloud) then
					tau=tau+1d4
				endif
				if(.not.tau.gt.1d-10) then
					tau=1d-10
				endif
				if(tau.gt.1d4) then
					tau=1d4
				endif
				if(ir.lt.nr) then
					tauR_nu(ir,ilam,ig)=tauR_nu(ir+1,ilam,ig)+tau
				else
					tauR_nu(ir,ilam,ig)=tau
				endif
			enddo
		enddo
	enddo

	call gauleg(0d0,1d0,nu,wnu,nnu)

	allocate(tauR(nr))
	allocate(Ij(nr))
	allocate(Itot(nr))
	allocate(Linv(nr,nr))
	allocate(Lmat(nr,nr))
	allocate(IWORKomp(10*nr*nr))
	do inu0=1,nnu0
		if(inu0.eq.nnu0.or..not.scattstar) then
			Jstar_nu=0d0
			if(inu0.ne.1) then
				Si(1:nlam,1:ng,1:nr,inu0)=Si(1:nlam,1:ng,1:nr,1)
				goto 1
			endif
		else
			must=(real(inu0)-0.5)/real(nnu0-1)
			do ilam=1,nlam-1
				do ig=1,ng
					tauRs(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(must)
					contr=(Fstar(ilam)/(2d0*pi*Dplanet**2))
					call SolveIjStar(tauRs,contr,Ijs,nr)
					Jstar_nu(1:nr,ilam,ig)=Ijs(1:nr)
				enddo
			enddo
		endif
	
		NRHS=1
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

c						tauR(1:nr)=tauR_nu(1:nr,ilam,ig)/abs(nu(inu))
c						do ir=1,nr
c							Itot=0d0
c							Itot(ir)=1d0
c							call SolveIj(tauR,Itot,Ij,nr)
c							Linv(ir,1:nr)=Linv(ir,1:nr)+wnu(inu)*Ij(1:nr)*Cs(ilam,1:nr)/(Ce(ilam,ig,1:nr))
c						enddo
					enddo
					Linv=-Linv
					do ir=1,nr
						Linv(ir,ir)=1d0+Linv(ir,ir)
					enddo
					do ir=1,nr
						contr=Jstar_nu(ir,ilam,ig)
						if(Ce(ilam,ig,ir).eq.0d0) then
							Itot(ir)=BBr(ilam,ir)
						else
							Itot(ir)=(BBr(ilam,ir)*Ca(ilam,ig,ir)+contr*Cs(ilam,ir))/(Ce(ilam,ig,ir))
						endif
					enddo
					call DGESV( nr, NRHS, Linv, nr, IWORKomp, Itot, nr, info )
					Si(ilam,ig,1:nr,inu0)=Itot(1:nr)

					do ir=1,nr
						if(Ca(ilam,ig,ir).eq.0d0.or.(.not.Si(ilam,ig,ir,inu0).gt.0d0)) then
							Si(ilam,ig,ir,inu0)=BBr(ilam,ir)
						else
							Si(ilam,ig,ir,inu0)=Si(ilam,ig,ir,inu0)*Ce(ilam,ig,ir)/Ca(ilam,ig,ir)
						endif
					enddo
				enddo
			else
				do ig=1,ng
					do ir=1,nr
						Si(ilam,ig,ir,inu0)=BBr(ilam,ir)
					enddo
				enddo
			endif
		enddo
1	continue
	enddo	
	deallocate(tauR)
	deallocate(Ij)
	deallocate(Itot)
	deallocate(Linv,Lmat,IWORKomp)
	
	return
	end
	


	subroutine InvertIj(tauR,Linv,nr)
	IMPLICIT NONE
	integer ir,nr,iir
	real*8 tauR(nr),Ij(nr),Si(nr),Linv(nr,nr),Lmat(nr,nr-2)
	real*8 x(nr),y(nr),fact,d
	real*8 MM(nr,3),MMal(nr,1),Ma(nr),Mb(nr),Mc(nr)
	integer indx(nr),info

	Ma=0d0
	Mb=0d0
	Mc=0d0
	do ir=2,nr-1
		fact=1d0/(0.5d0*(tauR(ir+1)+tauR(ir))-0.5d0*(tauR(ir)+tauR(ir-1)))
		Mb(ir)=1d0+fact*(1d0/(tauR(ir+1)-tauR(ir))+1d0/(tauR(ir)-tauR(ir-1)))
		Ma(ir)=-fact*1d0/(tauR(ir)-tauR(ir-1))
		Mc(ir)=-fact*1d0/(tauR(ir+1)-tauR(ir))
	enddo
	Mb(1)=1d0/(tauR(1)-tauR(2))
	Mc(1)=-1d0/(tauR(1)-tauR(2))
	Ma(nr)=1d0/(tauR(nr-1)-tauR(nr))
	Mb(nr)=-1d0-1d0/(tauR(nr-1)-tauR(nr))
	Linv=0d0

	Lmat=0d0
	do iir=1,nr-2
		Lmat(iir+1,iir)=1d0
	enddo
	call dgtsv(nr,nr-2,Ma(2:nr),Mb(1:nr),Mc(1:nr-1),Lmat,nr,info)
	do ir=2,nr-1
		do iir=1,nr
			Linv(ir,iir)=Lmat(iir,ir-1)
		enddo
	enddo

	return

	do iir=2,nr-1
		x=0d0
		x(iir)=1d0
		info=0
		call tridag(Ma,Mb,Mc,x,y,nr,info)
		Ij(1:nr)=y(1:nr)
		do ir=1,nr
			if(Ij(ir).lt.0d0) info=1
		enddo
		if(info.ne.0) then
			do ir=1,nr
				MM(ir,1)=Ma(ir)
				MM(ir,2)=Mb(ir)
				MM(ir,3)=Mc(ir)
			enddo
			call bandec(MM,nr,1,1,nr,3,MMal,1,indx,d)
			x=0d0
			x(iir)=1d0
			call banbks(MM,nr,1,1,nr,3,MMal,1,indx,x)
			Ij(1:nr)=x(1:nr)
		endif
		do ir=1,nr
			if(Ij(ir).lt.0d0) then
				Ij(ir)=0d0
			endif
		enddo
		Linv(iir,1:nr)=Ij(1:nr)
	enddo


	return
	end

