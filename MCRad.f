	subroutine MCRad(ilam,fluxg,phase0,docloud0)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iphase
	real*8 z,dz,E,Ca(nr),Cs(nr),Ce(nr),tau,Planck,random,v,fluxg,dx,dy,wphase(nphase)
	real*8 vR1,vR2,b,rr,R1,R2,tau_v,x,y,phase0(nphase),theta,fstop,albedo,ct1,ct2,E0
	real*8 tot,g(nr),Eabs,EJv(nr),Crw(nr)
	integer iphot,ir,jr,Nphot,ilam,ig,nscat,jrnext,NphotStar,NphotPlanet,irdark
	logical docloud0(nclouds),goingup,hitR,onedge,hitR1,hitR2,dorw(nr)
	type(Mueller) M(nr)
	
	NphotPlanet=1000/real(ng)+10
	NphotStar=nphase*5000/real(ng)+10

	EJv=0d0

	fluxg=0d0
	phase0=0d0
	do iphase=1,nphase
		ct1=1d0-2d0*real(iphase-1)/real(nphase)
		ct2=1d0-2d0*real(iphase)/real(nphase)
		wphase(iphase)=2d0/abs(ct1-ct2)
	enddo

	do ir=1,nr
		call GetMatrix(ir,ilam,M(ir),docloud0)
		g(ir)=0d0
		tot=0d0
		do iphase=1,180
			g(ir)=g(ir)+M(ir)%F11(iphase)*costheta(iphase)*sintheta(iphase)
			tot=tot+M(ir)%F11(iphase)*sintheta(iphase)
		enddo
		g(ir)=g(ir)/tot
	enddo
	
	do ig=1,ng
		do ir=1,nr
			call Crossections(ir,ilam,ig,Ca(ir),Cs(ir),docloud0)
			Ce(ir)=Ca(ir)+Cs(ir)
			dorw(ir)=.false.
			Crw(ir)=Ca(ir)+Cs(ir)*(1d0-g(ir))
			if((R(ir+1)-R(ir))*Crw(ir).gt.factRW) dorw(ir)=.true.
		enddo

		irdark=0
		do ir=nr,1,-1
			tau=0d0
			do jr=ir+1,nr
				tau=tau+Ca(jr)*(R(jr+1)-R(jr))
			enddo
			if(tau.gt.maxtau) then
				irdark=ir
				exit
			endif
		enddo

		do ir=irdark+1,nr
			E0=4d0*pi*(R(ir+1)**3-R(ir)**3)*Planck(T(ir),freq(ilam))*Ca(ir)/3d0
			if(E0.gt.0d0) then
				Nphot=NphotPlanet
				E0=E0/real(Nphot)
				do iphot=1,Nphot
					call randomdirection(x,y,z)
					rr=R(ir)+(R(ir+1)-R(ir))*random(idum)
					x=x*rr
					y=y*rr
					z=z*rr
					call randomdirection(dx,dy,dz)
					jr=ir
					onedge=.false.
					goingup=((x*dx+y*dy+z*dz).gt.0d0)
					call travel(x,y,z,dx,dy,dz,jr,onedge,goingup,E,nscat,Ce,Ca,Cs,Crw,g,M,dorw,0.5d0,Eabs,EJv)
					if(jr.gt.nr.and.nscat.gt.0) then
						fluxg=fluxg+E*E0/real(ng)
					endif
				enddo
			endif
		enddo

		if(scattstar) then
		do ir=0,nr
			if(ir.ne.0) then
				E0=Fstar(ilam)*(R(ir+1)**2-R(ir)**2)/Dplanet**2
				Nphot=NphotStar/real(nr)
			else
				E0=Fstar(ilam)*R(1)**2/Dplanet**2
				Nphot=NphotStar
			endif
			E0=E0/real(Nphot)
			do iphot=1,Nphot
				if(ir.eq.0) then
					call randomdisk(x,y,0d0,R(1))
				else
					call randomdisk(x,y,R(ir),R(ir+1))
				endif
				z=sqrt(R(nr+1)**2-x**2-y**2)

				dz=-1d0
				dx=0d0
				dy=0d0
				goingup=.false.
				onedge=.true.
				jr=nr
				call travel(x,y,z,dx,dy,dz,jr,onedge,goingup,E,nscat,Ce,Ca,Cs,Crw,g,M,dorw,0.5d0,Eabs,EJv)
				if(jr.gt.nr.and.nscat.gt.0) then
					iphase=real(nphase)*(1d0-dz)/2d0+1
					if(iphase.lt.1) iphase=1
					if(iphase.gt.nphase) iphase=nphase
					phase0(iphase)=phase0(iphase)+wphase(iphase)*E*E0/real(ng)
				endif
			enddo
		enddo		
		endif
	enddo

	return
	end
	
	
	subroutine travel(x,y,z,dx,dy,dz,jr,onedge,goingup,E,nscat,Ce,Ca,Cs,Crw,g,M,dorw,powstop,Eabs,EJv)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,y,z,dx,dy,dz,E,Ce(nr),Ca(nr),Cs(nr),g(nr),powstop,EJv(nr),Eabs,Crw(nr)
	integer jr,nscat,jrnext,nscat0
	logical onedge,goingup,hitR1,hitR2,hitR,RandomWalk,absorbed,dorw(nr)
	real*8 tau,R1,R2,b,vR1,vR2,rr,v,tau_v,albedo,fstop,random,rho
	type(Mueller) M(nr)

	Eabs=1d0
	E=1d0
	nscat=0
	nscat0=0
	absorbed=.false.

1	continue
	tau=-log(random(idum))
2	continue
	if(absorbed) return
	rr=x**2+y**2+z**2
	rho=sqrt(rr)
	if(dorw(jr)) then
		if(RandomWalk(x,y,z,dx,dy,dz,rho,E,Cs(jr),Ce(jr),Crw(jr),g(jr),jr,nscat,absorbed,0.5d0,Eabs,EJv(jr))) goto 2
	endif
	R1=R(jr)**2
	R2=R(jr+1)**2
	b=2d0*(x*dx+y*dy+z*dz)
	if(onedge) then
		if(goingup) then
			hitR1=.false.
			vR1=1d200
			hitR2=hitR(R2,rr,b,vR2)
		else
			hitR1=hitR(R1,rr,b,vR1)
			hitR2=.true.
			vR2=-b
		endif
	else
		hitR1=hitR(R1,rr,b,vR1)
		hitR2=hitR(R2,rr,b,vR2)
	endif

	v=vR2
	goingup=.true.
	jrnext=jr+1
	if(hitR1.and.vR1.lt.v.and.vR1.gt.0d0) then
		v=vR1
		goingup=.false.
		jrnext=jr-1
	endif
	tau_v=v*Ce(jr)
	albedo=Cs(jr)/(Ca(jr)+Cs(jr))
	if(tau_v.lt.tau) then
		x=x+v*dx
		y=y+v*dy
		z=z+v*dz
		tau=tau-tau_v
		EJv(jr)=EJv(jr)+tau_v*(1d0-albedo)
		jr=jrnext
		if(jr.gt.nr) return
		if(jr.lt.1) then
			jr=1
			call reflectsurface(x,y,z,dx,dy,dz)
			goingup=.true.
			nscat=nscat+1
		endif
		onedge=.true.
		goto 1
	endif
	v=tau/Ce(jr)
	x=x+v*dx
	y=y+v*dy
	z=z+v*dz
	EJv(jr)=EJv(jr)+tau*(1d0-albedo)

	fstop=1d0-albedo**powstop
	if(random(idum).lt.fstop) return

	call scattangle(M(jr),dx,dy,dz)
	nscat=nscat+1
	onedge=.false.
	E=E*albedo/(1d0-fstop)
	goto 1

	return
	end
	
	
	subroutine Crossections(ir,ilam,ig,Ca,Cs,docloud0)
	use GlobalSetup
	IMPLICIT NONE
	integer ir,ilam,ig,icloud,isize
	real*8 Ca,Cs
	logical docloud0(nclouds)

	Ca=Cabs(ir,ilam,ig)*Ndens(ir)
	Cs=Csca(ir,ilam)*Ndens(ir)
	do icloud=1,nclouds
		if(docloud0(icloud)) then
			do isize=1,Cloud(icloud)%nsize
				Ca=Ca+
     &		Cloud(icloud)%Kabs(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
				Cs=Cs+
     &		Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
			enddo
		endif
	enddo

	return
	end
	
	
	subroutine GetMatrix(ir,ilam,M,docloud0)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ir,i,ilam,icloud,isize
	type(Mueller) M
	logical docloud0(nclouds)
	
	M%F11=Rayleigh%F11*Csca(ir,ilam)*Ndens(ir)
	M%IF11=Rayleigh%IF11*Csca(ir,ilam)*Ndens(ir)
	do icloud=1,nclouds
		if(docloud0(icloud)) then
			do isize=1,Cloud(icloud)%nsize
				M%F11=M%F11+Cloud(icloud)%F(isize,ilam)%F11*
     &				Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
				M%IF11=M%IF11+Cloud(icloud)%F(isize,ilam)%IF11*
     &				Cloud(icloud)%Ksca(isize,ilam)*Cloud(icloud)%w(isize)*cloud_dens(ir,icloud)
			enddo
		endif
	enddo

	return
	end



	subroutine scattangle(M,dx,dy,dz)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,it
	real*8 dz,theta,Fr,Fi,random,dx,dy,x,y,z,rr,u,v,w
	type(Mueller) M

	Fr=random(idum)*M%IF11(180)
	it=0
	call hunt(M%IF11,180,Fr,it)
	if(it.lt.1) it=1
	if(it.gt.180) it=180

	x=dx
	y=dy
	z=dz

	u=0d0
	v=-dz
	w=dy
	rr=sqrt(u*u+v*v+w*w)
	u=u/rr
	v=v/rr
	w=w/rr

	call rotate(dx,dy,dz,u,v,w,costheta(it),sintheta(it))
	rr=sqrt(dx*dx+dy*dy+dz*dz)
	dx=dx/rr
	dy=dy/rr
	dz=dz/rr
	it=random(idum)*360d0
	it=it+1
	call rotate(dx,dy,dz,x,y,z,costheta(it),sintheta(it))
	rr=sqrt(dx*dx+dy*dy+dz*dz)
	dx=dx/rr
	dy=dy/rr
	dz=dz/rr

	return
	end


	subroutine reflectsurface(x,y,z,dx,dy,dz)
	use GlobalSetup
	IMPLICIT NONE
	real*8 x,y,z,dx,dy,dz
c for now Lambert surface with albedo 1
	
	call randomdirection(dx,dy,dz)
	if((x*dx+y*dy+z*dz).lt.0d0) then
		dx=-dx
		dy=-dy
		dz=-dz
	endif
	
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotate(x,y,z,u,v,w,cost,sint)
	IMPLICIT NONE
	real*8 x,y,z,u,v,w,yy(3),theta,inp
	real*8 cost,sint,u2,v2,w2
c	cost=cos(theta)
c	sint=sin(theta)
	u2=u*u
	v2=v*v
	w2=w*w

	inp=x*u+y*v+z*w
	yy(1)=u*inp
     & +(x*(v2+w2)-u*(v*y+w*z))*cost
     & +(v*z-w*y)*sint
	yy(2)=v*inp
     & +(y*(u2+w2)-v*(u*x+w*z))*cost
     & +(w*x-u*z)*sint
	yy(3)=w*inp
     & +(z*(u2+v2)-w*(u*x+v*y))*cost
     & +(u*y-v*x)*sint
	x=yy(1)
	y=yy(2)
	z=yy(3)
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotateY(x,y,z,cost,sint)
	IMPLICIT NONE
	real*8 x,y,z,xx,zz,r
	real*8 cost,sint

	xx=x*cost-z*sint
	zz=z*cost+x*sint
	x=xx
	z=zz

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	

	subroutine randomdirection(x,y,z)
	use GlobalSetup
	IMPLICIT NONE
	real*8 x,y,z,s,random
	
1	continue
	x=2d0*random(idum)-1d0
	y=2d0*random(idum)-1d0
	z=2d0*random(idum)-1d0
	s=x**2+y**2+z**2
	if(s.gt.1d0) goto 1
	s=sqrt(s)
	x=x/s
	y=y/s
	z=z/s
	
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	

	subroutine randomdisk(x,y,r1,r2)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,y,s,random,r1,r2,theta,sr1,sr2
	
	sr1=sqrt(r1)
	sr2=sqrt(r2)
	s=sr1+(sr2-sr1)*random(idum)
	theta=random(idum)*2d0*pi
	s=s*s
	
	x=s*cos(theta)
	y=s*sin(theta)

	return
	end


	logical function hitR(Rad,r,b,v)
	IMPLICIT NONE
	real*8 Rad,r,b,cc,discr,vr1,vr2,v,q
	
	hitR=.false.
	v=1d200

	cc=r-Rad
	discr=(b**2-4d0*cc)
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(b.gt.0d0) then
			q=-0.5d0*(b+discr)
		else
			q=-0.5d0*(b-discr)
		endif
		vr1=q
		vr2=cc/q
		if(vr1.gt.0d0) then
			v=vr1
			hitR=.true.
		endif
		if(vr2.gt.0d0.and.vr2.lt.v) then
			v=vr2
			hitR=.true.
		endif
	endif
	return
	end


	module RandomWalkModule
	IMPLICIT NONE
	integer NY
	parameter(NY=1000)
	real*8 phi(NY),yy(NY)
	end module RandomWalkModule

	
	logical function RandomWalk(x,y,z,dx,dy,dz,rr,E,Cs,Ce,Crw,g,jr,nscat,absorbed,powstop,Eabs,EJv)
	use GlobalSetup
	use Constants
	use RandomWalkModule
	IMPLICIT NONE
	real*8 dmin,v,ry,random,lr,x,y,z,dx,dy,dz,E,Crw,Kappa,g,Cs,Ce
	real*8 rr,albedo,d1,d2,fstop,powstop,EJv,Eabs
	integer i,jr,nscat,iy
	logical absorbed

	RandomWalk=.false.

	lr=1d0/Crw

	d1=abs(rr-R(jr))
	d2=abs(rr-R(jr+1))
	dmin=d1
	if(d2.lt.dmin) dmin=d2

	if(dmin.le.factRW*lr) return

	RandomWalk=.true.

	ry=random(idum)
	iy=1
	call hunt(phi,NY,ry,iy)

	v=-3d0*log(yy(iy))*dmin**2/(lr**2*pi**2)

	call randomdirection(dx,dy,dz)
	x=x+dmin*dx
	y=y+dmin*dy
	z=z+dmin*dz

	albedo=Cs/Ce
	fstop=1d0-albedo**powstop
	fstop=fstop**v
	absorbed=.false.

	E=E*(albedo/(1d0-fstop))**v

	EJv=EJv+v*(1d0-albedo)
	Eabs=Eabs+v*(1d0-albedo)

	nscat=nscat+v

	if(random(idum).gt.fstop) absorbed=.true.

	return
	end
	



	subroutine InitRandomWalk()
	use RandomWalkModule
	IMPLICIT NONE
	integer i,n,nmax
	nmax=1000

	yy(1)=1d0
	phi(1)=1d0
	do i=2,NY
		yy(i)=real(NY-i+1)/real(NY-1)*0.75d0
		phi(i)=0d0
		do n=1,nmax
			phi(i)=phi(i)+(-1d0)**(n+1)*yy(i)**(n**2)
		enddo
		phi(i)=phi(i)*2d0
	enddo

	return
	end
	


