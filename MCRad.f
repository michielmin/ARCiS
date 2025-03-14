	subroutine MCRad(ilam,fluxg,phase0,docloud0)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iphase
	real*8 z,dz,E,Ca(nr),Cs(nr),Ce(nr),tau,Planck,random,v,fluxg,dx,dy,wphase(nphase)
	real*8 vR1,vR2,b,rr,R1,R2,tau_v,x,y,phase0(nphase),theta,fstop,albedo,ct1,ct2,E0
	real*8 tot,g(nr),Eabs,EJv(nr),Crw(nr),F11(180),Gsca
	real*8 Eemit(nr),Femit(nr),xs,ys,zs,inp,Pb(nr+1)
	integer iphot,ir,jr,Nphot,ilam,ig,nscat,jrnext,NphotStar,NphotPlanet,irdark,j
	logical docloud0(nclouds),goingup,hitR,onedge,hitR1,hitR2,dorw(nr)
	real*8 starttime,stoptime
	logical hitplanet
	
	NphotPlanet=Nphot0/real(ng)+5
	NphotStar=nphase*Nphot0/real(ng*45)+5

	EJv=0d0

	fluxg=0d0
	phase0=0d0
	theta=0d0
	ct1=cos(theta)
	do iphase=1,nphase
c		ct1=1d0-2d0*real(iphase-1)/real(nphase)
c		ct2=1d0-2d0*real(iphase)/real(nphase)
		theta=real(iphase)*pi/real(nphase)
		ct2=cos(theta)
		wphase(iphase)=2d0/abs(ct1-ct2)
		if(IsNaN(wphase(iphase))) print*,wphase(iphase)
		ct1=ct2
	enddo

	do ir=1,nr
		g(ir)=0d0
	enddo
	
	Pb(1)=P(1)
	do ir=2,nr
		Pb(ir)=sqrt(P(ir-1)*P(ir))
	enddo
	Pb(nr+1)=P(nr)
	do ig=1,ng
		do ir=1,nr
			call Crossections(ir,ilam,ig,Ca(ir),Cs(ir),docloud0,0,F11,g(ir),.true.)
			Ce(ir)=Ca(ir)+Cs(ir)
			dorw(ir)=.false.
			Crw(ir)=Ca(ir)+Cs(ir)*(1d0-g(ir))
			if((R(ir+1)-R(ir))*Crw(ir).gt.factRW) dorw(ir)=.true.
		enddo


		tau=0d0
		Eemit=0d0
		Femit=0d0
		do ir=nr,1,-1
			Eemit(ir)=(4d0*pi/3d0)*(R(ir+1)**3-R(ir)**3)*Planck(T(ir),freq(ilam))*Ca(ir)
			Femit(ir)=Eemit(ir)*exp(-tau)
			tau=tau+Ca(ir)*(R(ir+1)-R(ir))
		enddo
		E0=sum(Femit(1:nr))
		Femit=Femit/E0
		
		Nphot=NphotPlanet
		do iphot=1,Nphot
			E0=random(idum)
			do ir=nr,1,-1
				E0=E0-Femit(ir)
				if(E0.lt.0d0) exit
			enddo
			if(ir.lt.1) ir=1
			if(ir.gt.nr) ir=nr
			E0=Eemit(ir)/(real(Nphot)*Femit(ir))
			call randomdirection(x,y,z)
			rr=(R(ir)**3+(R(ir+1)**3-R(ir)**3)*random(idum))**(1d0/3d0)
			x=x*rr
			y=y*rr
			z=z*rr
			call randomdirection(dx,dy,dz)
			jr=ir
			onedge=.false.
			goingup=((x*dx+y*dy+z*dz).gt.0d0)
			call travel(x,y,z,dx,dy,dz,jr,onedge,goingup,E,nscat,Ce,Ca,Cs,Crw,g,dorw,0.5d0,Eabs,EJv)
			if(jr.gt.nr.and.nscat.gt.0) then
				fluxg=fluxg+E*E0*wgg(ig)
			endif
		enddo


c		Nphot=NphotPlanet
c		do iphot=1,Nphot
c			ir=1
c			E0=4d0*pi*R(ir)**2*Planck(T(ir),freq(ilam))/real(Nphot)
c			call randomdirection(x,y,z)
c			rr=R(ir)**0.99*R(ir+1)**0.01
c			x=x*rr
c			y=y*rr
c			z=z*rr
c			call randomdirection(dx,dy,dz)
c			jr=ir
c			onedge=.false.
c			goingup=((x*dx+y*dy+z*dz).gt.0d0)
c			call travel(x,y,z,dx,dy,dz,jr,onedge,goingup,E,nscat,Ce,Ca,Cs,Crw,g,M,dorw,0.5d0,Eabs,EJv)
c			if(jr.gt.nr) then
c				fluxg=fluxg+E*E0*wgg(ig)
c			endif
c		enddo


		if(scattstar) then
c		do ir=0,nr
c			if(ir.ne.0) then
c				E0=Fstar(ilam)*(R(ir+1)**2-R(ir)**2)/Dplanet**2
c				Nphot=NphotStar/real(nr)
c			else
c				E0=Fstar(ilam)*R(1)**2/Dplanet**2
c				Nphot=NphotStar
c			endif
c			E0=E0/real(Nphot)
c			do iphot=1,Nphot
c				if(ir.eq.0) then
c					call randomdisk(x,y,0d0,R(1))
c				else
c					call randomdisk(x,y,R(ir),R(ir+1))
c				endif
c				z=sqrt(R(nr+1)**2-x**2-y**2)
c
c				dz=-1d0
c				dx=0d0
c				dy=0d0
c				goingup=.false.
c				onedge=.true.
c				jr=nr
c				call travel(x,y,z,dx,dy,dz,jr,onedge,goingup,E,nscat,Ce,Ca,Cs,Crw,g,M,dorw,0.5d0,Eabs,EJv)
c				if(jr.gt.nr.and.nscat.gt.0) then
c					iphase=real(nphase)*(1d0-dz)/2d0+1
c					if(dz.gt.1d0) dz=1d0
c					if(dz.lt.-1d0) dz=-1d0
c					iphase=real(nphase)*(acos(dz)/pi)+1
c					if(iphase.lt.1) iphase=1
c					if(iphase.gt.nphase) iphase=nphase
c					phase0(iphase)=phase0(iphase)+wphase(iphase)*E*E0/real(ng)
c				endif
c			enddo
c		enddo		


		do ir=0,nr
			if(ir.ne.0) then
				E0=Fstar(ilam)*(R(ir+1)**2-R(ir)**2)/(Dplanet**2)
				Nphot=NphotStar/real(nr)
			else
				E0=Fstar(ilam)*R(1)**2/(Dplanet**2)
				Nphot=NphotStar
			endif
			E0=E0/real(Nphot)
			do iphot=1,Nphot
1				if(ir.eq.0) then
					call randomdisk(x,y,0d0,R(1))
				else
					call randomdisk(x,y,R(ir),R(ir+1))
				endif
				z=sqrt(R(nr+1)**2-x**2-y**2)

				call randomdirection(xs,ys,zs)
				xs=xs*Rstar
				ys=ys*Rstar
				zs=zs*Rstar+Dplanet
			
				dx=x-xs
				dy=y-ys
				dz=z-zs

				inp=x*dx+y*dy+z*dz
				if(inp.gt.0d0) goto 1
				inp=xs*dx+ys*dy+(zs-Dplanet)*dz
				if(inp.lt.0d0) goto 1
				
				tot=sqrt(dx**2+dy**2+dz**2)
				dx=dx/tot
				dy=dy/tot
				dz=dz/tot

				goingup=.false.
				onedge=.true.
				jr=nr
				call travel(x,y,z,dx,dy,dz,jr,onedge,goingup,E,nscat,Ce,Ca,Cs,Crw,g,dorw,0.5d0,Eabs,EJv)
				if(jr.gt.nr.and.nscat.gt.0) then
					iphase=real(nphase)*(1d0-dz)/2d0+1
					if(dz.gt.1d0) dz=1d0
					if(dz.lt.-1d0) dz=-1d0
					iphase=real(nphase)*(acos(dz)/pi)+1
					if(iphase.lt.1) iphase=1
					if(iphase.gt.nphase) iphase=nphase
					phase0(iphase)=phase0(iphase)+wphase(iphase)*E*E0*wgg(ig)
				endif
			enddo
		enddo
		endif
	enddo

	return
	end

	logical function hitplanet(x,y,z,dx,dy,dz)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,y,z,dx,dy,dz,R1,R2,aa,b,cc
	real*8 v,vr1,vr2,qq,discr,xt,yt,zt,vxt,vyt,vzt

	xt=x
	yt=y
	zt=z+Dplanet
	vxt=dx
	vyt=dy
	vzt=dz

	R1=R(nr+1)**2
	R2=xt**2+yt**2+zt**2

	aa=vxt**2+vyt**2+vzt**2
	b=2d0*(xt*vxt+yt*vyt+zt*vzt)
	cc=R2-R1
	discr=(b**2-4d0*cc*aa)
	hitplanet=.false.
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(b.gt.0d0) then
			qq=-0.5d0*(b+discr)
		else
			qq=-0.5d0*(b-discr)
		endif
		vr1=qq/aa
		vr2=cc/qq
		v=vr1
		hitplanet=.true.
		if(vr2.lt.v) then
			v=vr2
			hitplanet=.true.
		endif
		x=x+v*dx
		y=y+v*dy
		z=z+v*dz
	endif

	return
	end
	

	
	
	subroutine travel(x,y,z,dx,dy,dz,jr,onedge,goingup,E,nscat,Ce,Ca,Cs,Crw,g,dorw,powstop,Eabs,EJv)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,y,z,dx,dy,dz,E,Ce(nr),Ca(nr),Cs(nr),g(nr),powstop,EJv(nr),Eabs,Crw(nr)
	integer jr,nscat,jrnext,nscat0
	logical onedge,goingup,hitR1,hitR2,hitR,RandomWalk,absorbed,dorw(nr)
	real*8 tau,R1,R2,b,vR1,vR2,rr,v,tau_v,albedo,fstop,random,rho

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
		if(RandomWalk(x,y,z,dx,dy,dz,rho,E,Cs(jr),Ce(jr),Crw(jr),g(jr),jr,nscat,absorbed,powstop,Eabs,EJv(jr))) goto 1
	endif
	R1=R(jr)**2
	R2=R(jr+1)**2
	b=2d0*(x*dx+y*dy+z*dz)

	goingup=((x*dx+y*dy+z*dz).gt.0d0)
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
		if(goingup) then
			hitR1=.false.
			vR1=1d200
			hitR2=hitR(R2,rr,b,vR2)
		else
			hitR1=hitR(R1,rr,b,vR1)
			hitR2=hitR(R2,rr,b,vR2)
		endif
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
	albedo=Cs(jr)/Ce(jr)
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
		goto 2
	endif
	v=tau/Ce(jr)
	x=x+v*dx
	y=y+v*dy
	z=z+v*dz
	EJv(jr)=EJv(jr)+tau*(1d0-albedo)

	fstop=max(1d0-albedo**powstop,1d-6)
	if(random(idum).lt.fstop.and.powstop.gt.0d0) return

	if(powstop.gt.0d0) then
		E=E*albedo/(1d0-fstop)
	else
		return
	endif
	call randomdirection(dx,dy,dz)
	nscat=nscat+1
	onedge=.false.
	goto 1

	return
	end
	
	
	subroutine Crossections(ir,ilam,ig,Ca,Cs,docloud0,ivel,F11,G,doF11)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ir,ilam,ig,icloud,isize,ivel,i
	real*8 Ca,Cs,F11(180),theta,CsAdd,G
	logical docloud0(nclouds),doF11

	Ca=Cabs(ir,ilam,ig,ivel)*Ndens(ir)
	Cs=Csca(ir,ilam)*Ndens(ir)
	if(.not.Ca.gt.0d0) Ca=0d0
	if(.not.Cs.gt.0d0) Cs=0d0
	if(doF11) then
		do i=1,180
			theta=real(i)*pi/180d0
			F11(i)=Cs*3d0*(1d0+cos(theta)**2)/4d0
		enddo
	endif

	G=0d0
	do icloud=1,nclouds
		if(docloud0(icloud)) then
			Ca=Ca+Cloud(icloud)%Kabs(ir,ilam)*cloud_dens(ir,icloud)
			CsAdd=Cloud(icloud)%Ksca(ir,ilam)*cloud_dens(ir,icloud)
			Cs=Cs+CsAdd
			if(doF11) then
				G=G+CsAdd*Cloud(icloud)%G(ir,ilam)
				do i=1,180
					F11(i)=F11(i)+CsAdd*Cloud(icloud)%F11(ir,ilam,i)
				enddo
			endif
		endif
	enddo
	if(.not.Ca.gt.0d0) Ca=0d0
	if(.not.Cs.gt.0d0) Cs=0d0
	if(doF11) then
		if(Cs.gt.0d0) then
			G=G/Cs
			F11=F11/Cs
		else
			G=0d0
			F11=1d0
		endif
	else
		G=0d0
		F11=1d0
	endif

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
	real*8 cost,sint,u2,v2,w2,ux,vy,wz
c	cost=cos(theta)
c	sint=sin(theta)
	u2=u*u
	v2=v*v
	w2=w*w

	ux=u*x
	vy=v*y
	wz=w*z

	inp=ux+vy+wz
	yy(1)=u*inp
     & +(x*(v2+w2)-u*(vy+wz))*cost
     & +(v*z-w*y)*sint
	yy(2)=v*inp
     & +(y*(u2+w2)-v*(ux+wz))*cost
     & +(w*x-u*z)*sint
	yy(3)=w*inp
     & +(z*(u2+v2)-w*(ux+vy))*cost
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
	
	sr1=r1**2
	sr2=r2**2
	s=sr1+(sr2-sr1)*random(idum)
	theta=random(idum)*2d0*pi
	s=sqrt(s)
	
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

	albedo=Cs/Ce
	fstop=max(1d0-albedo**powstop,1d-4)
	if(powstop.le.0d0) fstop=0d0
	if(powstop.gt.0d0) E=E*(albedo/(1d0-fstop))**v

	fstop=1d0-(1d0-fstop)**v
	absorbed=.false.

	if(random(idum).lt.fstop) then
		absorbed=.true.
		return
	endif

	call randomdirection(dx,dy,dz)
	x=x+dmin*dx
	y=y+dmin*dy
	z=z+dmin*dz

	EJv=EJv+v*(1d0-albedo)
	Eabs=Eabs+v*(1d0-albedo)

	nscat=nscat+v


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
	


