	subroutine MCRad(ilam,flux,docloud)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 z,dz,E,Ca(nr),Cs(nr),Ce(nr),tau,Planck,random,v,flux,dx,dy
	real*8 vR1,vR2,b,rr,R1,R2,tau_v,x,y
	integer iphot,ir,jr,Nphot,ilam,ig,nscat,jrnext
	logical docloud(nclouds),goingup,hitR,onedge,hitR1,hitR2
	type(Mueller) M(nr)
	
	Nphot=50
	flux=0d0
	
	do ir=1,nr
		call GetMatrix(ir,ilam,M(ir),docloud)
	enddo
	do ig=1,ng
		do ir=1,nr
			call Crossections(ir,ilam,ig,Ca(ir),Cs(ir),docloud)
			Ce(ir)=Ca(ir)+Cs(ir)
		enddo

		do ir=1,nr
		tau=0d0
		do jr=ir+1,nr
			tau=tau+Ca(jr)*(R(jr+1)-R(jr))
		enddo
		if(tau.lt.maxtau) then
		E=4d0*pi*(R(ir+1)**3-R(ir)**3)*Planck(T(ir),freq(ilam))*Ca(ir)/3d0
		E=E/real(Nphot)
		do iphot=1,Nphot
			x=0d0
			y=0d0
			z=R(ir)+(R(ir+1)-R(ir))*random(idum)
			call randomdirection(dx,dy,dz)
			jr=ir
			nscat=0
			onedge=.false.
			goingup=(dz.gt.0d0)
1			continue
			tau=-log(random(idum))
2			continue
			rr=x**2+y**2+z**2
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
			v=1d200
			if(hitR1) then
				v=vR1
				goingup=.false.
				jrnext=jr-1
			endif
			if(hitR2.and.vR2.lt.v) then
				v=vR2
				goingup=.true.
				jrnext=jr+1
			endif
			tau_v=v*Ce(jr)
			if(tau_v.lt.tau) then
				x=x+v*dx
				y=y+v*dy
				z=z+v*dz
				tau=tau-tau_v
				jr=jrnext
				if(jr.gt.nr.or.jr.lt.1) goto 3
				onedge=.true.
				goto 2
			endif
			v=tau/Ce(jr)
			x=x+v*dx
			y=y+v*dy
			z=z+v*dz
			call scattering(M(jr),dx,dy,dz)
			nscat=nscat+1
			onedge=.false.
			if(random(idum).lt.(Cs(jr)/Ce(jr))) goto 1
3			continue
			if(jr.gt.nr.and.nscat.gt.0) flux=flux+E/real(ng)
		enddo
		endif
		enddo		
	enddo
	
	return
	end
	
	subroutine Crossections(ir,ilam,ig,Ca,Cs,docloud)
	use GlobalSetup
	IMPLICIT NONE
	integer ir,ilam,ig,icloud,isize
	real*8 Ca,Cs
	logical docloud(nclouds)

	Ca=Cabs(ir,ilam,ig)*Ndens(ir)
	Cs=Csca(ir,ilam)*Ndens(ir)
	do icloud=1,nclouds
		if(docloud(icloud)) then
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
	
	
	subroutine GetMatrix(ir,ilam,M,docloud)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ir,i,ilam,icloud,isize
	type(Mueller) M
	logical docloud(nclouds)
	
	M%F11=Rayleigh%F11*Csca(ir,ilam)*Ndens(ir)
	M%IF11=Rayleigh%IF11*Csca(ir,ilam)*Ndens(ir)
	do icloud=1,nclouds
		if(docloud(icloud)) then
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



	subroutine scattering(M,dx,dy,dz)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,it
	real*8 dz,theta,Fr,Fi,random,dx,dy,x,y,z,rr,u,v,w
	type(Mueller) M
	
	Fr=180d0*random(idum)*M%IF11/pi
	Fi=0d0
	do i=1,180
		Fi=Fi+sintheta(i)*M%F11(i)
		if(Fi.gt.Fr) then
			it=i
			goto 1
		endif
	enddo
1	continue

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
	it=random(idum)*360d0+1
	call rotate(dx,dy,dz,x,y,z,costheta(it),sintheta(it))

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


