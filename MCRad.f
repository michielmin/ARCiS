	subroutine MCRad(ilam,flux,docloud)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 z,dz,E,Ca(nr),Cs(nr),Ce(nr),tau,Planck,random,d,flux,dx,dy
	integer iphot,ir,jr,Nphot,ilam,ig,nscat
	logical docloud(nclouds)
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
			z=R(ir)+(R(ir+1)-R(ir))*random(idum)
			call randomdirection(dx,dy,dz)
			jr=ir
			nscat=0
1			continue
			tau=-log(random(idum))
2			d=dz*tau/(Ce(jr))
			if(d.lt.0d0) then
				if((z+d).lt.R(jr)) then
					tau=tau-abs(z-R(jr))*Ce(jr)/dz
					z=R(jr)
					jr=jr-1
					if(jr.eq.0) goto 3
					goto 2
				endif
			else
				if((z+d).gt.R(jr+1)) then
					tau=tau-abs(z-R(jr+1))*Ce(jr)/dz
					z=R(jr+1)
					jr=jr+1
					if(jr.gt.nr) goto 3
					goto 2
				endif
			endif
			z=z+d
			call scattering(M(jr),dz)
			nscat=nscat+1
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



	subroutine scattering(M,dz)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,it
	real*8 dz,theta,Fr,Fi,random,dx,dy,x,y,z
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

	dx=sqrt(1d0-dz**2)
	dy=0d0
	x=dx
	y=dy
	z=dz

	call rotateY(dx,dy,dz,costheta(it),sintheta(it))
	theta=2d0*pi*random(idum)
	call rotate(dx,dy,dz,x,y,z,theta)

	return
	end



c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotate(x,y,z,u,v,w,theta)
	IMPLICIT NONE
	real*8 x,y,z,u,v,w,yy(3),theta,inp
	real*8 cost,sint,u2,v2,w2
	cost=cos(theta)
	sint=sin(theta)
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

