	subroutine MCRad(ilam,flux,docloud)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 z,dz,E,Ca(nr),Cs(nr),Ce(nr),tau,Planck,random,d,flux
	integer iphot,ir,jr,Nphot,ilam,ig,nscat
	logical docloud(nclouds)
	
	Nphot=50
	flux=0d0
	
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
		E=2d0*pi*(R(ir+1)**3-R(ir)**3)*Planck(T(ir),freq(ilam))*Ca(ir)
		E=E/real(Nphot)
		do iphot=1,Nphot
			z=R(ir)+(R(ir+1)-R(ir))*random(idum)
			dz=2d0*asin(2d0*random(idum)-1d0)/pi
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
			dz=acos(2d0*random(idum)-1d0)
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
	