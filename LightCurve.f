	subroutine LightCurveRetrieval_Fit(rtrace,nrtrace)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nrtrace
	real*8 rtrace(nrtrace)
	integer nc,n,i,j,k,ilam,iobs,im,in,iter,ndata,nl,ilam0
	parameter(nc=3)
	parameter(n=100)
	real*8 Ifunc,x(n),A1(n),A2(n),Ftot(0:nc-1),Intersect,d,Atot(n),F(0:nc-1,n),mu
	external Ifunc
	real*8 Ma,nu,Ma0,Ftot0(0:nc-1),beta,sumt,chi1,chi2,w
	real*8,allocatable :: LC(:,:,:),cLD(:),b(:),Matrix(:,:),ll(:),cLD0(:)

	integer MEQ,MAP,MGT,MODE,MDW
	integer,allocatable :: IP(:)
	real*8 PRGOPT(10),RNORME,RNORML
	real*8,allocatable :: WS(:)
	
	do i=1,n
		x(i)=2d0*acos(1d0-real(i-1)/real(n-1))/pi
	enddo
	
	Ftot0=0d0
	do i=1,n-1
		mu=sqrt(1d0-x(i)**2)
		Atot(i)=pi*(x(i+1)**2-x(i)**2)
		F(0,i)=1d0
		F(1,i)=-(1d0-mu)
		F(2,i)=-(1d0-mu)**2
		Ftot0(0:2)=Ftot0(0:2)+F(0:2,i)*Atot(i)

c		F(0,i)=1d0
c		F(1,i)=-(1d0-sqrt(mu))
c		F(2,i)=-(1d0-mu)
c		F(3,i)=-(1d0-mu**(3d0/2d0))
c		F(4,i)=-(1d0-mu**2)
c		Ftot0(0:4)=Ftot0(0:4)+F(0:4,i)*Atot(i)
	enddo

	nLightCurve=10000

	x=x*Rstar
	Atot=Atot*Rstar**2
	Ftot0=Ftot0*Rstar**2
	call ComputeTc(orbit_omega,orbit_e,Ma0)
	ndata=0
	nl=0
	do iobs=1,nobs
		if(ObsSpec(iobs)%type.eq.'lightcurve') then
			ndata=ndata+ObsSpec(iobs)%nt*ObsSpec(iobs)%nlam
			do ilam=1,ObsSpec(iobs)%nlam
				nl=nl+1
			enddo
		endif
	enddo
	allocate(ll(nl))
	nl=0
	do iobs=1,nobs
		if(ObsSpec(iobs)%type.eq.'lightcurve') then
			do ilam=1,ObsSpec(iobs)%nlam
				nl=nl+1
				ll(nl)=ObsSpec(iobs)%lam(ilam)
			enddo
		endif
	enddo

				MEQ=0
				MAP=ndata
				MGT=2*n-1
				MDW=MEQ+MAP+MGT

			PRGOPT(1)=1
			PRGOPT(2)=0
			PRGOPT(3)=0
			PRGOPT(4)=0
			PRGOPT(5)=0
			PRGOPT(6)=0
			PRGOPT(7)=0

			allocate(WS(2*(MEQ+nc)+max(MAP+MGT,nc)+(MGT+2)*(nc+7)))
			allocate(IP(MGT+2*nc+2))
			IP(1)=2*(MEQ+nc)+max(MAP+MGT,nc)+(MGT+2)*(nc+7)
			IP(2)=MGT+2*nc+2

	allocate(Matrix(MDW,nc+1))
	allocate(b(nc))
	allocate(cLD(nl*nc))
	allocate(cLD0(nl*nc))

	allocate(LC(0:nc-1,ndata,max(nlam,nl)))
	LC=0d0
	ilam0=1
	Matrix=0d0
	do iobs=1,nobs
		if(ObsSpec(iobs)%type.eq.'lightcurve') then
			in=0
			do j=1,ObsSpec(iobs)%nt
				Ma=Ma0+2d0*pi*ObsSpec(iobs)%t(j)/orbit_P
				call TrueAnomaly(Ma,orbit_e,nu)

				d=Dplanet*(1d0-orbit_e**2)*sqrt(1d0-(sin(orbit_omega+nu)*sin(orbit_inc))**2)/(1d0+orbit_e*cos(nu))
				A1(1:n)=0d0
				do k=1,nrtrace-1
					if(d.le.Rstar+rtrace(k)) then
						do i=1,n
							A2(i)=Intersect(x(i),rtrace(k+1),d)
						enddo
						Ftot=0d0
						do i=1,n-1
							Ftot(0:nc-1)=Ftot(0:nc-1)+F(0:nc-1,i)*(A2(i+1)-A1(i+1)-A2(i)+A1(i))
						enddo
						do ilam=ilam0,ilam0+ObsSpec(iobs)%nlam-1
							LC(0:nc-1,j,ilam)=LC(0:nc-1,j,ilam)+obsA_LC(k,ilam)*Ftot(0:nc-1)
						enddo
						A1=A2
					endif
				enddo
				do ilam=ilam0,ilam0+ObsSpec(iobs)%nlam-1
					LC(0:nc-1,j,ilam)=(Ftot0(0:nc-1)-LC(0:nc-1,j,ilam))/Ftot0(0)
				enddo
			enddo

			do ilam=ilam0,ilam0+ObsSpec(iobs)%nlam-1
				do j=1,ObsSpec(iobs)%nt
					in=(j-1)*ObsSpec(iobs)%nlam+ilam-ilam0+1
					Matrix(j,1:nc)=LC(0:nc-1,j,ilam)/ObsSpec(iobs)%dy(in)
					Matrix(j,nc+1)=ObsSpec(iobs)%y(in)/ObsSpec(iobs)%dy(in)
				enddo
				do j=1,n
					Matrix(ObsSpec(iobs)%nt+j,1:nc)=F(0:nc-1,j)
					Matrix(ObsSpec(iobs)%nt+j,nc+1)=0d0
				enddo
				do j=1,n-1
					Matrix(ObsSpec(iobs)%nt+j+n,1:nc)=F(0:nc-1,j)-F(0:nc-1,j+1)
					Matrix(ObsSpec(iobs)%nt+j+n,nc+1)=0d0
				enddo
				MAP=ObsSpec(iobs)%nt
				MGT=2*n-1
				IP(1)=2*(MEQ+nc)+max(MAP+MGT,nc)+(MGT+2)*(nc+7)
				IP(2)=MGT+2*nc+2
				PRGOPT(1)=1
				MODE=0
				call dlsei(Matrix, MDW, MEQ, MAP, MGT, nc, PRGOPT, b, RNORME, RNORML, MODE, WS, IP)
				if(MODE.eq.0) then
	     			cLD((ilam-1)*nc+1:(ilam-1)*nc+nc)=b(1:nc)
	     		else
	     			print*,MODE
	     			cLD((ilam-1)*nc+1)=1d0
	     			cLD((ilam-1)*nc+2:(ilam-1)*nc+nc)=0d0
				endif
			enddo
		endif
		ilam0=ilam0+ObsSpec(iobs)%nlam
	enddo

	cLD0=cLD
	cLD=0d0
	do ilam=1,nl
		sumt=0d0
		do i=1,nl
			w=exp(-((1d0/ll(ilam)-1d0/ll(i))/1000d0)**2)
			cLD((ilam-1)*nc+1:(ilam-1)*nc+nc)=cLD((ilam-1)*nc+1:(ilam-1)*nc+nc)+
     &					w*cLD0((i-1)*nc+1:(i-1)*nc+nc)
			sumt=sumt+w
		enddo
		cLD((ilam-1)*nc+1:(ilam-1)*nc+nc)=cLD((ilam-1)*nc+1:(ilam-1)*nc+nc)/sumt
	enddo

	im=0
	ilam0=1
	do iobs=1,nobs
		if(ObsSpec(iobs)%type.eq.'lightcurve') then
			in=0
			do j=1,ObsSpec(iobs)%nt
				do ilam=ilam0,ilam0+ObsSpec(iobs)%nlam-1
					im=im+1
					in=in+1
					ObsSpec(iobs)%model(in)=sum(LC(0:nc-1,j,ilam)*cLD((ilam-1)*nc+1:(ilam-1)*nc+nc))
					if(.not.ObsSpec(iobs)%model(in).ge.0d0) then
						ObsSpec(iobs)%model(in)=1d0
					endif
				enddo
			enddo
		endif
		ilam0=ilam0+ObsSpec(iobs)%nlam
	enddo

c	open(unit=30,file='LDcoeff.dat',FORM="FORMATTED")
c	do ilam=1,nlam-1
c		write(30,*) lam(ilam),(cLD((ilam-1)*nc+i),i=1,nc)
c	enddo
c	close(unit=30)
c
c	do iobs=1,nobs
c		im=0
c		do ilam=1,ObsSpec(iobs)%nlam
c			open(unit=20,file=trim(outputdir) // "lightcurve" // trim(int2string(iobs,'(i0.4)')) // "_"  // trim(int2string(ilam,'(i0.4)')) // ".dat",FORM="FORMATTED")
c			do j=1,ObsSpec(iobs)%nt
c				im=(j-1)*ObsSpec(iobs)%nlam+ilam
c				write(20,*) ObsSpec(iobs)%t(j),ObsSpec(iobs)%model(im),ObsSpec(iobs)%y(im),(LC(i,j,ilam),i=0,nc-1)
c			enddo
c		enddo
c		close(unit=20)
c	enddo

	deallocate(Matrix)
	deallocate(b)
	deallocate(cLD)
	deallocate(cLD0)
	deallocate(LC)
	deallocate(ll)
	deallocate(WS)
	deallocate(IP)
	
	return
	end
	



	subroutine LightCurveRetrieval(rtrace,nrtrace)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nrtrace
	real*8 rtrace(nrtrace)
	integer nc,n,i,j,k,ilam,iobs,im,in,iter,ndata,nl,ilam0
	parameter(n=100)
	real*8 Ifunc,x(n),A1(n),A2(n),Ftot(0:4),Intersect,d,Atot(n),F(0:4,n),mu
	external Ifunc
	real*8 Ma,nu,Ma0,Ftot0(0:4),beta,sum,chi1,chi2,betamin,betamax
	real*8,allocatable :: LC(:,:,:),cLD(:),b1(:,:),b2(:,:),Matrix(:,:),Matrix0(:,:),ll(:)

	character Trans
	integer Md,Nd,NRHS,LWORK,INFO
	real*8,allocatable :: WORK(:)	
	
	do i=1,n
		x(i)=2d0*acos(1d0-real(i-1)/real(n-1))/pi
	enddo
	
	Ftot0=0d0
	do i=1,n-1
		mu=sqrt(1d0-x(i)**2)
		Atot(i)=pi*(x(i+1)**2-x(i)**2)
		F(0,i)=1d0
		F(1,i)=-(1d0-sqrt(mu))
		F(2,i)=-(1d0-mu)
		F(3,i)=-(1d0-mu**(3d0/2d0))
		F(4,i)=-(1d0-mu**2)
		Ftot0(0:4)=Ftot0(0:4)+F(0:4,i)*Atot(i)
	enddo

	nLightCurve=10000

	x=x*Rstar
	Atot=Atot*Rstar**2
	Ftot0=Ftot0*Rstar**2
	call ComputeTc(orbit_omega,orbit_e,Ma0)
	ndata=0
	nl=0
	do iobs=1,nobs
		if(ObsSpec(iobs)%type.eq.'lightcurve') then
			ndata=ndata+ObsSpec(iobs)%nt*ObsSpec(iobs)%nlam
			do ilam=1,ObsSpec(iobs)%nlam
				nl=nl+1
			enddo
		endif
	enddo
	allocate(ll(nl))
	nl=0
	do iobs=1,nobs
		if(ObsSpec(iobs)%type.eq.'lightcurve') then
			do ilam=1,ObsSpec(iobs)%nlam
				nl=nl+1
				ll(nl)=ObsSpec(iobs)%lam(ilam)
			enddo
		endif
	enddo
	allocate(Matrix(ndata+(nl-1)*5,nl*5))
	allocate(Matrix0(ndata+(nl-1)*5,nl*5))
	allocate(b1(ndata+(nl-1)*5,1))
	allocate(b2(ndata+(nl-1)*5,1))
	allocate(cLD(nl*5))

	allocate(LC(0:4,ndata,max(nlam,nl)))
	LC=0d0
	ilam0=1
	im=0
	b1=0d0
	Matrix=0d0
	do iobs=1,nobs
		if(ObsSpec(iobs)%type.eq.'lightcurve') then
			in=0
			do j=1,ObsSpec(iobs)%nt
				Ma=Ma0+2d0*pi*ObsSpec(iobs)%t(j)/orbit_P
				call TrueAnomaly(Ma,orbit_e,nu)

				d=Dplanet*(1d0-orbit_e**2)*sqrt(1d0-(sin(orbit_omega+nu)*sin(orbit_inc))**2)/(1d0+orbit_e*cos(nu))
				A1(1:n)=0d0
				do k=1,nrtrace-1
					if(d.le.Rstar+rtrace(k)) then
						do i=1,n
							A2(i)=Intersect(x(i),rtrace(k+1),d)
						enddo
						Ftot=0d0
						do i=1,n-1
							Ftot(0:4)=Ftot(0:4)+F(0:4,i)*(A2(i+1)-A1(i+1)-A2(i)+A1(i))
						enddo
						do ilam=ilam0,ilam0+ObsSpec(iobs)%nlam
							LC(0:4,j,ilam)=LC(0:4,j,ilam)+obsA_LC(k,ilam)*Ftot(0:4)
						enddo
						A1=A2
					endif
				enddo
				do ilam=ilam0,ilam0+ObsSpec(iobs)%nlam-1
					LC(0:4,j,ilam)=(Ftot0(0:4)-LC(0:4,j,ilam))/Ftot0(0)
				enddo
				do ilam=ilam0,ilam0+ObsSpec(iobs)%nlam-1
					im=im+1
					in=in+1
					Matrix(im,(ilam-1)*5+1)=LC(0,j,ilam)/ObsSpec(iobs)%dy(in)
					Matrix(im,(ilam-1)*5+2)=LC(1,j,ilam)/ObsSpec(iobs)%dy(in)
					Matrix(im,(ilam-1)*5+3)=LC(2,j,ilam)/ObsSpec(iobs)%dy(in)
					Matrix(im,(ilam-1)*5+4)=LC(3,j,ilam)/ObsSpec(iobs)%dy(in)
					Matrix(im,(ilam-1)*5+5)=LC(4,j,ilam)/ObsSpec(iobs)%dy(in)
					b1(im,1)=ObsSpec(iobs)%y(in)/ObsSpec(iobs)%dy(in)
				enddo
			enddo
		endif
		ilam0=ilam0+ObsSpec(iobs)%nlam
	enddo

	Trans='N'
	Md=ndata+(nl-1)*5
	Nd=nl*5
	NRHS=1
	LWORK = max( 1, min(Md,Nd) + max( min(Md,Nd), NRHS )*2 )
	allocate(WORK(LWORK))

	beta=1d0
	betamin=-1d0
	betamax=-1d0
	do iter=1,3
		k=0
		do ilam=1,nl-1
			k=k+1
			Matrix(ndata+k,(ilam-1)*5+1)=beta/abs(ll(ilam)-ll(ilam+1))
			Matrix(ndata+k,(ilam  )*5+1)=-beta/abs(ll(ilam)-ll(ilam+1))
			k=k+1
			Matrix(ndata+k,(ilam-1)*5+2)=beta/abs(ll(ilam)-ll(ilam+1))
			Matrix(ndata+k,(ilam  )*5+2)=-beta/abs(ll(ilam)-ll(ilam+1))
			k=k+1
			Matrix(ndata+k,(ilam-1)*5+3)=beta/abs(ll(ilam)-ll(ilam+1))
			Matrix(ndata+k,(ilam  )*5+3)=-beta/abs(ll(ilam)-ll(ilam+1))
			k=k+1
			Matrix(ndata+k,(ilam-1)*5+4)=beta/abs(ll(ilam)-ll(ilam+1))
			Matrix(ndata+k,(ilam  )*5+4)=-beta/abs(ll(ilam)-ll(ilam+1))
			k=k+1
			Matrix(ndata+k,(ilam-1)*5+5)=beta/abs(ll(ilam)-ll(ilam+1))
			Matrix(ndata+k,(ilam  )*5+5)=-beta/abs(ll(ilam)-ll(ilam+1))
		enddo
		Matrix0=Matrix
		b2=b1
		call DGELS(Trans,Md,Nd,NRHS,Matrix,Md,b2,Md,WORK,LWORK,INFO)
		cLD(1:Nd)=b2(1:Nd,1)
		Matrix=Matrix0

		chi1=0d0
		do i=1,ndata
			sum=0d0
			do j=1,Nd
				sum=sum+Matrix(i,j)*cLD(j)
			enddo
			chi1=chi1+(sum-b1(i,1))**2
		enddo
		chi2=0d0
		do i=ndata+1,Md
			sum=0d0
			do j=1,Nd
				sum=sum+Matrix(i,j)*cLD(j)
			enddo
			chi2=chi2+(sum-b1(i,1))**2
		enddo
		if(chi2.gt.chi1*0.1) then
			if(betamax.lt.0d0) then
				betamin=beta
				beta=beta*2d0
			else
				beta=(beta+betamax)/2d0
			endif
		else
			if(betamin.lt.0d0) then
				betamax=beta
				beta=beta/2d0
			else
				beta=(beta+betamin)/2d0
			endif
		endif
	enddo

	im=0
	ilam0=1
	do iobs=1,nobs
		if(ObsSpec(iobs)%type.eq.'lightcurve') then
			in=0
			do j=1,ObsSpec(iobs)%nt
				do ilam=ilam0,ilam0+ObsSpec(iobs)%nlam-1
					im=im+1
					in=in+1
					ObsSpec(iobs)%model(in)=LC(0,j,ilam)*cLD((ilam-1)*5+1)+
     &										LC(1,j,ilam)*cLD((ilam-1)*5+2)+
     &										LC(2,j,ilam)*cLD((ilam-1)*5+3)+
     &										LC(3,j,ilam)*cLD((ilam-1)*5+4)+
     &										LC(4,j,ilam)*cLD((ilam-1)*5+5)
				enddo
			enddo
		endif
		ilam0=ilam0+ObsSpec(iobs)%nlam
	enddo

c	open(unit=30,file='LDcoeff.dat',FORM="FORMATTED")
c	do ilam=1,nlam-1
c		write(30,*) lam(ilam),(cLD((ilam-1)*5+i),i=1,5)
c	enddo
c	close(unit=30)
c
c	do iobs=1,nobs
c		im=0
c		do ilam=1,ObsSpec(iobs)%nlam
c			open(unit=20,file=trim(outputdir) // "lightcurve" // trim(int2string(iobs,'(i0.4)')) // "_"  // trim(int2string(ilam,'(i0.4)')) // ".dat",FORM="FORMATTED")
c			do j=1,ObsSpec(iobs)%nt
c				im=(j-1)*ObsSpec(iobs)%nlam+ilam
c				write(20,*) ObsSpec(iobs)%t(j),ObsSpec(iobs)%model(im),ObsSpec(iobs)%y(im),LC(0,j,ilam),LC(1,j,ilam),LC(2,j,ilam),LC(3,j,ilam),LC(4,j,ilam)
c			enddo
c		enddo
c		close(unit=20)
c	enddo

	deallocate(Matrix)
	deallocate(Matrix0)
	deallocate(b1)
	deallocate(b2)
	deallocate(cLD)
	deallocate(LC)
	deallocate(WORK)
	deallocate(ll)
	
	return
	end
	





	subroutine LightCurve(rtrace,nrtrace)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nrtrace
	real*8 cLD(0:4,nlam),rtrace(nrtrace)
	integer nc,n,i,j,k,ilam
	parameter(n=100)
	real*8 Ifunc,x(n),A1(n),A2(n),Ftot(0:4),Intersect,d,Atot(n),F(0:4,n),mu
	external Ifunc
	real*8 Ma,nu,Ma0,Ftot0(0:4)
	real*8,allocatable :: LC(:,:,:)
	
	do i=1,n
		x(i)=2d0*acos(1d0-real(i-1)/real(n-1))/pi
	enddo
	
	Ftot0=0d0
	do i=1,n-1
		mu=sqrt(1d0-x(i)**2)
		Atot(i)=pi*(x(i+1)**2-x(i)**2)
		F(0,i)=1d0
		F(1,i)=-(1d0-sqrt(mu))
		F(2,i)=-(1d0-mu)
		F(3,i)=-(1d0-mu**(3d0/2d0))
		F(4,i)=-(1d0-mu**2)
		Ftot0(0:4)=Ftot0(0:4)+F(0:4,i)*Atot(i)
	enddo

	nLightCurve=10000
	allocate(timeLightCurve(nLightCurve))

	x=x*Rstar
	Atot=Atot*Rstar**2
	Ftot0=Ftot0*Rstar**2
	call ComputeTc(orbit_omega,orbit_e,Ma0)
	allocate(LC(0:4,nLightCurve,nlam))
	LC=0d0
	do j=1,nLightCurve
		timeLightCurve(j)=orbit_P*(real(j-1)/real(nLightCurve-1)-0.5d0)/2d0
		Ma=Ma0+2d0*pi*timeLightCurve(j)/orbit_P
		call TrueAnomaly(Ma,orbit_e,nu)

		d=Dplanet*(1d0-orbit_e**2)*sqrt(1d0-(sin(orbit_omega+nu)*sin(orbit_inc))**2)/(1d0+orbit_e*cos(nu))
		A1(1:n)=0d0
		do k=1,nrtrace-1
			if(d.le.Rstar+rtrace(k)) then
				do i=1,n
					A2(i)=Intersect(x(i),rtrace(k+1),d)
				enddo
				Ftot=0d0
				do i=1,n-1
					Ftot(0:4)=Ftot(0:4)+F(0:4,i)*(A2(i+1)-A1(i+1)-A2(i)+A1(i))
				enddo
				do ilam=1,nlam
					LC(0:4,j,ilam)=LC(0:4,j,ilam)+obsA_LC(k,ilam)*Ftot(0:4)
				enddo
				A1=A2
			endif
		enddo
		do ilam=1,nlam
			LC(0:4,j,ilam)=Ftot0(0:4)-LC(0:4,j,ilam)
		enddo
	enddo

	do ilam=1,nlam
		cLD(0,ilam)=1d0
		cLD(1,ilam)=1d-1*lam(1)/lam(ilam)
		cLD(2,ilam)=1d-1*lam(1)/lam(ilam)
		cLD(3,ilam)=1d-3*lam(1)/lam(ilam)
		cLD(4,ilam)=-1d-2*lam(1)/lam(ilam)
	enddo
	
	open(unit=20,file=trim(outputdir) // "lightcurve.dat",FORM="FORMATTED")
	do j=1,nLightCurve
		write(20,*) timeLightCurve(j),(sum(cLD(0:4,ilam)*LC(0:4,j,ilam))/sum(cLD(0:4,ilam)*Ftot0(0:4)),ilam=1,nlam-1)
	enddo
	close(unit=20)
	
	return
	end
	
	
	real*8 function Ifunc(x,c,nc)
	IMPLICIT NONE
	integer nc
	real*8 x,c(nc),mu
	mu=sqrt(1d0-x**2)
	Ifunc=1d0-c(1)*(1d0-sqrt(mu))-c(2)*(1d0-mu)-c(3)*(1d0-mu**(3d0/2d0))-c(4)*(1d0-mu**2)
	return
	end



	
	real*8 function Intersect(R1,R2,dist)
	IMPLICIT NONE
	real*8 R1,R2,Ar,Br,d,x,y,pi,asinA,asinB,t1,t2,dist
	parameter(pi=3.1415926536)
	d=abs(dist)
	if(R1.gt.R2) then
		Ar=R1
		Br=R2
	else
		Ar=R2
		Br=R1
	endif
	
	x=(Ar**2-Br**2+d**2)/(2d0*d)
	if(x.gt.Ar) then
		if(Ar.gt.d+Br) then
			Intersect=pi*Br**2
		else
			Intersect=0d0
		endif
		return
	else
		y=sqrt(Ar**2-x**2)
	endif
	asinB=asin(y/Br)
	asinA=asin(y/Ar)
	t1=y*sqrt(Ar**2-y**2)
	if(x.lt.d) then
		t2=y*sqrt(Br**2-y**2)
	else
		t2=-y*sqrt(Br**2-y**2)
		asinB=pi-asinB
	endif		
	Intersect=Ar**2*asinA+Br**2*asinB-t1-t2
	
	return
	end



	subroutine TrueAnomaly(Ma,e,nu)
	IMPLICIT NONE
	real*8 Ma,Ea,nu,e,err,sinEa,eps
	eps=1d-4
	
c	Ma=Ea-e*sin(Ea)
	Ea=Ma
	sinEa=sin(Ea)
	err=abs(Ma-Ea+e*sinEa)
	do while(err.gt.eps)
		Ea=Ma+e*sinEa
		sinEa=sin(Ea)
		err=abs(Ma-Ea+e*sinEa)
	enddo
	nu=2d0*atan(sqrt((1d0+e)/(1d0-e))*tan(Ea/2d0))
	if(nu.lt.0d0) nu=nu+2d0*3.1415926536
	
	return
	end	

	subroutine ComputeTc(omega,e,Ma)
	IMPLICIT NONE
	real*8 omega,Ma,Ea,e,nu,pi
	parameter(pi=3.1415926536)

	nu=pi/2d0-omega
	Ea=2d0*atan(sqrt((1d0-e)/(1d0+e))*tan(nu/2d0))
	Ma=Ea-e*sin(Ea)
	
	return
	end
	
