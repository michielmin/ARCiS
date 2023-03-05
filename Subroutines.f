
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END



	subroutine Smooth(y,dy,z,n,idum)
	IMPLICIT NONE
	integer m,n,MAXM,MAXITER
	parameter(MAXM=1000,MAXITER=100)
	real*8 y(n),dy(n),x(n),a(MAXM,n),b(MAXM),min,z(n),maxsigma
	real*8 rnorm,w(MAXM*3),beta(MAXITER),error(MAXITER),eps
	integer info,i,j,niter,ii,idum
	real*8 gasdev
	m=2*n-1

	eps=1d-2
	beta(1)=0.1d0

	maxsigma=1d0

	niter=1
1	continue
	do i=1,n
		do j=1,n
			a(i,j)=0d0
			a(i+n,j)=0d0
		enddo
		a(i,i)=1d0/dy(i)
		b(i)=(y(i)+gasdev(idum)*dy(i))/dy(i)
		if(i.lt.n.and.i.gt.1) then
			a(i+n,i-1)=-beta(niter)
			a(i+n,i)=2d0*beta(niter)
			a(i+n,i+1)=-beta(niter)
			b(i+n)=0d0
		endif
	enddo

	call dgels('N',m,n,1,a,MAXM,b,MAXM,w,MAXM*3,info)

	error(niter)=0d0
	do i=1,n
		error(niter)=error(niter)+((y(i)-b(i))/dy(i))**2
	enddo
	error(niter)=error(niter)/real(n-1)
	
	if(abs(error(niter)-maxsigma).gt.eps) then
		niter=niter+1
		if(niter.gt.MAXITER) goto 2
		beta(niter)=beta(niter-1)/sqrt(error(niter-1))
		goto 1
	endif

	goto 3

2	min=0d0
	ii=MAXITER
	do i=1,MAXITER
		if(error(i).gt.min.and.error(i).lt.1d0) then
			min=error(i)
			ii=i
		endif
	enddo
	
	do i=1,n
		do j=1,n
			a(i,j)=0d0
			a(i+n,j)=0d0
		enddo
		a(i,i)=1d0/dy(i)
		b(i)=y(i)/dy(i)
		if(i.lt.n.and.i.gt.1) then
			a(i+n,i-1)=-beta(ii)
			a(i+n,i)=2d0*beta(ii)
			a(i+n,i+1)=-beta(ii)
			b(i+n)=0d0
		endif
	enddo
	call dgels('N',m,n,1,a,MAXM,b,MAXM,w,MAXM*3,info)
	
3	continue
	do i=1,n
		z(i)=b(i)
	enddo

	return
	end



	real*8 function random(idum)
	IMPLICIT NONE
	real*8 ran1
	integer idum

	random=ran1(idum)

	return
	end
	
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
!$OMP THREADPRIVATE(iv,iy)
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END



      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END




      FUNCTION gasdev(idum)
      INTEGER idum
      REAL*8 gasdev
CU    USES ran1
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,random
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0.or.idum.lt.0) then
1       v1=2.*random(idum)-1.
        v2=2.*random(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END



c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine regrid(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 grid(n),y(n),x0,y0,x1,y1
	character*500 input
	logical truefalse
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,FORM="FORMATTED",ACCESS="STREAM")
	i=1
1	read(20,*,end=102,err=1) x0,y0
103	if(x0.ge.grid(i)) then
		y(i)=y0
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) x1,y1
101	if(grid(i).le.x1.and.grid(i).ge.x0) then
		y(i)=y1+(grid(i)-x1)*(y0-y1)/(x0-x1)
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	goto 100
102	continue
	do j=i,n
		y(j)=y(i-1)*grid(i-1)/grid(j)
	enddo
	close(unit=20)
	return
	end



c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine regridN(input,grid_in,y_in,n,i1,i2,nn,nskip,invert,extrapol)
	IMPLICIT NONE
	integer i,j,n,i1,i2,nn,nskip
	real*8 grid(n),y(n,nn),x0,y0(nn),x1,y1(nn),dummy(max(i1,i2+nn))
	real*8 grid_in(n),y_in(n,nn)
	character*500 input
	logical truefalse,invert,extrapol
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	if(invert) then
		do i=1,n
			grid(i)=grid_in(n-i+1)
			y(i,1:nn)=y_in(n-i+1,1:nn)
		enddo
	else
		grid=grid_in
		y=y_in
	endif

	y0=1d-60
	y1=1d-60
	y=1d-60
	open(unit=20,file=input,FORM="FORMATTED",ACCESS="STREAM")
	do j=1,nskip
		read(20,*)
	enddo
	i=1
1	read(20,*,end=102,err=1) dummy(1:max(i1,i2+nn-1))
	x0=dummy(i1)
	y0(1:nn)=dummy(i2:i2+nn-1)
c	do j=1,nn
c		if(y0(j).le.0d0) y0(j)=1d-50
c	enddo
103	if((invert.and.x0.ge.grid(i)).or.(.not.invert.and.x0.le.grid(i))) then
		if(extrapol) y(i,1:nn)=y0(1:nn)
		i=i+1
		goto 103
	endif
100	read(20,*,end=102,err=100) dummy(1:max(i1,i2+nn-1))
	x1=dummy(i1)
	y1(1:nn)=dummy(i2:i2+nn-1)
c	do j=1,nn
c		if(y1(j).le.0d0) y1(j)=y0(j)*1d-50
c	enddo
101	if((invert.and.(grid(i).le.x1.and.grid(i).ge.x0)).or.
     &		(.not.invert.and.(grid(i).ge.x1.and.grid(i).le.x0))) then
c		y(i,1:nn)=10d0**(log10(y1(1:nn))+(log10(grid(i)/x1))*(log10(y0(1:nn)/y1(1:nn)))/(log10(x0/x1)))
		y(i,1:nn)=y1(1:nn)+(grid(i)-x1)*(y0(1:nn)-y1(1:nn))/(x0-x1)
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	goto 100
102	continue
	if(extrapol) then
		do j=i,n
			y(j,1:nn)=y0(1:nn)
		enddo
	endif
	close(unit=20)
c	do i=1,n
c		do j=1,nn
c			if(y(i,j).le.0d0) y(i,j)=1d-60
c		enddo
c	enddo

	if(invert) then
		do i=1,n
			y_in(i,1:nn)=y(n-i+1,1:nn)
		enddo
	else
		y_in=y
	endif

	return
	end



c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine regridlog(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 grid(n),y(n),x0,y0,x1,y1
	real*8 lx0,ly0,lx1,ly1,lx,ly
	character*500 input
	logical truefalse
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,FORM="FORMATTED",ACCESS="STREAM")
	i=1
1	read(20,*,end=102,err=1) x0,y0
103	if(x0.ge.grid(i)) then
		y(i)=y0
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) x1,y1
	if(y1.le.0d0) y1=y0*1d-50
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		y(i)=10d0**(log10(y1)+(log10(grid(i)/x1))*(log10(y0/y1))/(log10(x0/x1)))
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	goto 100
102	continue
	do j=i,n
		y(j)=y0*x0/grid(j)
	enddo
	close(unit=20)
	do i=1,n
		if(y(i).le.0d0) y(i)=y0*1d-60
	enddo
	return
	end


c-----------------------------------------------------------------------
c The new readstar subroutine uses a boxcar filtering to read in 
c high resolution spectra.
c-----------------------------------------------------------------------

	subroutine regridstar(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n,nls
	real*8 grid(n),y(n),x0,y0,xedge(n+1)
	real*8 grid2(n),y2(n),tot(n)
	real*8,allocatable :: ls(:),Fs(:),dls(:)
	character*500 input
	logical truefalse,done(n)

	do i=1,n-1
		xedge(i+1)=sqrt(grid(i)*grid(i+1))
	enddo
	xedge(1)=grid(1)**2/xedge(2)
	xedge(n)=grid(n)**2/xedge(n-1)

	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,FORM="FORMATTED",ACCESS="STREAM")
	j=1
1	read(20,*,end=2,err=1) x0,y0
	j=j+1
	goto 1
2	continue
	nls=j-1
	close(unit=20)
	allocate(ls(nls))
	allocate(dls(nls))
	allocate(Fs(nls))

	open(unit=20,file=input,FORM="FORMATTED",ACCESS="STREAM")
	do i=1,nls
3		read(20,*,end=4,err=3) ls(i),Fs(i)
	enddo
4	close(unit=20)
	do i=2,nls-1
		dls(i)=(1d0/ls(i-1)-1d0/ls(i+1))/2d0
	enddo
	dls(1)=dls(2)
	dls(nls)=dls(nls-1)
	
	j=1
	tot=0d0
	y=0d0
	done=.false.
	do i=1,nls
5		continue
		if(ls(i).lt.xedge(j)) goto 6
		if(ls(i).gt.xedge(j+1)) then
			j=j+1
			if(j.gt.n) goto 7
			goto 5
		endif
		y(j)=y(j)+Fs(i)*dls(i)
		tot(j)=tot(j)+dls(i)
		done(j)=.true.
6		continue
	enddo
7	continue
	j=0
	do i=1,n
		if(done(i)) then
			y(i)=y(i)/tot(i)
		else
			j=j+1
			grid2(j)=grid(i)
		endif
	enddo
	if(j.ne.0) then
		call readstar_interpol(input,grid2,y2,j)
		j=0
		do i=1,n
			if(.not.done(i)) then
				j=j+1
				y(i)=y2(j)
			endif
		enddo
	endif
	
	return
	end
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine readstar_interpol(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 grid(n),y(n),x0,y0,x1,y1
	character*500 input
	logical truefalse
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,FORM="FORMATTED",ACCESS="STREAM")
	i=1
1	read(20,*,end=102,err=1) x0,y0
103	if(x0.ge.grid(i)) then
		y(i)=y0
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) x1,y1
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		if(y1.gt.1d-60.and.y0.gt.1d-60) then
			y(i)=10d0**(log10(y1)+(log10(grid(i))-log10(x1))*(log10(y0)-log10(y1))/(log10(x0)-log10(x1)))
		else
			y(i)=y1+(grid(i)-x1)*(y0-y1)/(x0-x1)
		endif
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	goto 100
102	continue
	do j=i,n
		y(j)=y1*x1**2/grid(j)**2
	enddo
	close(unit=20)
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
      SUBROUTINE sort(arr,n)
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,temp
	if(n.le.1) return
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=l-1
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END


	subroutine sortw(array,warray,n)
	IMPLICIT NONE
	integer n
	real*8 array(n),warray(n)
	integer i,j
	real*8 temp,wtemp
	
	if(n.le.1) return

	do i=2,n
		temp=array(i)
		wtemp=warray(i)
		do j=i-1,1,-1
			if (array(j).le.temp) exit
			array(j+1)=array(j)
			warray(j+1)=warray(j)
		enddo
		array(j+1)=temp
		warray(j+1)=wtemp
	enddo
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine sortw_new(array,warray,last)
	IMPLICIT NONE
	integer i,j,left,right,last
	real*8 array(last),warray(last)
	real*8 temp,p,next,wtemp,wnext

	p=0.5*(array(1)+array(last))
	if (array(1).gt.array(last)) then
		temp=array(last)
		array(last)=array(1)
		array(1)=temp
		wtemp=warray(last)
		warray(last)=warray(1)
		warray(1)=wtemp
	endif

	left=1
	right=last
	temp=array(2)
	wtemp=warray(2)

	do i=2,last-1
		if (temp.lt.p) then
			do j=left,1,-1
				if (array(j).le.temp) exit
				array(j+1)=array(j)
				warray(j+1)=warray(j)
			enddo
			array(j+1)=temp
			warray(j+1)=wtemp
			temp=array(left+2)
			wtemp=warray(left+2)
			left=left+1
		else
			next=array(right-1)
			wnext=warray(right-1)
			do j=right,last
				if (array(j).ge.temp) exit
				array(j-1)=array(j)
				warray(j-1)=warray(j)
			enddo
			array(j-1)=temp
			warray(j-1)=wtemp
			temp=next
			wtemp=wnext
			right=right-1	
		endif
	enddo

	return
	end



	
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END


      FUNCTION expint(n,x)
	IMPLICIT NONE
      INTEGER n,MAXIT
      REAL*8 expint,x,EPS,FPMIN,EULER
      PARAMETER (MAXIT=100,EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649)
      INTEGER i,ii,nm1
      REAL*8 a,b,c,d,del,fact,h,psi
      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0d0.and.(n.eq.0.or.n.eq.1)))then
        print*,'bad arguments in expint'
        expint=1./real(nm1)
      else if(n.eq.0)then
        expint=exp(-x)/x
      else if(x.eq.0d0)then
        expint=1./real(nm1)
      else if(x.gt.1d0)then
        b=x+n
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=1,MAXIT
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.).lt.EPS)then
            expint=h*exp(-x)
            return
          endif
11      continue
c        stop 'continued fraction failed in expint'
c        print*,'continued fraction failed in expint'
		expint=0d0
      else
        if(nm1.ne.0)then
          expint=1./nm1
        else
          expint=-log(x)-EULER
        endif
        fact=1.
        do 13 i=1,MAXIT
          fact=-fact*x/i
          if(i.ne.nm1)then
            del=-fact/(i-nm1)
          else
            psi=-EULER
            do 12 ii=1,nm1
              psi=psi+1./ii
12          continue
            del=fact*(-log(x)+psi)
          endif
          expint=expint+del
          if(abs(del).lt.abs(expint)*EPS) return
13      continue
c        stop 'series failed in expint'
c        print*,'series failed in expint'
		expint=0d0
      endif
      return
      END


	
	subroutine quadint(x1,x2,x3,y1,y2,y3,a,b,c)
	IMPLICIT NONE
	real*8 x1,x2,x3,y1,y2,y3,a,b,c
	real*8 y12,y23,x12,x23
	
	y12=(y1-y2)/(x1-x2)
	y23=(y2-y3)/(x2-x3)
	x12=(x1**2-x2**2)/(x1-x2)
	x23=(x2**2-x3**2)/(x2-x3)
	
	a=(y12-y23)/(x12-x23)
	b=y12-a*x12
	c=y1-b*x1-a*x1**2
	
	return
	end


      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
	
	
	
	subroutine stats(dat,n,aver,var)
	IMPLICIT NONE
	integer n,i
	real*8 dat(n),aver,var
	aver=0d0
	var=0d0
	do i=1,n
		aver=aver+dat(i)
	enddo
	aver=aver/real(n)
	do i=1,n
		var=var+(aver-dat(i))**2
	enddo
	var=sqrt(var/real(n-1))
	return
	end


	subroutine stati(dat,n,aver,var1,var2)
	IMPLICIT NONE
	integer n,i,n1,n2
	real*8 dat(n),aver,var1,var2
	aver=0d0
	var1=0d0
	var2=0d0
	n1=0
	n2=0
	do i=1,n
		aver=aver+dat(i)
	enddo
	aver=aver/real(n)
	do i=1,n
		if(dat(i).lt.aver) then
			var1=var1+(aver-dat(i))**2
			n1=n1+1
		else
			var2=var2+(aver-dat(i))**2
			n2=n2+1
		endif
	enddo
	var1=sqrt(var1/real(n1-1))
	var2=sqrt(var2/real(n2-1))
	if(n1.le.1) var1=0d0
	if(n2.le.1) var2=0d0
	return
	end


	subroutine writeppmfile(filename,im,nlam,n)
	IMPLICIT NONE
	character*500 filename
	integer n,nlam,cdepth,i,j,k,c
	real*8 im(n,n,nlam)

	cdepth=256

	open(unit=35,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	write(35,'("P3")')
	write(35,'("# ",a)') trim(filename)
	write(35,'(i5,i5)') n,n
	write(35,'(i5)') cdepth
	do i=1,n
		do j=1,n
			do k=1,nlam
				c=im(j,i,k)*real(cdepth)
				write(35,*) c
			enddo
		enddo
	enddo
	close(unit=35)
	return
	end
	
		
	subroutine computemedian(x,n,xm)
	IMPLICIT NONE
	integer i,n
	real*8 x(n),xm

	if(n.eq.1) then
		xm=x(1)
		return
	endif

	call sort(x,n)
	if((n/2)*2.eq.n) then
		i=n/2
		xm=(x(i)+x(i+1))/2d0
	else
		i=n/2+1
		xm=x(i)
	endif
	
	return
	end
	
	subroutine computeav50(x,n,xm)
	IMPLICIT NONE
	integer i,n,i1,i2
	real*8 x(n),xm,tot,w,sig,c
	
	if(n.eq.1) then
		xm=x(1)
		return
	endif

	call sort(x,n)
	xm=0d0
	tot=0d0
	c=real(n)/2d0+0.5
	sig=real(n)/8d0
	if(n.lt.5) then
		i1=1
		i2=n
	else if(n.lt.9) then
		i1=2
		i2=n-1
	else if(n.lt.11) then
		i1=3
		i2=n-2
	else
		i1=4
		i2=n-3
	endif
	do i=i1,i2
		w=exp(-((c-real(i))/sig)**2)
		tot=tot+w
		xm=xm+w*x(i)
	enddo
	xm=xm/tot
	
	return
	end

	subroutine printstats(x,n)
	IMPLICIT NONE
	integer i,n
	real*8 x(n),dx,x0

	if(n.le.1) return	
	x0=sum(x(1:n))/real(n)
	dx=0d0
	do i=1,n
		dx=dx+(x(n)-x0)**2
	enddo
	dx=sqrt(dx/(real(n-1)))
	print*,dx/x0
	
	return
	end
	
