
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
      if (iset.eq.0) then
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
	open(unit=20,file=input,RECL=100)
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
	open(unit=20,file=input,RECL=100)
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
c-----------------------------------------------------------------------
	
	subroutine sort(x,n)
	IMPLICIT NONE
	integer n,i,j,imin
	real*8 x(n),min
	
	do j=1,n-1
	min=x(j)
	imin=j
	do i=j,n
		if(x(i).lt.min) then
			min=x(i)
			imin=i
		endif
	enddo
	min=x(j)
	x(j)=x(imin)
	x(imin)=min
	enddo
	
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	
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
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
        print*,'bad arguments in expint'
        expint=1./nm1
      else if(n.eq.0)then
        expint=exp(-x)/x
      else if(x.eq.0.)then
        expint=1./nm1
      else if(x.gt.1.)then
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
      endif
      return
      END


