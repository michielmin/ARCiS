      subroutine voigt(a, v, ans)
      implicit none

! computes the voight function per zaghloul mnras, 375, 1043, 2007
! this retouine also returns the derivaties of the voigt functions 
! with respect to its arguments

! input: 
! a - the ratio of the natural width to the doppler width 
!     usually called "y" in the literature
! v - distance from line center in units of the doppler width
!     usually called "x" in the literature


! output:
! ans - the voigt function, the amplitude of the line profile
! dansda - derivative of the voigt function with respect to a
! dansdv - derivative of the voigt function with respect to v


! declare the pass
      double precision a,v,ans,dansda,dansdv


! common block to the integrands
      double precision  a_ratio,vat
      common /voigtcom/ a_ratio,vat


! local variables
      integer           i,imax
      parameter         (imax = 10000)
      external          zaghloul,dzaghloulda,dzaghlouldv
      double precision  zaghloul,dzaghloulda,dzaghlouldv
      double precision  term1,term2,ulo,uhi,
     1                  t1,dt1da,ff,gg,hh,t2,dt2da,dt2dv,
     2                  dterm1da,dterm1dv,dterm2da,dterm2dv,
     3                  sum,half_period,zero_cross


! for the expansion when a is greater than 26.6
      double precision  ainv2,c2,c4,c6,c8,c10,c12
      parameter         (c2  = 1.0d0/2.0d0,    
     1                   c4  = 3.0d0/4.0d0,    
     2                   c6  = 15.0d0/8.0d0,    
     3                   c8  = 105.0d0/16.0d0,  
     4                   c10 = 945.0d0/32.0d0, 
     5                   c12 = 10395.0d0/64.0d0)


! note rtpi = sqrt(pi), f
! whm = 2*sqrt(ln(2)) = full width half max of gaussian with variance 1/rtpi
      double precision  pi,rtpi,factor,fwhm,zero,eps,eps2
      parameter         (pi     = 3.1415926535897932384d0,
     1                   rtpi   = 1.7724538509055159d0,    
     2                   factor = 2.0d0/rtpi,
     3                   fwhm   = 1.6651d0,
     4                   zero   = 0.0d0,                    
     5                   eps    = 1.0d-12,
     6                   eps2   = 1.0d-16)


! formats statements
 78   format(1x,i6,1p6e14.6)

	if(abs(v).gt.8d0) then
		ans=a/(v**2+a**2)/sqrt(3.1415926536)
		return
	endif

! load the common block with the passed values
      a_ratio = a
      vat     = v


! the first term of zaghloul's equation 4, using equation 6 as needed,
! with derivatives
! here i rely on erfc being an intrinsic f90 function

      if (a_ratio .lt. 26.6) then
       t1    = exp(a_ratio*a_ratio) * erfc(a_ratio) 
       dt1da = 2.0d0*a_ratio*t1 - factor

      else
       ainv2  = 1.0d0/(a_ratio * a_ratio)
       t1     = 1.0d0/(rtpi * a_ratio) * 
     &          (1.0d0 - ainv2*(c2 - ainv2*(c4 - ainv2*(c6  
     &          - ainv2*(c8 - ainv2*(c10 - ainv2*c12))))))
       dt1da  = -ainv2/rtpi * 
     &          (1.0d0 - ainv2*(3.0d0*c2 - ainv2*(5.0d0*c4 
     &           - ainv2*(7.0d0*c6 - ainv2*(9.0d0*c8 
     &          - ainv2*(11.0d0*c10 - ainv2*13.0d0*c12))))))

      end if

      ff       = exp(-(v*v))
      gg       = cos(2.0d0*a_ratio*v)
      hh       = sin(2.0d0*a_ratio*v)

      term1    = t1 * ff * gg
!      write(6,78) 0,term1,dterm1da,dterm1dv



! the second term of zaghloul's eq 4
! the sinusoids in the integrands have a period of pi/a_ratio.
! we'll be integrating over each half-period so that all
! contributions have the same sign if the half-period.
! we'll also be intergrating over the full-width-half-max
! of the gaussian. we'll take the smallest length scale
! as the one to integrate over.


      half_period = min(0.5d0 * pi/a_ratio, fwhm)


! first the voigt function itself
! initialize the integral sum and the upper integration limit
      sum  = 0.0d0
      uhi  = vat

! form a lower integration bound, 
! allow for positive of negative values of v

      do i=1,imax
       ulo = uhi - sign(1.0d0,uhi)*half_period
       if (uhi .eq. 0.0d0) then
        ulo = 0.0d0
       else if (uhi .lt. 0.0) then
        ulo = min(zero,ulo)
       else if (uhi .gt. 0.0) then
        ulo = max(zero,ulo)
       end if

! pop each integral
       call qromb_zag(zaghloul,ulo,uhi,eps,t2)
!       write(6,78) i,ulo,uhi,t2

! add in these contributions to the total
       sum = sum + t2
       
! convergence check
       if (ulo .eq. 0.0 .or. abs(t2/sum) .le. eps2) goto 100

! swap the limits for the next half period, and end of loop
       uhi = ulo
      enddo

! write a warning if we excedded the maximum number of half-periods
      write(6,*) 'voigt integral did not converge'

! arrive here after convergence
 100  continue
      term2    = sum * factor
!      read(5,*)



! voigt function and its derivatives are the sum of the two terms

      ans    = term1 + term2  

      return
      end





      double precision function zaghloul(u)
      implicit none
      save

! this routine forms the integrand of zaghloul equation 4

! declare the pass
      double precision u

! common block to the voigt integrand
      double precision  a_ratio,vat
      common /voigtcom/ a_ratio,vat

! go
      zaghloul = exp(-((vat*vat) - (u*u))) * sin(2.0d0*a_ratio*(vat-u))

!      write(6,77) a_ratio,vat,u,zaghloul
! 77   format(1x,1p4e14.6)


      return
      end




      double precision function dzaghloulda(u)
      implicit none
      save

! this routine forms the integrand of zaghloul equation 4

! declare the pass
      double precision u

! common block to the voigt integrand
      double precision  a_ratio,vat
      common /voigtcom/ a_ratio,vat

! go
      dzaghloulda = exp(-((vat*vat) - (u*u))) 
     1              * cos(2.0d0*a_ratio*(vat-u)) * 2.0d0*(vat-u)

      return
      end




      double precision function dzaghlouldv(u)
      implicit none
      save

! this routine forms the integrand of zaghloul equation 4

! declare the pass
      double precision u

! local variables
      double precision ff,gg

! common block to the voigt integrand
      double precision  a_ratio,vat
      common /voigtcom/ a_ratio,vat

! go
      ff = exp(-((vat*vat) - (u*u)))
      gg = 2.0d0*a_ratio*(vat - u)
      dzaghlouldv = 2.0d0* ff * (cos(gg) * a_ratio - sin(gg) * vat)

      return
      end




      subroutine qromb_zag(func,a,b,eps,ss)  
      implicit none
      save

! returns as ss the integral of the function func from a to b with fractional 
! accuracy eps. integration by romberg's method of order 2k where e.g k=2 is  
! simpson's rule.  
!   
! jmax limits the total number of steps; k is the 
! the number of points used in the extrapolation; arrays s and h store the  
! trapazoidal approximations and their relative step sizes. 
!  
! declare 
      external          func 
      integer           j,jmax,jmaxp,k,km    
      parameter         (jmax=20, jmaxp=jmax+1, k=5, km=k-1) 
      double precision  a,b,ss,s(jmaxp),h(jmaxp),eps,dss,func 

      h(1) = 1.0d0 
      do j=1,jmax 
       call trapzd_zag(func,a,b,s(j),j)  
       if (j .ge. k) then    
        call polint_zag(h(j-km),s(j-km),k,0.0d0,ss,dss)    
        if (abs(dss) .le. eps*abs(ss)) return    
       end if    
       s(j+1) = s(j) 
       h(j+1) = 0.25d0 * h(j)  
      enddo
!      write(6,*) ' after ',jmax,' iterations ' 
!      write(6,*) ' of trying to integrate between ',a,' and ',b 
!      write(6,*) ' with fractional accuracy ',eps 
!      write(6,*) ' the integral is ',ss 
!      write(6,*) ' and error estimate ',dss 
!      write(6,*) ' so abs(dss) ',abs(dss),'> eps*abs(ss)',eps*abs(ss) 
!      stop       'too many steps in qromb'  
      return 
      end    




      subroutine trapzd_zag(func,a,b,s,n)    
      implicit none
      save

! this routine computes the n'th stage of refinement of an extended 
! trapazoidal rule. func is input as the name of a function to be   
! integrated between limits a and b. when n=1 the routine returns as s  
! the crudest estimate of the integral of func(x)dx from a to b.    
! subsequent calls with n=2,3... will improve the accuracy of s by adding   
! 2**(n-2) additional interior points. s should not be modified between 
! sequential calls. 
!  
! this routine is the workhorse of all the following closed formula 
! integration routines. 
!   
! local it  is the number of points to be added on the next call    
! local del is the step size.   
!  
! declare  
      external          func 
      integer           n,it,j   
      double precision  func,a,b,s,del,x,sum,tnm 

! go 
      if (n.eq.1) then   
       s  = 0.5d0 * (b-a) * ( func(a) + func(b) )   
      else   
       it  = 2**(n-2) 
       tnm = it  
       del = (b-a)/tnm   
       x   = a + (0.5d0 *del)    
       sum = 0.0d0 
       do j=1,it  
        sum = sum + func(x)  
        x   = x + del  
       enddo
       s  = 0.5d0 * (s + (b-a)*sum/tnm) 
      end if 
      return     
      end    





      subroutine polint_zag(xa,ya,n,x,y,dy)
      implicit none
      save

! given arrays xa and ya of length n and a value x, this routine returns a 
! value y and an error estimate dy. if p(x) is the polynomial of degree n-1
! such that ya = p(xa) ya then the returned value is y = p(x) 
! 
! declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=10)
      double precision xa(n),ya(n),x,y,dy,c(nmax),d(nmax),dif,dift,
     1                 ho,hp,w,den

! find the index ns of the closest table entry; initialize the c and d tables
      ns  = 1
      dif = abs(x - xa(1))
      do i=1,n
       dift = abs(x - xa(i))
       if (dift .lt. dif) then
        ns  = i
        dif = dift
       end if
       c(i)  = ya(i)
       d(i)  = ya(i)
      enddo

! first guess for y
      y = ya(ns)

! for each column of the table, loop over the c's and d's and update them
      ns = ns - 1
      do m=1,n-1
       do i=1,n-m
        ho   = xa(i) - x
        hp   = xa(i+m) - x
        w    = c(i+1) - d(i)
        den  = ho - hp
        if (den .eq. 0.0) stop ' 2 xa entries are the same in polint'
        den  = w/den
        d(i) = hp * den
        c(i) = ho * den
       enddo

! after each column is completed, decide which correction c or d, to add
! to the accumulating value of y, that is, which path to take in the table
! by forking up or down. ns is updated as we go to keep track of where we
! are. the last dy added is the error indicator.
       if (2*ns .lt. n-m) then
        dy = c(ns+1)
       else
        dy = d(ns)
        ns = ns - 1
       end if
       y = y + dy
      enddo
      return
      end

