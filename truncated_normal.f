      function i4_uniform_ab ( a, b, seed )

c*********************************************************************72
c
cc I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer I4_UNIFORM_AB, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer i4_uniform_ab
      integer k
      real r
      integer seed
      integer value

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop 1
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if

      r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
      r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &  +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
      value = nint ( r )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      i4_uniform_ab = value

      return
      end
      subroutine normal_01_cdf ( x, cdf )

c*********************************************************************72
c
cc NORMAL_01_CDF evaluates the Normal 01 CDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    AG Adams,
c    Algorithm 39,
c    Areas Under the Normal Curve,
c    Computer Journal,
c    Volume 12, pages 197-198, 1969.
c
c  Parameters:
c
c    Input, double precision X, the argument of the CDF.
c
c    Output, double precision CDF, the value of the CDF.
c
      implicit none

      double precision a1
      parameter ( a1 = 0.398942280444D+00 )
      double precision a2
      parameter ( a2 = 0.399903438504D+00 )
      double precision, parameter :: a3 = 5.75885480458D+00
      double precision, parameter :: a4 = 29.8213557808D+00
      double precision, parameter :: a5 = 2.62433121679D+00
      double precision, parameter :: a6 = 48.6959930692D+00
      double precision, parameter :: a7 = 5.92885724438D+00
      double precision, parameter :: b0 = 0.398942280385D+00
      double precision, parameter :: b1 = 3.8052D-08
      double precision, parameter :: b2 = 1.00000615302D+00
      double precision, parameter :: b3 = 3.98064794D-04
      double precision, parameter :: b4 = 1.98615381364D+00
      double precision, parameter :: b5 = 0.151679116635D+00
      double precision, parameter :: b6 = 5.29330324926D+00
      double precision, parameter :: b7 = 4.8385912808D+00
      double precision, parameter :: b8 = 15.1508972451D+00
      double precision, parameter :: b9 = 0.742380924027D+00
      double precision, parameter :: b10 = 30.789933034D+00
      double precision b11
      parameter ( b11 = 3.99019417011D+00 )
      double precision cdf
      double precision q
      double precision x
      double precision y
c
c  |X| .le. 1.28.
c
      if ( abs ( x ) .le. 1.28D+00 ) then

        y = 0.5D+00 * x * x

        q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 
     &    / ( y + a5 + a6 / ( y + a7 ) ) ) )
c
c  1.28 .lt. |X| .le. 12.7
c
      else if ( abs ( x ) .le. 12.7D+00 ) then

        y = 0.5D+00 * x * x

        q = exp ( - y ) * b0 / ( abs ( x ) - b1
     &    + b2  / ( abs ( x ) + b3
     &    + b4  / ( abs ( x ) - b5
     &    + b6  / ( abs ( x ) + b7 
     &    - b8  / ( abs ( x ) + b9
     &    + b10 / ( abs ( x ) + b11 ) ) 
     &) ) ) )
c
c  12.7 .lt. |X|
c
      else

        q = 0.0D+00

      end if
c
c  Take account of negative X.
c
      if ( x .lt. 0.0D+00 ) then
        cdf = q
      else
        cdf = 1.0D+00 - q
      end if

      return
      end
      subroutine normal_01_cdf_inv ( p, x )

c*********************************************************************72
c
cc NORMAL_01_CDF_INV inverts the standard normal CDF.
c
c  Discussion:
c
c    The result is accurate to about 1 part in 10**16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 February 2015
c
c  Author:
c
c    Original FORTRAN77 version by Michael Wichura.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Michael Wichura,
c    Algorithm AS241:
c    The Percentage Points of the Normal Distribution,
c    Applied Statistics,
c    Volume 37, Number 3, pages 477-484, 1988.
c
c  Parameters:
c
c    Input, double precision P, the value of the cumulative probability
c    densitity function.  0 .lt. P .lt. 1.  If P is outside this range, an
c    "infinite" value will be returned.
c
c    Output, double precision X, the normal deviate value
c    with the property that the probability of a standard normal deviate being
c    less than or equal to the value is P.
c
      implicit none

      double precision a(8)
      double precision b(8)
      double precision c(8)
      double precision const1
      parameter ( const1 = 0.180625D+00 )
      double precision const2
      parameter ( const2 = 1.6D+00 )
      double precision d(8)
      double precision e(8)
      double precision f(8)
      double precision p
      double precision q
      double precision r
      double precision r8_huge
      double precision r8poly_value_horner
      double precision split1
      parameter ( split1 = 0.425D+00 )
      double precision split2
      parameter ( split2 = 5.0D+00 )
      double precision x

      save a
      save b
      save c
      save d
      save e
      save f

      data a /
     &  3.3871328727963666080D+00,
     &  1.3314166789178437745D+02,
     &  1.9715909503065514427D+03,
     &  1.3731693765509461125D+04,
     &  4.5921953931549871457D+04,     
     &  6.7265770927008700853D+04,     
     &  3.3430575583588128105D+04,     
     &  2.5090809287301226727D+03 /
      data b /
     &  1.0D+00,
     &  4.2313330701600911252D+01,
     &  6.8718700749205790830D+02,  
     &  5.3941960214247511077D+03,
     &  2.1213794301586595867D+04,
     &  3.9307895800092710610D+04,
     &  2.8729085735721942674D+04,
     &  5.2264952788528545610D+03 /
      data c /
     &  1.42343711074968357734D+00,
     &  4.63033784615654529590D+00,
     &  5.76949722146069140550D+00,
     &  3.64784832476320460504D+00,
     &  1.27045825245236838258D+00,
     &  2.41780725177450611770D-01,
     &  2.27238449892691845833D-02,
     &  7.74545014278341407640D-04 /
      data d /
     &  1.0D+00,
     &  2.05319162663775882187D+00,
     &  1.67638483018380384940D+00,
     &  6.89767334985100004550D-01,
     &  1.48103976427480074590D-01,  
     &  1.51986665636164571966D-02,
     &  5.47593808499534494600D-04,    
     &  1.05075007164441684324D-09 /
      data e /
     &  6.65790464350110377720D+00,
     &  5.46378491116411436990D+00,
     &  1.78482653991729133580D+00,
     &  2.96560571828504891230D-01,
     &  2.65321895265761230930D-02,
     &  1.24266094738807843860D-03,
     &  2.71155556874348757815D-05,
     &  2.01033439929228813265D-07 /
      data f /
     &  1.0D+00,
     &  5.99832206555887937690D-01,
     &  1.36929880922735805310D-01,
     &  1.48753612908506148525D-02,
     &  7.86869131145613259100D-04,  
     &  1.84631831751005468180D-05,
     &  1.42151175831644588870D-07,    
     &  2.04426310338993978564D-15 /

      if ( p .le. 0.0D+00 ) then
        x = - r8_huge ( )
        return
      end if

      if ( 1.0D+00 .le. p ) then
        x = r8_huge ( )
        return
      end if

      q = p - 0.5D+00

      if ( abs ( q ) .le. split1 ) then

        r = const1 - q * q
        x = q * r8poly_value_horner ( 7, a, r ) 
     &        / r8poly_value_horner ( 7, b, r )

      else

        if ( q .lt. 0.0D+00 ) then
          r = p
        else
          r = 1.0D+00 - p
        end if

        if ( r .le. 0.0D+00 ) then

          x = r8_huge ( )

        else

          r = sqrt ( - log ( r ) )

          if ( r .le. split2 ) then

            r = r - const2
            x = r8poly_value_horner ( 7, c, r ) 
     &        / r8poly_value_horner ( 7, d, r )

          else

            r = r - split2
            x = r8poly_value_horner ( 7, e, r ) 
     &        / r8poly_value_horner ( 7, f, r )

          end if

        end if

        if ( q .lt. 0.0D+00 ) then
          x = -x
        end if

      end if

      return
      end
      subroutine normal_01_cdf_values ( n_data, x, fx )

c*********************************************************************72
c
cc NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      Needs["Statistics`ContinuousDistributions`"]
c      dist = NormalDistribution [ 0, 1 ]
c      CDF [ dist, x ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 17 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  0.5000000000000000D+00,
     &  0.5398278372770290D+00,
     &  0.5792597094391030D+00,
     &  0.6179114221889526D+00,
     &  0.6554217416103242D+00,
     &  0.6914624612740131D+00,
     &  0.7257468822499270D+00,
     &  0.7580363477769270D+00,
     &  0.7881446014166033D+00,
     &  0.8159398746532405D+00,
     &  0.8413447460685429D+00,
     &  0.9331927987311419D+00,
     &  0.9772498680518208D+00,
     &  0.9937903346742239D+00,
     &  0.9986501019683699D+00,
     &  0.9997673709209645D+00,
     &  0.9999683287581669D+00 /
      data x_vec /
     &  0.0000000000000000D+00,
     &  0.1000000000000000D+00,
     &  0.2000000000000000D+00,
     &  0.3000000000000000D+00,
     &  0.4000000000000000D+00,
     &  0.5000000000000000D+00,
     &  0.6000000000000000D+00,
     &  0.7000000000000000D+00,
     &  0.8000000000000000D+00,
     &  0.9000000000000000D+00,
     &  0.1000000000000000D+01,
     &  0.1500000000000000D+01,
     &  0.2000000000000000D+01,
     &  0.2500000000000000D+01,
     &  0.3000000000000000D+01,
     &  0.3500000000000000D+01,
     &  0.4000000000000000D+01 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine normal_01_mean ( mean )

c*********************************************************************72
c
cc NORMAL_01_MEAN returns the mean of the Normal 01 PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 December 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision MEAN, the mean of the PDF.
c
      implicit none

      double precision mean

      mean = 0.0D+00

      return
      end
      subroutine normal_01_moment ( order, value )

c*********************************************************************72
c
cc NORMAL_01_MOMENT evaluates moments of the Normal 01 PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ORDER, the order of the moment.
c    0 <= ORDER.
c
c    Output, double precision VALUE, the value of the moment.
c
      implicit none

      integer order
      double precision r8_factorial2
      double precision value

      if ( mod ( order, 2 ) .eq. 0 ) then
        value = r8_factorial2 ( order - 1 )
      else
        value = 0.0D+00
      end if

      return
      end
      subroutine normal_01_pdf ( x, pdf )

c*********************************************************************72
c
cc NORMAL_01_PDF evaluates the Normal 01 PDF.
c
c  Discussion:
c
c    The Normal 01 PDF is also called the "Standard Normal" PDF, or
c    the Normal PDF with 0 mean and standard deviation 1.
c
c    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 December 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the PDF.
c
c    Output, double precision PDF, the value of the PDF.
c
      implicit none

      double precision pdf
      double precision r8_pi
      parameter ( r8_pi = 3.141592653589793D+00 )
      double precision x

      pdf = exp ( -0.5D+00 * x * x ) / sqrt ( 2.0D+00 * r8_pi )

      return
      end
      subroutine normal_01_sample ( seed, x )

c*********************************************************************72
c
cc NORMAL_01_SAMPLE samples the standard normal probability distribution.
c
c  Discussion:
c
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c    The Box-Muller method is used, which is efficient, but
c    generates two values at a time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 February 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X, a sample of the standard normal PDF.
c
      implicit none

      double precision r1
      double precision r2
      double precision r8_pi
      parameter ( r8_pi = 3.141592653589793D+00 )
      double precision r8_uniform_01
      integer seed
      integer used
      double precision x
      double precision y

      save used
      save y

      data used / -1 /
      data y / 0.0D+00 /

      if ( used .eq. -1 ) then
        used = 0
      end if
c
c  If we've used an even number of values so far, generate two more,
c  return one and save one.
c
      if ( mod ( used, 2 ) .eq. 0 ) then

10      continue

          r1 = r8_uniform_01 ( seed )

          if ( r1 .ne. 0.0D+00 ) then
            go to 20
          end if

        go to 10

20      continue

        r2 = r8_uniform_01 ( seed )

        x = sqrt ( - 2.0D+00 * log ( r1 ) ) 
     &    * cos ( 2.0D+00 * r8_pi * r2 )
        y = sqrt ( - 2.0D+00 * log ( r1 ) ) 
     &    * sin ( 2.0D+00 * r8_pi * r2 )
c
c  Otherwise, return the second, saved, value.
c
      else

        x = y

      end if

      used = used + 1

      return
      end
      subroutine normal_01_variance ( variance )

c*********************************************************************72
c
cc NORMAL_01_VARIANCE returns the variance of the Normal 01 PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 December 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision VARIANCE, the variance of the PDF.
c
      implicit none

      double precision variance

      variance = 1.0D+00

      return
      end
      subroutine normal_ms_cdf ( x, mu, sigma, cdf )

c*********************************************************************72
c
cc NORMAL_MS_CDF evaluates the Normal CDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 February 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the CDF.
c
c    Input, double precision MU, SIGMA, the parameters of the PDF.
c    0.0D+00 .lt. SIGMA.
c
c    Output, double precision CDF, the value of the CDF.
c
      implicit none

      double precision cdf
      double precision mu
      double precision sigma
      double precision x
      double precision y

      y = ( x - mu ) / sigma

      call normal_01_cdf ( y, cdf )

      return
      end
      subroutine normal_ms_cdf_inv ( cdf, mu, sigma, x )

c*********************************************************************72
c
cc NORMAL_MS_CDF_INV inverts the Normal CDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 February 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision CDF, the value of the CDF.
c    0.0D+00 .le. CDF .le. 1.0.
c
c    Input, double precision MU, SIGMA, the parameters of the PDF.
c    0.0D+00 .lt. SIGMA.
c
c    Output, double precision X, the corresponding argument.
c
      implicit none

      double precision cdf
      double precision mu
      double precision sigma
      double precision x
      double precision x2

      if ( cdf .lt. 0.0D+00 .or. 1.0D+00 .lt. cdf ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NORMAL_MS_CDF_INV - Fatal error!'
        write ( *, '(a)' ) '  CDF .lt. 0 or 1 .lt. CDF.'
        stop 1
      end if

      call normal_01_cdf_inv ( cdf, x2 )

      x = mu + sigma * x2

      return
      end
      subroutine normal_ms_mean ( mu, sigma, mean )

c*********************************************************************72
c
cc NORMAL_MS_MEAN returns the mean of the Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the parameters of the PDF.
c    0.0D+00 .lt. SIGMA.
c
c    Output, double precision MEAN, the mean of the PDF.
c
      implicit none

      double precision mean
      double precision mu
      double precision sigma

      mean = mu

      return
      end
      subroutine normal_ms_moment ( order, mu, sigma, value )

c*********************************************************************72
c
cc NORMAL_MS_MOMENT evaluates moments of the Normal PDF.
c
c  Discussion:
c
c    The formula was posted by John D Cook.
c
c    Order  Moment
c    -----  ------
c      0    1
c      1    mu
c      2    mu^2 +         sigma^2
c      3    mu^3 +  3 mu   sigma^2
c      4    mu^4 +  6 mu^2 sigma^2 +   3      sigma^4
c      5    mu^5 + 10 mu^3 sigma^2 +  15 mu   sigma^4
c      6    mu^6 + 15 mu^4 sigma^2 +  45 mu^2 sigma^4 +  15      sigma^6
c      7    mu^7 + 21 mu^5 sigma^2 + 105 mu^3 sigma^4 + 105 mu   sigma^6
c      8    mu^8 + 28 mu^6 sigma^2 + 210 mu^4 sigma^4 + 420 mu^2 sigma^6 
c           + 105 sigma^8
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ORDER, the order of the moment.
c    0 <= ORDER.
c
c    Input, double precision MU, the mean of the distribution.
c
c    Input, double precision SIGMA, the standard deviation of the distribution.
c
c    Output, double precision VALUE, the value of the central moment.
c
      implicit none

      integer j
      integer j_hi
      double precision mu
      integer order
      double precision r8_choose
      double precision r8_factorial2
      double precision sigma
      double precision value

      j_hi = ( order / 2 )

      value = 0.0D+00
      do j = 0, j_hi
        value = value 
     &    + r8_choose ( order, 2 * j ) 
     &    * r8_factorial2 ( 2 * j - 1 ) 
     &    * mu ** ( order - 2 * j ) * sigma ** ( 2 * j )
      end do

      return
      end
      subroutine normal_ms_moment_central ( order, mu, sigma, value )

c*********************************************************************72
c
cc NORMAL_MS_MOMENT_CENTRAL evaluates central moments of the Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ORDER, the order of the moment.
c    0 <= ORDER.
c
c    Input, double precision MU, the mean of the distribution.
c
c    Input, double precision SIGMA, the standard deviation of the distribution.
c
c    Output, double precision VALUE, the value of the central moment.
c
      implicit none

      double precision mu
      integer order
      double precision r8_factorial2
      double precision sigma
      double precision value

      if ( mod ( order, 2 ) .eq. 0 ) then
        value = r8_factorial2 ( order - 1 ) * sigma ** order
      else
        value = 0.0D+00
      end if

      return
      end
      subroutine normal_ms_moment_central_values ( order, mu, sigma, 
     &  value )

c*********************************************************************72
c
cc NORMAL_MS_MOMENT_CENTRAL_VALUES: moments 0 through 10 of the Normal PDF.
c
c  Discussion:
c
c    The formula was posted by John D Cook.
c
c    Order  Moment
c    -----  ------
c      0    1
c      1    0
c      2    sigma^2
c      3    0
c      4    3 sigma^4
c      5    0
c      6    15 sigma^6
c      7    0
c      8    105 sigma^8
c      9    0
c     10    945 sigma^10
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ORDER, the order of the moment.
c    0 <= ORDER <= 10.
c
c    Input, double precision MU, the mean of the distribution.
c
c    Input, double precision SIGMA, the standard deviation of the distribution.
c
c    Output, double precision VALUE, the value of the central moment.
c
      implicit none

      double precision mu
      integer order
      double precision sigma
      double precision value

      if ( order .eq. 0 ) then
        value = 1.0D+00
      else if ( order .eq. 1 ) then
        value = 0.0D+00
      else if ( order .eq. 2 ) then
        value = sigma ** 2
      else if ( order .eq. 3 ) then
        value = 0.0D+00
      else if ( order .eq. 4 ) then
        value = 3.0D+00 * sigma ** 4
      else if ( order .eq. 5 ) then
        value = 0.0D+00
      else if ( order .eq. 6 ) then
        value = 15.0D+00 * sigma ** 6
      else if ( order .eq. 7 ) then
        value = 0.0D+00
      else if ( order .eq. 8 ) then
        value = 105.0D+00 * sigma ** 8
      else if ( order .eq. 9 ) then
        value = 0.0D+00
      else if ( order .eq. 10 ) then
        value = 945.0D+00 * sigma ** 10
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    'NORMAL_MS_MOMENT_CENTRAL_VALUES - Fatal error!'
        write ( *, '(a)' ) '  Only ORDERS 0 through 10 are available.'
        stop 1
      end if

      return
      end
      subroutine normal_ms_moment_values ( order, mu, sigma, value )

c*********************************************************************72
c
cc NORMAL_MS_MOMENT_VALUES evaluates moments 0 through 8 of the Normal PDF.
c
c  Discussion:
c
c    The formula was posted by John D Cook.
c
c    Order  Moment
c    -----  ------
c      0    1
c      1    mu
c      2    mu^2 +         sigma^2
c      3    mu^3 +  3 mu   sigma^2
c      4    mu^4 +  6 mu^2 sigma^2 +   3      sigma^4
c      5    mu^5 + 10 mu^3 sigma^2 +  15 mu   sigma^4
c      6    mu^6 + 15 mu^4 sigma^2 +  45 mu^2 sigma^4 +  15      sigma^6
c      7    mu^7 + 21 mu^5 sigma^2 + 105 mu^3 sigma^4 + 105 mu   sigma^6
c      8    mu^8 + 28 mu^6 sigma^2 + 210 mu^4 sigma^4 + 420 mu^2 sigma^6 
c           + 105 sigma^8
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ORDER, the order of the moment.
c    0 <= ORDER <= 8.
c
c    Input, double precision MU, the mean of the distribution.
c
c    Input, double precision SIGMA, the standard deviation of the distribution.
c
c    Output, double precision VALUE, the value of the central moment.
c
      implicit none

      double precision mu
      integer order
      double precision sigma
      double precision value

      if ( order .eq. 0 ) then
        value = 1.0D+00
      else if ( order .eq. 1 ) then
        value = mu
      else if ( order .eq. 2 ) then
        value = mu ** 2 + sigma ** 2
      else if ( order .eq. 3 ) then
        value = mu ** 3 + 3.0D+00 * mu * sigma ** 2
      else if ( order .eq. 4 ) then
        value = mu ** 4 + 6.0D+00 * mu ** 2 * sigma ** 2 
     &    + 3.0D+00 * sigma ** 4
      else if ( order .eq. 5 ) then
        value = mu ** 5 + 10.0D+00 * mu ** 3 * sigma ** 2 
     &    + 15.0D+00 * mu * sigma ** 4
      else if ( order .eq. 6 ) then
        value = mu ** 6 + 15.0D+00 * mu ** 4 * sigma ** 2 
     &    + 45.0D+00 * mu ** 2 * sigma ** 4 
     &    + 15.0D+00 * sigma ** 6
      else if ( order .eq. 7 ) then
        value = mu ** 7 + 21.0D+00 * mu ** 5 * sigma ** 2 
     &    + 105.0D+00 * mu ** 3 * sigma ** 4 
     &    + 105.0D+00 * mu * sigma ** 6
      else if ( order .eq. 8 ) then
        value = mu ** 8 + 28.0D+00 * mu ** 6 * sigma ** 2 
     &    + 210.0D+00 * mu ** 4 * sigma ** 4 
     &    + 420.0D+00 * mu ** 2 * sigma ** 6 + 105.0D+00 * sigma ** 8
      else
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'NORMAL_MS_MOMENT_VALUES - Fatal error!'
        write ( *, '(a)' ) '  Only ORDERS 0 through 8 are available.'
        stop 1
      end if

      return
      end
      subroutine normal_ms_pdf ( x, mu, sigma, pdf )

c*********************************************************************72
c
cc NORMAL_MS_PDF evaluates the Normal PDF.
c
c  Discussion:
c
c    PDF(MU,SIGMA;X)
c      = exp ( - 0.5 * ( ( X - MU ) / SIGMA )^2 ) / ( SIGMA * sqrt ( 2 * PI ) )
c
c    The normal PDF is also known as the Gaussian PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the PDF.
c
c    Input, double precision MU, SIGMA, the parameters of the PDF.
c    0.0D+00 .lt. SIGMA.
c
c    Output, double precision PDF, the value of the PDF.
c
      implicit none

      double precision mu
      double precision pdf
      double precision r8_pi
      parameter ( r8_pi = 3.141592653589793D+00 )
      double precision sigma
      double precision x
      double precision y

      y = ( x - mu ) / sigma

      pdf = exp ( - 0.5D+00 * y * y ) 
     &  / ( sigma * sqrt ( 2.0D+00 * r8_pi ) )

      return
      end
      subroutine normal_ms_sample ( mu, sigma, seed, x )

c*********************************************************************72
c
cc NORMAL_MS_SAMPLE samples the Normal PDF.
c
c  Discussion:
c
c    The Box-Muller method is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 October 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the parameters of the PDF.
c    0.0D+00 .lt. SIGMA.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X, a sample of the PDF.
c
      implicit none

      double precision mu
      integer seed
      double precision sigma
      double precision x

      call normal_01_sample ( seed, x )

      x = mu + sigma * x

      return
      end
      subroutine normal_ms_variance ( mu, sigma, variance )

c*********************************************************************72
c
cc NORMAL_MS_VARIANCE returns the variance of the Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the parameters of the PDF.
c    0.0D+00 .lt. SIGMA.
c
c    Output, double precision VARIANCE, the variance of the PDF.
c
      implicit none

      double precision mu
      double precision sigma
      double precision variance

      variance = sigma * sigma

      return
      end
      function r8_choose ( n, k )

c*********************************************************************72
c
cc R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
c
c  Discussion:
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in R8 arithmetic.
c
c    The formula used is:
c
c      C(N,K) = N! / ( K! * (N-K)! )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    ML Wolfson, HV Wright,
c    Algorithm 160:
c    Combinatorial of M Things Taken N at a Time,
c    Communications of the ACM,
c    Volume 6, Number 4, April 1963, page 161.
c
c  Parameters:
c
c    Input, integer N, K, are the values of N and K.
c
c    Output, double precision R8_CHOOSE, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer k
      integer mn
      integer mx
      integer n
      double precision r8_choose
      double precision value

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        value = 0.0D+00

      else if ( mn .eq. 0 ) then

        value = 1.0D+00

      else

        mx = max ( k, n - k )
        value = dble ( mx + 1 )

        do i = 2, mn
          value = ( value * dble ( mx + i ) ) / dble ( i )
        end do

      end if

      r8_choose = value

      return
      end
      function r8_factorial2 ( n )

c*********************************************************************72
c
cc R8_FACTORIAL2 computes the double factorial function.
c
c  Discussion:
c
c    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
c                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
c
c  Example:
c
c     N Value
c
c     0     1
c     1     1
c     2     2
c     3     3
c     4     8
c     5    15
c     6    48
c     7   105
c     8   384
c     9   945
c    10  3840
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the double factorial
c    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
c
c    Output, double precision R8_FACTORIAL2, the value.
c
      implicit none

      integer n
      double precision r8_factorial2
      double precision r8_n

      if ( n .lt. 1 ) then
        r8_factorial2 = 1.0D+00
        return
      end if

      r8_n = dble ( n )
      r8_factorial2 = 1.0D+00

10    continue

      if ( 1.0D+00 .lt. r8_n ) then
        r8_factorial2 = r8_factorial2 * r8_n
        r8_n = r8_n - 2.0D+00
        go to 10
      end if

      return
      end
      subroutine r8_factorial2_values ( n_data, n, f )

c*********************************************************************72
c
cc R8_FACTORIAL2_VALUES returns values of the double factorial function.
c
c  Discussion:
c
c    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
c                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
c
c    In Mathematica, the function can be evaluated by:
c
c      n!!
c
c  Example:
c
c     N    N!!
c
c     0     1
c     1     1
c     2     2
c     3     3
c     4     8
c     5    15
c     6    48
c     7   105
c     8   384
c     9   945
c    10  3840
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2015
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996,
c    ISBN: 0-8493-2479-3,
c    LC: QA47.M315.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer N, the argument of the function.
c
c    Output, double precision F, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 16 )

      double precision f_vec(n_max)
      double precision f
      integer n_data
      integer n
      integer n_vec(n_max)

      save f_vec
      save n_vec

      data f_vec /
     &        1.0D+00,
     &        1.0D+00,
     &        2.0D+00,
     &        3.0D+00,
     &        8.0D+00,
     &       15.0D+00,
     &       48.0D+00,
     &      105.0D+00,
     &      384.0D+00,
     &      945.0D+00,
     &     3840.0D+00,
     &    10395.0D+00,
     &    46080.0D+00,
     &   135135.0D+00,
     &   645120.0D+00,
     &  2027025.0D+00 /
      data n_vec /
     &   0,
     &   1,  2,  3,  4,  5,
     &   6,  7,  8,  9, 10,
     &  11, 12, 13, 14, 15 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        f = 0.0D+00
      else
        n = n_vec(n_data)
        f = f_vec(n_data)
      end if

      return
      end
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
c
c  Discussion:
c
c    The value returned by this function is NOT required to be the
c    maximum representable R8.  This value varies from machine to machine,
c    from compiler to compiler, and may cause problems when being printed.
c    We simply want a "very large" but non-infinite number.
c
c    FORTRAN90 provides a built-in routine HUGE ( X ) that
c    can return the maximum representable number of the same datatype
c    as X, if that is what is really desired.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_HUGE, a huge number.
c
      implicit none

      double precision r8_huge

      r8_huge = 1.0D+30

      return
      end
      function r8_log_2 ( x )

c*********************************************************************72
c
cc R8_LOG_2 returns the logarithm base 2 of an R8.
c
c  Discussion:
c
c    value = Log ( |X| ) / Log ( 2.0 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number whose base 2 logarithm is desired.
c    X should not be 0.
c
c    Output, double precision R8_LOG_2, the logarithm base 2 of the absolute
c    value of X.  It should be true that |X| = 2**D_LOG_2.
c
      implicit none

      double precision r8_huge
      double precision r8_log_2
      double precision x

      if ( x .eq. 0.0D+00 ) then
        r8_log_2 = -r8_huge ( )
      else
        r8_log_2 = log ( abs ( x ) ) / log ( 2.0D+00 )
      end if

      return
      end
      function r8_mop ( i )

c*********************************************************************72
c
cc R8_MOP returns the I-th power of -1 as an R8.
c
c  Discussion:
c
c    An R8 is a double precision real value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the power of -1.
c
c    Output, double precision R8_MOP, the I-th power of -1.
c
      implicit none

      integer i
      double precision r8_mop

      if ( mod ( i, 2 ) .eq. 0 ) then
        r8_mop = + 1.0D+00
      else
        r8_mop = - 1.0D+00
      end if

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a pseudorandom R8 scaled to [0,1].
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r8_uniform_01 = seed / ( 2^31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop 1
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8poly_print ( n, a, title )

c*********************************************************************72
c
cc R8POLY_PRINT prints out a polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 February 2015
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input, double precision A(0:N), the polynomial coefficients.
c    A(0) is the constant term and
c    A(N) is the coefficient of X**N.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(0:n)
      integer i
      double precision mag
      character plus_minus
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      if ( n .le. 0 ) then
        write ( *, '( ''  p(x) = 0'' )' )
        return
      end if

      if ( a(n) .lt. 0.0D+00 ) then
        plus_minus = '-'
      else
        plus_minus = ' '
      end if

      mag = abs ( a(n) )

      if ( 2 .le. n ) then
        write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' )
     &    plus_minus, mag, n
      else if ( n .eq. 1 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' )
     &    plus_minus, mag
      else if ( n .eq. 0 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
      end if

      do i = n - 1, 0, -1

        if ( a(i) .lt. 0.0D+00 ) then
          plus_minus = '-'
        else
          plus_minus = '+'
        end if

        mag = abs ( a(i) )

        if ( mag .ne. 0.0D+00 ) then

          if ( 2 .le. i ) then
            write ( *,
     &        ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' )
     &        plus_minus, mag, i
          else if ( i .eq. 1 ) then
            write ( *,
     &        ' ( ''         '', a1, g14.6, '' * x'' )' )
     &        plus_minus, mag
          else if ( i .eq. 0 ) then
            write ( *, ' ( ''         '', a1, g14.6 )' )
     &        plus_minus, mag
          end if
        end if

      end do

      return
      end
      function r8poly_value_horner ( m, c, x )

c*********************************************************************72
c
cc R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
c
c  Discussion:
c
c    The polynomial 
c
c      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
c
c    is to be evaluated at X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 January 2015
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the degree.
c
c    Input, double precision C(0:M), the polynomial coefficients.  
c    C(I) is the coefficient of X^I.
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision R8POLY_VALUE_HORNER, the polynomial value.
c
      implicit none

      integer m

      double precision c(0:m)
      integer i
      double precision r8poly_value_horner
      double precision value
      double precision x

      value = c(m)
      do i = m - 1, 0, -1
        value = value * x + c(i)
      end do

      r8poly_value_horner = value

      return
      end
      subroutine r8vec_linspace ( n, a, b, x )

c*********************************************************************72
c
cc R8VEC_LINSPACE creates a vector of linearly spaced values.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
c
c    In other words, the interval is divided into N-1 even subintervals,
c    and the endpoints of intervals are used as the points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A, B, the first and last entries.
c
c    Output, double precision X(N), a vector of linearly spaced data.
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      double precision x(n)

      if ( n .eq. 1 ) then

        x(1) = ( a + b ) / 2.0D+00

      else

        do i = 1, n
          x(i) = ( dble ( n - i     ) * a
     &           + dble (     i - 1 ) * b )
     &           / dble ( n     - 1 )
        end do

      end if

      return
      end
      subroutine r8vec_max ( n, a, amax )

c*********************************************************************72
c
cc R8VEC_MAX returns the maximum value in an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, double precision A(N), the array.
c
c    Output, double precision AMAX, the value of the largest entry.
c
      implicit none

      integer n

      double precision a(n)
      double precision amax
      integer i

      amax = a(1)
      do i = 2, n
        amax = max ( amax, a(i) )
      end do

      return
      end
      subroutine r8vec_mean ( n, a, mean )

c*********************************************************************72
c
cc R8VEC_MEAN returns the mean of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A(N), the vector whose mean is desired.
c
c    Output, double precision MEAN, the mean of the vector entries.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision mean

      mean = 0.0D+00
      do i = 1, n
        mean = mean + a(i)
      end do
      mean = mean / dble ( n )

      return
      end
      subroutine r8vec_min ( n, a, amin )

c*********************************************************************72
c
cc R8VEC_MIN returns the minimum value in an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, double precision A(N), the array.
c
c    Output, double precision AMIN, the value of the smallest entry.
c
      implicit none

      integer n

      double precision a(n)
      double precision amin
      integer i

      amin = a(1)
      do i = 2, n
        amin = min ( amin, a(i) )
      end do

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      subroutine r8vec_variance ( n, a, variance )

c*********************************************************************72
c
cc R8VEC_VARIANCE returns the variance of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The variance of a vector X of length N is defined as
c
c      mean ( X(1:n) ) = sum ( X(1:n) ) / n
c
c      var ( X(1:n) ) = sum ( ( X(1:n) - mean )^2 ) / ( n - 1 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c    N should be at least 2.
c
c    Input, double precision A(N), the vector.
c
c    Output, double precision VARIANCE, the variance of the vector.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision mean
      double precision variance

      if ( n .lt. 2 ) then

        variance = 0.0D+00

      else

        mean = 0.0D+00
        do i = 1, n
          mean = mean + a(i)
        end do
        mean = mean / dble ( n )

        variance = 0.0D+00
        do i = 1, n
          variance = variance + ( a(i) - mean ) ** 2
        end do
        variance = variance / dble ( n - 1 )

      end if

      return
      end
      subroutine timestamp_TN ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
      subroutine truncated_normal_ab_cdf ( x, mu, sigma, a, b, cdf )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_AB_CDF evaluates the truncated Normal CDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the CDF.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, B, the lower and upper truncation limits.
c
c    Output, double precision CDF, the value of the CDF.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision b
      double precision beta
      double precision beta_cdf
      double precision cdf
      double precision mu
      double precision sigma
      double precision x
      double precision xi
      double precision xi_cdf

      alpha = ( a - mu ) / sigma
      beta = ( b - mu ) / sigma
      xi = ( x - mu ) / sigma

      call normal_01_cdf ( alpha, alpha_cdf )
      call normal_01_cdf ( beta, beta_cdf )
      call normal_01_cdf ( xi, xi_cdf )

      cdf = ( xi_cdf - alpha_cdf ) / ( beta_cdf - alpha_cdf )

      return
      end
      subroutine truncated_normal_ab_cdf_values ( n_data, mu, sigma, 
     &  a, b, x, fx )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_AB_CDF_VALUES: values of the Truncated Normal CDF.
c
c  Discussion:
c
c    The Normal distribution, with mean Mu and standard deviation Sigma,
c    is truncated to the interval [A,B].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 September 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision MU, the mean of the distribution.
c
c    Output, double precision SIGMA, the standard deviation of the distribution.
c
c    Output, double precision A, B, the lower and upper truncation limits.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      double precision a
      double precision a_vec(n_max)
      double precision b
      double precision b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      double precision mu
      double precision mu_vec(n_max)
      integer n_data
      double precision sigma
      double precision sigma_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save b_vec
      save fx_vec
      save mu_vec
      save sigma_vec
      save x_vec

      data a_vec /
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00 /
      data b_vec /
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00 /
      data fx_vec /
     & 0.3371694242213513D+00,
     & 0.3685009225506048D+00,
     & 0.4006444233448185D+00,
     & 0.4334107066903040D+00,
     & 0.4665988676496338D+00,
     & 0.5000000000000000D+00,
     & 0.5334011323503662D+00,
     & 0.5665892933096960D+00,
     & 0.5993555766551815D+00,
     & 0.6314990774493952D+00,
     & 0.6628305757786487D+00 /
      data mu_vec /
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00 /
      data sigma_vec /
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00 /
      data x_vec /
     &  90.0D+00,
     &  92.0D+00,
     &  94.0D+00,
     &  96.0D+00,
     &  98.0D+00,
     & 100.0D+00,
     & 102.0D+00,
     & 104.0D+00,
     & 106.0D+00,
     & 108.0D+00,
     & 110.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        mu = 0.0D+00
        sigma = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        mu = mu_vec(n_data)
        sigma = sigma_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine truncated_normal_ab_cdf_inv ( cdf, mu, sigma, a, b, x )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_AB_CDF_INV inverts the truncated Normal CDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision CDF, the value of the CDF.
c    0.0D+00 <= CDF <= 1.0.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, B, the lower and upper truncation limits.
c
c    Output, double precision X, the corresponding argument.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision b
      double precision beta
      double precision beta_cdf
      double precision cdf
      double precision mu
      double precision sigma
      double precision x
      double precision xi
      double precision xi_cdf

      if ( cdf .lt. 0.0D+00 .or. 1.0D+00 .lt. cdf ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_CDF_INV - Fatal error!'
        write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
        stop 1
      end if

      alpha = ( a - mu ) / sigma
      beta = ( b - mu ) / sigma

      call normal_01_cdf ( alpha, alpha_cdf )
      call normal_01_cdf ( beta, beta_cdf )

      xi_cdf = ( beta_cdf - alpha_cdf ) * cdf + alpha_cdf
      call normal_01_cdf_inv ( xi_cdf, xi )

      x = mu + sigma * xi

      return
      end
      subroutine truncated_normal_ab_mean ( mu, sigma, a, b, mean )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_AB_MEAN returns the mean of the truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, B, the lower and upper truncation limits.
c
c    Output, double precision MEAN, the mean of the PDF.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision alpha_pdf
      double precision b
      double precision beta
      double precision beta_cdf
      double precision beta_pdf
      double precision mean
      double precision mu
      double precision sigma

      alpha = ( a - mu ) / sigma
      beta = ( b - mu ) / sigma

      call normal_01_cdf ( alpha, alpha_cdf )
      call normal_01_cdf ( beta, beta_cdf )

      call normal_01_pdf ( alpha, alpha_pdf )
      call normal_01_pdf ( beta, beta_pdf )

      mean = mu + sigma * ( alpha_pdf - beta_pdf ) 
     &  / ( beta_cdf - alpha_cdf )

      return
      end
      subroutine truncated_normal_ab_moment ( order, mu, sigma, a, b, 
     &  moment )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_AB_MOMENT: moments of the truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Phoebus Dhrymes,
c    Moments of Truncated Normal Distributions,
c    May 2005.
c
c  Parameters:
c
c    Input, integer ORDER, the order of the moment.
c    0 <= ORDER.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c    0.0 < SIGMA.
c
c    Input, double precision A, B, the lower and upper truncation limits.
c
c    Output, double precision MOMENT, the moment of the PDF.
c
      implicit none

      double precision a
      double precision a_h
      double precision a_cdf
      double precision a_pdf
      double precision b
      double precision b_h
      double precision b_cdf
      double precision b_pdf
      double precision ir
      double precision irm1
      double precision irm2
      double precision moment
      double precision mu
      integer order
      integer r
      double precision r8_choose
      double precision sigma

      if ( order .lt. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_MOMENT - Fatal error!'
        write ( *, '(a)' ) '  ORDER < 0.'
        stop 1
      end if

      if ( sigma .le. 0.0D+00 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_MOMENT - Fatal error!'
        write ( *, '(a)' ) '  SIGMA <= 0.0.'
        stop 1
      end if

      if ( b .le. a ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_MOMENT - Fatal error!'
        write ( *, '(a)' ) '  B <= A.'
        stop 1
      end if

      a_h = ( a - mu ) / sigma
      call normal_01_pdf ( a_h, a_pdf )
      call normal_01_cdf ( a_h, a_cdf )

      if ( a_cdf .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_MOMENT - Fatal error!'
        write ( *, '(a)' ) 
     &    '  PDF/CDF ratio fails, because A_CDF is too small.'
        write ( *, '(a,g14.6)' ) '  A_PDF = %g\n', a_pdf
        write ( *, '(a,g14.6)' ) '  A_CDF = %g\n', a_cdf
        stop 1
      end if

      b_h = ( b - mu ) / sigma
      call normal_01_pdf ( b_h, b_pdf )
      call normal_01_cdf ( b_h, b_cdf )

      if ( b_cdf .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_MOMENT - Fatal error!'
        write ( *, '(a)' ) 
     &    '  PDF/CDF ratio fails, because B_CDF is too small.'
        write ( *, '(a,g14.6)' ) '  B_PDF = %g\n', b_pdf
        write ( *, '(a,g14.6)' ) '  B_CDF = %g\n', b_cdf
        stop 1
      end if

      moment = 0.0D+00
      irm2 = 0.0D+00
      irm1 = 0.0D+00

      do r = 0, order

        if ( r .eq. 0 ) then
          ir = 1.0D+00
        else if ( r .eq. 1 ) then
          ir = - ( b_pdf - a_pdf ) / ( b_cdf - a_cdf )
        else
          ir = dble ( r - 1 ) * irm2 
     &      - ( b_h ** ( r - 1 ) * b_pdf - a_h ** ( r - 1 ) * a_pdf ) 
     &      / ( b_cdf - a_cdf )
        end if

        moment = moment + r8_choose ( order, r ) * mu ** ( order - r ) 
     &    * ( sigma ** r ) * ir

        irm2 = irm1
        irm1 = ir

      end do

      return
      end
      subroutine truncated_normal_ab_pdf ( x, mu, sigma, a, b, pdf )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_AB_PDF evaluates the truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the PDF.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, B, the lower and upper truncation limits.
c
c    Output, double precision PDF, the value of the PDF.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision b
      double precision beta
      double precision beta_cdf
      double precision mu
      double precision pdf
      double precision sigma
      double precision x
      double precision xi
      double precision xi_pdf

      alpha = ( a - mu ) / sigma
      beta = ( b - mu ) / sigma
      xi = ( x - mu ) / sigma

      call normal_01_cdf ( alpha, alpha_cdf )
      call normal_01_cdf ( beta, beta_cdf )
      call normal_01_pdf ( xi, xi_pdf )

      pdf = xi_pdf / ( beta_cdf - alpha_cdf ) / sigma

      return
      end
      subroutine truncated_normal_ab_pdf_values ( n_data, mu, sigma, 
     &  a, b, x, fx )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_AB_PDF_VALUES: values of the Truncated Normal PDF.
c
c  Discussion:
c
c    The Normal distribution, with mean Mu and standard deviation Sigma,
c    is truncated to the interval [A,B].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 September 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision MU, the mean of the distribution.
c
c    Output, double precision SIGMA, the variance of the distribution.
c
c    Output, double precision A, B, the lower and upper truncation limits.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      double precision a
      double precision a_vec(n_max)
      double precision b
      double precision b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      double precision mu
      double precision mu_vec(n_max)
      integer n_data
      double precision sigma
      double precision sigma_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save b_vec
      save fx_vec
      save mu_vec
      save sigma_vec
      save x_vec

      data a_vec /
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00 /
      data b_vec /
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00 /
      data fx_vec /
     &  0.01543301171801836D+00,
     &  0.01588394472270638D+00,
     &  0.01624375997031919D+00,
     &  0.01650575046469259D+00,
     &  0.01666496869385951D+00,
     &  0.01671838200940538D+00,
     &  0.01666496869385951D+00,
     &  0.01650575046469259D+00,
     &  0.01624375997031919D+00,
     &  0.01588394472270638D+00,
     &  0.01543301171801836D+00 /
      data mu_vec /
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00 /
      data sigma_vec /
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00 /
      data x_vec /
     &  90.0D+00,
     &  92.0D+00,
     &  94.0D+00,
     &  96.0D+00,
     &  98.0D+00,
     & 100.0D+00,
     & 102.0D+00,
     & 104.0D+00,
     & 106.0D+00,
     & 108.0D+00,
     & 110.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        mu = 0.0D+00
        sigma = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        mu = mu_vec(n_data)
        sigma = sigma_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine truncated_normal_ab_sample ( mu, sigma, a, b, seed, x )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_AB_SAMPLE samples the truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, B, the lower and upper truncation limits.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X, a sample of the PDF.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision b
      double precision beta
      double precision beta_cdf
      double precision mu
      double precision r8_uniform_01
      integer seed
      double precision sigma
      double precision u
      double precision x
      double precision xi
      double precision xi_cdf

      alpha = ( a - mu ) / sigma
      beta = ( b - mu ) / sigma

      call normal_01_cdf ( alpha, alpha_cdf )
      call normal_01_cdf ( beta, beta_cdf )

      u = r8_uniform_01 ( seed )
      xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf )
      call normal_01_cdf_inv ( xi_cdf, xi )

      x = mu + sigma * xi

      return
      end
      subroutine truncated_normal_ab_variance ( mu, sigma, a, b, 
     &  variance )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_AB_VARIANCE returns the variance of the truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, B, the lower and upper truncation limits.
c
c    Output, double precision VARIANCE, the variance of the PDF.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision alpha_pdf
      double precision b
      double precision beta
      double precision beta_cdf
      double precision beta_pdf
      double precision mu
      double precision sigma
      double precision variance

      alpha = ( a - mu ) / sigma
      beta = ( b - mu ) / sigma

      call normal_01_pdf ( alpha, alpha_pdf )
      call normal_01_pdf ( beta, beta_pdf )

      call normal_01_cdf ( alpha, alpha_cdf )
      call normal_01_cdf ( beta, beta_cdf )

      variance = sigma * sigma * ( 1.0D+00 
     &  + ( alpha * alpha_pdf - beta * beta_pdf ) 
     &  / ( beta_cdf - alpha_cdf ) 
     &  - ( ( alpha_pdf - beta_pdf ) / ( beta_cdf - alpha_cdf ) ) ** 2 )

      return
      end
      subroutine truncated_normal_a_cdf ( x, mu, sigma, a, cdf )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_A_CDF evaluates the lower truncated Normal CDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the CDF.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, the lower truncation limit.
c
c    Output, double precision CDF, the value of the CDF.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision cdf
      double precision mu
      double precision sigma
      double precision x
      double precision xi
      double precision xi_cdf

      alpha = ( a - mu ) / sigma
      xi = ( x - mu ) / sigma

      call normal_01_cdf ( alpha, alpha_cdf )
      call normal_01_cdf ( xi, xi_cdf )

      cdf = ( xi_cdf - alpha_cdf ) / ( 1.0D+00 - alpha_cdf )

      return
      end
      subroutine truncated_normal_a_cdf_values ( n_data, mu, sigma, 
     &  a, x, fx )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_A_CDF_VALUES: values of the lower Truncated Normal CDF.
c
c  Discussion:
c
c    The Normal distribution, with mean Mu and standard deviation Sigma,
c    is truncated to the interval [A,B].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 September 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision MU, the mean of the distribution.
c
c    Output, double precision SIGMA, the standard deviation of the distribution.
c
c    Output, double precision A, the lower truncation limit.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      double precision a
      double precision a_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      double precision mu
      double precision mu_vec(n_max)
      integer n_data
      double precision sigma
      double precision sigma_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save fx_vec
      save mu_vec
      save sigma_vec
      save x_vec

      data a_vec /
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00 /
      data fx_vec /
     &  0.3293202045481688D+00, 
     &  0.3599223134505957D+00, 
     &  0.3913175216041539D+00, 
     &  0.4233210140873113D+00, 
     &  0.4557365629792204D+00, 
     &  0.4883601253415709D+00, 
     &  0.5209836877039214D+00, 
     &  0.5533992365958304D+00, 
     &  0.5854027290789878D+00, 
     &  0.6167979372325460D+00, 
     &  0.6474000461349729D+00 /
      data mu_vec /
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00 /
      data sigma_vec /
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00 /
      data x_vec /
     &  90.0D+00,
     &  92.0D+00,
     &  94.0D+00,
     &  96.0D+00,
     &  98.0D+00,
     & 100.0D+00,
     & 102.0D+00,
     & 104.0D+00,
     & 106.0D+00,
     & 108.0D+00,
     & 110.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        mu = 0.0D+00
        sigma = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        mu = mu_vec(n_data)
        sigma = sigma_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine truncated_normal_a_cdf_inv ( cdf, mu, sigma, a, x )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_A_CDF_INV inverts the lower truncated Normal CDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision CDF, the value of the CDF.
c    0.0D+00 <= CDF <= 1.0.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, the lower truncation limit.
c
c    Output, double precision X, the corresponding argument.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision cdf
      double precision mu
      double precision sigma
      double precision x
      double precision xi
      double precision xi_cdf

      if ( cdf .lt. 0.0D+00 .or. 1.0D+00 .lt. cdf ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRUNCATED_NORMAL_A_CDF_INV - Fatal error!'
        write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
        stop 1
      end if

      alpha = ( a - mu ) / sigma

      call normal_01_cdf ( alpha, alpha_cdf )

      xi_cdf = ( 1.0D+00 - alpha_cdf ) * cdf + alpha_cdf
      call normal_01_cdf_inv ( xi_cdf, xi )

      x = mu + sigma * xi

      return
      end
      subroutine truncated_normal_a_mean ( mu, sigma, a, mean )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_A_MEAN returns the mean of the lower truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, the lower truncation limit.
c
c    Output, double precision MEAN, the mean of the PDF.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision alpha_pdf
      double precision mean
      double precision mu
      double precision sigma

      alpha = ( a - mu ) / sigma

      call normal_01_cdf ( alpha, alpha_cdf )

      call normal_01_pdf ( alpha, alpha_pdf )

      mean = mu + sigma * alpha_pdf
     &  / ( 1.0D+00 - alpha_cdf )

      return
      end
      subroutine truncated_normal_a_moment ( order, mu, sigma, a, 
     &  moment )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_A_MOMENT: moments of the lower truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Phoebus Dhrymes,
c    Moments of Truncated Normal Distributions,
c    May 2005.
c
c  Parameters:
c
c    Input, integer ORDER, the order of the moment.
c    0 <= ORDER.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c    0.0 < SIGMA.
c
c    Input, double precision A, the lower truncation limit.
c
c    Output, double precision MOMENT, the moment of the PDF.
c
      implicit none

      double precision a
      double precision moment
      double precision mu
      integer order
      double precision r8_mop
      double precision sigma

      call truncated_normal_b_moment ( order, - mu, sigma, - a, moment )

      moment = r8_mop ( order ) * moment

      return
      end
      subroutine truncated_normal_a_pdf ( x, mu, sigma, a, pdf )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_A_PDF evaluates the lower truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the PDF.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, the lower truncation limit.
c
c    Output, double precision PDF, the value of the PDF.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision mu
      double precision pdf
      double precision sigma
      double precision x
      double precision xi
      double precision xi_pdf

      alpha = ( a - mu ) / sigma
      xi = ( x - mu ) / sigma

      call normal_01_cdf ( alpha, alpha_cdf )
      call normal_01_pdf ( xi, xi_pdf )

      pdf = xi_pdf / ( 1.0D+00 - alpha_cdf ) / sigma

      return
      end
      subroutine truncated_normal_a_pdf_values ( n_data, mu, sigma, 
     &  a, x, fx )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_A_PDF_VALUES: values of the lower Truncated Normal PDF.
c
c  Discussion:
c
c    The Normal distribution, with mean Mu and standard deviation Sigma,
c    is truncated to the interval [A,+oo).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 September 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision MU, the mean of the distribution.
c
c    Output, double precision SIGMA, the standard deviation of the distribution.
c
c    Output, double precision A, the lower truncation limit.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      double precision a
      double precision a_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      double precision mu
      double precision mu_vec(n_max)
      integer n_data
      double precision sigma
      double precision sigma_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save fx_vec
      save mu_vec
      save sigma_vec
      save x_vec

      data a_vec /
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00,
     &  50.0D+00 /
      data fx_vec /
     &  0.01507373507401876D+00,
     &  0.01551417047139894D+00,
     &  0.01586560931024694D+00,
     &  0.01612150073158793D+00,
     &  0.01627701240029317D+00,
     &  0.01632918226724295D+00,
     &  0.01627701240029317D+00,
     &  0.01612150073158793D+00,
     &  0.01586560931024694D+00,
     &  0.01551417047139894D+00,
     &  0.01507373507401876D+00 /
      data mu_vec /
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00 /
      data sigma_vec /
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00 /
      data x_vec /
     &  90.0D+00,
     &  92.0D+00,
     &  94.0D+00,
     &  96.0D+00,
     &  98.0D+00,
     & 100.0D+00,
     & 102.0D+00,
     & 104.0D+00,
     & 106.0D+00,
     & 108.0D+00,
     & 110.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        mu = 0.0D+00
        sigma = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        mu = mu_vec(n_data)
        sigma = sigma_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine truncated_normal_a_sample ( mu, sigma, a, seed, x )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_A_SAMPLE samples the lower truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, the lower truncation limit.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X, a sample of the PDF.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision mu
      double precision r8_uniform_01
      integer seed
      double precision sigma
      double precision u
      double precision x
      double precision xi
      double precision xi_cdf

      alpha = ( a - mu ) / sigma

      call normal_01_cdf ( alpha, alpha_cdf )

      u = r8_uniform_01 ( seed )
      xi_cdf = alpha_cdf + u * ( 1.0D+00 - alpha_cdf )
      call normal_01_cdf_inv ( xi_cdf, xi )

      x = mu + sigma * xi

      return
      end
      subroutine truncated_normal_a_variance ( mu, sigma, a, variance )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_A_VARIANCE: variance of the lower truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision A, the lower truncation limit.
c
c    Output, double precision VARIANCE, the variance of the PDF.
c
      implicit none

      double precision a
      double precision alpha
      double precision alpha_cdf
      double precision alpha_pdf
      double precision mu
      double precision sigma
      double precision variance

      alpha = ( a - mu ) / sigma

      call normal_01_pdf ( alpha, alpha_pdf )

      call normal_01_cdf ( alpha, alpha_cdf )

      variance = sigma * sigma * ( 1.0D+00 
     &  + ( alpha * alpha_pdf ) 
     &  / ( 1.0D+00 - alpha_cdf ) 
     &  - ( alpha_pdf / ( 1.0D+00 - alpha_cdf ) ) ** 2 )

      return
      end
      subroutine truncated_normal_b_cdf_values ( n_data, mu, sigma, 
     &  b, x, fx )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_B_CDF_VALUES: values of the Upper Truncated Normal CDF.
c
c  Discussion:
c
c    The Normal distribution, with mean Mu and standard deviation Sigma,
c    is truncated to the interval (-oo,B].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 September 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision MU, the mean of the distribution.
c
c    Output, double precision SIGMA, the standard deviation of the distribution.
c
c    Output, double precision B, the upper truncation limit.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      double precision b
      double precision b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      double precision mu
      double precision mu_vec(n_max)
      integer n_data
      double precision sigma
      double precision sigma_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save b_vec
      save fx_vec
      save mu_vec
      save sigma_vec
      save x_vec

      data b_vec /
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00 /
      data fx_vec /
     &  0.3525999538650271D+00, 
     &  0.3832020627674540D+00, 
     &  0.4145972709210122D+00, 
     &  0.4466007634041696D+00, 
     &  0.4790163122960786D+00, 
     &  0.5116398746584291D+00, 
     &  0.5442634370207796D+00, 
     &  0.5766789859126887D+00,
     &  0.6086824783958461D+00, 
     &  0.6400776865494043D+00, 
     &  0.6706797954518312D+00 /
      data mu_vec /
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00 /
      data sigma_vec /
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00 /
      data x_vec /
     &  90.0D+00,
     &  92.0D+00,
     &  94.0D+00,
     &  96.0D+00,
     &  98.0D+00,
     & 100.0D+00,
     & 102.0D+00,
     & 104.0D+00,
     & 106.0D+00,
     & 108.0D+00,
     & 110.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        b = 0.0D+00
        mu = 0.0D+00
        sigma = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        b = b_vec(n_data)
        mu = mu_vec(n_data)
        sigma = sigma_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine truncated_normal_b_cdf ( x, mu, sigma, b, cdf )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_B_CDF evaluates the upper truncated Normal CDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the CDF.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision B, the upper truncation limit.
c
c    Output, double precision CDF, the value of the CDF.
c
      implicit none

      double precision b
      double precision beta
      double precision beta_cdf
      double precision cdf
      double precision mu
      double precision sigma
      double precision x
      double precision xi
      double precision xi_cdf

      beta = ( b - mu ) / sigma
      xi = ( x - mu ) / sigma

      call normal_01_cdf ( beta, beta_cdf )
      call normal_01_cdf ( xi, xi_cdf )

      cdf = xi_cdf / beta_cdf

      return
      end
      subroutine truncated_normal_b_cdf_inv ( cdf, mu, sigma, b, x )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_B_CDF_INV inverts the upper truncated Normal CDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision CDF, the value of the CDF.
c    0.0D+00 <= CDF <= 1.0.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision B, the upper truncation limit.
c
c    Output, double precision X, the corresponding argument.
c
      implicit none

      double precision b
      double precision beta
      double precision beta_cdf
      double precision cdf
      double precision mu
      double precision sigma
      double precision x
      double precision xi
      double precision xi_cdf

      if ( cdf .lt. 0.0D+00 .or. 1.0D+00 .lt. cdf ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRUNCATED_NORMAL_B_CDF_INV - Fatal error!'
        write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
        stop 1
      end if

      beta = ( b - mu ) / sigma

      call normal_01_cdf ( beta, beta_cdf )

      xi_cdf = beta_cdf * cdf
      call normal_01_cdf_inv ( xi_cdf, xi )

      x = mu + sigma * xi

      return
      end
      subroutine truncated_normal_b_mean ( mu, sigma, b, mean )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_B_MEAN returns the mean of the upper truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision B, the upper truncation limit.
c
c    Output, double precision MEAN, the mean of the PDF.
c
      implicit none

      double precision b
      double precision beta
      double precision beta_cdf
      double precision beta_pdf
      double precision mean
      double precision mu
      double precision sigma

      beta = ( b - mu ) / sigma

      call normal_01_cdf ( beta, beta_cdf )

      call normal_01_pdf ( beta, beta_pdf )

      mean = mu - sigma * beta_pdf / beta_cdf

      return
      end
      subroutine truncated_normal_b_moment ( order, mu, sigma, b, 
     &  moment )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_B_MOMENT: moments of the upper truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Phoebus Dhrymes,
c    Moments of Truncated Normal Distributions,
c    May 2005.
c
c  Parameters:
c
c    Input, integer ORDER, the order of the moment.
c    0 <= ORDER.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c    0.0 < SIGMA.
c
c    Input, double precision B, the upper truncation limit.
c
c    Output, double precision MOMENT, the moment of the PDF.
c
      implicit none

      double precision b
      double precision f
      double precision h
      double precision h_cdf
      double precision h_pdf
      double precision ir
      double precision irm1
      double precision irm2
      double precision moment
      double precision mu
      integer order
      integer r
      double precision r8_choose
      double precision sigma

      if ( order .lt. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'TRUNCATED_NORMAL_B_MOMENT - Fatal error!'
        write ( *, '(a)' ) '  ORDER < 0.'
        stop 1
      end if

      h = ( b - mu ) / sigma
      call normal_01_pdf ( h, h_pdf )
      call normal_01_cdf ( h, h_cdf )

      if ( h_cdf .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'TRUNCATED_NORMAL_B_MOMENT - Fatal error!'
        write ( *, '(a)' ) '  CDF((B-MU)/SIGMA) = 0.'
        stop 1
      end if

      f = h_pdf / h_cdf

      moment = 0.0D+00
      irm2 = 0.0D+00
      irm1 = 0.0D+00

      do r = 0, order

        if ( r .eq. 0 ) then
          ir = 1.0D+00
        else if ( r .eq. 1 ) then
          ir = - f
        else
          ir = - h ** ( r - 1 ) * f + dble ( r - 1 ) * irm2
        end if

        moment = moment + r8_choose ( order, r ) * mu ** ( order - r ) 
     &    * ( sigma ** r ) * ir

        irm2 = irm1
        irm1 = ir

      end do

      return
      end
      subroutine truncated_normal_b_pdf ( x, mu, sigma, b, pdf )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_B_PDF evaluates the upper truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the PDF.
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision B, the upper truncation limit.
c
c    Output, double precision PDF, the value of the PDF.
c
      implicit none

      double precision b
      double precision beta
      double precision beta_cdf
      double precision mu
      double precision pdf
      double precision sigma
      double precision x
      double precision xi
      double precision xi_pdf

      beta = ( b - mu ) / sigma
      xi = ( x - mu ) / sigma

      call normal_01_cdf ( beta, beta_cdf )
      call normal_01_pdf ( xi, xi_pdf )

      pdf = xi_pdf / beta_cdf / sigma

      return
      end
      subroutine truncated_normal_b_pdf_values ( n_data, mu, sigma, 
     &  b, x, fx )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_B_PDF_VALUES: values of the upper Truncated Normal PDF.
c
c  Discussion:
c
c    The Normal distribution, with mean Mu and standard deviation Sigma,
c    is truncated to the interval (-oo,B].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 September 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision MU, the mean of the distribution.
c
c    Output, double precision SIGMA, the standard deviation of the distribution.
c
c    Output, double precision B, the upper truncation limit.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      double precision b
      double precision b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      double precision mu
      double precision mu_vec(n_max)
      integer n_data
      double precision sigma
      double precision sigma_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save b_vec
      save fx_vec
      save mu_vec
      save sigma_vec
      save x_vec

      data b_vec /
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00,
     &  150.0D+00 /
      data fx_vec /
     &  0.01507373507401876D+00, 
     &  0.01551417047139894D+00, 
     &  0.01586560931024694D+00, 
     &  0.01612150073158793D+00, 
     &  0.01627701240029317D+00, 
     &  0.01632918226724295D+00, 
     &  0.01627701240029317D+00, 
     &  0.01612150073158793D+00, 
     &  0.01586560931024694D+00, 
     &  0.01551417047139894D+00, 
     &  0.01507373507401876D+00 /
      data mu_vec /
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00,
     &  100.0D+00 /
      data sigma_vec /
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00,
     &  25.0D+00 /
      data x_vec /
     &  90.0D+00,
     &  92.0D+00,
     &  94.0D+00,
     &  96.0D+00,
     &  98.0D+00,
     & 100.0D+00,
     & 102.0D+00,
     & 104.0D+00,
     & 106.0D+00,
     & 108.0D+00,
     & 110.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        b = 0.0D+00
        mu = 0.0D+00
        sigma = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        b = b_vec(n_data)
        mu = mu_vec(n_data)
        sigma = sigma_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine truncated_normal_b_sample ( mu, sigma, b, seed, x )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_B_SAMPLE samples the upper truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision B, the upper truncation limit.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X, a sample of the PDF.
c
      implicit none

      double precision b
      double precision beta
      double precision beta_cdf
      double precision mu
      double precision r8_uniform_01
      integer seed
      double precision sigma
      double precision u
      double precision x
      double precision xi
      double precision xi_cdf

      beta = ( b - mu ) / sigma

      call normal_01_cdf ( beta, beta_cdf )

      u = r8_uniform_01 ( seed )
      xi_cdf = u * beta_cdf
      call normal_01_cdf_inv ( xi_cdf, xi )

      x = mu + sigma * xi

      return
      end
      subroutine truncated_normal_b_variance ( mu, sigma, b, variance )

c*********************************************************************72
c
cc TRUNCATED_NORMAL_B_VARIANCE: variance of the upper truncated Normal PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision MU, SIGMA, the mean and standard deviation of the
c    parent Normal distribution.
c
c    Input, double precision B, the upper truncation limit.
c
c    Output, double precision VARIANCE, the variance of the PDF.
c
      implicit none

      double precision b
      double precision beta
      double precision beta_cdf
      double precision beta_pdf
      double precision mu
      double precision sigma
      double precision variance

      beta = ( b - mu ) / sigma

      call normal_01_pdf ( beta, beta_pdf )

      call normal_01_cdf ( beta, beta_cdf )

      variance = sigma * sigma * ( 1.0D+00 
     &  - beta * beta_pdf / beta_cdf 
     &  - ( beta_pdf / beta_cdf ) ** 2 )

      return
      end

