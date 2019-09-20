subroutine rate(nr, T, nt, k)
  implicit none
  integer(4), intent(in) :: nr
  real(8), intent(in) :: T(nr), nt(nr)
  real(8), intent(out) :: k(18, nr)
  real(8) :: k_0(nr), k_inf(nr)

  
  !Tsai et al. (2018)
  k_0(1:nr) = 1.93d3 * T(1:nr)**(-9.88d0) * exp(-7544d0/T(1:nr)) &
       + 5.11d-11 * T(1:nr)**(-6.25d0) * exp(-1433d0/T(1:nr))
  k_inf(1:nr) = 1.03d-10 * T(1:nr)**(-0.018d0) * exp(16.74d0/T(1:nr))
  k(1, 1:nr) = (k_0(1:nr) * nt(1:nr)) / (1.0d0 + k_0(1:nr) * nt(1:nr) / k_inf(1:nr))
  
  k(2, 1:nr) = 1.60d-10

  k_0(1:nr) = 1.66d-10 * exp(-12630d0/T(1:nr))
  k_inf(1:nr) = 3d9 * exp(-14600d0/T(1:nr))
  k(3, 1:nr) = (k_0(1:nr) * nt(1:nr)) / (1.0d0 + k_0(1:nr) * nt(1:nr) / k_inf(1:nr))
  

  k(4, 1:nr) = 1.4d-10
  k(5, 1:nr) = 1.05d-12 * T(1:nr)**0.5d0
  k(6, 1:nr) = 2.20d-20 * T(1:nr)**3.0d0 * exp(-4040d0/T(1:nr))  
  k(7, 1:nr) = 6.82d-20 * T(1:nr)**2.685d0 * exp(-4643d0/T(1:nr))  
  k(8, 1:nr) = 6.78d-16 * T(1:nr)**1.5d0 * exp(-854d0/T(1:nr))  
  k(9, 1:nr) = 4.91d-19 * T(1:nr)**2.485d0 * exp(-10380d0/T(1:nr))  
  k(10, 1:nr) = 8.04d-28 * T(1:nr)**4.0d0 * exp(1010d0/T(1:nr))  
  k(11, 1:nr) = 2.89d-16 * T(1:nr)**1.02d0 * exp(-5930d0/T(1:nr))  
 
  k_0(1:nr) = 3.49d38 * T(1:nr)**(-13.13d0) * exp(-36825d0/T(1:nr))
  k_inf(1:nr) = 7.95d13 * exp(-27463d0/T(1:nr))
  k(12, 1:nr) = (k_0(1:nr) * nt(1:nr)) / (1.0d0 + k_0(1:nr) * nt(1:nr) / k_inf(1:nr))
  
  k(13, 1:nr) = 6.98d-10 * T(1:nr)**(-0.27d0) * exp(39d0/T(1:nr))  
  k(14, 1:nr) = 7.9d-9 * T(1:nr)**(-1.1d0) * exp(-98d0/T(1:nr))  
  k(15, 1:nr) = 3.7d-11
  k(16, 1:nr) = 1.17d-11 * exp(-1260d0/T(1:nr))

  k_0(1:nr) = 2.7d-31 * T(1:nr)**(-0.6d0)
  k_inf(1:nr) = 3.31d-6 * T(1:nr)**(-1.0d0)
  k(17, 1:nr) = (k_0(1:nr) * nt(1:nr)) / (1.0d0 + k_0(1:nr) * nt(1:nr) / k_inf(1:nr))
  
  k(18, 1:nr) = 1.05d-17 * T(1:nr)**(1.5d0) * exp(259d0/T(1:nr))  
  
end subroutine rate
