subroutine diseq_calc(nr, R, P, T, nmol, molname, mixrat_r, COratio, Kzz)
  implicit none
  integer(4), intent(in) :: nr, nmol
  real(8), intent(in) :: R(nr+1), P(nr), T(nr), COratio, Kzz
  character(10), intent(in) :: molname(nmol)
  real(8), intent(inout) :: mixrat_r(nr, nmol)
  
  real(8) :: nt(nr), n_eq(nmol, nr), n(nmol, nr), tau(nmol, nr)
  real(8) :: tau_chem(nr, nr), eddy(nr, nr), e(nr, nr), a(nr, nr)

  integer(4) :: i, j

  integer(4) :: LDA, IPIV(nr), INFO, LWORK
  real(8) :: WORK(nr)

  !Boltzmann constant [erg/K]
  real(8), parameter :: kb = 1.3806504d-16


  do j = 1, nr
     nt(j) = P(j) / kb / T(j) * 1.0d6
     n_eq(1:nmol, j) = mixrat_r(j, 1:nmol) * nt(j)
  end do

  call timescale(nr, T(1:nr), nt(1:nr), nmol, molname(1:nmol), n_eq(1:nmol, 1:nr), COratio, tau(1:nmol, 1:nr))
  
  call diffusion(nr, R(1:nr), nt(1:nr), Kzz, eddy(1:nr, 1:nr))

  
  e(1:nr, 1:nr) = 0.0d0
  do j = 1, nr
     e(j, j) = 1.0d0
  end do

  
  n(1:nmol, 1:nr) = n_eq(1:nmol, 1:nr)
  do i = 1, nmol
     if(i == 6 .or. i == 5 .or. i == 1 .or. i == 11 .or. i == 22) then
        tau_chem(1:nr, 1:nr) = 0.0d0
        do j = 1, nr
           tau_chem(j, j) = tau(i, j)
        end do
        a(1:nr, 1:nr) = e(1:nr, 1:nr) + matmul(tau_chem(1:nr, 1:nr), eddy(1:nr, 1:nr))

        LDA = nr
        call DGETRF(nr, nr, a(1:nr, 1:nr), LDA, IPIV, INFO)

        LWORK = nr
        call DGETRI(nr, a(1:nr, 1:nr), LDA, IPIV, WORK, LWORK, INFO)

        n(i, 1:nr) = matmul(a(1:nr, 1:nr), n_eq(i, 1:nr))

        do j = 1, nr
           mixrat_r(j, i) = n(i, j) / nt(j)
        end do
     end if
  end do

  !CO2 calculated with pseudo-eq value
  n(2, 1:nr) = n(5, 1:nr) * n(1, 1:nr) * n_eq(64, 1:nr) / n_eq(5, 1:nr) / n_eq(1, 1:nr) / n(64, 1:nr) * n_eq(2, 1:nr)

end subroutine diseq_calc
