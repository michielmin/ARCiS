subroutine timescale(nr, T, nt, nmol, molname, n_eq, COratio, tau)
  implicit none
  integer(4), intent(in) :: nr, nmol
  real(8), intent(in) :: T(nr), nt(nr), n_eq(nmol, nr), COratio
  character(10), intent(in) :: molname(nmol)
  real(8), intent(out) :: tau(nmol, nr)

  real(8) :: k(18, nr)
  integer(4) :: j

  
  call rate(nr, T(1:nr), nt(1:nr), k(1:18, 1:nr))

  
  !Tsai et al. (2018)
  !H2
  tau(45, 1:nr) = n_eq(45, 1:nr) / (k(17, 1:nr) * n_eq(64, 1:nr) * n_eq(64, 1:nr))

  if(COratio <= 1.0d0) then
     !CH4
     tau(6, 1:nr) = n_eq(6, 1:nr) / (k(1, 1:nr) * n_eq(63, 1:nr) * n_eq(13, 1:nr) &
          + max(min(k(2, 1:nr) * n_eq(62, 1:nr) * n_eq(64, 1:nr), k(3, 1:nr) * n_eq(62, 1:nr)), &
          k(4, 1:nr) * n_eq(63, 1:nr) * n_eq(34, 1:nr) + k(9, 1:nr) * n_eq(39, 1:nr) * n_eq(64, 1:nr) &
          + min(k(5, 1:nr) * n_eq(61, 1:nr) * n_eq(13, 1:nr), k(6, 1:nr) * n_eq(6, 1:nr) * n_eq(64, 1:nr)))) &
          + tau(45, 1:nr) * 3.0d0 * n_eq(5, 1:nr) / n_eq(45, 1:nr)

    !CO
     tau(5, 1:nr) = n_eq(5, 1:nr) / (k(1, 1:nr) * n_eq(63, 1:nr) * n_eq(13, 1:nr) &
          + max(min(k(2, 1:nr) * n_eq(62, 1:nr) * n_eq(64, 1:nr), k(3, 1:nr) * n_eq(62, 1:nr)), &
          k(4, 1:nr) * n_eq(63, 1:nr) * n_eq(34, 1:nr) + k(9, 1:nr) * n_eq(39, 1:nr) * n_eq(64, 1:nr) &
          + min(k(5, 1:nr) * n_eq(61, 1:nr) * n_eq(13, 1:nr), k(6, 1:nr) * n_eq(6, 1:nr) * n_eq(64, 1:nr)))) &
          + tau(45, 1:nr) * 3.0d0 * n_eq(5, 1:nr) / n_eq(45, 1:nr)

    !H2O
     tau(1, 1:nr) = n_eq(1, 1:nr) / (k(1, 1:nr) * n_eq(63, 1:nr) * n_eq(13, 1:nr) &
          + max(min(k(2, 1:nr) * n_eq(62, 1:nr) * n_eq(64, 1:nr), k(3, 1:nr) * n_eq(62, 1:nr)), &
          k(4, 1:nr) * n_eq(63, 1:nr) * n_eq(34, 1:nr) + k(9, 1:nr) * n_eq(39, 1:nr) * n_eq(64, 1:nr) &
          + min(k(5, 1:nr) * n_eq(61, 1:nr) * n_eq(13, 1:nr), k(6, 1:nr) * n_eq(6, 1:nr) * n_eq(64, 1:nr)))) &
          + tau(45, 1:nr) * 3.0d0 * n_eq(5, 1:nr) / n_eq(45, 1:nr)
     
  else
     !CH4
     tau(6, 1:nr) = n_eq(6, 1:nr) / (k(1, 1:nr) * n_eq(63, 1:nr) * n_eq(13, 1:nr) &
          + min(k(2, 1:nr) * n_eq(62, 1:nr) * n_eq(64, 1:nr), k(3, 1:nr) * n_eq(62, 1:nr)) &
          + k(9, 1:nr) * n_eq(39, 1:nr) * n_eq(64, 1:nr) &
          + max(k(8, 1:nr) * n_eq(26, 1:nr) * n_eq(34, 1:nr), k(10, 1:nr) * n_eq(26, 1:nr) * n_eq(13, 1:nr)))

     !CO
     tau(5, 1:nr) = n_eq(5, 1:nr) / (k(1, 1:nr) * n_eq(63, 1:nr) * n_eq(13, 1:nr) &
          + min(k(2, 1:nr) * n_eq(62, 1:nr) * n_eq(64, 1:nr), k(3, 1:nr) * n_eq(62, 1:nr)) &
          + k(9, 1:nr) * n_eq(39, 1:nr) * n_eq(64, 1:nr) &
          + max(k(8, 1:nr) * n_eq(26, 1:nr) * n_eq(34, 1:nr), k(10, 1:nr) * n_eq(26, 1:nr) * n_eq(13, 1:nr)))

     !H2O
     tau(1, 1:nr) = n_eq(1, 1:nr) / (k(1, 1:nr) * n_eq(63, 1:nr) * n_eq(13, 1:nr) &
          + min(k(2, 1:nr) * n_eq(62, 1:nr) * n_eq(64, 1:nr), k(3, 1:nr) * n_eq(62, 1:nr)) &
          + k(9, 1:nr) * n_eq(39, 1:nr) * n_eq(64, 1:nr) &
          + max(k(8, 1:nr) * n_eq(26, 1:nr) * n_eq(34, 1:nr), k(10, 1:nr) * n_eq(26, 1:nr) * n_eq(13, 1:nr)))
  end if
  
  !NH3
  tau(11, 1:nr) = 0.5d0 * (n_eq(11, 1:nr) &
       /(max(k(11, 1:nr) * n_eq(67, 1:nr) * n_eq(67, 1:nr), k(12, 1:nr) * n_eq(68, 1:nr)) &
       + k(13, 1:nr) * n_eq(66, 1:nr) * n_eq(67, 1:nr) + k(14, 1:nr) * n_eq(8, 1:nr) * n_eq(67, 1:nr) &
       + k(15, 1:nr) * n_eq(65, 1:nr) * n_eq(8, 1:nr)) &
       + tau(45, 1:nr) * 3.0d0 * n_eq(22, 1:nr) / n_eq(45, 1:nr))

  !N2
  tau(22, 1:nr) = n_eq(22, 1:nr) &
       /(max(k(11, 1:nr) * n_eq(67, 1:nr) * n_eq(67, 1:nr), k(12, 1:nr) * n_eq(68, 1:nr)) &
       + k(13, 1:nr) * n_eq(66, 1:nr) * n_eq(67, 1:nr) + k(14, 1:nr) * n_eq(8, 1:nr) * n_eq(67, 1:nr) &
       + k(15, 1:nr) * n_eq(65, 1:nr) * n_eq(8, 1:nr)) &
       + tau(45, 1:nr) * 3.0d0 * n_eq(22, 1:nr) / n_eq(45, 1:nr)

end subroutine timescale
