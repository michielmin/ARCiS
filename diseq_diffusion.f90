subroutine diffusion(nr, R, nt, Kzz, eddy)
  implicit none
  integer(4), intent(in) :: nr
  real(8), intent(in) :: R(nr), nt(nr), Kzz
  real(8), intent(out) :: eddy(nr, nr)

  integer(4) :: j
  real(8) :: dr
  
  
  eddy(1:nr, 1:nr) = 0.0d0
  do j = 1, nr
     if(j == 1) then
     elseif(j == nr) then
        dr = R(j) - R(j-1)
        eddy(j, j-1) = 0.5d0 * (nt(j-1) + nt(j)) / nt(j-1) / dr**2.0d0
        eddy(j, j) = - 0.5d0 * (nt(j-1) + nt(j)) / nt(j) / dr**2.0d0
     else
        dr = 0.5d0 * (R(j+1) - R(j-1))
        eddy(j, j-1) = 0.5d0 * (nt(j-1) + nt(j)) / nt(j-1) / dr**2.0d0
        eddy(j, j) = - (0.5d0 * (nt(j-1) + nt(j)) + 0.5d0 * (nt(j) + nt(j+1))) / nt(j) / dr**2.0d0
        eddy(j, j+1) = 0.5d0 * (nt(j) + nt(j+1)) / nt(j+1) / dr**2.0d0
     end if
  end do
  eddy(1:nr, 1:nr) = - Kzz * eddy(1:nr, 1:nr)
   
end subroutine diffusion
