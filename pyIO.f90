subroutine ask_python(x,nx,y,ny,init,command,info)
  use iso_c_binding
  implicit none

  interface
     function py_popen(cmd, mode) bind(C, name="py_popen") result(fp)
       import :: c_ptr, c_char
       character(kind=c_char), dimension(*) :: cmd
       character(kind=c_char), dimension(*) :: mode
       type(c_ptr) :: fp
     end function

     function py_pclose(fp) bind(C, name="py_pclose") result(rc)
       import :: c_ptr, c_int
       type(c_ptr), value :: fp
       integer(c_int) :: rc
     end function

     function py_write_int(fp, v) bind(C, name="py_write_int") result(rc)
       import :: c_ptr, c_int
       type(c_ptr), value :: fp
       integer(c_int), value :: v
       integer(c_int) :: rc
     end function

     function py_write_double(fp, v) bind(C, name="py_write_double") result(rc)
       import :: c_ptr, c_double, c_int
       type(c_ptr), value :: fp
       real(c_double), value :: v
       integer(c_int) :: rc
     end function

     function py_flush(fp) bind(C, name="py_flush") result(rc)
       import :: c_ptr, c_int
       type(c_ptr), value :: fp
       integer(c_int) :: rc
     end function

     function py_read_int(fp, v) bind(C, name="py_read_int") result(rc)
       import :: c_ptr, c_int
       type(c_ptr), value :: fp
       integer(c_int) :: v
       integer(c_int) :: rc
     end function

     function py_read_double(fp, v) bind(C, name="py_read_double") result(rc)
       import :: c_ptr, c_double, c_int
       type(c_ptr), value :: fp
       real(c_double) :: v
       integer(c_int) :: rc
     end function
  end interface

  integer :: init
! init=0 -> open python interface
! init=1 -> do interface
! init=2 -> close interface
  integer(c_int) :: nx, ny
  real(c_double) :: x(nx), y(ny)
  integer :: i,j
  integer(c_int) :: nret, rc, info
  type(c_ptr),save :: fp
  character*500 :: command

  if(init.eq.0) then
  ! Start persistent python process (it will loop until pipe closes)
  	fp = py_popen(command, "r+")
  	if (.not. c_associated(fp)) then
     	print *, "ERROR: popen failed"
     	stop 1
  	end if
  endif

  if(init.eq.1) then
  ! Send request
!	rc = py_write_int(fp, nx)
  	do i = 1, nx
     	rc = py_write_double(fp, x(i))
  	end do
  	rc = py_flush(fp)

  ! Read response
!  	rc = py_read_int(fp, nret)
!  	if (nret /= ny) then
!     	print *, "ERROR: returned N differs: ", nret
!     	rc = py_pclose(fp)
!     	stop 1
!  	end if

  	do i = 1, ny
     	rc = py_read_double(fp, y(i))
  	end do
	rc = py_read_int(fp, info)
  endif

  if(init.eq.2) then
  	rc = py_pclose(fp)
  endif

end subroutine ask_python

