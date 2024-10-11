	subroutine outputform(string,form)
	IMPLICIT NONE
	character string*(*)
	character,intent(in),optional :: form*(*)

	if(present(form)) then
		write(*,form) string
		write(9,form) string
	else
		write(*,'(a)') string
		write(9,'(a)') string
	endif
	
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine flushoutput()
	IMPLICIT NONE
	
	call flush(6)
	call flush(9)

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	module OutputModeModule
	IMPLICIT NONE
	logical do_output
	end module OutputModeModule

	subroutine SetOutputMode(doit)
	use OutputModeModule
	IMPLICIT NONE
	logical doit
	do_output=doit
	return
	end

	subroutine output(string)
	use OutputModeModule
	IMPLICIT NONE
	character(len=*),intent(in) :: string

	if(.not.do_output) return

	write(*,'(a)') string
	write(9,'(a)') string
	
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine output_erase(string)
	use OutputModeModule
	IMPLICIT NONE
	character string*(*)
	if(.not.do_output) return

	write(*,'(1a1,a,$)') char(13),string
	
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine ignorestar(un)
	IMPLICIT NONE
	integer un
	character c
1	read(un,fmt=3,end=2) c
	if(c.eq.'*') goto 1
	if(c.eq.'#') goto 1
	backspace(unit=un)
2	continue
3	format(a1)
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine tellertje(i,n)
	use GlobalSetup
	use OutputModeModule
	IMPLICIT NONE
	integer i,n,f

	if(.not.do_output) return

	if(i.eq.1) call output("....................")

	f=int(20d0*dble(i)/dble(n))
	
	if(20d0*dble(i-1)/dble(n).lt.dble(f)
     &   .and.20d0*dble(i+1)/dble(n).gt.dble(f)) then
		call outputform(".",'(a1,$)')
		call flushoutput
	endif
	
	if(i.eq.n) call output("")

	return
	end

	logical function checktellertje(i,n)
	IMPLICIT NONE
	integer i,n,f

	checktellertje=.false.
	f=int(20d0*dble(i)/dble(n))
	
	if(20d0*dble(i-1)/dble(n).lt.dble(f)
     &   .and.20d0*dble(i+1)/dble(n).gt.dble(f)) then
		checktellertje=.true.
	endif

	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine tellertje_perc(i,n)
	use dbl2str_mod
	IMPLICIT NONE
	integer i,n,f
	real*8 xx
	
	f=int(100d0*dble(i)/dble(n))
	
	xx=100d0*dble(i)/dble(n)
	if(i.eq.1) then
		call output_erase(trim(dbl2string(xx,'(f5.1)')) // " %")
	else if(100d0*dble(i-1)/dble(n).lt.dble(f)
     &   .and.100d0*dble(i+1)/dble(n).gt.dble(f)) then
		call output_erase(trim(dbl2string(xx,'(f5.1)')) // " %")
	endif
	if(i.eq.n) then
		call output("")
	endif

	return
	end



