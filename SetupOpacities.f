	subroutine SetupOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer imol
	
	do imol=1,nmol
		select case(Mol(imol)%filetype)
			case("LAMBDA")
				call ReadLambdaFiles(imol)
			case default
				call output("Wrong filetype")
		end select
	enddo
	
	return
	end
	
	
