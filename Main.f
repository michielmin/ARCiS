	program SPARC
	IMPLICIT NONE
	
	call Init()
	
	call SetupStructure()
	
	call SetupLTE()
	
	call Raytrace()
	
	call output()
	
	end
	
