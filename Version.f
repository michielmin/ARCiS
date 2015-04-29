	subroutine VersionDateTime(string)
	character*500 version,string
	parameter(version=trim(__DATE__)//' '//trim(__TIME__))

	string=version
	
	return
	end
	
	
	subroutine VersionGIT(string)
	character*500 string
#include "gitversion.h"

	string=gitversion
	
	return
	end
	
