	character*500 function VersionGIT()
	IMPLICIT NONE
	character*500 gitversion
#include "gitversion.h"

	VersionGIT=trim(gitversion)
	
	return
	end
	
