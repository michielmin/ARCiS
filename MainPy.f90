	module pyEx
	IMPLICIT NONE
	integer nlam,nr,nret,nobs
	end module pyEx
	

	subroutine pysetoutputmode(tf)
	IMPLICIT NONE
	logical tf

	call SetOutputMode(tf)
	
	return
	end
	

	subroutine pyInit(pyinputfile,pyoutputdir,cline)
	use pyEx, py_nlam => nlam, py_nr => nr, py_nret => nret, py_nobs => nobs
	IMPLICIT NONE
	character(*),optional :: cline
	character(*) pyinputfile,pyoutputdir
	integer i,nlam,nr,n_ret,nobs

	call pyInitSupport(pyinputfile,pyoutputdir,cline,nlam,nr,n_ret,nobs)

	py_nlam=nlam
	py_nr=nr
	py_nret=n_ret
	py_nobs=nobs

	return
	end

	subroutine pyRunARCiS()
	IMPLICIT NONE

	call pyRunARCiSSupport()

	return
	end
	
	subroutine pyComputeModel()
	IMPLICIT NONE

	call InitDens()
	call ComputeModel(.true.)
	
	return
	end
	
	subroutine pyWriteFiles()
	IMPLICIT NONE

	call WriteStructure()
	call WriteOutput()

	return
	end
	
	subroutine pyWriteBestfit()
	IMPLICIT NONE

	call pyWriteBestfitSupport()

	return
	end
	
	
	subroutine pySetKeyword(pykey,pyval)
	IMPLICIT NONE
	character(*) pykey,pyval
	character*1000 readline

	readline=trim(pykey) // "=" // trim(pyval)
	call GetAndSetKeyPy(readline)

	return
	end

	
	subroutine pySetValue(pykey,pyval)
	IMPLICIT NONE
	character(*) pykey
	real*8 pyval
	character*1000 readline,dbstr
	
	write(dbstr,'(es14.7)') pyval

	readline=trim(pykey) // "=" // trim(dbstr)
	call GetAndSetKeyPy(readline)

	return
	end

	function pyGetLike() result(lnew)
	IMPLICIT NONE
	real*8 lnew

	call ComputeLike(lnew)
 	if(.not.lnew.gt.-1d100) lnew=-1d100

	return
	end


	function pyGetRetrievalNames(n) result(names)
	IMPLICIT NONE
	integer,intent(in) :: n
	character*500 :: names,pyGetRetrievalNamesSupport

	names=pyGetRetrievalNamesSupport(n)	

	return
	end

	subroutine pyGetLam(l)
	use pyEx
	IMPLICIT NONE
	real*8 l(:)

	call pyGetLamSupport(l,nlam)

	return
	end

	subroutine pyGetTrans(l,trans)
	use pyEx
	IMPLICIT NONE
	real*8 l(:),trans(:)

	call pyGetTransSupport(l,trans,nlam)

	return
	end

	subroutine pyGetEmis(l,emis,ip)
	use pyEx
	IMPLICIT NONE
	real*8 l(:),emis(:)
	integer ip
	
	call pyGetEmisSupport(l,emis,ip,nlam)

	return
	end

	subroutine pyGetStar(l,star)
	use pyEx
	IMPLICIT NONE
	real*8 l(:),star(:)
	
	call pyGetStarSupport(l,star,nlam)

	return
	end


	subroutine pyGetPT(Ppy,Tpy)
	IMPLICIT NONE
	real*8 Ppy(:),Tpy(:)

	call pyGetPTSupport(Ppy,Tpy)

	return
	end



	function pyTransformRetPar(var,i) result(x)
	IMPLICIT NONE
	real*8 var,x,pyTransformRetParSupport
	integer i,j

	x=pyTransformRetParSupport(var,i)

	return
	end


	subroutine pyMapRetrieval(var,nret)
	IMPLICIT NONE
	integer nret
	real*8 var(nret),dvar(2,nret)
	
	call MapRetrievalMN(var,dvar)
	
	return
	end
	

	function pyGetObsN(i) result(n)
	integer n,pyGetObsNSupport
	n=pyGetObsNSupport(i)
	return
	end

	subroutine pyGetObs(i,l,obs,dobs,model)
	IMPLICIT NONE
	real*8 l(:),obs(:),model(:),dobs(:)
	integer i

	call pyGetObsSupport(i,l,obs,dobs,model)

	return
	end

	subroutine pyVerbose(switch)
	IMPLICIT NONE
	logical switch

	call pyVerboseSupport(switch)
	
	return
	end
	

