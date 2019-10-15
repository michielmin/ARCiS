! Include file for example MultiNest program 'Gaussian' (see arXiv:1001.0719)

module params_multinest
implicit none


! Toy Model Parameters

	!dimensionality
      	integer sdim
      !	parameter( sdim = 4 )
      

! Parameters for MultiNest
	
      	!whether to do use Nested Importance Sampling
	logical nest_IS
 	parameter(nest_IS=.true.)
	
      	!whether to do multimodal sampling
	logical nest_mmodal 
 	parameter(nest_mmodal=.true.)
	
      	!sample with constant efficiency
	logical nest_ceff
 	parameter(nest_ceff=.true.)
	
      	!max no. of live points
      	integer nest_nlive
	!parameter(nest_nlive=100)
      
      
      	!seed for MultiNest, -ve means take it from sys clock
	integer nest_rseed 
	parameter(nest_rseed=-1)
      
      	!evidence tolerance factor
      	double precision nest_tol 
!      	parameter(nest_tol=0.5)
      
      	!enlargement factor reduction parameter
      	double precision nest_efr
!      	parameter(nest_efr=0.3d0)
      
      	!root for saving posterior files
      	character*1000 nest_root
	!parameter(nest_root='./')
	
	!after how many iterations feedback is required & the output files should be updated
	!note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	integer nest_updInt
	parameter(nest_updInt=100)
	
	!null evidence (set it to very high negative no. if null evidence is unknown)
	double precision nest_Ztol
	parameter(nest_Ztol=-1.d90)
      
      	!max modes expected, for memory allocation
      	integer nest_maxModes 
      	parameter(nest_maxModes=10)
      
      	!no. of parameters to cluster (for mode detection)
      	integer nest_nClsPar
      	!parameter(nest_nClsPar=2)
      
      	!whether to resume from a previous run
      	logical nest_resume
      	!parameter(nest_resume=.false.)
      
      	!whether to write output files
      	logical nest_outfile
      	parameter(nest_outfile=.true.)
      
      	!initialize MPI routines?, relevant only if compiling with MPI
	!set it to F if you want your main program to handle MPI initialization
      	logical nest_initMPI
      	parameter(nest_initMPI=.false.)
      
      	!points with loglike < nest_logZero will be ignored by MultiNest
      	double precision nest_logZero
      	parameter(nest_logZero=-huge(1d0))
      
      	!max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
	!has done max no. of iterations or convergence criterion (defined through nest_tol) has been satisfied
      	integer nest_maxIter
      	parameter(nest_maxIter=0)
	
	
      	!feedback on the sampling progress?
      	logical nest_fb 
      	parameter(nest_fb=.true.)
!=======================================================================


end module params_multinest
