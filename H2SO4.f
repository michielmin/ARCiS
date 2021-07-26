c Data from Palmer (1975)
	subroutine H2SO4(l_lnk,n_lnk,k_lnk,nlam)
	IMPLICIT NONE
	integer nlam,j
	real*8 l_lnk(*),n_lnk(*),k_lnk(*)
	real*8 l0(227),n0(227),k0(227)
	data (l0(j),j=1,227) /
     &  0.35971E+00, 0.40816E+00, 0.44944E+00, 0.55556E+00, 0.70175E+00,
     &  0.71429E+00, 0.72464E+00, 0.73529E+00, 0.74627E+00, 0.75758E+00,
     &  0.76923E+00, 0.78125E+00, 0.79365E+00, 0.80645E+00, 0.81967E+00,
     &  0.83333E+00, 0.84746E+00, 0.86207E+00, 0.87719E+00, 0.89286E+00,
     &  0.90909E+00, 0.92593E+00, 0.94340E+00, 0.96154E+00, 0.98039E+00,
     &  0.10000E+01, 0.10204E+01, 0.10417E+01, 0.10638E+01, 0.10870E+01,
     &  0.11111E+01, 0.11364E+01, 0.11628E+01, 0.11905E+01, 0.12195E+01,
     &  0.12500E+01, 0.12658E+01, 0.12821E+01, 0.12987E+01, 0.13158E+01,
     &  0.13333E+01, 0.13514E+01, 0.13699E+01, 0.13889E+01, 0.14085E+01,
     &  0.14286E+01, 0.14493E+01, 0.14706E+01, 0.14925E+01, 0.15152E+01,
     &  0.15385E+01, 0.15625E+01, 0.15873E+01, 0.16129E+01, 0.16393E+01,
     &  0.16667E+01, 0.16949E+01, 0.17241E+01, 0.17544E+01, 0.17857E+01,
     &  0.18182E+01, 0.18519E+01, 0.18868E+01, 0.19231E+01, 0.19608E+01,
     &  0.20000E+01, 0.20408E+01, 0.20833E+01, 0.21277E+01, 0.21739E+01,
     &  0.22222E+01, 0.22727E+01, 0.23256E+01, 0.23810E+01, 0.24390E+01,
     &  0.25000E+01, 0.25641E+01, 0.26316E+01, 0.26882E+01, 0.27248E+01,
     &  0.27624E+01, 0.27701E+01, 0.28329E+01, 0.28409E+01, 0.28818E+01,
     &  0.29155E+01, 0.29412E+01, 0.29851E+01, 0.30211E+01, 0.30769E+01,
     &  0.31746E+01, 0.32787E+01, 0.33445E+01, 0.34130E+01, 0.34722E+01,
     &  0.35587E+01, 0.36232E+01, 0.36900E+01, 0.37594E+01, 0.38168E+01,
     &  0.38462E+01, 0.38610E+01, 0.39062E+01, 0.39526E+01, 0.40000E+01,
     &  0.40816E+01, 0.41494E+01, 0.42194E+01, 0.42735E+01, 0.42918E+01,
     &  0.43668E+01, 0.44643E+01, 0.45872E+01, 0.47170E+01, 0.48077E+01,
     &  0.49505E+01, 0.51020E+01, 0.51813E+01, 0.52632E+01, 0.53191E+01,
     &  0.53763E+01, 0.54348E+01, 0.54945E+01, 0.55556E+01, 0.56180E+01,
     &  0.56497E+01, 0.56818E+01, 0.57471E+01, 0.58140E+01, 0.58824E+01,
     &  0.59524E+01, 0.60241E+01, 0.60606E+01, 0.60976E+01, 0.61350E+01,
     &  0.61728E+01, 0.62112E+01, 0.62500E+01, 0.63291E+01, 0.64103E+01,
     &  0.64935E+01, 0.65359E+01, 0.65789E+01, 0.66225E+01, 0.66667E+01,
     &  0.67114E+01, 0.68027E+01, 0.68966E+01, 0.69930E+01, 0.70922E+01,
     &  0.71942E+01, 0.72993E+01, 0.73529E+01, 0.74627E+01, 0.75758E+01,
     &  0.76336E+01, 0.77519E+01, 0.78740E+01, 0.80000E+01, 0.80645E+01,
     &  0.81301E+01, 0.82645E+01, 0.84034E+01, 0.85470E+01, 0.86207E+01,
     &  0.86957E+01, 0.88496E+01, 0.89286E+01, 0.90090E+01, 0.90909E+01,
     &  0.91743E+01, 0.92593E+01, 0.93458E+01, 0.94340E+01, 0.95238E+01,
     &  0.96154E+01, 0.97087E+01, 0.98039E+01, 0.99010E+01, 0.10101E+02,
     &  0.10204E+02, 0.10309E+02, 0.10384E+02, 0.10417E+02, 0.10526E+02,
     &  0.10638E+02, 0.10753E+02, 0.10870E+02, 0.10989E+02, 0.11111E+02,
     &  0.11236E+02, 0.11364E+02, 0.11442E+02, 0.11494E+02, 0.11765E+02,
     &  0.11905E+02, 0.12195E+02, 0.12500E+02, 0.12658E+02, 0.12821E+02,
     &  0.13158E+02, 0.13514E+02, 0.13889E+02, 0.14286E+02, 0.14706E+02,
     &  0.14925E+02, 0.15385E+02, 0.15873E+02, 0.16129E+02, 0.16667E+02,
     &  0.16949E+02, 0.17241E+02, 0.17544E+02, 0.17857E+02, 0.18182E+02,
     &  0.18519E+02, 0.18868E+02, 0.19608E+02, 0.20000E+02, 0.20408E+02,
     &  0.20833E+02, 0.21277E+02, 0.22222E+02, 0.22727E+02, 0.23256E+02,
     &  0.24390E+02, 0.25000E+02 /
	data (n0(j),j=1,227) /
     &  0.14590E+01, 0.14430E+01, 0.14380E+01, 0.14340E+01, 0.14320E+01,
     &  0.14320E+01, 0.14320E+01, 0.14310E+01, 0.14310E+01, 0.14310E+01,
     &  0.14310E+01, 0.14310E+01, 0.14300E+01, 0.14300E+01, 0.14300E+01,
     &  0.14300E+01, 0.14300E+01, 0.14290E+01, 0.14290E+01, 0.14280E+01,
     &  0.14280E+01, 0.14280E+01, 0.14270E+01, 0.14270E+01, 0.14270E+01,
     &  0.14270E+01, 0.14270E+01, 0.14260E+01, 0.14260E+01, 0.14250E+01,
     &  0.14250E+01, 0.14240E+01, 0.14230E+01, 0.14230E+01, 0.14220E+01,
     &  0.14220E+01, 0.14220E+01, 0.14210E+01, 0.14200E+01, 0.14200E+01,
     &  0.14190E+01, 0.14180E+01, 0.14170E+01, 0.14160E+01, 0.14160E+01,
     &  0.14160E+01, 0.14150E+01, 0.14140E+01, 0.14130E+01, 0.14120E+01,
     &  0.14110E+01, 0.14100E+01, 0.14100E+01, 0.14090E+01, 0.14080E+01,
     &  0.14070E+01, 0.14060E+01, 0.14050E+01, 0.14040E+01, 0.14040E+01,
     &  0.14030E+01, 0.14010E+01, 0.14000E+01, 0.13990E+01, 0.13980E+01,
     &  0.13960E+01, 0.13940E+01, 0.13930E+01, 0.13910E+01, 0.13880E+01,
     &  0.13850E+01, 0.13830E+01, 0.13790E+01, 0.13760E+01, 0.13720E+01,
     &  0.13680E+01, 0.13600E+01, 0.13530E+01, 0.13450E+01, 0.13390E+01,
     &  0.13290E+01, 0.13270E+01, 0.13050E+01, 0.13030E+01, 0.12910E+01,
     &  0.12830E+01, 0.12760E+01, 0.12670E+01, 0.12630E+01, 0.12600E+01,
     &  0.12740E+01, 0.13150E+01, 0.13470E+01, 0.13770E+01, 0.14030E+01,
     &  0.14310E+01, 0.14440E+01, 0.14480E+01, 0.14450E+01, 0.14380E+01,
     &  0.14370E+01, 0.14370E+01, 0.14380E+01, 0.14360E+01, 0.14350E+01,
     &  0.14380E+01, 0.14460E+01, 0.14610E+01, 0.14720E+01, 0.14700E+01,
     &  0.14570E+01, 0.14460E+01, 0.14410E+01, 0.14380E+01, 0.14340E+01,
     &  0.14260E+01, 0.14160E+01, 0.14090E+01, 0.14030E+01, 0.13970E+01,
     &  0.13930E+01, 0.13870E+01, 0.13840E+01, 0.13790E+01, 0.13750E+01,
     &  0.13730E+01, 0.13710E+01, 0.13690E+01, 0.13650E+01, 0.13570E+01,
     &  0.13470E+01, 0.13410E+01, 0.13390E+01, 0.13360E+01, 0.13330E+01,
     &  0.13310E+01, 0.13290E+01, 0.13260E+01, 0.13220E+01, 0.13180E+01,
     &  0.13050E+01, 0.12970E+01, 0.12850E+01, 0.12700E+01, 0.12520E+01,
     &  0.12370E+01, 0.12090E+01, 0.11590E+01, 0.10920E+01, 0.10240E+01,
     &  0.10570E+01, 0.11800E+01, 0.12450E+01, 0.13140E+01, 0.12920E+01,
     &  0.12760E+01, 0.12440E+01, 0.12250E+01, 0.12180E+01, 0.12190E+01,
     &  0.12170E+01, 0.12300E+01, 0.13000E+01, 0.14470E+01, 0.15160E+01,
     &  0.15720E+01, 0.16320E+01, 0.16450E+01, 0.16540E+01, 0.16390E+01,
     &  0.16170E+01, 0.15780E+01, 0.15480E+01, 0.15370E+01, 0.15480E+01,
     &  0.15480E+01, 0.15270E+01, 0.14570E+01, 0.13720E+01, 0.13010E+01,
     &  0.14150E+01, 0.16340E+01, 0.17890E+01, 0.18560E+01, 0.19670E+01,
     &  0.19530E+01, 0.18700E+01, 0.18080E+01, 0.18410E+01, 0.19370E+01,
     &  0.20070E+01, 0.19780E+01, 0.19560E+01, 0.19400E+01, 0.18480E+01,
     &  0.18120E+01, 0.17510E+01, 0.17100E+01, 0.16930E+01, 0.16630E+01,
     &  0.16280E+01, 0.16040E+01, 0.15840E+01, 0.15670E+01, 0.15520E+01,
     &  0.15430E+01, 0.15200E+01, 0.14880E+01, 0.14660E+01, 0.14100E+01,
     &  0.13820E+01, 0.13890E+01, 0.14770E+01, 0.16800E+01, 0.19120E+01,
     &  0.20450E+01, 0.20570E+01, 0.19610E+01, 0.19130E+01, 0.18740E+01,
     &  0.18480E+01, 0.18260E+01, 0.17850E+01, 0.17810E+01, 0.18220E+01,
     &  0.18800E+01, 0.18960E+01 /
	data (k0(j),j=1,227) /
     &  0.14816E-17, 0.31946E-15, 0.12461E-13, 0.12621E-10, 0.56625E-08,
     &  0.85075E-08, 0.11782E-07, 0.16318E-07, 0.22600E-07, 0.31300E-07,
     &  0.39200E-07, 0.48500E-07, 0.58700E-07, 0.70600E-07, 0.83500E-07,
     &  0.99500E-07, 0.11700E-06, 0.13400E-06, 0.15200E-06, 0.17100E-06,
     &  0.19700E-06, 0.23900E-06, 0.32700E-06, 0.45300E-06, 0.62000E-06,
     &  0.86700E-06, 0.11100E-05, 0.13200E-05, 0.14000E-05, 0.14400E-05,
     &  0.15100E-05, 0.16900E-05, 0.19400E-05, 0.22900E-05, 0.29100E-05,
     &  0.39500E-05, 0.46800E-05, 0.56100E-05, 0.66700E-05, 0.78500E-05,
     &  0.94400E-05, 0.11500E-04, 0.13800E-04, 0.16500E-04, 0.19800E-04,
     &  0.25400E-04, 0.34000E-04, 0.45300E-04, 0.62100E-04, 0.78500E-04,
     &  0.98800E-04, 0.12400E-03, 0.15200E-03, 0.18100E-03, 0.21400E-03,
     &  0.25500E-03, 0.29900E-03, 0.34800E-03, 0.40200E-03, 0.45200E-03,
     &  0.49300E-03, 0.53300E-03, 0.58600E-03, 0.65800E-03, 0.74100E-03,
     &  0.83700E-03, 0.89300E-03, 0.94700E-03, 0.10200E-02, 0.11500E-02,
     &  0.12900E-02, 0.14100E-02, 0.16200E-02, 0.17100E-02, 0.18600E-02,
     &  0.21100E-02, 0.29290E-02, 0.40658E-02, 0.52856E-02, 0.62275E-02,
     &  0.73372E-02, 0.75818E-02, 0.98565E-02, 0.10185E-01, 0.12000E-01,
     &  0.22000E-01, 0.30000E-01, 0.48000E-01, 0.64000E-01, 0.92000E-01,
     &  0.14300E+00, 0.18200E+00, 0.19500E+00, 0.19700E+00, 0.19200E+00,
     &  0.17300E+00, 0.15500E+00, 0.13700E+00, 0.12200E+00, 0.11600E+00,
     &  0.11700E+00, 0.11800E+00, 0.11500E+00, 0.11300E+00, 0.11400E+00,
     &  0.11800E+00, 0.12000E+00, 0.11600E+00, 0.99000E-01, 0.92000E-01,
     &  0.77000E-01, 0.75000E-01, 0.74000E-01, 0.69000E-01, 0.64000E-01,
     &  0.58000E-01, 0.53000E-01, 0.52000E-01, 0.50000E-01, 0.51000E-01,
     &  0.51000E-01, 0.52000E-01, 0.53000E-01, 0.55000E-01, 0.56000E-01,
     &  0.57000E-01, 0.58000E-01, 0.58000E-01, 0.57000E-01, 0.54000E-01,
     &  0.57000E-01, 0.62000E-01, 0.64000E-01, 0.65000E-01, 0.67000E-01,
     &  0.68000E-01, 0.69000E-01, 0.69000E-01, 0.70000E-01, 0.65000E-01,
     &  0.56000E-01, 0.52000E-01, 0.48000E-01, 0.45000E-01, 0.47000E-01,
     &  0.52000E-01, 0.58000E-01, 0.67000E-01, 0.10000E+00, 0.22600E+00,
     &  0.33400E+00, 0.43700E+00, 0.43100E+00, 0.36200E+00, 0.31300E+00,
     &  0.30800E+00, 0.32100E+00, 0.35800E+00, 0.40000E+00, 0.42700E+00,
     &  0.45300E+00, 0.53200E+00, 0.64300E+00, 0.66900E+00, 0.65100E+00,
     &  0.61500E+00, 0.52400E+00, 0.48100E+00, 0.43900E+00, 0.39100E+00,
     &  0.36000E+00, 0.34200E+00, 0.34700E+00, 0.36600E+00, 0.37300E+00,
     &  0.35200E+00, 0.32600E+00, 0.31100E+00, 0.35000E+00, 0.62400E+00,
     &  0.80900E+00, 0.86900E+00, 0.83300E+00, 0.79500E+00, 0.62700E+00,
     &  0.46800E+00, 0.39100E+00, 0.40700E+00, 0.46000E+00, 0.45800E+00,
     &  0.35600E+00, 0.24300E+00, 0.20600E+00, 0.18100E+00, 0.12100E+00,
     &  0.10700E+00, 0.96000E-01, 0.94000E-01, 0.90000E-01, 0.90000E-01,
     &  0.11000E+00, 0.12600E+00, 0.14300E+00, 0.16000E+00, 0.17600E+00,
     &  0.18300E+00, 0.20300E+00, 0.23200E+00, 0.25300E+00, 0.34000E+00,
     &  0.41500E+00, 0.54000E+00, 0.69900E+00, 0.80200E+00, 0.74000E+00,
     &  0.56900E+00, 0.40500E+00, 0.24100E+00, 0.20900E+00, 0.19700E+00,
     &  0.19400E+00, 0.19300E+00, 0.21800E+00, 0.24800E+00, 0.27400E+00,
     &  0.24500E+00, 0.21200E+00 /
	nlam=227
	l_lnk(1:nlam)=l0(1:nlam)
	n_lnk(1:nlam)=n0(1:nlam)
	k_lnk(1:nlam)=k0(1:nlam)
	return
	end
