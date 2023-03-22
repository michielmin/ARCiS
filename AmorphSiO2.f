cRefractive index of amorphous SiO2[s] at 300K; Henning&Mutschke (1997), A&A Vol. 327; Philipp (1985) in Palik: "Handbook of Optical Constants of Solids"; Database of Optical Constants for Cosmic Dust, Laboratory Astrophysics Group of the AIU Jena
cData points: 434
cWavelength (micron) n  k
	subroutine AmorphSiO2(l_lnk,n_lnk,k_lnk,nlam)
	IMPLICIT NONE
	integer nlam,j
	real*8 l_lnk(*),n_lnk(*),k_lnk(*)
	real*8 l0(         434),n0(         434),k0(         434)
	data (l0(j),j=1,         434) /
     &  0.47690E-01, 0.48620E-01, 0.49590E-01, 0.50610E-01, 0.51660E-01, 
     &  0.52760E-01, 0.53910E-01, 0.55100E-01, 0.56360E-01, 0.57670E-01, 
     &  0.59040E-01, 0.60480E-01, 0.61990E-01, 0.63580E-01, 0.65260E-01, 
     &  0.67020E-01, 0.68880E-01, 0.69850E-01, 0.70850E-01, 0.71250E-01, 
     &  0.71870E-01, 0.72930E-01, 0.74020E-01, 0.75140E-01, 0.76290E-01, 
     &  0.77490E-01, 0.78720E-01, 0.79990E-01, 0.81300E-01, 0.82660E-01, 
     &  0.83770E-01, 0.84920E-01, 0.86100E-01, 0.86700E-01, 0.87310E-01, 
     &  0.88560E-01, 0.89840E-01, 0.91160E-01, 0.92520E-01, 0.93920E-01, 
     &  0.95370E-01, 0.96860E-01, 0.98390E-01, 0.99980E-01, 0.10160E+00, 
     &  0.10330E+00, 0.10420E+00, 0.10510E+00, 0.10600E+00, 0.10690E+00, 
     &  0.10780E+00, 0.10880E+00, 0.10970E+00, 0.11070E+00, 0.11170E+00, 
     &  0.11270E+00, 0.11370E+00, 0.11430E+00, 0.11450E+00, 0.11480E+00, 
     &  0.11530E+00, 0.11590E+00, 0.11640E+00, 0.11700E+00, 0.11750E+00, 
     &  0.11810E+00, 0.11860E+00, 0.11920E+00, 0.11960E+00, 0.11980E+00, 
     &  0.12040E+00, 0.12100E+00, 0.12150E+00, 0.12270E+00, 0.12400E+00, 
     &  0.12520E+00, 0.12650E+00, 0.12780E+00, 0.12910E+00, 0.13050E+00, 
     &  0.13190E+00, 0.13330E+00, 0.13480E+00, 0.13620E+00, 0.13780E+00, 
     &  0.13930E+00, 0.14090E+00, 0.14250E+00, 0.14420E+00, 0.14590E+00, 
     &  0.14760E+00, 0.14940E+00, 0.15120E+00, 0.15310E+00, 0.15500E+00, 
     &  0.15900E+00, 0.16310E+00, 0.16750E+00, 0.17220E+00, 0.17710E+00, 
     &  0.18370E+00, 0.18510E+00, 0.18547E+00, 0.19358E+00, 0.19800E+00, 
     &  0.20006E+00, 0.20255E+00, 0.20445E+00, 0.21107E+00, 0.21444E+00, 
     &  0.21946E+00, 0.22650E+00, 0.23100E+00, 0.23129E+00, 0.24280E+00, 
     &  0.25033E+00, 0.25730E+00, 0.26315E+00, 0.27487E+00, 0.29136E+00, 
     &  0.30341E+00, 0.31228E+00, 0.32525E+00, 0.34000E+00, 0.34036E+00, 
     &  0.35868E+00, 0.39400E+00, 0.39685E+00, 0.40466E+00, 0.41017E+00, 
     &  0.43400E+00, 0.43405E+00, 0.43583E+00, 0.46781E+00, 0.47999E+00, 
     &  0.48613E+00, 0.50800E+00, 0.50858E+00, 0.51836E+00, 0.53385E+00, 
     &  0.54607E+00, 0.57907E+00, 0.58756E+00, 0.58929E+00, 0.58930E+00, 
     &  0.62782E+00, 0.64385E+00, 0.65628E+00, 0.66782E+00, 0.67079E+00, 
     &  0.70652E+00, 0.72813E+00, 0.76649E+00, 0.76800E+00, 0.79476E+00, 
     &  0.83250E+00, 0.84467E+00, 0.99140E+00, 0.10000E+01, 0.10141E+01, 
     &  0.10830E+01, 0.11592E+01, 0.12000E+01, 0.13000E+01, 0.13070E+01, 
     &  0.13958E+01, 0.14000E+01, 0.14792E+01, 0.15296E+01, 0.15414E+01, 
     &  0.16000E+01, 0.16815E+01, 0.17614E+01, 0.18000E+01, 0.19457E+01, 
     &  0.20531E+01, 0.20582E+01, 0.23000E+01, 0.25000E+01, 0.26000E+01, 
     &  0.30000E+01, 0.35000E+01, 0.40000E+01, 0.42000E+01, 0.44000E+01, 
     &  0.46000E+01, 0.48000E+01, 0.50000E+01, 0.52000E+01, 0.54000E+01, 
     &  0.56000E+01, 0.58000E+01, 0.60000E+01, 0.61000E+01, 0.62000E+01, 
     &  0.63000E+01, 0.64000E+01, 0.65000E+01, 0.66000E+01, 0.67000E+01, 
     &  0.68000E+01, 0.69000E+01, 0.70000E+01, 0.70500E+01, 0.71000E+01, 
     &  0.71500E+01, 0.72000E+01, 0.72500E+01, 0.73000E+01, 0.73500E+01, 
     &  0.74000E+01, 0.74500E+01, 0.75000E+01, 0.75500E+01, 0.76000E+01, 
     &  0.76500E+01, 0.77000E+01, 0.77500E+01, 0.78000E+01, 0.78500E+01, 
     &  0.79000E+01, 0.79500E+01, 0.79526E+01, 0.79898E+01, 0.80274E+01, 
     &  0.80654E+01, 0.81037E+01, 0.81424E+01, 0.81814E+01, 0.82209E+01, 
     &  0.82607E+01, 0.83009E+01, 0.83415E+01, 0.83824E+01, 0.84239E+01, 
     &  0.84657E+01, 0.85079E+01, 0.85505E+01, 0.85937E+01, 0.86372E+01, 
     &  0.86811E+01, 0.87255E+01, 0.87703E+01, 0.88157E+01, 0.88615E+01, 
     &  0.89077E+01, 0.89545E+01, 0.90018E+01, 0.90495E+01, 0.90978E+01, 
     &  0.91465E+01, 0.91958E+01, 0.92457E+01, 0.92961E+01, 0.93470E+01, 
     &  0.93985E+01, 0.94506E+01, 0.95033E+01, 0.95565E+01, 0.96103E+01, 
     &  0.96648E+01, 0.97199E+01, 0.97756E+01, 0.98319E+01, 0.98888E+01, 
     &  0.99466E+01, 0.10005E+02, 0.10064E+02, 0.10124E+02, 0.10184E+02, 
     &  0.10245E+02, 0.10307E+02, 0.10370E+02, 0.10433E+02, 0.10497E+02, 
     &  0.10562E+02, 0.10628E+02, 0.10695E+02, 0.10762E+02, 0.10831E+02, 
     &  0.10900E+02, 0.10970E+02, 0.11041E+02, 0.11113E+02, 0.11186E+02, 
     &  0.11260E+02, 0.11334E+02, 0.11410E+02, 0.11487E+02, 0.11565E+02, 
     &  0.11644E+02, 0.11724E+02, 0.11805E+02, 0.11887E+02, 0.11971E+02, 
     &  0.12055E+02, 0.12141E+02, 0.12228E+02, 0.12316E+02, 0.12406E+02, 
     &  0.12497E+02, 0.12589E+02, 0.12683E+02, 0.12778E+02, 0.12874E+02, 
     &  0.12972E+02, 0.13072E+02, 0.13172E+02, 0.13275E+02, 0.13379E+02, 
     &  0.13485E+02, 0.13592E+02, 0.13701E+02, 0.13812E+02, 0.13925E+02, 
     &  0.14040E+02, 0.14156E+02, 0.14275E+02, 0.14395E+02, 0.14518E+02, 
     &  0.14642E+02, 0.14769E+02, 0.14898E+02, 0.15029E+02, 0.15163E+02, 
     &  0.15299E+02, 0.15437E+02, 0.15578E+02, 0.15722E+02, 0.15868E+02, 
     &  0.16017E+02, 0.16169E+02, 0.16324E+02, 0.16482E+02, 0.16642E+02, 
     &  0.16806E+02, 0.16974E+02, 0.17144E+02, 0.17318E+02, 0.17496E+02, 
     &  0.17677E+02, 0.17862E+02, 0.18051E+02, 0.18244E+02, 0.18441E+02, 
     &  0.18643E+02, 0.18849E+02, 0.19060E+02, 0.19275E+02, 0.19495E+02, 
     &  0.19721E+02, 0.19951E+02, 0.20187E+02, 0.20429E+02, 0.20677E+02, 
     &  0.20930E+02, 0.21191E+02, 0.21457E+02, 0.21730E+02, 0.22011E+02, 
     &  0.22299E+02, 0.22594E+02, 0.22897E+02, 0.23209E+02, 0.23529E+02, 
     &  0.23858E+02, 0.24196E+02, 0.24544E+02, 0.24903E+02, 0.25272E+02, 
     &  0.25652E+02, 0.26043E+02, 0.26447E+02, 0.26863E+02, 0.27293E+02, 
     &  0.27737E+02, 0.28196E+02, 0.28669E+02, 0.29160E+02, 0.29667E+02, 
     &  0.30192E+02, 0.30736E+02, 0.31300E+02, 0.31885E+02, 0.32492E+02, 
     &  0.33123E+02, 0.33779E+02, 0.34461E+02, 0.35172E+02, 0.35912E+02, 
     &  0.36685E+02, 0.37491E+02, 0.38333E+02, 0.39215E+02, 0.40137E+02, 
     &  0.41104E+02, 0.42119E+02, 0.43186E+02, 0.44307E+02, 0.45489E+02, 
     &  0.46735E+02, 0.48052E+02, 0.49444E+02, 0.50920E+02, 0.52487E+02, 
     &  0.54154E+02, 0.55929E+02, 0.57825E+02, 0.59854E+02, 0.62030E+02, 
     &  0.64371E+02, 0.66895E+02, 0.69626E+02, 0.72588E+02, 0.75815E+02, 
     &  0.79341E+02, 0.83211E+02, 0.87478E+02, 0.92207E+02, 0.97476E+02, 
     &  0.10338E+03, 0.11005E+03, 0.11764E+03, 0.12636E+03, 0.13647E+03, 
     &  0.14833E+03, 0.16246E+03, 0.17956E+03, 0.20069E+03, 0.22744E+03, 
     &  0.26244E+03, 0.31015E+03, 0.37907E+03, 0.48738E+03 /
	data (n0(j),j=1,         434) /
     &  0.62850E+00, 0.64320E+00, 0.65170E+00, 0.65910E+00, 0.67230E+00, 
     &  0.68620E+00, 0.69980E+00, 0.71380E+00, 0.73040E+00, 0.74700E+00, 
     &  0.76640E+00, 0.79000E+00, 0.81310E+00, 0.83710E+00, 0.86420E+00, 
     &  0.89490E+00, 0.92890E+00, 0.94990E+00, 0.98940E+00, 0.10130E+01, 
     &  0.10520E+01, 0.11180E+01, 0.11690E+01, 0.12050E+01, 0.12240E+01, 
     &  0.12210E+01, 0.12100E+01, 0.12020E+01, 0.11910E+01, 0.11840E+01, 
     &  0.11870E+01, 0.12170E+01, 0.12850E+01, 0.13350E+01, 0.13860E+01, 
     &  0.14570E+01, 0.14840E+01, 0.14690E+01, 0.14470E+01, 0.14300E+01, 
     &  0.14100E+01, 0.13940E+01, 0.13900E+01, 0.14000E+01, 0.14350E+01, 
     &  0.15040E+01, 0.15510E+01, 0.16020E+01, 0.16610E+01, 0.17220E+01, 
     &  0.17820E+01, 0.18390E+01, 0.18920E+01, 0.19240E+01, 0.19250E+01, 
     &  0.18910E+01, 0.18050E+01, 0.17380E+01, 0.16880E+01, 0.16260E+01, 
     &  0.15690E+01, 0.15280E+01, 0.15060E+01, 0.15210E+01, 0.16070E+01, 
     &  0.17690E+01, 0.19590E+01, 0.21590E+01, 0.22890E+01, 0.23730E+01, 
     &  0.25150E+01, 0.25900E+01, 0.26190E+01, 0.26210E+01, 0.25810E+01, 
     &  0.25290E+01, 0.24640E+01, 0.23970E+01, 0.23350E+01, 0.22750E+01, 
     &  0.22220E+01, 0.21740E+01, 0.21330E+01, 0.20950E+01, 0.20600E+01, 
     &  0.20280E+01, 0.19990E+01, 0.19710E+01, 0.19460E+01, 0.19220E+01, 
     &  0.18990E+01, 0.18770E+01, 0.18560E+01, 0.18380E+01, 0.18200E+01, 
     &  0.17890E+01, 0.17610E+01, 0.17370E+01, 0.17160E+01, 0.16980E+01, 
     &  0.16790E+01, 0.16775E+01, 0.16758E+01, 0.16600E+01, 0.16509E+01, 
     &  0.16493E+01, 0.16456E+01, 0.16429E+01, 0.16343E+01, 0.16304E+01, 
     &  0.16250E+01, 0.16182E+01, 0.16139E+01, 0.16140E+01, 0.16053E+01, 
     &  0.16003E+01, 0.15962E+01, 0.15931E+01, 0.15875E+01, 0.15810E+01, 
     &  0.15796E+01, 0.15743E+01, 0.15709E+01, 0.15675E+01, 0.15675E+01, 
     &  0.15639E+01, 0.15585E+01, 0.15581E+01, 0.15572E+01, 0.15565E+01, 
     &  0.15540E+01, 0.15540E+01, 0.15538E+01, 0.15510E+01, 0.15501E+01, 
     &  0.15497E+01, 0.15482E+01, 0.15482E+01, 0.15477E+01, 0.15468E+01, 
     &  0.15462E+01, 0.15447E+01, 0.15443E+01, 0.15442E+01, 0.15442E+01, 
     &  0.15428E+01, 0.15423E+01, 0.15419E+01, 0.15416E+01, 0.15415E+01, 
     &  0.15405E+01, 0.15399E+01, 0.15391E+01, 0.15390E+01, 0.15385E+01, 
     &  0.15377E+01, 0.15375E+01, 0.15351E+01, 0.15350E+01, 0.15348E+01, 
     &  0.15339E+01, 0.15328E+01, 0.15323E+01, 0.15310E+01, 0.15309E+01, 
     &  0.15298E+01, 0.15297E+01, 0.15287E+01, 0.15280E+01, 0.15278E+01, 
     &  0.15270E+01, 0.15258E+01, 0.15247E+01, 0.15241E+01, 0.15218E+01, 
     &  0.15200E+01, 0.15200E+01, 0.15156E+01, 0.15116E+01, 0.15099E+01, 
     &  0.14995E+01, 0.14845E+01, 0.14662E+01, 0.14570E+01, 0.14480E+01, 
     &  0.14370E+01, 0.14250E+01, 0.14120E+01, 0.13980E+01, 0.13820E+01, 
     &  0.13640E+01, 0.13430E+01, 0.13200E+01, 0.13070E+01, 0.12940E+01, 
     &  0.12790E+01, 0.12630E+01, 0.12450E+01, 0.12260E+01, 0.12060E+01, 
     &  0.11830E+01, 0.11580E+01, 0.11300E+01, 0.11150E+01, 0.10990E+01, 
     &  0.10830E+01, 0.10650E+01, 0.10460E+01, 0.10250E+01, 0.10040E+01, 
     &  0.98010E+00, 0.95490E+00, 0.92770E+00, 0.89820E+00, 0.86600E+00, 
     &  0.83080E+00, 0.79210E+00, 0.74920E+00, 0.70130E+00, 0.64720E+00, 
     &  0.59530E+00, 0.51320E+00, 0.42378E+00, 0.39963E+00, 0.39185E+00, 
     &  0.39461E+00, 0.40300E+00, 0.41541E+00, 0.42963E+00, 0.44294E+00, 
     &  0.45575E+00, 0.46706E+00, 0.47370E+00, 0.47837E+00, 0.48130E+00, 
     &  0.48035E+00, 0.47706E+00, 0.47287E+00, 0.46379E+00, 0.45171E+00, 
     &  0.43931E+00, 0.42327E+00, 0.40994E+00, 0.41375E+00, 0.44297E+00, 
     &  0.52178E+00, 0.68759E+00, 0.93926E+00, 0.12512E+01, 0.15842E+01, 
     &  0.18916E+01, 0.21557E+01, 0.23871E+01, 0.25714E+01, 0.27163E+01, 
     &  0.28404E+01, 0.29278E+01, 0.29756E+01, 0.30044E+01, 0.29977E+01, 
     &  0.29606E+01, 0.29073E+01, 0.28239E+01, 0.27239E+01, 0.26289E+01, 
     &  0.25325E+01, 0.24387E+01, 0.23653E+01, 0.22972E+01, 0.22338E+01, 
     &  0.21868E+01, 0.21433E+01, 0.20974E+01, 0.20625E+01, 0.20305E+01, 
     &  0.19939E+01, 0.19655E+01, 0.19401E+01, 0.19117E+01, 0.18856E+01, 
     &  0.18659E+01, 0.18410E+01, 0.18177E+01, 0.17979E+01, 0.17728E+01, 
     &  0.17516E+01, 0.17315E+01, 0.17053E+01, 0.16770E+01, 0.16508E+01, 
     &  0.16220E+01, 0.15964E+01, 0.15819E+01, 0.15740E+01, 0.15675E+01, 
     &  0.15725E+01, 0.15876E+01, 0.16060E+01, 0.16332E+01, 0.16716E+01, 
     &  0.17074E+01, 0.17419E+01, 0.17705E+01, 0.17900E+01, 0.17865E+01, 
     &  0.17729E+01, 0.17438E+01, 0.17209E+01, 0.17026E+01, 0.16846E+01, 
     &  0.16647E+01, 0.16490E+01, 0.16355E+01, 0.16236E+01, 0.16117E+01, 
     &  0.15996E+01, 0.15874E+01, 0.15748E+01, 0.15662E+01, 0.15585E+01, 
     &  0.15498E+01, 0.15362E+01, 0.15262E+01, 0.15131E+01, 0.14985E+01, 
     &  0.14860E+01, 0.14716E+01, 0.14576E+01, 0.14438E+01, 0.14291E+01, 
     &  0.14136E+01, 0.14013E+01, 0.13872E+01, 0.13659E+01, 0.13491E+01, 
     &  0.13325E+01, 0.13157E+01, 0.13027E+01, 0.12870E+01, 0.12650E+01, 
     &  0.12437E+01, 0.12222E+01, 0.12017E+01, 0.11720E+01, 0.11393E+01, 
     &  0.11011E+01, 0.10392E+01, 0.96123E+00, 0.87050E+00, 0.76737E+00, 
     &  0.66298E+00, 0.57067E+00, 0.51942E+00, 0.53220E+00, 0.63417E+00, 
     &  0.87537E+00, 0.12768E+01, 0.17699E+01, 0.22249E+01, 0.25482E+01, 
     &  0.27308E+01, 0.28115E+01, 0.28193E+01, 0.27852E+01, 0.27419E+01, 
     &  0.26852E+01, 0.26247E+01, 0.25743E+01, 0.25183E+01, 0.24660E+01, 
     &  0.24331E+01, 0.23950E+01, 0.23590E+01, 0.23346E+01, 0.23104E+01, 
     &  0.22838E+01, 0.22672E+01, 0.22494E+01, 0.22240E+01, 0.22039E+01, 
     &  0.21883E+01, 0.21620E+01, 0.21481E+01, 0.21387E+01, 0.21258E+01, 
     &  0.21168E+01, 0.21140E+01, 0.21048E+01, 0.20956E+01, 0.20899E+01, 
     &  0.20795E+01, 0.20653E+01, 0.20568E+01, 0.20406E+01, 0.20265E+01, 
     &  0.20229E+01, 0.20161E+01, 0.20069E+01, 0.19999E+01, 0.19987E+01, 
     &  0.20024E+01, 0.19986E+01, 0.19983E+01, 0.19875E+01, 0.19867E+01, 
     &  0.19800E+01, 0.19689E+01, 0.19645E+01, 0.19592E+01, 0.19495E+01, 
     &  0.19481E+01, 0.19519E+01, 0.19449E+01, 0.19426E+01, 0.19445E+01, 
     &  0.19364E+01, 0.19322E+01, 0.19326E+01, 0.19258E+01, 0.19186E+01, 
     &  0.19204E+01, 0.19144E+01, 0.19060E+01, 0.19184E+01, 0.19329E+01, 
     &  0.19254E+01, 0.19191E+01, 0.19255E+01, 0.19209E+01, 0.19035E+01, 
     &  0.19034E+01, 0.19033E+01, 0.19033E+01, 0.19034E+01 /
	data (k0(j),j=1,         434) /
     &  0.32430E+00, 0.36090E+00, 0.39070E+00, 0.42310E+00, 0.45400E+00, 
     &  0.48220E+00, 0.50970E+00, 0.53870E+00, 0.56670E+00, 0.59590E+00, 
     &  0.62720E+00, 0.65590E+00, 0.68440E+00, 0.71540E+00, 0.74880E+00, 
     &  0.78470E+00, 0.82730E+00, 0.86190E+00, 0.89870E+00, 0.91080E+00, 
     &  0.92360E+00, 0.91530E+00, 0.89170E+00, 0.85990E+00, 0.82230E+00, 
     &  0.78960E+00, 0.78480E+00, 0.78410E+00, 0.79680E+00, 0.82650E+00, 
     &  0.86940E+00, 0.92900E+00, 0.97240E+00, 0.98290E+00, 0.97090E+00, 
     &  0.91870E+00, 0.84620E+00, 0.79100E+00, 0.77610E+00, 0.77150E+00, 
     &  0.77530E+00, 0.80170E+00, 0.84060E+00, 0.89080E+00, 0.95730E+00, 
     &  0.10150E+01, 0.10350E+01, 0.10480E+01, 0.10560E+01, 0.10460E+01, 
     &  0.10270E+01, 0.99390E+00, 0.94520E+00, 0.87300E+00, 0.79580E+00, 
     &  0.72100E+00, 0.65970E+00, 0.64180E+00, 0.63560E+00, 0.67320E+00, 
     &  0.75940E+00, 0.84630E+00, 0.96660E+00, 0.11140E+01, 0.12890E+01, 
     &  0.14000E+01, 0.14540E+01, 0.14520E+01, 0.14230E+01, 0.13780E+01, 
     &  0.12290E+01, 0.10720E+01, 0.93260E+00, 0.71230E+00, 0.54410E+00, 
     &  0.41610E+00, 0.31480E+00, 0.24020E+00, 0.17000E+00, 0.12000E+00, 
     &  0.82000E-01, 0.56000E-01, 0.36000E-01, 0.23000E-01, 0.13000E-01, 
     &  0.70000E-02, 0.31000E-02, 0.10000E-02, 0.23000E-03, 0.35000E-04, 
     &  0.27000E-05, 0.27002E-05, 0.27008E-05, 0.27018E-05, 0.27033E-05, 
     &  0.27076E-05, 0.27137E-05, 0.27220E-05, 0.27328E-05, 0.27459E-05, 
     &  0.27665E-05, 0.27712E-05, 0.27725E-05, 0.28027E-05, 0.28209E-05, 
     &  0.28298E-05, 0.28408E-05, 0.28494E-05, 0.28810E-05, 0.28980E-05, 
     &  0.29242E-05, 0.29630E-05, 0.29890E-05, 0.29907E-05, 0.30610E-05, 
     &  0.31098E-05, 0.31570E-05, 0.31979E-05, 0.32833E-05, 0.34113E-05, 
     &  0.35103E-05, 0.35859E-05, 0.37008E-05, 0.38374E-05, 0.38408E-05, 
     &  0.40192E-05, 0.43897E-05, 0.44210E-05, 0.45082E-05, 0.45708E-05, 
     &  0.48508E-05, 0.48513E-05, 0.48729E-05, 0.52750E-05, 0.54356E-05, 
     &  0.55182E-05, 0.58210E-05, 0.58292E-05, 0.59692E-05, 0.61966E-05, 
     &  0.63809E-05, 0.69008E-05, 0.70399E-05, 0.70685E-05, 0.70687E-05, 
     &  0.77296E-05, 0.80184E-05, 0.82480E-05, 0.84656E-05, 0.85223E-05, 
     &  0.92278E-05, 0.96755E-05, 0.10510E-04, 0.10544E-04, 0.11159E-04, 
     &  0.12071E-04, 0.12377E-04, 0.16533E-04, 0.16805E-04, 0.17257E-04, 
     &  0.19604E-04, 0.22462E-04, 0.24114E-04, 0.28542E-04, 0.28873E-04, 
     &  0.33323E-04, 0.33545E-04, 0.37948E-04, 0.40965E-04, 0.41695E-04, 
     &  0.45472E-04, 0.51141E-04, 0.57196E-04, 0.60306E-04, 0.73201E-04, 
     &  0.83956E-04, 0.84494E-04, 0.11318E-03, 0.14201E-03, 0.15837E-03, 
     &  0.23860E-03, 0.37902E-03, 0.57606E-03, 0.67402E-03, 0.78450E-03, 
     &  0.90862E-03, 0.10476E-02, 0.12027E-02, 0.13753E-02, 0.15668E-02, 
     &  0.17786E-02, 0.20125E-02, 0.22700E-02, 0.25450E-02, 0.28620E-02, 
     &  0.32290E-02, 0.36570E-02, 0.41600E-02, 0.47550E-02, 0.54670E-02, 
     &  0.63280E-02, 0.73800E-02, 0.86830E-02, 0.94550E-02, 0.10320E-01, 
     &  0.11310E-01, 0.12430E-01, 0.13700E-01, 0.15170E-01, 0.16870E-01, 
     &  0.18850E-01, 0.21170E-01, 0.23910E-01, 0.27170E-01, 0.31080E-01, 
     &  0.35830E-01, 0.41630E-01, 0.48830E-01, 0.57880E-01, 0.69490E-01, 
     &  0.84870E-01, 0.10630E+00, 0.24887E+00, 0.32537E+00, 0.40042E+00, 
     &  0.46837E+00, 0.53014E+00, 0.58613E+00, 0.63592E+00, 0.68142E+00, 
     &  0.72415E+00, 0.76326E+00, 0.80161E+00, 0.84201E+00, 0.88309E+00, 
     &  0.92619E+00, 0.97485E+00, 0.10268E+01, 0.10844E+01, 0.11530E+01, 
     &  0.12330E+01, 0.13284E+01, 0.14522E+01, 0.16060E+01, 0.17911E+01, 
     &  0.20144E+01, 0.22473E+01, 0.24389E+01, 0.25562E+01, 0.25772E+01, 
     &  0.25136E+01, 0.24021E+01, 0.22524E+01, 0.20742E+01, 0.18949E+01, 
     &  0.17024E+01, 0.14957E+01, 0.12974E+01, 0.10995E+01, 0.90342E+00, 
     &  0.72745E+00, 0.56471E+00, 0.41700E+00, 0.31012E+00, 0.22741E+00, 
     &  0.16166E+00, 0.12549E+00, 0.10124E+00, 0.79340E-01, 0.70724E-01, 
     &  0.64838E-01, 0.52878E-01, 0.47674E-01, 0.46450E-01, 0.38904E-01, 
     &  0.35116E-01, 0.37200E-01, 0.32984E-01, 0.29962E-01, 0.32340E-01, 
     &  0.30172E-01, 0.27368E-01, 0.28612E-01, 0.28284E-01, 0.27712E-01, 
     &  0.30826E-01, 0.32530E-01, 0.30548E-01, 0.37498E-01, 0.48242E-01, 
     &  0.62304E-01, 0.88676E-01, 0.11917E+00, 0.14704E+00, 0.17700E+00, 
     &  0.21244E+00, 0.23983E+00, 0.26463E+00, 0.28627E+00, 0.29644E+00, 
     &  0.28964E+00, 0.27544E+00, 0.24604E+00, 0.20576E+00, 0.16296E+00, 
     &  0.12179E+00, 0.10209E+00, 0.89536E-01, 0.83262E-01, 0.75268E-01, 
     &  0.74116E-01, 0.71756E-01, 0.73404E-01, 0.73198E-01, 0.72882E-01, 
     &  0.72884E-01, 0.73498E-01, 0.75270E-01, 0.80838E-01, 0.80502E-01, 
     &  0.78668E-01, 0.80346E-01, 0.79604E-01, 0.80606E-01, 0.81386E-01, 
     &  0.84862E-01, 0.86650E-01, 0.93126E-01, 0.95176E-01, 0.10170E+00, 
     &  0.10975E+00, 0.11560E+00, 0.11861E+00, 0.12706E+00, 0.13776E+00, 
     &  0.14961E+00, 0.16178E+00, 0.17467E+00, 0.18339E+00, 0.18973E+00, 
     &  0.20719E+00, 0.22175E+00, 0.23172E+00, 0.24523E+00, 0.25735E+00, 
     &  0.26788E+00, 0.28169E+00, 0.30983E+00, 0.36783E+00, 0.45176E+00, 
     &  0.58202E+00, 0.76736E+00, 0.10108E+01, 0.13024E+01, 0.16324E+01, 
     &  0.19716E+01, 0.22228E+01, 0.22867E+01, 0.21471E+01, 0.18774E+01, 
     &  0.15797E+01, 0.12983E+01, 0.10589E+01, 0.87170E+00, 0.72331E+00, 
     &  0.60039E+00, 0.51195E+00, 0.43912E+00, 0.37780E+00, 0.33936E+00, 
     &  0.30788E+00, 0.27242E+00, 0.24929E+00, 0.23313E+00, 0.20780E+00, 
     &  0.19239E+00, 0.17958E+00, 0.15706E+00, 0.14208E+00, 0.13307E+00, 
     &  0.11907E+00, 0.11236E+00, 0.11177E+00, 0.10791E+00, 0.10074E+00, 
     &  0.10108E+00, 0.92454E-01, 0.82412E-01, 0.76994E-01, 0.67468E-01, 
     &  0.56968E-01, 0.51742E-01, 0.43990E-01, 0.38784E-01, 0.41376E-01, 
     &  0.42664E-01, 0.39344E-01, 0.36706E-01, 0.46074E-01, 0.41226E-01, 
     &  0.46460E-01, 0.36008E-01, 0.28570E-01, 0.25838E-01, 0.22268E-01, 
     &  0.15082E-01, 0.12618E-01, 0.14152E-01, 0.96960E-02, 0.11570E-01, 
     &  0.18468E-01, 0.14476E-01, 0.10788E-01, 0.13330E-01, 0.12478E-01, 
     &  0.71040E-02, 0.11900E-01, 0.14176E-01, 0.98360E-02, 0.12138E-01, 
     &  0.11296E-01, 0.95220E-02, 0.16706E-01, 0.26996E-01, 0.25082E-01, 
     &  0.17284E-01, 0.18146E-01, 0.22936E-01, 0.18238E-01, 0.12304E-01, 
     &  0.50200E-02, 0.24520E-02, 0.12740E-02, 0.75800E-03 /
	nlam=         434
	l_lnk(1:nlam)=l0(1:nlam)
	n_lnk(1:nlam)=n0(1:nlam)
	k_lnk(1:nlam)=k0(1:nlam)
	return
	end
