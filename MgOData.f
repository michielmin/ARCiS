	subroutine MgO(l_lnk,n_lnk,k_lnk,nlam)
	IMPLICIT NONE
	integer nlam,j
	real*8 l_lnk(*),n_lnk(*),k_lnk(*)
	real*8 l0(         500),n0(         500),k0(         500)
	data (l0(j),j=1,         500) /
     &  0.17000E-01, 0.17354E-01, 0.17716E-01, 0.18085E-01, 0.18462E-01, 
     &  0.18846E-01, 0.19239E-01, 0.19640E-01, 0.20049E-01, 0.20466E-01, 
     &  0.20893E-01, 0.21328E-01, 0.21773E-01, 0.22226E-01, 0.22689E-01, 
     &  0.23162E-01, 0.23644E-01, 0.24137E-01, 0.24640E-01, 0.25153E-01, 
     &  0.25677E-01, 0.26212E-01, 0.26758E-01, 0.27316E-01, 0.27885E-01, 
     &  0.28466E-01, 0.29059E-01, 0.29664E-01, 0.30282E-01, 0.30913E-01, 
     &  0.31557E-01, 0.32215E-01, 0.32886E-01, 0.33571E-01, 0.34270E-01, 
     &  0.34984E-01, 0.35713E-01, 0.36457E-01, 0.37217E-01, 0.37992E-01, 
     &  0.38784E-01, 0.39592E-01, 0.40416E-01, 0.41258E-01, 0.42118E-01, 
     &  0.42995E-01, 0.43891E-01, 0.44806E-01, 0.45739E-01, 0.46692E-01, 
     &  0.47665E-01, 0.48658E-01, 0.49671E-01, 0.50706E-01, 0.51763E-01, 
     &  0.52841E-01, 0.53942E-01, 0.55066E-01, 0.56213E-01, 0.57384E-01, 
     &  0.58580E-01, 0.59800E-01, 0.61046E-01, 0.62318E-01, 0.63616E-01, 
     &  0.64941E-01, 0.66294E-01, 0.67675E-01, 0.69085E-01, 0.70525E-01, 
     &  0.71994E-01, 0.73494E-01, 0.75025E-01, 0.76588E-01, 0.78184E-01, 
     &  0.79812E-01, 0.81475E-01, 0.83173E-01, 0.84905E-01, 0.86674E-01, 
     &  0.88480E-01, 0.90323E-01, 0.92205E-01, 0.94126E-01, 0.96087E-01, 
     &  0.98089E-01, 0.10013E+00, 0.10222E+00, 0.10435E+00, 0.10652E+00, 
     &  0.10874E+00, 0.11101E+00, 0.11332E+00, 0.11568E+00, 0.11809E+00, 
     &  0.12055E+00, 0.12306E+00, 0.12563E+00, 0.12824E+00, 0.13091E+00, 
     &  0.13364E+00, 0.13643E+00, 0.13927E+00, 0.14217E+00, 0.14513E+00, 
     &  0.14816E+00, 0.15124E+00, 0.15439E+00, 0.15761E+00, 0.16089E+00, 
     &  0.16425E+00, 0.16767E+00, 0.17116E+00, 0.17473E+00, 0.17837E+00, 
     &  0.18208E+00, 0.18588E+00, 0.18975E+00, 0.19370E+00, 0.19774E+00, 
     &  0.20186E+00, 0.20606E+00, 0.21036E+00, 0.21474E+00, 0.21921E+00, 
     &  0.22378E+00, 0.22844E+00, 0.23320E+00, 0.23806E+00, 0.24302E+00, 
     &  0.24808E+00, 0.25325E+00, 0.25853E+00, 0.26391E+00, 0.26941E+00, 
     &  0.27502E+00, 0.28075E+00, 0.28660E+00, 0.29257E+00, 0.29867E+00, 
     &  0.30489E+00, 0.31124E+00, 0.31773E+00, 0.32434E+00, 0.33110E+00, 
     &  0.33800E+00, 0.34504E+00, 0.35223E+00, 0.35957E+00, 0.36706E+00, 
     &  0.37471E+00, 0.38251E+00, 0.39048E+00, 0.39862E+00, 0.40692E+00, 
     &  0.41540E+00, 0.42405E+00, 0.43289E+00, 0.44191E+00, 0.45111E+00, 
     &  0.46051E+00, 0.47011E+00, 0.47990E+00, 0.48990E+00, 0.50010E+00, 
     &  0.51052E+00, 0.52116E+00, 0.53202E+00, 0.54310E+00, 0.55442E+00, 
     &  0.56597E+00, 0.57776E+00, 0.58979E+00, 0.60208E+00, 0.61463E+00, 
     &  0.62743E+00, 0.64050E+00, 0.65385E+00, 0.66747E+00, 0.68137E+00, 
     &  0.69557E+00, 0.71006E+00, 0.72485E+00, 0.73995E+00, 0.75537E+00, 
     &  0.77111E+00, 0.78717E+00, 0.80357E+00, 0.82031E+00, 0.83740E+00, 
     &  0.85485E+00, 0.87266E+00, 0.89084E+00, 0.90940E+00, 0.92835E+00, 
     &  0.94769E+00, 0.96743E+00, 0.98759E+00, 0.10082E+01, 0.10292E+01, 
     &  0.10506E+01, 0.10725E+01, 0.10948E+01, 0.11176E+01, 0.11409E+01, 
     &  0.11647E+01, 0.11890E+01, 0.12137E+01, 0.12390E+01, 0.12648E+01, 
     &  0.12912E+01, 0.13181E+01, 0.13455E+01, 0.13736E+01, 0.14022E+01, 
     &  0.14314E+01, 0.14612E+01, 0.14917E+01, 0.15227E+01, 0.15545E+01, 
     &  0.15869E+01, 0.16199E+01, 0.16537E+01, 0.16881E+01, 0.17233E+01, 
     &  0.17592E+01, 0.17958E+01, 0.18333E+01, 0.18714E+01, 0.19104E+01, 
     &  0.19502E+01, 0.19909E+01, 0.20323E+01, 0.20747E+01, 0.21179E+01, 
     &  0.21620E+01, 0.22071E+01, 0.22531E+01, 0.23000E+01, 0.23479E+01, 
     &  0.23968E+01, 0.24468E+01, 0.24977E+01, 0.25498E+01, 0.26029E+01, 
     &  0.26571E+01, 0.27125E+01, 0.27690E+01, 0.28267E+01, 0.28856E+01, 
     &  0.29457E+01, 0.30071E+01, 0.30697E+01, 0.31337E+01, 0.31989E+01, 
     &  0.32656E+01, 0.33336E+01, 0.34031E+01, 0.34740E+01, 0.35463E+01, 
     &  0.36202E+01, 0.36957E+01, 0.37726E+01, 0.38512E+01, 0.39315E+01, 
     &  0.40134E+01, 0.40970E+01, 0.41824E+01, 0.42695E+01, 0.43584E+01, 
     &  0.44492E+01, 0.45419E+01, 0.46366E+01, 0.47332E+01, 0.48318E+01, 
     &  0.49324E+01, 0.50352E+01, 0.51401E+01, 0.52472E+01, 0.53565E+01, 
     &  0.54681E+01, 0.55820E+01, 0.56983E+01, 0.58170E+01, 0.59382E+01, 
     &  0.60619E+01, 0.61882E+01, 0.63171E+01, 0.64487E+01, 0.65831E+01, 
     &  0.67202E+01, 0.68603E+01, 0.70032E+01, 0.71491E+01, 0.72980E+01, 
     &  0.74501E+01, 0.76053E+01, 0.77637E+01, 0.79255E+01, 0.80906E+01, 
     &  0.82591E+01, 0.84312E+01, 0.86069E+01, 0.87862E+01, 0.89692E+01, 
     &  0.91561E+01, 0.93468E+01, 0.95416E+01, 0.97403E+01, 0.99433E+01, 
     &  0.10150E+02, 0.10362E+02, 0.10578E+02, 0.10798E+02, 0.11023E+02, 
     &  0.11253E+02, 0.11487E+02, 0.11727E+02, 0.11971E+02, 0.12220E+02, 
     &  0.12475E+02, 0.12735E+02, 0.13000E+02, 0.13271E+02, 0.13547E+02, 
     &  0.13830E+02, 0.14118E+02, 0.14412E+02, 0.14712E+02, 0.15019E+02, 
     &  0.15331E+02, 0.15651E+02, 0.15977E+02, 0.16310E+02, 0.16650E+02, 
     &  0.16996E+02, 0.17351E+02, 0.17712E+02, 0.18081E+02, 0.18458E+02, 
     &  0.18842E+02, 0.19235E+02, 0.19636E+02, 0.20045E+02, 0.20462E+02, 
     &  0.20889E+02, 0.21324E+02, 0.21768E+02, 0.22221E+02, 0.22684E+02, 
     &  0.23157E+02, 0.23639E+02, 0.24132E+02, 0.24635E+02, 0.25148E+02, 
     &  0.25672E+02, 0.26207E+02, 0.26753E+02, 0.27310E+02, 0.27879E+02, 
     &  0.28460E+02, 0.29053E+02, 0.29658E+02, 0.30276E+02, 0.30907E+02, 
     &  0.31551E+02, 0.32208E+02, 0.32879E+02, 0.33564E+02, 0.34263E+02, 
     &  0.34977E+02, 0.35706E+02, 0.36449E+02, 0.37209E+02, 0.37984E+02, 
     &  0.38775E+02, 0.39583E+02, 0.40408E+02, 0.41250E+02, 0.42109E+02, 
     &  0.42986E+02, 0.43882E+02, 0.44796E+02, 0.45729E+02, 0.46682E+02, 
     &  0.47655E+02, 0.48647E+02, 0.49661E+02, 0.50696E+02, 0.51752E+02, 
     &  0.52830E+02, 0.53931E+02, 0.55054E+02, 0.56201E+02, 0.57372E+02, 
     &  0.58567E+02, 0.59787E+02, 0.61033E+02, 0.62305E+02, 0.63603E+02, 
     &  0.64928E+02, 0.66280E+02, 0.67661E+02, 0.69071E+02, 0.70510E+02, 
     &  0.71979E+02, 0.73478E+02, 0.75009E+02, 0.76572E+02, 0.78167E+02, 
     &  0.79796E+02, 0.81458E+02, 0.83155E+02, 0.84888E+02, 0.86656E+02, 
     &  0.88461E+02, 0.90304E+02, 0.92186E+02, 0.94106E+02, 0.96067E+02, 
     &  0.98068E+02, 0.10011E+03, 0.10220E+03, 0.10433E+03, 0.10650E+03, 
     &  0.10872E+03, 0.11098E+03, 0.11330E+03, 0.11566E+03, 0.11807E+03, 
     &  0.12053E+03, 0.12304E+03, 0.12560E+03, 0.12822E+03, 0.13089E+03, 
     &  0.13361E+03, 0.13640E+03, 0.13924E+03, 0.14214E+03, 0.14510E+03, 
     &  0.14812E+03, 0.15121E+03, 0.15436E+03, 0.15758E+03, 0.16086E+03, 
     &  0.16421E+03, 0.16763E+03, 0.17112E+03, 0.17469E+03, 0.17833E+03, 
     &  0.18204E+03, 0.18584E+03, 0.18971E+03, 0.19366E+03, 0.19770E+03, 
     &  0.20181E+03, 0.20602E+03, 0.21031E+03, 0.21469E+03, 0.21917E+03, 
     &  0.22373E+03, 0.22839E+03, 0.23315E+03, 0.23801E+03, 0.24297E+03, 
     &  0.24803E+03, 0.25320E+03, 0.25847E+03, 0.26386E+03, 0.26935E+03, 
     &  0.27496E+03, 0.28069E+03, 0.28654E+03, 0.29251E+03, 0.29860E+03, 
     &  0.30483E+03, 0.31118E+03, 0.31766E+03, 0.32428E+03, 0.33103E+03, 
     &  0.33793E+03, 0.34497E+03, 0.35216E+03, 0.35949E+03, 0.36698E+03, 
     &  0.37463E+03, 0.38243E+03, 0.39040E+03, 0.39853E+03, 0.40684E+03, 
     &  0.41531E+03, 0.42397E+03, 0.43280E+03, 0.44181E+03, 0.45102E+03, 
     &  0.46042E+03, 0.47001E+03, 0.47980E+03, 0.48980E+03, 0.50000E+03 /
	data (n0(j),j=1,         500) /
     &  0.97704E+00, 0.97792E+00, 0.97881E+00, 0.97970E+00, 0.98059E+00, 
     &  0.98147E+00, 0.98236E+00, 0.98326E+00, 0.98415E+00, 0.98504E+00, 
     &  0.98593E+00, 0.98683E+00, 0.98772E+00, 0.98575E+00, 0.98190E+00, 
     &  0.97806E+00, 0.97423E+00, 0.97042E+00, 0.96663E+00, 0.96285E+00, 
     &  0.95908E+00, 0.95533E+00, 0.95160E+00, 0.94787E+00, 0.94417E+00, 
     &  0.94048E+00, 0.93680E+00, 0.93329E+00, 0.93278E+00, 0.93227E+00, 
     &  0.93177E+00, 0.93126E+00, 0.93075E+00, 0.93024E+00, 0.92974E+00, 
     &  0.92923E+00, 0.92873E+00, 0.92822E+00, 0.92771E+00, 0.92721E+00, 
     &  0.91412E+00, 0.89390E+00, 0.87412E+00, 0.85477E+00, 0.83586E+00, 
     &  0.81736E+00, 0.79928E+00, 0.78159E+00, 0.76354E+00, 0.74586E+00, 
     &  0.72859E+00, 0.71172E+00, 0.69524E+00, 0.67914E+00, 0.66342E+00, 
     &  0.64621E+00, 0.62446E+00, 0.60345E+00, 0.58314E+00, 0.56352E+00, 
     &  0.59685E+00, 0.63547E+00, 0.67660E+00, 0.72038E+00, 0.80890E+00, 
     &  0.92946E+00, 0.96061E+00, 0.99282E+00, 0.10261E+01, 0.10605E+01, 
     &  0.11705E+01, 0.12925E+01, 0.12510E+01, 0.12108E+01, 0.11719E+01, 
     &  0.11343E+01, 0.10999E+01, 0.10707E+01, 0.10423E+01, 0.10147E+01, 
     &  0.11038E+01, 0.13698E+01, 0.14967E+01, 0.16086E+01, 0.17030E+01, 
     &  0.17326E+01, 0.15056E+01, 0.14796E+01, 0.14547E+01, 0.14301E+01, 
     &  0.14099E+01, 0.15585E+01, 0.18231E+01, 0.19841E+01, 0.20095E+01, 
     &  0.20352E+01, 0.20612E+01, 0.20875E+01, 0.21142E+01, 0.21413E+01, 
     &  0.21686E+01, 0.21754E+01, 0.21810E+01, 0.21866E+01, 0.21922E+01, 
     &  0.21979E+01, 0.22036E+01, 0.22092E+01, 0.22985E+01, 0.27390E+01, 
     &  0.30391E+01, 0.25648E+01, 0.23479E+01, 0.22968E+01, 0.22468E+01, 
     &  0.21978E+01, 0.21500E+01, 0.21185E+01, 0.20991E+01, 0.20800E+01, 
     &  0.20610E+01, 0.20421E+01, 0.20235E+01, 0.20050E+01, 0.19867E+01, 
     &  0.19685E+01, 0.19562E+01, 0.19477E+01, 0.19392E+01, 0.19307E+01, 
     &  0.19222E+01, 0.19139E+01, 0.19055E+01, 0.18972E+01, 0.18889E+01, 
     &  0.18806E+01, 0.18724E+01, 0.18642E+01, 0.18561E+01, 0.18480E+01, 
     &  0.18399E+01, 0.18327E+01, 0.18330E+01, 0.18332E+01, 0.18334E+01, 
     &  0.18336E+01, 0.18338E+01, 0.18340E+01, 0.18342E+01, 0.18344E+01, 
     &  0.18346E+01, 0.18348E+01, 0.18350E+01, 0.18352E+01, 0.18354E+01, 
     &  0.18356E+01, 0.18332E+01, 0.18291E+01, 0.18251E+01, 0.18210E+01, 
     &  0.18170E+01, 0.18129E+01, 0.18089E+01, 0.18049E+01, 0.18009E+01, 
     &  0.17969E+01, 0.17929E+01, 0.17889E+01, 0.17849E+01, 0.17810E+01, 
     &  0.17798E+01, 0.17786E+01, 0.17773E+01, 0.17761E+01, 0.17749E+01, 
     &  0.17736E+01, 0.17724E+01, 0.17712E+01, 0.17699E+01, 0.17687E+01, 
     &  0.17675E+01, 0.17662E+01, 0.17650E+01, 0.17638E+01, 0.17625E+01, 
     &  0.17611E+01, 0.17596E+01, 0.17582E+01, 0.17567E+01, 0.17553E+01, 
     &  0.17538E+01, 0.17524E+01, 0.17510E+01, 0.17495E+01, 0.17481E+01, 
     &  0.17467E+01, 0.17452E+01, 0.17438E+01, 0.17430E+01, 0.17429E+01, 
     &  0.17429E+01, 0.17428E+01, 0.17427E+01, 0.17427E+01, 0.17426E+01, 
     &  0.17425E+01, 0.17424E+01, 0.17424E+01, 0.17423E+01, 0.17422E+01, 
     &  0.17421E+01, 0.17421E+01, 0.17420E+01, 0.17419E+01, 0.17419E+01, 
     &  0.17423E+01, 0.17444E+01, 0.17465E+01, 0.17486E+01, 0.17506E+01, 
     &  0.17527E+01, 0.17548E+01, 0.17569E+01, 0.17590E+01, 0.17611E+01, 
     &  0.17632E+01, 0.17653E+01, 0.17675E+01, 0.17696E+01, 0.17717E+01, 
     &  0.17738E+01, 0.17759E+01, 0.17757E+01, 0.17726E+01, 0.17696E+01, 
     &  0.17665E+01, 0.17634E+01, 0.17603E+01, 0.17572E+01, 0.17542E+01, 
     &  0.17511E+01, 0.17480E+01, 0.17450E+01, 0.17419E+01, 0.17389E+01, 
     &  0.17359E+01, 0.17328E+01, 0.17298E+01, 0.17268E+01, 0.17238E+01, 
     &  0.17208E+01, 0.17213E+01, 0.17220E+01, 0.17226E+01, 0.17232E+01, 
     &  0.17238E+01, 0.17245E+01, 0.17251E+01, 0.17257E+01, 0.17264E+01, 
     &  0.17270E+01, 0.17270E+01, 0.17241E+01, 0.17211E+01, 0.17182E+01, 
     &  0.17152E+01, 0.17123E+01, 0.17093E+01, 0.17064E+01, 0.17035E+01, 
     &  0.17006E+01, 0.16976E+01, 0.16947E+01, 0.16918E+01, 0.16889E+01, 
     &  0.16845E+01, 0.16797E+01, 0.16749E+01, 0.16700E+01, 0.16652E+01, 
     &  0.16604E+01, 0.16557E+01, 0.16509E+01, 0.16462E+01, 0.16414E+01, 
     &  0.16348E+01, 0.16230E+01, 0.16113E+01, 0.15997E+01, 0.15881E+01, 
     &  0.15767E+01, 0.15653E+01, 0.15540E+01, 0.15428E+01, 0.15316E+01, 
     &  0.15206E+01, 0.15070E+01, 0.14923E+01, 0.14777E+01, 0.14633E+01, 
     &  0.14490E+01, 0.14349E+01, 0.14209E+01, 0.14070E+01, 0.13805E+01, 
     &  0.13502E+01, 0.13206E+01, 0.12916E+01, 0.12633E+01, 0.12356E+01, 
     &  0.12085E+01, 0.11719E+01, 0.11127E+01, 0.10565E+01, 0.10031E+01, 
     &  0.95237E+00, 0.90425E+00, 0.85386E+00, 0.78376E+00, 0.71941E+00, 
     &  0.66035E+00, 0.53575E+00, 0.42196E+00, 0.21380E+00, 0.15751E+00, 
     &  0.14666E+00, 0.15088E+00, 0.15523E+00, 0.15135E+00, 0.13840E+00, 
     &  0.12217E+00, 0.10448E+00, 0.92709E-01, 0.90513E-01, 0.91190E-01, 
     &  0.91871E-01, 0.92558E-01, 0.93250E-01, 0.96380E-01, 0.10025E+00, 
     &  0.10428E+00, 0.11134E+00, 0.12215E+00, 0.13506E+00, 0.14997E+00, 
     &  0.17472E+00, 0.20777E+00, 0.27081E+00, 0.37355E+00, 0.99919E+00, 
     &  0.18908E+01, 0.77600E+01, 0.12169E+02, 0.96006E+01, 0.81113E+01, 
     &  0.71202E+01, 0.66406E+01, 0.62573E+01, 0.58961E+01, 0.55677E+01, 
     &  0.53128E+01, 0.50697E+01, 0.48376E+01, 0.46162E+01, 0.44802E+01, 
     &  0.43717E+01, 0.42659E+01, 0.41626E+01, 0.40618E+01, 0.40115E+01, 
     &  0.39676E+01, 0.39242E+01, 0.38812E+01, 0.38388E+01, 0.37967E+01, 
     &  0.37552E+01, 0.37283E+01, 0.37044E+01, 0.36808E+01, 0.36572E+01, 
     &  0.36339E+01, 0.36106E+01, 0.35876E+01, 0.35646E+01, 0.35451E+01, 
     &  0.35299E+01, 0.35148E+01, 0.34998E+01, 0.34848E+01, 0.34699E+01, 
     &  0.34551E+01, 0.34403E+01, 0.34256E+01, 0.34109E+01, 0.33963E+01, 
     &  0.33865E+01, 0.33784E+01, 0.33702E+01, 0.33622E+01, 0.33541E+01, 
     &  0.33460E+01, 0.33380E+01, 0.33300E+01, 0.33220E+01, 0.33140E+01, 
     &  0.33060E+01, 0.32981E+01, 0.32904E+01, 0.32830E+01, 0.32757E+01, 
     &  0.32684E+01, 0.32611E+01, 0.32538E+01, 0.32466E+01, 0.32393E+01, 
     &  0.32321E+01, 0.32249E+01, 0.32177E+01, 0.32105E+01, 0.32034E+01, 
     &  0.31965E+01, 0.31980E+01, 0.31994E+01, 0.32009E+01, 0.32023E+01, 
     &  0.32038E+01, 0.32053E+01, 0.32067E+01, 0.32082E+01, 0.32096E+01, 
     &  0.32111E+01, 0.32126E+01, 0.32140E+01, 0.32155E+01, 0.32170E+01, 
     &  0.32184E+01, 0.32199E+01, 0.32205E+01, 0.32196E+01, 0.32187E+01, 
     &  0.32178E+01, 0.32169E+01, 0.32161E+01, 0.32152E+01, 0.32143E+01, 
     &  0.32134E+01, 0.32125E+01, 0.32116E+01, 0.32107E+01, 0.32098E+01, 
     &  0.32090E+01, 0.32081E+01, 0.32104E+01, 0.32152E+01, 0.32201E+01, 
     &  0.32249E+01, 0.32298E+01, 0.32347E+01, 0.32395E+01, 0.32444E+01, 
     &  0.32493E+01, 0.32542E+01, 0.32513E+01, 0.32384E+01, 0.32256E+01, 
     &  0.32128E+01, 0.32001E+01, 0.31910E+01, 0.32036E+01, 0.32162E+01, 
     &  0.32289E+01, 0.32416E+01, 0.32544E+01, 0.32673E+01, 0.32802E+01, 
     &  0.32685E+01, 0.32557E+01, 0.32430E+01, 0.32303E+01, 0.32177E+01, 
     &  0.32051E+01, 0.31926E+01, 0.31801E+01, 0.31677E+01, 0.31663E+01, 
     &  0.31653E+01, 0.31643E+01, 0.31633E+01, 0.31623E+01, 0.31613E+01, 
     &  0.31603E+01, 0.31593E+01, 0.31606E+01, 0.31623E+01, 0.31639E+01, 
     &  0.31655E+01, 0.31671E+01, 0.31687E+01, 0.31703E+01, 0.31720E+01 /
	data (k0(j),j=1,         500) /
     &  0.71834E-01, 0.73224E-01, 0.74640E-01, 0.76084E-01, 0.74265E-01, 
     &  0.72171E-01, 0.70135E-01, 0.68157E-01, 0.68330E-01, 0.81208E-01, 
     &  0.97592E-01, 0.11345E+00, 0.11101E+00, 0.10862E+00, 0.10038E+00, 
     &  0.88247E-01, 0.85450E-01, 0.86808E-01, 0.88188E-01, 0.89590E-01, 
     &  0.91779E-01, 0.94495E-01, 0.97292E-01, 0.10052E+00, 0.10483E+00, 
     &  0.10933E+00, 0.11402E+00, 0.11892E+00, 0.12348E+00, 0.12743E+00, 
     &  0.13149E+00, 0.13569E+00, 0.14002E+00, 0.14449E+00, 0.15025E+00, 
     &  0.15625E+00, 0.16249E+00, 0.16898E+00, 0.16818E+00, 0.16466E+00, 
     &  0.16122E+00, 0.15786E+00, 0.15456E+00, 0.15133E+00, 0.14838E+00, 
     &  0.15404E+00, 0.15993E+00, 0.16764E+00, 0.19041E+00, 0.21565E+00, 
     &  0.23943E+00, 0.26582E+00, 0.29513E+00, 0.32729E+00, 0.35739E+00, 
     &  0.39025E+00, 0.42614E+00, 0.47435E+00, 0.54489E+00, 0.62280E+00, 
     &  0.68632E+00, 0.75633E+00, 0.82312E+00, 0.86274E+00, 0.90427E+00, 
     &  0.94780E+00, 0.97667E+00, 0.10044E+01, 0.98535E+00, 0.91847E+00, 
     &  0.85583E+00, 0.79067E+00, 0.78170E+00, 0.81070E+00, 0.84078E+00, 
     &  0.86800E+00, 0.88325E+00, 0.89876E+00, 0.91455E+00, 0.93061E+00, 
     &  0.94696E+00, 0.10459E+01, 0.11708E+01, 0.12200E+01, 0.10529E+01, 
     &  0.93634E+00, 0.88198E+00, 0.89345E+00, 0.10273E+01, 0.11152E+01, 
     &  0.11671E+01, 0.12214E+01, 0.11724E+01, 0.10720E+01, 0.98713E+00, 
     &  0.91769E+00, 0.85313E+00, 0.79311E+00, 0.74736E+00, 0.71665E+00, 
     &  0.68720E+00, 0.65897E+00, 0.63569E+00, 0.61387E+00, 0.59280E+00, 
     &  0.58345E+00, 0.64665E+00, 0.73144E+00, 0.88772E+00, 0.59235E+00, 
     &  0.21665E-01, 0.30965E-03, 0.49648E-04, 0.37587E-04, 0.28456E-04, 
     &  0.22768E-04, 0.18976E-04, 0.15815E-04, 0.13181E-04, 0.11204E-04, 
     &  0.99021E-05, 0.87511E-05, 0.77340E-05, 0.68350E-05, 0.60406E-05, 
     &  0.53385E-05, 0.47180E-05, 0.41696E-05, 0.36850E-05, 0.32567E-05, 
     &  0.29654E-05, 0.27171E-05, 0.24897E-05, 0.22812E-05, 0.20903E-05, 
     &  0.19153E-05, 0.17550E-05, 0.16080E-05, 0.14734E-05, 0.13501E-05, 
     &  0.12371E-05, 0.11335E-05, 0.10386E-05, 0.95166E-06, 0.87199E-06, 
     &  0.79899E-06, 0.66834E-06, 0.55108E-06, 0.45440E-06, 0.37468E-06, 
     &  0.29458E-06, 0.22963E-06, 0.20299E-06, 0.19102E-06, 0.17976E-06, 
     &  0.16916E-06, 0.15918E-06, 0.14980E-06, 0.14097E-06, 0.13266E-06, 
     &  0.12445E-06, 0.10722E-06, 0.92371E-07, 0.79580E-07, 0.68560E-07, 
     &  0.62539E-07, 0.61805E-07, 0.61080E-07, 0.60363E-07, 0.59655E-07, 
     &  0.58955E-07, 0.58263E-07, 0.57580E-07, 0.56904E-07, 0.56236E-07, 
     &  0.55577E-07, 0.54925E-07, 0.54280E-07, 0.53643E-07, 0.53014E-07, 
     &  0.52392E-07, 0.51777E-07, 0.51170E-07, 0.50569E-07, 0.49976E-07, 
     &  0.49390E-07, 0.48810E-07, 0.48237E-07, 0.47671E-07, 0.47112E-07, 
     &  0.46559E-07, 0.46013E-07, 0.45473E-07, 0.44940E-07, 0.44412E-07, 
     &  0.43891E-07, 0.43376E-07, 0.42867E-07, 0.42364E-07, 0.41867E-07, 
     &  0.41376E-07, 0.40891E-07, 0.40411E-07, 0.39937E-07, 0.39468E-07, 
     &  0.39005E-07, 0.38547E-07, 0.38095E-07, 0.37648E-07, 0.37207E-07, 
     &  0.36770E-07, 0.36339E-07, 0.35912E-07, 0.35491E-07, 0.35074E-07, 
     &  0.34663E-07, 0.34256E-07, 0.33854E-07, 0.33457E-07, 0.33065E-07, 
     &  0.32677E-07, 0.32293E-07, 0.31914E-07, 0.31540E-07, 0.31170E-07, 
     &  0.30804E-07, 0.30443E-07, 0.30085E-07, 0.29732E-07, 0.29384E-07, 
     &  0.29039E-07, 0.28698E-07, 0.28361E-07, 0.28029E-07, 0.27700E-07, 
     &  0.27375E-07, 0.27054E-07, 0.26736E-07, 0.26423E-07, 0.26112E-07, 
     &  0.25806E-07, 0.25503E-07, 0.25204E-07, 0.24908E-07, 0.24616E-07, 
     &  0.24327E-07, 0.24042E-07, 0.23760E-07, 0.23481E-07, 0.23206E-07, 
     &  0.22933E-07, 0.22664E-07, 0.22398E-07, 0.22135E-07, 0.21876E-07, 
     &  0.21619E-07, 0.21365E-07, 0.21115E-07, 0.20867E-07, 0.20622E-07, 
     &  0.20380E-07, 0.24967E-07, 0.32809E-07, 0.43813E-07, 0.58507E-07, 
     &  0.78128E-07, 0.10433E-06, 0.13932E-06, 0.18604E-06, 0.24844E-06, 
     &  0.33322E-06, 0.45197E-06, 0.61303E-06, 0.83148E-06, 0.11278E-05, 
     &  0.15297E-05, 0.20748E-05, 0.28141E-05, 0.36682E-05, 0.46738E-05, 
     &  0.59551E-05, 0.75876E-05, 0.96677E-05, 0.12318E-04, 0.15695E-04, 
     &  0.19998E-04, 0.25480E-04, 0.32465E-04, 0.41365E-04, 0.51529E-04, 
     &  0.63363E-04, 0.77916E-04, 0.95811E-04, 0.11782E-03, 0.14487E-03, 
     &  0.17280E-03, 0.20183E-03, 0.23573E-03, 0.27532E-03, 0.32156E-03, 
     &  0.37558E-03, 0.48132E-03, 0.65170E-03, 0.88240E-03, 0.12381E-02, 
     &  0.21617E-02, 0.37744E-02, 0.44439E-02, 0.49483E-02, 0.55099E-02, 
     &  0.61353E-02, 0.70023E-02, 0.84804E-02, 0.10271E-01, 0.12439E-01, 
     &  0.15064E-01, 0.19473E-01, 0.26488E-01, 0.36029E-01, 0.49008E-01, 
     &  0.66661E-01, 0.95706E-01, 0.19926E+00, 0.32280E+00, 0.43192E+00, 
     &  0.57793E+00, 0.77329E+00, 0.98711E+00, 0.10815E+01, 0.11848E+01, 
     &  0.12981E+01, 0.14222E+01, 0.15524E+01, 0.16580E+01, 0.17708E+01, 
     &  0.18912E+01, 0.20198E+01, 0.21572E+01, 0.23039E+01, 0.24606E+01, 
     &  0.26279E+01, 0.28769E+01, 0.31753E+01, 0.35046E+01, 0.38681E+01, 
     &  0.42712E+01, 0.49044E+01, 0.56315E+01, 0.64663E+01, 0.76933E+01, 
     &  0.92416E+01, 0.11535E+02, 0.97581E+01, 0.49050E+01, 0.12710E+01, 
     &  0.72422E+00, 0.48267E+00, 0.35367E+00, 0.29306E+00, 0.24283E+00, 
     &  0.20122E+00, 0.16673E+00, 0.13816E+00, 0.12910E+00, 0.12158E+00, 
     &  0.11360E+00, 0.87758E-01, 0.67795E-01, 0.52373E-01, 0.40459E-01, 
     &  0.31256E-01, 0.28534E-01, 0.27830E-01, 0.27144E-01, 0.26474E-01, 
     &  0.25821E-01, 0.25184E-01, 0.24563E-01, 0.23957E-01, 0.23366E-01, 
     &  0.22790E-01, 0.22228E-01, 0.22215E-01, 0.22543E-01, 0.22876E-01, 
     &  0.23214E-01, 0.23558E-01, 0.23906E-01, 0.24259E-01, 0.24617E-01, 
     &  0.24981E-01, 0.25350E-01, 0.25726E-01, 0.26111E-01, 0.26501E-01, 
     &  0.26896E-01, 0.27298E-01, 0.27706E-01, 0.28120E-01, 0.28540E-01, 
     &  0.28966E-01, 0.29399E-01, 0.29838E-01, 0.30284E-01, 0.30736E-01, 
     &  0.31196E-01, 0.31662E-01, 0.32135E-01, 0.33201E-01, 0.34326E-01, 
     &  0.35489E-01, 0.36692E-01, 0.37935E-01, 0.39220E-01, 0.40549E-01, 
     &  0.41923E-01, 0.43344E-01, 0.41230E-01, 0.34854E-01, 0.29464E-01, 
     &  0.24907E-01, 0.21055E-01, 0.16458E-01, 0.12105E-01, 0.95365E-02, 
     &  0.77531E-02, 0.63032E-02, 0.52901E-02, 0.50495E-02, 0.48199E-02, 
     &  0.46007E-02, 0.43915E-02, 0.41918E-02, 0.40012E-02, 0.38192E-02, 
     &  0.36455E-02, 0.34798E-02, 0.33215E-02, 0.31705E-02, 0.30263E-02, 
     &  0.29352E-02, 0.28604E-02, 0.27874E-02, 0.27163E-02, 0.26470E-02, 
     &  0.25795E-02, 0.25137E-02, 0.24495E-02, 0.23870E-02, 0.23261E-02, 
     &  0.22668E-02, 0.22090E-02, 0.21526E-02, 0.20977E-02, 0.20444E-02, 
     &  0.19932E-02, 0.19432E-02, 0.18945E-02, 0.18471E-02, 0.18008E-02, 
     &  0.17556E-02, 0.17117E-02, 0.16688E-02, 0.16269E-02, 0.15862E-02, 
     &  0.15464E-02, 0.15077E-02, 0.14713E-02, 0.14367E-02, 0.14029E-02, 
     &  0.13699E-02, 0.13377E-02, 0.13062E-02, 0.12755E-02, 0.12455E-02, 
     &  0.12162E-02, 0.11875E-02, 0.11596E-02, 0.11323E-02, 0.11057E-02, 
     &  0.10797E-02, 0.10543E-02, 0.10279E-02, 0.10014E-02, 0.97556E-03, 
     &  0.95038E-03, 0.92585E-03, 0.90196E-03, 0.87868E-03, 0.85600E-03, 
     &  0.83391E-03, 0.81239E-03, 0.79142E-03, 0.77100E-03, 0.75110E-03, 
     &  0.73171E-03, 0.71283E-03, 0.69443E-03, 0.67651E-03, 0.65905E-03 /
	nlam=         500
	l_lnk(1:nlam)=l0(1:nlam)
	n_lnk(1:nlam)=n0(1:nlam)
	k_lnk(1:nlam)=k0(1:nlam)
	return
	end
