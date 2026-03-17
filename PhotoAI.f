	subroutine PhotoAI()
	use GlobalSetup
	use AtomsModule
	use Constants
	IMPLICIT NONE
	integer nPvulc,i,j,k,nmol_ex
	parameter(nmol_ex=11)
	parameter(nPvulc=150)
	real*8 x(10),N_H,n_H_Solar,Z,Pvulc(nPvulc),y(nPvulc*nmol_ex)
	real*8 C_O,S_O,S_O_Solar,C_O_Solar
	real*8 mixrat_ex(nPvulc)
	integer nx,ny,init,imol
	integer imol_ex(nmol_ex)
	character*4 cmol_ex(nmol_ex)
	parameter(cmol_ex = (/ "H2O ","CO  ","CO2 ","CH4 ","NH3 ","H2S ","HCN ","C2H2","NS  ","OCS ","SO2 " /))
	character*500 command

	do i=1,nmol_ex
		do j=1,nmol
			if(trim(molname(j)).eq.trim(cmol_ex(i))) then
				imol_ex(i)=j
				exit
			endif
		enddo
	enddo

	command=' '
	init=1
	
	x(1)=log10(Mplanet)
	x(2)=loggPlanet
	x(3)=log10(sqrt(Rstar/(2d0*Dplanet))*Tstar)
	x(4)=log10(Dplanet/AU)
	Z=10d0**metallicity
	x(5)=Z
	C_O=molfracs_atoms(3)/molfracs_atoms(5)
	C_O_Solar=0.0002478241/0.0004509658
	S_O=molfracs_atoms(11)/molfracs_atoms(5)
	S_O_Solar=1.2137900734604e-05/0.0004509658
	x(6)=C_O/C_O_Solar	! C/O compared to Solar
	N_H=molfracs_atoms(4)/molfracs_atoms(1)
	N_H_Solar=6.22506056949881e-05/0.9207539305
	x(7)=(N_H/N_H_Solar-1d0)/(Z-1d0)
	x(8)=S_O/S_O_Solar	! C/O compared to Solar
	x(9)=Tstar
	x(10)=Rstar/Rsun

	nx=10
	ny=nPvulc*nmol_ex
	call ask_python(x,nx,y,ny,init,command)

	do i=1,nPvulc
		Pvulc(i)=10d0**(3d0-11d0*real(i-1)/real(nPvulc-1))
	enddo

	k=0
	do i=1,nmol_ex
		do j=1,nPvulc
			k=k+1
			mixrat_ex(j)=y(k)
		enddo
		call regridarray(-log(Pvulc(1:nPvulc)),log(mixrat_ex(1:nPvulc)),nPvulc,-log(P(1:nr)),mixrat_r(1:nr,imol_ex(i)),nr)
		mixrat_r(1:nr,imol_ex(i))=exp(mixrat_r(1:nr,imol_ex(i)))
	enddo

	if(useEOS) then
		do i=1,nr
			call GetNablaEOS(P(i),T(i),molfracs_atoms(1)/(molfracs_atoms(1)+molfracs_atoms(2)),nabla_ad(i))
		enddo
	endif

	do i=1,nr
		MMW(i)=0d0
		do imol=1,nmol
			if(includemol(imol)) MMW(i)=MMW(i)+mixrat_r(i,imol)*Mmol(imol)
		enddo
	enddo

	return
	end
	
	
	subroutine InitPhotoAI()
	use GlobalSetup
	use AtomsModule
	IMPLICIT NONE
	real*8 x(10),y(10)
	integer nx,ny,init
	character*500 command,homedir

	call getenv('HOME',homedir) 
	command='python ' // trim(homedir) // '/ARCiS/Data/PhotoAI/PhotoAI.py'
	
	nx=10
	ny=10
	init=0
	call ask_python(x,nx,y,ny,init,command)

	return
	end

	subroutine ClosePhotoAI()
	IMPLICIT NONE
	real*8 x(10),y(10)
	integer nx,ny,init
	character*500 command,homedir

	command=' '
	
	nx=10
	ny=10
	init=2
	call ask_python(x,nx,y,ny,init,command)

	return
	end

