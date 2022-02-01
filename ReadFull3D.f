c long and latt are part of the Struct3D module.
c They contain the longitude and latitude boundaries of a cell 
c Note that nlong and nlatt are the number of longitude and latitude boundaries!
c This means there are (nlong-1) (nlatt-1) cells in longitude and latitude.
c The cell with number ilong,ilatt has boundaries:
c     Longitude: long(ilong) to long(ilong+1)
c     Latitude:  latt(ilatt) to latt(ilatt+1)
c Longitude cells run from the nightside where the first cell boundary is defined to be long(1) = 0
c Other cell boundaries are distributed equally over long = 0 to 2*pi
c This means that:
c - for odd values of nlong there is no cell directly centered at the substellar point.
c - for even values of nlong the cell with ilong=nlong/2 is centered at the substellar point.
c Latitudinal cells run from equator to pole where the first cell boundary is the equator and the final one is the pole.
c Assumption is that North and South halves of the planet are the same.
c The exact content of the modules can be found in Modules.f
	subroutine DoReadFull3D(i,ilong,ilatt)
	use GlobalSetup
	use Constants
	use Struct3D
	IMPLICIT NONE
	integer i,ilong,ilatt
	
c Parameters to redefine are:
	TPfile='whatever filename is now needed as a function of long and latt'
	Cloud(1)%file='whatever filename is now needed as a function of long and latt'

	return
	end
	