	subroutine WriteOutput()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i
	character*500 filename
	character*6000 form
	real*8,allocatable :: theta(:)
	real*8 Fp1,Fp2,ApAs,wr(nr)
	logical,allocatable :: docloud0(:,:)
	real*8,allocatable :: spec(:,:),specR(:),lamR(:),specRexp(:),specErr(:),Fstar_obs(:)
	real*8 x,specres_obs,expspecres_obs,gasdev,tot,Dmirror,f_phot,noisefloor,molweight(nmol),Tweight,Pweight
	real*8 lam_out(nlam)
	integer nlam_out
	integer ilam,j,nj,nlamR,i_instr,k,ir
	character*1000 line,instr_add
	character*10 side
	
	if(i2d.eq.0) then
		side=" "
	else
		write(side,'("_",i0.2)') i2d
	endif

	lam_out(1:nlam)=lam(1:nlam)/micron
	nlam_out=nlam

	allocate(docloud0(max(nclouds,1),ncc))
	allocate(theta(nphase))
	do i=1,nphase
		theta(i)=theta_phase(i)
	enddo
	do i=1,nclouds
		docloud0(i,:)=docloud(:,i)
	enddo
	call output("==================================================================")

	filename=trim(outputdir) // "star" // trim(side)
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	write(30,'("#",a13,a19)') "lambda [mu]","flux_star[Jy]"
	form='(f14.6,es19.7E3)'
	do i=1,nlam_out
		if(computelam(i)) then
			write(30,form) lam_out(i),Fstar(i)*1d23/distance**2
		endif
	enddo
	close(unit=30)

	if(emisspec) then
	
	filename=trim(outputdir) // "emis" // trim(side)
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	if(nclouds.gt.0) then
		form='("#",a13,a19,' // trim(int2string(ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]","flux [Jy]",docloud0(1:nclouds,1:ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","flux [Jy]"
	endif
	form='(f14.6,' // int2string(ncc+1,'(i3)') // 'es19.7E3)'
	do i=1,nlam_out
		if(lamemis(i).and.computelam(i)) then
		write(30,form) lam_out(i),
c     &					flux(0:ncc,i)
c     &					4d0*pi*1d-34*(phase(1,0,i)+flux(0,i))*clight*distance**2/(lam(i)*lam(i+1))
     &					(phase(1,j,i)+flux(j,i),j=0,ncc)
		endif
	enddo
	close(unit=30)

	filename=trim(outputdir) // "emisR" // trim(side)
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	if(nclouds.gt.0) then
		form='("#",a13,a19,' // trim(int2string(ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]","flux/flux_star",docloud0(1:nclouds,1:ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","flux/flux_star"
	endif
	form='(f14.6,' // int2string(ncc+1,'(i3)') // 'es19.7E3)'
	do i=1,nlam_out
		if(lamemis(i).and.computelam(i)) then
		write(30,form) lam_out(i),
c     &					flux(0:ncc,i)/(Fstar(i)*1d23/distance**2)
     &					((phase(1,j,i)+flux(j,i))/(Fstar(i)*1d23/distance**2),j=0,ncc)
		endif
	enddo
	close(unit=30)

	endif
	
	if(transspec) then

	filename=trim(outputdir) // "trans" // trim(side)
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	if(nclouds.gt.0) then
		form='("#",a13,a19,' // trim(int2string(ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]","Rp^2/Rstar^2",docloud0(1:nclouds,1:ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","Rp^2/Rstar^2"
	endif
	form='(f14.6,' // int2string(ncc+1,'(i3)') // 'es19.7E3)'
	do i=1,nlam_out
		if(computelam(i)) then
			write(30,form) lam_out(i),
     &					obsA(0:ncc,i)/(pi*Rstar**2)
    	endif
	enddo
	close(unit=30)

	if(do3D) then
	filename=trim(outputdir) // "trans_split" // trim(side)
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	write(30,'("#",a13,3a19)') "lambda [mu]","Morning","Evening","Rp^2/Rstar^2"
	form='(f14.6,3es19.7E3)'
	do i=1,nlam_out
		if(computelam(i)) then
			write(30,form) lam_out(i),obsA_split(i,1:2)/(pi*Rstar**2),
     &					obsA(0,i)/(pi*Rstar**2)
     	endif
	enddo
	close(unit=30)
	endif

	endif

	if(emisspec) then

	if(nphase.le.310) then
		filename=trim(outputdir) // "phase" // trim(side)
		call output("Writing spectrum to: " // trim(filename))
		open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
		form='("#",a13,' // trim(int2string(nphase,'(i4)')) // 
     &				 '("   flux(",f5.1,") [Jy]"),"         fstar [Jy]")'
		write(30,form) "lambda [mu]",theta(1:nphase)
		form='(f14.6,' // int2string(nphase+2,'(i3)') // 'es19.7E3)'
		do i=1,nlam_out
			if(lamemis(i).and.computelam(i)) then
			write(30,form) lam_out(i),
     &					phase(1:nphase,0,i)+flux(0,i),
     &					Fstar(i)*1d23/distance**2,
     &					(pi*Rplanet**2)*Fstar(i)*1d23/distance**2/(4d0*Dplanet**2)
			endif
		enddo
		close(unit=30)

		if(computealbedo) then
			filename=trim(outputdir) // "albedo" // trim(side)
			call output("Writing albedo to: " // trim(filename))
			open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
			form='("#",a13,' // trim(int2string(nphase,'(i4)')) // 
     &				 '("      albedo(",f5.1,")"),"         fstar [Jy]")'
			write(30,form) "lambda [mu]",theta(1:nphase)
			form='(f14.6,' // int2string(nphase+1,'(i3)') // 'es19.7E3)'
			do i=1,nlam_out
				if(lamemis(i).and.computelam(i)) then
				write(30,form) lam_out(i),planet_albedo(1:nphase,i),Fstar(i)*1d23/distance**2
				endif
			enddo
			close(unit=30)
		endif
	endif

	nj=0
	do i=1,nlam_out
		if(lamemis(i).and.computelam(i))nj=nj+1
	enddo
	if(nj.lt.350) then
		allocate(specR(nlam))
		filename=trim(outputdir) // "phasecurve" // trim(side)
		call output("Writing phasecurve to: " // trim(filename))
		open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
		form='("#",a13,' // trim(int2string(nj,'(i4)')) // 
     &				 '("      F(",es8.1E3,")"))'
		specR(1:nlam_out)=lam_out(1:nlam_out)
		nj=0
		do i=1,nlam_out
			if(lamemis(i).and.computelam(i)) then
				nj=nj+1
				specR(nj)=specR(i)
			endif
		enddo
		write(30,form) "phase [degrees]",specR(1:nj)
		form='(f14.6,' // int2string(nj,'(i3)') // 'es17.9E3)'
		do i=1,nphase
			specR(1:nlam_out)=Fstar(1:nlam_out)*1d23/distance**2
			specR(1:nlam_out)=specR(1:nlam_out)+phase(i,0,1:nlam_out)+flux(0,1:nlam_out)
			specR(1:nlam_out)=specR(1:nlam_out)/(Fstar(1:nlam_out)*1d23/distance**2)
			nj=0
			do j=1,nlam_out
				if(lamemis(j).and.computelam(j)) then
					nj=nj+1
					specR(nj)=specR(j)
				endif
			enddo
			write(30,form) theta(i),specR(1:nj)
		enddo
		close(unit=30)
		deallocate(specR)
	endif

	endif
	
	filename=trim(outputdir) // "tau1depth" // trim(side)
	call output("Writing tau1depth to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	if(nclouds.gt.0) then
		form='("#",a13,' // trim(int2string(ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]",docloud0(1:nclouds,1:ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","P [bar]"
	endif
	form='(f14.6,' // int2string(ncc,'(i3)') // 'es19.7E3)'
	do i=1,nlam_out
		if(computelam(i)) then
		write(30,form) lam_out(i),tau1depth(1:ncc,i)
		endif
	enddo
	close(unit=30)


	filename=trim(outputdir) // "cloudtau" // trim(side)
	call output("Writing cloud optical depth to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	if(nclouds.gt.0) then
		form='("#",a13,' // trim(int2string(ncc,'(i3)')) // 
     &				 '(' // trim(int2string(19-nclouds,'(i3)')) // '(" "),' // 
     &				trim(int2string(nclouds,'(i3)')) // 'l1))'
		write(30,form) "lambda [mu]",docloud0(1:nclouds,1:ncc)
	else
		write(30,'("#",a13,a19)') "lambda [mu]","optical depth"
	endif
	form='(f14.6,' // int2string(ncc,'(i3)') // 'es19.7E3)'
	do i=1,nlam_out
		if(computelam(i)) then
		write(30,form) lam_out(i),cloudtau(1:ncc,i)
		endif
	enddo
	close(unit=30)

	if(n_instr.gt.0) then

	do i_instr=1,n_instr
	
	instr_add=instrument(i_instr)
	select case(instrument(i_instr))
		case("ARIEL")
			nlamR=52
			allocate(lamR(nlamR))
			allocate(specR(nlamR))
			allocate(specRexp(nlamR))
			allocate(specErr(nlamR))
			call ARIELspecres(lamR,specR,specRexp)
			Dmirror=1d0
			f_phot=1d0/1.3d0
			noisefloor=2d-6
		case("JWST")
			nlamR=470
			allocate(lamR(nlamR))
			allocate(specR(nlamR))
			allocate(specRexp(nlamR))
			allocate(specErr(nlamR))
			call JWSTspecres(lamR,specR,specRexp)
			Dmirror=6.5d0
			f_phot=0.25d0
			noisefloor=20d-6
		case("MIRI")
			nlamR=308
			allocate(lamR(nlamR))
			allocate(specR(nlamR))
			allocate(specRexp(nlamR))
			allocate(specErr(nlamR))
			call MIRIspecres(lamR,specR,specRexp)
			Dmirror=6.5d0
			f_phot=0.25d0
			noisefloor=20d-6
		case("NIRSPEC")
			nlamR=162
			allocate(lamR(nlamR))
			allocate(specR(nlamR))
			allocate(specRexp(nlamR))
			allocate(specErr(nlamR))
			call NIRSPECspecres(lamR,specR,specRexp)
			Dmirror=6.5d0
			f_phot=0.25d0
			noisefloor=20d-6
		case("WFC3")
			nlamR=13
			allocate(lamR(nlamR))
			allocate(specR(nlamR))
			allocate(specRexp(nlamR))
			allocate(specErr(nlamR))
			call WFC3specres(lamR,specR,specRexp)
			Dmirror=2.4d0
			f_phot=0.1d0
			noisefloor=100d-6
		case("obs","OBS")
			instr_add="obs" // trim(int2string(instr_nobs(i_instr),'(i0.3)'))
			instr_ntrans(i_instr)=1d0
			noisefloor=0d0
			nlamR=ObsSpec(instr_nobs(i_instr))%ndata
			allocate(lamR(nlamR))
			allocate(specR(nlamR))
			allocate(specRexp(nlamR))
			allocate(specErr(nlamR))
			lamR(1:nlamR)=ObsSpec(instr_nobs(i_instr))%lam(1:nlamR)
			specRexp(1:nlamR)=ObsSpec(instr_nobs(i_instr))%Rexp(1:nlamR)
			specR(1:nlamR)=ObsSpec(instr_nobs(i_instr))%R(1:nlamR)
			specErr(1:nlamR)=ObsSpec(instr_nobs(i_instr))%dy(1:nlamR)
		case default
			instr_add="simulated"
			instr_ntrans(i_instr)=1d0
			noisefloor=0d0
			call CountSimInstrument(instrument(i_instr),nlamR)
			allocate(lamR(nlamR))
			allocate(specR(nlamR))
			allocate(specRexp(nlamR))
			allocate(specErr(nlamR))
			call ReadSimInstrument(instrument(i_instr),lamR,specR,specRexp,specErr,nlamR)
	end select
	allocate(spec(nphase,nlamR))
	allocate(Fstar_obs(nlamR))
	call regridspecres(lam,Fstar(1:nlam_out),nlam_out,
     &						lamR,Fstar_obs(1:nlamR),specR,specRexp,nlamR)
	if(instr_add.ne."simulated".and.instr_add(1:3).ne."obs") then
		do i=1,nlamR
			tot=1.51d7*(Fstar_obs(i)*1d23/distance**2)
			tot=tot*(pi*(Dmirror/2d0)**2)
			tot=tot*2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet)
			tot=tot*max(1d0,instr_ntrans(i_instr))*f_phot/specR(i)
			specErr(i)=1d0/sqrt(tot)
			if(specErr(i).lt.noisefloor) specErr(i)=noisefloor
		enddo
	endif

	if(computecontrib.and.allocated(obsA_contr)) then
		molweight=0d0
		Tweight=0d0
		Pweight=0d0
		tot=0d0
		do ir=1,nr
			wr(ir)=0d0
			if(sum(obsA_contr(ir,1:nlam_out)).gt.0d0) then
				call regridspecres(lam,obsA_contr(ir,1:nlam_out),nlam_out,lamR,spec(1,1:nlamR),specR,specRexp,nlamR)
				do j=1,nlamR
					x=(spec(1,j)/specErr(j))**2
					tot=tot+x
					wr(ir)=wr(ir)+x
					Tweight=Tweight+x*T(ir)
					Pweight=Pweight+x*P(ir)
					do i=1,nmol
						molweight(i)=molweight(i)+x*mixrat_r(ir,i)
					enddo
				enddo
			endif
		enddo
		wr=wr/sum(wr)
		Tweight=Tweight/tot
		Pweight=Pweight/tot
		molweight=molweight/tot
		call output("T average : " // dbl2string(Tweight,'(es10.3)'))
		call output("P average : " // dbl2string(Pweight,'(es10.3)'))
		do i=1,nmol
			if(includemol(i)) then
				call output(molname(i) // ": " // dbl2string(molweight(i),'(es10.3)'))
			endif
		enddo
		filename=trim(outputdir) // "contr_trans_" // trim(instr_add) // trim(side)
		open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
		write(30,'("#",a13,f10.3)') "T average ",Tweight
		write(30,'("#",a13,es10.3)') "P average ",Pweight
		do i=1,nmol
			if(includemol(i)) then
				write(30,'("#",a13,es10.3)') molname(i),molweight(i)
			endif
		enddo
		do ir=1,nr
			write(30,*) P(ir),wr(ir)
		enddo
		close(unit=30)



		molweight=0d0
		Tweight=0d0
		Pweight=0d0
		tot=0d0
		do ir=1,nr
			wr(ir)=0d0
			if(sum(obsA_contr(ir,1:nlam_out)).gt.0d0) then
				call regridspecres(lam,flux_contr(ir,1:nlam_out),nlam_out,lamR,spec(1,1:nlamR),specR,specRexp,nlamR)
				do j=1,nlamR
					x=(spec(1,j)/specErr(j))**2
					tot=tot+x
					wr(ir)=wr(ir)+x
					Tweight=Tweight+x*T(ir)
					Pweight=Pweight+x*P(ir)
					do i=1,nmol
						molweight(i)=molweight(i)+x*mixrat_r(ir,i)
					enddo
				enddo
			endif
		enddo
		wr=wr/sum(wr)
		Tweight=Tweight/tot
		Pweight=Pweight/tot
		molweight=molweight/tot
		call output("T average : " // dbl2string(Tweight,'(es10.3)'))
		call output("P average : " // dbl2string(Pweight,'(es10.3)'))
		do i=1,nmol
			if(includemol(i)) then
				call output(molname(i) // ": " // dbl2string(molweight(i),'(es10.3)'))
			endif
		enddo
		filename=trim(outputdir) // "contr_emisR_" // trim(instr_add) // trim(side)
		open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
		write(30,'("#",a13,f10.3)') "T average ",Tweight
		write(30,'("#",a13,es10.3)') "P average ",Pweight
		do i=1,nmol
			if(includemol(i)) then
				write(30,'("#",a13,es10.3)') molname(i),molweight(i)
			endif
		enddo
		do ir=1,nr
			write(30,*) P(ir),wr(ir)
		enddo
		close(unit=30)
	endif

	call regridspecres(lam,obsA(0,1:nlam_out),nlam_out,
     &					lamR,spec(1,1:nlamR),specR,specRexp,nlamR)
	spec=spec/(pi*Rstar**2)

	k=1
	if(instr_ntrans(i_instr).lt.1d0) then
		do ir=1,nr
			if(P(ir).gt.1d0.and.P(ir+1).le.1d0) exit
		enddo
		x=(sqrt(Rstar/(2d0*Dplanet))*Tstar*kb)/((Ggrav*Mplanet/(Rplanet**2))*mp*2.3d0)
		do i=1,nlamR
			do while((sqrt(real(k))*(5d0*x*Rplanet/(Rstar**2))/specErr(i)).lt.(7d0*0.9d0))
				k=k+1
			enddo
		enddo
		call output("assuming " // trim(int2string(k,'(i5)')) // " orbits")
	endif
	filename=trim(outputdir) // "obs_trans_" // trim(instr_add) // trim(side)
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	if(instr_add.ne."simulated") then
		write(30,'("# transit time       : ",f10.3," sec")') 2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet)
		write(30,'("# number of transits : ",f10.3)') instr_ntrans(i_instr)
		write(30,'("# integration time   : ",f10.3," hours")') 
     &			(instr_ntrans(i_instr)*2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet))/3600d0
	endif
	write(30,'("#",a13,4a19)') "lambda [mu]","Rp^2/Rstar^2","error","R"
	form='(f14.6,4es19.7E3)'
	do i=1,nlamR
		write(30,form) lamR(i)/micron,spec(1,i),specErr(i)/sqrt(real(k)),specR(i)
	enddo
	close(unit=30)
	do i=1,nlamR
		spec(1,i)=spec(1,i)+gasdev(idum)*specErr(i)/sqrt(real(k))
	enddo
	filename=trim(outputdir) // "obs_trans_noise_" // trim(instr_add) // trim(side)
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	if(instr_add.ne."simulated") then
		write(30,'("# transit time       : ",f10.3," sec")') 2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet)
		write(30,'("# number of transits : ",f10.3)') instr_ntrans(i_instr)
		write(30,'("# integration time   : ",f10.3," hours")') 
     &			(instr_ntrans(i_instr)*2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet))/3600d0
	endif
	write(30,'("#",a13,4a19)') "lambda [mu]","Rp^2/Rstar^2","error","R"
	form='(f14.6,4es19.7E3)'
	do i=1,nlamR
		write(30,form) lamR(i)/micron,spec(1,i),specErr(i)/sqrt(real(k)),specR(i)
	enddo
	close(unit=30)

	do i=1,nphase
		call regridspecres(lam,phase(i,0,1:nlam_out)+flux(0,1:nlam_out),nlam_out,
     &						lamR,spec(i,1:nlamR),specR,specRexp,nlamR)
	enddo
	call regridspecres(lam,Fstar(1:nlam_out),nlam_out,
     &						lamR,Fstar_obs(1:nlamR),specR,specRexp,nlamR)
	filename=trim(outputdir) // "obs_emis_" // trim(instr_add) // trim(side)
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	if(instr_add.ne."simulated") then
		write(30,'("# transit time       : ",es19.7E3, "sec")') 2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet)
		write(30,'("# number of transits : ",es19.7E3)') instr_ntrans(i_instr)
		write(30,'("# integration time   : ",es19.7E3,"hours")') 
     &		2d0*(instr_ntrans(i_instr)*2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet))/3600d0
	endif
	form='("#",a13,"        flux [Jy]","         error","             R")'
	write(30,form) "lambda [mu]"
	form='(f14.6,3es19.7E3)'
	do i=1,nlamR
c		write(30,form) lamR(i)/micron,4d0*pi*1d-34*spec(1,i)*clight*distance**2/lamR(i)**2,specErr(i)/sqrt(real(k)),specR(i)
		write(30,form) lamR(i)/micron,spec(1,i),specErr(i)/sqrt(real(k)),specR(i)
	enddo
	close(unit=30)
	do j=1,nphase
		spec(j,1:nlamR)=spec(j,1:nlamR)/(Fstar_obs(1:nlamR)*1d23/distance**2)
	enddo
	filename=trim(outputdir) // "obs_emisR_" // trim(instr_add) // trim(side)
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	if(instr_add.ne."simulated") then
		write(30,'("# transit time       : ",es19.7E3, "sec")') 2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet)
		write(30,'("# number of transits : ",es19.7E3)') instr_ntrans(i_instr)
		write(30,'("# integration time   : ",es19.7E3,"hours")') 
     &		2d0*(instr_ntrans(i_instr)*2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet))/3600d0
	endif
	form='("#",a13,' // trim(int2string(nphase,'(i4)')) // 
     &				 '("   flux(",f5.1,") [Jy]"),"         error","             R")'
	write(30,form) "lambda [mu]",theta(1:nphase)
	form='(f14.6,' // trim(int2string(nphase+2,'(i4)')) // 'es19.7E3)'
	do i=1,nlamR
		write(30,form) lamR(i)/micron,spec(1:nphase,i),specErr(i)/sqrt(real(k)),specR(i)
	enddo
	close(unit=30)
	do j=1,nphase
		do i=1,nlamR
			spec(j,i)=spec(j,i)+gasdev(idum)*specErr(i)/sqrt(real(k))
		enddo
	enddo
	filename=trim(outputdir) // "obs_emisR_noise_" // trim(instr_add) // trim(side)
	call output("Writing spectrum to: " // trim(filename))
	open(unit=30,file=filename,FORM="FORMATTED",ACCESS="STREAM")
	if(instr_add.ne."simulated") then
		write(30,'("# transit time       : ",es19.7E3, "sec")') 2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet)
		write(30,'("# number of transits : ",es19.7E3)') instr_ntrans(i_instr)
		write(30,'("# integration time   : ",es19.7E3,"hours")') 
     &		2d0*(instr_ntrans(i_instr)*2d0*pi*sqrt(Dplanet**3/(Ggrav*Mstar))*Rstar/(pi*Dplanet))/3600d0
	endif
	form='("#",a13,' // trim(int2string(nphase,'(i4)')) // 
     &				 '("   flux(",f5.1,") [Jy]"),"         error","             R")'
	write(30,form) "lambda [mu]",theta(1:nphase)
	form='(f14.6,' // trim(int2string(nphase+2,'(i4)')) // 'es19.7E3)'
	do i=1,nlamR
		write(30,form) lamR(i)/micron,spec(1:nphase,i),specErr(i)/sqrt(real(k)),specR(i)
	enddo
	close(unit=30)

	deallocate(lamR)
	deallocate(specR)
	deallocate(specRexp)
	deallocate(specErr)
	deallocate(spec)
	deallocate(Fstar_obs)

	enddo

	endif
	

	deallocate(docloud0)
	deallocate(theta)
	call output("==================================================================")
	
	
	return
	end
	
	


	subroutine WriteOpacity(ir,flag,nu0,kappa0,nnu0,ng0)
	use GlobalSetup
	IMPLICIT NONE
	character*500 file
	character*4 flag
	integer nnu0,i,ir,ng0,j
	real*8 nu0(nnu0),kappa0(nnu0,ng0)
	
	file=trim(outputdir) // "opacity_" // trim(flag) // "_" // trim(int2string(ir,'(i0.4)')) // ".dat"
	open(unit=30,file=file,FORM="FORMATTED",ACCESS="STREAM")
	write(30,'("# Pressure:    ",es10.3E3," bar")') P(ir)
	write(30,'("# Temperature: ",f10.3," K")') T(ir)
	write(30,'("#",a13,a19)') "lambda [mu]","kappa [cm^2/mol]"
	do i=1,nnu0
		do j=1,ng0
			write(30,'(f12.6,es19.7E3)') 1d4/nu0(i),kappa0(i,j)
		enddo
	enddo
	close(unit=30)
	
	return
	end

	subroutine ARIELspecres(lam,R,Rexp)
	IMPLICIT NONE
	integer nlam,j
	parameter(nlam=52)
	real*8 lam(*),R(*),Rexp(*)
	real*8 l0(          52),R0(          52),e0(          52)
	data (l0(j),j=1,          52) /
     &  0.55000E+00, 0.70500E+00, 0.95500E+00, 0.11564E+01, 0.12749E+01, 
     &  0.14056E+01, 0.15497E+01, 0.17085E+01, 0.18836E+01, 0.19696E+01, 
     &  0.20092E+01, 0.20496E+01, 0.20908E+01, 0.21328E+01, 0.21757E+01, 
     &  0.22194E+01, 0.22640E+01, 0.23095E+01, 0.23559E+01, 0.24033E+01, 
     &  0.24516E+01, 0.25009E+01, 0.25511E+01, 0.26024E+01, 0.26547E+01, 
     &  0.27081E+01, 0.27625E+01, 0.28180E+01, 0.28747E+01, 0.29325E+01, 
     &  0.29914E+01, 0.30515E+01, 0.31129E+01, 0.31754E+01, 0.32393E+01, 
     &  0.33044E+01, 0.33708E+01, 0.34385E+01, 0.35077E+01, 0.35782E+01, 
     &  0.36501E+01, 0.37234E+01, 0.40322E+01, 0.43055E+01, 0.45973E+01, 
     &  0.49089E+01, 0.52416E+01, 0.55968E+01, 0.59762E+01, 0.63812E+01, 
     &  0.68137E+01, 0.72756E+01 /
	data (R0(j),j=1,          52) /
     &  0.55000E+01, 0.33571E+01, 0.32931E+01, 0.10256E+02, 0.10256E+02, 
     &  0.10256E+02, 0.10256E+02, 0.10256E+02, 0.10256E+02, 0.50251E+02, 
     &  0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 
     &  0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 
     &  0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 
     &  0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 
     &  0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 
     &  0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 0.50251E+02, 
     &  0.50251E+02, 0.50251E+02, 0.15254E+02, 0.15254E+02, 0.15254E+02, 
     &  0.15254E+02, 0.15254E+02, 0.15254E+02, 0.15254E+02, 0.15254E+02, 
     &  0.15254E+02, 0.15254E+02 /
	data (e0(j),j=1,          52) /
     &  0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 
     &  0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 
     &  0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 
     &  0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 
     &  0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 
     &  0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 
     &  0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 
     &  0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 
     &  0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 
     &  0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 0.20000E+02, 
     &  0.20000E+02, 0.20000E+02 /
	lam(1:nlam)=l0(1:nlam)/1d4
	R(1:nlam)=R0(1:nlam)
	Rexp(1:nlam)=e0(1:nlam)
	return
	end

	subroutine JWSTspecres(lam,R,Rexp)
	IMPLICIT NONE
	integer nlam,j
	parameter(nlam=470)
	real*8 lam(*),R(*),Rexp(*)
	real*8 l0(         470),R0(         470),e0(         470)
	data (l0(j),j=1,         470) /
     &  0.10000E+01, 0.10154E+01, 0.10493E+01, 0.10846E+01, 0.11213E+01, 
     &  0.11594E+01, 0.11991E+01, 0.12404E+01, 0.12822E+01, 0.13245E+01, 
     &  0.13671E+01, 0.14101E+01, 0.14535E+01, 0.14972E+01, 0.15413E+01, 
     &  0.15852E+01, 0.16287E+01, 0.16721E+01, 0.17151E+01, 0.17580E+01, 
     &  0.18006E+01, 0.18431E+01, 0.18853E+01, 0.19274E+01, 0.19694E+01, 
     &  0.20111E+01, 0.20526E+01, 0.20934E+01, 0.21336E+01, 0.21732E+01, 
     &  0.22123E+01, 0.22508E+01, 0.22890E+01, 0.23267E+01, 0.23640E+01, 
     &  0.24009E+01, 0.24375E+01, 0.24737E+01, 0.25097E+01, 0.25452E+01, 
     &  0.25802E+01, 0.26147E+01, 0.26487E+01, 0.26823E+01, 0.27154E+01, 
     &  0.27481E+01, 0.27805E+01, 0.28125E+01, 0.28442E+01, 0.28756E+01, 
     &  0.29067E+01, 0.29375E+01, 0.29680E+01, 0.29983E+01, 0.30283E+01, 
     &  0.30579E+01, 0.30871E+01, 0.31160E+01, 0.31445E+01, 0.31727E+01, 
     &  0.32005E+01, 0.32281E+01, 0.32554E+01, 0.32825E+01, 0.33092E+01, 
     &  0.33358E+01, 0.33621E+01, 0.33881E+01, 0.34140E+01, 0.34396E+01, 
     &  0.34651E+01, 0.34903E+01, 0.35154E+01, 0.35403E+01, 0.35650E+01, 
     &  0.35896E+01, 0.36139E+01, 0.36382E+01, 0.36623E+01, 0.36862E+01, 
     &  0.37100E+01, 0.37337E+01, 0.37572E+01, 0.37806E+01, 0.38039E+01, 
     &  0.38270E+01, 0.38500E+01, 0.38730E+01, 0.38958E+01, 0.39185E+01, 
     &  0.39411E+01, 0.39635E+01, 0.39859E+01, 0.40082E+01, 0.40304E+01, 
     &  0.40524E+01, 0.40742E+01, 0.40959E+01, 0.41175E+01, 0.41389E+01, 
     &  0.41602E+01, 0.41814E+01, 0.42024E+01, 0.42233E+01, 0.42441E+01, 
     &  0.42647E+01, 0.42853E+01, 0.43057E+01, 0.43260E+01, 0.43463E+01, 
     &  0.43664E+01, 0.43864E+01, 0.44063E+01, 0.44261E+01, 0.44459E+01, 
     &  0.44655E+01, 0.44851E+01, 0.45045E+01, 0.45239E+01, 0.45432E+01, 
     &  0.45624E+01, 0.45815E+01, 0.45970E+01, 0.46006E+01, 0.46195E+01, 
     &  0.46384E+01, 0.46572E+01, 0.46759E+01, 0.46935E+01, 0.46946E+01, 
     &  0.47132E+01, 0.47317E+01, 0.47501E+01, 0.47685E+01, 0.47821E+01, 
     &  0.47868E+01, 0.48050E+01, 0.48232E+01, 0.48414E+01, 0.48594E+01, 
     &  0.48650E+01, 0.48774E+01, 0.48953E+01, 0.49132E+01, 0.49310E+01, 
     &  0.49435E+01, 0.49488E+01, 0.49665E+01, 0.49841E+01, 0.50017E+01, 
     &  0.50185E+01, 0.50192E+01, 0.50367E+01, 0.50540E+01, 0.50713E+01, 
     &  0.50884E+01, 0.50901E+01, 0.51055E+01, 0.51224E+01, 0.51393E+01, 
     &  0.51561E+01, 0.51583E+01, 0.51728E+01, 0.51895E+01, 0.52060E+01, 
     &  0.52225E+01, 0.52248E+01, 0.52389E+01, 0.52553E+01, 0.52715E+01, 
     &  0.52877E+01, 0.52891E+01, 0.53511E+01, 0.54120E+01, 0.54708E+01, 
     &  0.55285E+01, 0.55850E+01, 0.56400E+01, 0.56945E+01, 0.57472E+01, 
     &  0.57997E+01, 0.58505E+01, 0.59012E+01, 0.59502E+01, 0.59992E+01, 
     &  0.60467E+01, 0.60942E+01, 0.61405E+01, 0.61866E+01, 0.62319E+01, 
     &  0.62768E+01, 0.63212E+01, 0.63651E+01, 0.64088E+01, 0.64517E+01, 
     &  0.64947E+01, 0.65369E+01, 0.65789E+01, 0.66206E+01, 0.66619E+01, 
     &  0.67031E+01, 0.67437E+01, 0.67842E+01, 0.68243E+01, 0.68642E+01, 
     &  0.69040E+01, 0.69432E+01, 0.69824E+01, 0.70213E+01, 0.70599E+01, 
     &  0.70985E+01, 0.71365E+01, 0.71745E+01, 0.72123E+01, 0.72498E+01, 
     &  0.72872E+01, 0.73243E+01, 0.73611E+01, 0.73980E+01, 0.74343E+01, 
     &  0.74706E+01, 0.75068E+01, 0.75426E+01, 0.75783E+01, 0.76138E+01, 
     &  0.76490E+01, 0.76842E+01, 0.77191E+01, 0.77537E+01, 0.77883E+01, 
     &  0.78226E+01, 0.78567E+01, 0.78908E+01, 0.79244E+01, 0.79579E+01, 
     &  0.79914E+01, 0.80244E+01, 0.80572E+01, 0.80901E+01, 0.81224E+01, 
     &  0.81545E+01, 0.81866E+01, 0.82183E+01, 0.82497E+01, 0.82811E+01, 
     &  0.83123E+01, 0.83430E+01, 0.83738E+01, 0.84044E+01, 0.84346E+01, 
     &  0.84647E+01, 0.84949E+01, 0.85246E+01, 0.85542E+01, 0.85838E+01, 
     &  0.86132E+01, 0.86423E+01, 0.86713E+01, 0.87004E+01, 0.87290E+01, 
     &  0.87575E+01, 0.87861E+01, 0.88144E+01, 0.88424E+01, 0.88705E+01, 
     &  0.88985E+01, 0.89261E+01, 0.89537E+01, 0.89812E+01, 0.90086E+01, 
     &  0.90357E+01, 0.90628E+01, 0.90899E+01, 0.91167E+01, 0.91433E+01, 
     &  0.91699E+01, 0.91965E+01, 0.92228E+01, 0.92489E+01, 0.92751E+01, 
     &  0.93013E+01, 0.93270E+01, 0.93527E+01, 0.93785E+01, 0.94041E+01, 
     &  0.94294E+01, 0.94547E+01, 0.94801E+01, 0.95053E+01, 0.95302E+01, 
     &  0.95551E+01, 0.95800E+01, 0.96048E+01, 0.96293E+01, 0.96538E+01, 
     &  0.96783E+01, 0.97027E+01, 0.97268E+01, 0.97509E+01, 0.97750E+01, 
     &  0.97991E+01, 0.98229E+01, 0.98466E+01, 0.98703E+01, 0.98940E+01, 
     &  0.99175E+01, 0.99408E+01, 0.99642E+01, 0.99875E+01, 0.10011E+02, 
     &  0.10034E+02, 0.10057E+02, 0.10080E+02, 0.10103E+02, 0.10125E+02, 
     &  0.10148E+02, 0.10171E+02, 0.10193E+02, 0.10216E+02, 0.10238E+02, 
     &  0.10260E+02, 0.10282E+02, 0.10305E+02, 0.10327E+02, 0.10349E+02, 
     &  0.10371E+02, 0.10392E+02, 0.10414E+02, 0.10436E+02, 0.10457E+02, 
     &  0.10479E+02, 0.10501E+02, 0.10522E+02, 0.10543E+02, 0.10565E+02, 
     &  0.10586E+02, 0.10607E+02, 0.10628E+02, 0.10649E+02, 0.10670E+02, 
     &  0.10691E+02, 0.10712E+02, 0.10733E+02, 0.10753E+02, 0.10774E+02, 
     &  0.10795E+02, 0.10815E+02, 0.10836E+02, 0.10856E+02, 0.10876E+02, 
     &  0.10897E+02, 0.10917E+02, 0.10937E+02, 0.10957E+02, 0.10977E+02, 
     &  0.10997E+02, 0.11017E+02, 0.11037E+02, 0.11057E+02, 0.11077E+02, 
     &  0.11097E+02, 0.11116E+02, 0.11136E+02, 0.11155E+02, 0.11175E+02, 
     &  0.11194E+02, 0.11214E+02, 0.11233E+02, 0.11252E+02, 0.11272E+02, 
     &  0.11291E+02, 0.11310E+02, 0.11329E+02, 0.11348E+02, 0.11367E+02, 
     &  0.11386E+02, 0.11405E+02, 0.11424E+02, 0.11443E+02, 0.11461E+02, 
     &  0.11480E+02, 0.11499E+02, 0.11517E+02, 0.11536E+02, 0.11554E+02, 
     &  0.11573E+02, 0.11591E+02, 0.11610E+02, 0.11628E+02, 0.11646E+02, 
     &  0.11665E+02, 0.11683E+02, 0.11701E+02, 0.11719E+02, 0.11737E+02, 
     &  0.11755E+02, 0.11773E+02, 0.11791E+02, 0.11809E+02, 0.11827E+02, 
     &  0.11845E+02, 0.11862E+02, 0.11880E+02, 0.11898E+02, 0.11915E+02, 
     &  0.11933E+02, 0.11951E+02, 0.11968E+02, 0.11986E+02, 0.12003E+02, 
     &  0.12020E+02, 0.12038E+02, 0.12055E+02, 0.12072E+02, 0.12090E+02, 
     &  0.12107E+02, 0.12124E+02, 0.12141E+02, 0.12158E+02, 0.12175E+02, 
     &  0.12193E+02, 0.12209E+02, 0.12226E+02, 0.12243E+02, 0.12260E+02, 
     &  0.12277E+02, 0.12294E+02, 0.12311E+02, 0.12327E+02, 0.12344E+02, 
     &  0.12361E+02, 0.12377E+02, 0.12394E+02, 0.12411E+02, 0.12427E+02, 
     &  0.12443E+02, 0.12460E+02, 0.12476E+02, 0.12493E+02, 0.12509E+02, 
     &  0.12525E+02, 0.12542E+02, 0.12558E+02, 0.12574E+02, 0.12590E+02, 
     &  0.12607E+02, 0.12623E+02, 0.12639E+02, 0.12655E+02, 0.12671E+02, 
     &  0.12687E+02, 0.12703E+02, 0.12719E+02, 0.12735E+02, 0.12750E+02, 
     &  0.12766E+02, 0.12782E+02, 0.12798E+02, 0.12814E+02, 0.12829E+02, 
     &  0.12845E+02, 0.12861E+02, 0.12876E+02, 0.12892E+02, 0.12908E+02, 
     &  0.12923E+02, 0.12939E+02, 0.12954E+02, 0.12970E+02, 0.12985E+02 /
	data (R0(j),j=1,         470) /
     &  0.32500E+02, 0.20589E+02, 0.15163E+02, 0.15077E+02, 0.14987E+02, 
     &  0.14893E+02, 0.14796E+02, 0.14919E+02, 0.15261E+02, 0.15611E+02, 
     &  0.15964E+02, 0.16321E+02, 0.16681E+02, 0.17044E+02, 0.17528E+02, 
     &  0.18135E+02, 0.18747E+02, 0.19355E+02, 0.19960E+02, 0.20562E+02, 
     &  0.21160E+02, 0.21756E+02, 0.22349E+02, 0.22939E+02, 0.23527E+02, 
     &  0.24158E+02, 0.24950E+02, 0.25855E+02, 0.26745E+02, 0.27623E+02, 
     &  0.28489E+02, 0.29343E+02, 0.30187E+02, 0.31021E+02, 0.31846E+02, 
     &  0.32663E+02, 0.33471E+02, 0.34272E+02, 0.35105E+02, 0.36071E+02, 
     &  0.37129E+02, 0.38172E+02, 0.39199E+02, 0.40213E+02, 0.41214E+02, 
     &  0.42202E+02, 0.43179E+02, 0.44145E+02, 0.45101E+02, 0.46047E+02, 
     &  0.46983E+02, 0.47912E+02, 0.48831E+02, 0.49743E+02, 0.50787E+02, 
     &  0.51971E+02, 0.53148E+02, 0.54310E+02, 0.55458E+02, 0.56592E+02, 
     &  0.57713E+02, 0.58822E+02, 0.59920E+02, 0.61006E+02, 0.62083E+02, 
     &  0.63149E+02, 0.64205E+02, 0.65252E+02, 0.66291E+02, 0.67321E+02, 
     &  0.68342E+02, 0.69357E+02, 0.70363E+02, 0.71363E+02, 0.72355E+02, 
     &  0.73341E+02, 0.74320E+02, 0.75292E+02, 0.76259E+02, 0.77219E+02, 
     &  0.78174E+02, 0.79123E+02, 0.80067E+02, 0.81005E+02, 0.81939E+02, 
     &  0.82867E+02, 0.83791E+02, 0.84710E+02, 0.85624E+02, 0.86534E+02, 
     &  0.87439E+02, 0.88342E+02, 0.89239E+02, 0.90181E+02, 0.91252E+02, 
     &  0.92399E+02, 0.93540E+02, 0.94673E+02, 0.95798E+02, 0.96915E+02, 
     &  0.98027E+02, 0.99131E+02, 0.10023E+03, 0.10132E+03, 0.10240E+03, 
     &  0.10348E+03, 0.10455E+03, 0.10562E+03, 0.10668E+03, 0.10773E+03, 
     &  0.10878E+03, 0.10982E+03, 0.11086E+03, 0.11189E+03, 0.11292E+03, 
     &  0.11395E+03, 0.11497E+03, 0.11598E+03, 0.11700E+03, 0.11803E+03, 
     &  0.11905E+03, 0.12006E+03, 0.23839E+02, 0.12108E+03, 0.12208E+03, 
     &  0.12308E+03, 0.12408E+03, 0.12508E+03, 0.25363E+02, 0.12607E+03, 
     &  0.12705E+03, 0.12804E+03, 0.12902E+03, 0.12999E+03, 0.27881E+02, 
     &  0.13096E+03, 0.13193E+03, 0.13290E+03, 0.13386E+03, 0.13482E+03, 
     &  0.30149E+02, 0.13577E+03, 0.13673E+03, 0.13767E+03, 0.13862E+03, 
     &  0.32211E+02, 0.13956E+03, 0.14051E+03, 0.14144E+03, 0.14239E+03, 
     &  0.34230E+02, 0.14353E+03, 0.14484E+03, 0.14615E+03, 0.14745E+03, 
     &  0.14874E+03, 0.36388E+02, 0.15002E+03, 0.15129E+03, 0.15256E+03, 
     &  0.15383E+03, 0.38294E+02, 0.15509E+03, 0.15634E+03, 0.15758E+03, 
     &  0.15882E+03, 0.39941E+02, 0.16005E+03, 0.16128E+03, 0.16250E+03, 
     &  0.16336E+03, 0.41875E+02, 0.43561E+02, 0.45193E+02, 0.46952E+02, 
     &  0.48409E+02, 0.50077E+02, 0.51531E+02, 0.53131E+02, 0.54599E+02, 
     &  0.56165E+02, 0.57692E+02, 0.59182E+02, 0.60676E+02, 0.62138E+02, 
     &  0.63647E+02, 0.64990E+02, 0.66458E+02, 0.67697E+02, 0.69090E+02, 
     &  0.70251E+02, 0.71587E+02, 0.72689E+02, 0.73988E+02, 0.75137E+02, 
     &  0.76306E+02, 0.77586E+02, 0.78553E+02, 0.79806E+02, 0.80744E+02, 
     &  0.81972E+02, 0.83157E+02, 0.84090E+02, 0.85312E+02, 0.86169E+02, 
     &  0.87380E+02, 0.88529E+02, 0.89418E+02, 0.90631E+02, 0.91455E+02, 
     &  0.92638E+02, 0.93859E+02, 0.94624E+02, 0.95853E+02, 0.96849E+02, 
     &  0.97835E+02, 0.99085E+02, 0.99844E+02, 0.10107E+03, 0.10235E+03, 
     &  0.10304E+03, 0.10435E+03, 0.10551E+03, 0.10633E+03, 0.10768E+03, 
     &  0.10870E+03, 0.10969E+03, 0.11108E+03, 0.11199E+03, 0.11306E+03, 
     &  0.11443E+03, 0.11523E+03, 0.11646E+03, 0.11798E+03, 0.11877E+03, 
     &  0.12018E+03, 0.12194E+03, 0.12277E+03, 0.12418E+03, 0.12600E+03, 
     &  0.12691E+03, 0.12826E+03, 0.13023E+03, 0.13136E+03, 0.13242E+03, 
     &  0.13435E+03, 0.13572E+03, 0.13641E+03, 0.13819E+03, 0.13980E+03, 
     &  0.14031E+03, 0.14189E+03, 0.14371E+03, 0.14444E+03, 0.14554E+03, 
     &  0.14738E+03, 0.14862E+03, 0.14914E+03, 0.15099E+03, 0.15285E+03, 
     &  0.15335E+03, 0.15455E+03, 0.15644E+03, 0.15762E+03, 0.15812E+03, 
     &  0.15995E+03, 0.16187E+03, 0.16245E+03, 0.16340E+03, 0.16533E+03, 
     &  0.16682E+03, 0.16732E+03, 0.16872E+03, 0.17068E+03, 0.17173E+03, 
     &  0.17223E+03, 0.17401E+03, 0.17600E+03, 0.17670E+03, 0.17727E+03, 
     &  0.17927E+03, 0.18122E+03, 0.18172E+03, 0.18247E+03, 0.18449E+03, 
     &  0.18628E+03, 0.18678E+03, 0.18760E+03, 0.18966E+03, 0.19140E+03, 
     &  0.19190E+03, 0.19270E+03, 0.19478E+03, 0.19656E+03, 0.19706E+03, 
     &  0.19774E+03, 0.19984E+03, 0.20178E+03, 0.20227E+03, 0.20278E+03, 
     &  0.20484E+03, 0.20699E+03, 0.20754E+03, 0.20805E+03, 0.20978E+03, 
     &  0.21194E+03, 0.21286E+03, 0.21336E+03, 0.21463E+03, 0.21682E+03, 
     &  0.21823E+03, 0.21873E+03, 0.21942E+03, 0.22162E+03, 0.22365E+03, 
     &  0.22415E+03, 0.22465E+03, 0.22636E+03, 0.22860E+03, 0.22962E+03, 
     &  0.23012E+03, 0.23099E+03, 0.23325E+03, 0.23515E+03, 0.23565E+03, 
     &  0.23615E+03, 0.23782E+03, 0.24011E+03, 0.24123E+03, 0.24173E+03, 
     &  0.24230E+03, 0.24460E+03, 0.24686E+03, 0.24736E+03, 0.24787E+03, 
     &  0.24898E+03, 0.25132E+03, 0.25305E+03, 0.25355E+03, 0.25405E+03, 
     &  0.25562E+03, 0.25799E+03, 0.25929E+03, 0.25979E+03, 0.26029E+03, 
     &  0.26220E+03, 0.26460E+03, 0.26558E+03, 0.26609E+03, 0.26659E+03, 
     &  0.26869E+03, 0.27113E+03, 0.27193E+03, 0.27243E+03, 0.27293E+03, 
     &  0.27512E+03, 0.27758E+03, 0.27834E+03, 0.27884E+03, 0.27934E+03, 
     &  0.28146E+03, 0.28396E+03, 0.28480E+03, 0.28530E+03, 0.28580E+03, 
     &  0.28772E+03, 0.29023E+03, 0.29131E+03, 0.29181E+03, 0.29231E+03, 
     &  0.29388E+03, 0.29641E+03, 0.29789E+03, 0.29839E+03, 0.29889E+03, 
     &  0.29993E+03, 0.30249E+03, 0.30451E+03, 0.30502E+03, 0.30552E+03, 
     &  0.30601E+03, 0.30846E+03, 0.31107E+03, 0.31170E+03, 0.31220E+03, 
     &  0.31269E+03, 0.31432E+03, 0.31695E+03, 0.31845E+03, 0.31894E+03, 
     &  0.31944E+03, 0.32007E+03, 0.32270E+03, 0.32524E+03, 0.32574E+03, 
     &  0.32625E+03, 0.32675E+03, 0.32834E+03, 0.33102E+03, 0.33260E+03, 
     &  0.33311E+03, 0.33361E+03, 0.33411E+03, 0.33654E+03, 0.33925E+03, 
     &  0.34002E+03, 0.34052E+03, 0.34102E+03, 0.34192E+03, 0.34465E+03, 
     &  0.34700E+03, 0.34750E+03, 0.34799E+03, 0.34849E+03, 0.34991E+03, 
     &  0.35268E+03, 0.35453E+03, 0.35503E+03, 0.35553E+03, 0.35603E+03, 
     &  0.35781E+03, 0.36061E+03, 0.36213E+03, 0.36262E+03, 0.36312E+03, 
     &  0.36362E+03, 0.36560E+03, 0.36844E+03, 0.36978E+03, 0.37028E+03, 
     &  0.37078E+03, 0.37128E+03, 0.37327E+03, 0.37613E+03, 0.37750E+03, 
     &  0.37799E+03, 0.37849E+03, 0.37899E+03, 0.38083E+03, 0.38372E+03, 
     &  0.38526E+03, 0.38576E+03, 0.38627E+03, 0.38677E+03, 0.38826E+03, 
     &  0.39118E+03, 0.39310E+03, 0.39360E+03, 0.39412E+03, 0.39462E+03, 
     &  0.39555E+03, 0.39848E+03, 0.40099E+03, 0.40150E+03, 0.40200E+03, 
     &  0.40250E+03, 0.40300E+03, 0.40566E+03, 0.40866E+03, 0.40945E+03, 
     &  0.40995E+03, 0.41047E+03, 0.41097E+03, 0.41268E+03, 0.41570E+03, 
     &  0.41749E+03, 0.41799E+03, 0.41849E+03, 0.41898E+03, 0.41946E+03 /
	data (e0(j),j=1,         470) /
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01 /
	lam(1:nlam)=l0(1:nlam)/1d4
	R(1:nlam)=R0(1:nlam)
	Rexp(1:nlam)=e0(1:nlam)
	return
	end

	subroutine MIRIspecres(lam,R,Rexp)
	IMPLICIT NONE
	integer nlam,j
	parameter(nlam=308)
	real*8 lam(*),R(*),Rexp(*)
	real*8 l0(         308),R0(         308),e0(         308)
	data (l0(j),j=1,         308) /
     &  0.45970E+01, 0.46935E+01, 0.47821E+01, 0.48650E+01, 0.49435E+01, 
     &  0.50185E+01, 0.50901E+01, 0.51583E+01, 0.52248E+01, 0.52891E+01, 
     &  0.53511E+01, 0.54120E+01, 0.54708E+01, 0.55285E+01, 0.55850E+01, 
     &  0.56400E+01, 0.56945E+01, 0.57472E+01, 0.57997E+01, 0.58505E+01, 
     &  0.59012E+01, 0.59502E+01, 0.59992E+01, 0.60467E+01, 0.60942E+01, 
     &  0.61405E+01, 0.61866E+01, 0.62319E+01, 0.62768E+01, 0.63212E+01, 
     &  0.63651E+01, 0.64088E+01, 0.64517E+01, 0.64947E+01, 0.65369E+01, 
     &  0.65789E+01, 0.66206E+01, 0.66619E+01, 0.67031E+01, 0.67437E+01, 
     &  0.67842E+01, 0.68243E+01, 0.68642E+01, 0.69040E+01, 0.69432E+01, 
     &  0.69824E+01, 0.70213E+01, 0.70599E+01, 0.70985E+01, 0.71365E+01, 
     &  0.71745E+01, 0.72123E+01, 0.72498E+01, 0.72872E+01, 0.73243E+01, 
     &  0.73611E+01, 0.73980E+01, 0.74343E+01, 0.74706E+01, 0.75068E+01, 
     &  0.75426E+01, 0.75783E+01, 0.76138E+01, 0.76490E+01, 0.76842E+01, 
     &  0.77191E+01, 0.77537E+01, 0.77883E+01, 0.78226E+01, 0.78567E+01, 
     &  0.78908E+01, 0.79244E+01, 0.79579E+01, 0.79914E+01, 0.80244E+01, 
     &  0.80572E+01, 0.80901E+01, 0.81224E+01, 0.81545E+01, 0.81866E+01, 
     &  0.82183E+01, 0.82497E+01, 0.82811E+01, 0.83123E+01, 0.83430E+01, 
     &  0.83738E+01, 0.84044E+01, 0.84346E+01, 0.84647E+01, 0.84949E+01, 
     &  0.85246E+01, 0.85542E+01, 0.85838E+01, 0.86132E+01, 0.86423E+01, 
     &  0.86713E+01, 0.87004E+01, 0.87290E+01, 0.87575E+01, 0.87861E+01, 
     &  0.88144E+01, 0.88424E+01, 0.88705E+01, 0.88985E+01, 0.89261E+01, 
     &  0.89537E+01, 0.89812E+01, 0.90086E+01, 0.90357E+01, 0.90628E+01, 
     &  0.90899E+01, 0.91167E+01, 0.91433E+01, 0.91699E+01, 0.91965E+01, 
     &  0.92228E+01, 0.92489E+01, 0.92751E+01, 0.93013E+01, 0.93270E+01, 
     &  0.93527E+01, 0.93785E+01, 0.94041E+01, 0.94294E+01, 0.94547E+01, 
     &  0.94801E+01, 0.95053E+01, 0.95302E+01, 0.95551E+01, 0.95800E+01, 
     &  0.96048E+01, 0.96293E+01, 0.96538E+01, 0.96783E+01, 0.97027E+01, 
     &  0.97268E+01, 0.97509E+01, 0.97750E+01, 0.97991E+01, 0.98229E+01, 
     &  0.98466E+01, 0.98703E+01, 0.98940E+01, 0.99175E+01, 0.99408E+01, 
     &  0.99642E+01, 0.99875E+01, 0.10011E+02, 0.10034E+02, 0.10057E+02, 
     &  0.10080E+02, 0.10103E+02, 0.10125E+02, 0.10148E+02, 0.10171E+02, 
     &  0.10193E+02, 0.10216E+02, 0.10238E+02, 0.10260E+02, 0.10282E+02, 
     &  0.10305E+02, 0.10327E+02, 0.10349E+02, 0.10371E+02, 0.10392E+02, 
     &  0.10414E+02, 0.10436E+02, 0.10457E+02, 0.10479E+02, 0.10501E+02, 
     &  0.10522E+02, 0.10543E+02, 0.10565E+02, 0.10586E+02, 0.10607E+02, 
     &  0.10628E+02, 0.10649E+02, 0.10670E+02, 0.10691E+02, 0.10712E+02, 
     &  0.10733E+02, 0.10753E+02, 0.10774E+02, 0.10795E+02, 0.10815E+02, 
     &  0.10836E+02, 0.10856E+02, 0.10876E+02, 0.10897E+02, 0.10917E+02, 
     &  0.10937E+02, 0.10957E+02, 0.10977E+02, 0.10997E+02, 0.11017E+02, 
     &  0.11037E+02, 0.11057E+02, 0.11077E+02, 0.11097E+02, 0.11116E+02, 
     &  0.11136E+02, 0.11155E+02, 0.11175E+02, 0.11194E+02, 0.11214E+02, 
     &  0.11233E+02, 0.11252E+02, 0.11272E+02, 0.11291E+02, 0.11310E+02, 
     &  0.11329E+02, 0.11348E+02, 0.11367E+02, 0.11386E+02, 0.11405E+02, 
     &  0.11424E+02, 0.11443E+02, 0.11461E+02, 0.11480E+02, 0.11499E+02, 
     &  0.11517E+02, 0.11536E+02, 0.11554E+02, 0.11573E+02, 0.11591E+02, 
     &  0.11610E+02, 0.11628E+02, 0.11646E+02, 0.11665E+02, 0.11683E+02, 
     &  0.11701E+02, 0.11719E+02, 0.11737E+02, 0.11755E+02, 0.11773E+02, 
     &  0.11791E+02, 0.11809E+02, 0.11827E+02, 0.11845E+02, 0.11862E+02, 
     &  0.11880E+02, 0.11898E+02, 0.11915E+02, 0.11933E+02, 0.11951E+02, 
     &  0.11968E+02, 0.11986E+02, 0.12003E+02, 0.12020E+02, 0.12038E+02, 
     &  0.12055E+02, 0.12072E+02, 0.12090E+02, 0.12107E+02, 0.12124E+02, 
     &  0.12141E+02, 0.12158E+02, 0.12175E+02, 0.12193E+02, 0.12209E+02, 
     &  0.12226E+02, 0.12243E+02, 0.12260E+02, 0.12277E+02, 0.12294E+02, 
     &  0.12311E+02, 0.12327E+02, 0.12344E+02, 0.12361E+02, 0.12377E+02, 
     &  0.12394E+02, 0.12411E+02, 0.12427E+02, 0.12443E+02, 0.12460E+02, 
     &  0.12476E+02, 0.12493E+02, 0.12509E+02, 0.12525E+02, 0.12542E+02, 
     &  0.12558E+02, 0.12574E+02, 0.12590E+02, 0.12607E+02, 0.12623E+02, 
     &  0.12639E+02, 0.12655E+02, 0.12671E+02, 0.12687E+02, 0.12703E+02, 
     &  0.12719E+02, 0.12735E+02, 0.12750E+02, 0.12766E+02, 0.12782E+02, 
     &  0.12798E+02, 0.12814E+02, 0.12829E+02, 0.12845E+02, 0.12861E+02, 
     &  0.12876E+02, 0.12892E+02, 0.12908E+02, 0.12923E+02, 0.12939E+02, 
     &  0.12954E+02, 0.12970E+02, 0.12985E+02 /
	data (R0(j),j=1,         308) /
     &  0.23839E+02, 0.25363E+02, 0.27881E+02, 0.30149E+02, 0.32211E+02, 
     &  0.34230E+02, 0.36388E+02, 0.38294E+02, 0.39941E+02, 0.41875E+02, 
     &  0.43561E+02, 0.45193E+02, 0.46952E+02, 0.48409E+02, 0.50077E+02, 
     &  0.51531E+02, 0.53131E+02, 0.54599E+02, 0.56165E+02, 0.57692E+02, 
     &  0.59182E+02, 0.60676E+02, 0.62138E+02, 0.63647E+02, 0.64990E+02, 
     &  0.66458E+02, 0.67697E+02, 0.69090E+02, 0.70251E+02, 0.71587E+02, 
     &  0.72689E+02, 0.73988E+02, 0.75137E+02, 0.76306E+02, 0.77586E+02, 
     &  0.78553E+02, 0.79806E+02, 0.80744E+02, 0.81972E+02, 0.83157E+02, 
     &  0.84090E+02, 0.85312E+02, 0.86169E+02, 0.87380E+02, 0.88529E+02, 
     &  0.89418E+02, 0.90631E+02, 0.91455E+02, 0.92638E+02, 0.93859E+02, 
     &  0.94624E+02, 0.95853E+02, 0.96849E+02, 0.97835E+02, 0.99085E+02, 
     &  0.99844E+02, 0.10107E+03, 0.10235E+03, 0.10304E+03, 0.10435E+03, 
     &  0.10551E+03, 0.10633E+03, 0.10768E+03, 0.10870E+03, 0.10969E+03, 
     &  0.11108E+03, 0.11199E+03, 0.11306E+03, 0.11443E+03, 0.11523E+03, 
     &  0.11646E+03, 0.11798E+03, 0.11877E+03, 0.12018E+03, 0.12194E+03, 
     &  0.12277E+03, 0.12418E+03, 0.12600E+03, 0.12691E+03, 0.12826E+03, 
     &  0.13023E+03, 0.13136E+03, 0.13242E+03, 0.13435E+03, 0.13572E+03, 
     &  0.13641E+03, 0.13819E+03, 0.13980E+03, 0.14031E+03, 0.14189E+03, 
     &  0.14371E+03, 0.14444E+03, 0.14554E+03, 0.14738E+03, 0.14862E+03, 
     &  0.14914E+03, 0.15099E+03, 0.15285E+03, 0.15335E+03, 0.15455E+03, 
     &  0.15644E+03, 0.15762E+03, 0.15812E+03, 0.15995E+03, 0.16187E+03, 
     &  0.16245E+03, 0.16340E+03, 0.16533E+03, 0.16682E+03, 0.16732E+03, 
     &  0.16872E+03, 0.17068E+03, 0.17173E+03, 0.17223E+03, 0.17401E+03, 
     &  0.17600E+03, 0.17670E+03, 0.17727E+03, 0.17927E+03, 0.18122E+03, 
     &  0.18172E+03, 0.18247E+03, 0.18449E+03, 0.18628E+03, 0.18678E+03, 
     &  0.18760E+03, 0.18966E+03, 0.19140E+03, 0.19190E+03, 0.19270E+03, 
     &  0.19478E+03, 0.19656E+03, 0.19706E+03, 0.19774E+03, 0.19984E+03, 
     &  0.20178E+03, 0.20227E+03, 0.20278E+03, 0.20484E+03, 0.20699E+03, 
     &  0.20754E+03, 0.20805E+03, 0.20978E+03, 0.21194E+03, 0.21286E+03, 
     &  0.21336E+03, 0.21463E+03, 0.21682E+03, 0.21823E+03, 0.21873E+03, 
     &  0.21942E+03, 0.22162E+03, 0.22365E+03, 0.22415E+03, 0.22465E+03, 
     &  0.22636E+03, 0.22860E+03, 0.22962E+03, 0.23012E+03, 0.23099E+03, 
     &  0.23325E+03, 0.23515E+03, 0.23565E+03, 0.23615E+03, 0.23782E+03, 
     &  0.24011E+03, 0.24123E+03, 0.24173E+03, 0.24230E+03, 0.24460E+03, 
     &  0.24686E+03, 0.24736E+03, 0.24787E+03, 0.24898E+03, 0.25132E+03, 
     &  0.25305E+03, 0.25355E+03, 0.25405E+03, 0.25562E+03, 0.25799E+03, 
     &  0.25929E+03, 0.25979E+03, 0.26029E+03, 0.26220E+03, 0.26460E+03, 
     &  0.26558E+03, 0.26609E+03, 0.26659E+03, 0.26869E+03, 0.27113E+03, 
     &  0.27193E+03, 0.27243E+03, 0.27293E+03, 0.27512E+03, 0.27758E+03, 
     &  0.27834E+03, 0.27884E+03, 0.27934E+03, 0.28146E+03, 0.28396E+03, 
     &  0.28480E+03, 0.28530E+03, 0.28580E+03, 0.28772E+03, 0.29023E+03, 
     &  0.29131E+03, 0.29181E+03, 0.29231E+03, 0.29388E+03, 0.29641E+03, 
     &  0.29789E+03, 0.29839E+03, 0.29889E+03, 0.29993E+03, 0.30249E+03, 
     &  0.30451E+03, 0.30502E+03, 0.30552E+03, 0.30601E+03, 0.30846E+03, 
     &  0.31107E+03, 0.31170E+03, 0.31220E+03, 0.31269E+03, 0.31432E+03, 
     &  0.31695E+03, 0.31845E+03, 0.31894E+03, 0.31944E+03, 0.32007E+03, 
     &  0.32270E+03, 0.32524E+03, 0.32574E+03, 0.32625E+03, 0.32675E+03, 
     &  0.32834E+03, 0.33102E+03, 0.33260E+03, 0.33311E+03, 0.33361E+03, 
     &  0.33411E+03, 0.33654E+03, 0.33925E+03, 0.34002E+03, 0.34052E+03, 
     &  0.34102E+03, 0.34192E+03, 0.34465E+03, 0.34700E+03, 0.34750E+03, 
     &  0.34799E+03, 0.34849E+03, 0.34991E+03, 0.35268E+03, 0.35453E+03, 
     &  0.35503E+03, 0.35553E+03, 0.35603E+03, 0.35781E+03, 0.36061E+03, 
     &  0.36213E+03, 0.36262E+03, 0.36312E+03, 0.36362E+03, 0.36560E+03, 
     &  0.36844E+03, 0.36978E+03, 0.37028E+03, 0.37078E+03, 0.37128E+03, 
     &  0.37327E+03, 0.37613E+03, 0.37750E+03, 0.37799E+03, 0.37849E+03, 
     &  0.37899E+03, 0.38083E+03, 0.38372E+03, 0.38526E+03, 0.38576E+03, 
     &  0.38627E+03, 0.38677E+03, 0.38826E+03, 0.39118E+03, 0.39310E+03, 
     &  0.39360E+03, 0.39412E+03, 0.39462E+03, 0.39555E+03, 0.39848E+03, 
     &  0.40099E+03, 0.40150E+03, 0.40200E+03, 0.40250E+03, 0.40300E+03, 
     &  0.40566E+03, 0.40866E+03, 0.40945E+03, 0.40995E+03, 0.41047E+03, 
     &  0.41097E+03, 0.41268E+03, 0.41570E+03, 0.41749E+03, 0.41799E+03, 
     &  0.41849E+03, 0.41898E+03, 0.41946E+03 /
	data (e0(j),j=1,         308) /
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01 /
	lam(1:nlam)=l0(1:nlam)/1d4
	R(1:nlam)=R0(1:nlam)
	Rexp(1:nlam)=e0(1:nlam)
	return
	end

	subroutine NIRSPECspecres(lam,R,Rexp)
	IMPLICIT NONE
	integer nlam,j
	parameter(nlam=162)
	real*8 lam(*),R(*),Rexp(*)
	real*8 l0(         162),R0(         162),e0(         162)
	data (l0(j),j=1,         162) /
     &  0.10000E+01, 0.10154E+01, 0.10493E+01, 0.10846E+01, 0.11213E+01, 
     &  0.11594E+01, 0.11991E+01, 0.12404E+01, 0.12822E+01, 0.13245E+01, 
     &  0.13671E+01, 0.14101E+01, 0.14535E+01, 0.14972E+01, 0.15413E+01, 
     &  0.15852E+01, 0.16287E+01, 0.16721E+01, 0.17151E+01, 0.17580E+01, 
     &  0.18006E+01, 0.18431E+01, 0.18853E+01, 0.19274E+01, 0.19694E+01, 
     &  0.20111E+01, 0.20526E+01, 0.20934E+01, 0.21336E+01, 0.21732E+01, 
     &  0.22123E+01, 0.22508E+01, 0.22890E+01, 0.23267E+01, 0.23640E+01, 
     &  0.24009E+01, 0.24375E+01, 0.24737E+01, 0.25097E+01, 0.25452E+01, 
     &  0.25802E+01, 0.26147E+01, 0.26487E+01, 0.26823E+01, 0.27154E+01, 
     &  0.27481E+01, 0.27805E+01, 0.28125E+01, 0.28442E+01, 0.28756E+01, 
     &  0.29067E+01, 0.29375E+01, 0.29680E+01, 0.29983E+01, 0.30283E+01, 
     &  0.30579E+01, 0.30871E+01, 0.31160E+01, 0.31445E+01, 0.31727E+01, 
     &  0.32005E+01, 0.32281E+01, 0.32554E+01, 0.32825E+01, 0.33092E+01, 
     &  0.33358E+01, 0.33621E+01, 0.33881E+01, 0.34140E+01, 0.34396E+01, 
     &  0.34651E+01, 0.34903E+01, 0.35154E+01, 0.35403E+01, 0.35650E+01, 
     &  0.35896E+01, 0.36139E+01, 0.36382E+01, 0.36623E+01, 0.36862E+01, 
     &  0.37100E+01, 0.37337E+01, 0.37572E+01, 0.37806E+01, 0.38039E+01, 
     &  0.38270E+01, 0.38500E+01, 0.38730E+01, 0.38958E+01, 0.39185E+01, 
     &  0.39411E+01, 0.39635E+01, 0.39859E+01, 0.40082E+01, 0.40304E+01, 
     &  0.40524E+01, 0.40742E+01, 0.40959E+01, 0.41175E+01, 0.41389E+01, 
     &  0.41602E+01, 0.41814E+01, 0.42024E+01, 0.42233E+01, 0.42441E+01, 
     &  0.42647E+01, 0.42853E+01, 0.43057E+01, 0.43260E+01, 0.43463E+01, 
     &  0.43664E+01, 0.43864E+01, 0.44063E+01, 0.44261E+01, 0.44459E+01, 
     &  0.44655E+01, 0.44851E+01, 0.45045E+01, 0.45239E+01, 0.45432E+01, 
     &  0.45624E+01, 0.45815E+01, 0.46006E+01, 0.46195E+01, 0.46384E+01, 
     &  0.46572E+01, 0.46759E+01, 0.46946E+01, 0.47132E+01, 0.47317E+01, 
     &  0.47501E+01, 0.47685E+01, 0.47868E+01, 0.48050E+01, 0.48232E+01, 
     &  0.48414E+01, 0.48594E+01, 0.48774E+01, 0.48953E+01, 0.49132E+01, 
     &  0.49310E+01, 0.49488E+01, 0.49665E+01, 0.49841E+01, 0.50017E+01, 
     &  0.50192E+01, 0.50367E+01, 0.50540E+01, 0.50713E+01, 0.50884E+01, 
     &  0.51055E+01, 0.51224E+01, 0.51393E+01, 0.51561E+01, 0.51728E+01, 
     &  0.51895E+01, 0.52060E+01, 0.52225E+01, 0.52389E+01, 0.52553E+01, 
     &  0.52715E+01, 0.52877E+01 /
	data (R0(j),j=1,         162) /
     &  0.32500E+02, 0.20589E+02, 0.15163E+02, 0.15077E+02, 0.14987E+02, 
     &  0.14893E+02, 0.14796E+02, 0.14919E+02, 0.15261E+02, 0.15611E+02, 
     &  0.15964E+02, 0.16321E+02, 0.16681E+02, 0.17044E+02, 0.17528E+02, 
     &  0.18135E+02, 0.18747E+02, 0.19355E+02, 0.19960E+02, 0.20562E+02, 
     &  0.21160E+02, 0.21756E+02, 0.22349E+02, 0.22939E+02, 0.23527E+02, 
     &  0.24158E+02, 0.24950E+02, 0.25855E+02, 0.26745E+02, 0.27623E+02, 
     &  0.28489E+02, 0.29343E+02, 0.30187E+02, 0.31021E+02, 0.31846E+02, 
     &  0.32663E+02, 0.33471E+02, 0.34272E+02, 0.35105E+02, 0.36071E+02, 
     &  0.37129E+02, 0.38172E+02, 0.39199E+02, 0.40213E+02, 0.41214E+02, 
     &  0.42202E+02, 0.43179E+02, 0.44145E+02, 0.45101E+02, 0.46047E+02, 
     &  0.46983E+02, 0.47912E+02, 0.48831E+02, 0.49743E+02, 0.50787E+02, 
     &  0.51971E+02, 0.53148E+02, 0.54310E+02, 0.55458E+02, 0.56592E+02, 
     &  0.57713E+02, 0.58822E+02, 0.59920E+02, 0.61006E+02, 0.62083E+02, 
     &  0.63149E+02, 0.64205E+02, 0.65252E+02, 0.66291E+02, 0.67321E+02, 
     &  0.68342E+02, 0.69357E+02, 0.70363E+02, 0.71363E+02, 0.72355E+02, 
     &  0.73341E+02, 0.74320E+02, 0.75292E+02, 0.76259E+02, 0.77219E+02, 
     &  0.78174E+02, 0.79123E+02, 0.80067E+02, 0.81005E+02, 0.81939E+02, 
     &  0.82867E+02, 0.83791E+02, 0.84710E+02, 0.85624E+02, 0.86534E+02, 
     &  0.87439E+02, 0.88342E+02, 0.89239E+02, 0.90181E+02, 0.91252E+02, 
     &  0.92399E+02, 0.93540E+02, 0.94673E+02, 0.95798E+02, 0.96915E+02, 
     &  0.98027E+02, 0.99131E+02, 0.10023E+03, 0.10132E+03, 0.10240E+03, 
     &  0.10348E+03, 0.10455E+03, 0.10562E+03, 0.10668E+03, 0.10773E+03, 
     &  0.10878E+03, 0.10982E+03, 0.11086E+03, 0.11189E+03, 0.11292E+03, 
     &  0.11395E+03, 0.11497E+03, 0.11598E+03, 0.11700E+03, 0.11803E+03, 
     &  0.11905E+03, 0.12006E+03, 0.12108E+03, 0.12208E+03, 0.12308E+03, 
     &  0.12408E+03, 0.12508E+03, 0.12607E+03, 0.12705E+03, 0.12804E+03, 
     &  0.12902E+03, 0.12999E+03, 0.13096E+03, 0.13193E+03, 0.13290E+03, 
     &  0.13386E+03, 0.13482E+03, 0.13577E+03, 0.13673E+03, 0.13767E+03, 
     &  0.13862E+03, 0.13956E+03, 0.14051E+03, 0.14144E+03, 0.14239E+03, 
     &  0.14353E+03, 0.14484E+03, 0.14615E+03, 0.14745E+03, 0.14874E+03, 
     &  0.15002E+03, 0.15129E+03, 0.15256E+03, 0.15383E+03, 0.15509E+03, 
     &  0.15634E+03, 0.15758E+03, 0.15882E+03, 0.16005E+03, 0.16128E+03, 
     &  0.16250E+03, 0.16336E+03 /
	data (e0(j),j=1,         162) /
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01 /
	lam(1:nlam)=l0(1:nlam)/1d4
	R(1:nlam)=R0(1:nlam)
	Rexp(1:nlam)=e0(1:nlam)
	return
	end

	subroutine WFC3specres(lam,R,Rexp)
	IMPLICIT NONE
	integer nlam,j
	parameter(nlam=          13)
	real*8 lam(*),R(*),Rexp(*)
	real*8 l0(          13),R0(          13),e0(          13)
	data (l0(j),j=1,          13) /
     &  0.86700E+00, 0.92500E+00, 0.98300E+00, 0.10410E+01, 0.10995E+01, 
     &  0.11470E+01, 0.12165E+01, 0.12855E+01, 0.13545E+01, 0.14235E+01, 
     &  0.14925E+01, 0.15620E+01, 0.16315E+01 /
	data (R0(j),j=1,          13) /
     &  0.14948E+02, 0.15948E+02, 0.16948E+02, 0.17948E+02, 0.18636E+02, 
     &  0.16386E+02, 0.17630E+02, 0.18630E+02, 0.19630E+02, 0.20630E+02, 
     &  0.21630E+02, 0.22314E+02, 0.23645E+02 /
	data (e0(j),j=1,          13) /
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 0.20000E+01, 
     &  0.20000E+01, 0.20000E+01, 0.20000E+01 /
	lam(1:nlam)=l0(1:nlam)/1d4
	R(1:nlam)=R0(1:nlam)/2d0
	Rexp(1:nlam)=e0(1:nlam)
	return
	end


	subroutine CountSimInstrument(filename,nlamR)
	IMPLICIT NONE
	character*500 filename
	integer nlamR,ilam,i,n
	character*6000 line
	character*1000 key(100),value(100)
	
	open(unit=20,file=filename,FORM="FORMATTED",ACCESS="STREAM")

1	read(20,'(a6000)',end=3) line
	call getkeys(line,key,n)
	if(line(1:1).eq.' '.or.key(1)(1:1).eq.'#') goto 1

	ilam=0
2	read(20,'(a6000)',end=3) line
	call getkeys(line,value,n)
	if(value(1).eq.' ') goto 2
	ilam=ilam+1
	goto 2
3	close(unit=20)
	nlamR=ilam

	return
	end



	subroutine ReadSimInstrument(filename,lamR,specR,specRexp,specErr,nlamR)
	IMPLICIT NONE
	character*500 filename
	integer nlamR,ilam,i,n
	real*8 lamR(nlamR),specR(nlamR),specRexp(nlamR),specErr(nlamR)
	real*8 TRSignal,TRSNR
	character*6000 line
	character*1000 key(100),value(100)
	
	specRexp=20d0

	open(unit=20,file=filename,FORM="FORMATTED",ACCESS="STREAM")

4	read(20,'(a6000)',end=6) line
	call getkeys(line,key,n)
	if(line(1:1).eq.' '.or.line(1:1).eq.'#') goto 4

	ilam=0
5	read(20,'(a6000)',end=6) line
	call getkeys(line,value,n)
	if(value(1).eq.' ') goto 5
	ilam=ilam+1
	specErr(ilam)=-1d0
	do i=1,n
		if(value(i).ne.' ') then
			select case(key(i))
				case('Wavelength')
					read(value(i),*) lamR(ilam)
				case('BandWidth')
					read(value(i),*) specR(ilam)
				case('TRSignal')
					read(value(i),*) TRSignal
				case('TRSNR')
					read(value(i),*) TRSNR
				case('NoiseOnTransitFloorStack')
					read(value(i),*) specErr(ilam)
			end select
		endif
	enddo
	if(specErr(ilam).le.0d0) specErr(ilam)=TRSignal/TRSNR
	specR(ilam)=lamR(ilam)/(2d0*specR(ilam))
	lamR(ilam)=1d-4*lamR(ilam)
	goto 5
6	close(unit=20)

	return
	end


