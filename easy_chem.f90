!!!!---------------------------
!!!!  EASY CHEM, A CEA CLONE
!!!!---------------------------
!!!! (C) BY PAUL MOLLIERE 2015

!!!!----
!!!! SOME GLOBAL DATA STORAGE STUFF
!!!!----

module thermo_data_block
  implicit none
  CHARACTER*80, parameter, public :: fpath = '~/SPARC/Data/thermo_easy_chem_simp_own.inp'                
  INTEGER, parameter, public      :: N_coeffs = 10, N_temps = 10, N_reac_save = 1000
  INTEGER, public                 :: thermo_data_n_coeffs(N_temps,N_reac_save), thermo_data_n_intervs(N_reac_save)
  DOUBLE PRECISION, public        :: thermo_data(N_coeffs,N_temps,N_reac_save), mol_weight(N_reac_save)
  DOUBLE PRECISION, public        :: thermo_data_temps(2,N_temps,N_reac_save), reac_stoich(5,N_reac_save), &
       form_heat_Jmol_298_15_K(N_reac_save), thermo_data_T_exps(8,N_temps,N_reac_save), &
       H_0_298_15_K_m_H_0_0_K(N_temps,N_reac_save)
  DOUBLE PRECISION, parameter, public :: R = 8.3144598d0
  CHARACTER*2, public             :: reac_atoms_names(5,N_reac_save)
  LOGICAL, public                 :: reac_condensed(N_reac_save), reac_ion(N_reac_save)
  LOGICAL, public                 :: verbose, ions, quick, remove_ions
  INTEGER, public                 :: iter_max
  INTEGER, public                 :: N_gas, N_ions
  CHARACTER*40                    :: names_reactants_reorg(N_reac_save)
  DOUBLE PRECISION, public, parameter :: amu = 1.660538921d-24, kB=1.3806488d-16
  DOUBLE PRECISION, public, parameter :: mol = 6.02214129d23
  INTEGER, public, parameter      :: N_atoms_save = 104
  ! The atomic mass data shown below was taken from http://www.science.co.il/PTelements.asp
  CHARACTER*2, public, parameter  :: names_atoms_save(N_atoms_save) = &
       (/ 'E ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na', &
       'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn', &
       'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr', &
       'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb', &
       'Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd', &
       'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir', &
       'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
       'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr' /)
  DOUBLE PRECISION, public, parameter  :: masses_atoms_save(N_atoms_save) = &
       amu*(/ 0.000548579909d0, 1.0079d0,4.0026d0,6.941d0,9.0122d0,10.811d0,12.0107d0,14.0067d0 &
       ,15.9994d0,18.9984d0,20.1797d0,22.9897d0,24.305d0,26.9815d0,28.0855d0,30.9738d0,32.065d0 &
       ,35.453d0,39.948d0,39.0983d0,40.078d0,44.9559d0,47.867d0,50.9415d0,51.9961d0,54.938d0,55.845d0 &
       ,58.9332d0,58.6934d0,63.546d0,65.39d0,69.723d0,72.64d0,74.9216d0,78.96d0,79.904d0,83.8d0,85.4678d0 &
       ,87.62d0,88.9059d0,91.224d0,92.9064d0,95.94d0,98d0,101.07d0,102.9055d0,106.42d0,107.8682d0 &
       ,112.411d0,114.818d0,118.71d0,121.76d0,127.6d0,126.9045d0,131.293d0,132.9055d0,137.327d0,138.9055d0 &
       ,140.116d0,140.9077d0,144.24d0,145d0,150.36d0,151.964d0,157.25d0,158.9253d0,162.5d0,164.9303d0 &
       ,167.259d0,168.9342d0,173.04d0,174.967d0,178.49d0,180.9479d0,183.84d0,186.207d0,190.23d0,192.217d0 &
       ,195.078d0,196.9665d0,200.59d0,204.3833d0,207.2d0,208.9804d0,209d0,210d0,222d0,223d0,226d0,227d0,232.0381d0 &
       ,231.0359d0,238.0289d0,237d0,244d0,243d0,247d0,247d0,251d0,252d0,257d0,258d0,259d0,262d0/)
  
end module thermo_data_block

!!!!----
!!!! MAIN SUBROUTINE
!!!!----

subroutine EASY_CHEM(N_atoms,N_reactants,names_atoms,names_reactants,molfracs_atoms, &
     molfracs_reactants,massfracs_reactants,temp,press,ini,nabla_ad,gamma2,MMW,rho,c_pe)

  use thermo_data_block
  implicit none

  !! I/O:
  INTEGER                      :: N_atoms, N_reactants
  CHARACTER*40                 :: names_atoms(N_atoms), names_reactants(N_reactants), &
       names_reactants_orig(N_reactants)
  DOUBLE PRECISION             :: molfracs_atoms(N_atoms), molfracs_reactants(N_reactants), &
       massfracs_reactants(N_reactants)
  DOUBLE PRECISION             :: temp, press
  DOUBLE PRECISION             :: thermo_quants,nabla_ad,gamma2,MMW,rho,c_pe
  LOGICAL                      :: ini

  !! Internal:
  DOUBLE PRECISION             :: C_P_0(N_reactants), H_0(N_reactants), S_0(N_reactants), &
       molfracs_atoms_ions(N_atoms+1), temp_use
  INTEGER                      :: i_reac, N_atoms_use, gamma_neg_try
  CHARACTER*40                 :: names_atoms_ions(N_atoms+1)

  verbose = .FALSE.
  quick = .FALSE.
  remove_ions = .FALSE.
  
  ! Contains the original order of reactant names
  names_reactants_orig = names_reactants

  call init_random_seed()
  
  ! If EASY_CHEM is called the first time it needs to read in the thermodynamic data for all
  ! the considered reactant species:

  if (ini) then
     call ec_READ_ALL_DATA(N_reactants,names_reactants)
     ! Reordered such that all condensed species are at the end
     names_reactants_reorg(1:N_reactants) = names_reactants
  else
     ! ec_READ_ALL_DATA is only called when ini == .TRUE.
     ! (must be done at least once!).
     ! For ini == .FALSE. restore the ordering where the
     ! condensed species are at the end, because this
     ! is how the thermodynamic data is stored
     names_reactants = names_reactants_reorg(1:N_reactants)
  end if

  names_atoms_ions(1:N_atoms) = names_atoms
  molfracs_atoms_ions(1:N_atoms) = molfracs_atoms
  names_atoms_ions(N_atoms+1) = 'E'
  molfracs_atoms_ions(N_atoms+1) = 0d0
  
  IF (ions .AND. temp > 750d0) THEN
     N_atoms_use = N_atoms+1
  ELSE IF (ions .AND. temp <= 750d0) THEN
     remove_ions = .TRUE.
     N_atoms_use = N_atoms
  ELSE
     N_atoms_use = N_atoms
  END IF
  
  call ec_CALC_THERMO_QUANTS(N_reactants,names_reactants,temp,C_P_0, H_0, S_0)
  gamma2 = 0d0
  temp_use = temp
  gamma_neg_try = 0d0
  do while (gamma2 < 1d0)
     call ec_CALC_EQU_CHEM(N_atoms_use,N_reactants,names_atoms_ions(1:N_atoms_use), &
          names_reactants,molfracs_atoms_ions(1:N_atoms_use), &
          molfracs_reactants,massfracs_reactants,temp_use,press,C_P_0, H_0, S_0, &
          names_reactants_orig,nabla_ad,gamma2,MMW,rho,c_pe)
     if (gamma2 < 1d0) then
        write(*,*) 'Gamma was < 1, redo! gamma2, temp, ', gamma2, temp
        gamma_neg_try = gamma_neg_try + 1
        if (gamma_neg_try > 10) then
           call random_number(temp_use)
           temp_use = temp*(1d0 + 0.01d0*temp_use)
           write(*,*) 'temp, temp_use', temp, temp_use
           call ec_CALC_THERMO_QUANTS(N_reactants,names_reactants,temp_use,C_P_0, H_0, S_0)
        end if
     end if
  end do

  ! Restore names etc. to original, i.e. user input order.
  names_reactants = names_reactants_orig

  c_pe = c_pe*1d7 ! J/(g K) to erg/(g K)
  
end subroutine EASY_CHEM

!!!!----
!!!! "SUB"-SUBROUTINES
!!!!----

!#####################################################

subroutine ec_READ_ALL_DATA(N_reactants,names_reactants)

  use thermo_data_block
  implicit none
  !! I/O:
  INTEGER                      :: N_reactants
  CHARACTER*40                 :: names_reactants(N_reactants)

  !! FILE I/O:
  CHARACTER*80                 :: file_line

  !! INTERNAL
  INTEGER                      :: i_reac, n_interv, i_interv, &
       i_stoich, stoich_start, reac_found

  reac_found = 0
  thermo_data_n_intervs = -1
  N_gas = 0
  N_ions = 0

  reac_ion = .FALSE.
  ions = .FALSE.

  OPEN(unit=17,file=fpath)
  DO WHILE (1>0)
     READ(17,'(A80)',end=122) file_line
     DO i_reac = 1, N_reactants
        IF (TRIM(ADJUSTL(file_line(1:18))) .EQ. &
             TRIM(ADJUSTL(names_reactants(i_reac)))) THEN
           reac_found = reac_found+1
           READ(17,'(A80)',end=122) file_line
           READ(file_line(1:3),'(I2)') n_interv
           thermo_data_n_intervs(i_reac) = n_interv
           IF (file_line(52:52) == '0') THEN
              reac_condensed(i_reac) = .FALSE.
              N_gas = N_gas + 1
           ELSE
              reac_condensed(i_reac) = .TRUE.
           END IF
        END IF
     END DO
  END DO
122 close(17)

  IF (reac_found .NE. N_reactants) THEN
     WRITE(*,*) 
     WRITE(*,*) 'EASY CHEM ERROR! For the following species no thermodynamical data was found:'
     DO i_reac = 1, N_reactants
        IF (thermo_data_n_intervs(i_reac) .EQ. -1) THEN
           WRITE(*,*) trim(adjustl(names_reactants(i_reac)))
        END IF
     END DO
     STOP
  END IF

  ! Put condensed species at the end of the reactants array
  call reorder_specs(N_reactants,names_reactants)
  
  ! BASED ON THE ~RIGHT DESCRIPTION GIVEN IN GORDON 1996, page 73
  ! AND THE APPEARANCE OF THERMO.INP.
  OPEN(unit=17,file=fpath)
  DO WHILE (1>0)
     READ(17,'(A80)',end=123) file_line
     DO i_reac = 1, N_reactants
        IF (TRIM(ADJUSTL(file_line(1:18))) .EQ. &
             TRIM(ADJUSTL(names_reactants(i_reac)))) THEN
           READ(17,'(A80)',end=123) file_line
           READ(file_line(1:3),'(I2)') n_interv
           thermo_data_n_intervs(i_reac) = n_interv
           stoich_start = 11
           DO i_stoich = 1, 5
              reac_atoms_names(i_stoich,i_reac) = file_line(stoich_start:stoich_start+1)
              ! Are there ions to be treated?
              IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                 ions = .TRUE.
                 reac_ion(i_reac) = .TRUE.
                 N_ions = N_ions + 1
              END IF
              READ(file_line(stoich_start+2:stoich_start+7),'(F6.2)') reac_stoich(i_stoich,i_reac)
              stoich_start = stoich_start+8
           END DO
           IF (file_line(52:52) == '0') THEN
              reac_condensed(i_reac) = .FALSE.
           ELSE
              reac_condensed(i_reac) = .TRUE.
           END IF
           READ(file_line(53:65),'(F13.5)') mol_weight(i_reac)
           READ(file_line(66:80),'(F13.5)') form_heat_Jmol_298_15_K(i_reac)
           DO i_interv = 1, 3*n_interv
              READ(17,'(A80)',end=123) file_line
              IF (MOD(i_interv,3) .EQ. 1) THEN
                 READ(file_line(2:22),'(F10.3,X,F10.3)') thermo_data_temps(1,i_interv/3+1,i_reac), &
                      thermo_data_temps(2,i_interv/3+1,i_reac)
                 READ(file_line(23:23),'(I1)') thermo_data_n_coeffs(i_interv/3+1,i_reac)
                 READ(file_line(24:63),'(8F5.1)') thermo_data_T_exps(1,i_interv/3+1,i_reac), &
                      thermo_data_T_exps(2,i_interv/3+1,i_reac), thermo_data_T_exps(3,i_interv/3+1,i_reac), &
                      thermo_data_T_exps(4,i_interv/3+1,i_reac), thermo_data_T_exps(5,i_interv/3+1,i_reac), &
                      thermo_data_T_exps(6,i_interv/3+1,i_reac), thermo_data_T_exps(7,i_interv/3+1,i_reac), &
                      thermo_data_T_exps(8,i_interv/3+1,i_reac)
                 READ(file_line(66:80),'(F15.3)') H_0_298_15_K_m_H_0_0_K(i_interv/3+1,i_reac)
              END IF
              IF (MOD(i_interv,3) .EQ. 2) THEN
                 READ(file_line(1:80),'(5D16.8)') thermo_data(1,i_interv/3+1,i_reac), &
                      thermo_data(2,i_interv/3+1,i_reac), thermo_data(3,i_interv/3+1,i_reac) , &
                      thermo_data(4,i_interv/3+1,i_reac), thermo_data(5,i_interv/3+1,i_reac)
              END IF
              IF (MOD(i_interv,3) .EQ. 0) THEN
                 READ(file_line(1:80),'(5D16.8)') thermo_data(6,i_interv/3,i_reac), &
                      thermo_data(7,i_interv/3,i_reac), thermo_data(8,i_interv/3,i_reac) , &
                      thermo_data(9,i_interv/3,i_reac), thermo_data(10,i_interv/3,i_reac)
              END IF
           END DO
        END IF
     END DO
  END DO
123 CLOSE(17)

end subroutine ec_READ_ALL_DATA

!#####################################################

subroutine reorder_specs(N_reactants,names_reactants)

  use thermo_data_block
  implicit none
  !! I/O:
  INTEGER                      :: N_reactants
  CHARACTER*40                 :: names_reactants(N_reactants), &
       names_reactants_buff(N_reactants)
  !! Internal:
  INTEGER                      :: i_reac, gas_offset, cond_offset
  
  names_reactants_buff = names_reactants
  gas_offset = 1
  cond_offset = 1
  DO i_reac = 1, N_reactants
     IF (reac_condensed(i_reac)) THEN
        names_reactants(N_gas+cond_offset) = names_reactants_buff(i_reac)
        cond_offset = cond_offset + 1
     ELSE
        names_reactants(gas_offset) = names_reactants_buff(i_reac)
        gas_offset = gas_offset + 1
     END IF
  END DO 

end subroutine reorder_specs

!#####################################################

subroutine ec_CALC_THERMO_QUANTS(N_reactants,names_reactants,temp,C_P_0, H_0, S_0)

  use thermo_data_block
  implicit none
  !! I/O
  INTEGER                      :: N_reactants
  CHARACTER*40                 :: names_reactants(N_reactants)
  DOUBLE PRECISION             :: temp
  DOUBLE PRECISION             :: C_P_0(N_reactants), H_0(N_reactants), &
       S_0(N_reactants)
  !! internal
  INTEGER                      :: i_reac, i_tempinv, tempinv_ind

  DO i_reac = 1, N_reactants

     ! Get temperature interpolation range
     IF (temp < thermo_data_temps(1,1,i_reac)) THEN
        tempinv_ind = 1
     ELSE IF (temp >= thermo_data_temps(2,thermo_data_n_intervs(i_reac),i_reac)) THEN
        tempinv_ind = thermo_data_n_intervs(i_reac)
     ELSE
        DO i_tempinv = 1, thermo_data_n_intervs(i_reac)
           IF ((temp >= thermo_data_temps(1,i_tempinv,i_reac)) .AND. &
                (temp < thermo_data_temps(2,i_tempinv,i_reac))) THEN
              tempinv_ind = i_tempinv
              EXIT
           END IF
        END DO
     END IF

     ! Calculate thermodynamic quantities as explained in Gordon 1996, page 74
     C_P_0(i_reac) = (thermo_data(1,tempinv_ind,i_reac)*temp**(-2)+ &
          thermo_data(2,tempinv_ind,i_reac)*temp**(-1)+ &
          thermo_data(3,tempinv_ind,i_reac)+thermo_data(4,tempinv_ind,i_reac)* &
          temp**(1)+thermo_data(5,tempinv_ind,i_reac)*temp**(2)+ &
          thermo_data(6,tempinv_ind,i_reac)*temp**(3)+thermo_data(7,tempinv_ind,i_reac)* &
          temp**(4))*R
     H_0(i_reac) = (-thermo_data(1,tempinv_ind,i_reac)*temp**(-2)+ &
          thermo_data(2,tempinv_ind,i_reac)*temp**(-1)*log(temp)+ &
          thermo_data(3,tempinv_ind,i_reac)+thermo_data(4,tempinv_ind,i_reac)*temp**(1)/2d0+ &
          thermo_data(5,tempinv_ind,i_reac)*temp**(2)/3d0+ &
          thermo_data(6,tempinv_ind,i_reac)*temp**(3)/4d0+thermo_data(7,tempinv_ind,i_reac)* &
          temp**(4)/5d0+thermo_data(9,tempinv_ind,i_reac)/temp)* &
          R*temp
     S_0(i_reac) = (-thermo_data(1,tempinv_ind,i_reac)*temp**(-2)/2d0- &
          thermo_data(2,tempinv_ind,i_reac)*temp**(-1)+ &
          thermo_data(3,tempinv_ind,i_reac)*log(temp)+ &
          thermo_data(4,tempinv_ind,i_reac)*temp**(1)+ &
          thermo_data(5,tempinv_ind,i_reac)*temp**(2)/2d0+ &
          thermo_data(6,tempinv_ind,i_reac)*temp**(3)/3d0+thermo_data(7,tempinv_ind,i_reac)* &
          temp**(4)/4d0+thermo_data(10,tempinv_ind,i_reac))*R
     
  END DO
  
end subroutine ec_CALC_THERMO_QUANTS

!#####################################################

recursive subroutine ec_CALC_EQU_CHEM(N_atoms,N_reactants,names_atoms,names_reactants,molfracs_atoms, &
     molfracs_reactants,massfracs_reactants,temp,press,C_P_0, H_0, S_0,names_reactants_orig, &
     nabla_ad,gamma2,MMW,rho,c_pe)

  use thermo_data_block
  implicit none

  !! I/O:
  INTEGER                      :: N_atoms, N_reactants
  CHARACTER*40                 :: names_atoms(N_atoms), names_reactants(N_reactants), &
       names_reactants_orig(N_reactants)
  DOUBLE PRECISION             :: molfracs_atoms(N_atoms), molfracs_reactants(N_reactants), &
       massfracs_reactants(N_reactants)
  DOUBLE PRECISION             :: temp, press
  DOUBLE PRECISION             :: C_P_0(N_reactants), H_0(N_reactants), S_0(N_reactants)

  !! CEA McBride 1994 style variables:
  DOUBLE PRECISION             :: n ! Moles of gas particles per total mass of mixture in kg
  DOUBLE PRECISION             :: n_spec(N_reactants) ! Moles of species per total mass of mixture in kg
  DOUBLE PRECISION             :: n_spec_old(N_reactants) ! Moles of species per total mass of mixture in kg
                                                          ! of previous iteration
  DOUBLE PRECISION             :: pi_atom(N_atoms) ! Lagrangian multipliers for the atomic species divided
                                                   ! by (R*T)
  DOUBLE PRECISION             :: matrix(N_reactants+N_atoms+1,N_reactants+N_atoms+1)
                                 ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
                                 ! condensed species, the pis and the delta log(n)
  DOUBLE PRECISION             :: vector(N_reactants+N_atoms+1), &
       solution_vector(N_reactants+N_atoms+1), nabla_ad,gamma2,MMW,rho
  
  !! Internal:
  INTEGER                       :: i_iter, i_reac, j_reac, inc_next, &
       current_solids_number, N_spec_eff, buffer_ind, i_atom
  LOGICAL                       :: converged, remove_cond, slowed
  LOGICAL, allocatable          :: solid_inc(:), neg_cond(:)
  DOUBLE PRECISION, allocatable :: dgdnj(:)
  INTEGER, allocatable          :: solid_indices(:), solid_indices_buff(:)
  DOUBLE PRECISION              :: nsum, mu_gas(N_gas), a_gas(N_gas,N_atoms), mass_species, atom_mass, &
       msum, c_pe

  converged = .FALSE.
  slowed = .FALSE.
  call ec_MAKE_ini_vals(N_reactants,N_atoms,n,n_spec,pi_atom)

  iter_max = 50000 + N_reactants/2
  current_solids_number = 0

  MMW = 0d0

  n_spec_old = n_spec

  ! FIRST: DO GAS ONLY!
  DO i_iter = 1, iter_max

     IF (quick) THEN
        call ec_MAKE_MATRIX_short(N_atoms,names_atoms,molfracs_atoms,N_gas,press,temp,&
             C_P_0,H_0,S_0,n,n_spec,pi_atom,matrix(1:N_atoms+1,1:N_atoms+1),vector(1:N_atoms+1),(/1,1,1,1,1/),&
             names_reactants, N_reactants, 5, mu_gas,a_gas)
        call ec_INV_MATRIX_short(N_atoms+1, &
             matrix(1:N_atoms+1,1:N_atoms+1),vector(1:N_atoms+1), &
             solution_vector(1:N_atoms+1))
        call ec_CHANGE_ABUNDS_short(N_atoms,N_gas,solution_vector(1:N_atoms+1),n_spec,pi_atom,n,converged,&
             (/1,1,1,1,1/),5,mu_gas,a_gas,temp,names_atoms,molfracs_atoms,N_reactants,n_spec_old)
     ELSE
        call ec_MAKE_MATRIX_long(N_atoms,names_atoms,molfracs_atoms,N_gas,&
             press,temp,C_P_0,H_0,S_0,n,n_spec,pi_atom, &
             matrix(1:N_gas+N_atoms+1,1:N_gas+N_atoms+1),vector(1:N_gas+N_atoms+1),&
             (/1,1,1,1,1/),names_reactants,N_reactants,5)
        call ec_INV_MATRIX_long(N_atoms+N_gas+1, &
             matrix(1:N_gas+N_atoms+1,1:N_gas+N_atoms+1),vector(1:N_gas+N_atoms+1), &
             solution_vector(1:N_gas+N_atoms+1))
        call ec_CHANGE_ABUNDS_long(N_atoms,N_gas,solution_vector(1:N_gas+N_atoms+1), &
             n_spec,pi_atom,n,converged,(/1,1,1,1,1/),5,names_atoms,molfracs_atoms,N_reactants,&
             n_spec_old)
     END IF

     n_spec_old = n_spec
     
     IF (verbose) THEN
        write(*,*)
        write(*,*)
        write(*,*) i_iter
        DO i_reac = 1, N_reactants
           write(*,*) names_reactants(i_reac), n_spec(i_reac)/SUM(n_spec)
        END DO
     END IF

     IF (converged) THEN
        EXIT
     END IF
  END DO

  IF (.NOT. converged) THEN
     WRITE(*,*) 
     WRITE(*,*) 'EASY CHEM WARNING: One ore more convergence criteria not satisfied!'
  END IF

  converged = .FALSE.
  remove_cond = .FALSE.

  IF (N_gas .EQ. N_reactants) THEN
     call ec_CALC_ADIABATIC_GRADIENT(N_atoms,N_gas,N_reactants,n_spec, &
          n,H_0,C_P_0,(/ 1,1,1,1,1 /),5,temp,names_atoms,nabla_ad,gamma2,c_pe)
  END IF
  
  ! THEN: INCLUDE CONDENSATES!
  IF (N_gas < N_reactants) THEN
     allocate(solid_inc(N_reactants-N_gas))
     allocate(dgdnj(N_reactants-N_gas))
     allocate(solid_indices(N_reactants-N_gas))
     allocate(solid_indices_buff(N_reactants-N_gas))
     allocate(neg_cond(N_reactants-N_gas))
     solid_inc = .FALSE.
     inc_next = 0
     neg_cond = .FALSE.

     N_spec_eff = N_gas
     
     DO WHILE (.NOT. (inc_next .EQ. -1))

        call ec_INCLUDE_SOLIDS_QUESTIONMARK(N_reactants,N_atoms,pi_atom,names_atoms,H_0,S_0,dgdnj,temp, &
             solid_inc,inc_next,neg_cond)

        IF (.NOT. (inc_next .EQ. 0)) THEN
           DO i_reac = N_gas+1, N_reactants
              IF (n_spec(i_reac) < 0d0) THEN
                 n_spec(i_reac) = 0d0
                 inc_next = i_reac
                 remove_cond = .TRUE.
              END IF
           END DO
        END IF
        
        IF (.NOT. (inc_next .EQ. -1)) THEN

           IF (remove_cond) THEN
              current_solids_number = current_solids_number - 1
              solid_indices_buff = 0
              buffer_ind = 1
              DO i_reac = 1, N_reactants-N_gas
                 IF (solid_indices(i_reac) .NE. inc_next) THEN
                    solid_indices_buff(buffer_ind) = solid_indices(i_reac)
                    buffer_ind = buffer_ind + 1
                 END IF
              END DO
              solid_indices = solid_indices_buff
              solid_inc(inc_next-N_gas) = .FALSE.
              neg_cond(inc_next-N_gas) = .TRUE.
           ELSE
              current_solids_number = current_solids_number + 1
              solid_indices(current_solids_number) = inc_next
              solid_inc(inc_next-N_gas) = .TRUE.
           END IF
           
           N_spec_eff = N_gas+current_solids_number
           DO i_iter = 1, iter_max

              IF (quick) THEN
                 call ec_MAKE_MATRIX_short(N_atoms,names_atoms,molfracs_atoms,N_spec_eff,press,temp,&
                      C_P_0,H_0,S_0,n,n_spec,pi_atom, &
                      matrix(1:N_atoms+1+N_spec_eff-N_gas,1:N_atoms+1+N_spec_eff-N_gas), &
                      vector(1:N_atoms+1+N_spec_eff-N_gas),solid_indices,&
                      names_reactants, N_reactants, N_spec_eff-N_gas, mu_gas,a_gas)
                 call ec_INV_MATRIX_short(N_atoms+1+N_spec_eff-N_gas, &
                      matrix(1:N_atoms+1+N_spec_eff-N_gas,1:N_atoms+1+N_spec_eff-N_gas), &
                      vector(1:N_atoms+1+N_spec_eff-N_gas), &
                      solution_vector(1:N_atoms+1+N_spec_eff-N_gas))
                 call ec_CHANGE_ABUNDS_short(N_atoms,N_spec_eff,solution_vector(1:N_atoms+1+N_spec_eff-N_gas), &
                      n_spec,pi_atom,n,converged,&
                      solid_indices,N_spec_eff-N_gas,mu_gas,a_gas,temp,names_atoms,molfracs_atoms,N_reactants, &
                      n_spec_old)
              ELSE
                 call ec_MAKE_MATRIX_long(N_atoms,names_atoms,molfracs_atoms,N_spec_eff,&
                      press,temp,C_P_0,H_0,S_0,n,n_spec,pi_atom, &
                      matrix(1:N_spec_eff+N_atoms+1,1:N_spec_eff+N_atoms+1),vector(1:N_spec_eff+N_atoms+1),&
                      solid_indices,names_reactants,N_reactants,N_spec_eff-N_gas)
                 call ec_INV_MATRIX_long(N_atoms+N_spec_eff+1, &
                      matrix(1:N_spec_eff+N_atoms+1,1:N_spec_eff+N_atoms+1),vector(1:N_spec_eff+N_atoms+1), &
                      solution_vector(1:N_spec_eff+N_atoms+1))
                 call ec_CHANGE_ABUNDS_long(N_atoms,N_spec_eff,solution_vector(1:N_spec_eff+N_atoms+1), &
                      n_spec,pi_atom,n,converged,&
                      solid_indices,N_spec_eff-N_gas,names_atoms,molfracs_atoms,N_reactants, &
                      n_spec_old)
              END IF

              n_spec_old = n_spec
              
              IF (verbose) THEN
                 write(*,*)
                 write(*,*)
                 write(*,*) i_iter
                 DO i_reac = 1, N_reactants
                    write(*,*) names_reactants(i_reac), n_spec(i_reac)/SUM(n_spec)
                 END DO
              END IF
              
              DO i_reac = N_gas+1, N_reactants
                 IF ((n_spec(i_reac) < 0d0) .AND. (i_iter > 30)) THEN
                    converged = .TRUE.
                    EXIT
                 END IF
              END DO

              IF (converged) THEN
                 EXIT
              END IF
              
           END DO

           IF (.NOT. converged) THEN
              WRITE(*,*)
              IF (quick) THEN
                 quick = .FALSE.
                 call ec_CALC_EQU_CHEM(N_atoms,N_reactants,names_atoms, &
                      names_reactants,molfracs_atoms, &
                      molfracs_reactants,massfracs_reactants, &
                      temp,press,C_P_0, H_0, S_0,names_reactants_orig, &
                      nabla_ad,gamma2,MMW,rho,c_pe)
                 quick = .TRUE.
                 slowed = .TRUE.
                 if (verbose) THEN
                    WRITE(*,*)
                    write(*,*) 'SLOW!'
                 END IF
                 EXIT
              ELSE
                 WRITE(*,*)
                 WRITE(*,*) 'EASY CHEM WARNING: One ore more convergence criteria not satisfied!'
              END IF
           END IF

           converged = .FALSE.
           
        END IF

        remove_cond = .FALSE.
        
     END DO

     ! Calc. nabla_ad
     IF (.NOT. slowed) THEN
        call ec_CALC_ADIABATIC_GRADIENT(N_atoms,N_spec_eff,N_reactants,n_spec, &
             n,H_0,C_P_0,solid_indices,N_spec_eff-N_gas,temp,names_atoms,nabla_ad,gamma2,c_pe)
     END IF
     
     deallocate(solid_inc)
     deallocate(dgdnj)
     deallocate(solid_indices)
     deallocate(solid_indices_buff)
     deallocate(neg_cond)
     
  END IF

  ! PREPARE FINAL OUTPUT
  
  IF (.NOT. slowed) THEN
     nsum = SUM(n_spec)
     DO i_reac = 1, N_reactants
        IF (n_spec(i_reac)/nsum < 1d-50) THEN
           n_spec(i_reac) = 0d0
        END IF
     END DO

     nsum = SUM(n_spec)
     DO i_reac = 1, N_reactants
        DO j_reac = 1, N_reactants
           IF (names_reactants(j_reac) .EQ. names_reactants_orig(i_reac)) THEN
              molfracs_reactants(i_reac) = n_spec(j_reac)/nsum
              EXIT
           END IF
        END DO
     END DO

     DO i_reac = 1, N_reactants
        mass_species = 0d0
        DO i_atom = 1, 5
           if (trim(adjustl(reac_atoms_names(i_atom,i_reac))) .NE. '') THEN
              call ec_ATOM_MASS(reac_atoms_names(i_atom,i_reac),atom_mass)
              mass_species = mass_species+atom_mass*DBLE(reac_stoich(i_atom,i_reac))
           END IF
        END DO
        DO j_reac = 1, N_reactants
           IF (names_reactants(i_reac) .EQ. names_reactants_orig(j_reac)) THEN
              massfracs_reactants(j_reac) = n_spec(i_reac)*mass_species
              MMW = MMW + massfracs_reactants(j_reac)/mass_species
              EXIT
           END IF
        END DO
     END DO
     msum = SUM(massfracs_reactants)
     massfracs_reactants = massfracs_reactants / msum

     MMW = MMW/msum
     MMW = 1d0/MMW

     rho = (press*1d6)/kB/temp*MMW
     MMW = MMW/amu
     
  END IF

end subroutine ec_CALC_EQU_CHEM

!#####################################################

subroutine ec_MAKE_ini_vals(N_reactants,N_atoms,n,n_spec,pi_atom)

  use thermo_data_block       
  implicit none
  !! I/O:
  INTEGER                      :: N_atoms, N_reactants
  DOUBLE PRECISION             :: n ! Moles of gas particles per total mass of mixture in kg
  DOUBLE PRECISION             :: n_spec(N_reactants) ! Moles of species per total mass of mixture in kg
  DOUBLE PRECISION             :: pi_atom(N_atoms) ! Lagrangian multipliers for the atomic species divided
                                                   ! by (R*T)

  !! Internal:
  INTEGER                      :: i_reac

  n = 0.1d0
  n_spec = 0d0
  pi_atom = 0d0
  DO i_reac = 1, N_gas
     n_spec(i_reac) = n/DBLE(N_gas)
     IF (remove_ions) THEN
        IF(reac_ion(i_reac)) THEN
           n_spec(i_reac) = 0d0
        END IF
     END IF
  END DO

end subroutine ec_MAKE_ini_vals

!#####################################################

subroutine ec_MAKE_MATRIX_long(N_atoms,names_atoms,molfracs_atoms,N_reactants,press,temp,&
     C_P_0,H_0,S_0,n,n_spec,pi_atom,matrix,vector,solid_indices,&
     names_reactants, N_reactants2, N_solids)

  use thermo_data_block       
  implicit none
  !! I/O:
  INTEGER                      :: N_atoms, N_reactants, N_reactants2, N_solids
  INTEGER                      :: solid_indices(N_solids)
  CHARACTER*40                 :: names_atoms(N_atoms), names_reactants(N_reactants2)
  DOUBLE PRECISION             :: molfracs_atoms(N_atoms), press, temp
  DOUBLE PRECISION             :: C_P_0(N_reactants2), H_0(N_reactants2), S_0(N_reactants2)
  DOUBLE PRECISION             :: n ! Moles of gas particles per total mass of mixture in kg
  DOUBLE PRECISION             :: n_spec(N_reactants2) ! Moles of species per total mass of mixture in kg
  DOUBLE PRECISION             :: pi_atom(N_atoms) ! Lagrangian multipliers for the atomic species divided
                                                   ! by (R*T)
  DOUBLE PRECISION             :: matrix(N_reactants+N_atoms+1,N_reactants+N_atoms+1)
                                 ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
                                 ! condensed species, the pis and the delta log(n)
  DOUBLE PRECISION             :: vector(N_reactants+N_atoms+1)
  
  !! Internal:
  DOUBLE PRECISION             :: b_0(N_atoms), b_0_norm, mass_atom, b(N_atoms)
  DOUBLE PRECISION             :: a(N_reactants,N_atoms), mu(N_reactants)
  INTEGER                      :: i_atom, i_reac, i_ratom
  CHARACTER*40                 :: upper_atom_name
  CHARACTER*2                  :: upper_ratom_name

!!$  write(*,*) N_reactants

  ! Set up b0
  b_0_norm = 0d0
  DO i_atom = 1, N_atoms
     call ec_ATOM_MASS(names_atoms(i_atom),mass_atom)
     b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
  END DO
  b_0 = molfracs_atoms/b_0_norm

  ! Set up a_ij
  a = 0d0
  DO i_atom = 1, N_atoms
     call To_upper(names_atoms(i_atom),upper_atom_name)
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              a(i_reac,1:N_atoms) = 0d0
              CYCLE
           END IF
        END IF
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
           END IF
        END DO
     END DO
     DO i_reac = N_gas+1, N_reactants
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,solid_indices(i_reac-N_gas)),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
           END IF
        END DO
     END DO
  END DO

  ! Set up mu_j
  DO i_reac = 1, N_reactants
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           mu(i_reac) = 0d0
           CYCLE
        END IF
     END IF
     ! Taken from Venot et al. (2012), in comparison with McBride 1996.
     IF (i_reac <= N_gas) THEN
        mu(i_reac) = H_0(i_reac) - temp*S_0(i_reac)

        IF (n_spec(i_reac) > 1d-290) THEN
           mu(i_reac) = mu(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*log(press)
        ELSE
           IF (verbose) THEN
              write(*,*) 'n_spec(i_reac) == 0 for '//trim(adjustl(names_reactants(i_reac)))// &
                   ' set to 1d-13 and try again.'
           END IF
           call RANDOM_NUMBER(n_spec(i_reac))
           n_spec(i_reac) = n_spec(i_reac)*1d-13
           mu(i_reac) = mu(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*log(press)
        END IF
        
     ELSE
        mu(i_reac) = H_0(solid_indices(i_reac-N_gas)) - temp*S_0(solid_indices(i_reac-N_gas))
     END IF
  END DO

  ! MATRIX SETUP
  matrix = 0d0
  ! Set up the matrix for the N_gas equations (Eq. 2.18)
  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     matrix(i_reac,i_reac) = 1d0
     DO i_atom = 1, N_atoms
        matrix(i_reac,N_reactants+i_atom) = -a(i_reac,i_atom)
     END DO
     matrix(i_reac,N_reactants+N_atoms+1) = -1d0
  END DO

  ! Set up the matrix for the N_reactants-N_gas equations (Eq. 2.19)
  IF (N_gas < N_reactants) THEN
     DO i_reac = N_gas+1, N_reactants
        DO i_atom = 1, N_atoms
           matrix(i_reac,N_reactants+i_atom) = -a(i_reac,i_atom)
        END DO
     END DO
  END IF

  ! Set up the matrix for the N_atom equations (Eq. 2.20)
  DO i_atom = 1, N_atoms
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              CYCLE
           END IF
        END IF
        matrix(N_reactants+i_atom,i_reac) = a(i_reac,i_atom)*n_spec(i_reac)
     END DO
     IF (N_gas < N_reactants) THEN
        DO i_reac = N_gas+1, N_reactants
           matrix(N_reactants+i_atom,i_reac) = a(i_reac,i_atom)
        END DO
     END IF
  END DO

  ! Set up the matrix for the last equation (Eq. 2.21)
  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     matrix(N_reactants+N_atoms+1,i_reac) = n_spec(i_reac)
  END DO
  matrix(N_reactants+N_atoms+1,N_reactants+N_atoms+1) = -n

  ! VECTOR SETUP
  !vector(N_reactants+N_atoms+1)
  vector = 0d0

  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     vector(i_reac) = -mu(i_reac)/R/temp ! (Eq. 2.18)
  END DO
  
  IF (N_gas < N_reactants) THEN
     vector(N_gas+1:N_reactants) = -mu(N_gas+1:N_reactants)/R/temp ! (Eq. 2.19)
  END IF

  b = 0d0
  DO i_atom = 1, N_atoms
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              CYCLE
           END IF
        END IF
        b(i_atom) = b(i_atom) + a(i_reac,i_atom)*n_spec(i_reac)
     END DO
     DO i_reac = N_gas+1, N_reactants
        b(i_atom) = b(i_atom) + a(i_reac,i_atom)*n_spec(solid_indices(i_reac-N_gas))
     END DO
  END DO
  vector(N_reactants+1:N_reactants+N_atoms) = b_0 - b ! (Eq. 2.20)

  vector(N_reactants+N_atoms+1) = n - SUM(n_spec(1:N_gas)) ! (Eq. 2.21)
  
end subroutine ec_MAKE_MATRIX_long

!#####################################################

subroutine ec_MAKE_MATRIX_short(N_atoms,names_atoms,molfracs_atoms,N_reactants,press,temp,&
     C_P_0,H_0,S_0,n,n_spec,pi_atom,matrix,vector,solid_indices,&
     names_reactants, N_reactants2, N_solids, mu_gas,a_gas)

  use thermo_data_block       
  implicit none
  !! I/O:
  INTEGER                      :: N_atoms, N_reactants, N_reactants2, N_solids
  INTEGER                      :: solid_indices(N_solids)
  CHARACTER*40                 :: names_atoms(N_atoms), names_reactants(N_reactants2)
  DOUBLE PRECISION             :: molfracs_atoms(N_atoms), press, temp
  DOUBLE PRECISION             :: C_P_0(N_reactants2), H_0(N_reactants2), S_0(N_reactants2)
  DOUBLE PRECISION             :: n ! Moles of gas particles per total mass of mixture in kg
  DOUBLE PRECISION             :: n_spec(N_reactants2) ! Moles of species per total mass of mixture in kg
  DOUBLE PRECISION             :: pi_atom(N_atoms) ! Lagrangian multipliers for the atomic species divided
                                                   ! by (R*T)
  DOUBLE PRECISION             :: matrix(N_atoms+1+(N_reactants-N_gas),N_atoms+1+(N_reactants-N_gas))
                                 ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
                                 ! condensed species, the pis and the delta log(n)
  DOUBLE PRECISION             :: vector(N_atoms+1+(N_reactants-N_gas))
  DOUBLE PRECISION             :: mu_gas(N_gas), a_gas(N_gas,N_atoms)
  
  !! Internal:
  DOUBLE PRECISION             :: b_0(N_atoms), b_0_norm, mass_atom, b(N_atoms)
  DOUBLE PRECISION             :: a(N_reactants,N_atoms), mu(N_reactants)
  INTEGER                      :: i_atom, i_reac, i_ratom, i_atom2
  CHARACTER*40                 :: upper_atom_name
  CHARACTER*2                  :: upper_ratom_name

  ! Set up b0
  b_0_norm = 0d0
  DO i_atom = 1, N_atoms
     call ec_ATOM_MASS(names_atoms(i_atom),mass_atom)
     b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
  END DO
  b_0 = molfracs_atoms/b_0_norm

  ! Set up a_ij
  a = 0d0
  DO i_atom = 1, N_atoms
     call To_upper(names_atoms(i_atom),upper_atom_name)
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              a(i_reac,1:N_atoms) = 0d0
              CYCLE
           END IF
        END IF
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
           END IF
        END DO
     END DO
     DO i_reac = N_gas+1, N_reactants
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,solid_indices(i_reac-N_gas)),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
           END IF
        END DO
     END DO
  END DO

  ! Set up mu_j
  DO i_reac = 1, N_reactants
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           mu(i_reac) = 0d0
           CYCLE
        END IF
     END IF
     ! Taken from Venot et al. (2012), in comparison with McBride 1996.
     IF (i_reac <= N_gas) THEN
        mu(i_reac) = H_0(i_reac) - temp*S_0(i_reac)

        IF (n_spec(i_reac) > 1d-290) THEN
           mu(i_reac) = mu(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*log(press)
        ELSE
           IF (verbose) THEN
              write(*,*) 'n_spec(i_reac) == 0 for '//trim(adjustl(names_reactants(i_reac)))// &
                   ' set to 1d-13 and try again.'
           END IF
           call RANDOM_NUMBER(n_spec(i_reac))
           n_spec(i_reac) = n_spec(i_reac)*1d-13
           mu(i_reac) = mu(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*log(press)
        END IF        
     ELSE
        mu(i_reac) = H_0(solid_indices(i_reac-N_gas)) - temp*S_0(solid_indices(i_reac-N_gas))
     END IF
  END DO

  a_gas = a(1:N_gas,1:N_atoms)
  mu_gas = mu(1:N_gas)

  ! MATRIX SETUP
  matrix = 0d0

  ! Set up the matrix for the N_atoms equations (Eq. 2.24)
  DO i_atom = 1, N_atoms
     DO i_atom2 = 1, N_atoms
        DO i_reac = 1, N_gas
           IF (remove_ions) THEN
              IF (reac_ion(i_reac)) THEN
                 CYCLE
              END IF
           END IF
           matrix(i_atom,i_atom2) = matrix(i_atom,i_atom2) + &
                a(i_reac,i_atom)*a(i_reac,i_atom2)*n_spec(i_reac)
        END DO
     END DO

     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              CYCLE
           END IF
        END IF
        matrix(i_atom,N_atoms+1) = matrix(i_atom,N_atoms+1) + &
             a(i_reac,i_atom)*n_spec(i_reac)
     END DO

     IF (N_gas < N_reactants) THEN
        DO i_reac = N_gas+1, N_reactants
           matrix(i_atom,N_atoms+1+i_reac-N_gas) = a(i_reac,i_atom)
        END DO
     END IF
     
  END DO

  ! Set up the matrix for the equation (Eq. 2.26)
  DO i_atom = 1, N_atoms
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              CYCLE
           END IF
        END IF
        matrix(N_atoms+1,i_atom) = matrix(N_atoms+1,i_atom) + &
             a(i_reac,i_atom)*n_spec(i_reac)
     END DO
  END DO
  
  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     matrix(N_atoms+1,N_atoms+1) = matrix(N_atoms+1,N_atoms+1) + n_spec(i_reac) !!
  END DO
  matrix(N_atoms+1,N_atoms+1) = matrix(N_atoms+1,N_atoms+1) - n

  ! Set up the matrix for the (N_reactants-N_gas) equations (Eq. 2.25)
  
  IF (N_gas < N_reactants) THEN
     DO i_reac = N_gas+1, N_reactants
        DO i_atom = 1, N_atoms
           matrix(N_atoms+1+i_reac-N_gas,i_atom) = a(i_reac,i_atom)
        END DO
     END DO
  END IF

  ! VECTOR SETUP
  !vector(N_atoms+1+(N_reactants-N_gas))
  vector = 0d0
  
  ! (Eq. 2.25)
  IF (N_gas < N_reactants) THEN
     vector(N_atoms+2:N_atoms+1+(N_reactants-N_gas)) = mu(N_gas+1:N_reactants)/R/temp
  END IF

  ! (Eq. 2.24)
  b = 0d0
  DO i_atom = 1, N_atoms
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              CYCLE
           END IF
        END IF
        b(i_atom) = b(i_atom) + a(i_reac,i_atom)*n_spec(i_reac)
     END DO
     DO i_reac = N_gas+1, N_reactants
        b(i_atom) = b(i_atom) + a(i_reac,i_atom)*n_spec(solid_indices(i_reac-N_gas))
     END DO
  END DO
  vector(1:N_atoms) = b_0 - b
  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     vector(1:N_atoms) = vector(1:N_atoms) + &
          a(i_reac,1:N_atoms)*n_spec(i_reac)*mu(i_reac)/R/temp 
  END DO

  ! (Eq. 2.26)
  vector(N_atoms+1) = n - SUM(n_spec(1:N_gas)) + SUM(n_spec(1:N_gas)*mu(1:N_gas))/R/temp
  
end subroutine ec_MAKE_MATRIX_short

!#####################################################

subroutine ec_INV_MATRIX_long(lens,matrix,vector,solution_vector)

  use thermo_data_block       
  implicit none
  !! I/O:
  INTEGER                      :: lens, write, comp
  DOUBLE PRECISION             :: matrix(lens,lens), matrix_nions(lens-N_ions,lens-N_ions) 
                                 ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
                                 ! condensed species, the pis and the delta log(n)
  DOUBLE PRECISION             :: vector(lens),solution_vector(lens), vector_nions(lens-N_ions), &
       solution_vector_nions(lens-N_ions)
  !! Internal:
  INTEGER                      :: index(lens), corrf_i, corrf_j, index_nions(lens-N_ions)
  DOUBLE PRECISION             :: d
  INTEGER                      :: i_mat, j_mat

  INTERFACE
     SUBROUTINE ludcmp(a,indx,d)
       USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: a
       INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
       DOUBLE PRECISION, INTENT(OUT) :: d
       DOUBLE PRECISION, DIMENSION(size(a,1)) :: vv
       DOUBLE PRECISION, PARAMETER :: TINY=1.0e-20_sp
       INTEGER(I4B) :: j,n,imax
     END SUBROUTINE ludcmp
     SUBROUTINE lubksb(a,indx,b)
       USE nrtype; USE nrutil, ONLY : assert_eq
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: a
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
       DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: b
       INTEGER(I4B) :: i,n,ii,ll
       DOUBLE PRECISION :: summ
     END SUBROUTINE lubksb
  END INTERFACE

  solution_vector = vector

  !matrix(N_reactants+N_atoms+1,N_reactants+N_atoms+1)
  IF (remove_ions) THEN
     vector_nions = 0d0
     matrix_nions = 0d0
     corrf_i = 0
     DO i_mat = 1, lens
        corrf_j = 0
        IF (i_mat <= N_gas) THEN
           IF (reac_ion(i_mat)) THEN
              corrf_i = corrf_i + 1
              cycle
           END IF
        END IF
        DO j_mat = 1, lens
           IF (j_mat <= N_gas) THEN
              IF (reac_ion(j_mat)) THEN
                 corrf_j = corrf_j + 1
                 cycle
              END IF
           END IF
           matrix_nions(j_mat-corrf_j,i_mat-corrf_i) = matrix(j_mat,i_mat)
        END DO
        vector_nions(i_mat-corrf_i) = vector(i_mat)
     END DO
     solution_vector_nions = vector_nions
     
     call ludcmp(matrix_nions,index_nions,d)
     call lubksb(matrix_nions,index_nions,solution_vector_nions)
     corrf_i = 0
     DO i_mat = 1, lens
        IF (i_mat <= N_gas) THEN
           IF (reac_ion(i_mat)) THEN
              corrf_i = corrf_i + 1
              cycle
           END IF
        END IF
        solution_vector(i_mat) = solution_vector_nions(i_mat-corrf_i)
     END DO
  ELSE
     call ludcmp(matrix,index,d)
     call lubksb(matrix,index,solution_vector)
  END IF

end subroutine ec_INV_MATRIX_long

!#####################################################

subroutine ec_INV_MATRIX_short(lens,matrix,vector,solution_vector)

  use thermo_data_block       
  implicit none
  !! I/O:
  INTEGER                      :: lens
  DOUBLE PRECISION             :: matrix(lens,lens) 
                                 ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
                                 ! condensed species, the pis and the delta log(n)
  DOUBLE PRECISION             :: vector(lens),solution_vector(lens)
  !! Internal:
  INTEGER                      :: index(lens)
  DOUBLE PRECISION             :: d

  INTERFACE
     SUBROUTINE ludcmp(a,indx,d)
       USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: a
       INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
       DOUBLE PRECISION, INTENT(OUT) :: d
       DOUBLE PRECISION, DIMENSION(size(a,1)) :: vv
       DOUBLE PRECISION, PARAMETER :: TINY=1.0e-20_sp
       INTEGER(I4B) :: j,n,imax
     END SUBROUTINE ludcmp
     SUBROUTINE lubksb(a,indx,b)
       USE nrtype; USE nrutil, ONLY : assert_eq
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: a
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
       DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: b
       INTEGER(I4B) :: i,n,ii,ll
       DOUBLE PRECISION :: summ
     END SUBROUTINE lubksb
  END INTERFACE

  solution_vector = vector

  call ludcmp(matrix,index,d)
  call lubksb(matrix,index,solution_vector)

end subroutine ec_INV_MATRIX_short

!#####################################################

subroutine ec_CHANGE_ABUNDS_long(N_atoms,N_reactants,solution_vector,n_spec,pi_atom,n,converged,&
     solid_indices,N_solids,names_atoms,molfracs_atoms,N_reactants2,n_spec_old)

  use thermo_data_block       
  implicit none
  !! I/O:
  INTEGER                      :: N_atoms, N_reactants,N_solids,N_reactants2
  INTEGER                      :: solid_indices(N_solids)
  DOUBLE PRECISION             :: solution_vector(N_reactants+N_atoms+1)
  DOUBLE PRECISION             :: n ! Moles of gas particles per total mass of mixture in kg
  DOUBLE PRECISION             :: n_spec(N_reactants2) ! Moles of species per total mass of mixture in kg
  DOUBLE PRECISION             :: n_spec_old(N_reactants2) ! Moles of species per total mass of mixture in kg
  DOUBLE PRECISION             :: pi_atom(N_atoms) ! Lagrangian multipliers for the atomic species divided
                                                   ! by (R*T)
  LOGICAL                      :: converged
  
  !! Internal:
  INTEGER                      :: i_reac
  INTEGER, save                :: n_done = 0
  DOUBLE PRECISION             :: lambda, lambda1, lambda2
  DOUBLE PRECISION, parameter  :: SIZE = 18.420681
  LOGICAL                      :: gas_good, solids_good, total_good

  ! IONS
  INTEGER                      :: i_ion, i_stoich
  DOUBLE PRECISION             :: pi_ion, pi_ion_norm

  ! MASS BALANCE CHECKS
  DOUBLE PRECISION             :: b_0(N_atoms), b_0_norm, mass_atom, b(N_atoms), pi_atom_old(N_atoms)
  DOUBLE PRECISION             :: a(N_reactants,N_atoms), mu(N_reactants), mval_mass_good
  INTEGER                      :: i_atom, i_ratom
  CHARACTER*40                 :: upper_atom_name
  CHARACTER*2                  :: upper_ratom_name
  LOGICAL                      :: mass_good, pi_good
  CHARACTER*40                 :: names_atoms(N_atoms)
  DOUBLE PRECISION             :: molfracs_atoms(N_atoms)


  ! Calculate correction factors as described in Section 3.3 of the McBride Manual
  lambda1 = 9d99
  lambda2 = 9d99
  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     IF (LOG(n_spec(i_reac)/n) > -SIZE) THEN
        lambda1 = MIN(lambda1,2d0/(MAX(5d0*ABS(solution_vector(N_reactants+N_atoms+1)), &
             ABS(solution_vector(i_reac)))))
     ELSE IF ((LOG(n_spec(i_reac)/n) <= -SIZE) .AND. &
          (solution_vector(i_reac) >= 0d0)) THEN
        lambda2 = MIN(lambda2,ABS((-LOG(n_spec(i_reac)/n)-9.2103404)/ &
             (solution_vector(i_reac)-solution_vector(N_reactants+N_atoms+1))))
     END IF
  END DO
  lambda = MIN(1d0,lambda1,lambda2)

  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     n_spec(i_reac) = n_spec(i_reac)*exp(lambda*solution_vector(i_reac))
  END DO
  
  IF (N_gas < N_reactants) THEN
     DO i_reac = N_gas+1, N_reactants
        n_spec(solid_indices(i_reac-N_gas)) = n_spec(solid_indices(i_reac-N_gas))+ &
             lambda*solution_vector(i_reac)
     END DO
  END IF
  pi_atom_old = pi_atom
  pi_atom = solution_vector(N_reactants+1:N_reactants+N_atoms)
  n = n*exp(lambda*solution_vector(N_reactants+N_atoms+1))

  gas_good = .TRUE.
  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     IF (n_spec(i_reac)*ABS(solution_vector(i_reac))/SUM(n_spec) > 0.5d-5) THEN
        gas_good = .FALSE.
     END IF
  END DO
  solids_good = .TRUE.
  IF (N_gas < N_reactants) THEN
     DO i_reac = N_gas+1, N_reactants
        IF (ABS(solution_vector(i_reac))/SUM(n_spec) > 0.5d-5) THEN
           solids_good = .FALSE.
        END IF
     END DO
  END IF
  total_good = .TRUE.
  IF (n*ABS(solution_vector(N_reactants+N_atoms+1))/SUM(n_spec) > 05.d-5) THEN
     total_good = .FALSE.
  END IF

  !!!-----------------------
  
  mass_good = .TRUE.
  pi_good = .TRUE.
  
  
  ! Set up b0
  b_0_norm = 0d0
  DO i_atom = 1, N_atoms
     call ec_ATOM_MASS(names_atoms(i_atom),mass_atom)
     b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
  END DO
  b_0 = molfracs_atoms/b_0_norm

  ! Set up a_ij
  a = 0d0
  DO i_atom = 1, N_atoms
     call To_upper(names_atoms(i_atom),upper_atom_name)
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              a(i_reac,1:N_atoms) = 0d0
              CYCLE
           END IF
        END IF
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
           END IF
        END DO
     END DO
     DO i_reac = N_gas+1, N_reactants
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,solid_indices(i_reac-N_gas)),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
           END IF
        END DO
     END DO
  END DO

  mval_mass_good = MAXVAL(b_0)*1d-2
  DO i_atom = 1, N_atoms
     IF ((abs(b_0(i_atom)-sum(a(1:N_reactants,i_atom)*n_spec(1:N_reactants))) > mval_mass_good) .AND. (b_0(i_atom) > 1d-6)) THEN
        mass_good = .FALSE.
     END IF
  END DO

  DO i_atom = 1, N_atoms
     IF (abs((pi_atom_old(i_atom)-pi_atom(i_atom))/pi_atom(i_atom)) > 1d-3) THEN
        pi_good = .FALSE.
     END IF
  END DO

  IF ((.NOT. mass_good) .OR. (.NOT. pi_good)) THEN
     mass_good = .TRUE.
     pi_good = .TRUE.   
     DO i_reac = 1, N_reactants2
        IF (ABS(n_spec(i_reac)-n_spec_old(i_reac)) > 1d-10) THEN
           mass_good = .FALSE.
           pi_good = .FALSE.
        END IF
     END DO
  END IF

  !!!-------------------

  ! ION CONVERGENCE?

  IF (ions .AND. (.NOT. remove_ions)) THEN

     ! DO THE MAGIC THEY DO IN SECT. 3.7 in McBride
     
     pi_ion = 0d0
     pi_ion_norm = 0d0
     DO i_reac = 1, N_reactants
        DO i_stoich = 1, 5
           IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
              pi_ion = pi_ion - n_spec(i_reac)*reac_stoich(i_stoich,i_reac)
              pi_ion_norm = pi_ion_norm + n_spec(i_reac)*reac_stoich(i_stoich,i_reac)**2d0
              EXIT
           END IF
        END DO
     END DO

     pi_ion = pi_ion / pi_ion_norm

     IF (ABS(pi_ion) > 1d-4) THEN
        DO i_ion = 1, 80
           DO i_reac = 1, N_reactants
              DO i_stoich = 1, 5
                 IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                    n_spec(i_reac) = n_spec(i_reac)*exp(reac_stoich(i_stoich,i_reac)*pi_ion)
                    EXIT
                 END IF
              END DO
           END DO

           pi_ion = 0d0
           pi_ion_norm = 0d0
           DO i_reac = 1, N_reactants
              DO i_stoich = 1, 5
                 IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                    pi_ion = pi_ion - n_spec(i_reac)*reac_stoich(i_stoich,i_reac)
                    pi_ion_norm = pi_ion_norm + n_spec(i_reac)*reac_stoich(i_stoich,i_reac)**2d0
                    EXIT
                 END IF
              END DO
           END DO

        END DO
     END IF

     IF (((gas_good .AND. solids_good) .AND. total_good) .AND. (ABS(pi_ion) <= 1d-4)) THEN
        IF (mass_good .AND. pi_good) THEN
           converged = .TRUE.
        END IF
     END IF

  ELSE

     IF ((gas_good .AND. solids_good) .AND. total_good) THEN
        IF (mass_good .AND. pi_good) THEN
           converged = .TRUE.
        END IF
     END IF

  END IF
     
  n_done = n_done + 1

end subroutine ec_CHANGE_ABUNDS_long

!#####################################################

subroutine ec_CHANGE_ABUNDS_short(N_atoms,N_reactants,solution_vector,n_spec,pi_atom,n,converged,&
     solid_indices,N_solids,mu_gas,a_gas,temp,names_atoms,molfracs_atoms,N_reactants2,n_spec_old)

  use thermo_data_block       
  implicit none
  !! I/O:
  INTEGER                      :: N_atoms, N_reactants,N_solids, N_reactants2
  INTEGER                      :: solid_indices(N_solids)
  DOUBLE PRECISION             :: solution_vector(N_atoms+1+(N_reactants-N_gas))
  DOUBLE PRECISION             :: n ! Moles of gas particles per total mass of mixture in kg
  DOUBLE PRECISION             :: n_spec(N_reactants2) ! Moles of species per total mass of mixture in kg
  DOUBLE PRECISION             :: n_spec_old(N_reactants2) ! Moles of species per total mass of mixture in kg
  DOUBLE PRECISION             :: pi_atom(N_atoms) ! Lagrangian multipliers for the atomic species divided
                                                   ! by (R*T)
  LOGICAL                      :: converged
  DOUBLE PRECISION             :: mu_gas(N_gas), a_gas(N_gas,N_atoms), temp
  
  !! Internal:
  INTEGER                      :: i_reac
  INTEGER, save                :: n_done = 0
  DOUBLE PRECISION             :: lambda, lambda1, lambda2
  DOUBLE PRECISION, parameter  :: SIZE = 18.420681
  LOGICAL                      :: gas_good, solids_good, total_good
  DOUBLE PRECISION             :: delta_n_gas(N_gas)

  ! IONS
  INTEGER                      :: i_ion, i_stoich
  DOUBLE PRECISION             :: pi_ion, pi_ion_norm

  ! MASS BALANCE CHECKS
  DOUBLE PRECISION             :: b_0(N_atoms), b_0_norm, mass_atom, b(N_atoms), pi_atom_old(N_atoms)
  DOUBLE PRECISION             :: a(N_reactants,N_atoms), mu(N_reactants), mval_mass_good
  INTEGER                      :: i_atom, i_ratom
  CHARACTER*40                 :: upper_atom_name
  CHARACTER*2                  :: upper_ratom_name
  LOGICAL                      :: mass_good, pi_good
  CHARACTER*40                 :: names_atoms(N_atoms)
  DOUBLE PRECISION             :: molfracs_atoms(N_atoms)

  ! Get delta_n_gas, following Eq. 2.18:
  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     delta_n_gas(i_reac) = SUM(a_gas(i_reac,1:N_atoms)*solution_vector(1:N_atoms)) + &
          solution_vector(N_atoms+1) - mu_gas(i_reac)/R/temp
  END DO

  ! Calculate correction factors as described in Section 3.3 of the McBride Manual
  lambda1 = 9d99
  lambda2 = 9d99
  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     IF (LOG(n_spec(i_reac)/n) > -SIZE) THEN
        lambda1 = MIN(lambda1,2d0/(MAX(5d0*ABS(solution_vector(N_atoms+1)), &
             ABS(delta_n_gas(i_reac)))))
     ELSE IF ((LOG(n_spec(i_reac)/n) <= -SIZE) .AND. &
          (delta_n_gas(i_reac) >= 0d0)) THEN
        lambda2 = MIN(lambda2,ABS((-LOG(n_spec(i_reac)/n)-9.2103404)/ &
             (delta_n_gas(i_reac)-solution_vector(N_atoms+1))))
     END IF
  END DO
  lambda = MIN(1d0,lambda1,lambda2)

  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     n_spec(i_reac) = n_spec(i_reac)*exp(lambda*delta_n_gas(i_reac))
  END DO
     
  IF (N_gas < N_reactants) THEN
     DO i_reac = N_gas+1, N_reactants
        n_spec(solid_indices(i_reac-N_gas)) = n_spec(solid_indices(i_reac-N_gas))+ &
             lambda*solution_vector(N_atoms+1+i_reac-N_gas)
     END DO
  END IF
  pi_atom_old = pi_atom
  pi_atom = solution_vector(1:N_atoms)
  n = n*exp(lambda*solution_vector(N_atoms+1))

  gas_good = .TRUE.
  DO i_reac = 1, N_gas
     IF (remove_ions) THEN
        IF (reac_ion(i_reac)) THEN
           CYCLE
        END IF
     END IF
     IF (n_spec(i_reac)*ABS(delta_n_gas(i_reac))/SUM(n_spec) > 0.5d-5) THEN
        gas_good = .FALSE.
     END IF
  END DO
  solids_good = .TRUE.
  IF (N_gas < N_reactants) THEN
     DO i_reac = N_gas+1, N_reactants
        IF (ABS(solution_vector(N_atoms+1+i_reac-N_gas))/SUM(n_spec) > 0.5d-5) THEN
           solids_good = .FALSE.
        END IF
     END DO
  END IF
  total_good = .TRUE.
  IF (n*ABS(solution_vector(N_atoms+1))/SUM(n_spec) > 05.d-5) THEN
     total_good = .FALSE.
  END IF

  !!!-----------------------
  
  mass_good = .TRUE.
  pi_good = .TRUE.
  
  ! Set up b0
  b_0_norm = 0d0
  DO i_atom = 1, N_atoms
     call ec_ATOM_MASS(names_atoms(i_atom),mass_atom)
     b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
  END DO
  b_0 = molfracs_atoms/b_0_norm

  ! Set up a_ij
  a = 0d0
  DO i_atom = 1, N_atoms
     call To_upper(names_atoms(i_atom),upper_atom_name)
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              a(i_reac,1:N_atoms) = 0d0
              CYCLE
           END IF
        END IF
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
           END IF
        END DO
     END DO
     DO i_reac = N_gas+1, N_reactants
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,solid_indices(i_reac-N_gas)),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
           END IF
        END DO
     END DO
  END DO

  mval_mass_good = MAXVAL(b_0)*1d-2
  DO i_atom = 1, N_atoms
     IF ((abs(b_0(i_atom)-sum(a(1:N_reactants,i_atom)*n_spec(1:N_reactants))) > mval_mass_good) .AND. (b_0(i_atom) > 1d-6)) THEN
        mass_good = .FALSE.
     END IF
  END DO

  DO i_atom = 1, N_atoms
     IF (abs((pi_atom_old(i_atom)-pi_atom(i_atom))/pi_atom(i_atom)) > 1d-3) THEN
        pi_good = .FALSE.
     END IF
  END DO

  IF ((.NOT. mass_good) .OR. (.NOT. pi_good)) THEN
     mass_good = .TRUE.
     pi_good = .TRUE.   
     DO i_reac = 1, N_reactants2
        IF (ABS(n_spec(i_reac)-n_spec_old(i_reac)) > 1d-10) THEN
           mass_good = .FALSE.
           pi_good = .FALSE.
        END IF
     END DO
  END IF
  
  !!!-------------------

  ! ION CONVERGENCE?

  IF (ions .AND. (.NOT. remove_ions)) THEN

     ! DO THE MAGIC THEY DO IN SECT. 3.7 in McBride
  
     pi_ion = 0d0
     pi_ion_norm = 0d0
     DO i_reac = 1, N_reactants
        DO i_stoich = 1, 5
           IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
              pi_ion = pi_ion - n_spec(i_reac)*reac_stoich(i_stoich,i_reac)
              pi_ion_norm = pi_ion_norm + n_spec(i_reac)*reac_stoich(i_stoich,i_reac)**2d0
              EXIT
           END IF
        END DO
     END DO

     pi_ion = pi_ion / pi_ion_norm

     IF (ABS(pi_ion) > 1d-4) THEN
        DO i_ion = 1, 80
           DO i_reac = 1, N_reactants
              DO i_stoich = 1, 5
                 IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                    n_spec(i_reac) = n_spec(i_reac)*exp(reac_stoich(i_stoich,i_reac)*pi_ion)
                    EXIT
                 END IF
              END DO
           END DO

           pi_ion = 0d0
           pi_ion_norm = 0d0
           DO i_reac = 1, N_reactants
              DO i_stoich = 1, 5
                 IF (trim(adjustl(reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                    pi_ion = pi_ion - n_spec(i_reac)*reac_stoich(i_stoich,i_reac)
                    pi_ion_norm = pi_ion_norm + n_spec(i_reac)*reac_stoich(i_stoich,i_reac)**2d0
                    EXIT
                 END IF
              END DO
           END DO

        END DO
     END IF

     IF (((gas_good .AND. solids_good) .AND. total_good) .AND. (ABS(pi_ion) <= 1d-4)) THEN
        IF (mass_good .AND. pi_good) THEN
           converged = .TRUE.
        END IF
     END IF

  ELSE

     IF ((gas_good .AND. solids_good) .AND. total_good) THEN
        IF (mass_good .AND. pi_good) THEN
           converged = .TRUE.
        END IF
     END IF

  END IF
     
  n_done = n_done + 1

end subroutine ec_CHANGE_ABUNDS_short

!#####################################################

subroutine ec_ATOM_MASS(atom_name,atom_mass)

  use thermo_data_block
  implicit none

  !! I/O:
  CHARACTER(*)                 :: atom_name
  DOUBLE PRECISION             :: atom_mass
  !! Internal:
  INTEGER                      :: i_atom
  CHARACTER(len(atom_name))    :: atom_name_UPPER
  CHARACTER*2                  :: names_atoms_save_UPPER

  atom_mass = -1d0

  DO i_atom = 1, N_atoms_save
     call To_upper(atom_name,atom_name_UPPER)
     call To_upper(names_atoms_save(i_atom),names_atoms_save_UPPER)
     IF (trim(adjustl(atom_name_UPPER)) .EQ. &
          trim(adjustl(names_atoms_save_UPPER))) THEN
        atom_mass = masses_atoms_save(i_atom)
        EXIT
     END IF
  END DO

  IF (atom_mass < 0d0) THEN
     WRITE(*,*) 'EASY CHEM ERROR: Mass of atom of species '// &
          trim(adjustl(atom_name))//' not found.'
     STOP
  END IF

end subroutine ec_ATOM_MASS

!#####################################################

subroutine ec_INCLUDE_SOLIDS_QUESTIONMARK(N_reactants,N_atoms,pi_atom,names_atoms,H_0,S_0,dgdnj,temp, &
     solid_inc,inc_next,neg_cond)

  use thermo_data_block
  implicit none
  !! I/O
  INTEGER                      :: N_reactants, N_atoms, inc_next
  DOUBLE PRECISION             :: pi_atom(N_atoms), temp
  CHARACTER*40                 :: names_atoms(N_atoms)
  DOUBLE PRECISION             :: H_0(N_reactants), S_0(N_reactants), dgdnj(N_reactants-N_gas)
  LOGICAL                      :: solid_inc(N_reactants-N_gas), all_needed, neg_cond(N_reactants-N_gas)

  !! Internal:
  DOUBLE PRECISION             :: a(N_reactants,N_atoms)
  DOUBLE PRECISION             :: mu(N_reactants), minval_inc
  INTEGER                      :: i_atom, i_reac, i_ratom
  CHARACTER*40                 :: upper_atom_name
  CHARACTER*2                  :: upper_ratom_name

  inc_next = -1
  minval_inc = 1d0
  
  ! Set up a_ij
  a = 0d0
  DO i_atom = 1, N_atoms
     call To_upper(names_atoms(i_atom),upper_atom_name)
     DO i_reac = 1, N_reactants
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              CYCLE
           END IF
        END IF
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
           END IF
        END DO
     END DO
  END DO

  ! EVAL Eq. 3.7 in McBride Manual
  DO i_reac = N_gas+1, N_reactants
     mu(i_reac) = H_0(i_reac) - temp*S_0(i_reac)
     dgdnj(i_reac-N_gas) = mu(i_reac)/R/temp - SUM(a(i_reac,1:N_atoms)*pi_atom)
  END DO

  DO i_reac = N_gas+1, N_reactants
     IF ((dgdnj(i_reac-N_gas) < 0d0) .AND. (.NOT. solid_inc(i_reac-N_gas))) THEN
        IF (((dgdnj(i_reac-N_gas) < minval_inc) .AND. (.NOT. neg_cond(i_reac-N_gas)) .AND. &
             temp <= thermo_data_temps(2,thermo_data_n_intervs(i_reac),i_reac)) .AND. &
             (temp >= thermo_data_temps(1,1,i_reac))) THEN
           minval_inc = dgdnj(i_reac-N_gas)
           inc_next = i_reac
        END IF
     END IF
  END DO

end subroutine ec_INCLUDE_SOLIDS_QUESTIONMARK

!#####################################################

subroutine ec_CALC_ADIABATIC_GRADIENT(N_atoms,N_spec_eff,N_reactants,n_spec, &
     n,H_0,C_P_0,solid_indices,N_solids,temp,names_atoms,nabla_ad,gamma2,c_pe)

  use thermo_data_block
  implicit none

    !! I/O:
  INTEGER                      :: N_atoms, N_spec_eff, N_reactants, N_solids
  INTEGER                      :: solid_indices(N_solids)
  CHARACTER*40                 :: names_atoms(N_atoms)
  DOUBLE PRECISION             :: temp
  DOUBLE PRECISION             :: C_P_0(N_reactants), H_0(N_reactants)
  DOUBLE PRECISION             :: n ! Moles of gas particles per total mass of mixture in kg
  DOUBLE PRECISION             :: n_spec(N_reactants) ! Moles of species per total mass of mixture in kg
  DOUBLE PRECISION             :: matrix(N_atoms+1+(N_spec_eff-N_gas),N_atoms+1+(N_spec_eff-N_gas))
                                 ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
                                 ! condensed species, the pis and the delta log(n)
  DOUBLE PRECISION             :: vector(N_atoms+1+(N_spec_eff-N_gas)), &
       solution_vector(N_atoms+1+(N_spec_eff-N_gas))
  DOUBLE PRECISION             :: nabla_ad, gamma2
  
  !! Internal:
  DOUBLE PRECISION             :: a(N_spec_eff,N_atoms), c_pe
  INTEGER                      :: i_atom, i_reac, i_ratom, i_atom2
  CHARACTER*40                 :: upper_atom_name
  CHARACTER*2                  :: upper_ratom_name

  ! Set up a_ij
  a = 0d0
  DO i_atom = 1, N_atoms
     call To_upper(names_atoms(i_atom),upper_atom_name)
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              a(i_reac,i_atom) = 0d0
              CYCLE
           END IF
        END IF
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,i_reac)*mol
           END IF
        END DO
     END DO
     DO i_reac = N_gas+1, N_spec_eff
        DO i_ratom = 1, 5
           call To_upper(reac_atoms_names(i_ratom,solid_indices(i_reac-N_gas)),upper_ratom_name)
           IF (trim(adjustl(upper_atom_name)) .EQ. &
                trim(adjustl(upper_ratom_name))) THEN
              a(i_reac,i_atom) = reac_stoich(i_ratom,solid_indices(i_reac-N_gas))*mol
           END IF
        END DO
     END DO
  END DO

  matrix = 0d0
  ! Setup matrix, following Eq. 2.56
  DO i_atom = 1, N_atoms
     ! First term, LHS
     DO i_atom2 = 1, N_atoms
        DO i_reac = 1, N_gas
           IF (remove_ions) THEN
              IF (reac_ion(i_reac)) THEN
                 CYCLE
              END IF
           END IF
           matrix(i_atom,i_atom2) = matrix(i_atom,i_atom2) + &
                n_spec(i_reac)*a(i_reac,i_atom2)*a(i_reac,i_atom)
        END DO
     END DO
     ! Second term, LHS
     DO i_reac = N_gas+1, N_spec_eff
        matrix(i_atom,N_atoms+1+i_reac-N_gas) = a(i_reac,i_atom)
     END DO
     ! Third term, LHS
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              CYCLE
           END IF
        END IF
        matrix(i_atom,N_atoms+1) = matrix(i_atom,N_atoms+1) + &
             a(i_reac,i_atom)*n_spec(i_reac)
     END DO
  END DO

  ! Setup matrix, following Eq. 2.58
  DO i_atom = 1, N_atoms
     DO i_reac = 1, N_gas
        IF (remove_ions) THEN
           IF (reac_ion(i_reac)) THEN
              CYCLE
           END IF
        END IF
        matrix(N_atoms+1,i_atom) = matrix(N_atoms+1,i_atom) + &
             a(i_reac,i_atom)*n_spec(i_reac)
     END DO
  END DO

  ! Setup matrix, following Eq. 2.57
  DO i_reac = N_gas+1, N_spec_eff
     DO i_atom = 1, N_atoms
        matrix(N_atoms+1+i_reac-N_gas,i_atom) = a(i_reac,i_atom)
     END DO
  END DO

  vector = 0d0
  ! Setup the vector, following Eq. 2.56
  DO i_atom = 1, N_atoms
     vector(i_atom) = -SUM(a(1:N_gas,i_atom)*n_spec(1:N_gas)*H_0(1:N_gas)) &
          /R/temp
  END DO

  ! Setup the vector, following Eq. 2.58
  vector(N_atoms+1) = -SUM(n_spec(1:N_gas)*H_0(1:N_gas))/R/temp

  ! Setup the vector, following Eq. 2.57
  DO i_reac = N_gas+1, N_spec_eff
     vector(N_atoms+1+i_reac-N_gas) = -H_0(solid_indices(i_reac-N_gas))/R/temp
  END DO

  ! Solve the system
  call ec_INV_MATRIX_short(N_atoms+1+N_spec_eff-N_gas,matrix,vector,solution_vector)

  ! Calculate c_pe, following Eq. 2.59
  c_pe = 0d0
  DO i_atom = 1, N_atoms
     c_pe = c_pe + SUM(a(1:N_gas,i_atom)*n_spec(1:N_gas)*H_0(1:N_gas)/R/temp) * &
          solution_vector(i_atom)
  END DO
  DO i_reac = N_gas+1, N_spec_eff
     c_pe = c_pe + H_0(solid_indices(i_reac-N_gas))/R/temp* &
          solution_vector(N_atoms+1+i_reac-N_gas) + &
          n_spec(solid_indices(i_reac-N_gas))* &
          C_P_0(solid_indices(i_reac-N_gas))/R
  END DO
  c_pe = c_pe + SUM(n_spec(1:N_gas)*C_P_0(1:N_gas)/R)
  c_pe = c_pe + SUM(n_spec(1:N_gas)*H_0(1:N_gas)/R/temp)* &
       solution_vector(N_atoms+1)
  c_pe = c_pe + SUM(n_spec(1:N_gas)*(H_0(1:N_gas)/R/temp)**2d0)
  c_pe = c_pe*R

  ! Calculate nabla_ad, using Eq. 2.50 and Eq. 2.75
  nabla_ad = c_pe/n/R/(1d0+solution_vector(N_atoms+1))
  nabla_ad = 1/nabla_ad
  gamma2 = 1d0/(1d0-nabla_ad)
  
end subroutine ec_CALC_ADIABATIC_GRADIENT

!#####################################################

!!!!!!!!!!!!!
! START FOREIGN CODE
!!!!!!!!!!!!!

!TAKEN FROM http://rosettacode.org/wiki/String_case#Fortran :

subroutine To_upper(strin,strout)
  character(*) :: strin,strout
  integer :: i

  strout = strin
  do i = 1, len(strout)
     select case(strout(i:i))
     case("a":"z")
        strout(i:i) = achar(iachar(strout(i:i))-32)
     end select
  end do
end subroutine To_upper

!TAKEN FROM NUMERICAL RECIPES :

SUBROUTINE ludcmp(a,indx,d)
  USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
  INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
  REAL(DP), INTENT(OUT) :: d
  REAL(DP), DIMENSION(size(a,1)) :: vv
  REAL(DP), PARAMETER :: TINY=1.0e-20_sp
  INTEGER(I4B) :: j,n,imax
  n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
  d=1.0
  vv=maxval(abs(a),dim=2)
  if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
  vv=1.0_sp/vv
  do j=1,n
     imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
     if (j /= imax) then
        call swap(a(imax,:),a(j,:))
        d=-d
        vv(imax)=vv(j)
     end if
     indx(j)=imax
     if (a(j,j) == 0.0) a(j,j)=TINY
     a(j+1:n,j)=a(j+1:n,j)/a(j,j)
     a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
  end do
END SUBROUTINE ludcmp

!TAKEN FROM NUMERICAL RECIPES :

SUBROUTINE lubksb(a,indx,b)
  USE nrtype; USE nrutil, ONLY : assert_eq
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
  INTEGER(I4B) :: i,n,ii,ll
  REAL(DP) :: summ
  n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
  ii=0
  do i=1,n
     ll=indx(i)
     summ=b(ll)
     b(ll)=b(i)
     if (ii /= 0) then
        summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
     else if (summ /= 0.0) then
        ii=i
     end if
     b(i)=summ
  end do
  do i=n,1,-1
     b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
  end do
END SUBROUTINE lubksb

!!!!!!!!!!!!!
! END FOREIGN CODE
!!!!!!!!!!!!!


!**********************************************************
! TAKEN FROM
! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
!**********************************************************

subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to OR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = 0!getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed
