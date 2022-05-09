! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019-2020
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains all cross sections for electrons inelastic scattering 
! References used:
! [1] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2008: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2009)



module CS_ion_inelastic
use Universal_constants
use Relativity, only: rest_energy
use Objects
use Initial_conditions, only : count_types_of_SHIs, repeated_type
use CDF_delta, only : integrated_delta_CDF_CS, P_prefactor, Find_linear_a_b, integral_CS_x_W, energy_loss_delta, CDF_total_CS_delta
use Read_input_data, only : m_folder_materials, m_input_folder, m_ion_CS, m_ion_inelast_CS
use Dealing_with_files, only : read_file
use SHI_charge_state, only : Equilibrium_charge_SHI
use CS_general_tools, only: MFP_from_sigma
use CS_integration_limits, only: find_Wmax_equal_Wmin, W_min, W_max, W_max_nonrel_free
use CS_electrons_inelastic, only: CDF_total_CS_nonrel
use Little_subroutines, only: interpolate_data_single

implicit none


real(8) :: m_Eeq_factor
! Used for testing purposes:
parameter (m_Eeq_factor = 50e6)    ! Ekin of ion [eV]

 contains
 
!-----------------------------------------------------------------------------
! Subroutines for inelastic scattering:

subroutine get_ion_IMFP(Material, numpar, bunch, MC, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material  !material parameters of each target that it's constructed of
   type(Num_par), intent(in), target :: numpar  ! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   type(MC_arrays), dimension(:), intent(in) :: MC	! all MC arrays for particles: photons, electrons and holes; size equals to number of iterations
   type(Error_handling), intent(inout) :: Err   ! error log
   !-------------------------------- 
   integer :: N_targets, i, j, k, m, ibunch, N_bunch, ipart, N_types_SHI, Ngrid, N_elements, N_shells, FN, Reason, count_lines, N_SHI, Nsiz
   real(8) :: sigma, Se, lambda,  sigma_cur, Se_cur, lambda_cur
   real(8) :: Total_Se, Total_lambda, N_elem, elem_contrib, ph_to_ion
   character(200) :: Path, Folder_with_CS, command, Path_valent, File_name
   character(20) :: Model_name, Ion_name, temp_c2, temp_c3
   logical :: file_exist, file_opened, read_well
   !-------------------------------- 
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   !-------------------------------- 

   ph_to_ion = 1.0d3    ! rescaling factor to set the grid [eV -> keV]
   
   ! Check if there is any SHI to calcualte:
   N_types_SHI = count_types_of_SHIs(bunch) ! module "Initial_conditions"
   
   ! Only if there are some ions, calculate MFPs for them
   if (N_types_SHI < 1) goto 9999   ! if there are none, skip the whole routine

   ! Find how many sources of radiation are used:
   N_bunch = size (bunch)
   ipart = 0
   ! Check for all of them, whether there are ions (do we need to include ions here or not):
   BNCH:do ibunch = 1, N_bunch
      if (bunch(ibunch)%KOP == 3) then ! it is an ion
         Ion_name = bunch(ibunch)%Name  ! this is the name of the SHI
         ipart = ipart + 1
         write(*, '(a)', advance='no') ' Obtaining ion inelastic scattering cross sections for '// trim(adjustl(Ion_name)) //'...'

         path_sep => numpar%path_sep
         Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_ion_CS))	! Ion CSs to be (or are) storred in this folder
         N_targets = size(Material)	! that's how many different targets user specified

         ! Get ion MFPs:
         TRGT:do i = 1, N_targets    ! for each target
            N_elem = dble(SUM(Material(i)%Elements(:)%percentage))	! number of elements in this compound material
            N_elements = size(Material(i)%Elements)   ! that's how many different elements are in this target
            LMNT:do j =1, N_elements  ! for each element
               Element => Material(i)%Elements(j) ! all information about this element
               if (.not.allocated(Element%SHI_inelastic)) allocate(Element%SHI_inelastic(N_types_SHI)) ! that's how many different types of ions we have in this calculation

               ! Check if a folder for this element already exists:
               Folder_with_CS = trim(adjustl(Path))//path_sep//trim(adjustl(Element%Name))
               inquire(DIRECTORY=trim(adjustl(Folder_with_CS)),exist=file_exist)    ! check if input file excists
               if (.not.file_exist) then  ! to make sure that such a folder is present (even if empty)
                  ! Create a new directory for output files:
                  command='mkdir '//trim(adjustl(Folder_with_CS)) ! to create a folder use this command
                  CALL system(command)    ! create the folder
               endif
               ! Check if a folder for this SHI already exists:
               Folder_with_CS = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(Ion_name))
               inquire(DIRECTORY=trim(adjustl(Folder_with_CS)),exist=file_exist)    ! check if input file excists
               if (.not.file_exist) then  ! to make sure that such a folder is present (even if empty)
                  ! Create a new directory for output files:
                  command='mkdir '//trim(adjustl(Folder_with_CS)) ! to create a folder use this command
                  CALL system(command)    ! create the folder
               endif

               !-------------------------------------
               write(temp_c2,'(f8.1)') numpar%CDF_Eeq_factor
               ! Set part of the file name corresponding to the model used for ion inelastic CS:
               select case (numpar%SHI_inelast)    ! so far, we can only have CDF-based models for SHI
               case (1:3)	! CDF with Delta oscillators
                  select case (numpar%CDF_model)
                  case default   ! Ritchie
                     write(temp_c3,'(a)') '_R_'    ! Ritchie
                  case (1)   ! Mermin
                     write(temp_c3,'(a)') '_M_'    ! Mermin
                  endselect
                  Model_name = 'CDF_Delta'//trim(adjustl(temp_c3))//'Eeq_'//trim(adjustl(temp_c2))
               case (4)	! CDF with Ritchi's oscillators
                  select case (numpar%CDF_model)
                  case default   ! Ritchie
                     write(temp_c3,'(a)') '_R'    ! Ritchie
                  case (1)   ! Mermin
                     write(temp_c3,'(a)') '_M'    ! Mermin
                  endselect
                  Model_name = 'CDF'//trim(adjustl(temp_c3))
               case (5)	! CDF with single pole delta-functions
                  select case (numpar%CDF_model)
                  case default   ! Ritchie
                     write(temp_c3,'(a)') '_R_'    ! Ritchie
                  case (1)   ! Mermin
                     write(temp_c3,'(a)') '_M_'    ! Mermin
                  endselect
                  Model_name = 'CDF_SP_delta'//trim(adjustl(temp_c3))//'Eeq_'//trim(adjustl(temp_c2))
               case default ! exclude
                  Model_name = 'NO'
               end select

               ! Path to files with the valence band and total (material specific) MFPs:
!                if (allocated(Material(i)%CDF_valence%A)) then    ! Valence band
                  Path_valent = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
                  Path_valent = trim(adjustl(Path_valent))//path_sep//trim(adjustl(Material(i)%Name))
!                endif

               ! For electron CSs use the same grid as for photons:
               Ngrid = size(Element%Phot_absorption%E)
               allocate(Element%SHI_inelastic(ipart)%E(Ngrid))
               Element%SHI_inelastic(ipart)%E = Element%Phot_absorption%E*ph_to_ion  ! [eV] -> [keV] energy grid for SHIs
               allocate(Element%SHI_inelastic(ipart)%Total(Ngrid))
               allocate(Element%SHI_inelastic(ipart)%Total_MFP(Ngrid))
               allocate(Element%SHI_inelastic(ipart)%Total_Se(Ngrid))
               allocate(Element%SHI_inelastic(ipart)%Per_shell(Element%N_shl,Ngrid))
               allocate(Element%SHI_inelastic(ipart)%MFP(Element%N_shl,Ngrid))
               allocate(Element%SHI_inelastic(ipart)%Se(Element%N_shl,Ngrid))
               N_shells = Element%N_shl
               
               ! Calculate total cross sections for all shells of this element:
              do k = 1, N_shells
                  ! Check valence band:
                  ! If this shell forms valent band, and data for the valence band exists, separate it from the atomic shells:
                  VAL:if ( (Element%valent(k)) .and. (allocated(Material(i)%CDF_valence%A)) ) then    ! Valence band
                     ! Allocate Valent cross sections for each type of SHI:
                     if (.not.allocated(Material(i)%SHI_inelastic_valent)) allocate(Material(i)%SHI_inelastic_valent(N_types_SHI))
                     ! Allocate parameters of the valence CSs:
                     if (.not.allocated(Material(i)%SHI_inelastic_valent(ipart)%E)) then   ! valence band has not been defined yet, so do that
                        Nsiz = size(Material(i)%Ph_absorption_valent%E)   ! valence band has a different grid from core shells
                        allocate(Material(i)%SHI_inelastic_valent(ipart)%E(Nsiz))
                        Material(i)%SHI_inelastic_valent(ipart)%E = Material(i)%Ph_absorption_valent%E*ph_to_ion  ! [eV] -> [keV] energy grid for SHIs
                        allocate(Material(i)%SHI_inelastic_valent(ipart)%Total(Nsiz))
                        allocate(Material(i)%SHI_inelastic_valent(ipart)%Total_MFP(Nsiz))
                        allocate(Material(i)%SHI_inelastic_valent(ipart)%Total_Se(Nsiz))
                        ! Check file with CSs:
                        File_name = trim(adjustl(Path_valent))//path_sep//trim(adjustl(Ion_name))//'_in_'//trim(adjustl(Material(i)%Name))//trim(adjustl(m_ion_inelast_CS))//'_valence_'//trim(adjustl(Model_name))//'.dat'
                        inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!                         if (file_exist) then	! just read from the file, no need to recalculate:
                        if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
                           open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
                           ! Get the CSs:
                           count_lines = 0
                           do m = 1, Nsiz	! for all energy grid points:
                              E => Material(i)%SHI_inelastic_valent(ipart)%E(m)	! electron energy [eV]
                              read(FN,'(es,es,es,es)', IOSTAT=Reason) E, Material(i)%SHI_inelastic_valent(ipart)%Total(m), Material(i)%SHI_inelastic_valent(ipart)%Total_MFP(m), Material(i)%SHI_inelastic_valent(ipart)%Total_Se(m)
                              call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                              if (.not.read_well) then
                                 close(FN)	! redo the file
                                 goto 8391
                              endif
                              ! Recalculate the MFP since material density might be different from the one saved in the file:
                              Material(i)%SHI_inelastic_valent(ipart)%Total_MFP(m) = MFP_from_sigma(Material(i)%SHI_inelastic_valent(ipart)%Total(m),  Material(i)%At_Dens) ! [A] module "CS_general_tools"
                           enddo
                        else	! no such file => create it
8391                    open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
                           do m = 1, Nsiz	! for all energy grid points:
                              E => Material(i)%SHI_inelastic_valent(ipart)%E(m)	! electron energy [eV]
                              sigma = 0.0d0
                              Se = 0.0d0
                              call get_SHI_inelastic_CS(E, bunch(ibunch)%Z, bunch(ibunch)%Meff, Material(i), Material(i)%DOS, numpar, bunch, ibunch, j, k, Material(i)%DOS%Egap, sigma, Se=Se)	! see below
                              Material(i)%SHI_inelastic_valent(ipart)%Total(m) = sigma   ! [A^2]
                              lambda = MFP_from_sigma(sigma,  Material(i)%At_Dens)    ! module "CS_general_tools"
                              Material(i)%SHI_inelastic_valent(ipart)%Total_MFP(m) = lambda  ! [1/A]
                              Material(i)%SHI_inelastic_valent(ipart)%Total_Se(m) = Se    ! [eV/A]
                              write(FN,'(es,es,es,es)') Material(i)%SHI_inelastic_valent(ipart)%E(m), sigma, lambda, Se
                           enddo ! m = 1, Nsiz
                        endif ! (file_exist)
                        close(FN)
                     endif ! (.not.allocated(Material(i)%SHI_inelastic_valent(ipart)%E)) 
                  else VAL    ! core shell
                     ! Check file with CSs:
                     File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(Ion_name))//'_in_'//trim(adjustl(Element%Name))//trim(adjustl(m_ion_inelast_CS))//'_'//trim(adjustl(Element%Shell_name(k)))//'_'//trim(adjustl(Model_name))//'.dat'
                     !File_name = trim(adjustl(Path_valent))//path_sep//trim(adjustl(Ion_name))//'_in_'//trim(adjustl(Element%Name))//trim(adjustl(m_ion_inelast_CS))//'_valence_'//trim(adjustl(Model_name))//'.dat'
                     inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!                      if (file_exist) then	! just read from the file, no need to recalculate:
                     if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
                        open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
                        ! Get the MFP and CDF points:
                        count_lines = 0
                        do m = 1, Ngrid	! for all energy grid points:
                           E => Element%SHI_inelastic(ipart)%E(m)	! electron energy [eV]
                           read(FN,'(es,es,es,es)', IOSTAT=Reason) E, Element%SHI_inelastic(ipart)%Per_shell(k,m), Element%SHI_inelastic(ipart)%MFP(k,m), Element%SHI_inelastic(ipart)%Se(k,m)
                           call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                           if (.not.read_well) then
                              close(FN)	! redo the file
                              goto 8390
                           endif
                        enddo
                     else	! no such file => create it
8390                 open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
                        do m = 1, Ngrid	! for all energy grid points:
                           E => Element%SHI_inelastic(ipart)%E(m)	! electron energy [eV]
                           sigma = 0.0d0
                           Se = 0.0d0
                           call get_SHI_inelastic_CS(E, bunch(ibunch)%Z, bunch(ibunch)%Meff, Material(i), Material(i)%DOS, numpar, bunch, ibunch, j, k, Material(i)%Elements(j)%Ip(k), sigma, Se=Se)	! see below
                           Element%SHI_inelastic(ipart)%Per_shell(k,m) = sigma	! [A^2]
                           Element%SHI_inelastic(ipart)%Se(k,m) = Se * 1.0d24/Material(i)%At_Dens    ! [eV/A]
                           Element%SHI_inelastic(ipart)%MFP(k,m) = MFP_from_sigma(Element%SHI_inelastic(ipart)%Per_shell(k,m),  1.0d24)    ! module "CS_general_tools"
                           write(FN,'(es,es,es,es)') Element%SHI_inelastic(ipart)%E(m), sigma, Element%SHI_inelastic(ipart)%MFP(k,m), Element%SHI_inelastic(ipart)%Se(k,m)
                        enddo ! m = 1, Ngrid
                     endif ! (file_exist)
                     close(FN)
                     ! Normalize to the material density:
                     ! Take into account the density of atoms of this particular kind:
                     elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
                     Element%SHI_inelastic(ipart)%MFP(k,:) = Element%SHI_inelastic(ipart)%MFP(k,:) * 1.0d24/(Material(i)%At_Dens * elem_contrib)
                     Element%SHI_inelastic(ipart)%Se(k,:) = Element%SHI_inelastic(ipart)%Se(k,:) * Material(i)%At_Dens/1.0d24 * elem_contrib
                  endif VAL
                  
               enddo ! k = 1, N_shells
            enddo LMNT
            
            ! And the total CS:
            if (.not.allocated(Material(i)%SHI_inelastic_total)) allocate(Material(i)%SHI_inelastic_total(N_types_SHI))
            Nsiz = size(Material(i)%Ph_absorption_valent%E)   ! valence band has a different grid from core shells
            allocate(Material(i)%SHI_inelastic_total(ipart)%E(Nsiz))
            Material(i)%SHI_inelastic_total(ipart)%E = Material(i)%Ph_absorption_valent%E*ph_to_ion   ! [eV] -> [keV] energy grid for SHIs
            allocate(Material(i)%SHI_inelastic_total(ipart)%Total(Nsiz))
            allocate(Material(i)%SHI_inelastic_total(ipart)%Total_MFP(Nsiz))
            allocate(Material(i)%SHI_inelastic_total(ipart)%Total_Se(Nsiz))
            Material(i)%SHI_inelastic_total(ipart)%Total(:) = 0.0d0
            Material(i)%SHI_inelastic_total(ipart)%Total_MFP(:) = 0.0d0
            Material(i)%SHI_inelastic_total(ipart)%Total_Se(:) = 0.0d0
            ! Sum up CSs from each element:
            do j =1, N_elements	! for each element
               Element => Material(i)%Elements(j)	! all information about this element
               elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
               do m = 1, Nsiz	! for all energy grid points:
                  E => Material(i)%SHI_inelastic_total(ipart)%E(m)    ! electron energy [eV]
                  ! Core shells:
                  sigma = 0.0d0 ! to start with
                  Se = 0.0d0
                  lambda = 0.0d0
                  ! Get the cross section for this element on this grid point (which may not coincide with the grid for this element due to element-specific grid points):
                  do k = 1, Element%N_shl		! for all shells in this elements
                     if (.not.Element%valent(k)) then ! core orbitals
                        ! For each shell, find its partial cross section:
                        call interpolate_data_single(Element%SHI_inelastic(ipart)%E,  Element%SHI_inelastic(ipart)%Per_shell(k,:), E, sigma_cur) ! module "Little_subroutines"
                        call interpolate_data_single(Element%SHI_inelastic(ipart)%E,  Element%SHI_inelastic(ipart)%Se(k,:), E, Se_cur) ! module "Little_subroutines"
                        call interpolate_data_single(Element%SHI_inelastic(ipart)%E,  Element%SHI_inelastic(ipart)%MFP(k,:), E, lambda_cur) ! module "Little_subroutines"
                     else    ! exclude valent shells
                        sigma_cur = 0.0d0
                        Se_cur = 0.0d0
                        lambda_cur = 1.0d30
                     endif
                     ! Summing over all shells:
                     sigma = sigma + sigma_cur*elem_contrib
                     Se = Se + Se_cur !*elem_contrib
                     lambda = lambda + 1.0d0/lambda_cur !*elem_contrib
                  enddo ! k
                  Material(i)%SHI_inelastic_total(ipart)%Total(m) = Material(i)%SHI_inelastic_total(ipart)%Total(m) + sigma   ! [A^2]
                  Material(i)%SHI_inelastic_total(ipart)%Total_MFP(m) = Material(i)%SHI_inelastic_total(ipart)%Total_MFP(m) + lambda  ! [1/A] to be inversed below
                  Material(i)%SHI_inelastic_total(ipart)%Total_Se(m) = Material(i)%SHI_inelastic_total(ipart)%Total_Se(m) + Se   ! [eV/A]

                  ! Add valence band / shells:
                  !if ((numpar%El_inelast /= 2) .and. (allocated(Material(i)%CDF_valence%A)) ) then    ! valence band is defined and used in the CS model (RBEB excluded)
                  if (allocated(Material(i)%CDF_valence%A)) then    ! valence band is defined and used in the CS model
                     if (j == 1) then ! add valence band only once - when studying the first element of the material
                        Material(i)%SHI_inelastic_total(ipart)%Total(m) = Material(i)%SHI_inelastic_total(ipart)%Total(m) + Material(i)%SHI_inelastic_valent(ipart)%Total(m)   ! [A^2]
                        Material(i)%SHI_inelastic_total(ipart)%Total_MFP(m) = Material(i)%SHI_inelastic_total(ipart)%Total_MFP(m) + 1.0d0/Material(i)%SHI_inelastic_valent(ipart)%Total_MFP(m)   ! [1/A]
                        Material(i)%SHI_inelastic_total(ipart)%Total_Se(m) = Material(i)%SHI_inelastic_total(ipart)%Total_Se(m) + Material(i)%SHI_inelastic_valent(ipart)%Total_Se(m)   ! [eV/A]
                     endif
                  else ! valent atomic levels
                     sigma = 0.0d0 ! to start with
                     Se = 0.0d0
                     lambda = 0.0d0
                     ! Get the cross section for this element on this grid point (which may not coincide with the grid for this element due to element-specific grid points):
                     do k = 1, Element%N_shl		! for all shells in this elements
                        if (Element%valent(k)) then ! core orbitals
                           ! For each shell, find its partial cross section:
                           call interpolate_data_single(Element%SHI_inelastic(ipart)%E,  Element%SHI_inelastic(ipart)%Per_shell(k,:), E, sigma_cur) ! module "Little_subroutines"
                           call interpolate_data_single(Element%SHI_inelastic(ipart)%E,  Element%SHI_inelastic(ipart)%Se(k,:), E, Se_cur) ! module "Little_subroutines"
                           call interpolate_data_single(Element%SHI_inelastic(ipart)%E,  Element%SHI_inelastic(ipart)%MFP(k,:), E, lambda_cur) ! module "Little_subroutines"
                        else    ! exclude valent shells
                           sigma_cur = 0.0d0
                           Se_cur = 0.0d0
                           lambda_cur = 1.0d30
                        endif
                        ! Summing over all shells:
                        sigma = sigma + sigma_cur*elem_contrib
                        Se = Se + Se_cur !*elem_contrib
                        lambda = lambda + 1.0d0/lambda_cur !*elem_contrib
                     enddo ! k
                     Material(i)%SHI_inelastic_total(ipart)%Total(m) = Material(i)%SHI_inelastic_total(ipart)%Total(m) + sigma   ! [A^2]
                     Material(i)%SHI_inelastic_total(ipart)%Total_MFP(m) = Material(i)%SHI_inelastic_total(ipart)%Total_MFP(m) + lambda  ! [1/A] to be inversed below
                     Material(i)%SHI_inelastic_total(ipart)%Total_Se(m) = Material(i)%SHI_inelastic_total(ipart)%Total_Se(m) + Se   ! [eV/A]
                  endif
               enddo ! m = 1, Nsiz
            enddo !  j =1, N_elements
            ! Inverse to get MFP:
            Material(i)%SHI_inelastic_total(ipart)%Total_MFP = 1.0d0/Material(i)%SHI_inelastic_total(ipart)%Total_MFP  ! [1/A] -> [A]
            File_name = trim(adjustl(Path_valent))//path_sep//trim(adjustl(Ion_name))//'_in_'//trim(adjustl(Material(i)%Name))//trim(adjustl(m_ion_inelast_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
            inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
            if (.not.file_exist .or. numpar%recalculate_MFPs) then	! only create it if file does not exist
               open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
               do m = 1, Ngrid	! for all energy grid points:
                  write(FN,'(es,es,es,es)') Material(i)%SHI_inelastic_total(ipart)%E(m), Material(i)%SHI_inelastic_total(ipart)%Total(m), &
                                               Material(i)%SHI_inelastic_total(ipart)%Total_MFP(m), Material(i)%SHI_inelastic_total(ipart)%Total_Se(m)
               enddo ! m = 1, Nsiz
               close(FN)
            endif
            
         enddo TRGT
        
        write(*, '(a)') ' Done.'
!          print*, 'Ion inelastic scattering cross sections for '// trim(adjustl(Ion_name)) //' are obtained.'
      endif ! (bunch(ibunch)%KOP == 3)
   enddo BNCH
   
   nullify(path_sep, E, Element)
9999 continue
end subroutine get_ion_IMFP


! Interface to select the model of SHI inelastic scattering cross section:
subroutine get_SHI_inelastic_CS(Ekin, ZSHI, MSHI_amu, Material, DOS, numpar, bunch, ibunch, j, k, Ip, sigma, Se)
   real(8), intent(in) :: Ekin ! [eV] SHI kinetic energy
   integer, intent(in) :: ZSHI  ! atomic number of SHI
   real(8), intent(in) :: MSHI_amu  ! [amu] SHI mass
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   integer, intent(in) :: ibunch    ! number of kinds of SHI
   integer, intent(in) :: j, k  ! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   real(8), intent(out) :: sigma	! [A^2] cross section
   real(8), intent(out), optional :: Se   ! [eV/A] energy loss
   !---------------------------------------
   real(8) :: Zeff, MSHI, max_E0, Ekinq, Mc2, mtc2, Se1
   integer :: SHI_inelast
   MSHI = MSHI_amu*g_amu    ! [kg]
   ! Effective ion charge:
   Zeff = Equilibrium_charge_SHI(Ekin, MSHI, dble(ZSHI), Material%Mean_Z, numpar%SHI_ch_st, bunch(ibunch)%Zeff) ! module "SHI_charge_state"
!    SHI_inelast = numpar%El_inelast
   SHI_inelast = numpar%SHI_inelast
   
   ! Check if the ion energy is too low to use delta-CDF:
   select case (SHI_inelast)    ! only if delta-CDF is used
   case (1:3,5) ! CDF with delta-functions
       VAL0:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         max_E0 = maxval(Material%CDF_valence%E0(:))
      else VAL0
         max_E0 = maxval(Material%Elements(j)%CDF(k)%E0(:))
      endif VAL0
      ! To find the limit of aplicability, we need the rest energies of SHI and electron:
      Mc2 = rest_energy(MSHI)   ! module "Relativity"
      mtc2 = rest_energy(g_me)   ! module "Relativity"
      call find_Wmax_equal_Wmin(Mc2, mtc2, .false., Ekin, Ip, max_E0, Ekinq)   ! module "CS_integration_limits"
      
      ! Check if delta-functional CDF works here:
!       if (Ekin < m_Eeq_factor) then   ! switch to nonrelativistic numerically-integrable CDF:
      if (Ekin < numpar%CDF_Eeq_factor * Ekinq) then   ! switch to nonrelativistic numerically-integrable CDF:
          SHI_inelast = 4
      endif
!       print*, 'm', Ekin, numpar%CDF_Eeq_factor*Ekinq, SHI_inelast
   endselect
   
   select case (SHI_inelast)
   case default ! exculde
      sigma = 0.0d0
      if (present(Se)) Se = 0.0d0
   
   case (1:3,5) ! CDF with delta-functions
      VAL2:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         sigma = CDF_total_CS_delta(SHI_inelast, Ekin, MSHI, Zeff, Ip, Material%At_Dens, Material%CDF_valence, g_me, .false.)    ! module "CDF_delta"
         if (present(Se)) Se = energy_loss_delta(SHI_inelast, Ekin, MSHI, dble(ZSHI), Zeff, Ip, Material%At_Dens, g_me, Material%CDF_valence, .false.) ! module "CDF_delta"
      else VAL2  ! core shells:
         sigma = CDF_total_CS_delta(SHI_inelast, Ekin, MSHI, Zeff, Ip, Material%At_Dens, Material%Elements(j)%CDF(k), g_me, .false.)    ! module "CDF_delta"
         if (present(Se)) Se = energy_loss_delta(SHI_inelast, Ekin, MSHI, dble(ZSHI), Zeff, Ip, Material%At_Dens, g_me, Material%Elements(j)%CDF(k), .false.) ! module "CDF_delta"
      endif VAL2
   
   case (4) ! non-relativistic CDF Ritchie
      VAL3:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers)  ! module "CS_electrons_inelastic"                     
         end select ! (numpar%CDF_dispers) 
      else VAL3 ! core shells
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .false.,0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers)  ! module "CS_electrons_inelastic"
         end select ! (numpar%CDF_dispers) 
      endif VAL3
      if (present(Se)) Se = Se1   
   end select
end subroutine  get_SHI_inelastic_CS



! Interface to select the model of electron inelastic scattering cross section:
function get_SHI_inelastic_energy_transfer(Ekin, ZSHI, MSHI_amu, Material, DOS, numpar, bunch, ibunch, j, k, Ip, CS_total, CDF_dispers, CDF_m_eff, hw_phonon) result(dE)
   real(8) :: dE   ! [eV] sampled transferred energy
   real(8), intent(in) :: Ekin	! [eV] SHI kinetic energy
   integer, intent(in) :: ZSHI  ! atomic number of SHI
   real(8), intent(in) :: MSHI_amu  ! [amu] SHI mass
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   integer, intent(in) :: ibunch    ! number of kinds of SHI
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   real(8), intent(in) :: CS_total  ! [A^2] precalculated total cross section
   integer, intent(in), optional :: CDF_dispers, CDF_m_eff	! dispersion relation and effective mass
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency, for particle-phonon scattering
   !---------------------------------------
   real(8) :: Zeff, MSHI, Mc2, mtc2, Se1
   real(8) :: eps, eps2, eps3, eps_cur
   real(8) :: max_E0, Ekinq, RN, CS_sampled, CS_cur, E_left, E_right, E_cur, sigma
   real(8) :: CS_tot	! [A^2] precalculated total cross section
   integer :: dispers, m_eff, SHI_inelast, i
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-2
   eps2 = 1.0d-4
   eps3 = 1.0d-6
   
   ! Set the model parameters:
   MSHI = MSHI_amu*g_amu    ! [kg]
   ! Effective ion charge:
   Zeff = Equilibrium_charge_SHI(Ekin, MSHI, dble(ZSHI), Material%Mean_Z, numpar%SHI_ch_st, bunch(ibunch)%Zeff) ! module "SHI_charge_state"
!    SHI_inelast = numpar%El_inelast  ! to chose the model of CS calculations below
   SHI_inelast = numpar%SHI_inelast  ! to chose the model of CS calculations below
   if (present(CDF_dispers)) then	! user provided alternative (e.g. for all but VB)
      dispers = CDF_dispers
   else	! for the VB
      dispers = numpar%CDF_dispers
   endif
   if (present(CDF_m_eff)) then	! user provided alternative (e.g. for all but VB)
      m_eff = CDF_m_eff
   else
      m_eff = numpar%CDF_m_eff
   endif
   
   Mc2 = rest_energy(MSHI)   ! module "Relativity"
   mtc2 = rest_energy(g_me)   ! module "Relativity"

   ! In case it is delta-CDF model, we have to make sure the energy is not below its applicability limit:
   select case (SHI_inelast)    ! only if delta-CDF is used
   case (1:3,5) ! CDF with delta-functions
      !VAL0:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
      VAL0:if ( ( k == 0 ) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         max_E0 = maxval(Material%CDF_valence%E0(:))
      else VAL0
         max_E0 = maxval(Material%Elements(j)%CDF(k)%E0(:))
      endif VAL0
      call find_Wmax_equal_Wmin(Mc2, mtc2, .false., Ekin, Ip, max_E0, Ekinq)   ! module "CS_integration_limits"
      ! Define minimal and maximal energies (integration limits):
      E_left = W_min(Ip, Mc2, mtc2, Ekin, Ip) ! module "CS_integration_limits"
      !E_right = W_max(Mc2, mtc2, .false., Ekin, Ip)  ! module "CS_integration_limits"
      if (present(hw_phonon)) then  ! scattering on phonons:
         E_right = W_max(Mc2, mtc2, .false., Ekin, Ip, hw_phonon)  ! module "CS_integration_limits"
      else  ! scattering on electrons or atoms:
         E_right = W_max(Mc2, mtc2, .false., Ekin, Ip)  ! module "CS_integration_limits"
      endif
      
      ! Check if delta-functional CDF works here:
!       if (Ekin < m_Eeq_factor) then   ! switch to nonrelativistic numerically-integrable CDF:
      if (Ekin < numpar%CDF_Eeq_factor * Ekinq) then   ! switch to nonrelativistic numerically-integrable CDF:
         SHI_inelast = 4
         ! Define minimal and maximal energies (integration limits):
!          E_left = Ip
         E_left = W_min(Ip, Mc2, mtc2, Ekin, Ip) ! module "CS_integration_limits"
         !E_right = W_max(Mc2, mtc2, .false., Ekin, Ip)  ! module "CS_integration_limits"
         E_right = W_max_nonrel_free(Mc2, mtc2, .false., Ekin, Ip)    ! module "CS_integration_limits" 
      endif

   case default
!       E_left = Ip
      E_left = W_min(Ip, Mc2, mtc2, Ekin, Ip) ! module "CS_integration_limits"
      E_right = W_max_nonrel_free(Mc2, mtc2, .false., Ekin, Ip)    ! module "CS_integration_limits" 
   end select
   
!    E_left = W_min(Ip, Mc2, mtc2, Ekin, 0.0d0) ! module "CS_integration_limits"
!    E_right = W_max(Mc2, mtc2, .false., Ekin, 0.0d0)  ! module "CS_integration_limits"
   ! Use of the maximal plasmon energy as the upper limit of integration is not included yet !!!
   
   ! Sample the cross section:
   call random_number(RN)
   
!    do i = 1, 200    ! Testing
!    RN = dble(i)/200 ! Testing
   
!    CS_tot = get_SHI_integral_inelastic_CS(Ekin, Material, DOS, numpar, j, k, Ip, SHI_inelast, MSHI, Zeff, dispers, m_eff, E_right)   ! below
!    print*, 'SHI', k, CS_tot, CS_total
   CS_tot = CS_total    ! use precalculated value
   
   CS_sampled = RN*CS_tot
   
   ! For analytucally integrable CDF, use bisection method:
   if (SHI_inelast /= 4) then
   
      ! Start finding CS:
      E_cur = (E_left + E_right)*0.5d0
      CS_cur = get_SHI_integral_inelastic_CS(Ekin, Material, DOS, numpar, j, k, Ip, SHI_inelast, MSHI, Zeff, dispers, m_eff, E_cur)   ! below
   
!        call get_SHI_inelastic_CS(Ekin, ZSHI, MSHI_amu, Material, DOS, numpar, bunch, ibunch, j, k, Ip, sigma)	! see below
!       write(*,'(f,f,f,f,f,f)') Ekin, E_cur, CS_cur, CS_sampled, CS_tot, sigma
   
      ! Sampled cross-sections closer to the total one require higher precision:
      if (RN > 0.99d0) then
         eps_cur = eps3
      elseif (RN > 0.8d0) then
         eps_cur = eps2
      else
         eps_cur = eps
      endif
   
      ! Search by bisection method:
      !do while (ABS(CS_cur - CS_sampled) >= eps_cur*CS_sampled)
      do while (abs(E_left - E_right) >= eps*E_left)
         if (CS_cur > CS_sampled) then
            E_right = E_cur
         else
            E_left = E_cur
         endif
         E_cur = (E_left + E_right)*0.5d0
         if (abs(E_left - E_right) < eps) exit  ! precise enough
         CS_cur = get_SHI_integral_inelastic_CS(Ekin, Material, DOS, numpar, j, k, Ip, SHI_inelast, MSHI, Zeff, dispers, m_eff, E_cur)   ! below
!          if (RN > 0.99d0) write(*,'(f,f,f,f,f,f)') Ekin, E_cur, CS_cur, CS_sampled, CS_tot, RN
      enddo

   else     ! for numerically integrable one, just integrate until found the sampled value:
   
      VAL3:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                    Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, Sigma_sampled=CS_sampled, &
                     E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"                     
         end select ! (numpar%CDF_dispers) 
      else VAL3 ! core shells
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                    Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .false.,0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         end select ! (numpar%CDF_dispers) 
      endif VAL3
!         print*, 'get_SHI_inelastic_energy_transfer 2', CS_sampled, CS_cur, E_cur, RN, CS_tot
    endif    ! (El_inelast /= 4)

   ! Output: sampled transferred energy:
   dE = E_cur
   
!    print*, RN, dE ! Testing
!    enddo          ! Testing
!    PAUSE 'TEST get_SHI_inelastic_energy_transfer'
!    pause 'get_SHI_inelastic_energy_transfer'
!    if (dE > 150.0d0) then
!       print*, 't1', Ekin, m_Eeq_factor*Ekinq, SHI_inelast
!       print*, 't2', CS_sampled, CS_tot, RN
!       print*, 't3', dE
!    endif
end function get_SHI_inelastic_energy_transfer


function get_SHI_integral_inelastic_CS(Ekin, Material, DOS, numpar, j, k, Ip, SHI_inelast, MSHI, Zeff, dispers, m_eff, Emax) result (sigma)
   real(8) :: sigma	! [A^2] cross section
   real(8), intent(in) :: Ekin	! [eV] electron kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   integer, intent(in) :: SHI_inelast    ! model for the inelastic cross section
   real(8), intent(in) :: MSHI  ! [kg] SHI mass
   real(8), intent(in) :: Zeff  ! effective charge of the incident particle
   integer, intent(in) :: dispers, m_eff	! dispersion relation and effective mass
   real(8), intent(in) :: Emax  ! [eV] upper integration limit
   !---------------------------------------
   real(8) :: Se1

! Now chose the model of CS:
   select case (SHI_inelast)  ! chose which model for electron inelastic cross section to use
   case (1:3,5)	! CDF with delta-functions
      VAL2: if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         sigma = CDF_total_CS_delta(SHI_inelast, Ekin, MSHI, Zeff, Ip, Material%At_Dens, Material%CDF_valence, g_me, .false., Emax_in = Emax)    ! module "CDF_delta"
      elseif ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band, another definition
         sigma = CDF_total_CS_delta(SHI_inelast, Ekin, MSHI, Zeff, Ip, Material%At_Dens, Material%CDF_valence, g_me, .false., Emax_in = Emax)    ! module "CDF_delta"
      else VAL2  ! core shells:
         sigma = CDF_total_CS_delta(SHI_inelast, Ekin, MSHI, Zeff, Ip, Material%At_Dens, Material%Elements(j)%CDF(k), g_me, .false., Emax_in = Emax)    ! module "CDF_delta"
      endif VAL2
   
   case (4) ! nonrelativistic Ritchie CDF
      VAL3:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band, another definition
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                    Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, Wmax_in=Emax) ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, & 
                     CDF_dispers=numpar%CDF_dispers, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         end select ! (numpar%CDF_dispers)
      else VAL3 ! core shells
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                    Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .false.,0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, Wmax_in=Emax) ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ekin, MSHI, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         end select ! (numpar%CDF_dispers) 
         
      endif VAL3
   case default	! exclude
      sigma = 0.0d0
   end select
end function get_SHI_integral_inelastic_CS




end module CS_ion_inelastic
