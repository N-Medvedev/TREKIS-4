! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018-2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains all cross sections for electrons inelastic scattering 
! References used:
! [1] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2008: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2009)
! [2] F. Salvat, J.M. Fernhdez-Varea, NIMB 63 (1992) 255-269
! [3] I. Plante, F,Cucinotta, New Journal of Physics 11 (2009) 063047
! [4] "Transport of Electrons and Photons"  Edited by Jenkins, Nelson, Rindi
! [5] Kim, Santos, Parente, PRB 62 (2000) 052710
! [6] R. Rymzhanov et al., Nucl. Instrum. and Methods B 388 (2016) 41


module CS_electrons_inelastic
use Universal_constants
use CDF_Ritchi, only: Diel_func, Find_Ritchie_peak_Q
use CDF_Mermin, only: Mermin_Diel_func
use CDF_delta, only: integrated_delta_CDF_CS, P_prefactor, Find_linear_a_b, energy_loss_delta, CDF_total_CS_delta
use Relativity
use Objects
use Read_input_data, only: m_electron_CS, m_input_folder, m_electron_inelast_CS, m_folder_materials
use Dealing_with_files, only: read_file
use CS_general_tools, only: MFP_from_sigma, plasmon_energy, remap_energy_step, temperature_factor, remap_dq
use Little_subroutines, only: interpolate_data_single, find_in_array_monoton, print_time_step
use CS_integration_limits, only: W_max, Q_min_nonrel, Q_max_nonrel, find_Wmax_equal_Wmin, W_min, W_max_nonrel_free
use SHI_charge_state, only: Equilibrium_charge_SHI

implicit none

real(8) :: m_Eeq_factor
! Used for testing purposes:
parameter (m_Eeq_factor = 3200.0d0)    ! Ekin of ion [eV] -- Testing

 contains
 
!-----------------------------------------------------------------------------
! Subroutines for inelastic scattering:

subroutine get_electron_IMFP(Material, numpar, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material	!material parameters of each target that it's constructed of
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------------------------
   integer :: i, j, k, m, Ngrid, N_targets, N_elements, N_shells, FN, Reason, count_lines, Nsiz
   real(8) :: sigma, Se, lambda, sigma_cur, Se_cur, lambda_cur
   real(8) :: Total_Se, Total_lambda, N_elem, elem_contrib
   character(200) :: Path, Folder_with_CS, command, File_name, Model_name
   character(200) :: Path_valent
   character(10) :: temp_c, temp_c2, temp_c3
   logical :: file_exist, read_well
   !--------------------------------
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   
   write(*, '(a)', advance='no') ' Obtaining electron inelastic scattering cross sections...'
   
   path_sep => numpar%path_sep
   Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_electron_CS))	! Electron CSs are storred in this folder
   
   N_targets = size(Material)	! that's how many different targets user specified
   
   ! Get electron MFPs:
   TRGT:do i = 1, N_targets	! for each target
      N_elem = dble(SUM(Material(i)%Elements(:)%percentage))	! number of elements in this compound material
      N_elements = size(Material(i)%Elements)	! that's how many different elements are in this target
      LMNT:do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         ! Check if folder for this element already exists:
         Folder_with_CS = trim(adjustl(Path))//path_sep//trim(adjustl(Element%Name))
         inquire(DIRECTORY=trim(adjustl(Folder_with_CS)),exist=file_exist)    ! check if input file excists
         if (.not.file_exist) then	! to make sure that such a folder is present (even if empty)
            ! Create a new directory for output files:
            command='mkdir '//trim(adjustl(Folder_with_CS))	! to create a folder use this command
            CALL system(command)	! create the folder
         endif
         !-------------------------------------
         ! Set part of the file name corresponding to the model used for electron inelastic CS:
         write(temp_c,'(f8.1)') Material(i)%T    ! target temperature
         write(temp_c2,'(f8.1)') numpar%CDF_Eeq_factor
         select case (numpar%CDF_model)
         case default   ! Ritchie
            write(temp_c3,'(a)') '_R_'    ! Ritchie
         case (1)   ! Mermin
            write(temp_c3,'(a)') '_M_'    ! Mermin
         endselect
         select case (numpar%El_inelast)
         case (1)	! CDF with Ritchi's oscillators
            Model_name = 'CDF_Ritchi_T_'//trim(adjustl(temp_c))//'K'
         case (2)	! RBEB
            Model_name = 'RBEB'
         case (3)	! CDF with delta-functions
            Model_name = 'CDF_delta_T_'//trim(adjustl(temp_c))//'K'//trim(adjustl(temp_c3))//'Eeq_'//trim(adjustl(temp_c2))
         case (5)	! CDF with single pole delta-functions
            Model_name = 'CDF_SP_delta_T_'//trim(adjustl(temp_c))//'K'//trim(adjustl(temp_c3))//'Eeq_'//trim(adjustl(temp_c2))
         case default	! exclude
            Model_name = 'NO'
         end select

         ! Path to files with the valence band and total (material specific) MFPs:
         !if ((numpar%El_inelast /= 2) .and. (allocated(Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
            Path_valent = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
            Path_valent = trim(adjustl(Path_valent))//path_sep//trim(adjustl(Material(i)%Name))
         !endif
         
         ! For electron CSs use the same grid as for photons:
         Ngrid = size(Element%Phot_absorption%E)
         allocate(Element%El_inelastic%E(Ngrid))
         Element%El_inelastic%E = Element%Phot_absorption%E
         allocate(Element%El_inelastic%Total(Ngrid))
         allocate(Element%El_inelastic%Total_MFP(Ngrid))
         allocate(Element%El_inelastic%Total_Se(Ngrid))
         allocate(Element%El_inelastic%Per_shell(Element%N_shl,Ngrid))
         allocate(Element%El_inelastic%MFP(Element%N_shl,Ngrid))
         allocate(Element%El_inelastic%Se(Element%N_shl,Ngrid))
         N_shells = Element%N_shl
         
         ! Calculate total cross sections for all shells of this element:
         do k = 1, N_shells
            ! Check valence band:
            ! If this shell forms valent band, and data for the valence band exists, separate it from the atomic shells:
            VAL:if ((numpar%El_inelast /= 2) .and. (Element%valent(k)) .and. (allocated(Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
               if (.not.allocated(Material(i)%El_inelastic_valent%E)) then   ! valence band has not been defined yet, so do that
                  Nsiz = size(Material(i)%Ph_absorption_valent%E)   ! valence band has a different grid from core shells
                  allocate(Material(i)%El_inelastic_valent%E(Nsiz))
                  Material(i)%El_inelastic_valent%E = Material(i)%Ph_absorption_valent%E    ! we already set this grid for photons, reuse it
                  allocate(Material(i)%El_inelastic_valent%Total(Nsiz))
                  allocate(Material(i)%El_inelastic_valent%Total_MFP(Nsiz))
                  allocate(Material(i)%El_inelastic_valent%Total_Se(Nsiz))
                  !allocate(Material(i)%El_inelastic_valent%Per_shell(Element%N_shl,Nsiz)) ! --- leave it not allocated !!!
                  ! Check file with CSs:
                  !File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_electron_inelast_CS))//'_valence_'//trim(adjustl(Model_name))//'.dat'
                  File_name = trim(adjustl(Path_valent))//path_sep//trim(adjustl(m_electron_inelast_CS))//'_valence_'//trim(adjustl(Model_name))//'.dat'
                  inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
                  if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then   ! just read from the file, no need to recalculate:
                     open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
                     ! Get the MFP and CDF points:
                     count_lines = 0
                     do m = 1, Nsiz	! for all energy grid points:
                        E => Material(i)%El_inelastic_valent%E(m)	! electron energy [eV]
                        read(FN,'(es,es,es,es)', IOSTAT=Reason) E, Material(i)%El_inelastic_valent%Total(m), Material(i)%El_inelastic_valent%Total_MFP(m), Material(i)%El_inelastic_valent%Total_Se(m)
                        call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                        if (.not.read_well) then
                           close(FN)	! redo the file
                           goto 8391
                        endif
                     enddo
                  else	! no such file (or user wants to recalculate it) => create it
8391                 open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
                     do m = 1, Nsiz	! for all energy grid points:
                        E => Material(i)%El_inelastic_valent%E(m)	! electron energy [eV]
                        !call get_el_inelastic_CS(E, Material(i), Material(i)%DOS, numpar, j, k, Material(i)%DOS%Egap, sigma, CDF_dispers=0, CDF_m_eff=1, Se=Se)	! see below
                        call get_el_inelastic_CS(E, Material(i), Material(i)%DOS, numpar, j, k, Material(i)%DOS%Egap, sigma, &
                                                          CDF_dispers=numpar%CDF_dispers, CDF_m_eff=numpar%CDF_m_eff, Se=Se)    ! see below
                        Material(i)%El_inelastic_valent%Total(m) = sigma   ! [A^2]
                        lambda = MFP_from_sigma(sigma,  Material(i)%At_Dens)    ! module "CS_general_tools"
                        Material(i)%El_inelastic_valent%Total_MFP(m) = lambda  ! [A]
                        Material(i)%El_inelastic_valent%Total_Se(m) = Se    ! [eV/A]
                        write(FN,'(es,es,es,es)') Material(i)%El_inelastic_valent%E(m), sigma, lambda, Se
!                         write(*,'(a,es,es,es,es)') 'EE:', E, sigma, lambda, Material(i)%DOS%Egap
                     enddo ! m = 1, Nsiz
                  endif ! (file_exist)
                  close(FN)
               endif ! (.not.allocated(Material(i)%El_inelastic_valent))
            else VAL    ! core shell
               ! Check file with CSs:
               File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_electron_inelast_CS))//'_'//trim(adjustl(Element%Shell_name(k)))//'_'//trim(adjustl(Model_name))//'.dat'
               inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
               if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then  ! just read from the file, no need to recalculate:
                  open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
                  ! Get the MFP and CDF points:
                  count_lines = 0
                  do m = 1, Ngrid	! for all energy grid points:
                     E => Element%El_inelastic%E(m)	! electron energy [eV]
                     read(FN,'(es,es,es,es)', IOSTAT=Reason) E, Element%El_inelastic%Per_shell(k,m), Element%El_inelastic%MFP(k,m),  Element%El_inelastic%Se(k,m)
                     call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                     if (.not.read_well) then
                        close(FN)	! redo the file
                        goto 8390
                     endif
!                       print*, k, m, Element%El_inelastic%E(m), Element%El_inelastic%Per_shell(k,m)
                  enddo
               else	! no such file => create it
8390           open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
                  do m = 1, Ngrid	! for all energy grid points:
                     E => Element%El_inelastic%E(m)	! electron energy [eV]
                     call get_el_inelastic_CS(E, Material(i), Material(i)%DOS, numpar, j, k, Material(i)%Elements(j)%Ip(k), sigma, &
                                                          CDF_dispers=numpar%CDF_dispers, CDF_m_eff=numpar%CDF_m_eff, Se=Se)    ! see below
                     Element%El_inelastic%Per_shell(k,m) = sigma	! [A^2]
                     ! MFP and Se are calculated per 1.0d24 1/cm^3 density; need to be rescaled below:
                     Element%El_inelastic%Se(k,m) = Se * 1.0d24/Material(i)%At_Dens    ! [eV/A]
                     Element%El_inelastic%MFP(k,m) = MFP_from_sigma(Element%El_inelastic%Per_shell(k,m), 1.0d24)    ! module "CS_general_tools"
                     write(FN,'(es,es,es,es)') Element%El_inelastic%E(m), sigma, Element%El_inelastic%MFP(k,m),  Element%El_inelastic%Se(k,m)
                  enddo ! m = 1, Ngrid
               endif ! (file_exist)
               close(FN)
               
               ! Normalize MFPs and Se for real material density:
               ! Take into account the density of atoms of this particular kind:
               elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
               Element%El_inelastic%MFP(k,:) = Element%El_inelastic%MFP(k,:) * 1.0d24/(Material(i)%At_Dens * elem_contrib)
               Element%El_inelastic%Se(k,:) = Element%El_inelastic%Se(k,:) * Material(i)%At_Dens/1.0d24 * elem_contrib
            endif VAL
         enddo ! k = 1, N_shells
      enddo LMNT
     
      ! And the total CS:
      Nsiz = size(Material(i)%Ph_absorption_total%E)   ! grid for total cross section is different from core-shells grid (but the same as valence grid)
      allocate(Material(i)%El_inelastic_total%E(Nsiz))
      Material(i)%El_inelastic_total%E = Material(i)%Ph_absorption_total%E
      allocate(Material(i)%El_inelastic_total%Total(Nsiz))
      allocate(Material(i)%El_inelastic_total%Total_MFP(Nsiz))
      allocate(Material(i)%El_inelastic_total%Total_Se(Nsiz))
      Material(i)%El_inelastic_total%Total(:) = 0.0d0
      Material(i)%El_inelastic_total%Total_MFP(:) = 0.0d0
      Material(i)%El_inelastic_total%Total_Se(:) = 0.0d0
      ! Sum up CSs from each element:
      do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         do m = 1, Nsiz	! for all energy grid points:
            E => Material(i)%El_inelastic_total%E(m)    ! electron energy [eV]
            ! Core shells:
            sigma = 0.0d0 ! to start with
            Se = 0.0d0
            lambda = 0.0d0
            ! Get the cross section for this element on this grid point (which may not coincide with the grid for this element due to element-specific grid points):
            do k = 1, Element%N_shl		! for all shells in this elements
               if (.not.Element%valent(k)) then ! core orbitals
                  ! For each shell, find its partial cross section:
                  call interpolate_data_single(Element%El_inelastic%E,  Element%El_inelastic%Per_shell(k,:), E, sigma_cur) ! module "Little_subroutines"
                  call interpolate_data_single(Element%El_inelastic%E,  Element%El_inelastic%Se(k,:), E, Se_cur) ! module "Little_subroutines"
                  call interpolate_data_single(Element%El_inelastic%E,  Element%El_inelastic%MFP(k,:), E, lambda_cur) ! module "Little_subroutines"
               else    ! exclude valent shells
                  sigma_cur = 0.0d0
                  Se_cur = 0.0d0
                  lambda_cur = 1.0d30
               endif
               ! Summing over all shells:
               sigma = sigma + sigma_cur*elem_contrib
               Se = Se + Se_cur !*elem_contrib
               lambda = lambda + 1.0d0/lambda_cur !*elem_contrib
!                write(*,'(a,i2,i2,i2,e,e)') 'MFP', j, k, Element%N_shl, E, sigma
            enddo ! k
            Material(i)%El_inelastic_total%Total(m) = Material(i)%El_inelastic_total%Total(m) + sigma   ! [A^2]
            Material(i)%El_inelastic_total%Total_Se(m) = Material(i)%El_inelastic_total%Total_Se(m) +  Se   ! [eV/A]
            Material(i)%El_inelastic_total%Total_MFP(m) = Material(i)%El_inelastic_total%Total_MFP(m) + lambda  ! [1/A] to be inversed below
            
            ! Add valence band / shells:
            if ((numpar%El_inelast /= 2) .and. (allocated(Material(i)%CDF_valence%A)) ) then    ! valence band is defined and used in the CS model (RBEB excluded)
               if (j == 1) then ! add valence band only once - when studying the first element of the material
                  Material(i)%El_inelastic_total%Total(m) = Material(i)%El_inelastic_total%Total(m) + &
                                                    Material(i)%El_inelastic_valent%Total(m)   ! [A^2]
                  Material(i)%El_inelastic_total%Total_MFP(m) = Material(i)%El_inelastic_total%Total_MFP(m) + &
                                                    1.0d0/Material(i)%El_inelastic_valent%Total_MFP(m)   ! [1/A]
                  Material(i)%El_inelastic_total%Total_Se(m) = Material(i)%El_inelastic_total%Total_Se(m) + &
                                                    Material(i)%El_inelastic_valent%Total_Se(m)   ! [eV/A]
               endif
            else ! valent atomic levels
               ! Valence shells:
               sigma = 0.0d0 ! to start with
               Se = 0.0d0
               lambda = 0.0d0
               ! Get the cross section for this element on this grid point (which may not coincide with the grid for this element due to element-specific grid points):
               do k = 1, Element%N_shl		! for all shells in this elements
                  if (Element%valent(k)) then ! core orbitals
                     ! For each shell, find its partial cross section:
                     call interpolate_data_single(Element%El_inelastic%E,  Element%El_inelastic%Per_shell(k,:), E, sigma_cur) ! module "Little_subroutines"
                     call interpolate_data_single(Element%El_inelastic%E,  Element%El_inelastic%Se(k,:), E, Se_cur) ! module "Little_subroutines"
                     call interpolate_data_single(Element%El_inelastic%E,  Element%El_inelastic%MFP(k,:), E, lambda_cur) ! module "Little_subroutines"
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
               Material(i)%El_inelastic_total%Total(m) = Material(i)%El_inelastic_total%Total(m) + sigma   ! [A^2]
               Material(i)%El_inelastic_total%Total_Se(m) = Material(i)%El_inelastic_total%Total_Se(m) +  Se   ! [eV/A]
               Material(i)%El_inelastic_total%Total_MFP(m) = Material(i)%El_inelastic_total%Total_MFP(m) + lambda  ! [1/A] to be inversed below
            endif
         enddo ! m = 1, Nsiz
      enddo !  j =1, N_elements
      ! Inverse it to get the MFP:
      Material(i)%El_inelastic_total%Total_MFP = 1.0d0/Material(i)%El_inelastic_total%Total_MFP  ! [1/A] -> [A]
      ! Save the total cross section into the file:
      File_name = trim(adjustl(Path_valent))//path_sep//trim(adjustl(m_electron_inelast_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!       if (.not.file_exist) then	! only create it if file does not exist
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            write(FN,'(es,es,es,es)') Material(i)%El_inelastic_total%E(m), Material(i)%El_inelastic_total%Total(m), &
                            Material(i)%El_inelastic_total%Total_MFP(m), Material(i)%El_inelastic_total%Total_Se(m)
         enddo ! m = 1, Nsiz
         close(FN)
!       endif
   enddo TRGT
   
!    PAUSE 'get_electron_IMFP'    ! Testing
   write(*,'(a)') ' Done.'
!    print*, 'Electron inelastic scattering cross sections are obtained.'
   
   nullify(path_sep, E, Element)
end subroutine get_electron_IMFP



! Interface to select the model of electron inelastic scattering cross section:
 subroutine get_el_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, sigma, CDF_dispers, CDF_m_eff, Se)
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   real(8), intent(out) :: sigma	! [A^2] cross section
   real(8), intent(out), optional :: Se   ! [eV/A] energy loss
   integer, intent(in), optional :: CDF_dispers, CDF_m_eff	! dispersion relation and effective mass
   !---------------------------------------
   real(8) :: Zeff, Se1, max_E0, Eeq
   integer :: dispers, m_eff, El_inelast
   Zeff = 1.0d0	! electron charge
   El_inelast = numpar%El_inelast   ! to chose the model of CS calculations below
   
   ! In case it is delta-CDF model, we have to make sure the energy is not below its applicability limit:
   if ((El_inelast == 3) .or. (El_inelast == 5)) then
      VAL0:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         max_E0 = maxval(Material%CDF_valence%E0(:))
      else VAL0
         max_E0 = maxval(Material%Elements(j)%CDF(k)%E0(:))
      endif VAL0
      call find_Wmax_equal_Wmin(0.0d0, 0.0d0, .true., Ee, Ip, max_E0, Eeq)   ! module "CS_integration_limits"
      ! Check if delta-functional CDF works here:
!       if (Ee < m_Eeq_factor) then   ! switch to nonrelativistic numerically-integrable CDF:
      if (Ee < numpar%CDF_Eeq_factor * Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
          El_inelast = 4
      endif
   endif ! (El_inelast == 3)
   
   ! Now chose the model of CS:
   select case (El_inelast)  ! chose which model for electron inelastic cross section to use
   case (1)	! relativistic CDF with Ritchi model
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
      VAL:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         sigma = CDF_total_CS(Ee, g_me, Zeff, Ip, Material%T_eV, Material%At_Dens, Material%CDF_valence, &
            g_me, dispers, m_eff, DOS%k, DOS%Eff_m, Material%DOS%v_f)	! see below -- relativistic version, too slow, do not use!
      else VAL  ! core shells:
         sigma = CDF_total_CS(Ee, g_me, Zeff, Ip, Material%T_eV, Material%At_Dens, Material%Elements(j)%CDF(k), &
            g_me, dispers, m_eff, DOS%k, DOS%Eff_m, Material%DOS%v_f)	! see below -- relativistic version, too slow, do not use!
      endif VAL
      
   case (2)	! RBEB
      sigma = RBEB_total_CS(Ee, Material%Elements(j)%Ne_shell(k), Ip, Material%Elements(j)%Ek(k))	! see below
   
   case (3,5)	! CDF with delta-functions
       VAL2:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         sigma = CDF_total_CS_delta(El_inelast, Ee, g_me, Zeff, Ip, Material%At_Dens, Material%CDF_valence, g_me, .true.)    ! module "CDF_delta"
         if (present(Se)) Se = energy_loss_delta(El_inelast, Ee, g_me, 1.0d0, Zeff, Ip, Material%At_Dens, &
            g_me, Material%CDF_valence, .true., 0) ! module "CDF_delta"
      else VAL2  ! core shells:
         sigma = CDF_total_CS_delta(El_inelast, Ee, g_me, Zeff, Ip, Material%At_Dens, Material%Elements(j)%CDF(k), &
         g_me, .true.)    ! module "CDF_delta"
         if (present(Se)) Se = energy_loss_delta(El_inelast, Ee, g_me, 1.0d0, Zeff, Ip, Material%At_Dens, &
            g_me, Material%Elements(j)%CDF(k), .true., 0) ! module "CDF_delta"
      endif VAL2
   
   case (4) ! nonrelativistic Ritchie CDF
      VAL3:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, &
                    DOS%k, DOS%Eff_m, .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model)  ! below
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f)  ! below
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, CDF_dispers=numpar%CDF_dispers)  ! below
         end select ! (numpar%CDF_dispers) 
      else VAL3 ! core shells
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, &
                    DOS%k, DOS%Eff_m, .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model)  ! below
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .true.,1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f)  ! below
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, CDF_dispers=numpar%CDF_dispers)  ! below
         end select ! (numpar%CDF_dispers) 
      endif VAL3
      if (present(Se)) Se = Se1
   
   case default	! exclude
      sigma = 0.0d0
   end select
end subroutine get_el_inelastic_CS



! Interface to select the model of electron inelastic scattering cross section:
function get_inelastic_energy_transfer(Ee, Material, DOS, numpar, j, k, Ip, CS_total, CDF_dispers, CDF_m_eff, hw_phonon) result(dE)
   real(8) :: dE   ! [eV] sampled transferred energy
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   real(8), intent(in) :: CS_total  ! [A^2] precalculated total cross section
   integer, intent(in), optional :: CDF_dispers, CDF_m_eff	! dispersion relation and effective mass
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency, for particle-phonon scattering
   !---------------------------------------
   real(8) :: eps, eps1, eps2, eps3
   real(8) :: Zeff, max_E0, Eeq, RN, CS_sampled, CS_cur, E_left, E_right, E_cur, Se1
   real(8) :: CS_tot	! [A^2] precalculated total cross section
   real(8) :: Mc2, mtc2
   integer :: dispers, m_eff, El_inelast
   ! Set accepteble margin of precision for the angle:
   eps1 = 1.0d-2
   eps2 = 1.0d-4
   eps3 = 1.0d-6
   
   ! Set the model parameters:
   Zeff = 1.0d0	! electron charge
   El_inelast = numpar%El_inelast   ! to chose the model of CS calculations below
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

   ! In case it is delta-CDF model, we have to make sure the energy is not below its applicability limit:
   if ( (El_inelast == 3) .or. (El_inelast == 5) ) then
      !VAL0:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
      VAL0:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         max_E0 = maxval(Material%CDF_valence%E0(:))
      else VAL0
         max_E0 = maxval(Material%Elements(j)%CDF(k)%E0(:))
      endif VAL0
      call find_Wmax_equal_Wmin(0.0d0, 0.0d0, .true., Ee, Ip, max_E0, Eeq)   ! module "CS_integration_limits"
      
      ! Check if delta-functional CDF works here:
!       if (Ee < m_Eeq_factor) then   ! switch to nonrelativistic numerically-integrable CDF:
      if (Ee < numpar%CDF_Eeq_factor * Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
         El_inelast = 4
         E_left = Ip ! [eV] minimal transferred energy
         E_right = (Ip + Ee)*0.5d0    ! [eV] maximal transferred energy
      else
         ! Define minimal and maximal energies (integration limits):
         Mc2 = rest_energy(g_me)   ! module "Relativity"
         mtc2 = Mc2 !rest_energy(g_me)   ! module "Relativity"
         E_left = W_min(Ip, Mc2, mtc2, Ee, Ip) ! module "CS_integration_limits"
         E_right = W_max(Mc2, mtc2, .true., Ee, Ip)  ! module "CS_integration_limits"      
      endif
   else
      E_left = Ip ! [eV] minimal transferred energy
      E_right = (Ip + Ee)*0.5d0    ! [eV] maximal transferred energy
   endif ! (El_inelast == 3)
   ! Use of the maximal plasmon energy as the upper limit of integration is not included yet !!!         
   
   ! Sample the cross section:
   call random_number(RN)
!    CS_tot = get_integral_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, dispers, m_eff, E_right)   ! below
!    print*, 'El_dE:', Ee, k, CS_tot, CS_total ! Testing  
   CS_tot = CS_total    ! use precalculated value
   
   CS_sampled = RN*CS_tot
        
   ! For analytucally integrable CDF, use bisection method:
   if (El_inelast /= 4) then
        !if ( (RN < 1.0d-8) .or. (CS_sampled < 1.0d-12) )then   ! no need to sample, it's just the lower limit
        if ( RN < 1.0d-8 )then   ! no need to sample, it's just the lower limit
            E_cur = E_left
            CS_cur = 0.0d0
        else ! sample it
           
           ! sampled cross-sections closer to the total one require higher precision:
           if (RN > 0.99d0) then    
              eps = eps3
           elseif (RN > 0.8d0) then
              eps = eps2
           else
              eps = eps1
           endif
        
            ! Start finding CS:
            E_cur = (E_left + E_right)*0.5d0
            CS_cur = get_integral_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, dispers, m_eff, E_cur)   ! below
            ! Search by bisection method:
!             do while (ABS(CS_cur - CS_sampled) > eps*CS_sampled)
            do while (abs(E_left - E_right) >= eps1*E_left)
                if (CS_cur > CS_sampled) then
                    E_right = E_cur
                else
                    E_left = E_cur
                endif
                E_cur = (E_left + E_right)/2.0d0
                if (abs(E_left - E_right) < eps1) exit  ! precise enough
                CS_cur = get_integral_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, dispers, m_eff, E_cur)   ! below
            enddo
        endif
        
!         print*, 'get_inelastic_energy_transfer 1', CS_sampled, CS_cur, E_cur, RN, CS_tot
    else     ! for numerically integrable one, just integrate until found the smapled value:
         VAL3:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
            select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
            case default  ! free electron
                call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, &
                    DOS%k, DOS%Eff_m, .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                    Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! below
            case (1)  ! Plasmon pole
                call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! below
            case (2)  ! Ritchie extended
                call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! below
            end select ! (numpar%CDF_dispers) 
        else VAL3 ! core shells
            select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
            case default  ! free electron
                call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, &
                    DOS%k, DOS%Eff_m, .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                    Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! below
            case (1)  ! Plasmon pole
                call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .true.,1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! below
            case (2)  ! Ritchie extended
                call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! below
            end select ! (numpar%CDF_dispers) 
        endif VAL3
!         print*, 'get_inelastic_energy_transfer 2', CS_sampled, CS_cur, E_cur, RN, CS_tot
    endif    ! (El_inelast /= 4)
   
   ! Output: sampled transferred energy:
   dE = E_cur

!    if (dE <= Ip) print*, 'get_inelastic_energy_transfer', Ee, dE, Ip, El_inelast, ( (k == 0) .and. (allocated(Material%CDF_valence%A)) )
end function get_inelastic_energy_transfer


! Interface to select the model of electron inelastic scattering cross section:
function get_inelastic_energy_transfer_OLD(Ee, Material, DOS, numpar, j, k, Ip, CDF_dispers, CDF_m_eff, hw_phonon) result(dE)
   real(8) :: dE   ! [eV] sampled transferred energy
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   integer, intent(in), optional :: CDF_dispers, CDF_m_eff	! dispersion relation and effective mass
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency, for particle-phonon scattering
   !---------------------------------------
   real(8) :: eps, Zeff, max_E0, Eeq, RN, CS_sampled, CS_cur, E_left, E_right, E_cur
   real(8) :: CS_tot	! [A^2] precalculated total cross section
   integer :: dispers, m_eff, El_inelast
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-3
   
   ! Set the model parameters:
   Zeff = 1.0d0	! electron charge
   El_inelast = numpar%El_inelast   ! to chose the model of CS calculations below
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

   ! In case it is delta-CDF model, we have to make sure the energy is not below its applicability limit:
   if (El_inelast == 3) then
      !VAL0:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
      VAL0:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         max_E0 = maxval(Material%CDF_valence%E0(:))
      else VAL0
         max_E0 = maxval(Material%Elements(j)%CDF(k)%E0(:))
      endif VAL0
      call find_Wmax_equal_Wmin(0.0d0, 0.0d0, .true., Ee, Ip, max_E0, Eeq)   ! module "CS_integration_limits"
      ! Check if delta-functional CDF works here:
!       if (Ee < 1.1d0*Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
      if (Ee < numpar%CDF_Eeq_factor * Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
         El_inelast = 4
      endif
   endif ! (El_inelast == 3)
   
   E_left = Ip ! [eV] minimal transferred energy
   E_right = (Ip + Ee)*0.5d0    ! [eV] maximal transferred energy
   ! Use of the maximal plasmon energy as the upper limit of integration is not included yet !!!         
   
   ! Sample the cross section:
   call random_number(RN)
   CS_tot = get_integral_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, dispers, m_eff, E_right)   ! below
   CS_sampled = RN*CS_tot
   
   if ( (RN < 1.0d-10) .or. (CS_sampled < 1.0d-12) )then   ! no need to sample, it's just the lower limit
      E_cur = (E_left + E_right)*0.5d0
!       print*, 'get_inelastic_energy_transfer', Ee, CS_sampled, RN, CS_tot, E_cur
   else ! sample it 
      ! Start finding CS:
      E_cur = (E_left + E_right)*0.5d0
      CS_cur = get_integral_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, dispers, m_eff, E_cur)   ! below
      ! Search by bisection method:
      do while (ABS(CS_cur - CS_sampled) > eps*CS_sampled)
         if (CS_cur > CS_sampled) then
            E_right = E_cur
         else
            E_left = E_cur
         endif
         E_cur = (E_left + E_right)/2.0d0
         if (abs(E_left - E_right) < eps) exit  ! precise enough
         CS_cur = get_integral_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, dispers, m_eff, E_cur)   ! below
      enddo
   endif
   
   ! Output: sampled transferred energy:
   dE = E_cur
end function get_inelastic_energy_transfer_OLD




subroutine renormalize_alpha_CDF(numpar, Material, switch_SHI) ! module "CDF_get_from_data"
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   type(Target_atoms), dimension(:), intent(inout), target :: Material  !material parameters of each target that it's constructed of
   logical, intent(in), optional :: switch_SHI  ! if it is SHI
   !---------------------------------------
   real(8) :: CS_Ritchie, CS_Delta, Ee
   real(8) :: Se1, Zeff, Ip, Eeq, max_E0, hw_ph_max
   integer :: i, j, k, N_targets, N_elements, N_shells
   type(Atom_kind), pointer :: Element
   logical :: it_is_SHI

   ! Identify if this subroutine is used for SHI:
   if (present(switch_SHI)) then    ! if user provided
      it_is_SHI = switch_SHI
   else ! by default
      it_is_SHI = .false.
   endif
   Ip = 0.0d0 ! unused
   
   select case (numpar%El_inelast)  ! chose which model for electron inelastic cross section to use
      case (3:5)   ! in case Delta-CDF is used
      Zeff = 1.0d0  ! electron charge
      N_targets = size(Material)	! that's how many different targets user specified
      
      ! Go through all the shells of all elements of all targets:
      TRGT:do i = 1, N_targets	! for each target
         
         ! Inelastic CDF:
         N_elements = size(Material(i)%Elements)	! that's how many different elements are in this target
         
         LMNT:do j =1, N_elements	! for each element
            Element => Material(i)%Elements(j)	! all information about this element
            N_shells = Element%N_shl
            
            ! Calculate total cross sections for all shells of this element:
            ALLSHL:do k = 1, N_shells
               
               ! Define the switching energy for this shell:
               VAL0:if ( (Material(i)%Elements(j)%valent(k)) .and. (allocated(Material(i)%CDF_valence%A)) ) then    ! Valence band
                  max_E0 = maxval(Material(i)%CDF_valence%E0(:))
               else VAL0
                  max_E0 = maxval(Material(i)%Elements(j)%CDF(k)%E0(:))
               endif VAL0
               call find_Wmax_equal_Wmin(0.0d0, 0.0d0, .true., 100.0d0, Ip, max_E0, Eeq)   ! module "CS_integration_limits"
               ! The energy where cross section switches from Ritchie to Delta:
               Ee = numpar%CDF_Eeq_factor * Eeq
               
               ! Get Ritchie CDF cross secion:
               VAL3:if ( (Material(i)%Elements(j)%valent(k)) .and. (allocated(Material(i)%CDF_valence%A)) ) then    ! Valence band
                  Ip = Material(i)%DOS%Egap ! [eV] band gap
                  
                  select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
                  case default  ! free electron
                     call CDF_total_CS_nonrel(numpar, CS_Ritchie, Se1, Ee, g_me, Zeff, Ip, Material(i)%T_eV, Material(i)%CDF_valence, g_me, &
                          Material(i)%DOS%k, Material(i)%DOS%Eff_m, .true., 1.0d0, Material(i)%At_Dens, Material(i)%DOS, &
                          numpar%CDF_model)  ! below
                  case (1)  ! Plasmon pole
                     call CDF_total_CS_nonrel(numpar, CS_Ritchie, Se1, Ee, g_me, Zeff, Ip, Material(i)%T_eV, Material(i)%CDF_valence, &
                            g_me, Material(i)%DOS%k, Material(i)%DOS%Eff_m, .true., 1.0d0, Material(i)%At_Dens, Material(i)%DOS, &
                            numpar%CDF_model, CDF_dispers=numpar%CDF_dispers, v_f=Material(i)%DOS%v_f)  ! below
                  case (2)  ! Ritchie extended
                     call CDF_total_CS_nonrel(numpar, CS_Ritchie, Se1, Ee, g_me, Zeff, Ip, Material(i)%T_eV, Material(i)%CDF_valence, &
                            g_me, Material(i)%DOS%k, Material(i)%DOS%Eff_m, .true., 1.0d0, Material(i)%At_Dens, Material(i)%DOS, &
                            numpar%CDF_model, CDF_dispers=numpar%CDF_dispers)  ! below
                  end select ! (numpar%CDF_dispers) 
                  
                  ! Get Delta CDF cross secion:
                  CS_Delta = CDF_total_CS_delta(numpar%El_inelast, Ee, g_me, Zeff, Ip, Material(i)%At_Dens, Material(i)%CDF_valence, &
                                g_me, .true.) ! module "CDF_delta"
                 
                  ! Rescale alpha coefficients for this shell to match the cross-sections:
                  if (it_is_SHI) then
                     Material(i)%CDF_valence%alpha_SHI = Material(i)%CDF_valence%alpha_SHI * CS_Ritchie/CS_Delta
                  else ! for electrons:
                     Material(i)%CDF_valence%alpha = Material(i)%CDF_valence%alpha * CS_Ritchie/CS_Delta
                  endif
                  
               else VAL3 ! core shells
                  Ip = Material(i)%Elements(j)%Ip(k)    ! [eV] ionization potential
               
                  select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
                  case default  ! free electron
                     call CDF_total_CS_nonrel(numpar, CS_Ritchie, Se1, Ee, g_me, Zeff, Ip, Material(i)%T_eV, Material(i)%Elements(j)%CDF(k), &
                            g_me, Material(i)%DOS%k, Material(i)%DOS%Eff_m, .true., 1.0d0, Material(i)%At_Dens, Material(i)%DOS, &
                            numpar%CDF_model)  ! below
                  case (1)  ! Plasmon pole
                     call CDF_total_CS_nonrel(numpar, CS_Ritchie, Se1, Ee, g_me, Zeff, Ip, Material(i)%T_eV, Material(i)%Elements(j)%CDF(k), &
                            g_me, Material(i)%DOS%k, Material(i)%DOS%Eff_m, .true.,1.0d0, Material(i)%At_Dens, Material(i)%DOS, &
                            numpar%CDF_model, CDF_dispers=numpar%CDF_dispers, v_f=Material(i)%DOS%v_f)  ! below
                  case (2)  ! Ritchie extended
                     call CDF_total_CS_nonrel(numpar, CS_Ritchie, Se1, Ee, g_me, Zeff, Ip, Material(i)%T_eV, Material(i)%Elements(j)%CDF(k), &
                            g_me, Material(i)%DOS%k, Material(i)%DOS%Eff_m, .true., 1.0d0, Material(i)%At_Dens, Material(i)%DOS, &
                            numpar%CDF_model, CDF_dispers=numpar%CDF_dispers)  ! below
                  end select ! (numpar%CDF_dispers) 
                  
                  ! Get Delta CDF cross secion:
                  CS_Delta = CDF_total_CS_delta(numpar%El_inelast, Ee, g_me, Zeff, Ip, Material(i)%At_Dens, &
                                Material(i)%Elements(j)%CDF(k), g_me, .true.) ! module "CDF_delta"
                  
                  ! Rescale alpha coefficients for this shell to match the cross-sections:
                  if (it_is_SHI) then
                     Material(i)%Elements(j)%CDF(k)%alpha_SHI = Material(i)%Elements(j)%CDF(k)%alpha_SHI * CS_Ritchie/CS_Delta
                  else ! for electrons:
                     Material(i)%Elements(j)%CDF(k)%alpha = Material(i)%Elements(j)%CDF(k)%alpha * CS_Ritchie/CS_Delta
                  endif

                  ! SKIP ALPHA RENORMALIZATION, RENORMALIZE A(:) INSTEAD: (Testing)
                  ! Rescale A coefficients for this shell to match delta-cross-sections: (Testing)
!                   Material(i)%Elements(j)%CDF(k)%A = Material(i)%Elements(j)%CDF(k)%A * CS_Delta/CS_Ritchie (Testing)
               endif VAL3

            enddo ALLSHL
         enddo LMNT
      enddo TRGT 
      nullify(Element)
   end select ! Inelastic

   ! Elastic CDF:
   select case (numpar%El_elastic)  ! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   case (1,5) ! CDF
      N_targets = size(Material)	! that's how many different targets user specified
      ! Go through all the shells of all elements of all targets:
      TRGT2:do i = 1, N_targets	! for each target
         
         ! Maximal phonon frequency is defined by the maximal phononic CDF peak:
         !hw_ph_max = maxval( Material(i)%CDF_phonon%E0(:) + Material(i)%CDF_phonon%Gamma(:) )
         hw_ph_max = maxval( Material(i)%CDF_phonon%E0(:) + 10.0d0*Material(i)%CDF_phonon%Gamma(:) )  ! Testing
         ! In case it is delta-CDF model, we have to make sure the energy is not below its applicability limit:
         max_E0 = maxval(Material(i)%CDF_phonon%E0(:))
         call find_Wmax_equal_Wmin(g_me, Material(i)%Mean_Mass, .false., 10.0d0, 1.0d-6, max_E0, &
                        Eeq, hw_ph_max) ! module "CS_integration_limits"
            
         Ee = numpar%CDF_Eeq_elast*Eeq
            
         ! Target mean atomic number:
         if (numpar%CDF_elast_Zeff /= 1) then
            Zeff = 1.0d0 + Equilibrium_charge_SHI(Ee, g_me, Material(i)%Mean_Z, (Material(i)%Mean_Z-1.0d0), 0, 1.0d0) ! module "SHI_charge_state"
         else
            Zeff = 1.0d0    ! electron charge
         endif
            
!          print*, 'renormalize_alpha_CDF Elastic', Ee, Eeq
         ! Ritchie CDF cross section:
         call CDF_total_CS_nonrel(numpar, CS_Ritchie, Se1, Ee, g_me, Zeff, 1.0d-20, Material(i)%T_eV, Material(i)%CDF_phonon, &
                 Material(i)%Mean_Mass, &
                 Material(i)%DOS%k, Material(i)%DOS%Eff_m, .false., 1.0d0, Material(i)%At_Dens, Material(i)%DOS, numpar%CDF_model, &
                 hw_phonon = hw_ph_max)  ! module "CS_electrons_inelastic"     
                    
         ! Detla CDF cross section:
         CS_Delta = CDF_total_CS_delta(numpar%El_elastic, Ee, g_me, Zeff, 1.0d-16, Material(i)%At_Dens, Material(i)%CDF_phonon, &
                 Material(i)%Mean_Mass, .false., hw_phonon=hw_ph_max)    ! module "CDF_delta"
                                 
         ! Rescale alpha coefficients for phonons to match the cross-sections:
         Material(i)%CDF_phonon%alpha = Material(i)%CDF_phonon%alpha * CS_Ritchie/CS_Delta

      enddo TRGT2
   end select ! Elastic
   
end subroutine renormalize_alpha_CDF



function get_integral_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, dispers, m_eff, Emax) result (sigma)
   real(8) :: sigma	! [A^2] cross section
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   integer, intent(in) :: El_inelast    ! model for the inelastic cross section
   real(8), intent(in) :: Zeff  ! effective charge of the incident particle
   integer, intent(in) :: dispers, m_eff	! dispersion relation and effective mass
   real(8), intent(in) :: Emax  ! [eV] upper integration limit
   !---------------------------------------
   real(8) :: Se1
   integer :: ish

   ! Now chose the model of CS:
   select case (El_inelast)  ! chose which model for electron inelastic cross section to use
   case (1)	! relativistic CDF with Ritchi model
      !VAL:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
      VAL:if ( ( k==0 ) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         sigma = CDF_total_CS(Ee, g_me, Zeff, Ip, Material%T_eV, Material%At_Dens, Material%CDF_valence, g_me, dispers, &
                 m_eff, DOS%k, DOS%Eff_m, Material%DOS%v_f, Wmax_in=Emax)    ! below -- relativistic version, too slow, do not use!
      else VAL  ! core shells:
         sigma = CDF_total_CS(Ee, g_me, Zeff, Ip, Material%T_eV, Material%At_Dens, Material%Elements(j)%CDF(k), g_me, dispers, &
                 m_eff, DOS%k, DOS%Eff_m, Material%DOS%v_f, Wmax_in=Emax) ! below -- relativistic version, too slow, do not use!
      endif VAL
      
   case (2)	! RBEB
      if ( k==0 ) then  ! valence band
         ish = size(Material%Elements(1)%Ne_shell)  ! fix proper shell chosing later!
         sigma = RBEB_integral_CS(Ee, Material%Elements(1)%Ne_shell(ish), Ip, Material%Elements(1)%Ek(ish), Emax)  ! below
      else
         sigma = RBEB_integral_CS(Ee, Material%Elements(j)%Ne_shell(k), Ip, Material%Elements(j)%Ek(k), Emax)  ! below
      endif
      
   case (3,5)	! CDF with delta-functions
       !VAL2:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
       VAL2:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         sigma = CDF_total_CS_delta(El_inelast, Ee, g_me, Zeff, Ip, Material%At_Dens, Material%CDF_valence, g_me, .true., Emax_in = Emax)    ! module "CDF_delta"
      else VAL2  ! core shells:
         sigma = CDF_total_CS_delta(El_inelast, Ee, g_me, Zeff, Ip, Material%At_Dens, Material%Elements(j)%CDF(k), g_me, .true., Emax_in = Emax)    ! module "CDF_delta"
      endif VAL2
   
   case (4) ! nonrelativistic CDF
      ! Select model for slow particles:
      !VAL3:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
      VAL3:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, &
                 DOS%k, DOS%Eff_m, .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, Wmax_in=Emax)  ! below
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                  .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                  CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, Wmax_in=Emax)  ! below
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                  .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, CDF_dispers=numpar%CDF_dispers, Wmax_in=Emax)  ! below
         end select ! (numpar%CDF_dispers) 
      else VAL3 ! core shells
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, &
                 DOS%k, DOS%Eff_m, .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, Wmax_in=Emax)  ! below
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                  .true.,1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                  CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, Wmax_in=Emax)  ! below
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                  .true., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, CDF_dispers=numpar%CDF_dispers, Wmax_in=Emax)  ! below
         end select ! (numpar%CDF_dispers) 
      endif VAL3
      
   case default ! exclude
      sigma = 0.0d0
   end select
end function get_integral_inelastic_CS


!-----------------------------------------------------------------------------
! Non-relativistic CDF cross section, Ritchie from [6], or Mermin:
subroutine CDF_total_CS_nonrel(numpar, sigma, Se, Ekin, Mass, Zeff, Ip, T_target, CDF, M_sc, &
                   k, Eff_m, identical, m_eff, At_Dens, DOS, CDF_model, &
                   v_f, Zmodel, hw_phonon, Nel, CDF_dispers, Wmax_in, Sigma_sampled, E_sampled)
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), intent(out) :: sigma	! [A^2] cross section
   real(8), intent(out) :: Se  ! [eV/A] energy loss
   real(8), intent(in) :: Ekin	! [eV] kinetic energy of incident electron
   real(8), intent(in) :: Mass	! [kg] mass of the incident particle
   real(8), intent(in) :: Zeff	! [e] effective charge of the incident particle
   real(8), intent(in) :: Ip	! [eV] ionization potential of the concidered shell
   real(8), intent(in) :: T_target	! [eV] target temeprature
   type(Ritchi_CDF), intent(in) :: CDF	! CDF coefficients
   real(8), intent(in) :: M_sc	! [kg] mass of the target scattering center (eletron, atom)
   real(8), dimension(:), intent(in) :: k, Eff_m	! arrays with the target dispersion curve from DOS
   logical, intent(in) :: identical ! are those particles identical or not (electron-electron vs hole-electron, etc.)
   real(8), intent(in) :: m_eff ! coefficient of the effective mass of the incident particle
   real(8), intent(in):: At_Dens ! atomic density [1/cm^3]
   type(Density_of_states), intent(in) :: DOS    ! DOS
   integer, intent(in) :: CDF_model ! which model to use for slow particles
   real(8), intent(in), optional :: v_f	! fermi velosity
   integer, intent(in), optional :: Zmodel	! sets model choice for the effective charge for SHI (variable not used for other particles)
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency, for particle-phonon scattering
   real(8), intent(in), optional :: Nel ! electrons per atom (to get plasmon energy)
   integer, intent(in), optional :: CDF_dispers     ! whilch model to use for the extension of the dielectric functions
   real(8), intent(in), optional :: Wmax_in ! [eV] maximal transfered energy
   real(8), intent(in), optional :: Sigma_sampled   ! to search for transfered energy, sampled cross section
   real(8), intent(inout), optional :: E_sampled    ! energy corresponding to sampled cross section
   !----------------------------------------
   real(8) :: eps, Wmax, Wmin, M_in, Mc2, mtc2, mass_temp, Eplasmon, dE, E_cur, E_cur0, E_cur1, ddEdx
   real(8) :: sigma1, sigma0, sigma2, sigma_cur, dEdx1, dEdx0, dEdx2, dEdx_cur
   real(8) :: dEdx, prefact, Tfact0, Tfact1, Tfact2, dsigma_max, dE_test
   integer :: n
   
!    if (present(hw_phonon)) then
!       print*, 'CDF_total_CS_nonrel', hw_phonon, Ekin, Ip
!    endif
   
   if (Ekin > Ip) then  ! ionization is possible
      eps = 1.0d-6 ! margin within which the effective mass is equal to zero
      ! Get the particle mass:
      if (m_eff > eps) then ! a constant coefficient * me
         M_in = m_eff * g_me
      elseif (abs(m_eff) < eps) then ! equals to user-provided mass (electron, ion, etc.)
         M_in = Mass   
      else  ! effective mass from DOS
         call interpolate_data_single(DOS%E, DOS%Eff_m, (DOS%E_VB_top - Ekin), mass_temp) ! module "Little_subroutines"
         M_in = mass_temp * g_me
         !print*, 'M_eff:', Ekin, mass_temp, DOS%E_VB_top - Ekin
      endif
      ! Incident particle mass in the energy units:
      Mc2 = rest_energy(M_in)   ! module "Relativity"
      ! Target particle mass in the energy units:
      mtc2 = rest_energy(M_sc)   ! module "Relativity"
      ! Upper integration limit [eV]:
      if (present(Wmax_in)) then
         Wmax = Wmax_in ! user defined
      elseif (present(hw_phonon)) then
         !Wmax = W_max(Mc2, mtc2, identical, Ekin, Ip, hw_phonon) ! module "CS_integration_limits"
         Wmax = W_max_nonrel_free(Mc2, mtc2, identical, Ekin, Ip, hw_phonon)   ! module "CS_integration_limits"
      else
         !Wmax = W_max(Mc2, mtc2, identical, Ekin, Ip) ! module "CS_integration_limits"
         Wmax = W_max_nonrel_free(Mc2, mtc2, identical, Ekin, Ip)   ! module "CS_integration_limits"
      endif
      ! Minimal transfered energy:
      Wmin = W_min(Ip, Mc2, mtc2, Ekin, Ip) ! module "CS_integration_limits"
!       print*, 'CDF_total_CS_2', Ekin, Wmax, Wmin

      WiP:if ( Wmax > Ip ) then  ! ionization is possible
!       WiP:if ( Wmax > max(Ip,Wmin) ) then  ! ionization is possible
         ! Get the constant prefactor for the cross section
         prefact = (Zeff*Zeff) * P_prefactor_nonrel(M_in, Ekin, At_Dens)  ! below
         
         ! Lower integration limit:
         Wmin = Ip ! [eV] minimal transferred energy
   
         ! Use maximal plasmon energy as the upper limit of integration
         if ( (present(Nel)) .and.  (.not.present(Wmax_in)) ) then       ! If included, only for total cross section
            Eplasmon = plasmon_energy(Nel, At_Dens, Ip)    ! module "CS_general_tools"
!             if (Eplasmon >= Wmax) Wmax = Eplasmon ! single atom vs plasmon
         endif
         if (Ekin < Wmax) Wmax = Ekin ! no more than the total electron energy
   
         ! Start integration by energy (non-relativistic case):
         if (present(hw_phonon)) then
            n = numpar%CDF_int_n_elastE    ! number of integration steps
         else
            n = numpar%CDF_int_n_inelastE    ! number of integration steps
         endif

         E_cur = Wmin    ! to start integration
         ! Just to start from the min transfered energy:
         if (present(Sigma_sampled) .and. present(E_sampled)) then   ! it is a searching expedition
            E_sampled = E_cur
         endif
         
         ! Differential cross section:
         if (present(v_f) .and. present(CDF_dispers)) then
            call Diff_cross_section(numpar,sigma0, Ekin, E_cur, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model, &
                        E0_model=CDF_dispers, v_f=v_f, Mass=M_sc)  ! below
         elseif (present(CDF_dispers)) then
            call Diff_cross_section(numpar,sigma0, Ekin, E_cur, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model, E0_model=CDF_dispers)  ! below
         elseif (present(hw_phonon)) then
            call Diff_cross_section(numpar,sigma0, Ekin, E_cur, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model, hw_phonon=hw_phonon)  ! below
         else ! by defualt, no additional parameters
            call Diff_cross_section(numpar,sigma0, Ekin, E_cur, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model)  ! below
         endif
         ! Temperature factor:
         Tfact0 = temperature_factor(E_cur, T_target)    ! module "CS_general_tools"
         sigma0 = sigma0*Tfact0   ! include temprature factor
         dsigma_max = sigma0    ! current definition of the maximal value, to be changed below
         dEdx = 0.0d0
         sigma = 0.0d0
         
!          print*, 'Before', E_cur, dE, Wmin, Wmax
         
         ! Integrate:
         SRCH:do while (E_cur <= Wmax)
            ! Change integration step as needed by change of the energy:
            dE = remap_energy_step(E_cur, n, G=CDF%Gamma, E0=CDF%E0)  ! module "CS_general_tools" CORRECT
            
!             if (.not.present(hw_phonon)) then ! do not do for phonons
!                dE = remap_energy_step(E_cur, n, minval(CDF%Gamma))  ! module "CS_general_tools"               
!             else
!                dE = remap_energy_step(E_cur, n)  ! module "CS_general_tools" CORRECT
!             endif
            
            ! Simpson integration:
            E_cur0 =  E_cur + dE/2.0d0    ! middle step [eV]
            
!             if (Ekin == 1000.0d0) then
!                E_cur0 = 0.1d0 ! Testing
!                write(*,'(a,es,es,es,es,es,es)') 'CS:', Ekin, dE, E_cur0, Wmin, Wmax ! Testing
!             endif
            
            
            if (E_cur0 > Wmax) exit
            ! Differential cross section at this point:
            if (present(v_f) .and. present(CDF_dispers)) then
               call Diff_cross_section(numpar,sigma1, Ekin, E_cur0, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model, &
                            E0_model=CDF_dispers, v_f=v_f, Mass=M_sc)  ! below
            elseif (present(CDF_dispers)) then
               call Diff_cross_section(numpar,sigma1, Ekin, E_cur0, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model, E0_model=CDF_dispers)  ! below
            elseif (present(hw_phonon)) then
               call Diff_cross_section(numpar,sigma1, Ekin, E_cur0, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model, hw_phonon=hw_phonon)  ! below
            else ! by defualt, no additional parameters
               call Diff_cross_section(numpar,sigma1, Ekin, E_cur0, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model)  ! below
            endif
            ! Temperature factor:
            Tfact1 = temperature_factor(E_cur0, T_target)    ! module "CS_general_tools"
            sigma1 = sigma1*Tfact1   ! include temprature factor
      
!             write(*,'(a,es,es,es,es,es,es)') 'CS:', Ekin, dE, E_cur0, Wmin, Wmax, sigma1 ! Testing
            
            ! End-point for this step:
            E_cur1 = E_cur + dE
            if (E_cur1 > Wmax) exit
            ! Differential cross section at this point:
            if (present(v_f) .and. present(CDF_dispers)) then
               call Diff_cross_section(numpar,sigma2, Ekin, E_cur1, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model, &
                            E0_model=CDF_dispers, v_f=v_f, Mass=M_sc)  ! below
            elseif (present(CDF_dispers)) then
               call Diff_cross_section(numpar,sigma2, Ekin, E_cur1, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model, E0_model=CDF_dispers)  ! below
            elseif (present(hw_phonon)) then
               call Diff_cross_section(numpar,sigma2, Ekin, E_cur1, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model, hw_phonon=hw_phonon)  ! below
            else ! by defualt, no additional parameters
               call Diff_cross_section(numpar,sigma2, Ekin, E_cur1, CDF%A, CDF%E0, CDF%Gamma, M_in, M_sc, CDF_model)  ! below
            endif
            ! Temperature factor:
            Tfact2 = temperature_factor(E_cur1, T_target)    ! module "CS_general_tools"
            sigma2 = sigma2*Tfact2    ! include temprature factor
         
            ! Integral of diff.cross section according to Simpson method:
            sigma_cur = dE/6.0d0*(sigma0 + 4.0d0*sigma1 + sigma2)
            !if (sigma_cur > dsigma_max) dsigma_max = sigma_cur  ! redefine the maximal value
            ! Sum it into the total cross section:
            sigma = sigma + sigma_cur     ! integrated cross section ([eV], to be converted into [A^2] below)
            ! Sum the energy loss (stopping power) function:
            dEdx = dEdx + E_cur*sigma_cur ! stopping power  ([eV^2], to be converted into [eV/A] below)
            sigma0 = sigma2   ! save for the next step 
            E_cur = E_cur1  ! [eV] save for the next step
            
            if (present(Sigma_sampled) .and. present(E_sampled)) then   ! it is a searching expedition
               E_sampled = E_cur     ! save corresponding energy
               if (sigma*prefact >= Sigma_sampled) then ! found the cross section that equals to the sampled one
!                   print*, 'CDF_total_CS_nonrel', sigma, Sigma_sampled, E_cur 
                  exit SRCH ! if found the energy, no need to continue, exit the cicle
               endif
            endif
            
            !if (sigma_cur/dsigma_max < 1.0d-3) exit   ! no need to continue integration, the function is already too small
!             if (abs(Ekin-1.0d0) < 1.0d-2) print*, 'CDF_total_CS_3', Ekin, E_cur, 1.0d0/( (At_Dens*1.0d-24) * sigma*(Zeff*Zeff)*P_prefactor_nonrel(M_in, Ekin, At_Dens) ),  (At_Dens*1.0d-24) 
         enddo SRCH
!          print*, 'CDF_total_CS_4', prefact, sigma
         sigma = prefact * sigma  ! [A^2]
!          print*, 'CDF_total_CS_5', sigma, E_cur
         Se = prefact * (At_Dens*1.0d-24) * dEdx    ! [eV/A]
      else WiP  ! ionization is impossible
         sigma = 0.0d0
         Se = 0.0d0
         if (present(Sigma_sampled) .and. present(E_sampled)) then   ! it is a searching expedition
            E_sampled = Ip
         endif
      endif WiP
      
!       write(*,'(a,f,f,f,es)') 'CS: ', Ekin, Wmax, Wmin, sigma
!       print*, 'M_in', M_in, Ekin
!       !call print_time_step('Time:', Ekin/1.0d6, msec=.true.)   ! module "Little_subroutines"
!       pause 'CDF_total_CS_nonrel'
   
   else ! ionization is impossible
      sigma = 0.0d0
      Se = 0.0d0
      if (present(Sigma_sampled) .and. present(E_sampled)) then   ! it is a searching expedition
         E_sampled = Ip
      endif
   endif
   
!    if (Ekin == 1000.0d0) pause 'CDF_total_CS_nonrel' ! Testing
   
end subroutine CDF_total_CS_nonrel



subroutine Diff_cross_section(numpar, dCS, Ekin, hw, A, E, Gamma, M_in, M_t, CDF_model, v_f, E0_model, Mass, hw_phonon)
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), intent(out) :: dCS  ! single-differential cross section (unnormalized, hence in [eV])
   real(8), intent(in) :: Ekin  ! kinetic energy of incident particle [eV]
   real(8), intent(in) :: hw   ! transferred energy [eV]
   real(8), dimension(:), intent(in) :: A, E, Gamma     ! CDF parameters
   real(8), intent(in) ::  M_in, M_t  ! [kg] masses of the incident particle and the target one
   integer, intent(in) :: CDF_model ! Which model for CDF: Ritchie, Mermin, etc.
   real(8), intent(in), optional ::  v_f    ! target fermi velosity [sqrt(eV/kg)]
   integer, intent(in), optional :: E0_model     ! whilch model to use for the extension of the dielectric functions
   real(8), intent(in), optional :: Mass           ! mass [kg]
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   !------------------------------------------
   integer :: i, n
   real(8) :: Qcenter, Qmin, Qmax, dsigma0, dsigma1, dsigma2, dsigma_cur
   real(8) :: hq, hq0, hq1, dq, eps, dq_min, q_max_test, Mr

   eps = 1.0d-6 ! margin within which masses are considered identical
   if (abs((M_in-M_t)/M_t) > eps) then   ! different masses of the target and incident particles
      ! Transferred energy corresponding to the minimum transfered momentum:
      Qmin = Q_min_nonrel(Ekin, hw, M_in, M_t)   ! module "CS_integration_limits"
      ! Transferred energy corresponding to the maximum transfered momentum:
      Qmax = Q_max_nonrel(Ekin, hw, M_in, M_t)   ! module "CS_integration_limits"
      ! Define number used for integration grid:
      n = numpar%CDF_int_n_elastQ    ! number of integration steps
   else ! electron by default
      ! Transferred energy corresponding to the minimum transfered momentum:
      Qmin = Q_min_nonrel(Ekin, hw)   ! module "CS_integration_limits"
      ! Transferred energy corresponding to the maximum transfered momentum:
      Qmax = Q_max_nonrel(Ekin, hw)   ! module "CS_integration_limits"
      ! Define number used for integration grid:
      n = numpar%CDF_int_n_inelastQ    ! number of integration steps
   endif

   ! Mass ratio:
!    Mr = g_me / M_t

   ! Considering that a CDF peak has a fixed width within Drude model,
   ! we don't have to integrate a lot of empty space too far from it,
   ! the peak of the oscillator along Q is (hw-E0):
   if (abs((M_in-M_t)/M_t) > eps) then   ! different masses of the target and incident particles
!       Qmin = max(1.0d-8, Qmin, minval(hw-E(:)-5.0d0*Gamma(:))) ! we have to also exclude zero
      Qmin = max(1.0d-10, Qmin, minval(hw-E(:)-10.0d0*Gamma(:))) ! we have to also exclude zero
   else     ! electron-electron scattering
      Qmin = max(1.0d-5, Qmin, minval(hw-E(:)-5.0d0*Gamma(:))) ! we have to also exclude zero
   endif
   q_max_test = maxval(hw-E(:)+10.0d0*Gamma(:))
   if (q_max_test > Qmin) then  ! in case of phonons, it may turn unphisycal, so exclude such cases
      Qmax = min(Qmax, q_max_test)
   endif
   
   dsigma_cur = 0.0d0 ! starting integration, diff.CS [1/eV]
   hq = Qmin    ! to start with
   ! Get the dielectric function:
   select case (CDF_model)
   case default ! Ritchie
      if (present(v_f) .and. present(Mass) .and. present(E0_model)) then
         dsigma0 = 1.0d0/hq * Diel_func(A, E, Gamma, hw, hq, E0_model=E0_model, v_f=v_f, Mass=Mass) ! module "CDF_Ritchi"
      elseif (present(E0_model)) then
         dsigma0 = 1.0d0/hq * Diel_func(A, E, Gamma, hw, hq, E0_model=E0_model) ! module "CDF_Ritchi"
      else ! be defualt, no additional parameters
         dsigma0 = 1.0d0/hq * Diel_func(A, E, Gamma, hw, hq) ! module "CDF_Ritchi"
      endif
   case (1) ! Mermin
      dsigma0 = 1.0d0/hq * Mermin_Diel_func(A, E, Gamma, hw, hq) ! module "CDF_Mermin"
   endselect
   
   ! integrate double-diff.CS numerically:
   do while (hq < Qmax) ! no matter how many points, go till the end
      ! Make the step small around the peak, but larger when far from it:
      dq = remap_dq(hq, n, hw, E, Gamma) ! module "CS_general_tools"
      
      ! Simpson integration:
      hq0 = hq + dq/2.0d0
      if (hq0 > Qmax) hq0=Qmax
      select case (CDF_model)
      case default ! Ritchie
         if (present(v_f) .and. present(Mass) .and. present(E0_model)) then
            dsigma1 = 1.0d0/hq0 * Diel_func(A, E, Gamma, hw, hq0 , E0_model=E0_model, v_f=v_f, Mass=Mass) ! module "CDF_Ritchi"
         elseif (present(E0_model)) then
            dsigma1 = 1.0d0/hq0 * Diel_func(A, E, Gamma, hw, hq0 , E0_model=E0_model) ! module "CDF_Ritchi"
         else ! be defualt, no additional parameters
            dsigma1 = 1.0d0/hq0 * Diel_func(A, E, Gamma, hw, hq0 ) ! module "CDF_Ritchi"
         endif
      case (1) ! Mermin
         dsigma1 = 1.0d0/hq0 * Mermin_Diel_func(A, E, Gamma, hw, hq0 ) ! module "CDF_Mermin"
      endselect

!       write(*,'(a,es,es,es,es,es,es,es)') 'DCS:', Ekin, hw, Qmin, Qmax, dq, hq, dsigma1    ! Testing
!       if ( abs(Ekin - 1000.0d0)/1000.0d0 < 0.01d0) write(*,'(a,es,es,es,es,es,es,es)') 'DCS:', Ekin, hw, Qmin, Qmax, dq, hq, dsigma1    ! Testing
!       if (Ekin > 7500.0) write(*,'(a,es,es,es,es,es,es,es)') 'Diff_cross_section:', Ekin, Qmin, Qmax, dq, hq, dsigma1
      
      hq1 = hq + dq
      if (hq1 > Qmax) hq1=Qmax
      select case (CDF_model)
      case default ! Ritchie
         if (present(v_f) .and. present(Mass) .and. present(E0_model)) then
            dsigma2 = 1.0d0/hq1 * Diel_func(A, E, Gamma, hw, hq1 , E0_model=E0_model, v_f=v_f, Mass=Mass) ! module "CDF_Ritchi"
         elseif (present(E0_model)) then
            dsigma2 = 1.0d0/hq1 * Diel_func(A, E, Gamma, hw, hq1 , E0_model=E0_model) ! module "CDF_Ritchi"
         else ! be defualt, no additional parameters
            dsigma2 = 1.0d0/hq1 * Diel_func(A, E, Gamma, hw, hq1 ) ! module "CDF_Ritchi"
         endif
      case (1) ! Mermin
         dsigma2 = 1.0d0/hq1 * Mermin_Diel_func(A, E, Gamma, hw, hq1 ) ! module "CDF_Mermin"
      endselect
      
      ! Assemble the integrated diff.CS:
      dsigma_cur = dsigma_cur + dq/6.0d0*(dsigma0 + 4.0d0*dsigma1 + dsigma2)
      
!       if ((abs(Ekin - 1.0) < 1.0d-2) .and. (abs(hw - 0.1d0) < 1.0d-2)) then
!           write(*,'(a,f,f,f,f,f)') 'Diff_cross_section', Ekin, hw, hq, dq, dsigma_cur
!       endif
!       
      dsigma0 = dsigma2 ! save for the next step
      hq0 = hq1 ! save for the next step
      hq = hq0
   enddo

   ! Output variable:
   ! Mass ratio comes from full relativistic DSF connection to CDF: (Q+M_t*c^2)/m_e*c^2
   dCS = dsigma_cur      ! single-diff.cross section; unnormalized, in [eV]
   
!    if ( abs(Ekin - 1000.0d0)/1000.0d0 < 0.01d0) pause 'Diff_cross_section' ! Testing
end subroutine Diff_cross_section



pure function P_prefactor_nonrel(M, E, nat) result(P)
   real(8) :: P
   real(8), intent(in) :: M, E  ! mass [kg] and kinetic energy [eV] of the projectile
   real(8), intent(in) :: nat   ! [1/cm^3] atomic density
   ! Factor of 0.5 comes from the fact that dQ ~ 1/2 q*dq
   P = 0.5d0*1.0d24/(g_Pi*g_a0*nat)*M/(g_me*E) ! converting into [A^2/eV]
end function P_prefactor_nonrel



!-----------------------------------------------------------------------------
! Relativistic CDF cross section from [1], generalization of Eq.(3.148) to finite temperatures of the target, according to [6]:
 function CDF_total_CS(Ekin, Mass, Zeff, Ip, T_target, n_target, CDF, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, v_f, Zmodel, Wmax_in) result(sigma)
   real(8) sigma	! [A^2] cross section
   real(8), intent(in) :: Ekin	! [eV] kinetic energy of incident electron
   real(8), intent(in) :: Mass	! [kg] mass of the incident particle
   real(8), intent(in) :: Zeff	! [e] effective charge of the incident particle
   real(8), intent(in) :: Ip	! [eV] ionization potential of the conciderred shell
   real(8), intent(in) :: T_target	! [eV] target temeprature
   real(8), intent(in) :: n_target	! [1/cm^3] atomic/molecular density of the target atoms (depending on normalization of the CDF coefficients)
   type(Ritchi_CDF), intent(in) :: CDF	! CDF coefficients
   real(8), intent(in) :: M_sc	! [kg] mass of the scattering center (eletron, SHI, hole)
   integer, intent(in) :: CDF_dispers, CDF_m_eff	! indices of the models for target dispersion relation, and effective mass
   real(8), dimension(:), intent(in) :: k, Eff_m	! arrays with the target dispersion curve from DOS
   real(8), intent(in) :: v_f	! fermi velosity
   integer, intent(in), optional :: Zmodel	! sets model choice for the effective charge for SHI (variable not used for other particles)
   real(8), intent(in), optional :: Wmax_in ! [eV] maximal transferred energy (used for energy transfer)
   !----------------------------------------
   real(8) :: Wmax
   if (present(Wmax_in)) then   ! user provided
      Wmax = Wmax_in
   else ! maximal possible transfered energy for the total CS:
      Wmax = Ekin	! [eV] upper integration limit
   endif
   
   if (present(Zmodel)) then
      sigma = CDF_integral_CS(Ekin, Mass, Zeff, Ip, T_target, n_target, CDF, Wmax, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, v_f, Zmodel)	! see below
   else
      sigma = CDF_integral_CS(Ekin, Mass, Zeff, Ip, T_target, n_target, CDF, Wmax, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, v_f)	! see below
   endif
end function CDF_total_CS


 function CDF_integral_CS(Ekin, Mass, Zeff, Ip, T_target, n_target, CDF, Wmax, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, v_f, Zmodel) result(sigma)
   real(8) sigma	! [A^2] cross section
   real(8), intent(in) :: Ekin	! [eV] kinetic energy of incident electron
   real(8), intent(in) :: Mass	! [kg] mass of the incident particle
   real(8), intent(in) :: Zeff	! [e] effective charge of the incident particle
   real(8), intent(in) :: Ip	! [eV] ionization potential of the conciderred shell
   real(8), intent(in) :: T_target	! [eV] target temeprature
   real(8), intent(in) :: n_target	! [1/cm^3] atomic/molecular density of the target atoms (depending on normalization of the CDF coefficients)
   type(Ritchi_CDF), intent(in) :: CDF	! CDF coefficients
   real(8), intent(in) :: Wmax		! [eV] upper inegration limit
   real(8), intent(in) :: M_sc	! [kg] mass of the scattering center (eletron, SHI, hole)
   integer, intent(in) :: CDF_dispers, CDF_m_eff	! indices of the models for target dispersion relation, and effective mass
   real(8), dimension(:), intent(in) :: k, Eff_m	! arrays with the target dispersion curve from DOS
   real(8), intent(in) :: v_f	! fermi velosity
   integer, intent(in), optional :: Zmodel	! sets model choice for the effective charge for SHI (variable not used for other particles)
   !----------------------------------------
   real(8) :: Mc2, beta, v, W_cur, W_cur0, Wmin, W_mid, prefac, sigma_cur, eps, dW, dW_min, dW_max, dW_half
   real(8) :: expfac, dS, dSigma, dSigma_mid, dSigma0
   integer :: i, Ngrid
 
    if (Ekin > Ip) then	! ionization of this shell is possible
      eps = 1.0d-12	! precision limit
      dS = 0.01d0	! maximal allowed change in dSigma per step [%]
      Mc2 = 2.0d0*rest_energy(Mass)	! [eV] module "Relativity"   
      v = velosity_from_kinetic_energy(Ekin, Mass,afs=.false.)	! module "Relativity"
      beta = beta_factor(v)	! module "Relativity"
      prefac = 2.0d24/(g_Pi*g_a0*n_target*Mc2*beta)	! converting CS into [A^2]
   
      ! to start integration:
      sigma_cur = 0.0d0
      Ngrid = 100	! grid point for integration over W_cur
      Wmin = Ip	! [eV] lower integration limit
      W_cur = Wmin
!       dW_max = (Wmax - Wmin)/dble(Ngrid)	! maximal allowed integration step
      dW_max = (Wmax - Wmin)/2.0d0		! maximal allowed integration step
      dW_min = (1.0d0/(W_cur+1.0d0) + W_cur)/dble(Ngrid)
      dW = dW_min		! start with it, and later reduce if needed
      dW_half = dW*0.50d0
      dSigma = 0.0d0
      if (present(Zmodel)) then
         dSigma0 = dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W_cur, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_f, Zmodel) 	! below
      else
         dSigma0 = dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W_cur, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_f) 	! below
      endif
      
      i = 0
      WINT:do while (W_cur < Wmax)
         dW = (1.0d0/(W_cur+1.0d0) + W_cur*W_cur)/dble(Ngrid)
         i = i + 1	! steps counter
         W_cur0 = W_cur
         W_cur = W_cur + dW
         if (W_cur > Wmax) then	! if by chance we exceeded the limit
            W_cur = Wmax
            dW = Wmax - W_cur0
            dW_half = 0.5d0*dW
         endif
         
         if (present(Zmodel)) then
            dSigma  = dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W_cur, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_f, Zmodel) 	! below
         else
            dSigma  = dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W_cur, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_f) 	! below
         endif
      
         ! Adaptive step: if it's too large, reduce it:
         if ((dSigma > eps) .and. (dSigma0 > eps)) then	! makes sense to take care of integration
            ! Check if it's too large change per step:
!             print*, 'Sigma', dSigma, dSigma0
            do while ((ABS(dSigma - dSigma0) > dS*dSigma0) .and. (dW > dW_min))
               dW = dW*0.5d0	! reduce timestep
               dW_half = 0.5d0*dW	! adjust half-step correspondingly
               W_cur = W_cur0 + dW
               if (W_cur > Wmax) exit
               
               if (present(Zmodel)) then
                  dSigma  = dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W_cur, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_f, Zmodel) 	! below
               else
                  dSigma  = dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W_cur, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_f) 	! below
               endif
               
               if (dW < dW_min) exit	! too smal step
            enddo
            ! Check if it's too small change per step (inefficient):
            !do while ((ABS(dSigma - dSigma0) < 0.5d0*dS*dSigma0) .and. (dW < dW_max))
            do while (ABS(dSigma - dSigma0) < 0.5d0*dS*dSigma0)
               dW = 1.5d0*dW		! increase timestep
               dW_half = 0.5d0*dW	! adjust half-step correspondingly
               W_cur = W_cur0 + dW
               if (W_cur > Wmax) exit
               
               if (present(Zmodel)) then
                  dSigma  = dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W_cur, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_f, Zmodel) 	! below
               else
                  dSigma  = dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W_cur, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_f) 	! below
               endif
               
               if (dW > dW_max) exit	! too large step
            enddo
         else	! reset the step
            dW = dW_max
            dW_half = dW*0.5d0
            W_cur = W_cur0 + dW
         endif
         ! Now proceed with the integration:
         if (W_cur > Wmax) then	! if by chance we exceeded the limit
            W_cur = Wmax
            dW = Wmax - W_cur0
            dW_half = dW*0.50d0
         endif
         W_mid =  W_cur0 + dW_half
         if (present(Zmodel)) then
            dSigma_mid = dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W_mid, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_f, Zmodel) 	! below
         else
            dSigma_mid = dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W_mid, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_f) 	! below
         endif
      
         expfac = 1.0d0 - exp(-W_cur/T_target)	! asymetry factor from DSF and fluctuation-dissipation theorem [6]
         
         ! Add up the contributions according to Simpson-3/8 scheme:
         sigma_cur = sigma_cur + expfac*dW/6.0d0*(dSigma0 + 4.0d0*dSigma_mid + dSigma)

         !write(*,'(a,es,es,es,es,es)')  'S', W_cur, dW, ABS(dSigma - dSigma0), dS*dSigma0
         
         ! Save the data for the next point of integration:
         dSigma0 = dSigma
      enddo WINT
      
      sigma = prefac*sigma_cur	! [A^2] cross section
   else	! below threshold, no ionization
      sigma = 0.0d0
   endif
!    print*, 'T', Ip, Ekin, sigma
!    PAUSE 'CDF_integral_CS'
end function CDF_integral_CS



 function dCDF_CS_dQ(CDF, Ekin, Mass, Mc2, beta, v, W, M_sc, CDF_dispers, CDF_m_eff, k, Eff_m, Zeff, v_fermi, Zmodel) result(dSigma)
   real(8) dSigma
   type(Ritchi_CDF), intent(in) :: CDF	! CDF coefficients
   real(8), intent(in) :: Ekin, v, beta	! [eV] kinetic energy of incident particle, its velosity [m/s], and relativistiv beta factor
   real(8), intent(in) :: Mass, Mc2	! mass [kg] and rest energy [eV] of the incident particle
   real(8), intent(in) :: W	! [eV] inegration energy as independent variable
   real(8), intent(in) :: M_sc	! [kg] mass of the scattering center (eletron, SHI, hole)
   integer, intent(in) :: CDF_dispers, CDF_m_eff	! indices of the models for target dispersion relation, and effective mass
   real(8), dimension(:), intent(in) :: k, Eff_m	! arrays with the target dispersion curve from DOS
   real(8), intent(in) :: Zeff	! [e] effective charge of the incident SHI
   real(8), intent(in) :: v_fermi	! fermi velosity
   integer, intent(in), optional :: Zmodel	! sets model choice for the effective charge for SHI (not used for other particles)
   !---------------------------
   real(8) :: Z, Q, Q0, Q_mid, Qmin, Qmax, Qmc2, dQ, dQ_half, dbleN, dQmin, dQmax, dq_small
   real(8) :: EW, Emc2, sqr1, sqr2, Mc4, sqr12, Mc2_sc, Mc2_sc2, AB, ImE, dSigma0, dSigma_cur, dSigma_sum, dSigma_mid, Mass1
   integer :: i, j, Ngrid

   ! Rest energy of the scattering senter:
   Mc2_sc = rest_energy(M_sc)	! [eV] targets M*c^2, assuming the target's effective mass does not depend on momentum (INCONSISTENT WITH MASS-FROM-DOS MODEL!)
   Mc2_sc2 = 2.0d0*Mc2_sc
   Emc2 = Ekin + Mc2_sc2
   ! Get the integration limits, [1] P. 280, Eq.(A.31):
   sqr1 = sqrt(Ekin*Emc2)
   EW = Ekin - W
   sqr2 = sqrt(EW*(EW + Mc2_sc2))
   Mc4 = Mc2*Mc2	! square of the rest energy of the incident particle [eV]
   sqr12 = sqr1 - sqr2
   Qmin = sqrt(sqr12*sqr12 + Mc4) - Mc2	! [eV]
   sqr12 = sqr1 + sqr2
   Qmax = sqrt(sqr12*sqr12 + Mc4) - Mc2	! [eV]
   
!     print*, 'Qlims', Qmin, Qmax, W
   
   ! to start integration:
   Ngrid = 100	! chose minimal number of points for integration
   dbleN = dble(Ngrid)
   Q = Qmin
   Q0 = Q
   dQmin = minval(CDF%A(:))/dbleN
   dQmax = (Qmax - Qmin)/dbleN
!    print*, 'dQlims', dQmin, dQmax
   dQ = dQmin
   dSigma_sum = 0.0d0
   dSigma_cur = 0.0d0
   dSigma0 = 0.0d0
   dSigma_mid = 0.0d0
   do while (Q < Qmax)
      dQ_half = dQ*0.5d0
      Q = Q0 + dQ_half	! half step forward
      AB = get_AB(Q, Mc2, W, beta, Emc2)	! below
      ! Effective charge:
      if (present(Zmodel)) then	! set effective charge of the incident SHI
         Z = Zeff	! Equilibrium_charge_SHI(Ekin, Mass, Zeff, Ztarget, perc, Kind_Zeff, fixed_Zeff)	! module "SHI_charge_state"
      else
         Z = Zeff
      endif
      ! Get momentum transfer from the recoil energy:
      dq_small = momentum_from_kinetic_energy(Q, Mass)	! module "Relativity"
      ! Get the effective mass of the scattering center:
      Mass1 = chose_mass(CDF_m_eff, M_sc, dq_small, k, Eff_m) 	! below
      select case (CDF_dispers)	! model for the dispertion relation of the target:
      case default	! free electron
         ImE = Diel_func(CDF%A(:), CDF%E0(:), CDF%Gamma(:), W, dq_small, Mass=Mass1) 	! module "CDF_Ritchi"
      case (1)	! plasmon pole
         ImE = Diel_func(CDF%A(:), CDF%E0(:), CDF%Gamma(:), W, dq_small, Mass=Mass1, v_f=v_fermi, E0_model=1) 	! module "CDF_Ritchi"
      case (2)	! Ritchi model
         ImE = Diel_func(CDF%A(:), CDF%E0(:), CDF%Gamma(:), W, dq_small, Mass=Mass1, E0_model=2) 	! module "CDF_Ritchi"
      end select
      dSigma_mid = (Zeff*Zeff)*AB*ImE
      
      ! And on the next half-step forward:
      Q = Q + dQ_half	! half step forward
      AB = get_AB(Q, Mc2, W, beta, Emc2)	! below
      ! Effective charge:
      if (present(Zmodel)) then	! set effective charge of the incident SHI
         Z = Zeff	! Equilibrium_charge_SHI(Ekin, Mass, Zeff, Ztarget, perc, Kind_Zeff, fixed_Zeff)	! module "SHI_charge_state"
      else
         Z = Zeff
      endif
      ! Get momentum transfer from the recoil energy:
      dq_small = momentum_from_kinetic_energy(Q, Mass)	! module "Relativity"
      ! Get the effective mass of the scattering center:
      Mass1 = chose_mass(CDF_m_eff, M_sc, dq_small, k, Eff_m) 	! below
      select case (CDF_dispers)	! model for the dispertion relation of the target:
      case default	! free electron
         ImE = Diel_func(CDF%A(:), CDF%E0(:), CDF%Gamma(:), W, dq_small, Mass=Mass1) 	! module "CDF_Ritchi"
      case (1)	! plasmon pole
         ImE = Diel_func(CDF%A(:), CDF%E0(:), CDF%Gamma(:), W, dq_small, Mass=Mass1, v_f=v_fermi, E0_model=1) 	! module "CDF_Ritchi"
      case (2)	! Ritchi model
         ImE = Diel_func(CDF%A(:), CDF%E0(:), CDF%Gamma(:), W, dq_small, Mass=Mass1, E0_model=2) 	! module "CDF_Ritchi"
      end select
      dSigma_cur = (Zeff*Zeff)*AB*ImE
      
!        write(*,'(a,es,es,es,es,es,es)') 'Im', Q, dQ, dq_small, AB, W, ImE
      
      dSigma_sum = dSigma_sum +  dQ/6.0d0*(dSigma_cur + 4.0d0*dSigma_mid + dSigma0)
      ! save for the next integration step:
      dSigma0 = dSigma_cur
      Q0 = Q
   enddo ! while (Q < Qmax)
   dSigma = dSigma_sum
!     PAUSE 'dCDF_CS_dQ'
end function dCDF_CS_dQ


pure function get_AB(Q, Mc2, W, beta, Emc2) result(AB)
   real(8) AB	! multiplicator in calculations of d2sigma/dWdQ
   real(8), intent(in) :: Q	! [eV] recoil energy (related to transferred momentum)
   real(8), intent(in) :: Mc2	! rest energy of incident particle [eV]
   real(8), intent(in) :: W	! energy loss [eV]
   real(8), intent(in) :: beta	! relativistic beta factor of incident particle
   real(8), intent(in) :: Emc2	! [eV] relativistic total energy of the scattering center
   real(8) :: Qmc2 , W2, beta2, QQmc2, cosbra, cos2teta, dem, A, B
   Qmc2 = Q + Mc2	! [eV]
   ! Get the cos^2(teta):
   W2 = W*W	! [eV^2]
   beta2 = beta*beta
   QQmc2 = Q*Qmc2	! [eV^3]
   cosbra = 1.0d0 + (QQmc2 - W2)/(2.0d0*W*Emc2)
   cos2teta = W2/beta2/QQmc2*cosbra*cosbra	! [1] P. 282, Eq.(A.42)
   dem = QQmc2 - W2
   A = 2.0d0*Mc2/(W*QQmc2) + 2.0d0*beta2*(1.0d0 - cos2teta)*W*Mc2 / (dem*dem)
   B =  W*Qmc2/Mc2
   AB = A*B
end function get_AB


function chose_mass(CDF_m_eff, M_sc, dq, k, Eff_m) result(Mass1)
   real(8) Mass1	! [kg]
   integer, intent(in) :: CDF_m_eff		! index of the model for effective mass
   real(8), intent(in) :: M_sc	! [kg] mass, used in case of some models
   real(8), intent(in) :: dq	! [sqqrt(J*kg)] transferred momentum
   real(8), dimension(:), intent(in) :: k, Eff_m	! dispersion curve from DOS
   real(8) :: qlim
   integer :: j
   ! Find effective mass of the scattering center:
   select case (CDF_m_eff)
   case default	! free particle (electron, atom, hole...)
      Mass1 = M_sc	! [kg]
   case (0)	! effective electron mass from DOS, Eq.(6) [6]
      ! Find the index of the mass from DOS:
      qlim = dq*sqrt(g_e)	! convert into [sqrt(eV*kg)]
      if (qlim <= k(size(k))) then	! only if the transferred energy is within our awailable DOS
         call find_in_array_monoton(k, qlim, j)	! module "Little_subroutines"
         Mass1 = Eff_m(j)*g_me	! [kg]
      else	! otherwise, assume free electron
         Mass1 = g_me	! [kg] free electron mass at rest
      endif
   end select
end function chose_mass



!-----------------------------------------------------------------------------
! Integral RBEB cross section from [5]
pure function RBEB_total_CS(Ekin, Nel, Ip, U) result(sigma)
   real(8) sigma	! [A^2] cross section of inelastic electron scattering on atomic shell
   real(8), intent(in) :: Ekin	! [eV] kinetic energy of incident electron
   real(8), intent(in) :: Nel	! number of electrons in the chosen shell
   real(8), intent(in) :: Ip	! [eV] ionization potential of this shell
   real(8), intent(in) :: U	! [eV] kinetic energy of electrons in this shell
   real(8) :: Erest, bp, t, tp, up, betat2, prefac, tempt, A, B, C, temp
   
   if (Ekin > Ip) then	! there is a cross section
      Erest = rest_energy(g_me)	! module "Relativity"
      t = Ekin/Ip	! [5] under Eq.(3)
      tp = Ekin/Erest	! [5] Eq.(12)
      bp = Ip/Erest	! [5] Eq.(13)
      up = U/Erest	! [5] Eq.(14)
      tempt = 1.0d0 + tp
      betat2 = 1.0d0 - 1.0d0/(tempt*tempt)
      prefac = RBEB_prefac(Nel, tp, up, bp) ! below
      A = 0.5d0 * (1.0d0 - 1.0d0/(t*t)) * (log((betat2)/(1.0d0 - betat2)) - betat2 - log(2.0d0*bp))
      temp = 1.0d0 + tp*0.5d0
      temp = temp*temp
      B = 1.0d0 - 1.0d0/t - log(t)/(t + 1.0d0)*(1.0d0 + 2.0d0*tp)/temp
      C = bp*bp/temp*(t - 1.0d0)*0.5d0
      sigma = prefac * (A + B + C)	! [A^2] from [5] Eq.(22)
   else	! no ionization for incident energies below Ip
      sigma = 0.0d0
   endif
end function RBEB_total_CS


pure function RBEB_integral_CS(Ekin, Nel, Ip, U, Emax) result(sigma)
   real(8) sigma	! [A^2] cross section of inelastic electron scattering on atomic shell
   real(8), intent(in) :: Ekin	! [eV] kinetic energy of incident electron
   real(8), intent(in) :: Nel	! number of electrons in the chosen shell
   real(8), intent(in) :: Ip	! [eV] ionization potential of this shell
   real(8), intent(in) :: U	! [eV] kinetic energy of electrons in this shell
   real(8), intent(in) :: Emax  ! [eV] maximal transferred energy
   real(8) :: tmax, Wmax, Erest, bp, t, tp, up, betat2, prefac, tempt, A, B, C, temp, tmax1, t_tmax
   
   if (Ekin > Ip) then	! there is a cross section
      t = Ekin/Ip	! [5] under Eq.(3)
      tmax = (t-1.0d0)*0.5d0   ! maximal transfered energy
      Wmax = Emax/Ip
      if (Wmax < tmax) then
         tmax = Wmax   ! maximal tranfered energy
      endif
      Erest = rest_energy(g_me)	! module "Relativity"
      t = Ekin/Ip  ! [5] under Eq.(3)
      tp = Ekin/Erest  ! [5] Eq.(12)
      bp = Ip/Erest     ! [5] Eq.(13)
      up = U/Erest     ! [5] Eq.(14)
      prefac = RBEB_prefac(Nel, tp, up, bp) ! below
      tempt = 1.0d0 + tp
      betat2 = 1.0d0 - 1.0d0/(tempt*tempt)
      tmax1 = tmax + 1.0d0
      t_tmax = t - tmax
      A = 0.5d0 * (log((betat2)/(1.0d0 - betat2)) - betat2 - log(2.0d0*bp)) * tmax * ( (tmax + 2.0d0)/(tmax1*tmax1) + (2.0d0*t - tmax)/(t_tmax*t_tmax) )
      temp = 1.0d0 + 0.5d0*tp
      temp = temp*temp
      B = tmax * ( bp*bp/(temp) + (t*t + tmax*(1.0d0 - t) + 1.0d0)/(t*t_tmax*tmax1) )
      C = -1.0d0/tmax1 * (1.0d0 + 2.0d0*tp)/temp * log(t*tmax1/t_tmax)
      sigma = prefac * (A + B + C)	! [A^2] from [5] Eq.(22)
   else ! no ionization for incident energies below Ip
      sigma = 0.0d0
   endif
end function RBEB_integral_CS


pure function RBEB_prefac(Nel, tp, up, bp) result (prefac)
   real(8) prefac   ! Eq.(22) from [5] 
   real(8), intent(in) :: Nel, tp, up, bp
   real(8) :: tempt, betat2, tempu, betau2, tempb, betab2
   tempt = 1.0d0 + tp
   betat2 = 1.0d0 - 1.0d0/(tempt*tempt)
   tempu = 1.0d0 + up
   betau2 =  1.0d0 - 1.0d0/(tempu*tempu)
   tempb = 1.0d0 + bp
   betab2 = 1.0d0 - 1.0d0/(tempb*tempb)
   prefac = 4.0d0*g_Pi*g_a0*g_a0*(g_alpha**4)*Nel/((betat2 + betau2 + betab2)*2.0d0*bp)
end function RBEB_prefac


end module CS_electrons_inelastic
