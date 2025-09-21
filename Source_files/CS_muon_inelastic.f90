! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2025
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains all cross sections for muon inelastic scattering
! References used:
! [1] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2008: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2009)

module CS_muon_inelastic
use Universal_constants
use CDF_delta, only : energy_loss_delta, CDF_total_CS_delta
use Objects
use Read_input_data, only: m_input_folder, m_folder_materials, m_muon_CS, m_muon_inelast_CS
use Dealing_with_files, only: read_file
use CS_general_tools, only: MFP_from_sigma
use CS_integration_limits, only: find_Wmax_equal_Wmin, W_max
use CS_electrons_inelastic, only: CDF_total_CS, CDF_total_CS_nonrel
use Little_subroutines, only: interpolate_data_single
use Relativity, only: rest_energy

implicit none

 contains

!-----------------------------------------------------------------------------
! Subroutines for inelastic scattering:

subroutine get_muon_IMFP(Material, numpar, Err)
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
   
   write(*, '(a)', advance='no') ' Obtaining muon inelastic scattering cross sections...'
   
   path_sep => numpar%path_sep
   Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_muon_CS))	! Muon CSs are storred in this folder
   
   N_targets = size(Material)	! that's how many different targets user specified
   
   ! Get muon MFPs:
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
         write(temp_c,'(f8.1)') Material(i)%T    ! target temperature
         write(temp_c2,'(f8.1)') numpar%CDF_Eeq_factor
         ! Set part of the file name corresponding to the model used for muon inelastic CS:
         select case (numpar%CDF_model)
         case default   ! Ritchie
            write(temp_c3,'(a)') '_R_'    ! Ritchie
         case (1)   ! Mermin
            write(temp_c3,'(a)') '_M_'    ! Mermin
         endselect
         
         select case (numpar%Mu_inelast)
         case (1:3)	! CDF with delta-functions : FOR NOW, ONLY THIS MODEL IS IMPLEMENTED
            Model_name = 'CDF_delta_T_'//trim(adjustl(temp_c))//'K'//trim(adjustl(temp_c3))//'Eeq_'//trim(adjustl(temp_c2))
         case (5)	! CDF with single pole delta-functions
            Model_name = 'CDF_SP_delta_T_'//trim(adjustl(temp_c))//'K'//trim(adjustl(temp_c3))//'Eeq_'//trim(adjustl(temp_c2))
         case default	! exclude
            Model_name = 'NO'
         end select
         
         ! Path to files with the valence band and total (material specific) MFPs:
         Path_valent = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
         Path_valent = trim(adjustl(Path_valent))//path_sep//trim(adjustl(Material(i)%Name))
         
         ! For muon CSs use the same grid as for photons:
         Ngrid = size(Element%Phot_absorption%E)
         allocate(Element%Muon_inelastic%E(Ngrid))
         Element%Muon_inelastic%E = Element%Phot_absorption%E
         allocate(Element%Muon_inelastic%Total(Ngrid))
         allocate(Element%Muon_inelastic%Total_MFP(Ngrid))
         allocate(Element%Muon_inelastic%Total_Se(Ngrid))
         allocate(Element%Muon_inelastic%Per_shell(Element%N_shl,Ngrid))
         allocate(Element%Muon_inelastic%MFP(Element%N_shl,Ngrid))
         allocate(Element%Muon_inelastic%Se(Element%N_shl,Ngrid))
         N_shells = Element%N_shl
         
         ! Calculate total cross sections for all shells of this element:
         do k = 1, N_shells
            ! Check valence band:
            ! If this shell forms valent band, and data for the valence band exists, separate it from the atomic shells:
            VAL:if ((numpar%Mu_inelast /= 2) .and.  (Element%valent(k)) .and. (allocated(Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
               if (.not.allocated(Material(i)%Muon_inelastic_valent%E)) then   ! valence band has not been defined yet, so do that
                  Nsiz = size(Material(i)%Ph_absorption_valent%E)   ! valence band has a different grid from core shells
                  allocate(Material(i)%Muon_inelastic_valent%E(Nsiz))
                  Material(i)%Muon_inelastic_valent%E = Material(i)%Ph_absorption_valent%E    ! we already set this grid for photons, reuse it
                  allocate(Material(i)%Muon_inelastic_valent%Total(Nsiz))
                  allocate(Material(i)%Muon_inelastic_valent%Total_MFP(Nsiz))
                  allocate(Material(i)%Muon_inelastic_valent%Total_Se(Nsiz))
                  !allocate(Material(i)%Muon_inelastic_valent%Per_shell(Element%N_shl,Nsiz)) ! --- leave it not allocated !!!
                  ! Check file with CSs:
                  !File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_muon_inelast_CS))//'_valence_'//trim(adjustl(Model_name))//'.dat'
                  File_name = trim(adjustl(Path_valent))//path_sep//trim(adjustl(m_muon_inelast_CS))//'_valence_'//trim(adjustl(Model_name))//'.dat'
                  inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!                   if (file_exist) then	! just read from the file, no need to recalculate:
                  if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
                     open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
                     ! Get the MFP and CDF points:
                     count_lines = 0
                     do m = 1, Nsiz	! for all energy grid points:
                        E => Material(i)%Muon_inelastic_valent%E(m)	! muon energy [eV]
                        read(FN,'(es,es,es,es)', IOSTAT=Reason) E, Material(i)%Muon_inelastic_valent%Total(m), &
                              Material(i)%Muon_inelastic_valent%Total_MFP(m), Material(i)%Muon_inelastic_valent%Total_Se(m)
                        call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                        if (.not.read_well) then
                           close(FN)	! redo the file
                           goto 8391
                        endif
                        ! Recalculate the MFP since material density might be different from the one saved in the file:
                        Material(i)%Muon_inelastic_valent%Total_MFP(m) = MFP_from_sigma(Material(i)%Muon_inelastic_valent%Total(m),  Material(i)%At_Dens) ! [A] module "CS_general_tools"
                     enddo
                  else	! no such file => create it
8391                 open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
                     do m = 1, Nsiz	! for all energy grid points:
                        E => Material(i)%Muon_inelastic_valent%E(m)	! muon energy [eV]
                        !call get_muon_inelastic_CS(E, Material(i), Material(i)%DOS, numpar, j, k, Material(i)%DOS%Egap, sigma, CDF_dispers=0, CDF_m_eff=1, Se=Se)	! see below
                        call get_muon_inelastic_CS(E, Material(i), Material(i)%DOS, numpar, j, k, Material(i)%DOS%Egap, sigma, &
                                                          CDF_dispers=numpar%CDF_dispers, CDF_m_eff=numpar%CDF_m_eff, Se=Se)    ! see below
                        Material(i)%Muon_inelastic_valent%Total(m) = sigma   ! [A^2]
                        lambda = MFP_from_sigma(sigma,  Material(i)%At_Dens)    ! module "CS_general_tools"
                        Material(i)%Muon_inelastic_valent%Total_MFP(m) = lambda  ! [A]
                        Material(i)%Muon_inelastic_valent%Total_Se(m) = Se    ! [eV/A]
                        write(FN,'(es,es,es,es)') Material(i)%Muon_inelastic_valent%E(m), sigma, lambda, Se
                     enddo ! m = 1, Nsiz
                  endif ! (file_exist)
                  close(FN)
               endif ! (.not.allocated(Material(i)%Muon_inelastic_valent))
            else VAL    ! core shell
               ! Check file with CSs:
               File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_muon_inelast_CS))//'_'//trim(adjustl(Element%Shell_name(k)))//'_'//trim(adjustl(Model_name))//'.dat'
               inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!                if (file_exist) then	! just read from the file, no need to recalculate:
               if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
                  open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
                  ! Get the MFP and CDF points:
                  count_lines = 0
                  do m = 1, Ngrid	! for all energy grid points:
                     E => Element%Muon_inelastic%E(m)	! muon energy [eV]
                     read(FN,'(es,es,es,es)', IOSTAT=Reason) E, Element%Muon_inelastic%Per_shell(k,m), Element%Muon_inelastic%MFP(k,m), Element%Muon_inelastic%Se(k,m)
                     call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                     if (.not.read_well) then
                        close(FN)	! redo the file
                        goto 8390
                     endif
!                       print*, k, m, Element%Muon_inelastic%E(m), Element%Muon_inelastic%Per_shell(k,m)
                  enddo
               else	! no such file => create it
8390           open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
                  do m = 1, Ngrid	! for all energy grid points:
                     E => Element%Muon_inelastic%E(m)	! muon energy [eV]
                     !call get_muon_inelastic_CS(E, Material(i), Material(i)%DOS, numpar, j, k, Material(i)%Elements(j)%Ip(k), sigma, CDF_dispers=0, CDF_m_eff=1, Se=Se)	! see below
                     call get_muon_inelastic_CS(E, Material(i), Material(i)%DOS, numpar, j, k, Material(i)%Elements(j)%Ip(k), sigma, &
                                                          CDF_dispers=numpar%CDF_dispers, CDF_m_eff=numpar%CDF_m_eff, Se=Se)    ! see below
                     Element%Muon_inelastic%Per_shell(k,m) = sigma	! [A^2]
                     Element%Muon_inelastic%Se(k,m) = Se * 1.0d24/Material(i)%At_Dens    ! [eV/A]
                     Element%Muon_inelastic%MFP(k,m) = MFP_from_sigma(Element%Muon_inelastic%Per_shell(k,m),  1.0d24)    ! module "CS_general_tools"
                     write(FN,'(es,es,es,es)') Element%Muon_inelastic%E(m), sigma, Element%Muon_inelastic%MFP(k,m), Element%Muon_inelastic%Se(k,m)
                  enddo ! m = 1, Ngrid
               endif ! (file_exist)
               close(FN)
               ! Normalize MFPs and Se to real material density:
               ! Take into account the density of atoms of this particular kind:
               elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
               Element%Muon_inelastic%MFP(k,:) = Element%Muon_inelastic%MFP(k,:) * 1.0d24/(Material(i)%At_Dens * elem_contrib)
               Element%Muon_inelastic%Se(k,:) = Element%Muon_inelastic%Se(k,:) * Material(i)%At_Dens/1.0d24 * elem_contrib
            endif VAL
         enddo ! k = 1, N_shells
     
      enddo LMNT
     
      ! And the total CS:
      Nsiz = size(Material(i)%Ph_absorption_total%E)   ! grid for total cross section is different from core-shells grid (but the same as valence grid)
      allocate(Material(i)%Muon_inelastic_total%E(Nsiz))
      Material(i)%Muon_inelastic_total%E = Material(i)%Ph_absorption_total%E
      allocate(Material(i)%Muon_inelastic_total%Total(Nsiz))
      allocate(Material(i)%Muon_inelastic_total%Total_MFP(Nsiz))
      allocate(Material(i)%Muon_inelastic_total%Total_Se(Nsiz))
      Material(i)%Muon_inelastic_total%Total(:) = 0.0d0
      Material(i)%Muon_inelastic_total%Total_MFP(:) = 0.0d0
      Material(i)%Muon_inelastic_total%Total_Se(:) = 0.0d0
      ! Sum up CSs from each element:
      do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         do m = 1, Nsiz	! for all energy grid points:
            E => Material(i)%Muon_inelastic_total%E(m)    ! muon energy [eV]
            ! Core shells:
            sigma = 0.0d0 ! to start with
            Se = 0.0d0
            lambda = 0.0d0
            ! Get the cross section for this element on this grid point (which may not coincide with the grid for this element due to element-specific grid points):
            do k = 1, Element%N_shl		! for all shells in this elements
               if (.not.Element%valent(k)) then ! core orbitals
                  ! For each shell, find its partial cross section:
                  call interpolate_data_single(Element%Muon_inelastic%E,  Element%Muon_inelastic%Per_shell(k,:), E, sigma_cur) ! module "Little_subroutines"
                  call interpolate_data_single(Element%Muon_inelastic%E,  Element%Muon_inelastic%Se(k,:), E, Se_cur) ! module "Little_subroutines"
                  call interpolate_data_single(Element%Muon_inelastic%E,  Element%Muon_inelastic%MFP(k,:), E, lambda_cur) ! module "Little_subroutines"
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
            Material(i)%Muon_inelastic_total%Total(m) = Material(i)%Muon_inelastic_total%Total(m) + sigma   ! [A^2]
            Material(i)%Muon_inelastic_total%Total_Se(m) = Material(i)%Muon_inelastic_total%Total_Se(m) +  Se   ! [eV/A]
            Material(i)%Muon_inelastic_total%Total_MFP(m) = Material(i)%Muon_inelastic_total%Total_MFP(m) + lambda  ! [1/A] to be inversed below
            
            ! Add valence band / shells:
            if ((numpar%Mu_inelast /= 2) .and. (allocated(Material(i)%CDF_valence%A)) ) then    ! valence band is defined and used in the CS model (RBEB excluded)
               if (j == 1) then ! add valence band only once - when studying the first element of the material
                  Material(i)%Muon_inelastic_total%Total(m) = Material(i)%Muon_inelastic_total%Total(m) + Material(i)%Muon_inelastic_valent%Total(m)   ! [A^2]
                  Material(i)%Muon_inelastic_total%Total_MFP(m) = Material(i)%Muon_inelastic_total%Total_MFP(m) + 1.0d0/Material(i)%Muon_inelastic_valent%Total_MFP(m)   ! [1/A]
                  Material(i)%Muon_inelastic_total%Total_Se(m) = Material(i)%Muon_inelastic_total%Total_Se(m) + Material(i)%Muon_inelastic_valent%Total_Se(m)   ! [eV/A]
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
                     call interpolate_data_single(Element%Muon_inelastic%E,  Element%Muon_inelastic%Per_shell(k,:), E, sigma_cur) ! module "Little_subroutines"
                     call interpolate_data_single(Element%Muon_inelastic%E,  Element%Muon_inelastic%Se(k,:), E, Se_cur) ! module "Little_subroutines"
                     call interpolate_data_single(Element%Muon_inelastic%E,  Element%Muon_inelastic%MFP(k,:), E, lambda_cur) ! module "Little_subroutines"
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
               Material(i)%Muon_inelastic_total%Total(m) = Material(i)%Muon_inelastic_total%Total(m) + sigma   ! [A^2]
               Material(i)%Muon_inelastic_total%Total_Se(m) = Material(i)%Muon_inelastic_total%Total_Se(m) +  Se   ! [eV/A]
               Material(i)%Muon_inelastic_total%Total_MFP(m) = Material(i)%Muon_inelastic_total%Total_MFP(m) + lambda  ! [1/A] to be inversed below
            endif
         enddo ! m = 1, Nsiz
      enddo !  j =1, N_elements
      ! Inverse it to get the MFP:
      Material(i)%Muon_inelastic_total%Total_MFP = 1.0d0/Material(i)%Muon_inelastic_total%Total_MFP  ! [1/A] -> [A]
      ! Save the total cross section into the file:
      File_name = trim(adjustl(Path_valent))//path_sep//trim(adjustl(m_muon_inelast_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
      if (.not.file_exist .or. numpar%recalculate_MFPs) then	! only create it if file does not exist
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            write(FN,'(es,es,es,es)') Material(i)%Muon_inelastic_total%E(m), Material(i)%Muon_inelastic_total%Total(m), &
                                                  Material(i)%Muon_inelastic_total%Total_MFP(m), Material(i)%Muon_inelastic_total%Total_Se(m)
         enddo ! m = 1, Nsiz
         close(FN)
      endif
   enddo TRGT
   
   write(*, '(a)') ' Done.'
!    print*, 'Muon inelastic scattering cross sections are obtained.'
   
   nullify(path_sep, E, Element)
end subroutine get_muon_IMFP




! Interface to select the model of muon inelastic scattering cross section:
function get_inelastic_energy_transfer_muon(Ee, Material, DOS, numpar, j, k, Ip, CDF_dispers, CDF_m_eff, hw_phonon) result(dE)
   real(8) :: dE   ! [eV] sampled transferred energy
   real(8), intent(in) :: Ee	! [eV] muon kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   integer, intent(in), optional :: CDF_dispers, CDF_m_eff	! dispersion relation and effective mass
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency, for particle-phonon scattering
   !---------------------------------------
   real(8) :: eps, eps1, eps2, eps3
   real(8) :: Zeff, max_E0, Eeq, RN, CS_sampled, CS_cur, E_left, E_right, E_cur, mtc2, Se1, m_in_c2
   real(8) :: CS_tot    ! [A^2] precalculated total cross section
   integer :: dispers, m_eff, Mu_inelast
   ! Set accepteble margin of precision for the angle:
   eps1 = 1.0d-2
   eps2 = 1.0d-4
   eps3 = 1.0d-6
   
   ! Set the model parameters:
   Zeff = -1.0d0	! muon charge
   Mu_inelast = numpar%Mu_inelast   ! to chose the model of CS calculations below
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
   
   m_in_c2 = rest_energy(g_M_muon)   ! module "Relativity"
   mtc2 = rest_energy(g_me)   ! module "Relativity"

   ! In case it is delta-CDF model, we have to make sure the energy is not below its applicability limit:
   if (Mu_inelast == 3) then
      !VAL0:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
      VAL0:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         max_E0 = maxval(Material%CDF_valence%E0(:))
      else VAL0
         max_E0 = maxval(Material%Elements(j)%CDF(k)%E0(:))
      endif VAL0

      call find_Wmax_equal_Wmin(m_in_c2, mtc2, .false., Ee, Ip, max_E0, Eeq)   ! module "CS_integration_limits"
      ! Check if delta-functional CDF works here:
!       if (Ee < 1.1d0*Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
      if (Ee < numpar%CDF_Eeq_factor * Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
         Mu_inelast = 4
      endif
   endif ! (Mu_inelast == 3)
   
   E_left = Ip ! [eV] minimal transferred energy
!    E_right = (Ip + Ee)*0.5d0    ! [eV] maximal transferred energy
   E_right = W_max(m_in_c2, mtc2, .false., Ee, Ip)  ! module "CS_integration_limits"
   ! * Note: in case of delta-CDF, the limits will automatically be adjusted to Wmin, Wmax
   ! Use of the maximal plasmon energy as the upper limit of integration is not included yet !!!         
   
   ! Sample the cross section:
   call random_number(RN)
   CS_tot = get_integral_inelastic_CS_muon(Ee, Material, DOS, numpar, j, k, Ip, Mu_inelast, Zeff, dispers, m_eff, E_right)   ! below
   CS_sampled = RN*CS_tot
   
   ! For analytucally integrable CDF, use bisection method:
   if (Mu_inelast /= 4) then
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
      CS_cur = get_integral_inelastic_CS_muon(Ee, Material, DOS, numpar, j, k, Ip, Mu_inelast, Zeff, dispers, m_eff, E_cur)   ! below
      ! Search by bisection method:
!       do while (ABS(CS_cur - CS_sampled)/CS_sampled > eps)
      do while (abs(E_left - E_right) >= eps1*E_left)
         if (CS_cur > CS_sampled) then
            E_right = E_cur
         else
            E_left = E_cur
         endif
         E_cur = (E_left + E_right)/2.0d0
         if (abs(E_left - E_right) < eps1) exit  ! precise enough
         CS_cur = get_integral_inelastic_CS_muon(Ee, Material, DOS, numpar, j, k, Ip, Mu_inelast, Zeff, dispers, m_eff, E_cur)   ! below
      enddo
!       print*, 'get_inelastic_energy_transfer_muon 1', Mu_inelast, CS_sampled, CS_cur, E_cur
   else     ! for numerically integrable one, just integrate until found the smapled value:
      VAL3:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                    Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model,&
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, &
                     Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, CDF_dispers=numpar%CDF_dispers, &
                     Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         end select ! (numpar%CDF_dispers) 
      else VAL3 ! core shells
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                    Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .false.,0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, &
                     Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, CDF_dispers=numpar%CDF_dispers, &
                     Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
         end select ! (numpar%CDF_dispers) 
      endif VAL3
!       print*, 'get_inelastic_energy_transfer_muon 2', Mu_inelast, CS_sampled, CS_cur, E_cur
   endif    ! (El_inelast /= 4)

   ! Output: sampled transferred energy:
   dE = E_cur
end function get_inelastic_energy_transfer_muon


function get_integral_inelastic_CS_muon(Ee, Material, DOS, numpar, j, k, Ip, Mu_inelast, Zeff, dispers, m_eff, Emax) result (sigma)
   real(8) :: sigma	! [A^2] cross section
   real(8), intent(in) :: Ee	! [eV] muon kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   integer, intent(in) :: Mu_inelast    ! model for the inelastic cross section
   real(8), intent(in) :: Zeff  ! effective charge of the incident particle
   integer, intent(in) :: dispers, m_eff	! dispersion relation and effective mass
   real(8), intent(in) :: Emax  ! [eV] upper integration limit
   !---------------------------------------
   real(8) :: Se1

   ! Now choose the model of CS:
   select case (Mu_inelast)  ! chose which model for muon inelastic cross section to use
   case (1:3,5)	! CDF with delta-functions
      !VAL2:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
      VAL2:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         sigma = CDF_total_CS_delta(Mu_inelast, Ee, g_M_muon, Zeff, Ip, Material%At_Dens, Material%CDF_valence, g_me, .false., Emax_in = Emax)    ! module "CDF_delta"
      else VAL2  ! core shells:
         sigma = CDF_total_CS_delta(Mu_inelast, Ee, g_M_muon, Zeff, Ip, Material%At_Dens, Material%Elements(j)%CDF(k), g_me, .false., Emax_in = Emax)    ! module "CDF_delta"
      endif VAL2
   
   case (4) ! nonrelativistic Ritchie CDF
     !VAL3:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
     VAL3:if ( (k == 0) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                    Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         end select ! (numpar%CDF_dispers) 
      else VAL3 ! core shells
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                    Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         end select ! (numpar%CDF_dispers) 
      endif VAL3
   case default	! exclude
      sigma = 0.0d0
   end select
end function get_integral_inelastic_CS_muon




! Interface to select the model of muon inelastic scattering cross section:
 subroutine get_muon_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, sigma, CDF_dispers, CDF_m_eff, Se)
   real(8), intent(in) :: Ee	! [eV] muon kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   real(8), intent(out) :: sigma	! [A^2] cross section
   real(8), intent(out), optional :: Se   ! [eV/A] energy loss
   integer, intent(in), optional :: CDF_dispers, CDF_m_eff	! dispersion relation and effective mass
   !---------------------------------------
   real(8) :: Zeff, Eeq, max_E0, Se1, m_in_c2, mtc2
   integer :: dispers, m_eff, Mu_inelast
   Zeff = 1.0d0   ! effective charge (electron)
   Mu_inelast = numpar%Mu_inelast   ! to chose the model of CS calculations below
   
   ! In case it is delta-CDF model, we have to make sure the energy is not below its applicability limit:
  select case (Mu_inelast)
   case (1:3,5)	! CDF with delta-functions
      VAL0:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         max_E0 = maxval(Material%CDF_valence%E0(:))
      else VAL0
         max_E0 = maxval(Material%Elements(j)%CDF(k)%E0(:))
      endif VAL0
      m_in_c2 = rest_energy(g_M_muon)   ! module "Relativity"
      mtc2 = rest_energy(g_me)   ! module "Relativity"
      call find_Wmax_equal_Wmin(m_in_c2, mtc2, .false., Ee, Ip, max_E0, Eeq)   ! module "CS_integration_limits"
      ! Check if delta-functional CDF works here:
      if (Ee < numpar%CDF_Eeq_factor * Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
         Mu_inelast = 4
      endif
   end select
   
   
   ! Select the appropriate model:
   select case (Mu_inelast)
   case (1:3,5)	! CDF with delta-functions
       VAL2:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         sigma = CDF_total_CS_delta(Mu_inelast,Ee, g_M_muon, Zeff, Ip, Material%At_Dens, Material%CDF_valence, g_me, .false.)    ! module "CDF_delta"
         if (present(Se)) Se = energy_loss_delta(Mu_inelast,Ee, g_M_muon, 1.0d0, Zeff, Ip, Material%At_Dens, g_me, &
            Material%CDF_valence, .false., 0) ! module "CDF_delta"
      else VAL2  ! core shells:
         sigma = CDF_total_CS_delta(Mu_inelast,Ee, g_M_muon, Zeff, Ip, Material%At_Dens, Material%Elements(j)%CDF(k), g_me, .false.)    ! module "CDF_delta"
         if (present(Se)) Se = energy_loss_delta(Mu_inelast,Ee, g_M_muon, 1.0d0, Zeff, Ip, Material%At_Dens, g_me, &
            Material%Elements(j)%CDF(k), .false., 0) ! module "CDF_delta"
      endif VAL2
   
   case (4)  ! Non-relativistic Ritchie CDF
      VAL3:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers)  ! module "CS_electrons_inelastic"
         end select ! (numpar%CDF_dispers) 
      else VAL3 ! core shells
         select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, &
                    DOS%k, DOS%Eff_m, .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model)  ! module "CS_electrons_inelastic"
         case (1)  ! Plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), &
                     g_me, DOS%k, DOS%Eff_m, &
                     .false.,0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers, v_f=Material%DOS%v_f)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie extended
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_M_muon, Zeff, Ip, Material%T_eV, Material%Elements(j)%CDF(k), g_me, &
                     DOS%k, DOS%Eff_m, &
                     .false., 0.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                     CDF_dispers=numpar%CDF_dispers)  ! module "CS_electrons_inelastic"
         end select ! (numpar%CDF_dispers) 
      endif VAL3
      if (present(Se)) Se = Se1   

   case default	! exclude
      sigma = 0.0d0
      if (present(Se)) Se = 0.0d0
   end select
end subroutine get_muon_inelastic_CS


end module CS_muon_inelastic
