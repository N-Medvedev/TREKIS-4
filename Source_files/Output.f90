! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! and R. Rymzhanov
! in 2018-2025
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to create output files:

MODULE Output

! Open_MP related modules from external libraries:
#ifdef _OPENMP
   USE IFLPORT
   USE OMP_LIB
#endif

use Objects
use Universal_constants
use Dealing_with_files, only: copy_file, close_file, get_file_stat, Count_lines_in_file, read_file
use Dealing_with_XYZ_files, only: write_XYZ
use Dealing_with_LAMMPS, only: Write_LAMMPS_input_file
use Gnuplotting
use Little_subroutines, only: print_time, find_order_of_number, print_energy, order_of_energy, m_starline, m_dashline
use Read_input_data, only: m_input_minimal, m_input_data, m_numerical_parameters, m_input_folder, m_databases, &
                           m_EADL, m_EPDL
use Periodic_table, only : Atomic_data, Read_from_periodic_table
use Initial_conditions, only: count_types_of_SHIs, repeated_type
use CS_general_tools, only: get_ranges_from_Se

implicit none

! In this module, all output file names are collected:
character(200) :: m_output_parameters, m_communication, m_output_DOS, m_output_DOS_k, m_output_DOS_effm
character(200) :: m_folder_MFP, m_output_MFP, m_output_IMFP, m_output_EMFP, m_output_Brems, m_output_annihil
character(200) :: m_output_Compton, m_output_Rayleigh, m_output_pair, m_output_absorb
character(200) :: m_output_Se, m_output_Range, m_output_Se_vs_range

character(200) :: m_output_total, m_output_N_gnu, m_output_E_gnu, m_output_total_cutoff
character(200) :: m_output_MD, m_output_MD_T_gnu, m_output_MD_MSD_gnu, m_output_MD_E_gnu
character(200) :: m_output_spectrum_ph, m_output_spectrum_e, m_output_spectrum_h, m_output_spectrum_p, m_output_spectrum_SHI, m_output_spectrum_mu
character(200) :: m_output_spectrum_ph_1d, m_output_spectrum_e_1d, m_output_spectrum_h_1d, m_output_spectrum_p_1d, m_output_spectrum_SHI_1d, &
                  m_output_spectrum_mu_1d
character(200) :: m_output_theta_ph_1d, m_output_theta_e_1d, m_output_theta_h_1d, m_output_theta_p_1d, m_output_theta_SHI_1d, m_output_theta_mu_1d

character(200) :: m_output_cartesian_1d_X_ph, m_output_cartesian_1d_Y_ph, m_output_cartesian_1d_Z_ph
character(200) :: m_output_cartesian_1d_X_e, m_output_cartesian_1d_Y_e, m_output_cartesian_1d_Z_e
character(200) :: m_output_cartesian_1d_X_p, m_output_cartesian_1d_Y_p, m_output_cartesian_1d_Z_p
character(200) :: m_output_cartesian_1d_X_h, m_output_cartesian_1d_Y_h, m_output_cartesian_1d_Z_h
character(200) :: m_output_cartesian_1d_X_SHI, m_output_cartesian_1d_Y_SHI, m_output_cartesian_1d_Z_SHI
character(200) :: m_output_cartesian_1d_X_a, m_output_cartesian_1d_Y_a, m_output_cartesian_1d_Z_a
character(200) :: m_output_cartesian_1d_X_mu, m_output_cartesian_1d_Y_mu, m_output_cartesian_1d_Z_mu
character(200) :: m_output_cartesian_1d_X_E_ph, m_output_cartesian_1d_Y_E_ph, m_output_cartesian_1d_Z_E_ph
character(200) :: m_output_cartesian_1d_X_E_e, m_output_cartesian_1d_Y_E_e, m_output_cartesian_1d_Z_E_e
character(200) :: m_output_cartesian_1d_X_E_p, m_output_cartesian_1d_Y_E_p, m_output_cartesian_1d_Z_E_p
character(200) :: m_output_cartesian_1d_X_E_h, m_output_cartesian_1d_Y_E_h, m_output_cartesian_1d_Z_E_h
character(200) :: m_output_cartesian_1d_X_E_SHI, m_output_cartesian_1d_Y_E_SHI, m_output_cartesian_1d_Z_E_SHI
character(200) :: m_output_cartesian_1d_X_E_a, m_output_cartesian_1d_Y_E_a, m_output_cartesian_1d_Z_E_a
character(200) :: m_output_cartesian_1d_X_E_mu, m_output_cartesian_1d_Y_E_mu, m_output_cartesian_1d_Z_E_mu
! Velosity theta distribution:
character(200) :: m_output_velocity_theta_distr_ph, m_output_velocity_theta_distr_e, &
            m_output_velocity_theta_distr_h, m_output_velocity_theta_distr_p, m_output_velocity_theta_distr_SHI, &
            m_output_velocity_theta_distr_mu

!Cylindrical
!1d
character(200) :: m_output_cylindric_1d_R_ph, m_output_cylindric_1d_R_e, m_output_cylindric_1d_R_p
character(200) :: m_output_cylindric_1d_R_h, m_output_cylindric_1d_R_SHI, m_output_cylindric_1d_R_a, m_output_cylindric_1d_R_mu
character(200) :: m_output_cylindric_1d_R_E_ph, m_output_cylindric_1d_R_E_e, m_output_cylindric_1d_R_E_p
character(200) :: m_output_cylindric_1d_R_E_h, m_output_cylindric_1d_R_E_SHI, m_output_cylindric_1d_R_E_a, m_output_cylindric_1d_R_E_mu
!2d
character(200) :: m_output_cylindric_2d_RL_ph, m_output_cylindric_2d_RL_e, m_output_cylindric_2d_RL_p, m_output_cylindric_2d_RL_mu
character(200) :: m_output_cylindric_2d_RL_h, m_output_cylindric_2d_RL_SHI, m_output_cylindric_2d_RL_a
character(200) :: m_output_cylindric_2d_RL_E_ph, m_output_cylindric_2d_RL_E_e, m_output_cylindric_2d_RL_E_p, m_output_cylindric_2d_RL_E_mu
character(200) :: m_output_cylindric_2d_RL_E_h, m_output_cylindric_2d_RL_E_SHI, m_output_cylindric_2d_RL_E_a

! MD output
character(200) :: m_output_MD_energies
character(200) :: m_output_MD_cell_params
character(200) :: m_output_MD_coordinates
character(200) :: m_output_MD_velocities
character(200) :: m_output_MCMD
character(200) :: m_output_MD_LAMMPS
character(200) :: m_output_MD_displacements


! code version:
character(30), parameter :: m_TREKIS_version = 'TREKIS-4 (version 24.11.2025)'


! All output file names:
parameter (m_output_parameters = '!OUTPUT_parameters.txt')
parameter (m_communication = 'Communication.txt')
!dddddddddddddddddddddddddddddddddddddddddddd
parameter (m_output_DOS = 'OUTPUT_DOS_of_')
parameter (m_output_DOS_k = 'OUTPUT_DOS_k_vector_of_')
parameter (m_output_DOS_effm = 'OUTPUT_DOS_effective_mass_of_')
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
parameter (m_output_MD_energies = 'OUTPUT_MD_energies.txt')
parameter (m_output_MD_cell_params = 'OUTPUT_MD_average_parameters.txt')
parameter (m_output_MD_displacements = 'OUTPUT_MD_mean_displacements.txt')
parameter (m_output_MD_coordinates = 'OUTPUT_MD_coordinates.xyz')
parameter (m_output_MD_velocities = 'OUTPUT_MD_velocities.xyz')
parameter (m_output_MCMD = 'OUTPUT_MCMD_energy_transfer.txt')
parameter (m_output_MD_LAMMPS = 'OUTPUT_MD_LAMMPS.lmp')
!ssssssssssssssssssssssssssssssssssssssssssss
parameter (m_output_spectrum_ph = 'OUTPUT_photon_spectrum_')
parameter (m_output_spectrum_e = 'OUTPUT_electron_spectrum_')
parameter (m_output_spectrum_h = 'OUTPUT_valence_hole_spectrum_')
parameter (m_output_spectrum_p = 'OUTPUT_positron_spectrum_')
parameter (m_output_spectrum_SHI = 'OUTPUT_SHI_spectrum_')
parameter (m_output_spectrum_mu = 'OUTPUT_muon_spectrum_')
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
parameter (m_output_velocity_theta_distr_ph = 'OUTPUT_photon_velocity_theta_')
parameter (m_output_velocity_theta_distr_e = 'OUTPUT_electron_velocity_theta_')
parameter (m_output_velocity_theta_distr_h = 'OUTPUT_hole_velocity_theta_')
parameter (m_output_velocity_theta_distr_p = 'OUTPUT_positron_velocity_theta_')
parameter (m_output_velocity_theta_distr_SHI = 'OUTPUT_SHI_velocity_theta_')
parameter (m_output_velocity_theta_distr_mu = 'OUTPUT_muon_velocity_theta_')
!ssssssssssssssssssssssssssssssssssssssssssss
parameter (m_output_spectrum_ph_1d = 'OUTPUT_photon_spectrum_1d_')
parameter (m_output_spectrum_e_1d = 'OUTPUT_electron_spectrum_1d_')
parameter (m_output_spectrum_h_1d = 'OUTPUT_valence_hole_spectrum_1d_')
parameter (m_output_spectrum_p_1d = 'OUTPUT_positron_spectrum_1d_')
parameter (m_output_spectrum_mu_1d = 'OUTPUT_muon_spectrum_1d_')
parameter (m_output_spectrum_SHI_1d = 'OUTPUT_SHI_spectrum_1d_')
!ssssssssssssssssssssssssssssssssssssssssssss
parameter (m_output_theta_ph_1d = 'OUTPUT_photon_velocity_theta_distr_1d_')
parameter (m_output_theta_e_1d = 'OUTPUT_electron_velocity_theta_distr_1d_')
parameter (m_output_theta_h_1d = 'OUTPUT_valence_hole_velocity_theta_distr_1d_')
parameter (m_output_theta_p_1d = 'OUTPUT_positron_velocity_theta_distr_1d_')
parameter (m_output_theta_SHI_1d = 'OUTPUT_SHI_velocity_theta_distr_1d_')
parameter (m_output_theta_mu_1d = 'OUTPUT_muon_velocity_theta_distr_1d_')
!dddddddddddddddddddddddddddddddddddddddddddd
! Cartesian:
parameter (m_output_cartesian_1d_X_ph = 'OUTPUT_photon_density_1d_X_')
parameter (m_output_cartesian_1d_Y_ph = 'OUTPUT_photon_density_1d_Y_')
parameter (m_output_cartesian_1d_Z_ph = 'OUTPUT_photon_density_1d_Z_')
parameter (m_output_cartesian_1d_X_e = 'OUTPUT_electron_density_1d_X_')
parameter (m_output_cartesian_1d_Y_e = 'OUTPUT_electron_density_1d_Y_')
parameter (m_output_cartesian_1d_Z_e = 'OUTPUT_electron_density_1d_Z_')
parameter (m_output_cartesian_1d_X_p = 'OUTPUT_positron_density_1d_X_')
parameter (m_output_cartesian_1d_Y_p = 'OUTPUT_positron_density_1d_Y_')
parameter (m_output_cartesian_1d_Z_p = 'OUTPUT_positron_density_1d_Z_')
parameter (m_output_cartesian_1d_X_h = 'OUTPUT_holes_density_1d_X_')
parameter (m_output_cartesian_1d_Y_h = 'OUTPUT_holes_density_1d_Y_')
parameter (m_output_cartesian_1d_Z_h = 'OUTPUT_holes_density_1d_Z_')
parameter (m_output_cartesian_1d_X_SHI = 'OUTPUT_SHI_density_1d_X_')
parameter (m_output_cartesian_1d_Y_SHI = 'OUTPUT_SHI_density_1d_Y_')
parameter (m_output_cartesian_1d_Z_SHI = 'OUTPUT_SHI_density_1d_Z_')
parameter (m_output_cartesian_1d_X_a = 'OUTPUT_atomic_density_1d_X_')
parameter (m_output_cartesian_1d_Y_a = 'OUTPUT_atomic_density_1d_Y_')
parameter (m_output_cartesian_1d_Z_a = 'OUTPUT_atomic_density_1d_Z_')
parameter (m_output_cartesian_1d_X_mu = 'OUTPUT_muon_density_1d_X_')
parameter (m_output_cartesian_1d_Y_mu = 'OUTPUT_muon_density_1d_Y_')
parameter (m_output_cartesian_1d_Z_mu = 'OUTPUT_muon_density_1d_Z_')
parameter (m_output_cartesian_1d_X_E_ph = 'OUTPUT_photon_energy_1d_X_')
parameter (m_output_cartesian_1d_Y_E_ph = 'OUTPUT_photon_energy_1d_Y_')
parameter (m_output_cartesian_1d_Z_E_ph = 'OUTPUT_photon_energy_1d_Z_')
parameter (m_output_cartesian_1d_X_E_e = 'OUTPUT_electron_energy_1d_X_')
parameter (m_output_cartesian_1d_Y_E_e = 'OUTPUT_electron_energy_1d_Y_')
parameter (m_output_cartesian_1d_Z_E_e = 'OUTPUT_electron_energy_1d_Z_')
parameter (m_output_cartesian_1d_X_E_p = 'OUTPUT_positron_energy_1d_X_')
parameter (m_output_cartesian_1d_Y_E_p = 'OUTPUT_positron_energy_1d_Y_')
parameter (m_output_cartesian_1d_Z_E_p = 'OUTPUT_positron_energy_1d_Z_')
parameter (m_output_cartesian_1d_X_E_h = 'OUTPUT_holes_energy_1d_X_')
parameter (m_output_cartesian_1d_Y_E_h = 'OUTPUT_holes_energy_1d_Y_')
parameter (m_output_cartesian_1d_Z_E_h = 'OUTPUT_holes_energy_1d_Z_')
parameter (m_output_cartesian_1d_X_E_SHI = 'OUTPUT_SHI_energy_1d_X_')
parameter (m_output_cartesian_1d_Y_E_SHI = 'OUTPUT_SHI_energy_1d_Y_')
parameter (m_output_cartesian_1d_Z_E_SHI = 'OUTPUT_SHI_energy_1d_Z_')
parameter (m_output_cartesian_1d_X_E_a = 'OUTPUT_atomic_energy_1d_X_')
parameter (m_output_cartesian_1d_Y_E_a = 'OUTPUT_atomic_energy_1d_Y_')
parameter (m_output_cartesian_1d_Z_E_a = 'OUTPUT_atomic_energy_1d_Z_')
parameter (m_output_cartesian_1d_X_E_mu = 'OUTPUT_muon_energy_1d_X_')
parameter (m_output_cartesian_1d_Y_E_mu = 'OUTPUT_muon_energy_1d_Y_')
parameter (m_output_cartesian_1d_Z_E_mu = 'OUTPUT_muon_energy_1d_Z_')
!dddddddddddddddddddddddddddddddddddddddddddd
! Cyllindrical:
! 1-dimensional
parameter (m_output_cylindric_1d_R_ph = 'OUTPUT_photon_density_1d_R_')
parameter (m_output_cylindric_1d_R_e = 'OUTPUT_electron_density_1d_R_')
parameter (m_output_cylindric_1d_R_p = 'OUTPUT_positron_density_1d_R_')
parameter (m_output_cylindric_1d_R_h = 'OUTPUT_holes_density_1d_R_')
parameter (m_output_cylindric_1d_R_SHI = 'OUTPUT_SHI_density_1d_R_')
parameter (m_output_cylindric_1d_R_a = 'OUTPUT_atomic_density_1d_R_')
parameter (m_output_cylindric_1d_R_mu = 'OUTPUT_muon_density_1d_R_')
parameter (m_output_cylindric_1d_R_E_ph = 'OUTPUT_photon_energy_1d_R_')
parameter (m_output_cylindric_1d_R_E_e = 'OUTPUT_electron_energy_1d_R_')
parameter (m_output_cylindric_1d_R_E_p = 'OUTPUT_positron_energy_1d_R_')
parameter (m_output_cylindric_1d_R_E_h = 'OUTPUT_holes_energy_1d_R_')
parameter (m_output_cylindric_1d_R_E_SHI = 'OUTPUT_SHI_energy_1d_R_')
parameter (m_output_cylindric_1d_R_E_a = 'OUTPUT_atomic_energy_1d_R_')
parameter (m_output_cylindric_1d_R_E_mu = 'OUTPUT_muon_energy_1d_R_')
!2-dimensional
parameter (m_output_cylindric_2d_RL_ph = 'OUTPUT_photon_density_2d_RL_')
parameter (m_output_cylindric_2d_RL_e = 'OUTPUT_electron_density_2d_RL_')
parameter (m_output_cylindric_2d_RL_p = 'OUTPUT_positron_density_2d_RL_')
parameter (m_output_cylindric_2d_RL_h = 'OUTPUT_holes_density_2d_RL_')
parameter (m_output_cylindric_2d_RL_SHI = 'OUTPUT_SHI_density_2d_RL_')
parameter (m_output_cylindric_2d_RL_a = 'OUTPUT_atomic_density_2d_RL_')
parameter (m_output_cylindric_2d_RL_mu = 'OUTPUT_muon_density_2d_RL_')
parameter (m_output_cylindric_2d_RL_E_ph = 'OUTPUT_photon_energy_2d_RL_')
parameter (m_output_cylindric_2d_RL_E_e = 'OUTPUT_electron_energy_2d_RL_')
parameter (m_output_cylindric_2d_RL_E_p = 'OUTPUT_positron_energy_2d_RL_')
parameter (m_output_cylindric_2d_RL_E_h = 'OUTPUT_holes_energy_2d_RL_')
parameter (m_output_cylindric_2d_RL_E_SHI = 'OUTPUT_SHI_energy_2d_RL_')
parameter (m_output_cylindric_2d_RL_E_a = 'OUTPUT_atomic_energy_2d_RL_')
parameter (m_output_cylindric_2d_RL_E_mu = 'OUTPUT_muon_energy_2d_RL_')
!ttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
parameter (m_output_total = 'OUTPUT_total_')
parameter (m_output_total_cutoff = 'OUTPUT_total_above_cutoff_')
parameter (m_output_N_gnu = 'OUTPUT_total_numbers_')
parameter (m_output_E_gnu = 'OUTPUT_total_energies_')
parameter (m_output_MD = 'OUTPUT_MD_')
parameter (m_output_MD_E_gnu = 'OUTPUT_MD_energies_')
parameter (m_output_MD_T_gnu = 'OUTPUT_MD_temperature_')
parameter (m_output_MD_MSD_gnu = 'OUTPUT_MD_displacements_')
!rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
parameter (m_folder_MFP = 'MFPs_and_Ranges_in_')
parameter (m_output_Compton = 'OUTPUT_Compton_')
parameter (m_output_Rayleigh = 'OUTPUT_Rayleigh_')
parameter (m_output_pair = 'OUTPUT_pair_')
parameter (m_output_absorb = 'OUTPUT_absorption_')
parameter (m_output_MFP = 'OUTPUT_MFPs_')
parameter (m_output_IMFP = 'OUTPUT_IMFPs_')
parameter (m_output_EMFP = 'OUTPUT_EMFPs_')
parameter (m_output_Brems = 'OUTPUT_Brems_MFPs_')
parameter (m_output_annihil = 'OUTPUT_Annihilation_MFPs_')
parameter (m_output_Range = 'OUTPUT_Ranges_')
parameter (m_output_Se = 'OUTPUT_Stopping_')
parameter (m_output_Se_vs_range = 'OUTPUT_Se_vs_range_')


 contains

!===============================================
! Write out all output files:
subroutine write_output_files(used_target, numpar, out_data, MD_atoms, MD_supce, MD_pots, tim, laststep)
    type(Matter), intent(in) :: used_target      ! parameters of the target
    type(Num_par), intent(inout) :: numpar    ! all numerical parameters
    type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
    type(Atom), dimension(:), intent(inout), allocatable :: MD_atoms ! all atoms in MD as objects
    type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
    type(MD_potential), dimension(:,:), intent(in), allocatable :: MD_pots    ! MD potentials for each kind of atom-atom interactions
    real(8), intent(in) :: tim  ! simulation time step [fs]
    logical, intent(in), optional :: laststep   ! if it is last time step of the simulation
    !---------------------------------
    logical :: last_step
    
    if (present(laststep)) then ! set a marker if it is last time step
       last_step = laststep
    else
       last_step = .false.
    endif

    ! Printout total values:
    call total_values(used_target, numpar, out_data, tim)   ! below
    
    ! Printout spectra:
    call energy_spectra(used_target, numpar, out_data, tim)   ! below

    ! Printout spatial distributions:
    call spatial_distributions(used_target, numpar, out_data, tim)   ! below
    
    ! Printout velocity distributions by theta:
    call velocity_theta_distributions(used_target, numpar, out_data, tim)   ! below
    
    ! Printout MD data:
    call MD_data_printout(used_target, numpar, out_data, MD_atoms, MD_supce, MD_pots, tim, last_step)  ! below
end subroutine write_output_files



subroutine MD_data_printout(used_target, numpar, out_data, MD_atoms, MD_supce, MD_pots, tim, last_step)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout) :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   type(Atom), dimension(:), intent(inout), allocatable :: MD_atoms ! all atoms in MD as objects
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   type(MD_potential), dimension(:,:), intent(in), allocatable :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   real(8), intent(in) :: tim  ! simulation time step [fs]
   logical, intent(in), optional :: last_step ! some options (like LAMMPS) are printed only on last step
   !--------------------------
   logical :: file_exist, its_last_step
   character(20) :: text_ch
   
   if (numpar%DO_MD) then   ! only if MD module is active it makes sense to analyze the data
      ! MD energy ballance:
      call MD_total_values(numpar, MD_supce, MD_pots, out_data, tim)    ! below
      
      ! Atomic coordinates:
      if (numpar%print_MD_R_xyz) then   ! if user requested atomic coordinates:
         ! Create a file, if does not exist yet:
         numpar%FILE_MD_XYZ = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_MD_coordinates))
         inquire(file=trim(adjustl(numpar%FILE_MD_XYZ)),exist=file_exist) ! check if input file is there
         if (.not. file_exist) then   ! it's the first time, create file and write the header
            open(newunit = numpar%FN_MD_XYZ, FILE = trim(adjustl(numpar%FILE_MD_XYZ)))
         endif
         write(text_ch,'(f16.2)') tim
         call write_XYZ(numpar%FN_MD_XYZ, MD_atoms, 'Time '//trim(adjustl(text_ch)))  ! module "Dealing_with_XYZ_files"
      endif ! (numpar%print_MD_R_xyz)
      
      ! Atomic velocities:
      if (numpar%print_MD_V_xyz) then   ! if user requested atomic velocities:
         ! Create a file, if does not exist yet:
         numpar%FILE_MD_V_XYZ = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_MD_velocities))
         inquire(file=trim(adjustl(numpar%FILE_MD_V_XYZ)),exist=file_exist) ! check if input file is there
         if (.not. file_exist) then   ! it's the first time, create file and write the header
            open(newunit = numpar%FN_MD_V_XYZ, FILE = trim(adjustl(numpar%FILE_MD_V_XYZ)))
         endif
         write(text_ch,'(f16.2)') tim
         call write_XYZ(numpar%FN_MD_V_XYZ, MD_atoms, 'Time '//trim(adjustl(text_ch)), print_velocities = .true.)  ! module "Dealing_with_XYZ_files"
      endif ! (numpar%print_MD_V_xyz)

      ! LAMMPS input file:
      if (present(last_step)) then
         its_last_step = last_step  ! user defined
      else
         its_last_step = .false.    ! default
      endif
      ! LAMMPS file is only created at the last timestep:
      if (numpar%print_MD_LAMMPS .and. its_last_step) then   ! if user requested LAMMPS input file:
         ! Create a file, if does not exist yet:
         numpar%FILE_MD_LAMMPS = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_MD_LAMMPS))
         inquire(file=trim(adjustl(numpar%FILE_MD_LAMMPS)),exist=file_exist) ! check if input file is there
         if (.not. file_exist) then   ! it's the first time, create file and write the header
            open(newunit = numpar%FN_MD_LAMMPS, FILE = trim(adjustl(numpar%FILE_MD_LAMMPS)))
         endif
         write(text_ch,'(f16.2)') tim
         call Write_LAMMPS_input_file(numpar%FN_MD_LAMMPS, used_target, numpar, MD_supce, MD_atoms, MD_pots, &
            'TREKIS-4 : '//trim(adjustl(numpar%output_path))//' at '//trim(adjustl(text_ch))//' fs' )  ! module "Dealing_with_LAMMPS"
      endif ! (numpar%print_MD_LAMMPS)
   endif

   ! Energy transferred from MC to MD:
   if (numpar%print_MC_MD_energy) then   ! if user requested
      call write_MCMD(numpar, MD_supce, tim)  ! below
   endif ! (numpar%print_MC_MD_energy)

end subroutine MD_data_printout


subroutine write_MCMD(numpar, MD_supce, tim)
   type(Num_par), intent(inout) :: numpar    ! all numerical parameters
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------
   logical :: file_exist
   integer :: i, j, k, Nsiz, Nsiz2, Nsiz3
   character(10) :: axis_name(3)

   ! Create a file with energies, if does not exist yet:
   numpar%FILE_MCMD_info = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_MCMD))
   inquire(file=trim(adjustl(numpar%FILE_MCMD_info)),exist=file_exist) ! check if input file is there

   if (.not. file_exist) then   ! it's the first time, create file and write the header
      open(newunit = numpar%FN_MCMD_info, FILE = trim(adjustl(numpar%FILE_MCMD_info)))
      select case (MD_supce%coord_type)
      case (1)   ! Cartesian
         do i = 1,3 ! set names of the coordinate axes
            select case (MD_supce%coord_axis(i))
            case (1)
               write(axis_name(i),'(a)') 'X'
            case (2)
               write(axis_name(i),'(a)') 'Y'
            case (3)
               write(axis_name(i),'(a)') 'Z'
            case default
               write(axis_name(i),'(a)') 'O'
            endselect
         enddo
         write(numpar%FN_MCMD_info,'(a)') '#Time    '//trim(adjustl(axis_name(1)))//'   '//&
                trim(adjustl(axis_name(2)))//'  '//trim(adjustl(axis_name(3)))//'   '//&
                'el_elast    h_elast     e_cut   h_cut  tot'
         write(numpar%FN_MCMD_info,'(a)') '#fs    A   A A eV    eV  eV  eV  eV  eV'
      case (2)   ! Spherical
         do i = 1,3 ! set names of the coordinate axes
            select case (MD_supce%coord_axis(i))
            case (1)
               write(axis_name(i),'(a)') 'R'
            case (2)
               write(axis_name(i),'(a)') 'Theta'
            case (3)
               write(axis_name(i),'(a)') 'Phi'
            case default
               write(axis_name(i),'(a)') 'O'
            endselect
         enddo
         write(numpar%FN_MCMD_info,'(a)') '#Time    '//trim(adjustl(axis_name(1)))//'   '//&
                trim(adjustl(axis_name(2)))//'  '//trim(adjustl(axis_name(3)))//'   '//&
                'el_elast    h_elast     e_cut   h_cut  tot'
         write(numpar%FN_MCMD_info,'(a)') '#fs    A   deg deg eV    eV  eV  eV  eV  eV'
      case (3)   ! Cylindrical
         do i = 1,3 ! set names of the coordinate axes
            select case (MD_supce%coord_axis(i))
            case (1)
               write(axis_name(i),'(a)') 'R'
            case (2)
               write(axis_name(i),'(a)') 'L'
            case (3)
               write(axis_name(i),'(a)') 'Theta'
            case default
               write(axis_name(i),'(a)') 'O'
            endselect
         enddo
         write(numpar%FN_MCMD_info,'(a)') '#Time    '//trim(adjustl(axis_name(1)))//'   '//&
                trim(adjustl(axis_name(2)))//'  '//trim(adjustl(axis_name(3)))//'   '//&
                'el_elast    h_elast     e_cut   h_cut  tot'
         write(numpar%FN_MCMD_info,'(a)') '#fs    A   A deg eV    eV  eV  eV  eV    eV'
      end select
   endif

   ! Get the energy transfer from MC to MD module to print out:
   Nsiz = size(MD_supce%MC2MD_dimen(1)%grid)
   Nsiz2 = size(MD_supce%MC2MD_dimen(2)%grid)
   Nsiz3 = size(MD_supce%MC2MD_dimen(3)%grid)

   do k = 1, Nsiz3
      do j = 1, Nsiz2
         do i = 1, Nsiz
             write(numpar%FN_MCMD_info,'(f,f,f,f,es,es,es,es,es)') tim, &
             MD_supce%MC2MD_dimen(1)%grid(i), MD_supce%MC2MD_dimen(2)%grid(j), MD_supce%MC2MD_dimen(3)%grid(k), &
             MD_supce%E_e_at_from_MC(i,j,k), MD_supce%E_h_at_from_MC(i,j,k), &
             MD_supce%E_e_from_MC(i,j,k), MD_supce%E_h_from_MC(i,j,k), &
             MD_supce%E_e_at_from_MC(i,j,k) + MD_supce%E_h_at_from_MC(i,j,k) + &
             MD_supce%E_e_from_MC(i,j,k) + MD_supce%E_h_from_MC(i,j,k)
         enddo ! i
      enddo ! j
   enddo ! k

   if (MD_supce%coord_dim > 0) then ! no need to skip lines for zero dimensional data
      write(numpar%FN_MCMD_info,'(a)') '   '   ! an empty line to separate the time-blocks
   endif
end subroutine write_MCMD


subroutine MD_total_values(numpar, MD_supce, MD_pots, out_data, tim)
   type(Num_par), intent(inout) :: numpar    ! all numerical parameters
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   type(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(output_data), intent(in) :: out_data ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !---------------------------------
   real(8) :: Etot, Ptot
   integer :: N_KOA, i
   logical :: file_exist
   character(10) :: chtemp

   ! Create a file with energies, if does not exist yet:
   numpar%FILE_MD_totals = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_MD_energies))
   inquire(file=trim(adjustl(numpar%FILE_MD_totals)),exist=file_exist) ! check if input file is there
   
   if (.not. file_exist) then   ! it's the first time, create file and write the header
      open(newunit = numpar%FN_MD_totals, FILE = trim(adjustl(numpar%FILE_MD_totals)))
      write(numpar%FN_MD_totals,'(a)') '#Time    E_kin  E_pot  E_tot '
      write(numpar%FN_MD_totals,'(a)') '#fs    eV/atom   eV/atom eV/atom '
   endif
   ! Get the total energy in the MD module to print out:
   Etot = out_data%MD_Ekin + out_data%MD_Epot
   write(numpar%FN_MD_totals,'(es24.16, $)') tim, out_data%MD_Ekin, out_data%MD_Epot, Etot
   write(numpar%FN_MD_totals,'(a)') ''
   
   
   ! Create a file with other parameters (average atomic temperature and pressure)
   numpar%FILE_MD_average = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_MD_cell_params))
   inquire(file=trim(adjustl(numpar%FILE_MD_average)),exist=file_exist) ! check if input file is there
   
   if (.not. file_exist) then   ! it's the first time, create file and write the header
      open(newunit = numpar%FN_MD_average, FILE = trim(adjustl(numpar%FILE_MD_average)))
      write(numpar%FN_MD_average,'(a)') '#Time    Temperature  Pressure  Pxx  Pxy    Pxz Pyx Pyy Pyz Pzx Pzy Pzz'
      write(numpar%FN_MD_average,'(a)') '#fs    K  GPa   GPa   GPa GPa GPa GPa GPa GPa GPa GPa'
   endif
   ! Get the total pressure in the MD module to print out:
   Ptot = 1.0d0/3.0d0 * (MD_supce%Pressure(1,1) + MD_supce%Pressure(2,2) + MD_supce%Pressure(3,3))
   write(numpar%FN_MD_average,'(es24.16, $)') tim, MD_supce%T_kin, Ptot, MD_supce%Pressure(1:3,1:3)
   write(numpar%FN_MD_average,'(a)') ''


   ! Number of different kinds of atoms (elements) in MD:
   N_KOA = size(out_data%MD_MSDP)

   ! Create a file with mean atomic displacements:
   numpar%FILE_MD_displacement = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_MD_displacements))
   inquire(file=trim(adjustl(numpar%FILE_MD_displacement)),exist=file_exist) ! check if input file is there
   if (.not. file_exist) then   ! it's the first time, create file and write the header
      open(newunit = numpar%FN_MD_displacement, FILE = trim(adjustl(numpar%FILE_MD_displacement)))
      if (numpar%n_MSD /= 1) then
         write(chtemp,'(i2)') numpar%n_MSD
         if (N_KOA > 1) then  ! many kinds of atoms
            write(numpar%FN_MD_displacement,'(a)', ADVANCE='no') '#Time    Mean_displacement  Displacement_'//trim(adjustl(MD_pots(1,1)%El1))//'  '
            do i = 2, N_KOA
               write(numpar%FN_MD_displacement,'(a)', ADVANCE='no') 'Displacement_'//trim(adjustl(MD_pots(i,i)%El1))//' '
            enddo
            write(numpar%FN_MD_displacement,'(a)') ''
            write(numpar%FN_MD_displacement,'(a)', ADVANCE='no') '#fs    A^'//trim(adjustl(chtemp))//' A^'//trim(adjustl(chtemp))//'   '
            do i = 2, N_KOA
               write(numpar%FN_MD_displacement,'(a)', ADVANCE='no') 'A^'//trim(adjustl(chtemp))//' '
            enddo
            write(numpar%FN_MD_displacement,'(a)') ''
         else  ! only one kind of atoms
            write(numpar%FN_MD_displacement,'(a)') '#Time    Mean_displacement'
            write(numpar%FN_MD_displacement,'(a)') '#fs    A^'//trim(adjustl(chtemp))
         endif
      else
         if (N_KOA > 1) then  ! many kinds of atoms
            write(numpar%FN_MD_displacement,'(a)', ADVANCE='no') '#Time    Mean_displacement  Displacement_'//trim(adjustl(MD_pots(1,1)%El1))//'  '
            do i = 2, N_KOA
               write(numpar%FN_MD_displacement,'(a)', ADVANCE='no') 'Displacement_'//trim(adjustl(MD_pots(i,i)%El1))//' '
            enddo
            write(numpar%FN_MD_displacement,'(a)') ''
            write(numpar%FN_MD_displacement,'(a)', ADVANCE='no') '#fs    A A  '
            do i = 2, N_KOA
               write(numpar%FN_MD_displacement,'(a)', ADVANCE='no') 'A  '
            enddo
            write(numpar%FN_MD_displacement,'(a)') ''
         else  ! only one kind of atoms
            write(numpar%FN_MD_displacement,'(a)') '#Time    Mean_displacement'
            write(numpar%FN_MD_displacement,'(a)') '#fs    A'
         endif
      endif
   endif
   ! Save the data into the file:
   if (N_KOA > 1) then ! many kinds of atoms
      write(numpar%FN_MD_displacement,'(es24.16, $)') tim, out_data%MD_MSD, out_data%MD_MSDP(:)
      write(numpar%FN_MD_displacement,'(a)') ''
   else ! only one kind of atoms
      write(numpar%FN_MD_displacement,'(es24.16, es24.16)') tim, out_data%MD_MSD
   endif



end subroutine MD_total_values





subroutine velocity_theta_distributions(used_target, numpar, out_data, tim) 
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout) :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !---------------------------------
   ! printout theta=acos(Vz/V):
   if (numpar%vel_theta_grid_par%along_axis) then
      ! Photons:
      call photon_velotheta(used_target, numpar, out_data, tim) ! below
      ! Electrons:
      call electron_velotheta(used_target, numpar, out_data, tim) ! below
      ! Positrons:
      call positron_velotheta(used_target, numpar, out_data, tim) ! below
      ! Holes:
      call hole_velotheta(used_target, numpar, out_data, tim) ! below
      ! SHIs:
      call SHI_velotheta(used_target, numpar, out_data, tim) ! below
      ! Muons:
      call muon_velotheta(used_target, numpar, out_data, tim) ! below
   endif

   if (ANY(numpar%Theta_grid_par(:)%along_axis)) then  ! Along axis
      call electron_theta_distr_1d(used_target, numpar, out_data, tim)  ! below
   endif
end subroutine velocity_theta_distributions



subroutine spatial_distributions(used_target, numpar, out_data, tim) 
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout) :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !---------------------------------
   
   !1111111111111111111111
   ! Spatial distributions in 1d:
   
   ! Cartesian:
      call printout_cartesian_1d(used_target, numpar, out_data, tim) ! below
   ! Cylindrical:          
      call printout_cylindric_1d(used_target, numpar, out_data, tim) ! below
   ! Spherical:
        ! NOT READY YET 
   
   !2222222222222222222222
   ! Spatial distributions in 2d:
   
   ! Cartesian:
        ! NOT READY YET
   ! Cylindrical:
      call printout_cylindric_2d(used_target, numpar, out_data, tim) ! below
   ! Spherical:
        ! NOT READY YET
   
   !3333333333333333333333
   ! Spatial distributions in 3d:

   ! Cartesian:
        ! NOT READY YET
   ! Cylindrical:
        ! NOT READY YET
   ! Spherical:
        ! NOT READY YET
   
end subroutine spatial_distributions




subroutine printout_cartesian_1d(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
    ! Along radius X:
   if (numpar%grid_par(1)%along_axis) then

      call printout_cartesian_1d_particle(numpar%FN_car_1d_X_ph, numpar%FN_car_1d_X_E_ph, m_output_cartesian_1d_X_ph, m_output_cartesian_1d_X_E_ph, &
                    used_target, numpar, out_data%Distr_ph_X(:),  out_data%E_Distr_ph_X(:), tim, 1, 1)  ! photons; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_X_e, numpar%FN_car_1d_X_E_e, m_output_cartesian_1d_X_e, m_output_cartesian_1d_X_E_e, &
                    used_target, numpar, out_data%Distr_e_X(:),  out_data%E_Distr_e_X(:), tim, 1, 1)  ! electrons; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_X_p, numpar%FN_car_1d_X_E_p, m_output_cartesian_1d_X_p, m_output_cartesian_1d_X_E_p, &
                    used_target, numpar, out_data%Distr_p_X(:),  out_data%E_Distr_p_X(:), tim, 1, 1)  ! positrons; below
      call printout_cartesian_1d_hole(used_target, numpar, out_data, tim, 1)  ! all holes; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_X_SHI, numpar%FN_car_1d_X_E_SHI, &
                    m_output_cartesian_1d_X_SHI, m_output_cartesian_1d_X_E_SHI, &
                    used_target, numpar, out_data%Distr_SHI_X(:),  out_data%E_Distr_SHI_X(:), tim, 1, 1)  ! shi; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_X_a, numpar%FN_car_1d_X_E_a, m_output_cartesian_1d_X_a, m_output_cartesian_1d_X_E_a, &
                    used_target, numpar, out_data%Distr_a_X(:),  out_data%E_Distr_a_X(:), tim, 1, 1)  ! atomic; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_X_mu, numpar%FN_car_1d_X_E_mu, m_output_cartesian_1d_X_mu, &
                    m_output_cartesian_1d_X_E_mu, used_target, numpar, out_data%Distr_mu_X(:),  out_data%E_Distr_mu_X(:), tim, 1, 1)  ! muons; below
   endif

   ! Along depth Y:
   if (numpar%grid_par(2)%along_axis) then
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Y_ph, numpar%FN_car_1d_Y_E_ph, m_output_cartesian_1d_Y_ph, m_output_cartesian_1d_Y_E_ph, &
                    used_target, numpar, out_data%Distr_ph_Y(:),  out_data%E_Distr_ph_Y(:), tim, 2, 2)  ! photons; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Y_e, numpar%FN_car_1d_Y_E_e, m_output_cartesian_1d_Y_e, m_output_cartesian_1d_Y_E_e, &
                    used_target, numpar, out_data%Distr_e_Y(:),  out_data%E_Distr_e_Y(:), tim, 2, 2)  ! electrons; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Y_p, numpar%FN_car_1d_Y_E_p, m_output_cartesian_1d_Y_p, m_output_cartesian_1d_Y_E_p, &
                    used_target, numpar, out_data%Distr_p_Y(:),  out_data%E_Distr_p_Y(:), tim, 2, 2)  ! positrons; below
      call printout_cartesian_1d_hole(used_target, numpar, out_data, tim, 2)  ! all holes; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Y_SHI, numpar%FN_car_1d_Y_E_SHI, m_output_cartesian_1d_Y_SHI, m_output_cartesian_1d_Y_E_SHI, &
                    used_target, numpar, out_data%Distr_SHI_Y(:),  out_data%E_Distr_SHI_Y(:), tim, 2, 2)  ! shi; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Y_a, numpar%FN_car_1d_Y_E_a, m_output_cartesian_1d_Y_a, m_output_cartesian_1d_Y_E_a, &
                    used_target, numpar, out_data%Distr_a_Y(:),  out_data%E_Distr_a_Y(:), tim, 2, 2)  ! atomic; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Y_mu, numpar%FN_car_1d_Y_E_mu, m_output_cartesian_1d_Y_mu, &
                    m_output_cartesian_1d_Y_E_mu, &
                    used_target, numpar, out_data%Distr_mu_Y(:),  out_data%E_Distr_mu_Y(:), tim, 2, 2)  ! muon; below
   endif

   ! Along Z:
   if (numpar%grid_par(3)%along_axis) then
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Z_ph, numpar%FN_car_1d_Z_E_ph, m_output_cartesian_1d_Z_ph, m_output_cartesian_1d_Z_E_ph, &
                    used_target, numpar, out_data%Distr_ph_Z(:),  out_data%E_Distr_ph_Z(:), tim, 3, 3)  ! photons; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Z_e, numpar%FN_car_1d_Z_E_e, m_output_cartesian_1d_Z_e, m_output_cartesian_1d_Z_E_e, &
                    used_target, numpar, out_data%Distr_e_Z(:),  out_data%E_Distr_e_Z(:), tim, 3, 3)  ! electrons; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Z_p, numpar%FN_car_1d_Z_E_p, m_output_cartesian_1d_Z_p, m_output_cartesian_1d_Z_E_p, &
                    used_target, numpar, out_data%Distr_p_Z(:),  out_data%E_Distr_p_Z(:), tim, 3, 3)  ! positrons; below
      call printout_cartesian_1d_hole(used_target, numpar, out_data, tim, 3)  ! all holes; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Z_SHI, numpar%FN_car_1d_Z_E_SHI, m_output_cartesian_1d_Z_SHI, m_output_cartesian_1d_Z_E_SHI, &
                    used_target, numpar, out_data%Distr_SHI_Z(:),  out_data%E_Distr_SHI_Z(:), tim, 3, 3)  ! shi; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Z_a, numpar%FN_car_1d_Z_E_a, m_output_cartesian_1d_Z_a, m_output_cartesian_1d_Z_E_a, &
                    used_target, numpar, out_data%Distr_a_Z(:),  out_data%E_Distr_a_Z(:), tim, 3, 3)  ! atomic; below
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Z_mu, numpar%FN_car_1d_Z_E_mu, m_output_cartesian_1d_Z_mu, &
                    m_output_cartesian_1d_Z_E_mu, &
                    used_target, numpar, out_data%Distr_mu_Z(:),  out_data%E_Distr_mu_Z(:), tim, 3, 3)  ! muons; below
   endif
end subroutine printout_cartesian_1d




subroutine printout_cartesian_1d_particle(FN, FN_E, Part_file_name_D, Part_file_name_E, used_target, numpar, Distr_coord, Distr_E, tim, grid_ind, axis_ind)  ! particles
   integer, intent(inout) :: FN, FN_E    ! file number to save data into, for density and energy
   character(*), intent(in) :: Part_file_name_D, Part_file_name_E   ! part of the name of the output file for a particle density and energy density
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(in), target :: numpar    ! all numerical parameters
   real(8), dimension(:), intent(in) :: Distr_coord, Distr_E  ! output distributions
   real(8), intent(in) :: tim  ! simulation time step [fs]
   integer, intent(in) :: grid_ind  ! index, which grid to use
   integer, intent(in) :: axis_ind  ! 1=X, 2=Y, 3=Z
   !----------------------------------------
   integer :: i_tar
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(Part_file_name_D))//trim(adjustl(used_target%Name))//'.dat'
                    !trim(adjustl(m_output_cylindric_1d_R_ph))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = FN, FILE = trim(adjustl(File_name)))
         select case (axis_ind)
         case (1)   ! X
            write(FN,'(a)') '#Time    X  Density'
         case (2)   ! Y
            write(FN,'(a)') '#Time    Y  Density'
         case (3)   ! Z
            write(FN,'(a)') '#Time    Z  Density'
         end select
         write(FN,'(a)') '#fs    A   1/A^3'
      endif
!       FN = numpar%FN_car_1d_R_ph ! just set a number

      ! Write data with densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(grid_ind)%spatial_grid1(:), Distr_coord(:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(Part_file_name_E))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_E_ph))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = FN_E, FILE = trim(adjustl(File_name)))
!          write(numpar%FN_car_1d_R_E_ph,'(a)') '#Time    R    Energy'
         select case (axis_ind)
         case (1)   ! X
            write(FN_E,'(a)') '#Time    X  Energy'
         case (2)   ! Y
            write(FN_E,'(a)') '#Time    Y  Energy'
         case (3)   ! Z
            write(FN_E,'(a)') '#Time    Z  Energy'
         end select
         write(FN_E,'(a)') '#fs    A   eV/A^3'
      endif

      ! Write data with energy densities into this file:
      call printout_data_on_1d_grid(FN_E, numpar%grids(grid_ind)%spatial_grid1(:), Distr_E(:), tim)    ! below
   enddo TRGT
end subroutine printout_cartesian_1d_particle



subroutine printout_cartesian_1d_hole(used_target, numpar, out_data, tim, x_ind)  ! valence holes
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   integer, intent(in) :: x_ind ! 1=X, 2=Y, 3=Z
   !----------------------------------------
   integer :: i_tar, i_KOA, i_NSH, i_arr
   character(50) :: File_name
   
   ! Make sure the arrays with files numbers for holes are allocated:
   select case (x_ind)
   case (1) ! X
      if (.not.allocated(numpar%FN_car_1d_X_h)) allocate(numpar%FN_car_1d_X_h(numpar%N_sh_tot))
      if (.not.allocated(numpar%FN_car_1d_X_E_h)) allocate(numpar%FN_car_1d_X_E_h(numpar%N_sh_tot))
      ! Printout hole distribution in each target:
      i_arr = 1    ! do valence band first
      call printout_cartesian_1d_particle(numpar%FN_car_1d_X_h(i_arr), numpar%FN_car_1d_X_E_h(i_arr), trim(adjustl(m_output_cartesian_1d_X_h))//'Valence_', &
                trim(adjustl(m_output_cartesian_1d_X_E_h))//'Valence_', used_target, numpar, out_data%Distr_h_X(i_arr,:),  out_data%E_Distr_h_X(i_arr,:), tim, 1, i_arr)  ! holes; below

      TRGT1:do i_tar = 1, used_target%NOC ! for all targets
         KOA1:do i_KOA = 1, size(used_target%Material(i_tar)%Elements)    ! for all elements within this target
            SH1:do i_NSH = 1, used_target%Material(i_tar)%Elements(i_KOA)%N_core_shl    ! for all core shells of this element
               i_arr = i_arr + 1    ! next array index
               write(File_name,'(a)') trim(adjustl(used_target%Material(i_tar)%Elements(i_KOA)%Name))//'_'// &
                                            trim(adjustl(used_target%Material(i_tar)%Elements(i_KOA)%Shell_name(i_NSH)))
!                call printout_cylindric_1d_h_single_shell(trim(adjustl(File_name)), used_target, numpar, out_data, i_arr, tim)  ! below
               call printout_cartesian_1d_particle(numpar%FN_car_1d_X_h(i_arr), numpar%FN_car_1d_X_E_h(i_arr), &
                trim(adjustl(m_output_cartesian_1d_X_h))//trim(adjustl(File_name))//'_', trim(adjustl(m_output_cartesian_1d_X_E_h))//trim(adjustl(File_name))//'_', &
                used_target, numpar, out_data%Distr_h_X(i_arr,:),  out_data%E_Distr_h_X(i_arr,:), tim, 1, i_arr)  ! holes; below
            enddo SH1
         enddo KOA1
      enddo TRGT1
      
   case (2) ! Y
      if (.not.allocated(numpar%FN_car_1d_Y_h)) allocate(numpar%FN_car_1d_Y_h(numpar%N_sh_tot))
      if (.not.allocated(numpar%FN_car_1d_Y_E_h)) allocate(numpar%FN_car_1d_Y_E_h(numpar%N_sh_tot))
      ! Printout hole distribution in each target:
      i_arr = 1    ! do valence band first
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Y_h(i_arr), numpar%FN_car_1d_Y_E_h(i_arr), trim(adjustl(m_output_cartesian_1d_Y_h))//'Valence_', &
                trim(adjustl(m_output_cartesian_1d_Y_E_h))//'Valence_', used_target, numpar, out_data%Distr_h_Y(i_arr,:),  out_data%E_Distr_h_Y(i_arr,:), tim, 2, i_arr)  ! holes; below
   
      TRGT2:do i_tar = 1, used_target%NOC ! for all targets
         KOA2:do i_KOA = 1, size(used_target%Material(i_tar)%Elements)    ! for all elements within this target
            SH2:do i_NSH = 1, used_target%Material(i_tar)%Elements(i_KOA)%N_core_shl    ! for all core shells of this element
               i_arr = i_arr + 1    ! next array index
               write(File_name,'(a)') trim(adjustl(used_target%Material(i_tar)%Elements(i_KOA)%Name))//'_'// &
                                            trim(adjustl(used_target%Material(i_tar)%Elements(i_KOA)%Shell_name(i_NSH)))
!                call printout_cylindric_1d_h_single_shell(trim(adjustl(File_name)), used_target, numpar, out_data, i_arr, tim)  ! below
               call printout_cartesian_1d_particle(numpar%FN_car_1d_Y_h(i_arr), numpar%FN_car_1d_Y_E_h(i_arr), &
                trim(adjustl(m_output_cartesian_1d_Y_h))//trim(adjustl(File_name))//'_', trim(adjustl(m_output_cartesian_1d_Y_E_h))//trim(adjustl(File_name))//'_', &
                used_target, numpar, out_data%Distr_h_Y(i_arr,:),  out_data%E_Distr_h_Y(i_arr,:), tim, 2, i_arr)  ! holes; below
            enddo SH2
         enddo KOA2
      enddo TRGT2
      
   case (3) ! Z
      if (.not.allocated(numpar%FN_car_1d_Z_h)) allocate(numpar%FN_car_1d_Z_h(numpar%N_sh_tot))
      if (.not.allocated(numpar%FN_car_1d_Z_E_h)) allocate(numpar%FN_car_1d_Z_E_h(numpar%N_sh_tot))
      ! Printout hole distribution in each target:
      i_arr = 1    ! do valence band first
      call printout_cartesian_1d_particle(numpar%FN_car_1d_Z_h(i_arr), numpar%FN_car_1d_Z_E_h(i_arr), trim(adjustl(m_output_cartesian_1d_Z_h))//'Valence_', &
                trim(adjustl(m_output_cartesian_1d_Z_E_h))//'Valence_', used_target, numpar, out_data%Distr_h_Z(i_arr,:),  out_data%E_Distr_h_Z(i_arr,:), tim, 3, i_arr)  ! holes; below
   
      TRGT3:do i_tar = 1, used_target%NOC ! for all targets
         KOA3:do i_KOA = 1, size(used_target%Material(i_tar)%Elements)    ! for all elements within this target
            SH3:do i_NSH = 1, used_target%Material(i_tar)%Elements(i_KOA)%N_core_shl    ! for all core shells of this element
               i_arr = i_arr + 1    ! next array index
               write(File_name,'(a)') trim(adjustl(used_target%Material(i_tar)%Elements(i_KOA)%Name))//'_'// &
                                            trim(adjustl(used_target%Material(i_tar)%Elements(i_KOA)%Shell_name(i_NSH)))
!                call printout_cylindric_1d_h_single_shell(trim(adjustl(File_name)), used_target, numpar, out_data, i_arr, tim)  ! below
               call printout_cartesian_1d_particle(numpar%FN_car_1d_Z_h(i_arr), numpar%FN_car_1d_Z_E_h(i_arr), &
                trim(adjustl(m_output_cartesian_1d_Z_h))//trim(adjustl(File_name))//'_', trim(adjustl(m_output_cartesian_1d_Z_E_h))//trim(adjustl(File_name))//'_', &
                used_target, numpar, out_data%Distr_h_Z(i_arr,:),  out_data%E_Distr_h_Z(i_arr,:), tim, 3, i_arr)  ! holes; below
            enddo SH3
         enddo KOA3
      enddo TRGT3
      
   end select
end subroutine printout_cartesian_1d_hole




subroutine printout_cylindric_1d(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
    ! Along radius R:
   if (numpar%grid_par(8)%along_axis) then
      call printout_cylindric_1d_ph(used_target, numpar, out_data, tim)   ! photons; below
      call printout_cylindric_1d_e(used_target, numpar, out_data, tim)    ! electrons; below
      call printout_cylindric_1d_p(used_target, numpar, out_data, tim)    ! positrons; below
      call printout_cylindric_1d_h(used_target, numpar, out_data, tim)    ! valence holes; below
      call printout_cylindric_1d_SHI(used_target, numpar, out_data, tim)  ! shi; below
      call printout_cylindric_1d_a(used_target, numpar, out_data, tim)    ! atomic; below
      call printout_cylindric_1d_mu(used_target, numpar, out_data, tim)   ! muon; below
   endif

   ! Along depth L:
   if (numpar%grid_par(9)%along_axis) then
      ! NOT READY YET
   endif

   ! Along angle Theta:
   if (numpar%grid_par(10)%along_axis) then
      ! NOT READY YET
   endif   
end subroutine printout_cylindric_1d



subroutine printout_cylindric_1d_ph(used_target, numpar, out_data, tim)  ! photons
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_1d_R_ph))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_ph))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_ph, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_ph,'(a)') '#Time    R  Density'
         write(numpar%FN_cyl_1d_R_ph,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_ph ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%Distr_ph_R(:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(m_output_cylindric_1d_R_E_ph))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_E_ph))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_E_ph, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_E_ph,'(a)') '#Time    R    Energy'
         write(numpar%FN_cyl_1d_R_E_ph,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_E_ph ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%E_Distr_ph_R(:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_1d_ph



subroutine printout_cylindric_1d_e(used_target, numpar, out_data, tim)  ! electrons
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_1d_R_e))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_e))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_e, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_e,'(a)') '#Time    R  Density'
         write(numpar%FN_cyl_1d_R_e,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_e ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%Distr_e_R(:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_1d_R_E_e))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_E_e))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_E_e, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_E_e,'(a)') '#Time    R    Energy'
         write(numpar%FN_cyl_1d_R_E_e,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_E_e ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%E_Distr_e_R(:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_1d_e



subroutine printout_cylindric_1d_p(used_target, numpar, out_data, tim)  ! positrons
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_1d_R_p))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_p))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_p, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_p,'(a)') '#Time    R  Density'
         write(numpar%FN_cyl_1d_R_p,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_p ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%Distr_p_R(:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_1d_R_E_p))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_E_p))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_E_p, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_E_p,'(a)') '#Time    R    Energy'
         write(numpar%FN_cyl_1d_R_E_p,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_E_p ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%E_Distr_p_R(:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_1d_p



subroutine printout_cylindric_1d_mu(used_target, numpar, out_data, tim)  ! muons
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_1d_R_mu))//trim(adjustl(used_target%Name))//'.dat'

      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_mu, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_mu,'(a)') '#Time    R  Density'
         write(numpar%FN_cyl_1d_R_mu,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_mu ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%Distr_mu_R(:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_1d_R_E_mu))//trim(adjustl(used_target%Name))//'.dat'

      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_E_mu, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_E_mu,'(a)') '#Time    R    Energy'
         write(numpar%FN_cyl_1d_R_E_mu,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_E_mu ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%E_Distr_mu_R(:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_1d_mu



subroutine printout_cylindric_1d_h(used_target, numpar, out_data, tim)  ! valence holes
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, i_KOA, i_NSH, i_arr
   character(50) :: File_name
   
   ! Make sure the arrays with files numbers for holes are allocated:
   if (.not.allocated(numpar%FN_cyl_1d_R_h)) allocate(numpar%FN_cyl_1d_R_h(numpar%N_sh_tot))
   if (.not.allocated(numpar%FN_cyl_1d_R_E_h)) allocate(numpar%FN_cyl_1d_R_E_h(numpar%N_sh_tot))
   
   ! Printout electron distribution in each target:
   i_arr = 1    ! do valence band first
   call printout_cylindric_1d_h_single_shell('Valence', used_target, numpar, out_data, i_arr, tim)  ! below
   
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      KOA:do i_KOA = 1, size(used_target%Material(i_tar)%Elements)    ! for all elements within this target
         SH:do i_NSH = 1, used_target%Material(i_tar)%Elements(i_KOA)%N_core_shl    ! for all core shells of this element
               i_arr = i_arr + 1    ! next array index
               write(File_name,'(a)') trim(adjustl(used_target%Material(i_tar)%Elements(i_KOA)%Name))//'_'// &
                                            trim(adjustl(used_target%Material(i_tar)%Elements(i_KOA)%Shell_name(i_NSH)))
               call printout_cylindric_1d_h_single_shell(trim(adjustl(File_name)), used_target, numpar, out_data, i_arr, tim)  ! below
         enddo SH
      enddo KOA
   enddo TRGT
end subroutine printout_cylindric_1d_h



subroutine printout_cylindric_1d_h_single_shell(File_name_part, used_target, numpar, out_data, Narr, tim)  ! valence holes
   character(*), intent(in) :: File_name_part
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   integer, intent(in) :: Narr  ! index of the array where to use the data from
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(300) :: File_name
   integer :: FN
   logical :: file_exist
   
   ! Create file with density:
   File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                 trim(adjustl(m_output_cylindric_1d_R_h))//trim(adjustl(File_name_part))//'.dat'
   inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
   if (.not. file_exist) then   ! it's the first time, create file and write the header
      open(newunit = numpar%FN_cyl_1d_R_h(Narr), FILE = trim(adjustl(File_name)))
      write(numpar%FN_cyl_1d_R_h(Narr),'(a)') '#Time    R  Density'
      write(numpar%FN_cyl_1d_R_h(Narr),'(a)') '#fs    A   1/A^3'
   endif
   FN = numpar%FN_cyl_1d_R_h(Narr) ! just set a number
   ! Write data with densities into this file:
   call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%Distr_h_R(Narr,:), tim)    ! below

   ! Create file with energy density:
   File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                 trim(adjustl(m_output_cylindric_1d_R_E_h))//trim(adjustl(File_name_part))//'.dat'
   inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
   if (.not. file_exist) then   ! it's the first time, create file and write the header
      open(newunit = numpar%FN_cyl_1d_R_E_h(Narr), FILE = trim(adjustl(File_name)))
      write(numpar%FN_cyl_1d_R_E_h(Narr),'(a)') '#Time    R    Energy'
      write(numpar%FN_cyl_1d_R_E_h(Narr),'(a)') '#fs    A   eV/A^3'
   endif
   FN = numpar%FN_cyl_1d_R_E_h(Narr) ! just set a number
   ! Write data with energy densities into this file:
   call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%E_Distr_h_R(Narr,:), tim)    ! below
end subroutine printout_cylindric_1d_h_single_shell



subroutine printout_cylindric_1d_SHI(used_target, numpar, out_data, tim)  ! SHI
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(m_output_cylindric_1d_R_SHI))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_SHI))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_SHI, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_SHI,'(a)') '#Time    R  Density'
         write(numpar%FN_cyl_1d_R_SHI,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_SHI ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%Distr_SHI_R(:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(m_output_cylindric_1d_R_E_SHI))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_E_e))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_E_SHI, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_E_SHI,'(a)') '#Time    R    Energy'
         write(numpar%FN_cyl_1d_R_E_SHI,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_E_SHI ! just set a number
      ! Write data with energy densities into this file:
!       call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%E_Distr_e_R(:), tim)    ! below
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%E_Distr_SHI_R(:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_1d_SHI



subroutine printout_cylindric_1d_a(used_target, numpar, out_data, tim)  ! atomic energy
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(m_output_cylindric_1d_R_a))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_a))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_a, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_a,'(a)') '#Time    R  Density'
         write(numpar%FN_cyl_1d_R_a,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_a ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%Distr_a_R(:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(m_output_cylindric_1d_R_E_a))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_1d_R_E_a))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_1d_R_E_a, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_1d_R_E_a,'(a)') '#Time    R    Energy'
         write(numpar%FN_cyl_1d_R_E_a,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_1d_R_E_a ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_1d_grid(FN, numpar%grids(8)%spatial_grid1(:), out_data%E_Distr_a_R(:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_1d_a


subroutine printout_cylindric_2d(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
    ! Along radius RL:
   if (numpar%grid_par(11)%along_axis) then
      call printout_cylindric_2d_RL_ph(used_target, numpar, out_data, tim)  ! photons; below
      call printout_cylindric_2d_RL_e(used_target, numpar, out_data, tim)  ! electrons; below
      call printout_cylindric_2d_RL_p(used_target, numpar, out_data, tim)  ! positrons; below
      call printout_cylindric_2d_RL_h(used_target, numpar, out_data, tim)  ! valence holes; below
      call printout_cylindric_2d_RL_SHI(used_target, numpar, out_data, tim)  ! shi; below
      call printout_cylindric_2d_RL_a(used_target, numpar, out_data, tim)  ! atomic; below
      call printout_cylindric_2d_RL_mu(used_target, numpar, out_data, tim)  ! muons; below
   endif

   ! Along depth RTheta:
   if (numpar%grid_par(12)%along_axis) then
      ! NOT READY YET
   endif

end subroutine printout_cylindric_2d



subroutine printout_cylindric_2d_RL_ph(used_target, numpar, out_data, tim)  ! photons
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_2d_RL_ph))//trim(adjustl(used_target%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_ph, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_ph,'(a)') '#R[A]:L[A]:Density[eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_ph,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_ph ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%Distr_ph_RL(:,:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(m_output_cylindric_2d_RL_E_ph))//trim(adjustl(used_target%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_E_ph, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_E_ph,'(a)') '#R[A]:L[A]:Energy [eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_E_ph,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_E_ph ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%E_Distr_ph_RL(:,:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_2d_RL_ph


subroutine printout_cylindric_2d_RL_e(used_target, numpar, out_data, tim)  ! electrons
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_2d_RL_e))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_2d_RL_e))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_e, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_e,'(a)') '#R[A]:L[A]:Density[eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_e,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_e ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%Distr_e_RL(:,:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_2d_RL_E_e))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_2d_RL_E_e))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_E_e, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_E_e,'(a)') '#R[A]:L[A]:Energy[eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_E_e,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_E_e ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%E_Distr_e_RL(:,:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_2d_RL_e


subroutine printout_cylindric_2d_RL_p(used_target, numpar, out_data, tim)  ! positrons
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_2d_RL_p))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_2d_RL_p))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_p, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_p,'(a)') '#R[A]:L[A]:Density[eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_p,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_p ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%Distr_p_RL(:,:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_2d_RL_E_p))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_2d_RL_E_p))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_E_p, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_E_p,'(a)') '#R[A]:L[A]:Energy [eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_E_p,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_E_p ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%E_Distr_p_RL(:,:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_2d_RL_p



subroutine printout_cylindric_2d_RL_mu(used_target, numpar, out_data, tim)  ! muons
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_2d_RL_mu))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_2d_RL_mu))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_mu, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_mu,'(a)') '#R[A]:L[A]:Density[eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_mu,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_mu ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%Distr_mu_RL(:,:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_cylindric_2d_RL_E_mu))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_2d_RL_E_mu))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_E_mu, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_E_mu,'(a)') '#R[A]:L[A]:Energy [eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_E_mu,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_E_mu ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%E_Distr_mu_RL(:,:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_2d_RL_mu


subroutine printout_cylindric_2d_RL_h(used_target, numpar, out_data, tim)  ! valence holes
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, i_KOA, i_NSH, i_arr
   character(50) :: File_name
   
   ! Make sure the arrays with files numbers for holes are allocated:
   if (.not.allocated(numpar%FN_cyl_2d_RL_h)) allocate(numpar%FN_cyl_2d_RL_h(numpar%N_sh_tot))
   if (.not.allocated(numpar%FN_cyl_2d_RL_E_h)) allocate(numpar%FN_cyl_2d_RL_E_h(numpar%N_sh_tot))
   
   ! Printout electron distribution in each target:
   i_arr = 1    ! do valence band first
   call printout_cylindric_2d_RL_h_single_shell('Valence', used_target, numpar, out_data, i_arr, tim)  ! below
   
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      KOA:do i_KOA = 1, size(used_target%Material(i_tar)%Elements)    ! for all elements within this target
         SH:do i_NSH = 1, used_target%Material(i_tar)%Elements(i_KOA)%N_core_shl    ! for all core shells of this element
               i_arr = i_arr + 1    ! next array index
               write(File_name,'(a)') trim(adjustl(used_target%Material(i_tar)%Elements(i_KOA)%Name))//'_'// &
                                            trim(adjustl(used_target%Material(i_tar)%Elements(i_KOA)%Shell_name(i_NSH)))
               call printout_cylindric_2d_RL_h_single_shell(trim(adjustl(File_name)), used_target, numpar, out_data, i_arr, tim)  ! below
         enddo SH
      enddo KOA
   enddo TRGT
end subroutine printout_cylindric_2d_RL_h



subroutine printout_cylindric_2d_RL_h_single_shell(File_name_part, used_target, numpar, out_data, Narr, tim)  ! valence holes
   character(*), intent(in) :: File_name_part
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   integer, intent(in) :: Narr  ! index of the array where to use the data from
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(300) :: File_name
   integer :: FN
   logical :: file_exist
   
   ! Create file with density:
   File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                 trim(adjustl(m_output_cylindric_2d_RL_h))//trim(adjustl(File_name_part))//'.dat'
   inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
   if (.not. file_exist) then   ! it's the first time, create file and write the header
      open(newunit = numpar%FN_cyl_2d_RL_h(Narr), FILE = trim(adjustl(File_name)))
      write(numpar%FN_cyl_2d_RL_h(Narr),'(a)') '#R[A]:L[A]:Density[eV/A^3]'
      !write(numpar%FN_cyl_2d_RL_h(Narr),'(a)') '#fs    A   1/A^3'
   endif
   FN = numpar%FN_cyl_2d_RL_h(Narr) ! just set a number
   ! Write data with densities into this file:
   call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                out_data%Distr_h_RL(Narr,:,:), tim)    ! below

   ! Create file with energy density:
   File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                 trim(adjustl(m_output_cylindric_2d_RL_E_h))//trim(adjustl(File_name_part))//'.dat'
   inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
   if (.not. file_exist) then   ! it's the first time, create file and write the header
      open(newunit = numpar%FN_cyl_2d_RL_E_h(Narr), FILE = trim(adjustl(File_name)))
      write(numpar%FN_cyl_2d_RL_E_h(Narr),'(a)') '#R[A]:L[A]:Energy[eV/A^3]'
      !write(numpar%FN_cyl_2d_RL_E_h(Narr),'(a)') '#fs    A   eV/A^3'
   endif
   FN = numpar%FN_cyl_2d_RL_E_h(Narr) ! just set a number
   ! Write data with energy densities into this file:
   call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:),&
                                out_data%E_Distr_h_RL(Narr,:,:), tim)    ! below
end subroutine printout_cylindric_2d_RL_h_single_shell


subroutine printout_cylindric_2d_RL_SHI(used_target, numpar, out_data, tim)  ! SHI
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(m_output_cylindric_2d_RL_SHI))//trim(adjustl(used_target%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_SHI, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_SHI,'(a)') '#R[A]:L[A]:Density[eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_SHI,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_SHI ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%Distr_SHI_RL(:,:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(m_output_cylindric_2d_RL_E_SHI))//trim(adjustl(used_target%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_E_SHI, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_E_SHI,'(a)') '#R[A]:L[A]:Energy[eV/A^3]'
!         write(numpar%FN_cyl_2d_RL_E_SHI,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_E_SHI ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%E_Distr_e_RL(:,:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_2d_RL_SHI



subroutine printout_cylindric_2d_RL_a(used_target, numpar, out_data, tim)  ! atomic energy
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   integer :: i_tar, FN
   character(300) :: File_name
   logical :: file_exist
   ! Printout electron distribution in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create file with density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(m_output_cylindric_2d_RL_a))//trim(adjustl(used_target%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_a, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_a,'(a)') '#R[A]:L[A]:Density[eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_a,'(a)') '#fs    A   1/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_a ! just set a number
      ! Write data with densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%Distr_a_RL(:,:), tim)    ! below

      ! Create file with energy density:
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                     trim(adjustl(m_output_cylindric_2d_RL_E_a))//trim(adjustl(used_target%Name))//'.dat'
                     !trim(adjustl(m_output_cylindric_2d_RL_E_a))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_cyl_2d_RL_E_a, FILE = trim(adjustl(File_name)))
         write(numpar%FN_cyl_2d_RL_E_a,'(a)') '#R[A]:L[A]:Energy[eV/A^3]'
         !write(numpar%FN_cyl_2d_RL_E_a,'(a)') '#fs    A   eV/A^3'
      endif
      FN = numpar%FN_cyl_2d_RL_E_a ! just set a number
      ! Write data with energy densities into this file:
      call printout_data_on_2d_grid(FN, numpar%grids(11)%spatial_grid1(:), numpar%grids(11)%spatial_grid2(:), &
                                    out_data%E_Distr_a_RL(:,:), tim)    ! below
   enddo TRGT
end subroutine printout_cylindric_2d_RL_a


subroutine printout_data_on_1d_grid(FN, grid_array, data_array, tim)
   integer, intent(in) :: FN    ! file number (file must be already open!)
   real(8), dimension(:), intent(in) :: grid_array  ! array with the grid
   real(8), dimension(:), intent(in) :: data_array  ! 1d array with the data to printout
   real(8), intent(in) :: tim   ! [fs] current time step to print out
   !------------------------------
   integer :: Nsiz, i
   Nsiz = size(grid_array)
   if (Nsiz == size(data_array)) then   ! make sure the sizes match
      do i = 1, Nsiz
          write(FN,'(f16.6,es24.16,es24.16)') tim, grid_array(i), data_array(i)
      enddo
      write(FN,'(a)') ! skip line between timesteps
   else ! something is wrong
      write(6,'(a)') 'Error in printout_data_on_1d_grid: size of grid does not match size of data array'
   endif
end subroutine printout_data_on_1d_grid

subroutine printout_data_on_2d_grid(FN, grid_array1, grid_array2, data_array, tim)
   integer, intent(in) :: FN    ! file number (file must be already open!)
   real(8), dimension(:), intent(in) :: grid_array1, grid_array2  ! array with the grid
   real(8), dimension(:,:), intent(in) :: data_array  ! 1d array with the data to printout
   real(8), intent(in) :: tim   ! [fs] current time step to print out
   !------------------------------
   integer :: Nsiz1, Nsiz2, i, j
   
   Nsiz1 = size(grid_array1)
   Nsiz2 = size(grid_array2)
   if (Nsiz1 == size(data_array, 1)) then   ! make sure the sizes match
        if (Nsiz2 == size(data_array, 2)) then   ! make sure the sizes match      
            write(FN,'(f16.6,A)') tim, ' fs' !printout header with time
            write(FN,'(A)', advance='no') ' Gr1/Gr2   '
            !Printout first row with grid_array1 vaues
            do i = 1, Nsiz2 
                write(FN,'(es24.16)', advance='no'), grid_array2(i)
            enddo
            write(FN,'(A)') ' ' !Jump to the next line after row
            do i = 1, Nsiz1
                write(FN,'(es24.16)', advance='no'), grid_array1(i)
                do j = 1, Nsiz2
                     write(FN,'(es24.16)', advance='no') data_array(i,j)
                enddo
                write(FN,'(A)') ' ' !Jump to the next line after row
            enddo
            write(FN,'(a)') ! skip line between timesteps
        else ! something is wrong
            write(6,'(a)') 'Error in printout_data_on_2d_grid: size of grid2 does not match size of data array'
        endif    
   else ! something is wrong
      write(6,'(a)') 'Error in printout_data_on_2d_grid: size of grid1 does not match size of data array'
   endif
end subroutine printout_data_on_2d_grid


subroutine energy_spectra(used_target, numpar, out_data, tim)
    type(Matter), intent(in) :: used_target      ! parameters of the target
    type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
    type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
    real(8), intent(in) :: tim  ! simulation time step [fs]
    !---------------------------------
    if (numpar%NRG_grid_par%along_axis) then ! if user requested energy data
       ! Printout photon spectrum:
       call photon_spectrum(used_target, numpar, out_data, tim)   ! below
       ! Printout electron spectrum:
       call electron_spectrum(used_target, numpar, out_data, tim)   ! below
       ! Printout valence hole spectrum:
       call hole_spectrum(used_target, numpar, out_data, tim)   ! below
       ! Printout positron spectrum:
       call positron_spectrum(used_target, numpar, out_data, tim)   ! below
       ! SHI spectum:
        ! NOT READY YET (because we typically don't need it...)
       ! Printout muon spectrum:
       call muon_spectrum(used_target, numpar, out_data, tim)   ! below
    endif ! (numpar%NRG_grid_par%along_axis)
    !---------------------------------
    ! Printout electron spectra vs space in 1d:    
    call electron_spectra_1d(used_target, numpar, out_data, tim)    ! below
end subroutine energy_spectra



subroutine total_values(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !---------------------------------
   logical :: file_exist, file_exists2, do_high
   integer :: i_tar
   real(8) :: Etot, Etot_high
   
   ! Allocate the file names if needed:
   file_exist = .true.  ! to start with
   if (.not.allocated(numpar%FILE_totals)) allocate(numpar%FILE_totals(used_target%NOC))
   if (.not.allocated(numpar%FN_totals)) then
      allocate(numpar%FN_totals(used_target%NOC))
      file_exist = .false.
   endif

   ! Set files for particels with energies above cut-off, only if cut-offs are used:
   file_exists2 = .true.    ! to start with
   do_high = .false.    ! to start with
   if ( (numpar%El_Cutoff > 0.0d0) .or. (numpar%H_Cutoff > 0.0d0) .or. (numpar%Ph_Cutoff > 0.0d0) ) then
      do_high = .true.  ! there are cut-offs used, save the output data
      if (.not.allocated(numpar%FILE_totals_cutoff)) allocate(numpar%FILE_totals_cutoff(used_target%NOC))
      if (.not.allocated(numpar%FN_totals_cutoff)) then
         allocate(numpar%FN_totals_cutoff(used_target%NOC))
         file_exists2 = .false.
      endif
   endif
   
   ! Printout total values in each target:
   !TRGT:do i_tar = 1, used_target%NOC ! for all targets
   ! Total nmumbers for the entire target, not different materials:
   i_tar = 1
      if (.not.file_exist) then ! Create a file:
         numpar%FILE_totals(i_tar) = trim(adjustl(numpar%output_path))//numpar%path_sep// &
         !           trim(adjustl(m_output_total))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
                    trim(adjustl(m_output_total))//'all.dat'
      endif
      inquire(file=trim(adjustl(numpar%FILE_totals(i_tar))),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_totals(i_tar), FILE = trim(adjustl(numpar%FILE_totals(i_tar))))
         write(numpar%FN_totals(i_tar),'(a)') '#Time    Nph Ne Nh  Np  Eph Ee  Eh_kin  Eh_pot  Ep  Eat  Etot'
         write(numpar%FN_totals(i_tar),'(a)') '#fs    a.u. a.u.    a.u.    a.u.    eV  eV  eV  eV  eV  eV  eV'
      endif
      ! Get the total energy in the target to print out:
      Etot = out_data%Eph + out_data%Ee + out_data%Eh_kin + out_data%Eh_pot + out_data%Ep + out_data%Eat
      write(numpar%FN_totals(i_tar),'(es24.16, $)') tim, out_data%Nph, out_data%Ne, out_data%Nh, out_data%Np, &
              out_data%Eph, out_data%Ee, out_data%Eh_kin, out_data%Eh_pot, out_data%Ep, out_data%Eat, Etot
      write(numpar%FN_totals(i_tar),'(a)') ''

      ! The same for high-energy particles (above cut-off):
      if (do_high) then
         if (.not.file_exists2) then   ! create files
            numpar%FILE_totals_cutoff(i_tar) = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    !trim(adjustl(m_output_total_cutoff))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
                    trim(adjustl(m_output_total_cutoff))//'all.dat'
         endif
         inquire(file=trim(adjustl(numpar%FILE_totals_cutoff(i_tar))),exist=file_exists2) ! check if input file is there
         if (.not. file_exists2) then   ! it's the first time, create file and write the header
            open(newunit = numpar%FN_totals_cutoff(i_tar), FILE = trim(adjustl(numpar%FILE_totals_cutoff(i_tar))))
            write(numpar%FN_totals_cutoff(i_tar),'(a)') '#Time    Nph Ne Nh  Np  Eph Ee  Eh_kin  Eh_pot  Ep  Etot'
            write(numpar%FN_totals_cutoff(i_tar),'(a)') '#fs    a.u. a.u.    a.u.    a.u.    eV  eV  eV  eV  eV  eV'
         endif
         ! Get the total energy of high-energy particles (above cut-off) in the target to print out:
         Etot_high = out_data%Eph_high + out_data%Ee_high + out_data%Eh_kin_high + out_data%Eh_pot_high + out_data%Ep_high
         write(numpar%FN_totals_cutoff(i_tar),'(es24.16, $)') tim, out_data%Nph_high, out_data%Ne_high, out_data%Nh_high, &
            out_data%Np_high, out_data%Eph_high, out_data%Ee_high, out_data%Eh_kin_high, out_data%Eh_pot_high, &
            out_data%Ep_high, Etot_high
         write(numpar%FN_totals_cutoff(i_tar),'(a)') ''
      endif ! (do_high)
   !enddo TRGT
end subroutine total_values



subroutine create_gnuplot_files(used_target, MD_pots, numpar)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(MD_potential), dimension(:,:), allocatable, intent(in) :: MD_pots    ! MD potentials
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   !---------------------------------
   ! Create gnuplot scripts for plotting total numbers and energies:
   call gnuplot_total_values(used_target, numpar)  ! below

   ! Create gnuplot for MD part:
   if (numpar%DO_MD) then
      call gnuplot_MD_values(used_target, MD_pots, numpar)  ! below
   endif
end subroutine create_gnuplot_files



subroutine gnuplot_MD_values(used_target, MD_pots, numpar)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   !---------------------------------
   integer :: i_tar
   !TRGT:do i_tar = 1, used_target%NOC ! for all targets
   i_tar = 1    ! for global target only, not each material
      ! Create a file for atomic temperature:
      call create_MD_energies_gnuplot(used_target%Material(i_tar), numpar, &
      trim(adjustl(m_output_MD_energies)), numpar%t_start, numpar%t_total, log_x = .false.)   ! below

      ! Create a file for atomic temperature:
      call create_MD_temperature_gnuplot(used_target%Material(i_tar), numpar, &
      trim(adjustl(m_output_MD_cell_params)), numpar%t_start, numpar%t_total, log_x = .false.)   ! below

      ! Create a file for atomic displacements:
      call create_MD_displacement_gnuplot(MD_pots, used_target%Material(i_tar), numpar, &
      trim(adjustl(m_output_MD_displacements)), numpar%t_start, numpar%t_total, log_x = .false.)   ! below

!    enddo TRGT
end subroutine gnuplot_MD_values



subroutine create_MD_energies_gnuplot(Material, numpar, Datafile, x_start, x_end, log_x)
   type(Target_atoms), intent(in), target :: Material ! parameters of this material
   type(Num_par), intent(in), optional :: numpar	! all numerical parameters
   character(*), intent(in) :: Datafile
   real(8), intent(in) :: x_start, x_end
   logical, intent(in), optional :: log_x
   !------------------------
   character(200) :: File_script, Out_file
   character(50) :: Title, temp, temp2
   character(10) :: units
   character(5) ::  call_slash, sh_cmd, col_y
   logical :: logx
   integer :: FN_gnu, i_first, Reason
   real(8) :: tics, ord

   if (present(log_x)) then
      if (log_x) then   ! user set it to make x-axis logscale
         logx = .true.
      else  ! x axis linear
         logx = .false.
      endif
   else ! x axis linear
      logx = .false.
   endif

   ! Set the grid step on the plots:
   if (logx) then
      tics = 10.0d0
   else
      ord = dble( find_order_of_number( abs(x_start-x_end) ) - 2 )   ! module "Little_subroutines"

      !tics = 10.0d0**ord
      write(temp2,'(es)')  abs(x_start-x_end) ! make it a string
      temp = trim(adjustl(temp2))
      read(temp(1:1),*,IOSTAT=Reason) i_first

      if (i_first < 3) then
         tics = 10.0d0**ord
      else
         tics = 10.0d0**(ord+1)
      endif
   endif

   ! Get the extension and slash in this OS:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"

   ! Printout total values in each target:
   ! Get the paths and file names:
   File_script = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_MD_E_gnu))//trim(adjustl(Material%Name))//trim(adjustl(sh_cmd))
   open(newunit = FN_gnu, FILE = trim(adjustl(File_script)))
   Out_file = trim(adjustl(m_output_MD_E_gnu))//'in_'//trim(adjustl(Material%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))

   ! Create the gnuplot-script header:
   call write_gnuplot_script_header_new(FN_gnu, 1, 3.0d0, tics, "Energies vs Time", "Time (fs)", "Energy (eV/atom)", &
      trim(adjustl(Out_file)), trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, setkey=0, &
      logx=logx, logy=.false.)  ! module "Gnuplotting"

   ! Create the plotting part:
   write(Title, '(a)') ' Total'
   write(col_y, '(i4)') 4   ! in this column there is Etot
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .true., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .true., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   write(Title, '(a)') ' Potential'
   write(col_y, '(i4)') 3   ! in this column there is Eph
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif

   ! Create the gnuplot-script ending:
   call  write_gnuplot_script_ending_new(FN_gnu, File_script, numpar%path_sep)  ! module "Gnuplotting"

   call close_file('save',FN=FN_gnu)
end subroutine create_MD_energies_gnuplot




subroutine create_MD_temperature_gnuplot(Material, numpar, Datafile, x_start, x_end, log_x)
   type(Target_atoms), intent(in), target :: Material ! parameters of this material
   type(Num_par), intent(in), optional :: numpar	! all numerical parameters
   character(*), intent(in) :: Datafile
   real(8), intent(in) :: x_start, x_end
   logical, intent(in), optional :: log_x
   !------------------------
   character(200) :: File_script, Out_file
   character(50) :: Title, temp, temp2
   character(10) :: units
   character(5) ::  call_slash, sh_cmd, col_y
   logical :: logx
   integer :: FN_gnu, i_first, Reason
   real(8) :: tics, ord

   if (present(log_x)) then
      if (log_x) then   ! user set it to make x-axis logscale
         logx = .true.
      else  ! x axis linear
         logx = .false.
      endif
   else ! x axis linear
      logx = .false.
   endif

   ! Set the grid step on the plots:
   if (logx) then
      tics = 10.0d0
   else
      ord = dble(find_order_of_number( abs(x_start-x_end) ) - 2)   ! module "Little_subroutines"
      write(temp2,'(es)')  abs(x_start-x_end) ! make it a string
      temp = trim(adjustl(temp2))
      read(temp(1:1),*,IOSTAT=Reason) i_first
      if (i_first < 3) then
         tics = 10.0d0**ord
      else
         tics = 10.0d0**(ord+1)
      endif
   endif

   ! Get the extension and slash in this OS:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"

   ! Printout total values in each target:
   ! Get the paths and file names:
   File_script = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_MD_T_gnu))//trim(adjustl(Material%Name))//trim(adjustl(sh_cmd))
   open(newunit = FN_gnu, FILE = trim(adjustl(File_script)))
   Out_file = trim(adjustl(m_output_MD_T_gnu))//'in_'//trim(adjustl(Material%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))

   ! Create the gnuplot-script header:
   call write_gnuplot_script_header_new(FN_gnu, 1, 3.0d0, tics, "Temperature vs Time", "Time (fs)", "Temperature (K)", &
                  trim(adjustl(Out_file)), trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, &
                  setkey=0, logx=logx, logy=.false.)  ! module "Gnuplotting"

   ! Create the plotting part:
   write(Title, '(a)') ' Atoms'
   write(col_y, '(i4)') 2   ! in this column there is mean atomic temperature
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .true., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .true., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif

   ! Create the gnuplot-script ending:
   call  write_gnuplot_script_ending_new(FN_gnu, File_script, numpar%path_sep)  ! module "Gnuplotting"

   call close_file('save',FN=FN_gnu)
end subroutine create_MD_temperature_gnuplot



subroutine create_MD_displacement_gnuplot(MD_pots, Material, numpar, Datafile, x_start, x_end, log_x)
   type(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(Target_atoms), intent(in), target :: Material ! parameters of this material
   type(Num_par), intent(in), optional :: numpar	! all numerical parameters
   character(*), intent(in) :: Datafile
   real(8), intent(in) :: x_start, x_end
   logical, intent(in), optional :: log_x
   !------------------------
   character(200) :: File_script, Out_file
   character(50) :: Title, temp, temp2
   character(10) :: units, chtemp
   character(5) ::  call_slash, sh_cmd, col_y
   logical :: logx
   integer :: FN_gnu, i_first, Reason, N_KOA, i
   real(8) :: tics, ord

   ! Number of different kinds of atoms (defined by different potentials):
   N_KOA = size(MD_pots,1)

   if (present(log_x)) then
      if (log_x) then   ! user set it to make x-axis logscale
         logx = .true.
      else  ! x axis linear
         logx = .false.
      endif
   else ! x axis linear
      logx = .false.
   endif

   ! Set the grid step on the plots:
   if (logx) then
      tics = 10.0d0
   else
      ord = dble(find_order_of_number( abs(x_start-x_end) ) - 2)   ! module "Little_subroutines"
      write(temp2,'(es)')  abs(x_start-x_end) ! make it a string
      temp = trim(adjustl(temp2))
      read(temp(1:1),*,IOSTAT=Reason) i_first
      if (i_first < 3) then
         tics = 10.0d0**ord
      else
         tics = 10.0d0**(ord+1)
      endif
   endif

   ! Get the extension and slash in this OS:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"

   ! Printout total values in each target:
   ! Get the paths and file names:
   File_script = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_MD_MSD_gnu))//trim(adjustl(Material%Name))//trim(adjustl(sh_cmd))
   open(newunit = FN_gnu, FILE = trim(adjustl(File_script)))
   Out_file = trim(adjustl(m_output_MD_MSD_gnu))//'in_'//trim(adjustl(Material%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))

   ! Create the gnuplot-script header:
   if (numpar%n_MSD /= 1) then
      write(chtemp,'(i2)') numpar%n_MSD
      call write_gnuplot_script_header_new(FN_gnu, 1, 3.0d0, tics, "Displacement vs Time", "Time (fs)", &
                  "Displacement (A^"//trim(adjustl(chtemp))//')', &
                  trim(adjustl(Out_file)), trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, &
                  setkey=0, logx=logx, logy=.false.)  ! module "Gnuplotting"
   else
      call write_gnuplot_script_header_new(FN_gnu, 1, 3.0d0, tics, "Displacement vs Time", "Time (fs)", "Displacement (A)", &
                  trim(adjustl(Out_file)), trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, &
                  setkey=0, logx=logx, logy=.false.)  ! module "Gnuplotting"
   endif

   ! Create the plotting part:
   if (N_KOA > 1) then ! many kinds of atoms
      write(Title, '(a)') ' Average'
      write(col_y, '(i4)') 2   ! in this column there is mean atomic temperature
      if (numpar%path_sep == '\') then	! if it is Windows
         call write_gnu_printout(FN_gnu, .true., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), &
            x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
      else
         call write_gnu_printout(FN_gnu, .true., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), &
            x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
      endif
      do i = 3, (2+N_KOA)-1 ! for all kinds of atoms
         write(Title, '(a)') trim(adjustl(MD_pots(i-2,i-2)%El1))
         write(col_y, '(i4)') i   ! in this column there is mean atomic temperature
         if (numpar%path_sep == '\') then	! if it is Windows
            call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), &
               x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
         else
            call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), &
               x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
         endif
      enddo
      write(Title, '(a)') trim(adjustl(MD_pots(N_KOA,N_KOA)%El1))
      write(col_y, '(i4)') 2+N_KOA   ! in this column there is mean atomic temperature
      if (numpar%path_sep == '\') then	! if it is Windows
         call write_gnu_printout(FN_gnu, .false., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)), &
            x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
      else
         call write_gnu_printout(FN_gnu, .false., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)), &
            x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
      endif
   else ! only one kind of atoms
      write(Title, '(a)') ' Average'
      write(col_y, '(i4)') 2   ! in this column there is mean atomic temperature
      if (numpar%path_sep == '\') then	! if it is Windows
         call write_gnu_printout(FN_gnu, .true., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)), &
            x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
      else
         call write_gnu_printout(FN_gnu, .true., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)), &
            x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
      endif
   endif ! (N_KOA > 1)

   ! Create the gnuplot-script ending:
   call  write_gnuplot_script_ending_new(FN_gnu, File_script, numpar%path_sep)  ! module "Gnuplotting"

   call close_file('save',FN=FN_gnu)
end subroutine create_MD_displacement_gnuplot




subroutine gnuplot_total_values(used_target, numpar)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   !---------------------------------
   integer :: i_tar
   !TRGT:do i_tar = 1, used_target%NOC ! for all targets
   i_tar = 1    ! for global target only, not each material
      ! Create a file for total numbers:
      call create_total_numbers_gnuplot(used_target%Material(i_tar), numpar, &
      trim(adjustl(m_output_total))//'all.dat', numpar%t_start, numpar%t_total, log_x = .false.)   ! below
      !trim(adjustl(m_output_total))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat', numpar%t_start, numpar%t_total, log_x = .false.)   ! below


      ! Create a file for total eneries:
      call create_total_energies_gnuplot(used_target%Material(i_tar), numpar, &
      trim(adjustl(m_output_total))//'all.dat', numpar%t_start, numpar%t_total, log_x = .false.)   ! below
      !trim(adjustl(m_output_total))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat', numpar%t_start, numpar%t_total, log_x = .false.)   ! below
!    enddo TRGT
end subroutine gnuplot_total_values




subroutine create_total_numbers_gnuplot(Material, numpar, Datafile, x_start, x_end, log_x)
   type(Target_atoms), intent(in), target :: Material ! parameters of this material
   type(Num_par), intent(in), optional :: numpar	! all numerical parameters
   character(*), intent(in) :: Datafile
   real(8), intent(in) :: x_start, x_end
   logical, intent(in), optional :: log_x
   !------------------------
   character(200) :: File_script, Out_file
   character(50) :: Title, temp, temp2
   character(10) :: units
   character(5) ::  call_slash, sh_cmd, col_y
   logical :: logx
   integer :: FN_gnu, i_first, Reason
   real(8) :: tics, ord
   
   if (present(log_x)) then
      if (log_x) then   ! user set it to make x-axis logscale
         logx = .true.
      else  ! x axis linear
         logx = .false.
      endif
   else ! x axis linear
      logx = .false.
   endif
   
   ! Set the grid step on the plots:
   if (logx) then
      tics = 10.0d0
   else
      ord = dble(find_order_of_number( abs(x_start-x_end) ) - 2)   ! module "Little_subroutines"
      write(temp2,'(es)')  abs(x_start-x_end) ! make it a string
      temp = trim(adjustl(temp2))
      read(temp(1:1),*,IOSTAT=Reason) i_first
!       print*, 'CHECK', i_first
!       pause 'Check'
      if (i_first < 3) then
         tics = 10.0d0**ord
      else
         tics = 10.0d0**(ord+1)
      endif
   endif
   
   ! Get the extension and slash in this OS:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"
   
   ! Printout total values in each target:
   ! Get the paths and file names:
   File_script = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_N_gnu))//trim(adjustl(Material%Name))//trim(adjustl(sh_cmd))
   open(newunit = FN_gnu, FILE = trim(adjustl(File_script)))
   Out_file = trim(adjustl(m_output_N_gnu))//'in_'//trim(adjustl(Material%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))
   
   ! Create the gnuplot-script header:
   call write_gnuplot_script_header_new(FN_gnu, 1, 3.0d0, tics, "Numbers vs Time", "Time (fs)", "Numbers (arb. units)", trim(adjustl(Out_file)), &
            trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, setkey=2, logx=logx, logy=.false.)  ! module "Gnuplotting"
   
   ! Create the plotting part:
   write(Title, '(a)') ' Photons'
   write(col_y, '(i4)') 2   ! in this column there is Nph
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .true., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .true., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif  
   write(Title, '(a)') ' Electrons'
   write(col_y, '(i4)') 3   ! in this column there is Nph
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, &
                 lw=3, title=trim(adjustl(Title)), additional_info = 'lt -1')  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, &
                 lw=3, title=trim(adjustl(Title)), additional_info = 'lt -1', linux_s =.true.)  ! module "Gnuplotting"
   endif
   write(Title, '(a)') ' All holes'
   write(col_y, '(i4)') 4   ! in this column there is Nph
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, &
               lw=1, title=trim(adjustl(Title)) )  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, &
               lw=1, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   write(Title, '(a)') ' Positrons'
   write(col_y, '(i4)') 5   ! in this column there is Nph
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   
   ! Create the gnuplot-script ending:
   call  write_gnuplot_script_ending_new(FN_gnu, File_script, numpar%path_sep)  ! module "Gnuplotting"
   
   call close_file('save',FN=FN_gnu)
end subroutine create_total_numbers_gnuplot




subroutine create_total_energies_gnuplot(Material, numpar, Datafile, x_start, x_end, log_x)
   type(Target_atoms), intent(in), target :: Material ! parameters of this material
   type(Num_par), intent(in), optional :: numpar	! all numerical parameters
   character(*), intent(in) :: Datafile
   real(8), intent(in) :: x_start, x_end
   logical, intent(in), optional :: log_x
   !------------------------
   character(200) :: File_script, Out_file
   character(50) :: Title, temp, temp2
   character(10) :: units
   character(5) ::  call_slash, sh_cmd, col_y
   logical :: logx
   integer :: FN_gnu, i_first, Reason
   real(8) :: tics, ord
   
   if (present(log_x)) then
      if (log_x) then   ! user set it to make x-axis logscale
         logx = .true.
      else  ! x axis linear
         logx = .false.
      endif
   else ! x axis linear
      logx = .false.
   endif
   
   ! Set the grid step on the plots:
   if (logx) then
      tics = 10.0d0
   else
      ord = dble( find_order_of_number( abs(x_start-x_end) ) - 2 )   ! module "Little_subroutines"
      write(temp2,'(es)')  abs(x_start-x_end) ! make it a string
      temp = trim(adjustl(temp2))
      read(temp(1:1),*,IOSTAT=Reason) i_first
      if (i_first < 3) then
         tics = 10.0d0**ord
      else
         tics = 10.0d0**(ord+1)
      endif
   endif
   
   ! Get the extension and slash in this OS:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"
   
   ! Printout total values in each target:
   ! Get the paths and file names:
   File_script = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_E_gnu))//trim(adjustl(Material%Name))//trim(adjustl(sh_cmd))
   open(newunit = FN_gnu, FILE = trim(adjustl(File_script)))
!    Out_file = trim(adjustl(m_output_E_gnu))//'in_'//trim(adjustl(Material%Name))//'.eps'
   Out_file = trim(adjustl(m_output_E_gnu))//'in_'//trim(adjustl(Material%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))
   
   ! Create the gnuplot-script header:
   call write_gnuplot_script_header_new(FN_gnu, 1, 3.0d0, tics, "Energies vs Time", "Time (fs)", "Energy (eV)", trim(adjustl(Out_file)), &
            trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, setkey=0, logx=logx, logy=.false.)  ! module "Gnuplotting"
   
   ! Create the plotting part:
   write(Title, '(a)') ' Total'
   write(col_y, '(i4)') 12   ! in this column there is Etot
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .true., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .true., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   write(Title, '(a)') ' Photons'
   write(col_y, '(i4)') 6   ! in this column there is Eph
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif  
   write(Title, '(a)') ' Electrons'
   write(col_y, '(i4)') 7   ! in this column there is Ee
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, &
                 lw=3, title=trim(adjustl(Title)), additional_info = 'lt -1')  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, &
                 lw=3, title=trim(adjustl(Title)), additional_info = 'lt -1', linux_s =.true.)  ! module "Gnuplotting"
   endif
   write(Title, '(a)') ' Holes (kin)'
   write(col_y, '(i4)') 8   ! in this column there is Eh_pot
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, &
               lw=3, title=trim(adjustl(Title)) )  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, &
               lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   write(Title, '(a)') ' Holes (pot)'
   write(col_y, '(i4)') 9   ! in this column there is Eh_pot
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, &
               lw=3, title=trim(adjustl(Title)) )  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, &
               lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   write(Title, '(a)') ' Positrons'
   write(col_y, '(i4)') 10   ! in this column there is Ep
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .false., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   write(Title, '(a)') ' Atoms'
   write(col_y, '(i4)') 11   ! in this column there is Ep
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif

   ! Create the gnuplot-script ending:
   call  write_gnuplot_script_ending_new(FN_gnu, File_script, numpar%path_sep)  ! module "Gnuplotting"
   
   call close_file('save',FN=FN_gnu)
end subroutine create_total_energies_gnuplot




subroutine electron_spectra_1d(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in), target :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   character(3) :: axis_name
   integer :: FN, i_tar, i, Nsiz, j, Nsiz2
   logical :: file_opened, file_exist
   real(8), dimension(:), pointer :: Spectrum_cur
   
   ! Create a directory:
   if (numpar%Spectr_grid_par(1)%along_axis) then  ! Along X
      axis_name = 'X'
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_spectrum_e_1d))//trim(adjustl(axis_name))//'.dat'                 
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_spectrum_e_X, FILE = trim(adjustl(File_name)))
         write(numpar%FN_spectrum_e_X,'(a)') '#Time    E Distribution_vs_space_X'
         write(numpar%FN_spectrum_e_X,'(a)') '#fs    eV 1/eV    etc.'
         write(numpar%FN_spectrum_e_X,'(a)', advance='no') '#   '
         write(numpar%FN_spectrum_e_X,'(es15.3,$)') numpar%Spectr_grid(1)%spatial_grid1(:)
         write(numpar%FN_spectrum_e_X,'(a)')    ! to start new line
      endif
      FN = numpar%FN_spectrum_e_X ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(numpar%NRG_grid)
       do i = 1, Nsiz
          Spectrum_cur => out_data%Spectra_e_X(i,:)
          Nsiz2 = size(Spectrum_cur)
          write(FN,'(es24.16, $)') tim, numpar%NRG_grid(i), ( Spectrum_cur(j), j=1,Nsiz2-1 )
          !write(FN,'(f16.6, es24.16, $)') tim, numpar%NRG_grid(i), out_data%Spectra_e_X(i,:)
          write(FN,'(a)')   ! start a new line
       enddo
       write(FN,'(a)') ! skip line between timesteps
   endif    ! X
   
   if (numpar%Spectr_grid_par(2)%along_axis) then  ! Along Y
      axis_name = 'Y'
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_spectrum_e_1d))//trim(adjustl(axis_name))//'.dat'                 
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_spectrum_e_Y, FILE = trim(adjustl(File_name)))
         write(numpar%FN_spectrum_e_Y,'(a)') '#Time    E Distribution_vs_space_Y'
         write(numpar%FN_spectrum_e_Y,'(a)') '#fs    eV 1/eV    etc.'
         write(numpar%FN_spectrum_e_Y,'(a)', advance='no') '#   '
         write(numpar%FN_spectrum_e_Y,'(es15.3,$)') numpar%Spectr_grid(2)%spatial_grid1(:)
         write(numpar%FN_spectrum_e_Y,'(a)')    ! to start new line
      endif
      FN = numpar%FN_spectrum_e_Y ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(numpar%NRG_grid)
       do i = 1, Nsiz
          Spectrum_cur => out_data%Spectra_e_Y(i,:)
          Nsiz2 = size(Spectrum_cur)
          write(FN,'(es24.16, $)') tim, numpar%NRG_grid(i), ( Spectrum_cur(j), j=1,Nsiz2-1 )
          !write(FN,'(f16.6, es24.16, $)') tim, numpar%NRG_grid(i), out_data%Spectra_e_Y(i,:)
          write(FN,'(a)')   ! start a new line
       enddo
       write(FN,'(a)') ! skip line between timesteps
   endif    ! Y
   
   if (numpar%Spectr_grid_par(3)%along_axis) then  ! Along Z
      axis_name = 'Z'
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_spectrum_e_1d))//trim(adjustl(axis_name))//'.dat'                 
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_spectrum_e_Z, FILE = trim(adjustl(File_name)))
         write(numpar%FN_spectrum_e_Z,'(a)') '#Time    E Distribution_vs_space_Z'
         write(numpar%FN_spectrum_e_Z,'(a)') '#fs    eV 1/eV    etc.'
         write(numpar%FN_spectrum_e_Z,'(a)', advance='no') '#   '
         write(numpar%FN_spectrum_e_Z,'(e15.3,$)') numpar%Spectr_grid(3)%spatial_grid1(:)
         write(numpar%FN_spectrum_e_Z,'(a)')    ! to start new line
      endif
      FN = numpar%FN_spectrum_e_Z ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(numpar%NRG_grid)
!        print*, 'SIZES:', size(numpar%NRG_grid), size(out_data%Spectra_e_Z, 1)
       do i = 1, Nsiz
          Spectrum_cur => out_data%Spectra_e_Z(i,:)
          Nsiz2 = size(Spectrum_cur)
          write(FN,'(es24.16, $)') tim, numpar%NRG_grid(i), ( Spectrum_cur(j), j=1,Nsiz2-1 )
          !write(FN,'(f16.6, es24.16, $)') tim, numpar%NRG_grid(i), out_data%Spectra_e_Z(i,:)
          write(FN,'(a)')   ! start a new line
       enddo
       write(FN,'(a)') ! skip line between timesteps
       
!        print*, 'electron_spectra_1d'
       
   endif    ! Z
   
   nullify (Spectrum_cur)
end subroutine electron_spectra_1d



subroutine electron_theta_distr_1d(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in), target :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   character(3) :: axis_name
   integer :: FN, i_tar, i, Nsiz, j, Nsiz2
   logical :: file_opened, file_exist
   real(8), dimension(:), pointer :: Theta_distr_cur

   ! Create a directory:
   if (numpar%Theta_grid_par(1)%along_axis) then  ! Along X
      axis_name = 'X'
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_theta_e_1d))//trim(adjustl(axis_name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_theta_e_X, FILE = trim(adjustl(File_name)))
         write(numpar%FN_theta_e_X,'(a)') '#Time    Theta Distribution_vs_space_X'
         write(numpar%FN_theta_e_X,'(a)') '#fs    deg 1/deg    etc.'
         write(numpar%FN_theta_e_X,'(a)', advance='no') '#   '
         write(numpar%FN_theta_e_X,'(es15.3,$)') numpar%Theta_grid(1)%spatial_grid1(:)
         write(numpar%FN_theta_e_X,'(a)')    ! to start new line
      endif
      FN = numpar%FN_theta_e_X ! just set a number
      ! Size of the array for output theta:
      Nsiz = size(out_data%Vel_theta_ph)

      do i = 1, Nsiz
          Theta_distr_cur => out_data%Theta_e_X(i,:)
          Nsiz2 = size(Theta_distr_cur)
          write(FN,'(es24.16, $)') tim, numpar%vel_theta_grid(i), ( Theta_distr_cur(j), j=1,Nsiz2-1 )
          !write(FN,'(f16.6, es24.16, $)') tim, numpar%vel_theta_grid(i), out_data%Theta_e_X(i,:)
          write(FN,'(a)')   ! start a new line
       enddo
       write(FN,'(a)') ! skip line between timesteps
   endif    ! X

   if (numpar%Theta_grid_par(2)%along_axis) then  ! Along Y
      axis_name = 'Y'
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_theta_e_1d))//trim(adjustl(axis_name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_theta_e_Y, FILE = trim(adjustl(File_name)))
         write(numpar%FN_theta_e_Y,'(a)') '#Time    Theta Distribution_vs_space_Y'
         write(numpar%FN_theta_e_Y,'(a)') '#fs    deg 1/deg    etc.'
         write(numpar%FN_theta_e_Y,'(a)', advance='no') '#   '
         write(numpar%FN_theta_e_Y,'(es15.3,$)') numpar%Theta_grid(2)%spatial_grid1(:)
         write(numpar%FN_theta_e_Y,'(a)')    ! to start new line
      endif
      FN = numpar%FN_theta_e_Y ! just set a number
      ! Size of the array for output theta:
      Nsiz = size(out_data%Vel_theta_ph)

      do i = 1, Nsiz
          Theta_distr_cur => out_data%Theta_e_Y(i,:)
          Nsiz2 = size(Theta_distr_cur)
          write(FN,'(es24.16, $)') tim, numpar%vel_theta_grid(i), ( Theta_distr_cur(j), j=1,Nsiz2-1 )
          !write(FN,'(f16.6, es24.16, $)') tim, numpar%vel_theta_grid(i), out_data%Theta_e_Y(i,:)
          write(FN,'(a)')   ! start a new line
       enddo
       write(FN,'(a)') ! skip line between timesteps
   endif    ! Y

   if (numpar%Theta_grid_par(3)%along_axis) then  ! Along Z
      axis_name = 'Z'
      File_name = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_theta_e_1d))//trim(adjustl(axis_name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_theta_e_Z, FILE = trim(adjustl(File_name)))
         write(numpar%FN_theta_e_Z,'(a)') '#Time    Theta Distribution_vs_space_Z'
         write(numpar%FN_theta_e_Z,'(a)') '#fs    deg 1/deg    etc.'
         write(numpar%FN_theta_e_Z,'(a)', advance='no') '#   '
         write(numpar%FN_theta_e_Z,'(e15.3,$)') numpar%Theta_grid(3)%spatial_grid1(:)
         write(numpar%FN_theta_e_Z,'(a)')    ! to start new line
      endif
      FN = numpar%FN_theta_e_Z ! just set a number
      ! Size of the array for output theta:
      Nsiz = size(out_data%Vel_theta_ph)
!        print*, 'SIZES:', size(numpar%vel_theta_grid), size(out_data%Theta_e_Z, 1)
      do i = 1, Nsiz
          Theta_distr_cur => out_data%Theta_e_Z(i,:)
          Nsiz2 = size(Theta_distr_cur)
          write(FN,'(es24.16, $)') tim, numpar%vel_theta_grid(i), ( Theta_distr_cur(j), j=1,Nsiz2-1 )
          !write(FN,'(f16.6, es24.16, $)') tim, numpar%vel_theta_grid(i), out_data%Theta_e_Z(i,:)
          write(FN,'(a)')   ! start a new line
       enddo
       write(FN,'(a)') ! skip line between timesteps
   endif    ! Z

   nullify (Theta_distr_cur)
end subroutine electron_theta_distr_1d




subroutine photon_spectrum(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist
   
   ! Printout electron spectrum in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_spectrum_ph = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_spectrum_ph))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_spectrum_ph)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_spectrum_ph, FILE = trim(adjustl(numpar%FILE_spectrum_ph)))
         write(numpar%FN_spectrum_ph,'(a)') '#Time    E Distribution'
         write(numpar%FN_spectrum_ph,'(a)') '#fs    eV 1/eV'
      endif
      FN = numpar%FN_spectrum_ph ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(out_data%Spectrum_ph)
!       print*, 'photon_spectrum', Nsiz, size(numpar%NRG_grid)
      do i = 1, Nsiz
         write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%NRG_grid(i), out_data%Spectrum_ph(i)
      enddo
      write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine photon_spectrum



subroutine electron_spectrum(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist
   
   ! Printout electron spectrum in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_spectrum_e = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_spectrum_e))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_spectrum_e)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_spectrum_e, FILE = trim(adjustl(numpar%FILE_spectrum_e)))
         write(numpar%FN_spectrum_e,'(a)') '#Time    E Distribution'
         write(numpar%FN_spectrum_e,'(a)') '#fs    eV 1/eV'
      endif
      FN = numpar%FN_spectrum_e ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(out_data%Spectrum_e)
       do i = 1, Nsiz
          write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%NRG_grid(i), out_data%Spectrum_e(i)
       enddo
       write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine electron_spectrum



subroutine hole_spectrum(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist
   
   if (.not.allocated(numpar%FILE_spectrum_h)) then
      allocate(numpar%FILE_spectrum_h(used_target%NOC))
      allocate(numpar%FN_spectrum_h(used_target%NOC))
   endif

   ! Printout electron spectrum in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_spectrum_h(i_tar) = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_spectrum_h))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_spectrum_h(i_tar))),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_spectrum_h(i_tar), FILE = trim(adjustl(numpar%FILE_spectrum_h(i_tar))))
         write(numpar%FN_spectrum_h(i_tar),'(a)') '#Time    E Distribution'
         write(numpar%FN_spectrum_h(i_tar),'(a)') '#fs    eV 1/eV'
      endif
      FN = numpar%FN_spectrum_h(i_tar) ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(numpar%NRG_grid_VB)
      do i = 1, Nsiz
         write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%NRG_grid_VB(i), out_data%Spectrum_h(i_tar,i)
      enddo
      write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine hole_spectrum



subroutine positron_spectrum(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist
   
   ! Printout electron spectrum in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_spectrum_p = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_spectrum_p))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_spectrum_p)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_spectrum_p, FILE = trim(adjustl(numpar%FILE_spectrum_p)))
         write(numpar%FN_spectrum_p,'(a)') '#Time    E Distribution'
         write(numpar%FN_spectrum_p,'(a)') '#fs    eV 1/eV'
      endif
      FN = numpar%FN_spectrum_p ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(out_data%Spectrum_p)
       do i = 1, Nsiz
          write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%NRG_grid(i), out_data%Spectrum_p(i)
       enddo
       write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine positron_spectrum



subroutine muon_spectrum(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist

   ! Printout electron spectrum in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_spectrum_mu = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_spectrum_mu))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_spectrum_mu)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_spectrum_mu, FILE = trim(adjustl(numpar%FILE_spectrum_mu)))
         write(numpar%FN_spectrum_mu,'(a)') '#Time    E Distribution'
         write(numpar%FN_spectrum_mu,'(a)') '#fs    eV 1/eV'
      endif
      FN = numpar%FN_spectrum_mu ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(out_data%Spectrum_mu)
       do i = 1, Nsiz
          write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%NRG_grid(i), out_data%Spectrum_mu(i)
       enddo
       write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine muon_spectrum



!===================================================
! Printing out velocity theta distribution:
subroutine photon_velotheta(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist
   
   ! Printout photon velocity theta in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_vel_theta_ph = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_velocity_theta_distr_ph))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_vel_theta_ph)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_vel_theta_ph, FILE = trim(adjustl(numpar%FILE_vel_theta_ph)))
         write(numpar%FN_vel_theta_ph,'(a)') '#Time    Theta Distribution'
         write(numpar%FN_vel_theta_ph,'(a)') '#fs    deg 1/deg'
      endif
      FN = numpar%FN_vel_theta_ph ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(out_data%Vel_theta_ph)
       do i = 1, Nsiz
          write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%vel_theta_grid(i), out_data%Vel_theta_ph(i)
       enddo
       write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine photon_velotheta


subroutine electron_velotheta(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist
   
   ! Printout electron velocity theta in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_vel_theta_e = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_velocity_theta_distr_e))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_vel_theta_e)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_vel_theta_e, FILE = trim(adjustl(numpar%FILE_vel_theta_e)))
         write(numpar%FN_vel_theta_e,'(a)') '#Time    Theta Distribution'
         write(numpar%FN_vel_theta_e,'(a)') '#fs    deg 1/deg'
      endif
      FN = numpar%FN_vel_theta_e ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(out_data%Vel_theta_e)
       do i = 1, Nsiz
          write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%vel_theta_grid(i), out_data%Vel_theta_e(i)
       enddo
       write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine electron_velotheta


subroutine positron_velotheta(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist
   
   ! Printout electron velocity theta in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_vel_theta_p = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_velocity_theta_distr_p))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_vel_theta_p)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_vel_theta_p, FILE = trim(adjustl(numpar%FILE_vel_theta_p)))
         write(numpar%FN_vel_theta_p,'(a)') '#Time    Theta Distribution'
         write(numpar%FN_vel_theta_p,'(a)') '#fs    deg 1/deg'
      endif
      FN = numpar%FN_vel_theta_p ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(out_data%Vel_theta_p)
       do i = 1, Nsiz
          write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%vel_theta_grid(i), out_data%Vel_theta_p(i)
       enddo
       write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine positron_velotheta



subroutine muon_velotheta(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist

   ! Printout electron velocity theta in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_vel_theta_mu = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_velocity_theta_distr_mu))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_vel_theta_mu)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_vel_theta_mu, FILE = trim(adjustl(numpar%FILE_vel_theta_mu)))
         write(numpar%FN_vel_theta_mu,'(a)') '#Time    Theta Distribution'
         write(numpar%FN_vel_theta_mu,'(a)') '#fs    deg 1/deg'
      endif
      FN = numpar%FN_vel_theta_mu ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(out_data%Vel_theta_mu)
       do i = 1, Nsiz
          write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%vel_theta_grid(i), out_data%Vel_theta_mu(i)
       enddo
       write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine muon_velotheta



subroutine hole_velotheta(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist
   
   ! Printout electron velocity theta in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_vel_theta_h = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_velocity_theta_distr_h))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_vel_theta_h)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_vel_theta_h, FILE = trim(adjustl(numpar%FILE_vel_theta_h)))
         write(numpar%FN_vel_theta_h,'(a)') '#Time    Theta Distribution'
         write(numpar%FN_vel_theta_h,'(a)') '#fs    deg 1/deg'
      endif
      FN = numpar%FN_vel_theta_h ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(out_data%Vel_theta_h)
       do i = 1, Nsiz
          write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%vel_theta_grid(i), out_data%Vel_theta_h(i)
       enddo
       write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine hole_velotheta


subroutine SHI_velotheta(used_target, numpar, out_data, tim)
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(inout), target :: numpar    ! all numerical parameters
   type(output_data), intent(in) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim  ! simulation time step [fs]
   !----------------------------------------
   character(200) :: File_name
   integer :: FN, i_tar, i, Nsiz
   logical :: file_opened, file_exist
   
   ! Printout electron velocity theta in each target:
   TRGT:do i_tar = 1, used_target%NOC ! for all targets
      ! Create a directory for MFPs:
      numpar%FILE_vel_theta_SHI = trim(adjustl(numpar%output_path))//numpar%path_sep// &
                    trim(adjustl(m_output_velocity_theta_distr_SHI))//trim(adjustl(used_target%Material(i_tar)%Name))//'.dat'
      inquire(file=trim(adjustl(numpar%FILE_vel_theta_SHI)),exist=file_exist) ! check if input file is there
      if (.not. file_exist) then   ! it's the first time, create file and write the header
         open(newunit = numpar%FN_vel_theta_SHI, FILE = trim(adjustl(numpar%FILE_vel_theta_SHI)))
         write(numpar%FN_vel_theta_SHI,'(a)') '#Time    Theta Distribution'
         write(numpar%FN_vel_theta_SHI,'(a)') '#fs    deg 1/deg'
      endif
      FN = numpar%FN_vel_theta_SHI ! just set a number
      ! Size of the array for output spectrum:
      Nsiz = size(out_data%Vel_theta_SHI)
       do i = 1, Nsiz
          write(FN,'(f16.6,es24.16,es24.16)') tim, numpar%vel_theta_grid(i), out_data%Vel_theta_SHI(i)
       enddo
       write(FN,'(a)') ! skip line between timesteps
   enddo TRGT
end subroutine SHI_velotheta


 
!===================================================
subroutine print_MFPs(used_target, numpar, bunch)    ! create all output files
   type(Matter), intent(in) :: used_target      ! parameters of the target
   type(Num_par), intent(in) :: numpar    ! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   ! All MFPs:
   call printout_MFPs(used_target, numpar, bunch)  ! below
   ! All Ranges:
   call printout_Se_and_ranges(used_target, numpar, bunch)  ! below
end subroutine print_MFPs


!===================================================
! Printing out stopping powers and ranges:
subroutine printout_Se_and_ranges(used_target, numpar, bunch)
   type(Matter), intent(in), target :: used_target      ! parameters of the target
   type(Num_par), intent(in), target :: numpar    ! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   !----------------------------------------------------------
   character(250) :: Path, command, Filename_gnu, Filename
   character(50) :: Model_name
   character(5) :: Ion_name
   character(250), dimension(:), allocatable :: File_name_a, File_for_gnu
   integer, dimension(:), allocatable :: FNA
   real(8) :: x_start, x_end, y_start, y_end, eps, E_start, E_start_ord
   integer :: i, j, k, FN, Nsiz, m, N_types_SHI, N_bunch, ipart, ibunch, iret, N_elements, N_shells
   logical :: file_exist
   type(Cross_section), pointer :: CS
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   real(8), pointer :: E
   
   if (numpar%printout_ranges) then    ! do it only if a user requested it
      path_sep => numpar%path_sep
      ! Prepare file with the data output:
      TRGT:do i = 1, used_target%NOC ! for all targets
         ! Create a directory for MFPs:
         Path = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_folder_MFP))//trim(adjustl(used_target%Material(i)%Name))
         inquire(DIRECTORY=trim(adjustl(Path)),exist=file_exist)    ! check if input file excists
         if (.not.file_exist) then ! to make sure that such a folder is present (even if empty)
            ! Create a new directory for output files:
            command='mkdir '//trim(adjustl(Path))  ! to create a folder use this command
#ifdef _OPENMP
            iret = system(command)   ! create the folder
#else
            call system(command)   ! create the folder
#endif
         endif
         
         ! Define arrays needed for each element of the target:
         N_elements = size(used_target%Material(i)%Elements)	! that's how many different elements are in this target
         if (allocated(File_name_a)) deallocate(File_name_a)
         allocate(File_name_a(N_elements))
         if (allocated(File_for_gnu)) deallocate(File_for_gnu)
         allocate(File_for_gnu(N_elements))
         if (allocated(FNA)) deallocate(FNA)
         allocate(FNA(N_elements))
         
         
         ! Calculate, printout, and plot stoppings and ranges:
         ! 1) Electrons:
         ! Pointer to the cross section:
         CS => used_target%Material(i)%El_inelastic_total
         
         ! Inelastic:
         select case (numpar%El_inelast)
         case (1)	! CDF with Ritchi's oscillators
            Model_name = '_CDF_Ritchie'
         case (2)	! RBEB
            Model_name = '_RBEB'
         case (3)	! CDF with delta functions
            Model_name = '_CDF_delta'
         case (4)	! CDF 
            Model_name = '_CDF_nonrel'
         case (5)	! CDF with delta functions
            Model_name = '_CDF_SPdelta'
         case default	! exclude
            Model_name = '_NO'
         end select
         
         ! Electron inelastic Se and range:
         Filename_gnu= trim(adjustl(m_output_Se))//trim(adjustl(used_target%Material(i)%Name))//'_electron'//trim(adjustl(Model_name))//'.dat'
         Filename = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(Filename_gnu))
         open(newunit = FN, FILE = trim(adjustl(Filename)))
         write(FN,'(a)') '#Energy(eV)    Se(eV/A)   Range(A)'
         Nsiz = size(CS%E)
         ! Assume the lower integration limit for electronic range calculation as:
         E_start = used_target%Material(i)%DOS%Egap + 10.0d0    ! [eV]
         !E_start = 10.0d0    ! [eV] testing standard definition
         ! Calculate Ranges:
         call get_ranges_from_Se(CS%E, CS%Total_Se, E_start, CS%Total_Range)    ! module "CS_general_tools"
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN,'(es,es,es)')  CS%E(m), CS%Total_Se(m), CS%Total_Range(m)
         enddo
         call close_file('save', FN=FN)  ! module "Dealing_with_files"
         
         ! Define axes for the plot of Se:
         x_start = 1.0d0
         x_end = 10.0**dble(-1 + find_order_of_number(CS%E(size(CS%E))) ) ! module "Little_subroutines"
         y_start = 0.0d0
         y_end = 10.0d0
         ! Prepare gnuplot script:
         call create_Se_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, 'electron') ! below
         
         ! Define axes for plotting Range:
         x_start = 10.0d0
         x_end = 10.0**dble(-1 + find_order_of_number(CS%E(size(CS%E))) ) ! module "Little_subroutines"
         y_start = 10.0d0
         y_end = 1.0d12
         ! Prepare gnuplot script:
         call create_Range_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, 'electron')    ! below
         x_start = 10.0d0
         x_end = 1.0d12
         call create_Se_vs_Range_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, 'electron')    ! below
         
         
         ! 2) Ions:
         ! Check if there is any SHI to calcualte:
         N_types_SHI = count_types_of_SHIs(bunch) ! module "Initial_conditions"
         ! Only if there are some ions, calculate MFPs for them
         if (N_types_SHI > 0) then  ! there are ions
            ! Find how many sources of radiation are used:
            N_bunch = size (bunch)
            ipart = 0
            ! Check for all of them, whether there are ions (do we need to include ions here or not):
            BNCH:do ibunch = 1, N_bunch
               if (bunch(ibunch)%KOP == 3) then ! it is an ion
                  if (.not.repeated_type(bunch, ibunch)) then   ! new type of SHI, create its plots:
                  
                     Ion_name = bunch(ibunch)%Name  ! this is the name of the SHI
                     ipart = ipart + 1
                     select case (numpar%El_inelast)    ! so far, we can only have CDF-based models for SHI
                     case (4)	! non-relativistic CDF with Ritchi's oscillators
                        Model_name = '_CDF_Ritchie'
                     case default	! delta-function CDF
                        Model_name = '_CDF'
                     end select
                  
                     ! Partial stopping:
                     ! Prepare the files for each element to be written into:
                     do j =1, N_elements	! for each element
                        Element => used_target%Material(i)%Elements(j)	! all information about this element
                        File_for_gnu(j) = trim(adjustl(m_output_Se))//trim(adjustl(used_target%Material(i)%Name))//'_' &
                        //trim(adjustl(Element%Name))//'_'//trim(adjustl(Ion_name))//trim(adjustl(Model_name))//'.dat'
                        
                        File_name_a(j) = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu(j)))
                        open(newunit = FNA(j), FILE = trim(adjustl(File_name_a(j))))
                        ! Create the header:
                        write(FNA(j),'(a)', advance='no') '#Energy(MeV) '
                        call write_shell_headers_output(FNA(j), used_target%Material(i), Element, numpar=numpar)  ! below        
                     enddo
                     ! Save core-shells data into the file:
                     do j =1, N_elements	! for each element
                        Element => used_target%Material(i)%Elements(j)	! all information about this element
                        Nsiz = size(Element%SHI_inelastic(ipart)%E)
                        do m = 1, Nsiz	! for all energy grid points:
                           E => Element%SHI_inelastic(ipart)%E(m)       ! Ion energy [eV]
                           write(FNA(j),'(es)', advance='no') E*1.0d-6   ! energy [MeV]
                           N_shells = Element%N_shl
                           do k = 1, N_shells
                              ! Check valence band
                              if ((numpar%Pos_inelast /= 2) .and. (Element%valent(k)) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
                                 ! Valence band will be added at the end
                              else    ! core shell
                                 write(FNA(j),'(es)', advance='no') Element%SHI_inelastic(ipart)%Se(k,m) ! Se [eV/A]
                              endif
                           enddo ! k = 1, N_shells
                           write(FNA(j),'(a)') ''  ! to go to the next line
                        enddo ! m = 1, Nsiz
                        call close_file('save', FN=FNA(j))  ! module "Dealing_with_files"
                     enddo !  j =1, N_elements	! for each element
                     ! Total and Valence stopping:
                     ! Pointer to the cross section:
                     CS => used_target%Material(i)%SHI_inelastic_total(ipart)
                     Filename_gnu= trim(adjustl(m_output_Se))//trim(adjustl(used_target%Material(i)%Name))//'_'//trim(adjustl(Ion_name))//trim(adjustl(Model_name))//'_Total.dat'
                     Filename = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(Filename_gnu))
                     open(newunit = FN, FILE = trim(adjustl(Filename)))
                     write(FN,'(a)') '#Energy(MeV)    Se(eV/A)   Range(A)   Valence'
                     ! Assume the lower integration limit for electronic range calculation as:
                     E_start_ord = 10.0d0
                     E_start = bunch(ibunch)%Meff*g_amu/(4.0d0*g_me)*(used_target%Material(i)%DOS%Egap + E_start_ord)    ! [MeV]
                     ! Calculate Ranges:
                     call get_ranges_from_Se(CS%E, CS%Total_Se, E_start, CS%Total_Range)    ! module "CS_general_tools"
                     do m = 1, Nsiz     ! for all energy grid points:
                        write(FN,'(es,es,es,es)')  CS%E(m)*1.0d-6, CS%Total_Se(m), CS%Total_Range(m), &
                                                used_target%Material(i)%SHI_inelastic_valent(ipart)%Total_Se(m)
                     enddo
                     call close_file('save', FN=FN)  ! module "Dealing_with_files"
                  
                     ! Define axes for the plot of Se:
                     ! (The plotting limits are chosen for no physical reason, just such that the curves fall inside the plot)
                     x_start = E_start*1.0d-6
                     x_end = (bunch(ibunch)%Z/500.0d0) * 1.0d-6*10.0**dble(-2 + find_order_of_number(CS%E(size(CS%E))) ) ! module "Little_subroutines"
                     y_start = 0.0d0
                     y_end = 5000.0d0
                     ! Prepare gnuplot script:
                     call create_Se_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, trim(adjustl(Ion_name)), y_units='(MeV)') ! below
         
                     ! Define axes for plotting Range:
                     ! (The plotting limits are chosen for no physical reason, just such that the curves fall inside the plot)
                     y_start = 10.0d0 * (100.0d0/E_start_ord)**2
                     y_end = 1.0d12 * (E_start_ord/100.0d0)**2
                     ! Prepare gnuplot script:
                     call create_Range_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, trim(adjustl(Ion_name)), y_units='(MeV)') ! below
                     ! (The plotting limits are chosen for no physical reason, just such that the curves fall inside the plot)
                     x_start = 1.0d6 * sqrt(bunch(ibunch)%Z/100.0d0)
                     x_end = sqrt(bunch(ibunch)%Z/100.0d0) * 10.0**dble(-1 + find_order_of_number(CS%Total_Range(size(CS%E))) ) ! module "Little_subroutines"
                     call create_Se_vs_Range_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, trim(adjustl(Ion_name)))    ! below
                  endif ! (.not.repeated_type(bunch, ibunch))
               endif ! (bunch(ibunch)%KOP == 3)
            enddo BNCH
         endif ! (N_types_SHI > 0)


         ! 3) Positrons:
         ! Pointer to the cross section:
         CS => used_target%Material(i)%Pos_inelastic_total

         ! Inelastic:
         select case (numpar%Pos_inelast)
         case (1)	! CDF with Ritchi's oscillators
            Model_name = '_CDF_Ritchie'
         case (2)	! RBEB
            Model_name = '_RBEB'
         case (3)	! CDF with delta functions
            Model_name = '_CDF_delta'
         case (4)	! CDF
            Model_name = '_CDF_nonrel'
         case (5)	! CDF with delta functions
            Model_name = '_CDF_SPdelta'
         case default	! exclude
            Model_name = '_NO'
         end select

         ! Positron inelastic Se and range:
         Filename_gnu= trim(adjustl(m_output_Se))//trim(adjustl(used_target%Material(i)%Name))//'_posiron'//trim(adjustl(Model_name))//'.dat'
         Filename = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(Filename_gnu))
         open(newunit = FN, FILE = trim(adjustl(Filename)))
         write(FN,'(a)') '#Energy(eV)    Se(eV/A)   Range(A)'
         Nsiz = size(CS%E)
         ! Assume the lower integration limit for electronic range calculation as:
         E_start = used_target%Material(i)%DOS%Egap + 10.0d0    ! [eV]
         !E_start = 10.0d0    ! [eV] testing standard definition
         ! Calculate Ranges:
         call get_ranges_from_Se(CS%E, CS%Total_Se, E_start, CS%Total_Range)    ! module "CS_general_tools"
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN,'(es,es,es)')  CS%E(m), CS%Total_Se(m), CS%Total_Range(m)
         enddo
         call close_file('save', FN=FN)  ! module "Dealing_with_files"

         if (numpar%Pos_inelast /= 0) then ! gnuplot only if used
            ! Define axes for the plot of Se:
            x_start = 1.0d0
            x_end = 10.0**dble(-1 + find_order_of_number(CS%E(size(CS%E))) ) ! module "Little_subroutines"
            y_start = 0.0d0
            y_end = 10.0d0
            ! Prepare gnuplot script:
            call create_Se_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, 'positron') ! below

            ! Define axes for plotting Range:
            x_start = 10.0d0
            x_end = 10.0**dble(-1 + find_order_of_number(CS%E(size(CS%E))) ) ! module "Little_subroutines"
            y_start = 10.0d0
            y_end = 1.0d12
            ! Prepare gnuplot script:
            call create_Range_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, 'positron')    ! below
            x_start = 10.0d0
            x_end = 1.0d12
            call create_Se_vs_Range_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, 'positron')    ! below
         endif

         ! 4) Muon:
         ! Pointer to the cross section:
         CS => used_target%Material(i)%Muon_inelastic_total

         ! Inelastic:
         select case (numpar%Mu_inelast)
         case (1)	! CDF with Ritchi's oscillators
            Model_name = '_CDF_Ritchie'
         case (2)	! RBEB
            Model_name = '_RBEB'
         case (3)	! CDF with delta functions
            Model_name = '_CDF_delta'
         case (4)	! CDF
            Model_name = '_CDF_nonrel'
         case (5)	! CDF with delta functions
            Model_name = '_CDF_SPdelta'
         case default	! exclude
            Model_name = '_NO'
         end select

         ! Positron inelastic Se and range:
         Filename_gnu= trim(adjustl(m_output_Se))//trim(adjustl(used_target%Material(i)%Name))//'_muon'//trim(adjustl(Model_name))//'.dat'
         Filename = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(Filename_gnu))
         open(newunit = FN, FILE = trim(adjustl(Filename)))
         write(FN,'(a)') '#Energy(eV)    Se(eV/A)   Range(A)'
         Nsiz = size(CS%E)
         ! Assume the lower integration limit for electronic range calculation as:
         E_start = used_target%Material(i)%DOS%Egap + 10.0d0    ! [eV]
         !E_start = 10.0d0    ! [eV] testing standard definition
         ! Calculate Ranges:
         call get_ranges_from_Se(CS%E, CS%Total_Se, E_start, CS%Total_Range)    ! module "CS_general_tools"
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN,'(es,es,es)')  CS%E(m), CS%Total_Se(m), CS%Total_Range(m)
         enddo
         call close_file('save', FN=FN)  ! module "Dealing_with_files"

         if (numpar%Mu_inelast /= 0) then ! gnuplot only if used
            ! Define axes for the plot of Se:
            x_start = 10.0d0
            x_end = 10.0**dble(-1 + find_order_of_number(CS%E(size(CS%E))) ) ! module "Little_subroutines"
            y_start = 0.0d0
            y_end = 10.0d0
            ! Prepare gnuplot script:
            call create_Se_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, 'muon') ! below

            ! Define axes for plotting Range:
            x_start = 10.0d0
            x_end = 10.0**dble(-1 + find_order_of_number(CS%E(size(CS%E))) ) ! module "Little_subroutines"
            y_start = 10.0d0
            y_end = 1.0d12
            ! Prepare gnuplot script:
            call create_Range_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, 'muon')    ! below
            x_start = 10.0d0
            x_end = 1.0d12
            call create_Se_vs_Range_gnuplot(used_target%Material(i), numpar, trim(adjustl(Filename_gnu)), x_start, x_end, y_start, y_end, 'muon')    ! below
         endif

         
         !=========================
         ! Execute all gnuplot scripts for MFPs and Ranges:
         if (numpar%gnupl%do_gnuplot) call execute_all_gnuplots(numpar, Path)    ! below
         
      enddo TRGT
   endif ! (numpar%printout_ranges) 
   nullify(CS)
end subroutine printout_Se_and_ranges


subroutine create_Se_vs_Range_gnuplot(Material, numpar, Datafile, x_start, x_end, y_start, y_end, particle_name, y_units)
   type(Target_atoms), intent(in), target :: Material ! parameters of this material
   type(Num_par), intent(in), optional :: numpar	! all numerical parameters
   character(*), intent(in) :: Datafile
   real(8), intent(in) :: x_start, x_end, y_start, y_end
   character(*), intent(in) :: particle_name
   character(*), intent(in), optional :: y_units
   !------------------------
   character(200) :: File_script, Out_file, Path
   character(50) :: Title
   character(10) :: units
   character(5) ::  call_slash, sh_cmd, col_y
   integer :: FN_gnu
   real(8) :: tics
   
  if (numpar%gnupl%do_gnuplot) then ! do only if user wants plots
   
   if (present(y_units)) then
      units = y_units  ! user provided units
   else
      units = '(A)'    ! by default, assume eV
   endif
   
   tics = 10.0d0
   
   ! Get the extension and slash in this OS:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"
   
    ! Get the paths and file names:
   Path = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_folder_MFP))//trim(adjustl(Material%Name))
   File_script = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(m_output_Se_vs_range))//trim(adjustl(Material%Name))//'_'// &
                      trim(adjustl(particle_name))//trim(adjustl(sh_cmd))
   open(newunit = FN_gnu, FILE = trim(adjustl(File_script)))
!    Out_file = trim(adjustl(m_output_Se_vs_range))//trim(adjustl(particle_name))//'_in_'//trim(adjustl(Material%Name))//'.eps'
   Out_file = trim(adjustl(m_output_Se_vs_range))//trim(adjustl(particle_name))//'_in_'//trim(adjustl(Material%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))
   
   ! Create the gnuplot-script header:
   call write_gnuplot_script_header_new(FN_gnu, 1, 3.0d0, tics, "Se vs Range", "Range (A)", "Stopping power Se (eV/A)", trim(adjustl(Out_file)), &
            trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, setkey=0, logx=.true., logy=.false.)  ! module "Gnuplotting"
   
   ! Create the plotting part:
   write(col_y, '(i4)') 2   ! in this column there is Se
   write(Title, '(a)') trim(adjustl(particle_name))//' inelastic Se'
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .true., .true., Datafile, col_x="3", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .true., .true., Datafile, col_x="3", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   
   ! Create the gnuplot-script ending:
   call  write_gnuplot_script_ending_new(FN_gnu, File_script, numpar%path_sep)  ! module "Gnuplotting"
   
   call close_file('save',FN=FN_gnu)
   
  endif ! (numpar%gnupl%do_gnuplot) then ! do only if user wants plots
end subroutine create_Se_vs_Range_gnuplot


subroutine create_Se_gnuplot(Material, numpar, Datafile, x_start, x_end, y_start, y_end, particle_name, y_units)
   type(Target_atoms), intent(in), target :: Material ! parameters of this material
   type(Num_par), intent(in), optional :: numpar	! all numerical parameters
   character(*), intent(in) :: Datafile
   real(8), intent(in) :: x_start, x_end, y_start, y_end
   character(*), intent(in) :: particle_name
   character(*), intent(in), optional :: y_units
   !------------------------
   character(200) :: File_script, Out_file, Path
   character(50) :: Title
   character(10) :: units
   character(5) ::  call_slash, sh_cmd, col_y
   integer :: FN_gnu
   real(8) :: tics
   
  if (numpar%gnupl%do_gnuplot) then ! do only if user wants plots
   
   if (present(y_units)) then
      units = y_units  ! user provided units
   else
      units = '(eV)'    ! by default, assume eV
   endif
   
   tics = 10.0d0
   
   ! Get the extension and slash in this OS:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"
   
    ! Get the paths and file names:
   Path = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_folder_MFP))//trim(adjustl(Material%Name))
   File_script = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(m_output_Se))//trim(adjustl(Material%Name))//'_'// &
                      trim(adjustl(particle_name))//trim(adjustl(sh_cmd))
   open(newunit = FN_gnu, FILE = trim(adjustl(File_script)))
   Out_file = trim(adjustl(m_output_Se))//trim(adjustl(particle_name))//'_in_'//trim(adjustl(Material%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))
   
   ! Create the gnuplot-script header:
   call write_gnuplot_script_header_new(FN_gnu, 1, 3.0d0, tics, "Se", "Energy "//trim(adjustl(units)), "Stopping power Se (eV/A)", trim(adjustl(Out_file)), &
            trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, setkey=0, logx=.true., logy=.false.)  ! module "Gnuplotting"
   
   ! Create the plotting part:
   write(col_y, '(i4)') 2   ! in this column there is Se
   write(Title, '(a)') trim(adjustl(particle_name))//' inelastic Se'
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .true., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .true., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   
   ! Create the gnuplot-script ending:
   call  write_gnuplot_script_ending_new(FN_gnu, File_script, numpar%path_sep)  ! module "Gnuplotting"
   
   call close_file('save',FN=FN_gnu)
   
 endif ! (numpar%gnupl%do_gnuplot) then ! do only if user wants plots
end subroutine create_Se_gnuplot



subroutine create_Range_gnuplot(Material, numpar, Datafile, x_start, x_end, y_start, y_end, particle_name, y_units)
   type(Target_atoms), intent(in), target :: Material ! parameters of this material
   type(Num_par), intent(in), optional :: numpar	! all numerical parameters
   character(*), intent(in) :: Datafile
   real(8), intent(in) :: x_start, x_end, y_start, y_end
   character(*), intent(in) :: particle_name
   character(*), intent(in), optional :: y_units
   !------------------------
   character(200) :: File_script, Out_file, Path
   character(50) :: Title
   character(10) :: units
   character(5) ::  call_slash, sh_cmd, col_y
   integer :: FN_gnu
   real(8) :: tics
   
  if (numpar%gnupl%do_gnuplot) then ! do only if user wants plots

   if (present(y_units)) then
      units = y_units  ! user provided units
   else
      units = '(eV)'    ! by default, assume eV
   endif
   
   tics = 10.0d0
   
   ! Get the extension and slash in this OS:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"
   
    ! Get the paths and file names:
   Path = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_folder_MFP))//trim(adjustl(Material%Name))
   File_script = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(m_output_Range))//trim(adjustl(Material%Name))//'_'// &
                      trim(adjustl(particle_name))//trim(adjustl(sh_cmd))
   open(newunit = FN_gnu, FILE = trim(adjustl(File_script)))
   Out_file = trim(adjustl(m_output_Range))//trim(adjustl(particle_name))//'_in_'//trim(adjustl(Material%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))
   
   ! Create the gnuplot-script header:
   call write_gnuplot_script_header_new(FN_gnu, 1, 3.0d0, tics, "Range", "Energy "//trim(adjustl(units)), "Range (A)", trim(adjustl(Out_file)), &
            trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, setkey=0, logx=.true., logy=.true.)  ! module "Gnuplotting"
   
   ! Create the plotting part:
   write(col_y, '(i4)') 3   ! in this column there is Range
   write(Title, '(a)') 'Inelastic '//trim(adjustl(particle_name))//' range'
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .true., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)),  x_start=x_start, x_end=x_end, y_start=y_start, y_end=y_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .true., .true., Datafile, col_x="1", col_y=trim(adjustl(col_y)), x_start=x_start, x_end=x_end, y_start=y_start, y_end=y_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   
   ! Create the gnuplot-script ending:
   call  write_gnuplot_script_ending_new(FN_gnu, File_script, numpar%path_sep)  ! module "Gnuplotting"
   
   call close_file('save',FN=FN_gnu)
   
  endif ! (numpar%gnupl%do_gnuplot) then ! do only if user wants plots
end subroutine create_Range_gnuplot



!===================================================
! Printing out mean free paths:
subroutine printout_MFPs(used_target, numpar, bunch)
   type(Matter), intent(in), target :: used_target	! parameters of the target
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   !-----------------------------------------------------
   character(250) :: Path, command, File_name2, File_EMFL, File_Brems, File_annihil, File_pair, File_Compton, File_Rayleigh
   character(250) :: File_for_gnu2, File_EMFL_gnu, File_Brems_gnu, File_annihil_gnu, File_pair_gnu, File_Compton_gnu, File_Rayleigh_gnu
   character(250), dimension(:), allocatable :: File_name, File_for_gnu
   character(50) :: Model_name, mass_model
   character(5) :: Ion_name
   real(8) :: x_start, x_end, y_start, y_end, eps
   integer :: i, j, k, m, FN2, N_elements, Nsiz, N_shells, N_types_SHI, N_bunch, ipart, ibunch, iret 
   integer, dimension(:), allocatable :: FN
   logical :: file_exist
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   
   if (numpar%printout_MFPs) then    ! do it only if a user requested it
      path_sep => numpar%path_sep
      ! Prepare file with the data output:
      TRGT:do i = 1, used_target%NOC ! for all targets
         ! Create a directory for MFPs:
         Path = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_folder_MFP))//trim(adjustl(used_target%Material(i)%Name))
         inquire(DIRECTORY=trim(adjustl(Path)),exist=file_exist)    ! check if input file excists
         if (.not.file_exist) then ! to make sure that such a folder is present (even if empty)
            ! Create a new directory for output files:
            command='mkdir '//trim(adjustl(Path))  ! to create a folder use this command
#ifdef _OPENMP
            iret = system(command)   ! create the folder
#else
            call system(command)   ! create the folder
#endif
         endif
         
         ! make sure names of the files with elements are available:
         if (allocated(File_name)) deallocate(File_name)
         if (allocated(File_for_gnu)) deallocate(File_for_gnu)
         if (allocated(FN)) deallocate(FN)
         N_elements = size(used_target%Material(i)%Elements)	! that's how many different elements are in this target
         allocate(File_name(N_elements))
         allocate(File_for_gnu(N_elements))
         allocate(FN(N_elements))
         !=========================
         ! Create the output files:
         
         !-------------------------------------------------
         ! 1) Photons:
         x_start = find_order_of_number(used_target%Material(i)%Ph_absorption_total%E(1)) ! module "Little_subroutines"
         x_end = 10.0**dble(-1 + find_order_of_number(used_target%Material(i)%Ph_absorption_total%E(size(used_target%Material(i)%Ph_absorption_total%E))) ) ! module "Little_subroutines"
         y_start = 10.0d0
         y_end = 1.0d12
         
         ! Absorption:
         select case (numpar%Ph_absorb) ! Include photoabsorption
         case (1)
            Model_name = '_CDF'
         case (2)
            Model_name = '_EPDL'
         case default	! no Compton effect included
            Model_name = '_NO'
         end select
         do j =1, N_elements	! for each element
            Element => used_target%Material(i)%Elements(j)	! all information about this element
            File_for_gnu(j) = trim(adjustl(m_output_absorb))//trim(adjustl(used_target%Material(i)%Name))//'_'//trim(adjustl(Element%Name))//'_photon'//trim(adjustl(Model_name))//'.dat'
            File_name(j) = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu(j)))
            open(newunit = FN(j), FILE = trim(adjustl(File_name(j))))
            ! Create the header:
            write(FN(j),'(a)', advance='no') '#Energy(eV) '
            call write_shell_headers_output(FN(j), used_target%Material(i), Element, numpar=numpar)  ! below        
         enddo
          ! Save core-shells data into the file for photons:
         do j =1, N_elements	! for each element
            Element => used_target%Material(i)%Elements(j)	! all information about this element
            Nsiz = size(Element%Phot_absorption%E)
            do m = 1, Nsiz	! for all energy grid points:
               E => Element%Phot_absorption%E(m)    ! photons energy [eV]
               write(FN(j),'(es)', advance='no') E   ! energy
               N_shells = Element%N_shl
               ! Calculate total cross sections for all shells of this element:
               do k = 1, N_shells
                  ! Check valence band
                  if (Element%valent(k) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band
                        ! Valence band will be added at the end
                  else ! core shell
                     write(FN(j),'(es)', advance='no') Element%Phot_absorption%MFP(k,m) ! MFP [A]
                  endif
               enddo ! k = 1, N_shells
               write(FN(j),'(a)') ''  ! to go to the next line
            enddo ! m = 1, Nsiz
            call close_file('save', FN=FN(j))  ! module "Dealing_with_files"
         enddo
         
         ! Write valence and total data for photons:
         File_for_gnu2 = trim(adjustl(m_output_absorb))//trim(adjustl(used_target%Material(i)%Name))//'_photon'//trim(adjustl(Model_name))//'_total.dat'
         File_name2 = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu2))
         open(newunit = FN2, FILE = trim(adjustl(File_name2)))
         write(FN2,'(a)', advance='no') '#Energy(eV) '
         if ( (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band
            write(FN2,'(a)', advance='no') 'Valence '
         endif
         write(FN2,'(a)') 'Total'
         Nsiz = size(used_target%Material(i)%Ph_absorption_valent%E)   ! valence band has a different grid from core shells
         do m = 1, Nsiz	! for all energy grid points:
            E => used_target%Material(i)%Ph_absorption_valent%E(m)	! photons energy [eV]
            write(FN2,'(es)', advance='no') E   ! energy
            if ( (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band
               write(FN2,'(es)', advance='no') used_target%Material(i)%Ph_absorption_valent%Total_MFP(m)   ! MFP [A]
            endif
            write(FN2,'(es)') used_target%Material(i)%Ph_absorption_total%Total_MFP(m)   ! total MFP [A]
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         
         ! Compton scattering for photons:
         select case (numpar%Ph_Compton)
         case (1)	! Include Compton:
            Model_name = '_PENELOPE'
         case default	! no Compton effect included
            Model_name = '_NO'
         end select
         File_Compton_gnu = trim(adjustl(m_output_Compton))//trim(adjustl(used_target%Material(i)%Name))//'_photon'//trim(adjustl(Model_name))//'_total.dat'
         File_Compton = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_Compton_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_Compton)))
         write(FN2,'(a)') '#Energy(eV)    Compton_MFP(A)'
         Nsiz = size(used_target%Material(i)%Ph_Compton_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%Ph_Compton_total%E(m), used_target%Material(i)%Ph_Compton_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         
         ! Coherent (Thomson or Rayleigh) scattering of photons:
         select case (numpar%Ph_Thomson)
         case (1)	! Include Rayleigh:
            Model_name = '_PENELOPE'
         case default	! no Compton effect included
            Model_name = '_NO'
         end select
         File_Rayleigh_gnu = trim(adjustl(m_output_Rayleigh))//trim(adjustl(used_target%Material(i)%Name))//'_photon'//trim(adjustl(Model_name))//'_total.dat'
         File_Rayleigh = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_Rayleigh_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_Rayleigh)))
         write(FN2,'(a)') '#Energy(eV)    Rayleigh_MFP(A)'
         Nsiz = size(used_target%Material(i)%Ph_coherent_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%Ph_coherent_total%E(m), used_target%Material(i)%Ph_coherent_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         
         ! Pair creation by photons:
         select case (numpar%Ph_Pairs)
         case (1)	! Include Compton:
            Model_name = '_BH'
         case default	! no Compton effect included
            Model_name = '_NO'
         end select
         File_pair_gnu= trim(adjustl(m_output_pair))//trim(adjustl(used_target%Material(i)%Name))//'_photon'//trim(adjustl(Model_name))//'_total.dat'
         File_pair = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_pair_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_pair)))
         write(FN2,'(a)') '#Energy(eV)    Pair_creation_MFP(A)'
         Nsiz = size(used_target%Material(i)%Ph_pair_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%Ph_pair_total%E(m), used_target%Material(i)%Ph_pair_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
   
         
         ! Create the plot with MFPs of photons:
         call create_MFPs_gnuplot(File_for_gnu, File_for_gnu2, used_target%Material(i), x_start, x_end, y_start, y_end, &
                'photon', numpar=numpar, File_Compton=File_Compton_gnu, File_Rayleigh=File_Rayleigh_gnu, File_pair=File_pair_gnu)    ! below
         !-------------------------------------------------
         ! 2) Electrons:
         
         ! Define axes for the plot:
         x_start = 10.0**dble( find_order_of_number(used_target%Material(i)%El_elastic_total%E(1)) )   ! module "Little_subroutines"
         x_end = 10.0**dble(-1 + find_order_of_number(used_target%Material(i)%El_inelastic_total%E(size(used_target%Material(i)%El_inelastic_total%E))) ) ! module "Little_subroutines"
         y_start = 1.0d0
         y_end = 1.0d6
         
         ! Inelastic:
         select case (numpar%El_inelast)
         case (1)	! CDF with Ritchi's oscillators
            Model_name = '_CDF_Ritchie'
         case (2)	! RBEB
            Model_name = '_RBEB'
         case (3)	! CDF with extended Liljequist model (delta-functions)
            Model_name = '_CDF_delta'
         case (4)	! CDF with Ritchi's oscillators
            Model_name = '_CDF_nonrel'
         case (5)	! CDF with Ritchi's oscillators
            Model_name = '_CDF_SPdelta'
         case default	! exclude
            Model_name = '_NO'
         end select
         
         ! Prepare the files for each element to be written into:
         N_elements = size(used_target%Material(i)%Elements)	! that's how many different elements are in this target
         LMNT:do j =1, N_elements	! for each element
            Element => used_target%Material(i)%Elements(j)	! all information about this element
            File_for_gnu(j) = trim(adjustl(m_output_IMFP))//trim(adjustl(used_target%Material(i)%Name))//'_'//trim(adjustl(Element%Name))//'_electron'//trim(adjustl(Model_name))//'.dat'
            File_name(j) = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu(j)))
            open(newunit = FN(j), FILE = trim(adjustl(File_name(j))))
            ! Create the header:
            write(FN(j),'(a)', advance='no') '#Energy(eV) '
            call write_shell_headers_output(FN(j), used_target%Material(i), Element, numpar=numpar)  ! below        
         enddo LMNT
         
         ! Save core-shells data into the file for electrons:
         LMNT2:do j =1, N_elements	! for each element
            Element => used_target%Material(i)%Elements(j)	! all information about this element
            Nsiz = size(Element%Phot_absorption%E)
            do m = 1, Nsiz	! for all energy grid points:
               E => Element%El_inelastic%E(m)    ! electron energy [eV]
               write(FN(j),'(es)', advance='no') E   ! energy
               N_shells = Element%N_shl
               ! Calculate total cross sections for all shells of this element:
               do k = 1, N_shells
                  ! Check valence band
                  VAL:if ((numpar%El_inelast /= 2) .and. (Element%valent(k)) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
                        ! Valence band will be added at the end
                  else VAL    ! core shell
                     write(FN(j),'(es)', advance='no') Element%El_inelastic%MFP(k,m) ! MFP [A]
                  endif VAL
               enddo ! k = 1, N_shells
               write(FN(j),'(a)') ''  ! to go to the next line
            enddo ! m = 1, Nsiz
            call close_file('save', FN=FN(j))  ! module "Dealing_with_files"
         enddo LMNT2
         
         ! Write valence and total data for electrons:
         File_for_gnu2 = trim(adjustl(m_output_IMFP))//trim(adjustl(used_target%Material(i)%Name))//'_electron'//trim(adjustl(Model_name))//'_total.dat'
         File_name2 = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu2))
         open(newunit = FN2, FILE = trim(adjustl(File_name2)))
         write(FN2,'(a)', advance='no') '#Energy(eV) '
         if ((numpar%El_inelast /= 2) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
            write(FN2,'(a)', advance='no') 'Valence '
         endif
         write(FN2,'(a)') 'Total'
         Nsiz = size(used_target%Material(i)%Ph_absorption_valent%E)   ! valence band has a different grid from core shells
         do m = 1, Nsiz	! for all energy grid points:
            E => used_target%Material(i)%El_inelastic_total%E(m)	! electron energy [eV]
            write(FN2,'(es)', advance='no') E   ! energy
            if ((numpar%El_inelast /= 2) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
               write(FN2,'(es)', advance='no') used_target%Material(i)%El_inelastic_valent%Total_MFP(m)   ! MFP [A]
            endif
            write(FN2,'(es)') used_target%Material(i)%El_inelastic_total%Total_MFP(m)   ! total MFP [A]
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         
         ! Elastic electron scattering MFP:
         select case (numpar%El_elastic)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
         case (1)	! CDF
            Model_name = '_CDF'
         case (2)	! Mott
            Model_name = '_Mott'
         case (3)	! DSF
            Model_name = '_DSF'
         case (5)	! CDF
            Model_name = '_CDF_SP'
         case default	! exclude
            Model_name = '_NO'
         end select
         File_EMFL_gnu = trim(adjustl(m_output_EMFP))//trim(adjustl(used_target%Material(i)%Name))//'_electron'//trim(adjustl(Model_name))//'_total.dat'
         File_EMFL = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_EMFL_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_EMFL)))
         write(FN2,'(a)') '#Energy(eV)    Elastic_MFP(A)'
         Nsiz = size(used_target%Material(i)%El_elastic_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%El_elastic_total%E(m), used_target%Material(i)%El_elastic_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         
         ! Electron Bremsstrahlung MFP:
         select case (numpar%El_Brems)	! elastic scattering: 0=excluded, 1=BHW
         case (1)	! BHW
            Model_name = '_BHW'
         case default	! excluded
            Model_name = '_NO'
         end select
         File_Brems_gnu = trim(adjustl(m_output_Brems))//trim(adjustl(used_target%Material(i)%Name))//'_electron'//trim(adjustl(Model_name))//'.dat'
         File_Brems = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_Brems_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_Brems)))
         write(FN2,'(a)') '#Energy(eV)    Brems_MFP(A)'
         Nsiz = size(used_target%Material(i)%El_Brems_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%El_Brems_total%E(m), used_target%Material(i)%El_Brems_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         

         ! Create the plot with MFPs of electrons:
         call create_MFPs_gnuplot(File_for_gnu, File_for_gnu2, used_target%Material(i), x_start, x_end, y_start, y_end, &
                'electron', numpar=numpar, File_EMFL=File_EMFL_gnu, File_Brems=File_Brems_gnu)    ! below
         
         !------------------------------------------------
         ! 3) Positrons:
         ! Define axes for the plot:
!          x_start = find_order_of_number(used_target%Material(i)%Pos_inelastic_total%E(1)) ! module "Little_subroutines"
         x_start = 10.0**dble( find_order_of_number(used_target%Material(i)%Pos_elastic_total%E(1)) )   ! module "Little_subroutines"
         x_end = 10.0**dble(-1 + find_order_of_number(used_target%Material(i)%Pos_inelastic_total%E(size(used_target%Material(i)%Pos_inelastic_total%E))) ) ! module "Little_subroutines"
         y_start = 1.0d0
         y_end = 1.0d7
         ! Inelastic
         select case (numpar%Pos_inelast)
         case (1:3)	! CDF with delta-functions
            Model_name = '_CDF_delta'
         case (5)	! CDF with delta-functions
            Model_name = '_CDF_SPdelta'
         case default	! exclude
            Model_name = '_NO'
         end select
         ! Prepare the files for each element to be written into:
         do j =1, N_elements	! for each element
            Element => used_target%Material(i)%Elements(j)	! all information about this element
            File_for_gnu(j) = trim(adjustl(m_output_IMFP))//trim(adjustl(used_target%Material(i)%Name))//'_'//trim(adjustl(Element%Name))//'_positron'//trim(adjustl(Model_name))//'.dat'
            File_name(j) = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu(j)))
            open(newunit = FN(j), FILE = trim(adjustl(File_name(j))))
            ! Create the header:
            write(FN(j),'(a)', advance='no') '#Energy(eV) '
            call write_shell_headers_output(FN(j), used_target%Material(i), Element, numpar=numpar)  ! below        
         enddo
         ! Save core-shells data into the file:
         do j =1, N_elements	! for each element
            Element => used_target%Material(i)%Elements(j)	! all information about this element
            Nsiz = size(Element%Phot_absorption%E)
            do m = 1, Nsiz	! for all energy grid points:
               E => Element%Pos_inelastic%E(m)    ! electron energy [eV]
               write(FN(j),'(es)', advance='no') E   ! energy
               N_shells = Element%N_shl
               ! Calculate total cross sections for all shells of this element:
               do k = 1, N_shells
                  ! Check valence band
                  if ((numpar%Pos_inelast /= 2) .and. (Element%valent(k)) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
                        ! Valence band will be added at the end
                  else    ! core shell
                     write(FN(j),'(es)', advance='no') Element%Pos_inelastic%MFP(k,m) ! MFP [A]
                  endif
               enddo ! k = 1, N_shells
               write(FN(j),'(a)') ''  ! to go to the next line
            enddo ! m = 1, Nsiz
            call close_file('save', FN=FN(j))  ! module "Dealing_with_files"
         enddo
         ! Write valence and total data for positions:
         File_for_gnu2 = trim(adjustl(m_output_IMFP))//trim(adjustl(used_target%Material(i)%Name))//'_positron'//trim(adjustl(Model_name))//'_total.dat'
         File_name2 = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu2))
         open(newunit = FN2, FILE = trim(adjustl(File_name2)))
         write(FN2,'(a)', advance='no') '#Energy(eV) '
         if ((numpar%Pos_inelast /= 2) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
            write(FN2,'(a)', advance='no') 'Valence '
         endif
         write(FN2,'(a)') 'Total'
         Nsiz = size(used_target%Material(i)%Pos_inelastic_total%E)   ! valence band has a different grid from core shells
         do m = 1, Nsiz	! for all energy grid points:
            E => used_target%Material(i)%Pos_inelastic_total%E(m)	! positron energy [eV]
            write(FN2,'(es)', advance='no') E   ! energy
            if ((numpar%Pos_inelast /= 2) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
               write(FN2,'(es)', advance='no') used_target%Material(i)%Pos_inelastic_valent%Total_MFP(m)   ! MFP [A]
            endif
            write(FN2,'(es)') used_target%Material(i)%Pos_inelastic_total%Total_MFP(m)   ! total MFP [A]
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         
         ! Elastic positron scattering MFP:
         select case (numpar%Pos_elastic)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
         case (1)	! CDF
            Model_name = '_CDF'
         case (2)	! Mott
            Model_name = '_Mott'
         case (3)	! DSF
            Model_name = '_DSF'
         case default	! exclude
            Model_name = '_NO'
         end select
         File_EMFL_gnu = trim(adjustl(m_output_EMFP))//trim(adjustl(used_target%Material(i)%Name))//'_positron'//trim(adjustl(Model_name))//'_total.dat'
         File_EMFL = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_EMFL_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_EMFL)))
         write(FN2,'(a)') '#Energy(eV)    Elastic_MFP(A)'
         Nsiz = size(used_target%Material(i)%Pos_elastic_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%Pos_elastic_total%E(m), used_target%Material(i)%Pos_elastic_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         
         ! Positron Bremsstrahlung MFP:
         select case (numpar%Pos_Brems)	! elastic scattering: 0=excluded, 1=BHW
         case (1)	! BHW
            Model_name = '_BHW'
         case default	! excluded
            Model_name = '_NO'
         end select
         File_Brems_gnu = trim(adjustl(m_output_Brems))//trim(adjustl(used_target%Material(i)%Name))//'_positron'//trim(adjustl(Model_name))//'.dat'
         File_Brems = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_Brems_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_Brems)))
         write(FN2,'(a)') '#Energy(eV)    Brems_MFP(A)'
         Nsiz = size(used_target%Material(i)%Pos_Brems_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%Pos_Brems_total%E(m), used_target%Material(i)%Pos_Brems_total%Total_MFP(m)
!             write(FN2,'(es,es)')  used_target%Material(i)%El_Brems_total%E(m), used_target%Material(i)%El_Brems_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
              
         ! Positron annihilation MFP:
         select case (numpar%Pos_annih)
         case (1)	! Heitler model
            Model_name = '_Heitler'
         case default	! exclude
            Model_name = '_NO'   
         end select
         File_annihil_gnu = trim(adjustl(m_output_annihil))//trim(adjustl(used_target%Material(i)%Name))//'_positron'//trim(adjustl(Model_name))//'.dat'
         File_annihil = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_annihil_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_annihil)))
         write(FN2,'(a)') '#Energy(eV)    Annihilation_MFP(A)'
         Nsiz = size(used_target%Material(i)%Pos_annihil_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%Pos_annihil_total%E(m), used_target%Material(i)%Pos_annihil_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         
         ! Create the plot with MFPs of positrons:
         call create_MFPs_gnuplot(File_for_gnu, File_for_gnu2, used_target%Material(i), x_start, x_end, y_start, y_end, &
                'positron', numpar=numpar, File_EMFL=File_EMFL_gnu, File_Brems=File_Brems_gnu, File_annihil=File_annihil_gnu)    ! below
         
         
         !------------------------------------------------
         ! 4) Ions:
         ! Check if there is any SHI to calcualte:
         N_types_SHI = count_types_of_SHIs(bunch) ! module "Initial_conditions"
         ! Only if there are some ions, calculate MFPs for them
         if (N_types_SHI > 0) then  ! there are ions
            ! Find how many sources of radiation are used:
            N_bunch = size (bunch)
            ipart = 0
            ! Check for all of them, whether there are ions (do we need to include ions here or not):
            BNCH:do ibunch = 1, N_bunch
               if (bunch(ibunch)%KOP == 3) then ! it is an ion
                  if (.not.repeated_type(bunch, ibunch)) then   ! new type of SHI, create its plots:
                     Ion_name = bunch(ibunch)%Name  ! this is the name of the SHI
                     ipart = ipart + 1
                  
                     ! Define axes for the plot:
                     x_start = (bunch(ibunch)%Z/100.0d0) * 1.0d-4 * 10.0**dble(find_order_of_number(used_target%Material(i)%SHI_inelastic_total(ipart)%E(1))) ! module "Little_subroutines"
                     x_end = (bunch(ibunch)%Z/100.0d0) *1.0d-6 * 10.0**dble(-2 + find_order_of_number(used_target%Material(i)%SHI_inelastic_total(ipart)%E( &
                                                 size(used_target%Material(i)%SHI_inelastic_total(ipart)%E))) ) ! module "Little_subroutines"
                     y_start = 0.01 * (100.0d0/bunch(ibunch)%Z)
                     y_end = 1.0d3 * (100.0d0/bunch(ibunch)%Z)
                     ! Inelastic
                     select case (numpar%SHI_inelast)    ! so far, we can only have CDF-based models for SHI
                     case (4)	! non-relativistic CDF with Ritchi's oscillators
                        Model_name = '_CDF_Ritchie'
                     case default	! delta-function CDF
                        Model_name = '_CDF'
                     end select
                     ! Prepare the files for each element to be written into:
                     do j =1, N_elements	! for each element
                        Element => used_target%Material(i)%Elements(j)	! all information about this element
                        File_for_gnu(j) = trim(adjustl(m_output_IMFP))//trim(adjustl(used_target%Material(i)%Name))//'_'//trim(adjustl(Element%Name)) &
                                              //'_'//trim(adjustl(Ion_name))//trim(adjustl(Model_name))//'.dat'
                        File_name(j) = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu(j)))
                        open(newunit = FN(j), FILE = trim(adjustl(File_name(j))))
                        ! Create the header:
                        write(FN(j),'(a)', advance='no') '#Energy(MeV) '
                        call write_shell_headers_output(FN(j), used_target%Material(i), Element, numpar=numpar)  ! below        
                     enddo
                     ! Save core-shells data into the file:
                     do j =1, N_elements	! for each element
                        Element => used_target%Material(i)%Elements(j)	! all information about this element
                        Nsiz = size(Element%SHI_inelastic(ipart)%E)
                        do m = 1, Nsiz	! for all energy grid points:
                           E => Element%SHI_inelastic(ipart)%E(m)    ! electron energy [eV]
                           write(FN(j),'(es)', advance='no') E*1.0d-6   ! energy [MeV]
                           N_shells = Element%N_shl
                           ! Calculate total cross sections for all shells of this element:
                           do k = 1, N_shells
                              ! Check valence band
                              if ((numpar%Pos_inelast /= 2) .and. (Element%valent(k)) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
                                 ! Valence band will be added at the end
                              else    ! core shell
                                 write(FN(j),'(es)', advance='no') Element%SHI_inelastic(ipart)%MFP(k,m) ! MFP [A]
                              endif
                           enddo ! k = 1, N_shells
                           write(FN(j),'(a)') ''  ! to go to the next line
                        enddo ! m = 1, Nsiz
                        call close_file('save', FN=FN(j))  ! module "Dealing_with_files"
                     enddo !  j =1, N_elements	! for each element
                     
                     ! Write valence and total data for positions:
                     File_for_gnu2 = trim(adjustl(m_output_IMFP))//trim(adjustl(used_target%Material(i)%Name)) &
                                          //'_'//trim(adjustl(Ion_name))//'_'//trim(adjustl(Model_name))//'_total.dat'
                     File_name2 = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu2))
                     open(newunit = FN2, FILE = trim(adjustl(File_name2)))
                     write(FN2,'(a)', advance='no') '#Energy(eV) '
                     if ((numpar%Pos_inelast /= 2) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
                        write(FN2,'(a)', advance='no') 'Valence '
                     endif
                     write(FN2,'(a)') 'Total'
                     Nsiz = size(used_target%Material(i)%SHI_inelastic_total(ipart)%E)   ! valence band has a different grid from core shells
                     do m = 1, Nsiz	! for all energy grid points:
                        E => used_target%Material(i)%SHI_inelastic_total(ipart)%E(m)	! positron energy [eV]
                        write(FN2,'(es)', advance='no') E*1.0d-6   ! energy [MeV]
                        if ((numpar%El_inelast /= 2) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
                           write(FN2,'(es)', advance='no') used_target%Material(i)%SHI_inelastic_valent(ipart)%Total_MFP(m)   ! MFP [A]
                        endif
                        write(FN2,'(es)') used_target%Material(i)%SHI_inelastic_total(ipart)%Total_MFP(m)   ! total MFP [A]
                     enddo ! m = 1, Nsiz
                     call close_file('save', FN=FN2)  ! module "Dealing_with_files"
                  
                     ! Create the plot with MFPs of SHI:
                     call create_MFPs_gnuplot(File_for_gnu, File_for_gnu2, used_target%Material(i), x_start, x_end, y_start, y_end, &
                           Ion_name, numpar=numpar, ch_units='(MeV)')    ! below
                  endif ! (.not.repeated_type(bunch, ibunch))
               endif ! (bunch(ibunch)%KOP == 3) then ! it is an ion
            enddo BNCH
         endif ! (N_types_SHI > 0) then  ! there are ions


         !------------------------------------------------
         ! 5) Muons:
         ! Define axes for the plot:
         x_start = 10.0**dble( find_order_of_number(used_target%Material(i)%Muon_elastic_total%E(1)) )   ! module "Little_subroutines"
         x_end = 10.0**dble(-1 + find_order_of_number(used_target%Material(i)%Muon_inelastic_total%E(size(used_target%Material(i)%Muon_inelastic_total%E))) ) ! module "Little_subroutines"
         y_start = 1.0d0
         y_end = 1.0d7
         ! Inelastic
         select case (numpar%Mu_inelast)
         case (1:3)	! CDF with delta-functions
            Model_name = '_CDF_delta'
         case (5)	! CDF with delta-functions
            Model_name = '_CDF_SPdelta'
         case default	! exclude
            Model_name = '_NO'
         end select
         ! Prepare the files for each element to be written into:
         do j =1, N_elements	! for each element
            Element => used_target%Material(i)%Elements(j)	! all information about this element
            File_for_gnu(j) = trim(adjustl(m_output_IMFP))//trim(adjustl(used_target%Material(i)%Name))//'_'//trim(adjustl(Element%Name))//'_muon'//trim(adjustl(Model_name))//'.dat'
            File_name(j) = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu(j)))
            open(newunit = FN(j), FILE = trim(adjustl(File_name(j))))
            ! Create the header:
            write(FN(j),'(a)', advance='no') '#Energy(eV) '
            call write_shell_headers_output(FN(j), used_target%Material(i), Element, numpar=numpar)  ! below
         enddo
         ! Save core-shells data into the file:
         do j =1, N_elements	! for each element
            Element => used_target%Material(i)%Elements(j)	! all information about this element
            Nsiz = size(Element%Phot_absorption%E)
            do m = 1, Nsiz	! for all energy grid points:
               E => Element%Muon_inelastic%E(m)    ! electron energy [eV]
               write(FN(j),'(es)', advance='no') E   ! energy
               N_shells = Element%N_shl
               ! Calculate total cross sections for all shells of this element:
               do k = 1, N_shells
                  ! Check valence band
                  if ((numpar%Mu_inelast /= 2) .and. (Element%valent(k)) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
                        ! Valence band will be added at the end
                  else    ! core shell
                     write(FN(j),'(es)', advance='no') Element%Muon_inelastic%MFP(k,m) ! MFP [A]
                  endif
               enddo ! k = 1, N_shells
               write(FN(j),'(a)') ''  ! to go to the next line
            enddo ! m = 1, Nsiz
            call close_file('save', FN=FN(j))  ! module "Dealing_with_files"
         enddo
         ! Write valence and total data for positions:
         File_for_gnu2 = trim(adjustl(m_output_IMFP))//trim(adjustl(used_target%Material(i)%Name))//'_muon'//trim(adjustl(Model_name))//'_total.dat'
         File_name2 = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu2))
         open(newunit = FN2, FILE = trim(adjustl(File_name2)))
         write(FN2,'(a)', advance='no') '#Energy(eV) '
         if ((numpar%Mu_inelast /= 2) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
            write(FN2,'(a)', advance='no') 'Valence '
         endif
         write(FN2,'(a)') 'Total'
         Nsiz = size(used_target%Material(i)%Muon_inelastic_total%E)   ! valence band has a different grid from core shells
         do m = 1, Nsiz	! for all energy grid points:
            E => used_target%Material(i)%Muon_inelastic_total%E(m)	! muon energy [eV]
            write(FN2,'(es)', advance='no') E   ! energy
            if ((numpar%Mu_inelast /= 2) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
               write(FN2,'(es)', advance='no') used_target%Material(i)%Muon_inelastic_valent%Total_MFP(m)   ! MFP [A]
            endif
            write(FN2,'(es)') used_target%Material(i)%Muon_inelastic_total%Total_MFP(m)   ! total MFP [A]
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"

         ! Elastic muon scattering MFP:
         select case (numpar%Mu_elastic)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
         case (1)	! CDF
            Model_name = '_CDF'
         case (2)	! Mott
            Model_name = '_Mott'
         case (3)	! DSF
            Model_name = '_DSF'
         case default	! exclude
            Model_name = '_NO'
         end select
         File_EMFL_gnu = trim(adjustl(m_output_EMFP))//trim(adjustl(used_target%Material(i)%Name))//'_muon'//trim(adjustl(Model_name))//'_total.dat'
         File_EMFL = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_EMFL_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_EMFL)))
         write(FN2,'(a)') '#Energy(eV)    Elastic_MFP(A)'
         Nsiz = size(used_target%Material(i)%Muon_elastic_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%Muon_elastic_total%E(m), used_target%Material(i)%Muon_elastic_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"

         ! Muon Bremsstrahlung MFP:
         select case (numpar%Mu_Brems)	! elastic scattering: 0=excluded, 1=BHW
         case (1)	! BHW
            Model_name = '_BHW'
         case default	! excluded
            Model_name = '_NO'
         end select
         File_Brems_gnu = trim(adjustl(m_output_Brems))//trim(adjustl(used_target%Material(i)%Name))//'_muon'//trim(adjustl(Model_name))//'.dat'
         File_Brems = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_Brems_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_Brems)))
         write(FN2,'(a)') '#Energy(eV)    Brems_MFP(A)'
         Nsiz = size(used_target%Material(i)%Muon_Brems_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%Muon_Brems_total%E(m), used_target%Material(i)%Muon_Brems_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"

         ! Muon Pair production MFP:
         ! NOT READY

         ! Create the plot with MFPs of muons:
         call create_MFPs_gnuplot(File_for_gnu, File_for_gnu2, used_target%Material(i), x_start, x_end, y_start, y_end, &
                'muon', numpar=numpar, File_EMFL=File_EMFL_gnu, File_Brems=File_Brems_gnu)    ! below

         !------------------------------------------------
         ! 6) Valence holes (MUST BE THE LAST TYPE OF PARTICLES; PLACE ALL OTHERS ABOVE):
         ! Define axes for the plot:
         x_start = 0.0d0
         if (size(used_target%Material(i)%H_inelastic_total%E) > 0) then
            x_end = used_target%Material(i)%H_inelastic_total%E(size(used_target%Material(i)%H_inelastic_total%E))  ! module "Little_subroutines"
         else
            x_end = 1.0d0
         endif
         y_start = 1.0d0
         y_end = 1.0d4
         eps = 1.0d-6 ! margin within which mass equals to zero
         if (numpar%H_m_eff > eps) then ! a constant coefficient * me
            if (abs(numpar%H_m_eff) < 1d6) then
               write(mass_model,'(f6.2)') numpar%H_m_eff
               mass_model = '_meff_'//trim(adjustl(mass_model))//'_'
            else  ! hole is too heave to move
               mass_model = '_frozen_'
            endif
         elseif (abs(numpar%H_m_eff) < eps) then ! equals to free electron mass
            mass_model = '_free_'
         else  ! effective mass from DOS
            mass_model = '_DOS_'   
         endif
         ! Inelastic
         select case (numpar%H_inelast)
         case (1,3)	! CDF with delta-functions
            Model_name = trim(adjustl(mass_model))//'CDF_Ritchie'
         case (5)	! CDF with delta-functions
            Model_name = trim(adjustl(mass_model))//'CDF_SPdelta'
         case default	! exclude
            Model_name = trim(adjustl(mass_model))//'NO'
         end select
         ! Write valence and total data for positions:
         File_for_gnu2 = trim(adjustl(m_output_IMFP))//trim(adjustl(used_target%Material(i)%Name))//'_hole'//trim(adjustl(Model_name))//'.dat'
         File_name2 = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_for_gnu2))
         open(newunit = FN2, FILE = trim(adjustl(File_name2)))
         write(FN2,'(a)', advance='no') '#Energy(eV) '
         if ((numpar%Pos_inelast /= 2) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
            write(FN2,'(a)', advance='no') 'Valence '
         endif
         write(FN2,'(a)') 'Total'
         Nsiz = size(used_target%Material(i)%H_inelastic_total%E)   ! valence band has a different grid from core shells
         do m = 1, Nsiz	! for all energy grid points:
            E => used_target%Material(i)%H_inelastic_total%E(m)	! valence hole energy [eV]
            write(FN2,'(es)', advance='no') E   ! energy
            if ((numpar%El_inelast /= 2) .and. (allocated(used_target%Material(i)%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
               write(FN2,'(es)', advance='no') used_target%Material(i)%H_inelastic_total%Total_MFP(m)   ! MFP [A]
            endif
            write(FN2,'(es)') used_target%Material(i)%H_inelastic_total%Total_MFP(m)   ! total MFP [A]
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         
         ! Elastic hole scattering MFP:
         select case (numpar%H_elast)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
         case (1)	! CDF
            Model_name =  trim(adjustl(mass_model))//'CDF'
         case (2)	! Mott
            Model_name =  trim(adjustl(mass_model))//'Mott'
         case (3)	! DSF
            Model_name =  trim(adjustl(mass_model))//'DSF'
         case default	! exclude
            Model_name =  trim(adjustl(mass_model))//'NO'
         end select
         File_EMFL_gnu = trim(adjustl(m_output_EMFP))//trim(adjustl(used_target%Material(i)%Name))//'_hole'//trim(adjustl(Model_name))//'_total.dat'
         File_EMFL = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(File_EMFL_gnu))
         open(newunit = FN2, FILE = trim(adjustl(File_EMFL)))
         write(FN2,'(a)') '#Energy(eV)    Elastic_MFP(A)'
         Nsiz = size(used_target%Material(i)%H_elastic_total%E)
         do m = 1, Nsiz     ! for all energy grid points:
            write(FN2,'(es,es)')  used_target%Material(i)%H_elastic_total%E(m), used_target%Material(i)%H_elastic_total%Total_MFP(m)
         enddo
         call close_file('save', FN=FN2)  ! module "Dealing_with_files"
         
         ! Create the plot with MFPs of valence holes:
         if (allocated(File_for_gnu)) deallocate(File_for_gnu)    ! no core shells ionization by valence hole is possible
         call create_MFPs_gnuplot(File_for_gnu, File_for_gnu2, used_target%Material(i), x_start, x_end, y_start, y_end, &
                'hole', numpar=numpar, File_EMFL=File_EMFL_gnu, logx_in=.false., logy_in=.true.)    ! below
         !------------------------------------------------

      enddo TRGT ! for all targets
   endif
   
   nullify(E, path_sep, Element)
end subroutine printout_MFPs




subroutine create_MFPs_gnuplot(File_name, File_name2, Material, x_start, x_end, y_start, y_end, particle_name, numpar, File_EMFL, File_Brems, File_annihil, File_pair, File_Compton, File_Rayleigh, logx_in, logy_in, ch_units)    ! below
   character(*), dimension(:), allocatable, intent(in) :: File_name      ! file with the core data to plot from
   character(*), intent(in) :: File_name2    ! file with the total data to plot from
   type(Target_atoms), intent(in), target :: Material ! parameters of this material
   character(*), intent(in) :: particle_name  ! name of the model and particle to include into the file names
   real(8), intent(in) :: x_start, x_end, y_start, y_end
   type(Num_par), intent(in), optional :: numpar	! all numerical parameters
   character(*), intent(in), optional :: File_EMFL, File_Brems, File_annihil, File_pair, File_Compton, File_Rayleigh ! files with the data to plot from: elastic, bremsstrahlung, annihilation, pair creation, Compton, Rayleigh
   logical, intent(in), optional :: logx_in, logy_in
   character(*), intent(in), optional :: ch_units
   !------------------------------------
   type(Atom_kind), pointer  :: Element
   real(8) :: tics
   integer :: j, k, N_elements, N_shells, FN_gnu, FN_eps, i
   character(200) :: File_script, Out_file, Path
   character(50) :: Title
   character(10) :: units
   character(5) ::  call_slash, sh_cmd, col_y
   logical :: logx, logy, first_line
   
  if (numpar%gnupl%do_gnuplot) then ! do only if user wants plots
   
   if (present(logx_in)) then
      logx = logx_in    ! follow what user set
   else
      logx = .true. ! set logscale x by default
   endif
   if (present(logy_in)) then
      logy = logy_in    ! follow what user set
   else
      logy = .true. ! set logscale y by default
   endif
   if (present(ch_units)) then
      units = ch_units  ! user provided units
   else
      units = '(eV)'    ! by default, assume eV
   endif

   ! Get the extension and slash in this OS:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"
   ! Get the paths and file names:
   Path = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_folder_MFP))//trim(adjustl(Material%Name))
   File_script = trim(adjustl(Path))//numpar%path_sep//trim(adjustl(m_output_MFP))//trim(adjustl(Material%Name))//'_'// &
                      trim(adjustl(particle_name))//trim(adjustl(sh_cmd))
   open(newunit = FN_gnu, FILE = trim(adjustl(File_script)))
   Out_file = trim(adjustl(m_output_MFP))//trim(adjustl(particle_name))//'_in_'//trim(adjustl(Material%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))
   
   ! Set the grid step on the plots:
   if (logx) then
      tics = 10.0d0
   else
      tics = dble(find_order_of_number( abs(x_start-x_end) )) ! module "Little_subroutines"
   endif
   
   ! Create the gnuplot-script header:
   call write_gnuplot_script_header_new(FN_gnu, 1, 3.0d0, tics, "MFP", "Energy "//trim(adjustl(units)), "Mean free path (A)", trim(adjustl(Out_file)), &
                trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, setkey=1, logx=logx, logy=logy)  ! module "Gnuplotting"
   
   if (allocated(File_name)) then
      N_elements = size(Material%Elements)	! that's how many different elements are in this target
      ! Core-shells:
      LMNT:do j =1, N_elements	! for each element
         Element => Material%Elements(j)	! all information about this element
         N_shells = Element%N_shl
         ! MFPs for all shells of this element:
         do k = 1, N_shells
            VAL:if ( (Element%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
                  ! Valence band will be added at the end
            else VAL    ! core shell
               write(col_y, '(i4)') 1 + k  ! number of column
               write(Title, '(a)') trim(adjustl(Element%Name))//' '//trim(adjustl(Element%Shell_name(k)))//'-shell inelastic'
               ! Create the plotting options:
               if ((j ==1) .and. (k==1)) then    ! first line in gnuplot script
                  if (numpar%path_sep == '\') then	! if it is Windows
                     call write_gnu_printout(FN_gnu, .true., .false., File_name(j), x_start=x_start, x_end=x_end, y_start=y_start, y_end=y_end, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
                  else
                     call write_gnu_printout(FN_gnu, .true., .false., File_name(j),  x_start=x_start, x_end=x_end, y_start=y_start, y_end=y_end, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
                  endif
               else
                  if (numpar%path_sep == '\') then	! if it is Windows
                     call write_gnu_printout(FN_gnu, .false., .false., File_name(j), col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
                  else
                     call write_gnu_printout(FN_gnu, .false., .false., File_name(j), col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
                  endif
               endif
            endif VAL
         enddo
      enddo LMNT
      first_line = .false.
   else ! Core shells are undefined, start from the valence band (e.g. for holes scattering)
      first_line = .true.
   endif
   ! Valence MFPs:
   i = 1
   if ( (allocated(Material%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
      i = i + 1 ! column with valence MFP
      write(col_y, '(i4)') i
      write(Title, '(a)') 'Valence inelastic'
      if (numpar%path_sep == '\') then	! if it is Windows
         call write_gnu_printout(FN_gnu, first_line, .false., File_name2, col_x="1", col_y=trim(adjustl(col_y)), &
            x_start=x_start, x_end=x_end, y_start=y_start, y_end=y_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
       else
          call write_gnu_printout(FN_gnu, first_line, .false., File_name2, col_x="1", col_y=trim(adjustl(col_y)), &
            x_start=x_start, x_end=x_end, y_start=y_start, y_end=y_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
       endif
       first_line = .false.
   endif
   
   if(present(File_EMFL)) then  ! plot also elastic MFP
      write(col_y, '(i4)') 2
      write(Title, '(a)') 'Elastic'
      if (numpar%path_sep == '\') then	! if it is Windows
         call write_gnu_printout(FN_gnu, first_line, .false., File_EMFL, col_x="1", col_y=trim(adjustl(col_y)), &
            x_start=x_start, x_end=x_end, y_start=y_start, y_end=y_end, lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
       else
          call write_gnu_printout(FN_gnu, first_line, .false., File_EMFL, col_x="1", col_y=trim(adjustl(col_y)), &
            x_start=x_start, x_end=x_end, y_start=y_start, y_end=y_end, lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
       endif
   endif ! present(File_EMFL)
   
   if(present(File_Brems)) then  ! plot also Bremsstrahlung
      write(col_y, '(i4)') 2
      write(Title, '(a)') 'Bremsstrahlung'
      if (numpar%path_sep == '\') then	! if it is Windows
         call write_gnu_printout(FN_gnu, .false., .false., File_Brems, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
       else
          call write_gnu_printout(FN_gnu, .false., .false., File_Brems, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
       endif
   endif ! present(File_Brems)
   
   if(present(File_annihil)) then  ! plot also Annihilation
      write(col_y, '(i4)') 2
      write(Title, '(a)') 'Annihilation'
      if (numpar%path_sep == '\') then	! if it is Windows
         call write_gnu_printout(FN_gnu, .false., .false., File_annihil, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
       else
          call write_gnu_printout(FN_gnu, .false., .false., File_annihil, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
       endif
   endif ! present(File_annihil)
   
   if(present(File_Rayleigh)) then  ! plot also Rayleigh
      write(col_y, '(i4)') 2
      write(Title, '(a)') 'Rayleigh'
      if (numpar%path_sep == '\') then	! if it is Windows
         call write_gnu_printout(FN_gnu, .false., .false., File_Rayleigh, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
       else
          call write_gnu_printout(FN_gnu, .false., .false., File_Rayleigh, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
       endif
   endif ! present(File_Rayleigh)
   
   if(present(File_Compton)) then  ! plot also Compton
      write(col_y, '(i4)') 2
      write(Title, '(a)') 'Compton'
      if (numpar%path_sep == '\') then	! if it is Windows
         call write_gnu_printout(FN_gnu, .false., .false., File_Compton, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
       else
          call write_gnu_printout(FN_gnu, .false., .false., File_Compton, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
       endif
   endif ! present(File_Compton)

   if(present(File_pair)) then  ! plot also e-e+ pair creation
      write(col_y, '(i4)') 2
      write(Title, '(a)') 'e-e+ pair creation'
      if (numpar%path_sep == '\') then	! if it is Windows
         call write_gnu_printout(FN_gnu, .false., .false., File_pair, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
       else
          call write_gnu_printout(FN_gnu, .false., .false., File_pair, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
       endif
   endif ! present(File_pair)
   
   
   ! Total inelastic MFPs:
   i = i + 1    ! column with total MFP
   write(col_y, '(i4)') i
   write(Title, '(a)') 'Total inelastic'
   if (numpar%path_sep == '\') then	! if it is Windows
      call write_gnu_printout(FN_gnu, .false., .true., File_name2, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)))  ! module "Gnuplotting"
   else
      call write_gnu_printout(FN_gnu, .false., .true., File_name2, col_x="1", col_y=trim(adjustl(col_y)), lw=3, title=trim(adjustl(Title)), linux_s=.true.)  ! module "Gnuplotting"
   endif
   
   ! Create the gnuplot-script ending:
   call  write_gnuplot_script_ending_new(FN_gnu, File_script, numpar%path_sep)  ! module "Gnuplotting"
   nullify(Element)
   
   call close_file('save', FN=FN_gnu)  ! module "Dealing_with_files"
   
  endif ! (g_numpar%gnupl%do_gnuplot)
end subroutine create_MFPs_gnuplot


subroutine write_shell_headers_output(FN, Material, Element, numpar, novalence)
   integer, intent(in) :: FN    ! file to write to
   type(Target_atoms), intent(in) :: Material ! parameters of this material
   type(Atom_kind), intent(in)  :: Element
   type(Num_par), intent(in), optional :: numpar	! all numerical parameters
   logical, optional :: novalence   ! if no need to include valence band, do all the shells instead
   integer :: k, N_elements, N_shells

   N_shells = Element%N_shl
   ! Calculate total cross sections for all shells of this element:
   do k = 1, N_shells
      ! Check valence band
      if (present(numpar)) then  ! model-dependent
         VAL:if ((numpar%El_inelast /= 2) .and.  (Element%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
            ! Valence band will be added at the end
         else VAL    ! core shell
            write(FN,'(a)', advance='no') trim(adjustl(Element%Name))//'_'//trim(adjustl(Element%Shell_name(k)))//'    '  
         endif VAL
      else ! not model-dependent
         VAL2:if ((Element%valent(k)) .and. (allocated(Material%CDF_valence%A))) then    ! Valence band
            ! Valence band will be added at the end
         else VAL2    ! core shell
            write(FN,'(a)', advance='no') trim(adjustl(Element%Name))//'_'//trim(adjustl(Element%Shell_name(k)))//'    '  
         endif VAL2   
      endif
   enddo
   write(FN,'(a)')
end subroutine write_shell_headers_output


!===================================================
! All subroutines for DOS output:
subroutine printout_DOS(used_target, numpar)
   type(Matter), intent(in) :: used_target	! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   character(250) :: File_name
   integer :: i, j, FN
   
   if (numpar%printout_DOS) then    ! do it only if a user requested it
      ! Prepare DOS file with the data:
      do i = 1, used_target%NOC
         ! For this material in the target:
         File_name = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_DOS))//trim(adjustl(used_target%Material(i)%Name))//'.dat'
         open(newunit = FN, FILE = trim(adjustl(File_name)))
      
         ! Titles:
         write(FN,'(a)') '#E    DOS k   me_eff'
         ! Data:
         do j = 1, size(used_target%Material(i)%DOS%E)
            write(FN, '(f, es, es, es)') used_target%Material(i)%DOS%E(j), used_target%Material(i)%DOS%DOS(j), used_target%Material(i)%DOS%k(j), used_target%Material(i)%DOS%Eff_m(j)
         enddo
      
         call close_file('save', FN=FN)  ! module "Dealing_with_files"
      enddo
      
      ! Plot the data:
      if (numpar%gnupl%do_gnuplot) then
         call gnuplot_DOS(used_target, numpar) ! below
         call gnuplot_DOS_k(used_target, numpar) ! below
         call gnuplot_DOS_m_eff(used_target, numpar) ! below
      endif
   endif
end subroutine printout_DOS


subroutine gnuplot_DOS(used_target, numpar)
   type(Matter), intent(in) :: used_target	! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   character(250) :: File_name, File_script, Out_file
   character(5) ::  call_slash, sh_cmd
   integer :: i, j, FN
   if (numpar%printout_DOS) then    ! do it only if a user requested it
      ! Get the extension and slash in this OS:
      call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"
      
      ! Prepare DOS file with the data:
      do i = 1, used_target%NOC
         ! Create gnuplot script file:
         !File_script = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_DOS_k))// &
         File_script = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_DOS))// &
                            trim(adjustl(used_target%Material(i)%Name))//trim(adjustl(sh_cmd))
         open(NEWUNIT=FN, FILE = trim(adjustl(File_script)), action="write", status="replace")
         Out_file = 'OUTPUT_DOS_in_'//trim(adjustl(used_target%Material(i)%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))
         
         ! For this material in the target:
         File_name = trim(adjustl(m_output_DOS))//trim(adjustl(used_target%Material(i)%Name))//'.dat'
         
         ! Create the gnuplot-script header:
         call write_gnuplot_script_header_new(FN, 1, 3.0d0, 2.0d0, "DOS", "Energy (eV)", "DOS (1/eV)", trim(adjustl(Out_file)), &
                    trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, setkey=0)  ! module "Gnuplotting"
         
         ! DOS:
         ! Create the plotting options:
         if (numpar%path_sep == '\') then	! if it is Windows
            call  write_gnu_printout(FN, .true., .true., File_name, col_x="1", col_y="2", lw=3, title='DOS')  ! module "Gnuplotting"
         else
            call  write_gnu_printout(FN, .true., .true., File_name, col_x="1", col_y="2", lw=3, title='DOS', linux_s=.true.)  ! module "Gnuplotting"
         endif
!          ! Effective mass:
!          ! Create the plotting options:
!          if (numpar%path_sep == '\') then	! if it is Windows
!             call  write_gnu_printout(FN, .false., .true., File_name, col_x="1", col_y="4", lw=3, title='Effective mass')  ! module "Gnuplotting"
!          else
!             call  write_gnu_printout(FN, .false., .true., File_name, col_x="1", col_y="4", lw=3, title='Effective mass', linux_s=.true.)  ! module "Gnuplotting"
!          endif
         
         ! Create the gnuplot-script ending:
         call  write_gnuplot_script_ending_new(FN, File_script, numpar%path_sep)  ! module "Gnuplotting"
         close(FN)
      enddo
   endif
end subroutine gnuplot_DOS


subroutine gnuplot_DOS_k(used_target, numpar)
   type(Matter), intent(in) :: used_target	! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   character(250) :: File_name, File_script, Out_file
   character(5) ::  call_slash, sh_cmd
   integer :: i, j, FN
   if (numpar%printout_DOS) then    ! do it only if a user requested it
      ! Get the extension and slash in this OS:
      call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"
      
      ! Prepare DOS file with the data:
      do i = 1, used_target%NOC
         ! Create gnuplot script file:
         !File_script = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_DOS_effm))//&
         File_script = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_DOS_k))//&
                            trim(adjustl(used_target%Material(i)%Name))//trim(adjustl(sh_cmd))
                            
         open(NEWUNIT=FN, FILE = trim(adjustl(File_script)), action="write", status="replace")
         Out_file = 'OUTPUT_DOS_k_vector_in_'//trim(adjustl(used_target%Material(i)%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))
         
         ! For this material in the target:
         File_name = trim(adjustl(m_output_DOS))//trim(adjustl(used_target%Material(i)%Name))//'.dat'
         
         ! Create the gnuplot-script header:
         call write_gnuplot_script_header_new(FN, 1, 3.0d0, 2.0d0, "k-vector", "Energy (eV)", "k-vector (1/m)", trim(adjustl(Out_file)), &
                    trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, setkey=0)  ! module "Gnuplotting"
         
         ! Create the plotting options:
         if (numpar%path_sep == '\') then	! if it is Windows
            call  write_gnu_printout(FN, .true., .true., File_name, col_x="1", col_y="3", lw=3, title='k-vector')  ! module "Gnuplotting"
         else
            call  write_gnu_printout(FN, .true., .true., File_name, col_x="1", col_y="3", lw=3, title='k-vector', linux_s=.true.)  ! module "Gnuplotting"
         endif
         
         ! Create the gnuplot-script ending:
         call  write_gnuplot_script_ending_new(FN, File_script, numpar%path_sep)  ! module "Gnuplotting"
         close(FN)
      enddo
   endif
end subroutine gnuplot_DOS_k



subroutine gnuplot_DOS_m_eff(used_target, numpar)
   type(Matter), intent(in) :: used_target	! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   character(250) :: File_name, File_script, Out_file
   character(5) ::  call_slash, sh_cmd
   integer :: i, j, FN
   if (numpar%printout_DOS) then    ! do it only if a user requested it
      ! Get the extension and slash in this OS:
      call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"
      
      ! Prepare DOS file with the data:
      do i = 1, used_target%NOC
         ! Create gnuplot script file:
         File_script = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_DOS_effm))// &
                            trim(adjustl(used_target%Material(i)%Name))//trim(adjustl(sh_cmd))
         open(NEWUNIT=FN, FILE = trim(adjustl(File_script)), action="write", status="replace")
         Out_file = 'OUTPUT_DOS_effective_mass_in_'//trim(adjustl(used_target%Material(i)%Name))//'.'//trim(adjustl(numpar%gnupl%gnu_extension))
         
         ! For this material in the target:
         File_name = trim(adjustl(m_output_DOS))//trim(adjustl(used_target%Material(i)%Name))//'.dat'
         
         ! Create the gnuplot-script header:
         call write_gnuplot_script_header_new(FN, 1, 3.0d0, 2.0d0, "Effective mass", "Energy (eV)", "Effective mass (me)", trim(adjustl(Out_file)), &
                    trim(adjustl(numpar%gnupl%gnu_terminal)), numpar%path_sep, setkey=0)  ! module "Gnuplotting"
         
         ! Create the plotting options:
         if (numpar%path_sep == '\') then	! if it is Windows
            call  write_gnu_printout(FN, .true., .true., File_name, col_x="1", col_y="4", y_end=5.0d0, lw=3, title='Effective mass')  ! module "Gnuplotting"
         else
            call  write_gnu_printout(FN, .true., .true., File_name, col_x="1", col_y="4", y_end=5.0d0, lw=3, title='Effective mass', linux_s=.true.)  ! module "Gnuplotting"
         endif
         
         ! Create the gnuplot-script ending:
         call  write_gnuplot_script_ending_new(FN, File_script, numpar%path_sep)  ! module "Gnuplotting"
         close(FN)
      enddo
   endif
end subroutine gnuplot_DOS_m_eff



!===================================================
! Gnuplotting all the scripts:
subroutine execute_all_gnuplots(numpar, out_path)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   character(*), intent(in), optional :: out_path    ! folder with the cmd-files
   character(200) :: File_name, command
   integer :: FN, N_f, i, n_slash
   integer :: open_status, iret, idir
   character(200), dimension(:), allocatable :: All_files
   character(300) :: output_path
   character(5) ::  call_slash, sh_cmd
   
   if (present(out_path)) then  ! if user provided the folder name:
      output_path = out_path
   else ! if not, by default use the main output folder:
      output_path = numpar%output_path
   endif
   
   ! Create a temporary file:
   File_name = trim(adjustl(output_path))//numpar%path_sep//'Temp.txt'   

   ! Find the extension of the gnuplot scripts:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"

   ! Save the names of all gnuplot scripts into this file:
   if (numpar%path_sep == '\') then	! if it is Windows
      command = 'dir '//trim(adjustl(output_path))//'\*'//trim(adjustl(sh_cmd))//' /b >'//trim(adjustl(File_name))
   else ! linux:
      command = "ls -t "//trim(adjustl(output_path))//" | grep '"//trim(adjustl(sh_cmd))//"' >"//trim(adjustl(File_name))
   endif
#ifdef _OPENMP
   iret = system(trim(adjustl(command)))   ! execute the command to save file names in the temp file
#else
   call system(trim(adjustl(command))) ! execute the command to save file names in the temp file
#endif
   
   ! Open the files with gnuplot script names:
   open(NEWUNIT=FN, file=trim(adjustl(File_name)), iostat=open_status, action='read')
   if ( open_status /= 0 ) then
      print *, 'Could not open ',trim(adjustl(File_name)),' for gnuplotting.', ' Unit = ', FN
   endif
   
   ! Find out how many there are:
   call Count_lines_in_file(FN, N_f) ! below
   
   ! Allocate array with them:
   allocate(All_files(N_f)) ! array with all relevant file names
   
   ! Read file names:
   do i = 1,N_f
      read(FN,*) All_files(i)
   enddo
   
   ! Delete the unneeded temporary file:
   call close_file('delete',FN=FN) ! temp file is not needed anymore, erase it
   
   ! Execute all the gnuplot scripts:
#ifdef _OPENMP
   idir = chdir(trim(adjustl(output_path))) ! go into the directory with output files
#else
   call chdir(trim(adjustl(output_path))) ! go into the directory with output files
#endif

   if (numpar%path_sep == '\') then	! if it is Windows
      do i = 1,N_f  ! execute all gnuplot scripts (for all targets)
#ifdef _OPENMP
         iret = system( '@echo off' )   ! create the folder
         iret = system(trim(adjustl(call_slash))//' '//trim(adjustl(All_files(i))))   ! create the folder
#else
         call system( '@echo off' )
         call system(trim(adjustl(call_slash))//' '//trim(adjustl(All_files(i))))
#endif
         
      enddo
   else ! linux:
      do i = 1,N_f  ! execute all gnuplot scripts (for all targets)
#ifdef _OPENMP
         iret = system( '#!/bin/bash' )
         iret = system(trim(adjustl(call_slash))//trim(adjustl(All_files(i))))
#else
         call system( '#!/bin/bash' )
         call system(trim(adjustl(call_slash))//trim(adjustl(All_files(i))))
#endif
         
      enddo
   endif
   
   ! Count how many times the system has to go out of the directory to get back into the original directory:
   ! Defined by the number of slashes in the path given:
   n_slash = count( (/ (trim(adjustl(output_path(i:i))), i=1,len_trim(output_path)) /) == trim(adjustl(numpar%path_sep)) )
   do i = 1, n_slash+1  ! go up in the directory tree as many times as needed
#ifdef _OPENMP
      idir = chdir("../")    ! exit the directory with output files
#else
      call chdir("../")    ! exit the directory with output files
#endif   
   enddo

end subroutine execute_all_gnuplots


!===================================================
! Creating output folder:
subroutine make_output_folder(used_target, numpar, bunch)
   type(Matter), intent(in) :: used_target	! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(Radiation_param), dimension(:), allocatable, intent(inout) :: bunch	! incomming radiation
   !--------------------------------------------------
   character(200) :: File_name, File_name2, command, bunch_name, NRG, temp, ch1, chtest, chtest1, error_message, Path
   integer :: i, Nsiz, iret, INFO
   real(8) :: E_reduced
   logical :: file_exist
   type(atomic_data), dimension(:), allocatable :: Periodic_table ! this is an internal module variable
   
   ! Output directory name will contain the target name:
   File_name = 'OUTPUT_'//trim(adjustl(used_target%Name))
   
   ! Output name will contain the details of radiation source:
   select case(bunch(1)%KOP)	! use the first radiation bunch to define the name
   case (0)  ! Kind of particle: 0=photon, 1=electron, 2=positron, 3=SHI, 4=hole, 5=muon
      bunch_name = 'photon'
   case (1)
      bunch_name = 'electron'
   case (2)
      bunch_name = 'positron'
   case (3)
      ! Read the periodic table, in case we need it:
      Path = trim(adjustl(m_input_folder))//numpar%path_sep//trim(adjustl(m_databases)) ! where to find the periodic table
      call Read_from_periodic_table(Path, numpar%path_sep, INFO, error_message, 1, Periodic_table_out=Periodic_table) ! module "Periodic_table"
      ! Get the name of the ion from its ztomic number:
      bunch(1)%Name = Periodic_table(bunch(1)%Z)%Name  ! name of the ion
      !print*, bunch(1)%Z
      !print*, 'Name: ', bunch(1)%Name
      bunch_name = trim(adjustl(bunch(1)%Name))//'_ion'
      !pause 'OUTPUT'
      deallocate(Periodic_table)    ! done with the periodic table, clean up
   case (4)
      bunch_name = 'hole'
   case (5)
      bunch_name = 'muon'
   endselect
   File_name = trim(adjustl(File_name))//'_'//trim(adjustl(bunch_name))
   
   ! Output will also contain the energy of incomming radiation particle:
   !write(NRG, '(f16.2)') bunch(1)%E
   call order_of_energy(bunch(1)%E, temp, E_reduced=E_reduced)   ! module "Little_subroutines"
   write(NRG, '(f16.2)') E_reduced
   NRG = trim(adjustl(NRG))//'_'//trim(adjustl(temp))

   File_name = trim(adjustl(File_name))//'_'//trim(adjustl(NRG))
   
   ! Output will also contain the number of bunches:
   Nsiz = size(bunch(:))
   if (Nsiz > 1) then
      write(temp, '(i8)') Nsiz
      File_name = trim(adjustl(File_name))//'_'//trim(adjustl(temp))//'_bunches'
   endif
   
   ! Check if directory with the same name already exists. If yes, add a number at the end:
   File_name2 = File_name
   i = 0
   inquire(DIRECTORY=trim(adjustl(File_name2)),exist=file_exist)    ! check if input file excists
   do while (file_exist)
      i = i + 1
      write(ch1,'(i6)') i
      write(File_name2,'(a,a,a)') trim(adjustl(File_name)), '_v', trim(adjustl(ch1))
      inquire(DIRECTORY=trim(adjustl(File_name2)),exist=file_exist)    ! check if input file excists
   enddo
   ! Save the output direcory name:
   numpar%output_path = File_name2
   
   ! Create a new directory for output files:
   command='mkdir '//trim(adjustl(File_name2))	! to create a folder use this command
#ifdef _OPENMP
   iret = system(trim(adjustl(command)))    ! create the folder
#else
   CALL system(command)	! create the folder
#endif
   
   ! Copy the input parameters files into this folder, to store all the details for reproducibility:
   if (numpar%new_input_format) then ! new format used
      chtest = trim(adjustl(numpar%input_path))//trim(adjustl(m_input_minimal))
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         call copy_file(trim(adjustl(chtest)),trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_files"
      else ! it is linux
         call copy_file(trim(adjustl(chtest)),trim(adjustl(numpar%output_path))) ! module "Dealing_with_files"
      endif
   else     ! old format used (2 files)
      chtest = trim(adjustl(numpar%input_path))//trim(adjustl(m_input_data))
      chtest1 = trim(adjustl(numpar%input_path))//trim(adjustl(m_numerical_parameters))
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         call copy_file(trim(adjustl(chtest)),trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_files"
         call copy_file(trim(adjustl(chtest1)),trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_files"
      else ! it is linux
         call copy_file(trim(adjustl(chtest)),trim(adjustl(numpar%output_path))) ! module "Dealing_with_files"
         call copy_file(trim(adjustl(chtest1)),trim(adjustl(numpar%output_path))) ! module "Dealing_with_files"
      endif
   endif

   ! Save name for the parameters file:
   numpar%FILE_parameters = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_output_parameters))

   ! Open communication channel for the user:
   ! Set the name of the file thru which user can communicate with the program:
   numpar%FILE_communication	= trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_communication))
   ! file thru which user can communicate with the program:
   numpar%FN_communication = 110    ! can't use newunit in this case, set number manually
   open(unit=numpar%FN_communication, FILE = trim(adjustl(numpar%FILE_communication)), status = 'replace')
   ! get the time when the file was last modified
   close(numpar%FN_communication)   ! close file
   call get_file_stat(trim(adjustl(File_name)), Last_modification_time=numpar%MOD_TIME)	! module "Dealing_with_files"
end subroutine make_output_folder

 
!===================================================
! Printing out the title of the program:
subroutine print_TREKIS4_label(print_to)
   integer, intent(in) :: print_to
   !------------------
   write(print_to, '(a)') trim(adjustl(m_starline))
   write(print_to, '(a)') '        _______   ____    _____   _   _   _    ___           '
   write(print_to, '(a)') '       |__   __| |  _ \  |  ___| | | / / | |  / __|     __   '
   write(print_to, '(a)') '          | |    | |_) ) | |___  | |/ /  | | ( (_      /  |  '
   write(print_to, '(a)') '          | |    |    /  |  ___| |   (   | |  \_ \    /   |  '
   write(print_to, '(a)') '          | |    | |\ \  | |___  | |\ \  | |  __) )  / /| |  '
   write(print_to, '(a)') '          |_|    |_| \_\ |_____| |_| \_\ |_| |___/  / /_| |  '
   write(print_to, '(a)') '                                                   |____   | '
   write(print_to, '(a)') '                                                        |_|  '
   write(print_to, '(a)') trim(adjustl(m_starline))
   write(print_to, '(a)') '       TREKIS: Time-REsolved Kinetics in Irradiated Solids'
   write(print_to, '(a)') '       '//trim(adjustl(m_TREKIS_version))
   write(print_to, '(a)') '       Available at https://github.com/N-Medvedev/TREKIS-4'

   write(print_to, '(a)') trim(adjustl(m_starline))
end subroutine print_TREKIS4_label



subroutine Print_title(print_to, used_target, numpar, bunch, MD_atoms, MD_supce, MD_pots, do_lable)
   integer, intent(in) :: print_to ! the screen, or file
   type(Matter), intent(in) :: used_target	! parameters of the target
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! radiation bunch parameters
   type(Atom), dimension(:), intent(inout), allocatable :: MD_atoms ! all atoms in MD as objects
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   type(MD_potential), dimension(:,:), intent(inout), allocatable :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   logical, intent(in) :: do_lable  ! printout lable or not
   !------------------------
   integer i , Nsiz, j, k
   character(100) :: text, text1, text2, text3
   logical :: MC_output, MD_output
   !------------------------

   if (do_lable) then   ! print TRKEIS lable
      call print_TREKIS4_label(print_to)  ! below
   else
      write(print_to,'(a)') trim(adjustl(m_starline))
   endif
   !write(print_to,'(a)') trim(adjustl(m_starline))
   !write(print_to,'(a)') '*      TREKIS: Time-Resolved Kinetics in Irradiated Solids       *'
   !write(print_to,'(a)') trim(adjustl(m_starline))
   !**********************************************
   if (numpar%new_input_format) then ! new format used
      write(print_to,'(a)') ' Input read from the file: '//trim(adjustl(m_input_minimal))
   else
      write(print_to,'(a)') ' Input read from files: '//trim(adjustl(m_input_data))//' and '//&
                                                       trim(adjustl(m_numerical_parameters))
   endif

   ! TARGET:
   write(print_to,'(a,a,a)') ' Performed for the following parameters:'
   write(print_to,'(a,a)') ' Material : ', trim(adjustl(used_target%Name))
   write(text,'(i10)') used_target%NOC
   if (used_target%NOC == 1) then
      text1 = 'target'
   else
      text1 = 'targets'
   endif
   write(print_to,'(a,a,a)') ' Contains ', trim(adjustl(text)), ' '//trim(adjustl(text1))
   do i = 1, used_target%NOC
      write(text,'(i10)') i
      write(print_to,'(a,a,a)') ' Target #', trim(adjustl(text)), ' is '//trim(adjustl(used_target%Material(i)%Name))
      write(print_to,'(a)') ' Chemical formula interpreted as: ' 
      do j = 1, used_target%Material(i)%N_Elements
         write(text,'(f12.6)') used_target%Material(i)%Elements(j)%percentage
         write(text1,'(i3)') used_target%Material(i)%Elements(j)%Zat
         write(print_to,'(a,a,a,a,a)') ' '//trim(adjustl(text)), ' of ', trim(adjustl(used_target%Material(i)%Elements(j)%Name)), ' (element #', trim(adjustl(text1))//')'
      enddo
      write(text,'(e12.6)') used_target%Material(i)%At_Dens
      write(print_to, '(a,a,a)') ' With the atomic density: ', trim(adjustl(text)), ' [1/cm^3]'
   enddo

   !**********************************************
   ! RADIATION:
   write(print_to,'(a)') trim(adjustl(m_starline))
   Nsiz = size(bunch)
   write(text,'(i10)') Nsiz
   select case (Nsiz)
   case (1)
      text1 = trim(adjustl(text))//' bunch'
      write(print_to,'(a,a)') ' The target is irradiated with ', trim(adjustl(text1))
   case (2:)   
      text1 = trim(adjustl(text))//' bunches'
      write(print_to,'(a,a)') ' The target is irradiated with ', trim(adjustl(text1))
   case default
      write(print_to,'(a)') 'No irradiation is modelled'
   end select
   
   do i = 1, Nsiz ! for each bunch:
      select case(bunch(i)%KOP)	! use the first radiation bunch to define the name
      case (0)  ! Kind of particle: 0=photon, 1=electron, 2=positron, 3=SHI, 4=hole, 5=muon
         text = 'photon'
      case (1)
         text = 'electron'
      case (2)
         text = 'positron'
      case (3)
         text = 'ion '//trim(adjustl(bunch(i)%Name))
      case (4)
         text = 'VB hole'
      case (5)
         text = 'muon'
      endselect
      write(text1,'(i10)') i
      write(print_to,'(a,a,a,a)') ' Bunch #', trim(adjustl(text1)), ' consists of ', trim(adjustl(text))
      !write(text1,'(f16.2)') bunch(i)%E
      !write(print_to,'(a,a,a)') ' with energy ', trim(adjustl(text1)), ' [eV]'
      call print_energy(bunch(i)%E, 'with energy', print_to) ! module "Little_subroutines"
   enddo
   
  !**********************************************
  ! NUMERICS:
   write(print_to,'(a)') trim(adjustl(m_starline))
       
   if (numpar%DO_MC) then
      write(print_to,'(a,a)') ' MC module is active  (v)'
   else 
      write(print_to,'(a,a)') ' MC module is off  (x)'
   endif
   
   if (numpar%DO_MD) then
      write(print_to,'(a,a)') ' MD module is active  (v)'
   else 
      write(print_to,'(a,a)') ' MD module is off  (x)'
   endif
   
   if (numpar%DO_TTM) then
      write(print_to,'(a,a)') ' TTM module is active (v)'
   else 
      write(print_to,'(a,a)') ' TTM module is off (x)'
   endif

#ifdef _OPENMP
   write(text1,'(i10)') numpar%NOMP
   write(print_to,'(a,a)') ' Number of threads for OPENMP: ', trim(adjustl(text1))
#else ! if you set to use OpenMP in compiling: 'make OMP=no'
   write(print_to,'(a)') ' The code is compiled without OPENMP'
#endif

   !**********************************************
   write(print_to,'(a)') trim(adjustl(m_starline))
   
   write(text1,'(i10)') numpar%NMC
   write(print_to,'(a,a)') ' Number of MC iterations: ', trim(adjustl(text1))
   
   if (numpar%recalculate_MFPs) then  ! Force program to recalculate cross sections and MFPs even if there are files containing precalculated ones
      write(print_to,'(a)') ' Cross sections are recalculated in this simulation run'
   else
      write(print_to,'(a)') ' Cross sections are read from files (if exist)'
   endif
   !**********************************************

   write(print_to,'(a)') ' Databases used: '//trim(adjustl(m_EADL))//' and '//trim(adjustl(m_EPDL))

   !**********************************************
   ! MODELS:
   write(print_to,'(a)') trim(adjustl(m_starline)) 
   write(print_to,'(a)') ' Models used:'

   ! CDF numerical parameters:
   !if (numpar%verbose) then
      select case (abs(numpar%CDF_CS_method))
      case (1)
         write(print_to,'(a)') ' Using TREKIS-3 default integration grid'
      case (2)
         write(print_to,'(a)') ' Using TREKIS-4 default integration grid'
      case default
         write(print_to,'(a)') ' Using accelerated TREKIS-4 integration grid'
      end select

      if (numpar%CDF_CS_method < 0) then
         write(print_to,'(a)') ' Integrating diff.CS every time required'
      else
         write(print_to,'(a)') ' Precalculated diff.CS are interpolated'
      endif

      select case (numpar%CDF_plasmon)
      case default
         write(print_to,'(a)') ' Plasmon integration limit unused'
      case (1)
         write(print_to,'(a)') ' Plasmon integration limit used'
      end select


      write(text1,'(f16.2)') numpar%CDF_Eeq_factor
      write(print_to,'(a)') ' Coefficient switching from Ritchie to Delta CDF: E = k*Wmin (inelastic): '//trim(adjustl(text1))

      write(text1,'(f16.2)') numpar%CDF_Eeq_elast
      write(print_to,'(a)') ' Coefficient switching from Ritchie to Delta CDF: E = k*Wmin (elastic): '//trim(adjustl(text1))

      select case (numpar%CDF_elast_Zeff)
      case default
         write(print_to,'(a)') ' Effective charge of target atoms: Barkas-like'
      case (1)
         write(print_to,'(a)') ' Effective charge of target atoms: Z=1'
      end select

      write(text1,'(i0)') numpar%CDF_int_n_inelastE
      write(print_to,'(a)') ' Effective number of grid points for inelastic CS (energy): '//trim(adjustl(text1))

      write(text1,'(i0)') numpar%CDF_int_n_inelastQ
      write(print_to,'(a)') ' Effective number of grid points for inelastic CS (momentum): '//trim(adjustl(text1))

      write(text1,'(i0)') numpar%CDF_int_n_elastE
      write(print_to,'(a)') ' Effective number of grid points for elastic CS (energy): '//trim(adjustl(text1))

      write(text1,'(i0)') numpar%CDF_int_n_elastQ
      write(print_to,'(a)') ' Effective number of grid points for elastic CS (momentum): '//trim(adjustl(text1))

      write(print_to,'(a)') trim(adjustl(m_dashline))
   !endif

   ! Photons:
   select case (numpar%Ph_absorb) ! Include photoabsorption
   case (1)
      write(print_to,'(a)') ' Photon absorption: CDF'
   case (2)
      write(print_to,'(a)') ' Photon absorption: EPDL database'
   case default	! no Compton effect included
      write(print_to,'(a)') ' Photon absorption: No'
   end select
   
   select case (numpar%Ph_Compton)
   case (1)	! Include Compton:
      write(print_to,'(a)') ' Photon Compton effect: PENELOPE'
   case default	! no Compton effect included
      write(print_to,'(a)') ' Photon Compton effect: No'
   end select
   
   select case (numpar%Ph_Pairs)
   case (1)	! Include Compton:
      write(print_to,'(a)') ' Photon pair creation: BH'
   case default	! no Compton effect included
      write(print_to,'(a)') ' Photon pair creation: No'
   end select
   
   select case (numpar%Ph_Thomson)
   case (1)	! Include Rayleigh:
      write(print_to,'(a)') ' Photon Rayleigh: PENELOPE'
   case default	! no Compton effect included
      write(print_to,'(a)') ' Photon Rayleigh: No'
   end select

   if (numpar%Ph_Cutoff > 0.0d0) then
      if (abs(numpar%Ph_Cutoff) < 1.0d7) then
         write(text,'(f12.3)') numpar%Ph_Cutoff
      else
         write(text,'(es20.3)') numpar%Ph_Cutoff
      endif

      write(print_to,'(a,a)') ' Photon energy cut-off: ', trim(adjustl(text))//' [eV]'
   else
      write(print_to,'(a)') ' No energy cut-off for photon transport'
   endif

   if (numpar%Ph_att_eff > 0.0d0) then
      write(text,'(f12.3)') numpar%Ph_att_eff
      write(print_to,'(a,a)') ' Photon attenuation length is set to: ', trim(adjustl(text))//' [A]'
   else
      write(print_to,'(a)') ' Photon attenuation length is taken from database'
   endif

   !**********************************************
   ! Electrons
   select case (numpar%El_inelast)
   case (1)	! CDF
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Ritchie CDF'
      case (1)   ! Mermin
         text = 'Mermin CDF'
      endselect
   case (2)	! RBEB
      text = 'RBEB'
   case (3)	! CDF
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Delta-CDF (with Ritchie)'
      case (1)   ! Mermin
         text = 'Delta-CDF (with Mermin)'
      endselect
   case (4)	! Nonrelativistic CDF
      text = 'Nonrel-CDF'
   case (5)	! Singe-pole delta CDF
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'SP-delta-CDF (with Ritchie)'
      case (1)   ! Mermin
         text = 'SP-delta-CDF (with Mermin)'
      endselect
   case default	! exclude
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Electron inelastic scattering: ', trim(adjustl(text))
   
   select case (numpar%El_elastic)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   case (1)	! CDF
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'CDF (with Ritchie)'
      case (1)   ! Mermin
         text = 'CDF (with Mermin)'
      endselect
   case (2)	! Mott
      text = 'Mott'
   case (3)	! DSF
      text = 'DSF'
   case (5)	! Singe-pole delta CDF
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'SP-delta-CDF (with Ritchie)'
      case (1)   ! Mermin
         text = 'SP-delta-CDF (with Mermin)'
      endselect
   case default	! exclude
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Electron elastic scattering: ', trim(adjustl(text))
   
   select case (numpar%El_brems)
   case (1)	! Bethe-Heitler-Wenzel; see Salvat et al., NIMB 1992
      text = 'BHW'
   case default	! exclude
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Electron Bremsstrahlung: ', trim(adjustl(text))

   if (numpar%El_Cutoff > 0.0d0) then
      if (abs(numpar%El_Cutoff) < 1.0d7) then
         write(text,'(f12.3)') numpar%El_Cutoff
      else
         write(text,'(es20.3)') numpar%El_Cutoff
      endif

      write(print_to,'(a,a)') ' Electron energy cut-off: ', trim(adjustl(text))//' [eV]'
   else
      write(print_to,'(a)') ' No energy cut-off for electron transport'
   endif
   
   !**********************************************
   ! Ions:
   select case (numpar%SHI_inelast)
   case (1:3)	! CDF
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Delta-CDF (with Ritchie)'
      case (1)   ! Mermin
         text = 'Delta-CDF (with Mermin)'
      endselect
   case (4)
      text = 'Ritchie CDF'
   case (5) ! SP
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'SP-delta-CDF (with Ritchie)'
      case (1)   ! Mermin
         text = 'SP-delta-CDF (with Mermin)'
      endselect
   case default
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Ion inelastic scattering: ', trim(adjustl(text))

   if (numpar%SHI_Cutoff > 0.0d0) then
      if (abs(numpar%SHI_Cutoff) < 1.0d7) then
         write(text,'(f12.3)') numpar%SHI_Cutoff
      else
         write(text,'(es20.3)') numpar%SHI_Cutoff
      endif

      write(print_to,'(a,a)') ' SHI energy cut-off: ', trim(adjustl(text))//' [eV]'
   else
      write(print_to,'(a)') ' No energy cut-off for SHI transport'
   endif
  
   !**********************************************
   ! Positrons:
   select case (numpar%Pos_inelast)
   case (1:3)	! CDF
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Delta-CDF (with Ritchie)'
      case (1)   ! Mermin
         text = 'Delta-CDF (with Mermin)'
      endselect
   case (4)
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Ritchie CDF'
      case (1)   ! Mermin
         text = 'Mermin CDF'
      endselect
   case (5) ! SP
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'SP-delta-CDF (with Ritchie)'
      case (1)   ! Mermin
         text = 'SP-delta-CDF (with Mermin)'
      endselect
   case default
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Positron inelastic scattering: ', trim(adjustl(text))

   select case (numpar%Pos_elastic)
   case (1)
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Ritchie CDF'
      case (1)   ! Mermin
         text = 'Mermin CDF'
      endselect
   case (2)
      text = 'Mott'
   case (3)
      text = 'DSF'
   case default
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Positron elastic scattering: ', trim(adjustl(text))
   
   select case (numpar%Pos_Brems)
   case (1)
      text = 'BHW'
   case default
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Positron Bremsstrahlung: ', trim(adjustl(text))
   
   select case (numpar%Pos_annih)
   case (1)
      text = 'Heitler'
   case default
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Positron annihilation: ', trim(adjustl(text))

   if (numpar%Pos_Cutoff > 0.0d0) then
      if (abs(numpar%Pos_Cutoff) < 1.0d7) then
         write(text,'(f12.3)') numpar%Pos_Cutoff
      else
         write(text,'(es20.3)') numpar%Pos_Cutoff
      endif

      write(print_to,'(a,a)') ' Positron energy cut-off: ', trim(adjustl(text))//' [eV]'
   else
      write(print_to,'(a)') ' No energy cut-off for positron transport'
   endif
   
   !**********************************************
   ! Holes
   select case (numpar%H_inelast)
   case (1:5)	! CDF with Ritchi's oscillators
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Ritchie CDF'
      case (1)   ! Mermin
         text = 'Mermin CDF'
      endselect
   case default	! exclude
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Valence hole inelastic scattering: ', trim(adjustl(text))
   
   select case (numpar%H_elast)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   case (1)	! CDF - only total cross section
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Ritchie CDF'
      case (1)   ! Mermin
         text = 'Mermin CDF'
      endselect
   case (2)	! Mott - cross sections for each type of atom
      text = 'Mott'
   case default	! exclude
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Valence hole elastic scattering: ', trim(adjustl(text))
   
   select case (numpar%H_Auger) !  Auger decays:  0=excluded, 1=EADL
   case (1)
      text = 'EADL'
   case default	! exclude
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Core hole Auger decay: ', trim(adjustl(text))
   
   select case (numpar%H_Radiat) !  Radiative decays: 0=excluded, 1=EADL
   case (1)
      text = 'EADL'
   case default	! exclude
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Core hole radiative decay: ', trim(adjustl(text))

   if (numpar%H_Cutoff > 0.0d0) then
      if (abs(numpar%H_Cutoff) < 1.0d7) then
         write(text,'(f12.3)') numpar%H_Cutoff
      else
         write(text,'(es20.3)') numpar%H_Cutoff
      endif
      write(print_to,'(a,a)') ' VB holes energy cut-off: ', trim(adjustl(text))//' [eV]'
   else
      write(print_to,'(a)') ' No energy cut-off for VB holes transport'
   endif
   

   !**********************************************
   ! Muons
   select case (numpar%Mu_inelast)
   case (1:3)	! CDF
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Delta-CDF (with Ritchie)'
      case (1)   ! Mermin
         text = 'Delta-CDF (with Mermin)'
      endselect
   case (4)
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Ritchie CDF'
      case (1)   ! Mermin
         text = 'Mermin CDF'
      endselect
   case (5) ! SP
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'SP-delta-CDF (with Ritchie)'
      case (1)   ! Mermin
         text = 'SP-delta-CDF (with Mermin)'
      endselect
   case default
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Muon inelastic scattering: ', trim(adjustl(text))

   select case (numpar%Mu_elastic)	! elastic scattering: 0=excluded, 1=CDF
   case (1)
      select case (numpar%CDF_model)
      case default   ! Ritchie
         text = 'Ritchie CDF'
      case (1)   ! Mermin
         text = 'Mermin CDF'
      endselect
   case (2)
      text = 'Mott'
   case (3)
      text = 'DSF'
   case default
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Muon elastic scattering: ', trim(adjustl(text))

   select case (numpar%Mu_Brems) ! Muons bremsstrahlung: 0=excluded, 1=BHW
   case (1)
      text = 'BHW'
   case default	! exclude
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Muon bremsstrahlung: ', trim(adjustl(text))

   select case (numpar%Mu_Pairs) !  Electron-positron pair creation: 0=excluded, 1=... (not ready)
   case (1)
      text = '(NOT READY)'
   case default	! exclude
      text = 'No'
   end select
   write(print_to,'(a,a)') ' Muon pair creation: ', trim(adjustl(text))

   if (numpar%Mu_Cutoff > 0.0d0) then
      if (abs(numpar%Mu_Cutoff) < 1.0d7) then
         write(text,'(f12.3)') numpar%Mu_Cutoff
      else
         write(text,'(es20.3)') numpar%Mu_Cutoff
      endif

      write(print_to,'(a,a)') ' Muon energy cut-off: ', trim(adjustl(text))//' [eV]'
   else
      write(print_to,'(a)') ' No energy cut-off for muon transport'
   endif


   !**********************************************
   ! MD model parameters:
   if (numpar%DO_MD) then
      write(print_to,'(a)') trim(adjustl(m_starline))
      write(text1,'(i16)') size(MD_atoms)
      write(print_to,'(a)') ' Number of atoms in MD is: '//trim(adjustl(text1))
      
      write(text1,'(f12.3)') MD_supce%x_start
      write(text2,'(f12.3)') MD_supce%x_end
      write(print_to,'(a)') ' Supercell size along X: '//trim(adjustl(text1))//' to '//trim(adjustl(text2))
      write(text1,'(f12.3)') MD_supce%y_start
      write(text2,'(f12.3)') MD_supce%y_end
      write(print_to,'(a)') ' Supercell size along Y: '//trim(adjustl(text1))//' to '//trim(adjustl(text2))
      write(text1,'(f12.3)') MD_supce%z_start
      write(text2,'(f12.3)') MD_supce%z_end
      write(print_to,'(a)') ' Supercell size along Z: '//trim(adjustl(text1))//' to '//trim(adjustl(text2))
      
      write(print_to,'(a)') ' The following MD potentials are used:'
      do i = 1, size(MD_pots,1)
         do j = i, size(MD_pots,2)
            write(print_to,'(a)') ' Interaction '//trim(adjustl(MD_pots(i,j)%El1))//'-'//   &
                                    trim(adjustl(MD_pots(i,j)%El2))//' is modeled with:'
            do k = 1, size(MD_pots(i,j)%Set)
               write(text,'(i3)') k
               write(print_to,'(a)') ' '//trim(adjustl(text))//') '//trim(adjustl(MD_pots(i,j)%Set(k)%Par%Name))
            enddo ! k
         enddo ! j
      enddo ! i

      if (numpar%use_thermostat) then    ! only if user included it
         select case (numpar%thermostat_inx)   ! which thermostat to use
         case (0)  ! Berendsen
            write(print_to,'(a)') ' Berendsen thermostat is used in MD'
            write(text1,'(f12.3)') MD_supce%Bath_T
            write(print_to,'(a)') ' with bath temperature: '//trim(adjustl(text1))//' [K]'
            write(text1,'(f12.3)') MD_supce%therm_tau
            write(print_to,'(a)') ' and characteristic time: '//trim(adjustl(text1))//' [fs]'
         case default
            write(print_to,'(a)') ' No thermostat is used in MD'
         end select
      else
         write(print_to,'(a)') ' No thermostat is used in MD'
      endif

      if (numpar%damp_pressure) then    ! only if user included it
         write(print_to,'(a)') ' Pressure dampening viscosity is used in MD'
         write(text1,'(f12.3)') MD_supce%press_tau
         write(print_to,'(a)') ' with characteristic time: '//trim(adjustl(text1))//' [fs]'
      else
         write(print_to,'(a)') ' No pressure dampening viscosity is used in MD'
      endif

      if (numpar%MD_force_ind == 1) then
         write(print_to,'(a)') ' Interatomic forces calcualted as numerical derivative of the potential'
      else
         write(print_to,'(a)') ' Interatomic forces calcualted analytically'
      endif

   endif

   !**********************************************
   MC_output = .false. ! to start with
   MD_output = .false. ! to start with
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') 'Optional output: '
   !---------------------
   write(print_to,'(a)') ' In MC module: '

   if (numpar%vel_theta_grid_par%along_axis) then
      write(print_to,'(a)') '  Angular velosity distribution: theta=acos(Vz/V)'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%NRG_grid_par%along_axis) then
      write(print_to,'(a)') '  Energy spectrum of particles'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%Spectr_grid_par(1)%along_axis) then
      write(print_to,'(a)') '  Spectrum resolved along X'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%Spectr_grid_par(2)%along_axis) then
      write(print_to,'(a)') '  Spectrum resolved along Y'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%Spectr_grid_par(3)%along_axis) then
      write(print_to,'(a)') '  Spectrum resolved along Z'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%Spectr_grid_par(8)%along_axis) then
      write(print_to,'(a)') '  Spectrum resolved along R (cylindric coordinates)'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%grid_par(1)%along_axis) then
      write(print_to,'(a)') '  Distribution along X'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%grid_par(2)%along_axis) then
      write(print_to,'(a)') '  Distribution along Y'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%grid_par(3)%along_axis) then
      write(print_to,'(a)') '  Distribution along Z'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%grid_par(8)%along_axis) then
      write(print_to,'(a)') '  Distribution along R (cylindric coordinates)'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%grid_par(11)%along_axis) then
      write(print_to,'(a)') '  2d-distribution along R and Z (cylindric coordinates)'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%Theta_grid_par(1)%along_axis) then
      write(print_to,'(a)') '  Theta-distribution resolved along X'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%Theta_grid_par(2)%along_axis) then
      write(print_to,'(a)') '  Theta-distribution resolved along Y'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (numpar%Theta_grid_par(3)%along_axis) then
      write(print_to,'(a)') '  Theta-distribution resolved along Z'
      MC_output = .true.   ! mark that there was any optional output
   endif

   if (.not.MC_output) then
      write(print_to,'(a)') '  none'
   endif

   !---------------------
   if (numpar%DO_MD) then
      write(print_to,'(a)') ' In MD module: '

      if (numpar%print_MC_MD_energy) then
         write(print_to,'(a)') '  MC-MD energy transfer'
         MD_output = .true.   ! mark that there was any optional output
      endif

      if (numpar%n_MSD > 0) then
         write(print_to,'(a)') '  Mean atomic displacements'
         MD_output = .true.   ! mark that there was any optional output
      endif

      if (numpar%print_MD_R_xyz) then
         write(print_to,'(a)') '  Atomic coordinates in XYZ format'
         MD_output = .true.   ! mark that there was any optional output
      endif

      if (numpar%print_MD_V_xyz) then
         write(print_to,'(a)') '  Atomic velosities in XYZ format'
         MD_output = .true.   ! mark that there was any optional output
      endif

      if (numpar%print_MD_LAMMPS) then
         write(print_to,'(a)') '  Input file for LAMMPS at the final time instant'
         MD_output = .true.   ! mark that there was any optional output
      endif

      if (.not.MD_output) then
         write(print_to,'(a)') '  none'
      endif
   endif
   
   !**********************************************
9999   write(print_to,'(a)') trim(adjustl(m_starline))
end subroutine Print_title


subroutine Print_atomic_parameters(print_to, used_target)
   integer, intent(in) :: print_to ! the screen, or file
   type(Matter), intent(in) :: used_target	! parameters of the target
   integer :: i, j, k
   character(100) :: text
!    write(print_to,'(a)') trim(adjustl(m_starline))
   do i = 1, used_target%NOC	! for all target constituents
      do j = 1, used_target%Material(i)%N_Elements	! for all elements this target material contains
         write(text, '(i2)') used_target%Material(i)%Elements(j)%N_shl
         if (used_target%Material(i)%Elements(j)%N_shl > 1) then
            write(print_to,'(a,a,a,a,a)') 'Element ', trim(adjustl(used_target%Material(i)%Elements(j)%Name)), ' has ', trim(adjustl(text)), ' shells:'
         else
            write(print_to,'(a,a,a,a,a)') 'Element ', trim(adjustl(used_target%Material(i)%Elements(j)%Name)), ' has ', trim(adjustl(text)), ' shell:'
         endif
         write(print_to, '(a)') ' #	Designator	Name	Valent  Ne	Ip[eV]	Ek[eV]	Auger[fs]	Radiative[fs]'
         do k = 1, used_target%Material(i)%Elements(j)%N_shl	! for all shell of this element from this target
            write(print_to, '(i3, i3, a6, L, f5.2, f12.4, f12.4, es12.4, es12.4)') k, used_target%Material(i)%Elements(j)%Shl_dsgnr(k), &
             ' '//trim(adjustl(used_target%Material(i)%Elements(j)%Shell_name(k)))//' ', &
             used_target%Material(i)%Elements(j)%valent(k), &
             used_target%Material(i)%Elements(j)%Ne_shell(k), &
             used_target%Material(i)%Elements(j)%Ip(k), used_target%Material(i)%Elements(j)%Ek(k), &
             used_target%Material(i)%Elements(j)%Auger(k), used_target%Material(i)%Elements(j)%Radiat(k)
! used_target%Material(i)%Elements(j)%f_rad, used_target%Material(i)%Elements(j)%f_auger
         enddo ! k
         write(print_to,'(a)') trim(adjustl(m_starline))
      enddo ! j

      write(text, '(f12.4)') used_target%Material(i)%DOS%Egap
      write(print_to,'(a)') ' Band gap in this material: '//trim(adjustl(text))//' [eV]'
      write(print_to,'(a)') trim(adjustl(m_starline))
  enddo ! i
end subroutine Print_atomic_parameters


subroutine Print_CDF_parameters(print_to, used_target)
   integer, intent(in) :: print_to ! the screen, or file
   type(Matter), intent(in), target :: used_target	! parameters of the target
   integer :: i, j, m, k
   real(8) :: elem_contrib, N_elem 
   type(Atom_kind), pointer :: Element
      
   TRGT:do i = 1, used_target%NOC		! number of different components of the target (layers, clusters, etc.)   
      N_elem = dble(SUM(used_target%Material(i)%Elements(:)%percentage))	! number of elements in this compound material
      
      write(print_to, '(a,a)') 'CDF coefficients (normalized per molecule) in ', trim(adjustl(used_target%Material(i)%Name))
      ELMNT:do j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent
         Element => used_target%Material(i)%Elements(j)	! all information about this element
         
         elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         
         write(print_to, '(a,a)') 'Element: ', trim(adjustl(Element%Name))
         ! Get CDFs for each shell:
         SHL:do m = 1, Element%N_shl		! for all shells in this elements
            write(print_to,'(a,a,a)') 'Shell: ', trim(adjustl(Element%Shell_name(m))), ', A E0 Gamma    alpha'
            do k = 1, size(Element%CDF(m)%A)	! all coefficients
               write(print_to,'(f12.4, f12.4, f12.4, f12.4)') Element%CDF(m)%A(k) * elem_contrib, Element%CDF(m)%E0(k), &
                                                Element%CDF(m)%Gamma(k), Element%CDF(m)%alpha(k)*elem_contrib
               write(print_to, '(i3, i3, a, f5.2, f12.4, f12.4, es12.4, es12.4)')
            enddo
         enddo SHL
         write(print_to,'(a)') trim(adjustl(m_starline))
         
      enddo ELMNT
      
      ! Valence:
      if (allocated(used_target%Material(i)%CDF_valence%A)) then
         write(print_to,'(a)') 'Valence band, A E0 Gamma    alpha'
         do k = 1, size(used_target%Material(i)%CDF_valence%A)	! all coefficients
            write(print_to,'(f12.4, f12.4, f12.4, f12.4)') used_target%Material(i)%CDF_valence%A(k),  used_target%Material(i)%CDF_valence%E0(k), &
                                            used_target%Material(i)%CDF_valence%Gamma(k), used_target%Material(i)%CDF_valence%alpha(k)
         enddo
      endif
      write(print_to,'(a)') trim(adjustl(m_starline))
      
      ! Phonons:
      if (allocated(used_target%Material(i)%CDF_phonon%A)) then
         write(print_to,'(a)') 'Phonons, A E0 Gamma alpha'
         do k = 1, size(used_target%Material(i)%CDF_phonon%A)	! all coefficients
            write(print_to,'(f16.6, f16.6, f16.6, f16.6)') used_target%Material(i)%CDF_phonon%A(k),  used_target%Material(i)%CDF_phonon%E0(k), &
                                            used_target%Material(i)%CDF_phonon%Gamma(k), used_target%Material(i)%CDF_phonon%alpha(k)
         enddo
      endif
      write(print_to,'(a)') trim(adjustl(m_starline))
      
   enddo TRGT
   
   nullify(Element)
end subroutine Print_CDF_parameters



subroutine print_duration(print_to, chtext, ctim)
   integer, intent(in) :: print_to	! into which unit to print (file number or screen)
   character(*), intent(in) :: chtext	! time duration to print out
   integer, intent(in) :: ctim(8)	! time stamps
   
   write(print_to,'(a,a)') ' Duration of execution of program: ', trim(adjustl(chtext))
   call print_time(' Started  at', ind=0, ctim=ctim, print_to=print_to)	! module "Little_subroutines"
   call print_time(' Finished at', ind=0, print_to=print_to)	! module "Little_subroutines"
   write(print_to,'(a)') trim(adjustl(m_starline))
end subroutine print_duration


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! Communication with the used via file:
subroutine communicate_with_user(FN, time, numpar)
   integer, intent(in) :: FN ! file number to read from
   real(8), intent(in) :: time ! current time [fs]
   type(Num_par), intent(inout), target :: numpar	! all numerical parameters
   !----------------------------------------------------------------------
   integer :: Reason, i, MOD_TIM, errnum, sz
   character(200) :: readline, given_line, File_name
   real(8) :: given_num
   logical :: read_well, read_well_2, file_opened
   
   readline = ''  ! to start with
   File_name = trim(adjustl(numpar%FILE_communication))

   inquire(UNIT=FN,opened=file_opened)
   if (file_opened) close(FN) ! for windows, we have to close the file to let the user write into it
   ! Check if the file was modified since the last time:
   call get_file_stat(trim(adjustl(File_name)), Last_modification_time=MOD_TIM) ! module 'Dealing_with_files'
   if (MOD_TIM /= numpar%MOD_TIME) then ! open file again only if it was modified by the user
      numpar%MOD_TIME = MOD_TIM ! save new time of the last modification
      open(UNIT=FN, FILE=trim(adjustl(File_name)),IOSTAT=errnum, SHARED, STATUS='old')
   endif
7777     continue ! in case if the program could not open the file

   inquire(UNIT=FN,opened=file_opened)
   COM_OPEN:if (file_opened) then ! read it
      rewind(FN)  ! to start reading from the start
      ! read the first line
      read(FN,'(a)', IOSTAT=Reason, SIZE=sz, ADVANCE='no') readline
      call read_file(Reason, i, read_well)

      if ( (.not.read_well) .and. (numpar%path_sep == '/') ) then ! if it is Linux
         rewind(FN)  ! to start reading from the start
         read(FN, '(a)', IOSTAT=Reason) readline(1:sz) ! read it again, now knowing the size
         call read_file(Reason, i, read_well)
      endif

      if (read_well) then
         ! Interpret what user wrote:
         call pars_communications(trim(adjustl(readline)), given_line, given_num, read_well_2)   ! below
         ! Act accordingly:
         call act_on_communication(read_well_2, given_line, given_num, numpar, time) ! below
         i = 1
         do while (Reason .EQ. 0) ! read all lines if there is more than one
            i = i + 1
            read(FN,'(a)',IOSTAT=Reason) readline
            call read_file(Reason, i, read_well)    ! module "Dealing_with_files"
            if (Reason .NE. 0) exit
            ! Interpret what user wrote:
            call pars_communications(trim(adjustl(readline)), given_line, given_num, read_well_2)    ! below
            ! Act accordingly:
            call act_on_communication(read_well_2, given_line, given_num, numpar, time)  ! below
         enddo
         rewind(FN)
         write(FN,'(a)') ''
         rewind(FN)
      endif

      call get_file_stat(trim(adjustl(File_name)), Last_modification_time=MOD_TIM) ! module 'Dealing_with_files'
      if (MOD_TIM /= numpar%MOD_TIME) then ! if it was modified by the user, then
         numpar%MOD_TIME = MOD_TIM         ! save new time of the last modification
      endif
      close(FN) ! for windows, we have to close the file to let the user write into it

   endif COM_OPEN
end subroutine communicate_with_user




subroutine pars_communications(readline, out_line, out_num, read_well)
   character(*), intent(in) :: readline
   character(*), intent(out) :: out_line
   real(8), intent(out) :: out_num
   logical, intent(out) :: read_well
   !---------------------------------
   integer :: Reason, i
   read_well = .false.
   out_line = ''
   out_num = 0.0d0

   i = 1    ! to start with
   read(readline,*, IOSTAT=Reason) out_line, out_num
   call read_file(Reason, i, read_well)  ! module "Dealing_with_files"
   if (Reason .LT. 0) then
      print*, 'No descriptor or value found in the communication file'
   else if (Reason .GT. 0) then
      print*, 'Given number interpreted as', out_num, ', it does not match the variable type'
   endif
   if (.not.read_well) then
      print*, 'Comunication format must be as follows:'
      print*, 'Two columns: 1) descriptor; 2) value'
      print*, 'Allowed descriptors: Time; dt; Save_dt; OMP'
   endif
end subroutine pars_communications




subroutine act_on_communication(read_well, given_line, given_num, numpar, time)
   logical, intent(in) :: read_well ! did we read something meaningful from the communication file?
   character(*), intent(in) :: given_line ! line read from the file
   real(8), intent(in) :: given_num  ! number read from the file
   type(Num_par), intent(inout), target :: numpar	! all numerical parameters   
   real(8), intent(in) :: time ! current time [fs]
   !-----------------------------------------------------------
   integer FN, noth
   logical file_opened
   character(200) :: File_name, temp1, temp2
   character(1) path_sep
   path_sep = trim(adjustl(numpar%path_sep))
   if (read_well) then
      FN = numpar%FN_parameters   ! parameters output file, to save info about the user-defined change of parameters
      
      select case(trim(adjustl(given_line)))
      case ('time', 'TIME', 'Time', 'TIme', 'TIMe', 'tIme')
         numpar%t_total = given_num ! total duration of simulation [fs]
         write(temp1,'(f12.3)') time
         write(temp2,'(f12.3)') given_num
         write(6,'(a,a)') 'Duration of simulation is changed to ', trim(adjustl(temp2))
         write(FN,'(a,a,a,a)') 'At time instant of ', trim(adjustl(temp1)), '[fs], duration of simulation is changed to ', trim(adjustl(temp2))
      
      case ('Step', 'STEP', 'step', 'STep', 'STEp', 'dt', 'Dt', 'DT')
         numpar%dt_MD = given_num ! Time step [fs]
!          numpar%halfdt = numpar%dt/2.0d0      ! dt/2, often used
!          numpar%dtsqare = numpar%dt*numpar%halfdt ! dt*dt/2, often used
         write(temp1,'(f12.3)') time
         write(temp2,'(f12.3)') given_num
         write(6,'(a,a)') 'Time-step of simulation is changed to ',  trim(adjustl(temp2))
         write(FN,'(a,a,a,a)') 'At time instant of ',  trim(adjustl(temp1)), '[fs], time-step of simulation is changed to ',  trim(adjustl(temp2)) 
      
      case ('SAVEdt', 'savedt', 'dtsave', 'dtSAVE', 'Savedt', 'SaveDT', 'SaveDt', 'save_dt', 'Save_dt', 'SAVE_dt', 'SAVE_DT')
         numpar%dt_printout = given_num ! save data into files every 'dt_printout' [fs]
         write(temp1,'(f12.3)') time
         write(temp2,'(f12.3)') given_num
         write(6,'(a,a)') 'Time-step of saving output files is changed to ',  trim(adjustl(temp2))
         write(FN,'(a,a,a,a)') 'At time instant of ',  trim(adjustl(temp1)), '[fs], time-step of saving output files is changed to ',  trim(adjustl(temp2))
      
      case ('OMP', 'omp', 'NOMP', 'nomp', 'Nomp', 'N_OMP', 'n_omp')
         ! Reset the OpenMP parallelization options:
         numpar%NOMP = given_num
         write(temp1,'(f12.3)') time
         write(temp2,'(i10)') INT(given_num)
         
#ifdef _OPENMP
         noth = OMP_GET_MAX_THREADS()   ! to chech if the function worked
         call set_OMP_number( numpar%NOMP, .true., 6, 'Reset number of threads in OpenMP to '//trim(adjustl(temp2)) )    ! below
         if ( noth /= OMP_GET_MAX_THREADS() ) then
            write(FN,'(a,a,a,a)') 'At time instant of ',  trim(adjustl(temp1)), '[fs], number of threads in OpenMP is changed to ',  trim(adjustl(temp2))
         else
            write(FN,'(a,a,a,a)') 'At time instant of ',  trim(adjustl(temp1)), '[fs]: unsuccessful attempt to change number of threads in OpenMP to ',  trim(adjustl(temp2))
            write(6,'(a)') 'Number of threads in OpenMP is unchanged: ',  trim(adjustl(temp2))
         endif
#else
         write(FN,'(a)') ' The code compiled without OpenMP, cannot set parallelization'
         write(6,'(a)') 'The code compiled without OpenMP, cannot set parallelization'
#endif

      case default
         print*, 'Could not interpret what is read from the file: ', trim(adjustl(given_line)), given_num
      end select
      
   endif
end subroutine act_on_communication


subroutine set_OMP_number(NOMP, prnt, FN, lin)
   integer, intent(inout) :: NOMP  ! number of threads to be set; negative means = maximal threads available
   logical, intent(in) :: prnt  ! do we want to print out anything?
   integer, intent(in) :: FN    ! file number to print into
   character(*), intent(in), optional :: lin    ! a line to print out
   !------------------------------------
   character(10) :: temp2
   
#ifdef _OPENMP
   call OMP_SET_DYNAMIC(0) ! standard openmp subroutine
   if (NOMP <= 0) then ! use all available processors / threads:
      NOMP = OMP_GET_MAX_THREADS() ! number of threads for openmp defined in INPUT_PARAMETERS.txt
   endif
   call OMP_SET_NUM_THREADS(NOMP) ! number of threads for openmp defined in INPUT_PARAMETERS.txt
   if (prnt) then
      if (present(lin)) then
         write(FN,'(a)') trim(adjustl(lin))
      else  ! printout default message
         write(temp2,'(i10)') INT(NOMP)
         write(FN,'(a,a)') ' The code was compiled with OpenMP parallelization, THREADS: ', trim(adjustl(temp2))
      endif
   endif
#else
   if (prnt) then
      if (present(lin)) then
         write(FN,'(a)') trim(adjustl(lin))
      else  ! printout default message
         write(FN,'(a)') ' The code was compiled without using OpenMP'
      endif
   endif
#endif
end subroutine set_OMP_number





subroutine close_all_output(numpar)
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   !----------------------
   integer :: Nsiz, i
   ! For all output files:

   ! Files with spectra:
   call close_file('close', FN=numpar%FN_spectrum_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_e)	! module "Dealing_with_files"
   if (allocated(numpar%FN_spectrum_h)) then
      Nsiz = size(numpar%FN_spectrum_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_spectrum_h(i))	! module "Dealing_with_files"
      enddo
   endif
   call close_file('close', FN=numpar%FN_spectrum_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_mu)	! module "Dealing_with_files"

   ! Files with theta distributions:
   call close_file('close', FN=numpar%FN_vel_theta_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_vel_theta_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_vel_theta_h)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_vel_theta_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_vel_theta_SHI)	! module "Dealing_with_files"

   ! Files with spectra vs space:
   call close_file('close', FN=numpar%FN_spectrum_ph_X)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_e_X)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_p_X)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_SHI_X)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_h_X)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_spectrum_ph_Y)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_e_Y)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_p_Y)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_SHI_Y)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_h_Y)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_spectrum_ph_Z)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_e_Z)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_p_Z)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_SHI_Z)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_spectrum_h_Z)	! module "Dealing_with_files"

   ! Files with theta distribution vs space:
   call close_file('close', FN=numpar%FN_theta_ph_X)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_e_X)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_p_X)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_SHI_X)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_h_X)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_theta_ph_Y)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_e_Y)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_p_Y)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_SHI_Y)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_h_Y)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_theta_ph_Z)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_e_Z)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_p_Z)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_SHI_Z)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_theta_h_Z)	! module "Dealing_with_files"

   ! File numbers with spatial distributions:
   ! Cartesian:
   call close_file('close', FN=numpar%FN_car_1d_X_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_X_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_X_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_X_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_X_a)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_car_1d_Y_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Y_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Y_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Y_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Y_a)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_car_1d_Z_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Z_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Z_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Z_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Z_a)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_car_1d_X_E_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_X_E_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_X_E_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_X_E_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_X_E_a)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_car_1d_Y_E_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Y_E_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Y_E_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Y_E_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Y_E_a)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_car_1d_Z_E_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Z_E_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Z_E_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Z_E_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_car_1d_Z_E_a)	! module "Dealing_with_files"

   if (allocated(numpar%FN_car_1d_X_h)) then
      Nsiz = size(numpar%FN_car_1d_X_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_car_1d_X_h(i))	! module "Dealing_with_files"
      enddo
   endif
   if (allocated(numpar%FN_car_1d_X_E_h)) then
      Nsiz = size(numpar%FN_car_1d_X_E_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_car_1d_X_E_h(i))	! module "Dealing_with_files"
      enddo
   endif
   if (allocated(numpar%FN_car_1d_Y_h)) then
      Nsiz = size(numpar%FN_car_1d_Y_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_car_1d_Y_h(i))	! module "Dealing_with_files"
      enddo
   endif
   if (allocated(numpar%FN_car_1d_Y_E_h)) then
      Nsiz = size(numpar%FN_car_1d_Y_E_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_car_1d_Y_E_h(i))	! module "Dealing_with_files"
      enddo
   endif
   if (allocated(numpar%FN_car_1d_Z_h)) then
      Nsiz = size(numpar%FN_car_1d_Z_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_car_1d_Z_h(i))	! module "Dealing_with_files"
      enddo
   endif
   if (allocated(numpar%FN_car_1d_Z_E_h)) then
      Nsiz = size(numpar%FN_car_1d_Z_E_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_car_1d_Z_E_h(i))	! module "Dealing_with_files"
      enddo
   endif

   call close_file('close', FN=numpar%FN_cyl_1d_R_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_1d_R_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_1d_R_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_1d_R_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_1d_R_a)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_cyl_1d_R_E_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_1d_R_E_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_1d_R_E_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_1d_R_E_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_1d_R_E_a)	! module "Dealing_with_files"

   if (allocated(numpar%FN_cyl_1d_R_h)) then
      Nsiz = size(numpar%FN_cyl_1d_R_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_cyl_1d_R_h(i))	! module "Dealing_with_files"
      enddo
   endif
   if (allocated(numpar%FN_cyl_1d_R_E_h)) then
      Nsiz = size(numpar%FN_cyl_1d_R_E_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_cyl_1d_R_E_h(i))	! module "Dealing_with_files"
      enddo
   endif

   call close_file('close', FN=numpar%FN_cyl_2d_RL_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_2d_RL_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_2d_RL_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_2d_RL_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_2d_RL_a)	! module "Dealing_with_files"

   call close_file('close', FN=numpar%FN_cyl_2d_RL_E_ph)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_2d_RL_E_e)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_2d_RL_E_p)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_2d_RL_E_SHI)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_cyl_2d_RL_E_a)	! module "Dealing_with_files"

   if (allocated(numpar%FN_cyl_2d_RL_h)) then
      Nsiz = size(numpar%FN_cyl_2d_RL_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_cyl_2d_RL_h(i))	! module "Dealing_with_files"
      enddo
   endif
   if (allocated(numpar%FN_cyl_2d_RL_E_h)) then
      Nsiz = size(numpar%FN_cyl_2d_RL_E_h)
      do i = 1, Nsiz
         call close_file('close', FN=numpar%FN_cyl_2d_RL_E_h(i))	! module "Dealing_with_files"
      enddo
   endif

   ! MD related files:
   call close_file('close', FN=numpar%FN_MD_totals)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_MD_average)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_MD_displacement)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_MD_XYZ)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_MD_V_XYZ)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_MCMD_info)	! module "Dealing_with_files"
   call close_file('close', FN=numpar%FN_MD_LAMMPS)	! module "Dealing_with_files"

end subroutine close_all_output



subroutine pars_communications_old(readline, out_line, out_num, read_well)
   character(*), intent(in) :: readline
   character(*), intent(out) :: out_line
   real(8), intent(out) :: out_num
   logical, intent(out) :: read_well
   character(LEN(readline)) partline, lastpart
   real(8) :: newnum
   integer :: i, leng, firstnum, coun, Reason
   character(*), parameter :: numbers = '0123456789'
   character(*), parameter :: plusminus = '-+'
   logical :: allowpoint, allowminus, allowed
   read_well = .false.
   out_line = ''
   out_num = 0.0d0
   allowpoint = .true.
   allowminus = .true.
   allowed = .true.
   leng = LEN(trim(adjustl(readline))) ! how many characters are in the line

   !print*, readline, trim(adjustl(readline)), leng
   if (leng .GT. 0) then
      partline = ''
      LENGT:do i = 1,leng ! compare all name character by character
         if (verify(trim(adjustl(readline(i:i))),trim(adjustl(numbers))) == 0) then ! it's an integer number
            firstnum = i
            exit LENGT
         endif
      enddo LENGT
      if (firstnum .LT. 2) goto 1111 !skip everything, there is no command to interpret
      partline(1:firstnum-1) = readline(1:firstnum-1) ! part of line with text
      if (firstnum .LT. leng) then ! there is some number, apparently
         lastpart = ''
         coun = 1
         LENGT2:do i = firstnum-2,leng ! compare all name character by character
            if (trim(adjustl(readline(i:i))) == ' ') goto 1112
             if (verify(trim(adjustl(readline(i:i))),trim(adjustl(numbers))) == 0) then ! it's a number
               lastpart(coun:coun) = readline(i:i) ! part of line with numbers
               coun = coun + 1
               if (.not.allowed) then
                  allowpoint = .false. ! only one is allowed
               endif
               allowminus = .false. ! only one is allowed
             else ! it's not a number, but may be a correct symbol:
               select case(trim(adjustl(readline(i:i))))
               case ('.') ! decimal point
                  if (allowpoint) then
                     lastpart(coun:coun) = readline(i:i) ! part of line with numbers
                     coun = coun + 1
                     allowpoint = .false. ! only one is allowed
                     allowminus = .false. ! only one is allowed
                  endif
               case ('-', '+') ! exp format
                  if (i .LT.leng) then ! after the minus/plus must be a number:
                     if ((allowminus) .AND. ((verify(trim(adjustl(readline(i+1:i+1))),trim(adjustl(numbers))) == 0) .OR. (trim(adjustl(readline(i+1:i+1))) .EQ. '.'))) then
                        lastpart(coun:coun) = readline(i:i) ! part of line with numbers
                        coun = coun + 1
                        allowminus = .false. ! only one is allowed
                     endif
                  endif
               case ('e', 'd', 'E', 'D') ! exp format
                  if ((i .LT. leng) .AND. (len(trim(adjustl(lastpart))) .GT. 0) ) then ! after the E must be a number, or a minus:
                     if ((allowed) .AND. (verify(trim(adjustl(readline(i+1:i+1))),trim(adjustl(numbers))) == 0)) then
                        lastpart(coun:coun) = '0' ! set 0 before any 'e' or 'd'
                        lastpart(coun+1:coun+1) = readline(i:i) ! part of line with numbers
                        coun = coun + 2
                        allowed = .false. ! only one is allowed
                     else if (i .LT. leng - 1) then ! may be first minus/plus, and then the number:
                        if ((allowed) .AND. ((verify(trim(adjustl(readline(i+1:i+1))),trim(adjustl(plusminus))) == 0) .AND. (verify(trim(adjustl(readline(i+2:i+2))),trim(adjustl(numbers))) == 0) )) then
                           lastpart(coun:coun) = '0' ! set 0 before any 'e' or 'd'
                           lastpart(coun+1:coun+1) = readline(i:i) ! part of line with numbers
                           coun = coun + 2
                           allowed = .false. ! only one is allowed
                           allowminus = .true. ! one more minus is allowed in the exponent
                        endif
                     endif
                  endif
               end select
             endif
            1112 continue
         enddo LENGT2
      endif

      if (LEN(trim(adjustl(lastpart))) <= 0) then ! no need to even try reading it
         Reason = -1
      else
         read(lastpart,'(es25.16)',IOSTAT=Reason) newnum
      endif
      if (Reason .EQ. 0) then
         read_well = .true. ! we read something meaningfull from the file
         out_line = trim(adjustl(partline)) ! get the meaningful line out
         if (ABS(newnum) <= 1.0d-10) then
            read(lastpart,*) newnum
         endif
         out_num = newnum ! get the number out
         !print*, 'The number read is', out_num
      else if (Reason .LT. 0) then
         print*, 'No number found in the communication file'
      else if (Reason .GT. 0) then
         print*, 'Given number interpreted as', trim(adjustl(lastpart)), ', it does not match the variable type'
      endif
   else
      !print*, 'EMPTY FILE'
   endif
1111 continue
end subroutine pars_communications_old



END MODULE Output
