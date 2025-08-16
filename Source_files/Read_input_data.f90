! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to read input files:

MODULE Read_input_data
use Universal_constants
use Objects
use Dealing_with_files
use Dealing_with_cdf_files
use Dealing_with_EADL
use Read_numerical_parameters, only: m_input_folder, read_num_pars, read_output_grid_coord, &
                set_MD_step_grid, set_time_grid
use Read_MD_parameters, only: read_supercell_parameters, set_atomic_coordinates, &
                set_atomic_velocities, set_MD_potential
use MD_general_tools, only: get_nearest_neighbors_list
use Periodic_table
use Little_subroutines


implicit none

! All paths to input data and databases are collected here within this module:
character(50) :: m_databases, m_input_minimal, m_input_data, m_numerical_parameters, m_EADL, m_EPDL
character(50) :: m_photon_CS, m_electron_CS, m_ion_CS, m_positron_CS, m_hole_CS
character(50) :: m_raw_EPDL, m_interpolated_Ph_CS, m_folder_materials
character(50) :: m_folder_normalized_CDF, m_folder_DOS, m_folder_MD
character(50) :: m_optical_data, m_fitted_coeffs, m_photon_Compton
character(50) :: m_file_Compton_profiles, m_file_pair_coeffs, m_file_form_factors
character(50) :: m_photon_pair, m_photon_Rayleigh, m_photon_absorption
character(50) :: m_electron_inelast_CS, m_electron_elast_CS, m_electron_Brems_CS
character(50) :: m_ion_inelast_CS, m_holes_inelast_CS, m_holes_elast_CS
character(50) :: m_positron_inelast_CS, m_positron_elast_CS, m_positron_Brems_CS, m_positron_annihilation


! All input folders / directories:
parameter(m_databases = 'Atomic_parameters')		! folder with the periodic table and ENDL databases
parameter(m_photon_CS = 'Photon_cross_sections')	! folder with all photons cross sections and MFPs
parameter(m_electron_CS = 'Electron_cross_sections')	! folder with all electron cross sections and MFPs
parameter(m_hole_CS = 'Valence_hole_cross_sections')	! folder with all valence holes cross sections and MFPs
parameter(m_ion_CS = 'Ion_cross_sections')	! folder with all ion cross sections and MFPs
parameter(m_positron_CS = 'Positron_cross_sections')	! folder with all positron cross sections and MFPs
parameter(m_folder_materials = 'Materials_parameters')	! folder with all material parameters files (chemical formalae, density, etc.)
parameter(m_folder_normalized_CDF = 'Normalized_CDF')	! folder with normalized CDFdata for all elements comvmerted from CSs extracted from EPDL
parameter(m_folder_DOS = 'DOS')	    ! folder with all material DOSes
parameter(m_folder_MD = 'MD_input') ! folder with all parameters for MD simulations
! Databases and input files:
! New format (minimal) input:
parameter(m_input_minimal = 'INPUT.txt')		! file frim INPUT parameters
! Old format (complete) input:
parameter(m_input_data = 'INPUT_DATA.txt')		! file frim INPUT_DATA parameters
parameter(m_numerical_parameters = 'NUMERICAL_PARAMETERS.txt')	! file with all numerical parameters
!parameter(m_EADL = 'EADL2017.all')				! file with EADL database OLD
!parameter(m_EPDL = 'EPDL2017.all')				! file with EPDL database OLD
!parameter(m_EADL = 'EADL2023.ALL')				! file with EADL database OLD
!parameter(m_EPDL = 'EPDL2023.ALL')				! file with EPDL database OLD
parameter(m_EADL = 'EADL2025.ALL')				! file with EADL database
parameter(m_EPDL = 'EPDL2025.ALL')				! file with EPDL database
parameter(m_file_Compton_profiles = 'Compton_profiles.dat')	! database of Compton profiles to be used for calculations of Compton cross sections
parameter(m_file_pair_coeffs = 'Pair_creation_coefficients.dat')	! database with coefficients to be used in photon pair creation cross sections
parameter(m_file_form_factors  = 'Atomic_form_factors.dat')		! database with coefficients used to construct atomic form factors
! Calculated data to be reused as input parameters:
parameter(m_raw_EPDL = 'Raw_absorption')		! part of name of files with photon cross sections extracted from EPDL
parameter(m_interpolated_Ph_CS = 'Interpolated_absorption')	! part of name of files with interpolated photoabsorption CSs, MFPs
parameter(m_optical_data = 'Optical_data')				! part of name of files with optical data extracted from EPDL
parameter(m_fitted_coeffs = 'Fifted_CDF_coefficients')		! part of name of files with fitted CDF coefficients
parameter(m_photon_Compton = 'Photon_Compton_CS')	! part of name of files with photon Compton CS
parameter(m_photon_pair = 'Photon_Pair_CS')			! part of name of files with photon pair creation CS
parameter(m_photon_Rayleigh = 'Photon_Rayleigh_CS')	! part of name of files with photon Rayleigh (elastic) CS
parameter(m_photon_absorption = 'Photon_absorption_CS')	! part of name of files with photon absorption CS
parameter(m_electron_inelast_CS = 'Electron_inelastic_CS')	! part of name of files with electron inelastic CS
parameter(m_electron_elast_CS = 'Electron_elastic_CS')	! part of name of files with electron elastic CS
parameter(m_electron_Brems_CS = 'Electron_Bremsstrahlung_CS')	! part of name of files with electron Bremsstrahlung CS
parameter(m_ion_inelast_CS = '_inelastic_CS')	! part of name of files with ion inelastic CS
parameter(m_positron_inelast_CS = 'Positron_inelastic_CS')	! part of name of files with positron inelastic CS
parameter(m_positron_elast_CS = 'Positron_elastic_CS')	! part of name of files with positron elastic CS
parameter(m_positron_Brems_CS = 'Positron_Bremsstrahlung_CS')	! part of name of files with positron Bremsstrahlung CS
parameter(m_positron_annihilation = 'Positron_annihilation_CS')	! part of name of files with positron annihilation
parameter(m_holes_inelast_CS = 'Holes_inelastic_CS')	! part of name of files with valence hole inelastic CS
parameter(m_holes_elast_CS = 'Holes_elastic_CS')	! part of name of files with valence hole elastic CS


 contains


subroutine Get_targets_parameters(used_target, numpar, Err)
   type(Matter), intent(inout), target :: used_target	! parameters of the target
   type(Num_par), intent(inout), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !---------------------------------------------------------------
   integer, dimension(:), allocatable :: Zat, at_NVB
   real(8), dimension(:), allocatable :: at_percentage, at_masses
   character(3), dimension(:), allocatable :: at_short_names   
   integer, dimension(:), allocatable :: Shl_dsgnr	! EADL shell designator
   character(11), dimension(:), allocatable :: Shell_name	! names of the shells
   real(8), dimension(:), allocatable :: Ip		! [eV] ionization potentials for all shells
   real(8), dimension(:), allocatable :: Ek		! [eV] mean kinetic energy of all shells
   real(8), dimension(:), allocatable :: Ne_shell	! number of electron in each shell
   real(8), dimension(:), allocatable :: Radiat	! [fs] radiative-decay times for all shells
   real(8), dimension(:,:), allocatable :: f_rad	! probability of radiative decay from a given shell
   real(8), dimension(:), allocatable :: Auger	! [fs] Auger-decay times for all shells
   real(8), dimension(:,:,:), allocatable :: f_auger	! probability of Auger decay between two given shells
   real(8), dimension(:,:), allocatable :: Phot_abs_CS_tot	! total cross sections of photoabsorption (to be extracted from EPDL)
   real(8), dimension(:,:,:), allocatable :: Phot_abs_CS_shl	! per shell cross sections of photoabsorption (to be extracted from EPDL)
   real(8), dimension(:), allocatable :: Ip_target  ! [eV] ionization potentials of all elements in the target
   real(8) :: Omega, ksum, fsum
   integer :: Reason, INFO, i, N_elements, j, FN, FN2, k, FN3, m, Npoints, FN4, count_lines, FN5, FN6, FN7, N_size, N_size2
   character(200) :: File_name, File_name2, File_name3, File_name4, File_name5
   character(200) :: chemical_formula, Error_descript, Path, folder_photons, file_photons, command, file_interpolated, File_DOS
   character(1200) :: Error_descript_extended
   character, pointer :: path_sep
   logical :: file_exist, file_opened, read_well
   integer, pointer :: Ntargets
   path_sep => numpar%path_sep
   
   ! 0) Check if the ENDL databases are present, and get their addresses:
   call Check_EPICS_database(path_sep, File_name, File_name2, INFO, Error_descript)	! below
   if (INFO /= 0) then	! there was some error with the ENDL databases
      ! print out the error into the log file:
      call Save_error_details(Err, 1, Error_descript)	! module "Objects"
      !print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
      Error_descript_extended = trim(adjustl(Error_descript))//' TREKIS terminates. '// &
      'Please download the EPICS databases (EADL2025.ALL and EPDL2025.ALL) from https://nds.iaea.org/epics/ '// &
      'or  http://redcullen1.net/HOMEPAGE.NEW/index.htm '// &
      'and place them into the directory: INPUT_DATA\Atomic_parameters'
      call print_error(Error_descript_extended)   ! module "Little_subroutines"
      goto 9992	! skip executing the program, exit the subroutine
   endif
   
   print*, 'Deciphering the chemical formula of the target, reading periodic table...'
   
   numpar%N_sh_tot = 0  ! to start counting number of total shells
   
   ! 1) Interprete the names and read parameters for all target componenets:
   TRGT:do i = 1, used_target%NOC		! number of different components of the target (layers, clusters, etc.)   
      ! Check if the folder with material parameters is present:
      Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
      inquire(DIRECTORY=trim(adjustl(Path)),exist=file_exist)    ! check if input file excists
      if (.not.file_exist) then
         Error_descript = 'Folder '//trim(adjustl(Path))//' not found'
         call Save_error_details(Err, 1, Error_descript)	! module "Objects"
         !print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         Error_descript_extended = trim(adjustl(Error_descript))//', TREKIS terminates'
         call print_error(Error_descript_extended)   ! module "Little_subroutines"
         goto 9992	! skip executing the program, exit the subroutin
      endif
      ! Folder with the particular material:
      Path = trim(adjustl(Path))//path_sep//trim(adjustl(used_target%Material(i)%Name))
      inquire(DIRECTORY=trim(adjustl(Path)),exist=file_exist)    ! check if input file excists
      if (.not.file_exist) then
         Error_descript = 'Folder '//trim(adjustl(Path))//' not found'
         call Save_error_details(Err, 1, Error_descript)	! module "Objects"
         !print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         Error_descript_extended = trim(adjustl(Error_descript))//', TREKIS terminates'
         call print_error(Error_descript_extended)   ! module "Little_subroutines"
         goto 9992	! skip executing the program, exit the subroutin
      endif
      ! Check if the file with the given name is present:
      File_name3 = trim(adjustl(Path))//path_sep//trim(adjustl(used_target%Material(i)%Name))//'.txt'
      inquire(file=trim(adjustl(File_name3)),exist=file_exist) ! check if input file is there
      if (.not.file_exist) then
         Error_descript = 'File '//trim(adjustl(File_name3))//' not found'
         call Save_error_details(Err, 1, Error_descript)	! module "Objects"
         !print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         Error_descript_extended = trim(adjustl(Error_descript))//', TREKIS terminates'
         call print_error(Error_descript_extended)   ! module "Little_subroutines"
         goto 9992	! skip executing the program, exit the subroutin
      endif
      open(newunit = FN, FILE = trim(adjustl(File_name3)), status = 'old', action='read')
      ! Read the chemical fomula from the file:
      call read_material_parameters(FN, File_name3, chemical_formula, used_target%Material(i), Err)	! below
      close (FN)	! done reading material parameters

      ! Make sure parameters are consistent:
      if ((numpar%El_inelast == 1) .or. (numpar%El_inelast == 3)) then  ! check if we have CDF parameteres:
         if (.not.allocated(used_target%Material(i)%CDF_valence%A)) then
            numpar%El_inelast = 5   ! use SPDelta instead
            numpar%Pos_inelast = 5   ! use SPDelta instead
            numpar%H_inelast = 5   ! use SPDelta instead
         endif
      endif

      ! Make sure parameters are consistent:
      if (numpar%El_elastic == 1) then ! for electron
         if (.not.allocated(used_target%Material(i)%CDF_phonon%A)) then    ! there is no CDF in the file
            numpar%El_elastic = 5   ! use SPDelta for electrons instead
         endif
      endif
      if (numpar%Pos_elastic == 1) then ! for positron
         if (.not.allocated(used_target%Material(i)%CDF_phonon%A)) then    ! there is no CDF in the file
            numpar%Pos_elastic = 5   ! use SPDelta for electrons instead
         endif
      endif
      if (numpar%H_elast == 1) then ! for valence hole
         if (.not.allocated(used_target%Material(i)%CDF_phonon%A)) then    ! there is no CDF in the file
            numpar%H_elast = 5   ! use SPDelta for electrons instead
         endif
      endif
      !print*, chemical_formula, used_target%Material(i)%Dens, used_target%Material(i)%CDF_valence%E0(:)
      
      ! Interprete the target chemical formula, and read data for each element from the periodic table:
      Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_databases))
      !call Decompose_compound(Path, used_target%Material(i)%Name, path_sep, INFO, Error_descript, at_num=N_elements, at_numbers=Zat, at_percentage=at_percentage, at_short_names=at_short_names, at_masses=at_masses, at_NVB=at_NVB)	! module "Periodic_table"
      call Decompose_compound(Path, chemical_formula, path_sep, INFO, Error_descript, at_num=N_elements, at_numbers=Zat, at_percentage=at_percentage, at_short_names=at_short_names, at_masses=at_masses, at_NVB=at_NVB)	! module "Periodic_table"
      if (INFO /= 0) then
         call Save_error_details(Err, 4, Error_descript)	! module "Objects"
         !print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         Error_descript_extended = trim(adjustl(Error_descript))//', TREKIS terminates'
         call print_error(Error_descript_extended)   ! module "Little_subroutines"
         goto 9992	! skip executing the program, exit the subroutin
      endif
      ! Now we know how many elements are in this material of the target componen:
      allocate(used_target%Material(i)%Elements(N_elements))	! allocate array of element details
      ! Copy the output into arrays:
      used_target%Material(i)%Elements(:)%Zat = Zat(:)
      used_target%Material(i)%N_Elements = N_Elements
      used_target%Material(i)%Elements(:)%percentage = at_percentage(:)
      used_target%Material(i)%Elements(:)%Name = at_short_names(:)
      used_target%Material(i)%Elements(:)%Mass = at_masses
      used_target%Material(i)%Elements(:)%M = used_target%Material(i)%Elements(:)%Mass * g_amu  ! [kg]
      used_target%Material(i)%Elements(:)%NVB = at_NVB
      ! Get atomic derived parameters:
      used_target%Material(i)%Mean_Mass = g_amu*SUM(used_target%Material(i)%Elements(:)%Mass*used_target%Material(i)%Elements(:)%percentage)/dble(SUM(used_target%Material(i)%Elements(:)%percentage))	! [kg] average atomic mass
      used_target%Material(i)%Mean_Z = SUM(used_target%Material(i)%Elements(:)%Zat*used_target%Material(i)%Elements(:)%percentage)/dble(SUM(used_target%Material(i)%Elements(:)%percentage))	! average atomic number
      used_target%Material(i)%At_Dens = 1.0d-3*used_target%Material(i)%Dens/used_target%Material(i)%Mean_Mass	! [1/cm^3] atomic density
      used_target%Material(i)%N_VB_el = dble(SUM(used_target%Material(i)%Elements(:)%NVB*used_target%Material(i)%Elements(:)%percentage))	! total number of VB electrons in this material [1/molecule]
      ! Debye energy [eV]:
      used_target%Material(i)%E_debye = Debye_energy(used_target%Material(i)%At_Dens, used_target%Material(i)%v_sound) ! module "Little_subroutines"
      ! Einstein energy [eV]:
      used_target%Material(i)%E_eistein = Einstein_energy(used_target%Material(i)%E_debye) ! module "Little_subroutines"
      
!       print*, 'NAMES: ', allocated(at_short_names), at_short_names(:)

      ! Clean up:
      if (allocated(Zat)) deallocate(Zat)
      if (allocated(at_percentage)) deallocate(at_percentage)
      if (allocated(at_masses)) deallocate(at_masses)
      if (allocated(at_NVB)) deallocate(at_NVB)
      if (allocated(at_short_names)) deallocate(at_short_names)
   
!       print*, 'TARGET', i, used_target%Material(i)%N_Elements, used_target%Material(i)%Elements(:)%Name, used_target%Material(i)%Elements(:)%Zat, used_target%Material(i)%Elements(:)%percentage, used_target%Material(i)%Elements(:)%Mass, used_target%Material(i)%Elements(:)%NVB
      
      write(*,'(a)', advance='no')  ' Periodic table read successfully for: '
      do j = 1, size(used_target%Material(i)%Elements(:))
         write(*,'(a,a)', advance='no') used_target%Material(i)%Elements(j)%Name, " "
      enddo
      write(*,'(a)')
!       print*, 'Z: ', used_target%Material(i)%Elements(:)%Zat
      print*, 'Reading EADL database with atomic data...'
      
      !--------------------------------------------------
      ! 2) Read EADL database with atomic parameters for the target atoms
      ! Open EADL database file:
      open(newunit = FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
      inquire(unit=FN,opened=file_opened)    ! check if this file is opened
      if (file_opened) then
         do j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent

            call Read_EADL_data(path_sep, FN, File_name, INFO, used_target%Material(i)%Elements(j)%Zat, Shl_dsgnr=Shl_dsgnr, Shell_name=Shell_name, Ip=Ip, Ek=Ek, Ne_shell=Ne_shell, Radiat=Radiat, f_rad=f_rad, Auger=Auger, f_auger=f_auger)	! module "Dealing_with_EADL"

            if (INFO /= 0) then	! there was some error with the ENDL databases
               ! print out the error into the log file:
               Error_descript = 'Error while reading EADL file: '//trim(adjustl(File_name))
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               !print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
               Error_descript_extended = trim(adjustl(Error_descript))//', TREKIS terminates'
               call print_error(Error_descript_extended)   ! module "Little_subroutines"
               goto 9992	! skip executing the program, exit the subroutin
            endif
            ! Save the data into arrays:
            used_target%Material(i)%Elements(j)%N_shl = size(Shl_dsgnr)	! number of shells in this element
            if (.not.allocated(used_target%Material(i)%Elements(j)%Shl_dsgnr)) allocate(used_target%Material(i)%Elements(j)%Shl_dsgnr(size(Shl_dsgnr)))
            used_target%Material(i)%Elements(j)%Shl_dsgnr = Shl_dsgnr
            if (.not.allocated(used_target%Material(i)%Elements(j)%Shell_name)) allocate(used_target%Material(i)%Elements(j)%Shell_name(size(Shell_name)))
            used_target%Material(i)%Elements(j)%Shell_name = Shell_name
            if (.not.allocated(used_target%Material(i)%Elements(j)%Ip)) allocate(used_target%Material(i)%Elements(j)%Ip(size(Ip)))
            used_target%Material(i)%Elements(j)%Ip= Ip
            if (.not.allocated(used_target%Material(i)%Elements(j)%Ek)) allocate(used_target%Material(i)%Elements(j)%Ek(size(Ek)))
            used_target%Material(i)%Elements(j)%Ek = Ek
            if (.not.allocated(used_target%Material(i)%Elements(j)%Ne_shell)) allocate(used_target%Material(i)%Elements(j)%Ne_shell(size(Ne_shell)))
            used_target%Material(i)%Elements(j)%Ne_shell = Ne_shell
            if (.not.allocated(used_target%Material(i)%Elements(j)%Radiat)) allocate(used_target%Material(i)%Elements(j)%Radiat(size(Radiat)))
            select case(numpar%H_Radiat)
            case default
               used_target%Material(i)%Elements(j)%Radiat = 3.0d23   ! [fs] set infinite, if user set to exclude it
            case(1) ! use EPICS database:
               used_target%Material(i)%Elements(j)%Radiat = Radiat
            end select
            if (.not.allocated(used_target%Material(i)%Elements(j)%f_rad)) allocate(used_target%Material(i)%Elements(j)%f_rad(size(f_rad,1),size(f_rad,2)))
            used_target%Material(i)%Elements(j)%f_rad = f_rad
            if (.not.allocated(used_target%Material(i)%Elements(j)%Auger)) allocate(used_target%Material(i)%Elements(j)%Auger(size(Auger)))
            select case(numpar%H_Auger)
            case default
               used_target%Material(i)%Elements(j)%Auger = 2.0d23   ! [fs] set infinite, if user set to exclude it
            case(1) ! use EPICS database:
               used_target%Material(i)%Elements(j)%Auger = Auger
            end select
            if (.not.allocated(used_target%Material(i)%Elements(j)%f_auger)) allocate(used_target%Material(i)%Elements(j)%f_auger(size(f_auger,1),size(f_auger,2),size(f_auger,3)))
            used_target%Material(i)%Elements(j)%f_auger = f_auger
            if (.not.allocated(used_target%Material(i)%Elements(j)%Compton)) then 
               allocate(used_target%Material(i)%Elements(j)%Compton(size(Auger)))
               used_target%Material(i)%Elements(j)%Compton = 0.0d0	! just to start
            endif
            
            ! Mark all shells as valent or core:
            call tick_valence(used_target%Material(i)%Elements(j)%Ne_shell, dble(used_target%Material(i)%Elements(j)%NVB), used_target%Material(i)%Elements(j)%valent)	! module "Little_subroutines"
             
            ! Count how many core shell there are in total, in all types of atoms:
            used_target%Material(i)%Elements(j)%N_core_shl = COUNT(.not.used_target%Material(i)%Elements(j)%valent(:))  ! number of core shells in this atom
            numpar%N_sh_tot = numpar%N_sh_tot + used_target%Material(i)%Elements(j)%N_core_shl  ! sum them up

            ! Clean up after dealing with this element:
            if (allocated(Shl_dsgnr)) deallocate(Shl_dsgnr)	! EADL shell designator
            if (allocated(Shell_name)) deallocate(Shell_name)	! names of the shells
            if (allocated(Ip)) deallocate(Ip)		! [eV] ionization potentials for all shells
            if (allocated(Ek)) deallocate(Ek)		! [eV] mean kinetic energy of all shells
            if (allocated(Ne_shell)) deallocate(Ne_shell)	! number of electron in each shell
            if (allocated(Radiat)) deallocate(Radiat)	! [fs] radiative-decay times for all shells
            if (allocated(f_rad)) deallocate(f_rad)	! probability of radiative decay from a given shell
            if (allocated(Auger)) deallocate(Auger)	! [fs] Auger-decay times for all shells
            if (allocated(f_auger)) deallocate(f_auger)	! probability of Auger decay for each pairs of shells
            rewind(FN)	! after we read everything for this element, start over for the new element
         enddo ! j = 1, used_target%Material(i)%N_Elements
      else ! could not open EADL
         INFO = 2
         Error_descript = 'Could not open EADL database from: '//trim(adjustl(File_name))
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         !print*, trim(adjustl(Error_descript))//', TREKIS terminates'
         Error_descript_extended = trim(adjustl(Error_descript))//', TREKIS terminates'
         call print_error(Error_descript_extended)   ! module "Little_subroutines"
         goto 9992	! skip executing the program, exit the subroutin
      endif
      
      ! Count how many core shell there are in total, in all types of atoms:
      numpar%N_sh_tot = numpar%N_sh_tot + 1 ! add also valence band as a separate shell
      
      print*, 'EADL database read successfully.'
      print*, 'Reading EPDL database for photon cross sections...'
      
      !--------------------------------------------------
      ! 3) Read EPDL database with atomic parameters for the target atoms
      ! Open EPDL database file:
      open(newunit = FN2, FILE = trim(adjustl(File_name2)), status = 'old', action='read')
      inquire(unit=FN2,opened=file_opened)    ! check if this file is opened
      if (file_opened) then
         do j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent
            ! Check if data for this element were already extracted from the database:
            folder_photons =  trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_photon_CS))
            folder_photons = trim(adjustl(folder_photons))//path_sep//trim(adjustl(used_target%Material(i)%Elements(j)%Name))

!             print*, 'TEST 0: ', i, j, trim(adjustl(used_target%Material(i)%Elements(j)%Name))

            inquire(DIRECTORY=trim(adjustl(folder_photons)),exist=file_exist)    ! check if input file excists
            if (.not.file_exist) then	! to make sure that such a folder is present (even if empty)
              ! Create a new directory for output files:
              command='mkdir '//trim(adjustl(folder_photons))	! to create a folder use this command
              CALL system(command)	! create the folder
            endif
            
            ! Files with raw data:
            file_photons = trim(adjustl(folder_photons))//path_sep//trim(adjustl(m_raw_EPDL))//'_total.dat'
            inquire(file=trim(adjustl(file_photons)),exist=file_exist) ! check if input file is there
            if (file_exist) then	! read from the file, no need to use EPDL
               open(newunit = FN3, FILE = trim(adjustl(file_photons)))
               ! Find out how many data points we have, to allocate the arrays:
               call Count_lines_in_file(FN3, Npoints)	! module "Dealing_with_files"
               ! Now allocate the arrays:
               allocate(Phot_abs_CS_tot(2,Npoints))
               Phot_abs_CS_tot(1,:) = 1.0d30
               Phot_abs_CS_tot(2,:) = -1.0d-12
               allocate(Phot_abs_CS_shl(used_target%Material(i)%Elements(j)%N_shl,2,Npoints))
               Phot_abs_CS_shl(:,1,:) = 1.0d30
               Phot_abs_CS_shl(:,2,:) = -1.0d-12
               ! and read data into them:
               do k = 1, size(Phot_abs_CS_tot,2)
                  read(FN3, '(es,es)',IOSTAT=Reason) Phot_abs_CS_tot(1,k), Phot_abs_CS_tot(2,k)
                  if (Reason < 0) exit	! end of file
               enddo
               close(FN3)
               
               ! Allocate and create energy grid as an object:
               call create_energy_grid(1.0d0, Phot_abs_CS_tot(1,size(Phot_abs_CS_tot,2)), used_target%Material(i)%Elements(j)%Phot_absorption%E, special_points=used_target%Material(i)%Elements(j)%Ip)	! module "Little_subroutines"
               ! Now interpolate the data onto the constructed energy grid:
               allocate(used_target%Material(i)%Elements(j)%Phot_absorption%Total(size(used_target%Material(i)%Elements(j)%Phot_absorption%E)))
               call interpolate_data_on_grid(Phot_abs_CS_tot(1,:), Phot_abs_CS_tot(2,:), used_target%Material(i)%Elements(j)%Phot_absorption%E, used_target%Material(i)%Elements(j)%Phot_absorption%Total) ! module "Little_subroutines"
               
               ! Also save interpolated data, for user to check in case if needed:
               file_interpolated =  trim(adjustl(folder_photons))//path_sep//trim(adjustl(m_interpolated_Ph_CS))//'_total.dat'
               open(newunit = FN4, FILE = trim(adjustl(file_interpolated)))
               do m = 1, size(used_target%Material(i)%Elements(j)%Phot_absorption%E)
                  write(FN4,'(es,es)') used_target%Material(i)%Elements(j)%Phot_absorption%E(m), used_target%Material(i)%Elements(j)%Phot_absorption%Total(m)
               enddo
               close (FN4)
               
               ! Allocate the object-arrays for all shell on the same energy grid:
               allocate(used_target%Material(i)%Elements(j)%Phot_absorption%Per_shell(used_target%Material(i)%Elements(j)%N_shl,size(used_target%Material(i)%Elements(j)%Phot_absorption%E)))
               allocate(used_target%Material(i)%Elements(j)%Phot_absorption%MFP(used_target%Material(i)%Elements(j)%N_shl,size(used_target%Material(i)%Elements(j)%Phot_absorption%E)))
               
               ! Get cross sections for each shell too:
               do m = 1, used_target%Material(i)%Elements(j)%N_shl		! for all shells in this element
                  file_photons = trim(adjustl(folder_photons))//path_sep//trim(adjustl(m_raw_EPDL))//'_'//trim(adjustl(used_target%Material(i)%Elements(j)%Shell_name(m)))//'.dat'
                  inquire(file=trim(adjustl(file_photons)),exist=file_exist) ! check if input file is there
                  if (.not.file_exist) then ! in an extremely unprobable case that some file is missing, just recreate that file (and all the other files) from EPDL
                     deallocate(Phot_abs_CS_tot)
                     deallocate(Phot_abs_CS_shl)
                     goto 9991
                  endif
                  open(newunit = FN3, FILE = trim(adjustl(file_photons)))
                  do k = 1, size(Phot_abs_CS_shl,3)
                     read(FN3, '(es,es)',IOSTAT=Reason) Phot_abs_CS_shl(m,1,k), Phot_abs_CS_shl(m,2,k)
                     if (Reason < 0) exit	! end of file
                  enddo
                  close(FN3)
                  ! Fill the photoabsorption cross sections arrays as objects for each shell:
                  call interpolate_data_on_grid(Phot_abs_CS_shl(m,1,:), Phot_abs_CS_shl(m,2,:), used_target%Material(i)%Elements(j)%Phot_absorption%E, used_target%Material(i)%Elements(j)%Phot_absorption%Per_shell(m,:)) ! module "Little_subroutines"
                  
                  ! Also save interpolated data, for user to check in case if needed:
                  file_interpolated = trim(adjustl(folder_photons))//path_sep//trim(adjustl(m_interpolated_Ph_CS))//'_'//trim(adjustl(used_target%Material(i)%Elements(j)%Shell_name(m)))//'.dat'
                  open(newunit = FN4, FILE = trim(adjustl(file_interpolated)))
                  do k = 1, size(used_target%Material(i)%Elements(j)%Phot_absorption%E)
                     write(FN4,'(es,es)') used_target%Material(i)%Elements(j)%Phot_absorption%E(k), used_target%Material(i)%Elements(j)%Phot_absorption%Per_shell(m,k)
                  enddo
                  close (FN4)
                  
               enddo ! m = 1, used_target%Material(i)%Elements(j)%N_shl		! for all shells in this element
               if (allocated(Phot_abs_CS_tot)) deallocate(Phot_abs_CS_tot)	! probability of Auger decay for each pairs of shells
               if (allocated(Phot_abs_CS_shl)) deallocate(Phot_abs_CS_shl)	! probability of Auger decay for each pairs of shells
            else	! read directly from EPDL
               ! Read cross sections for photons from EPDL database:
9991           call Read_EPDL_rata(path_sep, FN2, File_name, INFO, used_target%Material(i)%Elements(j)%Zat, &
                    used_target%Material(i)%Elements(j)%Shl_dsgnr, Phot_abs_CS_tot, Phot_abs_CS_shl)	! module "Dealing_with_EADL"
               Phot_abs_CS_tot(1,:) = Phot_abs_CS_tot(1,:)*1.0d6	! Convert [MeV] -> [eV]
               Phot_abs_CS_tot(2,:) = Phot_abs_CS_tot(2,:)*1.0d-8	! Convert [b] -> [A^2]
               Phot_abs_CS_shl(:,1,:) = Phot_abs_CS_shl(:,1,:)*1.0d6	! Convert [MeV] -> [eV]
               Phot_abs_CS_shl(:,2,:) = Phot_abs_CS_shl(:,2,:)*1.0d-8	! Convert [b] -> [A^2]
               
               ! Save the data into files:
               file_photons = trim(adjustl(folder_photons))//path_sep//trim(adjustl(m_raw_EPDL))//'_total.dat'

!                print*, 'TEST: ', trim(adjustl(file_photons))
!                print*, 'TEST 2: ', path_sep, ' : ', trim(adjustl(m_raw_EPDL))//'_total.dat'
!                print*, 'TEST 3: ', trim(adjustl(folder_photons))

               open(newunit = FN3, FILE = trim(adjustl(file_photons)))
               do k = 1, size(Phot_abs_CS_tot,2)
                  write(FN3, '(es,es)') Phot_abs_CS_tot(1,k), Phot_abs_CS_tot(2,k)
               enddo
               close(FN3)
               
               ! Allocate and create energy grid as an object:
               call create_energy_grid(1.0d0, Phot_abs_CS_tot(1,size(Phot_abs_CS_tot,2)), used_target%Material(i)%Elements(j)%Phot_absorption%E, special_points=used_target%Material(i)%Elements(j)%Ip)	! module "Little_subroutines"
               allocate(used_target%Material(i)%Elements(j)%Phot_absorption%Total(size(used_target%Material(i)%Elements(j)%Phot_absorption%E)))
               ! Now interpolate the data onto the constructed energy grid:
               call interpolate_data_on_grid(Phot_abs_CS_tot(1,:), Phot_abs_CS_tot(2,:), used_target%Material(i)%Elements(j)%Phot_absorption%E, used_target%Material(i)%Elements(j)%Phot_absorption%Total) ! module "Little_subroutines"
               
               ! Also save interpolated data, for user to check in case if needed:
               file_interpolated =  trim(adjustl(folder_photons))//path_sep//trim(adjustl(m_interpolated_Ph_CS))//'_total.dat'
               open(newunit = FN4, FILE = trim(adjustl(file_interpolated)))
               do m = 1, size(used_target%Material(i)%Elements(j)%Phot_absorption%E)
                   write(FN4,'(es,es)') used_target%Material(i)%Elements(j)%Phot_absorption%E(m), used_target%Material(i)%Elements(j)%Phot_absorption%Total(m)
               enddo
               close (FN4)
               
               ! Allocate the object-arrays for all shell on the same energy grid:
               allocate(used_target%Material(i)%Elements(j)%Phot_absorption%Per_shell(used_target%Material(i)%Elements(j)%N_shl,size(used_target%Material(i)%Elements(j)%Phot_absorption%E)))
               
               ! Save cross sections for each shell too:
               do m = 1, used_target%Material(i)%Elements(j)%N_shl		! for all shells in this element
                  file_photons = trim(adjustl(folder_photons))//path_sep//trim(adjustl(m_raw_EPDL))//'_'//trim(adjustl(used_target%Material(i)%Elements(j)%Shell_name(m)))//'.dat'
                  open(newunit = FN3, FILE = trim(adjustl(file_photons)))
                  do k = 1, size(Phot_abs_CS_shl,3)
                     if (Phot_abs_CS_shl(m,2,k) < 0.0d0) exit
                     write(FN3, '(es,es)') Phot_abs_CS_shl(m,1,k), Phot_abs_CS_shl(m,2,k)
                  enddo
                  close(FN3)
                  ! Fill the photoabsorption cross sections arrays as objects for each shell:
                  !call interpolate_data_on_grid(Phot_abs_CS_tot(1,:), Phot_abs_CS_tot(2,:), used_target%Material(i)%Elements(j)%Phot_absorption%E, used_target%Material(i)%Elements(j)%Phot_absorption%Per_shell(m,:)) ! module "Little_subroutines"
                  call interpolate_data_on_grid(Phot_abs_CS_shl(m,1,:), Phot_abs_CS_shl(m,2,:), used_target%Material(i)%Elements(j)%Phot_absorption%E, used_target%Material(i)%Elements(j)%Phot_absorption%Per_shell(m,:)) ! module "Little_subroutines"
                  
                  ! Also save interpolated data, for user to check in case if needed:
                  file_interpolated = trim(adjustl(folder_photons))//path_sep//trim(adjustl(m_interpolated_Ph_CS))//'_'//trim(adjustl(used_target%Material(i)%Elements(j)%Shell_name(m)))//'.dat'
                  open(newunit = FN4, FILE = trim(adjustl(file_interpolated)))
                  do k = 1, size(used_target%Material(i)%Elements(j)%Phot_absorption%E)
                     write(FN4,'(es,es)') used_target%Material(i)%Elements(j)%Phot_absorption%E(k), used_target%Material(i)%Elements(j)%Phot_absorption%Per_shell(m,k)
                  enddo
                  close (FN4)
               enddo ! m = 1, used_target%Material(i)%Elements(j)%N_shl		! for all shells in this element

               if (allocated(Phot_abs_CS_tot)) deallocate(Phot_abs_CS_tot)	! probability of Auger decay for each pairs of shells
               if (allocated(Phot_abs_CS_shl)) deallocate(Phot_abs_CS_shl)	! probability of Auger decay for each pairs of shells
               rewind(FN2)	! after we read everything for this element, start over for the new element
            endif
         enddo ! j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent
      else ! could not open EPDL
         INFO = 2
         Error_descript = 'Could not open EPDL database from: '//trim(adjustl(File_name2))
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         !print*, trim(adjustl(Error_descript))//', TREKIS terminates'
         Error_descript_extended = trim(adjustl(Error_descript))//', TREKIS terminates'
         call print_error(Error_descript_extended)   ! module "Little_subroutines"
         goto 9992	! skip executing the program, exit the subroutin
      endif ! (file_opened)
      print*, 'EPDL database read successfully.'
      
      
      print*, 'Constructing valence and total grids...'
      !--------------------------------------------------
      ! 4) Since energy grids have special points corresponding to the ionization potential of the element,
      ! different elements have to be combined to construct the grids for the valence band and the total cross sections:
      ! Create a temporary array with all ionization potentials for all elements in the target:
      ! Determine the size of the array
      N_size = 0
      do j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent
         N_size = N_size + size(used_target%Material(i)%Elements(j)%Ip)
      enddo
      if (allocated(Ip_target)) deallocate(Ip_target)
      allocate (Ip_target(N_size))
      Ip_target = 0.0d0
      ! Concotinate all Ip arrays into one:
      N_size = 1
      do j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent
         N_size2 = (N_size-1) + size(used_target%Material(i)%Elements(j)%Ip)    ! numbers of elements are to be added
         Ip_target(N_size:N_size2) = used_target%Material(i)%Elements(j)%Ip(:)
         N_size = N_size2+1   ! for the next loop of the cycle
      enddo
      ! Create energy grid for valence and total cross photoabsorption sections:
      call create_energy_grid(1.0d0, used_target%Material(i)%Elements(1)%Phot_absorption%E(size(used_target%Material(i)%Elements(1)%Phot_absorption%E)), used_target%Material(i)%Ph_absorption_valent%E, special_points=Ip_target)	! module "Little_subroutines"
      allocate(used_target%Material(i)%Ph_absorption_total%E(size(used_target%Material(i)%Ph_absorption_valent%E)))
      used_target%Material(i)%Ph_absorption_total%E = used_target%Material(i)%Ph_absorption_valent%E
      print*, 'Valence and total grids constructed successfully.'

      !print*, 'Val=', used_target%Material(i)%Ph_absorption_valent%E(:)
      !print*, 'Tot=', used_target%Material(i)%Ph_absorption_total%E(:)
      !pause 'PHOTON_VAL'

      
      print*, 'Reading Compton profiles database...'      
      !--------------------------------------------------
      ! 5) Read Compton profiles for the target material:
      ! The database was extracted from:
      ! [Biggs et al., At. Data Nuclear Data Tables, 16 (1975) 201]
      select case (numpar%Ph_Compton)
      case (1)	! Include Compton:
         ! Open a file with Compton profiles:
         Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_databases))
         inquire(DIRECTORY=trim(adjustl(Path)),exist=file_exist)    ! check if input file excists
         if (.not.file_exist) then
            Error_descript = 'Folder '//trim(adjustl(Path))//' not found'
            call Save_error_details(Err, 1, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9992	! skip executing the program, exit the subroutin
         endif
         ! Check if the file with the given name is present:
         File_name4 =  trim(adjustl(Path))//path_sep//trim(adjustl(m_file_Compton_profiles))
         inquire(file=trim(adjustl(File_name4)),exist=file_exist) ! check if input file is there
         if (.not.file_exist) then
            Error_descript = 'File '//trim(adjustl(File_name4))//' not found'
            call Save_error_details(Err, 1, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9992	! skip executing the program, exit the subroutin
         endif
         open(newunit = FN6, FILE = trim(adjustl(File_name4)), status = 'old', action='read')
         ! Read Compton profiles from the database:
         do j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent
            call read_Compton_profiles(FN6, used_target%Material(i)%Elements(j)%Zat, used_target%Material(i)%Elements(j)%Shell_name, used_target%Material(i)%Elements(j)%Compton, Err, read_well)	! below
            if (.not.read_well) goto 9992	! skip executing the program, exit the subroutine
         enddo
         close(FN6)
      case default	! no Compton effect included
      end select
      print*, 'Compton profiles read successfully.'
      
      
      print*, 'Reading pair creation coefficients database...'
      !--------------------------------------------------
      ! 6) Read Pair creation coefficients for the target material:
      ! The database was extracted from:
      ! [Biggs et al., At. Data Nuclear Data Tables, 16 (1975) 201]
      if ( (numpar%Ph_Pairs == 1) .or. (numpar%El_Brems == 1) ) then    ! read Pair creation coefficients:
         ! Open a file with Pair creation coeffs:
         Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_databases))
         inquire(DIRECTORY=trim(adjustl(Path)),exist=file_exist)    ! check if input file excists
         if (.not.file_exist) then
            Error_descript = 'Folder '//trim(adjustl(Path))//' not found'
            call Save_error_details(Err, 1, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9992	! skip executing the program, exit the subroutin
         endif
         ! Check if the file with the given name is present:
         File_name5 =  trim(adjustl(Path))//path_sep//trim(adjustl(m_file_pair_coeffs))
         inquire(file=trim(adjustl(File_name5)),exist=file_exist) ! check if input file is there
         if (.not.file_exist) then
            Error_descript = 'File '//trim(adjustl(File_name5))//' not found'
            call Save_error_details(Err, 1, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9992	! skip executing the program, exit the subroutin
         endif
         open(newunit = FN7, FILE = trim(adjustl(File_name5)), status = 'old', action='read')
         ! Read pair creation coeffs from the database:
         do j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent
            call read_Pair_creation_coeffs(FN7, used_target%Material(i)%Elements(j)%Zat, used_target%Material(i)%Elements(j)%Pair_R, used_target%Material(i)%Elements(j)%Pair_nu_inf, Err, read_well)	! below
            if (.not.read_well) goto 9992	! skip executing the program, exit the subroutine
         enddo
         close(FN7)
      else  ! no pair creation included
      end if
      print*, 'Pair creation coefficients  read successfully.'
      
      
      print*, 'Reading form factors coefficients database...'
      !--------------------------------------------------
      ! 7) Read form factors coefficients for the target material:
      ! The database was extracted from PyPENELOPE source code files (ElementPhysicalData.py):
      ! http://pypenelope.sourceforge.net/models.html
      select case (numpar%Ph_Thomson)
      case (1)	! Include Pair creation:
         ! Open a file with Pair creation coeffs:
         Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_databases))
         inquire(DIRECTORY=trim(adjustl(Path)),exist=file_exist)    ! check if input file excists
         if (.not.file_exist) then
            Error_descript = 'Folder '//trim(adjustl(Path))//' not found'
            call Save_error_details(Err, 1, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9992	! skip executing the program, exit the subroutin
         endif
         ! Check if the file with the given name is present:
         File_name5 =  trim(adjustl(Path))//path_sep//trim(adjustl(m_file_form_factors))
         inquire(file=trim(adjustl(File_name5)),exist=file_exist) ! check if input file is there
         if (.not.file_exist) then
            Error_descript = 'File '//trim(adjustl(File_name5))//' not found'
            call Save_error_details(Err, 1, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9992	! skip executing the program, exit the subroutin
         endif
         open(newunit = FN7, FILE = trim(adjustl(File_name5)), status = 'old', action='read')
         ! Read pair creation coeffs from the database:
         do j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent
            call read_form_factors(FN7, used_target%Material(i)%Elements(j)%Zat, used_target%Material(i)%Elements(j)%form_a, Err, read_well)	! below
            if (.not.read_well) goto 9992	! skip executing the program, exit the subroutine
!             print*, used_target%Material(i)%Elements(j)%Zat, used_target%Material(i)%Elements(j)%form_a
         enddo
         close(FN7)
      case default	! no pair creation included
      end select
      print*, 'Form factors coefficients  read successfully.'
      
   enddo TRGT
   
9992 continue 
   ! Close EADL and EPDL database files (if they are still opened):
   call close_file('close',FN=FN)	! module "Dealing_with_files"
   call close_file('close',FN=FN2)	! module "Dealing_with_files"
   call close_file('close',FN=FN3)
   call close_file('close',FN=FN4)
   call close_file('close',FN=FN5)
   call close_file('close',FN=FN6)
   call close_file('close',FN=FN7)
   ! Clean up:
   nullify(path_sep)
end subroutine Get_targets_parameters



subroutine read_form_factors(FN, Zat, a, Err, read_well)
   integer, intent(in) :: FN, Zat	! file number with the database; atomic number
   real(8), dimension(5), intent(out) :: a	! coefficients used to construct form factors
   type(Error_handling), intent(inout) :: Err	! error log
   logical, intent(inout) :: read_well
   real(8), dimension(:,:), allocatable :: read_data	! to read the coeffs
   real(8), dimension(:), allocatable :: read_vec	! to read the coeffs
   integer :: N_line, N_col, i, j, Reason, count_lines
   character(200) :: Error_descript, temp
   N_col = size(a)
   ! Count how many lines are in the file:
   call Count_lines_in_file(FN, N_line, skip_lines=1)	! module "Dealing_with_files"
   if (Zat > N_line) then	! we don't have this element in our database
      read_well = .false.
      write(temp,'(i4)') Zat
      Error_descript = 'Element #'// trim(adjustl(temp)) //'in read_form_factors not found'
      call Save_error_details(Err, 3, Error_descript)	! module "Objects"
      print*, trim(adjustl(Error_descript)), ', TREKIS terminates'  
   else	! we have this element
      allocate(read_data(N_line, N_col))
      allocate(read_vec(N_col))
      count_lines = 1
      read(FN,*)	! skip the first line with description
      do i = 1, N_line	! read the pair creation coeffs
         !read(FN,*,IOSTAT=Reason) read_data(i,:)
         read(FN,*,IOSTAT=Reason) read_vec(:)   ! read into temporary array
         read_data(i,:) = read_vec(:)   ! save into working array
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i3)') 'In read_form_factors could not read line ', count_lines
            call Save_error_details(Err, 2, Error_descript)	! module "Objects"
            goto 9989
         endif
      enddo
      rewind(FN)	! for future use of the file
      ! Get the coefficients:
      a(:) = read_data(Zat,:)
      ! Clean up at the end:
      deallocate(read_data, read_vec)
   endif
9989 continue
end subroutine read_form_factors



subroutine read_Pair_creation_coeffs(FN, Zat, R, nu_inf, Err, read_well)
   integer, intent(in) :: FN, Zat	! file number with the database; atomic number
   real(8), intent(out) :: R, nu_inf	! coefficients used to construct pair creation cross section
   type(Error_handling), intent(inout) :: Err	! error log
   logical, intent(inout) :: read_well
   real(8), dimension(:,:), allocatable :: read_data	! to read the coeffs
   integer :: N_line, N_col, i, j, Reason, count_lines
   character(200) :: Error_descript, temp
   N_col = 2
   ! Count how many lines are in the file:
   call Count_lines_in_file(FN, N_line, skip_lines=1)	! module "Dealing_with_files"
   if (Zat > N_line) then	! we don't have this element in our database
      read_well = .false.
      write(temp,'(i4)') Zat
      Error_descript = 'Element #'//trim(adjustl(temp)) //'in read_Pair_creation_coeffs not found'
      call Save_error_details(Err, 3, Error_descript)	! module "Objects"
      print*, trim(adjustl(Error_descript)), ', TREKIS terminates'  
   else	! we have this element
      allocate(read_data(N_line, N_col))
      read(FN,*)	! skip the first line with description
      do i = 1, N_line	! read the pair creation coeffs
         !read(FN,*,IOSTAT=Reason) read_data(i,:)
         read(FN,*,IOSTAT=Reason) read_data(i,1), read_data(i,2)
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i3)') 'In read_Pair_creation_coeffs could not read line ', count_lines
            call Save_error_details(Err, 2, Error_descript)	! module "Objects"
            goto 9988
         endif
      enddo
      rewind(FN)	! for future use of the file
      ! Get the coefficients:
      R = read_data(Zat,1)
      nu_inf = read_data(Zat,2)
      ! Clean up at the end:
      deallocate(read_data)
   endif
9988 continue
! print*, 'read_Pair_creation_coeffs', R
end subroutine read_Pair_creation_coeffs



subroutine read_Compton_profiles(FN, Z, Shells, Compton, Err, read_well)
   integer, intent(in) :: FN	! file number of the database (must be already opened)
   integer, intent(in) :: Z	! element number
   character(*), dimension(:), intent(in) :: Shells	! shell names for the given element
   real(8), dimension(:), intent(inout) :: Compton	! compton profiles for each atomic shell
   type(Error_handling), intent(inout) :: Err	! error log
   logical, intent(inout) :: read_well
   real(8), dimension(:,:), allocatable :: read_data	! to read the compton profiles database
   real(8), dimension(:), allocatable :: read_vec	! to read the compton profiles database
   character(3), dimension(:), allocatable :: shell_names	! to identify the shells that are listed in the database
   integer :: N_line, N_col, i, j, Reason, count_lines
   character(200) :: Error_descript, temp
   ! Count how many lines are in the file:
   call Count_lines_in_file(FN, N_line, skip_lines=2)	! module "Dealing_with_files"
    if (Z > N_line) then	! we don't have this element in our database
      read_well = .false.
      write(temp,'(i4)') Z
      Error_descript = 'Element #'//trim(adjustl(temp)) //'in read_Compton_profiles not found'
      call Save_error_details(Err, 3, Error_descript)	! module "Objects"
      print*, trim(adjustl(Error_descript)), ', TREKIS terminates'  
   else	! we have this element
      ! Count how many columns are in the file:
      call Count_columns_in_file(FN, N_col, skip_lines=2)	! module "Dealing_with_files"
      allocate(shell_names(N_col))
      allocate(read_data(N_line, N_col))
      allocate(read_vec(N_col))
      read(FN,*) shell_names(:)	! first line contains shell names
      read(FN,*)	! skip the second line with shell names in another convention, not used here
      do i = 1, N_line	! read the compton databases
         !read(FN,*,IOSTAT=Reason) read_data(i,1:N_col)
         read(FN,*,IOSTAT=Reason) read_vec  ! read into temporary array
         read_data(i,:) = read_vec(:)   ! save into the working array
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i3)') 'In read_Compton_profiles could not read line ', count_lines
            call Save_error_details(Err, 2, Error_descript)	! module "Objects"
            goto 9987
         endif
      enddo
      rewind(FN)	! for future use of the file
      ! Identify the shells:
      do i = 1, size(Shells)
         do j = 1, size(shell_names)	! find the same shell:
            if (trim(adjustl(shell_names(j)(1:2))) == trim(adjustl(Shells(i)(1:2))) ) then
               Compton(i) = read_data(Z,j)	! save the Compton profile for this shell of this element
               exit
            endif
         enddo
      enddo
      ! Clean up at the end:
      deallocate(shell_names, read_data, read_vec)
   endif
9987 continue
end subroutine read_Compton_profiles




subroutine read_material_parameters(FN, File_name, chemical_formula, Material, Err)
   integer, intent(in) :: FN	! file from where to read
   character(*), intent(in) :: File_name	! file name
   character(*), intent(out) :: chemical_formula	! chemical formula of a target material
   type(Target_atoms), intent(inout):: Material	!material parameters of each target that it's constructed of
   type(Error_handling), intent(inout) :: Err	! error log
   !-----------------------------
   integer :: count_lines, Reason, N
   character(200) :: Error_descript, temp_ch
   logical :: read_well, valence_present, phonon_present
   
   count_lines = 0	! to start counting lines in the file

   read(FN,*,IOSTAT=Reason) chemical_formula
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9990
   endif
   
   read(FN,*,IOSTAT=Reason) Material%Dens	! material density [g/cm^3]
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9990
   endif
   
   read(FN,*,IOSTAT=Reason) Material%DOS%Egap	! band gap [eV]
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9990
   endif
   
   read(FN,*,IOSTAT=Reason) Material%DOS%E_f		! [eV] Fermi energy
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9990
   endif
   
   read(FN,*,IOSTAT=Reason) Material%v_sound    ! [m/s] speed of sound
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9990
   endif
   
   ! Default values to start with:
   Material%me_eff = 1.0d0  ! default CB effective mass (for single-pole CDF)
   valence_present = .false.
   phonon_present= .false.
   Material%DOS%E_VB_bottom = 0.0d0
   Material%DOS%E_VB_top = 0.0d0
   Material%DOS%E_CB_bottom = 0.0d0
   Material%DOS%E_CB_top = 0.0d0
   ! Set default barrier values:
   call set_default_barrier(Material%Surface_barrier%Work_func, Material%Surface_barrier%Surf_bar, Material%Surface_barrier%Bar_height)    ! below
   
   ! Proceed reading the file to change default values for the user-defined ones:
   do while (read_well)	! check all possible entries (valence band and phonons)
      read(FN,*,IOSTAT=Reason) temp_ch  ! check whether there are CDF parameters provided for valence or phonon system
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         exit	! nothing more to read from the file
      else
         select case (trim(adjustl(temp_ch)))
         
         case ('EFFECTIVE_MASS', 'Effective_mass', 'me_eff', 'ME_EFF', 'Eff_mass')
            ! Next line defines the effective electron mass in CB in units of free-electron mass [me]
            read(FN,*,IOSTAT=Reason) Material%me_eff
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9990
            endif
         case ('VALENCE', 'Valence', 'valence', 'VAL', 'Val', 'val', 'Valent', 'valent', 'VALENT')
            read(FN,*,IOSTAT=Reason) N  ! how many CDF oscillators are used
            ! Next lines will define Ritchie CDF parameters for the valence band:
            if (N > 0) then
               valence_present = .true.
               allocate(Material%CDF_valence%A(N))
               allocate(Material%CDF_valence%E0(N))
               allocate(Material%CDF_valence%Gamma(N))
               call read_CDF_from_file(FN, Material%CDF_valence, N, read_well)	! below
               if (.not.read_well) then
                  deallocate(Material%CDF_valence%A, Material%CDF_valence%E0, Material%CDF_valence%Gamma)
                  print*, 'Could not read CDF parameters for Valence band from file '//trim(adjustl(File_name))
                  print*, 'Using atomic approximation'
                  read_well = .true.    ! to continue reading from the file
               else
                  print*, 'Valence band CDF parameters read from file successfully'
               endif
            endif
         case ('PHONONS', 'PHONON', 'Phonons', 'Phonon', 'phonons', 'phonon', 'PHON', 'Phon', 'phon')
            read(FN,*,IOSTAT=Reason) N  ! how many CDF oscillators are used
            ! Next lines will define Ritchie CDF parameters for the phonons:
            if (N > 0) then
               phonon_present = .true.
               allocate(Material%CDF_phonon%A(N))
               allocate(Material%CDF_phonon%E0(N))
               allocate(Material%CDF_phonon%Gamma(N))
               call read_CDF_from_file(FN, Material%CDF_phonon, N, read_well)   ! below
               if (.not.read_well) then
                  deallocate(Material%CDF_phonon%A, Material%CDF_phonon%E0, Material%CDF_phonon%Gamma)
                  print*, 'Could not read CDF parameters for Phonons from file '//trim(adjustl(File_name))
                  print*, 'Using atomic approximation'
                  read_well = .true.    ! to continue reading from the file
               else
                  print*, 'Phonons CDF parameters read from file successfully'
               endif
            endif
         case ('DOS', 'DOs', 'Dos', 'dos')
            ! Next line defines the name of the file containing DOS for this material:
            read(FN,*,IOSTAT=Reason) Material%DOS_file   ! DOS file name
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9990
            endif
         case ('BANDS', 'Bands', 'bands')
            ! Next line defines parameters of the free-electron DOS to be constructed,
            ! in the following order:
            ! bottom and top of the valence band, bottom and top of the conduction band [eV]
            ! (must be in accord with gap and Fermi energies)
            read(FN,*,IOSTAT=Reason) Material%DOS%E_VB_bottom, Material%DOS%E_VB_top, Material%DOS%E_CB_bottom, Material%DOS%E_CB_top
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9990
            endif
         case ('BARRIER', 'Barrier', 'barrier', 'SURFACE', 'Surface', 'surface')
            ! Next line defines parameters of the surface barrier for particle emission:
            ! Must be in this order:  Work function [eV] ; surface barrier length [A] ; barrier height for electron emission [eV]
            read(FN,*,IOSTAT=Reason) Material%Surface_barrier%Work_func, Material%Surface_barrier%Surf_bar, Material%Surface_barrier%Bar_height
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9990
            endif
         
         end select
      endif
   enddo
   ! Print out warning messages:
   if (.not.valence_present) then
      print*, 'Could not find CDF parameters for Valence band in the file '//trim(adjustl(File_name))
   endif
   if (.not.phonon_present) then
      print*, 'Could not find CDF parameters for Phonons in the file '//trim(adjustl(File_name))
   endif
   
9990 continue
end subroutine read_material_parameters



pure subroutine set_default_barrier(Work_func, Surf_bar, Bar_height)
   real(8), intent(inout) :: Work_func, Surf_bar, Bar_height
   ! Default values for mateirals, for which parameters are unknown:
   Work_func = 4.0d0   ! [eV]
   Surf_bar = 10.0d0    ! [A]
   Bar_height = 6.0d0   ! [eV]
end subroutine set_default_barrier


subroutine read_CDF_from_file(FN2, CDF, N_CDF, read_well)
   integer, intent(in) :: FN2	! file number to read CDF coeffs from
   type(Ritchi_CDF), intent(inout) :: CDF		! parameters entering Ritchi-CDF
   integer, intent(in) :: N_CDF	! number of CDF oscillators
   logical, intent(inout) :: read_well	! check for errors during file reading
   integer :: k, count_lines, Reason
   count_lines = 0
   do k = 1, N_CDF	! all coefficients
      read(FN2, * ,IOSTAT=Reason) CDF%A(k), CDF%E0(k), CDF%Gamma(k)
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
!          print*, 'Error in read_CDF_from_file: problem in line', count_lines
!          pause 'read_CDF_from_file'
         exit	! semething wen wrong while reading the file
      endif
   enddo
end subroutine read_CDF_from_file




subroutine Check_EPICS_database(path_sep, File_name, File_name2, INFO, Error_descript)
    character(1), intent(in) :: path_sep
    character(200), intent(inout) :: File_name
    character(200), intent(inout) :: File_name2
    character(200), intent(inout) :: Error_descript
    integer, intent(out) :: INFO
    character(200) :: Folder_name
    logical file_exist
    Error_descript = ''	! no error at the beginning
    Folder_name = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_databases))  ! here we keep databases
    File_name = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//trim(adjustl(m_EADL))
    File_name2 = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//trim(adjustl(m_EPDL))
    inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
    if (.not. file_exist) then
       write(Error_descript,'(a,a,a)') 'File ',trim(adjustl(File_name)),' is not found.'
       INFO = 1
    else
       INFO = 0
    endif
    inquire(file=trim(adjustl(File_name2)),exist=file_exist) ! check if input file is there
    if (.not. file_exist) then
       write(Error_descript,'(a,a,a)') 'File ',trim(adjustl(File_name2)),' is not found.'
       INFO = 2
    endif
end subroutine Check_EPICS_database


!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
! Call subroutines to read MD parameters:
 
subroutine Read_MD_input(used_target, numpar, MD_atoms, MD_supce, MD_pots, Err)
   type(Matter), intent(in) :: used_target	! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(Atom), dimension(:), allocatable, intent(inout) :: MD_atoms     ! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   type(MD_potential), dimension(:,:), allocatable, intent(inout) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(Error_handling), intent(inout) :: Err	! error log
   !---------------------------
   character(250) :: Folder_name, Path_to_MD, Path_to_Periodic
   
   ! Read all MD-required parameters, if user requested an MD run:
   if (numpar%DO_MD) then

      numpar%MD_integrator = 0  ! default index of MD integrator (velocity Verlet)

      ! Define the path to MD parameter:
      Path_to_MD = trim(adjustl(m_input_folder))//numpar%path_sep//trim(adjustl(m_folder_MD))// &
                    numpar%path_sep//trim(adjustl(used_target%Name))//numpar%path_sep
      Path_to_Periodic = trim(adjustl(m_input_folder))//numpar%path_sep//trim(adjustl(m_databases))
      
      ! Each parameters set from different files:
      ! 1) Supercell size and boundary conditions:
      call read_supercell_parameters(Path_to_MD, numpar, MD_supce, Err)  ! module "Read_MD_parameters"
      
      ! 2) Read or set atomic coordinates:
      call set_atomic_coordinates(numpar, Path_to_Periodic, Path_to_MD, MD_supce, MD_atoms, Err)  ! module "Read_MD_parameters"
      
      ! 3) Read or set atomic velocities:
      call set_atomic_velocities(used_target, numpar, Path_to_MD, MD_atoms, Err)  ! module "Read_MD_parameters"
      
      ! 4) Read parameters of interatomic potentials (assume target #1 for now, may be changed later):
      call set_MD_potential(Path_to_MD, numpar, used_target, 1, MD_pots, Err)  ! module "Read_MD_parameters"

      ! 5) Set the list of nearest neighbors for all atoms:
      call get_nearest_neighbors_list(MD_atoms, MD_supce, MD_pots, numpar) ! module "MD_general_tools"
      
      ! 6) If requested, conform simulation box and MD supercell:
      call conform_box_and_MDsupercell(numpar, MD_supce)    ! below
   endif

   ! Whether there is MD or not, MC-MD info-exchange arrays must be allocated:
   if (.not.allocated(MD_supce%E_e_at_from_MC)) then
      allocate(MD_supce%E_e_at_from_MC(1,1,1), source = 0.0d0)  ! energy elastically transferred from e to atoms from within MC module
   endif
   if (.not.allocated(MD_supce%E_h_at_from_MC)) then
      allocate(MD_supce%E_h_at_from_MC(1,1,1), source = 0.0d0)  ! energy elastically transferred from h to atoms from within MC module
   endif
   if (.not.allocated(MD_supce%E_p_at_from_MC)) then
      allocate(MD_supce%E_p_at_from_MC(1,1,1), source = 0.0d0)  ! energy elastically transferred from p to atoms from within MC module
   endif
   if (.not.allocated(MD_supce%E_e_from_MC)) then
      allocate(MD_supce%E_e_from_MC(1,1,1), source = 0.0d0)   ! energy of electrons below cut-off from within MC module
   endif
   if (.not.allocated(MD_supce%E_h_from_MC)) then
      allocate(MD_supce%E_h_from_MC(1,1,1), source = 0.0d0)   ! energy of holes below cut-off from within MC module
   endif
end subroutine Read_MD_input



subroutine conform_box_and_MDsupercell(numpar, MD_supce)
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters

   ! if the simulation box boundaries should be taken MD supercell:
   if (numpar%reset_from_MD(1)) then   ! [A] coordinates of the left and right ends of the simulation box along X
      numpar%box_start_x = MD_supce%x_start
      numpar%box_end_x = MD_supce%x_end
   endif

   if (numpar%reset_from_MD(2)) then   ! [A] coordinates of the left and right ends of the simulation box along Y
      numpar%box_start_y = MD_supce%y_start
      numpar%box_end_y = MD_supce%y_end
   endif

   if (numpar%reset_from_MD(3)) then   ! [A] coordinates of the left and right ends of the simulation box along Z
      numpar%box_start_z = MD_supce%z_start
      numpar%box_end_z = MD_supce%z_end
   endif
end subroutine conform_box_and_MDsupercell



!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! Read the main input file:
subroutine Read_input(used_target, numpar, bunch, Err)
   type(Matter), intent(inout) :: used_target	! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(Radiation_param), dimension(:), allocatable, intent(inout) :: bunch	! incomming radiation
   type(Error_handling), intent(inout) :: Err	! error log
   !---------------------------
   integer :: FN, FN2, FN3
   character(200) :: Error_descript, read_line, Folder_name, File_name, File_name_2, File_name_new
   logical :: file_exists, file_opened
   Error_descript = ''	! start with no error
   !---------------------------
   ! Read path separator from environment (setting your OS at the same time):
    call Path_separator(numpar%path_sep)	! module "Dealing_with_files"
   !---------------------------
   ! Define initial files:
   Folder_name = trim(adjustl(m_input_folder))//numpar%path_sep	! folder with all input files
   numpar%input_path = Folder_name	! save the address with input files



   !---------------------------
   ! Reading new format of input file using Fortran Namelists (still in the testing mode!):
   call set_defaults(used_target, bunch, numpar)  ! to start with, set all default flags; below

   numpar%new_input_format = .false. ! to start with, assumed new format wasn't used
   File_name_new = trim(adjustl(Folder_name))//trim(adjustl(m_input_minimal))
   print*, 'Reading file: '//trim(adjustl(File_name_new))
   FN3 = 103
   ! Check if files exist and open them:
   call open_file(FN3, File_name_new, Error_descript, status = 'old', action='read')	! module "Dealing_with_files"
   if (LEN(trim(adjustl(Error_descript))) > 0) then	! it means, some error occured
      Error_descript = ''    ! renew it, and just skip the file
      print*, 'Could not read it, checking if old format of input exists...'
   else
      call Read_single_file_input(FN3, File_name_new, used_target, bunch, numpar, Err)   ! below
      numpar%new_input_format = .true. ! new format was used
   endif
   call close_file('close', FN=FN3)	! module "Dealing_with_files"


   !---------------------------
   ! Check old format of input:
   if (.not.numpar%new_input_format) then ! only if didn't read the data in new format, use the old one:
      FN = 101
      FN2 = 102
      File_name = trim(adjustl(Folder_name))//trim(adjustl(m_input_data))	!'INPUT_DATA.txt'
      File_name_2 = trim(adjustl(Folder_name))//trim(adjustl(m_numerical_parameters))	!'NUMERICAL_PARAMETERS.txt'
      !---------------------------
      ! Check if files exist and open them:
      call open_file(FN, File_name, Error_descript, status = 'old', action='read')	! module "Dealing_with_files"
      if (LEN(trim(adjustl(Error_descript))) > 0) then	! it means, some error occured
         ! print out the error into the log file:
         call Save_error_details(Err, 1, Error_descript)	! module "Objects"
         print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         goto 9999	! skip executing the program, exit the subroutine
      endif
   
      call open_file(FN2, File_name_2, Error_descript, status = 'old', action='read')	! module "Dealing_with_files"
      if (LEN(trim(adjustl(Error_descript))) > 0) then	! it means, some error occured
         ! print out the error into the log file:
         call Save_error_details(Err, 1, Error_descript)	! module "Objects"
         print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         goto 9999	! skip executing the program, exit the subroutine
      endif
      !---------------------------
      ! Read the file with input data:
      print*, 'Reading file: '//trim(adjustl(File_name))
      call read_input_parameters(FN, File_name, used_target, bunch, Err)	! see below
      !---------------------------
      ! Read the file with all numerical parameters:
      print*, 'Reading file: '//trim(adjustl(File_name_2))
      call read_num_pars(FN2, File_name_2, numpar, Err)	! module "Read_numerical_parameters"
   endif ! (.not.numpar%new_input_format)
   
   !---------------------------
   ! close opened input files:
9999 continue
   call close_file('close', FN=FN)	! module "Dealing_with_files"
   call close_file('close', FN=FN2)	! module "Dealing_with_files"
end subroutine Read_input



! If all input data are provided in a single file (new format, not all data needed to be provided):
subroutine Read_single_file_input(FN, File_name, used_target, bunch, numpar, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Matter), intent(inout) :: used_target	! parameters of the target
   type(Radiation_param), dimension(:), allocatable, intent(inout) :: bunch	! incomming radiation
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !----------------------
   integer :: Reason, count_lines, i, N_rp, i_rp
   character(200) :: Error_descript, read_line
   logical :: read_well, stay_in_the_loop

   Error_descript = ''  ! to start with
   count_lines = 0	! to start counting lines in the file

   !----------------------
   ! Read all line in file, until find a block keyword:

   !----------------------
   ! I) Mandatory sections:
   !----------------------

   stay_in_the_loop = .true.  ! to start with
   ML:do while (stay_in_the_loop)   ! check all possible entries
      read(FN,*,IOSTAT=Reason) read_line  ! read line to interprete it
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then     ! probably, end of the file reached
         stay_in_the_loop = .false. ! time to leave
         exit ML  ! exit the main loop
      endif


      ! Interprete the line:
      select case (trim(adjustl(read_line)))
      !----------------------
      case ('TARGET', 'Target', 'target', '::: TARGET :::', '::: Target :::', '::: target :::')
         ! 1) Material definition
         call read_target_definition(FN, File_name, used_target, count_lines, read_well, Err)    ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('RADIATION', 'Radiation', 'radiation', '::: RADIATION :::', '::: Radiation :::', '::: radiation :::')
         ! 2) Radiation definition
         call read_radiation_definition(FN, File_name, bunch, count_lines, read_well, Err) ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('::: NUMERICAL PARAMETERS :::', '::: NUMERICS:::', 'NUMERICAL PARAMETERS', 'numerical parameters', 'NUMERICS', 'numerics')
         ! 3) Mandatory numerical parameters
         call read_numerics_definition(FN, File_name, numpar, count_lines, read_well, Err)    ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('OUTPUT', 'Output', 'output', '::: OUTPUT DATA :::', '::: output data :::', '::: OUTPUT :::', '::: output:::')
         ! 4) Mandatory output
         call read_mandatory_output(FN, File_name, numpar, count_lines, read_well, Err)    ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('ELECTRONS', 'Electrons', 'electrons', '::: MODELS FOR ELECTRONS :::')
         ! 5) Optional numerical details for electrons
         call read_optional_electrons_numerics(FN, File_name, numpar, count_lines, read_well, Err)    ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('PHOTONS', 'Photons', 'photons', '::: MODELS FOR PHOTONS :::')
         ! 6) Optional numerical details for photons
         call read_optional_photons_numerics(FN, File_name, numpar, count_lines, read_well, Err)    ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('SHI', 'SHIs', 'Ions', 'IONS', 'ions', '::: MODELS FOR SHI :::')
         ! 7) Optional numerical details for SHIs
         call read_optional_SHI_numerics(FN, File_name, numpar, count_lines, read_well, Err)    ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('Positrons', 'POSITRONS', 'positrons', '::: MODELS FOR POSITRONS :::')
         ! 8) Optional numerical details for positrons
         call read_optional_positrons_numerics(FN, File_name, numpar, count_lines, read_well, Err)    ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('Holes', 'HOLES', 'holes', '::: MODELS FOR HOLES :::')
         ! 9) Optional numerical details for positrons
         call read_optional_holes_numerics(FN, File_name, numpar, count_lines, read_well, Err)    ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('CDF', 'cdf', '::: MODELS FOR CDF :::')
         ! 10) Optional numerical details for CDF
         call read_optional_CDF_numerics(FN, File_name, numpar, count_lines, read_well, Err)    ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('Quenching', 'QUENCHING', 'quenching', 'MD', '::: MODELS FOR MD :::')
         ! 11) Optional numerical details for CDF
         call read_optional_MD_numerics(FN, File_name, numpar, count_lines, read_well, Err)    ! below
         if (.not. read_well) then
            stay_in_the_loop = .false. ! time to leave
            exit ML  ! exit the main loop
         endif

      !----------------------
      case ('Optional', 'OPTIONAL', 'optional')
         ! All mandatory optional were read, proceed to the optional
         exit ML  ! exit the main loop

      end select

   enddo ML
   if (Err%Err) return    ! if an error occured while reading input file, terminate the program

   !----------------------
   ! II) Optional optput section:
   !----------------------
   ! Optional output
   call read_output_grid_coord(FN, File_name, numpar, Err, count_lines) ! module "Read_numerical_parameters"

end subroutine Read_single_file_input



subroutine read_optional_MD_numerics(FN, File_name, numpar, count_lines, read_well, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout) :: numpar	! all numerical parametersn
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   integer :: Reason
   character(200) :: Error_descript

   ! Use quenching in MD (T=true; F=false), when to start [fs], how often nullify velocities [fs]:
   read(FN,*,IOSTAT=Reason) numpar%do_quenching, numpar%t_quench_start, numpar%dt_quench
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
   numpar%t_quench_run = 0.0d0  ! starting

end subroutine read_optional_MD_numerics



subroutine read_optional_CDF_numerics(FN, File_name, numpar, count_lines, read_well, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout) :: numpar	! all numerical parametersn
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   integer :: Reason
   character(200) :: Error_descript

   ! CDF model: 0 = Drude / Lindhard CDF, 1=Mermin CDF, 2=Full conserving CDF (not ready)
   read(FN,*,IOSTAT=Reason) numpar%CDF_model
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! target dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie; effective mass [in me] (0=effective mass from DOS of VB; -1=free-electron):
   read(FN,*,IOSTAT=Reason) numpar%CDF_dispers, numpar%CDF_m_eff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Include plasmon integration limit (0=no, 1=yes):
   read(FN,*,IOSTAT=Reason) numpar%CDF_plasmon
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Coefficient where to switch from Ritchie to Delta CDF: E = k * Wmin (INELASTIC):
   read(FN,*,IOSTAT=Reason) numpar%CDF_Eeq_factor
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Coefficient where to switch from Ritchie to Delta CDF: E = k * Wmin (ELASTIC):
   read(FN,*,IOSTAT=Reason) numpar%CDF_Eeq_elast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

    ! Flag to use for target atoms Zeff (set 0), or Z=1 (set 1):
   read(FN,*,IOSTAT=Reason) numpar%CDF_elast_Zeff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! effective number of grid points for inelastic cross section integration over energy
   read(FN,*,IOSTAT=Reason) numpar%CDF_int_n_inelastE
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! effective number of grid points for inelastic cross section integration over momentum
   read(FN,*,IOSTAT=Reason) numpar%CDF_int_n_inelastQ
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! effective number of grid points for elastic cross section integration over energy
   read(FN,*,IOSTAT=Reason) numpar%CDF_int_n_elastE
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! effective number of grid points for elastic cross section integration over momentum
   read(FN,*,IOSTAT=Reason) numpar%CDF_int_n_elastQ
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
end subroutine read_optional_CDF_numerics


subroutine read_optional_holes_numerics(FN, File_name, numpar, count_lines, read_well, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout) :: numpar	! all numerical parametersn
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   integer :: Reason
   character(200) :: Error_descript

   ! Auger decays:  0=excluded, 1=EADL
   read(FN,*,IOSTAT=Reason) numpar%H_Auger
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Radiative decays: 0=excluded, 1=EADL
   read(FN,*,IOSTAT=Reason) numpar%H_Radiat
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [me] effective valence hole mass in units of electron mass
   read(FN,*,IOSTAT=Reason) numpar%H_m_eff ! [me] effective valence hole mass (-1=from DOS; 0=free electron; >0=fixed mass in m_e)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Valence hole inelastic scattering: 0=excluded, 1=CDF, 2=BEB, 3=Delta
   read(FN,*,IOSTAT=Reason) numpar%H_inelast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Valence hole elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF (NOT READY YET!), 5=SP-CDF
   read(FN,*,IOSTAT=Reason) numpar%H_elast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [eV] Cut-off energy (holes with lower energies are excluded from calculation):
   read(FN,*,IOSTAT=Reason) numpar%H_Cutoff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
end subroutine read_optional_holes_numerics




subroutine read_optional_positrons_numerics(FN, File_name, numpar, count_lines, read_well, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout) :: numpar	! all numerical parametersn
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   integer :: Reason
   character(200) :: Error_descript

   ! Positron inelastic scattering: 0=excluded, 1=delta
   read(FN,*,IOSTAT=Reason) numpar%Pos_inelast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Positron elastic scattering: 0=excluded, 1=delta
   read(FN,*,IOSTAT=Reason) numpar%Pos_elastic
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Positron Bremsstrahlung scattering: 0=excluded, 1=delta
   read(FN,*,IOSTAT=Reason) numpar%Pos_Brems
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Positron annihilation: 0=excluded, 1= Heitler
   read(FN,*,IOSTAT=Reason) numpar%Pos_annih
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [eV] Cut-off energy (electrons with lower energies are excluded from calculation):
   read(FN,*,IOSTAT=Reason) numpar%Pos_Cutoff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
end subroutine read_optional_positrons_numerics



subroutine read_optional_SHI_numerics(FN, File_name, numpar, count_lines, read_well, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout) :: numpar	! all numerical parametersn
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   integer :: Reason
   character(200) :: Error_descript

   ! SHI inelastic scattering: 0=excluded, 1:3=delta, 4=nonrelativ.CDF (DO NOT USE!), 5=SPdelta
   read(FN,*,IOSTAT=Reason) numpar%SHI_inelast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Charge state: 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande, 4=fixed Zeff, 5=charge exchange:
   read(FN,*,IOSTAT=Reason) numpar%SHI_ch_st
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! SHI charge shape: 0=point-like charge; 1=Brandt-Kitagawa ion:
   read(FN,*,IOSTAT=Reason) numpar%SHI_ch_shape
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [eV] Cut-off energy (SHIs with lower energies are excluded from calculation):
   read(FN,*,IOSTAT=Reason) numpar%SHI_Cutoff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
end subroutine read_optional_SHI_numerics



subroutine read_optional_photons_numerics(FN, File_name, numpar, count_lines, read_well, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout) :: numpar	! all numerical parametersn
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   integer :: Reason
   character(200) :: Error_descript

   ! Photoabsorption CSs: 0=excluded, 1=CDF, 2=EPDL:
   read(FN,*,IOSTAT=Reason) numpar%Ph_absorb
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Compton effect: 0=excluded, 1=PENELOPE
   read(FN,*,IOSTAT=Reason) numpar%Ph_Compton
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Thomson / Rayleigh scattering: 0=excluded, 1=PENELOPE
   read(FN,*,IOSTAT=Reason) numpar%Ph_Thomson
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Electron-positron pair creation: 0=excluded, 1=included, 2=included with 6) Landau-Pomeranchuk-Migdal suppression effect:
   read(FN,*,IOSTAT=Reason) numpar%Ph_Pairs
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Photonuclear Physics: 0=excluded, 1=included:
   read(FN,*,IOSTAT=Reason) numpar%Ph_Nucl
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [eV] Cut-off energy (photons with lower energies are excluded from calculation):
   read(FN,*,IOSTAT=Reason) numpar%Ph_Cutoff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [A] Effective photon attenuation length (<0 means do not use effective one, use real one from EPDL):
   read(FN,*,IOSTAT=Reason) numpar%Ph_att_eff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
end subroutine read_optional_photons_numerics




subroutine read_optional_electrons_numerics(FN, File_name, numpar, count_lines, read_well, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout) :: numpar	! all numerical parametersn
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   integer :: Reason
   character(200) :: Error_descript


   ! MC or MD target model: 0=MC, 1=MD
   read(FN,*,IOSTAT=Reason) numpar%MC_vs_MD
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Include forces and fields among electrons: 0=exclude, 1=include
   read(FN,*,IOSTAT=Reason) numpar%El_forces
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! inelastic scattering: 0=excluded, 1=relativ.CDF, 2=RBEB, 3=delta, 4=nonrelativ.CDF (DO NOT USE!), 5=SPdelta
   read(FN,*,IOSTAT=Reason) numpar%El_inelast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   read(FN,*,IOSTAT=Reason) numpar%El_elastic
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Bremsstrahlung: 0=excluded, 1=BHW
   read(FN,*,IOSTAT=Reason) numpar%El_Brems
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Cherenkov radiation: 0=excluded, 1= ...
   read(FN,*,IOSTAT=Reason) numpar%El_Cheren
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [eV] Cut-off energy (electrons with lower energies are excluded from calculation):
   read(FN,*,IOSTAT=Reason) numpar%El_Cutoff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

end subroutine read_optional_electrons_numerics



subroutine read_mandatory_output(FN, File_name, numpar, count_lines, read_well, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout) :: numpar	! all numerical parametersn
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   integer :: temp, Reason
   character(200) :: Error_descript

   ! Which format to use for gnuplot figures:
   read(FN,*,IOSTAT=Reason) numpar%gnupl%gnu_extension
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Does user require to printout DOS of materials?
   read(FN,*,IOSTAT=Reason) numpar%printout_DOS
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Does user require to printout MFPs in materials?
   read(FN,*,IOSTAT=Reason) numpar%printout_MFPs
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Does user require to printout Ranges in materials?
   read(FN,*,IOSTAT=Reason) numpar%printout_ranges
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
end subroutine read_mandatory_output



subroutine read_numerics_definition(FN, File_name, numpar, count_lines, read_well, Err)    ! below
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout) :: numpar	! all numerical parametersn
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   integer :: Reason, temp
   character(200) :: temp_ch, Error_descript

   ! number of MC iterations
   read(FN,*,IOSTAT=Reason) numpar%NMC
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! number of threads for parallel calculations with OpenMP (1 if nonparrelelized)
   read(FN,*,IOSTAT=Reason) numpar%NOMP
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [fs] Time-step for MD
   read(FN,*,IOSTAT=Reason) temp_ch  ! check whether there is time grid provided
   !read(FN,*,IOSTAT=Reason) numpar%dt_MD
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   else
      call set_MD_step_grid(temp_ch, numpar, read_well, Error_descript)  ! module "Read_numerical_parameters"
      if (.not. read_well) then
         call Save_error_details(Err, 2, Error_descript)    ! module "Objects"
         return
      endif
   endif

   ! [fs] How often to print out the data into files
   read(FN,*,IOSTAT=Reason) temp_ch  ! check whether there is time grid provided
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      write(Error_descript,'(a)') trim(adjustl(Error_descript))//'. Setting printout timestep equal to simulation step.'
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      numpar%dt_printout = numpar%dt_MD ! printout data every timestep by default
   else
      call set_time_grid(temp_ch, numpar)   ! module "Read_numerical_parameters"
   endif

   ! [fs] when to start simulation
   read(FN,*,IOSTAT=Reason) numpar%t_start
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [fs] when to stop simulation
   read(FN,*,IOSTAT=Reason) numpar%t_total
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Activate MC module or not:
   read(FN,*,IOSTAT=Reason) temp
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
   if (temp /= 0) then	! perform MC calculations:
      numpar%DO_MC = .true.
   else	! do NOT perform MC calculations:
      numpar%DO_MC = .false.
   endif

   ! Activate MD module or not:
   read(FN,*,IOSTAT=Reason) temp
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
   if (temp /= 0) then	! perform MD calculations:
      numpar%DO_MD = .true.
   else	! do NOT perform MD calculations:
      numpar%DO_MD = .false.
   endif

   ! Activate TTM module or not:
   read(FN,*,IOSTAT=Reason) temp
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
   if (temp /= 0) then	! perform MD calculations:
      numpar%DO_TTM = .true.
   else	! do NOT perform MD calculations:
      numpar%DO_TTM = .false.
   endif

   ! Recalculate cross sections and MFPs (T=true; F=false):
   read(FN,*,IOSTAT=Reason) numpar%recalculate_MFPs
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [A] coordinates of the left and right ends of the simulation window along X, and whether reset it via MD:
   read(FN,*,IOSTAT=Reason) numpar%box_start_x, numpar%box_end_x, numpar%reset_from_MD(1)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [A] coordinates of the left and right ends of the simulation window along Y, and whether reset it via MD:
   read(FN,*,IOSTAT=Reason) numpar%box_start_y, numpar%box_end_y, numpar%reset_from_MD(2)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! [A] coordinates of the left and right ends of the simulation window along Z, and whether reset it via MD:
   read(FN,*,IOSTAT=Reason) numpar%box_start_z, numpar%box_end_z, numpar%reset_from_MD(3)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Bondary conditions along x,y,z: 0=free, 1=periodic
   read(FN,*,IOSTAT=Reason) numpar%periodic_x, numpar%periodic_y, numpar%periodic_z
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

end subroutine read_numerics_definition



subroutine read_radiation_definition(FN, File_name, bunch, count_lines, read_well, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Radiation_param), dimension(:), allocatable, intent(inout) :: bunch   ! incomming radiation
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   integer :: N_rp, i_rp, Reason
   character(200) :: Error_descript

   ! Number of different types of incomming particles:
   read(FN,*,IOSTAT=Reason) N_rp	! How many incoming particles/pulses to model
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! Now we know how many incomming particles we have, allocate arrays:
   if (allocated(bunch)) then
      if (size(bunch) /= N_rp) deallocate(bunch)      ! if default doesn't coincide with defined
   endif
   if (.not. allocated(bunch)) allocate(bunch(N_rp))

   ! Read the parameters of each of the incoming particles / pulses:
   PULS:do i_rp = 1, N_rp	! for each of them
      ! Type of radiation: 0=single particles, 1=pulse / bunch
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%NOP
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif
      if (bunch(i_rp)%NOP < 1) bunch(i_rp)%NOP = 1  ! by default use 1 particle, not less

      ! Kind of particle: 0=photon, 1=electron, 2=positron, 3=SHI, 4=hole
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%KOP
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif

      ! [A] coordinates of impact:
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%R(:)
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif

      ! [A] spread (or uncertainty) of coordinates of impact
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%R_spread(:)
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif

!       ! Spatial shape: 0 = rectangular, 1 = Gaussian, 2 = SASE
!       read(FN,*,IOSTAT=Reason) bunch(i_rp)%Space_shape
!       call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
!       if (.not. read_well) then
!          write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
!          call Save_error_details(Err, 2, Error_descript)	! module "Objects"
!          return
!       endif

      ! [degrees] theta and phi angles of impact:
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%theta, bunch(i_rp)%phi
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif
      ! Convert from degrees to radians:
      bunch(i_rp)%theta = bunch(i_rp)%theta * g_deg2rad
      bunch(i_rp)%phi = bunch(i_rp)%phi * g_deg2rad

      ! [eV] energy of the incomming particle / pulse:
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%E
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif

      ! [eV] spread of energies (or energy uncertainty):
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%E_spread
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif

!       ! Spectral shape: 0 = rectangular, 1 = Gaussian, 2 = SASE:
!       read(FN,*,IOSTAT=Reason) bunch(i_rp)%Spectr_shape
!       call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
!       if (.not. read_well) then
!          write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
!          call Save_error_details(Err, 2, Error_descript)	! module "Objects"
!          return
!       endif

      ! [fs] arrival time of the incomming particle / center of the pulse:
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%t
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif

      ! [fs] FWHM-duration of the pulse (ignorred for single particles):
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%FWHM
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif

!       ! Temporal shape: 0 = rectangular, 1 = Gaussian, 2 = SASE:
!       read(FN,*,IOSTAT=Reason) bunch(i_rp)%Time_shape
!       call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
!       if (.not. read_well) then
!          write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
!          call Save_error_details(Err, 2, Error_descript)	! module "Objects"
!          return
!       endif

      ! In case of SHI, there are additional lines specifying
      if (bunch(i_rp)%KOP == 3) then
         ! Z, atomic number of ion in periodic table:
         read(FN,*,IOSTAT=Reason) bunch(i_rp)%Z
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            call Save_error_details(Err, 2, Error_descript)	! module "Objects"
            return
         endif

         ! User-provided fixed charge [e], and SHI mass [a.m.u]:
         read(FN,*,IOSTAT=Reason) bunch(i_rp)%Zeff,  bunch(i_rp)%Meff
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            call Save_error_details(Err, 2, Error_descript)	! module "Objects"
            return
         endif
      endif

   enddo PULS

end subroutine read_radiation_definition



subroutine read_target_definition(FN, File_name, used_target, count_lines, read_well, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Matter), intent(inout) :: used_target	! parameters of the target
   integer, intent(inout) :: count_lines  ! number of the read line
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------------
   integer :: Reason, i
   character(200) :: Error_descript

   ! Target name:
   read(FN,*,IOSTAT=Reason) used_target%Name
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif

   ! how many target components there are:
   read(FN,*,IOSTAT=Reason) used_target%NOC
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      return
   endif
   ! Now we know how many different types of the targets we have, allocate arrays for them:
   If (allocated(used_target%Material)) then
      if (size(used_target%Material) /= used_target%NOC) deallocate(used_target%Material) ! if default doesn't coincide with defined
   endif
   If (allocated(used_target%Mat_types)) then
      if (size(used_target%Mat_types) /= used_target%NOC) deallocate(used_target%Mat_types) ! if default doesn't coincide with defined
   endif
   If (allocated(used_target%Geom)) then
      if (size(used_target%Geom) /= used_target%NOC) deallocate(used_target%Geom) ! if default doesn't coincide with defined
   endif
   if (.not. allocated(used_target%Material)) allocate(used_target%Material(used_target%NOC))
   if (.not. allocated(used_target%Mat_types)) allocate(used_target%Mat_types(used_target%NOC))
   if (.not. allocated(used_target%Geom)) allocate(used_target%Geom(used_target%NOC))

   ! Read parameters for all the types of the target, all components one by one:
   TRGT:do i = 1, used_target%NOC
      ! Read the chemical formula of this target component:
      read(FN,*,IOSTAT=Reason) used_target%Material(i)%Name
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif

      ! Read the temperature of this target component:
      read(FN,*,IOSTAT=Reason) used_target%Material(i)%T		! [K]
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif
      used_target%Material(i)%T_eV = used_target%Material(i)%T*g_kb_EV	! convert to [eV]

      ! Read the type of the target:
      read(FN,*,IOSTAT=Reason) used_target%Mat_types(i)
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif

      ! Once we know the index, lets set the type for the geometry:
      select case (used_target%Mat_types(i))
      case default	! Rectangle:
         allocate(Rectangle::used_target%Geom(i)%Set)
      case (1)	! Sphere
         allocate(Sphere::used_target%Geom(i)%Set)
      case (2)	! Sphere_segment
         allocate(Sphere_segment::used_target%Geom(i)%Set)
      case (3)	! Cylinder
         allocate(Cylinder::used_target%Geom(i)%Set)
      case (4)	! Cylinder_segment
         allocate(Cylinder_segment::used_target%Geom(i)%Set)
      endselect

      ! Now we can read the parameters of this given type:
      ASSOCIATE (ARRAY => used_target%Geom(i)%Set)	! that's the syntax to use when passing polimorphic arrays into subroutines
       select type (ARRAY)	! Depending on the type of used_target%Geom(i)%Set
         type is (Rectangle)
            call read_param_Rectangle(FN, File_name, Reason, count_lines, array, Err)	! see below
         type is (Sphere)
            call read_param_Sphere(FN, File_name, Reason, count_lines, array, Err)	! see below
         type is (Sphere_segment)
            call read_param_Sphere_segment(FN, File_name, Reason, count_lines, array, Err)	! see below
         type is (Cylinder)
            call read_param_Cylinder(FN, File_name, Reason, count_lines, array, Err)		! see below
         type is (Cylinder_segment)
            call read_param_Cylinder_segment(FN, File_name, Reason, count_lines, array, Err)	! see below
       endselect
      END ASSOCIATE

      ! Read type of potential barrier:
      read(FN,*,IOSTAT=Reason) used_target%Material(i)%Surface_barrier%barr_type	! Type of emission barrier: 0=step, 1=Eckart
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         return
      endif
   enddo TRGT

end subroutine read_target_definition




subroutine set_defaults(used_target, bunch, numpar)
   type(Matter), intent(inout) :: used_target   ! parameters of the target
   type(Radiation_param), dimension(:), allocatable, intent(inout) :: bunch	! incomming radiation
   type(Num_par), intent(inout) :: numpar       ! all numerical parameters
   !---------------------------
   ! Set default parameters:

   ! Target:
   used_target%Name = 'Test_target'
   used_target%NOC = 1
   allocate(used_target%Material(used_target%NOC))
   allocate(used_target%Mat_types(used_target%NOC))
   allocate(used_target%Geom(used_target%NOC))
   used_target%Material(1)%Name = 'Null'    ! No target by default
   used_target%Material(1)%T = 300.0d0      ! Temperature [K]
   used_target%Material(1)%T_eV = used_target%Material(1)%T*g_kb_EV ! [eV]
   used_target%Mat_types(1) = 0     ! rectangle
   used_target%Material(1)%Surface_barrier%barr_type = 1    ! Eckart barrier

   ! Radiation:
   allocate(bunch(1))
   bunch(1)%NOP = 1  ! by default use 1 particle
   bunch(1)%KOP = 0  ! Kind of particle: 0=photon, 1=electron, 2=positron, 3=SHI, 4=hole
   bunch(1)%R(:) = 0.0d0   ! [A] coordinates of impact
   bunch(1)%R_spread(:) = 0.0d0 ! [A] spread (or uncertainty) of coordinates of impact
   bunch(1)%theta = 0.0d0 ! [degrees] theta and phi angles of impact:
   bunch(1)%phi = 0.0d0
   bunch(1)%E = 100.0d0   ! [eV] energy of the incomming particle / pulse
   bunch(1)%E_spread = 0.0d0 ! [eV] spread of energies (or energy uncertainty)
   bunch(1)%t = 0.0d0  ! [fs] arrival time of the incomming particle / center of the pulse
   bunch(1)%FWHM = 0.0d0 ! [fs] FWHM-duration of the pulse (ignorred for single particles)

   ! NUMERICAL PARAMETERS:
   numpar%NMC = 1       ! number of MC iterations
   numpar%NOMP = 1      ! number of threads for parallel calculations with OpenMP (1 if nonparrelelized)
   numpar%dt_MD = 1.0d0 ! time step
   numpar%i_dt = -1     ! to mark that the reset option is unused
   numpar%dt_printout = numpar%dt_MD ! printout data every timestep by default
   numpar%t_start = 0.0d0     ! [fs] when to start simulation
   numpar%t_total = 100.0d0   ! [fs] when to stop simulation
   numpar%DO_MC = .false.     ! do NOT perform MC calculations
   numpar%DO_MD = .false.     ! do NOT perform MD calculations
   numpar%DO_TTM = .false.    ! do NOT perform TTN calculations
   numpar%recalculate_MFPs = .false. ! do NOT recalculate MFPs
   ! [A] coordinates of the left and right ends of the simulation box along X, and whether reset it via MD:
   numpar%box_start_x = -10.0d10
   numpar%box_end_x = 10.0d10
   numpar%reset_from_MD(1) = .false.
   ! [A] coordinates of the left and right ends of the simulation box along Y, and whether reset it via MD:
   numpar%box_start_y = -10.0d10
   numpar%box_end_y  = 10.0d10
   numpar%reset_from_MD(2) = .false.
   ! [A] coordinates of the left and right ends of the simulation box along Z, and whether reset it via MD:
   numpar%box_start_z = 0.0d0
   numpar%box_end_z = 10.0d0
   numpar%reset_from_MD(3) = .false.
   ! Kind of boundary along X, Y, Z (0=absorbing; 1=periodic; 2=reflective; 3=white):
   numpar%periodic_x = 1
   numpar%periodic_y = 1
   numpar%periodic_z = 1

   !::: MODELS FOR ELECTRONS :::
   numpar%MC_vs_MD = 0  ! MC or MD target model: 0=MC, 1=MD
   numpar%El_forces = 0 ! Include forces and fields among electrons: 0=exclude, 1=include
   numpar%El_inelast = 3 ! inelastic scattering: 0=excluded, 1=numerical relativ.CDF, 2=RBEB, 3=delta, 4=nonrelativ.CDF (DO NOT USE!), 5=SPdelta
   numpar%El_elastic = 1 ! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   numpar%El_Brems = 1        ! Bremsstrahlung: 0=excluded, 1=BHW.
   numpar%El_Cheren = 0       ! Cherenkov radiation: 0=excluded, 1= NOT READY
   numpar%El_Cutoff = 0.1d0 ! [eV] Cut-off energy (electrons with lower energies are excluded from calculation)

   !::: MODELS FOR PHOTONS :::
   numpar%Ph_absorb = 2  ! Photoabsorption CSs: 0=excluded, 1=CDF, 2=EPDL
   numpar%Ph_Compton = 1 ! Compton effect: 0=excluded, 1=PENELOPE
   numpar%Ph_Thomson = 1 ! Thomson / Rayleigh scattering: 0=excluded, 1=PENELOPE
   numpar%Ph_Pairs = 1   ! Electron-positron pair creation: 0=excluded, 1=included, 2=included with Landau-Pomeranchuk-Migdal (not ready)
   numpar%Ph_Nucl = 0  ! ! Photonuclear Physics: 0=excluded, 1=included (not ready)
   numpar%Ph_Cutoff = -1.0d0    ! [eV] Cut-off energy (photons with lower energies are excluded from calculation)
   numpar%Ph_att_eff = -1.0d0   ! [A] Effective photon attenuation length (<0 means do not use effective one, use real one from EPDL)

   !::: MODELS FOR SHI :::
   numpar%SHI_inelast = 1     ! SHI inelastic scattering: 0=excluded, 1:3=delta, 4=nonrelativ.CDF (DO NOT USE!), 5=SPdelta
   numpar%SHI_ch_st = 0       ! Charge state: 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande, 4=fixed Zeff, 5=charge exchange
   numpar%SHI_ch_shape = 0    ! SHI charge shape: 0=point-like charge; 1=Brandt-Kitagawa ion:
   numpar%SHI_Cutoff = -1.0d0 ! [eV] Cut-off energy (SHIs with lower energies are excluded from calculation)

   !::: MODELS FOR POSITRONS :::
   numpar%Pos_inelast = 1 ! Positron inelastic scattering: 0=excluded, 1=delta
   numpar%Pos_elastic = 1 ! Positron elastic scattering: 0=excluded, 1=delta
   numpar%Pos_Brems = 1 ! Positron Bremsstrahlung scattering: 0=excluded, 1=delta
   numpar%Pos_annih = 1 ! Positron annihilation: 0=excluded, 1= Heitler
   numpar%Pos_Cutoff = 0.1d0  ! [eV] Cut-off energy (positrons with lower energies are excluded from calculation)

   !::: MODEL PARAMETERS FOR CDF :::
   numpar%CDF_model = 0 ! CDF model: 0 = Drude / Lindhard CDF, 1=Mermin CDF, 2=Full conserving CDF (not ready)
   ! target dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie; effective mass [in me] (0=effective mass from DOS of VB; -1=free-electron):
   numpar%CDF_dispers = 0
   numpar%CDF_m_eff = 0
   numpar%CDF_plasmon = 0     ! Include plasmon integration limit (0=no, 1=yes)
   numpar%CDF_Eeq_factor = 10.0d0 ! Coefficient where to switch from Ritchie to Delta CDF: E = k * Wmin (INELASTIC)
   numpar%CDF_Eeq_elast= 10.0d0   ! coeff.k where to switch from nonrelativistic to Delta CDF for ELASTIC scattering: E = k * Wmin
   numpar%CDF_elast_Zeff = 0      ! use for target atoms Zeff (set 0), or Z=1 (set 1)
   numpar%CDF_int_n_inelastE = 50  ! n grid points for INELASTIC cross section integration over energy (E): dE = max((E - E0(:)), G(:))/n
   numpar%CDF_int_n_inelastQ = 100 ! n grid points for INELASTIC cross section integration over momentum (Q): dQ = max((Q - (W-E0(:))), G(:))/n
   numpar%CDF_int_n_elastE = 10  ! n grid points for ELASTIC cross section integration over energy (E): dE = max((E - E0(:)), G(:))/n
   numpar%CDF_int_n_elastQ = 100 ! n grid points for ELASTIC cross section integration over momentum (Q): dQ = max((Q - (W-E0(:))), G(:))/n

   !::: MODELS FOR HOLES :::
   numpar%H_Auger = 1   ! Auger decays:  0=excluded, 1=EADL
   numpar%H_Radiat = 1  ! Radiative decays: 0=excluded, 1=EADL
   numpar%H_m_eff = -1  ! [me] effective valence hole mass (-1=from DOS; 0=free electron; >0=fixed mass in m_e)
   numpar%H_inelast = 1 ! Valence hole inelastic scattering: 0=excluded, 1=CDF, 2=BEB, 3=Delta
   numpar%H_elast = 1   ! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF (NOT READY YET!), 5=SP-CDF
   numpar%H_Cutoff = 0.1d0       ! [eV] Cut-off energy (holes with lower energies are excluded from calculation)

   ! ::: MD MODEL PARAMETERS :::
   ! Use quenching in MD (T=true; F=false), when to start [fs], how often nullify velocities [fs]:
   numpar%do_quenching = .false.
   numpar%t_quench_start = 0.0d0
   numpar%dt_quench = 1.0d0
   numpar%t_quench_run = 0.0d0

   ! ::: OUTPUT DATA :::
   numpar%gnupl%gnu_extension = 'png'
   numpar%printout_DOS = .true.     ! printout DOS
   numpar%printout_MFPs = .true.    ! printout MFPs
   numpar%printout_ranges = .true.  ! printout Ranges
   ! Optional:
   numpar%MD_force_ind = 0 ! do NOT calculate forces as numerical derivative of the potential
   numpar%do_cohesive = .false.    ! calculate cohesive energy (instead of full TREKIS run)
   numpar%print_MD_R_xyz = .false.    ! printout MD atomic coordinates in XYZ
   numpar%print_MD_V_xyz = .false.    ! printout MD atomic velosities in XYZ
   numpar%print_MD_LAMMPS = .false.   ! create input file for LAMMPS at the final time instant
   numpar%print_MC_MD_energy = .false.    ! printout MC-MD energy transfer
   numpar%vel_theta_grid_par%along_axis = .false.     ! printout particles theta-distribution
   numpar%print_each_step = .false.    ! printout each MD step timing
   numpar%NRG_grid_par%along_axis = .false.     ! spectra printout
   numpar%Spectr_grid_par(:)%along_axis = .false.     ! printout spectra along various axes
   numpar%grid_par(:)%along_axis = .false.
end subroutine set_defaults



subroutine read_new_format_input(FN, File_name, used_target, bunch, Err)      ! UNFINISHED
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Matter), intent(inout) :: used_target	! parameters of the target
   type(Radiation_param), dimension(:), allocatable, intent(inout) :: bunch	! incomming radiation
   type(Error_handling), intent(inout) :: Err	! error log
   !----------------------
   integer :: N_targets, N_bunches
   character(200) :: Error_descript

   ! Read global target parameters:
   call read_global_target(FN, used_target) ! below

end subroutine read_new_format_input


subroutine read_global_target(FN, used_target) ! UNFINISHED
   integer, intent(in) :: FN    ! file number to read from
   type(Matter), intent(inout) :: used_target	! parameters of the target
   !----------------------
   integer :: n_targets
   character(100) :: target_name
   namelist / global_target / n_targets, target_name
   read(FN, nml = global_target)
   ! Save the read variables into the global ones:
   used_target%NOC = n_targets  ! save number of targets
   used_target%Name = target_name
end subroutine read_global_target


!----------------------------
! If two input files are used (old format, all data given):
subroutine read_input_parameters(FN, File_name, used_target, bunch, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Matter), intent(inout) :: used_target	! parameters of the target
   type(Radiation_param), dimension(:), allocatable, intent(inout) :: bunch	! incomming radiation
   type(Error_handling), intent(inout) :: Err	! error log
   !----------------------
   integer :: Reason, count_lines, i, N_rp, i_rp
   character(200) :: Error_descript
   logical :: read_well
   
   count_lines = 0	! to start counting lines in the file

   ! Skip the first line: it is only an indicator. May be used for comments.
   read(FN,*,IOSTAT=Reason)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Target name:
   read(FN,*,IOSTAT=Reason) used_target%Name
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! how many target components there are:
   read(FN,*,IOSTAT=Reason) used_target%NOC
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   ! Now we know how many different types of the targets we have, allocate arrays for them:
   If (allocated(used_target%Material)) then
      if (size(used_target%Material) /= used_target%NOC) deallocate(used_target%Material) ! if default doesn't coincide with defined
   endif
   If (allocated(used_target%Mat_types)) then
      if (size(used_target%Mat_types) /= used_target%NOC) deallocate(used_target%Mat_types) ! if default doesn't coincide with defined
   endif
   If (allocated(used_target%Geom)) then
      if (size(used_target%Geom) /= used_target%NOC) deallocate(used_target%Geom) ! if default doesn't coincide with defined
   endif
   if (.not. allocated(used_target%Material)) allocate(used_target%Material(used_target%NOC))
   if (.not. allocated(used_target%Mat_types)) allocate(used_target%Mat_types(used_target%NOC))
   if (.not. allocated(used_target%Geom)) allocate(used_target%Geom(used_target%NOC))

   ! Read parameters for all the types of the target, all components one by one:
   TRGT:do i = 1, used_target%NOC
      ! Read the chemical formula of this target component:
      read(FN,*,IOSTAT=Reason) used_target%Material(i)%Name
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
     
      ! Read the temperature of this target component:
      read(FN,*,IOSTAT=Reason) used_target%Material(i)%T		! [K]
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      used_target%Material(i)%T_eV = used_target%Material(i)%T*g_kb_EV	! convert to [eV]
      
      ! Read the type of the target:
      read(FN,*,IOSTAT=Reason) used_target%Mat_types(i)
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      
      ! Once we know the index, lets set the type for the geometry:
      select case (used_target%Mat_types(i))
      case default	! Rectangle:
         allocate(Rectangle::used_target%Geom(i)%Set)
      case (1)	! Sphere
         allocate(Sphere::used_target%Geom(i)%Set)
      case (2)	! Sphere_segment
         allocate(Sphere_segment::used_target%Geom(i)%Set) 
      case (3)	! Cylinder
         allocate(Cylinder::used_target%Geom(i)%Set)
      case (4)	! Cylinder_segment
         allocate(Cylinder_segment::used_target%Geom(i)%Set)
      endselect
      
      ! Now we can read the parameters of this given type:
      ASSOCIATE (ARRAY => used_target%Geom(i)%Set)	! that's the syntax to use when passing polimorphic arrays into subroutines
       select type (ARRAY)	! Depending on the type of used_target%Geom(i)%Set
         type is (Rectangle)
            call read_param_Rectangle(FN, File_name, Reason, count_lines, array, Err)	! see below
         type is (Sphere)
            call read_param_Sphere(FN, File_name, Reason, count_lines, array, Err)	! see below
         type is (Sphere_segment)
            call read_param_Sphere_segment(FN, File_name, Reason, count_lines, array, Err)	! see below
         type is (Cylinder)
            call read_param_Cylinder(FN, File_name, Reason, count_lines, array, Err)		! see below
         type is (Cylinder_segment)
            call read_param_Cylinder_segment(FN, File_name, Reason, count_lines, array, Err)	! see below
       endselect
      END ASSOCIATE
      
      ! Read type of potential barrier:
      read(FN,*,IOSTAT=Reason) used_target%Material(i)%Surface_barrier%barr_type	! Type of emission barrier: 0=step, 1=Eckart
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      
   enddo TRGT
   
   !00000000000000000000000000
   ! Now, define the radiation type:
   ! Skip line
   read(FN,*,IOSTAT=Reason) !    ::: RADIATION :::
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Number of different types of incomming particles:
   read(FN,*,IOSTAT=Reason) N_rp	! How many incoming particles/pulses to model
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Now we know how many incomming particles we have, allocate arrays:
   if (allocated(bunch)) then
      if (size(bunch) /= N_rp) deallocate(bunch)      ! if default doesn't coincide with defined
   endif
   if (.not. allocated(bunch)) allocate(bunch(N_rp))
   
   ! Read the parameters of each of the incoming particles / pulses:
   PULS:do i_rp = 1, N_rp	! for each of them
      ! Type of radiation: 0=single particles, 1=pulse / bunch
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%NOP
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      if (bunch(i_rp)%NOP < 1) bunch(i_rp)%NOP = 1  ! by default use 1 particle, not less
      
      ! Kind of particle: 0=photon, 1=electron, 2=positron, 3=SHI, 4=hole
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%KOP
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      
      ! [A] coordinates of impact:
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%R(:)
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
   
      ! [A] spread (or uncertainty) of coordinates of impact
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%R_spread(:)
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      
!       ! Spatial shape: 0 = rectangular, 1 = Gaussian, 2 = SASE
!       read(FN,*,IOSTAT=Reason) bunch(i_rp)%Space_shape
!       call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
!       if (.not. read_well) then
!          write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
!          call Save_error_details(Err, 2, Error_descript)	! module "Objects"
!          goto 9998
!       endif
      
      ! [degrees] theta and phi angles of impact:
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%theta, bunch(i_rp)%phi
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      ! Convert from degrees to radians:
      bunch(i_rp)%theta = bunch(i_rp)%theta * g_deg2rad
      bunch(i_rp)%phi = bunch(i_rp)%phi * g_deg2rad
      
      ! [eV] energy of the incomming particle / pulse:
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%E
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      
      ! [eV] spread of energies (or energy uncertainty):
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%E_spread
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      
!       ! Spectral shape: 0 = rectangular, 1 = Gaussian, 2 = SASE:
!       read(FN,*,IOSTAT=Reason) bunch(i_rp)%Spectr_shape
!       call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
!       if (.not. read_well) then
!          write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
!          call Save_error_details(Err, 2, Error_descript)	! module "Objects"
!          goto 9998
!       endif
      
      ! [fs] arrival time of the incomming particle / center of the pulse:
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%t
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      
      ! [fs] FWHM-duration of the pulse (ignorred for single particles):
      read(FN,*,IOSTAT=Reason) bunch(i_rp)%FWHM
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9998
      endif
      
!       ! Temporal shape: 0 = rectangular, 1 = Gaussian, 2 = SASE:
!       read(FN,*,IOSTAT=Reason) bunch(i_rp)%Time_shape
!       call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
!       if (.not. read_well) then
!          write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
!          call Save_error_details(Err, 2, Error_descript)	! module "Objects"
!          goto 9998
!       endif
      
      ! In case of SHI, there are additional lines specifying
      if (bunch(i_rp)%KOP == 3) then
         ! Z, atomic number of ion in periodic table:
         read(FN,*,IOSTAT=Reason) bunch(i_rp)%Z
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            call Save_error_details(Err, 2, Error_descript)	! module "Objects"
            goto 9998
         endif
         
         ! User-provided fixed charge [e], and SHI mass [a.m.u]:
         read(FN,*,IOSTAT=Reason) bunch(i_rp)%Zeff,  bunch(i_rp)%Meff
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            call Save_error_details(Err, 2, Error_descript)	! module "Objects"
            goto 9998
         endif
      endif
      
   enddo PULS
9998 continue
end subroutine read_input_parameters


!==============================
! Read different geometries from input files:
subroutine read_param_Rectangle(FN, File_name, Reason, count_lines, array, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   integer, intent(inout) :: Reason	! index for reading file
   integer, intent(inout) :: count_lines	! which line of the file we are reading now
   type(Rectangle), intent(inout) :: array	! the target geometry array for the case of rectangle
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------
   character(200) :: Error_descript
   logical :: read_well
   
   read(FN,*,IOSTAT=Reason) array%X, array%Y, array%Z	! coordinates of its center
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9997
   endif
   
   read(FN,*,IOSTAT=Reason) array%Xstart, array%Xend	! Length: its beginning and its end along X axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9997
   endif
   
   read(FN,*,IOSTAT=Reason) array%Ystart, array%Yend	! Length: its beginning and its end along Y axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9997
   endif
   
   read(FN,*,IOSTAT=Reason) array%Zstart, array%Zend	! Length: its beginning and its end along Z axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9997
   endif
   
   read(FN,*,IOSTAT=Reason) array%angle_x	! its rotation around X-axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9997
   endif
   array%angle_x = array%angle_x * g_deg2rad    ! [deg] -> [rad]
   
   read(FN,*,IOSTAT=Reason) array%angle_y	! its rotation around Y-axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9997
   endif
   array%angle_y = array%angle_y * g_deg2rad    ! [deg] -> [rad]
   
   read(FN,*,IOSTAT=Reason) array%angle_z	! its rotation around Z-axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9997
   endif
   array%angle_z = array%angle_z * g_deg2rad    ! [deg] -> [rad]

9997 continue
end subroutine read_param_Rectangle



subroutine read_param_Sphere(FN, File_name, Reason, count_lines, array, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   integer, intent(inout) :: Reason	! index for reading file
   integer, intent(inout) :: count_lines	! which line of the file we are reading now
   type(Sphere), intent(inout) :: array	! the target geometry array for the case of sphere
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------
   character(200) :: Error_descript
   logical :: read_well
   
   read(FN,*,IOSTAT=Reason) array%X, array%Y, array%Z	! coordinates of its center
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9996
   endif
   
   read(FN,*,IOSTAT=Reason) array%R	! its radius
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9996
   endif
   
9996 continue
end subroutine read_param_Sphere


subroutine read_param_Sphere_segment(FN, File_name, Reason, count_lines, array, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   integer, intent(inout) :: Reason	! index for reading file
   integer, intent(inout) :: count_lines	! which line of the file we are reading now
   type(Sphere_segment), intent(inout) :: array	! the target geometry array for the case of Sphere_segment
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------
   character(200) :: Error_descript
   logical :: read_well
   
   read(FN,*,IOSTAT=Reason) array%X, array%Y, array%Z	! coordinates of its center
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9995
   endif
   
   read(FN,*,IOSTAT=Reason) array%R_start, array%R_end		! starting and ending points by radius (if it's cut as a spherical layer)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9995
   endif
   
   read(FN,*,IOSTAT=Reason) array%phi_start, array%phi_end	! starting and ending points by phi angle within [-Pi/2..Pi/2] (if it's cut as excluded cone)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9995
   endif
   
   read(FN,*,IOSTAT=Reason) array%theta_start, array%theta_end	! starting and ending points by theta angle within [0..2*Pi) (if it's cut as a segment)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9995
   endif
   
9995 continue
end subroutine read_param_Sphere_segment


subroutine read_param_Cylinder(FN, File_name, Reason, count_lines, array, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   integer, intent(inout) :: Reason	! index for reading file
   integer, intent(inout) :: count_lines	! which line of the file we are reading now
   type(Cylinder), intent(inout) :: array	! the target geometry array for the case of Cylinder
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------
   character(200) :: Error_descript
   logical :: read_well
   
   read(FN,*,IOSTAT=Reason) array%X, array%Y, array%Z	! coordinates of its center
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9994
   endif
   
   read(FN,*,IOSTAT=Reason) array%R	! radius
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9994
   endif
   
   read(FN,*,IOSTAT=Reason) array%L_start, array%L_end	! starting and ending points of the cylinder defining its length
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9994
   endif
   
   read(FN,*,IOSTAT=Reason) array%angle_x	! its rotation around X-axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9994
   endif
   array%angle_x = array%angle_x * g_deg2rad    ! [deg] -> [rad]
   
   read(FN,*,IOSTAT=Reason) array%angle_y	! its rotation around Y-axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9994
   endif
   array%angle_y = array%angle_y * g_deg2rad    ! [deg] -> [rad]
   
   read(FN,*,IOSTAT=Reason) array%angle_z	! its rotation around Z-axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9994
   endif
   array%angle_z = array%angle_z * g_deg2rad    ! [deg] -> [rad]
   
9994 continue
end subroutine read_param_Cylinder


subroutine read_param_Cylinder_segment(FN, File_name, Reason, count_lines, array, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   integer, intent(inout) :: Reason	! index for reading file
   integer, intent(inout) :: count_lines	! which line of the file we are reading now
   type(Cylinder_segment), intent(inout) :: array	! the target geometry array for the case of Cylinder_segment
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------
   character(200) :: Error_descript
   logical :: read_well
   
   read(FN,*,IOSTAT=Reason) array%X, array%Y, array%Z	! coordinates of its center
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9993
   endif
   
   read(FN,*,IOSTAT=Reason) array%R_start, array%R_end	! starting and ending radius (if it's cut as cylindrical layer)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9993
   endif
   
   read(FN,*,IOSTAT=Reason) array%L_start, array%L_end	! starting and ending points of the cylinder defining its length
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9993
   endif
   
   read(FN,*,IOSTAT=Reason) array%theta_start, array%theta_end	! starting and ending points by theta angle within [0..2*Pi) (if it's cut as a segment)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9993
   endif
   
   read(FN,*,IOSTAT=Reason) array%angle_x	! its rotation around X-axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9993
   endif
   array%angle_x = array%angle_x * g_deg2rad    ! [deg] -> [rad]
   
   read(FN,*,IOSTAT=Reason) array%angle_y	! its rotation around Y-axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9993
   endif
   array%angle_y = array%angle_y * g_deg2rad    ! [deg] -> [rad]
   
   read(FN,*,IOSTAT=Reason) array%angle_z	! its rotation around Z-axis
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9993
   endif
   array%angle_z = array%angle_z * g_deg2rad    ! [deg] -> [rad]
   
9993 continue
end subroutine read_param_Cylinder_segment


END MODULE Read_input_data
