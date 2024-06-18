!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TREKIS: Time-REsolved Kinetics in Irradiated Solids
! The code is written in 2018 - 2021 by
! N. Medvedev, R. Rymzhanov, F. Akhmetov, R. Voronkov...
! 
!
! Should you have any questions, contact the authors: 
! (1) nikita.medvedev@fzu.cz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DISCLAMER:
! Although we endeavour to ensure that the code TREKIS and results delivered are correct, no warranty is given as to its accuracy. 
! We assume no responsibility for possible errors or omissions. We shall not be liable for any damage arising from the use of this 
! code or its parts or any results produced with it, or from any action or decision taken as a result of using this code or any related material.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONVENTIONS OF PROGRAMMING:
! 1) All global variables start with "g_", e.g. g_numpar, and all defined in the module "Variables"
! 2) All modular variable names are defined starting as "m_", e.g. "m_number"
! 3) All local variables used within subrounies should NOT start with "g_" or "m_"
! 4) Add a comment after each subroutine and function specifying in which module it can be found
! 5) Leave comments describing what EACH LINE of the code is doing
! 6) Each end(smth) statement should be commented to which block it belongs, e.g.: if (i<k) then ... endif ! (i<k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONVENTIONS OF FILE AND MODULE NAMES:
! 1) Modules dealing with different file formats must start with "Dealing_with_[format]", e.g. Dealing_with_cdf
! 2) Modules with cross sections calculations must start with "CS_[particle]_[model]", e.g. CS_electron_CDF, etc.
! 3) Modules related to Monte Carlo subroutines must start with "MC_"
! 4) Modules related to Molecular Dynamics must start with "MD_"
! 5) Modules related to any other model must start with the name or abbreviation of this model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM TREKIS
!MMMMMMMMMMMMMMMMMMMMMMM
! Initiate modules with all the 'use' statements collected in a separate file:
include 'Use_statements.f90'    ! include parts of the code from an external file
! Open_MP related modules from external libraries:
#ifdef _OPENMP
   USE IFLPORT
   USE OMP_LIB
#endif

implicit none

!--------------------------------------------------------------
! Print the title on the screen:
call print_TRKEIS4_lable(6) ! module "Output"
!--------------------------------------------------------------


! The code is executed in multiple steps:
! 0) Define the basic variables used, initiate random seed and get the current time:
call random_seed()  ! standard FORTRAN seeding of random numbers
call date_and_time(values=g_c1)     ! standard FORTRAN time and date
g_ctim=g_c1     ! save the timestamp
g_Err%Err = .false.         ! no errors yet
g_Err%File_Name = 'Error_log.txt'   ! set the name of the error log file
g_Err%Err_descript = ''         ! no errors yet
open(newunit=g_Err%File_Num, FILE = trim(adjustl(g_Err%File_Name)))     ! Open error log file
! Print out the current time:
call print_time('Start at', ind=0)  ! module "Little_subroutines"
!--------------------------------------------------------------
! 1) Read input files:
call Read_input(g_target, g_numpar, g_bunch, g_Err)     ! module "Read_input_data"
if (g_Err%Err) goto 9999    ! if an error occured while reading input files, terminate the program

! In case user provided time-grid, conform the MD-step parameters to it:
call reset_dt(g_numpar, 0.0d0)  ! below

! Set gnuplot parameters:
call process_user_gnu_parameters(g_numpar%gnupl%gnu_extension, g_numpar%gnupl%gnu_terminal, g_numpar%gnupl%do_gnuplot)   ! module "Gnuplotting"

! Set the OpenMP parallelization options:
call set_OMP_number(g_numpar%NOMP, .true., 6)    ! module "Output"

! Set default values:
call Set_defaults(g_numpar, g_bunch, g_MC)  ! module "Initial_conditions"
! Interprete the target names, and get the material parameters accordingly:
call Get_targets_parameters(g_target, g_numpar, g_Err)  ! module "Read_input_data"
if (g_Err%Err) goto 9999    ! if an error occured while reading input files, terminate the program
! Read DOS for the target material:
call read_DOS_files(g_target, g_numpar, g_Err)    ! module "Dealing_with_DOS"
if (g_Err%Err) goto 9999    ! if an error occured while reading input files, terminate the program
! Set the surface barrier parameters:
call set_all_barriers_parameters(g_target)  ! module "MC_general_tools"

! Get MD parameters, if user runs MD simulation, such as
! supercell parameters, atomic coordinates and velocities, MD potentials:
call Read_MD_input(g_target, g_numpar, g_MD_atoms, g_MD_supce, g_MD_pots, g_Err) ! module "Read_input_data"
if (g_Err%Err) goto 9999    ! if an error occured while reading input files, terminate the program

!--------------------------------------------------------------
! 2) Create output directory:
call make_output_folder(g_target, g_numpar, g_bunch)    ! module "Output"

! Save DOS parameters if required:
call printout_DOS(g_target, g_numpar)  ! module "Output"

!--------------------------------------------------------------
! 3) Set initial parameters:
call Set_initial_parameters(g_numpar, g_target, g_bunch, g_MC)  ! module "Initial_conditions"

! 3.a) Find within which target each incident particle is:
call Find_starting_targets(g_target, g_numpar, g_bunch, g_MC)  ! module "MC_general_tools"

!--------------------------------------------------------------
! Printout the title of the program on the screen:
call Print_title(6, g_target, g_numpar, g_bunch, g_MD_atoms, g_MD_supce, g_MD_pots, do_lable=.false.)    ! module "Output"
! and also printout the title into an output file:
open(newunit = g_numpar%FN_parameters, FILE = trim(adjustl(g_numpar%FILE_parameters )))
call Print_title(g_numpar%FN_parameters, g_target, g_numpar, g_bunch, g_MD_atoms, g_MD_supce, g_MD_pots, do_lable=.true.)   ! module "Output"
! Save the atomic parameters from EADL to the Parameters file:
call Print_atomic_parameters(g_numpar%FN_parameters, g_target)  ! module "Output"

! Printout atomic data on the screen, if needed:
! call Print_atomic_parameters(0, g_target)	! module "Output"

! 3a) Prepare arrays for output distributions from MC:
call allocate_output(g_target, g_numpar, g_output)  ! module "MC_data_analysis"

!--------------------------------------------------------------
! 4) Prepare cross sections (and mean free paths):
! a) photons:
! a.1) cross sections of photoabsorption (together with construction of CDF functions):
call get_CDF(g_numpar, g_target, g_Err) ! module "CDF_get_from_data"
! a.2) renormalize alpha parameters in Delta-CDF to match Ritchie-CDF cross-section at the switching point:
! Commented out for Testing (uncomment for release):
call renormalize_alpha_CDF(g_numpar, g_target%Material) ! module "CS_electrons_inelastic"

!--------------------------------------------------------------
! Save CDF parameters into the output file:
call Print_CDF_parameters(g_numpar%FN_parameters, g_target)  ! module "Output"
!--------------------------------------------------------------

! a.2) cross sections of Compton scattering:
call get_photon_Compton(g_target%Material, g_numpar, g_Err)     ! module "CS_photons"
! a.3) cross sections of pair creation:
call get_photon_pair_creation(g_target%Material, g_numpar, g_Err)   ! module "CS_photons"
! a.4) cross sections of coherent (aka elastic, aka Rayleigh, aka Thomson) scattering:
call get_photon_Rayleigh(g_target%Material, g_numpar, g_Err)    ! module "CS_photons"

! b) electrons:
! b.1) Inelastic:
call get_electron_IMFP(g_target%Material, g_numpar, g_Err)  ! module "CS_electrons_inelastic"
! b.2) Elastic:
call get_electron_EMFP(g_target%Material, g_numpar, g_Err)  ! module "CS_electrons_elastic"
! b.3) Bremsstrahlung:
call get_electron_Brems(g_target%Material, g_numpar, g_Err) ! module "CS_electrons_Bremsstrahlung"

! c) SHI:
! c.1) Inelastic scattering:
call get_ion_IMFP(g_target%Material, g_numpar, g_bunch, g_MC, g_Err)	! module "CS_ion_inelastic"
! c.2) Elastic scattering:
! Not ready yet...

! d) Positron:
! d.1) Inelastic:
call get_positron_IMFP(g_target%Material, g_numpar, g_Err)	! module "CS_positrons_inelastic"
! d.2) Elastic:
call get_positron_EMFP(g_target%Material, g_numpar, g_Err)	! module "CS_positrons_elastic"
! d.3) Bremsstrahlung:
call get_positron_Brems(g_target%Material, g_numpar, g_Err)	! module "CS_positrons_Bremsstrahlung"
! d.4) Annihilation:
call get_positron_annihilation(g_target%Material, g_numpar, g_Err)	! module "CS_positrons_annihilation"

! e) Valence hole:
! e.1) Inelastic scattering:
call get_hole_IMFP(g_target%Material, g_numpar, g_Err)  ! module "CS_holes_inelastic"
! e.2) Elastic scattering:
call get_holes_EMFP(g_target%Material, g_numpar, g_Err)  ! module "CS_holes_elastic"

!--------------------------------------------------------------
! Prepare output files with mean free paths, if user requested:
call print_MFPs(g_target, g_numpar, g_bunch)  ! module "Output"

! Create gnuplot scripts:
if (g_numpar%gnupl%do_gnuplot) call create_gnuplot_files(g_target, g_MD_pots, g_numpar)  ! module "Output"

!--------------------------------------------------------------
! Simulation run:
! write(6,'(a)') '******************************************************************'
write(6,'(a)') trim(adjustl(m_starline))    ! module "Output"

! Sample the MC particles first free flight distance from the surface to the first event:
call prepare_MC_run(g_target, g_MC, g_numpar, g_MD_supce, g_MD_supce%E_e_from_MC, g_MD_supce%E_h_from_MC)   ! module "MC"

! Prepare for MD run - potential and forces before the start:
call prepare_MD_run(g_MD_atoms, g_MD_supce, g_MD_pots, g_numpar)    ! module "MD"

! If the user requested to calculated cohesive energy instead of full run:
if (g_numpar%do_cohesive) then
   call calculate_cohesive_energy(g_numpar, g_MD_atoms, g_MD_supce, g_MD_pots)  ! module "MD"
   goto 9999
endif


! Reset starting time of simulation, if needed:
g_numpar%t_start = set_starting_time(g_bunch, g_numpar%t_start)   ! module "Initial_conditions"
g_numpar%i_dt_printout = 1  ! to start with
g_dt_save = 0.0d0    ! just to start
g_time = g_numpar%t_start         ! just to start
! Check the time of printout:
call reset_dt_out(g_numpar, g_dt_out)   ! below
call print_time_step('Simulation time start:', g_time, msec=.true.)   ! module "Little_subroutines"

! print*, 'TEST -1'

! Collect printable output data from raw MC and MD data:
call analyze_MC_output_data(g_target, g_numpar, g_MC, g_output, g_time)    ! module "MC_data_analysis"

! print*, 'TEST -0.5'

call analyze_MD_output_data(g_numpar, g_MD_supce, g_MD_pots, g_MD_atoms, g_output)       ! module "MD_data_analysis"

! print*, 'TEST 0'

! Save the first output data:
call write_output_files(g_target, g_numpar, g_output, g_MD_atoms, g_MD_supce, g_MD_pots, g_time)    ! module "Output"

! print*, 'TEST 1'

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
! Start time propagation:
TP:do while (g_time+g_numpar%dt_MD <= g_numpar%t_total)
   ! In case the user provided time-grid, reset the time-step if required:
   call reset_dt(g_numpar, g_time)  ! below

   ! 1) High-energy electron kinetics (Monte Carlo):
   call MC_step(g_time, g_numpar%dt_MD, g_MC, g_target, g_numpar, g_bunch, g_MD_supce) ! module "MC"

!    print*, 'TEST 2'

   !--------------------------------------------------------------
   ! 2) Low-energy electron kinetics (TTM, Boltzmann, etc):
   ! NOT INCLUDED YET
   
   !--------------------------------------------------------------
   ! 3) Atomic system (Molecular Dynamics):
   call MD_step(g_MD_atoms, g_MD_supce, g_MD_pots, g_numpar, g_time, g_numpar%dt_MD)  ! module "MD"

   !--------------------------------------------------------------
   ! Prepare for the next timestep:
   g_time = g_time + g_numpar%dt_MD     ! [fs] next time-step
   g_dt_save = g_dt_save + g_numpar%dt_MD   ! [fs] for tracing when to save the output data
   !--------------------------------------------------------------

   !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   ! Write current data into output files:
   ! If it's time to printout data (or, always printout the last step):
   if ( (g_dt_save >= g_dt_out) .or. (g_time+g_numpar%dt_MD >= g_numpar%t_total) ) then
      ! Print out the curent time-step
      call print_time_step('Simulation time (out):', g_time, msec=.true.)   ! module "Little_subroutines"
      ! Collect printable output data from raw MC data:
      call analyze_MC_output_data(g_target, g_numpar, g_MC, g_output, g_time)    ! module "MC_data_analysis"
      call analyze_MD_output_data(g_numpar, g_MD_supce, g_MD_pots, g_MD_atoms, g_output)       ! module "MD_data_analysis"

!       call print_time_step('Before printout:', g_time, msec=.true.)   ! module "Little_subroutines"
      ! Save current output data:
      if (g_time+g_numpar%dt_MD > g_numpar%t_total) then   ! last timestep is special
         call write_output_files(g_target, g_numpar, g_output, g_MD_atoms, g_MD_supce, g_MD_pots, g_time, laststep=.true.) ! module "Output"
      else ! other steps are regular
         call write_output_files(g_target, g_numpar, g_output, g_MD_atoms, g_MD_supce, g_MD_pots, g_time)    ! module "Output"
      endif
!       call print_time_step('After printout:', g_time, msec=.true.)   ! module "Little_subroutines"
      ! Allow user to comunicate with the program (program reads your commands from the comunication file):
      call communicate_with_user(g_numpar%FN_communication, g_time, g_numpar) ! module "Output"
      g_dt_save = 0.0d0 ! reset the counter
      ! Check if the time of printout changes:
      call reset_dt_out(g_numpar, g_dt_out)   ! below
   else
      if (g_numpar%print_each_step) call print_time_step('Simulation time step :', g_time, msec=.true.)   ! module "Little_subroutines"
   endif

!    print*, 'TEST 3'

   !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
enddo TP ! while (g_time .LT. g_numpar%t_total)


!'******************************************************************'
write(6,'(a)') trim(adjustl(m_starline))    ! module "Output"

! Execute all gnuplot scripts to plot the data user requested to:
if (g_numpar%gnupl%do_gnuplot) call execute_all_gnuplots(g_numpar) ! module "Output"

!--------------------------------------------------------------
! Finilize simulation run:
9999 continue
! Printing out the duration of the program, starting and ending time and date:
call date_and_time(values=g_c1)     ! For calculation of the time of execution of the program
g_as1=dble(24.0d0*60.0d0*60.0d0*(g_c1(3)-g_ctim(3))+3600.0d0*(g_c1(5)-g_ctim(5))+ &
    60.0d0*(g_c1(6)-g_ctim(6))+(g_c1(7)-g_ctim(7))+(g_c1(8)-g_ctim(8))*0.001d0)    ! sec
call parse_time(g_as1,g_text)   ! module "Little_subroutines"

! Print duration of execution of the program on the screen:
call print_duration(6, trim(adjustl(g_text)), g_ctim) ! module "Output"
! And into the Parameters file:
call print_duration(g_numpar%FN_parameters, trim(adjustl(g_text)), g_ctim) ! module "Output"

! Close remaining files:
call close_file('close', FN=g_numpar%FN_parameters)	! module "Dealing_with_files"
if (g_Err%Err) then
   call close_file('close', FN=g_Err%File_Num) ! module "Dealing_with_files"
else ! if there was no error, no need to keep the file, delete it
   call close_file('delete', FN=g_Err%File_Num) ! module "Dealing_with_files"
endif
! call close_file('close', File_name=g_numpar%FILE_communication)    ! module "Dealing_with_files"
STOP

!-----------------------------------------------------------

contains



pure subroutine reset_dt(numpar, tim_cur)
   real(8), intent(in) :: tim_cur        ! current time step of the simulation
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   real(8) :: est
   est = 1.0d-6 ! precision
   if ((numpar%i_dt > 0) .and. (numpar%i_dt <= size(numpar%dt_MD_grid))) then   ! only if there is an option to change dt
      if (tim_cur >= numpar%dt_MD_reset_grid(numpar%i_dt)-est) then ! time to change dt
         numpar%dt_MD = numpar%dt_MD_grid(numpar%i_dt)              ! to this value
         numpar%i_dt = numpar%i_dt + 1 ! next step to read from
      endif
   elseif (numpar%i_dt == 0) then   ! its before the simulation start, reset the starting time
      numpar%i_dt = numpar%i_dt + 1 ! next step to read from
      numpar%t_start = numpar%dt_MD_reset_grid(1)   ! to start from
      numpar%dt_MD = numpar%dt_MD_grid(1)           ! to start from
   endif
end subroutine reset_dt



pure subroutine reset_dt_out(numpar, dt_out)
   real(8), intent(inout) :: dt_out ! time to print out the data [fs]
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   real(8) :: est
   est = 1.0d-6 ! precision

   ! Check, if it's time to printout data into file:
   if (allocated(numpar%dt_out_grid)) then  ! only if a file with time grid was provided
      if (numpar%i_dt_printout <= size(numpar%dt_out_grid)) then    ! if didn't reach the end of file
         if (numpar%i_dt_printout > 1) then   ! not the first step:
            dt_out = numpar%dt_out_grid(numpar%i_dt_printout) - numpar%dt_out_grid(numpar%i_dt_printout-1) - est ! next step
         else   ! the first one
            dt_out = numpar%dt_out_grid(numpar%i_dt_printout) - est
         endif
         ! Check if the user requested to change the simulation step:
         if (allocated(numpar%dt_reset_grid)) then
            numpar%dt_MD = numpar%dt_reset_grid(numpar%i_dt_printout)   ! reset the MD time step according to the one from the file
            if (numpar%dt_MD > dt_out) numpar%dt_MD = dt_out    ! in case of inconsistency, set it equal to the printout step
         endif
         numpar%i_dt_printout = numpar%i_dt_printout + 1    ! next point on the time grid for next time
      else  ! end of the grid reached
         dt_out = numpar%t_total - numpar%dt_out_grid(size(numpar%dt_out_grid)) - est ! print the last point
      endif
   else ! given number:
      dt_out = numpar%dt_printout - est  ! next point on the equidistant grid with the constant step provided in NUMERICAL_PARAMETERS file
   endif
!    print*, 'dt_out', dt_out
end subroutine reset_dt_out


 
END PROGRAM TREKIS
