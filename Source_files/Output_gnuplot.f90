! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2026
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to create gnuplot shell scripts:

MODULE Output_gnuplot

! Open_MP related modules from external libraries:
#ifdef _OPENMP
   USE IFLPORT
   USE OMP_LIB
#endif

use Objects
use Universal_constants
use Gnuplotting
use Little_subroutines, only: print_time, find_order_of_number
use Dealing_with_files, only: copy_file, close_file, get_file_stat, Count_lines_in_file, read_file


implicit none

character(100) :: m_output_MD_energies
character(100) :: m_output_MD_cell_params
character(100) :: m_output_MD_displacements
character(100) :: m_output_MD_E_gnu
character(100) :: m_output_MD_T_gnu
character(100) :: m_output_MD_MSD_gnu
character(100) :: m_output_total
character(100) :: m_output_total_cutoff
character(100) :: m_output_N_gnu
character(100) :: m_output_E_gnu
character(100) :: m_output_MD

character(100) :: m_folder_MFP, m_output_Compton, m_output_Rayleigh, m_output_pair, m_output_absorb, m_output_MFP, m_output_IMFP, &
                  m_output_EMFP, m_output_Brems, m_output_annihil, m_output_Range, m_output_Se, m_output_Se_vs_range, &
                  m_output_DOS, m_output_DOS_k, m_output_DOS_effm

parameter (m_output_total = 'OUTPUT_total_')
parameter (m_output_total_cutoff = 'OUTPUT_total_above_cutoff_')
parameter (m_output_N_gnu = 'OUTPUT_total_numbers_')
parameter (m_output_E_gnu = 'OUTPUT_total_energies_')

parameter (m_output_MD = 'OUTPUT_MD_')
parameter (m_output_MD_energies = 'OUTPUT_MD_energies.txt')
parameter (m_output_MD_cell_params = 'OUTPUT_MD_average_parameters.txt')
parameter (m_output_MD_displacements = 'OUTPUT_MD_mean_displacements.txt')
parameter (m_output_MD_E_gnu = 'OUTPUT_MD_energies_')
parameter (m_output_MD_T_gnu = 'OUTPUT_MD_temperature_')
parameter (m_output_MD_MSD_gnu = 'OUTPUT_MD_displacements_')

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

parameter (m_output_DOS = 'OUTPUT_DOS_of_')
parameter (m_output_DOS_k = 'OUTPUT_DOS_k_vector_of_')
parameter (m_output_DOS_effm = 'OUTPUT_DOS_effective_mass_of_')


contains



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







END MODULE Output_gnuplot
