! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2021
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to read input files for MD:

MODULE Read_MD_parameters
use Universal_constants
use Objects
use Periodic_table, only: atomic_data, Read_from_periodic_table, find_atomic_number
use Dealing_with_files, only: open_file, read_file, close_file !, Count_lines_in_file, Count_columns_in_file
use Dealing_with_XYZ_files, only: read_XYZ
use MD_general_tools, only: Sample_Maxwelian, remove_CoM_velocity, Check_MD_periodic_boundaries
use Little_subroutines, only: create_grid_1d


implicit none


! Names of the files with MD parameters:
character(50) :: m_supercell, m_coords_XYZ, m_vel_XYZ, m_unit_cell, m_potential

! All names are fixed here:
parameter(m_supercell   = 'Supercell.txt')       ! Parameters with supercell sizes and boundary conditions
parameter(m_coords_XYZ  = 'Coordinates.xyz')     ! Atomic coordinates [A]
parameter(m_vel_XYZ     = 'Velocities.xyz')      ! Atomic velocities [A/fs]
parameter(m_unit_cell   = 'Unit_cell.txt')        ! Atomic coordinates in unit cell, and unit cell parameters
parameter(m_potential   = 'Potential')           ! PArameters of all interatomic potentials / part of the name


 contains

subroutine read_supercell_parameters(Path_to_MD, numpar, MD_supce, Err)
   character(*), intent(in) :: Path_to_MD   ! path to the directory with MD parameters
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !-------------------
   character(250) :: File_supce, Error_descript
   character(50) :: temp_ch, temp_ch2
   integer :: FN, Reason, count_lines
   logical :: file_exists, read_well, MC2MD_descript_done

   Error_descript = ''  ! start with no error
   read_well = .true.   ! to start with
   MC2MD_descript_done = .false.    ! to start with
   ! Default values:
   numpar%use_thermostat = .false.  ! no thermostat to use for the supercell cooling
   numpar%thermostat_inx = -1       ! no thermostat
   numpar%damp_pressure = .false.   ! no pressure dampener to use
   MD_supce%V = -1.0d0  ! to start with (no supercell is defined yet)
   numpar%recenter = .false.    ! to start with
   MD_supce%coord_type = 0  ! to start with
   MD_supce%coord_dim = 0   ! to start with

   ! file with supercell parameters
   File_supce = trim(adjustl(Path_to_MD))//trim(adjustl(m_supercell))  ! Supercell.txt
   FN = 201
   call open_file(FN, File_supce, Error_descript, status = 'old', action='read')	! module "Dealing_with_files"
   if (LEN(trim(adjustl(Error_descript))) > 0) then	! it means, some error occured
      ! print out the error into the log file:
      call Save_error_details(Err, 1, Error_descript)	! module "Objects"
      print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
      goto 9999	! skip executing the program, exit the subroutine
   endif

   ! Read supercell parameters that user defined:
   count_lines = 0  ! to start with
   RSC:do while (read_well)	! check all possible entries
      read(FN,*,IOSTAT=Reason) temp_ch
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"

      if (.not. read_well) then
         exit RSC     ! nothing more to read from the file
      else  !  (.not. read_well)
         select case (trim(adjustl(temp_ch)))   ! block of parametres name
         case ('Supercell', 'SUPERCELL', 'supercell')  ! size of the supercell
            ! Read file with MD supercell parameters:
            count_lines = 0  ! to start with
            ! Supercell parameters along X-axis:
            read(FN,*,IOSTAT=Reason) MD_supce%x_start, MD_supce%x_end, MD_supce%boundary(1)
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9999
            endif

            ! Supercell parameters along Y-axis:
            read(FN,*,IOSTAT=Reason) MD_supce%y_start, MD_supce%y_end, MD_supce%boundary(2)
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9999
            endif
   
            ! Supercell parameters along Z-axis:
            read(FN,*,IOSTAT=Reason) MD_supce%z_start, MD_supce%z_end, MD_supce%boundary(3)
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9999
            endif
   
            ! Get the supercell sizes:
            MD_supce%vect(1) = abs(MD_supce%x_end - MD_supce%x_start)    ! [A]
            MD_supce%vect(2) = abs(MD_supce%y_end - MD_supce%y_start)    ! [A]
            MD_supce%vect(3) = abs(MD_supce%z_end - MD_supce%z_start)    ! [A]
            ! And volume:
            MD_supce%V = PRODUCT(MD_supce%vect)  ! [A^3]

         case ('Recenter', 'recenter', 'RECENTER')  ! recenter the supercell
            read(FN,*,IOSTAT=Reason) temp_ch2
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9999
            endif
            select case (temp_ch2)  ! define along which zxes to recenter the supercell
            case ('X', 'x')
               numpar%recenter(1) = .true.
            case ('Y', 'y')
               numpar%recenter(2) = .true.
            case ('Z', 'z')
               numpar%recenter(3) = .true.
            case ('XY', 'Xy', 'xY', 'xy', 'YX', 'Yx', 'yX', 'yx')
               numpar%recenter(1) = .true.
               numpar%recenter(2) = .true.
            case ('XZ', 'Xz', 'xZ', 'xz', 'ZX', 'zX', 'Zx', 'zx')
               numpar%recenter(1) = .true.
               numpar%recenter(3) = .true.
            case ('YZ', 'Yz', 'yZ', 'yz', 'ZY', 'Zy', 'zY', 'zy')
               numpar%recenter(2) = .true.
               numpar%recenter(3) = .true.
            case ('XYZ', 'YZX', 'ZXY', 'XZY', 'YXZ', 'ZYX', 'xyz', 'yzx', 'zxy', 'xzy', 'yxz', 'zyx')
               numpar%recenter(:) = .true.
            endselect

         case ('MC2MD', 'MC_MD_coupling', 'mc_md_coupling', 'MC_MD_Coupling', 'MC_MD_COUPLING')  ! grid for energy transfer between MC and MD

            ! Define the grid(s) for exchange of info between MC and MD:
            call read_MC2MD_energy_transfer(MD_supce, FN, File_supce, count_lines, read_well, Err)  ! below
            if (.not. read_well) goto 9999

         case ('Thermostat', 'THERMOSTAT', 'thermostat')   ! parameters of the thermostat
            read(FN,*,IOSTAT=Reason) temp_ch2
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9999
            endif
            select case (temp_ch2)
            case ('Berendsen', 'berendsen', 'BERENDSEN')
               numpar%use_thermostat = .true.  ! use thermostat in the MD simulation
               numpar%thermostat_inx = 0       ! Berendsen thermostat
               ! Thermostat (bath) temeprature and characteristic time:
               read(FN,*,IOSTAT=Reason) MD_supce%Bath_T, MD_supce%therm_tau  ! [K] bath temperature, [fs] characteristic time of cooling
               call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
               if (.not. read_well) then
                  write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
                  call Save_error_details(Err, 2, Error_descript)	! module "Objects"
                  goto 9999
               endif

               ! Thermostat boundaries along X-axis counted from the border of the supercell:
               read(FN,*,IOSTAT=Reason) MD_supce%thermo_dx_start, MD_supce%thermo_dx_end
               call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
               if (.not. read_well) then
                  write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
                  call Save_error_details(Err, 2, Error_descript)	! module "Objects"
                  goto 9999
               endif

               ! Thermostat boundaries along Y-axis counted from the border of the supercell:
               read(FN,*,IOSTAT=Reason) MD_supce%thermo_dy_start, MD_supce%thermo_dy_end
               call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
               if (.not. read_well) then
                  write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
                  call Save_error_details(Err, 2, Error_descript)	! module "Objects"
                  goto 9999
               endif

               ! Thermostat boundaries along Z-axis counted from the border of the supercell:
               read(FN,*,IOSTAT=Reason) MD_supce%thermo_dz_start, MD_supce%thermo_dz_end
               call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
               if (.not. read_well) then
                  write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
                  call Save_error_details(Err, 2, Error_descript)	! module "Objects"
                  goto 9999
               endif
            end select

         case ('Pressure_damp', 'Pressure_Damp', 'PRESSURE_DAMP', 'pressure_damp', & ! parameters of the pressure wave dampener
            'Pressure_dampener', 'Pressure_Dampener', 'PRESSURE_DAMPENER', 'pressure_dampener')
            numpar%damp_pressure = .true.

            ! Pressure dampening characteristic time [fs]:
            read(FN,*,IOSTAT=Reason) MD_supce%press_tau
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9999
            endif

            ! Pressure wave dampening boundaries along X-axis counted from the border of the supercell:
            read(FN,*,IOSTAT=Reason) MD_supce%pressdamp_dx_start, MD_supce%pressdamp_dx_end
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9999
            endif

            ! Pressure wave dampening boundaries along Y-axis counted from the border of the supercell:
            read(FN,*,IOSTAT=Reason) MD_supce%pressdamp_dy_start, MD_supce%pressdamp_dy_end
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9999
            endif

            ! Pressure wave dampening boundaries along Z-axis counted from the border of the supercell:
            read(FN,*,IOSTAT=Reason) MD_supce%pressdamp_dz_start, MD_supce%pressdamp_dz_end
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
               call Save_error_details(Err, 2, Error_descript)	! module "Objects"
               goto 9999
            endif

         end select
      endif ! (.not. read_well)
   enddo RSC
   !---------------------------

   ! Now we have supercell defined (probably):
   if (MD_supce%V <= 0.0d0) then ! supercell size is zero or negative, cannot proceed
      write(temp_ch,'(f20.3)') MD_supce%V
      write(Error_descript,'(a,i3)') 'Supercell is ill-defined (volume='//trim(adjustl(temp_ch))//&
                                      ') in the file '//trim(adjustl(File_supce))
      call Save_error_details(Err, 6, Error_descript)	! module "Objects"
      print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
      goto 9999
   endif

   ! close opened input files:
9999 continue
   call close_file('close', FN=FN)	! module "Dealing_with_files"

!    print*, MD_supce%x_start, MD_supce%x_end, MD_supce%boundary(1)
!    print*, MD_supce%y_start, MD_supce%y_end, MD_supce%boundary(2)
!    print*, MD_supce%z_start, MD_supce%z_end, MD_supce%boundary(3)
!    print*, MD_supce%V, MD_supce%vect(1)*MD_supce%vect(2)*MD_supce%vect(3)
!    pause 'read_supercell_parameters'
end subroutine read_supercell_parameters


subroutine read_MC2MD_energy_transfer(MD_supce, FN, File_supce, count_lines, read_well, Err)  ! below
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   integer, intent(in) :: FN
   character(*), intent(in) :: File_supce
   integer, intent(inout) :: count_lines
   logical, intent(inout) :: read_well
   type(Error_handling), intent(inout) :: Err	! error log
   !----------------------
   real(8) :: step(3), R_max, x_max, y_max, z_max
   integer :: FN2, Reason, i, axis(3), Nsiz, Nsiz2, Nsiz3
   character(15) :: coord_type
   character(250) :: Error_descript

   ! Which coordinate system to use, and what dimensionality:
   read(FN,*,IOSTAT=Reason) coord_type, MD_supce%coord_dim
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      print*, trim(adjustl(Error_descript)), ' TREKIS terminates'
      goto 9995
   endif

   ! Identify and save index of the coordinate system:
   select case (trim(adjustl(coord_type)))
   case ('Cartesian', 'cartesian', 'CARTESIAN')
      MD_supce%coord_type = 1
   case ('Spherical', 'spherical', 'SPHERICAL')
      MD_supce%coord_type = 2
   case ('Cylindrical', 'cylindrical', 'CYLINDRICAL')
      MD_supce%coord_type = 3
   case default
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' descriptor '// &
            trim(adjustl(coord_type))//' could not be interpreted'
      call Save_error_details(Err, 6, Error_descript)	! module "Objects"
      print*, trim(adjustl(Error_descript)), ' TREKIS terminates'
      goto 9995
   end select

   ! Get the grid parameters for MC-MD info exchange:
   if (MD_supce%coord_dim > 0) then ! at least 1d grid:
      do i = 1, MD_supce%coord_dim ! for all dimensions
         read(FN,*,IOSTAT=Reason) MD_supce%coord_axis(i), MD_supce%step(i)      ! along which axis and with which step [A] (or [rad] for angles)
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' could not read line ', count_lines
            call Save_error_details(Err, 2, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ' TREKIS terminates'
            goto 9995
         endif
         if (i > 1) then   ! check consistency:
            if (MD_supce%coord_axis(i) == MD_supce%coord_axis(i-1)) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' repeated line ', count_lines
               call Save_error_details(Err, 6, Error_descript)	! module "Objects"
               print*, trim(adjustl(Error_descript)), ' TREKIS terminates'
               goto 9995
            endif
         endif ! (i > 1)
         if (i > 2) then   ! check consistency:
            if (MD_supce%coord_axis(i) == MD_supce%coord_axis(i-2)) then
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_supce))//' repeated line ', count_lines
               call Save_error_details(Err, 6, Error_descript)	! module "Objects"
               print*, trim(adjustl(Error_descript)), ' TREKIS terminates'
               goto 9995
            endif
         endif ! (i > 2)
      enddo ! i
   endif ! (MD_supce%coord_dim > 0)

9995 continue
end subroutine read_MC2MD_energy_transfer


subroutine set_MC2MD_arrays(MD_supce)
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   !----------------------
   integer :: i, Nsiz, Nsiz2, Nsiz3
   real(8) :: x_max, y_max, z_max, R_max

   ! Get the grid parameters for MC-MD info exchange:
   if (MD_supce%coord_dim > 0) then ! at least 1d grid:
      do i = 1, MD_supce%coord_dim ! for all dimensions
         select case (MD_supce%coord_axis(i))  ! along which axis (for Cart, Spher, or Cyl)
         case (1)  ! 1=X or R or R
            select case (MD_supce%coord_type)
            case (1)   ! Cartesian: X
               call create_grid_1d(MD_supce%x_start, MD_supce%x_end, MD_supce%step(i), &
                                    MD_supce%MC2MD_dimen(i)%grid, 0)   ! module "Little_subroutines"
            case (2)   ! Spherical: R
               ! To encapsulate entire box inside of a shpere, use maximal distance:
               x_max = max(MD_supce%x_start, MD_supce%x_end)
               y_max = max(MD_supce%y_start, MD_supce%y_end)
               z_max = max(MD_supce%z_start, MD_supce%z_end)
               R_max = sqrt(x_max*x_max + y_max*y_max + z_max*z_max)
               call create_grid_1d(0.0d0, R_max, MD_supce%step(i), MD_supce%MC2MD_dimen(i)%grid, 0)   ! module "Little_subroutines"
            case (3)   ! Cylindrical: R
               ! To encapsulate entire box inside of a cylinder, use maximal distance:
               x_max = max(MD_supce%x_start, MD_supce%x_end)
               y_max = max(MD_supce%y_start, MD_supce%y_end)
               R_max = sqrt(x_max*x_max + y_max*y_max)
               call create_grid_1d(0.0d0, R_max, MD_supce%step(i), MD_supce%MC2MD_dimen(i)%grid, 0)   ! module "Little_subroutines"
            endselect

         case (2)  ! 2=Y or Theta or L
            select case (MD_supce%coord_type)
            case (1)   ! Cartesian: Y
               call create_grid_1d(MD_supce%y_start, MD_supce%y_end, MD_supce%step(i), &
                                    MD_supce%MC2MD_dimen(i)%grid, 0)   ! module "Little_subroutines"
            case (2)   ! Spherical: theta
               call create_grid_1d(0.0d0, g_2Pi, MD_supce%step(i), MD_supce%MC2MD_dimen(i)%grid, 0)   ! module "Little_subroutines"
            case (3)   ! Cylindrical: L
               call create_grid_1d(MD_supce%z_start, MD_supce%z_end, MD_supce%step(i), &
                                    MD_supce%MC2MD_dimen(i)%grid, 0)   ! module "Little_subroutines"
            endselect

         case (3)  ! 3=Z or Phi or Phi
            select case (MD_supce%coord_type)
            case (1)   ! Cartesian: Z
               call create_grid_1d(MD_supce%z_start, MD_supce%z_end, MD_supce%step(i), &
                                    MD_supce%MC2MD_dimen(i)%grid, 0)   ! module "Little_subroutines"
            case (2)   ! Spherical: phi
               call create_grid_1d(0.0d0, g_half_Pi, MD_supce%step(i), MD_supce%MC2MD_dimen(i)%grid, 0)   ! module "Little_subroutines"
            case (3)   ! Cylindrical: phi
               call create_grid_1d(0.0d0, g_2Pi, MD_supce%step(i), MD_supce%MC2MD_dimen(i)%grid, 0)   ! module "Little_subroutines"
            endselect
         end select
      enddo ! i
   endif ! (MD_supce%coord_dim > 0)

   ! Set the arrays with transferred energies between MC and MD:
   if (allocated(MD_supce%MC2MD_dimen(1)%grid)) then
      Nsiz = size(MD_supce%MC2MD_dimen(1)%grid)
   else
      Nsiz = 1
      allocate(MD_supce%MC2MD_dimen(1)%grid(1), source = 0.0d0)
   endif
   if (allocated(MD_supce%MC2MD_dimen(2)%grid)) then
      Nsiz2 = size(MD_supce%MC2MD_dimen(2)%grid)
   else
      Nsiz2 = 1
      allocate(MD_supce%MC2MD_dimen(2)%grid(1), source = 0.0d0)
   endif
   if (allocated(MD_supce%MC2MD_dimen(3)%grid)) then
      Nsiz3 = size(MD_supce%MC2MD_dimen(3)%grid)
   else
      Nsiz3 = 1
      allocate(MD_supce%MC2MD_dimen(3)%grid(1), source = 0.0d0)
   endif
   allocate(MD_supce%E_e_at_from_MC(Nsiz,Nsiz2,Nsiz3), source = 0.0d0)  ! energy elastically transferred e-to-atoms from within MC module
   allocate(MD_supce%E_h_at_from_MC(Nsiz,Nsiz2,Nsiz3), source = 0.0d0)  ! energy elastically transferred h-to-atoms from within MC module
   allocate(MD_supce%E_p_at_from_MC(Nsiz,Nsiz2,Nsiz3), source = 0.0d0)  ! energy elastically transferred p-to-atoms from within MC module
   allocate(MD_supce%E_e_from_MC(Nsiz,Nsiz2,Nsiz3), source = 0.0d0)   ! energy of electrons below cut-off from within MC module
   allocate(MD_supce%E_h_from_MC(Nsiz,Nsiz2,Nsiz3), source = 0.0d0)   ! energy of holes below cut-off from within MC module
end subroutine set_MC2MD_arrays



subroutine set_atomic_coordinates(numpar, Path_to_Periodic, Path_to_MD, MD_supce, MD_atoms, Err)
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   character(*), intent(in) :: Path_to_Periodic ! path to the directory with periodic table
   character(*), intent(in) :: Path_to_MD   ! path to the directory with MD parameters
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   type(Atom), dimension(:), allocatable, intent(inout) :: MD_atoms     ! all atoms in MD as objects
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------
   character(250) :: File_XYZ, File_unitcell, Error_descript
   character(30) :: temp_ch
   character(15) :: Full_name
   character(3), dimension(:), allocatable :: At_names
   real(8) :: coord_shift
   real(8), dimension(:,:), allocatable :: At_coords
   type(atomic_data), dimension(:), allocatable :: Periodic_table_out   ! periodic table data, if user requested
   real(8), dimension(3) :: Unit_cell_sizes
   integer, dimension(3) :: How_many_cells
   integer :: FN, FN2, Reason, count_lines, Nat, Nat_tot, i, j, k, at, count_at, INFO, NVB
   logical :: file_exists, file_exists2, read_well, do_xyz, do_unitcell, done_atomic, done_supce
   
   Error_descript = ''  ! start with no error
   ! file with supercell parameters
   File_XYZ = trim(adjustl(Path_to_MD))//trim(adjustl(m_coords_XYZ))  ! Coordinates.xyz
   FN = 202
   ! 1) check if file with XYZ coordinates exists (it has priority):
   inquire(file=trim(adjustl(File_XYZ)),exist=file_exists)
   if (file_exists) then    ! if we have XYZ file, use it
      do_xyz = .true.
   else ! Only if there is no XYZ data, check if there are data on unit cell:
      do_xyz = .false.
      File_unitcell = trim(adjustl(Path_to_MD))//trim(adjustl(m_unit_cell))  ! Uni_cell.txt
      FN2 = 203
      inquire(file=trim(adjustl(File_unitcell)),exist=file_exists2)
      if (file_exists2) then    ! if we have file with unit cell parameters, use it
         do_unitcell = .true.
      else
         Error_descript = 'Neither file '//trim(adjustl(File_XYZ))//', nor file '//trim(adjustl(File_XYZ))//' found'
         ! print out the error into the log file:
         call Save_error_details(Err, 1, Error_descript)	! module "Objects"
         print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         goto 9998	! skip executing the program, exit the subroutine
      endif
   endif
   
   ! Having found a file with coordinates, read it:
   if (do_xyz) then ! XYZ file with coordinates for all atoms in the supercell:
      call open_file(FN, File_XYZ, Error_descript, status = 'old', action='read')	! module "Dealing_with_files"
      if (LEN(trim(adjustl(Error_descript))) > 0) then	! it means, some error occured
         ! print out the error into the log file:
         call Save_error_details(Err, 1, Error_descript)	! module "Objects"
         print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         goto 9998	! skip executing the program, exit the subroutine
      endif
      
      ! Read coordinates:
      call read_XYZ(FN, File_XYZ, Path_to_Periodic, numpar%path_sep, Err, atoms=MD_atoms, do_atomic_prprts=.true.) ! module "Dealing_with_XYZ_files"
   
   elseif(do_unitcell) then ! file with parameters of a unit cell, to construct supercell:
      call open_file(FN2, File_unitcell, Error_descript, status = 'old', action='read')	! module "Dealing_with_files"
      if (LEN(trim(adjustl(Error_descript))) > 0) then	! it means, some error occured
         ! print out the error into the log file:
         call Save_error_details(Err, 1, Error_descript)	! module "Objects"
         print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         goto 9998	! skip executing the program, exit the subroutine
      endif
      
      ! Read parameters of unit cell and coordinates within this unit cell:
      ! to start with:
      done_atomic = .false.
      done_supce = .false.
      temp_ch = ''
      count_lines = 0
      read_well = .true. 
      ! Check what user defined:
      UC:do while (read_well)	! check all possible entries
         read(FN2,*,IOSTAT=Reason) temp_ch  ! check whether there are grids for printout
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            exit UC     ! nothing more to read from the file
         else  !  (.not. read_well)
            select case (trim(adjustl(temp_ch)))
            ! Atomic coordinates inside of a unit (primitive) cell:
            case ('Coordinates', 'COORDINATES', 'coordinates', 'COORDS', 'coords', 'Coords')
               read(FN2,*,IOSTAT=Reason) Nat ! number of atoms in the unit cell
               call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
               if (.not. read_well) then
                  write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_unitcell))//' could not read line ', count_lines
                  call Save_error_details(Err, 2, Error_descript)	! module "Objects"
                  goto 9998
               endif
               ! Knowing the number of atoms, proceed dealing with them:
               allocate(At_names(Nat))
               allocate(At_coords(Nat,3))
               ! For all atoms in the unit cell, read their names and coords:
               do i = 1, Nat
                  read(FN2,*,IOSTAT=Reason) At_names(i), At_coords(i,1), At_coords(i,2), At_coords(i,3)
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_unitcell))//' could not read line ', count_lines
                     call Save_error_details(Err, 2, Error_descript)	! module "Objects"
                     goto 9998
                  endif
               enddo
               done_atomic = .true.  ! read atomic data well
               
            ! Dimensions of the unit cell:
            case ('Unit_cell', 'UNIT_CELL', 'unit_cell', 'unitcell', 'UNITCELL', 'Unitcell')
               ! Size of the unit cell along X, Y, Z in [A]
               read(FN2,*,IOSTAT=Reason) Unit_cell_sizes(:)
               call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
               if (.not. read_well) then
                  write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_unitcell))//' could not read line ', count_lines
                  call Save_error_details(Err, 2, Error_descript)	! module "Objects"
                  goto 9998
               endif
               ! How many unit cells should be used to construct the supercell:
               read(FN2,*,IOSTAT=Reason) How_many_cells(:)
               call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
               if (.not. read_well) then
                  write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_unitcell))//' could not read line ', count_lines
                  call Save_error_details(Err, 2, Error_descript)	! module "Objects"
                  goto 9998
               endif
               done_supce = .true. ! read supercell parameters well
               
            end select
         endif ! (.not. read_well)
      enddo UC
      
      ! Make sure all necessary parameters have been read:
      if (.not.done_atomic .or. .not.done_supce) then
         Error_descript = 'Did not find correct MD parameters in the file '//trim(adjustl(File_unitcell))
         call Save_error_details(Err, 6, Error_descript)	! module "Objects"
         print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         print*, 'Please read the manual on how to set MD parameters for UNIT_CELL'
         goto 9998	! skip executing the program, exit the subroutine
      else  ! if all ok, set atomic coordinates for all atoms in the supercell:
         ! Get data from the Periodic table to set atomic parameters below:
         call Read_from_periodic_table(Path_to_Periodic, numpar%path_sep, INFO, Error_descript, &
                            1, Periodic_table_out=Periodic_table_out)    ! module "Periodic_table"
         if (INFO /= 0) then
            Error_descript = 'Subroutine set_atomic_coordinates: '//trim(adjustl(Error_descript))
            call Save_error_details(Err, 4, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9998	! skip executing the program, exit the subroutin
         endif
         
         ! Total number of atoms in the supercell, accounting for the number of unit cells:
         Nat_tot = Nat * How_many_cells(1) * How_many_cells(2) * How_many_cells(3)
         if (allocated(MD_atoms)) deallocate(MD_atoms)
         allocate(MD_atoms(Nat_tot))
         count_at = 0   ! to start with
         ! Now, set absolute coordinates of all atoms:
!          !$omp parallel private (i, j, k, at, count_at, NVB, Full_name, INFO, Error_descript)
!          !$omp do schedule(dynamic)
         do i = 1, How_many_cells(1)        ! X
            do j = 1, How_many_cells(2)     ! Y
               do k = 1, How_many_cells(3)  ! Z
                  do at = 1, Nat            ! atoms in a unit cell
                     ! To count atoms one by one:
                     !count_at = count_at + 1    ! count atoms
                     count_at = (Nat*(((i-1)*How_many_cells(2) + (j-1))*How_many_cells(3) + (k-1)) + at)

                     MD_atoms(count_at)%Name = At_names(at)
                     call find_atomic_number(MD_atoms(count_at)%Name, MD_atoms(count_at)%Z, NVB, Full_name, &
                        MD_atoms(count_at)%Mass, Periodic_table_out, INFO, Error_descript)   ! module "Periodic_table"
                     if (INFO /= 0) then
                        Error_descript = 'Subroutine set_atomic_coordinates #2: '//trim(adjustl(Error_descript))
                        call Save_error_details(Err, 4, Error_descript)	! module "Objects"
                        print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
                        goto 9998	! skip executing the program, exit the subroutin
                     endif
                     MD_atoms(count_at)%Mass = MD_atoms(count_at)%Mass * g_amu  ! [a.m.u] -> [kg]
                     MD_atoms(count_at)%R(1) = MD_supce%x_start + Unit_cell_sizes(1)*(dble(i - 1) + At_coords(at,1)) ! X [A]
                     MD_atoms(count_at)%R(2) = MD_supce%y_start + Unit_cell_sizes(2)*(dble(j - 1) + At_coords(at,2)) ! Y [A]
                     MD_atoms(count_at)%R(3) = MD_supce%z_start + Unit_cell_sizes(3)*(dble(k - 1) + At_coords(at,3)) ! Z [A]
                     MD_atoms(count_at)%R0(:) = MD_atoms(count_at)%R(:) ! last timestep
                  enddo ! at = 1, Nat
               enddo ! k = 1, How_many_cells(3)
            enddo ! j = 1, How_many_cells(2)
         enddo ! i = 1, How_many_cells(1)
!          !$omp enddo
!          !$omp end parallel

         ! Adjust the supercell size according to unit-cell parameters:
         MD_supce%x_end = MD_supce%x_start + Unit_cell_sizes(1)*How_many_cells(1)   ! X [A]
         MD_supce%y_end = MD_supce%y_start + Unit_cell_sizes(2)*How_many_cells(2)   ! Y [A]
         MD_supce%z_end = MD_supce%z_start + Unit_cell_sizes(3)*How_many_cells(3)   ! Z [A]
         ! Update the supercell sizes:
         MD_supce%vect(1) = abs(MD_supce%x_end - MD_supce%x_start)    ! [A]
         MD_supce%vect(2) = abs(MD_supce%y_end - MD_supce%y_start)    ! [A]
         MD_supce%vect(3) = abs(MD_supce%z_end - MD_supce%z_start)    ! [A]
         ! And volume:
         MD_supce%V = PRODUCT(MD_supce%vect)  ! [A^3]
      endif
   endif

   ! Check if the user requested to recenter the supercell:
   Nat = size(MD_atoms) ! total number of atoms in the supercell
   if (numpar%recenter(1)) then  ! recenter the supercell and shift coordinates along X
      coord_shift = -MD_supce%vect(1)*0.5d0 - MD_supce%x_start
      ! Place supercell symmetric around X=0:
      MD_supce%x_start = -MD_supce%vect(1)*0.5d0
      MD_supce%x_end   =  MD_supce%vect(1)*0.5d0
      ! And shift the atomic coordinates:
      !$omp parallel private (i)
      !$omp do schedule(dynamic)
      do i = 1, Nat     ! for all atoms
         MD_atoms(i)%R(1) = MD_atoms(i)%R(1) + coord_shift
         MD_atoms(i)%R0(1) = MD_atoms(i)%R(1)
      enddo
      !$omp enddo
      !$omp end parallel
   endif
   if (numpar%recenter(2)) then  ! recenter the supercell and shift coordinates along Y
      coord_shift = -MD_supce%vect(2)*0.5d0 - MD_supce%y_start
      ! Place supercell symmetric around Y=0:
      MD_supce%y_start = -MD_supce%vect(2)*0.5d0
      MD_supce%y_end   =  MD_supce%vect(2)*0.5d0
      ! And shift the atomic coordinates:
      !$omp parallel private (i)
      !$omp do schedule(dynamic)
      do i = 1, Nat     ! for all atoms
         MD_atoms(i)%R(2) = MD_atoms(i)%R(2) + coord_shift
         MD_atoms(i)%R0(2) = MD_atoms(i)%R(2)
      enddo
      !$omp enddo
      !$omp end parallel
   endif
   if (numpar%recenter(3)) then  ! recenter the supercell and shift coordinates along Z
      coord_shift = -MD_supce%vect(3)*0.5d0 - MD_supce%z_start
      ! Place supercell symmetric around Z=0:
      MD_supce%z_start = -MD_supce%vect(3)*0.5d0
      MD_supce%z_end   =  MD_supce%vect(3)*0.5d0
      ! And shift the atomic coordinates:
      !$omp parallel private (i)
      !$omp do schedule(dynamic)
      do i = 1, Nat     ! for all atoms
         MD_atoms(i)%R(3) = MD_atoms(i)%R(3) + coord_shift
         MD_atoms(i)%R0(3) = MD_atoms(i)%R(3)
      enddo
      !$omp enddo
      !$omp end parallel
   endif

   ! Make sure atoms are within the supercell:
   ! Place atoms back inside the supercell with periodic boundaries, if needed:
   call Check_MD_periodic_boundaries(MD_atoms, MD_supce)    ! module "MD_general_tools"

   ! Set arrays for energy exchange between MC and MD modules:
   call set_MC2MD_arrays(MD_supce)  ! below

   !---------------------------
   ! close opened input files:
9998 continue
   call close_file('close', FN=FN)	! module "Dealing_with_files"
   call close_file('close', FN=FN2)	! module "Dealing_with_files"
!    print*, '--------------'
!    print*, MD_supce%x_start, MD_supce%x_end, MD_supce%boundary(1)
!    print*, MD_supce%y_start, MD_supce%y_end, MD_supce%boundary(2)
!    print*, MD_supce%z_start, MD_supce%z_end, MD_supce%boundary(3)
!    print*, '--------------'
!    do i = 1, size(MD_atoms)
!       write(*,'(a,i4,a,f,f,f)') 'At:', i, MD_atoms(i)%Name, MD_atoms(i)%R(:)
!    enddo
!    pause 'set_atomic_coordinates'
end subroutine set_atomic_coordinates


subroutine set_atomic_velocities(used_target, numpar, Path_to_MD, MD_atoms, Err)
   type(Matter), intent(in) :: used_target	! parameters of the target
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   character(*), intent(in) :: Path_to_MD   ! path to the directory with MD parameters
   type(Atom), dimension(:), allocatable, intent(inout) :: MD_atoms     ! all atoms in MD as objects
   type(Error_handling), intent(inout) :: Err	! error log
   !-------------------------
   character(250) :: File_XYZ, Error_descript
   real(8) :: Ta_eV
   integer :: FN, Nat, i
   logical :: file_exists, do_xyz
   
   Error_descript = ''  ! start with no error
   ! file with supercell parameters
   File_XYZ = trim(adjustl(Path_to_MD))//trim(adjustl(m_vel_XYZ))  ! Velocities.xyz
   FN = 203
   ! 1) check if file with XYZ coordinates exists (it has priority):
   inquire(file=trim(adjustl(File_XYZ)),exist=file_exists)
   if (file_exists) then    ! if we have XYZ file, use it
      do_xyz = .true.
   else ! Only if there is no XYZ data, check if there are data on unit cell:
      do_xyz = .false.
   endif
   
   ! Having found a file with coordinates, read it:
   if (do_xyz) then ! XYZ file with coordinates for all atoms in the supercell:
      call open_file(FN, File_XYZ, Error_descript, status = 'old', action='read')	! module "Dealing_with_files"
      if (LEN(trim(adjustl(Error_descript))) > 0) then	! it means, some error occured
         ! print out the error into the log file:
         call Save_error_details(Err, 1, Error_descript)	! module "Objects"
         print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         goto 9998  ! skip executing the program, exit the subroutine
      endif
      
      ! Read coordinates:
      call read_XYZ(FN, File_XYZ, File_XYZ, numpar%path_sep, Err, atoms=MD_atoms, do_velosities=.true.) ! module "Dealing_with_XYZ_files"
   else ! Set velocities randomly accodring to Maxwellian distribution:
      Nat = size(MD_atoms)  ! total number of atoms in the supercell      
      ! Using temperature of the FIRST target:
      ! Sample atomic velosities:
      !$omp parallel private (i)
      !$omp do schedule(dynamic)
      do i = 1, Nat     ! for all atoms
         call Sample_Maxwelian(used_target%Material(1)%T_eV, MD_atoms(i)%Mass, &
                            MD_atoms(i)%V(1), MD_atoms(i)%V(2), MD_atoms(i)%V(3))    ! module "MD_general_tools"
      enddo
      !$omp enddo
      !$omp end parallel
   endif
   
   ! Set "last" time step velosities:
   !$omp parallel private (i)
   !$omp do schedule(dynamic)
   do i = 1, Nat     ! for all atoms
      MD_atoms(i)%V0(:) = MD_atoms(i)%V(:)
   enddo
   !$omp enddo
   !$omp end parallel
   
   ! Get rid of the center-of-mass motion:
   call remove_CoM_velocity(MD_atoms)   ! module "MD_general_tools"
   
   ! close opened input files:
9998 continue
   call close_file('close', FN=FN)	! module "Dealing_with_files"
!    print*, '--------------'
!    do i = 1, size(MD_atoms)
!       write(*,'(a,i4,a,i,es,f,f,f)') 'V:', i, MD_atoms(i)%Name, MD_atoms(i)%Z, MD_atoms(i)%Mass, MD_atoms(i)%V(:)
!    enddo
!    pause 'set_atomic_velocities'
end subroutine set_atomic_velocities


subroutine set_MD_potential(Path_to_MD, numpar, used_target, ind_target, MD_pots, Err)
   character(*), intent(in) :: Path_to_MD   ! path to the directory with MD parameters
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(Matter), intent(in), target :: used_target	! parameters of the target
   integer, intent(in) :: ind_target        ! target index
   type(MD_potential), dimension(:,:), allocatable, intent(inout) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(Error_handling), intent(inout) :: Err	! error log
   !-------------------------
   character(500) :: File_Pot, Error_descript, File_Pot1, File_Pot2
   integer :: FN, NoElem, i, j, k
   logical :: file_exists
   character(3), pointer :: Elem_name1, Elem_name2
   
   ! to start with:
   Error_descript = ''
   FN = 5000

   ! Number of different chemical elements the target is composed of:
   NoElem = size(used_target%Material(ind_target)%Elements(:))
   ! Allocate array of pair-wise potentials:
   allocate(MD_pots(NoElem,NoElem))
   
   ! For each pair of elements, find files with parameters defining the potential:
   do i = 1, NoElem ! for all elements
      Elem_name1 => used_target%Material(ind_target)%Elements(i)%Name   ! name of the first element
      do j = i, NoElem  ! for elements that were not checked yet
         ! Find file with parameters for this pair of elements:
         Elem_name2 => used_target%Material(ind_target)%Elements(j)%Name    ! name of the second element
         ! Construct the name of the file with the potential (elements could be in any order):
         File_Pot1 = trim(adjustl(Path_to_MD))//trim(adjustl(Elem_name1))//'_'// &
                    trim(adjustl(Elem_name2))//'_'//trim(adjustl(m_potential))//'.txt'  ! element1_element2
         File_Pot2 = trim(adjustl(Path_to_MD))//trim(adjustl(Elem_name2))//'_'// &
                    trim(adjustl(Elem_name1))//'_'//trim(adjustl(m_potential))//'.txt'  ! element2_element1
         ! Check if there is such a file:
         inquire(file=trim(adjustl(File_Pot1)),exist=file_exists)
         if (file_exists) then
            File_Pot = File_Pot1    ! file found, proceed
         else ! check if the other order:
            ! Check if there is such a file:
            inquire(file=trim(adjustl(File_Pot2)),exist=file_exists)
            if (file_exists) then 
               File_Pot = File_Pot2    ! file found, proceed
            else ! check if there is a file with a general name:
               Error_descript = 'Could not find file with MD potential: '//trim(adjustl(File_Pot))
               print*, trim(adjustl(Error_descript)), ', checking file '//trim(adjustl(m_potential))//'.txt'
               ! Try the general file:
               File_Pot = trim(adjustl(Path_to_MD))//trim(adjustl(m_potential))//'.txt'
               inquire(file=trim(adjustl(File_Pot)),exist=file_exists)
               if (.not.file_exists) then  ! no file with the potential
                  ! print out the error into the log file:
                  Error_descript = trim(adjustl(Error_descript))//' or '//trim(adjustl(File_Pot))
                  call Save_error_details(Err, 1, Error_descript)	! module "Objects"
                  print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
                  goto 9997  ! skip executing the program, exit the subroutine
               endif
            endif
         endif
         
         ! If file with the potential found, proceed with reading it:
         FN = 5000+i    ! file number
         call open_file(FN, File_Pot, Error_descript, status = 'old', action='read')	! module "Dealing_with_files"
         if (LEN(trim(adjustl(Error_descript))) > 0) then	! it means, some error occured
            ! print out the error into the log file:
            call Save_error_details(Err, 1, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9997  ! skip executing the program, exit the subroutine
         endif
         
         ! Save names of elements for this potential:
         MD_pots(i,j)%El1 = Elem_name1
         MD_pots(i,j)%El2 = Elem_name2
         
         call read_MD_potential(FN, File_Pot, numpar, MD_pots(i, j), Error_descript)   ! below
         if (LEN(trim(adjustl(Error_descript))) > 0) then	! it means, some error occured
            ! print out the error into the log file:
            call Save_error_details(Err, 1, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9997  ! skip executing the program, exit the subroutine
         endif   
         
         ! When done with the file, close it:
         call close_file('close', FN=FN)	! module "Dealing_with_files"

         
         ! Lower triangle is the same as the upper one:
         if (j > i) then
            MD_pots(j,i)%El1 = MD_pots(i,j)%El2
            MD_pots(j,i)%El2 = MD_pots(i,j)%El1
            allocate(MD_pots(j,i)%Set(size(MD_pots(i,j)%Set)))
            do k = 1, size(MD_pots(i,j)%Set)
               ASSOCIATE (ARRAY => MD_pots(i,j)%Set(k)%Par) ! this is the sintax we have to use to check the class of defined types
                  select type(ARRAY)
                  type is (LJ) ! set all the same parameters of the potential
                     allocate(LJ::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (LJ)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%C12 = ARRAY%C12
                           ARRAY2%C6 = ARRAY%C6
                        endselect
                     ENDASSOCIATE
                  type is (Buck)
                     allocate(Buck::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (Buck)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%A = ARRAY%A
                           ARRAY2%B = ARRAY%B
                           ARRAY2%C = ARRAY%C
                        endselect
                     ENDASSOCIATE
                  type is (Power)
                     allocate(Power::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (Power)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%n = ARRAY%n
                           ARRAY2%C = ARRAY%C
                        endselect
                     ENDASSOCIATE
                  type is (Exp_law)
                     allocate(Exp_law::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (Exp_law)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%k = ARRAY%k
                           ARRAY2%C = ARRAY%C
                        endselect
                     ENDASSOCIATE
                  type is (Matsui)
                     allocate(Matsui::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (Matsui)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%A = ARRAY%A
                           ARRAY2%B = ARRAY%B
                           ARRAY2%C = ARRAY%C
                        endselect
                     ENDASSOCIATE
                  type is (SW)
                     allocate(SW::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (SW)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%E = ARRAY%E
                           ARRAY2%sigma = ARRAY%sigma
                           ARRAY2%A = ARRAY%A
                           ARRAY2%B = ARRAY%B
                           ARRAY2%p = ARRAY%p
                           ARRAY2%q = ARRAY%q
                           ARRAY2%aa = ARRAY%aa
                           ARRAY2%ll = ARRAY%ll
                           ARRAY2%gam = ARRAY%gam
                        endselect
                     ENDASSOCIATE
                  type is (Morse)
                     allocate(Morse::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (Morse)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%De = ARRAY%De
                           ARRAY2%re = ARRAY%re
                           ARRAY2%a = ARRAY%a
                        endselect
                     ENDASSOCIATE
                  type is (Coulomb)
                     allocate(Coulomb::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (Coulomb)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%Z2 = ARRAY%Z1
                           ARRAY2%Z1 = ARRAY%Z2
                        endselect
                     ENDASSOCIATE
                  type is (Coulomb_screened)
                     allocate(Coulomb_screened::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (Coulomb_screened)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%Z2 = ARRAY%Z1
                           ARRAY2%Z1 = ARRAY%Z2
                        endselect
                     ENDASSOCIATE
                  type is (Soft_Coulomb)
                     allocate(Soft_Coulomb::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (Soft_Coulomb)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%Z2 = ARRAY%Z1
                           ARRAY2%Z1 = ARRAY%Z2
                           ARRAY2%r0 = ARRAY%r0
                        endselect
                     ENDASSOCIATE
                  type is (Coulomb_Ewald)
                     allocate(Coulomb_Ewald::MD_pots(j,i)%Set(k)%Par)
                     ASSOCIATE (ARRAY2 => MD_pots(j,i)%Set(k)%Par)
                        select type(ARRAY2)
                        type is (Coulomb_Ewald)
                           ARRAY2%d_cut = ARRAY%d_cut
                           ARRAY2%dd = ARRAY%dd
                           ARRAY2%Z2 = ARRAY%Z1
                           ARRAY2%Z1 = ARRAY%Z2
                        endselect
                     ENDASSOCIATE
                  type is (ZBL)
                     allocate(ZBL::MD_pots(j,i)%Set(k)%Par)
                     MD_pots(j,i)%Set(k)%Par%d_cut = ARRAY%d_cut
                     MD_pots(j,i)%Set(k)%Par%dd = ARRAY%dd
                  endselect ! select type(ARRAY)
               ENDASSOCIATE ! (ARRAY => MD_pots(i,j)%Set(k)%Par)
            enddo ! k
         endif ! (j > i)
      enddo ! j = i, NoElem
   enddo ! i

   ! close opened input files:
9997 continue
   call close_file('close', FN=FN)	! module "Dealing_with_files"
   nullify(Elem_name1, Elem_name2)
end subroutine set_MD_potential



subroutine read_MD_potential(FN, File_name, numpar, MD_pots, Error_descript)
   integer, intent(in) :: FN    ! file with potential, must be opened
   character(*), intent(inout) :: File_name
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   type(MD_potential), intent(inout) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   character(*), intent(inout) :: Error_descript
   !---------------------
   integer :: NoP, Reason, count_lines, count_pots
   character(50) :: temp_ch
   logical :: read_well, read_potential
   
   count_lines = 0	! to start counting lines in the file

   ! First line defines number of potentials:
   read(FN,*,IOSTAT=Reason) NoP
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      goto 9996
   endif
   
   ! When we know how many potentials we have, allocate the arrays:
   allocate(MD_pots%Set(NoP))

   ! Read potentials that user defined:
   count_pots = 1   ! to start with
   PF:do while (read_well)	! check all possible entries
      read(FN,*,IOSTAT=Reason) temp_ch  ! check whether there are grids for printout
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         exit PF     ! nothing more to read from the file
      else  !  (.not. read_well)
         read_potential = .false. ! did not read any potential yet
         select case (trim(adjustl(temp_ch)))   ! potential name
         case ('LJ', 'Lj', 'lj', 'Lennard_Jones', 'Lennard-Jones', 'LennardJones')    ! Lennard-Jones
            allocate(LJ::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (LJ)
                  ARRAY%Name = 'Lennard-Jones'
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%C12, ARRAY%C6
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
!                   print*, 'Potential:', trim(adjustl(ARRAY%Name)), ARRAY%d_cut, ARRAY%dd, ARRAY%C12, ARRAY%C6
!                   pause 'read_MD_potential'
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential

         case ('Power_law', 'Power', 'Power_Law', 'POWER', 'POWER_LAW', 'power', 'power_law')  ! power law
            allocate(Power::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (Power)
                  ARRAY%Name = 'Power_law'

                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif

                  read(FN,*,IOSTAT=Reason) ARRAY%C, ARRAY%n
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential

         case ('Exp_law', 'Exponential', 'Exp_Law', 'EXPONENTIAL', 'EXP_LAW', 'exponential', 'exp_law')  ! exponential law
            allocate(Exp_law::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (Exp_law)
                  ARRAY%Name = 'Exponential_law'

                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif

                  read(FN,*,IOSTAT=Reason) ARRAY%C, ARRAY%k
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential

         case ('Buck', 'BUCK', 'buck', 'Buckingham', 'BUCKINGHAM', 'buckingham')  ! Buckingham
            allocate(Buck::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (Buck)
                  ARRAY%Name = 'Buckingham'
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%A, ARRAY%B, ARRAY%C
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential
         
         case ('Matsui', 'MATSUI', 'matsui')  ! Matsui
            allocate(Matsui::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (Matsui)
                  ARRAY%Name = 'Matsui'
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%A, ARRAY%B, ARRAY%C
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential

         case ('SW', 'sw', 'Stillinger-Weber', 'StillingerWeber', 'Stillinger_Weber')  ! Stillinger-Weber
            allocate(SW::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (SW)
                  ARRAY%Name = 'Stillinger-Weber'

                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif

                  read(FN,*,IOSTAT=Reason) ARRAY%E, ARRAY%sigma, ARRAY%aa
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif

                  read(FN,*,IOSTAT=Reason) ARRAY%A, ARRAY%B
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif

                  read(FN,*,IOSTAT=Reason) ARRAY%p, ARRAY%q
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif

                  read(FN,*,IOSTAT=Reason) ARRAY%ll, ARRAY%gam
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential

            ! For now, only numberical derivative of SW is implemented, so mark it:
!              numpar%MD_force_ind = 1

         case ('Morse', 'MORSE', 'morse')  ! Morse
            allocate(Morse::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (Morse)
                  ARRAY%Name = 'Morse'
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%De, ARRAY%re, ARRAY%a
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential

         case ('Coulomb', 'COULOMB', 'coulomb')  ! Coulomb
            allocate(Coulomb::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (Coulomb)
                  ARRAY%Name = 'Coulomb'
               
                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
                  ARRAY%dd = 0.0d0  ! for Wolf's method, hard cut-off must be used
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%Z1, ARRAY%Z2
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential

         case ('Coulomb_Ewald', 'COULOMB_EWALD', 'coulomb_ewald', 'Coulomb_Ewalds', 'EOULOMB_EWALDS', 'eoulomb_ewalds')  ! Coulomb Ewald
            allocate(Coulomb_Ewald::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (Coulomb_Ewald)
                  ARRAY%Name = 'Coulomb_Ewald'

                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif

                  read(FN,*,IOSTAT=Reason) ARRAY%Z1, ARRAY%Z2
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential


         case ('Coulomb_screen', 'COULOMB_SCREEN', 'coulomb_screen', 'Coulomb_screened', 'COULOMB_SCREENED', 'coulomb_screened', 'Screened_Coulomb', 'SCREENED_COULOMB', 'screened_coulomb')  ! SHort-range screened Coulomb
            allocate(Coulomb_screened::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (Coulomb_screened)
                  ARRAY%Name = 'Screened Coulomb'
               
                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%Z1, ARRAY%Z2
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential
            
         case ('Soft_Coulomb', 'Soft_coulomb', 'soft_soulomb', 'SOFT_COULOMB')  ! Soft Coulomb
            allocate(Soft_Coulomb::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (Soft_Coulomb)
                  ARRAY%Name = 'Soft Coulomb'
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
                  
                  read(FN,*,IOSTAT=Reason) ARRAY%Z1, ARRAY%Z2, ARRAY%r0
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential

         case ('ZBL', 'zbl')  ! Universal ZBL
            allocate(ZBL::MD_pots%Set(count_pots)%Par)    ! set the type of potential
            ASSOCIATE (ARRAY => MD_pots%Set(count_pots)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(ARRAY)
               type is (ZBL)
                  ARRAY%Name = 'ZBL'
                  read(FN,*,IOSTAT=Reason) ARRAY%d_cut, ARRAY%dd
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not. read_well) then
                     write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
                     goto 9996
                  endif
               endselect
            ENDASSOCIATE
            read_potential = .true. ! read one more potential
            
         end select
         ! If we've read one more potenital from the file:
         if (read_potential) count_pots = count_pots + 1    ! next potential in the file
         
      endif ! (.not. read_well)
   enddo PF

9996 continue
end subroutine read_MD_potential

END MODULE Read_MD_parameters
