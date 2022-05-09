! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2021
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with files in XYZ format:
! [1] https://en.wikipedia.org/wiki/XYZ_file_format

MODULE Dealing_with_XYZ_files

use Universal_constants
use Objects
use Periodic_table, only: atomic_data, Read_from_periodic_table, find_atomic_number
use Dealing_with_files

implicit none

 contains


subroutine read_XYZ(FN, File_name, Path_to_Periodic, path_sep, Err, atoms, R, do_velosities, do_atomic_prprts)
   integer, intent(in) :: FN	! file number
   character(*), intent(in) :: File_name   ! file name
   character(*), intent(in) :: Path_to_Periodic ! path to the directory with periodic table:
   ! Path_to_Periodic = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_databases))
   character, intent(in) :: path_sep    ! path separator in the OS (windows vs linux)
   type(Error_handling), intent(inout) :: Err	! error log
   type(Atom), dimension(:), allocatable, intent(inout), optional :: atoms ! atomic parameters
   real(8), dimension(:,:), allocatable, intent(inout), optional :: R ! atomic coordinates in an array, if required
   logical, intent(in), optional :: do_velosities     ! reading atomic velosities instead of coordinates
   logical, intent(in), optional :: do_atomic_prprts  ! do we want to assign atomic properties from the Periodic Table?
   !---------------------------
   character(250) :: text, Path, Error_descript
   character(3) :: Name
   character(15) :: Full_name
   integer i, Nat, INFO, NVB, Reason, Z, count_lines
   logical :: assign_atomic_prprts, set_velosities, read_well
   type(atomic_data), dimension(:), allocatable :: Periodic_table_out   ! periodic table data, if user requested
   
   if (.not.present(atoms) .and. .not.present(R)) then
      Error_descript = 'Subroutine read_XYZ, user provided neither Atoms nor R arrays'
      call Save_error_details(Err, 5, Error_descript)	! module "Objects"
      print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
      goto 9999	! skip executing the program, exit the subroutin
   endif
   
   ! If user specified if we are reading coords or velosities:
   if (present(do_velosities)) then
      set_velosities = do_velosities
   else ! ba default, reading coordinates
      set_velosities = .false.
   endif
   
   ! If user specified if we need to set atomic parameters:
   if (present(atoms) .and. present(do_atomic_prprts)) then
      assign_atomic_prprts = do_atomic_prprts
   else ! ba default, no need to read properties from the Periodic Table
      assign_atomic_prprts = .false.
   endif
   
   
   ! Get the periodic table data to assign atomic numbers from names of elements:
   if (assign_atomic_prprts) then   ! only if user requested
      call Read_from_periodic_table(Path_to_Periodic, path_sep, INFO, Error_descript, 1, Periodic_table_out=Periodic_table_out)    ! module "Periodic_table"
      if (INFO /= 0) then
         Error_descript = 'Subroutine read_XYZ: '//trim(adjustl(Error_descript))
         call Save_error_details(Err, 4, Error_descript)	! module "Objects"
         print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
         goto 9999	! skip executing the program, exit the subroutin
      endif
   endif
   
   ! Reading first two line of XYZ file:
   count_lines = 0
   read_well = .true.
   read(FN,*,IOSTAT=Reason) Nat
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
      goto 9999
   endif
   
   read(FN,*,IOSTAT=Reason) text
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
      goto 9999
   endif
   
   ! Make sure atoms array has the same dimension as used provided:
   if (present(atoms)) then ! array of atoms as objects
      if (allocated(atoms)) deallocate(atoms)
      allocate(atoms(Nat))
      ! Read atomic names and coordinates in all further lines:
      do i = 1, Nat
         ! What to read from the file:
         if (set_velosities) then   ! velocties:
            read(FN,*,IOSTAT=Reason) Name, atoms(i)%V(1), atoms(i)%V(2), atoms(i)%V(3)
            if (trim(adjustl(Name)) /= trim(adjustl(atoms(i)%Name))) then
               Error_descript = 'Subroutine read_XYZ #3: inconsistency between atomic coordinates and velosities (different atom type of the same atom)'
               call Save_error_details(Err, 4, Error_descript)	! module "Objects"
               print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
               goto 9999	! skip executing the program, exit the subroutin
            endif
            atoms(i)%V0(:) = atoms(i)%V(:)  ! last timestep
         else   ! coordinates:
            read(FN,*,IOSTAT=Reason) atoms(i)%Name, atoms(i)%R(1), atoms(i)%R(2), atoms(i)%R(3)
            atoms(i)%R0(:) = atoms(i)%R(:)  ! last timestep
         endif
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i10)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            call Save_error_details(Err, 2, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9999
         endif
         
         ! Assign atomic properties for this element from its name:
         if (assign_atomic_prprts) then
            call find_atomic_number(atoms(i)%Name, atoms(i)%Z, NVB, Full_name, atoms(i)%Mass, Periodic_table_out, &
                                    INFO, Error_descript)   ! module "Periodic_table"
            if (INFO /= 0) then
               Error_descript = 'Subroutine read_XYZ #2: '//trim(adjustl(Error_descript))
               call Save_error_details(Err, 4, Error_descript)	! module "Objects"
               print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
               goto 9999	! skip executing the program, exit the subroutin
            endif
            atoms(i)%Mass = atoms(i)%Mass * g_amu  ! [a.m.u] -> [kg]
            
         endif ! (assign_atomic_prprts)
      enddo ! i = 1, Nat
   elseif (present(R)) then ! array of coordinatess
      if (allocated(R)) deallocate(R)
      allocate(R(Nat,3))
      ! Read atomic names and coordinates in all further lines:
      do i = 1, Nat
         read(FN,*,IOSTAT=Reason) Name, R(i,1), R(i,2), R(i,3)
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i10)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            call Save_error_details(Err, 2, Error_descript)	! module "Objects"
            print*, trim(adjustl(Error_descript)), ', TREKIS terminates'
            goto 9999
         endif
      enddo !  i = 1, Nat
   endif ! (present(atoms))
   
9999 continue   
end subroutine read_XYZ


subroutine write_XYZ(FN, atoms, text_to_write, print_velocities)
   integer, intent(in) :: FN	! file number
   type(Atom), dimension(:), intent(in) :: atoms	! atomic parameters
   character(*), intent(in) :: text_to_write    ! message to save in the XYZ block
   logical, intent(in), optional :: print_velocities    ! write velocities instead of coordinates
   integer i
   character(10) :: Numb_out
   logical :: print_coords
   
   print_coords = .true.    ! to start with
   if (present(print_velocities)) then  ! check if user wants to write out coordinates or velocities
      if (print_velocities) print_coords = .false.
   endif
   
   write(Numb_out, '(i10)') size(atoms)
   write(FN, '(a)') trim(adjustl(Numb_out))
   write(FN, '(a)') trim(adjustl(text_to_write))
   if (print_coords) then   ! write out coordinates
      do i = 1, size(atoms)
         write(FN, '(a,es25.16,es25.16,es25.16)') trim(adjustl(atoms(i)%Name)), atoms(i)%R(1), atoms(i)%R(2), atoms(i)%R(3)
      enddo
   else ! write out velocities:
      do i = 1, size(atoms)
         write(FN, '(a,es25.16,es25.16,es25.16)') trim(adjustl(atoms(i)%Name)), atoms(i)%V(1), atoms(i)%V(2), atoms(i)%V(3)
      enddo
   endif
end subroutine write_XYZ


END MODULE Dealing_with_XYZ_files
