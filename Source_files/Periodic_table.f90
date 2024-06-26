! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2015-2021
! 1111111111111111111111111111111111111111111111111111111111111
MODULE Periodic_table

implicit none

!==============================================
! For reading atomic data from our periodic table:
type Atomic_data    ! our internal atomic database "INPUT_atomic_data.dat"
   integer :: Z     ! atomic number
   real(8) :: Mass  ! atomic mass [atomic units]
   character(15) :: Full_Name   ! Full atomic name
   character(3) :: Name         ! element name 
   real(8) :: Nvb   ! number of valence electrons
   real(8) :: cov_r ! covalent radius
endtype Atomic_data
!==============================================

 contains
 

! Find data on the given element from the periodic table:
subroutine Read_from_periodic_table(Path, path_sep, INFO, error_message, Z, name, elem_name, mass, NVB, cov_R, Periodic_table_out)
   character(*), intent(in) :: Path ! path to the folder with 'INPUT_atomic_data.dat'
   character(*), intent(in) :: path_sep ! path separator
   integer, intent(inout) :: INFO  ! 0=file read well; 1=no file; 2=couldn't open; 3=error while reading
   character(100), intent(inout) :: error_message
   integer, intent(in) :: Z     ! element we are looking data for
   character(*), intent(out), optional :: name    ! name of the element
   character(*), intent(out), optional :: elem_name ! full name of the element
   real(8), intent(out), optional :: mass     ! mass of the element [a.m.u.]
   integer, intent(out), optional :: NVB ! number of valence electrons
   real(8), intent(out), optional :: cov_R  ! covalent radius [A]
   type(atomic_data), dimension(:), allocatable, intent(inout), optional :: Periodic_table_out   ! the entire periodic table to print out for the user
   !==============================================
   type(atomic_data), dimension(:), allocatable :: Periodic_table ! this is an internal module variable
   integer FN, Reason, i, N
   character(100) :: Folder_name, File_name
   logical :: file_exists, file_opened
   INFO = 0 ! start with no error
   Folder_name = trim(adjustl(Path))    ! where to find the periodic table
   File_name = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//'INPUT_atomic_data.dat' ! fixed name of the database
   
   inquire(file=trim(adjustl(File_name)),exist=file_exists) ! check if input file is there
   exists:if (file_exists) then   
      FN = 101
      open (unit=FN, file=trim(adjustl(File_name)))
      inquire(unit=FN,opened=file_opened)    ! check if this file is opened
      if (file_opened) then
         call Count_lines(FN, N)    ! below
         if (.not.allocated(Periodic_table)) allocate(Periodic_table(N-1))
         if (present(Periodic_table_out)) then
            if (.not.allocated(Periodic_table_out)) allocate(Periodic_table_out(N-1))
         endif
         do i = 1, N-1
            read(FN,*,IOSTAT=Reason) Periodic_table(i)%Z, Periodic_table(i)%Full_name, Periodic_table(i)%Name, Periodic_table(i)%Mass, Periodic_table(i)%Nvb, Periodic_table(i)%cov_r
            !print*, Periodic_table(i)%Z, Periodic_table(i)%Full_name, Periodic_table(i)%Name, Periodic_table(i)%Mass, Periodic_table(i)%Nvb
            if (Reason /= 0) exit
         enddo
         close(FN)             ! and if it is, close it
         if (Reason /= 0) goto 912
      else
         write(error_message,*) 'Could not open the file: ', trim(adjustl(File_name))
         INFO = 2 ! could not open file
         goto 912
      endif
      ! Printout the requested parameters:
      if (present(name)) name =  Periodic_table(Z)%Full_name
      if (present(elem_name)) elem_name =  Periodic_table(Z)%Name
      if (present(mass)) mass =  Periodic_table(Z)%Mass
      if (present(NVB)) NVB =  Periodic_table(Z)%NVB
      if (present(cov_R)) cov_R =  Periodic_table(Z)%cov_r
      if (present(Periodic_table_out)) Periodic_table_out = Periodic_table

   else exists
      write(error_message,*) 'Could not find the file: ', trim(adjustl(File_name))
      Reason = 1
      INFO = 1
   endif exists

912 if (allocated(Periodic_table)) deallocate(Periodic_table)
   if (Reason /= 0) then
      write(error_message,*) 'Could not read input file: ', trim(adjustl(File_name))
      INFO = 3
   endif
end subroutine Read_from_periodic_table


 

! Find which elements constitute given chemical formula:
subroutine Decompose_compound(Path, Name, path_sep, INFO, error_message, at_num, at_numbers, &
                                at_percentage, at_short_names, at_names, at_masses, at_NVB)
   character(*), intent(in) :: Path ! path to the folder with 'INPUT_atomic_data.dat'
   character(*), intent(in) :: Name ! compound name (SiO2, Al2O3, etc.)
   character(*), intent(in) :: path_sep ! path separator
   integer, intent(inout) :: INFO  ! 0=file read well; 1=no file; 2=couldn't open; 3=error while reading; 4=user-requested exit
   character(100), intent(inout) :: error_message
   integer, intent(out) :: at_num ! how many different elements are in this compound
   integer, dimension(:), allocatable, intent(out), optional :: at_numbers
   real(8), dimension(:), allocatable, intent(out), optional :: at_percentage
   !character(*), dimension(:), allocatable, intent(out), optional :: at_short_names ! name of the element
   character(3), dimension(:), allocatable, intent(out), optional :: at_short_names ! name of the element
   !character(*), dimension(:), allocatable, intent(out), optional :: at_names ! full name of the element
   character(3), dimension(:), allocatable, intent(out), optional :: at_names ! full name of the element
   real(8), dimension(:), allocatable, intent(out), optional :: at_masses ! mass of each element [Mp]
   integer, dimension(:), allocatable, intent(out), optional :: at_NVB ! number of valence electrons
   !==============================================
   type(atomic_data), dimension(:), allocatable :: Periodic_table ! this is an internal module variable
   integer FN, Reason, leng, i, cur, N, num, coun, NVB
   real(8) C
   character(3), dimension(100) :: ElNames ! all elements in the compound, assuming they are not more then 100
   real(8), dimension(100) :: ElPersent     ! persentage of this element in the compound
   integer, dimension(100) :: ElNumbers     ! numbers of elements in the compound
   character(100) :: Folder_name, File_name
   character(15) :: Full_name
   real(8) :: M
   character(3) :: El
   character(*), parameter :: numbers = '0123456789'
   CHARACTER(*), PARAMETER :: LowCase = 'abcdefghijklmnopqrstuvwxyz'
   CHARACTER(*), PARAMETER :: UpCase  = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   logical :: file_exists, num_vs_char, found_atom, file_opened, devide
   INFO = 0 ! start with no error
   !Folder_name = 'INPUT_DATA'//trim(adjustl(path_sep))//'EADL_parameters'  ! here we keep databases
   Folder_name = trim(adjustl(Path))
   File_name = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//'INPUT_atomic_data.dat' ! fixed name of the database
   inquire(file=trim(adjustl(File_name)),exist=file_exists) ! check if input file is there
   exists:if (file_exists) then   
      FN = 101
      open (unit=FN, file=trim(adjustl(File_name)))
      !open (newunit=FN, file=trim(adjustl(File_name)))
      inquire(unit=FN,opened=file_opened)    ! check if this file is opened
      if (file_opened) then
         call Count_lines(FN, N)    ! below
         if (.not.allocated(Periodic_table)) allocate(Periodic_table(N-1))
!          read(FN,*,IOSTAT=Reason) ! skip first line with names
         do i = 1, N-1
            read(FN,*,IOSTAT=Reason) Periodic_table(i)%Z, Periodic_table(i)%Full_name, Periodic_table(i)%Name, Periodic_table(i)%Mass, Periodic_table(i)%Nvb
            !print*, Periodic_table(i)%Z, Periodic_table(i)%Full_name, Periodic_table(i)%Name, Periodic_table(i)%Mass, Periodic_table(i)%Nvb
            if (Reason /= 0) exit
         enddo
         close(FN)             ! and if it is, close it
         if (Reason /= 0) goto 911
      else
         INFO = 2 ! could not open file
         goto 911
      endif

      devide = .false.
      coun = 0
      leng = LEN(trim(adjustl(Name))) ! how many characters are in the name

      num = 0 ! how many different elements are in this compound
      El = ' ' ! start a new name
      do i = 1,leng ! compare all name character by character
         if (verify(trim(adjustl(Name(i:i))),trim(adjustl(numbers))) == 0) then ! it's an integer number
         ! it tells you how many of these atoms are in the compound
            if (i == 1) then ! it's the first symbol, it must be an element!
               write(error_message,*) 'Symbol ', trim(adjustl(Name(i:i))), ' in the compound formula could not be identified'
               INFO = 4
               goto 911
            endif
            read(Name(i:i),*) cur ! read it into a real number
            if (devide) then
               coun = coun + 1
            endif
            C = C*10.0d0 + real(cur)
         else if (verify(trim(adjustl(Name(i:i))),trim(adjustl(UpCase))) == 0) then ! it's an upper-case latine letter
         ! the element name starts

            if (num > 0) then ! at least one element was found:
                if (C <= 0) then
                   ElPersent(num) = 1  ! there is at least one
                else
                   C = C/(10.0d0**coun) ! to account for non-integer element percentage
                   ElPersent(num) = C  ! that's how many of this element in the compound
                endif
                ElNames(num) = trim(adjustl(El))
            endif
            num = num + 1   ! new element found in the compound
            ! Start a new element:
            C = 0            ! to start
            coun = 0         ! to start
            devide = .false. ! to start
            El = '   '       ! start a new name
            El = trim(adjustl(El))//Name(i:i)   ! write it's name 
         else if (verify(trim(adjustl(Name(i:i))),trim(adjustl(LowCase))) == 0) then ! it's a lower-case latine letter
         ! the element name still goes on
            if (i == 1) then ! it's the first symbol, it must be an element!
               write(error_message,*) 'Symbol ', trim(adjustl(Name(i:i))), ' in the compound formula could not be identified'
               INFO = 4
               goto 911
            endif
            El = trim(adjustl(El))//Name(i:i)   ! continue writing it's name
         else if ((trim(adjustl(Name(i:i))) == '.')) then ! decimal point
            if (i == 1) then ! it's the first symbol, it must be an element!
               write(error_message,*) 'Symbol ', trim(adjustl(Name(i:i))), ' in the compound formula could not be identified'
               INFO = 4
               goto 911
            endif
            devide =.true. 
         else ! it's another symbol
         ! no idea what that might be...
            write(error_message,*) 'Symbol ', trim(adjustl(Name(i:i))), ' in the compound formula could not be identified'
            INFO = 4
            goto 911
         endif
      enddo
      if (num > 0) then ! the last element:
        C = C/(10.0d0**coun) ! to account for non-integer element percentage
        if (C <= 0) then
           ElPersent(num) = 1  ! there is at least one
        else
           ElPersent(num) = C  ! that's how many of this element in the compound
        endif
        ElNames(num) = trim(adjustl(El))
      endif

      at_num = num ! how many different elements are in this compound
      if (present(at_numbers)) then
         if (.not.allocated(at_numbers)) then
            allocate(at_numbers(at_num))
         else
            deallocate(at_numbers)
            allocate(at_numbers(at_num))
         endif
      endif
      if (present(at_percentage)) then
         if (.not.allocated(at_percentage)) then
            allocate(at_percentage(at_num))
         else
            deallocate(at_percentage)
            allocate(at_percentage(at_num))
         endif
      endif
      if (present(at_short_names)) then
         if (.not.allocated(at_short_names)) then
            !allocate(character(3) :: at_short_names(at_num))
            allocate(at_short_names(at_num))
         else
            deallocate(at_short_names)
            !allocate(character(3) :: at_short_names(at_num))
            allocate(at_short_names(at_num))
         endif
      endif
      if (present(at_names)) then
         if (.not.allocated(at_names)) then
            !allocate(character(15) :: at_names(at_num))
            allocate(at_names(at_num))
         else
            deallocate(at_names)
            !allocate(character(15) :: at_names(at_num))
            allocate(at_names(at_num))
         endif
      endif
      if (present(at_masses)) then
         if (.not.allocated(at_masses)) then
            allocate(at_masses(at_num))
         else
            deallocate(at_masses)
            allocate(at_masses(at_num))
         endif
      endif
      if (present(at_NVB)) then
         if (.not.allocated(at_NVB)) then
            allocate(at_NVB(at_num))
         else
            deallocate(at_NVB)
            allocate(at_NVB(at_num))
         endif
      endif

      do i = 1, num
         call find_atomic_number(ElNames(i), ElNumbers(i), NVB, Full_name, M, Periodic_table, INFO, error_message)  ! below
         if (INFO /= 0) goto 911
         if (present(at_numbers)) at_numbers(i) = ElNumbers(i)  ! element number, Z
         if (present(at_percentage)) at_percentage(i) = ElPersent(i) ! percentage
         if (present(at_short_names)) then
            write(at_short_names(i),'(a)') ElNames(i) ! name of the element
         endif
         if (present(at_names)) then
            write(at_names(i),'(a)') Full_name  ! full name of the element
         endif
         if (present(at_masses)) at_masses(i) = M  ! mass of the element in the proton-mass units
         if (present(at_NVB)) at_NVB(i) = NVB      ! number of valence electrons
!          print*, 'Check: ', ElNames(i)
!          print*, 'LEN: ', LEN(at_short_names(i))
!          print*, 'Decompose_compound: ', at_short_names(i)
      enddo

   else exists
      write(error_message,*) 'Could not find the file: ', trim(adjustl(File_name))
      Reason = 1
      INFO = 1
   endif exists

911 if (allocated(Periodic_table)) deallocate(Periodic_table)
   if (Reason /= 0) then
      write(error_message,*) 'Could not read input file: ', trim(adjustl(File_name))
      INFO = 3
   endif
end subroutine Decompose_compound


subroutine find_atomic_number(Name, Z, NVB, Full_name, M, Periodic_table, INFO, error_message)
   character(*), intent(in) :: Name ! element abbreviation
   integer, intent(out) :: Z        ! atomic number
   integer, intent(out) :: NVB      ! number of valence electrons
   type(atomic_data), dimension(:), intent(in) :: Periodic_table ! this is an internal module variable
   character(15), intent(out), optional :: Full_name ! element name
   real(8), intent(out), optional  :: M ! mass
   integer, intent(inout) :: INFO  ! 0=file read well; 1=no file; 2=couldn't open; 3=error while reading
   character(100) :: error_message
   integer N, i
   N = size(Periodic_table) ! that's how many elements we have in our table
   do i = 1, N ! compare names
      if (trim(adjustl(Name)) == trim(adjustl(Periodic_table(i)%Name))) then
         Z = Periodic_table(i)%Z
         NVB = Periodic_table(i)%Nvb
         if (present(Full_name)) Full_name = Periodic_table(i)%Full_name
         if (present(M)) M = Periodic_table(i)%Mass
         exit
      endif
      if (i .GE. N) then
         write(error_message,*) 'Element ', trim(adjustl(Name)), ' was not found in the database: INPUT_atomic_data.dat'
         INFO = 4
      endif
   enddo
end subroutine find_atomic_number


subroutine Count_lines(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    integer i
    if (present(skip_lines)) then
       do i=1,skip_lines
          read(File_num,*, end=604) 
       enddo
       604 continue
    endif
    i = 0
    do
        read(File_num,*, end=603)
        i = i + 1
    enddo
    603 continue
    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    N = i
end subroutine Count_lines


END MODULE Periodic_table
