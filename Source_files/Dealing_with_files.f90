! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2016-2018
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with files:
MODULE Dealing_with_files

! USE IFPORT ! standart FORTRAN module for dealing with files

implicit none

 contains


subroutine get_file_stat(File_name, device_ID, Inode_number, File_mode, Number_of_links, O_uid, O_gid, where_located, &
                         File_size, Last_access_time, Last_modification_time, Last_status_change, blocks_allocated)
! See description here: https://software.intel.com/en-us/node/526830
   character(*), intent(in) :: File_name ! which file we are checking?
   integer, intent(out), optional :: device_ID ! Device the file resides on
   integer, intent(out), optional :: Inode_number ! File inode number
   integer, intent(out), optional :: File_mode ! Access mode of the file
   integer, intent(out), optional :: Number_of_links ! Number of hard links to the file
   integer, intent(out), optional :: O_uid ! User ID of owner
   integer, intent(out), optional :: O_gid ! Group ID of owner
   integer, intent(out), optional :: where_located ! Raw device the file resides on
   integer, intent(out), optional :: File_size ! Size of the file
   integer, intent(out), optional :: Last_access_time ! Time when the file was last accessed (*)
   integer, intent(out), optional :: Last_modification_time ! Time when the file was last modified(*)
   integer, intent(out), optional :: Last_status_change ! Time of last file status change (*)
   integer, intent(out), optional :: blocks_allocated ! Blocksize for file system I/O operations
   !(*) Times are in the same format returned by the TIME function (number of seconds since 00:00:00 Greenwich mean time, January 1, 1970).
   !=====================
   INTEGER :: info_array(12)
   
   ! Get the statistics on the file:
   call STAT(trim(adjustl(File_name)), info_array) ! intrinsic fortran subroutine
   
   if (present(device_ID)) device_ID = info_array(1)  ! Device the file resides on
   if (present(Inode_number)) Inode_number = info_array(2) ! File inode number
   if (present(File_mode)) File_mode = info_array(3) ! Access mode of the file
   if (present(Number_of_links)) Number_of_links = info_array(4) ! Number of hard links to the file
   if (present(O_uid)) O_uid = info_array(5) ! User ID of owner
   if (present(O_gid)) O_gid = info_array(6) ! Group ID of owner
   if (present(where_located)) where_located = info_array(7) ! Raw device the file resides on
   if (present(File_size)) File_size = info_array(8) ! Size of the file
   if (present(Last_access_time)) Last_access_time = info_array(9) ! Time when the file was last accessed (*)
   if (present(Last_modification_time)) Last_modification_time = info_array(10) ! Time when the file was last modified(*)
   if (present(Last_status_change)) Last_status_change = info_array(11) ! Time of last file status change (*)
   if (present(blocks_allocated)) blocks_allocated = info_array(12) ! Blocksize for file system I/O operations
end subroutine get_file_stat


subroutine open_file(FN, File_name, Error_descript, status, action)
   integer, intent(inout) :: FN	! file number
   character(*), intent(inout) :: File_name	! file name
   character(*), intent(inout) :: Error_descript	! message about an error in opening file
   character(*), intent(in), optional:: status, action	! status and action variables for the fortraen OPEN statement
   !-------------------------------
   logical :: file_exists, file_opened
   ! 1) check if file exists:
   inquire(file=trim(adjustl(File_name)),exist=file_exists)
   if (.not.file_exists) then		! print out a message for the user:
      Error_descript = 'File '//trim(adjustl(File_name))//' could not be found'
      goto 9999	! exit the subroutine with an error message
   endif
   ! 2) try to open the file:
   if (present(status)) then	! use the user-defined status:
      if (present(action)) then	! both options are present, e.g. status = 'old', action='read'
         open(UNIT=FN, FILE = trim(adjustl(File_name)), status = trim(adjustl(status)), action=trim(adjustl(action)))
      else	! only action is present, e.g. status = 'old'
         open(UNIT=FN, FILE = trim(adjustl(File_name)), status = trim(adjustl(status)))
      endif
   else
      if (present(action)) then	! only action is present, e.g. action='read'
         open(UNIT=FN, FILE = trim(adjustl(File_name)), action=trim(adjustl(action)))
      else	! no additional options are present
         open(UNIT=FN, FILE = trim(adjustl(File_name)))
      endif
   endif
   ! 3) check if file was opened successfully:
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (.not.file_opened) then	! print out a message for the user:
      Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened'
      goto 9999	! exit the subroutine with an error message
   endif
9999 continue
end subroutine open_file


subroutine close_file(safe, File_name, FN)
   character(*), intent(in) :: safe ! close or delete
   character(*), intent(in), optional :: File_name ! path to file named
   integer, intent(in), optional :: FN   ! file number to be referred to
   logical file_exist, file_opened
   integer FN1
   character(200) :: command
   
   if (present(FN)) then	! if the user specified the file number, use it:
      inquire(UNIT = FN, opened=file_opened)
      FN1 = FN
   endif
   if (present(File_name)) then	! if the user specified the file name, use it instead:
      inquire(file=trim(adjustl(File_name)),opened=file_opened, number=FN1)
   endif
   if (file_opened) then	! If the defined file is still opened:
      select case (trim(adjustl(safe)))	! if the user specified special options for closing the file, use them:
      case ('delete')	! delete the file if not needed:
         close(FN1, status=trim(adjustl(safe)))
      case default	! just close it as is:
         close(FN1)
      end select
   else
      select case (trim(adjustl(safe)))	! if the user specified special options for closing the file, use them:
      case ('delete')	! delete the file if not needed:
        if (present(File_name)) then	! if the user specified the file name, use it instead:
            FN1 = 200
            open(UNIT = FN1, FILE = trim(adjustl(File_name)))   ! reopen it to close properly
            inquire(file=trim(adjustl(File_name)),opened=file_opened, number=FN1)
            close(FN1, status=trim(adjustl(safe)))  ! close and delete
         endif
       end select
   endif
end subroutine close_file



subroutine read_file(Reason, i, read_well, Error_descript)
   integer, intent(in) :: Reason    ! description of error
   integer, intent(inout) :: i      ! line number
   logical, intent(inout) :: read_well  ! did we read ok?
   character(*), intent(out), optional :: Error_descript ! what is wrong
   character(16) chi
   i = i + 1    ! it's next line
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       if (present(Error_descript)) then
          write(chi, '(i12)') i
          write(Error_descript,'(a,a,a)') 'Problem reading input file in line ', trim(adjustl(chi)), ' ,wrong type of variable'
       endif
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       if (present(Error_descript)) then
          write(chi, '(i12)') i
          write(Error_descript,'(a,a,a)') 'Problem reading input file in line ', trim(adjustl(chi)), ' ,unexpected END of file'
       endif
       read_well = .false.
   ELSE   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   END IF
end subroutine read_file




subroutine Count_columns_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of columns in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    real(8) temp
    character(1000) temp_ch
    integer i, Reason
    integer :: temp_i
    if (present(skip_lines)) then
       do i=1,skip_lines
          read(File_num,*, end=601) 
       enddo
       601 continue
    endif

    read(File_num,'(a)', IOSTAT=Reason) temp_ch ! count columns in this line
    N = number_of_columns(trim(adjustl(temp_ch))) ! see below

    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
end subroutine Count_columns_in_file


pure function number_of_columns(line)
   integer :: number_of_columns
   character(*), intent(in) :: line
   integer i, n, toks
   logical :: same_space
   same_space = .false.
   i = 0
   n = len(line)
   number_of_columns = 0
   do while(i < n) ! scan through all the line
      i = i + 1
      selectcase (line(i:I))
      case (' ', '	') ! space or tab can be a separator between the columns
         if (.not.same_space) number_of_columns = number_of_columns + 1
         same_space = .true. ! in case columns are separated by more than one space or tab
      case default ! column data themselves, not a space inbetween
         same_space = .false.
      endselect
   enddo
   number_of_columns = number_of_columns + 1	! number of columns is by 1 more than number of spaces inbetween
end function number_of_columns 


subroutine Count_columns_in_file_OLD(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of columns in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    real(8) temp
    integer i, Reason
    if (present(skip_lines)) then	! If the user specified to start counting columns not from the first one:
       do i=1,skip_lines
          read(File_num,*, end=601) 
       enddo
       601 continue
    endif
    i = 0
    do	! count all the columns in the file from the specified one to the end:
        read(File_num, '(e25.16)', advance='no', IOSTAT=Reason) temp
        if (Reason .NE. 0) exit	! found the end of file
        i = i + 1
    enddo
    602 continue
    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    N = i	! print out the number of columns from the specified one to the last one
end subroutine Count_columns_in_file_OLD


subroutine Count_lines_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    integer i
    if (present(skip_lines)) then	! If the user specified to start counting lines not from the first one:
       do i=1,skip_lines
          read(File_num,*, end=604) 
       enddo
       604 continue
    endif
    i = 0
    do	! count all the lines in the file from the specified one to the end:
        read(File_num,*, end=603)
        i = i + 1
    enddo
    603 continue
    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    N = i	! print out the number of lines from the specified one to the last one
end subroutine Count_lines_in_file



subroutine copy_file(file_to_copy, folder_copy_to, OS_ind, add_com)
   character(len=*), intent(in) :: file_to_copy, folder_copy_to
   integer, intent(in), optional :: OS_ind ! windows or linux
   character(len=*), intent(in), optional :: add_com  ! additional options provided (such as /Y /F)
   character(200) command, add_option

   if (present(add_com)) then
      add_option = add_com
   else
      add_option = ''   ! no additional options
   endif

   if (present(OS_ind)) then
      select case (OS_ind)
      case (0) ! linux
         command='/cp '//trim(adjustl(file_to_copy))//' '//trim(adjustl(folder_copy_to))//' '//trim(adjustl(add_option))
      case default ! assume windows
         command='xcopy '//trim(adjustl(file_to_copy))//' '//trim(adjustl(folder_copy_to))//' '//trim(adjustl(add_option))
      end select
   else ! assume linux
      command='/cp '//trim(adjustl(file_to_copy))//' '//trim(adjustl(folder_copy_to))//' '//trim(adjustl(add_option))
   endif
   CALL system(command)
end subroutine copy_file


subroutine unix_unlimiting()	! set unlimited memory in UNIX systems, and set specific paths to MKL library (environment  specific routine!)
   character(200) command
   command='ulimit -s unlimited'
   CALL system(command)
   command='limit stacksize unlimited'
   CALL system(command)
   command='export LD_LIBRARY_PATH=/opt/intel/2011/lib/intel64:$LD_LIBRARY_PATH'
   CALL system(command)
   command='export LD_LIBRARY_PATH=/opt/products/mkl/11.0/mkl/lib/em64t:$LD_LIBRARY_PATH'
   CALL system(command)
end subroutine unix_unlimiting


subroutine Path_separator(path_sep)	! fine which path separator is used in your environment, it defines your OS 	(Windows vs Linux)
   CHARACTER(len=1), intent(out) :: path_sep
   CHARACTER(len = 100) :: path
   CALL get_environment_variable("PATH",path)
   if (path(1:1) .EQ. '/') then	!unix based OS
       path_sep = '/'
   else if (path(3:3) .EQ. '\') then	!Windows OS
       path_sep = '\'
   else
       path_sep = ''
       print*, 'Path separator is not defined'	!Unknown OS
   endif 
end subroutine Path_separator


END MODULE Dealing_with_files
