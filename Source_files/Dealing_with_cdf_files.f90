! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2014-2018
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with files:
MODULE Dealing_with_cdf_files
use Universal_constants
use Objects
use Periodic_table
use Dealing_with_files
use Dealing_with_EADL

implicit none

 contains


subroutine read_cdf_file(FN, File_name, Material, Numpar, Err)
   integer, intent(in) :: FN	! file number to read cdf from
   character(*), intent(in) :: File_name	! name of the file
   type(Target_atoms), intent(inout) :: Material	!material parameters of this one target of all
   type(Num_par) :: Numpar	! numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !---------------
   integer :: i, j, k, l, Reason, N, INFO, N_sh, CDF_coef, N_sh_num, PQN
   character(200) :: Path, temp, error_message	! where to find periodic table files
   character(30) :: Full_Name, Shell_name
   integer, dimension(:), allocatable :: at_numbers
   real(8), dimension(:), allocatable :: at_percentage
   character(3), dimension(:), allocatable :: at_short_names ! name of the element
   real(8), dimension(:), allocatable :: at_masses ! mass of each element [Mp]
   integer, dimension(:), allocatable :: at_NVB ! number of valence electrons
   logical :: read_well
   i = 0	! start counting lines
   
   Path = Numpar%input_path//'Atomic_parameters'	! where to find pariodic table file
   
   READ(FN,'(a100)',IOSTAT=Reason) Material%Name ! first line is the full material name of this target
   
   READ(FN,*,IOSTAT=Reason) N   ! number of elements in this compound / or name of the compound in "short cdf" format
   
   SHRT:IF (Reason .GT. 0)  THEN	! if it's a formula:
      backspace(FN) ! read this line again:
      READ(FN,*,IOSTAT=Reason) temp   ! chemical formula of this compound
      call read_file(Reason, i, read_well) ! reports if everything read well
      if (.not. read_well) then
         write(error_message,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', i
         call Save_error_details(Err, 2, error_message)	! module "Objects"
         goto 2014
      endif
      ! Extract the element names from the string passed, and get all the parameters from the periodic table for all elements:
      call Decompose_compound(Path, Material%Name, Numpar%path_sep, INFO, error_message, Material%N_Elements, at_numbers, at_percentage, at_short_names, at_masses=at_masses, at_NVB=at_NVB) ! module "Periodic_table"
      if (INFO == 1) then	!  INFO: 0=file read well; 1=no file; 2=couldn't open; 3=error while reading
         write(error_message,'(a)') 'File '//trim(adjustl(Path))//Numpar%path_sep//'INPUT_atomic_data.dat  not found'
         call Save_error_details(Err, 1, error_message)	! module "Objects"
         goto 2014
      else if (INFO == 2) then
         write(error_message,'(a)') 'File '//trim(adjustl(Path))//Numpar%path_sep//'INPUT_atomic_data.dat  could not be opened'
         call Save_error_details(Err, 1, error_message)	! module "Objects"
         goto 2014
      else if (INFO == 3) then
         write(error_message,'(a)') 'File '//trim(adjustl(Path))//Numpar%path_sep//'INPUT_atomic_data.dat  could not be read'
         call Save_error_details(Err, 1, error_message)	! module "Objects"
         goto 2014
      endif
      N = Material%N_Elements	! that's how many kinds of atoms we have in this target
      if (.not.allocated(Material%Elements)) allocate(Material%Elements(N))	! allocate objects with parameters
      ! Save all the parameters into those objects:
      do j = 1, N  ! read for each element it's basic data:
         Material%Elements(j)%Name = at_short_names(j)	! chemical element name
         Material%Elements(j)%Zat = at_numbers(j)	 ! atomic number [electron charge]
         Material%Elements(j)%Mass = at_masses(j)	! [a.m.u]
         Material%Elements(j)%M = Material%Elements(j)%Mass*g_amu		! [kg]
         Material%Elements(j)%percentage = at_percentage(j)	! contribution to the compound
         Material%Elements(j)%NVB = at_NVB(j)	! number of valence electrons
      enddo
   else SHRT	! if it is a number:
      call read_file(Reason, i, read_well) ! reports if everything read well
      if (.not. read_well) then
         write(error_message,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', i
         call Save_error_details(Err, 2, error_message)	! module "Objects"
         goto 2014
      endif
      if (.not. allocated(Material%Elements)) allocate(Material%Elements(N)) ! that's how many atom kinds we have
      do j = 1, N  ! read for each element it's basic data:
         READ(FN,*,IOSTAT=Reason) Material%Elements(j)%Zat, Material%Elements(j)%percentage ! atomic number and its persentage in the compound's molecule
         call read_file(Reason, i, read_well) ! reports if everything read well
         if (.not. read_well) goto 2014
!          call Find_element_name(INT(Material%Elements(j)%Zat), Material%Elements(j)%Name, Full_Name, Material%Elements(j)%Mass) ! from module 'Dealing_with_EADL'
         Material%Elements(j)%M = Material%Elements(j)%Mass*g_amu		! [kg]
      enddo
   endif SHRT

   READ(FN,*,IOSTAT=Reason) Material%Dens, Material%v_sound, Material%DOS%E_f   ! material density [g/cm^3], speed of sound in the material [m/s], and Fermi energy [eV]
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) then
      write(error_message,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', i
      call Save_error_details(Err, 2, error_message)	! module "Objects"
      goto 2014
   endif
   Material%At_Dens = 1.0d-3*Material%Dens/(SUM(Material%Elements(:)%M)/dble(N))	! [1/cm^3] atomic density
   Material%DOS%v_f = sqrt(2.0d0*Material%DOS%E_f/g_me)	! units as used in one-pole approximation (not SI!)
  
   do j = 1, N	! read for each element its shells data:
      READ(FN,*,IOSTAT=Reason) N_sh	! number of shells in the first atom kind:
      call read_file(Reason, i, read_well) ! reports if everything read well
      if (.not. read_well) then
         write(error_message,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', i
         call Save_error_details(Err, 2, error_message)	! module "Objects"
         goto 2014
      endif
   
      ! Now we know how many shells are in this atom:
      Material%Elements(j)%N_shl = N_sh   ! save it
      if (.not. allocated(Material%Elements(j)%Shell_name)) allocate(Material%Elements(j)%Shell_name(N_sh)) ! allocate shell-names for each shell
      if (.not. allocated(Material%Elements(j)%Shl_dsgnr)) allocate(Material%Elements(j)%Shl_dsgnr(N_sh)) ! allocate shell disignator for each shell
      if (.not. allocated(Material%Elements(j)%Ne_shell)) then
         allocate(Material%Elements(j)%Ne_shell(N_sh)) ! allocate numbers of electrons for each shell
         Material%Elements(j)%Ne_shell(:) = 0.0d0
      endif
      if (.not. allocated(Material%Elements(j)%Ip)) then
         allocate(Material%Elements(j)%Ip(N_sh)) ! allocate ionization potentials for each shell
         Material%Elements(j)%Ip(:) = -1.0d15
      endif
      if (.not. allocated(Material%Elements(j)%Ek)) then
         allocate(Material%Elements(j)%Ek(N_sh)) ! allocate mean kinetic energies for each shell
         Material%Elements(j)%Ek(:) = -1.0d-15
      endif
      if (.not. allocated(Material%Elements(j)%Auger)) then
         allocate(Material%Elements(j)%Auger(N_sh)) ! allocate auger-times for each shell
         Material%Elements(j)%Auger(:) = 1.0d31
      endif
      if (.not. allocated(Material%Elements(j)%Radiat)) then
         allocate(Material%Elements(j)%Radiat(N_sh)) ! allocate radiative-times for each shell
         Material%Elements(j)%Radiat(:) = 2.0d31
      endif
!       if (.not. allocated(Material%Elements(j)%PQN)) then
!          allocate(Material%Elements(j)%PQN(N_sh)) ! allocate principle quantum numbers for each shell
!          Material%Elements(j)%PQN(:) = 0
!       endif
!       if (.not. allocated(Material%Elements(j)%KOCS)) then
!          allocate(Material%Elements(j)%KOCS(N_sh)) ! allocate kind of inelastic cross sections
!          Material%Elements(j)%KOCS(:) = 0
!       endif
!       if (.not. allocated(Material%Elements(j)%KOCS_SHI)) then
!          allocate(Material%Elements(j)%KOCS_SHI(N_sh)) ! allocate kind of inelastic cross sections
!          Material%Elements(j)%KOCS_SHI(:) = 0
!       endif
      if (.not. allocated(Material%Elements(j)%CDF)) then
         allocate(Material%Elements(j)%CDF(N_sh)) ! allocate CDF-functions' coefficiants for each shell
      endif

      do k = 1, N_sh ! read all the data for each shell:
          READ(FN,*,IOSTAT=Reason) CDF_coef, N_sh_num, Material%Elements(j)%Ip(k), Material%Elements(j)%Ne_shell(k), Material%Elements(j)%Auger(k)	! number of CDF functions, shell-designator, ionization potential, number of electrons, Auger-time
          call read_file(Reason, i, read_well)	! reports if everything read well
          if (.not. read_well) then
             write(error_message,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', i
             call Save_error_details(Err, 2, error_message)	! module "Objects"
             goto 2014
          endif
          Material%Elements(j)%Shl_dsgnr(k) = ABS(N_sh_num)  ! save shell designator
!           call define_PQN(Material%Elements(j)%Shl_dsgnr(k), Shell_name, PQN)	! from module 'Dealing_with_EADL'
          if (Material%Elements(j)%Shl_dsgnr(k) .GE. 63) Material%N_VB_el = Material%Elements(j)%Ne_shell(k)		! then it is the valence band, save number of VB electrons
          Material%Elements(j)%Shell_name(k) = Shell_name	! save the name of this shell

          ! If some data are missing in the input file, get it from EADL database:
!           call check_atomic_parameters(NumPar, Material%Elements, j, k, N_sh, Err, read_well)	! see below

          if (Material%Elements(j)%Ip(k) .LT. 1.0d-1) Material%Elements(j)%Ip(k) = 1.0d-1	! [eV], introduce at least a minimum "band gap"
          if (CDF_coef .GT. 0) then	! there is CDF to be used in a cross section
!              Material%Elements(j)%KOCS_SHI(k) = 1	! CDF cross section for this shell
!              if (N_sh_num .GT. 0) then	! use normal CDF cross sections for electrons:
!                 Material%Elements(j)%KOCS(k) = 1	! CDF cross section for this shell
!              else ! use BEB for electrons (but still CDF for the SHI):
!                 Material%Elements(j)%KOCS(k) = 2	! BEB cross section for this shell
!              endif
             if (.not. allocated(Material%Elements(j)%CDF(k)%E0)) allocate(Material%Elements(j)%CDF(k)%E0(CDF_coef))   ! that's how many CDF function for this shell of this atom
             if (.not. allocated(Material%Elements(j)%CDF(k)%A)) allocate(Material%Elements(j)%CDF(k)%A(CDF_coef))   ! that's how many CDF function for this shell of this atom
             if (.not. allocated(Material%Elements(j)%CDF(k)%Gamma)) allocate(Material%Elements(j)%CDF(k)%Gamma(CDF_coef))   ! that's how many CDF function for this shell of this atom
             do l = 1, CDF_coef    ! read all the CDF parameters:
                READ(FN,*,IOSTAT=Reason) Material%Elements(j)%CDF(k)%E0(l), Material%Elements(j)%CDF(k)%A(l), Material%Elements(j)%CDF(k)%Gamma(l)
                call read_file(Reason, i, read_well) ! reports if everything read well
                if (.not. read_well) then
                   write(error_message,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', i
                   call Save_error_details(Err, 2, error_message)	! module "Objects"
                   goto 2014
                endif
             enddo
          else ! there is no CDF, so let's use atomic BEB cross section:
!              Material%Elements(j)%KOCS(k) = 2 ! BEB cross section for this shell
!              Material%Elements(j)%KOCS_SHI(k) = 2 ! BEB cross section for this shell
          endif
      enddo
   enddo
   
   READ(FN,*,IOSTAT=Reason) CDF_coef
   i = i + 1
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       write(*,'(a,i3,a)') 'Line ', i, ', END of file; no CDF data for phonon-peak' 
       write(*,'(a,i3,a)') 'Using atomic cross-sections for elastic scattering'
       read_well = .true.
   ELSE   ! normal reading
       read_well = .true.  ! it read well, nothing to report
       if (.not. allocated(Material%CDF_phonon%A)) allocate(Material%CDF_phonon%A(CDF_coef))
       if (.not. allocated(Material%CDF_phonon%E0)) allocate(Material%CDF_phonon%E0(CDF_coef))
       if (.not. allocated(Material%CDF_phonon%Gamma)) allocate(Material%CDF_phonon%Gamma(CDF_coef))
       do l = 1, CDF_coef    ! read all the CDF parameters for phonon peak:
          READ(FN,*,IOSTAT=Reason) Material%CDF_phonon%E0(l), Material%CDF_phonon%A(l), Material%CDF_phonon%Gamma(l)
          call read_file(Reason, i, read_well) ! reports if everything read well
          if (.not. read_well) goto 2014
       enddo
   END IF

   
2014 continue
   ! Clean up at the end:
   if (allocated(at_numbers)) deallocate(at_numbers)
   if (allocated(at_percentage)) deallocate(at_percentage)
   if (allocated(at_short_names)) deallocate(at_short_names)
   if (allocated(at_masses)) deallocate(at_masses)
   if (allocated(at_NVB)) deallocate(at_NVB)
end subroutine read_cdf_file



! 
! subroutine check_atomic_parameters(NumPar, Target_atoms, N_at, cur_shl, shl, Error_message, read_well) ! from module 'Dealing_with_EADL'
!    type(Num_par), intent(inout) :: NumPar ! numerical parameters
!    type(Atom_kind), dimension(:), allocatable, intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
!    integer, intent(in), optional :: N_at
!    integer, intent(in), optional :: cur_shl, shl  ! current atoms, shell, and total number of shells
!    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
!    logical, intent(inout) :: read_well  ! did we read the file well?
!    
!    integer INFO, i, j, k, FN, FN1, Z, Shl_dsgtr, Nat
!    character(100) :: File_name_EADL, Folder_name, Error_descript
!    character(100) :: File_name_EPDL97, File_name, File_name2
!    logical :: File_opened
!    
!    ! check if these database-files exist:
!    call Check_EADP_EPDL97(NumPar%path_sep, File_name_EADL, File_name_EPDL97, INFO) ! see below
!    if (INFO .EQ. 0) then ! files are there, reading them:
!       if (present(cur_shl)) then ! for this shell only:
!          Z = Target_atoms(N_at)%Zat ! atomic number
! 
!          if (Target_atoms(N_at)%Ne_shell(cur_shl) .LE. 0) then ! nondefined in the input file, find the default value from EADL:
! !             call READ_EADL_TYPE_FILE_int(FN1, File_name_EADL, Z, 912, Target_atoms, N_at, cur_shl, shl) ! read PQNs and shell-names
!             call READ_EADL_TYPE_FILE_int(FN1, File_name_EADL, Z, 912, INFO, Error_descript, N_shl=Target_atoms(N_at)%N_shl, Nel=Target_atoms(N_at)%Ne_shell, Shell_name=Target_atoms(N_at)%Shell_name, Shl_num=Target_atoms(N_at)%Shl_dsgnr)
!          endif
!          
!          if (Target_atoms(N_at)%Ip(cur_shl) .LE. -1.0d-14) then ! nondefined in the input file, find the default value from EADL:
! !             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 913, Target_atoms(N_at)%Ip, cur_shl, shl, Target_atoms(N_at)%Shl_num(cur_shl)) ! read binding energies
! 
!             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 913, Target_atoms(N_at)%Ip, cur_shl, INFO=INFO, error_message=Error_descript)
!          endif
!          
!          if (Target_atoms(N_at)%Ek(cur_shl) .LE. 0.0d0) then ! nondefined in the input file, find the default value from EADL:
!             if (Target_atoms(N_at)%Shl_dsgnr(cur_shl) .GE. 62) then ! find approximate Ekin for the VB:
!                call next_designator(Shl_dsgtr, Target_atoms(N_at)%Shl_dsgnr(cur_shl-1)) ! find the designator for the VB (next shell after last one)
!             else
!                Shl_dsgtr = Target_atoms(N_at)%Shl_dsgnr(cur_shl)
!             endif
!             !print*, Target_atoms(N_at)%Ek(cur_shl), Target_atoms(N_at)%Ip(cur_shl), Shl_dsgtr
! !             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 914, Target_atoms(N_at)%Ek, cur_shl, shl, Shl_dsgtr) ! read kinetic energies
!             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 914, Target_atoms(N_at)%Ek, cur_shl, INFO=INFO, error_message=Error_descript)
!          endif
!          
!          if (Target_atoms(N_at)%Shl_dsgnr(cur_shl) .GE. 63) then
!             Target_atoms(N_at)%Radiat(cur_shl)=1d23 ! [fs] not possible process (or not included) => infinite time
!          else
! !             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 921, Target_atoms(N_at)%Radiat, cur_shl, shl, Target_atoms(N_at)%Shl_dsgnr(cur_shl)) ! read radiative decay-times
!             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 921, Target_atoms(N_at)%Radiat, cur_shl, INFO=INFO, error_message=Error_descript)
!             ! Convert into fs:
!             Target_atoms(N_at)%Radiat(cur_shl)=1d15*g_h/(g_e*Target_atoms(N_at)%Radiat(cur_shl)) ! [fs]
! !             print*, N_at, cur_shl, Target_atoms(N_at)%Radiat(cur_shl)
!          endif
!          
!          if (Target_atoms(N_at)%Shl_dsgnr(cur_shl) .GE. 63) then
!             Target_atoms(N_at)%Auger(cur_shl)=1.0d23 ! [fs] not possible process => infinite time
!          else if ((Target_atoms(N_at)%Auger(cur_shl) .LE. 0.0d0) .OR. (Target_atoms(N_at)%Auger(cur_shl) .GT. 1.0d30)) then
! !             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 922, Target_atoms(N_at)%Auger, cur_shl, shl, Target_atoms(N_at)%Shl_dsgnr(cur_shl)) ! read auger-times
!             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 922, Target_atoms(N_at)%Auger, cur_shl, INFO=INFO, error_message=Error_descript) ! read auger-times
!             ! Convert into fs:
!             Target_atoms(N_at)%Auger(cur_shl)=1d15*g_h/(g_e*Target_atoms(N_at)%Auger(cur_shl)) ! [fs]
!          endif
!       else ! for all shells of this atom:
!          Nat = size(Target_atoms)
!          do i = 1, Nat
!             Z = Target_atoms(i)%Zat ! atomic number
!             !call READ_EADL_TYPE_FILE_int (FN1, File_name_EADL, Z, 912, Target_atoms, i) ! read PQNs and shell-names
!              call READ_EADL_TYPE_FILE_int(FN1, File_name_EADL, Z, 912, INFO, Error_descript, N_shl=Target_atoms(N_at)%N_shl, Nel=Target_atoms(N_at)%Ne_shell, Shell_name=Target_atoms(N_at)%Shell_name, Shl_num=Target_atoms(N_at)%Shl_dsgnr)
!             
!             !call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 913, Target_atoms(i)%Ip) ! read binding energies
!             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 913, Target_atoms(N_at)%Ip, INFO=INFO, error_message=Error_descript)
!             
!             !call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 914, Target_atoms(i)%Ek) ! read kinetic energies
!              call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 914, Target_atoms(N_at)%Ek, INFO=INFO, error_message=Error_descript)
!             
! !             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 921, Target_atoms(i)%Radiat) ! read radiative decay-times
!             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 921, Target_atoms(N_at)%Radiat, INFO=INFO, error_message=Error_descript)
!             
! !             call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 922, Target_atoms(i)%Auger) ! read auger-times
!              call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 922, Target_atoms(N_at)%Auger, INFO=INFO, error_message=Error_descript) ! read auger-times
! 
!             !print*, 'Atom:',i, Nat, size(Target_atoms(i)%Ip), size(Target_atoms(i)%Ek), size(Target_atoms(i)%Radiat), size(Target_atoms(i)%Auger)
!             do j = 1,size(Target_atoms(i)%Ip)
!                !print*, 'Shell',j, size(Target_atoms(i)%Ip)
!                Target_atoms(i)%Auger(j)=1d15*g_h/(g_e*Target_atoms(i)%Auger(j)) ! [fs]
!                if (Target_atoms(i)%Shl_dsgnr(j) .GE. 63) then
!                   Target_atoms(i)%Radiat(j)=1d23 ! [fs] not possible process (or not included) => infinite time
!                else
!                   Target_atoms(i)%Radiat(j)=1d15*g_h/(g_e*Target_atoms(i)%Radiat(j)) ! [fs]
!                endif
!             enddo ! j
!          enddo ! i
!       endif
!    else ! couldn't find databases:
!       Folder_name = 'INPUT_EADL'  ! here we keep databases
!       if (INFO .EQ. 1) then
!          File_name = trim(adjustl(Folder_name))//trim(adjustl(NumPar%path_sep))//'eadl.all'
!          Error_descript = 'File '//trim(adjustl(File_name))//' is not found!'    ! description of an error
!          call Save_error_details(Error_message, 21, Error_descript) ! write it into the error-log file
!          print*, trim(adjustl(Error_descript)) ! print it also on the sreen
!       endif
!       if (INFO .EQ. 2) then
!          File_name2 = trim(adjustl(Folder_name))//trim(adjustl(NumPar%path_sep))//'epdl97.all'
!          Error_descript = 'File '//trim(adjustl(File_name2))//' is not found!'    ! description of an error
!          call Save_error_details(Error_message, 22, Error_descript) ! write it into the error-log file
!          print*, trim(adjustl(Error_descript)) ! print it also on the sreen
!       endif
!       read_well = .false.   ! it didn't read well the input file...
!    endif
!    inquire(file=trim(adjustl(File_name_EADL)),opened=File_opened)
!    if (File_opened) close(FN1)   ! close the EADL database
! end subroutine check_atomic_parameters




! subroutine Check_EADP_EPDL97(path_sep, File_name, File_name2, INFO)
!     character(1), intent(in) :: path_sep
!     character(100), intent(inout) :: File_name
!     character(100), intent(inout) :: File_name2
!     integer INFO
!     character(100) :: Folder_name
!     logical file_exist
!     Folder_name = 'INPUT_EADL'  ! here we keep databases
!     File_name = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//'eadl.all'
!     File_name2 = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//'epdl97.all'
!     inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
!     if (.not. file_exist) then
!        write(*,'(a,a,a)') 'File ',trim(adjustl(File_name)),' is not found.'
!        INFO = 1
!        !print*, 'The program cannot continue.'
!        !goto 9999
!     else
!        INFO = 0
!     endif
!     inquire(file=trim(adjustl(File_name2)),exist=file_exist) ! check if input file is there
!     if (.not. file_exist) then
!        write(*,'(a,a,a)') 'File ',trim(adjustl(File_name2)),' is not found.'
!        !print*, 'The program cannot continue.'
!        !goto 9999
!        INFO = 2
!     endif
! end subroutine Check_EADP_EPDL97
 

END MODULE Dealing_with_cdf_files
