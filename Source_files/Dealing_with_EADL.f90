! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2014-2018
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines for dealing with EADL and EPDL databases
! The format description is taken from:
! "ENDL Type Formats for the LLNL Evaluated Atomic Data Library (EADL), Evaluated Electron Data Library (EEDL), and Evaluated Photon Data Library (EPDL)"
! by S. T. Perkins and D. E. Cullen, 2002
! 1111111111111111111111111111111111111111111111111111111111111
MODULE Dealing_with_EADL

use Universal_constants
use Little_subroutines

implicit none

!==============================================
! For reading atomic data from our periodic table:
type Shell_designator    ! our internal database "ENDL_shell_designators.dat"
   integer :: Shl	! shell designator
   character(7) :: Shell_spect	! spectroscopic shell name
   character(6) :: Shell_at	! atomic shell name
endtype Shell_designator
!==============================================

 contains



subroutine Read_EPDL_rata(path_sep, FN, File_name, INFO, Z, Shl_dsgnr, Phot_abs_CS_tot, Phot_abs_CS_shl)
   character(1), intent(in) :: path_sep	! path separator
   integer, intent(in) :: FN	! file number of ENDL database
   character(100), intent(in) :: File_name	! file name of ENDL database
   integer, intent(inout) :: INFO	! info whether file read well
   integer, intent(in) :: Z	! atomic number
   integer, dimension(:), allocatable, intent(in) :: Shl_dsgnr	! array of EADL shell designators
   real(8), dimension(:,:), allocatable, intent(out) :: Phot_abs_CS_tot	! Total photoabsorption cross section:  Photon energy [eV], Cross section [A^2]
   real(8), dimension(:,:,:), allocatable, intent(out) :: Phot_abs_CS_shl	! per shell cross sections of photoabsorption (to be extracted from EPDL): shell, Cross section [A^2]
   !-----------------------------
   integer :: EndCheck, i, counter, j, Nsiz, shl_num_1, Nshl
   integer :: Z0, A0, Yi0, Yo0, Date0, C0, I0, S0, X10, Iflag0
   real(8) :: AW0, temp
   
   Nshl = size(Shl_dsgnr)	! number of shells in this atom
   
   ! Start by reading the header of the first block
   call Read_ENDL_header(FN, File_name, INFO, Z0, A0, Yi0, Yo0, AW0, Date0, C0, I0, S0, X10, Iflag0)	! see below
   if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine

   EndCheck = 0	! to start with
   do while (Z0 < Z)	! skip lines until you find the element you need
      do while (EndCheck /= 1)	! until end of the block
         read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
         if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
      enddo
      EndCheck = 0	! restart checker
      ! Read the header of the next block:
      call Read_ENDL_header(FN, File_name, INFO, Z0, A0, Yi0, Yo0, AW0, Date0, C0, I0, S0, X10, Iflag0)	! see below
      if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
!       print*, 'EPDL:', Z0
   enddo
   
   ! We found the element we need:
   Z0Z:do while (Z0 == Z)	! read all the data for this element
      ! Read data from the blocks of this elements:
      EndCheck = 0
      do while (EndCheck /= 1)	! read the entire block
         read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
         if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
         if (EndCheck /= 1) then	! if it is not the end of a block, read the line
            backspace(FN)	! get back to read this line again
            select case (C0)	! which process is it?
            case (71)	! Coherent scattering
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case (72)	! Incoherent scattering
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case (73)	! Photoabsorption
               select case (I0)	! which kind of data is it?
               case (0)	! cross section
                  select case (S0)
                  case (0)	! total cross section
!                      print*, 'Iflag0 tot', Iflag0
                     ! Allocate array if needed:
                     ! First, count how many data points we have here:
                     if (.not. allocated(Phot_abs_CS_tot)) then
                        counter = 0
                        do while (EndCheck /= 1)	! until end of the block
                           counter = counter + 1		! count lines
                           read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
                           if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
                        enddo
                        ! Then, rewind the data back to read:
                        do j = 1, counter	! rewind back to the start of the block
                           backspace(FN)	! get back to read this line again
                        enddo
                        ! Finaly, allocate arrays to the given number:
                        Nsiz = counter-1
                        allocate(Phot_abs_CS_tot(2,Nsiz))
                        Phot_abs_CS_tot = 0.0d0
                        allocate(Phot_abs_CS_shl(Nshl,2,Nsiz))
                        Phot_abs_CS_shl(:,1,:) = 1d24
                        Phot_abs_CS_shl(:,2,:) = -1.0d-12
                     else	! in case we know the size, for some reason...
                         Nsiz = size(Phot_abs_CS_tot,2)
                     endif
                     
                     EndCheck = 0
                     counter = 0
                     ! Now read the total cross section into this array:
                     do while (EndCheck /= 1)	! until end of the block
                        counter = counter + 1
                        read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
                        if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
                        if (EndCheck == 1) then
                           exit	! done with this block
                        else
                           backspace(FN)
                           read(FN,*,IOSTAT=INFO) Phot_abs_CS_tot(1,counter), Phot_abs_CS_tot(2,counter)
                           if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
                        endif
                     enddo ! while (EndCheck /= 1)
                  case (91)	! cross section per shell
!                      print*, 'Iflag0 shl', Iflag0
                     counter = 0
                     do while (EndCheck /= 1)	! until end of the block
                        counter = counter + 1		! count lines
                        read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
                        if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
                        if (EndCheck == 1) then
                           exit	! done with this block
                        else
                           backspace(FN)
                           call find_shell_by_designator(X10, Shl_dsgnr, shl_num_1)
                           read(FN,*,IOSTAT=INFO) Phot_abs_CS_shl(shl_num_1,1,counter), Phot_abs_CS_shl(shl_num_1,2,counter)
                           if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
                        endif
                     enddo 
                  case default ! skip line
                     read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
                  endselect
               case (10)	! average energy of emitted particle
                  read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
               case (11)	! average energy of atom
                  read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
               case default ! skip line
                  read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
               end select
            case (74)	! Pair production
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case (75)	! Triplet production
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case default ! skip line
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            end select
         endif	! it is the end of the block, go to the next iteration of the cycle
         
      enddo ! while (EndCheck /= 1)
      ! After the end of the block, read the header of the next one:
      call Read_ENDL_header(FN, File_name, INFO, Z0, A0, Yi0, Yo0, AW0, Date0, C0, I0, S0, X10, Iflag0)	! see below
      if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
   enddo Z0Z ! (Z0 == Z)
      
9998 continue
9997   format(71X,I1)			! block separator
end subroutine Read_EPDL_rata



 
subroutine Read_EADL_data(path_sep, FN, File_name, INFO, Z, Mass, Shl_dsgnr, Shell_name, Ip, Ek, Ne_shell, Radiat, f_rad, Auger, f_auger)
   character(1), intent(in) :: path_sep	! path separator
   integer, intent(in) :: FN	! file number of ENDL database
   character(100), intent(in) :: File_name	! file name of ENDL database
   integer, intent(inout) :: INFO	! info whether file read well
   integer, intent(in) :: Z	! atomic number
   integer, intent(inout), optional :: Mass	! atomic mass
   integer, dimension(:), allocatable, intent(inout), optional :: Shl_dsgnr	! EADL shell designator
   character(11), dimension(:), allocatable, intent(inout), optional :: Shell_name	! names of the shells
   real(8), dimension(:), allocatable, intent(inout), optional :: Ip		! [eV] ionization potentials for all shells
   real(8), dimension(:), allocatable, intent(inout), optional :: Ek		! [eV] mean kinetic energy of all shells
   real(8), dimension(:), allocatable, intent(inout), optional :: Ne_shell	! number of electron in each shell
   real(8), dimension(:), allocatable, intent(inout), optional :: Radiat	! [fs] radiative-decay times for all shells
   real(8), dimension(:,:), allocatable, intent(inout), optional :: f_rad	! probability of radiative decay from a given shell
   real(8), dimension(:), allocatable, intent(inout), optional :: Auger	! [fs] Auger-decay times for all shells
   real(8), dimension(:,:,:), allocatable, intent(inout), optional :: f_auger	! probability of Auger decay between two given shells
   !------------------------------------------------------------------------
   integer, dimension(:), allocatable :: Shl_dsgnr0	! EADL shell designator
   character(11), dimension(:), allocatable :: Shell_name0	! names of the shells
   real(8), dimension(:), allocatable :: Ip0		! [eV] ionization potentials for all shells
   real(8), dimension(:), allocatable :: Ek0		! [eV] mean kinetic energy of all shells
   real(8), dimension(:), allocatable :: Ne_shell0	! number of electron in each shell
   real(8), dimension(:), allocatable :: Radiat0	! [fs] radiative-decay times for all shells
   real(8), dimension(:,:), allocatable :: f_rad0	! probability of radiative decay from a given shell
   real(8), dimension(:), allocatable :: Auger0	! [fs] Auger-decay times for all shells
   real(8), dimension(:,:,:), allocatable :: f_auger0	! probability of Auger decay between two given shells
   integer :: EndCheck, i, counter, j, i_temp, i_temp2, shl_num_1, shl_num_2, shl_num_3, Nsiz
   integer :: Z0, A0, Yi0, Yo0, Date0, C0, I0, S0, X10
   real(8) :: AW0, temp
   character(72) test_line

   ! Start by reading the header of the first block
   call Read_ENDL_header(FN, File_name, INFO, Z0, A0, Yi0, Yo0, AW0, Date0, C0, I0, S0, X10)	! see below
   if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine

   EndCheck = 0	! to start with
   do while (Z0 < Z)	! skip lines until you find the element you need
      do while (EndCheck /= 1)	! until end of the block
         read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
         if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
      enddo
      EndCheck = 0	! restart checker
      ! Read the header of the next block:
      call Read_ENDL_header(FN, File_name, INFO, Z0, A0, Yi0, Yo0, AW0, Date0, C0, I0, S0, X10)	! see below
      if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
   enddo
   ! Now we found the element we need, proceed
   if (present(Mass)) Mass = AW0
   i = 0	! just to start
   
!    print*, 'Z', Z, Z0,  INFO, Z0, A0, Yi0, Yo0, AW0, Date0, C0, I0, S0, X10
   
   do while (Z0 == Z)	! read all the data for this element
      i = i + 1	! count blocks

      if (i == 1) then	! first time => allocate arrays
         ! First, count how many shells we have here:
         counter = 0
         do while (EndCheck /= 1)	! until end of the block
            counter = counter + 1		! count lines
            read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
            if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
         enddo
         ! Then, rewind the data back to read:
         do j = 1, counter	! rewind back to the start of the block
            backspace(FN)	! get back to read this line again
         enddo
         ! Finaly, allocate arrays to the given number:
         Nsiz = counter-1
         call allocate_EADL_data_arrays(Nsiz, Shl_dsgnr0, Shell_name0, Ip0, Ek0, Ne_shell0, Radiat0, f_rad0, Auger0, f_auger0)	! see below

         ! First, read all the shell designators to be used below:
         ! Now we have all arrays allocated, read data into them:
         EndCheck = 0
         counter = 0
         do while (EndCheck /= 1)	! read the entire block
            counter = counter + 1
            read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
            if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
            if (EndCheck /= 1) then	! if it is not the end of a block, read the line
               backspace(FN)	! get back to read this line again
               select case (I0)	! which data are in this block:
               case (912)	! number of electrons
                  read(FN,*,IOSTAT=INFO) Shl_dsgnr0(counter), Ne_shell0(counter)
                  ! Knowing the designator, fine shell name:
                  call find_shell_names(path_sep, INFO, Shl_dsgnr0(counter), 0, Shell_name0(counter))	! see below
                  if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
               case default
                  read(FN,*,IOSTAT=INFO) ! We already read it, so just skip such lines
               end select
            endif	! it is the end of the block, go to the next iteration of the cycle
         enddo ! while (EndCheck /= 1)
         ! Then, rewind the data back to read:
         do j = 1, counter	! rewind back to the start of the block
            backspace(FN)	! get back to read this line again
         enddo
      endif

      
      ! Now we have all arrays allocated, read data into them:
      EndCheck = 0
      counter = 0
      do while (EndCheck /= 1)	! read the entire block
         counter = counter + 1
         read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
         if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
         if (EndCheck /= 1) then	! if it is not the end of a block, read the line
            backspace(FN)	! get back to read this line again
            select case (I0)	! which data are in this block:
            case (912)	! number of electrons
               read(FN,*,IOSTAT=INFO) ! We already read it, so just skip such lines
               !read(FN,*,IOSTAT=INFO) Shl_dsgnr0(counter), Ne_shell0(counter)
               ! Knowing the designator, fine shell name:
               !call find_shell_names(path_sep, INFO, Shl_dsgnr0(counter), 0, Shell_name0(counter))	! see below
               if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
            case (913)	! binding energy
               read(FN,*,IOSTAT=INFO) temp, Ip0(counter)
               if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
               Ip0(counter) = Ip0(counter)*1.0d6	! [MeV] -> [eV]
            case (914)	! kinetic energy
               read(FN,*,IOSTAT=INFO) temp, Ek0(counter)
               if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
               Ek0(counter) = Ek0(counter)*1.0d6	! [MeV] -> [eV]
            case (915)	! average radius
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case (921)	! radiative level width
               read(FN,*,IOSTAT=INFO) temp, Radiat0(counter)
               if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
               Radiat0(counter) = 1d15*g_h/(g_e*Radiat0(counter)*1.0d6)	! [MeV] -> [fs]
            case (922)	! nonradiative level width
               read(FN,*,IOSTAT=INFO) temp, Auger0(counter)
               if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
               Auger0(counter) = 1d15*g_h/(g_e*Auger0(counter)*1.0d6)	! [MeV] -> [fs]
            case (931)	! radiative transition probability
               read(FN,*,IOSTAT=INFO) i_temp, temp
               if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
               call find_shell_by_designator(X10, Shl_dsgnr0, shl_num_1)
               call find_shell_by_designator(i_temp, Shl_dsgnr0, shl_num_2)
               f_rad0(shl_num_1,shl_num_2) = temp
            case (932)	! nonradiative transition probability
               read(FN,*,IOSTAT=INFO) i_temp, i_temp2, temp
               if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
               call find_shell_by_designator(X10, Shl_dsgnr0, shl_num_1)
               call find_shell_by_designator(i_temp, Shl_dsgnr0, shl_num_2)
               call find_shell_by_designator(i_temp2, Shl_dsgnr0, shl_num_3)
               f_auger0(shl_num_1,shl_num_2,shl_num_3) = temp
            case (933)	! particles per initial vacancy
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case (934)	! energy of particles per initial vacancy
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case (935)	! average energy to the residual atom, i.e., local deposition, per initial vacancy
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            end select
         endif	! it is the end of the block, go to the next iteration of the cycle
      enddo ! while (EndCheck /= 1)
      ! After the end of the block, read the header of the next one:
      call Read_ENDL_header(FN, File_name, INFO, Z0, A0, Yi0, Yo0, AW0, Date0, C0, I0, S0, X10)	! see below
      if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
   enddo ! (Z0 == Z)
   
   ! Save parameteres required by the user:
   if (present(Shl_dsgnr)) then
      if(.not.allocated(Shl_dsgnr)) allocate(Shl_dsgnr(Nsiz))	! EADL shell designator
      Shl_dsgnr = Shl_dsgnr0
   endif
   if (present(Shell_name)) then
      if(.not.allocated(Shell_name)) allocate(Shell_name(Nsiz))	! names of the shells
      Shell_name = Shell_name0
   endif
   if (present(Ip)) then
      if(.not.allocated(Ip)) allocate(Ip(Nsiz))		! [eV] ionization potentials for all shells
      Ip = Ip0
   endif
   if (present(Ek)) then
      if(.not.allocated(Ek)) allocate(Ek(Nsiz))		! [eV] mean kinetic energy of all shells
      Ek = Ek0
   endif
   if (present(Ne_shell)) then 
      if(.not.allocated(Ne_shell)) allocate(Ne_shell(Nsiz))	! number of electron in each shell
      Ne_shell = Ne_shell0
   endif
   if (present(Radiat)) then
      if(.not.allocated(Radiat)) allocate(Radiat(Nsiz))	! [fs] radiative-decay times for all shells
      Radiat = Radiat0
   endif
   if (present(f_rad)) then
      if(.not.allocated(f_rad)) allocate(f_rad(Nsiz,Nsiz))	! probability of radiative decay from a given shell
      f_rad = f_rad0
   endif
   if (present(Auger)) then
      if(.not.allocated(Auger)) allocate(Auger(Nsiz))	! [fs] Auger-decay times for all shells
      Auger = Auger0
   endif
   if (present(f_auger)) then
      if(.not.allocated(f_auger)) allocate(f_auger(Nsiz,Nsiz,Nsiz))	! probability of Auger radiative decay between a given shells
      f_auger = f_auger0
   endif
   
   ! Clean up at the end:
   if (allocated(Shl_dsgnr0)) deallocate(Shl_dsgnr0, Shell_name0, Ip0, Ek0, Ne_shell0, Radiat0, f_rad0, Auger0, f_auger0)
9998 continue
9997   format(71X,I1)			! block separator
end subroutine Read_EADL_data


subroutine allocate_EADL_data_arrays(Nsiz, Shl_dsgnr, Shell_name, Ip, Ek, Ne_shell, Radiat, f_rad, Auger, f_auger)
   integer, intent(in) :: Nsiz	! known array size
   integer, dimension(:), allocatable, intent(inout) :: Shl_dsgnr	! EADL shell designator
   character(11), dimension(:), allocatable, intent(inout) :: Shell_name	! names of the shells
   real(8), dimension(:), allocatable, intent(inout) :: Ip		! [eV] ionization potentials for all shells
   real(8), dimension(:), allocatable, intent(inout) :: Ek		! [eV] mean kinetic energy of all shells
   real(8), dimension(:), allocatable, intent(inout) :: Ne_shell	! number of electron in each shell
   real(8), dimension(:), allocatable, intent(inout) :: Radiat	! [fs] radiative-decay times for all shells
   real(8), dimension(:,:), allocatable, intent(inout) :: f_rad	! probability of radiative decay from a given shell
   real(8), dimension(:), allocatable, intent(inout) :: Auger	! [fs] Auger-decay times for all shells
   real(8), dimension(:,:,:), allocatable, intent(inout) :: f_auger	! probability of Auger decay between two given shells
   ! allocate arrays:
   allocate(Shl_dsgnr(Nsiz))
   allocate(Shell_name(Nsiz))
   allocate(Ip(Nsiz))
   allocate(Ek(Nsiz))
   allocate(Ne_shell(Nsiz))
   allocate(Radiat(Nsiz))
   allocate(f_rad(Nsiz,Nsiz))
   allocate(Auger(Nsiz))
   allocate(f_auger(Nsiz,Nsiz,Nsiz))
   ! set default values:
   Shl_dsgnr = 0
   Shell_name = ''
   Ip = 1.0d23	! [eV]
   Ek = 0.0d0	! [eV]
   Ne_shell = 0
   Radiat = 1.0d24	! [fs]
   f_rad = 0.0d0
   Auger = 1.0d23		! [fs]
   f_auger = 0.0d0
end subroutine allocate_EADL_data_arrays



subroutine Read_ENDL_header(FN, File_name, INFO, Z, A, Yi, Yo, AW, Date, C, I, S, X1, Iflag)
   integer, intent(in) :: FN	! file number of ENDL database
   character(100), intent(in) :: File_name	! file name of ENDL database
   integer, intent(inout) :: INFO	! info whether file read well
   integer, intent(out), optional :: Z 	! atomic number
   integer, intent(out), optional :: A	! mass number (in all cases=0, for elemental data)
   integer, intent(out), optional :: Yi	! incident particle designator (see Table II)
   integer, intent(out), optional :: Yo	! outgoing particle designator (see Table II)
   real(8), intent(out), optional :: AW	! atomic mass (amu)
   integer, intent(out), optional :: Date	! date of evaluation (YYMMDD)
   integer, intent(out), optional :: Iflag	! interpolation flag (used in EPDL, not EADL)
   integer, intent(out), optional :: C	! reaction descriptor (see Table II)
   integer, intent(out), optional :: I		! reaction property (see Table II)
   integer, intent(out), optional :: S	! reaction modifier (see Table II)
   integer, intent(out), optional :: X1	! subshell designator (see Table VI)
   !----------------------------------------------
   integer :: Z0, A0, Yi0, Yo0, Date0, C0, I0, S0, Iflag0
   real(8) :: AW0, X10
!    character(72) :: test
   ! First header line:
   
!    read(FN,'(a)') test
!    print*, ':::: ', test
!    backspace(FN)
   if (present(Iflag)) then ! EPDL database:
      read(FN,'(i3, i3, 1X, i2, 1X, i2, 1X, E11.4, 1x, i6, i1)',IOSTAT=INFO) Z0, A0, Yi0, Yo0, AW0, Date0, Iflag0
   else ! EADL database:
      read(FN,'(i3, i3, 1X, i2, 1X, i2, 1X, E11.4, 1x, i6)',IOSTAT=INFO) Z0, A0, Yi0, Yo0, AW0, Date0
   endif
!    print*, '1) Read_ENDL_header', INFO
!    print*, Z0, A0, Yi0, Yo0, AW0, Date0
   if (INFO /= 0) then	! something went wrong while reading this line
      goto 9998	! exit the subroutine
   endif
   
   ! Second header line:
   read(FN,'(i2, i3, i3, 14X, E11.4)',IOSTAT=INFO) C0, I0, S0, X10
!    print*, '2) Read_ENDL_header', INFO
!    print*, C0, I0, S0, X10
   if (INFO /= 0) then	! something went wrong while reading this line
      goto 9998	! exit the subroutine
   endif
   
   ! Save all required output:
   if (present(Z)) Z = Z0
   if (present(A)) A = A0
   if (present(Yi)) Yi = Yi0
   if (present(Yo)) Yo = Yo0
   if (present(AW)) AW = AW0
   if (present(Date)) Date = Date0
   if (present(Iflag)) Iflag = Iflag0
   if (present(C)) C = C0
   if (present(I)) I = I0
   if (present(S)) S = S0
   if (present(X1)) X1 = NINT(X10)
   
9998 continue
end subroutine Read_ENDL_header
 
 

subroutine select_imin_imax(imin, imax, shl_dsgntr)
   integer, intent(in) :: shl_dsgntr
   integer, intent(out) :: imin, imax ! according to EADL database, those are corresponding shells:
   select case (shl_dsgntr) ! sum subshells:
   case (2)
      imin = 3
      imax = 6
   case (4)
      imin = 5
      imax = 6
   case (7)
      imin = 8
      imax = 14
   case (9)
      imin = 10
      imax = 11
   case (12)
      imin = 13
      imax = 14
   case (15)
      imin = 16
      imax = 25
   case (26)
      imin = 27
      imax = 39
   case (40)
      imin = 41
      imax = 56
   case (57)
      imin = 58
      imax = 61
   case default
      imin = 0
      imax = 0
   endselect
end subroutine select_imin_imax


subroutine next_designator(Shl_dsgtr, Last_shl) ! find the designator for the VB (next shell after last one)
   integer, intent(in) :: Last_shl
   integer, intent(out) :: Shl_dsgtr
   select case (Last_shl) ! sum subshells:
   case (:1)
      Shl_dsgtr = 2
   case (2:6)
      Shl_dsgtr = 7
   case (7:14)
      Shl_dsgtr = 15
   case (15:25)
      Shl_dsgtr = 26
   case (26:39)
      Shl_dsgtr = 40
   case (40:)
      Shl_dsgtr = 57
   endselect
end subroutine next_designator


! Interpolation of photoabsorbtion cross-section according to EADL database:
subroutine Interpolate_EPDL(Iflag, E1, E2, Sigma1, Sigma2, E_needed, OUT_value)
   integer, intent(in) :: Iflag
   real(8), intent(in) :: E1, E2, Sigma1, Sigma2, E_needed
   real(8), intent(out) :: OUT_value
   real(8) E2log, E1log, E_needed_log, Sigma1log, Sigma2log
   select case(Iflag) ! what interpolation to use:
      case(0,2) ! linear x and y
         OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2 - E1)*(E_needed - E1)
      case(3)	! logarithmic x, linear y
         E2log = log(E2)
         E1log = log(E1)
         E_needed_log = log(E_needed)
         OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2log - E1log)*(E_needed_log - E1log)
      case(4)	! linear x, logarithmic y
         Sigma1log = log(Sigma1)
         Sigma2log = log(Sigma2)
         OUT_value = Sigma1log + (Sigma2log - Sigma1log)/(E2 - E1)*(E_needed - E1)
         OUT_value = exp(OUT_value)
      case(5)	! logarithmic x and y
         E2log = log(E2)
         E1log = log(E1)
         E_needed_log = log(E_needed)
         Sigma1log = log(Sigma1)
         Sigma2log = log(Sigma2)
         OUT_value = Sigma1log + (Sigma2log - Sigma1log)/(E2log - E1log)*(E_needed_log - E1log)
         OUT_value = exp(OUT_value)
      case default ! linear x and y
         OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2 - E1)*(E_needed - E1) 
   end select
end subroutine Interpolate_EPDL


subroutine find_shell_by_designator(shl_dsgntr, sch_array, shl_num)
   integer, intent(in) :: shl_dsgntr	! shell designator
   integer, dimension(:), intent(in) :: sch_array	! array of all designators
   integer, intent(out) :: shl_num	! shell number in arrays
   call Find_in_array_monoton(dble(sch_array), dble(shl_dsgntr), shl_num)	! module "Little_subroutines"
end subroutine find_shell_by_designator


subroutine find_shell_names(path_sep, INFO, Shl_dsgnr, Ind, Shell_name, PQN)
   character(1), intent(in) :: path_sep	! path separator in the system
   integer, intent(inout) :: INFO		! INFO=0 file found and read well, INFO/=0 trouble reading file
   integer, intent(in) :: Shl_dsgnr	! shell designator
   integer, intent(in) :: Ind	! which format of shell name to use: 0 = spectroscopic (e.g. K-shell, L-shell...); 1 = atomic (e.g. 1s1/2, 2p3/2...)
   character(11), intent(out), optional :: Shell_name	! names
   integer, intent(out), optional :: PQN	!Principal quantum number
   !------------------
   type(Shell_designator) :: Shell_ENDL
   integer :: FN, i, N
   character(100) :: Folder_name, File_name, temp
   logical file_exist, file_open
   ! Check if file with the database exists:
   Folder_name = 'INPUT_DATA'//path_sep//'Atomic_parameters'  ! here we keep databases
   File_name = trim(adjustl(Folder_name))//path_sep//'ENDL_shell_designators.dat'	! this name is fixed
   inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
   if (.not. file_exist) then
      write(*,'(a,a,a)') 'File ',trim(adjustl(File_name)),' is not found.'
      INFO = 1
      goto 9996
   else
      INFO = 0
   endif
   ! If file found, read it:
   open(newunit=FN, FILE=trim(adjustl(File_name)), ACTION='READ')		! open database
   ! Read all lines in the file:
   SRCH:do while (INFO == 0)
      read(FN,*,IOSTAT=INFO) Shell_ENDL%Shl, Shell_ENDL%Shell_spect, Shell_ENDL%Shell_at
!       print*, 'SRCH', INFO, Shell_ENDL%Shl, Shell_ENDL%Shell_spect, Shell_ENDL%Shell_at
      if (INFO /= 0 ) then	! end of file or error while reading
         write(temp,'(i3)') Shl_dsgnr
         write(*,'(a,a,a)') 'Shell designator #',trim(adjustl(temp)),' is not found.'
         exit SRCH
      elseif (Shell_ENDL%Shl == Shl_dsgnr) then	! we found the element we needed
         select case (Ind)	! in which format shell name should be used:
         case (1)	! atomic shell name
            if (present(Shell_name)) Shell_name = Shell_ENDL%Shell_at
         case default	! spectroscopic shell name
            if (present(Shell_name)) Shell_name = Shell_ENDL%Shell_spect
         endselect
         ! Get principal quantum number if required:
         if (present(PQN)) then
            Shell_ENDL%Shell_at = trim(adjustl(Shell_ENDL%Shell_at))
            read(Shell_ENDL%Shell_at,'(i1)') PQN
         endif
         exit SRCH	! once we found the required values, exit the search cycle
      endif
   enddo SRCH

9996 continue
   inquire(file=trim(adjustl(File_name)),opened=file_open) ! check if input file is there
   if (file_open) close(FN)
end subroutine find_shell_names


END MODULE Dealing_with_EADL
