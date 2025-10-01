! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2014-2020
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with DOS files:
MODULE Dealing_with_DOS
use Universal_constants
use Objects
use Dealing_with_files, only : Count_lines_in_file, read_file, close_file
use Read_input_data, only : m_input_folder, m_folder_DOS, m_folder_materials
use Little_subroutines, only : interpolate_data_on_grid, Find_in_array_monoton, Fermi_function

implicit none

real(8) :: m_dE

parameter(m_dE = 0.05d0) ! [eV] energy grid step for DOS

 contains
 
subroutine read_DOS_files(used_target, numpar, Err)
   type(Matter), intent(inout), target :: used_target	! parameters of the target
   type(Num_par), intent(inout), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !------------------------------
   real(8) :: Int_DOS, E_cur, eps
   integer :: FN, i, j
   character(250) :: Path, Error_descript, File_name
   logical :: file_exist
   character, pointer :: path_sep
   type(Density_of_states), pointer :: DOS
   
   if (numpar%verbose) print*, 'Getting DOS of materials...'
   
   ! The limit for smallest effective hole allowed:
   eps = 1.0d-2     ! [me]
   
   path_sep => numpar%path_sep  ! path separator
   
   ! For all target materials:
   TRGT:do i = 1, used_target%NOC   ! number of different components of the target (layers, clusters, etc.)   
      file_exist = .false.   ! to start with
      ! Get the DOS:
      ! Check if file with DOS is in the material folder:
      Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
      ! Folder with the particular material:
      Path = trim(adjustl(Path))//path_sep//trim(adjustl(used_target%Material(i)%Name))
      File_name = trim(adjustl(Path))//path_sep//trim(adjustl(used_target%Material(i)%DOS_file))
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there

      ! If not, check if the file with DOS is present in DOS folder:
      if (.not.file_exist) then
         ! Check if the folder with material parameters is present:
         Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_DOS))
         File_name = trim(adjustl(Path))//path_sep//trim(adjustl(used_target%Material(i)%DOS_file))
         inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
      endif

      if (file_exist) then ! use DOS from the file:
         open(newunit = FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
         if (numpar%verbose) print*, 'Reading DOS from file '//trim(adjustl(File_name))
         ! Read DOS from it into the arrays:
         call read_DOS_from_file(FN, File_name, used_target%Material(i)%DOS, Err)    ! below
         
         ! Close the file after reading DOS from it:
         call close_file('close',FN=FN)
      else ! construct free-electron DOS:
         if (numpar%verbose) print*, 'Creating free-electron DOS...'
         call construct_free_electron_DOS(used_target%Material(i)%DOS)
      endif
      
      ! Make sure we didn't get any negative values during interpolation onto the grid:
      do j = 1, size(used_target%Material(i)%DOS%E)
         if (used_target%Material(i)%DOS%DOS(j) < 0.0d0) used_target%Material(i)%DOS%DOS(j) = 0.0d0
      enddo

      ! Set electron distribution function on the given grid (the same grid as used for DOS):
      call set_electron_distribution(used_target%Material(i)%DOS%E, used_target%Material(i)%fe, used_target%Material(i)%DOS%E_f, used_target%Material(i)%T_eV) ! below
      
      ! Normalize DOS to the given number of electrons:
      ! get number of electrons as an integral of the unnormalized DOS
      Int_DOS = SUM(used_target%Material(i)%DOS%DOS(:) * used_target%Material(i)%fe(:)) * m_dE  ! for equal-spaced grid
      ! normalized DOS to the number of electrons in the VB (or CB) PER MOLECULE (not per atom):
      used_target%Material(i)%DOS%DOS(:) = used_target%Material(i)%DOS%DOS(:) * used_target%Material(i)%N_VB_el / Int_DOS
      
      DOS => used_target%Material(i)%DOS
      
      ! Calculate the effective mass for holes:
      ! VB is filled from top to bottom:
      Int_DOS = 0.0d0  ! start to integrate:
      do j = DOS%N_VB_top, 1, -1
         Int_DOS = Int_DOS + DOS%DOS(j) * m_dE
         DOS%k(j) = (3.0d0*2.0d0*g_Pi*g_Pi/2.0d0*Int_DOS*used_target%Material(i)%At_Dens/dble(SUM(used_target%Material(i)%Elements(:)%percentage))*1d6)**(1.0d0/3.0d0)  ! [1/m]
         E_cur = abs(DOS%E(j) - DOS%E(DOS%N_VB_top)) ! counting down from the top of VB
         if  (E_cur > 1.0d-6) then
            DOS%Eff_m(j) = g_h*g_h*DOS%k(j)*DOS%k(j)/(2.0d0*E_cur*g_e)/g_me ! [electron mass]
            if (DOS%Eff_m(j) < eps) DOS%Eff_m(j) = 1.0d10   ! make it immobile, if for some reason the mass is too low
         else
            DOS%Eff_m(j) = 1.0d10   ! immobile electron/hole at the bottom
         endif
      enddo
      ! CB is filled from bottom to top:
      Int_DOS = 0.0d0  ! restart to integrate:
      do j = DOS%N_CB_bottom, size(DOS%E)
         Int_DOS = Int_DOS + DOS%DOS(j) * m_dE
         DOS%k(j) = (3.0d0*2.0d0*g_Pi*g_Pi/2.0d0*Int_DOS*used_target%Material(i)%At_Dens / &
                             dble(SUM(used_target%Material(i)%Elements(:)%percentage))*1d6)**(1.0d0/3.0d0)  ! [1/m]
         E_cur = abs(DOS%E(j) - DOS%E(DOS%N_CB_bottom)) ! counting from the bottom of CB
         if  (E_cur > 1.0d-6 ) then
            DOS%Eff_m(j) = g_h*g_h*DOS%k(j)*DOS%k(j)/(2.0d0*E_cur*g_e)/g_me ! [electron mass]
            if (DOS%Eff_m(j) < eps) DOS%Eff_m(j) = 1.0d10   ! make it immobile, if for some reason the mass is too low
         else
            DOS%Eff_m(j) = 1.0d10   ! immobile electron/hole at the bottom
         endif
      enddo

      ! Get the coefficient that joins free-electron DOS to the given DOS at its top:
      DOS%alpha_CB = set_alpha_CB(DOS%DOS(:), DOS%E_CB_bottom, DOS%E_CB_top) ! below
      
!       print*, 'read_DOS_files', DOS%N_CB_bottom, DOS%E_CB_bottom, DOS%E(DOS%N_CB_bottom) 
      ! Shift the DOS to the top of VB:
      DOS%E(:) = DOS%E(:) - DOS%E_VB_top
      ! Adjust the other parameters accrodingly:
      DOS%E_VB_bottom = DOS%E_VB_bottom - DOS%E_VB_top
      DOS%E_CB_bottom = DOS%E_CB_bottom - DOS%E_VB_top
      DOS%E_CB_top = DOS%E_CB_top - DOS%E_VB_top
      DOS%E_VB_top = 0.0d0
      
      ! Make sure band gap value is consistent with the gap in the DOS:
      E_cur = DOS%E_CB_bottom - (DOS%E_VB_top + DOS%Egap)  ! error to be corrected
      DOS%E_CB_bottom = DOS%E_VB_top + DOS%Egap    ! correct the CB starting point for the numerical precision shift
      DOS%E(DOS%N_CB_bottom:size(DOS%E)) = DOS%E(DOS%N_CB_bottom:size(DOS%E)) - E_cur  ! correct all CB for the numerical precision shift
!       print*, 'read_DOS_files', DOS%N_CB_bottom, DOS%E_CB_bottom, DOS%E(DOS%N_CB_bottom) 
!       pause
 
      ! Get the integrals to be used later in MC:
      call set_integral_DOS_fe(DOS%E, used_target%Material(i)%fe, DOS%DOS(:), used_target%Material(i)%Integral_DOS_fe) ! below
   
      nullify (DOS)
   enddo TRGT
   
   if (numpar%verbose) print*, 'DOS of all materials read successfully.'
9992 continue
end subroutine read_DOS_files


pure function set_alpha_CB(DOS, E_CB_bottom, E_CB_top) result(alpha)
   real(8) alpha    ! [eV]
   real(8), dimension(:), intent(in) :: DOS ! DOS
   real(8), intent(in) ::  E_CB_bottom, E_CB_top
   integer :: Nsiz
   Nsiz = size(DOS)
   alpha = DOS(Nsiz)/sqrt(E_CB_top - E_CB_bottom)
end function set_alpha_CB



pure subroutine set_electron_distribution(Egrid, fe, mu, Te)
   real(8), dimension(:), intent(in) :: Egrid
   real(8), dimension(:), allocatable, intent(inout) :: fe  ! electron distribution function
   real(8), intent(in) :: mu, Te  ! chemical potential [eV] and temperature [eV]
   integer :: Nsiz, i
   Nsiz = size(Egrid)
   if (.not.allocated(fe)) allocate(fe(Nsiz))
   do i = 1, Nsiz
      fe(i) = Fermi_function(mu, Te, Egrid(i))    ! module "Little_subroutines"
   enddo
end subroutine set_electron_distribution


pure subroutine set_integral_DOS_fe(Egrid, fe, DOS, Integral_DOS_fe)
   real(8), dimension(:), intent(in) :: Egrid
   real(8), dimension(:), intent(in) :: fe  ! electron distribution function
   real(8), dimension(:), intent(in) :: DOS ! DOS
   real(8), dimension(:), allocatable, intent(inout) :: Integral_DOS_fe
   integer :: Nsiz, i
   Nsiz = size(Egrid)
   if (.not.allocated(Integral_DOS_fe)) allocate(Integral_DOS_fe(Nsiz))
   do i = 1, Nsiz
      Integral_DOS_fe(i) = DOS(i) * fe(i)
   enddo
end subroutine set_integral_DOS_fe



subroutine select_energy_DOS(Egrid, DOS, Integral_DOS_fe, Egap, dE, alpha_CB, E_out)
   real(8), dimension(:), intent(in) :: Egrid     ! energy grid in DOS [eV]
   real(8), dimension(:), intent(in) :: DOS     ! DOS
   real(8), dimension(:), intent(in) :: Integral_DOS_fe ! fe(E)*DOS(E)
   real(8), intent(in) :: Egap  ! [eV] band gap
   real(8), intent(in) :: dE    ! transferred energy [eV]
   real(8), intent(in) :: alpha_CB  ! coefficient for transfering DOS into free-electron DOS
   real(8), intent(out) :: E_out    ! the energy level from where an electron is being excited (counted from the TOP OF VB!) [eV]
   !-----------------------------
   real(8), dimension(:), allocatable :: Int_Boltzmann
   real(8) :: RN, eps
   integer :: i, Nsiz, i_plus, i2, i_out, i_start
   ! Set accepteble margin of precision for the energy level:
   eps = 1.0d-3

   if (isnan(dE)) print*, 'Error in select_energy_DOS: dE=', dE

   STBG:if (abs(dE - Egap) < eps) then   ! it's the first grid point from the band gap
!       E_out = -eps   ! [eV]
      E_out = -abs(dE - Egap)    ! [eV]
!       print*, 'select_energy_DOS 0:', E_out, dE, Egap
   else STBG
      ! The maximal energy level on the grid:
      Nsiz = size(Egrid)
      allocate(Int_Boltzmann(Nsiz))
   
      ! In case not all of the VB can be ionized with this energy
      if ( dE < abs(Egrid(1))+Egap ) then  ! find which states can be ionized
         call Find_in_array_monoton(Egrid, -dE+Egap, i_start)   ! module "Little_subroutines"
         i_start = i_start + 1
         !print*, 'select_energy_DOS', dE, i_start, Egrid(i_start)
      else ! all the VB can be ionized
         i_start = 1
      endif
   
      ! Find the index for the shift by dE:
      call Find_in_array_monoton(Egrid, Egrid(1) + dE, i_plus)   ! module "Little_subroutines"

      ! Form the Boltzmann-integral-like function to be used for sampling probability:
      Int_Boltzmann = 0.0d0    ! to start with
      do i = i_start, Nsiz
         i2 = i + i_plus
         if (i2 <= Nsiz) then  ! defined DOS
            if (i == 1) then
               Int_Boltzmann(i) = Integral_DOS_fe(i) * DOS(i2)
            else
               Int_Boltzmann(i) = Int_Boltzmann(i-1) + Integral_DOS_fe(i) * DOS(i2)
            endif
         else  ! free-electron DOS state an electron ends up in
            if (i == 1) then
               Int_Boltzmann(i) = Integral_DOS_fe(i) * alpha_CB*sqrt( Egrid(i) + dE )
            else
               Int_Boltzmann(i) = Int_Boltzmann(i-1) + Integral_DOS_fe(i)* alpha_CB*sqrt( Egrid(i) + dE )
            endif
         endif
!          write(*,'(a,i3,f,f,f)') 'select_energy_DOS', i, Int_Boltzmann(i),  Integral_DOS_fe(i), dE
      enddo

      ! In case DOS is too small within a too small interval around the band gap
      if (Int_Boltzmann(Nsiz) < eps) then   ! just set it close to zero
!          E_out = -eps   ! [eV]
         E_out = -abs(dE - Egap)    ! [eV]
!          print*, 'select_energy_DOS 1:', E_out, dE, Egap
      else  ! sample a proper energy interval
         ! Get the random energy level from DOS according to DOS and fe:
         call random_number(RN)
   
         if (isnan(RN*Int_Boltzmann(Nsiz))) then
            print*, 'Error in select_energy_DOS #2:'
            print*, RN, Int_Boltzmann(Nsiz), Nsiz
         endif
         
         ! Find the sampled energy level:
         call Find_in_array_monoton(Int_Boltzmann, RN*Int_Boltzmann(Nsiz), i_out)   ! module "Little_subroutines"
   
         ! Now set the energy level an electron is ionized from:
!          if (Egap > eps) then   ! semiconduction or inslator
            E_out = Egrid(i_out+1)
!          else   ! metal
!             E_out = Egrid(i_out)
!          endif
!          if (abs(E_out - Egap) < eps) print*, 'select_energy_DOS', E_out, i_out, Egrid(i_out+1), Egrid(i_out), RN
      endif ! (Int_Boltzmann(Nsiz) < eps)
      
      deallocate(Int_Boltzmann)
   endif STBG
   
   if (isnan(E_out)) then
      print*, 'Error in select_energy_DOS #3:', E_out, dE
   endif
   
!    print*, dE, E_out, Egrid(i_out+1), RN*Int_Boltzmann(Nsiz)
!    if ( dE < abs(Egrid(1))+Egap ) then 
!       do i = 1, Nsiz
!          write(*,'(i3, f, f, f, f)') i, Egrid(i), Integral_DOS_fe(i), Int_Boltzmann(i), RN*Int_Boltzmann(Nsiz)
!       enddo
!    endif
end subroutine select_energy_DOS

 
subroutine construct_free_electron_DOS(DOS)
   type(Density_of_states), intent(inout) :: DOS    ! DOS
   real(8) :: eps, E_VB_bottom, E_VB_top, E_CB_bottom, E_CB_top
   integer :: N_lines, i
   
   eps = 1.0d-6 ! margin of error
   if (abs(DOS%E_VB_bottom-DOS%E_VB_top) < eps) then ! if they coinside, they are undefined, use default values:
      E_VB_bottom =  DOS%E_f - 10.0d0 ! [eV]
      E_VB_top = DOS%E_f  ! [eV]
      E_CB_bottom = DOS%E_f + DOS%Egap ! [eV]
      E_CB_top = E_CB_bottom + 10.0d0   ! [eV]
      ! And save these values to be used later for valence holes:
      DOS%E_VB_bottom = E_VB_bottom
      DOS%E_VB_top = E_VB_top
      DOS%E_CB_bottom = E_CB_bottom
      DOS%E_CB_top = E_CB_top
   else
      E_VB_bottom = DOS%E_VB_bottom
      E_VB_top = DOS%E_VB_top
      E_CB_bottom = DOS%E_CB_bottom
      E_CB_top = DOS%E_CB_top
   endif
   
   ! Create the energy grid for DOS to be used with the given spacing:
   call create_DOS_grid(DOS%E, E_VB_bottom, E_CB_top) ! below
   N_lines = size(DOS%E)
   allocate(DOS%DOS(N_lines), DOS%k(N_lines), DOS%Eff_m(N_lines))
   
   ! Construct DOS by parts: valence and conduction (not normalized!):
   if (DOS%Egap > m_dE) then ! it is a dielectric
      DOS%DOS = 0.0d0   ! default values are zero (and stay zero within the gap)
      ! Find the point on the grid where Fermi energy lies:
      call Find_in_array_monoton(DOS%E, DOS%E_f, DOS%N_VB_top)  ! module "Little_subroutines"
      do i = DOS%N_VB_top, 1, -1    ! VB is filled from top to bottom
         DOS%DOS(i) = sqrt(abs(DOS%E(i) - DOS%E_f))
      enddo
      ! Find the point on the grid where the conduction band starts:
      call Find_in_array_monoton(DOS%E, E_CB_bottom, DOS%N_CB_bottom, from_above=.true.)  ! module "Little_subroutines"
      do i = DOS%N_CB_bottom, size(DOS%E)    ! CB is filled from bottom to top
         DOS%DOS(i) = sqrt(abs(DOS%E(i) - E_CB_bottom))
      enddo
   else ! it is a metal
      DOS%N_VB_top = 1  ! conduction band starts from the first element of the DOS array
      DOS%N_CB_bottom = 1  ! conduction band starts from the first element of the DOS array
      do i = 1, size(DOS%E)    ! CB is filled from bottom to top
         DOS%DOS(i) = sqrt(abs(DOS%E(i) - E_VB_bottom))
      enddo
   endif
end subroutine construct_free_electron_DOS



subroutine read_DOS_from_file(FN, File_name, DOS, Err)
   integer, intent(in) :: FN    ! file number (file must be opened prior to calling this function)
   character(*), intent(in) :: File_name    ! file with DOS
   type(Density_of_states), intent(inout) :: DOS    ! DOS
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------------
   integer :: N_lines, i, Nread, Reason, count_lines
   real(8), dimension(:), allocatable :: E_temp, DOS_temp
   logical :: read_well
   character(200) :: Error_descript
   
   ! Count how many energy grid points are in the file for this DOS:
   call Count_lines_in_file(FN, Nread)  ! module "Dealing_with_files"

   ! Read DOS from the file:
   ! allocate temporary arrays for a raw DOS:
   allocate(E_temp(Nread), DOS_temp(Nread))
   count_lines = 0  ! to start counting lines from
   do i = 1, Nread  ! read the file
      read(FN,*,IOSTAT=Reason)  E_temp(i), DOS_temp(i)    ! energy [eV], DOS [1/eV]
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         call Save_error_details(Err, 2, Error_descript)	! module "Objects"
         goto 9990
      endif
!       print*, i, E_temp(i), DOS_temp(i)
   enddo
!    pause 'DOS_temp'
   
   ! Create the energy grid for DOS to be used with the given spacing:
   call create_DOS_grid(DOS%E, E_temp(1), E_temp(size(E_temp))) ! below
   ! Allocate other parameters of the DOS object:
   N_lines = size(DOS%E)
   allocate(DOS%DOS(N_lines), DOS%k(N_lines), DOS%Eff_m(N_lines))
   ! Interpolate DOS on the grid to be used everywhere:
   call interpolate_data_on_grid(E_temp, DOS_temp, DOS%E, DOS%DOS)  ! module "Little_subroutines"
   ! Clean up the temporary arrays:
   deallocate(E_temp, DOS_temp)

   ! Find the point on the grid where Fermi energy lies:
   call Find_in_array_monoton(DOS%E, DOS%E_f, DOS%N_VB_top)  ! module "Little_subroutines"
   ! Identify the bottom and top indices for VB and CB:
   if (DOS%Egap > m_dE) then ! it is a dielectric
      ! Find the point on the grid where the conduction band starts:
      call Find_in_array_monoton(DOS%E, (DOS%E_f + DOS%Egap), DOS%N_CB_bottom)  ! module "Little_subroutines"
   else ! it is a metal
!       DOS%N_VB_top = 1  ! conduction band starts from the first element of the DOS array
!       DOS%N_CB_bottom = 1  ! conduction band starts from the first element of the DOS array
      DOS%N_CB_bottom = DOS%N_VB_top
   endif
   ! Define bottom and top of the bands:
   DOS%E_VB_bottom = DOS%E(1)
   DOS%E_VB_top = DOS%E(DOS%N_VB_top)
   DOS%E_CB_bottom = DOS%E(DOS%N_CB_bottom)
   DOS%E_CB_top = DOS%E(size(DOS%E))
   
!    print*, 'DOS:', DOS%E_VB_bottom , DOS%E_VB_top
!    do i = 1, N_lines
!       print*, DOS%E(i), DOS%DOS(i)
!    enddo
!    pause 'read_DOS_from_file'

9990 continue
end subroutine read_DOS_from_file


pure subroutine create_DOS_grid(E, E_VB_bottom, E_CB_top) ! creates DOS grid in a free form (without predefined array from a file)
   real(8), dimension(:), allocatable, intent(inout) :: E  ! energy grid for DOS [eV]
   real(8), intent(in) :: E_VB_bottom, E_CB_top ! bottom and top of the DOS
   integer :: Nsiz, i
   
   Nsiz = FLOOR(abs(E_VB_bottom - E_CB_top)/m_dE)    ! grid size for the given grid step m_dE
   
   ! Knowing the size, allocate the energy grid for DOS array:
   allocate(E(Nsiz))
   
   ! Start exactly from the given starting point (the last point may shift though within the amount of m_dE):
   E(1) = E_VB_bottom
   do i = 2, Nsiz 
      E(i) = E(i-1) + m_dE     ! fill the energy grid with equal spacing of m_dE
   enddo
end subroutine create_DOS_grid

 
END MODULE Dealing_with_DOS
