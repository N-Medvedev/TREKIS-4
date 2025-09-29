! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2025
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains all cross sections for positrons Bremsstrahlung photon emission, etc.
! References used:
! [1] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2008: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2009)
! [2] F. Salvat, J.M. Fernhdez-Varea, NIMB 63 (1992) 255-269


module CS_muon_Bremsstrahlung
use Universal_constants
use CDF_Ritchi
use Relativity
use Objects
use Read_input_data, only: m_muon_CS, m_input_folder, m_muon_Brems_CS, m_folder_materials
use Dealing_with_files, only: read_file
use CS_general_tools, only: MFP_from_sigma
use Little_subroutines, only: interpolate_data_single
use CS_electrons_Bremsstrahlung, only: m_Wth, Bremsstrahlung_total_CS

implicit none


 contains
 

!-----------------------------------------------------------------------------
! Subroutines for Bremsstrahlung:
subroutine get_muon_Brems(Material, numpar, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material	!material parameters of each target that it's constructed of
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------------------------
   integer :: i, j, k, m, Ngrid, N_targets, N_elements, N_shells, FN, Reason, count_lines, Nsiz
   real(8) :: sigma, N_elem, elem_contrib, MFP
   character(200) :: Path, Folder_with_CS, command, File_name, Model_name, Path_total
   logical :: file_exist, read_well
   !--------------------------------
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   
   write(*, '(a)', advance='no') ' Obtaining muon Bremsstrahlung cross sections...'
   
   path_sep => numpar%path_sep
   Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_muon_CS))	! Electron CSs are storred in this folder
   
   N_targets = size(Material)	! that's how many different targets user specified
   
   ! Get muon MFPs:
   TRGT:do i = 1, N_targets	! for each target
      N_elem = dble(SUM(Material(i)%Elements(:)%percentage))	! number of elements in this compound material
      N_elements = size(Material(i)%Elements)	! that's how many different elements are in this target
      LMNT:do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         ! Check if folder for this element already exists:
         Folder_with_CS = trim(adjustl(Path))//path_sep//trim(adjustl(Element%Name))
         inquire(DIRECTORY=trim(adjustl(Folder_with_CS)),exist=file_exist)    ! check if input file excists
         if (.not.file_exist) then	! to make sure that such a folder is present (even if empty)
            ! Create a new directory for output files:
            command='mkdir '//trim(adjustl(Folder_with_CS))	! to create a folder use this command
            CALL system(command)	! create the folder
         endif
         !-------------------------------------
         ! Set part of the file name corresponding to the model used for muon Bremsstrahlung CS:
         select case (numpar%Mu_Brems)	! elastic scattering: 0=excluded, 1=BHW
         case (1)	! BHW
            Model_name = 'BHW'
         case default	! exclude
            Model_name = 'NO'
         end select
         
         ! For muon CSs use the same grid as for photons:
         Ngrid = size(Element%Phot_absorption%E)
         allocate(Element%Muon_Brems%E(Ngrid))
         Element%Muon_Brems%E = Element%Phot_absorption%E
         allocate(Element%Muon_Brems%Total(Ngrid))
         allocate(Element%Muon_Brems%Total_MFP(Ngrid))
         
         ! Calculate total elastic cross sections for this element:
         File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_muon_Brems_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!          if (file_exist) then	! only create it if file does not exist
         if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
            open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
            count_lines = 0
            do m = 1, Ngrid	! for all energy grid points:
               read(FN,'(es,es,es)', IOSTAT=Reason) Element%Muon_Brems%E(m), Element%Muon_Brems%Total(m), Element%Muon_Brems%Total_MFP(m)
               call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
               if (.not.read_well) then
                  close(FN)	! redo the file
                  goto 8392
               endif
            enddo ! m = 1, Ngrid
         else 
8392     open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
            do m = 1, Ngrid	! for all energy grid points:
               call get_muon_Brems_CS(Element%Muon_Brems%E(m), Element, numpar, Element%Muon_Brems%Total(m))	! below
               Element%Muon_Brems%Total_MFP(m) = MFP_from_sigma(Element%Muon_Brems%Total(m),  1.0d24)    ! module "CS_general_tools"
               write(FN,'(es,es,es)') Element%Muon_Brems%E(m), Element%Muon_Brems%Total(m), Element%Muon_Brems%Total_MFP(m)
            enddo ! m = 1, Ngrid
         endif
         close(FN)
         ! Normalize MFPs and Se for real material density:
         ! Take into account the density of atoms of this particular kind:
         elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         Element%Muon_Brems%Total_MFP(:) = Element%Muon_Brems%Total_MFP(:) * 1.0d24/(Material(i)%At_Dens * elem_contrib)
      enddo LMNT
      
      ! And total cross sections:
      ! Where to save:
      Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
      Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
      File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_muon_Brems_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
      ! Get the cross sections to save:
      Nsiz = size(Material(i)%Ph_absorption_total%E)   ! grid for total cross section is different from core-shells grid
      allocate(Material(i)%Muon_Brems_total%E(Nsiz))
      Material(i)%Muon_Brems_total%E = Material(i)%Ph_absorption_total%E
      allocate(Material(i)%Muon_Brems_total%Total(Nsiz))
      allocate(Material(i)%Muon_Brems_total%Total_MFP(Nsiz))
      Material(i)%Muon_Brems_total%Total(:) = 0.0d0
      Material(i)%Muon_Brems_total%Total_MFP(:) = 0.0d0
       ! Sum up CSs from each element:
      do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         do m = 1, Nsiz	! for all energy grid points:
            E => Material(i)%Muon_Brems_total%E(m)    ! photon energy [eV]
            call interpolate_data_single(Element%Muon_Brems%E,  Element%Muon_Brems%Total(:), E, sigma) ! module "Little_subroutines"
            call interpolate_data_single(Element%Muon_Brems%E,  Element%Muon_Brems%Total_MFP(:), E, MFP) ! module "Little_subroutines"
            ! Add them into the arrays:
            Material(i)%Muon_Brems_total%Total(m) = Material(i)%Muon_Brems_total%Total(m) + sigma*elem_contrib
            Material(i)%Muon_Brems_total%Total_MFP(m) = Material(i)%Muon_Brems_total%Total_MFP(m) + 1.0d0/MFP !*elem_contrib ! to be inverted
         enddo ! m
      enddo ! j
      ! Invert to get MFP:
      Material(i)%Muon_Brems_total%Total_MFP(:) = 1.0d0/Material(i)%Muon_Brems_total%Total_MFP(:) ! [A]
      ! Save the total cross section into the file:
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
      if (.not.file_exist .or. numpar%recalculate_MFPs) then	! only create it if file does not exist
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            write(FN,'(es,es,es)') Material(i)%Muon_Brems_total%E(m), Material(i)%Muon_Brems_total%Total(m), Material(i)%Muon_Brems_total%Total_MFP(m)
         enddo ! m = 1, Nsiz
         close(FN)
      endif
      
   enddo TRGT
   
   write(*, '(a)') ' Done.'
!    print*, 'Electron Bremsstrahlung cross sections are obtained.'
   nullify(path_sep, E, Element)
end subroutine get_muon_Brems


subroutine get_muon_Brems_CS(Ee, Element, numpar, sigma, E_max_in)
   real(8), intent(in) :: Ee	! [eV] muon energy
   type(Atom_kind), intent(in), target :: Element	! data for this element
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   real(8), intent(out) :: sigma	! [A^2] cross section of Bremsstrahlung
   real(8), intent(in), optional :: E_max_in    ! upper integration limit [eV]
   real(8) :: Fp
   select case (numpar%Mu_Brems)	! Bremsstrahlung scattering: 0=excluded, 1=BHW
   case (1)	! BHW
      if (present(E_max_in)) then
         sigma = Bremsstrahlung_total_CS(Ee, Element, E_max_in)	! module "CS_electrons_Bremsstrahlung"
      else
         sigma = Bremsstrahlung_total_CS(Ee, Element)	! module "CS_electrons_Bremsstrahlung"
      endif
      ! Add the muon correction to the electron CS, Eq.(58) [2]:
      !Fp = muon_Brems_correction(dble(Element%Zat), Ee, g_me_eV)  ! below
      Fp = muon_Brems_correction(dble(Element%Zat), Ee, (g_M_muon_MeV*1.0d6) )  ! below
      sigma = Fp * sigma
   case default	! exclude
      sigma = 0.0d0
   end select
end subroutine get_muon_Brems_CS


pure function muon_Brems_correction(Z, Ee, mc2) result(Fp)
   real(8) Fp   ! Eq.(3.154) [1]
   real(8), intent(in) :: Z, Ee, mc2
   real(8) :: t, arg, t2, t3, t6
   t = log(1.0d0 + 1.0d6/(Z*Z)*Ee/mc2) ! Eq.(3.155) [1]
   t2 = t*t
   t3 = t2*t
   t6 = t3*t3
   arg = -1.2359d-1*t + 6.1274d-2*t2 - 3.1516d-2*t3 + &
         7.7446d-3*t2*t2 - 1.0595d-3*t2*t3 +7.0568d-5 * t6 - &
         1.8080d-6*t6*t
   Fp = 1.0d0 - exp(arg)  ! Eq.(3.154) [1]
end function muon_Brems_correction



function get_energy_transfer_Bremsstrahlung_muon(Ee, Element, numpar, j, CS_tot) result(dE)
   real(8) :: dE   ! [eV] sampled transferred energy
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Atom_kind), intent(in) :: Element    ! all information about this element
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), intent(in) :: CS_tot    ! [A^2] precalculated total cross section
   integer, intent(in) :: j     ! number or element
   !---------------------------------------
   real(8) :: eps, RN, CS_sampled, CS_cur, E_left, E_right, E_cur
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-3
   
   E_left = m_Wth ! [eV] minimal photon energy
   !E_right = Ee/(Ee + g_me_eV)   ! [eV] maximal transferred energy
   E_right = Ee/(Ee + (g_M_muon_MeV*1.0d6) )   ! [eV] maximal transferred energy

   ! Sample the cross section:
   call random_number(RN)
   CS_sampled = RN*CS_tot
   
   ! Start finding CS:
   E_cur = (E_left + E_right)*0.5d0
!    CS_cur = Bremsstrahlung_total_CS(Ee, Element, E_cur)  ! module "CS_electrons_Bremsstrahlung"
   call get_muon_Brems_CS(Ee, Element, numpar, CS_cur, E_cur)    ! above
   
   ! Search by bisection method:
   do while (ABS(CS_cur - CS_sampled)/CS_sampled > eps)
      if (CS_cur > CS_sampled) then
         E_right = E_cur
      else
         E_left = E_cur
      endif
      E_cur = (E_left + E_right)/2.0d0
      if (abs(E_left - E_right) < eps) exit  ! precise enough
!       CS_cur =  Bremsstrahlung_total_CS(Ee, Element, E_cur) ! module "CS_electrons_Bremsstrahlung"
      call get_muon_Brems_CS(Ee, Element, numpar, CS_cur, E_cur)    ! above
   enddo
   
   ! Output: sampled transferred energy:
   dE = E_cur
end function get_energy_transfer_Bremsstrahlung_muon


end module CS_muon_Bremsstrahlung
