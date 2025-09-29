! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains all cross sections for valence holes elastic scattering
! Since valence holes are always non-relativistic, one can reuse subroutines from TREKIS-3


module CS_holes_elastic
use Universal_constants
use CDF_Ritchi, only: m_two_third
use CDF_delta, only: CDF_total_CS_delta, energy_loss_delta
use Relativity
use Objects
use Little_subroutines, only: interpolate_data_single
use Read_input_data, only: m_input_folder, m_folder_materials, m_holes_elast_CS, m_hole_CS
use Dealing_with_files, only: read_file
use CS_general_tools, only: MFP_from_sigma, find_valence_hole_mass
use CS_electrons_elastic, only: Mott_total_CS
use CS_electrons_inelastic,only: CDF_total_CS_nonrel
use CS_integration_limits, only: find_Wmax_equal_Wmin
use SHI_charge_state, only: Equilibrium_charge_SHI


implicit none

 contains

!-----------------------------------------------------------------------------
! Subroutines for elastic scattering:
subroutine get_holes_EMFP(Material, numpar, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material  !material parameters of each target that it's constructed of
   type(Num_par), intent(inout), target :: numpar  ! all numerical parameters
   type(Error_handling), intent(inout) :: Err   ! error log
   !--------------------------------
   integer :: i, j, k, m, Ngrid, N_targets, N_elements, N_shells, FN, Reason, count_lines, Nsiz
   real(8) :: sigma, MFP, N_elem, elem_contrib, eps
   character(200) :: Path, Folder_with_CS, command, File_name, Model_name, Path_total
   character(25) :: mass_model, temp_c, temp_c2, temp_c3, temp_c4
   logical :: file_exist, read_well
   !--------------------------------
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   
   eps = 1.0d-6 ! margin for considering effective mass equal to free electron mass
   
   !write(*, '(a)') ' Obtaining valence hole elastic scattering cross sections...'
   write(*, '(a)', advance='no') ' Obtaining velence hole elastic scattering cross sections...'
   
   path_sep => numpar%path_sep
   Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_hole_CS))	! Electron CSs are storred in this folder
   
   if (numpar%H_m_eff > eps) then ! a constant coefficient * me
      if (abs(numpar%H_m_eff) < 1d6) then
         write(mass_model,'(f6.2)') numpar%H_m_eff
         mass_model = '_meff_'//trim(adjustl(mass_model))//'_'
      else  ! hole is too heave to move
         mass_model = '_frozen_'
         numpar%H_elast = -1    ! no scattering for a frozen hole
      endif
   elseif (abs(numpar%H_m_eff) < eps) then ! equals to free electron mass
      mass_model = '_free_'
   else  ! effective mass from DOS
      mass_model = '_DOS_'   
   endif
   
   N_targets = size(Material)	! that's how many different targets user specified
   
   ! Get electron MFPs:
   TRGT:do i = 1, N_targets	! for each target
      ! Set total grid:
      call set_grid_valence_holes_CS(Material(i)%H_elastic_total%E, abs(Material(i)%DOS%E_VB_bottom - Material(i)%DOS%E_VB_top) )   ! below
      Nsiz = size(Material(i)%H_elastic_total%E)    ! grid for total cross section
      allocate(Material(i)%H_elastic_total%Total(Nsiz))
      allocate(Material(i)%H_elastic_total%Total_MFP(Nsiz))
      Material(i)%H_elastic_total%Total(:) = 0.0d0
      Material(i)%H_elastic_total%Total_MFP(:) = 0.0d0
   
   
      ! Set part of the file name corresponding to the model used for electron elastic CS:
      select case (numpar%H_elast)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
      case (1,5)	! CDF - only total cross section
         ! Construct the model name:
         write(temp_c,'(f8.1)') Material(i)%T    ! target temperature
         write(temp_c2,'(f8.1)') numpar%CDF_Eeq_elast
         if (numpar%CDF_elast_Zeff /= 1) then
            write(temp_c3,'(a)') '_Zeff'
         else
            write(temp_c3,'(a)') '_Z1'
         endif

         select case (numpar%CDF_model)
         case default   ! Ritchie
            write(temp_c4,'(a)') '_R_'    ! Ritchie
         case (1)   ! Mermin
            write(temp_c4,'(a)') '_M_'    ! Mermin
         case (5)   ! Single pole
            write(temp_c4,'(a)') '_SP_'    ! Single-pole
         endselect

         Model_name = 'CDF_delta_T_'//trim(adjustl(temp_c))//'K'//trim(adjustl(temp_c4))//'Eeq_'//trim(adjustl(temp_c2))//trim(adjustl(temp_c3))
!          Model_name = 'CDF_delta_T_'//trim(adjustl(temp_c))//'K_Eeq_'//trim(adjustl(temp_c2))
          
         ! Check file with CSs:
         ! Where to save:
         Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
         Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
         File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_holes_elast_CS))//'_total_'//trim(adjustl(Model_name))//trim(adjustl(mass_model))//'.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!          if (file_exist) then	! just read from the file, no need to recalculate:
         if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
            open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
            ! Get the MFP and CDF points:
            count_lines = 0
            do m = 1, Nsiz  ! for all energy grid points:
               read(FN,'(es,es,es)', IOSTAT=Reason) Material(i)%H_elastic_total%E(m), Material(i)%H_elastic_total%Total(m), Material(i)%H_elastic_total%Total_MFP(m)
               call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
               if (.not.read_well) then
                  close(FN)	! redo the file
                  goto 8392
               endif
            enddo
         else	! no such file => create it
8392     open(newunit = FN, FILE = trim(adjustl(File_name)), action='write')
            Element => Material(i)%Elements(1)  ! unused here, so just set an "empty" value
            do m = 1, Nsiz  ! for all energy grid points:
               call get_h_elastic_CS(Material(i)%H_elastic_total%E(m), Material(i), Element, numpar, Material(i)%H_elastic_total%Total(m))   ! [A^2] below
               Material(i)%H_elastic_total%Total_MFP(m) = MFP_from_sigma(Material(i)%H_elastic_total%Total(m), Material(i)%At_Dens)    ! module "CS_general_tools"
!                Material(i)%El_inelastic_valent%Total_Se(m) = Se    ! [eV/A]
               write(FN,'(es,es,es)') Material(i)%H_elastic_total%E(m), Material(i)%H_elastic_total%Total(m), Material(i)%H_elastic_total%Total_MFP(m)
!                write(6,'(es,es,es)') Material(i)%H_elastic_total%E(m), Material(i)%H_elastic_total%Total(m), Material(i)%H_elastic_total%Total_MFP(m)
            enddo ! m = 1, Nsiz
         endif ! (file_exist)
         close(FN)
         
      case (2)	! Mott - cross sections for each type of atom
         Model_name = 'Mott'
         N_elem = dble(SUM(Material(i)%Elements(:)%percentage))	! number of elements in this compound material
         N_elements = size(Material(i)%Elements)	! that's how many different elements are in this target
         LMNT:do j =1, N_elements	! for each element
            Element => Material(i)%Elements(j)	! all information about this element
            ! Check if folder for this element already exists:
!             Folder_with_CS = trim(adjustl(Path))//path_sep//trim(adjustl(Element%Name))
!             inquire(DIRECTORY=trim(adjustl(Folder_with_CS)),exist=file_exist)    ! check if input file excists
!             if (.not.file_exist) then	! to make sure that such a folder is present (even if empty)
!                ! Create a new directory for output files:
!                command='mkdir '//trim(adjustl(Folder_with_CS))	! to create a folder use this command
!                CALL system(command)	! create the folder
!             endif
            !-------------------------------------
      
            ! For valence hole CSs use specific grid:
            call set_grid_valence_holes_CS(Element%H_elastic%E, abs(Material(i)%DOS%E_VB_bottom - Material(i)%DOS%E_VB_top) )   ! below
            Ngrid = size(Element%H_elastic%E)
            allocate(Element%H_elastic%Total(Ngrid))
            allocate(Element%H_elastic%Total_MFP(Ngrid))
            
            ! Where to save:
            Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
            Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
            ! Calculate elastic cross sections for this element:
            File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_holes_elast_CS))//'_'//trim(adjustl(Element%Name))//'_'// &
                                trim(adjustl(Model_name))//trim(adjustl(mass_model))//'.dat'
            inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!             if (file_exist) then	! only create it if file does not exist
            if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
               open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
               count_lines = 0
               do m = 1, Ngrid	! for all energy grid points:
                  read(FN,'(es,es,es)', IOSTAT=Reason) Element%H_elastic%E(m), Element%H_elastic%Total(m), Element%H_elastic%Total_MFP(m)
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not.read_well) then
                     close(FN)	! redo the file
                     goto 8391
                  endif
               enddo ! m = 1, Ngrid
            else
8391        open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
               do m = 1, Ngrid	! for all energy grid points:
                  call get_h_elastic_CS(Element%H_elastic%E(m), Material(i), Element, numpar, Element%H_elastic%Total(m))   ! below
                  Element%H_elastic%Total_MFP(m) = MFP_from_sigma(Element%H_elastic%Total(m),  1.0d24)    ! module "CS_general_tools"
                  write(FN,'(es,es,es)') Element%H_elastic%E(m), Element%H_elastic%Total(m), Element%H_elastic%Total_MFP(m)
               enddo ! m = 1, Ngrid
            endif
            close(FN)
            ! Normalize to material density:
            ! Normalize to the material density:
            elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
            Element%H_elastic%Total_MFP(:) = Element%H_elastic%Total_MFP(:) * 1.0d24/(Material(i)%At_Dens * elem_contrib)
         enddo LMNT

         ! Sum up CSs from each element:
         do j =1, N_elements	! for each element
            Element => Material(i)%Elements(j)	! all information about this element
            elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
            do m = 1, Nsiz  ! for all energy grid points:
               E => Material(i)%H_elastic_total%E(m)    ! photon energy [eV]
               call interpolate_data_single(Element%H_elastic%E,  Element%H_elastic%Total(:), E, sigma) ! module "Little_subroutines"
               call interpolate_data_single(Element%H_elastic%E,  Element%H_elastic%Total_MFP(:), E, MFP) ! module "Little_subroutines"
               ! Add them into the arrays:
               Material(i)%H_elastic_total%Total(m) = Material(i)%H_elastic_total%Total(m) + sigma*elem_contrib
               Material(i)%H_elastic_total%Total_MFP(m) = Material(i)%H_elastic_total%Total_MFP(m) + 1.0d0/MFP !*elem_contrib ! to be inverted
            enddo ! m
         enddo ! j
         ! Invert to get MFP:
         Material(i)%H_elastic_total%Total_MFP(:) = 1.0d0/Material(i)%H_elastic_total%Total_MFP(:) ! [A]
         
      case (3)	! DSF - total cross section
         Model_name = 'DSF'
      case default	! exclude
         Model_name = 'NO'
         ! All cross sections are zero:
         do m = 1, Nsiz ! for all energy grid points:
            Material(i)%H_elastic_total%Total(m) = 0.0d0
            Material(i)%H_elastic_total%Total_MFP(m) = 1.0d25
         enddo ! m
      end select

      ! And total cross sections:
      ! Where to save:
      Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
      Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
      File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_holes_elast_CS))//'_total_'//trim(adjustl(Model_name))//trim(adjustl(mass_model))//'.dat'
      
      ! Save the total cross section into the file:
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
      if (.not.file_exist .or. numpar%recalculate_MFPs) then	! only create it if file does not exist or user wants it recalculated
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            write(FN,'(es,es,es)') Material(i)%H_elastic_total%E(m), Material(i)%H_elastic_total%Total(m), Material(i)%H_elastic_total%Total_MFP(m)
         enddo ! m = 1, Nsiz
         close(FN)
      endif
      
   enddo TRGT

   write(*,'(a)') ' Done.'

   nullify(path_sep, E, Element)
end subroutine get_holes_EMFP


pure subroutine set_grid_valence_holes_CS(Egrid, E_VB)
   real(8), dimension(:), allocatable, intent(inout) :: Egrid  ! valence holes energy grid to be defined
   real(8), intent(in) :: E_VB ! [eV] width of the valence band
   real(8) :: dE
   integer :: Nsiz, i
   dE = 0.1d0   ! [eV] grid step size
   Nsiz = CEILING(E_VB/dE)    ! size of the grid
   allocate(Egrid(Nsiz))
   if (Nsiz >= 1) then  ! only if this array is non-empty
      Egrid(1) = 0.0d0 ! starting point
      do i = 2, Nsiz   ! for all grid points
         Egrid(i) = Egrid(i-1) + dE
      enddo
   endif
end subroutine set_grid_valence_holes_CS


! Interface to select the model of electron elastic scattering cross section:
subroutine get_h_elastic_CS(Ee, Material, Element, numpar, sigma, E_max, Se)
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Atom_kind), intent(in) :: Element	! data for this element
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), intent(out) :: sigma	! [A^2] cross section
   real(8), intent(in), optional :: E_max   ! [eV], integration limit
   real(8), intent(out), optional :: Se ! [eV/A]
   !-----------------------------------------
   real(8) :: Zeff, mass, eps, Eeq, Se1, max_E0, hw_ph_max
   integer :: H_elast
   
   eps = 1.0d-6
   
   select case (numpar%H_elast)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   case (1,5) ! CDF
      ! Find the mass of the valence hole:
      mass = find_valence_hole_mass(numpar, Material%DOS, Ee ) ! module "MC_general_tools"

      ! Target mean atomic number:
      if (numpar%CDF_elast_Zeff /= 1) then
         Zeff = 1.0d0 + Equilibrium_charge_SHI(Ee, mass, Material%Mean_Z, (Material%Mean_Z-1.0d0), 0, 1.0d0) ! module "SHI_charge_state"
      else
         Zeff = 1.0d0    ! electron charge
      endif
      
      
      ! Maximal phonon frequency is defined by the maximal phononic CDF peak:
      hw_ph_max = maxval( Material%CDF_phonon%E0(:) + Material%CDF_phonon%Gamma(:) )
      !hw_ph_max = maxval( Material%CDF_phonon%E0(:) + 10.0d0*Material%CDF_phonon%Gamma(:) )
      
      H_elast = 4   ! for holes, always use Ritchie CDF
!       H_elast = numpar%H_elast   ! to chose the model of CS calculations below
!       ! In case it is delta-CDF model, we have to make sure the energy is not below its applicability limit:
!       if (H_elast == 1) then
!          max_E0 = maxval(Material%CDF_phonon%E0(:))
!          call find_Wmax_equal_Wmin(mass, Material%Mean_Mass, .false., Ee, 1.0d-6, max_E0, Eeq, hw_ph_max)   ! module "CS_integration_limits"
!          ! Check if delta-functional CDF works here,  and apply for electrons above chosen 100 eV:
!          !if (Ee < max(numpar%CDF_Eeq_elast*Eeq, 100.0d0)) then   ! switch to nonrelativistic numerically-integrable CDF:
!          if (Ee < numpar%CDF_Eeq_elast*Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
!              H_elast = 4
!          endif
!       endif ! (H_elast == 1)
      
      select case (H_elast)  ! chose which model for electron inelastic cross section to use
      case (1,5)  ! Delta relativistic CDF
         if (present (E_max)) then
            sigma = CDF_total_CS_delta(H_elast, Ee, mass, Zeff, 1.0d-8, Material%At_Dens, Material%CDF_phonon, &
                                                     Material%Mean_Mass, .false., hw_phonon=hw_ph_max, Emax_in=E_max)    ! module "CDF_delta"
         else
            sigma = CDF_total_CS_delta(H_elast, Ee, mass, Zeff, 1.0d-8, Material%At_Dens, Material%CDF_phonon, &
                                                        Material%Mean_Mass, .false., hw_phonon=hw_ph_max)    ! module "CDF_delta"
         endif
         if (present(Se)) Se = energy_loss_delta(H_elast, Ee, mass, 1.0d0, Zeff, 1.0d-8, Material%At_Dens, Material%Mean_Mass, &
                                          Material%CDF_phonon, .false., 0, hw_ph_max) ! module "CDF_delta"
      case (4) ! nonrelativistic Ritchie CDF
         if (present (E_max)) then
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, mass, Zeff, 1.0d-10, Material%T_eV, Material%CDF_phonon, Material%Mean_Mass, &
                        Material%DOS%k, Material%DOS%Eff_m, .false., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                        hw_phonon = hw_ph_max, Wmax_in=E_max)  ! module "CS_electrons_inelastic"
         else
               call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, mass, Zeff, 1.0d-10, Material%T_eV, Material%CDF_phonon, Material%Mean_Mass, &
                        Material%DOS%k, Material%DOS%Eff_m, .false., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                        hw_phonon = hw_ph_max)  ! module "CS_electrons_inelastic"
         endif
         if (present(Se)) Se = Se1
      end select

!       print*, 'get_h_elastic_CS', Ee, sigma, mass
      
   case (2) ! Mott
      if (numpar%H_m_eff > eps) then ! a constant coefficient * me
         sigma = Mott_total_CS(Ee, Element%Zat, numpar%H_m_eff)     ! module "CS_electrons_elastic"
      elseif (abs(numpar%H_m_eff) < eps) then ! equals to free electron mass
         sigma = Mott_total_CS(Ee, Element%Zat)     ! module "CS_electrons_elastic"
      else  ! effective mass from DOS
         call interpolate_data_single(Material%DOS%E, Material%DOS%Eff_m, (Material%DOS%E_VB_top - Ee), mass) ! module "Little_subroutines"
!          print*, 'get_h_elastic_CS', Ee, mass
         sigma = Mott_total_CS(Ee, Element%Zat, mass=mass)       ! module "CS_electrons_elastic"
      endif
      
   case (3) ! DSF
      ! NOT READY YET
      
   case default ! exclude
      sigma = 0.0d0
      if (present(Se)) Se = 0.0d0
      
   end select
end subroutine get_h_elastic_CS



end module CS_holes_elastic
