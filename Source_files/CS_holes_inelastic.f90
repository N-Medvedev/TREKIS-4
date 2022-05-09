! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018-2019
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains all cross sections for electrons inelastic scattering 
! References used:
! [1] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2008: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2009)
! [2] F. Salvat, J.M. Fernhdez-Varea, NIMB 63 (1992) 255-269
! [3] I. Plante, F,Cucinotta, New Journal of Physics 11 (2009) 063047
! [4] "Transport of Electrons and Photons"  Edited by Jenkins, Nelson, Rindi
! [5] Kim, Santos, Parente, PRB 62 (2000) 052710
! [6] R. Rymzhanov et al., Nucl. Instrum. and Methods B 388 (2016) 41


module CS_holes_inelastic
use Universal_constants
use CDF_Ritchi
use Relativity
use Objects
use Read_input_data, only: m_input_folder, m_holes_inelast_CS, m_folder_materials
use Dealing_with_files, only: read_file
use CS_general_tools, only: MFP_from_sigma, find_valence_hole_mass
use CS_electrons_inelastic, only: CDF_total_CS_nonrel
use CS_holes_elastic, only: set_grid_valence_holes_CS
use CS_integration_limits, only: W_min, W_max, W_max_nonrel_free
use Little_subroutines, only: interpolate_data_single

implicit none

 contains
 
!-----------------------------------------------------------------------------
! Subroutines for inelastic scattering:
subroutine get_hole_IMFP(Material, numpar, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material	!material parameters of each target that it's constructed of
   type(Num_par), intent(inout), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------------------------
   integer :: i, j, k, m, Ngrid, N_targets, N_elements, N_shells, FN, Reason, count_lines
   real(8) :: eps, sigma, Se, lambda, sigma_cur, Se_cur, lambda_cur
   real(8) :: Total_Se, Total_lambda, N_elem, elem_contrib
   character(200) :: command, File_name, Path_valent
   character(25) :: Model_name, mass_model
   logical :: file_exist, read_well
   !--------------------------------
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   character(10) :: temp_c

   write(*, '(a)') ' Obtaining valence hole inelastic scattering cross sections...'

   path_sep => numpar%path_sep
   N_targets = size(Material)	! that's how many different targets user specified
   
   eps = 1.0d-6 ! margin within which mass equals to zero
   if (numpar%H_m_eff > eps) then ! a constant coefficient * me
      if (abs(numpar%H_m_eff) < 1d6) then
         write(mass_model,'(f6.2)') numpar%H_m_eff
         mass_model = '_meff_'//trim(adjustl(mass_model))//'_'
      else  ! hole is too heave to move
         mass_model = '_frozen_'
         numpar%H_inelast = -1    ! no scattering for a frozen hole
      endif
   elseif (abs(numpar%H_m_eff) < eps) then ! equals to free electron mass
      mass_model = '_free_'
   else  ! effective mass from DOS
      mass_model = '_DOS_'   
   endif
   
   ! Get valence hole MFPs:
   TRGT:do i = 1, N_targets	! for each target
      ! Set part of the file name corresponding to the model used for electron inelastic CS:
      select case (numpar%H_inelast)
      case (1:5)	! CDF with Ritchi's oscillators 
         select case (numpar%CDF_model)
         case default   ! Ritchie
            write(temp_c,'(a)') '_R'    ! Ritchie
         case (1)   ! Mermin
            write(temp_c,'(a)') '_M'    ! Mermin
         endselect
         
         Model_name = 'CDF'//trim(adjustl(temp_c))
      case default	! exclude
         Model_name = 'NO'
      end select

      ! Path to files with the valence band and total MFPs (which are the same, so only use total one!):
      Path_valent = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
      Path_valent = trim(adjustl(Path_valent))//path_sep//trim(adjustl(Material(i)%Name))
         
      ! For valence holes inelastic CSs use the same grid as for elastic ones:
      call set_grid_valence_holes_CS(Material(i)%H_inelastic_total%E, abs(Material(i)%DOS%E_VB_bottom - Material(i)%DOS%E_VB_top) )   ! module "CS_holes_elastic"
      Ngrid = size(Material(i)%H_inelastic_total%E)
      allocate(Material(i)%H_inelastic_total%Total(Ngrid))
      allocate(Material(i)%H_inelastic_total%Total_MFP(Ngrid))
      allocate(Material(i)%H_inelastic_total%Total_Se(Ngrid))

      ! Calculate valence band scattering CS:
      VAL:if (allocated(Material(i)%CDF_valence%A)) then    ! if valence band CDF exists
         ! Check file with CSs:
         File_name = trim(adjustl(Path_valent))//path_sep//trim(adjustl(m_holes_inelast_CS))//'_total_'//trim(adjustl(Model_name))//trim(adjustl(mass_model))//'.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!          if (file_exist) then	! just read from the file, no need to recalculate:
         if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
            open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
            ! Get the MFP and CDF points:
            count_lines = 0
            do m = 1, Ngrid	! for all energy grid points:
               E => Material(i)%H_inelastic_total%E(m)	! electron energy [eV]
               read(FN,'(es,es,es,es)', IOSTAT=Reason) E, Material(i)%H_inelastic_total%Total(m), &
                            Material(i)%H_inelastic_total%Total_MFP(m), Material(i)%H_inelastic_total%Total_Se(m)
              call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
              if (.not.read_well) then
                 close(FN)	! redo the file
                 goto 8391
              endif
              ! Recalculate the MFP since material density might be different from the one saved in the file:
              Material(i)%H_inelastic_total%Total_MFP(m) = MFP_from_sigma(Material(i)%H_inelastic_total%Total(m),  Material(i)%At_Dens) ! [A] module "CS_general_tools"
            enddo
         else	! no such file => create it
8391     open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
            do m = 1, Ngrid	! for all energy grid points:
               E => Material(i)%H_inelastic_total%E(m)	! electron energy [eV]
               call get_h_inelastic_CS(E, Material(i), Material(i)%DOS, numpar, j, k, Material(i)%DOS%Egap, sigma, CDF_dispers=numpar%CDF_dispers, Se=Se)	! see below
               Material(i)%H_inelastic_total%Total(m) = sigma   ! [A^2]
               lambda = MFP_from_sigma(sigma,  Material(i)%At_Dens)    ! module "CS_general_tools"
               Material(i)%H_inelastic_total%Total_MFP(m) = lambda  ! [A]
               Material(i)%H_inelastic_total%Total_Se(m) = Se    ! [eV/A]
               write(FN,'(es,es,es,es)') Material(i)%H_inelastic_total%E(m), sigma, lambda, Se
            enddo ! m = 1, Ngrid
         endif ! (file_exist)
         close(FN)
      else VAL  ! No scattering possible
         Material(i)%H_inelastic_total%Total(:) = 0.0d0            ! [A^2]
         Material(i)%H_inelastic_total%Total_MFP(:) = 1.0d25 ! [A]
         Material(i)%H_inelastic_total%Total_Se(:) = 0.0d0      ! [eV/A]
      endif VAL
   enddo TRGT
   
!    write(*,'(a)') ' Done.'
   
   nullify(path_sep, E, Element)
end subroutine get_hole_IMFP



! Interface to select the model of electron inelastic scattering cross section:
 subroutine get_h_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, sigma, CDF_dispers, Se, Emax, Sigma_sampled, E_sampled)
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   real(8), intent(out) :: sigma	! [A^2] cross section
   real(8), intent(out), optional :: Se   ! [eV/A] energy loss
   integer, intent(in), optional :: CDF_dispers ! dispersion relation
   real(8), intent(in), optional :: Emax  ! [eV] upper integration limit (maximal transfered energy)
   real(8), intent(in), optional :: Sigma_sampled   ! to search for transfered energy, sampled cross section
   real(8), intent(inout), optional :: E_sampled    ! energy corresponding to sampled cross section
   !---------------------------------------
   real(8) :: Zeff,  Se1
   integer :: dispers
   Zeff = 1.0d0	! hole charge
   select case (numpar%H_inelast)
   case (1:5)	! CDF with Ritchie model
      if (present(CDF_dispers)) then	! user provided alternative (e.g. for all but VB)
         dispers = CDF_dispers
      else	! for the VB
         dispers = numpar%CDF_dispers
      endif

      ! Integrals as requested by user:
      if (present(Emax)) then  ! integral up to given limit:
         select case (dispers) ! 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         case (1)  ! plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                   v_f=Material%DOS%v_f, CDF_dispers=dispers, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                   CDF_dispers=dispers, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
         end select
         
      elseif (present(Sigma_sampled) .and. present(E_sampled)) then ! Integral to find sampled energy according to sampled CS:
         select case (dispers) ! 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                   Sigma_sampled=Sigma_sampled, E_sampled=E_sampled)  ! module "CS_electrons_inelastic"
         case (1)  ! plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                   v_f=Material%DOS%v_f, CDF_dispers=dispers, Sigma_sampled=Sigma_sampled, E_sampled=E_sampled)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, CDF_dispers=dispers, &
                   Sigma_sampled=Sigma_sampled, E_sampled=E_sampled)  ! module "CS_electrons_inelastic"
         end select
      
      else ! Total cross section:
         select case (dispers) ! 0=free electron, 1=plasmon-pole, 2=Ritchie
         case default  ! free electron
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model)  ! module "CS_electrons_inelastic"
         case (1)  ! plasmon pole
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                   v_f=Material%DOS%v_f, CDF_dispers=dispers)  ! module "CS_electrons_inelastic"
         case (2)  ! Ritchie
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                   CDF_dispers=dispers)  ! module "CS_electrons_inelastic"
         end select
      endif

   case default	! exclude
      sigma = 0.0d0
      Se1 = 0.0d0
      if (present(E_sampled)) then
         E_sampled = 0.0d0
      endif
   end select
   if (present(Se)) Se = Se1 ! output if required
end subroutine get_h_inelastic_CS



! Interface to select the model of electron inelastic scattering cross section:
function get_inelastic_energy_transfer_h(Ee, Material, DOS, numpar, j, k, Ip, CDF_dispers, CDF_m_eff, hw_phonon) result(dE)
   real(8) :: dE   ! [eV] sampled transferred energy
   real(8), intent(in) :: Ee	! [eV] valence hole kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   integer, intent(in), optional :: CDF_dispers, CDF_m_eff	! dispersion relation and effective mass
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency, for particle-phonon scattering
   !---------------------------------------
   real(8) :: eps, Zeff, max_E0, Eeq, RN, CS_sampled, CS_cur, E_left, E_right, E_cur, M_h, Se1, Mc2, mtc2
   real(8) :: CS_tot    ! [A^2] precalculated total cross section
   integer :: dispers, m_eff, i
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-3
   
   ! Set the model parameters:
   Zeff = 1.0d0 ! hole charge
   M_h = find_valence_hole_mass(numpar, DOS, Ee)     ! module "CS_general_tools"
   ! Check just in case some problems with discret DOS:
   if (M_h <= 0.0d0) then   ! use free electron mass
      M_h = g_me
   endif
   
   if (present(CDF_dispers)) then	! user provided alternative (e.g. for all but VB)
      dispers = CDF_dispers
   else	! for the VB
      dispers = numpar%CDF_dispers
   endif
   if (present(CDF_m_eff)) then	! user provided alternative (e.g. for all but VB)
      m_eff = CDF_m_eff
   else
      if (numpar%CDF_m_eff <= 0) then
         m_eff = 1
      else
         m_eff = numpar%CDF_m_eff
      endif
   endif

   Mc2 = rest_energy(M_h)   ! module "Relativity"
   mtc2 = rest_energy(m_eff*g_me)   ! module "Relativity"

   E_left = W_min(Ip, Mc2, mtc2, Ee, Ip) ! module "CS_integration_limits"
   ! Upper integration limit [eV]:
   if (present(hw_phonon)) then
      E_right = W_max_nonrel_free(Mc2, mtc2, .false., Ee, Ip, hw_phonon)   ! module "CS_integration_limits"
   else
      E_right = W_max_nonrel_free(Mc2, mtc2, .false., Ee, Ip)   ! module "CS_integration_limits"
   endif

   ! Sample the cross section:
   call random_number(RN)
   
!    do i = 1, 100    ! testing
!    RN = dble(i)/100.0   ! testing
!    print*, 'E:', E_left, E_right
   
   
   ! Get total cross sections (recalculate):
!    call get_h_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, CS_tot)   ! above
!    print*, 'Recalc:', Ee, CS_tot
   ! Getting precalculated cross section (linear interpolation may be too crude approximation!):
   call interpolate_data_single(Material%H_inelastic_total%E, Material%H_inelastic_total%Total, Ee, CS_tot) ! module "Little_subroutines"
!   print*, 'Pre:', Ee, CS_tot
   
!    CS_tot = get_integral_inelastic_CS_h(Ee, Material, DOS, numpar, j, k, Ip, 4, Zeff, M_h, dispers, m_eff, E_right)   ! below
   CS_sampled = RN*CS_tot
   
   if ( (RN < 1.0d-10) .or. (CS_sampled < 1.0d-12) )then   ! no need to sample, it's just the lower limit
      E_cur = E_left
   else ! sample it
      ! Find the sampled CS:
      call get_h_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, CS_cur, Sigma_sampled=CS_sampled, E_sampled = E_cur)   ! above 
   endif

!    write(*,'(f,es,f,f,f)') Ee, E_cur, CS_sampled, RN, CS_tot
!    enddo ! testing
!    pause 'get_inelastic_energy_transfer_h'
   
   ! Output: sampled transferred energy:
   dE = E_cur
end function get_inelastic_energy_transfer_h



function get_integral_inelastic_CS_h(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, M_h, dispers, m_eff, Emax) result (sigma)
   real(8) :: sigma	! [A^2] cross section
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Density_of_states), intent(in):: DOS	! DOS of the material
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   integer, intent(in) :: El_inelast    ! model for the inelastic cross section
   real(8), intent(in) :: Zeff  ! effective charge of the incident particle
   real(8), intent(in) :: M_h   ! [kg] incident particle mass
   integer, intent(in) :: dispers, m_eff	! dispersion relation and effective mass
   real(8), intent(in) :: Emax  ! [eV] upper integration limit
   !---------------------------------------
   real(8) :: Se1
   integer :: ish

   ! Now chose the model of CS:
   select case (El_inelast)  ! chose which model for electron inelastic cross section to use
   case (4) ! nonrelativistic Ritchie CDF, holes can only ionize VB, not core shells
      select case (numpar%CDF_dispers)  ! dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
!       case default  ! free electron
!          call CDF_total_CS_nonrel(sigma, Se1, Ee, M_h, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, &
!                  DOS%k, DOS%Eff_m, .false., 1.0d0, Material%At_Dens, Material%DOS, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
!       case (1)  ! Plasmon pole
!          call CDF_total_CS_nonrel(sigma, Se1, Ee, M_h, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
!                   .false., 1.0d0, Material%At_Dens, Material%DOS, CDF_dispers=numpar%CDF_dispers, &
!                   v_f=Material%DOS%v_f, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
!       case (2)  ! Ritchie extended
!          call CDF_total_CS_nonrel(sigma, Se1, Ee, M_h, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
!                   .false., 1.0d0, Material%At_Dens, Material%DOS, CDF_dispers=numpar%CDF_dispers, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
!       end select ! (numpar%CDF_dispers)

      case default  ! free electron
         call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
      case (1)  ! plasmon pole
         call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, v_f=Material%DOS%v_f, &
                   CDF_dispers=numpar%CDF_dispers, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
      case (2)  ! Ritchie
         call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, Ip, Material%T_eV, Material%CDF_valence, g_me, DOS%k, DOS%Eff_m, &
                   .false., numpar%H_m_eff, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                   CDF_dispers=numpar%CDF_dispers, Wmax_in=Emax)  ! module "CS_electrons_inelastic"
      end select
   case default	! exclude
      sigma = 0.0d0
   end select
end function get_integral_inelastic_CS_h


end module CS_holes_inelastic
