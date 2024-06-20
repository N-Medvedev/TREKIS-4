! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018-2020
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains all cross sections for electrons elastic scattering
! References used:
! [1] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2008: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2009)
! [2] F. Salvat, J.M. Fernhdez-Varea, NIMB 63 (1992) 255-269
! [3] I. Plante, F,Cucinotta, New Journal of Physics 11 (2009) 063047
! [4] "Transport of Electrons and Photons"  Edited by Jenkins, Nelson, Rindi
! [5] Kim, Santos, Parente, PRB 62 (2000) 052710


module CS_electrons_elastic
use Universal_constants
use CDF_Ritchi, only: m_two_third
use Relativity
use Objects
use Read_input_data, only: m_electron_CS, m_input_folder, m_electron_elast_CS, m_folder_materials
use Dealing_with_files, only: read_file
use CS_integration_limits, only: find_Wmax_equal_Wmin
use CS_general_tools, only: MFP_from_sigma
use CS_electrons_inelastic, only: CDF_total_CS_nonrel
use CDF_delta, only: CDF_total_CS_delta, energy_loss_delta
use Little_subroutines, only: interpolate_data_single, create_energy_grid
use SHI_charge_state, only: Equilibrium_charge_SHI


implicit none

 contains
 
!-----------------------------------------------------------------------------
! Subroutines for elastic scattering:
subroutine get_electron_EMFP(Material, numpar, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material	!material parameters of each target that it's constructed of
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------------------------
   integer :: i, j, k, m, Ngrid, N_targets, N_elements, N_shells, FN, Reason, count_lines, Nsiz
   real(8) :: sigma, MFP, N_elem, elem_contrib, grid_renorm
   character(500) :: Path, Folder_with_CS, command, File_name, Model_name, Path_total
   character(10) :: temp_c, temp_c2, temp_c3, temp_c4
   logical :: file_exist, read_well
   !--------------------------------
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   
   ! That's the factor to use to set the grid, with respect to the photon grid
   grid_renorm = 10.0d0
   
   write(*, '(a)', advance='no') ' Obtaining electron elastic scattering cross sections...'
   
   path_sep => numpar%path_sep
   Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_electron_CS))	! Electron CSs are storred in this folder
   
   N_targets = size(Material)	! that's how many different targets user specified

   ! Get electron MFPs:
   TRGT:do i = 1, N_targets	! for each target
   
      ! Allocate the cross sections to save:
      Nsiz = size(Material(i)%Ph_absorption_total%E)  ! to find the highest energy
      call create_energy_grid(0.01d0, Material(i)%Ph_absorption_total%E(Nsiz), Material(i)%El_elastic_total%E)    ! module "Little_subroutines"
      Nsiz = size(Material(i)%El_elastic_total%E)  ! grid for total cross section is different from core-shells grid
!       allocate(Material(i)%El_elastic_total%E(Nsiz))
      Material(i)%El_elastic_total%E = Material(i)%El_elastic_total%E / grid_renorm    ! smaller grid to account for phonons
      allocate(Material(i)%El_elastic_total%Total(Nsiz))
      allocate(Material(i)%El_elastic_total%Total_MFP(Nsiz))
      Material(i)%El_elastic_total%Total(:) = 0.0d0
      Material(i)%El_elastic_total%Total_MFP(:) = 0.0d0

      !-------------------------------------
      ! Set part of the file name corresponding to the model used for electron elastic CS:
      select case (numpar%El_elastic)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
      case (1,5)  ! CDF
         ! Construct the model name:
         write(temp_c,'(f8.1)') Material(i)%T    ! target temperature
         write(temp_c2,'(f8.1)') numpar%CDF_Eeq_elast
         select case (numpar%CDF_model)
         case default   ! Ritchie
            write(temp_c4,'(a)') '_R_'    ! Ritchie
         case (1)   ! Mermin
            write(temp_c4,'(a)') '_M_'    ! Mermin
         case (5)   ! Single pole
            write(temp_c4,'(a)') '_SP_'    ! Single-pole
         endselect
         if (numpar%CDF_elast_Zeff /= 1) then
            write(temp_c3,'(a)') '_Zeff'
         else
            write(temp_c3,'(a)') '_Z1'
         endif
         
         Model_name = 'CDF_delta_T_'//trim(adjustl(temp_c))//'K'//trim(adjustl(temp_c4))//'Eeq_'//trim(adjustl(temp_c2))//trim(adjustl(temp_c3))

         ! Check file with CSs:
         ! Where to save:
         Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
         Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
         File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_electron_elast_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
         !if (file_exist) then  ! just read from the file, no need to recalculate:
         if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
            open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
            ! Get the MFP and CDF points:
            count_lines = 0
            do m = 1, Nsiz  ! for all energy grid points:
               read(FN,'(es,es,es)', IOSTAT=Reason) Material(i)%El_elastic_total%E(m), Material(i)%El_elastic_total%Total(m), Material(i)%El_elastic_total%Total_MFP(m)
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
               call get_el_elastic_CS(Material(i)%El_elastic_total%E(m), Material(i), Element, numpar, Material(i)%El_elastic_total%Total(m))   ! [A^2] below
               Material(i)%El_elastic_total%Total_MFP(m) = MFP_from_sigma(Material(i)%El_elastic_total%Total(m), Material(i)%At_Dens)    ! module "CS_general_tools"
!                Material(i)%El_inelastic_valent%Total_Se(m) = Se    ! [eV/A]
               write(FN,'(es,es,es)') Material(i)%El_elastic_total%E(m), Material(i)%El_elastic_total%Total(m), Material(i)%El_elastic_total%Total_MFP(m)
!                write(6,'(es,es,es)') Material(i)%El_elastic_total%E(m), Material(i)%El_elastic_total%Total(m), Material(i)%El_elastic_total%Total_MFP(m)
            enddo ! m = 1, Nsiz
         endif ! (file_exist)
         close(FN)

      case (2)	! Mott
         Model_name = 'Mott'
         
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
         
            ! For electron CSs use the same grid as for photons:
            Ngrid = size(Material(i)%El_elastic_total%E)
            allocate(Element%El_elastic%E(Ngrid))
            Element%El_elastic%E = Material(i)%El_elastic_total%E
            allocate(Element%El_elastic%Total(Ngrid))
            allocate(Element%El_elastic%Total_MFP(Ngrid))
         
            ! Calculate total elastic cross sections for this element:
            File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_electron_elast_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
            inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!             if (file_exist) then	! only create it if file does not exist
            if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
               open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
               count_lines = 0
               do m = 1, Ngrid	! for all energy grid points:
                  read(FN,'(es,es,es)', IOSTAT=Reason) Element%El_elastic%E(m), Element%El_elastic%Total(m), Element%El_elastic%Total_MFP(m)
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not.read_well) then
                     close(FN)	! redo the file
                     goto 8391
                  endif
               enddo ! m = 1, Ngrid
            else
8391        open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
               do m = 1, Ngrid	! for all energy grid points:
                  call get_el_elastic_CS(Element%El_elastic%E(m), Material(i), Element, numpar, Element%El_elastic%Total(m))	! below
                  Element%El_elastic%Total_MFP(m) = MFP_from_sigma(Element%El_elastic%Total(m),  1.0d24)    ! module "CS_general_tools"
                  write(FN,'(es,es,es)') Element%El_elastic%E(m), Element%El_elastic%Total(m), Element%El_elastic%Total_MFP(m)
               enddo ! m = 1, Ngrid
            endif
            close(FN)
            ! Normalize to the material density:
            elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
            Element%El_elastic%Total_MFP(:) = Element%El_elastic%Total_MFP(:) * 1.0d24/(Material(i)%At_Dens * elem_contrib)
         enddo LMNT
      
         ! And total cross sections:
         ! Where to save:
         Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
         Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
         File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_electron_elast_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
         ! Sum up CSs from each element:
         do j =1, N_elements	! for each element
            Element => Material(i)%Elements(j)	! all information about this element
            elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
            do m = 1, Nsiz	! for all energy grid points:
               E => Material(i)%El_elastic_total%E(m)    ! photon energy [eV]
               call interpolate_data_single(Element%El_elastic%E,  Element%El_elastic%Total(:), E, sigma) ! module "Little_subroutines"
               call interpolate_data_single(Element%El_elastic%E,  Element%El_elastic%Total_MFP(:), E, MFP) ! module "Little_subroutines"
               ! Add them into the arrays:
               Material(i)%El_elastic_total%Total(m) = Material(i)%El_elastic_total%Total(m) + sigma*elem_contrib
               Material(i)%El_elastic_total%Total_MFP(m) = Material(i)%El_elastic_total%Total_MFP(m) + 1.0d0/MFP !*elem_contrib ! to be inverted
            enddo ! m
         enddo ! j
         ! Invert to get MFP:
         Material(i)%El_elastic_total%Total_MFP(:) = 1.0d0/Material(i)%El_elastic_total%Total_MFP(:) ! [A]
      
      case (3)	! DSF
         Model_name = 'DSF'
         ! NOT READY YET
      case default	! exclude
         Model_name = 'NO'
         Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
         Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
         File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_electron_elast_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
         ! All cross sections are zero:
         do m = 1, Nsiz	! for all energy grid points:
            Material(i)%El_elastic_total%Total(m) = 0.0d0
            Material(i)%El_elastic_total%Total_MFP(m) = 1.0d25
         enddo ! m
      end select
      
      ! Save the total cross section into the file:
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
      if (.not.file_exist .or. numpar%recalculate_MFPs) then	! only create it if file does not exist or user wants to recalculate it
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            write(FN,'(es,es,es)') Material(i)%El_elastic_total%E(m), Material(i)%El_elastic_total%Total(m), Material(i)%El_elastic_total%Total_MFP(m)
         enddo ! m = 1, Nsiz
         close(FN)
      endif
      
   enddo TRGT
   
   write(*, '(a)') ' Done.'
!    print*, 'Electron elastic scattering cross sections are obtained.'
   
   nullify(path_sep, E, Element)
end subroutine get_electron_EMFP



! Interface to select the model of electron elastic scattering cross section:
subroutine get_el_elastic_CS(Ee, Material, Element, numpar, sigma, mu_max_in, E_max, Se)
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Target_atoms), intent(in) :: Material	!material parameters of each target that it's constructed of
   type(Atom_kind), intent(in) :: Element	! data for this element
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), intent(out) :: sigma	! [A^2] cross section
   real(8), intent(in), optional :: mu_max_in  ! mu=cos(theta), integration limit
   real(8), intent(in), optional :: E_max   ! [eV], integration limit
   real(8), intent(out), optional :: Se ! [eV/A]
   !-----------------------------------------
   real(8) :: Zeff, Se1, max_E0, Eeq, hw_ph_max
   integer :: dispers, m_eff, El_elast
   
   select case (numpar%El_elastic)  ! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF, 4=nonrelativistic, 5=single-pole
   case (1,5) ! CDF
      ! Target mean atomic number:
      if (numpar%CDF_elast_Zeff /= 1) then
         Zeff = 1.0d0 + Equilibrium_charge_SHI(Ee, g_me, Material%Mean_Z, (Material%Mean_Z-1.0d0), 0, 1.0d0) ! module "SHI_charge_state"
      else
         Zeff = 1.0d0    ! electron charge
      endif
      
      ! Maximal phonon frequency is defined by the maximal phononic CDF peak:
      hw_ph_max = maxval( Material%CDF_phonon%E0(:) + 10.0d0*Material%CDF_phonon%Gamma(:) )  ! Testing
!       hw_ph_max = maxval( Material%CDF_phonon%E0(:) + Material%CDF_phonon%Gamma(:) )
      
      El_elast = numpar%El_elastic   ! to chose the model of CS calculations below
      ! In case it is delta-CDF model, we have to make sure the energy is not below its applicability limit:
      if (El_elast /= 4) then
         max_E0 = maxval(Material%CDF_phonon%E0(:))
         call find_Wmax_equal_Wmin(g_me, Material%Mean_Mass, .false., Ee, 1.0d-16, max_E0, Eeq, hw_ph_max)   ! module "CS_integration_limits"
         ! Check if delta-functional CDF works here,  and apply for electrons above chosen 100 eV:
         !if (Ee < max(numpar%CDF_Eeq_elast*Eeq, 100.0d0)) then   ! switch to nonrelativistic numerically-integrable CDF:
         if (Ee < numpar%CDF_Eeq_elast*Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
             El_elast = 4
         endif
      endif ! (El_elast == 1)
!       print*, 'E_el:', Ee, Eeq, El_elast

      
!       call CDF_total_CS_nonrel(numpar, sigma, Se1, 10000.0d0, g_me, Zeff, 1.0d-20, Material%T_eV, Material%CDF_phonon, Material%Mean_Mass, &
!                         Material%DOS%k, Material%DOS%Eff_m, .false., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
!                         hw_phonon = hw_ph_max)  ! module "CS_electrons_inelastic"
!       print*, 'get_el_elastic_CS', 10000.0d0, sigma
!       pause 'get_el_elastic_CS'
      
      select case (El_elast)  ! chose which model for electron inelastic cross section to use
      case (1,5)  ! Delta relativistic CDF
         if (present (E_max)) then
            sigma = CDF_total_CS_delta(El_elast, Ee, g_me, Zeff, 1.0d-16, Material%At_Dens, Material%CDF_phonon, &
                                                     Material%Mean_Mass, .false., hw_phonon=hw_ph_max, Emax_in=E_max)    ! module "CDF_delta"
         else
            sigma = CDF_total_CS_delta(El_elast, Ee, g_me, Zeff, 1.0d-16, Material%At_Dens, Material%CDF_phonon, &
                                                        Material%Mean_Mass, .false., hw_phonon=hw_ph_max)    ! module "CDF_delta"
         endif
         if (present(Se)) Se = energy_loss_delta(El_elast, Ee, g_me, 1.0d0, Zeff, 1.0d-16, Material%At_Dens, Material%Mean_Mass, &
                                          Material%CDF_phonon, .false., hw_ph_max) ! module "CDF_delta"
      case (4) ! nonrelativistic Ritchie or Mermin CDF
         if (present (E_max)) then
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, 1.0d-20, Material%T_eV, Material%CDF_phonon, Material%Mean_Mass, &
                        Material%DOS%k, Material%DOS%Eff_m, .false., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                        hw_phonon = hw_ph_max, Wmax_in=E_max)  ! module "CS_electrons_inelastic"
         else
            !call CDF_total_CS_nonrel(numpar, sigma, Se1, 1000.0d0, g_me, Zeff, 1.0d-20, Material%T_eV, Material%CDF_phonon, Material%Mean_Mass, & ! Testing
            call CDF_total_CS_nonrel(numpar, sigma, Se1, Ee, g_me, Zeff, 1.0d-20, Material%T_eV, Material%CDF_phonon, Material%Mean_Mass, &
                        Material%DOS%k, Material%DOS%Eff_m, .false., 1.0d0, Material%At_Dens, Material%DOS, numpar%CDF_model, &
                        hw_phonon = hw_ph_max)  ! module "CS_electrons_inelastic"
         endif
         if (present(Se)) Se = Se1
      end select
!       
!       print*, 'get_el_elastic_CS', Ee, Eeq, El_elast, sigma
!       pause 'get_el_elastic_CS'

   case (2) ! Mott
      if (present(mu_max_in)) then
         sigma = Mott_total_CS(Ee, Element%Zat, mu_max_in=mu_max_in)	! below
      else
         sigma = Mott_total_CS(Ee, Element%Zat)	! below
      endif
      
   case (3) ! DSF (Relativistic)
      ! NOT READY YET
   
   case default ! exclude
      sigma = 0.0d0
      if (present(Se)) Se = 0.0d0
   end select
end subroutine get_el_elastic_CS



!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! Mott's cross section according to [3]
pure function Mott_total_CS(Ee, Z, mass, mu_max_in) result(sigma)
   real(8) sigma	! [A^2] cross section
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   integer, intent(in) :: Z	! Ion atomic number
   real(8), intent(in), optional :: mass    ! [me]
   real(8), intent(in), optional :: mu_max_in  ! mu=cos(theta), integration limit
   real(8) :: nu, beta, v, dbleZ, beta2, me, mu_max
   if (present(mass)) then
      me = mass*g_me
   else
      me = g_me
   endif
   if (present(mu_max_in)) then
      if (mu_max_in > 1.0d0) then ! total CS
         mu_max = 1.0d0
      else if (mu_max_in < -1.0d0) then ! zero
         mu_max = -1.0d0
      else  ! partially integrated CS, used for definition of the transfered energy
         mu_max = mu_max_in
      endif
   else ! total CS
      mu_max = 1.0d0
   endif
   v = velosity_from_kinetic_energy(Ee, me, afs=.false.)    ! [m/s] module "Relativity"
   if (v < 1.0d-6) then ! immobile particle does not scatter
      sigma = 0.0d0
   else
      beta = beta_factor(v)	! module "Relativity"
      beta2 = beta*beta
      dbleZ = dble(Z)
      if (present(mass)) then
         nu = screening_parameter(beta, dbleZ, Ee, me)	! below
      else
         nu = screening_parameter(beta, dbleZ, Ee)	! below
      endif
      sigma = g_Pi*g_r0*g_r0*dbleZ*(dbleZ + 1.0d0) * (mu_max + 1.0d0) / ((2.0d0*nu + 1.0d0 - mu_max)*(nu + 1.0d0)) * (1.0d0 - beta2)/(beta2*beta2)	! [3] int of Eq.(39)

!       ! Just for testing Vova's cross section:
!       sigma = Lipp_Mott_test(Ee, 1.0d0, dbleZ)  ! below
!
!       print*, Ee, nu
!       print*, '--------------------------'

   endif
end function Mott_total_CS


 function Lipp_Mott_test(E_ev, meff, z_real) result(Elastic_cross_section)
   real(8) Elastic_cross_section    ! [A^2]
   real(8), intent(in) :: E_ev, meff, z_real
   real(8) :: Pi, clight, alpha, me, el, re, mass, beta, tau, eta_prime, beta_sq, eta_prime_test
	Pi = 3.1415926535897932d0
	clight = 299792458.d0		! speed of light (m/s)
	alpha = 0.0072973525664d0	! fine-structure constant, ~1/137
	me = 9.10938356d-31		! electron mass
	el = 1.6021766208d-19		! electron charge
	re = 2.8179403267d-15		! classical electron radius

	mass = meff*me

	beta =      2.0d0*E_ev*el/(mass*clight**2)	! speed of the electron in units of c, or simply v/c,  squared
	tau = E_ev*el/(mass*clight**2)			! kinetic energy of the electron in units of mc^2

	eta_prime_test = 0.25d0*(alpha/(0.885d0))**2 / beta * (z_real)**(2.d0/3.d0) * ( 1.13d0 + 3.76d0*(alpha*z_real)**2/beta**2 * dsqrt(tau/(tau+1.d0)) )	! Epmirically modified Moliere's screening

	!Elastic_cross_section = Pi*z_real*(z_real+1.d0)*re**2 *(1.d0-beta**2)/ ( eta_prime*(eta_prime+1.d0) * beta**4 ) * 1.d20	! m^2 -> A^2

	beta_sq = 2.0d0*E_ev*el/(mass*clight**2)	! speed of the electron in units of c, or simply v/c,  squared
	eta_prime = 0.25d0*(alpha/(0.885d0))**2 / beta_sq * (z_real)**(2.d0/3.d0) * ( 1.13d0 + 3.76d0*(alpha*z_real)**2/beta_sq * dsqrt(tau/(tau+1.d0)) )	! Epmirically modified Moliere's screening

	Elastic_cross_section = Pi*z_real*(z_real+1.d0)*re**2 *(1.d0-beta_sq)/ ( eta_prime*(eta_prime+1.d0) * beta_sq**2 )  * 1.d20	! m^2 -> A^2

    print*, E_ev, eta_prime, eta_prime_test

end function Lipp_Mott_test



! Solution of diff.Mott's cross section:
pure function Mott_sample_mu(Ee, Z, RN, mass) result(theta)
   real(8) theta	! deflection angle
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   integer, intent(in) :: Z	! Ion atomic number
   real(8), intent(in) :: RN    ! random number [0,1]
   real(8), intent(in), optional :: mass    ! in units of [me]
   real(8) :: nu, beta, v, dbleZ, beta2, me, mu
   if (present(mass)) then
      me = mass*g_me
   else
      me = g_me
   endif
   v = velosity_from_kinetic_energy(Ee, me, afs=.false.)    ! [m/s] module "Relativity"
   if (v < 1.0d-6) then ! immobile particle does not scatter
      theta = 0.0d0
   else
      beta = beta_factor(v)	! module "Relativity"
      beta2 = beta*beta
      dbleZ = dble(Z)
      if (present(mass)) then
         nu = screening_parameter(beta, dbleZ, Ee, me)	! below
      else
         nu = screening_parameter(beta, dbleZ, Ee)	! below
      endif
      mu = (RN*(2.0d0*nu + 1.0d0) - nu) / (RN + nu)
      theta = ACOS(mu)
   endif
end function Mott_sample_mu




pure function screening_parameter(beta, Z, Ee, me) result(nu)
   real(8) nu	! screening parameter
   real(8), intent(in) :: beta	! relativistic beta
   real(8), intent(in) :: Z	! atomic number
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   real(8), intent(in), optional :: me  ! [kg] mass of the particle
   real(8) :: beta2, tau, Erest
   beta2 = beta*beta
   if (present(me)) then
      Erest = rest_energy(me)	! module "Relativity"
   else
      Erest = rest_energy(g_me)	! module "Relativity"
   endif
   tau = Ee / Erest
   nu = 1.7d-5*Z**m_two_third*(1.0d0 - beta2)/beta2 * ( 1.13d0 + 3.76d0 * g_alpha*g_alpha/beta2 * Z*Z * sqrt(tau/(1.0d0 + tau)) )	! [4] Page 160, Eq. (7.8)
end function screening_parameter


end module CS_electrons_elastic
