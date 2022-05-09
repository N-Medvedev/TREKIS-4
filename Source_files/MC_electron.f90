! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019-2021
! 1111111111111111111111111111111111111111111111111111111111111
! [1] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2014: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2015)
! Module contains Monte Carlo routines to simulate electron events
module MC_electron
use Universal_constants
use Objects
use Geometries ! ,only: m_tollerance_eps
use MC_general_tools, only: flight_time_to_boundary, out_of_simulation_box, particle_cross_box_boundary, check_shell, check_element, &
                                deflect_velosity, get_positron_flight_time, get_hole_flight_time, get_electron_flight_time, &
                                sample_phi, sample_theta, transfered_E_from_theta, cos_theta_from_W, cos_recoil_from_W, check_phi, &
                                get_photon_flight_time, remove_particle_from_MC_array, extend_MC_array, find_the_target, &
                                define_normal_to_surface, electron_transmission_probability, reflection_from_surface, &
                                electron_transmission_probability_step, add_energy_into_MD_array
use CS_electrons_elastic, only: get_el_elastic_CS, Mott_sample_mu
use CS_electrons_inelastic, only: get_inelastic_energy_transfer, get_integral_inelastic_CS, CDF_total_CS_nonrel
use CS_electrons_Bremsstrahlung, only: sample_Bremsstrahlung_theta, get_energy_transfer_Bremsstrahlung, sample_Bremsstrahlung_electron_angles
use CS_general_tools, only: total_CS_from_chennels, MFP_from_sigma, Time_from_MFP, find_type_of_scattering, find_valence_hole_mass
use CS_integration_limits, only: W_max, find_Wmax_equal_Wmin
use Relativity, only: velosity_from_kinetic_energy, beta_factor, rest_energy, kinetic_energy_from_velosity
use Little_subroutines, only: interpolate_data_single, print_time_step
use Dealing_with_DOS, only: select_energy_dos
use SHI_charge_state, only: Equilibrium_charge_SHI

implicit none

 contains

subroutine MC_electron_event(used_target, numpar, N_e, Prtcl, NOP, MC, MD_supce, E_e_at, E_e, E_h)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(inout) :: N_e   ! number of electrons
   type(Electron), dimension(:), allocatable, intent(inout) :: Prtcl  ! the electron to perform some event
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout) :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e_at, E_e, E_h ! data to pass to MD later
   !-----------------------------------
   integer :: i_type   ! type of event: -1=box crossing; 0=boundary; 1=inelastic; 2=elastic; 3=Bremsstrahlung
   integer :: INFO  ! to check if procidure executed correctly
   real(8) :: t0, R0(3)    ! save for checking purposes
   
   if ( (isnan(Prtcl(NOP)%R(1))) .or. (isnan(Prtcl(NOP)%R(2))) .or. (isnan(Prtcl(NOP)%R(3))) .or. &
       (SUM(ABS(Prtcl(NOP)%V(:))) < 1.0d-12) .or. (Prtcl(NOP)%Ekin < 0.0d0) ) then
      print*, 'Error in MC_electron_event:',  Prtcl(NOP)%active
      print*, 'N=', NOP, MC%N_e
      print*, Prtcl(NOP)%R(:)
      print*, Prtcl(NOP)%R0(:) 
      print*, Prtcl(NOP)%V(:)
      print*, Prtcl(NOP)%V0(:)
      print*, Prtcl(NOP)%Ekin
      print*, Prtcl(NOP)%ti - Prtcl(NOP)%t0, Prtcl(NOP)%ti, Prtcl(NOP)%t0
   endif
   
   R0 = Prtcl(NOP)%R0   ! save for testing
   ! 0) Advance the last-step parameters:
   Prtcl(NOP)%R0(:) = Prtcl(NOP)%R(:)
   Prtcl(NOP)%V0(:) = Prtcl(NOP)%V(:)
   
   ! 1) Move the particle to the end-point:
   Prtcl(NOP)%R(:) = Prtcl(NOP)%R0(:) + (Prtcl(NOP)%V(:)) * (Prtcl(NOP)%ti - Prtcl(NOP)%t0) !* g_ms2Afs  
   t0 = Prtcl(NOP)%t0   ! save for testing
   Prtcl(NOP)%t0 = Prtcl(NOP)%ti
      
   ! 2) Choose an electron event:
   call find_type_of_electron_event(used_target, numpar, Prtcl(NOP), i_type)  ! below

   ! 3) Perform the event according to the chosen type:
   select case(i_type)
   case (0) ! target boundary crossing
      call event_electron_target_boundary(used_target, numpar, Prtcl(NOP), NOP, MD_supce, E_e, INFO)   ! below
      if (INFO /= 0) then   ! somthing went wrong in the boundary scattering:
         print*, 'Error in MC_electron_event'
         print*, 'problem with event_electron_target_boundary'
         print*, 't0=', t0
         print*, 'R0=', R0
         print*, '****************************************'
      endif
   case (1) ! inelastic
      call event_electron_inelastic(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)   ! below
   case (2) ! elastic
      call event_electron_elastic(used_target, numpar, MC, NOP, MD_supce, E_e, E_e_at)   ! below
   case (3) ! Bremsstrahlung
      call event_electron_Bremsstrahlung(used_target, numpar, MC, NOP, MD_supce, E_e)   ! below
   case default ! simulation box boundary crossing
      call particle_cross_box_boundary(used_target, numpar, N_e, NOP, Prtcl)    ! module "MC_general_tools"
   endselect
end subroutine MC_electron_event


!ССССССССССССССССССССССССССССССССССССССССССС
! Electron reflecting from or crossing target boundary:
subroutine event_electron_target_boundary(used_target, numpar, Prtcl, NOP, MD_supce, E_e, INFO)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Electron), intent(inout) :: Prtcl        ! electron as an object
   integer, intent(in) :: NOP   ! number of electron
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) ::E_e  ! data to pass to MD later
   integer, intent(inout) :: INFO  ! info about errors to pass to main subroutine
   !------------------------------------------------------
   real(8), dimension(3)  :: R_shift, norm_to_surf, V_perp, V
   real(8) :: Vnorm, Ekin_norm, T, RN, Vabs, Z
   integer :: target_ind 
   logical :: transmit
   type(Emission_barrier), pointer :: Em_Barr
   
   ! Find whether an electron crosses the boundary or reflects back:
   ! 1) Kinetic energy towards crossing the boundary:
   call define_normal_to_surface(used_target,  Prtcl, norm_to_surf, 'event_electron_target_boundary', INFO)    ! module "MC_general_tools"
   if (INFO /=0) then   !some error occured
      print*, 'Electron: define_normal_to_surface failed for Prtcl:', NOP
      print*, 'R=', Prtcl%R
      print*, 'V=', Prtcl%V
      print*, 'E=', Prtcl%Ekin
      print*, 'ti=', Prtcl%ti
      print*, 't_sc=', Prtcl%t_sc
      print*, 'Inside target #', Prtcl%in_target
      print*, 'n to surface:', norm_to_surf
   endif
   
   ! 2) Find the velosity projection onto the normal to the surface:
   Vnorm = scalar_projection(Prtcl%V, norm_to_surf)     ! module "Geometries"
   Vnorm = Vnorm * g_Afs2ms ! [fs/A] -> [m/s]
   
   ! 3) Get the normal component of the kinetic energy:
   Ekin_norm = kinetic_energy_from_velosity(Vnorm, g_me)    ! module "Relativity"
   
   ! 4) Get the probability of boundary crossing:
   Em_Barr => used_target%Material(Prtcl%in_target)%Surface_barrier
   ! Projection of energy on the normal to use to calculate transmission probability:
!    T = electron_transmission_probability(Em_Barr%Work_func, Em_Barr%gamma, Em_Barr%E1, Ekin_norm)     ! module "MC_general_tools"
   ! Total kinetic energy to use to calculate transmission probability:
   select case(Em_Barr%barr_type)
   case (1) ! Eckart-type barrier
      T = electron_transmission_probability(Em_Barr%Work_func, Em_Barr%gamma, Em_Barr%E1, Prtcl%Ekin)     ! module "MC_general_tools"
   case default ! step barrier
      T = electron_transmission_probability_step(Em_Barr%Work_func, Ekin_norm) ! module "MC_general_tools"
   end select
   
   ! 5) Sample transmission vs. reflection:
   transmit = .true.   ! to start with, assume transmission
   RN = 0.0d0
   if (T < 0.99999999d0) then   ! reflection is plausible:
      call random_number(RN)
      if (RN > T) transmit = .false.    ! it will be reflection
   endif

!     print*, 'Barr', Em_Barr%gamma, Em_Barr%E1, Ekin_norm
!     print*, 'T=', T, RN, Ekin_norm, Prtcl%Ekin
!      print*, 'V0=', Prtcl%V(:)
!      print*, 'R0=', Prtcl%R(:)

   ! 6) Model reflection or transmission:
   if (transmit) then   ! transmission into another target:
      ! Find into which target this electron enters:
      ! Find which target's boundary the particle is crossing:
      !Vabs = SQRT( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
      !R_shift = m_tollerance_eps * Prtcl%V(:)/Vabs     ! to place particle inside of the material
      R_shift = 10.0d0*m_tollerance_eps * Prtcl%V(:)/abs(Prtcl%V(:))     ! to place particle inside of the material
      ! Update particle's material index according to the new material it enters:
      call find_the_target(used_target, Prtcl, R_shift) ! module "MC_general_tools"
!        print*, 'Transmission'
         
      ! 6.a) Shift electron just across the border (along the direction of velocity):
      !Vabs = SQRT( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
      !Prtcl%R(:) = Prtcl%R(:) + 1.0d-6 * Prtcl%V(:)/Vabs
      Prtcl%R(:) = Prtcl%R(:) + R_shift(:)
      
      if (INFO /=0) then   !some error occured
         print*, 'After Transmission:'
         print*, 'R=', Prtcl%R(:)
         print*, 'V=', Prtcl%V(:)
         print*, 'Vnorm=', Vnorm
         print*, 'Ekin_norm=', Ekin_norm
         print*, 'Shift:', R_shift(:)
      endif

      ! 6.b) Get the next flight inside the new target (transmitted) or old target (reflected):
      call get_electron_flight_time(used_target, numpar, Prtcl, MD_supce, E_e)  ! module "MC_general_tools"
   else ! reflection back
      ! Change velosity according to reflection from the surface with given normal:
      call reflection_from_surface(Prtcl%V, norm_to_surf)   ! module "MC_general_tools"
      if (INFO /=0) then   !some error occured
         print*, 'After Reflection:'
         print*, 'R=', Prtcl%R(:)
         print*, 'V=', Prtcl%V
         print*, 'Vnorm=', Vnorm
         print*, 'Ekin_norm=', Ekin_norm
      endif
      
      ! 6.c) Shift electron just across the border (along the direction of velocity):
!       Vabs = SQRT( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
!       Prtcl%R(:) = Prtcl%R(:) + m_tollerance_eps * Prtcl%V(:)/Vabs
      Prtcl%R(:) = Prtcl%R(:) + 10.0d0*m_tollerance_eps * Prtcl%V(:)/abs(Prtcl%V(:))     ! to place particle inside of the material
      
      ! To check exceptional cases of particles hitting a corner of a target and emitting through another wall,
      ! update particle's material index according to the new material it may enter:
      call find_the_target(used_target, Prtcl) ! module "MC_general_tools"

      ! 6.d) Get the next flight inside the new target (transmitted) or old target (reflected):
      call get_electron_flight_time(used_target, numpar, Prtcl, MD_supce, E_e, .true.)  ! module "MC_general_tools"
   endif
   
!      print*, 'V1=', Prtcl%V(:)
!      print*, 'R1=', Prtcl%R(:)
!     print*, 'Flight time:', Prtcl%ti, Prtcl%t_sc
! !    pause
   
   nullify(Em_Barr)
end subroutine event_electron_target_boundary



!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
! Electron Bremsstrahlung emission subroutines:
subroutine event_electron_Bremsstrahlung(used_target, numpar, MC, NOP, MD_supce, E_e)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e  ! data to pass to MD later
   !------------------------------------------------------
   type(Electron), pointer :: Prtcl  ! the electron to perform some event
   type(Target_atoms), pointer :: matter
   type(Atom_kind), pointer :: Element
   real(8) :: RN, phi, theta, CS_tot, CS_sampled, CS_sum, CS_elem, N_elem, elem_contrib, CS_cur, Vtot0, Vtot, a0(3), a(3), dE, V_ph(3)
   real(8) :: phi_ph, theta_ph, beta, A_coef, B_coef
   integer :: j, KOA, KOA1
   logical :: found_shl
   
   ! 0) Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Electrons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! 1) Select the element the electron scatters on:
   ! Get the total cross section:
   call interpolate_data_single(matter%El_Brems_total%E(:), matter%El_Brems_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
   ! Get the random number, to sample partial cross section:
   call random_number(RN)   ! intrinsic FORTRAN subroutine
   ! Realized CS:
   CS_sampled = RN*CS_tot
   ! Find which of the partial cross sections correspond to the sampled one:
   CS_sum = 0.0d0   ! to start with
   N_elem = dble(SUM(matter%Elements(:)%percentage))   ! number of elements in this compound material
   ! Scan through all the atoms' shells to find the one:
   found_shl = .false.  ! to start with
   KOA = matter%N_Elements  ! by default, it is the last element
   ! (no need to check, because true by exclusion if all the other elements are not true):
   if (KOA > 1) then    ! find element
      SRCH:do j =1, matter%N_Elements-1     ! for each element, expect for the last one
         Element => matter%Elements(j)          ! all information about this element
         elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         ! For each shell, find its partial cross section:
         call interpolate_data_single(Element%El_brems%E,  Element%El_brems%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
         CS_sum = CS_sum + CS_cur*elem_contrib
         CS_elem = CS_cur  ! save for later use
         ! Check if this shell is the sampled one according to the integral CS:
         call check_element(CS_sum, CS_sampled, j, KOA1, found_shl)    ! module "MC_general_tools"
         if (found_shl) then
            KOA = KOA1
            exit SRCH  ! no need to continue the search if we already found the KOA
         endif
      enddo SRCH
   else ! only one element, no neet to search
      Element => matter%Elements(KOA)          ! all information about this element
      ! For each shell, find its partial cross section:
      call interpolate_data_single(Element%El_brems%E,  Element%El_brems%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
      CS_elem = CS_cur  ! save for later use
   endif
   
   if ( (KOA > size(matter%Elements(:))) .or. (KOA < 1) ) then
      print*, "Error 1 in event_electron_Bremsstrahlung:", CS_sampled, RN*CS_tot, CS_tot
   endif
   
   if (isnan(CS_elem) .or. (CS_elem <= 0.0d0)) then
      print*, "Error 2 in event_electron_Bremsstrahlung:", Prtcl%Ekin, CS_elem
   endif
   
   ! Now we know which element the electron scatters on:
   Element => matter%Elements(KOA)  ! all information about this element
   
   ! 2) Get the energy of emitted photon:
   dE = get_energy_transfer_Bremsstrahlung(Prtcl%Ekin, Element, numpar, KOA, CS_elem)    ! module "CS_electrons_Bremsstrahlung"
   
   ! 3) Get the cosines of the electron velosity:
   Vtot = sqrt( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
   a0(:) = Prtcl%V(:)/Vtot
   Vtot0 = Vtot ! save it for later setting photon parameters
   
   ! 4) Get the electron deflection angles:
   call sample_Bremsstrahlung_electron_angles(phi,theta)    ! module "CS_electrons_Bremsstrahlung"
   
   ! 5a) Get the deflection angle (angles of scattering of electron in the laboratory system):
!    call  deflect_velosity(a0(1), a0(2), a0(3), phi, theta, a(1), a(2), a(3))  ! module "MC_general_tools"
   call  deflect_velosity(a0(1), a0(2), a0(3), theta, phi, a(1), a(2), a(3))  ! module "MC_general_tools"
   ! 5b) Change electron energy accordingly:
   Prtcl%Ekin = Prtcl%Ekin - dE ! [eV]
   ! 5c) Change the velosity accordingly:
   Vtot = velosity_from_kinetic_energy(Prtcl%Ekin, g_me)  ! [A/fs] module "Relativity"
   Prtcl%V(:) = Vtot * a(:)  ! components of the electron velosity
   
   ! 6) Now update the electron parameters:
   call get_electron_flight_time(used_target, numpar, Prtcl, MD_supce, E_e)  ! module "MC_general_tools"
   
   ! 7) Create emitted photon:
   MC%N_ph = MC%N_ph + 1
   ! 7a) Get the photon emission angles:
   phi_ph = sample_phi()   ! module "MC_general_tools"
   beta = beta_factor(Vtot0, afs_in = .true.)    ! module "Relativity"
   ! coefficients we don't know yet [1], so set the default ones:
   A_coef = 1.0d0
   B_coef = 0.0d0
   theta_ph = sample_Bremsstrahlung_theta(beta, A_coef, B_coef)   ! module "CS_electrons_Bremsstrahlung"
   ! 7b) Get the cosines of the emitting electron velosity:
!    call deflect_velosity(a0(1), a0(2), a0(3), phi_ph, theta_ph, a(1), a(2), a(3))  ! module "MC_general_tools"
   call deflect_velosity(a0(1), a0(2), a0(3), theta_ph, phi_ph, a(1), a(2), a(3))  ! module "MC_general_tools"
!    V_ph(:) =g_cvel * a(:)   ! derection of photon motion [m/s]
   V_ph(:) =g_c_Afs * a(:)   ! derection of photon motion [A/fs]
   ! 7c) Save the photon data into the particle properties:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_ph > size(MC%MC_Photons)) call extend_MC_array(MC%MC_Photons)   ! module "MC_general_tools"
   call make_new_particle(MC%MC_Photons(MC%N_ph), Ekin=dE, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_ph)    ! module "Objects"
   ! 7d) Define the next photon scattering event:
   call get_photon_flight_time(used_target, numpar, MC%MC_Photons(MC%N_ph))  ! module "MC_general_tools"
   
   nullify(Prtcl, matter, Element)
end subroutine event_electron_Bremsstrahlung




!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Electron inelastic scattering subroutines:
subroutine event_electron_inelastic(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h  ! data to pass to MD later
   !------------------------------------------------------
   integer :: KOA, NSH, i
   real(8) :: Ekin, Ekin_new, dE, theta, phi, theta_r, phi_r, eps, Ekin_SAVE
   real(8) :: V0(3), V(3), V_2(3), V_tot, V_e_abs, V_e(3), E_DOS, a0(3), a(3)
   real(8) :: polar_angle_h, azimuth_angle_h, Mh, mass_temp, V_h_abs, V_h(3)
   real(8) :: CS_total
   type(Electron), pointer :: Prtcl  ! the electron to perform some event
   type(Target_atoms), pointer :: matter
!    type(Atom_kind), pointer :: Element
   logical :: valent
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Electrons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
!   do i = 1, 1000   ! Testing
!   Prtcl%Ekin = 1000.0d0    ! Testing
!   Prtcl%Ekin = matter%DOS%Egap + dble(i) ! Testing
   
   if ( (Prtcl%Ekin < matter%DOS%Egap) .and. (allocated(matter%CDF_valence%A)) ) then   ! it cannot perform ionization, just let it fly:
!       print*, 'Too small Ekin e:', Prtcl%Ekin, matter%DOS%Egap
      ! Get the old electron next sampled free path:
      call get_electron_flight_time(used_target, numpar, Prtcl, MD_supce, E_e)  ! module "MC_general_tools"
      goto 9999 ! skip all the scattering stuff
   endif
   
   ! 1) Which shell of which element is being ionized:
   call select_shell_ionization(used_target, numpar, MC, NOP, KOA, NSH)  ! below
   
   ! 2) Sample transfered energy in the impact ionization event:
   ! Sample the transfered energy:
   if (KOA == 0) then   ! valence band:
      valent = .true.
      ! Find the precalculated total cross section:
      call interpolate_data_single(matter%El_inelastic_valent%E, matter%El_inelastic_valent%Total, &
                        Prtcl%Ekin, CS_total) ! module "Little_subroutines"
      ! Pass it into energy transfer finding:
      dE = get_inelastic_energy_transfer(Prtcl%Ekin, matter, matter%DOS, numpar, KOA, NSH, matter%DOS%Egap, CS_total)    ! module "CS_electrons_inelastic"
   else ! core shell:
      valent = .false.
      ! Find the precalculated total cross section:
      call interpolate_data_single(matter%Elements(KOA)%El_inelastic%E, matter%Elements(KOA)%El_inelastic%Per_shell(NSH,:), &
                        Prtcl%Ekin, CS_total) ! module "Little_subroutines"
      ! Pass it into energy transfer finding:
      dE = get_inelastic_energy_transfer(Prtcl%Ekin, matter, matter%DOS, numpar, KOA, NSH, &
                        matter%Elements(KOA)%Ip(NSH), CS_total)    ! module "CS_electrons_inelastic"
   endif
!    dE = matter%DOS%Egap + dble(i) ! Testing
   
   ! Just in case of numerical inaccuracies, make sure transfered energy is not above the electron energy:
   if (dE > Prtcl%Ekin) then
      print*, 'Error in event_electron_inelastic 0', dE, Prtcl%Ekin, valent
      dE = Prtcl%Ekin
!    else if (dE < 0.0d0) then
   else if (dE < matter%DOS%Egap) then
      print*, 'Error in event_electron_inelastic 1', dE, Prtcl%Ekin, valent
      dE = matter%DOS%Egap
   endif
   
   ! 2a) Scattering angles according to the transferred energy:
   call sample_angles_inelastic(Prtcl%Ekin, dE, theta, phi) ! below
   
   ! And the angles of emission of a new electron:
   call set_angles_of_new_electron(Prtcl%Ekin, dE, phi, theta_r, phi_r)  ! below
   
   ! 3) Energy an electron is emitted with:
   if (KOA == 0) then   ! valence band:
      ! Sample from where within the valence band it is ionized according to DOS:
      call select_energy_DOS(matter%DOS%E, matter%DOS%DOS, matter%Integral_DOS_fe, &
                                matter%DOS%Egap, dE, matter%DOS%alpha_CB, E_DOS)   ! module "Dealing_with_DOS"
      Ekin = dE + (E_DOS - matter%DOS%Egap)
   else ! core shell:
      Ekin = dE - matter%Elements(KOA)%Ip(NSH)
      E_DOS = 0.0d0
   endif
   ! Make sure the transfered energy is sufficient:
   if (Ekin < 0.0d0) then   ! insufficient energy is replaced with the minimal allowed:
      print*, "Error in event_electron_inelastic 2:", Ekin, E_DOS, Prtcl%Ekin, dE
!       pause 
      Ekin = 0.0d0
      if (KOA == 0) then   ! valence band:
         dE = -E_DOS + matter%DOS%Egap
      else
         dE = matter%Elements(KOA)%Ip(NSH)
      endif
   endif
   
   ! 4) Change the incident electron direction of motion:
   ! Get the cosines of the electron velosity:
   V_tot = sqrt(SUM(Prtcl%V(:)*Prtcl%V(:)))
   V0(:) = Prtcl%V(:)/V_tot
   
   
   ! Get the deflection angle (angles of emission of photoelectron in the laboratory system):
   call deflect_velosity(V0(1), V0(2), V0(3), theta, phi, V(1), V(2), V(3))  ! module "MC_general_tools"
   ! Get the NEW absolute electron velosity, to set its components according to the cosines:
   Ekin_SAVE = Prtcl%Ekin   ! Used to check up below
   
   Prtcl%Ekin = Prtcl%Ekin - dE ! update old electron energy
   
   if (Prtcl%Ekin < 0.0d0) then
      print*, 'Error in event_electron_inelastic (Ekin<0)', Prtcl%Ekin, dE
      print*, E_DOS, KOA, NSH
      print*, Ekin, matter%DOS%Egap
      ! Fix the energy to the original value:
      Prtcl%Ekin = Prtcl%Ekin + dE ! update old electron energy
      print*, 'Energy patched to the original value:', Prtcl%Ekin
      ! Get the old electron next sampled free path:
      call get_electron_flight_time(used_target, numpar, Prtcl, MD_supce, E_e)  ! module "MC_general_tools"
      goto 9999 ! skip all the scattering stuff if it's impossible
   endif
   V_e_abs = velosity_from_kinetic_energy(Prtcl%Ekin, g_me)     !  [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the electron velosity

   ! Update new velosity:
   Prtcl%V(:) = V_e(:)
   
   ! Get the old electron next sampled free path:
   call get_electron_flight_time(used_target, numpar, Prtcl, MD_supce, E_e)  ! module "MC_general_tools"
   
   ! 4a) Get the direction of the new electron:
   call deflect_velosity(V0(1), V0(2), V0(3), theta_r, phi_r, V(1), V(2), V(3))  ! module "MC_general_tools"
   V_e_abs = velosity_from_kinetic_energy(Ekin, g_me)     !  [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the electron velosity

!    write(*,'(a,f,f,f,e,e,e)') 'e:', dE, Ekin, theta_r, phi, V_e(:)
!    enddo ! Testing
!    PAUSE 'event_electron_inelastic'
   
   ! 5) Make a new electron:
   MC%N_e = MC%N_e + 1
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_e > size(MC%MC_Electrons)) then
      nullify(Prtcl, matter)    ! nullify pointers
      call extend_MC_array(MC%MC_Electrons)   ! module "MC_general_tools"
      ! Pointers for easier access to scattering particle and target properties:
      Prtcl => MC%MC_Electrons(NOP)
      matter => used_target%Material(Prtcl%in_target)
   endif
   call make_new_particle(MC%MC_Electrons(MC%N_e), Ekin=Ekin, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_e)    ! module "Objects"
   
   ! 6) Get new electrons time of the next event:
   call get_electron_flight_time(used_target, numpar, MC%MC_Electrons(MC%N_e), MD_supce, E_e)  ! module "MC_general_tools"
   
   if (MC%MC_Electrons(MC%N_e)%ti < Prtcl%t0) then
      print*, '(event_electron_inelastic)', MC%N_e, MC%MC_Electrons(MC%N_e)%ti, Prtcl%ti, Prtcl%t0
   endif
   
   
   ! 7) Create a new hole after electron removal:
   MC%N_h = MC%N_h + 1
   ! for valent hole, sample its velosity:
   if (valent) then
      ! Find hole's mass:
      Mh =  find_valence_hole_mass(numpar, used_target%Material(Prtcl%in_target)%DOS, abs(E_DOS)) ! module "MC_general_tools"
      ! Set the velosity of the new hole:
      V_h_abs = velosity_from_kinetic_energy(abs(E_DOS), Mh)    !  [A/fs] module "Relativity"
      polar_angle_h = sample_phi() ! module "MC_general_tools"
      azimuth_angle_h = sample_theta()  ! module "MC_general_tools"
      V0(1) = 1.0d0
      V0(2:3) = 0.0d0
      ! Get the holes velosity in the sampled direction:
!       call  deflect_velosity(V0(1), V0(2), V0(3), polar_angle_h, azimuth_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      call  deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle_h, polar_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      V_h(:) = V_h_abs * V_h(:)   ! set also the absolute value
   else
      E_DOS = 0.0d0
      V_h(:) = 0.0d0   ! core holes don't fly
   endif
   ! Save parameters of the new particles into array:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_h > size(MC%MC_Holes)) call extend_MC_array(MC%MC_Holes)    ! module "MC_general_tools"
   call make_new_particle(MC%MC_Holes(MC%N_h), Ekin=-(E_DOS), t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_h, &
                                      KOA = KOA, Sh = NSH, valent = valent )    ! module "Objects"
   
   if (.not.valent .and. (KOA == 0)) then ! inconsistency found
      print*, 'event_electron_inelastic: valency of hole error'
      print*, valent, KOA, NSH, MC%N_h
   endif
   
!    print*, 'Ekin after:', Prtcl%Ekin, dE, Prtcl%Ekin+dE
   if (valent) then
      if (abs(Ekin_SAVE - (Prtcl%Ekin + MC%MC_Holes(MC%N_h)%Ekin + matter%DOS%Egap + MC%MC_Electrons(MC%N_e)%Ekin)) > 1.0d-12*Ekin_SAVE) then
         print*, 'Error in event_electron_inelastic 3: conservation'
         print*, 'E_in', Ekin_SAVE, abs(Ekin_SAVE - (Prtcl%Ekin + MC%MC_Holes(MC%N_h)%Ekin + matter%DOS%Egap + MC%MC_Electrons(MC%N_e)%Ekin)) 
         print*, 'Eoth after:', MC%MC_Holes(MC%N_h)%Ekin, matter%DOS%Egap, MC%MC_Electrons(MC%N_e)%Ekin
         print*, 'Etot :', Prtcl%Ekin + MC%MC_Holes(MC%N_h)%Ekin + matter%DOS%Egap + MC%MC_Electrons(MC%N_e)%Ekin
         print*, 'dE:',  Ekin, dE , E_DOS
      endif
   else
!       write(*,'(a,f,f,f)') 'event_electron_inelastic:', MC%MC_Holes(MC%N_h)%R
      if (abs(Ekin_SAVE - (Prtcl%Ekin + MC%MC_Holes(MC%N_h)%Ekin + matter%Elements(KOA)%Ip(NSH) + MC%MC_Electrons(MC%N_e)%Ekin)) &
                > 1.0d-12*Ekin_SAVE) then
         print*, 'Error in event_electron_inelastic 4: conservation'
         print*, 'E_in', Ekin_SAVE, abs(Ekin_SAVE - (Prtcl%Ekin + MC%MC_Holes(MC%N_h)%Ekin + matter%Elements(KOA)%Ip(NSH) + MC%MC_Electrons(MC%N_e)%Ekin)) 
         print*, 'Eoth after:', MC%MC_Holes(MC%N_h)%Ekin, matter%Elements(KOA)%Ip(NSH), MC%MC_Electrons(MC%N_e)%Ekin
         print*, 'Etot :', Prtcl%Ekin + MC%MC_Holes(MC%N_h)%Ekin + matter%Elements(KOA)%Ip(NSH) + MC%MC_Electrons(MC%N_e)%Ekin
         print*, 'dE:',  Ekin, dE , E_DOS
      endif
   endif
   

   ! 8) Get new holes time of the next event:
   call get_hole_flight_time(used_target, numpar, MC%MC_Holes(MC%N_h), MD_supce, E_h)  ! module "MC_general_tools"
   
9999 nullify(Prtcl, matter)
end subroutine event_electron_inelastic


subroutine select_shell_ionization(used_target, numpar, MC, NOP, KOA, NSH)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   integer, intent(out) :: KOA, NSH ! kind of atom, and the shell number, which absorbs the photon
   !------------------------------------------------------
   type(Target_atoms), pointer :: matter
   type(Electron), pointer :: Prtcl  ! the photon to perform some event
   type(Atom_kind), pointer :: Element
   real(8) :: RN, CS_tot, CS_sampled, CS_cur, CS_sum,  elem_contrib, N_elem
   integer :: i_mat, j, k
   logical :: found_shl, valence_done
   
   ! Starting default values (valence band):
   KOA = 0
   NSH = 0
   valence_done = .false.   ! to start, we did not calculate valence band yet
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Electrons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! Get the total cross section:
!    Prtcl%Ekin = 3000.0d0    ! Testing
   call interpolate_data_single(matter%El_inelastic_total%E(:), matter%El_inelastic_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
   ! Get the random number, to sample partial cross section:
   call random_number(RN)   ! intrinsic FORTRAN subroutine
   ! Realized CS:
!    RN = 0.999d0 ! Testing
   CS_sampled = RN*CS_tot

   ! Find which of the partial cross sections correspond to the sampled one:
   CS_sum = 0.0d0   ! to start with
   N_elem = dble(SUM(matter%Elements(:)%percentage))   ! number of elements in this compound material
   ! Scan through all the atoms' shells to find the one:
   SRCH:do j =1, matter%N_Elements	! for each element
      Element => matter%Elements(j)	! all information about this element
      elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
!       KOA = j
      ! Check the atomic CSs:
      do k = 1, Element%N_shl    ! for all shells in this elements
!          NSH = k
         
         VAL:if ((numpar%El_inelast /= 2) .and. (Element%valent(k)) .and. (allocated(matter%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
            if (.not.valence_done) then ! add VB (but only once):
               call interpolate_data_single(matter%El_inelastic_valent%E,  matter%El_inelastic_valent%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
               CS_sum = CS_sum + CS_cur  ! contribution from the band
               valence_done = .true.     ! we added VB once (not do double-count it)
            endif
         else VAL
            call interpolate_data_single(Element%El_inelastic%E,  Element%El_inelastic%Per_shell(k,:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"      
            CS_sum = CS_sum + CS_cur*elem_contrib   ! contribution from the element
         endif VAL
         
!          CS_sum = CS_sum + CS_cur*elem_contrib
         
         if (.not.Element%valent(k)) then ! core orbitals
            call check_shell(CS_sum, CS_sampled, j, k, KOA, NSH, found_shl)    ! module "MC_general_tools"
         else  ! valence orbitals
            if (allocated(matter%CDF_valence%A)) then    ! valence band
               call check_shell(CS_sum, CS_sampled, 0, 0, KOA, NSH, found_shl)    ! module "MC_general_tools"
            else    ! valence shells
               call check_shell(CS_sum, CS_sampled, j, k, KOA, NSH, found_shl)    ! module "MC_general_tools"
            endif
         endif ! (Element%valent(k))
         
!          write(*,'(a,i2,i2,e,e,e,e,L2)') 'select_shell_ionization', j, k, CS_sum, CS_cur, CS_sampled, CS_tot, found_shl ! Testing
         
         if (found_shl) exit SRCH  ! no need to continue the search if we already found the KOA nad NSH
      enddo ! k
   enddo SRCH !  j =1, N_elements

   if ( KOA>matter%N_Elements ) then
      print*, 'Error in select_shell_ionization', KOA, NSH, Prtcl%Ekin, CS_sampled, CS_tot, RN
   endif
   
   ! Clean up the memory:
   nullify(matter, Prtcl, Element)
end subroutine select_shell_ionization


subroutine sample_angles_inelastic(Ekin, dE, theta, phi)
   real(8), intent(in) :: Ekin, dE  ! [eV] incident and transfered energy
   real(8), intent(out) :: theta, phi   ! scattering angles
   real(8) :: mu
   ! Assume uniform scattering probability into 2Pi for phi:
   phi = sample_phi()    ! module "MC_general_tools"
   ! Assume absolutely elastic scattering for theta (no recoil):
   mu = cos_theta_from_W(Ekin, dE, g_me, g_me)    ! module "MC_general_tools"
   theta = acos(mu)
end subroutine sample_angles_inelastic


pure subroutine set_angles_of_new_electron(Ekin, dE, phi, theta_r, phi_r)  
   real(8), intent(in) :: Ekin, dE  ! [eV] incident and transfered energy
   real(8), intent(in) :: phi   ! angle of deflection of the original particle
   real(8), intent(out) :: theta_r, phi_r   ! emission angles of a new electron
   real(8) :: mu
   ! Assume phis are in the same plane:
   phi_r = phi + g_Pi
   call check_phi(phi_r) ! module "MC_general_tools"
   ! Assume absolutely elastic scattering for theta:
   mu = cos_recoil_from_W(Ekin, dE, g_me, g_me)    ! module "MC_general_tools"
!    if ((mu < -1.0d0) .or. (mu > 1.0d0)) print*, mu, Ekin, dE
   theta_r = acos(mu)
end subroutine set_angles_of_new_electron





!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! Electron elastic scattering subroutines:

subroutine event_electron_elastic(used_target, numpar, MC, NOP, MD_supce, E_e, E_e_at)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_e_at   ! data to pass to MD later
   !------------------------------------------------------
   type(Electron), pointer :: Prtcl  ! the electron to perform some event
   type(Target_atoms), pointer :: matter
   real(8) :: RN, phi, theta, CS_tot, CS_sampled, CS_sum, N_elem, elem_contrib, CS_cur, Vtot, a0(3), a(3), dE, Ekin
   integer :: j, KOA, KOA1
   logical :: found_at
   
   ! 0) Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Electrons(NOP)
   matter => used_target%Material(Prtcl%in_target)

   ! 1) Set the deflection angles:
   ! Sample azimuthal angle randomly within 2Pi:
   phi = sample_phi()   ! module "MC_general_tools"
   
   ! 2,3) Find out transfered energy and deflecting angle:
   call elastic_energy_transfer(numpar, Prtcl, matter, KOA, dE, theta)   ! below
   ! Just in case, check that the transfered energy is not too high:
   if (dE > Prtcl%Ekin) dE = Prtcl%Ekin ! [eV]
   
   ! 4) Get the cosines of the electron velosity:
   Vtot = sqrt( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
   a0(:) = Prtcl%V(:)/Vtot
   
   ! 5) Get the deflection angle (angles of scattering of electron in the laboratory system):
   call deflect_velosity(a0(1), a0(2), a0(3), theta, phi, a(1), a(2), a(3))  ! module "MC_general_tools"

!    print*, 'event_electron_elastic1', Prtcl%Ekin, dE, Vtot, Prtcl%V(:)
!    print*, 'a0(:)', a0(:), sqrt(SUM(a0(:)*a0(:)))
   Ekin = Prtcl%Ekin    ! save for testing
   Prtcl%Ekin = Prtcl%Ekin - dE ! [eV]
   
   ! 5c) Change the velosity accordingly:
   Vtot = velosity_from_kinetic_energy(Prtcl%Ekin, g_me)  !  [A/fs] module "Relativity"
   ! Update velosity:
   Prtcl%V(:) = Vtot * a(:)  ! components of the electron velosity
   
   ! 6) Now update the electron parameters:
   call get_electron_flight_time(used_target, numpar, Prtcl, MD_supce, E_e)  ! module "MC_general_tools"
   
   ! 7) Save the atomic in this collision parameters:
   MC%N_at_nrg = MC%N_at_nrg + 1
   if (MC%N_at_nrg > size(MC%MC_Atoms_events)) call extend_MC_array(MC%MC_Atoms_events)    ! module "MC_general_tools"
   ! Save the parameters of this collision (ONLY ENERGY TRANSFER IS SAVED, MOMENTUM NOT DONE YET!)
   call make_new_particle(MC%MC_Atoms_events(MC%N_at_nrg), Ekin=dE, t0=Prtcl%t0, KOA = KOA, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R)    ! module "Objects"
   
   ! Save the energy to pass to MD module, if needed:
   if (numpar%DO_MD) then   ! if user requested MD at all
      call add_energy_into_MD_array(numpar, dE, Prtcl%R, MD_supce, E_e_at)   ! module "MC_general_tools"
!       print*, 'event_electron_elastic 2', dE, MD_supce%coord_dim
   endif

   if (abs(Ekin - (Prtcl%Ekin+MC%MC_Atoms_events(MC%N_at_nrg)%Ekin)) > 1.0d-10) then
      print*, 'Error in event_electron_elastic', Ekin
      print*, Prtcl%Ekin, MC%MC_Atoms_events(MC%N_at_nrg)%Ekin, Prtcl%Ekin+MC%MC_Atoms_events(MC%N_at_nrg)%Ekin
      print*, 'N', MC%N_at_nrg, MC%MC_Atoms_events(MC%N_at_nrg)%active, COUNT(MC%MC_Atoms_events(:)%active), &
              SUM(MC%MC_Atoms_events(:)%Ekin, MASK = MC%MC_Atoms_events(:)%active)
   endif
!     print*, 'event_electron_elastic2', Prtcl%Ekin, Vtot, Prtcl%ti, Prtcl%t_sc, Prtcl%V(:)
!     print*, 'a(:)', a(:), sqrt(SUM(a(:)*a(:)))
   
   nullify(Prtcl, matter)
end subroutine event_electron_elastic



subroutine elastic_energy_transfer(numpar, Prtcl, matter, KOA, dE, theta)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Electron), intent(in) :: Prtcl  ! the electron to perform some event
   type(Target_atoms), intent(in), target :: matter ! material parameters
   integer, intent(out) :: KOA    ! kind of atom
   real(8), intent(out) :: dE, theta    ! [eV] transfered energy, and scattering angle
   !--------------------------------
   real(8) :: CS_tot, RN, CS_sampled, CS_sum, N_elem, elem_contrib, CS_cur, Eeq, max_E0, Zeff
   real(8) :: E_left, E_right, E_cur, mu, hw_phonon, mtc2, Mc2, Se1
   real(8) :: eps, eps1, eps2, eps3, CS_tot_test
   logical :: found_at
   integer :: j, KOA1, El_inelast, El_elast
   type(Atom_kind), pointer :: Element
   
   ! Set accepteble margin of precision for the angle:
   eps1 = 1.0d-2
   eps2 = 1.0d-4
   eps3 = 1.0d-6
   
   select case (numpar%El_elastic)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF, 5=SP-CDF
   case (1,5)  ! CDF
!       eps = 1.0d-3
      E_left = 0.0d0! [eV] minimal transferred energy
!       E_right = matter%E_debye    ! [eV] maximal transferred energy
      hw_phonon = maxval( matter%CDF_phonon%E0(:) + matter%CDF_phonon%Gamma(:) )
      mtc2 = rest_energy(matter%Mean_Mass)   ! target rest energy; module "Relativity"
      Mc2 = rest_energy(g_me)   ! incident electron rest energy; module "Relativity"
      !E_right = W_max(g_me, mtc2, .false., Prtcl%Ekin, 1.0d-8, hw_phonon) ! module "CS_integration_limits"
      E_right = W_max(Mc2, mtc2, .false., Prtcl%Ekin, 1.0d-8, hw_phonon) ! module "CS_integration_limits"
      
      El_elast = numpar%El_elastic
      if ((El_elast == 1) .or. (El_elast == 5)) then
         max_E0 = maxval(matter%CDF_phonon%E0(:))
         call find_Wmax_equal_Wmin(g_me, matter%Mean_Mass, .false., Prtcl%Ekin, 1.0d-6, max_E0, Eeq, hw_phonon)   ! module "CS_integration_limits"
         ! Check if delta-functional CDF works here,  and apply for electrons above chosen threshold:
         if (Prtcl%Ekin < numpar%CDF_Eeq_elast*Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
             El_elast = 4
         endif
      endif ! (El_elast == 1)
      
      ! Sample the cross section:
      call random_number(RN)
      Element => matter%Elements(1) ! unused here, so just to pass something into the subroutine
!       call get_el_elastic_CS(Prtcl%Ekin, matter, Element, numpar, CS_tot_test)   ! module "CS_electrons_elastic"
      ! Use precalculated total cross section:
      call interpolate_data_single(matter%El_elastic_total%E(:), matter%El_elastic_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
!       print*, 'Elast:', Prtcl%Ekin, CS_tot_test, CS_tot   ! Testing
      
      ! Sample the cross section:
      CS_sampled = RN*CS_tot
   
      if ( (RN < 1.0d-10) .or. (CS_sampled < 1.0d-12) ) then   ! no need to sample, it's just the lower limit
         E_cur = E_left
      else ! sample it
         if (El_elast /= 4) then    ! for analytical CSs, use bisection search algorithm:
            ! sampled cross-sections closer to the total one require higher precision:
            if (RN > 0.99d0) then    
               eps = eps3
            elseif (RN > 0.8d0) then
               eps = eps2
            else
               eps = eps1
            endif
      
            ! Start finding CS:
            E_cur = (E_left + E_right)*0.5d0
            call get_el_elastic_CS(Prtcl%Ekin, matter, Element, numpar, CS_cur, E_max=E_cur)   ! module "CS_electrons_elastic"
            ! Search by bisection method:
!             do while (ABS(CS_cur - CS_sampled) > eps*CS_sampled)
            do while (abs(E_left - E_right) >= eps2*E_left)
                if (CS_cur > CS_sampled) then
                    E_right = E_cur
                else
                    E_left = E_cur
                endif
                E_cur = (E_left + E_right)/2.0d0
                if (abs(E_left - E_right) < eps2) exit  ! precise enough
                call get_el_elastic_CS(Prtcl%Ekin, matter, Element, numpar, CS_cur, E_max=E_cur)   ! module "CS_electrons_elastic"
            enddo
!             print*, 'elastic_energy_transfer 1', CS_sampled, CS_cur, E_cur
         else   ! for numerically integrable CS, use direct search
             ! Target mean atomic number:
             if (numpar%CDF_elast_Zeff /= 1) then
                Zeff = 1.0d0 + Equilibrium_charge_SHI(Prtcl%Ekin, g_me, matter%Mean_Z, (matter%Mean_Z-1.0d0), 0, 1.0d0) ! module "SHI_charge_state"
             else
                Zeff = 1.0d0    ! electron charge
             endif
             
             call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Prtcl%Ekin, g_me, Zeff, 1.0d-10, matter%T_eV, matter%CDF_phonon, matter%Mean_Mass, &
                        matter%DOS%k, matter%DOS%Eff_m, .false., 1.0d0, matter%At_Dens, matter%DOS, numpar%CDF_model, &
                        hw_phonon = hw_phonon, Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
!             print*, 'elastic_energy_transfer 2', CS_sampled, CS_cur, E_cur
         endif  ! (El_elast /= 4)
      endif ! ( (RN < 1.0d-10) .or. (CS_sampled < 1.0d-12) )
   
      ! Output: sampled transferred energy:
      dE = E_cur
   
      ! Get the corresponding scattering angle:
      mu = cos_theta_from_W(Prtcl%Ekin, dE, g_me, matter%Mean_Mass) ! module "MC_general_tools"
!       print*, 'mu', mu, Prtcl%Ekin, dE, E_left, CS_sampled, CS_tot
      theta = acos(mu)

   case (2)  ! Mott
   
      ! 1) Select the element the electron scatters on:
      ! Get the total cross section:
      call interpolate_data_single(matter%El_elastic_total%E(:), matter%El_elastic_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
      ! Get the random number, to sample partial cross section:
      call random_number(RN)   ! intrinsic FORTRAN subroutine
      ! Realized CS:
      CS_sampled = RN*CS_tot
      ! Find which of the partial cross sections correspond to the sampled one:
      CS_sum = 0.0d0   ! to start with
      N_elem = dble(SUM(matter%Elements(:)%percentage))   ! number of elements in this compound material
      ! Scan through all the atoms to find the one:
      found_at = .false.   ! to start with
      KOA = matter%N_Elements  ! by default, it is the last element
      ! (no need to check, because true by exclusion if all the other elements are not true):
      SRCH:do j =1, matter%N_Elements-1     ! for each element, expect for the last one
         Element => matter%Elements(j)          ! all information about this element
         elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         ! For each shell, find its partial cross section:
         call interpolate_data_single(Element%El_elastic%E,  Element%El_elastic%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
         CS_sum = CS_sum + CS_cur*elem_contrib
         ! Check if this shell is the sampled one according to the integral CS:
         call check_element(CS_sum, CS_sampled, j, KOA1, found_at)    ! module "MC_general_tools"
!          if (found_at) exit SRCH  ! no need to continue the search if we already found the KOA
         if (found_at) then
            KOA = KOA1
            exit SRCH  ! no need to continue the search if we already found the KOA
         endif
      enddo SRCH
   
      if ( (KOA > size(matter%Elements(:))) .or. (KOA < 1) ) then
         print*, 'Elast', CS_sampled, RN*CS_tot, CS_tot
      endif
   
      ! Now we know which element the electron scatters on:
      Element => matter%Elements(KOA)  ! all information about this element
   
      ! 2) Sample polar angle according to the differential cross section
      call  get_electron_elastic_polar_angle(Prtcl%Ekin, Element, theta) ! below
         
      ! 3) Change electron energy accordingly:
      dE = transfered_E_from_theta(Prtcl%Ekin, theta, g_me, Element%M) ! module "MC_general_tools"

   case (3) ! DSF
      ! NOT READY YET
   case default ! no scattering
      dE = 0.0d0
      theta = 0.0d0
   end select

   nullify(Element)
end subroutine elastic_energy_transfer



subroutine get_electron_elastic_polar_angle(Ekin, Element, theta)
   real(8), intent(in) :: Ekin  ! [eV] photon energy
   type(Atom_kind), intent(in) :: Element   ! parameters of the element the photon scatters on
   real(8), intent(out) :: theta   ! polar angle of photoelectron emission
   real(8) :: RN
   ! Sample the angle:
   call random_number(RN)
   ! The scattering angle is:
   theta = Mott_sample_mu(Ekin, Element%Zat, RN)  ! module "CS_electrons_elastic"
end subroutine get_electron_elastic_polar_angle



!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

subroutine find_type_of_electron_event(used_target, numpar, Prtcl, i_type)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Electron), intent(in) :: Prtcl  ! the electron to perform some event
   integer, intent(out) :: i_type   ! type of event: -1=box crossing; 0=boundary; 1=inelastic; 2=elastic; 3=Bremsstrahlung
   !----------------------------------------------
   real(8) :: eps, MFP, V
   logical :: out
   type(Target_atoms), pointer :: matter
   
   !eps = 1.0d-12
   eps = 1.0d-10

!    print*, 'Elbe', Prtcl%Ekin, Prtcl%ti, Prtcl%t_sc, Prtcl%R(:), Prtcl%V(:)
   
   ! 1) Check if it is a boundary crossing:
   if (abs(Prtcl%ti - Prtcl%t_sc) > eps) then   ! it is a boundary crossing
      ! 2) Check whether it is simulation box crossing or a target crossing:
      out = out_of_simulation_box(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), numpar%box_start_x, numpar%box_end_x, &
               numpar%box_start_y, numpar%box_end_y, numpar%box_start_z, numpar%box_end_z) ! module "MC_general_tools"
      if (out) then ! it is crossing the simulation box boundary
         i_type = -1
      else  ! it is crossing the target boundary
         i_type = 0
      endif
      
!        print*, 'El3', Prtcl%Ekin, i_type, Prtcl%R(:)
!        print*, 'Box:', numpar%box_start_x, numpar%box_end_x, &
!                 numpar%box_start_y, numpar%box_end_y, numpar%box_start_z, numpar%box_end_z
   else ! it is a scattering-type event
      ! Properties of the target material, inside of which the particle is:
      matter => used_target%Material(Prtcl%in_target)
      
!       print*, 'El2', Prtcl%Ekin, Prtcl%ti, Prtcl%t_sc, Prtcl%R(:), Prtcl%V(:)
      
      ! 3) Find out what kind of scattering event it is:
      call find_type_of_scattering(i_type, Prtcl%Ekin,  matter%El_inelastic_total%E, matter%El_inelastic_total%Total, &
                    E_array2=matter%El_elastic_total%E, CS_array2=matter%El_elastic_total%Total, &
                    E_array3=matter%El_Brems_total%E, CS_array3=matter%El_Brems_total%Total) ! module "CS_general_tools"
   endif
   
!    if (Prtcl%Ekin < 8.9) then
!       if (i_type == 1) print*,'(impossible electron scattering)', Prtcl%Ekin, i_type
!    endif
!    print*, 'El1', i_type
   
!    print*, 'Elaf', Prtcl%Ekin, i_type, Prtcl%ti, Prtcl%t_sc
   
   nullify(matter)
end subroutine find_type_of_electron_event


   
   
end module MC_electron
