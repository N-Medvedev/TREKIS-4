! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019-2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains Monte Carlo routines to simulate photon events
! References used in the module:
! [1]  F. Salvat, J. M. Fernandez-Varea, E. Acosta, J. Sempau "PENELOPE -2015 A Code System for Monte Carlo Simulation of Electron and Photon Transport", OECD (2015)

module MC_photon
use Universal_constants
use Objects
use Geometries ! ,only: m_tollerance_eps
use MC_general_tools, only: flight_time_to_boundary, out_of_simulation_box, particle_cross_box_boundary, check_shell, check_element, &
                                            deflect_velosity, get_electron_flight_time, get_positron_flight_time, get_hole_flight_time, get_photon_flight_time, &
                                            sample_phi, sample_theta, remove_particle_from_MC_array, extend_MC_array, &
                                            put_back_into_box, find_the_target
use CS_general_tools, only: total_CS_from_chennels, MFP_from_sigma, Time_from_MFP, find_type_of_scattering, find_valence_hole_mass
use CS_photons_Rayleigh, only: get_Rayleigh_CS
use CS_photons_Compton, only: Compton_CS, Compton_Ec, Compton_mu_from_Ec
use CS_photons_pair_creation, only: Pair_CS
use Little_subroutines, only: sample_Poisson, interpolate_data_single
use Relativity, only: velosity_from_kinetic_energy, rest_energy, beta_factor, gamma_factor
use Dealing_with_DOS, only: select_energy_DOS

implicit none

 contains

subroutine MC_photon_event(used_target, numpar, N_ph, Prtcl, NOP, MC, MD_supce, E_e, E_h)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(inout) :: N_ph   ! number of photons
   type(Photon), dimension(:), allocatable, intent(inout) :: Prtcl  ! the photon to perform some event
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout) :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h ! data to pass to MD later
   !---------------------------------------
   integer :: i_type    ! type of event: -1=crossing box boundary; 0=crossing target boundary; 1=absorption; 2=elastic; 3=Compton; 4=pair creation
   
   ! 0) Advance the last-step parameters:
   Prtcl(NOP)%R0(:) = Prtcl(NOP)%R(:)
   Prtcl(NOP)%V0(:) = Prtcl(NOP)%V(:)
   
   ! 1) Move the particle to the end-point:
   Prtcl(NOP)%R(:) = Prtcl(NOP)%R0(:) + Prtcl(NOP)%V(:) * (Prtcl(NOP)%ti - Prtcl(NOP)%t0) !* g_ms2Afs
   Prtcl(NOP)%t0 = Prtcl(NOP)%ti
   
   ! 1.5) For photons (neutral particles), we use simplified periodic boundaries:
   ! place the photon back into the simulation box, just in case it escaped it:
   call put_back_into_box(numpar, Prtcl(NOP))    ! module "MC_general_tools"

   ! 2) Choose photon event:
   call find_type_of_photon_event(used_target, numpar, Prtcl(NOP), i_type)  ! below
   
!    if (Prtcl(NOP)%R(3) < 0.0d0) then
!       print*, 'ph:', NOP, Prtcl(NOP)%R(3), Prtcl(NOP)%Ekin, i_type, Prtcl(NOP)%in_target
!    endif
   
   ! 3) Knowing the type of event, act accordingly:
   select case(i_type)
   case (0) ! target boundary crossing
!       print*, 'event_photon_target_boundary'
      call event_photon_target_boundary(used_target, numpar, Prtcl(NOP))   ! below
      
   case (1) ! absorption
      ! Create new electron by photon:
      call event_photoabsorption(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)  ! below
      ! Photon is absorbed and disappears from the MC array:
      call remove_particle_from_MC_array(N_ph, NOP, Prtcl)  ! module "MC_general_tools"
      
   case (2) ! coherent
!       print*, 'MC_photon_event photon elastic', N_ph, i_type
      call event_photon_elastic(used_target, numpar, MC, NOP)   ! below
      
   case (3) ! Compton
!        print*, 'MC_photon_event photon Compton', N_ph, i_type
      call event_photon_Compton(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)  ! below
      
   case (4) ! e-e+ pair creation
      ! Pair production by photon:
      call event_photon_pair_production(used_target, numpar, MC, NOP, MD_supce, E_e)   ! below
      ! Photon is absorbed and disappears from the MC array:
      call remove_particle_from_MC_array(N_ph, NOP, Prtcl)   ! module "MC_general_tools"
      
   case default ! simulation box boundary crossing
!       print*, 'particle_cross_box_boundary', Prtcl(NOP)%Ekin, sqrt(SUM(Prtcl(NOP)%V(:)*Prtcl(NOP)%V(:))), Prtcl(NOP)%R(:), N_ph
      call particle_cross_box_boundary(used_target, numpar, N_ph, NOP, Prtcl)    ! module "MC_general_tools"
!       print*, 'MC_photon_event photon removed', N_ph, i_type

   endselect
end subroutine MC_photon_event


!ССССССССССССССССССССССССССССССССССССССССССС
! Photon crossing target boundary:
subroutine event_photon_target_boundary(used_target, numpar, Prtcl)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Photon), intent(inout) :: Prtcl        ! electron as an object
   !------------------------------------------------------
   real(8), dimension(3)  :: R_shift
   real(8) :: PrtclV
   integer :: target_ind
   
   ! Find into which target this photon enters:
   ! 1) Find which target's boundary the particle is crossing:
   !R_shift = 1.0d-7 * Prtcl%V(:)     ! to place particle inside of the material
   PrtclV = max( SQRT( SUM( Prtcl%V(:)*Prtcl%V(:) ) ), m_tollerance_eps) ! exclude zero
   R_shift = m_tollerance_eps * Prtcl%V(:)/PrtclV     ! to place particle inside of the material
   ! Update particle's material index according to the new material it enters:
   call find_the_target(used_target, Prtcl, R_shift) ! module "MC_general_tools"
   
!    print*, 'ph target_ind', Prtcl%in_target
   
   ! Test case (assume vacuum):
!    Prtcl%in_target = 0   ! vacuum
   call get_photon_flight_time(used_target, numpar, Prtcl)  ! module "MC_general_tools"

end subroutine event_photon_target_boundary


!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
! Photon e-e+ pair-production subroutines:
subroutine event_photon_pair_production(used_target, numpar, MC, NOP, MD_supce, E_e)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e  ! data to pass to MD later
   !------------------------------------------------------
   type(Photon), pointer :: Prtcl  ! the photon to perform some event
   type(Target_atoms), pointer :: matter
   type(Atom_kind), pointer :: Element
   real(8) :: Ee, Epos, theta_e, phi_e, theta_p, phi_p, N_elem, elem_contrib, CS_cur, CS_sum, CS_sampled
   real(8) :: V0(3), V(3), V_e(3), V_e_abs
   integer :: j, KOA, KOA1
   logical :: found_shl
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Photons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! 1) Choose the element the photon scatters on:
   N_elem = dble(SUM(matter%Elements(:)%percentage))   ! number of elements in this compound material
   ! Scan through all the atoms' shells to find the one:
   found_shl = .false.  ! to start with
   KOA = matter%N_Elements  ! by default, it is the last element
   ! (no need to check, because true by exclusion if all the other elements are not true):
   SRCH:do j =1, matter%N_Elements-1     ! for each element, expect for the last one
      Element => matter%Elements(j)	! all information about this element
      elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
      ! For each shell, find its partial cross section:
      call interpolate_data_single(Element%Phot_pair%E,  Element%Phot_pair%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
      CS_sum = CS_sum + CS_cur*elem_contrib
      ! Check if this shell is the sampled one according to the integral CS:
      call check_element(CS_sum, CS_sampled, j, KOA1, found_shl)    ! module "MC_general_tools"
!       if (found_shl) exit SRCH  ! no need to continue the search if we already found the KOA
      if (found_shl) then
         KOA = KOA1
         exit SRCH  ! no need to continue the search if we already found the KOA
      endif
   enddo SRCH
   ! Now we know which element the photon scatters on:
   Element => matter%Elements(KOA)  ! all information about this element
   
   ! 2) Get energy of the emitted electron:
   call emitted_pair_electron_energy(Element, Prtcl%Ekin, Ee) ! below
   
   ! 3) Make a new electron:
   MC%N_e = MC%N_e + 1
   ! 3a) Sample emitted electron angles:
   call get_pair_electron_angles(Ee, theta_e, phi_e) ! below   
   
   ! 3b) Set the direction of motion:
!    V0(:) = MC%MC_Photons(NOP)%V(:)/g_cvel   ! [m/s]
   V0(:) = Prtcl%V(:)/ sqrt(SUM(Prtcl%V(:)*Prtcl%V(:)))   ! [A/fs]
   ! Get the deflection angle (angles of emission of photoelectron in the laboratory system):
   call  deflect_velosity(V0(1), V0(2), V0(3), theta_e, phi_e, V(1), V(2), V(3))  ! module "MC_general_tools"
   ! Get the absolute electron velosity, to set its components according to the cosines:
   V_e_abs = velosity_from_kinetic_energy(Ee, g_me)     ! [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the electron velosity

   ! 3c) Set the parameters of the electron into MC array:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_e > size(MC%MC_Electrons)) then 
      call extend_MC_array(MC%MC_Electrons)   ! module "MC_general_tools"
   endif
   call make_new_particle(MC%MC_Electrons(MC%N_e), Ekin=Ee, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_e)    ! module "Objects"
   
   ! 4) Get new electrons time of the next event:
   call get_electron_flight_time(used_target, numpar, MC%MC_Electrons(MC%N_e), MD_supce, E_e)  ! module "MC_general_tools"
   
   ! 5) Make a new positron:
   MC%N_p = MC%N_p + 1
   ! Its energy:
   Epos = Prtcl%Ekin - Ee - 2.0d0*g_me_eV
   ! 5a) Sample emitted positron angles:
   call get_pair_electron_angles(Epos, theta_p, phi_p) ! below
   
   ! 5b) Get the deflection angle (angles of emission of photoelectron in the laboratory system):
   call  deflect_velosity(V0(1), V0(2), V0(3), theta_p, phi_p, V(1), V(2), V(3))  ! module "MC_general_tools"
   ! Get the absolute positron velosity, to set its components according to the cosines:
   V_e_abs = velosity_from_kinetic_energy(Epos, g_me)     ! [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the positron velosity
   
   ! 5c) Set the parameters of the positron into MC array:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_p > size(MC%MC_Positrons)) call extend_MC_array(MC%MC_Positrons)   ! module "MC_general_tools"
   call make_new_particle(MC%MC_Positrons(MC%N_p), Ekin=Epos, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_e)    ! module "Objects"
   
   ! 6) Get new positrons time of the next event:
   call get_positron_flight_time(used_target, numpar, MC%MC_Positrons(MC%N_p))  ! module "MC_general_tools"
   
   nullify(Prtcl, matter, Element)
end subroutine event_photon_pair_production


subroutine emitted_pair_electron_energy(Element, Ekin, Ee)
   type(Atom_kind), intent(in) :: Element   ! element properties
   real(8), intent(in) :: Ekin  ! [eV] incomming photon energy
   real(8), intent(out) :: Ee   ! [eV] emitted electron energy
   real(8) :: RN, CS_tot, CS_sampled, eps
   real(8) :: CS_start, CS_cur, E_left, E_right, E_cur, k
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-3
   
   ! Sample the angle:
   call random_number(RN)
   ! Get the starting integration point of the cross section:
   k = Ekin/g_me_eV
   E_left = 1.0d0/k
   CS_start = Pair_CS(Ekin, Element, E_left) ! module "CS_photons_pair_creation"
   E_right = 1.0d0 - 1.0d0/k
   CS_tot = Pair_CS(Ekin, Element, E_right) - CS_start ! module "CS_photons_pair_creation"
   CS_sampled = RN * CS_tot ! sampled cross section kernel
   ! Find the angle that corresponds to the sampled cross section kernel:
   E_cur = (E_left + E_right)/2.0d0
   CS_cur = Pair_CS(Ekin, Element, E_cur) - CS_start ! module "CS_photons_pair_creation"
   ! Search by bisection method:
   do while (ABS(CS_cur - CS_sampled)/CS_sampled > eps)
      if (CS_cur > CS_sampled) then
         E_right = E_cur
      else
         E_left = E_cur
      endif
      E_cur = (E_left + E_right)/2.0d0
      if (abs(E_left - E_right) < eps) exit  ! precise enough
      CS_cur =  Pair_CS(Ekin, Element, E_cur) - CS_start ! module "CS_photons_pair_creation"
   enddo
   
   Ee = E_cur   ! [eV] electron energy sampled
end subroutine emitted_pair_electron_energy


subroutine get_pair_electron_angles(Ekin, theta, phi)
   real(8), intent(in) :: Ekin  ! [eV] emitted electron energy
   real(8), intent(out) :: phi, theta   ! electron scattering angles, azymiuthal and polar
   
   ! Sample azimuthal angle randomly within 2Pi:
   phi = sample_phi()   ! module "MC_general_tools"
   
   ! Sample polar angle according to the differential cross section, Eq(2.49) [1]
   call get_pair_polar_angle(Ekin, theta)  ! below
end subroutine get_pair_electron_angles


subroutine get_pair_polar_angle(Ekin, theta)
   real(8), intent(in) :: Ekin  ! [eV] emitted electron or positron energy
   real(8), intent(out) :: theta    ! polar angle
   real(8) :: RN, beta, v, ksi
   v = velosity_from_kinetic_energy(Ekin, g_me, afs =.false.) ! [m/s] module "Relativity"
   beta = beta_factor(v)    ! modul "Relativity"
   call random_number(RN)
   ksi = 2.0d0*RN - 1.0d0
   theta = (ksi + beta)/(ksi*beta + 1.0d0)  ! Eq.(2.99) in 2015-edition of [1]
end subroutine get_pair_polar_angle


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Photon Compton scattering subroutines:
subroutine event_photon_Compton(used_target, numpar, MC, NOP,  MD_supce, E_e, E_h)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h ! data to pass to MD later
   !------------------------------------------------------
    integer :: KOA, NSH
   real(8) :: Ekin, Ekin_new, dE, polar_angle, azimuth_angle   ! electron energy and angles of emittion
   real(8) :: polar_angle_ph, azimuth_angle_ph     ! photon deflection angle
   real(8) :: polar_angle_h, azimuth_angle_h        ! hole deflection angle
   real(8) :: V0(3), V(3), V_tot, V_e_abs, V_e(3), E_DOS, a0(3), a(3), V_h(3), V_h_abs, Mh, mass_temp, eps, Epot
   type(Photon), pointer :: Prtcl  ! the photon to perform some event
   logical :: valent
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Photons(NOP)
   !matter => used_target%Material(Prtcl%in_target)
   
   ! 1) Compton by which shell of which element:
   call select_shell_Compton(used_target, MC, NOP, KOA, NSH)  ! below
   
   ! 2) Sample emitted electron angles:
   ! We don't have a model for valence band Compton, use atomic shells instead:
   if (KOA == 0) then ! assume uppermost valence shell of the first element
      call get_Compton_photon_angles(used_target%Material(Prtcl%in_target)%Elements(1), &   ! element #1
            size(used_target%Material(Prtcl%in_target)%Elements(1)%Ip), & ! last shell
            Prtcl%Ekin, azimuth_angle, polar_angle)  ! below
   else ! atomic shell as defined:
      call get_Compton_photon_angles(used_target%Material(Prtcl%in_target)%Elements(KOA), NSH, Prtcl%Ekin, azimuth_angle, polar_angle)  ! below
   endif

   ! 2a) Energy of the photon after Compton scattering:
   Ekin_new = Compton_Ec(Prtcl%Ekin, cos(azimuth_angle))   ! module "CS_photons_Compton"
   
   ! 3) Get emitted electron energy:
   ! Transferred energy from the photon:
   dE = Prtcl%Ekin - Ekin_new   ! [eV] energy transfered to the electron
   ! additional check of consistency:
   if (KOA > 0) then    ! for atomic shell scattering
      if (dE < used_target%Material(Prtcl%in_target)%Elements(KOA)%Ip(NSH)) then
         !print*, "Problem in event_photon_Compton: dE < Ip"
         !print*, dE, used_target%Material(Prtcl%in_target)%Elements(KOA)%Ip(NSH), azimuth_angle, cos(azimuth_angle)
         ! Redefine it to conserve energy:
         dE = used_target%Material(Prtcl%in_target)%Elements(KOA)%Ip(NSH)
         Ekin_new = Prtcl%Ekin - dE
      endif
   else ! for VB scattering
      if (dE < used_target%Material(Prtcl%in_target)%DOS%Egap) then
         !print*, "Problem in event_photon_Compton: dE < Egap"
         !print*, dE, used_target%Material(Prtcl%in_target)%DOS%Egap, azimuth_angle, cos(azimuth_angle)
         ! Redefine it to conserve energy:
         dE = used_target%Material(Prtcl%in_target)%DOS%Egap
         Ekin_new = Prtcl%Ekin - dE
      endif
   endif

   if (isnan(dE)) then
      print*, 'Error in event_photon_Compton:'
      print*, 'Transferred energy is dE = ', dE
      print*, 'E=', Prtcl%Ekin, Ekin_new
      print*, 'KOA', KOA, NSH
      print*, 'theta:', azimuth_angle, cos(azimuth_angle)
      print*, 'Ip', used_target%Material(Prtcl%in_target)%Elements(KOA)%Ip(NSH)
      print*, 'Egap', used_target%Material(Prtcl%in_target)%DOS%Egap
   endif

   ! Energy an electron is emitted with:
   if (KOA == 0) then   ! valence band:
   !if (used_target%Material(Prtcl%in_target)%Elements(KOA)%valent(NSH)) then
      valent = .true.
      ! Sample from where within the valence band it is ionized according to DOS:
      call select_energy_DOS(used_target%Material(Prtcl%in_target)%DOS%E, used_target%Material(Prtcl%in_target)%DOS%DOS, &
                                           used_target%Material(Prtcl%in_target)%Integral_DOS_fe, &
                                           used_target%Material(Prtcl%in_target)%DOS%Egap, &
                                           dE, used_target%Material(Prtcl%in_target)%DOS%alpha_CB, E_DOS)   ! module "Dealing_with_DOS"
      Ekin = dE + (E_DOS - used_target%Material(Prtcl%in_target)%DOS%Egap)
      Epot = used_target%Material(Prtcl%in_target)%DOS%Egap ! potential energy of hole, used for testing only!
   else ! core shell:
      valent = .false.
      E_DOS = 0.0d0
      Ekin = dE - used_target%Material(Prtcl%in_target)%Elements(KOA)%Ip(NSH)
      Epot = used_target%Material(Prtcl%in_target)%Elements(KOA)%Ip(NSH)    ! potential energy of hole, used for testing only!
   endif
   ! Make sure the transfered energy is sufficient:
   if (Ekin < 0.0d0) then   ! insufficient energy is replaced with the minimal allowed:
      print*, "Problem in event_photon_Compton: "
      print*, 'Compton Ekin:', Ekin
      Ekin = 0.0d0
      if (KOA == 0) then   ! valence band:
         Ekin_new = -E_DOS + used_target%Material(Prtcl%in_target)%DOS%Egap
      else
         Ekin_new = used_target%Material(Prtcl%in_target)%Elements(KOA)%Ip(NSH)
      endif
      ! Change the transferred energy accordingly to updated value:
      dE = Prtcl%Ekin - Ekin_new   ! [eV] energy transfered to the electron
   endif

   
   ! 4) Get the photoelectron direction of emission:
   ! Get the cosines of the photon velosity:
   V_tot = sqrt(SUM(Prtcl%V(:)*Prtcl%V(:)))
   V0(:) = Prtcl%V(:)/V_tot
   ! Get the deflection angle (angles of emission of photoelectron in the laboratory system):
!    call  deflect_velosity(V0(1), V0(2), V0(3), polar_angle, azimuth_angle, V(1), V(2), V(3))  ! module "MC_general_tools"
   call  deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle, polar_angle, V(1), V(2), V(3))  ! module "MC_general_tools"
   ! Get the absolute electron velosity, to set its components according to the cosines:
   V_e_abs = velosity_from_kinetic_energy(Ekin, g_me)     ! [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the electron velosity

   ! 5) Make a new electron:
   MC%N_e = MC%N_e + 1
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_e > size(MC%MC_Electrons)) call extend_MC_array(MC%MC_Electrons)    ! module "MC_general_tools"
   call make_new_particle(MC%MC_Electrons(MC%N_e), Ekin=Ekin, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_e)    ! module "Objects"
   
   ! 6) Get new electrons time of the next event:
   call get_electron_flight_time(used_target, numpar, MC%MC_Electrons(MC%N_e), MD_supce, E_e)  ! module "MC_general_tools"
   
   ! 7) Create a new hole after electron removal:
    MC%N_h = MC%N_h + 1
    ! for valent hole, sample its velosity:
    if (valent) then
       ! Find hole's mass:
       Mh =  find_valence_hole_mass(numpar, used_target%Material(Prtcl%in_target)%DOS, E_DOS) ! module "MC_general_tools"
       ! Set the velosity of the new hole:
       V_h_abs = velosity_from_kinetic_energy(abs(E_DOS), Mh)    ! [A/fs] module "Relativity"
       polar_angle_h = sample_phi() ! module "MC_general_tools"
       azimuth_angle_h = sample_theta()  ! module "MC_general_tools"
       V0(3) = 1.0d0
       V0(1:2) = 0.0d0
       ! Ger the holes velosity in the sampled direction:   
!        call  deflect_velosity(V0(1), V0(2), V0(3), polar_angle, azimuth_angle, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
       call  deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle, polar_angle, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
       V_h(:) = V_h_abs * V_h(:)   ! set also the absolute value
       KOA = 0
       NSH = 0
    else
       V_h(:) = 0.0d0   ! core holes don't fly
    endif
    ! in case we have more particles than spaces in the array, extend the array:
    if (MC%N_h > size(MC%MC_Holes)) call extend_MC_array(MC%MC_Holes)   ! module "MC_general_tools"
    call make_new_particle(MC%MC_Holes(MC%N_h), Ekin=abs(E_DOS), t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_h, &
                                      KOA = KOA, Sh = NSH, valent = valent )    ! module "Objects"
   
   if ( Prtcl%Ekin - (Ekin_new+Ekin+abs(E_DOS)+Epot) > 1.0d-10) then
      print*, "Problem in event_photon_Compton: "
      print*, 'event_photon_Compton0', Prtcl%Ekin, Ekin_new, dE
      print*, 'event_photon_Compton1', Ekin, E_DOS, Epot
      print*, 'event_photon_Compton2', Ekin_new+Ekin+abs(E_DOS)+Epot
      print*, 'Fin', KOA, NSH, Prtcl%Ekin - (Ekin_new+Ekin+abs(E_DOS)+Epot)
      if (MC%MC_Holes(MC%N_h)%valent) then
         print*, used_target%Material(MC%MC_Holes(MC%N_h)%in_target)%DOS%Egap
      else
         print*, used_target%Material(Prtcl%in_target)%Elements(KOA)%Ip(NSH)
      endif
   endif
   
   
   ! 8) Get new hole's time of the next event:
   call get_hole_flight_time(used_target, numpar, MC%MC_Holes(MC%N_h), MD_supce, E_h)  ! module "MC_general_tools"
   
   ! 9) Update the photon parameters:
   ! 9a) Get the cosines of the photon velosity:
   a0(:) = Prtcl%V(:)/ sqrt(SUM(Prtcl%V(:)*Prtcl%V(:)))  ! using untis of [A/fs]
!    a0(:) = Prtcl%V(:)/g_cvel  ! using untis of [m/s]
   ! photon deflection angle:
   polar_angle_ph = sample_theta()  ! module "MC_general_tools"
   azimuth_angle_ph = azimuth_angle + g_Pi  ! in the same plane as electron
   ! 9b) Get the deflection angle (angles of scattering of photon in the laboratory system):
!    call deflect_velosity(a0(1), a0(2), a0(3), polar_angle_ph, azimuth_angle_ph, a(1), a(2), a(3))  ! module "MC_general_tools"
   call deflect_velosity(a0(1), a0(2), a0(3), azimuth_angle_ph, polar_angle_ph, a(1), a(2), a(3))  ! module "MC_general_tools"
   
   ! 9c) New photon energy:
   Prtcl%Ekin = Ekin_new   ! [eV] photon energy after Compton scattering
   ! Update velosity:
   Prtcl%V(:) = g_c_Afs * a(:)  ! components of the photon velosity
!    Prtcl%V(:) = g_cvel * a(:)  ! components of the photon velosity
   
   ! 9d) Define the next photon scattering event:
   call get_photon_flight_time(used_target, numpar, Prtcl)  ! module "MC_general_tools"
   
   nullify(Prtcl)
end subroutine event_photon_Compton


subroutine select_shell_Compton(used_target, MC, NOP, KOA, NSH)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   integer, intent(out) :: KOA, NSH ! kind of atom, and the shell number, which absorbs the photon
   !------------------------------------------------------
   type(Target_atoms), pointer :: matter
   type(Photon), pointer :: Prtcl  ! the photon to perform some event
   type(Atom_kind), pointer :: Element
   real(8) :: RN, CS_tot, CS_sampled, CS_cur, CS_sum,  elem_contrib, N_elem
   integer :: i_mat, j, k
   logical :: found_shl
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Photons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! Get the total cross section:
   call interpolate_data_single(matter%Ph_Compton_total%E(:), matter%Ph_Compton_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
   ! Get the random number, to sample partial cross section:
   call random_number(RN)   ! intrinsic FORTRAN subroutine
   ! Realized CS:
   CS_sampled = RN*CS_tot

   ! Find which of the partial cross sections correspond to the sampled one:
   CS_sum = 0.0d0   ! to start with
   N_elem = dble(SUM(matter%Elements(:)%percentage))   ! number of elements in this compound material
   ! Scan through all the atoms' shells to find the one:
   SRCH:do j =1, matter%N_Elements	! for each element
      Element => matter%Elements(j)	! all information about this element
      elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
      ! Check the atomic CSs:
      do k = 1, Element%N_shl    ! for all shells in this elements
         call interpolate_data_single(Element%Phot_compton%E,  Element%Phot_compton%Per_shell(k,:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
         CS_sum = CS_sum + CS_cur*elem_contrib
!          if (.not.Element%valent(k)) then ! core orbitals
            call check_shell(CS_sum, CS_sampled, j, k, KOA, NSH, found_shl)    ! module "MC_general_tools"
            if (Element%valent(k) .and. allocated(matter%CDF_valence%A) ) then ! valence band
               KOA = 0
               NSH = 0
            endif
            
!          else  ! valence orbitals
!             if (allocated(matter%CDF_valence%A)) then    ! valence band
!                call check_shell(CS_sum, CS_sampled, 0, 0, KOA, NSH, found_shl)    ! module "MC_general_tools"
!             else    ! valence shells
!                call check_shell(CS_sum, CS_sampled, j, k, KOA, NSH, found_shl)    ! module "MC_general_tools"
!             endif
!          endif ! (Element%valent(k)) 
         if (found_shl) exit SRCH  ! no need to continue the search if we already found the KOA nad NSH
      enddo ! k
   enddo SRCH !  j =1, N_elements

   ! Clean up the memory:
   nullify(matter, Prtcl, Element)
end subroutine select_shell_Compton



subroutine get_Compton_photon_angles(Element, k, Ekin, theta, phi)
   type(Atom_kind), intent(in) :: Element   ! element properties
   integer, intent(in) :: k ! index of the shell
   real(8), intent(in) :: Ekin  ! [eV] photon energy
   real(8), intent(out) :: phi, theta   ! electron scattering angles, azymiuthal and polar
   
   ! Sample azimuthal angle randomly within 2Pi:
   phi = sample_phi()   ! module "MC_general_tools"
   
   ! Sample polar angle according to the differential cross section, Eq(2.48) [1]
   call get_Compton_photon_polar_angle(Element, k, Ekin, theta)  ! below
end subroutine get_Compton_photon_angles



subroutine get_Compton_photon_polar_angle(Element, NSH, Ekin, theta)
   type(Atom_kind), intent(in) :: Element   ! element properties
   integer, intent(in) :: NSH ! index of the shell
   real(8), intent(in) :: Ekin  ! [eV] photon energy
   real(8), intent(out) :: theta   ! polar angle of photoelectron emission
   real(8) :: RN, CS_start, CS_end, CS_tot, CS_sampled, eps
   real(8) :: CS_cur, mu_left, mu_right, mu_cur, mu_Ip, Ip
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-3
   
   ! Sample the angle:
   call random_number(RN)
   ! Get the starting integration point of the cross section:
   Ip = Element%Ip(NSH)   ! ionization potential of this shell
   ! Min possible mu, accounting for Ip:
   mu_Ip = Compton_mu_from_Ec(Ekin, Ip)    ! module "CS_photons_Compton"
!    if (mu_Ip > -1.0d0) print*, 'get_Compton_photon_polar_angle', Ekin, NSH
   mu_left = max(-1.0d0, mu_Ip)  ! to start with
   mu_right = 1.0d0 ! to start with
   CS_start = Compton_CS(Ekin, Element, NSH, mu_left) ! module "CS_photons_Compton"
   CS_end = Compton_CS(Ekin, Element, NSH, mu_right) ! module "CS_photons_Compton"
   CS_tot = CS_end - CS_start   ! total cross section kernel
   CS_sampled = RN * CS_tot ! sampled cross section kernel   
   ! Find the angle that corresponds to the sampled cross section kernel:
   mu_cur = (mu_left + mu_right)/2.0d0
   CS_cur = Compton_CS(Ekin, Element, NSH, mu_cur) - CS_start  ! module "CS_photons_Compton"
   ! Search by bisection method:
   do while (ABS(CS_cur - CS_sampled) > eps * CS_sampled)
      if (CS_cur > CS_sampled) then
         mu_right = mu_cur
      else
         mu_left = mu_cur
      endif
      mu_cur = (mu_left + mu_right)/2.0d0
      if (abs(mu_left - mu_right) < eps) exit  ! precise enough
      CS_cur = Compton_CS(Ekin, Element, NSH, mu_cur) - CS_start  ! module "CS_photons_Compton"
   enddo
   
   ! The scattering angle is:
   theta = ACOS(mu_cur)
end subroutine get_Compton_photon_polar_angle



!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! Photon elastic scattering subroutines:

subroutine event_photon_elastic(used_target, numpar, MC, NOP)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   !------------------------------------------------------
   type(Photon), pointer :: Prtcl  ! the photon to perform some event
   type(Target_atoms), pointer :: matter
   type(Atom_kind), pointer :: Element
   real(8) :: a0(3), a(3), RN, phi, theta
   real(8) :: CS_sampled, CS_sum, CS_tot, CS_cur, N_elem, elem_contrib
   integer :: j, KOA, KOA1
   logical :: found_shl
   
   ! 0) Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Photons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! 1) Set the deflection angles:
   ! Sample azimuthal angle randomly within 2Pi:
   phi = sample_phi()   ! module "MC_general_tools"
   
   ! 2) Select the element the photon scatters on:
   ! Get the total cross section:
   call interpolate_data_single(matter%Ph_coherent_total%E(:), matter%Ph_coherent_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
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
   SRCH:do j =1, matter%N_Elements-1     ! for each element, expect for the last one
      Element => matter%Elements(j)	! all information about this element
      elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
      ! For each shell, find its partial cross section:
      call interpolate_data_single(Element%Phot_coherent%E,  Element%Phot_coherent%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
      CS_sum = CS_sum + CS_cur*elem_contrib
      ! Check if this shell is the sampled one according to the integral CS:
      call check_element(CS_sum, CS_sampled, j, KOA1, found_shl)    ! module "MC_general_tools"
!       if (found_shl) exit SRCH  ! no need to continue the search if we already found the KOA
      if (found_shl) then
         KOA = KOA1
         exit SRCH  ! no need to continue the search if we already found the KOA
      endif
   enddo SRCH
   ! Now we know which element the photon scatters on:
   Element => matter%Elements(KOA)  ! all information about this element
   
   ! 3) Sample polar angle according to the differential cross section, Eq(2.17) [1]
   call get_photon_elastic_polar_angle(CS_cur, Prtcl%Ekin, Element, numpar, theta)  ! below
   
   ! 4) Get the cosines of the photon velosity:
   a0(:) = Prtcl%V(:)/ sqrt(SUM(Prtcl%V(:)*Prtcl%V(:)))  ! using untis of [A/fs]
!    a0(:) = Prtcl%V(:)/g_cvel  ! using untis of [m/s]
   
   ! 5) Get the deflection angle (angles of scattering of photon in the laboratory system):
!    call  deflect_velosity(a0(1), a0(2), a0(3), phi, theta, a(1), a(2), a(3))  ! module "MC_general_tools"
   call  deflect_velosity(a0(1), a0(2), a0(3), theta, phi, a(1), a(2), a(3))  ! module "MC_general_tools"
   ! Update velosity:
   Prtcl%V(:) = g_c_Afs * a(:)  ! components of the photon velosity [A/fs]
!    Prtcl%V(:) = g_cvel * a(:)  ! components of the photon velosity [m/s]
   
   ! 6) Now update the photon parameters:
   call get_photon_flight_time(used_target, numpar, Prtcl)  ! module "MC_general_tools"
   
   nullify(Prtcl, matter, Element)
end subroutine event_photon_elastic



subroutine get_photon_elastic_polar_angle(CS_tot, Eph, Element, numpar, theta)
   real(8), intent(in) :: CS_tot    ! [A^2] total cross section
   real(8), intent(in) :: Eph  ! [eV] photon energy
   type(Atom_kind), intent(in) :: Element   ! parameters of the element the photon scatters on
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(out) :: theta   ! polar angle of photoelectron emission
   real(8) :: RN, CS_start, CS_end, CS_sampled, eps
   real(8) :: CS_cur, mu_left, mu_right, mu_cur
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-3
   
   ! Sample the angle:
   call random_number(RN)
   ! Get the starting integration point of the cross section:
   mu_left = -1.0d0  ! to start with
   mu_right = 1.0d0 ! to start with

   CS_sampled = RN * CS_tot ! sampled cross section
   
   call get_Rayleigh_CS(Eph, Element, numpar, CS_start, mu_left)   ! module "CS_photons_Rayleigh"
   ! Find the angle that corresponds to the sampled cross section:
   mu_cur = (mu_left + mu_right)/2.0d0
   call get_Rayleigh_CS(Eph, Element, numpar, CS_cur, mu_cur)    ! module "CS_photons_Rayleigh"
   CS_cur = CS_cur - CS_start
   ! Search by bisection method:
   do while (ABS(CS_cur - CS_sampled) > eps*CS_sampled)
      if (CS_cur > CS_sampled) then
         mu_right = mu_cur
      else
         mu_left = mu_cur
      endif
      mu_cur = (mu_left + mu_right)/2.0d0
      if (abs(mu_left - mu_right) < eps) exit  ! precise enough
      call get_Rayleigh_CS(Eph, Element, numpar, CS_cur, mu_cur)   ! module "CS_photons_Rayleigh"
      CS_cur = CS_cur - CS_start
   enddo
   
   ! The scattering angle is:
   theta = ACOS(mu_cur)
end subroutine get_photon_elastic_polar_angle



!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Photoabsorption subroutines:

subroutine event_photoabsorption(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h ! data to pass to MD later
   !------------------------------------------------------
   integer :: KOA, NSH
   real(8) :: Ekin, polar_angle, azimuth_angle   ! electron energy and angles of emittion
   real(8) :: V0(3), V(3), V_tot, V_e_abs, V_e(3), V_h(3), E_DOS, Mh, polar_angle_h, azimuth_angle_h, V_h_abs
   type(Photon), pointer :: Prtcl  ! the photon to perform some event
   logical :: valent
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Photons(NOP)
   !matter => used_target%Material(Prtcl%in_target)
   
   ! 1) Photoabsorption by which shell of which element:
   call select_shell_photoabsorption(used_target, MC, NOP, KOA, NSH)  ! below
   
   ! 2) Get emitted electron energy: (hw - Ip)
   if (KOA == 0) then   ! valence band:
      valent = .true.
      ! Sample from where within the valence band it is ionized according to DOS:
      call select_energy_DOS(used_target%Material(Prtcl%in_target)%DOS%E, used_target%Material(Prtcl%in_target)%DOS%DOS, &
                                           used_target%Material(Prtcl%in_target)%Integral_DOS_fe, &
                                           used_target%Material(Prtcl%in_target)%DOS%Egap, &
                                           Prtcl%Ekin, used_target%Material(Prtcl%in_target)%DOS%alpha_CB, E_DOS)   ! module "Dealing_with_DOS"
      Ekin = Prtcl%Ekin + (E_DOS - used_target%Material(Prtcl%in_target)%DOS%Egap)
   else ! core shell:
      valent = .false.
      Ekin = Prtcl%Ekin - used_target%Material(Prtcl%in_target)%Elements(KOA)%Ip(NSH)
   endif
   
   ! 3) Sample electron emission angles:
   !call get_photoelectron_emission_angles(Ekin, polar_angle, azimuth_angle)  ! below
   call get_photoelectron_emission_angles(Ekin, azimuth_angle, polar_angle)  ! below
   
   ! 4) Get the photoelectron direction of emission:
   ! Get the cosines of the photon velosity:
   V_tot = sqrt(SUM(Prtcl%V(:)*Prtcl%V(:)))
   if (V_tot > 0.0d0) then
      V0(:) = Prtcl%V(:)/V_tot
   else
      V0(:) = 0.0d0
   endif
   ! Get the deflection angle (angles of emission of photoelectron in the laboratory system):
!    call  deflect_velosity(V0(1), V0(2), V0(3), polar_angle, azimuth_angle, V(1), V(2), V(3))  ! module "MC_general_tools"
   call  deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle, polar_angle, V(1), V(2), V(3))  ! module "MC_general_tools"
   ! Get the absolute electron velosity, to set its components according to the cosines:
   V_e_abs = velosity_from_kinetic_energy(Ekin, g_me)     ! [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the electron velosity

   ! 5) Make a new electron:
   MC%N_e = MC%N_e + 1
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_e > size(MC%MC_Electrons)) call extend_MC_array(MC%MC_Electrons)   ! module "MC_general_tools"
   call make_new_particle(MC%MC_Electrons(MC%N_e), Ekin=Ekin, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_e)    ! module "Objects"
   
   ! 6) Get new electron's time of the next event:
   call get_electron_flight_time(used_target, numpar, MC%MC_Electrons(MC%N_e), MD_supce, E_e)  ! module "MC_general_tools"
   
!     print*, 'event_photoabsorption 0:', MC%MC_Electrons(MC%N_e)%R(:)
!     print*, 'event_photoabsorption 1:', MC%MC_Electrons(MC%N_e)%V(:)
!     print*, 'event_photoabsorption 2:', MC%MC_Electrons(MC%N_e)%ti, MC%MC_Electrons(MC%N_e)%t_sc
!     print*, 'event_photoabsorption 3:', Ekin, abs(E_DOS)
!     print*, '--------------------------------------'
   
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
    ! in case we have more particles than spaces in the array, extend the array:
    if (MC%N_h > size(MC%MC_Holes)) call extend_MC_array(MC%MC_Holes)   ! module "MC_general_tools"
    call make_new_particle(MC%MC_Holes(MC%N_h), Ekin=abs(E_DOS), t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_e, &
                                      KOA = KOA, Sh = NSH, valent = valent )    ! module "Objects"

   ! 8) Get new hole's time of the next event:
   call get_hole_flight_time(used_target, numpar, MC%MC_Holes(MC%N_h), MD_supce, E_h)  ! module "MC_general_tools"
   
!    print*, 'event_photoabsorption:', MC%N_e, Ekin, MC%N_h, E_DOS
   
   nullify(Prtcl)
end subroutine event_photoabsorption


subroutine get_photoelectron_emission_angles(Ekin, theta, phi)
   real(8), intent(in) :: Ekin  ! [eV] emitted electron energy
   real(8), intent(out) :: phi, theta   ! electron scattering angles, azymiuthal and polar
   
   ! Sample azimuthal angle randomly within 2Pi:
   phi = sample_phi()   ! module "MC_general_tools"
   
   ! Sample polar angle according to the differential cross section, Eq(2.4) [1]
   call get_photoelectron_polar_angle(Ekin, theta)  ! below

   ! TESTing:
!    theta = acos(0.0d0)
!    print*, 'P:', Ekin, phi, theta
end subroutine get_photoelectron_emission_angles


subroutine get_photoelectron_polar_angle(Ekin, theta)
   real(8), intent(in) :: Ekin  ! [eV] emitted electron energy
   real(8), intent(out) :: theta   ! polar angle of photoelectron emission
   real(8) :: RN, CS_start, CS_end, CS_tot, CS_sampled, eps
   real(8) :: CS_cur, mu_left, mu_right, mu_cur
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-3
   
   if (Ekin < eps) then ! too low energy, assume perpendicular
      mu_cur = 0.0d0
   else ! sample for this energy
      ! Sample the angle:
      call random_number(RN)
      ! Get the starting integration point of the cross section (kernel, without prefactor):
      mu_left = -1.0d0  ! to start with
      mu_right = 1.0d0 ! to start with
      CS_start = Integral_Sauter_Kernel(Ekin, mu_left) ! below
      CS_end = Integral_Sauter_Kernel(Ekin, mu_right)   ! below
      CS_tot = CS_end - CS_start   ! total cross section kernel
      CS_sampled = RN * CS_tot ! sampled cross section kernel
      ! Find the angle that corresponds to the sampled cross section kernel:
      mu_cur = (mu_left + mu_right)/2.0d0
      CS_cur = Integral_Sauter_Kernel(Ekin, mu_cur) - CS_start  ! below
      ! Search by bisection method:
      do while (ABS(CS_cur - CS_sampled)/CS_sampled > eps)
         if (CS_cur > CS_sampled) then
            mu_right = mu_cur
         else
            mu_left = mu_cur
         endif
         mu_cur = (mu_left + mu_right)/2.0d0
         if (abs(mu_left - mu_right) < eps) exit  ! precise enough
         CS_cur = Integral_Sauter_Kernel(Ekin, mu_cur) - CS_start  ! below
      enddo
   endif ! (Ekin < eps)
   
   ! The scattering angle is:
   theta = ACOS(mu_cur)
end subroutine get_photoelectron_polar_angle



pure function Integral_Sauter_Kernel(Ee, mu) result(M)
   real(8) M    ! Kernel of the integral of the diff.CS of photoabsorption by Sauter, Eq.(2.4) p.54 [1]
   real(8), intent(in) :: Ee    ! [eV] energy of emitted photon
   real(8), intent(in) :: mu    ! cos(theta) electron emission angle with respect to the photon directino of motion
   real(8) :: beta, gamma, B, onebmu, b3, mc2
   mc2 = g_me_eV    ! electron mass in [eV]
   beta = sqrt(Ee*(Ee+2.0d0*mc2))/(Ee+mc2)
   b3 = beta*beta*beta
   gamma = 1.0d0 + Ee/mc2
   B = 0.5d0*gamma*(gamma - 1.0d0)*(gamma - 2.0d0)
   onebmu = 1.0d0 - beta*mu
   M = (3.0d0*B*onebmu+2.0d0)/(6.0d0*beta*onebmu*onebmu*onebmu) + &
          (2.0d0*B - 1.0d0 + (1.0d0 - 0.5d0*B)/onebmu - 1.0d0/(3.0d0*onebmu*onebmu))/(b3*onebmu) + &
          B*log(abs(onebmu))/b3  ! integral of Eq.(2.4) [1]
end function Integral_Sauter_Kernel



pure function Sauter_Kernel(Ee, mu) result(M)
   real(8) M    ! Kernel of the integral of the diff.CS of photoabsorption by Sauter, Eq.(2.4) p.54 [1]
   real(8), intent(in) :: Ee    ! [eV] energy of emitted photon
   real(8), intent(in) :: mu    ! cos(theta) electron emission angle with respect to the photon directino of motion
   real(8) :: V, beta, gamma, B, mu2, onebmu
   V =  velosity_from_kinetic_energy(Ee, g_me, afs=.false.)  ! [m/s] module "Relativity"
   beta = beta_factor(V)  ! module "Relativity"
   gamma = gamma_factor(V)  ! module "Relativity"
   B = 0.5d0*gamma*(gamma - 1.0d0)*(gamma - 2.0d0)
   mu2 = mu*mu
   onebmu = 1.0d0 - beta*mu
   M = (1.0d0 - mu2)/(onebmu*onebmu*onebmu*onebmu)*(1.0d0 + B*onebmu)   ! from Eq.(2.4) [1]
end function Sauter_Kernel



subroutine select_shell_photoabsorption(used_target, MC, NOP, KOA, NSH)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   integer, intent(out) :: KOA, NSH ! kind of atom, and the shell number, which absorbs the photon
   !------------------------------------------------------
   type(Target_atoms), pointer :: matter
   type(Photon), pointer :: Prtcl  ! the photon to perform some event
   type(Atom_kind), pointer :: Element
   real(8) :: RN, CS_tot, CS_sampled, CS_cur, CS_sum,  elem_contrib, N_elem
   integer :: i_mat, j, k
   logical :: found_shl
   
   ! Starting default values (valence band):
   KOA = 0
   NSH = 0
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Photons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! Get the total cross section:
   call interpolate_data_single(matter%Ph_absorption_total%E(:), matter%Ph_absorption_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
   ! Get the random number, to sample partial cross section:
   call random_number(RN)   ! intrinsic FORTRAN subroutine
   ! Realized CS:
   CS_sampled = RN*CS_tot

   ! Find which of the partial cross sections correspond to the sampled one:
   CS_sum = 0.0d0   ! to start with
   N_elem = dble(SUM(matter%Elements(:)%percentage))   ! number of elements in this compound material
   ! Scan through all the atoms' shells to find the one:
   SRCH:do j =1, matter%N_Elements	! for each element
      Element => matter%Elements(j)	! all information about this element
      elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
      ! Check the atomic CSs:
      do k = 1, Element%N_shl    ! for all shells in this elements
         if (.not.Element%valent(k)) then ! core orbitals
            ! For each shell, find its partial cross section:
            call interpolate_data_single(Element%Phot_absorption%E,  Element%Phot_absorption%Per_shell(k,:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
            CS_sum = CS_sum + CS_cur*elem_contrib
            ! Check if this shell is the sampled one according to the integral CS:
            call check_shell(CS_sum, CS_sampled, j, k, KOA, NSH, found_shl)    ! module "MC_general_tools"
            if (found_shl) exit SRCH  ! no need to continue the search if we already found the KOA nad NSH
         endif ! (.not.Element%valent(k)) 
      enddo ! k
      ! Check valence band / shells:
      if (allocated(matter%CDF_valence%A)) then    ! valence band
         if (j == 1) then ! add valence band only once - when studying the first element of the material
            call interpolate_data_single(matter%Ph_absorption_valent%E, matter%Ph_absorption_valent%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
            CS_sum = CS_sum + CS_cur
            ! Check if this shell is the sampled one according to the integral CS:
            call check_shell(CS_sum, CS_sampled, 0, 0, KOA, NSH, found_shl)    ! module "MC_general_tools"
            if (found_shl) exit SRCH  ! no need to continue the search if we already found the KOA nad NSH
         endif ! (j == 1)
      else ! valent atomic levels
         do k = 1, Element%N_shl     ! for all shells in this elements
            if (Element%valent(k)) then ! valence orbitals
               call interpolate_data_single(Element%Phot_absorption%E,  Element%Phot_absorption%Per_shell(k,:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
               CS_sum = CS_sum + CS_cur*elem_contrib
               ! Check if this shell is the sampled one according to the integral CS:
               call check_shell(CS_sum, CS_sampled, j, k, KOA, NSH, found_shl)    ! module "MC_general_tools"
               if (found_shl) exit SRCH  ! no need to continue the search if we already found the KOA nad NSH
            endif ! (Element%valent(k)) 
         enddo ! k
      endif ! (allocated(used_target%Material(i)%CDF_valence%A)) 
   enddo SRCH !  j =1, N_elements

   ! Clean up the memory:
   nullify(matter, Prtcl, Element)
end subroutine select_shell_photoabsorption




subroutine find_type_of_photon_event(used_target, numpar, Prtcl, i_type)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Photon), intent(in) :: Prtcl  ! the photon to perform some event
   integer, intent(out) :: i_type   ! type of event: -1=box crossing; 0=boundary; 1=absorption; 2=coherent; 3=Compton;4=e-e+pair creation
   !----------------------------------------------
   real(8) :: eps, MFP, V
   logical :: out
   type(Target_atoms), pointer :: matter
   
   eps = 1.0d-12

!    print*, 'Ph:', Prtcl%Ekin, Prtcl%ti, Prtcl%t_sc, Prtcl%R(:)

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
   else ! it is a scattering-type event
      ! Properties of the target material, inside of which the particle is:
      matter => used_target%Material(Prtcl%in_target)
      
      ! 3) Find out what kind of scattering event it is:
      call find_type_of_scattering(i_type, Prtcl%Ekin,  matter%Ph_absorption_total%E, matter%Ph_absorption_total%Total, &
                    E_array2=matter%Ph_coherent_total%E, CS_array2=matter%Ph_coherent_total%Total, &
                    E_array3=matter%Ph_Compton_total%E, CS_array3=matter%Ph_Compton_total%Total, &
                    E_array4=matter%Ph_pair_total%E, CS_array4=matter%Ph_pair_total%Total) ! module "CS_general_tools"
   endif
   
!    print*, 'Ph2', Prtcl%Ekin, Prtcl%ti, Prtcl%t_sc, Prtcl%R(:), i_type
   
   nullify(matter)
end subroutine find_type_of_photon_event


   
   
end module MC_photon
