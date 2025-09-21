! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018 - 2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains encapsulating routines for Monte Carlo simulations:
module MC
use Universal_constants
use Objects
use Geometries, only: shift_coordinate_system, Cartesian_to_cylindrical, Cartesian_to_spherical
use MC_photon, only: MC_photon_event
use MC_electron, only: MC_electron_event
use MC_positron, only: MC_positron_event
use MC_SHI, only: MC_SHI_event
use MC_hole, only: MC_hole_event
use MC_muon, only: MC_muon_event
use MC_general_tools, only: find_the_target, get_photon_flight_time, get_hole_flight_time, get_SHI_flight_time, &
                            get_positron_flight_time, get_electron_flight_time, get_muon_flight_time, renew_atomic_arrays, &
                            remove_particle_from_MC_array
! use Little_subroutines, only: Find_in_array_monoton


implicit none

 contains



subroutine MC_step(t_cur, dt, MC, used_target, numpar, bunch, MD_supce)
   real(8), intent(in) :: t_cur, dt     ! [fs] current time, and timestep
   type(MC_arrays), dimension(:), intent(inout) :: MC	! all MC arrays for all particles; size = number of iterations
   type(Matter), intent(in) :: used_target  ! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   !---------------------------------
   real(8) :: t_lim
   integer :: iter
   logical :: anything_to_do(numpar%NMC)

   ! Check if there is anything to do within MC module:
   if (.not.numpar%DO_MC) goto 9999    ! if there is no active particle in MC, no need to even start it
   anything_to_do(:) = .false.  ! just to start with
   t_lim = t_cur + dt   ! [fs] allow particles run up to this time (next time step)
   
   ! Perform all MC iterations:
   call Iterate_MC(MC, t_lim, used_target, numpar, bunch, anything_to_do, MD_supce)    ! below
   
   ! Save for the next run information about any active particle within any MC iteration:
   numpar%DO_MC = any(anything_to_do(:))
9999 continue
end subroutine MC_step



subroutine Iterate_MC(MC, t_lim, used_target, numpar, bunch, anything_to_do, MD_supce)
   real(8), intent(in) :: t_lim     ! [fs] end of this timestep
   type(MC_arrays), dimension(:), intent(inout) :: MC	! all MC arrays for all particles; size = number of iterations
   type(Matter), intent(in) :: used_target  ! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   logical, dimension(:), intent(inout) :: anything_to_do   ! are there any active particles?
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   !---------------------------------
   integer :: iter
   real(8), dimension(size(MD_supce%E_e_at_from_MC,1),size(MD_supce%E_e_at_from_MC,2),size(MD_supce%E_e_at_from_MC,3)) :: E_e_at
   real(8), dimension(size(MD_supce%E_h_at_from_MC,1),size(MD_supce%E_h_at_from_MC,2),size(MD_supce%E_h_at_from_MC,3)) :: E_h_at
   real(8), dimension(size(MD_supce%E_p_at_from_MC,1),size(MD_supce%E_p_at_from_MC,2),size(MD_supce%E_p_at_from_MC,3)) :: E_p_at
   real(8), dimension(size(MD_supce%E_mu_at_from_MC,1),size(MD_supce%E_mu_at_from_MC,2),size(MD_supce%E_mu_at_from_MC,3)) :: E_mu_at
   real(8), dimension(size(MD_supce%E_e_from_MC,1),size(MD_supce%E_e_from_MC,2),size(MD_supce%E_e_from_MC,3)) :: E_e
   real(8), dimension(size(MD_supce%E_h_from_MC,1),size(MD_supce%E_h_from_MC,2),size(MD_supce%E_h_from_MC,3)) :: E_h
   
   ! To start with:
   E_e_at = 0.0d0
   E_h_at = 0.0d0
   E_p_at = 0.0d0
   E_mu_at = 0.0d0
   E_e = 0.0d0
   E_h = 0.0d0

   !$omp parallel &
   !$omp private (iter)
   !$omp do schedule(dynamic) reduction(+: E_e_at, E_h_at, E_p_at, E_e, E_h, E_mu_at )
   MC_ITER:do iter = 1, numpar%NMC
      ! Find out if there is any particle that is to be simulated:
      anything_to_do(iter) = ( any(MC(iter)%MC_Photons(:)%active) .or. any(MC(iter)%MC_Electrons(:)%active) &
                                  .or. any(MC(iter)%MC_Positrons(:)%active) .or. any(MC(iter)%MC_Holes(:)%active) &
                                  .or. any(MC(iter)%MC_SHIs(:)%active) .or. any(MC(iter)%MC_Muons(:)%active) )
      ! If there is, then run MC simulation:
      if (anything_to_do(iter)) then    ! do
         call MC_single_iteration(MC(iter), t_lim, used_target, numpar, bunch, MD_supce, E_e_at, E_h_at, E_p_at, E_e, E_h, E_mu_at)    ! below
      endif
   enddo MC_ITER
   !$omp end do
   !$omp end parallel

   ! Save the info exchange between MC and MD into supercell object:
   if (numpar%DO_MD .or. numpar%print_MC_MD_energy) then   ! if user requested MD at all
      MD_supce%E_e_at_from_MC = E_e_at/dble(numpar%NMC)
      MD_supce%E_h_at_from_MC = E_h_at/dble(numpar%NMC)
      MD_supce%E_p_at_from_MC = E_p_at/dble(numpar%NMC)
      MD_supce%E_mu_at_from_MC = E_mu_at/dble(numpar%NMC)
      MD_supce%E_e_from_MC = E_e/dble(numpar%NMC)
      MD_supce%E_h_from_MC = E_h/dble(numpar%NMC)
   endif
end subroutine Iterate_MC


subroutine MC_single_iteration(MC, t_lim, used_target, numpar, bunch, MD_supce, E_e_at, E_h_at, E_p_at, E_e, E_h, E_mu_at)
   type(MC_arrays), intent(inout) :: MC      ! elements of MC array for all particles in one iteration
   real(8), intent(in) :: t_lim                        ! [fs] end of the time step
   type(Matter), intent(in) :: used_target     ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e_at, E_h_at, E_p_at, E_mu_at, E_e, E_h  ! data to pass to MD later
   !----------------------------------------------------
   real(8) :: t_shortest    ! [fs] shortest time of the next event
   real(8) :: t_previous
   integer :: KOP   ! kind of particle
   integer :: NOP, NOP_previous, Ntot_previous   ! index, i.e. number of particle in the array
   integer :: iter_count    ! counter of interations with the same shortest time (for checking progress of the oop)
   
   ! 1) Find the shortest time until the next event, to figure out by which particle it is performed:
   call find_the_shortest_time(MC, t_shortest, KOP, NOP) ! below
   
   t_previous = t_shortest  ! save it for checking next step
   NOP_previous = 0 ! to start with
   iter_count = 0   ! no getting stuck in a loop so far
   ! 2) Run the MC simulation, starting from the particles with the shortest next event time:
   TIM:do while (t_shortest <= t_lim)   ! simulate as long as there are events within this timestep

      select case (KOP)    ! this particle is going to undergo an event
      case (0)  ! photon
         call MC_photon_event(used_target, numpar, MC%N_ph, MC%MC_Photons, NOP, MC, &
                                MD_supce, E_e, E_h)  ! module "MC_photon"
      case (1)  ! electron
         call MC_electron_event(used_target, numpar, MC%N_e, MC%MC_Electrons, NOP, MC, &
                                MD_supce, E_e_at, E_e, E_h)  ! module "MC_electron"
      case (2)  ! positron
         call MC_positron_event(used_target, numpar, MC%N_p, MC%MC_Positrons, NOP, MC, &
                                MD_supce, E_p_at, E_e, E_h)  ! module "MC_positron"
      case (3)  ! ion
         call MC_SHI_event(used_target, numpar, bunch, MC%N_SHI, MC%MC_SHIs, NOP, MC, &
                                MD_supce, E_e, E_h)  ! module "MC_SHI"
      case (4)  ! hole
         call MC_hole_event(used_target, numpar, MC%N_h, MC%MC_Holes, NOP, MC, &
                                MD_supce, E_h_at, E_e, E_h)  ! module "MC_hole"
      case (5)  ! muon
         call MC_muon_event(used_target, numpar, MC%N_mu, MC%MC_Muons, NOP, MC, &
                                MD_supce, E_mu_at, E_e, E_h)  ! module "MC_muon"
      endselect
      ! Check if there are still particles to perform an event within this timestep:
      call find_the_shortest_time(MC, t_shortest, KOP, NOP) ! below



      ! Check for consistency:
      if ( (t_shortest <= t_previous ) ) then ! somehow, the next step is not advancing after the previous one
         if ((Ntot_previous == MC%N_ph) .and. ((NOP_previous == NOP) .or. (t_shortest <= t_previous-1.0d-3))) then ! either it's the same particle, or moved back in time
            iter_count = iter_count + 1    ! count repeated attempts
         endif

         if ((iter_count > 1) .and. .not.( (numpar%Ph_att_eff >= 0.0d0) .and. (KOP == 0) ) ) then   ! inform about an Error
            write(*,'(a,i6,es,es,i6,i6)') 'MC step error', iter_count, t_previous, t_shortest, KOP, NOP
            write(*,'(i3, es, i6)') MC%MC_Electrons(NOP)%generation, MC%MC_Electrons(NOP)%Ekin, MC%N_e
            if (iter_count > 10) then
               write(*,'(a)') 'We got stuck in a time loop. Emergency exit.'
               !exit TIM
               select case (KOP)    ! this particle needs to be removed:
               case (0)  ! photon
                  call remove_particle_from_MC_array(MC%N_ph, NOP, MC%MC_Photons)  ! module "MC_general_tools"
               case (1)  ! electron
                  call remove_particle_from_MC_array(MC%N_e, NOP, MC%MC_Electrons)  ! module "MC_general_tools"
               case (2)  ! positron
                  call remove_particle_from_MC_array(MC%N_p, NOP, MC%MC_Positrons)  ! module "MC_general_tools"
               case (3)  ! ion
                  call remove_particle_from_MC_array(MC%N_SHI, NOP, MC%MC_SHIs)  ! module "MC_general_tools"
               case (4)  ! hole
                  call remove_particle_from_MC_array(MC%N_h, NOP, MC%MC_Holes)  ! module "MC_general_tools"
               case (5)  ! muon
                  call remove_particle_from_MC_array(MC%N_mu, NOP, MC%MC_Muons)  ! module "MC_general_tools"
               endselect
            endif
         endif
      else
         t_previous = t_shortest  ! save it for checking next step
         iter_count = 0
         NOP_previous = NOP
         Ntot_previous = MC%N_ph
      endif
   enddo TIM ! while (t_shortest <= t_lim)
end subroutine MC_single_iteration


subroutine find_the_shortest_time(MC, t_shortest, KOP, NOP)
   type(MC_arrays), intent(in) :: MC	! elements of all MC arrays: photons, electrons and holes; size equals to number of iterations
   real(8), intent(out) :: t_shortest  ! [fs] the shortest next collision time
   integer, intent(out) :: KOP  ! index of the kind of the particle that is next to perform an event
   integer, intent(out) :: NOP  ! index of the particle
   ! Reminder: particle kinds: 0=photon, 1=electron, 2=positron, 3=SHI, 4=hole, 5=muon
   real(8) :: t_e, t_i, t_ph, t_pos, t_h, t_muon    ! for each kind of particles defined in the code
   integer :: NOP_ph, NOP_e, NOP_p, NOP_SHI, NOP_h, NOP_mu
   
   ! The shortest time till next collision among only active photons:
   NOP_ph = transfer( (/ minloc(MC%MC_Photons(:)%ti, MASK=MC%MC_Photons(:)%active) /), NOP)
   if (NOP_ph > 0) then
      t_ph = MC%MC_Photons(NOP_ph)%ti
   else
      t_ph = 1.0d25  ! no collisions if no active photons
   endif
   NOP = NOP_ph ! start with assuming photon
   
   ! The shortest time till next collision among only active electrons:
   NOP_e = transfer( (/ minloc(MC%MC_Electrons(:)%ti, MASK=MC%MC_Electrons(:)%active) /), NOP)
   if (NOP_e > 0) then
      t_e = MC%MC_Electrons(NOP_e)%ti
   else
      t_e = 1.0d25  ! no collisions if no active electrons
   endif
!    print*, "find_the_shortest_time", NOP_e, MC%N_e, t_e, size(MC%MC_Electrons(:)), MC%MC_Electrons(NOP_e)%active, MC%MC_Electrons(NOP_e)%ti
   
   ! The shortest time till next collision among only active positrons:
   NOP_p = transfer( (/ minloc(MC%MC_Positrons(:)%ti, MASK=MC%MC_Positrons(:)%active) /), NOP)
   if (NOP_p > 0) then
      t_pos = MC%MC_Positrons(NOP_p)%ti
   else
      t_pos = 1.0d25  ! no collisions if no active positrons
   endif
   
   ! The shortest time till next collision among only active holes:
   NOP_h = transfer( (/ minloc(MC%MC_Holes(:)%ti, MASK=MC%MC_Holes(:)%active) /), NOP)
   if (NOP_h > 0) then
      t_h = MC%MC_Holes(NOP_h)%ti
   else
      t_h = 1.0d25  ! no collisions if no active holes
   endif
   
   ! The shortest time till next collision among only active ions:
   NOP_SHI = transfer( (/ minloc(MC%MC_SHIs(:)%ti, MASK=MC%MC_SHIs(:)%active) /), NOP)
   if (NOP_SHI > 0) then
      t_i = MC%MC_SHIs(NOP_SHI)%ti
   else
      t_i = 1.0d25  ! no collisions if no active SHI
   endif

   ! The shortest time till next collision among only active mouns:
   NOP_mu = transfer( (/ minloc(MC%MC_Muons(:)%ti, MASK=MC%MC_Muons(:)%active) /), NOP)
   if (NOP_mu > 0) then
      t_muon = MC%MC_Muons(NOP_mu)%ti
   else
      t_muon = 1.0d25  ! no collisions if no active muons
   endif
   
   ! Chose the shortest time among the types of particles:
   ! To start with:
   KOP = 0  ! photon
   t_shortest = t_ph
   
   if (t_e < t_shortest) then   ! electron time is shorter
      KOP = 1   ! electron
      t_shortest = t_e
      NOP = NOP_e
   endif
   
   if (t_pos < t_shortest) then   ! positron time is shorter
      KOP = 2   ! positron
      t_shortest = t_pos
      NOP = NOP_p
   endif
   
   if (t_i < t_shortest) then   ! SHI time is shorter
      KOP = 3   ! ion
      t_shortest = t_i
      NOP = NOP_SHI
   endif
   
   if (t_h < t_shortest) then   ! hole time is shorter
      KOP = 4   ! hole
      t_shortest = t_h
      NOP = NOP_h
   endif

   if (t_muon < t_shortest) then   ! muon time is shorter
      KOP = 5   ! muon
      t_shortest = t_muon
      NOP = NOP_mu
   endif
end subroutine find_the_shortest_time


!-------------------------------------------------------------
! Preparing before the first MC run:
subroutine prepare_MC_run(used_target, MC, numpar, MD_supce, E_e, E_h)
   type(Matter), intent(in) :: used_target  ! parameters of the target
   type(MC_arrays), dimension(:), intent(inout) :: MC   ! all MC arrays for particles: photons, electrons and holes; size equals to number of iterations
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h ! data to pass to MD later
   !----------------
   integer :: iter
   logical :: anything_to_do(numpar%NMC)
   
   ! Chech if there is anything to do within MC module:
   if (.not.numpar%DO_MC) goto 9998    ! if there is no active particle in MC, no need to even start it
   anything_to_do(:) = .false.  ! just to start with
   
   !$omp parallel &
   !$omp private (iter)
   !$omp do
   MC_ITER:do iter = 1, numpar%NMC
      ! Find out if there is any particle that is to be simulated:
      anything_to_do(iter) = ( any(MC(iter)%MC_Photons(:)%active) .or. any(MC(iter)%MC_Electrons(:)%active) &
                                  .or. any(MC(iter)%MC_Positrons(:)%active) .or. any(MC(iter)%MC_Holes(:)%active) &
                                  .or. any(MC(iter)%MC_SHIs(:)%active) .or. any(MC(iter)%MC_Muons(:)%active) )
      ! If there is, then set the free flight distances in this iteration:
      if (anything_to_do(iter)) then
         call set_free_flights(used_target, numpar, MC(iter), MD_supce, E_e, E_h)    ! below
      endif
   enddo MC_ITER
   !$omp end do
   !$omp end parallel
   
   ! Save for the next run information about any active particle within any MC iteration:
   numpar%DO_MC = any(anything_to_do(:))
9998 continue
end subroutine prepare_MC_run


subroutine set_free_flights(used_target, numpar, MC, MD_supce, E_e, E_h)
   type(Matter), intent(in) :: used_target  ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(MC_arrays), intent(inout) :: MC   ! MC arrays for particles: photons, electrons and holes within one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h ! data to pass to MD later
   integer :: i
   ! Do for all photons:
   if (any(MC%MC_Photons(:)%active) ) then    ! if there is at least one
      do i = 1, MC%N_ph
         ! Find inside of which target this particle enters:
         call find_the_target(used_target, MC%MC_Photons(i)) ! module "MC_general_tools"
!           print*, 'set_free_flights 0:', i, MC%MC_Photons(i)%ti, MC%MC_Photons(i)%t0, MC%MC_Photons(i)%in_target 
         ! Get its flight time within this target:
         call get_photon_flight_time(used_target, numpar, MC%MC_Photons(i))  ! module "MC_general_tools"
!          print*, 'set_free_flights', i, MC%MC_Photons(i)%ti, MC%MC_Photons(i)%t0
      enddo
   endif
   ! Do for all electrons:
   if (any(MC%MC_Electrons(:)%active)  ) then    ! if there is at least one
      do i = 1, MC%N_e
         ! Find inside of which target this particle enters:
         call find_the_target(used_target, MC%MC_Electrons(i))  ! module "MC_general_tools"
         ! Get its flight time within this target:
         call get_electron_flight_time(used_target, numpar, MC%MC_Electrons(i), MD_supce, E_e)  ! module "MC_electron"
      enddo
   endif
   ! Do for all positrons:
   if (any(MC%MC_Positrons(:)%active) ) then    ! if there is at least one
      do i = 1, MC%N_p
         ! Find inside of which target this particle enters:
         call find_the_target(used_target, MC%MC_Positrons(i))  ! module "MC_general_tools"
         ! Get its flight time within this target:
         call get_positron_flight_time(used_target, numpar, MC%MC_Positrons(i))  ! module "MC_positron"
      enddo
   endif
   ! Do for all holes:
   if (any(MC%MC_Holes(:)%active) ) then    ! if there is at least one
      do i = 1, MC%N_h
         ! Find inside of which target this particle enters:
         call find_the_target(used_target, MC%MC_Holes(i))  ! module "MC_general_tools"
         ! Get its flight time within this target:
         call get_hole_flight_time(used_target, numpar, MC%MC_Holes(i), MD_supce, E_h)  ! module "MC_hole"
      enddo
   endif
   ! Do for all SHIs:
   if (any(MC%MC_SHIs(:)%active) ) then    ! if there is at least one
      do i = 1, MC%N_SHI
         ! Find inside of which target this particle enters:
         call find_the_target(used_target, MC%MC_SHIs(i))  ! module "MC_general_tools"
         ! Get its flight time within this target:
         call get_SHI_flight_time(used_target, numpar, MC%MC_SHIs(i))  ! module "MC_SHI"
      enddo
   endif
   ! Do for all muons:
   if (any(MC%MC_Muons(:)%active) ) then    ! if there is at least one
      do i = 1, MC%N_mu
         ! Find inside of which target this particle enters:
         call find_the_target(used_target, MC%MC_Muons(i))  ! module "MC_general_tools"
         ! Get its flight time within this target:
         call get_muon_flight_time(used_target, numpar, MC%MC_Muons(i))  ! module "MC_muons"
      enddo
   endif
end subroutine set_free_flights


end module MC
