! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018-2019
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains definition of default variables and initial conditions

MODULE Initial_conditions
use Universal_constants
use Objects
use Periodic_table, only : Atomic_data, Read_from_periodic_table
use Little_subroutines, only : sample_Gaussian
use Relativity, only : velosity_from_kinetic_energy
use Geometries, only : Spherical_to_cartesian
use Read_input_data, only : m_input_folder, m_databases
use SHI_charge_state, only : Equilibrium_charge_SHI

 implicit none
 
integer :: m_MC_size

parameter (m_MC_size = 1)    ! default starting value of particles arrays for MC
 
 contains


! Set the right starting time, depending on whether we use the pulse or not:
pure function set_starting_time(bunch, t_start) result (tim)
   real(8) :: tim ! starting time step [fs]
   type(Radiation_param), dimension(:), intent(in) :: bunch ! incomming radiation
   real(8), intent(in) :: t_start   ! [fs] predefined starting time
   tim = dble(CEILING(min(minval(bunch(:)%t-bunch(:)%FWHM*5.0d0), t_start)))  ! [fs]
end function set_starting_time


subroutine Set_initial_parameters(numpar, used_target, bunch, MC)
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   type(Matter), intent(in):: used_target	! parameters of the target
   type(Radiation_param), dimension(:), allocatable, intent(inout) :: bunch	! incomming radiation
   type(MC_arrays), dimension(:), intent(inout) :: MC	! all MC arrays for particles: photons, electrons and holes; size equals to number of iterations
!    type(Atom), dimension(:), allocatable, intent(inout) :: MD_atoms	! all atoms in MD as objects
   !----------------------------------------
   type(atomic_data), dimension(:), allocatable :: Periodic_table ! this is an internal module variable
   integer :: N_bunch, iter, ibunch, ipart, i_SHI
   integer :: INFO  ! 0=file read well; 1=no file; 2=couldn't open; 3=error while reading
   character(100) :: error_message , Path
   
   print*, 'Setting initial conditions...'
   
   ! Read the periodic table, in case we need it:
   Path = trim(adjustl(m_input_folder))//numpar%path_sep//trim(adjustl(m_databases)) ! where to find the periodic table
   call Read_from_periodic_table(Path, numpar%path_sep, INFO, error_message, 1, Periodic_table_out=Periodic_table) ! module "Periodic_table"
   
   N_bunch = size(bunch(:)) ! that's how many different radiation types are incoming on the target
   ! Define them within MC:
   if (numpar%NMC>0) then
      !$omp parallel private (iter, ibunch, ipart, i_SHI)
      !$omp do
      do iter = 1, numpar%NMC  ! in each MC iteration
         ! To start counting active particles:
         MC(iter)%N_ph = 0  ! number of active photons
         MC(iter)%N_e = 0    ! number of active electrons
         MC(iter)%N_p = 0    ! number of active positrons
         MC(iter)%N_h = 0    ! number of active holes
         MC(iter)%N_SHI = 0  ! number of active SHIs
         MC(iter)%N_at_nrg = 0  ! number of active atomic events of energy transfer
         i_SHI = 0 ! to count typies of SHI
      
         ! Start setting the parameters into MC arrays:
         do ibunch = 1, N_bunch
            ! For all particles in the bunch/pulse:
            do ipart = 1, bunch(ibunch)%NOP
               call  define_projectile(numpar, used_target, bunch, MC, iter, ibunch, ipart, i_SHI, Periodic_table) ! below
            enddo
         enddo
      enddo
      !$omp end do
      !$omp end parallel
   else ! just to get MFPs and parameters:
      do ibunch = 1, N_bunch
         ! For all particles in the bunch/pulse:
         do ipart = 1, bunch(ibunch)%NOP
            call  define_projectile(numpar, used_target, bunch, MC, 1, ibunch, ipart, i_SHI, Periodic_table) ! below
         enddo
      enddo
   endif
   ! Once we read all the information into arrays, no need to keep the projectiles info:
   !deallocate(bunch, Periodic_table)
   deallocate(Periodic_table)
   
   print*, 'Initial conditions all set, ready to start simulation.'
end subroutine Set_initial_parameters


subroutine define_projectile(numpar, used_target, bunch, MC, iter, ibunch, iparticle, i_SHI, Periodic_table)
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   type(Matter), intent(in), target :: used_target	! parameters of the target
   type(Radiation_param), dimension(:), intent(inout), target :: bunch	! incomming radiation
   type(MC_arrays), dimension(:), intent(inout) :: MC	! all MC arrays for particles: photons, electrons and holes; size equals to number of iterations
   integer, intent(in) :: iter, ibunch, iparticle    ! nuber of iteration, number of the bunch, number of incoming particle in the bunch
   integer, intent(inout) :: i_SHI  ! typies of SHI
   type(atomic_data), dimension(:), intent(in) :: Periodic_table ! this is an internal module variable
   !----------------------
   integer, pointer :: KOP
   real(8) :: V
   
   ! Kind of particle: 0=photon, 1=electron, 2=positron, 3=SHI, 4=hole
   KOP => bunch(ibunch)%KOP
   ! Sort data from the projectile/bunch into the MC arrays to start MC modelling:
   select case (KOP)
   case default ! photon
      MC(iter)%N_ph = MC(iter)%N_ph + 1 ! one more photon among incoming particles
      MC(iter)%MC_Photons(MC(iter)%N_ph)%active = .true.
      MC(iter)%MC_Photons(MC(iter)%N_ph)%generation = 0 ! incomming particle
      ! Sample energy around the given one according to gaussian shape:
      MC(iter)%MC_Photons(MC(iter)%N_ph)%Ekin = sample_Gaussian(bunch(ibunch)%E, bunch(ibunch)%E_spread)    ! module "Little_subroutines"
      ! Sample arrival time around the given one according to gaussian shape:
      MC(iter)%MC_Photons(MC(iter)%N_ph)%t0 = sample_Gaussian(bunch(ibunch)%t, bunch(ibunch)%FWHM+1.0d-6)    ! module "Little_subroutines"
      ! Sample coordinates of impact around the given one according to gaussian shape:
      MC(iter)%MC_Photons(MC(iter)%N_ph)%R(1) = sample_Gaussian(bunch(ibunch)%R(1)+1.0d-6, bunch(ibunch)%R_spread(1))    ! module "Little_subroutines"
      MC(iter)%MC_Photons(MC(iter)%N_ph)%R(2) = sample_Gaussian(bunch(ibunch)%R(2)+1.0d-6, bunch(ibunch)%R_spread(2))    ! module "Little_subroutines"
      MC(iter)%MC_Photons(MC(iter)%N_ph)%R(3) = sample_Gaussian(bunch(ibunch)%R(3)+1.0d-6, bunch(ibunch)%R_spread(3))    ! module "Little_subroutines"
      ! Absolute value of velosity:
      V = velosity_from_kinetic_energy(MC(iter)%MC_Photons(MC(iter)%N_ph)%Ekin, 0.0d0)  ! [A/fs] module "Relativity" [m/s]
!       V = V * g_ms2Afs    ! [m/s] -> [A/fs]
      ! Velosity projections in cartesian coordinates:
      call Spherical_to_cartesian(V, bunch(ibunch)%theta, bunch(ibunch)%phi, MC(iter)%MC_Photons(MC(iter)%N_ph)%V(1), MC(iter)%MC_Photons(MC(iter)%N_ph)%V(2), MC(iter)%MC_Photons(MC(iter)%N_ph)%V(3))   ! module "Geometries"
      ! Data on the previous time step:
      MC(iter)%MC_Photons(MC(iter)%N_ph)%V0(:) = MC(iter)%MC_Photons(MC(iter)%N_ph)%V(:)
      MC(iter)%MC_Photons(MC(iter)%N_ph)%R0(:) = MC(iter)%MC_Photons(MC(iter)%N_ph)%R(:) - numpar%dt_MD * MC(iter)%MC_Photons(MC(iter)%N_ph)%V(:)   ! making one step back in time
   case (1) ! electron
      MC(iter)%N_e = MC(iter)%N_e + 1 ! one more electron among incoming particles
      MC(iter)%MC_Electrons(MC(iter)%N_e)%active = .true.
      MC(iter)%MC_Electrons(MC(iter)%N_e)%generation = 0 ! incomming particle
      ! Sample energy around the given one according to gaussian shape:
      MC(iter)%MC_Electrons(MC(iter)%N_e)%Ekin = sample_Gaussian(bunch(ibunch)%E, bunch(ibunch)%E_spread)    ! module "Little_subroutines"
      ! Sample arrival time around the given one according to gaussian shape:
      MC(iter)%MC_Electrons(MC(iter)%N_e)%t0 = sample_Gaussian(bunch(ibunch)%t, bunch(ibunch)%FWHM+1.0d-6)    ! module "Little_subroutines"
      ! Sample coordinates of impact around the given one according to gaussian shape:
      MC(iter)%MC_Electrons(MC(iter)%N_e)%R(1) = sample_Gaussian(bunch(ibunch)%R(1), bunch(ibunch)%R_spread(1))    ! module "Little_subroutines"
      MC(iter)%MC_Electrons(MC(iter)%N_e)%R(2) = sample_Gaussian(bunch(ibunch)%R(2), bunch(ibunch)%R_spread(2))    ! module "Little_subroutines"
      MC(iter)%MC_Electrons(MC(iter)%N_e)%R(3) = sample_Gaussian(bunch(ibunch)%R(3), bunch(ibunch)%R_spread(3))    ! module "Little_subroutines"
      ! Absolute value of velosity:
      V = velosity_from_kinetic_energy(MC(iter)%MC_Electrons(MC(iter)%N_e)%Ekin, g_me)  ! [A/fs] module "Relativity"
!       V = V * g_ms2Afs    ! [m/s] -> [A/fs]
      ! Velosity projections in cartesian coordinates:
      call Spherical_to_cartesian(V, bunch(ibunch)%theta, bunch(ibunch)%phi, MC(iter)%MC_Electrons(MC(iter)%N_e)%V(1), MC(iter)%MC_Electrons(MC(iter)%N_e)%V(2), MC(iter)%MC_Electrons(MC(iter)%N_e)%V(3))   ! module "Geometries"
      ! Data on the previous time step:
      MC(iter)%MC_Electrons(MC(iter)%N_e)%V0(:) = MC(iter)%MC_Electrons(MC(iter)%N_e)%V(:)
      MC(iter)%MC_Electrons(MC(iter)%N_e)%R0(:) = MC(iter)%MC_Electrons(MC(iter)%N_e)%R(:) - numpar%dt_MD * MC(iter)%MC_Electrons(MC(iter)%N_e)%V(:)   ! making one step back in time
   
   case (2) ! positron
      MC(iter)%N_p = MC(iter)%N_p + 1 ! one more positron among incoming particles
      MC(iter)%MC_Positrons(MC(iter)%N_p)%active = .true.
      MC(iter)%MC_Positrons(MC(iter)%N_p)%generation = 0 ! incomming particle
      ! Sample energy around the given one according to gaussian shape:
      MC(iter)%MC_Positrons(MC(iter)%N_p)%Ekin = sample_Gaussian(bunch(ibunch)%E, bunch(ibunch)%E_spread)    ! module "Little_subroutines"
      ! Sample arrival time around the given one according to gaussian shape:
      MC(iter)%MC_Positrons(MC(iter)%N_p)%t0 = sample_Gaussian(bunch(ibunch)%t, bunch(ibunch)%FWHM+1.0d-6)    ! module "Little_subroutines"
      ! Sample coordinates of impact around the given one according to gaussian shape:
      MC(iter)%MC_Positrons(MC(iter)%N_p)%R(1) = sample_Gaussian(bunch(ibunch)%R(1), bunch(ibunch)%R_spread(1))    ! module "Little_subroutines"
      MC(iter)%MC_Positrons(MC(iter)%N_p)%R(2) = sample_Gaussian(bunch(ibunch)%R(2), bunch(ibunch)%R_spread(2))    ! module "Little_subroutines"
      MC(iter)%MC_Positrons(MC(iter)%N_p)%R(3) = sample_Gaussian(bunch(ibunch)%R(3), bunch(ibunch)%R_spread(3))    ! module "Little_subroutines"
      ! Absolute value of velosity:
      V = velosity_from_kinetic_energy(MC(iter)%MC_Positrons(MC(iter)%N_p)%Ekin, g_me)  ! [A/fs] module "Relativity"
!       V = V * g_ms2Afs    ! [m/s] -> [A/fs]
      ! Velosity projections in cartesian coordinates:
      call Spherical_to_cartesian(V, bunch(ibunch)%theta, bunch(ibunch)%phi, MC(iter)%MC_Positrons(MC(iter)%N_p)%V(1), MC(iter)%MC_Positrons(MC(iter)%N_p)%V(2), MC(iter)%MC_Positrons(MC(iter)%N_p)%V(3))   ! module "Geometries"
      ! Data on the previous time step:
      MC(iter)%MC_Positrons(MC(iter)%N_p)%V0(:) = MC(iter)%MC_Positrons(MC(iter)%N_p)%V(:)
      MC(iter)%MC_Positrons(MC(iter)%N_p)%R0(:) = MC(iter)%MC_Positrons(MC(iter)%N_p)%R(:) - numpar%dt_MD * MC(iter)%MC_Positrons(MC(iter)%N_p)%V(:)   ! making one step back in time
   
   case (3) ! SHI
      if (.not.repeated_type(bunch, ibunch)) i_SHI = i_SHI + 1   ! new type of SHI
      MC(iter)%N_SHI = MC(iter)%N_SHI + 1 ! one more SHI among incoming particles
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%active = .true.
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%generation = 0 ! incomming particle
      ! Sample energy around the given one according to gaussian shape:
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Ekin = sample_Gaussian(bunch(ibunch)%E, bunch(ibunch)%E_spread)    ! module "Little_subroutines"
      ! Sample arrival time around the given one according to gaussian shape:
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%t0 = sample_Gaussian(bunch(ibunch)%t, bunch(ibunch)%FWHM+1.0d-6)    ! module "Little_subroutines"
      ! Sample coordinates of impact around the given one according to gaussian shape:
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%R(1) = sample_Gaussian(bunch(ibunch)%R(1), bunch(ibunch)%R_spread(1))    ! module "Little_subroutines"
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%R(2) = sample_Gaussian(bunch(ibunch)%R(2), bunch(ibunch)%R_spread(2))    ! module "Little_subroutines"
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%R(3) = sample_Gaussian(bunch(ibunch)%R(3), bunch(ibunch)%R_spread(3))    ! module "Little_subroutines"
      ! For SHI, it's atomic number, effective charge etc:
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Z = bunch(ibunch)%Z
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Name = Periodic_table(MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Z)%Name
      bunch(ibunch)%Name = MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Name ! also save the name here
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Mass = Periodic_table(MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Z)%Mass
      
!       print*, 'bunch(ibunch)%Name ', bunch(ibunch)%Name
!       pause 
      
      if (bunch(ibunch)%Meff <= 0.0d0) then
         MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Meff = MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Mass    ! use default value
         bunch(ibunch)%Meff = MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Mass    ! also replace with a default value
      else
         MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Meff = bunch(ibunch)%Meff ! [amu]
      endif
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Mass = MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Meff * g_amu   ! [amu] -> [kg]
      ! Absolute value of velosity:
      V = velosity_from_kinetic_energy(MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Ekin, MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Mass)  ! [Afs/] module "Relativity"
!       V = V * g_ms2Afs    ! [m/s] -> [A/fs]
      ! Velosity projections in cartesian coordinates:
      call Spherical_to_cartesian(V, bunch(ibunch)%theta, bunch(ibunch)%phi, MC(iter)%MC_SHIs(MC(iter)%N_SHI)%V(1), &
                MC(iter)%MC_SHIs(MC(iter)%N_SHI)%V(2), MC(iter)%MC_SHIs(MC(iter)%N_SHI)%V(3))   ! module "Geometries"
                
      ! Data on the previous time step:
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%V0(:) = MC(iter)%MC_SHIs(MC(iter)%N_SHI)%V(:)
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%R0(:) = MC(iter)%MC_SHIs(MC(iter)%N_SHI)%R(:) - numpar%dt_MD * MC(iter)%MC_SHIs(MC(iter)%N_SHI)%V(:)   ! making one step back in time
      ! Assumed that the SHI is within the target number 1 (TO CHECK AND CHANGE LATER) : 
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Zeff = Equilibrium_charge_SHI(MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Ekin, &
                MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Mass, dble(MC(iter)%MC_SHIs(MC(iter)%N_SHI)%Z), &
                used_target%Material(1)%Mean_Z, numpar%SHI_ch_st, bunch(ibunch)%Zeff) ! module "SHI_charge_state"
      
      ! Index of the kind of atom of this SHI:
      MC(iter)%MC_SHIs(MC(iter)%N_SHI)%KOA = i_SHI
   case (4) ! hole
      MC(iter)%N_h = MC(iter)%N_h + 1 ! one more hole among incoming particles
      MC(iter)%MC_Holes(MC(iter)%N_h)%active = .true.
      MC(iter)%MC_Holes(MC(iter)%N_h)%generation = 0 ! incomming particle
      ! Sample energy around the given one according to gaussian shape:
      MC(iter)%MC_Holes(MC(iter)%N_h)%Ekin = sample_Gaussian(bunch(ibunch)%E, bunch(ibunch)%E_spread)    ! module "Little_subroutines"
      ! Sample arrival time around the given one according to gaussian shape:
      MC(iter)%MC_Holes(MC(iter)%N_h)%t0 = sample_Gaussian(bunch(ibunch)%t, bunch(ibunch)%FWHM+1.0d-6)    ! module "Little_subroutines"
      ! Sample coordinates of impact around the given one according to gaussian shape:
      MC(iter)%MC_Holes(MC(iter)%N_h)%R(1) = sample_Gaussian(bunch(ibunch)%R(1), bunch(ibunch)%R_spread(1))    ! module "Little_subroutines"
      MC(iter)%MC_Holes(MC(iter)%N_h)%R(2) = sample_Gaussian(bunch(ibunch)%R(2), bunch(ibunch)%R_spread(2))    ! module "Little_subroutines"
      MC(iter)%MC_Holes(MC(iter)%N_h)%R(3) = sample_Gaussian(bunch(ibunch)%R(3), bunch(ibunch)%R_spread(3))    ! module "Little_subroutines"
      ! Absolute value of velosity:
      V = velosity_from_kinetic_energy(MC(iter)%MC_Holes(MC(iter)%N_h)%Ekin, g_me)  ! [A/fs] module "Relativity"
!       V = V * g_ms2Afs    ! [m/s] -> [A/fs]
      ! Velosity projections in cartesian coordinates:
      call Spherical_to_cartesian(V, bunch(ibunch)%theta, bunch(ibunch)%phi, MC(iter)%MC_Holes(MC(iter)%N_h)%V(1), MC(iter)%MC_Holes(MC(iter)%N_h)%V(2), MC(iter)%MC_Holes(MC(iter)%N_h)%V(3))   ! module "Geometries"
      ! Data on the previous time step:
      MC(iter)%MC_Holes(MC(iter)%N_h)%V0(:) = MC(iter)%MC_Holes(MC(iter)%N_h)%V(:)
      MC(iter)%MC_Holes(MC(iter)%N_h)%R0(:) = MC(iter)%MC_Holes(MC(iter)%N_h)%R(:) - numpar%dt_MD * MC(iter)%MC_Holes(MC(iter)%N_h)%V(:)   ! making one step back in time
      ! For hole, kind of atom and shell it sits in:
      MC(iter)%MC_Holes(MC(iter)%N_h)%KOA = bunch(ibunch)%KOA
      MC(iter)%MC_Holes(MC(iter)%N_h)%Sh = bunch(ibunch)%Sh

   case (5) ! muon
      MC(iter)%N_mu = MC(iter)%N_mu + 1 ! one more muon among incoming particles
      MC(iter)%MC_Muons(MC(iter)%N_mu)%active = .true.
      MC(iter)%MC_Muons(MC(iter)%N_mu)%generation = 0 ! incomming particle
      ! Sample energy around the given one according to gaussian shape:
      MC(iter)%MC_Muons(MC(iter)%N_mu)%Ekin = sample_Gaussian(bunch(ibunch)%E, bunch(ibunch)%E_spread)    ! module "Little_subroutines"
      ! Sample arrival time around the given one according to gaussian shape:
      MC(iter)%MC_Muons(MC(iter)%N_mu)%t0 = sample_Gaussian(bunch(ibunch)%t, bunch(ibunch)%FWHM+1.0d-6)    ! module "Little_subroutines"
      ! Sample coordinates of impact around the given one according to gaussian shape:
      MC(iter)%MC_Muons(MC(iter)%N_mu)%R(1) = sample_Gaussian(bunch(ibunch)%R(1), bunch(ibunch)%R_spread(1))    ! module "Little_subroutines"
      MC(iter)%MC_Muons(MC(iter)%N_mu)%R(2) = sample_Gaussian(bunch(ibunch)%R(2), bunch(ibunch)%R_spread(2))    ! module "Little_subroutines"
      MC(iter)%MC_Muons(MC(iter)%N_mu)%R(3) = sample_Gaussian(bunch(ibunch)%R(3), bunch(ibunch)%R_spread(3))    ! module "Little_subroutines"
      ! Absolute value of velosity:
      V = velosity_from_kinetic_energy(MC(iter)%MC_Muons(MC(iter)%N_mu)%Ekin, g_M_muon)  ! [A/fs] module "Relativity"
!       V = V * g_ms2Afs    ! [m/s] -> [A/fs]
      ! Velosity projections in cartesian coordinates:
      call Spherical_to_cartesian(V, bunch(ibunch)%theta, bunch(ibunch)%phi, &
           MC(iter)%MC_Muons(MC(iter)%N_mu)%V(1), MC(iter)%MC_Muons(MC(iter)%N_mu)%V(2), MC(iter)%MC_Muons(MC(iter)%N_mu)%V(3))   ! module "Geometries"
      ! Data on the previous time step:
      MC(iter)%MC_Muons(MC(iter)%N_mu)%V0(:) = MC(iter)%MC_Muons(MC(iter)%N_mu)%V(:)
      MC(iter)%MC_Muons(MC(iter)%N_mu)%R0(:) = MC(iter)%MC_Muons(MC(iter)%N_mu)%R(:) - numpar%dt_MD * MC(iter)%MC_Muons(MC(iter)%N_mu)%V(:)   ! making one step back in time

   end select
   
   nullify(KOP)
end subroutine define_projectile



function count_types_of_SHIs(bunch) result(N_SHI)
   integer :: N_SHI ! how many different types of ions
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   integer :: N_bunch, i
   N_bunch = size(bunch)
   N_SHI = 0 ! no ions to be analysed
   if (N_bunch >= 1) then
      do i = 1, N_bunch
         if (bunch(i)%KOP == 3) N_SHI = N_SHI + 1   ! it is an SHI
      enddo
   endif
end function count_types_of_SHIs



pure function count_types_of_SHIs_OLD(bunch) result(N_SHI)
   integer :: N_SHI ! how many different types of ions
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   integer :: ibunch, N_bunch, i, j
   N_bunch = size(bunch)
   if (N_bunch < 1) then
      N_SHI = 0 ! no ions to be analysed
   elseif (N_bunch == 1) then
      if (bunch(1)%KOP == 3) N_SHI = 1
   else ! count how many are different
      if (bunch(1)%KOP == 3) N_SHI = 1
      do i = 2, N_bunch
         if (bunch(i)%KOP == 3) N_SHI = N_SHI + 1   ! it is an SHI
         if ( repeated_type(bunch,i) ) N_SHI = N_SHI - 1  ! the same type of Ion -> exclude
      enddo ! i
   endif
end function count_types_of_SHIs_OLD


pure function repeated_type(bunch, i) result(the_same)   ! check if the ion number "i" is the same element as was already in the bunch
   logical :: the_same
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   integer, intent(in) :: i
   integer :: j, KOP
   the_same = .false.
   if (i > 1) then
      SM:do j = 1, i-1
         if (bunch(j)%KOP == 3) then ! it was an ion
            if (bunch(i)%Z == bunch(j)%Z) then
               the_same = .true.  ! found the same ion
               exit SM
            endif
         endif
      enddo SM
   endif
end function repeated_type

 
subroutine Set_defaults(numpar, bunch, MC)
!    type(Matter), intent(inout) :: used_target	! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(Radiation_param), dimension(:), allocatable, intent(inout) :: bunch	! incomming radiation
   type(MC_arrays), dimension(:), allocatable, intent(inout) :: MC	! all MC arrays for particles: photons, electrons and holes; size equals to number of iterations
!    type(Atom), dimension(:), allocatable, intent(inout) :: MD_atoms	! all atoms in MD as objects
   !----------------------------------------
   integer :: i, N_size, j
   N_size = max( m_MC_size, size(bunch), maxval(bunch(:)%NOP) )        ! initial size of all MC particle arrays

   !-------------------------------------------------
   ! All defaults about MC routines:
   
!    numpar%DO_MC = .true.       ! by default, run MC simulation
!    numpar%DO_MD = .false.      ! Activate MD module or not
!    numpar%DO_TTM = .false.     ! Activate TTM module or not
!    numpar%MC_vs_MD = 0         ! MC or MD target model: 0=MC, 1=MD
!    numpar%recalculate_MFPs = .true. ! force program to recalculated MFPs
   
   ! Allocate the array of objects containing the data for all MC iterations:
   if (.not.allocated(MC)) allocate(MC(max(numpar%NMC,1)))
   ! Arrays of particles inside of MC iterations:
   do i = 1, max(numpar%NMC,1)
      if (.not.allocated(MC(i)%MC_Photons)) then
         allocate(MC(i)%MC_Photons(N_size))	! all photons as objects
         do j = 1, N_size
            call set_default_particle(MC(i)%MC_Photons(j))  ! module "Objects"
         enddo
      endif
      
      if (.not.allocated(MC(i)%MC_Electrons)) then
         allocate(MC(i)%MC_Electrons(N_size))	! all electrons as objects
         do j = 1, N_size
            call set_default_particle(MC(i)%MC_Electrons(j))  ! module "Objects"
         enddo
      endif
      
      if (.not.allocated(MC(i)%MC_Positrons)) then
         allocate(MC(i)%MC_Positrons(N_size))	! all positrons as objects
         do j = 1, N_size
            call set_default_particle(MC(i)%MC_Positrons(j))  ! module "Objects"
         enddo
      endif
      
      if (.not.allocated(MC(i)%MC_Holes)) then
         allocate(MC(i)%MC_Holes(N_size))	! all holes as objects
         do j = 1, N_size
            call set_default_particle(MC(i)%MC_Holes(j))  ! module "Objects"
         enddo
      endif
      
      if (.not.allocated(MC(i)%MC_SHIs)) then
         allocate(MC(i)%MC_SHIs(N_size))	! all SHIs as objects
         do j = 1, N_size
            call set_default_particle(MC(i)%MC_SHIs(j))  ! module "Objects"
         enddo
      endif
      
      if (.not.allocated(MC(i)%MC_Atoms_events)) then
         allocate(MC(i)%MC_Atoms_events(N_size))    ! all elastic scattering events as objects
         do j = 1, N_size
            call set_default_particle(MC(i)%MC_Atoms_events(j))  ! module "Objects"
         enddo
      endif

      if (.not.allocated(MC(i)%MC_Muons)) then
         allocate(MC(i)%MC_Muons(N_size))	! all positrons as objects
         do j = 1, N_size
            call set_default_particle(MC(i)%MC_Muons(j))  ! module "Objects"
         enddo
      endif

   enddo
end subroutine Set_defaults

 
 
END MODULE Initial_conditions
