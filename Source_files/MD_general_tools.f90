! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains generally used MD routines
! [1] Thompson et al., J. Chem. Phys. 131, 154107 (2009); https://doi.org/10.1063/1.3245303
! [2] Namilae et al., PRB 76, 144111 (2007); https://doi.org/10.1103/PhysRevB.76.144111
! [3] Mohamed, J. Stat. Phys. 145, 1653 (2011); https://doi.org/10.1007/s10955-011-0364-y
! [4] Fan et al., PRB 92, 094301 (2015); http://dx.doi.org/10.1103/PhysRevB.92.094301

module MD_general_tools
use Universal_constants
use Objects
use Relativity, only: kinetic_energy_from_velosity
use Little_subroutines, only: solve_rate_equation
use Geometries, only: Spherical_R

implicit none

real(8) :: m_kin_press_factor, m_eVA_2_GPa, m_Afs2kg_2_eVA

parameter(m_kin_press_factor = (g_Afs2ms*g_Afs2ms)/g_e)     ! [A^2/fs^2/kg] -> [J] -> [eV]
parameter(m_eVA_2_GPa = g_e*1d30*1d-9)                      ! [eV/A^3] -> [J/m^3] -> [GPa]
parameter(m_Afs2kg_2_eVA = 1.0d10/g_e)                      ! [A/fs^2/kg] -> [eV/A]


 contains


pure function three_body_cos_theta(R_ij, R_ik) result(cos_theta)
   real(8) cos_theta  ! cosine of the angle between atoms j,i,k (centered on i)
   real(8), dimension(3), intent(in) :: R_ij, R_ik
   real(8) rij, rik
   ! Get absolute values:
   rij = Spherical_R(R_ij(1),R_ij(2),R_ij(3))   ! module "Geometries"
   rik = Spherical_R(R_ik(1),R_ik(2),R_ik(3))   ! module "Geometries"
   if ((rij > 0.0d0) .and. (rik > 0.0d0)) then
      cos_theta = DOT_PRODUCT(R_ij, R_ik) / (rij*rik)  ! Eq.(25) in [4]
   else  ! undefined, set default:
      cos_theta = 0.0d0
   endif
end function three_body_cos_theta


pure function d_three_body_cos_theta(R_ij, R_ik, alpha) result(cos_theta)
   real(8) cos_theta  ! derivative of cosine of the angle between atoms j,i,k (centered on i) by R_ij,alpha
   real(8), dimension(3), intent(in) :: R_ij, R_ik
   integer, intent(in) :: alpha  ! companent of derivative (1=x, 2=y, 3=z)
   real(8) rij, rik
   ! Get absolute values:
   rij = Spherical_R(R_ij(1),R_ij(2),R_ij(3))   ! module "Geometries"
   rik = Spherical_R(R_ik(1),R_ik(2),R_ik(3))   ! module "Geometries"
   ! Get cosine:
   cos_theta = three_body_cos_theta(R_ij, R_ik)   ! above
   ! Construct the derivative:
   if ((rij > 0.0d0) .and. (rik > 0.0d0)) then
      cos_theta = 1.0d0/rij*(R_ik(alpha)/rik - R_ij(alpha)/rij*cos_theta) ! Eq.(A4) in [4]
   else  ! undefined, set default:
      cos_theta = 0.0d0
   endif
end function d_three_body_cos_theta



subroutine Thermostat(numpar, MD_atoms, MD_supce, dt)
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   real(8), intent(in) :: dt ! [fs] timestep
   !------------------------------
   real(8) :: therm_border_start(3), therm_border_end(3), Ekin, E_new, Epot, Ta, dE
   integer :: i, k, Nat, Na_thermos
   logical, dimension(size(MD_atoms)) :: atom_in_thermos    ! mask of atoms subject to action of thermostat

   if (numpar%use_thermostat) then    ! only if user included it
      ! Check which atoms in the supercell are subject to the thermostat:
      if ( (MD_supce%thermo_dx_start+MD_supce%thermo_dx_end >= MD_supce%vect(1)) .or. &
           (MD_supce%thermo_dy_start+MD_supce%thermo_dy_end >= MD_supce%vect(2)) .or. &
           (MD_supce%thermo_dz_start+MD_supce%thermo_dz_end >= MD_supce%vect(3)) ) then ! all atoms
         atom_in_thermos = .true.   ! all atoms are in the thermostat
      else  ! not all atoms are in the thermostat
         Nat = size(MD_atoms)   ! number of atoms
         ! Define borders of the thermostat:
         therm_border_start(1) = MD_supce%x_start + MD_supce%thermo_dx_start  ! [A]
         therm_border_start(2) = MD_supce%y_start + MD_supce%thermo_dy_start  ! [A]
         therm_border_start(3) = MD_supce%z_start + MD_supce%thermo_dz_start  ! [A]
         therm_border_end(1) = MD_supce%x_end - MD_supce%thermo_dx_end  ! [A]
         therm_border_end(2) = MD_supce%y_end - MD_supce%thermo_dy_end  ! [A]
         therm_border_end(3) = MD_supce%z_end - MD_supce%thermo_dz_end  ! [A]
         ! Check for each atom if it is inside the thermostat:
         do i = 1, Nat
            if ( ANY(MD_atoms(i)%R(:) < therm_border_start(:)) .or. &
                 ANY(MD_atoms(i)%R(:) > therm_border_end(:)) ) then ! atom in the thermostat
               atom_in_thermos(i) = .true.
            else    ! atom is not in the thermostat
               atom_in_thermos(i) = .false.
            endif
         enddo  ! i
      endif

      ! Apply thermostat on the selected atoms:
      select case (numpar%thermostat_inx)   ! which thermostat to use
      case (0)  ! Berendsen
         ! Total number of atoms inside of the thermostat:
         Na_thermos = COUNT(atom_in_thermos)
         ! Find total energy of the subset of atoms:
         call get_total_energy(MD_atoms, Ekin, Epot, atom_in_thermos) ! below
         ! Get temperature (Ta) of these atoms:
         call get_temperature_from_energy(dble(Na_thermos), Ekin, Ta, .true., .false.)    ! below
         ! Get change of temperature due to interaction with the bath (thermostat):
         call solve_rate_equation(Ta, MD_supce%Bath_T, dt, MD_supce%therm_tau, 2) ! module "Little_subroutines"
         ! Get energy corresonding to the updated temperature:
         call get_energy_from_temperature(dble(Na_thermos), Ta, E_new, .true., .false.)  ! below
         ! Get change of the total energy:
         dE = E_new - Ekin  ! [eV]
         ! Update atomic velosities correspondingly, to deliver this energy to atoms:
         call rescale_velosities(MD_atoms, Ekin, dE, atom_in_thermos) ! below

      end select
   endif ! (numpar%use_thermostat)
end subroutine Thermostat


subroutine rescale_velosities(MD_atoms, Ekin, dE, atoms_mask)
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   real(8), intent(in) :: Ekin  ! [eV] total kinetic energy of all considerred atoms
   real(8), intent(in) :: dE    ! [eV] change of the kinetic energy of atoms
   logical, dimension(:), intent(in), optional :: atoms_mask  ! mask for choosing which atoms to include
   !-------------------
   real(8) :: eps, b, alpha, Enew, Epot
   integer :: i, Nat

   eps = 1.0d-8   ! acceptable change of energy

   Nat = size(MD_atoms)

   if ((Ekin > eps) .and. (ABS(dE) > (Ekin*eps)) ) then ! makes sense to change it
      ! Get the scaling coefficient:
      b = dE/Ekin   ! ratio entering the scaling factor
      if (b < -1.0d0) then
         alpha = 0.0d0  ! not possible to remove more energy than present in the system
      else
         alpha = dsqrt(1.0d0 + b)    ! scaling factor for velosities, providing correct energy change
      endif
      ! Rescale atomic velocities:
      if (present(atoms_mask)) then ! subset of atoms
         !$omp parallel private (i)
         !$omp do schedule(dynamic)
         do i = 1, Nat
            if (atoms_mask(i)) then    ! atom i belongs to the subset
               MD_atoms(i)%V(:) = MD_atoms(i)%V(:) * alpha
               MD_atoms(i)%V0(:) = MD_atoms(i)%V0(:) * alpha
            endif
         enddo
         !$omp enddo
         !$omp end parallel
      else  ! all atoms
         !$omp parallel private (i)
         !$omp do schedule(dynamic)
         do i = 1, Nat
            MD_atoms(i)%V(:) = MD_atoms(i)%V(:) * alpha
            MD_atoms(i)%V0(:) = MD_atoms(i)%V0(:) * alpha
         enddo
         !$omp enddo
         !$omp end parallel
      endif

      ! Check that everything is correct:
      if (present(atoms_mask)) then ! subset of atoms
         call get_total_energy(MD_atoms, Enew, Epot, atoms_mask)   ! below
      else  ! all atoms
         call get_total_energy(MD_atoms, Enew, Epot)   ! below
      endif
      if ( abs(Enew - (Ekin + dE)) > 1.0d-4*(Ekin + dE)) then ! energy is not concerved due to rescaling
         write(*,'(a,f,f,f,f)') 'Possible error in rescale_velosities: ', dE, Enew-(Ekin+dE), Enew, (Ekin+dE)
      endif
   endif
end subroutine rescale_velosities


subroutine get_total_energy(MD_atoms, Ekin, Epot, atoms_mask)
   type(Atom), dimension(:), intent(in) :: MD_atoms	! all atoms in MD as objects
   real(8), intent(out) :: Ekin, Epot   ! [eV] total kinetic and potential energy of atoms
   logical, dimension(:), intent(in), optional :: atoms_mask  ! mask for chising which atoms to include
   !--------------
   real(8) :: Vtot, Ekin_atom
   integer :: i, Nat

   Nat = size(MD_atoms)  ! total number of atoms
   Ekin = 0.0d0 ! to start with

   if (.not.present(atoms_mask)) then   ! all atoms:
      call get_total_energies(MD_atoms, Ekin, Epot)    ! below
      Ekin = Ekin * dble(Nat)   ! average energy -> total energy
      Epot = Epot * dble(Nat)   ! average energy -> total energy
   else ! subset of atoms chosen by the mask:
      !$omp parallel private (i, Vtot, Ekin_atom)
      !$omp do reduction(+:Ekin) schedule(dynamic)
      do i = 1, Nat    ! for all atoms
         if (atoms_mask(i)) then    ! chose from subset
            Vtot = dsqrt(SUM( MD_atoms(i)%V(:)*MD_atoms(i)%V(:) )) * 1.0d5       ! [A/fs] -> [m/s]
            Ekin_atom = kinetic_energy_from_velosity(Vtot, MD_atoms(i)%Mass)    ! module "Relativity"
            ! Total kinetic energy:
            Ekin = Ekin + Ekin_atom   ! [eV]
         endif
      enddo
      !$omp enddo
      !$omp end parallel
      Epot = 0.5d0 * SUM(MD_atoms(:)%U, MASK=atoms_mask)  ! [eV] total potential energy: 1/2*sum_i,j(U)
   endif
end subroutine get_total_energy



subroutine Damp_pressure(numpar, MD_atoms, MD_supce)    ! following methodology from Ref.[2]
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   !------------------------------
   integer :: Nat, i
   real(8) :: visc_border_start(3), visc_border_end(3)  ! boundaries of viscous regions around the borders
   real(8) :: Force(3), viscos_factor

   if (numpar%damp_pressure) then    ! only if user included it
      Nat = size(MD_atoms)
      ! Define borders of the thermostat:
      visc_border_start(1) = MD_supce%x_start + MD_supce%pressdamp_dx_start  ! [A]
      visc_border_start(2) = MD_supce%y_start + MD_supce%pressdamp_dy_start  ! [A]
      visc_border_start(3) = MD_supce%z_start + MD_supce%pressdamp_dz_start  ! [A]
      visc_border_end(1) = MD_supce%x_end - MD_supce%pressdamp_dx_end  ! [A]
      visc_border_end(2) = MD_supce%y_end - MD_supce%pressdamp_dy_end  ! [A]
      visc_border_end(3) = MD_supce%z_end - MD_supce%pressdamp_dz_end  ! [A]

      !$omp parallel private (i, Force, viscos_factor)
      !$omp do schedule(dynamic)
      do i = 1, Nat
         Force = 0.0d0 ! to start with
         ! Get the viscose forces at each border:
         if ( MD_atoms(i)%R(1) < visc_border_start(1) ) then
            viscos_factor = (abs(MD_atoms(i)%R(1)-visc_border_start(1))/MD_supce%pressdamp_dx_start)**3
            viscos_factor = viscos_factor / MD_supce%press_tau  ! [1/fs]
            Force(:) = Force(:) + MD_atoms(i)%V(:) * viscos_factor  ! [A/fs^2]
         endif
         if ( MD_atoms(i)%R(2) < visc_border_start(2) ) then
            viscos_factor = (abs(MD_atoms(i)%R(2)-visc_border_start(2))/MD_supce%pressdamp_dy_start)**3
            viscos_factor = viscos_factor / MD_supce%press_tau  ! [1/fs]
            Force(:) = Force(:) + MD_atoms(i)%V(:) * viscos_factor  ! [A/fs^2]
         endif
         if ( MD_atoms(i)%R(3) < visc_border_start(3) ) then
            viscos_factor = (abs(MD_atoms(i)%R(3)-visc_border_start(3))/MD_supce%pressdamp_dz_start)**3
            viscos_factor = viscos_factor / MD_supce%press_tau  ! [1/fs]
            Force(:) = Force(:) + MD_atoms(i)%V(:) * viscos_factor  ! [A/fs^2]
         endif
         if ( MD_atoms(i)%R(1) > visc_border_end(1) ) then
            viscos_factor = (abs(MD_atoms(i)%R(1)-visc_border_end(1))/MD_supce%pressdamp_dx_end)**3
            viscos_factor = viscos_factor / MD_supce%press_tau  ! [1/fs]
            Force(:) = Force(:) + MD_atoms(i)%V(:) * viscos_factor  ! [A/fs^2]
         endif
         if ( MD_atoms(i)%R(2) > visc_border_end(2) ) then
            viscos_factor = (abs(MD_atoms(i)%R(2)-visc_border_end(2))/MD_supce%pressdamp_dy_end)**3
            viscos_factor = viscos_factor / MD_supce%press_tau  ! [1/fs]
            Force(:) = Force(:) + MD_atoms(i)%V(:) * viscos_factor  ! [A/fs^2]
         endif
         if ( MD_atoms(i)%R(3) > visc_border_end(3) ) then
            viscos_factor = (abs(MD_atoms(i)%R(3)-visc_border_end(3))/MD_supce%pressdamp_dz_end)**3
            viscos_factor = viscos_factor / MD_supce%press_tau  ! [1/fs]
            Force(:) = Force(:) + MD_atoms(i)%V(:) * viscos_factor  ! [A/fs^2]
         endif
         ! Update forces acting on the atom:
         MD_atoms(i)%Force(:) = MD_atoms(i)%Force(:) + Force(:) / MD_atoms(i)%Mass * m_Afs2kg_2_eVA ! [eV/A]
      enddo
      !$omp enddo
      !$omp end parallel
   endif
end subroutine Damp_pressure



! Check periodic boundary crossing:
subroutine Check_MD_periodic_boundaries(MD_atoms, MD_supce)
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   !---------------------------
   integer :: x_shift, y_shift, z_shift
   integer :: Nat, i
   Nat = size(MD_atoms) ! number of atoms
   
   if ( ANY(MD_supce%boundary(:) == 1) ) then    ! there are some boundaries that are periodic:
      
      !$omp parallel private (i, x_shift, y_shift, z_shift)
      !$omp do
      do i = 1, Nat    ! check all atoms for crossing a periodic boundary
         
         if (MD_supce%boundary(1) == 1) then ! periodic along X
            if (MD_atoms(i)%R(1) < MD_supce%x_start) then ! cross "left" boundary, place atom back
               x_shift = CEILING ( ABS(MD_atoms(i)%R(1)-MD_supce%x_start)/MD_supce%vect(1) )
               MD_atoms(i)%R(1) = MD_atoms(i)%R(1) + MD_supce%vect(1) * dble(x_shift)
               MD_atoms(i)%R0(1) = MD_atoms(i)%R0(1) + MD_supce%vect(1) * dble(x_shift)
            endif
            if (MD_atoms(i)%R(1) > MD_supce%x_end) then ! cross "right" boundary, place atom back
               x_shift = CEILING ( ABS(MD_atoms(i)%R(1)-MD_supce%x_end)/MD_supce%vect(1) )
               MD_atoms(i)%R(1) = MD_atoms(i)%R(1) - MD_supce%vect(1) * dble(x_shift)
               MD_atoms(i)%R0(1) = MD_atoms(i)%R0(1) - MD_supce%vect(1) * dble(x_shift)
            endif
         endif
         
         if (MD_supce%boundary(2) == 1) then ! periodic along Y
            if (MD_atoms(i)%R(2) < MD_supce%y_start) then ! cross "left" boundary, place atom back
               y_shift = CEILING ( ABS(MD_atoms(i)%R(2)-MD_supce%y_start)/MD_supce%vect(2) )
               MD_atoms(i)%R(2) = MD_atoms(i)%R(2) + MD_supce%vect(2) * dble(y_shift)
               MD_atoms(i)%R0(2) = MD_atoms(i)%R0(2) + MD_supce%vect(2) * dble(y_shift)
            endif
            if (MD_atoms(i)%R(2) > MD_supce%y_end) then ! cross "right" boundary, place atom back
               y_shift = CEILING ( ABS(MD_atoms(i)%R(2)-MD_supce%y_end)/MD_supce%vect(2) )
               MD_atoms(i)%R(2) = MD_atoms(i)%R(2) - MD_supce%vect(2) * dble(y_shift)
               MD_atoms(i)%R0(2) = MD_atoms(i)%R0(2) - MD_supce%vect(2) * dble(y_shift)
            endif
         endif
         
         if (MD_supce%boundary(3) == 1) then ! periodic along Z
            if (MD_atoms(i)%R(3) < MD_supce%z_start) then ! cross "left" boundary, place atom back
               z_shift = CEILING ( ABS(MD_atoms(i)%R(3)-MD_supce%z_start)/MD_supce%vect(3) )
               MD_atoms(i)%R(3) = MD_atoms(i)%R(3) + MD_supce%vect(3) * dble(z_shift)
               MD_atoms(i)%R0(3) = MD_atoms(i)%R0(3) + MD_supce%vect(3) * dble(z_shift)
            endif
            if (MD_atoms(i)%R(3) > MD_supce%z_end) then ! cross "right" boundary, place atom back
               z_shift = CEILING ( ABS(MD_atoms(i)%R(3)-MD_supce%z_end)/MD_supce%vect(3) )
               MD_atoms(i)%R(3) = MD_atoms(i)%R(3) - MD_supce%vect(3) * dble(z_shift)
               MD_atoms(i)%R0(3) = MD_atoms(i)%R0(3) - MD_supce%vect(3) * dble(z_shift)
            endif
         endif
      
      enddo
      !$omp enddo
      !$omp end parallel
         
   endif
end subroutine Check_MD_periodic_boundaries


 

! Sample velosity according to Maxwell distribution:
subroutine Sample_Maxwelian(T, Mass, Vx, vy, Vz)
   real(8), intent(in) :: T     ! [eV] Temperature to set the velocities accordingly
   real(8), intent(in) :: Mass  ! [kg] mass of the atom
   real(8), intent(out) :: Vx, Vy, Vz   ! velocities [A/fs]
   real(8) RN(6) ! random numbers
   real(8) E, V, theta, phi, cos_phi, cos_rn, cos2
   integer i

   E = T*25.0d0  ! just to start
   do while (E > T*15.0d0) ! exclude too high energies
      do i = 1,size(RN)
         call random_number(RN(i))
      enddo
      cos_rn = cos(RN(3)*g_half_Pi)
      cos2 = cos_rn*cos_rn
      E = T*(-log(RN(1)) - log(RN(2))*cos2) ! [eV] Using Eq.(9) Ref.[3]
      E = 2.0d0*E ! half of the energy is kinetic, half potential, so double it to get the right temperature
   enddo
   V = sqrt(E*2.0d0*g_e/Mass)*1d-5 ! [A/fs] absolute value of velocity
   theta = 2.0d0*g_Pi*RN(4) ! angle
   phi = -g_half_Pi + g_Pi*RN(5)  ! second angle
   cos_phi = cos(phi)
   ! To exclude artifacts due to choice of coordinate system, sample them randomly:
   if (RN(6) <= 0.33) then  ! Vz along Z
      Vx = V*cos_phi*cos(theta)
      Vy = V*cos_phi*sin(theta)
      Vz = V*sin(phi)
   elseif(RN(6) <=0.67) then    ! Vy along Z
      Vz = V*cos_phi*cos(theta)
      Vx = V*cos_phi*sin(theta)
      Vy = V*sin(phi)
   else ! Vx along Z
      Vy = V*cos_phi*cos(theta)
      Vz = V*cos_phi*sin(theta)
      Vx = V*sin(phi)
   endif
end subroutine Sample_Maxwelian
 


subroutine remove_CoM_velocity(atoms, print_out)
   type(Atom), dimension(:), intent(inout) :: atoms ! all atoms in MD as objects
   logical, optional, intent(in) :: print_out ! print our c-o-m momemntum
   real(8) :: vx, vy, vz, Masstot
   integer :: i
   
   Masstot = SUM(atoms(:)%Mass) ! net-mass of all atoms
   
   ! Center of mass velocities:
   vx = SUM(atoms(:)%Mass * atoms(:)%V(1))/Masstot
   vy = SUM(atoms(:)%Mass * atoms(:)%V(2))/Masstot
   vz = SUM(atoms(:)%Mass * atoms(:)%V(3))/Masstot

   ! Subtract the velocity of the center of mass:
   !$omp parallel private (i)
   !$omp do schedule(dynamic)
   do i = 1, size(atoms)     ! for all atoms
      atoms(i)%V(1) = atoms(i)%V(1) - vx
      atoms(i)%V(2) = atoms(i)%V(2) - vy
      atoms(i)%V(3) = atoms(i)%V(3) - vz
   enddo
   !$omp enddo
   !$omp end parallel

   ! In case user wants to see CoM velosities:
   if (present(print_out)) then
      if (print_out) write(*,'(a,es25.16,es25.16,es25.16)') 'CoM:', vx, vy, vz
   endif
end subroutine remove_CoM_velocity


subroutine Quenching(numpar, MD_atoms, t_time)
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   real(8), intent(in) :: t_time ! [fs] current timestep
   !------------------------------
   integer k, Nat
   if (numpar%do_quenching) then    ! only if user included it
      if ((t_time >= numpar%t_quench_start) .and. (numpar%t_quench_run > numpar%dt_quench)) then ! it's time to cool atoms down:
         numpar%t_quench_run = 0.0d0    ! reset counter
         Nat = size(MD_atoms)   ! total number of atoms
         ! Cooling atoms instantly to Ta=0:
         !$omp parallel private (k)
         !$omp do
         do k = 1, Nat ! for all atoms, instantaneous cool-down
            MD_atoms(k)%V(:) = 0.0d0
            MD_atoms(k)%V0(:) = 0.0d0
         enddo
         !$omp enddo
         !$omp end parallel
         
         !call Atomic_kinetic_energies()
      else ! it's not the time to quench, keep counting:
         numpar%t_quench_run = numpar%t_quench_run + numpar%dt_MD ! [fs]
      endif
   endif
end subroutine Quenching


subroutine rescale_supercell(MD_atoms_in, MD_atoms, MD_supce_in, MD_supce, rescal)
   type(Atom), dimension(:), intent(in) :: MD_atoms_in
   type(Atom), dimension(:), intent(inout) :: MD_atoms
   type(MD_supcell), intent(in) :: MD_supce_in
   type(MD_supcell), intent(inout) :: MD_supce
   real(8), intent(in) :: rescal    ! rescaling coefficient
   integer :: i
   ! Rescale supercell size:
   MD_supce%x_start = MD_supce_in%x_start * rescal
   MD_supce%x_end = MD_supce_in%x_end * rescal
   MD_supce%y_start = MD_supce_in%y_start * rescal
   MD_supce%y_end = MD_supce_in%y_end * rescal
   MD_supce%z_start = MD_supce_in%z_start * rescal
   MD_supce%z_end = MD_supce_in%z_end * rescal
   MD_supce%vect(:) = MD_supce_in%vect(:) * rescal
   MD_supce%V = PRODUCT(MD_supce%vect)  ! [A^3]
   ! Rescale atomic coordinates:
   do i = 1, size(MD_atoms)
      MD_atoms(i)%R(:) = MD_atoms_in(i)%R(:) * rescal
   enddo
end subroutine rescale_supercell



subroutine get_nearest_neighbors_list(MD_atoms, MD_supce, MD_pots, numpar)
   ! Define parameters of nearest neighbors (NN) for all atoms
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   !---------------------
   real(8) :: r_cut, r(1), d(1), a_r, x, y, z, r_max
   integer :: Nat, i, j, k, Npot
   integer, dimension(size(MD_atoms)) :: coun
   
   Nat = size(MD_atoms) ! total number of atoms in the simulation box
   Npot = size(MD_pots,1)   ! number of different elements
   coun = 0 ! to start with
   
   ! No interaction between atoms at distance larger than this:
   r_cut = 0.0d0   ! to start with
   do i = 1, Npot
      do j = 1, Npot
         do k = 1, size(MD_pots(i,j)%Set)   ! for all potentials defined
            ASSOCIATE (ARRAY => MD_pots(i,j)%Set(k)%Par) ! this is the sintax we have to use to check the class of defined types
               r_max = ARRAY%d_cut + 10.0d0*ARRAY%dd   ! [A] cut off distance
               if (r_cut < r_max) then    ! maximal length
                  r_cut = r_max   ! this is the maximal length a potential can reach
               endif
            ENDASSOCIATE
         enddo ! k
      enddo ! j
   enddo ! i
   
   ! Make sure NN lists are allocated:
   if (.not.allocated(numpar%Neighbors_Num)) then   ! number of nearest neighbors for each atom
      allocate(numpar%Neighbors_Num(Nat), source = 0)
      allocate(numpar%Neighbors_list(Nat,Nat), source = 0)
      allocate(numpar%Neighbors_Rij(Nat,Nat), source = 1.0d10)
      allocate(numpar%Neighbors_Xij(Nat,Nat,3), source = 1.0d10)
   endif
   
   ! Define the parameters of nearest neighbors:
   !$omp parallel private (i, j, a_r, x, y, z)
   !$omp do schedule(dynamic)
   do i = 1, Nat    ! for all atoms
      ! to reset counter:
      coun(i) = 0
      numpar%Neighbors_list(i,:) = 0
      numpar%Neighbors_Rij(i,:) = 1.0d10
      numpar%Neighbors_Xij(i,:,:) = 1.0d10
      do j = 1, Nat    ! for all pair-atoms
         if (i /= j) then   ! only neighbors, no self-interaction
            ! Get the interaction distance, accounting for periodic boundaries:
            call shortest_distance(i, j, MD_atoms, MD_supce, a_r, x1=x, y1=y, z1=z) ! below
            
            if (a_r < 1.0d-6) then
               print*, 'Error in get_nearest_neighbors_list:'
               print*, 'a_r=', a_r
               print*, x, y, z
               print*, 'Atoms:', i, j
               print*, MD_atoms(i)%R(:)
               print*, MD_atoms(j)%R(:)
            endif
            
            ! If atoms are interacting, save their parameters:
            if (a_r < r_cut) then   ! these atoms do interact
               coun(i) = coun(i) + 1
               numpar%Neighbors_list(i, coun(i)) = j    ! atom i interacts with this atom
               numpar%Neighbors_Xij(i, coun(i), 1) = x  ! at this distance, X
               numpar%Neighbors_Xij(i, coun(i), 2) = y  ! at this distance, Y
               numpar%Neighbors_Xij(i, coun(i), 3) = z  ! at this distance, Z
               numpar%Neighbors_Rij(i, coun(i)) = a_r   ! at this distance, R
            endif ! (a_r < r_cut)
         endif ! (i /= j)
      enddo ! j = 1, Nat
      ! Save the parameters of NN:
      numpar%Neighbors_Num(i) = coun(i) ! that's how many nearest neighbours there are for the atom i
   enddo ! i = 1, Nat
   !$omp enddo
   !$omp end parallel
end subroutine get_nearest_neighbors_list


pure subroutine shortest_distance(i, j, MD_atoms, MD_supce, r_shortest, x1, y1, z1, disp, ind_disp)
   integer, intent(in) :: i, j  ! indices of atoms interacting
   type(Atom), dimension(:), intent(in) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   real(8), intent(out) :: r_shortest   ! [A] shortest distance between two atoms, accounting for periodic boundaries
   real(8), intent(out), optional :: x1, y1, z1 ! [A] distance along X, Y, Z
   real(8), intent(in), optional :: disp  ! [A] displacement of atom i
   integer, intent(in), optional :: ind_disp ! along which axis is the displacement
   !---------------------
   integer :: ix, iy, iz
   real(8) :: r_min, x_min, y_min, z_min, x, y, z, r, x_min_0, y_min_0, z_min_0
   
   ! To start with:
   x_min = MD_atoms(i)%R(1) - MD_atoms(j)%R(1)
   y_min = MD_atoms(i)%R(2) - MD_atoms(j)%R(2)
   z_min = MD_atoms(i)%R(3) - MD_atoms(j)%R(3)
   if (present(disp) .and. present(ind_disp)) then ! atom is displaeced:
      select case (ind_disp)
      case (1)
         x_min = x_min + disp
      case (2)
         y_min = y_min + disp
      case (3)
         z_min = z_min + disp
      end select
   endif
   r_min = sqrt(x_min*x_min + y_min*y_min + z_min*z_min)
   ! Save to reuse below:
   x_min_0 = x_min
   y_min_0 = y_min
   z_min_0 = z_min

   ! And only if the boundaries are periodic, check if atoms in the image cells are closer:
   if (any((MD_supce%boundary(:) == 1))) then
      ! To start with:
      x = x_min
      y = y_min
      z = z_min
      r = r_min
      ! Check all image cells around the supercell:
      do ix = -1,1
         do iy = -1,1
            do iz = -1,1
               if (MD_supce%boundary(1) == 1) then ! periodic along X
                  !x = MD_atoms(i)%R(1) - MD_atoms(j)%R(1) + dble(ix)*MD_supce%vect(1)
                  x = x_min_0 + dble(ix)*MD_supce%vect(1)
               endif
               if (MD_supce%boundary(2) == 1) then ! periodic along Y
                  !y = MD_atoms(i)%R(2) - MD_atoms(j)%R(2) + dble(iy)*MD_supce%vect(2)
                  y = y_min_0 + dble(iy)*MD_supce%vect(2)
               endif
               if (MD_supce%boundary(3) == 1) then ! periodic along Z
                  !z = MD_atoms(i)%R(3) - MD_atoms(j)%R(3) + dble(iz)*MD_supce%vect(3)
                  z = z_min_0 + dble(iz)*MD_supce%vect(3)
               endif
               r = sqrt(x*x + y*y + z*z)
               ! If the distance is shorter to this image copy of the atom, save it:
               if (r < r_min) then
                  r_min = r
                  x_min = x
                  y_min = y
                  z_min = z
               endif
            enddo ! iz = -1,1
         enddo ! iy = -1,1
      enddo ! ix = -1,1
   endif ! (any((MD_supce%boundary(:) == 1)))
   ! Out:
   r_shortest = r_min
   if (present(x1)) x1 = x_min
   if (present(y1)) y1 = y_min
   if (present(z1)) z1 = z_min
end subroutine shortest_distance



pure subroutine Verlet_step(R, R0, V, V0, A, dt, dt_half, first_step)
   real(8), dimension(3), intent(inout) :: R    ! [A] coordinates
   real(8), dimension(3), intent(inout) :: R0   ! [A] coordinates on last timestep
   real(8), dimension(3), intent(inout) :: V    ! [A/fs] velosities
   real(8), dimension(3), intent(inout) :: V0   ! [A/fs] velosities on last timestep
   real(8), dimension(3), intent(in) :: A   ! [A/fs^2] acceleration
   real(8), intent(in) :: dt                ! [fs] timestep
   real(8), intent(in) :: dt_half           ! dt/2, for speeding up calculations
   logical, intent(in) :: first_step        ! first or second step of Verlet

   ! update velocities:
   V(:) = V0(:) + A(:)*dt_half
   ! On the second step, update velocities for the next half-timestep:
   V0(:) = V(:)
   
   if (first_step) then    ! update coordinates too:
      R(:) = R0(:) + dt*V(:) ! new coordinates
      ! Update for the next timestep:
      R0(:) = R(:)
   endif
end subroutine Verlet_step




pure subroutine Verlet_step_OLD(R, R0, V, V0, A, dt, dt2, dt_half, do_coords)
   real(8), dimension(3), intent(inout) :: R    ! [A] coordinates
   real(8), dimension(3), intent(inout) :: R0   ! [A] coordinates on last timestep
   real(8), dimension(3), intent(inout) :: V    ! [A/fs] velosities
   real(8), dimension(3), intent(inout) :: V0   ! [A/fs] velosities on last timestep
   real(8), dimension(3), intent(in) :: A   ! [A/fs^2] acceleration
   real(8), intent(in) :: dt                ! [fs] timestep
   real(8), intent(in) :: dt2, dt_half      ! dt^2/2, dt/2, for speeding up calculations
   logical, intent(in) :: do_coords         ! change coordinates or not (first or second step of Verlet)

   ! Verlet velocities (half-step):
   V(:) = V0(:) + A(:)*dt_half

   ! Verlet step of coordinates (full step):
   if (do_coords) then
      R(:) = R0(:) + dt*V0(:) + A(:)*dt2 ! new coordinates
      ! Update for the next timestep:
      R0(:) = R(:)
      ! On the second step, update velocities for the next half-timestep:
      V0(:) = V(:)
   endif
end subroutine Verlet_step_OLD



subroutine find_which_potential(MD_atoms, atom_1, atom_2, MD_pots, KOP1, KOP2)
   type(Atom), dimension(:), intent(in) :: MD_atoms	! all atoms in MD as objects
   integer, intent(in) :: atom_1, atom_2    ! indices of the atoms in the MD array
   type(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   integer, intent(out) :: KOP1, KOP2    ! indices of the potentials in the Pot array
   !-------------------------------
   integer :: i, j
   logical :: found_pot
   
   found_pot = .false. ! to start with
   
   i = 0    ! to start
   j = 0    ! to start
   FP:do i = 1, size(MD_pots,1)
      do j = 1, size(MD_pots,2)
!          print*, 'P:', trim(adjustl(MD_pots(i,j)%El1))//' '//trim(adjustl(MD_pots(i,j)%El2))
         if ( (trim(adjustl(MD_atoms(atom_1)%Name)) == trim(adjustl(MD_pots(i,j)%El1)) ) .and. &
              (trim(adjustl(MD_atoms(atom_2)%Name)) == trim(adjustl(MD_pots(i,j)%El2)) ) ) then
            ! save indices of the potential:
            KOP1 = i
            KOP2 = j
            found_pot = .true.   ! mark the potential as found
            exit FP ! potential indices found, stop the subroutine
         endif
      enddo ! j
   enddo FP

   if (.not.found_pot) then   ! there is no potential given for this pair of atoms:
      print*, 'ERROR in find_which_potential: no potential found for atoms '// &
               trim(adjustl(MD_atoms(atom_1)%Name))//'-'// &
               trim(adjustl(MD_atoms(atom_2)%Name))
   endif
!    pause 'Pause find_which_potential'
end subroutine find_which_potential



subroutine get_pressure(MD_atoms, MD_supce, Press_tens, done_Press_Pot) ! Following Eq.(28) in Ref.[1]
   type(Atom), dimension(:), intent(in), target :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   real(8), dimension(3,3), intent(out) :: Press_tens   ! pressure tensor
   logical, intent(in) :: done_Press_Pot   ! if potential contribution is precalculated, use it
   !----------------------
   real(8), dimension(3,3) :: Press_Kin   ! Kinetic part of the pressure tensor
   real(8), dimension(3,3) :: Press_Pot   ! Potential part of the pressure tensor
   integer :: i, j, i_at, Nat
   real(8), pointer :: Mass
   
   Nat = size(MD_atoms) ! total number of atoms
   
   Press_Kin(:,:) = 0.0d0   ! to start with
   Press_Pot(:,:) = 0.0d0   ! to start with
   !$omp parallel private (i_at, Mass, i, j)
   !$omp do 
   do i_at = 1, Nat
      Mass => MD_atoms(i_at)%Mass
      do i = 1, 3
         do j = 1, 3
            ! Kinetic part:
            Press_Kin(i,j) = Press_Kin(i,j) + MD_atoms(i_at)%V(i)*MD_atoms(i_at)%V(j)*Mass*m_kin_press_factor   ! [eV]
         enddo   
      enddo
   enddo
   !$omp enddo
   ! Potential part:
   if (.not.done_Press_Pot) then
      !$omp do
      do i_at = 1, Nat
         do i = 1, 3
            do j = 1, 3
            ! Potential part:
            Press_Pot(i,j) = Press_Pot(i,j) + MD_atoms(i_at)%r(i)*MD_atoms(i_at)%Force(j)   ! [eV]
            enddo   
         enddo
      enddo
      !$omp enddo
   else ! use precalculated potential part:
      do i = 1, 3
         do j = 1, 3
            Press_Pot(i,j) = MD_supce%P_pot(i,j)    ! [eV]
         enddo
      enddo
   endif
   !$omp end parallel
   nullify(Mass)
   
   ! Total pressure tensor:
   Press_tens = (Press_Kin + 0.5d0*Press_Pot) / MD_supce%V * m_eVA_2_GPa  ! [GPa]
!     print*, 'Ptot:', Press_tens
!     print*, 'Pkin:', Press_Kin
!     print*, 'Ppot:', Press_Pot
!     print*, 'V=', MD_supce%V
end subroutine get_pressure



subroutine get_total_energies(MD_atoms, Ekin, Epot)
   type(Atom), dimension(:), intent(in) :: MD_atoms	! all atoms in MD as objects
   real(8), intent(out) :: Ekin, Epot   ! [eV/atom] average kinetic and potential energy of atoms
   !-----------------------
   integer :: Nat, i
   real(8) :: Ekin_atom, Vtot
   
   Nat = size(MD_atoms) ! total number of atoms
   Ekin = 0.0d0 ! to start with
   
   !$omp parallel private (i, Vtot, Ekin_atom)
   !$omp do reduction(+:Ekin)
   do i = 1, Nat    ! for all atoms
      Vtot = dsqrt(SUM( MD_atoms(i)%V(:)*MD_atoms(i)%V(:) )) * g_Afs2ms ! [A/fs] -> [m/s]
      Ekin_atom = kinetic_energy_from_velosity(Vtot, MD_atoms(i)%Mass)  ! module "Relativity"
      ! Total kinetic energy:
      Ekin = Ekin + Ekin_atom   ! [eV]
   enddo
   !$omp enddo
   !$omp end parallel
   Ekin = Ekin/dble(Nat)   ! [eV/atom] kinetic energy
   Epot = 0.5d0 * SUM(MD_atoms(:)%U)/dble(Nat)  ! [eV/atom] total potential energy: 1/2*sum_i,j(U)
end subroutine get_total_energies



pure subroutine get_energy_from_temperature(Na, Ta, Ekin, in_K, periodic)
   real(8), intent(in) :: Na ! number of atoms
   real(8), intent(in) :: Ta ! temperature of atoms [K] or [eV]
   real(8), intent(out) :: Ekin ! total kinetic energy of atoms [eV]
   logical, intent(in) :: in_K  ! is temperature is [K]? (if not, it must be in [eV])
   logical, intent(in) :: periodic  ! periodic boundary conditions used or not?
   real(8) :: T_convet, DoF
   if (in_K) then   ! temperature is in [K]
      T_convet = g_kb_EV
   else ! temperature is in [eV]
      T_convet = 1.0d0  ! no conversion factor is needed
   endif

   if (periodic) then
      DoF = 6.0d0   ! those degrees of freedom need to be subtracted
   else
      DoF = 0.0d0   ! all degrees of freedom are unconstrained
   endif

   Ekin = 0.5d0*(3.0d0*Na - DoF)*(Ta*T_convet) ! kinetic energy in a box with periodic boundary
end subroutine get_energy_from_temperature



pure subroutine get_temperature_from_energy(Na, Ekin, Ta, in_K, periodic)
   real(8), intent(in) :: Na    ! number of atoms
   real(8), intent(in) :: Ekin  ! total kinetic energy of atoms [eV]
   real(8), intent(out) :: Ta   ! temperature of atoms [K] or [eV]
   logical, intent(in) :: in_K  ! temperature output is [K]? (if not, it will be in [eV])
   logical, intent(in) :: periodic  ! periodic boundary conditions used or not?
   real(8) :: T_convet, DoF

   if (in_K) then   ! temperature required in [K]
      T_convet = g_kb
   else ! temperature required in [eV]
      T_convet = 1.0d0  ! no conversion factor is needed
   endif

   if (periodic) then
      DoF = 6.0d0   ! those degrees of freedom need to be subtracted
   else
      DoF = 0.0d0   ! all degrees of freedom are unconstrained
   endif

   Ta = 2.0d0/(3.0d0*Na - DoF)*Ekin*T_convet ! temperature in a box with periodic boundary
end subroutine get_temperature_from_energy



pure function Fermi_cut_off(a_r, r_L, d_L) result(f_cut) ! cut-off function (Fermi shape)
   real(8) :: f_cut
   real(8), intent(in) :: a_r, r_L, d_L
   real(8) :: eps
   eps = 1.0d-6
   if (d_L <= eps) then ! zero smearing, use step function
      if (a_r <= r_L) then
         f_cut = 1.0d0
      else
         f_cut = 0.0d0
      endif
   else ! calculte it
      f_cut = 1.0d0/(1.0d0 + dexp((a_r - r_L)/d_L))
   endif
end function Fermi_cut_off


pure function d_Fermi_cut_off(a_r, r_L, d_L) result(d_f_cut) ! derivative of cut-off function (Fermi shape)
   real(8) :: d_f_cut
   real(8), intent(in) :: a_r, r_L, d_L
   real(8) :: exp_r, exp_r2, arg
   real(8) :: eps
   eps = 1.0d-6
   if (d_L <= eps) then
      arg = 1.0d20
   else
      arg = (a_r - r_L)/d_L
   endif

   if (arg >= log(HUGE(a_r))) then
      d_f_cut = 0.0d0
   else
      exp_r = dexp(arg)
      exp_r2 = 1.0d0 + exp_r
      d_f_cut = -exp_r/(d_L*exp_r2*exp_r2)
   endif
end function d_Fermi_cut_off


   
end module MD_general_tools
