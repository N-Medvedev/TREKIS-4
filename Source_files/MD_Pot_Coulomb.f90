! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2021
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with Coulomb potential
! For Ewalds summation, it follows the example of Allen and Tildesley:
! [1] https://github.com/Allen-Tildesley/examples/blob/master/ewald_module.f90
! Which accompanies the book "Computer Simulation of Liquids", second edition, 2017
! [2] http://micro.stanford.edu/mediawiki/images/4/46/Ewald_notes.pdf
! For Wolf's et al. treatement of truncated Coulomb, see:
! [3] D. Wolf, et al., J. Chem. Phys. 110, 8254 (1999); https://doi.org/10.1063/1.478738
! [4] C. J. Fennell and J. D. Gezelter, J. Chem. Phys. 124, 234104 (2006)



MODULE MD_Pot_Coulomb
use Universal_constants
use Objects
use MD_general_tools, only : shortest_distance, find_which_potential
use Little_subroutines, only: print_time_step

implicit none

real(8) :: m_k, m_2Pi2, m_sqrtPi

parameter(m_k = g_ke * g_e * 1.0d10)  ! Constant in the Coulomb law, converting potential into [eV]
parameter(m_2Pi2 = g_2Pi*g_2Pi)     ! 2*Pi^2
parameter(m_sqrtPi = sqrt(g_Pi))    ! sqrt(Pi)

 contains


pure function Coulomb_Wolf_pot(q1, q2, r, Rc, alpha) result(WPot) ! truncated Coulomb potential [4]
   real(8) :: WPot  ! [eV]
   real(8), intent(in) :: q1, q2    ! [e] charges
   real(8), intent(in) :: r, Rc ! [A] interatomic distance, truncation distance
   real(8), intent(in) :: alpha ! truncation parameter
   real(8) :: Rc2, a2, term
   if (r < Rc) then
      Rc2 = Rc*Rc
      a2 = alpha*alpha
      term = erfc(alpha*Rc)/Rc
      WPot = m_k*q1*q2 * (erfc(alpha*r)/r - term + &
         (term/Rc + 2.0d0*alpha/m_sqrtPi*exp(-a2*Rc2)/Rc)*(r-Rc) )  ! [eV]
   else
      WPot = 0.0d0
   endif
end function Coulomb_Wolf_pot


pure function Coulomb_Wolf_self_term(q1, Rc, alpha) result(SelfPot)   ! Self-term, Eq.(5.13) [3]
   real(8) :: SelfPot  ! [eV]
   real(8), intent(in) :: q1    ! [e] charges
   real(8), intent(in) :: Rc    ! [A] interatomic distance, truncation distance
   real(8), intent(in) :: alpha ! truncation parameter
   SelfPot = m_k*q1*q1*(erfc(alpha*Rc)/Rc + 2.0d0*alpha/m_sqrtPi)    ! [eV]
end function Coulomb_Wolf_self_term


pure function d_Coulomb_Wolf_pot(q1, q2, r, Rc, alpha) result(dWPot) ! derivative truncated Coulomb potential, Eq.(5.22) [3]
   real(8) :: dWPot ! [eV/A]
   real(8), intent(in) :: q1, q2    ! [e] charges
   real(8), intent(in) :: r, Rc ! [A] interatomic distance, truncation distance
   real(8), intent(in) :: alpha ! truncation parameter
   real(8) :: r2, Rc2, a2
   if (r < Rc) then
      r2 = r*r
      Rc2 = Rc*Rc
      a2 = alpha*alpha
      dWPot = -m_k*q1*q2*( erfc(alpha*r)/r2 - erfc(alpha*Rc)/Rc2 + &
               2.0d0*alpha/m_sqrtPi*(exp(-a2*r2)/r - exp(-a2*Rc2)/Rc) )  ! [eV/A]
   else
      dWPot = 0.0d0
   endif
end function d_Coulomb_Wolf_pot




pure function Coulomb_Wolf_pot_simple(q1, q2, r, Rc) result(WPot) ! truncated Coulomb potential, Eq.(3.9) in [3]
   real(8) :: WPot  ! [eV]
   real(8), intent(in) :: q1, q2    ! [e] charges
   real(8), intent(in) :: r, Rc ! [A] interatomic distance, truncation distance
   if (r < Rc) then
      WPot = m_k*q1*q2*(1.0d0/r - 1.0d0/Rc)  ! [eV]
   else
      WPot = 0.0d0
   endif
end function Coulomb_Wolf_pot_simple


pure function Coulomb_Wolf_self_term_simple(q1, Rc) result(SelfPot)   ! Self-term, Eq.(3.15) [3]
   real(8) :: SelfPot  ! [eV]
   real(8), intent(in) :: q1    ! [e] charges
   real(8), intent(in) :: Rc    ! [A] interatomic distance, truncation distance
   SelfPot = m_k*q1*q1/Rc    ! [eV]
end function Coulomb_Wolf_self_term_simple


pure function d_Coulomb_Wolf_pot_simple(q1, q2, r, Rc) result(dWPot) ! derivative truncated Coulomb potential, Eq.(3.17) [3]
   real(8) :: dWPot ! [eV/A]
   real(8), intent(in) :: q1, q2    ! [e] charges
   real(8), intent(in) :: r, Rc ! [A] interatomic distance, truncation distance
   if (r < Rc) then
      dWPot = -m_k*q1*q2*(1.0d0/(r*r) - 1.0d0/(Rc*Rc))  ! [eV/A]
   else
      dWPot = 0.0d0
   endif
end function d_Coulomb_Wolf_pot_simple




subroutine Ewalds_Coulomb(numpar, MD_atoms, MD_supce, MD_pots)
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   class(MD_potential), dimension(:,:), intent(in), target :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   !-----------------------------
   integer :: nk, i, j, Nat, i_pot, N_pot, KOP1, KOP2, k_sq_max, kx, ky, kz, k_sq
   real(8) :: x, y, z, r, cos_a, cos_b, cos_c, Pot_part, d_Pot_part, kappa, b, kr_sq, Force(3), kap_r, erfc_kap
   real(8) :: prefac, prefac_f, k_r(3), tot_derive, eps, U_cur
   real(8), dimension(size(MD_atoms)) :: U_array
   real(8), dimension(size(MD_atoms),3) :: F_array
   ! Consider k-vectors within a sphere bounded by the cube (-nk,+nk) in each component:
   real(8), dimension(:), allocatable :: kfac     ! k-dept quantities (k_sq_max)
   real(8), dimension(:), allocatable :: q     ! charges of all atoms used in the Coulomb potential
   complex, dimension(:,:), allocatable :: eikx, eiky, eikz
   complex, dimension(:,:,:), allocatable :: Sk
   complex :: term

   eps = 1.0d-12    ! precision
   Nat = size(MD_atoms) ! all atoms
   allocate(q(Nat))

   ! Define Ewalds parameter:
   kappa = 6.0d0 / maxval(MD_supce%vect(:))    ! [1/A] Taken from [1]
   b = 0.25d0/(kappa*kappa) ! [A^2] parameter in the exponent of Ewalds k-part
   nk = 8           ! Taken from [1]
   k_sq_max = nk**2         ! maximal k for sum in Ewalds k-space part
   prefac = 0.5d0*g_e*1.0d10/(MD_supce%V*g_e0) ! factor in k-space Ewalds contribution
   prefac_f = g_e*1.0d20/(MD_supce%V*g_e0) ! factor in k-space Ewalds contribution

!    call print_time_step('Ewalds_Coulomb:', 0.0, msec=.true.)   ! module "Little_subroutines"

   ! Precalculate some parameters:
   if (.not.allocated(numpar%Ewalds_kfac)) then ! first use, define parameters for integratino of k-space Ewalds
      allocate( numpar%Ewalds_kfac(k_sq_max), source=0.0d0 ) ! Allocate module data array
      do kx = 0, nk
         do ky = 0, nk
            do kz = 0, nk
               k_sq = kx**2 + ky**2 + kz**2
               if ( ( k_sq <= k_sq_max ) .and. ( k_sq /= 0 ) ) then ! to ensure it is within range
                   kr_sq = m_2Pi2*( (dble(kx)/MD_supce%vect(1))**2 + (dble(ky)/MD_supce%vect(2))**2 + (dble(kz)/MD_supce%vect(3))**2 ) ! [1/A^2]
                   numpar%Ewalds_kfac(k_sq) = exp( -b * kr_sq ) / kr_sq ! [A^2] stored expression for later reuse
               endif
            enddo ! kz
         enddo ! ky
      enddo ! kx
   endif
   ! Allocate arrays for k-space Ewalds calculations:
   allocate(eikx(Nat,-nk:nk), source=cmplx(0.0d0,0.0d0))
   allocate(eiky(Nat,-nk:nk), source=cmplx(0.0d0,0.0d0))
   allocate(eikz(Nat,-nk:nk), source=cmplx(0.0d0,0.0d0))

!    call print_time_step('Ewalds_Coulomb:', 1.0, msec=.true.)   ! module "Little_subroutines"

   ! Fill the arrays of exponents for Sk:
   call get_eik_Ewalds(MD_atoms, MD_supce%vect(:), nk, eikx, eiky, eikz)  ! below

!    call print_time_step('Ewalds_Coulomb:', 2.0, msec=.true.)   ! module "Little_subroutines"

   U_array = 0.0d0  ! to start with
   F_array = 0.0d0  ! to start with
   !$omp parallel private (i, j, r, x, y, z, kap_r, erfc_kap, cos_a, cos_b, cos_c, KOP1, KOP2, Force)
   ! 1) Real part of the Ewalds sum:
   !$omp do schedule(dynamic)
   do i = 1, Nat    ! all atoms in the supercell
      do j = 1, Nat ! upper triangle of pairs of atoms
         if (i /= j) then   ! no self-interaction
            ! Get the distance between the atoms, accounting for periodic boundaries:
            call shortest_distance(i, j, MD_atoms, MD_supce, r, x1=x, y1=y, z1=z)    ! module "MD_general_tools"
            if (r < 1.0d-6) print*, 'Error#2 in get_MD_pot_n_force: distance too small:', r, x, y, z

            kap_r = kappa*r    ! argumet in the erfc-function for real part of Ewalds sum
            if (kap_r > eps) then ! large enough energy to consider
               ! Direction cosines between the two atoms (to calculate forces):
               cos_a = x/r
               cos_b = y/r
               cos_c = z/r
               erfc_kap = ERFC(kap_r)   ! cut-off ERFC function of real Ewalds part

               ! Find indices of the potential(s) used for this pair of atoms:
               call find_which_potential(MD_atoms, i, j, MD_pots, KOP1, KOP2)    ! module "MD_general_tools"

               Force = 0.0d0
               ! Get potential and force for real part of the Ewald sum:
               call add_potential_and_force(U_array(i), Force, MD_pots(KOP1,KOP2)%Set, &
                                            r, erfc_kap, kappa, kap_r, cos_a, cos_b, cos_c)    ! below
               F_array(i,:) = F_array(i,:) + Force

            endif ! (kap_r > eps)
         endif ! (i /= j)
      enddo ! j
      ! Save the data into the atomic arrays:
      MD_atoms(i)%U = MD_atoms(i)%U + U_array(i)
      MD_atoms(i)%Force(:) = MD_atoms(i)%Force(:) + F_array(i,:)
   enddo ! i
   !$omp enddo
   !$omp end parallel

!    call print_time_step('Ewalds_Coulomb:', 3.0, msec=.true.)   ! module "Little_subroutines"

   !-----------------------------------
   ! 2) k-space part of the Ewalds sum:
   ! a) construct S(k)=exp(ik.r):
   ! Get the array of charges first:
   call define_charge_array(MD_atoms, MD_pots, q)   ! below
   ! Calcualte Sk itself now:
   call get_Sk_Ewalds(Nat, nk, eikx, eiky, eikz, q, k_sq_max, Sk)    ! below

!    call print_time_step('Ewalds_Coulomb:', 4.0, msec=.true.)   ! module "Little_subroutines"

   ! b) get k-space contribution into Ewalds summation:
   !$omp parallel private (i, Force, U_cur, kx, ky, kz, k_sq, Pot_part, k_r, tot_derive)
   !$omp do schedule(dynamic)
   do i = 1, Nat    ! all atoms in the supercell
      Force = 0.0d0 ! to start summing up total forces
      U_cur = 0.0d0 ! to start
      do kx = -nk, nk
         do ky = -nk, nk
            do kz = -nk, nk
               k_sq = kx**2 + ky**2 + kz**2
               if ( ( k_sq <= k_sq_max ) .and. ( k_sq /= 0 ) ) then ! to ensure it is within range
                  Pot_part = numpar%Ewalds_kfac(k_sq) * &
                             dble(Sk(kx,ky,kz)*eikx(i,-kx)*eiky(i,-ky)*eikz(i,-kz))
                  ! Get the total energy:
                  !MD_atoms(i)%U = MD_atoms(i)%U + Pot_part * 2.0d0
                  U_cur = U_cur + Pot_part  ! [eV]
!                   write(*,'(i5,i3,i3,i3,es,es,es)') i, kx, ky, kz, Pot_part, prefac, MD_supce%V

                  ! k-vector:
                  k_r(1) = g_2Pi*dble(kx)/MD_supce%vect(1)  ! [1/A]
                  k_r(2) = g_2Pi*dble(ky)/MD_supce%vect(2)  ! [1/A]
                  k_r(3) = g_2Pi*dble(kz)/MD_supce%vect(3)  ! [1/A]

                  ! Get the interaction distance, accounting for periodic boundaries:
                  tot_derive = dble( g_CI * ( Sk(kx,ky,kz)*eikx(i,-kx)*eiky(i,-ky)*eikz(i,-kz) - &
                                          CONJG(Sk(kx,ky,kz))*eikx(i,kx)*eiky(i,ky)*eikz(i,kz) ) )   ! Eq.(42) in [2]
                  ! Add to total force contribution:
                  Force(:) = Force(:) - tot_derive*numpar%Ewalds_kfac(k_sq)*k_r(:)    ! add to total force along all axes
               endif ! ( ( k_sq <= k_sq_max ) .and. ( k_sq /= 0 ) )
            enddo ! kz
         enddo ! ky
      enddo ! kx
      ! Add to the total potential and force on atom i:
      MD_atoms(i)%U = MD_atoms(i)%U + prefac*q(i)*U_cur * 2.0d0 ! [eV] factor of 2 is to double-count, because total U includes 1/2
      MD_atoms(i)%Force(:) = MD_atoms(i)%Force(:) + prefac*q(i)*Force(:)  ! [eV/A]

      ! Subtract self part of k-space sum
      MD_atoms(i)%U = MD_atoms(i)%U - kappa*q(i)**2/m_sqrtPi * 2.0d0    ! [eV] same double-counting
   enddo ! i
   !$omp enddo
   !$omp end parallel

!    call print_time_step('Ewalds_Coulomb:', 5.0, msec=.true.)   ! module "Little_subroutines"

   if (allocated(q)) deallocate(q)
   if (allocated(kfac)) deallocate(kfac)
   if (allocated(eikx)) deallocate(eikx)
   if (allocated(eiky)) deallocate(eiky)
   if (allocated(eikz)) deallocate(eikz)
   if (allocated(Sk)) deallocate(Sk)
end subroutine Ewalds_Coulomb



subroutine add_potential_and_force(Upot, Force, MD_pots_Set, r, erfc_kap, kappa, kap_r, cos_a, cos_b, cos_c)
   real(8), intent(inout) :: Upot
   real(8), dimension(3), intent(inout) :: Force
   type(MD_pot_set), dimension(:), intent(in) :: MD_pots_Set   ! set of MD parameters
   real(8), intent(in) :: r, erfc_kap, kappa, kap_r, cos_a, cos_b, cos_c
   real(8) :: d_Pot_part, Pot_part, Force_one(3), tot_derive
   integer :: i_pot, N_pot
   !*/ Remark: ASSOCIATE construct inside of nested OMP do-loop has problems,
   !*/ that is why ASSOCIATE construct must be in a separate subroutine
   !*/ called from within a nested OMP do-loop.

   ! Potential and force for pure Coulomb:
   Force_one = 0.0d0  ! to start summing up total forces
   N_pot = size(MD_pots_Set)
   do i_pot = 1, N_pot
      ASSOCIATE (MDPot => MD_pots_Set(i_pot)%Par) ! this is the sintax we have to use to check the type of defined class
         select type(MDPot)
         type is (Coulomb_Ewald)   ! Coulomb potential as a long range one
            ! Get potential:
            Pot_part = Bare_Coulomb_pot(MDPot%Z1, MDPot%Z2, r)        ! [eV] below
            d_Pot_part = d_Bare_Coulomb_pot(MDPot%Z1, MDPot%Z2, r)    ! [eV/A] below
            ! Get the total energy:
            Upot = Upot + Pot_part*erfc_kap  ! [eV] Real part of Ewalds method for Coulomb potential
            ! Get projection of forces (F=-dU(r_ij)/dr_ij):
            tot_derive = -d_Pot_part*erfc_kap + Pot_part*(2.0d0*kappa/m_sqrtPi*exp(-kap_r*kap_r))
            Force_one(1) = Force_one(1) + tot_derive * cos_a    ! [eV/A] add to total force along X
            Force_one(2) = Force_one(2) + tot_derive * cos_b    ! [eV/A] add to total force along Y
            Force_one(3) = Force_one(3) + tot_derive * cos_c    ! [eV/A] add to total force along Z
         endselect
      END ASSOCIATE
   enddo ! i_pot
   ! Add to the total force on atom i:
   Force(:) = Force(:) + Force_one(:)
end subroutine add_potential_and_force


subroutine define_charge_array(MD_atoms, MD_pots, q)
   type(Atom), dimension(:), intent(in) :: MD_atoms	! all atoms in MD as objects
   type(MD_potential), dimension(:,:), intent(in), target :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   real(8), dimension(:), intent(inout) :: q    ! charges of all atoms used in Coulomb potential
   integer :: i, Nat, KOP1, KOP2, N_pot, i_pot
   class(MD_pot), pointer :: MDPot

   Nat = size(MD_atoms) ! total number of atoms

   !$omp parallel private (i, KOP1, KOP2, N_pot, i_pot, MDPot)
   !$omp do schedule(dynamic)
   do i = 2, Nat
      ! Find indices of the potential(s) is used for this pair of atoms:
      call find_which_potential(MD_atoms, i, 1, MD_pots, KOP1, KOP2)    ! module "MD_general_tools"

      N_pot = size(MD_pots(KOP1,KOP2)%Set)
      do i_pot = 1, N_pot
         ASSOCIATE (MDPot => MD_pots(KOP1,KOP2)%Set(i_pot)%Par) ! this is the sintax we have to use to check the type of defined class
            select type(MDPot)
            type is (Coulomb_Ewald)   ! Coulomb potential is long range
               ! Fill array of charges:
               q(i) = MDPot%Z1    ! we will use it below for k-space contribution of Ewalds
               if (i == 2) q(1) = MDPot%Z2
!                write(*,'(i6,f,i3,i3,a4)') i, q(i), KOP1, KOP2, MD_atoms(i)%Name
            endselect
         END ASSOCIATE
      enddo ! i_pot
   enddo ! i
   !$omp enddo
   !$omp end parallel
end subroutine define_charge_array



subroutine get_Sk_Ewalds(Nat, nk, eikx, eiky, eikz, q, k_sq_max, Sk)
   integer, intent(in) :: Nat, nk    ! number of atoms and summation points
   complex, dimension(Nat,-nk:nk), intent(in) :: eikx, eiky, eikz   ! exp(i*k*r)
   real(8), dimension(:), intent(in) :: q     ! charges of all atoms used in the Coulomb potential
   integer, intent(in) :: k_sq_max  ! max k for integration in Ewalds
   complex, dimension(:,:,:), allocatable, intent(inout) :: Sk
   !------------------------
   integer :: kx, ky, kz, k_sq

   allocate(Sk(-nk:nk,-nk:nk,-nk:nk))

   !$omp parallel private (kx,ky,kz,k_sq)
   !$omp do schedule(dynamic)
   do kx = -nk, nk
      do ky = -nk, nk
         do kz = -nk, nk
            k_sq = kx**2 + ky**2 + kz**2
            if ( ( k_sq <= k_sq_max ) .and. ( k_sq /= 0 ) ) then ! to ensure it is within range
               Sk(kx, ky, kz) = SUM(q(:)*eikx(:,kx)*eiky(:,ky)*eikz(:,kz))  ! sum over all atoms j
            endif
         enddo ! kz
      enddo ! ky
   enddo ! kx
   !$omp enddo
   !$omp end parallel
end subroutine get_Sk_Ewalds



subroutine get_eik_Ewalds(MD_atoms, MD_supce_vec, nk, eikx, eiky, eikz) ! Following method from [1]
   type(Atom), dimension(:), intent(in) :: MD_atoms	! all atoms in MD as objects
   real(8), dimension(3), intent(in) :: MD_supce_vec    ! [A] orthogonal supercell sizes along X, Y, Z
   integer, intent(in) :: nk    ! number of summation points
   complex, dimension(size(MD_atoms),-nk:nk), intent(inout) :: eikx, eiky, eikz   ! exp(i*k*r)
   real(8), dimension(size(MD_atoms)) :: x2pi, y2pi, z2pi
   integer :: i, j, k

   ! Calculate kx, ky, kz = 0, 1 explicitly
   eikx(:,0) = cmplx(1.0d0, 0.0d0)
   eiky(:,0) = cmplx(1.0d0, 0.0d0)
   eikz(:,0) = cmplx(1.0d0, 0.0d0)
   x2pi = -g_2Pi * MD_atoms(:)%R(1) / MD_supce_vec(1)
   y2pi = -g_2Pi * MD_atoms(:)%R(2) / MD_supce_vec(2)
   z2pi = -g_2Pi * MD_atoms(:)%R(3) / MD_supce_vec(3)
   eikx(:,1) = exp(g_CI*x2pi)   ! exp(-i.k_x.r_x)
   eiky(:,1) = exp(g_CI*y2pi)   ! exp(-i.k_y.r_y)
   eikz(:,1) = exp(g_CI*z2pi)   ! exp(-i.k_z.r_z)
   ! Calculate remaining positive kx, ky and kz by recurrence
   do k = 2, nk
      eikx(:,k) = eikx(:,k-1) * eikx(:,1)
      eiky(:,k) = eiky(:,k-1) * eiky(:,1)
      eikz(:,k) = eikz(:,k-1) * eikz(:,1)
   enddo

   ! Negative k values are complex conjugates of positive ones
   ! We do not need negative values of kx
   eikx(:,-nk:-1) = CONJG ( eikx(:,nk:1:-1) )
   eiky(:,-nk:-1) = CONJG ( eiky(:,nk:1:-1) )
   eikz(:,-nk:-1) = CONJG ( eikz(:,nk:1:-1) )
end subroutine get_eik_Ewalds




pure function Bare_Coulomb_pot(q1, q2, r) result(CPot) ! Coulomb potential
   real(8) :: CPot  ! [eV]
   real(8), intent(in) :: q1, q2    ! [e] charges
   real(8), intent(in) :: r ! [A] interatomic distance
   CPot = m_k*q1*q2/r  ! [eV]
end function Bare_Coulomb_pot


pure function d_Bare_Coulomb_pot(q1, q2, r) result(dCPot) ! derivative of Coulomb potential
   real(8) :: dCPot ! [eV/A]
   real(8), intent(in) :: q1, q2    ! [e] charges
   real(8), intent(in) :: r ! [A] interatomic distance
   dCPot = -m_k*q1*q2/(r*r)  ! [eV/A]
end function d_Bare_Coulomb_pot



END MODULE MD_Pot_Coulomb
