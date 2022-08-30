! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2022
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with Stillinger-Weber potential
! [1] Stillinger, Weber, PRB 31, 5262 (1985)
! [2] Fan et al., PRB 92, 094301 (2015)
! [3] Zhou et al., PRB 88, 085309 (2013)

MODULE MD_Pot_SW

use Universal_constants
use Objects
use Geometries, only: Spherical_R
use MD_general_tools, only: three_body_cos_theta, find_which_potential, Fermi_cut_off, d_Fermi_cut_off, shortest_distance


implicit none

real(8) :: m_eps, m_one_third

parameter(m_eps = 1.0d-10) ! presision
parameter(m_one_third = 1.0d0/3.0d0)

 contains


subroutine SW_potential(MDPar, r, i, j, R_ij, numpar, U_ij, dU_ij, MD_atoms, MD_supce, Pot, Pot_3bdy, disp, ind_disp) ! Stillinger-Weber potential [1]
   type(SW), intent(in) :: MDPar  ! All prameters of SW potential
   real(8), intent(in) :: r   ! [A] distance between atoms i and j
   integer, intent(in) :: i, j  ! indices of the pair of interacting atoms
   real(8), dimension(3), intent(in) :: R_ij ! [A] projections of distance between atoms i and j
   type(Num_par), intent(in), target :: numpar   ! all numerical parameters
   real(8), dimension(:,:), intent(in) :: U_ij, dU_ij ! precalculated repeated contributions of SW potential
   type(Atom), dimension(:), intent(in) :: MD_atoms
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   real(8), intent(out) :: Pot, Pot_3bdy  ! [eV]
   real(8), intent(in), optional :: disp  ! [A] displacement of atom i
   integer, intent(in), optional :: ind_disp  ! along which axis
   !-----------------------
   integer :: k, KOP1, KOP2, KOP3, kk, i_pot
   real(8) :: v2, v3, v3_cur, f2, h, R_ik(3), cos_theta, Rs, ra, F3bdy_1, F3bdy_2(3), cos3, dcos_theta(3), cos_ij(3), cos_ik(3)
   real(8) :: r2, x, y, z, f_cut
   integer, pointer :: atom_3

   ! Get two-body potential:
   Rs = r/MDPar%sigma ! normalized value
   f2 = SW_f2(MDPar, Rs) ! below
   v2 = MDPar%E * f2 ! [eV]
   Pot = v2

   !--------------------------
   ! Get three-body potential:
   v3_cur = 0.0d0 ! potential to start with
   ! Find nearest neighbors it is interacting with:
   do k = 1, numpar%Neighbors_Num(i)    ! that's how many atoms the atom "i" interacts with
      atom_3 => numpar%Neighbors_list(i, k)     ! that's index of the third atom
      if (atom_3 /= j) then   ! the third atom cannot equal to the second atom
         ! And the distance between atoms i and k (atom_3):
         if (present(disp) .and. present(ind_disp)) then ! atom i is displaced:
            call shortest_distance(i, atom_3, MD_atoms, MD_supce, r2, x1=R_ik(1), y1=R_ik(2), z1=R_ik(3), &
                        disp = disp, ind_disp = ind_disp) ! module "MD_general_tools"
         else  ! no displacement, regular coordinates:
            call shortest_distance(i, atom_3, MD_atoms, MD_supce, r2, x1=R_ik(1), y1=R_ik(2), z1=R_ik(3)) ! module "MD_general_tools"
         endif
         if (r2 < 1.0d-6) print*, 'Error in SW_potential: distance too small:', r2, R_ik(:)

         ! Get the three-body cosine:
         cos_theta = three_body_cos_theta(R_ij, R_ik)   ! module "MD_general_tools"
         cos3 = cos_theta + m_one_third ! to reuse below


         ! Testing:
!           if ((i == 5) .and. (j == 1) .and. (atom_3 == 2)) then
!          if ((i == 5) .and. (j == 1)) then
!             print*, 'Rij', R_ij(:)
!             print*, 'Rik', R_ik(:)
!             print*, 'cos', cos_theta
!             print*, atom_3, cos_theta, U_ij(i,atom_3) * cos3**2
!          endif

         v3_cur = v3_cur + U_ij(i,atom_3) * cos3**2  ! part of Eq.(1) [3]
         ! Testing:
!          f_cut = Fermi_cut_off(r2, 3.5d0, 0.1d0) ! module "MD_general_tools"
!          v3_cur = v3_cur + U_ij(i,atom_3) * f_cut * cos3**2  ! part of Eq.(1) [3]

      endif ! (k /= j)
   enddo ! k = 1, numpar%Neighbors_Num(i)
   ! Save three-body part [3]:
   Pot_3bdy = U_ij(i,j) * v3_cur  ! Eq.(1) [3]

   nullify(atom_3)
end subroutine SW_potential




subroutine d_SW_potential(MDPar, r, i, j, R_ij, numpar, U_ij, dU_ij, Pot, Pot_3bdy, &
                           d_Pot_part, Force_3bdy, MD_atoms, MD_pots, MD_supce) ! Stillinger-Weber potential and force [1]
   type(SW), intent(in) :: MDPar  ! All prameters of SW potential
   real(8), intent(in) :: r   ! [A] distance between atoms i and j
   integer, intent(in) :: i, j  ! indices of the pair of interacting atoms
   real(8), dimension(3), intent(in) :: R_ij ! [A] projections of distance between atoms i and j
   type(Num_par), intent(in), target :: numpar   ! all numerical parameters
   real(8), dimension(:,:), intent(in) :: U_ij, dU_ij ! precalculated repeated contributions of SW potential
   real(8), intent(out) :: Pot, Pot_3bdy  ! [eV]
   real(8), intent(out) :: d_Pot_part
   real(8), dimension(3), intent(out) :: Force_3bdy
   type(Atom), dimension(:), intent(in) :: MD_atoms
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   !-----------------------
   integer :: k, KOP1, KOP2, KOP3, kk, i_pot
   real(8) :: v2, v3, v3_cur, f2, h, R_ik(3), Rs, ra, F3bdy_1, F3bdy_2(3), cos3, dcos_theta(3), cos_ij(3), cos_ik(3)
   real(8) :: f_cut, d_f_cut, Pot_3, Forc_3(3), R_jk(3), Rij, Rik, Rjk, cos_jk(3)
   real(8) :: cos_theta_ijk, cos_theta_kji, cos_theta_jik
   real(8) :: cos3_ijk, cos3_kji, cos3_jik, A, C(3)
   real(8) :: dF_ijk_dr_ij(3), dF_jik_dr_ji(3), dF_ijk_dr_ik(3), dF_kji_dr_ki(3)
   real(8), pointer :: r2, x, y, z
   integer, pointer :: atom_3

   ! Get two-body potential:
   Rs = r/MDPar%sigma ! normalized value
   f2 = SW_f2(MDPar, Rs) ! below
   v2 = MDPar%E * f2 ! [eV]
   Pot = v2

   ! Get two-body force:
   f2 = SW_d_f2(MDPar, Rs) ! below
   d_Pot_part = (MDPar%E * MDPar%A * f2)     ! [eV/atom]


   !--------------------------
   ! Get three-body potential:
   v3_cur = 0.0d0 ! potential to start with
   F3bdy_2(:) = 0.0d0   ! force to start with
   A = 0.0d0
   C = 0.0d0
   ! cosine direction for atoms i and j:
   cos_ij(:) = R_ij(:)/r
   Rij = r

   ! Find nearest neighbors it is interacting with:
   do k = 1, numpar%Neighbors_Num(i)    ! that's how many atoms the atom "i" interacts with
      atom_3 => numpar%Neighbors_list(i, k)     ! that's index of the third atom
      if (atom_3 /= j) then   ! the third atom cannot equal to the second atom
         r2 => numpar%Neighbors_Rij(i, k)     ! [A] shortest distance between the atoms "i" and "atom_3" (which has number k in the list)
         x => numpar%Neighbors_Xij(i, k, 1)   ! [A] shortest distance between the atoms "i" and "atom_3" along X
         y => numpar%Neighbors_Xij(i, k, 2)   ! [A] shortest distance between the atoms "i" and "atom_3" along Y
         z => numpar%Neighbors_Xij(i, k, 3)   ! [A] shortest distance between the atoms "i" and "atom_3" along Z
         if (r2 < 1.0d-6) print*, 'Error in SW_potential: distance too small:', r2, x, y, z
         R_ik = (/ x, y, z /)
         Rik = r2
         ! Direction cosines between atoms i and atom_3 (neighbor "k"):
         !cos_ik(:) = R_ik(:)/r2

         ! And the distance between atoms j (atom_2) and k (atom_3):
         call shortest_distance(j, atom_3, MD_atoms, MD_supce, Rjk, x1=R_jk(1), y1=R_jk(2), z1=R_jk(3)) ! module "MD_general_tools"

         ! Get the three-body cosine:
         cos_theta_ijk = three_body_cos_theta( R_ij,  R_ik)   ! module "MD_general_tools"
         cos_theta_jik = three_body_cos_theta(-R_ij,  R_jk)   ! module "MD_general_tools"

         cos3_ijk = cos_theta_ijk + m_one_third ! to reuse below
         cos3_jik = cos_theta_jik + m_one_third ! to reuse below

         Pot_3 = U_ij(i,atom_3)*cos3_ijk**2 ! part of Eq.(1) [3]
         v3_cur = v3_cur + Pot_3

         !------------
         ! Force:
         ! The derivative of the cos(teta) by R_ij,{x,y,z}:
         dF_ijk_dr_ij = d_cos_theta_dr(cos_theta_ijk, Rij, R_ij, Rik, R_ik)   ! below
         dF_jik_dr_ji = d_cos_theta_dr(cos_theta_jik, Rij, -R_ij, Rjk, R_jk)   ! below

         ! Derivative of V^(3), collect the terms under the sum sign (part of Eq.(1) [3]):
         !A = A + U_ij(i,atom_3)*cos3_ijk**2 + U_ij(j,atom_3)*cos3_jik**2
         A = A + Pot_3 + U_ij(j,atom_3)*cos3_jik**2

         C(:) = C(:) + U_ij(i, atom_3)*2.0d0*cos3_ijk*dF_ijk_dr_ij(:) - U_ij(j, atom_3)*2.0d0*cos3_jik*dF_jik_dr_ji(:)

      endif ! (k /= j)
   enddo ! k = 1, numpar%Neighbors_Num(i)
   ! Save three-body part [3]:
   Pot_3 = U_ij(i,j) * v3_cur
   Pot_3bdy = Pot_3 ! Eq.(1) [3]

   ! Collect all the V^(3) terms, including those outside of the sum sign:
   Force_3bdy(:) = dU_ij(i,j)*cos_ij(:)*A + U_ij(i,j) * C(:)

   nullify(r2, x, y, z, atom_3)
end subroutine d_SW_potential


pure function d_cos_theta_dr(cos_theta, Rij, r_ij, Rik, r_ik) result(dFdr)
   real(8), dimension(3) :: dFdr
   real(8), intent(in) :: cos_theta, Rij, Rik, r_ij(3), r_ik(3)
   if (Rik >= 0.0d0 .and. Rik >= 0.0d0) then
      !dFdr(:) = 1/Rij * (r_ik(:)/Rik - r_ij(:)/Rij * cos_theta) + 1/Rik * (r_ij(:)/Rij - r_ik(:)/Rik * cos_theta)
      dFdr(:) = 1/Rij * (r_ik(:)/Rik - r_ij(:)/Rij * cos_theta)
   else
      dFdr = 0.0d0
   endif
end function d_cos_theta_dr



pure function SW_f2(MDPar, r) result(f2)
   real(8) f2  ! two-body potential function [1]
   type(SW), intent(in) :: MDPar  ! All prameters of SW potential
   real(8), intent(in) :: r   ! R/sigma [dimensionless], distance between atoms i and j
   real(8) :: ra
   ra = r-MDPar%aa
   ! Get f_2 function via Eq.(2.3) [1]:
   if (ra >= -m_eps) then ! atoms too far apart
      f2 = 0.0d0
   else
      f2 = MDPar%A*(MDPar%B*(r**(-MDPar%p)) - (r**(-MDPar%q)))*exp(1.0d0/ra)
   endif
end function SW_f2


pure function SW_d_f2(MDPar, Rs) result(f2)
   real(8) f2  ! two-body potential function [1]
   type(SW), intent(in) :: MDPar  ! All prameters of SW potential
   real(8), intent(in) :: Rs   ! R/sigma [dimensionless], distance between atoms i and j
   real(8) :: ra
   ra = Rs - MDPar%aa
   if (ra >= -m_eps) then ! atoms too far apart
      f2 = 0.0d0
   else
      f2 = ( -MDPar%p*MDPar%B*(Rs**(-MDPar%p-1.0d0)) + MDPar%q*(Rs**(-MDPar%q-1.0d0)) - &
           (MDPar%B*(Rs**(-MDPar%p)) - Rs**(-MDPar%q))/ra**2 ) * exp(1.0d0/ra) / MDPar%sigma
   endif
end function SW_d_f2




pure function SW_U_ij(MDPar, R) result(U_ij)
   real(8) U_ij   ! three-body potential function [3]
   type(SW), intent(in) :: MDPar  ! All prameters of SW potential
   real(8), intent(in) :: R ! distance between the atoms i and j
   real(8) :: rij, rija

   ! Get normalized value:
   rij = R/MDPar%sigma
   rija = rij-MDPar%aa
   ! Construct the U_ij-function via Eq.(4) [3]:
   if (rija >= -m_eps) then   ! atoms too far apart
      U_ij = 0.0d0
   else
      U_ij = sqrt(MDPar%ll * MDPar%E) * exp(MDPar%gam/rija)
   endif
end function SW_U_ij



pure function SW_dU_ij(MDPar, R) result(dU_ij)
   real(8) dU_ij   ! derivative of three-body potential function [3]
   type(SW), intent(in) :: MDPar  ! All prameters of SW potential
   real(8), intent(in) :: R ! distance between the atoms i and j
   real(8) :: rij, rija

   ! Get normalized value:
   rij = R/MDPar%sigma
   rija = rij-MDPar%aa
   ! Construct the U_ij-function via Eq.(4) [3]:
   if (rija >= -m_eps) then   ! atoms too far apart
      dU_ij = 0.0d0
   else
      dU_ij = -sqrt(MDPar%ll * MDPar%E)/MDPar%sigma * MDPar%gam/(rija**2) * exp(MDPar%gam/rija)
   endif
end function SW_dU_ij





END MODULE MD_Pot_SW
