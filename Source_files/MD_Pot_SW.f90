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
use MD_general_tools, only: three_body_cos_teta, find_which_potential


implicit none

real(8) :: m_eps, m_one_third

parameter(m_eps = 1.0d-10) ! presision
parameter(m_one_third = 1.0d0/3.0d0)

 contains


function SW_potential(MDPar, r, i, j, R_ij, numpar, U_ij) result(Pot) ! Stillinger-Weber potential [1]
   real(8) :: Pot  ! [eV]
   type(SW), intent(in) :: MDPar  ! All prameters of SW potential
   real(8), intent(in) :: r   ! [A] distance between atoms i and j
   integer, intent(in) :: i, j  ! indices of the interacting atoms
   real(8), dimension(3), intent(in) :: R_ij ! [A] projections of distance between atoms i and j
   type(Num_par), intent(in), target :: numpar   ! all numerical parameters
   real(8), dimension(:,:), intent(in) :: U_ij ! repeated contributions for three-body potentials (e.g. SW)
   !-----------------------
   integer :: k, KOP2, KOP3
   real(8) :: v2, v3, f2, h, R_ik(3), cos_teta
   real(8), pointer :: r2, x, y, z
   integer, pointer :: atom_3

   ! Get two-body potential:
   f2 = SW_f2(MDPar, r/MDPar%sigma) ! below
   v2 = MDPar%E * f2 ! [eV]

   ! Get three-body potential:
   v3 = 0.0d0  ! to start with
   ! Find nearest neighbors it is interacting with:
   do k = 1, numpar%Neighbors_Num(i)    ! that's how many atoms the atom "i" interacts with
      atom_3 => numpar%Neighbors_list(i, k)     ! that's index of the atom, which "i" is interacting with
      if (atom_3 /= j) then
         r2 => numpar%Neighbors_Rij(i, k)     ! [A] shortest distance between the atoms "i" and "atom_2" (which has number j in the list)
         x => numpar%Neighbors_Xij(i, k, 1)   ! [A] shortest distance between the atoms "i" and "atom_2" along X
         y => numpar%Neighbors_Xij(i, k, 2)   ! [A] shortest distance between the atoms "i" and "atom_2" along Y
         z => numpar%Neighbors_Xij(i, k, 3)   ! [A] shortest distance between the atoms "i" and "atom_2" along Z
         if (r2 < 1.0d-6) print*, 'Error in SW_potential: distance too small:', r, x, y, z
         R_ik = (/ x,y,z /)

         ! Get the three-body cosine:
         cos_teta = three_body_cos_teta(R_ij, R_ik)   ! module "MD_general_tools"

         v3 = v3 + U_ij(i,atom_3)*(cos_teta + m_one_third)**2  ! part of Eq.(1) [3]
      endif ! (k /= j)
   enddo ! k = 1, numpar%Neighbors_Num(i)
   v3 = U_ij(i,j) * v3  ! Eq.(1) [3]

   ! Add two-body and three-body parts [3]:
   Pot = v2 + v3

   nullify(r2, x, y, z, atom_3)
end function SW_potential



subroutine d_SW_potential(MDPar, r, i, j, R_ij, numpar, U_ij, dU_ij, d_Pot_part, Force_3bdy) ! Forces for Stillinger-Weber potential [1]
   type(SW), intent(in) :: MDPar  ! All prameters of SW potential
   real(8), intent(in) :: r   ! [A] distance between atoms i and j
   integer, intent(in) :: i, j  ! indices of the interacting atoms
   real(8), dimension(3), intent(in) :: R_ij ! [A] projections of distance between atoms i and j
   type(Num_par), intent(in), target :: numpar   ! all numerical parameters
   real(8), dimension(:,:), intent(in) :: U_ij, dU_ij ! repeated contributions for three-body potentials (e.g. SW)
   real(8), intent(out) :: d_Pot_part
   real(8), dimension(3), intent(out) :: Force_3bdy
   !-----------------------------------------------
   integer :: k
   real(8) :: Rs, ra, f2, cos_ij(3), cos_ik(3), dcos_ik(3), R_ik(3), F3bdy_1(3), F3bdy_2(3), cos_teta, cos3
   real(8), pointer :: rk, x, y, z
   integer, pointer :: atom_3

   ! Get normalized value:
   Rs = r/MDPar%sigma
   ra = Rs - MDPar%aa

   ! Two-body contribution:
   if (ra >= -m_eps) then ! atoms too far apart
      f2 = 0.0d0
   else
      f2 = ( -MDPar%p*MDPar%B*Rs**(-MDPar%p-1.0d0) + MDPar%q*Rs**(-MDPar%q-1.0d0) - &
           (MDPar%B*Rs**(-MDPar%p) - Rs**(-MDPar%q))/ra**2 ) * exp(1.0d0/ra) / MDPar%sigma
   endif
   d_Pot_part = MDPar%E * MDPar%A * f2 ! [eV/atom]


   ! Three-body contribution:
   Force_3bdy(:) = 0.0d0   ! to start with
   F3bdy_1(:) = 0.0d0
   F3bdy_2(:) = 0.0d0

   ! cosine direction for atoms i and j:
   cos_ij(:) = R_ij(:)/r

   do k = 1, numpar%Neighbors_Num(i)    ! that's how many atoms the atom "i" interacts with
      atom_3 => numpar%Neighbors_list(i, k)     ! that's index of the atom, which "i" is interacting with
      if (atom_3 /= j) then
         rk => numpar%Neighbors_Rij(i, k)      ! [A] shortest distance between the atoms "i" and "atom_3" (which has number j in the list)
         x => numpar%Neighbors_Xij(i, k, 1)   ! [A] shortest distance between the atoms "i" and "atom_3" along X
         y => numpar%Neighbors_Xij(i, k, 2)   ! [A] shortest distance between the atoms "i" and "atom_3" along Y
         z => numpar%Neighbors_Xij(i, k, 3)   ! [A] shortest distance between the atoms "i" and "atom_3" along Z
         if (rk < 1.0d-6) print*, 'Error in d_SW_potential: distance too small:', rk, x, y, z
         R_ik = (/ x,y,z /)

         ! Direction cosines between atoms i and k:
         cos_ik(:) = R_ik(:)/rk

         ! Get the three-body cosine:
         cos_teta = three_body_cos_teta(R_ij, R_ik)   ! module "MD_general_tools"
         cos3 = cos_teta + m_one_third ! to reuse below

         ! The derivative of the cos(teta) by R_i,{x,y,z}:
         dcos_ik(:) = 1.0d0/r * (cos_ik(:) - cos_ij(:)*cos_teta) + 1.0d0/rk*(cos_ij(:) - cos_ik(:)*cos_teta)

         ! Derivative of V^(3), collect the terms under the sum sign:
         F3bdy_1(:) = F3bdy_1(:) + U_ij(i,atom_3)*cos3**2  ! part of Eq.(1) [3]
         F3bdy_2(:) = F3bdy_2(:) + dU_ij(i,atom_3)*cos3**2*cos_ik(:) + U_ij(i,atom_3)*2.0d0*cos3*dcos_ik(:)
      endif ! (atom_3 /= j)
   enddo ! k
   ! Collect all the V^(3) terms, including those outside of the sum sign:
   Force_3bdy(:) = F3bdy_1(:)*dU_ij(i,j)*cos_ij(:) + F3bdy_2(:)*U_ij(i,j)

   nullify(rk, x, y, z, atom_3)
end subroutine d_SW_potential



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
      f2 = MDPar%A*(MDPar%B*r**(-MDPar%p) - r**(-MDPar%q))*exp(1.0d0/ra)
   endif
end function SW_f2


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
      dU_ij = - sqrt(MDPar%ll * MDPar%E) * MDPar%gam/rija**2 * exp(MDPar%gam/rija)/MDPar%sigma
   endif
end function SW_dU_ij





END MODULE MD_Pot_SW
