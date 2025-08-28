! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018-2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains definitions of different geometric objects used for targets, and coordinate transformations between different systems
! References:
! [1] T.M. Jenkins, W.R. Nelson, A. Rindi, "Monte Carlo Transport of Electrons and Photons" (1988) Plenum Press, New York


module Geometries
use Universal_constants

implicit none

real(8) :: m_tollerance_eps

parameter (m_tollerance_eps = 1.0d-6)    ! precision of any boundary crossing
 
 
 contains

 

pure function Spherical_R(x,y,z) result(R)
   real(8) R
   real(8), intent(in) :: x, y, z
   R = sqrt(x*x + y*y + z*z)
end function Spherical_R

pure function Cylindrical_R(x,y) result(R)
   real(8) R
   real(8), intent(in) :: x, y
   R = sqrt(x*x + y*y)
end function Cylindrical_R




subroutine fit_parabola_to_3points(x1, y1, x2, y2, x3, y3, A, B, C)  ! find coefficients of parabola by given three points
   real(8), intent(in) :: x1, x2, x3    ! given three points x coordinates
   real(8), intent(in) :: y1, y2, y3    ! function values for the given three points
   real(8), intent(out) :: A, B, C      ! parabola coefficients
   real(8) :: denom, x1x2, x2x3, x1x3, x3y1y2, x2y1y3, x1y2y3
   x1x2 = x1 - x2
   x2x3 = x2 - x3
   x1x3 = x1 - x3
   x3y1y2 = x3*(y1 - y2)
   x2y1y3 = x2*(y1 - y3)
   x1y2y3 = x1*(y2 - y3)
   denom = x1x2*x1x3*x2x3
   A = (-x3y1y2 + x2y1y3 - x1y2y3) / denom
   B = (x3*x3y1y2 - x2*x2y1y3 + x1*x1y2y3) / denom
   C = (x2*x3*x2x3*y1 - x3*x1*x1x3*y2 + x1*x2*x1x2*y3) / denom
end subroutine fit_parabola_to_3points



subroutine fit_line_to_2points(x1, y1, x2, y2, A, B)  ! find coefficients of a line by given two points
   real(8), intent(in) :: x1, x2    ! given 2 points x coordinates
   real(8), intent(in) :: y1, y2    ! function values for the given 2 points
   real(8), intent(out) :: A, B      ! line coefficients
   real(8) :: x1x2
   x1x2 = x2 - x1
   A = (x2*y1 - x1*y2) / x1x2
   B = (y2 - y1) / x1x2
end subroutine fit_line_to_2points




!------------------------------------------------------------------------------------ 
!Find the velosity projection onto the normal to the surface:
pure function scalar_projection(V_to_project, V_normal) result(Vnorm)
   real(8) Vnorm    ! component of V_to_project projected onto V_normal
   real(8), dimension(:), intent(in) :: V_to_project    ! vector to be projected
   real(8), dimension(:), intent(in) :: V_normal        ! vector to project on
   real(8) abs_V_norm
   abs_V_norm = SQRT( SUM(V_normal(:)*V_normal(:) ) )   ! normalization of vector we are projecting on
   ! Scalar projection:
   Vnorm = DOT_PRODUCT(V_to_project, V_normal) / abs_V_norm ! intrinsic fortran subroutine
end function scalar_projection

 
!------------------------------------------------------------------------------------
! Points of intersection of a vector with a surface:

pure subroutine where_crossing_cylinder(R, V, Cylinder_bottom, Cylinder_top, Cylinder_R, T)
   real(8), dimension(3), intent(in) :: R, V    ! [A] coordinates (with respect to "center" of cylinder) and [A/fs] velosity of the particle
   real(8), intent(in) :: Cylinder_bottom, Cylinder_top, Cylinder_R    ! dimensions of the cylinder, aligned along Z-axis
   real(8), intent(out) :: T ! [fs] time to the crossing-boundary point
   !========================
   real(8) :: t_cyl, L_coord, t_cur, V_rad, R_coord
   real(8) :: Cx, Cy, Cz, n1, n2, n3
   ! Find intersection with the infinite cylliner:
   call intersection_with_cylinder(R(1), R(2), V(1), V(2), Cylinder_R, t_cyl)   ! below
   L_coord = V(3)*t_cyl
   ! Check if this coordinate along Z is within the cylinder:
   if ( (L_coord >= Cylinder_bottom) .and. (L_coord <= Cylinder_top) ) then ! hit
      T = t_cyl
   else ! possible miss, check the bases of the cylinder
      ! Check bottom plane:
      Cx = 0.0d0
      Cy = 0.0d0
      Cz = Cylinder_bottom
      call define_normal_to_plane(Cx, Cy, Cz, n1, n2, n3)   ! below
      call intersection_with_plane(R(1), R(2), R(3), V(1), V(2), V(3), n1, n2, n3, Cx, Cy, Cz, t_cur) ! below
      V_rad = sqrt(V(1)*V(1) + V(2)*V(2))
      R_coord = V_rad * t_cur
      if (R_coord <= Cylinder_R) then   ! hit
         T = t_cur
      else  ! possible miss
         T = 1.0d25
      endif
      ! Check top plane:
      Cx = 0.0d0
      Cy = 0.0d0
      Cz = Cylinder_top
      call define_normal_to_plane(Cx, Cy, Cz, n1, n2, n3)   ! below
      call intersection_with_plane(R(1), R(2), R(3), V(1), V(2), V(3), n1, n2, n3, Cx, Cy, Cz, t_cur) ! below
      V_rad = sqrt(V(1)*V(1) + V(2)*V(2))
      R_coord = V_rad * t_cur
      if (R_coord <= Cylinder_R) then   ! hit
         if (T > t_cur) then    ! hits this plane first
            T = t_cur   ! save the time to this plane
         endif
      endif
   endif
end subroutine where_crossing_cylinder



pure subroutine where_crossing_box_boundary(R, V, box_start_x, box_end_x, box_start_y, box_end_y, box_start_z, box_end_z, &
   generation, periodic_x, periodic_y, periodic_z, T, neutral)
   real(8), dimension(3), intent(in) :: R, V    ! [A] coordinates and [A/fs] velosity of the particle
   real(8), intent(in) :: box_start_x, box_end_x, box_start_y, box_end_y, box_start_z, box_end_z    ! dimensions of the simulation box
   integer, intent(in) :: generation, periodic_x, periodic_y, periodic_z   ! generation of prticles; periodic or not boundary along x, y, z
   real(8), intent(out) :: T ! [fs] time to the crossing-boundary point
   logical, intent(in), optional :: neutral ! if it's neutral particle, periodic conditions can be simplified
   !========================
   real(8) :: t_cur, L_cur
   real(8) :: Cx, Cy, Cz, n1, n2, n3
   logical :: simple_cond
   ! Simulation box is always parallelepiped, so crossing its planes are here:
   
   ! Along Z:
   if (present(neutral)) then
      if ( (neutral) .and. (generation > 0) .and. (periodic_z==1) ) then
         simple_cond =.true.
      else
         simple_cond =.false.
      endif
   else     ! normal particle and boundary conditions:
      simple_cond =.false.
   endif
   ! Top surface:
   Cx = 0.0d0
   Cy = 0.0d0
   Cz = box_start_z
   if (simple_cond) then
      T = 1.1d25    ! don't account for the periodic boundary crossing
   else     ! normal particle and boundary conditions:
      call define_normal_to_plane(Cx, Cy, Cz, n1, n2, n3)   ! below
      call intersection_with_plane(R(1), R(2), R(3), V(1), V(2), V(3), n1, n2, n3, Cx, Cy, Cz, t_cur) ! below
      T = t_cur
!       print*, 'z1', t_cur
   endif
   ! Bottom surface:
   Cx = 0.0d0
   Cy = 0.0d0
   Cz = box_end_z
   if (.not.simple_cond) then
      call define_normal_to_plane(Cx, Cy, Cz, n1, n2, n3)   ! below
      call intersection_with_plane(R(1), R(2), R(3), V(1), V(2), V(3), n1, n2, n3, Cx, Cy, Cz, t_cur) ! below
      if ((t_cur > 0.0d0) .and. (T>t_cur)) T = t_cur
!       print*, 'z2', t_cur
   endif 
   
   ! Along X:
   if (present(neutral)) then
      if ( (neutral) .and. (generation > 0) .and. (periodic_x==1) ) then
         simple_cond =.true.
      else
         simple_cond =.false.
      endif
   else     ! normal particle and boundary conditions:
      simple_cond =.false.
   endif 
   ! Left surface:
   Cx = box_start_x
   Cy = 0.0d0
   Cz = 0.0d0
   if (.not.simple_cond) then
      call define_normal_to_plane(Cx, Cy, Cz, n1, n2, n3, ind=1)   ! below
      call intersection_with_plane(R(1), R(2), R(3), V(1), V(2), V(3), n1, n2, n3, Cx, Cy, Cz, t_cur) ! below
      if ((t_cur > 0.0d0) .and. (T>t_cur)) T = t_cur
!       print*, '--------------------'
!       print*, R(1), R(2), R(3)
!       print*, V(1), V(2), V(3)
!       print*, n1, n2, n3
!       print*, Cx, Cy, Cz
!       print*, 'x1', t_cur
!       print*, '--------------------'
   endif
   ! Right surface:
   Cx = box_end_x
   Cy = 0.0d0
   Cz = 0.0d0
   if (.not.simple_cond) then
      call define_normal_to_plane(Cx, Cy, Cz, n1, n2, n3, ind=1)   ! below
      call intersection_with_plane(R(1), R(2), R(3), V(1), V(2), V(3), n1, n2, n3, Cx, Cy, Cz, t_cur) ! below
      if ((t_cur > 0.0d0) .and. (T>t_cur)) T = t_cur
!       print*, 'x2', t_cur
   endif
   
   ! Along Y:
   if (present(neutral)) then
      if ( (neutral) .and. (generation > 0) .and. (periodic_y==1) ) then
         simple_cond =.true.
      else
         simple_cond =.false.
      endif
   else     ! normal particle and boundary conditions:
      simple_cond =.false.
   endif 
   ! Back surface:
   Cx = 0.0d0
   Cy = box_start_y
   Cz = 0.0d0
   if (.not.simple_cond) then
      call define_normal_to_plane(Cx, Cy, Cz, n1, n2, n3, ind=2)   ! below
      call intersection_with_plane(R(1), R(2), R(3), V(1), V(2), V(3), n1, n2, n3, Cx, Cy, Cz, t_cur) ! below
      if ((t_cur > 0.0d0) .and. (T>t_cur)) T = t_cur
!       print*, 'y1', t_cur
   endif
   ! Front surface:
   Cx = 0.0d0
   Cy = box_end_y
   Cz = 0.0d0
   if (.not.simple_cond) then
      call define_normal_to_plane(Cx, Cy, Cz, n1, n2, n3, ind=2)   ! below
      call intersection_with_plane(R(1), R(2), R(3), V(1), V(2), V(3), n1, n2, n3, Cx, Cy, Cz, t_cur) ! below
      if ((t_cur > 0.0d0) .and. (T>t_cur)) T = t_cur
!       print*, 'y2', t_cur
   endif
end subroutine where_crossing_box_boundary



pure subroutine define_normal_to_plane(Cx,Cy,Cz,n1,n2,n3,ind)
   real(8), intent(in) :: Cx,Cy,Cz  ! coordinates of the plane point closest to the coordinates origin
   real(8), intent(out) :: n1,n2,n3
   integer, optional, intent(in) :: ind
   real(8) :: Cnorm, Cx1, Cy1, Cz1, eps
   !eps = 1.0d-6
   eps = m_tollerance_eps
   if ( (ABS(Cx) < eps) .and. (ABS(Cy) < eps) .and. (ABS(Cz) < eps) ) then
      if (present(ind)) then
         select case (ind)
         case (1)       ! x-axis
            Cx1 = 1.0d0
            Cy1 = 0.0d0
            Cz1 = 0.0d0
         case (2)       ! y-axis
            Cx1 = 0.0d0
            Cy1 = 1.0d0
            Cz1 = 0.0d0
         case default   ! z-axis
            Cx1 = 0.0d0
            Cy1 = 0.0d0
            Cz1 = 1.0d0
         end select
      else  ! if nothing is known, asume z-axis 
         Cx1 = 0.0d0
         Cy1 = 0.0d0
         Cz1 = 1.0d0
      endif
   else
      Cx1 = Cx
      Cy1 = Cy
      Cz1 = Cz
   endif
   Cnorm = sqrt(Cx1*Cx1 + Cy1*Cy1 + Cz1*Cz1)
   n1 = Cx1/Cnorm
   n2 = Cy1/Cnorm
   n3 = Cz1/Cnorm
end subroutine define_normal_to_plane



 function which_box_wall_crossing(R, box_start_x, box_end_x, box_start_y, box_end_y, box_start_z, box_end_z) result(iw)
   real(8), dimension(3), intent(in) :: R   ! [A] coordinates
   real(8), intent(in) :: box_start_x, box_end_x, box_start_y, box_end_y, box_start_z, box_end_z
   integer :: iw    ! 1=top (Z_start), 2= bottom(Z_end), 3=left (Y_start), 4=right (Y_end), 5=back (X_start), 6=front (X_end)
   real(8) :: eps   ! precision
   !eps = 1.0d-6 ! how close to the boundary, to consider that it is crossing
   eps = m_tollerance_eps
   iw = 0   ! to start with
   if (R(1) >= box_end_x - eps) iw = 6
   if (R(1) <= box_start_x + eps) iw = 5
   if (R(2) >= box_end_y - eps) iw = 4
   if (R(2) <= box_start_y + eps) iw = 3
   if (R(3) >= box_end_z - eps) iw = 2
   if (R(3) <= box_start_z + eps) iw = 1
   if (iw == 0) then
      print*, 'Problem in which_box_wall_crossing:'
      print*, R
      print*, box_start_x, box_start_y, box_start_z
      print*, box_end_x, box_end_y, box_end_z
      pause 'DO NOT KNOW WHAT TO DO :-('
   endif
end function which_box_wall_crossing


pure subroutine intersection_with_plane(x, y, z, vx, vy, vz, n1, n2, n3, Cx, Cy, Cz, t)
   ! Subroutine follows the Section (17.1.1) from [1]
   real(8), intent(in) :: x, y, z   ! [A] coordinates of the particle
   real(8), intent(in) :: vx, vy, vz  ! [A/fs] velosity of the particle
   real(8), intent(in) :: n1, n2, n3  ! [A] coordinates of a point on the plane 
   real(8), intent(in) :: Cx, Cy, Cz  ! [A] normal vector to the plane
   real(8), intent(out) :: t    ! time when the intersection will take place (intersection point)
   ! t>0 means there is an intersection point; t<=0 means no intersection(in the future)
   !--------------------------
   real(8) :: nom, den, eps
   !eps = 1.0d-8 ! tollerance for parallel definition
   eps = m_tollerance_eps
   ! Define intersection point according to Eq.(17.8) [1]:
   nom = (Cx-x)*n1 + (Cy-y)*n2 + (Cz-z)*n3
   den = vx*n1 + vy*n2 + vz*n3
   if (abs(nom) < eps*1.0d-2) then  ! it just crossed the border
      t = eps
!       t = 1.0d20
   else
      !if ((abs(den) <= eps) .or. (abs(nom) <= eps)) then ! parallel to the plane, or just crossed the plane
      if (abs(den) <= eps) then ! parallel to the plane
         t = 1.0d20
      else
         t = nom/den
         if (t < 0.0d0) t = 2.0d20
      endif
   endif
end subroutine intersection_with_plane


pure subroutine intersection_with_cylinder(x, y, vx, vy, R, t)
   ! Subroutine follows the Section (17.1.3) from [1]
   real(8), intent(in) :: x, y   ! [A] coordinates of the particle
   real(8), intent(in) :: vx, vy  ! [A/fs] velosity of the particle
   real(8), intent(in) :: R  ! [A] radius of the cylinder located at (0,0,0) along Z-axis
   real(8), intent(out) :: t    ! time when the intersection will take place (intersection point)
   ! t>0 means there is an intersection point; t<=0 means no intersection(in the future)
   !--------------------------
   real(8) :: a, b, c, eps, t_cur(2)
   !eps = 1.0d-6 ! tolerance
   eps = m_tollerance_eps
   ! Define intersection point according to Eq.(17.13) [1]:
   ! Coefficients according Eq.(17.14) [1]:
   a = vx*vx + vy*vy
   b = x*vx + y*vy
   c = x*x + y*y - R*R
   t = -1.0d0    ! no crossing possible, to start with
   if (abs(a) >  eps) then  ! crossing is possible
      call solution_of_quadratic_equation(a, b, c, t_cur(1), t_cur(2))   ! below
      if (t_cur(1) > 0.0d0) then
         if (t_cur(2) > 0.0d0) then
            t = minval(t_cur(:))
         else
            t = t_cur(1)
         endif
      else
         if (t_cur(2) > 0.0d0) t = t_cur(2)
      endif
   endif ! (abs(a) >  eps)
end subroutine intersection_with_cylinder



pure subroutine intersection_with_sphere(x, y, z, vx, vy, vz, R, t)
   ! Subroutine derived equations following similar logic to the Section (17.1.3) from [1]
   real(8), intent(in) :: x, y, z   ! [A] coordinates of the particle
   real(8), intent(in) :: vx, vy, vz  ! [A/fs] velosity of the particle
   real(8), intent(in) :: R  ! [A] radius of the sphere centered at (0,0,0)
   real(8), intent(out) :: t    ! time when the intersection will take place (intersection point)
   ! t>0 means there is an intersection point; t<=0 means no intersection(in the future)
   !--------------------------
   real(8) :: a, b, c, t_cur(2)
   ! Define intersection point:
   ! Coefficients:
   a = vx*vx + vy*vy + vz*vz
   b = x*vx + y*vy + z*vz
   c = x*x + y*y + z*z - R*R
   t = -1.0d0    ! no crossing possible, to start with
   call solution_of_quadratic_equation(a, b, c, t_cur(1), t_cur(2))   ! below
   if (t_cur(1) > 0.0d0) then
      if (t_cur(2) > 0.0d0) then
         t = minval(t_cur(:))
      else
         t = t_cur(1)
      endif
   else
      if (t_cur(2) > 0.0d0) t = t_cur(2)
   endif
end subroutine intersection_with_sphere


pure subroutine solution_of_quadratic_equation(a, b, c, t1, t2)
   real(8), intent(out) :: t1, t2   ! solution of equation: a*t^2 + 2*t*b + c^2 = 0
   real(8), intent(in) :: a, b, c   ! coefficients
   real(8) :: ba, sa, s
   s = b*b - a*c
   if (s < 0.0d0) then  ! no solution
      t1 = -1.0d20
      t2 = -1.0d20
   else
      ba = -b/a
      sa = sqrt(s)/a
      t1 = -ba + sa
      t2 = -ba - sa
   endif
end subroutine solution_of_quadratic_equation


!------------------------------------------------------------------------------------
! Coordinate systems convertions:
 
pure subroutine shift_coordinate_system(X0, Y0, Z0, Xs, Ys, Zs, X, Y, Z)
   real(8), intent(in) :: X0, Y0, Z0	! old coordinates
   real(8), intent(in) :: Xs, Ys, Zs	! vector to shift by
   real(8), intent(out) :: X, Y, Z	! new coordinates
   X = X0 + Xs
   Y = Y0 + Ys
   Z = Z0 + Zs
end subroutine shift_coordinate_system


pure subroutine rotate_coordinate_system(X0, Y0, Z0, X, Y, Z, angle_x, angle_y, angle_z)
   real(8), intent(in) :: X0, Y0, Z0	! old coordinates
   real(8), intent(out) :: X, Y, Z	! new coordinates
   real(8), intent(in), optional :: angle_x, angle_y, angle_z	! roration angles around X, Y, Z axes, [0..2*Pi)
   real(8) :: R, eps
   !eps = 1.0d-6  ! tolerance for considering zngle to be zero
   eps = m_tollerance_eps
   ! To start with:
   X = X0
   Y = Y0
   Z = Z0
   if (present(angle_x)) then
      if (abs(angle_x) > eps) then
         R = sqrt(Y*Y + Z*Z)
         Z = R*cos(angle_x)
         Y = R*sin(angle_x)
      endif
   endif
   if (present(angle_y)) then
      if (abs(angle_y) > eps) then
         R = sqrt(X*X + Z*Z)
         Z = R*cos(angle_y)
         X = R*sin(angle_y)
      endif
   endif
   if (present(angle_z)) then
      if (abs(angle_z) > eps) then
         R = sqrt(Y*Y + X*X)
         X = R*cos(angle_z)
         Y = R*sin(angle_z)
      endif
   endif
end subroutine rotate_coordinate_system


! Since rotation is a non-commutative procidure, when we need to rotate back, we have to change the order of rotations:
pure subroutine rotate_coordinate_system_inverse(X0, Y0, Z0, X, Y, Z, angle_x, angle_y, angle_z)
   real(8), intent(in) :: X0, Y0, Z0	! old coordinates
   real(8), intent(out) :: X, Y, Z	! new coordinates
   real(8), intent(in), optional :: angle_x, angle_y, angle_z	! roration angles around X, Y, Z axes, [0..2*Pi)
   real(8) :: R, eps
   !eps = 1.0d-6  ! tolerance for considering zngle to be zero
   eps = m_tollerance_eps
   ! To start with:
   X = X0
   Y = Y0
   Z = Z0
   if (present(angle_z)) then
      if (abs(angle_z) > eps) then
         R = sqrt(Y*Y + X*X)
         X = R*cos(angle_z)
         Y = R*sin(angle_z)
      endif
   endif
   if (present(angle_y)) then
      if (abs(angle_y) > eps) then
         R = sqrt(X*X + Z*Z)
         Z = R*cos(angle_y)
         X = R*sin(angle_y)
      endif
   endif
   if (present(angle_x)) then
      if (abs(angle_x) > eps) then
         R = sqrt(Y*Y + Z*Z)
         Z = R*cos(angle_x)
         Y = R*sin(angle_x)
      endif
   endif
end subroutine rotate_coordinate_system_inverse


pure function get_v_theta(V, Vasb_in) result(theta)
   real(8) theta    ! [deg] angle with respect to the local Z-axis shifted to the particle position
   real(8), dimension(3), intent(in) :: V   ! particle velosity
   real(8), optional, intent(in) :: Vasb_in ! absolute value of velosity
   real(8) :: Vabs
   if (present(Vasb_in)) then
      Vabs = Vasb_in
   else
      Vabs = sqrt( SUM(V(:)*V(:)) )
   endif
   if (Vabs > m_tollerance_eps) then   ! moving particle
      theta = acos(V(3)/Vabs) * g_rad2deg  ! [deg]
   else ! still particle
      theta = 0.0d0
   endif
end function get_v_theta



pure subroutine Cartesian_to_spherical(X, Y, Z, R, theta, phi)
   real(8), intent(in) :: X, Y, Z		! cartesian coordinates
   real(8), intent(out) :: R, theta, phi	! spherical coordinates
   real(8) :: X2, Y2, XY2
   X2 = X*X
   Y2 = Y*Y
   XY2 = sqrt(X2 + Y2)
   R = sqrt(XY2*XY2 + Z*Z)
   theta = acos(Y/XY2)
   phi = acos(Z/XY2)
end subroutine Cartesian_to_spherical


pure subroutine Spherical_to_cartesian(R, theta, phi, X, Y, Z)
   real(8), intent(in) :: R, theta, phi	! spherical coordinates: theta=[0,Pi] counted from Z; phi=[0,2*Pi] within (X,Y) plane
   real(8), intent(out) :: X, Y, Z		! cartesian coordinates
   real(8) :: sinphi, Rsin
   sinphi = sin(theta)
   Rsin = R*sinphi
   X = Rsin*sin(phi)
   Y = Rsin*cos(phi)
   Z = R*cos(theta)
end subroutine Spherical_to_cartesian


pure subroutine Spherical_to_cartesian_OLD(R, theta, phi, X, Y, Z)      ! Different definition of angles
   real(8), intent(in) :: R, theta, phi	! spherical coordinates: phi=[0,Pi] counted from Z; theta=[0,2*Pi] within (X,Y) plane
   real(8), intent(out) :: X, Y, Z		! cartesian coordinates
   real(8) :: sinphi, Rsin
   sinphi = sin(phi)
   Rsin = R*sinphi
   X = Rsin*sin(theta)
   Y = Rsin*cos(theta)
   Z = R*cos(phi)
end subroutine Spherical_to_cartesian_OLD





pure subroutine Cartesian_to_cylindrical(X, Y, Z, R, L, theta)
   real(8), intent(in) :: X, Y, Z		! cartesian coordinates
   real(8), intent(out) :: R, L, theta		! cylindrical coordinates: R [A], L [A], Phi [rad]
   real(8) :: X2, Y2
   X2 = X*X
   Y2 = Y*Y
   R = sqrt(X2 + Y2)
   L = Z
   theta = acos(Y/R)
end subroutine Cartesian_to_cylindrical


pure subroutine Cylindrical_to_cartesian(R, L, theta, X, Y, Z)
   real(8), intent(in) :: R, L, theta		! cylindrical coordinates
   real(8), intent(out) :: X, Y, Z		! cartesian coordinates
   X = R*sin(theta)
   Y = R*cos(theta)
   Z = L
end subroutine Cylindrical_to_cartesian

end module Geometries
