! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2015-2021
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains useful universal-purpuse subroutines

module Little_subroutines

! Other module:
use Universal_constants, only : g_Pi, g_h, g_e

implicit none

real(8) :: m_Pi_6

parameter (m_Pi_6 = (g_Pi/6.0d0)**(1.0d0/3.0d0))

! this interface finds by itself which of the two subroutine to use depending on the parameters passed:
interface extend_array_size ! extend array size by factor of 2
    module procedure extend_array_size_int
    module procedure extend_array_size_real
end interface extend_array_size

! this interface finds by itself which of the two subroutine to use depending on the parameters passed:
interface extend_array_size_by_one ! extend array size by 1
    module procedure extend_array_size_by_one_int
    module procedure extend_array_size_by_one_real
end interface extend_array_size_by_one

! this interface finds by itself which of the two subroutine to use depending on the dimensions of the array passed:
interface Find_in_array ! search cheaking one by one
    module procedure Find_in_1D_array
    module procedure Find_in_2D_array
end interface Find_in_array
! this interface finds by itself which of the two subroutine to use depending on the dimensions of the array passed:
interface Find_in_array_monoton ! search with bisection method
    module procedure Find_in_monotonous_1D_array
    module procedure Find_in_monotonous_2D_array
end interface Find_in_array_monoton

! This interface allows to sort real or integer arrays:
interface sort_array_bubble
    module procedure sort_array_bubble_real
    module procedure sort_array_bubble_int
end interface sort_array_bubble

! this is a function that returns the order of a passed number:
interface find_order_of_number
   module procedure find_order_of_number_real
   module procedure find_order_of_number_int
end interface find_order_of_number   

! this is a function that creates a spatial grid:
interface create_grid
   module procedure create_grid_1d
   module procedure create_grid_2d
end interface create_grid

! this is a function that set a spatial grid:
interface grid_count
   module procedure grid_count_1d
   module procedure grid_count_2d
end interface grid_count


! private :: ! hides items not listed on public statement 
public :: Find_in_array, Find_in_array_monoton, extend_array_size, find_order_of_number, sort_array_bubble, extend_array_size_by_one

 contains



subroutine solve_rate_equation(X, X0, dt, tau, ind) ! solving dX/dt = -(X - X0)/tau
   real(8), intent(inout) :: X ! variable in the rate equation
   real(8), intent(in) :: X0   ! the equilibrium value that X approaches
   real(8), intent(in) :: dt   ! the time-step
   real(8), intent(in) :: tau  ! the characteristic time of the rate equation
   integer, intent(in) :: ind  ! which scheme to use
   real(8) dttau
   dttau = dt/tau
   select case (ind)
   case (1) ! implicit
      X = (X + dttau*X0)/(1.0d0 + dttau)
   case (2) ! explicit
      X = X - dttau*(X - X0)
   case default ! analytical exact solution
      X = X0 + (X - X0)*dexp(-dttau)
   end select
end subroutine solve_rate_equation
 

pure function Debye_energy(At_Dens, Vsound) result(Edebye)
   real(8) Edebye   ! [eV] maximal phonon energy within Debye approximation
   real(8), intent(in) :: At_Dens   ! [1/cm^3] atomic density
   real(8), intent(in) :: Vsound    ! [m/s] speed of sound in the material
   real(8) :: qdebye
   qdebye = (6.0d0*g_Pi*g_Pi*(At_Dens*1d6))**(0.33333333d0)   ! Debye momentum [1/m]
   Edebye = g_h*Vsound*qdebye/g_e  ! Debye energy [eV]
end function Debye_energy


pure function Einstein_energy(Edebye) result(Eeinstein) ! https://en.wikipedia.org/wiki/Debye_model#Debye_versus_Einstein
   real(8) Eeinstein   ! [eV] maximal phonon energy within Debye approximation
   real(8), intent(in) :: Edebye   ! [eV] maximal phonon energy within Debye approximation
   Eeinstein = Edebye * m_Pi_6  ! Einstein energy [eV]
end function Einstein_energy




! Mark all shells as core or valent:
pure subroutine tick_valence(Ne_per_shell, NVB, is_valence)
   real(8), dimension(:), intent(in) :: Ne_per_shell	! number fo electrons in a shell
   real(8), intent(in) :: NVB	! total number of valence electrons
   logical, dimension(:), allocatable, intent(inout) :: is_valence	! mark whether this shell is valent or not, for each shell
   real(8) :: Ne
   integer :: i, N
   N = size(Ne_per_shell)	! total number of shells
   if (.not.allocated(is_valence)) allocate(is_valence(N))
   is_valence = .false.	! by default, the shell is not valent
   Ne = 0.0d0	! to start counting
   do i = N, 1, -1	! count from top to bottom
      Ne = Ne + Ne_per_shell(i)
      if (Ne <= NVB) then	! it is a valent shell
         is_valence(i) = .true.
      else	! it is a deep shell, as well as all deeper shells
         exit	! so no need to check them
      endif
   enddo
end subroutine tick_valence
 

pure subroutine create_energy_grid(x_start, x_end, array1d, special_points)
   real(8), intent(in) :: x_start, x_end	! starting and ending points of the grid
   real(8), dimension(:), allocatable, intent(inout) :: array1d	! array into which the grid will be saved
   real(8), dimension(:), intent(in), optional :: special_points	! make finer grid around those points
   !------------------------------------------
   real(8) :: deg0, deg1, x_0, x_1
   integer :: i, j, siz
   deg0 = find_order_of_number(x_start)	! function below
   x_0 = 10.0d0**(deg0-1)	! starting point of the grid
   deg1 = find_order_of_number(x_end) - 1	! function below
   x_1 = 10.0d0**deg1	! ending point of the grid
   ! Count how many grid points will be in such a grid:
   if (present(special_points)) then	! add spectial points to the grid
      call energy_grid_count(x_0, x_1, deg0, deg1, siz=siz, special_points=special_points)	! below
   else ! no special points, smooth grid
      call energy_grid_count(x_0, x_1, deg0, deg1, siz=siz)	! below
   endif
   ! Allocate array correspondingly:
   allocate(array1d(siz))
   ! Fill array with grid points:
   if (present(special_points)) then	! add spectial points to the grid
      call energy_grid_count(x_0, x_1, deg0, deg1, array1d=array1d, special_points=special_points)	! below
   else ! no special points, smooth grid
      call energy_grid_count(x_0, x_1, deg0, deg1, array1d=array1d)	! below
   endif
end subroutine create_energy_grid



pure subroutine energy_grid_count(x_0, x_1, deg0, deg1, array1d, siz, special_points)
   real(8), intent(in) :: deg0, deg1, x_0, x_1	! degrees and numbers of the starting and ending grid point
   real(8), dimension(:), intent(out), optional :: array1d	! array into which the grid will be saved
   integer, intent(out), optional :: siz	! calculate only the grid size
   real(8), dimension(:), intent(in), optional :: special_points	! make finer grid around those points
   integer :: i, j, k, i_sp, N_sp
   real(8) :: x_cur, dx, x_sp_eps
   integer, dimension(:), allocatable :: array_sp
   real(8), dimension(:), allocatable :: copy_special_points
   
   if (present(special_points)) then
      x_sp_eps = 1.0d-3	! how close a grid point should be around a special point
      allocate(copy_special_points(size(special_points)))
      copy_special_points = special_points	! to operate with
      ! Make sure all special points are in accending order:
      call sort_array_bubble(copy_special_points)	! below
   endif
   
   k = 1
   x_cur = x_0
   i = deg0
   if (present(array1d)) array1d(k) = x_cur
   MN:do while (x_cur < x_1)
      k = k + 1
      if (x_cur < 0.01d0) then
         dx = 10.0d0**(find_order_of_number(x_cur)-2)	! function below
      else if (x_cur < 1.0d0) then
         dx = 0.1
      else if (x_cur < 100.0d0) then
         dx = 1.0d0
      else
         dx = 10.0d0**(find_order_of_number(x_cur)-2)	! function below
      endif
      
      ! Check if we need a finer grid around special points:
      if (present(special_points)) then
         call check_special_point(x_cur, dx, copy_special_points, N_sp, array_sp)	! see below
         if (N_sp > 0) then	! there are special points within this interval:
            do i_sp = 1, N_sp	! for all special points, make a grid around them
               if (i_sp > 1) then ! make sure the point does not repeat - there is no need to refine grid around the same point twice
                  if (ABS(copy_special_points(array_sp(i_sp))-copy_special_points(array_sp(i_sp-1))) > x_sp_eps ) then
!                      print*, 'SP',   N_sp, i_sp, array_sp(i_sp) ,copy_special_points(array_sp(i_sp)) 
                     if (present(array1d)) array1d(k) =  copy_special_points(array_sp(i_sp)) - x_sp_eps
                     k = k + 1
                     if (present(array1d)) array1d(k) =  copy_special_points(array_sp(i_sp)) + x_sp_eps
                     k = k + 1
                  endif
               else ! it's the first point
                  ! If special point councides with a grid point, exlude it:
                  if (ABS(copy_special_points(array_sp(i_sp))-x_cur) > x_sp_eps ) then ! it does not coincide
!                      print*, 'SP',   N_sp, i_sp, array_sp(i_sp) ,copy_special_points(array_sp(i_sp)) 
                     if (present(array1d)) array1d(k) =  copy_special_points(array_sp(i_sp)) - x_sp_eps
                     k = k + 1
                     if (present(array1d)) array1d(k) =  copy_special_points(array_sp(i_sp)) + x_sp_eps
                     k = k + 1
                  else ! it coincides
!                      print*, 'SP',   N_sp, i_sp, array_sp(i_sp) ,copy_special_points(array_sp(i_sp)) 
                     k = k - 1 ! to overwrite the point
                     if (present(array1d)) array1d(k) =  copy_special_points(array_sp(i_sp)) - x_sp_eps
                     k = k + 1
                     if (present(array1d)) array1d(k) =  copy_special_points(array_sp(i_sp)) + x_sp_eps
                     k = k + 1
                  endif
               endif ! (i_sp > 1)
            enddo ! i_sp = 1, N_sp
         endif ! (N_sp > 0) 
      endif ! (present(special_points))

      x_cur = x_cur + dx
      if (present(array1d)) array1d(k) = x_cur
   enddo MN
   if (present(siz)) siz = k
   if (allocated(copy_special_points)) deallocate(copy_special_points)
end subroutine energy_grid_count


pure subroutine check_special_point(x_cur, dx, special_points, N_sp, array_sp)
   real(8), intent(in) :: x_cur, dx
   real(8), dimension(:), intent(in) :: special_points	! make finer grid around those points
   integer, intent(out) :: N_sp	! how many special points are inside of this interval: [x_cur: x_cur+dx]
   integer, dimension(:), allocatable, intent(out) :: array_sp	! locations of special points within the given interval
   integer :: i, j
   real(8) :: x_end
   x_end = x_cur+dx
   ! Count how many special points are within this interval:
   N_sp = count( ((special_points(:) >= x_cur) .and. (special_points(:) < x_end)) )
   if (N_sp > 0) then
      ! Alloate new array for the locations of special points:
      if (allocated(array_sp)) deallocate(array_sp)
      allocate(array_sp(N_sp))
      array_sp = 0
      ! Save locations of special points into this array:
      j = 0	! to start with
      do i = 1, size(special_points)
         if ((special_points(i) >= x_cur) .and. (special_points(i) < x_end)) then
            j = j + 1
            array_sp(j) = i
         endif
      enddo
      if (N_sp > 1) then	! make sure they are in accending order:
         call sort_array_bubble(array_sp)	! below
      endif
   endif
end subroutine check_special_point


 subroutine create_grid_1d(x_start, x_end, dx, array1d, ind)
   real(8), intent(in) :: x_start, x_end, dx	! starting and ending points of the grid, and its step
   real(8), dimension(:), allocatable, intent(inout) :: array1d	! array into which the grid will be saved
   integer, intent(in) :: ind	! kind of grid: 0=linear, 1=log
   !------------------------------------------
   real(8) :: deg0, deg1, x_0, x_1
   integer :: i, j, siz
   
!    print*, 'create_grid_1d', x_start
   
   if (ind == 1) then   ! log scale
      deg0 = find_order_of_number(x_start)	! function below
      x_0 = 10.0d0**(deg0-1)	! starting point of the grid
      x_0 = SIGN(x_0, x_start)
   else ! linear scale
      x_0 = FLOOR(x_start)
   endif
!    print*, 'create_grid_1d', deg0, x_0
   
   if (ind == 1) then   ! log scale
      deg1 = find_order_of_number(x_end-1.0d-6)   ! function below
      x_1 = 10.0d0**deg1	! ending point of the grid
      x_1 = SIGN(x_1, x_end)
   else
      x_1 = CEILING(x_end)
   endif
   
   ! Count how many grid points will be in such a grid:
   call grid_count(x_0, x_1, deg0, deg1, dx, ind, siz=siz)	! below
   ! Allocate array correspondingly:
   allocate(array1d(siz))
   ! Fill array with grid points:
   call grid_count(x_0, x_1, deg0, deg1, dx, ind, array1d=array1d)	! below
end subroutine create_grid_1d


 subroutine create_grid_2d(x_start, x_end, dx, array2d, ind, dim_siz, dim_ind)
   real(8), intent(in) :: x_start, x_end, dx	! starting and ending points of the grid, and its step
   real(8), dimension(:,:), allocatable, intent(inout) :: array2d	! array into which the grid will be saved
   integer, intent(in) :: ind          ! kind of grid: 0=linear, 1=log
   integer, intent(in) :: dim_siz   ! size of the first dimension of the array
   integer, intent(in) :: dim_ind   ! which dimension to set
   !------------------------------------------
   real(8) :: deg0, deg1, x_0, x_1
   integer :: i, j, siz
   
!    print*, 'create_grid_2d', x_start
   
   deg0 = find_order_of_number(x_start)	! function below
   x_0 = 10.0d0**(deg0-1)	! starting point of the grid
   
!    print*, 'create_grid_2d', deg0, x_0
   
!    deg1 = find_order_of_number(x_end) - 1	! function below
   deg1 = find_order_of_number(x_end)   ! function below
   x_1 = 10.0d0**deg1	! ending point of the grid
   ! Count how many grid points will be in such a grid:
   call grid_count(x_0, x_1, deg0, deg1, dx, ind, siz=siz)	! below
   ! Allocate array correspondingly:
   allocate(array2d(dim_siz,siz))
   ! Fill array with grid points:
   call grid_count(x_0, x_1, deg0, deg1, dx, ind, dim_ind, array2d=array2d)	! below
end subroutine create_grid_2d



pure subroutine grid_count_1d(x_0, x_1, deg0, deg1, dx, ind, array1d, siz)
   real(8), intent(in) :: deg0, deg1, x_0, x_1, dx	! degrees and numbers of the starting and ending grid point and its step
   integer, intent(in) :: ind	! kind of the grid
   real(8), dimension(:), intent(out), optional :: array1d	! array into which the grid will be saved
   integer, intent(out), optional :: siz
   integer :: i, j, k
   real(8) :: x_cur, dx_add
   
   k = 1
   x_cur = x_0
   select case (ind)	! kind of grid
   case (1)	! log
      i = deg0
   case default
      i = 0
   endselect
   
!    print*, x_0, x_1
   MN:do while (x_cur < x_1)
      select case (ind)	! kind of grid
      case default	! linear
         x_cur = x_0
         if (present(array1d)) array1d(k) = x_cur
         do while (x_cur < x_1)
            k = k +1
            x_cur = x_cur + dx
            if (present(array1d)) array1d(k) = x_cur
            if (x_cur >= x_1) exit MN
         enddo
      case (1)	! log
         dx_add = 10.0d0**(i-1) * dx
         do j = 1, CEILING(9.0/dx)
            if (present(array1d)) array1d(k) = x_cur
            if (x_cur >= x_1) exit MN
            k = k +1
            x_cur = x_cur + dx_add
         enddo  !  j = 1, 9
      endselect
      i = i + 1	! increase the order by one
   enddo MN
   if (present(array1d) .and. (ind == 1)) array1d(k) = x_cur
   if (present(siz)) siz = k
end subroutine grid_count_1d


 subroutine grid_count_2d(x_0, x_1, deg0, deg1, dx, ind, dim_ind, array2d, siz)
   real(8), intent(in) :: deg0, deg1, x_0, x_1, dx	! degrees and numbers of the starting and ending grid point and its step
   integer, intent(in) :: ind	! kind of the grid
   integer, intent(in) :: dim_ind   ! dimension to fill
   real(8), dimension(:,:), intent(out), optional :: array2d	! array into which the grid will be saved
   integer, intent(out), optional :: siz
   integer :: i, j, k
   real(8) :: x_cur, dx_add
   
   k = 1
   x_cur = x_0
   i = deg0
!     print*, 'grid_count_2d', x_0, x_1
   MN:do while (x_cur < x_1)
      select case (ind)	! kind of grid
      case default	! linear
         x_cur = x_0
         if (present(array2d)) array2d(dim_ind,k) = x_cur
         do while (x_cur < x_1)
            k = k +1
            x_cur = x_cur + dx
            if (present(array2d)) array2d(dim_ind,k) = x_cur
            if (x_cur >= x_1) exit MN
         enddo
      case (1)	! log
          dx_add = 10.0d0**(i-1) * dx
         do j = 1, CEILING(9.0/dx)
            if (present(array2d)) array2d(dim_ind,k) = x_cur
            if (x_cur >= x_1) exit MN
            k = k +1
            x_cur = x_cur + dx_add
         enddo  !  j = 1, 9
      endselect
      i = i + 1	! increase the order by one
   enddo MN
   if (present(array2d) .and. (ind == 1)) array2d(dim_ind,k) = x_cur
   if (present(siz)) siz = k
end subroutine grid_count_2d


pure subroutine grid_count_1d_OLD(x_0, x_1, deg0, deg1, dx, ind, array1d, siz)
   real(8), intent(in) :: deg0, deg1, x_0, x_1, dx	! degrees and numbers of the starting and ending grid point and its step
   integer, intent(in) :: ind	! kind of the grid
   real(8), dimension(:), intent(out), optional :: array1d	! array into which the grid will be saved
   integer, intent(out), optional :: siz
   integer :: i, j, k
   real(8) :: x_cur
   
   k = 1
   x_cur = x_0
   i = deg0
!    print*, x_0, x_1
   MN:do while (x_cur < x_1)
      select case (ind)	! kind of grid
      case default	! linear
         x_cur = x_0
         if (present(array1d)) array1d(k) = x_cur
         do while (x_cur < x_1)
            k = k +1
            x_cur = x_cur + dx
            if (present(array1d)) array1d(k) = x_cur
            if (x_cur >= x_1) exit MN
         enddo
      case (1)	! log
         do j = 1, 9
!             x_cur = (dble(j)*dx)*10.0d0**(i-1)
            x_cur = (dble(j))*10.0d0**(i-1) * dx
            if (present(array1d)) array1d(k) = x_cur
            if (x_cur >= x_1) exit MN
            if (j == 1) then
               k = k +1
!                x_cur = (1.5d0*dx)*10.0d0**(i-1)
               x_cur = (1.5d0)*10.0d0**(i-1) * dx
               if (present(array1d)) array1d(k) = x_cur
               if (x_cur >= x_1) exit MN
            endif
            k = k +1
         enddo  !  j = 1, 9
      endselect
      i = i + 1	! increase the order by one
   enddo MN
   if (present(siz)) siz = k
end subroutine grid_count_1d_OLD


 subroutine grid_count_2d_OLD(x_0, x_1, deg0, deg1, dx, ind, dim_ind, array2d, siz)
   real(8), intent(in) :: deg0, deg1, x_0, x_1, dx	! degrees and numbers of the starting and ending grid point and its step
   integer, intent(in) :: ind	! kind of the grid
   integer, intent(in) :: dim_ind   ! dimension to fill
   real(8), dimension(:,:), intent(out), optional :: array2d	! array into which the grid will be saved
   integer, intent(out), optional :: siz
   integer :: i, j, k
   real(8) :: x_cur
   
   k = 1
   x_cur = x_0
   i = deg0
!     print*, 'grid_count_2d', x_0, x_1
   MN:do while (x_cur < x_1)
      select case (ind)	! kind of grid
      case default	! linear
         x_cur = x_0
         if (present(array2d)) array2d(dim_ind,k) = x_cur
         do while (x_cur < x_1)
            k = k +1
            x_cur = x_cur + dx
            if (present(array2d)) array2d(dim_ind,k) = x_cur
            if (x_cur >= x_1) exit MN
         enddo
      case (1)	! log
         do j = 1, 9
!             x_cur = (dble(j)*dx)*10.0d0**(i-1)
            x_cur = (dble(j))*10.0d0**(i-1) * dx
            if (present(array2d)) array2d(dim_ind,k) = x_cur
            if (x_cur >= x_1) exit MN
            if (j == 1) then
               k = k +1
!                x_cur = (1.5d0*dx)*10.0d0**(i-1)
               x_cur = (1.5d0)*10.0d0**(i-1) * dx
               if (present(array2d)) array2d(dim_ind,k) = x_cur
               if (x_cur >= x_1) exit MN
            endif
            k = k +1
         enddo  !  j = 1, 9
      endselect
      i = i + 1	! increase the order by one
   enddo MN
   if (present(siz)) siz = k
end subroutine grid_count_2d_OLD


pure subroutine order_of_time(tim, text, gnu_text, x_tics)
   real(8), intent(in) :: tim ! time to find its order
   character(*), intent(out) :: text ! fs, ps, ns, mks, ms, s
   character(*), intent(out), optional :: gnu_text ! culomn to set in gnuplot
   real(8), intent(out), optional :: x_tics ! tics for gnuplot
   integer :: time_ord
   time_ord = find_order_of_number(tim) ! module "Little_subroutines"
   if (present(x_tics)) then
      x_tics = 10**(time_ord) ! set tics for gnuplot
      if (tim/dble(x_tics) > 0.5) then
         x_tics = 10**(time_ord-1) ! set tics for gnuplot
      else if (tim/dble(x_tics) > 0.2) then
         x_tics = 0.5*10**(time_ord-1) ! set tics for gnuplot
      else
         x_tics = 10**(time_ord-2) ! set tics for gnuplot
      endif
   endif

   if (time_ord > 1e15) then ! s
      text = '(s)'
      if (present(gnu_text)) gnu_text = '($1/1e15)'
   else if (time_ord > 1e12) then ! ms
      text = '(ms)'
      if (present(gnu_text)) gnu_text = '($1/1e12)'
   else if (time_ord > 1e9) then ! mks
      text = '(mks)'
      if (present(gnu_text)) gnu_text = '($1/1e9)'
   else if (time_ord > 1e6) then ! ns
      text = '(ns)'
      if (present(gnu_text)) gnu_text = '($1/1e6)'
   else if (time_ord > 1e3) then ! ps
      text = '(ps)'
      if (present(gnu_text)) gnu_text = '($1/1e3)'
   else ! fs
      text = '(fs)'
      if (present(gnu_text)) gnu_text = '($1)'
   endif
end subroutine order_of_time



pure subroutine order_of_energy(E, text, gnu_text, x_tics, E_reduced)
   real(8), intent(in) :: E ! energy to find its order
   character(*), intent(out) :: text ! meV, eV, keV, MeV, TeV, GeV
   character(*), intent(out), optional :: gnu_text ! culomn to set in gnuplot
   real(8), intent(out), optional :: x_tics ! tics for gnuplot
   real(8), intent(out), optional :: E_reduced  ! energy in the found units
   integer :: time_ord

   if (present(x_tics)) then
      time_ord = find_order_of_number(E) ! module "Little_subroutines"
      x_tics = 10**(time_ord) ! set tics for gnuplot
      if (E/dble(x_tics) > 0.5) then
         x_tics = 10**(time_ord-1) ! set tics for gnuplot
      else if (E/dble(x_tics) > 0.2) then
         x_tics = 0.5*10**(time_ord-1) ! set tics for gnuplot
      else
         x_tics = 10**(time_ord-2) ! set tics for gnuplot
      endif
   endif

   if (E >= 1e12) then ! TeV
      text = 'TeV'
      if (present(E_reduced)) E_reduced = E/1.0d12
      if (present(gnu_text)) gnu_text = '($1/1e12)'
   else if (E >= 1e9) then ! GeV
      text = 'GeV'
      if (present(E_reduced)) E_reduced = E/1.0d9
      if (present(gnu_text)) gnu_text = '($1/1e9)'
   else if (E >= 1e6) then ! MeV
      text = 'MeV'
      if (present(E_reduced)) E_reduced = E/1.0d6
      if (present(gnu_text)) gnu_text = '($1/1e6)'
   else if (E >= 1e3) then ! keV
      text = 'keV'
      if (present(E_reduced)) E_reduced = E/1.0d3
      if (present(gnu_text)) gnu_text = '($1/1e3)'
   else if (E >= 1.0d0) then ! eV
      text = 'eV'
      if (present(E_reduced)) E_reduced = E
      if (present(gnu_text)) gnu_text = '($1)'
   else if (E >= 1.0d-3) then ! meV
      text = 'meV'
      if (present(E_reduced)) E_reduced = E/1.0d-3
      if (present(gnu_text)) gnu_text = '($1*1e3)'
   else ! default
      text = 'eV'
      if (present(E_reduced)) E_reduced = E
      if (present(gnu_text)) gnu_text = '($1)'
   endif
end subroutine order_of_energy



pure function find_order_of_number_real_OLD(num) result(find_order_of_number_real)
   integer find_order_of_number_real
   real(8), intent(in) :: num
   character(64) :: temp, temp2
   real(8) num0
   integer i, counter
   num0 = ABS(num)
   if (num0 >= 1) then
      write(temp,'(e)') num0	! make it a string
      do i = 1,LEN(TRIM(adjustl(temp)))
         if (TRIM(adjustl(temp(i:i))) == 'E') then
            counter = i+1
            exit
         endif
      enddo
      temp2 = TRIM(adjustl(temp(counter:)))
      read(temp2,*)  num0
      find_order_of_number_real = num0
   else
      counter = -1
      write(temp,'(f)') num0	! make it a string
      do i = 3,LEN(TRIM(adjustl(temp)))
         if (TRIM(adjustl(temp(i:i))) /= '0') exit
         counter = counter -1
      enddo
      find_order_of_number_real = counter
   endif
end function find_order_of_number_real_OLD


pure function find_order_of_number_real(num)
   integer find_order_of_number_real
   real(8), intent(in) :: num
   character(64) :: temp, temp2
   real(8) num0
   integer i, counter, i_first
   num0 = ABS(num)
   write(temp,'(e)') num0	! make it a string
   do i = 1,LEN(TRIM(adjustl(temp)))
      if (TRIM(adjustl(temp(i:i))) == 'E') then
         counter = i+1
         exit
      endif
   enddo
   temp2 = TRIM(adjustl(temp(counter:)))
   read(temp2,*)  num0
   find_order_of_number_real = num0
end function find_order_of_number_real

pure function find_order_of_number_int(num)
   integer find_order_of_number_int
   integer, intent(in) :: num
   character(64) :: temp
!    integer :: i_first
   write(temp,'(i)') num ! make it a string
   find_order_of_number_int = LEN(TRIM(adjustl(temp))) ! find how many characters in this string
end function find_order_of_number_int


pure subroutine parse_time(sec,chtest)
   real(8), intent(inout) :: sec ! time interval in [sec]
   character(*), intent(out) :: chtest ! split it into mituns, hours, days...
   character(100) temp
   real(8) days, hours, mins
   days = 0.0d0
   hours = 0.0d0
   mins = 0.0d0
   if (sec .GE. 60.0d0) then
      mins = FLOOR(sec/60.0d0)
      sec = sec - mins*60.0d0
      if (mins .GT. 60.0d0) then
         hours = FLOOR(mins/60.0d0)
         mins = mins - hours*60.0d0
         if (hours .GT. 24.0d0) then
            days = FLOOR(hours/24.0d0)
            hours = hours - days*24.0d0
         endif
      endif
   endif
   chtest = ''
   temp = ''
   if (days .GT. 1.0d0) then
      write(temp, '(i9)') int(days)
      write(chtest, '(a,a)') trim(adjustl(temp)), ' days'
   else if (days .GT. 0.0d0) then
      write(temp, '(i9)') int(days)
      write(chtest, '(a,a)') trim(adjustl(temp)), ' day'
   endif
   if (hours .GT. 1.0d0) then
      write(temp, '(i9)') int(hours)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' hours'
   else if (hours .GT. 0.0d0) then
      write(temp, '(i9)') int(hours)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' hour'
   endif
   if (mins .GT. 1.0d0) then
      write(temp, '(i9)') int(mins)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' mins'
   else if (mins .GT. 0.0d0) then
      write(temp, '(i9)') int(mins)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' min'
   endif
   write(temp, '(f7.3)') sec
   write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' sec'
end subroutine parse_time


subroutine print_progress(string,ndone,ntotal)
    implicit none
    character*(*) string
    character(255) prog,oldprog
    integer ndone,ntotal,i
    if (100.0*ndone/ntotal .GE. 100.0d0) then
        write(0,'(a,$)') '                                                                                   ',char(13)
    else
        write(prog,'(a25,1x,''['')') string
        do i=1,40
            prog(27+i:27+i)=' '
        enddo
        write(prog(44:51),'(f7.2,''%'')') 100.0*ndone/ntotal
        do i=1,40
            if ((1.0*ndone/ntotal).gt.(1.0*i/40)) then
                if (prog(27+i:27+i).eq.' ') prog(27+i:27+i)='-'
            endif
        enddo
        prog(67:67)=']'
        write(0,'(a,a,$)') prog(1:77),char(13)
        return
    endif
end subroutine print_progress



function sample_Poisson(l0) result(l_sampled)
   real(8) l_sampled    ! Sample point according to Poisson distribution
   real(8), intent(in) :: l0    ! mean value to sample around (e.g. MFP, decay time, etc.)
   real(8) :: RN
   call random_number(RN)
   l_sampled = -l0 * log(RN)
end function sample_Poisson



function sample_Gaussian(x0, FWHM) result(x)
   real(8) :: x ! sampled point according to Gaussian distribution with given parameters
   real(8), intent(in) :: x0, FWHM  ! central point and full width at half maximum of a gaussian
   real(8) :: sigma, R1, R2
   if (FWHM < 1.0d-6) then  ! assume delta function
      x = x0
   else ! sample gaussian
      sigma = FWHM/2.35482d0  ! gaussian sigma parameter
      call random_number(R1)
      call random_number(R2)
      x = x0 + sigma*dsqrt(-2.0d0*dlog(R1))*dcos(2.0d0*g_Pi*R2)    ! Box-Muller method of sampling a gaussian function
   endif
end function sample_Gaussian


pure subroutine Gaussian(mu, sigma, x, Gaus, normalized_max) ! at the time x according to Gaussian shape
   real(8), intent(in) :: mu, sigma, x
   real(8), intent(out) :: Gaus
   real(8), parameter :: Pi = 3.1415926535897932384626433832795d0
   real(8), intent(in), optional :: normalized_max
   if (present(normalized_max)) then
      Gaus = exp(-(x-mu)*(x-mu)/(2.0d0*sigma*sigma))*normalized_max ! it will be equal to "normalized_max" at the maximum
   else ! normalized as integral(Gaus) = 1
      Gaus = 1.0d0/(sqrt(2.0d0*Pi)*sigma)*exp(-(x-mu)*(x-mu)/(2.0d0*sigma*sigma)) ! it will be normalized to integral=1
   endif
end subroutine Gaussian


pure function Rect_pulse(mu, sigma, x) ! number of particles at timestep x, according to flat-top pulse
   real(8), intent(in) :: mu, sigma, x
   real(8) Rect_pulse
   real(8) Gaus_sigma
   Gaus_sigma = sigma !*2.35482d0
   if ((x .GT. (mu - Gaus_sigma/2.0d0)) .AND. (x .LE. (mu + Gaus_sigma/2.0d0)) ) then
      Rect_pulse = 1.0d0/Gaus_sigma
   else
      Rect_pulse = 0.0d0
   endif
end function Rect_pulse

pure function Gaus_pulse(mu, sigma, x) ! number of particles at the time x according to Gaussian shape
   real(8), intent(in) :: mu, sigma, x
   real(8) Gaus_pulse
   if ((x .GT. (mu - 2.5d0*sigma)) .AND. (x .LE. (mu + 2.5d0*sigma)) ) then
      Gaus_pulse = 1.0d0/(sqrt(2.0d0*g_Pi)*sigma)*dexp(-(x-mu)*(x-mu)/(2.0d0*sigma*sigma))
   else
      Gaus_pulse = 0.0d0
   endif
end function Gaus_pulse

pure function SASE_pulse(mu, sigma, x) ! number of particles at the time x according to schematic SASE shape
   real(8), intent(in) :: mu, sigma, x
   real(8) SASE_pulse
   integer i
   real(8) RN(4), Co(4), y, f, SASE
   if ((x .GT. (mu - sigma)) .AND. (x .LE. (mu + sigma)) ) then
      RN(1) = 0.4563449303d0
      RN(2) = 0.1271999433d0
      RN(3) = 1.0d0 - RN(1)
      RN(4) = 1.0d0 - RN(2)
      Co(1) = 2.0d0
      Co(2) = 4.0d0
      Co(3) = 5.0d0
      Co(4) = 9.0d0
      SASE = 0.0d0
      do i = 1,4
         y = Co(i)*g_Pi*(x+sigma-mu)/sigma
         f = sin(y)
         SASE = SASE + f*f/sigma*RN(i)
      enddo
      SASE_pulse = SASE/2.0d0 ! normalization to 1
   else
      SASE_pulse = 0.0d0
   endif
end function SASE_pulse



subroutine interpolate_data_single(given_array_x, given_array_y, x_value_in, y_value_out, x0, y0)
   real(8), dimension(:), intent(in) :: given_array_x, given_array_y
   real(8), intent(in) :: x_value_in    ! value of x, corresponging to which the value of y should be interpolated
   real(8), intent(out) :: y_value_out  ! interpolated value
   real(8), intent(in), optional :: x0, y0 ! assume value for extrapolation outside of the array
   integer :: i_closest
   y_value_out = 0.0d0  ! initiate to start
   call Find_in_array_monoton(given_array_x,  x_value_in, i_closest)	! see below   
   if (present(x0) .and. present(y0)) then
      call linear_interpolation(given_array_x, given_array_y, x_value_in, y_value_out, i_closest, x0, y0)	! see below
   elseif (present(x0)) then
      call linear_interpolation(given_array_x, given_array_y, x_value_in, y_value_out, i_closest, x0)	! see below
   elseif (present(y0)) then
      call linear_interpolation(given_array_x, given_array_y, x_value_in, y_value_out, i_closest, y0)	! see below
   else
      call linear_interpolation(given_array_x, given_array_y, x_value_in, y_value_out, i_closest)	! see below
   endif
end subroutine interpolate_data_single



subroutine interpolate_data_on_grid(given_array_x, given_array_y, given_grid, array_to_fill)
   real(8), dimension(:), intent(in) :: given_array_x, given_array_y, given_grid
   real(8), dimension(:), intent(out) :: array_to_fill
   !----------------------------------------------
   integer :: i, N, i_closest, i_found
   real(8) :: eps
   logical :: repeated
   repeated = .false.
   eps = 1.0d-7
   !N = min(size(given_grid), size(array_to_fill), size(given_array_x))
   N = min(size(given_grid), size(array_to_fill))

!    print*, 'interpolate_data_on_grid', N, size(given_array_x)

   do i = 1, N
      call Find_in_array_monoton(given_array_x,  given_grid(i), i_closest)	! see below
      if (given_grid(i) >= given_array_x(i_closest)) then
         if (i_closest > 1) then
            do while (ABS(given_array_x(i_closest) - given_array_x(i_closest-1)) < eps) ! exclude the same value
               i_closest = i_closest + 1
               !if (i_closest == N) exit
               if (i_closest == size(given_array_x)) exit
               repeated = .true.
            enddo
         endif
         if (repeated .or. (i_closest >= N)) then
            repeated = .false.
            call linear_interpolation(given_array_x, given_array_y, given_grid(i), array_to_fill(i), i_closest) ! see below
         else
            call linear_interpolation(given_array_x, given_array_y, given_grid(i), array_to_fill(i), i_closest+1)	! see below
         endif
!             if ((i_closest > 1) .and. (i_closest < size(given_array_x) )) then
!              write(*,'(a,i4,i5,es,es,es,es)') 'a1', i, i_closest+1, given_grid(i), given_array_x(i_closest-1), given_array_x(i_closest), given_array_x(i_closest+1)
!              write(*,'(a,i4,i5,es,es,es,es)') 'a2', i, i_closest+1, array_to_fill(i), given_array_y(i_closest-1), given_array_y(i_closest), given_array_y(i_closest+1)
!             endif 
      else ! no photoabsorption below the edge:
!          print*, i, i_closest, array_to_fill(i), array_to_fill(i_closest)
         array_to_fill(i) = 0.0d0
      endif
   enddo
!    pause "interpolate_data_on_grid"
end subroutine interpolate_data_on_grid


subroutine linear_interpolation(xarray, yarray, x, y, i, x0, y0, replac)
   real(8), dimension(:), intent(in) :: xarray, yarray	! x-array, y-array
   real(8), intent(in) :: x	! input
   real(8), intent(out) :: y	! output
   integer, intent(in) :: i	! index for x-array
   real(8), intent(in), optional :: x0, y0 ! assume initial value for i = 1
   logical, intent(in), optional :: replac ! replace x0 by the given one independantly on which i is it?
   real(8) eps, eps_inv, dx
   eps = 1.0d10 ! maximal allowed change from y(i-1) to y(i) in case y(i) << y(i-1)
   eps_inv = 1.0d-10
   REDO: if (.not.present(replac)) then
    if (i .GT. 1) then
       dx = xarray(i) - xarray(i-1)
       if (dx < eps_inv) then    ! the values are too close, no need to interpolate
          y = yarray(i)
       else
         if (x - xarray(i-1) .GE. 0.0d0) then
            if ((yarray(i-1) > eps*yarray(i)) .and. (ABS(yarray(i)) > eps_inv) ) then  !  in case y(i) << y(i-1)
               !y = yarray(i) - yarray(i)*(1.0d0-eps)/(xarray(i) - xarray(i-1))*(x - xarray(i-1))
               y = yarray(i) - yarray(i)*(1.0d0-eps)/(dx)*(xarray(i) - x)
!                print*, 'TEST linear_interpolation',  yarray(i),  yarray(i-1)
            else
!                print*, i, xarray(i), xarray(i-1)            
               y = yarray(i-1) + (yarray(i) - yarray(i-1))/(dx)*(x - xarray(i-1))
            endif  ! (yarray(i-1) > eps*yarray(i))
         else !  (x - xarray(i-1) .GE. 0.0d0) 
            if (present(y0) .and. present(x0)) then
               y = y0 + (yarray(i) - y0)/(xarray(i) - x0)*(x - x0)
            else
               if (present(x0)) then
                  y = (yarray(i) - 0.0)/(xarray(i) - x0)*(x - x0)
               else
                  y = (yarray(i) - 0.0)/(xarray(i) - 0.0)*(x - 0.0)
               endif
            endif
         
         endif !  (x - xarray(i-1) .GE. 0.0d0) 
       endif ! (dx < eps_inv)
     else ! (i .GT. 1)
      if (present(y0) .and. present(x0)) then
         if (abs(xarray(i) - x0) < 1.0d-15) then
            y = y0 + (yarray(i) - y0)/(xarray(i+1) - xarray(i)) * (x - x0)
         else
            y = y0 + (yarray(i) - y0)/(xarray(i) - x0) * (x - x0)
         endif
      else
         if (present(x0)) then
            y = (yarray(i) - 0.0)/(xarray(i) - x0)*(x - x0)
         else
!             y = (yarray(i) - 0)/(xarray(i) - 0)*(x - 0)
            y = yarray(i) + (yarray(i+1) - yarray(i))/(xarray(i+1) - xarray(i))*(x - xarray(i))            
         endif
      endif
    endif
   else REDO ! in this case use just the given values:
    if ((replac) .AND. present(x0) .AND. present(y0)) then
       y = y0 + (yarray(i) - y0)/(xarray(i) - x0)*(x - x0)
    else
      if (i .GT. 1) then
         if (x - xarray(i-1) .GE. 0.0d0) then
            y = yarray(i-1) + (yarray(i) - yarray(i-1))/(xarray(i) - xarray(i-1))*(x - xarray(i-1))
         else
            if (present(y0) .and. present(x0)) then
               y = y0 + (yarray(i) - y0)/(xarray(i) - x0)*(x - x0)
            else
               if (present(x0)) then
                  y = (yarray(i) - 0.0)/(xarray(i) - x0)*(x - x0)
               else
                  y = (yarray(i) - 0.0)/(xarray(i) - 0.0)*(x - 0.0)
               endif
            endif
         endif
      endif
    endif
   endif REDO
end subroutine linear_interpolation



 subroutine Fermi_interpolation(xarray, yarray, x, y, i)
   real(8), dimension(:), intent(in) :: xarray, yarray	! x-array, y-array
   real(8), intent(in) :: x	! input
   real(8), intent(out) :: y	! output
   integer, intent(in) :: i	! index for x-array
   real(8) :: mu, T ! parameters of the fermi function to be found and used
   real(8) :: temp
   if (i > 1) then
      if ((yarray(i) == yarray(i-1)) .or. (xarray(i) == xarray(i-1))) then ! no need for doing anything
         y = yarray(i)
      else
         if ((yarray(i-1) < tiny(x)) .or. (yarray(i) < tiny(x)) .or. (yarray(i-1) == 2.0d0) .or. (yarray(i) == 2.0d0)) then ! just in case it's 0
            call linear_interpolation(xarray, yarray, x, y, i)
         else
            temp = log(2.0d0/yarray(i)-1.0d0)
            T = (temp - log(2.0d0/yarray(i-1)-1.0d0))/(xarray(i)-xarray(i-1))
            mu = xarray(i) - T*temp
            y = 2.0d0/(1.0d0 + exp((x - mu)/T))
         endif
      endif
   else ! use default value
      y = yarray(1)
   endif   
end subroutine Fermi_interpolation


pure function Fermi_function(mu, Te, E) result(fe)
   real(8) :: fe 
   real(8), intent(in) :: mu, Te    ! chemical potential and temperature, both in [eV]
   real(8), intent(in) :: E ! energy point [eV]
   real(8) :: eps, arg
   eps = 1.0d-4
   if (Te < eps) then   ! zero temperature, no influence
      if (E <= mu) then
         fe = 1.0d0
      else
         fe = 0.0d0
      endif
   else
      arg = (E - mu)/Te
      if ( arg >= log(HUGE(mu)) ) then   ! too small temperature
         if (E <= mu) then
            fe = 1.0d0
         else
            fe = 0.0d0
         endif
      else ! non-zero temperature
         fe = 1.0d0/(1.0d0 + exp(arg))
      endif
   endif
end function Fermi_function



pure subroutine Find_in_1D_array(Array, Value, Number)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i
   i = 1
   do while (Array(i) .LT. Value-1.0d-10)
      i = i + 1
   enddo
   Number = i
end subroutine Find_in_1D_array

pure subroutine Find_in_2D_array(Array, Value, Indx, Number)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value   ! to be found in the array as near as possible
   integer, INTENT(in) :: Indx    ! index of the array, showing in which colonm we search
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i
   i = 1
   do while (Array(Indx,i) .LT. Value)
      i = i + 1
   enddo
   Number = i
end subroutine Find_in_2D_array

subroutine Find_in_monotonous_1D_array(Array, Value0, Number)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2

   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_1D_array'
        print*, 'cannot proceed, the value of Value0 is', Value0
        pause 'STOPPED WORKING...'
   else
       N = size(Array)
       if (Value0 .LE. Array(1)) then ! it's the first value, no need to search
           i_cur = 1
       else if (Value0 .GE. Array(N)) then ! it's the last value, no need to search
           i_cur = N
       else
           i_1 = 1
           i_2 = N
           i_cur = FLOOR((i_1+i_2)/2.0)
           temp_val = Array(i_cur)
           coun = 0
           do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if ((Value0 .GE. Array(i_cur)) .AND. (Value0 .LT. Array(i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if (temp_val .LE. Value0) then
                   i_1 = i_cur
                else
                   i_2 = i_cur
                endif
                i_cur = FLOOR((i_1+i_2)/2.0)
                temp_val = Array(i_cur)
                coun = coun + 1
                if (coun .GT. 1e3) then
                    print*, 'PROBLEM WITH CONVERGANCE IN'
                    print*, 'Find_in_monotonous_1D_array', coun
                    write(*, '(es25.16,es25.16,es25.16,es25.16)') Value0, Array(i_cur), Array(i_1), Array(i_2)
                    pause 'STOPPED WORKING...'
                endif
           enddo
       endif
   endif    ! isnan
   Number = i_cur
end subroutine Find_in_monotonous_1D_array

subroutine Find_in_monotonous_2D_array(Array, Value0, Indx, Number)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(in) :: Indx    ! index of the array, showing in which colonm we search
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2

   N = size(Array,2)
   i_1 = 1
   val_1 = Array(Indx,i_1)
   i_2 = N
   val_2 = Array(Indx,i_2)
   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(Indx,i_cur)

   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_2D_array'
        print*, 'cannot proceed, the value of Value0 is', Value0
        write(*, '(f25.16,f25.16,f25.16,f25.16)') Value0, Array(Indx,i_cur), Array(Indx,i_1), Array(Indx,i_2)
        pause 'STOPPED WORKING...'
   else
       if (Value0 .LT. Array(Indx,1)) then ! it's the first value, no need to search
           i_cur = 0
       else if (Value0 .GE. Array(Indx,N)) then ! it's the last value, no need to search
           i_cur = N-1
       else
           coun = 0
           do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if ((Value0 .GE. Array(Indx,i_cur)) .AND. (Value0 .LE. Array(Indx,i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if (temp_val .LE. Value0) then
                   i_1 = i_cur
                   val_1 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(Indx,i_cur)
                else
                   i_2 = i_cur
                   val_2 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(Indx,i_cur)
                endif
                coun = coun + 1
                if (coun .GT. 1e3) then
                    print*, 'PROBLEM WITH CONVERGANCE IN'
                    print*, 'Find_in_monotonous_2D_array', coun
                    write(*, '(f25.16,f25.16,f25.16,f25.16)') Value0, Array(Indx,i_cur), Array(Indx,i_1), Array(Indx,i_2)
                    pause 'STOPPED WORKING...'
                endif
           enddo
       endif
   endif    ! isnan
   Number = i_cur+1
end subroutine Find_in_monotonous_2D_array


pure subroutine sort_array_bubble_real(array_to_sort)	! sort accending
   real(8), dimension(:), intent(inout) :: array_to_sort
   real(8) :: temp
   integer :: N,i,j
   logical :: swapped
   N = size(array_to_sort)
   swapped = .false.
   do j = N-1, 1, -1
      do i = 1, j
         if (array_to_sort(i) > array_to_sort(i+1)) then	! swap them
            temp = array_to_sort(i)
            array_to_sort(i) = array_to_sort(i+1)
            array_to_sort(i+1) = temp
            swapped = .true.
         end if
      enddo
      if (.not. swapped) exit	! already sorted, no need to continue
   enddo
end subroutine sort_array_bubble_real

pure subroutine sort_array_bubble_int(array_to_sort)	! sort accending
   integer, dimension(:), intent(inout) :: array_to_sort
   integer :: temp
   integer :: N,i,j
   logical :: swapped
   N = size(array_to_sort)
   swapped = .false.
   do j = N-1, 1, -1
      do i = 1, j
         if (array_to_sort(i) > array_to_sort(i+1)) then	! swap them
            temp = array_to_sort(i)
            array_to_sort(i) = array_to_sort(i+1)
            array_to_sort(i+1) = temp
            swapped = .true.
         end if
      enddo
      if (.not. swapped) exit	! already sorted, no need to continue
   enddo
end subroutine sort_array_bubble_int


pure subroutine extend_array_size_by_one_int(array1)
   integer, dimension(:), allocatable, intent(inout) :: array1
   integer, dimension(:), allocatable :: tmp
   integer :: N_old, N_new
   if (.not.allocated(array1)) then
      allocate(array1(1))
   else
      N_old = size(array1)	! size of the present array
      N_new = N_old + 1	! new size increased by factor of 2 
      allocate( tmp( N_new ) )	! array tp be used temporary storing data
      tmp(1:N_old) = array1	! store the data from old array
      ! internal function in FORTRAN-2008 format:
      call move_alloc( tmp, array1 )	! shift data from old array into the new one
   endif
end subroutine extend_array_size_by_one_int

pure subroutine extend_array_size_by_one_real(array1)
   real(8), dimension(:), allocatable, intent(inout) :: array1
   real(8), dimension(:), allocatable :: tmp
   integer :: N_old, N_new
   if (.not.allocated(array1)) then
      allocate(array1(1))
   else
      N_old = size(array1)	! size of the present array
      N_new = N_old + 1	! new size increased by factor of 2 
      allocate( tmp( N_new ) )	! array tp be used temporary storing data
      tmp(1:N_old) = array1	! store the data from old array
      ! internal function in FORTRAN-2008 format:
      call move_alloc( tmp, array1 )	! shift data from old array into the new one
   endif
end subroutine extend_array_size_by_one_real


pure subroutine extend_array_size_int(array1)
   integer, dimension(:), allocatable, intent(inout) :: array1
   integer, dimension(:), allocatable :: tmp
   integer :: N_old, N_new
   N_old = size(array1)	! size of the present array
   N_new = 2*N_old		! new size increased by factor of 2 
   allocate( tmp( N_new ) )	! array tp be used temporary storing data
   tmp(1:N_old) = array1	! store the data from old array
!    deallocate( array1 )	! free memory from old array
   ! internal function in FORTRAN-2008 format:
   call move_alloc( tmp, array1 )	! shift data from old array into the new one
end subroutine extend_array_size_int

pure subroutine extend_array_size_real(array1)
   real(8), dimension(:), allocatable, intent(inout) :: array1
   real(8), dimension(:), allocatable :: tmp
   integer :: N_old, N_new
   N_old = size(array1)	! size of the present array
   N_new = 2*N_old		! new size increased by factor of 2 
   allocate( tmp( N_new ) )	! array tp be used temporary storing data
   tmp(1:N_old) = array1	! store the data from old array
!    deallocate( array1 )	! free memory from old array
   ! internal function in FORTRAN-2008 format:
   call move_alloc( tmp, array1 )	! shift data from old array into the new one
end subroutine extend_array_size_real


pure subroutine extend_array_size_int_old(array1)
   integer, dimension(:), allocatable, intent(inout) :: array1
   integer N
   integer, dimension(:), allocatable :: array2
   N = size(array1)
   allocate(array2(N))
   array2 = array1
   deallocate(array1)
   allocate(array1(2*N))
   array1(1:N) = array2(1:N)
   deallocate(array2)
end subroutine extend_array_size_int_old

pure subroutine extend_array_size_real_old(array1)
   real(8), dimension(:), allocatable, intent(inout) :: array1
   integer N
   real(8), dimension(:), allocatable :: array2
   N = size(array1)
   allocate(array2(N))
   array2 = array1
   deallocate(array1)
   allocate(array1(2*N))
   array1(1:N) = array2(1:N)
   deallocate(array2)
end subroutine extend_array_size_real_old


subroutine print_time_step(text, num, msec, print_to)
   CHARACTER(len=*) :: text	! text to print out
   real(8), intent(in), optional :: num	! to print out this number
   logical, intent(in), optional :: msec ! print msec or not?
   integer, intent(in), optional :: print_to	! which unit to print to
   character(len=100) :: var 
   integer c1(8) ! time stamps
   integer :: print_to_unit
   if (present(print_to)) then
      print_to_unit = print_to	! print into theis unit (file or screen)
   else
      print_to_unit = 6	! print on the screen by default
   endif
   
   call date_and_time(values=c1) ! standard FORTRAN time and date
   if (present(num) .and. present(msec)) then
      write(var,'(f16.2)') num
      write(print_to_unit, 1002) trim(adjustl(text))//' '//trim(adjustl(var))//' fs at ', c1(5), c1(6), c1(7), c1(8), c1(3), c1(2), c1(1)
   elseif (present(msec)) then
      write(print_to_unit, 1002) trim(adjustl(text))//' fs at ', c1(5), c1(6), c1(7), c1(8), c1(3), c1(2), c1(1)
   elseif (present(num)) then
      write(var,'(f16.2)') num
      write(print_to_unit, 1001) trim(adjustl(text))//' '//trim(adjustl(var))//' fs at ', c1(5), c1(6), c1(7), c1(3), c1(2), c1(1)
   else
      write(print_to_unit, 1001) trim(adjustl(text))//' fs at ', c1(5), c1(6), c1(7), c1(3), c1(2), c1(1)
   endif
   ! Formats for printing:
   1001 format (a, i2.2, ':', i2.2, ':', i2.2, ' on ', i2.2, '/', i2.2, '/', i4.4)
   1002 format (a, i2.2, ':', i2.2, ':', i2.2, ':', i3.3, ' on ', i2.2, '/', i2.2, '/', i4.4)
end subroutine print_time_step


pure subroutine trim_zeros(in_line)
   character(len=*), intent(inout) :: in_line
   integer :: i, j, leng
   i = 0
   leng = LEN(trim(in_line))
   j = leng
   do while (in_line(j:j) == '0') ! trim extra zeros from the end of the line
      in_line(j:j) = ' ' ! delete last characters
      i = i + 1
      j = leng-i
   enddo
   if (in_line(j:j) == '.') in_line(j+1:j+1) = '00'
end subroutine trim_zeros


subroutine print_time(text, ind, iter, ctim, print_to)
   CHARACTER(len=*) :: text	! text to print out
   integer, intent(in), optional :: ind, iter	! ind = to print out miliseconds or not; iter = integer number to print out
   integer, intent(in), optional ::  ctim(8)	! given time to print it out
   integer, intent(in), optional :: print_to	! which unit to print to
   integer c1(8) ! time stamps
   integer :: print_to_unit
   character(len=100) :: var 
   if (present(print_to)) then
      print_to_unit = print_to	! print into theis unit (file or screen)
   else
      print_to_unit = 6	! print on the screen by default
   endif
   
   call date_and_time(values=c1) ! standard FORTRAN time and date
   if (present(ind)) then
      if (present(ctim)) then
         if (present(iter)) then
            write(var,'(i12)') iter
            write(print_to_unit, 1002) ' '//trim(adjustl(text))//' '//trim(adjustl(var))//' ', ctim(5), ctim(6), ctim(7), ctim(8), ctim(3), ctim(2), ctim(1)
         else
            write(print_to_unit, 1002) ' '//trim(adjustl(text))//' ', ctim(5), ctim(6), ctim(7), ctim(8), ctim(3), ctim(2), ctim(1)
         endif
      else
         if (present(iter)) then
            write(var,'(i12)') iter
            write(print_to_unit, 1002) ' '//trim(adjustl(text))//' '//trim(adjustl(var))//' ', c1(5), c1(6), c1(7), c1(8), c1(3), c1(2), c1(1)
         else
            write(print_to_unit, 1002) ' '//trim(adjustl(text))//' ', c1(5), c1(6), c1(7), c1(8), c1(3), c1(2), c1(1)
         endif
      endif
   else
      if (present(ctim)) then
         if (present(iter)) then
            write(var,'(i12)') iter
            write(print_to_unit, 1001) ' '//trim(adjustl(text))//' '//trim(adjustl(var))//' ', ctim(5), ctim(6), ctim(7), ctim(8), ctim(3), ctim(2), ctim(1)
         else
            write(print_to_unit, 1001) ' '//trim(adjustl(text))//' ', ctim(5), ctim(6), ctim(7), ctim(3), ctim(2), ctim(1)
         endif
      else
         if (present(iter)) then
            write(var,'(i12)') iter
            write(print_to_unit, 1001) ' '//trim(adjustl(text))//' '//trim(adjustl(var))//' ', c1(5), c1(6), c1(7), c1(3), c1(2), c1(1)
         else
            write(print_to_unit, 1001) ' '//trim(adjustl(text))//' ', c1(5), c1(6), c1(7), c1(3), c1(2), c1(1)
         endif
      endif
   endif
   ! Formats for printing:
   1001 format (a, i2.2, ':', i2.2, ':', i2.2, ' on ', i2.2, '/', i2.2, '/', i4.4)
   1002 format (a, i2.2, ':', i2.2, ':', i2.2, ':', i3.3, ' on ', i2.2, '/', i2.2, '/', i4.4)
end subroutine print_time


subroutine print_energy(E, print_txt, print_to)
   real(8), intent(in) :: E ! energy in eV
   character(*), intent(in) :: print_txt
   integer, intent(in), optional :: print_to    ! which unit to print to
   integer :: print_to_unit
   real(8) :: E_converted
   character(6) :: E_units
   character(30) :: E_txt
   
   if (present(print_to)) then
      print_to_unit = print_to	! print into theis unit (file or screen)
   else
      print_to_unit = 6	! print on the screen by default
   endif
   
   if (E <= 1.0d-3) then ! meV
      E_converted = E*1.0d3
      E_units = '[meV]'
   elseif (E >= 1.0d12 ) then ! TeV
      E_converted = E*1.0d-12
      E_units = '[TeV]'
   elseif (E >= 1.0d9 ) then ! GeV
      E_converted = E*1.0d-9
      E_units = '[GeV]'
   elseif (E >= 1.0d6 ) then ! MeV
      E_converted = E*1.0d-6
      E_units = '[MeV]'
   elseif (E >= 1.0d3 ) then ! keV
      E_converted = E*1.0d-3
      E_units = '[keV]'
   else ! eV
      E_converted = E
      E_units = '[eV]'
   endif
   
   write(E_txt, '(f10.3)') E_converted
   write(print_to_unit, '(a,a,a)') ' '//trim(adjustl(print_txt))//' ', trim(adjustl(E_txt)), ' '//trim(adjustl(E_units))
   
end subroutine print_energy



pure subroutine parse_yes_no(string, yesno)
   character(*), intent(in) :: string ! figure out whether it is 'yes' or 'no'
   logical, intent(out) :: yesno ! yes = true, no = false
   select case (trim(adjustl(string))) 
   case ('y', 'Y', 'yes', 'YES', 'Yeah', 'Yes', 'yEs', 'yeS', 'YEs', 'yES', 'YeS', '1', 'true', 'Da', 'Ja', 'da', 'aha', 'ja', 'Aha')
      yesno = .true.
   case default
      yesno = .false.
   end select
end subroutine parse_yes_no


function fast_pow(x,n)
! this function calculates integer power of a variable x for powers up to 13 in a fast manner:
   real(8) :: fast_pow
   real(8), intent(in) :: x
   integer, intent(in) :: n
   real(8) x2, x3, x4, x5, x6, x12
   select case(n)
   case(0)
      fast_pow = 1.0d0
   case (1)
      fast_pow = x
   case (2)
      fast_pow = x*x
   case (3)
      x2 = x*x
      fast_pow = x2*x
   case (4)
      x2 = x*x
      fast_pow = x2*x2
   case (5)
      x2 = x*x
      x3 = x2*x
      fast_pow = x3*x2
   case (6)
      x2 = x*x
      x3 = x2*x
      fast_pow = x3*x3
   case (7)
      x2 = x*x
      x3 = x2*x
      x4 = x2*x2
      fast_pow = x4*x3
   case (8)
      x2 = x*x
      x3 = x2*x
      x4 = x2*x2
      fast_pow = x4*x4
   case (9)
      x2 = x*x
      x3 = x2*x
      x4 = x2*x2
      x5 = x2*x3
      fast_pow = x5*x4
   case (10)
      x2 = x*x
      x3 = x2*x
      !x4 = x2*x2
      x5 = x2*x3
      fast_pow = x5*x5
   case (11)
      x2 = x*x
      x3 = x2*x
      !x4 = x2*x2
      x5 = x2*x3
      x6 = x3*x3
      fast_pow = x6*x5
   case (12)
      x2 = x*x
      x3 = x2*x
      !x4 = x2*x2
      !x5 = x2*x3
      x6 = x3*x3
      fast_pow = x6*x6
   case (13)
      x2 = x*x
      x3 = x2*x
      !x4 = x2*x2
      !x5 = x2*x3
      x6 = x3*x3
      fast_pow = x6*x6*x
   case (24)
      x2 = x*x
      x3 = x2*x
      x6 = x3*x3
      x12 = x6*x6
      fast_pow = x12*x12
   case default
      print*, 'This power is not supported in function "fast_pow"'
      pause 'What now?'
   end select
end function fast_pow


end module Little_subroutines
