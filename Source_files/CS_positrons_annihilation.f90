! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains all cross sections for positron annihilatin
! [1] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2014: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2015)


module CS_positrons_annihilation
use Universal_constants
use CDF_Ritchi
use Relativity
use Objects
use Read_input_data, only: m_positron_CS, m_input_folder, m_folder_materials, m_positron_annihilation
use Dealing_with_files, only: read_file
use CS_general_tools, only: MFP_from_sigma

implicit none

 contains

!-----------------------------------------------------------------------------
! Subroutines for position annihilation:
subroutine get_positron_annihilation(Material, numpar, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material	!material parameters of each target that it's constructed of
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------------------------
   integer :: i, j, k, m, Nsiz, N_targets, N_elements, N_shells, FN, Reason, count_lines
   real(8) :: sigma, MFP, N_elem, elem_contrib, Ne_tot
   character(200) :: Path, Folder_with_CS, command, File_name, Model_name, Path_total
   logical :: file_exist, read_well
   !--------------------------------
   character, pointer :: path_sep

   write(*, '(a)', advance='no') ' Obtaining positron annihilation cross sections...'

   path_sep => numpar%path_sep
   Folder_with_CS = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_positron_CS))   ! Electron CSs are storred in this folder   

   N_targets = size(Material)	! that's how many different targets user specified
   ! Get positron MFPs:
   TRGT:do i = 1, N_targets	! for each target
      ! Set part of the file name corresponding to the model used for positron annihilation CS:
      select case (numpar%Pos_annih)
      case (1)	! Heitler model
         Model_name = 'Heitler'
      case default	! exclude
         Model_name = 'NO'   
      end select

      inquire(DIRECTORY=trim(adjustl(Folder_with_CS)),exist=file_exist)    ! check if input file excists
      if (.not.file_exist) then	! to make sure that such a folder is present (even if empty)
         ! Create a new directory for output files:
         command='mkdir '//trim(adjustl(Folder_with_CS))	! to create a folder use this command
         CALL system(command)	! create the folder
      endif

      ! For positron CSs use the same grid as for photons:
      Nsiz = size(Material(i)%El_elastic_total%E)   ! grid for total cross section
      allocate(Material(i)%Pos_annihil_total%E(Nsiz))
      Material(i)%Pos_annihil_total%E = Material(i)%El_elastic_total%E
      allocate(Material(i)%Pos_annihil_total%Total(Nsiz))
      allocate(Material(i)%Pos_annihil_total%Total_MFP(Nsiz))
      Material(i)%Pos_annihil_total%Total(:) = 0.0d0
      Material(i)%Pos_annihil_total%Total_MFP(:) = 0.0d0

      ! Get annihilation cross sections per electron:
      File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_positron_annihilation))//'_per_electron_'//trim(adjustl(Model_name))//'.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!       if (file_exist) then	! only create it if file does not exist
      if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
         count_lines = 0
         do m = 1, Nsiz	! for all energy grid points:
            read(FN,'(es,es)', IOSTAT=Reason) Material(i)%Pos_annihil_total%E(m), Material(i)%Pos_annihil_total%Total(m)
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not.read_well) then
               close(FN)	! redo the file
               goto 8391
            endif
         enddo ! m = 1, Nsiz
      else
8391     open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            call get_pos_annihilation_CS(Material(i)%Pos_annihil_total%E(m), numpar, Material(i)%Pos_annihil_total%Total(m))   ! below
            write(FN,'(es,es)') Material(i)%Pos_annihil_total%E(m), Material(i)%Pos_annihil_total%Total(m)
         enddo ! m = 1, Nsiz
      endif
      close(FN)

      ! Get how many electrons we have per atom in each element:
      N_elem = dble(SUM(Material(i)%Elements(:)%percentage))	! number of elements in this compound material
      Ne_tot = 0.0d0
      N_elements = size(Material(i)%Elements)	! that's how many different elements are in this target
      do j =1, N_elements	! for each element in this target
         elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         Ne_tot = Ne_tot + dble(Material(i)%Elements(j)%Zat)*elem_contrib
      enddo ! j
      ! Convert CSs to those per atom:
      Material(i)%Pos_annihil_total%Total(:) = Material(i)%Pos_annihil_total%Total(:) * Ne_tot  ! per atom

      ! Also get MFPs:
      do m = 1, Nsiz	! for all energy grid points:
         Material(i)%Pos_annihil_total%Total_MFP(m) = MFP_from_sigma(Material(i)%Pos_annihil_total%Total(m),  Material(i)%At_Dens)    ! module "CS_general_tools"
      enddo

   enddo TRGT

   write(*, '(a)') ' Done.'
!    print*, 'Positron annihilation cross sections are obtained.'

   nullify(path_sep)
end subroutine get_positron_annihilation

 

! Interface to select the model of positron elastic scattering cross section:
pure subroutine get_pos_annihilation_CS(Ee, numpar, sigma)
   real(8), intent(in) :: Ee	! [eV] positron kinetic energy
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), intent(out) :: sigma	! [A^2] cross section (per electron, not per atom)
   select case (numpar%Pos_annih)	! annihilation scattering: 0=excluded, 1=Heitler
   case (1)	! Heitler model
      sigma = get_Heitler_CS(Ee)  ! below
   case default	! exclude
      sigma = 0.0d0
   end select
end subroutine get_pos_annihilation_CS


pure function get_Heitler_CS(Ee) result(CS)
   real(8) CS   ! [A^2] cross section (per electron, not per atom)
   real(8), intent(in) :: Ee    ! [eV] kinetic energy of positron
   real(8) :: v, g, g1
   v = velosity_from_kinetic_energy(Ee, g_me, afs=.false.)   ! [m/s] module "Relativity"
   g = gamma_factor(v)  ! module "Relativity"
   g1 = (g*g-1.0d0)
   ! Eq. (3.189) p.182  [1]:
   CS = g_Pi*g_r0*g_r0/((g+1.0d0)*g1) * ( (g*g + 4.0d0*g + 1.0d0)*log(g + sqrt(g1)) - (3.0d0 + g)*sqrt(g1) )
end function get_Heitler_CS


pure function get_Heitler_intCS(Ee, ksi) result(CS)
   real(8) CS   ! [A^2] cross section (per electron, not per atom)
   real(8), intent(in) :: ksi   ! emitted photon energy normalized to positron total energy, we integrate until
   real(8), intent(in) :: Ee    ! [eV] kinetic enegry of positron
   real(8) :: v, g, g1, ksi1
   v = velosity_from_kinetic_energy(Ee, g_me, afs=.false.)   ! [m/s] module "Relativity"
   g = gamma_factor(v)  ! module "Relativity"
   g1 = g + 1.0d0
   ksi1 = 1.0d0 - ksi
   ! Analytical integral of Eq.(3.187) p.182 [1]:
   CS = -2.0d0 * g1*g1 * ksi + 1.0d0/ksi - 1.0d0/ksi1 + log(ksi/ksi1) * (g*g+4.0d0*g+1.0d0)
   CS = CS * g_Pi*g_r0*g_r0/((g-1.0d0)*g1*g1)
end function get_Heitler_intCS



pure function get_Heitler_dCS(Ee, ksi) result(CS)
   real(8) CS   ! [A^2] cross section (per electron, not per atom)
   real(8), intent(in) :: ksi   ! emitted photon energy normalized to positron total energy
   real(8), intent(in) :: Ee    ! [eV] kinetic enegry of positron
   real(8) :: v, g, g1
   v = velosity_from_kinetic_energy(Ee, g_me, afs=.false.)   ! [m/s] module "Relativity"
   g = gamma_factor(v)  ! module "Relativity"
   g1 = (g*g-1.0d0)
   ! Eq.(3.187) p.182  [1]:
   CS = g_Pi*g_r0*g_r0/((g+1)*g1) * ( get_zeta(g,ksi) + get_zeta(g,1.0d0 - ksi) )    ! below
end function get_Heitler_dCS


pure function get_zeta(g, ksi) result(zeta)
   real(8) zeta
   real(8), intent(in) :: ksi   ! Eq. (3.183)
   real(8), intent(in) :: g ! gamma
   real(8) :: g1
   g1 = g + 1.0d0
   zeta = -g1*g1 + (g*g + 4.0d0*g + 1.0d0)/ksi - 1.0d0/(ksi * ksi)  ! Eq.(3.188) [1]
end function get_zeta


pure function ksi_min(Ekin)
   real(8) ksi_min  ! Eq.(3.186) [1]
   real(8), intent(in) :: Ekin
   real(8) :: gam
   gam = 1.0d0 + Ekin/g_me
   ksi_min = 1.0d0/( 1.0d0 +  gam + sqrt(gam*gam - 1.0d0) )
end function ksi_min


pure function theta_annih(Epos, Eph) result(theta)
   real(8) theta
   real(8), intent(in) :: Epos, Eph
   real(8) :: g, ksi, cos_th
   g = gamma_factor_from_Ekin(Epos, g_me)   ! module "Relativity"
   ksi = ksi_annih(Eph, Epos)   ! below
   cos_th = sqrt(g*g - 1.0d0) * (g + 1.0d0 - 1.0d0/ksi) ! Eq.(3.184) [1]
   theta = ACOS(cos_th)
end function theta_annih


pure function theta_plus_annih(Epos, Eph) result(theta)
   real(8) theta
   real(8), intent(in) :: Epos, Eph
   real(8) :: g, ksi, cos_th
   g = gamma_factor_from_Ekin(Epos, g_me)   ! module "Relativity"
   ksi = ksi_annih(Eph, Epos)   ! below
   cos_th = sqrt(g*g - 1.0d0) * (g + 1.0d0 - 1.0d0/(1.0d0-ksi)) ! Eq.(3.185) [1]
   theta = ACOS(cos_th)
end function theta_plus_annih


pure function ksi_annih(Eph, Epos) result(ksi)
   real(8) ksi
   real(8), intent(in) :: Eph, Epos ! [eV] photon energy, and positron kinetic energy
   ksi = Eph/(Epos + 2.0d0*g_me_eV)    ! Eq.(3.183) [1]
end function ksi_annih

pure function Eph_from_ksi_annih(Epos, ksi) result(Eph)
   real(8) Eph  ! photon energy [eV]
   real(8), intent(in) :: ksi, Epos ! normalized photon energy, and positron kinetic energy [eV] 
   Eph = ksi/(Epos + 2.0d0*g_me_eV)    ! inverse of Eq.(3.183) [1]
end function Eph_from_ksi_annih



end module CS_positrons_annihilation
