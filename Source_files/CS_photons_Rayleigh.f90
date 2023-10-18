! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains cross sections for phoon scattering: Rayleigh / Thomson
! References used in the module:
! [1]  F. Salvat, J. M. Fernandez-Varea, E. Acosta, J. Sempau "PENELOPE – A Code System for Monte Carlo Simulation of Electron and Photon Transport", OECD (2001)


module CS_photons_Rayleigh
use Universal_constants
use Objects
use Read_input_data, only: m_input_folder, m_photon_CS, m_photon_Rayleigh, m_folder_materials
use Dealing_with_files, only: read_file
use CS_general_tools, only: MFP_from_sigma
use Little_subroutines, only: interpolate_data_single

implicit none

real(8) :: m_Thom_factor, m_five_sixteenth

parameter (m_Thom_factor = 20.6074d0)	! [1] Constant for Eq.(2.5)
parameter (m_five_sixteenth = 5.0d0/16.0d0)

 contains

 
!---------------------------------------------------------------------------------
! Rayleigh (photon elastic scattering) subroutines:

subroutine get_photon_Rayleigh(Material, numpar, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material      !material parameters of each target that it's constructed of
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------------------------
   real(8) :: sigma, MFP, N_elem, elem_contrib
   integer :: i, j, k, N_targets, N_elements, Ngrid, N_shells, FN, m
   integer :: Reason, count_lines, Nsiz
   character(200) :: Folder_with_CS, Path, command, Model_name, File_name, Path_total
   logical :: file_exist, read_well
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   
   write(*, '(a)', advance='no') ' Obtaining photon Rayleigh scattering cross sections...'
   
   path_sep => numpar%path_sep
   Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_photon_CS))	! Electron CSs are storred in this folder
   
   N_targets = size(Material)	! that's how many different targets user specified
   
   ! Get electron MFPs:
   TRGT:do i = 1, N_targets	! for each target
      N_elem = dble(SUM(Material(i)%Elements(:)%percentage))	! number of elements in this compound material
      N_elements = size(Material(i)%Elements)	! that's how many different elements are in this target
      LMNT:do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         ! Check if folder for this element already exists:
         Folder_with_CS = trim(adjustl(Path))//path_sep//trim(adjustl(Element%Name))
         inquire(DIRECTORY=trim(adjustl(Folder_with_CS)),exist=file_exist)    ! check if input file excists
         if (.not.file_exist) then	! to make sure that such a folder is present (even if empty)
            ! Create a new directory for output files:
            command='mkdir '//trim(adjustl(Folder_with_CS))	! to create a folder use this command
            CALL system(command)	! create the folder
         endif
         !-------------------------------------
          ! Set part of the file name corresponding to the model used for pair creation CS:
         select case (numpar%Ph_Thomson)
         case (1)	! PENELOPE
            Model_name = 'PENELOPE'
         case default	! exclude
            Model_name = 'NO'
         end select
         
         ! For pair creation CSs use the same grid as for photon absorption:
         Ngrid = size(Element%Phot_absorption%E)
         allocate(Element%Phot_coherent%E(Ngrid))
         Element%Phot_coherent%E = Element%Phot_absorption%E
         allocate(Element%Phot_coherent%Total(Ngrid))
         allocate(Element%Phot_coherent%Total_MFP(Ngrid))
         
         File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_photon_Rayleigh))//'_total_'//trim(adjustl(Model_name))//'.dat'
!          print*, trim(adjustl(File_name))
         
         inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!          if (file_exist) then	! read from the file
         if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
            open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
            do m = 1, Ngrid	! for all energy grid points:
               E => Element%Phot_coherent%E(m)
               read(FN,'(es,es,es)', IOSTAT=Reason) E, Element%Phot_coherent%Total(m), Element%Phot_coherent%Total_MFP(m) 
               call read_file(Reason, count_lines, read_well)	! module "Dealing_wil_files"
               if (.not.read_well) then
                  close(FN)	! redo the file
                  goto 8391
               endif
            enddo ! m = 1, Ngrid
         else	! only create it if file does not exist 
8391     open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
            do m = 1, Ngrid	! for all energy grid points:
               E => Element%Phot_coherent%E(m)
               call get_Rayleigh_CS(E, Element, numpar, sigma)	! below
               Element%Phot_coherent%Total(m) = sigma   ! [A^2]
               Element%Phot_coherent%Total_MFP(m) = MFP_from_sigma(Element%Phot_coherent%Total(m),  1.0d24)    ! module "CS_general_tools"
               write(FN,'(es,es,es)') E, sigma, Element%Phot_coherent%Total_MFP(m) 
            enddo ! m = 1, Ngrid
         endif
         close(FN)
         ! Normalize MFPs and Se for real material density:
         ! Take into account the density of atoms of this particular kind:
         elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         Element%Phot_coherent%Total_MFP(:) = Element%Phot_coherent%Total_MFP(:) * 1.0d24/(Material(i)%At_Dens * elem_contrib)
      enddo LMNT
      
      ! Total CS:
      ! Where to save:
      Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
      Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
      File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_photon_Rayleigh))//'_total_'//trim(adjustl(Model_name))//'.dat'
      ! Get the cross sections to save:
      Nsiz = size(Material(i)%Ph_absorption_total%E)   ! grid for total cross section is different from core-shells grid
      allocate(Material(i)%Ph_coherent_total%E(Nsiz))
      Material(i)%Ph_coherent_total%E = Material(i)%Ph_absorption_total%E
      allocate(Material(i)%Ph_coherent_total%Total(Nsiz))
      allocate(Material(i)%Ph_coherent_total%Total_MFP(Nsiz))
      Material(i)%Ph_coherent_total%Total(:) = 0.0d0
      Material(i)%Ph_coherent_total%Total_MFP(:) = 0.0d0
       ! Sum up CSs from each element:
      do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         do m = 1, Nsiz	! for all energy grid points:
            E => Material(i)%Ph_coherent_total%E(m)    ! photon energy [eV]
            call interpolate_data_single(Element%Phot_coherent%E,  Element%Phot_coherent%Total(:), E, sigma) ! module "Little_subroutines"
            call interpolate_data_single(Element%Phot_coherent%E,  Element%Phot_coherent%Total_MFP(:), E, MFP) ! module "Little_subroutines"
            ! Add them into the arrays:
            Material(i)%Ph_coherent_total%Total(m) = Material(i)%Ph_coherent_total%Total(m) + sigma*elem_contrib
            Material(i)%Ph_coherent_total%Total_MFP(m) = Material(i)%Ph_coherent_total%Total_MFP(m) + 1.0d0/MFP !*elem_contrib ! to be inverted
         enddo ! m
      enddo ! j
      ! Invert to get MFP:
      Material(i)%Ph_coherent_total%Total_MFP(:) = 1.0d0/Material(i)%Ph_coherent_total%Total_MFP(:) ! [A]
      ! Save the total cross section into the file:
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
      if (.not.file_exist .or. numpar%recalculate_MFPs) then	! only create it if file does not exist
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            write(FN,'(es,es,es)') Material(i)%Ph_coherent_total%E(m), Material(i)%Ph_coherent_total%Total(m), Material(i)%Ph_coherent_total%Total_MFP(m)
         enddo ! m = 1, Nsiz
         close(FN)
      endif
      
   enddo TRGT
   
   write(*, '(a)') ' Done.'
!    print*, 'Photon Rayleigh scattering cross sections are obtained.'
   nullify(path_sep, Element, E)
end subroutine get_photon_Rayleigh

 
! Rayleigh scattering cross section per element:
subroutine get_Rayleigh_CS(Eph, Element, numpar, sigma, mu_max_in)
   real(8), intent(in) :: Eph	! [eV] photon energy
   type(Atom_kind), intent(in) :: Element	! data for this element
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), intent(out) :: sigma	! [A^2] cross section
   real(8), intent(in), optional :: mu_max_in   ! maximal angle to integrate up to
   select case (numpar%Ph_Thomson)
   case (1)	! PENELOPE
      if (present(mu_max_in)) then
         sigma = Rayleigh_CS(Eph, Element, mu_max_in)   ! see below
      else
         sigma = Rayleigh_CS(Eph, Element)	! see below
      endif
   case default	! exclude
      sigma = 0.0d0
   end select
end subroutine get_Rayleigh_CS


function Rayleigh_CS(Eph, Element, mu_max_in) result(sigma)
   real(8) :: sigma	! [A^2] cross section
   real(8), intent(in) :: Eph	! [eV] photon energy
   type(Atom_kind), intent(in), target :: Element	! data for this element
   real(8), intent(in), optional :: mu_max_in   ! maximal angle to integrate up to
   !------------------------------------
   real(8) :: dSigma, dSigma0, dSigma_mid, dS, eps, sig_step, dblZ	! d Sigma / s Omega
   real(8) :: mu, mu0, dmu, dmu_half, mu_mid, mu_min, mu_max, dmu_max, dmu_min
   integer :: i, Ngrid
   
   NEGLEC:if (Eph <= 1.0d6) then	! and only then calculated the cross sections; above this energy, the cross sections calculation is unstable! But there Rayleigh can always be neglected
      dblZ = dble(Element%Zat)
      eps = 1.0d-17	! precision limit
      dS = 0.01d0	! maximal allowed change in dSigma per step
      Ngrid = 100	! grid point for integration over mu
   
      mu_min = -1.0d0	! cos(-Pi)
      if (present(mu_max_in)) then
         if (mu_max_in >= 1.0d0) then
            mu_max = 1.0d0  ! cos(0)
         else
            mu_max = mu_max_in
         endif
      else
         mu_max = 1.0d0 ! cos(0)
      endif
      mu = mu_min		! starting from cos(-Pi)
      dmu_max = (1.0d0 - mu_min)/dble(Ngrid)	! maximal allowed integration step
      dmu_min = dmu_max/100.0d0
      dmu = dmu_max	! start with it, and later reduce if needed
      dmu_half = dmu*0.50d0
      ! Start integration of differential cross section to obtain the total one:
      sigma = 0.0d0
      dSigma0 = Rayleigh_DCS(Eph, mu, dblZ, Element%form_a) 	! below
      i = 0
      do while (mu < mu_max)
         i = i + 1	! steps counter
         ! Simpson 3/8 method of integration:
         mu0 = mu
         mu = mu + dmu
         if (mu > mu_max) then	! if by chance we exceeded the limit
            mu = mu_max
            dmu = mu_max - mu0
            dmu_half = 0.5d0*dmu
         endif
         dSigma = Rayleigh_DCS(Eph, mu, dblZ, Element%form_a) 	! below
         ! Adaptive step: if it's too large, reduce it:
         if ((dSigma > eps) .and. (dSigma0 > eps)) then	! makes sense to take care of integration
            ! Check if it's too large change per step:
!           print*, 'Sigma', dSigma, dSigma0
            do while ((ABS(dSigma - dSigma0) > dS*dSigma0) .and. (dmu > dmu_min))
               dmu = dmu*0.5d0	! reduce timestep
               dmu_half = 0.5d0*dmu	! adjust half-step correspondingly
               mu = mu0 + dmu
               if (mu > mu_max) exit
               dSigma = Rayleigh_DCS(Eph, mu, dblZ, Element%form_a) 	! below
               if (dmu < dmu_min) exit	! too smal step
            enddo
            ! Check if it's too small change per step (inefficient):
            do while ((ABS(dSigma - dSigma0) < 0.5d0*dS*dSigma0) .and. (dmu < dmu_max))
               dmu = dmu + dmu_min	! increase timestep
               dmu_half = 0.5d0*dmu	! adjust half-step correspondingly
               mu = mu0 + dmu
               if (mu > mu_max) exit
               dSigma = Rayleigh_DCS(Eph, mu, dblZ, Element%form_a) 	! below
               if (dmu > dmu_max) exit	! too large step
            enddo
         else	! reset the step
            dmu = dmu_max
            dmu_half = dmu*0.5d0
            mu = mu0 + dmu
         endif
         ! Now proceed with the integration:
         if (mu > mu_max) then	! if by chance we exceeded the limit
            mu = mu_max
            dmu = mu_max - mu0
            dmu_half = dmu*0.50d0
         endif
         mu_mid =  mu0 + dmu_half
         dSigma_mid = Rayleigh_DCS(Eph, mu_mid, dblZ, Element%form_a) 	! below
         ! Add up the contributions according to Simpson-3/8 scheme:
         sig_step = dmu/6.0d0*(dSigma0 + 4.0d0*dSigma_mid + dSigma)
         sigma = sigma + sig_step
!         write(*,'(a,f,f,es,es)') 'S', mu, dmu, dSigma, sigma
         ! Save the data for the next point of integration:
         dSigma0 = dSigma
      enddo
      sigma = 2.0d0*g_Pi*sigma
   else NEGLEC
      sigma = 0.0d0	! ignor negligibly small cross sections
   endif NEGLEC
end function Rayleigh_CS


pure function Rayleigh_DCS(Eph, mu, Z, a) result(dSigma)
   real(8) dSigma
   real(8), intent(in) :: Eph, mu, Z, a(5)
   real(8) :: q,  FF, SigmTh
   SigmTh = g_r0*g_r0*0.5d0*(1.0d0 + mu*mu)	! Thompson cross section, [1] Eq.(2.3)
   q = 2.0d0*Eph*g_e/g_cvel*sqrt(2.0d0*(1.0d0 - mu))	! [1] Eq.(2.4)
   FF = form_factor(q, a, Z)	! below
   dSigma = SigmTh * FF* FF	! [1] Eq.(2.2)
end function Rayleigh_DCS


pure function form_factor(q, a, Z) result(FF)
   real(8) FF
   real(8), intent(in) :: q, Z, a(5)
   real(8) fxz, Fk, x, x2, x4, demf, b, CapQ, al, mc
   mc = g_me*g_cvel
   x = q/mc*m_Thom_factor     ! [1] Eq.(2.5)
   x2 = x*x
   x4 = x2*x2
   demf = 1.0d0 + a(4)*x2 + a(5)*x4
   fxz = Z*(1.0d0 + a(1)*x2 + a(2)*x*x2 + a(3)*x4) / (demf*demf)     ! [1] Eq.(2.7)
   FF = fxz
   if ((Z > 10.0d0) .and. (fxz < 2.0d0)) then
      al = g_alpha*(Z - m_five_sixteenth)    ! [1] Eq.(2.9)
      b = sqrt(1.0d0 - al*al)    ! [1] Eq.(2.9)
      CapQ = q/(2.0d0*mc*al)     ! [1] Eq.(2.9)
      Fk = sin(2.0d0*b*atan(CapQ))/(b*CapQ*(1.0d0 + CapQ*CapQ)**b)   ! [1] Eq.(2.8)
      if (FF < Fk) FF = Fk
   endif
end function form_factor



end module CS_photons_Rayleigh
