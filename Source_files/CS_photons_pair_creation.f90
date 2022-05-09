! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains cross sections for phoon scattering: pair production
! References used in the module:
! [1]  F. Salvat, J. M. Fernandez-Varea, E. Acosta, J. Sempau "PENELOPE – A Code System for Monte Carlo Simulation of Electron and Photon Transport", OECD (2001)


module CS_photons_pair_creation
use Universal_constants
use Objects
use Read_input_data, only: m_input_folder, m_photon_CS, m_photon_pair, m_folder_materials
use Dealing_with_files, only: read_file
use CS_general_tools, only: MFP_from_sigma
use Little_subroutines, only: interpolate_data_single, Find_in_array_monoton

implicit none

real(8) :: m_Cr, m_two_third

parameter (m_Cr = 1.0093d0)	! [1], Page 59, under Eq.(2.92)
parameter (m_two_third = 2.0d0/3.0d0)

 contains

 
 
!---------------------------------------------------------------------------------
! Pair production subroutines:

subroutine get_photon_pair_creation(Material, numpar, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material	!material parameters of each target that it's constructed of
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------------------------
   real(8) :: sigma, N_elem, elem_contrib, MFP
   integer :: i, j, k, N_targets, N_elements, Ngrid, N_shells, FN, m
   integer :: Reason, count_lines, Nsiz, i_closest
   character(200) :: Folder_with_CS, Path, command, Model_name, File_name, Path_total
   logical :: file_exist, read_well
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   
   write(*, '(a)', advance='no') ' Obtaining photon pair creation scattering cross sections...'
   !write(*, '(a)') ' Obtaining photon pair creation scattering cross sections...'
   
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
         select case (numpar%Ph_Pairs)
         case (1)	! PENELOPE
            Model_name = 'PENELOPE'
         case default	! exclude
            Model_name = 'NO'
         end select
         
         ! For pair creation CSs use the same grid as for photon absorption:
         Ngrid = size(Element%Phot_absorption%E)
         allocate(Element%Phot_pair%E(Ngrid))
         Element%Phot_pair%E = Element%Phot_absorption%E
         allocate(Element%Phot_pair%Total(Ngrid))
         allocate(Element%Phot_pair%Total_MFP(Ngrid))
         
         File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_photon_pair))//'_total_'//trim(adjustl(Model_name))//'.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!          if (file_exist) then	! read from the file
         if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
            open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
            do m = 1, Ngrid	! for all energy grid points:
               read(FN,'(es,es,es)', IOSTAT=Reason) Element%Phot_pair%E(m), Element%Phot_pair%Total(m), Element%Phot_pair%Total_MFP(m)
               call read_file(Reason, count_lines, read_well)	! module "Dealing_wil_files"
               if (.not.read_well) then
                  close(FN)	! redo the file
                  goto 8391
               endif
!                write(*,'(i5,es,es,es)')  m, Element%Phot_pair%E(m), Element%Phot_pair%Total(m), Element%Phot_pair%Total_MFP(m)
            enddo ! m = 1, Ngrid
         else	! only create it if file does not exist 
8391     open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
            do m = 1, Ngrid	! for all energy grid points:
               E => Element%Phot_pair%E(m)
               call get_pair_CS(E, Element, numpar, sigma)	! below
               Element%Phot_pair%Total(m) = sigma
               Element%Phot_pair%Total_MFP(m) = MFP_from_sigma(Element%Phot_pair%Total(m),  1.0d24)    ! module "CS_general_tools"
               write(FN,'(es,es,es)') E, sigma, Element%Phot_pair%Total_MFP(m)
            enddo ! m = 1, Ngrid
         endif
         close(FN)
         ! Normalize MFPs and Se for real material density:
         ! Take into account the density of atoms of this particular kind:
         elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         Element%Phot_pair%Total_MFP(:) = Element%Phot_pair%Total_MFP(:) * 1.0d24/(Material(i)%At_Dens * elem_contrib)
      enddo LMNT
      
      ! Total CS:
      ! Where to save:
      Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
      Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
      File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_photon_pair))//'_total_'//trim(adjustl(Model_name))//'.dat'
      ! Get the cross sections to save:
      Nsiz = size(Material(i)%Ph_absorption_total%E)   ! grid for total cross section is different from core-shells grid
      allocate(Material(i)%Ph_pair_total%E(Nsiz))
      Material(i)%Ph_pair_total%E = Material(i)%Ph_absorption_total%E
      allocate(Material(i)%Ph_pair_total%Total(Nsiz))
      allocate(Material(i)%Ph_pair_total%Total_MFP(Nsiz))
      Material(i)%Ph_pair_total%Total(:) = 0.0d0
      Material(i)%Ph_pair_total%Total_MFP(:) = 0.0d0
      ! Sum up CSs from each element:
      do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         do m = 1, Nsiz	! for all energy grid points:
            E => Material(i)%Ph_pair_total%E(m)    ! photon energy [eV]
!             call Find_in_array_monoton(Element%Phot_pair%E,  E, i_closest)	! see below
!             pause
!             print*, 'i_closest', i_closest, E
!             print*, 'X', Element%Phot_pair%E
            
            call interpolate_data_single(Element%Phot_pair%E,  Element%Phot_pair%Total(:), E, sigma) ! module "Little_subroutines"
            call interpolate_data_single(Element%Phot_pair%E,  Element%Phot_pair%Total_MFP(:), E, MFP) ! module "Little_subroutines"
            ! Add them into the arrays:
            Material(i)%Ph_pair_total%Total(m) = Material(i)%Ph_pair_total%Total(m) + sigma*elem_contrib
            Material(i)%Ph_pair_total%Total_MFP(m) = Material(i)%Ph_pair_total%Total_MFP(m) + 1.0d0/MFP !*elem_contrib ! to be inverted
         enddo ! m
      enddo ! j
      ! Invert to get MFP:
      Material(i)%Ph_pair_total%Total_MFP(:) = 1.0d0/Material(i)%Ph_pair_total%Total_MFP(:) ! [A]
      ! Save the total cross section into the file:
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
      if (.not.file_exist .or. numpar%recalculate_MFPs) then	! only create it if file does not exist
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            write(FN,'(es,es,es)') Material(i)%Ph_pair_total%E(m), Material(i)%Ph_pair_total%Total(m), Material(i)%Ph_pair_total%Total_MFP(m)
         enddo ! m = 1, Nsiz
         close(FN)
      endif
      
   enddo TRGT
   
   write(*, '(a)') ' Done.'
!    print*, 'Photon pair creation scattering cross sections are obtained.'
   nullify(path_sep, Element, E)
end subroutine get_photon_pair_creation



! Pair creation scattering cross section per element:
subroutine get_pair_CS(Ee, Element, numpar, sigma)
   real(8), intent(in) :: Ee	! [eV] photon energy
   type(Atom_kind), intent(in) :: Element	! data for this element
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), intent(out) :: sigma	! [A^2] cross section
   select case (numpar%Ph_Pairs)
   case (1)	! PENELOPE
      sigma = Pair_CS(Ee, Element)	! see below
   case default	! exclude
      sigma = 0.0d0
   end select
end subroutine get_pair_CS


function Pair_CS(Ee, Element, Emax_in) result(sigma)
   real(8) :: sigma	! [A^2] cross section
   real(8), intent(in) :: Ee	! [eV] photon energy
   type(Atom_kind), intent(in), target :: Element	! data for this element
   real(8), intent(in), optional :: Emax_in
   !------------------------------------
   real(8) :: dSigma, dSigma0, dSigma_mid, dS, eps, sig_step
   real(8) :: Emin, Emax, k, k_in, Ecur, Ecur0, dE, dE_max, dE_min, dE_half, E_mid
   integer :: i, Ngrid

   if (Ee > 2.0d0*g_me_eV) then	! Photon with such energy can create a (e-, e+) pair
      eps = 1.0d-8	! precision limit
      dS = 0.01d0	! maximal allowed change in dSigma per step
      Ngrid = 100	! grid point for integration over Enorm
      
      k = Ee/g_me_eV
      Emin = 1.0d0/k
      Emax = 1.0d0 - 1.0d0/k
      if (present(Emax_in)) then
         if (Emax_in < Emin) then
            Emax = Emin
         elseif (Emax_in < Emax) then
            Emax = Emax_in
         endif
      endif
      Ecur = Emin	! starting
      dE_max = (Emax - Emin)/dble(Ngrid)	! maximal allowed integration step
      dE_min = dE_max/100.0d0
      dE = dE_max	! start with it, and later reduce if needed
      dE_half = dE*0.50d0
      ! Start integration of differential cross section to obtain the total one:
      sigma = 0.0d0
      dSigma0 = Pair_dSigma(Ee, Ecur, Element%Zat, Element%Pair_R, Element%Pair_nu_inf)	! below
      !(Ee, Element, k, mu)	! below
      i = 0   
      do while (Ecur < Emax)
         i = i + 1	! steps counter
         ! Simpson 3/8 method of integration:
         Ecur0 = Ecur
         Ecur = Ecur + dE
         if (Ecur > Emax) then	! if by chance we exceeded the limit
            Ecur = Emax
            dE = Emax - Ecur0
            dE_half = 0.5d0*dE
         endif
         dSigma = Pair_dSigma(Ee, Ecur, Element%Zat, Element%Pair_R, Element%Pair_nu_inf)	! below
         ! Adaptive step: if it's too large, reduce it:
         if ((dSigma > eps) .and. (dSigma0 > eps)) then	! makes sense to take care of integration
            ! Check if it's too large change per step:
!             print*, 'Sigma', dSigma, dSigma0
            do while ((ABS(dSigma - dSigma0) > dS*dSigma0) .and. (dE > dE_min))
               dE = dE*0.5d0	! reduce timestep
               dE_half = 0.5d0*dE	! adjust half-step correspondingly
               Ecur = Ecur0 + dE
               if (Ecur > Emax) exit
               dSigma = Pair_dSigma(Ee, Ecur, Element%Zat, Element%Pair_R, Element%Pair_nu_inf)	! below
               if (dE < dE_min) exit	! too smal step
            enddo
            ! Check if it's too small change per step (inefficient):
            do while ((ABS(dSigma - dSigma0) < 0.5d0*dS*dSigma0) .and. (dE < dE_max))
               dE = dE + dE_min	! increase timestep
               dE_half = 0.5d0*dE	! adjust half-step correspondingly
               Ecur = Ecur0 + dE
               if (Ecur > Emax) exit
               dSigma = Pair_dSigma(Ee, Ecur, Element%Zat, Element%Pair_R, Element%Pair_nu_inf)	! below
               if (dE > dE_max) exit	! too large step
            enddo
         else	! reset the step
            dE = dE_max
            dE_half = dE*0.5d0
            Ecur = Ecur0 + dE
         endif
         ! Now proceed with the integration:
         if (Ecur > Emax) then	! if by chance we exceeded the limit
            Ecur = Emax
            dE = Emax - Ecur0
            dE_half = dE*0.50d0
         endif
         E_mid =  Ecur0 + dE_half
         dSigma_mid = Pair_dSigma(Ee, E_mid, Element%Zat, Element%Pair_R, Element%Pair_nu_inf)	! below
      
         ! Add up the contributions according to Simpson-3/8 scheme:
         sig_step = dE/6.0d0*(dSigma0 + 4.0d0*dSigma_mid + dSigma)
         sigma = sigma + sig_step
!          write(*,'(a,f,f,f,es,es)') 'S', Ee*1e-6, Ecur, dE, dSigma, sigma
         ! Save the data for the next point of integration:
         dSigma0 = dSigma
      enddo
!       PAUSE 'PAIR'
   else	! such photon energy is insufficient to Compton-scatter
      sigma = 0.0d0
   endif
end function Pair_CS



pure function Pair_dSigma(E, dE, Z, R, nu_inf) result(dSigma)
   real(8) dSigma 	! [A^2]/rad
   real(8), intent(in) :: E, dE	! [eV] photon energy, transferred energy
   integer, intent(in) :: Z	! atomic number
   real(8), intent(in) :: R, nu_inf	! coefficients
   real(8) :: Eps, tempE, Phi1, Phi2, k, nu, b, dblZ
   Eps = dE   !/E
   k = E/g_me_eV
   tempE = 0.5d0 - Eps
   dblZ = dble(Z)
   nu = Pair_nu(k, nu_inf, dblZ)	! below
   b = R/(2.0d0*k*Eps*(1.0d0 - Eps))	! [1] Eq.(2.82) // (2.79) in 2015-edition
   call Pair_Phis(k, b, R, dblZ, Phi1, Phi2)	! below
!    write(*,'(a,f,f,es,f,es, es, es)') 'dS', Phi1, Phi2, b, k, Eps, E, dE
   dSigma = g_r0*g_r0*g_alpha*dblZ*(dblZ + nu)*m_Cr*m_two_third* (2.0d0*tempE*tempE*Phi1 + Phi2)	! [1] Eq.(2.90)
end function Pair_dSigma


pure subroutine Pair_Phis(k, b, R, Z, Phi1, Phi2)
   real(8), intent(in) :: k, b, R, Z
   real(8), intent(out) :: Phi1, Phi2
   real(8) :: CapPhi1, CapPhi2
   real(8) :: g0, g1, g2, lnR, CP3, fc, F0
   call Pair_Capital_Phis(b, R, CapPhi1, CapPhi2)	! below
   CP3 = 3.0d0*CapPhi1
   g1 = 0.5d0*(CP3 - CapPhi2)	! [1] Eq.(2.92)
   g2 = 0.25d0*(CP3 + CapPhi2)	! [1] Eq.(2.92)
   fc = Pair_fc(Z)	! below
   F0 = Pair_F0(k, Z)	! below
   g0 = -4.0d0*fc + F0	! [1] Eq.(2.92)
   Phi1 = g1 + g0	! [1] Eq.(2.91)
   Phi2 = g2 + g0	! [1] Eq.(2.91)
!    write(*,'(a,f,f,f,f,f)') 'P', Phi1, Phi2, g0, g1, g2
   if (Phi1 < 0.0d0) Phi1 = 0.0d0
   if (Phi2 < 0.0d0) Phi2 = 0.0d0
end subroutine Pair_Phis


pure function Pair_fc(Z) result(fc)
   real(8) fc
   real(8), intent(in) :: Z
   real(8) :: a, a2, a4, a6
   a = g_alpha*Z
   a2 = a*a
   a4 = a2*a2
   a6 = a2*a4
   fc = a2*( 1.0d0/(1.0d0 + a2) + 0.202059d0 - 0.0393d0*a2 + 0.00835d0*a4 - 0.00201d0*a6 + 0.00049d0*a4*a4 - 0.00012d0*a6*a4 + 0.00003d0*a6*a6)	! [1] Eq.(2.83)
end function Pair_fc


pure function Pair_F0(k, Z) result(F0)
   real(8) F0
   real(8), intent(in) :: k, Z
   real(8) :: k2, sqrk2, a, a2
   a = g_alpha*Z
   a2 = a*a
   k2 = 2.0d0/k
   sqrk2 = sqrt(k2)
   F0 = (-0.1774d0-12.1d0*a + 11.18*a2)*sqrk2 + (8.523d0 + 73.26d0*a - 44.41d0*a2)*k2 - (13.52d0 + 121.1d0*a - 96.41d0*a2)*k2*sqrk2 + (8.946d0 + 62.05d0*a - 63.41d0*a2)*k2*k2	! [1] Eq.(2.93)
end function Pair_F0



pure subroutine Pair_Capital_Phis(b, R, Phi1, Phi2)
   real(8), intent(in) :: b, R
   real(8), intent(out) :: Phi1, Phi2
   real(8) lnb, artb, ln4, binv, b2
   b2 = b*b
   lnb = 2.0d0*log(1.0d0 + b2)
   binv = 1.0d0/b
   artb = 4.0d0*b*atan(binv)
   ln4 = 4.0d0*log(R)
   Phi1 = 2.0d0 - lnb - artb + ln4	! [1] Eq.(2.81)
   Phi2 = 2.0d0*m_two_third - lnb + 2.0d0*b2*( 4.0d0 - artb - 3.0d0*log(1.0d0 + binv*binv) ) + ln4	! [1] Eq.(2.81)
end subroutine Pair_Capital_Phis



pure function Pair_nu(k, nu_inf, Z) result(nu)	! [1] Eq.(2.88)
   real(8) nu
   real(8), intent(in) :: k, nu_inf, Z
   real(8) :: v
   v = Pair_v(k, Z)	! below
   nu = (1.0d0 - exp(-v))*nu_inf
end function Pair_nu

pure function Pair_v(k, Z) result(v)	! [1] Eq.(2.89)
   real(8) v
   real(8), intent(in) :: k, Z
   real(8) :: k4, a, k42
   k4 = log(4.0d0/k)
   k42 = k4*k4
   a = g_alpha*Z
   v = (0.284d0 - 0.1909d0*a)*k4 + (0.1095d0 + 0.2206d0*a)*k42 + (0.02888d0 - 0.0469d0*a)*k4*k42 + (0.002527d0 + 0.002623d0*a)*k42*k42
end function Pair_v

 
 
end module CS_photons_pair_creation
