! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains all cross sections for electrons Bremsstrahlung photon emission, etc.
! References used:
! [1] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2008: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2009)
! [2] F. Salvat, J.M. Fernhdez-Varea, NIMB 63 (1992) 255-269
! [3] MCNP — A General Monte Carlo N-Particle Transport Code, Version 5. April 24, 2003 (Revised 10/3/05)


module CS_electrons_Bremsstrahlung
use Universal_constants
use CDF_Ritchi
use Relativity
use Objects
use Read_input_data, only: m_electron_CS, m_input_folder, m_electron_Brems_CS, m_folder_materials
use Dealing_with_files, only: read_file
use CS_photons_pair_creation, only: Pair_Phis, Pair_fc, Pair_Capital_Phis, Pair_nu, Pair_v
use CS_general_tools, only: MFP_from_sigma
use Little_subroutines, only: interpolate_data_single

implicit none

real(8) :: m_1_3, m_Wth

parameter (m_1_3 = 1.0d0/3.0d0)
parameter (m_Wth = 1.0d0)   ! [eV] Threshold for photon energy allowed to be emitted via Bremsstrahlung; [1] Eq.(3.145)

 contains
 

!-----------------------------------------------------------------------------
! Subroutines for Bremsstrahlung:
subroutine get_electron_Brems(Material, numpar, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material	!material parameters of each target that it's constructed of
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------------------------
   integer :: i, j, k, m, Ngrid, N_targets, N_elements, N_shells, FN, Reason, count_lines, Nsiz
   real(8) :: sigma, N_elem, elem_contrib, MFP
   character(200) :: Path, Folder_with_CS, command, File_name, Model_name, Path_total
   logical :: file_exist, read_well
   !--------------------------------
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   
   write(*, '(a)', advance='no') ' Obtaining electron Bremsstrahlung cross sections...'
   
   path_sep => numpar%path_sep
   Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_electron_CS))	! Electron CSs are storred in this folder
   
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
         ! Set part of the file name corresponding to the model used for electron Bremsstrahlung CS:
         select case (numpar%El_Brems)	! elastic scattering: 0=excluded, 1=BHW
         case (1)	! BHW
            Model_name = 'BHW'
         case default	! exclude
            Model_name = 'NO'
         end select
         
         ! For electron CSs use the same grid as for photons:
         Ngrid = size(Element%Phot_absorption%E)
         allocate(Element%El_brems%E(Ngrid))
         Element%El_brems%E = Element%Phot_absorption%E
         allocate(Element%El_brems%Total(Ngrid))
         allocate(Element%El_brems%Total_MFP(Ngrid))
         
         ! Calculate total elastic cross sections for this element:
         File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_electron_Brems_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
         if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! only create it if file does not exist
            open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
            count_lines = 0
            do m = 1, Ngrid	! for all energy grid points:
               read(FN,'(es,es,es)', IOSTAT=Reason) Element%El_brems%E(m), Element%El_brems%Total(m), Element%El_brems%Total_MFP(m)
               call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
               if (.not.read_well) then
                  close(FN)	! redo the file
                  goto 8392
               endif
            enddo ! m = 1, Ngrid
         else 
8392     open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
            do m = 1, Ngrid	! for all energy grid points:
               call get_el_Brems_CS(Element%El_brems%E(m), Element, numpar, Element%El_brems%Total(m))	! below
               Element%El_brems%Total_MFP(m) = MFP_from_sigma(Element%El_brems%Total(m),  1.0d24)    ! module "CS_general_tools"
               write(FN,'(es,es,es)') Element%El_brems%E(m), Element%El_brems%Total(m), Element%El_brems%Total_MFP(m)
            enddo ! m = 1, Ngrid
         endif
         close(FN)
         ! Normalize MFP to the material density:
         ! Take into account the density of atoms of this particular kind:
         elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         Element%El_brems%Total_MFP(:) = Element%El_brems%Total_MFP(:) * 1.0d24/(Material(i)%At_Dens * elem_contrib)
      enddo LMNT
      
      ! And total cross sections:
      ! Where to save:
      Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
      Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
      File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_electron_Brems_CS))//'_total_'//trim(adjustl(Model_name))//'.dat'
      ! Get the cross sections to save:
      Nsiz = size(Material(i)%Ph_absorption_total%E)   ! grid for total cross section is different from core-shells grid
      allocate(Material(i)%El_Brems_total%E(Nsiz))
      Material(i)%El_Brems_total%E = Material(i)%Ph_absorption_total%E
      allocate(Material(i)%El_Brems_total%Total(Nsiz))
      allocate(Material(i)%El_Brems_total%Total_MFP(Nsiz))
      Material(i)%El_Brems_total%Total(:) = 0.0d0
      Material(i)%El_Brems_total%Total_MFP(:) = 0.0d0
       ! Sum up CSs from each element:
      do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         do m = 1, Nsiz	! for all energy grid points:
            E => Material(i)%El_Brems_total%E(m)    ! photon energy [eV]
            call interpolate_data_single(Element%El_brems%E,  Element%El_brems%Total(:), E, sigma) ! module "Little_subroutines"
            call interpolate_data_single(Element%El_brems%E,  Element%El_brems%Total_MFP(:), E, MFP) ! module "Little_subroutines"
            ! Add them into the arrays:
            Material(i)%El_Brems_total%Total(m) = Material(i)%El_Brems_total%Total(m) + sigma*elem_contrib
            Material(i)%El_Brems_total%Total_MFP(m) = Material(i)%El_Brems_total%Total_MFP(m) + 1.0d0/MFP !*elem_contrib ! to be inverted
         enddo ! m
      enddo ! j
      ! Invert to get MFP:
      Material(i)%El_Brems_total%Total_MFP(:) = 1.0d0/Material(i)%El_Brems_total%Total_MFP(:) ! [A]
      ! Save the total cross section into the file:
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
      if (.not.file_exist .or. numpar%recalculate_MFPs) then   ! only create it if file does not exist or user wants to recalculate it
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            write(FN,'(es,es,es)') Material(i)%El_Brems_total%E(m), Material(i)%El_Brems_total%Total(m), Material(i)%El_Brems_total%Total_MFP(m)
         enddo ! m = 1, Nsiz
         close(FN)
      endif
      
   enddo TRGT
   
   write(*, '(a)') ' Done.'
!    print*, 'Electron Bremsstrahlung cross sections are obtained.'
   nullify(path_sep, E, Element)
end subroutine get_electron_Brems


subroutine get_el_Brems_CS(Ee, Element, numpar, sigma, Emax_in)
   real(8), intent(in) :: Ee	! [eV] electron energy
   type(Atom_kind), intent(in), target :: Element	! data for this element
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   real(8), intent(out) :: sigma	! [A^2] cross section of Bremsstrahlung
   real(8), intent(in), optional :: Emax_in ! user provided upper integration limit (used for calculation of transferred energy)
   select case (numpar%El_Brems)	! Bremsstrahlung scattering: 0=excluded, 1=BHW
   case (1)	! BHW
      if (present(Emax_in)) then
         sigma = Bremsstrahlung_total_CS(Ee, Element, Emax_in)	! below
      else
         sigma = Bremsstrahlung_total_CS(Ee, Element)	! below
      endif
   case default	! exclude
      sigma = 0.0d0
   end select
end subroutine get_el_Brems_CS



function Bremsstrahlung_total_CS(Ee, Element, Emax_in) result(sigma)
   real(8) :: sigma	! [A^2] cross section
   real(8), intent(in) :: Ee	! [eV] electron energy
   type(Atom_kind), intent(in), target :: Element	! data for this element
   real(8), intent(in), optional :: Emax_in ! user provided upper integration limit (used for calculation of transferred energy)
   !------------------------------------
   real(8) :: dSigma, dSigma0, dSigma_mid, dS, eps, sig_step, Wth
   real(8) :: Emin, Emax, Ecur, Ecur0, dE, dE_max, dE_min, dE_half, E_mid
   integer :: i, Ngrid

   Wth = m_Wth  ! [eV] Threshold for photon energy allowed to be emitted via Bremsstrahlung; [1] Eq.(3.145)
   if (Ee > 2.0d0*Wth) then	! only for energies larger than threshold
      eps = 1.0d-8	! precision limit
      dS = 0.01d0	! maximal allowed change in dSigma per step
      Ngrid = 100	! grid point for integration over Enorm
      Emin = Wth/(Ee + g_me_eV)
      Emax = Ee/(Ee + g_me_eV)
      Ecur = Emin	! starting
      dE_max = (Emax - Emin)/dble(Ngrid)	! maximal allowed integration step
      ! In a case if user provided the upper integration limit:
      if (present(Emax_in)) then
         ! make sure user did not provide unphysical value:
         if (Emax_in < Emax) Emax = Emax_in
         if (Emax_in < Emin) Emax = Emin
      endif
      dE_min = dE_max/100.0d0
      dE = dE_max	! start with it, and later reduce if needed
      dE_half = dE*0.50d0
      ! Start integration of differential cross section to obtain the total one:
      sigma = 0.0d0
      dSigma0 = Brems_dSigma(Ee, Ecur, Element%Zat, Element%Pair_R, Element%Pair_nu_inf)	! below
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
         dSigma = Brems_dSigma(Ee, Ecur, Element%Zat, Element%Pair_R, Element%Pair_nu_inf)	! below
         ! Adaptive step: if it's too large, reduce it:
         if ((dSigma > eps) .and. (dSigma0 > eps)) then	! makes sense to take care of integration
            ! Check if it's too large change per step:
!             print*, 'Sigma', dSigma, dSigma0
            do while ((ABS(dSigma - dSigma0) > dS*dSigma0) .and. (dE > dE_min))
               dE = dE*0.5d0	! reduce timestep
               dE_half = 0.5d0*dE	! adjust half-step correspondingly
               Ecur = Ecur0 + dE
               if (Ecur > Emax) exit
               dSigma = Brems_dSigma(Ee, Ecur, Element%Zat, Element%Pair_R, Element%Pair_nu_inf)	! below
               if (dE < dE_min) exit	! too smal step
            enddo
            ! Check if it's too small change per step (inefficient):
            do while ((ABS(dSigma - dSigma0) < 0.5d0*dS*dSigma0) .and. (dE < dE_max))
               dE = dE + dE_min	! increase timestep
               dE_half = 0.5d0*dE	! adjust half-step correspondingly
               Ecur = Ecur0 + dE
               if (Ecur > Emax) exit
               dSigma = Brems_dSigma(Ee, Ecur, Element%Zat, Element%Pair_R, Element%Pair_nu_inf)	! below
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
         dSigma_mid = Brems_dSigma(Ee, E_mid, Element%Zat, Element%Pair_R, Element%Pair_nu_inf)	! below
      
         ! Add up the contributions according to Simpson-3/8 scheme:
         sig_step = dE/6.0d0*(dSigma0 + 4.0d0*dSigma_mid + dSigma)
         sigma = sigma + sig_step
!       write(*,'(a,f,f,f,es,es)') 'S', Ee*1e-6, Ecur, dE, dSigma, sigma
         ! Save the data for the next point of integration:
         dSigma0 = dSigma
      enddo
   else	! below threshold, no Bremsstrahlung
      sigma = 0.0d0
   endif
end function Bremsstrahlung_total_CS


function  Brems_dSigma(Ee, Ecur, Z, R, nu_inf) result(dSigma)
   real(8) dSigma	! [A^2]/rad
   real(8), intent(in) :: Ee, Ecur	! [eV] photon energy, transferred energy
   integer, intent(in) :: Z	! atomic number
   real(8), intent(in) :: R, nu_inf	! coefficients
   real(8) :: Phi1, Phi2, k, nu, b, dblZ, gamm
   k = Ee/g_me_eV
   dblZ = dble(Z)
   nu = Pair_nu(k, nu_inf, dblZ)	! module "CS_photons_pair_creation"
   gamm = (Ee+g_me_eV)/g_me_eV
   b = R*Ecur/(2.0d0*gamm*(1.0d0 - Ecur))	! [1] Eq.(3.133)
   call Brems_Phis(Ee, Ecur, b, R, dblZ, Phi1, Phi2)	! below
!    dSigma = g_r0*g_r0*g_alpha*dblZ*(dblZ + nu)/(Ee+g_me_eV) * (Ecur*Phi1 + 4.0d0/(3.0d0*Ecur)*(1.0d0 - Ecur)*Phi2)	! [1] Eq.(3.132), or [2] Eq.(49)
   dSigma = g_r0*g_r0*g_alpha*dblZ*(dblZ + nu) * (Ecur*Phi1 + 4.0d0/(3.0d0*Ecur)*(1.0d0 - Ecur)*Phi2)	! [1] Eq.(3.132), or [2] Eq.(49)
   if (dSigma < 0.0d0) dSigma = 0.0d0
end function  Brems_dSigma


subroutine Brems_Phis(Ee, Ecur, b, R, Z, Phi1, Phi2)
   real(8), intent(in) :: Ee, Ecur, b, R, Z
   real(8), intent(out) :: Phi1, Phi2
   real(8) :: CapPhi1, CapPhi2
   real(8) :: f0, f1, f2, fc, F2ZE
!    print*, 'Brems_Phis', Ee, Ecur, R, b
   call Pair_Capital_Phis(b, R, CapPhi1, CapPhi2)	! module "CS_photons_pair_creation"
   f1 = CapPhi1		! [2] Eq.(52)
   f2 = 0.5d0*(3.0d0*CapPhi1 - CapPhi2)	! [2] Eq.(53)
   if (Ecur >= (Ee - 5.0d0*g_me_eV)/(Ee + g_me_eV) ) then	! [2] Eq.(40)
      fc = 0.0d0
   else
      fc = Pair_fc(Z)	! module "CS_photons_pair_creation"
   endif
   F2ZE = Brems_F2(Ee, Z)	! below
   f0 = -4.0d0*fc + F2ZE	! [2] Eq.(51)
   Phi1 = f1 + f0		! [2] Eq.(50)
   Phi2 = f2 + f0		! [2] Eq.(50)
!    if (Phi1 < 0.0d0) Phi1 = 0.0d0
!    if (Phi2 < 0.0d0) Phi2 = 0.0d0
end subroutine Brems_Phis


pure function Brems_F2(E, Z) result(F2)
   real(8) F2
   real(8), intent(in) :: E, Z
   real(8) :: a
   a = g_alpha*Z
   F2 = (2.04d0 + 9.09d0*a)*(g_me_eV*g_me_eV/(E*(E + g_me_eV)))**(1.26d0 - 0.93d0*a)	! [2] Eq.(54)
end function Brems_F2



function sample_Bremsstrahlung_theta(beta, A_in, B_in) result(theta)
   real(8) theta    ! photoemission angle in Bremsstrahlung, according to Eq.(3.175) [1]
   real(8), intent(in) :: beta    ! relative velosity beta=v/c
   real(8), intent(in), optional :: A_in, B_in    ! coefficients A, B
   real(8) :: A, B, RN, C, D, mu, mu_prime, beta_prime, A2, A32, E
   
   if (present(A_in)) then
      A = A_in
   else ! use default value, paragraph after Eq.(3.175) [1]
      A = 1.0d0
   endif
   
   if (present(B_in)) then
      B = B_in
   else ! use default value, paragraph after Eq.(3.175) [1]
      B = 0.0d0
   endif
   
   ! Sample random number:
   call random_number(RN)

   ! Analytical solution of Eq.(3.172) (not in [1], solved by me):
   A2 = A*A
   D = A*A2 - 6.0d0*A2 + 32.0d0*RN*(RN-1.0d0) + 48.0d0*A*RN*(1.0d0 - RN)
   A32 = 3.0d0*A - 2.0d0
   E = (sqrt(-D/A32) + 4.0d0*RN - 2)*A32*A32 
   C = E**m_1_3
   ! Combine it all into the analytical solution:
   mu_prime = C/A32 + (A-2.0d0)/C
   
   beta_prime = beta*(1.0d0 + B)    ! after Eq.(3.175)
   
   ! Convert into solution of EQ.(3.175), following point (iv) in section 3.3.4.2 [1]
   mu = (mu_prime + beta_prime)/(1.0d0 + mu_prime*beta_prime)
   
   theta = ACOS(mu)
end function sample_Bremsstrahlung_theta



function get_energy_transfer_Bremsstrahlung(Ee, Element, numpar, j, CS_tot) result(dE)
   real(8) :: dE   ! [eV] sampled transferred energy
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Atom_kind), intent(in) :: Element    ! all information about this element
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), intent(in) :: CS_tot    ! [A^2] precalculated total cross section
   integer, intent(in) :: j     ! number or element
   !---------------------------------------
   real(8) :: eps, RN, CS_sampled, CS_cur, E_left, E_right, E_cur
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-3
   
   E_left = m_Wth ! [eV] minimal photon energy
   E_right = Ee/(Ee + g_me_eV)   ! [eV] maximal transferred energy

   ! Sample the cross section:
   call random_number(RN)
   CS_sampled = RN*CS_tot
   
   ! Start finding CS:
   E_cur = (E_left + E_right)*0.5d0
   CS_cur = Bremsstrahlung_total_CS(Ee, Element, E_cur)  ! above
   
   if (CS_cur < 1.0d-15) then   ! in case CS is zero, Bremsstrahlung is impossible, so, how are we here??
      print*, 'Error in get_energy_transfer_Bremsstrahlung', Ee, CS_cur, CS_sampled
   endif
   
   ! Search by bisection method:
   do while (ABS(CS_cur - CS_sampled) > eps*CS_sampled)
      if (CS_cur > CS_sampled) then
         E_right = E_cur
      else
         E_left = E_cur
      endif
      E_cur = (E_left + E_right)/2.0d0
      if (abs(E_left - E_right) < eps) exit  ! precise enough
      CS_cur =  Bremsstrahlung_total_CS(Ee, Element, E_cur) ! below
   enddo
   
   ! Output: sampled transferred energy:
   dE = E_cur
end function get_energy_transfer_Bremsstrahlung


pure subroutine sample_Bremsstrahlung_electron_angles(phi, theta)
   real(8), intent(out) :: phi, theta
   ! Currently, we assume no electron deflection during Bremsstrahlung photon emission [3],
   ! Chapter 2, section 7. Bremsstrahlung (Page 112 in pdf)
   phi = 0.0d0
   theta = 0.0d0
end subroutine sample_Bremsstrahlung_electron_angles



end module CS_electrons_Bremsstrahlung
