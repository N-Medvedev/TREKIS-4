! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains all subroutines for dealing with relativistic particles
module Relativity
use Universal_constants

implicit none

 contains


pure function velosity_from_kinetic_energy(Ekin, M0, afs) result(v)
   real(8) v        ! velosity  [A/fs], unless afs=.true. then [m/s]
   real(8), intent(in) :: Ekin	! kinetic energy [eV]
   real(8), intent(in) :: M0	! rest mass [kg]
   logical, intent(in), optional :: afs ! output velosity in [A/fs] (true), or in [m/s] (false)
   real(8) :: fact, Erest
   if (M0 < 1.0d-10*g_me) then   ! assume massless particle
      v = g_cvel    ! [m/s]
   else
      if (Ekin < 1.0d-15) then
         v = 0.0d0
      else
         Erest = rest_energy(M0)	! [eV] rest energy
         fact = Ekin/Erest + 1.0d0
         v = g_cvel*sqrt(1.0d0 - 1.0d0/(fact*fact))	! [m/s]
      endif
   endif
   if (present(afs)) then
      if (afs) v = v * g_ms2Afs ! [m/s] -> [A/fs]
   else ! by default, use [A/fs]
      v = v * g_ms2Afs ! [m/s] -> [A/fs]
   endif
end function velosity_from_kinetic_energy



pure function kinetic_energy_from_momentum(p, M) result(Ekin)
   real(8) Ekin	! kinetic energy [eV]
   real(8), intent(in) :: p	! [kg*m/s] velosity
   real(8), intent(in) :: M	! [kg] rest mass of the particle
   real(8) :: Erest, Etot
   Erest = rest_energy(M)	! see below
   Etot = total_energy_from_momentum(p, M)	! see below
   Ekin = Etot - Erest	! [eV]
end function kinetic_energy_from_momentum



pure function total_energy_from_momentum(p, M) result(Etot)
   real(8) Etot	! total energy [eV]
   real(8), intent(in) :: p	! [kg*m/s] velosity
   real(8), intent(in) :: M	! [kg] rest mass of the particle
   real(8) :: Erest
   Erest = rest_energy(M)	! see below
   Etot = sqrt(p*p*g_cvel*g_cvel/(g_e*g_e) + Erest*Erest)	! [eV]
end function total_energy_from_momentum


pure function momentum_from_velosity(v, M0) result(p)
   real(8) p	! [kg*m/s] momentum
   real(8), intent(in) :: v	! [m/s] velosity
   real(8), intent(in) :: M0	! [kg] rest mass of the particle
   real(8) :: gamma
   gamma = gamma_factor(v)	! see below
   p = gamma*M0*v	! [kg*m/s]
end function momentum_from_velosity


pure function momentum_from_energy(E) result(p)	! (also applicable to photons)
   real(8) p	! [kg*m/s] momentum
   real(8), intent(in) :: E	! [eV] total energy
   p = E/g_cvel/g_e	! [kg*m/s]
end function momentum_from_energy


pure function momentum_from_kinetic_energy(Ekin, M0) result(p)	! (also applicable to photons)
   real(8) p	! [kg*m/s] momentum
   real(8), intent(in) :: Ekin	! [eV] kinetic energy
   real(8), intent(in) :: M0	! [kg] rest mass
   real(8) :: Erest
   Erest = rest_energy(M0)	! see below
   p = Ekin/g_cvel*sqrt(1.0d0 + 2.0d0*Erest/Ekin)*g_e	! [kg*m/s]
end function momentum_from_kinetic_energy


pure function kinetic_energy_from_velosity(v, M0) result(Ekin)
   real(8) Ekin	! kinetic energy [eV]
   real(8), intent(in) :: v     ! [m/s] velosity
   real(8), intent(in) :: M0    ! [kg] rest mass of the particle
   real(8) :: gamma, Etot
   gamma = gamma_factor(v)	! see below
   Etot = rest_energy(M0)	! see below
   Ekin = Etot*(gamma - 1.0d0)	! [eV]
end function kinetic_energy_from_velosity


pure function relativistic_mass_from_energy(E) result(M)
   real(8) M		! [kg] mass of the particle
   real(8), intent(in) :: E	! [eV] total energy
   M = E/(g_cvel*g_cvel)/g_e	! [kg]
end function relativistic_mass_from_energy



pure function relativistic_mass_from_velosity(v, M0) result(M)
   real(8) M		! [kg] mass of the particle
   real(8), intent(in) :: v	! [m/s] velosity
   real(8), intent(in) :: M0	! [kg] rest mass of the particle
   real(8) :: gamma
   gamma = gamma_factor(v)	! see below
   M = gamma*M0	! [kg]
end function relativistic_mass_from_velosity

 
pure function rest_energy(M0) result(E)
   real(8) E	! [eV] total energy
   real(8), intent(in) :: M0	! [kg] rest mass of the particle
   E = M0*g_cvel*g_cvel/g_e	! [eV]
end function rest_energy


 
pure function gamma_factor(v, afs_in) result(gamma)
   real(8) gamma
   real(8), intent(in) :: v	! [m/s] velosity, or in [A/fs] if afs_in = .true.
   logical, intent(in), optional :: afs_in ! v provided in [A/fs] if true, [m/s] otherwise
   real(8) :: beta
   logical :: afs
   if (present(afs_in)) then
      afs = afs_in
   else
      afs = .false. ! by default, v is in [m/s]
   endif
   beta = beta_factor(v, afs)	! see below
   gamma = 1.0d0/dsqrt(1.0d0 - beta*beta)
end function gamma_factor


pure function gamma_factor_from_Ekin(Ekin, M) result(gamma)
   real(8) gamma
   real(8), intent(in) :: Ekin, M   ! kinetic energy [eV], particle mass [kg]
   real(8) :: Mc2
   ! Rest energy [eV]:
   Mc2 = rest_energy(M)  ! above
   gamma = (Ekin + Mc2)/Mc2
end function gamma_factor_from_Ekin

 
pure function beta_factor(v, afs_in) result(beta)
   real(8) beta
   real(8), intent(in) :: v	! [m/s] velosity
   logical, intent(in), optional :: afs_in ! v provided in [A/fs] if true, [m/s] otherwise
   real(8) :: conv
   if (present(afs_in)) then
      if (afs_in) then
         conv = g_Afs2ms ! [A/fs] -> [m/s]
      else
         conv = 1.0d0 ! v is in [m/s] by default
      endif
   else
      conv = 1.0d0  ! v is in [m/s] by default
   endif
   beta = v/g_cvel * conv
end function beta_factor
 
 
end module Relativity
