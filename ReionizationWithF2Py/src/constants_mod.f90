module constants_mod
  use kinds_mod, only: dp
  implicit none
  private
  public :: c_light, pi, two_pi, kboltz, mprot, hPlanck, hPlanck_eV, &
            Mpcbycm, yrbysec, mp_by_mHe, T_HI, T_HeI, T_HeII, &
            E_HII, E_HeII, E_HeIII, nu_HII, nu_HeII, nu_HeIII, &
            delta_c_z_0, nu_alpha, nu_beta, nu_gamma, nu_delta, Y_He, rho_crit_0, Msun

  ! Mathematical constants
  real(dp), parameter :: pi       = 3.14159265358979324_dp
  real(dp), parameter :: two_pi   = 2.0_dp * pi

  ! Physical constants
  real(dp), parameter :: c_light  = 3.0e10_dp            ! speed of light [cm/s]
  real(dp), parameter :: kboltz   = 1.38e-16_dp          ! Boltzmann constant [erg/K]
  real(dp), parameter :: mprot    = 1.67e-24_dp          ! proton mass [g]
  real(dp), parameter :: hPlanck  = 6.6260755e-27_dp     ! Planck constant [erg.s]
  real(dp), parameter :: hPlanck_eV = 4.1356692e-15_dp   ! Planck constant [eV.s]
  real(dp), parameter :: Msun = 1.989e33_dp              !solar mass in cgs
  real(dp), parameter :: Y_He = 0.2453
  real(dp), parameter :: rho_crit_0 = 2.775e11_dp  ! h^2 Msun / Mpc^3
  real(dp), parameter :: delta_c_z_0 =1.69_dp
  ! Unit conversions
  real(dp), parameter :: Mpcbycm  = 3.0857e24_dp         ! Mpc in cm
  real(dp), parameter :: yrbysec  = 365.0_dp * 24.0_dp * 3600.0_dp  ! year in seconds

  ! Abundance
  real(dp), parameter :: mp_by_mHe = 0.2518_dp           ! proton to Helium mass ratio

  ! Ionization temperatures [K]
  real(dp), parameter :: T_HI    = 157807.0_dp
  real(dp), parameter :: T_HeI   = 285335.0_dp
  real(dp), parameter :: T_HeII  = 631515.0_dp

  ! Ionization energies [eV]
  real(dp), parameter :: E_HII   = 13.5984_dp
  real(dp), parameter :: E_HeII  = 24.5874_dp
  real(dp), parameter :: E_HeIII = 54.416_dp

  ! Ionization frequencies [Hz]
  real(dp), parameter :: nu_HII   = E_HII   / hPlanck_eV
  real(dp), parameter :: nu_HeII  = E_HeII  / hPlanck_eV
  real(dp), parameter :: nu_HeIII = E_HeIII / hPlanck_eV

  ! Lyman series frequencies [Hz]
  real(dp), parameter :: nu_alpha = c_light / 1215.67_dp * 1.0e8_dp
  real(dp), parameter :: nu_beta  = c_light / 1025.72_dp * 1.0e8_dp
  real(dp), parameter :: nu_gamma = c_light /  972.54_dp * 1.0e8_dp
  real(dp), parameter :: nu_delta = c_light /  949.74_dp * 1.0e8_dp


end module constants_mod
