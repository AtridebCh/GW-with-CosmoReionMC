module recrates_mod
  use kinds_mod,     only: dp
  use constants_mod, only: kboltz, T_HI, T_HeI, T_HeII
  implicit none
  private
  public :: R_HII_e_A,   RC_HII_A,   &
            R_HeII_e_A,  RC_HeII_A,  &
            R_HeIII_e_A, RC_HeIII_A, &
            R_HII_e_B,   RC_HII_B,   &
            R_HeII_e_B,  RC_HeII_B,  &
            R_HeIII_e_B, RC_HeIII_B, &
            CI_HI,  CC_HI,           &
            CI_HeI, CC_HeI,          &
            CI_HeII, CC_HeII,        &
            DI_HeII, DC_HeII,        &
            EC_HI, EC_HeII

contains

  ! ---------------------------------------------------------------------------
  ! Private helper: reduced temperature variables
  ! ---------------------------------------------------------------------------
  pure function lambda(T_ion, T) result(res)
    real(dp), intent(in) :: T_ion, T
    real(dp) :: res

    res = 2.0_dp * T_ion / T
  end function lambda


  ! ---------------------------------------------------------------------------
  ! Case A recombination rates [cm^3 s^-1]
  ! NOTE: functions 1-6 are currently zeroed out (switched off)
  ! ---------------------------------------------------------------------------

  !!! 1 - HII case A recombination
  pure function R_HII_e_A(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HI, T)
    res = 1.269e-13_dp * lam**1.503_dp / (1.0_dp + (lam / 0.522_dp)**0.470_dp)**1.923_dp
    res = 0.0_dp  ! switched off
  end function R_HII_e_A

  !!! 2 - HII case A cooling rate [erg cm^3 s^-1]
  pure function RC_HII_A(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HI, T)
    res = 1.778e-29_dp * T * lam**1.965_dp / (1.0_dp + (lam / 0.541_dp)**0.502_dp)**2.697_dp
    res = 0.0_dp  ! switched off
  end function RC_HII_A

  !!! 3 - HeII case A recombination
  pure function R_HeII_e_A(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HeI, T)
    res = 3.0e-14_dp * lam**0.654_dp
    res = 0.0_dp  ! switched off
  end function R_HeII_e_A

  !!! 4 - HeII case A cooling rate
  pure function RC_HeII_A(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res

    res = kboltz * T * R_HeII_e_A(T)
    res = 0.0_dp  ! switched off
  end function RC_HeII_A

  !!! 5 - HeIII case A recombination
  pure function R_HeIII_e_A(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HeII, T)
    res = 2.0_dp * 1.269e-13_dp * lam**1.503_dp / (1.0_dp + (lam / 0.522_dp)**0.470_dp)**1.923_dp
    res = 0.0_dp  ! switched off
  end function R_HeIII_e_A

  !!! 6 - HeIII case A cooling rate
  pure function RC_HeIII_A(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HeII, T)
    res = 8.0_dp * 1.778e-29_dp * T * lam**1.965_dp / (1.0_dp + (lam / 0.541_dp)**0.502_dp)**2.697_dp
    res = 0.0_dp  ! switched off
  end function RC_HeIII_A


  ! ---------------------------------------------------------------------------
  ! Case B recombination rates [cm^3 s^-1]
  ! ---------------------------------------------------------------------------

  !!! 7 - HII case B recombination
  pure function R_HII_e_B(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HI, T)
    res = 2.753e-14_dp * lam**1.5_dp / (1.0_dp + (lam / 2.74_dp)**0.407_dp)**2.242_dp
  end function R_HII_e_B

  !!! 8 - HII case B cooling rate
  pure function RC_HII_B(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HI, T)
    res = 3.435e-30_dp * T * lam**1.97_dp / (1.0_dp + (lam / 2.25_dp)**0.376_dp)**3.72_dp
  end function RC_HII_B

  !!! 9 - HeII case B recombination
  pure function R_HeII_e_B(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HeI, T)
    res = 1.26e-14_dp * lam**0.75_dp
  end function R_HeII_e_B

  !!! 10 - HeII case B cooling rate
  pure function RC_HeII_B(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res

    res = kboltz * T * R_HeII_e_B(T)
  end function RC_HeII_B

  !!! 11 - HeIII case B recombination
  pure function R_HeIII_e_B(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HeII, T)
    res = 2.0_dp * 2.753e-14_dp * lam**1.5_dp / (1.0_dp + (lam / 2.74_dp)**0.407_dp)**2.242_dp
  end function R_HeIII_e_B

  !!! 12 - HeIII case B cooling rate
  pure function RC_HeIII_B(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HeII, T)
    res = 8.0_dp * 3.435e-30_dp * T * lam**1.97_dp / (1.0_dp + (lam / 2.25_dp)**0.376_dp)**3.72_dp
  end function RC_HeIII_B


  ! ---------------------------------------------------------------------------
  ! Collisional ionization rates [cm^3 s^-1]
  ! ---------------------------------------------------------------------------

  !!! 13 - HI collisional ionization
  pure function CI_HI(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HI, T)
    res = 21.11_dp * T**(-1.5_dp) * exp(-0.5_dp * lam) * &
          lam**(-1.089_dp) / (1.0_dp + (lam / 0.354_dp)**0.874_dp)**1.101_dp
  end function CI_HI

  !!! 14 - HI collisional ionization cooling
  pure function CC_HI(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res

    res = kboltz * T_HI * CI_HI(T)
  end function CC_HI

  !!! 15 - HeI collisional ionization
  pure function CI_HeI(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HeI, T)
    res = 32.38_dp * T**(-1.5_dp) * exp(-0.5_dp * lam) * &
          lam**(-1.146_dp) / (1.0_dp + (lam / 0.416_dp)**0.987_dp)**1.056_dp
  end function CI_HeI

  !!! 16 - HeI collisional ionization cooling
  pure function CC_HeI(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res

    res = kboltz * T_HeI * CI_HeI(T)
  end function CC_HeI

  !!! 17 - HeII collisional ionization
  pure function CI_HeII(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HeII, T)
    res = 19.95_dp * T**(-1.5_dp) * exp(-0.5_dp * lam) * &
          lam**(-1.089_dp) / (1.0_dp + (lam / 0.553_dp)**0.735_dp)**1.275_dp
  end function CI_HeII

  !!! 18 - HeII collisional ionization cooling
  pure function CC_HeII(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res

    res = kboltz * T_HeII * CI_HeII(T)
  end function CC_HeII


  ! ---------------------------------------------------------------------------
  ! Dielectric recombination rates
  ! ---------------------------------------------------------------------------

  !!! 19 - HeII dielectric recombination
  pure function DI_HeII(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HeII, T)
    res = 1.9e-3_dp * T**(-1.5_dp) * exp(-0.5_dp * 0.75_dp * lam) * &
          (1.0_dp + 0.3_dp * exp(-0.5_dp * 0.15_dp * lam))
  end function DI_HeII

  !!! 20 - HeII dielectric recombination cooling
  pure function DC_HeII(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res

    res = 0.75_dp * kboltz * T_HeII * DI_HeII(T)
  end function DC_HeII


  ! ---------------------------------------------------------------------------
  ! Excitation cooling rates [erg cm^3 s^-1]
  ! ---------------------------------------------------------------------------

  !!! 21 - HI excitation cooling
  pure function EC_HI(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HI, T)
    res = 7.5e-19_dp * exp(-0.75_dp * 0.5_dp * lam) / (1.0_dp + sqrt(T / 1.0e5_dp))
  end function EC_HI

  !!! 22 - HeII excitation cooling
  pure function EC_HeII(T) result(res)
    real(dp), intent(in) :: T
    real(dp) :: res, lam

    lam = lambda(T_HeII, T)
    res = 5.54e-17_dp * T**(-0.397_dp) * exp(-0.75_dp * 0.5_dp * lam) / &
          (1.0_dp + sqrt(T / 1.0e5_dp))
  end function EC_HeII

end module recrates_mod
