module crossSection_mod
  use kinds_mod,     only: dp
  use constants_mod, only: hPlanck_eV
  use adaptint_mod,  only: d01arf
  implicit none
  private
  public :: sigma_HI, sigma_HeI, sigma_HeII, &
            sigma_integral_HI, sigma_integral_HeI, sigma_integral_HeII

  ! shared working variable for integrand callbacks
  real(dp) :: dumindex

contains

  ! ---------------------------------------------------------------------------
  ! Photoionization cross sections
  ! ---------------------------------------------------------------------------
  function sigma_HI(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res
    real(dp), parameter :: sigma_0 = 5.475e-14_dp
    real(dp), parameter :: P       = 2.963_dp
    real(dp) :: nu_0, x

    nu_0 = 4.298e-1_dp / hPlanck_eV
    x    = nu / nu_0
    res  = sigma_0 * (x - 1.0_dp)**2 * x**(0.5_dp * P - 5.5_dp) / &
           (1.0_dp + sqrt(x / 32.88_dp))**P
  end function sigma_HI


  function sigma_HeI(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res
    real(dp), parameter :: sigma_0 = 9.492e-16_dp
    real(dp), parameter :: P       = 3.188_dp
    real(dp), parameter :: yw      = 2.039_dp
    real(dp), parameter :: y0      = 0.4434_dp
    real(dp), parameter :: y1      = 2.136_dp
    real(dp) :: nu_0, x, yy

    nu_0 = 1.361e1_dp / hPlanck_eV
    x    = nu / nu_0 - y0
    yy   = sqrt(x * x + y1 * y1)
    res  = sigma_0 * (yw**2 + (x - 1.0_dp)**2) * yy**(0.5_dp * P - 5.5_dp) / &
           (1.0_dp + sqrt(yy / 1.469_dp))**P
  end function sigma_HeI


  function sigma_HeII(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res
    real(dp), parameter :: sigma_0 = 1.369e-14_dp
    real(dp), parameter :: P       = 2.963_dp
    real(dp) :: nu_0, x

    nu_0 = 1.72_dp / hPlanck_eV
    x    = nu / nu_0
    res  = sigma_0 * (x - 1.0_dp)**2 * x**(0.5_dp * P - 5.5_dp) / &
           (1.0_dp + sqrt(x / 32.88_dp))**P
  end function sigma_HeII


  ! ---------------------------------------------------------------------------
  ! Photoionization cross section integrals: integral of nu^index * sigma(nu)
  ! Uses a shared integrator to avoid code duplication
  ! ---------------------------------------------------------------------------
  function sigma_integral_HI(index, numin, numax) result(res)
    real(dp), intent(in) :: index, numin, numax
    real(dp) :: res

    call integrate_sigma(sigma_integrand_HI, index, numin, numax, res)
  end function sigma_integral_HI


  function sigma_integral_HeI(index, numin, numax) result(res)
    real(dp), intent(in) :: index, numin, numax
    real(dp) :: res

    call integrate_sigma(sigma_integrand_HeI, index, numin, numax, res)
  end function sigma_integral_HeI


  function sigma_integral_HeII(index, numin, numax) result(res)
    real(dp), intent(in) :: index, numin, numax
    real(dp) :: res

    call integrate_sigma(sigma_integrand_HeII, index, numin, numax, res)
  end function sigma_integral_HeII


  ! ---------------------------------------------------------------------------
  ! Private helpers
  ! ---------------------------------------------------------------------------
  subroutine integrate_sigma(integrand, index, numin, numax, res)
    real(dp), intent(in)  :: index, numin, numax
    real(dp), intent(out) :: res
    integer  :: ifail, maxrul, nn, iparm
    real(dp) :: epsabs, epsrel, abserr
    real(dp), dimension(390) :: alpha
    interface
      function integrand(nu) result(val)
        import dp
        real(dp), intent(in) :: nu
        real(dp) :: val
      end function integrand
    end interface

    dumindex = index
    maxrul   = 9
    iparm    = 0
    epsabs   = 0.0_dp
    epsrel   = 1.0e-8_dp
    ifail    = -1
    call d01arf(numin, numax, integrand, epsrel, epsabs, maxrul, iparm, &
                abserr, res, nn, alpha, ifail)
  end subroutine integrate_sigma


  function sigma_integrand_HI(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res

    res = nu**dumindex * sigma_HI(nu)
  end function sigma_integrand_HI


  function sigma_integrand_HeI(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res

    res = nu**dumindex * sigma_HeI(nu)
  end function sigma_integrand_HeI


  function sigma_integrand_HeII(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res

    res = nu**dumindex * sigma_HeII(nu)
  end function sigma_integrand_HeII

end module crossSection_mod
