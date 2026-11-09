module backgroundCosmology_mod
  use kinds_mod,      only: dp
  use constants_mod,  only: pi, two_pi, kboltz, mprot, yrbysec, c_light
  use parameters_mod, only: h, omega_m, omega_l, omega_r, &
                           omega_k,omega_b, ombh2, ns, sigma_8, &
                           rho_c, gamma, dn_dlnk
  use variables_mod,  only: delta_c, massfunc_name, logm, logsig, &
                            logx_arr, dlogsig_dlogx, coeffspl, coeffspl2, coeffspl_deriv
  use adaptint_mod ! got d01amf, d01arf
  use onedspline_mod, only: spline, splevl
  implicit none
  private
  

  public :: sigma, logsigma, dfdnu_PS, numdenm, massfrac,          &
            nu_parameter, tdyn, delvir, pspec, age, angdist, s_k,  &
            comoving, hubbledist, norm, lumdist, omega_z, d,       &
            sigmasq, sigmasq_b, xbsq, f_integrand, set_sigma8_norm,&
            sigma8_norm, generic_dndM, dfdnu_ST,                   &
            differential_comoving_volume, setspline_sigma

  ! module-level working variables (replacing implicit globals)
  real(dp), save, protected :: sigma8_norm   ! add here
  real(dp) :: r      ! used in sigmasq/integrand_sigmasq
  real(dp), private :: cb_x, cb_z, cb_T, cb_mu

contains


  ! ---------------------------------------------------------------------------
  ! Transfer function (Bardeen et al.)
  ! ---------------------------------------------------------------------------
  function transfun(q) result(res)
    real(dp), intent(in) :: q
    real(dp) :: res
    real(dp), parameter :: aa = 6.4_dp, bb = 3.0_dp, cc = 1.7_dp, pow = 1.13_dp

    res = 1.0_dp / ((1.0_dp + (aa*q + (bb*q)**1.5_dp + (cc*q)**2)**pow)**(1.0_dp / pow))
  end function transfun


  ! ---------------------------------------------------------------------------
  ! sigma8 normalisation
  ! ---------------------------------------------------------------------------
  
  subroutine set_sigma8_norm()
    sigma8_norm = norm(sigma_8)   ! sigma_8 already available from parameters_mod
  end subroutine set_sigma8_norm
  
  function norm(sig) result(res)
    real(dp), intent(in) :: sig
    real(dp) :: res
    integer,  parameter :: lw = 2000, liw = lw / 4
    integer  :: inf, ifail
    real(dp) :: epsabs, epsrel, abserr
    integer,  dimension(liw) :: iw
    real(dp), dimension(lw)  :: w

    inf    = 1
    epsabs = 0.0_dp
    epsrel = 1.0e-10_dp
    ifail  = -1
    call d01amf(pswin, 0.0_dp, inf, epsabs, epsrel, res, abserr, w, lw, iw, liw, ifail)
    res = sig * sig * two_pi * pi / res
  end function norm

  function pswin(k) result(res)
    real(dp), intent(in) :: k
    real(dp) :: res
    real(dp) :: r_8

    r_8 = 8.0_dp / h
    res = pspec(k) * k * k * window(k * r_8)
  end function pswin

  function window(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res

    res = 3.0_dp * (sin(x) / x**3 - cos(x) / x**2)
    res = res * res
  end function window


  ! ---------------------------------------------------------------------------
  ! Press-Schechter mass function
  ! ---------------------------------------------------------------------------
  function sigma(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res

    res = sqrt(sigmasq(x))
  end function sigma

  function logsigma(logx) result(res)
    real(dp), intent(in) :: logx
    real(dp) :: res

    res = log(sigma(exp(logx)))
  end function logsigma

  function probdist_PS(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res

    res = sqrt(2.0_dp / pi) * nu * exp(-nu * nu / 2.0_dp)
  end function probdist_PS
  
  function dfridr(func, x, step, err) result(res)
    real(dp), intent(in)  :: x, step
    real(dp), intent(out) :: err
    real(dp) :: res
    integer,  parameter :: NTAB = 10
    real(dp), parameter :: CON = 1.4_dp, CON2 = CON * CON, BIG = 1.0e30_dp, SAFE = 2.0_dp
    integer  :: ierrmin, i, j
    integer,  dimension(1) :: imin
    real(dp) :: hh
    real(dp), dimension(NTAB - 1)    :: errt, fac
    real(dp), dimension(NTAB, NTAB)  :: a
    interface
      function func(x) result(res)
        import dp
        real(dp), intent(in) :: x
        real(dp) :: res
      end function func
    end interface

    if (step == 0.0_dp) stop 'step must be nonzero in dfridr'
    hh      = step
    a(1, 1) = (func(x + hh) - func(x - hh)) / (2.0_dp * hh)
    err     = BIG
    fac(1)  = CON
    do i = 2, NTAB - 1
      fac(i) = fac(i - 1) * CON
    end do
    do i = 2, NTAB
      hh       = hh / CON
      a(1, i)  = (func(x + hh) - func(x - hh)) / (2.0_dp * hh)
      do j = 2, i
        a(j, i) = (a(j-1, i) * fac(j-1) - a(j-1, i-1)) / (fac(j-1) - 1.0_dp)
      end do
      errt(1:i-1) = max(abs(a(2:i, i) - a(1:i-1, i)), abs(a(2:i, i) - a(1:i-1, i-1)))
      imin        = minloc(errt(1:i-1))
      ierrmin     = imin(1)
      if (errt(ierrmin) <= err) then
        err = errt(ierrmin)
        res = a(1 + ierrmin, i)
      end if
      if (abs(a(i, i) - a(i-1, i-1)) >= SAFE * err) return
    end do
  end function dfridr


  function numdenx(x, z) result(res)
    real(dp), intent(in) :: x, z
    real(dp) :: res
    real(dp) :: sig, dlogsigmadlogx, step, err, nu

    sig  = sigma(x)
    step = 0.1_dp
    do
      dlogsigmadlogx = dfridr(logsigma, log(x), step, err)
      if (abs(err) < 1.0e-5_dp) exit
      step = step / 2.0_dp
    end do
    nu  = delta_c / (d(z) * sig)
    res = -probdist_PS(nu) * (3.0_dp / (4.0_dp * pi * x**3)) * dlogsigmadlogx / x
  end function numdenx

  function d(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res, dinv

    if (omega_l == 0.0_dp .and. omega_m == 1.0_dp) then
      dinv = 1.0_dp + z
    else if (omega_l == 0.0_dp .and. omega_m /= 1.0_dp) then
      dinv = 1.0_dp + 2.5_dp * omega_m * z / (1.0_dp + 1.5_dp * omega_m)
    else
      dinv = (1.0_dp + ((1.0_dp + z)**3 - 1.0_dp) / &
              (1.0_dp + 0.4545_dp * (omega_l / omega_m)))**(1.0_dp / 3.0_dp)
    end if
    res = 1.0_dp / dinv
  end function d

  function numdenm(m, z) result(res)
    real(dp), intent(in) :: m, z
    real(dp) :: res, x, rho_0

    rho_0 = rho_c * omega_m
    x     = (3.0_dp * m / (4.0_dp * pi * rho_0))**(1.0_dp / 3.0_dp)
    res   = numdenx(x, z) / (4.0_dp * pi * x * x * rho_0)
  end function numdenm
  
  function generic_dndM(m, z) result(res)
    real(dp), intent(in):: m, z
    real(dp)        :: res 
    
    real(dp)        :: x, rho_0, sig, step, err, nu
    real(dp)        :: dlogsigmadlogx, dum
    integer  :: ier
    
    rho_0 = rho_c * omega_m
    
    ! Lagrangian radius
    x     = (3.0_dp * m / (4.0_dp * pi * rho_0))**(1.0_dp / 3.0_dp)
    sig   = 10.0_dp**splevl(log10(m), logm, logsig, coeffspl2, dum, dum, ier)

    dlogsigmadlogx = splevl(log10(x), logx_arr, dlogsig_dlogx, coeffspl_deriv, dum, dum, ier)
    nu  = delta_c / (d(z) * sig)
    
    
    ! dn/dM = - rho_0/(3M^2) * f(nu) * d ln sigma / d ln x
    select case (trim(massfunc_name))
      case ('PS')
        res = -rho_0 / (3.0_dp * m**2) * probdist_PS(nu) * dlogsigmadlogx
      case ('ST')
        res = -rho_0 / (3.0_dp * m**2) * probdist_ST(nu) * dlogsigmadlogx
      case default
        write(*,*) 'Error: unknown mass function: ', trim(massfunc_name)
        res = 0.0_dp
    end select
  end function generic_dndM

  function massfrac(m, z) result(res)
    real(dp), intent(in) :: m, z
    real(dp) :: res, rho_0

    rho_0 = rho_c * omega_m
    res   = m * numdenm(m, z) / rho_0
  end function massfrac
  
  function probdist_ST(nu) result(res)
    real(dp), intent(in) :: nu
    double precision, parameter   :: small_a_st = 0.707_dp
    double precision, parameter   :: p_st = 0.300_dp
    double precision, parameter   :: cap_A_st = 0.3222_dp
    real(dp) :: nu2, res
    
    
    nu2 = nu * nu
    res = cap_A_st * sqrt(2.0_dp * small_a_st / pi)* (1+(small_a_st *nu2)**(-p_st)) * nu *  exp(-small_a_st *nu2 / 2.0_dp)
  end function probdist_ST


  function dfdnu_ST(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res

    real(dp), parameter :: cap_A_st = 0.3222_dp
    real(dp), parameter :: small_a_st = 0.707_dp
    real(dp), parameter :: p = 0.3_dp

    real(dp) :: C, anu2, power_term, exp_term

    C          = cap_A_st * sqrt(2.0_dp * small_a_st / pi)
    anu2       = small_a_st * nu**2
    power_term = anu2**(-p)
    exp_term   = exp(-anu2 / 2.0_dp)

    res = C * exp_term * ( &
          (1.0_dp + power_term) * (anu2-1) &
           + 2.0_dp * p * power_term ) !we absorbed the negative sign while calculating the mass integral

  end function dfdnu_ST
  
  function dfdnu_PS(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res

    res = sqrt(2.0_dp / pi) * (nu**2 -1.0_dp) * exp(-nu * nu / 2.0_dp)

  end function dfdnu_PS
  
  ! ---------------------------------------------------------------------------
  ! Power spectrum and variance
  ! ---------------------------------------------------------------------------
  function sigmasq(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res
    integer,  parameter :: lw = 2000, liw = lw / 4
    integer  :: ifail, inf
    real(dp) :: epsabs, epsrel, abserr
    real(dp), dimension(lw)  :: w
    integer,  dimension(liw) :: iw

    r      = x
    inf    = 1
    epsabs = 0.0_dp
    epsrel = 1.0e-3_dp
    ifail  = -1
    call d01amf(integrand_sigmasq, 0.0_dp, inf, epsabs, epsrel, res, abserr, w, lw, iw, liw, ifail)
    res = res / (2.0_dp * pi * pi)
  end function sigmasq

  function integrand_sigmasq(k) result(res)
    real(dp), intent(in) :: k
    real(dp) :: res

    if (r < 1.0e-7_dp) then
      res = (1.0_dp - k * k * r * r / 5.0_dp) * k * k * pspec(k)
    else
      res = (3.0_dp * (sin(k * r) - k * r * cos(k * r)))**2 * pspec(k) / (k * k * r**3)**2
    end if
  end function integrand_sigmasq

  function nu_parameter(mass, z) result(res)
    real(dp), intent(in) :: mass, z
    real(dp) :: res, rho_0, x

    rho_0 = rho_c * omega_m
    x     = (3.0_dp * mass / (4.0_dp * pi * rho_0))**(1.0_dp / 3.0_dp)
    res   = delta_c / (sigma(x) * d(z))
  end function nu_parameter

  function pspec(k) result(res)
    real(dp), intent(in) :: k
    real(dp) :: k_pivot, res, gam_h, n0
    
    k_pivot = 0.05
    gam_h = gamma * h
    if (abs(k) < 1.0e-40_dp) then
      res = 0.0_dp
    else
      n0  = ns + 0.5_dp * dn_dlnk * log(k / k_pivot)
      res = k**n0 * transfun(k / gam_h) * transfun(k / gam_h)
    end if
  end function pspec


  ! ---------------------------------------------------------------------------
  ! Dynamical time and virial overdensity
  ! ---------------------------------------------------------------------------
  function tdyn(zcoll) result(res)
    real(dp), intent(in) :: zcoll
    real(dp) :: res, grho

    grho = h * h * delvir(zcoll) * 6.6742e-8_dp * 1.8791e-29_dp * yrbysec**2 / &
           (hubbledist(zcoll) * hubbledist(zcoll))
    res  = sqrt(3.0_dp * pi / (32.0_dp * grho))
  end function tdyn

  function delvir(zcoll) result(res)
    real(dp), intent(in) :: zcoll
    real(dp) :: res, x

    x   = omega_z(zcoll) - 1.0_dp
    res = 18.0_dp * pi * pi + 82.0_dp * x - 39.0_dp * x * x
  end function delvir


  ! ---------------------------------------------------------------------------
  ! Cosmological distance and time measures
  ! ---------------------------------------------------------------------------
  function age(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res
    integer,  parameter :: lw = 2000, liw = lw / 4
    integer  :: inf, ifail
    real(dp) :: epsabs, epsrel, abserr
    real(dp), dimension(lw)  :: w
    integer,  dimension(liw) :: iw

    inf    = 1
    epsabs = 0.0_dp
    epsrel = 1.0e-6_dp
    ifail  = -1
    call d01amf(f_integrand, z, inf, epsabs, epsrel, res, abserr, w, lw, iw, liw, ifail)
  end function age

  function f_integrand(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res

    res = hubbledist(z) / (1.0_dp + z)
  end function f_integrand

  function hubbledist(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res, omega_k

    omega_k = 1.0_dp - (omega_l + omega_m + omega_r)
    res     = 1.0_dp / sqrt(omega_k * (1.0_dp + z)**2 + omega_m * (1.0_dp + z)**3 + &
                             omega_r * (1.0_dp + z)**4 + omega_l)
  end function hubbledist

  function omega_z(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res

    res = hubbledist(z) * hubbledist(z) * (1.0_dp + z)**3 * omega_m
  end function omega_z

  function comoving(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res, abserr
    integer  :: ifail, maxrul, nn, iparm
    real(dp), dimension(390) :: alpha

    maxrul = 9
    iparm  = 0
    ifail  = -1
    call d01arf(0.0_dp, z, d_h, 1.0e-6_dp, 0.0_dp, maxrul, iparm, abserr, res, nn, alpha, ifail)
  end function comoving

  function d_h(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res

    res = hubbledist(z)
  end function d_h

  function s_k(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res, k, rv, omega_k

    omega_k = 1.0_dp - omega_m - omega_l - omega_r
    k       = sqrt(abs(omega_k))
    rv      = comoving(z)
    if     (omega_k < 0.0_dp) then
      res = sin(rv * k) / k
    else if (omega_k == 0.0_dp) then
      res = rv
    else
      res = sinh(rv * k) / k
    end if
  end function s_k

  function angdist(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res

    res = s_k(z) / (1.0_dp + z)
  end function angdist

  function lumdist(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res

    res = angdist(z) * (1.0_dp + z) * (1.0_dp + z)
  end function lumdist
  
  function differential_comoving_volume(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res, chi, c_by_H0

    c_by_H0 = c_light / (100 * h * 1e05_dp)   ! in Mpc; Hubble parameter in cm/Mpc/sec km_to_cm = 1e05
    chi = c_by_H0 * comoving(z)            ! chi(z) in Mpc
    res = 4.0_dp * pi * c_by_H0 * chi**2 * hubbledist(z) ! Mpc^3/sr

  end function differential_comoving_volume


  ! ---------------------------------------------------------------------------
  ! Baryonic power spectrum
  ! ---------------------------------------------------------------------------
  
  function sigmasq_b(x, z, T, mu) result(res)
    real(dp), intent(in) :: x, z, T, mu
    real(dp) :: res

    integer,  parameter :: lw = 20000, liw = lw / 4
    integer  :: ifail, inf
    real(dp) :: abserr
    real(dp), dimension(lw)  :: w
    integer,  dimension(liw) :: iw

    cb_x  = x
    cb_z  = z
    cb_T  = T
    cb_mu = mu

    inf   = 1
    ifail = -1
    call d01amf(integrand_sigmasq_b, 0.0_dp, inf, 0.0_dp, 1.0e-8_dp, &
              res, abserr, w, lw, iw, liw, ifail)

    if (ifail /= 0) write(*, '(a, 4f12.4)') &
    'sigmasq_b: integration failed x, z, T, mu =', x, z, T, mu

    res = sigma8_norm * d(z)**2 * res / (2.0_dp * pi**2)
  end function sigmasq_b

  
  function integrand_sigmasq_b(k) result(res)
    real(dp), intent(in) :: k
    real(dp) :: res

    if (cb_x < 1.0e-7_dp) then
      res = (1.0_dp - k**2 * cb_x**2 / 5.0_dp) * k**2 * pspec_b(k)
    else
      res = (3.0_dp * (sin(k * cb_x) - k * cb_x * cos(k * cb_x)))**2 * &
          pspec_b(k) / (k**2 * cb_x**3)**2
    end if
  end function integrand_sigmasq_b

  function pspec_b(k) result(res)
    real(dp), intent(in) :: k
    real(dp) :: res

    res = pspec(k) / (1.0_dp + xbsq(cb_T, cb_z, cb_mu) * k**2)**2
  end function pspec_b

  function xbsq(T, z, mu) result(res) !CF05 eqn 40 and 14
    real(dp), intent(in) :: T, z, mu
    real(dp) :: vc_squared, res
    
    vc_squared= 2.0_dp * kboltz * T / &
        (mu * mprot)
    res = vc_squared/(3.0_dp * omega_m * h**2 * (1.0_dp + z)*1e14_dp) !the factor 1e14_dp comes from unit conversion, 
    !H0= 100*h*1e05 cm/sec/mpc, so square of it will give 1e14 
    !(mpc_cm)^2 factor from H0 conversion will cancel with the (mpc_cm)^2 from the v_c as finally x_b is in mpc unit
  end function xbsq
  
  subroutine setspline_sigma()
    integer  :: i, ier
    real(dp) :: x, rho_0, err

    rho_0 = rho_c * omega_m
    do i = 1, 100
      logm(i)   = 0.2_dp * i
      x         = (3.0_dp * 10.0_dp**logm(i) / (4.0_dp * pi * rho_0))**(1.0_dp / 3.0_dp)
      logx_arr(i) = log10(x)
      logsig(i) = log10(sigma(x))
      
      ! d log10(sigma) / d log10(x) via Ridders method
      dlogsig_dlogx(i) = dfridr(logsig_of_logx, log10(x), 0.01_dp, err)
    end do

    call spline(logsig, logm,   coeffspl,  ier)
    call spline(logm,   logsig, coeffspl2, ier)
    call spline(logx_arr, dlogsig_dlogx, coeffspl_deriv, ier)
  end subroutine setspline_sigma
  
  ! Helper — sigma as function of log10(x), needed by dfridr interface
  function logsig_of_logx(logx) result(res)
    real(dp), intent(in) :: logx
    real(dp) :: res
    res = log10(sigma(10.0_dp**logx))
  end function logsig_of_logx
  

end module backgroundCosmology_mod
