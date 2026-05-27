module stellar_mod
  use kinds_mod,        only: dp
  use constants_mod,    only: kboltz, mprot, yrbysec, hPlanck, Mpcbycm, &
                               nu_HII, nu_HeII, nu_HeIII, pi
  use state_mod,        only: ionstate_t, ionsource_t, lumsource_t, fillingfactor_t
  use parameters_mod,   only: parameters_t, h, omega_m, omega_l, omega_r, &
                               omega_k, omega_b, ombh2, ns, sigma_8, &
                               rho_c,  gamma, &
                               fzero, alpha_lo, alpha_hi, e_sf_III, e_QSO,        &
                               esc_PopII, esc_PopIII, lambda_0, Delta_H_overlap,   &
                               betaindex, vc_min
                           
  use variables_mod,    only: z, dz, dtimedz, k, delta_c, &
                               dfcolldt_pop2_ion, dfcolldt_pop2_neut, &
                               dfcolldt_pop3_ion, dfcolldt_pop3_neut, &
                               mass_integral_pop2_ion, mass_integral_pop2_neut, &
                               mass_integral_pop3_ion, mass_integral_pop3_neut, &
                               sfr_pop2_ion, sfr_pop2_neut, massfunc_name, &
                               sfr_pop3_ion, sfr_pop3_neut, &
                               Gamma_PI, Gamma_PH, rho_b, rho_b_cgs, &
                               dnphotdz_ion, dnphotdz_neut, dnphotdz_H, dnphotdz_He, &
                               QH, QHe, HII, HeIII, HII_0, HeIII_0, &
                               dnphotdm, sigma_PI, sigma_PH, escfrac, &
                               lumfun_integral_qso, f_coll_method, &
                               t_H_array, dz_t_ff_array, f_starII, f_starIII, &
                               logm, logsig, dlogsig_dlogx, coeffspl, &
                               coeffspl2, coeffspl_deriv
  use backgroundCosmology_mod, only: omega_z, hubbledist, delvir, &
                               xbsq, dfdnu_PS, d, sigma, nu_parameter, &
                               dfdnu_ST, generic_dndM, numdenm, growth_dynamical
                               
  use adaptint_mod,            only: d01amf, d01arf
  use onedspline_mod,          only: spline, splevl
  implicit none
  private
  public :: get_sfr, f_star, get_ionflux, lumfun_integral, sum_dnphotdz, QH

  ! QSO spectral index and normalisation
  real(dp), parameter :: alpha_qso   = 1.57_dp
  real(dp), parameter :: qso_lognorm = 18.05_dp
 

  ! callback state for QSO luminosity function
  real(dp), private :: cb_logphistar, cb_loglstar, cb_gamma1, cb_gamma2

  ! vc-mass conversion coefficient
  real(dp), parameter :: vc_mass_coeff = 3.504_dp / 1.0e8_dp**(1.0_dp / 3.0_dp)

  ! callback state for mass integrand routines
  real(dp) :: zmass
  !redshift for fpop3 chatterjee
  real(dp) :: z_pop3=16.0_dp

contains

  ! ---------------------------------------------------------------------------
  ! Star formation rates and photon production
  ! ---------------------------------------------------------------------------
  subroutine get_sfr()

    real(dp) :: prefactor, sfrfactor, fQ, f1mQ, qso_factor, dtz, rho_m

    prefactor  = omega_z(z(k))**0.55_dp
    fQ         = QH(k-1)%Q
    f1mQ       = 1.0_dp - fQ
    sfrfactor  = rho_b !in Msun/mpc^3 so that rho_sfr comes out in this unit
    rho_m = omega_m *rho_b/omega_b !remember, in MdndM_ST/MdndM_PS already has rho_m multiplied
    dtz = dtimedz(k)


    !new method suggested by Barun
    if (trim(f_coll_method) == 'new') then
      dfcolldt_pop2_ion(k) = max( &
         -1.0_dp/(t_H_array(k)*dz) * (1.0_dp + z(k)) * &
         ( mass_integral_pop2_new( &
             z(k), &
             vcmin_ion(z(k), HII(k)%T, HII(k)%X), &
             vcmax() ) &
         - mass_integral_pop2_new( &
             z(k-1), &
             vcmin_ion(z(k-1), HII(k-1)%T, HII(k-1)%X), &
             vcmax() ) ), &
         0.0_dp )


      dfcolldt_pop2_neut(k) = max( &
         -1.0_dp/(t_H_array(k)*dz) * (1.0_dp + z(k)) * &
         ( mass_integral_pop2_new( &
             z(k), &
             vcmin_neut(z(k)), &
             vcmax() ) &
         - mass_integral_pop2_new( &
             z(k-1), &
             vcmin_neut(z(k-1)), &
             vcmax() ) ), &
         0.0_dp )


      dfcolldt_pop3_ion(k) = max( &
         -1.0_dp/(t_H_array(k)*dz) * (1.0_dp + z(k)) * &
         ( mass_integral_pop3_new( &
             z(k), &
             vcmin_ion(z(k), HeIII(k)%T, HeIII(k)%X), &
             vcmax() ) &
         - mass_integral_pop3_new( &
             z(k-1), &
             vcmin_ion(z(k), HeIII(k-1)%T, HeIII(k-1)%X), &
             vcmax() ) ), &
         0.0_dp )


      dfcolldt_pop3_neut(k) = max( &
         -1.0_dp/(t_H_array(k)*dz) * (1.0_dp + z(k)) * &
         ( mass_integral_pop3_new( &
             z(k), &
             vcmin_neut(z(k)), &
             vcmax() ) &
         - mass_integral_pop3_new( &
             z(k-1), &
             vcmin_neut(z(k-1)), &
             vcmax() ) ), &
         0.0_dp )
    
    else   
      dfcolldt_pop2_ion(k)  = (prefactor / t_H_array(k)) * mass_integral_pop2(z(k), vcmin_ion(z(k), HII(k-1)%T,   HII(k-1)%X), vcmax())
      dfcolldt_pop2_neut(k) = (prefactor / t_H_array(k)) * mass_integral_pop2(z(k), vcmin_neut(z(k)), vcmax())
      dfcolldt_pop3_ion(k)  = (prefactor / t_H_array(k)) * mass_integral_pop3(z(k), vcmin_ion(z(k), HeIII(k-1)%T, HeIII(k-1)%X), vcmax())
      dfcolldt_pop3_neut(k) = (prefactor / t_H_array(k)) * mass_integral_pop3(z(k), vcmin_neut(z(k)), vcmax())   
    endif
    
    !open(unit=10, file='output_new_method.txt', status='unknown', position='append')
      !write(10,*) z(k), dfcolldt_pop2_ion(k), dfcolldt_pop2_neut(k)
    !close(10)
    
    mass_integral_pop2_neut(k) = dfcolldt_pop2_neut(k)
    mass_integral_pop2_ion(k)  = dfcolldt_pop2_ion(k)
    mass_integral_pop3_neut(k) = dfcolldt_pop3_neut(k)
    mass_integral_pop3_ion(k)  = dfcolldt_pop3_ion(k)

    ! stellar photon rates
    call set_stellar_dnphotdz(k, fQ, f1mQ, sfrfactor)

    ! QSO photon rates
    qso_factor = 10.0_dp**qso_lognorm * lumfun_integral_qso(k) / (alpha_qso * hPlanck)
    call set_qso_dnphotdz(k, fQ, f1mQ, qso_factor)

    ! convert dz -> dt
    call rescale_dnphotdz(k, dtz)

    ! photon sums
    call sum_dnphotdz(k)


    ! total photon rates
    dnphotdz_H(k) = sum_species_H(dnphotdz_neut(k)) + sum_species_H(dnphotdz_ion(k))
    dnphotdz_He(k) = sum_species_He(dnphotdz_neut(k)) + sum_species_He(dnphotdz_ion(k))

  end subroutine get_sfr


  ! ---------------------------------------------------------------------------
  ! Ionizing flux (mean free path and photoionization rates)
  ! ---------------------------------------------------------------------------
  subroutine get_ionflux()

    real(dp) :: lambda_0_H, lambda_0_He, zk, dtz, mfp_H, mfp_He
    real(dp) :: fac_H, fac_He

    zk  = z(k)
    dtz = dtimedz(k)

    lambda_0_H  = lambda_0 * &
                  sqrt(xbsq(HII_0(k-1)%T,  zk, 1.0_dp / HII_0(k-1)%X))  / (1.0_dp + zk)
    lambda_0_He = lambda_0 * &
                  sqrt(xbsq(HeIII_0(k-1)%T, zk, 1.0_dp / HeIII_0(k-1)%X)) / (1.0_dp + zk)

    if ((1.0_dp - QH(k-1)%F_V) >= 0.0_dp) then
      QH(k)%meanfreepath  = QH(k-1)%Q**(1.0_dp / 3.0_dp) * lambda_0_H  * Mpcbycm / &
                            (1.0_dp - QH(k-1)%F_V)**(2.0_dp / 3.0_dp)
    else
      QH(k)%meanfreepath  = QH(k-1)%meanfreepath
    end if

    if ((1.0_dp - QHe(k-1)%F_V) >= 0.0_dp) then
      QHe(k)%meanfreepath = QHe(k-1)%Q**(1.0_dp / 3.0_dp) * lambda_0_He * Mpcbycm / &
                            (1.0_dp - QHe(k-1)%F_V)**(2.0_dp / 3.0_dp)
    else
      QHe(k)%meanfreepath = QHe(k-1)%meanfreepath
    end if

    mfp_H  = QH(k)%meanfreepath
    mfp_He = QHe(k)%meanfreepath

    ! common factors: mfp * (1+z)^3 / (dtimedz * Mpcbycm^3)
    fac_H  = mfp_H  * (1.0_dp + zk)**3 / (dtz * Mpcbycm**3)
    fac_He = mfp_He * (1.0_dp + zk)**3 / (dtz * Mpcbycm**3)

    call set_gamma(k, fac_H, fac_He)
  end subroutine get_ionflux


  ! ---------------------------------------------------------------------------
  ! Minimum and maximum circular velocities
  ! ---------------------------------------------------------------------------
  
  function vc(mass,zcoll) result(res)
    !mass in h^{-1} M_sun, vc in km/s
    real(dp),        intent(in) :: mass, zcoll
    real(dp) :: res

    res=coeff_vc_mass(zcoll)*mass**(1.d0/3.d0)
  end function vc
  
  function vcmin_neut(z_in) result(res)
    real(dp),           intent(in) :: z_in
    real(dp) :: res, mmin

    res  = vc_min
  end function vcmin_neut

  function vcmin_ion(z_in, T, X) result(res)
    real(dp),           intent(in) :: z_in, T, X
    real(dp) :: res

    res = sqrt(2.0_dp * kboltz * T * X / mprot) * 1.0e-5_dp
    if (res < vcmin_neut(z_in)) res = vcmin_neut(z_in)
  end function vcmin_ion

  pure function vcmax() result(res)
    real(dp) :: res
    res = 1.0e9_dp
  end function vcmax


  ! ---------------------------------------------------------------------------
  ! Halo mass and circular velocity conversions
  ! ---------------------------------------------------------------------------
  function mass_from_vc(vc_in, zcoll) result(res)
    real(dp), intent(in) :: vc_in, zcoll
    real(dp) :: res

    res = (vc_in / coeff_vc_mass(zcoll))**3
  end function mass_from_vc

  function coeff_vc_mass(zcoll) result(res)
    real(dp), intent(in) :: zcoll
    real(dp) :: res, hsq

    hsq = 1.0_dp / hubbledist(zcoll)**2
    res = vc_mass_coeff * (0.5_dp * hsq * delvir(zcoll))**(1.0_dp / 6.0_dp) * &
          (1.0_dp - (2.0_dp * omega_l) / (3.0_dp * hsq * delvir(zcoll)))
  end function coeff_vc_mass


  ! ---------------------------------------------------------------------------
  ! Feedback suppression functions
  ! ---------------------------------------------------------------------------


  function fpop3(zcoll, m_log10) result(res)
    real(dp),           intent(in) :: zcoll, m_log10
    real(dp) :: res
    real(dp) :: mass_m, mass_min, mass_form, m_form, m_min
    real(dp) :: sigma_m, sigma_min, sigma_form, arg, dum
    integer  :: ier
    real(dp), parameter :: factor = 1.0e6_dp

    if (e_sf_III < 1.0e-6_dp) then
      res = 0.0_dp
      return
    end if

    mass_m   = 10.0_dp**m_log10
    mass_min = mass_from_vc(vcmin_neut(zcoll), zcoll) / h
    mass_form = mass_m - factor
    m_form    = log10(mass_form)
    m_min     = log10(mass_min)

    if (mass_m <= mass_min) then
      res = 0.0_dp
    else if (mass_form <= mass_min) then
      res = 1.0_dp
    else
      sigma_m    = 10.0_dp**splevl(m_log10, logm, logsig, coeffspl2, dum, dum, ier)
      sigma_min  = 10.0_dp**splevl(m_min,   logm, logsig, coeffspl2, dum, dum, ier)
      sigma_form = 10.0_dp**splevl(m_form,  logm, logsig, coeffspl2, dum, dum, ier)

      arg = sqrt((sigma_form - sigma_m) / (sigma_min - sigma_form))
      if (arg < 1.0e-5_dp) then
        res = arg
      else
        res = 2.0_dp * atan(arg) / pi
      end if
    end if
  end function fpop3
 


  ! ---------------------------------------------------------------------------
  ! Mass function integrals shared implementation for PS (original CF05)
  ! ---------------------------------------------------------------------------
  function mass_integral_pop2(zcoll, vc_min, vc_max) result(res)
    real(dp), intent(in) :: zcoll, vc_min, vc_max
    real(dp) :: res

    res = mass_integral_generic(zcoll, vc_min, vc_max, mass_integrand_pop2)
  end function mass_integral_pop2

  function mass_integral_pop3(zcoll, vc_min, vc_max) result(res)
    real(dp), intent(in) :: zcoll, vc_min, vc_max
    real(dp) :: res

    res = mass_integral_generic(zcoll, vc_min, vc_max, mass_integrand_pop3)
  end function mass_integral_pop3

  function mass_integral_generic(zcoll, vc_min, vc_max, integrand) result(res)
    real(dp), intent(in) :: zcoll, vc_min, vc_max
    real(dp) :: res
    integer,  parameter :: lw = 2000, liw = lw / 4
    integer  :: inf, ifail, maxrul, nn, iparm
    real(dp) :: epsabs, epsrel, abserr, mmin, numin, mmax, numax, boundary_term
    real(dp), dimension(lw)  :: w
    integer,  dimension(liw) :: iw
    real(dp), dimension(390) :: alpha
    interface
      function integrand(nu) result(val)
        import dp
        real(dp), intent(in) :: nu
        real(dp) :: val
      end function integrand
    end interface

    if (abs(vc_min - vc_max) < 1.0e-8_dp .or. vc_min > vc_max) then
      res = 0.0_dp
      return
    end if

    zmass  = zcoll
    mmin   = mass_from_vc(vc_min, zcoll) / h
    numin  = nu_parameter(mmin, zcoll)

    epsabs = 0.0_dp
    epsrel = 1.0e-3_dp
    ifail  = -1
    !if f_star(M)
    !boundary_term = f_star(mmin) * probdist_ST(numin) * numin

    if (vc_max > 1.0e5_dp) then
      inf = 1
      call d01amf(integrand, numin, inf, epsabs, epsrel, res, abserr, w, lw, iw, liw, ifail)      
    else
      mmax   = mass_from_vc(vc_max, zcoll) / h
      numax  = nu_parameter(mmax, zcoll)
      maxrul = 9
      iparm  = 0
      call d01arf(numin, numax, integrand, epsrel, epsabs, maxrul, iparm, abserr, res, nn, alpha, ifail)
      !if f_star(M)
      !res = boundary_term + res
    end if

    if (res < 0.0_dp) res = 0.0_dp
  end function mass_integral_generic

  function mass_integrand_pop2(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res, sig, m, dum
    integer  :: ier

    !sig = delta_c / (d(zmass) * nu)
    ! for dynamical dark energy
    sig = delta_c / (growth_dynamical(zmass) * nu)
    m   = splevl(log10(sig), logsig, logm, coeffspl, dum, dum, ier)
    
    select case (trim(massfunc_name))
      case ('PS')
        res = dfdnu_PS(nu) * (1.0_dp - fpop3(zmass, m)) *f_star(10**m, fzero, alpha_lo, alpha_hi)
      case ('ST')
        res = dfdnu_ST(nu) * (1.0_dp - fpop3(zmass, m)) *f_star(10**m, fzero, alpha_lo, alpha_hi)
    end select
  end function mass_integrand_pop2

  function mass_integrand_pop3(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res, sig, m, dum
    integer  :: ier

    !sig = delta_c / (d(zmass) * nu)
    ! for dynamical dark energy
    sig = delta_c / (growth_dynamical(zmass) * nu)
    m   = splevl(log10(sig), logsig, logm, coeffspl, dum, dum, ier)
    
    select case (trim(massfunc_name))
      case ('PS')
        res = dfdnu_PS(nu) * fpop3(zmass, m)*f_star(10**m, fzero, alpha_lo, alpha_hi)
      case ('ST')
        res = dfdnu_ST(nu) * fpop3(zmass, m) *f_star(10**m, fzero, alpha_lo, alpha_hi)
    end select
  end function mass_integrand_pop3


! ---------------------------------------------------------------------------
  ! Mass function integrals shared implementation for ST (added by Atri) and f_star(M)
  ! ---------------------------------------------------------------------------
  
  function fnu_PS(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res

    res = sqrt(2.0_dp/pi) * &
        exp(-0.5_dp * nu * nu)

   end function fnu_PS
  
  function f_star(Mh, fzero, alpha_lo, alpha_hi ) result(res)
    real(dp), intent(in) :: Mh, fzero, alpha_lo, alpha_hi   ! log10 of halo mass
    real(dp) :: res
    real(dp), parameter :: Mp = 10.0**11.7_dp          ! log10(Mp) = 11.7
    real(dp) :: ratio

    ratio = Mh/Mp

    res = 2*fzero / (ratio**(alpha_lo) + ratio**(-alpha_hi))

  end function f_star
  
  function fpop3_chatterjee(zcoll, dz_t_ff) result(res)
    real(dp), intent(in) :: zcoll, dz_t_ff
    real(dp) :: res
    
    res = 0.5*(1+tanh(abs((zcoll-z_pop3)/dz_t_ff)))
  end function fpop3_chatterjee
  
  function mass_integral_pop2_new(zcoll, vc_min, vc_max) result(res)
    real(dp), intent(in) :: zcoll, vc_min, vc_max
    real(dp) :: res

    res = mass_integral_generic(zcoll, vc_min, vc_max, mass_integrand_pop2_new)
  end function mass_integral_pop2_new

  function mass_integral_pop3_new(zcoll, vc_min, vc_max) result(res)
    real(dp), intent(in) :: zcoll, vc_min, vc_max
    real(dp) :: res

    res = mass_integral_generic(zcoll, vc_min, vc_max, mass_integrand_pop3_new)
  end function mass_integral_pop3_new 

  function mass_integrand_pop2_new(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res, sig, m, dum
    integer  :: ier

    !sig = delta_c / (d(zmass) * nu)
    !for dynamical dark energy
    sig = delta_c / (growth_dynamical(zmass) * nu)
    m   = splevl(log10(sig), logsig, logm, coeffspl, dum, dum, ier)  
     
    res = (1.0_dp- fpop3(zmass, m)) * fnu_PS(nu) *f_star(10**m, fzero, alpha_lo, alpha_hi)
  end function mass_integrand_pop2_new

  function mass_integrand_pop3_new(nu) result(res)
    real(dp), intent(in) :: nu
    real(dp) :: res, sig, m, dum
    integer  :: ier

    !sig = delta_c / (d(zmass) * nu)
    ! for dynamical dark energy
    sig = delta_c / (growth_dynamical(zmass) * nu)
    m   = splevl(log10(sig), logsig, logm, coeffspl, dum, dum, ier) 
    res = fpop3(zmass, m) * fnu_PS(nu)*f_star(10**m, fzero, alpha_lo, alpha_hi)
  end function mass_integrand_pop3_new

  ! ---------------------------------------------------------------------------
  ! Private helpers to reduce repetition in get_sfr and get_ionflux
  ! ---------------------------------------------------------------------------
  subroutine set_stellar_dnphotdz(ik, fQ, f1mQ, sfrfactor)
    integer,            intent(in) :: ik
    real(dp),           intent(in) :: fQ, f1mQ, sfrfactor

    dnphotdz_ion(ik)%pop2%HII   = -fQ   * dnphotdm%pop2%HII   * mass_integral_pop2_ion(ik)  * sfrfactor * escfrac%pop2%HII   * f_starII(ik)
    dnphotdz_ion(ik)%pop2%HeII  = -fQ   * dnphotdm%pop2%HeII  * mass_integral_pop2_ion(ik)  * sfrfactor * escfrac%pop2%HeII  * f_starII(ik)
    dnphotdz_ion(ik)%pop2%HeIII = -fQ   * dnphotdm%pop2%HeIII * mass_integral_pop2_ion(ik)  * sfrfactor * escfrac%pop2%HeIII * f_starII(ik)
    sfr_pop2_ion(ik)            = -fQ   * mass_integral_pop2_ion(ik) * sfrfactor * f_starII(ik)

    dnphotdz_neut(ik)%pop2%HII   = -f1mQ * dnphotdm%pop2%HII   * mass_integral_pop2_neut(ik) * sfrfactor * escfrac%pop2%HII   * f_starII(ik)
    dnphotdz_neut(ik)%pop2%HeII  = -f1mQ * dnphotdm%pop2%HeII  * mass_integral_pop2_neut(ik) * sfrfactor * escfrac%pop2%HeII  * f_starII(ik)
    dnphotdz_neut(ik)%pop2%HeIII = -f1mQ * dnphotdm%pop2%HeIII * mass_integral_pop2_neut(ik) * sfrfactor * escfrac%pop2%HeIII * f_starII(ik)
    sfr_pop2_neut(ik)            = -f1mQ * mass_integral_pop2_neut(ik) * sfrfactor * f_starII(ik)

    dnphotdz_ion(ik)%pop3%HII   = -fQ   * dnphotdm%pop3%HII   * mass_integral_pop3_ion(ik)  * sfrfactor * escfrac%pop3%HII   * f_starIII(ik)
    dnphotdz_ion(ik)%pop3%HeII  = -fQ   * dnphotdm%pop3%HeII  * mass_integral_pop3_ion(ik)  * sfrfactor * escfrac%pop3%HeII  * f_starIII(ik)
    dnphotdz_ion(ik)%pop3%HeIII = -fQ   * dnphotdm%pop3%HeIII * mass_integral_pop3_ion(ik)  * sfrfactor * escfrac%pop3%HeIII * f_starIII(ik)
    sfr_pop3_ion(ik)            = -fQ   * mass_integral_pop3_ion(ik)  * sfrfactor * f_starIII(ik)

    dnphotdz_neut(ik)%pop3%HII   = -f1mQ * dnphotdm%pop3%HII   * mass_integral_pop3_neut(ik) * sfrfactor * escfrac%pop3%HII   * f_starIII(ik)
    dnphotdz_neut(ik)%pop3%HeII  = -f1mQ * dnphotdm%pop3%HeII  * mass_integral_pop3_neut(ik) * sfrfactor * escfrac%pop3%HeII  * f_starIII(ik)
    dnphotdz_neut(ik)%pop3%HeIII = -f1mQ * dnphotdm%pop3%HeIII * mass_integral_pop3_neut(ik) * sfrfactor * escfrac%pop3%HeIII * f_starIII(ik)
    sfr_pop3_neut(ik)            = -f1mQ * mass_integral_pop3_neut(ik) * sfrfactor * f_starIII(ik)
  end subroutine set_stellar_dnphotdz

  subroutine set_qso_dnphotdz(ik, fQ, f1mQ, qso_factor)
    integer,  intent(in) :: ik
    real(dp), intent(in) :: fQ, f1mQ, qso_factor

    real(dp) :: frac_H, frac_He2, frac_He3

    frac_H   = 1.0_dp - (nu_HII / nu_HeII)**alpha_qso
    frac_He2 = (nu_HII / nu_HeII)**alpha_qso - (nu_HII / nu_HeIII)**alpha_qso
    frac_He3 = (nu_HII / nu_HeIII)**alpha_qso

    dnphotdz_neut(ik)%QSO%HII   = -fQ   * qso_factor * frac_H
    dnphotdz_neut(ik)%QSO%HeII  = -fQ   * qso_factor * frac_He2
    dnphotdz_neut(ik)%QSO%HeIII = -fQ   * qso_factor * frac_He3
    dnphotdz_ion(ik)%QSO%HII    = -f1mQ * qso_factor * frac_H
    dnphotdz_ion(ik)%QSO%HeII   = -f1mQ * qso_factor * frac_He2
    dnphotdz_ion(ik)%QSO%HeIII  = -f1mQ * qso_factor * frac_He3
    
  end subroutine set_qso_dnphotdz

  subroutine rescale_dnphotdz(ik, dtimedz_k)
    integer, intent(in) :: ik
    real(dp), intent(in):: dtimedz_k
    real(dp) :: scale

    scale = -dtimedz_k !sec

    dnphotdz_ion(ik)%pop2%HII   = dnphotdz_ion(ik)%pop2%HII   * scale
    dnphotdz_ion(ik)%pop2%HeII  = dnphotdz_ion(ik)%pop2%HeII  * scale
    dnphotdz_ion(ik)%pop2%HeIII = dnphotdz_ion(ik)%pop2%HeIII * scale
    dnphotdz_neut(ik)%pop2%HII  = dnphotdz_neut(ik)%pop2%HII  * scale
    dnphotdz_neut(ik)%pop2%HeII = dnphotdz_neut(ik)%pop2%HeII * scale
    dnphotdz_neut(ik)%pop2%HeIII= dnphotdz_neut(ik)%pop2%HeIII* scale

    dnphotdz_ion(ik)%pop3%HII   = dnphotdz_ion(ik)%pop3%HII   * scale
    dnphotdz_ion(ik)%pop3%HeII  = dnphotdz_ion(ik)%pop3%HeII  * scale
    dnphotdz_ion(ik)%pop3%HeIII = dnphotdz_ion(ik)%pop3%HeIII * scale
    dnphotdz_neut(ik)%pop3%HII  = dnphotdz_neut(ik)%pop3%HII  * scale
    dnphotdz_neut(ik)%pop3%HeII = dnphotdz_neut(ik)%pop3%HeII * scale
    dnphotdz_neut(ik)%pop3%HeIII= dnphotdz_neut(ik)%pop3%HeIII* scale

    dnphotdz_ion(ik)%QSO%HII    = dnphotdz_ion(ik)%QSO%HII    * scale
    dnphotdz_ion(ik)%QSO%HeII   = dnphotdz_ion(ik)%QSO%HeII   * scale
    dnphotdz_ion(ik)%QSO%HeIII  = dnphotdz_ion(ik)%QSO%HeIII  * scale
    dnphotdz_neut(ik)%QSO%HII   = dnphotdz_neut(ik)%QSO%HII   * scale
    dnphotdz_neut(ik)%QSO%HeII  = dnphotdz_neut(ik)%QSO%HeII  * scale
    dnphotdz_neut(ik)%QSO%HeIII = dnphotdz_neut(ik)%QSO%HeIII * scale
  end subroutine rescale_dnphotdz

  pure function sum_species_H(src) result(res)
    type(ionsource_t), intent(in) :: src
    real(dp) :: res

    res = src%pop2%HII  + src%pop3%HII  + src%QSO%HII  + &
          src%pop2%HeII + src%pop3%HeII + src%QSO%HeII
  end function sum_species_H

  pure function sum_species_He(src) result(res)
    type(ionsource_t), intent(in) :: src
    real(dp) :: res

    res = src%pop2%HeIII + src%pop3%HeIII + src%QSO%HeIII
  end function sum_species_He
  
  subroutine sum_dnphotdz(ik)
    integer, intent(in) :: ik
    dnphotdz_H(ik)  = sum_species_H(dnphotdz_neut(ik))  + sum_species_H(dnphotdz_ion(ik))
    dnphotdz_He(ik) = sum_species_He(dnphotdz_neut(ik)) + sum_species_He(dnphotdz_ion(ik))
  end subroutine sum_dnphotdz

  subroutine set_gamma(ik, fac_H, fac_He)
    integer,  intent(in) :: ik
    real(dp), intent(in) :: fac_H, fac_He

    ! photoionization rates
    Gamma_PI(ik)%pop2%HII   = fac_H  * sigma_PI%pop2%HII   * (dnphotdz_neut(ik)%pop2%HII   + dnphotdz_ion(ik)%pop2%HII)
    Gamma_PI(ik)%pop3%HII   = fac_H  * sigma_PI%pop3%HII   * (dnphotdz_neut(ik)%pop3%HII   + dnphotdz_ion(ik)%pop3%HII)
    Gamma_PI(ik)%QSO%HII    = fac_H  * sigma_PI%QSO%HII    * (dnphotdz_neut(ik)%QSO%HII    + dnphotdz_ion(ik)%QSO%HII)
    Gamma_PI(ik)%pop2%HeII  = fac_H  * sigma_PI%pop2%HeII  * (dnphotdz_neut(ik)%pop2%HeII  + dnphotdz_ion(ik)%pop2%HeII)
    Gamma_PI(ik)%pop3%HeII  = fac_H  * sigma_PI%pop3%HeII  * (dnphotdz_neut(ik)%pop3%HeII  + dnphotdz_ion(ik)%pop3%HeII)
    Gamma_PI(ik)%QSO%HeII   = fac_H  * sigma_PI%QSO%HeII   * (dnphotdz_neut(ik)%QSO%HeII   + dnphotdz_ion(ik)%QSO%HeII)
    Gamma_PI(ik)%pop2%HeIII = fac_He * sigma_PI%pop2%HeIII * (dnphotdz_neut(ik)%pop2%HeIII + dnphotdz_ion(ik)%pop2%HeIII)
    Gamma_PI(ik)%pop3%HeIII = fac_He * sigma_PI%pop3%HeIII * (dnphotdz_neut(ik)%pop3%HeIII + dnphotdz_ion(ik)%pop3%HeIII)
    Gamma_PI(ik)%QSO%HeIII  = fac_He * sigma_PI%QSO%HeIII  * (dnphotdz_neut(ik)%QSO%HeIII  + dnphotdz_ion(ik)%QSO%HeIII)

    ! photoheating rates
    Gamma_PH(ik)%pop2%HII   = hPlanck * nu_HII   * fac_H  * sigma_PH%pop2%HII   * (dnphotdz_neut(ik)%pop2%HII   + dnphotdz_ion(ik)%pop2%HII)
    Gamma_PH(ik)%pop3%HII   = hPlanck * nu_HII   * fac_H  * sigma_PH%pop3%HII   * (dnphotdz_neut(ik)%pop3%HII   + dnphotdz_ion(ik)%pop3%HII)
    Gamma_PH(ik)%QSO%HII    = hPlanck * nu_HII   * fac_H  * sigma_PH%QSO%HII    * (dnphotdz_neut(ik)%QSO%HII    + dnphotdz_ion(ik)%QSO%HII)
    Gamma_PH(ik)%pop2%HeII  = hPlanck * nu_HeII  * fac_H  * sigma_PH%pop2%HeII  * (dnphotdz_neut(ik)%pop2%HeII  + dnphotdz_ion(ik)%pop2%HeII)
    Gamma_PH(ik)%pop3%HeII  = hPlanck * nu_HeII  * fac_H  * sigma_PH%pop3%HeII  * (dnphotdz_neut(ik)%pop3%HeII  + dnphotdz_ion(ik)%pop3%HeII)
    Gamma_PH(ik)%QSO%HeII   = hPlanck * nu_HeII  * fac_H  * sigma_PH%QSO%HeII   * (dnphotdz_neut(ik)%QSO%HeII   + dnphotdz_ion(ik)%QSO%HeII)
    Gamma_PH(ik)%pop2%HeIII = hPlanck * nu_HeIII * fac_He * sigma_PH%pop2%HeIII * (dnphotdz_neut(ik)%pop2%HeIII + dnphotdz_ion(ik)%pop2%HeIII)
    Gamma_PH(ik)%pop3%HeIII = hPlanck * nu_HeIII * fac_He * sigma_PH%pop3%HeIII * (dnphotdz_neut(ik)%pop3%HeIII + dnphotdz_ion(ik)%pop3%HeIII)
    Gamma_PH(ik)%QSO%HeIII  = hPlanck * nu_HeIII * fac_He * sigma_PH%QSO%HeIII  * (dnphotdz_neut(ik)%QSO%HeIII  + dnphotdz_ion(ik)%QSO%HeIII)
  end subroutine set_gamma
  
  function lumfun_integral(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res

    integer,  parameter :: lw = 2000, liw = lw / 4
    real(dp), parameter :: P0 = -4.8250643_dp,  P1 = 13.035753_dp,  &
                           P2 =  0.63150872_dp, P3 = -11.763560_dp, &
                           P4 = -14.249833_dp,  P5 =  0.41698725_dp, &
                           P6 = -0.62298947_dp, P7 =  2.1744386_dp,  &
                           P8 =  1.4599393_dp,  P9 = -0.79280099_dp
    integer  :: ifail, inf
    real(dp) :: abserr, xi
    real(dp), dimension(lw)  :: w
    integer,  dimension(liw) :: iw

    if (z > 12.0_dp) then
      res = 0.0_dp
      return
    end if

    xi = log10((1.0_dp + z) / 3.0_dp)

    cb_logphistar = P0
    cb_loglstar   = P1 + P2*xi + P3*xi**2 + P4*xi**3
    cb_gamma1     = P5 * 10.0_dp**(P6 * xi)
    cb_gamma2     = max(1.3_dp, 2.0_dp * P7 / (10.0_dp**(P8*xi) + 10.0_dp**(P9*xi)))

    inf    = 1
    ifail  = -1
    call d01amf(lumfun_integrand, 0.0_dp, inf, 0.0_dp, 1.0e-4_dp, &
                res, abserr, w, lw, iw, liw, ifail)
  end function lumfun_integral

  function lumfun_integrand(lb) result(res)
    real(dp), intent(in) :: lb
    real(dp) :: res

    real(dp), parameter :: Q0 =  8.99833_dp,  Q1 =  6.24800_dp, &
                           Q2 = -0.370587_dp, Q3 = -0.0115970_dp
    real(dp) :: psi, quasar_lumfun, lstarby10, lband

    psi           = 10.0_dp**cb_logphistar / (lb**cb_gamma1 + lb**cb_gamma2)
    quasar_lumfun = 10.0_dp**cb_loglstar * psi / (log(10.0_dp) * lb)
    lstarby10     = 10.0_dp**cb_loglstar * 1.0e-10_dp
    lband         = Q0 * (lb * lstarby10)**Q3 + Q1 * (lb * lstarby10)**Q2
    res           = quasar_lumfun * lb / lband
  end function lumfun_integrand
  

end module stellar_mod
