module reionization_mod
  use kinds_mod,        only: dp
  use constants_mod,    only: mprot, yrbysec, pi, Mpcbycm, Y_He, mp_by_mHe, Msun, c_light
  use state_mod,        only: ionstate_t, ionsource_t, lumsource_t, species_t
  use parameters_mod,only: parameters_t, init_cosmo, h, omega_m, omega_l, omega_r, &
                               omega_k,omega_b, ombh2, ns, sigma_8, &
                               delta_c_z0, rho_c, gamma, init_params, &
                               esc_PopII, esc_PopIII, lambda_0, Delta_H_overlap,   &
                               e_sf_II, e_sf_III, e_QSO, betaindex, vc_min
  use variables_mod
  use thermal_mod,      only: compute_species, update_ionstate, get_recrates, get_ionrates
  use stellar_mod,      only: get_sfr, get_ionflux, lumfun_integral, setspline_sigma
  use recrates_mod,     only: R_HII_e_A, R_HII_e_B, R_HeIII_e_A, R_HeIII_e_B
  use backgroundCosmology_mod, only: hubbledist, xbsq, sigmasq_b, f_integrand, &
                                    set_sigma8_norm, sigma8_norm, tdyn,   &
                                    differential_comoving_volume, lumdist, age
  use inhomoReion_mod
  use utils_mod,        only: summarize, write_summary
  use SEDreader_mod,   only: get_sed
  implicit none
  private
  public :: initialize, finalize, filling, validate_params, set_escfrac

contains

  ! ---------------------------------------------------------------------------
  ! Parameter validation
  ! ---------------------------------------------------------------------------
  function validate_params() result(ok)
    logical :: ok

    ok = .true.
    if (h < 0.2_dp .or. h > 1.0_dp) then
      ok = .false.
      write(*, *) '  Warning: h has units of 100 km/s/Mpc. You have:', h
    end if

    if (Y_He < 0.2_dp .or. Y_He > 0.8_dp) then
      ok = .false.
      write(*, *) '  Warning: Y_He is the helium fraction of baryons. You have:', Y_He
    end if

    if (omega_b < 0.001_dp .or. omega_m < 0.0_dp .or. &
        omega_b > 1.0_dp   .or. omega_m > 3.0_dp) then
      ok = .false.
      write(*, *) '  Warning: matter densities are out of expected range'
    end if
  end function validate_params


  ! ---------------------------------------------------------------------------
  ! Initialization
  ! ---------------------------------------------------------------------------
  subroutine initialize(params)
    type(parameters_t), intent(inout) :: params

    integer :: ik
    call init_cosmo(params%cosmo)
    call init_params(params)            ! copy to protected module vars
    
    call set_sigma8_norm()          ! normalise to sigma_8, fixing sigma8_norm
    
    delta_c   = delta_c_z0 / sqrt(sigma8_norm)
    rho_b     = rho_c * omega_b !Msun/Mpc^3
    rho_b_cgs = rho_b*Msun/Mpcbycm**3

    if (.not. validate_params()) stop 'Stopped due to parameter error'

    call setspline_sigma()
    call get_SED()

    dz = -abs(dz)
    n  = nint(abs((zend - zstart) / dz))

    call allocate_arrays(n)
    call zero_arrays()

    ! redshift grid and time arrays
    do ik = 0, n
      z(ik)          = ik * dz + zstart
      dtimedz(ik)    = - f_integrand(z(ik)) / (100.0_dp*h*1e05_dp/Mpcbycm)     !sec
      tau_factor(ik) = dtimedz(ik)*(1+z(ik))**3
      t_H_array(ik)  = hubbledist(z(ik)) / (100.0_dp*h*1e05_dp/Mpcbycm)        !sec
      dvc_dz(ik)     = differential_comoving_volume(z(ik))
      D_L(ik)        = lumdist(z(ik))*c_light / (100 * h * 1e05_dp) !mpc
      age_Gyr(ik)    = age(z(ik))/(100.0_dp*h*1e05_dp/Mpcbycm)/(1e09_dp * yrbysec)
      
      dz_t_ff_array(ik)  = tdyn(z(ik))*yrbysec/dtimedz(ik) !tdyn=t_ff is in year so first multiplying wiht yrbysec to convert to sec
      lumfun_integral_qso(ik) = e_QSO * lumfun_integral(z(ik))
    end do

    ! initial ionization states
    call set_initial_ionstate()

    ! initial filling factor state
    sigma(0) = sqrt(sigmasq_b(0.0_dp, z(0), global(0)%T, 1.0_dp))
    call get_LN_norm(sigma(0))
    

    QH(0)%Q            = 1.0e-8_dp
    QH(0)%Delta        = Delta_H_overlap
    QH(0)%F_V          = F_V(sigma(0), QH(0)%Delta)
    QH(0)%F_M          = F_M(sigma(0), QH(0)%Delta)
    QH(0)%R            = R(sigma(0),   QH(0)%Delta)
    QH(0)%meanfreepath = 1.0e-16_dp

    QHe(0)%Q            = 1.0e-10_dp
    QHe(0)%Delta        = Delta_H_overlap
    QHe(0)%F_V          = F_V(sigma(0), QHe(0)%Delta)
    QHe(0)%F_M          = F_M(sigma(0), QHe(0)%Delta)
    QHe(0)%R            = R(sigma(0),   QHe(0)%Delta)
    QHe(0)%meanfreepath = 0.0_dp
    
    call compute_all_species(0)
    
  end subroutine initialize


  ! ---------------------------------------------------------------------------
  ! Main time loop
  ! ---------------------------------------------------------------------------
  subroutine filling(params, ierr)
    type(parameters_t), intent(in)  :: params
    integer,            intent(out) :: ierr

    type(lumsource_t) :: totGamma_PI, totGamma_PH
    real(dp) :: n_e_H, n_e_He, n_H, n_He, n_e, n_HI
    real(dp) :: conv_factor
    integer  :: ik

    ierr = 0

    conv_factor = rho_b_cgs / mprot !basically n_b in cgs

    do ik = 1, n
      k = ik

      call set_escfrac()
      call get_sfr()

      sigma(ik) = sqrt(sigmasq_b(0.0_dp, z(ik), global(ik-1)%T, 1.0_dp / global(ik-1)%X))
      
      
      call get_LN_norm(sigma(ik))
      call get_ionflux()
      

      ! copy previous state forward for all regions
      call copy_ionstate_frac_T(neutral(ik),   neutral(ik-1))
      call copy_ionstate_frac_T(HII(ik),       HII(ik-1))
      call copy_ionstate_frac_T(HeIII(ik),     HeIII(ik-1))

      call get_recrates(z(ik), neutral(ik), dtimedz(ik), 1.0_dp, 1.0_dp)
      !write(*,*)'before recrate for HII', HII(ik)%recrate, HII(ik)%coolrate
      call get_recrates(z(ik), HII(ik), dtimedz(ik), QH(ik-1)%R,  1.0_dp)
      !write(*,*)'after recrate for HII', HII(ik)%recrate, HII(ik)%coolrate
      call get_recrates(z(ik), HeIII(ik), dtimedz(ik), QHe(ik-1)%R, 1.0_dp)
      
      
      ! --- ionized/neutral regions ---
      neutral(ik)%ionrate  = species_t(0.0_dp, 0.0_dp, 0.0_dp)

      call build_totgamma(ik, totGamma_PI, totGamma_PH, include_HeIII=.false.)
      !write(*,*)'totGamma_PI%HII', totGamma_PI
      !write(*,*)'before ionrate for HII', HII(ik)%ionrate, HII(ik)%heatrate
      call get_ionrates(z(ik), HII(ik), dtimedz(ik), totGamma_PI, totGamma_PH, QH(ik-1)%Q)
      !write(*,*)'after ionrate for HII', HII(ik)%ionrate, HII(ik)%heatrate
      
      totGamma_PI%HeIII  = sum_HeIII(Gamma_PI(ik))  * (QHe(ik-1)%meanfreepath / QH(ik-1)%meanfreepath) * (QH(ik-1)%Q / QHe(ik-1)%Q)
      totGamma_PH%HeIII  = sum_HeIII(Gamma_PH(ik))  * (QHe(ik-1)%meanfreepath / QH(ik-1)%meanfreepath) * (QH(ik-1)%Q / QHe(ik-1)%Q)
      call get_ionrates(z(ik), HeIII(ik), dtimedz(ik), totGamma_PI, totGamma_PH, QH(ik-1)%Q)
      
      neutral(ik)%heatrate = species_t(0.0_dp, 0.0_dp, 0.0_dp)

      call update_ionstate(z(ik), dz, dtimedz(ik), neutral(ik))
      !write(*,*)'HII region before ionstate', HII(ik)%frac%HI, HII(ik)%frac%HeI, HII(ik)%frac%HeIII, HII(ik)%X_HII, HII(ik)%X_HeII, HII(ik)%X_e, HII(ik)%X
      call update_ionstate(z(ik), dz, dtimedz(ik), HII(ik))
      !write(*,*)'HII region after ionstate', HII(ik)%frac%HI, HII(ik)%frac%HeI, HII(ik)%frac%HeIII, HII(ik)%X_HII, HII(ik)%X_HeII, HII(ik)%X_e, HII(ik)%X
      call update_ionstate(z(ik), dz, dtimedz(ik), HeIII(ik))
      
      call copy_ionstate_frac_T(neutral_0(ik), neutral_0(ik-1))
      call copy_ionstate_frac_T(HII_0(ik),     HII_0(ik-1))
      call copy_ionstate_frac_T(HeIII_0(ik),   HeIII_0(ik-1))

      ! --- Delta = 1 regions (uniform density) ---
      call get_recrates(z(ik), neutral_0(ik), dtimedz(ik), 1.0_dp, 1.0_dp)
      call get_recrates(z(ik), HII_0(ik), dtimedz(ik), 1.0_dp, 1.0_dp)
      call get_recrates(z(ik), HeIII_0(ik), dtimedz(ik), 1.0_dp, 1.0_dp)
      

      call copy_ionheat_rates(neutral_0(ik), neutral(ik))
      call copy_ionheat_rates(HII_0(ik),     HII(ik))
      call copy_ionheat_rates(HeIII_0(ik),   HeIII(ik))

      call update_ionstate(z(ik), dz, dtimedz(ik), neutral_0(ik))
      call update_ionstate(z(ik), dz, dtimedz(ik), HII_0(ik))
      call update_ionstate(z(ik), dz, dtimedz(ik), HeIII_0(ik))
      

      ! --- Delta = 0.1 regions (underdense) ---
      call copy_ionstate_frac_T(neutral_m1(ik),neutral_m1(ik-1))
      call copy_ionstate_frac_T(HII_m1(ik),    HII_m1(ik-1))
      call copy_ionstate_frac_T(HeIII_m1(ik),  HeIII_m1(ik-1))
      
      call get_recrates(z(ik), neutral_m1(ik), dtimedz(ik), 1.0_dp, 0.1_dp)
      call get_recrates(z(ik), HII_m1(ik), dtimedz(ik), 1.0_dp, 0.1_dp)
      call get_recrates(z(ik), HeIII_m1(ik), dtimedz(ik), 1.0_dp, 0.1_dp)

      call copy_ionheat_rates(neutral_m1(ik), neutral_0(ik))
      call copy_ionheat_rates(HII_m1(ik),     HII_0(ik))
      call copy_ionheat_rates(HeIII_m1(ik),   HeIII_0(ik))

      call update_ionstate(z(ik), dz, dtimedz(ik), neutral_m1(ik))
      call update_ionstate(z(ik), dz, dtimedz(ik), HII_m1(ik))
      call update_ionstate(z(ik), dz, dtimedz(ik), HeIII_m1(ik))


      ! --- Delta = 10 regions (overdense) ---
      call copy_ionstate_frac_T(neutral_p1(ik),neutral_p1(ik-1))
      call copy_ionstate_frac_T(HII_p1(ik),    HII_p1(ik-1))
      call copy_ionstate_frac_T(HeIII_p1(ik),  HeIII_p1(ik-1))  
      
      call get_recrates(z(ik), neutral_p1(ik), dtimedz(ik), 1.0_dp, 10.0_dp)
      call get_recrates(z(ik), HII_p1(ik), dtimedz(ik), 1.0_dp, 10.0_dp)
      call get_recrates(z(ik), HeIII_p1(ik), dtimedz(ik), 1.0_dp, 10.0_dp)

      call copy_ionheat_rates(neutral_p1(ik), neutral_0(ik))
      call copy_ionheat_rates(HII_p1(ik),     HII_0(ik))
      call copy_ionheat_rates(HeIII_p1(ik),   HeIII_0(ik))

      call update_ionstate(z(ik), dz, dtimedz(ik), neutral_p1(ik))
      call update_ionstate(z(ik), dz, dtimedz(ik), HII_p1(ik))
      call update_ionstate(z(ik), dz, dtimedz(ik), HeIII_p1(ik))

      ! --- ionization evolution ---
      n_e_H = (HII(ik)%X_e + (QHe(ik-1)%Q / QH(ik-1)%Q) * HeIII(ik)%X_e) * conv_factor
      n_H   = (HII(ik)%X_HII + HII(ik)%frac%HI) * conv_factor
      call update_Q(z(ik), dz, dtimedz(ik), QH(ik), QH(ik-1), sigma(ik), sigma(ik-1), &
                    dnphotdz_H(ik), n_e_H, n_H, R_HII_e_A(HII(ik)%T) + R_HII_e_B(HII(ik)%T), ierr)
      if (ierr /= 0) return
      
      n_e_He = HeIII(ik)%X_e * conv_factor
      n_He   = (HeIII(ik)%frac%HeIII + HeIII(ik)%X_HeII + HeIII(ik)%frac%HeI) * conv_factor
      call update_Q(z(ik), dz, dtimedz(ik), QHe(ik), QHe(ik-1), sigma(ik), sigma(ik-1), &
                    dnphotdz_He(ik), n_e_He, n_He, R_HeIII_e_A(HeIII(ik)%T) + R_HeIII_e_B(HeIII(ik)%T), ierr)
      if (ierr /= 0) return
      

      ! --- global state averages ---
      call mix_global(ik, global(ik),    neutral(ik),    HII(ik),    HeIII(ik),    use_F_M=.true.)
      call compute_species(global(ik))
      
      call mix_global(ik, global_0(ik),  neutral_0(ik),  HII_0(ik),  HeIII_0(ik),  use_F_M=.false.)
      call compute_species(global_0(ik))
      
      call mix_global(ik, global_m1(ik), neutral_m1(ik), HII_m1(ik), HeIII_m1(ik), use_F_M=.false.)
      call compute_species(global_m1(ik))
      
      call mix_global(ik, global_p1(ik), neutral_p1(ik), HII_p1(ik), HeIII_p1(ik), use_F_M=.false.)
      call compute_species(global_p1(ik))


      cpbycv(ik) = 0.5_dp * log10(global_p1(ik)%T / global_m1(ik)%T) + 1.0_dp

      ! --- mean free path and optical depth ---
      !dNLLdz(ik) = 1.0_dp / ((1.0_dp + z(ik)) * &
      !             sqrt(pi) * QH(ik)%meanfreepath / (Mpcbycm * hubbledist(z(ik)) * 3.0e3_dp / params%cosmo%h))
      dNLLdz(ik)=sqrt(pi)*QH(ik)%meanfreepath/(Mpcbycm*hubbledist(z(ik))*3.d3/h)
      dNLLdz(ik)=1.0_dp/((1.0_dp+z(ik))*dNLLdz(k))

      n_HI = global_0(ik)%frac%HI * conv_factor 
      n_e  = global_0(ik)%X_e     * conv_factor
      gammaHI(ik) = totGamma_PI%HII
      tau_elsc(ik) = tau_elsc(ik-1) + dz * dtimedz(ik) * n_e * 3.0e10_dp * 6.652e-25_dp * (1.0_dp + z(ik))**3
      !write(*,*) 'QH', z(k), QH(ik)%Q
    end do

    tau_elsc_today    = tau_elsc(n)
    tau_elsc(0:n)     = tau_elsc_today - tau_elsc(0:n)
    !if (ifprint) write(*, *) "tau_elsc_today =", tau_elsc_today
    tau_factor = tau_factor * conv_factor  * c_light *6.652e-25_dp !sigma_T = 6.652e-25_dp

    call summarize(params)
  end subroutine filling


  ! ---------------------------------------------------------------------------
  ! Finalization
  ! ---------------------------------------------------------------------------
  subroutine finalize()
    call deallocate_arrays()
  end subroutine finalize


  ! ---------------------------------------------------------------------------
  ! Escape fraction
  ! ---------------------------------------------------------------------------
  subroutine set_escfrac()

    character(len=1) :: flag

    flag = adjustl(ifPopII)

    esc_II_param  = esc_II(k)
    esc_III_param = esc_III(k)

    ! simple case: ifPopII is not T/t/F/f
    if (flag /= 'T' .and. flag /= 't' .and. flag /= 'F' .and. flag /= 'f') then
      escfrac%pop2%HII   = esc_II_param
      escfrac%pop2%HeII  = esc_II_param
      escfrac%pop2%HeIII = esc_II_param
      escfrac%pop3%HII   = esc_III_param
      escfrac%pop3%HeII  = esc_III_param
      escfrac%pop3%HeIII = esc_III_param
      !write(*,*) "z, esc_II, esc_III",z(k),escfrac%pop2%HII,escfrac%pop3%HII
      return
    end if
  end subroutine set_escfrac

  subroutine copy_ionstate_frac_T(dst, src)
    type(ionstate_t), intent(inout) :: dst
    type(ionstate_t), intent(in)    :: src

    dst%frac = src%frac
    dst%T    = src%T
  end subroutine copy_ionstate_frac_T

  subroutine copy_ionheat_rates(dst, src)
    type(ionstate_t), intent(inout) :: dst
    type(ionstate_t), intent(in)    :: src

    dst%ionrate  = src%ionrate
    dst%heatrate  = src%heatrate
  end subroutine copy_ionheat_rates

  subroutine build_totgamma(ik, tot_PI, tot_PH, include_HeIII)
    integer,           intent(in)  :: ik
    type(lumsource_t), intent(out) :: tot_PI, tot_PH
    logical,           intent(in)  :: include_HeIII

    tot_PI%HII  = Gamma_PI(ik)%pop2%HII  + Gamma_PI(ik)%pop3%HII  + Gamma_PI(ik)%QSO%HII
    tot_PI%HeII = Gamma_PI(ik)%pop2%HeII + Gamma_PI(ik)%pop3%HeII + Gamma_PI(ik)%QSO%HeII
    tot_PH%HII  = Gamma_PH(ik)%pop2%HII  + Gamma_PH(ik)%pop3%HII  + Gamma_PH(ik)%QSO%HII
    tot_PH%HeII = Gamma_PH(ik)%pop2%HeII + Gamma_PH(ik)%pop3%HeII + Gamma_PH(ik)%QSO%HeII

    if (include_HeIII) then
      tot_PI%HeIII = sum_HeIII(Gamma_PI(ik))
      tot_PH%HeIII = sum_HeIII(Gamma_PH(ik))
    else
      tot_PI%HeIII = 0.0_dp
      tot_PH%HeIII = 0.0_dp
    end if
  end subroutine build_totgamma

  pure function sum_HeIII(src) result(res)
    type(ionsource_t), intent(in) :: src
    real(dp) :: res

    res = src%pop2%HeIII + src%pop3%HeIII + src%QSO%HeIII
  end function sum_HeIII

  subroutine mix_global(ik, glob, reg_neut, reg_HII, reg_HeIII, use_F_M)
    integer,          intent(in)    :: ik
    type(ionstate_t), intent(inout) :: glob
    type(ionstate_t), intent(in)    :: reg_neut, reg_HII, reg_HeIII
    logical,          intent(in)    :: use_F_M

    real(dp) :: w_neut, w_HII, w_HeIII

    if (use_F_M) then
      w_neut  = 1.0_dp - QH(ik)%Q  * QH(ik)%F_M
      w_HII   = QH(ik)%Q  * QH(ik)%F_M - QHe(ik)%Q * QHe(ik)%F_M
      w_HeIII = QHe(ik)%Q * QHe(ik)%F_M
    else
      w_neut  = 1.0_dp - QH(ik)%Q
      w_HII   = QH(ik)%Q - QHe(ik)%Q
      w_HeIII = QHe(ik)%Q
    end if

    glob%T          = w_neut * reg_neut%T          + w_HII * reg_HII%T          + w_HeIII * reg_HeIII%T
    glob%frac%HI    = w_neut * reg_neut%frac%HI    + w_HII * reg_HII%frac%HI    + w_HeIII * reg_HeIII%frac%HI
    glob%frac%HeI   = w_neut * reg_neut%frac%HeI   + w_HII * reg_HII%frac%HeI   + w_HeIII * reg_HeIII%frac%HeI
    glob%frac%HeIII = w_neut * reg_neut%frac%HeIII + w_HII * reg_HII%frac%HeIII + w_HeIII * reg_HeIII%frac%HeIII
  end subroutine mix_global

  subroutine set_initial_ionstate()

    real(dp) :: HI_ini, HeI_ini, T_ini
    type(ionstate_t) :: ini

    HI_ini  = 1.0_dp - Y_He - 1.2e-5_dp * sqrt(omega_m) / &
              (h * omega_b)
    HeI_ini = mp_by_mHe * Y_He
    T_ini   = 0.0099_dp * 2.726_dp * (1.0_dp + z(0))**2

    ini = ionstate_t(species_t(HI_ini, HeI_ini, 0.0_dp), &
                     species_t(0.0_dp, 0.0_dp,  0.0_dp), &
                     species_t(0.0_dp, 0.0_dp,  0.0_dp), &
                     species_t(0.0_dp, 0.0_dp,  0.0_dp), &
                     species_t(0.0_dp, 0.0_dp,  0.0_dp), &
                     T_ini, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)

    
    neutral(0)   = ini
    
    HII(0)    = ini;  HeIII(0)   = ini;  global(0)   = ini
    neutral_0(0) = ini;  HII_0(0)  = ini;  HeIII_0(0) = ini;  global_0(0) = ini
    neutral_m1(0)= ini;  HII_m1(0) = ini;  HeIII_m1(0)= ini;  global_m1(0)= ini
    neutral_p1(0)= ini;  HII_p1(0) = ini;  HeIII_p1(0)= ini;  global_p1(0)= ini
  end subroutine set_initial_ionstate

  subroutine compute_all_species(ik)
    integer,            intent(in) :: ik

    call compute_species(neutral(ik))
    call compute_species(HII(ik))
    call compute_species(HeIII(ik))
    call compute_species(global(ik))
    call compute_species(neutral_0(ik))
    call compute_species(HII_0(ik))
    call compute_species(HeIII_0(ik))
    call compute_species(global_0(ik))
    call compute_species(neutral_m1(ik))
    call compute_species(HII_m1(ik))
    call compute_species(HeIII_m1(ik))
    call compute_species(global_m1(ik))
    call compute_species(neutral_p1(ik))
    call compute_species(HII_p1(ik))
    call compute_species(HeIII_p1(ik))
    call compute_species(global_p1(ik))
  end subroutine compute_all_species

  subroutine allocate_arrays(nn)
    integer, intent(in) :: nn

    allocate(z(0:nn), QH(0:nn), QHe(0:nn), &
             dnphotdz_neut(0:nn), dnphotdz_ion(0:nn), &
             tau_elsc(0:nn), cpbycv(0:nn), gammaHI(0:nn), &
             dtimedz(0:nn), t_H_array(0:nn), dz_t_ff_array(0:nn), &
             dvc_dz(0:nn), D_L(0:nn), neutral(0:nn), HII(0:nn), &
             HeIII(0:nn), global(0:nn), age_Gyr(0:nn), &
             neutral_0(0:nn), HII_0(0:nn), HeIII_0(0:nn), global_0(0:nn), &
             neutral_m1(0:nn), HII_m1(0:nn), HeIII_m1(0:nn), global_m1(0:nn), &
             neutral_p1(0:nn), HII_p1(0:nn), HeIII_p1(0:nn), global_p1(0:nn), &
             Gamma_PH(0:nn), Gamma_PI(0:nn), &
             dnphotdz_H(0:nn), dnphotdz_He(0:nn), sigma(0:nn), &
             dNLLdz(0:nn), lumfun_integral_qso(0:nn), &
             mass_integral_pop2_ion(0:nn), mass_integral_pop2_neut(0:nn), &
             mass_integral_pop3_ion(0:nn), mass_integral_pop3_neut(0:nn), &
             dfcolldt_pop2_ion(0:nn), dfcolldt_pop2_neut(0:nn), &
             dfcolldt_pop3_ion(0:nn), dfcolldt_pop3_neut(0:nn), &
             sfr_pop2_ion(0:nn), sfr_pop2_neut(0:nn), &
             sfr_pop3_ion(0:nn), sfr_pop3_neut(0:nn), &
             esc_II(0:nn), esc_III(0:nn), tau_factor(0:nn))
  end subroutine allocate_arrays

  subroutine zero_arrays()
    sfr_pop2_ion(:)     = 0.0_dp;  sfr_pop2_neut(:)     = 0.0_dp
    sfr_pop3_ion(:)     = 0.0_dp;  sfr_pop3_neut(:)     = 0.0_dp
    tau_elsc(:)         = 0.0_dp
    dfcolldt_pop2_ion(0) = 0.0_dp;  dfcolldt_pop3_ion(0) = 0.0_dp
    gammaHI(:)  = 0.0_dp
  end subroutine zero_arrays

  subroutine deallocate_arrays()
    deallocate(z, QH, QHe, dnphotdz_neut, dnphotdz_ion, &
               tau_elsc, cpbycv, gammaHI, dtimedz, t_H_array, &
               age_Gyr, dvc_dz, D_L, dz_t_ff_array, neutral, HII, HeIII, global, &
               neutral_0, HII_0, HeIII_0, global_0, &
               neutral_m1, HII_m1, HeIII_m1, global_m1, &
               neutral_p1, HII_p1, HeIII_p1, global_p1, &
               Gamma_PH, Gamma_PI, dnphotdz_H, dnphotdz_He, &
               sigma, dNLLdz, lumfun_integral_qso, &
               mass_integral_pop2_ion, mass_integral_pop2_neut, &
               mass_integral_pop3_ion, mass_integral_pop3_neut, &
               dfcolldt_pop2_ion, dfcolldt_pop2_neut, &
               dfcolldt_pop3_ion, dfcolldt_pop3_neut, &
               sfr_pop2_ion, sfr_pop2_neut, sfr_pop3_ion, sfr_pop3_neut, &
               esc_II, esc_III, tau_factor)

  end subroutine deallocate_arrays
  

end module reionization_mod
