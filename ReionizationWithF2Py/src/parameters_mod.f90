module parameters_mod
  use kinds_mod, only: dp
  use constants_mod, only: rho_crit_0, delta_c_z_0
  implicit none
  private
  
  public :: cosmo_params_t, reion_params_t, parameters_t,     &
            init_cosmo, init_params,     &
            h, omega_m, omega_l, omega_r, omega_k, gamma,     &
            omega_b, ombh2, ns, dn_dlnk, sigma_8, delta_c_z0, &
            rho_c, fzero, alpha_lo, alpha_hi, alpha_z,        &
            lambda_0, esc_PopII, esc_PopIII, e_sf_III, e_QSO, &
            Delta_H_overlap, betaindex, vc_min

  type :: cosmo_params_t
    ! Primary parameters
    real(dp) :: H0
    real(dp) :: ombh2
    real(dp) :: omch2
    !real(dp) :: As
    real(dp) :: ns
    real(dp) :: dn_dlnk
    real(dp) :: sigma_8
    real(dp) :: m_wdm    ! set to -1 if not used
    real(dp) :: delta_c_z0 = delta_c_z_0
    real(dp) :: rho_c_0 = rho_crit_0

    ! Derived quantities unlike Python properties,
    ! these must be set explicitly after construction
    real(dp) :: h
    real(dp) :: omb
    real(dp) :: omega_m
    real(dp) :: omega_l
    real(dp) :: omega_r
    real(dp) :: omega_k
    real(dp) :: gamma
    real(dp) :: rho_c
  end type cosmo_params_t


  type :: reion_params_t
    real(dp)          :: fzero
    real(dp)          :: alpha_lo
    real(dp)          :: alpha_hi
    real(dp)          :: alpha_z
    real(dp)          :: lambda_0
    real(dp)          :: esc_PopII
    real(dp)          :: esc_PopIII
    real(dp)          :: e_sf_III
    real(dp)          :: e_QSO
    real(dp)          :: Delta_H_overlap
    real(dp)          :: betaindex
    real(dp)          :: vc_min
  end type reion_params_t


  type :: parameters_t
    type(cosmo_params_t) :: cosmo
    type(reion_params_t) :: reion
  end type parameters_t
  
  real(dp), save, protected :: h, omega_m, omega_l, omega_r, omega_k
  real(dp), save, protected :: omega_b, ombh2, ns, sigma_8, dn_dlnk, delta_c_z0, rho_c, gamma
  real(dp), save, protected :: fzero, alpha_lo, alpha_hi, alpha_z
  real(dp), save, protected :: lambda_0, esc_PopII, esc_PopIII
  real(dp), save, protected :: e_sf_III, e_QSO
  real(dp), save, protected :: Delta_H_overlap, betaindex, vc_min

contains

  ! Equivalent of Python's @property derived quantities
  ! Call this after setting primary parameters
  subroutine init_cosmo(cosmo)
    type(cosmo_params_t), intent(inout) :: cosmo

    cosmo%h       = cosmo%H0 / 100.0_dp
    cosmo%omb     = cosmo%ombh2 / cosmo%h**2
    cosmo%omega_m = (cosmo%ombh2 + cosmo%omch2) / cosmo%h**2
    cosmo%omega_r = 0.0_dp
    cosmo%omega_k = 0.0_dp
    cosmo%omega_l = 1.0_dp - (cosmo%omega_m + cosmo%omega_r + cosmo%omega_k)
    cosmo%gamma   = cosmo%omega_m*cosmo%h
    cosmo%rho_c   = cosmo%rho_c_0*(cosmo%h)**2
  end subroutine init_cosmo

  
  subroutine init_params(params)
    type(parameters_t), intent(in) :: params

    h        = params%cosmo%h
    omega_m  = params%cosmo%omega_m
    omega_l  = params%cosmo%omega_l
    omega_r  = params%cosmo%omega_r
    omega_k  = params%cosmo%omega_k
    omega_b  = params%cosmo%omb
    ombh2    = params%cosmo%ombh2
    ns       = params%cosmo%ns
    sigma_8  = params%cosmo%sigma_8
    delta_c_z0  = params%cosmo%delta_c_z0
    gamma    = params%cosmo%gamma
    rho_c    = params%cosmo%rho_c_0*h**2

    fzero           = params%reion%fzero
    alpha_lo        = params%reion%alpha_lo
    alpha_hi        = params%reion%alpha_hi
    alpha_z         = params%reion%alpha_z
    lambda_0        = params%reion%lambda_0
    esc_PopII       = params%reion%esc_PopII
    esc_PopIII      = params%reion%esc_PopIII
    e_sf_III        = params%reion%e_sf_III
    e_QSO           = params%reion%e_QSO
    Delta_H_overlap = params%reion%Delta_H_overlap
    betaindex       = params%reion%betaindex
    vc_min          = params%reion%vc_min
  end subroutine init_params

end module parameters_mod
