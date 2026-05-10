module variables_mod
  use kinds_mod, only:dp
  use state_mod
  implicit none
  private

  ! Loop/control integers
  integer, public :: i, j, k, n 
  

  ! Logical flags
  logical, public ::  ifsummary

  character(len=1), public :: ifPopII
  
  logical, public :: ifprint  = .false.
  character(len=100), public:: outroot  = ''

  ! Redshift and density parameters
  real(dp), public :: zend, zstart, dz, mfilt, delta_c !delta_c =delta_c_z/norm()
  real(dp), public :: rho_b, rho_b_cgs, tau_elsc_today, lambda_0_min, lambda_0_max
  real(dp), public :: esc_II_param, esc_III_param
  real(dp), public :: Delta_V, betaindex, mumean, P_V_norm

  ! Input parameter
  character(len=100), public :: Pop3IMF='salpeter'
  character(len=100), public :: massfunc_name = 'PS'
  logical, public :: f_star_mass = .true.
  character(len=100), public :: f_coll_method = 'new'
  

  ! 1D allocatable arrays
  real(dp), allocatable, public :: z(:), tau_elsc(:), dtimedz(:), sigma(:), cpbycv(:), gammaHI(:)
  real(dp), allocatable, public :: t_array(:), t_H_array(:), dz_t_ff_array(:), &
                                   dvc_dz(:), D_L(:), age_Gyr(:)
  real(dp), allocatable, public :: esc_II(:), esc_III(:), f_starII(:), f_starIII(:), tau_factor(:)
  
  ! spline arrays for sigma
  real(dp), public ::  logm(100), logsig(100), logx_arr(100), dlogsig_dlogx(100)
  real(dp), public ::  coeffspl(3,100,100), coeffspl2(3,100,100), coeffspl_deriv(3, 100, 100)

  ! Derived type arrays
  type(fillingfactor_t), allocatable, public :: QH(:), QHe(:)

  type(ionsource_t), public :: dnphotdm, sigma_PI, sigma_PH, escfrac
  type(ionsource_t), allocatable, public :: dnphotdz_neut(:), dnphotdz_ion(:), &
                                            Gamma_PH(:), Gamma_PI(:)

  type(ionstate_t), allocatable, public :: neutral(:), HII(:), HeIII(:), global(:)
  type(ionstate_t), allocatable, public :: neutral_0(:), HII_0(:), HeIII_0(:), global_0(:)

  ! Star formation and photon rate arrays
  real(dp), allocatable, public :: dfcolldt_pop2_ion(:), dfcolldt_pop2_neut(:), &
                                   dfcolldt_pop3_ion(:), dfcolldt_pop3_neut(:)

  real(dp), allocatable, public :: mass_integral_pop3_ion(:), mass_integral_pop3_neut(:), &
                                   mass_integral_pop2_ion(:), mass_integral_pop2_neut(:), &
                                   lumfun_integral_qso(:)

  real(dp), allocatable, public :: sfr_pop3_ion(:), sfr_pop3_neut(:), &
                                   sfr_pop2_ion(:), sfr_pop2_neut(:)

  real(dp), allocatable, public :: dnphotdz_H(:), dnphotdz_He(:)
  real(dp), allocatable, public :: dNLLdz(:)

  real(dp), public :: zmin, zmax, zmean, fmin, fmax, fmean

end module variables_mod
