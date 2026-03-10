program driver
  use kinds_mod,        only: dp
  use parameters_mod,   only: parameters_t, cosmo_params_t, reion_params_t, &
                               init_cosmo, init_reion_defaults, init_params
  use variables_mod,    only: zstart, zend, dz, ifprint, outroot, &
                               esc_II, esc_III, n
  use reionization_mod, only: initialize, filling, finalize, validate_params

  implicit none

  type(parameters_t) :: params
  integer :: ierr

  ! set parameters
  params%cosmo%H0      = 67.74_dp
  params%cosmo%ombh2   = 0.0223_dp
  params%cosmo%omch2   = 0.11945_dp
  params%cosmo%ns      = 0.9667_dp
  params%cosmo%dn_dlnk = 0.0_dp
  !derived cosmological quantities
  params%cosmo%sigma_8 = 0.8159_dp
  params%cosmo%m_wdm   = -1.0_dp !not using any WDM

  call init_reion_defaults(params%reion)
  params%reion%esc_PopII       = 0.0036_dp
  params%reion%esc_PopIII      = 0.0_dp
  params%reion%lambda_0        = 2.44_dp
  params%reion%Delta_H_overlap = 59.21_dp
  params%reion%e_sf_II         = 0.01_dp
  params%reion%e_sf_III        = 0.0_dp
  params%reion%e_QSO           = 0.36_dp
  params%reion%betaindex       = -2.2_dp
  params%reion%vc_min          = 13.3_dp

  zstart   = 30.0_dp
  zend     = 0.0_dp
  dz       = 0.2_dp
  ifprint  = .true.
  outroot  = 'output/run1_'

  ! run
  call initialize(params)

  esc_II(0:n)  = params%reion%esc_PopII
  esc_III(0:n) = params%reion%esc_PopIII

  call filling(params, ierr)
  if (ierr /= 0) stop 'Error in filling()'

  if (ierr == 0) call finalize()

end program driver
