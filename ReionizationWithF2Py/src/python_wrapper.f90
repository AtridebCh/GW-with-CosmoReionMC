subroutine run_model(H0, ombh2, omch2, ns, sigma_8, &
                     esc_popII, esc_popIII, lambda0, &
                     zstart_in, zend_in, dz_in, Z_arraySize, &
                     Z_out, QH_Q_out, dNLLdz_out, gammaHI_out, &
                     sfr_pop2_out, sfr_pop3_out, dvc_out, &
                     D_L_out, age_out, tau_factor_out, ierr)
  use kinds_mod,        only: dp
  use constants_mod,    only: yrbysec
  use parameters_mod,   only: parameters_t, init_reion_defaults
  use variables_mod,    only: esc_II, esc_III, n, zstart, &
                              zend, dz, dNLLdz, QH, gammaHI, &
                              sfr_pop2_ion, sfr_pop2_neut, &
                              sfr_pop3_ion, sfr_pop3_neut, &
                              dvc_dz, D_L, age_Gyr, z, tau_factor
                              
                              
  use reionization_mod, only: initialize, filling, finalize
  implicit none


  ! inputs
  real(kind=8), intent(in)  :: H0, ombh2, omch2, ns, sigma_8
  real(kind=8), intent(in)  :: esc_popII, esc_popIII, lambda0
  real(kind=8), intent(in)  :: zstart_in, zend_in, dz_in
  integer, intent(in)       :: Z_arraySize
  
  ! outputs
  real(kind=8), intent(out) :: Z_out(0:Z_arraySize)
  real(kind=8), intent(out) :: QH_Q_out(0:Z_arraySize)
  real(kind=8), intent(out) :: dNLLdz_out(0:Z_arraySize)
  real(kind=8), intent(out) :: gammaHI_out(0:Z_arraySize)
  real(kind=8), intent(out) :: sfr_pop2_out(0:Z_arraySize)
  real(kind=8), intent(out) :: sfr_pop3_out(0:Z_arraySize)
  real(kind=8), intent(out) :: dvc_out(0:Z_arraySize)
  real(kind=8), intent(out) :: D_L_out(0:Z_arraySize)
  real(kind=8), intent(out) :: age_out(0:Z_arraySize)
  real(kind=8), intent(out) :: tau_factor_out(0:Z_arraySize)
  integer,      intent(out) :: ierr
  
  integer                   :: ik

  type(parameters_t) :: params

  ! set grid from Python
  zstart  = zstart_in
  zend    = zend_in
  dz      = dz_in

  ! set parameters
  params%cosmo%H0      = H0
  params%cosmo%ombh2   = ombh2
  params%cosmo%omch2   = omch2
  params%cosmo%ns      = ns
  params%cosmo%dn_dlnk = 0.0_dp
  params%cosmo%sigma_8 = sigma_8
  params%cosmo%m_wdm   = -1.0_dp

  call init_reion_defaults(params%reion)
  params%reion%esc_PopII       = esc_popII
  params%reion%esc_PopIII      = esc_popIII
  params%reion%lambda_0        = lambda0
  params%reion%Delta_H_overlap = 59.21_dp
  params%reion%e_sf_II         = 0.01_dp
  params%reion%e_sf_III        = 0.0_dp
  params%reion%e_QSO           = 0.36_dp
  params%reion%betaindex       = -2.2_dp
  params%reion%vc_min          = 13.3_dp

  ! ... fill params as before ...

  call initialize(params)
  !add dvc_dz in initialize
  esc_II(0:n)  = esc_popII
  esc_III(0:n) = esc_popIII
  call filling(params, ierr)
  
  do ik = 0, Z_arraySize
    Z_out(ik)        = z(ik)
    QH_Q_out(ik)     = QH(ik)%Q
    dNLLdz_out(ik)   = dNLLdz(ik)
    gammaHI_out(ik)  = gammaHI(ik)
    sfr_pop2_out(ik) = -(sfr_pop2_ion(ik)+sfr_pop2_neut(ik))*yrbysec
    sfr_pop3_out(ik) = -(sfr_pop3_ion(ik)+sfr_pop3_neut(ik))*yrbysec
    dvc_out(ik)      = dvc_dz(ik)
    D_L_out(ik)      = D_L(ik)
    age_out(ik)      = age_Gyr(ik)
    tau_factor_out(ik) = tau_factor(ik)
  end do
  
  call finalize()
  
end subroutine run_model

!add the mass function
