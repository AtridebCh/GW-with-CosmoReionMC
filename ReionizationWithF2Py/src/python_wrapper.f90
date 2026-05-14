subroutine run_model(H0, ombh2, omch2, ns, sigma_8, omega_zero, omega_a, &
                     fzero, alpha_lo, alpha_hi, alpha_z, esc_PopII, lambda0, &
                     zstart_in, zend_in, dz_in, Z_arraySize, &
                     Z_out, QH_Q_out, dNLLdz_out, gammaHI_out, &
                     sfr_pop2_out, sfr_pop3_out, dvc_out, &
                     D_L_out, age_out, tau_factor_out, &
                     omega_dyn_out, omega_de_out, ierr)
  
  use kinds_mod,        only: dp
  use constants_mod,    only: yrbysec
  use parameters_mod,   only: parameters_t
  use variables_mod,    only: n, zstart, zend, dz, &
                              dNLLdz, QH, gammaHI, &
                              sfr_pop2_ion, sfr_pop2_neut, &
                              sfr_pop3_ion, sfr_pop3_neut, &
                              dvc_dz, D_L, age_Gyr, z, &
                              tau_factor, omega_dyn, omega_de
                              
                              
  use reionization_mod, only: initialize, filling, finalize
  implicit none


  ! inputs
  real(kind=8), intent(in)  :: H0, ombh2, omch2, ns, sigma_8, omega_zero, omega_a
  real(kind=8), intent(in)  :: fzero, alpha_lo, alpha_hi, alpha_z
  real(kind=8), intent(in)  :: esc_PopII, lambda0
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
  real(kind=8), intent(out) :: omega_dyn_out(0:Z_arraySize)
  real(kind=8), intent(out) :: omega_de_out(0:Z_arraySize)
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
  params%cosmo%omega_zero = omega_zero
  params%cosmo%omega_a = omega_a
  params%cosmo%m_wdm   = -1.0_dp

  params%reion%fzero           = fzero
  params%reion%alpha_lo        = alpha_lo
  params%reion%alpha_hi        = alpha_hi
  params%reion%alpha_z         = alpha_z
  params%reion%esc_PopII       = esc_PopII
  params%reion%lambda_0        = lambda0
  params%reion%esc_PopIII      = 0.0_dp
  params%reion%e_sf_III        = 0.0_dp
  params%reion%e_QSO           = 0.36_dp
  params%reion%Delta_H_overlap = 59.21_dp
  params%reion%betaindex       = -2.2_dp
  params%reion%vc_min          = 13.3_dp

  ! ... fill params as before ...

  call initialize(params) !e_sf and e_esc are done inside initialize
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
    omega_dyn_out(ik)  = omega_dyn(ik)
    omega_de_out(ik)   = omega_de(ik)
  end do
  
  call finalize()
  
end subroutine run_model



!add the mass function

subroutine run_dndm(m_arr, z_arr, nm, nz, dndm_out, ierr)
  use kinds_mod,            only: dp
  use backgroundCosmology_mod,  only:  generic_dndM  !numdenm
  implicit none

  ! inputs
  integer,      intent(in)  :: nm, nz
  !f2py intent(hide) :: nm, nz
  real(kind=8), intent(in)  :: m_arr(nm)
  real(kind=8), intent(in)  :: z_arr(nz)

  ! outputs
  real(kind=8), intent(out) :: dndm_out(nm, nz)
  integer,      intent(out) :: ierr

  integer :: i, j

  ierr = 0
  do j = 1, nz
    do i = 1, nm
      dndm_out(i, j) =  generic_dndM(m_arr(i), z_arr(j)) !numdenm(m_arr(i), z_arr(j))
    end do
  end do

end subroutine run_dndm

subroutine get_sfe(m_arr, nm, fzero, alpha_lo, alpha_hi, sfe_out, ierr)
  use kinds_mod,            only: dp
  use stellar_mod,          only: f_star
  implicit none
  
  ! inputs
  real(kind=8), intent(in)  :: fzero, alpha_lo, alpha_hi
  integer,      intent(in)  :: nm
  !f2py intent(hide) :: nm
  real(kind=8), intent(in)  :: m_arr(nm)
  
  ! outputs
  real(kind=8), intent(out) :: sfe_out(nm)
  integer,      intent(out) :: ierr
  
  integer :: i
  
  ierr = 0
  do i = 1, nm
    sfe_out(i) =  f_star(m_arr(i), fzero, alpha_lo, alpha_hi)
  end do
end subroutine get_sfe

