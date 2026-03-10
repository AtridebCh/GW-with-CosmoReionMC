subroutine get_outputs(nn, z_out, QH_Q_out, dNLLdz_out, gammaHIOut)
  use kinds_mod,        only: dp
  use variables_mod,    only: dNLLdz, QH, z, gammaHI
  !still need to add, rho_sfr
  implicit none

  integer,      intent(in)  :: nn
  real(kind=8), intent(out) :: z_out(0:nn)
  real(kind=8), intent(out) :: QH_Q_out(0:nn)
  real(kind=8), intent(out) :: dNLLdz_out(0:nn)
  real(kind=8), intent(out) :: gammaHI_out(0:nn)
  integer :: ik
  

  do ik = 0, nn
    z_out(ik)        = z(ik)
    QH_Q_out(ik)     = QH(ik)%Q
    dNLLdz_out(ik)   = dNLLdz(ik)
    gammaHI_out(ik)  = gammaHI(ik)
  end do

  tau_today = tau_elsc_today

end subroutine get_outputs
