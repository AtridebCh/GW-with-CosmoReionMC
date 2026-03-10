module inhomoReion_mod
  use kinds_mod,     only: dp
  use constants_mod, only: pi, two_pi, Mpcbycm
  use state_mod, only:fillingfactor_t
  use parameters_mod,only: h, omega_m, omega_l, omega_r, &
                           omega_k,omega_b, ombh2, ns, sigma_8, &
                           rho_c, gamma, esc_PopII, esc_PopIII, &
                           lambda_0, Delta_H_overlap,   &
                           e_sf_II, e_sf_III, e_QSO, betaindex, vc_min
  use variables_mod,  only: Delta_V, mumean, P_V_norm
  use backgroundCosmology_mod, only: f_integrand, hubbledist
  implicit none
  private
  public :: update_Q, get_LN_norm, P_V, F_V, F_M, R, Delta_i_F_M, Delta_i_F_V

contains

  ! ---------------------------------------------------------------------------
  ! Log-normal PDF normalisation
  ! ---------------------------------------------------------------------------
  subroutine get_LN_norm(sigma)
    real(dp), intent(in) :: sigma

    real(dp) :: b1, b2

    b1 = betaindex + 1.0_dp
    b2 = betaindex + 2.0_dp

    Delta_V = (0.5_dp * (1.0_dp - erf(sigma * b1 / sqrt(2.0_dp))) - &
               exp(-0.5_dp * sigma**2 * b1**2) / (sqrt(two_pi) * b1 * sigma)) / &
              (0.5_dp * exp(sigma**2 * (betaindex + 1.5_dp)) * &
               (1.0_dp - erf(sigma * b2 / sqrt(2.0_dp))) - &
               exp(-0.5_dp * sigma**2 * b1**2) / (sqrt(two_pi) * b2 * sigma))

    mumean  = log(Delta_V) + sigma**2 * b1

    P_V_norm = 0.5_dp * (1.0_dp + erf((log(Delta_V) - mumean) / (sigma * sqrt(2.0_dp)))) - &
               exp(-0.5_dp * ((log(Delta_V) - mumean) / sigma)**2) / (sqrt(two_pi) * b1 * sigma)
    P_V_norm = 1.0_dp / P_V_norm
  end subroutine get_LN_norm


  ! ---------------------------------------------------------------------------
  ! Log-normal + power-law PDF
  ! ---------------------------------------------------------------------------
  pure function P_V(sigma, Delta) result(res)
    real(dp), intent(in) :: sigma, Delta
    real(dp) :: res

    if (Delta <= Delta_V) then
      res = exp(-0.5_dp * ((log(Delta) - mumean) / sigma)**2) / &
            (sqrt(two_pi) * Delta * sigma)
    else
      res = (Delta / Delta_V)**betaindex * &
            exp(-0.5_dp * ((log(Delta_V) - mumean) / sigma)**2) / &
            (sqrt(two_pi) * Delta_V * sigma)
    end if
    res = res * P_V_norm
  end function P_V


  ! ---------------------------------------------------------------------------
  ! Volume filling factor
  ! ---------------------------------------------------------------------------
  pure function F_V(sigma, Delta_i) result(res)
    real(dp), intent(in) :: sigma, Delta_i
    real(dp) :: res

    if (Delta_i <= Delta_V) then
      res = 0.5_dp * (1.0_dp + erf((log(Delta_i) - mumean) / (sigma * sqrt(2.0_dp))))
    else
      res = 0.5_dp * (1.0_dp + erf((log(Delta_V) - mumean) / (sigma * sqrt(2.0_dp)))) + &
            exp(-0.5_dp * ((log(Delta_V) - mumean) / sigma)**2) / &
            (sqrt(two_pi) * Delta_V**(betaindex + 1.0_dp) * sigma) * &
            (Delta_i**(betaindex + 1.0_dp) - Delta_V**(betaindex + 1.0_dp)) / (betaindex + 1.0_dp)
    end if
    res = P_V_norm * res
  end function F_V


  ! ---------------------------------------------------------------------------
  ! Mass filling factor
  ! ---------------------------------------------------------------------------
  pure function F_M(sigma, Delta_i) result(res)
    real(dp), intent(in) :: sigma, Delta_i
    real(dp) :: res

    if (Delta_i <= Delta_V) then
      res = 0.5_dp * exp(mumean + 0.5_dp * sigma**2) * &
            (1.0_dp + erf((log(Delta_i) - mumean - sigma**2) / (sigma * sqrt(2.0_dp))))
    else
      res = 0.5_dp * exp(mumean + 0.5_dp * sigma**2) * &
            (1.0_dp + erf((log(Delta_V) - mumean - sigma**2) / (sigma * sqrt(2.0_dp)))) + &
            exp(-0.5_dp * ((log(Delta_V) - mumean) / sigma)**2) / &
            (sqrt(two_pi) * Delta_V**(betaindex + 1.0_dp) * sigma) * &
            (Delta_i**(betaindex + 2.0_dp) - Delta_V**(betaindex + 2.0_dp)) / (betaindex + 2.0_dp)
    end if
    res = P_V_norm * res
  end function F_M


  ! ---------------------------------------------------------------------------
  ! Clumping factor integral
  ! ---------------------------------------------------------------------------
  pure function R(sigma, Delta_i) result(res)
    real(dp), intent(in) :: sigma, Delta_i
    real(dp) :: res

    if (Delta_i <= Delta_V) then
      res = 0.5_dp * exp(2.0_dp * (mumean + sigma**2)) * &
            (1.0_dp + erf((log(Delta_i) - mumean - 2.0_dp * sigma**2) / (sigma * sqrt(2.0_dp))))
    else
      res = 0.5_dp * exp(2.0_dp * (mumean + sigma**2)) * &
            (1.0_dp + erf((log(Delta_V) - mumean - 2.0_dp * sigma**2) / (sigma * sqrt(2.0_dp)))) + &
            exp(-0.5_dp * ((log(Delta_V) - mumean) / sigma)**2) / &
            (sqrt(two_pi) * Delta_V**(betaindex + 1.0_dp) * sigma) * &
            (Delta_i**(betaindex + 3.0_dp) - Delta_V**(betaindex + 3.0_dp)) / (betaindex + 3.0_dp)
    end if
    res = P_V_norm * res
  end function R


  ! ---------------------------------------------------------------------------
  ! Invert F_M to find Delta_i (Newton-safe bisection)
  ! ---------------------------------------------------------------------------
  function Delta_i_F_M(sigma, F_M_in, Delta_initial, ierr) result(res)
    real(dp), intent(in)    :: sigma, F_M_in, Delta_initial
    integer,  intent(inout) :: ierr
    real(dp) :: res

    integer,  parameter :: maxit = 1000000
    real(dp), parameter :: tol   = 1.0e-8_dp
    integer  :: j
    real(dp) :: dx, flo, fhi, func, dfunc, dxold

    ierr = 0
    res  = Delta_initial

    if (Delta_initial > 1.0e6_dp) then
      ierr = 1
      return
    end if

    func  = F_M(sigma, res) - F_M_in
    dfunc = P_V(sigma, res) * res

    if (func < 0.0_dp) then
      flo = res
      fhi = 1.0e6_dp
    else
      flo = 0.0_dp
      fhi = res
    end if

    dxold = fhi - flo
    dx    = dxold

    do j = 1, maxit
      if (((res - fhi) * dfunc - func) * ((res - flo) * dfunc - func) >= 0.0_dp .or. &
           abs(2.0_dp * func) > abs(dxold * dfunc)) then
        dxold = dx
        dx    = 0.5_dp * (fhi - flo)
        res   = flo + dx
      else
        dxold = dx
        dx    = func / dfunc
        res   = res - dx
      end if

      if (abs(dx) < tol) return

      func  = F_M(sigma, res) - F_M_in
      dfunc = P_V(sigma, res) * res

      if (func < 0.0_dp) then
        flo = res
      else
        fhi = res
      end if
    end do

    ierr = 1
    write(*, '(a, 4e14.4)') 'F_M: Delta_i exceeded max iterations Delta_initial, Delta_i, F_M(Delta_i), F_M_in: ', &
                              Delta_initial, res, F_M(sigma, res), F_M_in
  end function Delta_i_F_M


  ! ---------------------------------------------------------------------------
  ! Invert F_V to find Delta_i (Newton-safe bisection)
  ! ---------------------------------------------------------------------------
  function Delta_i_F_V(sigma, F_V_in, Delta_initial) result(res)
    real(dp), intent(in) :: sigma, F_V_in, Delta_initial
    real(dp) :: res

    integer,  parameter :: maxit = 100
    real(dp), parameter :: tol   = 1.0e-8_dp
    integer  :: j
    real(dp) :: dx, flo, fhi, func, dfunc, dxold

    res   = Delta_initial
    func  = F_V(sigma, res) - F_V_in
    dfunc = P_V(sigma, res)

    if (func < 0.0_dp) then
      flo = res
      fhi = 1.0e6_dp
    else
      flo = 0.0_dp
      fhi = res
    end if

    dxold = fhi - flo
    dx    = dxold

    do j = 1, maxit
      if (((res - fhi) * dfunc - func) * ((res - flo) * dfunc - func) >= 0.0_dp .or. &
           abs(2.0_dp * func) > abs(dxold * dfunc)) then
        dxold = dx
        dx    = 0.5_dp * (fhi - flo)
        res   = flo + dx
      else
        dxold = dx
        dx    = func / dfunc
        res   = res - dx
      end if

      if (abs(dx) < tol) return

      func  = F_V(sigma, res) - F_V_in
      dfunc = P_V(sigma, res)

      if (func < 0.0_dp) then
        flo = res
      else
        fhi = res
      end if
    end do

    stop 'F_V: Delta_i exceeded maximum iterations'
  end function Delta_i_F_V


  ! ---------------------------------------------------------------------------
  ! Update ionized bubble filling factor
  ! ---------------------------------------------------------------------------
  subroutine update_Q(z, dz, dtimedz_k, Q, Qold, sigma, sigmaold, dnphotdz, &
                      n_e, n_H, R_e_B, ierr)
    real(dp),             intent(in)    :: z, dz, dtimedz_k, sigma, sigmaold, dnphotdz, n_e, n_H, R_e_B
    type(fillingfactor_t), intent(inout) :: Q, Qold
    integer,              intent(inout) :: ierr

    integer  :: j
    real(dp) :: dF_M_dz, flo, fhi, dx, dxold, func, dfunc
    real(dp), parameter :: tol = 1.0e-8_dp

    ierr    = 0

    if (Qold%Q < 1.0_dp .or. Qold%Delta < Delta_H_overlap) then

      Q%Delta   = Delta_H_overlap
      dF_M_dz   = (-log(Q%Delta) + 0.5_dp * sigma**2) / sigma * &
                   P_V(sigma, Q%Delta) * Q%Delta**2 * (sigma - sigmaold) / dz
      Q%F_M     = F_M(sigma, Q%Delta)
      Q%F_V     = F_V(sigma, Q%Delta)
      Q%R       = R(sigma, min(Q%Delta, 1.0e4_dp))
      Q%Q       = (Qold%Q + dz * dnphotdz / (n_H * Mpcbycm**3 * Q%F_M)) / &
                  (1.0_dp + dz * (R_e_B * dtimedz_k * Q%R * n_e * (1.0_dp + z)**3 + dF_M_dz) / Q%F_M)

    else

      Q%Q    = Qold%Q
      Q%F_M  = Qold%F_M
      Q%Delta = Delta_i_F_M(sigma, Q%F_M, Qold%Delta, ierr)
      if (ierr /= 0) return
      Q%R    = R(sigma, min(Q%Delta, 1.0e4_dp))

      flo    = 0.0_dp
      fhi    = 1.0_dp
      dxold  = fhi - flo
      dx     = dxold

      func  = Q%F_M - Qold%F_M + dz * (Q%Q * R_e_B * dtimedz_k * Q%R * n_e * (1.0_dp + z)**3 &
              - dnphotdz / (n_H * Mpcbycm**3))
      dfunc = 1.0_dp + dz * R_e_B * dtimedz_k * n_e * (1.0_dp + z)**3 * Q%Delta

      do j = 1, 100
        if (((Q%F_M - fhi) * dfunc - func) * ((Q%F_M - flo) * dfunc - func) >= 0.0_dp .or. &
             abs(2.0_dp * func) > abs(dxold * dfunc)) then
          dxold = dx
          dx    = 0.5_dp * (fhi - flo)
          Q%F_M = flo + dx
        else
          dxold = dx
          dx    = func / dfunc
          Q%F_M = Q%F_M - dx
        end if

        if (abs(dx) < tol) exit

        Q%Delta = Delta_i_F_M(sigma, Q%F_M, Qold%Delta, ierr)
        if (ierr /= 0) return
        Q%R     = R(sigma, min(Q%Delta, 1.0e4_dp))

        func  = Q%F_M - Qold%F_M + dz * (Q%Q * R_e_B * dtimedz_k * Q%R * n_e * (1.0_dp + z)**3 &
                - dnphotdz / (n_H * Mpcbycm**3))
        dfunc = 1.0_dp + dz * R_e_B * dtimedz_k * n_e * (1.0_dp + z)**3 * Q%Delta

        if (func < 0.0_dp) then
          flo = Q%F_M
        else
          fhi = Q%F_M
        end if
      end do

      Q%Delta = Delta_i_F_M(sigma, Q%F_M, Qold%Delta, ierr)
      if (ierr /= 0) return
      if (Q%Delta < Delta_H_overlap) Q%Delta = Delta_H_overlap - 1.0e-3_dp
      Q%F_M = F_M(sigma, Q%Delta)
      Q%F_V = F_V(sigma, Q%Delta)
      Q%R   = R(sigma, min(Q%Delta, 1.0e4_dp))

    end if

    if (Q%Q > 1.0_dp) Q%Q = 1.0_dp
  end subroutine update_Q

end module inhomoReion_mod
