module thermal_mod
  use kinds_mod,      only: dp
  use constants_mod,  only: kboltz, mprot, yrbysec, mp_by_mHe, Y_He
  use state_mod,      only: ionstate_t, lumsource_t
  use parameters_mod, only: parameters_t, h, omega_m, omega_l, omega_r, &
                            omega_k,omega_b, ombh2, ns, sigma_8, &
                            rho_c
  use variables_mod,  only: delta_c, rho_b_cgs
  use recrates_mod    ! R_HII_e_A, R_HII_e_B, RC_HII_A etc.
  use backgroundCosmology_mod , only: hubbledist, f_integrand
  use findroot_mod,   only: newt, newton_solve
  implicit none
  private
  public :: compute_species, update_ionstate, get_recrates, get_ionrates, &
                      dQ_compton_dt, dumionstate, dumdz, xold   ! add these

  ! working variables for newt callback
  type(ionstate_t)   :: dumionstate
  real(dp)           :: dumdz
  real(dp)           :: xold(3)

contains

  ! ---------------------------------------------------------------------------
  ! Compute derived ionization fractions from primary species fractions
  ! ---------------------------------------------------------------------------
  subroutine compute_species(region)
    type(ionstate_t),  intent(inout) :: region

    region%X_HII  = 1.0_dp - Y_He - region%frac%HI
    region%X_HeII = mp_by_mHe * Y_He - region%frac%HeI - region%frac%HeIII
    region%X_e    = region%X_HII + region%X_HeII + 2.0_dp * region%frac%HeIII
    region%X      = region%frac%HI  + region%X_HII  + &
                    region%frac%HeI + region%X_HeII  + &
                    region%frac%HeIII + region%X_e
  end subroutine compute_species 
  
  
  function dQ_compton_dt(T,z) result(res)
    real(dp), intent(in) :: T, z
    real(dp) :: res                            ! local variable this is fine

    res=6.35e-41_dp*ombh2*(1.0_dp+z)**7*(2.726_dp*(1.0_dp+z)-T)
  end function dQ_compton_dt


  ! ---------------------------------------------------------------------------
  ! Implicit update of ionization state over redshift step dz
  ! ---------------------------------------------------------------------------
  subroutine update_ionstate(z, dz, dtimedz_k, region)
    real(dp),           intent(in)    :: z, dz
    type(ionstate_t),   intent(inout) :: region

    real(dp) :: xp(3), dQdz, dtimedz_k, x1, x2
    logical  :: check


    ! set working variables for funcv callback
    dumdz      = dz
    dumionstate = region

    xp(1) = region%frac%HI
    xp(2) = region%frac%HeI
    xp(3) = region%frac%HeIII
    xold(:)  = xp(:)

    !write(*,*) 'xp before newt', xp
    call newt(xp, check)
    !call newton_solve(xp, funcv, tol=1.0e-10_dp, maxits=50, converged=check)
    !write(*,*)'xp after newton', xp
    
    region%frac%HI    = xp(1)
    region%frac%HeI   = xp(2)
    region%frac%HeIII = xp(3)
    !write(*,*)'Now inside get_ionstate after newton', region%frac%HI, region%frac%HeI, region%frac%HeIII, region%X_HII, region%X_HeII, region%X_e, region%X

    
    
    call compute_species(region)
    

    ! heating and cooling contributions
    dQdz = region%heatrate%HI    * region%frac%HI    + &
           region%heatrate%HeI   * region%frac%HeI   + &
           region%heatrate%HeIII * region%X_HeII     + &
           region%coolrate%HI    * region%X_HII    * region%X_e + &
           region%coolrate%HeI   * region%X_HeII   * region%X_e + &
           region%coolrate%HeIII * region%frac%HeIII * region%X_e


    if (z > 0.0_dp) dQdz = dQdz + dQ_compton_dt(region%T, z) * region%X_e * dtimedz_k
    !write(*,*)'region%T inside update_ionstate', region%T, region%X, region%X_e, region%X_HeII
    x1 = 2.0_dp / (1.0_dp + z) + &
         (1.0_dp / region%X) * (xp(1) - xold(1) + xp(2) - xold(2) - xp(3) + xold(3)) / dz
    x2 = 2.0_dp * mprot / (3.0_dp * kboltz * rho_b_cgs * (1.0_dp + z)**3 * region%X) * dQdz

    region%T = (region%T + dz * x2) / (1.0_dp - dz * x1)
    !write(*,*)'region%T inside update_ionstate', x2, region%X, dQdz
    
  end subroutine update_ionstate


  ! ---------------------------------------------------------------------------
  ! Recombination and cooling rates
  ! ---------------------------------------------------------------------------
  subroutine get_recrates(z, region, dtimedz_k, clumping, Delta_g)
    real(dp),           intent(in)    :: z, dtimedz_k, clumping, Delta_g
    type(ionstate_t),   intent(inout) :: region

    real(dp) ::  n_b

    n_b     = (rho_b_cgs * Delta_g / mprot) * (1.0_dp + z)**3

    region%recrate%HI    =  clumping * (R_HII_e_A(region%T)   + R_HII_e_B(region%T))   * n_b * dtimedz_k
    region%recrate%HeI   =  clumping * (R_HeII_e_A(region%T)  + R_HeII_e_B(region%T))  * n_b * dtimedz_k
    region%recrate%HeIII = -clumping * (R_HeIII_e_A(region%T) + R_HeIII_e_B(region%T)) * n_b * dtimedz_k

    region%coolrate%HI    = -clumping * (RC_HII_A(region%T)   + RC_HII_B(region%T))   * n_b**2 * dtimedz_k / Delta_g
    region%coolrate%HeI   = -clumping * (RC_HeII_A(region%T)  + RC_HeII_B(region%T))  * n_b**2 * dtimedz_k / Delta_g
    region%coolrate%HeIII = -clumping * (RC_HeIII_A(region%T) + RC_HeIII_B(region%T)) * n_b**2 * dtimedz_k / Delta_g
  end subroutine get_recrates


  ! ---------------------------------------------------------------------------
  ! Photoionization and photoheating rates
  ! ---------------------------------------------------------------------------
  subroutine get_ionrates(z, region, dtimedz_k, Gamma_PI, Gamma_PH, Q)
    real(dp),            intent(in)    :: z, Q, dtimedz_k
    type(ionstate_t),    intent(inout) :: region
    type(lumsource_t),   intent(in)    :: Gamma_PI, Gamma_PH

    real(dp) ::  n_b


    n_b     = (rho_b_cgs / mprot) * (1.0_dp + z)**3

    region%ionrate%HI    = -(Gamma_PI%HII   / Q) * dtimedz_k
    region%ionrate%HeI   = -(Gamma_PI%HeII  / Q) * dtimedz_k
    region%ionrate%HeIII =  (Gamma_PI%HeIII / Q) * dtimedz_k

    region%heatrate%HI    = (Gamma_PH%HII   / Q) * n_b * dtimedz_k
    region%heatrate%HeI   = (Gamma_PH%HeII  / Q) * n_b * dtimedz_k
    region%heatrate%HeIII = (Gamma_PH%HeIII / Q) * n_b * dtimedz_k
  end subroutine get_ionrates
  

end module thermal_mod

function funcv(x) result(res)
  use kinds_mod,   only: dp
  use thermal_mod, only: dumionstate, dumdz, xold, compute_species
  use state_mod,   only: ionstate_t
  implicit none
  real(dp), dimension(:), intent(in) :: x
  real(dp), dimension(size(x))       :: res, dxdz 
  
  dumionstate%frac%HI=x(1)
  dumionstate%frac%HeI=x(2)
  dumionstate%frac%HeIII=x(3)
  call compute_species(dumionstate)
  dxdz(1)=dumionstate%ionrate%HI*dumionstate%frac%HI + dumionstate%recrate%HI*dumionstate%X_HII*dumionstate%X_e 

  dxdz(2)=dumionstate%ionrate%HeI*dumionstate%frac%HeI + dumionstate%recrate%HeI*dumionstate%X_HeII*dumionstate%X_e 

  dxdz(3)=dumionstate%ionrate%HeIII*dumionstate%X_HeII + dumionstate%recrate%HeIII*dumionstate%frac%HeIII*dumionstate%X_e 
  
  res(:)=x(:)-xold(:)-dumdz*dxdz(:)
end function funcv


subroutine fdjac(x, fvec, df)
  use kinds_mod, only: dp
  implicit none
  real(dp), dimension(:),   intent(inout) :: x
  real(dp), dimension(:),   intent(in)    :: fvec
  real(dp), dimension(:,:), intent(out)   :: df
  
  interface
    function funcv(x) result(res)
      import dp
      real(dp), dimension(:), intent(in) :: x
      real(dp), dimension(size(x))       :: res
    end function funcv
  end interface
  real(dp), parameter :: EPS=1.0d-4
  integer :: j,n
  real(dp),  dimension(size(x)) :: xsav,xph,h
  
  n=SIZE(x)
  xsav=x
  h=EPS*abs(xsav)
  where (h == 0.0) h=EPS
  xph=xsav+h
  h=xph-xsav
  !write(*,*) 'inside fdjac before the operation', x, funcv(x)
  do j=1,n
     x(j)=xph(j)
     df(:,j)=(funcv(x)-fvec(:))/h(j)
     x(j)=xsav(j)
  end do
  !write(*,*) 'inside fdjac after the operation', x, funcv(x)
end subroutine fdjac

