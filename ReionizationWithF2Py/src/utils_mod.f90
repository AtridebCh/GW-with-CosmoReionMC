module utils_mod
  use parameters_mod, only: parameters_t
  use variables_mod,  only: tau_elsc_today, escfrac, ifprint
  implicit none
  private
  public :: summarize, write_summary

contains

  subroutine summarize(params)
    type(parameters_t), intent(in) :: params

    if (ifprint) call write_summary(params, unit=6)
  end subroutine summarize


  subroutine write_summary(params, unit)
    type(parameters_t), intent(in) :: params
    integer,            intent(in) :: unit

    write(unit, '(a)')    'Input parameters'
    write(unit, '(a)')    '----------------'
    write(unit, '(a,f8.4)') ' Omega_m       = ', params%cosmo%omega_m
    write(unit, '(a,f8.4)') ' Omega_lambda  = ', params%cosmo%omega_l
    write(unit, '(a,f8.4)') ' Omega_b h^2   = ', params%cosmo%ombh2
    write(unit, '(a,f8.4)') ' h             = ', params%cosmo%h
    write(unit, '(a,f8.4)') ' n_s           = ', params%cosmo%ns
    write(unit, '(a,f8.4)') ' e_sf_II       = ', params%reion%e_sf_II
    write(unit, '(a,f8.4)') ' escfrac_II    = ', escfrac%pop2%HII
    write(unit, '(a,f8.4)') ' e_sf_III      = ', params%reion%e_sf_III
    write(unit, '(a,f8.4)') ' escfrac_III   = ', escfrac%pop3%HII
    write(unit, '(a,f8.4)') ' lambda_0      = ', params%reion%lambda_0
    write(unit, *)
    write(unit, '(a)')    'Output parameters'
    write(unit, '(a)')    '----------------'
    write(unit, '(a,f8.4)') ' tau_el        = ', tau_elsc_today
  end subroutine write_summary

end module utils_mod
