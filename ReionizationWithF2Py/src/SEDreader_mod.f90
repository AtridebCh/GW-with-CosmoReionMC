module SEDreader_mod
  use kinds_mod,     only: dp
  use constants_mod, only: hPlanck, nu_HII, nu_HeII, nu_HeIII
  use state_mod,     only: lumsource_t, ionsource_t
  use variables_mod, only: dnphotdm, sigma_PH, sigma_PI, ifprint, Pop3IMF
  use crossSection_mod !sigma_HI(nu), sigma_HeII(nu), etc
  implicit none
  private
  public :: get_SED

contains
  
  subroutine get_SED()
    type(lumsource_t) :: dndm2, dndm3, sPI2, sPI3, sPH2, sPH3

    ! zero all
    dnphotdm = ionsource_t(lumsource_t(0.0_dp, 0.0_dp, 0.0_dp), &
                         lumsource_t(0.0_dp, 0.0_dp, 0.0_dp), &
                         lumsource_t(0.0_dp, 0.0_dp, 0.0_dp))
    sigma_PH = ionsource_t(lumsource_t(0.0_dp, 0.0_dp, 0.0_dp), &
                         lumsource_t(0.0_dp, 0.0_dp, 0.0_dp), &
                         lumsource_t(0.0_dp, 0.0_dp, 0.0_dp))
    sigma_PI = ionsource_t(lumsource_t(0.0_dp, 0.0_dp, 0.0_dp), &
                         lumsource_t(0.0_dp, 0.0_dp, 0.0_dp), &
                         lumsource_t(0.0_dp, 0.0_dp, 0.0_dp))

    call get_stellar_SED("/Users/users/achatterjee/Data/GW_with_CosmoReionMC/ReionizationWithF2Py/SpecData/dndm_PopII_salpeter_Z0.004", dndm2, sPI2, sPH2)
    call get_stellar_SED("/Users/users/achatterjee/Data/GW_with_CosmoReionMC/ReionizationWithF2Py/SpecData/dndm_PopIII_" // trim(Pop3IMF) // "_star", dndm3, sPI3, sPH3)

    dnphotdm%pop2 = dndm2;  sigma_PI%pop2 = sPI2;  sigma_PH%pop2 = sPH2
    dnphotdm%pop3 = dndm3;  sigma_PI%pop3 = sPI3;  sigma_PH%pop3 = sPH3

    call get_quasar_SED(sigma_PI%QSO, sigma_PH%QSO)
  end subroutine get_SED


  subroutine get_stellar_SED(sed_file, dnphotdm, sigma_PI, sigma_PH)
    character(len=*),  intent(in)  :: sed_file
    type(lumsource_t), intent(out) :: dnphotdm, sigma_PI, sigma_PH

    real(dp) :: nu, dndm, nu_old, dndm_old
    integer  :: io_unit, io_stat
    
    
    dnphotdm = lumsource_t(0.0_dp, 0.0_dp, 0.0_dp)
    sigma_PI = lumsource_t(0.0_dp, 0.0_dp, 0.0_dp)
    sigma_PH = lumsource_t(0.0_dp, 0.0_dp, 0.0_dp)
    

    open(newunit=io_unit, file=sed_file, status='old', action='read')
    !print *, "Reading file:", trim(sed_file)
    read(io_unit, *, iostat=io_stat)
    read(io_unit, *, iostat=io_stat)
    read(io_unit, *, iostat=io_stat) nu_old, dndm_old
    
    if (io_stat /= 0) then
      print *, "Error reading first data line in", trim(sed_file)
      stop
    end if

    do
      read(io_unit, *, iostat=io_stat) nu, dndm
      if (io_stat /= 0) exit

      ! photon count integrals
      if (nu >= nu_HII  .and. nu < nu_HeII)  dnphotdm%HII   = dnphotdm%HII   + (nu_old - nu) * 0.5_dp * (dndm + dndm_old)
      if (nu >= nu_HeII .and. nu < nu_HeIII) dnphotdm%HeII  = dnphotdm%HeII  + (nu_old - nu) * 0.5_dp * (dndm + dndm_old)
      if (nu >= nu_HeIII)                    dnphotdm%HeIII = dnphotdm%HeIII + (nu_old - nu) * 0.5_dp * (dndm + dndm_old)

      ! HI cross section integrals
      if (nu >= nu_HII) then
        if (nu >= nu_HeIII) then
          sigma_PI%HII = sigma_PI%HII + sigma_HI(nu) * 0.5_dp * (dndm + dndm_old) * (nu_old - nu) * (nu / nu_HeIII)**1.5_dp
          sigma_PH%HII = sigma_PH%HII + sigma_HI(nu) * hPlanck * (nu - nu_HII) * 0.5_dp * (dndm + dndm_old) * (nu_old - nu) * (nu / nu_HeIII)**1.5_dp
        else
          sigma_PI%HII = sigma_PI%HII + sigma_HI(nu) * 0.5_dp * (dndm + dndm_old) * (nu_old - nu)
          sigma_PH%HII = sigma_PH%HII + sigma_HI(nu) * hPlanck * (nu - nu_HII) * 0.5_dp * (dndm + dndm_old) * (nu_old - nu)
        end if
      end if

      ! HeI cross section integrals
      if (nu >= nu_HeII) then
        if (nu >= nu_HeIII) then
          sigma_PI%HeII = sigma_PI%HeII + sigma_HeI(nu) * 0.5_dp * (dndm + dndm_old) * (nu_old - nu) * (nu / nu_HeIII)**1.5_dp
          sigma_PH%HeII = sigma_PH%HeII + sigma_HeI(nu) * hPlanck * (nu - nu_HeII) * 0.5_dp * (dndm + dndm_old) * (nu_old - nu) * (nu / nu_HeIII)**1.5_dp
        else
          sigma_PI%HeII = sigma_PI%HeII + sigma_HeI(nu) * 0.5_dp * (dndm + dndm_old) * (nu_old - nu)
          sigma_PH%HeII = sigma_PH%HeII + sigma_HeI(nu) * hPlanck * (nu - nu_HeII) * 0.5_dp * (dndm + dndm_old) * (nu_old - nu)
        end if
      end if

      ! HeII cross section integrals
      if (nu >= nu_HeIII) then
        sigma_PI%HeIII = sigma_PI%HeIII + sigma_HeII(nu) * 0.5_dp * (dndm + dndm_old) * (nu_old - nu) * (nu / nu_HeIII)**1.5_dp
        sigma_PH%HeIII = sigma_PH%HeIII + sigma_HeII(nu) * hPlanck * (nu - nu_HeIII) * 0.5_dp * (dndm + dndm_old) * (nu_old - nu) * (nu / nu_HeIII)**1.5_dp
      end if

      nu_old   = nu
      dndm_old = dndm
    end do

    close(io_unit)

    sigma_PH%HII   = sigma_PH%HII   / (dnphotdm%HII   * hPlanck * nu_HII)
    sigma_PH%HeII  = sigma_PH%HeII  / (dnphotdm%HeII  * hPlanck * nu_HeII)
    sigma_PH%HeIII = sigma_PH%HeIII / (dnphotdm%HeIII * hPlanck * nu_HeIII)

    sigma_PI%HII   = sigma_PI%HII   / dnphotdm%HII
    sigma_PI%HeII  = sigma_PI%HeII  / dnphotdm%HeII
    sigma_PI%HeIII = sigma_PI%HeIII / dnphotdm%HeIII

  end subroutine get_stellar_SED


  subroutine get_quasar_SED(sigma_PI, sigma_PH)
    type(lumsource_t), intent(inout) :: sigma_PI, sigma_PH

    real(dp), parameter :: alpha = 1.57_dp

    sigma_PH%HII = (alpha / nu_HII) * ( &
      sigma_integral_HI(-alpha,        nu_HII,         nu_HeIII)        * nu_HII**alpha - &
      sigma_integral_HI(-(alpha+1.0_dp), nu_HII,       nu_HeIII)        * nu_HII**(alpha+1.0_dp) + &
      sigma_integral_HI(-(alpha-1.5_dp), nu_HeIII,     nu_HeIII*40.0_dp) * nu_HII**alpha * nu_HeIII**(-1.5_dp) - &
      sigma_integral_HI(-1.07_dp,          nu_HeIII,     nu_HeIII*40.0_dp) * nu_HII**(alpha+1.0_dp) * nu_HeIII**(-1.5_dp)) / &
      (1.0_dp - (nu_HII / nu_HeII)**alpha)

    sigma_PI%HII = alpha * ( &
      sigma_integral_HI(-(alpha+1.0_dp), nu_HII,       nu_HeIII)        * nu_HII**alpha + &
      sigma_integral_HI(-1.07_dp,          nu_HeIII,     nu_HeIII*40.0_dp) * nu_HII**alpha * nu_HeIII**(-1.5_dp)) / &
      (1.0_dp - (nu_HII / nu_HeII)**alpha)
         
         
    sigma_PH%HeII = (alpha / nu_HeII) * ( &
      sigma_integral_HeI(-alpha,           nu_HeII,    nu_HeIII)         * nu_HII**alpha - &
      sigma_integral_HeI(-(alpha+1.0_dp),  nu_HeII,    nu_HeIII)         * nu_HII**alpha * nu_HeII + &
      sigma_integral_HeI(-0.07_dp, nu_HeIII,   nu_HeIII*40.0_dp) * nu_HII**alpha * nu_HeIII**(-1.5_dp) - &
      sigma_integral_HeI(-1.07_dp, nu_HeIII,   nu_HeIII*40.0_dp) * nu_HII**alpha * nu_HeII * nu_HeIII**(-1.5_dp)) / &
      ((nu_HII / nu_HeII)**alpha - (nu_HII / nu_HeIII)**alpha)



    sigma_PI%HeII = alpha * ( &
      sigma_integral_HeI(-(alpha+1.0_dp),  nu_HeII,    nu_HeIII)         * nu_HII**alpha + &
      sigma_integral_HeI(-alpha,           nu_HeIII,   nu_HeIII*40.0_dp) * nu_HII**1.07 * nu_HeIII**(-1.5_dp)) / &
      ((nu_HII / nu_HeII)**alpha - (nu_HII / nu_HeIII)**alpha)

    sigma_PH%HeIII = (alpha / nu_HeIII) * ( &
      sigma_integral_HeII(-0.07_dp, nu_HeIII,   nu_HeIII*40.0_dp) * nu_HII**alpha * nu_HeIII**(-1.5_dp) - &
      sigma_integral_HeII(-1.07_dp, nu_HeIII,   nu_HeIII*40.0_dp) * nu_HII**alpha * nu_HeIII**(-0.5_dp)) / &
      (nu_HII / nu_HeIII)**alpha

    sigma_PI%HeIII = alpha * &
      sigma_integral_HeII(-1.07_dp, nu_HeIII, nu_HeIII*40.0_dp) * nu_HII**alpha * nu_HeIII**(-1.5_dp) / &
      (nu_HII / nu_HeIII)**alpha

  end subroutine get_quasar_SED


end module SEDreader_mod
