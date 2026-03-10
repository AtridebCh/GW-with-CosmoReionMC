module state_mod
  use kinds_mod, only:dp
  implicit none
  private

  type, public :: species_t
    real(dp) :: HI, HeI, HeIII
  end type species_t

  type, public :: ionstate_t
    type(species_t) :: frac, recrate, ionrate, coolrate, heatrate
    real(dp) :: T
    real(dp) :: X_HII, X_HeII, X_e, X
  end type ionstate_t

  type, public :: lumsource_t
    real(dp) :: HII, HeII, HeIII
  end type lumsource_t

  type, public :: ionsource_t
    type(lumsource_t) :: pop2, pop3, QSO
  end type ionsource_t

  type, public :: fillingfactor_t
    real(dp) :: Q, Delta, F_V, F_M, R, meanfreepath
  end type fillingfactor_t

end module state_mod

