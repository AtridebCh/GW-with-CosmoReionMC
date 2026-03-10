
    module Reionization
    use Precision
    use MiscUtils
    use classes
    use results
    implicit none
    private
    !This module puts smooth tanh reionization of specified mid-point (z_{re}) and width
    !The tanh function is in the variable (1+z)**Rionization_zexp

    !Rionization_zexp=1.5 has the property that for the same z_{re}
    !the optical depth agrees with infinitely sharp model for matter domination
    !So tau and zre can be mapped into each other easily (for any symmetric window)
    !However for generality the module maps tau into z_{re} using a binary search
    !so could be easily modified for other monatonic parameterizations.

    !The ionization history must be twice differentiable.

    !AL March 2008
    !AL July 2008 - added trap for setting optical depth without use_optical_depth

    !See CAMB notes for further discussion: http://cosmologist.info/notes/CAMB.pdf

    real(dl), parameter :: TTanhReionization_DefFraction = -1._dl
    !if -1 set from YHe assuming Hydrogen and first ionization of Helium follow each other

    real(dl) :: Tanh_zexp = 1.5_dl


    type, extends(TReionizationModel) :: TTanhReionization
        logical    :: use_optical_depth = .false.
!rpc
        logical    :: UsePCReion=.false.
        integer    :: array_length = 251
!end_rpc
        real(dl)   :: redshift = 10._dl
        real(dl)   :: optical_depth = 0._dl
        real(dl)   :: delta_redshift = 0.5_dl
        real(dl)   :: fraction = TTanhReionization_DefFraction
        !Parameters for the second reionization of Helium
        logical    :: include_helium_fullreion  = .true.
        real(dl)   :: helium_redshift  = 3.5_dl
        real(dl)   :: helium_delta_redshift  = 0.4_dl
        real(dl)   :: helium_redshiftstart  = 5.5_dl
        real(dl)   :: tau_solve_accuracy_boost = 1._dl
        real(dl)   :: timestep_boost =  1._dl
        real(dl)   :: max_redshift = 50._dl
        !The rest are internal to this module.
        real(dl), private ::  fHe, WindowVarMid, WindowVarDelta
        class(CAMBdata), pointer :: State
        real(dl), private, allocatable :: reion_array(:),red_array(:) !changed here
        
        
    contains
    procedure :: ReadParams => TTanhReionization_ReadParams
    procedure :: Init => TTanhReionization_Init
    procedure :: x_e => TTanhReionization_xe
    procedure :: get_timesteps => TTanhReionization_get_timesteps
    procedure, nopass :: SelfPointer => TTanhReionization_SelfPointer
    procedure, nopass ::  GetZreFromTau => TTanhReionization_GetZreFromTau
    procedure, private :: SetParamsForZre => TTanhReionization_SetParamsForZre
    procedure, private :: zreFromOptDepth => TTanhReionization_zreFromOptDepth

    end type TTanhReionization

    public TTanhReionization
    contains


    function TTanhReionization_xe(this, z, tau, xe_recomb)
    !a and time tau and redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    class(TTanhReionization) :: this
    real(dl), intent(in) :: z
    real(dl), intent(in), optional :: tau, xe_recomb
    real(dl) TTanhReionization_xe
    real(dl):: tgh, xod, dum
    real(dl) :: xstart,A
    !RPC
  integer :: rsbin
!END RPC
    integer :: ier,i
    xstart = PresentDefault( 0._dl, xe_recomb)

    xod = (this%WindowVarMid - (1+z)**Tanh_zexp)/this%WindowVarDelta
    if (xod > 100) then
        tgh=1.d0
    else
        tgh=tanh(xod)
    end if
    TTanhReionization_xe =(this%fraction-xstart)*(tgh+1._dl)/2._dl+xstart
    
    if (this%UsePCReion) then
       if (.not. allocated(this%red_array)) then
        TTanhReionization_xe = xstart
       return
       end if
    end if
    

    !RPC!Use linear interpolation between redshift bins
        if (this%UsePCReion) then
           do rsbin = 1, this%array_length
              if (this%reion_array(rsbin).eq.0._dl) then
                 this%reion_array(rsbin) = max(xstart,1.d-10)
                endif
           end do
            !write(*,*)'here', this%reion_array(1)
            !write(*,*) 'here',2.0*this%red_array(this%array_length)- this%red_array(this%array_length-1))
           rsbin = 1
           if (z < (2.*this%red_array(1)-this%red_array(2)) .and. &
                z > (2.*this%red_array(this%array_length)- &
                this%red_array(this%array_length-1))) then
              do while (this%red_array(rsbin)>z .and. rsbin/=this%array_length)
                 rsbin = rsbin+1
              end do
                
              if (rsbin==1) then
                 TTanhReionization_xe= xstart+(this%reion_array(1)-xstart)* &
                      (2.*this%red_array(1)-this%red_array(2)-z)/ &
                      (this%red_array(1)-this%red_array(2))+this%fHe*((this%reion_array(1)-xstart)* &
                      (2.*this%red_array(1)-this%red_array(2)-z)/ &
                      (this%red_array(1)-this%red_array(2))) !changes made after Kyungjin's suggestion ; For the original version see the download folder
                  
              else if (rsbin==this%array_length .and. this%red_array(rsbin)>z) then
                 TTanhReionization_xe= this%reion_array(rsbin) + &
                      (1._dl+this%fHe-this%reion_array(rsbin))* &
                      (this%red_array(rsbin)-z)/(this%red_array(this%array_length-1)- &
                      this%red_array(this%array_length))+this%fHe*(this%reion_array(rsbin) + &
                      (1._dl+this%fHe-this%reion_array(rsbin))* &
                      (this%red_array(rsbin)-z)/(this%red_array(this%array_length-1)- &
                      this%red_array(this%array_length))) !changes made after Kyungjin's suggestion ; For the original version see the download folder
                   
              else
                 TTanhReionization_xe = this%reion_array(rsbin-1) + (this%reion_array(rsbin)-this%reion_array(rsbin-1))* &
                      (this%red_array(rsbin-1)-z)/(this%red_array(rsbin-1)-this%red_array(rsbin))+this%fHe*(this%reion_array(rsbin-1) + (this%reion_array(rsbin)-this%reion_array(rsbin-1))* &
                      (this%red_array(rsbin-1)-z)/(this%red_array(rsbin-1)-this%red_array(rsbin))) !changes made after Kyungjin's suggestion ; For the original version see the download folder
                 
              end if
           else if (z<=(2.*this%red_array(this%array_length)- &
                this%red_array(this%array_length-1))) then
              TTanhReionization_xe = 1._dl + this%fHe
              
           end if
               
        end if
!END RPC



    if (this%include_helium_fullreion .and. z < this%helium_redshiftstart) then

        !Effect of Helium becoming fully ionized is small so details not important
        xod = (this%helium_redshift - z)/this%helium_delta_redshift
        if (xod > 100) then
            tgh=1.d0
        else
            tgh=tanh(xod)
        end if

        TTanhReionization_xe =TTanhReionization_xe+this%fHe*(tgh+1._dl)/2._dl
        !write(*,*)TTanhReionization_xe
    end if
    
    !if (TTanhReionization_xe>1.0) then
     !   TTanhReionization_xe=1.0
    !endif
    if (TTanhReionization_xe<0.0) then
        TTanhReionization_xe=0.0
    endif
    if (TTanhReionization_xe>1.18) then
    	TTanhReionization_xe=1._dl+2.0*this%fHe !changes made after Kyungjin's suggestion ; For the original version see the download folder
    endif
    end function TTanhReionization_xe
    
    FUNCTION splevl(xb,x,f,c,dfb,ddfb,ier)

    use Precision
    IMPLICIT NONE

    REAL(dl), INTENT(in) :: xb
    REAL(dl), DIMENSION(:), INTENT(in) :: x,f
    REAL(dl), DIMENSION(3,SIZE(x)), INTENT(in) :: c
    INTEGER, INTENT(inout) :: ier
    REAL(dl), INTENT(out) :: dfb,ddfb
    REAL(dl) :: splevl
    INTEGER :: n,high,low,nhigh,mid
    REAL(dl) :: dx
    LOGICAL :: qascnd

    n=SIZE(x)
    low=0
    IF(n <= 1) THEN
       ier=111
       splevl=0.d0
       RETURN
    ENDIF

    qascnd=x(n) > x(1)
    ier=0

    IF((low < 1).OR.(low >= n)) THEN
       low=1
       high=n
    ENDIF

1   IF((xb < x(low)).EQV.(xb < x(high))) THEN
       IF((xb > x(low)).EQV.qascnd) THEN
          IF(high >= n) THEN
             ier=11
             low=n-1
          ELSE
             nhigh=MIN(n,high+2*(high-low))
             low=high
             high=nhigh
             GOTO 1
          ENDIF
       ELSE
          IF(low <= 1) THEN
             ier=12
          ELSE
             nhigh=low
             low=MAX(1,low-2*(high-low))
             high=nhigh
             GOTO 1
          ENDIF
       ENDIF
    ELSE
2      IF((high-low) > 1) THEN
          mid=(low+high)/2
          IF((xb < x(mid)).EQV.(xb < x(low))) THEN
             low=mid
          ELSE
             high=mid
          ENDIF
          GOTO 2
       ENDIF
    ENDIF

    dx=xb-x(low)
    splevl=((c(3,low)*dx+c(2,low))*dx+c(1,low))*dx+f(low)
    dfb=(3.d0*c(3,low)*dx+2.d0*c(2,low))*dx+c(1,low)
    ddfb=6.d0*c(3,low)*dx+2.d0*c(2,low)

  END FUNCTION splevl

    subroutine TTanhReionization_get_timesteps(this, n_steps, z_start, z_complete)
    !minimum number of time steps to use between tau_start and tau_complete
    !Scaled by AccuracyBoost later
    !steps may be set smaller than this anyway
    class(TTanhReionization) :: this
    !class(TCAMBdata), target :: State
    integer, intent(out) :: n_steps
    real(dl), intent(out):: z_start, z_Complete

    n_steps = nint(50.0)
    z_start = this%redshift + this%delta_redshift*8
    z_complete = max(0.d0,this%redshift-this%delta_redshift*8)

    end subroutine TTanhReionization_get_timesteps

    subroutine TTanhReionization_ReadParams(this, Ini)
    use IniObjects
    class(TTanhReionization) :: this
    class(TIniFile), intent(in) :: Ini

    this%Reionization = Ini%Read_Logical('reionization')
    if (this%Reionization) then

        this%use_optical_depth = Ini%Read_Logical('re_use_optical_depth')

        if (this%use_optical_depth) then
            this%optical_depth = Ini%Read_Double('re_optical_depth')
        else
            this%redshift = Ini%Read_Double('re_redshift')
        end if

        call Ini%Read('re_delta_redshift',this%delta_redshift)
        call Ini%Read('re_ionization_frac',this%fraction)
        call Ini%Read('re_helium_redshift',this%helium_redshift)
        call Ini%Read('re_helium_delta_redshift',this%helium_delta_redshift)

        this%helium_redshiftstart  = Ini%Read_Double('re_helium_redshiftstart', &
            this%helium_redshift + 5*this%helium_delta_redshift)

    end if

    end subroutine TTanhReionization_ReadParams

    subroutine TTanhReionization_SetParamsForZre(this)
    class(TTanhReionization) :: this

    this%WindowVarMid = (1._dl+this%redshift)**Tanh_zexp
    this%WindowVarDelta = Tanh_zexp*(1._dl+this%redshift)**(Tanh_zexp-1._dl)*this%delta_redshift

    end subroutine TTanhReionization_SetParamsForZre

    subroutine TTanhReionization_Init(this, State)
    use constants
    use MathUtils
    
    real(dl), allocatable :: reionization_array(:)
    real(dl), allocatable :: redshift_array(:)
    integer :: reion_array_length, n

    class(TTanhReionization) :: this
    
    class(TCAMBdata), target :: State
    procedure(obj_function) :: dtauda

    select type (State)
    class is (CAMBdata)
        this%State => State
        
        ! ====================================================
        ! External reionization injection (your custom block)
        ! ====================================================
        if (this%State%CP%ReionExternal) then

           this%array_length = this%State%CP%length_array

           if (allocated(this%red_array)) deallocate(this%red_array)
           if (allocated(this%reion_array)) deallocate(this%reion_array)

           allocate(this%red_array(this%array_length))
           allocate(this%reion_array(this%array_length))
           
           !I need to check it for taking care of the seg fault

           this%red_array = this%State%CP%redshift_array_external
           this%reion_array = this%State%CP%reionization_history

           this%UsePCReion = .true.
        end if
        
        this%fHe =  State%CP%YHe/(mass_ratio_He_H*(1.d0-State%CP%YHe))
!rpc
        if (this%Reionization.and.(.not.this%UsePCReion)) then
!end rpc
            if (this%optical_depth /= 0._dl .and. .not. this%use_optical_depth) &
                write (*,*) 'WARNING: You seem to have set the optical depth, but use_optical_depth = F'

            if (this%use_optical_depth.and.this%optical_depth<0.001 &
                .or. .not.this%use_optical_depth .and. this%Redshift<0.001) then
                this%Reionization = .false.
            end if

        end if

        if (this%Reionization) then
			this%redshift = 0
!rpc
			if (this%UsePCReion) then
            	this%redshift=this%red_array(1)
!end rpc
			else 
            	if (this%fraction==TTanhReionization_DefFraction) &
                	this%fraction = 1._dl + this%fHe  !H + singly ionized He

            	if (this%use_optical_depth) then
                	
                	call this%zreFromOptDepth()
				endif
                if (global_error_flag/=0) return
                !if (FeedbackLevel > 0) write(*,'("Reion redshift       =  ",f6.3)') this%redshift
                !RPC
!rpc        if (FeedbackLevel > 0) write(*,'("Reion redshift       =  ",f6.3)') Reion%redshift
			
        if ((.not.this%UsePCReion).and.(FeedbackLevel > 0)) &
             write(*,'("Reion redshift       =  ",f6.3)') this%redshift
!END RPC
            end if

            call this%SetParamsForZre()

            !this is a check, agrees very well in default parameterization
            if (FeedbackLevel > 1) write(*,'("Integrated opt depth = ",f7.4)') this%State%GetReionizationOptDepth()

        end if
    end select
    end subroutine TTanhReionization_Init

    subroutine TTanhReionization_Validate(this, OK)
    class(TTanhReionization),intent(in) :: this
    logical, intent(inout) :: OK
!rpc
   if (this%Reionization.and.(.not.this%UsePCReion)) then
!end rpc
        if (this%use_optical_depth) then
            if (this%optical_depth<0 .or. this%optical_depth > 0.9  .or. &
                this%include_helium_fullreion .and. this%optical_depth<0.01) then
                OK = .false.
                write(*,*) 'Optical depth is strange. You have:', this%optical_depth
            end if
        else
            if (this%redshift < 0 .or. this%Redshift +this%delta_redshift*3 > this%max_redshift .or. &
                this%include_helium_fullreion .and. this%redshift < this%helium_redshift) then
                OK = .false.
                write(*,*) 'Reionization redshift strange. You have: ',this%Redshift
            end if
        end if
        if (this%fraction/= TTanhReionization_DefFraction .and. (this%fraction < 0 .or. this%fraction > 1.5)) then
            OK = .false.
            write(*,*) 'Reionization fraction strange. You have: ',this%fraction
        end if
        if (this%delta_redshift > 3 .or. this%delta_redshift<0.1 ) then
            !Very narrow windows likely to cause problems in interpolation etc.
            !Very broad likely to conflict with quasar data at z=6
            OK = .false.
            write(*,*) 'Reionization delta_redshift is strange. You have: ',this%delta_redshift
        end if
    end if

    end subroutine TTanhReionization_Validate

    subroutine TTanhReionization_zreFromOptDepth(this)
    !General routine to find zre parameter given optical depth
    class(TTanhReionization) :: this
    real(dl) try_b, try_t
    real(dl) tau, last_top, last_bot
    integer i

    try_b = 0
    try_t = this%max_redshift
    i=0
    do
        i=i+1
        this%redshift = (try_t + try_b)/2
        call this%SetParamsForZre()
        tau = this%State%GetReionizationOptDepth()

        if (tau > this%optical_depth) then
            try_t = this%redshift
            last_top = tau
        else
            try_b = this%redshift
            last_bot = tau
        end if
        if (abs(try_b - try_t) < 1e-2_dl/this%tau_solve_accuracy_boost) then
            if (try_b==0) last_bot = 0
            if (try_t/=this%max_redshift) this%redshift  = &
                (try_t*(this%optical_depth-last_bot) + try_b*(last_top-this%optical_depth))/(last_top-last_bot)
            exit
        end if
        if (i>100) call GlobalError('TTanhReionization_zreFromOptDepth: failed to converge',error_reionization)
    end do

    if (abs(tau - this%optical_depth) > 0.002 .and. global_error_flag==0) then
        write (*,*) 'TTanhReionization_zreFromOptDepth: Did not converge to optical depth'
        write (*,*) 'tau =',tau, 'optical_depth = ', this%optical_depth
        write (*,*) try_t, try_b
        write (*,*) '(If running a chain, have you put a constraint on tau?)'
        call GlobalError('Reionization did not converge to optical depth',error_reionization)
    end if

    end subroutine TTanhReionization_zreFromOptDepth

    real(dl) function TTanhReionization_GetZreFromTau(P, tau)
    type(CAMBparams) :: P, P2
    real(dl) tau
    integer error
    type(CAMBdata) :: State

    P2 = P

    select type(Reion=>P2%Reion)
    class is (TTanhReionization)
        Reion%Reionization = .true.
        Reion%use_optical_depth = .true.
        Reion%optical_depth = tau
    end select
    call State%SetParams(P2,error)
    if (error/=0)  then
        TTanhReionization_GetZreFromTau = -1
    else
        select type(Reion=>State%CP%Reion)
        class is (TTanhReionization)
            TTanhReionization_GetZreFromTau = Reion%redshift
        end select
    end if

    end function TTanhReionization_GetZreFromTau

    subroutine TTanhReionization_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TTanhReionization), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TTanhReionization_SelfPointer

    end module Reionization
