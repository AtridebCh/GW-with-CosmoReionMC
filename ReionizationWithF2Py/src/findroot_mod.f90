MODULE findroot_mod
implicit none
  DOUBLE PRECISION, DIMENSION(:), POINTER :: fmin_fvecp
CONTAINS
  !BL
  FUNCTION fmin(x)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: x
    DOUBLE PRECISION :: fmin
    INTERFACE
       FUNCTION funcv(x)
         IMPLICIT NONE
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: x
         DOUBLE PRECISION, DIMENSION(size(x)) :: funcv
       END FUNCTION funcv
    END INTERFACE
    IF (.NOT. ASSOCIATED(fmin_fvecp)) THEN 
       PRINT *,'fmin: problem with pointer for returned values'
       STOP
    END IF
    fmin_fvecp=funcv(x)
    fmin=0.5d0*dot_product(fmin_fvecp,fmin_fvecp)
  END FUNCTION fmin

  SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: xold,g
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: p
    DOUBLE PRECISION, INTENT(IN) :: fold,stpmax
    DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: x
    DOUBLE PRECISION, INTENT(OUT) :: f
    LOGICAL, INTENT(OUT) :: check
    INTERFACE
       FUNCTION func(x)
         IMPLICIT NONE
         DOUBLE PRECISION :: func
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE
    DOUBLE PRECISION, PARAMETER :: ALF=1.0d-4
!,TOLX=epsilon(x)
    INTEGER :: ndum
    DOUBLE PRECISION :: a,alam,alam2,alamin,b,disc,f2,fold2,pabs,rhs1,rhs2,slope,&
         tmplam,TOLX
    TOLX=epsilon(x)
    ndum=SIZE(g)
    check=.false.
    !pabs=vabs(p(:))
    pabs=SQRT(DOT_PRODUCT(p(:),p(:)))
    if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
    slope=dot_product(g,p)
    alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0d0))
    alam=1.0
    do
       x(:)=xold(:)+alam*p(:)
       f=func(x)
       if (alam < alamin) then
          x(:)=xold(:)
          check=.true.
          RETURN
       else if (f <= fold+ALF*alam*slope) then
          RETURN
       else
          if (alam == 1.0) then
             tmplam=-slope/(2.0d0*(f-fold-slope))
          else
             rhs1=f-fold-alam*slope
             rhs2=f2-fold2-alam2*slope
             a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
             b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                  (alam-alam2)
             if (a == 0.0) then
                tmplam=-slope/(2.0d0*b)
             else
                disc=b*b-3.0d0*a*slope
                IF (disc < 0.0) THEN 
                   PRINT *,'roundoff problem in lnsrch'
                   STOP
                END IF
                tmplam=(-b+sqrt(disc))/(3.0d0*a)
             end if
             if (tmplam > 0.5d0*alam) tmplam=0.5d0*alam
          end if
       end if
       alam2=alam
       f2=f
       fold2=fold
       alam=max(tmplam,0.1d0*alam)
    end do
  END SUBROUTINE lnsrch

  subroutine newton_solve(x, func, tol, maxits, converged)
  use kinds_mod, only: dp
  implicit none
    !------------------------------------------------------------------
    ! Simple Newton-Raphson solver for a system F(x) = 0
    ! Uses finite-difference Jacobian and Gaussian elimination
    !
    ! x        : on input,  initial guess
    !            on output, solution
    ! tol      : convergence tolerance on |F(x)|
    ! maxits   : maximum number of iterations
    ! converged: .true. if converged, .false. otherwise
    !------------------------------------------------------------------
    real(dp), dimension(:), intent(inout) :: x
    real(dp),               intent(in)    :: tol
    integer,                intent(in)    :: maxits
    logical,                intent(out)   :: converged
    
    interface
    function func(x) result(res)
      import dp
      real(dp), dimension(:), intent(in) :: x
      real(dp), dimension(size(x))       :: res
    end function func
  end interface

    integer  :: n, it, j
    real(dp), dimension(size(x))           :: f, dx, xsav, fph
    real(dp), dimension(size(x),size(x))  :: Jac
    real(dp) :: eps, h

    n   = size(x)
    eps = 1.0e-5_dp

    converged = .false.

    do it = 1, maxits

      ! evaluate residual at current x
      f = func(x)

      ! check convergence
      if (maxval(abs(f)) < tol) then
        converged = .true.
        return
      end if

      ! build finite-difference Jacobian  J(:,j) = dF/dx_j
      xsav = x
      do j = 1, n
        h       = eps * max(abs(xsav(j)), 1.0e-8_dp)
        x       = xsav
        x(j)    = xsav(j) + h
        fph     = func(x)
        Jac(:,j)  = (fph - f) / h
      end do
      x = xsav   ! restore

      ! solve J * dx = -F  via Gaussian elimination with partial pivoting
      dx = -f
      call gauss_solve(Jac, dx, n)

      ! update
      x = x + dx

      ! step-size convergence guard
      if (maxval(abs(dx) / max(abs(x), 1.0e-10_dp)) < epsilon(1.0_dp)) then
        converged = .true.
        return
      end if

    end do

    write(*,*) 'newton_solve: did not converge in', maxits, 'iterations'
    write(*,*) 'Residual: ', maxval(abs(func(x)))

  end subroutine newton_solve


  subroutine gauss_solve(A_, rhs, n)
    use kinds_mod, only: dp
    implicit none
    !------------------------------------------------------------------
    ! Solve A*x = b in-place (b holds solution on exit)
    ! Simple Gaussian elimination with partial pivoting
    !------------------------------------------------------------------
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: A_(n,n), rhs(n)

    integer  :: i, k, ipiv
    real(dp) :: fac, tmp, Atmp(n)

    do k = 1, n
      ! find pivot
      ipiv = k
      do i = k+1, n
        if (abs(A_(i,k)) > abs(A_(ipiv,k))) ipiv = i
      end do

      ! swap rows
      if (ipiv /= k) then
        Atmp     = A_(k,:);  A_(k,:)    = A_(ipiv,:);  A_(ipiv,:)  = Atmp
        tmp      = rhs(k);   rhs(k)      = rhs(ipiv);     rhs(ipiv)    = tmp
      end if

      if (abs(A_(k,k)) < 1.0e-30_dp) then
        write(*,*) 'gauss_solve: singular or near-singular matrix at row', k
        return
      end if

      ! eliminate
      do i = k+1, n
        fac    = A_(i,k) / A_(k,k)
        A_(i,:) = A_(i,:) - fac * A_(k,:)
        rhs(i)   = rhs(i)   - fac * rhs(k)
      end do
    end do

    ! back substitution
    do i = n, 1, -1
      rhs(i) = rhs(i) - dot_product(A_(i,i+1:n), rhs(i+1:n))
      rhs(i) = rhs(i) / A_(i,i)
    end do

  end subroutine gauss_solve


  SUBROUTINE newt(x,check)
    !USE fminln
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: x
    LOGICAL, INTENT(OUT) :: check
    
    INTEGER, PARAMETER :: MAXITS=20000
    DOUBLE PRECISION, PARAMETER :: TOLF=1.0d-8,TOLMIN=1.0d-12,STPMX=100.0
!TOLX=epsilon(x),&
         !STPMX=100.0
    DOUBLE PRECISION :: TOLX
    INTEGER :: its
    INTEGER, DIMENSION(size(x)) :: indx
    DOUBLE PRECISION :: d,f,fold,stpmax
    DOUBLE PRECISION, DIMENSION(size(x)) :: g,p,xold
    DOUBLE PRECISION, DIMENSION(size(x)), TARGET :: fvec
    DOUBLE PRECISION, DIMENSION(size(x),size(x)) :: fjac
    INTERFACE
       SUBROUTINE fdjac(x,fvec,df)
         DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: x
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: fvec
         DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: df
       END SUBROUTINE fdjac
    END INTERFACE

    TOLX=EPSILON(x)
    fmin_fvecp=>fvec
    f=fmin(x)
    if (maxval(abs(fvec(:))) < 0.01d0*TOLF) then
       check=.false.
       RETURN
    end if
    stpmax=STPMX*MAX(SQRT(DOT_PRODUCT(x(:),x(:))),DBLE(SIZE(x)))
    do its=1,MAXITS
       call fdjac(x,fvec,fjac)
       g(:)=matmul(fvec(:),fjac(:,:))
       xold(:)=x(:)
       fold=f
       p(:)=-fvec(:)
       call ludcmp(fjac,indx,d)
       call lubksb(fjac,indx,p)
       call lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin)
       if (maxval(abs(fvec(:))) < TOLF) then
          check=.false.
          RETURN
       end if
       if (check) then
          check=(maxval(abs(g(:))*max(abs(x(:)),1.0d0) / &
               max(f,0.5d0*size(x))) < TOLMIN)
          RETURN
       end if
       if (maxval(abs(x(:)-xold(:))/max(abs(x(:)),1.0d0)) < TOLX) &
            RETURN
    end do
    PRINT *,'MAXITS exceeded in newt'
    STOP
  END SUBROUTINE newt

  SUBROUTINE ludcmp(a,indx,d)

    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER, DIMENSION(:), INTENT(OUT) :: indx
    DOUBLE PRECISION, INTENT(OUT) :: d
    DOUBLE PRECISION, DIMENSION(size(a,1)) :: vv
    DOUBLE PRECISION, PARAMETER :: TINY=1.0d-20
    INTEGER :: j,n,imax
    n=size(a,1)
    d=1.0
    vv=maxval(abs(a),dim=2)
    IF (ANY(vv == 0.d0)) THEN 
       PRINT *,'singular matrix in ludcmp'
       STOP
    END IF
    vv=1.0d0/vv
    do j=1,n
       imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
       if (j /= imax) then
          call swap(a(imax,:),a(j,:))
          d=-d
          vv(imax)=vv(j)
       end if
       indx(j)=imax
       if (a(j,j) == 0.0) a(j,j)=TINY
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
    end do
  END SUBROUTINE ludcmp

  SUBROUTINE lubksb(a,indx,b)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: a
    INTEGER, DIMENSION(:), INTENT(IN) :: indx
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: b
    INTEGER :: i,n,ii,ll
    DOUBLE PRECISION :: summ
    n=size(a,1)
    ii=0
    do i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0) then
          ii=i
       end if
       b(i)=summ
    end do
    do i=n,1,-1
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
  END SUBROUTINE lubksb

  FUNCTION imaxloc(arr)
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: arr
    INTEGER :: imaxloc
    INTEGER, DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc=imax(1)
  END FUNCTION imaxloc
  FUNCTION outerprod(a,b)
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: a,b
    DOUBLE PRECISION, DIMENSION(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod
  SUBROUTINE swap(a,b)
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: a,b
    DOUBLE PRECISION, DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap


END MODULE findroot_mod
