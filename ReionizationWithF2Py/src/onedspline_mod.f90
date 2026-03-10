MODULE onedspline_mod
public :: spline, splevl

CONTAINS

  SUBROUTINE spline(x,f,c,ier)

    USE kinds_mod
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(in) :: x,f
    INTEGER, INTENT(inout) ::  ier
    REAL(DP), DIMENSION(3,SIZE(x)), INTENT(out) :: c
    INTEGER :: i,j,n
    REAL(DP) :: div12,div23,c1,cn,g

    ier=0
    n=SIZE(x)
    IF(n <= 1) THEN
       ier=111
       RETURN

    ELSE IF(n == 2) THEN
       c(1,1)=(f(2)-f(1))/(x(2)-x(1))
       c(2,1)=0.d0
       c(3,1)=0.d0
       RETURN

    ELSE IF(n == 3) THEN
       div12=(f(2)-f(1))/(x(2)-x(1))
       div23=(f(3)-f(2))/(x(3)-x(2))
       c(3,1)=0.d0
       c(3,2)=0.d0
       c(2,1)=(div23-div12)/(x(3)-x(1))
       c(2,2)=c(2,1)
       c(1,1)=div12+c(2,1)*(x(1)-x(2))
       c(1,2)=div23+c(2,1)*(x(2)-x(3))
       RETURN

    ELSE
       c(3,n)=(f(n)-f(n-1))/(x(n)-x(n-1))
       DO i=n-1,2,-1
          c(3,i)=(f(i)-f(i-1))/(x(i)-x(i-1))
          c(2,i)=2.d0*(x(i+1)-x(i-1))
          c(1,i)=3.d0*(c(3,i)*(x(i+1)-x(i))+c(3,i+1)*(x(i)-x(i-1)))
       ENDDO
       c1=x(3)-x(1)
       c(2,1)=x(3)-x(2)
       c(1,1)=c(3,2)*c(2,1)*(2.d0*c1+x(2)-x(1))+c(3,3)*(x(2)-x(1))**2
       c(1,1)=c(1,1)/c1
       cn=x(n)-x(n-2)
       c(2,n)=x(n-1)-x(n-2)
       c(1,n)=c(3,n)*c(2,n)*(2.d0*cn+x(n)-x(n-1))
       c(1,n)=(c(1,n)+c(3,n-1)*(x(n)-x(n-1))**2)/cn

       g=(x(3)-x(2))/c(2,1)
       c(2,2)=c(2,2)-g*c1
       c(1,2)=c(1,2)-g*c(1,1)
       DO j=2,n-2
          g=(x(j+2)-x(j+1))/c(2,j)
          c(2,j+1)=c(2,j+1)-g*(x(j)-x(j-1))
          c(1,j+1)=c(1,j+1)-g*c(1,j)
       ENDDO
       g=cn/c(2,n-1)
       c(2,n)=c(2,n)-g*(x(n-1)-x(n-2))
       c(1,n)=c(1,n)-g*c(1,n-1)

       c(1,n)=c(1,n)/c(2,n)
       DO i=n-1,2,-1
          c(1,i)=(c(1,i)-c(1,i+1)*(x(i)-x(i-1)))/c(2,i)
       ENDDO
       c(1,1)=(c(1,1)-c(1,2)*c1)/c(2,1)

       DO i=1,n-1
          c(2,i)=(3.d0*c(3,i+1)-2.d0*c(1,i)-c(1,i+1))/(x(i+1)-x(i))
          c(3,i)=(c(1,i)+c(1,i+1)-2.d0*c(3,i+1))/(x(i+1)-x(i))**2
       ENDDO
       c(2,n)=0.d0
       c(3,n)=0.d0
    ENDIF

    RETURN
  END SUBROUTINE spline


  FUNCTION splevl(xb,x,f,c,dfb,ddfb,ier)

    USE kinds_mod
    IMPLICIT NONE

    REAL(DP), INTENT(in) :: xb
    REAL(DP), DIMENSION(:), INTENT(in) :: x,f
    REAL(DP), DIMENSION(3,SIZE(x)), INTENT(in) :: c
    INTEGER, INTENT(inout) :: ier
    REAL(DP), INTENT(out) :: dfb,ddfb
    REAL(DP) :: splevl
    INTEGER :: n,high,low,nhigh,mid
    REAL(DP) :: dx
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
END MODULE onedspline_mod
