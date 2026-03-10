PROGRAM binned

  IMPLICIT NONE

  INTEGER :: i,j,n,iend,nbin,pgopen
  REAL, DIMENSION(1000) :: z,flux=0.,fluxmin=0.,fluxmax=0.
  REAL :: dzbin,zmin,zmax
  REAL :: dum
  REAL, DIMENSION(:), ALLOCATABLE :: zmean,zhi,zlo,fluxmean,err1max,err2max,err1min,err2min
  INTEGER, DIMENSION(:), ALLOCATABLE :: num

  dzbin=0.2
  zmin=2.
  zmax=6.
  nbin=INT((zmax-zmin)/dzbin+1.5)
  ALLOCATE(zmean(nbin),zhi(nbin),zlo(nbin),fluxmean(nbin),err1max(nbin),err2max(nbin),err1min(nbin),err2min(nbin),num(nbin))
  fluxmean=0.
  err1max=0.0
  err2max=0.0
  err1min=0.0
  err2min=0.0
  DO i=1,nbin
     zmean(i)=(i-1)*dzbin+zmin
     zhi(i)=zmean(i)+0.5*dzbin
     zlo(i)=zmean(i)-0.5*dzbin
     PRINT *,zmean(i),zhi(i),zlo(i)
  END DO

  IF (pgopen('binnedflux.ps/cps').LT.1) STOP
  CALL pgscf(2)
  CALL pgpap(0.0,1.0)

  CALL pgenv(2.0,6.5,-1.0,1.5,0,20)

  i=1
  OPEN(2,file='lyaflux_fan',status='old',action='read')
  DO
     READ(2,*,iostat=iend) z(i),flux(i),dum
     IF (iend /= 0) EXIT
     fluxmin(i)=flux(i)-dum
     fluxmax(i)=flux(i)+dum
!!$     IF (flux(i) > 0.0) THEN
!!$        flux(i)=-LOG(flux(i))
!!$        i=i+1
!!$     END IF
     i=i+1
  END DO
  CLOSE(2)

  OPEN(2,file='lyaflux_songaila',status='old',action='read')
  DO
     READ(2,*,iostat=iend) z(i),flux(i),fluxmax(i),fluxmin(i)
     IF (iend /= 0) EXIT
     if (z(i) < 15.0) i=i+1
  END DO
  CLOSE(2)

  n=i-1
  PRINT *,n,MINVAL(z(1:n)),MAXVAL(z(1:n))

  CALL pgpt(n,z(1:n),LOG10(-LOG(MAX(flux(1:n),1.e-20))),-3)

  num(:)=0.0
  fluxmean(:)=0.0
  err1max(:)=-100.0
  err1min(:)=100.0
  err2max(:)=-100.0
  err2min(:)=100.0
  DO i=1,n
     DO j=1,nbin
        IF ((z(i) < zhi(j)) .AND. (z(i) >= zlo(j))) THEN
           num(j)=num(j)+1
           fluxmean(j)=fluxmean(j)+flux(i)
           !err(j)=err(j)+flux(i)*flux(i)
           err1max(j)=MAX(err1max(j),flux(i))
           err1min(j)=MIN(err1min(j),flux(i))
           err2max(j)=MAX(err2max(j),fluxmax(i))
           err2min(j)=MIN(err2min(j),fluxmin(i))           
           !err2max(j)=err2max(j)+(fluxmax(i)-flux(i))**2
           !err2min(j)=err2min(j)+(fluxmin(i)-flux(i))**2
        END IF
     END DO
  END DO
  WHERE(num > 0) fluxmean=fluxmean/REAL(num)
  WHERE(num > 0) err2max=(err2max-fluxmean)/SQRT(REAL(num))+fluxmean
  WHERE(num > 0) err2min=(err2min-fluxmean)/SQRT(REAL(num))+fluxmean
  !WHERE(num > 0) err2max=SQRT(err2max/REAL(num))+fluxmean
  !WHERE(num > 0) err2min=-SQRT(err2min/REAL(num))+fluxmean


  OPEN(2,file='lyabinned',status='unknown')
  DO i=1,nbin
     IF (num(i) > 0) WRITE(2,'(5f10.4)') zmean(i),fluxmean(i),err1max(i),err1min(i)
     !IF (num(i) > 0) WRITE(2,'(5f10.4)') zmean(i),fluxmean(i),MAX(err1max(i),err2max(i)),MIN(err1min(i),err2min(i))
  END DO
  CLOSE(2)

  CALL pgsci(2)
  CALL pgpt(nbin,zmean(1:nbin),LOG10(-LOG(MAX(fluxmean(1:nbin),1.e-20))),-4)
  CALL pgerry(nbin,zmean(1:nbin),LOG10(-LOG(MAX(err1max(1:nbin),1.e-20))),LOG10(-LOG(MAX(1.e-20,MIN(err1min(1:nbin),0.9999)))),1.0)
  call pgsci(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL pgenv(2.0,6.5,-1.0,1.5,0,20)

  i=1
  OPEN(2,file='lybflux_fan',status='old',action='read')
  DO
     READ(2,*,iostat=iend) z(i),flux(i),dum
     IF (iend /= 0) EXIT
     fluxmax(i)=flux(i)+dum
     fluxmin(i)=flux(i)-dum
!!$     IF (flux(i) > 0.0) THEN
!!$        flux(i)=-LOG(flux(i))
!!$        i=i+1
!!$     END IF
     i=i+1
  END DO
  CLOSE(2)

  OPEN(2,file='lybflux_songaila',status='old',action='read')
  DO
     READ(2,*,iostat=iend) z(i),flux(i),fluxmax(i),fluxmin(i)
     IF (iend /= 0) EXIT
     if (z(i) < 15.0) i=i+1
  END DO
  CLOSE(2)

  n=i-1
  PRINT *,n,MINVAL(z(1:n)),MAXVAL(z(1:n))

  CALL pgpt(n,z(1:n),LOG10(-LOG(MAX(flux(1:n),1.e-20))),-3)

  OPEN(2,file='lybbinned',status='unknown')
  DO i=1,n
     IF (z(i) <= 5.5) WRITE(2,'(5f10.4)') z(i),flux(i),fluxmax(i),fluxmin(i)
  END DO

  num(:)=0.0
  fluxmean(:)=0.0
  err1max(:)=-100.0
  err1min(:)=100.0
  err2max(:)=-100.0
  err2min(:)=100.0
  DO i=1,n
     DO j=1,nbin
        IF ((z(i) < zhi(j)) .AND. (z(i) >= zlo(j))) THEN
           num(j)=num(j)+1
           fluxmean(j)=fluxmean(j)+flux(i)
           !err(j)=err(j)+flux(i)*flux(i)
           err1max(j)=MAX(err1max(j),flux(i))
           err1min(j)=MIN(err1min(j),flux(i))
           err2max(j)=MAX(err2max(j),fluxmax(i))
           err2min(j)=MIN(err2min(j),fluxmin(i))           
           !err2max(j)=err2max(j)+(fluxmax(i)-flux(i))**2
           !err2min(j)=err2min(j)+(fluxmin(i)-flux(i))**2
        END IF
     END DO
  END DO
  WHERE(num > 0) fluxmean=fluxmean/REAL(num)
  WHERE(num > 0) err2max=(err2max-fluxmean)/SQRT(REAL(num))+fluxmean
  WHERE(num > 0) err2min=(err2min-fluxmean)/SQRT(REAL(num))+fluxmean
  !WHERE(num > 0) err2max=SQRT(err2max/REAL(num))+fluxmean
  !WHERE(num > 0) err2min=-SQRT(err2min/REAL(num))+fluxmean


  
  DO i=1,nbin
     IF ((num(i) > 0) .AND. (zmean(i) > 5.5)) WRITE(2,'(5f10.4)') zmean(i),fluxmean(i),err1max(i),err1min(i)
  END DO
  CLOSE(2)

  CALL pgsci(2)
  CALL pgpt(nbin,zmean(1:nbin),LOG10(-LOG(MAX(fluxmean(1:nbin),1.e-20))),-4)
  DO i=1,nbin
     PRINT *,zmean(i),LOG10(-LOG(MAX(err1max(i),err2max(i),1.e-20))),LOG10(-LOG(MAX(1.e-20,MIN(err1min(i),err2min(i),0.9999))))
  END DO
  CALL pgerry(nbin,zmean(1:nbin),LOG10(-LOG(MAX(err1max(1:nbin),1.e-20))),LOG10(-LOG(MAX(1.e-20,MIN(err1min(1:nbin),0.9999)))),1.0)
  call pgsci(1)

  call pgclos

END PROGRAM binned
