subroutine calcoccup(occtmp, occout)
use globle
implicit none
!
complex*16 :: occtmp(lmat, lmat)
complex*16 :: cc1, cc2
real*8  :: occout
integer :: n, m
cc1 = (0.d0, 0.d0)
cc2 = (0.d0, 0.d0)
!write(16,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!do n=1, lmat
!  do m=1, lmat
!    write(16,*) occtmp(n,m)
!  end do
!end do
call zgemm('c', 'n', lmat, lmat, lmat, cunity, vr, lmat, occtmp, lmat, czero, cmtmp4, lmat)
call zgemm('n', 'n', lmat, lmat, lmat, cunity, cmtmp4, lmat, vr, lmat, czero, cmtmp5, lmat)
!write(16,*)'~~~~~~~~~~~~~~~~~~~~~~~~~'
!do n=1, lmat
!  write(16,*) cmtmp5(n,n)
!end do
!write(16,*) '~~~~~~~~~~~~~~~~~~~~~~'
do n=1, lmat
  cc1 = cc1 + cmtmp5(n,n) * exp(-w(n)/temp)
end do
!
do n=1, lmat
  cc2 = cc2 + exp(-w(n)/temp)
end do
occout = real(cc1) / real(cc2)
end subroutine calcoccup

