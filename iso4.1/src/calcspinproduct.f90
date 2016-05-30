subroutine calcspinproduct(iorbs, iorbs2, dout)
use globle
implicit none
integer :: iorbs, iorbs2, n, m
real*8  :: dout
complex*16 :: rr1, rr2
!
dout = 0.d0
rr1  = (0.d0, 0.d0)
rr2  = (0.d0, 0.d0)
call zgemm('n', 'n', lmat, lmat, lmat, cunity, sopr(1:lmat,1:lmat,1,iorbs), lmat,  &
            sopr(1:lmat,1:lmat,1,iorbs2), lmat, czero,  cmtmp1, lmat)
call zgemm('n', 'n', lmat, lmat, lmat, cunity, sopr(1:lmat,1:lmat,2,iorbs), lmat,  &
            sopr(1:lmat,1:lmat,2,iorbs2), lmat, cunity, cmtmp1, lmat)
call zgemm('n', 'n', lmat, lmat, lmat, cunity, sopr(1:lmat,1:lmat,3,iorbs), lmat,  &
            sopr(1:lmat,1:lmat,3,iorbs2), lmat, cunity, cmtmp1, lmat)
!
call zgemm('c', 'n', lmat, lmat, lmat, cunity, vr, lmat, cmtmp1, lmat, czero, cmtmp2, lmat)
call zgemm('n', 'n', lmat, lmat, lmat, cunity, cmtmp2, lmat, vr, lmat, czero, cmtmp3, lmat)
write(16,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
do n=1, lmat
  write(16,*) cmtmp3(n,n)
end do
do n=1, lmat
  rr1 = rr1 + cmtmp3(n,n) * exp(-w(n)/temp)
end do
!
do n=1, lmat
  rr2 = rr2 + exp(-w(n)/temp)
end do
dout =  real(rr1) / real(rr2)
end subroutine calcspinproduct

