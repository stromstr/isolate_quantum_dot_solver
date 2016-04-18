subroutine buildspin
use globle
implicit none
integer :: ni, nj, istat, ixyz, iorbs
real*8 :: dtmp1
!
!allocate(dmtmp1(lmat, lmat), STAT=istat)
!allocate(dmtmp2(lmat, lmat), STAT=istat)
!allocate(dmtmp3(lmat, lmat), STAT=istat)
!allocate(dmtmp4(lmat, lmat), STAT=istat)
!allocate(dmtmp5(lmat, lmat), STAT=istat)
allocate(pauli(2,2,3), STAT=istat)
allocate(sopr(lmat,lmat,3,norbs), STAT=istat)
!  x
pauli(1,1,1) = czero                            ! 1/2 * ( 0  1 )
pauli(1,2,1) = dcmplx(0.5d0, 0.d0)              !       ( 1  0 )
pauli(2,1,1) = dcmplx(0.5d0, 0.d0)
pauli(2,2,1) = czero
!  y
pauli(1,1,2) = czero                            ! 1/2 * ( 0 -i )
pauli(1,2,2) = dcmplx(0.d0, -0.5d0)             !       ( i  0 )
pauli(2,1,2) = dcmplx(0.d0,  0.5d0)
pauli(2,2,2) = czero
!  z
pauli(1,1,3) = dcmplx(0.5d0, 0.d0)              ! 1/2 * ( 1  0 )
pauli(1,2,3) = czero                            !       ( 0 -1 ) 
pauli(2,1,3) = czero
pauli(2,2,3) = dcmplx(-0.5d0, 0.d0)
!
do iorbs=1,norbs
   do ixyz=1,3
      cmtmp1(1:lmat,1:lmat) = czero
      do ni=1,2
        do nj=1,2
!             if ( cdabs(pauli(ni,nj, ixyz)) .lt. dpico) cycle
             cmtmp2(1:lmat, 1:lmat) = dcmplx(amsall(1:lmat, 1:lmat, iorbs, ni), 0.d0)
             !call calcams(norbs, nspin, iorbs, nj, lmat, tmp)
             cmtmp3(1:lmat, 1:lmat) = dcmplx(amsall(1:lmat, 1:lmat, iorbs, nj), 0.d0)
             call zgemm('c', 'n', lmat, lmat, lmat, pauli(ni,nj,ixyz), cmtmp2, lmat,  &
                        cmtmp3, lmat, cunity, cmtmp1, lmat)
         end do
      end do
      sopr(1:lmat,1:lmat,ixyz,iorbs) = cmtmp1(1:lmat,1:lmat)
   end do
end do
write(10,*) 'check commutation relation [Sx, Sy] = iSz, [Sy, Sz] = iSx, [Sz, Sx] = iSy'
do iorbs=1,norbs
   cmtmp1(1:lmat,1:lmat) = sopr(1:lmat,1:lmat,1,iorbs)
   cmtmp2(1:lmat,1:lmat) = sopr(1:lmat,1:lmat,2,iorbs)
   cmtmp3(1:lmat,1:lmat) = sopr(1:lmat,1:lmat,3,iorbs)
! [x, y]
   call zgemm('n', 'n', lmat, lmat, lmat,  cunity, cmtmp1, lmat, cmtmp2, lmat,  &
               czero, cmtmp4, lmat)
   call zgemm('n', 'n', lmat, lmat, lmat, -cunity, cmtmp2, lmat, cmtmp1, lmat,  &
              cunity, cmtmp4, lmat)
   cmtmp5(1:lmat,1:lmat) = cmtmp4(1:lmat,1:lmat) - eye * cmtmp3(1:lmat,1:lmat)
   call cmaxmat(lmat, lmat, cmtmp5, lmat, dtmp1)
   if (dtmp1 .ge. dpico) then
      write(10,*)'buildspin3d: error! [Sx, Sy] = iSz violated for iorbs = ', iorbs
      write(10,*)'buildspin3d: largest deviation = ', dtmp1
      stop
   end if
! [y, z]
   call zgemm('n', 'n', lmat, lmat, lmat,  cunity, cmtmp2, lmat, cmtmp3, lmat,  &
               czero, cmtmp4, lmat)
   call zgemm('n', 'n', lmat, lmat, lmat, -cunity, cmtmp3, lmat, cmtmp2, lmat,  &
              cunity, cmtmp4, lmat)
   cmtmp5(1:lmat,1:lmat) = cmtmp4(1:lmat,1:lmat) - eye * cmtmp1(1:lmat,1:lmat)
   call cmaxmat(lmat, lmat, cmtmp5, lmat, dtmp1)
   if (dtmp1 .ge. dpico) then
      write(10,*)'buildspin3d: error! [Sy, Sz] = iSx violated for iorbs = ', iorbs
      write(10,*)'buildspin3d: largest deviation = ', dtmp1
      stop
   end if
! [z, x]
   call zgemm('n', 'n', lmat, lmat, lmat,  cunity, cmtmp3, lmat, cmtmp1, lmat,  &
               czero, cmtmp4, lmat)
   call zgemm('n', 'n', lmat, lmat, lmat, -cunity, cmtmp1, lmat, cmtmp3, lmat,  &
              cunity, cmtmp4, lmat)
   cmtmp5(1:lmat,1:lmat) = cmtmp4(1:lmat,1:lmat) - eye * cmtmp2(1:lmat,1:lmat)
   call cmaxmat(lmat, lmat, cmtmp5, lmat, dtmp1)
   if (dtmp1 .ge. dpico) then
      write(10,*)'buildspin3d: error! [Sz, Sx] = iSy violated for iorbs = ', iorbs
      write(10,*)'buildspin3d: largest deviation = ', dtmp1
      stop
   end if
   write(10,*) 'commutation checked for orbs:', iorbs
end do
end subroutine
!
