subroutine buildoperator
use globle
implicit none
!include '../include/sizes'
!include '../include/common'
!
integer :: nunity, istat
integer :: ni, nj, nk, ispin
integer :: irowtmp, icoltmp, isgntmp
integer, allocatable :: irow(:,:), icol(:,:), isgn(:,:)
real*8 :: cpu1, cpu2
!
!write(6,*)
!write(6,*)' entering buildoperator '
!call cpu_time(cpu1)
!
call calcams(norbs, nspin, 1, 1, lmat, ams)
amsall(1:lmat, 1:lmat, 1, 1) = ams(1:lmat, 1:lmat)
!
!  
 nunity = 2 * 4**(norbs - 1) 
 allocate(irow(nunity,nspin), icol(nunity,nspin), isgn(nunity,nspin), STAT=istat)
 nk = 0
 do ni=1,lmat 
  do nj=1,lmat
   if (dabs(ams(ni,nj) - 1.d0) .le. dnano) then
    nk = nk + 1
    irow(nk,1) = ni
    icol(nk,1) = nj
    isgn(nk,1) = int(ams(ni,nj))
   end if
  end do
 end do
 if (nk .ne. nunity) then
  write(6,*)
  write(6,*)' error! wrong c1u '
!  call amatout(lmat, lmat, ams)
  stop
 end if
!
 dmtmp1(1:lmat, 1:lmat) = 0.d0
 do nj=1,nunity
  call swapspin(lmat, norbs, irow(nj,1), icol(nj,1), 1, irowtmp, &
                icoltmp, isgntmp)
  irow(nj,2) = irowtmp
  icol(nj,2) = icoltmp
  isgn(nj,2) = isgntmp * isgn(nj,1)
  dmtmp1(irowtmp, icoltmp) = dble(isgn(nj,2))
 end do
 amsall(1:lmat, 1:lmat, 1, 2) = dmtmp1(1:lmat, 1:lmat)
!
 ispin=2
 do ni=2,norbs
   dmtmp1(1:lmat, 1:lmat) = 0.d0
   do nj=1,nunity
     call swapms(lmat, norbs, nspin, irow(nj,ispin), icol(nj,ispin), ni,         &
                 ispin, irowtmp, icoltmp, isgntmp)
     dmtmp1(irowtmp, icoltmp) = dble(isgn(nj,ispin)) * dble(isgntmp)
     irow(nj,ispin) = irowtmp
     icol(nj,ispin) = icoltmp 
     isgn(nj,ispin) = int(dmtmp1(irowtmp, icoltmp))
   end do
   amsall(1:lmat, 1:lmat, ni, ispin) = dmtmp1(1:lmat, 1:lmat) 
!
   ispin = 3 - ispin
!
   dmtmp1(1:lmat, 1:lmat) = 0.d0
   do nj=1,nunity 
     call swapspin(lmat, norbs, irow(nj,3-ispin), icol(nj,3-ispin), ni, irowtmp, &
                   icoltmp, isgntmp)
     dmtmp1(irowtmp, icoltmp) = dble(isgn(nj,3-ispin)) * dble(isgntmp)
     irow(nj,ispin) = irowtmp
     icol(nj,ispin) = icoltmp
     isgn(nj,ispin) = int(dmtmp1(irowtmp, icoltmp))
   end do
   amsall(1:lmat, 1:lmat, ni, ispin) = dmtmp1(1:lmat, 1:lmat) 
 end do

 deallocate(irow, icol, isgn, STAT=istat)
end subroutine buildoperator
