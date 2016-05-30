subroutine cmaxmat(dim1, dim2, cmat, ldm, dmax)
implicit none
!
integer,    intent(in)  :: dim1, dim2, ldm
complex*16, intent(in)  :: cmat(ldm,*)
real*8,     intent(out) :: dmax
!
integer                 :: ni,nj
!
dmax = 0.d0
do nj=1,dim2
 do ni=1,dim1
  dmax = max(dmax, cdabs(cmat(ni,nj)))
 enddo
enddo
!
end subroutine cmaxmat

