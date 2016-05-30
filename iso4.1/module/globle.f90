module globle
implicit none
complex*16, save, allocatable :: cmtmp1(:,:), cmtmp2(:,:), cmtmp3(:,:), cmtmp4(:,:), cmtmp5(:,:), cmtmp6(:,:), cmtmp7(:,:), cmtmp8(:,:)
complex*16, save, allocatable :: cmtmp9(:,:), cmtmp10(:,:), cmtmp11(:,:), cmtmp12(:,:), cmtmp13(:,:), cmtmp14(:,:), cmtmp15(:,:)
real*8, save, allocatable :: dmtmp1(:,:), dmtmp2(:,:), dmtmp3(:,:), dmtmp4(:,:), dmtmp5(:,:)
complex*16, save, allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:), tmp4(:,:), tmp5(:,:)
complex*16, allocatable :: hs(:,:), hss(:,:), w(:), vl(:,:), vr(:,:), work(:), rwork(:)
complex*16, save, allocatable :: sopr(:,:,:,:)
complex*16, save, allocatable :: pauli(:,:,:)            ! Pauli matrices
complex*16, save :: czero = (0.d0, 0.d0)
complex*16, save :: cunity = (1.d0, 0.d0)
complex*16, save :: eye = (0.d0, 1.d0)
complex*16, save :: ccompl = (0.d0, 1.d0)
real*8, save :: hbar = 0.65821
!real*8, save :: dpico = 1.d-12
real*8, parameter :: dpico = 1.d-12, dnano = 1.d-9
real*8, save, allocatable :: tmp(:,:), ams(:,:)
real*8, save, allocatable :: amsall(:,:,:,:)
real*8 :: temp, e1up, e1down, e2up, e2down, e3up, e3down, uu1, uu2, uu3, t12, t13, t23, u12, u13, u23, j12, j23, rho
integer, save :: nspin, lwmax, norbs, lmat
end module

