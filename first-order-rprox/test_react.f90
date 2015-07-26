program main

  use network
  use bl_types

  implicit none

  integer, parameter :: NPTS = 8192
  real(kind=dp_t) :: Xin(nspec, NPTS), Xout(nspec, NPTS)
  real(kind=dp_t) :: rho(NPTS), temp(NPTS)
  real(kind=dp_t) :: tmax, strt, stp
  integer :: ierr

  interface
     subroutine react(Xin, T, rho, tmax, Xout, ierr)
       import dp_t
       real(kind=dp_t), intent(in) :: Xin(:,:), T(:), rho(:), tmax
       real(kind=dp_t), intent(out) :: Xout(:,:)
       integer,         intent(out) :: ierr
     end subroutine react
  end interface

  ! initialize the network
  call network_init()

  ! setup the composition
  Xin(:,:) = 0
  Xin(ih1,:) = 0.7d0
  Xin(ihe4,:) = 0.28d0
  Xin(ic12,:) = 0.02d0

  ! burn
  temp(:) = 9.e8
  rho(:) = 2.e6
  tmax = 0.00001

  call cpu_time(strt)
  call react(Xin, temp, rho, tmax, Xout, ierr)
  call cpu_time(stp)

  print *, Xout(:,1)
  print *, ierr
  print *, 'runtime: ', stp - strt, ' seconds'

end program main
