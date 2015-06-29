program main

  use network

  implicit none

  double precision :: Xin(nspec), Xout(nspec)
  double precision :: rho, temp
  double precision :: tmax
  double precision enucdot

  ! initialize the network
  call network_init()

  ! setup the composition
  Xin(:) = 0
  Xin(ih1) = 0.7d0
  Xin(ihe4) = 0.28d0
  Xin(ic12) = 0.02d0
  
  ! burn
  temp = 9.e8
  rho = 2.e6
  tmax = 0.001

  call react(Xin, temp, rho, tmax, Xout, enucdot)

  print *, Xout

end program main
