program main

  use network

  implicit none

  double precision :: Xin(nspec), Xout(nspec)
  double precision :: rho, temp
  double precision :: tmax
  double precision enucdot

  integer :: ic12, io16

  ! initialize the network
  !call network_init()
  call network_initialize()

  ! setup the composition
  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  Xin(:) = 0
  Xin(ic12) = 0.5d0
  Xin(io16) = 0.5d0

  ! burn
  temp = 8.e8
  rho = 2.6e9
  tmax = 1000.0

  call react(Xin, temp, rho, tmax, Xout, enucdot)

  print *, Xout

end program main
