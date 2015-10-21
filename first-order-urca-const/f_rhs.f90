! For Example:
! The change in number density of C12 is
! d(n12)/dt = - 2 * 1/2 (n12)**2 <sigma v>
!
! where <sigma v> is the average of the relative velocity times the cross
! section for the reaction, and the factor accounting for the total number
! of particle pairs has a 1/2 because we are considering a reaction involving 
! identical particles (see Clayton p. 293).  Finally, the -2 means that for
! each reaction, we lose 2 carbon nuclei.
!
! The corresponding Mg24 change is
! d(n24)/dt = + 1/2 (n12)**2 <sigma v>
!
! note that no factor of 2 appears here, because we create only 1 Mg nuclei.
!
! Switching over to mass fractions, using n = rho X N_A/A, where N_A is
! Avagadro's number, and A is the mass number of the nucleon, we get
!
! d(X12)/dt = -2 *1/2 (X12)**2 rho N_A <sigma v> / A12
!
! d(X24)/dt = + 1/2 (X12)**2 rho N_A <sigma v> (A24/A12**2)
!
! these are equal and opposite.
!
! The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.

! we will always refer to the species by integer indices that come from
! the network module -- this makes things robust to a shuffling of the 
! species ordering

subroutine rhs(n, t, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use network_indices
  use rpar_indices

  implicit none
  !$acc routine seq
  !$acc routine(screenz) seq

  ! our convention is that y(1:nspec) are the species (in the same
  ! order as defined in network.f90, and y(nspec+1) is the temperature
  integer,         intent(in   ) :: n, ipar
  real(kind=dp_t), intent(in   ) :: y(n), t
  real(kind=dp_t), intent(  out) :: ydot(n)
  real(kind=dp_t), intent(inout) :: rpar(:)

  integer :: k
!  real(kind=dp_t) :: ymass(nspec)

  !real(kind=dp_t) :: dens
  !real(kind=dp_t) :: temp, T9, T9a, dT9dt, dT9adt

  real(kind=dp_t) :: capture_rate, decay_rate

  real(kind=dp_t), PARAMETER :: &
                     one_twelvth = 1.0d0/12.0d0, &
                     five_sixths = 5.0d0/ 6.0d0, &
                       one_third = 1.0d0/ 3.0d0, &
                      two_thirds = 2.0d0/ 3.0d0

  !dens = rpar(irp_dens)
  !temp = rpar(irp_temp)

  ! compute the molar fractions -- needed for the screening
  !ymass(ic12_) = y(1)/aion(ic12_)
  !ymass(io16_) = X_O16/aion(io16_)
  !ymass(img24_) = (ONE - y(1) - X_O16)/aion(img24_)


  ! compute some often used temperature constants
  !T9     = temp/1.e9_dp_t
  !dT9dt  = ONE/1.e9_dp_t
  !T9a    = T9/(1.0e0_dp_t + 0.0396e0_dp_t*T9)
  !dT9adt = (T9a / T9 - (T9a / (1.0e0_dp_t + 0.0396e0_dp_t*T9)) * 0.0396e0_dp_t) * dT9dt

  capture_rate = 0.1_dp_t
  decay_rate   = 0.2_dp_t

  ydot(ine23_) = -decay_rate*y(ine23_) + capture_rate*y(ina23_) 
  ydot(ina23_) = -capture_rate*y(ina23_) + decay_rate*y(ine23_) 

  return

end subroutine rhs


