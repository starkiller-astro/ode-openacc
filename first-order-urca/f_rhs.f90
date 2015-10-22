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
  !!$acc routine(screenz) seq ! For now, I don't use screenz

  ! our convention is that y(1:nspec) are the species (in the same
  ! order as defined in network.f90, and y(nspec+1) is the temperature
  integer,         intent(in   ) :: n, ipar
  real(kind=dp_t), intent(in   ) :: y(n), t
  real(kind=dp_t), intent(  out) :: ydot(n)
  real(kind=dp_t), intent(inout) :: rpar(:)

!  real(kind=dp_t) :: ymass(nspec)

  !real(kind=dp_t) :: dens
  !real(kind=dp_t) :: temp, T9, T9a, dT9dt, dT9adt

  real(kind=dp_t) ::  mu_elec
  real(kind=dp_t) ::  capture_rate, emission_rate

  real(kind=dp_t), PARAMETER :: &
                     one_twelvth = 1.0d0/12.0d0, &
                     five_sixths = 5.0d0/ 6.0d0, &
                       one_third = 1.0d0/ 3.0d0, &
                      two_thirds = 2.0d0/ 3.0d0

  interface 
    function rate_capture_na23(rpar) result(lambda)
      !$acc routine seq
      use bl_types
      real(kind=dp_t)                :: lambda
      real(kind=dp_t), dimension(*), intent(in) :: rpar
    end function rate_capture_na23
  end interface

  interface 
    function rate_emission_ne23(rpar) result(lambda)
      !$acc routine seq
      use bl_types
      real(kind=dp_t)                :: lambda
      real(kind=dp_t), dimension(*), intent(in) :: rpar
    end function rate_emission_ne23
  end interface

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

  ! mu_elec = 1.0/(sum_i z_i * x_i * y_i / a_i)
  ! y_i: ionization fraction (assume y_i = 1 below)
  mu_elec = 1.0d0/(zion(ine23_)*y(ine23_)/aion(ine23_) + &
                  zion(ina23_)*y(ina23_)/aion(ina23_))
  rpar(irp_mu_elec) = mu_elec

  capture_rate  = rate_capture_na23(rpar)
  emission_rate = rate_emission_ne23(rpar)

  ydot(ine23_) = -emission_rate*y(ine23_) + capture_rate*y(ina23_) 
  ydot(ina23_) = -capture_rate*y(ina23_) + emission_rate*y(ine23_) 

  return

end subroutine rhs

function rate_emission_ne23(rpar) result (lambda)
  !$acc routine seq
  ! Compute the beta emission rate of Ne-23 given density, temp, and electron
  ! fermi energy. For now, assumes gs->gs transition.
  use bl_constants_module
  use physical_constants_module
  use rpar_indices 
  use network_indices
  use phase_par_indices

  implicit none

  real(kind=dp_t)                :: lambda
  real(kind=dp_t), dimension(*), intent(in)    :: rpar
  real(kind=dp_t), dimension(n_phase_par_inds)  :: fpar
  real(kind=dp_t), parameter     :: ft = 1.0d0 ! placeholder value
  real(kind=dp_t)                :: ufermi, kt
 
  interface 
    function gauss_legendre_5pt_emission(fpar) result(igral)
      !$acc routine seq
      import dp_t
      real(kind=dp_t) :: igral
      real(kind=dp_t), dimension(*), intent(in) :: fpar
    end function gauss_legendre_5pt_emission
  end interface

  ! assuming gs->gs so E_i = E_j = 0
  fpar(iqn)   = (MASS_NE23_CONST-MASS_NA23_CONST)*MEV_PER_AMU_CONST/MELEC_CONST
  fpar(itemp) = rpar(irp_temp)
  
  ! placeholder approx eta for T->0, need to get it from EOS (Eq. 4c, FFN 1980)
  ufermi = 0.5d-2*(rpar(irp_dens)/rpar(irp_mu_elec))**THIRD
  kt     = KBOLT_CONST*rpar(irp_temp)
  fpar(ieta)  = ufermi/kt
  lambda = log(2.0d0)*gauss_legendre_5pt_emission(fpar)/ft

end function rate_emission_ne23

function phase_emission_ne23(v, fpar) result (phi)
  !$acc routine seq
  ! w is the mass+kinetic energy in units of mc^2
  ! v is related to w via a linear transformation for emission phase
  ! v is the integration variable to scale the limits to [-1,1]
  ! w = 1 + (v + 1)*(qn - 1)/2
  use phase_par_indices
  use network_indices
  use physical_constants_module

  implicit none

  real(kind=dp_t), intent(in) :: v
  real(kind=dp_t), dimension(*), intent(in) :: fpar
  real(kind=dp_t) ::  w, g, s, z
  real(kind=dp_t) ::  phi
 
  interface
    function gzw(inuc,w) result(g)
      !$acc routine seq
      import dp_t
      real(kind=dp_t)               :: g
      real(kind=dp_t), intent(in)   :: w
      integer, intent(in)           :: inuc
    end function gzw
  end interface

  w = 1.0d0 + (v + 1.0d0)*(fpar(iqn) - 1.0d0)/2.0d0
  g = gzw(ina23_,w)
  z = MELEC_CONST/(KBOLT_CONST*fpar(itemp))
  s = 1.0d0/(1.0d0 + exp((w-1.0d0)*z-fpar(ieta)))
  phi = w**2*(fpar(iqn) - w)**2*g*(1.0d0-s)

end function phase_emission_ne23

function gauss_legendre_5pt_emission(fpar) result(igral)
  !$acc routine seq
  ! Do 5-pt Gauss Legendre integration
  ! Weight function is w(x) = 1
  use bl_types
  use integration

  real(kind=dp_t) :: igral, xkj, wkj
  real(kind=dp_t), dimension(*), intent(in) :: fpar
  integer :: j
  integer, parameter  ::  N=5

  interface 
    function phase_emission_ne23(v, fpar) result(phi)
      !$acc routine seq
      import dp_t
      real(kind=dp_t), intent(in) :: v
      real(kind=dp_t) ::  phi
      real(kind=dp_t), dimension(*), intent(in) :: fpar
    end function phase_emission_ne23
  end interface

  igral = 0.0d0
  !$acc data present(gauss_legendre_5pt_xk, gauss_legendre_5pt_wk)
  do j=1,N
    wkj = gauss_legendre_5pt_wk(j)
    xkj = gauss_legendre_5pt_xk(j)
    igral = igral + wkj*phase_emission_ne23(xkj, fpar)
  end do
  !$acc end data
  return
end function gauss_legendre_5pt_emission

function rate_capture_na23(rpar) result (lambda)
  !$acc routine seq
  ! Compute the beta capture rate of Na-23 given density, temp, and electron
  ! fermi energy. For now, assumes gs->gs transition.
  use physical_constants_module
  use bl_constants_module
  use rpar_indices 
  use network_indices
  use phase_par_indices

  implicit none

  real(kind=dp_t)                :: lambda
  real(kind=dp_t), dimension(*), intent(in)    :: rpar
  real(kind=dp_t), dimension(n_phase_par_inds)  :: fpar
  real(kind=dp_t), parameter     :: ft = 1.0d0 ! placeholder value
  real(kind=dp_t)                :: ufermi, kt
 
  interface 
    function gauss_laguerre_5pt_capture(fpar) result(igral)
      !$acc routine seq
      import dp_t
      real(kind=dp_t) :: igral
      real(kind=dp_t), dimension(*), intent(in) :: fpar
    end function gauss_laguerre_5pt_capture
  end interface

  ! assuming gs->gs so E_i = E_j = 0
  fpar(iqn)   = (MASS_NA23_CONST-MASS_NE23_CONST)*MEV_PER_AMU_CONST/MELEC_CONST
  
  fpar(itemp) = rpar(irp_temp)
  
  ! placeholder approx eta for T->0, need to get it from EOS (Eq. 4c, FFN 1980)
  ufermi = 0.5d-2*(rpar(irp_dens)/rpar(irp_mu_elec))**THIRD
  kt     = KBOLT_CONST*rpar(irp_temp)
  fpar(ieta)  = ufermi/kt
  lambda = log(2.0d0)*gauss_laguerre_5pt_capture(fpar)/ft

end function rate_capture_na23

function phase_capture_na23(v, fpar) result (phi)
  !$acc routine seq
  ! w is the mass+kinetic energy in units of mc^2
  ! v is related to w via a linear transformation for capture phase
  ! v is the integration variable to scale the limits to [0,Infinity]
  ! w = v + wl, wl = 1 (for qn > -1), wl = abs(qn) (for qn < -1)
  use bl_types
  use phase_par_indices
  use network_indices
  use physical_constants_module

  implicit none

  real(kind=dp_t), intent(in) ::  v 
  real(kind=dp_t), dimension(*), intent(in) :: fpar
  real(kind=dp_t) ::  w, wl, g, s, z
  real(kind=dp_t) ::  phi

  interface
    function gzw(inuc,w) result(g)
      !$acc routine seq
      import dp_t
      real(kind=dp_t)               :: g
      real(kind=dp_t), intent(in)   :: w
      integer, intent(in)           :: inuc
    end function gzw
  end interface

  if (fpar(iqn) <= -1.0d0) then
    wl = abs(fpar(iqn))
  else
    wl = 1.0d0
  end if
  w = v + wl
  g = gzw(ina23_,w)
  z = MELEC_CONST/(KBOLT_CONST*fpar(itemp))
  s = 1.0d0/(1.0d0 + exp((w-1.0d0)*z-fpar(ieta)))
  phi = w**2*(fpar(iqn) + w)**2*g*s

end function phase_capture_na23

function gauss_laguerre_5pt_capture(fpar) result(igral)
  !$acc routine seq
  ! Do 5-pt Gauss Laguerre integration
  ! Weight function is w(x) = x^0 * exp(-x)
  use bl_types
  use integration
  
  integer, parameter :: N=5
  real(kind=dp_t) :: igral, wkj, xkj
  real(kind=dp_t), dimension(*), intent(in) :: fpar
  integer :: j

  interface 
    function phase_capture_na23(v, fpar) result(phi)
      !$acc routine seq
      import dp_t
      real(kind=dp_t), intent(in) :: v
      real(kind=dp_t) ::  phi
      real(kind=dp_t), dimension(*), intent(in) :: fpar
    end function phase_capture_na23
  end interface

  igral = 0.0d0
  !$acc data present(gauss_laguerre_5pt_xk, gauss_laguerre_5pt_wk)
  do j=1,N
    wkj = gauss_laguerre_5pt_wk(j)
    xkj = gauss_laguerre_5pt_xk(j)
    igral = igral + wkj*phase_capture_na23(xkj, fpar)*exp(xkj)
  end do
  !$acc end data
  return
end function gauss_laguerre_5pt_capture

function gzw(inuc,w) result(g)
  !$acc routine seq
  ! Calculate the G(Z,w) factor (Eq. 5a-5b of FFN 1980)
  ! Signs are chosen for electron emission and capture
  use network
  use bl_constants_module
  use physical_constants_module

  implicit none

  real(kind=dp_t)               :: g
  real(kind=dp_t), intent(in)   :: w
  integer, intent(in)           :: inuc
  real(kind=dp_t)               :: p, s, x, z, r, a

  ! Parameters for RMS ground state nuclear charge radii
  ! From ADNDT 2013, Eq. 8
  real(kind=dp_t), parameter    ::  r0 = 0.9071 ! fm
  real(kind=dp_t), parameter    ::  r1 = 1.105  ! fm
  real(kind=dp_t), parameter    ::  r2 = -0.548 ! fm
  ! r for FFN 1980 is in electron compton wavelengths
  ! From PDG 2014
  real(kind=dp_t), parameter    ::  lambda_elec_div_2pi = 3.8615926800d2  ! fm
  
  p = sqrt(w**2 - 1.0d0)
  z = zion(inuc)
  s = sqrt(1.0d0 - (ALPHA_CONST*z)**2) ! I'll need this for the general case
  x = ALPHA_CONST*z*w/p
  a = aion(inuc)

  !print *,'a: ',a
  !r = (2.908d-3)*a**(THIRD) - (2.437d0)*a**(-THIRD) ! From FFN 1980, this gives a negative r
  !print *,'r: ',r

  ! Calculate rms nuclear charge radius in units of electron compton wavelength
  r = ((r0 + r1/a**(TWO3RD) + r2/a**(FOUR3RD))*a**(THIRD))/(2.0d0*M_PI*lambda_elec_div_2pi)
  !print *,'r: ',r

  
  ! For now, use the non-relativistic approximation for gzw (doesn't use s)
  g = exp(-2.0d0*M_PI*abs(x))*2.0d0*M_PI*ALPHA_CONST*z*(2.0d0*ALPHA_CONST*z*r)**(-(ALPHA_CONST*z)**2)
  !print *,'g: ',g
end function gzw
