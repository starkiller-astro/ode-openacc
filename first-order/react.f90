subroutine react(Xin, T, rho, tmax, Xout, enucdot)

  use network
  use rpar_indices

  implicit none

  ! do a simple first-order backward Euler method for integrating
  ! a stiff-system of ODEs.

  ! Note: this should not been seen as a science code, but rather
  ! just used to test ideas about putting reactions on GPUs

  ! this version works on a single-zone only

  double precision, dimension(nspec), intent(in) :: Xin
  double precision, dimension(nspec), intent(out) :: Xout
  double precision, intent(in) :: T, rho, tmax
  double precision, intent(out) :: enucdot

  double precision :: time, dt, I

  integer :: m, n
    
  double precision :: dXdt(nspec)
  double precision :: X1(nspec), X2(nspec), dX1dt(nspec), dX2dt(nspec)
  double precision :: A(nspec, nspec), J(nspec, nspec), b(nspec)
  double precision :: dfdX
  integer :: info
  
  double precision, parameter :: eps = 1.d-8

  double precision :: rpar(n_rpar_comps)
  integer :: ipar

  integer :: ipiv(nspec)

  integer :: ic12, io16

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")

  ! get an estimate of the timestep by looking at the RHS
  ! dt = min{X/(dX/dt)}
  rpar(irp_dens) = rho
  rpar(irp_temp) = T
  rpar(irp_o16) = Xin(io16)

  call rhs(nspec, 0.0d0, Xin, dXdt, rpar, ipar)
  print *, dXdt

  dt = 1.d33
  dt = min(dt, abs(Xin(ic12)/dXdt(ic12)))

  print *, "dt = ", dt

  Xout(:) = xin(:)
  
  ! do the integration
  time = 0.0d0
  do while (time < tmax)

     ! construct a numerical Jacobian via simple differencing
     ! Lapack matrix ordering is described here:
     ! http://www.netlib.org/lapack/lug/node116.html#subsecarrayargs
     ! in general, A(i,j) holds matrix element a_{i,j}

     ! J_{m,n} = df_m/dX_n
     do n = 1, nspec
        X1(:) = Xout(:)
        X1(n) = X1(n) + eps

        X2(:) = Xout(:)
        X2(n) = X2(n) - eps

        call rhs(nspec, 0.0d0, X1, dX1dt, rpar, ipar)
        call rhs(nspec, 0.0d0, X2, dX2dt, rpar, ipar)
        
        do m = 1, nspec
           dfdX = 0.5d0*(dX2dt(m) - dX1dt(m))/(X2(n) - X1(n))
           J(m,n) = dfdX
        enddo
     enddo

     ! construct the matrix for the linear system
     ! (I - dt A) X^{n+1} = X^n
     do n = 1, nspec
        do m = 1, nspec
           if (m == n) then
              I = 1.0d0
           else
              I = 0.0d0
           endif
           
           A(m,n) = I - dt*J(m,n)
        enddo

        b(n) = xout(n)
     enddo
     
     ! solve the linear system -- we'll just use Lapack.  Note that the
     ! input X is overwritten on output
     call dgesv(nspec, 1, A, nspec, ipiv, b, nspec, info)

     do n = 1, nspec
        xout(n) = b(n)
     enddo
          
     if (time + dt > tmax) then
        dt = tmax - time
     endif

     time = time + dt
  enddo


end subroutine react
