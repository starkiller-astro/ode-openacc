subroutine react(nspec, Xin, T, rho, tmax, Xout, enucdot)

  implicit none

  ! do a simple first-order backward Euler method for integrating
  ! a stiff-system of ODEs.

  ! Note: this should not been seen as a science code, but rather
  ! just used to test ideas about putting reactions on GPUs

  ! this version works on a single-zone only

  integer :: nspec
  double precision, dimension(nspec), intent(in) :: Xin
  double precision, dimension(nspeC), intent(out) :: Xout
  double precision, intent(in) :: T, rho, tmax
  double precision, intent(out) :: enucdot

  double precision :: J(nspec, nspec)

  double precision, parameter :: eps = 1.d-8
  

  ! get an estimate of the timestep by looking at the RHS
  ! dt = min{X/(dX/dt)}
  call rhs(rho, T, Xin, dXdt)

  dt = 1.d33
  do n = 1, nspec
     dt = min(dt, Xin(n)/dXdt(n))
  enddo

  ! do the integration
  t = 0.0d0
  do while (t < tmax)

     ! construct a numerical Jacobian via simple differencing
     ! Lapack matrix ordering is described here:
     ! http://www.netlib.org/lapack/lug/node116.html#subsecarrayargs
     ! in general, A(i,j) holds matrix element a_{i,j}

     ! J_{m,n} = df_m/dX_n
     do n = 1, nspec
        X1(:) = Xin(:)
        X1(n) = X1(n) + eps

        X2(:) = Xin(:)
        X2(n) = Xin(n) - eps

        call rhs(rho, T, X1, dX1dt)
        call rhs(rho, T, X2, dX2dt)
        
        do m = 1, nspec
           dfdX = (dX2dt(m) - dX1dt(m))/(X2(n) - X1(n))
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
     enddo
     
     ! solve the linear system -- we'll just use Lapack.  Note that the
     ! input X is overwritten on output


     
     if (t + dt > tmax) then
        dt = tmax - t
     endif

     t = t + dt
  enddo


end subroutine react
