subroutine react(Xin, T, rho, tmax, Xout, ierr)

  use network
  use rpar_indices
  use bl_types
  use bl_constants_module

  implicit none
  !$acc routine(dgesv) seq

  ! do a simple first-order backward Euler method for integrating
  ! a stiff-system of ODEs.

  ! Note: this should not been seen as a science code, but rather
  ! just used to test ideas about putting reactions on GPUs

  ! this version works on a single-zone only

  interface
     subroutine rhs(n, t, y, ydot, rpar, ipar)
       !$acc routine seq
       import dp_t, n_rpar_comps
       integer,         intent(in   ) :: n, ipar
       real(kind=dp_t), intent(in   ) :: y(n), t
       real(kind=dp_t), intent(  out) :: ydot(n)
       real(kind=dp_t), intent(inout) :: rpar(n_rpar_comps)
     end subroutine rhs
  end interface
  real(kind=dp_t), intent(in)  :: Xin(:,:), T(:), rho(:), tmax
  real(kind=dp_t), intent(out) :: Xout(:,:)
  integer,         intent(out) :: ierr

  real(kind=dp_t) :: time, dt, I

  integer :: m, n

  real(kind=dp_t) :: dXdt(nspec)
  real(kind=dp_t) :: X1(nspec), X2(nspec), dX1dt(nspec), dX2dt(nspec)
  real(kind=dp_t) :: X_n(nspec), X_np1(nspec), dX(nspec), dfdX(nspec)
  real(kind=dp_t) :: A(nspec, nspec), J(nspec, nspec), b(nspec)
  integer :: info

  real(kind=dp_t) :: eps = 1.d-8

  real(kind=dp_t) :: rpar(n_rpar_comps)
  integer :: ipar

  integer :: ipiv(nspec)

  real(kind=dp_t) :: tol = 1.d-4
  integer :: max_iter = 20
  integer :: iter, npts, p
  logical :: converged

  npts = size(T)
  ierr = 0

  !$acc data copyin(Xin, T, rho, tmax, eps, tol, max_iter)        &
  !$acc copyout(Xout)                                             &
  !$acc copy(ierr)

  !$acc kernels

  ierr = 0
  !$acc loop gang vector reduction(+:ierr)                            &
  !$acc private(m,n,iter,I,time,dt,converged)                         &
  !$acc private(dXdt, X1, X2, dX1dt, dX2dt, X_n, X_np1, dX, A, J, b)  &
  !$acc private(info, ipiv, rpar, ipar)

  !$omp parallel do private(p,m,n,iter,rpar,X_n,dt,converged,time) &
  !$omp private(X_np1,dXdt,ipar,X1,X2,dX1dt,dX2dt,J,I,A,b,dX,info,ipiv) &
  !$omp reduction(+:ierr)
  do p = 1, npts
     ! get an estimate of the timestep by looking at the RHS
     ! dt = min{X/(dX/dt)}
     rpar(:) = ZERO
     rpar(irp_dens) = rho(p)
     rpar(irp_temp) = T(p)
     X_n(:) = Xin(:,p)

     call rhs(nspec, ZERO, X_n, dXdt, rpar, ipar)

     dt = 1.d33
     dt = 0.01*min(dt, abs(Xin(ih1,p)/dXdt(ih1)))

     ! do the integration
     time = ZERO
     do while (time < tmax)
        converged = .false.

        ! initial guess
        X_np1(:) = X_n(:)

        do iter = 1, max_iter

           ! construct a numerical Jacobian via simple differencing
           ! Lapack matrix ordering is described here:
           ! http://www.netlib.org/lapack/lug/node116.html#subsecarrayargs
           ! in general, A(i,j) holds matrix element a_{i,j}

           ! J_{m,n} = df_m/dX_n
           do n = 1, nspec
              X1(:) = X_np1(:)
              X1(n) = X1(n) + eps

              X2(:) = X_np1(:)
              X2(n) = X2(n) - eps

              call rhs(nspec, ZERO, X1, dX1dt, rpar, ipar)
              call rhs(nspec, ZERO, X2, dX2dt, rpar, ipar)

              do m = 1, nspec
                 J(m,n) = 0.5d0*(dX2dt(m) - dX1dt(m))/(X2(n) - X1(n))
              enddo
           enddo

           ! construct the matrix for the linear system
           ! (I - dt J) dX^{n+1} = rhs
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

           ! construct the RHS
           call rhs(nspec, ZERO, X_np1, dXdt, rpar, ipar)
           b(:) = X_n(:) - X_np1(:) + dt*dXdt(:)

           ! solve the linear system -- we'll just use Lapack.  Note that the
           ! input b is overwritten on output
           call dgesv(nspec, 1, A, nspec, ipiv, b, nspec, info)

           dX(:) = b(:)

           X_np1(:) = X_np1(:) + dX(:)

           ! compare to see if we converged
           if (abs(dX(ih1)) < tol*X_n(ih1)) then
              converged = .true.
              exit
           endif

        enddo   ! dX iteration loop

        if (.not. converged) then
           ierr = ierr + 1
           exit
        endif

        if (time + dt > tmax) then
           dt = tmax - time
        endif

        X_n(:) = X_np1(:)
        time = time + dt
     enddo

     Xout(:,p) = X_n(:)
  enddo
  !$acc end kernels
  !$acc end data

  !$omp end parallel do
end subroutine react
