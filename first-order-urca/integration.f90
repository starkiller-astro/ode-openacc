module integration
  use bl_types 
  implicit none

  contains 

    function gauss_legendre_5(f, fpars) result(igral)
      ! Do 5-pt Gauss Legendre integration
      ! Weight function is w(x) = 1

      integer, parameter :: N=5
      real(kind=dp_t) :: igral
      real(kind=dp_t), external :: f
      real(kind=dp_t), intent(in) :: fpars(:)
      real(kind=dp_t), parameter, dimension(N) :: xk = (/ &
        -0.906179845938664d0, &
        -0.538469310105683d0, &
        0.0d0, &
        0.538469310105683d0, &
        0.906179845938664d0 /)
      real(kind=dp_t), parameter, dimension(N) :: wk = (/ &
        0.236926885056189d0, &
        0.478628670499366d0, &
        0.568888888888889d0, &
        0.478628670499366d0, &
        0.236926885056189d0 /)
      integer :: j

      igral = 0.0d0
      do j=1,N
        igral = igral + wk(j)*f(xk(j), fpars)
      end do
      return
    end function gauss_legendre_5

    function gauss_laguerre_5(f, fpars) result(igral)
      ! Do 5-pt Gauss Laguerre integration
      ! Weight function is w(x) = x^0 * exp(-x)

      integer, parameter :: N=5
      real(kind=dp_t) :: igral
      real(kind=dp_t), external :: f
      real(kind=dp_t), intent(in) :: fpars(:)
      real(kind=dp_t), parameter, dimension(N) :: xk = (/ &
        0.263560319718141d0, &
        0.141340305910652d1, &
        0.359642577104072d1, &
        0.708581000585884d1, &
        0.126408008442758d2 /)
      real(kind=dp_t), parameter, dimension(N) :: wk = (/ &
        0.521755610582809d0, &
        0.398666811083176d0, &
        0.759424496817076d-1, &
        0.361175867992205d-2, &
        0.233699723857762d-4 /)
      integer :: j

      igral = 0.0d0
      do j=1,N
        igral = igral + wk(j)*f(xk(j), fpars)*exp(xk(j))
      end do
      return
    end function gauss_laguerre_5

end module integration
