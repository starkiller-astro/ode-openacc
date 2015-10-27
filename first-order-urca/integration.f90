module integration

  use bl_types 
  implicit none

  real(kind=dp_t), allocatable :: gauss_legendre_5pt_xk(:)
  real(kind=dp_t), allocatable :: gauss_legendre_5pt_wk(:)
  !$acc declare create(gauss_legendre_5pt_xk, gauss_legendre_5pt_wk)
  real(kind=dp_t), allocatable :: gauss_laguerre_5pt_xk(:)
  real(kind=dp_t), allocatable :: gauss_laguerre_5pt_wk(:)
  !$acc declare create(gauss_laguerre_5pt_xk, gauss_laguerre_5pt_wk)

  contains
    subroutine integration_init()
      use bl_types
      implicit none

      allocate(gauss_legendre_5pt_xk(5))
      gauss_legendre_5pt_xk = (/ &
        -0.906179845938664d0, &
        -0.538469310105683d0, &
        0.0d0, &
        0.538469310105683d0, &
        0.906179845938664d0 /)
      allocate(gauss_legendre_5pt_wk(5))
      gauss_legendre_5pt_wk = (/ &
        0.236926885056189d0, &
        0.478628670499366d0, &
        0.568888888888889d0, &
        0.478628670499366d0, &
        0.236926885056189d0 /)
      !$acc update device(gauss_legendre_5pt_xk, gauss_legendre_5pt_wk)
      allocate(gauss_laguerre_5pt_xk(5))
      gauss_laguerre_5pt_xk = (/ &
        0.263560319718141d0, &
        0.141340305910652d1, &
        0.359642577104072d1, &
        0.708581000585884d1, &
        0.126408008442758d2 /)
      allocate(gauss_laguerre_5pt_wk(5))
      gauss_laguerre_5pt_wk = (/ &
        0.521755610582809d0, &
        0.398666811083176d0, &
        0.759424496817076d-1, &
        0.361175867992205d-2, &
        0.233699723857762d-4 /)
      !$acc update device(gauss_laguerre_5pt_xk, gauss_laguerre_5pt_wk)
    end subroutine integration_init

end module integration
