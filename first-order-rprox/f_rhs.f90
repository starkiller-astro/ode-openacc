! Note: we expect the input to be mass fractions -- this is different
! than the version in Maestro.  We convert to molar fractions here.
!
! Y = X/A

subroutine rhs(n, t, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use rates_module
  use rpar_indices

  implicit none
  !$acc routine seq
  !$acc routine(make_rates) seq
  !$acc routine(make_ydots) seq
  integer,         intent(IN   ) :: n, ipar
  real(kind=dp_t), intent(IN   ) :: t, y(n)
  real(kind=dp_t), intent(INOUT) :: rpar(n_rpar_comps)
  real(kind=dp_t), intent(  OUT) :: ydot(n)

  integer :: k

  real(kind=dp_t), parameter :: T2T9 = 1.0e-9_dp_t

  real(kind=dp_t) :: dens, t9
  real(kind=dp_t) :: ymol(nspec)

  ydot = ZERO

  ! several thermodynamic quantities come in via rpar
  dens = rpar(irp_dens)
  t9 = rpar(irp_temp)*T2T9

  ymol = y(1:nspec)/aion(1:nspec)

  ! build the rates; weak rates are the wk* variables
  call make_rates(t9, dens, ymol, rpar)

  ! set up the ODEs
  call make_ydots(ymol,t9,rpar,ydot(1:nspec))

  ! make them in terms of mass fractions
  ! dX/dt = dX/dY dY/dt
  ydot(1:nspec) = aion(1:nspec)*ydot(1:nspec)

  return

end subroutine rhs


subroutine make_rates(t9,dens,y,rpar)

  use bl_types
  use bl_constants_module
  use rates_module
  use network
  use rpar_indices
  
  implicit none
  !$acc routine seq
  !$acc routine( rate_p_c12_to_n13)
  !$acc routine( rate_f17_to_p_o16)
  !$acc routine( rate_f17_to_o17)
  !$acc routine( rate_p_f17_to_ne18)
  !$acc routine( rate_he4_he4_he4_to_c12)
  !$acc routine( rate_p_n14_to_o15)
  !$acc routine( rate_he4_ne18_to_p_na21)
  !$acc routine( rate_ne18_to_f18)
  !$acc routine( rate_ne19_to_f19)
  !$acc routine( rate_p_ne19_to_na20)
  !$acc routine( rate_he4_o14_to_p_f17)
  !$acc routine( rate_o14_to_n14)
  !$acc routine( rate_he4_o15_to_ne19)
  !$acc routine( rate_o15_to_n15)
  !$acc routine( rate_he4_o16_to_ne20)
  !$acc routine( rate_p_o16_to_f17)
  !$acc routine( rate_he4_si26_to_p_p29)
  !$acc routine( rate_he4_ti44_to_p_v47)
  !$acc routine( calc_tfactors)
  
  real(kind=dp_t), intent(in   ) :: t9, dens, y(nspec)
  real(kind=dp_t), intent(inout) :: rpar(n_rpar_comps)

  ! locally used rates
  real(kind=dp_t) :: rate,wk18ne,wk19ne
  real(kind=dp_t) :: r56pg,cutoni,r57decay,r56eff
  real(kind=dp_t) :: t9i32

  ! some numbers from appendix C in WW81; these should probably be
  ! updated with current rates
  real(kind=dp_t), parameter :: Lweak = 1.05d0, & ! this is for NS
!                                Lweak = 0.107d0, & ! this is for lower 
                                                    ! densities
                                la2 = ONE/FIFTEEN ! mean rate from 30s to 56ni
                                                  ! from p-capture and beta 
                                                  ! decays
                                
  type (temp_t) :: tfactors
  
  rpar(irp_dlambCNOdh1:n_rpar_comps) = ZERO ! other constants

  call calc_tfactors(t9, tfactors)

  ! some common parameters
  rpar(irp_rates-1+irLweak) = Lweak
  rpar(irp_rates-1+irla2)   = la2

  ! weak rates first
  !
  ! 14o(beta nu)14n
  call rate_o14_to_n14(tfactors,rate)
  rpar(irp_rates-1+irwk14o) = rate

  ! 15o(beta nu)15n
  call rate_o15_to_n15(tfactors,rate)
  rpar(irp_rates-1+irwk15o) = rate
  
  ! 17f(beta nu)17o
  call rate_f17_to_o17(tfactors,rate)
  rpar(irp_rates-1+irwk17f) = rate

  ! these weak rates aren't needed outside of this routine
  ! 18ne(beta nu)18f
  call rate_ne18_to_f18(tfactors,wk18ne) 
  ! 19ne(beta nu)19f
  call rate_ne19_to_f19(tfactors,wk19ne)

  ! 12c(p,g)13n
  call rate_p_c12_to_n13(tfactors,rate)
  rpar(irp_rates-1+irpg12c) = dens*rate

  ! triple alpha
  call rate_he4_he4_he4_to_c12(tfactors,rate)
  rpar(irp_rates-1+ir3a) = dens*dens*rate

  ! 17f(p,g)18ne
  call rate_p_f17_to_ne18(tfactors,rate)
  rpar(irp_rates-1+irpg17f) = dens*rate

  ! 17f(g,p)16o
  call rate_f17_to_p_o16(tfactors,rate)
  rpar(irp_rates-1+irgp17f) = rate

  ! 15o(a,g)19ne
  call rate_he4_o15_to_ne19(tfactors,rate)
  rpar(irp_rates-1+irag15o) = dens*rate

  ! 16o(a,g)20ne
  call rate_he4_o16_to_ne20(tfactors,rate)
  rpar(irp_rates-1+irag16o) = dens*rate

  ! 16o(p,g)17f
  call rate_p_o16_to_f17(tfactors,rate)
  rpar(irp_rates-1+irpg16o) = dens*rate

  ! 14o(a,p)17f
  call rate_he4_o14_to_p_f17(tfactors,rate)
  rpar(irp_rates-1+irap14o) = dens*rate

  ! limit CNO as minimum between 14n(p,g)15o and 15o(beta nu)15n
  ! we store the limited rate in irlambCNO; this is lambda_CNO in WW81
  call rate_p_n14_to_o15(tfactors,rate)
  rpar(irp_rates-1+irlambCNO) = min(rpar(irp_rates-1+irwk15o),rate*dens*y(ih1))
  if (rpar(irp_rates-1+irlambCNO) < rpar(irp_rates-1+irwk15o)) then
     rpar(irp_dlambCNOdh1) = rate*dens
  endif

  ! 22mg(...)30s
  ! check if this proceeds via p-captures or (a,p) reactions
  ! the Lweak is from WW81, eqn C15
  ! we store the rate in irlambda1; this is the lambda1 in WW81
  call rate_he4_si26_to_p_p29(tfactors,rate)
  rpar(irp_rates-1+irlambda1) = max(rpar(irp_rates-1+irLweak),dens*y(ihe4)*rate)
  if (rpar(irp_rates-1+irlambda1) > rpar(irp_rates-1+irLweak)) then
       rpar(irp_dlambda1dhe4) = dens*rate
  ! use the sign of rpar(irp_rates-1+irlambda1) to indicate the value of delta1 in WW81
  ! if delta1 = 1, then we multiply the rate by -1
       rpar(irp_rates-1+irlambda1) = -ONE*rpar(irp_rates-1+irlambda1)
  endif

  ! 30s(...) 56ni 
  ! check if this proceeds via p-captures or (a,p) reactions
  ! use 44ti(a,p)v47 as a typical limiting rate for the (a,p) process
  ! store this in irlambda2; this is lambda2 in WW81
  call rate_he4_ti44_to_p_v47(tfactors,rate)
  rpar(irp_rates-1+irlambda2) = max(rpar(irp_rates-1+irla2),dens*y(ihe4)*rate)
  if (rpar(irp_rates-1+irlambda2) > rpar(irp_rates-1+irla2)) then
       rpar(irp_dlambda2dhe4) = dens*rate
  ! use the sign of rpar(irp_rates-1+irlambda2) to indicate th value of delta2
  ! if delta2 = 1, then we multiply the rate by -1
       rpar(irp_rates-1+irlambda2) = -ONE*rpar(irp_rates-1+irlambda2)
  endif

  ! form s1 from WW81; branching ratio for 18ne beta decay (wk18ne) vs (a,p)
  ! store result in irs1
  ! 18ne(a,p)21na
  call rate_he4_ne18_to_p_na21(tfactors,rate)
  rpar(irp_rates-1+irs1) = wk18ne / (wk18ne + dens*y(ihe4)*rate)
  rpar(irp_drs1dhe4) = -rpar(irp_rates-1+irs1)*dens*rate &
       / (wk18ne + dens*y(ihe4)*rate)

  ! form r1 from WW81; ranching ratio for 19ne beta decay (wk19ne) vs (p,g)
  ! store result in irr1
  ! 19ne(p,g)20na
  call rate_p_ne19_to_na20(tfactors,rate)
  rpar(irp_rates-1+irr1) = wk19ne / (wk19ne + dens*y(ih1)*rate)
  rpar(irp_drr1dh1) = -rpar(irp_rates-1+irr1)*dens*rate &
       / (wk19ne + dens*y(ih1)*rate)


  !....
  !....  additional coding for proton capture on 56ni to heavier elements
  !....   kludge    56ni+56p -> 2 (56ni) at a rate given by min
  !....   of 56ni(pg) and 57cu decay rate
  !....
  !....  use 56ni rate from wallace and woosley 1981
  t9i32=tfactors%t9i*sqrt(tfactors%t9i)
  r56pg=dens*(1.29e-02_dp_t*exp(-4.897_dp_t*tfactors%t9i) &
       +7.065e+03_dp_t*exp(-20.33_dp_t*tfactors%t9i))*t9i32

  !....  use generic proton separation energy of 400 kev
  !....  8.02 -> 4.64
  !      cutoni=2.08d-10*dens*exp(8.02*t9m1)/t932
  cutoni=2.08e-10_dp_t*dens*exp(4.642_dp_t*tfactors%t9i)*t9i32
  r57decay=3.54_dp_t
  r56eff=min(r56pg,cutoni*r57decay)
!   rpar(irp_r56eff) = r56eff
!   if (r56eff < r56pg) rpar(irp_dr56effdt) = r57decay*dcutonidt
  rpar(irp_r56eff) = ZERO

end subroutine make_rates


subroutine make_ydots(ymol,t9,rpar,dydt)

  use bl_types
  use bl_constants_module
  use network
  use rpar_indices
  
  implicit none
  !$acc routine seq

  real(kind=dp_t), intent(IN   ) :: ymol(nspec), t9
  real(kind=dp_t), intent(INOUT) :: rpar(n_rpar_comps)
  real(kind=dp_t), intent(  OUT) :: dydt(nspec)
  
  integer :: irp_start
  real(kind=dp_t) :: delta1, delta2
  real(kind=dp_t) :: dens

  ! initialize
  dydt = ZERO
  dens = rpar(irp_dens)

  ! check to see if we are doing this with the t-derivatives
  ! if so, offset our starting index in rpar
  irp_start = irp_rates

  delta1 = ZERO; delta2 = ZERO
  ! figure out the delta's; we used negative rates to indicate delta=1
  if (rpar(irp_rates-1+irlambda1) < ZERO) then
     delta1 = ONE
     rpar(irp_rates-1+irlambda1) = -ONE*rpar(irp_rates-1+irlambda1)
  endif
  if (rpar(irp_rates-1+irlambda2) < ZERO) then
     delta2 = ONE
     rpar(irp_rates-1+irlambda2) = -ONE*rpar(irp_rates-1+irlambda2)
  endif
  rpar(irp_delta1) = delta1
  rpar(irp_delta2) = delta2

! setup ODEs
!
!....
!.... 12c = 1
!....
      dydt(ic12)=-ymol(ic12)*ymol(ih1)*rpar(irp_start-1+irpg12c) &
           +ymol(ihe4)**3*rpar(irp_start-1+ir3a)/SIX &
           +ymol(io15)*rpar(irp_start-1+irlambCNO)
!....
!.... 14o = 2
!....
      dydt(io14)=-ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o) &
           -ymol(io14)*rpar(irp_start-1+irwk14o) &
           +ymol(ic12)*ymol(ih1)*rpar(irp_start-1+irpg12c)
!....
!.... 15o = 3
!....
      dydt(io15)=ymol(io14)*rpar(irp_start-1+irwk14o) &
           -ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o) &
           -ymol(io15)*rpar(irp_start-1+irlambCNO) &
           +ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*rpar(irp_start-1+irs1) &
           +ymol(if17)*rpar(irp_start-1+irwk17f)
!....
!.... 16o = 4
!....
      dydt(io16) = ymol(if17)*rpar(irp_start-1+irgp17f) &
           -ymol(io16)*ymol(ih1)*rpar(irp_start-1+irpg16o) &
           +ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irr1)*rpar(irp_start-1+irag15o) &
           -ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o)
!....
!.... 17f = 5
!....
      dydt(if17)=ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o) &
           +ymol(io16)*ymol(ih1)*rpar(irp_start-1+irpg16o) &
           -ymol(if17)*rpar(irp_start-1+irgp17f) &
           -ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f) &
           -ymol(if17)*rpar(irp_start-1+irwk17f)
!....
!.... 22mg = 6
!....
      dydt(img22)=ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o) &
           +ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*(ONE-rpar(irp_start-1+irs1)) &
           +ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o)*(ONE-rpar(irp_start-1+irr1)) &
           -ymol(img22)*rpar(irp_start-1+irlambda1)
!....
!.... 30s = 7
!....
      dydt(is30)=ymol(img22)*rpar(irp_start-1+irlambda1) &
           -ymol(is30)*rpar(irp_start-1+irlambda2)
!....
!.... amax (56ni) = 8  (note that WW81 have a typo -- they write lambda1 here)
!....
      dydt(ini56)=ymol(is30)*rpar(irp_start-1+irlambda2)
!....
!.... 4he (alpha) = 9
!....
      dydt(ihe4)=-ymol(ihe4)**3*HALF*rpar(irp_start-1+ir3a) &
           +ymol(io15)*rpar(irp_start-1+irlambCNO) &
           -ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o) &
           +ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*rpar(irp_start-1+irs1) &
           -ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o)*(ONE-rpar(irp_start-1+irr1)) &
           -ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o) &
           -ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*(ONE-rpar(irp_start-1+irs1)) &
           +ymol(if17)*rpar(irp_start-1+irwk17f) &
           -TWO*ymol(img22)*rpar(irp_start-1+irlambda1)*rpar(irp_delta1) &
           -6.5d0*ymol(is30)*rpar(irp_start-1+irlambda2)*rpar(irp_delta2)
!....
!.... 1h (p) = 10
!....
      dydt(ih1)=-ymol(io14)*rpar(irp_start-1+irwk14o) &
           -ymol(io15)*rpar(irp_start-1+irlambCNO) &
           -TWO*ymol(ic12)*ymol(ih1)*rpar(irp_start-1+irpg12c) &
           +ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o) &
           -TWO*ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*rpar(irp_start-1+irs1) &
           +ymol(if17)*rpar(irp_start-1+irgp17f) &
           -ymol(io16)*ymol(ih1)*rpar(irp_start-1+irpg16o) &
           -ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o)*rpar(irp_start-1+irr1) &
           -TWO*ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o) &
           -THREE*ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o)*(ONE-rpar(irp_start-1+irr1)) &
           -ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*(ONE-rpar(irp_start-1+irs1)) &
           -TWO*ymol(if17)*rpar(irp_start-1+irwk17f) &
           -EIGHT*ymol(img22)*rpar(irp_start-1+irlambda1)*(ONE-rpar(irp_delta1)) &
           -26.e0_dp_t*ymol(is30)*rpar(irp_start-1+irlambda2)*(ONE-rpar(irp_delta2))


      dydt(ini56) = dydt(ini56)+ymol(ini56)*ymol(ih1)*rpar(irp_r56eff)
      dydt(ih1) = dydt(ih1)-56.0d0*ymol(ini56)*ymol(ih1)*rpar(irp_r56eff)

end subroutine make_ydots


