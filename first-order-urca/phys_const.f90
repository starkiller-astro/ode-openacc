module physical_constants_module
  use bl_types

  implicit none

  ! Below are from PDG 2014
  real(kind=dp_t), parameter  ::  MELEC_CONST = 0.510998928d0 ! MeV
  real(kind=dp_t), parameter  ::  KBOLT_CONST = 8.6173324d-11 ! MeV/K
  real(kind=dp_t), parameter  ::  MEV_PER_AMU_CONST = 931.494061d0 ! MeV/amu
  real(kind=dp_t), parameter  ::  ALPHA_CONST = 7.2973525698d-3

  ! Below are from AME 2012
  real(kind=dp_t), parameter  ::  MASS_NE23_CONST   = 22994466.91d-6 ! amu
  real(kind=dp_t), parameter  ::  MASS_NA23_CONST   = 22989769.2820d-6 ! amu
end module physical_constants_module
