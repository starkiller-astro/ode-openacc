! the network module provides the information about the species we are
! advecting:
!
! nspec      -- the number of species in the network
! nrat       -- the number of reactions in this network
!
! aion       -- atomic number
! zion       -- proton number
! ebin       -- nuclear binding energy
!
! spec_names -- the name of the isotope
! short_spec_names -- abbreviated names
!
! This module contains three routines:
!
!  network_init()         -- initialize the isotope properties
!
!  network_species_index  -- return the index of the species given its name
!

module network

  use bl_types

  implicit none

  character (len=*), parameter :: network_name = "rprox"

  integer, parameter :: nspec = 10, nrat = 18 !13
  integer, parameter :: naux  = 0
  !$acc declare copyin(nspec, nrat, naux)

  character (len=16), allocatable :: spec_names(:)
  character (len= 5), allocatable :: short_spec_names(:)
  !character (len= 5), save :: short_aux_names(naux)

  real(kind=dp_t), allocatable :: aion(:), zion(:), ebin(:)
  !$acc declare create(aion, zion, ebin)

  ! set the indices; ordering base on rprox.f
  integer, parameter :: ic12 = 1
  integer, parameter :: io14 = 2
  integer, parameter :: io15 = 3
  integer, parameter :: io16 = 4
  integer, parameter :: if17 = 5
  integer, parameter :: img22 = 6
  integer, parameter :: is30 = 7
  integer, parameter :: ini56 = 8
  integer, parameter :: ihe4 = 9
  integer, parameter :: ih1 = 10

  integer, parameter :: irlambCNO = 1
  integer, parameter :: irag15o = 2
  integer, parameter :: irr1 = 3
  integer, parameter :: irag16o = 4
  integer, parameter :: irpg16o = 5
  integer, parameter :: irpg17f = 6
  integer, parameter :: irgp17f = 7
  integer, parameter :: irlambda2 = 8
  integer, parameter :: irap14o = 9
  integer, parameter :: irs1 = 10
  integer, parameter :: irlambda1 = 11
  integer, parameter :: ir3a = 12
  integer, parameter :: irpg12c = 13
  integer, parameter :: irwk14o = 14
  integer, parameter :: irwk17f = 15
  integer, parameter :: irwk15o = 16
  integer, parameter :: irLweak = 17
  integer, parameter :: irla2 = 18

  !$acc declare copyin(ih1, ihe4, ic12, io14, io15, io16, if17, img22, is30, ini56)
  !$acc declare copyin(irlambCNO, irag15o, irr1, irag16o, irpg16o, irpg17f)
  !$acc declare copyin(irgp17f, irlambda2, irap14o, irs1, irlambda1)
  !$acc declare copyin(ir3a, irpg12c, irwk14o, irwk17f, irwk15o, irLweak, irla2)

contains

  subroutine network_init()

    use bl_constants_module

    real(kind=dp_t), parameter :: MeV2erg = 1.60217646e-6, &
                                  N_A = 6.0221415e23

    allocate(spec_names(nspec))
    spec_names(ic12) = "carbon-12"
    spec_names(io14) = "oxygen-14"
    spec_names(io15) = "oxygen-15"
    spec_names(io16) = "oxygen-16"
    spec_names(if17) = "flourine-17"
    spec_names(img22) = "magnesium-22"
    spec_names(is30) = "sulfur-30"
    spec_names(ini56) = "nickel-56"
    spec_names(ihe4) = "helium-4"
    spec_names(ih1) = "hydrogen-1"

    allocate(short_spec_names(nspec))
    short_spec_names(ic12) = "C12"
    short_spec_names(io14) = "O14"
    short_spec_names(io15) = "O15"
    short_spec_names(io16) = "O16"
    short_spec_names(if17) = "F17"
    short_spec_names(img22) = "Mg22"
    short_spec_names(is30) = "S30"
    short_spec_names(ini56) = "Ni56"
    short_spec_names(ihe4) = "He4"
    short_spec_names(ih1) = "H1"

    ! set the species properties
    allocate(aion(nspec))
    aion(ic12) = TWELVE
    aion(io14) = 14.0_dp_t
    aion(io15) = 15.0_dp_t
    aion(io16) = 16.0_dp_t
    aion(if17) = 17.0_dp_t
    aion(img22) = 22.0_dp_t
    aion(is30) = 30.0_dp_t
    aion(ini56) = 56.0_dp_t
    aion(ihe4) = FOUR
    aion(ih1) = ONE
    !$acc update device(aion)

    allocate(zion(nspec))
    zion(ic12) = SIX
    zion(io14) = EIGHT
    zion(io15) = EIGHT
    zion(io16) = EIGHT
    zion(if17) = NINE
    zion(img22) = TWELVE
    zion(is30) = 16.0_dp_t
    zion(ini56) = 28.0_dp_t
    zion(ihe4) = TWO
    zion(ih1) = ONE
    !$acc update device(zion)

    ! our convention is that binding energy is negative.  The
    ! following are the binding energy (MeV).
    allocate(ebin(nspec))
    ebin(ic12) =  92.16279_dp_t
    ebin(io14) =  98.7325_dp_t
    ebin(io15) = 111.9569_dp_t
    ebin(io16) = 127.6207_dp_t
    ebin(if17) = 128.2211_dp_t
    ebin(img22) = 168.5768_dp_t
    ebin(is30) = 243.6866_dp_t
    ebin(ini56) = 483.995_dp_t
    ebin(ihe4) = 28.29599_dp_t
    ebin(ih1) = ZERO

    ! convert to erg / g by multiplying by N_A / aion and converting to erg
    ebin = -ebin * N_A * MeV2erg / aion
    !$acc update device(ebin)

  end subroutine network_init


  function network_species_index(name)

    character (len=*) :: name
    integer :: network_species_index, n

    network_species_index = -1

    do n = 1, nspec
       if (name == spec_names(n) .or. name == short_spec_names(n)) then
          network_species_index = n
          exit
       endif
    enddo

    return
  end function network_species_index

  subroutine network_finalize()

  end subroutine network_finalize

end module network
