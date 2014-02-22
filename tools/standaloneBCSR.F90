!> This program converts the old style BCSR format, which does not contain basis
!! information, and can only be read in conjunction with the basis set, into the
!! new standalone BCSR format, which contains all necessary information, and can
!! be read from other other languages.
!!
!! The program is called as a standin replacement for SP2. The command line
!! arguments are identical to the ones used for the SP2 call.

program standaloneBCSR

  use DerivedTypes
  use InOut
  use Macros

  implicit none

  type(ARGMT) :: Args
  type(BCSR)  :: F, P

  call StartUp(Args, "standalone")

  call New(F)
  call Get(F, TrixFile("OrthoF", Args, 0), standalone_O = .false.)
  call Put(F, TrixFile("OrthoF-standalone", Args, 0), standalone_O = .true.)
  call Delete(F)

  call New(P)
  call Get(P, TrixFile("OrthoD", Args, 1), standalone_O = .false.)
  call Put(P, TrixFile("OrthoD-standalone", Args, 1), standalone_O = .true.)
  call Delete(P)

end program standaloneBCSR
