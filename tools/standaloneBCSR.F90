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

  type(ARGMT) :: args
  type(BCSR)  :: A

  call StartUp(Args, "standalone")

  call New(A)
  call Get(A, TrixFile("OrthoF"), standalone_O = .false.)
  call Put(A, TrixFile("OrthoF-standalone"), standalone_O = .true.)

  call Delete(A)

end program standaloneBCSR
