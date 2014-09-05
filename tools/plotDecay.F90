program plot

  use derivedtypes
  use inout
  use linalg
  use macros

  implicit none

  character(len = *), parameter :: prog = "plotDecay"

  type(argmt) :: args
  type(crds) :: GM
  type(bcsr) :: A, X2

  call startup(args, prog)

  call get(GM, "1")

  call get(A, trixfile("S", args_o = args))
  call plotdecay(A, GM, "decay_S")

  call get(A, trixfile("X", args_o = args))
  call plotdecay(A, GM, "decay_X")

  call multiply(A, A, X2)
  call plotdecay(X2, GM, "decay_X2")

  call delete(A)
  call delete(X2)

  call shutdown(prog)

end program plot
