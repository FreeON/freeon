program BCSRInfo

  use DerivedTypes
  use InOut

  type(BCSR) :: A

  call Get(A, "A.OrthoF")

  write(*,*) "NSMat = ", A%NSMat
  write(*,*) "NAtms = ", A%NAtms
  write(*,*) "NBlks = ", A%NBlks
  write(*,*) "NNon0 = ", A%NNon0

end program BCSRInfo
