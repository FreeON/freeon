program BCSRInfo

  use DerivedTypes
  use InOut

  character(len = 2000) :: filename
  integer :: i, j
  type(BCSR) :: A
  type(DBL_RNK2) :: ADense

  if(command_argument_count() /= 1) then
    write(*,*) "need a filename"
    stop
  endif

  call get_command_argument(1, filename)
  call Get(A, filename)
  call SetEq(ADense, A)

  do i = 1, size(ADense%D, 1)
    do j = 1, size(ADense%D, 2)
      write(*,*) i, j, ADense%D(i, j)
    enddo
  enddo

end program BCSRInfo
