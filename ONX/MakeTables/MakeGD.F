      program MakeGD
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      integer itmp(50000)
      integer id(3)

      id(1) = 101
      id(2) = 201
      id(3) = 202
      do i=1,3
         do j=1,3
            ib = id(i)
            ik = id(j)
            call gcmaker(itmp,ib,ik,i1,i2,i3,i4,i5)
         enddo
      enddo
      stop
      end

      
