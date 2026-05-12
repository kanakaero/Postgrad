! Fortran program to implement matrix multiplication using OpenACC and compare the results with a sequential version.

program q2
   implicit none
   integer, parameter :: N = 3
   integer :: i, j, k
   real :: a(N, N), b(N, N), c(N, N), d(N, N)

   ! Initialize matrices a and b with random values
   call random_seed()
   call random_number(a)
   call random_number(b)

   write (*, *) "Matrix A:"
   do i = 1, N
      write (*, *) (a(i, j), j=1, N)
   end do

   write (*, *) "Matrix B:"
   do i = 1, N
      write (*, *) (b(i, j), j=1, N)
   end do

   ! Sequential matrix multiplication
   do i = 1, N
      do j = 1, N
         c(i, j) = 0.0
         do k = 1, N
            c(i, j) = c(i, j) + a(i, k)*b(k, j)
         end do
      end do
   end do

   write (*, *) "Sequential matrix multiplication Result:"
   do i = 1, N
      write (*, *) (c(i, j), j=1, N)
   end do

   ! Parallel matrix multiplication using OpenACC
   !$acc parallel loop collapse(2)
   do i = 1, N
      do j = 1, N
         d(i, j) = 0.0
         do k = 1, N
            d(i, j) = d(i, j) + a(i, k)*b(k, j)
         end do
      end do
   end do

   write (*, *) "Parallel matrix multiplication Result:"
   do i = 1, N
      write (*, *) (d(i, j), j=1, N)
   end do

   ! Compare the results of the sequential and parallel versions
   if (all(c == d)) then
      print *, "The results are the same."
   else
      print *, "The results are different."
   end if

end program q2
