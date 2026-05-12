! Simpson's Rule for numerical integration

program simpsons_rule
   use omp_lib
   implicit none
   integer :: i, size_index, nthreads
   real :: a, b, h, sum, integral
   character(100) :: argument
   integer, parameter :: sizes(5) = [32, 64, 128, 256, 512]

   ! Get the limits of integration and number of subintervals
   print *, "Enter the lower limit (a):"
   read *, a
   print *, "Enter the upper limit (b):"
   read *, b

   call GET_COMMAND_ARGUMENT(1, argument)
   read (argument, *) nthreads

   do size_index = 1, size(sizes)
      h = (b - a)/sizes(size_index) ! Calculate the width of each subinterval

      sum = f(a) + f(b) ! Initialize the sum with the first and last terms
      !$OMP PARALLEL DO REDUCTION(+:sum) NUM_THREADS(nthreads)
      do i = 1, sizes(size_index) - 1
         if (mod(i, 2) == 0) then
            sum = sum + 2*f(a + i*h)
         else
            sum = sum + 4*f(a + i*h)
         end if
      end do
      !$OMP END PARALLEL DO
      integral = (h/3)*sum
      print *, "The integral for size ", sizes(size_index), " is approximately:", integral
   end do

contains
   real function f(x)
      real :: x
      ! Define the function to integrate here
      f = sin(x)  ! Example: integrating sin(x)
   end function f
end program simpsons_rule
