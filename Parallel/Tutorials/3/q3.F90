! Find Prime numbers between two given numbers using OpenMP and Sieve of Eratosthenes

program PrimeNumbers
   use omp_lib
   implicit none
   integer :: lower, upper, i, j, nthreads
   logical, allocatable :: is_prime(:)
   character(100) :: argument

   ! Define the range
   lower = 1
   upper = 100

   call GET_COMMAND_ARGUMENT(1, argument)
   read (argument, *) nthreads

   ! Allocate the is_prime array
   allocate (is_prime(1:upper))

   ! Initialize the is_prime array
   is_prime = .true.
   is_prime(1) = .false. ! 1 is not a prime number

   !$OMP PARALLEL DO NUM_THREADS(nthreads)
   do i = 2, upper
      if (is_prime(i)) then
         do j = i*2, upper, i

            print *, j
            is_prime(j) = .false.
         end do
      end if
   end do
   !$OMP END PARALLEL DO

   ! Print the prime numbers in the specified range
   print *, "Prime numbers between ", lower, " and ", upper, ":"
   do i = lower, upper
      if (is_prime(i)) then
         print *, i
      end if
   end do

end program PrimeNumbers
