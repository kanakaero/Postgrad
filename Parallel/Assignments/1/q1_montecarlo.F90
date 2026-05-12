! Monte Carlo Integration in Fortran
! Integral = ∫_0^4 ∫_0^3 ∫_0^2 (4x^3 + xy^2 + 5y + yz +6z) dx dy dz

program montecarlo_integral
   use omp_lib
   implicit none

   integer :: i, N
   integer :: nthreads
   character(len=32) :: arg1, arg2

   real(8) :: x, y, z
   real(8) :: r1, r2, r3
   real(8) :: sum, sum2
   real(8) :: favg, f2, frms
   real(8) :: volume
   real(8) :: integral
   real(8) :: start, finish, runtime
   real(8) :: fval

   ! Read command line arguments
   call get_command_argument(1, arg1)
   call get_command_argument(2, arg2)

   read (arg1, *) N
   read (arg2, *) nthreads

   call omp_set_num_threads(nthreads)

   sum = 0.0d0
   sum2 = 0.0d0

   volume = 24.0d0

   start = omp_get_wtime()

   call random_seed()

!$omp parallel private(i,x,y,z,r1,r2,r3,fval) reduction(+:sum,sum2)
!$omp do
   do i = 1, N

      call random_number(r1)
      call random_number(r2)
      call random_number(r3)

      x = 4.0d0*r1
      y = 3.0d0*r2
      z = 2.0d0*r3

      fval = 4d0*x**3 + x*y**2 + 5d0*y + y*z + 6d0*z

      sum = sum + fval
      sum2 = sum2 + fval*fval

   end do
!$omp end do
!$omp end parallel

   finish = omp_get_wtime()
   runtime = finish - start

   favg = sum/N
   f2 = sum2/N

   frms = sqrt(f2 - favg**2)

   integral = volume*favg

   print *
   print *, "Samples =", N
   print *, "Threads =", nthreads
   print *, "Integral =", integral
   print *, "Error estimate =", volume*frms/sqrt(real(N))
   print *, "Runtime =", runtime
   print *

   ! Append timing results
   open (unit=20, file="montecarlo_timing.dat", status="unknown", position="append")
   write (20, *) N, nthreads, runtime
   close (20)

   ! Append Monte Carlo statistics
    open(unit=30,file="montecarlo_stats.dat",status="unknown",position="append")
    write(30,'(I12,3E20.10)') N, favg, f2, frms
    close(30)

end program montecarlo_integral
