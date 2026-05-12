! Riemann Integration in Fortran
! Integral = ∫_0^4 ∫_0^3 ∫_0^2 (4x^3 + xy^2 + 5y + yz + 6z) dz dy dx

program riemann_integral
   use omp_lib
   implicit none

   integer :: i, j, k
   integer :: Nx, Ny, Nz
   integer :: nthreads
   character(len=10) :: arg

   real(8) :: dx, dy, dz
   real(8) :: x, y, z
   real(8) :: sum, integral
   real(8) :: exact, error
   real(8) :: start, finish, runtime

    ! Grid resolution
   Nx = 48
   Ny = 36
   Nz = 24

    ! Read number of threads from command line
   call get_command_argument(1, arg)
   read (arg, *) nthreads
   call omp_set_num_threads(nthreads)

    ! Grid spacing
   dx = 4.0d0/Nx
   dy = 3.0d0/Ny
   dz = 2.0d0/Nz

   sum = 0.0d0

   start = omp_get_wtime()

    ! Parallel Riemann sum
!$omp parallel do collapse(3) reduction(+:sum)
   do i = 1, Nx
      x = (i - 0.5d0)*dx
      do j = 1, Ny
         y = (j - 0.5d0)*dy
         do k = 1, Nz
            z = (k - 0.5d0)*dz

            sum = sum + (4d0*x**3 + x*y**2 + 5d0*y + y*z + 6d0*z)

         end do
      end do
   end do
!$omp end parallel do

   integral = sum*dx*dy*dz

   finish = omp_get_wtime()
   runtime = finish - start

    ! Analytical solution
   exact = 2040.0d0
   error = abs(integral - exact)/exact

   print *
   print *, "Threads =", nthreads
   print *, "Numerical Integral =", integral
   print *, "Analytical Integral =", exact
   print *, "Relative Error =", error
   print *, "Execution Time =", runtime
   print *

    ! Append results to a file
   open (unit=10, file="riemann_timing.dat", status="unknown", position="append")
   write (10, *) nthreads, runtime
   close (10)

end program riemann_integral
