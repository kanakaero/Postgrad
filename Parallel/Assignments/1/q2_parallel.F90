! Poisson Equation Solver using Jacobi Iteration

program jacobi_parallel
   use omp_lib
   implicit none

   integer :: i,j,k,N,threads
   real(8) :: dx,x,y,error,maxerr
   real(8), allocatable :: phi(:,:),phi_new(:,:),exact(:,:),q(:,:)
   character(len=32) :: arg1,arg2
   real(8) :: num,den
   real(8) :: t1,t2

! ---- Read command line arguments ----
   call get_command_argument(1,arg1)
   call get_command_argument(2,arg2)

   read(arg1,*) dx
   read(arg2,*) threads

   call omp_set_num_threads(threads)

! ---- Grid size ----
   N = int(2.0/dx) + 1

   allocate(phi(N,N),phi_new(N,N),exact(N,N),q(N,N))

   phi = 0.0d0

! ---- Initialize RHS and exact solution ----
!$omp parallel do private(i,x,y)
   do j=1,N
      y = -1 + (j-1)*dx
      do i=1,N
         x = -1 + (i-1)*dx
         q(i,j) = 2*(2 - x*x - y*y)
         exact(i,j) = (x*x - 1)*(y*y - 1)
      end do
   end do
!$omp end parallel do

   k = 0
   t1 = omp_get_wtime()

! ---- Jacobi iteration ----
   do

! ---- Jacobi update ----
!$omp parallel do private(i)
      do j=2,N-1
         do i=2,N-1

            phi_new(i,j) = 0.25d0 * (phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1)) + (dx*dx/4.0d0)*q(i,j)

         end do
      end do
!$omp end parallel do

      phi(2:N-1,2:N-1) = phi_new(2:N-1,2:N-1)

! ---- compute normalized L2 error ----
      num = 0.0d0
      den = 0.0d0

!$omp parallel do private(i) reduction(+:num,den)
      do j=2,N-1
         do i=2,N-1
            num = num + (phi(i,j)-exact(i,j))**2
            den = den + exact(i,j)**2
         end do
      end do
!$omp end parallel do

      maxerr = sqrt(num/den)

      k = k + 1

      if(mod(k,50)==0) then
         print*,"Iter =",k," Relative Error =",maxerr*100.0d0,"%"
      end if

      if(maxerr < 0.01d0) then
         print*,"Converged (within 1%)"
         exit
      end if

   end do

   t2 = omp_get_wtime()

   print*,"Iterations =",k
   print*,"Runtime =",t2-t1," seconds"

! ---- Save runtime vs threads ----
   open(20,file="runtime_vs_threads.dat",status="unknown",position="append")
   write(20,'(4E16.8)') dx, threads, t2-t1, k
   close(20)

! ---- Output slice y=0.5 ----
   open(10,file="phi_slice_all.dat",status="unknown",position="append")

   j = nint((0.5 + 1)/dx) + 1

   write(10,*) "# dx =",dx," threads =",threads

   do i=1,N
      x = -1 + (i-1)*dx
      write(10,'(3E20.10)') x,phi(i,j),exact(i,j)
   end do

   write(10,*)

   close(10)

! ---- Compute final max absolute error ----
   maxerr = 0.0d0

!$omp parallel do private(i,error) reduction(max:maxerr)
   do j=2,N-1
      do i=2,N-1
         error = abs(phi(i,j)-exact(i,j))
         maxerr = max(maxerr,error)
      end do
   end do
!$omp end parallel do

   print*,"Maximum Absolute Error =",maxerr

end program jacobi_parallel