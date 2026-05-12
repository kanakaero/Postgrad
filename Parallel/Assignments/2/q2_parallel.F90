program q2p
   use mpi
   implicit none

   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: Nx = 201, Ny = 201
   real(dp), parameter :: L = 2.0_dp
   real(dp), parameter :: pi = 3.141592653589793_dp
   real(dp), parameter :: tol = 1.0e-4_dp
   integer, parameter :: max_iter = 20000

   integer :: rank, nprocs, ierr
   integer :: i, j, iter
   integer :: nlocal, jstart
   real(dp) :: dx, dy, err, local_err

   real(dp), allocatable :: phi(:,:), phi_new(:,:)
   real(dp), allocatable :: phi_global(:,:)
   real(dp) :: x(Nx), y(Ny)

   integer :: midx, midy
   character(len=32) :: fname

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

   dx = L/(Nx-1)
   dy = dx

   do i = 1, Nx
      x(i) = -1.0_dp + (i-1)*dx
   end do
   do j = 1, Ny
      y(j) = -1.0_dp + (j-1)*dy
   end do

   ! -------- Row decomposition --------
   nlocal = Ny / nprocs
   jstart = rank*nlocal + 1

   allocate(phi(Nx, nlocal+2))
   allocate(phi_new(Nx, nlocal+2))

   phi = 0.0_dp
   phi_new = 0.0_dp

   ! Left boundary
   do j = 1, nlocal
      phi(1,j+1) = sin(2.0_dp*pi*y(jstart+j-1))
      phi_new(1,j+1) = phi(1,j+1)
   end do

   ! ------------- Iteration -------------
   do iter = 1, max_iter
      local_err = 0.0_dp

      ! Exchange halo rows
      if (rank .gt. 0) then
         call MPI_Sendrecv(phi(:,2),Nx,MPI_DOUBLE_PRECISION,rank-1,0, &
                           phi(:,1),Nx,MPI_DOUBLE_PRECISION,rank-1,1, &
                           MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      end if

      if (rank .lt. nprocs-1) then
         call MPI_Sendrecv(phi(:,nlocal+1),Nx,MPI_DOUBLE_PRECISION,rank+1,1, &
                           phi(:,nlocal+2),Nx,MPI_DOUBLE_PRECISION,rank+1,0, &
                           MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      end if

      ! Interior update
      do j = 2, nlocal+1
         do i = 2, Nx-1
            phi_new(i,j) = 0.25_dp*( &
               phi(i+1,j) + phi(i-1,j) + &
               phi(i,j+1) + phi(i,j-1) + &
               dx*dx*(x(i)**2 + y(jstart+j-2)**2) )
         end do
      end do

      ! Left boundary (Dirichlet)
      do j = 2, nlocal+1
         phi_new(1,j) = sin(2.0_dp*pi*y(jstart+j-2))
      end do

      ! Top & bottom (Dirichlet)
      if (rank .eq. 0) phi_new(:,2) = 0.0_dp
      if (rank .eq. nprocs-1) phi_new(:,nlocal+1) = 0.0_dp

      ! Right boundary (Neumann)
      do j = 2, nlocal+1
         phi_new(Nx,j) = (4.0_dp*phi_new(Nx-1,j) - phi_new(Nx-2,j))/3.0_dp
      end do

      ! Error
      do j = 2, nlocal+1
         do i = 2, Nx-1
            local_err = max(local_err, abs(phi_new(i,j)-phi(i,j)))
         end do
      end do

      call MPI_Allreduce(local_err,err,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                         MPI_COMM_WORLD,ierr)

      phi = phi_new
      if (err .lt. tol) exit
   end do

   if (rank .eq. 0) print *, "Iterations:", iter, " Error:", err

   ! -------- Gather full solution --------
   if (rank .eq. 0) allocate(phi_global(Nx,Ny))

   call MPI_Gather(phi(:,2), Nx*nlocal, MPI_DOUBLE_PRECISION, &
                   phi_global, Nx*nlocal, MPI_DOUBLE_PRECISION, &
                   0, MPI_COMM_WORLD, ierr)

   ! -------- Output --------
   if (rank .eq. 0) then
      midx = Nx/2 + 1
      midy = Ny/2 + 1

      write(fname,'("phi_y0_p",I0,".dat")') nprocs
      open(10,file=fname,status="replace")
      do i = 1, Nx
         write(10,*) x(i), phi_global(i,midy)
      end do
      close(10)

      write(fname,'("phi_x0_p",I0,".dat")') nprocs
      open(20,file=fname,status="replace")
      do j = 1, Ny
         write(20,*) y(j), phi_global(midx,j)
      end do
      close(20)
   end if

   call MPI_Finalize(ierr)
end program