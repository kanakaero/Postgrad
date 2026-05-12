program q1p
   use mpi
   implicit none

   integer :: ierr, rank, nprocs
   integer :: nx_global, nx_local, nt
   integer :: i, n, left, right
   integer, allocatable :: counts(:), displs(:)
   integer :: base, rem, offset

   real(8), parameter :: c = 1.0d0, L = 2.0d0
   real(8), parameter :: dx = 0.002d0, dt = 0.0001d0
   real(8), parameter :: t_end = 1.0d0
   real(8), parameter :: pi = 3.141592653589793d0
   real(8) :: courant, t

   real(8), allocatable :: x(:)
   real(8), allocatable :: u_up(:), uo_up(:)
   real(8), allocatable :: u_q(:), uo_q(:)

   real(8), allocatable :: xg(:), upg(:), qg(:), uex(:)

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

   nx_global = int(L/dx) + 1
   nt = int(t_end/dt)
   courant = c*dt/dx

   allocate (counts(nprocs), displs(nprocs))

   base = nx_global/nprocs
   rem = mod(nx_global, nprocs)

   do i = 0, nprocs - 1
      counts(i + 1) = base
      if (i .lt. rem) counts(i + 1) = base + 1
   end do

   displs(1) = 0
   do i = 2, nprocs
      displs(i) = displs(i - 1) + counts(i - 1)
   end do

   nx_local = counts(rank + 1)
   offset = displs(rank + 1)

   allocate (x(nx_local))
   allocate (u_up(0:nx_local + 1), uo_up(0:nx_local + 1))
   allocate (u_q(0:nx_local + 1), uo_q(0:nx_local + 1))

   if (rank .eq. 0) then
      allocate (xg(nx_global), upg(nx_global), qg(nx_global), uex(nx_global))
   end if

   left = rank - 1
   right = rank + 1
   if (left .lt. 0) left = MPI_PROC_NULL
   if (right .ge. nprocs) right = MPI_PROC_NULL

   do i = 1, nx_local
      x(i) = (offset + i - 1)*dx
   end do

   do i = 1, nx_local
      if (x(i) .le. 0.5d0) then
         u_up(i) = sin(4d0*pi*x(i))
      else
         u_up(i) = 0d0
      end if
      u_q(i) = u_up(i)
   end do

   u_up(0) = 0d0; u_up(nx_local + 1) = 0d0
   u_q(0) = 0d0; u_q(nx_local + 1) = 0d0

   call gather_and_write("t0p.dat", 0d0)

   do n = 1, nt
      t = n*dt
      uo_up = u_up
      uo_q = u_q

      call MPI_Sendrecv(uo_up(1), 1, MPI_DOUBLE_PRECISION, left, 1, &
                        uo_up(nx_local + 1), 1, MPI_DOUBLE_PRECISION, right, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call MPI_Sendrecv(uo_up(nx_local), 1, MPI_DOUBLE_PRECISION, right, 2, &
                        uo_up(0), 1, MPI_DOUBLE_PRECISION, left, 2, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call MPI_Sendrecv(uo_q(1), 1, MPI_DOUBLE_PRECISION, left, 3, &
                        uo_q(nx_local + 1), 1, MPI_DOUBLE_PRECISION, right, 3, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call MPI_Sendrecv(uo_q(nx_local), 1, MPI_DOUBLE_PRECISION, right, 4, &
                        uo_q(0), 1, MPI_DOUBLE_PRECISION, left, 4, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      do i = 1, nx_local
         u_up(i) = uo_up(i) - courant*(uo_up(i) - uo_up(i - 1))
      end do

      u_q(1) = uo_q(1) - courant*(uo_q(1) - uo_q(0))
      if (nx_local .ge. 2) u_q(2) = uo_q(2) - courant*(uo_q(2) - uo_q(1))

      do i = 3, nx_local - 1
         u_q(i) = uo_q(i) - courant*( &
                  0.375d0*uo_q(i + 1) + 0.375d0*uo_q(i) &
                  - 0.875d0*uo_q(i - 1) + 0.125d0*uo_q(i - 2))
      end do

      if (nx_local .ge. 2) u_q(nx_local) = uo_q(nx_local) - courant*(uo_q(nx_local) - uo_q(nx_local - 1))

      if (rank .eq. 0) then
         u_up(1) = 0d0
         u_q(1) = 0d0
      end if

      if (rank .eq. nprocs - 1) then
         u_up(nx_local) = 0d0
         u_q(nx_local) = 0d0
      end if

      if (abs(t - 0.5d0) .lt. 0.5d0*dt) call gather_and_write("t05p.dat", t)
      if (abs(t - 1.0d0) .lt. 0.5d0*dt) call gather_and_write("t1p.dat", t)

   end do

   call MPI_Finalize(ierr)

contains

   subroutine gather_and_write(fname, t)
      character(len=*), intent(in)::fname
      real(8), intent(in)::t
      integer::i, unit

      call MPI_Gatherv(x, nx_local, MPI_DOUBLE_PRECISION, &
                       xg, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_Gatherv(u_up(1), nx_local, MPI_DOUBLE_PRECISION, &
                       upg, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_Gatherv(u_q(1), nx_local, MPI_DOUBLE_PRECISION, &
                       qg, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      if (rank .eq. 0) then
         call compute_exact(uex, xg, t, nx_global, c)
         open (newunit=unit, file=fname, status="unknown", position="append")
         do i = 1, nx_global
            write (unit, '(i4,4f15.6)') nprocs, xg(i), upg(i), qg(i), uex(i)
         end do
         close (unit)
      end if
   end subroutine

   subroutine compute_exact(u_exact, x, t, n, c)
      real(8), intent(out)::u_exact(n)
      real(8), intent(in)::x(n), t, c
      integer, intent(in)::n
      integer::i
      real(8)::x0
      do i = 1, n
         x0 = x(i) - c*t
         if (x0 .ge. 0d0 .and. x0 .le. 0.5d0) then
            u_exact(i) = sin(4d0*pi*x0)
         else
            u_exact(i) = 0d0
         end if
      end do
   end subroutine

end program
