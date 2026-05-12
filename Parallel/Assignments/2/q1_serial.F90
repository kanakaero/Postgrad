program q1s
   implicit none

   ! Parameters
   real(8), parameter :: c = 1.0d0
   real(8), parameter :: L = 2.0d0
   real(8), parameter :: dx = 0.002d0
   real(8), parameter :: dt = 0.0001d0
   real(8), parameter :: t_end = 1.0d0
   integer, parameter :: nx = int(L/dx) + 1
   integer, parameter :: nt = int(t_end/dt)

   ! Arrays
   real(8) :: x(nx)
   real(8) :: u_up(nx), u_quick(nx)
   real(8) :: uo_up(nx), uo_quick(nx)
   real(8) :: u_exact(nx)

   ! Variables
   integer :: i, n
   real(8) :: courant, t

   courant = c*dt/dx

   ! Grid
   do i = 1, nx
      x(i) = (i - 1)*dx
   end do

   ! Initial condition
   do i = 1, nx
      if (x(i) .le. 0.5d0) then
         u_up(i) = sin(4.0d0*3.14159265358979d0*x(i))
      else
         u_up(i) = 0.0d0
      end if
      u_quick(i) = u_up(i)
   end do

   ! initial condition
   call compute_exact(u_exact, x, 0.0d0, nx, c)
   call write_output("t0.dat", x, u_up, u_quick, u_exact, nx)

   ! Time stepping
   do n = 1, nt
      t = n*dt

      uo_up = u_up
      uo_quick = u_quick

      ! Upwind scheme
      do i = 2, nx
         u_up(i) = uo_up(i) - courant*(uo_up(i) - uo_up(i - 1))
      end do
      u_up(1) = 0.0d0

      ! QUICK scheme
      do i = 3, nx - 1
         u_quick(i) = uo_quick(i) - courant* &
                      (0.375d0*uo_quick(i + 1) + 0.375d0*uo_quick(i) &
                       - 0.875d0*uo_quick(i - 1) + 0.125d0*uo_quick(i - 2))
      end do

      ! Upwind where QUICK cannot be applied
      u_quick(2) = uo_quick(2) - courant*(uo_quick(2) - uo_quick(1))
      u_quick(1) = 0.0d0
      u_quick(nx) = 0.0d0

      ! Output
      if (abs(t - 0.5d0) .lt. 0.5d0*dt) then
         call compute_exact(u_exact, x, t, nx, c)
         call write_output("t05.dat", x, u_up, u_quick, u_exact, nx)
      end if

      if (abs(t - 1.0d0) .lt. 0.5d0*dt) then
         call compute_exact(u_exact, x, t, nx, c)
         call write_output("t1.dat", x, u_up, u_quick, u_exact, nx)
      end if

   end do

contains

   subroutine write_output(fname, x, u1, u2, uex, n)
      character(len=*), intent(in) :: fname
      real(8), intent(in) :: x(n), u1(n), u2(n), uex(n)
      integer, intent(in) :: n
      integer :: i, unit

      open (newunit=unit, file=fname, status='replace')
      do i = 1, n
         write (unit, '(4f15.6)') x(i), u1(i), u2(i), uex(i)
      end do
      close (unit)
   end subroutine

   subroutine compute_exact(u_exact, x, t, n, c)
      real(8), intent(out) :: u_exact(n)
      real(8), intent(in) :: x(n), t, c
      integer, intent(in) :: n
      integer :: i
      real(8) :: x0

      do i = 1, n
         x0 = x(i) - c*t
         if (x0 >= 0.0d0 .and. x0 .le. 0.5d0) then
            u_exact(i) = sin(4.0d0*3.14159265358979d0*x0)
         else
            u_exact(i) = 0.0d0
         end if
      end do
   end subroutine

end program
