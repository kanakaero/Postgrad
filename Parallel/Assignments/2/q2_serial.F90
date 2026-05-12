program q2s
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: Nx = 21, Ny = 21
  real(dp), parameter :: L = 2.0_dp
  real(dp), parameter :: pi = 3.141592653589793_dp
  real(dp), parameter :: tol = 1.0e-4_dp
  integer, parameter :: max_iter = 20000

  real(dp) :: dx, dy, err
  real(dp) :: phi(Nx,Ny), phi_new(Nx,Ny)
  real(dp) :: x(Nx), y(Ny)
  integer :: i, j, iter
  integer :: midx, midy

  ! Grid spacing
  dx = L/(Nx-1)
  dy = dx

  ! Coordinates: from -1 to 1
  do i = 1, Nx
     x(i) = -1.0_dp + (i-1)*dx
  end do

  do j = 1, Ny
     y(j) = -1.0_dp + (j-1)*dy
  end do

  ! Initial condition
  phi = 0.0_dp
  phi_new = 0.0_dp

  ! Left boundary: phi(-1,y) = sin(2πy)
  do j = 1, Ny
     phi(1,j) = sin(2.0_dp*pi*y(j))
     phi_new(1,j) = phi(1,j)
  end do

  ! Iteration
  do iter = 1, max_iter
     err = 0.0_dp

     ! Interior update
     do i = 2, Nx-1
        do j = 2, Ny-1
           phi_new(i,j) = 0.25_dp * ( &
                phi(i+1,j) + phi(i-1,j) + &
                phi(i,j+1) + phi(i,j-1) + &
                dx*dx*(x(i)**2 + y(j)**2) )
        end do
     end do

     ! Boundary conditions

     ! Left (Dirichlet)
     do j = 1, Ny
        phi_new(1,j) = sin(2.0_dp*pi*y(j))
     end do

     ! Bottom and Top (Dirichlet)
     phi_new(:,1)  = 0.0_dp
     phi_new(:,Ny) = 0.0_dp

     ! Right (Neumann: dφ/dx = 0)
     do j = 2, Ny-1
        phi_new(Nx,j) = (4.0_dp*phi_new(Nx-1,j) - phi_new(Nx-2,j)) / 3.0_dp
     end do

     ! Compute L-infinity norm
     do i = 2, Nx-1
        do j = 2, Ny-1
           err = max(err, abs(phi_new(i,j) - phi(i,j)))
        end do
     end do

     phi = phi_new
     if (err .lt. tol) exit
  end do

  print *, "Converged in iterations:", iter
  print *, "Final error:", err

  ! Midpoints for centerline output
  midx = Nx/2 + 1
  midy = Ny/2 + 1

  ! φ(x, y=0)
  open(10,file="phi_y0.dat",status="replace")
  do i = 1, Nx
     write(10,*) x(i), phi(i,midy)
  end do
  close(10)

  ! φ(x=0, y)
  open(20,file="phi_x0.dat",status="replace")
  do j = 1, Ny
     write(20,*) y(j), phi(midx,j)
  end do
  close(20)

  print *, "Output written: phi_y0.dat and phi_x0.dat"

end program q2s