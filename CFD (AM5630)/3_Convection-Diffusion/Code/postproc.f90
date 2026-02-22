! ----------------------------------------------------------------------
! Subroutine to post-process the results
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 21st November 2025
! ----------------------------------------------------------------------

subroutine postproc()
    use var_dec
    implicit none

    real(dp) :: max_diff
    real(dp) :: rho, k_by_cp  ! temporary variables
    real(dp) :: dTdx_face, dTdy_face
    real(dp) :: Fw, Fe, Fs, Fn
    real(dp) :: net_flux_gsd, net_flux_tdma

    rho      = 1.0_dp
    k_by_cp  = 1.0_dp/50.0_dp

    max_diff = 0.0_dp
    do j=1,jmax
        do i=1,imax
            max_diff = max(max_diff, abs(T_gsd(i,j)-T_tdma(i,j)))
        end do
    end do
    print *, "max difference GS-TDMA =", max_diff
    print *, " "

    ! Net flux - Gauss Seidel
    
    Fw = 0.0_dp
    Fe = 0.0_dp
    Fs = 0.0_dp
    Fn = 0.0_dp

    ! Inlet (i = 1)
    do j = 1, jmax
        dTdx_face = (T_gsd(2,j) - T_gsd(1,j)) / (xc(2)-xc(1))
        Fw = Fw + (rho*u(1,j)*T_gsd(1,j) - k_by_cp * dTdx_face) * dy(j)
    end do

    ! Outlet (i = imax)
    do j = 1, jmax
        dTdx_face = (T_gsd(imax,j) - T_gsd(imax-1,j)) / (xc(imax)-xc(imax-1))
        Fe = Fe + (rho*u(imax,j)*T_gsd(imax,j) - k_by_cp * dTdx_face) * dy(j)
    end do

    ! Bottom (j = 1)
    do i = 1, imax
        dTdy_face = (T_gsd(i,2) - T_gsd(i,1)) / (yc(2)-yc(1))
        Fs = Fs + (rho*v(i,1)*T_gsd(i,1) - k_by_cp * dTdy_face) * dx(i)
    end do

    ! Top (j = jmax)
    do i = 1, imax
        dTdy_face = (T_gsd(i,jmax) - T_gsd(i,jmax-1)) / (yc(jmax)-yc(jmax-1))
        Fn = Fn + (rho*v(i,jmax)*T_gsd(i,jmax) - k_by_cp * dTdy_face) * dx(i)
    end do

    net_flux_gsd = Fw + Fe + Fs + Fn

    print *, 'Fw (west)  =', Fw
    print *, 'Fe (east)  =', Fe
    print *, 'Fs (south) =', Fs
    print *, 'Fn (north) =', Fn
    print *, 'Sum        =', Fw+Fe+Fs+Fn
    print *, " "

    print *, "Net convective + conductive flux = ", net_flux_gsd, " (Gauss Seidel)"
    print *, " "

    ! Net flux - TDMA

    Fw = 0.0_dp
    Fe = 0.0_dp
    Fs = 0.0_dp
    Fn = 0.0_dp

    ! Inlet (i = 1)
    do j = 1, jmax
        dTdx_face = (T_tdma(2,j) - T_tdma(1,j)) / (xc(2)-xc(1))
        Fw = Fw + (rho*u(1,j)*T_tdma(1,j) - k_by_cp * dTdx_face) * dy(j)
    end do

    ! Outlet (i = imax)
    do j = 1, jmax
        dTdx_face = (T_tdma(imax,j) - T_tdma(imax-1,j)) / (xc(imax)-xc(imax-1))
        Fe = Fe + (rho*u(imax,j)*T_tdma(imax,j) - k_by_cp * dTdx_face) * dy(j)
    end do

    ! Bottom (j = 1)
    do i = 1, imax
        dTdy_face = (T_tdma(i,2) - T_tdma(i,1)) / (yc(2)-yc(1))
        Fs = Fs + (rho*v(i,1)*T_tdma(i,1) - k_by_cp * dTdy_face) * dx(i)
    end do

    ! Top (j = jmax)
    do i = 1, imax
        dTdy_face = (T_tdma(i,jmax) - T_tdma(i,jmax-1)) / (yc(jmax)-yc(jmax-1))
        Fn = Fn + (rho*v(i,jmax)*T_tdma(i,jmax) - k_by_cp * dTdy_face) * dx(i)
    end do

    net_flux_tdma = Fw + Fe + Fs + Fn

    print *, 'Fw (west)  =', Fw
    print *, 'Fe (east)  =', Fe
    print *, 'Fs (south) =', Fs
    print *, 'Fn (north) =', Fn
    print *, 'Sum        =', Fw+Fe+Fs+Fn
    print *, " "

    print *, "Net convective + conductive flux = ", net_flux_tdma, " (TDMA)"
    print *, " "

end subroutine postproc