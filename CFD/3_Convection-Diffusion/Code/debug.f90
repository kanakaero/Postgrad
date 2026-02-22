subroutine debug()

    use var_dec
    implicit none

    real(dp), allocatable :: div(:,:)
    real(dp) :: dudx, dvdy
    real(dp) :: m_w, m_e, m_s, m_n, face_width, u_face, v_face
    real(dp) :: rho
    m_w = 0.0_dp; m_e = 0.0_dp; m_s = 0.0_dp; m_n = 0.0_dp
    rho = 1.0_dp

    ! West (i=1)
    do j = 1, jmax-1
        u_face = u(1,j)
        face_width = yg(j+1) - yg(j)
        m_w = m_w + rho * u_face * face_width
    end do

    ! East (i=imax)
    do j = 1, jmax-1
        u_face = u(imax,j)
        face_width = yg(j+1) - yg(j)
        m_e = m_e + rho * u_face * face_width
    end do

    ! South (j=1)
    do i = 1, imax-1
    v_face = v(i,1)
        face_width = xg(i+1) - xg(i)
        m_s = m_s + rho * v_face * face_width
    end do

    ! North (j=jmax)
    do i = 1, imax-1
        v_face = v(i,jmax)
        face_width = xg(i+1) - xg(i)
        m_n = m_n + rho * v_face * face_width
    end do

    print *, 'mass fluxes (m_w, m_e, m_s, m_n) =', m_w, m_e, m_s, m_n
    print *, 'net mass flux =', m_w + m_e + m_s + m_n
    print *, ' '

    allocate(div(imax, jmax))
    div = 0.0_dp

    ! central-difference divergence on interior cells
    do j = 2, jmax-1
        do i = 2, imax-1
            dudx = (u(i+1,j) - u(i-1,j))/(xc(i+1) - xc(i-1))
            dvdy = (v(i,j+1) - v(i,j-1))/(yc(j+1) - yc(j-1))
            div(i,j) = dudx + dvdy
        end do
    end do

    ! West and East boundaries (i=1 and i=imax)
    do j = 2, jmax-1
        div(1,j) = ((u(2,j) - u(1,j))/(xc(2)-xc(1))) + ((v(1,j+1) - v(1,j-1))/(yc(j+1)-yc(j-1)))
        div(imax,j) = ((u(imax,j) - u(imax-1,j))/(xc(imax)-xc(imax-1))) + ((v(imax,j+1) - v(imax,j-1))/( yc(j+1)-yc(j-1)))
    end do

    ! South and North boundaries (j=1 and j=jmax)
    do i = 2, imax-1
        div(i,1) = ((u(i+1,1)-u(i-1,1))/(xc(i+1)-xc(i-1))) + ((v(i,2)-v(i,1))/(yc(2)-yc(1)))
        div(i,jmax) = ((u(i+1,jmax)-u(i-1,jmax)) / (xc(i+1)-xc(i-1))) + ((v(i,jmax)-v(i,jmax-1))/( yc(jmax)-yc(jmax-1)))
    end do

    print *, 'Divergence diagnostics:'
    print *, ' min div =', minval(div)
    print *, ' max div =', maxval(div)

    deallocate(div)

end subroutine debug