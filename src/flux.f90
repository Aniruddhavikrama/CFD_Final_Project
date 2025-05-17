module flux
    use set_precision
    use fluid_constants
    use set_inputs
    use set_constants
    use grid_type
    use soln_type
    use variable_conversion
    use limiters
    implicit none

    procedure(calc_flux), pointer :: flux_fun

    abstract interface
        subroutine calc_flux(VL,VR, nx, ny, F)
            import :: prec, neq
            real(prec), dimension(neq), intent(in) :: VL,VR
            real(prec), intent(in) :: nx, ny
            real(prec), dimension(neq), intent(out) :: F
        end subroutine calc_flux
    end interface

    private
    public :: compute_fluxes, vanleer_flux, roe_flux, compute_lr_states
    
contains
subroutine select_flux
    use set_inputs, only: flux_scheme
    flux_fun => null()
    select case (flux_scheme)
    ! case (0)
    !     flux_fun => central_flux
    case (1)
        flux_fun => vanleer_flux
    case (2)
        flux_fun => roe_flux
    case default
        write(*,*) 'ERROR: Invalid flux_scheme value'
        stop
    end select
end subroutine select_flux

subroutine compute_lr_states(U, Lxi, Rxi, Leta, Reta,grid)
    type(grid_t),intent(in)::grid
    real(prec), dimension(neq,grid%ig_low:grid%ig_high,grid%jg_low:grid%jg_high), intent(in) :: U
    real(prec), dimension(neq,grid%i_low:grid%i_high,grid%j_low:grid%j_high-1), intent(out) :: Lxi, Rxi
    real(prec), dimension(neq,grid%i_low:grid%i_high-1,grid%j_low:grid%j_high), intent(out) :: Leta, Reta
    real(prec), dimension(neq,grid%ig_low:grid%ig_high,grid%jg_low:grid%jg_high) :: V
    real(prec), dimension(neq,grid%ig_low:grid%ig_high,grid%jg_low:grid%jg_high) :: psi_plus_xi, psi_minus_xi
    real(prec), dimension(neq,grid%ig_low:grid%ig_high,grid%jg_low:grid%jg_high) :: psi_plus_eta, psi_minus_eta
    real(prec) :: epsilon
    integer :: i, j

    ! Set epsilon based on order
    if (order == 1) then
        epsilon = 0.0_prec
    else if (order == 2) then
        epsilon = 1.0_prec
    else
        write(*,*) 'ERROR! order must be 1 or 2!!!'
        stop
    end if

    ! Convert conservative to primitive variables
    do j = grid%jg_low, grid%jg_high
        do i = grid%ig_low, grid%ig_high
            V(1,i,j) = U(1,i,j)  ! rho
            V(2,i,j) = U(2,i,j) / U(1,i,j)  ! u
            V(3,i,j) = U(3,i,j) / U(1,i,j)  ! v
            V(4,i,j) = (gamma - one) * (U(4,i,j) - half * (U(2,i,j)**2 + U(3,i,j)**2) / U(1,i,j))  ! p
        end do
    end do
    ! call limit_primitives(V)

    ! Initialize limiter arrays
    psi_plus_xi = one
    psi_minus_xi = one
    psi_plus_eta = one
    psi_minus_eta = one

    ! Compute limiters for second-order
    if (order == 2) then
        call select_limiter(limiter_scheme)
        call calculate_limiters(grid, V, psi_plus_xi, psi_minus_xi)
        call calculate_limiters(grid, V, psi_plus_eta, psi_minus_eta, direction='eta') 
    end if

    ! do j = grid%j_low, grid%j_high-1
    !     Lxi(:,grid%i_low,j) = V(:,grid%i_low-1,j) + 0.5_prec * (V(:,grid%i_low,j) - V(:,grid%i_low-1,j))
    !     Rxi(:,grid%i_low,j) = V(:,grid%i_low,j) - 0.5_prec * (V(:,grid%i_low+1,j) - V(:,grid%i_low,j))
    ! end do

    ! xsi-direction extrapolation (faces at i_low:i_high)
    do j = grid%j_low, grid%j_high-1
        do i = grid%i_low, grid%i_high
            Lxi(:,i,j) = V(:,i-1,j) + fourth * epsilon * ( &
                (one - kappa) * psi_plus_xi(:,i-2,j) * (V(:,i-1,j) - V(:,i-2,j)) + &
                (one + kappa) * psi_minus_xi(:,i-1,j) * (V(:,i,j) - V(:,i-1,j)))
            Rxi(:,i,j) = V(:,i,j) - fourth * epsilon * ( &
                (one - kappa) * psi_plus_xi(:,i,j) * (V(:,i+1,j) - V(:,i,j)) + &
                (one + kappa) * psi_minus_xi(:,i-1,j) * (V(:,i,j) - V(:,i-1,j)))
        end do
    end do

    !! Boundary at i=grid%i_high
    ! do j = grid%j_low, grid%j_high-1
    !     Lxi(:,grid%i_high,j) = V(:,grid%i_high-1,j) + 0.5_prec * (V(:,grid%i_high,j) - V(:,grid%i_high-1,j))
    !     Rxi(:,grid%i_high,j) = V(:,grid%i_high,j) - 0.5_prec * (V(:,grid%i_high,j) - V(:,grid%i_high-1,j))
    ! end do

    ! ! etadirection extrapolation (faces at j_low:j_high)

    ! ! Boundary at j=grid%j_low
    ! do i = grid%i_low, grid%i_high-1
    !     Leta(:,i,grid%j_low) = V(:,i,grid%j_low-1) + 0.5_prec * (V(:,i,grid%j_low) - V(:,i,grid%j_low-1))
    !     Reta(:,i,grid%j_low) = V(:,i,grid%j_low) - 0.5_prec * (V(:,i,grid%j_low+1) - V(:,i,grid%j_low))
    ! end do

    !interior
    do j = grid%j_low, grid%j_high
        do i = grid%i_low, grid%i_high-1
            Leta(:,i,j) = V(:,i,j-1) + fourth * epsilon * ( &
                (one - kappa) * psi_plus_eta(:,i,j-2) * (V(:,i,j-1) - V(:,i,j-2)) + &
                (one + kappa) * psi_minus_eta(:,i,j-1) * (V(:,i,j) - V(:,i,j-1)))
            Reta(:,i,j) = V(:,i,j) - fourth * epsilon * ( &
                (one - kappa) * psi_plus_eta(:,i,j) * (V(:,i,j+1) - V(:,i,j)) + &
                (one + kappa) * psi_minus_eta(:,i,j-1) * (V(:,i,j) - V(:,i,j-1)))
        end do
    end do

    ! Boundary at j=grid%j_high
    ! do i = grid%i_low, grid%i_high-1
    !     Leta(:,i,grid%j_high) = V(:,i,grid%j_high-1) + 0.5_prec * (V(:,i,grid%j_high) - V(:,i,grid%j_high-1))
    !     Reta(:,i,grid%j_high) = V(:,i,grid%j_high) - 0.5_prec * (V(:,i,grid%j_high) - V(:,i,grid%j_high-1))
    ! end do
end subroutine compute_lr_states
! subroutine compute_fluxes(grid, soln)
!     use set_precision, only: prec
!     type(grid_t), intent(in)    :: grid
!     type(soln_t), intent(inout) :: soln

!     real(prec), dimension(neq,grid%i_low:grid%i_high,grid%j_low:grid%j_high-1) :: Lxi, Rxi
!     real(prec), dimension(neq,grid%i_low:grid%i_high-1,grid%j_low:grid%j_high) :: Leta, Reta
!     real(prec) :: nx, ny
!     integer :: i, j

!     ! Debug indices: print fluxes at this location
!     integer, parameter :: debug_i = 5, debug_j = 5

!     call select_flux()
!     call compute_lr_states(soln%U, Lxi, Rxi, Leta, Reta, grid)

!     ! ξ-direction fluxes
!     do j = grid%j_low, grid%j_high-1
!         do i = grid%i_low, grid%i_high
!             nx = grid%n_xi(i,j,1)
!             ny = grid%n_xi(i,j,2)

!             call flux_fun(Lxi(:,i,j), Rxi(:,i,j), nx, ny, soln%Fxi(:,i,j))

!             if (i == debug_i .and. j == debug_j) then
!                 print *, '--- ξ-direction DEBUG ---'
!                 print *, 'i, j       =', i, j
!                 print *, 'Fxi        =', soln%Fxi(:,i,j)
!                 print *, 'n_xi       =', nx, ny
!                 print *, 'A_xi       =', grid%A_xi(i,j)
!                 print *, '---------------------------'
!             end if
!         end do
!     end do

!     ! η-direction fluxes
!     do j = grid%j_low, grid%j_high
!         do i = grid%i_low, grid%i_high-1
!             nx = grid%n_eta(i,j,1)
!             ny = grid%n_eta(i,j,2)

!             call flux_fun(Leta(:,i,j), Reta(:,i,j), nx, ny, soln%Feta(:,i,j))

!             if (i == debug_i .and. j == debug_j) then
!                 print *, '--- η-direction DEBUG ---'
!                 print *, 'i, j       =', i, j
!                 print *, 'Feta       =', soln%Feta(:,i,j)
!                 print *, 'n_eta      =', nx, ny
!                 print *, 'A_eta      =', grid%A_eta(i,j)
!                 print *, '---------------------------'
!             end if
!         end do
!     end do
! end subroutine compute_fluxes


subroutine compute_fluxes(grid, soln)
    type(grid_t), intent(in)    :: grid
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq,grid%i_low:grid%i_high,grid%j_low:grid%j_high-1) :: Lxi, Rxi
    real(prec), dimension(neq,grid%i_low:grid%i_high-1,grid%j_low:grid%j_high) :: Leta, Reta
    real(prec) :: nx, ny
    integer :: i, j

    ! Select flux scheme
    call select_flux()

    ! Compute left and right states at face indices
    call compute_lr_states(soln%U, Lxi, Rxi, Leta, Reta,grid)

    ! Limit primitive variables
    ! call limit_primitives(soln%U, Lxi)
    ! call limit_primitives(soln%U, Rxi)
    ! call limit_primitives(soln%U, Leta)
    ! call limit_primitives(soln%U, Reta)
    ! call limit_primitives(Lxi)
    ! call limit_primitives( Rxi)
    ! call limit_primitives( Leta)
    ! call limit_primitives( Reta)

    ! Compute ξ-direction fluxes
    do j = grid%j_low, grid%j_high-1
        do i = grid%i_low, grid%i_high
            nx = grid%n_xi(i,j,1)
            ny = grid%n_xi(i,j,2)
            call flux_fun(Lxi(:,i,j), Rxi(:,i,j), nx, ny, soln%Fxi(:,i,j))
        end do
    end do

    ! Compute η-direction fluxes
    do j = grid%j_low, grid%j_high
        do i = grid%i_low, grid%i_high-1
            nx = grid%n_eta(i,j,1)
            ny = grid%n_eta(i,j,2)
            ! print *, "nx and ny", nx,ny
            call flux_fun( Leta(:,i,j), Reta(:,i,j), nx, ny, soln%Feta(:,i,j) ) 
        end do
    end do
    end subroutine compute_fluxes


subroutine vanleer_flux(VL, VR,nx,ny,F)
    real(prec), dimension(neq), intent(in) :: VL, VR  ! Primitive variables: [rho, u,v, p]
    real(prec), dimension(neq), intent(out) :: F
    real(prec),intent(in):: nx, ny
    real(prec) :: rhoL, uL1, pL, aL, ML, rhoR, uR1, pR, aR, MR,vL1,vR1
    real(prec) :: Fc(neq), Fp(neq)
    real(prec) :: alpha_plus, alpha_minus, beta_plus, beta_minus
    real(prec) :: m_plus, m_minus, C_plus, C_minus, D_plus, D_minus
    real(prec) :: p_plus, p_minus
    real(prec) :: hL, hR
    real(prec) :: unL,unR
  
    ! Extract primitive variables from VL (left state)
    rhoL = VL(1)  ! Density
    uL1  = VL(2)  ! u-Velocity
    vl1 =  VL(3)   !v-velocity
    pL   = VL(4)  ! Pressure
  
    ! Compute speed of sound and Mach number for left state
    aL   = sqrt(gamma * pL / rhoL)  ! Speed of sound
    ! ML   = uL1 / aL  ! Mach number
  
    ! Extract primitive variables from VR (right state)
    rhoR = VR(1)  ! Density
    uR1  = VR(2)  ! u-Velocity
    vR1  = VR(3)  ! v-velocity
    pR   = VR(4)  ! Pressure
  
    ! Compute speed of sound and Mach number for right state
    aR   = sqrt(gamma * pR / rhoR)  ! Speed of sound
    ! MR   = uR1 / aR  ! Mach number
    unL= uL1*nx+vL1*ny
    unR = uR1*nx+vR1*ny


    ML = unL/aL
    MR = unR/aR
  
    ! Split Mach numbers
    m_plus  = ((ML + 1.0_prec)**2 )/ 4.0_prec
    m_minus = (-(MR - 1.0_prec)**2) / 4.0_prec
  
    ! Flow direction sensor
    alpha_plus  = (1.0_prec + sign(1.0_prec, ML)) / 2.0_prec
    alpha_minus = (1.0_prec - sign(1.0_prec, MR)) / 2.0_prec
  
    ! Supersonic/subsonic sensor
    beta_plus  = -real(max(0, 1 - int(abs(ML))), prec)
    beta_minus = -real(max(0, 1 - int(abs(MR))), prec)
  
    ! Convective flux coefficients
    C_plus  = (alpha_plus * (1.0_prec + beta_plus) * ML) - beta_plus * m_plus
    C_minus = (alpha_minus * (1.0_prec + beta_minus) * MR) - beta_minus * m_minus
  
    ! Pressure flux coefficients
    p_plus  = m_plus * (-ML + 2.0_prec)
    p_minus = m_minus * (-MR - 2.0_prec)
  
    D_plus  = alpha_plus * (1 + beta_plus) - beta_plus * p_plus
    D_minus = alpha_minus * (1 + beta_minus) - beta_minus * p_minus
  
    ! Compute total enthalpy
    hL = (pL / rhoL) * (gamma / (gamma - 1.0_prec)) + 0.5_prec * (uL1**2+vL1**2)
    hR = (pR / rhoR) * (gamma / (gamma - 1.0_prec)) + 0.5_prec * (uR1**2+vR1**2)
    ! hL = (aL**2) * (gamma / (gamma - 1.0_prec) + 0.5_prec * (uL1**2+vL1**2))
    ! hR = (aR**2)/ (gamma - 1.0_prec) + 0.5_prec * (uR1**2+vR1**2)
  
    ! Pressure flux
    Fp(1) = 0.0_prec
    Fp(2) = (D_plus * pL + D_minus * pR)*nx
    Fp(3) = (D_plus * pL + D_minus * pR)*ny
    Fp(4) = 0.0_prec

    ! Convective flux
    Fc(1) = (rhoL * aL * C_plus) + (rhoR * aR * C_minus)
    Fc(2) = (rhoL * aL * C_plus * uL1) + (rhoR * aR * C_minus * uR1)
    Fc(3) = (rhoL * aL * C_plus * vL1) + (rhoR * aR * C_minus * vR1)
    Fc(4) = (rhoL * aL * C_plus * hL) + (rhoR * aR * C_minus * hR)
  
    ! Total flux (in conservative form)
    F = Fc + Fp
  end subroutine vanleer_flux

 subroutine roe_flux(VL, VR, nx, ny, F)
    real(prec), dimension(neq), intent(in) :: VL, VR  ! Primitive variables: [rho, u, v, p]
    real(prec), dimension(neq), intent(out) :: F
    real(prec), intent(in) :: nx, ny
    real(prec) :: rhoL, uL1, pL, aL, hL, rhoR, uR1, pR, aR, hR, vL1, vR1
    real(prec) :: rho_avg, u_avg, h_avg, a_avg, Ri, unL, unR, v_avg, un_avg
    real(prec), dimension(neq) :: FL, FR, rvec1, rvec2, rvec3, rvec4, lambda
    real(prec) :: dw1, dw2, dw3, dw4, epsilon1
    integer :: i
  
    ! Extract primitive variables from VL (left state)
    rhoL = VL(1)  ! Density
    uL1  = VL(2)  ! x-velocity
    vL1  = VL(3)  ! y-velocity
    pL   = VL(4)  ! Pressure
    aL   = sqrt(gamma * pL / rhoL)  ! Speed of sound
    hL   = aL**2 / (gamma - one) + half * (uL1**2 + vL1**2)  ! Total enthalpy
  
    ! Extract primitive variables from VR (right state)
    rhoR = VR(1)  ! Density
    uR1  = VR(2)  ! x-velocity
    vR1  = VR(3)  ! y-velocity
    pR   = VR(4)  ! Pressure
    aR   = sqrt(gamma * pR / rhoR)  ! Speed of sound
    hR   = aR**2 / (gamma - one) + half * (uR1**2 + vR1**2)  ! Total enthalpy
    
    ! Calculate normal velocities
    unL = uL1*nx + vL1*ny
    unR = uR1*nx + vR1*ny
  
    ! Compute Roe averages
    Ri = sqrt(rhoR / rhoL)
    rho_avg = Ri * rhoL
    u_avg = (uL1 + Ri * uR1) / (one + Ri)
    v_avg = (vL1 + Ri * vR1) / (one + Ri)
    h_avg = (hL + Ri * hR) / (one + Ri)
    a_avg = sqrt((gamma - one) * (h_avg - half * (u_avg**2 + v_avg**2)))
    un_avg = u_avg*nx + v_avg*ny
  
    ! Define right eigenvectors
    rvec1 = [one, u_avg, v_avg, half*(u_avg**2 + v_avg**2)]
    rvec2 = [zero, ny*rho_avg, -nx*rho_avg, rho_avg*(ny*u_avg - nx*v_avg)]
    rvec3 = half*(rho_avg/a_avg)*[one, u_avg + nx*a_avg, v_avg + ny*a_avg, h_avg + un_avg*a_avg]
    rvec4 = -half*(rho_avg/a_avg)*[one, u_avg - nx*a_avg, v_avg - ny*a_avg, h_avg - un_avg*a_avg]
    
    ! Define eigenvalues
    lambda = [un_avg, un_avg, un_avg + a_avg, un_avg - a_avg]
    
    ! Apply entropy fix
    epsilon1 = 0.01_prec  ! Typically 0.1 or smaller
    do i = 1, 4
      lambda(i) = abs(lambda(i))
      if (lambda(i) < 2.0_prec * epsilon1 * a_avg) then
        lambda(i) = half * (lambda(i)**2 / (epsilon1 * a_avg) + epsilon1 * a_avg)
      end if
    end do
  
    ! Wave strengths (differences in primitive variables)
    dw1 = (rhoR - rhoL) - (pR - pL) / a_avg**2                  ! Entropy wave
    dw2 = ny*(uR1 - uL1) - nx*(vR1 - vL1)                       ! Shear wave
    dw3 = nx*(uR1 - uL1) + ny*(vR1 - vL1) + (pR - pL)/(rho_avg*a_avg)  ! Right acoustic wave
    dw4 = nx*(uR1 - uL1) + ny*(vR1 - vL1) - (pR - pL)/(rho_avg*a_avg)  ! Left acoustic wave
  
    ! Compute physical fluxes
    FL = [rhoL*unL, rhoL*uL1*unL + pL*nx, rhoL*vL1*unL + pL*ny, rhoL*hL*unL]
    FR = [rhoR*unR, rhoR*uR1*unR + pR*nx, rhoR*vR1*unR + pR*ny, rhoR*hR*unR]
  
    ! Roe flux: average flux minus dissipation term
    F = half * (FL + FR) - half * (lambda(1)*dw1*rvec1 + lambda(2)*dw2*rvec2 + &
                                   lambda(3)*dw3*rvec3 + lambda(4)*dw4*rvec4)
    
end subroutine roe_flux
    
end module flux