module flux
    use set_precision
    use fluid_constants
    use set_inputs
    use set_constants
    use grid_type
    use soln_type
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
    real(prec), dimension(grid%ig_low:grid%ig_high,grid%jg_low:grid%jg_high,neq), intent(in) :: U
    real(prec), dimension(neq,grid%i_low:grid%i_high,grid%j_low:grid%j_high), intent(out) :: Lxi, Rxi
    real(prec), dimension(neq,grid%i_low:grid%i_high,grid%j_low:grid%j_high), intent(out) :: Leta, Reta
    real(prec), dimension(grid%ig_low:grid%ig_high,grid%jg_low:grid%jg_high,neq) :: V
    real(prec), dimension(grid%ig_low:grid%ig_high,grid%jg_low:grid%jg_high,neq) :: psi_plus_xi, psi_minus_xi
    real(prec), dimension(grid%ig_low:grid%ig_high,grid%jg_low:grid%jg_high,neq) :: psi_plus_eta, psi_minus_eta
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
            V(i,j,1) = U(i,j,1)  ! rho
            V(i,j,2) = U(i,j,2) / U(i,j,1)  ! u
            V(i,j,3) = U(i,j,3) / U(i,j,1)  ! v
            V(i,j,4) = (gamma - one) * (U(i,j,4) - half * (U(i,j,2)**2 + U(i,j,3)**2) / U(i,j,1))  ! p
        end do
    end do
    ! call limit_primitives(U, V)

    ! Initialize limiter arrays
    psi_plus_xi = one
    psi_minus_xi = one
    psi_plus_eta = one
    psi_minus_eta = one

    ! Compute limiters for second-order
    if (order == 2) then
        ! call select_limiter(limiter_scheme)
        ! call calculate_limiters(V, psi_plus_xi, psi_minus_xi)  ! ξ-direction
        ! call calculate_limiters(V, psi_plus_eta, psi_minus_eta, direction='eta')  ! η-direction
    end if

    ! xsi-direction extrapolation (faces at i_low:i_high)
    do j = grid%j_low, grid%j_high-1
        do i = grid%i_low, grid%i_high
            Lxi(:,i,j) = V(i-1,j,:) + fourth * epsilon * ( &
                (one - kappa) * psi_plus_xi(i-2,j,:) * (V(i-1,j,:) - V(i-2,j,:)) + &
                (one + kappa) * psi_minus_xi(i-1,j,:) * (V(i,j,:) - V(i-1,j,:)))
            Rxi(:,i,j) = V(i,j,:) - fourth * epsilon * ( &
                (one - kappa) * psi_plus_xi(i,j,:) * (V(i+1,j,:) - V(i,j,:)) + &
                (one + kappa) * psi_minus_xi(i-1,j,:) * (V(i,j,:) - V(i-1,j,:)))
        end do
    end do

    ! etadirection extrapolation (faces at j_low:j_high)
    do j = grid%j_low, grid%j_high
        do i = grid%i_low, grid%i_high-1
            Leta(:,i,j) = V(i,j-1,:) + fourth * epsilon * ( &
                (one - kappa) * psi_plus_eta(i,j-2,:) * (V(i,j-1,:) - V(i,j-2,:)) + &
                (one + kappa) * psi_minus_eta(i,j-1,:) * (V(i,j,:) - V(i,j-1,:)))
            Reta(:,i,j) = V(i,j,:) - fourth * epsilon * ( &
                (one - kappa) * psi_plus_eta(i,j,:) * (V(i,j+1,:) - V(i,j,:)) + &
                (one + kappa) * psi_minus_eta(i,j-1,:) * (V(i,j,:) - V(i,j-1,:)))
        end do
    end do
end subroutine compute_lr_states

subroutine compute_fluxes(grid, soln)
    type(grid_t), intent(in)    :: grid
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq,grid%i_low:grid%i_high,grid%j_low:grid%j_high) :: Lxi, Rxi
    real(prec), dimension(neq,grid%i_low:grid%i_high,grid%j_low:grid%j_high) :: Leta, Reta
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
            call flux_fun( Leta(:,i,j), Reta(:,i,j), nx, ny, soln%Feta(:,i,j) ) 
        end do
    end do
end subroutine compute_fluxes

! subroutine calculate_fluxes_normals(grid,soln)

!     type(grid_t), intent(in)    :: grid
!     type(soln_t), intent(inout) :: soln
!     real(prec), dimension(neq,grid%i_low:grid%i_high,grid%j_low:grid%j_high) :: Lxi, Rxi
!     real(prec), dimension(neq,grid%i_low:grid%i_high,grid%j_low:grid%j_high) :: Leta, Reta
!     real(prec) :: nx, ny
!     integer :: i,j,k

!     do j = grid%j_low,grid%j_high-1
!         do i = grid%i_low,grid%i_high
!           nx = grid%n_xi(i,j,1)
!           ny = grid%n_xi(i,j,2)
!         end do
!     end do

!     do j = grid%j_low,grid%j_high
!         do i = grid%i_low,grid%i_high-1
!           nx = grid%n_eta(i,j,1)
!           ny = grid%n_eta(i,j,2)
!         end do
!     end do


! end subroutine calculate_fluxes_normals


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
    m_plus  = (ML + 1.0_prec)**2 / 4.0_prec
    m_minus = -(MR - 1.0_prec)**2 / 4.0_prec
  
    ! Flow direction sensor
    alpha_plus  = (1.0_prec + sign(1.0_prec, ML)) / 2.0_prec
    alpha_minus = (1.0_prec - sign(1.0_prec, MR)) / 2.0_prec
  
    ! Supersonic/subsonic sensor
    beta_plus  = -real(max(0, 1 - int(abs(ML))), prec)
    beta_minus = -real(max(0, 1 - int(abs(MR))), prec)
  
    ! Convective flux coefficients
    C_plus  = alpha_plus * (1.0_prec + beta_plus) * ML - beta_plus * m_plus
    C_minus = alpha_minus * (1.0_prec + beta_minus) * MR - beta_minus * m_minus
  
    ! Pressure flux coefficients
    p_plus  = m_plus * (-ML + 2.0_prec)
    p_minus = m_minus * (-MR - 2.0_prec)
  
    D_plus  = alpha_plus * (1 + beta_plus) - beta_plus * p_plus
    D_minus = alpha_minus * (1 + beta_minus) - beta_minus * p_minus
  
    ! Compute total enthalpy
    hL = (pL / rhoL) * (gamma / (gamma - 1.0_prec)) + 0.5_prec * (uL1**2+vL1**2)
    hR = (pR / rhoR) * (gamma / (gamma - 1.0_prec)) + 0.5_prec * (uR1**2+vR1**2)
  
    ! Pressure flux
    Fp(1) = 0.0_prec
    Fp(2) = D_plus * pL + D_minus * pR
    Fp(3) = 0.0_prec

    ! Convective flux
    Fc(1) = (rhoL * aL * C_plus) + (rhoR * aR * C_minus)
    Fc(2) = (rhoL * aL * C_plus * uL1) + (rhoR * aR * C_minus * uR1) + Fp(2)*nx
    Fc(3) = (rhoL * aL * C_plus * vL1) + (rhoR * aR * C_minus * vR1) + Fp(2)*ny
    Fc(4) = (rhoL * aL * C_plus * hL) + (rhoR * aR * C_minus * hR)
  
    ! Total flux (in conservative form)
    F = Fc + Fp
  end subroutine vanleer_flux

  subroutine roe_flux(VL, VR,nx,ny,F)
    real(prec), dimension(neq), intent(in) :: VL, VR  ! Primitive variables: [rho, u, p]
    real(prec), dimension(neq), intent(out) :: F
    real(prec),intent(in):: nx, ny
    real(prec) :: rhoL, uL1, pL, aL, hL, rhoR, uR1, pR, aR, hR
    real(prec) :: rho_avg, u_avg, h_avg, a_avg, Ri
    real(prec), dimension(neq) :: FL, FR, rvec1, rvec2, rvec3, lambda
    real(prec) :: dw1, dw2, dw3,epsilon1
    integer :: i
  
    ! Extract primitive variables from VL (left state)
    rhoL = VL(1)  ! Density
    uL1  = VL(2)  ! Velocity
    pL   = VL(3)  ! Pressure
    aL   = sqrt(gamma * pL / rhoL)  ! Speed of sound
    ! hL   = aL**2 / (gamma - one) + half * uL1**2  ! Total enthalpy
  
    ! Extract primitive variables from VR (right state)
    rhoR = VR(1)  ! Density
    uR1  = VR(2)  ! Velocity
    pR   = VR(3)  ! Pressure
    aR   = sqrt(gamma * pR / rhoR)  ! Speed of sound
    ! hR   = aR**2 / (gamma - one) + half * uR1**2  ! Total enthalpy
    hL = (pL / rhoL) * (gamma / (gamma - 1.0_prec)) + 0.5_prec * uL1**2
    hR = (pR / rhoR) * (gamma / (gamma - 1.0_prec)) + 0.5_prec * uR1**2
  
    ! Compute Roe averages
    Ri = sqrt(rhoR / rhoL)
    rho_avg = Ri * rhoL
    u_avg   = (Ri * uR1 + uL1) / (Ri+ one)
    h_avg   = (Ri * hR + hL) / (Ri + one)
    a_avg   = sqrt((gamma - one) * (h_avg - half * u_avg**2))
  
    ! Define right eigenvectors
    rvec1 = [one, u_avg, half * u_avg**2]  ! Acoustic wave (u)
    rvec2 = (half*(rho_avg/a_avg)*([one, u_avg + a_avg, h_avg + u_avg * a_avg]))  ! Right-moving wave (u + a)
    rvec3 = (half*(-rho_avg/a_avg)*([one, u_avg - a_avg, h_avg - u_avg * a_avg]))  ! Left-moving wave (u - a)

    lambda = [(u_avg), (u_avg + a_avg), (u_avg - a_avg)]
    
    ! lambda_old = lambda
    
    epsilon1 = 0.05_prec  ! Typically 0.1 or smaller
    do i = 1, 3
      if (abs(lambda(i)) < 2.0_prec * epsilon1 * a_avg) then
        lambda(i) = (abs(lambda(i)**2)) / ((4.0_prec * epsilon1 * a_avg) + epsilon1 * a_avg)
      end if
    end do
  
    ! Wave strengths (differences in primitive variables)
    dw1 = (rhoR - rhoL) - ((pR - pL) / a_avg**2)  ! Entropy wave
    dw2 = (uR1 - uL1)+(pR - pL) / (rho_avg * a_avg)   ! Right-moving acoustic wave
    dw3 = (uR1 - uL1)-((pR - pL) / (rho_avg * a_avg))  ! Left-moving acoustic wave
  
    ! Compute fluxes in conservative form from primitive variables
    FL = [rhoL * uL1, rhoL * uL1**2 + pL, rhoL * hL * uL1]  ! Left flux
    FR = [rhoR * uR1, rhoR * uR1**2 + pR, rhoR * hR * uR1]  ! Right flux

  
    ! Roe flux: average flux minus dissipation term
    F = half * (FL + FR) - (half *abs(lambda(1)) * dw1 * rvec1 + &
                                   abs(lambda(2)) * dw2 * rvec2 + &
                                   abs(lambda(3)) * dw3 * rvec3)
    
  end subroutine roe_flux
    
end module flux