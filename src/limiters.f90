module limiters
    use set_precision, only: prec
    use set_constants, only: zero, one, two, half, fourth
    use set_inputs
    use grid_type, only: grid_t
    implicit none
    private
    public :: calculate_limiters, van_leer_limiter, minmod_limiter, select_limiter
    public :: init_limiters, finalize_limiters,venkatakrishnan_limiter,beta_limiter
  
    procedure(calc_limiter), pointer :: limiter_fun => null()
    logical :: limiters_frozen = .false.
    real(prec), allocatable, dimension(:,:,:) :: frozen_psi_p_xi, frozen_psi_m_xi
    real(prec), allocatable, dimension(:,:,:) :: frozen_psi_p_eta, frozen_psi_m_eta

    abstract interface
        subroutine calc_limiter(r, psi)
            import :: prec
            real(prec), dimension(:,:), intent(in) :: r
            real(prec), dimension(:,:), intent(out) :: psi
        end subroutine calc_limiter
    end interface
  
contains
    subroutine init_limiters(grid)
        type(grid_t), intent(in) :: grid
        if (.not. allocated(frozen_psi_p_xi)) then
            allocate(frozen_psi_p_xi(neq, grid%i_low:grid%i_high, grid%j_low:grid%j_high))
            allocate(frozen_psi_m_xi(neq, grid%i_low:grid%i_high, grid%j_low:grid%j_high))
            allocate(frozen_psi_p_eta(neq, grid%i_low:grid%i_high, grid%j_low:grid%j_high))
            allocate(frozen_psi_m_eta(neq, grid%i_low:grid%i_high, grid%j_low:grid%j_high))
            frozen_psi_p_xi = one
            frozen_psi_m_xi = one
            frozen_psi_p_eta = one
            frozen_psi_m_eta = one
        end if
    end subroutine init_limiters

    subroutine finalize_limiters
        if (allocated(frozen_psi_p_xi)) then
            deallocate(frozen_psi_p_xi, frozen_psi_m_xi, frozen_psi_p_eta, frozen_psi_m_eta)
        end if
    end subroutine finalize_limiters
  
    subroutine select_limiter(limiter_scheme)
        integer, intent(in) :: limiter_scheme
        select case (limiter_scheme)
        case (0)
            limiter_fun => null_limiter
        case (1)
            limiter_fun => van_leer_limiter
        case (2)
            limiter_fun => minmod_limiter
        case (3)
            limiter_fun => venkatakrishnan_limiter
        case (4)
            limiter_fun => beta_limiter
        case default
            write(*,*) 'ERROR: Invalid limiter_scheme value. Use 0 (none), 1 (Van Leer), or 2 (Minmod)'
            stop
        end select
    end subroutine select_limiter
  
    subroutine calculate_limiters(grid, V, psi_plus, psi_minus, direction)
        type(grid_t), intent(in) :: grid
        real(prec), dimension(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high), intent(in) :: V
        real(prec), dimension(neq, grid%i_low:grid%i_high, grid%j_low:grid%j_high), intent(out) :: psi_plus, psi_minus
        character(len=*), intent(in), optional :: direction
        real(prec), dimension(neq, grid%i_low:grid%i_high, grid%j_low:grid%j_high) :: den_xi, r_plus_xi, r_minus_xi
        real(prec), dimension(neq, grid%i_low:grid%i_high, grid%j_low:grid%j_high) :: den_eta, r_plus_eta, r_minus_eta

        integer :: i, j
    
        if (limiters_frozen) then
            if (present(direction) .and. direction == 'eta') then
                psi_plus = frozen_psi_p_eta
                psi_minus = frozen_psi_m_eta
            else
                psi_plus = frozen_psi_p_xi
                psi_minus = frozen_psi_m_xi
            end if
            return
        end if

        psi_plus = one
        psi_minus = one
    

            do j = grid%j_low, grid%j_high
                do i = grid%i_low, grid%i_high-1
                    den_eta(:,i,j) = V(:,i,j) - V(:,i,j-1)
                    den_eta(:,i,j) = sign(one, den_eta(:,i,j)) * max(abs(den_eta(:,i,j)), 1.0e-6_prec)
                    r_plus_eta(:,i,j) = (V(:,i,j+1) - V(:,i,j)) / den_eta(:,i,j)
                    r_minus_eta(:,i,j) = (V(:,i,j-1) - V(:,i,j-2)) / den_eta(:,i,j)
                end do
            end do
            do j = grid%j_low, grid%j_high
                call limiter_fun(r_plus_eta(:,grid%i_low:grid%i_high-1,j), psi_plus(:,grid%i_low:grid%i_high-1,j))
                call limiter_fun(r_minus_eta(:,grid%i_low:grid%i_high-1,j), psi_minus(:,grid%i_low:grid%i_high-1,j))
            end do
            frozen_psi_p_eta = psi_plus
            frozen_psi_m_eta = psi_minus
        
            ! xsi-direction limiting (faces at i_low:i_high)
            do j = grid%j_low, grid%j_high-1
                do i = grid%i_low, grid%i_high
                    den_xi(:,i,j) = V(:,i,j) - V(:,i-1,j)
                    den_xi(:,i,j) = sign(one, den_xi(:,i,j)) * max(abs(den_xi(:,i,j)), 1.0e-6_prec)
                    r_plus_xi(:,i,j) = (V(:,i+1,j) - V(:,i,j)) / den_xi(:,i,j)
                    r_minus_xi(:,i,j) = (V(:,i-1,j) - V(:,i-2,j)) / den_xi(:,i,j)
                end do
            end do
            do j = grid%j_low, grid%j_high-1
                call limiter_fun(r_plus_xi(:,grid%i_low:grid%i_high,j), psi_plus(:,grid%i_low:grid%i_high,j))
                call limiter_fun(r_minus_xi(:,grid%i_low:grid%i_high,j), psi_minus(:,grid%i_low:grid%i_high,j))
            end do
            frozen_psi_p_xi = psi_plus
            frozen_psi_m_xi = psi_minus
     
    end subroutine calculate_limiters
  
    subroutine null_limiter(r, psi)
        real(prec), dimension(:,:), intent(in) :: r
        real(prec), dimension(:,:), intent(out) :: psi
        psi = one
    end subroutine null_limiter


    subroutine venkatakrishnan_limiter(r, psi)
    real(prec), dimension(:,:), intent(in) :: r
    real(prec), dimension(:,:), intent(out) :: psi
    real(prec), parameter :: K = 0.3_prec  ! Venkatakrishnan parameter
    real(prec) :: K_sq
    
    K_sq = K * K
    where (r > zero)
        psi = (r * r + two * r + K_sq) / (r * r + r + two + K_sq)
    elsewhere
        psi = zero
    end where
    end subroutine venkatakrishnan_limiter

    subroutine beta_limiter(r, psi)
        real(prec), dimension(:,:), intent(in) :: r
        real(prec), dimension(:,:), intent(out) :: psi
        real(prec), parameter :: beta = 2.0_prec  ! Beta parameter (1 ≤ β ≤ 2)
        
        where (r > zero)
            psi = max(zero, min(beta * r, one, min(r, beta)))
        elsewhere
            psi = zero
        end where
    end subroutine beta_limiter
  
    subroutine van_leer_limiter(r, psi)
        real(prec), dimension(:,:), intent(in) :: r
        real(prec), dimension(:,:), intent(out) :: psi
        psi = sign(one, one + r) * max(abs(one + r), 1.0e-12_prec)
        psi = (r + abs(r)) / psi
        psi = half * (one + sign(one, r)) * psi
    end subroutine van_leer_limiter
  
    subroutine minmod_limiter(r, psi)
        real(prec), dimension(:,:), intent(in) :: r
        real(prec), dimension(:,:), intent(out) :: psi
        psi = half * (one + sign(one, r)) * min(r, one)
    end subroutine minmod_limiter
end module limiters