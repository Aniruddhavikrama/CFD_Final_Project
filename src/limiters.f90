module limiters
    use set_precision, only: prec
    use set_constants, only: zero, one, two, half, fourth
    use set_inputs
    use grid_type, only: grid_t
    implicit none
    private
    public :: calculate_limiters, van_leer_limiter, minmod_limiter, select_limiter
    public :: init_limiters, finalize_limiters, venkatakrishnan_limiter, beta_limiter
  
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
            allocate(frozen_psi_p_xi(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high))
            allocate(frozen_psi_m_xi(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high))
            allocate(frozen_psi_p_eta(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high))
            allocate(frozen_psi_m_eta(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high))
        end if
        
        ! Initialize to unity (no limiting)
        frozen_psi_p_xi = one
        frozen_psi_m_xi = one
        frozen_psi_p_eta = one
        frozen_psi_m_eta = one
        
        limiters_frozen = .false.
    end subroutine init_limiters

    subroutine finalize_limiters()
        if (allocated(frozen_psi_p_xi)) deallocate(frozen_psi_p_xi)
        if (allocated(frozen_psi_m_xi)) deallocate(frozen_psi_m_xi)
        if (allocated(frozen_psi_p_eta)) deallocate(frozen_psi_p_eta)
        if (allocated(frozen_psi_m_eta)) deallocate(frozen_psi_m_eta)
        limiters_frozen = .false.
    end subroutine finalize_limiters

    subroutine select_limiter(scheme)
        integer, intent(in) :: scheme
        
        limiter_fun => null()
        select case (scheme)
        case (0)
            limiter_fun => no_limiter
        case (1)
            limiter_fun => minmod_limiter
        case (2)
            limiter_fun => van_leer_limiter
        case (3)
            limiter_fun => venkatakrishnan_limiter
        case (4)
            limiter_fun => beta_limiter
        case default
            write(*,*) 'ERROR: Invalid limiter_scheme value:', scheme
            stop
        end select
    end subroutine select_limiter

    subroutine calculate_limiters(grid, V, psi_plus, psi_minus, direction)
        type(grid_t), intent(in) :: grid
        real(prec), dimension(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high), intent(in) :: V
        real(prec), dimension(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high), intent(out) :: psi_plus, psi_minus
        character(len=*), intent(in), optional :: direction
        
        real(prec), dimension(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high) :: r_plus, r_minus
        character(len=10) :: dir
        integer :: i, j, k
        
        ! Default direction is xi
        dir = 'xi'
        if (present(direction)) dir = trim(direction)
        
        ! If limiters are frozen, use stored values
        if (limiters_frozen) then
            select case (trim(dir))
            case ('xi')
                psi_plus = frozen_psi_p_xi
                psi_minus = frozen_psi_m_xi
            case ('eta')
                psi_plus = frozen_psi_p_eta
                psi_minus = frozen_psi_m_eta
            end select
            return
        end if
        
        ! Initialize
        psi_plus = one
        psi_minus = one
        r_plus = zero
        r_minus = zero
        
        ! Compute smoothness indicators
        select case (trim(dir))
        case ('xi')
            call compute_smoothness_xi(grid, V, r_plus, r_minus)
        case ('eta')
            call compute_smoothness_eta(grid, V, r_plus, r_minus)
        end select
        
        ! Apply limiter function
        if (associated(limiter_fun)) then
            do k = 1, neq
                call limiter_fun(r_plus(k,:,:), psi_plus(k,:,:))
                call limiter_fun(r_minus(k,:,:), psi_minus(k,:,:))
            end do
        end if
        
        ! Store frozen values if needed
        if (limiters_frozen) then
            select case (trim(dir))
            case ('xi')
                frozen_psi_p_xi = psi_plus
                frozen_psi_m_xi = psi_minus
            case ('eta')
                frozen_psi_p_eta = psi_plus
                frozen_psi_m_eta = psi_minus
            end select
        end if
    end subroutine calculate_limiters

    subroutine compute_smoothness_xi(grid, V, r_plus, r_minus)
        type(grid_t), intent(in) :: grid
        real(prec), dimension(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high), intent(in) :: V
        real(prec), dimension(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high), intent(out) :: r_plus, r_minus
        
        integer :: i, j, k
        real(prec) :: eps = 1.0e-12_prec
        
        do k = 1, neq
            do j = grid%jg_low, grid%jg_high
                do i = grid%ig_low+1, grid%ig_high-1
                    ! Forward smoothness indicator
                    if (abs(V(k,i+1,j) - V(k,i,j)) > eps) then
                        r_plus(k,i,j) = (V(k,i,j) - V(k,i-1,j)) / (V(k,i+1,j) - V(k,i,j))
                    else
                        r_plus(k,i,j) = one
                    end if
                    
                    ! Backward smoothness indicator
                    if (abs(V(k,i,j) - V(k,i-1,j)) > eps) then
                        r_minus(k,i,j) = (V(k,i+1,j) - V(k,i,j)) / (V(k,i,j) - V(k,i-1,j))
                    else
                        r_minus(k,i,j) = one
                    end if
                end do
            end do
        end do
    end subroutine compute_smoothness_xi

    subroutine compute_smoothness_eta(grid, V, r_plus, r_minus)
        type(grid_t), intent(in) :: grid
        real(prec), dimension(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high), intent(in) :: V
        real(prec), dimension(neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high), intent(out) :: r_plus, r_minus
        
        integer :: i, j, k
        real(prec) :: eps = 1.0e-12_prec
        
        do k = 1, neq
            do j = grid%jg_low+1, grid%jg_high-1
                do i = grid%ig_low, grid%ig_high
                    ! Forward smoothness indicator
                    if (abs(V(k,i,j+1) - V(k,i,j)) > eps) then
                        r_plus(k,i,j) = (V(k,i,j) - V(k,i,j-1)) / (V(k,i,j+1) - V(k,i,j))
                    else
                        r_plus(k,i,j) = one
                    end if
                    
                    ! Backward smoothness indicator
                    if (abs(V(k,i,j) - V(k,i,j-1)) > eps) then
                        r_minus(k,i,j) = (V(k,i,j+1) - V(k,i,j)) / (V(k,i,j) - V(k,i,j-1))
                    else
                        r_minus(k,i,j) = one
                    end if
                end do
            end do
        end do
    end subroutine compute_smoothness_eta

    subroutine no_limiter(r, psi)
        real(prec), dimension(:,:), intent(in) :: r
        real(prec), dimension(:,:), intent(out) :: psi
        psi = one
    end subroutine no_limiter

    subroutine minmod_limiter(r, psi)
        real(prec), dimension(:,:), intent(in) :: r
        real(prec), dimension(:,:), intent(out) :: psi
        
        psi = max(zero, min(one, r))
    end subroutine minmod_limiter

    subroutine van_leer_limiter(r, psi)
        real(prec), dimension(:,:), intent(in) :: r
        real(prec), dimension(:,:), intent(out) :: psi
        
        where (r > zero)
            psi = (r + abs(r)) / (one + abs(r))
        elsewhere
            psi = zero
        end where
    end subroutine van_leer_limiter

    subroutine venkatakrishnan_limiter(r, psi)
        real(prec), dimension(:,:), intent(in) :: r
        real(prec), dimension(:,:), intent(out) :: psi
        
        real(prec), parameter :: k_param = 0.3_prec
        real(prec) :: r2, k2
        integer :: i, j
        
        k2 = k_param**2
        
        do j = lbound(r,2), ubound(r,2)
            do i = lbound(r,1), ubound(r,1)
                r2 = r(i,j)**2
                if (r(i,j) > zero) then
                    psi(i,j) = (r2 + two*r(i,j) + k2) / (r2 + r(i,j) + two*k2)
                else
                    psi(i,j) = zero
                end if
            end do
        end do
    end subroutine venkatakrishnan_limiter

    subroutine beta_limiter(r, psi)
        real(prec), dimension(:,:), intent(in) :: r
        real(prec), dimension(:,:), intent(out) :: psi
        
        real(prec), parameter :: beta = 2.0_prec
        
        where (r > zero)
            psi = max(zero, min(beta*r, one), min(r, beta))
        elsewhere
            psi = zero
        end where
    end subroutine beta_limiter

end module limiters