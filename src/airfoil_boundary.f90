module airfoil_boundary
    use set_precision, only: prec
    use set_constants, only: zero, one, two, half,pi
    use set_inputs, only: neq, i_low, i_high, j_low, j_high
    use set_inputs, only: ig_low, ig_high, jg_low, jg_high
    use soln_type, only: soln_t
    use grid_type, only: grid_t
    use variable_conversion, only: prim2cons, cons2prim, limit_primitives
    use fluid_constants, only: gamma
    
    implicit none
    private
    
    public :: initialize_airfoil_flow
    public :: apply_wake_cut_bc
    public :: apply_airfoil_slip_wall
    public :: apply_farfield_bc
    public :: apply_airfoil_outflow_bc
    ! public :: apply_all_airfoil_bc
    public :: set_airfoil_extent
    
    ! Farfield flow conditions
    real(prec), parameter :: P_inf = 65855.8_prec    ! Freestream pressure (Pa)
    real(prec), parameter :: T_inf = 300.0_prec      ! Freestream temperature (K)  
    real(prec), parameter :: M_inf = 0.84_prec         ! Freestream Mach number
    real(prec), parameter :: alpha_deg = 0.0_prec     ! Angle of attack (degrees)
    real(prec), parameter :: R_gas = 287.0_prec       ! Gas constant for air (J/kg-K)
    
    ! Airfoil surface configuration - will be set by user
    integer, save :: i_airfoil_start = 33   ! Default values
    integer, save :: i_airfoil_end = 160
    integer, save :: j_airfoil_start = 1
    integer, save :: j_airfoil_end = 1
    
contains

    ! Subroutine to set the airfoil surface extent
    subroutine set_airfoil_extent(i_start, i_end, j_start, j_end)
        integer, intent(in) :: i_start, i_end, j_start, j_end
        i_airfoil_start = i_start
        i_airfoil_end = i_end
        j_airfoil_start = j_start
        j_airfoil_end = j_end
    end subroutine set_airfoil_extent

    subroutine initialize_airfoil_flow(grid, soln)
        type(grid_t), intent(in) :: grid
        type(soln_t), intent(inout) :: soln
        
        integer :: i, j
        real(prec) :: rho_init, u_init, v_init, p_init
        real(prec) :: u_inf, v_inf, alpha_rad
        
        ! Convert angle of attack from degrees to radians
        alpha_rad = (alpha_deg * pi) / 180.0_prec
        
        ! Calculate freestream velocity components with angle of attack
        u_inf = M_inf * sqrt(gamma * R_gas * T_inf) * cos(alpha_rad)
        v_inf = M_inf * sqrt(gamma * R_gas * T_inf) * sin(alpha_rad)
        
        ! Initialize with freestream conditions everywhere
        rho_init = P_inf / (R_gas * T_inf)
        u_init = u_inf
        v_init = v_inf
        p_init = P_inf
        
        ! Set initial conditions in all cells
        do j = grid%j_cell_low, grid%j_cell_high
            do i = grid%i_cell_low, grid%i_cell_high 
                soln%V(1, i, j) = rho_init
                soln%V(2, i, j) = u_init
                soln%V(3, i, j) = v_init
                soln%V(4, i, j) = p_init
            end do
        end do
        
        call prim2cons(soln%U, soln%V)
    end subroutine initialize_airfoil_flow

subroutine apply_farfield_bc(grid, soln)
    type(grid_t), intent(in) :: grid
    type(soln_t), intent(inout) :: soln
    
    integer :: i, j
    real(prec) :: rho_inf, u_inf, v_inf, p_inf_local, alpha_rad
    
    ! Convert angle of attack from degrees to radians
    alpha_rad = (alpha_deg * pi) / 180.0_prec
    
    ! Calculate farfield conditions
    rho_inf = P_inf / (R_gas * T_inf)
    u_inf = M_inf * sqrt(gamma * R_gas * T_inf) * cos(alpha_rad)
    v_inf = M_inf * sqrt(gamma * R_gas * T_inf) * sin(alpha_rad)
    p_inf_local = P_inf
    
    ! Apply farfield BC to outer boundary (j = j_cell_high)
    do j = grid%jg_high-1, grid%jg_high  ! Top ghost cells
        do i = 1, grid%imax
            soln%V(1, i, j) = rho_inf
            soln%V(2, i, j) = u_inf
            soln%V(3, i, j) = v_inf
            soln%V(4, i, j) = p_inf_local
        end do
    end do
    
    call prim2cons(soln%U, soln%V)
    
end subroutine apply_farfield_bc

subroutine apply_airfoil_outflow_bc(grid, soln)
    type(grid_t), intent(in) :: grid
    type(soln_t), intent(inout) :: soln
    
    integer :: i, j, n
    
    ! First loop: Apply 2nd order extrapolation to ig_low to ig_low-1 ghost cells
    do i = grid%ig_low+1, grid%ig_low,-1 ! Left ghost cells
        do j = 1, grid%jmax
            do n = 1, neq
                ! 2nd order linear extrapolation: f(i) = 2*f(i+1) - f(i+2)
                soln%V(n, i, j) = two * soln%V(n, i+1, j) - soln%V(n, i+2, j)
            end do
        end do
    end do
    
    ! Second loop: Apply 2nd order extrapolation to ig_high-1 to ig_high ghost cells
    do i = grid%ig_high-1, grid%ig_high  ! Right ghost cells
        do j = 1, grid%jmax
            do n = 1, neq
                soln%V(n, i, j) = two * soln%V(n, i-1, j) - soln%V(n, i-2, j)
            end do
        end do
    end do
    
    call limit_primitives(soln%V)
    call prim2cons(soln%U, soln%V)
    
end subroutine apply_airfoil_outflow_bc

subroutine apply_airfoil_slip_wall(grid, soln)
    type(grid_t), intent(in) :: grid
    type(soln_t), intent(inout) :: soln
    
    integer :: i
    real(prec) :: u_interior, v_interior, u_reflected, v_reflected
    real(prec) :: nx, ny, dot_product
    real(prec) :: p_ghost1, p_ghost2, T_interior
    
    ! Apply slip wall BC to airfoil surface at j=1
    do i = i_airfoil_start, i_airfoil_end
        
        ! Get wall normal vector at j=1
        nx = grid%n_eta(i, 1, 1)  
        ny = grid%n_eta(i, 1, 2)
        
        ! First ghost cell (jg_low-1 = 0) - get interior velocity from 1st interior cell
        u_interior = soln%V(2, i, 1)
        v_interior = soln%V(3, i, 1)
        
        ! Reflect velocity: v_reflected = v_interior - 2*(v_interior·n)*n
        dot_product = u_interior * nx + v_interior * ny
        u_reflected = u_interior - two * dot_product * nx
        v_reflected = v_interior - two * dot_product * ny
        
        ! Set first ghost cell values
        soln%V(2, i, grid%jg_low+1) = u_reflected  ! u-velocity
        soln%V(3, i, grid%jg_low+1) = v_reflected  ! v-velocity
        
        ! Copy temperature from 1st interior cell
        T_interior = soln%V(4, i, 1) / (R_gas * soln%V(1, i, 1))
        
        ! Pressure: 2nd order extrapolation
        p_ghost1 = two * soln%V(4, i, 1) - soln%V(4, i, 2)
        soln%V(4, i, grid%jg_low+1) = p_ghost1
        
        ! Density from ideal gas law
        soln%V(1, i, grid%jg_low+1) = p_ghost1 / (R_gas * T_interior)
        
        ! Second ghost cell (jg_low = 0) - get interior velocity from 2nd interior cell
        u_interior = soln%V(2, i, 2)
        v_interior = soln%V(3, i, 2)
        
        ! Reflect velocity
        dot_product = u_interior * nx + v_interior * ny
        u_reflected = u_interior - two * dot_product * nx
        v_reflected = v_interior - two * dot_product * ny
        
        ! Set second ghost cell values
        soln%V(2, i, grid%jg_low) = u_reflected
        soln%V(3, i, grid%jg_low) = v_reflected
        
        ! Copy temperature from 2nd interior cell
        T_interior = soln%V(4, i, 2) / (R_gas * soln%V(1, i, 2))
        
        ! Pressure: 2nd order extrapolation
        p_ghost2 = two * p_ghost1 - soln%V(4, i, 1)
        soln%V(4, i, grid%jg_low) = p_ghost2
        
        ! Density from ideal gas law
        soln%V(1, i, grid%jg_low) = p_ghost2 / (R_gas * T_interior)
        
    end do
    
    ! call limit_primitives(soln%V)
    call prim2cons(soln%U, soln%V)
    
end subroutine apply_airfoil_slip_wall

subroutine apply_wake_cut_bc(grid, soln)
    type(grid_t), intent(in) :: grid
    type(soln_t), intent(inout) :: soln

    integer :: i, n
    integer :: i2

    ! First segment: i = 1 to i_airfoil_start-1
    i2 = grid%imax+1
    do i = 1, i_airfoil_start-1
        i2 = i2 - 1
        do n = 1, neq
            soln%V(n, i, grid%jg_low+1) = soln%V(n, i2, 1)
            soln%V(n, i, grid%jg_low)   = soln%V(n, i2, 2)
        end do
    end do

    ! Second segment: i = i_airfoil_end+1 to imax
    i2 = i_airfoil_start
    do i = i_airfoil_end+1, grid%imax
        i2 = i2 - 1
        !  i_airfoil_start, 1, -1
        do n = 1, neq
            soln%V(n, i, grid%jg_low+1) = soln%V(n, i2, 1)
            soln%V(n, i, grid%jg_low)   = soln%V(n, i2, 2)
        end do
    end do

    call prim2cons(soln%U, soln%V)
end subroutine apply_wake_cut_bc


end module airfoil_boundary