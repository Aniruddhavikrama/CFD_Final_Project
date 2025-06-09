module inlet_boundary
    use set_precision, only: prec
    use set_constants, only: zero, one, two, half
    use set_inputs, only: neq, i_low, i_high, j_low, j_high
    use set_inputs, only: ig_low, ig_high, jg_low, jg_high
    use soln_type, only: soln_t
    use grid_type, only: grid_t
    use variable_conversion, only: prim2cons, cons2prim,limit_primitives
    use fluid_constants, only: gamma
    
    implicit none
    private
    
    public :: initialize_flow
    public :: apply_inlet_bc
    public :: apply_slip_walls
    public :: apply_outflow_bc
    public :: set_inlet_extent
    
    ! Inlet flow conditions
    real(prec), parameter :: P_inf = 12270.0_prec   ! Freestream pressure (Pa) - from table
    real(prec), parameter :: T_inf = 217.0_prec     ! Freestream temperature (K) - from table  
    real(prec), parameter :: M_inf = 4.0_prec       ! Freestream Mach number - from table
    real(prec), parameter :: R_gas = 287.0_prec     ! Gas constant for air (J/kg-K)
    real(prec), parameter :: y_inlet_max = 0.6_prec ! Maximum y-coordinate for inlet
    
    ! Boundary configuration - will be set by user
    integer, save :: i_inlet_max = 20  ! Default value, can be changed by user
    
contains

    ! Subroutine to set the inlet extent
    subroutine set_inlet_extent(inlet_max)
        integer, intent(in) :: inlet_max
        i_inlet_max = inlet_max
    end subroutine set_inlet_extent

    subroutine initialize_flow(grid, soln)
        type(grid_t), intent(in) :: grid
        type(soln_t), intent(inout) :: soln
        
        integer :: i, j
        real(prec) :: rho_init, u_init, v_init, p_init
        
        ! Initialize with freestream conditions everywhere
        rho_init = P_inf / (R_gas * T_inf)
        u_init = M_inf * sqrt(gamma * R_gas * T_inf)
        v_init = zero
        p_init = P_inf
        
        ! Set initial conditions in all cells
        do j = grid%j_cell_low, grid%j_cell_high
            do i = grid%i_cell_low, grid%i_cell_high 
                soln%V(1, i, j) = rho_init
                soln%V(2, i, j) = u_init
                soln%V(3, i, j) = v_init
                soln%V(4, i, j) = p_init
                
                call prim2cons(soln%U, soln%V)
            end do
        end do
    end subroutine initialize_flow
    
    subroutine apply_inlet_bc(grid, soln)
        type(grid_t), intent(in) :: grid
        type(soln_t), intent(inout) :: soln
        
        integer :: i, j
        real(prec) :: rho_inlet, u_inlet, v_inlet, p_inlet
        
        ! Calculate inlet conditions
        rho_inlet = P_inf / (R_gas * T_inf)
        u_inlet = M_inf * sqrt(gamma * R_gas * T_inf)
        v_inlet = zero
        p_inlet = P_inf
        
        ! Apply inlet BC from i=1 to i_inlet_max at top boundary
        do j = grid%jg_high-1, grid%jg_high  ! Top ghost cells
            do i = 1, i_inlet_max
                soln%V(1, i, j) = rho_inlet
                soln%V(2, i, j) = u_inlet
                soln%V(3, i, j) = v_inlet
                soln%V(4, i, j) = p_inlet
            end do
        end do
        
        call prim2cons(soln%U, soln%V)
        
    end subroutine apply_inlet_bc
    
    subroutine apply_slip_walls(grid, soln)
        type(grid_t), intent(in) :: grid
        type(soln_t), intent(inout) :: soln
        
        integer :: i, j
        real(prec) :: u_interior, v_interior, u_reflected, v_reflected
        real(prec) :: nx, ny, dot_product, u_normal, v_normal, u_tangent, v_tangent
        real(prec) :: p_ghost1, p_ghost2, T_interior, rho_ghost
        
        ! Bottom wall: from i=1 to grid%imax - 2 ghost cells below
        do i = 1, grid%imax
            ! Get wall normal vector (pointing into domain, upward for bottom wall)
            nx = grid%n_eta(i, 1, 1)  ! Normal x-component at bottom face
            ny = grid%n_eta(i, 1, 2)  ! Normal y-component at bottom face
            
            ! First ghost cell (j = grid%jg_low)
            ! Get interior velocity from j = 1
            u_interior = soln%V(2, i, 1)
            v_interior = soln%V(3, i, 1)
            
            ! Reflect velocity: v_reflected = v_interior - 2*(v_interior·n)*n
            dot_product = u_interior * nx + v_interior * ny
            u_reflected = u_interior - two * dot_product * nx
            v_reflected = v_interior - two * dot_product * ny
            
            ! Set ghost cell values
            soln%V(2, i, grid%jg_low+1) = u_reflected  ! u-velocity
            soln%V(3, i, grid%jg_low+1) = v_reflected  ! v-velocity
            
            ! Copy temperature (via ideal gas law with copied density calculation)
            T_interior = soln%V(4, i, 1) / (R_gas * soln%V(1, i, 1))  ! T = p/(rho*R)
            
            ! Pressure: 2nd order extrapolation
            p_ghost1 = two * soln%V(4, i, 1) - soln%V(4, i, 2)
            soln%V(4, i, grid%jg_low+1) = p_ghost1
            
            ! Density from ideal gas law with copied temperature
            soln%V(1, i, grid%jg_low+1) = p_ghost1 / (R_gas * T_interior)
            
            ! Second ghost cell (j = grid%jg_low)
            ! Get interior velocity from j = 2
            u_interior = soln%V(2, i, 2)
            v_interior = soln%V(3, i, 2)
            
            ! Reflect velocity
            dot_product = u_interior * nx + v_interior * ny
            u_reflected = u_interior - two * dot_product * nx
            v_reflected = v_interior - two * dot_product * ny
            
            soln%V(2, i, grid%jg_low) = u_reflected
            soln%V(3, i, grid%jg_low) = v_reflected
            
            ! Copy temperature from j = 2
            T_interior = soln%V(4, i, 2) / (R_gas * soln%V(1, i, 2))
            
            ! Pressure: 2nd order extrapolation
            p_ghost2 = two * p_ghost1 - soln%V(4, i, 1)
            soln%V(4, i, grid%jg_low) = p_ghost2
            
            ! Density from ideal gas law
            soln%V(1, i, grid%jg_low) = p_ghost2 / (R_gas * T_interior)
        end do
        
        ! Top wall: from i=i_inlet_max+1 to grid%imax - 2 ghost cells above
        do i = i_inlet_max+1, grid%imax
            ! Get wall normal vector (pointing into domain, downward for top wall)
            nx = grid%n_eta(i, grid%jmax, 1)  ! Normal x-component at top face
            ny = grid%n_eta(i, grid%jmax, 2)  ! Normal y-component at top face
            
            ! First ghost cell (j = grid%jg_high-1)
            ! Get interior velocity from j = grid%jmax
            u_interior = soln%V(2, i, grid%jmax)
            v_interior = soln%V(3, i, grid%jmax)
            
            ! Reflect velocity
            dot_product = u_interior * nx + v_interior * ny
            u_reflected = u_interior - two * dot_product * nx
            v_reflected = v_interior - two * dot_product * ny
            
            soln%V(2, i, grid%jg_high-1) = u_reflected
            soln%V(3, i, grid%jg_high-1) = v_reflected
            
            ! Copy temperature from j = grid%jmax
            T_interior = soln%V(4, i, grid%jmax) / (R_gas * soln%V(1, i, grid%jmax))
            
            ! Pressure: 2nd order extrapolation
            p_ghost1 = two * soln%V(4, i, grid%jmax) - soln%V(4, i, grid%jmax-1)
            soln%V(4, i, grid%jg_high-1) = p_ghost1
            
            ! Density from ideal gas law
            soln%V(1, i, grid%jg_high-1) = p_ghost1 / (R_gas * T_interior)
            
            ! Second ghost cell (j = grid%jg_high)
            ! Get interior velocity from j = grid%jmax-1
            u_interior = soln%V(2, i, grid%jmax-1)
            v_interior = soln%V(3, i, grid%jmax-1)
            
            ! Reflect velocity
            dot_product = u_interior * nx + v_interior * ny
            u_reflected = u_interior - two * dot_product * nx
            v_reflected = v_interior - two * dot_product * ny
            
            soln%V(2, i, grid%jg_high) = u_reflected
            soln%V(3, i, grid%jg_high) = v_reflected
            
            ! Copy temperature from j = grid%jmax-1
            T_interior = soln%V(4, i, grid%jmax-1) / (R_gas * soln%V(1, i, grid%jmax-1))
            
            ! Pressure: 2nd order extrapolation
            p_ghost2 = two * p_ghost1 - soln%V(4, i, grid%jmax)
            soln%V(4, i, grid%jg_high) = p_ghost2
            
            ! Density from ideal gas law
            soln%V(1, i, grid%jg_high) = p_ghost2 / (R_gas * T_interior)
        end do
        
        ! Left wall: from j=1 to grid%jmax - 2 ghost cells to the left
        do j = 1, grid%jmax
            ! Get wall normal vector (pointing into domain, rightward for left wall)
            nx = grid%n_xi(1, j, 1)  ! Normal x-component at left face
            ny = grid%n_xi(1, j, 2)  ! Normal y-component at left face
            
            ! First ghost cell (i = grid%ig_low+1)
            ! Get interior velocity from i = 1
            u_interior = soln%V(2, 1, j)
            v_interior = soln%V(3, 1, j)
            
            ! Reflect velocity
            dot_product = u_interior * nx + v_interior * ny
            u_reflected = u_interior - two * dot_product * nx
            v_reflected = v_interior - two * dot_product * ny
            
            soln%V(2, grid%ig_low+1, j) = u_reflected
            soln%V(3, grid%ig_low+1, j) = v_reflected
            
            ! Copy temperature from i = 1
            T_interior = soln%V(4, 1, j) / (R_gas * soln%V(1, 1, j))
            
            ! Pressure: 2nd order extrapolation
            p_ghost1 = two * soln%V(4, 1, j) - soln%V(4, 2, j)
            soln%V(4, grid%ig_low+1, j) = p_ghost1
            
            ! Density from ideal gas law
            soln%V(1, grid%ig_low+1, j) = p_ghost1 / (R_gas * T_interior)
            
            ! Second ghost cell (i = grid%ig_low)
            ! Get interior velocity from i = 2
            u_interior = soln%V(2, 2, j)
            v_interior = soln%V(3, 2, j)
            
            ! Reflect velocity
            dot_product = u_interior * nx + v_interior * ny
            u_reflected = u_interior - two * dot_product * nx
            v_reflected = v_interior - two * dot_product * ny
            
            soln%V(2, grid%ig_low, j) = u_reflected
            soln%V(3, grid%ig_low, j) = v_reflected
            
            ! Copy temperature from i = 2
            T_interior = soln%V(4, 2, j) / (R_gas * soln%V(1, 2, j))
            
            ! Pressure: 2nd order extrapolation
            p_ghost2 = two * p_ghost1 - soln%V(4, 1, j)
            soln%V(4, grid%ig_low, j) = p_ghost2
            
            ! Density from ideal gas law
            soln%V(1, grid%ig_low, j) = p_ghost2 / (R_gas * T_interior)
        end do
        
        ! Convert primitive to conservative variables
        call limit_primitives(soln%V)
        call prim2cons(soln%U, soln%V)
        
    end subroutine apply_slip_walls

    subroutine apply_outflow_bc(grid, soln)
        type(grid_t), intent(in) :: grid
        type(soln_t), intent(inout) :: soln
        
        integer :: i, j, n
        
        ! Apply 2nd order extrapolation to right boundary ghost cells
        do i = grid%ig_high-1, grid%ig_high  ! Right ghost cells
            do j = 1, grid%jmax
                do n = 1, neq
                    ! 2nd order linear extrapolation: f(i) = 2*f(i-1) - f(i-2)
                    soln%V(n, i, j) = two * soln%V(n, i-1, j) - soln%V(n, i-2, j)
                end do
            end do
        end do
        
        ! Convert to conserved variables
        call limit_primitives(soln%V)
        call prim2cons(soln%U, soln%V)
        
    end subroutine apply_outflow_bc
    
end module inlet_boundary