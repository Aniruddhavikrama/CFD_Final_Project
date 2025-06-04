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
    
    ! Inlet flow conditions
    real(prec), parameter :: P_inf = 12270.0_prec   ! Freestream pressure (Pa) - from table
    real(prec), parameter :: T_inf = 217.0_prec     ! Freestream temperature (K) - from table  
    real(prec), parameter :: M_inf = 4.0_prec       ! Freestream Mach number - from table
    real(prec), parameter :: R_gas = 287.0_prec     ! Gas constant for air (J/kg-K)
    real(prec), parameter :: y_inlet_max = 0.6_prec ! Maximum y-coordinate for inlet
    
contains

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
    do j = grid%j_cell_low, grid%j_cell_high  ! 1 to 16
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
        
        ! Apply inlet BC to cells (1,16) to (20,16)
        ! j = 16
        ! do i = 1, 20
        !     soln%V(1, i, j) = rho_inlet
        !     soln%V(2, i, j) = u_inlet
        !     soln%V(3, i, j) = v_inlet
        !     soln%V(4, i, j) = p_inlet 
            
        !     ! Convert to conserved variables
        !    call prim2cons(soln%U, soln%V)
        ! end do
        
        ! Apply same values to ghost cells j=17 and j=18 (2 ghost cells above inlet)
        do j = 17, 18
            do i = 1, 20
                soln%V(1, i, j) = rho_inlet
                soln%V(2, i, j) = u_inlet
                soln%V(3, i, j) = v_inlet
                soln%V(4, i, j) = p_inlet
                
                ! Convert to conserved variables
                ! call prim2cons(soln%U, soln%V)
            end do
        end do
! left wall inflow(for debugging)
        do j = 1, 16
            do i = -1,0
                soln%V(1, i, j) = rho_inlet
                soln%V(2, i, j) = u_inlet
                soln%V(3, i, j) = v_inlet
                soln%V(4, i, j) = p_inlet
                
                ! Convert to conserved variables

            end do
        end do
        call prim2cons(soln%U, soln%V)
        
    end subroutine apply_inlet_bc
    
subroutine apply_slip_walls(grid, soln)
    type(grid_t), intent(in) :: grid
    type(soln_t), intent(inout) :: soln
    
    integer :: i, j
    real(prec) :: p_ghost1, p_ghost2, T_ghost1, T_ghost2
    real(prec) :: vel(2), n_wall(2), vel_dot_n, n_dot_n
    real(prec) :: rho_ghost1, rho_ghost2
    
    ! Top wall: cells (21,16) to (52,16), ghost cells at j=17,18
    do i = 21, 52
        n_wall(:) = grid%n_xi(i,16,:)
        n_dot_n = sum(n_wall(:)**2)

        n_wall(:) = n_wall(:) / sqrt(n_dot_n)
   
        
        ! First ghost cell (j=17)
        p_ghost1 = 2.0_prec * soln%V(4,i,16) - soln%V(4,i,15)  ! 2nd-order extrapolation
        T_ghost1 = soln%V(4,i,16) / (soln%V(1,i,16) * R_gas)  ! Copy T from 1st interior
        rho_ghost1 = p_ghost1 / (R_gas * T_ghost1)
        
        vel(1) = soln%V(2,i,16)  ! u from 1st interior
        vel(2) = soln%V(3,i,16)  ! v from 1st interior
        vel_dot_n = sum(vel(:) * n_wall(:))
        
        soln%V(2,i,17) = vel(1) - 2.0_prec * vel_dot_n * n_wall(1)  ! Reflect normal velocity
        soln%V(3,i,17) = vel(2) - 2.0_prec * vel_dot_n * n_wall(2)
        soln%V(1,i,17) = rho_ghost1
        soln%V(4,i,17) = p_ghost1
        
        ! Second ghost cell (j=18)
        p_ghost2 = 2.0_prec * p_ghost1 - soln%V(4,i,16)  ! Specified formula
        T_ghost2 = soln%V(4,i,15) / (soln%V(1,i,15) * R_gas)  ! Copy T from 2nd interior
        rho_ghost2 = p_ghost2 / (R_gas * T_ghost2)
        
        vel(1) = soln%V(2,i,15)  ! u from 2nd interior
        vel(2) = soln%V(3,i,15)  ! v from 2nd interior
        vel_dot_n = sum(vel(:) * n_wall(:))
        
        soln%V(2,i,18) = vel(1) - 2.0_prec * vel_dot_n * n_wall(1)  ! Reflect normal velocity
        soln%V(3,i,18) = vel(2) - 2.0_prec * vel_dot_n * n_wall(2)
        soln%V(1,i,18) = rho_ghost2
        soln%V(4,i,18) = p_ghost2
        
        call prim2cons(soln%U, soln%V)
    end do
    
    ! Bottom wall: cells (1,1) to (52,1), ghost cells at j=0,-1
    do i = 1, 52
        n_wall(:) = grid%n_xi(i,1,:)
        n_dot_n = sum(n_wall(:)**2)

        n_wall(:) = n_wall(:) / sqrt(n_dot_n)

        
        ! First ghost cell (j=0)
        p_ghost1 = 2.0_prec * soln%V(4,i,1) - soln%V(4,i,2)
        T_ghost1 = soln%V(4,i,1) / (soln%V(1,i,1) * R_gas)  ! Copy T from 1st interior
        rho_ghost1 = p_ghost1 / (R_gas * T_ghost1)
        
        vel(1) = soln%V(2,i,1)
        vel(2) = soln%V(3,i,1)
        vel_dot_n = sum(vel(:) * n_wall(:))
        
        soln%V(2,i,0) = vel(1) - 2.0_prec * vel_dot_n * n_wall(1)
        soln%V(3,i,0) = vel(2) - 2.0_prec * vel_dot_n * n_wall(2)
        soln%V(1,i,0) = rho_ghost1
        soln%V(4,i,0) = p_ghost1
        
        ! Second ghost cell (j=-1)
        p_ghost2 = 2.0_prec * p_ghost1 - soln%V(4,i,1)
        T_ghost2 = soln%V(4,i,2) / (soln%V(1,i,2) * R_gas)  ! Copy T from 2nd interior
        rho_ghost2 = p_ghost2 / (R_gas * T_ghost2)
        
        vel(1) = soln%V(2,i,2)
        vel(2) = soln%V(3,i,2)
        vel_dot_n = sum(vel(:) * n_wall(:))
        
        soln%V(2,i,-1) = vel(1) - 2.0_prec * vel_dot_n * n_wall(1)
        soln%V(3,i,-1) = vel(2) - 2.0_prec * vel_dot_n * n_wall(2)
        soln%V(1,i,-1) = rho_ghost2
        soln%V(4,i,-1) = p_ghost2
        
        call prim2cons(soln%U, soln%V)
    end do
    
    ! Left wall: cells (1,1) to (1,16), ghost cells at i=0,-1
    ! do j = 1, 16
    !     n_wall(:) = grid%n_eta(1,j,:)
    !     n_dot_n = sum(n_wall(:)**2)
    !     n_wall(:) = n_wall(:) / sqrt(n_dot_n)
        
    !     ! First ghost cell (i=0)
    !     p_ghost1 = 2.0_prec * soln%V(4,1,j) - soln%V(4,2,j)
    !     T_ghost1 = soln%V(4,1,j) / (soln%V(1,1,j) * R_gas)  ! Copy T from 1st interior
    !     rho_ghost1 = p_ghost1 / (R_gas * T_ghost1)
        
    !     vel(1) = soln%V(2,1,j)
    !     vel(2) = soln%V(3,1,j)
    !     vel_dot_n = sum(vel(:) * n_wall(:))
        
    !     soln%V(2,0,j) = vel(1) - 2.0_prec * vel_dot_n * n_wall(1)
    !     soln%V(3,0,j) = vel(2) - 2.0_prec * vel_dot_n * n_wall(2)
    !     soln%V(1,0,j) = rho_ghost1
    !     soln%V(4,0,j) = p_ghost1
        
    !     ! Second ghost cell (i=-1)
    !     p_ghost2 = 2.0_prec * p_ghost1 - soln%V(4,1,j)
    !     T_ghost2 = soln%V(4,2,j) / (soln%V(1,2,j) * R_gas)  ! Copy T from 2nd interior
    !     rho_ghost2 = p_ghost2 / (R_gas * T_ghost2)
        
    !     vel(1) = soln%V(2,2,j)
    !     vel(2) = soln%V(3,2,j)
    !     vel_dot_n = sum(vel(:) * n_wall(:))
        
    !     soln%V(2,-1,j) = vel(1) - 2.0_prec * vel_dot_n * n_wall(1)
    !     soln%V(3,-1,j) = vel(2) - 2.0_prec * vel_dot_n * n_wall(2)
    !     soln%V(1,-1,j) = rho_ghost2
    !     soln%V(4,-1,j) = p_ghost2
        
    !     call prim2cons(soln%U, soln%V)
    ! end do
end subroutine apply_slip_walls

    
    subroutine apply_outflow_bc(grid, soln)
        type(grid_t), intent(in) :: grid
        type(soln_t), intent(inout) :: soln
        
        integer :: i, j, n
        
        ! Apply 2nd order extrapolation to ghost cells i=53 and i=54
        do i = 53, 54
            do j = 1, 16
                do n = 1, neq
                    ! 2nd order linear extrapolation: f(i) = 2*f(i-1) - f(i-2)
                    soln%V(n, i, j) = two * soln%V(n, i-1, j) - soln%V(n, i-2, j)
                end do
                
                ! Convert to conserved variables
                call limit_primitives(soln%V)
                call prim2cons(soln%U, soln%V)
            end do
        end do
        
    end subroutine apply_outflow_bc
    
end module inlet_boundary