program main
    use set_precision, only : prec
    use set_constants
    use set_inputs
    use grid_type
    use geometry, only : write_tecplot_file, cartesian_grid, read_grid, write_solution_dat,write_solution_tecplot,&
    compute_discretization_error, write_grid_and_cells_tecplot,write_solution_tecplot_iterations, write_soln_header
    use init, only : initialize_mms
    use mms_boundary
    use soln_type
    use flux
    use time_module, only : explicit_RK, calc_residual, residual_norms, calc_time_step,evaluate_mms_source
    use variable_conversion
    use limiters
    use residual_io
    use inlet_boundary
  
    implicit none
  
    type(grid_t) :: grid
    type(soln_t) :: soln
    integer ::  iter
    character(*), parameter :: grid_file = 'Inlet_fixed.53x17.grd '
    integer :: ierr
    real(prec), dimension(4) :: Rnorm
  
    ! Initialize constants and inputs
    call set_derived_constants
    call set_derived_inputs
  
    ! Set up grid
    if (cartesian_grid_flag == 1) then
      call cartesian_grid(grid)
    else
      call read_grid(grid, trim(grid_file), ierr)
    end if
    call cell_geometry(grid)
  
    ! Write grid output
    call write_tecplot_file(grid, "grid_output.dat")

  
    ! Allocate solution
    call allocate_soln(soln, grid)
  
  
    ! Initialize MMS solutions
    call init_limiters(grid)

    ! call initialize_mms(grid, soln) ! Sets soln%V to MMS, then soln%U from soln%V
    call initialize_flow(grid,soln)
    ! call apply_mms_boundary(grid,soln) ! Extrapolates soln%V to ghost, updates soln%U in ghost
    call apply_inlet_bc(grid,soln)
    call apply_slip_walls(grid,soln)
    call apply_outflow_bc(grid,soln)
    call update_states   (soln, grid) ! Syncs U,V, limits V, updates U, calc asnd

    ! Calculate MMS source terms for the initial residual calculation
    ! call evaluate_mms_source(grid, soln) ! Calculates soln%Smms
    ! soln%S = soln%Smms                   ! Assign Smms to S

    ! Calculate fluxes based on the initial state
    call compute_fluxes(grid,soln)    ! Calculates soln%Fxi, soln%Feta
  
    ! Calculate initial residual for normalization
    call calc_time_step(grid, soln) ! Calculates soln%dt
    call calc_residual(grid, soln)    ! Calculates soln%R using current S, Fxi, Feta
    call residual_norms(soln%R,soln%rinit, [one,one,one,one], grid)


  
    ! Write initial solution
    call write_solution_dat(grid, soln, "mms_initial.dat")
    call write_grid_and_cells_tecplot(grid, 'grid_with_cells.dat')

    ! create file for solution evolution
    call write_soln_header('solution_evolution.dat')
    ! Main time-stepping loop
    write(*, '(A)') 'Iteration  Continuity    X-Momentum    Y-Momentum    Energy'
! Main loop
    do iter = 1, max_iter
        if (iter == second_order_switch_iterations .and. order == 1) then
            write(*, '(A, I8)') 'Switching to second-order at iteration ', iter
            order = 2
            ! Optionally, you might want to recalculate some initialization here
            call update_states(soln, grid) ! If needed to recalculate limiters
        end if
        call explicit_RK(grid, soln)
        ! call apply_mms_boundary(grid, soln)
          ! 2) reconstruct face‐values & compute fluxes on those faces
        ! call compute_fluxes(grid, soln)
        call calc_time_step(grid, soln)
        call calc_residual(grid, soln)
        call residual_norms(soln%R, Rnorm, soln%rinit, grid)
        call write_solution_tecplot_iterations(grid, soln, 'solution_evolution.dat', iter)
    
        if (mod(iter, residual_out_freq) == 0) then
            write(*, '(I8, 4ES14.6)') iter, Rnorm(1), Rnorm(2), Rnorm(3), Rnorm(4)
            call write_residuals(iter, Rnorm) ! Write residuals to residual.dat
        end if
    
        if (maxval(Rnorm) < tol) then
            write(*, '(A, I8)') 'Converged at iteration ', iter
            call write_solution_tecplot(grid, soln, 'converged_solution.dat')
            call compute_discretization_error(grid, soln)

            exit
        end if
    end do
  
    ! Write final solution
    call write_solution_dat(grid, soln, "mms_final.dat")
    call write_solution_tecplot(grid, soln, 'converged_solution.dat')
  
    ! Clean up
    call deallocate_soln(soln)
    call deallocate_grid(grid)
  
  end program main