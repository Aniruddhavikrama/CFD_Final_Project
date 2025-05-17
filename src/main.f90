program main
    use set_precision, only : prec
    use set_constants
    use set_inputs
    use grid_type
    use geometry, only : write_tecplot_file, cartesian_grid, read_grid, write_solution_dat,write_solution_tecplot,&
    compute_discretization_error
    use init, only : initialize_mms
    use mms_boundary
    use soln_type
    use flux
    use time_module, only : explicit_RK, calc_residual, residual_norms, calc_time_step,evaluate_mms_source
    use variable_conversion
    use limiters
  
    implicit none
  
    type(grid_t) :: grid
    type(soln_t) :: soln
    integer ::  iter, max_iter
    character(*), parameter :: grid_file = 'curv2d33_edited.x'
    integer :: ierr
    real(prec), dimension(4) :: Rnorm
    ! real(prec), parameter :: CFL = 0.5_prec  adjust in set_inputs if needed
  
    ! Set maximum iterations for MMS convergence
    max_iter = 100000 ! Adjust as needed for convergence
  
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

    call initialize_mms(grid, soln) ! Sets soln%V to MMS, then soln%U from soln%V
    call apply_mms_boundary(grid,soln) ! Extrapolates soln%V to ghost, updates soln%U in ghost
    call update_states   (soln, grid) ! Syncs U,V, limits V, updates U, calc asnd

    ! Calculate MMS source terms for the initial residual calculation
    call evaluate_mms_source(grid, soln) ! Calculates soln%Smms
    soln%S = soln%Smms                   ! Assign Smms to S

    ! Calculate fluxes based on the initial state
    call compute_fluxes(grid,soln)    ! Calculates soln%Fxi, soln%Feta
  
    ! Calculate initial residual for normalization
    call calc_time_step(grid, soln) ! Calculates soln%dt
    call calc_residual(grid, soln)    ! Calculates soln%R using current S, Fxi, Feta
    ! call residual_norms(soln%R, soln%rinit, one) ! Set initial norms
    ! call residual_norms(soln%R,Rnorm, soln%rinit, grid)

    call residual_norms(soln%R,soln%rinit, [one,one,one,one], grid)


  
    ! Write initial solution
    call write_solution_dat(grid, soln, "mms_initial.dat")
  
    ! Main time-stepping loop
    write(*, '(A)') 'Iteration  Continuity    X-Momentum    Y-Momentum    Energy'
! Main loop
    do iter = 1, max_iter
        call explicit_RK(grid, soln)
        ! call apply_mms_boundary(grid, soln)
          ! 2) reconstruct face‐values & compute fluxes on those faces
        ! call compute_fluxes(grid, soln)
        call calc_time_step(grid, soln)
        call calc_residual(grid, soln)
        call residual_norms(soln%R, Rnorm, soln%rinit, grid)
    
        if (mod(iter, residual_out_freq) == 0) then
            write(*, '(I8, 4ES14.6)') iter, Rnorm(1), Rnorm(2), Rnorm(3), Rnorm(4)
        end if
    
        if (maxval(Rnorm) < 1.0e-8_prec) then
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