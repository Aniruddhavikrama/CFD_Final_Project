program main
    use set_precision, only : prec
    use set_constants, only : set_derived_constants,one
    use set_inputs, only : set_derived_inputs, cartesian_grid_flag, Lmms, residual_out_freq
    use grid_type
    use geometry, only : write_tecplot_file, cartesian_grid, read_grid, write_solution_dat
    use init, only : initialize_mms
    use soln_type, only : soln_t, allocate_soln, deallocate_soln
    use time_module, only : explicit_RK, calc_residual, residual_norms, calc_time_step
  
    implicit none
  
    type(grid_t) :: grid
    type(soln_t) :: soln
    integer :: i, j, iter, max_iter
    character(*), parameter :: grid_file = 'curv2d9_1.x'
    integer :: ierr
    real(prec), dimension(4) :: Rnorm
    ! real(prec), parameter :: CFL = 0.5_prec ! Example CFL number, adjust in set_inputs if needed
  
    ! Set maximum iterations for MMS convergence
    max_iter = 10000 ! Adjust as needed for convergence
  
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
  
    ! Initialize MMS solution
    call initialize_mms(grid, soln)
  
    ! Initialize initial residual norms
    call calc_time_step(grid, soln)
    call calc_residual(grid, soln)
    ! call residual_norms(soln%R, soln%rinit, one) ! Set initial norms
    call residual_norms(soln%R, soln%rinit, Rnorm, grid)


  
    ! Write initial solution
    call write_solution_dat(grid, soln, "mms_initial.dat")
  
    ! Main time-stepping loop
    write(*, '(A)') 'Iteration  Continuity    X-Momentum    Y-Momentum    Energy'
! Main loop
    do iter = 1, max_iter
        call explicit_RK(grid, soln)
        call calc_time_step(grid, soln)
        call calc_residual(grid, soln)
        call residual_norms(soln%R, Rnorm, soln%rinit, grid)
    
        if (mod(iter, residual_out_freq) == 0) then
            write(*, '(I8, 4ES14.6)') iter, Rnorm(1), Rnorm(2), Rnorm(3), Rnorm(4)
        end if
    
        if (maxval(Rnorm) < 1.0e-10_prec) then
            write(*, '(A, I8)') 'Converged at iteration ', iter
            exit
        end if
    end do
  
    ! Write final solution
    call write_solution_dat(grid, soln, "mms_final.dat")
  
    ! Clean up
    call deallocate_soln(soln)
    call deallocate_grid(grid)
  
  end program main