program main
    use set_precision, only : prec
    use set_constants, only : set_derived_constants
    use set_inputs,only : set_derived_inputs, cartesian_grid_flag,Lmms
    use grid_type
    use geometry, only : write_tecplot_file, cartesian_grid, read_grid,write_solution_dat
    use init,only : initialize_mms
    use soln_type,only : soln_t,allocate_soln,deallocate_soln
    
    implicit none
    
    type(grid_t) :: grid
    type(soln_t) :: soln
    integer :: i, j
    character(*), parameter :: grid_file = 'curv2d9_1.x'
    integer :: ierr
    
    call set_derived_constants
    call set_derived_inputs
    
    if (cartesian_grid_flag == 1) then
      call cartesian_grid(grid)
    else
      call read_grid(grid, trim(grid_file), ierr)
    end if
    ! call ghost_cells(grid)
    call cell_geometry(grid)
    
    call write_tecplot_file(grid, "grid_output.dat")
    call allocate_soln(soln,grid)

    call initialize_mms(grid,soln)

    call write_solution_dat(grid,soln,"mms_out.dat")

    call deallocate_soln(soln)
    call deallocate_grid(grid)
    
  end program main