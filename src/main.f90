program main
    use set_precision, only : prec
    use set_constants, only : set_derived_constants
    use set_inputs,only : set_derived_inputs, cartesian_grid_flag
    use grid_type
    use geometry, only : write_tecplot_file, cartesian_grid, read_grid
    
    implicit none
    
    type(grid_t) :: grid
    integer :: i, j
    character(*), parameter :: grid_file = 'curv2d9_1.x'
    integer :: ierr
    
    ! Initialize constants and derived inputs
    call set_derived_constants
    call set_derived_inputs
    
    if (cartesian_grid_flag == 1) then
      ! build a uniform Cartesian grid
      call cartesian_grid(grid)
  
    else
      ! read a curvilinear grid from file
      ! grid_file = 'grid_input.dat'   ! <-- change to your actual filename
      call read_grid(grid, trim(grid_file), ierr)
      if (ierr /= 0) then
        print *, '*** Error: failed to read grid from ', trim(grid_file)
        stop 1
      end if
    end if
    
    ! Calculate cell geometry (with fixed calculations)
    call cell_geometry(grid)
    
    ! Print out the area vectors and normals for each face
    print *, "Grid Area Vectors and Normals:"
    print *, "------------------------------"
    
    ! CHANGE: Print only for interior and boundary faces
    ! A_xi, n_xi for i_low:i_high+1, j_low:j_high
    ! Print for interior and right boundary faces
   ! Loop over faces for A_xi and n_xi
    do j = grid%j_low, grid%j_high-1
      do i = grid%i_low, grid%i_high
          print *, "Vertical Face (xi) at (", i, ",", j, ")"
          print*, "  A_xi  =", grid%A_xi(i,j)
          print *, "  n_xi  = (", grid%n_xi(i,j,1), ",", grid%n_xi(i,j,2), ")"
      end do
  end do
  
  ! Correct loop for eta-direction faces
  do j = grid%j_low, grid%j_high
      do i = grid%i_low, grid%i_high-1
          print*, "Horizontal Face (eta) at (", i, ",", j, ")"
          print *, "  A_eta =", grid%A_eta(i,j)
          print*, "  n_eta = (", grid%n_eta(i,j,1), ",", grid%n_eta(i,j,2), ")"
      end do
  end do
  
  ! Special printing for top boundary faces (j=grid%j_high)
  do i = grid%i_low, grid%i_high-1
      print *, "Top Boundary Face (eta) at (", i, ",", grid%j_high, ")"
      print*, "  A_eta =", grid%A_eta(i,grid%j_high)
      print *, "  n_eta = (", grid%n_eta(i,grid%j_high,1), ",", grid%n_eta(i,grid%j_high,2), ")"
  end do
  call write_tecplot_file(grid, "grid_output.dat")
    
    ! Clean up
  
    call deallocate_grid(grid)
    
  end program main