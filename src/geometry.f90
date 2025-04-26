module geometry

  use set_inputs
  use grid_type
  
  implicit none
  
contains

subroutine cartesian_grid(grid)
  use set_precision, only : prec
  use set_constants, only : one
  use set_inputs, only : imax, jmax, xmin, xmax, ymin, ymax, n_ghost
  use set_inputs, only : i_low, i_high, j_low, j_high
  use set_inputs, only : ig_low, ig_high, jg_low, jg_high
  use set_inputs, only : i_cell_low, i_cell_high, j_cell_low, j_cell_high
  integer :: i,j

  type(grid_t),intent(inout) ::grid


! Set grid indices before allocation
  grid%i_low   = i_low
  grid%j_low   = j_low
  grid%i_high  = i_high
  grid%j_high  = j_high
  grid%ig_low  = ig_low
  grid%jg_low  = jg_low
  grid%ig_high = ig_high
  grid%jg_high = jg_high
  grid%i_cell_low  = i_cell_low
  grid%i_cell_high = i_cell_high
  grid%j_cell_low  = j_cell_low
  grid%j_cell_high = j_cell_high
  grid%imax = imax
  grid%jmax = jmax
  grid%n_ghost = n_ghost

  call allocate_grid(grid)




  do j = grid%j_low, grid%j_high
    do i = grid%i_low, grid%i_high
        ! For imax=3, this gives vertices at x=0, 0.5, 1.0
        grid%x(i,j) = xmin + real(i-grid%i_low, prec) * (xmax-xmin) / real(grid%i_high-grid%i_low, prec)
        grid%y(i,j) = ymin + real(j-grid%j_low, prec) * (ymax-ymin) / real(grid%j_high-grid%j_low, prec)
    end do
end do

end subroutine cartesian_grid    

subroutine write_tecplot_file(grid, filename)
  use set_precision, only : prec
  use set_constants, only : zero
  type(grid_t), intent(in) :: grid
  character(len=*), intent(in) :: filename
  integer :: i, j, unit
  real(prec) :: a_xi, n_xi_x, n_xi_y, a_eta, n_eta_x, n_eta_y

  ! Open the file
  open(newunit=unit, file=trim(filename), status='replace', action='write')

  ! Write the Tecplot header
  write(unit, '(A)') 'TITLE = "2D Cartesian Grid with Face Areas and Normals"'
  write(unit, '(A)') 'VARIABLES = "X", "Y", "A_xi", "n_xi_x", "n_xi_y", "A_eta", "n_eta_x", "n_eta_y"'
  write(unit, '(A,I0,A,I0,A)') 'ZONE I=', grid%i_high - grid%i_low + 1, &
                               ', J=', grid%j_high - grid%j_low + 1, &
                               ', DATAPACKING=POINT'

  ! Write the data for each vertex
  do j = grid%j_low, grid%j_high
      do i = grid%i_low, grid%i_high
          a_xi = grid%A_xi(i,j)
          n_xi_x = grid%n_xi(i,j,1)
          n_xi_y = grid%n_xi(i,j,2)
          a_eta = grid%A_eta(i,j)
          n_eta_x = grid%n_eta(i,j,1)
          n_eta_y = grid%n_eta(i,j,2)

          ! Write the data
          write(unit, '(8E15.7)') grid%x(i,j), grid%y(i,j), &
                                  a_xi, n_xi_x, n_xi_y, a_eta, n_eta_x, n_eta_y
      end do
  end do

  close(unit)
end subroutine write_tecplot_file


end module geometry