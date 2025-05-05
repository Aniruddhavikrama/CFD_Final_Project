module mms_boundary
    use set_precision
    use set_constants
    use soln_type,only : soln_t
    use grid_type,only : grid_t
    implicit none
    
contains
subroutine apply_mms_boundary(grid, soln)
    !> Apply Dirichlet MMS boundary conditions using analytic solution on ghost cells
    type(grid_t), intent(in)      :: grid
    type(soln_t), intent(inout)   :: soln
    integer :: i, j

! Left ghost columns: ig_low and ig_low+1
    do j = grid%jg_low, grid%jg_high
        soln%V(:, grid%ig_low+1, j) = 2 * soln%V(:, grid%ig_low+2, j) - soln%V(:, grid%ig_low+3, j)
        soln%V(:, grid%ig_low  , j) = 2 * soln%V(:, grid%ig_low+1, j) - soln%V(:, grid%ig_low+2, j)
      end do
  
      ! Right ghost columns: ig_high-1 and ig_high
      do j = grid%jg_low, grid%jg_high
        soln%V(:, grid%ig_high-1, j) = 2 * soln%V(:, grid%ig_high-2, j) - soln%V(:, grid%ig_high-3, j)
        soln%V(:, grid%ig_high  , j) = 2 * soln%V(:, grid%ig_high-1, j) - soln%V(:, grid%ig_high-2, j)
      end do
  
      ! Bottom ghost rows: jg_low and jg_low+1
      do i = grid%ig_low, grid%ig_high
        soln%V(:, i, grid%jg_low+1) = 2 * soln%V(:, i, grid%jg_low+2) - soln%V(:, i, grid%jg_low+3)
        soln%V(:, i, grid%jg_low  ) = 2 * soln%V(:, i, grid%jg_low+1) - soln%V(:, i, grid%jg_low+2)
      end do
  
      ! Top ghost rows: jg_high-1 and jg_high
      do i = grid%ig_low, grid%ig_high
        soln%V(:, i, grid%jg_high-1) = 2 * soln%V(:, i, grid%jg_high-2) - soln%V(:, i, grid%jg_high-3)
        soln%V(:, i, grid%jg_high  ) = 2 * soln%V(:, i, grid%jg_high-1) - soln%V(:, i, grid%jg_high-2)
      end do

  end subroutine apply_mms_boundary
    
end module mms_boundary