module time_module

    use set_precision,only: prec
    use set_constants,only : zero,one,half,third,two,fourth
    use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
    use set_inputs,    only : j_low, j_high, jg_low, jg_high,CFL
    use grid_type, only : grid_t
    use soln_type, only : soln_t
    use variable_conversion, only : speed_of_sound

    implicit none

    contains
    subroutine calc_time_step(grid, soln)
        type(grid_t), intent(in)    :: grid
        type(soln_t), intent(inout) :: soln
        ! real(prec), dimension(grid%i_low:grid%i_high, grid%j_low:grid%j_high-1) :: Lam_xi
        ! real(prec), dimension(grid%i_low:grid%i_high-1, grid%j_low:grid%j_high) :: Lam_eta
        real(prec), dimension(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) :: Lam_xi, Lam_eta
        integer :: i,j


    
        do j = grid%j_cell_low, grid%j_cell_high
            do i = grid%i_cell_low, grid%i_cell_high
              call speed_of_sound(soln%V(4,i,j), soln%V(1,i,j), soln%asnd(i,j))
            end do
          end do
    
        Lam_xi  = abs(soln%V(2,:,:) * grid%n_xi_avg(:,:,1) + &
                      soln%V(3,:,:) * grid%n_xi_avg(:,:,2)) + soln%asnd
        Lam_eta = abs(soln%V(2,:,:) * grid%n_eta_avg(:,:,1) + &
                      soln%V(3,:,:) * grid%n_eta_avg(:,:,2)) + soln%asnd
        soln%dt = CFL * grid%V / ( &
                   Lam_xi * grid%A_xi(grid%i_low:grid%i_high, grid%j_low:grid%j_high-1) + &
                   Lam_eta * grid%A_eta(grid%i_low:grid%i_high-1, grid%j_low:grid%j_high) )
      end subroutine calc_time_step
    
      subroutine explicit_RK(grid, soln)
        type(grid_t), intent(inout) :: grid
        type(soln_t), intent(inout) :: soln
        real(prec), dimension(4) :: k
        integer :: i, j
    
        k = (/ fourth, third, half, one /)
        do j = 1, 4
          ! Apply MMS Dirichlet boundary conditions
          call apply_mms_boundary(grid, soln)
          ! Set source term for MMS
          soln%S = soln%Smms
          ! Convert primitive to conserved variables
          call prim2cons(soln%U, soln%V)
          ! Compute fluxes
          call calc_flux_2D(grid, soln)
          ! Compute residual
          call calc_residual(grid, soln)
          ! Update conserved variables
          do i = 1, neq
            soln%U(i, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) = &
                 soln%U(i, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) - &
               k(j) * (soln%R(i,  grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) * &
                       soln%dt( grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high)) / &
               grid%V( grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high)
          end do
          ! Update primitive variables
          call update_states(soln)
        end do
      end subroutine explicit_RK
    
      subroutine calc_residual(grid, soln)
        type(grid_t), intent(in) :: grid
        type(soln_t), intent(inout) :: soln
        integer :: i, j
    
        do j = grid%j_cell_low, grid%j_cell_high
          do i = grid%i_cell_low, grid%i_cell_high
            soln%R(:, i, j) = grid%A_xi(i+1, j) * soln%Fxi(:, i, j) &
                            - grid%A_xi(i, j) * soln%Fxi(:, i-1, j) &
                            + grid%A_eta(i, j+1) * soln%Feta(:, i, j) &
                            - grid%A_eta(i, j) * soln%Feta(:, i, j-1) &
                            - grid%V(i, j) * soln%S(:, i, j)
          end do
        end do
      end subroutine calc_residual
    
      subroutine residual_norms(R, Rnorm, rinit,grid)
        type(grid_t),intent(in) :: grid
        real(prec), dimension(neq, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high), intent(in) :: R
        real(prec), dimension(neq), intent(in)  :: rinit
        real(prec), dimension(neq), intent(out) :: Rnorm
        real(prec) :: Linv
        integer :: i
    
        Linv = one / real(size(R(1,:,:)), prec)
        do i = 1, neq
          Rnorm(i) = sqrt(Linv * sum(R(i,:,:)**2))
        end do
        Rnorm = Rnorm / rinit
      end subroutine residual_norms
  
end module time_module    