module time_module

    use set_precision,only: prec
    use set_constants,only : zero,one,half,third,two,fourth
    use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
    use set_inputs,    only : j_low, j_high, jg_low, jg_high,CFL
    use grid_type, only : grid_t
    use soln_type, only : soln_t
    use variable_conversion, only : speed_of_sound,prim2cons,update_states
    use mms_boundary
    use mms_functions, only: wrap_rmassconv, wrap_xmtmconv, wrap_ymtmconv, wrap_energyconv
    use flux, only : compute_fluxes

    implicit none

    contains

    subroutine evaluate_mms_source(grid, soln)
        type(grid_t), intent(in) :: grid
        type(soln_t), intent(inout) :: soln
        integer :: i, j
        real(prec) :: x, y
      
        do j = grid%j_cell_low, grid%j_cell_high
          do i = grid%i_cell_low, grid%i_cell_high
            x = grid%xc(i,j)
            y = grid%yc(i,j)
            soln%Smms(1,i,j) = wrap_rmassconv(x, y)
            soln%Smms(2,i,j) = wrap_xmtmconv(x, y)
            soln%Smms(3,i,j) = wrap_ymtmconv(x, y)
            soln%Smms(4,i,j) = wrap_energyconv(x, y)

                    ! Optional: Debug printout (uncomment if needed)
            ! print *, 'Cell (i,j)=', i, j, 'x=', x, 'y=', y
            ! print *, '  Smms mass   =', soln%Smms(1,i,j)
            ! print *, '  Smms x-mom  =', soln%Smms(2,i,j)
            ! print *, '  Smms y-mom  =', soln%Smms(3,i,j)
            ! print *, '  Smms energy =', soln%Smms(4,i,j)
          end do
        end do
      end subroutine evaluate_mms_source
      
      subroutine calc_time_step(grid, soln)
        type(grid_t), intent(in)    :: grid
        type(soln_t), intent(inout) :: soln
        real(prec), allocatable :: Lam_xi(:,:), Lam_eta(:,:)
        real(prec) :: A_xi_avg, A_eta_avg
        integer :: i, j
    
        allocate(Lam_xi(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high))
        allocate(Lam_eta(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high))
    
        ! Compute local speed of sound
        do j = grid%j_cell_low, grid%j_cell_high
            do i = grid%i_cell_low, grid%i_cell_high
                call speed_of_sound(soln%V(4,i,j), soln%V(1,i,j), soln%asnd(i,j))
            end do
        end do
    
        ! Compute local wave speeds in xi and eta directions
        Lam_xi(:,:) = abs(soln%V(2, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) * &
                          grid%n_xi_avg(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high,1) + &
                          soln%V(3, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) * &
                          grid%n_xi_avg(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high,2)) + &
                          soln%asnd(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high)
    
        Lam_eta(:,:) = abs(soln%V(2, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) * &
                           grid%n_eta_avg(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high,1) + &
                           soln%V(3, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) * &
                           grid%n_eta_avg(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high,2)) + &
                           soln%asnd(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high)
    
        ! Compute local time step by averaging face areas to cell centers
        do j = grid%j_cell_low, grid%j_cell_high
            do i = grid%i_cell_low, grid%i_cell_high
                A_xi_avg  = 0.5_prec * (grid%A_xi(i,j) + grid%A_xi(i+1,j))
                A_eta_avg = 0.5_prec * (grid%A_eta(i,j) + grid%A_eta(i,j+1))
                soln%dt(i,j) = CFL * grid%V(i,j) / (Lam_xi(i,j) * A_xi_avg + Lam_eta(i,j) * A_eta_avg)
            end do
        end do
    
        deallocate(Lam_xi, Lam_eta)
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

        call update_states(soln, grid)
        ! Set source term for MMS
        call evaluate_mms_source(grid,soln)
        soln%S = soln%Smms
        ! Convert primitive to conserved variables
        call prim2cons(soln%U, soln%V)
        ! Compute fluxes
        call compute_fluxes(grid, soln)

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
        call update_states(soln,grid)
      end do
    end subroutine explicit_RK
  
    subroutine calc_residual(grid, soln)
      type(grid_t), intent(in) :: grid
      type(soln_t), intent(inout) :: soln
      integer :: i, j
      ! do j = grid%j_cell_low, grid%j_cell_high
      !     do i = grid%i_cell_low, grid%i_cell_high
      !         soln%R(:, i, j) = grid%A_xi(i+1, j) * soln%Fxi(:, i+1, j) - grid%A_xi(i, j) * soln%Fxi(:, i, j) &
      !                           + grid%A_eta(i, j+1) * soln%Feta(:, i, j+1) - grid%A_eta(i, j) * soln%Feta(:, i, j)- grid%V(i, j) * soln%S(:, i, j) 
      !     end do
      ! end do
!!separate
      ! do j = grid%j_cell_low, grid%j_cell_high
      !     do i = grid%i_cell_low, grid%i_cell_high
      !         soln%R(:, i, j) = -(grid%A_xi(i+1, j) * soln%Fxi(:, i+1, j) - grid%A_xi(i, j) * soln%Fxi(:, i, j) &
      !                           + grid%A_eta(i, j+1) * soln%Feta(:, i, j+1) - grid%A_eta(i, j) * soln%Feta(:, i, j))- grid%V(i, j) * soln%S(:, i, j) 
      !     end do
      ! end do
      
!!separate
        do j = grid%j_cell_low, grid%j_cell_high
            do i = grid%i_cell_low, grid%i_cell_high
                soln%R(:, i, j) = - grid%V(i, j) * soln%S(:, i, j)
            end do
        end do

        do j = grid%j_low,grid%j_high-1
            do i = grid%i_low,grid%i_high-1
                soln%R(:, i, j) = soln%R(:, i, j) &
                                + grid%A_xi(i+1, j) * soln%Fxi(:, i+1, j) &
                                - grid%A_xi(i,   j) * soln%Fxi(:, i,   j)
            end do
        end do

        do j = grid%j_low,grid%j_high-1
            do i = grid%i_low,grid%i_high-1
                soln%R(:, i, j) = soln%R(:, i, j)                         &
                                + grid%A_eta(i, j+1) * soln%Feta(:, i, j+1 ) &
                                - grid%A_eta(i, j  ) * soln%Feta(:, i, j   ) 
            end do
        end do

      !!!seperate

      
      ! do j = grid%j_low,grid%j_high
      !   do i = grid%i_low,grid%i_high
      !     soln%R(:, i, j) = soln%R(:, i, j) &
      !                     + grid%A_xi(i+1, j) * soln%Fxi(:, i+1, j) &
      !                     - grid%A_xi(i, j) * soln%Fxi(:, i, j) &
      !                     + grid%A_eta(i, j+1) * soln%Feta(:, i, j) &
      !                     - grid%A_eta(i, j) * soln%Feta(:, i, j-1) &
      !                     - grid%V(i, j) * soln%S(:, i, j)
      !   end do
      ! end do
    end subroutine calc_residual
  
    subroutine residual_norms(R, rnorm, rinit,grid)
      type(grid_t),intent(in) :: grid
      real(prec), dimension(neq, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high), intent(in) :: R
      real(prec), dimension(neq), intent(in)  :: rinit
      real(prec), dimension(neq), intent(out) :: rnorm
      real(prec) :: Linv
      integer :: i
  
      Linv = one / real(size(R(1,:,:)), prec)
      do i = 1, neq
        Rnorm(i) = sqrt(Linv * sum(R(i,:,:)**2))
      end do
      Rnorm = Rnorm / rinit
    end subroutine residual_norms
  
end module time_module    