module time_module

    use set_precision,only: prec
    use set_constants,only : zero,one,half,third,two,fourth
    use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
    use set_inputs,    only : j_low, j_high, jg_low, jg_high
    use grid_type, only : grid_t
    use soln_type, only : soln_t
    use variable_conversion, only : speed_of_sound

    implicit none

    contains

    subroutine calc_time_step(grid,soln)
        use set_inputs, only : CFL
        type(grid_t),intent(in):: grid
        type(soln_t),intent(inout):: soln
        real(prec), dimension(ig_low:ig_high,jg_low:jg_high) :: Lam_xi,Lam_eta

        call speed_of_sound(soln%V(4,:,:),soln%V(1,:,:),soln%asnd)

        Lam_xi  = abs(soln%V(2,:,:)*grid%n_xi_avg(:,:,1)   + &
        soln%V(3,:,:)*grid%n_xi_avg(:,:,2))  + soln%asnd
        Lam_eta = abs(soln%V(2,:,:)*grid%n_eta_avg(:,:,1)  + &
        soln%V(3,:,:)*grid%n_eta_avg(:,:,2)) + soln%asnd
        soln%dt = CFL*grid%V/( &
        Lam_xi*grid%A_xi(ig_low:ig_high,jg_low:jg_high) + &
        Lam_eta*grid%A_eta(ig_low:ig_high,jg_low:jg_high) )
    end subroutine calc_time_step    
  
end module time_module    