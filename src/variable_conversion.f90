module variable_conversion

    use set_precision,only : prec
    use set_constants,only: zero,one,two,half
    use fluid_constants, only : R_gas,gamma
    use soln_type
    use grid_type


    implicit none

    private
  
    public :: speed_of_sound, prim2cons, cons2prim
    public :: update_states,limit_primitives
contains
subroutine update_states(soln, grid)
    type(soln_t), intent(inout) :: soln
    type(grid_t), intent(in)    :: grid
    integer :: i, j

    call cons2prim(soln%U, soln%V)
    call limit_primitives(soln%V)
    call prim2cons(soln%U, soln%V)

    do j = grid%j_cell_low, grid%j_cell_high
        do i = grid%i_cell_low, grid%i_cell_high
            call speed_of_sound(soln%V(4,i,j), soln%V(1,i,j), soln%asnd(i,j))
        end do
    end do

    ! soln%mach(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) = &
    !     sqrt(soln%V(2, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high)**2 + &
    !          soln%V(3, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high)**2) / &
    !     soln%asnd(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high)

    ! soln%temp(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) = &
    !     soln%V(4, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) / &
    !     (soln%V(1, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high) * R_gas)
end subroutine update_states



subroutine speed_of_sound(pressure,rho,sound_speed)
    real(prec), intent(in)  :: pressure, rho
    real(prec), intent(out) :: sound_speed
    sound_speed = sqrt(gamma*pressure/rho)

end subroutine speed_of_sound

subroutine prim2cons( U, V )
    
    real(prec), dimension(:,:,:), intent(out)  :: U
    real(prec), dimension(:,:,:), intent(in)   :: V
    
    U(1,:,:) = V(1,:,:)
    U(2,:,:) = V(1,:,:)*V(2,:,:)
    U(3,:,:) = V(1,:,:)*V(3,:,:)
    U(4,:,:) = V(4,:,:)/( gamma - one ) + half*V(1,:,:)*&
             & ( V(2,:,:)**2 + V(3,:,:)**2 )
    
  end subroutine prim2cons

  subroutine cons2prim( U, V )
    
    real(prec), dimension(:,:,:), intent(in) :: U
    real(prec), dimension(:,:,:), intent(out) :: V
    
    V(1,:,:) = U(1,:,:)
    V(2,:,:) = U(2,:,:)/U(1,:,:)
    V(3,:,:) = U(3,:,:)/U(1,:,:)
    V(4,:,:) = (gamma - one)*( U(4,:,:) - half*&
             & ( U(2,:,:)**2 + U(3,:,:)**2 )/U(1,:,:) )
    
  end subroutine cons2prim

  subroutine limit_primitives(V)
    
    real(prec), dimension(:,:,:), intent(inout) :: V
    integer :: i, j
    
    do j = lbound(V,3),ubound(V,3)
      do i = lbound(V,2),ubound(V,2)
        V(1,i,j) = max(0.0001_prec,V(1,i,j))
        V(4,i,j) = max(1.0_prec,V(4,i,j))
      end do
    end do
    
  end subroutine limit_primitives
    
end module variable_conversion