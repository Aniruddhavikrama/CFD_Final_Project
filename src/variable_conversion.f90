module variable_conversion

    use set_precision,only : prec
    use set_constants,only: zero,one,two,half
    use fluid_constants, only : R_gas,gamma


    implicit none

    private
  
    public :: speed_of_sound, prim2cons, cons2prim
    ! public :: update_states
contains

! subroutine update_states(soln)
!     type(soln_t), intent(inout) :: soln
!     call cons2prim(soln%U,soln%V)
!     !call limit_primitives(soln%V)
!     call prim2cons(soln%U, soln%V)
!     call speed_of_sound(soln%V(4,:,:),soln%V(1,:,:),soln%asnd)
    
!     soln%mach = sqrt(soln%V(2,:,:)**2 + soln%V(3,:,:)**2)/soln%asnd
!     soln%temp = soln%V(4,:,:)/( soln%V(1,:,:)*R_gas )

! end subroutine update_states  

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
    
end module variable_conversion