module init
    use set_precision        , only: prec
    use grid_type            , only: grid_t
    use soln_type            , only: soln_t
    use set_inputs           , only: Lmms
    use fluid_constants      , only: gamma
    use set_constants        , only: zero, one, two, half
    use mms_functions        , only: wrap_rho_mms, wrap_uvel_mms, wrap_vvel_mms, wrap_press_mms
    use variable_conversion, only : prim2cons
    implicit none

    private
    public :: initialize_mms

contains

    subroutine initialize_mms(grid, soln)
        type(grid_t), intent(in   ) :: grid
        type(soln_t), intent(inout) :: soln
        integer :: i, j
        real(prec) :: x, y, rho, u, v, p, E
        !— Loop over physical cells only
        do j = grid%j_cell_low, grid%j_cell_high
            do i = grid%i_cell_low, grid%i_cell_high
                x = grid%xc(i,j)  ! Cell-center coordinates
                y = grid%yc(i,j)

                ! print *, "i=", i, "j=", j, "x=", x, "y=", y
                soln%V(1,i,j) = wrap_rho_mms(x, y)
                soln%V(2,i,j) = wrap_uvel_mms(x, y)
                soln%V(3,i,j) = wrap_vvel_mms(x, y)
                soln%V(4,i,j) = wrap_press_mms(x, y)
                ! print *, "rho=", soln%V(1,i,j), "u=", soln%V(2,i,j), &
                !          "v=", soln%V(3,i,j), "p=", soln%V(4,i,j)
            end do
        end do

        !— Convert to conservative variables 
        call prim2cons(soln%U, soln%V)

    end subroutine initialize_mms

end module init