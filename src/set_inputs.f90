module set_inputs
    use set_precision, only : prec
    use set_constants, only : zero,one,two,half,pi
    implicit none

    private

    public :: imax, jmax, neq, xmin, xmax, ymin, ymax, n_ghost
    public :: i_high, i_low, ig_high, ig_low
    public :: j_high, j_low, jg_high, jg_low




    integer :: imax           = 10
    integer :: jmax           = 10
    integer :: i_low          = 0
    integer :: i_high         = 0
    integer :: ig_low         = 0
    integer :: ig_high        = 0
    integer :: j_low          = 0
    integer :: j_high         = 0
    integer :: jg_low         = 0
    integer :: jg_high        = 0
    integer :: neq            = 4
    integer :: n_ghost        = 2

    real(prec) :: xmin       = zero
    real(prec) :: xmax       = one
    real(prec) :: ymin       = zero
    real(prec) :: ymax       = one
  
    
contains

subroutine set_derived_inputs
    
    ! a_inf   = sqrt(gamma*R_gas*T_inf)
    ! u_inf   = M_inf*a_inf
    ! u0 = u_inf*cos((pi/180.0_prec)*alpha)
    ! v0 = u_inf*sin((pi/180.0_prec)*alpha)
    ! rho_inf = p_inf/(R_gas*T_inf)
    i_low = 1
    j_low = 1
    i_high = imax-1
    j_high = jmax-1
    ig_low  = i_low - n_ghost
    jg_low  = j_low - n_ghost
    ig_high = i_high + n_ghost
    jg_high = j_high + n_ghost
    
  end subroutine set_derived_inputs
    
end module set_inputs
