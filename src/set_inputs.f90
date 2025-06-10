module set_inputs
  use set_precision, only : prec
  use set_constants, only : zero,one,two,half,pi
  implicit none

  private

  public :: imax, jmax, neq, xmin, xmax, ymin, ymax, n_ghost
  public :: i_high, i_low, ig_high, ig_low
  public :: j_high, j_low, jg_high, jg_low
  public :: i_cell_high, i_cell_low, j_cell_high, j_cell_low
  public :: ig_cell_high, ig_cell_low, jg_cell_high, jg_cell_low
  public :: set_derived_inputs,CFL
  public :: cartesian_grid_flag,Lmms,kappa



  integer :: imax           = 2
  integer :: jmax           = 2
  integer :: i_low          = 0  ! Now will represent face index lower bound
  integer :: i_high         = 0  ! Now will represent face index upper bound
  integer :: i_cell_low     = 0  ! New: cell index lower bound
  integer :: i_cell_high    = 0  ! New: cell index upper bound
  integer :: ig_low         = 0
  integer :: ig_high        = 0
  integer :: j_low          = 0  ! Now will represent face index lower bound
  integer :: j_high         = 0  ! Now will represent face index upper bound
  integer :: j_cell_low     = 0  ! New: cell index lower bound
  integer :: j_cell_high    = 0  ! New: cell index upper bound
  integer :: jg_cell_low     = 0  ! New: cell index lower bound
  integer :: jg_cell_high    = 0  ! New: cell index upper bound
  integer :: ig_cell_low     = 0  ! New: cell index lower bound
  integer :: ig_cell_high    = 0  ! New: cell index upper bound
  integer :: jg_low         = 0
  integer :: jg_high        = 0
  integer :: neq            = 4
  integer :: n_ghost        = 2
  integer :: cartesian_grid_flag = 0  ! Switch: 1 for Cartesian, 0 for curvilinear from file
  integer,public :: order =2
  integer,public :: flux_scheme =1
  integer,public :: residual_out_freq= 1
  integer,public :: limiter_scheme = 1  ! 0: none, 1: Van Leer, 2: Minmod ,3: Venkatakrishnan, 4: Beta
  integer,public :: freeze_limiter_iter =100000 ! Freeze limiters after this iteration
  integer,public :: max_iter = 300000
  integer,public :: soln_out_freq =1000
  real(prec), public :: tol = 1.0e-10_prec
  integer,public :: second_order_switch_iterations = 1000




  real(prec) :: xmin       = zero 
  real(prec) :: xmax       = one
  real(prec) :: ymin       = zero
  real(prec) :: ymax       = one
  real(prec) :: CFL =       0.1_prec
  real(prec) :: Lmms     = 1.0_prec
  real(prec) :: kappa   = -1.0_prec !!!MUSCL Scheme
   
  
contains

subroutine set_derived_inputs
  
  ! a_inf   = sqrt(gamma*R_gas*T_inf)
  ! u_inf   = M_inf*a_inf
  ! u0 = u_inf*cos((pi/180.0_prec)*alpha)
  ! v0 = u_inf*sin((pi/180.0_prec)*alpha)
  ! rho_inf = p_inf/(R_gas*T_inf)
! Define cell indices
i_cell_low  = 1
i_cell_high = imax  ! 3 cells: i = 1, 2, 3
j_cell_low  = 1
j_cell_high = jmax  ! 3 cells: j = 1, 2, 3

ig_cell_low  = 1-n_ghost
ig_cell_high = imax+n_ghost  ! 3 cells: i = 1, 2, 3
jg_cell_low  = 1-n_ghost
jg_cell_high = jmax+n_ghost  ! 3 cells: j = 1, 2, 3

! Define face indices (i_low, i_high are now face indices)
i_low  = i_cell_low      ! First face at i = 1
i_high = i_cell_high + 1 ! Last face at i = 4 (for 3 cells)
j_low  = j_cell_low      ! First face at j = 1
j_high = j_cell_high + 1 ! Last face at j = 4

! Define vertex indices (including ghost cells)
ig_low  = i_low - n_ghost   ! -1
jg_low  = j_low - n_ghost   ! -1
ig_high = i_high + n_ghost  ! 5
jg_high = j_high + n_ghost  ! 5
  
end subroutine set_derived_inputs
  
end module set_inputs