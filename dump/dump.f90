module set_precision
    use iso_c_binding, only : c_float, c_double, c_long_double, c_float_complex, &
                              c_double_complex, c_long_double_complex,          &
                              c_short, c_int, c_long_long
    implicit none
    private
    public :: r4, r8, r16, rd, dp, prec
  
    integer, parameter :: r4  = c_float
    integer, parameter :: r8  = c_double
    integer, parameter :: r16 = c_long_double
    integer, parameter :: dp  = r8
    integer, parameter :: rd  = r8
    integer, parameter :: prec  = r8
  end module set_precision

  module set_constants

    use set_precision, only : prec
  
    implicit none
  
    private
  
    public :: zero, one, two, three, four, five, six, seven, eight, nine
    public :: half, third, fourth, fifth, sixth, tenth
    public :: pi, set_derived_constants
  
    real(prec), parameter :: zero   = 0.0_prec
    real(prec), parameter :: tenth  = 0.1_prec
    real(prec), parameter :: sixth  = 1.0_prec/6.0_prec
    real(prec), parameter :: fifth  = 0.2_prec
    real(prec), parameter :: fourth = 0.25_prec
    real(prec), parameter :: third  = 1.0_prec/3.0_prec
    real(prec), parameter :: half   = 0.5_prec
    real(prec), parameter :: one    = 1.0_prec
    real(prec), parameter :: two    = 2.0_prec
    real(prec), parameter :: three  = 3.0_prec
    real(prec), parameter :: four   = 4.0_prec
    real(prec), parameter :: five   = 5.0_prec
    real(prec), parameter :: six    = 6.0_prec
    real(prec), parameter :: seven  = 7.0_prec
    real(prec), parameter :: eight  = 8.0_prec
    real(prec), parameter :: nine   = 9.0_prec
    real(prec)            :: pi     = 3.0_prec
  
    contains
  
      subroutine set_derived_constants
  
        implicit none
  
        pi = acos(-1.0_prec)
  
      end subroutine set_derived_constants
  
  end module set_constants

  module set_inputs
    use set_precision, only : prec
    use set_constants, only : zero,one,two,half,pi
    implicit none

    private

    public :: imax, jmax, neq, xmin, xmax, ymin, ymax, n_ghost
    public :: i_high, i_low, ig_high, ig_low
    public :: j_high, j_low, jg_high, jg_low
    public :: i_cell_high, i_cell_low, j_cell_high, j_cell_low
    public :: set_derived_inputs



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
! Define cell indices
  i_cell_low  = 1
  i_cell_high = imax  ! 3 cells: i = 1, 2, 3
  j_cell_low  = 1
  j_cell_high = jmax  ! 3 cells: j = 1, 2, 3

  ! Define face indices (i_low, i_high are now face indices)
  i_low  = i_cell_low      ! First face at i = 1
  i_high = i_cell_high + 1 ! Last face at i = 4 (for 3 cells)
  j_low  = j_cell_low      ! First face at j = 1
  j_high = j_cell_high + 1 ! Last face at j = 4

  ! Define vertex indices (including ghost cells)
  ig_low  = i_cell_low - n_ghost   ! -1
  jg_low  = j_cell_low - n_ghost   ! -1
  ig_high = i_cell_high + n_ghost  ! 5
  jg_high = j_cell_high + n_ghost  ! 5
    
  end subroutine set_derived_inputs
    
end module set_inputs

module grid_type
    use set_precision, only :prec
    use set_constants, only : zero, one, two, half,fourth
    use set_inputs,    only : i_low, i_high, ig_low, ig_high
    use set_inputs,    only : j_low, j_high, jg_low, jg_high
    use set_inputs,    only : i_cell_low, i_cell_high, j_cell_low, j_cell_high
    
    implicit none
    

    private
  
    public :: grid_t, allocate_grid, deallocate_grid,cell_geometry

    type grid_t

    integer :: imax, jmax, n_ghost
    integer :: i_low, i_high, ig_low, ig_high
    integer :: j_low, j_high, jg_low, jg_high
    integer :: i_cell_low, i_cell_high  ! New
    integer :: j_cell_low, j_cell_high

    real(prec),allocatable,dimension(:,:) :: x,y
    real(prec),allocatable,dimension(:,:) :: xc,yc
    real(prec), allocatable,dimension(:,:) :: A_xi, A_eta
    real(prec), allocatable, dimension(:,:,:) :: n_xi,n_eta
    real (prec), allocatable, dimension(:,:,:) :: n_xi_avg
    real(prec),allocatable,dimension(:,:,:)  :: n_eta_avg
    real(prec),allocatable, dimension(:,:) :: V


    end type grid_t

contains
subroutine allocate_grid(grid)
  type(grid_t), intent(inout) :: grid

  ! Allocate vertex coordinates
  allocate(grid%x(grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high))
  allocate(grid%y(grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high))

  ! Allocate cell center coordinates (use cell indices)
  allocate(grid%xc(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high))
  allocate(grid%yc(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high))

  ! Allocate face areas (now directly use i_low:i_high, j_low:j_high)
  allocate(grid%A_xi(grid%i_low:grid%i_high, grid%j_low:grid%j_high))
  allocate(grid%A_eta(grid%i_low:grid%i_high, grid%j_low:grid%j_high))

  ! Allocate normal vectors at faces
  allocate(grid%n_xi(grid%i_low:grid%i_high, grid%j_low:grid%j_high, 2))
  allocate(grid%n_eta(grid%i_low:grid%i_high, grid%j_low:grid%j_high, 2))

  ! Allocate averaged normal vectors
  allocate(grid%n_xi_avg(grid%i_low:grid%i_high, grid%j_low:grid%j_high, 2))
  allocate(grid%n_eta_avg(grid%i_low:grid%i_high, grid%j_low:grid%j_high, 2))

  ! Allocate cell volumes (use cell indices)
  allocate(grid%V(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high))

  ! Initialize all arrays to zero
  grid%x = zero
  grid%y = zero
  grid%xc = zero
  grid%yc = zero
  grid%A_xi = zero
  grid%A_eta = zero
  grid%n_xi = zero
  grid%n_eta = zero
  grid%n_xi_avg = zero
  grid%n_eta_avg = zero
  grid%V = zero
end subroutine allocate_grid   

subroutine cell_geometry(grid)
  type(grid_t), intent(inout) :: grid
  integer :: i, j
  real(prec) :: dx, dy

  ! Compute cell centers (xc, yc) as the average of the four vertices
  do j = grid%j_cell_low, grid%j_cell_high
      do i = grid%i_cell_low, grid%i_cell_high
          grid%xc(i,j) = fourth * (grid%x(i,j) + grid%x(i+1,j) + grid%x(i,j+1) + grid%x(i+1,j+1))
          grid%yc(i,j) = fourth * (grid%y(i,j) + grid%y(i+1,j) + grid%y(i,j+1) + grid%y(i+1,j+1))
      end do
  end do

  ! Compute face areas and normal vectors (xi-direction, vertical faces)
  do j = grid%j_low, grid%j_high-1
      do i = grid%i_low, grid%i_high
          dx = grid%x(i,j+1) - grid%x(i,j)
          dy = grid%y(i,j+1) - grid%y(i,j)
          grid%n_xi(i,j,1) = dy
          grid%n_xi(i,j,2) = -dx
          grid%A_xi(i,j) = sqrt(grid%n_xi(i,j,1)**2 + grid%n_xi(i,j,2)**2)
          if (grid%A_xi(i,j) > zero) then
              grid%n_xi(i,j,1) = grid%n_xi(i,j,1) / grid%A_xi(i,j)
              grid%n_xi(i,j,2) = grid%n_xi(i,j,2) / grid%A_xi(i,j)
          else
              grid%n_xi(i,j,1) = zero
              grid%n_xi(i,j,2) = zero
          end if
          grid%n_xi_avg(i,j,:) = grid%n_xi(i,j,:)
      end do
  end do

  !  Compute face areas and normal vectors (eta-direction, horizontal faces)
  do j = grid%j_low, grid%j_high
      do i = grid%i_low, grid%i_high-1
          dx = grid%x(i+1,j) - grid%x(i,j)
          dy = grid%y(i+1,j) - grid%y(i,j)
          grid%n_eta(i,j,1) = dy
          grid%n_eta(i,j,2) = dx
          grid%A_eta(i,j) = sqrt(grid%n_eta(i,j,1)**2 + grid%n_eta(i,j,2)**2)
          if (grid%A_eta(i,j) > zero) then
              grid%n_eta(i,j,1) = -grid%n_eta(i,j,1) / grid%A_eta(i,j)
              grid%n_eta(i,j,2) = grid%n_eta(i,j,2) / grid%A_eta(i,j)
          else
              grid%n_eta(i,j,1) = zero
              grid%n_eta(i,j,2) = zero
          end if
          grid%n_eta_avg(i,j,:) = grid%n_eta(i,j,:)
      end do
  end do

  ! Compute cell volumes (area in 2D) using the shoelace formula
  do j = grid%j_cell_low, grid%j_cell_high
      do i = grid%i_cell_low, grid%i_cell_high
          grid%V(i,j) = half * abs( &
              (grid%x(i,j) * grid%y(i+1,j) + grid%x(i+1,j) * grid%y(i+1,j+1) + &
               grid%x(i+1,j+1) * grid%y(i,j+1) + grid%x(i,j+1) * grid%y(i,j)) - &
              (grid%y(i,j) * grid%x(i+1,j) + grid%y(i+1,j) * grid%x(i+1,j+1) + &
               grid%y(i+1,j+1) * grid%x(i,j+1) + grid%y(i,j+1) * grid%x(i,j)) &
          )
      end do
  end do

end subroutine cell_geometry

subroutine deallocate_grid( grid )
    
    type( grid_t ), intent( inout ) :: grid
    
    deallocate( grid%x, grid%y, grid%xc, grid%yc, grid%A_xi, grid%A_eta, &
                grid%n_xi, grid%n_eta,  grid%V )
    
  end subroutine deallocate_grid
    
end module grid_type


module geometry

    use set_inputs
    use grid_type
    
    implicit none
    private
  
    public :: cartesian_grid,write_tecplot_file
    
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

program main
  use set_precision, only : prec
  use set_constants, only : set_derived_constants
  use set_inputs,only : set_derived_inputs
  use grid_type
  use geometry, only : write_tecplot_file,cartesian_grid
  
  implicit none
  
  type(grid_t) :: grid
  integer :: i, j
  
  ! Initialize constants and derived inputs
  call set_derived_constants
  call set_derived_inputs
  
  ! Create the cartesian grid
  call cartesian_grid(grid)
  
  ! Calculate cell geometry (with fixed calculations)
  call cell_geometry(grid)
  
  ! Print out the area vectors and normals for each face
  print *, "Grid Area Vectors and Normals:"
  print *, "------------------------------"
  
  ! CHANGE: Print only for interior and boundary faces
  ! A_xi, n_xi for i_low:i_high+1, j_low:j_high
  ! Print for interior and right boundary faces
 ! Loop over faces for A_xi and n_xi
  do j = grid%j_low, grid%j_high-1
    do i = grid%i_low, grid%i_high
        print *, "Vertical Face (xi) at (", i, ",", j, ")"
        print*, "  A_xi  =", grid%A_xi(i,j)
        print *, "  n_xi  = (", grid%n_xi(i,j,1), ",", grid%n_xi(i,j,2), ")"
    end do
end do

! Correct loop for eta-direction faces
do j = grid%j_low, grid%j_high
    do i = grid%i_low, grid%i_high-1
        print*, "Horizontal Face (eta) at (", i, ",", j, ")"
        print *, "  A_eta =", grid%A_eta(i,j)
        print*, "  n_eta = (", grid%n_eta(i,j,1), ",", grid%n_eta(i,j,2), ")"
    end do
end do

! Special printing for top boundary faces (j=grid%j_high)
do i = grid%i_low, grid%i_high-1
    print *, "Top Boundary Face (eta) at (", i, ",", grid%j_high, ")"
    print*, "  A_eta =", grid%A_eta(i,grid%j_high)
    print *, "  n_eta = (", grid%n_eta(i,grid%j_high,1), ",", grid%n_eta(i,grid%j_high,2), ")"
end do
call write_tecplot_file(grid, "grid_output.dat")
  
  ! Clean up

  call deallocate_grid(grid)
  
end program main



