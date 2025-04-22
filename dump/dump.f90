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
    public :: set_derived_inputs




    integer :: imax           = 11
    integer :: jmax           = 11
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

module grid_type
    use set_precision, only :prec
    use set_constants, only : zero, one, two, half
    use set_inputs,    only : i_low, i_high, ig_low, ig_high
    use set_inputs,    only : j_low, j_high, jg_low, jg_high
    
    implicit none
    

    private
  
    public :: grid_t, allocate_grid, deallocate_grid,cell_geometry

    type grid_t

    integer :: imax, jmax, n_ghost
    integer :: i_low, i_high, ig_low, ig_high
    integer :: j_low, j_high, jg_low, jg_high

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
    ! use set_inputs

    type(grid_t),intent(inout) ::grid

    allocate (              grid%x(ig_low:ig_high+1,jg_low:jg_high+1), &
                            grid%y(ig_low:ig_high+1,jg_low:jg_high+1), &
                            grid%xc(ig_low:ig_high+1,jg_low:jg_high+1), &
                            grid%yc(ig_low:ig_high,jg_low:jg_high), &
                            grid%A_xi(ig_low:ig_high,jg_low:jg_high), &
                            grid%A_eta(ig_low:ig_high,jg_low:jg_high), &
                            grid%n_xi(ig_low:ig_high,jg_low:jg_high,2),   &
                            grid%n_eta(ig_low:ig_high,jg_low:jg_high,2),   &
                            grid%n_xi_avg(ig_low:ig_high,jg_low:jg_high,2), &
                            grid%n_eta_avg(ig_low:ig_high,jg_low:jg_high,2), &
                            grid%V(ig_low:ig_high,jg_low:jg_high)    )
    
    grid%x= zero
    grid%y= zero


end subroutine allocate_grid    

subroutine cell_geometry(grid)
  use set_precision, only : prec
  use set_constants, only : zero
  
  type(grid_t), intent(inout) :: grid
  integer :: i, j
  real(prec) :: dx, dy,area,cx,cy
  
  ! do j = grid%jg_low, grid%jg_high
  !   do i = grid%ig_low, grid%ig_high
  !     ! Calculate xi face (j-direction face)
  !     ! Vector from (i,j) to (i,j+1)
  !     dx = grid%x(i,j+1) - grid%x(i,j)
  !     dy = grid%y(i,j+1) - grid%y(i,j)
      
  !     ! Area magnitude (length of the face)
  !     grid%A_xi(i,j) = sqrt(dx**2 + dy**2)
      
  !     ! Normal vector (ensure proper normalization)
  !     if (grid%A_xi(i,j) > zero) then
  !       ! For a cartesian grid, normal to xi face points in x direction
  !       grid%n_xi(i,j,1) = dy / grid%A_xi(i,j)  ! x-component
  !       grid%n_xi(i,j,2) = -dx / grid%A_xi(i,j) ! y-component
  !     else
  !       grid%n_xi(i,j,1) = zero
  !       grid%n_xi(i,j,2) = zero
  !     end if
      
  !     ! Calculate eta face (i-direction face)
  !     ! Vector from (i,j) to (i+1,j)
  !     dx = grid%x(i+1,j) - grid%x(i,j)
  !     dy = grid%y(i+1,j) - grid%y(i,j)
      
  !     ! Area magnitude (length of the face)
  !     grid%A_eta(i,j) = sqrt(dx**2 + dy**2)
      
  !     ! Normal vector (ensure proper normalization)
  !     if (grid%A_eta(i,j) > zero) then
  !       ! For a cartesian grid, normal to eta face points in y direction
  !       grid%n_eta(i,j,1) = -dy / grid%A_eta(i,j) ! x-component
  !       grid%n_eta(i,j,2) = dx / grid%A_eta(i,j)  ! y-component
  !     else
  !       grid%n_eta(i,j,1) = zero
  !       grid%n_eta(i,j,2) = zero
  !     end if
      
  !     ! Calculate cell centers
  !     grid%xc(i,j) = 0.25_prec * (grid%x(i,j) + grid%x(i+1,j) + &
  !                                grid%x(i,j+1) + grid%x(i+1,j+1))
  !     grid%yc(i,j) = 0.25_prec * (grid%y(i,j) + grid%y(i+1,j) + &
  !                                grid%y(i,j+1) + grid%y(i+1,j+1))
      
  !     ! Calculate cell volume (area for 2D)
  !     ! For cartesian grid, this is simply dx*dy
  !     ! grid%V(i,j) = (grid%x(i+1,j) - grid%x(i,j)) * (grid%y(i,j+1) - grid%y(i,j))
  !   end do
  ! end do
  
  
  ! do j = grid%jg_low, grid%jg_high
  !   do i = grid%ig_low, grid%ig_high
  !     ! Average xi face normals (if needed for your scheme)
  !     if (i > grid%ig_low) then
  !       grid%n_xi_avg(i,j,1) = 0.5_prec * (grid%n_xi(i,j,1) + grid%n_xi(i-1,j,1))
  !       grid%n_xi_avg(i,j,2) = 0.5_prec * (grid%n_xi(i,j,2) + grid%n_xi(i-1,j,2))
  !     else
  !       grid%n_xi_avg(i,j,1) = grid%n_xi(i,j,1)
  !       grid%n_xi_avg(i,j,2) = grid%n_xi(i,j,2)
  !     end if
      
  !     ! Average eta face normals (if needed for your scheme)
  !     if (j > grid%jg_low) then
  !       grid%n_eta_avg(i,j,1) = 0.5_prec * (grid%n_eta(i,j,1) + grid%n_eta(i,j-1,1))
  !       grid%n_eta_avg(i,j,2) = 0.5_prec * (grid%n_eta(i,j,2) + grid%n_eta(i,j-1,2))
  !     else
  !       grid%n_eta_avg(i,j,1) = grid%n_eta(i,j,1)
  !       grid%n_eta_avg(i,j,2) = grid%n_eta(i,j,2)
  !     end if
  !   end do
  ! end do

    ! Interior calculations
  do j = grid%jg_low, grid%jg_high
    do i = grid%ig_low, grid%ig_high
      ! Calculate xi face (j-direction face)
      ! Vector from (i,j) to (i,j+1)
      dx = grid%x(i,j+1) - grid%x(i,j)
      dy = grid%y(i,j+1) - grid%y(i,j)
      
      ! Area magnitude (length of the face)
      grid%A_xi(i,j) = sqrt(dx**2 + dy**2)
      
      ! Normal vector (ensure proper normalization)
      if (grid%A_xi(i,j) > zero) then
        grid%n_xi(i,j,1) = dy / grid%A_xi(i,j)  ! x-component
        grid%n_xi(i,j,2) = -dx / grid%A_xi(i,j) ! y-component
      else
        grid%n_xi(i,j,1) = zero
        grid%n_xi(i,j,2) = zero
      end if
      
      ! Calculate eta face (i-direction face)
      ! Vector from (i,j) to (i+1,j)
      dx = grid%x(i+1,j) - grid%x(i,j)
      dy = grid%y(i+1,j) - grid%y(i,j)
      
      ! Area magnitude (length of the face)
      grid%A_eta(i,j) = sqrt(dx**2 + dy**2)
      
      ! Normal vector (ensure proper normalization)
      if (grid%A_eta(i,j) > zero) then
        grid%n_eta(i,j,1) = -dy / grid%A_eta(i,j) ! x-component
        grid%n_eta(i,j,2) = dx / grid%A_eta(i,j)  ! y-component
      else
        grid%n_eta(i,j,1) = zero
        grid%n_eta(i,j,2) = zero
      end if
      
      ! Calculate cell volume (area for 2D) using shoelace formula
      area = 0.5_prec * abs( &
             (grid%x(i,j) * grid%y(i+1,j) - grid%y(i,j) * grid%x(i+1,j)) + &
             (grid%x(i+1,j) * grid%y(i+1,j+1) - grid%y(i+1,j) * grid%x(i+1,j+1)) + &
             (grid%x(i+1,j+1) * grid%y(i,j+1) - grid%y(i+1,j+1) * grid%x(i,j+1)) + &
             (grid%x(i,j+1) * grid%y(i,j) - grid%y(i,j+1) * grid%x(i,j)) )
      grid%V(i,j) = area
      
      ! Calculate cell centers using centroid formula
      cx = zero
      cy = zero
      cx = cx + (grid%x(i,j) + grid%x(i+1,j)) * &
                (grid%x(i,j) * grid%y(i+1,j) - grid%x(i+1,j) * grid%y(i,j))
      cx = cx + (grid%x(i+1,j) + grid%x(i+1,j+1)) * &
                (grid%x(i+1,j) * grid%y(i+1,j+1) - grid%x(i+1,j+1) * grid%y(i+1,j))
      cx = cx + (grid%x(i+1,j+1) + grid%x(i,j+1)) * &
                (grid%x(i+1,j+1) * grid%y(i,j+1) - grid%x(i,j+1) * grid%y(i+1,j+1))
      cx = cx + (grid%x(i,j+1) + grid%x(i,j)) * &
                (grid%x(i,j+1) * grid%y(i,j) - grid%x(i,j) * grid%y(i,j+1))
      cy = cy + (grid%y(i,j) + grid%y(i+1,j)) * &
                (grid%x(i,j) * grid%y(i+1,j) - grid%x(i+1,j) * grid%y(i,j))
      cy = cy + (grid%y(i+1,j) + grid%y(i+1,j+1)) * &
                (grid%x(i+1,j) * grid%y(i+1,j+1) - grid%x(i+1,j+1) * grid%y(i+1,j))
      cy = cy + (grid%y(i+1,j+1) + grid%y(i,j+1)) * &
                (grid%x(i+1,j+1) * grid%y(i,j+1) - grid%x(i,j+1) * grid%y(i+1,j+1))
      cy = cy + (grid%y(i,j+1) + grid%y(i,j)) * &
                (grid%x(i,j+1) * grid%y(i,j) - grid%x(i,j) * grid%y(i,j+1))
      grid%xc(i,j) = cx / (6.0_prec * area)
      grid%yc(i,j) = cy / (6.0_prec * area)
    end do
  end do
  
end subroutine cell_geometry

subroutine deallocate_grid( grid )
    
    type( grid_t ), intent( inout ) :: grid
    
    deallocate( grid%x, grid%y, grid%xc, grid%yc, grid%A_xi, grid%A_eta, &
                grid%n_xi, grid%n_eta, grid%n_xi_avg, grid%n_eta_avg, grid%V )
    
  end subroutine deallocate_grid
    
end module grid_type


module geometry

    use set_inputs
    use grid_type
    
    implicit none
    
contains

subroutine cartesian_grid(grid)
    use set_precision, only : prec
    use set_constants, only : one
    use set_inputs, only : imax, jmax, xmin, xmax, ymin, ymax, n_ghost
    use set_inputs, only : i_low, i_high, j_low, j_high
    use set_inputs, only : ig_low, ig_high, jg_low, jg_high
    integer :: i,j

    type(grid_t),intent(inout) ::grid

    call allocate_grid(grid)
    grid%i_low   = i_low
    grid%j_low   = j_low
    grid%i_high  = i_high
    grid%j_high  = j_high
    grid%ig_low  = ig_low
    grid%jg_low  = jg_low
    grid%ig_high = ig_high
    grid%jg_high = jg_high



    do j = j_low,j_high+1
        do i = i_low,i_high+1
          grid%x(i,j) = xmin + real(i-1,prec)/real(imax-1,prec)*(xmax-xmin)
          grid%y(i,j) = ymin + real(j-1,prec)/real(jmax-1,prec)*(ymax-ymin)
        end do
      end do
  
end subroutine cartesian_grid    
end module geometry

program main
  use set_precision, only : prec
  use set_constants, only : set_derived_constants
  use set_inputs,only : set_derived_inputs
  use grid_type
  use geometry
  
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
  
  do j = grid%j_low, grid%j_high
    do i = grid%i_low, grid%i_high
      print *, "Cell (", i, ",", j, "):"
      print *, "  A_xi  =", grid%A_xi(i,j)
      print *, "  A_eta =", grid%A_eta(i,j)
      print *, "  n_xi  = (", grid%n_xi(i,j,1), ",", grid%n_xi(i,j,2), ")"
      print *, "  n_eta = (", grid%n_eta(i,j,1), ",", grid%n_eta(i,j,2), ")"
      print *, "  Volume =", grid%V(i,j)
      print *, "  Center = (", grid%xc(i,j), ",", grid%yc(i,j), ")"
    end do
  end do
  
  ! Clean up
  call deallocate_grid(grid)
  
end program main



