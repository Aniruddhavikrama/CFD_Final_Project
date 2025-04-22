module grid_type
    use set_precision, only :prec
    use set_constants, only : zero, one, two, half
    use set_inputs,    only : i_low, i_high, ig_low, ig_high
    use set_inputs,    only : j_low, j_high, jg_low, jg_high
    
    implicit none
    

    private
  
    public :: grid_t, allocate_grid, deallocate_grid

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
    type(grid_t), intent(inout) :: grid
    real(prec):: dx,dy,area,cx,cy
    integer:: i,j

    ! do j = grid%jg_low, grid%jg_high
    !     do i = grid%ig_low, grid%ig_high
    !       ! Calculate xi face (j-direction face)
    !       ! Vector from (i,j) to (i,j+1)
    !       dx = grid%x(i,j+1) - grid%x(i,j)
    !       dy = grid%y(i,j+1) - grid%y(i,j)
          
    !       ! Area magnitude (length of the face)
    !       grid%A_xi(i,j) = sqrt(dx**2 + dy**2)
          
    !       ! Normal vector (ensure proper normalization)
    !       if (grid%A_xi(i,j) > zero) then
    !         ! For a cartesian grid, normal to xi face points in x direction
    !         grid%n_xi(i,j,1) = dy / grid%A_xi(i,j)  ! x-component
    !         grid%n_xi(i,j,2) = -dx / grid%A_xi(i,j) ! y-component
    !       else
    !         grid%n_xi(i,j,1) = zero
    !         grid%n_xi(i,j,2) = zero
    !       end if
          
    !       ! Calculate eta face (i-direction face)
    !       ! Vector from (i,j) to (i+1,j)
    !       dx = grid%x(i+1,j) - grid%x(i,j)
    !       dy = grid%y(i+1,j) - grid%y(i,j)
          
    !       ! Area magnitude (length of the face)
    !       grid%A_eta(i,j) = sqrt(dx**2 + dy**2)
          
    !       ! Normal vector (ensure proper normalization)
    !       if (grid%A_eta(i,j) > zero) then
    !         ! For a cartesian grid, normal to eta face points in y direction
    !         grid%n_eta(i,j,1) = -dy / grid%A_eta(i,j) ! x-component
    !         grid%n_eta(i,j,2) = dx / grid%A_eta(i,j)  ! y-component
    !       else
    !         grid%n_eta(i,j,1) = zero
    !         grid%n_eta(i,j,2) = zero
    !       end if
          
    !       ! Calculate cell centers
    !       grid%xc(i,j) = 0.25_prec * (grid%x(i,j) + grid%x(i+1,j) + &
    !                                  grid%x(i,j+1) + grid%x(i+1,j+1))
    !       grid%yc(i,j) = 0.25_prec * (grid%y(i,j) + grid%y(i+1,j) + &
    !                                  grid%y(i,j+1) + grid%y(i+1,j+1))
          
    !       ! Calculate cell volume (area for 2D)
    !       ! For cartesian grid, this is simply dx*dy
    !       ! grid%V(i,j) = (grid%x(i+1,j) - grid%x(i,j)) * (grid%y(i,j+1) - grid%y(i,j))
    !     end do
    !   end do
      
      
    !   do j = grid%jg_low, grid%jg_high
    !     do i = grid%ig_low, grid%ig_high
    !       ! Average xi face normals (if needed for your scheme)
    !       if (i > grid%ig_low) then
    !         grid%n_xi_avg(i,j,1) = 0.5_prec * (grid%n_xi(i,j,1) + grid%n_xi(i-1,j,1))
    !         grid%n_xi_avg(i,j,2) = 0.5_prec * (grid%n_xi(i,j,2) + grid%n_xi(i-1,j,2))
    !       else
    !         grid%n_xi_avg(i,j,1) = grid%n_xi(i,j,1)
    !         grid%n_xi_avg(i,j,2) = grid%n_xi(i,j,2)
    !       end if
          
    !       ! Average eta face normals (if needed for your scheme)
    !       if (j > grid%jg_low) then
    !         grid%n_eta_avg(i,j,1) = 0.5_prec * (grid%n_eta(i,j,1) + grid%n_eta(i,j-1,1))
    !         grid%n_eta_avg(i,j,2) = 0.5_prec * (grid%n_eta(i,j,2) + grid%n_eta(i,j-1,2))
    !       else
    !         grid%n_eta_avg(i,j,1) = grid%n_eta(i,j,1)
    !         grid%n_eta_avg(i,j,2) = grid%n_eta(i,j,2)
    !       end if

    !       call volume_calc( (/ grid%x(i,j), grid%x(i+1,j), grid%x(i+1,j+1), &
    !                             grid%x(i,j+1),grid%x(i,j) /), &
    !                                 (/ grid%y(i,j), grid%y(i+1,j), grid%y(i+1,j+1), &
    !                                         grid%y(i,j+1),grid%y(i,j) /), &
    !                                 grid%xc(i,j), grid%yc(i,j), grid%V(i,j) )
    !     end do
    !   end do
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

!  Boundary face calculations for ghost cells 
! do i = grid%ig_low, grid%ig_high
!   j = grid%jg_high + 1
!   dx = grid%x(i+1,j) - grid%x(i,j)
!   dy = grid%y(i+1,j) - grid%y(i,j)
!   grid%A_eta(i,j) = sqrt(dx**2 + dy**2)
!   if (grid%A_eta(i,j) > zero) then
!     grid%n_eta(i,j,1) = -dy / grid%A_eta(i,j)
!     grid%n_eta(i,j,2) = dx / grid%A_eta(i,j)
!   else
!     grid%n_eta(i,j,1) = zero
!     grid%n_eta(i,j,2) = zero
!   end if
! end do
!
! do j = grid%jg_low, grid%jg_high
!   i = grid%ig_high + 1
!   dx = grid%x(i,j+1) - grid%x(i,j)
!   dy = grid%y(i,j+1) - grid%y(i,j)
!   grid%A_xi(i,j) = sqrt(dx**2 + dy**2)
!   if (grid%A_xi(i,j) > zero) then
!     grid%n_xi(i,j,1) = dy / grid%A_xi(i,j)
!     grid%n_xi(i,j,2) = -dx / grid%A_xi(i,j)
!   else
!     grid%n_xi(i,j,1) = zero
!     grid%n_xi(i,j,2) = zero
!   end if
! end do



end subroutine cell_geometry



subroutine deallocate_grid( grid )
    
    type( grid_t ), intent( inout ) :: grid
    
    deallocate( grid%x, grid%y, grid%xc, grid%yc, grid%A_xi, grid%A_eta, &
                grid%n_xi, grid%n_eta, grid%n_xi_avg, grid%n_eta_avg, grid%V )
    
  end subroutine deallocate_grid
    
end module grid_type