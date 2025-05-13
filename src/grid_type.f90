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

! Allocate face areas (use i_low:i_high, j_low:j_high)
allocate(grid%A_xi(grid%i_low:grid%i_high, grid%j_low:grid%j_high-1))
allocate(grid%A_eta(grid%i_low:grid%i_high-1, grid%j_low:grid%j_high))

! Allocate normal vectors at faces
allocate(grid%n_xi(grid%i_low:grid%i_high, grid%j_low:grid%j_high-1, 2))
allocate(grid%n_eta(grid%i_low:grid%i_high-1, grid%j_low:grid%j_high, 2))

! Allocate averaged normal vectors
allocate(grid%n_xi_avg(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high, 2))
allocate(grid%n_eta_avg(grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high, 2))

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

! subroutine ghost_cells(grid)
!     type(grid_t), intent(inout) :: grid
!     integer :: i, j

!     ! bottom
!     do j = grid%j_low-1, grid%jg_low,-1
!       do i = grid%i_low, grid%i_high+1
!         grid%x(i,j) = two*grid%x(i,j+1) - grid%x(i,j+2)
!         grid%y(i,j) = two*grid%y(i,j+1) - grid%y(i,j+2)
!       end do
!     end do

!     ! top
!     do j = grid%j_high+1, grid%jg_high
!       do i = grid%i_low, grid%i_high+1
!         grid%x(i,j) = two*grid%x(i,j-1) - grid%x(i,j-2)
!         grid%y(i,j) = two*grid%y(i,j-1) - grid%y(i,j-2)
!       end do
!     end do

!     ! left
!     do j = grid%j_low, grid%j_high+1
!        do i = grid%i_low-1, grid%ig_low,-1
      
!         grid%x(i,j) = two*grid%x(i+1,j) - grid%x(i+2,j)
!         grid%y(i,j) = two*grid%y(i+1,j) - grid%y(i+2,j)
!       end do
!     end do

!     ! right
!     do j = grid%j_low, grid%j_high+1
!         do i = grid%i_high+1, grid%ig_high
!         grid%x(i,j) = two*grid%x(i-1,j) - grid%x(i-2,j)
!         grid%y(i,j) = two*grid%y(i-1,j) - grid%y(i-2,j)
!       end do
!     end do


!   end subroutine ghost_cells

subroutine cell_geometry(grid)
    type(grid_t), intent(inout) :: grid
    integer :: i, j,k
    real(prec) :: dx, dy
  
    ! Compute (approximate) cell centers (xc, yc) as the average of the four vertices
    ! do j = grid%j_cell_low, grid%j_cell_high
    !     do i = grid%i_cell_low, grid%i_cell_high
    !         grid%xc(i,j) = fourth * (grid%x(i,j) + grid%x(i+1,j) + grid%x(i,j+1) + grid%x(i+1,j+1))
    !         grid%yc(i,j) = fourth * (grid%y(i,j) + grid%y(i+1,j) + grid%y(i,j+1) + grid%y(i+1,j+1))
    !     end do
    ! end do
    real(prec) :: area, cx, cy, term
    real(kind=prec), dimension(4) :: x, y
    
    ! Compute cell centers (xc, yc) and volumes (V) using the shoelace formula
    do j = grid%j_cell_low, grid%j_cell_high
      do i = grid%i_cell_low, grid%i_cell_high
        ! Define the four vertices of the quadrilateral cell (counterclockwise)
        x(1) = grid%x(i,j)
        y(1) = grid%y(i,j)
        x(2) = grid%x(i+1,j)
        y(2) = grid%y(i+1,j)
        x(3) = grid%x(i+1,j+1)
        y(3) = grid%y(i+1,j+1)
        x(4) = grid%x(i,j+1)
        y(4) = grid%y(i,j+1)

        ! Initialize accumulators for area and centroid
        area = zero
        cx = zero
        cy = zero

        ! Compute signed area and centroid contributions using shoelace formula
        do k = 1, 3
          term = (x(k) * y(k+1) - x(k+1) * y(k))
          area = area + term
          cx = cx + (x(k) + x(k+1)) * term
          cy = cy + (y(k) + y(k+1)) * term
        end do
        term = (x(4) * y(1) - x(1) * y(4))
        area = area + term
        cx = cx + (x(4) + x(1)) * term
        cy = cy + (y(4) + y(1)) * term

        ! Finalize area (take absolute value for volume)
        area = half * area
        grid%V(i,j) = abs(area)

        ! Finalize centroid (divide by 6 * area)

          grid%xc(i,j) = cx / (6.0_prec * area)
          grid%yc(i,j) = cy / (6.0_prec * area)

      end do
    end do
    
    
    ! Compute face areas and normal vectors (xi-direction, vertical faces)
    do j = grid%j_low, grid%j_high-1
        do i = grid%i_low, grid%i_high
            dx = grid%x(i,j+1) - grid%x(i,j)
            dy = grid%y(i,j+1) - grid%y(i,j)
            grid%n_xi(i,j,1) = dy
            grid%n_xi(i,j,2) = -dx
            ! grid%A_xi(i,j) = sqrt(grid%n_xi(i,j,1)**2 + grid%n_xi(i,j,2)**2)
            grid%A_xi(i,j) = sqrt( sum( grid%n_xi(i,j,:)**2 ) )
            grid%n_xi(i,j,:) = grid%n_xi(i,j,:) / grid%A_xi(i,j)
        end do
        do i = grid%i_low, grid%i_high -1
            grid%n_xi_avg(i,j,:) = half * ( grid%n_xi(i,j,:) + grid%n_xi(i+1,j,:) )
            grid%n_xi_avg(i,j,:) = grid%n_xi_avg(i,j,:) / sqrt( sum( grid%n_xi_avg(i,j,:)**2 ) )
        end do
    end do
  
    ! Compute face areas and normal vectors for top boundary (xi-direction, j = grid%j_high)
    ! do i = grid%i_low, grid%i_high
    !     dx = grid%x(i,grid%j_high) - grid%x(i,grid%j_high-1)
    !     dy = grid%y(i,grid%j_high) - grid%y(i,grid%j_high-1)
    !     grid%n_xi(i,grid%j_high,1) = dy
    !     grid%n_xi(i,grid%j_high,2) = -dx
    !     grid%A_xi(i,grid%j_high) = sqrt(grid%n_xi(i,grid%j_high,1)**2 + grid%n_xi(i,grid%j_high,2)**2)
    !     grid%n_xi(i,grid%j_high,1) = grid%n_xi(i,grid%j_high,1) / grid%A_xi(i,grid%j_high)
    !     grid%n_xi(i,grid%j_high,2) = grid%n_xi(i,grid%j_high,2) / grid%A_xi(i,grid%j_high)
    !     grid%n_xi_avg(i,grid%j_high,:) = grid%n_xi(i,grid%j_high,:)
    ! end do
  
    ! Compute face areas and normal vectors (eta-direction, horizontal faces)
    do j = grid%j_low, grid%j_high
        do i = grid%i_low, grid%i_high-1
            dx = grid%x(i+1,j) - grid%x(i,j)
            dy = grid%y(i+1,j) - grid%y(i,j)
            grid%n_eta(i,j,1) = dy
            grid%n_eta(i,j,2) = -dx
            ! grid%A_eta(i,j) = sqrt(grid%n_eta(i,j,1)**2 + grid%n_eta(i,j,2)**2)
            grid%A_eta(i,j) = sqrt( sum( grid%n_eta(i,j,:)**2 ) )
            ! grid%n_eta(i,j,1) = -grid%n_eta(i,j,1) / grid%A_eta(i,j)
            ! grid%n_eta(i,j,2) =  grid%n_eta(i,j,2) / grid%A_eta(i,j)
            grid%n_eta(i,j,:) =  -grid%n_eta(i,j,:) / grid%A_eta(i,j)
            ! grid%n_eta(i,j,:) =  grid%n_eta(i,j,:) / grid%A_eta(i,j)
        end do        
    end do

    do j = grid%j_low, grid%j_high-1
        do i = grid%i_low, grid%i_high-1
            grid%n_eta_avg(i,j,:) = half * ( grid%n_eta(i,j,:) + grid%n_eta(i,j+1,:) )
            grid%n_eta_avg(i,j,:) = grid%n_eta_avg(i,j,:) / sqrt( sum( grid%n_eta_avg(i,j,:)**2 ) )
        end do
    end do

  
    ! Compute face areas and normal vectors for right boundary (eta-direction, i = grid%i_high)
    ! do j = grid%j_low, grid%j_high
    !     dx = grid%x(grid%i_high,j) - grid%x(grid%i_high-1,j)
    !     dy = grid%y(grid%i_high,j) - grid%y(grid%i_high-1,j)
    !     grid%n_eta(grid%i_high,j,1) = dy
    !     grid%n_eta(grid%i_high,j,2) = -dx
    !     grid%A_eta(grid%i_high,j) = sqrt(grid%n_eta(grid%i_high,j,1)**2 + grid%n_eta(grid%i_high,j,2)**2)
    !     if (grid%A_eta(grid%i_high,j) > zero) then
    !         grid%n_eta(grid%i_high,j,1) = -grid%n_eta(grid%i_high,j,1) / grid%A_eta(grid%i_high,j)
    !         grid%n_eta(grid%i_high,j,2) = grid%n_eta(grid%i_high,j,2) / grid%A_eta(grid%i_high,j)
    !     else
    !         grid%n_eta(grid%i_high,j,1) = zero
    !         grid%n_eta(grid%i_high,j,2) = zero
    !     end if
    !     grid%n_eta_avg(grid%i_high,j,:) = grid%n_eta(grid%i_high,j,:)
    ! end do
  
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
              grid%n_xi, grid%n_eta,  grid%V,grid%n_xi_avg,grid%n_eta_avg )
  
end subroutine deallocate_grid
  
end module grid_type