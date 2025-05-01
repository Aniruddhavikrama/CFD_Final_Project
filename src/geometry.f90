module geometry

  use set_inputs
  use grid_type
  use soln_type
  
  implicit none
  private

  public :: cartesian_grid,write_tecplot_file,read_grid,write_solution_dat
  
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


subroutine read_grid(grid, grid_file, ierr)
  type(grid_t),    intent(inout) :: grid
  character(*),    intent(in)    :: grid_file
  integer,         intent(out)   :: ierr

  integer :: ios, unit_grid
  integer :: nblocks, blockID
  integer :: ni, nj
  integer :: i, j

  ierr = 0

  !––– open with a free unit number –––
  open(newunit=unit_grid, file=trim(grid_file), status='old', form='formatted', &
       iostat=ios)
  if (ios /= 0) then
    ierr = ios
    print *, '*** Error: cannot open grid file ', trim(grid_file)
    return
  end if

  !––– skip the block‐count line –––
  read(unit_grid, *, iostat=ios) nblocks
  if (ios /= 0) then
    ierr = ios
    print *, '*** Error: reading block count from ', trim(grid_file)
    close(unit_grid)
    return
  end if

  !––– read ni, nj, (and discard any 3rd integer) –––
  read(unit_grid, *, iostat=ios) ni, nj, blockID
  if (ios /= 0) then
    ierr = ios
    print *, '*** Error: reading dims from ', trim(grid_file)
    close(unit_grid)
    return
  end if

  if (ni < 2 .or. nj < 2) then
    ierr = 1
    print *, '*** Error: invalid grid dimensions: ni=', ni, ' nj=', nj
    close(unit_grid)
    return
  end if

  !––– convert to cell‐counts and set all our index bounds –––
  grid%imax = ni - 1
  grid%jmax = nj - 1

  grid%i_cell_low  = 1
  grid%i_cell_high = grid%imax
  grid%j_cell_low  = 1
  grid%j_cell_high = grid%jmax

  grid%i_low  = grid%i_cell_low
  grid%i_high = grid%i_cell_high + 1
  grid%j_low  = grid%j_cell_low
  grid%j_high = grid%j_cell_high + 1

  grid%ig_low  = grid%i_cell_low - n_ghost
  grid%ig_high = grid%i_cell_high + n_ghost
  grid%jg_low  = grid%j_cell_low - n_ghost
  grid%jg_high = grid%j_cell_high + n_ghost

  !––– allocate all arrays (vertices include ghosts) –––
  call allocate_grid(grid)

  !––– read X‐coordinates into the 1:ni × 1:nj interior –––
  read(unit_grid, *, iostat=ios) (( grid%x(i,j), i = 1, ni ), j = 1, nj)
  if (ios /= 0) then
    ierr = ios
    print *, '*** Error: reading X coords from ', trim(grid_file)
    call deallocate_grid(grid)
    close(unit_grid)
    return
  end if

  !––– read Y‐coordinates –––
  read(unit_grid, *, iostat=ios) (( grid%y(i,j), i = 1, ni ), j = 1, nj)
  if (ios /= 0) then
    ierr = ios
    print *, '*** Error: reading Y coords from ', trim(grid_file)
    call deallocate_grid(grid)
    close(unit_grid)
    return
  end if

  close(unit_grid)
end subroutine read_grid


subroutine write_tecplot_file_old(grid, filename)
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
end subroutine write_tecplot_file_old


subroutine write_tecplot_file(grid, filename)
  use set_precision, only : prec
  use set_constants, only : zero
  type(grid_t), intent(in) :: grid
  character(len=*), intent(in) :: filename
  integer :: i, j, unit, cnt
  real(prec) :: a_xi, n_xi_x, n_xi_y, a_eta, n_eta_x, n_eta_y
  
  ! Open the file
  open(newunit=unit, file=trim(filename), status='replace', action='write')
  
  ! Write the Tecplot header
  write(unit, '(A)') 'TITLE = "2D Cartesian Grid with Face Areas and Normals"'
  write(unit, '(A)') 'VARIABLES = "X", "Y", "A_xi", "n_xi_x", "n_xi_y", "A_eta", "n_eta_x", "n_eta_y"'
  write(unit, '(A,I0,A,I0,A)') 'ZONE I=', grid%i_high - grid%i_low + 1, &
                               ', J=', grid%j_high - grid%j_low + 1, &
                               ', DATAPACKING=BLOCK'
  ! write(unit, '(A)') 'Varlocation=([3-8]=CELLCENTERED)'
  write(unit, '(A)') 'DT=(DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE)'
  
  ! Write nodal data

  write(unit,'(E24.16)') ( ( grid%x(i,j), i = grid%i_low, grid%i_high ), &
                                          j = grid%j_low, grid%j_high )
  write(unit,'(E24.16)') ( ( grid%y(i,j), i = grid%i_low, grid%i_high ), &
                                          j = grid%j_low, grid%j_high )
  ! write cell-centered data
  ! write(unit,'(E24.16)') ( ( grid%A_xi(i,j), i = grid%i_low, grid%i_high-1 ), &
  !                                            j = grid%j_low, grid%j_high-1 )
  ! write(unit,'(E24.16)') ( ( grid%n_xi(i,j,1), i = grid%i_low, grid%i_high-1 ), &
  !                                            j = grid%j_low, grid%j_high-1 )
  ! write(unit,'(E24.16)') ( ( grid%n_xi(i,j,2), i = grid%i_low, grid%i_high-1 ), &
  !                                            j = grid%j_low, grid%j_high-1 )
  ! write(unit,'(E24.16)') ( ( grid%A_eta(i,j), i = grid%i_low, grid%i_high-1 ), &
  !                                            j = grid%j_low, grid%j_high-1 )
  ! write(unit,'(E24.16)') ( ( grid%n_eta(i,j,1), i = grid%i_low, grid%i_high-1 ), &
  !                                            j = grid%j_low, grid%j_high-1 )
  ! write(unit,'(E24.16)') ( ( grid%n_eta(i,j,2), i = grid%i_low, grid%i_high-1 ), &
  !                                            j = grid%j_low, grid%j_high-1 )

  cnt = 0
  do j = grid%j_low, grid%j_high
    do i = grid%i_low, grid%i_high
      if (j==grid%j_high) then
        write(unit,'(E24.16)') zero
      else
        write(unit,'(E24.16)') grid%A_xi(i,j)
      end if
    end do
  end do

  do j = grid%j_low, grid%j_high
    do i = grid%i_low, grid%i_high
      if (j==grid%j_high) then
        write(unit,'(E24.16)') zero
      else
        write(unit,'(E24.16)') grid%n_xi(i,j,1)
      end if
    end do
  end do

  do j = grid%j_low, grid%j_high
    do i = grid%i_low, grid%i_high
      if (j==grid%j_high) then
        write(unit,'(E24.16)') zero
      else
          write(unit,'(E24.16)') grid%n_xi(i,j,2)
      end if
    end do
  end do



  do j = grid%j_low, grid%j_high
    do i = grid%i_low, grid%i_high
      if (i==grid%i_high) then
        write(unit,'(E24.16)') zero
      else
        write(unit,'(E24.16)') grid%A_eta(i,j)
      end if
    end do
  end do

  do j = grid%j_low, grid%j_high
    do i = grid%i_low, grid%i_high
      if (i==grid%i_high) then
        write(unit,'(E24.16)') zero
      else
        write(unit,'(E24.16)') grid%n_eta(i,j,1)
      end if
    end do
  end do

  do j = grid%j_low, grid%j_high
    do i = grid%i_low, grid%i_high
      if (i==grid%i_high) then
        write(unit,'(E24.16)') zero
      else
        write(unit,'(E24.16)') grid%n_eta(i,j,2)
      end if
    end do
  end do

  
  close(unit)
  end subroutine write_tecplot_file

  subroutine write_solution_dat(grid, soln, filename)
    type(grid_t),       intent(in)    :: grid
    type(soln_t),       intent(in)    :: soln
    character(len=*),   intent(in)    :: filename

    integer, parameter :: io_unit = 20
    integer :: i,j,unit

    !—Open (old file replaced automatically)—
    open(newunit=unit, file=trim(filename), status='replace', action='write')

    !—Tecplot header—
    write(io_unit,'(A)') 'TITLE = "MMS Solution"'
    write(io_unit,'(A)') 'VARIABLES = "Xcenter","Ycenter","rho","u","v","p","rho","rho_u","rho_v","E"'
    write(io_unit,'(A,I0,A,I0,A)') 'ZONE I=', &
         grid%i_cell_high - grid%i_cell_low + 1, &
         ', J=', grid%j_cell_high - grid%j_cell_low + 1, &
         ', DATAPACKING=BLOCK'
    write(io_unit,'(A)') 'DT=(DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE)'

    !—Block-writes of each field at cell centers—
    write(io_unit,'(E24.16)') ((grid%xc(i,j), i=grid%i_cell_low,grid%i_cell_high), &
                                j=grid%j_cell_low,grid%j_cell_high)
    write(io_unit,'(E24.16)') ((grid%yc(i,j), i=grid%i_cell_low,grid%i_cell_high), &
                                j=grid%j_cell_low,grid%j_cell_high)

    write(io_unit,'(E24.16)') ((soln%V(1,i,j), i=grid%i_cell_low,grid%i_cell_high), &
                                j=grid%j_cell_low,grid%j_cell_high)
    write(io_unit,'(E24.16)') ((soln%V(2,i,j), i=grid%i_cell_low,grid%i_cell_high), &
                                j=grid%j_cell_low,grid%j_cell_high)
    write(io_unit,'(E24.16)') ((soln%V(3,i,j), i=grid%i_cell_low,grid%i_cell_high), &
                                j=grid%j_cell_low,grid%j_cell_high)
    write(io_unit,'(E24.16)') ((soln%V(4,i,j), i=grid%i_cell_low,grid%i_cell_high), &
                                j=grid%j_cell_low,grid%j_cell_high)

    write(io_unit,'(E24.16)') ((soln%U(1,i,j), i=grid%i_cell_low,grid%i_cell_high), &
                                j=grid%j_cell_low,grid%j_cell_high)
    write(io_unit,'(E24.16)') ((soln%U(2,i,j), i=grid%i_cell_low,grid%i_cell_high), &
                                j=grid%j_cell_low,grid%j_cell_high)
    write(io_unit,'(E24.16)') ((soln%U(3,i,j), i=grid%i_cell_low,grid%i_cell_high), &
                                j=grid%j_cell_low,grid%j_cell_high)
    write(io_unit,'(E24.16)') ((soln%U(4,i,j), i=grid%i_cell_low,grid%i_cell_high), &
                                j=grid%j_cell_low,grid%j_cell_high)

    close(io_unit)
  end subroutine write_solution_dat


end module geometry