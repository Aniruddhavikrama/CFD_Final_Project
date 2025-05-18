module residual_io
  use set_precision, only : prec
  implicit none
contains

  subroutine write_residuals(iter, Rnorm)
    integer,       intent(in)    :: iter
    real(prec),    intent(in)    :: Rnorm(4)
    logical, save :: first = .true.
    integer        :: unit

    if (first) then
      open(newunit=unit,                              &
           file   = "residual_history.dat",           &
           status = "replace",                       &
           action = "write",                         &
           form   = "formatted")
      write(unit,'(A)') '#  iter     Res1           Res2           Res3           Res4'
      first = .false.
    else
      open(newunit=unit,                              &
           file   = "residual_history.dat",           &
           status = "old",                           &
           action = "write",                         &
           position = "append",                      &
           form   = "formatted")
    end if

    write(unit,'(I8,4E15.7)') iter,                    &
         Rnorm(1), Rnorm(2), Rnorm(3), Rnorm(4)

    close(unit)
  end subroutine write_residuals

end module residual_io
