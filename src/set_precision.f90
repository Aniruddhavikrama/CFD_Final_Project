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