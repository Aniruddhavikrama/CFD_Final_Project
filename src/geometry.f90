module geometry

    use set_inputs
    
    implicit none
    
contains

subroutine cartesian_grid
    use set_precision, only : prec
    use set_constants, only : one
    use set_inputs, only : imax, jmax, xmin, xmax, ymin, ymax, n_ghost
    use set_inputs, only : i_low, i_high, j_low, j_high
    use set_inputs, only : ig_low, ig_high, jg_low, jg_high



    do j = j_low,j_high+1
        do i = i_low,i_high+1
          grid%x(i,j) = xmin + real(i-1,prec)/real(imax-1,prec)*(xmax-xmin)
          grid%y(i,j) = ymin + real(j-1,prec)/real(jmax-1,prec)*(ymax-ymin)
        end do
      end do
  
end subroutine cartesian_grid    
end module geometry