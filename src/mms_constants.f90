module mms_constants
  

    use set_precision,only : prec   
  
    implicit none
  
   ! NOTE: These are currently set up to run the supersonic manufactured solution
    real(prec), parameter :: rho0   = 1.0_prec
    real(prec), parameter :: rhox   = 0.15_prec
    real(prec), parameter :: rhoy   = -0.1_prec
    real(prec), parameter :: uvel0  = 800.0_prec
    real(prec), parameter :: uvelx  = 50.0_prec
    real(prec), parameter :: uvely  = -30.0_prec
    real(prec), parameter :: vvel0  = 800.0_prec
    real(prec), parameter :: vvelx  = -75.0_prec
    real(prec), parameter :: vvely  = 40.0_prec
    real(prec), parameter :: wvel0  = 0.0_prec
    real(prec), parameter :: wvelx  = 0.0_prec
    real(prec), parameter :: wvely  = 0.0_prec
    real(prec), parameter :: press0 = 100000.0_prec
    real(prec), parameter :: pressx = 20000.0_prec
    real(prec), parameter :: pressy = 50000.0_prec


    !subsonic manufactued solution(should update the values from the paper)
    ! real(prec), parameter :: rho0   = 1.0_prec
    ! real(prec), parameter :: rhox   = 0.15_prec
    ! real(prec), parameter :: rhoy   = -0.1_prec
    ! real(prec), parameter :: uvel0  = 800.0_prec
    ! real(prec), parameter :: uvelx  = 50.0_prec
    ! real(prec), parameter :: uvely  = -30.0_prec
    ! real(prec), parameter :: vvel0  = 800.0_prec
    ! real(prec), parameter :: vvelx  = -75.0_prec
    ! real(prec), parameter :: vvely  = 40.0_prec
    ! real(prec), parameter :: wvel0  = 0.0_prec
    ! real(prec), parameter :: wvelx  = 0.0_prec
    ! real(prec), parameter :: wvely  = 0.0_prec
    ! real(prec), parameter :: press0 = 100000.0_prec
    ! real(prec), parameter :: pressx = 20000.0_prec
    ! real(prec), parameter :: pressy = 50000.0_prec
    
  end module mms_constants