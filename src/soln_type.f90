module soln_type
    use set_precision,only : prec
    use set_inputs,only: i_low,i_high,ig_low,ig_high,j_low,j_high,jg_low,jg_high
    use fluid_constants,only: gamma
    use set_constants,only : zero,one,two,half
    use set_inputs,only : i_low,i_high,ig_high,ig_low
    use set_inputs,only : j_low,j_high,jg_high,jg_low
    use set_inputs,only : neq
    
    
    implicit none
    private
    public :: soln_t,allocate_soln,deallocate_soln

    type soln_t

    real(prec),allocatable,dimension(:,:,:) :: U !Conserved
    real(prec),allocatable,dimension(:,:,:) :: V !primitive
    real(prec),allocatable,dimension(:,:,:) :: S !Source terms
    real(prec),allocatable,dimension(:,:,:) :: R ! Residual
    real(prec),allocatable,dimension(:) :: rinit
    real(prec), allocatable, dimension(:,:)   :: asnd
    real(prec), allocatable, dimension(:,:)   :: mach
    real(prec), allocatable, dimension(:,:)   :: temp
    real(prec), allocatable, dimension(:,:)   :: dt
    real(prec), allocatable, dimension(:,:)     :: DEnorm
    real(prec), allocatable, dimension(:)     :: rnorm
    real(prec), allocatable, dimension(:,:,:) :: DE ! discretization error
    real(prec), allocatable, dimension(:,:,:) :: Fxi  ! normal fluxes
    real(prec), allocatable, dimension(:,:,:) :: Feta ! normal fluxes
    real(prec), allocatable, dimension(:,:,:) :: psi_p_xi  ! limiters
    real(prec), allocatable, dimension(:,:,:) :: psi_p_eta  ! limiters
    real(prec), allocatable, dimension(:,:,:) :: psi_m_xi  ! limiters
    real(prec), allocatable, dimension(:,:,:) :: psi_m_eta  ! limiters
    

    end type soln_t

contains

subroutine allocate_soln(soln)
    type(soln_t),intent(inout) :: soln

    allocate(soln%U( neq, ig_low:ig_high, jg_low:jg_high ), &
            soln%V( neq, ig_low:ig_high, jg_low:jg_high ), &
            soln%S( neq, ig_low:ig_high, jg_low:jg_high ), &
            soln%R(  neq,i_low:i_high,   j_low:j_high )  , &                           
            soln%asnd( ig_low:ig_high, jg_low:jg_high ),   &
            soln%mach( ig_low:ig_high, jg_low:jg_high ),   &
            soln%temp( ig_low:ig_high, jg_low:jg_high ),   &
            soln%dt(   ig_low:ig_high, jg_low:jg_high ),   &
            soln%Fxi(  neq, i_low-1:i_high, j_low:j_high ), &
            soln%Feta( neq, i_low:i_high, j_low-1:j_high ), &
            soln%psi_p_xi(  neq, ig_low-1:ig_high, j_low:j_high ), &
            soln%psi_m_xi(  neq, ig_low-1:ig_high, j_low:j_high ), &
            soln%psi_p_eta( neq, i_low:i_high, jg_low-1:jg_high ), &
            soln%psi_m_eta( neq, i_low:i_high, jg_low-1:jg_high ) , &
            soln%rinit( neq )  )


            soln%U     = zero
            soln%Fxi   = zero
            soln%Feta  = zero
            soln%S     = zero
            soln%V     = zero
            soln%R     = one
            soln%asnd  = one
            soln%mach  = one
            soln%temp  = one
            soln%dt    = one
            soln%rnorm = one
            soln%rinit = one
            
            soln%psi_p_xi  = one
            soln%psi_m_xi  = one
            soln%psi_p_eta = one
            soln%psi_m_eta = one
end subroutine allocate_soln

subroutine deallocate_soln(soln)
    implicit none

    type(soln_t),intent(inout) :: soln

    deallocate(     soln%U,     &
                    soln%Fxi,   &
                    soln%Feta,  &
                    soln%S,     &
                    soln%V,     &
                    soln%R,     &
                    soln%asnd,  &
                    soln%mach,  &
                    soln%temp,  &
                    soln%dt,    &
                    soln%rnorm, &
                    soln%rinit,&
                    soln%psi_p_xi,&
                    soln%psi_m_xi, &
                    soln%psi_p_eta,&
                    soln%psi_m_eta ) 
                end subroutine deallocate_soln    
    
end module soln_type