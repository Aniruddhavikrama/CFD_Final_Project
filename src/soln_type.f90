module soln_type
    use set_precision,only : prec
    use set_inputs
    use fluid_constants,only: gamma
    use set_constants,only : zero,one,two,half
    use grid_type,only : grid_t
    
    
    implicit none
    private
    public :: soln_t,allocate_soln,deallocate_soln

    type soln_t
    

    real(prec),allocatable,dimension(:,:,:) :: U !Conserved
    real(prec),allocatable,dimension(:,:,:) :: V !primitive
    real(prec),allocatable,dimension(:,:,:) :: S !Source terms
    real(prec),allocatable,dimension(:,:,:) :: R ! Residual
    ! real(prec),allocatable,dimension(:) :: rinit
    real(prec), allocatable, dimension(:,:)   :: asnd
    ! real(prec), allocatable, dimension(:,:)   :: mach
    ! real(prec), allocatable, dimension(:,:)   :: temp
    real(prec), allocatable, dimension(:,:)   :: dt
    ! real(prec), allocatable, dimension(:,:)     :: DEnorm
    ! ! real(prec), allocatable, dimension(:)     :: rnorm
    ! real(prec), allocatable, dimension(:,:,:) :: DE ! discretization error
    real(prec), allocatable, dimension(:,:,:) :: Fxi  ! normal fluxes
    real(prec), allocatable, dimension(:,:,:) :: Feta ! normal fluxes
    ! real(prec), allocatable, dimension(:,:,:) :: psi_p_xi  ! limiters
    ! real(prec), allocatable, dimension(:,:,:) :: psi_p_eta  ! limiters
    ! real(prec), allocatable, dimension(:,:,:) :: psi_m_xi  ! limiters
    ! real(prec), allocatable, dimension(:,:,:) :: psi_m_eta  ! limiters
    real(prec), allocatable, dimension(:,:,:) :: Umms ! MMS conserved variables
    real(prec), allocatable, dimension(:,:,:) :: Vmms ! MMS primitive variables
    real(prec), allocatable, dimension(:,:,:) :: Smms ! MMS source terms
    
    

    end type soln_t

contains

subroutine allocate_soln(soln,grid)
    type(soln_t),intent(inout) :: soln
    type(grid_t), intent(in) :: grid

    allocate(soln%U( neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high ), &
            soln%V( neq, grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high ), &
        
            ! soln%U( neq, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high ), &
            ! soln%V( neq, grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high ), &
            soln%S( neq, ig_low:ig_high, jg_low:jg_high ), &
            soln%R(  neq,grid%i_low:grid%i_high,   grid%j_low:grid%j_high )  , &                           
            soln%asnd( ig_low:ig_high, jg_low:jg_high ),   &
            ! soln%mach( ig_low:ig_high, jg_low:jg_high ),   &
            ! soln%temp( ig_low:ig_high, jg_low:jg_high ),   &
            soln%dt(   i_cell_low:i_cell_high, j_cell_low:j_cell_high ),   &
            soln%Fxi(  neq, grid%i_low:grid%i_high, grid%j_low:grid%j_high ), &
            soln%Feta( neq, grid%i_low:grid%i_high, grid%j_low:grid%j_high ), &
            ! soln%psi_p_xi(  neq, ig_low-1:ig_high, j_low:j_high ), &
            ! soln%psi_m_xi(  neq, ig_low-1:ig_high, j_low:j_high ), &
            ! soln%psi_p_eta( neq, i_low:i_high, jg_low-1:jg_high ), &
            ! soln%psi_m_eta( neq, i_low:i_high, jg_low-1:jg_high ) , &
            ! soln%rinit( neq )
            soln%Vmms( neq, grid%ig_low:grid%ig_high,  grid%jg_low:grid%jg_high ),&
            soln%Umms( neq, grid%ig_low:grid%ig_high,  grid%jg_low:grid%jg_high ),&
            soln%Smms( neq, grid%ig_low:grid%ig_high,  grid%jg_low:grid%jg_high ),&  
            )


            soln%U     = zero
            soln%Fxi   = zero
            soln%Feta  = zero
            ! soln%S     = zero
            soln%V     = zero
            soln%R     = one
            soln%asnd  = one
            ! soln%mach  = one
            ! soln%temp  = one
            soln%dt    = one
            ! soln%rnorm = one
            ! soln%rinit = one
            
            ! soln%psi_p_xi  = one
            ! soln%psi_m_xi  = one
            ! soln%psi_p_eta = one
            ! soln%psi_m_eta = one
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
                    ! soln%mach,  &
                    ! soln%temp,  &
                    soln%dt,    &
                    ! soln%rnorm, &
                    ! soln%rinit,&
                    ! soln%psi_p_xi,&
                    ! soln%psi_m_xi, &
                    ! soln%psi_p_eta,&
                    ! soln%psi_m_eta 
                    soln%Umms ,&
                    soln%Vmms, &
                    soln%Smms ) 
                end subroutine deallocate_soln    
    
end module soln_type