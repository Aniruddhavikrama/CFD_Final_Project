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
    real(prec),allocatable,dimension(:) :: rinit
    real(prec), allocatable, dimension(:,:)   :: asnd
    ! real(prec), allocatable, dimension(:,:,:) :: L ! eigenvalues
    real(prec), allocatable, dimension(:,:)   :: mach
    real(prec), allocatable, dimension(:,:)   :: temp
    real(prec), allocatable, dimension(:,:)   :: dt
    ! real(prec), allocatable, dimension(:,:)     :: DEnorm
    real(prec), allocatable, dimension(:)     :: rnorm
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
            ! soln%S( neq, ig_low:ig_high, jg_low:jg_high ), &
            soln%S( neq,grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high ), &
            soln%R(  neq,grid%i_cell_low:grid%i_cell_high,   grid%j_cell_low:grid%j_cell_high )  , &                           
            soln%asnd( grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high ),   &
            ! soln%L( neq,  grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high ), &
            soln%mach( grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high ),   &
            soln%temp( grid%ig_low:grid%ig_high, grid%jg_low:grid%jg_high ),   &
            soln%dt(   grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high ),   &
            soln%Fxi(  neq, grid%i_low:grid%i_high, grid%j_low:grid%j_high-1 ), &
            soln%Feta( neq, grid%i_low:grid%i_high-1, grid%j_low:grid%j_high ), &
            ! soln%psi_p_xi(  neq, grid%i_low:grid%i_high,   grid%j_low:grid%j_high ), &
            ! soln%psi_m_xi(  neq, grid%i_low:grid%i_high,   grid%j_low:grid%j_high), &
            ! soln%psi_p_eta( neq, grid%i_low:grid%i_high,   grid%j_low:grid%j_high ), &
            ! soln%psi_m_eta( neq, grid%i_low:grid%i_high,   grid%j_low:grid%j_high ) , &
            soln%rinit( neq ) ,&
            soln%rnorm(neq) ,&
            soln%Vmms( neq, grid%ig_low:grid%ig_high,  grid%jg_low:grid%jg_high ),&
            soln%Umms( neq, grid%ig_low:grid%ig_high,  grid%jg_low:grid%jg_high ),&
            soln%Smms( neq,grid%i_cell_low:grid%i_cell_high, grid%j_cell_low:grid%j_cell_high ),&  
            )


            soln%U     = zero
            soln%Fxi   = zero
            soln%Feta  = zero
            soln%S     = zero
            soln%V     = zero
            soln%R     = one
            soln%asnd  = one
            soln%rinit = one
            ! soln%mach  = one
            ! soln%temp  = one
            soln%dt    = one
            soln%rnorm = one
         
            
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
                    soln%mach,  &
                    soln%temp,  &
                    soln%dt,    &
                    soln%rnorm, &
                    soln%rinit,&
                    ! soln%psi_p_xi,&
                    ! soln%psi_m_xi, &
                    ! soln%psi_p_eta,&
                    ! soln%psi_m_eta 
                    soln%Umms ,&
                    soln%Vmms, &
                    soln%Smms ,&
                    ) 
                end subroutine deallocate_soln    
    
end module soln_type