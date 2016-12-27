!--------------------------------------------------------------------------------------------
! Routine to calculate skin mass matrix.
subroutine SkinMassMatrix(k,M)
    
    use PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: Skin,SkinGeometry1,SkinGeometry2
    use IntegrationVariables, only: Mg,Ng 
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints,TrialPolynomials,MultiplyVectorByTransposedVector
    use GeneralVariables, only: PrintMatrix
  
    implicit none

    !Input variables
    integer, intent(in) ::	    k					! number of polynomials
    !Output variables  
    real(8), intent(out) :: 	M(5*k*k,5*k*k)		! mass matrix
    !Local variables
    integer ::				    i,j,n,m_g,n_g		! iteration variables
    real(8), allocatable  ::    BB(:)				! Legendre polynimial multiplication vector
    real(8), allocatable  ::    HH(:,:)			    ! Diagonal matrix of multiplication of BB by BB
    real(8) ::				    x,y,z				! point coordinates
    real(8) ::				    zeta				! zeta-coordinate
    real(8) ::				    z1,z2				! limits of integration
    real(8) ::				    I1,I2,I3
    real(8) ::			        zetaskF,zetaskA 
    real(8) ::                  slope               ! Skin dz/dy
    real(8) ::        		    DJ,J1(2,2)          ! jacobinan determinant and inverse Jacobian matrix
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)  
    character(60) ::			name				! matrix name

    M = 0.0D0

    if(Skin(Ntr)%Nsk.EQ.0) return
    allocate(BB(k*k),HH(k*k,k*k))
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))

    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Skin(Ntr)%Nsk
                call SkinGeometry1(n,etaj(n_g),zetaskF,zetaskA)									! compute zetaskF, zetaskA
                zeta = 0.5D0*(1.0D0-zetai(m_g))*zetaskF+0.5D0*(1.0D0+zetai(m_g))*zetaskA
                call AxesTranformation(1,zeta,etaj(n_g),x,y)					! compute x, y
                call AxesTranformation(2,zeta,etaj(n_g),x,y,DJ,J1)              ! compute DJ
                call SkinGeometry2(3,n,x,y,z,z1,z2,slope)						! compute z1, z2
                call TrialPolynomials(1,k,zeta,etaj(n_g),BB)
                call MultiplyVectorByTransposedVector(BB,HH)
                !Integration constants
                I1=0.5D0*(zetaskA-zetaskF)*(z2-z1)*DJ
                I2=0.5D0*(zetaskA-zetaskF)*(z2*z2-z1*z1)/2.0D0*DJ
                I3=0.5D0*(zetaskA-zetaskF)*(z2*z2*z2-z1*z1*z1)/3.0D0*DJ
                forall(i=1:k*k,j=1:k*k)
                    ! M11. 
                    M(i,j)	= M(i,j)+Skin(Ntr)%rhosk(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M14.
                    M(i,j+3*k*k)	= M(i,j+3*k*k)+Skin(Ntr)%rhosk(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                    ! M22.
                    M(i+k*k,j+k*k)	= M(i+k*k,j+k*k)+Skin(Ntr)%rhosk(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M25.
                    M(i+k*k,j+4*k*k)	= M(i+k*k,j+4*k*k)+Skin(Ntr)%rhosk(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                    ! M33.
                    M(i+2*k*k,j+2*k*k)	= M(i+2*k*k,j+2*k*k)+Skin(Ntr)%rhosk(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M41.
                    M(i+4*k*k,j)		= M(j,i+4*k*k)
                    ! M44.
                    M(i+3*k*k,j+3*k*k)	= M(i+3*k*k,j+3*k*k)+Skin(Ntr)%rhosk(n)*gi(m_g)*gj(n_g)*HH(i,j)*I3
                    ! M52.
                    M(i+4*k*k,j+k*k)	= M(j+k*k,i+4*k*k)
                    ! M55.
                    M(i+4*k*k,j+4*k*k)	= M(i+4*k*k,j+4*k*k)+Skin(Ntr)%rhosk(n)*gi(m_g)*gj(n_g)*HH(i,j)*I3
                end forall
            end do
        end do
    end do

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'skin mass matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if
    deallocate(BB,HH)
    deallocate(gi,gj,zetai,etaj)

end subroutine SkinMassMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate spar web mass matrix.
subroutine SparWebMassMatrix(k,M)
    use PlanformVariables, only: Ntr
    use SparVariablesRoutines, only: Spar,SparGeometry,SparWebGeometry
    use IntegrationVariables, only: Mg,Ng 
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints,TrialPolynomials,MultiplyVectorByTransposedVector
    use GeneralVariables, only: PrintMatrix
 
    implicit none

    !Input variables
    integer, intent(in) ::	    k					! number of polynomials
    !Outputt variables
    real(8), intent(out) ::	    M(5*k*k,5*k*k)		! mass matrix
    !Local variables
    integer ::				    i,j,n,m_g,n_g		! iteration variables
    real(8), allocatable  ::    BB(:)				! Legendre polynimial multiplication vector
    real(8), allocatable  ::    HH(:,:)			    ! Diagonal matrix of multiplication of BB by BB
    real(8) ::				    x,y					! point coordinates
    real(8) ::				    zeta				! zeta-coordinate
    real(8) ::				    z1,z2				! limits of integration
    real(8) ::				    I1,I2,I3
    real(8) ::        		    DJ,J1(2,2)          ! jacobinan determinant and inverse Jacobian matrix  
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)  
    character(60) ::			name				! matrix name

    M = 0.0D0

    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(BB(k*k),HH(k*k,k*k))
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Spar(Ntr)%Nsp
                call SparWebGeometry(1,n,zetai(m_g),etaj(n_g),zeta,x,y,z1,z2)	! compute Spar(Ntr)%tcsw, zeta
                call AxesTranformation(1,zeta,etaj(n_g),x,y,DJ,J1)					! compute x, y
                call AxesTranformation(2,zeta,etaj(n_g),x,y,DJ,J1)					! compute DJ
                call SparWebGeometry(2,n,zetai(m_g),etaj(n_g),zeta,x,y,z1,z2)	! compute z1, z2
                call TrialPolynomials(1,k,zeta,etaj(n_g),BB)
                call MultiplyVectorByTransposedVector(BB,HH)
                ! integration constants		
                I1=Spar(Ntr)%tcsw*(z2-z1)*DJ
                I2=Spar(Ntr)%tcsw*(z2*z2-z1*z1)/2.0D0*DJ
                I3=Spar(Ntr)%tcsw*(z2*z2*z2-z1*z1*z1)/3.0D0*DJ
                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j) = M(i,j)+Spar(Ntr)%rhosw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M14.
                    M(i,j+3*k*k) = M(i,j+3*k*k)+Spar(Ntr)%rhosw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                    ! M22.
                    M(i+k*k,j+k*k) = M(i+k*k,j+k*k)+Spar(Ntr)%rhosw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M25.
                    M(i+k*k,j+4*k*k) = M(i+k*k,j+4*k*k)+Spar(Ntr)%rhosw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                    ! M33.
                    M(i+2*k*k,j+2*k*k) = M(i+2*k*k,j+2*k*k)+Spar(Ntr)%rhosw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M41.
                    M(i+3*k*k,j) = M(j,i+3*k*k)
                    ! M44.
                    M(i+3*k*k,j+3*k*k) = M(i+3*k*k,j+3*k*k)+Spar(Ntr)%rhosw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I3
                    ! M52.
                    M(i+4*k*k,j+k*k) = M(j+k*k,i+4*k*k)
                    ! M55.
                    M(i+4*k*k,j+4*k*k) = M(i+4*k*k,j+4*k*k)+Spar(Ntr)%rhosw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I3
                end forall
            end do
        end do
    end do

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'spar web mass matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if
    deallocate(BB,HH)
    deallocate(gi,gj,zetai,etaj)

end subroutine SparWebMassMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate spar cap mass matrix.
subroutine SparCapMassMatrix(k,M)
    use PlanformVariables, only: Ntr
    use SparVariablesRoutines, only: Spar,SparGeometry,SparCapGeometry
    use IntegrationVariables, only: Mg,Ng 
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints,TrialPolynomials,MultiplyVectorByTransposedVector
    use GeneralVariables, only: PrintMatrix

    implicit none

    !Input variables
    integer, intent(in) ::            k					! number of polynomials
    !Output variables
    real(8), intent(out) ::           M(5*k*k,5*k*k)		! mass matrix
    !Local variables
    integer ::                        i,j,n,m_g,n_g		! iteration variables
    real(8), allocatable  ::          BB(:)				! Legendre polynimial multiplication vector
    real(8), allocatable  ::          HH(:,:)			    ! Diagonal matrix of multiplication of BB by BB
    real(8) ::                        x,y					! point coordinates
    real(8) ::                        zeta				! zeta-coordinate
    real(8) ::                        z1U,z2U,z1L,z2L		! limits of integration
    real(8) ::                        I1,I2,I3
    real(8) ::                        DJ,J1(2,2)          ! jacobinan determinant and inverse Jacobian matrix  
    real(8),allocatable ::            gi(:),gj(:),zetai(:),etaj(:)  
    character(60) ::                  name				! matrix name

    M = 0.0D0

    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(BB(k*k),HH(k*k,k*k))
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Spar(Ntr)%Nsp
                call SparCapGeometry(1,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute Spar(Ntr)%lcsc, zeta
                call AxesTranformation(1,zeta,etaj(n_g),x,y,DJ,J1)							! compute x, y
                call AxesTranformation(2,zeta,etaj(n_g),x,y,DJ,J1)							! compute DJ
                call SparCapGeometry(2,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
                call TrialPolynomials(1,k,zeta,etaj(n_g),BB)
                call MultiplyVectorByTransposedVector(BB,HH)
                !Integration constants
                I1=Spar(Ntr)%lcsc*((z2U-z1U)+(z2L-z1L))*DJ
                I2=Spar(Ntr)%lcsc*((z2U*z2U-z1U*z1U)+(z2L*z2L-z1L*z1L))/2.0D0*DJ
                I3=Spar(Ntr)%lcsc*((z2U*z2U*z2U-z1U*z1U*z1U) + (z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0*DJ
                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j) = M(i,j)+Spar(Ntr)%rhosc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M14.
                    M(i,j+3*k*k) = M(i,j+3*k*k)+Spar(Ntr)%rhosc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                    ! M22.
                    M(i+k*k,j+k*k) = M(i+k*k,j+k*k)+Spar(Ntr)%rhosc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M25.
                    M(i+k*k,j+4*k*k) = M(i+k*k,j+4*k*k)+Spar(Ntr)%rhosc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                    ! M33.
                    M(i+2*k*k,j+2*k*k) = M(i+2*k*k,j+2*k*k)+Spar(Ntr)%rhosc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M41.
                    M(i+3*k*k,j) = M(j,i+3*k*k)
                    ! M44.
                    M(i+3*k*k,j+3*k*k) = M(i+3*k*k,j+3*k*k)+Spar(Ntr)%rhosc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I3
                    ! M52.
                    M(i+4*k*k,j+k*k) = M(j+k*k,i+4*k*k)
                    ! M55.
                    M(i+4*k*k,j+4*k*k)	= M(i+4*k*k,j+4*k*k)+Spar(Ntr)%rhosc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                end forall
            end do
        end do
    end do

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'spar cap mass matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if
    deallocate(BB,HH)
    deallocate(gi,gj,zetai,etaj)

end subroutine SparCapMassMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate spar mass matrix.
subroutine SparMassMatrix(k,M)

  use GeneralVariables, only: PrintMatrix

  implicit none

  !Input variables
  integer, intent(in) ::	k					! number of polynomials
  !Output variables
  real(8), intent(out) ::	M(5*k*k,5*k*k)		! mass matrix
  !Local variables
  real(8), allocatable  ::  Msw(:,:)	        ! spar web mass matrix
  real(8), allocatable  ::  Msc(:,:)	        ! spar cap mass matrix
  character(60) ::			name				! matrix name

  allocate(Msw(5*k*k,5*k*k),Msc(5*k*k,5*k*k))
  call SparWebMassMatrix(k,Msw)
  call SparCapMassMatrix(k,Msc)
  
  M = Msw + Msc 

! Write matrix into a file.
  if(PrintMatrix) then
    name	= 'spar mass matrix'
    call WriteMatrix(M,5*k*k,5*k*k,name)
  end if
  deallocate(Msw,Msc)

end subroutine SparMassMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate rib web mass matrix.
subroutine RibWebMassMatrix(k,M)
    use PlanformVariables, only: Ntr
    use RibVariablesRoutines, only: Rib,RibGeometry,RibWebGeometry
    use IntegrationVariables, only: Mg,Ng 
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints,TrialPolynomials,MultiplyVectorByTransposedVector
    use GeneralVariables, only: PrintMatrix
  

    implicit none
    !Input variables
    integer, intent(in) ::	    k					! number of polynomials
    !Output variables
    real(8), intent(out) ::	    M(5*k*k,5*k*k)		! mass matrix
    !Local Variables
    integer ::				    i,j,n,m_g,n_g		! iteration variables
    real(8), allocatable  ::    BB(:)				! Legendre polynimial multiplication vector
    real(8), allocatable  ::    HH(:,:)			    ! Diagonal matrix of multiplication of BB by BB
    real(8) ::				    x,y					! point coordinates
    real(8) ::				    eta					! eta-coordinate
    real(8) ::				    z1,z2				! limits of integration
    real(8) ::				    I1,I2,I3
    real(8) ::        		    DJ,J1(2,2)          ! jacobinan determinant and inverse Jacobian matrix
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)
    logical                     error    
    character(60) ::			name				! matrix name

    M	= 0.0D0

    if(Rib(Ntr)%Nrb.EQ.0) return
    allocate(BB(k*k),HH(k*k,k*k))
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Rib(Ntr)%Nrb
                call RibWebGeometry(1,n,zetai(m_g),etaj(n_g),eta,x,y,z1,z2,error)		! compute Rib(Ntr)%tcrw, eta
                call AxesTranformation(1,zetai(m_g),eta,x,y)					! compute x, y
                call AxesTranformation(2,zetai(m_g),eta,x,y,DJ,J1)					! compute DJ
                call RibWebGeometry(2,n,zetai(m_g),etaj(n_g),eta,x,y,z1,z2,error)		! compute z1, z2
                call TrialPolynomials(1,k,zetai(m_g),eta,BB)
                call MultiplyVectorByTransposedVector(BB,HH)
                !Integration constants
                I1=Rib(Ntr)%tcrw*(z2-z1)*DJ
                I2=Rib(Ntr)%tcrw*(z2*z2-z1*z1)/2.0D0*DJ
                I3=Rib(Ntr)%tcrw*(z2*z2*z2-z1*z1*z1)/3.0D0*DJ
                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j) = M(i,j)+Rib(Ntr)%rhorw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M14.
                    M(i,j+3*k*k) = M(i,j+3*k*k)+Rib(Ntr)%rhorw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                    ! M22.
                    M(i+k*k,j+k*k) = M(i+k*k,j+k*k)+Rib(Ntr)%rhorw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M25.
                    M(i+k*k,j+4*k*k) = M(i+k*k,j+4*k*k)+Rib(Ntr)%rhorw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                    ! M33.
                    M(i+2*k*k,j+2*k*k) = M(i+2*k*k,j+2*k*k)+Rib(Ntr)%rhorw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M41.
                    M(i+3*k*k,j) = M(j,i+3*k*k)
                    ! M44.
                    M(i+3*k*k,j+3*k*k) = M(i+3*k*k,j+3*k*k)+Rib(Ntr)%rhorw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I3
                    ! M52.
                    M(i+4*k*k,j+k*k) = M(j+k*k,i+4*k*k)
                    ! M55.
                    M(i+4*k*k,j+4*k*k) = M(i+4*k*k,j+4*k*k)+Rib(Ntr)%rhorw(n)*gi(m_g)*gj(n_g)*HH(i,j)*I3
                end forall
            end do
        end do
    end do

! Write matrix into a file.
  if(PrintMatrix) then
    name	= 'rib web mass matrix'
    call WriteMatrix(M,5*k*k,5*k*k,name)
  end if
  deallocate(BB,HH)
  deallocate(gi,gj,zetai,etaj)
end subroutine RibWebMassMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate rib cap mass matrix.
subroutine RibCapMassMatrix(k,M)
    use PlanformVariables, only: Ntr
    use RibVariablesRoutines, only: Rib,RibGeometry,RibCapGeometry
    use IntegrationVariables, only: Mg,Ng
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints,TrialPolynomials,MultiplyVectorByTransposedVector
    use GeneralVariables, only: PrintMatrix
  
    implicit none
    !Input variables
    integer, intent(in) ::	    k					! number of polynomials
    !Output variables
    real(8), intent(out) ::	    M(5*k*k,5*k*k)		! mass matrix
    !Local Variables
    integer ::				    i,j,n,m_g,n_g		! iteration variables
    real(8), allocatable  ::    BB(:)				! Legendre polynimial multiplication vector
    real(8), allocatable  ::    HH(:,:)			    ! Diagonal matrix of multiplication of BB by BB
    real(8) ::				    x,y					! point coordinates
    real(8) ::				    eta					! eta-coordinate
    real(8) ::				    z1U,z2U,z1L,z2L		! limits of integration
    real(8) ::				    I1,I2,I3
    real(8) ::        		    DJ,J1(2,2)          ! jacobinan determinant and inverse Jacobian matrix
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    
    character(60) ::			name				! matrix name

    M = 0.0D0

    if(Rib(Ntr)%Nrb.EQ.0) return
    allocate(BB(k*k),HH(k*k,k*k))
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))

    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Rib(Ntr)%Nrb
                call RibCapGeometry(1,n,zetai(m_g),etaj(n_g),eta,x,y,z1U,z2U,z1L,z2L)	! compute Rib(Ntr)%lcrc, eta
                call AxesTranformation(1,zetai(m_g),eta,x,y)							! compute x, y
                call AxesTranformation(2,zetai(m_g),eta,x,y,DJ,J1)                      ! compute DJ
                call RibCapGeometry(2,n,zetai(m_g),etaj(n_g),eta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
                call TrialPolynomials(1,k,zetai(m_g),eta,BB)
                call MultiplyVectorByTransposedVector(BB,HH)
                ! integration constants
                I1=Rib(Ntr)%lcrc*((z2U-z1U)+(z2L-z1L))*DJ
                I2=Rib(Ntr)%lcrc*((z2U*z2U-z1U*z1U)+(z2L*z2L-z1L*z1L))/2.0D0*DJ
                I3=Rib(Ntr)%lcrc*((z2U*z2U*z2U-z1U*z1U*z1U)+(z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0*DJ
                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j) = M(i,j)+Rib(Ntr)%rhorc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M14.
                    M(i,j+3*k*k) = M(i,j+3*k*k)+Rib(Ntr)%rhorc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                    ! M22.
                    M(i+k*k,j+k*k) = M(i+k*k,j+k*k)+Rib(Ntr)%rhorc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M25.
                    M(i+k*k,j+4*k*k) = M(i+k*k,j+4*k*k)+Rib(Ntr)%rhorc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I2
                    ! M33.
                    M(i+2*k*k,j+2*k*k) = M(i+2*k*k,j+2*k*k)+Rib(Ntr)%rhorc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I1
                    ! M41.
                    M(i+3*k*k,j) = M(j,i+3*k*k)
                    ! M44.
                    M(i+3*k*k,j+3*k*k) = M(i+3*k*k,j+3*k*k)+Rib(Ntr)%rhorc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I3
                    ! M52.
                    M(i+4*k*k,j+k*k) = M(j+k*k,i+4*k*k)
                    ! M55.
                    M(i+4*k*k,j+4*k*k) = M(i+4*k*k,j+4*k*k)+Rib(Ntr)%rhorc(n)*gi(m_g)*gj(n_g)*HH(i,j)*I3
                end forall
            end do
        end do
    end do

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'rib cap mass matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if
    deallocate(BB,HH)
    deallocate(gi,gj,zetai,etaj)

end subroutine RibCapMassMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate rib mass matrix.
subroutine RibMassMatrix(k,M)

    use GeneralVariables, only: PrintMatrix

    implicit none
    !Input variables
    integer, intent(in) ::	    k					! number of polynomials
    !Output variables
    real(8), intent(out) ::	    M(5*k*k,5*k*k)		! mass matrix
    !Local variables
    real(8), allocatable  ::    Mrw(:,:)	        ! rib web mass matrix
    real(8), allocatable  ::    Mrc(:,:)	        ! rib cap mass matrix
    character(60) ::			name                ! matrix name

    allocate(Mrw(5*k*k,5*k*k),Mrc(5*k*k,5*k*k))

    call RibWebMassMatrix(k,Mrw)
    call RibCapMassMatrix(k,Mrc)

    M = Mrw + Mrc 

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'rib mass matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if
    deallocate(Mrw,Mrc)
end subroutine RibMassMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate stringer mass matrix.
subroutine StringerMassMatrix(k,M)
    use PlanformVariables, only: Ntr
    use StringerVariablesRoutines, only: Stringer,StringerGeometry1,StringerGeometry2
    use IntegrationVariables, only: Mg,Ng 
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints,TrialPolynomials,MultiplyVectorByTransposedVector
    use GeneralVariables, only: PrintMatrix
  

    implicit none
    !Input variables
    integer, intent(in) ::	    k					! number of polynomials
    !Output variables
    real(8), intent(out) ::	    M(5*k*k,5*k*k)		! mass matrix
    !Local variables
    integer ::				    i,j,n,m_g,n_g		! iteration variables
    real(8), allocatable  ::    BB(:)				! Legendre polynimial multiplication vector
    real(8), allocatable  ::    HH(:,:)			    ! Diagonal matrix of multiplication of BB by BB
    real(8) ::				    x,y					! point coordinates
    real(8) ::				    zeta				! zeta-coordinate
    real(8) ::				    z1U,z2U,z1L,z2L		! limits of integration
    real(8) ::        		    DJ,J1(2,2)          ! jacobinan determinant and inverse Jacobian matrix 
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)   
    character(60) ::			name				! matrix name

    M = 0.0D0

    if(Stringer(Ntr)%Nst.EQ.0) return
    allocate(BB(k*k),HH(k*k,k*k))
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Stringer(Ntr)%Nst
                call StringerGeometry2(1,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute Stringer(Ntr)%lcst, zeta
                call AxesTranformation(1,zeta,etaj(n_g),x,y,DJ,J1)							! compute x, y
                call AxesTranformation(2,zeta,etaj(n_g),x,y,DJ,J1)							! compute DJ
                call StringerGeometry2(2,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
                call TrialPolynomials(1,k,zeta,etaj(n_g),BB)
                call MultiplyVectorByTransposedVector(BB,HH)
                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j)	= M(i,j)+Stringer(Ntr)%rhost(n)*gi(m_g)*gj(n_g)*Stringer(Ntr)%lcst*HH(i,j)*((z2U-z1U)+(z2L-z1L))*DJ
                    ! M14.
                    M(i,j+3*k*k) = M(i,j+3*k*k)+Stringer(Ntr)%rhost(n)*gi(m_g)*gj(n_g)*Stringer(Ntr)%lcst*HH(i,j)*((z2U*z2U-z1U*z1U)+(z2L*z2L-z1L*z1L))/2.0D0*DJ
                    ! M22.
                    M(i+k*k,j+k*k) = M(i+k*k,j+k*k)+Stringer(Ntr)%rhost(n)*gi(m_g)*gj(n_g)*Stringer(Ntr)%lcst*HH(i,j)*((z2U-z1U)+(z2L-z1L))*DJ
                    ! M25.
                    M(i+k*k,j+4*k*k) = M(i+k*k,j+4*k*k)+Stringer(Ntr)%rhost(n)*gi(m_g)*gj(n_g)*Stringer(Ntr)%lcst*HH(i,j)*((z2U*z2U-z1U*z1U)+(z2L*z2L-z1L*z1L))/2.0D0*DJ
                    ! M33.
                    M(i+2*k*k,j+2*k*k) = M(i+2*k*k,j+2*k*k)+Stringer(Ntr)%rhost(n)*gi(m_g)*gj(n_g)*Stringer(Ntr)%lcst*HH(i,j)*((z2U-z1U)+(z2L-z1L))*DJ
                    ! M44.
                    M(i+3*k*k,j+3*k*k) = M(i+3*k*k,j+3*k*k)+Stringer(Ntr)%rhost(n)*gi(m_g)*gj(n_g)*Stringer(Ntr)%lcst*HH(i,j)*((z2U*z2U*z2U-z1U*z1U*z1U)+(z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0*DJ
                    ! M55.
                    M(i+4*k*k,j+4*k*k) = M(i+4*k*k,j+4*k*k)+Stringer(Ntr)%rhost(n)*gi(m_g)*gj(n_g)*Stringer(Ntr)%lcst*HH(i,j)*((z2U*z2U*z2U-z1U*z1U*z1U)+(z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0*DJ
                    ! M41.
                    M(i+3*k*k,j) = M(j,i+3*k*k)
                    ! M52.
                    M(i+4*k*k,j+k*k) = M(j+k*k,i+4*k*k)
                end forall
            end do
        end do
    end do

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'stringer mass matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if
    deallocate(BB,HH)
    deallocate(gi,gj,zetai,etaj)
end subroutine StringerMassMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate planform mass matrix.
subroutine PlanformMassMatrix(k)
  use PlanformVariables, only: Ntr
  use DeformationVariables, only: Deformation
  use GeneralVariables, only: PrintMatrix

  implicit none
  !Input variables
  integer, intent(in) ::    k
  !Local variables
  real(8), allocatable  ::  Msk(:,:)	    ! skin mass matrix
  real(8), allocatable  ::  Msp(:,:)	    ! spar mass matrix
  real(8), allocatable  ::  Mrb(:,:)	    ! rib mass matrix
  real(8), allocatable  ::  Mst(:,:)	    ! stringer mass matrix
  character(60) ::			name			! matrix name

  allocate(Msk(5*k*k,5*k*k),Msp(5*k*k,5*k*k),Mrb(5*k*k,5*k*k),Mst(5*k*k,5*k*k))

  write(*,'(1X,A,1X,I1,1X,A)') 'Computing trapezoid', Ntr, 'mass matrix...'

!$OMP PARALLEL SECTIONS
!$OMP SECTION
  write(*,*) '    - skin'
  call SkinMassMatrix(k,Msk)
!$OMP SECTION
  write(*,*) '    - spar'
  call SparMassMatrix(k,Msp)
!$OMP SECTION
  write(*,*) '    - rib'
  call RibMassMatrix(k,Mrb)
!$OMP SECTION
  write(*,*) '    - stringer'
  call StringerMassMatrix(k,Mst)
!$OMP END PARALLEL SECTIONS

  Deformation(Ntr)%MM = Msk + Msp + Mrb + Mst 

! Write matrix into a file.
  if(PrintMatrix) then
    name	= 'planform mass matrix'
    call WriteMatrix(Deformation(Ntr)%MM,5*k*k,5*k*k,name)
  end if
  deallocate(Msk,Msp,Mrb,Mst)

end subroutine PlanformMassMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate overall mass matrix.
subroutine OverallMassMatrix(k)
    use PlanformVariables, only: Npl
    use DeformationVariables, only: Deformation,MM_Global
    use GeneralVariables, only: PrintMatrix

    implicit none
    !Input variables
    integer, intent(in) ::  k
    !Local variables
    character(60) ::        name			! matrix name
    integer ::              i,j             ! loop var

    write(*,*) 'Computing overall mass matrix...'
    
    MM_Global = 0.0D0
    forall(i=1:Npl,j=1:5*k*k)
        MM_Global((i-1)*5*k*k + 1:i*5*k*k,(i-1)*5*k*k + j ) = Deformation(i)%MM(:,j)
    end forall  
    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'overall mass matrix'
        call WriteMatrix(MM_Global,size(MM_Global,DIM=1),size(MM_Global,DIM=2),name)
    end if

end subroutine OverallMassMatrix