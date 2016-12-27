!--------------------------------------------------------------------------------------------
! Routine to calculate skin moment od inertia.
subroutine SkinMomentOfInertia(yc,Ixx)
    USE PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: Skin,SkinGeometry1,SkinGeometry2
    use IntegrationVariables, only: Mg
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints
    

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::     Ixx					! inertia moment
    !Local Variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y,z				! point coordinates
    real(8) ::                  zeta,eta			! point coordinates
    real(8) ::                  x1,x2				! limits of integration
    real(8) ::                  z1,z2				! limits of integration
    integer ::                  n,m_g				! iteration variables
    real(8) ::                  zetaskF,zetaskA 
    real(8) ::                  slope 
    real(8),allocatable ::      gi(:),zetai(:)   

    Ixx	= 0.0D0
    if(Skin(Ntr)%Nsk.EQ.0) return
    allocate(gi(Mg),zetai(Mg))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Skin(Ntr)%Nsk
            call SkinGeometry1(n,zetai(m_g),zetaskF,zetaskA)	! compute zetaskF, zetaskA !!!!!!!!!!!!!Changed :etaj(Ng) to zetai(m_g)!!!!!!
            zeta = 0.5D0*(1.0D0-zetai(m_g))*zetaskF+0.5D0*(1.0D0+zetai(m_g))*zetaskA
            call AxesTranformation(1,zeta,eta,x,y)		! compute x, y
            call AxesTranformation(1,zetaskF,eta,x1,y)	! compute x1, y
            call AxesTranformation(1,zetaskA,eta,x2,y)	! compute x2, y
            call SkinGeometry2(3,n,x,y,z,z1,z2,slope)	! compute z1, z2
            G	= (x2-x1)/2.0D0*(z2*z2*z2-z1*z1*z1)/3.0D0
            Ixx	= Ixx+Skin(Ntr)%Eysk(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine SkinMomentOfInertia

!--------------------------------------------------------------------------------------------
! Routine to calculate spar web moment of inertia.
subroutine SparWebMomentOfInertia(yc,Ixx)

    use SparVariablesRoutines, only: Spar,SparGeometry,SparWebGeometry
    use IntegrationVariables, only: Mg
    use PlanformVariables, only: Ntr
    
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::     Ixx					! inertia moment
    !Local Variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y					! point coordinates
    real(8) ::                  eta,zeta			! point coordinates
    real(8) ::                  x1,x2				! limits of integration
    real(8) ::                  z1,z2				! limits of integration
    integer ::                  n,m_g				! iteration variables
    real(8),allocatable ::      gi(:),zetai(:)    

    Ixx	= 0.0D0
    if(Spar(Ntr)%Nsp.EQ.0) return
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    allocate(gi(Mg),zetai(Mg))
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Spar(Ntr)%Nsp
            !call SparGeometry(n,eta)							! compute Spar(Ntr)%zetasp
            call SparWebGeometry(1,n,zetai(m_g),eta,zeta,x,y,z1,z2)			! compute Spar(Ntr)%tcsw, zeta
            call AxesTranformation(1,zeta,eta,x,y)					! compute x, y
            x1	= x-(Spar(Ntr)%tsw1(n)-Spar(Ntr)%tsw2(n)*eta)/2.0D0
            x2	= x+(Spar(Ntr)%tsw1(n)-Spar(Ntr)%tsw2(n)*eta)/2.0D0
            call SparWebGeometry(2,n,zetai(m_g),eta,zeta,x,y,z1,z2)			! compute z1, z2
            G	= (x2-x1)/2.0D0*(z2*z2*z2-z1*z1*z1)/3.0D0
            Ixx	= Ixx+Spar(Ntr)%Eesc(n)*gi(m_g)*G						! corrigir Spar(Ntr)%Eesc(n)
        end do
    end do
    deallocate(gi,zetai)

end subroutine SparWebMomentOfInertia

!--------------------------------------------------------------------------------------------
! Routine to calculate spar cap moment of inertia.
subroutine SparCapMomentOfInertia(yc,Ixx)

    use SparVariablesRoutines, only:Spar,SparGeometry,SparCapGeometry
    use IntegrationVariables, only: Mg
    use PlanformVariables, only: Ntr
    
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::     Ixx					! inertia moment
    !Local Variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y					! point coordinates
    real(8) ::                  eta,zeta			! point coordinates
    real(8) ::                  x1,x2				! limits of integration
    real(8) ::                  z1U,z2U,z1L,z2L		! limits of integration
    integer ::                  n,m_g				! iteration variables
    real(8),allocatable ::      gi(:),zetai(:)    

    Ixx	= 0.0D0
    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),zetai(Mg))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Spar(Ntr)%Nsp
            !call SparGeometry(n,eta)											! compute Spar(Ntr)%zetasp
            call SparCapGeometry(1,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute Spar(Ntr)%lcsc, zeta
            call AxesTranformation(1,zeta,eta,x,y)			! compute x, y
            x1	= x-(Spar(Ntr)%lsc1(n)-Spar(Ntr)%lsc2(n)*eta)/2.0D0
            x2	= x+(Spar(Ntr)%lsc1(n)-Spar(Ntr)%lsc2(n)*eta)/2.0D0
            call SparCapGeometry(2,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
            G	= (x2-x1)/2.0D0*((z2U*z2U*z2U-z1U*z1U*z1U)+(z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0
            Ixx	= Ixx+Spar(Ntr)%Eesc(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine SparCapMomentOfInertia

!--------------------------------------------------------------------------------------------
! Routine to calculate spar moment of inertia.
subroutine SparMomentOfInertia(yc,Ixx)

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::     Ixx					! inertia moment
    !Local Variables
    real(8) ::                  Ixxsw				! spar web moment of inertia
    real(8) ::                  Ixxsc				! spar cap moment of inertia

    call SparWebMomentOfInertia(yc,Ixxsw)
    call SparCapMomentOfInertia(yc,Ixxsc)
    Ixx	= Ixxsw+Ixxsc

end subroutine SparMomentOfInertia

!--------------------------------------------------------------------------------------------
! Routine to calculate rib web moment of inertia.
subroutine RibWebMomentOfInertia(yc,Ixx)

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::     Ixx					! inertia moment
    !Local Variables

    !k		= k
    Ixx	= yc*0.0D0

end subroutine RibWebMomentOfInertia

!--------------------------------------------------------------------------------------------
! Routine to calculate rib cap moment of inertia.
subroutine RibCapMomentOfInertia(yc,Ixx)

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::     Ixx					! inertia moment
    !Local Variables

    !k	= k
    Ixx	= yc*0.0D0

end subroutine RibCapMomentOfInertia

!--------------------------------------------------------------------------------------------
! Routine to calculate rib moment of inertia.
subroutine RibMomentOfInertia(yc,Ixx)

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::     Ixx					! inertia moment
    !Local Variables
    real(8) ::                  Ixxrw				! rib web moment of inertia
    real(8) ::                  Ixxrc				! rib cap moment of inertia

    call RibWebMomentOfInertia(yc,Ixxrw)
    call RibCapMomentOfInertia(yc,Ixxrc)
    Ixx	= Ixxrw+Ixxrc

end subroutine RibMomentOfInertia

!--------------------------------------------------------------------------------------------
! Routine to calculate stringer moment of inertia.
subroutine StringerMomentOfInertia(yc,Ixx)
    use PlanformVariables, only: Ntr
    use StringerVariablesRoutines, only: Stringer,StringerGeometry1,StringerGeometry2
    use IntegrationVariables, only: Mg
    
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::     Ixx					! inertia moment
    !Local Variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y					! point coordinates
    real(8) ::                  eta,zeta			! point coordinates
    real(8) ::                  x1,x2				! limits of integration
    real(8) ::                  z1U,z2U,z1L,z2L		! limits of integration
    integer ::                  n,m_g				! iteration variables
    real(8),allocatable ::      gi(:),zetai(:)    

    Ixx	= 0.0
    if(Stringer(Ntr)%Nst.EQ.0) return
    allocate(gi(Mg),zetai(Mg))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0+2.0*yc
    do m_g=1,Mg
        do n=1,Stringer(Ntr)%Nst
            !call StringerGeometry1(n,eta)						! compute Stringer(Ntr)%zetast
            call StringerGeometry2(1,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute Stringer(Ntr)%lcst, zeta
            call AxesTranformation(1,zeta,eta,x,y)					! compute x, y
            x1	= x-(Stringer(Ntr)%lst1(n)-Stringer(Ntr)%lst2(n)*eta)/2.0D0
            x2	= x+(Stringer(Ntr)%lst1(n)-Stringer(Ntr)%lst2(n)*eta)/2.0D0
            call StringerGeometry2(2,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
            G	= (x2-x1)/2.0D0*((z2U*z2U*z2U-z1U*z1U*z1U)+(z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0
            Ixx	= Ixx+Stringer(Ntr)%Eest(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine StringerMomentOfInertia

!--------------------------------------------------------------------------------------------
! Routine to calculate overall moment of inertia.
subroutine MomentOfInertia(yc,E0,Ixx)
    USE PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: Skin

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    real(8), intent(in) ::      E0					! reference elastic modulus
    !Output variables
    real(8), intent(out) ::     Ixx					! inertia moment
    !Local Variables
    real(8) ::                  Ixxsk				! skin moment of inertia
    real(8) ::                  Ixxsp				! spar moment of inertia
    real(8) ::                  Ixxrb				! rib moment of inertia
    real(8) ::                  Ixxst				! stringer moment of inertia
    

    write(*,*) 'Computing moment of inertia...'

    !E0	= Skin(Ntr)%Eysk(1)

    write(*,*) '    - skin'
    call SkinMomentOfInertia(yc,Ixxsk)
    !  write(*,'(1X,A8,F12.9,A3)') 'Ixxsk = ',Ixxsk/E0,'m^4'

    write(*,*) '    - spar'
    call SparMomentOfInertia(yc,Ixxsp)
    !  write(*,'(1X,A8,F12.9,A3)') 'Ixxsp = ',Ixxsp/E0,'m^4'

    write(*,*) '    - rib'
    call RibMomentOfInertia(yc,Ixxrb)
    !  write(*,'(1X,A8,F12.9,A3)') 'Ixxrb = ',Ixxrb/E0,'m^4'

    write(*,*) '    - stringer'
    call StringerMomentOfInertia(yc,Ixxst)
    !  write(*,'(1X,A8,F12.9,A3)') 'Ixxst = ',Ixxst/E0,'m^4'

    Ixx	= (Ixxsk+Ixxsp+Ixxrb+Ixxst)/E0 
    write(*,'(1X,A7,F12.9,A3)') 'Ixx  = ',Ixx,' [m^4]'

    !  Pz	= 0.0
    !  do i=1,CLoad(Ntr)%Ncl
    !    Pz	= Pz+CLoad(Ntr)%Pclz(i)
    !  end do
    !  write(*,'(1X,A6,F12.9)')	  'k   = ',Pz*(b2(Ntr)/cos(Lamdapl(Ntr)))**2/(E0*Ixx)	   !/cos(Lamdapl(Ntr))
    !  w	= Pz*(b2(Ntr)/cos(Lamdapl(Ntr)))**3/(3.0*E0*Ixx)
    !  write(*,'(1X,A4,F12.6,A1,A11,F12.6)') 'w = ',w,'m','     w/L = ',w/(b2(Ntr)/cos(Lamdapl(Ntr)))

end subroutine MomentOfInertia

!--------------------------------------------------------------------------------------------
! Routine to calculate skin mass.
subroutine SkinMass(M)
    USE PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: Skin,SkinGeometry1,SkinGeometry2
    use IntegrationVariables, only: Mg,Ng
    
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Output variables
    real(8), intent(out) ::     M					! mass
    !Local variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y,z				! point coordinates
    real(8) ::                  zeta				! zeta-coordinate
    real(8) ::                  z1,z2				! limits of integration
    integer ::                  n,m_g,n_g			! iteration variables	  
    real(8) ::                  zetaskF,zetaskA
    real(8) ::                  DJ,J1(2,2)          ! jacobinan determinant and inverse Jacobian matrix  
    real(8) ::                  slope
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    
   
    M	= 0.0D0
    if(Skin(Ntr)%Nsk.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Skin(Ntr)%Nsk
                call SkinGeometry1(n,etaj(Ng),zetaskF,zetaskA)		! compute zetaskF, zetaskA
                zeta	= 0.5D0*(1.0D0-zetai(m_g))*zetaskF+0.5D0*(1.0D0+zetai(m_g))*zetaskA
                call AxesTranformation(1,zeta,etaj(n_g),x,y)					! compute x, y
                call AxesTranformation(2,zeta,etaj(n_g),x,y,DJ,J1)				! compute DJ
                call SkinGeometry2(3,n,x,y,z,z1,z2,slope)						! compute z1, z2
                G = 0.5D0*(zetaskA-zetaskF)*(z2-z1)*DJ
                M = M+Skin(Ntr)%rhosk(n)*gi(m_g)*gj(n_g)*G
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)

end subroutine SkinMass

!--------------------------------------------------------------------------------------------
! Routine to calculate spar mass.
subroutine SparWebMass(M)

    use SparVariablesRoutines, only:Spar,SparGeometry,SparWebGeometry
    use IntegrationVariables, only: Mg,Ng
    use PlanformVariables, only: Ntr
    
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Output variables
    real(8), intent(out) ::     M					! mass
    !Local variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y					! point coordinates
    real(8) ::                  zeta				! point coordinates
    real(8) ::                  z1,z2				! limits of integration
    real(8) ::                  DJ,J1(2,2)          ! jacobinan determinant and inverse Jacobian matrix  
    integer ::                  n,m_g,n_g			! iteration variables
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    

    M	= 0.0D0
    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Spar(Ntr)%Nsp
                !call SparGeometry(n,etaj(n_g))					! compute Spar(Ntr)%zetasp
                call SparWebGeometry(1,n,zetai(m_g),etaj(n_g),zeta,x,y,z1,z2)	! compute Spar(Ntr)%tcsw, zeta
                call AxesTranformation(1,zeta,etaj(n_g),x,y,DJ,J1)			! compute x, y
                call AxesTranformation(3,zeta,etaj(n_g),x,y,DJ,J1)			! compute DJ, J1
                call SparWebGeometry(2,n,zetai(m_g),etaj(n_g),zeta,x,y,z1,z2)	! compute z1, z2
                G = Spar(Ntr)%tcsw*(z2-z1)*DJ
                M = M+Spar(Ntr)%rhosw(n)*gi(m_g)*gj(n_g)*G
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)

end subroutine SparWebMass

!--------------------------------------------------------------------------------------------
! Routine to calculate spar cap mass.
subroutine SparCapMass(M)

    use SparVariablesRoutines, only: Spar,SparGeometry,SparCapGeometry
    use IntegrationVariables, only: Mg,Ng
    use PlanformVariables, only: Ntr
    
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Output variables
    real(8), intent(out) ::     M					! mass
    !Local variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y					! point coordinates
    real(8) ::                  zeta                                ! point coordinates
    real(8) ::                  z1U,z2U,z1L,z2L                     ! limits of integration
    real(8) ::                  DJ,J1(2,2)                          ! jacobinan determinant and inverse Jacobian matrix  
    integer ::                  n,m_g,n_g                           ! iteration variables
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    

    M	= 0.0D0
    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Spar(Ntr)%Nsp
                !call SparGeometry(n,etaj(n_g))						! compute Spar(Ntr)%zetasp
                call SparCapGeometry(1,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute Spar(Ntr)%lcsc, zeta
                call AxesTranformation(1,zeta,etaj(n_g),x,y,DJ,J1)				! compute x, y
                call AxesTranformation(3,zeta,etaj(n_g),x,y,DJ,J1)				! compute DJ, J1
                call SparCapGeometry(2,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
                G = Spar(Ntr)%lcsc*((z2U-z1U)+(z2L-z1L))*DJ
                M = M+Spar(Ntr)%rhosc(n)*gi(m_g)*gj(n_g)*G
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)

end subroutine SparCapMass

!--------------------------------------------------------------------------------------------
! Routine to calculate spar mass.
subroutine SparMass(M)

    implicit none

    !Output variables
    real(8), intent(out) ::     M					! mass
    !Local variables
    real(8) ::                  Masssw				! spar web mass
    real(8) ::                  Masssc				! spar cap mass

    call SparWebMass(Masssw)
    call SparCapMass(Masssc)
    M	= Masssw+Masssc

end subroutine SparMass

!--------------------------------------------------------------------------------------------
! Routine to calculate rib web mass.
subroutine RibWebMass(M)
    use PlanformVariables, only: Ntr
    use RibVariablesRoutines, only:Rib,RibGeometry,RibWebGeometry
    use IntegrationVariables, only: Mg,Ng
    
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Output variables
    real(8), intent(out) ::     M					! mass
    !Local variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y					! point coordinates
    real(8) ::                  eta					! eta-coordinate
    real(8) ::                  z1,z2                   ! limits of integration
    real(8) ::                  DJ,J1(2,2)              ! jacobinan determinant and inverse Jacobian matrix  
    integer ::                  n,m_g,n_g               ! iteration variables
    logical ::                  error
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    

    M	= 0.0D0
    if(Rib(Ntr)%Nrb.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Rib(Ntr)%Nrb
                !call RibGeometry(n,zetai(m_g))									! compute Rib(Ntr)%etarb
                call RibWebGeometry(1,n,zetai(m_g),etaj(n_g),eta,x,y,z1,z2,error)		! compute Rib(Ntr)%tcrw, eta
                call AxesTranformation(1,zetai(m_g),eta,x,y,DJ,J1)					! compute x, y
                call AxesTranformation(3,zetai(m_g),eta,x,y,DJ,J1)					! compute DJ, J1
                call RibWebGeometry(2,n,zetai(m_g),etaj(n_g),eta,x,y,z1,z2,error)		! compute z1, z2
                G = Rib(Ntr)%tcrw*(z2-z1)*DJ
                M = M+Rib(Ntr)%rhorw(n)*gi(m_g)*gj(n_g)*G
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)

end subroutine RibWebMass

!--------------------------------------------------------------------------------------------
! Routine to calculate rib cap mass.
subroutine RibCapMass(M)
    use PlanformVariables, only: Ntr
    use RibVariablesRoutines, only:Rib,RibGeometry,RibCapGeometry
    use IntegrationVariables, only: Mg,Ng
    
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Output variables
    real(8), intent(out) ::     M					! mass
    !Local variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y					! point coordinates
    real(8) ::                  eta					! eta-coordinate
    real(8) ::                  z1U,z2U,z1L,z2L		        ! limits of integration
    real(8) ::                  DJ,J1(2,2)                          ! jacobinan determinant and inverse Jacobian matrix  
    integer ::                  n,m_g,n_g			        ! iteration variables
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    

    M	= 0.0D0
    if(Rib(Ntr)%Nrb.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Rib(Ntr)%Nrb
                !call RibGeometry(n,zetai(m_g))											! compute Rib(Ntr)%etarb
                call RibCapGeometry(1,n,zetai(m_g),etaj(n_g),eta,x,y,z1U,z2U,z1L,z2L)	! compute Rib(Ntr)%lcrc, eta
                call AxesTranformation(1,zetai(m_g),eta,x,y,DJ,J1)							! compute x, y
                call AxesTranformation(3,zetai(m_g),eta,x,y,DJ,J1)							! compute DJ, J1
                call RibCapGeometry(2,n,zetai(m_g),etaj(n_g),eta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
                G = Rib(Ntr)%lcrc*((z2U-z1U)+(z2L-z1L))*DJ
                M = M+Rib(Ntr)%rhorc(n)*gi(m_g)*gj(n_g)*G
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)

end subroutine RibCapMass

!--------------------------------------------------------------------------------------------
! Routine to calculate rib mass.
subroutine RibMass(M)

    implicit none

    !Output variables
    real(8), intent(out) ::     M					! mass
    !Local variables
    real(8) ::                  Massrw				! rib web mass
    real(8) ::                  Massrc				! rib cap mass

    call RibWebMass(Massrw)
    call RibCapMass(Massrc)
    M	= Massrw+Massrc

end subroutine RibMass

!--------------------------------------------------------------------------------------------
! Routine to calculate stringer mass.
subroutine StringerMass(M)
    use PlanformVariables, only: Ntr
    use StringerVariablesRoutines, only:Stringer, StringerGeometry1,StringerGeometry2
    use IntegrationVariables, only: Mg,Ng
    
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Output variables
    real(8), intent(out) ::     M					! mass
    !Local variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y					! point coordinates
    real(8) ::                  zeta				! point coordinates
    real(8) ::                  z1U,z2U,z1L,z2L		! limits of integration
    real(8) ::                  DJ,J1(2,2)          ! jacobinan determinant and inverse Jacobian matrix  
    integer ::                  n,m_g,n_g			! iteration variables
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    

    M	= 0.0D0
    if(Stringer(Ntr)%Nst.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Stringer(Ntr)%Nst
                !call StringerGeometry1(n,etaj(n_g))						! compute Stringer(Ntr)%zetast
                call StringerGeometry2(1,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute Stringer(Ntr)%lcst, zeta
                call AxesTranformation(1,zeta,etaj(n_g),x,y,DJ,J1)				! compute x, y
                call AxesTranformation(3,zeta,etaj(n_g),x,y,DJ,J1)				! compute DJ, J1
                call StringerGeometry2(2,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
                G = Stringer(Ntr)%lcst*((z2U-z1U)+(z2L-z1L))*DJ
                M = M+Stringer(Ntr)%rhost(n)*gi(m_g)*gj(n_g)*G
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)

end subroutine StringerMass

!--------------------------------------------------------------------------------------------
! Routine to calculate overall mass.
subroutine Mass(M)

    implicit none

    !Output variables
    real(8), intent(out) ::		M				! wing mass
    !Local variables
    real(8) ::				Masssk				! skin mass
    real(8) ::				Masssp				! spar mass
    real(8) ::				Massrb				! rib mass
    real(8) ::				Massst				! stringer mass


    write(*,*) 'Computing mass...'

    write(*,*) '    - skin'
    call SkinMass(Masssk)
    !  write(*,'(1X,A9,F12.3,A2)') 'Masssk = ',Masssk,'kg'

    write(*,*) '    - spar'
    call SparMass(Masssp)
    !  write(*,'(1X,A9,F12.3,A2)') 'Masssp = ',Masssp,'kg'

    write(*,*) '    - rib'
    call RibMass(Massrb)
    !  write(*,'(1X,A9,F12.3,A2)') 'Massrb = ',Massrb,'kg'

    write(*,*) '    - stringer'
    call StringerMass(Massst)
    !  write(*,'(1X,A9,F12.3,A2)') 'Massst = ',Massst,'kg'

    M	= Masssk+Masssp+Massrb+Massst 
    write(*,'(1X,A7,F12.3,A2)') 'Mass = ',M,'kg'

end subroutine Mass

