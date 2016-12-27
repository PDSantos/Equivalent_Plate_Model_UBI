
!--------------------------------------------------------------------------------------------
! Routine to calculate skin section area.
subroutine SkinSectionArea(yc,A)
    use PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: Skin,SkinGeometry1,SkinGeometry2
    use IntegrationVariables, only: Mg,Ng 
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none
  
    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::     A					! section area
    !Local Variables
    real(8) ::                  G					! function to integrate
    real(8) ::                  x,y,z				! point coordinates
    real(8) ::                  zeta,eta			! point coordinates
    real(8) ::                  x1,x2				! limits of integration
    real(8) ::                  z1,z2				! limits of integration
    real(8) ::                  zetaskF,zetaskA
    real(8) ::                  slope                       ! Skin dz/dy
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    
    integer ::                  n,m_g				! iteration variables	  

    A	= 0.0D0
    if(Skin(Ntr)%Nsk.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Skin(Ntr)%Nsk
            call SkinGeometry1(n,etaj(Ng),zetaskF,zetaskA)                                  ! compute zetaskF, zetaskA
            zeta = 0.5*(1.0-zetai(m_g))*zetaskF+0.5*(1.0+zetai(m_g))*zetaskA
            call AxesTranformation(1,zetai(m_g),eta,x,y)                                    ! compute x, y
            call AxesTranformation(1,zetaskF,eta,x1,y)                                      ! compute x1, y
            call AxesTranformation(1,zetaskA,eta,x2,y)                                      ! compute x2, y
            call SkinGeometry2(3,n,x,y,z,z1,z2,slope)                                       ! compute z1, z2
            G	= (x2-x1)/2.0*(z2-z1)
            A	= A+Skin(Ntr)%Eysk(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,gj,zetai,etaj)

end subroutine SkinSectionArea

!--------------------------------------------------------------------------------------------
! Routine to calculate spar web section area.
subroutine SparWebSectionArea(yc,A)
    use PlanformVariables, only: Ntr
    use SparVariablesRoutines
    use IntegrationVariables, only: Mg
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Input variables
    real(8), intent(in) ::          yc				! fraction of semi-span
    !Output variables
    real(8), intent(out) ::         A					! section area
    !Local Variables
    real(8) ::				G				! function to integrate
    real(8) ::				x,y				! point coordinates
    real(8) ::				eta,zeta			! point coordinates
    real(8) ::				x1,x2				! limits of integration
    real(8) ::				z1,z2				! limits of integration
    integer ::				n,m_g				! iteration variables
    real(8),allocatable ::          gi(:),zetai(:)    

    A	= 0.0D0
    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),zetai(Mg))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Spar(Ntr)%Nsp
            !call SparGeometryn,eta)											! compute Spar(Ntr)%zetasp
            call SparWebGeometry(1,n,zetai(m_g),eta,zeta,x,y,z1,z2)			! compute Spar(Ntr)%tcsw, zeta
            call AxesTranformation(1,zeta,eta,x,y)							! compute x, y
            x1	= x-(Spar(Ntr)%tsw1(n)-Spar(Ntr)%tsw2(n)*eta)/2.0
            x2	= x+(Spar(Ntr)%tsw1(n)-Spar(Ntr)%tsw2(n)*eta)/2.0
            call SparWebGeometry(2,n,zetai(m_g),eta,zeta,x,y,z1,z2)			! compute z1, z2
            G	= (x2-x1)/2.0*(z2-z1)
            A	= A+Spar(Ntr)%Eesw(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)
end subroutine SparWebSectionArea

!--------------------------------------------------------------------------------------------
! Routine to calculate spar cap section area.
subroutine SparCapSectionArea(yc,A)
    use PlanformVariables, only: Ntr
    use SparVariablesRoutines
    use IntegrationVariables, only: Mg
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::     A					! section area
    !Local Variables
    real(8) ::				G					! function to integrate
    real(8) ::				x,y					! point coordinates
    real(8) ::				eta,zeta			! point coordinates
    real(8) ::				x1,x2				! limits of integration
    real(8) ::				z1U,z2U,z1L,z2L		! limits of integration
    integer ::				n,m_g				! iteration variables
    real(8),allocatable ::          gi(:),zetai(:)    

    A	= 0.0D0
    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),zetai(Mg))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Spar(Ntr)%Nsp
            !call SparGeometryn,eta)											! compute Spar(Ntr)%zetasp
            call SparCapGeometry(1,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute Spar(Ntr)%lcsc, zeta
            call AxesTranformation(1,zeta,eta,x,y)							! compute x, y
            x1	= x-(Spar(Ntr)%lsc1(n)-Spar(Ntr)%lsc2(n)*eta)/2.0D0
            x2	= x+(Spar(Ntr)%lsc1(n)-Spar(Ntr)%lsc2(n)*eta)/2.0D0
            call SparCapGeometry(2,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
            G	= (x2-x1)/2.0D0*((z2U-z1U)+(z2L-z1L))
            A	= A+Spar(Ntr)%Eesc(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine SparCapSectionArea

!--------------------------------------------------------------------------------------------
! Routine to calculate spar section area.
subroutine SparSectionArea(yc,A)

    implicit none

    !Input variables
    real(8), intent(in) ::    yc					! fraction of semi-span!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Output variables
    real(8), intent(out) ::   A					! section area
    !Local Variables
    real(8) ::			Asw					! spar web section area
    real(8) ::			Asc					! spar cap section area

    call SparWebSectionArea(yc,Asw)
    call SparCapSectionArea(yc,Asc)
    A = Asw+Asc

end subroutine SparSectionArea

!--------------------------------------------------------------------------------------------
! Routine to calculate rib web section area.
subroutine RibWebSectionArea(yc,A)

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Output variables
    real(8), intent(out) ::     A					! section area
    !Local Variables

    !k	= k
    A	= yc*0.0D0

end subroutine RibWebSectionArea

!--------------------------------------------------------------------------------------------
! Routine to calculate rib cap section area.
subroutine RibCapSectionArea(yc,A)

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Output variables
    real(8), intent(out) ::     A					! section area
    !Local Variables

   !k	= k
    A	= yc*0.0D0

end subroutine RibCapSectionArea

!--------------------------------------------------------------------------------------------
! Routine to calculate rib section area.
subroutine RibSectionArea(yc,A)

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Output variables
    real(8), intent(out) ::     A					! section area
    !Local Variables
    real(8) ::				Arw					! rib web section area
    real(8) ::				Arc					! rib cap section area

    !  call RibWebMomentOfInertia(k,Arw)
    !  call RibCapMomentOfInertia(k,Arc)
    call RibWebSectionArea(yc,Arw)
    call RibCapSectionArea(yc,Arc)
    A	= Arw+Arc

end subroutine RibSectionArea

!--------------------------------------------------------------------------------------------
! Routine to calculate stringer section area.
subroutine StringerSectionArea(yc,A)
    use PlanformVariables, only: Ntr
    use StringerVariablesRoutines
    use IntegrationVariables, only: Mg
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints
    
    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Output variables
    real(8), intent(out) ::             A					! section area
    !Local Variables
    real(8) ::                          G					! function to integrate
    real(8) ::				x,y					! point coordinates
    real(8) ::				eta,zeta			! point coordinates
    real(8) ::				x1,x2				! limits of integration
    real(8) ::				z1U,z2U,z1L,z2L		! limits of integration
    integer ::				n,m_g				! iteration variables
    real(8),allocatable ::          gi(:),zetai(:)    

    A	= 0.0D0
    if(Stringer(Ntr)%Nst.EQ.0) return
    allocate(gi(Mg),zetai(Mg))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Stringer(Ntr)%Nst
            !call StringerGeometry1(n,eta)											! compute Stringer(Ntr)%zetast
            call StringerGeometry2(1,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute Stringer(Ntr)%lcst, zeta
            call AxesTranformation(1,zeta,eta,x,y)								! compute x, y
            x1	= x-(Stringer(Ntr)%lst1(n)-Stringer(Ntr)%lst2(n)*eta)/2.0D0
            x2	= x+(Stringer(Ntr)%lst1(n)-Stringer(Ntr)%lst2(n)*eta)/2.0D0
            call StringerGeometry2(2,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
            G	= (x2-x1)/2.0D0*((z2U-z1U)+(z2L-z1L))
            A	= A+Stringer(Ntr)%Eest(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine StringerSectionArea

!--------------------------------------------------------------------------------------------
! Routine to calculate overall section area.
subroutine SectionArea(screen,yc,E0,Area)

    implicit none

    !Input variables
    logical, intent(in) ::      screen				! TRUE: write info to scree; False otherwise
    real(8), intent(in) ::      yc					! fraction of semi-span
    real(8), intent(in) ::      E0					! reference elastic modulus
    !Output variables
    real(8), intent(out) ::     Area                ! section area
    !Local Variables
    real(8) ::    				Ask					! skin section area
    real(8) ::    				Asp					! spar section area
    real(8) ::    				Arb					! rib section area
    real(8) ::    				Ast					! stringer section area

    if(screen) write(*,*) 'Computing section area...'

    if(screen) write(*,*) '    - skin'
    call SkinSectionArea(yc,Ask)
    !  write(*,'(1X,A8,F12.9,A3)') 'Ask   = ',Ask/E0,'m^2'

    if(screen) write(*,*) '    - spar'
    call SparSectionArea(yc,Asp)
    !  write(*,'(1X,A8,F12.9,A3)') 'Asp   = ',Asp/E0,'m^2'

    if(screen) write(*,*) '    - rib'
    call RibSectionArea(yc,Arb)
    !  write(*,'(1X,A8,F12.9,A3)') 'Arb   = ',Arb/E0,'m^2'

    if(screen) write(*,*) '    - stringer'
    call StringerSectionArea(yc,Ast)
    !  write(*,'(1X,A8,F12.9,A3)') 'Ast   = ',Ast/E0,'m^2'

    Area = (Ask+Asp+Arb+Ast)/E0 
    if(screen) write(*,'(1X,A7,F12.9,A3)') 'Area = ',Area,'m^2'

end subroutine SectionArea

!--------------------------------------------------------------------------------------------
! Routine to calculate skin section first moment of area.
subroutine SkinSectionFirstMomentOfArea(yc,Am)
    use PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: Skin,SkinGeometry1,SkinGeometry2
    use IntegrationVariables, only: Mg
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints
    
    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Output variables
    real(8), intent(out) ::     Am					! section area
    !Local Variables
    real(8) ::                    G					! function to integrate
    real(8) ::                    x,y,z				! point coordinates
    real(8) ::                    zeta,eta			! point coordinates
    real(8) ::                    x1,x2				! limits of integration
    real(8) ::                    z1,z2				! limits of integration
    real(8) ::                    zetaskF,zetaskA
    real(8) ::                    slope             ! Skin dz/dy
    integer ::                    n,m_g				! iteration variables	  
    real(8),allocatable ::        gi(:),zetai(:) 


    Am	= 0.0D0
    if(Skin(Ntr)%Nsk.EQ.0) return
    allocate(gi(Mg),zetai(Mg)) 
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Skin(Ntr)%Nsk
            call SkinGeometry1(n,zetai(m_g),zetaskF,zetaskA)			! compute zetaskF, zetaskA ---BEFORE etaj -ERROR???????
            zeta = 0.5D0*(1.0D0-zetai(m_g))*zetaskF+0.5D0*(1.0D0+zetai(m_g))*zetaskA
            call AxesTranformation(1,zetai(m_g),eta,x,y)			! compute x, y
            call AxesTranformation(1,zetaskF,eta,x1,y)			! compute x1, y
            call AxesTranformation(1,zetaskA,eta,x2,y)			! compute x2, y
            call SkinGeometry2(3,n,x,y,z,z1,z2,slope)			! compute z1, z2
            !	  write(*,*) z,z1,z2
            G	= (x2-x1)/2.0D0*(z2*z2-z1*z1)/2.0D0
            Am	= Am+Skin(Ntr)%Eysk(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine SkinSectionFirstMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate spar web section first moment of area.
subroutine SparWebSectionFirstMomentOfArea(yc,Am)

    use PlanformVariables, only: Ntr
    use SparVariablesRoutines
    use IntegrationVariables, only: Mg
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Output variables
    real(8), intent(out) ::     Am					! section area
    !Local Variables
    real(8) ::				G					! function to integrate
    real(8) ::				x,y					! point coordinates
    real(8) ::				eta,zeta			! point coordinates
    real(8) ::				x1,x2				! limits of integration
    real(8) ::				z1,z2				! limits of integration
    integer ::				n,m_g				! iteration variables
    real(8),allocatable ::            gi(:),zetai(:) 

    Am	= 0.0D0
    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),zetai(Mg)) 
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Spar(Ntr)%Nsp
            !call SparGeometryn,eta)						! compute Spar(Ntr)%zetasp
            call SparWebGeometry(1,n,zetai(m_g),eta,zeta,x,y,z1,z2)		! compute Spar(Ntr)%tcsw, zeta
            call AxesTranformation(1,zeta,eta,x,y)				! compute x, y
            x1	= x-(Spar(Ntr)%tsw1(n)-Spar(Ntr)%tsw2(n)*eta)/2.0D0
            x2	= x+(Spar(Ntr)%tsw1(n)-Spar(Ntr)%tsw2(n)*eta)/2.0D0
            call SparWebGeometry(2,n,zetai(m_g),eta,zeta,x,y,z1,z2)		! compute z1, z2
            G	= (x2-x1)/2.0*(z2*z2-z1*z1)/2.0D0
            Am	= Am+Spar(Ntr)%Eesw(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine SparWebSectionFirstMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate spar cap section first moment of area.
subroutine SparCapSectionFirstMomentOfArea(yc,Am)
    
    use PlanformVariables, only: Ntr
    use SparVariablesRoutines
    use IntegrationVariables, only: Mg
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints
    
    implicit none

    !Input variables
    real(8), intent(in) ::      yc					! fraction of semi-span!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Output variables
    real(8), intent(out) ::     Am					! section area
    !Local Variables
    real(8) ::				G					! function to integrate
    real(8) ::				x,y					! point coordinates
    real(8) ::				eta,zeta			! point coordinates
    real(8) ::				x1,x2				! limits of integration
    real(8) ::				z1U,z2U,z1L,z2L		! limits of integration
    integer ::				n,m_g				! iteration variables
    real(8),allocatable ::            gi(:),zetai(:) 

    Am	= 0.0D0
    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),zetai(Mg)) 
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Spar(Ntr)%Nsp
            !call SparGeometryn,eta)							! compute Spar(Ntr)%zetasp
            call SparCapGeometry(1,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute Spar(Ntr)%lcsc, zeta
            call AxesTranformation(1,zeta,eta,x,y)					! compute x, y
            x1	= x-(Spar(Ntr)%lsc1(n)-Spar(Ntr)%lsc2(n)*eta)/2.0D0
            x2	= x+(Spar(Ntr)%lsc1(n)-Spar(Ntr)%lsc2(n)*eta)/2.0D0
            call SparCapGeometry(2,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
            !write(*,'(I5,4F12.6)') m_g,z1U,z2U,z1L,z2L
            G	= (x2-x1)/2.0D0*((z2U*z2U-z1U*z1U)+(z2L*z2L-z1L*z1L))/2.0D0
            Am	= Am+Spar(Ntr)%Eesc(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine SparCapSectionFirstMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate spar section first moment of area.
subroutine SparSectionFirstMomentOfArea(yc,Am)

    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Am					! section area
    !Local Variables
    real(8) ::				Amsw				! spar web section first moment of area
    real(8) ::				Amsc				! spar cap section first moment of area

    call SparWebSectionFirstMomentOfArea(yc,Amsw)
    call SparCapSectionFirstMomentOfArea(yc,Amsc)
    Am	= Amsw+Amsc

end subroutine SparSectionFirstMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate rib web section first moment of area.
subroutine RibWebSectionFirstMomentOfArea(yc,Am)

    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Am					! section area
    !Local Variables

    !k	= k
    Am	= yc*0.0D0

end subroutine RibWebSectionFirstMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate rib cap section first moment of area.
subroutine RibCapSectionFirstMomentOfArea(yc,Am)

    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Am					! section area
    !Local Variables

    !k	= k
    Am	= yc*0.0D0

end subroutine RibCapSectionFirstMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate rib section first moment of area.
subroutine RibSectionFirstMomentOfArea(yc,Am)

    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Am					! section area
    !Local Variables
    real(8) ::				Amrw				! rib web section first moment of area
    real(8) ::				Amrc				! rib cap section first moment of area

    !  call RibWebMomentOfInertia(k,Amrw)
    !  call RibCapMomentOfInertia(k,Amrc)
    call RibWebSectionFirstMomentOfArea(yc,Amrw)
    call RibCapSectionFirstMomentOfArea(yc,Amrc)
    Am	= Amrw+Amrc

end subroutine RibSectionFirstMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate stringer section first moment of area.
subroutine StringerSectionFirstMomentOfArea(yc,Am)
    use PlanformVariables, only: Ntr
    use StringerVariablesRoutines
    use IntegrationVariables, only: Mg
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Am					! section area
    !Local Variables
    real(8) ::				G					! function to integrate
    real(8) ::				x,y					! point coordinates
    real(8) ::				eta,zeta			! point coordinates
    real(8) ::				x1,x2				! limits of integration
    real(8) ::				z1U,z2U,z1L,z2L		! limits of integration
    integer ::				n,m_g				! iteration variables
    real(8),allocatable ::            gi(:),zetai(:) 

    Am	= 0.0D0
    if(Stringer(Ntr)%Nst.EQ.0) return
    allocate(gi(Mg),zetai(Mg)) 
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Stringer(Ntr)%Nst
            !call StringerGeometry1(n,eta)			! compute Stringer(Ntr)%zetast
            call StringerGeometry2(1,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute Stringer(Ntr)%lcst, zeta
            call AxesTranformation(1,zeta,eta,x,y)	! compute x, y
            x1	= x-(Stringer(Ntr)%lst1(n)-Stringer(Ntr)%lst2(n)*eta)/2.0D0
            x2	= x+(Stringer(Ntr)%lst1(n)-Stringer(Ntr)%lst2(n)*eta)/2.0D0
            call StringerGeometry2(2,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
            !  write(*,'(I5,4F12.6)') m_g,z1U,z2U,z1L,z2L
            G = (x2-x1)/2.0D0*((z2U*z2U-z1U*z1U)+(z2L*z2L-z1L*z1L))/2.0D0
            Am	= Am+Stringer(Ntr)%Eest(n)*gi(m_g)*G
        end do
    end do
    deallocate(zetai,gi)

end subroutine StringerSectionFirstMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate overall section area.
subroutine SectionFirstMomentOfAreaZ(screen,yc,E0,Amz)

    use PlanformVariables, only:b2

    implicit none

    !Input variables
    logical, intent(in) ::              screen				! TRUE: write info to scree; False otherwise
    real(8), intent(in) ::              yc					! fraction of semi-span
    real(8), intent(in) ::              E0					! reference elastic modulus
    !Output variables
    real(8), intent(out) ::             Amz					! section area
    !Local Variables
    real(8) ::				Amsk				! skin section first moment of area
    real(8) ::				Amsp				! spar section first moment of area
    real(8) ::				Amrb				! rib section first moment of area
    real(8) ::				Amst				! stringer section first moment of area

    if(screen) write(*,*) 'Computing section first moment of area at y = ',yc*b2(1),' [m]'    
  !!$OMP PARALLEL SECTIONS
  !!$OMP SECTION
    if(screen) write(*,*) '    - skin'
    call SkinSectionFirstMomentOfArea(yc,Amsk)
    !  write(*,'(1X,A8,F12.9,A3)') 'Amsk  = ',Amsk/E0,'m^3'
  !!$OMP SECTION
    if(screen) write(*,*) '    - spar'
    call SparSectionFirstMomentOfArea(yc,Amsp)
    !  write(*,'(1X,A8,F12.9,A3)') 'Amsp  = ',Amsp/E0,'m^3'
  !!$OMP SECTION
    if(screen) write(*,*) '    - rib'
    call RibSectionFirstMomentOfArea(yc,Amrb)
    !  write(*,'(1X,A8,F12.9,A3)') 'Amrb  = ',Amrb/E0,'m^3'
  !!$OMP SECTION
    if(screen) write(*,*) '    - stringer'
    call StringerSectionFirstMomentOfArea(yc,Amst)
    !  write(*,'(1X,A8,F12.9,A3)') 'Amst  = ',Amst/E0,'m^3'
  !!$OMP END PARALLEL SECTIONS
    Amz	= (Amsk+Amsp+Amrb+Amst)/E0 
    if(screen) write(*,'(1X,A7,F12.9,A)') 'Am   = ',Amz,' [m^3]'

end subroutine SectionFirstMomentOfAreaZ

!--------------------------------------------------------------------------------------------
! Routine to calculate overall section centroid.
subroutine SectionCentroid(screen,yc,zcentroid)

    use PlanformVariables, only: Ntr,b2,yroot
    use SkinVariablesRoutines, only: Skin

    implicit none

    !Input variables
    logical, intent(in) ::              screen				! TRUE: write info to scree; False otherwise
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             zcentroid           ! section z centroid position
    !Local Variables
    real(8) ::                          Area                ! section Area
    real(8) ::                          Amz					! section First Moment of Area
    real(8) ::                          E0					! reference elastic modulus
    
    if(screen) write(*,*) 'Computing section centroid at y = ',yroot(Ntr)+yc*b2(Ntr),' [m]'    
    E0	= Skin(Ntr)%Eysk(1)
    !Compute Section area
    call SectionArea(screen,yc,E0,Area)
    !Compute Section First Moment of Inertia 
    call SectionFirstMomentOfAreaZ(screen,yc,E0,Amz)
    !Compute section centroid
    zcentroid	= Amz/Area 
    if(screen) write(*,'(1X,A7,F12.9,A3)') 'z    = ',zcentroid,'[m]  '

end subroutine SectionCentroid

!--------------------------------------------------------------------------------------------
! Routine to calculate skin moment od inertia.
subroutine SkinSectionSecondMomentOfArea(yc,Ixx)
    use PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: Skin,SkinGeometry1,SkinGeometry2
    use IntegrationVariables
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Ixx					! section area
    !Local Variables
    real(8) ::                          G					! function to integrate
    real(8) ::                          x,y,z				! point coordinates
    real(8) ::                          zeta,eta            ! point coordinates
    real(8) ::                          x1,x2				! limits of integration
    real(8) ::                          z1,z2				! limits of integration
    real(8) ::                          zetaskF,zetaskA
    real(8) ::                          slope               ! Skin dz/dy
    integer ::                          n,m_g				! iteration variables	  
    real(8),allocatable ::              gi(:),zetai(:)   

    Ixx	= 0.0D0
    if(Skin(Ntr)%Nsk.EQ.0) return
    allocate(gi(Mg),zetai(Mg))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Skin(Ntr)%Nsk
            call SkinGeometry1(n,zetai(m_g),zetaskF,zetaskA)			! compute zetaskF, zetaskA    !!!!!!!!!!Changed: etaj(Ng) to zetai(m_g)
            zeta = 0.5D0*(1.0D0-zetai(m_g))*zetaskF+0.5D0*(1.0D0+zetai(m_g))*zetaskA
            call AxesTranformation(1,zetai(m_g),eta,x,y)						! compute x, y
            call AxesTranformation(1,zetaskF,eta,x1,y)						! compute x1, y
            call AxesTranformation(1,zetaskA,eta,x2,y)						! compute x2, y
            call SkinGeometry2(3,n,x,y,z,z1,z2,slope)							! compute z1, z2
            G	= (x2-x1)/2.0D0*(z2*z2*z2-z1*z1*z1)/3.0D0
            Ixx	= Ixx+Skin(Ntr)%Eysk(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine SkinSectionSecondMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate spar web section second moment of area.
subroutine SparWebSectionSecondMomentOfArea(yc,Ixx)

    use PlanformVariables, only: Ntr
    use SparVariablesRoutines
    use IntegrationVariables
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Input variables
    real(8), intent(in) ::              yc				! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Ixx				! section area
    !Local Variables
    real(8) ::                          G				! function to integrate
    real(8) ::                          x,y				! point coordinates
    real(8) ::                          eta,zeta			! point coordinates
    real(8) ::                          x1,x2				! limits of integration
    real(8) ::                          z1,z2				! limits of integration
    integer ::                          n,m_g				! iteration variables
    real(8),allocatable ::              gi(:),zetai(:)   

    Ixx	= 0.0D0
    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),zetai(Mg))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Spar(Ntr)%Nsp
            !call SparGeometryn,eta)							! compute Spar(Ntr)%zetasp
            call SparWebGeometry(1,n,zetai(m_g),eta,zeta,x,y,z1,z2)			! compute Spar(Ntr)%tcsw, zeta
            call AxesTranformation(1,zeta,eta,x,y)					! compute x, y
            x1	= x-(Spar(Ntr)%tsw1(n)-Spar(Ntr)%tsw2(n)*eta)/2.0D0
            x2	= x+(Spar(Ntr)%tsw1(n)-Spar(Ntr)%tsw2(n)*eta)/2.0D0
            call SparWebGeometry(2,n,zetai(m_g),eta,zeta,x,y,z1,z2)			! compute z1, z2
            G = (x2-x1)/2.0D0*(z2*z2*z2-z1*z1*z1)/3.0D0
            Ixx	= Ixx+Spar(Ntr)%Eesw(n)*gi(m_g)*G						! corrigir Spar(Ntr)%Eesc(n)
        end do
    end do
    deallocate(gi,zetai)

end subroutine SparWebSectionSecondMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate spar cap section second moment of area.
subroutine SparCapSectionSecondMomentOfArea(yc,Ixx)

    use PlanformVariables, only: Ntr
    use SparVariablesRoutines
    use IntegrationVariables
    
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints
    
    implicit none

    !Input variables
    real(8), intent(in) ::              yc				! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Ixx				! section area
    !Local Variables
    real(8) ::                          G				! function to integrate
    real(8) ::                          x,y				! point coordinates
    real(8) ::                          eta,zeta			! point coordinates
    real(8) ::                          x1,x2				! limits of integration
    real(8) ::                          z1U,z2U,z1L,z2L		! limits of integration
    integer ::                          n,m_g				! iteration variables
    real(8),allocatable ::              gi(:),zetai(:)   

    Ixx	= 0.0D0
    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),zetai(Mg))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta = -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Spar(Ntr)%Nsp
            !call SparGeometry(n,eta)							! compute Spar(Ntr)%zetasp
            call SparCapGeometry(1,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute Spar(Ntr)%lcsc, zeta
            call AxesTranformation(1,zeta,eta,x,y)					! compute x, y
            x1	= x-(Spar(Ntr)%lsc1(n)-Spar(Ntr)%lsc2(n)*eta)/2.0D0
            x2	= x+(Spar(Ntr)%lsc1(n)-Spar(Ntr)%lsc2(n)*eta)/2.0D0
            call SparCapGeometry(2,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
            G	= (x2-x1)/2.0D0*((z2U*z2U*z2U-z1U*z1U*z1U)+(z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0
            Ixx	= Ixx+Spar(Ntr)%Eesc(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine SparCapSectionSecondMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate spar section second moment of area.
subroutine SparSectionSecondMomentOfArea(yc,Ixx)

    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Ixx					! section second moment of inertia
    !Local Variables
    real(8) ::				Ixxsw				! spar web section second moment of area
    real(8) ::				Ixxsc				! spar cap section second moment of area

    call SparWebSectionSecondMomentOfArea(yc,Ixxsw)
    call SparCapSectionSecondMomentOfArea(yc,Ixxsc)
    Ixx	= Ixxsw+Ixxsc

end subroutine SparSectionSecondMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate rib web section second moment of area.
subroutine RibWebSectionSecondMomentOfArea(yc,Ixx)

    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Ixx					! section second moment of inertia
    !Local Variables

    !k = k
    Ixx	= yc*0.0D0

end subroutine RibWebSectionSecondMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate rib cap section second moment of area.
subroutine RibCapSectionSecondMomentOfArea(yc,Ixx)

    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Ixx					! section second moment of inertia
    !Local Variables

    !k	= k
    Ixx	= yc*0.0D0

end subroutine RibCapSectionSecondMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate rib section second moment of area.
subroutine RibSectionSecondMomentOfArea(yc,Ixx)

    implicit none

    !Input variables
    real(8), intent(in) ::              yc					! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Ixx					! section second moment of inertia
    !Local Variables
    real(8) ::				Ixxrw				! rib web section second moment of area
    real(8) ::				Ixxrc				! rib cap section second moment of area

    call RibWebSectionSecondMomentOfArea(yc,Ixxrw)
    call RibCapSectionSecondMomentOfArea(yc,Ixxrc)
    Ixx	= Ixxrw+Ixxrc

end subroutine RibSectionSecondMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate stringer section second moment of area.
subroutine StringerSectionSecondMomentOfArea(yc,Ixx)
    use PlanformVariables, only: Ntr
    use StringerVariablesRoutines
    use IntegrationVariables
    USE maths, only: AxesTranformation,GaussQuadratureWeightsSamplingPoints

    implicit none

    !Input variables
    real(8), intent(in) ::              yc				! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Ixx				! section second moment of inertia
    !Local Variables
    real(8) ::                          G				! function to integrate
    real(8) ::                          x,y				! point coordinates
    real(8) ::                          eta,zeta			! point coordinates
    real(8) ::                          x1,x2				! limits of integration
    real(8) ::                          z1U,z2U,z1L,z2L		! limits of integration
    integer ::                          n,m_g				! iteration variables
    real(8),allocatable ::              gi(:),zetai(:)   

    Ixx	= 0.0D0
    if(Stringer(Ntr)%Nst.EQ.0) return
    allocate(gi(Mg),zetai(Mg))
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0+2.0D0*yc
    do m_g=1,Mg
        do n=1,Stringer(Ntr)%Nst
            !call StringerGeometry1(n,eta)											! compute Stringer(Ntr)%zetast
            call StringerGeometry2(1,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute Stringer(Ntr)%lcst, zeta
            call AxesTranformation(1,zeta,eta,x,y)								! compute x, y
            x1	= x-(Stringer(Ntr)%lst1(n)-Stringer(Ntr)%lst2(n)*eta)/2.0D0
            x2	= x+(Stringer(Ntr)%lst1(n)-Stringer(Ntr)%lst2(n)*eta)/2.0D0
            call StringerGeometry2(2,n,zetai(m_g),eta,zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
            G	= (x2-x1)/2.0D0*((z2U*z2U*z2U-z1U*z1U*z1U)+(z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0
            Ixx	= Ixx+Stringer(Ntr)%Eest(n)*gi(m_g)*G
        end do
    end do
    deallocate(gi,zetai)

end subroutine StringerSectionSecondMomentOfArea

!--------------------------------------------------------------------------------------------
! Routine to calculate overall section second moment of area.
subroutine SectionSecondMomentOfArea(screen,yc,Ixx)

    use PlanformVariables, only: Ntr,b2,yroot
    use SkinVariablesRoutines, only: Skin

    implicit none

    !Input variables
    logical, intent(in) ::              screen          ! TRUE: write info to scree; False otherwise
    real(8), intent(in) ::              yc				! fraction of semi-span
    !Output variables
    real(8), intent(out) ::             Ixx				! section second moment of inertia
    !Local Variables
    real(8) ::				            Ixxsk           ! skin section second moment of area
    real(8) ::				            Ixxsp           ! spar section second moment of area
    real(8) ::				            Ixxrb           ! rib section second moment of area
    real(8) ::				            Ixxst           ! stringer section second moment of area
    real(8) ::                          E0              ! reference elastic modulus
    
    if(screen) write(*,*) 'Computing section second moment of area at y= ',yroot(Ntr) + yc*b2(Ntr),' [m]'
    
    E0	= Skin(Ntr)%Eysk(1)
    write(*,*) '    - skin'
    call SkinSectionSecondMomentOfArea(yc,Ixxsk)
    if(screen) write(*,'(1X,A8,F12.9,A3)') 'Ixxsk = ',Ixxsk/E0,' [m^4]'

    write(*,*) '    - spar'
    call SparSectionSecondMomentOfArea(yc,Ixxsp)
    if(screen) write(*,'(1X,A8,F12.9,A3)') 'Ixxsp = ',Ixxsp/E0,' [m^4]'

    write(*,*) '    - rib'
    call RibSectionSecondMomentOfArea(yc,Ixxrb)
    if(screen) write(*,'(1X,A8,F12.9,A3)') 'Ixxrb = ',Ixxrb/E0,' [m^4]'

    write(*,*) '    - stringer'
    call StringerSectionSecondMomentOfArea(yc,Ixxst)
    if(screen) write(*,'(1X,A8,F12.9,A3)') 'Ixxst = ',Ixxst/E0,' [m^4]'

    Ixx	= (Ixxsk+Ixxsp+Ixxrb+Ixxst)/E0 
    write(*,'(1X,A7,F12.9,A3)') 'Ixx  = ',Ixx,' [m^4]'

end subroutine SectionSecondMomentOfArea

