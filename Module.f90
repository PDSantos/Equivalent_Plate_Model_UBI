!--------------------------------------------------------------------------------------------
! Module with planform shape variables.
module PlanformVariables
    implicit none
    save

    integer ::                Npl		            ! number of wing trapezoids
    real(8),allocatable ::    yroot(:)		    ! y-coordinate of root chord
    real(8),allocatable ::    ytip(:)		        ! y-coordinate of tip chord
    real(8),allocatable ::    xLEroot(:)	        ! x-coordinate of root chord leading edge
    real(8),allocatable ::    xLEtip(:)		    ! x-coordinate of tip chord leading edge
    real(8),allocatable ::    zroot_pl(:)	        ! z-coordinate of root chord
    real(8),allocatable ::    ztip_pl(:)	        ! z-coordinate of tip chord
    real(8),allocatable ::    thetaroot(:)     	! incidence of root chord
    real(8),allocatable ::    thetatip(:)	        ! incidence of tip chord
    real(8),allocatable ::    croot(:)		    ! root chord
    real(8),allocatable ::    ctip(:)		        ! tip chord
    logical ::                boundary	        ! boudary conditions .TRUE.=built in, .FALSE.=simply supported
    real(8),allocatable ::    b2(:)		        ! semi-span
    real(8),allocatable ::    Lamdapl(:)	        ! leading edge sweep angle
    real(8),allocatable ::    mcpl(:)		        ! chord variation in y
    integer ::                Ntr		            ! number of the wing trapezoid beeing computed
    TYPE :: SectionVarsType
        !Spline section variables
        integer ::              n_sections
        real(8),allocatable ::  x_sections(:)           ! section x array
        real(8),allocatable ::  x_sections_centroid(:)  ! section centroid array      
        real(8),allocatable ::  y_sections(:)           ! section y array
        real(8),allocatable ::  z_sections(:)           ! section z array
        real(8),allocatable ::  z_sections_centroid(:)  ! section centroid array      
        real(8),allocatable ::  c_sections(:)           ! section chord array
        integer,allocatable ::  n_sec_airfoil(:)
    END TYPE SectionVarsType
    TYPE(SectionVarsType), ALLOCATABLE :: Section(:)
    TYPE :: TransMatrixType
        real(8) ::              alpha = 0.0D0       ! 
        real(8) ::              betta = 0.0D0       ! 
        real(8) ::              gamma = 0.0D0       ! 
        real(8) ::              Rot(3,3)            !Rotation matrix                                                      
    END TYPE TransMatrixType
    TYPE(TransMatrixType), ALLOCATABLE :: TransMatrix(:) !
contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all planform variables
    subroutine AllocatePlanformVariables()
      
        allocate(yroot(npl),ytip(npl),xLEroot(npl),xLEtip(npl),zroot_pl(npl),ztip_pl(npl))
        allocate(thetaroot(npl),thetatip(npl),croot(npl),ctip(npl))
        allocate(b2(npl),Lamdapl(npl),mcpl(npl))
        !Tranformation variables
        allocate(TransMatrix(npl))
	
    end subroutine AllocatePlanformVariables
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all planform variables
    subroutine AllocateSectionVariables(Nsk,nl_section)
        !Input vars
        integer, intent(in) ::		Nsk				! number of skin(s) for current planform
        integer, intent(in) ::		nl_section      ! number of sections for current planform
        if(.not. allocated(Section)) then !First planform
            allocate(Section(Npl))             
        end if
        Section(Ntr)%n_sections = nl_section
        allocate(Section(Ntr)%n_sec_airfoil(Nsk),Section(Ntr)%x_sections(nl_section),Section(Ntr)%x_sections_centroid(nl_section),Section(Ntr)%z_sections(nl_section),Section(Ntr)%y_sections(nl_section),Section(Ntr)%c_sections(nl_section),Section(Ntr)%z_sections_centroid(nl_section))
	
    end subroutine AllocateSectionVariables    
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all planform variables
    subroutine ComputeRotationMatrix(Ntr)
        !Input vars
        integer, intent(in) ::		Ntr				! current planform number
        
        TransMatrix(Ntr)%Rot(1,1) = cos(TransMatrix(Ntr)%gamma) *cos(TransMatrix(Ntr)%alpha) - cos(TransMatrix(Ntr)%betta)*sin(TransMatrix(Ntr)%alpha)*sin(TransMatrix(Ntr)%gamma)
        TransMatrix(Ntr)%Rot(1,2) = -cos(TransMatrix(Ntr)%alpha) *sin(TransMatrix(Ntr)%gamma) - sin(TransMatrix(Ntr)%alpha)*cos(TransMatrix(Ntr)%betta)*cos(TransMatrix(Ntr)%gamma)
        TransMatrix(Ntr)%Rot(1,3) = sin(TransMatrix(Ntr)%alpha)*sin(TransMatrix(Ntr)%betta)
        TransMatrix(Ntr)%Rot(2,1) = sin(TransMatrix(Ntr)%alpha) *cos(TransMatrix(Ntr)%gamma) + cos(TransMatrix(Ntr)%alpha)*cos(TransMatrix(Ntr)%betta)*sin(TransMatrix(Ntr)%gamma)
        TransMatrix(Ntr)%Rot(2,2) = cos(TransMatrix(Ntr)%alpha) *cos(TransMatrix(Ntr)%betta)*cos(TransMatrix(Ntr)%gamma) - sin(TransMatrix(Ntr)%alpha)*sin(TransMatrix(Ntr)%gamma)
        TransMatrix(Ntr)%Rot(2,3) = -cos(TransMatrix(Ntr)%alpha)*sin(TransMatrix(Ntr)%betta)
        TransMatrix(Ntr)%Rot(3,1) = sin(TransMatrix(Ntr)%betta)*sin(TransMatrix(Ntr)%gamma) 
        TransMatrix(Ntr)%Rot(3,2) = sin(TransMatrix(Ntr)%betta)*cos(TransMatrix(Ntr)%gamma)
        TransMatrix(Ntr)%Rot(3,3) = cos(TransMatrix(Ntr)%betta)
	
    end subroutine ComputeRotationMatrix  
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all planform variables   
    subroutine ApplyScalePlanform(scale)
        !Input vars
        real(8), intent(in) ::		scale				! current planform number
        !Local vars
        integer ::                  n
        forall(n=1:Npl)
            yroot(n) = yroot(n) * scale 
            xLEroot(n)= xLEroot(n) * scale 
            zroot_pl(n)= zroot_pl(n) * scale 
            croot(n)= croot(n) * scale 
            ytip(n) = ytip(n) * scale 
            xLEtip(n)= xLEtip(n) * scale 
            ztip_pl(n)= ztip_pl(n) * scale 
            ctip(n)= ctip(n) * scale 
        end forall
        
    end subroutine ApplyScalePlanform

    
end module PlanformVariables

!--------------------------------------------------------------------------------------------
! Module with skin variables.
module SkinVariablesRoutines
    implicit none
    save
    TYPE :: SkinPolVarsType
        integer ::		            Nt			    ! no. of polynomial terms defining skin thickness
        real(8), allocatable ::		Tz(:)           ! coefficients of polynomial defining thickness
        integer, allocatable ::		n_t(:)	        ! powers of polynomial defining thickness
    END TYPE SkinpolVarsType
    TYPE :: SkinSectionVarsType !xcai_section(Skin(Ntr)%Nsk,n_max,Section(Ntr)%n_sections)
        real(4),allocatable ::		XS_COEF(:,:)
        real(4),allocatable ::		xcai(:,:)
        real(4),allocatable ::		zcai(:,:) 
    END TYPE SkinSectionVarsType
    TYPE :: SkinVarsType
        integer ::            		Nsk				! number of skin panels
        real(8), allocatable ::		kskFroot(:)	    ! chord fraction of front line at root
        real(8), allocatable ::		kskFtip(:)		! chord fraction of front line at tip
        real(8), allocatable ::		kskAroot(:)	    ! chord fraction of aft line at root
        real(8), allocatable ::		kskAtip(:)		! chord fraction of aft line at tip
        real(8), allocatable ::		kskF1(:)		! front line equation coefficient
        real(8), allocatable ::		kskF2(:)		! front line equation coefficient
        real(8), allocatable ::		kskA1(:)		! aft line equation coefficient
        real(8), allocatable ::		kskA2(:)		! aft line equation coefficient

        TYPE(SkinPolVarsType), ALLOCATABLE :: Poly(:)
        
        real(8), allocatable ::		t0sk(:)		        ! panel thickness
        real(8), allocatable ::		ttiptroot(:)        ! ratio of tip thickness to root thickness
        real(8), allocatable ::		mtsk(:)		        ! thickness variation in y

        real(8), allocatable ::		rhosk(:)	        ! material density

        real(8), allocatable ::		Exsk(:)		    ! longitudinal elastic modulus in x
        real(8), allocatable ::		Eysk(:)		    ! longitudinal elastic modulus in y
        real(8), allocatable ::		Gxysk(:)		! transverse elastic modulus
        real(8), allocatable ::		miuxysk(:)		! poisson ratio xy
        real(8), allocatable ::		miuyxsk(:)		! poisson ratio yx
        real(8), allocatable ::		sk11(:)		    ! constitutive matrix element
        real(8), allocatable ::		sk12(:)		    ! constitutive matrix element
        real(8), allocatable ::		sk21(:)		    ! constitutive matrix element
        real(8), allocatable ::		sk22(:)		    ! constitutive matrix element
        real(8), allocatable ::		sk44(:)		    ! constitutive matrix element

        real(8), allocatable ::		Gxzsk(:)		! transverse elastic modulus
        real(8), allocatable ::		Gyzsk(:)		! transverse elastic modulus
        real(8), allocatable ::		sk55(:)		    ! constitutive matrix element
        real(8), allocatable ::		sk66(:)		    ! constitutive matrix element
        
        real(8),allocatable ::		D11sk(:),D12sk(:),D14sk(:),D15sk(:),D16sk(:)		! constitutive matrix element
        real(8),allocatable ::		D21sk(:),D22sk(:),D24sk(:),D25sk(:),D26sk(:)		! constitutive matrix element
        real(8),allocatable ::		D41sk(:),D42sk(:),D44sk(:),D45sk(:),D46sk(:)		! constitutive matrix element
        real(8),allocatable ::		D51sk(:),D52sk(:),D54sk(:),D55sk(:),D56sk(:)		! constitutive matrix element
        real(8),allocatable ::		D61sk(:),D62sk(:),D64sk(:),D65sk(:),D66sk(:)		! constitutive matrix element
        TYPE(SkinSectionVarsType), ALLOCATABLE :: Section(:)
  END TYPE SkinVarsType
  TYPE(SkinVarsType), ALLOCATABLE :: Skin(:)
  
contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all skin variables
    subroutine AllocateSkinVariables(Npl,Ntr,Nsk)
        !Input vars
        integer, intent(in) ::		Npl				! number of planforms
        integer, intent(in) ::		Ntr				! current planform number
        integer, intent(in) ::		Nsk				! number of skin(s) for current planform
        if(.not. allocated(Skin)) then !First planform
            allocate(Skin(Npl))             
        end if
        Skin(Ntr)%Nsk = Nsk
        If ( Nsk > 0) then
            allocate(Skin(Ntr)%kskFroot(Nsk),Skin(Ntr)%kskFtip(Nsk),Skin(Ntr)%kskAroot(Nsk),Skin(Ntr)%kskAtip(Nsk),Skin(Ntr)%kskF1(Nsk),Skin(Ntr)%kskF2(Nsk),Skin(Ntr)%kskA1(Nsk),Skin(Ntr)%kskA2(Nsk))
            allocate(Skin(Ntr)%Exsk(Nsk),Skin(Ntr)%Eysk(Nsk),Skin(Ntr)%Gxysk(Nsk),Skin(Ntr)%miuxysk(Nsk),Skin(Ntr)%miuyxsk(Nsk),Skin(Ntr)%sk11(Nsk),Skin(Ntr)%sk12(Nsk),Skin(Ntr)%sk21(Nsk))
            allocate(Skin(Ntr)%sk22(Nsk),Skin(Ntr)%sk44(Nsk),Skin(Ntr)%Gxzsk(Nsk),Skin(Ntr)%Gyzsk(Nsk),Skin(Ntr)%sk55(Nsk),Skin(Ntr)%sk66(Nsk))               

            allocate(Skin(Ntr)%t0sk(nsk),Skin(Ntr)%ttiptroot(Nsk),Skin(Ntr)%mtsk(Nsk),Skin(Ntr)%rhosk(nsk))
            allocate(Skin(Ntr)%Poly(nsk))

            allocate(Skin(Ntr)%D11sk(Nsk),Skin(Ntr)%D12sk(Nsk),Skin(Ntr)%D14sk(Nsk),Skin(Ntr)%D15sk(Nsk),Skin(Ntr)%D16sk(Nsk),Skin(Ntr)%D21sk(Nsk),Skin(Ntr)%D22sk(Nsk),Skin(Ntr)%D24sk(Nsk),Skin(Ntr)%D25sk(Nsk),Skin(Ntr)%D26sk(Nsk))
            allocate(Skin(Ntr)%D41sk(Nsk),Skin(Ntr)%D42sk(Nsk),Skin(Ntr)%D44sk(Nsk),Skin(Ntr)%D45sk(Nsk),Skin(Ntr)%D46sk(Nsk),Skin(Ntr)%D51sk(Nsk),Skin(Ntr)%D52sk(Nsk),Skin(Ntr)%D54sk(Nsk),Skin(Ntr)%D55sk(Nsk),Skin(Ntr)%D56sk(Nsk))
            allocate(Skin(Ntr)%D61sk(Nsk),Skin(Ntr)%D62sk(Nsk),Skin(Ntr)%D64sk(Nsk),Skin(Ntr)%D65sk(Nsk),Skin(Ntr)%D66sk(Nsk))
            
            allocate(Skin(Ntr)%Section(Nsk))
        end if  
    end subroutine AllocateSkinVariables
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all skin variables    
    subroutine DeallocateSkinVariables() 
        
        deallocate(Skin)
   
    end subroutine DeallocateSkinVariables
    !--------------------------------------------------------------------------------------------
    ! Routine to calculate skin constitutive matrix.
    subroutine SkinConstitutiveMatrix(flag,theta,Ntr,nskin)
        ! Input variables.
        integer, intent(in) ::	            flag    ! calculation option: 0 = initialize matrix;
                                                    !                     1 = calculate matrix
        real(8), intent(in) ::	            theta   ! skin inclination angle
        integer, intent(in) ::	            Ntr     ! current planform
        integer, intent(in), optional ::	nskin   ! current computed skin

        ! Local variables.
        real(8) ::              denominator         !
        integer ::              n                   ! iteration variables

        if(flag.EQ.0) then    
            do n=1,Skin(Ntr)%Nsk
                denominator	= 1.0D0-Skin(Ntr)%miuxysk(n)*Skin(Ntr)%miuyxsk(n)
                Skin(Ntr)%sk11(n)	= Skin(Ntr)%Exsk(n)/denominator
                Skin(Ntr)%sk12(n)	= Skin(Ntr)%miuyxsk(n)*Skin(Ntr)%Exsk(n)/denominator
                Skin(Ntr)%sk21(n)	= Skin(Ntr)%miuxysk(n)*Skin(Ntr)%Eysk(n)/denominator
                Skin(Ntr)%sk22(n)	= Skin(Ntr)%Eysk(n)/denominator
                Skin(Ntr)%sk44(n)	= Skin(Ntr)%Gxysk(n)
                Skin(Ntr)%sk55(n)	= 0.0D0 !5/6*Skin(Ntr)%Gxzsk(n)
                Skin(Ntr)%sk66(n)	= 0.0D0 !5/6*Skin(Ntr)%Gyzsk(n)
            end do
        else if(flag.EQ.1) then
            n = nskin
            Skin(Ntr)%D11sk(n)	= Skin(Ntr)%sk11(n)*dcos(theta)*dcos(theta)*dcos(theta)*dcos(theta)+Skin(Ntr)%sk66(n)*dsin(2*theta)*dsin(2*theta)
            Skin(Ntr)%D12sk(n)	= Skin(Ntr)%sk12(n)*dcos(theta)*dcos(theta)
            Skin(Ntr)%D14sk(n)	= 0.0D0
            Skin(Ntr)%D15sk(n)	= 0.0D0
            Skin(Ntr)%D16sk(n)	= Skin(Ntr)%sk11(n)*dsin(theta)*dcos(theta)*dcos(theta)*dcos(theta)-Skin(Ntr)%sk66(n)*dsin(2.0D0*theta)*dcos(2.0D0*theta)
            Skin(Ntr)%D21sk(n)	= Skin(Ntr)%D12sk(n)
            Skin(Ntr)%D22sk(n)	= Skin(Ntr)%sk22(n)
            Skin(Ntr)%D24sk(n)	= 0.0D0
            Skin(Ntr)%D25sk(n)	= 0.0D0 
            Skin(Ntr)%D26sk(n)	= Skin(Ntr)%sk21(n)*dsin(theta)*dcos(theta)
            Skin(Ntr)%D41sk(n)	= Skin(Ntr)%D14sk(n)
            Skin(Ntr)%D42sk(n)	= Skin(Ntr)%D24sk(n)
            Skin(Ntr)%D44sk(n)	= Skin(Ntr)%sk44(n)*dcos(theta)*dcos(theta)+Skin(Ntr)%sk55(n)*dsin(theta)*dsin(theta)
            Skin(Ntr)%D45sk(n)	= Skin(Ntr)%sk44(n)*dcos(theta)*dsin(theta)-Skin(Ntr)%sk55(n)*dcos(theta)*dsin(theta) ! ao rui dá + a ultima parcela e - a primeira
            Skin(Ntr)%D46sk(n)	= 0.0D0 
            Skin(Ntr)%D51sk(n)	= Skin(Ntr)%D15sk(n)
            Skin(Ntr)%D52sk(n)	= Skin(Ntr)%D25sk(n)
            Skin(Ntr)%D54sk(n)	= Skin(Ntr)%D45sk(n)
            Skin(Ntr)%D55sk(n)	= Skin(Ntr)%sk44(n)*dsin(theta)*dsin(theta)+Skin(Ntr)%sk55(n)*dcos(theta)*dcos(theta)
            Skin(Ntr)%D56sk(n)	= 0.0D0
            Skin(Ntr)%D61sk(n)	= Skin(Ntr)%D16sk(n)
            Skin(Ntr)%D62sk(n)	= Skin(Ntr)%D26sk(n)
            Skin(Ntr)%D64sk(n)	= Skin(Ntr)%D46sk(n)
            Skin(Ntr)%D65sk(n)	= Skin(Ntr)%D56sk(n)
            Skin(Ntr)%D66sk(n)	= Skin(Ntr)%sk11(n)*dsin(theta)*dsin(theta)*dcos(theta)*dcos(theta)-Skin(Ntr)%sk66(n)*dcos(2.0D0*theta)*dcos(2.0D0*theta) !Rui tem +
        end if
	
    end subroutine SkinConstitutiveMatrix
    !--------------------------------------------------------------------------------------------
    ! Routine to calculate zeta coodinate for geometry export to tecplot. - cossine distribution
    elemental real(8) function SkinZetaTecplot(n,ni,i) result(zeta) 
        USE maths, only: Pi
        !Input variables  
        integer, intent(in) ::		n					! skin number
        integer, intent(in) ::		ni					! total number of chordwise points
        integer, intent(in) ::		i					! current chordwise point
        
        if(n.EQ.1) then
            zeta	= -cos(Pi*dfloat(i)/dfloat(ni))
        else
            zeta	= cos(Pi*dfloat(i)/dfloat(ni))
        end if
    end function SkinZetaTecplot
    !--------------------------------------------------------------------------------------------
    ! Routine to calculate eta coodinate for geometry export to tecplot. - linear varying
    elemental real(8) function SkinEtaTecplot(nj,j) result(eta) 
        !Input variables  
        integer, intent(in) ::		nj					! total number of chordwise points
        integer, intent(in) ::		j					! current spanwise point
        
        eta = -1.0D0+2.0D0*float(j)/float(nj)
        
    end function SkinEtaTecplot
    !--------------------------------------------------------------------------------------------
    ! Routine to calculate skin geometry.
    subroutine SkinGeometry1(n,etaj,zetaskF,zetaskA)
        use PlanformVariables, only: Ntr
        !Input variables  
        integer, intent(in) ::		n					! skin number
        real(8), intent(in) ::		etaj				! input eta-coordinate
        !Output variables
        real(8), intent(out) ::		zetaskF,zetaskA
        
        zetaskF	= -1.0D0+2.0D0*(Skin(Ntr)%kskF1(n)-Skin(Ntr)%kskF2(n)*etaj)
        zetaskA	= -1.0D0+2.0D0*(Skin(Ntr)%kskA1(n)-Skin(Ntr)%kskA2(n)*etaj)

    end subroutine SkinGeometry1

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate skin geometry.
    subroutine SkinGeometry2(flag,n,x,y,z,z1,z2,slope)

        use PlanformVariables, only: Ntr,xLEroot,Lamdapl,yroot,croot,mcpl
  
        !Input variables
        integer, intent(in) ::			    flag	! flag to choose calculation
										            !   1 = calculate z given (x,y)
										            !   2 = calculate t given (x,y)
										            !   3 = calculate limits of integration (and slope)
        integer, intent(in) ::			    n		! number of skin panel
        real(8), intent(in) ::			    x,y 	! point coordinates
        !Output variables
        real(8), intent(out) ::			    z	    ! point z coordinate
        real(8),optional, intent(out) ::	z1,z2	! limits of integration 
        real(8),optional, intent(out) ::	slope	! slope of skinheight dz/dx
        !Local variables
        real(8) ::        				    xLE		! leading edge x-coordinate
        real(8) ::        				    c		! chord
        real(8) ::        				    xc		! x to chord ratio
        real(8) ::        				    tc		! thickness to chord ratio
        real(8) ::        				    t		! vertical thickness of skin
        integer ::        				    i		! iteration variable
        real(8) ::        				    tspan	! spanwise thickness correction
  
        if(flag.EQ.1.OR.flag.EQ.2.OR.flag.EQ.3) then
            xLE		= xLEroot(Ntr)+tan(Lamdapl(Ntr))*(y-yroot(Ntr))
            c		= croot(Ntr)+mcpl(Ntr)*(y-yroot(Ntr))
            xc		= (x-xLE)/c
            tc	    = 1.0D0
            z		= SkinHeightSlope(n,xc,y,1) 
        end if
        if(flag.EQ.2.OR.flag.EQ.3) then
            Skin(Ntr)%t0sk(n)	= 0.0D0
            do i=1,Skin(Ntr)%Poly(n)%Nt
	            Skin(Ntr)%t0sk(n)	= Skin(Ntr)%t0sk(n)+Skin(Ntr)%Poly(n)%Tz(i)*xc**Skin(Ntr)%Poly(n)%n_t(i)
            end do
            slope	= SkinHeightSlope(n,xc,y,2) !Skin slope
        end if
        if(flag.EQ.3) then
            tspan	= 1.0D0 + Skin(Ntr)%mtsk(n)*(y-yroot(Ntr)) !Skin thickeness variation
            t		= Skin(Ntr)%t0sk(n)*dsqrt(1.0D0+slope*slope)*tspan ! ???
            z1		= z-0.5D0*t
            z2		= z+0.5D0*t
        end if

    end subroutine SkinGeometry2

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate skin height or slope. Modified by Pedro Santos (17/11/2010) - ONLY DATA SET IS ACCEPTED!!!
    real(8) function SkinHeightSlope(n,xc,y,flag)

        use PlanformVariables, only: Ntr,b2,yroot,Section,yroot
        use SplineMaths, only: evaluate_spline
  
        !Input variables
        integer, intent(in) ::		n				! number of skin panel
        real(8), intent(in) ::		xc				! x to chord ratio
        real(8), intent(in) ::		y				! spanwise position
        integer, intent(in) ::		flag			! flag to choose calculation:
											        !   1 = height
											        !   2 = slope
        !Local variables
        integer ::    				i,iref					            ! iteration variable
        real(4) ::                  z_root_s(1),z_tip_s(1),xc_s(1)      ! single precision version of above mentioned vars
        real(8) ::    				z_root(1),z_tip(1)                  ! point from spline interpolation
        real(8) ::                  eps=epsilon(y)*10.0D0

        !Search y stations to interpolate
        iref=-1 !just to help debug
        do i=2,Section(Ntr)%n_sections
            IF(y - yroot(Ntr) >= Section(Ntr)%y_sections(i-1)-eps .AND. y - yroot(Ntr) <= Section(Ntr)%y_sections(i)+eps) Then !sum with eps to avoid numerical error when evaluating logical <=
                iref=i
                exit
            end if
        end do
        if(iref == -1) then
            continue
            iref = Section(Ntr)%n_sections
        end if

        xc_s(1)=real(xc)
        Call evaluate_spline(Section(Ntr)%n_sec_airfoil(n),real(Skin(Ntr)%Section(n)%xcai(1:Section(Ntr)%n_sec_airfoil(n),iref-1)),real(Skin(Ntr)%Section(n)%zcai(1:Section(Ntr)%n_sec_airfoil(n),iref-1)),Skin(Ntr)%Section(n)%XS_Coef(1:Section(Ntr)%n_sec_airfoil(n),iref-1),flag-1,1,xc_s,z_root_s)
        z_root(1)=z_root_s(1)
        Call evaluate_spline(Section(Ntr)%n_sec_airfoil(n),real(Skin(Ntr)%Section(n)%xcai(1:Section(Ntr)%n_sec_airfoil(n),iref)),real(Skin(Ntr)%Section(n)%zcai(1:Section(Ntr)%n_sec_airfoil(n),iref)),Skin(Ntr)%Section(n)%XS_Coef(1:Section(Ntr)%n_sec_airfoil(n),iref),flag-1,1,xc_s,z_tip_s)
        z_tip(1)=z_tip_s(1)
        SkinHeightSlope	= z_root(1)+(z_tip(1)-z_root(1))/(Section(Ntr)%y_sections(iref)-Section(Ntr)%y_sections(iref-1))*(y - yroot(Ntr) - Section(Ntr)%y_sections(iref-1))
    end function SkinHeightSlope
    !--------------------------------------------------------------------------------------------
    ! Routine to apply scale to skin   
    subroutine ApplyScaleSkin(scale)
        use PlanformVariables, only: Ntr
        !Input vars
        real(8), intent(in) ::		scale				! current planform number
        !Local vars
        integer ::                  n
        
        forall(n=1:Skin(Ntr)%Nsk)
            Skin(Ntr)%Poly(n)%Tz(1) = Skin(Ntr)%Poly(n)%Tz(1) * scale 
            Skin(Ntr)%rhosk(n) = Skin(Ntr)%rhosk(n) / scale**3
            Skin(Ntr)%Exsk(n) = Skin(Ntr)%Exsk(n) / scale
            Skin(Ntr)%Eysk(n) = Skin(Ntr)%Eysk(n) / scale
            Skin(Ntr)%Gxysk(n) = Skin(Ntr)%Gxysk(n) / scale
            Skin(Ntr)%Gxzsk(n) = Skin(Ntr)%Gxzsk(n) / scale
            Skin(Ntr)%Gyzsk(n) = Skin(Ntr)%Gyzsk(n) / scale
        end forall
        
    end subroutine ApplyScaleSkin

end module SkinVariablesRoutines

!--------------------------------------------------------------------------------------------
! Module with spar variables.
module SparVariablesRoutines
  implicit none
  save

    TYPE :: SparVarsType
      integer ::            				Nsp					! number of spars
      real(8), allocatable ::				ksproot(:)		    ! chord fraction of spar position at root
      real(8), allocatable ::				ksptip(:)		    ! chord fraction of spar position at tip
      real(8), allocatable ::				tswroot(:)		    ! spar web thickness at root
      real(8), allocatable ::				tswtip(:)		    ! spar web thickness at tip
      real(8), allocatable ::				hscroot(:)		    ! spar cap height at root
      real(8), allocatable ::				hsctip(:)		    ! spar cap height at tip
      real(8), allocatable ::				lscroot(:)		    ! spar cap width at root
      real(8), allocatable ::				lsctip(:)		    ! spar cap width at tip
      real(8), allocatable ::				ksp1(:)		        ! spar line equation coefficient
      real(8), allocatable ::				ksp2(:)		        ! spar line equation coefficient
      real(8), allocatable ::				tsw1(:)		        ! spar web thickness equation coefficient
      real(8), allocatable ::				tsw2(:)		        ! spar web thickness equation coefficient
      real(8), allocatable ::				hsc1(:)		        ! spar cap height equation coefficient
      real(8), allocatable ::				hsc2(:)		        ! spar cap height equation coefficient
      real(8), allocatable ::				lsc1(:)		        ! spar cap width equation coefficient
      real(8), allocatable ::				lsc2(:)		        ! spar cap width equation coefficient
      real(8), allocatable ::				Lamdasp(:)		    ! spar sweep angle

      real(8), allocatable ::				rhosw(:)		    ! spar web material density
      real(8), allocatable ::				rhosc(:)		    ! spar cap material density

      real(8) ::            				tcsw				! spar web thickness to wing chord ratio
      real(8) ::            				lcsc				! spar cap width to wing chord ratio

      real(8), allocatable ::				Eesc(:)		        ! longitudinal elastic modulus in x
      real(8), allocatable ::				miuezsc(:)		    ! poisson ratio xy
      real(8), allocatable ::				miuzesc(:)		    ! poisson ratio yx
      real(8), allocatable ::				D11sc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D12sc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D13sc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D21sc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D22sc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D23sc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D31sc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D32sc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D33sc(:)		    ! constitutive matrix element

      real(8), allocatable ::				Gezsw(:)		    ! transverse elastic modulus
      real(8), allocatable ::				Eesw(:)		        ! longitudinal elastic modulus in x
      real(8), allocatable ::				miuezsw(:)		    ! poisson ratio xy
      real(8), allocatable ::				miuzesw(:)		    ! poisson ratio yx
      real(8), allocatable ::				D11sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D12sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D13sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D21sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D22sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D23sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D31sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D32sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D33sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D44sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D45sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D54sw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D55sw(:)		    ! constitutive matrix element
  END TYPE SparVarsType
  TYPE(SparVarsType), ALLOCATABLE :: Spar(:)
  
contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all spar variables
    subroutine AllocateSparVariables(Npl,Ntr,Nsp)
        !Input vars
        integer, intent(in) ::		Npl				! number of planforms
        integer, intent(in) ::		Ntr				! current planform number
        integer, intent(in) ::		Nsp				! number spar(s) for current planform
        if(.not. allocated(Spar)) then !First planform
            allocate(Spar(Npl))             
        end if
        Spar(Ntr)%Nsp = Nsp
        If ( Nsp > 0) then
            allocate(Spar(Ntr)%ksproot(Nsp),Spar(Ntr)%ksptip(Nsp),Spar(Ntr)%tswroot(Nsp),Spar(Ntr)%tswtip(Nsp),Spar(Ntr)%hscroot(Nsp),Spar(Ntr)%hsctip(Nsp),Spar(Ntr)%lscroot(Nsp),Spar(Ntr)%lsctip(Nsp))
            allocate(Spar(Ntr)%ksp1(Spar(Ntr)%Nsp),Spar(Ntr)%ksp2(Nsp),Spar(Ntr)%tsw1(Nsp),Spar(Ntr)%tsw2(Nsp),Spar(Ntr)%hsc1(Nsp),Spar(Ntr)%hsc2(Nsp),Spar(Ntr)%lsc1(Nsp),Spar(Ntr)%lsc2(Nsp),Spar(Ntr)%Lamdasp(Nsp),Spar(Ntr)%rhosw(Nsp),Spar(Ntr)%rhosc(Nsp))
            allocate(Spar(Ntr)%Eesc(Nsp),Spar(Ntr)%miuezsc(Nsp),Spar(Ntr)%miuzesc(Nsp),Spar(Ntr)%D11sc(Nsp),Spar(Ntr)%D12sc(Nsp),Spar(Ntr)%D13sc(Nsp),Spar(Ntr)%D21sc(Nsp),Spar(Ntr)%D22sc(Nsp),Spar(Ntr)%D23sc(Nsp),Spar(Ntr)%D31sc(Nsp),Spar(Ntr)%D32sc(Nsp),Spar(Ntr)%D33sc(Nsp))
            allocate(Spar(Ntr)%Gezsw(Nsp),Spar(Ntr)%Eesw(Nsp),Spar(Ntr)%miuezsw(Nsp),Spar(Ntr)%miuzesw(Nsp),Spar(Ntr)%D11sw(Nsp),Spar(Ntr)%D12sw(Nsp),Spar(Ntr)%D13sw(Nsp),Spar(Ntr)%D21sw(Nsp),Spar(Ntr)%D22sw(Nsp),Spar(Ntr)%D23sw(Nsp),Spar(Ntr)%D31sw(Nsp))
            allocate(Spar(Ntr)%D32sw(Nsp),Spar(Ntr)%D33sw(Nsp),Spar(Ntr)%D44sw(Nsp),Spar(Ntr)%D45sw(Nsp),Spar(Ntr)%D54sw(Nsp),Spar(Ntr)%D55sw(Nsp))
        end if  
        
    end subroutine AllocateSparVariables    
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all spar variables
    subroutine DeallocateSparVariables()
    
        deallocate(Spar)

    end subroutine DeallocateSparVariables

!--------------------------------------------------------------------------------------------
    ! Routine to calculate spar constitutive matrix.
    subroutine SparConstitutiveMatrix(Ntr)
        !Input vars
        integer, intent(in) ::	Ntr				    ! current planform
        !local vars
        real(8) ::				L11,L12,L13,L24,L25	! direction cosines
        real(8) ::				E11            	    ! stiffness element
        integer ::				n					! iteration variables

        do n=1,Spar(Ntr)%Nsp
	        L11			= dsin(Spar(Ntr)%Lamdasp(n))*dsin(Spar(Ntr)%Lamdasp(n))
	        L12			= dcos(Spar(Ntr)%Lamdasp(n))*dcos(Spar(Ntr)%Lamdasp(n))
	        L13			= dsin(Spar(Ntr)%Lamdasp(n))*dcos(Spar(Ntr)%Lamdasp(n))				
	        L24			= dsin(Spar(Ntr)%Lamdasp(n))								
	        L25			= dcos(Spar(Ntr)%Lamdasp(n))								
	        ! Spar web.
	        E11 = Spar(Ntr)%Eesw(n)/(1.0D0-Spar(Ntr)%miuezsw(n)*Spar(Ntr)%miuzesw(n))
	        Spar(Ntr)%D11sw(n)	= L11*E11*L11
	        Spar(Ntr)%D12sw(n)	= L11*E11*L12
	        Spar(Ntr)%D13sw(n)	= -L11*E11*L13
	        Spar(Ntr)%D21sw(n)	= L12*E11*L11
	        Spar(Ntr)%D22sw(n)	= L12*E11*L12
	        Spar(Ntr)%D23sw(n)	= -L12*E11*L13
	        Spar(Ntr)%D31sw(n)	= -L13*E11*L11
	        Spar(Ntr)%D32sw(n)	= -L13*E11*L12
	        Spar(Ntr)%D33sw(n)	= L13*E11*L13
	        Spar(Ntr)%D44sw(n)	= L25*Spar(Ntr)%Gezsw(n)*L25
	        Spar(Ntr)%D45sw(n)	= L24*Spar(Ntr)%Gezsw(n)*L25
	        Spar(Ntr)%D54sw(n)	= L25*Spar(Ntr)%Gezsw(n)*L24
	        Spar(Ntr)%D55sw(n)	= L24*Spar(Ntr)%Gezsw(n)*L24
            ! Spar cap.
	        E11 = Spar(Ntr)%Eesc(n)/(1.0D0-Spar(Ntr)%miuezsc(n)*Spar(Ntr)%miuzesc(n))
	        Spar(Ntr)%D11sc(n)	= L11*E11*L11
	        Spar(Ntr)%D12sc(n)	= L11*E11*L12
	        Spar(Ntr)%D13sc(n)	= -L11*E11*L13
	        Spar(Ntr)%D21sc(n)	= L12*E11*L11
	        Spar(Ntr)%D22sc(n)	= L12*E11*L12
	        Spar(Ntr)%D23sc(n)	= -L12*E11*L13
	        Spar(Ntr)%D31sc(n)	= -L13*E11*L11
	        Spar(Ntr)%D32sc(n)	= -L13*E11*L12
	        Spar(Ntr)%D33sc(n)	= L13*E11*L13
        end do

    end subroutine SparConstitutiveMatrix
    !--------------------------------------------------------------------------------------------
    ! Routine to calculate spar geometry.
    elemental real(8) function SparGeometry(n,Ntr,etaj) result(zetasp_f)

      integer, intent(in) ::		n				! number of spar
      integer, intent(in) ::		Ntr				! planform number
      real(8), intent(in) ::		etaj            ! input eta-coordinate
      
      zetasp_f = -1.0D0+2.0D0*(Spar(Ntr)%ksp1(n)-Spar(Ntr)%ksp2(n)*etaj)

    end function SparGeometry

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate spar web geometry.
    subroutine SparWebGeometry(flag,n,zetai,etaj,zeta,x,y,z1sw,z2sw)

        use PlanformVariables, only : croot, ctip, Ntr
        use SkinVariablesRoutines, only: SkinGeometry2

        !Input variables
        integer, intent(in) ::		flag			! flag to choose calculation
								                    !   1 = calculate zeta
								                    !   2 = calculate limits of integration
        integer, intent(in) ::		n				! spar number
        real(8), intent(in) ::		zetai,etaj		! zeta,eta-coordinates of integration point
        real(8), intent(in) ::		x,y      		! point coordinates
        !Output variables
        real(8), intent(out) ::		zeta			! zeta-coordinate
        real(8), intent(out) ::		z1sw,z2sw		! limits of integration
        !Local variables
        real(8) ::				    z,slope 		! point coordinates
        real(8) ::				    z1,z2		    ! limits of integration
        real(8) ::				    tsw				! spar web thickness
        real(8) ::				    hsc				! spar cap height
        real(8) ::				    hsct			! spar cap height normal to skin		!**********************************
        real(8) ::				    c				! wing chord

        if(flag.EQ.1) then
            tsw		= Spar(Ntr)%tsw1(n)-Spar(Ntr)%tsw2(n)*etaj
            c		= 0.5D0*croot(Ntr)*(1.0D0-etaj)+0.5D0*ctip(Ntr)*(1.0D0+etaj)
            Spar(Ntr)%tcsw	= tsw/c
            zeta	= Spar(Ntr)%tcsw*zetai+SparGeometry(n,Ntr,etaj) !!!!!!!!!!!!!!!!!!Spar(Ntr)%zetasp
        end if
        if(flag.EQ.2) then
            hsct	= Spar(Ntr)%hsc1(n)-Spar(Ntr)%hsc2(n)*etaj					!***********************************   
            call SkinGeometry2(3,2,x,y,z,z1,z2,slope)
            hsc		= hsct*sqrt(1+slope*slope)				!************************************
            z1sw	= z2+hsc
            call SkinGeometry2(3,1,x,y,z,z1,z2,slope)
            hsc		= hsct*sqrt(1+slope*slope)				!************************************
            z2sw	= z1-hsc
        end if

    end subroutine SparWebGeometry

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate spar cap geometry.
    subroutine SparCapGeometry(flag,n,zetai,etaj,zeta,x,y,z1scU,z2scU,z1scL,z2scL)

        use PlanformVariables, only: croot, ctip, Ntr
        use SkinVariablesRoutines, only: SkinGeometry2

        !Input variables
        integer, intent(in) ::		flag			! flag to choose calculation
								                    !   1 = calculate zeta
								                    !   2 = calculate limits of integration
        integer, intent(in) ::		n				! spar number
        real(8), intent(in) ::		zetai,etaj		! zeta,eta-coordinates of integration point
        real(8), intent(in) ::		x,y      		! point coordinates
        !Output variables
        real(8), intent(out) ::		zeta			! zeta-coordinate
        real(8), intent(out) ::		z1scU,z2scU,z1scL,z2scL		! limits of integration
        !Local variables
        real(8) ::				    z,slope 		! point coordinates
        real(8) ::				    z1,z2		    ! limits of integration
        real(8) ::				    hsc				! spar cap height
        real(8) ::				    lsc             ! spar cap width
        real(8) ::				    hsct			! spar cap height normal to skin		!**********************************
        real(8) ::				    c				! wing chord

        if(flag.EQ.1) then
            lsc		= Spar(Ntr)%lsc1(n)-Spar(Ntr)%lsc2(n)*etaj
            c		= 0.5D0*croot(Ntr)*(1.0D0-etaj)+0.5D0*ctip(Ntr)*(1.0D0+etaj)
            Spar(Ntr)%lcsc	= lsc/c
            zeta	= Spar(Ntr)%lcsc*zetai+SparGeometry(n,Ntr,etaj) !!!!!!!!!!!!!!!!!!Spar(Ntr)%zetasp
        end if
        if(flag.EQ.2) then
            hsct	= Spar(Ntr)%hsc1(n)-Spar(Ntr)%hsc2(n)*etaj 
            call SkinGeometry2(3,1,x,y,z,z1,z2,slope)
            hsc		= hsct*sqrt(1+slope*slope)
            z1scU	= z1-hsc
            z2scU	= z1
            call SkinGeometry2(3,2,x,y,z,z1,z2,slope)
            hsc		= hsct*sqrt(1+slope*slope)
            z1scL	= z2
            z2scL	= z2+hsc
        end if

    end subroutine SparCapGeometry
    !--------------------------------------------------------------------------------------------
    ! Routine to apply scale to spar   
    subroutine ApplyScaleSpar(scale)
        use PlanformVariables, only: Ntr
        !Input vars
        real(8), intent(in) ::		scale				! current planform number
        !Local vars
        integer ::                  n
        
        forall(n=1:Spar(Ntr)%Nsp)
            Spar(Ntr)%tswroot(n) = Spar(Ntr)%tswroot(n) * scale 
            Spar(Ntr)%tswtip(n) = Spar(Ntr)%tswtip(n) * scale 
            
            Spar(Ntr)%hscroot(n) = Spar(Ntr)%hscroot(n) * scale 
            Spar(Ntr)%hsctip(n)  = Spar(Ntr)%hsctip(n) * scale  

            Spar(Ntr)%lscroot(n) = Spar(Ntr)%lscroot(n) * scale 
            Spar(Ntr)%lsctip(n) = Spar(Ntr)%lsctip(n) * scale    
            
            Spar(Ntr)%rhosw(n) = Spar(Ntr)%rhosw(n) / scale**3                     
            Spar(Ntr)%rhosc(n) = Spar(Ntr)%rhosc(n) / scale**3                         
            Spar(Ntr)%Gezsw(n) = Spar(Ntr)%Gezsw(n) / scale
            Spar(Ntr)%Eesw(n) = Spar(Ntr)%Eesw(n) / scale      

            Spar(Ntr)%Eesc(n) = Spar(Ntr)%Eesc(n) / scale                       
        end forall
        
    end subroutine ApplyScaleSpar

  
end module SparVariablesRoutines

!--------------------------------------------------------------------------------------------
! Module with rib variables.
module RibVariablesRoutines
  implicit none
  save
  
  TYPE :: RibVarsType
      integer ::            				Nrb					! number of ribs
      real(8), allocatable ::				krble(:)		    ! semi-span fraction of rib position at leading edge
      real(8), allocatable ::				krbte(:)		    ! semi-span fraction of rib position at trailing edge
      real(8), allocatable ::				trwle(:)		    ! rib web thickness at leading edge
      real(8), allocatable ::				trwte(:)		    ! rib web thickness at trailing edge
      real(8), allocatable ::				hrcle(:)		    ! rib cap height at leading edge
      real(8), allocatable ::				hrcte(:)		    ! rib cap height at trailing edge
      real(8), allocatable ::				lrcle(:)		    ! rib cap width at leading edge
      real(8), allocatable ::				lrcte(:)		    ! rib cap width at trailing edge
      real(8), allocatable ::				krb1(:)		        ! rib line equation coefficient
      real(8), allocatable ::				krb2(:)		        ! rib line equation coefficient
      real(8), allocatable ::				trw1(:)		        ! rib web thickness equation coefficient
      real(8), allocatable ::				trw2(:)		        ! rib web thickness equation coefficient
      real(8), allocatable ::				hrc1(:)		        ! rib cap height equation coefficient
      real(8), allocatable ::				hrc2(:)		        ! rib cap height equation coefficient
      real(8), allocatable ::				lrc1(:)		        ! rib cap width equation coefficient
      real(8), allocatable ::				lrc2(:)		        ! rib cap width equation coefficient
      real(8), allocatable ::				Lamdarb(:)		    ! rib sweep angle

      real(8), allocatable ::				rhorw(:)		    ! rib web material density
      real(8), allocatable ::				rhorc(:)		    ! rib cap material density

      real(8) ::            				tcrw				! rib web thickness to semi-span ratio
      real(8) ::            				lcrc				! rib cap width to semi-span ratio

      real(8), allocatable ::				Eerc(:)		        ! longitudinal elastic modulus in x
      real(8), allocatable ::				miuezrc(:)		    ! poisson ratio xy
      real(8), allocatable ::				miuzerc(:)		    ! poisson ratio yx
      real(8), allocatable ::				D11rc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D12rc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D13rc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D21rc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D22rc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D23rc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D31rc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D32rc(:)		    ! constitutive matrix element
      real(8), allocatable ::				D33rc(:)		    ! constitutive matrix element

      real(8), allocatable ::				Eerw(:)		        ! longitudinal elastic modulus in x
      real(8), allocatable ::				miuezrw(:)		    ! poisson ratio xy
      real(8), allocatable ::				miuzerw(:)		    ! poisson ratio yx
      real(8), allocatable ::				Gezrw(:)		    ! transverse elastic modulus
      real(8), allocatable ::				D11rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D12rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D13rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D21rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D22rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D23rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D31rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D32rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D33rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D44rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D45rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D54rw(:)		    ! constitutive matrix element
      real(8), allocatable ::				D55rw(:)		    ! constitutive matrix element
  END TYPE RibVarsType
  TYPE(RibVarsType), ALLOCATABLE :: Rib(:)

contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all rib variables
    subroutine AllocateRibVariables(Npl,Ntr,Nrb)
        !Input vars
        integer, intent(in) ::		Npl				! number of planforms
        integer, intent(in) ::		Ntr				! current planform number
        integer, intent(in) ::		Nrb				! number rib(s) for current planform
        if(.not. allocated(Rib)) then !First planform
            allocate(Rib(Npl))             
        end if
        Rib(Ntr)%Nrb = Nrb
        If ( Nrb > 0) then       
            allocate(Rib(Ntr)%krble(Nrb),Rib(Ntr)%krbte(Nrb),Rib(Ntr)%trwle(Nrb),Rib(Ntr)%trwte(Nrb),Rib(Ntr)%hrcle(Nrb),Rib(Ntr)%hrcte(Nrb),Rib(Ntr)%lrcle(Nrb),Rib(Ntr)%lrcte(Nrb),Rib(Ntr)%krb1(Nrb),Rib(Ntr)%krb2(Nrb),Rib(Ntr)%trw1(Nrb),Rib(Ntr)%trw2(Nrb))
            allocate(Rib(Ntr)%hrc1(Nrb),Rib(Ntr)%hrc2(Nrb),Rib(Ntr)%lrc1(Nrb),Rib(Ntr)%lrc2(Nrb),Rib(Ntr)%Lamdarb(Nrb),Rib(Ntr)%rhorw(Nrb),Rib(Ntr)%rhorc(Nrb))
            allocate(Rib(Ntr)%Eerc(Nrb),Rib(Ntr)%miuezrc(Nrb),Rib(Ntr)%miuzerc(Nrb),Rib(Ntr)%D11rc(Nrb),Rib(Ntr)%D12rc(Nrb),Rib(Ntr)%D13rc(Nrb),Rib(Ntr)%D21rc(Nrb),Rib(Ntr)%D22rc(Nrb))
            allocate(Rib(Ntr)%D23rc(Nrb),Rib(Ntr)%D31rc(Nrb),Rib(Ntr)%D32rc(Nrb),Rib(Ntr)%D33rc(Nrb),Rib(Ntr)%Eerw(Nrb),Rib(Ntr)%miuezrw(Nrb),Rib(Ntr)%miuzerw(Nrb),Rib(Ntr)%Gezrw(Nrb))
            allocate(Rib(Ntr)%D11rw(Nrb),Rib(Ntr)%D12rw(Nrb),Rib(Ntr)%D13rw(Nrb),Rib(Ntr)%D21rw(Nrb),Rib(Ntr)%D22rw(Nrb),Rib(Ntr)%D23rw(Nrb),Rib(Ntr)%D31rw(Nrb),Rib(Ntr)%D32rw(Nrb),Rib(Ntr)%D33rw(Nrb))
            allocate(Rib(Ntr)%D44rw(Nrb),Rib(Ntr)%D45rw(Nrb),Rib(Ntr)%D54rw(Nrb),Rib(Ntr)%D55rw(Nrb))
        end if

    end subroutine AllocateRibVariables    
    !--------------------------------------------------------------------------------------------
    ! Routine to deallocate all rib variables
    subroutine DeallocateRibVariables()
    
        deallocate(Rib)

    end subroutine DeallocateRibVariables
    !--------------------------------------------------------------------------------------------
    ! Routine to calculate rib constitutive matrix.
    subroutine RibConstitutiveMatrix(Ntr) 
        !Input Vars
        integer, intent(in) ::  Ntr		            ! number of the wing trapezoid beeing computed
        !Local Vars
        real(8) ::				L11,L12,L13,L24,L25	! direction cosines
        real(8) ::				E11             	! stiffness element
        integer ::				n					! iteration variables

        do n=1,Rib(Ntr)%Nrb
            L11			= dcos(Rib(Ntr)%Lamdarb(n))*dcos(Rib(Ntr)%Lamdarb(n))
            L12			= dsin(Rib(Ntr)%Lamdarb(n))*dsin(Rib(Ntr)%Lamdarb(n))
            L13			= dsin(Rib(Ntr)%Lamdarb(n))*dcos(Rib(Ntr)%Lamdarb(n))
            L24			= dcos(Rib(Ntr)%Lamdarb(n))				
            L25			= dsin(Rib(Ntr)%Lamdarb(n))	
					
            ! Rib web.
            E11			= Rib(Ntr)%Eerw(n)/(1.0D0-Rib(Ntr)%miuezrw(n)*Rib(Ntr)%miuzerw(n))
            Rib(Ntr)%D11rw(n)	= L11*E11*L11
            Rib(Ntr)%D12rw(n)	= L11*E11*L12
            Rib(Ntr)%D13rw(n)	= L11*E11*L13
            Rib(Ntr)%D21rw(n)	= L12*E11*L11
            Rib(Ntr)%D22rw(n)	= L12*E11*L12
            Rib(Ntr)%D23rw(n)	= L12*E11*L13
            Rib(Ntr)%D31rw(n)	= L13*E11*L11
            Rib(Ntr)%D32rw(n)	= L13*E11*L12
            Rib(Ntr)%D33rw(n)	= L13*E11*L13
            Rib(Ntr)%D44rw(n)	= L25*Rib(Ntr)%Gezrw(n)*L25
            Rib(Ntr)%D45rw(n)	= -L24*Rib(Ntr)%Gezrw(n)*L25
            Rib(Ntr)%D54rw(n)	= -L25*Rib(Ntr)%Gezrw(n)*L24
            Rib(Ntr)%D55rw(n)	= L24*Rib(Ntr)%Gezrw(n)*L24

            ! Rib cap.
            E11			= Rib(Ntr)%Eerc(n)/(1.0D0-Rib(Ntr)%miuezrc(n)*Rib(Ntr)%miuzerc(n))
            Rib(Ntr)%D11rc(n)	= L11*E11*L11
            Rib(Ntr)%D12rc(n)	= L11*E11*L12
            Rib(Ntr)%D13rc(n)	= L11*E11*L13
            Rib(Ntr)%D21rc(n)	= L12*E11*L11
            Rib(Ntr)%D22rc(n)	= L12*E11*L12
            Rib(Ntr)%D23rc(n)	= L12*E11*L13
            Rib(Ntr)%D31rc(n)	= L13*E11*L11
            Rib(Ntr)%D32rc(n)	= L13*E11*L12
            Rib(Ntr)%D33rc(n)	= L13*E11*L13
        end do
            
    end subroutine RibConstitutiveMatrix
    !--------------------------------------------------------------------------------------------
    ! Routine to calculate rib geometry.    
    elemental real(8) function RibGeometry(n,Ntr,zetai) result(etarb_f)
        !Input vars
        integer, intent(in) ::			n				! number of rib
        integer, intent(in) ::			Ntr				! trapezoid number
        real(8), intent(in) ::			zetai			! input zeta-coordinate

        etarb_f = -1.0D0+2.0D0*(Rib(Ntr)%krb1(n)-Rib(Ntr)%krb2(n)*zetai)

    end function RibGeometry

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate rib web geometry.
    subroutine RibWebGeometry(flag,n,zetai,etaj,eta,x,y,z1rw,z2rw,error)

        use SkinVariablesRoutines, only:SkinGeometry2
        use PlanformVariables, only: Ntr,b2,Section,yroot

        !Input Variables
        integer, intent(in) ::                      flag		! flag to choose calculation
									                            !  1 = calculate zeta
									                            !  2 = calculate limits of integration
        integer, intent(in) ::                      n			! rib number
        real(8), intent(in) ::                      zetai,etaj  ! zeta,eta-coordinates of integration point
        real(8), intent(in) ::                      x,y         ! point coordinates
        !Output Variables
        real(8), intent(out) ::                     z1rw,z2rw   ! limits of integration
        logical, optional, intent(out) ::           error
        !In/output variables
        real(8), intent(inout) ::                   eta			! eta-coordinate
        !Local Variables
        real(8) ::                                  z,slope				! point coordinates, skin slope
        real(8) ::                                  z1,z2		        ! limits of integration
        real(8) ::                                  trw					! spar web thickness
        real(8) ::                                  hrc					! spar cap height
        real(8) ::                                  hrct				! rib cap height normal to skin		!**********************************
        real(8) ::                                  z_ref
        real(8) ::                                  y_find
        real(8) ::                                  eps=epsilon(z_ref)*10.0D0
        integer ::                                  i,iref

        If(Present(error)) error=.False.
        if(flag.EQ.1) then
            trw		= Rib(Ntr)%trw1(n)-Rib(Ntr)%trw2(n)*zetai
            Rib(Ntr)%tcrw	= trw/b2(Ntr)
            eta		= Rib(Ntr)%tcrw*etaj+RibGeometry(n,Ntr,zetai) !!!!!!!!!!!!!!!!!!! Rib(Ntr)%etarb
            if(eta > 1.0D0) then
                eta=1.0D0
            elseif(eta < -1.0D0) then
                eta=-1.0D0
            end if
        end if
        if(flag.EQ.2) then
            y_find = y - yroot(Ntr) 
            iref = -1
            do i=2,Section(Ntr)%n_sections
                IF(y_find - Section(Ntr)%y_sections(i-1) >= -eps  .AND. y_find - Section(Ntr)%y_sections(i) <= eps ) Then !sum with eps to avoid numerical error when evaluating logical <=
                    iref=i
                    exit
                end if
            end do
            if(iref == -1) then
                iref= Section(Ntr)%n_sections
            end if
            z_ref=(Section(Ntr)%z_sections(iref)-Section(Ntr)%z_sections(iref-1))/(Section(Ntr)%y_sections(iref)-Section(Ntr)%y_sections(iref-1))*(y - yroot(Ntr) - Section(Ntr)%y_sections(iref-1))

            hrct	= Rib(Ntr)%hrc1(n)-Rib(Ntr)%hrc2(n)*zetai				!*********************************
            call SkinGeometry2(3,2,x,y,z,z1,z2,slope)
            hrc = hrct*dsqrt(1.0D0+slope*slope)			!**********************************
            z1rw	= z2+hrc
            if (z1rw >= z_ref) then
                z1rw = z_ref
                If(Present(error)) error=.True.
            end if
            call SkinGeometry2(3,1,x,y,z,z1,z2,slope) !top 
            hrc	= hrct*dsqrt(1.0D0+slope*slope)			!**********************************
            z2rw = z1-hrc
            if (z2rw <= z_ref) then
	            z2rw = z_ref
	            If(Present(error)) error=.True.
            end if  
        end if

    end subroutine RibWebGeometry

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate rib cap geometry.
    subroutine RibCapGeometry(flag,n,zetai,etaj,eta,x,y,z1rcU,z2rcU,z1rcL,z2rcL,error)

        use PlanformVariables, only: ctip,croot,Ntr,b2,Section,yroot
        use SkinVariablesRoutines, only: SkinGeometry2

        !Input Variables
        integer, intent(in) ::              flag		! flag to choose calculation:
                                                        !  1 = calculate zeta
                                                        !  2 = calculate limits of integration
        integer, intent(in) ::              n			! rib number
        real(8), intent(in) ::              zetai,etaj  ! zeta,eta-coordinates of integration point
        real(8), intent(in) ::              x,y         ! point coordinates
        !In/output variables
        real(8), intent(inout) ::           eta			! eta-coordinate    
        !Output Variables
        real(8), intent(out) ::             z1rcU,z2rcU,z1rcL,z2rcL ! limits of integration   
        logical, optional, intent(out) ::   error
        !Local vars
        real(8) ::				            z,slope				    ! point coordinates
        real(8) ::				            z1,z2				    ! limits of integration
        real(8) ::				            lrc					    ! rib cap width
        real(8) ::				            hrc					    ! rib cap height
        real(8) ::				            hrct				    ! rib cap height normal to skin		!**********************************
        real(8) ::                          z_ref
        real(8) ::                          eps=epsilon(z_ref)*10.0D0
        integer ::                          i,iref
        
        If(Present(error)) error=.False.
        if(flag.EQ.1) then
            lrc		= Rib(Ntr)%lrc1(n)-Rib(Ntr)%lrc2(n)*zetai
            Rib(Ntr)%lcrc	= lrc/b2(Ntr)
            eta		= Rib(Ntr)%lcrc*etaj+RibGeometry(n,Ntr,zetai) !!!!!!!!!!!!!!!!!!! Rib(Ntr)%etarb
            if(eta > 1.0D0) then
                eta=1.0D0
            elseif(eta < -1.0D0) then
                eta=-1.0D0
            end if
        end if
        if(flag.EQ.2) then
            iref = -1
            do i=2,Section(Ntr)%n_sections
                IF(y - yroot(Ntr) >= Section(Ntr)%y_sections(i-1)-eps .AND. y - yroot(Ntr) <= Section(Ntr)%y_sections(i)+eps) Then !sum with eps to avoid numerical error when evaluating logical <=
                    iref=i
                    exit
                end if
            end do
            if(iref == -1) then
                iref= Section(Ntr)%n_sections
            end if
            z_ref=(Section(Ntr)%z_sections(iref)-Section(Ntr)%z_sections(iref-1))/(Section(Ntr)%y_sections(iref)-Section(Ntr)%y_sections(iref-1))*(y - yroot(Ntr) - Section(Ntr)%y_sections(iref-1))

            hrct	= Rib(Ntr)%hrc1(n)-Rib(Ntr)%hrc2(n)*zetai								!*********************************
            call SkinGeometry2(3,1,x,y,z,z1,z2,slope)
            hrc		= hrct*dsqrt(1+slope*slope)							!**********************************
            z1rcU	= z1-hrc
            if (z1rcU <= z_ref) then
                z1rcU = z_ref
                If(Present(error)) error=.True.
            end if
            !if(z1rcU < 0) z1rcU=0.0D0
            z2rcU	= z1
            call SkinGeometry2(3,2,x,y,z,z1,z2,slope)
            hrc		= hrct*dsqrt(1.0D0+slope*slope)							!**********************************
            z1rcL	= z2
            z2rcL	= z2+hrc
            if (z2rcL >= z_ref) then
	            z2rcL = z_ref
	            If(Present(error)) error=.True.
            end if  
            !if(z2rcL > 0) z1rcU=0.0D0
        end if

    end subroutine RibCapGeometry
    !--------------------------------------------------------------------------------------------
    ! Routine to apply scale to spar   
    subroutine ApplyScaleRib(scale)
        use PlanformVariables, only: Ntr
        !Input vars
        real(8), intent(in) ::		scale				! current planform number
        !Local vars
        integer ::                  n
        
        forall(n=1:Rib(Ntr)%Nrb)
            Rib(Ntr)%trwle(n) = Rib(Ntr)%trwle(n) * scale
            Rib(Ntr)%trwte(n) = Rib(Ntr)%trwte(n) * scale       ! rib web thickness at leading and trailing edges
            Rib(Ntr)%hrcle(n) = Rib(Ntr)%hrcle(n) * scale
            Rib(Ntr)%hrcte(n) = Rib(Ntr)%hrcte(n) * scale       ! rib cap height at leading and trailing edges
            Rib(Ntr)%lrcle(n) = Rib(Ntr)%lrcle(n) * scale
            Rib(Ntr)%lrcte(n) = Rib(Ntr)%lrcte(n) * scale       ! rib cap width at leading and trailing edges
            Rib(Ntr)%rhorw(n) = Rib(Ntr)%rhorw(n) / scale**3                       ! material density of rib web
            Rib(Ntr)%rhorc(n) = Rib(Ntr)%rhorc(n) / scale**3                       ! material density of rib cap
            Rib(Ntr)%Gezrw(n) = Rib(Ntr)%Gezrw(n) / scale
            Rib(Ntr)%Eerw(n) = Rib(Ntr)%Eerw(n) / scale          ! transversal and logitudinal elastic moduli of rib web
            Rib(Ntr)%Eerc(n) = Rib(Ntr)%Eerc(n) / scale                        ! logitudinal elastic modulus of rib cap
        end forall
        
    end subroutine ApplyScaleRib

end module RibVariablesRoutines

!--------------------------------------------------------------------------------------------
! Module with stringer variables.
module StringerVariablesRoutines
    
  implicit none
  save
  TYPE :: StringerVarsType
      integer ::            				Nst					! number of stringers
      real(8), allocatable ::				kstroot(:)		    ! chord fraction of stringer position at root
      real(8), allocatable ::				ksttip(:)		    ! chord fraction of stringer position at tip
      real(8), allocatable ::				hstroot(:)		    ! stringer height at root
      real(8), allocatable ::				hsttip(:)		    ! stringer height at tip
      real(8), allocatable ::				lstroot(:)		    ! stringer width at root
      real(8), allocatable ::				lsttip(:)		    ! stringer width at tip
      real(8), allocatable ::				kst1(:)		        ! stringer line equation coefficient
      real(8), allocatable ::				kst2(:)		        ! stringer line equation coefficient
      real(8), allocatable ::				hst1(:)		        ! stringer height equation coefficient
      real(8), allocatable ::				hst2(:)		        ! stringer height equation coefficient
      real(8), allocatable ::				lst1(:)		        ! stringer width equation coefficient
      real(8), allocatable ::				lst2(:)		        ! stringer width equation coefficient
      real(8), allocatable ::				Lamdast(:)		    ! stringer sweep angle

      real(8), allocatable ::				rhost(:)		    ! stringer material density

      real(8) ::               				lcst				! stringer width to wing chord ratio

      real(8), allocatable ::				Eest(:)		        ! longitudinal elastic modulus in x
      real(8), allocatable ::				miuezst(:)		    ! poisson ratio xy
      real(8), allocatable ::				miuzest(:)		    ! poisson ratio yx
      real(8), allocatable ::				D11st(:)		    ! constitutive matrix element
      real(8), allocatable ::				D12st(:)		    ! constitutive matrix element
      real(8), allocatable ::				D13st(:)		    ! constitutive matrix element
      real(8), allocatable ::				D21st(:)		    ! constitutive matrix element
      real(8), allocatable ::				D22st(:)		    ! constitutive matrix element
      real(8), allocatable ::				D23st(:)		    ! constitutive matrix element
      real(8), allocatable ::				D31st(:)		    ! constitutive matrix element
      real(8), allocatable ::				D32st(:)		    ! constitutive matrix element
      real(8), allocatable ::				D33st(:)		    ! constitutive matrix element
  END TYPE StringerVarsType
  TYPE(StringerVarsType), ALLOCATABLE :: Stringer(:)
contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all stringer variables
    subroutine AllocateStringerVariables(Npl,Ntr,Nst)
        !Input vars
        integer, intent(in) ::		Npl				! number of planforms
        integer, intent(in) ::		Ntr				! current planform number
        integer, intent(in) ::		Nst				! number stringer(s) for current planform
        if( .not. allocated(Stringer)) then !First planform
            allocate(Stringer(Npl))             
        end if
        Stringer(Ntr)%Nst = Nst
        if( Nst > 0 ) then
            allocate(Stringer(Ntr)%kstroot(Nst),Stringer(Ntr)%ksttip(Nst),Stringer(Ntr)%hstroot(Nst),Stringer(Ntr)%hsttip(Nst),Stringer(Ntr)%lstroot(Nst),Stringer(Ntr)%lsttip(Nst),Stringer(Ntr)%kst1(Nst))
            allocate(Stringer(Ntr)%kst2(Nst),Stringer(Ntr)%hst1(Nst),Stringer(Ntr)%hst2(Nst),Stringer(Ntr)%lst1(Nst),Stringer(Ntr)%lst2(Nst),Stringer(Ntr)%Lamdast(Nst),Stringer(Ntr)%rhost(Nst))
            allocate(Stringer(Ntr)%Eest(Nst),Stringer(Ntr)%miuezst(Nst),Stringer(Ntr)%miuzest(Nst),Stringer(Ntr)%D11st(Nst),Stringer(Ntr)%D12st(Nst),Stringer(Ntr)%D13st(Nst),Stringer(Ntr)%D21st(Nst),Stringer(Ntr)%D22st(Nst))
            allocate(Stringer(Ntr)%D23st(Nst),Stringer(Ntr)%D31st(Nst),Stringer(Ntr)%D32st(Nst),Stringer(Ntr)%D33st(Nst))
        end if       
       
    end subroutine AllocateStringerVariables    
    !--------------------------------------------------------------------------------------------
    ! Routine to deallocate all stringer variables
    subroutine DeallocateStringerVariables ()
        
        deallocate(Stringer)        
        !deallocate(Stringer(Ntr)%kstroot,Stringer(Ntr)%ksttip,Stringer(Ntr)%hstroot,Stringer(Ntr)%hsttip,Stringer(Ntr)%lstroot,Stringer(Ntr)%lsttip,Stringer(Ntr)%kst1)
        !deallocate(Stringer(Ntr)%kst2,Stringer(Ntr)%hst1,Stringer(Ntr)%hst2,Stringer(Ntr)%lst1,Stringer(Ntr)%lst2,Stringer(Ntr)%Lamdast,Stringer(Ntr)%rhost)
        !deallocate(Stringer(Ntr)%Eest,Stringer(Ntr)%miuezst,Stringer(Ntr)%miuzest,Stringer(Ntr)%D11st,Stringer(Ntr)%D12st,Stringer(Ntr)%D13st,Stringer(Ntr)%D21st,Stringer(Ntr)%D22st)
        !deallocate(Stringer(Ntr)%D23st,Stringer(Ntr)%D31st,Stringer(Ntr)%D32st,Stringer(Ntr)%D33st)
    
    end subroutine DeallocateStringerVariables 
    !--------------------------------------------------------------------------------------------
    ! Routine to calculate stringer constitutive matrix.
    subroutine StringerConstitutiveMatrix(Ntr)
        
        !Input variables
        integer, intent(in) ::		Ntr				! current planform number
        !Local variables
        real(8) ::				    L11,L12,L13			! direction cosines
        real(8) ::				    E11					! stiffness element
        integer ::				    n					! iteration variables

        do n=1,Stringer(Ntr)%Nst
            L11			= dsin(Stringer(Ntr)%Lamdast(n))*dsin(Stringer(Ntr)%Lamdast(n))
            L12			= dcos(Stringer(Ntr)%Lamdast(n))*dcos(Stringer(Ntr)%Lamdast(n))
            L13			= dsin(Stringer(Ntr)%Lamdast(n))*dcos(Stringer(Ntr)%Lamdast(n))
            E11			= Stringer(Ntr)%Eest(n)/(1.0D0-Stringer(Ntr)%miuezst(n)*Stringer(Ntr)%miuzest(n))
            Stringer(Ntr)%D11st(n)	= L11*E11*L11
            Stringer(Ntr)%D12st(n)	= L11*E11*L12
            Stringer(Ntr)%D13st(n)	= -L11*E11*L13
            Stringer(Ntr)%D21st(n)	= L12*E11*L11
            Stringer(Ntr)%D22st(n)	= L12*E11*L12
            Stringer(Ntr)%D23st(n)	= -L12*E11*L13
            Stringer(Ntr)%D31st(n)	= -L13*E11*L11
            Stringer(Ntr)%D32st(n)	= -L13*E11*L12
            Stringer(Ntr)%D33st(n)	= L13*E11*L13
        end do

    end subroutine StringerConstitutiveMatrix
    !--------------------------------------------------------------------------------------------
    ! Routine to calculate stringer geometry.
    elemental real(8) function StringerGeometry1(n,Ntr,etaj) result(zetast_f)
        !Input variables
        integer, intent(in) ::		n				! number of stringer
        integer, intent(in) ::		Ntr				! current planform
        real(8), intent(in) ::		etaj            ! input eta-coordinate

        zetast_f = -1.0D0+2.0D0*(Stringer(Ntr)%kst1(n)-Stringer(Ntr)%kst2(n)*etaj)

    end function StringerGeometry1
    !--------------------------------------------------------------------------------------------
    ! Routine to calculate stringer geometry.
    subroutine StringerGeometry2(flag,n,zetai,etaj,zeta,x,y,z1stU,z2stU,z1stL,z2stL)

        use PlanformVariables, only: croot, ctip,ntr
        use SkinVariablesRoutines, only: SkinGeometry2
       
        !Input Variables
        integer, intent(in) ::              flag		! flag to choose calculation:
                                                        !  1 = calculate zeta
                                                        !  2 = calculate limits of integration
        integer, intent(in) ::              n			! stringer number
        real(8), intent(in) ::              zetai,etaj  ! zeta,eta-coordinates of integration point
        real(8), intent(in) ::              x,y         ! point coordinates
        !In/output variables   
        real(8), intent(inout) ::           zeta        ! point 
        !Output Variables
        real(8), intent(out) ::             z1stU,z2stU,z1stL,z2stL ! limits of integration   
        !Local vars
        real(8) ::				            z,slope				    ! point coordinates
        real(8) ::				            z1,z2				    ! limits of integration
        real(8) ::				            lst					    ! stringer width
        real(8) ::				            hst					    ! stringer height
        real(8) ::				            hstt				    ! stringer height normal to skin
        real(8) ::				            c					    ! wing chord

        if(flag.EQ.1) then
            lst		= Stringer(Ntr)%lst1(n)-Stringer(Ntr)%lst2(n)*etaj
            c		= 0.5D0*croot(Ntr)*(1.0D0-etaj)+0.5D0*ctip(Ntr)*(1.0D0+etaj)
            Stringer(Ntr)%lcst	= lst/c
            zeta	= Stringer(Ntr)%lcst*zetai+StringerGeometry1(n,Ntr,etaj) !!!!!!!!!Stringer(Ntr)%zetast
        end if
        if(flag.EQ.2) then
            hstt	= Stringer(Ntr)%hst1(n)-Stringer(Ntr)%hst2(n)*etaj
            call SkinGeometry2(3,1,x,y,z,z1,z2,slope)
            hst		= hstt*(1.0D0+slope*slope) 
            z1stU	= z1-hst
            z2stU	= z1
            call SkinGeometry2(3,2,x,y,z,z1,z2,slope)
            hst		= hstt*(1.0D0+slope*slope) 
            z1stL	= z2
            z2stL	= z2+hst
        end if       

    end subroutine StringerGeometry2
    !--------------------------------------------------------------------------------------------
    ! Routine to apply scale to stringer  
    subroutine ApplyScaleStringer(scale)
        use PlanformVariables, only: Ntr
        !Input vars
        real(8), intent(in) ::		scale				! current planform number
        !Local vars
        integer ::                  n
        
        forall(n=1:Stringer(Ntr)%Nst)
            Stringer(Ntr)%hstroot(n) = Stringer(Ntr)%hstroot(n) * scale
            Stringer(Ntr)%hsttip(n) = Stringer(Ntr)%hsttip(n) * scale    ! stringer height at root and tip
            Stringer(Ntr)%lstroot(n) = Stringer(Ntr)%lstroot(n) * scale
            Stringer(Ntr)%lsttip(n) = Stringer(Ntr)%lsttip(n) * scale    ! stringer width at root and tip
            Stringer(Ntr)%rhost(n) = Stringer(Ntr)%rhost(n) / scale**3                ! material density
            Stringer(Ntr)%Eest(n) = Stringer(Ntr)%Eest(n) / scale
        end forall
        
    end subroutine ApplyScaleStringer
      
end module StringerVariablesRoutines

! ----------------------------------------------------------------------------------------------
! Module for spring variables.
module SpringVariables
    implicit none
    save
    TYPE SpringVarsType
        integer ::                      Nspr				! number of springs
        real(8), allocatable ::			ku(:)			    ! linear spring stiffness in x direction (Nspr)
        real(8), allocatable ::			kv(:)			    ! linear stiffness in y direction (Nspr)
        real(8), allocatable ::			kw(:)			    ! linear spring stiffness in z direction (Nspr)
        real(8), allocatable ::			krx(:)		        ! rotation spring stiffness in x direction (Nspr)
        real(8), allocatable ::			kry(:)		        ! rotation spring stiffness in y direction (Nspr)
        real(8), allocatable ::			xspr(:)		        ! x coordinate of spring (Nspr)
        real(8), allocatable ::			yspr(:)		        ! y coordinate of spring (Nspr)
    END TYPE SpringVarsType
    TYPE(SpringVarsType), ALLOCATABLE :: Spring(:)
contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all Spring variables
    subroutine AllocateSpringVariables(Npl,Ntr,Nspr)
        !Input vars
        integer, intent(in) ::		Npl				! number of planforms
        integer, intent(in) ::		Ntr				! current planform number
        integer, intent(in) ::		Nspr				! number stringer(s) for current planform
        if( .not. allocated(Spring)) then !First planform
            allocate(Spring(Npl))             
        end if
        Spring(Ntr)%Nspr = Nspr
        if( Nspr > 0 ) then
            allocate(Spring(Ntr)%ku(Nspr),Spring(Ntr)%kv(Nspr),Spring(Ntr)%kw(Nspr),Spring(Ntr)%krx(Nspr),Spring(Ntr)%kry(Nspr),Spring(Ntr)%xspr(Nspr),Spring(Ntr)%yspr(Nspr))
        end if       
       
    end subroutine AllocateSpringVariables    
    !--------------------------------------------------------------------------------------------
    ! Routine to deallocate all Spring variables
    subroutine DeallocateSpringVariables ()
        
        deallocate(Spring)    
        !deallocate(Spring(Ntr)%ku,kv,kw,krx,kry,xspr,yspr)
    
    end subroutine DeallocateSpringVariables 

    
end module SpringVariables

!--------------------------------------------------------------------------------------------
! Module with coefficients of Legendre polynomials.
module PolynomialCoefficients
    implicit none
    save
    
    logical ::                  polynomial              ! .TRUE.=Legendre, .FALSE.=Ritz
    integer ::                  k                       ! number of polynomials 
    real(8), allocatable ::     Lcoef(:,:)              ! coefficients' matrix for the Legendre Polynomials
    real(8), allocatable ::     weights(:),nodes(:)     ! gauss-legendre quadrature weights and sampling points
    !Ritz Polynomial variables
    integer ::				    nx(15)
    integer ::				    ny(15)

    data nx/0,1,0,2,1,0,3,2,1,0,4,3,2,1,0/
    data ny/1,1,2,1,2,3,1,2,3,4,1,2,3,4,5/
contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all polynomial variables
    subroutine AllocatePolynomialCoefficients(ik)
        !Input vars
        integer, intent(in) ::		ik				! number of planforms

        allocate(weights(ik),nodes(ik),Lcoef(k,k))
       
    end subroutine AllocatePolynomialCoefficients  
    !--------------------------------------------------------------------------------------------
    ! Routine to deallocate all polynomial variables
    subroutine DeallocateSpringVariables ()
        
        deallocate(weights,nodes,Lcoef)    
    
    end subroutine DeallocateSpringVariables     
end module PolynomialCoefficients

!--------------------------------------------------------------------------------------------
! Module with deformation variables.
module DeformationVariables
    implicit none
    save
    
    TYPE DeformationVarsType
        real(8), allocatable  ::		MM(:,:)         ! mass matrix
        real(8), allocatable  ::		KM(:,:)	        ! stiffness matrix
        real(8), allocatable  ::		PV(:)           ! load vector
        real(8), allocatable  ::		q(:)			! deformation coefficients vector
    END TYPE DeformationVarsType
    TYPE(DeformationVarsType), ALLOCATABLE :: Deformation(:) 
    
    real(8), allocatable  ::		MM_Global(:,:)	    ! complete wing mass matrix
    real(8), allocatable  ::		KM_Global(:,:)	    ! complete wing stiffness matrix
    real(8), allocatable  ::		PV_Global(:)        ! complete wing load vector
    real(8), allocatable  ::		q_Global(:)        ! complete wing deformation coefficients vector
contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all deformation variables
    subroutine AllocateLocalDeformationVariables(Npl,k)

        !Input vars
        integer, intent(in) ::		Npl				! number of planforms
        integer ::                  k               ! number of polynomials 
        !local vars
        integer ::                  i               ! iteration var 
        
        if( .not. allocated(Deformation)) then !First planform
            allocate(Deformation(Npl))             
        end if
        
        do i = 1, Npl
            allocate(Deformation(i)%MM(5*k*k,5*k*k),Deformation(i)%KM(5*k*k,5*k*k),Deformation(i)%PV(5*k*k),Deformation(i)%q(5*k*k))
        end do
        allocate(KM_Global(Npl*5*k*k,Npl*5*k*k),MM_Global(Npl*5*k*k,Npl*5*k*k),PV_Global(Npl*5*k*k),q_Global(Npl*5*k*k))
        
    end subroutine AllocateLocalDeformationVariables
    !--------------------------------------------------------------------------------------------
    ! Routine to deallocate all deformation variables    
    subroutine DeallocateLocalDeformationVariables()
        
        deallocate(Deformation)
    
    end subroutine DeallocateLocalDeformationVariables
    
end module DeformationVariables


! ----------------------------------------------------------------------------------------------
! Module for frequencies.
module FrequenciesVariables

    save

    real(8), allocatable  ::		freq0(:)	    ! complete frequency matrix
    
contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all deformation variables
    subroutine AllocateFrequenciesVariables(k)

        !Input vars
        integer, intent(in) ::		k				! number of planforms
        
        allocate(freq0(5*k*k))
       
    end subroutine AllocateFrequenciesVariables
    !--------------------------------------------------------------------------------------------
    ! Routine to deallocate all deformation variables    
    subroutine DeallocateFrequenciesVariables()
        
        deallocate(freq0)
    
    end subroutine DeallocateFrequenciesVariables
    
end module FrequenciesVariables

! ----------------------------------------------------------------------------------------------
! Module for concentrated load variables.
module ConcentratedLoadVariables
    implicit none
    save

    TYPE :: CLoadVarsType
        integer ::                      Ncl				    ! number of loads
        real(8), allocatable ::	        Pclx(:)				! load in x direction
        real(8), allocatable ::	        Pcly(:)				! load in y direction
        real(8), allocatable ::	        Pclz(:)				! load in z direction
        real(8), allocatable ::	        xcl(:)				! load x-coordinate
        real(8), allocatable ::	        ycl(:)				! load y-coordinate
        real(8), allocatable ::	        zcl(:)				! load z-coordinate
        real(8), allocatable ::	        xcl0(:)				! load x-coordinate
        real(8), allocatable ::	        ycl0(:)				! load y-coordinate
        real(8), allocatable ::	        zcl0(:)				! load z-coordinate
    END TYPE CLoadVarsType
    TYPE(CLoadVarsType), ALLOCATABLE :: CLoad(:)

contains

    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all Concentrated Load variables
    subroutine AllocateConcentratedLoadVariables(Npl,Ntr,Ncl)
        !Input vars
        integer, intent(in) ::		Npl				! number of planforms
        integer, intent(in) ::		Ntr				! current planform number
        integer, intent(in) ::		Ncl				! number Concentrated Load for current planform
        if( .not. allocated(CLoad)) then !First planform
            allocate(CLoad(Npl))             
        end if
        CLoad(Ntr)%Ncl = Ncl
        if( Ncl > 0 ) then
            allocate(CLoad(Ntr)%Pclx(Ncl),CLoad(Ntr)%Pcly(Ncl),CLoad(Ntr)%Pclz(Ncl),CLoad(Ntr)%xcl(Ncl),CLoad(Ntr)%ycl(Ncl),CLoad(Ntr)%zcl(Ncl),CLoad(Ntr)%xcl0(Ncl),CLoad(Ntr)%ycl0(Ncl),CLoad(Ntr)%zcl0(Ncl))
        end if
    end subroutine AllocateConcentratedLoadVariables 
    !--------------------------------------------------------------------------------------------
    ! Routine to deallocate all Concentrated Load variables    
    subroutine DeallocateConcentratedLoadVariables()
        
        deallocate(CLoad)
    
    end subroutine DeallocateConcentratedLoadVariables
    
end module ConcentratedLoadVariables

! ----------------------------------------------------------------------------------------------
! Module for actuator load variables.
module ActuatorLoadVariables
    implicit none
    save
    TYPE ALoadVarsType
        integer ::            				Nac					! number of actuators
        real(8), allocatable ::				Pac(:)				! load
        real(8), allocatable ::				kxac1(:)			! chord fraction of end load position
        real(8), allocatable ::				kyac1(:)			! span fraction	of end load position
        real(8), allocatable ::				xac1(:)			    ! load x-coordinate
        real(8), allocatable ::				yac1(:)			    ! load y-coordinate
        real(8), allocatable ::				zac1(:)			    ! load z-coordinate
        real(8), allocatable ::				kxac2(:)			! chord fraction of end load position
        real(8), allocatable ::				kyac2(:)			! span fraction	of end load position
        real(8), allocatable ::				xac2(:)			    ! load x-coordinate
        real(8), allocatable ::				yac2(:)			    ! load y-coordinate
        real(8), allocatable ::				zac2(:)			    ! load z-coordinate
        real(8), allocatable ::				Dxac(:)			    ! end point x-coordinate difference
        real(8), allocatable ::				Dyac(:)			    ! end point y-coordinate difference
        real(8), allocatable ::				Dzac(:)			    ! end point z-coordinate difference
        real(8), allocatable ::				lac(:)				! actuator length

        integer ::                          Nal					! number of loads
        real(8), allocatable ::	            Palx(:)			    ! load in x direction
        real(8), allocatable ::	            Paly(:)			    ! load in y direction
        real(8), allocatable ::	            Palz(:)			    ! load in z direction
        real(8), allocatable ::	            xal(:)				! load x-coordinate
        real(8), allocatable ::	            yal(:)				! load y-coordinate
        real(8), allocatable ::	            zal(:)				! load z-coordinate
        real(8), allocatable ::	            xal0(:)			    ! load x-coordinate
        real(8), allocatable ::	            yal0(:)			    ! load y-coordinate
        real(8), allocatable ::	            zal0(:)			    ! load z-coordinate
    END TYPE ALoadVarsType
    TYPE(ALoadVarsType), ALLOCATABLE :: ALoad(:)
  
  contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all Actuator Load variables
    subroutine AllocateActuatorLoadVariables(Npl,Ntr,Nac,Nal)
        !Input vars
        integer, intent(in) ::		Npl				! number of planforms
        integer, intent(in) ::		Ntr				! current planform number
        integer, intent(in) ::		Nac			! number Concentrated Load for current planform
        integer, intent(in) ::		Nal			! number Concentrated Load for current planform
        if( .not. allocated(ALoad)) then !First planform
            allocate(ALoad(Npl))             
        end if
        ALoad(Ntr)%Nac = Nac
        ALoad(Ntr)%Nal = Nal
        if( Nac > 0 ) then
            allocate(ALoad(Ntr)%Pac(Nac),ALoad(Ntr)%kxac1(Nac),ALoad(Ntr)%kyac1(Nac),ALoad(Ntr)%xac1(Nac),ALoad(Ntr)%yac1(Nac),ALoad(Ntr)%kxac2(Nac),ALoad(Ntr)%kyac2(Nac),ALoad(Ntr)%xac2(Nac),ALoad(Ntr)%yac2(Nac),ALoad(Ntr)%zac2(Nac),ALoad(Ntr)%Dxac(Nac),ALoad(Ntr)%Dyac(Nac),ALoad(Ntr)%Dzac(Nac),ALoad(Ntr)%lac(Nac))
            allocate(ALoad(Ntr)%Palx(Nal),ALoad(Ntr)%Paly(Nal),ALoad(Ntr)%Palz(Nal),ALoad(Ntr)%xal(Nal),ALoad(Ntr)%yal(Nal),ALoad(Ntr)%zal(Nal),ALoad(Ntr)%xal0(Nal),ALoad(Ntr)%yal0(Nal),ALoad(Ntr)%zal0(Nal))        
        end if
    end subroutine AllocateActuatorLoadVariables 
    !--------------------------------------------------------------------------------------------
    ! Routine to deallocate all Actuator Load variables
    subroutine DeallocateActuatorLoadVars()
        
        deallocate(ALoad)
        !deallocate(Pac,kxac1,kyac1,xac1,yac1,kxac2,kyac2,xac2,yac2,zac2,Dxac,Dyac,Dzac,lac)
        !deallocate(Palx,Paly,Palz,xal,yal,zal,xal0,yal0,zal0)
    
    end subroutine DeallocateActuatorLoadVars
    
end module ActuatorLoadVariables

!--------------------------------------------------------------------------------------------
! Module with axes transformation variables.
module AxesTransformationVariables
    implicit none
    save

    TYPE :: AxesTransformationVarsType
        real(8) ::			xT(4)					! x-coordinates of trapezoid vertices
        real(8) ::			yT(4)					! y-coordinates of trapezoid vertices
    END TYPE AxesTransformationVarsType
    TYPE(AxesTransformationVarsType), ALLOCATABLE :: ATranformation(:)

contains
    !--------------------------------------------------------------------------------------------
    ! Routine to allocate all axes transformationVars variables
    subroutine AllocatePlanformTransformationVariables(Npl)
        !Input vars
        integer, intent(in) ::		Npl				! number of planforms
        
        if( .not. allocated(ATranformation)) then !First planform
            allocate(ATranformation(Npl))             
        end if

    end subroutine AllocatePlanformTransformationVariables 
    !--------------------------------------------------------------------------------------------
    ! Routine to compute planform axes tranformation coefficients.
    subroutine PlanformTransformationData(Ntr)

        use PlanformVariables, only:xLEroot,xLEroot,xLEtip,croot,ctip,yroot,ytip
        !Input variables
        integer, intent(in) ::          Ntr		! number of the wing trapezoid beeing computed
        
        ATranformation(Ntr)%xT(1)		= xLEroot(Ntr)
        ATranformation(Ntr)%xT(2)		= xLEroot(Ntr)+croot(Ntr)	
        ATranformation(Ntr)%xT(3)		= xLEtip(Ntr)+ctip(Ntr)
        ATranformation(Ntr)%xT(4)		= xLEtip(Ntr)
        ATranformation(Ntr)%yT(1)		= yroot(Ntr)
        ATranformation(Ntr)%yT(2)		= yroot(Ntr)
        ATranformation(Ntr)%yT(3)		= ytip(Ntr)
        ATranformation(Ntr)%yT(4)		= ytip(Ntr)

    end subroutine PlanformTransformationData

end module AxesTransformationVariables

!--------------------------------------------------------------------------------------------
! Module with integration variables.
module IntegrationVariables
    implicit none
    save
    integer ::                  ik          ! number of gauss-legendre integration points
    integer ::                  Mg,Ng		! number of sampling points in zeta,eta directions

end module IntegrationVariables

!--------------------------------------------------------------------------------------------
! Module with section variables.
module StressStrainVariables
    implicit none
    save

    real(8) ::					epsilonsk(5)			! skin strain
    real(8) ::					epsilonsk1(5)			! skin strain
    real(8) ::					epsilonsw(5)			! sparweb strain
    real(8) ::					epsilonsw1(5)			! sparweb strain
    real(8) ::					epsilonsc(5)			! sparcap strain
    real(8) ::					epsilonsc1(5)			! sparcap strain
    real(8) ::					epsilonrw(5)			! ribweb strain
    real(8) ::					epsilonrw1(5)			! ribweb strain
    real(8) ::					epsilonrc(5)			! ribcap strain
    real(8) ::					epsilonrc1(5)			! ribcap strain
    real(8) ::					epsilonst(5)			! stringer strain
    real(8) ::					epsilonst1(5)			! stringer strain

    real(8) ::					sigmavonMisessk			! skin vonMises stress
    real(8) ::					sigmavonMisesskmax      ! max skin vonMises stress
    real(8) ::					sigmavonMisessw			! sparweb vonMises stress
    real(8) ::					sigmavonMisesswmax      ! max sparweb vonMises stress
    real(8) ::					sigmavonMisessc			! sparcap vonMises stress
    real(8) ::					sigmavonMisesscmax      ! max sparcap vonMises stress
    real(8) ::					sigmavonMisesrw			! ribweb vonMises stress
    real(8) ::					sigmavonMisesrwmax      ! max ribweb vonMises stress
    real(8) ::					sigmavonMisesrc			! ribcap vonMises stress
    real(8) ::					sigmavonMisesrcmax      ! max ribcap vonMises stress  
    real(8) ::					sigmavonMisesst			! stringer vonMises stress
    real(8) ::					sigmavonMisesstmax      ! max stringer vonMises stress
  
    real(8) ::					sigmask(5)				! skin stress
    real(8) ::					sigmasw(5)				! sparweb stress
    real(8) ::					sigmasc(5)				! sparcap stress
    real(8) ::					sigmarw(5)				! ribweb stress
    real(8) ::					sigmarc(5)				! ribcap stress
    real(8) ::					sigmast(5)				! stringer stress
contains
    !--------------------------------------------------------------------------------------------
    ! Routine to compute planform axes tranformation coefficients.
    subroutine zero_sigmavonMises()
        sigmavonMisessk = 0.0D0			! skin vonMises stress
        sigmavonMisesskmax = 0.0D0      ! max skin vonMises stress
        sigmavonMisessw = 0.0D0			! sparweb vonMises stress
        sigmavonMisesswmax = 0.0D0      ! max sparweb vonMises stress
        sigmavonMisessc = 0.0D0			! sparcap vonMises stress
        sigmavonMisesscmax = 0.0D0      ! max sparcap vonMises stress
        sigmavonMisesrw = 0.0D0		! ribweb vonMises stress
        sigmavonMisesrwmax  = 0.0D0     ! max ribweb vonMises stress
        sigmavonMisesrc = 0.0D0			! ribcap vonMises stress
        sigmavonMisesrcmax = 0.0D0      ! max ribcap vonMises stress  
        sigmavonMisesst = 0.0D0			! stringer vonMises stress
        sigmavonMisesstmax = 0.0D0      ! max stringer vonMises stress
    end subroutine zero_sigmavonMises

end module StressStrainVariables

!--------------------------------------------------------------------------------------------
! Module with general variables.
module GeneralVariables
    implicit none
    save

    logical ::              PrintMatrix = .FALSE.
    logical ::              PrintVector = .FALSE.
    logical ::              OutputAirfoil = .True.
    logical ::              OutputStress = .True.
    logical ::              OutputDeformation = .TRUE.
    logical ::              wing_definition = .False.
    logical ::              OutputModes = .TRUE.
    real(8), parameter ::   scale = 1.0D0

end module GeneralVariables

!--------------------------------------------------------------------------------------------
!Module with file path variables
module FilePathModule
    implicit none
    save
  
    integer, parameter ::         in_unit = 16        ! input file unit number
    integer, parameter ::         airfoil_unit = 15        ! input file unit number
    integer, parameter ::         static_unit = 3     ! static output file unit number
    integer, parameter ::         wingdef_unit = 4    ! wing definition output file unit number
    integer, parameter ::         LE_TE_unit = 5      ! LE_TE output file unit number
    integer, parameter ::         mshapes_unit = 6    ! mode shapes output file unit number
    integer, parameter ::         skin_stress_unit = 7    ! mode shapes output file unit number
    integer, parameter ::         sparweb_stress_unit = 8    ! mode shapes output file unit number
    integer, parameter ::         sparcap_stress_unit = 9    ! mode shapes output file unit number
    integer, parameter ::         ribweb_stress_unit = 10    ! mode shapes output file unit number
    integer, parameter ::         ribcap_stress_unit = 11    ! mode shapes output file unit number
    integer, parameter ::         stringer_stress_unit = 12    ! mode shapes output file unit number
  
    integer ::                    ierror            ! error number durig opening file
    character(256) ::             file			    ! file name
    character(256) ::             dirproject		! project directory path name
    character(500) ::             output_dir		! project directory path name  
    character(250) ::             msg

contains
    !--------------------------------------------------------------------------------------------
    ! Routine to find next block of data in a file
    subroutine evaluate_next_input()
        !Variable declaration
        character(5) ::			    string			! string

        string=" " 
        do while(string /= "-----") 
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string 
        end do

    end subroutine evaluate_next_input
    !--------------------------------------------------------------------------------------------
    ! Routine to abort program
    subroutine abort_program()
        USE IFCORE, only:PEEKCHARQQ
        !Variable declaration
        LOGICAL(4)                        pressed / .FALSE. /
        
        write(*,'(300A)') 'Error: ',trim(msg)
        write(*,'(A)') 'Press any key to close ...'
        DO WHILE (.NOT. pressed)   
            pressed = PEEKCHARQQ ( ) 
        END DO
        STOP
        end subroutine abort_program
     
end module FilePathModule
