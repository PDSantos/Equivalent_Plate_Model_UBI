!-----------------------------------------------------------------------------------------------
! Program to calculate deformations of a wing.
subroutine Analysis(choice)

    use PlanformVariables, only: Npl, Ntr
    use PolynomialCoefficients, only: k
  
    implicit none
    !Input Vars
    integer,intent(in) ::               choice
    !Local Vars
    real(8) ::                          yc				! fraction of semi-span
    real(8) ::                          Ixx				! section second moment of inertia
    real(8) ::                          zcentroid		! section z centroid position
    real(8) ::                          wing_mass		! wing mass
    real(8) ::                          flutter_speed
    logical ::                          screen=.True.	! TRUE: write info to screen; False otherwise
    real(8) ::                          FlutterSpeedCalc

    !Compute MassMatrix, StiffnessMatrix or LoadVector for each planform (given the user choice)
    do Ntr=1,Npl
        if (choice.EQ.1) then !Compute Mass, Section Centroid, Deformation, Stresses, Strains and Frequencies.
            yc	= 0.0D0
            ! Mass.
            call Mass(wing_mass)
            ! Section centroid.
            call SectionCentroid(screen,yc,zcentroid)
            ! Section second moment of area.
            call SectionSecondMomentOfArea(screen,yc,Ixx)
            ! Mass matrix.
            call PlanformMassMatrix(k)
            !! Stiffness matrix.
            call PlanformStiffnessMatrix(k)
            ! Load vector.
            call PlanformLoadVector(k)
        else if (choice.EQ.2) then ! Compute Mass.
            ! Mass.
            call Mass(wing_mass)
        else if (choice.EQ.3) then ! Compute Section Centroid.
            ! Section centroid.
            call SectionCentroid(screen,yc,zcentroid)
            ! Section second moment of area.
            call SectionSecondMomentOfArea(screen,yc,Ixx)
        else if (choice.EQ.4) then ! Compute Deformation.
            ! Stiffness matrix.
            call PlanformStiffnessMatrix(k)
            ! Load vector.
            call PlanformLoadVector(k)
        else if (choice.EQ.5) then ! Compute Stresses and Strains.
            ! Stiffness matrix.
            call PlanformStiffnessMatrix(k)
            ! Load vector.
            call PlanformLoadVector(k)
        else if (choice.EQ.6) then ! Compute Wing Frequencies.
            ! Planform Mass matrix.
            call PlanformMassMatrix(k)
            ! Planform Stiffness matrix.
            call PlanformStiffnessMatrix(k)
        else if (choice.EQ.7) then !Flutter Speed Calculation.
            ! Planform Mass matrix.
            call PlanformMassMatrix(k)
            ! Planform Stiffness matrix.
            call PlanformStiffnessMatrix(k)
        end if
    end do
    !Solve Problem
    if (choice.EQ.1) then !Stresses, Strains and Frequencies.
        ! Overall Load vector.
        call OverallLoadVector(k)
        ! Overall Stiffness matrix.
        call OverallStiffnessMatrix(k)
        ! Compute static deformation coefficients.
        call SolveStatic(k)
        ! Compute stresses.
        call Stress()
        ! Compute natural frequencies and mode shapes.
        call SolveEigenvalueProblem(k)
    else if (choice.EQ.2) then ! Compute Mass.
        continue
    else if (choice.EQ.3) then ! Compute Section Centroid.
        continue
    else if (choice.EQ.4) then ! Compute Deformation.
        ! Overall Load vector.
        call OverallLoadVector(k)
        ! Overall Stiffness matrix.
        call OverallStiffnessMatrix(k)
        ! Compute static deformation coefficients.
        call SolveStatic(k)
    else if (choice.EQ.5) then ! Compute Stresses and Strains.
        ! Overall Load vector.
        call OverallLoadVector(k)
        ! Overall Stiffness matrix.
        call OverallStiffnessMatrix(k)
        ! Compute static deformation coefficients.
        call SolveStatic(k)
        ! Compute stresses.
        do Ntr=1,Npl !!Global NTR!!
            call Stress() 
        end do
    else if (choice.EQ.6) then ! Compute Wing Frequencies.
        ! Overall Stiffness matrix.
        call OverallStiffnessMatrix(k)
        ! Overall Mass matrix.
        Call OverallMassMatrix(k)
        ! Compute natural frequencies and mode shapes.
        call SolveEigenvalueProblem(k)
    else if (choice.EQ.7) then !Flutter Speed Calculation.
        ! Overall Stiffness matrix.
        call OverallStiffnessMatrix(k)
        ! Overall Mass matrix.
        Call OverallMassMatrix(k)
        ! Compute natural frequencies and mode shapes.
        call SolveEigenvalueProblem(k)
        ! Compute flutter speed
        flutter_speed=FlutterSpeedCalc(wing_mass)
    end if

end subroutine Analysis

subroutine ReadOverallGeometry()

    use PlanformVariables, only: Npl, Ntr
    use AxesTransformationVariables, only: PlanformTransformationData
    use ReadInputGeometry, only: InputGeometryData
  
    implicit none
    !Input Vars
    !Local Vars

    do Ntr=1,Npl
        !Read all the geometry for planform i
        call InputGeometryData()
        ! Planform transformation coefficients.
        call PlanformTransformationData(Ntr)
    end do

end subroutine ReadOverallGeometry

! ----------------------------------------------------------------------------------------------
! Routine to compute flutter speed based on wing firt torsional frequency
real(8) function FlutterSpeedCalc(wing_mass) result(flutter_speed)
    use FrequenciesVariables, only: freq0
    use PlanformVariables, only: Ntr,b2

    implicit none

    real(8), intent(in) ::              wing_mass
    !local variables
    real(8) ::                          r_alfa=0.5D0,eps_flutter=0.25D0,air_density=1.225D0,dCl_dAlfa=0.0D0


    flutter_speed=freq0(2)*r_alfa/sqrt(2.0D0*dCl_dAlfa*air_density*eps_flutter*b2(Ntr)/wing_mass)
    write(*,'(A,F10.4,A)') ' Flutter speed = ',flutter_speed,' m/s'


end function FlutterSpeedCalc

! ----------------------------------------------------------------------------------------------
! Routine to solve the system of equations K.q-P=0.
subroutine SolveStatic(k)
    USE mkl95_LAPACK, ONLY: GETRF,GETRS,GERFS

    use PlanformVariables, only: Ntr,Npl
    use DeformationVariables, only: Deformation,KM_Global,PV_Global,q_Global
    use SkinVariablesRoutines, only: Skin
    use GeneralVariables, only: OutputDeformation,wing_definition
    USE write_routines, only: Write_static_tecplot,Write_wing_definition

    implicit none
  
    real(8), parameter ::	    errorbound	= 1.0D-8
    !Input variables
    integer, intent(in) ::      k
    !Local variables
    real(8), allocatable  ::	A(:,:)	            ! system of equations
    real(8), allocatable  ::	F(:)			    ! function
    integer ::                  info
    integer ::				    i,j,n			        ! iteration variables
    integer, allocatable  ::    IPVT(:)
    real(8), allocatable  ::    FACT(:,:)
    real(8) ::				    error				! error in moving load coordinates

    allocate(A(Npl*5*k*k,Npl*5*k*k+1),F(Npl*5*k*k),FACT(Npl*5*k*k,Npl*5*k*k),IPVT(Npl*5*k*k))
    !allocate(ferr(1),berr(1),work(3*5*k*k),iwork(5*k*k))
    write(*,*) 'Solving system of linear equations...'

    ! Calculate K.q-P=0.
    n		= 0
    error	= 1.0D0
    !Perform LU Factorization
    FACT = KM_Global !Transfer coefficient matrix (stifeness matrix)
    call GETRF(FACT,ipiv=IPVT,info=INFO) ! Compute LU Factorization
    do i=1,1000
        n = n+1
        q_Global = PV_Global !Transfer RHS data
        call GETRS(FACT,IPVT,q_Global,info=INFO) !Solve Linear System given LU Factorization     
        call GERFS(KM_Global,FACT,IPVT,PV_Global,q_Global,info=INFO) !refine solution
        !Update individual coefficient deformation vector 
        forall( j=1: Npl) Deformation(j)%q = q_Global((j-1)*5*k*k+1 : j*5*k*k )
        !Compute error
        call DetermineError(error)
        IF(error < errorbound) exit
        !Recalculate load vector.
        do Ntr=1, Npl
            call PlanformLoadVector(k)
        end do
        call OverallLoadVector(k)
    end do
    write(*,'(1X,A,I6,1X,A,F13.10)') 'Total number of iterations is: ',n,'Final error: ',error
    write(*,'(1X,A)') '-'
    ! Write static_tecplot.txt file for Tecplot.
    if(OutputDeformation) Call Write_static_tecplot()
    if(wing_definition) Call Write_wing_definition()
    deallocate(A,F,FACT,IPVT)

end subroutine SolveStatic

! ----------------------------------------------------------------------------------------------
! Routine to determine the deformation at a given point - Version without side effects.
!subroutine Deformation(q,zeta,eta,z)
subroutine Deformation_calc(k,zeta,eta,z,u,v,w)
    use PlanformVariables, only: Ntr
    use DeformationVariables, only: Deformation 
    USE maths, only: MultiplyLegendrePolynomials
    implicit none

    !Input Vars
    integer, intent(in) ::      k
    real(8), intent(in) ::      zeta,eta,z			! point coordinates
    !Output Vars
    real(8), intent(out) ::     u,v,w   			! point deformations
    !Local Vars
    integer ::				    i					! iteration variable
    real(8) ::				    BB(k*k)				! Legendre polynomial multiplication vector
    real(8) ::                  phix,phiy
    real(8) ::				    u0,v0,w0

    call MultiplyLegendrePolynomials(1,k,zeta,eta,BB)
    u0	= 0.0D0
    v0	= 0.0D0
    w0	= 0.0D0
    phix = 0.0D0
    phiy = 0.0D0
    do i=1,k*k
        u0		= u0+Deformation(Ntr)%q(i)*BB(i)
        v0		= v0+Deformation(Ntr)%q(i+k*k)*BB(i)
        w0		= w0+Deformation(Ntr)%q(i+2*k*k)*BB(i)
        phix	= phix+Deformation(Ntr)%q(i+3*k*k)*BB(i)
        phiy	= phiy+Deformation(Ntr)%q(i+4*k*k)*BB(i)
    end do
    u   = u0+z*phix
    v   = v0+z*phiy
    w   = w0

end subroutine Deformation_calc

! ----------------------------------------------------------------------------------------------
! Routine to determine the error for the moving load.
!subroutine DetermineError(q,error)
subroutine DetermineError(error)

    use PlanformVariables, only: Npl,Ntr
    use ConcentratedLoadVariables, only: CLoad
    use PolynomialCoefficients, only: k
    USE maths, only: CoordinatesTransformation

    implicit none
  
    !Out Vars
    real(8), intent(out) ::	error				! error bound
    !Local Vars
    real(8) ::				x,y,z				! point coordinates
    real(8) ::				zeta,eta			! point coordinates
    real(8) ::              u,v,w   			! point deformations
    integer ::				i					! iteration variable

    error	= 0.0D0
    do Ntr=1,Npl
        do i=1,CLoad(Ntr)%Ncl
            ! Determine initial reference coordinates.
            x	= CLoad(Ntr)%xcl(i)
            y	= CLoad(Ntr)%ycl(i)
            z	= CLoad(Ntr)%zcl(i)
            ! Determine deformation at load position.
            call CoordinatesTransformation(CLoad(Ntr)%xcl0(i),CLoad(Ntr)%ycl0(i),zeta,eta)
            call Deformation_calc(k,zeta,eta,CLoad(Ntr)%zcl0(i),u,v,w)
            CLoad(Ntr)%xcl(i)	= CLoad(Ntr)%xcl0(i)+u
            CLoad(Ntr)%ycl(i)	= CLoad(Ntr)%ycl0(i)+v
            CLoad(Ntr)%zcl(i)	= CLoad(Ntr)%zcl0(i)+w
            ! Determine error in load coordinates.
	        error = max(error,abs((CLoad(Ntr)%xcl(i)-x)/x),abs((CLoad(Ntr)%ycl(i)-y)/y),abs((CLoad(Ntr)%zcl(i)-z)/z))
        end do
    end do

end subroutine DetermineError

! ----------------------------------------------------------------------------------------------
! Routine to solve the eigenvalue problem (K-lM)q=0.
subroutine SolveEigenvalueProblem(k)
    USE mkl95_LAPACK, ONLY: GGEV
    
    USE PlanformVariables, only: Ntr,Npl  
    use DeformationVariables, only: Deformation,KM_Global,MM_Global
    use GeneralVariables, only: OutputModes
    USE maths, only: PERMUTATION 
    USE write_routines, only: WriteModeShapesTecplot

    implicit none
    !Input vars
    integer, intent(in) ::                k
    !Local Variables 
    real(8),allocatable ::                alphar(:),alphai(:),beta_mkl(:),VR(:,:),eig(:)
    integer,allocatable  ::               key(:)
    integer ::                            INFO = 0
    !Local variables
    integer ::                            i

    write(*,*) 'Solving eigenvalue problem...'
    allocate(alphar(Npl*5*k*k),alphai(Npl*5*k*k),beta_mkl(Npl*5*k*k),VR(Npl*5*k*k,Npl*5*k*k),eig(Npl*5*k*k),key(Npl*5*k*k))
    CALL GGEV(A=KM_Global,B=MM_Global,ALPHAR=ALPHAR, ALPHAI=ALPHAI, BETA=BETA_mkl, VR=VR,INFO=INFO)
    eig=ALPHAR/BETA_mkl
    forall(i=1:Npl*5*k**2)
        key(i)=i
    end forall
    call dlasrt2('D',Npl*5*k**2,eig,key,info) !sorts "eig"  in decreasing order  
    call PERMUTATION (BETA_mkl,key,1)
    call PERMUTATION (ALPHAR,key,1)
    call PERMUTATION (ALPHAI,key,1)

    if(OutputModes) Call WriteModeShapesTecplot(key,alphar,beta_mkl,VR) !output 
    deallocate(alphar,alphai,beta_mkl,VR,eig,key)

end subroutine SolveEigenvalueProblem

