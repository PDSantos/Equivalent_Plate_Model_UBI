module ReadInputGeometry
    implicit none
    !Local variables 
    integer, allocatable ::		Nairoot(:)		        ! number of root airfoil data points - local
    real(8), allocatable ::		xcairoot(:,:)	        ! abcissas of root airfoil data (x/c) - local
    real(8), allocatable ::		zcairoot(:,:)	        ! ordinates of root airfoil data (z/c) - local
    integer, allocatable ::		Naitip(:)	            ! number of tip airfoil data points - local
    real(8), allocatable ::		xcaitip(:,:)	        ! abcissas of tip airfoil data (x/c) - local
    real(8), allocatable ::		zcaitip(:,:)	        ! ordinates of tip airfoil data (z/c) - local

contains
    !-----------------------------------------------------------------------------------------------
    ! Routine to read input data files
    subroutine ReadInputData(choice)

        use FilePathModule, only: abort_program,in_unit,ierror,msg,dirproject,output_dir,file
        use PolynomialCoefficients, only: polynomial,k,AllocatePolynomialCoefficients
        use IntegrationVariables, only: ik
        use IFPORT, only: MAKEDIRQQ
  
        implicit none
        !Output Variables
        integer, intent(out)  ::			choice
        !Local Variables  
        character(90)         ::		    string,string2				! string
        logical                             create_Result
	
        open(1,FILE='project.txt',STATUS='OLD',ACTION='READ',IOSTAT=ierror,IOMSG=msg)
        IF(ierror.NE.0) Call abort_program()
        read(1,'(A)',IOSTAT=ierror,IOMSG=msg) dirproject
        close(1)
        IF(ierror.NE.0) Call abort_program()
  
        !!!!!!!!!!!!!!!!!!!!!!!!!!Creat Output dir!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        output_dir=trim(adjustl(dirproject))//"equivalent_folder\"
        create_Result = MAKEDIRQQ (output_dir)

        file = trim(adjustl(dirproject))//'equivalent_plate_input.txt'
        open(in_unit,FILE=file,STATUS='OLD',ACTION='READ',IOSTAT=ierror,IOMSG=msg)
        IF(ierror.NE.0) Call abort_program()
  
        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) choice
        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) polynomial
        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) k
        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) ik

        Call AllocatePolynomialCoefficients(ik)
	
        select case (choice)
        case (1)
            write(*,*) 'Compute Mass, Section Centroid, Deformation, Stresses, Strains and Frequencies.' 
        case (2)
            write(*,*) 'Compute Mass.' 
        case (3)
            write(*,*) 'Compute Section Centroid.' 
        case (4)
            write(*,*) 'Compute Deformation.' 
        case (5)
            write(*,*) 'Compute Stresses and Strains.' 
        case (6)
            write(*,*) 'Compute Wing Frequencies.'
        case (7)
            write(*,*) 'Flutter Speed Calculation.'
        case default
            write(*,*) 'Invalid calculation selected!'
            Call abort_program()
        end select 
    
	    write(string,*) k
	    write(string2,*) ik
        If(polynomial) Then
           write(*,'(7A)') ' Use Legendre Polynomial of ',trim(adjustl(string)),'th order and ',trim(adjustl(string2)),'X',trim(adjustl(string2)),' integration points' 
        else
           write(*,*) 'Use Ritz Polynomial.' 
        End if 
        write(*,*) '-'
    
    end subroutine ReadInputData



    !--------------------------------------------------------------------------------------------
    ! Routine to read planform shape data.
    subroutine ReadPlanformData(choice)

        use PolynomialCoefficients, only: k
        use DeformationVariables, only: AllocateLocalDeformationVariables
        use FilePathModule, only: in_unit,evaluate_next_input,abort_program,ierror,msg
        use PlanformVariables
        use AxesTransformationVariables
        use FrequenciesVariables, only: AllocateFrequenciesVariables
        use GeneralVariables, only: scale

        implicit none
        !Input Vars
        integer,intent(in) ::       choice
        !Local vars
        integer ::			        n           ! iteration variable
        character(20) ::		    string      ! string

        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Npl
        Npl=Npl-1 ! Npl=Wing sections-1 !!!!
        !Allocate all planform variables
        Call AllocatePlanformVariables()
        Call AllocateLocalDeformationVariables(Npl,k)

        IF(choice == 1 .OR. choice == 7 .OR. choice == 6) Call AllocateFrequenciesVariables(k)

        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) boundary
        do n=1,Npl+1
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            If( n==1 ) then
                read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
                read(in_unit,*,IOSTAT=ierror,IOMSG=msg) yroot(n) 
                read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
                read(in_unit,*,IOSTAT=ierror,IOMSG=msg) xLEroot(n)
                read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
                read(in_unit,*,IOSTAT=ierror,IOMSG=msg) zroot_pl(n)
                read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
                read(in_unit,*,IOSTAT=ierror,IOMSG=msg) croot(n)
                read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
                read(in_unit,*,IOSTAT=ierror,IOMSG=msg) thetaroot(n)
            else
                read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
                read(in_unit,*,IOSTAT=ierror,IOMSG=msg) ytip(n-1)
                read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
                read(in_unit,*,IOSTAT=ierror,IOMSG=msg) xLEtip(n-1)
                read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
                read(in_unit,*,IOSTAT=ierror,IOMSG=msg) ztip_pl(n-1)
                read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
                read(in_unit,*,IOSTAT=ierror,IOMSG=msg) ctip(n-1)
                read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
                read(in_unit,*,IOSTAT=ierror,IOMSG=msg) thetatip(n-1)
                if(n /= Npl+1) then !test for last wing section
                    yroot(n) = ytip(n-1)
                    xLEroot(n) = xLEtip(n-1)
                    zroot_pl(n) = ztip_pl(n-1)
                    croot(n) = ctip(n-1)
                    thetaroot(n) = thetatip(n-1)
                end if
            end if
        end do
        Call ApplyScalePlanform(scale)
        !Test for read error
        IF(ierror.NE.0) Call abort_program()

        !Search for next block of inputs
        Call evaluate_next_input()

        do n=1,Npl
            !Correct xLE coordinates => input x relative to c/4 (and not referenced to leading edge!!) 
            xLEroot(n) = xLEroot(n) - croot(n)/4.0D0
            xLEtip(n) = xLEtip(n) - ctip(n)/4.0D0
            ! Determine semi-span.
            b2(n)		= ytip(n)-yroot(n)
            ! Determine trapezoid leading edge sweep angle.
            Lamdapl(n)	= atan((xLEtip(n)-xLEroot(n))/b2(n))
            ! Determine chord slope.
            mcpl(n)		= (ctip(n)-croot(n))/b2(n)
            Call ComputeRotationMatrix(n)
        end do
        Call AllocatePlanformTransformationVariables(Npl)

    end subroutine ReadPlanformData

    !--------------------------------------------------------------------------------------------
    ! Routine to read skin data.
    subroutine ReadSkinData()

        use PlanformVariables
        use SkinVariablesRoutines
        use FilePathModule, only: in_unit,evaluate_next_input,abort_program,ierror,msg
        use SplineMaths, only: functionLinearInter,SPLINE_function,evaluate_spline,LINEARINTERPOLATION,spline
        use GeneralVariables, only: scale

        implicit none

        integer ::                      n,i 			    ! iteration variable
        integer ::                      Nsk 			    ! number of skin(s) of current planform
        integer ::                      n_max,nl_section
        real(8) ::                      dy                  ! y increment between sections
        real(8) ::                      y1(2),x1(2),z1(2),c1(2)
        real(8), allocatable ::         zaux(:,:)
        character(20) ::                string			    ! string
        logical ::                      usetransverseshear  ! if .TRUE. use transverse shear
        real(4),allocatable ::          XS_COEF(:)
        real(4), allocatable ::		    zcairoot_new(:)     ! ordinates of root airfoil data (z/c) - local
        real(4), allocatable ::         zcaitip_new(:)      ! ordinates of tip airfoil data (z/c) - local
        real(8) ::                      offset(2)           !Top and bottom ofset values

        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string              ! Skin(s) Geometry Definition
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Nsk                     ! no. of skin panels (maximum 10)
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) usetransverseshear      ! T to use transverse shear or F otherwise
        !Allocate memory for skin variables
        Call AllocateSkinVariables(Npl,Ntr,Nsk)

        do n=1,Skin(Ntr)%Nsk
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Skin(Ntr)%kskFroot(n),Skin(Ntr)%kskFtip(n)  ! root and tip chord fraction of skin front line
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Skin(Ntr)%kskAroot(n),Skin(Ntr)%kskAtip(n)  ! root and tip chord fraction of skin aft line

            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Skin(Ntr)%Poly(n)%Nt                        ! no. of polynomial terms defining skin thickness
            allocate(Skin(Ntr)%Poly(n)%Tz(Skin(Ntr)%Poly(n)%Nt))
            allocate(Skin(Ntr)%Poly(n)%n_t(Skin(Ntr)%Poly(n)%Nt))

            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) (Skin(Ntr)%Poly(n)%Tz(i),i=1,Skin(Ntr)%Poly(n)%Nt)          ! polynonial coefficients defining skin thickness
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) (Skin(Ntr)%Poly(n)%n_t(i),i=1,Skin(Ntr)%Poly(n)%Nt)         ! polynomial powers defining skin thickness
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Skin(Ntr)%ttiptroot(n)                                      ! ratio of tip thickness to root thickness

            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Skin(Ntr)%rhosk(n)                                          ! material density
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Skin(Ntr)%Exsk(n),Skin(Ntr)%Eysk(n),Skin(Ntr)%Gxysk(n)      ! in-plane elastic moduli (Ex,Ey,Gxy)
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Skin(Ntr)%miuxysk(n),Skin(Ntr)%miuyxsk(n)                   ! poisson ratio (miuxy,miuyx)
	
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Skin(Ntr)%Gxzsk(n),Skin(Ntr)%Gyzsk(n)                       ! normal elastic moduli (Gxz,Gyz)
        end do
        !Test for read error
        IF(ierror.NE.0) Call abort_program()
        !Search for next block of inputs
        Call evaluate_next_input()
        ! Apply scale
        Call ApplyScaleSkin(scale)     

        ! Read airfoil data points.
        call ReadAirfoilData()

        ! Initialize skin constitutive matrix.
        call SkinConstitutiveMatrix(0,0.0D0,Ntr)
 
        !Generate data for all sections - Created and Modified by Pedro Santos (15/11/2010)
        !generate y coordinate for all sections
        nl_section =  int((ytip(Ntr)-yroot(Ntr)) / 0.01D0)
        if( nl_section <= 3) nl_section = 4
        dy = (ytip(Ntr)-yroot(Ntr)) / dfloat(nl_section-1)
        Call AllocateSectionVariables(Skin(Ntr)%Nsk,nl_section)
        !Compute y increment
        do i=1,Section(Ntr)%n_sections
            Section(Ntr)%y_sections(i)	= dfloat(i-1)*dy
        end do
        y1(1)=Section(Ntr)%y_sections(1)
        y1(2)=Section(Ntr)%y_sections(Section(Ntr)%n_sections)  
        x1(1)=xLEroot(Ntr)
        x1(2)=xLEtip(Ntr)
        z1(1)=zroot_pl(Ntr)
        z1(2)=ztip_pl(Ntr)
        c1(1)=croot(Ntr)
        c1(2)=ctip(Ntr)
        !Compute chords for all sections
        call LinearInterpolation(2,y1,c1,Section(Ntr)%n_sections,Section(Ntr)%y_sections,Section(Ntr)%c_sections)
        !Compute z
        call LinearInterpolation(2,y1,z1,Section(Ntr)%n_sections,Section(Ntr)%y_sections,Section(Ntr)%z_sections)
        !Compute theta
        !call LinearInterpolation(nsections,y1,theta1,jb1,y,theta)
        !Compute x
        call LinearInterpolation(2,y1,x1,Section(Ntr)%n_sections,Section(Ntr)%y_sections,Section(Ntr)%x_sections)

        !Generate interpolant data given xroot or xtip (select x data set with maximum size)
        do n=1,Skin(Ntr)%Nsk
            if(Naitip(n) >= Nairoot(n)) then !tip aifoil has more points
                allocate(XS_COEF(Nairoot(n)),zcairoot_new(1:Naitip(n)))
                Call SPLINE(real(zcairoot(n,1:Nairoot(n))),XS_Coef,real(xcairoot(n,1:Nairoot(n))),Nairoot(n))
                !Generate z data given xtip
                Call evaluate_spline(Nairoot(n),real(xcairoot(n,1:Nairoot(n))),real(zcairoot(n,1:Nairoot(n))),XS_Coef,0,Naitip(n),real(xcaitip(n,1:Naitip(n))),zcairoot_new)
                zcairoot(n,1:Naitip(n))=dble(zcairoot_new(1:Naitip(n)))
                xcairoot(n,1:Naitip(n))=xcaitip(n,1:Naitip(n))
                deallocate(XS_COEF,zcairoot_new)
                Nairoot(n) = Naitip(n)
                Section(Ntr)%n_sec_airfoil(n)=Naitip(n)
            else !root aifoil has more points
                allocate(XS_COEF(Naitip(n)),zcaitip_new(1:Nairoot(n)))
                Call SPLINE(real(zcaitip(n,1:Naitip(n))),XS_Coef,real(xcairoot(n,1:Naitip(n))),Naitip(n))
                !Generate z data given xtip
                Call evaluate_spline(Naitip(n),real(xcaitip(n,1:Naitip(n))),real(zcaitip(n,1:Naitip(n))),XS_Coef,0,Nairoot(n),real(xcairoot(n,1:Nairoot(n))),zcaitip_new)
                zcaitip(n,1:Nairoot(n))=dble(zcaitip_new(1:Nairoot(n)))
                xcaitip(n,1:Nairoot(n))=xcairoot(n,1:Nairoot(n))
                deallocate(XS_COEF,zcaitip_new)
                Naitip(n) = Nairoot(n)
                Section(Ntr)%n_sec_airfoil(n)=Nairoot(n)
            end if
        end do
        deallocate(Nairoot,Naitip)

        !Offset outer airfoil given skin thickness 
        do n=3,Skin(Ntr)%Nsk,2
            !Planform root offset
            offset(1) = Skin(Ntr)%Poly(1)%Tz(1)/2.0D0 !top outer shell thickness
            offset(2) = Skin(Ntr)%Poly(2)%Tz(1)/2.0D0 !bottom outer shell thickness
            do i=3, n,2
                if(i  == n) then
                    offset(1) = offset(1) + Skin(Ntr)%Poly(i)%Tz(1)/2.0D0
                    offset(2) = offset(2) + Skin(Ntr)%Poly(i+1)%Tz(1)/2.0D0
                    offset(1) = offset(1)/croot(Ntr) !normalize
                    offset(2) = offset(2)/croot(Ntr) !normalize
                else
                    offset(1) = offset(1) + Skin(Ntr)%Poly(i)%Tz(1)
                    offset(2) = offset(2) + Skin(Ntr)%Poly(i+1)%Tz(1)
                end if
            end do
            Call OffsetAirfoilData(Section(Ntr)%n_sec_airfoil(1:2),xcairoot(1:2,:),zcairoot(1:2,:),offset,Section(Ntr)%n_sec_airfoil(n:n+1),xcairoot(n:n+1,:),zcairoot(n:n+1,:)) !Root airfoil
            !Update skin start / end position
            Skin(Ntr)%kskFroot(n) = minval(xcairoot(n,1:Section(Ntr)%n_sec_airfoil(n))) !Upper skin Front
            Skin(Ntr)%kskAroot(n) = maxval(xcairoot(n,1:Section(Ntr)%n_sec_airfoil(n))) !Upper skin Back
            Skin(Ntr)%kskFroot(n+1) = minval(xcairoot(n+1,1:Section(Ntr)%n_sec_airfoil(n+1))) !lower skin Front
            Skin(Ntr)%kskAroot(n+1) = maxval(xcairoot(n+1,1:Section(Ntr)%n_sec_airfoil(n+1))) !lower skin Back
            !!!!Planform tip offset
            offset(1) = Skin(Ntr)%Poly(1)%Tz(1)/2.0D0 * Skin(Ntr)%ttiptroot(1) !top outer shell thickness
            offset(2) = Skin(Ntr)%Poly(2)%Tz(1)/2.0D0 * Skin(Ntr)%ttiptroot(2) !bottom outer shell thickness
            do i=3, n,2
                if(i == n) then
                    offset(1) = offset(1) + Skin(Ntr)%Poly(i)%Tz(1)/2.0D0 * Skin(Ntr)%ttiptroot(i)
                    offset(2) = offset(2) + Skin(Ntr)%Poly(i+1)%Tz(1)/2.0D0 * Skin(Ntr)%ttiptroot(i+1)
                    offset(1) = offset(1)/ctip(Ntr) !normalize
                    offset(2) = offset(2)/ctip(Ntr) !normalize
                else
                    offset(1) = offset(1) + Skin(Ntr)%Poly(i)%Tz(1) * Skin(Ntr)%ttiptroot(i)
                    offset(2) = offset(2) + Skin(Ntr)%Poly(i+1)%Tz(1) * Skin(Ntr)%ttiptroot(i+1)
                end if
            end do
            Call OffsetAirfoilData(Section(Ntr)%n_sec_airfoil(1:2),xcaitip(1:2,:),zcaitip(1:2,:),offset,Section(Ntr)%n_sec_airfoil(n:n+1),xcaitip(n:n+1,:),zcaitip(n:n+1,:)) !Tip airfoil
            !Update skin start / end position
            Skin(Ntr)%kskFtip(n) = minval(xcaitip(n,1:Section(Ntr)%n_sec_airfoil(n))) !Upper skin Front
            Skin(Ntr)%kskAtip(n) = maxval(xcaitip(n,1:Section(Ntr)%n_sec_airfoil(n))) !Upper skin Back
            Skin(Ntr)%kskFtip(n+1) = minval(xcaitip(n+1,1:Section(Ntr)%n_sec_airfoil(n+1))) !lower skin Front
            Skin(Ntr)%kskAtip(n+1) = maxval(xcaitip(n+1,1:Section(Ntr)%n_sec_airfoil(n+1))) !lower skin Back
        end do
        !Calculate skin start position and end position coeficents      
        do n=1,Skin(Ntr)%Nsk
            ! Determine thickness slope.
            Skin(Ntr)%mtsk(n)	= (Skin(Ntr)%ttiptroot(n)-1.0D0)/b2(Ntr)
            ! Determine skin front line equation coefficients.
            Skin(Ntr)%kskF1(n)	= (Skin(Ntr)%kskFroot(n)+Skin(Ntr)%kskFtip(n))/2.0D0
            Skin(Ntr)%kskF2(n)	= (Skin(Ntr)%kskFroot(n)-Skin(Ntr)%kskFtip(n))/2.0D0
            ! Determine skin aft line equation coefficients.
            Skin(Ntr)%kskA1(n)	= (Skin(Ntr)%kskAroot(n)+Skin(Ntr)%kskAtip(n))/2.0D0
            Skin(Ntr)%kskA2(n)	= (Skin(Ntr)%kskAroot(n)-Skin(Ntr)%kskAtip(n))/2.0D0
        end do
                
        n_max=maxval(Section(Ntr)%n_sec_airfoil) 
        allocate(zaux(2,n_max))
        !Generate z data for all sections using linear interpolation
        do n=1, Skin(Ntr)%Nsk
            allocate(Skin(Ntr)%Section(n)%xcai(n_max,Section(Ntr)%n_sections))
            allocate(Skin(Ntr)%Section(n)%zcai(n_max,Section(Ntr)%n_sections))
            allocate(Skin(Ntr)%Section(n)%XS_COEF(n_max,Section(Ntr)%n_sections))
            forall(i=1:Section(Ntr)%n_sec_airfoil(n))
                zaux(1,i)=zcairoot(n,i)
                zaux(2,i)=zcaitip(n,i)
                Skin(Ntr)%Section(n)%zcai(i,:) = functionLinearInter(2,y1,zaux(:,i),Section(Ntr)%n_sections,Section(Ntr)%y_sections)
            end forall
        end do
        !Generate x data array for all sections and rescale data set due to chord variation
        forall(i=1:Section(Ntr)%n_sections,n=1:nsk)
            Skin(Ntr)%Section(n)%xcai(1:Section(Ntr)%n_sec_airfoil(n),i) = xcaitip(n,1:Section(Ntr)%n_sec_airfoil(n))
            Skin(Ntr)%Section(n)%zcai(1:Section(Ntr)%n_sec_airfoil(n),i) = Section(Ntr)%z_sections(i) + &
                                                                            & Skin(Ntr)%Section(n)%zcai(1:Section(Ntr)%n_sec_airfoil(n),i)* &
                                                                            & Section(Ntr)%c_sections(i)
        end forall
   
        !Generate Spline coefficients to interpolate airfoil data for all sections
        forall(i=1:Section(Ntr)%n_sections,n=1:nsk)
            Skin(Ntr)%Section(n)%XS_Coef(1:Section(Ntr)%n_sec_airfoil(n),i)=SPLINE_function(real(Skin(Ntr)%Section(n)%zcai(1:Section(Ntr)%n_sec_airfoil(n),i)),real(Skin(Ntr)%Section(n)%xcai(1:Section(Ntr)%n_sec_airfoil(n),i)),Section(Ntr)%n_sec_airfoil(n))
        end forall

        deallocate(xcairoot,zcairoot)
        deallocate(xcaitip,zcaitip)
        deallocate(zaux)

    end subroutine ReadSkinData
    !--------------------------------------------------------------------------------------------
    ! Routine to offset
    subroutine OffsetSectionCentroid()
        use PlanformVariables
        use SkinVariablesRoutines

        implicit none  
        integer ::          i
        
        
    end subroutine OffsetSectionCentroid

    !--------------------------------------------------------------------------------------------
    ! Routine to read airfoil data - Modified by Pedro Santos in 03-11-2010:
    !                                => New airfoil scheme.
    !                               - Modified by Pedro Santos in 28-04-2012:
    !                                => read airfoil per planform
    subroutine ReadAirfoilData()

        use SkinVariablesRoutines
        use PlanformVariables, only: Ntr
        use FilePathModule, only: airfoil_unit,in_unit,evaluate_next_input,abort_program,ierror,msg,file,dirproject
        use GeneralVariables, only: OutputAirfoil
        USE maths, only: dReverseData
        use write_routines, only: write_airfoil

        implicit none

        integer ::                      n,i,j				            !iteration variable
        !!New vars
        integer ::                      maxfoil_new
        integer ::                      max,min,int_aux(1)
        integer ::                      n_airfoil(2),iend(2),ixmin(2)
        real ::                         dummy                           !Dummy variable to find number of airfoil points
        real(8), allocatable ::         xc(:,:),zc(:,:)
        real(8), allocatable ::         xc_aux(:,:),zc_aux(:,:) 
        real(8) ::                      scale
        character(2) ::                 string			                !string

        !Find airfoil number of points
        do j=1,2
            write(string,'(I2.2)') j-1+Ntr
            file = trim(adjustl(dirproject))//'airfoil_coordinates_'//string//'.txt'
            open(unit=airfoil_unit,file=file,STATUS='OLD',ACTION='READ',IOSTAT=ierror,IOMSG=msg)
            read(airfoil_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string     !Airfoil Name
            do i=1,1000
                read(airfoil_unit,*,IOSTAT=ierror) dummy,dummy
                if(ierror < 0) then !end of file detected
                    n_airfoil(j) = i-1 
                    exit
                end if
            end do
            close(airfoil_unit,IOSTAT=ierror,IOMSG=msg)
        end do  
        max = maxval(n_airfoil) !max number of points
        min = minval(n_airfoil) !min number of points
        !nairfoil=max
        allocate(xc(max,2),zc(max,2))
        xc=100.0
        zc=100.0
        !Read airfoil coordinates.
        do j=1,2
            write(string,'(I2.2)') j-1+Ntr
            file = trim(adjustl(dirproject))//'airfoil_coordinates_'//string//'.txt'
            open(unit=airfoil_unit,file=file,STATUS='OLD',ACTION='READ',IOSTAT=ierror,IOMSG=msg)
            read(airfoil_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            do i=1,n_airfoil(j)
                read(airfoil_unit,*,IOSTAT=ierror) xc(i,j),zc(i,j)
            end do
            close(unit=airfoil_unit,IOSTAT=ierror,IOMSG=msg)
            iend(j)  = n_airfoil(j) !index of trailing edge :: before resizing (if any) 
            int_aux  = minloc(xc(:,j)) !index of leading edge :: before resizing (if needed)
            ixmin(j) = int_aux(1) !Transfer value
        end do
        !Test for read error
        IF(ierror.NE.0) Call abort_program()
        !Allocate aifoil variables
        maxfoil_new= max0(maxval(ixmin),maxval(iend-ixmin)+1)
        allocate(Nairoot(Skin(Ntr)%Nsk),xcairoot(Skin(Ntr)%Nsk,maxfoil_new),zcairoot(Skin(Ntr)%Nsk,maxfoil_new))
        allocate(Naitip(Skin(Ntr)%Nsk),xcaitip(Skin(Ntr)%Nsk,maxfoil_new),zcaitip(Skin(Ntr)%Nsk,maxfoil_new))
        !Different number of airfoil points verified in readskindata!!!!!
        !Transfer airfoil data to skinvars
        allocate(xc_aux(max,2),zc_aux(max,2))
        do n=1,Skin(Ntr)%Nsk
            If(n==1) then
                Call dReverseData(xc(1:ixmin(1),1),xc_aux(1:ixmin(1),1))
                Call dReverseData(zc(1:ixmin(1),1),zc_aux(1:ixmin(1),1))
                xcairoot(n,1:ixmin(1))=xc_aux(1:ixmin(1),1)
                zcairoot(n,1:ixmin(1))=zc_aux(1:ixmin(1),1)
                Nairoot(n) = ixmin(1)
            else
                xcairoot(n,1:(iend(1)-ixmin(1)+1))=xc(ixmin(1):iend(1),1)
                zcairoot(n,1:(iend(1)-ixmin(1)+1))=zc(ixmin(1):iend(1),1)
                Nairoot(n) = iend(1)-ixmin(1)+1
            end if
            !Unit chord / LE at x=0
            xcairoot(n,:) = xcairoot(n,:) + abs(xcairoot(n,1))
            scale = xcairoot(n,Nairoot(n))
            xcairoot(n,:) = xcairoot(n,:)/scale
            zcairoot(n,:) = zcairoot(n,:)/scale
        end do
        j=2
        do n=1,Skin(Ntr)%Nsk
            If(n==1) then
                Call dReverseData(xc(1:ixmin(j),j),xc_aux(1:ixmin(j),j))
                Call dReverseData(zc(1:ixmin(j),j),zc_aux(1:ixmin(j),j))
                xcaitip(n,1:ixmin(j))=xc_aux(1:ixmin(j),j)
                zcaitip(n,1:ixmin(j))=zc_aux(1:ixmin(j),j)
                Naitip(n) = ixmin(j)
            else
                xcaitip(n,1:iend(j)-ixmin(j)+1)=xc(ixmin(j):iend(j),j)
                zcaitip(n,1:iend(j)-ixmin(j)+1)=zc(ixmin(j):iend(j),j)
                Naitip(n) = iend(j)-ixmin(j)+1
            end if
            !Unit chord / LE at x=0
            xcaitip(n,:) = xcaitip(n,:) + abs(xcaitip(n,1))
            scale = xcaitip(n,Naitip(n))
            xcaitip(n,:) = xcaitip(n,:)/scale
            zcaitip(n,:) = zcaitip(n,:)/scale
        end do  
        
        ! Write coordinates of airfoil.
        if(OutputAirfoil) Call write_airfoil("airfoil_tecplot.dat",2,Nairoot,xcairoot,zcairoot,Naitip,xcaitip,zcaitip)

        deallocate(xc_aux,zc_aux)

    end subroutine ReadAirfoilData

    !--------------------------------------------------------------------------------------------
    ! Routine to offset airfoil data given geometry and offset constant.    
    subroutine OffsetAirfoilData(Na,Xc,Zc,offset,Na_new,Xc_new,Zc_new)
        USE FilePathModule, only: output_dir,ierror,msg

        !Input vars
        integer, intent(in) ::                  Na(2)
        real(8), intent(in) ::                  Xc(:,:)
        real(8), intent(in) ::                  Zc(:,:)
        real(8), intent(in) ::                  offset(2)
        !Output vars
        integer, intent(out) ::                 Na_new(2)
        real(8), intent(out) ::                 Xc_new(:,:)
        real(8), intent(out) ::                 Zc_new(:,:)
        !In/out vars
        !Local vars
        real(8) ::                      xint,yint
        real(8) ::                      theta
        real(8) ::                      tx,tz,b(2),m(2)
        integer ::                      i,j,n
        integer ::                      i_ref,j_ref
        integer, save ::                count = 0
        
        count = count +1

        do n=1,2
            do i=1,Na(n)  
                !Compute derivative at point
                if( i == Na(n) ) then
                    theta = atan2(-Zc(n,i),Xc(n,i))
                else
                    theta = atan2(-(Zc(n,i+1) - Zc(n,i)),Xc(n,i+1) - Xc(n,i))
                end if
                tx = sin(theta) * offset(n)
                tz = cos(theta) * offset(n)
                if( n==1) then
                    Xc_new(n,i) = Xc(n,i) - tx
                    Zc_new(n,i) = Zc(n,i) - tz
                else
                    Xc_new(n,i) = Xc(n,i) + tx
                    Zc_new(n,i) = Zc(n,i) + tz
                end if
            end do
        end do
        !!! TE correction
        outerTE: do i= Na(1)/2, Na(1) - 1  !Top airfoil  
            do j= Na(2)/2, Na(2) - 1  !Bottom airfoil
                if(Zc_new(1,i) >= Zc_new(2,j) .AND. Zc_new(1,i+1) <= Zc_new(2,j+1)) then
                    if(Xc_new(1,i) <= Xc_new(2,j+1) .AND. Xc_new(1,i+1) >= Xc_new(2,j)) Then
                        i_ref = i 
                        j_ref = j 
                        exit outerTE
                    end if
                end if     
            end do
        end do outerTE
        m(1) = (Zc_new(2,j_ref+1) - Zc_new(2,j_ref)) / (Xc_new(2,j_ref+1) - Xc_new(2,j_ref)) !Bottom
        m(2) = (Zc_new(1,i_ref+1) - Zc_new(1,i_ref)) / (Xc_new(1,i_ref+1) - Xc_new(1,i_ref)) !Top
        b(1) = Zc_new(2,j_ref) - m(1)*Xc_new(2,j_ref)
        b(2) = Zc_new(1,i_ref) - m(2)*Xc_new(1,i_ref)
        xint = -(b(2)-b(1))/(m(2)-m(1))
        yint = m(2)*xint + b(2)
        !! Change TE coordinate
        Na_new(1) = i_ref + 1 
        Na_new(2) = j_ref + 1 
        Xc_new(1,i_ref+1) = xint 
        Xc_new(2,j_ref+1) = xint
        Zc_new(1,i_ref+1) = yint 
        Zc_new(2,j_ref+1) = yint
        
        !!!!! LE correction
        outerLE: do i= Na(1)/2, 2,-1  !Top airfoil  
            do j= Na(2)/2, 2,-1  !Bottom airfoil
                if(Zc_new(1,i) >= Zc_new(2,j) .AND. Zc_new(1,i-1) <= Zc_new(2,j-1)) then
                    if(Xc_new(1,i) >= Xc_new(2,j-1) .And. Xc_new(1,i-1) <= Xc_new(2,j)) Then
                        i_ref = i  !?????
                        j_ref = j + 1
                        exit outerLE
                    end if
                end if        
            end do
        end do outerLE
        m(1) = (Zc_new(2,j_ref-1) - Zc_new(2,j_ref)) / (Xc_new(2,j_ref-1) - Xc_new(2,j_ref)) !Bottom
        m(2) = (Zc_new(1,i_ref-1) - Zc_new(1,i_ref)) / (Xc_new(1,i_ref-1) - Xc_new(1,i_ref)) !Top
        b(1) = Zc_new(2,j_ref) - m(1)*Xc_new(2,j_ref)
        b(2) = Zc_new(1,i_ref) - m(2)*Xc_new(1,i_ref)
        xint = -(b(2)-b(1))/(m(2)-m(1))
        yint = m(2)*xint + b(2)
        !! Change LE coordinate
        Na_new(1) = Na_new(1) - (i_ref - 1) + 1
        Na_new(2) = Na_new(2) - (j_ref - 1) + 1
        Xc_new(1,i_ref-1) = xint 
        Xc_new(2,j_ref-1) = xint
        Zc_new(1,i_ref-1) = yint 
        Zc_new(2,j_ref-1) = yint  
        do i=1, Na_new(1)
            Xc_new(1,i) = Xc_new(1,i + i_ref-2)
            Zc_new(1,i) = Zc_new(1,i + i_ref-2)
        end do
        do i=1, Na_new(2)
            Xc_new(2,i) = Xc_new(2,i + j_ref-2)
            Zc_new(2,i) = Zc_new(2,i + j_ref-2)
        end do 
        if(count == 1) then
            open(2,file=trim(adjustl(output_dir))//"airfoil_tecplot_new.dat",ACTION='Write',IOSTAT=ierror,IOMSG=msg)
            write(2,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "Airfoil"'
            write(2,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES = x,y'
            do n=1,2
                write(2,'(A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="skin ',count,'" I=',Na(n)
                do i=1,Na(n) 
                    write(2,'(5F12.7)',IOSTAT=ierror,IOMSG=msg) Xc(n,i),Zc(n,i)
	            end do
            end do
        else 
            open(2,file=trim(adjustl(output_dir))//"airfoil_tecplot_new.dat",ACTION='Write',ACCESS = 'APPEND',IOSTAT=ierror,IOMSG=msg)
        end if
        do n=1,2
            write(2,'(A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="skin new',count,'" I=',Na_new(n)
            do i=1,Na_new(n) 
                write(2,'(5F12.7)',IOSTAT=ierror,IOMSG=msg) Xc_new(n,i),Zc_new(n,i)
	        end do
        end do
        close(2,IOSTAT=ierror,IOMSG=msg)
    end subroutine OffsetAirfoilData

    !--------------------------------------------------------------------------------------------
    ! Routine to read spar data.
    subroutine ReadSparData()

        use FilePathModule, only: in_unit,evaluate_next_input,abort_program,ierror,msg
        use PlanformVariables, only: Npl,Ntr,b2,croot,ctip,xLEtip,xLEroot
        use SparVariablesRoutines
        use GeneralVariables, only: scale
        implicit none

        integer ::				n				        ! iteration variable
        integer ::            	Nsp                     ! number of spars
        character(20) ::        string			        ! string
        logical				    useinplanenormalshear   ! if .TRUE. use inplane normal shear

        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Nsp                       ! no. of spars (maximum is 10)
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) useinplanenormalshear     ! T to use in-plane normal shear or F otherwise

        Call AllocateSparVariables(Npl,Ntr,Nsp)
 
        do n=1,Spar(Ntr)%Nsp
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spar(Ntr)%ksproot(n),Spar(Ntr)%ksptip(n)    ! chord fraction of spar position at root and tip
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spar(Ntr)%tswroot(n),Spar(Ntr)%tswtip(n)    ! spar web thickness at root and tip

            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spar(Ntr)%hscroot(n),Spar(Ntr)%hsctip(n)    ! spar cap height at root and tip

            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spar(Ntr)%lscroot(n),Spar(Ntr)%lsctip(n)    ! spar cap width at root and tip

            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spar(Ntr)%rhosw(n)                          ! material density of spar web
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spar(Ntr)%rhosc(n)                          ! material density of spar cap
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spar(Ntr)%Gezsw(n),Spar(Ntr)%Eesw(n)        ! transversal and logitudinal elastic moduli of spar web
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spar(Ntr)%miuezsw(n),Spar(Ntr)%miuzesw(n)   ! poisson ratio of spar web (miuez,miuze) - e is direction along web

            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spar(Ntr)%Eesc(n)                           ! logitudinal elastic modulus of spar cap
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spar(Ntr)%miuezsc(n),Spar(Ntr)%miuzesc(n)   ! poisson ratio of spar cap (miuez,miuze) - e is direction along cap
        end do
        !Test for read error
        IF(ierror.NE.0) Call abort_program()
        !Search for next block of inputs
        Call evaluate_next_input()
        !Apply scale
        Call ApplyScaleSpar(scale)
        
        !!
        if(.NOT.useinplanenormalshear) then
            do n=1,Spar(Ntr)%Nsp
                Spar(Ntr)%Eesw(n)	= 0.0D0	
            end do
        end if

        do n=1,Spar(Ntr)%Nsp
            ! Determine spar line sweep angle.
            Spar(Ntr)%Lamdasp(n)	= datan((xLEtip(Ntr)-xLEroot(Ntr)+Spar(Ntr)%ksptip(n)*ctip(Ntr)-Spar(Ntr)%ksproot(n)*croot(Ntr))/b2(Ntr))
            ! Determine spar line equation coefficients.
            Spar(Ntr)%ksp1(n)		= (Spar(Ntr)%ksproot(n)+Spar(Ntr)%ksptip(n))/2.0D0
            Spar(Ntr)%ksp2(n)		= (Spar(Ntr)%ksproot(n)-Spar(Ntr)%ksptip(n))/2.0D0
            ! Determine spar web thickness equation coefficients.
            Spar(Ntr)%tsw1(n)		= (Spar(Ntr)%tswroot(n)+Spar(Ntr)%tswtip(n))/2.0D0
            Spar(Ntr)%tsw2(n)		= (Spar(Ntr)%tswroot(n)-Spar(Ntr)%tswtip(n))/2.0D0
            ! Determine spar cap height equation coefficients.
            Spar(Ntr)%hsc1(n)		= (Spar(Ntr)%hscroot(n)+Spar(Ntr)%hsctip(n))/2.0D0
            Spar(Ntr)%hsc2(n)		= (Spar(Ntr)%hscroot(n)-Spar(Ntr)%hsctip(n))/2.0D0
            ! Determine spar cap width equation coefficients.
            Spar(Ntr)%lsc1(n)		= (Spar(Ntr)%lscroot(n)+Spar(Ntr)%lsctip(n))/2.0D0
            Spar(Ntr)%lsc2(n)		= (Spar(Ntr)%lscroot(n)-Spar(Ntr)%lsctip(n))/2.0D0
        end do

        call SparConstitutiveMatrix(Ntr)

    end subroutine ReadSparData

    !--------------------------------------------------------------------------------------------
    ! Routine to read rib data.
    subroutine ReadRibData()
  
        use FilePathModule, only: in_unit,evaluate_next_input,abort_program,ierror,msg
        use PlanformVariables, only: Npl,Ntr,b2,yroot,croot,ctip,xLEtip,xLEroot,Lamdapl,mcpl
        use RibVariablesRoutines
        use GeneralVariables, only: scale

        implicit none

        integer ::				n				        ! iteration variable
        integer ::            	Nrb                     ! number of rib(s)
        real(8) ::				xle,yle,xte,yte	        ! leading endge and trailing edge rib coordinates
        character(20) ::        string			        ! string
        logical				    useinplanenormalshear   ! if .TRUE. use inplane normal shear

        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*) Nrb                                                 ! no. of ribs          
        read(in_unit,*) useinplanenormalshear                               ! T to use in-plane shear or F otherwise
        Call AllocateRibVariables(Npl,Ntr,Nrb)

        do n=1,Rib(Ntr)%Nrb
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Rib(Ntr)%krble(n),Rib(Ntr)%krbte(n)     ! semi-span fraction of rib position at leading and trailing edges
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Rib(Ntr)%trwle(n),Rib(Ntr)%trwte(n)     ! rib web thickness at leading and trailing edges
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Rib(Ntr)%hrcle(n),Rib(Ntr)%hrcte(n)     ! rib cap height at leading and trailing edges
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Rib(Ntr)%lrcle(n),Rib(Ntr)%lrcte(n)     ! rib cap width at leading and trailing edges
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Rib(Ntr)%rhorw(n)                       ! material density of rib web
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Rib(Ntr)%rhorc(n)                       ! material density of rib cap
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Rib(Ntr)%Gezrw(n),Rib(Ntr)%Eerw(n)      ! transversal and logitudinal elastic moduli of rib web
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Rib(Ntr)%miuezrw(n),Rib(Ntr)%miuzerw(n) ! poisson ratio of rib web (miuez,miuze) - e is direction along rib

            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Rib(Ntr)%Eerc(n)                        ! logitudinal elastic modulus of rib cap
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Rib(Ntr)%miuezrc(n),Rib(Ntr)%miuzerc(n) ! poisson ratio of rib cap (miuez,miuze) - e is direction along rib
        end do
        !Test for read error
        IF(ierror.NE.0) Call abort_program()
        !Search for next block of inputs
        Call evaluate_next_input()
        
        !Apply scale
        Call ApplyScaleRib(scale)
        
        !!
        if(.NOT.useinplanenormalshear) then
            do n=1,Rib(Ntr)%Nrb
                Rib(Ntr)%Eerw(n) = 0.0D0	
            end do
        end if

        do n=1,Rib(Ntr)%Nrb
            ! Determine rib line sweep angle.
            yle			= yroot(Ntr)+Rib(Ntr)%krble(n)*b2(Ntr)
            xle			= xLEroot(Ntr)+Rib(Ntr)%krble(n)*b2(Ntr)*dtan(Lamdapl(Ntr))
            yte			= yroot(Ntr)+Rib(Ntr)%krbte(n)*b2(Ntr)
            xte			= Rib(Ntr)%krbte(n)*b2(Ntr)*dtan(Lamdapl(Ntr))+(croot(Ntr)-mcpl(Ntr)*(yte-yroot(Ntr)))
            Rib(Ntr)%Lamdarb(n)	= datan((yle-yte)/(xte-xle))
            ! Determine rib line equation coefficients.
            Rib(Ntr)%krb1(n)		= (Rib(Ntr)%krble(n)+Rib(Ntr)%krbte(n))/2.0D0
            Rib(Ntr)%krb2(n)		= (Rib(Ntr)%krble(n)-Rib(Ntr)%krbte(n))/2.0D0
            ! Determine rib web thickness equation coefficients.
            Rib(Ntr)%trw1(n)		= (Rib(Ntr)%trwle(n)+Rib(Ntr)%trwte(n))/2.0D0
            Rib(Ntr)%trw2(n)		= (Rib(Ntr)%trwle(n)-Rib(Ntr)%trwte(n))/2.0D0
            ! Determine rib cap height equation coefficients.
            Rib(Ntr)%hrc1(n)		= (Rib(Ntr)%hrcle(n)+Rib(Ntr)%hrcte(n))/2.0D0
            Rib(Ntr)%hrc2(n)		= (Rib(Ntr)%hrcle(n)-Rib(Ntr)%hrcte(n))/2.0D0
            ! Determine rib cap width equation coefficients.
            Rib(Ntr)%lrc1(n)		= (Rib(Ntr)%lrcle(n)+Rib(Ntr)%lrcte(n))/2.0D0
            Rib(Ntr)%lrc2(n)		= (Rib(Ntr)%lrcle(n)-Rib(Ntr)%lrcte(n))/2.0D0
        end do

        call RibConstitutiveMatrix(Ntr)

    end subroutine ReadRibData

    !--------------------------------------------------------------------------------------------
    ! Routine to read stringer data.
    subroutine ReadStringerData()

        use FilePathModule, only: in_unit,evaluate_next_input,abort_program,ierror,msg
        use PlanformVariables, only: Npl,Ntr,b2,croot,ctip,xLEtip,xLEroot
        use StringerVariablesRoutines
        use GeneralVariables, only: scale

        implicit none

        integer ::				n				! iteration variable
        integer ::				Nst				! iteration variable
        character(20) ::        string			! string

        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Nst                         ! no. of stringers (maximum is 20)
        Call AllocateStringerVariables(Npl,Ntr,Nst)
        do n=1,Stringer(Ntr)%Nst
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Stringer(Ntr)%kstroot(n),Stringer(Ntr)%ksttip(n)    ! chord fraction of stringer position at root and tip
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Stringer(Ntr)%hstroot(n),Stringer(Ntr)%hsttip(n)    ! stringer height at root and tip
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Stringer(Ntr)%lstroot(n),Stringer(Ntr)%lsttip(n)    ! stringer width at root and tip
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Stringer(Ntr)%rhost(n)                ! material density
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Stringer(Ntr)%Eest(n)                 ! longitudinal elastic modulus
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Stringer(Ntr)%miuezst(n),Stringer(Ntr)%miuzest(n)   ! poisson ratio of stringer (miuez,miuze) - e is direction along stringer
        end do
        !Test for read error
        IF(ierror.NE.0) Call abort_program()
        !Search for next block of inputs
        Call evaluate_next_input()
        Call ApplyScaleStringer(scale)        
        call StringerConstitutiveMatrix(Ntr)

        do n=1,Stringer(Ntr)%Nst
            ! Determine spar line sweep angle.
            Stringer(Ntr)%Lamdast(n)	= datan((xLEtip(Ntr)-xLEroot(Ntr)+Stringer(Ntr)%ksttip(n)*ctip(Ntr)-Stringer(Ntr)%kstroot(n)*croot(Ntr))/b2(Ntr))
            ! Determine stringer line equation coefficients.
            Stringer(Ntr)%kst1(n)		= (Stringer(Ntr)%kstroot(n)+Stringer(Ntr)%ksttip(n))/2.0D0
            Stringer(Ntr)%kst2(n)		= (Stringer(Ntr)%kstroot(n)-Stringer(Ntr)%ksttip(n))/2.0D0
            ! Determine stringer height equation coefficients.
            Stringer(Ntr)%hst1(n)		= (Stringer(Ntr)%hstroot(n)+Stringer(Ntr)%hsttip(n))/2.0D0
            Stringer(Ntr)%hst2(n)		= (Stringer(Ntr)%hstroot(n)-Stringer(Ntr)%hsttip(n))/2.0D0
            ! Determine stringer width equation coefficients.
            Stringer(Ntr)%lst1(n)		= (Stringer(Ntr)%lstroot(n)+Stringer(Ntr)%lsttip(n))/2.0D0
            Stringer(Ntr)%lst2(n)		= (Stringer(Ntr)%lstroot(n)-Stringer(Ntr)%lsttip(n))/2.0D0
        end do

    end subroutine ReadStringerData

    ! ----------------------------------------------------------------------------------------------
    ! Routine to read boundary spring data.
    subroutine ReadSpringData()

        use FilePathModule, only: in_unit,evaluate_next_input,abort_program,ierror,msg
        use PlanformVariables, only: Npl, Ntr
        use SpringVariables
  
        implicit none

        integer ::				n					! iteration variables
        integer ::				Nspr				! 
        character(20) ::        string				! string

        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        ! Read boundary spring data.
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Nspr
        Call AllocateSpringVariables(Npl,Ntr,Nspr)
        do n=1,Spring(Ntr)%Nspr
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spring(Ntr)%ku(n)
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spring(Ntr)%kv(n)
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spring(Ntr)%kw(n)
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spring(Ntr)%krx(n)
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spring(Ntr)%kry(n)
	        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Spring(Ntr)%xspr(n),Spring(Ntr)%yspr(n)
        end do
        !Test for read error
        IF(ierror.NE.0) Call abort_program()
        !Search for next block of inputs
        Call evaluate_next_input()

    end subroutine ReadSpringData
    ! ----------------------------------------------------------------------------------------------
    ! Routine to read concentrated load data.
    subroutine ReadConcentratedLoadData()

        use FilePathModule, only: in_unit,evaluate_next_input,abort_program,ierror,msg
        use PlanformVariables, only: Npl,Ntr,b2,croot,yroot,xleroot,lamdapl,mcpl
        use ConcentratedLoadVariables

        implicit none

        integer ::				n					! iteration variables
        integer ::				Ncl					! number of concentrated loads (in the current planform)
        character(20) ::        string              ! string
        logical ::              typecoordinate		! .TRUE. = absolute coordinates (x,y,z)
											        ! .FALSE. = relative coordinates [(x-xLE)/c,(y-yROOT)/b/2,z/c]
        ! Read conventrated load data.
        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Ncl
        !Allocate variables
        Call AllocateConcentratedLoadVariables(Npl,Ntr,Ncl)
        if(Ncl > 0) then
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) typecoordinate
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            do n=1,CLoad(Ntr)%Ncl
	            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) CLoad(Ntr)%xcl(n),CLoad(Ntr)%ycl(n),CLoad(Ntr)%zcl(n),CLoad(Ntr)%Pclx(n),CLoad(Ntr)%Pcly(n),CLoad(Ntr)%Pclz(n)
            end do
        else
           write(*,'(3X,A)') '      No load(s) applied to trapezoid!'
        end if
        !Test for read error
        IF(ierror.NE.0) Call abort_program()
        !Search for next block of inputs
        Call evaluate_next_input()
        if(Ncl > 0) then
            if(.NOT.typecoordinate) then
                do n=1,CLoad(Ntr)%Ncl
	                CLoad(Ntr)%ycl(n)	= CLoad(Ntr)%ycl(n)*b2(Ntr)+yroot(Ntr)
	                CLoad(Ntr)%xcl(n)	= CLoad(Ntr)%xcl(n)*(croot(Ntr)+mcpl(Ntr)*(CLoad(Ntr)%ycl(n)-yroot(Ntr)))+xLEroot(Ntr)+(CLoad(Ntr)%ycl(n)-yroot(Ntr))*dtan(Lamdapl(Ntr))
	                CLoad(Ntr)%zcl(n)	= CLoad(Ntr)%zcl(n)*(croot(Ntr)+mcpl(Ntr)*(CLoad(Ntr)%ycl(n)-yroot(Ntr)))
                end do
            end if
            do n=1,CLoad(Ntr)%Ncl
	            CLoad(Ntr)%xcl0(n)	= CLoad(Ntr)%xcl(n)
	            CLoad(Ntr)%ycl0(n)	= CLoad(Ntr)%ycl(n)
	            CLoad(Ntr)%zcl0(n)	= CLoad(Ntr)%zcl(n)
            end do
        end if

    end subroutine ReadConcentratedLoadData
    ! ----------------------------------------------------------------------------------------------
    ! Routine to read actuator load data.
    subroutine ReadActuatorLoadData()

        use FilePathModule, only: in_unit,evaluate_next_input,abort_program,ierror,msg
        use PlanformVariables
        use SkinVariablesRoutines, only:SkinGeometry2
        use ActuatorLoadVariables
        USE maths, only: AxesTranformation

        implicit none

        real(8) ::				zeta,eta			! actuator load coordinates
        real(8) ::				z					! z-coordinate
        real(8) ::				z1,z2				! limits of integration
        integer ::				n					! iteration variables
        integer ::				Nac, Nal            !
        character(20) ::        string				! string

    ! Read actuator load data.
        read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) Nac
        Nal=2*Nac
        call AllocateActuatorLoadVariables(Npl,Ntr,Nac,Nal)
        !read all actuator loads
        do n=1,ALoad(Ntr)%Nac
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,'(A)',IOSTAT=ierror,IOMSG=msg) string
            read(in_unit,*,IOSTAT=ierror,IOMSG=msg) ALoad(Ntr)%Pac(n)
	        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) ALoad(Ntr)%kxac1(n),ALoad(Ntr)%kyac1(n),ALoad(Ntr)%zac1(n)
	        read(in_unit,*,IOSTAT=ierror,IOMSG=msg) ALoad(Ntr)%kxac2(n),ALoad(Ntr)%kyac2(n),ALoad(Ntr)%zac2(n)
        end do
        !Test for read error
        IF(ierror.NE.0) Then
            Call abort_program
        End IF
        if(Ntr == Npl)  then
            close(in_unit,IOSTAT=ierror,IOMSG=msg) !last plaform - close file
        else
            !Search for next block of inputs
            Call evaluate_next_input()
        end if
        
        do n=1,ALoad(Ntr)%Nac
            zeta	= -1.0D0+2.0D0*ALoad(Ntr)%kxac1(n)
            eta		= -1.0D0+2.0D0*ALoad(Ntr)%kyac1(n)
            call AxesTranformation(1,zeta,eta,ALoad(Ntr)%xac1(n),ALoad(Ntr)%yac1(n))
            zeta	= -1.0D0+2.0D0*ALoad(Ntr)%kxac2(n)
            eta		= -1.0D0+2.0D0*ALoad(Ntr)%kyac2(n)
            call AxesTranformation(1,zeta,eta,ALoad(Ntr)%xac2(n),ALoad(Ntr)%yac2(n))
        end do

        do n=1,ALoad(Ntr)%Nac
            z	= ALoad(Ntr)%zac1(n)
            if(z.EQ.1000.0) then
                ! Assume z-coordinate on upper skin.
                call SkinGeometry2(1,1,ALoad(Ntr)%xac1(n),ALoad(Ntr)%yac1(n),ALoad(Ntr)%zac1(n),z1,z2)
	        else if(z.EQ.-1000.0) then
                ! Assume z-coordinate on lower skin.
                call SkinGeometry2(1,2,ALoad(Ntr)%xac1(n),ALoad(Ntr)%yac1(n),ALoad(Ntr)%zac1(n),z1,z2)
	        end if
            z	= ALoad(Ntr)%zac2(n)
            if(z.EQ.1000.0) then
                ! Assume z-coordinate on upper skin.
                call SkinGeometry2(1,1,ALoad(Ntr)%xac2(n),ALoad(Ntr)%yac2(n),ALoad(Ntr)%zac2(n),z1,z2)
	        else if(z.EQ.-1000.0) then
                ! Assume z-coordinate on lower skin.
                call SkinGeometry2(1,2,ALoad(Ntr)%xac2(n),ALoad(Ntr)%yac2(n),ALoad(Ntr)%zac2(n),z1,z2)
	        end if
        end do
        ! Actuator length.
        do n=1,ALoad(Ntr)%Nac
            ALoad(Ntr)%Dxac(n) = ALoad(Ntr)%xac2(n)-ALoad(Ntr)%xac1(n)
            ALoad(Ntr)%Dyac(n) = ALoad(Ntr)%yac2(n)-ALoad(Ntr)%yac1(n)
            ALoad(Ntr)%Dzac(n) = ALoad(Ntr)%zac2(n)-ALoad(Ntr)%zac1(n)
	        ALoad(Ntr)%lac(n)	= dsqrt(ALoad(Ntr)%Dxac(n)*ALoad(Ntr)%Dxac(n)+ALoad(Ntr)%Dyac(n)*ALoad(Ntr)%Dyac(n)+ALoad(Ntr)%Dzac(n)*ALoad(Ntr)%Dzac(n))
        end do
        ! End point coordinates.
        do n=1,ALoad(Ntr)%Nac
	        ALoad(Ntr)%xal(2*n-1)	= ALoad(Ntr)%xac1(n)
	        ALoad(Ntr)%yal(2*n-1)	= ALoad(Ntr)%yac1(n)
	        ALoad(Ntr)%zal(2*n-1)	= ALoad(Ntr)%zac1(n)
	        ALoad(Ntr)%xal(2*n)	= ALoad(Ntr)%xac2(n)
	        ALoad(Ntr)%yal(2*n)	= ALoad(Ntr)%yac2(n)
	        ALoad(Ntr)%zal(2*n)	= ALoad(Ntr)%zac2(n)
        end do
        ! End point loads.
        do n=1,ALoad(Ntr)%Nac
            ALoad(Ntr)%Palx(2*n-1)	= -ALoad(Ntr)%Dxac(n)/ALoad(Ntr)%lac(n)*ALoad(Ntr)%Pac(n)
	        ALoad(Ntr)%Paly(2*n-1)	= -ALoad(Ntr)%Dyac(n)/ALoad(Ntr)%lac(n)*ALoad(Ntr)%Pac(n)
	        ALoad(Ntr)%Palz(2*n-1)	= -ALoad(Ntr)%Dzac(n)/ALoad(Ntr)%lac(n)*ALoad(Ntr)%Pac(n)
            ALoad(Ntr)%Palx(2*n)	=  ALoad(Ntr)%Dxac(n)/ALoad(Ntr)%lac(n)*ALoad(Ntr)%Pac(n)
	        ALoad(Ntr)%Paly(2*n)	=  ALoad(Ntr)%Dyac(n)/ALoad(Ntr)%lac(n)*ALoad(Ntr)%Pac(n)
	        ALoad(Ntr)%Palz(2*n)	=  ALoad(Ntr)%Dzac(n)/ALoad(Ntr)%lac(n)*ALoad(Ntr)%Pac(n)
        end do

        do n=1,ALoad(Ntr)%Nal
	        ALoad(Ntr)%xal0(n)	= ALoad(Ntr)%xal(n)
	        ALoad(Ntr)%yal0(n)	= ALoad(Ntr)%yal(n)
	        ALoad(Ntr)%zal0(n)	= ALoad(Ntr)%zal(n)
	        !write(*,'(1X,3F9.3,3F9.1)') ALoad(Ntr)%xal(n),ALoad(Ntr)%yal(n),ALoad(Ntr)%zal(n),ALoad(Ntr)%Palx(n),ALoad(Ntr)%Paly(n),ALoad(Ntr)%Palz(n)
        end do

    end subroutine ReadActuatorLoadData
    !--------------------------------------------------------------------------------------------
    ! Routine to read integration data.
    subroutine ReadIntegrationData()

        use IntegrationVariables, only: ik,Mg,Ng
        use PolynomialCoefficients, only: k,polynomial
 
        implicit none

        if(polynomial) then
            !Mg		= min(ik,Mgmax)
            !Ng		= min(ik,Ngmax)
            Mg = ik
            Ng = ik
        else
            if(k*k.LT.6) then
                Mg = 6
                Ng = 6
            else if(k*k.GT.10) then
                Mg = 10
                Ng = 10
            else
                Mg = k*k
                Ng = k*k
            end if
        end if

    end subroutine ReadIntegrationData

    !--------------------------------------------------------------------------------------------
    ! Routine to read common data.
    subroutine InputCommonData(choice)
        USE maths, only: QuadratureWeights,LegendrePolynomialsCoefficients
        implicit none

        !Input Vars
        integer,intent(in) ::               choice
        write(*,*) 'Reading input files and computing input data...'

        call ReadPlanformData(choice)
 
        call ReadIntegrationData()
        call LegendrePolynomialsCoefficients()
        call QuadratureWeights()
  

    end subroutine InputCommonData

    !--------------------------------------------------------------------------------------------
    ! Routine to read input geometry data.
    subroutine InputGeometryData()

      use PlanformVariables, only:Ntr

      implicit none

      write(*,'(1X,A,I1)') '    - trapezoid ',Ntr

      call ReadSkinData() !***0)
      call ReadSparData() !***0)
      call ReadRibData() !***0)
      call ReadStringerData() !***0)
      call ReadSpringData()
      call ReadConcentratedLoadData()
      call ReadActuatorLoadData()
      
      !Call OffsetSectionCentroid()

    end subroutine InputGeometryData
    
end module ReadInputGeometry