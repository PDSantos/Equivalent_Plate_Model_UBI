!--------------------------------------------------------------------------------------------
! Routine to read skin data.
subroutine update_skin_splines()

    use PlanformVariables, only: Ntr,Section,b2
    use SkinVariablesRoutines, only: Skin
    use maths, only: CoordinatesTransformation
    use SplineMaths, only: SPLINE_function
    implicit none
    
    integer ::                          n_max           ! iteration variable 
    integer ::                          n,i,j           ! iteration variable 
    real(8) ::                          zeta,eta,z      ! point coordinates
    real(8) ::                          x,y
    real(8) ::                          u,v,w         ! point displacements
    
    !Update section z coordinate due to plate deformation
    Ntr = 1
    do n=1,Skin(Ntr)%Nsk !up or lower surface
        do j=1,Section(Ntr)%n_sections !section y position
            y=Section(Ntr)%y_sections(j)
            !!!!!Update section
            do i=1,Section(Ntr)%n_sec_airfoil(n) !airfoil x position (upper or lower surface)
                !x=Skin(Ntr)%xcai_section(n,i,j)
                !x_rel= (x-x_sections(j))/c_sections(j)
                x=Section(Ntr)%x_sections(j)+Skin(Ntr)%Section(n)%xcai(i,j)*Section(Ntr)%c_sections(j) !Transform relative x in absolute coordinate
                Call CoordinatesTransformation(x,y,zeta,eta)  !Compute zeta,eta corresponding to x,y   
                !Call SkinGeometry2(1,n,x,y,z,z1,z2) !calculate z given (x,y)
                z=Skin(Ntr)%Section(n)%zcai(i,j)
                Call compute_displacements(zeta,eta,z,u,v,w) ! Compute point displacements u,v,w
                !Skin(Ntr)%xcai_section(n,1:n_sec_airfoil(n),i)=xcaitip(n,1:n_sec_airfoil(n))
                Skin(Ntr)%Section(n)%zcai(i,j)=Skin(Ntr)%Section(n)%zcai(i,j)+w
            end do
            !!!!!!!!Update z_sections refrence due to deformation
            x=Section(Ntr)%x_sections(j)+0.25D0*Section(Ntr)%c_sections(j) !Take 1/4 chord has wing reference
            z=Section(Ntr)%z_sections(j)
            Call CoordinatesTransformation(x,y,zeta,eta)  !Compute zeta,eta corresponding to x,y   
            Call compute_displacements(zeta,eta,z,u,v,w) ! Compute point displacements u,v,w
            Section(Ntr)%z_sections(j)=Section(Ntr)%z_sections(j)+w
            !!!!!!!Compute centroid point deformation
            x=Section(Ntr)%x_sections_centroid(j) !x centroid in absolute coordinate!!!
            Call CoordinatesTransformation(x,y,zeta,eta)    !Compute zeta,eta corresponding to x,y   
            z=Section(Ntr)%z_sections_centroid(j)
            Call compute_displacements(zeta,eta,z,u,v,w)    ! Compute point displacements u,v,w
            !<<<<<<<subtract centroid point deformation>>>>>>>>>>>>
            Skin(Ntr)%Section(n)%zcai(:,j) = Skin(Ntr)%Section(n)%zcai(:,j)-w
            !-Subtract to z_section reference
            Section(Ntr)%z_sections(j)=Section(Ntr)%z_sections(j)-w
        end do
    end do
    !Compute residual array deformation

    !Update Spline coefficients of modified airfoil data for all sections
    !deallocate(Skin(Ntr)%Section%XS_COEF)
    n_max=maxval(Section(Ntr)%n_sec_airfoil)
    !allocate(Skin(Ntr)%Section(Skin(Ntr)%Nsk)%XS_COEF(n_max,Section(Ntr)%n_sections))
    forall(i=1:Section(Ntr)%n_sections,n=1:Skin(Ntr)%nsk)
        Skin(Ntr)%Section(n)%XS_Coef(1:Section(Ntr)%n_sec_airfoil(n),i)=SPLINE_function(real(Skin(Ntr)%Section(n)%zcai(1:Section(Ntr)%n_sec_airfoil(n),i)),real(Skin(Ntr)%Section(n)%xcai(1:Section(Ntr)%n_sec_airfoil(n),i)),Section(Ntr)%n_sec_airfoil(n))
    end forall
    !!!!Update centroid position (with new splines)
    do i=1,Section(Ntr)%n_sections
        Call SectionCentroid(.False.,Section(Ntr)%y_sections(i)/b2(Ntr),Section(Ntr)%z_sections_centroid(i))
        !write(*,'(F10.6)')  z_sections_centroid(i)     
    end do
  

end subroutine update_skin_splines

! ----------------------------------------------------------------------------------------------
!Subroutine to compute displacements at point zeta,eta
subroutine compute_displacements(zeta,eta,z,u,v,w)

    use PolynomialCoefficients, only: k
    USE maths, only: AxesTranformation

    implicit none

    !Input variables
    real(8), intent(in) ::		zeta,eta	! point transformed coordinates 
    real(8), intent(in) ::		z			! cartesian z coordinate
    !Output variables
    real(8), intent(out) ::		u, v, w		! point displacements 
    real(8) ::    			    x,y			! point cartesian coordinates 
 

    call AxesTranformation(1,zeta,eta,x,y) !x,y matching zeta,eta
    call Deformation_calc(k,zeta,eta,z,u,v,w)

end subroutine compute_displacements