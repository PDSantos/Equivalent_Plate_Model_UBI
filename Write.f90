! ----------------------------------------------------------------------------------------------
! Routine to write the matrices.
subroutine WriteMatrix(K,m,n,name)
    USE FilePathModule, only: output_dir
  implicit none

  integer		m					! number of rows
  integer		n					! number of columns
  real(8)		K(m,n)				! matrix coefficients 

  integer		i,j					! iteration variables

  character(60)	name				! matrix name

  open(unit=1,file=trim(adjustl(output_dir))//"matrix.txt",ACCESS='APPEND')

  write(1,*) ' '
  write(1,*) name
  do i=1,m
    write(1,'(<n>F30.7)') (K(i,j),j=1,n)
  end do

  close(1)

end subroutine WriteMatrix

! ----------------------------------------------------------------------------------------------
! Routine to write the vectors.
subroutine WriteVector(P,n,name)

  implicit none

  integer		n					! number elements
  real(8)		P(n)				! vector coefficients 

  integer		i					! iteration variables

  character(60)	name				! vector name

  open(unit=2,file="vector.txt",ACCESS='APPEND')

  write(2,*) ' '
  write(2,*) name
  do i=1,n
    write(2,'(F11.2)') P(i)*1e-3
  end do

  close(2)

end subroutine WriteVector

! ----------------------------------------------------------------------------------------------
! Routine to write point coordinates and displacements.
!subroutine WriteCoordinatesAndDisplacements2(q,zeta,eta,z)
subroutine WriteCoordinatesAndDisplacements2(unit,zeta,eta,z)

    use DeformationVariables, only: Deformation
    use PolynomialCoefficients, only: k
    use maths, only: AxesTranformation
    use FilePathModule, only: ierror,msg
    

    implicit none
  
    integer, intent(in) ::      unit
    real(8) ::				    x,y,z				! point coordinates 
    real(8) ::				    zeta,eta			! point coordinates 
    real(8) ::                  u,v,w               ! point deformations

    call AxesTranformation(1,zeta,eta,x,y)
    call Deformation_calc(k,zeta,eta,z,u,v,w)
    write(unit,100,IOSTAT=ierror,IOMSG=msg) x,y,z,u,v,w

    100 format(1X,6F15.9)

end subroutine WriteCoordinatesAndDisplacements2

! ----------------------------------------------------------------------------------------------
! Routine to write final point coordinates.
!subroutine WriteCoordinatesAndDisplacements3(q,zeta,eta,z)
subroutine WriteCoordinatesAndDisplacements3(unit,zeta,eta,z)

    use DeformationVariables, only: Deformation
    use PolynomialCoefficients, only: k
    USE maths, only: AxesTranformation
    use FilePathModule, only: ierror,msg

  implicit none

    integer, intent(in) ::      unit
    real(8) ::				    x,y,z				! point coordinates 
    real(8) ::				    zeta,eta			! point coordinates 
    real(8) ::                  u,v,w               ! point deformations

  call AxesTranformation(1,zeta,eta,x,y)
  call Deformation_calc(k,zeta,eta,z,u, v, w)
  write(unit,100,IOSTAT=ierror,IOMSG=msg) x+u,y+v,z+w,u,v,w

  100 format(1X,6F15.9)

end subroutine WriteCoordinatesAndDisplacements3

module write_routines
    implicit none
contains
    ! ----------------------------------------------------------------------------------------------
    ! Routine to write original wing surface and deformed surface.
    subroutine Write_static_tecplot()
    
        use PlanformVariables, only: Ntr,Npl
        use SkinVariablesRoutines, only: SkinZetaTecplot,Skin,SkinGeometry1,SkinGeometry2,SkinEtaTecplot
        USE FilePathModule, only: static_unit,output_dir,ierror,msg
        use maths, only: AxesTranformation
  
        integer ::                  ni,nj				! number of output points along chord and span
        integer ::                  i,j,n
        integer ::                  plate_run_count=1
        real(8) ::                  zeta,zetai,eta,z    ! point coordinates
        real(8) ::                  zetaskF,zetaskA
        real(8) ::                  x,y
        real(8) ::                  z1,z2				! limits of integration
        character(5) ::             string

        ni		= 70
        nj		= 30
        write(string,'(I2.1)') plate_run_count
        if(plate_run_count == 1) then
            open(unit=static_unit,file=trim(adjustl(output_dir))//"static_tecplot.dat",IOSTAT=ierror,IOMSG=msg)
            !Write undeformed surface - only first run
            write(static_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "Static Deformation"'
            write(static_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES = x,y,z,"u, m","v, m","w, m"'
            do Ntr=1, Npl
                write(static_unit,'(A,I2.2,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="undeformed Ntr=',Ntr,'" I=',Skin(Ntr)%Nsk*(ni+1),' J=',nj+1
                do j=0,nj
                    eta = SkinEtaTecplot(nj,j)
                    do n=1,Skin(Ntr)%Nsk
                        call SkinGeometry1(n,eta,zetaskF,zetaskA)		    ! compute zetaskF, zetaskA   
                        do i=0,ni
                            zetai = SkinZetaTecplot(n,ni,i)
                            zeta = 0.5D0*(1.0D0-zetai)*zetaskF+0.5D0*(1.0D0+zetai)*zetaskA
                            call AxesTranformation(1,zeta,eta,x,y)
                            call SkinGeometry2(1,n,x,y,z,z1,z2)
                            call WriteCoordinatesAndDisplacements2(static_unit,zeta,eta,z)
                        end do
                    end do
                end do
            end do
        else
            open(unit=static_unit,file=trim(adjustl(output_dir))//"static_tecplot.dat",position='append',action='write',IOSTAT=ierror,IOMSG=msg)
        end if
        !Output deformed surface - all runs
        do Ntr=1, Npl
            write(static_unit,'(A,I2.2,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="deformed Ntr=',Ntr,'" I=',Skin(Ntr)%Nsk*(ni+1),' J=',nj+1
            !write(static_unit,'(A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="deformed '//trim(string)//'" I=',Skin(Ntr)%Nsk*(ni+1),' J=',nj+1
            do j=0,nj
                eta = SkinEtaTecplot(nj,j)
                do n=1,Skin(Ntr)%Nsk
                    call SkinGeometry1(n,eta,zetaskF,zetaskA)		    ! compute zetaskF, zetaskA   
                    do i=0,ni
                        zetai = SkinZetaTecplot(n,ni,i)
                        zeta = 0.5D0*(1.0D0-zetai)*zetaskF+0.5D0*(1.0D0+zetai)*zetaskA
                        call AxesTranformation(1,zeta,eta,x,y)
                        call SkinGeometry2(1,n,x,y,z,z1,z2)
                        call WriteCoordinatesAndDisplacements3(static_unit,zeta,eta,z)
                    end do
                end do
            end do
        end do
        close(static_unit,IOSTAT=ierror,IOMSG=msg)
    end subroutine Write_static_tecplot

    ! ----------------------------------------------------------------------------------------------
    ! Routine to write wing_surface_definition.txt for MorphingWing.
    subroutine Write_wing_definition()

        use PlanformVariables, only: Ntr
        use SkinVariablesRoutines, only: Skin,SkinZetaTecplot,SkinGeometry1,SkinGeometry2,SkinEtaTecplot
        USE FilePathModule, only: wingdef_unit,LE_TE_unit,output_dir,ierror,msg
        use maths, only: AxesTranformation

        integer ::				ni,nj				! number of output points along chord and span
        integer ::              i,j,n
        real(8) ::				zetai,zeta,eta,z			! point coordinates
        real(8) ::              zetaskF,zetaskA
        real(8) ::				x,y
        real(8) ::				z1,z2				! limits of integration
        
        ni		= 75
        nj		= 20
        open(unit=wingdef_unit,file=trim(adjustl(output_dir))//"wing_surface_definition.dat",IOSTAT=ierror,IOMSG=msg)
        do j=0,nj
            eta = SkinEtaTecplot(nj,j)
            do n=1,Skin(Ntr)%Nsk
                call SkinGeometry1(n,eta,zetaskF,zetaskA)		    ! compute zetaskF, zetaskA   
                do i=0,ni
                    zetai = SkinZetaTecplot(n,ni,i)
                    zeta = 0.5D0*(1.0D0-zetai)*zetaskF+0.5D0*(1.0D0+zetai)*zetaskA
                    call AxesTranformation(1,zeta,eta,x,y)
                    call SkinGeometry2(1,n,x,y,z,z1,z2)
                    call WriteCoordinatesAndDisplacements2(wingdef_unit,zeta,eta,z)
                    write(wingdef_unit,'(4F20.14)') x,y,z,y
                end do
            end do
        end do
        close(wingdef_unit,IOSTAT=ierror,IOMSG=msg)

        ! Write LE_TE_static_tecplot.txt file.
        nj  = 20
        open(unit=LE_TE_unit,file=trim(adjustl(output_dir))//"LE_TE_static_tecplot.dat",IOSTAT=ierror,IOMSG=msg)
        write(LE_TE_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "LE & TE Static Deformation"'
        write(LE_TE_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES = x,y,z,u,v,w'
        write(LE_TE_unit,'(A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="LE" I=',nj+1 !,' J=',nj+1
        do j=0, nj
            zeta = -1.0D0 
            eta	= SkinEtaTecplot(nj,j)
            z = 0.0D0
            call WriteCoordinatesAndDisplacements2(LE_TE_unit,zeta,eta,z)
        end do
        write(LE_TE_unit,'(A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="TE" I=',nj+1,' J=',nj+1
        do j=0, nj
            zeta	= -1.0D0+2.0D0
            eta	= SkinEtaTecplot(nj,j)
            z = 0.0D0
            call WriteCoordinatesAndDisplacements3(LE_TE_unit,zeta,eta,z)
        end do
        close(LE_TE_unit,IOSTAT=ierror,IOMSG=msg)

    end subroutine Write_wing_definition

    ! ----------------------------------------------------------------------------------------------
    ! Routine to write original wing surface and deformed surface.
    subroutine WriteModeShapesTecplot(key,alphar,beta_mkl,VR)

        USE PlanformVariables, only: Npl,Ntr
        use DeformationVariables, only: Deformation
        use SkinVariablesRoutines, only: Skin,SkinGeometry1,SkinGeometry2,SkinZetaTecplot,SkinEtaTecplot
        use FrequenciesVariables, only: freq0
        use FilePathModule, only: mshapes_unit,output_dir,ierror,msg
        USE maths, only: ToHz,AxesTranformation 
        !Input vars
        integer, intent(in) ::              key(:)
        real(8), intent(in) ::              alphar(:),beta_mkl(:),VR(:,:)
        !Local Vars
        integer ::                          k 		                ! k = Npl*5*k*k    
        integer ::                          ni,nj					! number of output points along chord and span
        integer ::                          i,j,l,m,n
        real(8) ::                          vr_vec(size(alphar)),norm
        real(8) ::                          zetai,zeta,eta                ! point coordinates
        real(8) ::                          zetaskF,zetaskA
        real(8) ::                          x,y,z                   ! point coordinates
        real(8) ::                          z1,z2					! limits of integration
    
        k = size(alphar)
        ! Open mode_shape_tecplot.txt file.
        ni		= 50
        nj		= 30
        open(unit=mshapes_unit,file=trim(adjustl(output_dir))//"mode_shapes_tecplot.dat",IOSTAT=ierror,IOMSG=msg)
        write(mshapes_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "Mode Shapes"'
        write(mshapes_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES = x,y,z,u,v,w'
        do Ntr=1, Npl
            write(mshapes_unit,'(A,I2.2,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="undeformed Ntr=',Ntr,'" I=',Skin(Ntr)%Nsk*(ni+1),' J=',nj+1
            Deformation(Ntr)%q	= 0.0D0
            do j=0,nj
                eta = SkinEtaTecplot(nj,j)
                do n=1,Skin(Ntr)%Nsk
                    call SkinGeometry1(n,eta,zetaskF,zetaskA)		    ! compute zetaskF, zetaskA   
                    do i=0,ni
                        zetai = SkinZetaTecplot(n,ni,i) 
                        zeta = 0.5D0*(1.0D0-zetai)*zetaskF+0.5D0*(1.0D0+zetai)*zetaskA
                        call AxesTranformation(1,zeta,eta,x,y)
                        call SkinGeometry2(1,n,x,y,z,z1,z2)
                        call WriteCoordinatesAndDisplacements2(mshapes_unit,zeta,eta,z)
                    end do
                end do
            end do
        end do
        ! Compute natural frequencies.
        ! IMPORTANT NOTE:
        ! EVAL --> Complex array of size N containing the eigenvalues of A in decreasing
        ! order of magnitude.
        ! EVEC --> Complex array containing the matrix of eigenvectors. (Output)
        ! The J-th eigenvector, corresponding to EVAL(J), is stored in the J-th column.
        ! Each vector is normalized to have Euclidean length equal to the value one.
        l = 0
        do i=1,k	       !!!WARNING!!! k_here = Npl*5*k*k 
            IF(ALPHAR(k+1-i)*BETA_mkl(k+1-i) > 0) then
            !IF(BETA_mkl(k+1-i).NE.0.0D0) then
                freq0(i) = sqrt(ALPHAR(k+1-i)/BETA_mkl(k+1-i))
            Else
                freq0(i) = 0.0D0
            End IF           
	        if(freq0(i).EQ.0.0D0) then 
	            l = l
	        else
	            l = l+1
                write(*,'(1X,A,I1.1,A,F9.3,A,F9.3,A)') 'freq(',l,') = ',freq0(i),'rad/s = ',ToHz(freq0(i)),'Hz'
                vr_vec=VR(:,key(k+1-i)) !Select appropriate eigen vector
                norm=sqrt(sum(vr_vec**2)) !Compute eigen vector euclidian norm
                vr_vec=vr_vec/norm !Normalize eigen vector to have euclidian norm=1
                Ntr_loop: do Ntr=1, Npl 
                    ! Write mode_shape_tecplot.txt file.
                    write(mshapes_unit,'(A,I2.2,A,F7.2,A,I2.2,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="Mode ',l,' F=',ToHz(freq0(i)),'Hz Ntr=',Ntr,'" I=',2*(ni+1),' J=',nj+1
                    ! Determine load vector
                    Deformation(Ntr)%q = vr_vec((Ntr - 1)*k/Npl + 1 : Ntr*k/Npl)
                    do j=0,nj
                        eta	= SkinEtaTecplot(nj,j)
                        do n=1, 2 !Only superior and inferior outer skins | Skin(Ntr)%Nsk
                            call SkinGeometry1(n,eta,zetaskF,zetaskA)		    ! compute zetaskF, zetaskA   
                            do m=0,ni
                                zetai = SkinZetaTecplot(n,ni,m)   
                                zeta = 0.5D0*(1.0D0-zetai)*zetaskF+0.5D0*(1.0D0+zetai)*zetaskA
                                call AxesTranformation(1,zeta,eta,x,y)
                                call SkinGeometry2(1,n,x,y,z,z1,z2)
                                call WriteCoordinatesAndDisplacements3(mshapes_unit,zeta,eta,z)
                            end do
                        end do
                    end do
                end do Ntr_loop
            end if
            if(l.GE.6) exit
        end do
        close(mshapes_unit,IOSTAT=ierror,IOMSG=msg)
    end subroutine WriteModeShapesTecplot
    ! ----------------------------------------------------------------------------------------------
    ! Routine to write original wing surface and deformed surface. IMSL VERSION
    subroutine WriteModeShapesTecplotIMSL(BETA,EVAL,EVEC)

        USE PlanformVariables, only: Npl,Ntr
        use DeformationVariables, only: Deformation
        use SkinVariablesRoutines, only: Skin,SkinGeometry1,SkinGeometry2,SkinZetaTecplot,SkinEtaTecplot
        use FrequenciesVariables, only: freq0
        use FilePathModule, only: mshapes_unit,output_dir,ierror,msg
        USE maths, only: ToHz,AxesTranformation 
        !Input vars
        real(8), intent(in)  ::             BETA(:) 
        complex(8), intent(in)  ::          EVAL(:),EVEC(:,:) 
        !Local Vars
        integer ::                          k 		                ! k = Npl*5*k*k    
        integer ::                          ni,nj					! number of output points along chord and span
        integer ::                          i,j,l,m,n
        real(8) ::                          zetai,zeta,eta          ! point coordinates
        real(8) ::                          zetaskF,zetaskA
        real(8) ::                          x,y,z                   ! point coordinates
        real(8) ::                          z1,z2					! limits of integration
    
        k = size(BETA)
        ! Open mode_shape_tecplot.txt file.
        ni		= 50
        nj		= 30
        open(unit=mshapes_unit,file=trim(adjustl(output_dir))//"mode_shapes_tecplot.dat",IOSTAT=ierror,IOMSG=msg)
        write(mshapes_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "Mode Shapes"'
        write(mshapes_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES = x,y,z,u,v,w'
        do Ntr=1, Npl
            write(mshapes_unit,'(A,I2.2,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="undeformed Ntr=',Ntr,'" I=',Skin(Ntr)%Nsk*(ni+1),' J=',nj+1
            Deformation(Ntr)%q	= 0.0D0
            do j=0,nj
                eta = SkinEtaTecplot(nj,j)
                do n=1,Skin(Ntr)%Nsk
                    call SkinGeometry1(n,eta,zetaskF,zetaskA)		    ! compute zetaskF, zetaskA   
                    do i=0,ni
                        zetai = SkinZetaTecplot(n,ni,i) 
                        zeta = 0.5D0*(1.0D0-zetai)*zetaskF+0.5D0*(1.0D0+zetai)*zetaskA
                        call AxesTranformation(1,zeta,eta,x,y)
                        call SkinGeometry2(1,n,x,y,z,z1,z2)
                        call WriteCoordinatesAndDisplacements2(mshapes_unit,zeta,eta,z)
                    end do
                end do
            end do
        end do
        ! Compute natural frequencies.
        ! IMPORTANT NOTE:
        ! EVAL --> Complex array of size N containing the eigenvalues of A in decreasing
        ! order of magnitude.
        ! EVEC --> Complex array containing the matrix of eigenvectors. (Output)
        ! The J-th eigenvector, corresponding to EVAL(J), is stored in the J-th column.
        ! Each vector is normalized to have Euclidean length equal to the value one.
        l = 0
        do i=1,k	       !!!WARNING!!! k_here = Npl*5*k*k 
            IF(BETA(k+1-i).NE.0.0) then
                freq0(i)	= sqrt(eval(k+1-i)/BETA(k+1-i))
            Else
                freq0(i)	= sqrt(eval(k+1-i))
            End IF           
	        if(freq0(i).EQ.0.0D0) then 
	            l = l
	        else
	            l = l+1
                write(*,'(1X,A,I1.1,A,F9.3,A,F9.3,A)') 'freq(',l,') = ',freq0(i),'rad/s = ',ToHz(freq0(i)),'Hz'
                Ntr_loop: do Ntr=1, Npl 
                    ! Write mode_shape_tecplot.txt file.
                    write(mshapes_unit,'(A,I2.2,A,F7.2,A,I2.2,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="Mode ',l,' F=',ToHz(freq0(i)),'Hz Ntr=',Ntr,'" I=',2*(ni+1),' J=',nj+1
                    ! Determine load vector
                    Deformation(Ntr)%q = evec((Ntr - 1)*k/Npl + 1 : Ntr*k/Npl,k+1-i)
                    do j=0,nj
                        eta	= SkinEtaTecplot(nj,j)
                        do n=1, 2 !Only superior and inferior outer skins | Skin(Ntr)%Nsk
                            call SkinGeometry1(n,eta,zetaskF,zetaskA)		    ! compute zetaskF, zetaskA   
                            do m=0,ni
                                zetai = SkinZetaTecplot(n,ni,m)   
                                zeta = 0.5D0*(1.0D0-zetai)*zetaskF+0.5D0*(1.0D0+zetai)*zetaskA
                                call AxesTranformation(1,zeta,eta,x,y)
                                call SkinGeometry2(1,n,x,y,z,z1,z2)
                                call WriteCoordinatesAndDisplacements3(mshapes_unit,zeta,eta,z)
                            end do
                        end do
                    end do
                end do Ntr_loop
            end if
            if(l.GE.6) exit
        end do
        close(mshapes_unit,IOSTAT=ierror,IOMSG=msg)
    end subroutine WriteModeShapesTecplotIMSL

    !***************************************************************************
    ! Subroutine to write data to a file.
    ! The file is to be read by tecplot.
    !***************************************************************************
    subroutine WriteTecplot3D(access,tecplot_file,tecplot_unit,title,variables,zone,x)

      ! Variables declaration.
      ! Input variables.
      integer, intent(in) ::				access			! access (1 = first time, 2 = following times)
      integer, intent(in) ::				tecplot_unit	! tecplot file unit
      character(*), intent(in) ::			tecplot_file	! file that contains the airfoil shape for Tecplot
      character(*), intent(in) ::			title			! graph title
      character(*), intent(in) ::			variables		! variables name
      character(*), intent(in) ::			zone			! zone name
      real, intent(in) ::	                x(:,:,:)        ! variables values

      ! Local variables.
      integer ::						    i,j,k           ! counters

      ! Write data to file, in Tecplot format.
      select case(access)
        case(1)
	      100 continue
          open(tecplot_unit,file=tecplot_file,ERR=100,action='write')
          write(tecplot_unit,'(A)') 'TITLE = "'//trim(title)//'"'
          write(tecplot_unit,'(A)') 'VARIABLES = '//trim(variables) !//'"'
        case(2)
	      200 continue
          open(tecplot_unit,file=tecplot_file,ERR=200,position='append',action='write')
          write(tecplot_unit,'(A,I4.4,A,I4.4)') 'ZONE T="'//trim(zone)//'", I= ',size(x,dim=1),' J= ',size(x,dim=2)
	      do j=1,size(x,dim=2)
	        do i=1,size(x,dim=1)
		      write(tecplot_unit,'(<size(x,dim=3)>E15.7)') (x(i,j,k),k=1,size(x,dim=3))
		    end do
          end do
      end select
      close(tecplot_unit)

    end subroutine WriteTecplot3D
    ! ----------------------------------------------------------------------------------------------
    !Subroutine to export
    subroutine export_displacement_line(tecplot_file,tecplot_unit,title,variables,zone_in,append_on,n,zeta_in,eta_in,factor_u,factor_v)
 
        use PlanformVariables, only: Ntr,Npl
        use IntegrationVariables, only: ik
        use PolynomialCoefficients, only: k
        USE maths, only: AxesTranformation
        use SkinVariablesRoutines, only: SkinZetaTecplot,SkinEtaTecplot,SkinGeometry2
        use FilePathModule, only: abort_program,ierror,msg

        implicit none

        !Input variables
        character(*), intent(in) ::         tecplot_file	! name of file to write data
        integer, intent(in) ::              tecplot_unit 
        character(*), intent(in) ::			title			! graph title
        character(*), intent(in) ::			variables		! variables name
        character(*), intent(in) ::   		zone_in			! zone name
        logical, intent(in) ::           	append_on		! zone name
        integer,intent(in) ::		        n	            ! number of points 
        real(8),optional, intent(in) ::     zeta_in,eta_in  ! point transformed coordinates 
        real(8), intent(in) ::              factor_u,factor_v
        !Local variables
        integer ::                          i                          
        real(8),dimension(n) ::             x,y,z,u,v,w 	! point cartesian coordinates 
        real(8) ::                          zeta,eta
        character(200) ::             		zone
        character(3) ::                		string,Ntrstring
        real(8) ::				            zup,zdown,facx 	! limits of integration
      
        IF(.not.append_on) then
            open(tecplot_unit,file=tecplot_file,ACTION = 'Write')
            write(tecplot_unit,'(A)') 'TITLE = "'//trim(title)//'"'
            write(tecplot_unit,'(A)') 'VARIABLES = '//trim(variables) 
        else
            open(tecplot_unit,file=tecplot_file,ACCESS = 'APPEND')
        end if
      
        zone=zone_in
        write(string,'(I3.2)') k
        If(Present(zeta_in)) then
            zeta = zeta_in
            do Ntr = 1, Npl 
                do i=1,n
                    eta	= SkinEtaTecplot(n-1,i-1)
                    call AxesTranformation(1,zeta,eta,x(i),y(i)) !x,y matching zeta,eta
                    Call SkinGeometry2(1,1,x(i),y(i),zdown) !calculate zdown given (x,y)
                    Call SkinGeometry2(1,2,x(i),y(i),zup) !calculate zup given (x,y)
                    !z(i)=(zup+zdown)/2.0D0
                    z(i)=zup
                    call Deformation_calc(k,zeta,eta,z(i),u(i),v(i),w(i))
                end do
                write(Ntrstring,'(I3.2)') Ntr
                !Write u
                write(tecplot_unit,'(A,I4.4)') 'ZONE T= "u -'//string//' terms Planform'//Ntrstring//'", I= ',n
                do i=1,n  
                    write(tecplot_unit,'(2E15.7)') y(i),u(i)*factor_u
                end do
                !Write v
                write(tecplot_unit,'(A,I4.4)') 'ZONE T= "v -'//string//' terms Planform'//Ntrstring//'", I= ',n
                do i=1,n  
                    write(tecplot_unit,'(2E15.7)') y(i),v(i)*factor_v
                end do
                !Write w
                write(tecplot_unit,'(A,I4.4)') 'ZONE T= "w -'//string//' terms Planform'//Ntrstring//'", I= ' ,n
                do i=1,n  
                    write(tecplot_unit,'(2E15.7)') y(i),w(i)
                end do
            end do !Ntr
        elseif(Present(eta_in)) then
            eta=eta_in
            Ntr = 1
            do i=0,n-1  
                zeta = SkinZetaTecplot(1,n-1,i)
                call AxesTranformation(1,zeta,eta,x(i+1),y(i+1)) !x,y matching zeta,eta
                Call SkinGeometry2(1,1,x(i+1),y(i+1),zdown) !calculate zdown given (x,y)
                Call SkinGeometry2(1,2,x(i+1),y(i+1),zup) !calculate zup given (x,y)
                !z(i+1)=(zup+zdown)/2.0D0
                z(i+1)=zup
                call Deformation_calc(k,zeta,eta,z(i+1),u(i+1),v(i+1),w(i+1))
            end do
            !Adimensionalize x and y
            facx=(x(n)-x(1))
            x=x-x(1)
            !Write u
            zone='ZONE T= "u - '//string//' terms", I= ' 
            write(tecplot_unit,'(A,I4.4)') trim(adjustl(zone)),n
            do i=1,n  
                write(tecplot_unit,'(2E15.7)') x(i),u(i)*factor_u
            end do
            !Write v
            zone='ZONE T= "v - '//string//' terms", I= ' 
            write(tecplot_unit,'(A,I4.4)') trim(adjustl(zone)),n
            do i=1,n  
                write(tecplot_unit,'(2E15.7)') x(i),v(i)*factor_v
            end do
            !Write w
            zone='ZONE T= "w - '//string//' terms", I= ' 
            write(tecplot_unit,'(A,I4.4)') trim(adjustl(zone)),n
            do i=1,n  
                write(tecplot_unit,'(2E15.7)') x(i),w(i)
            end do           
        else
            write(*,*) "Error in export_displacement_line"
        end if
        write(string,'(I2.2)') ik
        close(tecplot_unit)
    IF(ierror.NE.0) Call abort_program()

    end subroutine export_displacement_line
    
    ! ----------------------------------------------------------------------------------------------
    !Subroutine to export 
    subroutine export_test_line(tecplot_file,tecplot_unit,title,variables,zone_in,append_on,n,zeta_in,eta_in)

        use IntegrationVariables, only: ik
        use PolynomialCoefficients, only: k
        USE maths, only: AxesTranformation
        use SkinVariablesRoutines, only: SkinZetaTecplot,SkinEtaTecplot,SkinGeometry2
        USE FilePathModule, only: abort_program,ierror,msg

        implicit none

        !Input variables
        character(*), intent(in) ::         tecplot_file	! name of file to write data
        integer, intent(in) ::              tecplot_unit 
        character(*), intent(in) ::			title			! graph title
        character(*), intent(in) ::			variables		! variables name
        character(*), intent(in) ::   		zone_in			! zone name
        logical, intent(in) ::              append_on		! zone name
        integer,intent(in) ::               n	            ! number of points 
        real(8),optional, intent(in) ::     zeta_in,eta_in  ! point transformed coordinates 
        !Local variables
        integer ::                          i,case_sel                            
        real(8),dimension(n) ::             x,y,z,u,v,w 	! point cartesian coordinates 
        real(8) ::                          zeta,eta
        character(100) ::             		zone
        character(2) ::                		string
        real(8) ::				            zup,zdown,facx 	! limits of integration

        IF(.not.append_on) then
            open(tecplot_unit,file=tecplot_file,ACTION = 'Write',IOSTAT=ierror,IOMSG=msg)
            write(tecplot_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "'//trim(title)//'"'
            write(tecplot_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES = '//trim(variables) 
        else
            open(tecplot_unit,file=tecplot_file,ACCESS = 'APPEND',IOSTAT=ierror,IOMSG=msg)
        end if

        zone=zone_in
        If(Present(zeta_in)) then
            case_sel=1
            zeta=zeta_in
            do i=0,n-1
                eta	= SkinEtaTecplot(n-1,i)
                call AxesTranformation(1,zeta,eta,x(i+1),y(i+1)) !x,y matching zeta,eta
                Call SkinGeometry2(1,1,x(i+1),y(i+1),zdown) !calculate zdown given (x,y)
                Call SkinGeometry2(1,2,x(i+1),y(i+1),zup) !calculate zup given (x,y)
                !z(i+1)=(zup+zdown)/2.0D0
                z(i+1)=zup
                call Deformation_calc(k,zeta,eta,z(i+1),u(i+1),v(i+1),w(i+1))
            end do
        elseif(Present(eta_in)) then
            case_sel=2
            eta=eta_in
            do i=0,n-1  
                zeta = SkinZetaTecplot(1,n-1,i)
                call AxesTranformation(1,zeta,eta,x(i+1),y(i+1)) !x,y matching zeta,eta
                Call SkinGeometry2(1,1,x(i+1),y(i+1),zdown) !calculate zdown given (x,y)
                Call SkinGeometry2(1,2,x(i+1),y(i+1),zup) !calculate zup given (x,y)
                !z(i+1)=(zup+zdown)/2.0D0
                z(i+1)=zup
                call Deformation_calc(k,zeta,eta,z(i+1),u(i+1),v(i+1),w(i+1))
            end do
            !Adimensionalize x and y
            facx=(x(n)-x(1))
            x=x-x(1)
            x=x/facx
        else
            write(*,*) "Error in export_displacement_line"
        end if
        write(string,'(I2.2)') ik
        !Write u
        zone='ZONE T= "u - '//string//' terms", I= ' 
        write(tecplot_unit,'(A,I4.4)') trim(adjustl(zone)),n
        !write(tecplot_unit,'(A,I4.4)') 'ZONE T= "u", I= ',n
        if(case_sel == 1) then
            do i=1,n  
                write(tecplot_unit,'(2E15.7)',IOSTAT=ierror,IOMSG=msg) y(i),z(i)
            end do
        else
            do i=1,n  
                write(tecplot_unit,'(2E15.7)',IOSTAT=ierror,IOMSG=msg) x(i),z(i)
            end do
        end if       
        close(tecplot_unit)
        IF(ierror.NE.0) Call abort_program()
    end subroutine export_test_line
    
    subroutine write_airfoil(tecplot_file,tecplot_unit,Nairoot,xcairoot,zcairoot,Naitip,xcaitip,zcaitip)
        
        use PlanformVariables, only: Ntr
        use SkinVariablesRoutines, only: Skin
        USE FilePathModule, only: abort_program,output_dir,ierror,msg
      
        implicit none

        !Input variables
        character(*), intent(in) ::     tecplot_file	! name of file to write data
        integer, intent(in) ::          tecplot_unit 
        integer, intent(in) ::          Nairoot(:)		        ! number of root airfoil data points - local
        real(8), intent(in) ::          xcairoot(:,:)	        ! abcissas of root airfoil data (x/c) - local
        real(8), intent(in) ::          zcairoot(:,:)	        ! ordinates of root airfoil data (z/c) - local
        integer, intent(in) ::          Naitip(:)	            ! number of tip airfoil data points - local
        real(8), intent(in) ::          xcaitip(:,:)	        ! abcissas of tip airfoil data (x/c) - local
        real(8), intent(in) ::          zcaitip(:,:)	        ! ordinates of tip airfoil data (z/c) - local
        !Local vars
        integer ::                      i,n
    
        open(tecplot_unit,file=trim(adjustl(output_dir))//tecplot_file,ACTION='Write',IOSTAT=ierror,IOMSG=msg)
        write(tecplot_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "Airfoil"'
        write(tecplot_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES = x,y'
        do n=1,Skin(Ntr)%Nsk
            write(tecplot_unit,'(A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="skin" I=',Nairoot(n)
            do i=1,Nairoot(n)
                write(tecplot_unit,'(5F12.7)',IOSTAT=ierror,IOMSG=msg) xcairoot(n,i),zcairoot(n,i)
	        end do
        end do
        do n=1,Skin(Ntr)%Nsk
            write(tecplot_unit,'(A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="skin" I=',Naitip(n)
            do i=1,Naitip(n)
                write(tecplot_unit,'(5F12.7)',IOSTAT=ierror,IOMSG=msg) xcaitip(n,i),zcaitip(n,i)
	        end do
        end do
        close(tecplot_unit,IOSTAT=ierror,IOMSG=msg)
        !Test for write error
        IF(ierror.NE.0) Call abort_program()
    end subroutine write_airfoil
    
end module write_routines


