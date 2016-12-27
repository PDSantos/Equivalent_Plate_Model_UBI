 
!--------------------------------------------------------------------------------------------
! Routine to calculate skin stiffness matrix.
subroutine SkinStiffnessMatrix(k,M)

    use PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: Skin,SkinGeometry1,SkinGeometry2,SkinConstitutiveMatrix
    use IntegrationVariables, only: Mg,Ng
    use GeneralVariables, only: PrintMatrix
    
    USE maths

    implicit none

    integer, intent(in) ::      k					! number of polynomials
    real(8), intent(out) ::     M(5*k*k,5*k*k)  	! stiffness matrix
    real(8) ::				    BB(k*k)				! Legendre polynimial multiplication vector
    real(8) ::				    BdB(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k)
    real(8) ::				    HH11(k*k,k*k),HH12(k*k,k*k),HH21(k*k,k*k),HH22(k*k,k*k),HH145(k*k,k*k)
    real(8) ::				    HH245(k*k,k*k),HH451(k*k,k*k),HH452(k*k,k*k),HH4545(k*k,k*k)
    real(8) ::				    N11(k*k,k*k),N12(k*k,k*k),N13(k*k,k*k),N14(k*k,k*k),N15(k*k,k*k),N21(k*k,k*k)
    real(8) ::				    N22(k*k,k*k),N23(k*k,k*k),N24(k*k,k*k),N25(k*k,k*k),N31(k*k,k*k),N32(k*k,k*k)
    real(8) ::				    N33(k*k,k*k),N34(k*k,k*k),N35(k*k,k*k),N41(k*k,k*k),N42(k*k,k*k),N43(k*k,k*k)
    real(8) ::				    N44(k*k,k*k),N51(k*k,k*k),N52(k*k,k*k),N53(k*k,k*k),N55(k*k,k*k)
    !This variables are now local
    real(8) ::				    zetaskF,zetaskA
    real(8) ::                  slope               ! Skin dz/dy
    real(8) ::				    DJ,J1(2,2)
    real(8) ::				    I1,I2,I3			! integration constants
    real(8) ::				    x,y,z				! point coordinates
    real(8) ::				    zeta				! zeta-coordinate
    real(8) ::				    z1,z2				! limits of integration
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    
    integer ::				    i,j,n,m_g,n_g		! iteration variables	  
    character(60) ::            name				! matrix name

    M	= 0.0D0
    if(Skin(Ntr)%Nsk.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng)) 
   
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Skin(Ntr)%Nsk
                call SkinGeometry1(n,etaj(n_g),zetaskF,zetaskA)		    ! compute zetaskF, zetaskA   
                zeta = 0.5D0*(1.0D0-zetai(m_g))*zetaskF+0.5D0*(1.0D0+zetai(m_g))*zetaskA
                call AxesTranformation(1,zeta,etaj(n_g),x,y)		    ! compute x, y
                call AxesTranformation(3,zeta,etaj(n_g),x,y,DJ,J1)		! compute DJ, J1
                call SkinGeometry2(3,n,x,y,z,z1,z2,slope)				! compute z1, z2 and slope												
                call TrialPolynomials(1,k,zeta,etaj(n_g),BB)
                call TrialPolynomials(2,k,zeta,etaj(n_g),BdB)
                call TrialPolynomials(3,k,zeta,etaj(n_g),BBd)
                TC11 = J1(1,1)*BdB+J1(2,1)*BBd
                TC22 = J1(1,2)*BdB+J1(2,2)*BBd
                call MultiplyVectorByTransposedVector(TC11,HH11)
                call MultiplyVectorByVector(TC11,TC22,HH12)
                call MultiplyVectorByVector(TC22,TC11,HH21)
                call MultiplyVectorByTransposedVector(TC22,HH22)
                call MultiplyVectorByVector(TC11,BB,HH145)
                call MultiplyVectorByVector(TC22,BB,HH245)
                call MultiplyVectorByVector(BB,TC11,HH451)
                call MultiplyVectorByVector(BB,TC22,HH452)  
                call MultiplyVectorByTransposedVector(BB,HH4545)
                call SkinConstitutiveMatrix(1,datan(slope),Ntr,n)
                N11=Skin(Ntr)%D11sk(n)*HH11+Skin(Ntr)%D44sk(n)*HH22
                N12=Skin(Ntr)%D12sk(n)*HH12+Skin(Ntr)%D44sk(n)*HH21
                N13=Skin(Ntr)%D16sk(n)*HH11+Skin(Ntr)%D45sk(n)*HH22
                N14=Skin(Ntr)%D16sk(n)*HH145
                N15=Skin(Ntr)%D45sk(n)*HH245
                
                N21=Skin(Ntr)%D21sk(n)*HH21+Skin(Ntr)%D44sk(n)*HH12
                N22=Skin(Ntr)%D22sk(n)*HH22+Skin(Ntr)%D44sk(n)*HH11
                N23=Skin(Ntr)%D26sk(n)*HH21+Skin(Ntr)%D45sk(n)*HH12	
                N24=Skin(Ntr)%D26sk(n)*HH245
                N25=Skin(Ntr)%D45sk(n)*HH145

                N31=Skin(Ntr)%D54sk(n)*HH22+Skin(Ntr)%D61sk(n)*HH11
                N32=Skin(Ntr)%D54sk(n)*HH21+Skin(Ntr)%D62sk(n)*HH12
                N33=Skin(Ntr)%D55sk(n)*HH22+Skin(Ntr)%D66sk(n)*HH11
                N34=Skin(Ntr)%D66sk(n)*HH145
                N35=Skin(Ntr)%D55sk(n)*HH245

                N41=Skin(Ntr)%D61sk(n)*HH451
                N42=Skin(Ntr)%D62sk(n)*HH452
                N43=Skin(Ntr)%D66sk(n)*HH451
                N44=Skin(Ntr)%D66sk(n)*HH4545

                N51=Skin(Ntr)%D54sk(n)*HH452
                N52=Skin(Ntr)%D54sk(n)*HH451
                N53=Skin(Ntr)%D55sk(n)*HH452
                N55=Skin(Ntr)%D55sk(n)*HH4545
                !Integration constants
                I1	= 0.5D0*(zetaskA-zetaskF)*(z2-z1)*DJ
                I2	= 0.5D0*(zetaskA-zetaskF)*(z2*z2-z1*z1)/2.0D0*DJ
                I3	= 0.5D0*(zetaskA-zetaskF)*(z2*z2*z2-z1*z1*z1)/3.0D0*DJ
                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j)				= M(i,j)+gi(m_g)*gj(n_g)*N11(i,j)*I1
                    ! M12.
                    M(i,j+k*k)			= M(i,j+k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I1
                    ! M13.
                    M(i,j+2*k*k)		= M(i,j+2*k*k)+gi(m_g)*gj(n_g)*N13(i,j)*I1
                    ! M14.
                    M(i,j+3*k*k)		= M(i,j+3*k*k)+gi(m_g)*gj(n_g)*(N11(i,j)*I2+N14(i,j)*I1)
                    ! M15.
                    M(i,j+4*k*k)		= M(i,j+4*k*k)+gi(m_g)*gj(n_g)*(N12(i,j)*I2+N15(i,j)*I1)
                    ! M21.
                    M(i+k*k,j)			= M(i+k*k,j)+gi(m_g)*gj(n_g)*N21(i,j)*I1
                    ! M22.
                    M(i+k*k,j+k*k)		= M(i+k*k,j+k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I1
                    ! M23.
                    M(i+k*k,j+2*k*k)	= M(i+k*k,j+2*k*k)+gi(m_g)*gj(n_g)*N23(i,j)*I1
                    ! M24.
                    M(i+k*k,j+3*k*k)	= M(i+k*k,j+3*k*k)+gi(m_g)*gj(n_g)*(N21(i,j)*I2+N24(i,j)*I1)
                    ! M25.
                    M(i+k*k,j+4*k*k)	= M(i+k*k,j+4*k*k)+gi(m_g)*gj(n_g)*(N22(i,j)*I2+N25(i,j)*I1)
                    ! M31.
                    M(i+2*k*k,j)		= M(i+2*k*k,j)+gi(m_g)*gj(n_g)*N31(i,j)*I1
                    ! M32.
                    M(i+2*k*k,j+k*k)	= M(i+2*k*k,j+k*k)+gi(m_g)*gj(n_g)*N32(i,j)*I1
                    ! M33.
                    M(i+2*k*k,j+2*k*k)	= M(i+2*k*k,j+2*k*k)+gi(m_g)*gj(n_g)*N33(i,j)*I1
                    ! M34.
                    M(i+2*k*k,j+3*k*k)	= M(i+2*k*k,j+3*k*k)+gi(m_g)*gj(n_g)*(N31(i,j)*I2+N34(i,j)*I1)
                    ! M35.
                    M(i+2*k*k,j+4*k*k)	= M(i+2*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*(N32(i,j)*I2+N35(i,j)*I1)
                    ! M41.
                    M(i+3*k*k,j)		= M(i+3*k*k,j)+gi(m_g)*gj(n_g)*(N11(i,j)*I2+N41(i,j)*I1)
                    ! M42.
                    M(i+3*k*k,j+k*k)	= M(i+3*k*k,j+k*k)+gi(m_g)*gj(n_g)*(N12(i,j)*I2+N42(i,j)*I1)
                    ! M43.
                    M(i+3*k*k,j+2*k*k)	= M(i+3*k*k,j+2*k*k)+gi(m_g)*gj(n_g)*(N13(i,j)*I2+N43(i,j)*I1)
                    ! M44.
                    M(i+3*k*k,j+3*k*k)	= M(i+3*k*k,j+3*k*k)+gi(m_g)*gj(n_g)*(N11(i,j)*I3+(N14(i,j)+N41(i,j))*I2+N44(i,j)*I1)
                    ! M45.
                    M(i+3*k*k,j+4*k*k)	= M(i+3*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*(N12(i,j)*I3+(N15(i,j)+N42(i,j))*I2)
                    ! M51.
                    M(i+4*k*k,j)		= M(i+4*k*k,j)+gi(m_g)*gj(n_g)*(N21(i,j)*I2+N51(i,j)*I1)
                    ! M52.
                    M(i+4*k*k,j+k*k)	= M(i+4*k*k,j+k*k)+gi(m_g)*gj(n_g)*(N22(i,j)*I2+N52(i,j)*I1)
                    ! M53.
                    M(i+4*k*k,j+2*k*k)	= M(i+4*k*k,j+2*k*k)+gi(m_g)*gj(n_g)*(N23(i,j)*I2+N53(i,j)*I1)
                    ! M54.
                    M(i+4*k*k,j+3*k*k)	= M(i+4*k*k,j+3*k*k)+gi(m_g)*gj(n_g)*(N21(i,j)*I3+(N24(i,j)+N51(i,j))*I2)
                    ! M55.
                    M(i+4*k*k,j+4*k*k)	= M(i+4*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*(N22(i,j)*I3+(N25(i,j)+N52(i,j))*I2+N55(i,j)*I1)
                end forall
            end do
        end do
    end do
    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'skin stiffness matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if
    deallocate(gi,gj,zetai,etaj)  
    return

end subroutine SkinStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate spar web stiffness matrix.
subroutine SparWebStiffnessMatrix(k,M)

    use PlanformVariables, only: Ntr
    use SparVariablesRoutines, only: Spar,SparGeometry,SparWebGeometry
    use IntegrationVariables, only: Mg,Ng
    use GeneralVariables, only: PrintMatrix
    USE maths

    implicit none

    integer, intent(in)   ::    k					! number of polynomials
    real(8), intent(out)  ::    M(5*k*k,5*k*k)  	! stiffness matrix
    real(8) ::				    BB(k*k)				! Legendre polynimial multiplication vector
    real(8) ::				    BdB(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k)
    real(8) ::				    HH11(k*k,k*k),HH12(k*k,k*k),HH145(k*k,k*k)
    real(8) ::				    HH21(k*k,k*k),HH22(k*k,k*k),HH245(k*k,k*k),HH4545(k*k,k*k)
    real(8) ::				    N11(k*k,k*k),N12(k*k,k*k),N22(k*k,k*k),N33(k*k,k*k),N34(k*k,k*k)
    real(8) ::				    N35(k*k,k*k),N44(k*k,k*k),N45(k*k,k*k),N55(k*k,k*k)
    real(8) ::				    I1,I2,I3			! integration constants
    real(8) ::				    x,y					! point coordinates
    real(8) ::				    zeta				! zeta-coordinate
    real(8) ::				    z1,z2				! limits of integration
    real(8) ::                  DJ,J1(2,2)
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    
    integer ::				    i,j,n,m_g,n_g		! iteration variables
    character(60) ::			name				! matrix name

    M	= 0.0D0

    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng)) 
   
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Spar(Ntr)%Nsp
                call SparWebGeometry(1,n,zetai(m_g),etaj(n_g),zeta,x,y,z1,z2)	! compute Spar(Ntr)%tcsw, zeta
                call AxesTranformation(1,zeta,etaj(n_g),x,y,DJ,J1)				! compute x, y
                call AxesTranformation(3,zeta,etaj(n_g),x,y,DJ,J1)				! compute DJ, J1
                call SparWebGeometry(2,n,zetai(m_g),etaj(n_g),zeta,x,y,z1,z2)	! compute z1, z2

                call TrialPolynomials(1,k,zeta,etaj(n_g),BB) !TC45
                call TrialPolynomials(2,k,zeta,etaj(n_g),BdB)
                call TrialPolynomials(3,k,zeta,etaj(n_g),BBd)
                TC11=J1(1,1)*BdB+J1(2,1)*BBd
                TC22=J1(1,2)*BdB+J1(2,2)*BBd

                call MultiplyVectorByTransposedVector(TC11,HH11)
                call MultiplyVectorByVector(TC11,TC22,HH12)
                call MultiplyVectorByVector(TC22,TC11,HH21)
                call MultiplyVectorByTransposedVector(TC22,HH22)
                call MultiplyVectorByVector(TC11,BB,HH145)
                call MultiplyVectorByVector(TC22,BB,HH245)
                call MultiplyVectorByTransposedVector(BB,HH4545)

                N11= Spar(Ntr)%D11sw(n)*HH11+Spar(Ntr)%D31sw(n)*HH21+Spar(Ntr)%D13sw(n)*HH12+Spar(Ntr)%D33sw(n)*HH22
                N12=Spar(Ntr)%D12sw(n)*HH12+Spar(Ntr)%D32sw(n)*HH22+Spar(Ntr)%D13sw(n)*HH11+Spar(Ntr)%D33sw(n)*HH21
                N22=Spar(Ntr)%D22sw(n)*HH22+Spar(Ntr)%D32sw(n)*HH12+Spar(Ntr)%D23sw(n)*HH21+Spar(Ntr)%D33sw(n)*HH11
                N33=Spar(Ntr)%D44sw(n)*HH22+Spar(Ntr)%D54sw(n)*HH12+Spar(Ntr)%D45sw(n)*HH21+Spar(Ntr)%D55sw(n)*HH11
                N34=Spar(Ntr)%D45sw(n)*HH245+Spar(Ntr)%D55sw(n)*HH145
                N35=Spar(Ntr)%D44sw(n)*HH245+Spar(Ntr)%D54sw(n)*HH145
                N44=Spar(Ntr)%D55sw(n)*HH4545
                N45=Spar(Ntr)%D54sw(n)*HH4545
                N55=Spar(Ntr)%D44sw(n)*HH4545

                ! Integration constants.		
                I1	= Spar(Ntr)%tcsw*(z2-z1)*DJ
                I2	= Spar(Ntr)%tcsw*(z2*z2-z1*z1)/2.0D0*DJ
                I3	= Spar(Ntr)%tcsw*(z2*z2*z2-z1*z1*z1)/3.0D0*DJ

                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j)= M(i,j)+gi(m_g)*gj(n_g)*N11(i,j)*I1
                    ! M12.
                    M(i,j+k*k)	= M(i,j+k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I1
                    ! M13.=0
                    ! M14.
                    M(i,j+3*k*k)= M(i,j+3*k*k)+gi(m_g)*gj(n_g)*N11(i,j)*I2
                    ! M15.
                    M(i,j+4*k*k)= M(i,j+4*k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I2
                    ! M21.
                    M(i+k*k,j)			= M(j,i+k*k)
                    ! M22.
                    M(i+k*k,j+k*k)		= M(i+k*k,j+k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I1		   
                    ! M23.=0
                    ! M24.
                    M(i+k*k,j+3*k*k)	= M(j,i+4*k*k)
                    ! M25.
                    M(i+k*k,j+4*k*k)	= M(i+k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I2
                    ! M31.=0
                    ! M32.=0
                    ! M33.
                    M(i+2*k*k,j+2*k*k)	= M(i+2*k*k,j+2*k*k)+gi(m_g)*gj(n_g)*N33(i,j)*I1	
                    ! M34.
                    M(i+2*k*k,j+3*k*k)	= M(i+2*k*k,j+3*k*k)+gi(m_g)*gj(n_g)*N34(i,j)*I1
                    ! M35.
                    M(i+2*k*k,j+4*k*k)	= M(i+2*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N35(i,j)*I1
                    ! M41.
                    M(i+3*k*k,j)		= M(j,i+3*k*k)
                    ! M42.
                    M(i+3*k*k,j+k*k)	= M(j+k*k,i+3*k*k)
                    ! M43.
                    M(i+3*k*k,j+2*k*k)	= M(j+2*k*k,i+3*k*k)
                    ! M44.
                    M(i+3*k*k,j+3*k*k)	= M(i+3*k*k,j+3*k*k)+gi(m_g)*gj(n_g)*(N11(i,j)*I3+N44(i,j)*I1)   
                    ! M45.
                    M(i+3*k*k,j+4*k*k)	= M(i+3*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*(N12(i,j)*I3+N45(i,j)*I1)
                    ! M51.
                    M(i+4*k*k,j)		= M(j,i+4*k*k)
                    ! M52.
                    M(i+4*k*k,j+k*k)	= M(j+k*k,i+4*k*k)
                    ! M53.
                    M(i+4*k*k,j+2*k*k)	= M(j+2*k*k,i+4*k*k)
                    ! M54.
                    M(i+4*k*k,j+3*k*k)	= M(j+3*k*k,i+4*k*k)
                    ! M55.
                    M(i+4*k*k,j+4*k*k)	= M(i+4*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*(N22(i,j)*I3+N55(i,j)*I1)
                end forall
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)  
    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'spar web stiffness matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if

end subroutine SparWebStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate spar cap stiffness matrix.
subroutine SparCapStiffnessMatrix(k,M)

    use PlanformVariables, only: Ntr
    use SparVariablesRoutines, only: Spar,SparGeometry,SparCapGeometry
    use IntegrationVariables, only: Mg,Ng
    use GeneralVariables, only: PrintMatrix
    USE maths

    implicit none

    integer, intent(in)   ::    k					! number of polynomials
    real(8), intent(out)  ::    M(5*k*k,5*k*k)  	! stiffness matrix
    real(8) ::				    BdB(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k)
    real(8) ::				    HH11(k*k,k*k),HH12(k*k,k*k),HH21(k*k,k*k),HH22(k*k,k*k)
    real(8) ::				    N11(k*k,k*k),N12(k*k,k*k),N22(k*k,k*k),N24(k*k,k*k)
    real(8) ::				    I1,I2,I3			! integration constants
    real(8) ::				    x,y					! point coordinates
    real(8) ::				    zeta				! zeta-coordinate
    real(8) ::				    z1U,z2U,z1L,z2L		! limits of integration
    real(8) ::                  DJ,J1(2,2)
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    
    integer ::				    i,j,n,m_g,n_g		! iteration variables
    character(60) ::			name				! matrix name

    M	= 0.0D0
  
    if(Spar(Ntr)%Nsp.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))  
   
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Spar(Ntr)%Nsp
                call SparCapGeometry(1,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute Spar(Ntr)%lcsc, zeta
                call AxesTranformation(1,zeta,etaj(n_g),x,y)						! compute x, y
                call AxesTranformation(3,zeta,etaj(n_g),x,y,DJ,J1)						! compute DJ, J1
                call SparCapGeometry(2,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L

                call TrialPolynomials(2,k,zeta,etaj(n_g),BdB)
                call TrialPolynomials(3,k,zeta,etaj(n_g),BBd)
                TC11=J1(1,1)*BdB+J1(2,1)*BBd
                TC22=J1(1,2)*BdB+J1(2,2)*BBd
                call MultiplyVectorByTransposedVector(TC11,HH11)
                call MultiplyVectorByVector(TC11,TC22,HH12)
                call MultiplyVectorByVector(TC22,TC11,HH21)
                call MultiplyVectorByTransposedVector(TC22,HH22)

                N11=Spar(Ntr)%D11sc(n)*HH11+Spar(Ntr)%D31sc(n)*HH21+Spar(Ntr)%D13sc(n)*HH12+Spar(Ntr)%D33sc(n)*HH22
                N12=Spar(Ntr)%D12sc(n)*HH12+Spar(Ntr)%D32sc(n)*HH22+Spar(Ntr)%D13sc(n)*HH11+Spar(Ntr)%D33sc(n)*HH21
                N22=Spar(Ntr)%D22sc(n)*HH22+Spar(Ntr)%D32sc(n)*HH12+Spar(Ntr)%D23sc(n)*HH21+Spar(Ntr)%D33sc(n)*HH11
                N24=Spar(Ntr)%D21sc(n)*HH21+Spar(Ntr)%D31sc(n)*HH11+Spar(Ntr)%D23sc(n)*HH22+Spar(Ntr)%D33sc(n)*HH12
                ! Integration constants.
                I1	= Spar(Ntr)%lcsc*((z2U-z1U)+(z2L-z1L))*DJ
                I2	= Spar(Ntr)%lcsc*((z2U*z2U-z1U*z1U)+(z2L*z2L-z1L*z1L))/2.0D0*DJ
                I3	= Spar(Ntr)%lcsc*((z2U*z2U*z2U-z1U*z1U*z1U) + (z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0*DJ
                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j) = M(i,j)+gi(m_g)*gj(n_g)*N11(i,j)*I1
                    ! M12.
                    M(i,j+k*k) = M(i,j+k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I1
                    ! M13.=0
                    ! M14.
                    M(i,j+3*k*k) = M(i,j+3*k*k)+gi(m_g)*gj(n_g)*N11(i,j)*I2
                    ! M15.
                    M(i,j+4*k*k) = M(i,j+4*k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I2
                    ! M21.
                    M(i+k*k,j) = M(j,i+k*k)
                    ! M22.
                    M(i+k*k,j+k*k) = M(i+k*k,j+k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I1
                    ! M23.=0
                    ! M24.
                    M(i+k*k,j+3*k*k) = M(i+k*k,j+3*k*k)+gi(m_g)*gj(n_g)*N24(i,j)*I2
                    ! M25.
                    M(i+k*k,j+4*k*k)	= M(i+k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I2
                    ! M31.=0
                    ! M32.=0
                    ! M33.=0
                    ! M34.=0
                    ! M35.=0
                    ! M41.
                    M(i+3*k*k,j) = M(j,i+3*k*k)
                    ! M42.
                    M(i+3*k*k,j+k*k) = M(j+k*k,i+3*k*k)
                    ! M43.=0
                    ! M44.
                    M(i+3*k*k,j+3*k*k) = M(i+3*k*k,j+3*k*k)+gi(m_g)*gj(n_g)*N11(i,j)*I3
                    ! M45.
                    M(i+3*k*k,j+4*k*k) = M(i+3*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I3
                    ! M51.
                    M(i+4*k*k,j) = M(j,i+4*k*k)
                    ! M52.
                    M(i+4*k*k,j+k*k) = M(j+k*k,i+4*k*k)
                    ! M53.=0
                    ! M54.
                    M(i+4*k*k,j+3*k*k) = M(j+3*k*k,i+4*k*k)
                    !M55
                    M(i+4*k*k,j+4*k*k) = M(i+4*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I3	
                end forall
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)  

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'spar cap stiffness matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if

end subroutine SparCapStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate spar stiffness matrix.
subroutine SparStiffnessMatrix(k,M)

  use GeneralVariables, only: PrintMatrix 

  implicit none

  integer, intent(in)   ::  k					! number of polynomials
  real(8), intent(out)  ::  M(5*k*k,5*k*k)  	! stiffness matrix
  real(8) ::                Msw(5*k*k,5*k*k)	! skin stiffness matrix
  real(8) ::                Msc(5*k*k,5*k*k)	! skin stiffness matrix
  character(60) ::          name				! matrix name

  call SparWebStiffnessMatrix(k,Msw)
  call SparCapStiffnessMatrix(k,Msc)
  M = Msw + Msc 

  ! Write matrix into a file.
  if(PrintMatrix) then
    name	= 'spar stiffness matrix'
    call WriteMatrix(M,5*k*k,5*k*k,name)
  end if

end subroutine SparStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate rib web stiffness matrix.
subroutine RibWebStiffnessMatrix(k,M)
    use PlanformVariables, only: Ntr
    use RibVariablesRoutines, only: Rib,RibGeometry,RibWebGeometry
    use IntegrationVariables, only: mg,ng 
    use GeneralVariables, only: PrintMatrix
    USE maths
    
    implicit none

    integer, intent(in)   ::    k					! number of polynomials
    real(8), intent(out)  ::    M(5*k*k,5*k*k)  	! stiffness matrix
    real(8) ::				    BB(k*k)				! Legendre polynimial multiplication vector
    real(8) ::				    BdB(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k)
    real(8) ::				    HH11(k*k,k*k),HH12(k*k,k*k),HH145(k*k,k*k)
    real(8) ::				    HH21(k*k,k*k),HH22(k*k,k*k),HH245(k*k,k*k),HH4545(k*k,k*k)
    real(8) ::				    N11(k*k,k*k),N12(k*k,k*k),N22(k*k,k*k),N33(k*k,k*k),N34(k*k,k*k)
    real(8) ::				    N35(k*k,k*k),N44(k*k,k*k),N45(k*k,k*k),N55(k*k,k*k)
    real(8) ::				    I1,I2,I3			! integration constants
    real(8) ::				    x,y					! point coordinates
    real(8) ::				    eta					! eta-coordinate
    real(8) ::				    z1,z2				! limits of integration
    integer ::				    i,j,n,m_g,n_g		! iteration variables
    real(8) ::                  DJ,J1(2,2)
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    
    character(60) ::			name				! matrix name
    logical                     error

    M = 0.0D0

    if(Rib(Ntr)%Nrb.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))  

    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Rib(Ntr)%Nrb
                call RibWebGeometry(1,n,zetai(m_g),etaj(n_g),eta,x,y,z1,z2,error)		! compute Rib(Ntr)%tcrw, eta
                call AxesTranformation(1,zetai(m_g),eta,x,y)				! compute x, y
                
                call AxesTranformation(3,zetai(m_g),eta,x,y,DJ,J1)				! compute DJ, J1
                call RibWebGeometry(2,n,zetai(m_g),etaj(n_g),eta,x,y,z1,z2,error)		! compute z1, z2

                call TrialPolynomials(1,k,zetai(m_g),eta,BB) !TC45
                call TrialPolynomials(2,k,zetai(m_g),eta,BdB)
                call TrialPolynomials(3,k,zetai(m_g),eta,BBd)
                TC11=J1(1,1)*BdB+J1(2,1)*BBd
                TC22=J1(1,2)*BdB+J1(2,2)*BBd
                call MultiplyVectorByTransposedVector(TC11,HH11)
                call MultiplyVectorByVector(TC11,TC22,HH12)
                call MultiplyVectorByVector(TC22,TC11,HH21)
                call MultiplyVectorByTransposedVector(TC22,HH22)
                call MultiplyVectorByVector(TC11,BB,HH145)
                call MultiplyVectorByVector(TC22,BB,HH245)
                call MultiplyVectorByTransposedVector(BB,HH4545)
                N11=Rib(Ntr)%D11rw(n)*HH11+Rib(Ntr)%D31rw(n)*HH21+Rib(Ntr)%D13rw(n)*HH12+Rib(Ntr)%D33rw(n)*HH22
                N12=Rib(Ntr)%D12rw(n)*HH12+Rib(Ntr)%D32rw(n)*HH22+Rib(Ntr)%D13rw(n)*HH11+Rib(Ntr)%D33rw(n)*HH21
                N22=Rib(Ntr)%D22rw(n)*HH22+Rib(Ntr)%D32rw(n)*HH12+Rib(Ntr)%D23rw(n)*HH21+Rib(Ntr)%D33rw(n)*HH11
                N33=Rib(Ntr)%D44rw(n)*HH22+Rib(Ntr)%D54rw(n)*HH12+Rib(Ntr)%D45rw(n)*HH21+Rib(Ntr)%D55rw(n)*HH11
                N34=Rib(Ntr)%D45rw(n)*HH245+Rib(Ntr)%D55rw(n)*HH145
                N35=Rib(Ntr)%D44rw(n)*HH245+Rib(Ntr)%D54rw(n)*HH145
                N44=Rib(Ntr)%D55rw(n)*HH4545
                N45=Rib(Ntr)%D54rw(n)*HH4545
                N55=Rib(Ntr)%D44rw(n)*HH4545
                ! Integration constants.		
                I1	= Rib(Ntr)%tcrw*(z2-z1)*DJ
                I2	= Rib(Ntr)%tcrw*(z2*z2-z1*z1)/2.0D0*DJ
                I3	= Rib(Ntr)%tcrw*(z2*z2*z2-z1*z1*z1)/3.0D0*DJ
                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j)				= M(i,j)+gi(m_g)*gj(n_g)*N11(i,j)*I1
                    ! M12.
                    M(i,j+k*k)			= M(i,j+k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I1
                    ! M13.=0
                    ! M14.
                    M(i,j+3*k*k)		= M(i,j+3*k*k)+gi(m_g)*gj(n_g)*N11(i,j)*I2
                    ! M15.
                    M(i,j+4*k*k)		= M(i,j+4*k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I2
                    ! M21.
                    M(i+k*k,j)			= M(j,i+k*k) !!!!!!!!!
                    ! M22.
                    M(i+k*k,j+k*k)		= M(i+k*k,j+k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I1		   
                    ! M23.=0
                    ! M24.
                    M(i+k*k,j+3*k*k)	= M(j,i+4*k*k) !!!!!!
                    ! M25.
                    M(i+k*k,j+4*k*k)	= M(i+k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I2
                    ! M31.=0
                    ! M32.=0
                    ! M33.
                    M(i+2*k*k,j+2*k*k)	= M(i+2*k*k,j+2*k*k)+gi(m_g)*gj(n_g)*N33(i,j)*I1	
                    ! M34.
                    M(i+2*k*k,j+3*k*k)	= M(i+2*k*k,j+3*k*k)+gi(m_g)*gj(n_g)*N34(i,j)*I1
                    ! M35.
                    M(i+2*k*k,j+4*k*k)	= M(i+2*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N35(i,j)*I1
                    ! M41.
                    M(i+3*k*k,j)		= M(j,i+3*k*k) !!!!
                    ! M42.
                    M(i+3*k*k,j+k*k)	= M(j+k*k,i+3*k*k) !!!!!
                    ! M43.
                    M(i+3*k*k,j+2*k*k)	= M(j+2*k*k,i+3*k*k) !!!!!!
                    ! M44.
                    M(i+3*k*k,j+3*k*k)	= M(i+3*k*k,j+3*k*k)+gi(m_g)*gj(n_g)*(N11(i,j)*I3+N44(i,j)*I1)   
                    ! M45.
                    M(i+3*k*k,j+4*k*k)	= M(i+3*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*(N12(i,j)*I3+N45(i,j)*I1)
                    ! M51.
                    M(i+4*k*k,j)		= M(j,i+4*k*k)!!!
                    ! M52.
                    M(i+4*k*k,j+k*k)	= M(j+k*k,i+4*k*k)!!!!
                    ! M53.
                    M(i+4*k*k,j+2*k*k)	= M(j+2*k*k,i+4*k*k)!!!
                    ! M54.
                    M(i+4*k*k,j+3*k*k)	= M(j+3*k*k,i+4*k*k) !!!
                    ! M55.
                    M(i+4*k*k,j+4*k*k)	= M(i+4*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*(N22(i,j)*I3+N55(i,j)*I1)
                end forall          	  
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)  

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'rib web stiffness matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if

end subroutine RibWebStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate rib cap stiffness matrix.
subroutine RibCapStiffnessMatrix(k,M)
    use PlanformVariables, only: Ntr
    use RibVariablesRoutines, only: Rib,RibGeometry,RibCapGeometry
    use IntegrationVariables, only: mg,ng 
    use GeneralVariables, only: PrintMatrix
    USE maths

    implicit none

    integer, intent(in)   ::    k					! number of polynomials
    real(8), intent(out)  ::    M(5*k*k,5*k*k)  	! stiffness matrix
    real(8) ::				    BdB(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k)
    real(8) ::				    HH11(k*k,k*k),HH12(k*k,k*k),HH21(k*k,k*k),HH22(k*k,k*k)
    real(8) ::				    N11(k*k,k*k),N12(k*k,k*k),N22(k*k,k*k),N24(k*k,k*k)
    real(8) ::				    I1,I2,I3			! integration constants
    real(8) ::				    x,y					! point coordinates
    real(8) ::				    eta					! eta-coordinate
    real(8) ::				    z1U,z2U,z1L,z2L		! limits of integration
    integer ::				    i,j,n,m_g,n_g		! iteration variables
    real(8) ::                  DJ,J1(2,2)
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    
    character(60) ::			name				! matrix name


    M = 0.0D0
    if(Rib(Ntr)%Nrb.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))  

    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Rib(Ntr)%Nrb
                call RibCapGeometry(1,n,zetai(m_g),etaj(n_g),eta,x,y,z1U,z2U,z1L,z2L)	! compute Rib(Ntr)%lcrc, eta
                call AxesTranformation(1,zetai(m_g),eta,x,y,DJ,J1)						! compute x, y
                call AxesTranformation(3,zetai(m_g),eta,x,y,DJ,J1)						! compute DJ, J1
                call RibCapGeometry(2,n,zetai(m_g),etaj(n_g),eta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L

                call TrialPolynomials(2,k,zetai(m_g),eta,BdB)
                call TrialPolynomials(3,k,zetai(m_g),eta,BBd)
                TC11=J1(1,1)*BdB+J1(2,1)*BBd
                TC22=J1(1,2)*BdB+J1(2,2)*BBd
   
                call MultiplyVectorByTransposedVector(TC11,HH11)
                call MultiplyVectorByVector(TC11,TC22,HH12)
                call MultiplyVectorByVector(TC22,TC11,HH21)
                call MultiplyVectorByTransposedVector(TC22,HH22)

                N11=Rib(Ntr)%D11rc(n)*HH11+Rib(Ntr)%D31rc(n)*HH21+Rib(Ntr)%D13rc(n)*HH12+Rib(Ntr)%D33rc(n)*HH22
                N12=Rib(Ntr)%D12rc(n)*HH12+Rib(Ntr)%D32rc(n)*HH22+Rib(Ntr)%D13rc(n)*HH11+Rib(Ntr)%D33rc(n)*HH21
                N22=Rib(Ntr)%D22rc(n)*HH22+Rib(Ntr)%D32rc(n)*HH12+Rib(Ntr)%D23rc(n)*HH21+Rib(Ntr)%D33rc(n)*HH11
                N24=Rib(Ntr)%D21rc(n)*HH21+Rib(Ntr)%D31rc(n)*HH11+Rib(Ntr)%D23rc(n)*HH22+Rib(Ntr)%D33rc(n)*HH12

                ! Integration constants.
                I1	= Rib(Ntr)%lcrc*((z2U-z1U)+(z2L-z1L))*DJ
                I2	= Rib(Ntr)%lcrc*((z2U*z2U-z1U*z1U)+(z2L*z2L-z1L*z1L))/2.0D0*DJ
                I3	= Rib(Ntr)%lcrc*((z2U*z2U*z2U-z1U*z1U*z1U)+(z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0*DJ

                !do i=1,k*k
                !do j=1,k*k
                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j)				= M(i,j)+gi(m_g)*gj(n_g)*N11(i,j)*I1
                    ! M12.
                    M(i,j+k*k)			= M(i,j+k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I1
                    ! M13.=0
                    ! M14.
                    M(i,j+3*k*k)		= M(i,j+3*k*k)+gi(m_g)*gj(n_g)*N11(i,j)*I2
                    ! M15.
                    M(i,j+4*k*k)		= M(i,j+4*k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I2
                    ! M21.
                    M(i+k*k,j)			= M(j,i+k*k)
                    ! M22.
                    M(i+k*k,j+k*k)		= M(i+k*k,j+k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I1
                    ! M23.=0
                    ! M24.
                    M(i+k*k,j+3*k*k)	= M(i+k*k,j+3*k*k)+gi(m_g)*gj(n_g)*N24(i,j)*I2
                    ! M25.
                    M(i+k*k,j+4*k*k)	= M(i+k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I2
                    ! M31.=0
                    ! M32.=0
                    ! M33.=0
                    ! M34.=0
                    ! M35.=0
                    ! M41.
                    M(i+3*k*k,j)		= M(j,i+3*k*k)
                    ! M42.
                    M(i+3*k*k,j+k*k)	= M(j+k*k,i+3*k*k)
                    ! M43.=0
                    ! M44.
                    M(i+3*k*k,j+3*k*k)	= M(i+3*k*k,j+3*k*k)+gi(m_g)*gj(n_g)*N11(i,j)*I3
                    ! M45.
                    M(i+3*k*k,j+4*k*k)	= M(i+3*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I3
                    ! M51.
                    M(i+4*k*k,j)		= M(j,i+4*k*k)
                    ! M52.
                    M(i+4*k*k,j+k*k)	= M(j+k*k,i+4*k*k)
                    ! M53.=0
                    ! M54.
                    M(i+4*k*k,j+3*k*k)	= M(j+3*k*k,i+4*k*k)
                    !M55
                    M(i+4*k*k,j+4*k*k)	= M(i+4*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I3
                end forall  
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)  

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'rib cap stiffness matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if

end subroutine RibCapStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate rib stiffness matrix.
subroutine RibStiffnessMatrix(k,M)

  use GeneralVariables, only: PrintMatrix

  implicit none

  integer, intent(in)   ::  k					! number of polynomials
  real(8), intent(out)  ::  M(5*k*k,5*k*k)  	! stiffness matrix
  real(8) ::                Mrw(5*k*k,5*k*k)	! skin stiffness matrix
  real(8) ::                Mrc(5*k*k,5*k*k)	! skin stiffness matrix
  character(60) ::          name				! matrix name

  call RibWebStiffnessMatrix(k,Mrw)
  call RibCapStiffnessMatrix(k,Mrc)
  M = Mrw +Mrc 

  ! Write matrix into a file.
  if(PrintMatrix) then
    name	= 'rib stiffness matrix'
    call WriteMatrix(M,5*k*k,5*k*k,name)
  end if

end subroutine RibStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate stringer stiffness matrix.
subroutine StringerStiffnessMatrix(k,M)
    use PlanformVariables, only: Ntr
    use StringerVariablesRoutines, only: Stringer,StringerGeometry1,StringerGeometry2
    use IntegrationVariables, only: Mg,Ng
    use GeneralVariables, only: PrintMatrix
    USE maths

    implicit none

    integer, intent(in)   ::    k					! number of polynomials
    real(8), intent(out)  ::    M(5*k*k,5*k*k)  	! stiffness matrix
    real(8) ::				    BdB(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)			! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k)
    real(8) ::				    HH11(k*k,k*k),HH12(k*k,k*k),HH21(k*k,k*k),HH22(k*k,k*k)
    real(8) ::				    N11(k*k,k*k),N12(k*k,k*k),N22(k*k,k*k),N24(k*k,k*k)
    real(8) ::				    I1,I2,I3			! integration constants
    real(8) ::				    x,y					! point coordinates
    real(8) ::				    zeta				! zeta-coordinate
    real(8) ::				    z1U,z2U,z1L,z2L		! limits of integration
    integer ::				    i,j,n,m_g,n_g		! iteration variables
    real(8) ::                  DJ,J1(2,2)
    real(8),allocatable ::      gi(:),gj(:),zetai(:),etaj(:)    
    character(60) ::			name				! matrix name

    M = 0.0D0
    if(Stringer(Ntr)%Nst.EQ.0) return
    allocate(gi(Mg),gj(Ng),zetai(Mg),etaj(Ng))  

    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    call GaussQuadratureWeightsSamplingPoints(gj,etaj)
    do m_g=1,Mg
        do n_g=1,Ng
            do n=1,Stringer(Ntr)%Nst
                call StringerGeometry2(1,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute Stringer(Ntr)%lcst, zeta
                call AxesTranformation(1,zeta,etaj(n_g),x,y)							    ! compute x, y
                call AxesTranformation(3,zeta,etaj(n_g),x,y,DJ,J1)							! compute DJ, J1
                call StringerGeometry2(2,n,zetai(m_g),etaj(n_g),zeta,x,y,z1U,z2U,z1L,z2L)	! compute z1U, z2U, z1L, z2L
                call TrialPolynomials(2,k,zeta,etaj(n_g),BdB)
                call TrialPolynomials(3,k,zeta,etaj(n_g),BBd)
                TC11=J1(1,1)*BdB+J1(2,1)*BBd
                TC22=J1(1,2)*BdB+J1(2,2)*BBd
                call MultiplyVectorByTransposedVector(TC11,HH11)
                call MultiplyVectorByVector(TC11,TC22,HH12)
                call MultiplyVectorByVector(TC22,TC11,HH21)
                call MultiplyVectorByTransposedVector(TC22,HH22)
                N11 = Stringer(Ntr)%D11st(n)*HH11+Stringer(Ntr)%D31st(n)*HH21+Stringer(Ntr)%D13st(n)*HH12+Stringer(Ntr)%D33st(n)*HH22
                N12 = Stringer(Ntr)%D12st(n)*HH12+Stringer(Ntr)%D32st(n)*HH22+Stringer(Ntr)%D13st(n)*HH11+Stringer(Ntr)%D33st(n)*HH21
                N22 = Stringer(Ntr)%D22st(n)*HH22+Stringer(Ntr)%D32st(n)*HH12+Stringer(Ntr)%D23st(n)*HH21+Stringer(Ntr)%D33st(n)*HH11
                N24 = Stringer(Ntr)%D21st(n)*HH21+Stringer(Ntr)%D31st(n)*HH11+Stringer(Ntr)%D23st(n)*HH22+Stringer(Ntr)%D33st(n)*HH12
                ! Integration constants.
                I1	= Stringer(Ntr)%lcst*((z2U-z1U)+(z2L-z1L))*DJ
                I2	= Stringer(Ntr)%lcst*((z2U*z2U-z1U*z1U)+(z2L*z2L-z1L*z1L))/2.0D0*DJ
                I3	= Stringer(Ntr)%lcst*((z2U*z2U*z2U-z1U*z1U*z1U)+(z2L*z2L*z2L-z1L*z1L*z1L))/3.0D0*DJ
                forall(i=1:k*k,j=1:k*k)
                    ! M11.
                    M(i,j) = M(i,j)+gi(m_g)*gj(n_g)*N11(i,j)*I1
                    ! M12.
                    M(i,j+k*k) = M(i,j+k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I1
                    ! M13.=0
                    ! M14.
                    M(i,j+3*k*k) = M(i,j+3*k*k)+gi(m_g)*gj(n_g)*N11(i,j)*I2
                    ! M15.
                    M(i,j+4*k*k) = M(i,j+4*k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I2
                    ! M21.
                    M(i+k*k,j) = M(j,i+k*k)
                    ! M22.
                    M(i+k*k,j+k*k) = M(i+k*k,j+k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I1
                    ! M23.=0
                    ! M24.
                    M(i+k*k,j+3*k*k) = M(i+k*k,j+3*k*k)+gi(m_g)*gj(n_g)*N24(i,j)*I2
                    ! M25.
                    M(i+k*k,j+4*k*k)	= M(i+k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I2
                    ! M31.=0
                    ! M32.=0
                    ! M33.=0
                    ! M34.=0
                    ! M35.=0
                    ! M41.
                    M(i+3*k*k,j) = M(j,i+3*k*k)
                    ! M42.
                    M(i+3*k*k,j+k*k) = M(j+k*k,i+3*k*k)
                    ! M43.=0
                    ! M44.
                    M(i+3*k*k,j+3*k*k) = M(i+3*k*k,j+3*k*k)+gi(m_g)*gj(n_g)*N11(i,j)*I3
                    ! M45.
                    M(i+3*k*k,j+4*k*k) = M(i+3*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N12(i,j)*I3
                    ! M51.
                    M(i+4*k*k,j) = M(j,i+4*k*k)
                    ! M52.
                    M(i+4*k*k,j+k*k) = M(j+k*k,i+4*k*k)
                    ! M53.=0
                    ! M54.
                    M(i+4*k*k,j+3*k*k) = M(j+3*k*k,i+4*k*k)
                    !M55
                    M(i+4*k*k,j+4*k*k) = M(i+4*k*k,j+4*k*k)+gi(m_g)*gj(n_g)*N22(i,j)*I3	                
                end forall
            end do
        end do
    end do
    deallocate(gi,gj,zetai,etaj)  

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'stringer stiffness matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if

end subroutine StringerStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate spring stiffness matrix.
subroutine SpringStiffnessMatrix(k,M)

    use PlanformVariables, only: ntr,croot,yroot,xLEroot,lamdapl,mcpl,b2 
    use SpringVariables, only: Spring
    use GeneralVariables, only: PrintMatrix
    USE maths, only: TrialPolynomials, MultiplyVectorByTransposedVector

    implicit none

    integer, intent(in)   ::    k					! number of polynomials
    real(8), intent(out)  ::    M(5*k*k,5*k*k)  	! stiffness matrix
    real(8) ::				    BB(k*k)				! Legendre polynimial multiplication vector
    real(8) ::				    HH(k*k,k*k)
    real(8) ::				    zeta,eta			! point coordinates
    integer ::				    i,j,n				! iteration variables
    character(60) ::			name				! matrix name

    M = 0.0D0
    if(Spring(Ntr)%Nspr.EQ.0) return
    do n=1,Spring(Ntr)%Nspr
        zeta = -1.0D0+2.0D0*(Spring(Ntr)%xspr(n)-xLEroot(Ntr)-(Spring(Ntr)%yspr(n)-yroot(Ntr))*dtan(Lamdapl(Ntr)))/(croot(Ntr)+mcpl(Ntr)*(Spring(Ntr)%yspr(n)-yroot(Ntr)))
        eta	 = -1.0D0+2.0D0*(Spring(Ntr)%yspr(n)-yroot(Ntr))/b2(Ntr)
        call TrialPolynomials(1,k,zeta,eta,BB)
        call MultiplyVectorByTransposedVector(BB,HH)
        forall(i=1:k*k,j=1:k*k)
            ! M11.
            M(i,j)	= M(i,j)+Spring(Ntr)%ku(n)*HH(i,j)
            ! M22.
            M(i+k*k,j+k*k) = M(i+k*k,j+k*k)+Spring(Ntr)%kv(n)*HH(i,j)
            ! M33.
            M(i+2*k*k,j+2*k*k) = M(i+2*k*k,j+2*k*k)+Spring(Ntr)%kw(n)*HH(i,j)
            ! M44.
            M(i+3*k*k,j+3*k*k) = M(i+3*k*k,j+3*k*k)+Spring(Ntr)%krx(n)*HH(i,j)        
            ! M55.
            M(i+4*k*k,j+4*k*k) = M(i+4*k*k,j+4*k*k)+Spring(Ntr)%kry(n)*HH(i,j)        
        end forall
    end do

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'spring stiffness matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if

end subroutine SpringStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate boundary stiffness matrix.
subroutine RootBoundaryStiffnessMatrix(k,M)

    use PlanformVariables, only: Ntr,xLEroot,cRoot,Boundary 
    use IntegrationVariables, only: Mg 
    use GeneralVariables, only: PrintMatrix
    USE maths, only: GaussQuadratureWeightsSamplingPoints,TrialPolynomials, MultiplyVectorByTransposedVector

    implicit none

    real(8), parameter ::	kgen= 1.0D9
    real(8), parameter ::	ku	= kgen
    real(8), parameter ::	kv	= kgen
    real(8), parameter ::	kw	= kgen
    real(8), parameter ::	krx	= kgen
    real(8), parameter ::	kry	= kgen

    integer, intent(in)   ::    k					! number of polynomials
    real(8), intent(out)  ::    M(5*k*k,5*k*k)  	! stiffness matrix
    real(8) ::				    BB(k*k)				! Legendre polynimial multiplication vector
    real(8) ::				    HH(5,5*k*k)
    real(8) ::				    Kspring(5,5)
    real(8) ::				    x1,x2				! limits of integration
    real(8) ::				    eta					! point coordinates
    real(8),allocatable ::      gi(:),zetai(:)   
    integer ::				    i,m_g				! iteration variables
    character(60) ::			name				! matrix name

    M = 0.0D0
    HH = 0.0D0
    allocate(gi(Mg),zetai(Mg))  
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    eta	= -1.0D0
    Ntr = 1
    x1	= xLEroot(Ntr)
    x2	= x1+croot(Ntr)
    Kspring = 0.0D0
    !Calculate Spring stiffness matrix
    Kspring(1,1) = ku
    Kspring(2,2) = kv
    Kspring(3,3) = kw
    if(boundary) then
        Kspring(4,4) = krx
        Kspring(5,5) = kry
    end if
    Kspring = Kspring*(x2-x1)/2.0D0
    do m_g=1,Mg
        call TrialPolynomials(1,k,zetai(m_g),eta,BB)
        forall(i=1:5) HH(i,k*k*(i-1)+1:k*k*i) = gi(m_g)*BB
        M = M + MatMul(MatMul(TRANSPOSE(HH),Kspring),HH)
    end do !Mg  
   
    deallocate(gi,zetai)
    ! Write matrix into a file
    if(PrintMatrix) then
        name	= 'boundary conditions stiffness matrix'
        call WriteMatrix(M,5*k*k,5*k*k,name)
    end if

end subroutine RootBoundaryStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate planform stiffness matrix.
subroutine PlanformStiffnessMatrix(k)
    use PlanformVariables, only: Ntr
    use DeformationVariables, only: Deformation
    use GeneralVariables, only: PrintMatrix

    implicit none

    integer, intent(in)   ::    k

    real(8), allocatable  ::    Msk(:,:)	! skin stiffness matrix
    real(8), allocatable  ::    Msp(:,:)	! spar stiffness matrix
    real(8), allocatable  ::    Mrb(:,:)	! rib stiffness matrix
    real(8), allocatable  ::    Mst(:,:)	! stringer stiffness matrix
    real(8), allocatable  ::    Mspr(:,:)	! spring stiffness matrix
    character(60) ::            name        ! matrix name

    allocate(Msk(5*k*k,5*k*k),Msp(5*k*k,5*k*k),Mrb(5*k*k,5*k*k),Mst(5*k*k,5*k*k),Mspr(5*k*k,5*k*k))
    write(*,'(1X,A,1X,I1,1X,A)') 'Computing trapezoid', Ntr, 'stiffness matrix...'
    !$OMP PARALLEL SECTIONS
    !$OMP SECTION
        write(*,*) '    - skin'
        call SkinStiffnessMatrix(k,Msk)
    !$OMP SECTION
        write(*,*) '    - spar'
        call SparStiffnessMatrix(k,Msp)
    !$OMP SECTION
        write(*,*) '    - rib'
        call RibStiffnessMatrix(k,Mrb)
    !$OMP SECTION
        write(*,*) '    - stringer'
        call StringerStiffnessMatrix(k,Mst)
    !$OMP SECTION
        write(*,*) '    - spring'
        call SpringStiffnessMatrix(k,Mspr)
    !$OMP END PARALLEL SECTIONS

    Deformation(Ntr)%KM = Msk + Msp + Mrb + Mst + Mspr  

    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'planform stiffness matrix'
        call WriteMatrix(Deformation(Ntr)%KM,5*k*k,5*k*k,name)
    end if
    deallocate(Msk,Msp,Mrb,Mst,Mspr)

end subroutine PlanformStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate boundary stiffness matrix.
subroutine JointStiffnessMatrix(k,M1,M2,M3)

    use PlanformVariables, only: Ntr,xLEroot,cRoot,TransMatrix 
    use IntegrationVariables, only: Mg 
    use GeneralVariables, only: PrintMatrix
    USE maths, only: GaussQuadratureWeightsSamplingPoints,TrialPolynomials, MultiplyVectorByTransposedVector,MultiplyVectorByTransposedVectorC

    implicit none

    real(8), parameter ::	kgen= 1.0D12
    real(8), parameter ::	ku	= kgen
    real(8), parameter ::	kv	= kgen
    real(8), parameter ::	kw	= kgen
    real(8), parameter ::	krx	= kgen !/ 35000.0D0 !/ 125000.0D0
    real(8), parameter ::	kry	= kgen !/ 35000.0D0 !/ 125000.0D0
    !Input vars
    integer, intent(in)   ::    k					! number of polynomials
    !Output vars
    real(8), intent(out)  ::    M1(5*k*k,5*k*k),M2(5*k*k,5*k*k),M3(5*k*k,5*k*k)  	! stiffness matrix
    !Local vars
    real(8) ::				    BB(k*k)				        ! Legendre polynomial multiplication vector
    real(8) ::				    HH1(5,5*k*k),HH2(5,5*k*k),M3_aux(5,5)
    real(8) ::				    Kspring(5,5)
    real(8) ::				    RR(5,5)                     ! Global rotation matrix
    real(8) ::				    x1,x2				        ! limits of integration
    real(8) ::				    eta					        ! point coordinates
    real(8),allocatable ::      gi(:),zetai(:)   
    integer ::				    i,m_g				        ! iteration variables
    character(60) ::			name				        ! matrix name

    M1 = 0.0D0
    M2 = 0.0D0
    M3 = 0.0D0
    HH1 = 0.0D0
    HH2 = 0.0D0
    Kspring = 0.0D0
    !Calculate Spring stiffness matrix
    x1 = xLEroot(Ntr) 
    x2 = x1+croot(Ntr)
    Kspring(1,1) = ku
    Kspring(2,2) = kv
    Kspring(3,3) = kw
    Kspring(4,4) = krx
    Kspring(5,5) = kry
    Kspring = Kspring*(x2-x1) /2.0D0
    !Build rotation matrix
    RR = 0.0D0
    forall(i=1:5) 
        RR(i,i) = 1.0D0 !TransMatrix(Ntr)%Rot(2,1)
    end forall
    !Calculate centre matrix in K12 term
    M3_aux = MATMUL(Kspring,TRANSPOSE(RR)) + MATMUL(TRANSPOSE(Kspring),RR)
  
    allocate(gi(Mg),zetai(Mg))  
    call GaussQuadratureWeightsSamplingPoints(gi,zetai)
    !Calculate 
    do m_g=1,Mg
        eta	= 1.0D0 !1st trapezoid joint location 
        call TrialPolynomials(1,k,zetai(m_g),eta,BB)
        forall(i=1:5) HH1(i,k*k*(i-1)+1:k*k*i) = gi(m_g)*BB
        eta	= -1.0D0 !2nd trapezoid joint location 
        call TrialPolynomials(1,k,zetai(m_g),eta,BB)
        forall(i=1:5) HH2(i,k*k*(i-1)+1:k*k*i) = gi(m_g)*BB
        !Compute [H1^T][K_JT][H1] term 
        M1 = M1 + MatMul(MatMul(TRANSPOSE(HH1),Kspring),HH1) !Only calculate on 
        !Compute [H2^T][R^T][K_JT][R][H2] term
        M2 = M2 + MATMUL(MATMUL(MATMUL(MATMUL(TRANSPOSE(HH2),TRANSPOSE(RR)),Kspring),RR),HH2)
        !Calculate -1/2 * [H1^T] ([K_JT][R^T] + [K_JT][R]) [H2] term
        M3 = M3 - 0.5D0 * MATMUL(MATMUL(TRANSPOSE(HH1),M3_aux),HH2)   
    end do !Mg
    deallocate(gi,zetai)
    ! Write matrix into a file
    if(PrintMatrix) then
        name	= 'joint stiffness matrix M1'
        call WriteMatrix(M1,5*k*k,5*k*k,name)
        name	= 'joint stiffness matrix M2'
        call WriteMatrix(M2,5*k*k,5*k*k,name)
        name	= 'joint stiffness matrix M3'
        call WriteMatrix(M3,5*k*k,5*k*k,name)
    end if

end subroutine JointStiffnessMatrix

!--------------------------------------------------------------------------------------------
! Routine to calculate overall stiffness matrix.
subroutine OverallStiffnessMatrix(k)
    use PlanformVariables, only: Npl,Ntr
    use DeformationVariables, only: Deformation,KM_Global
    use GeneralVariables, only: PrintMatrix

    implicit none
    !Input vars
    integer, intent(in)   ::    k
    !Local vars
    character(60) ::			name        ! matrix name
    real(8), allocatable  ::    Mbc(:,:)	! root boundary conditions stiffness matrix
    real(8), allocatable  ::    KM1(:,:),KM2(:,:),KM3(:,:)	! individual stiffness matrix groups

    write(*,'(1X,A)') 'Computing overall stiffness matrix...'

    allocate(Mbc(5*k*k,5*k*k)) !root boundary condition
     
    !Apply root boundary condion (only root section of first planform)
    call RootBoundaryStiffnessMatrix(k,Mbc)
    if(Npl > 1) then
        allocate(KM1(5*k*k,5*k*k),KM2(5*k*k,5*k*k),KM3(5*k*k,5*k*k))
        KM_Global = 0.0D0
        do Ntr=2, Npl
            call JointStiffnessMatrix(k,KM1,KM2,KM3)      
            !Build complete stifness matrix
            if(Ntr == 2) then
                KM_Global(1:5*k*k,1:5*k*k) = Deformation(1)%KM + KM1 + Mbc
            else
                KM_Global((Ntr-2)*5*k*k+1:(Ntr-1)*5*k*k,(Ntr-2)*5*k*k+1:(Ntr-1)*5*k*k) = KM_Global((Ntr-2)*5*k*k+1:(Ntr-1)*5*k*k,(Ntr-2)*5*k*k+1:(Ntr-1)*5*k*k) + KM1
            end if
            KM_Global((Ntr-1)*5*k*k+1:Ntr*5*k*k,(Ntr-1)*5*k*k+1:Ntr*5*k*k) = Deformation(Ntr)%KM + KM2 !Second Diagonal matrix
            KM_Global((Ntr-2)*5*k*k+1:(Ntr-1)*5*k*k,(Ntr-1)*5*k*k+1:Ntr*5*k*k) = KM3 !Upper Cross term 
            KM_Global((Ntr-1)*5*k*k+1:Ntr*5*k*k,(Ntr-2)*5*k*k+1:(Ntr-1)*5*k*k) = Transpose(KM3) !Lower Cross term 
            if(Ntr == Npl) deallocate(KM1,KM2,KM3) !Last planform: deallocate vars
        end do
    else 
        Ntr=1 !only one planform
        KM_Global = Deformation(Ntr)%KM + Mbc  
    end if
    deallocate(Mbc)
    ! Write matrix into a file.
    if(PrintMatrix) then
        name	= 'overall planform stiffness matrix'
        call WriteMatrix(KM_Global,size(KM_Global,DIM=1),size(KM_Global,DIM=2),name)
    end if

end subroutine OverallStiffnessMatrix