!--------------------------------------------------------------------------------------------
! Routine to calculate Skin Stresses
subroutine SkinStresses(k)
    
    use SkinVariablesRoutines
    use DeformationVariables, only: Deformation
    use PlanformVariables, only: Ntr
    use StressStrainVariables
    USE FilePathModule, only: output_dir,skin_stress_unit,ierror,msg,abort_program
    USE maths, only: AxesTranformation,TrialPolynomials,AddTwoVectors

    implicit none

    integer, intent(in) ::              k
    real(8) ::                          BdB(k*k)						! Legendre polynimial multiplication vector
    real(8) ::                          BBd(k*k)						! Legendre polynimial multiplication vector
    real(8) ::                          TC11(k*k),TC22(k*k)
    real(8) ::                          TC45(k*k)
    real(8) ::                          x,y,z,z1,z2,zetai,zeta,eta            ! point coordinates
    real(8) ::                          zetaskF,zetaskA
    real(8) ::                          u, v, w                         ! point displacements - needed only in coupled mode
    integer ::                          l,n,xi,yi,xtotal,ytotal         ! iteration variables	  
    real(8) ::                          du0dx,du0dy,dv0dx,dv0dy,dw0dx,dw0dy,dphixdx,dphixdy,dphiydx,dphiydy,phixx,phiyy
    real(8) ::                          DJ,J1(2,2),slope

    !xtotal=int(croot(Npl)*30.0D0) !40 !numero de pontos em x
    !ytotal=int(ytip(Npl)*11.0D0) !50 !numero de pontos em y
    xtotal=int(1.0D0/0.025D0)
    ytotal=int(1.0D0/0.025D0) !50 !numero de pontos em y
    !1 ponto em z
    if(Ntr == 1) then !first access
        open(unit=skin_stress_unit,file=trim(adjustl(output_dir))//"skinstresses.dat",STATUS='UNKNOWN',ACTION='Write',IOSTAT=ierror,IOMSG=msg)
        write(skin_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "SKIN STRESSES"'
        write(skin_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES =x,y,z,vonMises_Stress,stressxx,stressyy,stressxy,stressyz,stresszx,strainxx,strainyy,strainxy,strainyz,strainzx'
    else
        open(unit=skin_stress_unit,file=trim(adjustl(output_dir))//"skinstresses.dat",ACTION='Write',ACCESS = 'APPEND',IOSTAT=ierror,IOMSG=msg)
    end if


    if(Skin(Ntr)%Nsk.EQ.0) then
        close(skin_stress_unit,IOSTAT=ierror,IOMSG=msg)
        return
    end if

    sigmavonMisesskmax = 0.0D0
    sigmavonMisessk = 0.0D0

    do n=1,Skin(Ntr)%Nsk
        write(skin_stress_unit,'(A,I2.2,A,I2.2,A,I4.4,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="Planform ',Ntr,' Skin Panel ',n,'" I=',xtotal,' J=',ytotal,' K=',1
        do yi=0,ytotal-1
            eta	= SkinEtaTecplot(ytotal-1,yi)
            call SkinGeometry1(n,eta,zetaskF,zetaskA)		    ! compute zetaskF, zetaskA   
            do xi=0,xtotal-1  
                zetai = SkinZetaTecplot(1,xtotal-1,xi)
                zeta = 0.5D0*(1.0D0-zetai)*zetaskF+0.5D0*(1.0D0+zetai)*zetaskA
                call AxesTranformation(1,zeta,eta,x,y)
                call AxesTranformation(3,zeta,eta,x,y,DJ,J1)
                call TrialPolynomials(1,k,zeta,eta,TC45)
                call TrialPolynomials(2,k,zeta,eta,BdB)
                call TrialPolynomials(3,k,zeta,eta,BBd)  				
                call AddTwoVectors(J1(1,1),BdB,J1(2,1),BBd,TC11)
                call AddTwoVectors(J1(1,2),BdB,J1(2,2),BBd,TC22)
                call SkinGeometry2(2,n,x,y,z,z1,z2,slope)

                epsilonsk=0.0D0
                epsilonsk1=0.0D0
                sigmask=0.0D0

                du0dx	= 0.000D0
                du0dy	= 0.000D0
                dv0dx	= 0.000D0
                dv0dy	= 0.000D0
                dw0dx	= 0.000D0
                dw0dy	= 0.000D0
                dphixdx	= 0.000D0
                dphixdy	= 0.000D0
                dphiydx	= 0.000D0
                dphiydy	= 0.000D0
                phixx	= 0.000D0
                phiyy	= 0.000D0

                do l=1,k*k
                    du0dx	= du0dx+Deformation(Ntr)%q(l)*TC11(l)
                    du0dy	= du0dy+Deformation(Ntr)%q(l)*TC22(l)
                    dv0dx	= dv0dx+Deformation(Ntr)%q(l+k*k)*TC11(l)
                    dv0dy	= dv0dy+Deformation(Ntr)%q(l+k*k)*TC22(l)
                    dw0dx	= dw0dx+Deformation(Ntr)%q(l+2*k*k)*TC11(l)
                    dw0dy	= dw0dy+Deformation(Ntr)%q(l+2*k*k)*TC22(l)
                    dphixdx	= dphixdx+Deformation(Ntr)%q(l+3*k*k)*TC11(l)
                    dphixdy	= dphixdy+Deformation(Ntr)%q(l+3*k*k)*TC22(l)
                    dphiydx	= dphiydx+Deformation(Ntr)%q(l+4*k*k)*TC11(l)
                    dphiydy	= dphiydy+Deformation(Ntr)%q(l+4*k*k)*TC22(l)
                    phixx	= phixx+Deformation(Ntr)%q(l+3*k*k)*TC45(l)
                    phiyy	= phiyy+Deformation(Ntr)%q(l+4*k*k)*TC45(l)
                end do

                epsilonsk(1)= du0dx+z*dphixdx													!epsilon xx global
                epsilonsk(2)= dv0dy+z*dphiydy													!epsilon yy global
                epsilonsk(3)= 1.0D0/2.0D0*(du0dy+dv0dx+z*(dphixdy+dphiydx))							!epsilon xy global
                epsilonsk(4)= 1.0D0/2.0D0*(phiyy+dw0dy)												!epsilon yz global
                epsilonsk(5)= 1.0D0/2.0D0*(phixx+dw0dx)												!epsilon zx global

                epsilonsk1(1)= dcos(datan(-slope))*dcos(datan(-slope))*epsilonsk(1)-dcos(datan(-slope))*dsin(datan(-slope))*epsilonsk(5)	!epsilon xx local
                epsilonsk1(2)= epsilonsk(2)																								!epsilon yy local
                epsilonsk1(3)= dcos(datan(-slope))*epsilonsk(3)-dsin(datan(-slope))*epsilonsk(4)											!epsilon xy local
                epsilonsk1(4)= dsin(datan(-slope))*epsilonsk(3)+dcos(datan(-slope))*epsilonsk(4)											!epsilon yz local
                epsilonsk1(5)= 2.0D0*dsin(datan(-slope))*dcos(datan(-slope))*epsilonsk(1)+dcos(2.0D0*datan(-slope))*epsilonsk(5)				!epsilon zx local

                sigmask(1)=(Skin(Ntr)%sk11(n)*epsilonsk1(1)+Skin(Ntr)%sk12(n)*epsilonsk1(2))/1000000.0D0			!sigma xx local
                sigmask(2)=(Skin(Ntr)%sk21(n)*epsilonsk1(1)+Skin(Ntr)%sk22(n)*epsilonsk1(2))/1000000.0D0			!sigma yy local
                sigmask(3)=(Skin(Ntr)%sk44(n)*epsilonsk1(3))/1000000.0D0								!sigma xy local
                sigmask(4)=(Skin(Ntr)%sk55(n)*epsilonsk1(4))/1000000.0D0								!sigma yz local
                sigmask(5)=(Skin(Ntr)%sk66(n)*epsilonsk1(5))/1000000.0D0								!sigma zx local	

                sigmavonMisessk=dsqrt(((sigmask(1)-sigmask(2))*(sigmask(1)-sigmask(2))+sigmask(2)*sigmask(2)+(-sigmask(1))*(-sigmask(1)) &
                & +6.0D0*(sigmask(3)*sigmask(3)+sigmask(4)*sigmask(4)+sigmask(5)*sigmask(5)))/2.0D0)

                if (sigmavonMisessk.GT.sigmavonMisesskmax) sigmavonMisesskmax=sigmavonMisessk

                write(skin_stress_unit,'(F10.5,6F14.6,7F20.10)',IOSTAT=ierror,IOMSG=msg) x,y,z,sigmavonMisessk,sigmask(1),sigmask(2),sigmask(3),sigmask(4),sigmask(5),epsilonsk1(1),epsilonsk1(2),epsilonsk1(3),epsilonsk1(4),epsilonsk1(5)
            end do
        end do
    end do

    close(skin_stress_unit,IOSTAT=ierror,IOMSG=msg)
    IF(ierror.NE.0) Call abort_program()
   
end subroutine SkinStresses

!--------------------------------------------------------------------------------------------
! Routine to calculate Spar Web Stresses
subroutine SparWebStresses(k)
    
    use PlanformVariables, only: Ntr
    use SparVariablesRoutines, only: Spar,SparGeometry,SparWebGeometry 
    use SkinVariablesRoutines, only: SkinEtaTecplot
    use DeformationVariables, only: Deformation
    use StressStrainVariables
    USE FilePathModule, only: output_dir,sparweb_stress_unit,ierror,msg,abort_program
    USE maths, only: AxesTranformation, TrialPolynomials,AddTwoVectors
    
    implicit none

    integer, intent(in)   ::    k
    real(8) ::				    E22																	! stiffness element
    real(8) ::				    BdB(k*k)															! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)															! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k),TC45(k*k)							
    real(8) ::				    x,y,z,z1sw,z2sw,zeta,eta							! point coordinates				
    integer ::				    l,n,yi,zi,ytotal,ztotal											! iteration variables	  
    real(8) ::				    du0dx,du0dy,dv0dx,dv0dy,dw0dx,dw0dy,dphixdx,dphixdy,dphiydx,dphiydy,phixx,phiyy
    real(8) ::                  DJ,J1(2,2)

    ytotal=45				!numero de pontos em y
    ztotal=10				!numero de pontos em z

    if(Ntr == 1) then !first access
        open(unit=sparweb_stress_unit,file=trim(adjustl(output_dir))//"sparwebstresses.dat",STATUS='UNKNOWN',ACTION='Write',IOSTAT=ierror,IOMSG=msg)
        write(sparweb_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "SPARWEB STRESSES"'
        write(sparweb_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES =x,y,z,vonMises_Stress,stressxx,stressyy,stressxy,stressyz,stresszx,strainxx,strainyy,strainxy,strainyz,strainzx'
    else
        open(unit=sparweb_stress_unit,file=trim(adjustl(output_dir))//"sparwebstresses.dat",ACTION='Write',ACCESS = 'APPEND',IOSTAT=ierror,IOMSG=msg)
    end if

    if(Spar(Ntr)%Nsp.EQ.0) then
        close(sparweb_stress_unit,IOSTAT=ierror,IOMSG=msg)
        return
    end if

    sigmavonMisesswmax = 0.0D0
    sigmavonMisessw = 0.0D0

    do n=1,Spar(Ntr)%Nsp	
        write(sparweb_stress_unit,'(A,I2.2,A,I2.2,A,I4.4,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="Planform ',Ntr,' SparWeb ',n,'" I=',1,' J=',ytotal,' K=',ztotal
        do zi=0,ztotal-1
            do yi=0,ytotal-1
                !eta	 = -1.0D0+2.0D0*float(yi)/float(ytotal-1)
                eta = SkinEtaTecplot(ytotal-1,yi)
                zeta = -1.0D0+2.0D0*(Spar(Ntr)%ksp1(n)-Spar(Ntr)%ksp2(n)*eta)

                call AxesTranformation(1,zeta,eta,x,y)
                call AxesTranformation(3,zeta,eta,x,y,DJ,J1)
                call SparWebGeometry(2,n,zeta,eta,zeta,x,y,z1sw,z2sw)
                call TrialPolynomials(2,k,zeta,eta,BdB)
                call TrialPolynomials(3,k,zeta,eta,BBd)
                call AddTwoVectors(J1(1,1),BdB,J1(2,1),BBd,TC11)
                call AddTwoVectors(J1(1,2),BdB,J1(2,2),BBd,TC22)
                call TrialPolynomials(1,k,zeta,eta,TC45)

                z= z2sw-(zi)*(z2sw-z1sw)/(ztotal-1)		
                epsilonsw=0.0D0
                sigmasw=0.0D0

                du0dx	= 0.000D0
                du0dy	= 0.000D0
                dv0dx	= 0.000D0
                dv0dy	= 0.000D0
                dw0dx	= 0.000D0
                dw0dy	= 0.000D0
                dphixdx	= 0.000D0
                dphixdy	= 0.000D0
                dphiydx	= 0.000D0
                dphiydy	= 0.000D0
                phixx	= 0.000D0
                phiyy	= 0.000D0

                do l=1,k*k
                    du0dx	= du0dx+Deformation(Ntr)%q(l)*TC11(l)
                    du0dy	= du0dy+Deformation(Ntr)%q(l)*TC22(l)
                    dv0dx	= dv0dx+Deformation(Ntr)%q(l+k*k)*TC11(l)
                    dv0dy	= dv0dy+Deformation(Ntr)%q(l+k*k)*TC22(l)
                    dw0dx	= dw0dx+Deformation(Ntr)%q(l+2*k*k)*TC11(l)
                    dw0dy	= dw0dy+Deformation(Ntr)%q(l+2*k*k)*TC22(l)
                    dphixdx	= dphixdx+Deformation(Ntr)%q(l+3*k*k)*TC11(l)
                    dphixdy	= dphixdy+Deformation(Ntr)%q(l+3*k*k)*TC22(l)
                    dphiydx	= dphiydx+Deformation(Ntr)%q(l+4*k*k)*TC11(l)
                    dphiydy	= dphiydy+Deformation(Ntr)%q(l+4*k*k)*TC22(l)
                    phixx	= phixx+Deformation(Ntr)%q(l+3*k*k)*TC45(l)
                    phiyy	= phiyy+Deformation(Ntr)%q(l+4*k*k)*TC45(l)
                end do

                epsilonsw(1)= du0dx+z*dphixdx														!epsilon xx global
                epsilonsw(2)= dv0dy+z*dphiydy														!epsilon yy global
                epsilonsw(3)= 1.0D0/2.0D0*(du0dy+dv0dx+z*(dphixdy+dphiydx))								!epsilon xy global
                epsilonsw(4)= 1.0D0/2.0D0*(phiyy+dw0dy)													!epsilon zx global
                epsilonsw(5)= 1.0D0/2.0D0*(phixx+dw0dx)													!epsilon yz global

                epsilonsw1(1)= cos(Spar(Ntr)%Lamdasp(n))*cos(Spar(Ntr)%Lamdasp(n))*epsilonsw(1)+sin(Spar(Ntr)%Lamdasp(n))*sin(Spar(Ntr)%Lamdasp(n))*epsilonsw(2)-sin(Spar(Ntr)%Lamdasp(n))*cos(Spar(Ntr)%Lamdasp(n))*epsilonsw(3)	!epsilon xx local
                epsilonsw1(2)= sin(Spar(Ntr)%Lamdasp(n))*sin(Spar(Ntr)%Lamdasp(n))*epsilonsw(1)+cos(Spar(Ntr)%Lamdasp(n))*cos(Spar(Ntr)%Lamdasp(n))*epsilonsw(2)+sin(Spar(Ntr)%Lamdasp(n))*cos(Spar(Ntr)%Lamdasp(n))*epsilonsw(3)	!epsilon yy local
                epsilonsw1(3)= 2.0D0*cos(Spar(Ntr)%Lamdasp(n))*sin(Spar(Ntr)%Lamdasp(n))*epsilonsw(1)-2.0D0*cos(Spar(Ntr)%Lamdasp(n))*sin(Spar(Ntr)%Lamdasp(n))*epsilonsw(2)+cos(2.0D0*Spar(Ntr)%Lamdasp(n))*epsilonsw(3)			!epsilon xy local
                epsilonsw1(4)= cos(Spar(Ntr)%Lamdasp(n))*epsilonsw(4)+sin(Spar(Ntr)%Lamdasp(n))*epsilonsw(5)																						!epsilon yz local
                epsilonsw1(5)= -sin(Spar(Ntr)%Lamdasp(n))*epsilonsw(4)+cos(Spar(Ntr)%Lamdasp(n))*epsilonsw(5)																						!epsilon zx local

                E22	= Spar(Ntr)%Eesw(n)/(1.0D0-Spar(Ntr)%miuezsw(n)*Spar(Ntr)%miuzesw(n))

                sigmasw(1)=(epsilonsw1(1))/1000000.0D0					!sigma xx local
                sigmasw(2)=(E22*epsilonsw1(2))/1000000.0D0				!sigma yy local
                sigmasw(3)=(epsilonsw1(3))/1000000.0D0					!sigma xy local
                sigmasw(4)=(Spar(Ntr)%Gezsw(n)*epsilonsw1(4))/1000000.0D0		!sigma yz local
                sigmasw(5)=(epsilonsw1(5))/1000000.0D0					!sigma zx local	

                sigmavonMisessw=sqrt(((sigmasw(1)-sigmasw(2))*(sigmasw(1)-sigmasw(2))+sigmasw(2)*sigmasw(2)+(-sigmasw(1))*(-sigmasw(1)) &
                +6.0D0*(sigmasw(3)*sigmasw(3)+sigmasw(4)*sigmasw(4)+sigmasw(5)*sigmasw(5)))/2.0D0)

                if (sigmavonMisessw.GT.sigmavonMisesswmax) sigmavonMisesswmax=sigmavonMisessw
                
                write(sparweb_stress_unit,'(F10.5,6F14.6,7F20.10)',IOSTAT=ierror,IOMSG=msg) x,y,z,sigmavonMisessw,sigmasw(1),sigmasw(2),sigmasw(3),sigmasw(4),sigmasw(5),epsilonsw1(1),epsilonsw1(2),epsilonsw1(3),epsilonsw1(4),epsilonsw1(5)
            end do
        end do
    end do
    close(sparweb_stress_unit,IOSTAT=ierror,IOMSG=msg)
    IF(ierror.NE.0) Call abort_program()
   
end subroutine SparWebStresses

!--------------------------------------------------------------------------------------------
! Routine to calculate Spar Cap Stresses
subroutine SparCapStresses(k)

    use PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: SkinEtaTecplot
    use SparVariablesRoutines, only: Spar,SparGeometry,SparCapGeometry  
    use DeformationVariables, only: Deformation
    use StressStrainVariables
    USE FilePathModule, only: output_dir,sparcap_stress_unit,ierror,msg,abort_program
    USE maths, only: AxesTranformation, TrialPolynomials,AddTwoVectors

    implicit none

    integer, intent(in)   ::    k
    real(8) ::				    E22																	! stiffness element
    real(8) ::				    BdB(k*k)															! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)															! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k),TC45(k*k)
    real(8) ::				    x,y,z,z1scU,z2scU,z1scL,z2scL,zeta,eta                      		! point coordinates	
    integer ::				    l,n,yi,zi,ytotal,ztotal         								    ! iteration variables	  
    real(8) ::				    du0dx,du0dy,dv0dx,dv0dy,dw0dx,dw0dy,dphixdx,dphixdy,dphiydx,dphiydy,phixx,phiyy
    real(8) ::                  DJ,J1(2,2)
    !1 ponto em x
    ytotal=45				!numero de pontos em y
    ztotal=3				!numero de pontos em z
    
    if(Ntr == 1) then !first access
        open(unit=sparcap_stress_unit,file=trim(adjustl(output_dir))//"sparcapstresses.dat",STATUS='UNKNOWN',ACTION='Write',IOSTAT=ierror,IOMSG=msg)
        write(sparcap_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "SPARCAP STRESSES"'
        write(sparcap_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES =x,y,z,vonMises_Stress,stressxx,stressyy,stressxy,stressyz,stresszx,strainxx,strainyy,strainxy,strainyz,strainzx'
    else
        open(unit=sparcap_stress_unit,file=trim(adjustl(output_dir))//"sparcapstresses.dat",ACTION='Write',ACCESS = 'APPEND',IOSTAT=ierror,IOMSG=msg)
    end if

    if(Spar(Ntr)%Nsp.EQ.0) then
        close(sparcap_stress_unit,IOSTAT=ierror,IOMSG=msg)
        return
    end if

    sigmavonMisesscmax = 0.0D0
    sigmavonMisessc = 0.0D0

    do n=1,Spar(Ntr)%Nsp
        write(sparcap_stress_unit,'(A,I2.2,A,I2.2,A,I4.4,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="Planform ',Ntr,' SparCap ',n,'" I=',1,' J=',ytotal,' K=',ztotal
        do zi=0,ztotal-1
            do yi=0,ytotal-1
                !eta = -1.0D0+2.0D0*float(yi)/float(ytotal-1)
                eta = SkinEtaTecplot(ytotal-1,yi)
                zeta = -1.0D0+2.0D0*(Spar(Ntr)%ksp1(n)-Spar(Ntr)%ksp2(n)*eta)		  
                call AxesTranformation(1,zeta,eta,x,y)
                call AxesTranformation(3,zeta,eta,x,y,DJ,J1)
                call SparCapGeometry(2,n,zeta,eta,zeta,x,y,z1scU,z2scU,z1scL,z2scL)
                call TrialPolynomials(2,k,zeta,eta,BdB)
                call TrialPolynomials(3,k,zeta,eta,BBd)
                call AddTwoVectors(J1(1,1),BdB,J1(2,1),BBd,TC11)
                call AddTwoVectors(J1(1,2),BdB,J1(2,2),BBd,TC22)
                call TrialPolynomials(1,k,zeta,eta,TC45)

                z = z1scU+zi*(z2scU-z1scU)/dfloat(ztotal-1)
                epsilonsc=0.0D0
                sigmasc = 0.0D0
                du0dx	= 0.000D0
                du0dy	= 0.000D0
                dv0dx	= 0.000D0
                dv0dy	= 0.000D0
                dw0dx	= 0.000D0
                dw0dy	= 0.000D0
                dphixdx	= 0.000D0
                dphixdy	= 0.000D0
                dphiydx	= 0.000D0
                dphiydy	= 0.000D0
                phixx	= 0.000D0
                phiyy	= 0.000D0

                do l=1,k*k
                    du0dx	= du0dx+Deformation(Ntr)%q(l)*TC11(l)
                    du0dy	= du0dy+Deformation(Ntr)%q(l)*TC22(l)
                    dv0dx	= dv0dx+Deformation(Ntr)%q(l+k*k)*TC11(l)
                    dv0dy	= dv0dy+Deformation(Ntr)%q(l+k*k)*TC22(l)
                    dw0dx	= dw0dx+Deformation(Ntr)%q(l+2*k*k)*TC11(l)
                    dw0dy	= dw0dy+Deformation(Ntr)%q(l+2*k*k)*TC22(l)
                    dphixdx	= dphixdx+Deformation(Ntr)%q(l+3*k*k)*TC11(l)
                    dphixdy	= dphixdy+Deformation(Ntr)%q(l+3*k*k)*TC22(l)
                    dphiydx	= dphiydx+Deformation(Ntr)%q(l+4*k*k)*TC11(l)
                    dphiydy	= dphiydy+Deformation(Ntr)%q(l+4*k*k)*TC22(l)
                    phixx	= phixx+Deformation(Ntr)%q(l+3*k*k)*TC45(l)
                    phiyy	= phiyy+Deformation(Ntr)%q(l+4*k*k)*TC45(l)
                end do

                epsilonsc(1)= du0dx+z*dphixdx														!epsilon xx global
                epsilonsc(2)= dv0dy+z*dphiydy														!epsilon yy global
                epsilonsc(3)= 1.0D0/2.0D0*(du0dy+dv0dx+z*(dphixdy+dphiydx))								!epsilon xy global
                epsilonsc(4)= 1.0D0/2.0D0*(phiyy+dw0dy)													!epsilon yz global
                epsilonsc(5)= 1.0D0/2.0D0*(phixx+dw0dx)													!epsilon zx global


                epsilonsc1(1)= cos(Spar(Ntr)%Lamdasp(n))*cos(Spar(Ntr)%Lamdasp(n))*epsilonsc(1)+sin(Spar(Ntr)%Lamdasp(n))*sin(Spar(Ntr)%Lamdasp(n))*epsilonsc(2)-sin(Spar(Ntr)%Lamdasp(n))*cos(Spar(Ntr)%Lamdasp(n))*epsilonsc(3)	!epsilon xx local
                epsilonsc1(2)= sin(Spar(Ntr)%Lamdasp(n))*sin(Spar(Ntr)%Lamdasp(n))*epsilonsc(1)+cos(Spar(Ntr)%Lamdasp(n))*cos(Spar(Ntr)%Lamdasp(n))*epsilonsc(2)+sin(Spar(Ntr)%Lamdasp(n))*cos(Spar(Ntr)%Lamdasp(n))*epsilonsc(3)	!epsilon yy local
                epsilonsc1(3)= 2.0D0*cos(Spar(Ntr)%Lamdasp(n))*sin(Spar(Ntr)%Lamdasp(n))*epsilonsc(1)-2.0D0*cos(Spar(Ntr)%Lamdasp(n))*sin(Spar(Ntr)%Lamdasp(n))*epsilonsc(2)+cos(2.0D0*Spar(Ntr)%Lamdasp(n))*epsilonsc(3)			!epsilon xy local
                epsilonsc1(4)= cos(Spar(Ntr)%Lamdasp(n))*epsilonsc(4)+sin(Spar(Ntr)%Lamdasp(n))*epsilonsc(5)																						!epsilon yz local
                epsilonsc1(5)= -sin(Spar(Ntr)%Lamdasp(n))*epsilonsc(4)+cos(Spar(Ntr)%Lamdasp(n))*epsilonsc(5)																						!epsilon zx local

                E22	= Spar(Ntr)%Eesc(n)/(1.0D0-Spar(Ntr)%miuezsc(n)*Spar(Ntr)%miuzesc(n))

                sigmasc(1)=0.0D0*epsilonsc1(1)/1000000.0D0	    					!sigma xx local
                sigmasc(2)=E22*epsilonsc1(2)/1000000.0D0							!sigma yy local
                sigmasc(3)=0.0D0*epsilonsc1(3)/1000000.0D0							!sigma xy local
                sigmasc(4)=0.0D0*epsilonsc1(4)/1000000.0D0							!sigma yz local
                sigmasc(5)=0.0D0*epsilonsc1(5)/1000000.0D0							!sigma zx local																				!sigma yz

                sigmavonMisessc=sqrt(((sigmasc(1)-sigmasc(2))*(sigmasc(1)-sigmasc(2))+sigmasc(2)*sigmasc(2)+(-sigmasc(1))*(-sigmasc(1)) &
                +6.0D0*(sigmasc(3)*sigmasc(3)+sigmasc(4)*sigmasc(4)+sigmasc(5)*sigmasc(5)))/2.0D0)

                if (sigmavonMisessc.GT.sigmavonMisesscmax) sigmavonMisesscmax=sigmavonMisessc

                write(sparcap_stress_unit,'(F10.5,6F14.6,7F20.10)',IOSTAT=ierror,IOMSG=msg) x,y,z,sigmavonMisessc,sigmasc(1),sigmasc(2),sigmasc(3),sigmasc(4),sigmasc(5),epsilonsc(1),epsilonsc(2),epsilonsc(3),epsilonsc(4),epsilonsc(5)
            end do
        end do
    end do	
    close(sparcap_stress_unit,IOSTAT=ierror,IOMSG=msg)
    IF(ierror.NE.0) Call abort_program()

end subroutine SparCapStresses

!--------------------------------------------------------------------------------------------
! Routine to calculate Rib Web Stresses
subroutine RibWebStresses(k)
    use PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: SkinZetaTecplot
    use RibVariablesRoutines, only: Rib,RibGeometry,RibWebGeometry 
    use DeformationVariables, only: Deformation
    use StressStrainVariables
    USE FilePathModule, only: output_dir,ribweb_stress_unit,ierror,msg,abort_program
    USE maths, only: AxesTranformation,TrialPolynomials,AddTwoVectors
    
    implicit none

    integer, intent(in)   ::    k
    real(8) ::				    E11																	! stiffness element
    real(8) ::				    BdB(k*k)															! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)															! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k),TC45(k*k)			
    real(8) ::				    x,y,z,z1rw,z2rw,zeta,eta				! point coordinates
    integer ::				    l,n,xi,zi,xtotal,ztotal											! iteration variables	  
    real(8) ::				    du0dx,du0dy,dv0dx,dv0dy,dw0dx,dw0dy,dphixdx,dphixdy,dphiydx,dphiydy,phixx,phiyy
    real(8) ::                  DJ,J1(2,2)
    !Geometry errors search
    integer ::                  xtotal_new,count
    real(8) ::                  zeta_next,eta_next
    logical ::                  error,error_next,TE
    logical, allocatable ::     accept_array(:)


    !xtotal=int(croot(Npl)*40.0D0) !25 !numero de pontos em x
    xtotal=int(1.0D0/0.01D0)
    !1 ponto em y
    ztotal=9				!numero de pontos em z
    allocate(accept_array(0:xtotal-1))

    if(Ntr == 1) then !first access
        open(unit=ribweb_stress_unit,file=trim(adjustl(output_dir))//"ribwebstresses.dat",STATUS='UNKNOWN',ACTION='Write',IOSTAT=ierror,IOMSG=msg)
        write(ribweb_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "RIBWEB STRESSES"'
        write(ribweb_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES =x,y,z,vonMises_Stress,stressxx,stressyy,stressxy,stressyz,stresszx,strainxx,strainyy,strainxy,strainyz,strainzx'
    else
        open(unit=ribweb_stress_unit,file=trim(adjustl(output_dir))//"ribwebstresses.dat",ACTION='Write',ACCESS = 'APPEND',IOSTAT=ierror,IOMSG=msg)
    end if

    if(Rib(Ntr)%Nrb.EQ.0) then
        close(ribweb_stress_unit,IOSTAT=ierror,IOMSG=msg)
        return
    end if

    sigmavonMisesrwmax = 0.0D0
    sigmavonMisesrw = 0.0D0	

    do n=1,Rib(Ntr)%Nrb
        !Search for invalid rib web geometry and reduce total number of points
        xtotal_new=xtotal
        TE=.False.
        count=0
        do xi=0,xtotal-1
            if(TE) then
                if(count >= 4) then 
                    TE=.False.
                else
                    count=count+1
                    accept_array(xi)=.True.
                    cycle
                end if
            end if   
            zeta = SkinZetaTecplot(1,xtotal,xi)
            eta  = -1.0D0+2.0D0*(Rib(Ntr)%krb1(n)-Rib(Ntr)%krb2(n)*zeta)
            call AxesTranformation(1,zeta,eta,x,y)
            call RibWebGeometry(2,n,zeta,eta,eta,x,y,z1rw,z2rw,error)
            zeta_next = SkinZetaTecplot(1,xtotal,xi+1)
            eta_next	= -1.0D0+2.0D0*(Rib(Ntr)%krb1(n)-Rib(Ntr)%krb2(n)*zeta_next)
            call AxesTranformation(1,zeta_next,eta_next,x,y)
            call RibWebGeometry(2,n,zeta_next,eta_next,eta_next,x,y,z1rw,z2rw,error_next)
            !Select point 
            if(error.AND.error_next) then !problem: remove point
                xtotal_new=xtotal_new-1
                accept_array(xi)=.False.
            elseif(error.AND..not.error_next) then !problem near leading edge: allow point
                !xtotal_new=xtotal_new
                accept_array(xi)=.True.
                !LE=.True.
            elseif(.not.error.AND.error_next) then !problem near trailing edge: allow sebsequent points depending on count! (on top)
                !xtotal_new=xtotal_new
                accept_array(xi)=.True.
                !if(xi+1 <= xtotal-1) accept_array(xi+1)=.True.
                TE=.True.
                cycle
            elseif(.not.error.AND..not.error_next) then !no problem
                !xtotal_new=xtotal_new
                accept_array(xi)=.True.
            end if
            !if(error) xtotal_new=xtotal_new-1
        end do
        write(ribweb_stress_unit,'(A,I2.2,A,I2.2,A,I4.4,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="Planform ',Ntr,' RibWeb ',n,'" I=',xtotal_new,' J=',1,' K=',ztotal
        do zi=1,ztotal
            do xi=0,xtotal-1
                IF(.not.accept_array(xi)) then
                    cycle
                end if
                zeta = SkinZetaTecplot(1,xtotal,xi)
                eta  = -1.0D0+2.0D0*(Rib(Ntr)%krb1(n)-Rib(Ntr)%krb2(n)*zeta)
                call AxesTranformation(1,zeta,eta,x,y) !Compute x,y
                call RibWebGeometry(2,n,zeta,eta,eta,x,y,z1rw,z2rw,error)
                call AxesTranformation(3,zeta,eta,x,y,DJ,J1) !Compute DJ,J1
                call TrialPolynomials(2,k,zeta,eta,BdB)
                call TrialPolynomials(3,k,zeta,eta,BBd)
                call AddTwoVectors(J1(1,1),BdB,J1(2,1),BBd,TC11)
                call AddTwoVectors(J1(1,2),BdB,J1(2,2),BBd,TC22)
                call TrialPolynomials(1,k,zeta,eta,TC45)

                z= z1rw+dfloat(zi-1)*(z2rw-z1rw)/dfloat(ztotal-1)

                epsilonrw=0.0D0
                sigmarw=0.0D0

                du0dx	= 0.000D0
                du0dy	= 0.000D0
                dv0dx	= 0.000D0
                dv0dy	= 0.000D0
                dw0dx	= 0.000D0
                dw0dy	= 0.000D0
                dphixdx = 0.000D0
                dphixdy = 0.000D0
                dphiydx = 0.000D0
                dphiydy = 0.000D0
                phixx	= 0.000D0
                phiyy	= 0.000D0

                do l=1,k*k
                    du0dx	= du0dx+Deformation(Ntr)%q(l)*TC11(l)
                    du0dy	= du0dy+Deformation(Ntr)%q(l)*TC22(l)
                    dv0dx	= dv0dx+Deformation(Ntr)%q(l+k*k)*TC11(l)
                    dv0dy	= dv0dy+Deformation(Ntr)%q(l+k*k)*TC22(l)
                    dw0dx	= dw0dx+Deformation(Ntr)%q(l+2*k*k)*TC11(l)
                    dw0dy	= dw0dy+Deformation(Ntr)%q(l+2*k*k)*TC22(l)
                    dphixdx	= dphixdx+Deformation(Ntr)%q(l+3*k*k)*TC11(l)
                    dphixdy	= dphixdy+Deformation(Ntr)%q(l+3*k*k)*TC22(l)
                    dphiydx	= dphiydx+Deformation(Ntr)%q(l+4*k*k)*TC11(l)
                    dphiydy	= dphiydy+Deformation(Ntr)%q(l+4*k*k)*TC22(l)
                    phixx	= phixx+Deformation(Ntr)%q(l+3*k*k)*TC45(l)
                    phiyy	= phiyy+Deformation(Ntr)%q(l+4*k*k)*TC45(l)
                end do

                epsilonrw(1)= du0dx+z*dphixdx												!epsilon xx global
                epsilonrw(2)= dv0dy+z*dphiydy												!epsilon yy global
                epsilonrw(3)= 1.0D0/2.0D0*(du0dy+dv0dx+z*(dphixdy+dphiydx))						!epsilon xy global
                epsilonrw(4)= 1.0D0/2.0D0*(phiyy+dw0dy)											!epsilon yz global
                epsilonrw(5)= 1.0D0/2.0D0*(phixx+dw0dx)											!epsilon zx global

                epsilonrw1(1)= cos(Rib(Ntr)%Lamdarb(n))*cos(Rib(Ntr)%Lamdarb(n))*epsilonrw(1)+sin(Rib(Ntr)%Lamdarb(n))*sin(Rib(Ntr)%Lamdarb(n))*epsilonrw(2)-sin(Rib(Ntr)%Lamdarb(n))*cos(Rib(Ntr)%Lamdarb(n))*epsilonrw(3)	!epsilon xx local
                epsilonrw1(2)= sin(Rib(Ntr)%Lamdarb(n))*sin(Rib(Ntr)%Lamdarb(n))*epsilonrw(1)+cos(Rib(Ntr)%Lamdarb(n))*cos(Rib(Ntr)%Lamdarb(n))*epsilonrw(2)+sin(Rib(Ntr)%Lamdarb(n))*cos(Rib(Ntr)%Lamdarb(n))*epsilonrw(3)	!epsilon yy local
                epsilonrw1(3)= 2.0D0*cos(Rib(Ntr)%Lamdarb(n))*sin(Rib(Ntr)%Lamdarb(n))*epsilonrw(1)-2.0D0*cos(Rib(Ntr)%Lamdarb(n))*sin(Rib(Ntr)%Lamdarb(n))*epsilonrw(2)+cos(2.0D0*Rib(Ntr)%Lamdarb(n))*epsilonrw(3)			!epsilon xy local
                epsilonrw1(4)= cos(Rib(Ntr)%Lamdarb(n))*epsilonrw(4)+sin(Rib(Ntr)%Lamdarb(n))*epsilonrw(5)																						!epsilon yz local
                epsilonrw1(5)= -sin(Rib(Ntr)%Lamdarb(n))*epsilonrw(4)+cos(Rib(Ntr)%Lamdarb(n))*epsilonrw(5)																						!epsilon zx local

                E11	 = Rib(Ntr)%Eerw(n)/(1.0-Rib(Ntr)%miuezrw(n)*Rib(Ntr)%miuzerw(n))

                sigmarw(1)=E11*epsilonrw1(1)/1000000.0D0					!sigma xx local
                sigmarw(2)=0.0D0*epsilonrw1(2)/1000000.0D0					!sigma yy local
                sigmarw(3)=0.0D0*epsilonrw1(3)/1000000.0D0					!sigma xy local
                sigmarw(4)=0.0D0*epsilonrw1(4)/1000000.0D0					!sigma yz local
                sigmarw(5)=Rib(Ntr)%Gezrw(n)*(epsilonrw1(5))/1000000.0D0    !sigma zx local	

                sigmavonMisesrw=sqrt(((sigmarw(1)-sigmarw(2))*(sigmarw(1)-sigmarw(2))+sigmarw(2)*sigmarw(2)+(-sigmarw(1))*(-sigmarw(1)) &
                +6.0D0*(sigmarw(3)*sigmarw(3)+sigmarw(4)*sigmarw(4)+sigmarw(5)*sigmarw(5)))/2.0D0)

                if (sigmavonMisesrw.GT.sigmavonMisesrwmax) sigmavonMisesrwmax=sigmavonMisesrw
                
                write(ribweb_stress_unit,'(F10.5,6F14.6,7F20.10)',IOSTAT=ierror,IOMSG=msg) x,y,z,sigmavonMisesrw,sigmarw(1),sigmarw(2),sigmarw(3),sigmarw(4),sigmarw(5),epsilonrw1(1),epsilonrw1(2),epsilonrw1(3),epsilonrw1(4),epsilonrw1(5)
            end do
        end do
    end do
    close(ribweb_stress_unit,IOSTAT=ierror,IOMSG=msg)
    IF(ierror.NE.0) Then
        Call abort_program
    End IF  
    deallocate(accept_array)
end subroutine RibWebStresses

!--------------------------------------------------------------------------------------------
! Routine to calculate Rib Cap Stresses
subroutine RibCapStresses(k)
    use PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: SkinZetaTecplot
    use RibVariablesRoutines, only: Rib,RibGeometry,RibCapGeometry 
    use DeformationVariables, only: Deformation
    use StressStrainVariables
    USE FilePathModule, only: output_dir,ribcap_stress_unit,ierror,msg,abort_program
    USE maths, only: AxesTranformation,TrialPolynomials,AddTwoVectors
    
    implicit none

    integer, intent(in)   ::    k
    real(8) ::				    E11																	! stiffness element
    real(8) ::				    BdB(k*k)															! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)															! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k)
    real(8) ::				    TC45(k*k)
    real(8) ::				    x,y,z,z1rcU,z2rcU,z1rcL,z2rcL,zeta,eta                          	! point coordinates
    integer ::				    l,n,xi,zi,xtotal,ztotal											! iteration variables	  
    real(8) ::				    du0dx,du0dy,dv0dx,dv0dy,dw0dx,dw0dy,dphixdx,dphixdy,dphiydx,dphiydy,phixx,phiyy
    real(8) ::                  DJ,J1(2,2)
    logical ::                  error !,error_next,TE


    !xtotal=int(croot(Npl)*40.0D0) !25 !numero de pontos em x
    xtotal=int(1.0D0/0.01D0)
    !1 ponto em y
    ztotal=7				!numero de pontos em z
    
    if(Ntr == 1) then !first access
        open(unit=ribcap_stress_unit,file=trim(adjustl(output_dir))//"ribcapstresses.dat",STATUS='UNKNOWN',ACTION='Write',IOSTAT=ierror,IOMSG=msg)
        write(ribcap_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "RIBCAP STRESSES"'
        write(ribcap_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES =x,y,z,vonMises_Stress,stressxx,stressyy,stressxy,stressyz,stresszx,strainxx,strainyy,strainxy,strainyz,strainzx'
    else
        open(unit=ribcap_stress_unit,file=trim(adjustl(output_dir))//"ribcapstresses.dat",ACTION='Write',ACCESS='APPEND',IOSTAT=ierror,IOMSG=msg)
    end if

    if(Rib(Ntr)%Nrb.EQ.0) then
        close(ribcap_stress_unit,IOSTAT=ierror,IOMSG=msg)
        return
    end if

    sigmavonMisesrcmax = 0.0D0
    sigmavonMisesrc = 0.0D0	
    do n=1,Rib(Ntr)%Nrb
        write(ribcap_stress_unit,'(A,I2.2,A,I2.2,A,I4.4,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="Planform ',Ntr,' RibCap ',n,'" I=',xtotal-6,' J=',1,' K=',ztotal
        do zi=0,ztotal-1
            do xi=2,xtotal-5 
                zeta = SkinZetaTecplot(1,xtotal,xi)
                eta = -1.0D0+2.0D0*(Rib(Ntr)%krb1(n)-Rib(Ntr)%krb2(n)*zeta)

                call AxesTranformation(1,zeta,eta,x,y)		!Compute x,y
                call AxesTranformation(3,zeta,eta,x,y,DJ,J1) !Compute DJ,J1	  
                call RibCapGeometry(2,n,zeta,eta,eta,x,y,z1rcU,z2rcU,z1rcL,z2rcL,error)
                !if(error) write(*,*) "error ribcap"
                call TrialPolynomials(2,k,zeta,eta,BdB)
                call TrialPolynomials(3,k,zeta,eta,BBd)
                call AddTwoVectors(J1(1,1),BdB,J1(2,1),BBd,TC11)
                call AddTwoVectors(J1(1,2),BdB,J1(2,2),BBd,TC22)
                call TrialPolynomials(1,k,zeta,eta,TC45)

                z = z1rcU+zi*(z2rcU-z1rcU)/dfloat(ztotal-1)
                epsilonrc=0.0D0
                sigmarc=0.0D0
                du0dx	= 0.000D0
                du0dy	= 0.000D0
                dv0dx	= 0.000D0
                dv0dy	= 0.000D0
                dw0dx	= 0.000D0
                dw0dy	= 0.000D0
                dphixdx	= 0.000D0
                dphixdy	= 0.000D0
                dphiydx	= 0.000D0
                dphiydy	= 0.000D0
                phixx	= 0.000D0
                phiyy	= 0.000D0

                do l=1,k*k
                    du0dx	= du0dx+Deformation(Ntr)%q(l)*TC11(l)
                    du0dy	= du0dy+Deformation(Ntr)%q(l)*TC22(l)
                    dv0dx	= dv0dx+Deformation(Ntr)%q(l+k*k)*TC11(l)
                    dv0dy	= dv0dy+Deformation(Ntr)%q(l+k*k)*TC22(l)
                    dw0dx	= dw0dx+Deformation(Ntr)%q(l+2*k*k)*TC11(l)
                    dw0dy	= dw0dy+Deformation(Ntr)%q(l+2*k*k)*TC22(l)
                    dphixdx	= dphixdx+Deformation(Ntr)%q(l+3*k*k)*TC11(l)
                    dphixdy	= dphixdy+Deformation(Ntr)%q(l+3*k*k)*TC22(l)
                    dphiydx	= dphiydx+Deformation(Ntr)%q(l+4*k*k)*TC11(l)
                    dphiydy	= dphiydy+Deformation(Ntr)%q(l+4*k*k)*TC22(l)
                    phixx	= phixx+Deformation(Ntr)%q(l+3*k*k)*TC45(l)
                    phiyy	= phiyy+Deformation(Ntr)%q(l+4*k*k)*TC45(l)
                end do
                epsilonrc(1)= du0dx+z*dphixdx												!epsilon xx global
                epsilonrc(2)= dv0dy+z*dphiydy												!epsilon yy global
                epsilonrc(3)= 1.0D0/2.0D0*(du0dy+dv0dx+z*(dphixdy+dphiydx))						!epsilon xy global
                epsilonrc(4)= 1.0D0/2.0D0*(phiyy+dw0dy)											!epsilon yz global
                epsilonrc(5)= 1.0D0/2.0D0*(phixx+dw0dx)											!epsilon zx global
                epsilonrc1(1)= cos(Rib(Ntr)%Lamdarb(n))*cos(Rib(Ntr)%Lamdarb(n))*epsilonrc(1)+sin(Rib(Ntr)%Lamdarb(n))*sin(Rib(Ntr)%Lamdarb(n))*epsilonrc(2)-sin(Rib(Ntr)%Lamdarb(n))*cos(Rib(Ntr)%Lamdarb(n))*epsilonrc(3)	!epsilon xx local
                epsilonrc1(2)= sin(Rib(Ntr)%Lamdarb(n))*sin(Rib(Ntr)%Lamdarb(n))*epsilonrc(1)+cos(Rib(Ntr)%Lamdarb(n))*cos(Rib(Ntr)%Lamdarb(n))*epsilonrc(2)+sin(Rib(Ntr)%Lamdarb(n))*cos(Rib(Ntr)%Lamdarb(n))*epsilonrc(3)	!epsilon yy local
                epsilonrc1(3)= 2.0D0*cos(Rib(Ntr)%Lamdarb(n))*sin(Rib(Ntr)%Lamdarb(n))*epsilonrc(1)-2.0D0*cos(Rib(Ntr)%Lamdarb(n))*sin(Rib(Ntr)%Lamdarb(n))*epsilonrc(2)+cos(2.0D0*Rib(Ntr)%Lamdarb(n))*epsilonrc(3)			!epsilon xy local
                epsilonrc1(4)= cos(Rib(Ntr)%Lamdarb(n))*epsilonrc(4)+sin(Rib(Ntr)%Lamdarb(n))*epsilonrc(5)																						!epsilon yz local
                epsilonrc1(5)= -sin(Rib(Ntr)%Lamdarb(n))*epsilonrc(4)+cos(Rib(Ntr)%Lamdarb(n))*epsilonrc(5)																						!epsilon zx local

                E11	= Rib(Ntr)%Eerw(n)/(1.0D0-Rib(Ntr)%miuezrw(n)*Rib(Ntr)%miuzerw(n))

                sigmarc(1)=E11*epsilonrc1(1)/1000000.0D0					!sigma xx local
                sigmarc(2)=0.0D0*epsilonrc1(2)/1000000.0D0					!sigma yy local
                sigmarc(3)=0.0D0*epsilonrc1(3)/1000000.0D0					!sigma xy local
                sigmarc(4)=0.0D0*epsilonrc1(4)/1000000.0D0					!sigma yz local
                sigmarc(5)=0.0D0*(epsilonrc1(5))/1000000.0D0				!sigma zx local	

                sigmavonMisesrc=sqrt(((sigmarc(1)-sigmarc(2))*(sigmarc(1)-sigmarc(2))+sigmarc(2)*sigmarc(2)+(-sigmarc(1))*(-sigmarc(1)) &
                +6.0D0*(sigmarc(3)*sigmarc(3)+sigmarc(4)*sigmarc(4)+sigmarc(5)*sigmarc(5)))/2.0D0)

                if (sigmavonMisesrc.GT.sigmavonMisesrcmax) sigmavonMisesrcmax=sigmavonMisesrc

                write(ribcap_stress_unit,'(F10.5,6F14.6,7F20.10)',IOSTAT=ierror,IOMSG=msg) x,y,z,sigmavonMisesrc,sigmarc(1),sigmarc(2),sigmarc(3),sigmarc(4),sigmarc(5),epsilonrc(1),epsilonrc(2),epsilonrc(3),epsilonrc(4),epsilonrc(5)
            end do
        end do
    end do
    close(ribcap_stress_unit,IOSTAT=ierror,IOMSG=msg)
    IF(ierror.NE.0) Call abort_program()

end subroutine RibCapStresses

!--------------------------------------------------------------------------------------------
! Routine to calculate Stringer Stresses
subroutine StringerStresses(k)

    use PlanformVariables, only: Ntr
    use SkinVariablesRoutines, only: SkinEtaTecplot
    use StringerVariablesRoutines, only: Stringer,StringerGeometry1,StringerGeometry2    
    use DeformationVariables, only: Deformation
    use StressStrainVariables
    USE FilePathModule, only: output_dir,stringer_stress_unit,ierror,msg,abort_program
    USE maths, only: AxesTranformation,TrialPolynomials,AddTwoVectors

    implicit none

    integer, intent(in) ::      k
    real(8) ::                  E22																	! stiffness element
    real(8) ::                  BdB(k*k)															! Legendre polynimial multiplication vector
    real(8) ::				    BBd(k*k)															! Legendre polynimial multiplication vector
    real(8) ::				    TC11(k*k),TC22(k*k),TC45(k*k)
    real(8) ::				    x,y,z,z1stU,z2stU,z1stL,z2stL,zeta,eta    		! point coordinates	
    integer ::				    l,n,xi,yi,xtotal,ytotal,ztotal								! iteration variables	  
    real(8) ::				    du0dx,du0dy,dv0dx,dv0dy,dw0dx,dw0dy,dphixdx,dphixdy,dphiydx,dphiydy,phixx,phiyy
    real(8) ::                  DJ,J1(2,2)

    xtotal=3				!numero de pontos em x
    ytotal=25				!numero de pontos em y
    ztotal=1				!1 ponto em z

    if(Ntr == 1) then !first access
        open(unit=stringer_stress_unit,file=trim(adjustl(output_dir))//"stringerstresses.dat",STATUS='UNKNOWN',ACTION='Write',IOSTAT=ierror,IOMSG=msg)
        write(stringer_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'TITLE = "STRINGER STRESSES"'
        write(stringer_stress_unit,'(A)',IOSTAT=ierror,IOMSG=msg) 'VARIABLES =x,y,z,vonMises_Stress,stressxx,stressyy,stressxy,stressyz,stresszx,strainxx,strainyy,strainxy,strainyz,strainzx'
    else
        open(unit=stringer_stress_unit,file=trim(adjustl(output_dir))//"stringerstresses.dat",ACTION='Write',ACCESS='APPEND',IOSTAT=ierror,IOMSG=msg)
    end if
    
    if(Stringer(Ntr)%Nst.EQ.0) then
        close(stringer_stress_unit,IOSTAT=ierror,IOMSG=msg)
        return
    end if

    sigmavonMisesstmax = 0.0D0
    sigmavonMisesst = 0.0D0

    do n=1,Stringer(Ntr)%Nst	
        write(stringer_stress_unit,'(A,I2.2,A,I2.2,A,I4.4,A,I4.4,A,I4.4)',IOSTAT=ierror,IOMSG=msg) 'ZONE T="Planform ',Ntr,' Stringer ',n,'" I=',xtotal,' J=',ytotal,' K=',ztotal
        do yi=0,ytotal-1	  		
            eta  = SkinEtaTecplot(ytotal-1,yi)
            zeta = -1.0D0+2.0D0*(Stringer(Ntr)%kst1(n)-Stringer(Ntr)%kst2(n)*eta)		  
            call AxesTranformation(1,zeta,eta,x,y)		!Compute x,y			!****************************************
            call AxesTranformation(3,zeta,eta,x,y,DJ,J1) !Compute DJ,J1	
            do xi=0,xtotal-1
                call StringerGeometry2(2,n,zeta,eta,zeta,x,y,z1stU,z2stU,z1stL,z2stL)
                call TrialPolynomials(2,k,zeta,eta,BdB)
                call TrialPolynomials(3,k,zeta,eta,BBd)
                call AddTwoVectors(J1(1,1),BdB,J1(2,1),BBd,TC11)
                call AddTwoVectors(J1(1,2),BdB,J1(2,2),BBd,TC22)
                call TrialPolynomials(1,k,zeta,eta,TC45)

                z = z1stL+(z2stL-z1stL)/2.0D0

                epsilonst=0.0D0
                sigmast=0.0D0

                du0dx	= 0.000D0
                du0dy	= 0.000D0
                dv0dx	= 0.000D0
                dv0dy	= 0.000D0
                dw0dx	= 0.000D0
                dw0dy	= 0.000D0
                dphixdx = 0.000D0
                dphixdy = 0.000D0
                dphiydx = 0.000D0
                dphiydy = 0.000D0
                phixx	= 0.000D0
                phiyy	= 0.000D0

                do l=1,k*k
                    du0dx	= du0dx+Deformation(Ntr)%q(l)*TC11(l)
                    du0dy	= du0dy+Deformation(Ntr)%q(l)*TC22(l)
                    dv0dx	= dv0dx+Deformation(Ntr)%q(l+k*k)*TC11(l)
                    dv0dy	= dv0dy+Deformation(Ntr)%q(l+k*k)*TC22(l)
                    dw0dx	= dw0dx+Deformation(Ntr)%q(l+2*k*k)*TC11(l)
                    dw0dy	= dw0dy+Deformation(Ntr)%q(l+2*k*k)*TC22(l)
                    dphixdx	= dphixdx+Deformation(Ntr)%q(l+3*k*k)*TC11(l)
                    dphixdy	= dphixdy+Deformation(Ntr)%q(l+3*k*k)*TC22(l)
                    dphiydx	= dphiydx+Deformation(Ntr)%q(l+4*k*k)*TC11(l)
                    dphiydy	= dphiydy+Deformation(Ntr)%q(l+4*k*k)*TC22(l)
                    phixx	= phixx+Deformation(Ntr)%q(l+3*k*k)*TC45(l)
                    phiyy	= phiyy+Deformation(Ntr)%q(l+4*k*k)*TC45(l)
                end do

                epsilonst(1)= du0dx+z*dphixdx												!epsilon xx global
                epsilonst(2)= dv0dy+z*dphiydy												!epsilon yy global
                epsilonst(3)= 1.0D0/2.0D0*(du0dy+dv0dx+z*(dphixdy+dphiydx))						!epsilon xy global
                epsilonst(4)= 1.0D0/2.0D0*(phiyy+dw0dy)											!epsilon yz global
                epsilonst(5)= 1.0D0/2.0D0*(phixx+dw0dx)											!epsilon zx global

                epsilonst1(1)= cos(Stringer(Ntr)%Lamdast(n))*cos(Stringer(Ntr)%Lamdast(n))*epsilonst(1)+sin(Stringer(Ntr)%Lamdast(n))*sin(Stringer(Ntr)%Lamdast(n))*epsilonst(2)-sin(Stringer(Ntr)%Lamdast(n))*cos(Stringer(Ntr)%Lamdast(n))*epsilonst(3)	!epsilon xx local
                epsilonst1(2)= sin(Stringer(Ntr)%Lamdast(n))*sin(Stringer(Ntr)%Lamdast(n))*epsilonst(1)+cos(Stringer(Ntr)%Lamdast(n))*cos(Stringer(Ntr)%Lamdast(n))*epsilonst(2)+sin(Stringer(Ntr)%Lamdast(n))*cos(Stringer(Ntr)%Lamdast(n))*epsilonst(3)	!epsilon yy local
                epsilonst1(3)= 2.0D0*cos(Stringer(Ntr)%Lamdast(n))*sin(Stringer(Ntr)%Lamdast(n))*epsilonst(1)-2.0D0*cos(Stringer(Ntr)%Lamdast(n))*sin(Stringer(Ntr)%Lamdast(n))*epsilonst(2)+cos(2.0D0*Stringer(Ntr)%Lamdast(n))*epsilonst(3)			!epsilon xy local
                epsilonst1(4)= cos(Stringer(Ntr)%Lamdast(n))*epsilonst(4)+sin(Stringer(Ntr)%Lamdast(n))*epsilonst(5)																						!epsilon yz local
                epsilonst1(5)= -sin(Stringer(Ntr)%Lamdast(n))*epsilonst(4)+cos(Stringer(Ntr)%Lamdast(n))*epsilonst(5)																						!epsilon zx local

                E22	= Stringer(Ntr)%Eest(n)/(1.0D0-Stringer(Ntr)%miuezst(n)*Stringer(Ntr)%miuzest(n))

                sigmast(1)=0.0D0*epsilonst1(1)/1000000.0D0							!sigma xx local
                sigmast(2)=E22*epsilonst1(2)/1000000.0D0							!sigma yy local
                sigmast(3)=0.0D0*epsilonst1(3)/1000000.0D0							!sigma xy local
                sigmast(4)=0.0D0*epsilonst1(4)/1000000.0D0							!sigma yz local
                sigmast(5)=0.0D0*epsilonst1(5)/1000000.0D0							!sigma zx local		

                sigmavonMisesst=sqrt(((sigmast(1)-sigmast(2))*(sigmast(1)-sigmast(2))+sigmast(2)*sigmast(2)+(-sigmast(1))*(-sigmast(1)) &
                +6.0D0*(sigmast(3)*sigmast(3)+sigmast(4)*sigmast(4)+sigmast(5)*sigmast(5)))/2.0D0)

                if (sigmavonMisesst.GT.sigmavonMisesstmax) sigmavonMisesstmax=sigmavonMisesst
                !Error detected and corrected - invalid data write 
                write(stringer_stress_unit,'(F10.5,6F14.6,7F20.10)',IOSTAT=ierror,IOMSG=msg) x,y,z,sigmavonMisesst,sigmast(1),sigmast(2),sigmast(3),sigmast(4),sigmast(5),epsilonst(1),epsilonst(2),epsilonst(3),epsilonst(4),epsilonst(5)
            end do
        end do
    end do
    close(stringer_stress_unit,IOSTAT=ierror,IOMSG=msg)
    IF(ierror.NE.0) Call abort_program()

end subroutine StringerStresses

!--------------------------------------------------------------------------------------------
! Routine to calculate Stresses
subroutine Stress()

    use PolynomialCoefficients, only: k
    !USE FilePathModule
    use StressStrainVariables, only: sigmavonMisesskmax,sigmavonMisesswmax,sigmavonMisesscmax,sigmavonMisesrwmax,sigmavonMisesrcmax,sigmavonMisesstmax,zero_sigmavonMises
    use SkinVariablesRoutines, only: Skin 
    use SparVariablesRoutines, only: Spar
    use RibVariablesRoutines, only: Rib
    use StringerVariablesRoutines, only: Stringer
    use GeneralVariables, only: OutputStress
    use PlanformVariables, only: Ntr

    implicit none

    if(.not.OutputStress) return
    write(*,'(1X,A,1X,I1,1X,A)') 'Computing trapezoid', Ntr, 'stresses...'
    Call zero_sigmavonMises()

    !$OMP PARALLEL SECTIONS
    !$OMP SECTION
    call SkinStresses(k)
    !$OMP SECTION
    call SparWebStresses(k)
    !$OMP SECTION
    call SparCapStresses(k)
    !$OMP SECTION
    call RibWebStresses(k)
    !$OMP SECTION
    call RibCapStresses(k)
    !$OMP SECTION
    call StringerStresses(k)
    !$OMP END PARALLEL SECTIONS
    
    !Print maximum Stress Values
    !Skin stress
    IF(Skin(Ntr)%Nsk > 0) write(*,*) 'Max Skin Von Mises Stress:    ',sigmavonMisesskmax,'[MPa]'
    !Sparweb stress
    IF(Spar(Ntr)%Nsp > 0) write(*,*) 'Max Spar Web Von Mises Stress:',sigmavonMisesswmax,'[MPa]'
    !Sparcap stress    
    IF(Spar(Ntr)%Nsp > 0) write(*,*) 'Max Spar Cap Von Mises Stress:',sigmavonMisesscmax,'[MPa]'
    !Ribweb stress
    IF(Rib(Ntr)%Nrb > 0) write(*,*) 'Max Rib Web Von Mises Stress: ',sigmavonMisesrwmax,'[MPa]'
    !Ribcap stress    
    IF(Rib(Ntr)%Nrb > 0) write(*,*) 'Max Rib Cap Von Mises Stress: ',sigmavonMisesrcmax,'[MPa]'
    !Stringer stress
    IF(Stringer(Ntr)%Nst > 0) write(*,*) 'Max Stringer Von Mises Stress:',sigmavonMisesstmax,'[MPa]'

end subroutine Stress 
 