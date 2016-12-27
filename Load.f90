!--------------------------------------------------------------------------------------------
! Routine to calculate concentrated load vector.
subroutine ConcentratedLoadVector(k,P)

    use PlanformVariables, only: Ntr,xLEroot,Lamdapl,yroot,croot,b2,mcpl
    use ConcentratedLoadVariables
    use maths, only: TrialPolynomials
    use GeneralVariables, only: PrintVector

  implicit none
  
  !Input Variables
  integer, intent(in) ::	k					! number of polynomials
  !Output Variables
  real(8), intent(out) ::	P(5*k*k)			! load vector
  !Local Variables
  integer ::				i,n					! iteration variables
  real(8) ::				BB(k*k)				! Legendre polynimial multiplication vector
  real(8) ::				zeta,eta			! point coordinates
  character(60) ::			name				! vector name

  P	= 0.0D0
  if(CLoad(Ntr)%Ncl.EQ.0) return
  do n=1,CLoad(Ntr)%Ncl
   zeta	= -1.0D0+2.0D0*(CLoad(Ntr)%xcl(n)-xLEroot(Ntr)-(CLoad(Ntr)%ycl(n)-yroot(Ntr))*dtan(Lamdapl(Ntr)))/(croot(Ntr)+mcpl(Ntr)*(CLoad(Ntr)%ycl(n)-yroot(Ntr)))
   eta	= -1.0D0+2.0D0*(CLoad(Ntr)%ycl(n)-yroot(Ntr))/b2(Ntr)
   call TrialPolynomials(1,k,zeta,eta,BB)
! Pu.
    do i=1,k*k
	  P(i)	= P(i)+CLoad(Ntr)%Pclx(n)*BB(i)
    end do
! Pv.
    do i=k*k+1,2*k*k
	  P(i)	= P(i)+CLoad(Ntr)%Pcly(n)*BB(i-k*k)
    end do
! Pw.
    do i=2*k*k+1,3*k*k
	  P(i)	= P(i)+CLoad(Ntr)%Pclz(n)*BB(i-2*k*k)
    end do
! Px.
    do i=3*k*k+1,4*k*k
	  P(i)	= P(i)+CLoad(Ntr)%zcl(n)*CLoad(Ntr)%Pclx(n)*BB(i-3*k*k)
    end do
! Py.
    do i=4*k*k+1,5*k*k
	  P(i)	= P(i)+CLoad(Ntr)%zcl(n)*CLoad(Ntr)%Pcly(n)*BB(i-4*k*k)
    end do
  end do

! Write vector into a file.
  if(PrintVector) then
    name	= 'concentrated load vector'
    call WriteVector(P,5*k*k,name)
  end if

end subroutine ConcentratedLoadVector

!--------------------------------------------------------------------------------------------
! Routine to calculate actuator load vector.
subroutine ActuatorLoadVector(k,P)

    use PlanformVariables, only: Ntr,xLEroot,Lamdapl,yroot,croot,b2,mcpl
    use ActuatorLoadVariables, only: ALoad
    use GeneralVariables
    use maths, only: CoordinatesTransformation,TrialPolynomials
  implicit none

  !Input Variables
  integer, intent(in) ::	k					! number of polynomials
  !Output Variables
  real(8), intent(out) ::	P(5*k*k)			! load vector
  !Local Variables
  integer ::				i,n					! iteration variables
  real(8) ::				BB(k*k)				! Legendre polynimial multiplication vector
  real(8) ::				zeta,eta			! point coordinates
  character(60) ::			name				! vector name


  P	= 0.0D0
  if(ALoad(Ntr)%Nal.EQ.0) return
  do n=1,ALoad(Ntr)%Nal
    call CoordinatesTransformation(ALoad(Ntr)%xal0(n),ALoad(Ntr)%yal0(n),zeta,eta)
    call TrialPolynomials(1,k,zeta,eta,BB)
! Pu.
    do i=1,k*k
	  P(i)	= P(i)+ALoad(Ntr)%Palx(n)*BB(i)
    end do
! Pv.
    do i=k*k+1,2*k*k
	  P(i)	= P(i)+ALoad(Ntr)%Paly(n)*BB(i-k*k)
    end do
! Pw.
    do i=2*k*k+1,3*k*k
	  P(i)	= P(i)+ALoad(Ntr)%Palz(n)*BB(i-2*k*k)
    end do
! Px.
    do i=3*k*k+1,4*k*k
	  P(i)	= P(i)+ALoad(Ntr)%zal(n)*ALoad(Ntr)%Palx(n)*BB(i-3*k*k)
    end do
! Py.
    do i=4*k*k+1,5*k*k
	  P(i)	= P(i)+ALoad(Ntr)%zal(n)*ALoad(Ntr)%Paly(n)*BB(i-4*k*k)
    end do
  end do

! Write vector into a file.
  if(PrintVector) then
    name	= 'actuator load vector'
    call WriteVector(P,5*k*k,name)
  end if

end subroutine ActuatorLoadVector

!--------------------------------------------------------------------------------------------
! Routine to calculate overall load vector.
subroutine PlanformLoadVector(k)
  use PlanformVariables, only: Ntr
  use DeformationVariables, only: Deformation
  use GeneralVariables, only: PrintVector

  implicit none
  !Input variables
  integer, intent(in)   ::  k
  !Local variables
  real(8) ::				Pcl(5*k*k)			! concentrated load vector
  real(8) ::				Pal(5*k*k)			! actuator load vector
  character(60) ::			name				! matrix name

!  write(*,*) 'Computing load vector...'		!************************************

!  write(*,*) '    - concentrated'				!************************************
  call ConcentratedLoadVector(k,Pcl)
!  write(*,*) '    - actuator'					!************************************
  call ActuatorLoadVector(k,Pal)
!  write(*,*) '    - distributed'
!  call DistributedLoadVector(k,Pdl)
  Deformation(Ntr)%PV = Pcl + Pal !+Pdl(i)

! Write matrix into a file.
  if(PrintVector) then
    name	= 'planform load vector'
    call WriteVector(Deformation(Ntr)%PV,5*k*k,name)
  end if

end subroutine PlanformLoadVector

!--------------------------------------------------------------------------------------------
! Routine to calculate overall load vector.
subroutine OverallLoadVector(k)
    use PlanformVariables, only: Npl,Ntr
    use DeformationVariables, only: Deformation,PV_Global,q_Global
    use GeneralVariables, only: PrintVector

    implicit none
    !Input vars
    integer, intent(in) ::  k
    
    character(60) ::		name        ! matrix name
    integer ::              i           ! loop var

    forall(i=1:Npl)
        PV_Global((i-1)*5*k*k+1 : i*5*k*k ) = Deformation(i)%PV
    end forall

    ! Write matrix into a file.
    if(PrintVector) then
        name	= 'overall load vector'
        call WriteVector(PV_Global,size(PV_Global),name)
    end if

end subroutine OverallLoadVector
