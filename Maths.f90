module maths
    implicit none
    
    real(8), parameter ::       Pi	= 4.0D0*atan(1.0D0)     ! Pi number
    real(8), parameter ::		eps	= epsilon(Pi)           ! limit imposed to the DO WHILE cicle

contains
    elemental real(8) function ToHz(rad_s)
        
        real(8), intent(in) ::          rad_s
        ToHz = rad_s/(2.0D0*Pi)
    
    end function ToHz

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate the Legendre polynomials coefficients.
    subroutine LegendrePolynomialsCoefficients()

    !  use LegendreVariables
      use PolynomialCoefficients, only: k,Lcoef
      implicit none

      real(8) ::				c1(k),c2(k)
      integer ::				i,j,m				! iteration variable
      real(8) ::				DBLi

      ! Calculation of the C's coefficients.
      do i=1,k
	    DBLi	= dble(i)-2.0D0
	    if(i.EQ.1) then
		    c1(i)	= 0.0D0
		    c2(i)	= 0.0D0		!Initial values
	      else if(i.EQ.2) then
		    c1(i)	= 0.0D0
		    c2(i)	= 0.0D0		!Initial values
	      else
		    c1(i)	= (2.0D0*DBLi+1.0D0)/(DBLi+1.0D0)
		    c2(i)	= DBLi/(DBLi+1.0D0)
	    end if
      end do
      ! Coefficients' matrix initialization.

      Lcoef = 0.0D0

      ! Calculation of the coefficient matrix.
      Lcoef(1,1)	= 1.0D0
      Lcoef(1,2)	= 0.0D0		!Initial values
      Lcoef(2,1)	= 0.0D0
      Lcoef(2,2)	= 1.0D0		!Initial values
      do i=3,k
	    m = i
	    do j=1,m
	      if(dmod(dble(i-2),2.0D0).EQ.0.0D0) then
		      if(dmod(dble(j),2.0D0).EQ.0.0D0) then	!if even
			      Lcoef(i,j) = c1(i)*Lcoef(i-1,j-1)-c2(i)*Lcoef(i-2,j)
		        else	!if odd
			      Lcoef(i,j) = 0.0D0
		      end if
          else
		      if(j.EQ.1) then
			      Lcoef(i,j) = -c2(i)*Lcoef(i-2,j)
			    else if(dmod(dble(j),2.0D0).EQ.0.0D0) then	!if even
			      Lcoef(i,j)	= 0.0D0
			    else	!if odd
			      Lcoef(i,j)	= c1(i)*Lcoef(i-1,j-1)-c2(i)*Lcoef(i-2,j)
		      end if
	      end if
	    end do
      end do

    ! Print matrix.
    !  name= 'LgndrPol'
    !  call WriteMatrix(Lcoef,k,k,name)

    end subroutine LegendrePolynomialsCoefficients

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate Legendre polynomials.
    elemental real(8) function functionLegendrePolynomials(n,x)  ! value of Legendre polynomials

      use PolynomialCoefficients, only: Lcoef

      !Input variables
      real(8), intent(in) ::	    x					! variable for which B is calculated
      integer, intent(in) ::		n					! number of polynomials
      !Local variables
      integer ::    				j 				    ! iteration variable

      functionLegendrePolynomials = 0.0D0
      do j=1,n
          functionLegendrePolynomials	= functionLegendrePolynomials+Lcoef(n,j)*x**(j-1)
      end do
      
    end function functionLegendrePolynomials

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate first derivative of Legendre polynomials.
    elemental real(8) function functionLegendrePolynomialDerivatives(n,x)

      use PolynomialCoefficients, only: Lcoef

      !Input variables
      real(8), intent(in) ::	    x					! variable for which Bd is calculated
      integer, intent(in) ::		n					! number of polynomials
      !Local variables
      integer ::    				j 				! iteration variable
  
      functionLegendrePolynomialDerivatives=0.0D0
      do j=2,n
        functionLegendrePolynomialDerivatives	= functionLegendrePolynomialDerivatives+(j-1)*Lcoef(n,j)*x**(j-2)
      end do

    end function functionLegendrePolynomialDerivatives

    !--------------------------------------------------------------------------------------------
    ! Function to calculate the factorial of a number.
    elemental real(8) function factorial(n)
    
      implicit none
  
      integer, intent(in)   ::		n					! integer value
      integer ::    				i					! iteration variable

      factorial = 1.0D0
      do i=2,n
        factorial = factorial*dble(i)
      end do

    end function factorial

    !--------------------------------------------------------------------------------------------
    ! Routine to multiply Legendre polynomials.
    subroutine MultiplyLegendrePolynomials(flag,k,x,y,A)

        implicit none

        !Input variables
        integer, intent(in)   ::	k               ! number of polynomials
        real(8), intent(in)   ::	x,y             ! point coordinates
        integer, intent(in)   ::	flag            ! flag to choose calculation
											        !   1 = determine BB
											        !   2 = determine BdB
  											        !   3 = determine BBd
        !Output variables
        real(8), intent(out) ::	    A(k*k)          ! Output vector
        !Local variables
        integer ::				    i,j,l           ! iteration variables
        real(8) ::				    A_for(k,k)      ! Legendre polynomial multiplication vector
    
        l	= 0
        if(flag.EQ.1) then
            forall(i=1:k,j=1:k)
                A_for(i,j)  = functionLegendrePolynomials(i,x)*functionLegendrePolynomials(j,y)
            end forall
        else if(flag.EQ.2) then
            forall(i=1:k,j=1:k)
                A_for(i,j)  = functionLegendrePolynomialDerivatives(i,x)*functionLegendrePolynomials(j,y)
            end forall    
        else if(flag.EQ.3) then
            forall(i=1:k,j=1:k)
                A_for(i,j)  = functionLegendrePolynomials(i,x)*functionLegendrePolynomialDerivatives(j,y)
            end forall
        end if
        A = pack (A_for,.true.) 

    end subroutine MultiplyLegendrePolynomials

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate Ritz polynomials.
    Pure subroutine RitzPolynomials(n,x,y,B)

      use PolynomialCoefficients, only: nx,ny

      implicit none

      !Input variables
      integer, intent(in) ::		n					! number of terms
      real(8), intent(in) ::    	x,y					! variables for which B is calculated
      !Output variables
      real(8), intent(out) ::		B(n*n)			! Ritz polynomial derivtive
      !Local variables
      integer ::    				i					! iteration variables

      forall(i=1:n*n)
        B(i) = x**nx(i)*y**ny(i)
      end forall

    end subroutine RitzPolynomials

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate Ritz polynomial derivatives in x.
    Pure subroutine RitzPolynomialXDerivatives(n,x,y,Bdx)

      use PolynomialCoefficients, only: nx,ny

      implicit none

      !Input variables
      integer, intent(in) ::		n					! number of terms
      real(8), intent(in) ::    	x,y					! variables for which B is calculated
      !Output variables
      real(8), intent(out) ::		Bdx(n*n)			! Ritz polynomial derivtive
      !Local variables

      where(nx.EQ.0)
        Bdx	= 0.0D0
      elsewhere
        Bdx	= nx*x**(nx-1)*y**ny
      end where

    end subroutine RitzPolynomialXDerivatives

    !--------------------------------------------------------------------------------------------
    ! Routine to calculate Ritz polynomial derivatives in y.
    Pure subroutine RitzPolynomialYDerivatives(n,x,y,Bdy)

      use PolynomialCoefficients, only: nx,ny

      implicit none

      !Input variables
      integer, intent(in) ::		n					! number of terms
      real(8), intent(in) ::    	x,y					! variables for which B is calculated
      !Output variables
      real(8), intent(out) ::		Bdy(n*n)			! Ritz polynomial derivtive

      where(nx.EQ.0)
        Bdy	= 0.0D0
      elsewhere
        Bdy	= x**nx*ny*y**(ny-1)
      end where

    end subroutine RitzPolynomialYDerivatives

    !--------------------------------------------------------------------------------------------
    ! Routine to multiply Ritz polynomials.
    subroutine MultiplyRitzPolynomials(flag,n,x,y,A)

      implicit none
  
      !Input variables
      integer, intent(in) ::		flag				! flag to choose calculation
											            !   1 = determine BB
											            !   2 = determine BdB
											            !   3 = determine BBd
      integer, intent(in) ::		n					! number of polynomials
      real(8), intent(in) ::		x,y					! point coordinates
      !Output variables
      real(8), intent(out) ::		A(n*n)				! Output vector
      !Local variables
      integer ::    				i					! iteration variable
      real(8) ::    				BB(n*n),BdB(n*n),BBd(n*n)	! Ritz polynomial vector


      if(flag.EQ.1) then
          do i=1,n*n
            call RitzPolynomials(n,x,y,BB)
		    A(i)	= BB(i)
          end do
        else if(flag.EQ.2) then
          do i=1,n*n
            call RitzPolynomialXDerivatives(n,x,y,BdB)
		    A(i)	= BdB(i)
          end do
        else if(flag.EQ.3) then
          do i=1,n*n
            call RitzPolynomialYDerivatives(n,x,y,BBd)
		    A(i)	= BBd(i)
          end do
      end if

    end subroutine MultiplyRitzPolynomials

    !--------------------------------------------------------------------------------------------
    ! Routine to call trial polynomials.
    subroutine TrialPolynomials(flag,n,x,y,A)

      use IntegrationVariables, only: ik
      use PolynomialCoefficients, only: polynomial,k

      implicit none
      !Input variables
      integer, intent(in) ::		flag				! flag to choose calculation
											            !   1 = determine BB
											            !   2 = determine BdB
											            !   3 = determine BBd
      integer, intent(in) ::		n					! number of polynomials
      real(8), intent(in) ::		x,y					! point coordinates
      !Output variables
      real(8), intent(out) ::		A(n*n)				! Output vector

      if(polynomial) then
          call MultiplyLegendrePolynomials(flag,n,x,y,A)
        else
          call MultiplyRitzPolynomials(flag,n,x,y,A)
      end if

    end subroutine TrialPolynomials

    !--------------------------------------------------------------------------------------------
    ! Routine to transform (zeta,eta) coordinates in (x,y) coordinates.
    Pure subroutine AxesTranformation(flag,zeta,eta,x,y,DJ,J1)
        use PlanformVariables, only: Ntr
        use AxesTransformationVariables, only: ATranformation
  
      !Input variables
      integer, intent(in) ::				flag	    ! flag to choose calculation:
											            !   1 = calculate (x,y) given (zeta,eta)
											            !   2 = calculate Jacobian matrix and its determinant
											            !   3 = calculate inverse Jacobian matrix
      real(8), intent(in) ::				zeta,eta	! coordinates in (zeta,eta) plane
      !Output variables
      real(8), intent(out) ::				x,y			! coordinates in (x,y) plane
      real(8), intent(out), optional ::		DJ,J1(2,2)  ! jacobinan determinant and inverse Jacobian matrix
      !Local variables
      integer ::				            i			! iteration variable
      real(8) ::				            NT(4)		! transformation coefficients
      real(8) ::				            J(2,2)		! jacobian matrix elements

      if(flag.Eq.1) then
        NT(1)	= 0.25D0*(1.0D0-zeta)*(1.0D0-eta)
        NT(2)	= 0.25D0*(1.0D0+zeta)*(1.0D0-eta)
        NT(3)	= 0.25D0*(1.0D0+zeta)*(1.0D0+eta)
        NT(4)	= 0.25D0*(1.0D0-zeta)*(1.0D0+eta)
	    x		= 0.0D0
	    y		= 0.0D0
	    do i=1,4
	      x		= x+NT(i)*ATranformation(Ntr)%xT(i)
	      y		= y+NT(i)*ATranformation(Ntr)%yT(i)
	    end do
      end if
      if(flag.Eq.2.OR.flag.EQ.3) then
        J(1,1)	= 0.25D0*(-(1.0D0-eta)*(ATranformation(Ntr)%xT(1)-ATranformation(Ntr)%xT(2))+(1.0D0+eta)*(ATranformation(Ntr)%xT(3)-ATranformation(Ntr)%xT(4)))
        J(1,2)	= 0.25D0*(-(1.0D0-zeta)*(ATranformation(Ntr)%xT(1)-ATranformation(Ntr)%xT(4))-(1.0D0+zeta)*(ATranformation(Ntr)%xT(2)-ATranformation(Ntr)%xT(3)))
        J(2,1)	= 0.25D0*(-(1.0D0-eta)*(ATranformation(Ntr)%yT(1)-ATranformation(Ntr)%yT(2))+(1.0D0+eta)*(ATranformation(Ntr)%yT(3)-ATranformation(Ntr)%yT(4)))
        J(2,2)	= 0.25D0*(-(1.0D0-zeta)*(ATranformation(Ntr)%yT(1)-ATranformation(Ntr)%yT(4))-(1.0D0+zeta)*(ATranformation(Ntr)%yT(2)-ATranformation(Ntr)%yT(3)))
	    DJ		= J(1,1)*J(2,2)-J(1,2)*J(2,1)
      end if
      if(flag.EQ.3) then
        J1(1,1)	= J(2,2)/DJ
        J1(1,2)	= -J(1,2)/DJ
        J1(2,1)	= -J(2,1)/DJ
        J1(2,2)	= J(1,1)/DJ
      end if

    end subroutine AxesTranformation
    
    !--------------------------------------------------------------------------------------------
    ! Routine to transform (x,y) coordinates to (zeta,eta) coordinates.
    Pure subroutine CoordinatesTransformation(x,y,zeta,eta)

      use PlanformVariables, only: Ntr,xLEroot,Lamdapl,yroot,croot,mcpl,b2
  
      !Input variables
      real(8), intent(in) ::				x,y					! coordinates in (x,y) plane
      !Output variables
      real(8), intent(out) ::				zeta,eta			! coordinates in (zeta,eta) plane

      zeta	= -1.0D0+2.0D0*(x-xLEroot(Ntr)-(y-yroot(Ntr))*dtan(Lamdapl(Ntr)))/(croot(Ntr)+mcpl(Ntr)*(y-yroot(Ntr)))
      eta	= -1.0D0+2.0D0*(y-yroot(Ntr))/b2(Ntr)

    end subroutine CoordinatesTransformation

    !--------------------------------------------------------------------------------------------
    ! Routine to multiply one vector by its tranposed vector (nx1 by 1xn).
    Pure subroutine MultiplyVectorByTransposedVector(x,y)

      implicit none
      !Input variables
      real(8), dimension(:), intent(in) ::		            x			! vector
      !Output variables
      real(8), dimension(size(x),size(x)), intent(out) ::   y			! matrix k by k
      !Local variables
      integer ::    				i,j					! iteration variables

        forall(i=1:size(x),j=1:size(x))
            y(i,j)	= x(i)*x(j)
        end forall

    end subroutine MultiplyVectorByTransposedVector

    !--------------------------------------------------------------------------------------------
    ! Routine to multiply one vector by its tranposed vector (nx1 by 1xn).
    Pure subroutine MultiplyVectorByTransposedVectorC(x1,x2,y)

      implicit none
      !Input variables
      real(8), dimension(:), intent(in) ::		            x1,x2		! vectors
      !Output variables
      real(8), dimension(size(x1),size(x1)), intent(out) ::   y			! matrix k by k
      !Local variables
      integer ::    				i,j					! iteration variables

        forall(i=1:size(x1),j=1:size(x2))
            y(i,j)	= x1(i)*x2(j)
        end forall

    end subroutine MultiplyVectorByTransposedVectorC   
    
    !--------------------------------------------------------------------------------------------
    ! Routine to multiply one vector by another vector (nx1 by 1xn).
    Pure subroutine MultiplyVectorByVector(x1,x2,y)

      implicit none
      !Input variables
      real(8), dimension(:), intent(in) ::		                x1,x2		! vector
      !Output variables
      real(8), dimension(size(x1),size(x1)), intent(out) ::		y			! matrix k by k
      !Local variables
      integer ::    				i,j				! iteration variables

      forall(i=1:size(x1),j=1:size(x1))
         y(i,j)	= x1(i)*x2(j)
      end forall

    end subroutine MultiplyVectorByVector

    !--------------------------------------------------------------------------------------------
    ! Routine to add two vectors each multiplied by a constant.
    Pure subroutine AddTwoVectors(k1,x1,k2,x2,y)

      implicit none
      !Input variables
      real(8), dimension(:), intent(in) ::		        x1,x2		! vector
      real(8), intent(in) ::		                    k1,k2       ! constants
      !Output variables
      real(8), dimension(size(x1)), intent(out) ::		y           ! matrix k by k
      !Local variables
      integer ::    				i				! iteration variables

      forall(i=1:size(x1))
        y(i)	= k1*x1(i)+k2*x2(i)
      end forall

    end subroutine AddTwoVectors

    !--------------------------------------------------------------------------------------------
    ! Routine to add two matrices each multiplied by a constant.
    Pure subroutine AddTwoMatrices(k1,x1,k2,x2,y)

      implicit none
      !Input variables
      real(8), dimension(:,:), intent(in) ::		        x1,x2		! matrices
      real(8), intent(in) ::		                        k1,k2       ! constants
      !Output variables
      real(8), dimension(size(x1),size(x1)), intent(out) :: y           ! resultant matrix
      !Local variables
      integer ::    				i,j					! iteration variables

      forall(i=1:size(x1),j=1:size(x1))
        y(i,j)	= k1*x1(i,j)+k2*x2(i,j)
      end forall

    end subroutine AddTwoMatrices

    !--------------------------------------------------------------------------------------------
    ! Routine to add four matrices each multiplied by a constant.
    Pure subroutine AddFourMatrices(k1,x1,k2,x2,k3,x3,k4,x4,y)

      implicit none
      !Input variables
      real(8), dimension(:,:), intent(in) ::		x1,x2,x3,x4		! matrices
      real(8), intent(in) ::		                k1,k2,k3,k4     ! constants
      !Output variables
      real(8), dimension(size(x1),size(x1)), intent(out) ::		y   ! resultant matrix
      !Local variables
      integer ::    				i,j					! iteration variables

      forall(i=1:size(x1),j=1:size(x1))
        y(i,j)	= k1*x1(i,j)+k2*x2(i,j)+k3*x3(i,j)+k4*x4(i,j)
      end forall

    end subroutine AddFourMatrices

    !------------------------------------------------------------------------------
    !The algorithm receives as input:
    !    * n - a required number of nodes.
    !The algorithm returns:
    !    * x - array of nodes. Its index ranges within 0 and n-1.
    !    * w - array of weighting coefficients. Its index ranges between 0 and n-1.
    subroutine QuadratureWeights()

      use PolynomialCoefficients, only: nodes,weights
      use IntegrationVariables, only: ik
      use PolynomialCoefficients, only: k
      implicit none
      !Local vars
      integer ::				n					! number of sampling points
      integer ::				i,j,Nterms
      real(8) ::				r,r1,p1,p2,p3,dp3,aux1,aux2

      ! Variables initialization.
      n		= ik
      r		= 0.0D0
      r1	= 0.0D0
      p1	= 0.0D0
      p2	= 0.0D0
      p3	= 0.0D0
      dp3	= 0.0D0

      Nterms	= anint((dble(n)+1)/2)-1
      do i=0+1,Nterms+1
	    r		= dcos(PI*(4.0D0*dfloat(i-1)+3.0D0)/(4.0D0*dfloat(n)+2.0D0))
	    aux1	= dabs(r1-r)
	    aux2	= eps*(1.0D0+dabs(r))*100.0D0
	    do while(aux1.GE.aux2)
	      p2	= 0.0D0
	      p3	= 1.0D0
	      do j=1,n
		    p1	= p2
		    p2	= p3
		    p3	= (dfloat((2*(j-1)+1))*r*p2-(dfloat(j-1))*p1)/(dfloat(j))
	      END DO
	      dp3	= dfloat(n)*(r*p3-p2)/(r*r-1.0D0);
	      r1	= r
	      r		= r-p3/dp3
	      aux1	= dabs(r1-r)
	      aux2	= eps*((1.0D0+dabs(r))*100.0D0)
	    end do
	    nodes(i)		= r
	    nodes(n-i+1)	= -r
	    weights(i)		= 2.0D0/((1.0D0-r*r)*dp3*dp3)
	    weights(n-i+1)	= 2.0D0/((1.0D0-r*r)*dp3*dp3)
      end do

    end subroutine QuadratureWeights

    !------------------------------------------------------------------------------
    Pure subroutine GaussQuadratureWeightsSamplingPoints(A,x)

        use PolynomialCoefficients, only: nodes,weights

        implicit none

        real(8),dimension(:), intent(out) ::		A,x			! weights and sampling points

        x	= nodes
        A	= weights
    end subroutine GaussQuadratureWeightsSamplingPoints

    ! --------------------------------------------------------------------------
    ! Subroutine to reverse array data set - data is kind=4 
      Pure subroutine ReverseData(x,xrev)

      implicit none

      ! Parameters.

      ! Input varables.
      real, dimension(:), intent(in) ::	            x			    ! original data array
      ! Output varables.
      real, dimension(size(x)), intent(out) ::      xrev			! reversed data array
      ! Local variables.
      integer ::                                    i,k 
      real ::                                       xcopy(size(x))  ! copy of original data array
  
      k = size(x)
      xcopy=x
      do i=1,size(x)
        xrev(i)=xcopy(k)
	    k=k-1
      end do

    end subroutine ReverseData
    ! --------------------------------------------------------------------------
    ! Subroutine to reverse array data set - data is kind=8 
    Pure subroutine dReverseData(x,xrev)

      implicit none

      ! Parameters.

      ! Input varables.
      real(8), dimension(:), intent(in) ::	        x			    ! original data array
      ! Output varables.
      real(8), dimension(size(x)), intent(out) ::   xrev			! reversed data array
      ! Local variables.
      integer ::                                    i,k 
      real(8) ::                                    xcopy(size(x))  ! copy of original data array
  
      k = size(x)
      xcopy=x
      do i=1,size(x)
        xrev(i)=xcopy(k)
	    k=k-1
      end do

    end subroutine dReverseData

    subroutine PERMUTATION (X,IPERMU,IPATH)
        ! Variables declaration.
        implicit none

        ! Parameters.
  
        ! Input variables.
        real(8), dimension(:), intent(inout) ::     X 			! y coordinates
        integer, dimension(size(X)), intent(in) ::  IPERMU    	! x interval
        integer, intent(in) ::                      IPATH       ! x interval
 
        ! Input/output variables.

        ! Output varables.
        real(8), dimension(size(X)) ::               XPERMU     ! integral

        IF(IPATH == 1) Then
            XPERMU=X(IPERMU)
        Else IF(IPATH == 2) Then
            XPERMU(IPERMU)=X
        End IF
        X=XPERMU    	
        
    end subroutine PERMUTATION
    
end module maths