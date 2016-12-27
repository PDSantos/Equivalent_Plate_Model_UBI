module SplineMaths

    implicit none
contains

    ! -------------------------------------------------------------------
    ! Routine to linearly interpolate curve between two points.
    pure function functionLinearInter(n1,x1,y1,n,x)

      ! Input variables.
      integer, intent(in) ::			n1					! nº de pontos da série
      real(8), intent(in) ::			x1(n1)				! abcissas da série
      real(8), intent(in) ::			y1(n1)				! ordenadas da série
      integer, intent(in) ::			n					! nº de pontos da nova série    
      real(8), intent(in) ::		    x(n)				! abcissa do ponto a interpolar
  
      ! Output variables.
      real(8) ::            			functionLinearInter(n)

      ! Local variables.
      integer ::       					i,j,k				! indices
 
      if(x1(1).LT.x1(2)) then
          do i=1,n
            do j=1,n1-1
              k	= j
              if(x(i).LT.x1(1)) exit
              if(x(i).GE.x1(j).AND.x(i).LT.x1(j+1)) exit
            end do
            functionLinearInter(i)	= y1(k)+(y1(k+1)-y1(k))/(x1(k+1)-x1(k))*(x(i)-x1(k))
          end do
	    else
          do i=1,n
            do j=1,n1-1
              k	= j
              if(x(i).GT.x1(1)) exit
              if(x(i).LE.x1(j).AND.x(i).GT.x1(j+1)) exit
            end do
            functionLinearInter(i)	= y1(k)+(y1(k+1)-y1(k))/(x1(k+1)-x1(k))*(x(i)-x1(k))
          end do
      end if

      return

    end function functionLinearInter

    ! -------------------------------------------------------------------
    ! Routine to linearly interpolate curve between two points.
    subroutine LinearInterpolation(n1,x1,y1,n,x,y)

      ! Input variables.
      integer, intent(in) ::			n1					! nº de pontos da série
      real(8), intent(in) ::			x1(n1)				! abcissas da série
      real(8), intent(in) ::			y1(n1)				! ordenadas da série
      integer, intent(in) ::			n					! nº de pontos da nova série    
      real(8), intent(in) ::		    x(n)				! abcissa do ponto a interpolar
  
      ! Output variables.    
      real(8), intent(out) ::			y(n)				! ordenada interpolada

      ! Local variables.
      integer ::       					i,j,k				! indices
 
      if(x1(1).LT.x1(2)) then
          do i=1,n
            do j=1,n1-1
              k	= j
              if(x(i).LT.x1(1)) exit
              if(x(i).GE.x1(j).AND.x(i).LT.x1(j+1)) exit
            end do
            y(i)	= y1(k)+(y1(k+1)-y1(k))/(x1(k+1)-x1(k))*(x(i)-x1(k))
          end do
	    else
          do i=1,n
            do j=1,n1-1
              k	= j
              if(x(i).GT.x1(1)) exit
              if(x(i).LE.x1(j).AND.x(i).GT.x1(j+1)) exit
            end do
            y(i)	= y1(k)+(y1(k+1)-y1(k))/(x1(k+1)-x1(k))*(x(i)-x1(k))
          end do
      end if

      return

    end subroutine LinearInterpolation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!XFOIL Spline!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Pure function SPLINE_function(X,S,N)
      
          integer, intent(in) ::  N
          real(4), intent(in) ::  X(N),S(N)  
          real(4) ::              SPLINE_function(N) 
      
          !Local Variables
          integer             ::  i
          real(4)             ::  A(N),B(N),C(N)
          real(4)             ::  DSM,DSP

    !-------------------------------------------------------
    !     Calculates spline coefficients for X(S).          |
    !     Zero 2nd derivative end conditions are used.      |
    !     To evaluate the spline at some value of S,        |
    !     use SEVAL and/or DEVAL.                           |
    !                                                       |
    !     S        independent variable array (input)       |
    !     X        dependent variable array   (input)       |
    !     XS       dX/dS array                (calculated)  |
    !     N        number of points           (input)       |
    !                                                       |
    !-------------------------------------------------------
          DO I=2, N-1
            DSM = S(I) - S(I-1)
            DSP = S(I+1) - S(I)
            B(I) = DSP
            A(I) = 2.0*(DSM+DSP)
            C(I) = DSM
            SPLINE_function(I) = 3.0*((X(I+1)-X(I))*DSM/DSP + (X(I)-X(I-1))*DSP/DSM)
          end do
    !---- set zero second derivative end conditions
          A(1) = 2.0
          C(1) = 1.0
          SPLINE_function(1) = 3.0*(X(2)-X(1)) / (S(2)-S(1))
          B(N) = 1.0
          A(N) = 2.0
          SPLINE_function(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
    !---- solve for derivative array XS
          CALL TRISOL(A,B,C,SPLINE_function,N)
      
    end function SPLINE_function

    SUBROUTINE SPLINE(X,XS,S,N)
      
          integer, intent(in)           ::  N
          real(4), intent(in)           ::  X(N),S(N)  
          real(4), intent(out)          ::  XS(N) 
      
          !Local Variables
          integer                       ::  i
          real(4)                       ::  A(N),B(N),C(N) !,du2(N-2)
          real(4)                       ::  DSM,DSP
    !-------------------------------------------------------
    !     Calculates spline coefficients for X(S).          |
    !     Zero 2nd derivative end conditions are used.      |
    !     To evaluate the spline at some value of S,        |
    !     use SEVAL and/or DEVAL.                           |
    !                                                       |
    !     S        independent variable array (input)       |
    !     X        dependent variable array   (input)       |
    !     XS       dX/dS array                (calculated)  |
    !     N        number of points           (input)       |
    !                                                       |
    !-------------------------------------------------------
          DO I=2, N-1
            DSM = S(I) - S(I-1)
            DSP = S(I+1) - S(I)
            B(I) = DSP
            A(I) = 2.0*(DSM+DSP)
            C(I) = DSM
            XS(I) = 3.0*((X(I+1)-X(I))*DSM/DSP + (X(I)-X(I-1))*DSP/DSM)
          end do
    !---- set zero second derivative end conditions
          A(1) = 2.0
          C(1) = 1.0
          XS(1) = 3.0*(X(2)-X(1)) / (S(2)-S(1))
          B(N) = 1.0
          A(N) = 2.0
          XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
    !---- solve for derivative array XS
          CALL TRISOL(A,B,C,XS,N)
      
          END subroutine SPLINE


      SUBROUTINE SPLIND(X,XS,S,N,XS1,XS2)
      
          integer, intent(in)           ::  N
          real(4), intent(in)           ::  X(N),S(N),XS1,XS2  
          real(4), intent(out)          ::  XS(N) 
      
          !Local Variables
          integer                       ::  i
          real(4)                       ::  A(N),B(N),C(N)
          real(4)                       ::  DSM,DSP
    !-------------------------------------------------------
    !     Calculates spline coefficients for X(S).          |
    !     Specified 1st derivative and/or usual zero 2nd    |
    !     derivative end conditions are used.               |
    !     To evaluate the spline at some value of S,        |
    !     use SEVAL and/or DEVAL.                           |
    !                                                       |
    !     S        independent variable array (input)       |
    !     X        dependent variable array   (input)       |
    !     XS       dX/dS array                (calculated)  |
    !     N        number of points           (input)       |
    !     XS1,XS2  endpoint derivatives       (input)       |
    !              If = 999.0, then usual zero second       |
    !              derivative end condition(s) are used     |
    !              If = -999.0, then zero third             |
    !              derivative end condition(s) are used     |
    !                                                       |
    !-------------------------------------------------------
          DO I=2, N-1
            DSM = S(I) - S(I-1)
            DSP = S(I+1) - S(I)
            B(I) = DSP
            A(I) = 2.0*(DSM+DSP)
            C(I) = DSM
            XS(I) = 3.0*((X(I+1)-X(I))*DSM/DSP + (X(I)-X(I-1))*DSP/DSM)
          end do
    !
          IF(XS1.EQ.999.0) THEN
    !----- set zero second derivative end condition
           A(1) = 2.0
           C(1) = 1.0
           XS(1) = 3.0*(X(2)-X(1)) / (S(2)-S(1))
          ELSE IF(XS1.EQ.-999.0) THEN
    !----- set zero third derivative end condition
           A(1) = 1.0
           C(1) = 1.0
           XS(1) = 2.0*(X(2)-X(1)) / (S(2)-S(1))
          ELSE
    !----- set specified first derivative end condition
           A(1) = 1.0
           C(1) = 0.
           XS(1) = XS1
          ENDIF
    !
          IF(XS2.EQ.999.0) THEN
           B(N) = 1.0
           A(N) = 2.0
           XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
          ELSE IF(XS2.EQ.-999.0) THEN
           B(N) = 1.0
           A(N) = 1.0
           XS(N) = 2.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
          ELSE
           A(N) = 1.0
           B(N) = 0.
           XS(N) = XS2
          ENDIF
    !
          IF(N.EQ.2 .AND. XS1.EQ.-999.0 .AND. XS2.EQ.-999.0) THEN
           B(N) = 1.0
           A(N) = 2.0
           XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
          ENDIF
    !
    !---- solve for derivative array XS
          CALL TRISOL(A,B,C,XS,N)
    !
          RETURN
          END subroutine SPLIND


          SUBROUTINE SPLINA(X,XS,S,N)
      
          integer, intent(in)           ::  N
          real(4), intent(in)           ::  X(N),S(N) 
          real(4), intent(out)          ::  XS(N) 
      
          !Local Variables
          integer                       ::  i
          real(4)                       ::  ds,XS1,XS2,dx
          LOGICAL                           LEND
    !-------------------------------------------------------
    !     Calculates spline coefficients for X(S).          |
    !     A simple averaging of adjacent segment slopes     |
    !     is used to achieve non-oscillatory curve          |
    !     End conditions are set by end segment slope       |
    !     To evaluate the spline at some value of S,        |
    !     use SEVAL and/or DEVAL.                           |
    !                                                       |
    !     S        independent variable array (input)       |
    !     X        dependent variable array   (input)       |
    !     XS       dX/dS array                (calculated)  |
    !     N        number of points           (input)       |
    !                                                       |
    !-------------------------------------------------------
    !
          LEND = .TRUE.
          DO I=1, N-1
            DS = S(I+1)-S(I)
            IF (DS.EQ.0.) THEN
              XS(I) = XS1
              LEND = .TRUE.
             ELSE
              DX = X(I+1)-X(I)
              XS2 = DX / DS
              IF (LEND) THEN
                XS(I) = XS2
                LEND = .FALSE.
               ELSE
                XS(I) = 0.5*(XS1 + XS2)
              ENDIF
            ENDIF
            XS1 = XS2
          end do
          XS(N) = XS1

          END SUBROUTINE SPLINA

          Pure SUBROUTINE TRISOL(A,B,C,D,KK)
      
          implicit none
      
          integer, intent(in)           ::  KK
          real(4), intent(in)           ::  B(KK)
          real(4), intent(inout)        ::  A(KK),C(KK),D(KK)
          !Local Variables
          integer                       ::  K,KM
    !-----------------------------------------
    !     Solves KK long, tri-diagonal system |
    !                                         |
    !             A C          D              |
    !             B A C        D              |
    !               B A .      .              |
    !                 . . C    .              |
    !                   B A    D              |
    !                                         |
    !     The righthand side D is replaced by |
    !     the solution.  A, C are destroyed.  |
    !-----------------------------------------
    !
          DO K=2, KK
            KM = K-1
            C(KM) = C(KM) / A(KM)
            D(KM) = D(KM) / A(KM)
            A(K) = A(K) - B(K)*C(KM)
            D(K) = D(K) - B(K)*D(KM)
          end do
          D(KK) = D(KK)/A(KK)
          DO K=KK-1, 1, -1
            D(K) = D(K) - C(K)*D(K+1)
          end do
          END subroutine TRISOL



    !--------------------------------------------------
    !     Calculates X(SS)                             |
    !     XS array must have been calculated by SPLINE |
    !--------------------------------------------------
          Pure Real Function SEVAL(SS,X,XS,S,N)
          implicit none
      
          integer, intent(in)       ::  N
          real(4), intent(in)       ::  X(N), XS(N), S(N)
          real(4), intent(in)       ::  SS
      
          !Local variables
          integer                   ::  i,ILOW,IMID
          real(4)                   ::  DS,T,CX1,CX2

          ILOW = 1
          I = N
    !
       10 IF(I-ILOW .LE. 1) GO TO 11
    !
          IMID = (I+ILOW)/2
          IF(SS .LT. S(IMID)) THEN
           I = IMID
          ELSE
           ILOW = IMID
          ENDIF
          GO TO 10
    !
       11 DS = S(I) - S(I-1)
          T = (SS - S(I-1)) / DS
          CX1 = DS*XS(I-1) - X(I) + X(I-1)
          CX2 = DS*XS(I)   - X(I) + X(I-1)
          SEVAL = T*X(I) + (1.0-T)*X(I-1) + (T-T*T)*((1.0-T)*CX1 - T*CX2)
          END Function SEVAL
    !--------------------------------------------------
    !     Calculates dX/dS(SS)                         |
    !     XS array must have been calculated by SPLINE |
    !--------------------------------------------------     
          Pure Real FUNCTION DEVAL(SS,X,XS,S,N)

          implicit none
      
          integer, intent(in)       ::  N
          real(4), intent(in)       ::  X(N), XS(N), S(N)
          real(4), intent(in)       ::  SS
      
          !Local variables
          integer                   ::  i,ILOW,IMID
          real(4)                   ::  DS,T,CX1,CX2

          ILOW = 1
          I = N
    !
       10 IF(I-ILOW .LE. 1) GO TO 11
    !
          IMID = (I+ILOW)/2
          IF(SS .LT. S(IMID)) THEN
           I = IMID
          ELSE
           ILOW = IMID
          ENDIF
          GO TO 10
    !
       11 DS = S(I) - S(I-1)
          T = (SS - S(I-1)) / DS
          CX1 = DS*XS(I-1) - X(I) + X(I-1)
          CX2 = DS*XS(I)   - X(I) + X(I-1)
          DEVAL = X(I) - X(I-1) + (1.-4.0*T+3.0*T*T)*CX1 + T*(3.0*T-2.)*CX2
          DEVAL = DEVAL/DS
          RETURN
          END Function DEVAL
    !--------------------------------------------------
    !     Calculates d2X/dS2(SS)                       |
    !     XS array must have been calculated by SPLINE |
    !--------------------------------------------------      
          Pure Real FUNCTION D2VAL(SS,X,XS,S,N)
          implicit none
      
          integer, intent(in)       ::  N
          real(4), intent(in)       ::  X(N), XS(N), S(N)
          real(4), intent(in)       ::  SS
      
          !Local variables
          integer                   ::  i,ILOW,IMID
          real(4)                   ::  DS,T,CX1,CX2

          ILOW = 1
          I = N
    !
       10 IF(I-ILOW .LE. 1) GO TO 11
    !
          IMID = (I+ILOW)/2
          IF(SS .LT. S(IMID)) THEN
           I = IMID
          ELSE
           ILOW = IMID
          ENDIF
          GO TO 10
    !
       11 DS = S(I) - S(I-1)
          T = (SS - S(I-1)) / DS
          CX1 = DS*XS(I-1) - X(I) + X(I-1)
          CX2 = DS*XS(I)   - X(I) + X(I-1)
          D2VAL = (6.*T-4.)*CX1 + (6.*T-2.0)*CX2
          D2VAL = D2VAL/DS**2
          RETURN
          END Function D2VAL
      
    !-----------------------------------------------------------------------------------------------
    ! Routine to evaluate spline value and derivatives up to 2nd order,given spline coefficients.
    subroutine evaluate_spline(NDATA,XDATA,YDATA,XS_Coef,IDER,N,X,Y)

      implicit none
      ! Input varables.
      integer, intent(in) ::			            ndata,ider		! no. of data points and derivative desired
      real(4), dimension(ndata), intent(in) ::	    xdata,ydata		! x,y coordinates
      real(4), dimension(ndata), intent(in) ::	    XS_Coef         ! coefficents matrix
      integer, intent(in) ::			            n   			! no. of points to interpolate
      real(4), dimension(n), intent(in) ::	        x				! x coordinate of points to be interpolated
      ! Output varables.
      real(4), dimension(n), intent(out) ::         y				! y coordinate of interpolated point
      ! Local variables.
      integer ::                                    i
      !real ::                                       SEVAL,DEVAL,D2VAL
  
      IF(IDER.EQ.0) Then
          do i=1,n
            y(i)=SEVAL(x(i),ydata,XS_COEF,xdata,ndata)
          end do
      Else IF(IDER.EQ.1) Then
          do i=1,n
            y(i)=DEVAL(x(i),ydata,XS_COEF,xdata,ndata)
          end do
      Else IF(IDER.EQ.2) Then
          do i=1,n
            y(i)=D2VAL(x(i),ydata,XS_COEF,xdata,ndata)
          end do
      Else
          write(*,*) 'Error'
          y=-1.0E9   
      End IF
  
      end subroutine evaluate_spline
  
    !-----------------------------------------------------------------------------------------------
    ! Function to evaluate spline value and derivatives up to 2nd order,given spline coefficients.
    pure function functionevaluate_spline(NDATA,XDATA,YDATA,XS_Coef,IDER,N,X)

      implicit none
      ! Input varables.
      integer, intent(in) ::			            ndata,ider		! no. of data points and derivative desired
      real(4), dimension(ndata), intent(in) ::	    xdata,ydata		! x,y coordinates
      real(4), dimension(ndata), intent(in) ::	    XS_Coef         ! coefficents matrix
      integer, intent(in) ::			            n   			! no. of points to interpolate
      real(4), dimension(n), intent(in) ::	        x				! x coordinate of points to be interpolated
      ! Output varables.
      real(4), dimension(n) ::                      functionevaluate_spline	! y coordinate of interpolated point
      ! Local variables.
      integer ::                                    i
  
      IF(IDER.EQ.0) Then
          forall(i=1:n)
            functionevaluate_spline(i)=SEVAL(x(i),ydata,XS_COEF,xdata,ndata)
          end forall
      Else IF(IDER.EQ.1) Then
          forall(i=1:n)
            functionevaluate_spline(i)=DEVAL(x(i),ydata,XS_COEF,xdata,ndata)
          end forall
      Else IF(IDER.EQ.2) Then
          forall(i=1:n)
            functionevaluate_spline(i)=D2VAL(x(i),ydata,XS_COEF,xdata,ndata)
          end forall 
      End IF
  
    end function functionevaluate_spline
    
end module SplineMaths