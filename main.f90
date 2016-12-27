! *******************************************************************
! Program to evaluate structural performance of a wing based on the "Equivalent Plate model"
! Based on the work of Kapania and Liu [1], in which a first order shear deformation theory (FSDT) was applied to 
! create an EPM of a wing. 
!
! by
! Pedro V. Gamboa, BEng, MSc, PhD and
! Pedro D.R. Santos, BEng, MSc
! Department of Aerospace Sciences
! University of Beira Interior
! Portugal
!
! May 2012

! [1] R. K. Kapania and Y. Liu, "Static and Vibration Analyses of General Wing Structures using
! Equivalent-Plate Models," AIAA JOURNAL, vol. 38, pp. 1269–1277, July 2000.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>

! *******************************************************************
program Main
    !DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
    USE IFPORT, only: SLEEPQQ,ETIME
    USE IFCORE, only: PEEKCHARQQ 
    USE FilePathModule, only: output_dir
    use write_routines, only: export_displacement_line,export_test_line
    use ReadInputGeometry, only: ReadInputData,InputCommonData

    implicit none
  
    ! Variables declaration.
    integer ::                  choice			!analysis type:  
								!1-All 
								!2-Mass 
								!3-SectionCentroid 
								!4-Deformation 
								!5-Stresses and Strains 
								!6-Frequencies
    real ::                     TIME22,TA(2),SSS	!CPU time counter variables
    integer ::                  HHH,MMM !hour,minute and second of time spent on calculations
    LOGICAL(4)                  pressed / .FALSE. /,append_on
    character(500) ::           comp_file
    real(8) ::                  eta,zeta
    character(30) ::            arg
    
    !Main Program Start
    Call SetTitle()
    ! Body of close_handle
    Call SetupExit()
    !Write Heading
    write(*,1030)
    
    !Read Problem conditions
    Call ReadInputData(choice) !
    !Read common data Geometry
    call InputCommonData(choice)
    !Read all the geometry
    call ReadOverallGeometry()
   
    !Analyse wing
    call Analysis(choice)
    
    !!Write specific outputs
    append_on=.False.
    comp_file=trim(adjustl(output_dir))//"tip_deflection.dat"
    eta=0.93D0
    Call export_displacement_line(tecplot_file=comp_file,tecplot_unit=1,title="Tip Deflection",variables="X,Y",zone_in=" ",append_on=append_on,n=50,eta_in=eta,factor_u=100.0D0,factor_v=100.0D0)
    comp_file=trim(adjustl(output_dir))//"le_deflection.dat"
    !zeta=-1.0D0
    zeta = 0.0D0
    Call export_displacement_line(tecplot_file=comp_file,tecplot_unit=1,title="LE Deflection",variables="X,Y",zone_in=" ",append_on=append_on,n=50,zeta_in=zeta,factor_u=10.0D0,factor_v=10.0D0)
    comp_file=trim(adjustl(output_dir))//"35_deflection.dat"
    zeta=-0.75D0
    Call export_displacement_line(tecplot_file=comp_file,tecplot_unit=1,title="35% Deflection",variables="X,Y",zone_in=" ",append_on=append_on,n=50,zeta_in=zeta,factor_u=10.0D0,factor_v=10.0D0)
   
    comp_file=trim(adjustl(output_dir))//"le_z_coordinate.dat"
    zeta=-1.0D0
    !Call export_test_line(tecplot_file=comp_file,tecplot_unit=1,title="Z coordinate",variables="y,z",zone_in=" ",append_on=append_on,n=50,zeta_in=zeta)
    
1212 continue    
    TIME22= ETIME(TA)
    HHH=TIME22/3600.0
    MMM=(TIME22/3600.0-HHH)*3600.0/60.0
    SSS=((TIME22/3600.0-HHH)*3600.0/60.0-MMM)*60.0
    write(*,'(A)') '-----'
    write(*,'(A,I3,A,I3,A,F6.2,A)') 'Program has used:',HHH,'h',MMM,'m',SSS,'s'
    write(*,'(A)') 'Press any key to close ...'
    pressed=.False.
    DO WHILE (.NOT. pressed)   
        pressed = PEEKCHARQQ ( )
        CALL SLEEPQQ (200) 
    END DO

!Heading format
1030 FORMAT('-----', &
        /,'EPM - Equivalent Plate - Structural Analysis Program V2.0 (May 2012)', &
		/,'Authors: Pedro V. Gamboa and Pedro D.R. Santos', &
		/,'Department of Aerospace Sciences', &
		/,'University of Beira Interior', &
		/,'Portugal', &
        /,' ')

end program Main


SUBROUTINE SetupExit()

    USE KERNEL32, ONLY: SetConsoleCtrlHandler, BOOL, TRUE
    implicit none
        
    INTEGER(BOOL), external :: CloseWindow
    !DEC$ ATTRIBUTES STDCALL :: CloseWindow        
    INTEGER(BOOL) :: ff
    ff = SetConsoleCtrlHandler(loc(CloseWindow), TRUE)
    if (ff == 0)  write(*,*) "SetConsoleCtrlHandler failed"
    continue
END SUBROUTINE SetupExit

SUBROUTINE SetTitle()

    USE KERNEL32, ONLY: SetConsoleTitle, BOOL, TRUE
    implicit none
    INTEGER(BOOL) :: statConsole
        
    statConsole = SetConsoleTitle("Equivalent Plate - Structure Analysis Program"C) !Change console title (C API)
    if (statConsole == 0)  write(*,*) "SetConsoleTitle failed"
    !handle_win=FindWindow_G1(NULL,"Equivalent Plate - Structure Analysis Program"C)
END SUBROUTINE SetTitle

INTEGER(BOOL) FUNCTION CloseWindow(number)
!DEC$ ATTRIBUTES STDCALL :: CloseWindow

    USE KERNEL32, ONLY: BOOL, DWORD, TRUE, FALSE
    implicit none
    INTEGER(DWORD) :: number

    SELECT CASE (number)
        CASE (0) ! Ctrl c
            PRINT "('Ctrl-C pressed. Exiting...!')"
            CloseWindow = FALSE
            STOP
        CASE (1) ! Ctrl break
            PRINT "('Ctrl-Break pressed. Exiting...!')"
            CloseWindow = TRUE
            STOP
        CASE (2) ! Console Window close
            PRINT "('Exiting...!')"
            CloseWindow = TRUE
            STOP
        CASE DEFAULT ! 5 for log off, 6 for shutdown -
            PRINT "('Exiting...!')"
            CloseWindow = FALSE
            STOP
    END SELECT
END FUNCTION CloseWindow 