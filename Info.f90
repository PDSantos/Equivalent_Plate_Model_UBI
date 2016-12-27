! This source code contains routines to give miscellaneous information.

!-----------------------------------------------------------------------------------------------
! Routine to determine program execution time.
subroutine CPUTime()

  !use PortLib
  USE IFPORT
  implicit none

  real  	seconds,TA(2)
  real		hour,min,sec,dec
 
  seconds	= ETIME(TA)
!  write(*,'(1X,A17,F8.2,A21)') 'Program has used ',seconds,' seconds of CPU time.'
  hour		= ifix(seconds/3600.0)
  min		= ifix((seconds-hour*3600.0)/60.0)
  sec		= ifix(seconds-hour*3600.0-min*60.0)
  dec		= (seconds-hour*3600.0-min*60.0-sec)*100
!  write(*,'(1X,A17,I2.2,A2,I2.2,A4,F5.2,A14)') 'Program has used ',int(hour),'h ',int(min),'min ',sec,'s of CPU time.'
!  write(*,'(1X,A17,I2.2,A1,I2.2,A1,F5.2,A12)') 'Program has used ',int(hour),':',int(min),':',sec,' of CPU time.'
  write(*,'(1X,A17,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A12)')	&
&			  'Program has used ',int(hour),':',int(min),':',int(sec),'.',int(dec),' of CPU time.'
! write(*,'(A15,F8.2,A26,F8.2,A24)') 'This includes ',TA(1),' seconds of user time and ',	&
! &								TA(2),' seconds of system time.'

  open(1,FILE="optimization_parameters.txt",ACCESS='APPEND')
  write(1,'(1X,A17,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A12)')	&
&			  'Program has used ',int(hour),':',int(min),':',int(sec),'.',int(dec),' of CPU time.'
  close(1)

end subroutine CPUTime

!-----------------------------------------------------------------------------------------------
! Routine to determine system time.
subroutine SystemTime()

  USE IFPORT

  implicit none

!   real			hour,min,sec
!   integer(4)	int_time
   character*8	char_time
 
  call TIME(char_time)
  write(*,'(1X,A12,A8)') 'System time ',char_time
!  int_time = TIME()
!  hour		= ifix(seconds/3600.0)
!  min		= ifix((seconds-hour*3600.0)/60.0)
!  sec		= seconds-hour*3600.0-min*60.0
!  write(*,'(1X,A17,I2.2,A2,I2.2,A4,F5.2,A14)') 'Program has used ',int(hour),'h ',int(min),'min ',sec,'s of CPU time.'

end subroutine SystemTime
