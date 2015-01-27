module global
use config
implicit none

public char
 
contains
  ! Converts anything into a character
  character function char(i,r)
    ! INPUT VARIABLES
    integer,optional     :: i
    real(wp),optional    :: r
    
    ! CONVERT THE VARIABLE
    if (present(i)) then
      write(char, '(i10)')i
    else if (present(r)) then
      write(char, '(i10)')r
    end if
  end function char







end module global