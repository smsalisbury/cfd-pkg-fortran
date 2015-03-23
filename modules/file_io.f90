module file_io
use config
implicit none

private new_unit

!DATA DICTIONARY

public file_open,create_dir,dir_exists,file_close

contains
  ! This function opens a file for reading, writing, or both
  ! integer file_open(character filename, <character act='READWRITE'>, <character sta='UNKNOWN'>)
  !   filename    character(*)    The filename to read/write to.
  !   act         character(*)    Optional. Must be READ, WRITE, or READWRITE. Defaults to READWRITE.
  !   sta         character(*)    Optional. Must be OLD, NEW, SCRATCH, REPLACE, or UNKNOWN. Defaults to UNKNOWN.
  integer function file_open(filename,iact,ista)
    ! INPUT VARIABLES
    character(*)                :: filename
    character(*),optional       :: iact,ista
    
    ! INTERNAL VARIABLES
    integer                     :: i_err,file_unit
    logical                     :: file_exists
    character(10)               :: act,sta
    
    ! DEFAULTS
    file_open = -1
    file_unit=new_unit()
    if (present(iact)) then
      act = iact
    else
      act = 'READWRITE'
    end if
    
    ! CHECK IF THE FILE EXISTS
    inquire(file=filename,exist=file_exists)
    if (present(ista)) then
      sta = ista
    else
      if (file_exists) then
        sta = 'OLD'
      else
        sta = 'NEW'
      end if
    endif
    
    ! OPEN THE FILE
    ! write(*,*)'Opening '//trim(adjustl(filename))
    ! write(*,*)'Options:'
    ! write(*,*)'unit:',file_unit
    ! write(*,*)'file:',trim(adjustl(filename))
    ! write(*,*)'action:',act
    ! write(*,*)'status:',sta
    open(unit=file_unit,file=trim(adjustl(filename)),action=act,iostat=i_err,status=sta)
    
    if (i_err == 128) then
      write (*,*) 'File: ',file_unit,' does not exist '
    else if (i_err == 48) then
      write (*,*) 'Unable to open file: ',file_unit
    else if (i_err /= 0) then
      write(*,*)'Error opening file, code: ',i_err
      stop
    endif
    
    file_open = file_unit
  end function file_open
  
  ! This function closes a previously opened file
  subroutine file_close(file_unit)
    ! INPUT VARIABLES
    integer           :: file_unit
    
    ! CLOSE THE FILE
    close(file_unit)
  end subroutine file_close
  
  ! This is a function to search for an available unit.
  ! UN_MIN and UN_MAX define the range of possible units to check.
  ! Returns -1 if there are no units available
  integer function new_unit()
    ! INTERNAL VARIABLES
    integer,parameter :: UN_MIN = 10, UN_MAX = 1000
    logical           :: is_opened
    integer           :: un
    
    ! BEGIN SEARCH
    new_unit = -1
    do un=UN_MIN,UN_MAX
      inquire(unit=un,opened=is_opened)
      if (.NOT. is_opened) then
        new_unit = un
        exit
      end if
    end do    
  end function
  
  ! This function tests whether a directory exists or not
  ! Returns TRUE or FALSE correspondingly
  logical function dir_exists(directory)
    ! INPUT VARIABLES
    character(*)    :: directory
    
    ! DEFAULTS
    dir_exists = .FALSE.
    
    ! TEST
    inquire(file=trim(directory)//'/.',exist=dir_exists)
  end function dir_exists
  
  ! This subroutine makes a new directory
  ! This does nothing if the directory already exists
  subroutine create_dir(directory)
    ! INPUT VARIABLES
    character(*)    :: directory
    
    ! INTERNAL VARIABLES
    logical         :: dirExists
    
    ! TEST DIRECTORY
    dirExists = dir_exists(directory)
    
    ! MAKE DIRECTORY
    if (.NOT. dirExists) then
      call system('mkdir -p '//trim(directory))
    end if
  end subroutine create_dir
  
  ! This subroutine deletes a file
  subroutine rm_file(filename)
    ! INPUT VARIABLES
    character(*)    :: filename
    
    ! DELETE FILE
	call system('rm '//trim(filename))
  end subroutine rm_file

end module file_io
