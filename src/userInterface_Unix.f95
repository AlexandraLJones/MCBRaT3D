! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 31 $, $Date: 2010-06-14 23:12:00 +0100 (Mon, 14 Jun 2010) $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Code/userInterface_Unix.f95 $
module UserInterface
  ! $Revision: 31 $, $Date: 2010-06-14 23:12:00 +0100 (Mon, 14 Jun 2010) $
  ! $Name: Cornish-Gilliflower $
  use ErrorMessages
  use MultipleProcesses
  implicit none
  private
  public :: printStatus, getOneArgument 
  
contains
! ---------------------------------  
  subroutine printStatus(status)
    type(ErrorMessage), intent(inout) :: status
    
    logical :: shouldStop

    shouldStop = .false. 
    if(stateIsWarning(status)) then
      if(len_trim(getCurrentMessage(status)) > 0) then
        print *, trim(getCurrentMessage(status)) 
      else 
        print *, "Status is warning"
      end if
    end if
    if(stateIsFailure(status)) then
      shouldStop = .true. 
      if(len_trim(getCurrentMessage(status)) > 0) then
        print *, trim(getCurrentMessage(status)) 
      else 
        print *, "Status is Failure"
      end if
    end if
    
    if(.not. stateIsSuccess(status)) then 
      print *, "History:"
      call firstMessage(status)
      historyLoop: do
        if (.not. moreMessagesExist(status)) exit historyLoop
        if(len_trim(getCurrentMessage(status)) > 0) print *, "  ", trim(getCurrentMessage(status))
        call nextMessage(status)
      end do historyLoop
    end if

    if(shouldStop) stop
  end subroutine printStatus
! ---------------------------------  
!  function GetOneArgument(message)
!    character (len = *),&
!      optional, intent(in) :: message
!    character(len = 256)   :: GetOneArgument
!    !
!    ! Gets the first argument from the command line
!    !   Use this version on systems that do not support iarg and getarc
!    !
!   
!    if(present(message)) print *, message
!    read *, GetOneArgument
!    print *, 'Using value ' // trim(GetOneArgument)
!  end function GetOneArgument
! ---------------------------------  
  function GetOneArgument(message)
	character (len = *), optional, intent(in) :: message
	character(len = 256)   :: GetOneArgument

    ! Gets the first argument from the command line
    !   Use this version with GNU
    !

	CALL get_command_argument(1, GetOneArgument)
	print *, 'Using value ' // trim(GetOneArgument)
  end function GetOneArgument
!----------------------------------
!  function GetOneArgument(message)
!    character (len = *),&
!      optional, intent(in) :: message
!    character(len = 256)   :: GetOneArgument
!    integer                :: ilen, m, ierr
!
!    EXTERNAL PXFGETARG
!    m = 1  ! retrieve 1st commandline argument
!
 !   !
 !   ! Use this version for cray compilers
 !   !
 !    if(present(message)) print *, message
 !    CALL PXFGETARG(m, GetOneArgument, ilen, ierr)
!     if(ilen < 1 .and. MasterProc) then
!        print *, "No file name supplied"
!        stop
!     end if
!    if(MasterProc) print *, 'Using value ' // trim(GetOneArgument)
!  end function GetOneArgument
!---------------------------------
!  function GetOneArgument(message)
!    character (len = *),&
!      optional, intent(in) :: message
!    character(len = 256)   :: GetOneArgument
!    !
!    ! Gets the first argument from the command line
!    !   Use this version on systems that support  getarg
!    !
!    interface 
!      subroutine getarg(pos, value)
!        integer,            intent(in)  :: pos
!        character(len = *), intent(out) :: value
!      end subroutine getarg
!    end interface
!    
!    if(present(message)) print *, message
!    if(iargc() < 1 .and. MasterProc) then 
!      print *, "No file name supplied." 
!      stop
!    end if 
!    call getarg(1, GetOneArgument)
!    if(MasterProc) print *, 'Using value ' // trim(GetOneArgument)
!  end function GetOneArgument
! ---------------------------------  

end module UserInterface
