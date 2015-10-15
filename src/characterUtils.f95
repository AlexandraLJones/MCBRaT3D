! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Code/characterUtils.f95 $
! touch

module CharacterUtils
  implicit none
  private
  
  integer :: maxStringLength = 25
  public :: CharToInt, IntToChar, CharToReal

  interface IntToChar
    module procedure Int4ToChar, Int8ToChar
  end interface

contains
  elemental function CharToInt(inputString)
    character(len = *), intent( in) :: inputString
    integer                         :: CharToInt
    !
    ! Reads a character string into an integer variable 
    !    No internal error checking is done - if an incorrect string
    !    is passed in the entire program will fail with a Fortran runtime error. 
    
    ! Local variables
    character(len = maxStringLength) tempString

    ! Remove trailing blanks, which might be read as 0s on some processors
    ! Use list-directed read. 
    tempString = trim(inputString)
    read(tempString, *) CharToInt
    
  end function CharToInt
  ! --------------------------------------------------------------------------------
  elemental function Int4ToChar(integerValue)
    integer(4),              intent( in) :: integerValue
    character(len = maxStringLength)  :: Int4ToChar
    !
    !   Creates the character representation of an integer.  
    !  

    write(Int4ToChar, '(I10)') integerValue
    Int4ToChar = AdjustL(Int4ToChar)
  end function Int4ToChar
  ! --------------------------------------------------------------------------------
  elemental function Int8ToChar(integerValue)
    integer(8),              intent( in) :: integerValue
    character(len = maxStringLength)  :: Int8ToChar
    !
    !   Creates the character representation of an integer.
    !

    write(Int8ToChar, '(I19)') integerValue
    Int8ToChar = AdjustL(Int8ToChar)
  end function Int8ToChar

  ! --------------------------------------------------------------------------------
  elemental function CharToReal(inputString)
    character(len = *), intent( in) :: inputString
    real                            :: CharToReal
    !
    !   Reads a character string into an real variable 
    !    No internal error checking is done - if an incorrect string
    !      is passed in the entire program will fail with a Fortran runtime error. 
    !
    ! !END
    
    ! Local variables
    character(len = maxStringLength) tempString

    ! Remove trailing blanks, which might be read as 0s on some processors
    ! Use list-directed read. 
    !
    tempString = trim(AdjustL(inputString))
    read(tempString, *) CharToReal
  end function CharToReal
  
end module CharacterUtils
