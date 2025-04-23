MODULE RMHDS_math
!
!  Mathematical functions using SSLII library
!

  use RMHDS_header

  implicit none

  private

  public   math_random



CONTAINS


  SUBROUTINE math_random( rr )
!
!     Random number in [0,1]
!

    real(kind=DP), intent(inout), dimension(:) :: rr

    real(kind=4), allocatable, dimension(:) :: ss
    integer(kind=4), save :: iseed 
    integer :: nr, ierr

    data iseed / 211501 /

      nr = size(rr)

      allocate( ss(nr) )

!        call ranu2( iseed, ss, nr, ierr )
        call random_number( ss )

        rr(:) = ss(:)

      deallocate( ss )

    return

  END SUBROUTINE math_random


END MODULE RMHDS_math
