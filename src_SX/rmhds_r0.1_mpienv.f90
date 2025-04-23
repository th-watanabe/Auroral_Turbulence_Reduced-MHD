MODULE RMHDS_mpienv
!-------------------------------------------------------------------------------
!
!      Header and settings for using MPI
!
!
!-------------------------------------------------------------------------------

!--------------------------------------
!  Variables for MPI parallilization
!--------------------------------------
  implicit none

  public

  include "mpif.h"

  integer :: rank, nproc
  integer :: sizedouble_c, ierr_mpi
  integer, dimension(MPI_STATUS_SIZE)  :: status

  integer :: rankz

  integer :: izup, izdn

  integer ::  sendzdn, sendzup, recvzdn, recvzup


CONTAINS

!--------------------------------------
  SUBROUTINE mpienv_init( nprocz )
!--------------------------------------

  integer, intent(in) :: nprocz


!--- begin MPI settings

    call MPI_Init ( ierr_mpi )

    call MPI_Comm_rank ( MPI_COMM_WORLD, rank, ierr_mpi )

    call MPI_Comm_size ( MPI_COMM_WORLD, nproc, ierr_mpi )

    call MPI_Type_size ( MPI_DOUBLE_COMPLEX, sizedouble_c, ierr_mpi )


! --- process allocation to domain

    rankz   = rank

                             ! rank of the targets
    izup  = rankz + 1
    izdn  = rankz - 1

    if ( rankz == 0       ) izdn  = MPI_PROC_NULL
    if ( rankz == nproc-1 ) izup  = MPI_PROC_NULL

      if ( nproc /= nprocz ) then
        write( 6, * ) &
           " # proccesor assigment is invalid, nproc = ", nproc
        call MPI_Finalize ( ierr_mpi )
        stop
      end if




  END SUBROUTINE mpienv_init


END MODULE RMHDS_mpienv
