MODULE RMHDS_clock
!-------------------------------------------------------------------------------
!
!    Elapsed time measurements
!
!-------------------------------------------------------------------------------

  use RMHDS_header
  use RMHDS_mpienv

  implicit none

  private


  public   clock_timer


CONTAINS


!--------------------------------------
  SUBROUTINE clock_timer( isw, iflg )
!--------------------------------------

    implicit none

    integer, intent(in)  :: isw
    integer, intent(out) :: iflg

    real(kind=DP), save  :: sss0, eee0
    real(kind=DP)        :: ttotl, tother


      if( isw == 0 ) then

        if( rank == 0 ) then
          call clock_etime ( sss0 )
        end if

        call MPI_Bcast( sss0, 1, MPI_DOUBLE_PRECISION, 0, &
                      MPI_COMM_WORLD, ierr_mpi )

        iflg = 0


      else if( isw == 1 ) then

        if( rank == 0 ) then
          call clock_etime ( eee0 )
        end if

        call MPI_Bcast( eee0, 1, MPI_DOUBLE_PRECISION, 0, &
                      MPI_COMM_WORLD, ierr_mpi )

        if ( eee0-sss0 > e_limit ) then
          write( olog, * ) " # elapsed time is closing to the limit", eee0-sss0
          iflg = 1
        else
          iflg = 0
        end if


      else if( isw == 2 ) then

        if( rank == 0 ) then
          call clock_etime ( eee0 )
        end if

        call MPI_Bcast( eee0, 1, MPI_DOUBLE_PRECISION, 0, &
                      MPI_COMM_WORLD, ierr_mpi )

        ttotl   = eee0 - sss0


        write( olog, * ) ""
        write( olog, * ) " ######### elapsed time summary (sec) #########"
        write( olog, * ) "   total   = ", ttotl
        write( olog, * ) " ##############################################"

        iflg = 0


      end if

    return


  END SUBROUTINE clock_timer


!--------------------------------------
  SUBROUTINE clock_etime( ttt )
!--------------------------------------

    implicit none

    real(kind=DP) :: ttt

      ttt = MPI_Wtime()

    return

  END SUBROUTINE clock_etime


END MODULE RMHDS_clock
