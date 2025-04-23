PROGRAM RMHDS_main
!-------------------------------------------------------------------------------
!
!    Nonlinear reduced MHD code in a flux tube geometry
!
!      Hierarchy of modules (The lower should be complied earlier)
!
!        main
!         |
!        set, diag, advnc
!         |
!        exb
!         |
!        fld
!         |
!        bndry, fft, clock, tips
!         |
!        header, mpienv, math
!
!-------------------------------------------------------------------------------

  use RMHDS_header
  use RMHDS_mpienv
  use RMHDS_set,   only: set_init, set_close
  use RMHDS_clock, only: clock_timer
  use RMHDS_diag,  only: diag_cntrl
  use RMHDS_advnc, only: advnc_rkgsteps


  implicit none


  complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: omg, psi, dns

  real(kind=DP) :: time

   integer :: iflg  ! the variable "loop" is defined as a global variable 


!   integer :: mx, my


    call mpienv_init( nprocz )

    call set_init( omg, psi, dns, time )
!        print *, "# init end: rankz = ", rankz


        write( olog, * ) " # simulation is started at t = ", time


!! debug
!            if ( rankz == 0 ) then
!              write(olog,fmt="(a)") "# zz, dns in main 2"
!              mx   =-2
!              my   = 1
!              write(olog,fmt="(1p,5e15.7)") zz(-nz), dns(mx,my,-nz)
!              write(olog,*)
!            end if
!! debug


    loop  = 0

    call clock_timer( 0, iflg )
!        print *, "# clock end: rankz = ", rankz


    call diag_cntrl( omg, psi, dns, time, 0 )
!        print *, "# diag0 end: rankz = ", rankz


    do

      if ( time > tend - eps ) exit

      time   = time + dt
      loop   = loop + 1

      call advnc_rkgsteps( omg, psi, dns )

      call diag_cntrl( omg, psi, dns, time, 1 )

      call clock_timer( 1, iflg )
      if( iflg == 1 ) exit

    end do


    call diag_cntrl( omg, psi, dns, time, 2 )
        write( olog, * ) " # simulation is stopped at t = ", time


    call clock_timer( 2, iflg )


    call set_close


    call MPI_Finalize ( ierr_mpi )


  stop


END PROGRAM RMHDS_main
