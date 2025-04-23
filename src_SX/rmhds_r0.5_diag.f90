MODULE RMHDS_diag
!-------------------------------------------------------------------------------
!
!    Diagnostics including entropy transfer analysis
!
!-------------------------------------------------------------------------------

  use RMHDS_header
  use RMHDS_mpienv
  use RMHDS_set,    only: set_close
  use RMHDS_fld,    only: fld_calfld
  use RMHDS_intgrl, only: intgrl_zeta
  use RMHDS_pssn,   only: pssn_derivative

  implicit none

  private

  integer, parameter :: nk = min( nx, ny )

  public   diag_cntrl


CONTAINS


!--------------------------------------
  SUBROUTINE diag_cntrl( omg, psi, dns, time, id )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: psi

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: omg, dns

    real(kind=DP), intent(in) :: time

    integer, intent(in) :: id

    complex(kind=DP), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr

    real(kind=DP), save :: tout_mag, tout_ion, tout_eng

    integer :: mx, my, iz


      call fld_calfld( omg, psi, dns, phi, cpr )


      call dt_check2 ( time, phi, psi )


      if( id == 0 ) then

        if ( dtout_mag /= 0._DP ) then 
          call wrt ( omg, psi, dns, phi, cpr, time, 0 )
        end if

! debug
!            mx   = nx*3/4
!            my   = ny*3/4
            mx   =-2
            my   = 1
            write(olog,fmt="(a,1p,e15.7)") "# time = ", time
            write(olog,fmt="(a,2i5)")      "# mx, my = ", mx, my
            write(olog,fmt="(a)")          "# zz, omg, psi "
            do iz = -nz, nz-1
              write(olog,fmt="(1p,5e15.7)")  zz(iz), omg(mx,my,iz), psi(mx,my,iz)
            end do
              write(olog,*)
! debug

!! debug
!            if ( rankz == 0 ) then
!              write(olog,fmt="(a)") "# zz, dns in diag 1"
!              mx   =-2
!              my   = 1
!              write(olog,fmt="(1p,5e15.7)") zz(-nz), dns(mx,my,-nz)
!              write(olog,*)
!            end if
!! debug


        if ( dtout_ion /= 0._DP ) then 
          call wrt ( omg, psi, dns, phi, cpr, time, 1 )
        end if  
    
        call wrt ( omg, psi, dns, phi, cpr, time, 2 )

        tout_mag  = ( int( ( time + eps )/dtout_mag ) + 1 ) * dtout_mag
        tout_ion  = ( int( ( time + eps )/dtout_ion ) + 1 ) * dtout_ion
        tout_eng  = ( int( ( time + eps )/dtout_eng ) + 1 ) * dtout_eng


      else if( id == 1 ) then

        if ( time >= tout_mag - eps ) then
          tout_mag   = tout_mag + dtout_mag
          if ( dtout_mag /= 0._DP ) then 
            call wrt ( omg, psi, dns, phi, cpr, time, 0 )
            write( olog, * ) &
              " # Magnetospheric data output at time = ", time

! debug
!            mx   = nx*3/4
!            my   = ny*3/4
            mx   =-2
            my   = 1
            write(olog,fmt="(a,1p,e15.7)") "# time = ", time
            write(olog,fmt="(a,2i5)")      "# mx, my = ", mx, my
            write(olog,fmt="(a)")          "# zz, omg, psi "
            do iz = -nz, nz-1
              write(olog,fmt="(1p,5e15.7)")  zz(iz), omg(mx,my,iz), psi(mx,my,iz)
            end do
              write(olog,*)
! debug

          end if 
!          call contnu ( ff, time )
        end if

        if ( time >= tout_ion - eps ) then
          tout_ion   = tout_ion + dtout_ion
          if ( dtout_ion /= 0._DP ) then 
            call wrt ( omg, psi, dns, phi, cpr, time, 1 )
            write( olog, * ) &
               " # Ionospheric data output at time = ", time
          end if

        end if

! ---

!! debug
!            if ( rankz == 0 ) then
!              write(olog,fmt="(a)") "# zz, dns in diag 2"
!              mx   =-2
!              my   = 1
!              write(olog,fmt="(1p,5e15.7)") zz(-nz), dns(mx,my,-nz)
!              write(olog,*)
!            end if
!! debug

        if ( time >= tout_eng - eps ) then
          tout_eng   = tout_eng + dtout_eng
          call wrt ( omg, psi, dns, phi, cpr, time, 2 )
        end if


      else if( id == 2 ) then

          call contnu ( omg, psi, dns, time )

      end if


  END SUBROUTINE diag_cntrl


!--------------------------------------
  SUBROUTINE contnu ( omg, psi, dns, time )
!--------------------------------------

    implicit none

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: omg, psi, dns

    real(kind=DP), intent(in) :: time


    integer ::mx, my, iz


      rewind ocnt
      write( unit=ocnt ) time, dt, omg, psi, dns


  END SUBROUTINE contnu



!--------------------------------------
  SUBROUTINE wrt ( omg, psi, dns, phi, cpr, time, id )
!--------------------------------------

! --- arguments

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: omg, psi, dns

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr

    real(kind=DP), intent(in) :: time

    integer, intent(in) :: id

    real(kind=DP), dimension(0:nx) :: phi_mode_x, psi_mode_x
    real(kind=DP), dimension(0:ny) :: phi_mode_y, psi_mode_y
    real(kind=DP), dimension(0:nk) :: phi_mode_k, psi_mode_k

    real(kind=DP) :: phi_total_e, psi_total_e

    real(kind=DP), dimension(0:nx) :: idns_mode_x, icpr_mode_x, iomg_mode_x
    real(kind=DP), dimension(0:ny) :: idns_mode_y, icpr_mode_y, iomg_mode_y
    real(kind=DP), dimension(0:nk) :: idns_mode_k, icpr_mode_k, iomg_mode_k
 
    real(kind=DP) :: idns_total_e, icpr_total_e, iomg_total_e

    integer        ::  mx, my, iz


      if( id == 0 ) then

          write( unit=omag ) time, omg, cpr 

      else if( id == 1 ) then

        if ( rankz == 0 ) then
          write( unit=oion ) time, dns(:,:,-nz), omg(:,:,-nz), cpr(:,:,-nz)
        end if

      else if( id == 2 ) then

        call mode_energy ( phi, phi_mode_x, phi_mode_y, phi_mode_k, phi_total_e )
        call mode_energy ( psi, psi_mode_x, psi_mode_y, psi_mode_k, psi_total_e )

!!!        call dt_check ( time, phi_total_e, psi_total_e )
!!!        call dt_check2 ( time, phi, psi )

!! debug
!            if ( rankz == 0 ) then
!              write(olog,fmt="(a)") "# zz, dns in wrt"
!              mx   =-2
!              my   = 1
!              write(olog,fmt="(1p,5e15.7)") zz(-nz), dns(mx,my,-nz)
!              write(olog,*)
!            end if
!! debug

        if ( rankz == 0 ) then
          call imode_energy ( dns(:,:,-nz), idns_mode_x, idns_mode_y, idns_mode_k, idns_total_e )
          call imode_energy ( cpr(:,:,-nz), icpr_mode_x, icpr_mode_y, icpr_mode_k, icpr_total_e )
          call imode_energy ( omg(:,:,-nz), iomg_mode_x, iomg_mode_y, iomg_mode_k, iomg_total_e )
        end if

        if ( rank == 0 ) then
          write( unit=omph, fmt="(f15.8, SP, 2050ES24.15e3)" ) &
               time, phi_total_e, phi_mode_x, phi_mode_y, phi_mode_k
          write( unit=omps, fmt="(f15.8, SP, 2050ES24.15e3)" ) &
               time, psi_total_e, psi_mode_x, psi_mode_y, psi_mode_k

          write( unit=oidn, fmt="(f15.8, SP, 2050ES24.15e3)" ) &
               time, idns_total_e, idns_mode_x, idns_mode_y, idns_mode_k
          write( unit=oicr, fmt="(f15.8, SP, 2050ES24.15e3)" ) &
               time, icpr_total_e, icpr_mode_x, icpr_mode_y, icpr_mode_k
          write( unit=oiog, fmt="(f15.8, SP, 2050ES24.15e3)" ) &
               time, iomg_total_e, iomg_mode_x, iomg_mode_y, iomg_mode_k
        end if

      end if


  END SUBROUTINE wrt


!--------------------------------------
  SUBROUTINE dt_check ( time, phi_total_e, psi_total_e )
!--------------------------------------

    real(kind=DP), intent(in) :: time
    real(kind=DP), intent(in) :: phi_total_e, psi_total_e

    real(kind=DP) :: cfl, dx

      dx = 2._DP * pi / kx(nx)

      cfl = sqrt( max( phi_total_e, psi_total_e ) / ( 2._DP * lz ) ) * dt / dx

      if( cfl > 0.25_DP ) then
!!!      if( cfl > 0.1_DP ) then
        write(olog,fmt="(a,1p,e15.7,a,e15.7)") "### now, cfl = ", cfl, " at t = ", time
        write(olog,fmt="(a,1p,e15.7)")         "### then, dt is changed to dt = ", dt*0.5_DP
        dt   = dt * 0.5_DP
      end if


  END SUBROUTINE dt_check
 

!--------------------------------------
  SUBROUTINE dt_check2 ( time, phi, psi )
!--------------------------------------

    real(kind=DP), intent(in) :: time
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, psi
!!! real(kind=DP), dimension(-nz:nz-1)        :: dpara, valf, b0
    real(kind=DP), dimension(0:nny-1,0:nnx-1) :: vx, vy

    real(kind=DP) :: cfl, cfl_z, cfl_v, cfl_g, dx

    integer :: ix, iy, iz


      dx = 2._DP * pi / kx(nx)

      cfl   = 0._DP
      cfl_z = 0._DP
      cfl_v = 0._DP

      do iz = -nz, nz-1
        cfl_z = max( cfl_z, abs(valf(iz)*dt/dpara(iz)) )
      end do

      do iz = -nz, nz-1
        call pssn_derivative( phi(:,:,iz), vy, vx )
        do ix = 0, nnx-1
          do iy = 0, nny-1
            cfl_v = max( cfl_v, abs(vx(ix,iy)), abs(vx(ix,iy)) )
          end do
        end do
      end do
        cfl = max( cfl_z, cfl_v*dt/dx )

        call MPI_Allreduce( cfl, cfl_g, 1, MPI_DOUBLE_PRECISION, &
                            MPI_MAX, MPI_COMM_WORLD, ierr_mpi )


      if( cfl_g > 0.7_DP ) then
!!!      if( cfl > 0.1_DP ) then
        write(olog,fmt="(a,1p,2e15.7,a,e15.7)") &
             "### now, cfl, cfl_g = ", cfl, cfl_g, " at t = ", time
        write(olog,fmt="(a,1p,e15.7)")         "### then, dt is changed to dt = ", dt*0.5_DP
        dt   = dt * 0.5_DP
      end if


  END SUBROUTINE dt_check2


!--------------------------------------
  SUBROUTINE mode_energy ( wrk, mode_x, mode_y, mode_k, total_e )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: wrk

    real(kind=DP), dimension(0:nx) :: mode_x
    real(kind=DP), dimension(0:ny) :: mode_y
    real(kind=DP), dimension(0:nk) :: mode_k
 
    real(kind=DP) :: total_e

! --- local variables

    real(kind=DP), dimension(:,:,:), allocatable :: wr3
    real(kind=DP), dimension(:,:),   allocatable :: wr2

    integer  ::  mx, my, iz

    integer, dimension(0:nk) :: ic

    real(kind=DP) :: dk
    integer       :: ikp


      allocate( wr3(-nx:nx,0:ny,-nz:nz-1) )
      allocate( wr2(-nx:nx,0:ny) )

      mode_x(:)  = 0._DP
      mode_y(:)  = 0._DP
      total_e    = 0._DP

      mode_k(:)  = 0._DP
      ic(:)      = 0


        do iz = -nz, nz-1
          do my = 0, ny
            do mx = -nx, nx
              wr3(mx,my,iz) = real( wrk(mx,my,iz) * conjg( wrk(mx,my,iz) )  &
                                   , kind=DP ) * ksq(mx,my,iz)
            end do
          end do
        end do


      call intgrl_zeta ( wr3, wr2 )


! --- ky-modes
        do my = 1, ny
          do mx = 1, nx
            mode_y(my) = mode_y(my) + wr2(mx,my) + wr2(-mx,my)
          end do
        end do
    
        mx = 0
          do my = 1, ny
            mode_y(my) = mode_y(my) + wr2(mx,my) 
          end do
    
        my = 0
          do mx = 1, nx
            mode_y(my) = mode_y(my) + wr2(mx,my) 
          end do
 
! --- kx-modes
      do mx = 1, nx
        do my = 1, ny
          mode_x(mx) = mode_x(mx) + wr2(mx,my) + wr2(-mx,my)
        end do
          mode_x(mx) = mode_x(mx) + wr2(mx,0)
      end do
  
      mx = 0
        do my = 1, ny
          mode_x(mx) = mode_x(mx) + wr2(mx,my)
        end do

  
! --- total in whole k-space
      do my = 0, ny
        total_e = total_e + mode_y(my)
      end do


! --- k_perp-modes
      dk = kx(1)
      do my = 1, ny
        do mx = 1, nx
          ikp = nint( sqrt( kx(mx)**2 + ky(my)**2 ) / dk )
          if( ikp <= nk ) then
            mode_k(ikp) = mode_k(ikp) + wr2(mx,my) + wr2(-mx,my)
            ic(ikp)  = ic (ikp) + 2
          end if
        end do
      end do

      mx = 0
        do my = 1, ny
          ikp = nint( sqrt( kx(mx)**2 + ky(my)**2 ) / dk )
          if( ikp <= nx ) then
            mode_k(ikp) = mode_k(ikp) + wr2(mx,my)
            ic(ikp)  = ic (ikp) + 1
          end if
        end do

      my = 0
        do mx = 1, nx
          ikp = nint( sqrt( kx(mx)**2 + ky(my)**2 ) / dk )
          if( ikp <= nx ) then
            mode_k(ikp) = mode_k(ikp) + wr2(mx,my)
            ic(ikp)  = ic (ikp) + 1
          end if
        end do


!        do ikp = 1, nk
!            mode_k(ikp) = mode_k(ikp) / real( ic(ikp), kind=DP )
!        end do




      deallocate( wr3 )
      deallocate( wr2 )


  END SUBROUTINE mode_energy


!--------------------------------------
  SUBROUTINE imode_energy ( wrk, mode_x, mode_y, mode_k, total_e )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny) :: wrk

    real(kind=DP), dimension(0:nx) :: mode_x
    real(kind=DP), dimension(0:ny) :: mode_y
    real(kind=DP), dimension(0:nk) :: mode_k
 
    real(kind=DP) :: total_e

! --- local variables

    real(kind=DP), dimension(:,:),   allocatable :: wr2

    integer  ::  mx, my, iz

    integer, dimension(0:nk) :: ic

    real(kind=DP) :: dk
    integer       :: ikp


      allocate( wr2(-nx:nx,0:ny) )

!! debug
!            if ( rankz == 0 ) then
!              write(olog,fmt="(a)") "# zz, wrk in imode"
!              mx   =-2
!              my   = 1
!              write(olog,fmt="(1p,5e15.7)") zz(-nz), wrk(mx,my)
!              write(olog,*)
!            end if
!! debug

      mode_x(:)  = 0._DP
      mode_y(:)  = 0._DP
      total_e    = 0._DP

      mode_k(:)  = 0._DP
      ic(:)      = 0

          do my = 0, ny
            do mx = -nx, nx
              wr2(mx,my) = real( wrk(mx,my) * conjg( wrk(mx,my) )  &
                               , kind=DP )
            end do
          end do


! --- ky-modes
        do my = 1, ny
          do mx = 1, nx
            mode_y(my) = mode_y(my) + wr2(mx,my) + wr2(-mx,my)
          end do
        end do
    
        mx = 0
          do my = 1, ny
            mode_y(my) = mode_y(my) + wr2(mx,my) 
          end do
    
        my = 0
          do mx = 1, nx
            mode_y(my) = mode_y(my) + wr2(mx,my) 
          end do
 
! --- kx-modes
      do mx = 1, nx
        do my = 1, ny
          mode_x(mx) = mode_x(mx) + wr2(mx,my) + wr2(-mx,my)
        end do
          mode_x(mx) = mode_x(mx) + wr2(mx,0)
      end do
  
      mx = 0
        do my = 1, ny
          mode_x(mx) = mode_x(mx) + wr2(mx,my)
        end do
  

! --- total in whole k-space
      do my = 0, ny
        total_e = total_e + mode_y(my)
      end do


! --- k_perp-modes
      dk = kx(1)
      do my = 1, ny
        do mx = 1, nx
          ikp = nint( sqrt( kx(mx)**2 + ky(my)**2 ) / dk )
          if( ikp <= nk ) then
            mode_k(ikp) = mode_k(ikp) + wr2(mx,my) + wr2(-mx,my)
            ic(ikp)  = ic (ikp) + 2
          end if
        end do
      end do

      mx = 0
        do my = 1, ny
          ikp = nint( sqrt( kx(mx)**2 + ky(my)**2 ) / dk )
          if( ikp <= nx ) then
            mode_k(ikp) = mode_k(ikp) + wr2(mx,my)
            ic(ikp)  = ic (ikp) + 1
          end if
        end do

      my = 0
        do mx = 1, nx
          ikp = nint( sqrt( kx(mx)**2 + ky(my)**2 ) / dk )
          if( ikp <= nx ) then
            mode_k(ikp) = mode_k(ikp) + wr2(mx,my)
            ic(ikp)  = ic (ikp) + 1
          end if
        end do


!        do ikp = 1, nk
!            mode_k(ikp) = mode_k(ikp) / real( ic(ikp), kind=DP )
!        end do



      deallocate( wr2 )


  END SUBROUTINE imode_energy


END MODULE RMHDS_diag
