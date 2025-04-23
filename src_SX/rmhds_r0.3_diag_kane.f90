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
  use RMHDS_pssn,   only: pssn_brackets
  use RMHDS_advnc,  only: dpdz, dudz_dddz, calxi

  implicit none

  private

  integer, parameter :: nk = min( nx, ny )

  public   diag_cntrl, dt_check


CONTAINS


!--------------------------------------
  SUBROUTINE diag_cntrl( omg, psi, dns, time, id )
!--------------------------------------

    complex(kind=DP), intent(in), &
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


      if( id == 0 ) then

        tout_mag  = ( int( ( time + eps )/dtout_mag ) + 1 ) * dtout_mag
        tout_ion  = ( int( ( time + eps )/dtout_ion ) + 1 ) * dtout_ion
        tout_eng  = ( int( ( time + eps )/dtout_eng ) + 1 ) * dtout_eng

       if ( dtout_mag /= 0 ) then 
!          call wrt ( omg, psi, dns, phi, cpr, time, 0 )
           write(olog,fmt="(a)")"# dtout_mag exist!!"
       end if

        if ( dtout_ion /= 0 ) then 
!          call wrt ( omg, psi, dns, phi, cpr, time, 1 )
           write(olog,fmt="(a)")"# dtout_ion exist!!"
        end if  

       if ( time == 0 ) then 
          call wrt ( omg, psi, dns, phi, cpr, time, 0 )
          write( olog, * ) &
            " # Magnetospheric data output at time = ", time
       end if


! debug
!            mx   = nx*3/4
!            my   = ny*3/4
            mx   =-2
            my   = 1
            write(olog,fmt="(a,1p,e15.7)") "# time = ", time
            write(olog,fmt="(a,2i5)")      "# mx, my = ", mx, my
!            write(olog,fmt="(a,2i5)")      "# mx, my = ", mxi, myi
            write(olog,fmt="(a)")          "# zz, omg, psi "
            do iz = -nz, nz-1
!              write(olog,fmt="(1p,5e15.7)")  zz(iz), omg(mxi,myi,iz), psi(mxi,myi,iz)
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


        if ( time== 0 ) then 
          call wrt ( omg, psi, dns, phi, cpr, time, 1 )
            write( olog, * ) &
               " # Ionospheric data output at time = ", time
        end if

    
        call wrt ( omg, psi, dns, phi, cpr, time, 2 )


!      if( id == 1 ) then
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
!            write(olog,fmt="(a,2i5)")      "# mx, my = ", mxi, myi
            write(olog,fmt="(a)")          "# zz, omg, psi "
            do iz = -nz, nz-1
!              write(olog,fmt="(1p,5e15.7)")  zz(iz), omg(mxi,myi,iz), psi(mxi,myi,iz)
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


!    integer ::mx, my, iz


      rewind ocnt
      write( unit=ocnt ) time, dt, omg, psi, dns


  END SUBROUTINE contnu



!--------------------------------------
  SUBROUTINE wrt ( omg, psi, dns, phi, cpr, time, id )
!--------------------------------------

! --- arguments
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg, psi, dns

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr

    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1)  :: Epr

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

    real(kind=DP) :: vis_total_e, EJp_total_e
 
    real(kind=DP), dimension(0:nx) :: vis_mode_x, EJp_mode_x
    real(kind=DP), dimension(0:ny) :: vis_mode_y, EJp_mode_y
    real(kind=DP), dimension(0:nk) :: vis_mode_k, EJp_mode_k

    real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: ptg_mode_x, ptg_mode_y, ptg_total, ptg_mode_k
    real(kind=DP), dimension(-nx:nx,0:ny) :: ptg_x, ptg_y, ptg_k, ptg_t

    real(kind=DP), save :: EJp_t_all, vis_t_all, ptg_t_all
 
!    integer        ::  mx, my, iz
      if( time == 0) then
         EJp_t_all = 0.d0
         vis_t_all = 0.d0
         ptg_t_all = 0.d0
      end if

      if( id == 0 ) then

          call calc_Epara( psi, phi, Epr )
          write( unit=omag ) time, omg, cpr, Epr

      else if( id == 1 ) then

        if ( rankz == 0 ) then
          write( unit=oion ) time, dns(:,:,-nz), omg(:,:,-nz), cpr(:,:,-nz)
        end if

      else if( id == 2 ) then

        call calc_Epara( psi, phi, Epr )
        call mode_energy ( phi, phi_mode_x, phi_mode_y, phi_mode_k, phi_total_e, 0 )
        call mode_energy ( psi, psi_mode_x, psi_mode_y, psi_mode_k, psi_total_e, 1 )
        call mode_energy ( omg, vis_mode_x, vis_mode_y, vis_mode_k, vis_total_e, 2 )
        call mode_EJpara ( Epr, cpr, EJp_mode_x, EJp_mode_y, EJp_mode_k, EJp_total_e )

        call dt_check ( time, phi, psi )


!! debug
!            if ( rankz == 0 ) then
!              write(olog,fmt="(a)") "# zz, dns in wrt"
!              mx   =-2
!              my   = 1
!              write(olog,fmt="(1p,5e15.7)") zz(-nz), dns(mx,my,-nz)
!              write(olog,*)
!            end if
!! debug

        ptg_mode_x(:,:,:) = 0.d0
        ptg_mode_y(:,:,:) = 0.d0
        ptg_mode_k(:,:,:) = 0.d0
        ptg_total (:,:,:) = 0.d0

        if ( rankz == nprocz-1 .and. nk == nx ) then
          call pointing_flux ( phi(:,:,nz-1), psi(:,:,nz-1), ptg_mode_x(0:nx,0,nz-1), &
                          ptg_mode_y(0,0:ny,nz-1), ptg_mode_k(0:nk,0,nz-1), ptg_total(0,0,nz-1) )
        else if ( rankz == nprocz-1 .and. nk == ny ) then
          call pointing_flux ( phi(:,:,nz-1), psi(:,:,nz-1), ptg_mode_x(0:nx,0,nz-1), &
                          ptg_mode_y(0,0:ny,nz-1), ptg_mode_k(0,0:nk,nz-1), ptg_total(0,0,nz-1) )
        end if

        if ( rankz == 0 .and. nk == nx ) then
          call imode_energy ( dns(:,:,-nz), idns_mode_x, idns_mode_y, idns_mode_k, idns_total_e )
          call imode_energy ( cpr(:,:,-nz), icpr_mode_x, icpr_mode_y, icpr_mode_k, icpr_total_e )
          call imode_energy ( omg(:,:,-nz), iomg_mode_x, iomg_mode_y, iomg_mode_k, iomg_total_e )
          call pointing_flux ( phi(:,:,-nz), psi(:,:,-nz), -ptg_mode_x(0:nx,0,-nz), &
                          -ptg_mode_y(0,0:ny,-nz), -ptg_mode_k(0:nx,0,-nz), -ptg_total(0,0,-nz) )
        else if ( rankz == 0 .and. nk == ny ) then
          call imode_energy ( dns(:,:,-nz), idns_mode_x, idns_mode_y, idns_mode_k, idns_total_e )
          call imode_energy ( cpr(:,:,-nz), icpr_mode_x, icpr_mode_y, icpr_mode_k, icpr_total_e )
          call imode_energy ( omg(:,:,-nz), iomg_mode_x, iomg_mode_y, iomg_mode_k, iomg_total_e )
          call pointing_flux ( phi(:,:,-nz), psi(:,:,-nz), -ptg_mode_x(0:nx,0,-nz), &
                          -ptg_mode_y(0,0:ny,-nz), -ptg_mode_k(0,0:nk,-nz), -ptg_total(0,0,-nz) )
        end if
        
         call intgrl_zeta ( ptg_mode_x, ptg_x )
         call intgrl_zeta ( ptg_mode_y, ptg_y )
         call intgrl_zeta ( ptg_mode_k, ptg_k )
         call intgrl_zeta ( ptg_total, ptg_t )

        if ( rank == 0 ) then
          write( unit=omph, fmt="(f15.8, SP, 1026ES24.15e3)" ) &
               time, phi_total_e, phi_mode_x, phi_mode_y, phi_mode_k
          write( unit=omps, fmt="(f15.8, SP, 1026ES24.15e3)" ) &
               time, psi_total_e, psi_mode_x, psi_mode_y, psi_mode_k

          write( unit=oidn, fmt="(f15.8, SP, 1026ES24.15e3)" ) &
               time, idns_total_e, idns_mode_x, idns_mode_y, idns_mode_k
          write( unit=oicr, fmt="(f15.8, SP, 1026ES24.15e3)" ) &
               time, icpr_total_e, icpr_mode_x, icpr_mode_y, icpr_mode_k
          write( unit=oiog, fmt="(f15.8, SP, 1026ES24.15e3)" ) &
               time, iomg_total_e, iomg_mode_x, iomg_mode_y, iomg_mode_k
          ptg_t_all = ptg_t_all + ptg_t(0,0)
          if ( rankz == 0 .and. nk == nx ) then
           write( unit=optg, fmt="(f15.8, SP, 1027ES24.15e3)" ) &
               time, -ptg_t_all, ptg_t(0,0), ptg_x(0:nx,0), ptg_y(0,:), ptg_k(0:nk,0)
          else if ( rankz == 0 .and. nk == ny ) then
           write( unit=optg, fmt="(f15.8, SP, 1027ES24.15e3)" ) &
               time, -ptg_t_all, ptg_t(0,0), ptg_x(0:nx,0), ptg_y(0,:), ptg_k(0,:)
          end if

          vis_t_all = vis_t_all + vis_total_e
          write( unit=ovis, fmt="(f15.8, SP, 1027ES24.15e3)" ) &
               time, vis_t_all, vis_total_e, vis_mode_x, vis_mode_y, vis_mode_k
          EJp_t_all = EJp_t_all + EJp_total_e
          write( unit=oejp, fmt="(f15.8, SP, 1027ES24.15e3)" ) &
               time, EJp_t_all, EJp_total_e, EJp_mode_x, EJp_mode_y, EJp_mode_k
        end if

      end if

  END SUBROUTINE wrt


!--------------------------------------
  SUBROUTINE dt_check ( time, phi, psi )
!--------------------------------------

    real(kind=DP), intent(in) :: time
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, psi

    real(kind=DP), dimension(0:nx) :: phi_mode_x, psi_mode_x
    real(kind=DP), dimension(0:ny) :: phi_mode_y, psi_mode_y
    real(kind=DP), dimension(0:nk) :: phi_mode_k, psi_mode_k

    real(kind=DP) :: phi_total_e, psi_total_e

    real(kind=DP) :: cfl, dx

      dx = 2._DP * pi / kx(nx)
        call mode_energy ( phi, phi_mode_x, phi_mode_y, phi_mode_k, phi_total_e, 0 )
        call mode_energy ( psi, psi_mode_x, psi_mode_y, psi_mode_k, psi_total_e, 0 )

      cfl = sqrt( max( phi_total_e, psi_total_e )/ ( 2._DP * lz ) ) * dt / dx

      if( cfl > 0.25_DP ) then
!!!      if( cfl > 0.1_DP ) then
        write(olog,fmt="(a,1p,e15.7,a,e15.7)") "### now, cfl = ", cfl, " at t = ", time
        write(olog,fmt="(a,1p,e15.7)")         "### then, dt is changed to dt = ", dt*0.5_DP
        dt   = dt * 0.5_DP
      end if


  END SUBROUTINE dt_check
 

!--------------------------------------
  SUBROUTINE mode_energy ( wrk, mode_x, mode_y, mode_k, total_e, id )
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
    integer, intent(in) :: id

      allocate( wr3(-nx:nx,0:ny,-nz:nz-1) )
      allocate( wr2(-nx:nx,0:ny) )

      mode_x(:)  = 0._DP
      mode_y(:)  = 0._DP
      total_e    = 0._DP

      mode_k(:)  = 0._DP
      ic(:)      = 0

      if ( id == 0 ) then

        do iz = -nz, nz-1
          do my = 0, ny
            do mx = -nx, nx
              wr3(mx,my,iz) = real( wrk(mx,my,iz) * conjg( wrk(mx,my,iz) )  &
                                   , kind=DP ) * ksq(mx,my,iz) * 0.5d0 / valf(iz)**2
            end do
          end do
        end do

      else if ( id  == 1 ) then

        do iz = -nz, nz-1
          do my = 0, ny
            do mx = -nx, nx
              wr3(mx,my,iz) = real( wrk(mx,my,iz) * conjg( wrk(mx,my,iz) )  &
                                   , kind=DP ) * ksq(mx,my,iz) * 0.5d0
            end do
          end do
        end do

      else
       
        do iz = -nz, nz-1
          do my = 0, ny
            do mx = -nx, nx
              wr3(mx,my,iz) =  real( wrk(mx,my,iz) * conjg( wrk(mx,my,iz) )  &
                                   , kind=DP ) * nu / valf(iz)**2
            end do
          end do
        end do

      end if
       
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
  SUBROUTINE mode_EJpara ( Epr, cpr, mode_x, mode_y, mode_k, total_e )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: Epr, cpr

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
              wr3(mx,my,iz) = real( Epr(mx,my,iz) * conjg( cpr(mx,my,iz) )  &
                                   , kind=DP ) 
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


  END SUBROUTINE mode_EJpara

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

    integer  ::  mx, my!, iz

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

!--------------------------------------
  SUBROUTINE pointing_flux ( phi, psi, mode_x, mode_y, mode_k, total_e )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny) :: phi, psi

    real(kind=DP), dimension(0:nx) :: mode_x
    real(kind=DP), dimension(0:ny) :: mode_y
    real(kind=DP), dimension(0:nk) :: mode_k
 
    real(kind=DP) :: total_e

! --- local variables

    real(kind=DP), dimension(:,:),   allocatable :: wr2

    integer  ::  mx, my

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
              wr2(mx,my) = real( ui*kx(mx)*phi(mx,my) * conjg( ui*kx(mx)*psi(mx,my)), kind=DP ) &
                         + real( ui*ky(my)*phi(mx,my) * conjg( ui*ky(my)*psi(mx,my)), kind=DP )
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


  END SUBROUTINE pointing_flux

!--------------------------------------
  SUBROUTINE calc_Epara( psi, phi, Epara )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, psi
    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1)  :: xi

    complex(kind=DP), dimension(:,:,:), allocatable :: psnb, psnb1
    complex(kind=DP), dimension(:,:,:), allocatable :: dphidz
    complex(kind=DP), dimension(:,:,:), allocatable :: upg, dwg
    complex(kind=DP), dimension(:,:,:), allocatable :: dudz, dddz

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: Epara

    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: dpsdt
    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: aphi, apsi

    integer  ::  mx, my, iz

      allocate( psnb(-nx:nx,0:ny,-nz:nz-1) )
      allocate( psnb1(-nx:nx,0:ny,-nz:nz-1) )

      allocate( dphidz(-nx:nx,0:ny,-nz:nz-1) )

      allocate( upg(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dwg(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dudz(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dddz(-nx:nx,0:ny,-nz:nz-1) )

      if ( trim(fd_type) == "central" ) then
       aphi(:,:,:) = phi(:,:,:)
       call dpdz( aphi, dphidz,  1 )
      else if ( trim(fd_type) == "upwind" ) then
        upg(:,:,:) = sqrt(ksqi_xi(:,:,:))*phi(:,:,:) - psi(:,:,:)
        dwg(:,:,:) = sqrt(ksqi_xi(:,:,:))*phi(:,:,:) + psi(:,:,:)
       call dudz_dddz( upg, dwg, dudz, dddz )
        dphidz(:,:,:) = 0.5_DP * ( dudz(:,:,:) + dddz(:,:,:) ) * sqrt(ksq_xi(:,:,:))
      end if

      apsi(:,:,:) = psi(:,:,:)
      call calxi(apsi, xi, 1)

      if( trim(calc_type) == "nonlinear" ) then
        call pssn_brackets( xi,  phi, psnb, 2*nz )
        call pssn_brackets( psi, phi, psnb1, 2*nz )
        if( rankz == 0 ) then
          psnb(:,:,-nz) = ( 0._DP, 0._DP )
!          psnb1(:,:,-nz) = ( 0._DP, 0._DP )
        end if
      else
        psnb(:,:,:) = ( 0._DP, 0._DP )
        psnb1(:,:,:) = ( 0._DP, 0._DP )
      end if

      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            dpsdt(mx,my,iz) =   ksqi_xi(mx,my,iz) * (      &
                   ( psnb(mx,my,iz) + dphidz(mx,my,iz) )   &
                   - eta * ksq(mx,my,iz) * psi(mx,my,iz)   &
                   - psigma*ksq(mx,my,iz)**2*psi(mx,my,iz) )
          end do
        end do
      end do

      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
           Epara(mx,my,iz) = - ( dphidz(mx,my,iz) + psnb1(mx,my,iz) ) &
                             + dpsdt(mx,my,iz)
          end do
        end do
      end do

      deallocate( psnb )
      deallocate( psnb1 )
      deallocate( dphidz )

      deallocate( upg, dwg )
      deallocate( dudz, dddz )

  END SUBROUTINE calc_Epara

END MODULE RMHDS_diag
