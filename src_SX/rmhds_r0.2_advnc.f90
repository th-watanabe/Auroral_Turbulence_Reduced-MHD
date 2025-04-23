MODULE RMHDS_advnc
!-------------------------------------------------------------------------------
!
!    Time integration of the reduced MHD equations
!
!-------------------------------------------------------------------------------

  use RMHDS_header
  use RMHDS_mpienv
  use RMHDS_bndry, only: bndry_bound
  use RMHDS_fld,   only: fld_calfld
  use RMHDS_pssn,  only: pssn_brackets

  implicit none

  private

  public   advnc_rkgsteps


CONTAINS


!--------------------------------------
  SUBROUTINE advnc_rkgsteps( omg, psi, dns )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg, psi, dns


    complex(kind=DP), dimension(:,:,:), allocatable :: domg, qomg
    complex(kind=DP), dimension(:,:,:), allocatable :: dpsi, qpsi
    complex(kind=DP), dimension(:,:,:), allocatable :: ddns, qdns

    integer :: istep, mx, my, iz

      allocate( domg(-nx:nx,0:ny,-nz:nz-1), qomg(-nx:nx,0:ny,-nz:nz-1) ) 
      allocate( dpsi(-nx:nx,0:ny,-nz:nz-1), qpsi(-nx:nx,0:ny,-nz:nz-1) ) 
      allocate( ddns(-nx:nx,0:ny,-nz:nz-1), qdns(-nx:nx,0:ny,-nz:nz-1) ) 


      domg(:,:,:) = ( 0._DP, 0._DP )
      qomg(:,:,:) = ( 0._DP, 0._DP )
      dpsi(:,:,:) = ( 0._DP, 0._DP )
      qpsi(:,:,:) = ( 0._DP, 0._DP )
      ddns(:,:,:) = ( 0._DP, 0._DP )
      qdns(:,:,:) = ( 0._DP, 0._DP )


      do istep = 1, 4

!          print *, "# caldlt begin:  rankz, istep = ", rankz, istep

        if     ( trim(fd_type) == "central" ) then
          call caldlt_c( omg, psi, dns, domg, dpsi, ddns )
        else if( trim(fd_type) == "upwind"  ) then
          call caldlt_u( omg, psi, dns, domg, dpsi, ddns )
        else
          print *, "# selection of fd_type is invalid"
          stop
        end if

!          print *, "# caldlt end:  rankz, istep = ", rankz, istep

        call rkg   ( omg, psi, dns, domg, dpsi, ddns,       &
                                    qomg, qpsi, qdns, istep )

!          print *, "# rkg end:  istep = ", istep

      end do 

      deallocate( domg, qomg ) 
      deallocate( dpsi, qpsi ) 
      deallocate( ddns, qdns ) 


  END SUBROUTINE advnc_rkgsteps


!--------------------------------------
  SUBROUTINE rkg( omg, psi, dns, domg, dpsi, ddns, qomg, qpsi, qdns, istep )
!--------------------------------------
!     Runge-Kutta-Gill

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg,  psi,  dns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: domg, dpsi, ddns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: qomg, qpsi, qdns

    integer, intent(in) :: istep

    real(kind=DP) :: c1, c2, cq, c0
    integer :: mx, my, iz


      if      ( istep == 1 ) then
        c1   =  0.5_DP
        c2   = -1._DP
        cq   = -2._DP
        c0   =  1._DP
      else if ( istep == 2 ) then
        c1   =  1._DP - sqrt( 0.5_DP )
        c2   = -c1
        cq   =  1._DP - 3._DP * c1
        c0   =  2._DP * c1
      else if ( istep == 3 ) then
        c1   =  1._DP + sqrt( 0.5_DP )
        c2   = -c1
        cq   =  1._DP - 3._DP * c1
        c0   =  2._DP * c1
      else if ( istep == 4 ) then
        c1   =  1._DP / 6._DP
        c2   = -1._DP / 3._DP
        cq   =  0._DP
        c0   =  0._DP
      end if


      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            omg (mx,my,iz) = omg(mx,my,iz) &
                           + c1 * domg(mx,my,iz) + c2 * qomg(mx,my,iz)
            psi (mx,my,iz) = psi(mx,my,iz) &
                           + c1 * dpsi(mx,my,iz) + c2 * qpsi(mx,my,iz)
            dns (mx,my,iz) = dns(mx,my,iz) &
                           + c1 * ddns(mx,my,iz) + c2 * qdns(mx,my,iz)

            qomg(mx,my,iz) = cq * qomg(mx,my,iz) + c0 * domg(mx,my,iz)
            qpsi(mx,my,iz) = cq * qpsi(mx,my,iz) + c0 * dpsi(mx,my,iz)
            qdns(mx,my,iz) = cq * qdns(mx,my,iz) + c0 * ddns(mx,my,iz)
          end do
        end do
      end do


  END SUBROUTINE rkg


!--------------------------------------
  SUBROUTINE caldlt_c( omg, psi, dns, domg, dpsi, ddns )
!--------------------------------------
!     increment within a time step

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg,  psi,  dns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: domg, dpsi, ddns

    complex(kind=DP), dimension(:,:,:), allocatable :: phi, cpr
    complex(kind=DP), dimension(:,:,:), allocatable :: psnb1, psnb2, psnb3
    complex(kind=DP), dimension(:,:,:), allocatable :: dphidz, dcprdz

    integer  ::  mx, my, iz


      allocate( phi(-nx:nx,0:ny,-nz:nz-1) )
      allocate( cpr(-nx:nx,0:ny,-nz:nz-1) )

      allocate( psnb1(-nx:nx,0:ny,-nz:nz-1) )
      allocate( psnb2(-nx:nx,0:ny,-nz:nz-1) )
      allocate( psnb3(-nx:nx,0:ny,-nz:nz-1) )

      allocate( dphidz(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dcprdz(-nx:nx,0:ny,-nz:nz-1) )

!          print *, "# calfld begin: rankz = ", rankz

      call fld_calfld( omg, psi, dns, phi, cpr )

!          print *, "# calfld end: rankz = ", rankz
     
      if( trim(calc_type) == "nonlinear" ) then

        call pssn_brackets( phi, omg, psnb1, 2*nz )
        call pssn_brackets( psi, cpr, psnb2, 2*nz )
        call pssn_brackets( psi, phi, psnb3, 2*nz )

! debug
!        psnb1(:,:,:) = ( 0._DP, 0._DP )
!        psnb2(:,:,:) = ( 0._DP, 0._DP )
!        psnb3(:,:,:) = ( 0._DP, 0._DP )
! debug

! debug
         if( rankz == 0 ) then
           psnb3(:,:,-nz) = ( 0._DP, 0._DP )
         end if
! debug

      else

        psnb1(:,:,:) = ( 0._DP, 0._DP )
        psnb2(:,:,:) = ( 0._DP, 0._DP )
        psnb3(:,:,:) = ( 0._DP, 0._DP )

      end if

!          print *, "# dpdz begin: rankz = ", rankz

      call dpdz( phi, dphidz,  1 )
      call dpdz( cpr, dcprdz, -1 )

!          print *, "# dpdz end: rankz = ", rankz


      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            domg(mx,my,iz) = (     -psnb1(mx,my,iz)                       &
                   + valf(iz)**2 * ( psnb2(mx,my,iz) + dcprdz(mx,my,iz) ) &
                   - nu  * ksq(mx,my,iz) * omg(mx,my,iz) ) * dt
            dpsi(mx,my,iz) = (                                            &
                   +               ( psnb3(mx,my,iz) + dphidz(mx,my,iz) ) &
                   - eta * ksq(mx,my,iz) * psi(mx,my,iz) ) * dt
            ddns(mx,my,iz) = ( 0._DP, 0._DP )
          end do
        end do
      end do


!  zero-zero
      do iz = -nz, nz-1
        omg(0,0,iz)   = ( 0._DP, 0._DP )
        psi(0,0,iz)   = ( 0._DP, 0._DP )
        dns(0,0,iz)   = ( 0._DP, 0._DP )
      end do


! ionospheric boundary
      if( rankz == 0 ) then
        call iono_dns( dns(:,:,-nz), cpr(:,:,-nz), phi(:,:,-nz), ddns(:,:,-nz) )
        domg(:,:,-nz) = ( 0._DP, 0._DP )
      end if


! --- reality condition
      my = 0
        do iz = -nz, nz-1
          do mx = 1, nx
            domg(mx,my,iz) = conjg( domg(-mx,-my,iz) )
            dpsi(mx,my,iz) = conjg( dpsi(-mx,-my,iz) )
            ddns(mx,my,iz) = conjg( ddns(-mx,-my,iz) )
          end do
        end do
         

      if ( trim(fd_filt) == "on" ) then
        call zfilter ( domg,  1 )
!!!        call zfilter ( dpsi, -1 )
      end if


      deallocate( phi, cpr )
      deallocate( psnb1, psnb2, psnb3 )
      deallocate( dphidz, dcprdz )


  END SUBROUTINE caldlt_c


!--------------------------------------
  SUBROUTINE caldlt_u( omg, psi, dns, domg, dpsi, ddns )
!--------------------------------------
!     increment within a time step
!     with up-wind for the parallel derivatives
!     Note: this is valid only for the uniform V_A case

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg,  psi,  dns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: domg, dpsi, ddns

    complex(kind=DP), dimension(:,:,:), allocatable :: phi, cpr
    complex(kind=DP), dimension(:,:,:), allocatable :: psnb1, psnb2, psnb3
    complex(kind=DP), dimension(:,:,:), allocatable :: dphidz, dcprdz

    complex(kind=DP), dimension(:,:,:), allocatable :: upg, dwg
    complex(kind=DP), dimension(:,:,:), allocatable :: dudz, dddz

    integer  ::  mx, my, iz


      allocate( phi(-nx:nx,0:ny,-nz:nz-1) )
      allocate( cpr(-nx:nx,0:ny,-nz:nz-1) )

      allocate( psnb1(-nx:nx,0:ny,-nz:nz-1) )
      allocate( psnb2(-nx:nx,0:ny,-nz:nz-1) )
      allocate( psnb3(-nx:nx,0:ny,-nz:nz-1) )

      allocate( dphidz(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dcprdz(-nx:nx,0:ny,-nz:nz-1) )

      allocate( upg(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dwg(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dudz(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dddz(-nx:nx,0:ny,-nz:nz-1) )


      call fld_calfld( omg, psi, dns, phi, cpr )

     
      if( trim(calc_type) == "nonlinear" ) then

        call pssn_brackets( phi, omg, psnb1, 2*nz )
        call pssn_brackets( psi, cpr, psnb2, 2*nz )
        call pssn_brackets( psi, phi, psnb3, 2*nz )

        if( rankz == 0 ) then
          psnb3(:,:,-nz) = ( 0._DP, 0._DP )
        end if

      else

        psnb1(:,:,:) = ( 0._DP, 0._DP )
        psnb2(:,:,:) = ( 0._DP, 0._DP )
        psnb3(:,:,:) = ( 0._DP, 0._DP )

      end if


!--- up-wind
        upg(:,:,:)   = phi(:,:,:) - psi(:,:,:)
        dwg(:,:,:)   = phi(:,:,:) + psi(:,:,:)

      call dudz_dddz( upg, dwg, dudz, dddz )

        dphidz(:,:,:) = 0.5_DP * ( dudz(:,:,:) + dddz(:,:,:) )
        dcprdz(:,:,:) = 0.5_DP * ( dudz(:,:,:) - dddz(:,:,:) ) * ksq(:,:,:)
!--- up-wind


      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            domg(mx,my,iz) = (     -psnb1(mx,my,iz)                       &
                   + valf(iz)**2 * ( psnb2(mx,my,iz) + dcprdz(mx,my,iz) ) &
                   - nu  * ksq(mx,my,iz) * omg(mx,my,iz) ) * dt
            dpsi(mx,my,iz) = (                                            &
                   +               ( psnb3(mx,my,iz) + dphidz(mx,my,iz) ) &
                   - eta * ksq(mx,my,iz) * psi(mx,my,iz) ) * dt
            ddns(mx,my,iz) = ( 0._DP, 0._DP )
          end do
        end do
      end do


!  zero-zero
      do iz = -nz, nz-1
        omg(0,0,iz)   = ( 0._DP, 0._DP )
        psi(0,0,iz)   = ( 0._DP, 0._DP )
        dns(0,0,iz)   = ( 0._DP, 0._DP )
      end do


! ionospheric boundary
      if( rankz == 0 ) then
        call iono_dns( dns(:,:,-nz), cpr(:,:,-nz), phi(:,:,-nz), ddns(:,:,-nz) )
        domg(:,:,-nz) = ( 0._DP, 0._DP )
      end if


! --- reality condition
      my = 0
        do iz = -nz, nz-1
          do mx = 1, nx
            domg(mx,my,iz) = conjg( domg(-mx,-my,iz) )
            dpsi(mx,my,iz) = conjg( dpsi(-mx,-my,iz) )
            ddns(mx,my,iz) = conjg( ddns(-mx,-my,iz) )
          end do
        end do
         

!      if ( trim(fd_filt) == "on" ) then
!        call zfilter ( domg,  1 )
!        call zfilter ( dpsi, -1 )
!      end if


      deallocate( phi, cpr )
      deallocate( psnb1, psnb2, psnb3 )
      deallocate( dphidz, dcprdz )

      deallocate( upg, dwg )
      deallocate( dudz, dddz )


  END SUBROUTINE caldlt_u


!--------------------------------------
  SUBROUTINE dpdz ( w, dwdz, isw )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: w, dwdz

    integer, intent(in)  ::  isw

    integer  ::  mx, my, iz

! --- local variables

    complex(kind=DP), dimension(:,:,:), allocatable :: w_l

    real(kind=DP),dimension(-nz:nz-1) :: cefz


      do iz = -nz, nz-1
        cefz(iz)   = 1._DP / ( 12._DP * dpara(iz) )
      end do


      allocate( w_l(-nx:nx,0:ny,-nz-nb:nz-1+nb) )


      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            w_l(mx,my,iz) = w(mx,my,iz)
          end do
        end do
      end do

!          print *, "# bndry begin: rankz = ", rankz

      call bndry_bound ( w_l )

!          print *, "# bndry end: rankz = ", rankz


      if( rankz == nprocz-1 ) then

        if( isw >= 0 ) then        ! symmetric

          do my = 0, ny
            do mx = -nx, nx
              w_l(mx,my,nz  ) = ( 4._DP*w_l(mx,my,nz-1) - w_l(mx,my,nz-2) ) / 3._DP
              w_l(mx,my,nz+1) =   w_l(mx,my,nz-1)
            end do
          end do

        else if( isw < 0 ) then    ! anti-symmetric

          do my = 0, ny
            do mx = -nx, nx
              w_l(mx,my,nz  ) = ( 0._DP, 0._DP)
              w_l(mx,my,nz+1) = - w_l(mx,my,nz-1)
            end do
          end do

        end if

      end if


! --- gradient
      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            dwdz(mx,my,iz) =                          &
                 cefz(iz) * ( -     w_l(mx,my,iz+2)   &
                          + 8._DP * w_l(mx,my,iz+1)   &
                          - 8._DP * w_l(mx,my,iz-1)   &
                          +         w_l(mx,my,iz-2) )
          end do
        end do
      end do


      if( rankz == 0 ) then

            iz = -nz
            dwdz(:,:,iz) = (       - w_l(:,:,iz+2) &
                           + 4._DP * w_l(:,:,iz+1) &
                           - 3._DP * w_l(:,:,iz  ) &
                           ) / ( dpara(iz) * 2._DP )
            iz = -nz+1
            dwdz(:,:,iz) = (  w_l(:,:,iz+1) -  w_l(:,:,iz-1) &
                           ) / ( dpara(iz) * 2._DP )

      end if



  END SUBROUTINE dpdz


!--------------------------------------
  SUBROUTINE dudz_dddz( upg, dwg, dudz, dddz )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: upg, dwg

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: dudz, dddz

    integer  ::  mx, my, iz

! --- local variables

    complex(kind=DP), dimension(:,:,:), allocatable :: u_l, d_l

    real(kind=DP),dimension(-nz:nz-1) :: cefz


      do iz = -nz, nz-1
        cefz(iz)   = 1._DP / ( 6._DP * dpara(iz) )
      end do


      allocate( u_l(-nx:nx,0:ny,-nz-nb:nz-1+nb) )
      allocate( d_l(-nx:nx,0:ny,-nz-nb:nz-1+nb) )


      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            u_l(mx,my,iz) = upg(mx,my,iz)
            d_l(mx,my,iz) = dwg(mx,my,iz)
          end do
        end do
      end do


      call bndry_bound ( u_l )
      call bndry_bound ( d_l )


      if( rankz == nprocz-1 ) then

          do my = 0, ny
            do mx = -nx, nx
              u_l(mx,my,nz  ) = ( 4._DP*u_l(mx,my,nz-1) - u_l(mx,my,nz-2) ) / 3._DP
              u_l(mx,my,nz+1) =   u_l(mx,my,nz-1)

              d_l(mx,my,nz  ) =   u_l(mx,my,nz  )
              d_l(mx,my,nz+1) =   u_l(mx,my,nz-1)
            end do
          end do

      end if


! --- gradient
      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            dudz(mx,my,iz) = cefz(iz) * (   2._DP * u_l(mx,my,iz+1) &
                                          + 3._DP * u_l(mx,my,iz  ) &
                                          - 6._DP * u_l(mx,my,iz-1) &
                                          +         u_l(mx,my,iz-2) )

            dddz(mx,my,iz) = cefz(iz) * ( - 2._DP * d_l(mx,my,iz-1) &
                                          - 3._DP * d_l(mx,my,iz  ) &
                                          + 6._DP * d_l(mx,my,iz+1) &
                                          -         d_l(mx,my,iz+2) )
          end do
        end do
      end do


      if( rankz == 0 ) then

        do my = 0, ny
          do mx = -nx, nx
            iz = -nz+1
            dudz(mx,my,iz) = ( u_l(mx,my,iz) - u_l(mx,my,iz-1) ) / dpara(iz)

            iz = -nz
            dddz(mx,my,iz) = (     - d_l(mx,my,iz+2) &
                           + 4._DP * d_l(mx,my,iz+1) &
                           - 3._DP * d_l(mx,my,iz  ) &
                           ) / ( dpara(iz) * 2._DP )
            dudz(mx,my,iz) = dudz(mx,my,iz+1)
          end do
        end do

      end if


  END SUBROUTINE dudz_dddz


!--------------------------------------
  SUBROUTINE zfilter ( vv, isw )
!--------------------------------------
!     z-derivative of f

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: vv

    integer, intent(in)  ::  isw

    complex(kind=DP), dimension(:,:,:), allocatable :: ww

    real(kind=DP) :: aaa

    integer  ::  mx, my, iz, im


      allocate( ww(-nx:nx,0:ny,-nz-nb:nz-1+nb) )


        aaa   = 1._DP


        do iz = -nz, nz-1
          do my = 0, ny
            do mx = -nx, nx
              ww(mx,my,iz) = vv(mx,my,iz)
            end do
          end do
        end do


        call bndry_bound ( ww )


        if( rankz == nprocz-1 ) then

          if( isw >= 0 ) then        ! symmetric

            do my = 0, ny
              do mx = -nx, nx
                ww(mx,my,nz  ) = ( 4._DP*ww(mx,my,nz-1) - ww(mx,my,nz-2) ) / 3._DP
                ww(mx,my,nz+1) =   ww(mx,my,nz-1)
              end do
            end do

          else if( isw < 0 ) then    ! anti-symmetric

            do my = 0, ny
              do mx = -nx, nx
                ww(mx,my,nz  ) = ( 0._DP, 0._DP)
                ww(mx,my,nz+1) = - ww(mx,my,nz-1)
              end do
            end do

          end if

        end if


        if( rankz /= 0 ) then

          do iz = -nz, nz-1
            do my = 0, ny
              do mx = -nx, nx
                vv(mx,my,iz) =                                &
                          ( 1._DP - aaa )  * ww(mx,my,iz  )   &
                        + aaa * ( -          ww(mx,my,iz+2)   &
                                  +  4._DP * ww(mx,my,iz+1)   &
                                  + 10._DP * ww(mx,my,iz  )   &
                                  +  4._DP * ww(mx,my,iz-1)   &
                                  -          ww(mx,my,iz-2) ) &
                               / 16._DP
              end do
            end do
          end do

        else

          do iz = -nz+2, nz-1
            do my = 0, ny
              do mx = -nx, nx
                vv(mx,my,iz) =                                &
                          ( 1._DP - aaa )  * ww(mx,my,iz  )   &
                        + aaa * ( -          ww(mx,my,iz+2)   &
                                  +  4._DP * ww(mx,my,iz+1)   &
                                  + 10._DP * ww(mx,my,iz  )   &
                                  +  4._DP * ww(mx,my,iz-1)   &
                                  -          ww(mx,my,iz-2) ) &
                               / 16._DP
              end do
            end do
          end do


        end if


       deallocate( ww )


  END SUBROUTINE zfilter


!--------------------------------------
  SUBROUTINE iono_dns( idns, icpr, iphi, didns )
!--------------------------------------
!     increment within a time step

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny)  :: idns, icpr, iphi, didns

    complex(kind=DP), &
      dimension(-nx:nx,0:ny)  :: ipsnb


      if( trim(calc_type) == "nonlinear" ) then

        call pssn_brackets( iphi, idns, ipsnb, 1 )

      else

        ipsnb(:,:) = ( 0._DP, 0._DP )

      end if


      didns(:,:) = ( icpr(:,:) - 2.d+0 * alpha * idns(:,:) &
                   - ipsnb(:,:) ) * dt


      didns(0,0) = ( 0._DP, 0._DP )


  END SUBROUTINE iono_dns


END MODULE RMHDS_advnc
