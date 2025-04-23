MODULE RMHDS_fld
!-------------------------------------------------------------------------------
!
!    Time integration of the reduced MHD equations
!
!-------------------------------------------------------------------------------

  use RMHDS_header
  use RMHDS_mpienv
  use RMHDS_pssn

  implicit none

  private

  public   fld_calfld


CONTAINS


!--------------------------------------
  SUBROUTINE fld_calfld( omg, psi, dns, phi, cpr )
!--------------------------------------
!     coompute phi and cpr

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: dns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg, psi
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr

    integer :: mx, my, iz

      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            phi(mx,my,iz) = - ksqi(mx,my,iz) * omg(mx,my,iz)
!!! The sign is the same as the previous RMHD analyses 
!!! but differes from the popular definition in GK
            cpr(mx,my,iz) = - ksq (mx,my,iz) * psi(mx,my,iz)
!!!            cpr(mx,my,iz) =   ksq (mx,my,iz) * psi(mx,my,iz)
          end do
        end do
      end do


! ionospheric boundary
      if( rankz == 0 ) then
!          print *, "# iono_phi begin: rankz = ", rankz

!!! j -> phi
!!        call iono_phi( cpr(:,:,-nz), dns(:,:,-nz), phi(:,:,-nz), omg(:,:,-nz) )

! phi -> j
        call iono_cpr( phi(:,:,-nz), dns(:,:,-nz), cpr(:,:,-nz), psi(:,:,-nz) )

!          print *, "# iono_phi end: rankz = ", rankz
      end if


  END SUBROUTINE fld_calfld


!--------------------------------------
  SUBROUTINE iono_phi( icpr, idns, iphi, iomg )
!--------------------------------------
!    compute the ionospheric potential

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny)  :: icpr, idns
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny)  :: iphi, iomg

    integer :: mx, my, iz

      iz = -nz

      do my = 0, ny
        do mx = -nx, nx
          iphi(mx,my) = ksqi(mx,my,iz) / ( mp * idns0 )            &
                      * ( ( -ui * e0 * ( mp*ky(my) - mh*kx(mx) ) &
                            - ksq(mx,my,iz) * dperp ) * idns(mx,my)  &
                          - icpr(mx,my) )
          iomg(mx,my) = - ksq(mx,my,iz) * iphi(mx,my)
        end do
      end do

     
  END SUBROUTINE iono_phi


!--------------------------------------
  SUBROUTINE iono_cpr( iphi, idns, icpr, ipsi )
!--------------------------------------
!    compute the ionospheric potential

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny)  :: iphi, idns
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny)  :: icpr, ipsi

    complex(kind=DP), &
      dimension(-nx:nx,0:ny) :: ipsnb, idvgd

    integer :: mx, my, iz

      iz = -nz

      if( trim(calc_type) == "nonlinear" ) then
        call pssn_brackets( iphi, idns, ipsnb, 1 )
        call pssn_divgrad ( idns, iphi, idvgd, 1 )
      else
        ipsnb(:,:) = ( 0._DP, 0._DP )
        idvgd(:,:) = ( 0._DP, 0._DP )
      end if
      

      do my = 0, ny
        do mx = -nx, nx
          icpr(mx,my) = - mp * idns0 * ksq(mx,my,iz) * iphi(mx,my)  &
                        + ( -ui * e0 * ( mp*ky(my) - mh*kx(mx) )    &
                            - ksq(mx,my,iz) * dperp ) * idns(mx,my) &
                        + mp * idvgd(mx,my) + ipsnb(mx,my)

          ipsi(mx,my) = - ksqi(mx,my,iz) * icpr(mx,my)
        end do
      end do

     
  END SUBROUTINE iono_cpr


END MODULE RMHDS_fld
