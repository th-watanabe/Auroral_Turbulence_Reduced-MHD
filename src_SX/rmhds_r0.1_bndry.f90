MODULE RMHDS_bndry
!-------------------------------------------------------------------------------
!
!    Some useful tools and tips
!
!      GKV-plus r0.3 ( T.-H.Watanabe, Jun 2011)
!      GKV-plus r1.0 for Multi-K species ( M. Nakata, November 2011)
!
!-------------------------------------------------------------------------------

  use RMHDS_header
  use RMHDS_mpienv

  implicit none

  private

  public   bndry_bound


CONTAINS


!--------------------------------------
  SUBROUTINE bndry_bound ( ww )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nb:nz-1+nb)   :: ww

    complex(kind=DP), dimension(:,:,:), allocatable :: zb1, zb2

! --- local variables

    integer  ::  mx, my, ib
    integer  ::  slngz


      slngz  = (2*nx+1)*(ny+1)*nb


      allocate( zb1(-nx:nx,0:ny,nb*2) )
      allocate( zb2(-nx:nx,0:ny,nb*2) )


      do ib = 1, nb
        do my = 0, ny
          do mx = -nx, nx
            zb1(mx,my,ib)    = ww(mx,my,-nz   +ib-1)
            zb1(mx,my,nb+ib) = ww(mx,my, nz-nb+ib-1)
          end do
        end do
      end do

      call MPI_sendrecv( zb1(-nx,0,1),    slngz, MPI_DOUBLE_COMPLEX, izdn, 1, &
                         zb2(-nx,0,nb+1), slngz, MPI_DOUBLE_COMPLEX, izup, 1, &
                         MPI_COMM_WORLD, status, ierr_mpi )

      call MPI_sendrecv( zb1(-nx,0,nb+1), slngz, MPI_DOUBLE_COMPLEX, izup, 2, &
                         zb2(-nx,0,1),    slngz, MPI_DOUBLE_COMPLEX, izdn, 2, &
                         MPI_COMM_WORLD, status, ierr_mpi )


      if( rankz /= 0 ) then

        do ib = 1, nb
          do my = 0, ny
            do mx = -nx, nx
              ww(mx,my,-nz-nb+ib-1) = zb2(mx,my,ib)
            end do
          end do
        end do

      else

        do ib = 1, nb
          do my = 0, ny
            do mx = -nx, nx
              ww(mx,my,-nz-nb+ib-1) = ( 0._DP, 0._DP )
            end do
          end do
        end do

      end if


      if( rankz /= nprocz-1 ) then

        do ib = 1, nb
          do my = 0, ny
            do mx = -nx, nx
              ww(mx,my, nz+ib-1) = zb2(mx,my,nb+ib)
            end do
          end do
        end do

      else

        do ib = 1, nb
          do my = 0, ny
            do mx = -nx, nx
              ww(mx,my, nz+ib-1) = ( 0._DP, 0._DP )
            end do
          end do
        end do

      end if



  END SUBROUTINE bndry_bound


END MODULE RMHDS_bndry
