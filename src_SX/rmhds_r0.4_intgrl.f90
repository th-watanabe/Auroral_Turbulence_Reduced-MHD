MODULE RMHDS_intgrl
!-------------------------------------------------------------------------------
!
!    Flux surface and field line averages
!
!-------------------------------------------------------------------------------

  use RMHDS_header
  use RMHDS_mpienv

  implicit none

  private

  public   intgrl_zeta

    INTERFACE intgrl_zeta
      module procedure  intgrl_zeta_r, intgrl_zeta_z
    END INTERFACE


CONTAINS


!--------------------------------------
  SUBROUTINE intgrl_zeta_r ( wn, dpara, wa )
!--------------------------------------
!     average of a complex variable wn in the z space

    real(kind=DP), intent(in), dimension(:,:,:)  :: wn
    real(kind=DP), intent(in), dimension(:)      :: dpara

    real(kind=DP), intent(out), dimension(:,:)   :: wa

    real(kind=DP), dimension(:,:), allocatable   :: ww

    integer  ::  mx, my, iz
    integer  ::  nx0, ny0, nz0

      nx0 = size(wn, 1) 
      ny0 = size(wn, 2) 
      nz0 = size(wn, 3) 

      allocate( ww(1:nx0,1:ny0) )

!      wa   = ( 0._DP, 0._DP )
!      ww   = ( 0._DP, 0._DP )
      wa   = 0._DP
      ww   = 0._DP

      do iz = 1, nz0
        do my = 1, ny0
          do mx = 1, nx0
            ww(mx,my)   = ww(mx,my) + wn(mx,my,iz) * dpara(iz)
          end do
        end do
      end do

      call MPI_Allreduce( ww, wa, nx0*ny0, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, MPI_COMM_WORLD, ierr_mpi )

      deallocate( ww )


  END SUBROUTINE intgrl_zeta_r


!--------------------------------------
  SUBROUTINE intgrl_zeta_z ( wn, dpara, wa )
!--------------------------------------
!     average of a complex variable wn in the z space

    complex(kind=DP), intent(in), dimension(:,:,:)  :: wn
    real(kind=DP), intent(in), dimension(:)      :: dpara

    complex(kind=DP), intent(out), dimension(:,:)   :: wa

    complex(kind=DP), dimension(:,:), allocatable   :: ww

    integer  ::  mx, my, iz
    integer  ::  nx0, ny0, nz0

      nx0 = size(wn, 1) 
      ny0 = size(wn, 2) 
      nz0 = size(wn, 3) 

      allocate( ww(1:nx0,1:ny0) )

      wa   = ( 0._DP, 0._DP )
      ww   = ( 0._DP, 0._DP )

      do iz = 1, nz0
        do my = 1, ny0
          do mx = 1, nx0
            ww(mx,my)   = ww(mx,my) + wn(mx,my,iz) * dpara(iz)
          end do
        end do
      end do

      call MPI_Allreduce( ww, wa, nx0*ny0, MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, MPI_COMM_WORLD, ierr_mpi )

      deallocate( ww )


  END SUBROUTINE intgrl_zeta_z


END MODULE RMHDS_intgrl
