MODULE RMHDS_pssn
!-------------------------------------------------------------------------------
!
!    Calculation of Poisson bracket
!
!      created by Maeyama
!      modified by THW for aurora spectral code
!
!
!      FFTW:   1D FFT => 2D FFT
!      OpenMP: kx,ky-loop parallelization
!
!-------------------------------------------------------------------------------
  use RMHDS_header

  implicit none
  include "aslfftw3.f"
!  include "fftw3.f"
  private

  public  pssn_brackets, pssn_divgrad, pssn_derivative

  integer(kind=DP), save :: plan_backward_x, plan_forward_x
  integer(kind=DP), save :: plan_backward_y, plan_forward_y
  complex(kind=DP), dimension(0:nny/2,0:nnx-1) :: wwkk

  integer(kind=DP), save :: plan_backward_2d, plan_forward_2d

  integer, save :: iflg
  data iflg / 0 /


 CONTAINS


SUBROUTINE pssn_brackets( fk, pk, pb, nm )
!-------------------------------------------------------------------------------
!
!    Calculate time-differential term
!
!-------------------------------------------------------------------------------
  use RMHDS_header

  implicit none
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(in) :: fk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(in) :: pk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(out) :: pb
  integer, intent(in) :: nm

  complex(kind=DP), dimension(-nx:nx,0:ny) :: ikxf, ikyf, ikxp, ikyp
  real(kind=DP), dimension(0:nny-1,0:nnx-1) :: dfdx, dfdy, dpdx, dpdy, pbxy
  integer :: mx, my, ix, iy, im


!%%% Initialize FFTW %%%
    if (iflg == 0) then
      iflg = 1
      call exb_fft_pre
    end if


!%%% Calculate Poisson bracket %%%
!$OMP parallel default(none) &
!$OMP shared(kx,ky,fk,pk,ikxf,ikyf,ikxp,ikyp,dfdx,dfdy,dpdx,dpdy,pbxy,pb,nm) &
!$OMP private(mx,my,im,ix,iy)
    do im = 0, nm-1

      call pssn_derivative( fk(:,:,im), dfdx, dfdy )
      call pssn_derivative( pk(:,:,im), dpdx, dpdy )

!!$OMP do
!      do my = 0, ny
!        do mx = -nx, nx
!          ikxf(mx,my) = ui * kx(mx) * fk(mx,my,im)
!          ikyf(mx,my) = ui * ky(my) * fk(mx,my,im)
!          ikxp(mx,my) = ui * kx(mx) * pk(mx,my,im)
!          ikyp(mx,my) = ui * ky(my) * pk(mx,my,im)
!        end do
!      end do
!!$OMP end do
!      call exb_fft_backward(ikxf, dfdx)
!      call exb_fft_backward(ikyf, dfdy)
!      call exb_fft_backward(ikxp, dpdx)
!      call exb_fft_backward(ikyp, dpdy)

!$OMP do
      do ix = 0, nnx-1
        do iy = 0, nny-1
          pbxy(iy,ix) = - dpdx(iy,ix) * dfdy(iy,ix) + dpdy(iy,ix) * dfdx(iy,ix)
        end do
      end do
!$OMP end do
      call exb_fft_forward(pbxy, pb(:,:,im))
    end do
!$OMP end parallel


END SUBROUTINE pssn_brackets



SUBROUTINE pssn_derivative( gk, dgdx, dgdy )
!-------------------------------------------------------------------------------
!   Compute spatial derivative of g in real space
!-------------------------------------------------------------------------------

  use RMHDS_header

  implicit none

  complex(kind=DP), dimension(-nx:nx,0:ny),  intent(in)  :: gk
  real(kind=DP), dimension(0:nny-1,0:nnx-1), intent(out) :: dgdx, dgdy

  complex(kind=DP), dimension(-nx:nx,0:ny) :: ikxg, ikyg
  integer :: mx, my, ix, iy


!%%% Initialize FFTW %%%
    if (iflg == 0) then
      iflg = 1
      call exb_fft_pre
    end if


!$OMP do
      do my = 0, ny
        do mx = -nx, nx
          ikxg(mx,my) = ui * kx(mx) * gk(mx,my)
          ikyg(mx,my) = ui * ky(my) * gk(mx,my)
        end do
      end do
!$OMP end do
      call exb_fft_backward( ikxg, dgdx )
      call exb_fft_backward( ikyg, dgdy )


END SUBROUTINE pssn_derivative


SUBROUTINE pssn_divgrad( fk, pk, dg, nm )
!-------------------------------------------------------------------------------
!
!    Calculate time-differential term
!
!-------------------------------------------------------------------------------
  use RMHDS_header

  implicit none
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(in) :: fk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(in) :: pk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(out) :: dg
  integer, intent(in) :: nm

  complex(kind=DP), dimension(-nx:nx,0:ny) :: ikxg, ikyg, ikxp, ikyp
  real(kind=DP), dimension(0:nny-1,0:nnx-1) :: fdpdx, fdpdy, dpdx, dpdy, reff
  integer :: mx, my, ix, iy, im


!%%% Initialize FFTW %%%
    if (iflg == 0) then
      iflg = 1
      call exb_fft_pre
    end if


!%%% Calculate Poisson bracket %%%
!$OMP parallel default(none) &
!$OMP shared(kx,ky,fk,pk,ikxg,ikyg,ikxp,ikyp,fdpdx,fdpdy,dpdx,dpdy,reff,dg,nm) &
!$OMP private(mx,my,im,ix,iy)
    do im = 0, nm-1

!$OMP do
      do my = 0, ny
        do mx = -nx, nx
          ikxp(mx,my) = ui * kx(mx) * pk(mx,my,im)
          ikyp(mx,my) = ui * ky(my) * pk(mx,my,im)
        end do
      end do
!$OMP end do
      call exb_fft_backward(fk,   reff)
      call exb_fft_backward(ikxp, dpdx)
      call exb_fft_backward(ikyp, dpdy)

!$OMP do
      do ix = 0, nnx-1
        do iy = 0, nny-1
          fdpdx(iy,ix) = reff(iy,ix) * dpdx(iy,ix)
          fdpdy(iy,ix) = reff(iy,ix) * dpdy(iy,ix)
        end do
      end do
!$OMP end do
      call exb_fft_forward(fdpdx, ikxg)
      call exb_fft_forward(fdpdy, ikyg)
!$OMP do
      do my = 0, ny
        do mx = -nx, nx
          dg(mx,my,im) = ui * ( kx(mx) * ikxg(mx,my) + ky(my) * ikyg(mx,my) )
        end do
      end do
!$OMP end do

    end do
!$OMP end parallel


END SUBROUTINE pssn_divgrad


SUBROUTINE exb_fft_pre
!--------------------------------------
!  Initialization of FFT
  complex(kind=DP), dimension(0:nnx-1) :: wkx1, wkx2
  complex(kind=DP), dimension(0:nny/2) :: wky1
  real(kind=DP), dimension(0:nny-1) :: wky2

!  complex(kind=DP), dimension(0:nnx-1,0:nny/2) :: wkxy1
!  real(kind=DP), dimension(0:nnx-1,0:nny-1) :: wkxy2

  complex(kind=DP), dimension(0:nny/2,0:nnx-1) :: wkxy1
  real(kind=DP), dimension(0:nny-1,0:nny-1) :: wkxy2

!    call dfftw_plan_dft_1d(plan_backward_x,  &
!                           nnx,               &
!                           wkx1, wkx2,       &
!                           FFTW_BACKWARD,    &
!                           FFTW_MEASURE)
!    call dfftw_plan_dft_c2r_1d(plan_backward_y,  &
!                               nny,               &
!                               wky1, wky2,       &
!                               FFTW_MEASURE)
!    call dfftw_plan_dft_r2c_1d(plan_forward_y,   &
!                               nny,               &
!                               wky2, wky1,       &
!                               FFTW_MEASURE)
!    call dfftw_plan_dft_1d(plan_forward_x,   &
!                           nnx,               &
!                           wkx2, wkx1,       &
!                           FFTW_FORWARD,     &
!                           FFTW_MEASURE)

    call dfftw_plan_dft_c2r_2d(plan_backward_2d, &
                               nny, nnx,         &
                               wkxy1, wkxy2,     &
                               FFTW_MEASURE)
    call dfftw_plan_dft_r2c_2d(plan_forward_2d,  &
                               nny, nnx,         &
                               wkxy1, wkxy2,     &
                               FFTW_MEASURE)

END SUBROUTINE exb_fft_pre
  

SUBROUTINE exb_fft_backward(ww, wwxy)
!--------------------------------------
!  Execution of FFT
  complex(kind=DP), dimension(-nx:nx,0:ny), intent(in)   :: ww
  real(kind=DP), dimension(0:nny-1,0:nnx-1), intent(out) :: wwxy

  complex(kind=DP), dimension(0:nny/2,0:nnx-1) :: w1
  integer :: mx, my

    w1(:,:) = (0._DP, 0._DP)

    do my = 0, ny
      w1(my,0:nx)         = ww(0:nx,my)
      w1(my,nnx-nx:nnx-1) = ww(-nx:-1,my)
    end do
      call dfftw_execute_dft_c2r(plan_backward_2d, w1, wwxy)

END SUBROUTINE exb_fft_backward


SUBROUTINE exb_fft_forward(wwxy, ww)
!--------------------------------------
!  Execution of FFT
  real(kind=DP), dimension(0:nny-1,0:nnx-1), intent(in) :: wwxy
  complex(kind=DP), dimension(-nx:nx,0:ny), intent(out) :: ww

  complex(kind=DP), dimension(0:nny/2,0:nnx-1) :: w1
  real(kind=DP) :: ceff_norm
  integer :: mx, my

    ceff_norm = 1._DP / real(nnx * nny, kind=DP)

      call dfftw_execute_dft_r2c(plan_forward_2d, wwxy, w1)

    do my = 0, ny
      ww(0:nx,my)   = w1(my,0:nx)         * ceff_norm
      ww(-nx:-1,my) = w1(my,nnx-nx:nnx-1) * ceff_norm
    end do

    my = 0
      do mx = 1, nx
        ww(-mx,-my) = conjg(ww(mx,my))
      end do

END SUBROUTINE exb_fft_forward


END MODULE RMHDS_pssn
