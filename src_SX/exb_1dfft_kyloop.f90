MODULE exb
!-------------------------------------------------------------------------------
!
!    Calculation of Poisson bracket
!
!      FFTW:   1D FFT
!      OpenMP: kx,ky-loop parallelization
!
!-------------------------------------------------------------------------------
  use parameters
  implicit none
  include "fftw3.f"
  private

  public  exb_poisson_bracket

  integer(kind=DP), save :: plan_backward_x, plan_forward_x
  integer(kind=DP), save :: plan_backward_y, plan_forward_y
  complex(kind=DP), dimension(0:ny/2,0:nx-1) :: wwkk

 CONTAINS


SUBROUTINE exb_poisson_bracket(fk, pk, pb)
!-------------------------------------------------------------------------------
!
!    Calculate time-differential term
!
!-------------------------------------------------------------------------------
  use parameters
  use geometry, only : kx, ky, j0
  implicit none
  complex(kind=DP), dimension(-nkx:nkx,0:nky,0:nm-1), intent(in) :: fk
  complex(kind=DP), dimension(-nkx:nkx,0:nky),        intent(in) :: pk
  complex(kind=DP), dimension(-nkx:nkx,0:nky,0:nm-1), intent(out) :: pb

  complex(kind=DP), dimension(-nkx:nkx,0:nky) :: ikxf, ikyf, ikxp, ikyp
  real(kind=DP), dimension(0:ny-1,0:nx-1) :: dfdx, dfdy, dpdx, dpdy, pbxy
  integer :: mx, my, ix, iy, im

  integer, save :: iflg
  data iflg / 0 /


!%%% Initialize FFTW %%%
    if (iflg == 0) then
      iflg = 1
      call exb_fft_pre
    end if


!%%% Calculate Poisson bracket %%%
!$OMP parallel default(none) &
!$OMP shared(kx,ky,j0,fk,pk,ikxf,ikyf,ikxp,ikyp,dfdx,dfdy,dpdx,dpdy,pbxy,pb) &
!$OMP private(mx,my,im,ix,iy)
    do im = 0, nm-1
!$OMP do
      do my = 0, nky
        do mx = -nkx, nkx
          ikxf(mx,my) = ci * kx(mx) * fk(mx,my,im)
          ikyf(mx,my) = ci * ky(my) * fk(mx,my,im)
          ikxp(mx,my) = ci * kx(mx) * j0(mx,my,im) * pk(mx,my)
          ikyp(mx,my) = ci * ky(my) * j0(mx,my,im) * pk(mx,my)
        end do
      end do
!$OMP end do
      call exb_fft_backward(ikxf, dfdx)
      call exb_fft_backward(ikyf, dfdy)
      call exb_fft_backward(ikxp, dpdx)
      call exb_fft_backward(ikyp, dpdy)

!$OMP do
      do ix = 0, nx-1
        do iy = 0, ny-1
          pbxy(iy,ix) = - dpdx(iy,ix) * dfdy(iy,ix) + dpdy(iy,ix) * dfdx(iy,ix)
        end do
      end do
!$OMP end do
      call exb_fft_forward(pbxy, pb(:,:,im))
    end do
!$OMP end parallel


END SUBROUTINE exb_poisson_bracket


SUBROUTINE exb_fft_pre
!--------------------------------------
!  Initialization of FFT
  complex(kind=DP), dimension(0:nx-1) :: wkx1, wkx2
  complex(kind=DP), dimension(0:ny/2) :: wky1
  real(kind=DP), dimension(0:ny-1) :: wky2

    call dfftw_plan_dft_1d(plan_backward_x,  &
                           nx,               &
                           wkx1, wkx2,       &
                           FFTW_BACKWARD,    &
                           FFTW_MEASURE)
    call dfftw_plan_dft_c2r_1d(plan_backward_y,  &
                               ny,               &
                               wky1, wky2,       &
                               FFTW_MEASURE)
    call dfftw_plan_dft_r2c_1d(plan_forward_y,   &
                               ny,               &
                               wky2, wky1,       &
                               FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan_forward_x,   &
                           nx,               &
                           wkx2, wkx1,       &
                           FFTW_FORWARD,     &
                           FFTW_MEASURE)

END SUBROUTINE exb_fft_pre
  

SUBROUTINE exb_fft_backward(ww, wwxy)
!--------------------------------------
!  Execution of FFT
  complex(kind=DP), dimension(-nkx:nkx,0:nky), intent(in) :: ww
  real(kind=DP), dimension(0:ny-1,0:nx-1), intent(out)    :: wwxy

  complex(kind=DP), dimension(0:nx-1) :: w1, w2
  integer :: mx, my

    w1(:) = (0._DP, 0._DP)
!$OMP workshare
    wwkk(:,:) = (0._DP, 0._DP)
!$OMP end workshare
!$OMP do
    do my = 0, nky
      w1(0:nkx) = ww(0:nkx,my)
      w1(nx-nkx:nx-1) = ww(-nkx:-1,my)
      call dfftw_execute_dft(plan_backward_x, w1, w2)
      wwkk(my,0:nx-1) = w2(0:nx-1)
    end do
!$OMP end do

!$OMP do
    do mx = 0, nx-1
      call dfftw_execute_dft_c2r(plan_backward_y, wwkk(:,mx), wwxy(:,mx))
    end do
!$OMP end do

END SUBROUTINE exb_fft_backward


SUBROUTINE exb_fft_forward(wwxy, ww)
!--------------------------------------
!  Execution of FFT
  real(kind=DP), dimension(0:ny-1,0:nx-1), intent(in)      :: wwxy
  complex(kind=DP), dimension(-nkx:nkx,0:nky), intent(out) :: ww

  complex(kind=DP), dimension(0:nx-1) :: w1, w2
  real(kind=DP) :: ceff_norm
  integer :: mx, my

    ceff_norm = 1._DP / real(nx * ny, kind=DP)

!$OMP do
    do mx = 0, nx-1
      call dfftw_execute_dft_r2c(plan_forward_y, wwxy(:,mx), wwkk(:,mx))
    end do
!$OMP end do
!$OMP do
    do my = 0, nky
      w2(0:nx-1) = wwkk(my,0:nx-1)
      call dfftw_execute_dft(plan_forward_x, w2, w1)
      ww(0:nkx,my) = w1(0:nkx) * ceff_norm
      ww(-nkx:-1,my) = w1(nx-nkx:nx-1) * ceff_norm
    end do
!$OMP end do
    my = 0
      do mx = 1, nkx
        ww(-mx,-my) = conjg(ww(mx,my))
      end do

END SUBROUTINE exb_fft_forward


END MODULE exb
