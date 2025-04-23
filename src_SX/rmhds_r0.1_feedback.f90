MODULE RMHDS_feedback
!---------------------------------------------------------------------
!     analysis of feedback instability
!           based on the correct treatment of the diffusion term
!
!       Extended for the dipole coordinates
!
!                            by tomo  Jul. 23, 2009.
!---------------------------------------------------------------------

  use RMHDS_header

  implicit none

  private

!------------
! old parameters
!      real*8    zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .          dmeta, dmnu, zet0, dzet
!      common /comprm/ zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .                dmeta, dmnu, zet0, dzet
!------------

  real(kind=DP) :: zl, zet0, dzet

  public feedback_linana


CONTAINS
 
!---------------------------
  SUBROUTINE feedback_linana( kx, ky, zom, egphi, egpsi, egdns, nzg )
!---------------------------

    integer          ::   nzg

    complex(kind=DP), dimension(0:nzg) :: egphi, egpsi
    complex(kind=DP) ::  zom, egdns
    real(kind=DP)    :: kx, ky

    complex(kind=DP) :: zom0


      zl    = lz * 2._DP
      zet0  = pi/9._DP
      dzet  = ( pi*0.5_DP - zet0 ) / real( nzg, kind=DP )


      zom0 = ( 2.d-3, 1.d-4 )

      call solv( zom0, kx, ky, nzg, zom, egphi, egpsi, egdns )


    return

  END SUBROUTINE feedback_linana


!----------------------------------------------
!     solver in Newton method
!
!----------------------------------------------
  SUBROUTINE solv( zom0, kx, ky, nzg, zom, egphi, egpsi, egdns )

      integer          :: nzg
      real(kind=DP)    :: kx, ky
      complex(kind=DP) :: zom, zom0, egphi(0:nzg), egpsi(0:nzg), egdns

      complex(kind=DP) :: f, df, z, znew

      real(kind=DP)    :: eps, eps0

      integer      isw, i


      data     eps0     /  1.d-8 /


        z          = zom0
        znew       = zom0


        i          = 0


   10   continue
          i          = i + 1
          call  fun ( z,  f, 0, kx, ky, nzg, egphi, egpsi, egdns )
          call dfunr( z, df, 0, kx, ky, nzg, egphi, egpsi, egdns )

!          call dfuni( z, df, 0, 
!     .              kx, ky, nzg, egphi, egpsi, egdns )
!                               ----  these two routines give the
!                                     same result
!                                     (Cauchy-Riemann condition) 

          znew       = z - f / df
          eps        = abs( znew - z ) / abs( z )
          z          = znew
          if( i .ge. 30 ) goto 999
        if( eps .gt. eps0 )   goto 10

  999   continue


        zom = z

        call  fun( z,  f, 1, kx, ky, nzg, egphi, egpsi, egdns )

!
!        write(6,*) '   *   solution is obtained ; ', z
!        write(6,*) '   *   convergence check ; ', i, eps
!

      return

  END SUBROUTINE solv


!----------------------------------------------
!     function to be solved
!
!----------------------------------------------
  SUBROUTINE fun( z, f, iflg, kx, ky, nzg, egphi, egpsi, egdns )

      complex(kind=DP) :: z, f
      integer          :: iflg

      integer          :: nzg
      real(kind=DP)    :: kx, ky
      complex(kind=DP), dimension(0:nzg) :: egphi, egpsi
      complex(kind=DP) :: egdns


!      real*8    zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .          dmeta, dmnu, zet0, dzet
!      common /comprm/ zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .                dmeta, dmnu, zet0, dzet

      complex(kind=DP) :: zres
      real(kind=DP)    :: omgr, omgi, theta, kperp
      integer          :: nout, k


        omgr       =  dble(z)
        omgi       =  imag(z)
        kperp      =  sqrt( kx**2 + ky**2 )
        theta      =  atan2( ky, kx )

!
! magnetospheric response Phi/J in the slab geom.
!
!        zres       = - ( cdexp(ui*kpara*zl) + cdexp(-ui*kpara*zl) )
!     .               / ( cdexp(ui*kpara*zl) - cdexp(-ui*kpara*zl) )
!


        call response( z, theta, kperp, nzg, zres, egphi, egpsi )

        egdns = kperp**2 * egpsi(0) / ( ui * z - 2.d+0 * alpha )


        f   = ( z + ui * 2.d+0 * alpha ) * ( 1.d+0 + zres * mp * idns0 ) &
            - kperp * e0 * ( mp * sin(theta) - cos(theta) )              &
            + ui * kperp**2 * dperp



        if( iflg .eq. 1 ) then
           nout    = 6

           write(nout,630) '#------  '
           write(nout,630) '# kx, ky, omgr, omgi, f '
           write(nout,600) '#   ', kx, ky, z, f
           write(nout,630) '#------  '

!           write(nout,610) '# egdns = ', egdns
!           write(nout,630) '# k, egphi, egpsi '
!           do k = 0, nzg
!             write(nout,620) k, egphi(k), egpsi(k)
!           end do

  600        format(  a, 1p, 6e15.7 )
  610        format(  a, 1p, 2e15.7 )
  620        format( i5, 1p, 4e15.7 )
  630        format( a )
        end if


        return

  END SUBROUTINE fun


!----------------------------------------------
!     finite diffence function of f
!
!----------------------------------------------
  SUBROUTINE dfunr( z, df, iflg, kx, ky, nzg, egphi, egpsi, egdns )

      complex(kind=DP) :: z, df
      real(kind=DP)    :: kx, ky
      integer          :: iflg, nzg

      complex(kind=DP), dimension(0:nzg) :: egphi, egpsi
      complex(kind=DP) :: egdns


      complex(kind=DP) :: ff1, ff2, z1, z2
      real(kind=DP)    :: dx, d

      data     dx  /  1.d-10 /


        d          = dble(z) * dx
        z1         = z + dcmplx( d*0.5d0, 0.d0 )
        z2         = z - dcmplx( d*0.5d0, 0.d0 )


        call fun( z1, ff1, iflg, kx, ky, nzg, egphi, egpsi, egdns )
        call fun( z2, ff2, iflg, kx, ky, nzg, egphi, egpsi, egdns )


        df         = ( ff1 - ff2 ) / d


      return


  END SUBROUTINE dfunr


!----------------------------------------------
!     finite diffence function of f
!
!----------------------------------------------
  SUBROUTINE dfuni( z, df, iflg, kx, ky, nzg, egphi, egpsi, egdns )

      complex(kind=DP) :: z, df
      real(kind=DP)    :: kx, ky
      integer          :: iflg, nzg

      complex(kind=DP), dimension(0:nzg) :: egphi, egpsi
      complex(kind=DP) :: egdns

      complex(kind=DP) :: ff1, ff2, z1, z2
      real(kind=DP)    :: dx, d

      data     dx  /  1.d-7 /


        d          = imag(z) * dx
        z1         = z + dcmplx( 0.d0, d*0.5d0 )
        z2         = z - dcmplx( 0.d0, d*0.5d0 )


        call fun( z1, ff1, iflg, kx, ky, nzg, egphi, egpsi, egdns )
        call fun( z2, ff2, iflg, kx, ky, nzg, egphi, egpsi, egdns )


        df         = ( ff1 - ff2 ) / dcmplx( 0.d0, d )


      return


  END SUBROUTINE dfuni


!----------------------------------------------
!     magnetospheric response Phi/J
!
!----------------------------------------------
  SUBROUTINE response( z, theta, kperp, nzg, zres, egphi, egpsi )

      complex(kind=DP) :: z, zres
      real(kind=DP)    :: theta, kperp
      integer          :: nzg

      complex(kind=DP), dimension(0:nzg) :: egphi, egpsi

!      real*8    zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .          dmeta, dmnu, zet0, dzet
!      common /comprm/ zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .                dmeta, dmnu, zet0, dzet

      complex(kind=DP) :: psi, chi, phi, omg, bph
      real(kind=DP)    :: kperp_sq, zet, hph, hnu, dzds, sss, b0, valf

      complex(kind=DP), dimension(4) :: chi_k, chi_w, bph_k, bph_w
      real(kind=DP),    dimension(4) :: ck(4), ci(4)

      real(kind=DP)    :: fc1, fc2
      integer          :: i, j, k, isw


!
! coefficients used in the RK method
        ck(1)   = 1.d0/6.d0
        ck(2)   = 1.d0/3.d0
        ck(3)   = 1.d0/3.d0
        ck(4)   = 1.d0/6.d0

        ci(1)   = 0.d0
        ci(2)   = 0.5d0
        ci(3)   = 0.5d0
        ci(4)   = 1.d0


        chi_k = ( 0.d0, 0.d0 )
        chi_w = ( 0.d0, 0.d0 )
        bph_k = ( 0.d0, 0.d0 )
        bph_w = ( 0.d0, 0.d0 )

! boundary condition on the equator
        chi   = ( 0.d0, 0.d0 )
        bph   = ( 1.d0, 0.d0 )


        k     = nzg
        zet   = zet0 + dzet*dble(k)
        call metric( zl, zet0, zet, hnu, hph, &
                     dzds, sss, b0, valf )

        kperp_sq = kperp**2                    &
                 * ( ( cos(theta) / hph )**2   &
                   + ( sin(theta) / hnu )**2 )

        psi   = - chi / kperp_sq
        omg   = - kperp_sq * bph / b0

        egpsi(nzg) = psi * b0
        egphi(nzg) = bph


        do k = nzg, 1, -1

! Runge-Kutta starts

          do j = 1, 4
            if( j .eq. 1 ) then
              chi_w(j) = chi
              bph_w(j) = bph
            else
              chi_w(j) = chi + ci(j) * chi_k(j-1) * dzet
              bph_w(j) = bph + ci(j) * bph_k(j-1) * dzet
            end if

            if( j .eq. 1 ) then
              zet   = zet0 + dzet*dble(k)
              call metric( zl, zet0, zet, hnu, hph, &
                           dzds, sss, b0, valf )
            else if ( j .eq. 4 ) then
              zet   = zet0 + dzet*dble(k-1)
              call metric( zl, zet0, zet, hnu, hph, &
                           dzds, sss, b0, valf )
            else
              zet   = zet0 + dzet*(dble(k)-0.5d+0)
              call metric( zl, zet0, zet, hnu, hph, &
                           dzds, sss, b0, valf )
            end if

            kperp_sq = kperp**2                    &
                     * ( ( cos(theta) / hph )**2   &
                       + ( sin(theta) / hnu )**2 )

            fc1  = - kperp_sq /( valf**2 * dzds * b0 )
            fc2  = - b0 / ( dzds * kperp_sq )

            chi_k(j) = ui * z * bph_w(j) * fc1
            bph_k(j) = ui * z * chi_w(j) * fc2
          end do

          do j = 1, 4
            chi   = chi + dzet * ck(j)*chi_k(j)
            bph   = bph + dzet * ck(j)*bph_k(j)
          end do

          psi  = - chi / kperp_sq
          omg  = - kperp_sq * bph / b0

          egpsi(k-1) = psi * b0
          egphi(k-1) = bph

        end do

        phi   = - omg / kperp**2

        zres = kperp**2 * phi / chi


      return


  END SUBROUTINE response


!----------------------------------------------
!     metric
!
!----------------------------------------------
!!! SUBROUTINE metric_slab( zl, zet0, zet, hnu, hph, dzds, sss, b0, valf )
  SUBROUTINE metric( zl, zet0, zet, hnu, hph, dzds, sss, b0, valf )

      real(kind=DP) :: zl, zet0, zet
      real(kind=DP) :: hnu, hph, dzds, sss, b0, valf

      real(kind=DP) :: hnu0, hph0, w10, w20, w1, w2, al


        hph   = 1.d0
        hnu   = 1.d0
        dzds  = ( 0.5d0 * pi - zet0 ) / zl
        sss   = zl * ( zet - zet0 ) / ( 0.5d0 * pi - zet0 )
        b0    = 1.d0
        valf  = 1.d0


      return


  END SUBROUTINE metric
!!!  END SUBROUTINE metric_slab


!----------------------------------------------
!     metric
!
!----------------------------------------------
  SUBROUTINE metric_dipole( zl, zet0, zet, hnu, hph, dzds, sss, b0, valf )

      real(kind=DP) :: zl, zet0, zet
      real(kind=DP) :: hnu, hph, dzds, sss, b0, valf

      real(kind=DP) :: hnu0, hph0, w10, w20, w1, w2, al


        w10   = sqrt( 3.d0 ) * cos( zet0 )
        w20   = sqrt( 1.d0 + w10**2 )
        hph0  = sin( zet0 )**3
        hnu0  = sin( zet0 )**3 / w20

        al    = zl * 2.d0 * sqrt( 3.d0 )  &
              / ( w10 * w20 + log( abs( w10 + w20 ) ) )

        w1    = sqrt( 3.d0 ) * cos( zet )
        w2    = sqrt( 1.d0 + w1**2 )

        hph   = sin( zet )**3      / hph0
        hnu   = sin( zet )**3 / w2 / hnu0
        dzds  = 1.d0 / ( al * sin( zet ) * w2 )
        sss   = al / ( 2.d0*sqrt(3.d0) )                     &
                   * ( w10 * w20 + log( abs( w10 + w20 ) )   &
                     - w1  * w2  - log( abs( w1  + w2  ) ) )
        b0    = 1.d0 / ( hnu * hph )
        valf  = 1.d0


      return


  END SUBROUTINE metric_dipole


END MODULE RMHDS_feedback
