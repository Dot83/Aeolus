      subroutine r8_normal_ab( a, b, seed,r8 )

!*****************************************************************************80
!
!! R8_NORMAL_AB returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the mean of the PDF.
!
!    Input, real ( kind = 8 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_AB, a sample of the normal PDF.
!
      implicit none

      real ( kind = 8 ) a,r8
      real ( kind = 8 ) b
      real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
      real ( kind = 8 ) r1
      real ( kind = 8 ) r2
!      real ( kind = 8 ) r8_normal_ab
!      real ( kind = 8 ) r8_uniform_01
      integer ( kind = 4 ) seed
      real ( kind = 8 ) x

      call r8_uniform_01( seed, r1 )
      call r8_uniform_01( seed, r2 )
      x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )

      r8 = a + b * x
!     r8_normal_ab = a + b * x
      return
      end subroutine
