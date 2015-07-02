*> \brief \b DLAMCH
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAMCH determines double precision machine parameters.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] CMACH
*> \verbatim
*>          Specifies the value to be returned by DLAMCH:
*>          = 1  --> 'E' or 'e',   DLAMCH := eps
*>          = 2  --> 'S' or 's ,   DLAMCH := sfmin
*>          = 3  --> 'B' or 'b',   DLAMCH := base
*>          = 4  --> 'P' or 'p',   DLAMCH := eps*base
*>          = 5  --> 'N' or 'n',   DLAMCH := t
*>          = 6  --> 'R' or 'r',   DLAMCH := rnd
*>          = 7  --> 'M' or 'm',   DLAMCH := emin
*>          = 8  --> 'U' or 'u',   DLAMCH := rmin
*>          = 9  --> 'L' or 'l',   DLAMCH := emax
*>          = 10 --> 'O' or 'o',   DLAMCH := rmax
*>          where
*>          eps   = relative machine precision
*>          sfmin = safe minimum, such that 1/sfmin does not overflow
*>          base  = base of the machine
*>          prec  = eps*base
*>          t     = number of (base) digits in the mantissa
*>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*>          emin  = minimum exponent before (gradual) underflow
*>          rmin  = underflow threshold - base**(emin-1)
*>          emax  = largest exponent before overflow
*>          rmax  = overflow threshold  - (base**emax)*(1-eps)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!$acc routine seq
*
*  -- LAPACK auxiliary routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
*      CHARACTER          CMACH
      INTEGER            CMACH
*     ..
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT,
     $                   MINEXPONENT, RADIX, TINY
*     ..
*     .. Executable Statements ..
*
*
*     Assume rounding, not chopping. Always.
*
      RND = ONE
*
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
*'
*
*!      IF( LSAME( CMACH, 'E' ) ) THEN
      IF( CMACH == 1 ) THEN
         RMACH = EPS
*!      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
      ELSE IF( CMACH == 2 ) THEN
         SFMIN = TINY(ZERO)
         SMALL = ONE / HUGE(ZERO)
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
         RMACH = SFMIN
*!      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
      ELSE IF( CMACH == 3 ) THEN
         RMACH = RADIX(ZERO)
*!      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
      ELSE IF( CMACH == 4 ) THEN
         RMACH = EPS * RADIX(ZERO)
*!      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
      ELSE IF( CMACH == 5 ) THEN
         RMACH = DIGITS(ZERO)
*!      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
      ELSE IF( CMACH == 6 ) THEN
         RMACH = RND
*!      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
      ELSE IF( CMACH == 7 ) THEN
         RMACH = MINEXPONENT(ZERO)
*!      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
      ELSE IF( CMACH == 8 ) THEN
         RMACH = tiny(zero)
*!      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
      ELSE IF( CMACH == 9 ) THEN
         RMACH = MAXEXPONENT(ZERO)
*!      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
      ELSE IF( CMACH == 10 ) THEN
         RMACH = HUGE(ZERO)
      ELSE
         RMACH = ZERO
      END IF
*
      DLAMCH = RMACH
      RETURN
*
*     End of DLAMCH
*
      END
************************************************************************
*> \brief \b DLAMC3
*> \details
*> \b Purpose:
*> \verbatim
*> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*> the addition of  A  and  B ,  for use in situations where optimizers
*> might hold one of these in a register.
*> \endverbatim
*> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*> \date November 2011
*> \ingroup auxOTHERauxiliary
*>
*> \param[in] A
*> \verbatim
*>          A is a DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is a DOUBLE PRECISION
*>          The values A and B.
*> \endverbatim
*>
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.4.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
