!!!
!
! -- Inria
! -- (C) Copyright 2012
!
! This software is a computer program whose purpose is to process
! Matrices Over Runtime Systems @ Exascale (MORSE). More information
! can be found on the following website: http://www.inria.fr/en/teams/morse.
! 
! This software is governed by the CeCILL-B license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL-B
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
! 
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
! 
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
! 
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-B license and that you accept its terms.
!
!!!

      SUBROUTINE CLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, UPLO
      INTEGER            LDA, N
      REAL               AMAX, SCOND
*     ..
*     .. Array Arguments ..
      REAL               S( * )
      COMPLEX            A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  CLAQSY equilibrates a symmetric matrix A using the scaling factors
*  in the vector S.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n by n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if EQUED = 'Y', the equilibrated matrix:
*          diag(S) * A * diag(S).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(N,1).
*
*  S       (input) REAL array, dimension (N)
*          The scale factors for A.
*
*  SCOND   (input) REAL
*          Ratio of the smallest S(i) to the largest S(i).
*
*  AMAX    (input) REAL
*          Absolute value of largest matrix entry.
*
*  EQUED   (output) CHARACTER*1
*          Specifies whether or not equilibration was done.
*          = 'N':  No equilibration.
*          = 'Y':  Equilibration was done, i.e., A has been replaced by
*                  diag(S) * A * diag(S).
*
*  Internal Parameters
*  ===================
*
*  THRESH is a threshold value used to decide if scaling should be done
*  based on the ratio of the scaling factors.  If SCOND < THRESH,
*  scaling is done.
*
*  LARGE and SMALL are threshold values used to decide if scaling should
*  be done based on the absolute size of the largest matrix element.
*  If AMAX > LARGE or AMAX < SMALL, scaling is done.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, THRESH
      PARAMETER          ( ONE = 1.0E+0, THRESH = 0.1E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               CJ, LARGE, SMALL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
*
*     Initialize LARGE and SMALL.
*
      SMALL = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' )
      LARGE = ONE / SMALL
*
      IF( SCOND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) THEN
*
*        No equilibration
*
         EQUED = 'N'
      ELSE
*
*        Replace A by diag(S) * A * diag(S).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Upper triangle of A is stored.
*
            DO 20 J = 1, N
               CJ = S( J )
               DO 10 I = 1, J
                  A( I, J ) = CJ*S( I )*A( I, J )
   10          CONTINUE
   20       CONTINUE
         ELSE
*
*           Lower triangle of A is stored.
*
            DO 40 J = 1, N
               CJ = S( J )
               DO 30 I = J, N
                  A( I, J ) = CJ*S( I )*A( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         EQUED = 'Y'
      END IF
*
      RETURN
*
*     End of CLAQSY
*
      END
