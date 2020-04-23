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

      SUBROUTINE ZERRQR( PATH, NUNIT )
*
      INCLUDE 'morse_fortran.h'
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  ZERRQR tests the error exits for the COMPLEX*16 routines
*  that use the QR decomposition of a general matrix.
*
*  Arguments
*  =========
*
*  PATH    (input) CHARACTER*3
*          The LAPACK path name for the routines to be tested.
*
*  NUNIT   (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX
      PARAMETER          ( NMAX = 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
*     ..
*     .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ),
     $                   W( NMAX ), X( NMAX )
      INTEGER            HT( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, ZGEQR2, ZGEQRF, ZUNG2R,
     $                   ZUNGQR, ZUNM2R, ZUNMQR
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
*     Disable MORSE warnings/errors
* 
      CALL MORSE_DISABLE( MORSE_WARNINGS, INFO )
      CALL MORSE_DISABLE( MORSE_ERRORS,   INFO )
*
*     Set the variables to innocuous values.
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ),
     $                  -1.D0 / DBLE( I+J ) )
            AF( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ),
     $                   -1.D0 / DBLE( I+J ) )
   10    CONTINUE
         B( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
   20 CONTINUE
      OK = .TRUE.
*
*     Allocate HT
*
      CALL MORSE_ALLOC_WORKSPACE_ZGEQRF( 2, 2, HT, INFO )
*
*
*     Error exits for QR factorization
*
*     ZGEQRF
*
      SRNAMT = 'ZGEQRF'
      INFOT = 1
      CALL MORSE_ZGEQRF( -1, 0, A, 1, HT, INFO )
      CALL CHKXER( 'ZGEQRF', INFOT, NOUT, INFO, OK )
      INFOT = 2
      CALL MORSE_ZGEQRF( 0, -1, A, 1, HT, INFO )
      CALL CHKXER( 'ZGEQRF', INFOT, NOUT, INFO, OK )
      INFOT = 4
      CALL MORSE_ZGEQRF( 2, 1, A, 1, HT, INFO )
      CALL CHKXER( 'ZGEQRF', INFOT, NOUT, INFO, OK )
*
*     ZGEQRS
*
      SRNAMT = 'ZGEQRS'
      INFOT = 1
      CALL MORSE_ZGEQRS( -1, 0, 0, A, 1, X, B, 1, INFO )
      CALL CHKXER( 'ZGEQRS', INFOT, NOUT, INFO, OK )
      INFOT = 2
      CALL MORSE_ZGEQRS( 0, -1, 0, A, 1, X, B, 1, INFO )
      CALL CHKXER( 'ZGEQRS', INFOT, NOUT, INFO, OK )
      INFOT = 2
      CALL MORSE_ZGEQRS( 1, 2, 0, A, 2, X, B, 2, INFO )
      CALL CHKXER( 'ZGEQRS', INFOT, NOUT, INFO, OK )
      INFOT = 3
      CALL MORSE_ZGEQRS( 0, 0, -1, A, 1, X, B, 1, INFO )
      CALL CHKXER( 'ZGEQRS', INFOT, NOUT, INFO, OK )
      INFOT = 5
      CALL MORSE_ZGEQRS( 2, 1, 0, A, 1, X, B, 2, INFO )
      CALL CHKXER( 'ZGEQRS', INFOT, NOUT, INFO, OK )
      INFOT = 8
      CALL MORSE_ZGEQRS( 2, 1, 0, A, 2, X, B, 1, INFO )
      CALL CHKXER( 'ZGEQRS', INFOT, NOUT, INFO, OK )
*
*     ZUNGQR
*
      SRNAMT = 'ZUNGQR'
      INFOT = 1
      CALL MORSE_ZUNGQR( -1, 0, 0, A, 1, HT, W, 1, INFO )
      CALL CHKXER( 'ZUNGQR', INFOT, NOUT, INFO, OK )
      INFOT = 2
      CALL MORSE_ZUNGQR( 0, -1, 0, A, 1, HT, W, 1, INFO )
      CALL CHKXER( 'ZUNGQR', INFOT, NOUT, INFO, OK )
      INFOT = 2
      CALL MORSE_ZUNGQR( 1, 2, 0, A, 1, HT, W, 2, INFO )
      CALL CHKXER( 'ZUNGQR', INFOT, NOUT, INFO, OK )
      INFOT = 3
      CALL MORSE_ZUNGQR( 0, 0, -1, A, 1, HT, W, 1, INFO )
      CALL CHKXER( 'ZUNGQR', INFOT, NOUT, INFO, OK )
      INFOT = 3
      CALL MORSE_ZUNGQR( 1, 1, 2, A, 1, HT, W, 1, INFO )
      CALL CHKXER( 'ZUNGQR', INFOT, NOUT, INFO, OK )
      INFOT = 5
      CALL MORSE_ZUNGQR( 2, 2, 0, A, 1, HT, W, 2, INFO )
      CALL CHKXER( 'ZUNGQR', INFOT, NOUT, INFO, OK )
      INFOT = 8
      CALL MORSE_ZUNGQR( 2, 2, 0, A, 2, HT, W, 1, INFO )
      CALL CHKXER( 'ZUNGQR', INFOT, NOUT, INFO, OK )
*
*     MORSE_ZUNMQR
*
      SRNAMT = 'ZUNMQR'
      INFOT = 1
      CALL MORSE_ZUNMQR( '/', MORSECONJTRANS, 0, 0, 0, A, 1, HT, AF,
     $                   1, INFO )
      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
      INFOT = 2
      CALL MORSE_ZUNMQR( MORSELEFT, '/', 0, 0, 0, A, 1, HT, AF, 1,
     $                   INFO )
      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
      INFOT = 3
      CALL MORSE_ZUNMQR( MORSELEFT, MORSECONJTRANS, -1, 0, 0, A, 1,
     $                   HT, AF, 1, INFO )
      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
      INFOT = 4
      CALL MORSE_ZUNMQR( MORSELEFT, MORSECONJTRANS, 0, -1, 0, A, 1,
     $                   HT, AF, 1, INFO )
      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
      INFOT = 5
      CALL MORSE_ZUNMQR( MORSELEFT, MORSECONJTRANS, 0, 0, -1, A, 1,
     $                   HT, AF, 1, INFO )
      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
*      INFOT = 5
*      CALL MORSE_ZUNMQR( MORSELEFT, MORSECONJTRANS, 0, 1, 1, A, 1, HT,
*     4                   AF, 1, INFO )
*      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
*      INFOT = 5
*      CALL MORSE_ZUNMQR( MORSELEFT, MORSECONJTRANS, 1, 0, 1, A, 1, HT,
*     4                   AF, 1, INFO )
*      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
*      INFOT = 7
*      CALL MORSE_ZUNMQR( MORSELEFT, MORSECONJTRANS, 2, 1, 0, A, 1, HT,
*     4                   AF, 2, INFO )
*      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
*      INFOT = 7
*      CALL MORSE_ZUNMQR( MORSELEFT, MORSECONJTRANS, 1, 2, 0, A, 1, HT,
*     4                   AF, 1, INFO )
*      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
*      INFOT = 10
*      CALL MORSE_ZUNMQR( MORSELEFT, MORSECONJTRANS, 1, 2, 0, A, 1, HT,
*     4                   AF, 1, INFO )
*      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
*      INFOT = 10
*      CALL MORSE_ZUNMQR( MORSELEFT, MORSECONJTRANS, 2, 1, 0, A, 1, HT,
*     4                   AF, 2, INFO )
*      CALL CHKXER( 'ZUNMQR', INFOT, NOUT, INFO, OK )
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
*     Deallocate HT
*
      CALL MORSE_DEALLOC_HANDLE( HT, INFO )
*
*     Enable MORSE warnings/errors
* 
      CALL MORSE_ENABLE( MORSE_WARNINGS, INFO )
      CALL MORSE_ENABLE( MORSE_ERRORS,   INFO )
*
      RETURN
*
*     End of ZERRQR
*
      END
