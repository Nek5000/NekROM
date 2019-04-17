c-----------------------------------------------------------------------
      subroutine dgetrf( m, n, a, lda, ipiv, info )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      integer            info, lda, m, n
*     ..
*     .. Array Arguments ..
      integer            ipiv( * )
      double precision   a( lda, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      integer            i, iinfo, j, jb, nb
*     ..
*     .. External Subroutines ..
      external           dgemm, dgetf2, dlaswp, dtrsm, xerbla
*     ..
*     .. External Functions ..
      integer            ilaenv
      external           ilaenv
*     ..
*     .. Intrinsic Functions ..
      intrinsic          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      if( m.lt.0 ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, m ) ) then
         info = -4
      end if
      if( info.ne.0 ) then
         call xerbla( 'DGETRF', -info )
         return
      end if
*
*     Quick return if possible
*
      if( m.eq.0 .or. n.eq.0 )
     $   return
*
*     Determine the block size for this environment.
*
      nb = ilaenv( 1, 'DGETRF', ' ', m, n, -1, -1 )
      if( nb.le.1 .or. nb.ge.min( m, n ) ) then
*
*        Use unblocked code.
*
         call dgetf2( m, n, a, lda, ipiv, info )
      else
*
*        Use blocked code.
*
         do 20 j = 1, min( m, n ), nb
            jb = min( min( m, n )-j+1, nb )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            call dgetf2( m-j+1, jb, a( j, j ), lda, ipiv( j ), iinfo )
*
*           Adjust INFO and the pivot indices.
*
            if( info.eq.0 .and. iinfo.gt.0 )
     $         info = iinfo + j - 1
            do 10 i = j, min( m, j+jb-1 )
               ipiv( i ) = j - 1 + ipiv( i )
   10       continue
*
*           Apply interchanges to columns 1:J-1.
*
            call dlaswp( j-1, a, lda, j, j+jb-1, ipiv, 1 )
*
            if( j+jb.le.n ) then
*
*              Apply interchanges to columns J+JB:N.
*
               call dlaswp( n-j-jb+1, a( 1, j+jb ), lda, j, j+jb-1,
     $                      ipiv, 1 )
*
*              Compute block row of U.
*
               call dtrsm( 'Left', 'Lower', 'No transpose', 'Unit', jb,
     $                     n-j-jb+1, one, a( j, j ), lda, a( j, j+jb ),
     $                     lda )
               if( j+jb.le.m ) then
*
*                 Update trailing submatrix.
*
                  call dgemm( 'No transpose', 'No transpose', m-j-jb+1,
     $                        n-j-jb+1, jb, -one, a( j+jb, j ), lda,
     $                        a( j, j+jb ), lda, one, a( j+jb, j+jb ),
     $                        lda )
               end if
            end if
   20    continue
      end if
      return
*
*     End of DGETRF
*
      end
c-----------------------------------------------------------------------
      subroutine dgetf2( m, n, a, lda, ipiv, info )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      integer            info, lda, m, n
*     ..
*     .. Array Arguments ..
      integer            ipiv( * )
      double precision   a( lda, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      double precision   sfmin 
      integer            i, j, jp
*     ..
*     .. External Functions ..
      double precision   dlamch      
      integer            idamax
      external           dlamch, idamax
*     ..
*     .. External Subroutines ..
      external           dger, dscal, dswap, xerbla
*     ..
*     .. Intrinsic Functions ..
      intrinsic          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      if( m.lt.0 ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, m ) ) then
         info = -4
      end if
      if( info.ne.0 ) then
         call xerbla( 'DGETF2', -info )
         return
      end if
*
*     Quick return if possible
*
      if( m.eq.0 .or. n.eq.0 )
     $   return
*
*     Compute machine safe minimum 
* 
      sfmin = dlamch('S')  
*
      do 10 j = 1, min( m, n )
*
*        Find pivot and test for singularity.
*
         jp = j - 1 + idamax( m-j+1, a( j, j ), 1 )
         ipiv( j ) = jp
         if( a( jp, j ).ne.zero ) then
*
*           Apply the interchange to columns 1:N.
*
            if( jp.ne.j )
     $         call dswap( n, a( j, 1 ), lda, a( jp, 1 ), lda )
*
*           Compute elements J+1:M of J-th column.
*
            if( j.lt.m ) then 
               if( abs(a( j, j )) .ge. sfmin ) then 
                  call dscal( m-j, one / a( j, j ), a( j+1, j ), 1 ) 
               else 
                 do 20 i = 1, m-j 
                    a( j+i, j ) = a( j+i, j ) / a( j, j ) 
   20            continue 
               end if 
            end if 
*
         else if( info.eq.0 ) then
*
            info = j
         end if
*
         if( j.lt.min( m, n ) ) then
*
*           Update trailing submatrix.
*
            call dger( m-j, n-j, -one, a( j+1, j ), 1, a( j, j+1 ), lda,
     $                 a( j+1, j+1 ), lda )
         end if
   10 continue
      return
*
*     End of DGETF2
*
      end
c-----------------------------------------------------------------------
      subroutine dlaswp( n, a, lda, k1, k2, ipiv, incx )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      integer            incx, k1, k2, lda, n
*     ..
*     .. Array Arguments ..
      integer            ipiv( * )
      double precision   a( lda, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
*  Further Details
*  ===============
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
* =====================================================================
*
*     .. Local Scalars ..
      integer            i, i1, i2, inc, ip, ix, ix0, j, k, n32
      double precision   temp
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      if( incx.gt.0 ) then
         ix0 = k1
         i1 = k1
         i2 = k2
         inc = 1
      else if( incx.lt.0 ) then
         ix0 = 1 + ( 1-k2 )*incx
         i1 = k2
         i2 = k1
         inc = -1
      else
         return
      end if
*
      n32 = ( n / 32 )*32
      if( n32.ne.0 ) then
         do 30 j = 1, n32, 32
            ix = ix0
            do 20 i = i1, i2, inc
               ip = ipiv( ix )
               if( ip.ne.i ) then
                  do 10 k = j, j + 31
                     temp = a( i, k )
                     a( i, k ) = a( ip, k )
                     a( ip, k ) = temp
   10             continue
               end if
               ix = ix + incx
   20       continue
   30    continue
      end if
      if( n32.ne.n ) then
         n32 = n32 + 1
         ix = ix0
         do 50 i = i1, i2, inc
            ip = ipiv( ix )
            if( ip.ne.i ) then
               do 40 k = n32, n
                  temp = a( i, k )
                  a( i, k ) = a( ip, k )
                  a( ip, k ) = temp
   40          continue
            end if
            ix = ix + incx
   50    continue
      end if
*
      return
*
*     End of DLASWP
*
      end
c-----------------------------------------------------------------------
      subroutine dgetrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info )
*
*  -- LAPACK routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      character          trans
      integer            info, lda, ldb, n, nrhs
*     ..
*     .. Array Arguments ..
      integer            ipiv( * )
      double precision   a( lda, * ), b( ldb, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A**T * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A**T* X = B  (Transpose)
*          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      double precision   one
      parameter          ( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      logical            notran
*     ..
*     .. External Functions ..
      logical            lsame
      external           lsame
*     ..
*     .. External Subroutines ..
      external           dlaswp, dtrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      intrinsic          max
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      notran = lsame( trans, 'N' )
      if( .not.notran .and. .not.lsame( trans, 'T' ) .and. .not.
     $    lsame( trans, 'C' ) ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( nrhs.lt.0 ) then
         info = -3
      else if( lda.lt.max( 1, n ) ) then
         info = -5
      else if( ldb.lt.max( 1, n ) ) then
         info = -8
      end if
      if( info.ne.0 ) then
         call xerbla( 'DGETRS', -info )
         return
      end if
*
*     Quick return if possible
*
      if( n.eq.0 .or. nrhs.eq.0 )
     $   return
*
      if( notran ) then
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         call dlaswp( nrhs, b, ldb, 1, n, ipiv, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         call dtrsm( 'Left', 'Lower', 'No transpose', 'Unit', n, nrhs,
     $               one, a, lda, b, ldb )
*
*        Solve U*X = B, overwriting B with X.
*
         call dtrsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n,
     $               nrhs, one, a, lda, b, ldb )
      else
*
*        Solve A**T * X = B.
*
*        Solve U**T *X = B, overwriting B with X.
*
         call dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs,
     $               one, a, lda, b, ldb )
*
*        Solve L**T *X = B, overwriting B with X.
*
         call dtrsm( 'Left', 'Lower', 'Transpose', 'Unit', n, nrhs, one,
     $               a, lda, b, ldb )
*
*        Apply row interchanges to the solution vectors.
*
         call dlaswp( nrhs, b, ldb, 1, n, ipiv, -1 )
      end if
*
      return
*
*     End of DGETRS
*
      end
c-----------------------------------------------------------------------
