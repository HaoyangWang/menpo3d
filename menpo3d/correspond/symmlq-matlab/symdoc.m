%  symdoc.m is documentation for symmlq.m.
%
%        [ x, istop, itn, anorm, acond, rnorm, xnorm ] = ...
%          symmlq( n, b, aprodname, msolvename, iw, rw,...
%                  precon, shift, show, check, itnlim, rtol )
%
%  SYMMLQ is designed to solve the system of linear equations Ax = b
%  where A is an n by n symmetric matrix and b is a given vector.
%  The matrix A is not required to be positive definite.
%  (If A is known to be definite, the method of conjugate gradients
%  might be preferred, since it will require about the same number of
%  iterations as SYMMLQ but slightly less work per iteration.)
%  A is intended to be large and sparse.  It is accessed
%  by means of a function call of the form
%
%                  y = aprod( n, x, iw, rw )
%
%  which must return the product y = Ax for any given vector x.
%
%
%  More generally, SYMMLQ is designed to solve the system
%
%                  (A - shift*I) x = b
%
%  where shift is a specified scalar value.  If shift and b
%  are suitably chosen, the computed vector x may approximate an
%  (unnormalized) eigenvector of A, as in the methods of
%  inverse iteration and/or Rayleigh-quotient iteration.
%  Again, the matrix (A - shift*I) need not be positive definite.
%  The work per iteration is very slightly less if shift = 0.
%
%  A further option is that of preconditioning, which may reduce
%  the number of iterations required.  If M = C C' is a positive
%  definite matrix that is known to approximate  (A - shift*I)
%  in some sense, and if systems of the form  My = x  can be
%  solved efficiently, the parameters precon and msolve may be
%  used (see below).  When  precon = true, SYMMLQ will
%  implicitly solve the system of equations
%
%           P(A - shift*I)P' xbar  =  P b,
%
%  i.e.                Abar xbar  =  bbar
%  where                       P  =  C^(-1),
%                           Abar  =  P(A - shift*I)P',
%                           bbar  =  P b,
%
%  and return the solution     x  =  P' xbar.
%
%  In the discussion below, eps refers to the machine precision.
%
%  Parameters
%  ----------
%  n          input     n, the dimension of the matrix A.
%
%  b(n,1)     input     The rhs vector b.
%
%  aprodname  string    Name of a matlab function 'aprod'
%                       defining the matrix A.
%                       For a given vector x, the statement
%
%                             y = aprod( n, x, iw, rw )
%
%                       must return the product y = Ax
%                       without altering the vector x.
%
%  msolvename string    Name of an optional matlab function 'msolve'
%                       defining a preconditioning matrix M,
%                       which should approximate (A - shift*I) in some sense.
%                       M must be positive definite.
%                       For a given vector x, the statement
%
%                             y = msolve( n, x, iw, rw )
%
%                       must solve the linear system My = x
%                       without altering the vector x.
%
%                       In general, M should be chosen so that Abar
%                       has clustered eigenvalues.  For example,
%                       if A is positive definite, Abar would ideally
%                       be close to a multiple of I.
%                       If A or A - shift*I is indefinite, Abar might
%                       be close to a multiple of diag( I  -I ).
%
%  WARNING:   The files containing the functions 'aprod' and 'msolve'
%             must not be called aprodname.m or msolvename.m !!!!
%
%  iw, rw     workspace These are not used by symmlq but are passed
%                       "as is" to aprod and msolve.  They are intended
%                       for global information that aprod and msolve
%                       may need.
%
%  precon     input     If precon = true, preconditioning will
%                       be invoked.  Otherwise, function msolve
%                       will not be referenced.
%
%  shift      input     Should be zero if the system Ax = b is to be
%                       solved.  Otherwise, it could be an
%                       approximation to an eigenvalue of A, such as
%                       the Rayleigh quotient b'Ab / (b'b)
%                       corresponding to the vector b.
%                       If b is sufficiently like an eigenvector
%                       corresponding to an eigenvalue near shift,
%                       then the computed x may have very large
%                       components.  When normalized, x may be
%                       closer to an eigenvector than b.
%
%  show       input     If show = true, a summary of the iterations
%                       will be printed.
%
%  check      input     If check = true, an extra call of aprod will
%                       be used to check if A is symmetric.  Also,
%                       if precon = true, an extra call of msolve
%                       will be used to check if M is symmetric.
%
%  itnlim     input     An upper limit on the number of iterations.
%
%  rtol       input     A user-specified tolerance. SYMMLQ terminates
%                       if it appears that norm(rbar) is smaller than
%                           rtol * norm(Abar) * norm(xbar),  where
%                       rbar is the transformed residual vector,
%                           rbar = bbar - Abar xbar.
%
%                       If shift = 0 and precon = false, SYMMLQ
%                       terminates if norm(b - A*x) is smaller than
%                           rtol * norm(A) * norm(x).
%
%  x(n,1)     output    Returns the computed solution x.
%
%  istop      output    An integer giving the reason for termination.
%
%             -1        The matrix Abar appears to be a multiple of
%                       the identity. The solution is a multiple of b.
%
%              0        b = 0, so the exact solution is x = 0.
%                       No iterations were performed.
%
%              1        Norm(rbar) appears to be less than
%                       the value  rtol * norm(Abar) * norm(xbar).
%                       The solution in  x  should be acceptable.
%
%              2        Norm(rbar) appears to be less than
%                       the value eps * norm(Abar) * norm(xbar).
%                       This means that the residual is as small as
%                       seems reasonable on this machine.
%
%              3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
%                       which should indicate that x has essentially
%                       converged to an eigenvector of A
%                       corresponding to the eigenvalue shift.
%
%              4        acond (see below) has exceeded 0.1/eps, so
%                       the matrix Abar must be very ill-conditioned.
%                       x may not contain an acceptable solution.
%
%              5        The iteration limit was reached before any of
%                       the previous criteria were satisfied.
%
%              6        The matrix defined by aprod does not appear
%                       to be symmetric.
%                       For certain vectors y = Av and r = Ay, the
%                       products y'y and r'v differ significantly.
%
%              7        The matrix defined by msolve does not appear
%                       to be symmetric.
%                       For vectors satisfying My = v and Mr = y, the
%                       products y'y and r'v differ significantly.
%
%              8        An inner product of the form x'(M\x) was
%                       not positive, so the preconditioned matrix
%                       M does not appear to be positive definite.
%
%                       If istop >= 5, the final  x  may not be an
%                       acceptable solution.
%
%  itn        output    The number of iterations performed.
%
%  anorm      output    Estimate of the norm of the matrix operator
%                       Abar = P(A - shift*I)P'.
%
%  acond      output    An estimate of the condition of Abar above.
%                       This will usually be a substantial
%                       under-estimate of the true condition.
%
%  rnorm      output    The norm of the final residual vector,
%                       b  -  (A - shift*I) x.
%
%  xnorm      output    The norm of the final solution vector x.
%
%
%  Internal workspace arrays:
%
%  r1(n,1)   r2(n,1)   v(n,1)    w(n,1)    y(n,1)
%  Before exit, r1 contains the final residual vector.
%
%
%  This routine is an implementation of the algorithm described in
%  the following references:
%
%  C.C. Paige and M.A. Saunders,  Solution of Sparse Indefinite
%       Systems of Linear Equations,
%       SIAM J. Numer. Anal. 12, 4, September 1975, pp. 617-629.
%
%  J.G. Lewis,  Algorithms for Sparse Matrix Eigenvalue Problems,
%       Report STAN-CS-77-595, Computer Science Department,
%       Stanford University, Stanford, California, March 1977.
%
%  Applications of SYMMLQ and the theory of preconditioning
%  are described in the following references:
%
%  D.B. Szyld and O.B. Widlund,  Applications of Conjugate Gradient
%       Type Methods to Eigenvalue Calculations,
%       in R. Vichnevetsky and R.S. Steplman (editors),
%       Advances in Computer Methods for Partial Differential
%       Equations -- III, IMACS, 1979, 167-173.
%
%  D.B. Szyld,  A Two-level Iterative Method for Large Sparse
%       Generalized Eigenvalue Calculations,
%       Ph. D. dissertation, Department of Mathematics,
%       New York University, New York, October 1983.
%  ------------------------------------------------------------------
%
%  SYMMLQ development:
%         1972: First version.
%         1975: John Lewis recommended modifications to help with
%               inverse iteration:
%               1. Reorthogonalize v1 and v2.
%               2. Regard the solution as x = x1  +  bstep * b,
%                  with x1 and bstep accumulated separately
%                  and bstep % b added at the end.
%                  (In inverse iteration, b might be close to the
%                  required x already, so x1 may be a lot smaller
%                  than the multiple of b.)
%         1978: Daniel Szyld and Olaf Widlund implemented the first
%               form of preconditioning.
%               This required both a solve and a multiply with M.
%         1979: Implemented present method for preconditioning.
%               This requires only a solve with M.
%         1984: Sven Hammarling noted corrections to tnorm and x1lq.
%               SYMMLQ added to NAG library.
%  15 Sep 1985: Final F66 version.  SYMMLQ sent to "misc" in netlib.
%  16 Feb 1989: First F77 version.
%
%  22 Feb 1989: Hans Mittelmann observed beta2 = 0 (hence failure)
%               if Abar = const*I.  istop = -1 added for this case.
%
%  01 Mar 1989: Hans Mittelmann observed premature termination on
%               ( 1  1  1 )     (   )                   ( 1  1    )
%               ( 1  1    ) x = ( 1 ),  for which  T3 = ( 1  1  1 ).
%               ( 1     1 )     (   )                   (    1  1 )
%               T2 is exactly singular, so estimating cond(A) from
%               the diagonals of Lbar is unsafe.  We now use
%               L       or  Lbar         depending on whether
%               LQNORM  or  CGNORM       is least.
%
%  03 Mar 1989: Eps computed internally instead of coming in as a
%               parameter.
%  07 Jun 1989: Check added as a parameter to say if A and M
%               should be checked for symmetry.
%               This is the date of the Fortran version referred to next.
%  15 May 1990: MATLAB m-file symmlq.m derived from Fortran version
%               by Gerald N. Miranda Jr, UCSD.
%  16 May 1990: aprodname and msolvename added as parameters
%               to allow aprod and msolve to be called via eval.
%               iw and rw added as workspace parameters.
%
%  Michael A. Saunders            (na.saunders @ NA-net.stanford.edu)
%  Department of Operations Research
%  Stanford University
%  Stanford, CA 94305-4022
%  ------------------------------------------------------------------
