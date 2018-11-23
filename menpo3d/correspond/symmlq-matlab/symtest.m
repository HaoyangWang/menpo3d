
function symtest( n, precon, shift, pertbn )
%        symtest( n, precon, shift, pertbn )
% SYMTEST sets up and solves a system (A - shift * I)x = b, using
% SYMAPROD to define A and SYMMSOLV to define a preconditioner.

true   = 1;              false  = 0;
show   = true;           check  = true;
iw     = 0;
rw     = [ shift; pertbn ];

disp( ' ' )
disp( '==============' )
disp( 'Test of SYMMLQ' )
disp( '==============' )
disp( sprintf( 'shift = %23.14e  pertbn = %23.14e', shift, pertbn ) )

%  Set the true solution and the rhs
%  so that  (A - shift*I) * xtrue = b.

xtrue = (n : -1 : 1)';
b = symaprod ( n, xtrue, iw, rw );
b = (- shift) * xtrue +  b;
                
%  Set other parameters and solve.

itnlim = n * 2;
rtol   = 1.0e-12;

[ x, istop, itn, anorm, acond, rnorm, xnorm ] = ...
           symmlq( n, b, 'symaprod', 'symmsolv', iw, rw, ...
                   precon, shift, show, check, itnlim, rtol );
                       
%  print the solution and some clue about whether it is ok.

disp(' '); disp( 'Solution:' )
for j = 1:n
   disp( sprintf( '%22.14e', x(j) ) )
end
disp(' ')

w        = x - xtrue;
enorm    = norm( w ) / norm( xtrue );
etol     = 1.0e-5;

if enorm < etol
   disp( 'SYMMLQ appears to have been successful' )
else
   disp( 'SYMMLQ appears to have failed' )
end
disp( sprintf( 'Relative error in x = %10.1e', enorm ) )
%  End of test.

