
function y = symmsolv( n, x, iw, rw )
%        y = symmsolv( n, x, iw, rw )
%  SYMMSOLV solves M*y = x for some symmetric
%  positive-definite matrix M.
%  This is a simple example for testing symmlq.
%  rw is used to input shift and pertbn.
%
%  If pertbn = 0, the preconditioner will be exact, so
%  symmlq should require either one or two iterations,
%  depending on whether (A - shift*I) is positive definite or not.
%  If pertbn ~= 0, somewhat more iterations will be required.
                           
shift  = rw(1);
pertbn = rw(2); 
d  = (1.1/n) * (1:n)';

if shift ~= 0
   d  = abs( d - shift );
end %if

if pertbn ~= 0
   for i = 10:10:n
      d(i) = d(i) + pertbn;
   end;
end %if

y  = x./d;

%=================
%End of symmsolv.m
%=================

