
function y = symaprod( n, x, iw, rw )
%  y = symaprod( n, x, iw, rw )
%  symaprod  computes  y = A*x  for some matrix  A.
%  This is a simple example for testing  SYMMLQ.

d  = (1.1/n) * (1:n)';
y  = d.*x;

%=================
%End of symaprod.m
%=================

