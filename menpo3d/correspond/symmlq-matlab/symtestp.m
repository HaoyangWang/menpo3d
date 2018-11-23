
%======================================
%     File symtestp.m   (test program)
%     Script file for testing symmlq.m.
%======================================
                                                    
      true   = 1;      false  = 0;
      normal = false;  precon = true;
      shift  = 0.25;   pertbn = 0.1;

%     Test the unlikely tiny cases that often trip us up.

      symtest( 1, normal, 0 , 0 );
      symtest( 2, normal, 0 , 0 );
      symtest( 1, precon, 0 , 0 );
      symtest( 2, precon, 0 , 0 );

%     Test small positive-definite and indefinite systems with
%     exact preconditioners.  SYMMLQ should take 1 and 2 iterations
%     respectively.

      n      = 5;
      symtest( n, precon,     0, 0 );
      symtest( n, precon, shift, 0 );

%     Test other combinations on larger systems.
%     With no preconditioning, SYMMLQ will take about N iterations.

      n      = 50
      symtest( n, normal,     0, 0 );
      symtest( n, normal, shift, 0 );

%     PERTBN makes the preconditioners incorrect in N/10 entries.
%     SYMMLQ should take about N/10 iterations.

      symtest( n, precon,     0, pertbn );
      symtest( n, precon, shift, pertbn );

%     End of Main script for testing SYMMLQ

