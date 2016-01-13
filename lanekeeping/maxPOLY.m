function [ lambda ] = maxPOLY( x, w, hX, degree )
%
%  x{i}     -- n_i-by-1 free msspoly ( \forall i \in {1,...,nModes} )
%  w{i}     -- n_i-by-1 msspoly in x{i}
%
%  Computes an approximation to:
%  x{i}(4) * x{i}(2) s.t. w(x{i}) >= 1

% define the program
prog = spotsosprog;

% make x{i} indeterminate
prog = prog.withIndeterminate( x );

% create the program variables in each modes
dumb = monomials( x, 0 );
[ prog, dumber, lambda ] = prog.newFreePoly( dumb );

% make the cost function a constraint in the dual program: x_2 * x_4 -
% \lambda >= 0 for all w >= 1
% prog = sosOnK( prog, x( 2 ) * x( 4 ) - lambda, x,  w - 1, degree );
prog = sosOnK( prog, x( 2 ) * x( 4 ) - lambda, x,  [ w; hX ], degree );

% set options
options = spot_sdp_default_options();
options.verbose = 1;

% actually solve the problem
sol = prog.minimize( -lambda, @spot_mosek_sos, options );

lambda = -sol.eval( dumber );
