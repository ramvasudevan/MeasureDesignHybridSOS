function w = measureDesignHybridSOSReduceFootprint(t,x,f,g,hX,hXT,sX,R,dl,degree,freeFinalTime)
%
%  t        -- 1-by-1 msspoly.
%  x{i}     -- n_i-by-1 free msspoly ( \forall i \in {1,...,nModes} )
%  f{i}     -- n_i-by-1 msspoly in x{i}
%  g{i}     -- n_i-by-1 msspoly in x{i}
%  hX{i}    -- m_i-by-1 msspoly in x{i}
%  hXT{i}   -- k_i-by-1 msspoly in x{i}
%  sX{j,i}  -- l_i-by-1 msspoly in x{j,i} (defines the guard from i to j)
%  R{j,i}   -- nModes-by-1 msspoly in x{j,i} (defines the reset from i to j)
%  dl{i}    -- function mapping n_i-by-1 msspoly's into n-by-1 double.
%  degree   -- positive scalar integer
%
%  Computes an outerapproximation to the BRS for:
%
%  xdot{i} = f_i(t,x{i}) + g_i(t,x{i})u,
%
%  on the time horizon [0,1] and set X, moving trajectories into
%  XT at t = 1, subject to u in [-1, 1], where:
%
%  X{i}  = { x | each hX{i}(x) >= 0 and sX{i} >= 0  }
%  XT{i} = { x | each hXT{i}(x) >= 0 }
%
%  dlambda(p) = int_X p(x) dx component-wise (lebesgue moments)
%
%  degree dictates the order of certain approximations.
%
%  w{i}(x) >= 1 is an outer approximation of x{i}(0) that could reach XT.
%          /some/ control law.

if nargin < 9, freeFinalTime = 0; end

% number of hybrid modes
nModes = length( x );

% define the time constraint
T = 1;
hT  = t*(T-t);

% define the program
prog = spotsosprog;
prog = prog.withIndeterminate( t );

% create the program variables in each modes
[ prog, p ] = prog.newFree( 1 );
for i = 1:nModes
    % make x{i} indeterminate
    prog = prog.withIndeterminate( x{ i } );
    
    % create the v{i} variable
    vmonom{ i } = monomials( [ t; x{ i } ], 0:degree );
    [ prog, v{ i }, ~ ] = prog.newFreePoly( vmonom{ i } );
    
    % create the w{i} variable
    wmonom{ i } = monomials( x{ i }, 0:degree );
    [ prog, w{ i }, wcoeff{ i } ] = prog.newFreePoly( wmonom{ i } );
    
    % create the q{ i } variable which describes the g(t,x) * u(t,x) term
    q{ i } = msspoly( zeros( size( g{ i }, 2 ) ) );
    for j = 1:size( g{ i }, 2 )
        qmonom{ i, j } = monomials( [ t; x{ i } ], 0:degree );
        [ prog, q{ i }( j ), ~ ] = prog.newFreePoly( qmonom{ i, j } );
    end
    
    % creating the variables that will be used later
    v0{ i } = subs( v{ i }, t, 0 );
    vT{ i } = subs( v{ i }, t, T );
    dvdt{ i } = diff( v{ i }, t );
    dvdx{ i } = diff( v{ i }, x{ i } );
    Lfv{ i } = dvdt{ i } + dvdx{ i } * f{ i };
    Lgv{ i } = dvdx{ i } * g{ i };
    
end

% creating the constraints and cost function
obj = 0;
for i = 1:nModes
    % enforcing that -( Lv{i} + 1' * q{ i } ) is SOS on domain
    e = ones( size( g{ i }, 2 ), 1 );
    prog = sosOnK( prog, -( Lfv{ i } + e' * q{ i } ), [ t; x{ i } ], [ hT; hX{ i } ], degree );
    for j = 1:size( g{ i }, 2 )
        prog = sosOnK( prog, q{ i }( j ) - Lgv{ i }( j ), [ t; x{ i } ], [ hT; hX{ i } ], degree );
        prog = sosOnK( prog, q{ i }( j ) + Lgv{ i }( j ), [ t; x{ i } ], [ hT; hX{ i } ], degree );
    end
    
    % enforcing that w{i} - v{i}(0,*) - p - 1 is SOS on X{i}
    prog = sosOnK( prog, w{ i } - v0{ i } - p - 1, x{ i }, hX{ i }, degree );
    
    % enforcing that w{i} is SOS on X{i}
    prog = sosOnK( prog, w{ i }, x{ i }, hX{ i }, degree );
    
    % enforcing that v{i}(T,*) + p is SOS on XT{i}
    if freeFinalTime
        disp( 'Free Final Time' );
        prog = sosOnK( prog, v{ i } + p, [ t; x{ i } ], [ hT; hXT{ i } ], degree );
    else
        if ( ~isempty( hXT{ i } ) )
            prog = sosOnK( prog, vT{ i } + p, x{ i }, hXT{ i }, degree );
        end
    end
    
    % enforcing that v{i} - v{j} \circ (1,R_{j,i}) is SOS on sX{i}
    for j = 1:nModes
        if ( ~isempty( sX{ j, i } ) ) % if its empty there isn't a guard between these
            %         if ~isscalar( sX{ j, i }) && ~isempty( sX{ j, i }) % if its -1 there isn't a guard between these two modes
            %         if( sX{ j, i } ~= -1 )
            vj_helper = subs( v{ j }, x{ j }, R{ j , i } );
            prog = sosOnK( prog, v{ i } - vj_helper, [ t; x{ i } ] , [ hT; sX{ j, i } ], degree );
        end
    end
    
    % generate the cost function in each mode
    obj = obj + dl{ i }( wmonom{ i } )' * wcoeff{ i };
end

% set options
options = spot_sdp_default_options();
options.verbose = 1;

% actually solve the problem
sol = prog.minimize( obj, @spot_mosek_sos, options );

% construct w{i}
for i = 1:nModes
    w{ i } = sol.eval( w{ i } );
end









