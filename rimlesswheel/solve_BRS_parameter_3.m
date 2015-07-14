function [w,v] = solve_BRS_parameter_3(t,x,theta,f,hX,hXT,hTheta,dl,mu_theta,degree,freeFinalTime)

% number of hybrid modes % modifying the code for the rimless wheel
nModes = 1;

% re-define the time constraint
T = 1;
hT  = t*(T-t);

% define the program
prog = spotsosprog;
prog = prog.withIndeterminate( t );

% create the program variables in each modes
% [ prog, p ] = prog.newFree( 1 );
for i = 1:nModes
    % make x{i} indeterminate
    prog = prog.withIndeterminate( x{ i } );   
    prog = prog.withIndeterminate(theta);
    vmonom{ i } = monomials( [ t; x{ i };theta ], 0:degree );
    [ prog, v{ i }, ~ ] = prog.newFreePoly( vmonom{ i } );
%     keyboard
    % create the w{i} variable
    wmonom{ i } = monomials( x{ i }, 0:degree );
    [ prog, w{ i }, wcoeff{ i } ] = prog.newFreePoly( wmonom{ i } );
    
    % creating the variables that will be used later
    v0{ i } = subs( v{ i }, t, 0 );
    vT{ i } = subs( v{ i }, t, T );
    dvdt{ i } = diff( v{ i }, t );
    dvdx{ i } = diff( v{ i }, [x{ i };theta] );
    Lfv{ i } = dvdt{ i } + dvdx{ i } * f{ i };
end
% keyboard
% creating the constraints and cost function
obj = 0;
for i = 1:nModes
    [prog] = sosOnK( prog, -Lfv{ i }, [ t; x{ i };theta ], [ hT; hX{ i }; hTheta ], degree );
    
    % enforcing that w{i} is SOS on X{i}
    prog = sosOnK( prog, w{ i }, x{ i }, hX{ i }, degree );
    
    % enforcing that w{i} - v{i}(0,*) - 1 is SOS on X{i}
    [~,~,Theta_mass]=decomp(mu_theta(msspoly(1)));
    Theta_mass=full(Theta_mass);
    
    prog = sosOnK( prog, w{ i } - mu_theta(v0{ i })/Theta_mass - 1, [x{ i }], hX{ i }, degree );
    
    % enforcing that v{i}(T,*) is SOS on XT{i}
    if freeFinalTime
        disp( 'Free Final Time' );
        prog = sosOnK( prog, v{ i }, [ t; x{ i }; theta], [ hT; hXT{ i }; hTheta ], degree );
    else
        if ( ~isempty( hXT{ i } ) )
            prog = sosOnK( prog, vT{ i }, [x{ i }; theta], [hXT{ i }; hTheta], degree );
        end
    end
%     keyboard
    % generate the cost function in each mode
    obj = obj + dl{ i }( wmonom{ i } )' * wcoeff{ i };
end


% set options
options = spot_sdp_default_options();
options.verbose = 1;

% actually solve the problem
sol = prog.minimize( obj, @spot_mosek_sos, options );
[flag,idx]=sol.checkSOSsolutions();

if ~flag
    disp('SOS constraint violation!')
end
% putvar(sol)

i=1;
w{ i } = sol.eval( w{ i } );
v{i}=sol.eval(v{i});
end