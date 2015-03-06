function [w,v,x,u] = DI_SOS(t,x,f,g,h,hX,hXT,dl,degree,freeFinalTime,prob_options)
%
%  t        -- 1-by-1 msspoly.
%  x{i}     -- n_i-by-1 free msspoly ( \forall i \in {1,...,nModes} )
%  f{i}     -- n_i-by-1 msspoly in x{i}
%  g{i}     -- n_i-by-1 msspoly in x{i}
%  h{i}     -- m_i msspoly in x{i} (output function)
%  dl{i}    -- function mapping n_i-by-1 msspoly's into n-by-1 double.
%  degree   -- positive scalar integer
%  prob_options -> Uconstant, Utimedep, num_added, Uscale, Xscale, T
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
%  u{m} is a feedback control law. In the current iteration, it is assumed
%  that the outputs are partial state representations

if nargin < 9
    freeFinalTime = 0; 
    prob_options=struct;
end

if ~isfield(prob_options,'Uconstant')
    prob_options.Uconstant=1;
end
if ~isfield(prob_options,'Utimedep')
    prob_options.Utimedep=1;
end
if ~isfield(prob_options,'num_added')
    prob_options.num_added=0; 
end
if ~isfield(prob_options,'T')
    prob_options.T=1;
end
if ~isfield(prob_options,'Uscale')
    prob_options.Uscale=1;
end
if ~isfield(prob_options,'Xscale')
    prob_options.Xscale=ones(length(x{1}),1);
end


T=prob_options.T;
num_added=prob_options.num_added;

%---
% x and T -scaling
Xscale=prob_options.Xscale;
XTscale=[T;Xscale];
if ~isdouble(f{1})
    [k_a,k_b,k_c]=decomp(f{1});
    k_o=match([t;x{1}],k_a);
    multip=power(repmat(XTscale(k_o)',size(k_b,1),1),k_b);
    k_cnew=bsxfun(@times,k_c,prod(multip,2)');
    f{1}=recomp(k_a,k_b,k_cnew);
    f{1}=inv(diag(Xscale))*f{1}; % final scaling by T is performed below
end

if ~isdouble(g{1})
    [k_a,k_b,k_c]=decomp(g{1}(:));
    k_o=match([t;x{1}],k_a);
    multip=power(repmat(XTscale(k_o)',size(k_b,1),1),k_b);
    k_cnew=bsxfun(@times,k_c,prod(multip,2)');
    g{1}=reshape(recomp(k_a,k_b,k_cnew),size(g{1},1),size(g{1},2));
    g{1}=inv(diag(Xscale))*g{1}; % final scaling by T is performed below
end

%---

% adding auxiliary variables
z{1} = msspoly('za',prob_options.num_added);
x{1} = [x{1};z{1}];

% adding the dynamics of auxiliary states
f{1}=T*[f{1};zeros(num_added,1)];
g{1}=prob_options.Uscale*T*[[g{1};zeros(num_added,1)], [zeros(length(x{1})-length(z{1}),num_added);diag(ones(num_added,1))]];

% adding the auxiliary states to measured output
h{1}=[h{1};z{1}];

% overwriting the box constraint
dl{ 1 }    = boxMoments( x{ 1 }, -ones(length(x{1}),1), ones(length(x{1}),1) );

% amending terminal constraint
if num_added
    hXT{ 1 } = -x{1}'*x{1};
end

% amending the state constraints to impose box constaint on z
if num_added
    hX{1}=hX{1}+z{1}'*z{1};
end

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
    
    % create the v{i} variable
    % (HACK) force the degree of v to be the smallest multiple of two such
    % that L_fv is of degree 2k
    
%     [k_a,k_b,k_c]=decomp(f{1});
    keyboard
%     k_d=match(k_a,[t;x{1}]);
%     k_e=find(k_d~=0);
%     max_k=max(sum(k_b(:,k_d(k_e)),2));
%     
%     v_degree=floor((degree-max_k+1)/2)*2;
%     
    v_degree=degree;
    
    vmonom{ i } = monomials( [ t; x{ i } ], 0:v_degree );
    [ prog, v{ i }, ~ ] = prog.newFreePoly( vmonom{ i } );
    
    % create the w{i} variable
    wmonom{ i } = monomials( x{ i }, 0:degree );
    [ prog, w{ i }, wcoeff{ i } ] = prog.newFreePoly( wmonom{ i } );
    
    % create the q{ i } variable which describes the g(t,x) * u(t,x) term
    q{ i } = msspoly( zeros( size( g{ i }, 2 ) ,1) );
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
    % enforcing that -( Lv{i} + 1' * q{ i } ) is SOS on domain and that q{i}>0
    e = ones( size( g{ i }, 2 ), 1 );
%     keyboard
    [prog] = sosOnK( prog, -( Lfv{ i } + e' * q{ i } ), [ t; x{ i } ], [ hT; hX{ i } ], degree );
%     keyboard
    for j = 1:size( g{ i }, 2 )
        prog = sosOnK( prog, q{ i }( j ) - Lgv{ i }( j ), [ t; x{ i } ], [ hT; hX{ i } ], degree );
        prog = sosOnK( prog, q{ i }( j ) + Lgv{ i }( j ), [ t; x{ i } ], [ hT; hX{ i } ], degree );
        prog = sosOnK( prog, q{ i } , [ t; x{ i } ], [ hT; hX{ i } ], degree );
        prog=prog.updateSigmaTokens(); % this must be placed after enforcing positivity of q
    end
    
    % enforcing that w{i} is SOS on X{i}
    prog = sosOnK( prog, w{ i }, x{ i }, hX{ i }, degree );
    
    % enforcing that w{i} - v{i}(0,*) - 1 is SOS on X{i}
    prog = sosOnK( prog, w{ i } - v0{ i } - 1, x{ i }, hX{ i }, degree );
    
    % enforcing that v{i}(T,*) is SOS on XT{i}
    if freeFinalTime
        disp( 'Free Final Time' );
        prog = sosOnK( prog, v{ i }, [ t; x{ i } ], [ hT; hXT{ i } ], degree );
    else
        if ( ~isempty( hXT{ i } ) )
            prog = sosOnK( prog, vT{ i }, x{ i }, hXT{ i }, degree );
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
[flag,idx]=sol.checkSOSsolutions();

if ~flag
    disp('SOS constraint violation!')
end
putvar(sol)
i=1;
% construct w{i}
    w{ i } = sol.eval( w{ i } );
    
% constuct v{i}
    v{i}=sol.eval(v{i});
% construct u{i}

M_mu=sol.sosDual(sol.prog.sosTokens(1));
mono_mu=sol.prog.gramMonomials{sol.prog.sosTokens(1)};

% keyboard

if prob_options.Utimedep
    non_zero_coeff_indices=find(msubs(mono_mu,[t;x{1}],+(full(match([t;h{1}],[t;x{1}]))>0))==1);
else
    non_zero_coeff_indices=find(msubs(mono_mu,[t;x{1}],+(full(match([h{1}],[t;x{1}]))>0))==1);
end

% to remove any constant term in the control
if ~prob_options.Uconstant
    non_zero_coeff_indices=non_zero_coeff_indices(find(non_zero_coeff_indices>1));
end

[~,k_b1,~]=decomp(mono_mu);

% keyboard
for i=1:size(g{1},2)

    % -- extracting only rows and columns that match with the monomials of u
    % going to assume that the msspoly free is the same in both the monomials
    % of \sigma_plus and \mu

    [~,k_b2,~]=decomp(sol.prog.gramMonomials{sol.prog.sosPositiveTokens(i)});
    [~,idx]=ismember(k_b2,k_b1,'rows');
    non_zero_coeff_indices_reduced=intersect(non_zero_coeff_indices,idx);
    % --

    non_zero_coeff_indices_reduced=sort(non_zero_coeff_indices_reduced);

    M_mu_reduced=M_mu(idx,non_zero_coeff_indices_reduced);
    
    M_sp=sol.sosDual(sol.prog.sosPositiveTokens(i));
    M_sn=sol.sosDual(sol.prog.sosNegativeTokens(i));

    mom_sigma_pos=M_sp(1,:)';
    mom_sigma_neg=M_sn(1,:)';
    u_coeffs=pinv(M_mu_reduced)*(mom_sigma_pos-mom_sigma_neg);

    ma_coeff=max(abs(u_coeffs));
    u_coeffs =diag((abs(u_coeffs)>1e-4*ma_coeff))*u_coeffs; % removing smaller coefficients to remove insignificant monomials
                                                            % this is fine
                                                            % since all
                                                            % variables are in
                                                            % [-1,1]

    u{i}=u_coeffs'*mono_mu(non_zero_coeff_indices_reduced);
end
end









