function a_mom = computeApproximateMoments( dim, rel, e_mom, lb, ub )
%  dim      -- positive scalar integer (dimension of vector field)
%  rel      -- positive scalar integer (relaxation order of moment matrix)
%  e_mom    -- (rel + 1)^dim real valued vector (exact empirically computed moments)
%  lb       -- dim real valued vector (lower bound on box constraints in each dimension)
%  ub       -- dim real valued vector (upped bound on box constraints in each dimension)
%
%  Computes an outerapproximation 
%  Computes an outerapproximation to the BRS for:
%
%  xdot{i} = f_i(t,x{i}) + g_i(t,x{i})u,
%   
%  on the time horizon [0,1] and set X, moving trajectories into
%  XT at t = T, subject to u in U, where:
%
%  hU    = { u | each hU(u) >= 0 }
%  X{i}  = { x | each hX{i}(x) >= 0 and sX{i} >= 0  }
%  XT{i} = { x | each hXT{i}(x) >= 0 }
%
%  dlambda(p) = int_X p(x) dx component-wise (lebesgue moments)
%
%  degree dictates the order of certain approximations.
%
%  w{i}(x) >= 1 is an outer approximation of x{i}(0) that could reach XT.
%          /some/ control law.