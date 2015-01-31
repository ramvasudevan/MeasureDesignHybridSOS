% add the dual folder
addpath('..');

% parameters used for everything
versionsos    = 'sos';
degree        = 10;
Rx            = [-1 1; -1 1];
Rt            = 0.0001;
Ru            = 0;
T             = 10;
freeFinalTime = 0;

% these are the variables and constraints for everything
t = msspoly( 't', 1 );
u = msspoly( 'u', 1 );
hU = Ru - u'*u;

% these are parameters which are model specific
LCcoeff = [-0.3153   -0.1905    0.9751    0.0055    0.5027];
l = 1;
alpha = 0.4;
gamma = 0.2;
a = [1; 1];
guardX  = alpha + gamma;

% mode 1
x{ 1 }   = msspoly( 'xa', 2 );
xA       = x{1};
f{ 1 }   = [T * (a(2)*xA(2)); ...
    T * 1 / l * ((a(1)*xA(1)) - (a(1)*xA(1))^3/6)];
g{ 1 }   = [0;0];
hX{ 1 }  = [1 - xA(1)' * xA(1); ...
    1 - xA(2)' * xA(2)];
hXT{ 1 } = [Rt - (xA(2) - LCcoeff*flipud(monomials(xA(1),0:numel(LCcoeff)-1))); ...
    (xA(2) - LCcoeff*flipud(monomials(xA(1),0:numel(LCcoeff)-1))) + Rt];
sX{ 1, 1 } = [0 - (xA(1)-guardX)'*(xA(1)-guardX); ...
    1 - xA(2)'*xA(2)];
R{ 1, 1 }  = [2 * gamma - (a(1)*xA(1)); (a(2)*xA(2)) * (1 - (2*alpha)^2 / 2)];
dl{ 1 }    = boxMoments( x{ 1 }, Rx( :, 1 ), Rx( :, 2 ) );

% run the code
switch versionsos
    case 'sdsos'
       w  = measureDesignHybridSDSOS(t,x,u,f,g,hX,hXT,sX,R,hU,dl,degree,freeFinalTime);
    case 'dsos'
       w = measureDesignHybridDSOS(t,x,u,f,g,hX,hXT,sX,R,hU,dl,degree,freeFinalTime);
    case 'sos'
       w = measureDesignHybridSOS(t,x,u,f,g,hX,hXT,sX,R,hU,dl,degree,freeFinalTime);
end

% plot the measures on the domain
plotMeasuresonDomain( 'W', w{ 1 }, x{ 1 }, Rx, 'b' );

% save the created files
% save([mfilename '_d' num2str(degree) '_t' num2str(T)]);

