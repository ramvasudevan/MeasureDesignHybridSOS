% add the dual folder
addpath('..');

% parameters used for everything
versionsos    = 'sos';
degree        = 6;
Rx = [ -ones( 5, 1 ) ones( 5, 1 ) ];
Rt            = 0.001;
Ru            = 0.1;
T             = 1;
freeFinalTime = 0;

% dynamic specific parameters
m =1370 ; %kg
Iz= 2315.3; % kgm^2
a=1.11; % m
b = 1.59; % m
L=a+b;
Caf= 1.3308e5; % N/rad 
Car= 9.882e4; % N/rad
vf = 100/3.6; % m/s
ay_max = 0.3 * 9.8; 
K = [ 0.0816    0.0208    2.2891    0.3812 ]; % FROM LQR CODE XIANGRU GAVE
gamma = -( Caf + Car )/( m * vf );
beta = ( b * Car - a * Caf )/( m * vf ) - vf;
omega = ( b * Car - a * Caf )/( Iz * vf );
epsilon = -( a^2 * Caf + b^2 * Car )/( Iz * vf );
nu = Caf/m;
sigma = a*Caf/Iz;
% rd = vf^2/( 0.9 * ay_max );
rd = 0;

% these are the variables and constraints for everything
t = msspoly( 't', 1 );

% mode 1
x{ 1 }   = msspoly( 'xa', 5 );
xA       = x{ 1 };
df = K * xA( 1:4 );
f{ 1 } = T * [ xA( 2 ) + vf * xA( 3 ); ...
               gamma * xA( 2 ) + beta * xA( 4 ) + nu * df; ...
               xA( 4 ) - rd; ...
               omega * xA( 2 ) + epsilon * xA( 4 ) + sigma * df; ...
               omega * xA( 2 )^2 + epsilon * xA( 2 ) * xA( 4 ) + sigma * xA( 2 ) * df ...
               + gamma * xA( 2 ) * xA( 4 ) + beta * xA( 4 )^2 + nu * xA( 4 ) * df ];
g{ 1 }   = zeros( 5, 1 )/Ru;
hX{ 1 }  = [ 1 - xA( 1 )' * xA( 1 ); ...
             1 - xA( 2 )' * xA( 2 ); 
             1 - xA( 3 )' * xA( 3 );
             1 - xA( 4 )' * xA( 4 );
             1 - xA( 5 )' * xA( 5 ); ];
hXT{ 1 } = [ ( 0.1^2 - xA( 1 )' * xA( 1 ) );
             ( 0.1^2 - xA( 2 )' * xA( 2 ) );
             ( 0.1^2 - xA( 3 )' * xA( 3 ) );
             ( 0.1^2 - xA( 4 )' * xA( 4 ) );
             ( 0.1^2 - xA( 5 )' * xA( 5 ) ); ];
sX{ 1, 1 } = 1;
R{ 1, 1 } = -100 * ones( 5, 1 );
dl{ 1 }    = boxMoments( x{ 1 }, Rx( :, 1 ), Rx( :, 2 ) );

% run the code
switch versionsos
    case 'sdsos'
       w  = measureDesignHybridSDSOSReduceFootprint(t,x,f,g,hX,hXT,sX,R,dl,degree,freeFinalTime);
    case 'dsos'
       w = measureDesignHybridDSOSReduceFootprint(t,x,f,g,hX,hXT,sX,R,dl,degree,freeFinalTime);
    case 'sos'
       [ w, v ] = measureDesignHybridSOSReduceFootprint(t,x,f,g,hX,hXT,sX,R,dl,degree,freeFinalTime);
end

% plot the measures on the domain
% plotMeasuresonDomain( 'W', w{ 1 }, x{ 1 }, Rx, 'b' );
% x1 = linspace(-1,1,100);
% x2 = linspace(-1,1,100);
% X1 = repmat(x1,110,1);
% X2 = repmat(x2,100,1)';
% Wrestr = msubs(w{1},x{1},[X1(:) X2(:)]');
% contour(x1,x2,reshape(Wrestr,size(X1)),[1 1],'r');


% save the created files
% save([mfilename '_d' num2str(degree) '_t' num2str(T)]);

