% add the dual folder
addpath('..');

% parameters used for everything
versionsos    = 'sos';
degree        = 6;
Rx = [ -ones( 4, 1 ) ones( 4, 1 ) ];
Rt            = 0.001;
Ru            = 0.1;
T             = -1;
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
% K = zeros( 1, 4 );
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
x{ 1 }   = msspoly( 'xa', 4 );
xA       = x{ 1 };
df = -[ 0; nu; 0; sigma ] * K * xA( 1:4 );
f{ 1 } = T * ( [ xA( 2 ) + vf * xA( 3 ); ...
               gamma * xA( 2 ) + beta * xA( 4 ); ...
               xA( 4 ) - rd; ...
               omega * xA( 2 ) + epsilon * xA( 4 ) ] ... 
               + df );
% f{ 1 } = T * [ xA( 1 ); xA( 2 ); xA( 3 ); xA( 4 ) ];
g{ 1 }   = zeros( 4, 1 )/Ru;
hX{ 1 }  = [ 1 - xA( 1 )' * xA( 1 ); ...
             1 - xA( 2 )' * xA( 2 ); 
             1 - xA( 3 )' * xA( 3 );
             1 - xA( 4 )' * xA( 4 ) ];
hXT{ 1 } = [ ( 0.1^2 - xA( 1 )' * xA( 1 ) );
             ( 0.1^2 - xA( 2 )' * xA( 2 ) );
             ( 0.1^2 - xA( 3 )' * xA( 3 ) );
             ( 0.1^2 - xA( 4 )' * xA( 4 ) )];
sX{ 1, 1 } = -10;
R{ 1, 1 } = 0 * ones( 4, 1 );
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

[ foo1, foo2, foo3 ] = decomp( w{ 1 } );
foo3( find( abs( foo3 ) < 1e-10 ) ) = 0;
q = recomp( foo1, foo2, foo3 )
lambda = maxPOLY( x{ 1 }, ( q - 1 ), hX{ 1 }, 8 ) 

% plot the measures on the domain
% plotMeasuresonDomain( 'W', w{ 1 }, x{ 1 }, Rx, 'b' );
% Num = 11;
% x1 = linspace( -1, 1, Num );
% x1 = 0;
% x3 = linspace( -1, 1, Num );
% % x3 = 0;
% x2 = linspace(-1,1,100);
% x4 = linspace(-1,1,100);
% X2 = repmat(x2,100,1);
% X4 = repmat(x4,100,1);
% Wrestr = zeros( size( X2( : )' ) );
% 
% for i = 1:Num
%     for j = 1:Num
%         foo = subs( w{ 1 }, [ x{ 1 }( 1 ); x{ 1 }( 3 ) ], [ x1( i ); x3( j ) ] );
%         Wrestr = max( Wrestr, msubs( foo,[x{ 1 }( 2 ); x{ 1 }( 4 )],[X2(:) X4(:)]') );
%     end
% end
% contour(x2,x4,reshape(Wrestr,size(X2)),[1 1],'r');


% save the created files
% save([mfilename '_d' num2str(degree) '_t' num2str(T)]);

