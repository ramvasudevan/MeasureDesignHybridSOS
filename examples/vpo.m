% Van der pol oscillator

% add the dual folder
addpath('..');
close all
clear all
% parameters used for everything
degree        = 16;
T             = 100;
freeFinalTime = 1;
Rx            = [-1 1];
num_points=5;

prob_options.Uconstant=1;
prob_options.Utimedep=1;
prob_options.num_added=0;
prob_options.T=T;
prob_options.Uscale=1;

Uscale=prob_options.Uscale;

num_added=prob_options.num_added;

% these are the variables and constraints for everything
t = msspoly( 't', 1 );
x{ 1 }   = msspoly( 'xa', 2);
Rx=repmat(Rx,size(x{1},1),1);

prob_options.Xscale=[1.2;1.2];

aa=1;
f{ 1 }   =  (1 / aa) * [-2 * ( x{1}(2) * aa );
    0.8 * ( x{1}(1) * aa) + 10 * ( ( x{1}(1) * aa)^2 - .21 ) * (x{1}(2) * aa)];

g{ 1 }   = [0;0];
h{ 1 }   = x{1};

% auxiliary variables
hX{ 1 }  = [5-x{1}'*x{1};1-x{1};1+x{1}];
% hX{1}=[1-x{1};x{1}+1];
% hX{ 1 }  = [-x(1)+Rx(1,2);xA(1)+Rx(1,2);xA(2)+Rx(2,2);-xA(2)+Rx(2,2)];


hXT{ 1 } = .01^2 - ( x{1}'*x{1} ) / (aa ^ 2);

dl{ 1 }    = boxMoments( x{ 1 }, Rx( :, 1 ), Rx( :, 2 ) );

% run the code
[w,v,x,u] = DI_SOS(t,x,f,g,h,hX,hXT,dl,degree,freeFinalTime,prob_options);

%% 
    x1=linspace(-1,1,200);
    [X1,X2]=meshgrid(x1);
    Wrestr = msubs(w{1},x{1},[X1(:) X2(:)]');
    Vrestr = msubs(v{1},[t;x{1}],[zeros(size(X1(:))) X1(:) X2(:)]');
    surf(x1,x1,reshape(Wrestr,size(X1)));