%% Sample 1000 Trajectories
warning off;
close all;
clear;
%
% Vehicle parameters
%
m =1370 ; %kg
Iz= 2315.3; % kgm^2
a=1.11; % m
b = 1.59; % m
L=a+b;
%
% Tire parameters
%
Caf= 1.3308e5; % N/rad 
Car= 9.882e4; % N/rad
%
% Nominal longitudinal speed
%
u=100; % kph
u = u/3.6; % m/s
fprintf('longitudinal speed is %f m/s \n',u);
%
params.m=m; params.Iz=Iz; params.a=a; params.b=b;
params.L=L; params.Caf=Caf; params.Car=Car; params.u=u;
% Variables
% 
% y = lateral deviation from lane center in meters
% v = lateral velocity in m/s
% del_psi = psi - psi_desired = delta yaw angle in radians
% r = yaw rate in rad/s
% 
% del_f = Front wheel steering angle in rad
% rd = desired yaw rate
%
A=[0 1 u 0; 
    0 -(Caf+Car)/m/u 0 ((b*Car-a*Caf)/m/u - u); 
    0 0 0 1;
    0 (b*Car-a*Caf)/Iz/u  0 -(a^2 * Caf + b^2 * Car)/Iz/u];
%
B=[0;Caf/m; 0; a*Caf/Iz];
%
E=[0;0;-1;0];
Cy=[1 0 0 0];

c=[1 0 20 0];

model.A=A;
model.B=B;
model.c=c;
model.E=E;
model.Cy=Cy;

Sys=ss(A,B,c,0);
 
%     [K,S,E] = LQR(SYS,Q,R,N) calculates the optimal gain matrix K
%     such that:
%  
%       * For a continuous-time state-space model SYS, the state-feedback
%         law u = -Kx  minimizes the cost function
%  
%               J = Integral {x'Qx + u'Ru + 2*x'Nu} dt
%  
%         subject to the system dynamics  dx/dt = Ax + Bu
%
Kd=.25; Kp=4;
Q=Kp*c'*c + Kd*A'*c'*c*A;
R=600;
[K,S,E]=lqr(Sys,Q,R);
[V,D]=eig(A-B*K);

opts=[];

%% Monte Carlo 
N = 1000;
% K = zeros( 4, 1 );
IC = rand( 4, N );
IC = -0.1 + 0.2 .* IC;
FC = rand( 4, N );
t_end = 1;

% Atilde = [ 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1 ];
% stupid = @( t, x, params, model, K, opts ) Atilde * x;
% [ xdot ] = @stupid( t, x, params, model, K, opts ) Atilde * x;

for i = 1:N
    x0 = IC( :, i )';
%     [t, x] = ode45(stupid, [0, t_end], x0, [], params, model, K, opts);
    [t, x] = ode45(@xdot_LaneKeeping_Barrier, [0, t_end], x0, [], params, model, K, opts);
%     FC( :, i ) = max( abs( x ), [], 1 )';
    FC( :, i ) = x( end, : )';
end





