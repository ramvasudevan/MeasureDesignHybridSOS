function [xdot,delta_f,Fy,Barrier,rd] = xdot_LaneKeeping_Barrier(t,x,params,model,K,opts)
%
% parameters
m=params.m;Iz=params.Iz; a=params.a; b=params.b;
L=params.L; Caf=params.Caf; Car=params.Car; u=params.u;

A=model.A;
B=model.B;
E=model.E;
Cy=model.Cy;


% road curvature   rd * R = u, where R = radius of curvature
%
% centrifigual acceleration = (vel)^2/radius
%
% at 0.2 g, we get radius = (u^2)/(0.2 g)
%
% and we have r dot_theta = u
g=9.81;
ay_max=0.3*g;

radiusMin =  (u^2)/(0.9*ay_max);

if t> 22
    radius=-radiusMin;
    rd = u/radius;
elseif t > 20
    alpha = (t-20)/2;
    rd= alpha * u/(-radiusMin) + (1-alpha)*u/radiusMin;
    
elseif t> 11
    radius=radiusMin;
    rd = u/radius;
elseif t>10
    radius = radiusMin;
    rd = (t-10)*u/radius;
else
    rd=0;
end

BarrierFlag = 0;
switch BarrierFlag
    case 1 % barrier funciton exists
        
        %
        % Barrier Function: (1/2)*(dy/dt)^2 / (max allowed lateral acceleration) +
        % sign(dy/dt)y - 0.9 < 0
        %
        ymax=0.9;
        y=Cy*x;
        dy=Cy*A*x;
        
        h=-((1/2)*(dy)^2/ay_max + signum(dy)*y - ymax);
        
        Barrier=h;
        
        
        xbar=x-[0;0;0;1]*rd;  % do feedforward
        
        v=x(2);
        r=x(4);
        
        %  X = [delta_f;p1] % variable of the QP
        H=diag([1 1e5]);
        f=[0 0];
        %     Aineq=[(dy/ay_max)*(1/(h+h^2))*Cy*A*B, 0;
        %         Caf 0; -Caf 0];
        %     bineq=[-1/(h+h^2)*( (dy/ay_max)*(Cy*A^2*x + Cy*A*E*rd) + signum(dy)*(dy) ) + 1/Barrier;
        %         m*ay_max+Caf*atan((v+a*r)/u)+Car*atan((v-b*r)/u);  m*ay_max-Caf*atan((v+a*r)/u)-Car*atan((v-b*r)/u)];
        Aineq=[(dy/ay_max)*Cy*A*B, 0;
            Caf, 0; -Caf, 0];
        bineq=[ -((dy/ay_max)*(Cy*A^2*x + Cy*A*E*rd -u*rd) + signum(dy)*(dy)) + h ;
            m*ay_max+Caf*((v+a*r)/u)+Car*((v-b*r)/u)+m*u*rd;
            m*ay_max-Caf*((v+a*r)/u)-Car*((v-b*r)/u)-m*u*rd];
        Aeq=[1 1];
        beq=-K*xbar;
        
        [X, Val, ExitFlag] = quadprog(H,f,Aineq,bineq,Aeq,beq,[],[],[],opts);
        if ExitFlag < 1
            fprintf('QP can not generate a legal force for LK!')
            keyboard
        end
        
        delta_f=X(1);
        
    case 0% no barrier funciton exists
        
        xbar=x-[0;0;0;1]*rd;
        
        % LQR feedback 
        delta_f = -K*xbar; 
        
        Barrier = 0;
end


v=x(2);r=x(4);

alpha_f = delta_f - (v+a*r)/u;
alpha_r = 0     - (v-b*r)/u;
%lateral force
Fy = Caf*alpha_f + Car*alpha_r;

%
%
xdot=A*x + B*delta_f + E*rd;
%
end

function y = signum(x)
eps=1e-2;
y=x/( eps+abs(x) );
end


