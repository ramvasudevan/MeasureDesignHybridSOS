% add the dual folder
addpath('..');
addpath('../functions')
addpath(genpath('C:\Users\shankar\Documents\MATLAB\spotless'))
clear all
% parameters used for everything
versionsos    = 'sos';
degree        = 14;
T             = 1;
freeFinalTime = 0;


% these are the variables and constraints for everything
t = msspoly( 't', 1 );

% mode 1
x{ 1 }   = msspoly( 'xa',1);
xA       = x{1};
Rx   = repmat([-1 1],length(xA),2);
theta = msspoly('p',1);


% 1D example
example_no = 2;

switch example_no
    case 1
        f{ 1 } = xA + 1 * theta;
    case 2
        f{ 1 } = -xA + 0.2 * theta;
    case 3
        f{ 1 } = -0.7 * xA + 0.2 * theta - 0.1;
end
        
f{1} = [f{1};0];

% 2D example
if length(xA) == 2
    xAb=1*xA;
    thee=1-theta*0.3;
    f{1} =(thee*xAb-(1-thee)*xAb^3+xAb^5/120-xAb^7);
    f{1} = [(thee^3*xAb(1)-(1-thee)*xAb(2)^3+xAb(1)^5/120-xAb(2)^7);...
        (.1-(1-thee)^2*xAb(2)+(thee)*xAb(1)^2+xAb(1)^6/12-xAb(2)^6);...
        0];
end
hX{1}=[1-xA;1+xA];

hXT{1}=.1^2-(xA-.3*ones(length(xA),1))'*(xA-.3*ones(length(xA),1));
dl{ 1 }    = boxMoments( x{ 1 }, Rx( :, 1 ), Rx( :, 2 ) );
hTheta = [theta;1-theta];
mu_theta = boxMoments(theta,0,1);

% run the code
[w,v] = solve_BRS_parameter_3(t,x,theta,f,hX,hXT,hTheta,dl,mu_theta,degree,freeFinalTime);
%%
hold on

x1=linspace(-1,1,301);
Wrestr = msubs(w{1},x{1},x1);
plot(x1,Wrestr,'r',[-1,1],[1,1],'k--',[-1,1],[0,0],'k--')
ylim([-.2,1.5])
legend('w')
disp('uncertain case done. Press a key to continue.')

% running base_case
pause(1)



[w_n0,v_n0] = solve_BRS_normal(t,x, {msubs(f{1}(1:end-1), theta, 0)} ,hX,hXT,dl,degree,freeFinalTime);

x1=linspace(-1,1,301);
Wrestr = msubs(w_n0{1},x{1},x1);
plot(x1,Wrestr,'m')
ylim([-.2,1.5])


[w_n1,v_n1] = solve_BRS_normal(t,x, {msubs(f{1}(1:end-1), theta, 1)}, hX,hXT,dl,degree,freeFinalTime);

x1=linspace(-1,1,301);
Wrestr = msubs(w_n1{1},x{1},x1);
plot(x1,Wrestr,'b')
ylim([-.2,1.5])


if length(xA)==1
    return
end
%%
% clf
hold on
x1 = linspace(-1,1,100);
x2 = linspace(-1,1,100);
X1 = repmat(x1,100,1);
X2 = repmat(x2,100,1)';
Wrestr = msubs(w{1},x{1},[X1(:) X2(:)]');
[xv,~]=contour(x1,x2,reshape(Wrestr,size(X1)),[1 1],'r-');
return