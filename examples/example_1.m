% example from automatica -- " A unifying point of view on output
% feedback designs for global asymptotic stabilization", Andrieu and Praly,
% 2009

% add the dual folder
addpath('..');
close all
clear all
% parameters used for everything
degree        = 6;
T             = 10;
freeFinalTime = 0;
Rx            = [-1 1; -1 1];
num_points=5;

prob_options.Uconstant=1;
prob_options.Utimedep=1;
prob_options.num_added=0;
prob_options.T=T;
prob_options.Uscale=1000;
Uscale=prob_options.Uscale;

num_added=prob_options.num_added;

% these are the variables and constraints for everything
t = msspoly( 't', 1 );

% mode 1
x{ 1 }   = msspoly( 'xa', 2);
f{ 1 }   = [-x{1}(1)-x{1}(1)^2*x{1}(2); -x{1}(2)+x{1}(1)^2];
g{ 1 }   = [ x{1}(1)^2; 0 ];
h{ 1 }   = x{1}(1);

% auxiliary variables
hX{ 1 }  = [1.6^2-x{1}'*x{1}];
% hX{1}=[1-x{1};x{1}+1];
% hX{ 1 }  = [-x(1)+Rx(1,2);xA(1)+Rx(1,2);xA(2)+Rx(2,2);-xA(2)+Rx(2,2)];


hXT{ 1 } = .01-x{1}'*x{1};

dl{ 1 }    = boxMoments( x{ 1 }, Rx( :, 1 ), Rx( :, 2 ) );

% run the code
[w,v,x,u] = DI_SOS(t,x,f,g,h,hX,hXT,dl,degree,freeFinalTime,prob_options);

%% 
if num_added==1
    
    [X1,X2]=meshgrid(linspace(Rx(1,1),Rx(1,2),101));
    Wrestr = msubs(w{1},x{1},[X1(:) X2(:) zeros(numel(X1)*num_added,num_added)]');
    [xv,v]=contour(X1(1,:),X1(1,:),reshape(Wrestr,size(X1)),[1 1],'r');
        
    [X1,X2,X3]=meshgrid(linspace(Rx(1,1),Rx(1,2),20));
    Wrestr = msubs(w{1},x{1},[X1(:) X2(:) X3(:)]');
    isosurface(X1,X2,X3,reshape(full(Wrestr),size(X1)),1);
    hold on
else
    x1=linspace(Rx(1,1),Rx(1,2),101);
    [X1,X2]=meshgrid(x1);
    Wrestr = msubs(w{1},x{1},[X1(:) X2(:) zeros(numel(X1)*num_added,num_added)]');
    figure(1);
    surf(linspace(Rx(1,1),Rx(1,2),101),x1,reshape(Wrestr,size(X1)));
    figure(2);
    [xv,v]=contour(X1(1,:),X1(1,:),reshape(Wrestr,size(X1)),[1 1],'r');
    hold on
end

rp=randperm(size(xv,2));
for i=1:min(num_points,length(rp))
    y=ode45(@(time,x_state) sim_ex1(time,x_state,u,[t;x{1}],T,Uscale),0:1e-1:1,[xv(:,rp(i));zeros(num_added,num_added)]);
    hold on
    if num_added==1
        plot3(y.y(1,:),y.y(2,:),y.y(3,:),'b');
    else
        plot(y.y(1,:),y.y(2,:),'b');
    end
end

