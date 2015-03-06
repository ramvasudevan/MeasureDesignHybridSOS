% add the dual folder
addpath('..');
close all
clear all
% parameters used for everything
degree        = 6;
T             = 1;
freeFinalTime = 0;
Rx            = [-1 1; -1 1];
num_points=50;

prob_options.Uconstant=1;
prob_options.Utimedep=1;
prob_options.num_added=0;
prob_options.T=T;

num_added=prob_options.num_added;

% these are the variables and constraints for everything
t = msspoly( 't', 1 );

% mode 1
x{ 1 }   = msspoly( 'xa', 2);
f{ 1 }   = [x{1}(2); 0];
g{ 1 }   = [ 0; 1 ];
h{ 1 }   = x{1};

% auxiliary variables
hX{ 1 }  = [1.6^2-x{1}'*x{1}];
% hX{1}=[1-x{1};x{1}+1];
% hX{ 1 }  = [-x(1)+Rx(1,2);xA(1)+Rx(1,2);xA(2)+Rx(2,2);-xA(2)+Rx(2,2)];


hXT{ 1 } = -x{1}'*x{1};

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
    [X1,X2]=meshgrid(linspace(Rx(1,1),Rx(1,2),101));
    Wrestr = msubs(w{1},x{1},[X1(:) X2(:) zeros(numel(X1)*num_added,num_added)]');
    [xv,v]=contour(X1(1,:),X1(1,:),reshape(Wrestr,size(X1)),[1 1],'r');
    hold on
end

% True ROA
[X,Y] = meshgrid(Rx(1,1):0.005:Rx(1,2),Rx(2,1):0.005:Rx(2,2));
Z = zeros(size(X));
for i = 1:size(X,1)
    for j = 1:size(X,2)
        Z(i,j) = doubleIntegValueFun([X(i,j),Y(i,j)]);
    end
end
contour(X,Y,Z,[1 1], '-r', 'linewidth',2);
xlabel('x_1') ; ylabel('x_2');
legend('Outer', 'True')

% return
hand=gcf;

rp=randperm(size(xv,2));
for i=1:num_points
    y=ode45(@(time,x_state) simulate_DI(time,x_state,u,[t;x{1}]),0:1e-1:1,[xv(:,rp(i));zeros(num_added,num_added)]);
    hold on
    if num_added==1
        plot3(y.y(1,:),y.y(2,:),y.y(3,:),'b');
    else
        plot(y.y(1,:),y.y(2,:),'b');
    end
end

