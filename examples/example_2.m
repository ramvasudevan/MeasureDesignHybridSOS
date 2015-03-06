% Dubin's car

% add the dual folder
addpath('..');
close all
clear all
% parameters used for everything
degree        = 10;
T             = 4;
freeFinalTime = 0;
Rx            = [-1 1; -1 1;-1 1];
num_points=5;

prob_options.Uconstant=1;
prob_options.Utimedep=1;
prob_options.num_added=0;
prob_options.T=T;
prob_options.Uscale=1;
prob_options.Xscale=[2;2;2];
Uscale=prob_options.Uscale;

num_added=prob_options.num_added;

% these are the variables and constraints for everything
t = msspoly( 't', 1 );

% mode 1
x{ 1 }   = msspoly( 'xa', 3);
f{ 1 }   = [0;0;0];
g{ 1 }   = [1 0;0 1;-x{1}(2) x{1}(1)];
h{ 1 }   = x{1}(3);

% auxiliary variables
hX{ 1 }  = [4-x{1}'*x{1}];
% hX{1}=[1-x{1};x{1}+1];
% hX{ 1 }  = [-x(1)+Rx(1,2);xA(1)+Rx(1,2);xA(2)+Rx(2,2);-xA(2)+Rx(2,2)];


hXT{ 1 } = .01-x{1}'*x{1};

dl{ 1 }    = boxMoments( x{ 1 }, Rx( :, 1 ), Rx( :, 2 ) );

% run the code
[w,v,x,u] = DI_SOS(t,x,f,g,h,hX,hXT,dl,degree,freeFinalTime,prob_options);

%% 
    [X1,X2,X3]=meshgrid(linspace(1*Rx(1,1),1*Rx(1,2),30));
    Wrestr = msubs(w{1},x{1},[X1(:) X2(:) X3(:)]');
    Vrestr = msubs(v{1},[t;x{1}],[zeros(size(X1(:))) X1(:) X2(:) X3(:)]');
    figure(1)
    isosurface(X1,X2,X3,reshape(full(Wrestr),size(X1)),1);
    figure(2)
    isosurface(X1,X2,X3,reshape(full(Vrestr),size(X1)),0);
    hold on
return
%%
figure
hold on
x1=X1(:);
x2=X2(:);
x3=X3(:);
for i=1:length(Vrestr)
    if (Vrestr(i)>=0)
        plot3(x1(i),x2(i),x3(i),'.b')
    end
end
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
%%

rp=randperm(size(xv,2));
for i=1:min(num_points,length(rp))
    y=ode45(@(time,x_state) sim_BI(time,x_state,u,[t;x{1}],T,Uscale),0:1e-1:1,[xv(:,rp(i));zeros(num_added,num_added)]);
    hold on
    plot3(y.y(1,:),y.y(2,:),y.y(3,:),'b');
end

