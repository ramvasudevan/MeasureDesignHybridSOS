function dotx=sim_ex1(t,x,u,xx,T,Uscale)
for i=1:length(u)
u_current(i,1)=Uscale*full(msubs(u{i},xx,[t*T;x]));
end
dotx=T*[-x(1)+(u_current(1)-x(2))*x(1)^2;-x(2)+x(1)^2;u_current(2:end)];
end
