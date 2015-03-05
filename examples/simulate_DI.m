function dotx=simulate_DI(t,x,u,xx)
for i=1:length(u)
u_current(i,1)=full(msubs(u{i},xx,[t;x]));
end
dotx=[x(2);u_current];
end
