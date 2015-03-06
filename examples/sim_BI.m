function dotx=sim_BI(t,x,u,xx,T,Uscale)
for i=1:length(u)
u_current(i,1)=Uscale*full(msubs(u{i},xx,[t*T;x]));
end
dotx=T*[u_current(1);
    u_current(2);
    x(1)*u_current(2)-x(2)*u_current(1);u_current(3:end)];
end
