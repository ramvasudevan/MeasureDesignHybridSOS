function degree=nextEvenMssDegree(expr,vars)
    monos=p2d_decomp(expr);
    b=0;
    for i=1:length(monos)
        out=factor_n_test(monos(i),vars,vars,zeros(1,length(vars)));
        b=max(out.degree.clean,b);
    end
    
    if size(b,1) == 0
        degree=2;
    else
        degree=(1-mod(b,2))*b+(mod(b,2))*(floor(b/2)+1)*2;
    end
end