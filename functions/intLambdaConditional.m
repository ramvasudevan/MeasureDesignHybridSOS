% lambda is the lebesgue measure on unov
function out=intLambdaConditional(expr,all_variables,unov,lambda)
    mono_break=p2d_decomp(expr);
    out=msspoly(0);
    for i=1:length(mono_break)
%         keyboard
        split=factor_n_test(mono_break(i),all_variables,unov,zeros(1,length(unov)));
        out=out+split.clean.coeff*split.factor*lambda(split.clean);
    end
end