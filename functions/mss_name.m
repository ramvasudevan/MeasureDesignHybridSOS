% works only on monomials and will strip any coefficients
function [nm] = mss_name(p)
    [a,b,~]=decomp(p);
    nm='';
    for i=1:length(a)
        nm=strcat(nm,strcat(name(a(i)), '^',num2str(b(i))));
    end
    if strcmp(nm,'')
        nm='*const*';
    else
        nm=nm{1};
    end
end