% monomial match only! (at least for now... Can be sped-up by accepting vectors(future))
% 1 if true 0 ow.
% must provide all possible variables as first entry
function f=match_mono_mss(variables,mss_a,mss_b)
    [a_a,b_a,c_a]=decomp(mss_a);
    [a_b,b_b,c_b]=decomp(mss_b);

    % compare variables
    idx1=full(match(a_a,variables));
    idx2=full(match(a_b,variables));

    if sum(abs((idx1==0)-(idx2==0)))
        f=0;
        return;
    end

    % re-arrange matrix elements, if need be
    idx=full(match(a_a,a_b));
    if size(idx,2)==0
        idx=zeros(0,1);
    end
%     a_a_new=a_a(idx);
    b_a_new=b_a(idx);

    % now check if the exponents are the same
    if sum(abs(b_a_new-b_b)) || sum(abs(c_a-c_b))
        f=0;
    else
        f=1;
    end
end