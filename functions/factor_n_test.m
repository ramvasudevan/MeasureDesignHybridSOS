% check if a monomial can be factorized by a collection of variables (vars). In addition,
% checks if the factor has the provided exponents
% out.flag (1 is vars.^var_exp == cleaned_mss)
% out.msg is present only when there are error
function out=factor_n_test(mss_a,all_variables,vars,vars_exp)
    if length(vars)~=length(vars_exp)
        disp('Error')
        out.msg='Sizes of exponents and varibles do not match';
        return
    end
    
    if ~isfree(vars) || ~isfree(all_variables)
        out.factor=mss_a;
        out.clean=msspoly(1);
        out.flag=0;
        out.msg='Either all_varibles or vars is not a free msspoly';
        return;
    end
    
    % first decide which variables are not important and factorize
    % the monomial
    [a_a,b_a,c_a]=decomp(mss_a);

    unimportant_vars_idx=match(vars,all_variables);

    f_idx= unimportant_vars_idx==0;

    alien_idx=match(all_variables,a_a);
    alien_idxdx= alien_idx==0;
    alien_vars = a_a(alien_idxdx);% variables which are not native to this measure. 
                                  % (Think of them as constants)

    unimportant_vars=[all_variables(f_idx);alien_vars];

    idx_a=match(unimportant_vars,a_a);

    g_idx=find(idx_a~=0);
    h_idx=find(idx_a==0);

    a_clean=a_a(setdiff(1:length(a_a),g_idx));
    b_clean=b_a;
    b_clean(g_idx)=[];

    a_factor=a_a(setdiff(1:length(a_a),h_idx));
    b_factor=b_a;
    b_factor(h_idx)=[];

    % the following happens when (mss_a) cannot be divided by any
    % element of (vars) or when both are (1)
    if isempty(a_clean)
        cleaned_mss=msspoly(1);
        b_clean=0;
    else
        cleaned_mss=recomp(a_clean,b_clean,1);
    end

    % the following happens when (mss_a) does not contain any
    % variable besides those in (vars)
    if isempty(a_factor)
        factor_mss=msspoly(1);
    else
        factor_mss=recomp(a_factor,b_factor,c_a);
    end

    out.factor=factor_mss;
    out.clean=cleaned_mss;
    out.degree.clean=sum(full(b_clean));

    % now to check if the powers of cleaned_mss are as provided
%             keyboard
    out.flag=match_mono_mss(all_variables,cleaned_mss,recomp(vars,vars_exp,1));
end