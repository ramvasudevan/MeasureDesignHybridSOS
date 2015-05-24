classdef mommeas
    properties 
        degree=[]; % this is the relaxation for each conditional
        mom_matrix=zeros(0,0);
        mom_monomials=msspoly(0);
        conditionals_indx=containers.Map('KeyType','char','ValueType','double');
        moments=containers.Map('KeyType','char','ValueType','double');
        variables=msspoly(0);
        ov=msspoly(0);
        unov=msspoly(0);
        conditionals=[];
    end
    
    % make some of the variables private
    properties (Access=private)
        
    end
    methods (Access=private)
        
        % check if the provided inputs are consistent
        function meas=checkConsistency(meas)
        end
        
        % imposes the same order in which variables appear in monomials to
        % ensure that multiple keys do not appear and that the hashtable is
        % complete
        function mono=imposeOrder(meas,mono)
        end
        
        % monomial match only! (at least for now... Can be sped-up by accepting vectors(future))
        % 1 if true 0 ow
        function f=match_mono_mss(meas,mss_a,mss_b)
            [a_a,b_a,c_a]=decomp(mss_a);
            [a_b,b_b,c_b]=decomp(mss_b);
            
            % compare variables
            idx1=full(match(a_a,meas.variables));
            idx2=full(match(a_b,meas.variables));
            
            if sum(abs((idx1==0)-(idx2==0)))
                f=0;
                return;
            end
            
            % re-arrange matrix elements, if need be
            idx=full(match(a_a,a_b));
            if size(idx,2)==0
                idx=zeros(0,1);
            end
            a_a_new=a_a(idx);
            b_a_new=b_a(idx);

            % now check if the exponents are the same
            if sum(abs(b_a_new-b_b)) || sum(abs(c_a-c_b))
                f=0;
            else
                f=1;
            end
        end
        
        
        % extract moments until degree
        % inefficient code!! 
        % could benefit from using unique to wean repeated entries
        function meas=extractMoments(meas,monos)
            mom_mss_matrix=monos*monos';
            mom_mss_entries=mss_s2v(mom_mss_matrix);
            mom_entries=mss_s2v(meas.mom_matrix);
            
            mom_mss_moments=monomials([meas.ov;meas.unov],0:meas.degree);
            collect_moments=zeros(size(mom_mss_moments));
            
            for i=1:length(mom_mss_moments)
                j=1;
                flag=1;
                    while flag
                        if meas.match_mono_mss(mom_mss_moments(i),mom_mss_entries(j))
                            flag=0;
                            collect_moments(i,1)=mom_entries(j,1);
                        end
                        j=j+1;
                    end
            end
%             meas.moments=collect_moments;
%             meas.mom_monomials=mom_mss_moments;

            for i=1:length(mom_mss_moments)
%                 keyboard
                mono_name=mss_name(mom_mss_moments(i));
                meas.moments(mono_name)=collect_moments(i);
            end
        end
        
        
        function meas=extractConditionals(meas)
            % limiting max ov monomial to degree/2 to ensure
            % that 
            monos_ov=monomials(meas.ov,0:(meas.degree/2));
            monos_unov=monomials(meas.unov,0:(meas.degree));
            
            for i=1:length(monos_unov)
                indx_key=mss_name(monos_unov(i));
                Phi=zeros(0,length(monos_ov));
                Y=zeros(0,1);
%                 keyboard
                for j=1:length(monos_ov)
                    reg=monos_ov*monos_ov(j);
                    
                    % going to ensure that we can hit monos_ov with
                    % element monos_ov(j)
                    y_key=mss_name(monos_unov(i)*monos_ov(j));
                    for l=1:length(reg)
                        key=mss_name(reg(l));
                        if ~meas.moments.isKey(key) || ~meas.moments.isKey(y_key)
                            if size(Phi,1)==j
                                Phi(j,:)=[];
                            end
                            if size(Phi,1)==j
                                Y(j,1)=[];
                            end
                            break
                        else
                            Phi(j,l)=meas.moments(key);
                            Y(j,1)=meas.moments(y_key);
                        end
                    end
                end
                meas.conditionals=[meas.conditionals;monos_ov'*(pinv(Phi)*Y)];
                meas.conditionals_indx(indx_key)=length(meas.conditionals);
            end
           
            % sanity check examples
            meas.doChecks();
        end
        
    end
    
    methods
        % ov: observed states, unov: unobserved states 
        function meas=mommeas(degree_relaxtion,M,monos,ov,unov)
            meas.mom_matrix=M;
            meas.mom_monomials=monos;
            meas.degree=degree_relaxtion;
            meas.variables=[ov;unov];
            meas.ov=ov;
            meas.unov=unov;
            meas.checkConsistency();

            meas=meas.extractMoments(monos);
            meas=meas.extractConditionals();
        end
        
        % function to integrate-out conditional
        function out=intConditional(meas,expr)
            mono_break=meas.p2d_decomp(expr);
            out=msspoly(0);
%             keyboard
            for i=1:length(mono_break)
                split=meas.factor_n_test(mono_break(i),meas.unov,zeros(1,length(meas.unov)));
                out=out+split.clean.coeff*split.factor*meas.conditionals(meas.conditionals_indx(mss_name(split.clean)));
            end
        end
        
        % check if a monomial can be factorized by a collection of variables (vars). In addition,
        % checks if the factor has the provided exponents
        % out.flag (1 is vars.^var_exp == cleaned_mss)
        % out.msg is present only when there are error
        function out=factor_n_test(meas,mss_a,vars,vars_exp)
            if length(vars)~=length(vars_exp)
                disp('Error')
                out.msg='Sizes of exponents and varibles do not match';
                return
            end
            
            % first decide which variables are not important and factorize
            % the monomial
            [a_a,b_a,c_a]=decomp(mss_a);
            
            unimportant_vars_idx=match(vars,meas.variables);
            
            f_idx= unimportant_vars_idx==0;
            
            alien_idx=match(meas.variables,a_a);
            alien_idxdx= alien_idx==0;
            alien_vars = a_a(alien_idxdx);% variables which are not native to this measure. 
                                          % (Think of them as constants)
                       
            unimportant_vars=[meas.variables(f_idx);alien_vars];
            
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
            
            % now to check if the powers of cleaned_mss are as provided
%             keyboard
            out.flag=meas.match_mono_mss(cleaned_mss,recomp(vars,vars_exp,1));
        end
        
        % feed it a function and this will integrate using the moments
        % inputs in sequence : expression (polynomial), mode [optional :
        % 0 (initialization) /1 (decomposition)]
        function out=integrate(meas,varargin)
            if isempty(varargin)
                disp('Provide an expression to integrate over this measure')
                return;
            elseif length(varargin)==1
                expr=varargin{1};
                mode=0;
            else 
                expr=varargin{1};
                mode=+(varargin{2}~=0);
            end
            
            out=0;
            if mode==0
                mono_list=meas.p2d_decomp(expr);
%                 keyboard
                for i=1:length(mono_list)
                    tmp=mono_list(i);
                    out=out+tmp.coeff*meas.moments(mss_name(tmp));  
                end
            else
%                 keyboard
                e_stage=meas.intConditional(expr);
                mono_list=meas.p2d_decomp(e_stage);
%                 keyboard
                for i=1:length(mono_list)
                    tmp=mono_list(i);
%                     keyboard
                    out=out+tmp.coeff*meas.moments(mss_name(tmp));  
                end
            end
%             keyboard
        end
        
        % polynomial to monomial decomposition. Will retain coefficients
        function monos=p2d_decomp(meas,poly)
           [a,b,c]=decomp(poly);
           monos=msspoly(ones(size(c,2),1));
           for i=1:size(c,2)
               monos(i,1)=recomp(a,b(i,:),c(1,i));
           end
        end
        
        % sanity checks for conditional decomposition
        %%% SOMETHING BROKEN HERE
        function flag=doChecks(meas)
            flag=1;
            idx=[];
            num_tests=10;
            monos=monomials(meas.variables,0:(meas.degree)/2);
            test_coeffs=rand(num_tests,length(monos));
            test_funs=test_coeffs*monos;
%             keyboard
            for i=1:num_tests
                given=meas.integrate(test_funs(1),0);
                built=meas.integrate(test_funs(1),1);
                if (given-built)/given*100>=1e-2
                    flag=0;
                    idx=[idx;i];
                end
            end
            if ~flag
                disp('Failed some tests')
                keyboard
            end
        end
        
        % deconstructor
        function delete(meas)
        end
        
    end
end