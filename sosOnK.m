% [ prgout, tk ] = sosOnK( prg, p, x, h, d );
%
% prg -- spotsosprg.
% p   -- 1-by-1 msspoly in x.
% x   -- n-by-1 free msspoly
% h   -- m-by-1 msspoly in x. (semialgebraic constraints defining K)
% d   -- scalar integer, d > deg(g).
%
% prgout --  new program with constraints
%       s(i) SOS, p - s'*h SOS
%       s(i) of maximal degree s.t. deg(s'*h) <= d.
% tk -- token associated with p - s'*h SOS.
%
function [ prg, tk ] = sosOnK( prg, p, x, h, d )
    m = size( h, 1 ); % total number of h's
    S = msspoly( zeros( m, 1 ) ); % create this many multiplier variables
    for i = 1:m
        
        % choosing the multipler polynomials to be the right degree
        [Vexp,Pexp,~]=decomp(h(i));
        tmp1=match(Vexp,x);
        f=find(tmp1~=0);
        Ptmp=Pexp(:,tmp1(f));
        d_tmp=floor((d-full(max(sum(Ptmp,2))))/2)*2;
        
        % uncomment if you want to use the provided degree
%         d_tmp=d;
        
        bases = monomials( x, 0:d_tmp);
        [ prg, S( i ) ] = prg.newFreePoly( bases ); % make a polynomial whose coefficients are free variables
        prg = prg.withSOS( S( i ) ); % ensure that the polynomial is SOS
    end
%     keyboard;
    [ prg, tk ] = prg.withSOS( p - S'*h ); % ensure that p is SOS on K  % Punitar's Positivstellansatz
%     prg.sosVariable=length(prg.)
    prg.sosTokens=[prg.sosTokens;tk]; % storing the tokens corresponding to the desires SOS constraint
end