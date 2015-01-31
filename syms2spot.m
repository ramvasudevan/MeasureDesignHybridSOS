% converts a syms polynomial to spot msspoly
% z - syms polynomial
% x - syms variables
% y - msspoly variables
function total = syms2spot(z,x,y)

numX = numel(x);

[coeff, mon] = coeffs(z, x);

total  = 0;
% degBas = zeros(numel(mon),numX);

for j=1:numel(mon)
    temp = double(coeff(j));
    for k=1:numX
        %degBas(j,k) = double(feval(symengine,'degree', mon(j), x(k)));
        degX = double(feval(symengine,'degree', mon(j), x(k)));
        temp = temp * y(k) ^ degX;
    end
    total = total + temp;
end

end

% [xN,p,M] = decomp(z);
% total = 0;
% i = 1;
% % for i=1:size(M,1)
% for k = 1:size(M,2)
%     temp = M(i,k);
%     for j=1:numel(x)
%         temp = temp * x(j)^p(k,j);
%     end
%     total = total + temp;
% end
