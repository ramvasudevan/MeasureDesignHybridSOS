% converts a spot polynomial to syms
% z - spot polynomial
% x - msspoly variables
% y - syms variables
function total = spot2syms(z,x,y)

addpath('../base');

[xN,p,M] = decomp_consistent(z,x);

total = 0;
i = 1;
% for i=1:size(M,1)
for k = 1:size(M,2)
    temp = M(i,k);
    for j=1:numel(y)
        temp = temp * y(j)^p(k,j);
    end
    total = total + temp;
end

end