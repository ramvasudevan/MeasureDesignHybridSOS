function p = polyn2spot(z,x)

if numel(x) < numel(z.VarNames)
    disp('ERROR: incorrect number x variables');
    return;
end

p = 0;
for i = 1:numel(z.Coefficients)
    m = z.Coefficients(i);
    for j = 1:size(z.ModelTerms,2)
        m = m * x(j) ^ z.ModelTerms(i,j);
    end
    p = p + m;
end

end