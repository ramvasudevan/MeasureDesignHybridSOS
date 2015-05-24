% polynomial to monomial decomposition. Will retain coefficients
function monos=p2d_decomp(poly)
   [a,b,c]=decomp(poly);
   monos=msspoly(ones(size(c,2),1));
   for i=1:size(c,2)
       monos(i,1)=recomp(a,b(i,:),c(1,i));
   end
end