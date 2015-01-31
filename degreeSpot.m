% equivalent to Yalmip's degree function
%  gets the degree of p w.r.t. x variables
function degs = degreeSpot(p,x)

degs = zeros(x.dim);
pNum = p.var;
xNum = x.var;
if ~isempty(pNum)
    degrees = p.pow;
    for i = 1:numel(xNum)
        if all((pNum == xNum(i))==0) == 1
            degs(i) = 0;
        else
            if numel(degrees) == 1
                if pNum == xNum(i)
                    degs(i) = degrees;
                else
                    degs(i) = 0;
                end
            else
                degs(i) = max(degrees(pNum == xNum(i)));
            end
        end
    end
end

degs = degs';
end