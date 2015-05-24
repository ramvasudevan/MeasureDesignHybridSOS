function y=mss_v2s_vert(x,e)
% with nargin==1, re-arrange the rows of x into block symmetric matrix
% example 3 dimensional matrix
%     [ x(1) x(2) x(3)    ]
% y=  [ x(2) x(4) x(5)    ]
%     [ x(3) x(5) x(6)    ]

% SM 3.03.15

x=x(:);
n=size(x,1);
m=(sqrt(1+8*n)-1)/2;
if m~=round(m), error('unable to convert: wrong dimension'); y=[];return; end

y=zeros(m);
count=1;
for i=1:m
    for j=i:m
        y(j,i)=x(count);
        count=count+1;
    end
end

y = y + tril(y,-1)';
if nargin >1
   [u,s,v]=svd(y);
   s=s*(s>=e);
   y=u*s*v';
end