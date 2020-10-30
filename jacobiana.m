function [jh] = jacobiana(hx,x0)
% Calcula la matriz jacobiana de hx:R^n --> R^m

n = length(x0);
h = feval(hx, x0);
m = length(h);
ep = 1e-05;
jh = zeros(m,n);

for j = 1:n
    y = x0; y(j)=y(j) + ep;
    hy = feval(hx,y);
    jh(:,j) = (hy-h)/ep;
end