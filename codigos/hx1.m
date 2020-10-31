function [h] = hx1(x)

h(1) = x(1)^2 + x(2)^2 + x(3)^2-1;
h(2) = x(1)*x(2)*x(3)-1;

end