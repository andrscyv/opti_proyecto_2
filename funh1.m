function [h] = funh1(x)
%Restricciones para probar programación cuadrática sucesiva 

h(1) = x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(5)^2 -10;
h(2) = x(2)*x(3) - 5*x(4)*x(5);
h(3) = x(1)^3 + x(2)^3 +1 ;

h=h';
end