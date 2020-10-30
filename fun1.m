function [f] = fun1(x)
%función objetivo para prueba de programación cuadrática sucesiva.

f=exp(prod(x)) - (1/2)*( x(1)^3 + x(2)^3 + 1)^2;
end