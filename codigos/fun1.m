function [f] = fun1(x)
%funci�n objetivo para prueba de programaci�n cuadr�tica sucesiva.

f=exp(prod(x)) - (1/2)*( x(1)^3 + x(2)^3 + 1)^2;
end