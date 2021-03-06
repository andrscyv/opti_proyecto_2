function [h] = hesfera(x)
% Función de restricciones del problema de np puntos en la esfera unitaria
% de dimensión tres.
%
%  Optimzación Numérica
% ITAM
% 20 de octubre de 2020

n = length(x);
np = floor(n/3);
h = zeros(np,1);
for j = 1:np
    uj = x(3*(j-1)+1:3*j);
    h(j) = uj'*uj-1;
end

h(np+1) = x(1)-1;
h(np+2) = x(2);
h(np+3) = x(3);

end