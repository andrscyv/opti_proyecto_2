function [x, lambda] = pc(Q, A, c, b)
%Metodo directo para resolver el problema
%Min (1/2)*x'*Q*x + c'*x
%S.A. a*x = b
%
%Q es nxn simétrica y positiva definidad en el espacio nulo de A
%A es mxn con rango(A) = m
%c.-vector columna de oren n
%b.-vector columna de orden m.
%Out
%x.-vector columna de orden n con la solución numérica del problema
%lambda.- vector columna de orden que representa el multiplicador
%de Lagrange

% Optimización Numérica
%ITAM
% 25 de agosto de 2020
%-------------------------------------------------
m= length(b); % # de restricciones
n= length(c); % # de variables
K= [Q A'; A zeros(m)]; ld=[-c; b];
%Resolver el sistema lineal k*W=ld
%w= [x; lambda]

w=K\ld;
x=w(1:n);
lambda= w(n+1:n+m);
end