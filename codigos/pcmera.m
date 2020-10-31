function [x, lambda] = pcmera(Q, A, c, b)
%Metodo directo para resolver el problema
%Min (1/2)*x'*Q*x + c'*x
%S.A. a*x = b
%
%Q es nxn sim�trica positiva definidad en R^n y cuya inversa es f�cil
%    de calcular.
%A matriz mxn de rango m.
%c.-vector columna de oren n.
%b.-vector columna de orden m.
%Out
%x.-vector columna de orden n con la soluci�n num�rica del problema
%lambda.- vector columna de orden que representa el multiplicador
%de Lagrange

% Optimizaci�n Num�rica
%ITAM
% 25 de agosto de 2020
%----------------------------------------
S = inv(Q);
B = A*S*A';
ld = -( b + A*S*c);

%Usamos el gradiente conjugado
lambda = B \ ld;

x = -S *( c + A'*lambda);

end