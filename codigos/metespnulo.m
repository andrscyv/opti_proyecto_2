function[x]=metespnulo(Q,A,c,b)
%Método del espacio nulo para el problema
%
%Min (1/2)*x'*Q*x + c'*x
%S.A.: A*x=b
%A es mxn con rango m (m<=n)
%c vector en R^n y b vector en R^m, todos los vectores son columnas
%Out
%x.- aproximación al mínimo del problema
%
%Optimización Numérica
%ITAM
% 27 de agosto de 2020
%--------------------------------------------------------------------

Z=null(A);
%A*Ain = I_m
% Min (1/2)*x'*x s.a. A*x=b
xp=A\b;
%cambio de variable: x= Z*y+ xp

%definimos matriz y lado derecho
B= Z'*Q*Z;
d=Z'*(Q*xp+c);

%Sistemas lineal: B*y=d

y=B\(-d);

x= Z*y+xp;%Duda

end