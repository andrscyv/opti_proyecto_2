function [x, lambda, kiter] = pcspowell(fx, hx, x0)
% M�todo de Programaci�n Cuadr�tica Sucesiva con actualizaci�n de Powell
% para el problema 
% Min fx
% s.a. hx = 0
%
% fx: R^n --> R , hx: R^n-->R^m, ambas son dos veces continuamente
% diferenciables.
% In 
% fx.- cadena de caracteres con la funci�n objetivo en Matlab.
% hx.- cadena de caracteres con la funci�n de restricciones.
% x0.- vector inicial de dimensi�n n.
% Out
% x.- aproximaci�n al m�nimo local del problema
% lambda.- multiplicador de Lagrange de la restricci�n: hx = 0.
% kiter.- n�mero de iteraciones:
%--------------------------------------------------------------------------
% Actualizaci�n de Powell para la matriz del Lagrangiano
%
% Optimizaci�n Num�rica
% ITAM
% 6 de octubre de 2020
%--------------------------------------------------------------------------
% Funciones que se necesitan:
% gradiente.m .- calcula el gradiente de fx en un punto en R^n
% jacobiana.m .- calcua la matriz jacobiana de hx en cualquier punto de R^n
% , matriz jacobiana es mxn, donde m es el n�mero de restricciones 

tol = 1e-05;    % tolerancia para la CNPO
maxiter = 25;   % m�ximo n�mero de iteraciones
kiter = 0;      % contador de iteraciones

n = length(x0);      % dimensi�n de la variable
h = feval(hx,x0);
m = length(h);      % n�mero de restricciones

lambda = zeros(m,1);    % multiplicador de Lagrange original
B = eye(n);

% CNPO
gf = gradiente(fx,x0);
jh = jacobiana(hx,x0);
vcnpo = [gf + jh'*lambda; h];
x=x0;

% iteraci�n 
while(norm(vcnpo)>tol && kiter<maxiter)
    % Resolver el subproblema cuadr�tico de igualdad
    K = [B jh'; jh zeros(m)];
    ld = -[gf;h];
    d = K\ld;
    px = d(1:n);
    % Actualizamos 
    xn = x + px;                % Nuevo punto
    lambda = d(n+1:n+m);        % Nuevo multiplicador de Lagrange
    
    % Calcular la nueva matriz B
    s = xn-x;
    gfn = gradiente(fx,xn);
    jhn = jacobiana(hx,xn);
    lag1 = gf + jh'*lambda;
    lag2 = gfn + jhn'*lambda;
    y = lag2 - lag1;
    
    aux = s'*B*s;
    
    if (s'*y >= (0.2)*aux)
        r = y;
    else
        theta = ((0.2)*aux - s'*y)/(aux-s'*y); % actualizar el vector r
        r = theta*(B*s) + (1-theta)*y;
    end
    
    B = B - ((B*s*s'*B)/aux) + (r*r')/(s'*r);
    
    %Si B es mal condicionada
    if(rcond(B)<1e-04)
        B = norm(vcnpo)*eye(n);
    end
    
    % Condiciones necesarias de primer orden
    x = xn;
    gf = gradiente(fx,x); jh = jacobiana(hx, x);
    vcnpo = [ gf + jh'*lambda; h];
    kiter= kiter+1;
    fprintf('%2.0f %2.8f %2.8f\n',kiter, norm(vcnpo),norm(px));
end
end
