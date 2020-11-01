function [x, lambda, k] = pcsglobal(fx, hx, x0)
% Método de programacion Cuadratica Sucesiva con busqueda de linea,
% usando la funcion de merito L-1 y actualizacion de la hessiana
% con la formula BFGS para el problema
% Min fx
% Sujeto a hx = 0
%
% fx y hx son cadenas de caracteres con las funciones en Matlab
% de la funcion objetivo y las restricciones del problema
% El vector x0 es el valor inicial
% Salida
% x.- aproximacion al minimo local
% lambda.- multiplicador de Lagrange asociado a x.
% k.- numero de iteraciones realizadas.
%
% Debe usar las funciones: gradiente.m y jacobiana.m para calcular
% las primeras derivadas.
%
% Javier Montiel Gonzalez 159216, Alexis Ayala Redon 156916
% Andres Cruz y Vera 155899

% Par�metros inciales
tol = 1e-5;     % Tolerancia de las cnpo
maxk = 100;     % M�ximo n�mero de iteraciones
c1 = 1e-2;      % Valor de la constante para la b�squeda de l�nea
C = 1;          % Par�metro para encontrar una derivada direccional negativa

h = feval(hx, x0); % El valor de las restricciones en el punto inicial
n = length(x0); %Dimension de la variable
m = length(h); % Dimension de las restricciones


lambda = zeros(m,1); %Inicializamos el multiplicador de Lagrange
B = eye(n); % Inicializamos la aproximacion a la Hessiana

% Condiciones necesarias de primer orden
gf = gradiente(fx,x0);
jh = jacobiana(hx, x0);

v = [ gf + jh'*lambda; h];
x = x0;
k = 1;

while( norm(v) >= tol && k <= maxk)
    [p, ~, ~, ~, lambda] = quadprog(B, gf, [],[],jh, -h);
    lambda = lambda.eqlin;
%   [p, lambda] = pc(B, jh, gf, -h);
    
    % Buscamos que la derivada direccional sea una direccion de descenso
    if (gf'*p - C*norm(h,1) >= 0)
        C = min( 1e5, abs(gf'*p)/norm(h,1) +1);
    end

    %Busqueda de linea
    alpha = 1;
%     phi = feval(fx,x)+C*norm(h,1);
    phi = funcion_merito(x,C, fx, hx);
    D = gf'*p - C*norm(h,1);
    while( funcion_merito(x + alpha*p,C, fx,hx) > phi + alpha*c1*D)
        alpha = alpha/2;
    end

     % Calcular la nueva matriz B
    xn = x + alpha*p; % x_k+1
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
        theta = ((0.8)*aux - s'*y)/(aux-s'*y); % actualizar el vector r
        r = theta*(B*s) + (1-theta)*y;
    end

    B = B - ((B*s*s'*B)/aux) + (r*r')/(s'*r);

    %Si B es mal condicionada
    if(rcond(B)<1e-04)
        B = norm(v)*eye(n);
    end

    lambda = (jhn * jhn')\(jhn*(-gfn));
    x = xn;
    k = k + 1;
    gf = gfn;
    jh = jhn;
    h = feval(hx,x);
    v = [ gf + jh'*lambda; h];

end
end