function [x, lambda, k,v] = pcsglobal(fx, hx, x0)
% Metodo de programacion Cuadratica Sucesiva con busqueda de linea,
% usando la funcion de merito L-1 y actualizacion de la hessiana
% con la formula BFGS para el problema
% Min fx
% Sujeto a hx = 0
%
% Entradas
% fx y hx son cadenas de caracteres con las funciones en Matlab
% de la funcion objetivo y las restricciones del problema
% El vector x0 es el valor inicial
%
% Salida
% x.- aproximacion al minimo local
% lambda.- multiplicador de Lagrange asociado a x.
% k.- numero de iteraciones realizadas.
%
% Debe usar las funciones: gradiente.m y jacobiana.m para calcular
% las primeras derivadas.
% -------------------------------------------------------------------------
% Equipo: 
% Javier Montiel Gonzalez 159216
% Alexis Ayala Redon 156916
% Andres Cruz y Vera 155899

% Parámetros inciales
tol = 1e-5;     % Tolerancia de las cnpo
maxk = 100;     % Máximo número de iteraciones
c1 = 1e-2;      % Valor de la constante para la búsqueda de línea
C = 1;          % Parámetro para encontrar una derivada direccional negativa

h = feval(hx, x0); % El valor de las restricciones en el punto inicial
n = length(x0); %Dimension de la variable
m = length(h); % Dimension de las restricciones


lambda = zeros(m,1); %Inicializamos el multiplicador de Lagrange
B = eye(n); % Inicializamos la aproximacion a la Hessiana

gf = gradiente(fx,x0);
jh = jacobiana(hx, x0);

% Condiciones necesarias de primer orden
v = [ gf + jh'*lambda; h];
x = x0;
k = 1;

while( norm(v) >= tol && k <= maxk)
    
    % Resolvemos el subproblema cuadrático con restricciones con la funcion
    % de puntos interiores de MatLab
    [p, ~, ~, ~, lambda] = quadprog(B, gf, [],[],jh, -h);
    lambda = lambda.eqlin;
    
    % Buscamos que la derivada direccional en p sea una direccion de
    % descenso
    norm1H=norm(h,1);
    if (gf'*p - C*norm1H >= 0)
        C = min( 1e5, abs(gf'*p)/norm1H +1);
    end

    %Busqueda de linea
    alpha = 1;
    phi = funcion_merito(x,C, fx, hx);
    D = gf'*p - C*norm1H;
    while( funcion_merito(x + alpha*p,C, fx,hx) > phi + alpha*c1*D)
        alpha = alpha/2;
    end

    % Calcular la nueva matriz B con el metodo BFGS con la actualizacion de
    % Powell
    xn = x + alpha*p; 
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
        B = eye(n);
    end

    % Calculamos el nuevo multiplicador de Lagrange
    lambda = (jhn * jhn')\(jhn*(-gfn));
    
    % Actualizamos los valores
    x = xn;
    k = k + 1;
    gf = gfn;
    jh = jhn;
    h = feval(hx,x);
    v = [ gf + jh'*lambda; h];
end
end