function [x, lambda, k] = pcsglobal(fx, hx, x0)
% Método de programación Cuadrática Sucesiva con búsqueda de linea,
% usando la función de mérito L-1 y actualización de la hessiana
% con la fórmula BFGS para el problema
% Min fx
% Sujeto a hx = 0
%
% fx y hx son cadenas de caracteres con las funciones en Matlab
% de la función objetivo y las restricciones del problema
% El vector x0 es el valor inicial
% Salida
% x.- aproximación al mínimo local
% λ.- multiplicador de Lagrange asociado a x.
% k.- número de iteraciones realizadas.
%
% Debe usar las funciones: gradiente.m y jacobiana.m para calcular
% las primeras derivadas.
%
% Javier Montiel Gonzalez 159216, Alexis Ayala Redón 156916
% Andrés Cruz y Vera 155899

    tol = 1e-5;
    maxk = 100;
    c1 = 1e-2;
    C = 1;
    h = feval(hx, x0); 
    n = length(x0); %Dimensión de la variable
    m = length(h); % Dimensión de las restricciones
    lambda = zeros(m,1);
    B = eye(n);
    % Condiciones necesarias de primer orden
    gf = gradiente(fx,x0);
    jh = jacobiana(hx, x0);
    v = [ gf + jh'*lambda; h];
    x = x0;
    while( norm(v) >= tol && k <= maxk )
        [p, ~, ~, ~, lambda] = quadprog(B, gf, [],[],jh, -h);
        if (gf'*p - C*norm(h,1) >= 0)
            C = min( 1e5, abs(gf'*p)/norm(h,1) +1);
        end
        alpha = 1;
        phi = feval(f,x)-C*norm(h,1);
        phi_paso = feval(f, x + p) - C*norm(h,1);
        D = gf'*p - C*norm(h,1);
        while( phi_paso > phi + alpha*c1*D)
            alpha = alpha/2;
            phi_paso = feval(f,x + alpha*p) - C*norm(h,1);
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
            theta = ((0.2)*aux - s'*y)/(aux-s'*y); % actualizar el vector r
            r = theta*(B*s) + (1-theta)*y;
        end

        B = B - ((B*s*s'*B)/aux) + (r*r')/(s'*r);

        %Si B es mal condicionada
        if(rcond(B)<1e-04)
            B = norm(vcnpo)*eye(n);
        end
        
        lambda = (jhn * jhn')\(jhn*(-gfn);
        x = xn;
        k = k + 1;
        gf = gfn;
        jh = jhn;
        h = feval(hx,x);
        v = [ gf + jh'*lambda; h];
        
    end
end