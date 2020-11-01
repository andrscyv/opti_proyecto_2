% Resolver el problema de colocar np puntos en la esfera unitaria de tal
% forma que se minimice la funcion de repulsion entre ellos

% Equipo:
% Javier Montiel Gonzalez 159216
% Alexis Ayala Redon 156916
% Andres Cruz y Vera 155899

clear all; clc;
format long e;
% Graficamos la esfera unitaria con 50 meridianos
sphere(50)
axis equal
hold on

% Definimos el numero de puntos que se colocaran en la esfera
np=21;


P = randn(3,np); %  P(:,j) es un vector en R^3
x = zeros(3*np,1);
% Fijamos el primer punto como (1, 0, 0)'
x(1)=1;
% Generamos de manera aleatoria los np-1 puntos restantes
for j = 2:np
    v = P(:,j); nv = norm(v);
    P(:,j) = v/nv;
    x(3*(j-1)+1:3*j) = P(:,j);
end

% Encontramos la aproximacion al minimo local con pcsglobal
inicio= cputime;
[x, lambda, k, v] = pcsglobal('fesfera', 'hesfera', x); 
tiempo1= cputime-inicio;
fprintf('Valor de las CNPO en el mínimo %2.8f', norm(v))

% Graficamos los puntos en la esfera unitaria
 n=length(x);
 for j = 1:3:n-2
    plot3(x(j), x(j+1), x(j+2),'rd','Linewidth',3);
 end
title('Aproximación al minimo local con pcsglobal')
 
% Encontramos la aproximacion al minimo local con el metodo de MATLAB
options.MaxFunctionEvaluations= 1e+06;
options = optimset('Algorithm','sqp');
inicio2= cputime;
[xm,fxm] = fmincon('fesfera',x,[],[],[],[],[],[],'hesferam', options);
tiempo2= cputime-inicio2;

%Tabla de compracion de metodos
metodo=["pcsglobal"; "fmincon"];
fx=[fesfera(x); fxm];
tiempo= [tiempo1; tiempo2];

comp=table(metodo, fx, tiempo);
disp(comp);