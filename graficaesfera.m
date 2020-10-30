% graficaesfera.m 
% 20 de octubre de 2020
close all

sphere(50)
axis equal
hold on

np = 15;

P = randn(3,np); %  P(:,j) es un vector en R^3
x = zeros(3*np,1);
for j = 1:np
    v = P(:,j); nv = norm(v);
    P(:,j) = v/nv;
    x(3*(j-1)+1:3*j) = P(:,j);
    plot3(P(1,j), P(2,j), P(3,j),'rd','Linewidth',3)
    hold on
end



% optimización con MATLAB
options.MaxFunctionEvaluations= 1e+06;
options = optimset('Algorithm','sqp', 'Display','iter');

[xf,fx] = fmincon('fesfera',x,[],[],[],[],[],[],'hesfera', options);







 
 


