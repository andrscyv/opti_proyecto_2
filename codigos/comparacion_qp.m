%script: comparaciónqp.m
%Comparación entre quadprog.m y qpintpoint.m
%Generar problemas aleatorios de dimensión mayor
%resolvemos con cada código/ calidad de la solución/ cputime
%
%Optimización Numérica
%ITAM
%24 de septiembre de 2020
%
%Q nxn simétrica positiva definida
%A mxn y rango(A)=m
%F pxn
%c en R^n, b en R^m, d en R^p.
%dimensiones n,m,p
%-----------------------------------------------------------
tau=0.75; bag1= []; bag2=[];
for n=20:30:1000
    m=floor(n/2);
    p=floor(n/3);
    [Q ,A, F, c, b, d]=Generapc1(n,m,p,tau);
    t= cputime;
    [x,y,mu,z]=qpintpoint(Q,A,F,c,b,d);
    s= cputime;
    bag1= [bag1; s-t];
    t=cputime;
    [x,fx,exitflag,output]=quadprog(Q,c,-F,-d,A,b);
    s=cputime-t;
    bag2=[bag2; s];
end

%Graficar 
w=[20:30:1000]';
plot(w,bag1,'dr',w,bag1,'b', w, bag2, 'sg', w, bag2, 'm','Linewidth',3)
title('cputime en QP/quadprog vs puntos interiores', 'Fontsize', 16)
xlabel('Dimensión en n', 'Fontsize',16)
ylabel('CPUTIME', 'Fontsize',16)
    