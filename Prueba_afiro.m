%script: Prueba_afiro

%afiro.mat
load afiro
A=full(A);
n=length(c);
Q=eye(n);

x = metespnulo(Q,A,c,b);

[x1, y] = pc(Q,A,c,b);

norm(x-x1)
