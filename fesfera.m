function [f] = fesfera(x)
 warning('off','all');
% Funcion de repulsion para (n/3) puntos en la esfera unitaria de 
% dimension tres.
%
% Optimizacion Numerica 
% ITAM
% 20 de octubre de 2020
%
n = length(x);
np = floor(n/3); % numero de puntos en la esfera
f = 0;           % valor inicial de f

for i = 1:np-1
    ui = x(3*(i-1)+1:3*i);         % anclados en el punto ui
   for j =i+1:np
        uj = x(3*(j-1)+1:3*j);      % punto uj
      f = f+(1/norm(ui-uj));    % la funcion crece en sumandos.
   end    
end

end
