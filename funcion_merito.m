function [phi] = funcion_merito(x,C, fx, hx)
%Función de mérito para min f(x) s.a h(x) = 0

phi = feval(fx,x) + C * norm(feval(hx,x),1)
end

