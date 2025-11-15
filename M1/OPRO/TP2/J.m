function ju = J(u)
% renvoie la valeur de la fonctionnelle J au point u
global numex
global a b c p
%
switch numex
 case 1
   ju = - 2*u(1)*u(2) - 2*u(2) + u(1)^2 + 2*u(2)^2;
 case 2
   ju = - 2*u(1)*u(2) - 4*u(1) + u(1)^2 + 2*u(2)^2;
 case 3
   ju = 5*u(1)^2 + (5/2)*u(2)^2 + 7*u(1)*u(2) + 2*u(1) + u(2);
 case 4
   ju = c*( u(1)^2 + u(2)^2 + u(1)*u(2) + a*(u(2) - u(1)) )^2 + (1 - b*(u(1)+u(2)))^2;
 case 5
   ju = u(1)^2 - u(2)^2;
 case 6
   ju = ( u(1) - 1 )^2 + p*( u(1)^2 - u(2) )^2;
 otherwise
   error('Numero de fonction invalide.') 
end
