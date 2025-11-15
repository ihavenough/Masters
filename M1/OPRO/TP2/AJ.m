function gu = AJ(u)
% dans le cas d une fonctionnelle quadratique J, 
% AJ renvoie le gradient de la partie quadratique de J au point u,
% soit la partie lineaire de son gradient,
% c''est a dire la multiplication par la matrice associee a J.
%
global numex
global a b c p
%
switch numex
 case 1
   gu = [- 2*u(2) + 2*u(1); - 2*u(1) + 4*u(2)];
 case 2
   gu = [- 2*u(2) + 2*u(1); - 2*u(1) + 4*u(2)];
 case 3
   gu = [10*u(1) + 7*u(2); 7*u(1) + 5*u(2)];
 case {4, 5, 6}
   error('Fonction non quadratique.') 
 otherwise
   error('Numero de fonction invalide.') 
end

