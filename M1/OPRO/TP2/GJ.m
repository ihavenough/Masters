function gu = GJ(u)
% renvoie la valeur du gradient de la fonctionnelle J au point u
global numex
global a b c p
%
switch numex
 case 1
   gu = [- 2*u(2) + 2*u(1); - 2*u(1) - 2 + 4*u(2)];
 case 2
   gu = [ -2*u(2) - 4 + 2*u(1) ; -2*u(1) + 4*u(2) ];
 case 3
   gu = [ 10*u(1) + 7*u(2) + 2 ; 7*u(1) + 5*u(2) + 1 ];
 case 4
   gu = 2*[ c*( u(1)^2 + u(2)^2 + u(1)*u(2) + a*(u(2) - u(1)) )*(2*u(1)+u(2)-a) - b*(1 - b*(u(1)+u(2))) ;
            c*( u(1)^2 + u(2)^2 + u(1)*u(2) + a*(u(2) - u(1)) )*(2*u(2)+u(1)+a) - b*(1 - b*(u(1)+u(2)))  ];
 case 5
   gu = [ 2*u(1); - 2*u(2)];
 case 6
   gu = 2*[ u(1) - 1 + p*( u(1)^2 - u(2) )*(2*u(1)) ; 
            p*( u(2) - u(1)^2 )];
 otherwise
   error('Numero de fonction invalide.') 
end
