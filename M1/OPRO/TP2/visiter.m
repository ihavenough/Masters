function visiter(u,Tabuk,iter,itermax)
%
% permet de visualiser les iterations de l'algorithme du gradient
%
% variables entree:
% u       = solution obtenue  argmin J(v)
% Tabuk    = tableau contenant la suite (u_k) des points intermediaires
%           u_k = (Tabuk(1,k),Tabuk(2,k))
% iter    = nb d'iterations effectuees
% itermax = nb max d'iterations autorisees


global isov   % tableau des isovaleurs retournees par la fonction visiso.m 
global numex  % numero de l'exemple traite

hold on
plot(Tabuk(1,:),Tabuk(2,:),'ko')
hold on
plot(Tabuk(1,:),Tabuk(2,:),'k-')
%clc,home
disp(['Pour l''exemple numero ',int2str(numex)])
disp(blanks(2)')
disp(['le minimum trouve est : ']),u

disp(['la fonctionnelle en ce point vaut : ',num2str(J(u))])
disp(['||grad J(u*)|| =  ',num2str(norm(GJ(u),2))])
if (iter < itermax)
  disp(['le nb d''iterations effectue est : ',num2str(iter)])
else
  disp(['le nb d''iterations maximum (',num2str(iter),') a ete atteint'])
end
disp(blanks(2)')

text(u(1),u(2),['  argmin =  (',num2str(u(1)),' , ', num2str(u(2)), ')'])
xlabel(['nb d''itérations = ',num2str(iter), ...
        '        ||grad J(u*)|| =  ',num2str(norm(GJ(u),2)), ...
        '          J(u*) =  ',num2str(J(u))])
rep = input ('Voulez-vous afficher les isovaleurs (o/n) : ','s');
if rep=='o'
  % Comment the line below and uncomment the next one if using Octave.
   clabel(isov,'manual'); % OK for Matlab but not for Octave
  % clabel(isov); % Use this line if using Octave.
end
