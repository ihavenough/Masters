% Ce script met en oeuvre la methode du gradient a pas optimal.
%
% Cas quadratique elliptique et cas general.
%
%
clear all; close all;
global numex
global a b c p
numex=4;    % (1 2 3 4 5 ou 6)
quadratique = 'NON'; %('OUI' ou 'NON')
a=4; b=4; c=10; p=10;
%
% Initialisation des variables u0, epsil, nitermax, niter utilisees ci-dessous :
%
switch numex
  case 1
    u0 = [0; 0];
  case 2
    u0 = [0; 0];
  case 3
    u0 = [0; 0];
  case 4
    u0 = [-1; 2];
  case 5
    u0 = [1; 0];   % essayer [1; 1]
  case 6
    u0 = [0; 1];
  otherwise
    error('numex doit etre compris entre 1 et 6')
end

epsil = 1e-6;      % tolérance
nitermax = 1000;   % nombre max d’itérations
niter = 0;

%
stopCritere = 1; % Servira pour le test d''arret.
uk=u0           % uk contiendra les iteres successifs
tslesuk = u0;   % stockera l''ensemble des iteres successifs
if (quadratique=='OUI')
%
  while((stopCritere > epsil) & (niter < nitermax))
    niter = niter + 1;
    wk = GJ(uk);
    Awk = AJ(wk);
    alphak = (wk'*wk) / (wk'*Awk);
    uk = uk - alphak * wk;
    tslesuk = [tslesuk, uk];
    stopCritere = norm(GJ(uk));
  end

elseif (quadratique=='NON')
  %
  nbevalinfmin = 0;
  while((stopCritere > epsil) & (niter < nitermax))
   niter = niter+1;
   wk = GJ(uk) ;
   myfunction = @(x)J(uk - x * wk);
   [alphak, fval, exitflag, output] = fminbnd(myfunction,0,4); %CHOISIR ...);
   nbevalinfmin = nbevalinfmin + output.iterations;
   uk = uk - alphak * wk;
   tslesuk = [tslesuk, uk];
   stopCritere = norm(GJ(uk));
  end
  nbevalinfmin,niter
  tslesuk
  %
else
 error('Specifiez par OUI ou NON le caractere quadratique de la fonctionnelle.')
end
xmin = min(min(tslesuk(1,:)),-3.3);
xmax = max(tslesuk(1,:));
xdiff = xmax-xmin;
ymin = min(tslesuk(2,:));
ymax = max(tslesuk(2,:));
ydiff = ymax-ymin;
%xdiff = max(xdiff,ydiff); ydiff = xdiff;
xmin = xmin - 0.5*xdiff; xmax = xmax + 0.5*xdiff;
ymin = ymin - 0.5*ydiff; ymax = ymax + 0.5*ydiff;

figure(1)

visiso(xmin, xmax, ymin, ymax);
visiter(uk,tslesuk,niter,nitermax);

axis equal

