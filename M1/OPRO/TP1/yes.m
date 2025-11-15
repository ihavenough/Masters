
k_values = 5:0.5:10; % instead of 5:0.5:10
tol_values = 10.^(-k_values);

neval1 = zeros(size(tol_values));
neval5 = zeros(size(tol_values));

% Evaluate for function num = 1 and 3
for i = 1:length(tol_values)
    [~, neval1(i)] = dichoto(@fonc, 1, 0, 3, tol_values(i));  % Function 1
    [~, neval5(i)] = dichoto(@fonc, 5, 0, 3, tol_values(i));  % Function 3
end

% Plot for function 1
figure;
plot(k_values, neval1);
xlabel('log_{10}(1/\epsilon)');
ylabel('Nombre d''évaluations');
title('Fonction 1 — Méthode de dichotomie');


% Plot for function 3
figure;
plot(k_values, neval3);
xlabel('log_{10}(1/\epsilon)');
ylabel('Nombre d''évaluations');
title('Fonction 5 — Méthode de dichotomie');
