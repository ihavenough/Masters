%% Problem P6:  f(x1,x2) = (x1-1)^2 + p*(x1^2 - x2)^2

clear; clc;

p = 10;
u0 = [0; 1];
alpha = 0.01;
tol = 1e-6;
maxIter = 1000;

f = @(u) (u(1)-1)^2 + p*(u(1)^2 - u(2))^2;

gradf = @(u) [ 2*(u(1)-1) + 4*p*u(1)*(u(1)^2 - u(2));
              -2*p*(u(1)^2 - u(2)) ];

u_star = [1; 1];
u = u0;

errors = zeros(1, maxIter+1);
errors(1) = norm(u - u_star);

for k = 1:maxIter
    g = gradf(u);
    u = u - alpha*g;
    errors(k+1) = norm(u - u_star);

    if errors(k+1) < tol
        errors = errors(1:k+1);  % truncate
        break;
    end
end

figure;
semilogy(0:length(errors)-1, errors, 'LineWidth', 1.5);
xlabel('k (itération)');
ylabel('log ||u_k - u^*||');
title('Convergence pour P6: (x1-1)^2 + p(x1^2-x2)^2   avec p=10');
grid on;
print('-dpng', '-r300', 'convergence_p6.2.png');

