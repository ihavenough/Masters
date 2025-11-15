clear; clc; close all;

% Define A and b
A1 = [2 -2; -2 4];
b1 = [0; 2];

A3 = [10 7; 7 5];
b3 = [-2; -1];

% Exact solutions
u_star1 = A1\b1;
u_star3 = A3\b3;

% Conditioning
r1 = cond(A1);
r3 = cond(A3);
fprintf('Conditionnement A1: %.4f\n', r1);
fprintf('Conditionnement A3: %.4f\n', r3);

% Initial guess
u0 = [0; 0];

maxIter = 20;

%% === PROBLEM 1 ===
A = A1; b = b1; u_star = u_star1;
uk = u0;
tslesuk1 = uk;

for k = 1:maxIter
    wk = A*uk - b;                       % gradient
    alphak = (wk'*wk) / (wk'*(A*wk));    % optimal step
    uk = uk - alphak * wk;               % update
    tslesuk1 = [tslesuk1 uk];
end

%% === PROBLEM 3 ===
A = A3; b = b3; u_star = u_star3;
uk = u0;
tslesuk3 = uk;

for k = 1:maxIter
    wk = A*uk - b;
    alphak = (wk'*wk) / (wk'*(A*wk));
    uk = uk - alphak * wk;
    tslesuk3 = [tslesuk3 uk];
end

%% === Compute A-norm errors ===
for k = 1:size(tslesuk1,2)
    e1 = tslesuk1(:,k) - u_star1;
    normA1(k) = sqrt(e1' * A1 * e1);
end

for k = 1:size(tslesuk3,2)
    e3 = tslesuk3(:,k) - u_star3;
    normA3(k) = sqrt(e3' * A3 * e3);
end

%% === Theoretical bounds ===
bound1 = ((r1-1)/(r1+1)).^(0:size(tslesuk1,2)-1) ...
    * sqrt((u0-u_star1)'*A1*(u0-u_star1));
bound3 = ((r3-1)/(r3+1)).^(0:size(tslesuk3,2)-1) ...
    * sqrt((u0-u_star3)'*A3*(u0-u_star3));

%% === Plot ===
figure;
hold on;
semilogy(0:size(tslesuk1,2)-1, normA1, 'b-o', 'DisplayName','||uk-u*||A P1');
semilogy(0:size(tslesuk3,2)-1, normA3, 'r-o', 'DisplayName','||uk-u*||A P3');
semilogy(0:size(tslesuk1,2)-1, bound1, 'b--', 'DisplayName','Borne théorique P1');
semilogy(0:size(tslesuk3,2)-1, bound3, 'r--', 'DisplayName','Borne théorique P3');
xlabel('k (itération)');
ylabel('log||uk - u*||_A');
legend show;
grid on;
title('Convergence du gradient à pas optimal et influence du conditionnement');
hold off;
