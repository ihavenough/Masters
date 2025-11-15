function [xs, neval] = dichoto(fonction, num, a, b, tol)
    epsilon = tol / 2;
    neval = 0;
    max_iter = 1e5;
    iter = 0;

    while (b - a) > 2 * tol && iter < max_iter
        iter = iter + 1;
        x1 = (a + b) / 2 - epsilon;
        x2 = (a + b) / 2 + epsilon;
        f1 = fonction(num, x1);
        f2 = fonction(num, x2);
        neval = neval + 2;

        if f1 < f2
            b = x2;
        else
            a = x1;
        end
    end

    if iter == max_iter
        warning('Maximum iterations reached');
    end

    xs = (a + b) / 2;
end
