function [xs, neval] = secdor(fonction, num, a, b, tol, gamma)
    delta = gamma / (1 - gamma);
    x1 = a + gamma * (b - a);
    x2 = a + (1 - gamma) * (b - a);
    f1 = fonction(num, x1);
    f2 = fonction(num, x2);
    neval = 2;

    while (b - a) > tol
        if f1 <= f2
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + delta * (b - a);
            f1 = fonction(num, x1);
        else
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + (1 - gamma) * (b - a);
            f2 = fonction(num, x2);
        end
        neval = neval + 1;
    end

    xs = (a + b) / 2;
end
