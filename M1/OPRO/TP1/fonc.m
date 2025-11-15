function f = fonc(num, x)
    switch num
        case 1
            f = (x - 1).^2;
        case 2
            f = zeros(size(x));
            f(x <= 1) = (x(x <= 1) - 1).^2 + 2;
            idx = (x > 1) & (x <= 2);
            f(idx) = 1 - x(idx)/2;
            f(x > 2) = 0.5 + (x(x > 2) - 2).^2;
        case 3
            f = abs(x - 1) .* (1.1 - sin(6 * x));
        case 4
            f = x.^2;
        case 5
            f = (x - 3).^2;
        otherwise
            error('Function number not recognized');
    end
end
