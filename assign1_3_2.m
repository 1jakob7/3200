% Jakob Horvath, u1092049
% Plots the error when using Chebyshev polynomial interpolation of the GelbTanner function for polynomials
% of degree 28, 36, 44, 52, and 60 with 1001 sample points. L2 and
% L-infinity norms are used.

sp = 1001;
L1 = zeros(5);
L2 = zeros(5);
Linf = zeros(5);

yTrue = ones(1001, 1) .* GelbTanner(linspace(0, 1, 1001)); % obtain true y values from GelbTanner graph

n = [29, 37, 45, 53, 61]; % polyinterp uses degrees of length(x)-1
for k=1:5
    x = linspace(-1, 1, n(k));
    y = GelbTanner(x);
    
    xEqual = zeros(1, n(k));
    for m=1:n(k)
        xEqual(m) = -cos((2.0*m-1.0)*pi/(2.0*n(k)))/cos(pi/(2.0*n(k)));
    end
    yEqual = GelbTanner(xEqual);
    
    u = linspace(-1, 1, sp);
    v = polyinterp(xEqual, yEqual, u);
    plot(u, v, '-');
    hold on;
    
    % obtain the error between y-true and poly estimate for each sample point
    err = zeros(1, sp);
    for i=1:length(err)
        err(i) = err(i) + (yTrue(i) - v(i));
    end
    %L2 norm
    L2(k) = norm(err, 2);
    %L3 norm
    Linf(k) = norm(err, inf);
end

% plot the error rates for each norm
figure, plot(n(1:5), L2(1:5), '-x');
hold on;
plot(n(1:5), Linf(1:5), '-*');
legend('L2', 'Linf');