% Jakob Horvath, u1092049
% Plots the error when using Lagrange polynomial interpolation of the hubbard function for polynomials
% of degree 28, 36, 44, 52, and 60 with 1001 sample points. L1, L2, and
% L-infinity norms are used.

sp = 1001;
t = 0.25;
L1 = zeros(5);
L2 = zeros(5);
Linf = zeros(5);

yTrue = ones(1001, 1) .* hubbard(linspace(0, 1, 1001), t); % obtain true y values from hubbard graph

n = [29, 37, 45, 53, 61]; % polyinterp uses degrees of length(x)-1
for k=1:5
    x = linspace(0, 1, n(k));
    y = hubbard(x, t);
    u = linspace(0, 1, sp);
    v = polyinterp(x, y, u);
    plot(u, v, '-');
    hold on;

    % obtain the error between y-true and poly estimate for each sample point
    err = zeros(1, sp);
    for i=1:length(err)
        err(i) = err(i) + (yTrue(i) - v(i));
    end
    %L1 norm
    L1(k) = norm(err, 1);
    %L2 norm
    L2(k) = norm(err, 2);
    %L3 norm
    Linf(k) = norm(err, inf);
end

% plot the error rates for each norm
figure, plot(n(1:5), L1(1:5), '-o');
hold on;
plot(n(1:5), L2(1:5), '-x');
hold on;
plot(n(1:5), Linf(1:5), '-*');
legend('L1', 'L2', 'Linf');
