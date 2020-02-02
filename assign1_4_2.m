% Jakob Horvath, u1092049
% Plots the errors when using cubic spline and PCHIP polynomial interpolations of the GelbTanner function for polynomials
% of degree 28, 36, 44, 52, and 60 with 1001 sample points. L1, L2 and
% L-infinity norms are used.

sp = 1001;
L1 = zeros(5);
L2 = zeros(5);
Linf = zeros(5);

yTrue = ones(1001, 1) .* GelbTanner(linspace(-1, 1, 1001)); % obtain true y values from GelbTanner graph

n = [29, 37, 45, 53, 61]; % polyinterp uses degrees of length(x)-1
for k=1:5
    x = linspace(-1, 1, n(k));
    y = GelbTanner(x);
    u = linspace(-1, 1, sp);
    s = spline(x, y, u);
    p = pchip(x, y, u);
    plot(u, s, '-');
    hold on;
    plot(u, p, '-');
    hold on;

    % obtain the errors between y-true and poly estimate for each sample point
    sErr = zeros(1, sp);
    pErr = zeros(1, sp);
    for i=1:length(sErr)
        sErr(i) = sErr(i) + (yTrue(i) - s(i));
        pErr(i) = pErr(i) + (yTrue(i) - p(i));
    end
    %L1 norms
    sL1(k) = norm(sErr, 1);
    pL1(k) = norm(pErr, 1);
    %L2 norms
    sL2(k) = norm(sErr, 2);
    pL2(k) = norm(pErr, 2);
    %L3 norms
    sLinf(k) = norm(sErr, inf);
    pLinf(k) = norm(pErr, inf);
end

% plot the L1 error rates for spline and pchip
figure, plot(n(1:5), sL1(1:5), '-o');
hold on;
plot(n(1:5), pL1(1:5), '-o');
legend('sL1', 'pL1');
% plot the L2 and L-infinity error rates for spline and pchip
figure, plot(n(1:5), sL2(1:5), '-x');
hold on;
plot(n(1:5), pL2(1:5), '-x');
hold on;
plot(n(1:5), sLinf(1:5), '-*');
hold on;
plot(n(1:5), pLinf(1:5), '-*');
legend('sL2', 'pL2', 'sLinf', 'pLinf');