% Jakob Horvath, u1092049
% Calls quadtx with the function, 'func2_1', to record the estimated
% integral, the number of calls made to the function, and the time spent in
% quadtx. Then calls quadtx2, which checks the error to be less than
% 'tol*(h/(b-a))' as opposed to just 'tol', and gathers similar data.

T = zeros(8, 1);
fprintf('   tol        Q                   fcount \n')
a = 0;
b = 3;
% regular quadtx
for  k = 7:14
    tol = 10^(-k);
    tic % start timer
    [Q, fcount] = quadtx(@func2_1, a, b, tol);
    T(k-6) = toc; % store timing data for previous function call
    fprintf('%8.0e %21.14f %7d \n', tol, Q, fcount)
end

fprintf('----------------------------------------- \n')
fprintf('   tol2       Q2                  fcount2 \n')
% modified quadtx w/ error being less than tol(h/(a-b))
for  k = 7:14
    tol2 = 10^(-k);
    [Q2, fcount2] = quadtx2(@func2_1, a, b, tol2);
    fprintf('%8.0e %21.14f %7d \n', tol2, Q2, fcount2)
end