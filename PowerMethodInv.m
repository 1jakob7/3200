function [x] = PowerMethodInv(A, x0, tolerance)
Iteration = 0;
x_prev = linsolve(A, x0);
success = 0;
for k=1:100    
    disp(['Power iteration ' int2str(k) ':']);
    % Perform power iteration
    y = linsolve(A, x_prev); % change from the original PowerMethod.m
    x = y / norm(y,inf);
    x = x*sign(x(1));
    Iteration = Iteration + 1;
    % Check for convergence
    if norm(x - x_prev, inf) < tolerance
        disp(['Power method converged in ' int2str(Iteration) ' iterations.']);
        success = 1;
        break
    else
        x_prev = x;
    end    
end
if success == 0
    disp(['Power method failed to converge in ' int2str(Iteration) ' iterations.']);
end

