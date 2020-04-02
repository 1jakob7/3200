% Jakob Horvath
% Finds an estimate for the vector 'b' in a linear equation of the
% form 'Ax=b' using the Gauss-Seidel iterative method with 
% under-relaxtion applied.

N = [16, 32, 64]; % rank of matrix A
tol = [1.e-3, 1.e-4, 1.e-5]; % range of tolerances to test against
X = zeros(N(length(N)), length(N)*length(tol)); % hold vector x estimate
B = zeros(N(length(N)), length(N)*length(tol)); % hold A*x estimate
I = zeros(1, length(N)*length(tol)); % hold the number of iterations 

omega = 0.45;

for i=1:length(N)
    % Construct matrix A %
    e = ones(N(i), 1);
    A = spdiags([e -4*e 6*e -4*e e], -2:2, N(i), N(i));
    % Adjust main diagonal entries
    A(1, 1) = 9;
    A(N(i)-1, N(i)-1) = 5;
    A(N(i), N(i)) = 1;
    % Adjust offset diagonal entries
    A(N(i)-1, N(i)) = 2;
    A(N(i), N(i)-1) = -2;
    
    for k=1:length(tol)
        itr = 0;
        normVal = Inf;
        
        % Construct vectors b and x %
        b = 1/(N(i)^4)*ones(N(i), 1);
        x = zeros(N(i), 1);
        x_new = zeros(N(i), 1); % needed for computing 'c' (convergence rate)
        
        while normVal > tol(k)
            for p=1:N(i)
                sigma = 0;
                for j = 1:p-1
                    sigma = sigma + A(p, j)*x_new(j);
                end
                for j = p+1:N(i)
                    sigma = sigma + A(p, j)*x_new(j);
                end
                x_new(p) = (1/A(p, p))*(b(p)-sigma);
            end
            itr = itr + 1;
            if (itr > 1)
                x_new = omega*x_new + (1-omega)*x; % apply under-relaxation
                c = norm(x_new - x)/norm(x - x_old);
                normVal = (c/(1-c))*norm(x_new - x);
            end
            x_old = x;
            x = x_new;
        end
        b_est = A*x_new;
        for q=1:N(i)
            X(q, (i-1)*length(tol)+k) = x_new(q);
            B(q, (i-1)*length(tol)+k) = b_est(q);
        end
        I(1, (i-1)*length(tol)+k) = itr;
    end
end