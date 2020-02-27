% Jakob Horvath, u1092049
% Uses Simpson's Rule to evaluate the function, 'func2_1', using many
% intervals. Then, error estimates are inserted into quadtx's 'tol'
% parameter to see how many steps are required.

% Simpson's
N = [16, 32, 64, 128, 256, 512]; % corresponds to 33,65,129,257,513,1025
B = zeros(6, 1);
a = 0;
b = 3;
for k=1:6 % loop thru N's
   areaS = 0; 
   for i=1:2*N(k)+1 % summation from i=1 to 2N+1
       dx = (b-a)/(2*N(k)); % delta-x ('h')
       xi = a+(i-1)*dx;
       if i == 1 || i == 2*N(k)+1 % first or last element
           wi = dx/3;
       elseif mod(i, 2) == 0 % even element
           wi = 4*dx/3;
       else
           wi = 2*dx/3;
       end
       areaS = areaS + wi*func2_1(xi); % (cos(x^3))^200
   end
   B(k, 1) = areaS; % add estimated area for the integral to a 1d array
end

% obtain error estimates for Simpson's Rule -- (1/15)(Simp(h/2)-Simp(h))
E = zeros(5, 1);
E(1) = abs((1/15)*(B(2)-B(1)));
E(2) = abs((1/15)*(B(3)-B(2)));
E(3) = abs((1/15)*(B(4)-B(3)));
E(4) = abs((1/15)*(B(5)-B(4)));
E(5) = abs((1/15)*(B(6)-B(5)));

% use quadtx when tol is equal to a given error estimate from Simpson
for k=1:5
    tol = E(k);
    [Q, fcount] = quadtx(@func2_1, a, b, tol);
    fprintf('%23.10f %7d \n', Q, fcount)
end
