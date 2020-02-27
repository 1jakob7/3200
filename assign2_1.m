% Jakob Horvath, u1092049
% Implements the composite Simpson formula for x^p and 1+sin(x)*cos(2x/3)*sin(4x).

A = zeros(6, 6);
N = [8, 16, 32, 64, 128, 256]; % corresponds to 17,33,65,129,257,513
p = [2, 3, 4, 5, 6, 8];
a = 0;
b = 1;
for k=1:6 % loop thru N's
    for v=1:6 % loop thru p's
        area = 0;
        for i=1:2*N(k)+1 % summation from i=1 to 2N+1
            dx = (b-a)/(2*N(k)); % delta-x ('h')
            xi = a+(i-1)*dx;
            if i == 1 || i == 2*N(k)+1 % first or last element in summation
                wi = dx/3;
            elseif mod(i, 2) == 0 % an even element
                wi = 4*dx/3;
            else % an odd element
                wi = 2*dx/3;
            end
            area = area + wi*xi^p(v); % x^p
        end 
        A(k, v) = area; % add estimated area for the integral to a 2d array
    end
end

B = zeros(6, 1);
a = 0;
b = 2*pi;
fprintf('   2*N(k)     area \n');
for k=1:6 % loop thru N's
   area = 0; 
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
       area = area + wi*(1+sin(xi)*cos(2*xi/3)*sin(4*xi)); % 1+sin(x)*cos(2x/3)*sin(4x)
   end
   B(k, 1) = area; % add estimated area for the integral to a 1d array
   fprintf('%8.0e %21.14f \n', 2*N(k), area)
end