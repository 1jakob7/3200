% Jakob Horvath, u1092049
% Uses the supplied Gradient Descent algorithm to compute a least squares
% approximation for the US Census data (var x and y). The current program 
% uses what was found to be the optimal polynomial degree and alpha value
% for making the most accurate / work reducing estimate.

%%% US Census data
x = [1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000]';
y = [75.995, 91.972, 105.711, 123.203, 131.669, 150.697, 179.323, 203.212, 226.505, 249.633, 281.422]';

% scale the year values to be in the range: [-1,1]
for s=1:length(x)
    x(s) = (x(s) - 1950) / 50;
end

re = zeros(1, 5); % store the residuals
n = 10;
np1 = n+1;
%m = 3;
for m=1:5 % best polynomial degree for accuracy and the residual
    mp1 = m+1;
    XX = zeros(np1,mp1);
    for i = 1:np1
        for j = 1:mp1
            XX(i,j) = x(i)^(j-1);
        end
    end
    %solution vectors
    a = zeros(mp1,1);
    aold = zeros(mp1,1);
    %gradient vector
    r = zeros(mp1,1);
    res = y - XX*a;
    normres = norm(res);
    normVal = Inf;
    itr = 0;
    tol = 1.e-5;
    fac = 2.0/n;
    alpha = 0.69; % optimal value for reducing the work
     
    %%% Algorithm: Steepest Descent
    while normVal>tol
        aold=a;
        res = y-XX*a;
        normres= norm(res);
        for i = 1:mp1
            r(i) = 0;
            for j = 1:np1
                r(i) = r(i)+res(j)*XX(j,i);
            end
            r(i)= r(i)*fac;
        end
        a = a + alpha*r;
        itr=itr+1;
        normVal=abs(aold-a);
    end
    
    fprintf('\n degree = %i \n', m);
    %%% Print and plot results
    fprintf('Solution of the system is %f %f %f %f %f %f \n  ',a);
    %
    fprintf('\n norm residual %f in %i iterations \n',normres,itr);
    re(m) = normres;
    z = linspace(x(1), 1.5, 251);
    fz = zeros(251,1);
    for i = 1:251
        fz(i)=0.0;
        for j = 1:mp1
            fz(i) = fz(i)+a(j)*z(i)^(j-1);
        end
    end
    % 221 = 2010 and 239 = 2019
    fprintf(' Estimate for 2010 pop: %f \n', fz(221));
    fprintf(' Estimate for 2010 pop: %f \n', fz(239));
    
    plot(x,y,z,fz)
end
