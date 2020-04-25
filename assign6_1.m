% Jakob Horvath, u1092049
% Using a data set from a supplied text file, utilizes either normal
% equations or QR decomposition to produce a polynomial estimation with
% 1001 sample points.

warning('off','MATLAB:nearlySingularMatrix'); % triggered by normal equations method
warning('off','MATLAB:rankDeficientMatrix'); % triggered by QR when not properly scaled

n = 82;
p = 11;

D = zeros(82, 2);
count = 0;
file = fopen('data.txt', 'r'); % open file for 'r'eading
while ~feof(file)
    count = count + 1;
    tline = fgetl(file);
    % format and store data points
    extracted = strsplit(tline);
    % assign values in reverse to unswap x and y
    D(count, 2) = str2double(extracted(1));
    D(count, 1) = str2double(extracted(2));
end
fclose(file);

figure(1)
D = sortrows(D, 1);
%plot(D(:,1), D(:,2), '+');

%figure(2)
rescale = 1;

if (rescale == 1) % 0 leaves data as is, 1 scales both x and y to fall within [-1, 1]
    Ma = max(D, [], 1); % max x == Ma(1), max y == Ma(2)
    Mi = min(D, [], 1); % min x == Mi(1), min y == Mi(2)
    for i=1:n
        D(i, 1) = (2/(Ma(1)-Mi(1)))*D(i,1)-(Ma(1)+Mi(1))/(Ma(1)-Mi(1));
        D(i, 2) = (2/(Ma(2)-Mi(2)))*D(i,2)-(Ma(2)+Mi(2))/(Ma(2)-Mi(2));
    end
end
plot(D(:,1), D(:,2), '+');
hold on;

%%% Monomial Vandermonde Matrix Approx.
yyy = zeros(82, 1);
nxx =1001;
for m=9:2:25
    
    imeth=1;  % 0 gives normal equations, 1 gives QR
    
    if (rescale == 1)
        x = linspace(-1,1,n);
        xx = linspace(-1,1,nxx);
    else
        x = linspace(-9,-3,n);
        xx = linspace(-9,-3,nxx);
    end
    
    vv=zeros(nxx,1);
    A = zeros(n,m);
    for jj = 1:n
        for ii = 1:m
            A(jj,ii) = x(jj)^(ii-1);
        end
    end
    if(imeth == 0)
        C = A';
        yy= C*D(:,2);
        B = C*A;
        v= B\yy;
        CO =cond(B);
        w = warning('query','last');
        id = w.identifier;
        %fprintf(' n= %i m = %i cond = %7.2e \n',n,m,CO)
    else
        [Q,R]=qr(A);
        yy = Q'*D(:,2);
        v = R\yy;
    end
    elsqsum=0.0;
    for ii = 1:nxx
        for jj=1:m
            vv(ii) = vv(ii)+v(jj)*xx(ii)^(jj-1);
        end
    end
    for i=1:n
        yy = interp1(xx,vv,D(i,1)); % get estimated value from poly data
        yyy(i) = yy;
        elsqsum = elsqsum + (D(i,2)-yy)^2;
    end
    elsqsum = sqrt(elsqsum);
    plot(xx,vv);
    hold on;
    fprintf('  imeth= %i n=%i m=%i l(errlsq) = %e \n',imeth,n,m,elsqsum)
    s = lastwarn;
end
