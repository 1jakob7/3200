% Jakob Horvath, u1092049
% Modifies the supplied Gatlin_image.m script to work with MATLAB's
% "mandrill" and "durer" images. Uses SVD factorization to plot
% approximations for both images with ranks of 2, 4, 8, 16, 32, 64, and
% 128. 
% * Warning -- produces a lot of figures! *

load mandrill
% load the "mandrill" image data, built-in to MATLAB
[U,S,V] = svd(X);
% "mandrill" stores the image as the variable "X"
figure(1),clf
% plot the singular values
semilogy(diag(S),'b.','markersize',20)
set(gca,'fontsize',16)
title('singular values of the "mandrill" image matrix')
xlabel('k'), ylabel('\sigma_k')
figure(2),clf
% plot the original image
image(X), colormap(gray)
% image: MATLAB command to display a matrix as image
axis equal, axis off
title('true image (rank 480)','fontsize',16)
% plot the rank-k approximations
k = [2,4,8,16,32,64,128];
for i=1:length(k)
    figure(i+2),clf
    Xk = U(:,1:k(i))*S(1:k(i),1:k(i))*V(:,1:k(i))';
    image(Xk), colormap(gray)
    axis equal, axis off
    title(sprintf('best rank-%d approximation',k(i)),'fontsize',16)
end

load durer
% load the "durer" image data, built-in to MATLAB
[U2,S2,V2] = svd(X);
% "durer" stores the image as the variable "X"
figure(10),clf
% plot the singular values
semilogy(diag(S2),'b.','markersize',20)
set(gca,'fontsize',16)
title('singular values of the "durer" image matrix')
xlabel('k'), ylabel('\sigma_k')
figure(11),clf
% plot the original image
image(X), colormap(gray)
% image: MATLAB command to display a matrix as image
axis equal, axis off
title('true image (rank 480)','fontsize',16)
% plot the rank-k approximations
k = [2,4,8,16,32,64,128];
for i=1:length(k)
    figure(i+11),clf
    Xk = U2(:,1:k(i))*S2(1:k(i),1:k(i))*V2(:,1:k(i))';
    image(Xk), colormap(gray)
    axis equal, axis off
    title(sprintf('best rank-%d approximation',k(i)),'fontsize',16)
end
