function p = psnr1(x,y)

m1 = max(abs(y(:)))^2;
m2 = sum(abs((y(:)-x(:))).^2);
[m,n] = size(x);
p = 10*log10((m*n*m1)/m2);