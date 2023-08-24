% 
clear all
close all
I=imread('cameraman.tif');
I=double(I);
[m,n]=size(I); 
% 
% 
% % Add noise  with the mean 1
%%Gamma noise
L=10;%%Look
N= gamrnd(L,1/L,m, n);
% % Gaussian noise
% N=1+0.3*randn(m,n);
% % Rayleigh noise
%N=rand(m,n);
%%sigma=1
%%N=5*sqrt(-2*log(1-N));
%%caculate the mean of w
mvw=sum(sum(1./N))/m/n;
% 
A = fspecial('gaussian',7,2);
Bn = imfilter(I,A,'circular','conv');
Bn = Bn.*N;
figure,imshow(I,[]);
figure,imshow(Bn,[]);

opts.beta1 = 0.1;
opts.beta2 = 0.1;
opts.alpha1 = 0.1;
opts.maxitr = 1000;
opts.beta1 = 0.05;
opts.beta2 = 0.0003;
opts.alpha1 = 0.004;

out=CMC(A,Bn,I,N,mvw,opts);

fprintf('psnr: %2.4f, iteration: %4f, time: %4.4f\n',out.psnr(end),out.itr);
out.sol(out.sol>255)=255;
sol = out.sol;
figure,imshow(sol,[])
figure,plot(out.psnr);