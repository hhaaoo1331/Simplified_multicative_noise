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

beta1 = [100 50 10 1 0.1 0.01 0.001];
beta2 = [100 50 10 1 0.1 0.01 0.001];
alpha1 = [100 50 10 1 0.1 0.01 0.001];
opts.maxitr = 5000;
tic
for i=1:length(beta1)
    for j=1:length(beta2)
        for k=1:length(alpha1)
            opts.beta1 = 0.05;
            opts.beta2 = 0.0003;
            opts.alpha1 = 0.004;
            out3=CMC(A,Bn,I,N,mvw,opts);
        end
    end
end
toc
fprintf('psnr: %2.4f, iteration: %4f, time: %4.4f\n',out3.psnr(end),out3.itr,t3);

out3.sol(out3.sol>255)=255;
sol = out3.sol;
figure,imshow(sol,[])
figure,plot(out3.psnr);
% 
% % t1=cputime;out1= AA(Bn,I,alpha,3000);t1=cputime-t1;
% % fprintf('psnr: %2.4f, iteration: %4f, time: %4.4f\n',out1.psnr(end),out1.itr,t1);
% out1.sol(out1.sol>255)=255;figure,imshow(out1.sol,[])

%  
%  
% opts.nitr=3;opts.beta=100;t2=cputime;out2= ADM(Bn,I,0.5,opts);t2=cputime-t2;
% % fprintf('psnr: %2.4f, iteration: %4f, time: %4.4f\n',out2.psnr(end),out2.itr,t2);
% out2.sol(out2.sol>255)=255;figure,imshow(out2.sol,[])
%          end
%      end
%  end



