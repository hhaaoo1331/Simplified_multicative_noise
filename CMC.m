function [out] = CMC(H,Bn,I,N,mvw,opts)

beta1 = opts.beta1;
beta2 = opts.beta2;
alpha1 = opts.alpha1;
[m n] = size(Bn);%G
[D,Dt]= defDDt;
eigsH = psf2otf(H,[m,n]);
eigsHTH = abs(eigsH).^2;
eigsDTD=abs(psf2otf([1,-1],[m n])).^2 + abs(psf2otf([1;-1],[m n])).^2;
% initi
Lam1 = zeros(2*m,n);
Lam2 = zeros(m,n);
% initi x=[f; s]
f=Bn;
w=zeros(m,n);
[D1X,D2X] = D(Bn);
%  history
out.ssim = [];
out.psnr   = [];
out.relchg = [];
out.relf=[];
out.reln=[];
out.relchgw=[];
out.fval=[];
fvalchg=1;
out.fvalchg=[];
relchg=1;
fval=1;
ii=0;
while relchg> 5*10^(-4) && ii<opts.maxitr
    %%%% Step 1 update y=[p;z;w]
    % ==================
    %  p subproblem
    % ==================
    Z1 = D1X + Lam1(1:m,:)/beta2;
    Z2 = D2X + Lam1(m+1:2*m,:)/beta2;
    V = Z1.^2 + Z2.^2;
    V = sqrt(V);
    V(V==0) = 1;
    V = max(V - alpha1/beta2, 0)./V;
    p(1:m,:) = Z1.*V;
    p(m+1:2*m,:) = Z2.*V;
    % ==================
    %  w subproblem
    % ==================
    wp=w;
    Hf = imfilter(f,H,'circular','conv');
    right=mvw*ones(m,n)+Bn.*(beta1*Hf+Lam2);
    w=right./(ones(m,n)+beta1*(Bn.^2));
    %%%% Step 2 update [f; s]
    % ==================
    %  f subproblem
    % ==================
%     fp=f;
%     b1_1 = Dt(p(1:m,:),p(m+1:2*m,:));
%     b1_2 = -Dt(Lam1(1:m,:),Lam1(m+1:2*m,:))/beta1;
%     b1_3 = imfilter(Bn.*w,H,"circular","corr");
%     b1_4 = -imfilter(Lam2/beta1,H,"circular","corr");
%     b1 = fft2(b1_1 + b1_2 + b1_3 + b1_4);
%     b2=(eigsHTH+eigsDTD);
%     f = b1./b2;
%     f = real(ifft2(f));
%     [D1X,D2X] = D(f);
    fp=f;
    b1_1 = beta2*Dt(p(1:m,:),p(m+1:2*m,:));
    b1_2 = -Dt(Lam1(1:m,:),Lam1(m+1:2*m,:));
    b1_3 = beta1*imfilter(Bn.*w,H,'circular','corr');
    b1_4 = -imfilter(Lam2,H,'circular','corr');
    b1 = fft2(b1_1 + b1_2 + b1_3 + b1_4);
    b2=(beta1*eigsHTH+beta2*eigsDTD);
    f = b1./b2;
    f = real(ifft2(f));
    [D1X,D2X] = D(f);
    %%%%metric
    psnrf=psnr1(f,I);
    ssimf = SSIMour(I,f);
    out.ssim = [out.ssim,ssimf];
    out.psnr = [out.psnr;  psnrf ];
    relchg=norm(f - fp,'fro')/norm(f,'fro');
    out.relchg = [out.relchg;relchg];
    out.relchgw=[out.relchgw,norm(w - wp,'fro')/norm(w,'fro')];
    relf=norm(f- I,'fro')/norm(I,'fro'); %relative error wrt true image
    reln=norm(1./w-N,'fro')/norm(N,'fro');%relative error wrt true noise
    out.relf=[out.relf,relf];
    out.reln=[out.reln,reln];
    fval_old=fval;
    S = Bn.*w;
    fval=0.5*(norm(w-mvw*ones(m,n),'fro')).^2+alpha1*sum(sum(sqrt(D1X.^2 + D2X.^2)));
    fvalchg=abs(fval-fval_old)/abs(fval);
    out.fvalchg=[out.fvalchg,fvalchg];
    out.fval=[out.fval,fval];
    ii=ii+1;
    % ==================
    %  Update Lam
    % ==================
    Hf = imfilter(f,H,'circular','conv');
    Lam1(1:m,:) = Lam1(1:m,:) + beta2*(D1X-p(1:m,:));
    Lam1(m+1:2*m,:) = Lam1(m+1:2*m,:) + beta2*(D2X-p(m+1:2*m,:));
    Lam2 = Lam2 + beta1*(Hf-S);
end
out.sol = f;
% out.w=w;
% out.S=S;
%out.z=z;
out.itr = ii;
% fprintf('Iter: %d, psnrf: %4.2f, relf: %4.2f,reln: %4.2f\n',ii,psnrf,relf,reln);
% mv=sum(w(:))/m/n;
% var=sum((w(:)-mv).^2)/m/n;
% varw=sum(sum((1./N-mvw).^2))/m/n;
% fprintf('mv: %4.2f, var: %4.2f, mvw: %4.2f, varw: %4.2f\n',mv,var,mvw,varw);
