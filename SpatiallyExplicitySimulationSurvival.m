% Simulate spatially explity seed shadow and survival probability
% and fit model of CNDD
%%
clear
L=500;
[X,Y]=meshgrid(1:L,1:L);
[m,n]=size(X);
scale = 5;
npuls_2   = floor((n-1)/2);
pulsx     = 2*pi/n*[ 0:npuls_2  (npuls_2-n+1):-1 ];
npuls_2   = floor((m-1)/2);
pulsy     = 2*pi/m*[ 0:npuls_2  (npuls_2-m+1):-1 ];
[kx,ky] = meshgrid(pulsx,pulsy);
K2=kx.^2+ky.^2;
kernel = (1+K2*scale.^2).^(-3/2);


N0=4000*1e-4;    
D0=20;
p0=0.5;
a=log(p0/(1-p0));
b=1;a=10;
mu0 = N0*D0;
sd0 = sqrt((2*N0.*D0^2)/(8*pi*scale^2));

c=[1 0.95 0.9];
su=0.;
M=20;
N=logspace(log10(30*1e-4),log10(N0/3),M);

cl1=nan(M,3,3);

v1=nan(M,1);
v2=nan(M,1);
q1=nan(M,1);
q2=nan(M,1);
mu1=nan(M,1);
mu2=nan(M,1);
parfor i=1:M
    i
    
    [S,A] = SeedShadow(N(i)*1e+4,100,10,L,1);
    use=A(:)>0;A(use)=gamrnd(A(use),D0); %sum of exponentials is a gamma
    f=fft2(A);W1=ifft2(f.*kernel);
    
    B = poissrnd(N0-N(i),L,L);
    use=B(:)>0;B(use)=gamrnd(B(use),D0);
    f=fft2(B);W2=ifft2(f.*kernel);

   
    use=S>0;
    n0=S(use);
      
    for j=1:3
    p = (1+exp(-a+b*(W1(use).^c(j)+W2(use).^c(j)))).^-1;
    y = binornd(n0,p);
    cl1(i,:,j) = glmfit([W1(use) W2(use)],[y n0],'binomial');
    end

end

%% plotting
red =[0.8500    0.3250    0.0980];
blue = [0    0.4470    0.7410];
ax = [-2.3 -0.4 -b-.5 -b+.5]; 

M1 = N*D0;
V1 = (2*N).*D0^2/(8*pi*scale^2);
M2 = (N0-N)*D0;
V2 = (2*(N0-N)).*D0^2/(8*pi*scale^2);

i=M;
[S,A] = SeedShadow(N(i)*1e+4,100,10,L,1);
use=A(:)>0;A(use)=gamrnd(A(use),D0); %sum of exponentials is a gamma
f=fft2(A);W1=ifft2(f.*kernel);

B = poissrnd(N0-N(i),L,L);
use=B(:)>0;B(use)=gamrnd(B(use),D0);
f=fft2(B);W2=ifft2(f.*kernel);
    
figure(1);clf

for i=1:3
subplot(2,3,i)
plot(W2(use),W2(use).^c(i),'.','color',red,'markersize',2);hold all
plot(W1(use),W1(use).^c(i),'.','color',blue,'markersize',2)
ylabel('{\itZ}_{true}')
xlabel('{\itZ}_{obs}')
axis([0 ceil(max(W2(use))) 0 ceil(max(W2(use)))]);
axis square
title(['{\itZ}_{true} = {\itZ}_{obs}^{ ' num2str(c(i)) '}'])
refline(1,0)

subplot(2,3,3+i)
plot(log10(N/N0),squeeze(cl1(:,2:3,i)),'o','markersize',4);hold all
plot(log10(N/N0),-b.*c(i)*M1.^(c(i)-1),'-','color',blue)
plot(log10(N/N0),-b.*c(i)*M2.^(c(i)-1),'-','color',red)
xlabel('relative abundance')
ylabel('{\itb}_{GLM}')
axis(ax)
refline(0,-b)
set(gca,'xtick',log10([0.01:0.01:0.1 0.2:0.1:0.4]),'xticklabel',...
    {'10^-^2','','','','','','','','','10^-^1','','','',''})
axis square
pause(.1)
end
legend('CNDD','HNDD','Location','SouthEast');legend('boxoff')
sublabel;

