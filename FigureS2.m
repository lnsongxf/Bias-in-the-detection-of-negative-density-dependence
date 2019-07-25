a=1;
b=1;
N=200;
se=0.3;
lambda = 3/4; % s2x/(s2x+s2u)

figure(1);clf
%% classic error  (OLS)
su=sqrt((1-lambda)/lambda);
x  = randn(N,1);
y = a + b*x + se*randn(N,1);
w = x + su*randn(N,1);

p1 = polyfit(w,y,1);

subplot(231)
plot(x,y,'ko','markersize',4,'MarkerFaceColor','k');hold all
plot(w,y,'ro','markersize',4,'MarkerFaceColor','r')
plot([-4 4],b*[-4 4]+a,'k-','linewidth',2)
plot([-4 4],p1(1)*[-4 4]+p1(2),'r-','linewidth',2)
axis([- 4 4 -4+a  4+a])
axis square
title('Classic error')
legend('real','observed with error','location','northwest')
legend('boxoff')
ylabel({'OLS','y'})
xlabel('x')

%% Berkson error (OLS)
su=sqrt((1-lambda)/lambda);
w  = sqrt(1-su.^2)*randn(N,1);
x = w + su*randn(N,1);
y = a + b*x + se*randn(N,1);

p1 = polyfit(w,y,1);

subplot(232)
plot(x,y,'ko','markersize',4,'MarkerFaceColor','k');hold all
plot(w,y,'ro','markersize',4,'MarkerFaceColor','r')
plot([-4 4],b*[-4 4]+a,'k-','linewidth',2)
plot([-4 4],p1(1)*[-4 4]+p1(2),'r-','linewidth',2)
axis([- 4 4 -4+a  4+1])
axis square
title('Berkson error')
ylabel('y')
xlabel('x')

%% differential error  (OLS)
su=sqrt((1-lambda)/lambda);
x  = randn(N,1);
y = a + b*x + se*randn(N,1);
R = 0.5*[1 -0.5;-0.5 1];
C = chol(R);
v = randn(N,2)*C;
k = y + v(:,1);
w = x + v(:,2);

p1 = polyfit(w,k,1);

subplot(233)
plot(x,y,'ko','markersize',4,'MarkerFaceColor','k');hold all
plot(w,k,'ro','markersize',4,'MarkerFaceColor','r')
plot([-4 4],b*[-4 4]+a,'k-','linewidth',2)
plot([-4 4],p1(1)*[-4 4]+p1(2),'r-','linewidth',2)
axis([- 4 4 -4+a  4+a])
axis square
title('Differential error')
ylabel('y')
xlabel('x')
%% classic error (GLM)
a=1;u=rand(N,1);p=0.1;
x0 = poissrnd(3,N,1);
x = x0+u;
y  = poissrnd(a*x);

h1 = binornd(x0,p);
h2 = binornd(x0,p);
w = x0-h1+h2+u;

b1 = glmfit(log(w),y,'poisson','offset',log(w));

subplot(234)
xi=linspace(0.001,1.1*max(w),200);
plot(x,y,'ko','markersize',4,'MarkerFaceColor','k');hold all
plot(w,y,'ro','markersize',4,'MarkerFaceColor','r');hold all
plot(xi,a*xi,'k-','linewidth',2)
plot(xi,exp(b1(1))*xi.^(b1(2)+1),'r-','linewidth',2)
axis([0 xi(end) 0 1.1*max(y)])
axis square
ylabel({'GLM','N_1'})
xlabel('x')

%% Berkson error (GLM)
su=0.5;N=2000;
u=-1.2*su^2+su*randn(N,1);
b0=[2;-.1];

% w =exprnd(1,N,1);
% x=w.*exp(u);
% p=glmval(b0,x,'logit');
N0 = poissrnd(20,N,1);
x = N0.*exp(u);
p=glmval(b0,x,'logit');

N1 = binornd(N0,p);
use=N0>0;
b1 = glmfit(N0(use),[N1(use) N0(use)],'binomial');
p1 = N1./N0;

% subplot(235)
clf
xi=linspace(0,1.1*max(x),200);
plot(x(use),p1(use),'ko','markersize',4,'MarkerFaceColor','k');hold all
plot(N0(use),p1(use),'ro','markersize',4,'MarkerFaceColor','r');hold all
plot(xi,glmval(b0,xi,'logit'),'k-','linewidth',2)
plot(xi,glmval(b1,xi,'logit'),'r-','linewidth',2)
axis([0 xi(end) 0 1])
axis square
ylabel('suvival probability')
xlabel('x')

%% differential error  (GLM)
a=1;u=rand(N,1);p=0.1;
x0 = poissrnd(3,N,1);
x = x0+u;
y  = poissrnd(a*x);

hx = binornd(x0,p);
hy = binornd(y ,p);

w = x0-hx+hy+u;
k = y-hy+hx;

b1 = glmfit(log(w),k,'poisson','offset',log(w));

subplot(236)
xi=linspace(0.001,1.1*max(w),200);
plot(x,y,'ko','markersize',4,'markerfacecolor','k');hold all
plot(w,k,'ro','markersize',4,'markerfacecolor','r');hold all
plot(xi,a*xi,'k-','linewidth',2)
plot(xi,exp(b1(1))*xi.^(b1(2)+1),'r-','linewidth',2)
axis([0 xi(end) 0 1.1*max(y)])
axis square
ylabel('N_1')
xlabel('x')

sublabel
