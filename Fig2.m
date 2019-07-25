% first case
dx=10;
f=0.1;
L=1000;
S0 = SeedShadow(5,3000,20,L,1);
S = BlockMatrix(S0,dx);
use=randperm((L/dx)^2,100);
[x0,y0]=ind2sub([L/dx L/dx],use);
[x1,y1]=ind2sub([L/dx L/dx],use+2);
x = S(use);
y = binornd(x,f);
w = S(use+2);
p = glmfit(log(w),y,'poisson');

clf
subplot(321)
pcolor(S0);shading flat
hold all
contour(smooth(S0,5),[1 2 4 8],'color','w');
% clabel(C,h)
h = colorbar;
ylabel(h, 'Seed density')

for i=1:100
    rectangle('Position',[x0(i)*dx y0(i)*dx dx dx],'EdgeColor','w','linewidth',1)
    rectangle('Position',[x1(i)*dx y1(i)*dx dx dx],'EdgeColor','r','linewidth',1)
end
daspect([1 1 1])
axis([1 250 1 250])

subplot(234)
plot(x,y,'.','markersize',15)
hold all
plot(w,y,'.','markersize',15)

xi=0:200;
plot(xi,f*xi,'b-','linewidth',2)
plot(xi,exp(p(1))*xi.^p(2),'r-','linewidth',2)
xlabel('number of seeds')
ylabel('number of seedlings')
legend('true','observed','location','NorthWest')
legend('boxoff')



%% second case
f = 3;dx=50;dispersal=10;L=1000;
[S0,A0,xn,yn,x,y] = SeedShadow(20,f,dispersal,L,1);
S1 = BlockMatrix(S0,dx);
A1 = smooth(A0,dispersal,'fun','fat','dx',dx);
A2 = BlockMatrix(A0,dx);
mean(A2(:)==0 & S1(:)>0)

subplot(232)
hold off;
plot(xn,yn,'.','markersize',5);
hold all
plot(x,y,'.','markersize',10);
axis([0 250 0 250])
daspect([1 1 1])
set(gca,'xtick',0:dx:L,'ytick',0:dx:L);
grid on
legend('Juveniles','Adults')

use = A2(:)>0;
p = glmfit(log(A2(use)),S1(use),'Poisson');

subplot(235)
hold off
plot(A1(:),S1(:),'.','markersize',15)
hold all
plot(A2(:),S1(:),'.','markersize',15)

xi=0:.1:15.5;
plot(xi,f*xi,'b-','linewidth',2)
plot(xi,exp(p(1))*xi.^p(2),'r-','linewidth',2)
xlabel('number of adults')
ylabel('number of saplings')
legend('true','observed','location','NorthWest')
legend('boxoff')


%% third case
%load allometry
dat = importdata('Data\DiameterHeightCrown20190104.csv');

D = dat.data(:,4);
C = dat.data(:,6);

use = D>0 & C>0;

[res,gof,stats] = fit(log(D(use)),log(C(use)),'poly1');

su = gof.rmse;

figure
subplot(4,3,3)
loglog(D,C,'.');hold all
xi = logspace(-1,3,1000);
loglog(xi,exp(res.p2)*xi.^res.p1,'r-')
ylabel('Crown area (m^2)')
xlabel('Diameter (cm)')
set(gca,'xticklabel',[])

subplot(4,3,6)
histogram(log10(C) - (res.p2+res.p1*log(D))/log(10),40)
xlabel('log_{10}(residual)')
set(gca,'yticklabel',[],'box','off','xtick',[-1 0 1],'xticklabel',{'10^-^1','1','10^1'},'xlim',[-1.5 1.5])
%% MAIN LOOP
%load BCI data
load('Data\bci7.mat');

% compute abundance and max DBH
species = unique(sp);
S = length(species);
N0 = zeros(S,1);
Dmax = nan(S,1);
for i=1:S
    use = strcmp(sp,species(i))&dbh>0;
    if sum(use)>0
    N0(i) = sum(use);
    Dmax(i) = max(dbh(use));
    end
end


% select species with N>250 and max DBH>200
I=find(N0>250 & Dmax>200);
M =length(I);
S = species(I);
BA=zeros(M,1);
N=zeros(M,1);
n=zeros(M+1,1);
Lx = 1000;
Ly = 500;
scale = 5;
for i=1:M
    use = strcmp(sp,species(I(i)))&dbh>0;
    BA(i)=pi/4*sum((dbh(use)*1e-3).^2)/50e+4;
    N(i)=sum(use);
    
     use =strcmp(sp,species(I(i))) & dbh>0 & dbh<=100 & ...
        gx>=2*scale & gx<=Lx-2*scale & ...
        gy>=2*scale & gy<=Ly-2*scale;
    
    n(i+1)=sum(use);
end

% parameters&constants
X=zeros(sum(n),1);
W1=zeros(sum(n),1);
W2=zeros(sum(n),1);
basal=zeros(sum(n),1);
group=zeros(sum(n),1);
for i=1:M
    i
    SP = species(I(i));
    
    % select focal saplings, exclude a buffer zone 2x scale
    c1 = find(strcmp(sp,SP) & dbh>0 & dbh<=100 & ...
        gx>=2*scale & gx<=Lx-2*scale & ...
        gy>=2*scale & gy<=Ly-2*scale);
    
    n1 = length(c1);
    % select conspecifics including the buffer zone
    c0 = strcmp(sp,SP) & dbh>0 & gx>=0;
    
    u = -1/2*su^2+su*randn(sum(c0),1);       % random error
    D2 = exp(res.p2)*dbh(c0).^res.p1;        % allometric crown area
    H2 = D2.*exp(u);                         % true crown area
    Y2 = pi/4*dbh(c0).^2;             % bias allometric crown area
    % loop within focal individuals to estimate local neghborhhod
    
    x=zeros(n1,1);
    w1=zeros(n1,1);
    w2=zeros(n1,1);
    parfor j=1:n1
        r = sqrt((gx(c0)-gx(c1(j))).^2 + (gy(c0)-gy(c1(j))).^2); %distance to focal sapling
        use=r>0 & r<scale;
        x(j)  = sum(H2(use));  %true neighborood density
        w1(j) = sum(D2(use));  %observed neighborood density
        w2(j) = sum(Y2(use));  %observed neighborood density
    end
    
    index = sum(n(1:i))+1:sum(n(1:i+1));
    X(index)  = x;
    W1(index) = w1;
    W2(index) = w2;
    basal(index) = log10(BA(i)*1e+8);
    group(index)=i;
    
end

save('Cartoon3.mat','X','W1','W2','basal','group','BA','M','n')
%% figure
% load('Cartoon3.mat')
%logistic suvival probability
X0=X*1e-4; %change to ha
b=-.5;
p0=0.8; % survival probability in absence of conspecifics
a0=log(p0/(1-p0));
p1 = (1+exp(-a0-b*X0)).^-1;

%bounded logistic suvival probability
p0=0.5;
a0=log(p0/(1-p0));
p2 = (1+exp(-a0-b*X0)).^-1 + 0.3;

% plot(X,[p1 p2],'.')
% pause

beta1 = zeros(M,2);
beta2 = zeros(M,2);
E = zeros(M,1);
W01=W1*1e-4;
W02=W2*1e-4;
for i=1:M
    i
    index = sum(n(1:i))+1:sum(n(1:i+1));
    beta1(i,:) = glmfit(W01(index),p2(index),'Binomial');
    beta2(i,:) = glmfit(W02(index),p1(index),'Binomial');
    E(i)=var(W02(index)-X0(index))/var(W02(index));
end

% subplot(2,3,6)
figure(3);clf
subplot(211)
plot(log10(BA),beta2(:,2),'.','markersize',15)
set(gca,'xticklabel',[])
ylabel('stronger  <--- CNDD')
lsline

subplot(413)
plot(log10(BA),log(E),'.','markersize',15)
lsline
xlabel('log species abundance (basal area)')
ylabel('error')
subplot(414)
plot(X0,[p1 p2],'.')
% h=refline(0,b);
% set(h,'linestyle','-','linewidth',2)
% legend('fitted','true','Location','NorthWest')
% legend('boxoff')

