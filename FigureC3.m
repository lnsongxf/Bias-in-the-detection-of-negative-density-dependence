%% MAIN LOOP
%load BCI data census 7
load('bci7.mat');

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
X1 = zeros(sum(n),1);
X2 = zeros(sum(n),1);
W  = zeros(sum(n),1);
basal=zeros(sum(n),1);
su1=0.7;
su2=1.4;
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
    
    u1 = -1/2*su1^2+su1*randn(sum(c0),1);       % random error
    u2 = -1/2*su2^2+su2*randn(sum(c0),1);       % random error

    D2 = 0.6*1e-4*dbh(c0).^2;                      % allometric equation (D^2)
    H1 = D2.*exp(u1);                         % true biomass
    H2 = D2.*exp(u2);                         % true biomass
    
    % loop within focal individuals to estimate local neghborhhod
    x1=zeros(n1,1);
    x2=zeros(n1,1);
    w=zeros(n1,1);
    parfor j=1:n1
        r = sqrt((gx(c0)-gx(c1(j))).^2 + (gy(c0)-gy(c1(j))).^2); %distance to focal sapling
        use=r>0 & r<scale;
        x1(j) = sum(H1(use));  %true neighborood density
        x2(j) = sum(H2(use));  %true neighborood density
        w(j)  = sum(D2(use));  %observed neighborood density (Basal allometry)
    end
    
    index = sum(n(1:i))+1:sum(n(1:i+1));
    X1(index) = x1;
    X2(index) = x2;
    W(index)  = w;
    basal(index) = log10(BA(i)*1e+8);
    
end

%% figure
%logistic suvival probability

b=-[0.05 0.1 0.5];
p0=0.8; % survival probability in absence of conspecifics
a0=log(p0/(1-p0));


beta1 = zeros(M,4);
beta2 = zeros(M,3);
for j=1:3
    p1 = (1+exp(-a0-b(j)*X1)).^-1;
    p2 = (1+exp(-a0-b(j)*X2)).^-1;
    for i=1:M
        index = sum(n(1:i))+1:sum(n(1:i+1));
        bhat = glmfit(W(index),p1(index),'Binomial');
        beta1(i,j)=bhat(2);
        bhat = glmfit(W(index),p2(index),'Binomial');
        beta2(i,j)=bhat(2);
        
    end
end


figure(1);clf
subplot(311)
xi=linspace(0,max(X2),1000);
for j=1:3;plot(xi,(1+exp(-a0-b(j)*xi)).^-1,'linewidth',2);hold all;end
legend(num2str(b'))
xlim([0 max(X2)])
xlabel('conspecific density')
ylabel('survival probability')

for i=1:3
    subplot(3,3,3+i)
    plot(log10(BA),beta1(:,i),'.','markersize',15)
    title(['{\it b} = ' num2str(b(i))])
    if i==1;ylabel('{\itb}_{GLM}');end
    lsline
    
    subplot(3,3,6+i)
    plot(log10(BA),beta2(:,i),'.','markersize',15)
    if i==2;xlabel('log species abundance (basal area)');end
    if i==1;ylabel('{\itb}_{GLM}');end
    lsline
end
