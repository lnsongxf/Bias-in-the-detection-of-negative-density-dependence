%% fit survival model for BCI
clear
load('SurvivalAnalysis.mat')
finfo = dir('SaplingSurvSP_*.mat');
M = length(finfo);
c=logspace(-2,0,100);
FitMethod='Laplace';

inter=zeros(M,length(c));
cndd=zeros(M,length(c));
hndd=zeros(M,length(c));
sei=zeros(M,length(c));
sec=zeros(M,length(c));
seh=zeros(M,length(c));
LogL=zeros(length(c),2);
for k=1:length(c)
    L1=0;L2=0;P=zeros(M,2);
    parfor i=1:M
        [k i]
        dat = load(finfo(i).name);
        tbl=table(dat.surv,dat.Z1.^c(k),dat.Z2.^c(k),dat.census,...
            'VariableNames',{'y','z1','z2','c'});
        N = length(dat.Z1);
        P0 = mean(tbl.y);
        glm = fitglme(tbl,'y ~ 1 + z1 + z2 + (1|c)','distribution','Binomial',...
            'FitMethod',FitMethod,'Mustart',P0*ones(N,1));

        [B0,~,stats] = fixedEffects(glm);
        inter(i,k)=B0(1);
        cndd(i,k)=B0(2);
        hndd(i,k)=B0(3);
        
        sei(i,k)=stats.SE(1);
        sec(i,k)=stats.SE(2);
        seh(i,k)=stats.SE(3);
        
        p = predict(glm,tbl);
        L1 = L1 + sum(log(p(tbl.y==1)))+sum(log(1-p(tbl.y==0)));
        L2 = L2 + glm.LogLikelihood;
    end
    LogL(k,1)=L1;
    LogL(k,2)=L2;
end

%% figure
figure(1);clf

subplot(311)
errorbar(log10(BAavg),cndd(:,end),sec(:,end),'.');hold all
errorbar(log10(BAavg),hndd(:,end),seh(:,end),'.');hold all
trend = fit(log10(BAavg),cndd(:,end),'exp1','weights',1./sec(:,end));
plot(trend,'-b')
ylabel('{\itb}_{GLMM}')
xlabel([])
xlim([-6.3 -3.5])
title('\eta = {\itb}_0 + {\itb}_1 {\itZ}_1 + {\itb}_2 {\itZ}_2')
legend(['CNDD (' char(177) '1 SE)'],['HNDD (' char(177) '1 SE)']);legend('boxoff')
xlabel('log(abundance)')


subplot(312)
semilogx(c,LogL)
xlabel('{\itc}')
ylabel('Log Likelyhood')
hold all
[~,cmax]=max(LogL(:,1));
plot(c(b),LogL(b,1),'ro')
plot(c(end),LogL(end,1),'ro')


subplot(313)
[~,cmax]=max(LogL(:,1));
errorbar(log10(BAavg),cndd(:,cmax),sec(:,cmax),'.');hold all
errorbar(log10(BAavg),hndd(:,cmax),seh(:,cmax),'.');hold all
trend = fit(log10(BAavg),cndd(:,cmax),'poly1','weights',1./sec(:,cmax));
plot(trend,'-b')
ylabel('{\itb}_{GLMM}')
xlabel([])
xlim([-6.3 -3.5])
title('\eta = {\itb}_0 + {\itb}_1 {\itZ}_1^{ 0.22} + {\itb}_2 {\itZ}_2^{ 0.22}')
xlabel('log(abundance)')


%% k-fold cross-validation
clear
load('SurvivalAnalysis.mat')
finfo = dir('SaplingSurvSP_*.mat');
M = length(finfo);
c=[1 .22];
nmodel = length(c);
K=10;
LogL=zeros(M,nmodel);
RMSE=zeros(M,nmodel);

parfor i=1:M
    i
    dat = load(finfo(i).name);
    N = length(dat.Z1);
    n=floor(N/K);
    perm=randperm(N);
    for k=1:nmodel
        tbl=table(dat.surv(perm),dat.Z1(perm).^c(k),dat.Z2(perm).^c(k),dat.census(perm),...
            'VariableNames',{'y','z1','z2','c'});
        SSE = 0;
        if n*K<N; k0=K+1; else; k0=K; end
        for h=1:k0 %k-fold loop
            use = false(N,1);
            use((h-1)*n+1:min(h*n,N)) = true;
            glm = fitglme(tbl,'y ~ 1 + z1 + z2 + (1|c)','distribution','Binomial',...
                'FitMethod','MPL','Exclude',use);
            p = predict(glm,tbl);
            LogL(i,k) = LogL(i,k) + sum(log(p(tbl.y(use)==1)))+sum(log(1-p(tbl.y(use)==0)));
            SSE = SSE + sum((p(use)-tbl.y(use)).^2)
        end
        RMSE(i,k) = sqrt(SSE);
    end
end

%% plot cross-validatio
figure(1)

subplot(211)
plot(log(BAavg),RMSE(:,1)-RMSE(:,2),'.')
refline(0,0)
ylabel('RMSE({\itc}=1) - RMSE({\itc}=0.22)')
load('SurvivalAnalysis.mat')

subplot(212)
plot(log(BAavg),-LogL(:,1)+LogL(:,2),'.')
refline(0,0)
ylabel('-LogL({\itc}=1) + LogL({\itc}=0.22)')
xlabel('log (abundance)')
