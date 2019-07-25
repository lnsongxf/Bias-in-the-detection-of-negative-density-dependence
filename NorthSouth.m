clear
DX=20;
%% BCI: load data
load('Data\bci7.mat')
use=dbh>0;
gx=gx(use);gy=gy(use);sp=sp(use);dbh=dbh(use);
species=unique(sp);
S=length(species);
Lx=1000;Ly=500;bci.Area=50;
%% BCI: analyze
bci.sl=nan(S,1);
bci.slr=nan(S,1);
bci.sl0=nan(S,1);
bci.BA=zeros(S,1);
bci.N=zeros(S,1);
bci.mu=nan(S,1);
bci.s2=nan(S,1);

for i=1:S
    i
    use = strcmp(sp,species(i));
    bci.BA(i)=pi/4*sum((dbh(use)/1000).^2)/bci.Area;
    bci.N(i)=sum(use);
    M=median(dbh(use));
    if bci.N(i)>25
        use1=use&dbh>0&dbh<M;
        use2=use&dbh>0&dbh>=M;
        S1=QuadratCount(gx(use1),gy(use1),[Lx Ly],DX);%Saplings
        A1=QuadratCount(gx(use2),gy(use2),[Lx Ly],DX);%Adults  
        p=polyfit(log(A1(:)+1),log(S1(:)+1),1);
        bci.sl(i)=p(1);
        
        map = QuadratCount(gx(use),gy(use),[Lx Ly],[DX DX/2]); 
        S1 = map(1:2:2*Ly/DX,:);
        A1 = map(2:2:2*Ly/DX,:);%North   
        p=polyfit(log(A1(:)+1),log(S1(:)+1),1);
        bci.slr(i)=p(1);
        
        %theoretical slope
        S0=QuadratCount(gx(use),gy(use),[Lx Ly],DX);%all
        bci.mu(i)=mean(S0(:));
        bci.s2(i)=var(S0(:));
        bci.sl0(i)=2./(bci.mu(i)./bci.s2(i)+1) - 1;
        
    end
end
use1=bci.BA>1e-4&bci.N>25;
[bci.res,bci.gof]=fit(bci.sl0(use1),bci.sl(use1),'poly1');

%% SERC: load data
load('Data\serc2.mat')
use=dbh>0;
gx=gx(use);gy=gy(use);sp=sp(use);dbh=dbh(use);
species=unique(sp);
S=length(species);
Lx=400;Ly=400;serc.Area=16;
%% SERC: analysis
serc.sl=nan(S,1);
serc.slr=nan(S,1);
serc.sl0=nan(S,1);
serc.mu=nan(S,1);
serc.s2=nan(S,1);
serc.BA=nan(S,1);
serc.N=zeros(S,1);

for i=1:S
    i
    use = strcmp(sp,species(i));
    serc.BA(i)=pi/4*sum((dbh(use)/1000).^2)/serc.Area;
    serc.N(i)=sum(use);
    M=median(dbh(use));
    if serc.N(i)>11
        use1=use&dbh<M;
        use2=use&dbh>=M;
        S1=QuadratCount(gx(use1),gy(use1),[Lx Ly],DX);
        A1=QuadratCount(gx(use2),gy(use2),[Lx Ly],DX);       
        p=polyfit(log(A1(:)+1),log(S1(:)+1),1);
        serc.sl(i)=p(1);
        
        map = QuadratCount(gx(use),gy(use),[Lx Ly],[DX DX/2]); 
        S1 = map(1:2:2*Ly/DX,:);
        A1 = map(2:2:2*Ly/DX,:);%Adults   
        p=polyfit(log(A1(:)+1),log(S1(:)+1),1);
        serc.slr(i)=p(1);
        
        %theoretical slope
        S0=QuadratCount(gx(use),gy(use),[Lx Ly],DX);%all
        serc.mu(i)=mean(S0(:));
        serc.s2(i)=var(S0(:));
        serc.sl0(i)=2./(serc.mu(i)./serc.s2(i)+1) - 1;
    end
end
use2=serc.BA>1e-4&serc.N>11;   
[serc.res,serc.gof]=fit(serc.sl0(use2),serc.sl(use2),'poly1');
%% figure 
ax=[-0.7 3.4 -0.4 1.4];
h = figure(1);clf

x1=log10(bci.N/bci.Area);
use1=bci.BA>1e-4&bci.N>25;

subplot(321)
plot(fitlm(x1(use1),bci.sl(use1)),'marker','o','markersize',4)
xlabel([])
ylabel('slope of offset-power model')
title({'Adult ~ Saplings','Tropical (lat = 9.2)'})
axis(ax);legend('off')

subplot(322)
plot(fitlm(x1(use1),bci.slr(use1)),'marker','o','markersize',4)
xlabel([])
title({'North ~ South','Tropical (lat = 9.2)'})
axis(ax);legend('off');set(gca,'ylabel',[])
pause(.1)

x2=log10(serc.N/serc.Area);
use2=serc.BA>1e-4&serc.N>11;
subplot(323)
plot(fitlm(x2(use2),serc.sl(use2)),'marker','o','markersize',4)
xlabel('Abundance (log number of stems per hectare)')
ylabel('slope of offset-power model')
title( 'Temperate (lat = 38.9)');legend('off')
axis(ax)

subplot(324)  
plot(fitlm(x2(use2),serc.slr(use2)),'marker','o','markersize',4)
xlabel('Abundance (log number of stems per hectare)')
title( 'Temperate (lat = 38.9)');legend('off');set(gca,'ylabel',[])
axis(ax)
pause(.1)

% median slopes
use1=bci.N>25;
use2=serc.N>11;
use3=bci.N>25&bci.BA<0.1;
use4=serc.N>11&serc.BA<0.1;

B=zeros(2,4);
B(1,1)=nanmedian(bci.sl(use1));
B(1,2)=nanmedian(bci.slr(use1));
B(2,1)=nanmedian(serc.sl(use2));
B(2,2)=nanmedian(serc.slr(use2));

B(1,3)=nanmedian(bci.sl(use3));
B(1,4)=nanmedian(bci.slr(use3));
B(2,3)=nanmedian(serc.sl(use4));
B(2,4)=nanmedian(serc.slr(use4));


C1(1,:)=bootci(1000,@nanmedian,bci.sl(use1));
C1(2,:)=bootci(1000,@nanmedian,bci.slr(use1));
C2(1,:)=bootci(1000,@nanmedian,bci.sl(use3));
C2(2,:)=bootci(1000,@nanmedian,bci.slr(use3));

C3(1,:)=bootci(1000,@nanmedian,serc.sl(use2));
C3(2,:)=bootci(1000,@nanmedian,serc.slr(use2));
C4(1,:)=bootci(1000,@nanmedian,serc.sl(use4));
C4(2,:)=bootci(1000,@nanmedian,serc.slr(use4));


lower=[C1(:,1)' C2(:,1)'; C3(:,1)' C4(:,1)'];
upper=[C1(:,2)' C2(:,2)'; C3(:,2)' C4(:,2)'];

subplot(313)
ctrs=1:4;
hBar=bar(1:2,B);
ylabel('median slope of offset-power model')
set(gca,'xticklabel',{'Tropical (lat = 9.2)'; 'Temperate (lat = 38.9)'})
clear ctr ydt
for k1 = 1:size(B,2)
    ctr(:,k1) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(:,k1) = hBar(k1).YData;
end
hold on
errorbar(ctr, ydt, ydt-lower,upper-ydt, '.r')
legend('Adult ~ Saplings (all)','North ~ South (all)',...
       'Adult ~ Saplings (rare)','North ~ South (rare)','location','northwest')
legend('boxoff')

sublabel;

saveas(h,['Figures/North~South_' num2str(DX) '.fig'])


