clear

% parameters
quant=0.5;
L=2*20/pi;%mean dispersal = pi/2*L
DX=20;
model='power';
%% BCI: load data
load('bci.mat')
species=unique(sp);
S=length(species);
Lx=1000;Ly=500;bci.Area=50;
tr=25;
%% BCI: analyze
bci.sl=nan(S,1);
bci.slr=nan(S,1);

[saplings,adultDistWeighted,bci.N,bci.BA]=DistWeighted(sp,dbh,gx,gy,Lx,Ly,DX,tr,quant,L,'original');
[saplingsR,adultDistWeightedR]=DistWeighted(sp,dbh,gx,gy,Lx,Ly,DX,tr,quant,L,'random');
for i=1:S
    i
    if bci.N(i)>25
        
     if strcmp(model,'Ricker')
     tbl=table(adultDistWeighted(:,i),saplings(:,i),adultDistWeightedR(:,i),saplingsR(:,i),'VariableNames',{'A','S','AR','SR'});
     glm = fitglm(tbl,'S ~ 1 + A','distribution','Poisson','offset',log(adultDistWeighted(:,i)));
     bci.sl(i) = exp(glm.Coefficients.Estimate(2));
     glm = fitglm(tbl,'SR ~ 1 + AR','distribution','Poisson','offset',log(adultDistWeightedR(:,i)));
     bci.slr(i) = exp(glm.Coefficients.Estimate(2));   
     elseif strcmp(model,'power')
     tbl=table(log(adultDistWeighted(:,i)),saplings(:,i),log(adultDistWeightedR(:,i)),saplingsR(:,i),'VariableNames',{'A','S','AR','SR'});
     glm = fitglm(tbl,'S ~ 1 + A','distribution','Poisson');
     bci.sl(i) = glm.Coefficients.Estimate(2);
     glm = fitglm(tbl,'SR ~ 1 + AR','distribution','Poisson');
     bci.slr(i) = glm.Coefficients.Estimate(2);   
     end
    end
end

%% SERC: load data
load('serc.mat')
species=unique(sp);
S=length(species);
Lx=400;Ly=400;serc.Area=16;
tr=11;
%% SERC: analysis
serc.sl=nan(S,1);
serc.slr=nan(S,1);

[saplings,adultDistWeighted,serc.N,serc.BA]=DistWeighted(sp,dbh,gx,gy,Lx,Ly,DX,tr,quant,L,'original');
[saplingsR,adultDistWeightedR]=DistWeighted(sp,dbh,gx,gy,Lx,Ly,DX,tr,quant,L,'random');
for i=1:S
    i

    if serc.N(i)>11
     if strcmp(model,'Ricker')
     tbl=table(adultDistWeighted(:,i),saplings(:,i),adultDistWeightedR(:,i),saplingsR(:,i),'VariableNames',{'A','S','AR','SR'});
     glm = fitglm(tbl,'S ~ 1 + A','distribution','Poisson','offset',log(adultDistWeighted(:,i)));
     serc.sl(i) = exp(glm.Coefficients.Estimate(2));
     glm = fitglm(tbl,'SR ~ 1 + AR','distribution','Poisson','offset',log(adultDistWeightedR(:,i)));
     serc.slr(i) = exp(glm.Coefficients.Estimate(2));   
     elseif strcmp(model,'power')
     tbl=table(log(adultDistWeighted(:,i)),saplings(:,i),log(adultDistWeightedR(:,i)),saplingsR(:,i),'VariableNames',{'A','S','AR','SR'});
     glm = fitglm(tbl,'S ~ 1 + A','distribution','Poisson');
     serc.sl(i) = glm.Coefficients.Estimate(2);
     glm = fitglm(tbl,'SR ~ 1 + AR','distribution','Poisson');
     serc.slr(i) = glm.Coefficients.Estimate(2);   
     end
    end
end

%% figure 
ax=[-0.7 3.4 -0.4 1.4];
h = figure(1);clf

x1=log10(bci.N/bci.Area);
use1=~isnan(bci.N);

subplot(321)
plot(fitlm(x1(use1),bci.sl(use1)),'marker','o','markersize',4)
xlabel([])
ylabel('slope of offset-power model')
title({['Adult ~ Saplings (' num2str(DX) 'x' num2str(DX) ' m)'],'Tropical (lat = 9.2)'})
axis(ax);legend('off')

subplot(322)
plot(fitlm(x1(use1),bci.slr(use1)),'marker','o','markersize',4)
xlabel([])
title({['Odd ~ Even (' num2str(DX) 'x' num2str(DX) ' m)'],'Tropical (lat = 9.2)'})
axis(ax);legend('off');set(gca,'ylabel',[])
pause(.1)

x2=log10(serc.N/serc.Area);
use2=~isnan(serc.N);
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

%% median slopes
use1=~isnan(bci.N);
use2=~isnan(serc.N);
use3=~isnan(bci.N)&bci.BA<0.1;
use4=~isnan(serc.N)&serc.BA<0.1;

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
ctr=zeros(2,4); ydt=zeros(2,4);
for k1 = 1:size(B,2)
    ctr(:,k1) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(:,k1) = hBar(k1).YData;
end
hold on
errorbar(ctr, ydt, ydt-lower,upper-ydt, '.r')
legend('Adult ~ Saplings (all)','Odd ~ Even (all)',...
       'Adult ~ Saplings (rare)','Odd ~ Even (rare)','location','best')
legend('boxoff')

