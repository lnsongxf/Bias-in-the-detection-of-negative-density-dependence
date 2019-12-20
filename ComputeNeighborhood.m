%% Compute neighborhood density for survival analysis
clear
%% load BCI 8 census data 
% the data can be dowloaded from
% Condit, Richard et al. (2019), Complete data from the Barro Colorado 50-ha plot: 
% 423617 trees, 35 years, v3, DataONE, Dataset, https://doi.org/10.15146/5xcp-0d46
% data needs to be converted to MAT-files with names bci*.mat, where the arsterics indicates the census number (see Rcode below)
% require(R.matlab)
% load("bci.tree1.rdata")
%
% dat = subset(bci.tree1,status=='A')
%
% writeMat("bci1.mat",treeID=dat$treeID,
%                     sp=dat$sp,
%                     dbh=dat$dbh,
%                     gx=dat$gx,
%                     gy=dat$gy)
%
% each file contains the following variables
% sp: 6 letter species code
% dbh: diamter at breast height or (mm)
% gx, gy: tree coordinates in the local grid
% treeID: unique tree identifier

b = cell(8,1);
for i=1:8
    b{i}=load(['Data\bci' num2str(i) '.mat']);
end

species = unique(b{7}.sp);
S = length(species);
N0 = zeros(S,8);
Dmax = nan(S,8);
for j=1:8
    for i=1:S
        use = strcmp(b{j}.sp,species(i))&b{j}.dbh>0;
        if sum(use)>0
        N0(i,j) = sum(use);
        Dmax(i,j) = max(b{j}.dbh(use));
        end
    end
end

% select species with >250 individuals/census and max DBH>200 mm
I=find(mean(N0,2)>250 & max(Dmax,[],2)>200);

BA=zeros(length(I),8);
N=zeros(length(I),8);
Dmax=zeros(length(I),8);
for j=1:8
    for i=1:length(I)
        use = strcmp(b{j}.sp,species(I(i)))&b{j}.dbh>0;
        BA(i,j)=pi/4*sum((b{j}.dbh(use)*1e-3).^2)/50e+4;
        N(i,j)=sum(use);
        Dmax(i,j)=max(b{j}.dbh(use));
    end
end

BAavg=mean(BA,2);
Navg=mean(N,2);
Dmax=max(Dmax,[],2);
S = species(I);

save('SurvivalAnalysis.mat','BAavg','Navg','Dmax','S')
%% MAIN LOOP
Lx=1000;Ly=500;
scale = 5;
allometry = 2;
nZ=zeros(length(I),1);
tr = 100; %size treshold (mm)
buffer = 2*scale;
% compute vector size for pre-allocations
for i=1:length(I)
    for h=1:7
        % select focal individuals in the initial census, exclude a buffer zone
        c1 = strcmp(b{h}.sp,species(I(i))) & b{h}.dbh>0 &  b{h}.dbh<=tr & ...
                    b{h}.gx>=buffer & b{h}.gx<=Lx-buffer & ...
                    b{h}.gy>=buffer & b{h}.gy<=Ly-buffer;
        nZ(i)=nZ(i)+sum(c1);
    end
end


for i=1:length(I)
    % pre-allocate variables
    SP     = species(I(i));
    D      = zeros(nZ(i),1);
    surv   = zeros(nZ(i),1);
    Z1     = zeros(nZ(i),1);
    Z2     = zeros(nZ(i),1);
    census = zeros(nZ(i),1);

    n0=0;
    for h=1:7
        [i h]
        pause(.1)
        % select focal individuals in the initial census, exclude a buffer zone
        c1 = find(strcmp(b{h}.sp,SP) & b{h}.dbh>0 & b{h}.dbh<=tr & ...
                         b{h}.gx>=buffer & b{h}.gx<=Lx-buffer & ...
                         b{h}.gy>=buffer & b{h}.gy<=Ly-buffer);
        % select focal individuals in the following census, exclude a buffer zone
        c2 = find(strcmp(b{h+1}.sp,SP) & b{h+1}.dbh>0 & ...
                         b{h+1}.gx>=buffer & b{h+1}.gx<=Lx-buffer & ...
                         b{h+1}.gy>=buffer & b{h+1}.gy<=Ly-buffer);
        
        [~,bbb] = intersect(b{h}.treeID(c1),b{h+1}.treeID(c2));
        temp = zeros(length(c1),1);
        n1 = length(c1);
        temp(bbb) = 1;
        use = n0+1:n0+n1;n0=n0+n1;
        surv(use) = temp;
        D(use) = b{h}.dbh(c1);
        census(use)=h;
            
        % conspecific and heterospecific basal area(including individuals in the buffer zone)
        c0 = strcmp(b{h}.sp,SP) & b{h}.dbh>0 & b{h}.gx>=0;
        BCON = pi/4*(b{h}.dbh(c0)/10).'.^allometry;
        h0 = ~strcmp(b{h}.sp,SP) & b{h}.dbh>0 & b{h}.gx>=0;
        BHET = pi/4*(b{h}.dbh(h0)/10).'.^allometry;

        % loop within focal individuals to estimate local neghborhhod
        z1=zeros(n1,1);z2=zeros(n1,1);
        parfor j=1:n1   
        d2 = (b{h}.gx(c0)-b{h}.gx(c1(j))).^2 + (b{h}.gy(c0)-b{h}.gy(c1(j))).^2;  %distance to focal sapling
        w = exp(-sqrt(d2)/scale)/scale;
        z1(j)  = BCON*w - pi/4/scale*(b{h}.dbh(c1(j))/10).^allometry;  %conspecific neighborood density
        
        d2 = (b{h}.gx(h0)-b{h}.gx(c1(j))).^2 + (b{h}.gy(h0)-b{h}.gy(c1(j))).^2;  %distance to focal sapling
        w = exp(-sqrt(d2)/scale)/scale;
        z2(j)  = BHET*w;  %heterspecific neighborood density

        end
        Z1(use)=z1;
        Z2(use)=z2; 
    end
    
      save(['SaplingSurvSP_' num2str(i,'%02d') '.mat'],'Z1','Z2','surv','D','scale','census','allometry','tr')
end
