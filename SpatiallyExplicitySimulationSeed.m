% Simulate spatially explity seed shadow and germination probability
% and fit model of CNDD
clear
%% parameters&constsnst
L=500;      %plot size
dx=1;       %trap size
f=0.1;      %germination rate
K=25;       %number of independent realization
dL=4;
M = 20;  % number of scenarios


fecundity=2*logspace(3,4,M);   % average number of seeds per tree
mothers=5;                      % number of reproductive individuals
N=mothers.*fecundity*1e-4;      % expected numbr of seed desnity
dispersal=linspace(5,20,M);    % dospersal paramter od the fat kernel

fname='SimSeeds_v03.mat';            % file name for saving results
%enviromental variability
% E=f/10;
% alfa = (f^2 -E*f -f^3)/E;
% beta = (f - E + f*E -2*f^2 +f^3)/E;
E=0;
% gemination simulations
% fecundity=10000;
% f = logspace(-2,log10(0.2),M);
%% sampling design
% Wright sampling design
trap=zeros(L,L);
trap(dL+1:2*dL+1:L-dL,dL+1:2*dL+1:L-dL)=1;

site=zeros(L,L);
site(dL+1:2*dL+1:L-dL,dL-1:2*dL+1:L-dL)=1;
site(dL-1:2*dL+1:L-dL,dL+1:2*dL+1:L-dL)=2;
site(dL+3:2*dL+1:L-dL,dL+1:2*dL+1:L-dL)=3;

% % Batghi sampling design
% dL=40;       %half spacing between traps
% trap=zeros(L,L);
% trap(dL+1:2*dL+1:L,dL+1:2*dL+1:L)=1;
% trap(dL+3:2*dL+1:L,dL-1:2*dL+1:L)=2;
% trap(dL+3:2*dL+1:L,dL+3:2*dL+1:L)=3;
% 
% site=zeros(L,L);
% site(dL+1:2*dL+1:L,dL-1:2*dL+1:L)=1;
% site(dL-1:2*dL+1:L,dL+1:2*dL+1:L)=2;
% site(dL+3:2*dL+1:L,dL+1:2*dL+1:L)=3;
% site(dL+1:2*dL+1:L,dL+3:2*dL+1:L)=4;

% map=zeros(L,L);
% map(trap>0)=1;
% map(site>0)=2;
%
% pcolor(map);shading flat;daspect([1 1 1])
% axis([1 50 1 50])
%% MAIN LOOP
ML0_M=zeros(2,M,M,5);
ML0_V=zeros(2,M,M,5);
parfor i=1:M
    disp(i)
    ml0=zeros(2,M,K,5);
    for j=1:M
        for k=1:K
            [j k]
            S = SeedShadow(mothers,fecundity(i),dispersal(j),L,dx);
            
            if E>0
                % enviromenatal variability
                A=betarnd(alfa,beta,size(S));
                R = binornd(S,A);
            else
                % stocastic demography
                R = binornd(S,f);
            end
            
            R1 = [R(site==1) R(site==2) R(site==3)];
            S1 = S(trap==1);
            S0 = [S(site==1) S(site==2) S(site==3)];
            R2 = mean(R1,2);
            
            % perfect collocation (binomial)
            use=S0>0;
            tbl=table(S0(use),R1(use),'VariableNames',{'S0','R1'});
            glm = fitglm(tbl,'R1 ~ 1 + S0','distribution','Binomial','BinomialSize',S0(use));
            ml0(:,j,k,1) = glm.Coefficients.Estimate;
            
            % collocation errror: power-offset model
            use=~(S1==0&R2==0);
            ml0([2 1],j,k,2) = polyfit(log(S1(use)+1),log(R2(use)+1),1);
            
            % collocation errror: GLM(Poisson)
            S2=repmat(S1,1,3);S2(R1>0&S2==0)=R1(R1>0&S2==0);
            use=S2>0;
            tbl=table(log(S2(use)),R1(use),'VariableNames',{'LS2','R1'});
            glm = fitglm(tbl,'R1 ~ 1 + LS2','distribution','Poisson');
            ml0(:,j,k,3) = glm.Coefficients.Estimate;

            S0 = mean([S(site==1) S(site==2) S(site==3)],2);
            R0 = R(trap==1);
            S0(R0>0&S0==0)=R0(R0>0&S0==0);
            use=S0>0;
            tbl=table(log(S0(use)),R0(use),'VariableNames',{'LS2','R1'});
            glm = fitglm(tbl,'R1 ~ 1 + LS2','distribution','Poisson');
            ml0(:,j,k,3) = glm.Coefficients.Estimate;

        end
    end
    
    ML0_M(:,:,i,:) = nanmean(ml0(:,:,:,:),3);
    ML0_V(:,:,i,:) =  nanvar(ml0(:,:,:,:),0,3);
    
end


metafile={'ML0 is a multidimentional array (pxmxnxq)',...
    'p intercept and slope of the power-law model S = a*A^b',...
    'm dispersal, n mothers, q=1,..4, fit and error types ',...
    'q = 1: Binomial-logit with no error',...
    'q = 2: power off-sett with collocation error',...
    'q = 3: GLM(Poisson) with collocation error',...
    'q = 4: empty',...
    'q = 4: empty',...
    'M: refer to average of 25 simulationfor each combination of the parameters',...
    'V: refer to variance of 25 simulationfor each combination of the parameters'}';

save(fname,'metafile','fecundity','mothers','E','f','dispersal','ML0_M','ML0_V');



%% save simulations in binary files
%    fname=['S05072018' ...
%           '_F' num2str(fecundity(i),'%.2f')...
%           '_D' num2str(dispersal(j),'%.2f')...
%           '_N' num2str(k,'%03d')...
%           '_L' num2str(L)...
%           '_M' num2str(mothers)...
%           '.bin'];
%
%             fileID = fopen(fname,'w');
%             fwrite(fileID,[S(trap==1) S(trap==2) S(trap==3) ...
%                            S(site==1) S(site==2) S(site==3) S(site==4)],'int');
%             fclose(fileID);
%
% %             fileID = fopen(fname);
% %             A = fread(fileID,[36 7],'uint');
%         end
%     end
% end


