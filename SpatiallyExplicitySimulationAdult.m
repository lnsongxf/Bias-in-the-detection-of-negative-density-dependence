% Simulate spatially explity juvenile-adult reletionships
% and fit model of CNDD
clear
%% parameters & constants
L=1000;                         % plot size
dx=20;                          % quadrat size
K=25;                           % number of independent realizations per scenario
M=20;                           % number of steps
MissClass = 0.1;                % percentage of missclassified individuals
fecundity=2;                    % average number of recruit per tree
mothers=logspace(0,1,M);        % number of reproductive individuals ha-1
dispersal=linspace(5,20,M);     % dispersal paramter od the fat kernel
fname='SimAdult_v3.1.mat';      % file name for saving results

%% MAIN LOOP
ML0_M=zeros(2,M,M,6);
ML0_V=zeros(2,M,M,6);

parfor i=1:M
    disp(i)
    pause(.1)
    ml0=zeros(2,M,K,4);
    for j=1:M
        for k=1:K
            [i j k]
            [S,A] = SeedShadow(mothers(i),fecundity,dispersal(j),L,1);
            
            Saplings0   = BlockMatrix(S,dx);
            Adults0     = BlockMatrix(A,dx);
            Adults_dw0 = smooth(A,dispersal(j),'fun','fat','dx',dx);
              
            % missclassfication
            errS = binornd(S,MissClass);
            errA = binornd(A,MissClass);
            Saplings1   = BlockMatrix(S-errS+errA,dx);
            Adults1     = BlockMatrix(A+errS-errA,dx);
            Adults_dw1 = smooth(A+errS-errA,dispersal(j),'fun','fat','dx',dx);            
            
            % power-offset model
            ml0([2 1],j,k,1) = polyfit(log(Adults0(:)+1),log(Saplings0(:)+1),1);
            ml0([2 1],j,k,2) = polyfit(log(Adults1(:)+1),log(Saplings1(:)+1),1);
            
            % GLM
            tbl=table(Saplings0(:),Adults_dw0(:),log(Adults_dw0(:)),...
                Saplings1(:),Adults_dw1(:),log(Adults_dw1(:)),...
                'VariableNames',{'S0','A0','LA0','S1','A1','LA1'});
             
            % weighted-distance: Poisson - Ricker
            glm = fitglm(tbl,'S0 ~ 1 + A0','distribution','Poisson','link','log','offset',log(Adults_dw0(:)),'Exclude',Adults_dw0(:)<=0);
            ml0(:,j,k,3) = glm.Coefficients.Estimate;
            
            glm = fitglm(tbl,'S1 ~ 1 + A1','distribution','Poisson','link','log','offset',log(Adults_dw1(:)),'Exclude',Adults_dw1(:)<=0);
            ml0(:,j,k,4) = glm.Coefficients.Estimate;
            
            % weighted-distance: Poisson - Powerlaw
            glm = fitglm(tbl,'S0 ~ 1 + LA0','distribution','Poisson','link','log','Exclude',Adults_dw0(:)<=0);
            ml0(:,j,k,5) = glm.Coefficients.Estimate;
            
            glm = fitglm(tbl,'S1 ~ 1 + LA1','distribution','Poisson','link','log','Exclude',Adults_dw1(:)<=0);
            ml0(:,j,k,6) = glm.Coefficients.Estimate;
            
            % figure(1);clf
            % subplot(221);pcolor(Adults0);shading flat
            % daspect([1 1 1]);title('Adults (original)');caxis([0 2])
            % subplot(222);pcolor(Adults1);shading flat
            % daspect([1 1 1]);title('Adults (10% classification error)');caxis([0 2])
            %
            % subplot(223);pcolor(Saplings0);shading flat
            % daspect([1 1 1]);title('Saplings (original)');caxis([0 4])
            % subplot(224);pcolor(Saplings1);shading flat;daspect([1 1 1]);
            % title('Saplings (10% classification error)');caxis([0 4])
            
        end
    end
    
    ML0_M(:,:,i,:) = nanmean(ml0(:,:,:,:),3);
    ML0_V(:,:,i,:) =  nanvar(ml0(:,:,:,:),0,3);
end



metafile={'ML0 is a multidimentional array (pxmxnxq)',...
    'p intercept and slope of the power-law model S = a*A^b',...
    'm dispersal, n mothers, q=1,..4, fit and error types ',...
    'q = 1: OLS with no error (quadrat-based)',...
    'q = 2: OLS with classification error (quadrat-based)',...
    'q = 3: GLM with no error (Ricker, distance-weighted)',...
    'q = 4: GLM with classification error (Ricker, distance-weighted)',...
    'q = 5: GLM with no error (Powerlaw, distance-weighted)',...
    'q = 6: GLM with classification error (Powerlaw, distance-weighted)',...
    ['M: refer to average of ' num2str(K) ' simulationfor each combination of the parameters'],...
    ['V: refer to variance of ' num2str(K) ' simulationfor each combination of the parameters']}';

save(fname,'metafile','fecundity','mothers','dispersal','MissClass','ML0_M','ML0_V','L','K');


