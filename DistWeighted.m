function [saplings,adultDistWeighted,N,BasalArea]=DistWeighted(sp,dbh,gx,gy,Lx,Ly,DX,tr,quant,L,type)

% parameters
% quant=0.5;
% mean dispersal = pi/2*L
% L=2*20/pi;%
% tr=25;

% quadrat grid
x=0:DX:Lx;
y=0:DX:Ly;
n=length(x)-1;
m=length(y)-1;

% inizialize variables
species=unique(sp);
S=length(species);
saplings=zeros(n*m,S);
adultDistWeighted=zeros(n*m,S);
BasalArea=nan(S,1);
N=nan(S,1);


for s=1:S
    S-s
    use=strcmp(sp,species(s))&dbh>0;
    N(s,1)=sum(use);
    BasalArea(s,1)=pi/4*sum((dbh(use)/1000).^2)/(Lx*Ly/10000);
    
    if N(s,1)>tr
        x0=gx(use);y0=gy(use);
        D=dbh(use);
        if strcmp(type,'random')
            D=D(randperm(length(D)));
        end

        D50=quantile(D,quant);
        
        % get adults and saplings position
        use1=D<=D50;%sap
        use2=D>D50;%adults
        x1=x0(use1);y1=y0(use1);%sap
        x2=x0(use2);y2=y0(use2);%adults
        
        for i=1:n
            use1=find(x1>=x(i)&x1<x(i)+DX);
            dx=x(i)-x2+DX/2;
            for j=1:m
                LinInd = sub2ind([n m],i,j);
                saplings(LinInd,s)=sum(y1(use1)>=y(j)&y1(use1)<y(j)+DX);
                
                dy=y(j)-y2+DX/2;
                r2=dx.*dx+dy.*dy;
                adultDistWeighted(LinInd,s) = sum(1./(pi*L^2).*(L^2./(r2+L^2)).^2)*DX^2;
                
            end
        end
        
    end
end
