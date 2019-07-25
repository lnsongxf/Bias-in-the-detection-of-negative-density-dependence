% simulate the seed shadow in a LxL square arena
%
% parmeters:
% N: number of reproductive adult per hectare
% f: fecundity per tree
% hd:  dispersal distance
% L: plot size
% dx: seed plot size
% VM: overdispersion parameter to to create variability in seed production across reproductive individuals
% (to be changed in the code at line 19)

function [S,M,xn,yn,x,y] = SeedShadow(N,f,hd,L,dx)

N = round(N*L^2/10000);
x = rand(N,1)*L;
y = rand(N,1)*L;


VM = 2;  %overdispersion parameter (ratio of variance over mean, e.g. 1 for a poisson process)
if VM==1
    S = poissrnd(f,N,1);%random seeds per parent
else
    r = f./(VM-1);
    p=1/VM;
    S = nbinrnd(r,p,N,1);%random seeds per parent
end

N1 = sum(S);  
zfat=-4;
%Fat kernel K = -heta+2/(2*pi*L^2)*(1+x^2/L^2)^-2  follow Rosendal 2009
r=rand(N1,1);
R=hd*sqrt((1-r).^(2/(2+zfat))-1);
a=2*rand(N1,1)-1;
b=2*rand(N1,1)-1;
h=a.*a+b.*b;

use=h>1;N0=sum(use);
while N0>0
    a(use)=2*rand(N0,1)-1;
    b(use)=2*rand(N0,1)-1;
    h(use)=a(use).*a(use)+b(use).*b(use);
    use=h>1;N0=sum(use);
end
h1=sqrt(h);
Z(:,1)=R.*a./h1;
Z(:,2)=R.*b./h1;

xn = zeros(N1,1);
yn = zeros(N1,1);
b=1;
for i=1:N
    use=b:b+S(i)-1;
    xn(use)=x(i)+Z(use,1);
    yn(use)=y(i)+Z(use,2);
    b=b+S(i);
end

%periodic conditions
xn(xn>L)=xn(xn>L)-L;
yn(yn>L)=yn(yn>L)-L;

xn(xn<0)=xn(xn<0)+L;
yn(yn<0)=yn(yn<0)+L;

%out of the boundaries
im=xn>L | xn<0 | yn>L | yn<0;
xn(im)=rand(sum(im),1)*L;
yn(im)=rand(sum(im),1)*L;

if nargout==0  %plot only
    plot(xn,yn,'.','markersize',1);
    hold all
    plot(x,y,'.','markersize',12);
    axis([0 L 0 L])
    daspect([1 1 1])
    pause(.1)
    title(['seed: ' num2str(N*f/L^2,'%2.1f') ' (m^{-2}), dispersal: ' num2str(round(hd)) ' (m)'])
    
elseif nargout==1   
    
    S = QuadratCount(xn,yn,[L L],dx);
elseif nargout>1

    S = QuadratCount(xn,yn,[L L],dx);
    M = QuadratCount(x,y,[L L],dx);
end

end



