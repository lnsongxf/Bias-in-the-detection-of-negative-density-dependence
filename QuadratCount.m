function N = QuadratCount(x,y,L,dx,varargin)

%create a raster from a point process
%INPUTS: x,y: points coordinates
%        L:   plot dimensions (e.g. [1000,500], or 1000)
%        dx: grid size (can be a 2 elemrnt vector)
%
% N = QuadratCount(x,y,L,dx) return a raster of counts for each quadrat



if nargin==4
    
    N=histcounts2(x,y,0:dx:L(1),0:dx:L(2)).';
    
    %% basal area
elseif nargin==5
    XI=0:dx:L(1);
    YI=0:dx:L(2);
    I=length(XI)-1;
    J=length(YI)-1;
    
    N=zeros(J,I);
    %basal area (in m2/m2)
    A=varargin{1};
    
    for i=1:I-1
        use1=find(x>=XI(i) & x<XI(i+1));
        for j=1:J-1
            use2=y(use1)>=YI(j) & y(use1)<YI(j+1);
            if sum(use2); N(j,i)=sum(A(use1(use2))); end
        end
    end
    
    use1=find(x>=XI(I) & x<=XI(I+1));
    for j=1:J-1
        use2=y(use1)>=YI(j) & y(use1)<YI(j+1);
        if sum(use2); N(j,I)=sum(A(use1(use2))); end
    end
    
    use1=find(y>=YI(J) & y<=YI(J+1));
    for i=1:I-1
        use2=x(use1)>=XI(i) & x(use1)<XI(i+1);
        if sum(use2); N(J,i)=sum(A(use1(use2))); end
    end
    
    use1=x>=XI(I) & x<=XI(I+1) & y>=YI(J) & y<=YI(J+1);
    if sum(use1); N(J,I)=sum(A(use1)); end

    
end

end



