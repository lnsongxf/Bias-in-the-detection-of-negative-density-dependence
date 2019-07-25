%Smooth a surface with a 2D filter using Fourier convolution
%
%
% Matteo Detto, PhD
% Biometeorology lab, Ecosystem Science Division
% Dept. of Environmental Science, Policy, and Management (ESPM)
% University of California, Berkeley, CA 94720
% Phone: 1-510-642 9048
% Fax: 1-510-643-5098
% last update: 5 October 2018
% -------------------------------------------------------------------------
%   Copyright (C) 2007-2019, Matteo Detto
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.


function w = smooth(s,scale,varargin)


%default options
opts=struct('fun','gaussian','dx',1,'graph','n');
opts=parseArgs(varargin,opts);

[m,n]=size(s);

f=fft2(s);

npuls_2   = floor((n-1)/2);
pulsx     = 2*pi/n*[ 0:npuls_2  (npuls_2-n+1):-1 ];
npuls_2   = floor((m-1)/2);
pulsy     = 2*pi/m*[ 0:npuls_2  (npuls_2-m+1):-1 ];

[kx,ky] = meshgrid(pulsx,pulsy);
if length(scale)==2
    h=exp(-1/2*(kx*scale(1)).^2).*exp(-1/2*(ky*scale(2)).^2);
else
    K=sqrt(kx.^2+ky.^2);
    if strcmp(opts.fun,'exp')
        h = (1+(K*scale).^2).^(-3/2);
    elseif strcmp(opts.fun,'fat')
        h = scale*K.*besselk(1, scale*K);h(1)=1;
    elseif strcmp(opts.fun,'gaussian')
        h = exp(-1/2*(K*scale(1)).^2);
    end
end
w=ifft2(f.*h);
if opts.dx>1
    w  = BlockMatrix(w,opts.dx);
end

if opts.graph=='y'
    pcolor(w);shading interp
    daspect([1 1 1])
end

%%
function    B = BlockMatrix(A,dx)

[M,N]=size(A);

n=M/dx;
m=N/dx;
B = zeros(n,m);
for i=1:n
    for j=1:m
        temp = A((i-1)*dx+1:i*dx,(j-1)*dx+1:j*dx);
        B(i,j)=sum(temp(:));
    end
end


end

function ArgStruct=parseArgs(args,ArgStruct,varargin)

persistent matlabver

if isempty(matlabver)
    matlabver=ver('MATLAB');
    matlabver=str2num(matlabver.Version);
end

Aliases={};
FlagTypeParams='';

if (~isempty(varargin)) 
    FlagTypeParams=lower(strvcat(varargin{1}));
    if length(varargin)>1
        Aliases=varargin{2};
    end
end
 

%---------------Get "numeric" arguments
NumArgCount=1;
while (NumArgCount<=size(args,2))&(~ischar(args{NumArgCount}))
    NumArgCount=NumArgCount+1;
end
NumArgCount=NumArgCount-1;
if (NumArgCount>0)
    ArgStruct.NumericArguments={args{1:NumArgCount}};
else
    ArgStruct.NumericArguments={};
end 


%--------------Make an accepted fieldname matrix (case insensitive)
Fnames=fieldnames(ArgStruct);
for i=1:length(Fnames)
    name=lower(Fnames{i,1});
    Fnames{i,2}=name; %col2=lower
    AbbrevIdx=find(Fnames{i,1}~=name);
    Fnames{i,3}=[name(AbbrevIdx) ' ']; %col3=abreviation letters (those that are uppercase in the ArgStruct) e.g. SpacingHoriz->sh
    %the space prevents strvcat from removing empty lines
    Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams)); %Does this parameter have a value?
end
FnamesFull=strvcat(Fnames{:,2});
FnamesAbbr=strvcat(Fnames{:,3});

if ~isempty(Aliases)  
    for i=1:length(Aliases)
        name=lower(Aliases{i,1});
        FieldIdx=strmatch(name,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(name,FnamesFull); %&??????? exact or not? 
        end
        Aliases{i,2}=FieldIdx;
        AbbrevIdx=find(Aliases{i,1}~=name);
        Aliases{i,3}=[name(AbbrevIdx) ' ']; %the space prevents strvcat from removing empty lines
        Aliases{i,1}=name; %dont need the name in uppercase anymore for aliases
    end
    %Append aliases to the end of FnamesFull and FnamesAbbr
    FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1})); 
    FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3}));
end

%--------------get parameters--------------------
l=NumArgCount+1; 
while (l<=length(args))
    a=args{l};
    if ischar(a)
        paramHasValue=1; % assume that the parameter has is of type 'param',value
        a=lower(a);
        FieldIdx=strmatch(a,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(a,FnamesFull); 
        end
        if (length(FieldIdx)>1) %shortest fieldname should win 
            [mx,mxi]=max(sum(FnamesFull(FieldIdx,:)==' ',2));
            FieldIdx=FieldIdx(mxi);
        end
        if FieldIdx>length(Fnames) %then it's an alias type.
            FieldIdx=Aliases{FieldIdx-length(Fnames),2}; 
        end
        
        if isempty(FieldIdx) 
            error(['Unknown named parameter: ' a])
        end
        for curField=FieldIdx' %if it is an alias it could be more than one.
            if (Fnames{curField,4})
                if (l+1>length(args))
                    error(['Expected a value for parameter: ' Fnames{curField,1}])
                end
                val=args{l+1};
            else %FLAG PARAMETER
                if (l<length(args)) %there might be a explicitly specified value for the flag
                    val=args{l+1};
                    if isnumeric(val)
                        if (numel(val)==1)
                            val=logical(val);
                        else
                            error(['Invalid value for flag-parameter: ' Fnames{curField,1}])
                        end
                    else
                        val=true;
                        paramHasValue=0; 
                    end
                else
                    val=true;
                    paramHasValue=0; 
                end
            end
            if matlabver>=6
                ArgStruct.(Fnames{curField,1})=val; %try the line below if you get an error here
            else
                ArgStruct=setfield(ArgStruct,Fnames{curField,1},val); %<-works in old matlab versions
            end
        end
        l=l+1+paramHasValue; %if a wildcard matches more than one
    else
        error(['Expected a named parameter: ' num2str(a)])
    end
end

end

end