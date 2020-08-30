gure%function BCell(n,mu,sig,gam,rp)
% the function simulate the population dynamics of BCell-virus dynamics
% viral antigen V'=sigma-V sum_i (a_i+sum_j A_ij B_j)B_i + noise?
% B-cell B'_i=r V(a_i+sum_j A_ij B_j)B_i-B_i sum_j B_j/K+lambda (division)
% A is random matrix with mean mu/n, variance sigma^2/n, symmetry gamma
% a is bimodal distribution with probability p at low value
% the program is a duplication of BCell.m limited to the two B cell
% species case, one will need to specify the interaction matrix A

h = 0.02; % delta time
T = 200;
nRec = 5000;
eps  = 1e-8;

n = 2; % number of competing epitopes
a0  = -.5:0.1:1;
a1  = -1:0.1:.5;

l0 = length(a0);
l1 = length(a1);

al = .5;
ah = 1;
beta = 0.13;
rtot = 100;
%
lambda = eps;
sigma  = 1;
K = 6;  % normalized capacity

totStep = ceil(T/h);
%Esize = ceil(nRec/replica);
nh = ceil(totStep/nRec);
phi0 = zeros(l0,l1);
phi1 = zeros(l0,l1);

s = rng('shuffle');

%{
intf = zeros(Esize*n,n+1);

xdata = zeros(Esize*replica,n);
Vdata = zeros(Esize*replica,1);
flags = zeros(1,Esize*replica);
%}
a = [al;ah];
for i = 1:l0
    for j = 1:l1
        tic;
%LVmat = -abs(randn(n));
%LVmat = LVmat - mean(mean(LVmat))*ones(n); %remove the strong mutualism
%genm = sig/sqrt(n)*randn(n);
Amat = -[1,a1(j);a0(i),beta]; %asymmetry characterized by gam
    
tic;
%    for r = 1:replica
x = rand(n,1)*1e-3;
V = 1;

xt = zeros(n,nRec);
Vt = zeros(1,nRec);
rt = zeros(n,nRec);

for r = 1:rtot
for t=1:totStep
    dV = sigma-sum(a.*(1+Amat*x).*x)*V;
    dx = V*a.*(1+(Amat*x)).*x-x*sum(x)/K+lambda;
            if any(isnan(x)) || any(isnan(dx)) || isnan(dV) || isnan(V)
                break;
            end
    if abs(dV)<eps && all(abs(dx)<eps)
        %flags((i-1)*replica+r) = 1;
        break
    end
    % midpoint method
    Vtemp = V+h/2*dV;
    xtemp = x+h/2*dx;
    Vxi   = (sigma-sum(a.*(1+Amat*xtemp).*xtemp)*Vtemp); %sqrt(Vtemp*eps)*randn(1);
    xxi   = (Vtemp*a.*(1+(Amat*xtemp)).*xtemp-xtemp*sum(xtemp)/K+lambda); %sqrt(xtemp*eps).*randn(n,1);
    V     = V+h*(sigma-sum(a.*(1+Amat*xtemp).*xtemp)*Vtemp); %+Vxi
    x     = x+h*(Vtemp*a.*(1+(Amat*xtemp)).*xtemp-xtemp*sum(xtemp)/K+lambda); %+xxi
    nid   = x<0;
    x(nid)= 0;
    if abs(Vxi*h/V)>1 || max(abs(xxi*h./x))>1
        break
    end
    % perturb
    %if mod(t,totStep/2)==0
    %    x(2) = x(2)+0.5;
    %end
    % record
    %
    if mod(t,nh)==0
        xt(:,t/nh) = x;
        Vt(t/nh)   = V;
        rt(:,t/nh) = a+Amat*x;
    end
    %}
end
%
if x(1)>0.01
phi0(i,j) = phi0(i,j)+1;
end
if x(2)>0.01
phi1(i,j) = phi1(i,j)+1;
end
%}
end
%xdata((i-1)*replica+r,:) = x';
%Vdata((i-1)*replica+r)   = V;
% participation ratio
toc;
    end
end
        %{
    end
    toc;
end

dirc = './';
xname = 'BVconcentration';
iname = 'Interference';
fname = 'Exitflag';
nname = sprintf('n%d',n);
mname = sprintf('mu%.2f',mu);
sname = sprintf('sig%.2f',sig);
gname = sprintf('gam%.2f',gam);
pname = sprintf('p%.2f',p);
aname = sprintf('a%.2f',al);
kname = sprintf('K%.2f',K);
rname = sprintf('%02d',rp);
dtype = '.dat';
conname = [dirc,xname,'_',nname,'_',mname,'_',sname,'_',gname,'_',pname,'_',aname,'_',kname,'_',rname,dtype];
intname = [dirc,iname,'_',nname,'_',mname,'_',sname,'_',gname,'_',pname,'_',aname,'_',kname,'_',rname,dtype];
flgname = [dirc,fname,'_',nname,'_',mname,'_',sname,'_',gname,'_',pname,'_',aname,'_',kname,'_',rname,dtype];

dlmwrite(conname,[Vdata,xdata]);
dlmwrite(intname,intf);
dlmwrite(flgname,flags);
        %}