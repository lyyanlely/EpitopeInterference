%function BCellnp1(n,mu,sig,gam,a,b,rp)
% the function simulate the population dynamics of BCell-virus dynamics
% viral antigen V'=sigma-V sum_i (a_i+sum_j A_ij B_j)B_i + noise?
% B-cell B'_i=r V(a_i+sum_j A_ij B_j)B_i-B_i sum_j B_j/K+lambda (division)
% A is random matrix with mean mu/n, variance sigma^2/n, symmetry gamma
% a is bimodal distribution with probability p at low value

n = 6; % number of competing fit epitopes
%mu  = 0;
%sig = 2.;
%gam = 0;
a0  = 2.;  % a-b
a1  = 0.;  % a+b

h = 0.02; % delta time
T = 500*n;
nRec = 1;
replica = 1;
nRc  = 1000;
eps  = 1e-8;

%p  = .5; %2/n;
al = 0.9;
ah = 1;
%
lambda = eps;
sigma  = 1;
K = 2;  % normalized capacity

totStep = ceil(T/h);
Esize = ceil(nRec/replica);
nh = ceil(totStep/nRc);

s = rng('shuffle');

intf = zeros(Esize*(n),n+1);

xdata = zeros(Esize*replica,n);
Vdata = zeros(Esize*replica,1);
flags = zeros(1,Esize*replica);

tic;
for i = 1:Esize
    
    a = zeros(n,1);
    nl = 1; %round(n*p); %min(poissrnd(p*n),n);
    a(1:nl) = al; %*(1+0.1*randn(1,nl));
    a(nl+1:n) = ah; %*(1+0.1*randn(1,n-nl));
    Amat = zeros(n);
    %LVmat = -abs(randn(n));
    %LVmat = LVmat - mean(mean(LVmat))*ones(n); %remove the strong mutualism
    genm = sig/sqrt(n-1)*randn(n-1);
    Amat(2:n,2:n) = -(1+mu/(n-1))*eye(n-1)+mu/(n-1)+sqrt((1+gam)/2)*triu(genm,1)+sqrt((1-gam)/2)*tril(genm,-1)'+sqrt((1+gam)/2)*triu(genm,1)'-sqrt((1-gam)/2)*tril(genm,-1); %asymmetry characterized by gam
    Amat(1,1) = -1;
    Amat(1,2:n) = -a1;
    Amat(2:n,1) = -a0;
    
    intf((i-1)*n+(1:n),:) = [a,Amat];
    
    for r = 1:replica
        %%
        x = rand(n,1); %*1e-6;
        V = 1;
        
        xt = zeros(n,nRc);
        Vt = zeros(1,nRc);
        rt = zeros(n,nRc);
        
        for t=1:totStep
            dV = sigma-sum(a.*(1+Amat*x).*x)*V;
            dx = V*a.*(1+(Amat*x)).*x-x*sum(x)/K+lambda;
            if any(isnan(x)) || any(isnan(dx)) || isnan(dV) || isnan(V)
                break;
            end
            if abs(dV)<eps && all(abs(dx)<eps)
                flags((i-1)*replica+r) = 1;
                break
            end
            % midpoint method
            Vtemp = V+h/2*dV;
            xtemp = x+h/2*dx;
            nid   = xtemp<0;
            xtemp(nid)= 0;
            V     = V+h*(sigma-sum(a.*(1+Amat*xtemp).*xtemp)*Vtemp);
            x     = x+h*(Vtemp*a.*(1+(Amat*xtemp)).*xtemp-xtemp*sum(xtemp)/K+lambda);
            nid   = x<0|abs(dx*h)>x;
            x(nid)= 0;
            % record
            %
            if mod(t,nh)==0
                xt(:,t/nh) = x;
                Vt(t/nh)   = V;
                rt(:,t/nh) = a+Amat*x;
            end
            %}
        end
        %%
        xdata((i-1)*replica+r,:) = x';
        Vdata((i-1)*replica+r)   = V;
        % participation ratio
        
    end
end
toc;
%% save data
rp = 1;
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