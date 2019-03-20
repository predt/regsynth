clear variables;
close all;
rng(0);
% parpool('reducedform',100);
% parpool('local',4);

% n1 = 100;
% n0 = 1000;
n1 = 10;
n0 = 100;
k = 10;
a = 0.10; b = 0.90;
r = 2;
h = 0.10;

m = 0; 
for t=0:r
    m = m + nchoosek(r,t)*(a^(r-t))*((b-a)^t)/(t+1);
end
m2 = 0; 
for t=0:2*r
    m2 = m2 + nchoosek(2*r,t)*(a^(2*r-t))*((b-a)^t)/(t+1);
end

v = sqrt(k*(m2-m^2));

Lambda = 0.01:0.01:3.5;

% options = optimoptions('quadprog','StepTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',2000,'Display','off');
 options = optimoptions('quadprog','TolFun',1e-10,'TolX',1e-10,'MaxIter',2000,'Display','off');

T = 100;
MSEp = zeros(T,1); MSEnp = zeros(T,1); MSEm = zeros(T,1); MSEmopt = zeros(T,1);
lambdavalues = zeros(T,1);
M = 20; % nb. of neighbors to consider in matching

mvalues = zeros(T,1);
for t =1:T
    tic
    sprintf('Iteration: %d',t)
    
    x1 = a+(b-a)*rand(k,n1);
    x0 = (a-h)+(b-a+2*h)*rand(k,n0);
    x0 = sqrt(x0);
    y1 = sum(x1.^r,1)'/v+randn(n1,1);
    y0 = sum(x0.^r,1)'/v+randn(n0,1);

    H = 2*(x0'*x0);
    
    minMSE = inf;
    Wp = zeros(n0,n1);
    optlambda = 0;
    
    % 1. pen synth w/ optimized lambda
    
    for l = 1:length(Lambda)
        lambda = Lambda(l);
        W = zeros(n0,n1);
        for i=1:n1
            x = x1(:,i);
            D = x0 - kron(ones(1,n0),x);
            delta = diag(D'*D); 
            f = lambda*delta-2*x0'*x;
            w = quadprog(H,f,[],[],ones(1,n0),1,zeros(n0,1),ones(n0,1),[],options);
            W(:,i) = w;
        end
        mse = (y1-W'*y0)'*(y1-W'*y0);
        if mse < minMSE
            minMSE = mse;
            optlambda = lambda;
            Wp = W;
        end
    end
    
    lambdavalues(t) = optlambda;
    
    % 2. matching
    
    % 2.1 collects the indices of the m closet points in 'matches'
    matches = zeros(M,n1);
    for i=1:n1
         x = x1(:,i);
         D = x0 - kron(ones(1,n0),x);
         delta = diag(D'*D); 
         [sorted,I]=sort(delta);
         matches(:,i) = I(1:M);
    end
    
    % 2.2 computes the corresponding MSE
    minMSEm = inf;
    Wmopt = zeros(n0,n1);
    optm = 0;
    dim = size(Wmopt); 
    for m = 1:M
        W = zeros(n0,n1);
        W(sub2ind(dim,matches(1:m,:),ones(size(matches(1:m,:),1),1)*(1:dim(2))))=1/m;
        mse = (y1-W'*y0)'*(y1-W'*y0);
        if mse < minMSEm
            minMSEm = mse;
            optm = m;
            Wmopt = W;
        end
    end
    
    mvalues(t) = optm;
    
    % 3. Evaluate performance on subsequent period
    
    Wm = zeros(n0,n1);
    Wm(sub2ind(dim,matches(1,:),ones(size(matches(1,:),1),1)*(1:dim(2))))=1;
    
    
    y1 = sum(x1.^r,1)'/v+randn(n1,1);
    y0 = sum(x0.^r,1)'/v+randn(n0,1);
    
    % 3.1 non pen synth
    Wnp = zeros(n0,n1);
    for i=1:n1
        x = x1(:,i);
        D = x0 - kron(ones(1,n0),x);
        delta = diag(D'*D); 
        w = quadprog(H,-2*x0'*x,[],[],ones(1,n0),1,zeros(n0,1),ones(n0,1),[],options);
        Wnp(:,i) = w;
    end
    MSEp(t) = (y1-Wp'*y0)'*(y1-Wp'*y0);
    MSEnp(t) = (y1-Wnp'*y0)'*(y1-Wnp'*y0);
    MSEm(t) = (y1-Wm'*y0)'*(y1-Wm'*y0);
    MSEmopt(t) = (y1-Wmopt'*y0)'*(y1-Wmopt'*y0);
    
    toc
end

filename = sprintf('n1_%d_n0_%d_k_%d_r_%d_%d%d_%d_T_%d_nested',n1,n0,k,r,100*a,100*b,100*h,T);
save(filename,'n1','n0','k','r','a','b','h','MSEp','MSEnp','MSEm','MSEmopt','Lambda','M','lambdavalues','mvalues','T');

delete(gcp('nocreate')); 

%sum(W>.0000000001,1)

%vsample=sqrt(var(sum(x1.^r,2)))
