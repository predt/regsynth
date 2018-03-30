clear variables;
close all;
rng(0);
parpool('local')

n1 = 10;
n0 = 100;
k = 4;
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

Lambda = 0.01:0.01:1;

% options = optimoptions('quadprog','StepTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',2000,'Display','off');
 options = optimoptions('quadprog','TolFun',1e-10,'TolX',1e-10,'MaxIter',2000,'Display','off');

T = 100;
MSEp = zeros(T,1); MSEnp = zeros(T,1); MSEm = zeros(T,1); MSEm5 = zeros(T,1);
lambdavalues = zeros(T,1);
parfor t =1:T
    tic
    sprintf('Iteration: %d',t)
    
    x1 = a+(b-a)*rand(k,n1);
    x0 = (a-h)+(b-a+2*h)*rand(k,n0);
    y1 = sum(x1.^r,1)'/v+randn(n1,1);
    y0 = sum(x0.^r,1)'/v+randn(n0,1);

    H = 2*(x0'*x0);
    
    minMSE = inf;
    Wp = zeros(n0,n1);
    optlambda = 0;
    
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
    
    y1 = sum(x1.^r,1)'/v+randn(n1,1);
    y0 = sum(x0.^r,1)'/v+randn(n0,1);
    
    Wnp = zeros(n0,n1);
    Wm =  zeros(n0,n1);
    Wm5 = zeros(n0,n1);
    for i=1:n1
        x = x1(:,i);
        D = x0 - kron(ones(1,n0),x);
        delta = diag(D'*D); 
        w = quadprog(H,-2*x0'*x,[],[],ones(1,n0),1,zeros(n0,1),ones(n0,1),[],options);
        Wnp(:,i) = w;
        [sorted,I]=sort(delta);
        Wm(I(1),i)=1;
        Wm5(I(1:5),i)=1/5;
    end
    MSEp(t) = (y1-Wp'*y0)'*(y1-Wp'*y0);
    MSEnp(t) = (y1-Wnp'*y0)'*(y1-Wnp'*y0);
    MSEm(t) = (y1-Wm'*y0)'*(y1-Wm'*y0);
    MSEm5(t) = (y1-Wm5'*y0)'*(y1-Wm5'*y0);
    toc
end

save('n1_100_n0_1000_k_4_r_2_1090_10_nested.mat','n1','n0','k','r','a','b','h','MSEp','MSEnp','MSEm','MSEm5','Lambda','lambdavalues');

% delete(gcp); 

%sum(W>.0000000001,1)

%vsample=sqrt(var(sum(x1.^r,2)))
